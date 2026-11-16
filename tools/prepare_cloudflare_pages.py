#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import os
import posixpath
import re
import shutil
from collections import deque
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
WEB_ROOT = REPO_ROOT / "gbdraw" / "web"
DEFAULT_OUTPUT_ROOT = REPO_ROOT / "dist" / "cloudflare-pages"
DEFAULT_GOOGLE_ANALYTICS_MEASUREMENT_ID = "G-GG6JMKM02Y"
GOOGLE_ANALYTICS_MEASUREMENT_ID_ENV = "GBDRAW_GOOGLE_ANALYTICS_MEASUREMENT_ID"
GALLERY_REMOTE_REF_ENV = "GBDRAW_GALLERY_REMOTE_REF"
GALLERY_REMOTE_BASE_ENV = "GBDRAW_GALLERY_REMOTE_BASE"
DEFAULT_GALLERY_REMOTE_REPOSITORY = "satoshikawato/gbdraw"
GALLERY_REMOTE_MANIFEST = "gallery/remote-assets.json"
GALLERY_REMOTE_ASSET_LIMIT_BYTES = 25 * 1024 * 1024
GALLERY_REMOTE_ASSET_PATTERNS = (
    "gallery/examples/*.svg",
    "gallery/sessions/*.gbdraw",
    "gallery/sources/*.svg",
    "gallery/media/**/*",
)
SCRIPT_MARKER = "<!-- GOOGLE_ANALYTICS_SCRIPT -->"
# The flag that tells the page it is the hosted site, which loads Google
# Analytics; the About section shows its analytics notice from it
# (gbdraw/web/js/utils/hosted-build.js reads it by this name).
HOSTED_ANALYTICS_META_NAME = "gbdraw-hosted-analytics"
HOSTED_ANALYTICS_FLAG = f'<meta name="{HOSTED_ANALYTICS_META_NAME}" content="google-analytics-4">'
# The hosted pages that load Google Analytics, each with SCRIPT_MARKER in its
# head, and the same-origin script, written only into the bundle, that
# configures the tag for all of them (the Gallery's CSP allows no inline script).
GOOGLE_ANALYTICS_PAGES = ("index.html", "gallery/index.html", "gallery/palettes/index.html")
GOOGLE_ANALYTICS_BOOTSTRAP = "hosted-analytics.js"
# The hosts the Google tag needs, by CSP directive.
GOOGLE_ANALYTICS_CSP_SOURCES = {
    "script-src": "https://*.googletagmanager.com",
    "img-src": "https://*.google-analytics.com https://*.googletagmanager.com",
    "connect-src": (
        "https://*.google-analytics.com https://*.analytics.google.com "
        "https://*.googletagmanager.com"
    ),
}
_CSP_META_RE = re.compile(r'<meta http-equiv="Content-Security-Policy" content="([^"]*)"')
FULL_COMMIT_SHA_RE = re.compile(r"^[0-9a-fA-F]{40}$")
ISOLATION_HEADERS = """/*
  Cross-Origin-Opener-Policy: same-origin
  Cross-Origin-Embedder-Policy: require-corp
  Cross-Origin-Resource-Policy: same-origin
  Content-Security-Policy: frame-ancestors 'none'

/gallery/examples/*
  ! Content-Security-Policy
  Content-Security-Policy: default-src 'self' data: blob:; script-src 'unsafe-inline'; style-src 'unsafe-inline'; img-src 'self' data: blob:; frame-ancestors 'self'
"""
# Vendored files whose URL never serves other bytes: browsers keep them a year
# without revalidating. The manifest pins each file's SHA-256, and the bundle
# build fails when a file under these paths changes without a new URL.
IMMUTABLE_ASSETS_MANIFEST = REPO_ROOT / "tools" / "immutable-web-assets.json"
IMMUTABLE_CACHE_CONTROL = "public, max-age=31536000, immutable"


def _read_immutable_assets_manifest() -> tuple[tuple[str, ...], dict[str, str]]:
    manifest = json.loads(IMMUTABLE_ASSETS_MANIFEST.read_text(encoding="utf-8"))
    return tuple(manifest["prefixes"]), dict(manifest["files"])


def render_bundle_headers() -> str:
    """The bundle's `_headers`: isolation for every path, a year for each pinned file.

    One rule per pinned file, not per folder: Workers static assets answer a
    missing path with the single-page-app `index.html`, and a folder rule would
    let that page be cached as immutable.
    """
    _, pinned = _read_immutable_assets_manifest()
    immutable_rules = "".join(
        f"\n/{name}\n  Cache-Control: {IMMUTABLE_CACHE_CONTROL}\n" for name in sorted(pinned)
    )
    return ISOLATION_HEADERS + immutable_rules


def verify_immutable_assets(web_root: Path) -> None:
    """Fail when a file under an immutable path is not the one its URL was pinned to."""
    prefixes, pinned = _read_immutable_assets_manifest()
    found = {
        path.relative_to(web_root).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
        for prefix in prefixes
        for path in sorted((web_root / prefix).rglob("*"))
        if path.is_file()
    }
    problems = [
        f"{name}: {'not in the manifest' if name not in pinned else 'bytes changed'} ({digest})"
        for name, digest in sorted(found.items())
        if pinned.get(name) != digest
    ] + [f"{name}: missing" for name in sorted(set(pinned) - set(found))]
    if problems:
        raise RuntimeError(
            "Files under immutable URLs differ from "
            f"{IMMUTABLE_ASSETS_MANIFEST.relative_to(REPO_ROOT)}:\n"
            + "\n".join(problems)
            + "\nA file whose bytes change needs a new URL (a new file name or folder); "
            "then record the new path and SHA-256 in the manifest."
        )


# The page's module entry and its static imports (`import ... from`, a bare
# `import '...'`, `export ... from`); dynamic `import()` is not preloaded.
# tests/web/module-preload-closure.test.mjs compares the closure with V8's parser.
MODULE_ENTRY = "js/app.js"
_STATIC_IMPORT_RE = re.compile(
    r"""(?:^|[;}])[ \t]*(?:import|export)[ \t]+(?:[\w$*{}\s,]+?\s+from\s*)?['"]([^'"\n]+)['"]""",
    re.MULTILINE,
)


def module_preload_paths(web_root: Path) -> list[str]:
    """The modules js/app.js imports statically, nearest first, without the entry."""
    seen = {MODULE_ENTRY}
    order: list[str] = []
    pending = deque([MODULE_ENTRY])
    while pending:
        module = pending.popleft()
        source = (web_root / module).read_text(encoding="utf-8")
        for specifier in _STATIC_IMPORT_RE.findall(source):
            if not specifier.startswith("."):
                raise RuntimeError(f"{module} imports {specifier!r}, which is not a relative path")
            target = posixpath.normpath(posixpath.join(posixpath.dirname(module), specifier))
            if target not in seen:
                seen.add(target)
                order.append(target)
                pending.append(target)
    return order


def add_module_preloads(index_html: str, web_root: Path) -> str:
    """Let the browser fetch the whole module graph at once instead of level by level."""
    links = "".join(
        f'    <link rel="modulepreload" href="./{path}">\n' for path in module_preload_paths(web_root)
    )
    return _replace_once(index_html, "</head>", f"{links}</head>")


def _load_repo_module(module_path: Path, module_name: str):
    spec = spec_from_file_location(module_name, module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load module from {module_path}")
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _load_prepare_browser_wheel_module():
    return _load_repo_module(REPO_ROOT / "tools" / "prepare_browser_wheel.py", "prepare_browser_wheel")


def _load_build_support_module():
    return _load_repo_module(REPO_ROOT / "gbdraw" / "_build_support.py", "gbdraw_build_support")


def _load_stamp_web_build_module():
    return _load_repo_module(REPO_ROOT / "tools" / "stamp_web_build.py", "stamp_web_build")


def _load_gallery_artifact_manifest_module():
    return _load_repo_module(
        REPO_ROOT / "tools" / "gallery_artifact_manifest.py",
        "gallery_artifact_manifest",
    )


def _replace_once(source: str, old: str, new: str) -> str:
    if old not in source:
        raise RuntimeError(f"Expected text not found while preparing Cloudflare bundle: {old}")
    return source.replace(old, new, 1)


def _render_google_analytics_bootstrap(measurement_id: str) -> str:
    return f"""// The Google tag's configuration on the hosted site (tools/prepare_cloudflare_pages.py).
window.dataLayer = window.dataLayer || [];
function gtag(){{dataLayer.push(arguments);}}
gtag('js', new Date());
gtag('config', '{measurement_id}');
"""


def _render_google_analytics_script(measurement_id: str, page: str) -> str:
    bootstrap = posixpath.relpath(GOOGLE_ANALYTICS_BOOTSTRAP, posixpath.dirname(page) or ".")
    return f"""    {HOSTED_ANALYTICS_FLAG}
    <!-- Google tag (gtag.js) -->
    <script async src="https://www.googletagmanager.com/gtag/js?id={measurement_id}"></script>
    <script async src="{bootstrap}"></script>"""


def _allow_google_analytics_in_csp(page_html: str) -> str:
    """Add the Google tag's hosts to the page's CSP; a missing directive starts from default-src."""
    match = _CSP_META_RE.search(page_html)
    if match is None:
        raise RuntimeError("Expected a Content-Security-Policy meta element while preparing Cloudflare bundle")
    policy = match.group(1)
    for directive, hosts in GOOGLE_ANALYTICS_CSP_SOURCES.items():
        found = re.search(rf"(?<![\w-]){directive}\s([^;]*);", policy)
        if found is not None:
            policy = f"{policy[:found.end(1)]} {hosts}{policy[found.end(1):]}"
            continue
        default = re.search(r"(?<![\w-])default-src\s([^;]*);", policy)
        if default is None:
            raise RuntimeError(f"The page's CSP has neither {directive} nor default-src")
        policy = f"{policy[:default.end()]} {directive} {default.group(1)} {hosts};{policy[default.end():]}"
    return page_html[: match.start(1)] + policy + page_html[match.end(1) :]


def _default_gallery_remote_base(*, commit_sha: str | None = None) -> str:
    explicit_base = os.environ.get(GALLERY_REMOTE_BASE_ENV)
    if explicit_base:
        return explicit_base.rstrip("/") + "/"

    ref = (
        os.environ.get(GALLERY_REMOTE_REF_ENV, "").strip()
        or str(commit_sha or "").strip()
        or _load_stamp_web_build_module().resolve_commit_sha().strip()
    )
    if not FULL_COMMIT_SHA_RE.fullmatch(ref):
        raise RuntimeError(
            "Cloudflare Gallery remote assets require a full 40-character commit SHA. "
            f"Set {GALLERY_REMOTE_REF_ENV} to an immutable commit SHA."
        )
    return (
        f"https://raw.githubusercontent.com/{DEFAULT_GALLERY_REMOTE_REPOSITORY}/"
        f"{ref}/gbdraw/web/"
    )


def _write_remote_gallery_manifest(output_root: Path, remote_base: str) -> dict[str, str]:
    remote_assets: dict[str, str] = {}
    for pattern in GALLERY_REMOTE_ASSET_PATTERNS:
        for path in sorted(output_root.glob(pattern)):
            if not path.is_file() or path.stat().st_size <= GALLERY_REMOTE_ASSET_LIMIT_BYTES:
                continue
            relative_path = path.relative_to(output_root).as_posix()
            remote_assets[relative_path] = f"{remote_base}{relative_path}"
            path.unlink()

    manifest_path = output_root / GALLERY_REMOTE_MANIFEST
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(
        json.dumps(remote_assets, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return remote_assets


def build_cloudflare_pages_bundle(
    *,
    output_root: Path = DEFAULT_OUTPUT_ROOT,
    google_analytics_measurement_id: str = DEFAULT_GOOGLE_ANALYTICS_MEASUREMENT_ID,
    gallery_remote_base: str | None = None,
    commit_sha: str | None = None,
) -> Path:
    verify_gallery_artifacts = (
        _load_gallery_artifact_manifest_module().verify_gallery_artifacts
    )

    verify_gallery_artifacts()
    if output_root.exists():
        shutil.rmtree(output_root)
    shutil.copytree(WEB_ROOT, output_root)

    (output_root / GOOGLE_ANALYTICS_BOOTSTRAP).write_text(
        _render_google_analytics_bootstrap(google_analytics_measurement_id), encoding="utf-8"
    )
    for page in GOOGLE_ANALYTICS_PAGES:
        page_path = output_root / page
        page_html = _replace_once(
            page_path.read_text(encoding="utf-8"),
            SCRIPT_MARKER,
            _render_google_analytics_script(google_analytics_measurement_id, page),
        )
        page_path.write_text(_allow_google_analytics_in_csp(page_html), encoding="utf-8")
    index_path = output_root / "index.html"
    index_path.write_text(
        add_module_preloads(index_path.read_text(encoding="utf-8"), output_root), encoding="utf-8"
    )
    stamp_web_build_module = _load_stamp_web_build_module()
    resolved_commit_sha = str(commit_sha or stamp_web_build_module.resolve_commit_sha()).strip()
    stamp_web_build_module.stamp_web_build(
        output_root,
        commit_sha=resolved_commit_sha,
    )
    build_support = _load_build_support_module()
    # The token hashes the wheel as served: any build-time repack of the wheel
    # must happen before this point.
    served_wheel = output_root / build_support.expected_browser_wheel_name()
    build_support.write_wheel_cache_bust(
        output_root / "js" / "config.js",
        build_support.browser_wheel_cache_bust_token(served_wheel),
    )
    verify_immutable_assets(output_root)
    (output_root / "_headers").write_text(render_bundle_headers(), encoding="utf-8")
    _write_remote_gallery_manifest(
        output_root=output_root,
        remote_base=gallery_remote_base
        or _default_gallery_remote_base(commit_sha=resolved_commit_sha),
    )
    verify_gallery_artifacts(
        package_root=output_root,
        commit_sha=resolved_commit_sha.lower(),
    )
    return output_root


def prepare_cloudflare_pages(
    *,
    google_analytics_measurement_id: str | None = None,
    output_root: Path = DEFAULT_OUTPUT_ROOT,
) -> Path:
    prepare_browser_wheel_module = _load_prepare_browser_wheel_module()
    prepare_browser_wheel_module.prepare_browser_wheel()
    build_support = _load_build_support_module()
    build_support.refresh_open_source_notices()
    measurement_id = (
        google_analytics_measurement_id
        or os.environ.get(GOOGLE_ANALYTICS_MEASUREMENT_ID_ENV)
        or DEFAULT_GOOGLE_ANALYTICS_MEASUREMENT_ID
    )
    return build_cloudflare_pages_bundle(
        output_root=output_root,
        google_analytics_measurement_id=measurement_id,
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare the Cloudflare Pages bundle in dist/cloudflare-pages with the browser wheel, "
            "Google Analytics 4 tag, the hosted-site analytics flag, build label, and the "
            "wheel's content-hash cache-bust token."
        )
    )
    parser.add_argument(
        "--google-analytics-id",
        help=(
            "Google Analytics measurement ID. Defaults to the value of "
            f"{GOOGLE_ANALYTICS_MEASUREMENT_ID_ENV} or the repository default."
        ),
    )
    args = parser.parse_args(argv)
    output_root = prepare_cloudflare_pages(
        google_analytics_measurement_id=args.google_analytics_id,
    )
    print(f"Prepared Cloudflare Pages bundle: {output_root.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
