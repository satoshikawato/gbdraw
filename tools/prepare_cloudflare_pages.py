#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import shutil
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
WEB_ROOT = REPO_ROOT / "gbdraw" / "web"
DEFAULT_OUTPUT_ROOT = REPO_ROOT / "dist" / "cloudflare-pages"
DEFAULT_GOOGLE_ANALYTICS_MEASUREMENT_ID = "G-GG6JMKM02Y"
GOOGLE_ANALYTICS_MEASUREMENT_ID_ENV = "GBDRAW_GOOGLE_ANALYTICS_MEASUREMENT_ID"
GALLERY_REMOTE_REF_ENV = "GBDRAW_GALLERY_REMOTE_REF"
GALLERY_REMOTE_BASE_ENV = "GBDRAW_GALLERY_REMOTE_BASE"
DEFAULT_GALLERY_REMOTE_REF = "main"
DEFAULT_GALLERY_REMOTE_REPOSITORY = "satoshikawato/gbdraw"
GALLERY_REMOTE_MANIFEST = "gallery/remote-assets.json"
GALLERY_REMOTE_ASSET_LIMIT_BYTES = 25 * 1024 * 1024
GALLERY_REMOTE_ASSET_PATTERNS = (
    "gallery/examples/*.svg",
    "gallery/sessions/*.gbdraw-session.json",
    "gallery/sources/*.svg",
    "gallery/media/**/*",
)
SCRIPT_MARKER = "<!-- GOOGLE_ANALYTICS_SCRIPT -->"
NOTICE_MARKER = "<!-- GOOGLE_ANALYTICS_NOTICE -->"
ISOLATION_HEADERS = """/*
  Cross-Origin-Opener-Policy: same-origin
  Cross-Origin-Embedder-Policy: require-corp
  Cross-Origin-Resource-Policy: same-origin
  Content-Security-Policy: frame-ancestors 'none'

/gallery/examples/*
  ! Content-Security-Policy
  Content-Security-Policy: default-src 'self' data: blob:; script-src 'unsafe-inline'; style-src 'unsafe-inline'; img-src 'self' data: blob:; frame-ancestors 'self'
"""


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


def _load_refresh_gallery_sessions_module():
    return _load_repo_module(REPO_ROOT / "tools" / "refresh_gallery_sessions.py", "refresh_gallery_sessions")


def _load_stamp_web_build_module():
    return _load_repo_module(REPO_ROOT / "tools" / "stamp_web_build.py", "stamp_web_build")


def _replace_once(source: str, old: str, new: str) -> str:
    if old not in source:
        raise RuntimeError(f"Expected text not found while preparing Cloudflare bundle: {old}")
    return source.replace(old, new, 1)


def _render_google_analytics_script(measurement_id: str) -> str:
    return f"""    <!-- Google tag (gtag.js) -->
    <script async src="https://www.googletagmanager.com/gtag/js?id={measurement_id}"></script>
    <script>
      window.dataLayer = window.dataLayer || [];
      function gtag(){{dataLayer.push(arguments);}}
      gtag('js', new Date());
      gtag('config', '{measurement_id}');
    </script>"""


def _render_google_analytics_notice() -> str:
    return """                            <div class="bg-amber-50 text-amber-900 p-2.5 rounded-lg border border-amber-200 leading-relaxed">
                                <strong class="flex items-center gap-1 mb-1"><i class="ph ph-chart-line"></i> Hosted Site Analytics</strong>
                                The hosted <strong>gbdraw.app</strong> deployment uses Google Analytics 4 for aggregate page-usage metrics. Uploaded genome files and generated diagrams are still processed locally in your browser and are not sent to Google Analytics by gbdraw.
                            </div>"""


def _default_gallery_remote_base() -> str:
    explicit_base = os.environ.get(GALLERY_REMOTE_BASE_ENV)
    if explicit_base:
        return explicit_base.rstrip("/") + "/"
    ref = (
        os.environ.get(GALLERY_REMOTE_REF_ENV)
        or os.environ.get("CF_PAGES_COMMIT_SHA")
        or os.environ.get("GITHUB_SHA")
        or DEFAULT_GALLERY_REMOTE_REF
    )
    return f"https://raw.githubusercontent.com/{DEFAULT_GALLERY_REMOTE_REPOSITORY}/{ref}/gbdraw/web/"


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
    if output_root.exists():
        shutil.rmtree(output_root)
    shutil.copytree(WEB_ROOT, output_root)

    index_path = output_root / "index.html"
    index_html = index_path.read_text(encoding="utf-8")
    index_html = _replace_once(
        index_html,
        SCRIPT_MARKER,
        _render_google_analytics_script(google_analytics_measurement_id),
    )
    index_html = _replace_once(index_html, NOTICE_MARKER, _render_google_analytics_notice())
    index_html = _replace_once(
        index_html,
        "script-src 'self' 'unsafe-inline' 'unsafe-eval';",
        "script-src 'self' 'unsafe-inline' 'unsafe-eval' https://*.googletagmanager.com;",
    )
    index_html = _replace_once(
        index_html,
        "img-src 'self' data: blob:;",
        "img-src 'self' data: blob: https://*.google-analytics.com https://*.googletagmanager.com;",
    )
    index_html = _replace_once(
        index_html,
        "connect-src 'self';",
        "connect-src 'self' https://*.google-analytics.com https://*.analytics.google.com https://*.googletagmanager.com;",
    )
    index_path.write_text(index_html, encoding="utf-8")
    _load_stamp_web_build_module().stamp_web_build(output_root, commit_sha=commit_sha)
    (output_root / "_headers").write_text(ISOLATION_HEADERS, encoding="utf-8")
    _write_remote_gallery_manifest(
        output_root=output_root,
        remote_base=gallery_remote_base or _default_gallery_remote_base(),
    )
    return output_root


def prepare_cloudflare_pages(
    *,
    refresh_cache_bust: bool = False,
    refresh_gallery_sessions: bool = False,
    google_analytics_measurement_id: str | None = None,
    output_root: Path = DEFAULT_OUTPUT_ROOT,
) -> Path:
    if refresh_gallery_sessions:
        refresh_gallery_sessions_module = _load_refresh_gallery_sessions_module()
        refresh_gallery_sessions_module.refresh_gallery_sessions()
        refresh_gallery_sessions_module.prepare_gallery_assets()
    prepare_browser_wheel_module = _load_prepare_browser_wheel_module()
    prepare_browser_wheel_module.prepare_browser_wheel(refresh_cache_bust=refresh_cache_bust)
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
            "Google Analytics 4 tag, hosted-site privacy notice, and build label."
        )
    )
    parser.add_argument(
        "--refresh-cache-bust",
        action="store_true",
        help="Refresh GBDRAW_WHEEL_CACHE_BUST while preparing the browser wheel.",
    )
    parser.add_argument(
        "--google-analytics-id",
        help=(
            "Google Analytics measurement ID. Defaults to the value of "
            f"{GOOGLE_ANALYTICS_MEASUREMENT_ID_ENV} or the repository default."
        ),
    )
    parser.add_argument(
        "--refresh-gallery",
        action="store_true",
        help=(
            "Refresh web gallery session JSON/SVG assets before copying the web bundle. "
            "This requires the optional gallery asset dependencies."
        ),
    )
    args = parser.parse_args(argv)
    output_root = prepare_cloudflare_pages(
        refresh_cache_bust=args.refresh_cache_bust,
        refresh_gallery_sessions=args.refresh_gallery,
        google_analytics_measurement_id=args.google_analytics_id,
    )
    print(f"Prepared Cloudflare Pages bundle: {output_root.relative_to(REPO_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
