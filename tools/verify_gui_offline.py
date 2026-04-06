#!/usr/bin/env python3
from __future__ import annotations

import argparse
import contextlib
import http.server
import re
import shutil
import socketserver
import subprocess
import sys
import tarfile
import tempfile
import threading
import urllib.error
import urllib.request
import zipfile
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
WEB_ROOT = REPO_ROOT / "gbdraw" / "web"
VENDOR_ROOT = WEB_ROOT / "vendor"
CONFIG_PATH = WEB_ROOT / "js" / "config.js"

ASSET_URLS = {
    "vue": "https://unpkg.com/vue@3.5.25/dist/vue.global.js",
    "tailwindcss": "https://cdn.tailwindcss.com",
    "pyodide_core": "https://github.com/pyodide/pyodide/releases/download/0.29.0/pyodide-core-0.29.0.tar.bz2",
    "browser_wasi_shim": "https://registry.npmjs.org/@bjorn3/browser_wasi_shim/-/browser_wasi_shim-0.4.2.tgz",
    "jspdf": "https://registry.npmjs.org/jspdf/-/jspdf-3.0.3.tgz",
    "svg2pdf": "https://registry.npmjs.org/svg2pdf.js/-/svg2pdf.js-2.6.0.tgz",
    "dompurify": "https://registry.npmjs.org/dompurify/-/dompurify-3.2.7.tgz",
    "phosphor": "https://registry.npmjs.org/@phosphor-icons/web/-/web-2.1.2.tgz",
    "inter": "https://registry.npmjs.org/@fontsource/inter/-/inter-5.2.8.tgz",
    "noto_sans_jp": "https://registry.npmjs.org/@fontsource/noto-sans-jp/-/noto-sans-jp-5.2.9.tgz",
}

PYODIDE_CORE_FILES = {
    "package.json",
    "pyodide.asm.js",
    "pyodide.asm.wasm",
    "pyodide.js",
    "pyodide.mjs",
    "pyodide-lock.json",
    "python_stdlib.zip",
}

UI_FONT_ASSETS = {
    "inter": (
        "inter-latin-300-normal.woff2",
        "inter-latin-400-normal.woff2",
        "inter-latin-500-normal.woff2",
        "inter-latin-600-normal.woff2",
        "inter-latin-700-normal.woff2",
    ),
    "noto-sans-jp": (
        "noto-sans-jp-japanese-300-normal.woff2",
        "noto-sans-jp-japanese-400-normal.woff2",
        "noto-sans-jp-japanese-500-normal.woff2",
        "noto-sans-jp-japanese-700-normal.woff2",
    ),
}

REQUIRED_UI_FONT_FILES = tuple(
    Path("vendor") / "fonts" / family / filename
    for family, filenames in UI_FONT_ASSETS.items()
    for filename in filenames
)
TOP_LEVEL_BROWSER_WHEEL_RE = re.compile(r"^gbdraw/web/[^/]+\.whl$")


def _parse_local_wheel_paths() -> tuple[Path, ...]:
    config_text = CONFIG_PATH.read_text(encoding="utf-8")
    match = re.search(r"export const PYODIDE_LOCAL_WHEELS\s*=\s*\[(.*?)\];", config_text, re.DOTALL)
    if match is None:
        raise RuntimeError("Could not determine PYODIDE_LOCAL_WHEELS from gbdraw/web/js/config.js")
    wheel_paths = tuple(
        Path(raw.lstrip("./"))
        for raw in re.findall(r"""["']([^"']+\.whl)["']""", match.group(1))
    )
    if not wheel_paths:
        raise RuntimeError("PYODIDE_LOCAL_WHEELS is empty in gbdraw/web/js/config.js")
    return wheel_paths


def _download(url: str, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": "gbdraw-offline-vendor/1.0",
        },
    )
    with urllib.request.urlopen(request, timeout=120) as response, target.open("wb") as dst:
        shutil.copyfileobj(response, dst)


def _extract_tgz_member(archive: Path, member_name: str, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive, "r:gz") as tar:
        member = tar.getmember(member_name)
        extracted = tar.extractfile(member)
        if extracted is None:
            raise FileNotFoundError(f"Missing member {member_name} in {archive.name}")
        with destination.open("wb") as dst:
            shutil.copyfileobj(extracted, dst)


def _extract_tgz_prefix(archive: Path, prefix: str, destination_dir: Path) -> None:
    destination_dir.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive, "r:gz") as tar:
        for member in tar.getmembers():
            if not member.isfile() or not member.name.startswith(prefix):
                continue
            relative = Path(member.name).relative_to(prefix)
            if relative.name == "":
                continue
            target = destination_dir / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            extracted = tar.extractfile(member)
            if extracted is None:
                continue
            with target.open("wb") as dst:
                shutil.copyfileobj(extracted, dst)


def _extract_tgz_files(archive: Path, members: tuple[str, ...], destination_dir: Path) -> None:
    for member_name in members:
        _extract_tgz_member(archive, f"package/files/{member_name}", destination_dir / member_name)


def vendor_assets() -> None:
    with tempfile.TemporaryDirectory(prefix="gbdraw-offline-assets-") as tmpdir:
        tmp = Path(tmpdir)
        for path in [
            VENDOR_ROOT / "vue",
            VENDOR_ROOT / "tailwindcss",
            VENDOR_ROOT / "pyodide" / "v0.29.0" / "full",
            VENDOR_ROOT / "browser_wasi_shim" / "dist",
            VENDOR_ROOT / "jspdf",
            VENDOR_ROOT / "svg2pdf.js",
            VENDOR_ROOT / "dompurify",
            VENDOR_ROOT / "phosphor-icons" / "regular",
            VENDOR_ROOT / "fonts" / "inter",
            VENDOR_ROOT / "fonts" / "noto-sans-jp",
        ]:
            path.mkdir(parents=True, exist_ok=True)

        vue_archive = tmp / "vue.global.js"
        _download(ASSET_URLS["vue"], vue_archive)
        shutil.copy2(vue_archive, VENDOR_ROOT / "vue" / "vue.global.js")

        tailwind_archive = tmp / "tailwindcss-play.js"
        _download(ASSET_URLS["tailwindcss"], tailwind_archive)
        shutil.copy2(tailwind_archive, VENDOR_ROOT / "tailwindcss" / "tailwindcss-play.js")

        pyodide_archive = tmp / "pyodide-core.tar.bz2"
        _download(ASSET_URLS["pyodide_core"], pyodide_archive)
        pyodide_target = VENDOR_ROOT / "pyodide" / "v0.29.0" / "full"
        pyodide_target.mkdir(parents=True, exist_ok=True)
        with tarfile.open(pyodide_archive, "r:bz2") as tar:
            for member in tar.getmembers():
                if not member.isfile():
                    continue
                name = Path(member.name).name
                if name not in PYODIDE_CORE_FILES:
                    continue
                extracted = tar.extractfile(member)
                if extracted is None:
                    continue
                with (pyodide_target / name).open("wb") as dst:
                    shutil.copyfileobj(extracted, dst)

        browser_wasi_archive = tmp / "browser_wasi_shim.tgz"
        _download(ASSET_URLS["browser_wasi_shim"], browser_wasi_archive)
        _extract_tgz_prefix(
            browser_wasi_archive,
            "package/dist",
            VENDOR_ROOT / "browser_wasi_shim" / "dist",
        )

        jspdf_archive = tmp / "jspdf.tgz"
        _download(ASSET_URLS["jspdf"], jspdf_archive)
        _extract_tgz_member(
            jspdf_archive,
            "package/dist/jspdf.umd.min.js",
            VENDOR_ROOT / "jspdf" / "jspdf.umd.min.js",
        )

        svg2pdf_archive = tmp / "svg2pdf.tgz"
        _download(ASSET_URLS["svg2pdf"], svg2pdf_archive)
        _extract_tgz_member(
            svg2pdf_archive,
            "package/dist/svg2pdf.umd.min.js",
            VENDOR_ROOT / "svg2pdf.js" / "svg2pdf.umd.min.js",
        )

        dompurify_archive = tmp / "dompurify.tgz"
        _download(ASSET_URLS["dompurify"], dompurify_archive)
        _extract_tgz_member(
            dompurify_archive,
            "package/dist/purify.min.js",
            VENDOR_ROOT / "dompurify" / "purify.min.js",
        )

        phosphor_archive = tmp / "phosphor-icons.tgz"
        _download(ASSET_URLS["phosphor"], phosphor_archive)
        _extract_tgz_prefix(
            phosphor_archive,
            "package/src/regular",
            VENDOR_ROOT / "phosphor-icons" / "regular",
        )

        inter_archive = tmp / "inter.tgz"
        _download(ASSET_URLS["inter"], inter_archive)
        _extract_tgz_files(inter_archive, UI_FONT_ASSETS["inter"], VENDOR_ROOT / "fonts" / "inter")

        noto_sans_jp_archive = tmp / "noto-sans-jp.tgz"
        _download(ASSET_URLS["noto_sans_jp"], noto_sans_jp_archive)
        _extract_tgz_files(
            noto_sans_jp_archive,
            UI_FONT_ASSETS["noto-sans-jp"],
            VENDOR_ROOT / "fonts" / "noto-sans-jp",
        )


def _parse_wheel_name() -> str:
    config_text = CONFIG_PATH.read_text(encoding="utf-8")
    for line in config_text.splitlines():
        line = line.strip()
        if not line.startswith("export const GBDRAW_WHEEL_NAME"):
            continue
        _, value = line.split("=", 1)
        return value.strip().strip(";").strip().strip('"').strip("'")
    raise RuntimeError("Could not determine GBDRAW_WHEEL_NAME from gbdraw/web/js/config.js")


def _source_top_level_browser_wheels() -> tuple[Path, ...]:
    return tuple(sorted(path for path in WEB_ROOT.glob("*.whl") if path.is_file()))


def _assert_source_browser_wheel_state() -> Path:
    expected_name = _parse_wheel_name()
    top_level_wheels = _source_top_level_browser_wheels()
    expected_path = WEB_ROOT / expected_name
    if expected_path not in top_level_wheels:
        raise FileNotFoundError(f"Configured browser wheel is missing: {expected_path.relative_to(REPO_ROOT)}")
    extras = [path for path in top_level_wheels if path.name != expected_name]
    if extras:
        raise RuntimeError(
            "Unreferenced top-level browser wheels remain in gbdraw/web:\n"
            + "\n".join(str(path.relative_to(REPO_ROOT)) for path in extras)
        )
    return expected_path


def _read_wheel_file_bytes(wheel_path: Path) -> dict[str, bytes]:
    with zipfile.ZipFile(wheel_path) as zf:
        return {
            info.filename: zf.read(info.filename)
            for info in zf.infolist()
            if not info.is_dir()
        }


def _repo_python_sources() -> dict[str, bytes]:
    return {
        str(path.relative_to(REPO_ROOT)).replace("\\", "/"): path.read_bytes()
        for path in sorted((REPO_ROOT / "gbdraw").rglob("*.py"))
        if "__pycache__" not in path.parts
    }


def _assert_wheel_python_sources_match_source_tree(wheel_path: Path) -> None:
    wheel_entries = _read_wheel_file_bytes(wheel_path)
    wheel_python_sources = {
        name: data
        for name, data in wheel_entries.items()
        if name.startswith("gbdraw/") and name.endswith(".py")
    }
    repo_python_sources = _repo_python_sources()

    missing = sorted(set(repo_python_sources) - set(wheel_python_sources))
    extra = sorted(set(wheel_python_sources) - set(repo_python_sources))
    mismatched = sorted(
        name
        for name in repo_python_sources
        if name in wheel_python_sources and repo_python_sources[name] != wheel_python_sources[name]
    )
    if missing or extra or mismatched:
        details: list[str] = []
        if missing:
            details.append("Missing python sources in wheel:\n" + "\n".join(missing))
        if extra:
            details.append("Unexpected python sources in wheel:\n" + "\n".join(extra))
        if mismatched:
            details.append("Mismatched python sources in wheel:\n" + "\n".join(mismatched))
        details.append(
            "Refresh the checked-in browser wheel with:\n"
            "python -m build --wheel --no-isolation\n"
            "python tools/sync_browser_wheel.py dist/*.whl"
        )
        raise RuntimeError("\n\n".join(details))


def _outer_wheel_top_level_browser_entries(wheel_path: Path) -> tuple[str, ...]:
    entries = _read_wheel_file_bytes(wheel_path)
    return tuple(sorted(name for name in entries if TOP_LEVEL_BROWSER_WHEEL_RE.match(name)))


def _extract_embedded_browser_wheel(outer_wheel_path: Path) -> bytes:
    browser_wheel_name = _parse_wheel_name()
    entry_name = f"gbdraw/web/{browser_wheel_name}"
    entries = _read_wheel_file_bytes(outer_wheel_path)
    if entry_name not in entries:
        raise FileNotFoundError(f"Embedded browser wheel missing from {outer_wheel_path}: {entry_name}")
    return entries[entry_name]


def _assert_outer_wheel_embeds_source_browser_wheel(outer_wheel_path: Path) -> None:
    source_browser_wheel = _assert_source_browser_wheel_state()
    expected_entries = (f"gbdraw/web/{source_browser_wheel.name}",)
    actual_entries = _outer_wheel_top_level_browser_entries(outer_wheel_path)
    if actual_entries != expected_entries:
        raise RuntimeError(
            "Outer wheel contains unexpected top-level browser wheels:\n"
            + "\n".join(actual_entries or ["(none)"])
        )
    embedded_browser_wheel = _extract_embedded_browser_wheel(outer_wheel_path)
    source_browser_wheel_bytes = source_browser_wheel.read_bytes()
    if embedded_browser_wheel != source_browser_wheel_bytes:
        raise RuntimeError(
            "Embedded browser wheel does not match gbdraw/web browser wheel:\n"
            f"outer={outer_wheel_path}\n"
            f"source={source_browser_wheel}"
        )


def _assert_outer_wheel_config_matches_source(outer_wheel_path: Path) -> None:
    entries = _read_wheel_file_bytes(outer_wheel_path)
    entry_name = "gbdraw/web/js/config.js"
    if entry_name not in entries:
        raise FileNotFoundError(f"Outer wheel is missing {entry_name}")
    source_config = CONFIG_PATH.read_bytes()
    if entries[entry_name] != source_config:
        raise RuntimeError(
            "Outer wheel config.js does not match gbdraw/web/js/config.js:\n"
            f"outer={outer_wheel_path}\n"
            f"source={CONFIG_PATH}"
        )


def _assert_packaged_assets() -> None:
    _assert_source_browser_wheel_state()
    required = [
        WEB_ROOT / "index.html",
        WEB_ROOT / "open-source-notices.html",
        WEB_ROOT / "js" / "app.js",
        WEB_ROOT / "js" / "services" / "losat.js",
        WEB_ROOT / "vendor" / "vue" / "vue.global.js",
        WEB_ROOT / "vendor" / "tailwindcss" / "tailwindcss-play.js",
        *(WEB_ROOT / path for path in REQUIRED_UI_FONT_FILES),
        WEB_ROOT / "vendor" / "pyodide" / "v0.29.0" / "full" / "pyodide.js",
        WEB_ROOT / "vendor" / "pyodide" / "v0.29.0" / "full" / "pyodide.asm.js",
        WEB_ROOT / "vendor" / "pyodide" / "v0.29.0" / "full" / "pyodide.asm.wasm",
        WEB_ROOT / "vendor" / "pyodide" / "v0.29.0" / "full" / "pyodide-lock.json",
        WEB_ROOT / "vendor" / "pyodide" / "v0.29.0" / "full" / "python_stdlib.zip",
        WEB_ROOT / "vendor" / "pyodide" / "v0.29.0" / "full" / "micropip-0.11.0-py3-none-any.whl",
        WEB_ROOT / "vendor" / "browser_wasi_shim" / "dist" / "index.js",
        WEB_ROOT / "vendor" / "browser_wasi_shim" / "dist" / "wasi.js",
        WEB_ROOT / "vendor" / "jspdf" / "jspdf.umd.min.js",
        WEB_ROOT / "vendor" / "svg2pdf.js" / "svg2pdf.umd.min.js",
        WEB_ROOT / "vendor" / "dompurify" / "purify.min.js",
        WEB_ROOT / "vendor" / "phosphor-icons" / "regular" / "style.css",
        WEB_ROOT / "wasm" / "losat" / "losat.wasm",
        WEB_ROOT / _parse_wheel_name(),
        *(WEB_ROOT / path for path in _parse_local_wheel_paths()),
    ]
    missing = [path for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(
            "Missing required offline assets:\n" + "\n".join(str(path.relative_to(REPO_ROOT)) for path in missing)
        )


class QuietSimpleHTTPRequestHandler(http.server.SimpleHTTPRequestHandler):
    def log_message(self, format: str, *args: object) -> None:
        return


@contextlib.contextmanager
def _serve_web_root():
    handler = lambda *args, **kwargs: QuietSimpleHTTPRequestHandler(*args, directory=str(WEB_ROOT), **kwargs)
    with socketserver.TCPServer(("127.0.0.1", 0), handler) as httpd:
        port = httpd.server_address[1]
        thread = threading.Thread(target=httpd.serve_forever, daemon=True)
        thread.start()
        try:
            yield f"http://127.0.0.1:{port}"
        finally:
            httpd.shutdown()
            thread.join(timeout=5)


def _ensure_playwright_available() -> None:
    try:
        import playwright  # noqa: F401
    except ImportError as exc:
        raise RuntimeError(
            "Playwright is not installed. Run `python -m pip install playwright` and "
            "`python -m playwright install chromium` first."
        ) from exc


def smoke_test() -> None:
    _assert_packaged_assets()
    _ensure_playwright_available()

    from playwright.sync_api import sync_playwright

    circular_gbk = (REPO_ROOT / "tests" / "test_inputs" / "HmmtDNA.gbk").read_text(encoding="utf-8")
    linear_left_gbk = (REPO_ROOT / "tests" / "test_inputs" / "MERS-CoV.gbk").read_text(encoding="utf-8")
    linear_right_gbk = (REPO_ROOT / "tests" / "test_inputs" / "SARS-CoV-1.gbk").read_text(encoding="utf-8")

    with _serve_web_root() as base_url, sync_playwright() as p:
        try:
            browser = p.chromium.launch(headless=True)
        except Exception as exc:  # pragma: no cover - environment dependent
            raise RuntimeError(
                "Could not launch Playwright Chromium for offline GUI verification. "
                f"Underlying error: {exc}"
            ) from exc
        context = browser.new_context(accept_downloads=True)
        external_requests: list[str] = []

        def route_handler(route):
            url = route.request.url
            if url.startswith(base_url) or url.startswith("blob:") or url.startswith("data:") or url.startswith("about:"):
                route.continue_()
                return
            external_requests.append(url)
            route.abort()

        context.route("**/*", route_handler)
        page = context.new_page()
        page.goto(base_url, wait_until="domcontentloaded")
        page.wait_for_function(
            """
            () => {
              const app = window.__GBDRAW_APP__;
              if (!app) return false;
              const status = String(app.loadingStatus || '');
              return app.pyodideReady === true || status.startsWith('Startup Error:');
            }
            """,
            timeout=120000,
        )
        startup_state = page.evaluate(
            "() => ({ pyodideReady: window.__GBDRAW_APP__.pyodideReady, loadingStatus: window.__GBDRAW_APP__.loadingStatus })"
        )
        if not startup_state["pyodideReady"]:
            raise RuntimeError(f"GUI startup failed offline: {startup_state['loadingStatus']}")

        page.evaluate(
            """
            async ({ gbText }) => {
              const app = window.__GBDRAW_APP__;
              app.mode = 'circular';
              app.cInputType = 'gb';
              app.files.c_gb = new File([gbText], 'HmmtDNA.gbk', { type: 'text/plain' });
              app.files.c_gff = null;
              app.files.c_fasta = null;
              await app.runAnalysis();
            }
            """,
            {"gbText": circular_gbk},
        )
        page.wait_for_function(
            """
            () => {
              const app = window.__GBDRAW_APP__;
              return (app.results && app.results.length > 0) || Boolean(app.errorLog);
            }
            """,
            timeout=120000,
        )
        circular_state = page.evaluate(
            """
            () => ({
              errorLog: window.__GBDRAW_APP__.errorLog,
              resultCount: window.__GBDRAW_APP__.results.length
            })
            """
        )
        if circular_state["errorLog"]:
            raise RuntimeError(f"Circular offline generation failed:\n{circular_state['errorLog']}")
        if circular_state["resultCount"] < 1:
            raise RuntimeError("Circular offline generation produced no SVG results.")

        download_dir = Path(tempfile.mkdtemp(prefix="gbdraw-offline-downloads-"))
        try:
            for method_name, extension in [
                ("downloadSVG", ".svg"),
                ("downloadPNG", ".png"),
                ("downloadPDF", ".pdf"),
            ]:
                with page.expect_download(timeout=120000) as download_info:
                    page.evaluate(f"() => window.__GBDRAW_APP__.{method_name}()")
                download = download_info.value
                target = download_dir / download.suggested_filename
                download.save_as(target)
                if target.suffix.lower() != extension:
                    raise RuntimeError(
                        f"{method_name} produced unexpected filename {download.suggested_filename}"
                    )
                if target.stat().st_size <= 0:
                    raise RuntimeError(f"{method_name} produced an empty {extension} file.")
        finally:
            shutil.rmtree(download_dir, ignore_errors=True)

        page.evaluate(
            """
            async ({ leftText, rightText }) => {
              const app = window.__GBDRAW_APP__;
              app.mode = 'linear';
              app.lInputType = 'gb';
              app.blastSource = 'losat';
              if (app.linearSeqs.length < 2) {
                app.addLinearSeq();
              }
              app.linearSeqs[0].gb = new File([leftText], 'MERS-CoV.gbk', { type: 'text/plain' });
              app.linearSeqs[0].gff = null;
              app.linearSeqs[0].fasta = null;
              app.linearSeqs[0].blast = null;
              app.linearSeqs[0].losat_filename = '';
              app.linearSeqs[1].gb = new File([rightText], 'SARS-CoV-1.gbk', { type: 'text/plain' });
              app.linearSeqs[1].gff = null;
              app.linearSeqs[1].fasta = null;
              await app.runAnalysis();
            }
            """,
            {"leftText": linear_left_gbk, "rightText": linear_right_gbk},
        )
        page.wait_for_function(
            """
            () => {
              const app = window.__GBDRAW_APP__;
              return Boolean(app.errorLog) || (
                app.results &&
                app.results.length > 0 &&
                app.losatCacheInfo &&
                app.losatCacheInfo.length > 0
              );
            }
            """,
            timeout=180000,
        )
        linear_state = page.evaluate(
            """
            () => ({
              errorLog: window.__GBDRAW_APP__.errorLog,
              resultCount: window.__GBDRAW_APP__.results.length,
              losatCacheEntries: window.__GBDRAW_APP__.losatCacheInfo.length
            })
            """
        )
        if linear_state["errorLog"]:
            raise RuntimeError(f"Linear LOSAT offline generation failed:\n{linear_state['errorLog']}")
        if linear_state["resultCount"] < 1 or linear_state["losatCacheEntries"] < 1:
            raise RuntimeError(
                f"Linear LOSAT offline generation did not populate expected outputs: {linear_state}"
            )

        if external_requests:
            raise RuntimeError(
                "External network requests were attempted during offline GUI startup:\n"
                + "\n".join(external_requests)
            )

        context.close()
        browser.close()


def inspect_wheel(wheel_path: Path) -> None:
    if not wheel_path.exists():
        raise FileNotFoundError(wheel_path)
    source_browser_wheel = _assert_source_browser_wheel_state()
    _assert_wheel_python_sources_match_source_tree(source_browser_wheel)
    _assert_outer_wheel_embeds_source_browser_wheel(wheel_path)
    _assert_outer_wheel_config_matches_source(wheel_path)
    with zipfile.ZipFile(wheel_path) as zf:
        names = set(zf.namelist())
    required = {
        "gbdraw/web/index.html",
        "gbdraw/web/open-source-notices.html",
        "gbdraw/web/js/app.js",
        "gbdraw/web/js/services/losat.js",
        "gbdraw/web/wasm/losat/losat.wasm",
        f"gbdraw/web/{_parse_wheel_name()}",
        "gbdraw/web/vendor/vue/vue.global.js",
        "gbdraw/web/vendor/tailwindcss/tailwindcss-play.js",
        "gbdraw/web/vendor/pyodide/v0.29.0/full/pyodide.js",
        "gbdraw/web/vendor/pyodide/v0.29.0/full/pyodide.asm.wasm",
        "gbdraw/web/vendor/jspdf/jspdf.umd.min.js",
        "gbdraw/web/vendor/svg2pdf.js/svg2pdf.umd.min.js",
        "gbdraw/web/vendor/dompurify/purify.min.js",
        "gbdraw/web/vendor/browser_wasi_shim/dist/index.js",
        "gbdraw/web/vendor/browser_wasi_shim/dist/wasi.js",
        "gbdraw/web/vendor/phosphor-icons/regular/style.css",
    }
    required.update(f"gbdraw/web/{path.as_posix()}" for path in REQUIRED_UI_FONT_FILES)
    required.update(f"gbdraw/web/{path.as_posix()}" for path in _parse_local_wheel_paths())
    missing = sorted(required - names)
    if missing:
        raise FileNotFoundError(
            "Wheel is missing required offline assets:\n" + "\n".join(missing)
        )


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Vendor and verify offline gbdraw GUI assets. "
            "Run `vendor-assets` for third-party browser assets, keep the local Pyodide dependency "
            "wheels under `gbdraw/web/vendor/pyodide-wheels/`, refresh the browser wheel with "
            "`python tools/sync_browser_wheel.py`, and keep `gbdraw/web/open-source-notices.html` committed."
        )
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    subparsers.add_parser("vendor-assets", help="Download and vendor third-party browser assets.")
    subparsers.add_parser("smoke-test", help="Run an offline GUI startup smoke test with Playwright.")

    inspect_parser = subparsers.add_parser("inspect-wheel", help="Inspect a built wheel for offline GUI assets.")
    inspect_parser.add_argument("wheel_path", type=Path)

    args = parser.parse_args()

    try:
        if args.command == "vendor-assets":
            vendor_assets()
        elif args.command == "smoke-test":
            smoke_test()
        elif args.command == "inspect-wheel":
            inspect_wheel(args.wheel_path)
        else:
            parser.error(f"Unknown command: {args.command}")
    except (FileNotFoundError, RuntimeError, urllib.error.URLError, subprocess.SubprocessError) as exc:
        print(str(exc), file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
