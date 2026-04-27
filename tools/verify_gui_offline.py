#!/usr/bin/env python3
from __future__ import annotations

import argparse
import contextlib
import http.server
import io
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
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
WEB_ROOT = REPO_ROOT / "gbdraw" / "web"
VENDOR_ROOT = WEB_ROOT / "vendor"

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


def _load_build_support_module():
    build_support_path = REPO_ROOT / "gbdraw" / "_build_support.py"
    spec = spec_from_file_location("gbdraw_build_support", build_support_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not load build support module from {build_support_path}")

    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


BUILD_SUPPORT = _load_build_support_module()


def _parse_local_wheel_paths() -> tuple[Path, ...]:
    config_text = (WEB_ROOT / "js" / "config.js").read_text(encoding="utf-8")
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
    return BUILD_SUPPORT.read_browser_wheel_name_from_config()


def check_assets() -> None:
    BUILD_SUPPORT.validate_browser_wheel_prepared()
    _assert_packaged_assets()


def _assert_packaged_assets() -> None:
    browser_wheel_path = BUILD_SUPPORT.validate_browser_wheel_prepared()
    required = [
        WEB_ROOT / "index.html",
        WEB_ROOT / "open-source-notices.html",
        WEB_ROOT / "js" / "app.js",
        WEB_ROOT / "js" / "services" / "losat.js",
        WEB_ROOT / "js" / "workers" / "losat-worker.js",
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
        browser_wheel_path,
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

        deferred_palette_state = page.evaluate(
            """
            async () => {
              const app = window.__GBDRAW_APP__;
              const svg = app.svgContainer?.querySelector('svg');
              if (!svg) {
                return { error: 'SVG preview not found after circular generation.' };
              }

              const getFeatureFill = (svgId) => {
                if (!svgId) return null;
                const element = svg.querySelector(`#${CSS.escape(svgId)}`);
                return element ? element.getAttribute('fill') : null;
              };
              const findFeatureFillByType = (featureType, preferredSvgId = '') => {
                const currentSvg = app.svgContainer?.querySelector('svg');
                if (!currentSvg) return null;
                if (preferredSvgId) {
                  const preferred = currentSvg.querySelector(`#${CSS.escape(preferredSvgId)}`);
                  if (preferred) {
                    return preferred.getAttribute('fill');
                  }
                }
                for (const feat of app.extractedFeatures || []) {
                  if (feat?.type !== featureType || !feat?.svg_id) continue;
                  const element = currentSvg.querySelector(`#${CSS.escape(feat.svg_id)}`);
                  if (element) return element.getAttribute('fill');
                }
                return null;
              };

              const pickFeatureForPaletteCheck = () => {
                for (const feat of app.extractedFeatures || []) {
                  if (!feat?.svg_id || !feat?.type) continue;
                  const fill = getFeatureFill(feat.svg_id);
                  if (!fill) continue;
                  const baseColor = app.currentColors?.[feat.type];
                  if (!baseColor) continue;
                  return { svgId: feat.svg_id, type: feat.type, beforeFill: fill };
                }
                return null;
              };

              const chosenFeature = pickFeatureForPaletteCheck();
              if (!chosenFeature) {
                return { error: 'Could not find a rendered feature for palette verification.' };
              }

              const originalPalette = String(app.selectedPalette || 'default');
              const originalDraftColor = String(app.currentColors?.[chosenFeature.type] || '');
              app.paletteInstantPreviewEnabled = false;

              let pendingPalette = '';
              let pendingDraftColor = '';
              for (const paletteName of app.paletteNames || []) {
                if (paletteName === originalPalette) continue;
                app.selectedPalette = paletteName;
                app.updatePalette();
                const candidateColor = String(app.currentColors?.[chosenFeature.type] || '');
                if (!candidateColor || candidateColor === originalDraftColor) continue;
                pendingPalette = paletteName;
                pendingDraftColor = candidateColor;
                break;
              }

              if (!pendingPalette) {
                return {
                  error: `Could not find a palette with a different ${chosenFeature.type} color for deferred-preview verification.`
                };
              }

              await new Promise((resolve) => requestAnimationFrame(() => resolve()));
              const afterPaletteFill = getFeatureFill(chosenFeature.svgId);
              const pendingName = String(app.pendingPaletteName || '');
              const appliedColorAfterPalette = String(app.appliedPaletteColors?.[chosenFeature.type] || '');
              const pendingColorAfterPalette = String(app.pendingPaletteColors?.[chosenFeature.type] || '');

              const manualDraftColor = [pendingDraftColor, originalDraftColor].includes('#123456') ? '#654321' : '#123456';
              app.currentColors[chosenFeature.type] = manualDraftColor;
              await new Promise((resolve) => requestAnimationFrame(() => resolve()));
              const afterDraftEditFill = getFeatureFill(chosenFeature.svgId);

              await app.runAnalysis();
              const afterGenerateFill = findFeatureFillByType(chosenFeature.type, chosenFeature.svgId);
              const pendingAfterGenerate = String(app.pendingPaletteName || '');
              const appliedAfterGenerate = String(app.appliedPaletteColors?.[chosenFeature.type] || '');

              let immediatePalette = '';
              let immediateColor = '';
              for (const paletteName of app.paletteNames || []) {
                if (paletteName === String(app.selectedPalette || '')) continue;
                app.selectedPalette = paletteName;
                app.paletteInstantPreviewEnabled = true;
                app.updatePalette();
                const candidateColor = String(app.currentColors?.[chosenFeature.type] || '');
                if (!candidateColor || candidateColor === appliedAfterGenerate) continue;
                immediatePalette = paletteName;
                immediateColor = candidateColor;
                break;
              }

              if (!immediatePalette) {
                return {
                  error: `Could not find a second palette with a different ${chosenFeature.type} color for instant-preview verification.`
                };
              }

              await new Promise((resolve) => requestAnimationFrame(() => resolve()));
              const afterImmediateFill = findFeatureFillByType(chosenFeature.type, chosenFeature.svgId);

              return {
                featureType: chosenFeature.type,
                beforeFill: chosenFeature.beforeFill,
                afterPaletteFill,
                afterDraftEditFill,
                afterGenerateFill,
                afterImmediateFill,
                pendingName,
                pendingAfterGenerate,
                appliedColorAfterPalette,
                pendingColorAfterPalette,
                appliedAfterGenerate,
                manualDraftColor,
                pendingPalette,
                pendingDraftColor,
                immediatePalette,
                immediateColor
              };
            }
            """
        )
        if deferred_palette_state.get("error"):
            raise RuntimeError(deferred_palette_state["error"])
        if deferred_palette_state["afterPaletteFill"] != deferred_palette_state["beforeFill"]:
            raise RuntimeError(
                "Palette change updated the SVG even though palette instant preview was disabled."
            )
        if deferred_palette_state["afterDraftEditFill"] != deferred_palette_state["beforeFill"]:
            raise RuntimeError(
                "Manual -d edit updated the SVG before Generate Diagram while a palette draft was pending."
            )
        if deferred_palette_state["pendingName"] != deferred_palette_state["pendingPalette"]:
            raise RuntimeError(
                "Pending palette state did not record the deferred palette selection."
            )
        if deferred_palette_state["pendingColorAfterPalette"] != deferred_palette_state["pendingDraftColor"]:
            raise RuntimeError(
                "Pending palette colors were not updated after selecting a deferred palette."
            )
        if deferred_palette_state["afterGenerateFill"] != deferred_palette_state["manualDraftColor"]:
            raise RuntimeError(
                "Generate Diagram did not apply the deferred palette draft color to the SVG."
            )
        if deferred_palette_state["pendingAfterGenerate"]:
            raise RuntimeError("Pending palette state was not cleared after Generate Diagram.")
        if deferred_palette_state["appliedAfterGenerate"] != deferred_palette_state["manualDraftColor"]:
            raise RuntimeError("Applied palette colors were not promoted after Generate Diagram.")
        if deferred_palette_state["afterImmediateFill"] != deferred_palette_state["immediateColor"]:
            raise RuntimeError(
                "Palette instant preview did not update the SVG immediately after being re-enabled."
            )

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
    with zipfile.ZipFile(wheel_path) as zf:
        names = set(zf.namelist())
    expected_browser_wheel = f"gbdraw/web/{_parse_wheel_name()}"
    browser_wheels = sorted(
        name for name in names if name.startswith("gbdraw/web/gbdraw-") and name.endswith(".whl")
    )
    if browser_wheels != [expected_browser_wheel]:
        raise RuntimeError(
            "Wheel contains unexpected embedded browser wheels:\n" + "\n".join(browser_wheels)
        )
    required = {
        "gbdraw/web/index.html",
        "gbdraw/web/open-source-notices.html",
        "gbdraw/web/js/app.js",
        "gbdraw/web/js/services/losat.js",
        "gbdraw/web/wasm/losat/losat.wasm",
        expected_browser_wheel,
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
    assert_embedded_browser_wheel_is_not_recursive(wheel_path)


def assert_browser_wheel_is_not_recursive(wheel_path: Path) -> None:
    browser_wheel_member = f"gbdraw/web/{_parse_wheel_name()}"
    with zipfile.ZipFile(wheel_path) as browser_wheel:
        nested_names = set(browser_wheel.namelist())
    _assert_browser_wheel_names_not_recursive(nested_names, browser_wheel_member)


def assert_embedded_browser_wheel_is_not_recursive(wheel_path: Path) -> None:
    browser_wheel_member = f"gbdraw/web/{_parse_wheel_name()}"
    with zipfile.ZipFile(wheel_path) as zf:
        try:
            browser_wheel_bytes = zf.read(browser_wheel_member)
        except KeyError as exc:
            raise FileNotFoundError(
                f"Wheel is missing embedded browser wheel: {browser_wheel_member}"
            ) from exc

    with zipfile.ZipFile(io.BytesIO(browser_wheel_bytes)) as browser_wheel:
        nested_names = set(browser_wheel.namelist())
    _assert_browser_wheel_names_not_recursive(nested_names, browser_wheel_member)


def _assert_browser_wheel_names_not_recursive(nested_names: set[str], browser_wheel_member: str) -> None:
    if browser_wheel_member in nested_names:
        raise RuntimeError(
            f"Browser wheel recursively contains itself: {browser_wheel_member}"
        )


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Vendor and verify offline gbdraw GUI assets. "
            "Run `vendor-assets` for third-party browser assets, keep the local Pyodide dependency "
            "wheels under `gbdraw/web/vendor/pyodide-wheels/`, run `python tools/prepare_browser_wheel.py` "
            "for the local browser wheel, and keep `gbdraw/web/open-source-notices.html` committed."
        )
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    subparsers.add_parser("vendor-assets", help="Download and vendor third-party browser assets.")
    subparsers.add_parser("check-assets", help="Validate the prepared browser wheel and required offline assets.")
    subparsers.add_parser("smoke-test", help="Run an offline GUI startup smoke test with Playwright.")

    inspect_parser = subparsers.add_parser("inspect-wheel", help="Inspect a built wheel for offline GUI assets.")
    inspect_parser.add_argument("wheel_path", type=Path)

    args = parser.parse_args()

    try:
        if args.command == "vendor-assets":
            vendor_assets()
        elif args.command == "check-assets":
            check_assets()
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
