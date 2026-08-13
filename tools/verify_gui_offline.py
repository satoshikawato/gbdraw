#!/usr/bin/env python3
from __future__ import annotations

import argparse
import contextlib
import http.server
import io
import json
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
from runpy import run_path
from types import SimpleNamespace

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from gbdraw._web_assets import (
    ASSET_URLS,
    PYODIDE_CORE_FILES,
    PYODIDE_LOCAL_WHEELS,
    PYODIDE_RUNTIME_DIR,
    PYODIDE_RUNTIME_PACKAGE_WHEELS,
    REQUIRED_UI_FONT_FILES,
    UI_FONT_ASSETS,
    VENDOR_ROOT,
    WEB_ROOT,
)


def _load_build_support_module():
    return SimpleNamespace(**run_path(str(REPO_ROOT / "gbdraw" / "_build_support.py")))


BUILD_SUPPORT = _load_build_support_module()

def _required_tutorial_data_files() -> tuple[Path, ...]:
    manifest_path = WEB_ROOT / "tutorial-data" / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    declared = sorted(
        Path("tutorial-data") / metadata["relativePath"]
        for metadata in manifest["files"].values()
    )
    return (Path("tutorial-data/manifest.json"), *declared)


REQUIRED_TUTORIAL_DATA_FILES = _required_tutorial_data_files()


def _parse_local_wheel_paths() -> tuple[Path, ...]:
    return PYODIDE_LOCAL_WHEELS


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
            PYODIDE_RUNTIME_DIR,
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
        pyodide_target = PYODIDE_RUNTIME_DIR
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

        for asset_key, filename in PYODIDE_RUNTIME_PACKAGE_WHEELS.items():
            _download(ASSET_URLS[asset_key], PYODIDE_RUNTIME_DIR / filename)

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
        WEB_ROOT / "gallery" / "index.html",
        WEB_ROOT / "gallery" / "gallery.css",
        WEB_ROOT / "gallery" / "gallery.js",
        WEB_ROOT / "gallery" / "examples.json",
        WEB_ROOT / "assets" / "favicon.ico",
        WEB_ROOT / "assets" / "gbdraw-logo.svg",
        WEB_ROOT / "assets" / "gbdraw-logo-title.svg",
        WEB_ROOT / "assets" / "gbdraw-logo-title.png",
        WEB_ROOT / "js" / "app.js",
        WEB_ROOT / "js" / "services" / "losat.js",
        WEB_ROOT / "js" / "workers" / "losat-worker.js",
        WEB_ROOT / "js" / "workers" / "losat-threaded-worker.js",
        WEB_ROOT / "js" / "workers" / "losat-wasi-thread-worker.js",
        *(WEB_ROOT / path for path in REQUIRED_TUTORIAL_DATA_FILES),
        WEB_ROOT / "vendor" / "vue" / "vue.global.js",
        WEB_ROOT / "vendor" / "tailwindcss" / "tailwindcss-play.js",
        *(WEB_ROOT / path for path in REQUIRED_UI_FONT_FILES),
        *(PYODIDE_RUNTIME_DIR / filename for filename in PYODIDE_CORE_FILES),
        *(PYODIDE_RUNTIME_DIR / filename for filename in PYODIDE_RUNTIME_PACKAGE_WHEELS.values()),
        WEB_ROOT / "vendor" / "browser_wasi_shim" / "dist" / "index.js",
        WEB_ROOT / "vendor" / "browser_wasi_shim" / "dist" / "wasi.js",
        WEB_ROOT / "vendor" / "jspdf" / "jspdf.umd.min.js",
        WEB_ROOT / "vendor" / "svg2pdf.js" / "svg2pdf.umd.min.js",
        WEB_ROOT / "vendor" / "dompurify" / "purify.min.js",
        WEB_ROOT / "vendor" / "phosphor-icons" / "regular" / "style.css",
        WEB_ROOT / "wasm" / "losat" / "losat.wasm",
        WEB_ROOT / "wasm" / "losat" / "losat-threaded.wasm",
        browser_wheel_path,
        *(WEB_ROOT / path for path in _parse_local_wheel_paths()),
    ]
    missing = [path for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(
            "Missing required offline assets:\n" + "\n".join(str(path.relative_to(REPO_ROOT)) for path in missing)
        )


class QuietSimpleHTTPRequestHandler(http.server.SimpleHTTPRequestHandler):
    def end_headers(self) -> None:
        self.send_header("Cross-Origin-Opener-Policy", "same-origin")
        self.send_header("Cross-Origin-Embedder-Policy", "require-corp")
        self.send_header("Cross-Origin-Resource-Policy", "same-origin")
        self.send_header("Content-Security-Policy", "frame-ancestors 'none'")
        super().end_headers()

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


OFFLINE_GUI_BROWSER_CONTRACTS = (
    "offline-initialization",
    "generate",
    "helper-before-render",
    "palette-preview",
    "exports",
    "linear-losat",
)


DIAGRAM_WORKER_PROBE_INIT_SCRIPT = """
(() => {
  const activity = { constructions: 0, instances: [] };
  window.__GBDRAW_DIAGRAM_WORKER_ACTIVITY__ = activity;
  const NativeWorker = window.Worker;
  window.Worker = new Proxy(NativeWorker, {
    construct(target, args) {
      const worker = Reflect.construct(target, args, target);
      const url = String(args[0] || '');
      if (!url.includes('diagram-generation-worker.js')) return worker;

      activity.constructions += 1;
      const instance = {
        id: activity.constructions,
        url,
        initializations: 0,
        helpers: [],
        runs: [],
        settlements: [],
        terminated: false
      };
      activity.instances.push(instance);
      worker.addEventListener('message', (event) => {
        const message = event.data || {};
        if (!['init', 'helper', 'run'].includes(message.type)) return;
        instance.settlements.push({
          type: message.type,
          ok: message.ok === true,
          operation: String(message.operation || '')
        });
      });
      const nativePostMessage = worker.postMessage.bind(worker);
      worker.postMessage = (message, transfer) => {
        if (message?.type === 'init') {
          instance.initializations += 1;
        } else if (message?.type === 'helper') {
          instance.helpers.push({
            operation: String(message.operation || ''),
            requestId: String(message.requestId ?? '')
          });
        } else if (message?.type === 'run') {
          instance.runs.push({ requestId: String(message.requestId ?? '') });
        }
        if (transfer === undefined) return nativePostMessage(message);
        return nativePostMessage(message, transfer);
      };
      const nativeTerminate = worker.terminate.bind(worker);
      worker.terminate = () => {
        instance.terminated = true;
        return nativeTerminate();
      };
      return worker;
    }
  });
})();
"""


def _diagram_worker_activity(page) -> dict[str, object]:
    return page.evaluate(
        """
        () => {
          const activity = window.__GBDRAW_DIAGRAM_WORKER_ACTIVITY__ || {};
          const instances = Array.isArray(activity.instances) ? activity.instances : [];
          const settlements = instances.flatMap((instance) => instance.settlements || []);
          return {
            constructions: Number(activity.constructions || 0),
            initializations: instances.reduce(
              (total, instance) => total + Number(instance.initializations || 0), 0
            ),
            helpers: instances.reduce(
              (total, instance) => total + (instance.helpers?.length || 0), 0
            ),
            runs: instances.reduce(
              (total, instance) => total + (instance.runs?.length || 0), 0
            ),
            settledInitializations: settlements.filter(({ type }) => type === 'init').length,
            settledHelpers: settlements.filter(({ type }) => type === 'helper').length,
            settledRuns: settlements.filter(({ type }) => type === 'run').length,
            instances
          };
        }
        """
    )


def _assert_diagram_worker_activity(
    page,
    expected: dict[str, int],
    *,
    label: str,
) -> dict[str, object]:
    activity = _diagram_worker_activity(page)
    mismatches = {
        key: {"expected": value, "actual": activity.get(key)}
        for key, value in expected.items()
        if activity.get(key) != value
    }
    if mismatches:
        raise RuntimeError(
            f"{label} diagram Worker lifecycle mismatch:\n"
            f"{json.dumps({'mismatches': mismatches, 'activity': activity}, indent=2, sort_keys=True)}"
        )
    return activity


def _wait_for_semantic_state(
    page,
    predicate: str,
    *,
    timeout: int,
    label: str,
    snapshot: str,
) -> None:
    from playwright.sync_api import TimeoutError as PlaywrightTimeoutError

    try:
        page.wait_for_function(predicate, timeout=timeout)
    except PlaywrightTimeoutError as exc:
        try:
            state = page.evaluate(snapshot)
        except Exception as snapshot_exc:  # pragma: no cover - diagnostic fallback
            state = {"snapshotError": str(snapshot_exc)}
        raise RuntimeError(
            f"{label} timed out after {timeout} ms:\n"
            f"{json.dumps(state, indent=2, sort_keys=True)}"
        ) from exc


def _finish_browser_contract(
    context,
    browser,
    external_requests: list[str],
    contract: str,
) -> None:
    context.close()
    browser.close()
    if external_requests:
        raise RuntimeError(
            f"External network requests were attempted during {contract}:\n"
            + "\n".join(external_requests)
        )


def _verify_exports(page) -> None:
    download_dir = Path(tempfile.mkdtemp(prefix="gbdraw-offline-downloads-"))
    try:
        for method_name, filename_suffix in [
            ("downloadSVG", ".svg"),
            ("downloadPNG", ".png"),
            ("downloadPDF", ".pdf"),
            ("downloadInteractiveSVG", ".interactive.svg"),
        ]:
            with page.expect_download(timeout=120000) as download_info:
                page.evaluate(f"() => window.__GBDRAW_APP__.{method_name}()")
            download = download_info.value
            target = download_dir / download.suggested_filename
            download.save_as(target)
            if not target.name.lower().endswith(filename_suffix):
                raise RuntimeError(
                    f"{method_name} produced unexpected filename {download.suggested_filename}"
                )
            if target.stat().st_size <= 0:
                raise RuntimeError(
                    f"{method_name} produced an empty {filename_suffix} file."
                )

            content = target.read_bytes()
            if method_name == "downloadPNG" and not content.startswith(b"\x89PNG\r\n\x1a\n"):
                raise RuntimeError("downloadPNG did not produce a PNG payload.")
            if method_name == "downloadPDF" and not content.startswith(b"%PDF-"):
                raise RuntimeError("downloadPDF did not produce a PDF payload.")
            if method_name in {"downloadSVG", "downloadInteractiveSVG"}:
                svg_text = content.decode("utf-8")
                if "<svg" not in svg_text:
                    raise RuntimeError(f"{method_name} did not produce an SVG payload.")
                metadata_marker = 'id="gbdraw-interactive-feature-metadata"'
                if method_name == "downloadSVG" and metadata_marker in svg_text:
                    raise RuntimeError("downloadSVG unexpectedly included interactive metadata.")
                if method_name == "downloadInteractiveSVG":
                    required_markers = (
                        metadata_marker,
                        'id="gbdraw-interactive-feature-style"',
                        'id="gbdraw-interactive-feature-script"',
                        'data-gbdraw-interactive-svg="true"',
                    )
                    missing = [value for value in required_markers if value not in svg_text]
                    if missing:
                        raise RuntimeError(
                            "downloadInteractiveSVG omitted standalone assets or metadata: "
                            + ", ".join(missing)
                        )
    finally:
        shutil.rmtree(download_dir, ignore_errors=True)


def _verify_linear_losat(page, left_gbk: str, right_gbk: str) -> None:
    setup_state = page.evaluate(
        """
        ({ leftText, rightText }) => {
          const app = window.__GBDRAW_APP__;
          app.mode = 'linear';
          app.lInputType = 'gb';
          if (app.linearSeqs.length < 2) {
            app.addLinearSeq();
          }
          app.setLinearSeqPrimaryFile(
            0,
            'gb',
            new File([leftText], 'MERS-CoV.gbk', { type: 'text/plain' })
          );
          app.setLinearSeqPrimaryFile(
            1,
            'gb',
            new File([rightText], 'SARS-CoV-1.gbk', { type: 'text/plain' })
          );
          app.setLinearComparisonGlobalAction('losat');
          if (!app.setLinearComparisonLosatMode('blastn')) {
            throw new Error('Could not select the current Linear LOSAT nucleotide mode.');
          }
          const run = {
            status: 'running',
            error: '',
            startedResultCount: Array.isArray(app.results) ? app.results.length : 0,
            startedCacheCount: Array.isArray(app.losatCacheInfo)
              ? app.losatCacheInfo.length
              : 0
          };
          window.__GBDRAW_OFFLINE_CONTRACT_RUN__ = run;
          void (async () => {
            try {
              await app.runAnalysis();
              run.status = 'fulfilled';
            } catch (error) {
              run.status = 'rejected';
              run.error = String(error?.stack || error?.message || error);
            }
          })();
          return {
            comparisonMode: app.linearComparisonPlan?.mode || '',
            comparisonDefaultSource: app.linearComparisonPlan?.defaultSource || '',
            comparisonGlobalAction: app.linearComparisonGlobalAction,
            hasActiveLosatIntent: app.hasActiveLinearLosatIntent === true,
            losatProgram: app.losatProgram
          };
        }
        """,
        {"leftText": left_gbk, "rightText": right_gbk},
    )
    if setup_state != {
        "comparisonMode": "adjacent",
        "comparisonDefaultSource": "losat",
        "comparisonGlobalAction": "losat",
        "hasActiveLosatIntent": True,
        "losatProgram": "blastn",
    }:
        raise RuntimeError(
            "Linear LOSAT comparison intent was not activated before generation: "
            f"{setup_state}"
        )
    _wait_for_semantic_state(
        page,
        """
        () => {
          const app = window.__GBDRAW_APP__;
          const run = window.__GBDRAW_OFFLINE_CONTRACT_RUN__;
          if (!app || !run) return false;
          return Boolean(app.errorLog) || (
            run.status === 'rejected'
          ) || (
            run.status === 'fulfilled' &&
            Array.isArray(app.results) &&
            app.results.length > run.startedResultCount &&
            Array.isArray(app.losatCacheInfo) &&
            app.losatCacheInfo.length > run.startedCacheCount
          );
        }
        """,
        timeout=120000,
        label="Linear LOSAT generation and cache readiness",
        snapshot="""
        () => {
          const app = window.__GBDRAW_APP__;
          if (!app) return { appMounted: false };
          return {
            appMounted: true,
            run: window.__GBDRAW_OFFLINE_CONTRACT_RUN__ || null,
            mode: app.mode,
            inputType: app.lInputType,
            errorLog: app.errorLog,
            resultCount: Array.isArray(app.results) ? app.results.length : 0,
            losatCacheEntries: Array.isArray(app.losatCacheInfo)
              ? app.losatCacheInfo.length
              : 0,
            comparisonMode: app.linearComparisonPlan?.mode || '',
            comparisonDefaultSource: app.linearComparisonPlan?.defaultSource || '',
            comparisonGlobalAction: app.linearComparisonGlobalAction,
            hasActiveLosatIntent: app.hasActiveLinearLosatIntent === true,
            comparisonEdgeCount: app.linearComparisonResolution?.edges?.length || 0,
            comparisonErrors: app.linearComparisonResolution?.errors || [],
            losatProgram: app.losatProgram,
            losatThreadingStatus: app.losatThreadingStatus
          };
        }
        """,
    )
    linear_state = page.evaluate(
        """
        () => ({
          errorLog: window.__GBDRAW_APP__.errorLog,
          runStatus: window.__GBDRAW_OFFLINE_CONTRACT_RUN__?.status || '',
          runError: window.__GBDRAW_OFFLINE_CONTRACT_RUN__?.error || '',
          resultCount: window.__GBDRAW_APP__.results.length,
          losatCacheEntries: window.__GBDRAW_APP__.losatCacheInfo.length,
          comparisonMode: window.__GBDRAW_APP__.linearComparisonPlan?.mode || '',
          comparisonDefaultSource:
            window.__GBDRAW_APP__.linearComparisonPlan?.defaultSource || '',
          hasActiveLosatIntent:
            window.__GBDRAW_APP__.hasActiveLinearLosatIntent === true
        })
        """
    )
    if linear_state["errorLog"]:
        raise RuntimeError(
            f"Linear LOSAT offline generation failed:\n{linear_state['errorLog']}"
        )
    if linear_state["runStatus"] != "fulfilled":
        raise RuntimeError(
            "Linear LOSAT generation did not settle successfully: "
            f"{linear_state}"
        )
    if (
        linear_state["comparisonMode"] != "adjacent"
        or linear_state["comparisonDefaultSource"] != "losat"
        or not linear_state["hasActiveLosatIntent"]
    ):
        raise RuntimeError(
            "Linear LOSAT comparison intent changed during generation: "
            f"{linear_state}"
        )
    if linear_state["resultCount"] < 1 or linear_state["losatCacheEntries"] < 1:
        raise RuntimeError(
            "Linear LOSAT offline generation did not populate expected outputs: "
            f"{linear_state}"
        )


def _run_browser_contract(contract: str) -> None:
    if contract not in OFFLINE_GUI_BROWSER_CONTRACTS:
        raise ValueError(f"Unknown offline GUI browser contract: {contract}")
    _assert_packaged_assets()
    _ensure_playwright_available()

    from playwright.sync_api import sync_playwright

    with _serve_web_root() as base_url, sync_playwright() as p:
        try:
            browser = p.chromium.launch(headless=True)
        except Exception as exc:  # pragma: no cover - environment dependent
            raise RuntimeError(
                "Could not launch Playwright Chromium for offline GUI verification. "
                f"Underlying error: {exc}"
            ) from exc
        context = browser.new_context(accept_downloads=True)
        context.add_init_script(script=DIAGRAM_WORKER_PROBE_INIT_SCRIPT)
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
        isolation_state = page.evaluate(
            "() => ({ isolated: crossOriginIsolated === true, sharedArrayBuffer: typeof SharedArrayBuffer === 'function' })"
        )
        if not isolation_state["isolated"] or not isolation_state["sharedArrayBuffer"]:
            raise RuntimeError(
                "GUI local server did not enable browser isolation for threaded LOSAT: "
                f"{isolation_state}"
            )
        startup_snapshot = """
            () => {
              const app = window.__GBDRAW_APP__;
              if (!app) return { appMounted: false };
              const runtimeFields = [
                'pyodide', 'pyodideReady', 'pyodideLoading', 'pyodideError', 'pyodideStatus'
              ];
              return {
                appMounted: true,
                mainLoaderPresent: typeof window.loadPyodide === 'function',
                mainRuntimeFields: runtimeFields.filter(
                  (field) => Object.prototype.hasOwnProperty.call(app, field)
                ),
                paletteDefinitionCount: Object.keys(app.paletteDefinitions || {}).length
              };
            }
        """
        _wait_for_semantic_state(
            page,
            """
            () => {
              const app = window.__GBDRAW_APP__;
              if (!app) return false;
              return Object.keys(app.paletteDefinitions || {}).length > 0;
            }
            """,
            timeout=120000,
            label="Offline GUI startup",
            snapshot=startup_snapshot,
        )
        startup_state = page.evaluate(startup_snapshot)
        if startup_state["mainLoaderPresent"] or startup_state["mainRuntimeFields"]:
            raise RuntimeError(
                "GUI exposed a forbidden main-thread Python runtime boundary: "
                f"{startup_state}"
            )
        if startup_state["paletteDefinitionCount"] == 0:
            raise RuntimeError(
                f"GUI palettes failed to load offline: {startup_state}"
            )
        _assert_diagram_worker_activity(
            page,
            {
                "constructions": 0,
                "initializations": 0,
                "helpers": 0,
                "runs": 0,
            },
            label="Offline app shell",
        )

        if contract == "offline-initialization":
            _finish_browser_contract(context, browser, external_requests, contract)
            return

        if contract == "linear-losat":
            linear_left_gbk = (
                REPO_ROOT / "tests" / "test_inputs" / "MERS-CoV.gbk"
            ).read_text(encoding="utf-8")
            linear_right_gbk = (
                REPO_ROOT / "tests" / "test_inputs" / "SARS-CoV-1.gbk"
            ).read_text(encoding="utf-8")
            _verify_linear_losat(page, linear_left_gbk, linear_right_gbk)
            _assert_diagram_worker_activity(
                page,
                {"constructions": 1, "initializations": 1, "runs": 1},
                label="Offline Linear LOSAT generation",
            )
            _finish_browser_contract(context, browser, external_requests, contract)
            return

        if contract == "helper-before-render":
            helper_result = page.evaluate(
                """
                async () => {
                  const {
                    DIAGRAM_HELPER_OPERATIONS,
                    runDiagramHelperOperation
                  } = await import('./js/services/diagram-generation.js');
                  const response = await runDiagramHelperOperation(
                    DIAGRAM_HELPER_OPERATIONS.MEASURE_LEGEND_TEXT,
                    { caption: 'offline helper probe', fontFamily: 'Arial', fontSize: 14 }
                  );
                  return response?.result || null;
                }
                """
            )
            if not isinstance(helper_result, dict) or not isinstance(
                helper_result.get("width"), (int, float)
            ):
                raise RuntimeError(
                    "Offline helper operation did not settle with a measured width: "
                    f"{helper_result!r}"
                )
            _assert_diagram_worker_activity(
                page,
                {
                    "constructions": 1,
                    "initializations": 1,
                    "helpers": 1,
                    "runs": 0,
                    "settledInitializations": 1,
                    "settledHelpers": 1,
                    "settledRuns": 0,
                },
                label="Offline helper-only operation",
            )

        circular_gbk = (
            REPO_ROOT / "tests" / "test_inputs" / "HmmtDNA.gbk"
        ).read_text(encoding="utf-8")
        page.evaluate(
            """
            ({ gbText }) => {
              const app = window.__GBDRAW_APP__;
              app.mode = 'circular';
              app.cInputType = 'gb';
              app.files.c_gb = new File([gbText], 'HmmtDNA.gbk', { type: 'text/plain' });
              app.files.c_gff = null;
              app.files.c_fasta = null;
              const run = { status: 'running', error: '' };
              window.__GBDRAW_OFFLINE_CONTRACT_RUN__ = run;
              void (async () => {
                try {
                  await app.runAnalysis();
                  run.status = 'fulfilled';
                } catch (error) {
                  run.status = 'rejected';
                  run.error = String(error?.stack || error?.message || error);
                }
              })();
            }
            """,
            {"gbText": circular_gbk},
        )
        _wait_for_semantic_state(
            page,
            """
            () => {
              const app = window.__GBDRAW_APP__;
              const run = window.__GBDRAW_OFFLINE_CONTRACT_RUN__;
              if (!app || !run) return false;
              return run.status === 'rejected' || Boolean(app.errorLog) || (
                run.status === 'fulfilled' &&
                Array.isArray(app.results) &&
                app.results.length > 0
              );
            }
            """,
            timeout=120000,
            label="Circular diagram generation",
            snapshot="""
            () => {
              const app = window.__GBDRAW_APP__;
              return {
                run: window.__GBDRAW_OFFLINE_CONTRACT_RUN__ || null,
                errorLog: app?.errorLog || '',
                resultCount: Array.isArray(app?.results) ? app.results.length : 0
              };
            }
            """,
        )
        circular_state = page.evaluate(
            """
            () => ({
              errorLog: window.__GBDRAW_APP__.errorLog,
              runStatus: window.__GBDRAW_OFFLINE_CONTRACT_RUN__?.status || '',
              runError: window.__GBDRAW_OFFLINE_CONTRACT_RUN__?.error || '',
              resultCount: window.__GBDRAW_APP__.results.length
            })
            """
        )
        if circular_state["errorLog"]:
            raise RuntimeError(f"Circular offline generation failed:\n{circular_state['errorLog']}")
        if circular_state["runStatus"] != "fulfilled":
            raise RuntimeError(
                "Circular offline generation did not settle successfully: "
                f"{circular_state}"
            )
        if circular_state["resultCount"] < 1:
            raise RuntimeError("Circular offline generation produced no SVG results.")
        activity = _assert_diagram_worker_activity(
            page,
            {"constructions": 1, "initializations": 1, "runs": 1, "settledRuns": 1},
            label="Offline Circular generation",
        )
        if contract == "helper-before-render" and activity["helpers"] < 1:
            raise RuntimeError(
                "The helper-before-render scenario lost its settled helper operation: "
                f"{activity}"
            )
        if contract in {"generate", "helper-before-render"}:
            _finish_browser_contract(context, browser, external_requests, contract)
            return

        if contract == "exports":
            _verify_exports(page)
            _finish_browser_contract(context, browser, external_requests, contract)
            return

        _wait_for_semantic_state(
            page,
            """
            () => {
              const app = window.__GBDRAW_APP__;
              if (!app) return false;
              return Boolean(app.featureExtractionError) ||
                (Array.isArray(app.extractedFeatures) && app.extractedFeatures.length > 0);
            }
            """,
            timeout=120000,
            label="Circular feature extraction",
            snapshot="""
            () => {
              const app = window.__GBDRAW_APP__;
              return {
                featureExtractionError: app?.featureExtractionError || '',
                extractedFeatureCount: Array.isArray(app?.extractedFeatures)
                  ? app.extractedFeatures.length
                  : 0,
                resultCount: Array.isArray(app?.results) ? app.results.length : 0
              };
            }
            """,
        )
        feature_extraction_state = page.evaluate(
            """
            () => ({
              featureExtractionError: window.__GBDRAW_APP__.featureExtractionError,
              extractedFeatureCount: Array.isArray(window.__GBDRAW_APP__.extractedFeatures)
                ? window.__GBDRAW_APP__.extractedFeatures.length
                : 0
            })
            """
        )
        if feature_extraction_state["featureExtractionError"]:
            raise RuntimeError(
                "Circular feature extraction failed offline:\n"
                f"{feature_extraction_state['featureExtractionError']}"
            )
        if feature_extraction_state["extractedFeatureCount"] < 1:
            raise RuntimeError("Circular feature extraction produced no editable features.")

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
              const waitForExtractedFeatures = async () => {
                for (let attempt = 0; attempt < 240; attempt += 1) {
                  if (app.featureExtractionError) return false;
                  if (Array.isArray(app.extractedFeatures) && app.extractedFeatures.length > 0) {
                    return true;
                  }
                  await new Promise((resolve) => setTimeout(resolve, 50));
                }
                return false;
              };
              const waitForFeatureFillByType = async (featureType, preferredSvgId, expectedFill) => {
                let lastFill = findFeatureFillByType(featureType, preferredSvgId);
                for (let attempt = 0; attempt < 20; attempt += 1) {
                  if (lastFill === expectedFill) return lastFill;
                  await new Promise((resolve) => requestAnimationFrame(() => resolve()));
                  lastFill = findFeatureFillByType(featureType, preferredSvgId);
                }
                return lastFill;
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
              if (!(await waitForExtractedFeatures())) {
                return { error: 'Feature extraction did not complete after regenerating the circular diagram.' };
              }
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
              const afterImmediateFill = await waitForFeatureFillByType(
                chosenFeature.type,
                chosenFeature.svgId,
                immediateColor
              );

              return {
                featureType: chosenFeature.type,
                featureId: chosenFeature.id,
                svgId: chosenFeature.svgId,
                svgElementTag: chosenFeature.svgId
                  ? app.svgContainer?.querySelector('svg')?.querySelector(`#${CSS.escape(chosenFeature.svgId)}`)?.tagName || ''
                  : '',
                extractedHasSvgAfterImmediate: (app.extractedFeatures || []).some(
                  (feat) => feat?.svg_id === chosenFeature.svgId
                ),
                extractedFeatureCountAfterImmediate: Array.isArray(app.extractedFeatures)
                  ? app.extractedFeatures.length
                  : 0,
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
                appliedAfterImmediate: String(app.appliedPaletteColors?.[chosenFeature.type] || ''),
                pendingAfterImmediate: String(app.pendingPaletteName || ''),
                overrideAfterImmediate: app.featureColorOverrides?.[chosenFeature.id] || null,
                specificRuleCountAfterImmediate: Array.isArray(app.manualSpecificRules)
                  ? app.manualSpecificRules.length
                  : 0,
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
                "Palette instant preview did not update the SVG immediately after being re-enabled: "
                f"{deferred_palette_state}"
            )

        _finish_browser_contract(context, browser, external_requests, contract)
        return


def smoke_test(contract: str = "all") -> None:
    contracts = (
        OFFLINE_GUI_BROWSER_CONTRACTS
        if contract == "all"
        else (contract,)
    )
    for browser_contract in contracts:
        _run_browser_contract(browser_contract)


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
        "gbdraw/web/assets/favicon.ico",
        "gbdraw/web/assets/gbdraw-logo.svg",
        "gbdraw/web/assets/gbdraw-logo-title.svg",
        "gbdraw/web/assets/gbdraw-logo-title.png",
        "gbdraw/web/js/app.js",
        "gbdraw/web/js/services/losat.js",
        "gbdraw/web/js/workers/losat-worker.js",
        "gbdraw/web/js/workers/losat-threaded-worker.js",
        "gbdraw/web/js/workers/losat-wasi-thread-worker.js",
        "gbdraw/web/wasm/losat/losat.wasm",
        "gbdraw/web/wasm/losat/losat-threaded.wasm",
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
    required.update(
        f"gbdraw/web/{path.as_posix()}" for path in REQUIRED_TUTORIAL_DATA_FILES
    )
    required.update(f"gbdraw/web/{path.as_posix()}" for path in _parse_local_wheel_paths())
    required.update(
        f"gbdraw/web/vendor/pyodide/v0.29.0/full/{filename}"
        for filename in PYODIDE_CORE_FILES
    )
    required.update(
        f"gbdraw/web/vendor/pyodide/v0.29.0/full/{filename}"
        for filename in PYODIDE_RUNTIME_PACKAGE_WHEELS.values()
    )
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
    smoke_parser = subparsers.add_parser(
        "smoke-test",
        help="Run one or all independent offline GUI browser contracts with Playwright.",
    )
    smoke_parser.add_argument(
        "--contract",
        choices=("all", *OFFLINE_GUI_BROWSER_CONTRACTS),
        default="all",
        help="Select one browser contract; the default runs every contract in a fresh context.",
    )

    inspect_parser = subparsers.add_parser("inspect-wheel", help="Inspect a built wheel for offline GUI assets.")
    inspect_parser.add_argument("wheel_path", type=Path)

    args = parser.parse_args()

    try:
        if args.command == "vendor-assets":
            vendor_assets()
        elif args.command == "check-assets":
            check_assets()
        elif args.command == "smoke-test":
            smoke_test(args.contract)
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
