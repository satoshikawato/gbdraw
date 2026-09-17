"""Run the existing real-helper/replay acceptance offline at two viewports.

No analysis or render implementation is replaced here. The reused acceptance
also injects cancellation/render failure and prohibits fresh LOSAT search while
reusing its verified saved raw fixture.
"""
from __future__ import annotations

import argparse
import json
import hashlib
import shutil
import time
from pathlib import Path
from unittest.mock import patch
from urllib.parse import urlparse

from playwright.sync_api import Browser

from tests import run_losat_cache_browser_acceptance as acceptance


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    wheels = list((root / "gbdraw/web").glob("gbdraw-*-py3-none-any.whl"))
    if len(wheels) != 1:
        raise RuntimeError("Acceptance requires one prepared browser wheel")
    artifact_paths = [wheels[0],
                      *(root / name for name in (
                          "gbdraw/analysis/ortholog_paths.py", "gbdraw/analysis/protein_colinearity.py",
                          "gbdraw/session_request_codec.py", "gbdraw/api/session_compat.py",
                          "gbdraw/web/js/app/run-analysis.js", "gbdraw/web/js/services/session-request.js"))]
    artifact_hashes = {str(path.relative_to(root)): hashlib.sha256(path.read_bytes()).hexdigest()
                       for path in artifact_paths}
    runs = []
    original = Browser.new_context
    for width, height in ((1280, 720), (390, 844)):
        run = {"viewport": {"width": width, "height": height}, "paths": set(),
               "blockedExternal": [], "pageErrors": []}

        def new_context(browser, *positional, **kwargs):
            run["browser"] = browser.version
            context = original(browser, *positional, **kwargs, viewport=run["viewport"])

            def route(request_route):
                url = urlparse(request_route.request.url)
                if url.hostname not in ("127.0.0.1", "localhost"):
                    run["blockedExternal"].append(request_route.request.url)
                    request_route.abort()
                else:
                    run["paths"].add(url.path)
                    request_route.continue_()

            context.route("**/*", route)
            context.on("page", lambda page: page.on("pageerror", lambda error: run["pageErrors"].append(str(error))))
            return context

        operations = []
        original_generate, original_save = acceptance._generate, acceptance._save_session
        artifact_dir = root / '.venv/s06-browser-sessions' / args.output.stem / str(width)
        artifact_dir.mkdir(parents=True, exist_ok=True)
        def generate(page):
            started = time.perf_counter()
            result = original_generate(page)
            operations.append({"operation": "Generate", "elapsedSeconds": time.perf_counter() - started,
                               "result": result.get("result"), "telemetry": result.get("telemetry")})
            return result
        def save(page, checks):
            started = time.perf_counter()
            result = original_save(page, checks)
            elapsed = time.perf_counter() - started
            data = result.read_bytes()
            suffix = '.json.gz' if data.startswith(b'\x1f\x8b') else '.json'
            target = artifact_dir / (f'session-{len(operations)}' + suffix)
            shutil.copyfile(result, target)
            operations.append({"operation": "Save", "elapsedSeconds": elapsed, "bytes": len(data),
                               "sha256": hashlib.sha256(data).hexdigest(), "artifact": str(target)})
            return result
        with patch.object(Browser, "new_context", new_context), \
             patch.object(acceptance, "_generate", generate), \
             patch.object(acceptance, "_save_session", save):
            code = acceptance._run_python_adapter()
        run["operationObservations"] = operations
        run["paths"] = sorted(run["paths"])
        run["exitCode"] = code
        runs.append(run)
        if code or run["blockedExternal"] or run["pageErrors"]:
            break
    if any(hashlib.sha256(path.read_bytes()).hexdigest() != artifact_hashes[str(path.relative_to(root))]
           for path in artifact_paths):
        raise RuntimeError("Browser artifact changed during acceptance")
    report = {"artifact": "source SPA with generated local browser wheel",
              "artifactSha256": artifact_hashes,
              "network": "fresh contexts; external requests blocked before navigation",
              "acceptance": "tests/run_losat_cache_browser_acceptance.py",
              "timingScope": "lifecycle observations, not a repeated performance gate; Generate includes render; Save includes download",
              "runs": runs}
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    return int(len(runs) != 2 or any(
        r["exitCode"] or r["blockedExternal"] or r["pageErrors"] for r in runs))


if __name__ == "__main__":
    raise SystemExit(main())
