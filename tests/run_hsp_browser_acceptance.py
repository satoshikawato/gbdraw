"""Run the existing real-helper/replay acceptance offline at two viewports.

No analysis or render implementation is replaced here. The reused acceptance
also injects cancellation/render failure and prohibits fresh LOSAT search while
reusing its verified saved raw fixture.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from unittest.mock import patch
from urllib.parse import urlparse

from playwright.sync_api import Browser

from tests.run_losat_cache_browser_acceptance import _run_python_adapter


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
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

        with patch.object(Browser, "new_context", new_context):
            code = _run_python_adapter()
        run["paths"] = sorted(run["paths"])
        run["exitCode"] = code
        runs.append(run)
        if code or run["blockedExternal"] or run["pageErrors"]:
            break
    report = {"artifact": "source SPA with generated local browser wheel",
              "network": "fresh contexts; external requests blocked before navigation",
              "acceptance": "tests/run_losat_cache_browser_acceptance.py",
              "runs": runs}
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    return int(len(runs) != 2 or any(
        r["exitCode"] or r["blockedExternal"] or r["pageErrors"] for r in runs))


if __name__ == "__main__":
    raise SystemExit(main())
