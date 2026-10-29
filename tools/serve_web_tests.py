"""Serve the repository for Playwright, optionally with a Pyodide snapshot for the diagram Worker.

The server is `python3 -m http.server` (GET and HEAD, from the current
directory) plus a COEP header on the LOSAT thread worker (NESTED_WORKER_PATH).
With --snapshot-dir and --snapshot-path it also serves two paths: the snapshot
directory that tools/build_pyodide_snapshot.mjs wrote, at --snapshot-path, and
gbdraw/web/js/config.js with PYODIDE_SNAPSHOT_URL set to --snapshot-path, so
the Worker restores that snapshot. playwright.config.js starts it, with the
snapshot when GBDRAW_TEST_PYODIDE_SNAPSHOT_DIR is set
(docs/internal/SELECTIVE_CI.md, "Pyodide snapshot").
"""

from __future__ import annotations

import argparse
import io
import json
import os
import posixpath
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from typing import BinaryIO
from urllib.parse import unquote, urlsplit

CONFIG_PATH = "/gbdraw/web/js/config.js"
SNAPSHOT_URL_NULL = "export const PYODIDE_SNAPSHOT_URL = null;"
# Inside a cross-origin-isolated page this worker needs a COEP header of its own.
# It is the one worker that another worker creates, and Chromium 156 fetches such
# a worker past Playwright's routes, so routeCrossOriginIsolated
# (tests/web/helpers/app-lifecycle.cjs) cannot add the header. It does nothing
# for a script that does not run as a worker, and a page that is not
# cross-origin isolated never creates this worker (threaded LOSAT needs it).
NESTED_WORKER_PATH = "/gbdraw/web/js/workers/losat-wasi-thread-worker.js"


def make_handler(snapshot_dir: str | None, snapshot_path: str | None) -> type[SimpleHTTPRequestHandler]:
    class Handler(SimpleHTTPRequestHandler):
        def log_message(self, format: str, *args: object) -> None:  # noqa: A002 - base signature
            pass

        def end_headers(self) -> None:
            if urlsplit(self.path).path == NESTED_WORKER_PATH:
                self.send_header("Cross-Origin-Embedder-Policy", "require-corp")
            super().end_headers()

        def translate_path(self, path: str) -> str:
            clean = urlsplit(path).path
            if snapshot_dir is None or snapshot_path is None or not clean.startswith(snapshot_path):
                return super().translate_path(path)
            relative = posixpath.normpath(unquote(clean[len(snapshot_path):]))
            if relative.startswith(("..", "/")):
                return os.path.join(snapshot_dir, "missing")
            return os.path.join(snapshot_dir, relative)

        def send_head(self) -> BinaryIO | None:
            if snapshot_path is None or urlsplit(self.path).path != CONFIG_PATH:
                return super().send_head()
            with open(os.path.join(self.directory, CONFIG_PATH.lstrip("/")), encoding="utf-8") as handle:
                text = handle.read()
            if SNAPSHOT_URL_NULL not in text:
                self.send_error(500, "config.js does not declare PYODIDE_SNAPSHOT_URL = null")
                return None
            body = text.replace(
                SNAPSHOT_URL_NULL, f"export const PYODIDE_SNAPSHOT_URL = {json.dumps(snapshot_path)};"
            ).encode("utf-8")
            self.send_response(200)
            self.send_header("Content-Type", "text/javascript")
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            return io.BytesIO(body)

    return Handler


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--port", type=int, required=True)
    parser.add_argument("--snapshot-dir", help="Directory of snapshot.bin, state.json and fs.bin.")
    parser.add_argument("--snapshot-path", help="URL path prefix of the snapshot, ending in '/'.")
    args = parser.parse_args()
    if (args.snapshot_dir is None) != (args.snapshot_path is None):
        parser.error("--snapshot-dir and --snapshot-path go together")
    snapshot_dir = None
    if args.snapshot_dir is not None:
        snapshot_dir = os.path.abspath(args.snapshot_dir)
        if not args.snapshot_path.startswith("/") or not args.snapshot_path.endswith("/"):
            parser.error("--snapshot-path must start and end with '/'")
        if not os.path.isfile(os.path.join(snapshot_dir, "state.json")):
            parser.error(f"{snapshot_dir} holds no state.json; run tools/build_pyodide_snapshot.mjs first")
    handler = partial(make_handler(snapshot_dir, args.snapshot_path), directory=os.getcwd())
    ThreadingHTTPServer(("127.0.0.1", args.port), handler).serve_forever()


if __name__ == "__main__":
    main()
