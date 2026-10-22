"""Serve the repository for Playwright with a Pyodide snapshot for the diagram Worker.

The server is `python3 -m http.server` (GET and HEAD, from the current
directory) plus two paths: the snapshot directory that
tools/build_pyodide_snapshot.mjs wrote, at --snapshot-path, and
gbdraw/web/js/config.js with PYODIDE_SNAPSHOT_URL set to --snapshot-path, so
the Worker restores that snapshot. playwright.config.js starts it when
GBDRAW_TEST_PYODIDE_SNAPSHOT_DIR is set (docs/internal/SELECTIVE_CI.md,
"Pyodide snapshot").
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


def make_handler(snapshot_dir: str, snapshot_path: str) -> type[SimpleHTTPRequestHandler]:
    class Handler(SimpleHTTPRequestHandler):
        def log_message(self, format: str, *args: object) -> None:  # noqa: A002 - base signature
            pass

        def translate_path(self, path: str) -> str:
            clean = urlsplit(path).path
            if not clean.startswith(snapshot_path):
                return super().translate_path(path)
            relative = posixpath.normpath(unquote(clean[len(snapshot_path):]))
            if relative.startswith(("..", "/")):
                return os.path.join(snapshot_dir, "missing")
            return os.path.join(snapshot_dir, relative)

        def send_head(self) -> BinaryIO | None:
            if urlsplit(self.path).path != CONFIG_PATH:
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
    parser.add_argument("--snapshot-dir", required=True, help="Directory of snapshot.bin, state.json and fs.bin.")
    parser.add_argument("--snapshot-path", required=True, help="URL path prefix of the snapshot, ending in '/'.")
    args = parser.parse_args()
    snapshot_dir = os.path.abspath(args.snapshot_dir)
    if not args.snapshot_path.startswith("/") or not args.snapshot_path.endswith("/"):
        parser.error("--snapshot-path must start and end with '/'")
    if not os.path.isfile(os.path.join(snapshot_dir, "state.json")):
        parser.error(f"{snapshot_dir} holds no state.json; run tools/build_pyodide_snapshot.mjs first")
    handler = partial(make_handler(snapshot_dir, args.snapshot_path), directory=os.getcwd())
    ThreadingHTTPServer(("127.0.0.1", args.port), handler).serve_forever()


if __name__ == "__main__":
    main()
