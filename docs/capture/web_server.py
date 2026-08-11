"""Loopback-only static server used by documentation capture."""

from __future__ import annotations

import mimetypes
import threading
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from typing import Any

from config import (
    ISOLATION_HEADERS,
    LOCAL_BASE_URL_TEMPLATE,
    LOCAL_HOST,
    WEB_ROOT,
    validate_capture_base_url,
)


for suffix, content_type in (
    (".mjs", "text/javascript"),
    (".wasm", "application/wasm"),
    (".whl", "application/zip"),
):
    mimetypes.add_type(content_type, suffix)


class _CaptureRequestHandler(SimpleHTTPRequestHandler):
    def end_headers(self) -> None:
        for name, value in ISOLATION_HEADERS.items():
            self.send_header(name, value)
        super().end_headers()

    def log_message(self, format: str, *args: Any) -> None:
        del format, args


class CaptureWebServer:
    """Serve the packaged web app from a random loopback port."""

    def __init__(self) -> None:
        handler = partial(_CaptureRequestHandler, directory=str(WEB_ROOT))
        self._server = ThreadingHTTPServer((LOCAL_HOST, 0), handler)
        self._server.daemon_threads = True
        self._thread = threading.Thread(
            target=self._server.serve_forever,
            name="gbdraw-documentation-capture-server",
            daemon=True,
        )

    @property
    def base_url(self) -> str:
        port = int(self._server.server_address[1])
        base_url = LOCAL_BASE_URL_TEMPLATE.format(port=port)
        validate_capture_base_url(base_url)
        return base_url

    def __enter__(self) -> "CaptureWebServer":
        self._thread.start()
        return self

    def __exit__(self, exc_type: Any, exc: Any, traceback: Any) -> None:
        del exc_type, exc, traceback
        self._server.shutdown()
        self._server.server_close()
        self._thread.join(timeout=5)

