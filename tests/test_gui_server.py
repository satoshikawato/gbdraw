"""The local server of ``gbdraw gui``: one fixed address across runs, a free port
when that one is in use, and headers that make the browser revalidate every file.

One address keeps the browser's storage for it across runs: the layout choice and
the temporary input copies (PD-OI-121), whose leftovers the next run deletes.
"""
from __future__ import annotations

import socket
import threading
import urllib.request
from functools import partial
from pathlib import Path

import pytest

from gbdraw import cli


def _free_port() -> int:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as probe:
        probe.bind(("127.0.0.1", 0))
        return int(probe.getsockname()[1])


def _handler(directory: Path) -> partial[cli.IsolatedSimpleHTTPRequestHandler]:
    return partial(cli.IsolatedSimpleHTTPRequestHandler, directory=str(directory))


def test_the_gui_port_is_fixed_and_not_a_common_development_port() -> None:
    assert cli.GUI_PORT not in {3000, 5000, 8000, 8080, 8501, 8888}
    # Below the ephemeral ranges (Linux 32768-60999, Windows 49152-65535), so no outgoing connection holds it.
    assert 1024 < cli.GUI_PORT < 32768


def test_the_server_binds_the_fixed_port_when_it_is_free(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    port = _free_port()
    monkeypatch.setattr(cli, "GUI_PORT", port)
    server, fell_back = cli.bind_gui_server(_handler(tmp_path))
    try:
        assert (server.server_address, fell_back) == (("127.0.0.1", port), False)
    finally:
        server.server_close()


def test_the_server_falls_back_to_a_free_port_when_the_fixed_one_is_in_use(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as taken:
        taken.bind(("127.0.0.1", 0))
        taken.listen()
        busy = int(taken.getsockname()[1])
        monkeypatch.setattr(cli, "GUI_PORT", busy)
        server, fell_back = cli.bind_gui_server(_handler(tmp_path))
        try:
            assert fell_back is True
            assert server.server_address[1] != busy
        finally:
            server.server_close()


def test_served_files_are_revalidated_and_isolated(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    (tmp_path / "index.html").write_text("<!doctype html>", encoding="utf-8")
    monkeypatch.setattr(cli, "GUI_PORT", _free_port())
    server, _ = cli.bind_gui_server(_handler(tmp_path))
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        with urllib.request.urlopen(f"http://127.0.0.1:{server.server_address[1]}/index.html", timeout=10) as reply:
            headers = reply.headers
        assert headers["Cache-Control"] == "no-cache"
        assert headers["Cross-Origin-Opener-Policy"] == "same-origin"
        assert headers["Cross-Origin-Embedder-Policy"] == "require-corp"
    finally:
        server.shutdown()
        server.server_close()
