"""The local server of ``gbdraw gui``: one fixed address across runs, a free port
when a server answers on it, and headers that make the browser revalidate every
file against an ETag that names the file and the gbdraw version.

One address keeps the browser's storage for it across runs: the layout choice and
the temporary input copies (PD-OI-121), whose leftovers the next run deletes.
"""
from __future__ import annotations

import os
import socket
import threading
import urllib.error
import urllib.request
from email.message import Message
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


def test_the_server_falls_back_when_a_server_answers_on_the_ipv6_loopback(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # The browser opens localhost, which may reach ::1 before 127.0.0.1.
    try:
        taken = socket.socket(socket.AF_INET6, socket.SOCK_STREAM)
        taken.bind(("::1", 0))
    except OSError:
        pytest.skip("no IPv6 loopback here")
    with taken:
        taken.listen()
        busy = int(taken.getsockname()[1])
        monkeypatch.setattr(cli, "GUI_PORT", busy)
        server, fell_back = cli.bind_gui_server(_handler(tmp_path))
        try:
            assert fell_back is True
            assert server.server_address[1] != busy
        finally:
            server.server_close()


def test_the_server_falls_back_when_a_server_listens_on_every_address(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # macOS lets a bind with the reuse option share a port another program holds on 0.0.0.0.
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as taken:
        taken.bind(("0.0.0.0", 0))
        taken.listen()
        busy = int(taken.getsockname()[1])
        monkeypatch.setattr(cli, "GUI_PORT", busy)
        server, fell_back = cli.bind_gui_server(_handler(tmp_path))
        try:
            assert fell_back is True
        finally:
            server.server_close()


def test_a_restart_right_after_the_server_stopped_binds_the_fixed_port_again(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(cli, "GUI_PORT", _free_port())
    server, _ = cli.bind_gui_server(_handler(tmp_path))
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    with urllib.request.urlopen(f"http://127.0.0.1:{cli.GUI_PORT}/", timeout=10):
        pass
    server.shutdown()
    server.server_close()
    again, fell_back = cli.bind_gui_server(_handler(tmp_path))
    try:
        assert (again.server_address[1], fell_back) == (cli.GUI_PORT, False)
    finally:
        again.server_close()


class _ProbeSocket:
    """A socket whose connect never completes: each probe times out, as a refused
    loopback connection does on Windows."""

    timeouts: list[float | None] = []

    def __init__(self, *_args: object) -> None:
        pass

    def __enter__(self) -> _ProbeSocket:
        return self

    def __exit__(self, *_args: object) -> None:
        pass

    def settimeout(self, value: float | None) -> None:
        _ProbeSocket.timeouts.append(value)

    def connect_ex(self, _address: object) -> int:
        raise TimeoutError("timed out")


def test_a_probe_that_times_out_finds_no_server_and_waits_briefly(monkeypatch: pytest.MonkeyPatch) -> None:
    _ProbeSocket.timeouts = []
    monkeypatch.setattr(cli.socket, "socket", _ProbeSocket)
    assert cli._port_answers(cli.GUI_PORT) is False
    # Both loopback addresses are probed, each with a short timeout, so the start waits at most about 0.2 s.
    assert len(_ProbeSocket.timeouts) == 2
    assert all(value is not None and 0 < value <= 0.1 for value in _ProbeSocket.timeouts)


def _get(port: int, path: str, headers: dict[str, str] | None = None) -> tuple[int, Message]:
    request = urllib.request.Request(f"http://127.0.0.1:{port}{path}", headers=headers or {})
    try:
        with urllib.request.urlopen(request, timeout=10) as reply:
            reply.read()
            return reply.status, reply.headers
    except urllib.error.HTTPError as error:
        return error.code, error.headers


def test_a_cached_file_is_reused_only_when_its_etag_names_the_same_file_and_gbdraw(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    module = tmp_path / "app.js"
    module.write_text("export const a = 1;", encoding="utf-8")
    monkeypatch.setattr(cli, "GUI_PORT", _free_port())
    server, _ = cli.bind_gui_server(_handler(tmp_path))
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    port = server.server_address[1]
    try:
        status, headers = _get(port, "/app.js")
        etag = headers["ETag"]
        assert status == 200 and etag
        # A reload revalidates and downloads nothing.
        assert _get(port, "/app.js", {"If-None-Match": etag})[0] == 304
        # A date never answers 304: an older gbdraw can serve older files.
        assert _get(port, "/app.js", {"If-Modified-Since": "Fri, 01 Jan 2100 00:00:00 GMT"})[0] == 200
        # Other bytes under the same modification time, or another gbdraw version, are downloaded again.
        stat = module.stat()
        module.write_text("export const a = 22;", encoding="utf-8")
        os.utime(module, ns=(stat.st_atime_ns, stat.st_mtime_ns))
        status, changed = _get(port, "/app.js", {"If-None-Match": etag})
        assert status == 200 and changed["ETag"] != etag
        monkeypatch.setattr(cli, "__version_display__", "0.0.0-other")
        assert _get(port, "/app.js", {"If-None-Match": changed["ETag"]})[0] == 200
    finally:
        server.shutdown()
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
