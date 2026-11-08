"""The test-data manifest, its .gitignore block and the order fetch looks in."""

from __future__ import annotations

import hashlib
import json
import subprocess
import urllib.error
import urllib.request
from email.message import Message
from pathlib import Path

import pytest

from tools import test_data
from tools.test_data import InputFile

PROJECT_ROOT = Path(__file__).resolve().parents[1]
RELEASE_PREFIX = (
    f"https://github.com/{test_data.RELEASE_REPOSITORY}/releases/download/"
)


def _git(*args: str, cwd: Path, data: bytes | None = None) -> str:
    return subprocess.run(
        ["git", *args], cwd=cwd, input=data, capture_output=True, check=True
    ).stdout.decode().strip()


def test_manifest_is_in_the_form_the_tool_writes() -> None:
    files = test_data.load_manifest()
    assert test_data.MANIFEST_PATH.read_text(encoding="utf-8") == test_data.manifest_text(files)
    assert len({item.path for item in files}) == len(files)
    for item in files:
        assert item.path.startswith(("tests/test_inputs/", "examples/"))
        assert item.url.startswith(RELEASE_PREFIX)
        assert item.url.endswith(f"/{item.sha256[:12]}-{Path(item.path).name}")


def test_gitignore_holds_the_generated_block() -> None:
    text = test_data.GITIGNORE_PATH.read_text(encoding="utf-8")
    assert test_data.gitignore_block(test_data.load_manifest()) in text


def test_git_tracks_no_listed_input() -> None:
    if not (PROJECT_ROOT / ".git").exists():
        pytest.skip("not a git checkout")
    paths = [item.path for item in test_data.load_manifest()]
    assert _git("ls-files", "--", *paths, cwd=PROJECT_ROOT) == ""


def test_the_block_is_replaced_in_place_or_appended() -> None:
    old = InputFile("examples/a.gb", 1, "0" * 64, "0" * 40, RELEASE_PREFIX + "t/a")
    new = InputFile("examples/b.gb", 1, "1" * 64, "1" * 40, RELEASE_PREFIX + "t/b")
    appended = test_data.with_gitignore_block("*.pyc\n", [old])
    assert appended == "*.pyc\n\n" + test_data.gitignore_block([old])
    replaced = test_data.with_gitignore_block(appended + "/after\n", [new])
    assert replaced == "*.pyc\n\n" + test_data.gitignore_block([new]) + "/after\n"


def _input(path: str, content: bytes, git_blob: str, url: str) -> InputFile:
    return InputFile(path, len(content), hashlib.sha256(content).hexdigest(), git_blob, url)


def test_fetch_places_inputs_from_the_cache_git_history_or_the_url(tmp_path: Path) -> None:
    root = tmp_path / "checkout"
    root.mkdir()
    _git("init", "-q", cwd=root)
    cache = tmp_path / "cache"
    cache.mkdir()
    served = tmp_path / "served.gb"
    served.write_bytes(b"from the release\n")
    absent_url = (tmp_path / "absent").as_uri()

    in_history = _input(
        "tests/test_inputs/history.gb",
        b"from git history\n",
        _git("hash-object", "-w", "--stdin", cwd=root, data=b"from git history\n"),
        absent_url,
    )
    in_release = _input("examples/release.gb", served.read_bytes(), "2" * 40, served.as_uri())
    in_cache = _input("tests/test_inputs/cache.gb", b"from the cache\n", "3" * 40, absent_url)
    nowhere = _input("tests/test_inputs/nowhere.gb", b"nowhere\n", "4" * 40, absent_url)
    (cache / in_cache.sha256).write_bytes(b"from the cache\n")
    (root / "tests/test_inputs").mkdir(parents=True)
    (root / in_cache.path).write_bytes(b"stale copy\n")

    messages: list[str] = []
    missing = test_data.fetch(
        [in_history, in_release, in_cache, nowhere], root=root, cache=cache, log=messages.append
    )

    assert missing == (nowhere,)
    assert (root / in_history.path).read_bytes() == b"from git history\n"
    assert (root / in_release.path).read_bytes() == b"from the release\n"
    assert (root / in_cache.path).read_bytes() == b"from the cache\n"
    assert not (root / nowhere.path).exists()
    assert sorted(path.name for path in cache.iterdir()) == sorted(
        item.sha256 for item in (in_history, in_release, in_cache)
    )
    assert any("nowhere.gb" in message for message in messages)
    with pytest.raises(FileNotFoundError, match="tests/test_inputs/nowhere.gb"):
        test_data.fetch([nowhere], root=root, cache=cache, require=True, log=messages.append)


def test_no_copy_that_does_not_match_is_placed(tmp_path: Path) -> None:
    served = tmp_path / "served.gb"
    served.write_bytes(b"changed upstream\n")
    root = tmp_path / "not-a-repository"
    root.mkdir()
    cache = tmp_path / "cache"
    cache.mkdir()
    expected = _input("examples/x.gb", b"what the manifest says\n", "5" * 40, served.as_uri())
    (cache / expected.sha256).write_bytes(b"corrupted in place\n")
    messages: list[str] = []

    missing = test_data.fetch([expected], root=root, cache=cache, log=messages.append)

    assert missing == (expected,)
    assert not (root / expected.path).exists()
    assert (cache / expected.sha256).read_bytes() == b"corrupted in place\n"
    assert sorted(path.name for path in cache.iterdir()) == [expected.sha256]
    assert any("does not match" in message for message in messages)


TAG = "test-data-2026-11"
API_RELEASE = f"https://api.github.com/repos/{test_data.RELEASE_REPOSITORY}/releases/tags/{TAG}"


class _Response:
    def __init__(self, body: bytes) -> None:
        self._body = body

    def __enter__(self) -> "_Response":
        return self

    def __exit__(self, *exc: object) -> None:
        return None

    def read(self, size: int = -1) -> bytes:
        size = len(self._body) if size < 0 else size
        chunk, self._body = self._body[:size], self._body[size:]
        return chunk


def _serve(monkeypatch: pytest.MonkeyPatch, content: bytes, name: str, *, lookup_ok: bool = True) -> list[urllib.request.Request]:
    """Make the browser route answer 504; the API serves the tag lookup and the asset."""

    seen: list[urllib.request.Request] = []

    def urlopen(request: str | urllib.request.Request, timeout: float = 0) -> _Response:
        if isinstance(request, str):
            request = urllib.request.Request(request)
        seen.append(request)
        url = request.full_url
        if url.startswith(RELEASE_PREFIX):
            raise urllib.error.HTTPError(url, 504, "Gateway Time-out", Message(), None)
        if url == API_RELEASE:
            if not lookup_ok:
                raise urllib.error.HTTPError(url, 500, "Server Error", Message(), None)
            body = {"assets": [{"name": name, "url": "https://api.github.com/assets/7"}]}
            return _Response(json.dumps(body).encode())
        if url == "https://api.github.com/assets/7":
            assert request.get_header("Accept") == "application/octet-stream"
            return _Response(content)
        raise AssertionError(url)

    monkeypatch.setattr(urllib.request, "urlopen", urlopen)
    return seen


def _release_input(content: bytes) -> InputFile:
    name = f"{hashlib.sha256(content).hexdigest()[:12]}-big.gb"
    return _input("examples/big.gb", content, "6" * 40, f"{RELEASE_PREFIX}{TAG}/{name}")


@pytest.mark.parametrize("token", [None, "secret"])
def test_the_api_route_delivers_when_the_download_route_fails(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, token: str | None
) -> None:
    for variable in ("GITHUB_TOKEN", "GH_TOKEN"):
        monkeypatch.delenv(variable, raising=False)
    if token:
        monkeypatch.setenv("GITHUB_TOKEN", token)
    item = _release_input(b"big file\n")
    seen = _serve(monkeypatch, b"big file\n", item.url.rsplit("/", 1)[1])
    messages: list[str] = []

    missing = test_data.fetch([item], root=tmp_path / "root", cache=tmp_path / "cache", require=True, log=messages.append)

    assert missing == ()
    assert (tmp_path / "root" / item.path).read_bytes() == b"big file\n"
    asset_request = seen[-1]
    assert asset_request.full_url == "https://api.github.com/assets/7"
    assert asset_request.unredirected_hdrs.get("Authorization") == (f"Bearer {token}" if token else None)
    assert "Authorization" not in asset_request.headers
    assert any("504" in message for message in messages)


def test_a_failed_tag_lookup_leaves_the_item_missing(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    item = _release_input(b"big file\n")
    _serve(monkeypatch, b"big file\n", item.url.rsplit("/", 1)[1], lookup_ok=False)
    messages: list[str] = []

    missing = test_data.fetch([item], root=tmp_path / "root", cache=tmp_path / "cache", log=messages.append)

    assert missing == (item,)
    assert sum("release lookup" in message for message in messages) == 1
