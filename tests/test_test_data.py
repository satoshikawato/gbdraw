"""The test-data manifest, its .gitignore block and the order fetch looks in."""

from __future__ import annotations

import hashlib
import subprocess
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
