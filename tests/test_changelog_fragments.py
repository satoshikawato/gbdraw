from __future__ import annotations

import os
import subprocess
from pathlib import Path
from runpy import run_path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
CHANGELOG_TOOL = run_path(str(REPO_ROOT / "tools/changelog.py"))
check = CHANGELOG_TOOL["check"]
fold = CHANGELOG_TOOL["fold"]

CHANGELOG = """# Changelog

Intro.

## [Unreleased]

- Fixed: an entry already in the list.

## [0.14.5](./docs/RELEASE_NOTES_0.14.5.md) - 2026-11-03

- Added: a released entry.
"""


def _write(root: Path, name: str, text: str) -> None:
    path = root / "changelog.d" / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _git(root: Path, *args: str, when: str | None = None) -> None:
    environment = {**os.environ, "GIT_AUTHOR_NAME": "t", "GIT_AUTHOR_EMAIL": "t@example.invalid",
                   "GIT_COMMITTER_NAME": "t", "GIT_COMMITTER_EMAIL": "t@example.invalid"}
    if when:
        environment.update(GIT_AUTHOR_DATE=when, GIT_COMMITTER_DATE=when)
    subprocess.run(["git", *args], cwd=root, env=environment, check=True, capture_output=True)


def test_fragments_in_the_repository_can_be_folded() -> None:
    assert check(REPO_ROOT) == []


@pytest.mark.parametrize("category", CHANGELOG_TOOL["CATEGORIES"])
def test_every_category_reads_as_an_entry(tmp_path: Path, category: str) -> None:
    _write(tmp_path, "fix.md", f"- {category}: x.\n  More of it.\n")
    assert check(tmp_path) == []


@pytest.mark.parametrize(
    ("name", "text", "message"),
    [
        ("Fix_Legend.md", "- Fixed: x.\n", "a fragment is a file named like"),
        ("notes.txt", "- Fixed: x.\n", "a fragment is a file named like"),
        ("fix.md", "- Fixed: x.", "end the file with exactly one newline"),
        ("fix.md", "- Fixed: x.\n\n", "end the file with exactly one newline"),
        ("fix.md", "Fixed: x.\n", "start with '- <Category>: '"),
        ("fix.md", "- Bugfix: x.\n", "start with '- <Category>: '"),
        ("fix.md", "- Fixed: x.\n- Note: y.\n", "fix.md:2: an entry starts with"),
        ("fix.md", "- Fixed: x.\n\n- Added: y.\n", "fix.md:2: continue an entry"),
        ("fix.md", "- Fixed: x\ny.\n", "fix.md:2: continue an entry"),
        ("fix.md", "- Fixed: x. \n", "fix.md:1: trailing whitespace"),
    ],
)
def test_a_fragment_that_does_not_read_as_entries_is_reported(
    tmp_path: Path, name: str, text: str, message: str
) -> None:
    _write(tmp_path, "README.md", "# Changelog fragments\n\nAnything.\n")
    _write(tmp_path, name, text)
    errors = check(tmp_path)
    assert any(message in error for error in errors), errors


def test_fold_puts_the_newest_merge_first_and_deletes_the_fragments(tmp_path: Path) -> None:
    (tmp_path / "CHANGELOG.md").write_text(CHANGELOG, encoding="utf-8")
    _write(tmp_path, "README.md", "# Changelog fragments\n")
    _git(tmp_path, "init", "--quiet")
    _git(tmp_path, "add", ".")
    _git(tmp_path, "commit", "--quiet", "-m", "base", when="2026-11-01T00:00:00Z")
    _write(tmp_path, "fix-older.md", "- Fixed: the older merge.\n  Second line.\n")
    _git(tmp_path, "add", ".")
    _git(tmp_path, "commit", "--quiet", "-m", "older", when="2026-11-02T00:00:00Z")
    _write(tmp_path, "add-newer.md", "- Added: the newer merge.\n- Changed: a second entry of it.\n")
    _git(tmp_path, "add", ".")
    _git(tmp_path, "commit", "--quiet", "-m", "newer", when="2026-11-03T00:00:00Z")
    _write(tmp_path, "local.md", "- Removed: not committed yet.\n")

    assert fold(tmp_path) == 3

    assert (tmp_path / "CHANGELOG.md").read_text(encoding="utf-8") == CHANGELOG.replace(
        "## [Unreleased]\n\n",
        "## [Unreleased]\n\n"
        "- Removed: not committed yet.\n"
        "- Added: the newer merge.\n"
        "- Changed: a second entry of it.\n"
        "- Fixed: the older merge.\n"
        "  Second line.\n",
    )
    assert sorted(path.name for path in (tmp_path / "changelog.d").iterdir()) == ["README.md"]
    assert fold(tmp_path) == 0


def test_fold_into_an_empty_unreleased_list_keeps_the_next_heading_apart(tmp_path: Path) -> None:
    changelog = CHANGELOG.replace("- Fixed: an entry already in the list.\n\n", "")
    (tmp_path / "CHANGELOG.md").write_text(changelog, encoding="utf-8")
    _git(tmp_path, "init", "--quiet")
    _git(tmp_path, "add", ".")
    _git(tmp_path, "commit", "--quiet", "-m", "base")
    _write(tmp_path, "fix.md", "- Fixed: the only entry.\n")

    assert fold(tmp_path) == 1

    assert (tmp_path / "CHANGELOG.md").read_text(encoding="utf-8") == changelog.replace(
        "## [Unreleased]\n\n", "## [Unreleased]\n\n- Fixed: the only entry.\n\n"
    )


def test_fold_refuses_fragments_that_do_not_read_as_entries(tmp_path: Path) -> None:
    (tmp_path / "CHANGELOG.md").write_text(CHANGELOG, encoding="utf-8")
    _git(tmp_path, "init", "--quiet")
    _write(tmp_path, "fix.md", "Fixed: x.\n")

    with pytest.raises(ValueError, match="start with '- <Category>: '"):
        fold(tmp_path)
    assert (tmp_path / "CHANGELOG.md").read_text(encoding="utf-8") == CHANGELOG
    assert (tmp_path / "changelog.d" / "fix.md").exists()
