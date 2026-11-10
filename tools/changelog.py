#!/usr/bin/env python3
"""Check changelog fragments, or fold them into the [Unreleased] list of CHANGELOG.md.

Each pull request writes its entries to its own file under changelog.d/, so two
pull requests never edit the same lines of CHANGELOG.md. The release session
runs ``fold`` once, before the promotion: it puts every fragment's entries at
the top of the [Unreleased] list, newest merge first (the order in which they
reached the current branch), and deletes the fragments.

    python tools/changelog.py check
    python tools/changelog.py fold
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
FRAGMENT_DIRECTORY = "changelog.d"
UNRELEASED_HEADING = "## [Unreleased]"
CATEGORIES = ("Added", "Changed", "Deprecated", "Removed", "Fixed", "Security", "Performance", "Known issue")
FRAGMENT_NAME = re.compile(r"[a-z0-9][a-z0-9._-]*\.md")
ENTRY_START = re.compile(rf"- (?:{'|'.join(CATEGORIES)}): \S")


def fragment_paths(root: Path) -> list[Path]:
    """The fragment files of a checkout; README.md describes them and is not one."""
    directory = root / FRAGMENT_DIRECTORY
    return sorted(path for path in directory.iterdir() if path.name != "README.md")


def fragment_errors(path: Path) -> list[str]:
    """Why a fragment cannot be folded as written, or nothing."""
    name = f"{FRAGMENT_DIRECTORY}/{path.name}"
    if not path.is_file() or not FRAGMENT_NAME.fullmatch(path.name):
        return [f"{name}: a fragment is a file named like fix-legend-order.md (lower case, digits, '.', '_', '-')"]
    text = path.read_text(encoding="utf-8")
    if not text.endswith("\n") or text.endswith("\n\n"):
        return [f"{name}: end the file with exactly one newline"]
    errors = []
    lines = text[:-1].split("\n")
    if not ENTRY_START.match(lines[0]):
        errors.append(f"{name}:1: start with '- <Category>: ' where Category is one of {', '.join(CATEGORIES)}")
    for number, line in enumerate(lines, start=1):
        if line != line.rstrip():
            errors.append(f"{name}:{number}: trailing whitespace")
        elif line.startswith("- "):
            if not ENTRY_START.match(line):
                errors.append(f"{name}:{number}: an entry starts with '- <Category>: ' ({', '.join(CATEGORIES)})")
        elif not line.startswith("  ") or not line.strip():
            errors.append(f"{name}:{number}: continue an entry on lines indented by two spaces, with no blank lines")
    return errors


def check(root: Path) -> list[str]:
    return [error for path in fragment_paths(root) for error in fragment_errors(path)]


def _merged_at(root: Path, path: Path) -> int:
    """When the fragment reached the checked-out branch: the commit time of the first-parent commit
    that added it. A fragment that is not committed yet counts as the newest."""
    result = subprocess.run(
        ["git", "log", "--first-parent", "--diff-filter=A", "--format=%ct", "-1", "--",
         path.relative_to(root).as_posix()],
        cwd=root, capture_output=True, text=True, check=True,
    )
    stamp = result.stdout.strip()
    return int(stamp) if stamp else sys.maxsize


def fold(root: Path) -> int:
    """Insert every fragment's entries at the top of the [Unreleased] list, newest first, and
    delete the fragments. Returns the number of fragments folded."""
    errors = check(root)
    if errors:
        raise ValueError("\n".join(errors))
    paths = sorted(fragment_paths(root), key=lambda path: (-_merged_at(root, path), path.name))
    if not paths:
        return 0
    changelog = root / "CHANGELOG.md"
    lines = changelog.read_text(encoding="utf-8").split("\n")
    if lines.count(UNRELEASED_HEADING) != 1:
        raise ValueError(f"CHANGELOG.md must have exactly one '{UNRELEASED_HEADING}' heading")
    position = lines.index(UNRELEASED_HEADING) + 1
    while position < len(lines) and not lines[position]:
        position += 1
    entries = [line for path in paths for line in path.read_text(encoding="utf-8")[:-1].split("\n")]
    following = lines[position:]
    if following and following[0].startswith("## "):
        entries.append("")
    lines[position:] = [*entries, *following]
    changelog.write_text("\n".join(lines), encoding="utf-8")
    for path in paths:
        path.unlink()
    return len(paths)


def main(argv: list[str]) -> int:
    if argv not in (["check"], ["fold"]):
        print(__doc__.strip(), file=sys.stderr)
        return 2
    try:
        if argv == ["check"]:
            errors = check(REPO_ROOT)
            if errors:
                print("\n".join(errors), file=sys.stderr)
                return 1
            print(f"{len(fragment_paths(REPO_ROOT))} changelog fragment(s) can be folded.")
            return 0
        count = fold(REPO_ROOT)
    except ValueError as error:
        print(error, file=sys.stderr)
        return 1
    print(f"Folded {count} changelog fragment(s) into the [Unreleased] list of CHANGELOG.md.")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
