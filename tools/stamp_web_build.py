#!/usr/bin/env python3
from __future__ import annotations

import argparse
import html
import os
import re
import subprocess
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
PYPROJECT_PATH = REPO_ROOT / "pyproject.toml"
BUILD_LABEL_MARKER = "<!-- GBDRAW_HOSTED_BUILD_LABEL -->"
COMMIT_ENV_NAMES = ("CF_PAGES_COMMIT_SHA", "GITHUB_SHA")
_PROJECT_VERSION_RE = re.compile(r'^version\s*=\s*"([^"]+)"\s*$', re.MULTILINE)


def read_project_version() -> str:
    pyproject_text = PYPROJECT_PATH.read_text(encoding="utf-8")
    match = _PROJECT_VERSION_RE.search(pyproject_text)
    if match is None:
        raise RuntimeError(f"Could not determine project version from {PYPROJECT_PATH}")
    return match.group(1)


def resolve_commit_sha() -> str:
    for env_name in COMMIT_ENV_NAMES:
        value = os.environ.get(env_name, "").strip()
        if value:
            return value
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"],
            cwd=REPO_ROOT,
            stderr=subprocess.DEVNULL,
            text=True,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def render_build_label(*, project_version: str | None = None, commit_sha: str | None = None) -> str:
    version = project_version or read_project_version()
    full_commit = (commit_sha or resolve_commit_sha()).strip() or "unknown"
    short_commit = "unknown" if full_commit == "unknown" else full_commit[:7]
    label = f"v{version}+{short_commit}"
    return (
        '<span class="hosted-build-label" '
        f'title="Commit {html.escape(full_commit, quote=True)}">'
        f"Version: {html.escape(label)}</span>"
    )


def stamp_web_build(
    web_root: Path,
    *,
    project_version: str | None = None,
    commit_sha: str | None = None,
) -> Path:
    index_path = web_root / "index.html"
    source = index_path.read_text(encoding="utf-8")
    if BUILD_LABEL_MARKER not in source:
        if 'class="hosted-build-label"' in source:
            return index_path
        raise RuntimeError(f"{BUILD_LABEL_MARKER} not found in {index_path}")
    updated = source.replace(
        BUILD_LABEL_MARKER,
        render_build_label(project_version=project_version, commit_sha=commit_sha),
        1,
    )
    index_path.write_text(updated, encoding="utf-8")
    return index_path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Stamp hosted gbdraw web builds with version and commit.")
    parser.add_argument(
        "web_root",
        nargs="?",
        type=Path,
        default=Path("gbdraw/web"),
        help="Directory containing the index.html to stamp.",
    )
    args = parser.parse_args(argv)
    stamped_path = stamp_web_build(args.web_root)
    print(f"Stamped hosted build label in {stamped_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
