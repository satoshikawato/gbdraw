#!/usr/bin/env python3
"""Admit only a published release's exact version-tagged main source."""

from __future__ import annotations

import json
import os
from pathlib import Path
import re
from runpy import run_path
import subprocess


REPO_ROOT = Path(__file__).resolve().parents[1]


def check_release_source() -> str:
    ref = os.environ["GITHUB_REF"]
    if (
        os.environ["GITHUB_EVENT_NAME"] != "release"
        or os.environ["GITHUB_REF_TYPE"] != "tag"
        or not re.fullmatch(r"refs/tags/v\d+\.\d+\.\d+(?:rc\d+)?", ref)
    ):
        raise ValueError("Release source requires a published final or RC GitHub Release")
    event = json.loads(Path(os.environ["GITHUB_EVENT_PATH"]).read_text())
    release = event.get("release") if isinstance(event, dict) else None
    if (
        not isinstance(release, dict)
        or event.get("action") != "published"
        or release.get("draft") is not False
        or ref != f"refs/tags/{release.get('tag_name')}"
    ):
        raise ValueError("Published GitHub Release and tag ref must agree")
    version = run_path(str(REPO_ROOT / "gbdraw/_build_support.py"))["read_project_version"]()
    if ref != f"refs/tags/v{version}":
        raise ValueError("Release tag does not exactly match the project version")
    event_sha = os.environ["GITHUB_SHA"]
    if not re.fullmatch(r"[0-9a-f]{40}", event_sha):
        raise ValueError("Release event SHA must be an exact commit or tag object ID")

    def git(*args: str) -> str:
        return subprocess.check_output(["git", *args], cwd=REPO_ROOT, text=True).strip()

    target = git("rev-parse", "--verify", f"{ref}^{{commit}}")
    tag_object = git("rev-parse", "--verify", ref)
    if event_sha not in {tag_object, target} or git("rev-parse", "HEAD") != target:
        raise ValueError("Checkout, event SHA and release tag target must agree")
    # A release source is a main integration commit, not a side-branch commit
    # merely reachable through a merge. Historical main release commits work too.
    if target not in git("rev-list", "--first-parent", "refs/remotes/origin/main").splitlines():
        raise ValueError("Release tag target is not on main's first-parent history")
    if git("status", "--porcelain", "--untracked-files=no"):
        raise ValueError("Release checkout has modified tracked source")
    print(f"Release source admitted: {ref} -> {target}; version={version}")
    return version


if __name__ == "__main__":
    check_release_source()
