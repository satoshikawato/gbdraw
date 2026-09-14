"""Release admission uses real Git refs; archive verification never publishes."""

from __future__ import annotations

import io
import json
from pathlib import Path
from runpy import run_path
import subprocess
import tarfile
import zipfile

import pytest

from tests.test_web_packaging import _load_verify_module


ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def release_source(tmp_path, monkeypatch):
    def git(*args):
        return subprocess.check_output(["git", *args], cwd=tmp_path, text=True).strip()

    git("init", "-b", "main")
    git("config", "user.name", "Release test")
    git("config", "user.email", "release-test@example.invalid")
    (tmp_path / "pyproject.toml").write_text('[project]\nversion = "0.14.0rc1"\n')
    git("add", "pyproject.toml")
    git("commit", "-m", "Release fixture")
    sha = git("rev-parse", "HEAD")
    git("update-ref", "refs/remotes/origin/main", sha)
    git("tag", "v0.14.0rc1")
    check = run_path(str(ROOT / "tools/check_release_source.py"))["check_release_source"]
    # Reuse the real version reader, directed at the isolated fixture's TOML.
    read_version = run_path(str(ROOT / "gbdraw/_build_support.py"))["read_project_version"]
    read_version.__globals__["PYPROJECT_PATH"] = tmp_path / "pyproject.toml"
    monkeypatch.setitem(check.__globals__, "REPO_ROOT", tmp_path)
    monkeypatch.setitem(check.__globals__, "run_path", lambda _: {"read_project_version": read_version})
    event_path = tmp_path / "release-event.json"
    event_path.write_text(json.dumps({
        "action": "published",
        "release": {"tag_name": "v0.14.0rc1", "draft": False, "prerelease": True},
    }))
    for name, value in {
        "GITHUB_EVENT_NAME": "release", "GITHUB_REF_TYPE": "tag",
        "GITHUB_REF": "refs/tags/v0.14.0rc1", "GITHUB_SHA": sha,
        "GITHUB_EVENT_PATH": str(event_path),
    }.items():
        monkeypatch.setenv(name, value)
    return check, git, tmp_path


@pytest.mark.parametrize("annotated", [False, True])
def test_release_admits_exact_main_tag(release_source, monkeypatch, annotated):
    check, git, _ = release_source
    if annotated:
        git("tag", "-f", "-a", "v0.14.0rc1", "-m", "Release candidate")
        monkeypatch.setenv("GITHUB_SHA", git("rev-parse", "refs/tags/v0.14.0rc1"))
    assert check() == "0.14.0rc1"


def test_release_admits_final_and_historical_main_source(release_source, monkeypatch):
    check, git, root = release_source
    (root / "pyproject.toml").write_text('[project]\nversion = "0.14.0"\n')
    git("commit", "-am", "Final version")
    final = git("rev-parse", "HEAD")
    git("tag", "v0.14.0")
    git("commit", "--allow-empty", "-m", "Later main progress")
    git("update-ref", "refs/remotes/origin/main", git("rev-parse", "HEAD"))
    git("checkout", "--detach", final)
    monkeypatch.setenv("GITHUB_SHA", final)
    monkeypatch.setenv("GITHUB_REF", "refs/tags/v0.14.0")
    (root / "release-event.json").write_text(json.dumps({
        "action": "published",
        "release": {"tag_name": "v0.14.0", "draft": False, "prerelease": False},
    }))
    assert check() == "0.14.0"


@pytest.mark.parametrize("ref", ["refs/remotes/origin/main", "refs/tags/v0.14.0rc1"])
def test_release_rejects_missing_source_refs(release_source, ref):
    check, git, _ = release_source
    git("update-ref", "-d", ref)
    with pytest.raises(subprocess.CalledProcessError):
        check()


@pytest.mark.parametrize("key,value", [
    ("GITHUB_REF", "refs/tags/v0.14.0"),
    ("GITHUB_REF", "refs/heads/main"),
    ("GITHUB_REF", "v0.14.0rc1"),
    ("GITHUB_REF", "refs/tags/v0.14.0b0"),
    ("GITHUB_REF", "refs/tags/v0.14.0rc1/extra"),
    ("GITHUB_REF_TYPE", "branch"),
    ("GITHUB_EVENT_NAME", "workflow_dispatch"),
    ("GITHUB_EVENT_NAME", "push"),
    ("GITHUB_SHA", "HEAD"),
    ("GITHUB_SHA", "0" * 40),
])
def test_release_rejects_wrong_event_identity(release_source, monkeypatch, key, value):
    check, _, _ = release_source
    monkeypatch.setenv(key, value)
    with pytest.raises(ValueError):
        check()


@pytest.mark.parametrize("action", ["created", "edited", "deleted", "unpublished", "released", "prereleased"])
def test_release_rejects_other_release_actions(release_source, action):
    check, _, root = release_source
    path = root / "release-event.json"
    event = json.loads(path.read_text())
    event["action"] = action
    path.write_text(json.dumps(event))
    with pytest.raises(ValueError, match="Published GitHub Release and tag ref must agree"):
        check()


@pytest.mark.parametrize("release", [
    None,
    {},
    {"tag_name": "v0.14.0rc1", "draft": True},
    {"tag_name": "v0.14.0rc1"},
    {"tag_name": "v0.13.0", "draft": False},
])
def test_release_rejects_unpublished_or_mismatched_release(release_source, release):
    check, _, root = release_source
    (root / "release-event.json").write_text(json.dumps({"action": "published", "release": release}))
    with pytest.raises(ValueError, match="Published GitHub Release and tag ref must agree"):
        check()


@pytest.mark.parametrize("state", ["non-main", "side-parent", "wrong-checkout", "moved-tag", "dirty"])
def test_release_rejects_wrong_source(release_source, monkeypatch, state):
    check, git, root = release_source
    if state == "dirty":
        (root / "pyproject.toml").write_text('[project]\nversion = "0.14.0rc1"\n# dirty\n')
    else:
        initial = git("rev-parse", "HEAD")
        git("switch", "-c", "candidate")
        git("commit", "--allow-empty", "-m", "Unpromoted candidate")
        candidate = git("rev-parse", "HEAD")
        if state != "wrong-checkout":
            git("tag", "-f", "v0.14.0rc1")
        if state in {"non-main", "side-parent"}:
            monkeypatch.setenv("GITHUB_SHA", candidate)
        if state == "side-parent":
            git("switch", "main")
            git("commit", "--allow-empty", "-m", "Main progress")
            git("merge", "--no-ff", "candidate", "-m", "Integrate candidate")
            git("update-ref", "refs/remotes/origin/main", git("rev-parse", "HEAD"))
            git("checkout", "--detach", candidate)
        if state == "moved-tag":
            git("checkout", "--detach", initial)
    with pytest.raises(ValueError):
        check()


@pytest.mark.parametrize("fault", ["missing-wheel", "missing-sdist", "extra", "wheel-version", "sdist-version", "missing-content", "missing-wheel-content", None])
def test_release_distribution_pair_reuses_package_inspectors(tmp_path, monkeypatch, fault):
    verify = _load_verify_module()
    version = verify.BUILD_SUPPORT.read_project_version()
    metadata = f"Name: gbdraw\nVersion: {version}\n".encode()
    wrong = b"Name: gbdraw\nVersion: 0.0.0\n"
    wheel = tmp_path / verify.BUILD_SUPPORT.expected_browser_wheel_name()
    sdist = tmp_path / f"gbdraw-{version}.tar.gz"
    if fault != "missing-wheel":
        with zipfile.ZipFile(wheel, "w") as archive:
            archive.writestr(f"gbdraw-{version}.dist-info/METADATA", wrong if fault == "wheel-version" else metadata)
    if fault != "missing-sdist":
        with tarfile.open(sdist, "w:gz") as archive:
            data = wrong if fault == "sdist-version" else metadata
            member = tarfile.TarInfo(f"gbdraw-{version}/PKG-INFO")
            member.size = len(data)
            archive.addfile(member, io.BytesIO(data))
    if fault == "extra":
        (tmp_path / "unexpected.txt").write_text("not a distribution")
    calls = []
    if fault != "missing-content":
        monkeypatch.setattr(verify, "inspect_sdist", lambda path: calls.append(path))
    if fault not in {"missing-content", "missing-wheel-content"}:
        monkeypatch.setattr(verify, "inspect_wheel", lambda path: calls.append(path))
    if fault:
        with pytest.raises((RuntimeError, FileNotFoundError)):
            verify.inspect_distributions(tmp_path)
    else:
        verify.inspect_distributions(tmp_path)
        assert calls == [sdist, wheel]
