from __future__ import annotations

from contextlib import ExitStack
import hashlib
import io
import json
import stat
import tarfile
from unittest.mock import Mock
import zipfile

import pytest

from gbdraw import losat_setup as setup
from gbdraw.exceptions import ValidationError


@pytest.mark.parametrize("system,machine,target", [
    ("Linux", "x86_64", "x86_64-unknown-linux-gnu"),
    ("Windows", "AMD64", "x86_64-pc-windows-msvc"),
    ("Windows", "x64", "x86_64-pc-windows-msvc"),
    ("Darwin", "arm64", "aarch64-apple-darwin"),
    ("Darwin", "aarch64", "aarch64-apple-darwin"),
    ("Darwin", "x86_64", "x86_64-apple-darwin"),  # Intel / Rosetta Python
])
def test_target_uses_running_python(monkeypatch, system, machine, target):
    monkeypatch.setattr(setup.platform, "system", lambda: system)
    monkeypatch.setattr(setup.platform, "machine", lambda: machine)
    monkeypatch.setattr(setup.platform, "libc_ver", lambda: ("glibc", "2.34"))
    assert setup.runtime_target() == target


@pytest.mark.parametrize("libc,version", [("musl", "1.2.5"), ("glibc", "2.31"), ("", ""), ("glibc", "unknown")])
def test_unsupported_linux_runtime(monkeypatch, libc, version):
    monkeypatch.setattr(setup.platform, "system", lambda: "Linux")
    monkeypatch.setattr(setup.platform, "machine", lambda: "x86_64")
    monkeypatch.setattr(setup.platform, "libc_ver", lambda: (libc, version))
    with pytest.raises(ValidationError, match="glibc 2.34"):
        setup.runtime_target()


@pytest.mark.parametrize("system,machine,bits", [("Linux", "arm64", 8), ("Windows", "arm64", 8), ("FreeBSD", "amd64", 8), ("Linux", "x86_64", 4)])
def test_unsupported_target_is_explicit(monkeypatch, system, machine, bits):
    monkeypatch.setattr(setup.platform, "system", lambda: system)
    monkeypatch.setattr(setup.platform, "machine", lambda: machine)
    monkeypatch.setattr(setup.struct, "calcsize", lambda _: bits)
    with pytest.raises(ValidationError, match="Unsupported LOSAT platform"):
        setup.runtime_target()


@pytest.fixture
def distribution(tmp_path, monkeypatch):
    binary = b"#!/bin/sh\nprintf 'losat 0.1.0\\n'\n"
    target = "x86_64-unknown-linux-gnu"
    entry = {"filename": f"LOSAT-0.1.0-{target}.tar.gz", "binary": "LOSAT", "binary_sha256": hashlib.sha256(binary).hexdigest(), "binary_size": len(binary)}
    lock = {"schema_version": 1, "version": "0.1.0", "candidate_sha": "a" * 40, "artifacts": {target: entry}}
    monkeypatch.setattr(setup, "read_release_lock", lambda: lock)
    monkeypatch.setattr(setup, "runtime_target", lambda: target)
    monkeypatch.setattr(setup, "cache_root", lambda: tmp_path / "cache")
    monkeypatch.setattr(setup.subprocess, "run", lambda *a, **k: Mock(returncode=0, stdout="losat 0.1.0\n", stderr=""))

    def serve(*, archive_format="tar", extra=None, link=False, corrupt=False):
        suffix = ".zip" if archive_format == "zip" else ".tar.gz"
        entry["filename"] = f"LOSAT-0.1.0-{target}{suffix}"
        root = entry["filename"].removesuffix(suffix)
        files = {
            f"{root}/LOSAT": binary,
            f"{root}/LICENSE": b"license",
            f"{root}/README.md": b"readme",
            f"{root}/RELEASE-METADATA.json": json.dumps({"candidate_sha": lock["candidate_sha"], "release": "v0.1.0", "artifact": {"target": target}, "binary": {"sha256": entry["binary_sha256"]}}).encode(),
        }
        if extra:
            files[extra] = b"unexpected"
        buffer = io.BytesIO()
        if archive_format == "zip":
            with zipfile.ZipFile(buffer, "w") as archive:
                for name, data in files.items():
                    info = zipfile.ZipInfo(name)
                    info.external_attr = ((stat.S_IFLNK if link else stat.S_IFREG) | 0o755) << 16
                    archive.writestr(info, data)
        else:
            with tarfile.open(fileobj=buffer, mode="w:gz") as archive:
                for name, data in files.items():
                    info = tarfile.TarInfo(name)
                    info.size = len(data)
                    if link:
                        info.type = tarfile.SYMTYPE
                        info.linkname = "/outside"
                    archive.addfile(info, io.BytesIO(data))
        payload = b"not an archive" if corrupt else buffer.getvalue()
        entry.update(size=len(payload), sha256=hashlib.sha256(payload).hexdigest())
        def download(url, timeout):
            assert url == f"{setup.RELEASE_URL}/v0.1.0/{entry['filename']}"
            response = io.BytesIO(payload)
            response.geturl = lambda: url
            return response
        monkeypatch.setattr(setup, "urlopen", download)
        return lock
    return serve


@pytest.mark.parametrize("archive_format", ["tar", "zip"])
def test_install_atomic_and_offline_reuse(distribution, monkeypatch, archive_format):
    distribution(archive_format=archive_format)
    binary = setup.setup_losat()
    assert binary.name == "LOSAT" and binary.is_absolute()
    assert binary.stat().st_mode & stat.S_IXUSR
    monkeypatch.setattr(setup, "urlopen", Mock(side_effect=AssertionError("network forbidden")))
    assert setup.setup_losat() == binary
    assert setup.managed_losat() == binary


@pytest.mark.parametrize("options", [
    {"extra": "../outside"}, {"extra": "/absolute"}, {"extra": "unexpected"},
    {"link": True}, {"archive_format": "zip", "link": True},
    {"archive_format": "zip", "extra": "C:\\outside"}, {"corrupt": True},
])
def test_bad_archive_never_creates_cache(distribution, options):
    distribution(**options)
    with pytest.raises(ValidationError, match="setup failed"):
        setup.setup_losat()
    assert setup.managed_losat() is None


def test_archive_hash_mismatch_before_extraction(distribution):
    lock = distribution()
    next(iter(lock["artifacts"].values()))["sha256"] = "0" * 64
    with pytest.raises(ValidationError, match="SHA-256 mismatch"):
        setup.setup_losat()
    assert setup.managed_losat() is None


def test_cache_corruption_does_not_fall_back(distribution, monkeypatch):
    distribution()
    binary = setup.setup_losat()
    binary.write_bytes(b"corrupt")
    monkeypatch.setattr(setup, "urlopen", Mock(side_effect=AssertionError("network forbidden")))
    with pytest.raises(ValidationError, match="binary checksum mismatch"):
        setup.setup_losat()
    from gbdraw.analysis import protein_colinearity as protein
    monkeypatch.setattr(protein, "_conda_losatp_runtime", lambda: None)
    with ExitStack() as stack, pytest.raises(ValidationError, match="binary checksum mismatch"):
        protein._resolve_protein_blastp_runtime("losat", None, stack)


def test_interrupted_download_can_be_retried(distribution, monkeypatch):
    distribution()
    monkeypatch.setattr(setup, "urlopen", Mock(side_effect=KeyboardInterrupt))
    with pytest.raises(KeyboardInterrupt):
        setup.setup_losat()
    assert setup.managed_losat() is None
    distribution()
    assert setup.setup_losat().is_file()


def test_concurrent_installer_fails_without_touching_cache(distribution):
    lock = distribution()
    parent = setup.cache_root() / lock["version"]
    parent.mkdir(parents=True)
    target = setup.runtime_target()
    with setup._installation_lock(parent / f"{target}.lock"):
        with pytest.raises(ValidationError, match="Another LOSAT setup"):
            setup.setup_losat()
        assert setup.managed_losat() is None
    assert setup.setup_losat().is_file()


def test_version_failure_never_installs(distribution, monkeypatch):
    distribution()
    monkeypatch.setattr(setup.subprocess, "run", lambda *a, **k: Mock(returncode=0, stdout="losat 9.9.9", stderr=""))
    with pytest.raises(ValidationError, match="version check failed"):
        setup.setup_losat()
    assert setup.managed_losat() is None


def test_managed_resolver_precedes_bundled_and_explicit_precedes_managed(distribution, monkeypatch):
    from gbdraw.analysis import protein_colinearity as protein
    distribution()
    binary = setup.setup_losat()
    monkeypatch.setattr(protein, "_conda_losatp_runtime", lambda: None)
    monkeypatch.setattr(protein, "_bundled_losatp_resource", Mock(side_effect=AssertionError("bundled discovery")))
    with ExitStack() as stack:
        assert protein._resolve_protein_blastp_runtime("losat", None, stack).executable == str(binary)
        assert protein._resolve_protein_blastp_runtime("custom", None, stack).executable == "custom"
        assert protein._resolve_protein_blastp_runtime("losat", "blastp", stack).executable == "blastp"


def test_unpublished_lock_never_downloads(tmp_path, monkeypatch):
    (tmp_path / "losat-release.json").write_text(json.dumps({
        "schema_version": 1, "version": "0.1.0", "candidate_sha": None, "artifacts": {},
    }))
    monkeypatch.setattr(setup.resources, "files", lambda _: tmp_path)
    monkeypatch.setattr(setup, "urlopen", Mock(side_effect=AssertionError("network forbidden")))
    assert setup.managed_losat() is None
    with pytest.raises(ValidationError, match="no verified public release lock"):
        setup.setup_losat()


def test_malformed_packaged_lock_fails(tmp_path, monkeypatch):
    (tmp_path / "losat-release.json").write_text('{"schema_version":1,"version":"0.1.0","candidate_sha":null,"artifacts":[]}')
    monkeypatch.setattr(setup.resources, "files", lambda _: tmp_path)
    with pytest.raises(ValidationError, match="Invalid packaged LOSAT release lock"):
        setup.read_release_lock()
