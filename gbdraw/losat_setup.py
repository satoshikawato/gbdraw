"""Explicit installation and offline resolution of a pinned LOSAT release."""

from __future__ import annotations

import argparse
from contextlib import contextmanager
import hashlib
from importlib import resources
import json
import os
from pathlib import Path
import platform
import re
import shutil
import stat
import struct
import subprocess
import tarfile
import tempfile
from urllib.request import urlopen
import zipfile

from .exceptions import ValidationError


RELEASE_URL = "https://github.com/satoshikawato/LOSAT/releases/download"
TARGETS = {
    ("Linux", "x86_64"): "x86_64-unknown-linux-gnu",
    ("Windows", "x86_64"): "x86_64-pc-windows-msvc",
    ("Darwin", "aarch64"): "aarch64-apple-darwin",
    ("Darwin", "x86_64"): "x86_64-apple-darwin",
}


def runtime_target() -> str:
    # platform.machine() follows the running interpreter under Rosetta.
    machine = platform.machine().lower()
    machine = {"amd64": "x86_64", "x64": "x86_64", "arm64": "aarch64"}.get(machine, machine)
    target = TARGETS.get((platform.system(), machine))
    if target is None or struct.calcsize("P") != 8:
        raise ValidationError(
            f"Unsupported LOSAT platform: {platform.system()} / {platform.machine()}. "
            "Supported: Linux x64 (glibc), Windows x64, macOS arm64 and x64 (64-bit Python)."
        )
    if platform.system() == "Linux":
        libc, version = platform.libc_ver()
        if libc != "glibc" or not re.fullmatch(r"\d+\.\d+(?:\.\d+)*", version) or tuple(map(int, version.split("."))) < (2, 34):
            raise ValidationError("Unsupported LOSAT runtime: Linux requires glibc 2.34 or later; musl/Alpine is not supported.")
    return target


def read_release_lock() -> dict:
    try:
        lock = json.loads(resources.files("gbdraw.data").joinpath("losat-release.json").read_text())
        validate_release_lock(lock)
        return lock
    except (OSError, ValueError, KeyError, TypeError) as exc:
        raise ValidationError(f"Invalid packaged LOSAT release lock: {exc}") from exc


def validate_release_lock(lock: dict) -> None:
    """Shared schema boundary for the packaged lock and its release generator."""
    if lock["schema_version"] != 1 or lock["version"] != "0.1.0":
        raise ValueError("unsupported schema or version")
    artifacts = lock["artifacts"]
    if not isinstance(artifacts, dict):
        raise ValueError("invalid artifacts")
    if not artifacts and lock["candidate_sha"] is None:
        return  # No public release has been verified yet.
    if set(artifacts) != set(TARGETS.values()) or not re.fullmatch(r"[0-9a-f]{40}", lock["candidate_sha"]):
        raise ValueError("incomplete native release identity")
    for target, entry in artifacts.items():
        binary = "LOSAT.exe" if target.endswith("windows-msvc") else "LOSAT"
        suffix = ".zip" if binary.endswith(".exe") else ".tar.gz"
        root = f"LOSAT-{lock['version']}-{target}"
        if entry["filename"] != root + suffix or entry["binary"] != binary:
            raise ValueError("invalid artifact name")
        for field in ("sha256", "binary_sha256"):
            if not re.fullmatch(r"[0-9a-f]{64}", entry[field]):
                raise ValueError(f"invalid {field}")
        for field in ("size", "binary_size"):
            if type(entry[field]) is not int or entry[field] <= 0:
                raise ValueError(f"invalid {field}")


def cache_root() -> Path:
    if platform.system() == "Windows":
        base = Path(os.environ.get("LOCALAPPDATA", Path.home() / "AppData/Local"))
    elif platform.system() == "Darwin":
        base = Path.home() / "Library/Caches"
    else:
        base = Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache"))
    return base.expanduser().absolute() / "gbdraw/losat"


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _identity(lock: dict, target: str) -> dict:
    entry = lock["artifacts"][target]
    return {
        "version": lock["version"], "target": target, "candidate_sha": lock["candidate_sha"],
        "url": f"{RELEASE_URL}/v{lock['version']}/{entry['filename']}", **entry,
    }


def _verify_cache(directory: Path, identity: dict) -> Path | None:
    if not directory.exists() and not directory.is_symlink():
        return None
    try:
        binary = directory / identity["binary"]
        receipt = directory / "INSTALL.json"
        if directory.is_symlink() or binary.is_symlink() or receipt.is_symlink():
            raise ValueError("symbolic links are not allowed")
        if json.loads(receipt.read_text()) != identity:
            raise ValueError("installation identity differs from the pinned release")
        if binary.stat().st_size != identity["binary_size"] or file_sha256(binary) != identity["binary_sha256"]:
            raise ValueError("binary checksum mismatch")
        result = subprocess.run([str(binary), "--version"], capture_output=True, text=True, timeout=15, check=False)
        if result.returncode or result.stdout.strip() != f"losat {identity['version']}":
            raise ValueError(f"version check failed: {result.stderr.strip() or result.stdout.strip()}")
        return binary
    except (OSError, ValueError, subprocess.SubprocessError) as exc:
        raise ValidationError(
            f"Invalid LOSAT cache at {directory}: {exc}. "
            "Remove this version/target directory and rerun gbdraw setup-losat."
        ) from exc


def managed_losat() -> Path | None:
    """Resolve only an already installed, verified binary; never use the network."""
    lock = read_release_lock()
    if not lock["artifacts"]:
        return None
    try:
        target = runtime_target()
    except ValidationError:
        return None  # No managed target; existing explicit/bundled/PATH policy applies.
    identity = _identity(lock, target)
    return _verify_cache(cache_root() / lock["version"] / target, identity)


@contextmanager
def _installation_lock(path: Path):
    # OS advisory locks are released even when an installer is killed.
    with path.open("a+b") as handle:
        handle.seek(0)
        try:
            if os.name == "nt":
                import msvcrt
                if path.stat().st_size == 0:
                    handle.write(b"\0")
                    handle.flush()
                    handle.seek(0)
                msvcrt.locking(handle.fileno(), msvcrt.LK_NBLCK, 1)
            else:
                import fcntl
                fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as exc:
            raise ValidationError("Another LOSAT setup is running for this target; retry after it finishes.") from exc
        try:
            yield
        finally:
            if os.name == "nt":
                handle.seek(0)
                msvcrt.locking(handle.fileno(), msvcrt.LK_UNLCK, 1)
            else:
                fcntl.flock(handle, fcntl.LOCK_UN)


def _unpack_binary(archive: Path, destination: Path, identity: dict) -> None:
    root = identity["filename"].removesuffix(".tar.gz").removesuffix(".zip")
    expected = {f"{root}/{name}" for name in (identity["binary"], "LICENSE", "README.md", "RELEASE-METADATA.json")}
    if archive.suffix == ".zip":
        handle = zipfile.ZipFile(archive)
        members = handle.infolist()
        names = [m.filename for m in members]
        regular = all(stat.S_IFMT(m.external_attr >> 16) in (0, stat.S_IFREG) and not m.is_dir() for m in members)
        sizes = {m.filename: m.file_size for m in members}
        open_member = handle.open
    else:
        handle = tarfile.open(archive, "r:gz")
        members = handle.getmembers()
        names = [m.name for m in members]
        regular = all(m.isfile() for m in members)
        sizes = {m.name: m.size for m in members}
        open_member = handle.extractfile
    with handle:
        if set(names) != expected or len(names) != len(expected) or not regular:
            raise ValueError("unexpected archive members or links")
        binary_name = f"{root}/{identity['binary']}"
        if sizes[binary_name] != identity["binary_size"]:
            raise ValueError("unexpected binary size")
        metadata_name = f"{root}/RELEASE-METADATA.json"
        if sizes[metadata_name] > 1024 * 1024:
            raise ValueError("oversized release metadata")
        with open_member(metadata_name) as source:
            metadata = json.load(source)
        if (metadata["candidate_sha"] != identity["candidate_sha"]
                or metadata["release"] != f"v{identity['version']}"
                or metadata["artifact"]["target"] != identity["target"]
                or metadata["binary"]["sha256"] != identity["binary_sha256"]):
            raise ValueError("archive metadata differs from the pinned release")
        # Never extract paths, permissions or links from the archive.
        with open_member(binary_name) as source, (destination / identity["binary"]).open("xb") as output:
            shutil.copyfileobj(source, output)
    (destination / identity["binary"]).chmod(0o755)


def setup_losat() -> Path:
    target = runtime_target()
    lock = read_release_lock()
    if not lock["artifacts"]:
        raise ValidationError("LOSAT v0.1.0 has no verified public release lock yet. Use --losatp_bin with an explicit executable.")
    identity = _identity(lock, target)
    parent = cache_root() / lock["version"]
    directory = parent / target
    try:
        parent.mkdir(parents=True, exist_ok=True)
        with _installation_lock(parent / f"{target}.lock"):
            installed = _verify_cache(directory, identity)
            if installed is not None:
                return installed
            with tempfile.TemporaryDirectory(prefix=f".{target}-", dir=parent) as temporary:
                temporary = Path(temporary)
                archive = temporary / identity["filename"]
                with urlopen(identity["url"], timeout=60) as response, archive.open("xb") as output:
                    if not response.geturl().startswith("https://"):
                        raise ValueError("release download redirected away from HTTPS")
                    remaining = identity["size"]
                    while remaining:
                        chunk = response.read(min(1024 * 1024, remaining))
                        if not chunk:
                            raise ValueError("incomplete archive download")
                        output.write(chunk)
                        remaining -= len(chunk)
                    if response.read(1):
                        raise ValueError("archive exceeds pinned size")
                if file_sha256(archive) != identity["sha256"]:
                    raise ValueError("archive SHA-256 mismatch")
                staged = temporary / "install"
                staged.mkdir()
                _unpack_binary(archive, staged, identity)
                (staged / "INSTALL.json").write_text(json.dumps(identity, indent=2) + "\n")
                _verify_cache(staged, identity)
                os.replace(staged, directory)
                return directory / identity["binary"]
    except (OSError, ValueError, KeyError, TypeError, EOFError, tarfile.TarError, zipfile.BadZipFile) as exc:
        raise ValidationError(f"LOSAT setup failed for {target}: {exc}") from exc


def setup_main(args: list[str]) -> None:
    parser = argparse.ArgumentParser(prog="gbdraw setup-losat", description=__doc__)
    parser.parse_args(args)
    print(setup_losat())
