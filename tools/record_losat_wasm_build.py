#!/usr/bin/env python3
"""Record, or check, the LOSAT build of the Web app's two Wasm files.

Record (after building both artifacts from one LOSAT checkout, see
``gbdraw/web/wasm/losat/README.md``)::

    python tools/record_losat_wasm_build.py --losat-src <LOSAT checkout> --target-dir <CARGO_TARGET_DIR>

copies the two artifacts into ``gbdraw/web/wasm/losat/``, writes the build
manifest ``losat-build.json`` beside them (LOSAT commit, build commands,
toolchain, each file's sha256 and size), and writes the build the Web app runs
into ``gbdraw/web/js/services/losat-build.generated.js``: the version it records
with each search and each file's sha256, which the app checks before it compiles
a file. The threaded file goes through LOSAT's shared-memory guard on the way:
cargo's output is not guarded, and the Web app does not guard a module itself.
The guard of the recorded commit is copied, byte for byte, to
``tools/losat/wasi_shared_memory.js``, so that the check proves the shipped
threaded file with LOSAT's own verifier.

``--check`` fails when a Wasm file, the copied guard or the generated module does
not match the manifest, or when LOSAT's verifier (``node
tools/losat/wasi_shared_memory.js verify``) refuses the threaded file or counts
other guarded instructions than the manifest records.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
from collections.abc import Iterator
from pathlib import Path
from typing import TypedDict

if sys.version_info >= (3, 11):
    import tomllib
else:  # pragma: no cover - Python 3.10
    import tomli as tomllib

REPO_ROOT = Path(__file__).resolve().parents[1]
WASM_DIR = REPO_ROOT / "gbdraw" / "web" / "wasm" / "losat"
MANIFEST_PATH = WASM_DIR / "losat-build.json"
MODULE_PATH = REPO_ROOT / "gbdraw" / "web" / "js" / "services" / "losat-build.generated.js"
REPOSITORY = "https://github.com/satoshikawato/LOSAT"
MANIFEST_SCHEMA = 1
GUARDED_ARTIFACT = "losat-threaded.wasm"
GUARD_TOOL = "tests/wasi_shared_memory.js"  # in LOSAT's crate
GUARD_TOOL_COPY = REPO_ROOT / "tools" / "losat" / "wasi_shared_memory.js"
MINIMUM_GUARD_VERSION = 3  # version 3 guards memory.grow and atomic instructions as well as memory.copy and memory.fill


class ArtifactSpec(TypedDict):
    target: str
    cargo_args: list[str]
    built: str
    shape: str


# Published file name -> how LOSAT builds it (its `.cargo/config.toml` aliases
# `build-web-serial` and `build-web-threaded`, without their per-alias target
# directories). The serial file is the library's cdylib, a WASI reactor with the
# direct `losat_web_*` API (docs/web/abi_v1.md); the threaded file is the
# `losat` command for wasm32-wasip1-threads.
ARTIFACTS: dict[str, ArtifactSpec] = {
    "losat.wasm": {
        "target": "wasm32-wasip1",
        "cargo_args": ["--lib", "--no-default-features"],
        "built": "LOSAT.wasm",
        "shape": "serial reactor, direct losat_web_* API (_initialize)",
    },
    "losat-threaded.wasm": {
        "target": "wasm32-wasip1-threads",
        "cargo_args": ["--bin", "losat", "--features", "wasm-threads"],
        "built": "losat.wasm",
        "shape": "threaded command (_start, wasi_thread_start; imports env.memory, wasi.thread-spawn)",
    },
}


class MemoryPages(TypedDict):
    initial: int
    maximum: int


class SharedMemoryGuard(TypedDict):
    """The record LOSAT's guard writes into the module it guarded (its custom section) and prints."""

    version: int
    fill: int
    copy: int
    grow: int
    atomic: int
    atomic_helpers: int


class GuardTool(TypedDict):
    """The guard that guarded the threaded file: its path in LOSAT's crate and its sha256, which is the sha256 of
    ``tools/losat/wasi_shared_memory.js``."""

    path: str
    sha256: str


class _FileRecordRequired(TypedDict):
    sha256: str
    size: int
    target: str
    shape: str
    command: str
    rustflags: list[str]


class FileRecord(_FileRecordRequired, total=False):
    memory: MemoryPages
    sharedMemoryGuard: SharedMemoryGuard
    guardTool: GuardTool


class BuildManifest(TypedDict):
    schema: int
    repository: str
    commit: str
    commitDate: str
    cargoVersion: str
    version: str
    toolchain: dict[str, str]
    files: dict[str, FileRecord]


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _run(argv: list[str], cwd: Path) -> str:
    return subprocess.run(argv, cwd=cwd, check=True, capture_output=True, text=True).stdout.strip()


def _leb128(data: bytes, offset: int) -> tuple[int, int]:
    value = shift = 0
    while True:
        byte = data[offset]
        offset += 1
        value |= (byte & 0x7F) << shift
        shift += 7
        if not byte & 0x80:
            return value, offset


def _limits(data: bytes, offset: int) -> tuple[int, int | None, int]:
    flags = data[offset]
    initial, offset = _leb128(data, offset + 1)
    maximum = None
    if flags & 0x01:
        maximum, offset = _leb128(data, offset)
    return initial, maximum, offset


def _sections(wasm: bytes) -> Iterator[tuple[int, int, int]]:
    """Each section of a Wasm module: its id and where its body starts and ends."""

    if wasm[:8] != b"\0asm\x01\0\0\0":
        raise ValueError("not a Wasm module")
    offset = 8
    while offset < len(wasm):
        section_id = wasm[offset]
        size, body = _leb128(wasm, offset + 1)
        offset = body + size
        yield section_id, body, offset


def imported_memory_pages(wasm: bytes) -> MemoryPages | None:
    """The page limits of the memory a Wasm module imports, or None when it imports none."""

    for section_id, body, _end in _sections(wasm):
        if section_id != 2:
            continue
        count, cursor = _leb128(wasm, body)
        for _ in range(count):
            for _name in range(2):
                length, cursor = _leb128(wasm, cursor)
                cursor += length
            kind = wasm[cursor]
            cursor += 1
            if kind == 0x02:
                initial, maximum, cursor = _limits(wasm, cursor)
                if maximum is None:
                    raise ValueError("the imported memory declares no maximum")
                return {"initial": initial, "maximum": maximum}
            if kind == 0x00:
                _index, cursor = _leb128(wasm, cursor)
            elif kind == 0x01:
                _initial, _maximum, cursor = _limits(wasm, cursor + 1)
            elif kind == 0x03:
                cursor += 2
            else:
                raise ValueError(f"unsupported import kind {kind}")
        return None
    return None


def build_version(cargo_version: str, commit: str) -> str:
    """The version the Web records with each search: the Cargo version and the commit."""

    return f"{cargo_version}+{commit[:12]}"


def render_module(manifest: BuildManifest) -> str:
    identity = {
        "commit": manifest["commit"],
        "version": manifest["version"],
        "files": {name: {"sha256": record["sha256"]} for name, record in sorted(manifest["files"].items())},
        "threadedMemoryPages": manifest["files"]["losat-threaded.wasm"]["memory"],
    }
    return (
        "// @ts-check\n"
        "// Generated by tools/record_losat_wasm_build.py from gbdraw/web/wasm/losat/losat-build.json;\n"
        "// do not edit manually.\n"
        "\n"
        "/**\n"
        " * The LOSAT build of the shipped Wasm files: the commit, the version the Web records with each\n"
        " * search (the Cargo version and the commit), each file's sha256, which the app checks before it\n"
        " * compiles the file, and the page limits of the memory the threaded file imports.\n"
        " * @type {Readonly<{ commit: string, version: string, files: Readonly<Record<LosatWasmFile, { sha256: string }>>, threadedMemoryPages: { initial: number, maximum: number } }>}\n"
        " */\n"
        f"export const LOSAT_WASM_BUILD = Object.freeze({json.dumps(identity, indent=2)});\n"
        "\n"
        f"/** @typedef {{{' | '.join(repr(name) for name in sorted(manifest['files']))}}} LosatWasmFile A shipped LOSAT Wasm file. */\n"
    )


def _guard_record(printed: dict[str, int]) -> SharedMemoryGuard:
    return SharedMemoryGuard(
        version=printed["version"],
        fill=printed["fill"],
        copy=printed["copy"],
        grow=printed["grow"],
        atomic=printed["atomic"],
        atomic_helpers=printed["atomic_helpers"],
    )


def verified_guard(wasm: Path, guard_copy: Path = GUARD_TOOL_COPY) -> SharedMemoryGuard | str:
    """The guard record of ``wasm`` as LOSAT's verifier counts it in the code, or why the verifier refuses it."""

    result = subprocess.run(
        ["node", str(guard_copy), "verify", str(wasm)], capture_output=True, text=True, check=False
    )
    if result.returncode != 0:
        return result.stderr.strip() or f"exit {result.returncode}"
    return _guard_record(json.loads(result.stdout))


def guard_threaded_file(crate: Path, built: Path, destination: Path) -> tuple[SharedMemoryGuard, GuardTool]:
    """Write ``built`` through LOSAT's guard to ``destination`` and copy the guard to ``tools/losat/``.

    The guard verifies what it writes and prints its record; this keeps the record, refuses a guard older than
    ``MINIMUM_GUARD_VERSION`` or one that guarded nothing, and requires that the copy is the guard that ran.
    """

    GUARD_TOOL_COPY.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(crate / GUARD_TOOL, GUARD_TOOL_COPY)
    printed = json.loads(_run(["node", GUARD_TOOL, "guard", str(built), str(destination)], crate))
    guard = _guard_record(printed)
    tool = GuardTool(path=GUARD_TOOL, sha256=_sha256(GUARD_TOOL_COPY))
    if printed["guard_sha256"] != tool["sha256"]:
        raise SystemExit(f"{GUARD_TOOL_COPY} is not the guard that ran (sha256 {printed['guard_sha256']})")
    if printed["output_sha256"] != _sha256(destination):
        raise SystemExit(f"{destination} is not the module the guard wrote")
    if guard["version"] < MINIMUM_GUARD_VERSION:
        raise SystemExit(f"guard version {guard['version']} is older than {MINIMUM_GUARD_VERSION}")
    if guard["fill"] + guard["copy"] + guard["grow"] + guard["atomic"] == 0:
        raise SystemExit(f"the guard guarded nothing in {built}; it would run unguarded in the Web app")
    return guard, tool


def record(losat_src: Path, target_dir: Path) -> BuildManifest:
    crate = losat_src / "LOSAT"
    if _run(["git", "status", "--porcelain"], losat_src):
        raise SystemExit(f"{losat_src} has uncommitted changes; build from a clean checkout")
    commit = _run(["git", "rev-parse", "HEAD"], losat_src)
    cargo_version = tomllib.loads((crate / "Cargo.toml").read_text())["package"]["version"]
    config = tomllib.loads((crate / ".cargo" / "config.toml").read_text())
    files: dict[str, FileRecord] = {}
    for name, spec in ARTIFACTS.items():
        built = target_dir / spec["target"] / "release" / spec["built"]
        destination = WASM_DIR / name
        command = ["cargo", "build", "--release", "--locked", "--target", spec["target"], *spec["cargo_args"]]
        guarded: tuple[SharedMemoryGuard, GuardTool] | None = None
        if name == GUARDED_ARTIFACT:
            guarded = guard_threaded_file(crate, built, destination)
            command += ["&&", "node", GUARD_TOOL, "guard", spec["built"], name]
        else:
            shutil.copyfile(built, destination)
        files[name] = {
            "sha256": _sha256(destination),
            "size": destination.stat().st_size,
            "target": spec["target"],
            "shape": spec["shape"],
            "command": " ".join(command),
            "rustflags": list(config["target"][spec["target"]]["rustflags"]),
        }
        wasm = destination.read_bytes()
        memory = imported_memory_pages(wasm)
        if memory is not None:
            files[name]["memory"] = memory
        if guarded is not None:
            files[name]["sharedMemoryGuard"], files[name]["guardTool"] = guarded
    manifest: BuildManifest = {
        "schema": MANIFEST_SCHEMA,
        "repository": REPOSITORY,
        "commit": commit,
        "commitDate": _run(["git", "log", "-1", "--format=%cI"], losat_src),
        "cargoVersion": cargo_version,
        "version": build_version(cargo_version, commit),
        "toolchain": {
            "rustc": _run(["rustc", "-vV"], crate),
            "cargo": _run(["cargo", "-V"], crate),
        },
        "files": files,
    }
    MANIFEST_PATH.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    MODULE_PATH.write_text(render_module(manifest), encoding="utf-8")
    return manifest


def check(wasm_dir: Path = WASM_DIR, module_path: Path = MODULE_PATH, guard_copy: Path = GUARD_TOOL_COPY) -> list[str]:
    """The differences between the Wasm files, the copied guard, the generated module and the manifest.

    Whether the threaded file is guarded is for LOSAT's verifier to say (the copied guard's ``verify``, which walks
    the code): the record in the file's custom section alone proves nothing (OV-2441).
    """

    manifest: BuildManifest = json.loads((wasm_dir / MANIFEST_PATH.name).read_text(encoding="utf-8"))
    problems = []
    if sorted(manifest["files"]) != sorted(ARTIFACTS):
        problems.append(f"the manifest lists {sorted(manifest['files'])}, expected {sorted(ARTIFACTS)}")
    if manifest["version"] != build_version(manifest["cargoVersion"], manifest["commit"]):
        problems.append("the manifest version is not its Cargo version and commit")
    for name, expected in manifest["files"].items():
        path = wasm_dir / name
        if not path.is_file():
            problems.append(f"{name} is missing")
            continue
        size, digest = path.stat().st_size, _sha256(path)
        wasm = path.read_bytes()
        if name == GUARDED_ARTIFACT:
            tool = expected.get("guardTool")
            if tool is None or "sharedMemoryGuard" not in expected:
                problems.append(f"the manifest records no guard of {name}")
            elif not guard_copy.is_file() or _sha256(guard_copy) != tool["sha256"]:
                problems.append(f"{guard_copy.name} is not the guard the manifest records (sha256 {tool['sha256']})")
            else:
                verified = verified_guard(path, guard_copy)
                if isinstance(verified, str):
                    problems.append(f"LOSAT's guard refuses {name}: {verified}")
                elif verified != expected["sharedMemoryGuard"]:
                    problems.append(
                        f"LOSAT's guard counts {verified} in {name}; the manifest records {expected['sharedMemoryGuard']}"
                    )
        memory = imported_memory_pages(wasm)
        if memory != expected.get("memory"):
            problems.append(f"{name} imports memory {memory}; the manifest records {expected.get('memory')}")
        if (size, digest) != (expected["size"], expected["sha256"]):
            problems.append(
                f"{name} is {size} bytes, sha256 {digest}; the manifest records "
                f"{expected['size']} bytes, sha256 {expected['sha256']}"
            )
    if not module_path.is_file() or module_path.read_text(encoding="utf-8") != render_module(manifest):
        problems.append(f"{module_path.name} does not match the manifest; rerun this tool")
    return problems


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--check", action="store_true", help="check the files against the manifest")
    parser.add_argument("--losat-src", type=Path, help="the LOSAT checkout the artifacts were built from")
    parser.add_argument("--target-dir", type=Path, help="the CARGO_TARGET_DIR of that build")
    args = parser.parse_args(argv)
    if args.check:
        problems = check()
        for problem in problems:
            print(problem, file=sys.stderr)
        return 1 if problems else 0
    if args.losat_src is None or args.target_dir is None:
        parser.error("--losat-src and --target-dir are required to record a build")
    manifest = record(args.losat_src.resolve(), args.target_dir.resolve())
    print(f"recorded LOSAT {manifest['version']} ({manifest['commit']})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
