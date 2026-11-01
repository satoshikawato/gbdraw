"""The shipped LOSAT Wasm files match their recorded build (``losat-build.json``)."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import shutil
import sys
import types
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
_SPEC = importlib.util.spec_from_file_location(
    "record_losat_wasm_build", REPO_ROOT / "tools" / "record_losat_wasm_build.py"
)
assert _SPEC is not None and _SPEC.loader is not None
build_tool = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = build_tool
_SPEC.loader.exec_module(build_tool)


def test_the_tool_loads_on_python_3_10_with_tomli(monkeypatch: pytest.MonkeyPatch) -> None:
    # Python 3.10 has no tomllib; the tool must fall back to tomli, which
    # pyproject.toml installs there (OV-2406).
    tomli = types.ModuleType("tomli")
    monkeypatch.setitem(sys.modules, "tomli", tomli)
    monkeypatch.setitem(sys.modules, "tomllib", None)
    monkeypatch.setattr(sys, "version_info", (3, 10, 0, "final", 0))
    spec = importlib.util.spec_from_file_location("record_losat_wasm_build_py310", build_tool.__file__)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert module.tomllib is tomli


def test_shipped_wasm_files_match_the_manifest() -> None:
    assert build_tool.check() == []


def test_a_changed_wasm_file_fails_the_check(tmp_path: Path) -> None:
    wasm_dir = tmp_path / "losat"
    shutil.copytree(build_tool.WASM_DIR, wasm_dir)
    with (wasm_dir / "losat.wasm").open("ab") as handle:
        handle.write(b"\0")
    problems = build_tool.check(wasm_dir)
    assert len(problems) == 1 and problems[0].startswith("losat.wasm is ")


def test_the_generated_module_follows_the_manifest(tmp_path: Path) -> None:
    module = tmp_path / "losat-build.generated.js"
    module.write_text(build_tool.MODULE_PATH.read_text(encoding="utf-8").replace("+", "-", 1), encoding="utf-8")
    assert build_tool.check(module_path=module) == [f"{module.name} does not match the manifest; rerun this tool"]


def test_the_threaded_memory_limits_are_the_ones_the_file_imports() -> None:
    threaded = (build_tool.WASM_DIR / "losat-threaded.wasm").read_bytes()
    pages = build_tool.imported_memory_pages(threaded)
    assert pages is not None and pages["initial"] > 0 and pages["maximum"] >= pages["initial"]
    assert build_tool.imported_memory_pages((build_tool.WASM_DIR / "losat.wasm").read_bytes()) is None
    module = build_tool.MODULE_PATH.read_text(encoding="utf-8")
    assert f'"initial": {pages["initial"]},' in module and f'"maximum": {pages["maximum"]}' in module


def _leb(value: int) -> bytes:
    out = bytearray()
    while True:
        byte = value & 0x7F
        value >>= 7
        out.append(byte | (0x80 if value else 0))
        if not value:
            return bytes(out)


def _custom_section(name: str, payload: bytes) -> bytes:
    body = _leb(len(name)) + name.encode() + payload
    return b"\0" + _leb(len(body)) + body


_GUARD = {"version": 2, "fill": 1, "copy": 2, "atomic": 3, "atomic_helpers": 1}
_HEADER = b"\0asm\x01\0\0\0"


def _guarded(record: object = _GUARD) -> bytes:
    return _HEADER + _custom_section("losat.shared_memory_guard", json.dumps(record).encode())


def test_the_guard_record_is_the_one_losat_wrote() -> None:
    assert build_tool.shared_memory_guard(_guarded()) == _GUARD
    assert build_tool.shared_memory_guard(_HEADER) is None
    assert build_tool.shared_memory_guard(_HEADER + _custom_section("name", b"x")) is None
    assert build_tool.shared_memory_guard(_HEADER + _custom_section("losat.shared_memory_guard", b"{")) is None


@pytest.mark.parametrize(
    ("wasm", "expected"),
    [
        (_guarded(), None),
        (_HEADER, "has no shared-memory guard record"),
        (_guarded({**_GUARD, "version": 1}), "guard version 1 is older than 2"),
        (_guarded({**_GUARD, "fill": 0, "copy": 0, "atomic": 0}), "guard record guards nothing"),
    ],
)
def test_a_threaded_file_needs_a_current_guard(wasm: bytes, expected: str | None) -> None:
    problem = build_tool.guard_problem(wasm)
    assert problem == expected


def test_the_check_refuses_an_unguarded_threaded_file(tmp_path: Path) -> None:
    wasm_dir = tmp_path / "losat"
    wasm_dir.mkdir()
    threaded = _HEADER
    (wasm_dir / "losat.wasm").write_bytes(_HEADER)
    (wasm_dir / "losat-threaded.wasm").write_bytes(threaded)
    record = {"sha256": hashlib.sha256(_HEADER).hexdigest(), "size": len(_HEADER)}
    manifest = {
        "commit": "a" * 40,
        "cargoVersion": "0.2.0",
        "version": build_tool.build_version("0.2.0", "a" * 40),
        "files": {"losat.wasm": record, "losat-threaded.wasm": {**record, "memory": None}},
    }
    (wasm_dir / "losat-build.json").write_text(json.dumps(manifest), encoding="utf-8")
    module = tmp_path / "module.js"
    module.write_text(build_tool.render_module(manifest), encoding="utf-8")
    assert build_tool.check(wasm_dir, module) == ["losat-threaded.wasm has no shared-memory guard record"]


def test_the_check_compares_the_guard_record_with_the_manifest(tmp_path: Path) -> None:
    wasm_dir = tmp_path / "losat"
    wasm_dir.mkdir()
    threaded = _guarded()
    (wasm_dir / "losat.wasm").write_bytes(_HEADER)
    (wasm_dir / "losat-threaded.wasm").write_bytes(threaded)
    manifest = {
        "commit": "a" * 40,
        "cargoVersion": "0.2.0",
        "version": build_tool.build_version("0.2.0", "a" * 40),
        "files": {
            "losat.wasm": {"sha256": hashlib.sha256(_HEADER).hexdigest(), "size": len(_HEADER)},
            "losat-threaded.wasm": {
                "sha256": hashlib.sha256(threaded).hexdigest(),
                "size": len(threaded),
                "memory": None,
                "sharedMemoryGuard": {**_GUARD, "copy": 99},
            },
        },
    }
    (wasm_dir / "losat-build.json").write_text(json.dumps(manifest), encoding="utf-8")
    module = tmp_path / "module.js"
    module.write_text(build_tool.render_module(manifest), encoding="utf-8")
    problems = build_tool.check(wasm_dir, module)
    assert len(problems) == 1 and problems[0].startswith("losat-threaded.wasm carries guard record")
