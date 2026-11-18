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


def test_the_generated_module_records_each_shipped_file_digest() -> None:
    # The Web app compiles a fetched Wasm file only when its sha256 is the one recorded here (OV-2440).
    manifest = json.loads(build_tool.MANIFEST_PATH.read_text(encoding="utf-8"))
    module = build_tool.MODULE_PATH.read_text(encoding="utf-8")
    for name, record in manifest["files"].items():
        assert record["sha256"] == hashlib.sha256((build_tool.WASM_DIR / name).read_bytes()).hexdigest()
        assert f'"{name}": {{\n      "sha256": "{record["sha256"]}"\n    }}' in module


def test_the_copied_guard_is_the_one_the_manifest_records(tmp_path: Path) -> None:
    # tests/web/losat-wasm-guard.test.mjs proves the threaded file with this copy of LOSAT's guard (OV-2441).
    tool = json.loads(build_tool.MANIFEST_PATH.read_text(encoding="utf-8"))["files"]["losat-threaded.wasm"]["guardTool"]
    assert tool == {"path": build_tool.GUARD_TOOL, "sha256": hashlib.sha256(build_tool.GUARD_TOOL_COPY.read_bytes()).hexdigest()}
    changed = tmp_path / "wasi_shared_memory.js"
    changed.write_bytes(build_tool.GUARD_TOOL_COPY.read_bytes() + b"\n")
    assert build_tool.check(guard_copy=changed) == [f"{changed.name} is not the guard the manifest records (sha256 {tool['sha256']})"]


def test_the_check_refuses_a_threaded_file_without_a_recorded_guard(tmp_path: Path) -> None:
    wasm_dir = tmp_path / "losat"
    shutil.copytree(build_tool.WASM_DIR, wasm_dir)
    manifest = json.loads((wasm_dir / "losat-build.json").read_text(encoding="utf-8"))
    del manifest["files"]["losat-threaded.wasm"]["guardTool"]
    (wasm_dir / "losat-build.json").write_text(json.dumps(manifest), encoding="utf-8")
    assert build_tool.check(wasm_dir) == ["the manifest records no guard of losat-threaded.wasm"]


@pytest.mark.skipif(shutil.which("node") is None, reason="LOSAT's guard runs on Node.js")
def test_recording_guards_the_threaded_file_with_the_copied_guard(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    # record() writes cargo's threaded output through LOSAT's guard, keeps the record the guard prints and the copy
    # of the guard that ran. The guard returns a module it already guarded as it is, so the shipped file stands in
    # for cargo's output here.
    crate = tmp_path / "LOSAT"
    (crate / "tests").mkdir(parents=True)
    shutil.copyfile(build_tool.GUARD_TOOL_COPY, crate / build_tool.GUARD_TOOL)
    copy = tmp_path / "copy" / "wasi_shared_memory.js"
    monkeypatch.setattr(build_tool, "GUARD_TOOL_COPY", copy)
    destination = tmp_path / "losat-threaded.wasm"
    shipped = build_tool.WASM_DIR / "losat-threaded.wasm"
    guard, tool = build_tool.guard_threaded_file(crate, shipped, destination)
    recorded = json.loads(build_tool.MANIFEST_PATH.read_text(encoding="utf-8"))["files"]["losat-threaded.wasm"]
    assert (guard, tool) == (recorded["sharedMemoryGuard"], recorded["guardTool"])
    assert destination.read_bytes() == shipped.read_bytes()
    assert copy.read_bytes() == build_tool.GUARD_TOOL_COPY.read_bytes()


def _forged_threaded_file() -> bytes:
    """The shipped threaded file's header, memory import and guard record, without its code: a file whose guard
    record claims the shipped counts although nothing in it is guarded."""

    shipped = (build_tool.WASM_DIR / "losat-threaded.wasm").read_bytes()
    kept = []
    for section_id, body, end in build_tool._sections(shipped):
        length, cursor = build_tool._leb128(shipped, body)
        named = section_id == 0 and shipped[cursor : cursor + length] == b"losat.shared_memory_guard"
        if section_id == 2 or named:
            kept.append(bytes([section_id]) + _leb(end - body) + shipped[body:end])
    assert len(kept) == 2
    return shipped[:8] + b"".join(kept)


def _leb(value: int) -> bytes:
    out = bytearray()
    while True:
        byte = value & 0x7F
        value >>= 7
        out.append(byte | (0x80 if value else 0))
        if not value:
            return bytes(out)


def test_the_check_refuses_a_threaded_file_whose_guard_record_is_not_its_code(tmp_path: Path) -> None:
    # OV-2441: the record alone (a custom section) passed the check; LOSAT's verifier walks the code and refuses it.
    wasm_dir = tmp_path / "losat"
    shutil.copytree(build_tool.WASM_DIR, wasm_dir)
    forged = _forged_threaded_file()
    (wasm_dir / "losat-threaded.wasm").write_bytes(forged)
    manifest = json.loads((wasm_dir / "losat-build.json").read_text(encoding="utf-8"))
    manifest["files"]["losat-threaded.wasm"].update(sha256=hashlib.sha256(forged).hexdigest(), size=len(forged))
    (wasm_dir / "losat-build.json").write_text(json.dumps(manifest), encoding="utf-8")
    module = tmp_path / "module.js"
    module.write_text(build_tool.render_module(manifest), encoding="utf-8")
    problems = build_tool.check(wasm_dir, module)
    assert len(problems) == 1 and problems[0].startswith("LOSAT's guard refuses losat-threaded.wasm:"), problems
