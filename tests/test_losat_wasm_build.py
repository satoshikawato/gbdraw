"""The shipped LOSAT Wasm files match their recorded build (``losat-build.json``)."""

from __future__ import annotations

import importlib.util
import shutil
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
_SPEC = importlib.util.spec_from_file_location(
    "record_losat_wasm_build", REPO_ROOT / "tools" / "record_losat_wasm_build.py"
)
assert _SPEC is not None and _SPEC.loader is not None
build_tool = importlib.util.module_from_spec(_SPEC)
sys.modules[_SPEC.name] = build_tool
_SPEC.loader.exec_module(build_tool)


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
