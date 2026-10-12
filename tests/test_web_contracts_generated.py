from __future__ import annotations

import dataclasses
import functools
import importlib
import importlib.util
import re
import subprocess
import sys
from pathlib import Path
from types import ModuleType
from typing import Literal, Optional, get_type_hints

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
GENERATOR = REPO_ROOT / "tools" / "generate_web_contracts.py"


@functools.cache
def _generator() -> ModuleType:
    spec = importlib.util.spec_from_file_location("generate_web_contracts", GENERATOR)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# The `browser` marker puts this test in the jobs a Web-only pull request runs.
@pytest.mark.browser
def test_generated_web_contracts_match_python_source() -> None:
    result = subprocess.run(
        [sys.executable, "tools/generate_web_contracts.py", "--check"],
        cwd=REPO_ROOT,
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.parametrize(
    ("web_name", "module", "name"),
    [
        ("SESSION_VERSION", "gbdraw.session_io", "CURRENT_SESSION_VERSION"),
        ("CANONICAL_SESSION_MIN_VERSION", "gbdraw.session_io", "CANONICAL_SESSION_MIN_VERSION"),
        ("DISPLAY_PLACEMENT_SCHEMA", "gbdraw.session_request_codec", "DISPLAY_PLACEMENT_SCHEMA"),
        ("SIMILARITY_ALIGNMENT_PLAN_SCHEMA", "gbdraw.layout.similarity_alignment", "SIMILARITY_ALIGNMENT_PLAN_SCHEMA"),
        ("LEGACY_PROTEIN_LOSAT_CACHE_SCHEMA", "gbdraw.analysis.protein_colinearity", "LEGACY_PROTEIN_LOSAT_CACHE_SCHEMA"),
        ("ORTHOGROUP_INFERENCE_VERSION", "gbdraw.analysis.protein_colinearity", "ORTHOGROUP_INFERENCE_VERSION"),
    ],
)
def test_a_changed_python_constant_makes_the_generated_module_stale(
    monkeypatch: pytest.MonkeyPatch, web_name: str, module: str, name: str
) -> None:
    owner = importlib.import_module(module)
    value = getattr(owner, name)
    generator = _generator()
    current = generator.render_web_contracts_module()
    assert f"export const {web_name} = {generator._scalar(value)};\n" in current

    monkeypatch.setattr(owner, name, value + (1 if isinstance(value, int) else "-next"))

    assert generator.render_web_contracts_module() != current


@pytest.mark.parametrize("render", ["render_web_contracts_module", "render_request_typedefs_module"])
def test_generated_web_contracts_module_imports_nothing(render: str) -> None:
    # `contracts/` is a leaf directory of the Web layers (R13).
    source = getattr(_generator(), render)()

    assert re.search(r"\bimport\b|\bfrom\s+['\"]", source) is None


def test_a_changed_codec_table_makes_the_request_typedefs_stale(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import gbdraw.session_request_codec as codec

    generator = _generator()
    current = generator.render_request_typedefs_module()
    assert " * @property {boolean} [collinearInferOrthogroups]\n" in current
    assert current.endswith("export {};\n")
    assert re.search(r"\bany\b", current) is None

    monkeypatch.setattr(codec, "_OPTIONAL_PROTEIN_SETTINGS", frozenset())

    assert " * @property {boolean} collinearInferOrthogroups\n" in generator.render_request_typedefs_module()


@dataclasses.dataclass(frozen=True)
class _Inner:
    value: float


@dataclasses.dataclass(frozen=True)
class _Outer:
    name: str
    count: int | None
    mode: Literal["a", "b"]
    items: list[str]
    pairs: tuple[tuple[int, bool], ...]
    lookup: dict[str, float]
    inner: Optional[_Inner]
    anything: object


def test_type_hints_render_as_jsdoc_types() -> None:
    from gbdraw.contract_shapes import obj, ref, shape_of_hint

    def special(hint: object) -> dict[str, object] | None:
        return ref("Inner") if hint is _Inner else None

    hints = get_type_hints(_Outer)
    rendered = {
        name: _generator()._jsdoc_type(shape_of_hint(hint, special)) for name, hint in hints.items()
    }

    assert rendered == {
        "name": "string",
        "count": "number | null",
        "mode": "'a' | 'b'",
        "items": "string[]",
        "pairs": "[number, boolean][]",
        "lookup": "Record<string, number>",
        "inner": "Inner | null",
        "anything": "unknown",
    }
    assert _generator()._jsdoc_type(shape_of_hint(list[str | None])) == "(string | null)[]"
    assert _generator()._jsdoc_type(obj({})) == "Record<string, never>"
    with pytest.raises(TypeError, match="No JSON shape"):
        shape_of_hint(_Inner)


def test_named_shapes_render_as_typedefs() -> None:
    from gbdraw.contract_shapes import STRING, INTEGER, literal, nullable, obj, union

    typedefs = _generator()._typedefs({
        "Thing": obj(
            {"id": STRING, "size": nullable(INTEGER), "kind": obj({"tag": literal("x", 1, False, None)})},
            optional=["size"],
        ),
        "Either": union(STRING, obj({"a": INTEGER}, optional=["a"])),
    })

    assert typedefs == (
        "/**\n"
        " * @typedef {object} Thing\n"
        " * @property {string} id\n"
        " * @property {number | null} [size]\n"
        " * @property {{ tag: 'x' | 1 | false | null }} kind\n"
        " */\n"
        "/** @typedef {string | { a?: number }} Either */\n"
    )
