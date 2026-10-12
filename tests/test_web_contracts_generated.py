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


@pytest.mark.browser
def test_the_contracts_directory_holds_only_the_generated_modules() -> None:
    # R14: `contracts/` holds the modules that `--check` compares and nothing else.
    generator = _generator()
    targets = {generator.TARGET, generator.REQUEST_TARGET, generator.SESSION_TARGET}

    assert set(generator.TARGET.parent.rglob("*")) == targets


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


def test_a_changed_session_field_set_makes_the_generated_module_stale(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import gbdraw.session_io as session_io

    generator = _generator()
    assert (
        "/** @type {readonly string[]} */\n"
        "export const OTHER_MODE_RESULT_LEGEND_FIELDS = Object.freeze(['originalColors', 'originalOrder']);\n"
    ) in generator.render_web_contracts_module()

    monkeypatch.setattr(
        session_io, "OTHER_MODE_RESULT_LEGEND_FIELDS", session_io.OTHER_MODE_RESULT_LEGEND_FIELDS | {"another"}
    )

    assert "Object.freeze(['another', 'originalColors', 'originalOrder'])" in generator.render_web_contracts_module()


def test_the_session_typedefs_import_only_the_request_typedefs() -> None:
    # A leaf of the Web layers may import a leaf (R13).
    source = _generator().render_session_typedefs_module()

    assert re.findall(r"\bimport\b.*|\bfrom\s+['\"].*", source) == [
        "import { CanonicalRenderRequest } from './request.generated.js' */"
    ]
    assert source.endswith("export {};\n")
    assert re.search(r"\bany\b", source.split("\n\n", 1)[1]) is None
    assert (
        "/** @typedef {SessionInput & { format: 'gbdraw-session', version: number }} SessionEnvelope */\n"
    ) in source


def test_a_changed_session_declaration_makes_the_session_typedefs_stale(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import gbdraw.session_contract as contract

    generator = _generator()
    current = generator.render_session_typedefs_module()
    assert " * @property {string[]} [originalOrder]\n" in current
    assert "another" not in current

    monkeypatch.setitem(contract.SessionOtherModeLegend.__annotations__, "another", "str")

    assert " * @property {string} another\n" in generator.render_session_typedefs_module()


def test_a_changed_registry_row_makes_the_session_typedefs_stale(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import gbdraw.session_contract as contract

    generator = _generator()
    assert " * @property {unknown} [canvasPadding]\n" in generator.render_session_typedefs_module()

    monkeypatch.setitem(contract.SLICE_CONTAINERS, "ui", contract.SLICE_CONTAINERS["ui"] | {"another"})

    assert " * @property {unknown} [another]\n" in generator.render_session_typedefs_module()


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


def test_an_open_object_renders_with_an_index_signature() -> None:
    from gbdraw.contract_shapes import INTEGER, STRING, obj

    generator = _generator()

    # `@property` cannot carry an index signature, so the typedef is a type literal.
    assert generator._typedefs({"Open": obj({"id": STRING, "size": obj({"n": INTEGER}, open=True)}, optional=["size"], open=True)}) == (
        "/**\n"
        " * @typedef {{\n"
        " *   id: string,\n"
        " *   size?: { n: number, [key: string]: unknown },\n"
        " *   [key: string]: unknown\n"
        " * }} Open\n"
        " */\n"
    )
    assert generator._jsdoc_type(obj({}, open=True)) == "{ [key: string]: unknown }"
