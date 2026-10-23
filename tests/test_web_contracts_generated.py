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
    targets = {
        generator.TARGET, generator.REQUEST_TARGET, generator.SESSION_TARGET, generator.REQUEST_SHAPE_TARGET,
        generator.SESSION_SHAPE_TARGET,
    }

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


def test_the_member_role_and_confidence_values_come_from_the_python_literals(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # OV-1390: the Web shows these values in the UI language; it lists no value itself.
    from gbdraw.analysis import protein_colinearity

    generator = _generator()
    current = generator.render_web_contracts_module()
    assert (
        "export const ORTHOGROUP_MEMBER_ROLE = Object.freeze({\n"
        "  values: /** @type {readonly ('anchor' | 'coortholog' | 'inparalog' | 'low_confidence' | 'local_paralog')[]} */ "
        "(Object.freeze(['anchor', 'coortholog', 'inparalog', 'low_confidence', 'local_paralog'])),\n"
        "  default: 'anchor'\n"
        "});\n"
    ) in current
    assert (
        "export const ORTHOGROUP_MEMBER_CONFIDENCE = Object.freeze({\n"
        "  values: /** @type {readonly ('high' | 'medium' | 'low')[]} */ (Object.freeze(['high', 'medium', 'low'])),\n"
        "  default: 'high'\n"
        "});\n"
    ) in current

    @dataclasses.dataclass(frozen=True)
    class NextMember:
        role: Literal["anchor", "coortholog", "inparalog", "low_confidence", "local_paralog", "next"] = "anchor"
        confidence: Literal["high", "medium", "low"] = "medium"

    monkeypatch.setattr(protein_colinearity, "OrthogroupMember", NextMember)
    following = generator.render_web_contracts_module()

    assert "'local_paralog', 'next']" in following
    assert "(Object.freeze(['high', 'medium', 'low'])),\n  default: 'medium'\n" in following


def test_a_member_field_default_outside_its_literal_is_refused() -> None:
    @dataclasses.dataclass(frozen=True)
    class Broken:
        role: Literal["anchor"] = "other"  # type: ignore[assignment]

    with pytest.raises(TypeError, match="Broken.role"):
        _generator()._literal_field_declaration("BROKEN", Broken, "role")


def test_the_legend_positions_of_a_mode_have_one_python_owner(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # OV-652: the output options check the mapping the Web reads.
    import gbdraw.api.options as options
    from gbdraw.exceptions import ValidationError

    generator = _generator()
    assert (
        "// gbdraw/api/options.py: LEGEND_POSITIONS\n"
        "/** @type {Readonly<Record<'circular' | 'linear', readonly string[]>>} */\n"
        "export const LEGEND_POSITIONS = Object.freeze({\n"
        "  circular: Object.freeze(['left', 'right', 'top', 'bottom', 'upper_left', 'upper_right',"
        " 'lower_left', 'lower_right', 'none']),\n"
        "  linear: Object.freeze(['left', 'right', 'top', 'bottom', 'none'])\n"
        "});\n"
    ) in generator.render_web_contracts_module()
    with pytest.raises(ValidationError) as circular:
        options.CircularOutputOptions(legend="center")
    with pytest.raises(ValidationError) as linear:
        options.LinearOutputOptions(legend="upper_left")
    assert [str(circular.value), str(linear.value)] == [
        "Circular legend must be one of: left, right, top, bottom, "
        "upper_left, upper_right, lower_left, lower_right, none.",
        "Linear legend must be one of: left, right, top, bottom, none.",
    ]
    monkeypatch.setattr(
        options,
        "LEGEND_POSITIONS",
        {**options.LEGEND_POSITIONS, "linear": (*options.LEGEND_POSITIONS["linear"], "upper_left")},
    )

    assert options.LinearOutputOptions(legend="upper_left").legend == "upper_left"
    assert (
        "  linear: Object.freeze(['left', 'right', 'top', 'bottom', 'none', 'upper_left'])\n"
    ) in generator.render_web_contracts_module()


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


def test_the_request_shape_module_is_the_description_as_frozen_data() -> None:
    from gbdraw.session_request_codec import describe_canonical_request

    source = _generator().render_request_shape_module()

    # A leaf imports a leaf (R13): only the type the checker declares.
    assert re.findall(r"\bimport\b.*|\bfrom\s+['\"].*", source) == [
        "import { ContractDescription } from '../utils/contract-shape.js' */"
    ]
    assert re.search(r"\bany\b", source) is None
    assert "defaults: Object.freeze({\n        collinearInferOrthogroups: true,\n        collinearityParams: null," in source
    assert describe_canonical_request()["types"]["CanonicalProteinComparisonSettings"]["defaults"] == {
        "collinearInferOrthogroups": True,
        # A required field that holds null takes null (OV-777).
        **{key: None for key in (
            "collinearityParams", "losatpBin", "losatpThreads", "ncbiBlastpBin", "orthogroupMemberMaxHits",
            "proteinBlastpCandidateLimit",
        )},
    }


def test_a_changed_option_default_makes_the_request_shape_module_stale(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import gbdraw.session_request_codec as codec
    from gbdraw.api.options import LinearDiagramOptions

    generator = _generator()
    assert 'pairwiseMatchStyle: "ribbon"' in generator.render_request_shape_module()

    monkeypatch.setattr(codec, "_DEFAULT_LINEAR_OPTIONS", LinearDiagramOptions(pairwise_match_style="curve"))

    assert 'pairwiseMatchStyle: "curve"' in generator.render_request_shape_module()


def test_description_data_renders_any_json_key_but_the_prototype() -> None:
    generator = _generator()

    assert generator._frozen_json({"a.b": [1, None, True, 0.5], "c": {}}) == (
        'Object.freeze({"a.b": Object.freeze([1, null, true, 0.5]), c: Object.freeze({})})'
    )
    with pytest.raises(TypeError, match="__proto__"):
        generator._frozen_json({"__proto__": 1})
    with pytest.raises(ValueError):
        generator._frozen_json(float("nan"))


def test_the_session_shape_module_is_the_admission_description_as_frozen_data() -> None:
    source = _generator().render_session_shape_module()
    description = _generator().web_session_admission()

    assert re.findall(r"\bimport\b.*|\bfrom\s+['\"].*", source) == [
        "import { BandedContractDescription } from '../utils/contract-shape.js' */"
    ]
    assert re.search(r"\bany\b", source) is None
    assert "export const SESSION_ADMISSION_DESCRIPTION = Object.freeze({" in source
    for band in description["bands"]:
        versions = ", ".join(map(str, band["versions"]))
        assert (
            f'Object.freeze({{minVersion: {band["minVersion"]}, versions: Object.freeze([{versions}]), '
            f'root: "{band["root"]}"}})'
        ) in source


def test_the_web_contracts_stop_at_the_current_session_version() -> None:
    """D-15-81: Python reads Session 49 before the flip; the Web contracts list
    no version newer than CURRENT_SESSION_VERSION until then."""
    from gbdraw import session_io
    from gbdraw.session_contract import describe_session_admission

    newer = {version for version in session_io.SUPPORTED_SESSION_VERSIONS if version > session_io.CURRENT_SESSION_VERSION}
    python_bands = describe_session_admission()["bands"]
    web_bands = _generator().web_session_admission()["bands"]
    constants = _generator().render_web_contracts_module()

    assert web_bands == [band for band in python_bands if not set(band["versions"]) <= newer]
    for version in newer:
        assert f"{version}" not in re.search(r"SUPPORTED_SESSION_VERSIONS = [^;]*;", constants).group(0)
        assert f"{version}:" not in re.search(r"FEATURE_CATALOG_SCHEMA_BY_SESSION_VERSION = [^;]*;", constants).group(0)


def test_the_session_typedefs_declare_each_band_and_their_union() -> None:
    source = _generator().render_session_typedefs_module()
    roots = [band["root"] for band in _generator().web_session_admission()["bands"]]

    assert f"/** @typedef {{{' | '.join(roots)}}} SessionChecked */" in source
    assert " *   otherModeResult?: never,\n" in source  # a literal with no value admits nothing
    assert " * @property {CheckedSessionResult[]} results\n" in source


def test_a_literal_without_values_renders_as_never() -> None:
    from gbdraw.contract_shapes import literal

    assert _generator()._jsdoc_type(literal()) == "never"


def test_a_changed_session_threshold_makes_the_session_shape_module_stale(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import gbdraw.session_io as session_io

    generator = _generator()
    current = generator.render_session_shape_module()

    monkeypatch.setattr(session_io, "CURRENT_ARTIFACT_SESSION_MIN_VERSION", 40)
    monkeypatch.setattr(session_io, "CURRENT_AUTHORITY_SESSION_MIN_VERSION", 41)

    assert generator.render_session_shape_module() != current
