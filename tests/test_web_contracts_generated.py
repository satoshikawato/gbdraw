from __future__ import annotations

import base64
import dataclasses
import functools
import importlib
import importlib.util
import json
import re
import subprocess
import sys
from pathlib import Path
from types import ModuleType
from typing import Any, Callable, Literal, Optional, get_type_hints

import pytest

from gbdraw.session_contract import describe_session_admission

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
        generator.SESSION_SHAPE_TARGET, generator.RULE_SYNTHESIS_TARGET, generator.LOCATE_TARGET, generator.DERIVED_RESULT_SHAPE_TARGET,
        generator.CARD_WINDOWS_TARGET,
    }

    assert set(generator.TARGET.parent.rglob("*")) == targets


@pytest.mark.parametrize(
    ("web_name", "module", "name"),
    [
        ("SESSION_VERSION", "gbdraw.session_io", "CURRENT_SESSION_VERSION"),
        ("CANONICAL_SESSION_MIN_VERSION", "gbdraw.session_io", "CANONICAL_SESSION_MIN_VERSION"),
        ("DISPLAY_PLACEMENT_SCHEMA", "gbdraw.session_request_codec", "DISPLAY_PLACEMENT_SCHEMA"),
        ("SIMILARITY_ALIGNMENT_PLAN_SCHEMA", "gbdraw.layout.similarity_alignment", "SIMILARITY_ALIGNMENT_PLAN_SCHEMA"),
        # D5: the catalog schema older Sessions stored sequences in, and the schema an exported Interactive SVG
        # embeds (the Web export writes it, not the stored FEATURE_CATALOG_SCHEMA).
        ("SEQUENCE_FEATURE_CATALOG_SCHEMA", "gbdraw.session_io", "SEQUENCE_FEATURE_CATALOG_SCHEMA"),
        ("INTERACTIVE_EMBED_SCHEMA", "gbdraw.web_support.feature_catalog", "INTERACTIVE_EMBED_SCHEMA"),
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


def test_the_circular_file_card_slots_come_from_the_python_card(monkeypatch: pytest.MonkeyPatch) -> None:
    # OV-2961: the Web reads a Circular File card's fields (`circularSources` row) from Python's
    # `CircularSourceCard`, each slot typed as itself; a new field makes the generated module stale.
    import gbdraw.session_io as session_io

    assert session_io.CIRCULAR_SOURCE_FIELDS == {"uid", *session_io.CIRCULAR_SOURCE_FILE_SLOTS}
    generator = _generator()
    assert (
        "/** @type {readonly ['gb', 'gff', 'fasta']} */\n"
        "export const CIRCULAR_SOURCE_FILE_SLOTS = Object.freeze(['gb', 'gff', 'fasta']);\n"
    ) in generator.render_web_contracts_module()

    monkeypatch.setattr(session_io, "CIRCULAR_SOURCE_FILE_SLOTS", (*session_io.CIRCULAR_SOURCE_FILE_SLOTS, "depth"))

    assert "Object.freeze(['gb', 'gff', 'fasta', 'depth'])" in generator.render_web_contracts_module()


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
    assert (
        "defaults: Object.freeze({\n        collinearInferOrthogroups: true,\n"
        "        collinearityAnchorMode: \"rbh\",\n        collinearityParams: Object.freeze({"
    ) in source
    assert describe_canonical_request()["types"]["CanonicalProteinComparisonSettings"]["defaults"] == {
        "collinearInferOrthogroups": True,
        # A stored anchor mode that was removed (one_to_one) reads as rbh (OV-2303).
        "collinearityAnchorMode": "rbh",
        # Unset collinearity parameters read as the defaults the encoder writes (OV-3382).
        "collinearityParams": {
            "kind": "lossless",
            "parameters": {
                "minAnchors": 1, "maxUnitGap": 0, "maxDiagonalDrift": 0, "maxConflicts": 1, "mergeOrientation": "strand",
            },
        },
        # A required field that holds null takes null (OV-777).
        **{key: None for key in (
            "losatpBin", "losatpThreads", "ncbiBlastpBin", "orthogroupMemberMaxHits", "proteinBlastpCandidateLimit",
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


def test_the_fresh_protein_comparison_defaults_come_from_the_python_options(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import gbdraw.session_request_codec as codec
    from gbdraw.api.options import LinearDiagramOptions, LosatSearchOptions

    generator = _generator()
    assert codec.describe_fresh_protein_comparison() == {
        "collinearInferOrthogroups": LinearDiagramOptions.collinear_infer_orthogroups,
        "hitLimits": {
            "candidateLimit": LosatSearchOptions.losatp_max_target_seqs,
            "orthogroupMemberMaxHits": LosatSearchOptions.losatp_member_max_hits,
        },
    }
    assert "FRESH_PROTEIN_COMPARISON_DEFAULTS = Object.freeze({\n  collinearInferOrthogroups: false," in (
        generator.render_web_contracts_module()
    )

    monkeypatch.setattr(
        codec, "_DEFAULT_LINEAR_OPTIONS", LinearDiagramOptions(collinear_infer_orthogroups=True)
    )
    monkeypatch.setattr(
        codec,
        "_UNSET_LOSAT_SEARCH",
        LosatSearchOptions(
            program="losatp", losatp_mode="none", losatp_max_target_seqs=9, losatp_member_max_hits=3
        ),
    )

    module = generator.render_web_contracts_module()
    assert "FRESH_PROTEIN_COMPARISON_DEFAULTS = Object.freeze({\n  collinearInferOrthogroups: true," in module
    assert "candidateLimit: 9,\n    orthogroupMemberMaxHits: 3" in module


def test_the_absent_circular_output_comes_from_the_python_options(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import gbdraw.session_request_codec as codec
    from gbdraw.api.options import CIRCULAR_DEFAULT_PLOT_TITLE_POSITION, CircularOutputOptions

    generator = _generator()
    assert codec.describe_absent_circular_output() == {
        "legend": CircularOutputOptions().legend,
        "plotTitlePosition": CIRCULAR_DEFAULT_PLOT_TITLE_POSITION,
    }
    assert "ABSENT_CIRCULAR_OUTPUT = Object.freeze({\n  legend: 'right',\n  plotTitlePosition: 'none'\n})" in (
        generator.render_web_contracts_module()
    )

    monkeypatch.setattr(
        codec, "CircularOutputOptions", functools.partial(CircularOutputOptions, legend="left", plot_title_position="top")
    )
    assert "ABSENT_CIRCULAR_OUTPUT = Object.freeze({\n  legend: 'left',\n  plotTitlePosition: 'top'\n})" in (
        generator.render_web_contracts_module()
    )


def test_the_draft_annotation_style_comes_from_the_python_style_defaults(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # The Web's new annotation style is Python's RegionAnnotationStyle with the highlight fill (OV-2813).
    import gbdraw.session_request_codec as codec
    from gbdraw.annotations import RegionAnnotationStyle
    from gbdraw.annotations.models import DEFAULT_HIGHLIGHT_FILL

    generator = _generator()
    assert codec.describe_draft_annotation_style() == {
        **codec._encode_annotation_style(RegionAnnotationStyle()), "fill": DEFAULT_HIGHLIGHT_FILL,
    }
    current = generator.render_web_contracts_module()
    assert "DRAFT_ANNOTATION_STYLE = Object.freeze({\n  stroke: '#404040',\n  strokeWidth: 1.5," in current
    assert "  labelOffset: 4\n})" in current

    monkeypatch.setattr(
        codec, "RegionAnnotationStyle", functools.partial(RegionAnnotationStyle, stroke="#123456", label_offset=6)
    )
    module = generator.render_web_contracts_module()
    assert "DRAFT_ANNOTATION_STYLE = Object.freeze({\n  stroke: '#123456'," in module
    assert "  labelOffset: 6\n})" in module


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
    description = describe_session_admission()

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


def test_the_session_typedefs_declare_each_band_and_their_union() -> None:
    source = _generator().render_session_typedefs_module()
    roots = [band["root"] for band in describe_session_admission()["bands"]]

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


@pytest.mark.browser
def test_the_circular_conservation_flow_has_one_python_owner() -> None:
    """The Web reads the raw-entry flow of a conservation search from the generator (no JS literal)."""
    from gbdraw.comparisons.circular_losat import CIRCULAR_CONSERVATION_FLOW

    generated = (REPO_ROOT / "gbdraw/web/js/contracts/schema.generated.js").read_text(encoding="utf-8")
    assert f"export const CIRCULAR_CONSERVATION_FLOW = '{CIRCULAR_CONSERVATION_FLOW}';" in generated
    web_js = REPO_ROOT / "gbdraw" / "web" / "js"
    holders = sorted(
        str(path.relative_to(web_js))
        for path in web_js.rglob("*.js")
        if path.parent.name != "contracts"
        and (
            f"'{CIRCULAR_CONSERVATION_FLOW}'" in (text := path.read_text(encoding="utf-8"))
            or "const CIRCULAR_CONSERVATION_FLOW" in text
        )
    )
    assert holders == []


@pytest.mark.browser
def test_the_gallery_drawn_drawing_count_has_one_python_owner() -> None:
    """The Web publication reads how many drawn drawings a Gallery Session shows from the generator."""
    from gbdraw.session_contract import GALLERY_SESSION_DRAWN_DRAWINGS

    generated = (REPO_ROOT / "gbdraw/web/js/contracts/schema.generated.js").read_text(encoding="utf-8")
    assert f"export const GALLERY_SESSION_DRAWN_DRAWINGS = {GALLERY_SESSION_DRAWN_DRAWINGS};" in generated
    publication = (REPO_ROOT / "gbdraw/web/js/services/gallery-session-publication.js").read_text(encoding="utf-8")
    assert "drawn.length !== GALLERY_SESSION_DRAWN_DRAWINGS" in publication


@pytest.mark.browser
def test_the_losatp_no_search_spelling_has_one_python_owner() -> None:
    """Whether a request runs a LOSATP search reads Python's spelling of "no search" from the generator."""
    from gbdraw.api.options import LOSATP_MODE_WIRE, LOSATP_NO_SEARCH_MODE

    assert LOSATP_NO_SEARCH_MODE == LOSATP_MODE_WIRE["none"]
    generated = (REPO_ROOT / "gbdraw/web/js/contracts/schema.generated.js").read_text(encoding="utf-8")
    assert f"export const LOSATP_NO_SEARCH_MODE = '{LOSATP_NO_SEARCH_MODE}';" in generated
    reader = (REPO_ROOT / "gbdraw/web/js/services/session-request.js").read_text(encoding="utf-8")
    assert "comparison.mode !== LOSATP_NO_SEARCH_MODE" in reader


@pytest.mark.browser
def test_the_viewer_base_complements_have_one_python_owner() -> None:
    """The exported viewer reverse-complements with Python's table: the asset holds a slot, not a copy."""
    from gbdraw.render.interactive_svg import STANDALONE_BASE_COMPLEMENTS_SLOT
    from gbdraw.web_support.feature_catalog import BASE_COMPLEMENTS, _reverse_complement_sequence

    assert _reverse_complement_sequence("".join(BASE_COMPLEMENTS)) == "".join(reversed(BASE_COMPLEMENTS.values()))
    asset = (REPO_ROOT / "gbdraw/web/js/services/standalone-interactivity-assets.js").read_text(encoding="utf-8")
    assert asset.count(f"var BASE_COMPLEMENTS = {STANDALONE_BASE_COMPLEMENTS_SLOT};") == 1
    assert "A:'T'" not in asset


def test_the_derived_result_shape_module_is_the_description_as_frozen_data() -> None:
    source = _generator().render_derived_result_shape_module()

    # A leaf imports a leaf (R13): only the type the checker declares.
    assert re.findall(r"\bimport\b.*|\bfrom\s+['\"].*", source) == [
        "import { RootedContractDescription } from '../utils/contract-shape.js' */"
    ]
    assert re.search(r"\bany\b", source) is None


def _typed_results() -> dict[str, list[tuple[str, object]]]:
    """Per derived mode, typed result resources: written ones, empty ones and near misses."""

    import copy
    import json

    from gbdraw.analysis.collinearity import CollinearityAnchor, CollinearityBlock, CollinearityResult
    from gbdraw.analysis.protein_colinearity import OrthogroupGraphResult, OrthogroupMember
    from gbdraw.session_request_codec import encode_canonical_typed_resource

    anchor = CollinearityAnchor(
        query_protein_id="q", subject_protein_id="s", query_record_index=0, subject_record_index=1,
        query_order=0, subject_order=0, query_start=10, query_end=20, subject_start=30, subject_end=40,
        identity=95.0, evalue=1e-10, bitscore=100.0, alignment_length=10, query_feature_svg_id="qf",
        subject_feature_svg_id="sf", source="test", query_unit_id="qu", subject_unit_id="su",
        query_unit_kind="cds", subject_unit_kind="cds", query_locus_id=None, subject_locus_id=None,
        query_display_name="Q", subject_display_name="S",
    )
    member = OrthogroupMember(
        orthogroup_id="og1", protein_id="p1", record_index=0, feature_index=0, record_id="r", label="L",
        start=1, end=9, strand=1, feature_svg_id=None, source_protein_id=None,
    )
    graph = OrthogroupGraphResult(orthogroups={"og1": [member]}, member_by_protein_id={"p1": member})
    block = CollinearityBlock(
        block_id="b1", query_record_index=0, subject_record_index=1, orientation="plus", score=1.0, anchors=(anchor,),
    )

    def encoded(kind: str, value: object) -> dict[str, Any]:
        return json.loads(encode_canonical_typed_resource(kind, value))

    def changed(resource: dict[str, Any], change: Callable[[dict[str, Any]], object]) -> dict[str, Any]:
        result = copy.deepcopy(resource)
        change(result)
        return result

    empty_graph = encoded("orthogroupResult", OrthogroupGraphResult(orthogroups={}, member_by_protein_id={}))
    empty_collinear = encoded("result", CollinearityResult(
        blocks=(), orthogroups=OrthogroupGraphResult(orthogroups={}, member_by_protein_id={}),
    ))
    return {
        "orthogroup": [
            ("written", encoded("orthogroupResult", graph)),
            ("empty", empty_graph),
            ("no fields", changed(empty_graph, lambda r: r["value"].update(fields={}))),
            ("an unknown field", changed(empty_graph, lambda r: r["value"]["fields"].update(extra={}))),
            ("another key", changed(empty_graph, lambda r: r.update(extra=1))),
            ("an older schema", changed(empty_graph, lambda r: r.update(schema=2))),
            ("a fractional schema", changed(empty_graph, lambda r: r.update(schema=3.0))),
            ("another kind", changed(empty_graph, lambda r: r.update(kind="result"))),
            ("an array mapping", changed(empty_graph, lambda r: r["value"]["fields"].update(orthogroups=[]))),
            ("a member without its fields", changed(
                encoded("orthogroupResult", graph), lambda r: r["value"]["fields"]["orthogroups"]["og1"][0].update(fields={})
            )),
        ],
        "collinear": [
            ("written", encoded("result", CollinearityResult(blocks=(block,), unblocked_anchors=(anchor,)))),
            ("empty", empty_collinear),
            ("no orthogroups", changed(empty_collinear, lambda r: r["value"]["fields"].update(orthogroups=None))),
            ("an older schema without orthogroups", changed(empty_collinear, lambda r: (
                r.update(schema=2), r["value"]["fields"].update(orthogroups=None)
            ))),
            ("an older schema with graph orthogroups", changed(empty_collinear, lambda r: r.update(schema=2))),
            ("a boolean schema", changed(empty_collinear, lambda r: (
                r.update(schema=True), r["value"]["fields"].update(orthogroups=None)
            ))),
            ("nested no fields", changed(empty_collinear, lambda r: r["value"]["fields"]["orthogroups"].update(fields={}))),
            ("nested legacy type", changed(
                empty_collinear, lambda r: r["value"]["fields"]["orthogroups"].update(type="OrthogroupResult")
            )),
            ("a missing field", changed(empty_collinear, lambda r: r["value"]["fields"].pop("unblockedAnchors"))),
        ],
    }


def test_the_derived_result_shape_accepts_what_python_decodes() -> None:
    # OV-1926: the Web's zero-hit check reads this shape in place of Python's decode, so both accept the same results.
    from gbdraw.analysis.protein_artifacts import _decodes_as_derived_result, describe_derived_analysis_results
    from gbdraw.contract_shapes import ref, shape_violations

    description = describe_derived_analysis_results()
    verdicts = {
        (mode, label): (
            not shape_violations(resource, ref(description["roots"][mode]), description["types"], "result"),
            _decodes_as_derived_result(mode, resource),
        )
        for mode, cases in _typed_results().items()
        for label, resource in cases
    }

    assert {case: shape for case, (shape, _decoded) in verdicts.items()} == {
        case: decoded for case, (_shape, decoded) in verdicts.items()
    }
    assert {case for case, (_shape, decoded) in verdicts.items() if decoded} == {
        ("orthogroup", "written"), ("orthogroup", "empty"), ("collinear", "written"), ("collinear", "empty"),
        ("collinear", "no orthogroups"), ("collinear", "an older schema without orthogroups"),
    }


_RECORD_SELECTOR_TEXTS = (
    # text, the record the renderer draws, N with 0 (tests/web/record-draw-selection.test.mjs reads the same texts).
    ("#1", 1, 1), ("#01", 1, 1), (" # 2 ", 2, 2), ("#0", None, 0), ("#00", None, 0),
    ("#x", None, None), ("#", None, None), ("#1e1", None, None), ("#+1", 1, 1), ("# +02", 2, 2), ("#+ 1", None, None),
    ("#++1", None, None), ("#-1", None, None), ("#+0", None, 0), ("+#1", None, None), ("#0x2", None, None),
    ("#1.0", None, None), ("BGC0000708", None, None), ("#\u00853", 3, 3), ("#\u001c3", 3, 3), ("#﻿3", None, None),
)


def test_the_web_reads_a_record_selector_number_with_the_renderer_s_pattern(monkeypatch: pytest.MonkeyPatch) -> None:
    # OV-2960: the Web's `#N` reader matches the pattern Python's renderer reads a record selector with, spelled
    # for ECMA-262 (Python's `\s` is not JavaScript's); the generator writes the Python owner's pattern.
    from gbdraw.io import record_select

    for text, record, number in _RECORD_SELECTOR_TEXTS:
        assert record_select.record_selector_number(text) == record, text
        assert record_select.record_index_number(text) == number, text
        parsed = record_select.parse_record_selector(text) if record is not None else None
        assert parsed is None or parsed.record_index == record - 1
    generator = _generator()
    monkeypatch.setattr(record_select, "RECORD_INDEX_SELECTOR_PATTERN", r"#([0-9]+)")

    assert "export const RECORD_INDEX_SELECTOR_PATTERN = '#([0-9]+)';\n" in generator.render_web_contracts_module()


# Each BMP character around and inside a `#N` (a space before, between and after; a digit), read by the Web owner.
_RECORD_INDEX_TEXTS_SCRIPT = """
const { recordIndexNumber } = await import('./gbdraw/web/js/services/record-draw-selection.js');
const read = [];
for (let code = 0; code < 0x10000; code += 1) {
  const c = String.fromCharCode(code);
  read.push(recordIndexNumber(`${c}#${c}1${c}`), recordIndexNumber(`#1${c}`));
}
process.stdout.write(JSON.stringify(read));
"""


@pytest.mark.browser
def test_the_web_owner_reads_every_bmp_character_in_a_record_selector_as_python_does() -> None:
    # Review 1 of #1778: `recordIndexNumber` (services/record-draw-selection.js), run by Node over the BMP, gives
    # the number Python's `record_index_number` gives for each text, so the generated pattern holds in JavaScript.
    from gbdraw.io import record_select

    completed = subprocess.run(
        ["node", "--input-type=module", "-e", _RECORD_INDEX_TEXTS_SCRIPT],
        cwd=REPO_ROOT, check=True, capture_output=True, text=True,
    )
    web = json.loads(completed.stdout)
    texts = [text for code in range(0x10000) for text in (f"{chr(code)}#{chr(code)}1{chr(code)}", f"#1{chr(code)}")]
    python = [record_select.record_index_number(text) for text in texts]
    assert len(web) == len(python)
    differing = [(text, web_number, number) for text, web_number, number in zip(texts, web, python) if web_number != number]
    assert differing == []
    # The texts reach both branches: the spaces Python's `\s` holds (ASCII and others) and the ten ASCII digits.
    assert [python[2 * ord(digit) + 1] for digit in "0123456789"] == list(range(10, 20))
    assert sum(number == 1 for number in python[0::2]) > 20


# Record Order tokens and the row the renderer reads (None where it refuses it); the Web test reads the same file.
_RECORD_ROW_VECTORS = json.loads((REPO_ROOT / "tests" / "fixtures" / "record_order_row_vectors.json").read_text())
_RECORD_ROW_TOKENS = tuple((case["token"], case["row"]) for case in _RECORD_ROW_VECTORS["tokens"])
_RECORD_ROW_LISTS = _RECORD_ROW_VECTORS["rowLists"]
_RECORD_ORDER_PLACEMENTS = _RECORD_ROW_VECTORS["placements"]
_RECORD_ORDER_WINDOW_PLACEMENTS = _RECORD_ROW_VECTORS["windowPlacements"]
_RECORD_ORDER_CIRCULAR_COUNTS = _RECORD_ROW_VECTORS["circularCounts"]
# The refusals of `resolve_record_row_positions` for a token or the count of tokens, not another error of the same call.
_PLACEMENT_REFUSAL = r"multi_record_position (entry|selector|must (provide|include))|Invalid record selector|Record index must be"


def test_the_web_reads_a_record_order_row_with_the_renderer_s_pattern(monkeypatch: pytest.MonkeyPatch) -> None:
    # OV-3313: the Web reads the row of a Record Order token (`<selector>@<row>`) with the pattern Python's
    # `parse_record_row_position` reads it with; the generator writes the Python owner's pattern.
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    from gbdraw.exceptions import ValidationError
    from gbdraw.layout import record_placement

    for token, row in _RECORD_ROW_TOKENS:
        if row is None:
            with pytest.raises(ValidationError, match=r"multi_record_position entry .* must (use a positive integer row|be in)"):
                record_placement.parse_record_row_position(token)
        else:
            assert record_placement.parse_record_row_position(token)[1] == row, token
    # Lists: the order the renderer draws their records in (`resolve_record_row_positions`), None where it refuses a row;
    # rows above 2^53 and up to `int()`'s 4300 digits are read exactly (review 1 L1). The Web test reads the same lists.
    for case in _RECORD_ROW_LISTS:
        records = [SeqRecord(Seq("A"), id=f"r{index}") for index in range(len(case["tokens"]))]
        if case["order"] is None:
            with pytest.raises(ValidationError, match=r"multi_record_position entry .* must use a positive integer row"):
                record_placement.resolve_record_row_positions(records, case["tokens"])
        else:
            assert list(record_placement.resolve_record_row_positions(records, case["tokens"])[0]) == case["order"]
    generator = _generator()
    monkeypatch.setattr(record_placement, "RECORD_ROW_PATTERN", r"([0-9]+)")

    assert "export const RECORD_ROW_PATTERN = '([0-9]+)';\n" in generator.render_web_contracts_module()


def test_the_web_reads_a_record_order_token_s_record_with_the_renderer_s_rule(monkeypatch: pytest.MonkeyPatch) -> None:
    # OV-3473: each list of the shared vectors draws its records in the order Python's renderer resolves
    # (`resolve_record_row_positions` with `record_order_selector`), or Python refuses it; the Web test reads the same
    # lists with `recordOrderPlaces`. The Web strips a selector with the space Python's `str.strip` removes, carried by
    # the generator (`PYTHON_SPACE_PATTERN`).
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    from gbdraw.exceptions import ValidationError
    from gbdraw.layout import record_placement

    assert len(_RECORD_ORDER_PLACEMENTS) >= 20
    for case in _RECORD_ORDER_PLACEMENTS:
        records = [SeqRecord(Seq("A"), id=record_id) for record_id in case["recordIds"]]
        mode = case.get("mode", "linear")
        if case["order"] is None:
            with pytest.raises(ValidationError, match=_PLACEMENT_REFUSAL):
                record_placement.resolve_record_row_positions(records, case["tokens"], _compatibility=mode)
        else:
            resolved = record_placement.resolve_record_row_positions(records, case["tokens"], _compatibility=mode)
            assert list(resolved[0]) == case["order"], case
    assert all(chr(code).isspace() == bool(re.fullmatch(record_placement.PYTHON_SPACE_PATTERN, chr(code)))
               for code in range(0x110000))
    generator = _generator()
    monkeypatch.setattr(record_placement, "PYTHON_SPACE_PATTERN", r" ")

    assert "export const PYTHON_SPACE_PATTERN = ' ';\n" in generator.render_web_contracts_module()


def _vector_record_inputs(case: dict) -> tuple:
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    from gbdraw.api import InMemoryRecordSource, RecordInput
    from gbdraw.io.regions import RecordWindow

    inputs = []
    for number, (record_id, windows) in enumerate(zip(case["recordIds"], case["windows"]), start=1):
        record = SeqRecord(Seq("A" * 4000), id=record_id)
        record.annotations["molecule_type"] = "DNA"
        inputs.append(RecordInput(InMemoryRecordSource(record), record_key=f"record-{number}",
                                  windows=tuple(RecordWindow(start, end) for start, end in windows)))
    return tuple(inputs)


def test_the_web_places_record_order_windows_and_counts_circular_records_as_the_renderer_does() -> None:
    # OV-3313 review 3 M2: the renderer draws a record's windows together, in window order, at its token's place
    # (`_unit_positions`), and a Circular canvas whose records each draw one record counts its tokens against the
    # request's records. The Web test loads the same cases (tests/web/stored-request-read.test.mjs).
    from gbdraw.api import CircularDiagramRequest, LinearDiagramRequest, LinearMultiRecordOptions, build_request_diagram
    from gbdraw.api.options import CircularMultiRecordOptions
    from gbdraw.exceptions import ValidationError

    assert len(_RECORD_ORDER_WINDOW_PLACEMENTS) >= 8 and len(_RECORD_ORDER_CIRCULAR_COUNTS) >= 5
    for case in _RECORD_ORDER_WINDOW_PLACEMENTS:
        request = LinearDiagramRequest(records=_vector_record_inputs(case),
                                       layout=LinearMultiRecordOptions(multi_record_positions=tuple(case["tokens"])))
        if case["cells"] is None:
            with pytest.raises(ValidationError, match=r"multi_record_position must provide exactly \d+ entry"):
                build_request_diagram(request)
            continue
        units = build_request_diagram(request).drawing._gbdraw_track_slot_geometry["records"]
        rows = sorted({unit["axisYpx"] for unit in units})
        cells = [
            unit["recordKey"] + ("" if "windowIndex" not in unit else f"~w{unit['windowIndex'] + 1}")
            + f"@{rows.index(unit['axisYpx']) + 1}"
            for unit in sorted(units, key=lambda unit: (unit["axisYpx"], unit["axisXpx"]))
        ]
        assert cells == case["cells"], case
    for case in _RECORD_ORDER_CIRCULAR_COUNTS:
        request = CircularDiagramRequest(records=_vector_record_inputs(case),
                                         layout=CircularMultiRecordOptions(multi_record_positions=tuple(case["tokens"])))
        if case["refused"]:
            with pytest.raises(ValidationError, match=r"multi_record_position (selector|must include)"):
                build_request_diagram(request)
        else:
            build_request_diagram(request)


_PYTHON_STRIP_SCRIPT = """
const { pythonStrip } = await import('./gbdraw/web/js/services/genbank-header.js');
const read = [];
for (let code = 0; code < 0x10000; code += 1) read.push(pythonStrip(`${String.fromCharCode(code)}a${String.fromCharCode(code)}`));
process.stdout.write(JSON.stringify(read));
"""


_TYPED_LIST_SCRIPT = """
const { TYPED_LIST_SEPARATORS } = await import('./gbdraw/web/js/contracts/schema.generated.js');
const { splitTypedList } = await import('./gbdraw/web/js/utils/typed-list.js');
const { parseBlacklistWords } = await import('./gbdraw/web/js/services/file-imports.js');
const { bytesToText } = await import('./gbdraw/web/js/services/byte-utils.js');
const { readFileSync } = await import('node:fs');
const { texts, files } = JSON.parse(readFileSync(0, 'utf8'));
process.stdout.write(JSON.stringify({
  separators: TYPED_LIST_SEPARATORS,
  items: texts.map(splitTypedList),
  fileWords: files.map((file) => parseBlacklistWords(bytesToText(Buffer.from(file, 'base64'))).words)
}));
"""

_BMP = [chr(code) for code in range(0x10000) if not 0xD800 <= code <= 0xDFFF]


@pytest.mark.browser
def test_the_web_python_strip_removes_every_bmp_space_python_removes() -> None:
    # OV-3760: `pythonStrip` (services/genbank-header.js), run by Node over the BMP, strips what `str.strip` strips.
    completed = subprocess.run(
        ["node", "--input-type=module", "-e", _PYTHON_STRIP_SCRIPT],
        cwd=REPO_ROOT, check=True, capture_output=True, text=True,
    )
    web = json.loads(completed.stdout)
    python = [f"{chr(code)}a{chr(code)}".strip() for code in range(0x10000)]
    assert [code for code in range(0x10000) if web[code] != python[code]] == []


# Each BMP character around and inside a token's row (a space before and after the row; a digit), read by the Web owner.
_RECORD_ROW_TEXTS_SCRIPT = """
const { recordOrderRows } = await import('./gbdraw/web/js/services/record-draw-selection.js');
const read = [];
for (let code = 0; code < 0x10000; code += 1) {
  const c = String.fromCharCode(code);
  read.push(...recordOrderRows([`#1@${c}1${c}`]), ...recordOrderRows([`#1@1${c}`]));
}
process.stdout.write(JSON.stringify(read));
"""


@pytest.mark.browser
def test_the_web_owner_reads_every_bmp_character_in_a_record_order_row_as_python_does() -> None:
    # OV-3313: `recordOrderRows` (services/record-draw-selection.js), run by Node over the BMP, gives the row Python's
    # `parse_record_row_position` gives for each token, or none where Python refuses it.
    from gbdraw.exceptions import ValidationError
    from gbdraw.layout.record_placement import parse_record_row_position

    def python_row(token: str) -> int | None:
        try:
            return parse_record_row_position(token)[1]
        except ValidationError:
            return None

    completed = subprocess.run(
        ["node", "--input-type=module", "-e", _RECORD_ROW_TEXTS_SCRIPT],
        cwd=REPO_ROOT, check=True, capture_output=True, text=True,
    )
    web = json.loads(completed.stdout)
    tokens = [token for code in range(0x10000) for token in (f"#1@{chr(code)}1{chr(code)}", f"#1@1{chr(code)}")]
    python = [python_row(token) for token in tokens]
    assert len(web) == len(python)
    differing = [(token, web_row, row) for token, web_row, row in zip(tokens, web, python) if web_row != row]
    assert differing == []
    # The tokens reach both branches: the spaces Python's `\s` holds and the ten ASCII digits.
    assert [python[2 * ord(digit) + 1] for digit in "0123456789"] == list(range(10, 20))
    assert sum(row == 1 for row in python[0::2]) > 20


@pytest.mark.browser
def test_the_web_cuts_and_trims_a_typed_list_as_python_owns_it() -> None:
    """OV-2792: the CLI and the Web cut a typed list (a label blacklist, ring names) at the same separators and
    trim the same characters from each item's ends, so a CLI Session draws the same in the Web."""

    from gbdraw.core.typed_list import TYPED_LIST_SEPARATORS, TYPED_LIST_TRIMMED, split_typed_list

    cut = [f"a{c}b" for c in _BMP]
    trimmed = [f"{c}a{c},{c}b{c}" for c in _BMP]
    texts = cut + trimmed + ["a, b\u060cc\nd\r\ne\rf", ",\u060c\n a ,, \n\u060cb,", "x\r\n\r\ny", ""]
    # A Blacklist file that Windows Notepad saved: a byte order mark, then words on CRLF lines.
    notepad = "\ufeffhypothetical\r\nputative\u060c unknown\r\n".encode("utf-8")
    completed = subprocess.run(
        ["node", "--input-type=module", "-e", _TYPED_LIST_SCRIPT],
        input=json.dumps({"texts": texts, "files": [base64.b64encode(notepad).decode("ascii")]}),
        cwd=REPO_ROOT, check=True, capture_output=True, text=True, encoding="utf-8",
    )
    web = json.loads(completed.stdout)
    assert tuple(web["separators"]) == TYPED_LIST_SEPARATORS
    differing = [(text, items) for text, items in zip(texts, web["items"]) if items != split_typed_list(text)]
    assert differing == []
    # Four characters cut "a?b" in two: CR, LF, the comma and the Arabic comma, and nothing else.
    assert sorted(text[1] for text in cut if len(split_typed_list(text)) == 2) == ["\n", "\r", ",", "\u060c"]
    # The characters taken from an item's ends are the owner's set, and nothing else (the two commas cut instead).
    taken = {text[0] for text in trimmed if split_typed_list(text) == ["a", "b"]}
    assert taken - {",", "،"} == set(TYPED_LIST_TRIMMED)
    # The same file read by the Web's file import and pasted into --label_blacklist gives the same words.
    expected = ["hypothetical", "putative", "unknown"]
    assert web["fileWords"] == [expected]
    assert split_typed_list(notepad.decode("utf-8")) == expected


_TRIM_SCRIPT = """
const { readFileSync } = await import('node:fs');
const codes = JSON.parse(readFileSync(0, 'utf8'));
process.stdout.write(JSON.stringify(codes.filter((code) => String.fromCharCode(code).trim() === '')));
"""


@pytest.mark.browser
def test_the_typed_list_trims_what_python_strip_or_the_web_trim_took_before() -> None:
    """OV-2792: ``TYPED_LIST_TRIMMED`` is exactly the characters Python's ``str.strip()`` took from an item's ends
    together with those JavaScript's ``trim()`` took, so no item keeps an end character it lost before on either
    side (dropping, say, U+2028 fails here), and no other character is trimmed."""

    from gbdraw.core.typed_list import TYPED_LIST_TRIMMED

    completed = subprocess.run(
        ["node", "--input-type=module", "-e", _TRIM_SCRIPT],
        input=json.dumps([ord(c) for c in _BMP]),
        cwd=REPO_ROOT, check=True, capture_output=True, text=True,
    )
    web_trimmed = {chr(code) for code in json.loads(completed.stdout)}
    python_stripped = {c for c in _BMP if c.isspace()}
    # Both readers take spaces and line breaks; each takes some the other keeps (U+FEFF; U+0085, U+001C-U+001F).
    assert {" ", "\t", "\n", "\u2028", "\u3000"} <= web_trimmed & python_stripped
    assert "\ufeff" in web_trimmed - python_stripped and "\x85" in python_stripped - web_trimmed
    assert len(TYPED_LIST_TRIMMED) == len(set(TYPED_LIST_TRIMMED))
    assert set(TYPED_LIST_TRIMMED) == python_stripped | web_trimmed
