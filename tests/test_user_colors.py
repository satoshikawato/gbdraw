"""User colors and SVG keywords are checked once, where gbdraw reads them (P10c, OV-245).

The check accepts what a renderer reads today (svgwrite's ``paint`` values and
the CSS color forms browsers read) and rejects the rest with a gbdraw
``ValidationError``, so the SVG elements need no svgwrite validation.
"""

from __future__ import annotations

import json
import re
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace

import pytest
from pandas import DataFrame
from svgwrite.data.colors import colornames as SVGWRITE_COLOR_NAMES
from svgwrite.data.typechecker import Full11TypeChecker

from gbdraw import cli
from gbdraw.analysis.conservation import conservation_track_gradient_colors
from gbdraw.analysis.depth_tracks import clone_depth_config
from gbdraw.api.prepared import resolve_feature_inputs
from gbdraw.api.request_render import _prepare_diagram_inputs
from gbdraw.config.modify import validate_config_overrides
from gbdraw.exceptions import ValidationError
from gbdraw.io import colors as color_io
from gbdraw.io.colors import is_user_color, load_default_colors, resolve_color_to_hex
from gbdraw.web_support.rule_matching import evaluate_rules_json
from gbdraw.session import load_session_document, materialize_session, session_to_request

REPO = Path(__file__).resolve().parents[1]
EXAMPLES = REPO / "examples"
GALLERY_SESSIONS = sorted((REPO / "gbdraw" / "web" / "gallery" / "sessions").iterdir())

# Values that render today somewhere: svgwrite ``paint`` values, the CSS forms
# browsers read, and color keywords in any letter case.
ACCEPTED = (
    "red", "Red", "RED", " red ", "rebeccapurple", "none", "None", "currentColor", "currentcolor",
    "inherit", "transparent", "#abc", "#AABBCC", "#abcd", "#11223344",
    "rgb(1,2,3)", "rgb( 10 , 20 , 30 )", "rgb(10%,20%,30%)", "RGB(1,2,3)", "rgba(1,2,3,0.5)",
    "rgb(1 2 3 / 50%)", "hsl(120,50%,50%)", "hsla(120deg 50% 50% / .5)",
    "url(#g)", "url(#g) red", "icc-color(x,1)",
    "", "   ",  # svgwrite's paint type takes an empty value (no paint)
)
# Values no renderer reads.
REJECTED = (
    "notacolor", "#ab", "#abcde", "#ggg", "#junk", "rgb(1,2)", "rgb(1,2,3,4,5)",
    "rgb(a,b,c)", "hsl(1,2)", "rgb(1,2,3)/", "url(", "nan",
)


def _run_cli(*args: str) -> None:
    argv = sys.argv
    sys.argv = ["gbdraw", *args]
    try:
        cli.main()
    finally:
        sys.argv = argv


def _tsv(tmp_path: Path, name: str, text: str) -> str:
    path = tmp_path / name
    path.write_text(text, encoding="utf-8")
    return str(path)


@pytest.mark.parametrize("line", ["gc_content\tnotacolor", "CDS\tnotacolor", "skew_high\t#ggg"])
def test_invalid_default_color_is_a_gbdraw_error(
    line: str, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    # OV-245: origin/dev stopped with svgwrite's TypeError traceback.
    with pytest.raises(SystemExit) as stopped:
        _run_cli(
            "circular", "--gbk", str(EXAMPLES / "MellatMJNV.gb"),
            "-d", _tsv(tmp_path, "colors.tsv", f"{line}\n"),
            "-o", str(tmp_path / "out"), "-f", "svg",
        )
    assert stopped.value.code == 1
    feature_type, color = line.split("\t")
    assert f"ERROR: Invalid color {color!r} for feature type {feature_type!r} in the default colors" in (
        capsys.readouterr().err
    )
    assert not (tmp_path / "out.svg").exists()


@pytest.mark.parametrize("header", ["feature_type\tcolor", "Feature_Type\tColor"])
def test_default_colors_header_row_is_not_an_override(header: str, tmp_path: Path) -> None:
    # The Web writes `feature_type<TAB>color` above its -d rows (Run Info recipe, Generate request).
    for load_comparison in (False, True):
        colors = load_default_colors(
            _tsv(tmp_path, "colors.tsv", f"{header}\nCDS\t#ff0000\n"), load_comparison=load_comparison
        )
        assert "feature_type" not in set(colors["feature_type"].str.lower())
        assert colors.set_index("feature_type").at["CDS", "color"] == "#ff0000"
        resolve_feature_inputs(color_table=None, default_colors=colors, feature_visibility_table=None)


def test_invalid_specific_table_color_stays_a_gbdraw_error(tmp_path: Path) -> None:
    table = _tsv(tmp_path, "table.tsv", "CDS\tproduct\tkinase\tnotacolor\n")
    with pytest.raises(ValidationError, match="Invalid color 'notacolor'"):
        color_io.read_color_table(table)


def test_feature_color_dataframes_are_checked_at_resolve_feature_inputs() -> None:
    # Python API and Web canonical tables reach here as DataFrames; origin/dev
    # wrote fill="notacolor" for a color table without captions.
    defaults = load_default_colors("")
    table = DataFrame(
        [["CDS", "product", "kinase", "notacolor"]],
        columns=["feature_type", "qualifier_key", "value", "color"],
    )
    with pytest.raises(ValidationError, match=r"'CDS' in the color table \(row 1\)") as raised:
        resolve_feature_inputs(color_table=table, default_colors=defaults, feature_visibility_table=None)
    assert raised.value.diagnostic == {"code": "TABLE_INVALID", "field": "color", "reason": "COLOR", "row": 1}

    bad_defaults = defaults.copy()
    bad_defaults.loc[bad_defaults["feature_type"] == "CDS", "color"] = "#12"
    with pytest.raises(ValidationError, match="'CDS' in the default colors"):
        resolve_feature_inputs(color_table=None, default_colors=bad_defaults, feature_visibility_table=None)


@pytest.mark.parametrize("value", ACCEPTED)
def test_accepted_colors_pass_every_entry_point(value: str) -> None:
    assert is_user_color(value)
    validate_config_overrides({"objects.gc_content.stroke_color": value})
    defaults = load_default_colors("")
    defaults.loc[defaults["feature_type"] == "CDS", "color"] = value
    # A color table without captions writes its colors as given; captions need
    # a name or hex color (normalize_specific_color_captions, unchanged).
    table = DataFrame(
        [["CDS", "product", "kinase", value]],
        columns=["feature_type", "qualifier_key", "value", "color"],
    )
    resolve_feature_inputs(color_table=table, default_colors=defaults, feature_visibility_table=None)
    expected = value if value else "#000000"  # an empty depth color keeps the configured one
    assert clone_depth_config(SimpleNamespace(fill_color="#000000"), fill_color=value).fill_color == expected


@pytest.mark.parametrize("color", ["rgb(1,2,3)", "transparent", "hsl(120,50%,50%)", "#junk", "#abcd", "notacolor"])
def test_captioned_color_tables_keep_their_hex_domain(color: str) -> None:
    # Every value outside none / name / #RGB / #RRGGBB carries one diagnostic,
    # also through the Web rule check (evaluate_rules_json).
    table = DataFrame(
        [["CDS", "product", "kinase", color, "Kinase"]],
        columns=["feature_type", "qualifier_key", "value", "color", "caption"],
    )
    expected = {"code": "TABLE_INVALID", "field": "color", "reason": "COLOR", "row": 1}
    if is_user_color(color):  # the others stop earlier, at the feature color check
        with pytest.raises(ValidationError, match=re.escape(f"{color!r} in the color table (row 1)")) as raised:
            resolve_feature_inputs(
                color_table=table, default_colors=load_default_colors(""), feature_visibility_table=None
            )
        assert raised.value.diagnostic == expected
    rule = {"feat": "CDS", "qual": "product", "val": "kinase", "color": color, "cap": "Kinase"}
    with pytest.raises(ValidationError) as web:
        evaluate_rules_json("[]", json.dumps([rule]), "color-captions")
    assert web.value.diagnostic == expected


@pytest.mark.parametrize("color", ["#11223344", "#abcd"])
def test_captioned_color_table_with_alpha_hex_is_a_gbdraw_error(color: str) -> None:
    # OV-260: origin/dev raised a raw ValueError from normalize_hex_color.
    table = DataFrame(
        [["CDS", "product", "kinase", "red", "Kinase"], ["CDS", "product", "ligase", color, "Ligase"]],
        columns=["feature_type", "qualifier_key", "value", "color", "caption"],
    )
    with pytest.raises(ValidationError, match=rf"{color!r} in the color table \(row 2\)") as raised:
        resolve_feature_inputs(color_table=table, default_colors=load_default_colors(""), feature_visibility_table=None)
    assert raised.value.diagnostic == {"code": "TABLE_INVALID", "field": "color", "reason": "COLOR", "row": 2}


def test_conservation_default_colors_take_names() -> None:
    # OV-261: origin/dev raised a raw ValueError for a color name in
    # objects.conservation.min_color/max_color.
    assert conservation_track_gradient_colors(
        None, default_min_color="Red", default_max_color="#8B9CC1"
    ) == ("#ff0000", "#8b9cc1")


@pytest.mark.parametrize("bad", ["none", "transparent", "rgba(1,2,3,0.5)", ""])
@pytest.mark.parametrize("leaf", ["min_color", "max_color"])
def test_conservation_default_color_error_names_its_config_path(leaf: str, bad: str) -> None:
    # F7: a COLOR leaf the config accepts (none, transparent, rgba()) but the
    # ring's hex domain rejects stops with a classified error naming the leaf.
    from gbdraw.web_support.error_adapter import serialize_web_error

    defaults = {"default_min_color": "#d6e2f0", "default_max_color": "#8b9cc1", f"default_{leaf}": bad}
    with pytest.raises(ValidationError) as raised:
        conservation_track_gradient_colors(None, **defaults)
    path = f"objects.conservation.{leaf}"
    assert raised.value.diagnostic == {"code": "INPUT_INVALID", "reason": "COLOR", "configPath": path}
    payload = serialize_web_error(raised.value, operation="generate", stage="render")
    assert payload["code"] == "INPUT_INVALID"
    assert payload["context"]["configPath"] == path


@pytest.mark.parametrize("value", REJECTED)
def test_rejected_colors_fail_every_entry_point(value: str) -> None:
    assert not is_user_color(value)
    with pytest.raises(ValidationError) as raised:
        validate_config_overrides({"objects.gc_content.stroke_color": value})
    assert raised.value.diagnostic == {
        "code": "INPUT_INVALID", "configPath": "objects.gc_content.stroke_color", "reason": "COLOR",
    }
    defaults = load_default_colors("")
    defaults.loc[defaults["feature_type"] == "CDS", "color"] = value
    with pytest.raises(ValidationError):
        resolve_feature_inputs(color_table=None, default_colors=defaults, feature_visibility_table=None)
    with pytest.raises(ValidationError):
        clone_depth_config(SimpleNamespace(fill_color="#000000"), fill_color=value)


def test_the_check_accepts_everything_the_older_checks_accepted() -> None:
    svg = Full11TypeChecker()
    for value in (*ACCEPTED, *REJECTED):
        if svg.is_paint(value):  # svgwrite's debug check let it through
            assert is_user_color(value), value
        if color_io._is_specific_table_color(value):  # the -t file domain
            assert is_user_color(value), value
    for name in SVGWRITE_COLOR_NAMES:
        assert is_user_color(name) and is_user_color(name.upper()), name


@pytest.mark.parametrize(
    ("value", "accepted"),
    [("red", True), ("Red", True), ("#ABC", True), ("#junk", False), ("#12", False)],
)
def test_hex_resolver_rejects_malformed_hex(value: str, accepted: bool) -> None:
    # Track-slot skew colors and comparison default colors resolve through here;
    # origin/dev passed any "#..." text through.
    if accepted:
        assert resolve_color_to_hex(value).startswith("#")
    else:
        with pytest.raises(ValidationError, match="Invalid color"):
            resolve_color_to_hex(value)


@pytest.mark.parametrize(
    ("path", "value"),
    [
        ("objects.scale.font_weight", "heavy"),
        ("objects.legends.font_weight", "1200"),
        ("objects.legends.text_anchor", "left"),
        ("objects.legends.dominant_baseline", "top"),
        ("objects.definition.linear.text_anchor", "left"),
        ("objects.text.font_family", " "),
        ("labels.stroke_color.label_stroke_color", "notacolor"),
        ("objects.definition.linear.fill", "notacolor"),
    ],
)
def test_invalid_svg_keywords_stay_errors(path: str, value: str) -> None:
    # origin/dev raised svgwrite's TypeError at render (or wrote the value
    # silently on the debug=False writers); now a gbdraw error names the path.
    with pytest.raises(ValidationError, match=f"{path!r}: {value!r}") as raised:
        validate_config_overrides({path: value})
    assert raised.value.diagnostic["configPath"] == path


@pytest.mark.parametrize(
    ("path", "value"),
    [
        ("objects.scale.font_weight", "650"),
        ("objects.legends.text_anchor", "middle"),
        ("objects.legends.dominant_baseline", "central"),
        ("objects.definition.linear.line_styles.name.fill", "rgba(0,0,0,0.5)"),
    ],
)
def test_svg_keywords_and_css_colors_are_accepted(path: str, value: str) -> None:
    validate_config_overrides({path: value})


@pytest.mark.parametrize(
    ("path", "value"),
    [
        ("objects.scale.font_weight", "Bold"),
        ("objects.scale.font_weight", "700.0"),
        ("objects.scale.font_weight", " 700"),
        ("objects.scale.font_weight", "7e2"),
        ("objects.scale.font_weight", "0"),
        ("objects.scale.font_weight", "1001"),
        ("objects.legends.text_anchor", "MIDDLE"),
        ("objects.legends.text_anchor", " middle"),
        ("objects.legends.dominant_baseline", "Central"),
    ],
)
def test_svg_keywords_keep_svgwrites_case_sensitive_check(path: str, value: str) -> None:
    # D-20: CairoSVG (PNG/PDF/EPS/PS) compares these keywords case-sensitively
    # and falls back to start/baseline for "MIDDLE", so a case variant stays an
    # error, as it was under svgwrite; a color keeps any letter case.
    with pytest.raises(ValidationError, match=f"{path!r}: {value!r}") as raised:
        validate_config_overrides({path: value})
    assert raised.value.diagnostic["configPath"] == path


def test_cli_color_option_error(tmp_path: Path, capsys: pytest.CaptureFixture[str]) -> None:
    with pytest.raises(SystemExit) as stopped:
        _run_cli(
            "circular", "--gbk", str(EXAMPLES / "MellatMJNV.gb"),
            "--block_stroke_color", "notacolor",
            "-o", str(tmp_path / "out"), "-f", "svg",
        )
    assert stopped.value.code == 1
    assert "ERROR: Invalid value for config override 'objects.features.block_stroke_color': 'notacolor'" in (
        capsys.readouterr().err
    )


def test_every_built_in_palette_passes() -> None:
    palettes = color_io.tomllib.loads(
        (REPO / "gbdraw" / "data" / "color_palettes.toml").read_text(encoding="utf-8")
    )
    for palette, colors in palettes.items():
        if not isinstance(colors, dict):  # the file's title
            continue
        resolve_feature_inputs(
            color_table=None,
            default_colors=load_default_colors("", palette=palette),
            feature_visibility_table=None,
        )


@pytest.mark.parametrize("session_path", GALLERY_SESSIONS, ids=lambda path: path.name)
def test_every_gallery_session_passes_request_preparation(session_path: Path) -> None:
    # Decode and input preparation only: the render-time checks (depth, skew
    # slot and conservation default colors) are not reached by this test.
    document = load_session_document(session_path)
    with tempfile.TemporaryDirectory() as output, materialize_session(document, output_directory=output) as live:
        _prepare_diagram_inputs(session_to_request(live))


@pytest.mark.parametrize(
    "color", ["transparent", "rgba(1,2,3,0.5)", "hsl(120,50%,50%)", "#11223344", "rgb(1 2 3 / 50%)", "Red"]
)
def test_accepted_feature_colors_log_no_error(color: str, caplog: pytest.LogCaptureFixture) -> None:
    # The legend compares rows by color identity; a color the check accepts but
    # the identity cannot resolve to #RRGGBB falls back to its text quietly.
    from gbdraw.legend.table import _legend_fill_identity

    with caplog.at_level("DEBUG"):
        assert _legend_fill_identity(color) == (
            "#ff0000" if color == "Red" else color.lower()
        )
    assert [record.message for record in caplog.records if record.levelname == "ERROR"] == []


def _linear_with_default_colors(tmp_path: Path, row: str, *, blast: bool) -> Path:
    args = ["linear", "--gbk", str(EXAMPLES / "LvMJNV.gb")]
    if blast:
        args += [str(EXAMPLES / "MeenMJNV.gb"), "-b", str(EXAMPLES / "LvMJNV.TrcuMJNV.tblastx.out")]
    _run_cli(*args, "-d", _tsv(tmp_path, "colors.tsv", f"{row}\n"), "-o", str(tmp_path / "out"), "-f", "svg")
    return tmp_path / "out.svg"


@pytest.mark.parametrize(
    ("row", "blast"),
    [
        ("pairwise_match_min\tred", False),
        ("collinear_block_plus\tRed", False),
        ("pairwise_match_min\t#11223344", False),
        ("collinear_block_plus\t#11223344", True),
    ],
)
def test_comparison_gradient_colors_that_are_never_interpolated_still_render(
    row: str, blast: bool, tmp_path: Path
) -> None:
    # P1: origin/dev renders these (the key is unused, or the name is read later);
    # the OV-254 error belongs to the interpolation, not to building the configurator.
    assert _linear_with_default_colors(tmp_path, row, blast=blast).exists()


@pytest.mark.parametrize(
    ("option", "row"), [("-t", "CDS\tproduct\tportal\tDarkGrey\tX"), ("-d", "CDS\tDarkGrey")]
)
def test_color_tables_draw_a_mixed_case_color_name(option: str, row: str, tmp_path: Path) -> None:
    # OV-270: the color check reads a name in any case, and so does drawing.
    _run_cli(
        "linear", "--gbk", str(REPO / "tests" / "test_inputs" / "NC_001416.gb"),
        option, _tsv(tmp_path, "colors.tsv", f"{row}\n"), "-o", str(tmp_path / "out"), "-f", "svg",
    )
    assert 'fill="DarkGrey"' in (tmp_path / "out.svg").read_text(encoding="utf-8")


def test_circular_conservation_ignores_the_linear_comparison_colors(tmp_path: Path) -> None:
    _run_cli(
        "circular", "--gbk", str(EXAMPLES / "LvMJNV.gb"),
        "--conservation_blast", str(EXAMPLES / "LvMJNV.TrcuMJNV.tblastx.out"),
        "-d", _tsv(tmp_path, "colors.tsv", "pairwise_match_min\tred\n"),
        "-o", str(tmp_path / "out"), "-f", "svg",
    )
    assert (tmp_path / "out.svg").exists()


def test_comparison_gradient_color_with_alpha_is_a_table_error_where_it_is_interpolated(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    # OV-254: origin/dev raised a raw ValueError from interpolate_color at render.
    with pytest.raises(SystemExit) as stopped:
        _linear_with_default_colors(tmp_path, "pairwise_match_min\t#11223344", blast=True)
    assert stopped.value.code == 1
    assert "ERROR: Invalid color '#11223344' for feature type 'pairwise_match_min'" in capsys.readouterr().err


def test_orientation_identity_gradient_names_the_collinear_color_key() -> None:
    from gbdraw.render.groups.linear.pairwise_match import PairWiseMatchGroup

    group = PairWiseMatchGroup.__new__(PairWiseMatchGroup)
    group.match_min_color = "#ffffff"
    group.collinearity_orientation_colors = {"plus": "#112233", "minus": "#445566"}
    group.collinearity_orientation_min_colors = {"plus": "#11223344", "minus": "#ffeeee"}
    row = SimpleNamespace(
        collinearity_block_id="b1",
        collinearity_color_mode="orientation_identity",
        collinearity_orientation="plus",
    )
    with pytest.raises(ValidationError, match=re.escape("'#11223344' for feature type 'collinear_block_plus_min'")) as raised:
        group.resolve_match_fill_color(row, 0.5, "#000000")
    assert raised.value.diagnostic == {"code": "TABLE_INVALID", "field": "color", "reason": "COLOR"}
