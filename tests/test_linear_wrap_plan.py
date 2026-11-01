"""The Linear wrap length: its config leaf, its automatic rule, and what it refuses."""

from __future__ import annotations

from dataclasses import replace
import json
from pathlib import Path
import subprocess
import sys
from xml.etree import ElementTree

from Bio import SeqIO
import pytest

from gbdraw.api import plan_request, render_request
from gbdraw.api.config import apply_config_overrides
from gbdraw.api.diagram import build_linear_diagram_result
from gbdraw.api.options import LinearDiagramOptions, LinearMultiRecordOptions, LosatSearchOptions
from gbdraw.api.requests import (
    GenBankInputSource,
    LinearDiagramRequest,
    RecordInput,
    RenderOutputRequest,
    SimilarityAlignmentReference,
)
from gbdraw.config.models import GbdrawConfig
from gbdraw.config.models.canvas import LinearCanvasConfig
from gbdraw.config.modify import modify_config_dict
from gbdraw.config.toml import load_config_toml
from gbdraw.exceptions import ValidationError
from gbdraw.layout.linear_wrap import (
    auto_wrap_length,
    resolve_wrap_length,
    wrap_blocking_inputs,
)
from gbdraw.mode_profiles import (
    effective_linear_track_leaves,
    linear_wrap_default_overrides,
    resolve_mode_profile_overrides,
)
from gbdraw.session_io import CURRENT_SESSION_VERSION
from gbdraw.web_support.mode_scoped_settings import SLICE_CONTAINERS


WRAP = "canvas.linear.wrap_length"
LAMBDA = Path(__file__).parent / "test_inputs" / "NC_001416.gb"
LAMBDA_LENGTH = 48_502


@pytest.mark.parametrize(
    ("longest", "expected"),
    [
        (1, None),
        (10_000, None),  # no shorter row than the floor, so a single row
        (10_001, 10_000),
        (48_502, 10_000),  # phage lambda
        (50_000, 10_000),
        (50_001, 20_000),
        (100_000, 20_000),
        (100_001, 50_000),
        (155_943, 50_000),  # Nicotiana tabacum plastid
        (250_000, 50_000),
        (280_000, 50_000),  # LvMJNV, past five rows at the cap
        (13_033_779, 50_000),
    ],
)
def test_auto_wrap_length_golden_table(longest: int, expected: int | None) -> None:
    assert auto_wrap_length(longest) == expected


@pytest.mark.parametrize(
    ("requested", "lengths", "expected"),
    [
        ("none", [100_000], None),
        (10_000, [100_000], 10_000),
        (10_000, [10_001, 5], 10_000),
        (10_000, [10_000], None),  # no longer than one row
        (100_000, [48_502, 90_000], None),  # draws as today
        ("auto", [48_502, 3_000], 10_000),  # one length, from the longest
        ("auto", [10_000], None),
    ],
)
def test_resolve_wrap_length(requested, lengths, expected) -> None:
    assert resolve_wrap_length(requested, lengths) == expected


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        (None, "none"),
        ("none", "none"),
        (" None ", "none"),
        ("auto", "auto"),
        ("AUTO", "auto"),
        (50_000, 50_000),
        ("50kb", 50_000),
        ("0.5Mb", 500_000),
    ],
)
def test_config_leaf_canonical_values(value, expected) -> None:
    raw = load_config_toml("gbdraw.data", "config.toml")["canvas"]["linear"]
    parsed = LinearCanvasConfig.from_dict({**raw, "wrap_length": value, "track_layout": "above"})
    assert parsed.wrap_length == expected


@pytest.mark.parametrize("value", [0, -5, True, False, 1.5, "0", "abc", "1e3", [], {}, 2**53, str(2**53)])
def test_config_leaf_refuses_other_values_and_names_its_path(value) -> None:
    raw = load_config_toml("gbdraw.data", "config.toml")["canvas"]["linear"]
    with pytest.raises(ValidationError, match="canvas.linear.wrap_length"):
        LinearCanvasConfig.from_dict({**raw, "wrap_length": value, "track_layout": "above"})


def test_packaged_config_does_not_wrap() -> None:
    assert apply_config_overrides(None, None).canvas.linear.wrap_length == "none"
    assert apply_config_overrides(None, None).canvas.linear.wrap_length == "none"


@pytest.mark.parametrize("value", ["auto", 50_000])
def test_a_config_holds_no_wrap_default_before_the_lengths_are_known(value) -> None:
    resolved = resolve_mode_profile_overrides("linear", {WRAP: value})
    assert "canvas.linear.track_layout" not in resolved
    assert "objects.scale.style" not in resolved
    assert "canvas.linear.ruler_on_axis" not in resolved
    assert resolved["objects.axis.linear.stroke_color"] == "lightgray"


@pytest.mark.parametrize("value", ["auto", 50_000])
def test_wrap_defaults_are_what_a_cut_record_draws_with(value) -> None:
    assert linear_wrap_default_overrides({WRAP: value}) == {
        "canvas.linear.track_layout": "above",
        "objects.scale.style": "ruler",
        "canvas.linear.ruler_on_axis": True,
        "objects.axis.linear.stroke_color": "dimgray",
    }


def test_wrap_defaults_leave_what_the_request_names() -> None:
    defaults = linear_wrap_default_overrides(
        {WRAP: "auto", "canvas.linear.track_layout": "below", "objects.axis.linear.stroke_color": "red"}
    )
    assert "canvas.linear.track_layout" not in defaults
    assert "objects.axis.linear.stroke_color" not in defaults


def test_wrap_defaults_restate_the_scale_colors_the_axis_rule_would_overwrite() -> None:
    defaults = linear_wrap_default_overrides({WRAP: "auto", "objects.scale.stroke_color": "red"})
    assert defaults["objects.scale.stroke_color"] == "red"


def test_the_axis_ruler_is_fixed_in_the_wrap_defaults() -> None:
    defaults = linear_wrap_default_overrides({WRAP: "auto", "canvas.linear.ruler_on_axis": False})
    assert defaults["canvas.linear.ruler_on_axis"] is True


def test_hidden_scale_keeps_the_plain_axis_color_in_the_wrap_defaults() -> None:
    defaults = linear_wrap_default_overrides({WRAP: "auto", "objects.scale.show": False})
    assert "objects.axis.linear.stroke_color" not in defaults


@pytest.mark.parametrize("explicit", [{}, {WRAP: "none"}, {WRAP: None}, {WRAP: "NONE"}])
def test_no_wrap_defaults_without_a_wrap_length(explicit) -> None:
    assert linear_wrap_default_overrides(explicit) is None


def test_effective_track_leaves_read_the_request_then_the_packaged_config() -> None:
    assert effective_linear_track_leaves({}) == ("middle", "bar")
    assert effective_linear_track_leaves({WRAP: 5}) == ("middle", "bar")
    assert effective_linear_track_leaves(
        {"canvas.linear.track_layout": "below", "objects.scale.style": "ruler"}
    ) == ("below", "ruler")


@pytest.mark.parametrize(
    "explicit",
    [
        {"canvas.linear.track_layout": "middle"},
        {"objects.scale.style": "bar"},
        {"canvas.linear.normalize_length": True},
    ],
)
def test_a_config_holding_a_wrap_length_loads_whatever_else_it_holds(explicit) -> None:
    raw = modify_config_dict(load_config_toml("gbdraw.data", "config.toml"), {WRAP: 10_000, **explicit})
    assert GbdrawConfig.from_dict(raw).canvas.linear.wrap_length == 10_000


@pytest.mark.parametrize(
    ("explicit", "reason", "words"),
    [
        ({"canvas.linear.track_layout": "middle"}, "WRAP_TRACK_LAYOUT", "above or below"),
        ({"objects.scale.style": "bar"}, "WRAP_SCALE_STYLE", "ruler"),
        ({"canvas.linear.normalize_length": True}, "WRAP_NORMALIZE_LENGTH", "normalize_length"),
    ],
)
def test_option_conflicts_are_refused_when_a_record_is_cut(tmp_path: Path, explicit, reason, words) -> None:
    with pytest.raises(ValidationError, match=words) as caught:
        plan_request(_request(tmp_path, config_overrides={WRAP: 10_000, **explicit}))
    assert caught.value.diagnostic == {"code": "INPUT_INVALID", "field": "wrap_length", "reason": reason}


@pytest.mark.parametrize(
    "explicit",
    [
        {"canvas.linear.track_layout": "middle"},
        {"objects.scale.style": "bar"},
        {"canvas.linear.normalize_length": True},
    ],
)
def test_option_conflicts_are_not_refused_when_no_record_is_cut(tmp_path: Path, explicit) -> None:
    plan_request(_request(tmp_path, config_overrides={WRAP: 100_000, **explicit}))


def test_a_stored_config_conflict_is_refused_when_a_record_is_cut(tmp_path: Path) -> None:
    """A stored config says what its leaves are: middle there was asked for, not defaulted."""

    raw = modify_config_dict(load_config_toml("gbdraw.data", "config.toml"), {WRAP: 10_000})
    with pytest.raises(ValidationError, match="above or below"):
        plan_request(_request(tmp_path, config=raw))
    plain = modify_config_dict(raw, {WRAP: 100_000})
    assert plan_request(_request(tmp_path, config=plain)).build().metadata.wrap is None


@pytest.mark.parametrize(
    ("explicit", "setting"),
    [
        ({"canvas.linear.track_layout": "middle"}, "track_layout middle"),
        ({"canvas.linear.track_layout": "above", "objects.scale.style": "bar"}, "scale_style bar"),
        (
            {
                "canvas.linear.track_layout": "above",
                "objects.scale.style": "ruler",
                "canvas.linear.normalize_length": True,
            },
            "normalize_length",
        ),
    ],
)
def test_a_refusal_names_the_setting_not_an_option_the_user_may_not_have_typed(tmp_path: Path, explicit, setting) -> None:
    """A Session saved uncut and replayed on a longer record reaches these without any command line."""

    raw = modify_config_dict(load_config_toml("gbdraw.data", "config.toml"), {WRAP: 10_000, **explicit})
    with pytest.raises(ValidationError) as caught:
        plan_request(_request(tmp_path, config=raw))
    assert setting in str(caught.value)
    assert "--" not in str(caught.value)


def test_hidden_scale_allows_any_stored_style(tmp_path: Path) -> None:
    plan_request(
        _request(tmp_path, config_overrides={WRAP: "auto", "objects.scale.show": False, "objects.scale.style": "bar"})
    )


def test_the_override_refuses_what_the_leaf_refuses() -> None:
    raw = load_config_toml("gbdraw.data", "config.toml")
    for value in ("fifty", True, 0, -5):
        with pytest.raises(ValidationError, match=WRAP):
            modify_config_dict(raw, {WRAP: value})


@pytest.mark.parametrize(
    ("value", "canonical"),
    [("50kb", 50_000), ("0.5Mb", 500_000), (50_000, 50_000), ("auto", "auto"), ("none", "none"), (None, "none")],
)
def test_the_override_accepts_what_the_leaf_accepts(value, canonical) -> None:
    options = LinearDiagramOptions(config_overrides={WRAP: value})
    assert apply_config_overrides(options.config, options.config_overrides).canvas.linear.wrap_length == canonical


@pytest.mark.parametrize(
    "name",
    [
        "blast_files",
        "comparison_table_file",
        "linear_comparisons",
        "protein_comparisons",
        "orthogroups",
        "collinearity_blocks",
        "losat_search",
        "similarity_alignment",
    ],
)
def test_wrap_blocking_inputs_names_each_comparison_input(name) -> None:
    assert wrap_blocking_inputs(**{name: object()}) == (name,)


def test_wrap_blocking_inputs_passes_empty_inputs() -> None:
    assert wrap_blocking_inputs(blast_files=(), linear_comparisons=None) == ()


def _request(tmp_path: Path, *, records=1, layout=None, **option_fields) -> LinearDiagramRequest:
    return LinearDiagramRequest(
        records=tuple(
            RecordInput(GenBankInputSource(str(LAMBDA)), record_key=f"r{index}")
            for index in range(records)
        ),
        options=LinearDiagramOptions(**option_fields),
        layout=layout,
        output=RenderOutputRequest(output_prefix="out", output_directory=tmp_path, formats=("svg",)),
    )


def _rendered_svg(tmp_path: Path) -> str:
    return (tmp_path / "out.svg").read_text(encoding="utf-8")


def test_comparisons_are_refused_before_any_search(tmp_path: Path, monkeypatch) -> None:
    def no_search(*_args: object, **_kwargs: object) -> None:
        raise AssertionError("a search ran")

    monkeypatch.setattr("gbdraw.api.request_render._resolve_nucleotide_losat", no_search)
    request = _request(
        tmp_path,
        records=2,
        blast_files=(str(tmp_path / "missing.tsv"),),
        config_overrides={WRAP: "auto"},
    )
    with pytest.raises(ValidationError, match="comparison") as caught:
        render_request(request)
    assert caught.value.diagnostic == {"code": "COMPARISON_INPUT", "field": "wrap_length", "reason": "WRAP"}


def test_a_similarity_analysis_is_refused_before_its_search(tmp_path: Path, monkeypatch) -> None:
    def no_search(*_args: object, **_kwargs: object) -> None:
        raise AssertionError("a search ran")

    monkeypatch.setattr("gbdraw.api.request_render.build_request_plan_diagram", no_search)
    request = replace(
        _request(
            tmp_path,
            records=2,
            losat_search=LosatSearchOptions(program="losatp", losatp_mode="similarity_groups"),
            config_overrides={WRAP: 10_000},
        ),
        similarity_alignment=SimilarityAlignmentReference(feature_id="f1"),
    )
    with pytest.raises(ValidationError, match="comparison"):
        render_request(request)


def test_comparisons_are_not_refused_when_no_record_is_cut(tmp_path: Path) -> None:
    empty = tmp_path / "empty.tsv"
    empty.write_text("", encoding="utf-8")
    plan_request(_request(tmp_path, records=2, blast_files=(str(empty),), config_overrides={WRAP: 100_000}))


def test_records_in_one_grid_row_are_refused_when_a_record_is_cut(tmp_path: Path) -> None:
    request = _request(
        tmp_path,
        records=2,
        layout=LinearMultiRecordOptions(multi_record_positions=("#1@1", "#2@1")),
        config_overrides={WRAP: 10_000},
    )
    with pytest.raises(ValidationError, match="one grid row") as caught:
        render_request(request)
    assert caught.value.diagnostic["reason"] == "WRAP_GRID_ROW"


def test_records_in_one_grid_row_are_not_refused_when_no_record_is_cut(tmp_path: Path) -> None:
    plan_request(
        _request(
            tmp_path,
            records=2,
            layout=LinearMultiRecordOptions(multi_record_positions=("#1@1", "#2@1")),
            config_overrides={WRAP: 100_000},
        )
    )


def test_the_builder_decides_as_the_plan_does() -> None:
    record = SeqIO.read(LAMBDA, "genbank")
    wrapped = build_linear_diagram_result([record], options=LinearDiagramOptions(config_overrides={WRAP: 10_000}))
    assert wrapped.metadata.wrap.row_count == 5
    assert 'data-gbdraw-wrap-row' in wrapped.drawing.tostring()
    longer = build_linear_diagram_result([record], options=LinearDiagramOptions(config_overrides={WRAP: 100_000}))
    plain = build_linear_diagram_result([record], options=LinearDiagramOptions())
    assert longer.drawing.tostring() == plain.drawing.tostring()


def _cli(*args: str, cwd: Path, gbk: tuple[Path, ...] = (LAMBDA,)) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-c", "from gbdraw.cli import main; main()", "linear", "--gbk", *map(str, gbk), *args],
        cwd=cwd,
        capture_output=True,
        text=True,
        env={"PYTHONPATH": str(Path(__file__).parents[1]), "PATH": "/usr/bin:/bin"},
    )


@pytest.mark.parametrize(
    ("args", "words"),
    [
        (("--wrap_length", "auto", "--track_layout", "middle"), "above or below"),
        (("--wrap_length", "auto", "--scale_style", "bar"), "ruler"),
        (("--wrap_length", "auto", "--normalize_length"), "normalize_length"),
    ],
)
def test_cli_refuses_with_a_clear_message(tmp_path: Path, args, words) -> None:
    result = _cli(*args, "-o", "out", cwd=tmp_path)
    assert result.returncode != 0
    assert words in result.stderr + result.stdout
    assert "Traceback" not in result.stderr
    assert not (tmp_path / "out.svg").exists()


def test_cli_wraps_only_when_asked_and_defaults_stay(tmp_path: Path) -> None:
    result = _cli("-o", "out", cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    assert (tmp_path / "out.svg").exists()


def _short_record(tmp_path: Path, length: int = 5_000) -> Path:
    record = SeqIO.read(LAMBDA, "genbank")[:length]
    record.annotations["molecule_type"] = "DNA"
    record.annotations["topology"] = "linear"
    path = tmp_path / "short.gb"
    SeqIO.write(record, path, "genbank")
    return path


def _assert_draws_as_without_wrap_length(
    tmp_path: Path, wrap_args: tuple[str, ...], args: tuple[str, ...], gbk: tuple[Path, ...] = (LAMBDA,)
) -> None:
    (tmp_path / "plain").mkdir()
    (tmp_path / "wrap").mkdir()
    plain = _cli(*args, "-o", "out", cwd=tmp_path / "plain", gbk=gbk)
    assert plain.returncode == 0, plain.stderr
    wrapped = _cli(*wrap_args, *args, "-o", "out", cwd=tmp_path / "wrap", gbk=gbk)
    assert wrapped.returncode == 0, wrapped.stderr
    plain_svg = (tmp_path / "plain" / "out.svg").read_bytes()
    assert (tmp_path / "wrap" / "out.svg").read_bytes() == plain_svg


@pytest.mark.parametrize(
    "args",
    [
        (),
        ("--gc", "--skew", "--show_labels", "all"),
        ("--track_layout", "middle"),
        ("--scale_style", "bar"),
        ("--track_layout", "below", "--scale_style", "ruler", "--ruler_on_axis"),
        ("--normalize_length",),
        ("--hide_scale",),
        ("--linear_track_order", "features,gc_content"),
        ("--linear_track_slot", "features:features@side=overlay", "--linear_track_slot", "gc_content:dinucleotide_content@side=below"),
    ],
    ids=["plain", "tracks", "middle", "bar", "axis_ruler", "normalize", "hidden", "order", "overlay"],
)
def test_nothing_to_cut_draws_exactly_the_plain_command(tmp_path: Path, args) -> None:
    _assert_draws_as_without_wrap_length(tmp_path, ("--wrap_length", "100kb"), args)


ORDER_WITH_AXIS = ("--linear_track_order", "features,gc_content", "--gc", "--linear_track_axis_index")


def test_nothing_to_cut_checks_a_track_order_against_the_plain_layout(tmp_path: Path) -> None:
    """Features at the axis boundary: valid for the plain Middle, a conflict only for rows, which are not drawn."""

    _assert_draws_as_without_wrap_length(tmp_path, ("--wrap_length", "100kb"), (*ORDER_WITH_AXIS, "0"))


def test_nothing_to_cut_refuses_what_the_plain_layout_refuses(tmp_path: Path) -> None:
    args = (*ORDER_WITH_AXIS, "1")
    plain = _cli(*args, "-o", "a", cwd=tmp_path)
    wrapped = _cli("--wrap_length", "100kb", *args, "-o", "b", cwd=tmp_path)
    assert plain.returncode != 0
    assert wrapped.returncode != 0
    assert "side='overlay'" in wrapped.stderr + wrapped.stdout
    assert "Traceback" not in wrapped.stderr
    assert not (tmp_path / "b.svg").exists()


def test_nothing_to_cut_under_the_automatic_floor_draws_the_plain_command(tmp_path: Path) -> None:
    _assert_draws_as_without_wrap_length(tmp_path, ("--wrap_length", "auto"), (), gbk=(_short_record(tmp_path),))


def test_nothing_to_cut_in_one_grid_row_draws_the_plain_command(tmp_path: Path) -> None:
    _assert_draws_as_without_wrap_length(
        tmp_path,
        ("--wrap_length", "100kb"),
        ("--multi_record_position", "#1@1", "--multi_record_position", "#2@1"),
        gbk=(LAMBDA, LAMBDA),
    )


def test_nothing_to_cut_with_a_comparison_draws_the_plain_command(tmp_path: Path) -> None:
    blast = tmp_path / "blast.tsv"
    blast.write_text("NC_001416.1\tNC_001416.1\t100\t1000\t0\t0\t1\t1000\t1\t1000\t0.0\t2000\n", encoding="utf-8")
    _assert_draws_as_without_wrap_length(tmp_path, ("--wrap_length", "100kb"), ("-b", str(blast)), gbk=(LAMBDA, LAMBDA))


def test_nothing_to_cut_keeps_the_exit_status_of_the_plain_command(tmp_path: Path) -> None:
    """The plain command's own refusal stands; the option adds none and removes none."""

    plain = _cli("--track_layout", "middle", "--ruler_on_axis", "-o", "a", cwd=tmp_path)
    wrapped = _cli("--wrap_length", "100kb", "--track_layout", "middle", "--ruler_on_axis", "-o", "b", cwd=tmp_path)
    assert (plain.returncode, wrapped.returncode) == (0, 0), (plain.stderr, wrapped.stderr)
    assert (tmp_path / "a.svg").read_bytes() == (tmp_path / "b.svg").read_bytes()


def test_cli_explicit_above_draws_rows_like_the_default(tmp_path: Path) -> None:
    default = _cli("--wrap_length", "10kb", "-o", "a", cwd=tmp_path)
    assert default.returncode == 0, default.stderr
    explicit = _cli("--wrap_length", "10kb", "--track_layout", "above", "--scale_style", "ruler", "-o", "b", cwd=tmp_path)
    assert explicit.returncode == 0, explicit.stderr
    assert (tmp_path / "a.svg").read_bytes() == (tmp_path / "b.svg").read_bytes()


def test_cli_track_order_follows_the_wrap_default_when_a_record_is_cut(tmp_path: Path) -> None:
    result = _cli("--wrap_length", "10kb", "--linear_track_order", "features,gc_content", "-o", "out", cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    assert "data-gbdraw-wrap-row" in (tmp_path / "out.svg").read_text(encoding="utf-8")


def test_cli_track_order_is_checked_against_the_layout_rows_use(tmp_path: Path) -> None:
    default = _cli("--wrap_length", "10kb", *ORDER_WITH_AXIS, "1", "-o", "a", cwd=tmp_path)
    assert default.returncode == 0, default.stderr
    explicit = _cli("--wrap_length", "10kb", "--track_layout", "above", *ORDER_WITH_AXIS, "1", "-o", "b", cwd=tmp_path)
    assert explicit.returncode == 0, explicit.stderr
    drawn = (tmp_path / "a.svg").read_text(encoding="utf-8")
    assert drawn == (tmp_path / "b.svg").read_text(encoding="utf-8")
    assert drawn.count('data-gbdraw-wrap-row="') >= 5


def test_cli_track_order_refused_for_rows_is_refused_when_a_record_is_cut(tmp_path: Path) -> None:
    result = _cli("--wrap_length", "10kb", *ORDER_WITH_AXIS, "0", "-o", "out", cwd=tmp_path)
    assert result.returncode != 0
    assert "side='above'" in result.stderr + result.stdout
    assert "Traceback" not in result.stderr
    assert not (tmp_path / "out.svg").exists()


def test_a_track_order_decided_by_the_record_lengths_is_stored_as_drawn(tmp_path: Path) -> None:
    saved = _cli("--wrap_length", "10kb", *ORDER_WITH_AXIS, "1", "--save_session", "-o", "out", cwd=tmp_path)
    assert saved.returncode == 0, saved.stderr
    replayed = subprocess.run(
        [sys.executable, "-c", "from gbdraw.cli import main; main()", "linear", "--session", "out.gbdraw-session.json",
         "-o", "again"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        env={"PYTHONPATH": str(Path(__file__).parents[1]), "PATH": "/usr/bin:/bin"},
    )
    assert replayed.returncode == 0, replayed.stderr
    assert (tmp_path / "again.svg").read_bytes() == (tmp_path / "out.svg").read_bytes()


def test_cli_ruler_on_axis_with_rows_does_not_warn(tmp_path: Path) -> None:
    result = _cli("--wrap_length", "10kb", "--ruler_on_axis", "-o", "out", cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    assert "ignored" not in result.stdout + result.stderr


OVERLAY = ("--linear_track_slot", "features:features@side=overlay", "--linear_track_slot", "gc:dinucleotide_content@side=below")


@pytest.mark.parametrize(
    ("args", "gbk"),
    [
        (OVERLAY, (LAMBDA,)),
        (OVERLAY, (LAMBDA, LAMBDA)),
        ((*OVERLAY, "--track_layout", "below"), (LAMBDA,)),
        ((*OVERLAY, "--hide_scale"), (LAMBDA,)),
    ],
    ids=["one_record", "two_records", "named_layout", "hidden_scale"],
)
def test_a_feature_slot_on_the_axis_is_refused_when_a_record_is_cut(tmp_path: Path, args, gbk) -> None:
    result = _cli("--wrap_length", "10kb", *args, "-o", "out", cwd=tmp_path, gbk=gbk)
    assert result.returncode != 0
    assert "side=overlay" in result.stderr + result.stdout
    assert "Traceback" not in result.stderr
    assert not (tmp_path / "out.svg").exists()


def test_a_hidden_scale_stays_allowed_with_a_record_cut(tmp_path: Path) -> None:
    result = _cli("--wrap_length", "10kb", "--hide_scale", "-o", "out", cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    assert "data-gbdraw-wrap-row" in (tmp_path / "out.svg").read_text(encoding="utf-8")


@pytest.mark.parametrize(
    "args",
    [
        (),
        ("--track_layout", "below"),
        ("--ruler_on_axis",),
        ("--linear_track_order", "features,gc_content"),
        ("--linear_track_slot", "features:features@side=above", "--linear_track_slot", "gc:dinucleotide_content@side=below"),
        ("--linear_track_slot", "features:features@side=below"),
        ("--linear_track_slot", "gc:dinucleotide_content@side=below"),
        ("--separate_strands",),
        ("--show_labels", "first"),
    ],
    ids=["default", "below", "ruler_flag", "order", "slots_above", "slots_below", "no_feature_slot", "strands", "labels_first"],
)
@pytest.mark.parametrize("gbk", [(LAMBDA,), (LAMBDA, LAMBDA)], ids=["one_record", "two_records"])
def test_every_row_of_every_cut_figure_with_a_shown_scale_has_ticks_and_no_length_bar(tmp_path: Path, args, gbk) -> None:
    result = _cli("--wrap_length", "10kb", *args, "-o", "out", cwd=tmp_path, gbk=gbk)
    assert result.returncode == 0, result.stderr
    root = ElementTree.parse(tmp_path / "out.svg").getroot()
    axes = [element for element in root.iter() if element.attrib.get("data-gbdraw-wrap-part") == "axis"]
    assert len(axes) == 5 * len(gbk)
    assert all(sum(1 for item in axis.iter() if item.tag.endswith("}text")) > 1 for axis in axes)
    assert "length_bar" not in (tmp_path / "out.svg").read_text(encoding="utf-8")


def _session(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_a_cut_figure_stores_what_it_drew_and_replays_byte_identical(tmp_path: Path) -> None:
    saved = _cli("--wrap_length", "10kb", "--gc", "--save_session", "-o", "out", cwd=tmp_path)
    assert saved.returncode == 0, saved.stderr
    linear = _session(tmp_path / "out.gbdraw-session.json")["renderRequest"]["diagramOptions"]
    assert "wrapDefaults" not in linear
    assert linear["config"]["canvas"]["linear"]["track_layout"] == "above"
    assert linear["config"]["canvas"]["linear"]["wrap_length"] == 10_000
    assert linear["config"]["objects"]["scale"]["style"] == "ruler"
    replay = subprocess.run(
        [sys.executable, "-c", "from gbdraw.cli import main; main()", "linear", "--session", "out.gbdraw-session.json",
         "-o", "replayed", "-f", "svg"],
        cwd=tmp_path, capture_output=True, text=True, env={"PYTHONPATH": str(Path(__file__).parents[1]), "PATH": "/usr/bin:/bin"},
    )
    assert replay.returncode == 0, replay.stderr
    assert (tmp_path / "replayed.svg").read_bytes() == (tmp_path / "out.svg").read_bytes()


def test_a_figure_with_nothing_to_cut_stores_the_plain_leaves_and_replays_byte_identical(tmp_path: Path) -> None:
    saved = _cli("--wrap_length", "100kb", "--save_session", "-o", "out", cwd=tmp_path)
    assert saved.returncode == 0, saved.stderr
    linear = _session(tmp_path / "out.gbdraw-session.json")["renderRequest"]["diagramOptions"]
    assert linear["config"]["canvas"]["linear"]["track_layout"] == "middle"
    assert linear["config"]["canvas"]["linear"]["wrap_length"] == 100_000
    assert "data-gbdraw-wrap-row" not in (tmp_path / "out.svg").read_text(encoding="utf-8")
    replay = subprocess.run(
        [sys.executable, "-c", "from gbdraw.cli import main; main()", "linear", "--session", "out.gbdraw-session.json",
         "-o", "replayed", "-f", "svg"],
        cwd=tmp_path, capture_output=True, text=True, env={"PYTHONPATH": str(Path(__file__).parents[1]), "PATH": "/usr/bin:/bin"},
    )
    assert replay.returncode == 0, replay.stderr
    assert (tmp_path / "replayed.svg").read_bytes() == (tmp_path / "out.svg").read_bytes()


def _diagram_options(value: object) -> dict | None:
    if isinstance(value, dict):
        if isinstance(value.get("diagramOptions"), dict):
            return value["diagramOptions"]
        value = list(value.values())
    for item in value if isinstance(value, list) else ():
        found = _diagram_options(item)
        if found is not None:
            return found
    return None


def test_the_field_the_planner_consumes_is_refused_in_a_stored_request(tmp_path: Path) -> None:
    saved = _cli("--wrap_length", "10kb", "--save_session", "-o", "out", cwd=tmp_path)
    assert saved.returncode == 0, saved.stderr
    session_path = tmp_path / "out.gbdraw-session.json"
    session = _session(session_path)
    options = _diagram_options(session)
    assert options is not None and "wrapDefaults" not in options
    options["wrapDefaults"] = {"trackLayout": "above"}
    session_path.write_text(json.dumps(session), encoding="utf-8")
    replay = subprocess.run(
        [sys.executable, "-c", "from gbdraw.cli import main; main()", "linear", "--session", session_path.name,
         "-o", "replayed", "-f", "svg"],
        cwd=tmp_path, capture_output=True, text=True, env={"PYTHONPATH": str(Path(__file__).parents[1]), "PATH": "/usr/bin:/bin"},
    )
    assert replay.returncode != 0
    assert "wrapDefaults" in replay.stderr + replay.stdout
    assert not (tmp_path / "replayed.svg").exists()


def test_session_form_key_reads_in_the_linear_slice_and_replays() -> None:
    from gbdraw.session_io import _append_linear_gui_args

    assert "linear_wrap_length" in SLICE_CONTAINERS["config.form"]
    assert CURRENT_SESSION_VERSION == 50

    def replayed(value: object) -> list[str]:
        run_args: list[str] = []
        invocation: list[str] = []
        try:
            _append_linear_gui_args(
                run_args, invocation, [], session={}, files={}, ui={}, config={},
                form={"linear_wrap_length": value}, adv={}, temp_dir=Path("."),
            )
        except ValidationError:  # no embedded files in this fixture
            pass
        return run_args

    assert replayed("auto")[-2:] == ["--wrap_length", "auto"]
    assert replayed(50_000)[-2:] == ["--wrap_length", "50000"]
    for absent in (None, "", "none"):
        assert "--wrap_length" not in replayed(absent)
