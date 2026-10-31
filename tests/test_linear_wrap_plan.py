"""The Linear wrap length: its config leaf, its automatic rule, and what it refuses."""

from __future__ import annotations

from pathlib import Path
import subprocess
import sys

import pytest

from gbdraw.api import render_request
from gbdraw.api.config import apply_config_overrides
from gbdraw.api.options import LinearDiagramOptions, LinearMultiRecordOptions
from gbdraw.api.requests import (
    GenBankInputSource,
    LinearDiagramRequest,
    RecordInput,
    RenderOutputRequest,
)
from gbdraw.config.models.canvas import LinearCanvasConfig
from gbdraw.config.modify import modify_config_dict
from gbdraw.config.toml import load_config_toml
from gbdraw.exceptions import ValidationError
from gbdraw.layout.linear_wrap import (
    auto_wrap_length,
    resolve_wrap_length,
    wrap_blocking_inputs,
)
from gbdraw.mode_profiles import effective_linear_track_leaves, resolve_mode_profile_overrides
from gbdraw.session_io import CURRENT_SESSION_VERSION
from gbdraw.web_support.mode_scoped_settings import SLICE_CONTAINERS


WRAP = "canvas.linear.wrap_length"
LAMBDA = Path(__file__).parent / "test_inputs" / "NC_001416.gb"
LAMBDA_LENGTH = 48_502


def _config(**overrides):
    return apply_config_overrides(None, resolve_mode_profile_overrides("linear", overrides))


@pytest.mark.parametrize(
    ("longest", "expected"),
    [
        (1, None),
        (10_000, None),  # W-21: no shorter row than the floor, so a single row
        (10_001, 10_000),
        (48_502, 10_000),  # phage lambda (W-07)
        (50_000, 10_000),
        (50_001, 20_000),
        (100_000, 20_000),
        (100_001, 50_000),
        (155_943, 50_000),  # Nicotiana tabacum plastid (W-07)
        (250_000, 50_000),
        (280_000, 50_000),  # LvMJNV (W-07), past five rows at the cap
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
        (10_000, [10_000], None),  # W-08: no longer than one row
        (100_000, [48_502, 90_000], None),  # W-27: draws as today
        ("auto", [48_502, 3_000], 10_000),  # W-09: one length, from the longest
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
    assert parsed.wrap_active is (expected != "none")


@pytest.mark.parametrize("value", [0, -5, True, False, 1.5, "0", "abc", "1e3", [], {}])
def test_config_leaf_refuses_other_values_and_names_its_path(value) -> None:
    raw = load_config_toml("gbdraw.data", "config.toml")["canvas"]["linear"]
    with pytest.raises(ValidationError, match="canvas.linear.wrap_length"):
        LinearCanvasConfig.from_dict({**raw, "wrap_length": value, "track_layout": "above"})


def test_packaged_config_does_not_wrap() -> None:
    assert _config().canvas.linear.wrap_active is False
    assert _config().canvas.linear.wrap_length == "none"


@pytest.mark.parametrize("value", ["auto", 50_000])
def test_defaults_follow_wrapping(value) -> None:
    resolved = resolve_mode_profile_overrides("linear", {WRAP: value})
    assert resolved["canvas.linear.track_layout"] == "above"
    assert resolved["objects.scale.style"] == "ruler"
    assert resolved["canvas.linear.ruler_on_axis"] is True
    assert resolved["objects.axis.linear.stroke_color"] == "dimgray"


def test_explicit_below_keeps_its_layout_while_wrapping() -> None:
    resolved = resolve_mode_profile_overrides(
        "linear", {WRAP: "auto", "canvas.linear.track_layout": "below"}
    )
    assert resolved["canvas.linear.track_layout"] == "below"


def test_the_axis_ruler_is_fixed_while_wrapping() -> None:
    resolved = resolve_mode_profile_overrides(
        "linear", {WRAP: "auto", "canvas.linear.ruler_on_axis": False}
    )
    assert resolved["canvas.linear.ruler_on_axis"] is True


def test_hidden_scale_keeps_the_plain_axis_color_while_wrapping() -> None:
    resolved = resolve_mode_profile_overrides("linear", {WRAP: "auto", "objects.scale.show": False})
    assert resolved["objects.axis.linear.stroke_color"] == "lightgray"
    assert _config(**{WRAP: "auto", "objects.scale.show": False}).canvas.linear.wrap_active


@pytest.mark.parametrize("explicit", [{}, {WRAP: "none"}])
def test_defaults_do_not_change_without_wrapping(explicit) -> None:
    resolved = resolve_mode_profile_overrides("linear", explicit)
    assert "canvas.linear.track_layout" not in resolved
    assert "objects.scale.style" not in resolved
    assert "canvas.linear.ruler_on_axis" not in resolved
    assert effective_linear_track_leaves(explicit) == ("middle", "bar")


def test_effective_track_leaves_read_one_rule() -> None:
    assert effective_linear_track_leaves({WRAP: 5}) == ("above", "ruler")
    assert effective_linear_track_leaves(
        {WRAP: 5, "canvas.linear.track_layout": "below", "objects.scale.style": "bar"}
    ) == ("below", "bar")
    assert effective_linear_track_leaves({"objects.scale.style": "ruler"}) == ("middle", "ruler")


@pytest.mark.parametrize(
    ("explicit", "reason", "words"),
    [
        ({"canvas.linear.track_layout": "middle"}, "WRAP_TRACK_LAYOUT", "above or below"),
        ({"objects.scale.style": "bar"}, "WRAP_SCALE_STYLE", "ruler"),
        ({"canvas.linear.normalize_length": True}, "WRAP_NORMALIZE_LENGTH", "normalize_length"),
    ],
)
def test_option_conflicts_are_refused_on_every_surface(explicit, reason, words) -> None:
    with pytest.raises(ValidationError, match=words) as caught:
        _config(**{WRAP: "auto", **explicit})
    assert caught.value.diagnostic == {"code": "INPUT_INVALID", "field": "wrap_length", "reason": reason}


def test_hidden_scale_allows_any_stored_style() -> None:
    assert _config(**{WRAP: "auto", "objects.scale.show": False, "objects.scale.style": "bar"})


def test_the_override_takes_canonical_values_only() -> None:
    raw = load_config_toml("gbdraw.data", "config.toml")
    for value in ("50kb", True):
        with pytest.raises(ValidationError, match=WRAP):
            modify_config_dict(raw, {WRAP: value})


def test_conflicts_need_wrapping() -> None:
    assert _config(**{"canvas.linear.normalize_length": True}).canvas.linear.normalize_length is True
    assert _config(**{"canvas.linear.track_layout": "middle", "objects.scale.style": "bar"})


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


def test_comparisons_are_refused_before_any_search(tmp_path: Path) -> None:
    request = _request(
        tmp_path,
        records=2,
        blast_files=(str(tmp_path / "missing.tsv"),),
        config_overrides={WRAP: "auto"},
    )
    with pytest.raises(ValidationError, match="comparison") as caught:
        render_request(request)
    assert caught.value.diagnostic == {"code": "COMPARISON_INPUT", "field": "wrap_length", "reason": "WRAP"}


def test_records_in_one_grid_row_are_refused(tmp_path: Path) -> None:
    request = _request(
        tmp_path,
        records=2,
        layout=LinearMultiRecordOptions(multi_record_positions=("#1@1", "#2@1")),
        config_overrides={WRAP: 10_000},
    )
    with pytest.raises(ValidationError, match="one grid row") as caught:
        render_request(request)
    assert caught.value.diagnostic["reason"] == "WRAP_GRID_ROW"


def test_a_record_shorter_than_the_row_draws_one_row(tmp_path: Path) -> None:
    render_request(_request(tmp_path, config_overrides={WRAP: 100_000}))
    assert "<svg" in _rendered_svg(tmp_path)


def test_a_longer_record_says_wrapping_is_not_drawn_yet(tmp_path: Path) -> None:
    with pytest.raises(ValidationError, match="not drawn yet") as caught:
        render_request(_request(tmp_path, config_overrides={WRAP: 10_000}))
    assert caught.value.diagnostic["reason"] == "WRAP_NOT_DRAWN"
    assert not (tmp_path / "out.svg").exists()


def _cli(*args: str, cwd: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-c", "from gbdraw.cli import main; main()", "linear", "--gbk", str(LAMBDA), *args],
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
        (("--wrap_length", "10kb"), "not drawn yet"),
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


def test_cli_longer_row_draws_one_row_with_the_axis_ruler(tmp_path: Path) -> None:
    result = _cli("--wrap_length", "100kb", "--hide_scale", "-o", "hidden", cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    result = _cli("--wrap_length", "100kb", "-o", "shown", cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    assert "WARNING" not in result.stderr


def test_session_form_key_reads_in_the_linear_slice_and_replays() -> None:
    from gbdraw.session_io import _append_linear_gui_args

    assert "linear_wrap_length" in SLICE_CONTAINERS["config.form"]
    assert CURRENT_SESSION_VERSION == 49

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
