"""Producer-owned failure meaning: ``diagnostic=`` and bounded Web locators."""

from __future__ import annotations

import json
import math
from pathlib import Path
from types import SimpleNamespace

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gbdraw.analysis.depth_tracks import normalize_depth_tracks
from gbdraw.api.options import CircularDiagramOptions, LinearDiagramOptions
from gbdraw.canvas import CircularCanvasConfigurator
from gbdraw.config.models import CircularRenderProfile, GbdrawConfig
from gbdraw.config.modify import modify_config_dict, validate_config_overrides
from gbdraw.config.toml import load_config_toml
from gbdraw.configurators.legend import _legend_half_stroke_width
from gbdraw.diagrams.circular.radial_layout import resolve_circular_radial_layout
from gbdraw.exceptions import GbdrawError, ParseError, ValidationError
from gbdraw.mode_profiles import (
    COMPARISON_THRESHOLD_DOMAINS,
    ComparisonThresholds,
    mode_profiles_payload,
)
from gbdraw.tracks import CircularTrackSlot, ScalarSpec
from gbdraw.web_support.error_adapter import serialize_web_error


def _web(error: BaseException) -> dict:
    return serialize_web_error(error, operation="generate", stage="render")


def test_diagnostic_is_preferred_over_message_and_survives_wrapping():
    native = ValidationError(
        "PRIVATE_MESSAGE must be a positive integer",
        diagnostic={"code": "INPUT_INVALID", "field": "window", "reason": "POSITIVE_INTEGER_OR_AUTO"},
    )
    wrapper = ValidationError("Canonical request could not be converted: PRIVATE")
    wrapper.__cause__ = native
    for error in (native, wrapper):
        payload = _web(error)
        assert payload["code"] == "INPUT_INVALID"
        assert payload["context"] == {"field": "window", "reason": "POSITIVE_INTEGER_OR_AUTO"}
        assert "PRIVATE" not in json.dumps(payload)


def test_diagnostic_keeps_only_allowlisted_identifiers_and_locators():
    error = ValidationError(
        "PRIVATE",
        diagnostic={
            "code": "TRACK_LAYOUT",
            "reason": "PRIVATE_REASON",
            "field": "PRIVATE_FIELD",
            "slotIndex": 2,
            "innerPx": 181,
            "outerPx": -1,
            "row": True,
            "secret": "PRIVATE_VALUE",
            "configPath": "PRIVATE.path",
        },
    )
    assert _web(error) == {
        "code": "TRACK_LAYOUT",
        "operation": "generate",
        "stage": "render",
        "context": {"slotIndex": 2, "innerPx": 181},
    }


def test_unknown_diagnostic_code_falls_back_to_native_classification():
    error = ValidationError(
        "A request requires at least one RecordInput.",
        diagnostic={"code": "PRIVATE_CODE", "reason": "NONNEGATIVE"},
    )
    assert _web(error)["code"] == "INPUT_REQUIRED"


def test_diagnostic_is_pickle_safe_and_optional():
    import pickle

    error = ValidationError("plain")
    assert error.diagnostic is None
    typed = ValidationError("typed", diagnostic={"code": "INPUT_INVALID"})
    restored = pickle.loads(pickle.dumps(typed))
    assert restored.args == ("typed",)
    assert restored.diagnostic == {"code": "INPUT_INVALID"}


@pytest.mark.parametrize(
    ("values", "field", "reason"),
    [
        ({"identity": 150}, "identity", "PERCENT"),
        ({"identity": -5}, "identity", "PERCENT"),
        ({"alignment_length": 12.5}, "alignment_length", "NONNEGATIVE_INTEGER"),
        ({"alignment_length": -10}, "alignment_length", "NONNEGATIVE_INTEGER"),
        ({"bitscore": -1}, "bitscore", "NONNEGATIVE"),
        ({"evalue": math.nan}, "evalue", "NONNEGATIVE"),
    ],
)
def test_comparison_threshold_domains_produce_typed_reasons(values, field, reason):
    arguments = {"evalue": 0.0, "bitscore": 0.0, "identity": 0.0, "alignment_length": 0}
    arguments.update(values)
    with pytest.raises(ValidationError) as caught:
        ComparisonThresholds(**arguments)
    assert _web(caught.value)["context"] == {"field": field, "reason": reason}


def test_comparison_threshold_domains_are_published_to_the_web_profile():
    domains = mode_profiles_payload()["comparisonDomains"]
    assert set(domains) == {"evalue", "bitscore", "identity", "alignmentLength"}
    assert domains["identity"] == {"minimum": 0, "maximum": 100, "integer": False, "reason": "PERCENT"}
    assert domains["alignmentLength"] == {"minimum": 0, "maximum": None, "integer": True, "reason": "NONNEGATIVE_INTEGER"}
    assert set(COMPARISON_THRESHOLD_DOMAINS) == {"evalue", "bitscore", "identity", "alignment_length"}


def _small_radial_canvas():
    config_dict = modify_config_dict(
        load_config_toml("gbdraw.data", "config.toml"),
        {
            "labels.circular.scope": "none",
            "canvas.show_gc": True,
            "canvas.show_skew": True,
            "canvas.circular.track_type": "tuckin",
        },
    )
    cfg = GbdrawConfig.from_dict(config_dict)
    canvas_config = CircularCanvasConfigurator(
        "test", CircularRenderProfile(cfg), "none", SimpleNamespace(seq="N" * 1000)
    )
    canvas_config.radius = 100.0
    return canvas_config


_NUMERIC_STACK = (
    CircularTrackSlot(id="gc_content", renderer="dinucleotide_content"),
    CircularTrackSlot(id="gc_skew", renderer="dinucleotide_skew"),
    CircularTrackSlot(id="PRIVATE_SLOT", renderer="dinucleotide_skew"),
)


@pytest.mark.parametrize(
    ("slots", "explicit_radius", "reason", "slot_index", "cause"),
    [
        # The numeric stack sits directly on the center definition band.
        (_NUMERIC_STACK, False, "DEFINITION_RESERVED", 2, "the center definition text reserves 70.0px"),
        # The same band given as an explicit center_reserved_radius.
        (_NUMERIC_STACK, True, "CENTER_RESERVED", 2, "center_reserved_radius reserves 70.0px"),
        # A pinned inner slot, not the center band, is the inner limit.
        (
            (
                CircularTrackSlot(id="PRIVATE_SLOT", renderer="dinucleotide_content"),
                CircularTrackSlot(id="gc_skew", renderer="dinucleotide_skew"),
                CircularTrackSlot(id="at_skew", renderer="dinucleotide_skew"),
                CircularTrackSlot(
                    id="inner_spacer",
                    renderer="spacer",
                    radius=ScalarSpec(85.0, "px"),
                    width=ScalarSpec(8.0, "px"),
                ),
            ),
            False,
            "CANNOT_FIT",
            0,
            None,
        ),
    ],
)
def test_track_fit_failure_reports_row_and_band_without_slot_id(
    slots, explicit_radius, reason, slot_index, cause
):
    with pytest.raises(ValidationError) as caught:
        resolve_circular_radial_layout(
            total_length=1000,
            canvas_config=_small_radial_canvas(),
            slots=list(slots),
            definition_reserved_radius_px=70.0 if cause is not None else 30.0,
            center_reserved_radius_explicit=explicit_radius,
        )
    message = str(caught.value)
    # The message names the slot that could not be placed and the limiting cause;
    # an explicit radius is not blamed on the species text.
    assert message.startswith("Circular track slot 'PRIVATE_SLOT' cannot fit inside between ")
    if cause is not None:
        assert cause in message
    assert ("species" in message) is (reason == "DEFINITION_RESERVED")
    payload = _web(caught.value)
    assert payload["code"] == "TRACK_LAYOUT"
    context = payload["context"]
    assert context["reason"] == reason
    assert context["slotIndex"] == slot_index
    assert context["innerPx"] >= 0 and context["outerPx"] >= 0
    if cause is not None:  # a pinned slot can leave no band at all (inner > outer)
        assert context["innerPx"] < context["outerPx"]
    assert "PRIVATE" not in json.dumps(payload)
    # A message that merely looks like a layout failure is not classified (B9:
    # in the render stage it is a render failure that names only its class).
    assert _web(ValidationError(message))["code"] == "RENDER_FAILED"
    assert _web(ValidationError(message))["context"] == {"exceptionType": "ValidationError"}


@pytest.mark.parametrize("surface", ["logical", "record-major"])
@pytest.mark.parametrize(
    ("content", "code", "context"),
    [
        ("chr\t1\t4\nchr\t2\tPRIVATE\n", "DEPTH_INVALID", {"reason": "DEPTH_VALUES"}),
        ("chr\t1\t4\nchr\t2\t-5\n", "DEPTH_INVALID", {"reason": "NONNEGATIVE"}),
        ("chr\t1\n", "TABLE_INVALID", {"reason": "THREE_COLUMNS"}),
        (None, "INPUT_UNREADABLE", {}),
    ],
)
def test_depth_failure_reports_series_locator(tmp_path: Path, surface, content, code, context):
    from gbdraw.api.options import DepthTrackInput

    good = tmp_path / "good.tsv"
    good.write_text("chr\t1\t4\nchr\t2\t5\n", encoding="utf-8")
    bad = tmp_path / "PRIVATE_bad.tsv"
    if content is not None:
        bad.write_text(content, encoding="utf-8")
    record = SeqRecord(Seq("ACGT"), id="chr")
    with pytest.raises(GbdrawError) as caught:
        if surface == "logical":
            normalize_depth_tracks(
                [record],
                depth_tracks=[DepthTrackInput(source=str(good)), DepthTrackInput(source=str(bad))],
            )
        else:
            normalize_depth_tracks([record], depth_track_files=[[str(good), str(bad)]])
    payload = _web(caught.value)
    assert payload["code"] == code
    assert payload["context"] == {**context, "seriesIndex": 1}
    assert "PRIVATE" not in json.dumps(payload)
    if content is not None and code == "DEPTH_INVALID" and context["reason"] == "DEPTH_VALUES":
        assert isinstance(caught.value, ParseError)


@pytest.mark.parametrize(
    ("overrides", "reason", "path"),
    [
        ({"objects.scale.interval": 1.5}, "INTEGER", "objects.scale.interval"),
        ({"labels.font_size.short": 0}, "POSITIVE", "labels.font_size.short"),
        ({"objects.definition.linear.line_styles.name.font_size": -1.0}, "POSITIVE", "objects.definition.linear.line_styles.name.font_size"),
        ({"labels.stroke_width.long": -0.5}, "NONNEGATIVE", "labels.stroke_width.long"),
        ({"objects.depth.tick_font_size": math.inf}, "POSITIVE", "objects.depth.tick_font_size"),
        ({"objects.depth.tick_font_size": 0.0}, "POSITIVE", "objects.depth.tick_font_size"),
        ({"objects.ticks.tick_width": -1.0}, "NONNEGATIVE", "objects.ticks.tick_width"),
    ],
)
def test_config_leaf_failure_reports_canonical_setting(overrides, reason, path):
    with pytest.raises(ValidationError) as caught:
        validate_config_overrides(overrides)
    payload = _web(caught.value)
    assert payload["code"] == "INPUT_INVALID"
    assert payload["context"] == {"reason": reason, "configPath": path}


@pytest.mark.parametrize(
    ("options_type", "overrides", "reason", "path"),
    [
        (LinearDiagramOptions, {"objects.ticks.tick_width": 4, "canvas.circular.radius": 1.2}, "CIRCULAR_SETTING", "canvas.circular.radius"),
        (CircularDiagramOptions, {"objects.blast_match.curve_tension": 0.3}, "LINEAR_SETTING", "objects.blast_match.curve_tension"),
    ],
)
def test_other_mode_config_override_reports_setting_and_owning_mode(options_type, overrides, reason, path):
    # OV-130: the Web names the first other-mode setting and the mode it belongs to.
    with pytest.raises(ValidationError, match="cannot target") as caught:
        options_type(config_overrides=overrides)
    payload = _web(caught.value)
    assert payload["code"] == "MODE_SETTING"
    assert payload["context"] == {"reason": reason, "configPath": path}


def test_style_domains_cover_full_configs_but_not_offsets():
    config = load_config_toml("gbdraw.data", "config.toml")
    config["objects"]["legends"]["font_size"]["short"] = -3
    with pytest.raises(ValidationError) as caught:
        CircularDiagramOptions(config=config)
    assert _web(caught.value)["context"] == {
        "reason": "POSITIVE",
        "configPath": "objects.legends.font_size.short",
    }
    # D-27 keeps the current domain of settings whose negative meaning is unverified.
    LinearDiagramOptions(
        config_overrides={
            "canvas.linear.track_axis_gap": -5.0,
            "labels.linear.rotation": -30.0,
            "labels.spacing.linear": -1.0,
            "canvas.linear.vertical_offset": -10.0,
        }
    )
    CircularDiagramOptions(
        config_overrides={
            "labels.unified_adjustment.outer_labels.x_radius_offset": -0.5,
            "objects.features.block_stroke_width.short": 0.0,
        }
    )


def test_legend_stroke_width_failure_is_typed_not_value_error():
    with pytest.raises(ValidationError) as caught:
        _legend_half_stroke_width({"entry": {"stroke": "black", "width": -1}})
    assert _web(caught.value)["context"] == {"reason": "NONNEGATIVE"}


def test_specific_color_table_failure_reports_row_without_value(tmp_path: Path):
    from gbdraw.io.colors import read_color_table

    table = tmp_path / "PRIVATE_colors.tsv"
    table.write_text(
        "CDS\tproduct\tNADH\t#ff0000\tRed\nCDS\tproduct\tPRIVATE\tnotacolor\tBad\n",
        encoding="utf-8",
    )
    with pytest.raises(ValidationError) as caught:
        read_color_table(str(table))
    payload = _web(caught.value)
    assert payload["code"] == "TABLE_INVALID"
    assert payload["context"] == {"field": "color", "reason": "COLOR", "row": 2}
    assert "PRIVATE" not in json.dumps(payload)


def test_display_start_beyond_the_record_is_a_typed_input_failure():
    from gbdraw.api.record_planning import resolve_record_display
    from gbdraw.api.requests import RecordDisplayOptions

    with pytest.raises(ValidationError) as caught:
        resolve_record_display(
            RecordDisplayOptions(start_coordinate=500),
            source_length=100,
            detected_topology="circular",
            source_base=1,
            source_step=1,
            has_input_region=False,
            has_collection_region=False,
            is_cropped=False,
        )
    payload = _web(caught.value)
    assert payload["code"] == "INPUT_INVALID"
    assert payload["context"] == {"field": "start", "reason": "DISPLAY_START_BOUNDS"}


def test_unparsable_genbank_is_unreadable_input_not_a_render_failure(tmp_path: Path):
    from gbdraw.io.genome import load_gbks

    broken = tmp_path / "PRIVATE_broken.gb"
    broken.write_text("LOCUS       PRIVATE_BROKEN\nFEATURES             Location/Qualifiers\n     CDS  PRIVATE\n", encoding="utf-8")
    with pytest.raises(ParseError) as caught:
        load_gbks([str(broken)])
    payload = _web(caught.value)
    assert payload["code"] == "INPUT_UNREADABLE"
    assert payload["context"] == {}
    assert "PRIVATE" not in json.dumps(payload)
