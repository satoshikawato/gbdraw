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


def _pinned(slot_id: str, renderer: str, radius_px: float, width_px: float = 20.0, **fields) -> CircularTrackSlot:
    return CircularTrackSlot(
        id=slot_id,
        renderer=renderer,
        radius=ScalarSpec(radius_px, "px"),
        width=ScalarSpec(width_px, "px"),
        **fields,
    )


@pytest.mark.parametrize(
    ("slots", "definition_px", "explicit_radius", "message", "context"),
    [
        # Two pinned rows overlap: the later row is named.
        (
            (
                _pinned("gc_content", "dinucleotide_content", 50.0),
                _pinned("PRIVATE_SLOT", "dinucleotide_skew", 55.0),
            ),
            None,
            False,
            "Pinned circular track slot 'PRIVATE_SLOT' overlaps reserved circular slot 'gc_content'.",
            {"reason": "CANNOT_FIT", "slotIndex": 1},
        ),
        # A pinned row overlaps the center definition text ...
        (
            (_pinned("PRIVATE_SLOT", "dinucleotide_skew", 35.0),),
            40.0,
            False,
            "Pinned circular track slot 'PRIVATE_SLOT' overlaps reserved circular slot 'definition'.",
            {"reason": "DEFINITION_RESERVED", "slotIndex": 0},
        ),
        # ... or an explicit center_reserved_radius.
        (
            (_pinned("PRIVATE_SLOT", "dinucleotide_skew", 35.0),),
            40.0,
            True,
            "Pinned circular track slot 'PRIVATE_SLOT' overlaps reserved circular slot 'definition'.",
            {"reason": "CENTER_RESERVED", "slotIndex": 0},
        ),
        # A row pinned farther out than an unpinned row listed before it; the
        # order check skips the pinned feature row between them.
        (
            (
                CircularTrackSlot(id="gc_content", renderer="dinucleotide_content"),
                _pinned("features", "features", 50.0),
                _pinned("PRIVATE_SLOT", "dinucleotide_skew", 90.0, 10.0),
            ),
            None,
            False,
            "Circular track slot order cannot be honored with the supplied pinned geometry: "
            "'PRIVATE_SLOT' would overlap or move outside 'gc_content'.",
            {"reason": "CANNOT_FIT", "slotIndex": 2},
        ),
        # An outside row has no room between the axis and the pinned row above it.
        (
            (
                _pinned("gc_content", "dinucleotide_content", 105.0, 4.0, side="outside"),
                CircularTrackSlot(id="PRIVATE_SLOT", renderer="dinucleotide_skew", side="outside"),
            ),
            None,
            False,
            "Circular track slot 'PRIVATE_SLOT' cannot be placed outside without overlap.",
            {"reason": "CANNOT_FIT", "slotIndex": 1, "innerPx": 101, "outerPx": 102},
        ),
    ],
)
def test_radial_layout_conflicts_report_track_row_without_slot_id(
    slots, definition_px, explicit_radius, message, context
):
    # TK-02: every radial-layout conflict is a TRACK_LAYOUT failure, not a render failure.
    with pytest.raises(ValidationError) as caught:
        resolve_circular_radial_layout(
            total_length=1000,
            canvas_config=_small_radial_canvas(),
            slots=list(slots),
            definition_reserved_radius_px=definition_px,
            center_reserved_radius_explicit=explicit_radius,
        )
    assert str(caught.value) == message
    payload = _web(caught.value)
    assert payload["code"] == "TRACK_LAYOUT"
    assert payload["context"] == context
    assert "PRIVATE" not in json.dumps(payload)


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


@pytest.mark.parametrize(
    ("prefix", "reason"),
    [
        ('PRIVATE:c*?"<>|', "FILENAME"),
        ("nested/PRIVATE", "FILENAME"),
        ("CON", "FILENAME"),
        ("nul.txt", "FILENAME"),
        ("PRIVATE b.", "FILENAME"),
        ("", "REQUIRED"),
        ("P" * 201, "FILENAME_LENGTH"),
        # 101 two-byte characters are 202 bytes in UTF-8.
        ("é" * 101, "FILENAME_LENGTH"),
    ],
    ids=["characters", "folder", "device", "device-extension", "trailing-dot", "empty", "201-bytes", "202-utf8-bytes"],
)
def test_output_prefix_failure_names_the_output_prefix_field(prefix: str, reason: str):
    from gbdraw.api.requests import RenderOutputRequest

    with pytest.raises(ValidationError) as caught:
        RenderOutputRequest(output_prefix=prefix)
    payload = serialize_web_error(caught.value, operation="generate", stage="request-validation")
    assert payload["code"] == "INPUT_INVALID"
    assert payload["context"] == {"field": "output_prefix", "reason": reason}
    assert "PRIVATE" not in json.dumps(payload)


def test_output_prefix_length_limit_counts_utf8_bytes():
    from gbdraw.api.requests import RenderOutputRequest

    assert RenderOutputRequest(output_prefix="x" * 200).output_prefix == "x" * 200
    assert RenderOutputRequest(output_prefix="é" * 100).output_prefix == "é" * 100


def test_radial_inner_label_fit_failure_reports_the_fixed_track_row(tmp_path: Path):
    from gbdraw.circular import _get_args, run_circular_from_namespace

    with pytest.raises(ValidationError) as caught:
        run_circular_from_namespace(_get_args([
            "--gbk", str(Path(__file__).parent / "test_inputs" / "HmmtDNA.gbk"),
            "-o", str(tmp_path / "radial"),
            "--labels", "both", "--label_placement", "radial",
            "--circular_track_slot", "features:features@r=330px,w=40px",
            "--circular_track_slot", "gc_content:dinucleotide_content@side=inside,r=290px,w=40px",
        ]))
    assert str(caught.value).startswith("radial inner labels cannot fit the fixed circular geometry")
    payload = _web(caught.value)
    assert payload["code"] == "TRACK_LAYOUT"
    assert payload["context"] == {"reason": "CANNOT_FIT", "slotIndex": 0}


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


def _record(record_id: str, length: int) -> SeqRecord:
    record = SeqRecord(Seq("A" * length), id=record_id)
    record.annotations["molecule_type"] = "DNA"
    return record


@pytest.mark.parametrize(
    ("start", "end", "message", "reason"),
    [
        (1000, 2000, r"Region 1000\.\.2000 starts after the end of record TINY\.1 \(60 bp\)", "RECORD_BOUNDS"),
        (50, 10, r"Region start \(50\) must not exceed the region end \(10\)", "ORDER"),
    ],
)
def test_region_outside_the_record_names_the_record_and_is_a_region_failure(start, end, message, reason):
    # CI-01: the end was clamped before the order check, so a region wholly
    # beyond the record read "Start position (1000) must be less than end
    # position (2000)" and reached the Web as RENDER_FAILED.
    from gbdraw.crop_genbank import check_start_end_coords

    with pytest.raises(ValidationError, match=message) as caught:
        check_start_end_coords(_record("TINY.1", 60), start, end)
    payload = _web(caught.value)
    assert payload["code"] == "REGION_INVALID"
    assert payload["context"] == {"field": "region", "reason": reason}


def test_linear_record_region_beyond_the_record_keeps_its_region_diagnostic():
    from gbdraw.api import InMemoryRecordSource, RecordInput
    from gbdraw.api.record_planning import resolve_record_inputs
    from gbdraw.io.regions import parse_region_spec

    with pytest.raises(ValidationError, match="TINY.1 \\(60 bp\\)") as caught:
        resolve_record_inputs(
            [RecordInput(source=InMemoryRecordSource(_record("TINY.1", 60)), region=parse_region_spec("1000-2000"))],
            gff_candidate_features=None,
            gff_keep_all_features=False,
        )
    assert _web(caught.value)["code"] == "REGION_INVALID"
    # A partial overrun is still clamped with a warning.
    resolved = resolve_record_inputs(
        [RecordInput(source=InMemoryRecordSource(_record("TINY.1", 60)), region=parse_region_spec("30-2000"))],
        gff_candidate_features=None,
        gff_keep_all_features=False,
    )
    assert len(resolved.records[0]) == 31
