from __future__ import annotations

import re
import xml.etree.ElementTree as ET

import pandas as pd
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from gbdraw.api import (
    GenBankInputSource,
    InMemoryRecordSource,
    LinearComparison,
    LinearDiagramRequest,
    LinearMultiRecordOptions,
    RecordCardinality,
    RecordInput,
    RecordPresentation,
    build_request_diagram,
)
from gbdraw.api.config import apply_config_overrides
from gbdraw.api.diagram import assemble_linear_diagram_from_records
from gbdraw.config.models import GbdrawConfig
from gbdraw.config.toml import load_config_toml
from gbdraw.exceptions import ValidationError
from gbdraw.io.comparisons import COMPARISON_COLUMNS
from gbdraw.layout.linear_multi_record import (
    LinearRecordMeasurement,
    RecordKey,
    record_pairs_between_adjacent_rows,
    resolve_record_row_positions,
    solve_linear_layout,
)
from gbdraw.canvas.linear import LinearCanvasConfigurator
from gbdraw.config.models import LinearRenderProfile
from gbdraw.diagrams.linear.assemble import (
    _RECORD_LOCAL_LINE_KINDS,
    _ROW_ELIGIBLE_LINE_KINDS,
    _split_definition_line_kinds,
)
from gbdraw.layout.linear import CollisionBand
from gbdraw.render.groups.linear.definition import DefinitionGroup
from gbdraw.layout.record_placement import parse_record_row_position


def _records(*lengths: int) -> list[SeqRecord]:
    records = [SeqRecord(Seq("A" * length), id=f"record_{index}") for index, length in enumerate(lengths, 1)]
    for record in records:
        record.annotations["molecule_type"] = "DNA"
    return records


def _comparison(query: int, subject: int) -> LinearComparison:
    row = ["q", "s", 90.0, 100, 0, 0, 10, 100, 20, 110, 1e-20, 200]
    return LinearComparison(
        query,
        subject,
        pd.DataFrame([row], columns=COMPARISON_COLUMNS),
    )


def _record_group(
    root: ET.Element,
    record_id: str,
    record_index: int,
) -> ET.Element:
    namespace = {"svg": "http://www.w3.org/2000/svg"}
    return next(
        group
        for group in root.findall(".//svg:g", namespace)
        if group.attrib.get("data-gbdraw-record-id") == record_id
        and group.attrib.get("data-gbdraw-record-index") == str(record_index)
    )


def _translate(group: ET.Element) -> tuple[float, float]:
    translations = re.findall(
        r"translate\(\s*([-+0-9.eE]+)[,\s]+([-+0-9.eE]+)\s*\)",
        group.attrib.get("transform", ""),
    )
    assert translations
    return (
        sum(float(x_value) for x_value, _y_value in translations),
        sum(float(y_value) for _x_value, y_value in translations),
    )


def test_shared_scale_and_fixed_gap_for_two_by_two_layout() -> None:
    measurements = tuple(
        LinearRecordMeasurement(index, RecordKey(f"key-{index}"), length)
        for index, length in enumerate((1000, 500, 800, 700))
    )
    plan = solve_linear_layout(
        measurements,
        (0, 0, 1, 1),
        available_width=1024,
        record_gap_px=24,
        align_center=False,
    )

    assert plan.px_per_bp == pytest.approx(1000 / 1500)
    assert plan.placement_for_index(0).sequence_width == pytest.approx(2000 / 3)
    assert plan.placement_for_index(1).x == pytest.approx((2000 / 3) + 24)
    assert plan.placement_for_index(2).sequence_width == pytest.approx(1600 / 3)


def test_non_contiguous_rows_and_token_column_order_are_normalized() -> None:
    records = _records(100, 100, 100)
    ordered, rows = resolve_record_row_positions(
        records,
        ("#2@10", "#1@10", "#3@30"),
    )
    assert ordered == (1, 0, 2)
    assert rows == (0, 0, 1)


def test_shared_record_position_grammar_preserves_mode_compatibility() -> None:
    assert parse_record_row_position(
        "record@segment@2",
        _compatibility="circular",
    ) == ("record@segment", 2)
    with pytest.raises(ValidationError, match="'<selector>@<row>' format"):
        parse_record_row_position("record@segment@2")


def test_linear_record_selector_error_text_is_preserved() -> None:
    with pytest.raises(ValidationError, match=r"#0.*out of range for 1 loaded record"):
        resolve_record_row_positions(_records(100), ("#0@1",))


def test_adjacent_row_pairs_exclude_same_row_records() -> None:
    assert record_pairs_between_adjacent_rows((0, 0, 1, 1)) == (
        (0, 2),
        (0, 3),
        (1, 2),
        (1, 3),
    )


def test_layout_rejects_gap_that_consumes_available_width() -> None:
    measurements = (
        LinearRecordMeasurement(0, RecordKey("a"), 10),
        LinearRecordMeasurement(1, RecordKey("b"), 10),
    )
    with pytest.raises(ValidationError, match="record_count=2"):
        solve_linear_layout(
            measurements,
            (0, 0),
            available_width=24,
            record_gap_px=24,
        )


def test_comparison_anchor_can_overlay_x_disjoint_record_header() -> None:
    measurements = (
        LinearRecordMeasurement(
            0,
            RecordKey("top"),
            100,
            bottom_extent=12,
            comparison_bottom_extent=12,
            collision_bands=(
                CollisionBand("body", 0, 100, 0, 12),
                CollisionBand("comparison", 0, 100, 0, 12),
            ),
        ),
        LinearRecordMeasurement(
            1,
            RecordKey("bottom"),
            100,
            top_extent=40,
            comparison_top_extent=18,
            collision_bands=(
                CollisionBand("body", 0, 100, -18, 0),
                CollisionBand("comparison", 0, 100, -18, 0),
                CollisionBand("definition", -100, -10, -40, -20),
            ),
        ),
    )
    plan = solve_linear_layout(
        measurements,
        (0, 1),
        available_width=100,
        first_axis_y=50,
        comparison_height=20,
    )

    bottom = plan.placement_for_index(1)
    assert bottom.axis_y == pytest.approx(100)
    assert bottom.comparison_top_y == pytest.approx(82)
    assert plan.content_top == pytest.approx(50)
    assert plan.content_bottom == pytest.approx(100)
    assert plan.row_gap_resolutions[0].axis_gap == pytest.approx(50)
    assert plan.row_gap_resolutions[0].current_band.kind == "comparison"


def test_multi_record_solver_reserves_comparison_only_on_active_boundary() -> None:
    measurements = tuple(
        LinearRecordMeasurement(
            index,
            RecordKey(f"record-{index}"),
            100,
            top_extent=10,
            bottom_extent=10,
            comparison_top_extent=5,
            comparison_bottom_extent=5,
            collision_bands=(
                CollisionBand("body", 0, 100, -10, 10),
                CollisionBand("comparison", 0, 100, -5, 5),
            ),
        )
        for index in range(3)
    )

    plan = solve_linear_layout(
        measurements,
        (0, 1, 2),
        available_width=100,
        row_gap_px=8,
        comparison_height=60,
        comparison_endpoint_gap_px=4,
        comparison_record_indices_by_boundary={0: (0, 1)},
    )

    first, second = plan.row_gap_resolutions
    assert first.axis_gap == pytest.approx(78)
    assert first.current_band is not None
    assert first.current_band.kind == "comparison"
    assert second.axis_gap == pytest.approx(28)
    assert second.current_band is not None
    assert second.current_band.kind == "body"
    top = plan.placement_for_index(0)
    bottom = plan.placement_for_index(1)
    assert bottom.comparison_top_y - top.comparison_bottom_y == pytest.approx(60)


def test_single_record_rows_publish_boundary_constraints_and_fit_canvas() -> None:
    canvas = assemble_linear_diagram_from_records(
        _records(1000, 1000, 1000),
        cfg=apply_config_overrides(
            None,
            {
                "labels.linear.scope": "none",
                "canvas.show_gc": False,
                "canvas.show_skew": False,
                "canvas.linear.comparison_height": 200,
            },
        ),
        linear_comparisons=[_comparison(0, 1)],
        legend="none",
    )

    geometry = canvas._gbdraw_track_slot_geometry
    first, second = geometry["axisGapConstraints"]
    assert first["currentKind"] == "comparison"
    assert first["nextKind"] == "comparison"
    assert first["clearGapPx"] == pytest.approx(208)
    assert second["currentKind"] != "comparison"
    assert second["nextKind"] != "comparison"
    assert second["clearGapPx"] != pytest.approx(208)

    records = geometry["records"]
    corridor = (
        records[1]["comparisonExclusionBand"]["absoluteTopPx"]
        - records[0]["comparisonExclusionBand"]["absoluteBottomPx"]
    )
    assert corridor == pytest.approx(208)
    assert all(record["collisionBands"] for record in records)

    viewbox_height = float(str(canvas.attribs["viewBox"]).split()[-1])
    painted_bottom = max(
        record["canvasBand"]["absoluteBottomPx"] for record in records
    )
    assert painted_bottom <= viewbox_height


def test_multi_record_local_header_clears_previous_row_body() -> None:
    records = _records(1000, 1000, 1000)
    for index, record in enumerate(records):
        record.annotations["gbdraw_record_label"] = (
            f"Very long record header {index}"
        )
    canvas = assemble_linear_diagram_from_records(
        records,
        cfg=apply_config_overrides(
            None,
            {
                "labels.linear.scope": "none",
                "canvas.show_gc": False,
                "canvas.show_skew": False,
                "objects.definition.linear.font_size.short": 80,
                "objects.definition.linear.font_size.long": 80,
            },
        ),
        layout=LinearMultiRecordOptions(
            multi_record_positions=("#1@1", "#2@2", "#3@2"),
        ),
        legend="none",
    )

    geometry = canvas._gbdraw_track_slot_geometry
    constraint = geometry["axisGapConstraints"][0]
    assert constraint["currentKind"] == "body"
    assert constraint["nextKind"] == "definition"

    top_record = geometry["records"][0]
    lower_definition = next(
        band
        for band in geometry["records"][1]["collisionBands"]
        if band["kind"] == "definition"
    )
    body_bottom = top_record["recordBodyBand"]["absoluteBottomPx"]
    definition_top = (
        geometry["records"][1]["axisYpx"] + lower_definition["topPx"]
    )
    assert definition_top - body_bottom == pytest.approx(
        constraint["clearGapPx"]
    )


def test_comparison_extent_cannot_escape_reserved_record_extent() -> None:
    with pytest.raises(ValidationError, match="cannot exceed top_extent"):
        LinearRecordMeasurement(
            0,
            RecordKey("record"),
            100,
            top_extent=20,
            comparison_top_extent=21,
        )


def test_api_renders_record_local_widths_and_grid_metadata() -> None:
    records = _records(1000, 500, 800, 700)
    canvas = assemble_linear_diagram_from_records(
        records,
        cfg=apply_config_overrides(
            None,
            {
                "labels.linear.scope": "none",
                "canvas.show_gc": False,
                "canvas.show_skew": False,
            },
        ),
        layout=LinearMultiRecordOptions(
            record_gap_px=24,
            multi_record_positions=("#1@1", "#2@1", "#3@2", "#4@2"),
        ),
        legend="none",
    )
    svg = canvas.tostring()
    assert svg.count("data-record-row=") == 4
    assert 'data-record-row="0"' in svg
    assert 'data-record-column="1"' in svg


def test_bottom_legend_follows_last_resolved_row() -> None:
    records = _records(*(1000 for _index in range(10)))
    records[0].features.append(
        SeqFeature(FeatureLocation(100, 300, strand=1), type="CDS")
    )
    svg = assemble_linear_diagram_from_records(
        records,
        cfg=apply_config_overrides(
            None,
            {
                "labels.linear.scope": "none",
                "canvas.show_gc": False,
                "canvas.show_skew": False,
            },
        ),
        layout=LinearMultiRecordOptions(
            multi_record_positions=tuple(
                f"#{index + 1}@{1 if index < 5 else 2}"
                for index in range(10)
            ),
        ),
        legend="bottom",
        plot_title="Resolved rows",
        plot_title_position="bottom",
    ).tostring()
    root = ET.fromstring(svg)
    namespace = {"svg": "http://www.w3.org/2000/svg"}
    groups = {
        group.attrib["id"]: group
        for group in root.findall(".//svg:g", namespace)
        if "id" in group.attrib
    }

    last_row_axis = _translate(_record_group(root, "record_10", 9))[1]
    legend_top = _translate(groups["legend"])[1]
    assert 0 < legend_top - last_row_axis < 120


@pytest.mark.parametrize("ruler_on_axis", [False, True])
def test_multi_record_above_layout_separates_row_definitions_and_record_labels(
    ruler_on_axis: bool,
) -> None:
    records = _records(1000, 800)
    for record in records:
        record.annotations["gbdraw_record_label"] = "TUMSAT-TG-2018"
    records[0].annotations["gbdraw_record_subtitle"] = "chromosome 1"
    records[1].annotations["gbdraw_record_subtitle"] = "chromosome 2"
    for record in records:
        record.features.append(
            SeqFeature(
                FeatureLocation(100, min(700, len(record.seq)), strand=1),
                type="CDS",
            )
        )

    config_dict = load_config_toml("gbdraw.data", "config.toml")
    config_dict["canvas"]["show_gc"] = False
    config_dict["canvas"]["show_skew"] = False
    config_dict["labels"]["linear"]["scope"] = "none"
    config_dict["canvas"]["linear"]["track_layout"] = "above"
    config_dict["canvas"]["linear"]["keep_definition_left_aligned"] = True
    config_dict["canvas"]["linear"]["ruler_on_axis"] = ruler_on_axis
    definition_cfg = config_dict["objects"]["definition"]["linear"]
    definition_cfg["show_replicon"] = False
    definition_cfg["show_accession"] = False
    definition_cfg["show_length"] = False

    svg = assemble_linear_diagram_from_records(
        records,
        cfg=GbdrawConfig.from_dict(config_dict),
        layout=LinearMultiRecordOptions(
            record_gap_px=24,
            multi_record_positions=("#1@1", "#2@1"),
        ),
        legend="bottom",
    ).tostring()
    root = ET.fromstring(svg)
    namespace = {"svg": "http://www.w3.org/2000/svg"}
    groups = {
        group.attrib["id"]: group
        for group in root.findall(".//svg:g", namespace)
        if "id" in group.attrib
    }

    first_record = _record_group(root, "record_1", 0)
    second_record = _record_group(root, "record_2", 1)
    first_local_definition = groups["record_1_definition_record_1"]
    first_row_definition = groups["record_1_definition_record_1_row"]
    second_local_definition = groups["record_2_definition_record_2"]

    def text_values(group: ET.Element) -> list[str]:
        return [
            "".join(text.itertext())
            for text in group.findall(".//svg:text", namespace)
        ]

    # The label describes the whole row; each replicon name stays above its record.
    assert text_values(first_row_definition) == ["TUMSAT-TG-2018"]
    assert text_values(first_local_definition) == ["chromosome 1"]
    assert text_values(second_local_definition) == ["chromosome 2"]
    assert 0 <= _translate(first_row_definition)[0] < _translate(first_record)[0]

    first_feature_y_values = [
        float(y_value)
        for path in first_record.findall(".//svg:path", namespace)
        for _x_value, y_value in re.findall(
            r"[ML]\s*([-+0-9.eE]+)\s*,?\s*([-+0-9.eE]+)",
            path.attrib.get("d", ""),
        )
    ]
    assert first_feature_y_values
    expected_definition_center_y = _translate(first_record)[1] + 0.5 * (
        min(first_feature_y_values) + max(first_feature_y_values)
    )
    assert _translate(first_row_definition)[1] == pytest.approx(
        expected_definition_center_y
    )

    feature_y_values = [
        float(y_value)
        for path in second_record.findall(".//svg:path", namespace)
        for _x_value, y_value in re.findall(
            r"[ML]\s*([-+0-9.eE]+)\s*,?\s*([-+0-9.eE]+)",
            path.attrib.get("d", ""),
        )
    ]
    assert feature_y_values
    feature_top = _translate(second_record)[1] + min(feature_y_values)
    local_text = second_local_definition.find(".//svg:text", namespace)
    assert local_text is not None
    local_font_size = float(local_text.attrib["font-size"])
    local_bottom = _translate(second_local_definition)[1] + (0.5 * local_font_size)
    assert local_bottom <= feature_top
    assert feature_top - local_bottom >= (
        float(config_dict["canvas"]["linear"]["vertical_padding"]) - 0.5
    )


@pytest.mark.parametrize("keep_definition_left_aligned", [True, False])
def test_multi_record_row_keeps_only_row_wide_definition_lines_on_the_left(
    keep_definition_left_aligned: bool,
) -> None:
    """The row heading carries only what every record of the row repeats."""

    def render(labels: tuple[str, str], subtitles: tuple[str, str]) -> dict[str, list[str]]:
        records = _records(1000, 800)
        for record, label, subtitle in zip(records, labels, subtitles, strict=True):
            record.annotations["gbdraw_record_label"] = label
            record.annotations["gbdraw_record_subtitle"] = subtitle
        config_dict = load_config_toml("gbdraw.data", "config.toml")
        config_dict["canvas"]["show_gc"] = False
        config_dict["canvas"]["show_skew"] = False
        config_dict["labels"]["linear"]["scope"] = "none"
        config_dict["canvas"]["linear"]["keep_definition_left_aligned"] = (
            keep_definition_left_aligned
        )
        definition_cfg = config_dict["objects"]["definition"]["linear"]
        definition_cfg["show_replicon"] = False
        definition_cfg["show_accession"] = False
        definition_cfg["show_length"] = False
        svg = assemble_linear_diagram_from_records(
            records,
            cfg=GbdrawConfig.from_dict(config_dict),
            layout=LinearMultiRecordOptions(
                record_gap_px=24,
                multi_record_positions=("#1@1", "#2@1"),
            ),
            legend="none",
        ).tostring()
        root = ET.fromstring(svg)
        namespace = {"svg": "http://www.w3.org/2000/svg"}
        return {
            group.attrib["id"]: [
                "".join(text.itertext())
                for text in group.findall(".//svg:text", namespace)
            ]
            for group in root.findall(".//svg:g", namespace)
            if "definition" in group.attrib.get("id", "")
        }

    shared = "<i>Vibrio harveyi</i>"
    # Italic markup renders as tspans, so itertext() yields the plain name.
    shared_plain = "Vibrio harveyi"

    # Both lines describe the whole row, so both are drawn once beside it.
    repeated = render((shared, shared), ("SB1", "SB1"))
    assert repeated["record_1_definition_record_1_row"] == [shared_plain, "SB1"]
    assert repeated["record_1_definition_record_1"] == []
    assert repeated["record_2_definition_record_2"] == []
    assert "record_2_definition_record_2_row" not in repeated

    # Nothing describes the row as a whole, so the row gets no heading at all.
    distinct = render((shared, "<i>Vibrio owensii</i>"), ("SB1", "XSBZ03"))
    assert "record_1_definition_record_1_row" not in distinct
    assert distinct["record_1_definition_record_1"] == [shared_plain, "SB1"]
    assert distinct["record_2_definition_record_2"] == ["Vibrio owensii", "XSBZ03"]

    # One organism over several replicons: the replicon name is per record, so it
    # must not be promoted to the row heading just because it leads the row.
    partly = render((shared, shared), ("Chromosome 1", "Plasmid pVh1"))
    assert partly["record_1_definition_record_1_row"] == [shared_plain]
    assert partly["record_1_definition_record_1"] == ["Chromosome 1"]
    assert partly["record_2_definition_record_2"] == ["Plasmid pVh1"]

    # A subtitle belongs under its own label, so it follows the label down.
    orphaned = render((shared, "<i>Vibrio owensii</i>"), ("Complete genome", "Complete genome"))
    assert "record_1_definition_record_1_row" not in orphaned
    assert orphaned["record_1_definition_record_1"] == [shared_plain, "Complete genome"]
    assert orphaned["record_2_definition_record_2"] == [
        "Vibrio owensii",
        "Complete genome",
    ]

    # A records table names a row once on its leading record and leaves the rest
    # blank; an empty value has nothing of its own to say, so it does not
    # contradict the row.
    inherited = render((shared, ""), ("SB1", ""))
    assert inherited["record_1_definition_record_1_row"] == [shared_plain, "SB1"]
    assert inherited["record_1_definition_record_1"] == []
    assert inherited["record_2_definition_record_2"] == []

    # The row part is drawn from the leading record, so a name only a follower
    # carries stays above that follower rather than disappearing.
    follower_only = render(("", "<i>Vibrio owensii</i>"), ("", "XSBZ03"))
    assert "record_1_definition_record_1_row" not in follower_only
    assert follower_only["record_1_definition_record_1"] == []
    assert follower_only["record_2_definition_record_2"] == [
        "Vibrio owensii",
        "XSBZ03",
    ]


def test_row_definition_split_covers_every_definition_line_kind() -> None:
    """No line DefinitionGroup can draw may fall outside the row/local split."""

    record = _records(1000)[0]
    record.annotations["gbdraw_record_label"] = "Label"
    record.annotations["gbdraw_record_subtitle"] = "Subtitle"
    record.features.append(
        SeqFeature(
            FeatureLocation(0, 10, strand=1),
            type="source",
            qualifiers={"chromosome": ["1"]},
        )
    )

    config_dict = load_config_toml("gbdraw.data", "config.toml")
    definition_cfg = config_dict["objects"]["definition"]["linear"]
    definition_cfg["show_replicon"] = True
    definition_cfg["show_accession"] = True
    definition_cfg["show_length"] = True
    cfg = GbdrawConfig.from_dict(config_dict)
    canvas_config = LinearCanvasConfigurator(
        num_of_entries=1,
        longest_genome=len(record.seq),
        profile=LinearRenderProfile(cfg),
        legend="none",
    )

    drawable_kinds = {
        line.kind
        for line in DefinitionGroup(record, canvas_config, cfg=cfg).definition_lines
    }
    assert drawable_kinds == (
        set(_ROW_ELIGIBLE_LINE_KINDS) | set(_RECORD_LOCAL_LINE_KINDS)
    )


def test_row_definition_split_is_decided_per_row() -> None:
    """A row that shares nothing loses its heading without affecting other rows."""

    records = _records(1000, 800, 600)
    labels = ("Lambda selected region", "Lividomycin cluster", "Ribostamycin cluster")
    for record, label in zip(records, labels, strict=True):
        record.annotations["gbdraw_record_label"] = label

    local_kinds, row_kinds = _split_definition_line_kinds(
        list(records),
        rows_by_record=(0, 1, 1),
        row_leading_indices={0, 2},
    )

    # Row 0 holds one record, so its label describes the whole row.
    assert row_kinds[0] == frozenset({"name"})
    assert "name" not in local_kinds[0]
    # Row 1 holds two differently labelled records, so neither leads a heading.
    assert row_kinds[1] == frozenset()
    assert row_kinds[2] == frozenset()
    assert "name" in local_kinds[1]
    assert "name" in local_kinds[2]


def test_multi_record_layout_rejects_normalize_length() -> None:
    with pytest.raises(ValidationError, match="normalize_length"):
        assemble_linear_diagram_from_records(
            _records(100, 100),
            cfg=apply_config_overrides(
                None,
                {"canvas.linear.normalize_length": True},
            ),
            layout=LinearMultiRecordOptions(
                multi_record_positions=("#1@1", "#2@1"),
            ),
            legend="none",
        )


@pytest.mark.linear
def test_final_width_label_band_recomputes_first_row_axis_y() -> None:
    records = _records(3000, 3000)
    records[0].features = [
        SeqFeature(
            FeatureLocation(start, start + 12, strand=1),
            type="CDS",
            qualifiers={"product": [f"Long final-width label number {index:02d}"]},
        )
        for index, start in enumerate(range(80, 2640, 80), start=1)
    ]
    config_dict = load_config_toml("gbdraw.data", "config.toml")
    config_dict["canvas"]["show_gc"] = False
    config_dict["canvas"]["show_skew"] = False
    config_dict["labels"]["linear"]["scope"] = "all"
    config_dict["labels"]["rendering"] = "external_only"
    definition_cfg = config_dict["objects"]["definition"]["linear"]
    definition_cfg["show_replicon"] = False
    definition_cfg["show_accession"] = False
    definition_cfg["show_length"] = False

    drawing = assemble_linear_diagram_from_records(
        records,
        cfg=GbdrawConfig.from_dict(config_dict),
        layout=LinearMultiRecordOptions(
            multi_record_positions=("#1@1", "#2@1"),
        ),
        selected_features_set=["CDS"],
        legend="none",
    )
    first_record = drawing._gbdraw_track_slot_geometry["records"][0]
    canvas_band = first_record["canvasBand"]

    assert canvas_band["topPx"] < -100.0
    assert canvas_band["absoluteTopPx"] >= -1e-6


def test_typed_request_preserves_stable_record_keys_in_svg_metadata() -> None:
    records = _records(100, 100)
    request = LinearDiagramRequest(
        records=tuple(
            RecordInput(
                source=InMemoryRecordSource(record),
                record_key=f"stable-{index}",
                presentation=RecordPresentation(grid_row=1, grid_column=index),
            )
            for index, record in enumerate(records, start=1)
        ),
        layout=LinearMultiRecordOptions(),
    )
    svg = build_request_diagram(request).drawing.tostring()
    assert 'data-record-key="stable-1"' in svg
    assert 'data-record-key="stable-2"' in svg


def test_linear_all_source_cards_inherit_rows_and_contiguous_columns(
    tmp_path,
) -> None:
    sources = []
    for source_index, row in enumerate((1, 2), start=1):
        source = tmp_path / f"source-{source_index}.gbk"
        records = _records(100, 120)
        for record_index, record in enumerate(records, start=1):
            record.id = f"source-{source_index}-record-{record_index}"
            record.name = record.id
        SeqIO.write(records, source, "genbank")
        sources.append(
            RecordInput(
                source=GenBankInputSource(source),
                cardinality=RecordCardinality.ALL,
                record_key=f"card-{source_index}",
                presentation=RecordPresentation(grid_row=row),
            )
        )

    svg = build_request_diagram(
        LinearDiagramRequest(
            records=tuple(sources),
            layout=LinearMultiRecordOptions(),
        )
    ).drawing.tostring()
    root = ET.fromstring(svg)
    namespace = {"svg": "http://www.w3.org/2000/svg"}
    placements = [
        (
            group.attrib["data-record-key"],
            int(group.attrib["data-record-row"]),
            int(group.attrib["data-record-column"]),
        )
        for group in root.findall(".//svg:g", namespace)
        if "data-record-row" in group.attrib
    ]

    assert placements == [
        ("card-1:1", 0, 0),
        ("card-1:2", 0, 1),
        ("card-2:1", 1, 0),
        ("card-2:2", 1, 1),
    ]


def test_linear_all_source_card_rejects_explicit_column() -> None:
    with pytest.raises(ValidationError, match="grid_column.*RecordCardinality.ALL"):
        LinearDiagramRequest(
            records=(
                RecordInput(
                    source=InMemoryRecordSource(_records(100)[0]),
                    cardinality=RecordCardinality.ALL,
                    presentation=RecordPresentation(grid_row=1, grid_column=1),
                ),
            ),
            layout=LinearMultiRecordOptions(),
        )
