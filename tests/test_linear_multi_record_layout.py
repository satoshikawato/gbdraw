from __future__ import annotations

import json
import re
import xml.etree.ElementTree as ET

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from gbdraw.api import (
    InMemoryRecordSource,
    LinearComparison,
    LinearDiagramRequest,
    LinearMultiRecordOptions,
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
from gbdraw.layout.linear import CollisionBand
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


def test_row_definition_gutter_preserves_shared_row_biological_geometry() -> None:
    lengths = (1000, 800, 600)
    namespace = {"svg": "http://www.w3.org/2000/svg"}
    cfg = apply_config_overrides(
        None,
        {
            "labels.linear.scope": "none",
            "canvas.show_gc": False,
            "canvas.show_skew": False,
            "canvas.linear.keep_definition_left_aligned": True,
        },
    )

    def render(*, labeled: bool) -> dict[str, object]:
        records = _records(*lengths)
        if labeled:
            records[0].annotations["gbdraw_record_label"] = "Species alpha"
            records[0].annotations["gbdraw_record_subtitle"] = "strain beta"
        root = ET.fromstring(
            assemble_linear_diagram_from_records(
                records,
                cfg=cfg,
                layout=LinearMultiRecordOptions(
                    record_gap_px=24,
                    multi_record_positions=("#1@1", "#2@1", "#3@1"),
                ),
                legend="none",
            ).tostring()
        )
        biological_groups = sorted(
            (
                group
                for group in root.findall(".//svg:g", namespace)
                if "data-record-index" in group.attrib
            ),
            key=lambda group: int(group.attrib["data-record-index"]),
        )
        widths = []
        for group in biological_groups:
            axis = next(
                line
                for line in group.findall("./svg:line", namespace)
                if float(line.attrib["y1"]) == float(line.attrib["y2"]) == 0.0
            )
            widths.append(float(axis.attrib["x2"]) - float(axis.attrib["x1"]))
        record_x = tuple(_translate(group)[0] for group in biological_groups)
        return {
            "root": root,
            "placements": tuple(
                (
                    int(group.attrib["data-record-row"]),
                    int(group.attrib["data-record-column"]),
                )
                for group in biological_groups
            ),
            "sequence_widths": tuple(widths),
            "px_per_bp": tuple(
                width / len(record.seq) for width, record in zip(widths, records)
            ),
            "record_x": record_x,
            "view_box": tuple(
                float(value) for value in root.attrib["viewBox"].split()
            ),
        }

    unlabeled = render(labeled=False)
    labeled = render(labeled=True)

    assert unlabeled["placements"] == labeled["placements"] == (
        (0, 0),
        (0, 1),
        (0, 2),
    )
    assert unlabeled["px_per_bp"][0] > 0
    assert unlabeled["px_per_bp"] == pytest.approx(
        (unlabeled["px_per_bp"][0],) * 3
    )
    assert labeled["px_per_bp"] == pytest.approx(unlabeled["px_per_bp"])
    assert labeled["sequence_widths"] == pytest.approx(
        unlabeled["sequence_widths"]
    )
    unlabeled_x = unlabeled["record_x"]
    labeled_x = labeled["record_x"]
    assert tuple(x - labeled_x[0] for x in labeled_x) == pytest.approx(
        tuple(x - unlabeled_x[0] for x in unlabeled_x)
    )

    labeled_root = labeled["root"]
    definition_groups = tuple(
        group
        for group in labeled_root.findall(".//svg:g", namespace)
        if group.attrib.get("data-gbdraw-role")
        in {"record-definition", "record-definition-row"}
    )
    assert definition_groups
    assert all(
        {"data-record-row", "data-record-column"}.isdisjoint(group.attrib)
        for group in definition_groups
    )
    row_definition = next(
        group
        for group in definition_groups
        if group.attrib.get("data-gbdraw-role") == "record-definition-row"
    )
    row_text = row_definition.findall("./svg:text", namespace)
    assert ["".join(text.itertext()) for text in row_text] == [
        "Species alpha",
        "strain beta",
    ]
    assert all(
        text.attrib.get("text-anchor") == "start" and float(text.attrib["x"]) == 0.0
        for text in row_text
    )

    composition = json.loads(labeled_root.attrib["data-gbdraw-composition"])
    primary_bounds = composition["primary"]["finalBounds"]
    edge_padding = composition["spacing"]["edgePaddingPx"]
    view_box_left, _view_box_top, view_box_width, _view_box_height = labeled[
        "view_box"
    ]
    view_box_right = view_box_left + view_box_width
    assert primary_bounds["x"] - view_box_left == pytest.approx(edge_padding)
    assert view_box_right - (
        primary_bounds["x"] + primary_bounds["width"]
    ) == pytest.approx(edge_padding)
    assert _translate(row_definition)[0] == pytest.approx(primary_bounds["x"])

    gutter_growth = view_box_width - unlabeled["view_box"][2]
    record_shift = labeled_x[0] - unlabeled_x[0]
    assert gutter_growth > 0
    assert gutter_growth == pytest.approx(record_shift)


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
    records[0].annotations["gbdraw_record_label"] = "TUMSAT-TG-2018"
    records[0].annotations["gbdraw_record_subtitle"] = "chromosome 1"
    records[1].annotations["gbdraw_record_label"] = "chromosome 2"
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

    assert text_values(first_row_definition) == ["TUMSAT-TG-2018", "chromosome 1"]
    assert text_values(first_local_definition) == []
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
