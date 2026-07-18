from __future__ import annotations

import re
import xml.etree.ElementTree as ET

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from gbdraw.api import (
    InMemoryRecordSource,
    LinearDiagramRequest,
    LinearMultiRecordOptions,
    RecordInput,
    RecordPresentation,
    assemble_linear_diagram_from_records,
    build_request_diagram,
)
from gbdraw.config.toml import load_config_toml
from gbdraw.exceptions import ValidationError
from gbdraw.layout.linear_multi_record import (
    LinearRecordMeasurement,
    RecordKey,
    record_pairs_between_adjacent_rows,
    resolve_record_row_positions,
    solve_linear_layout,
)


def _records(*lengths: int) -> list[SeqRecord]:
    records = [SeqRecord(Seq("A" * length), id=f"record_{index}") for index, length in enumerate(lengths, 1)]
    for record in records:
        record.annotations["molecule_type"] = "DNA"
    return records


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


def test_comparison_anchor_can_overlay_reserved_record_header() -> None:
    measurements = (
        LinearRecordMeasurement(
            0,
            RecordKey("top"),
            100,
            bottom_extent=12,
        ),
        LinearRecordMeasurement(
            1,
            RecordKey("bottom"),
            100,
            top_extent=40,
            comparison_top_extent=18,
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
    assert bottom.axis_y == pytest.approx(122)
    assert bottom.comparison_top_y == pytest.approx(104)
    assert plan.content_top == pytest.approx(50)
    assert plan.content_bottom == pytest.approx(122)


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
        layout=LinearMultiRecordOptions(
            record_gap_px=24,
            multi_record_positions=("#1@1", "#2@1", "#3@2", "#4@2"),
        ),
        legend="none",
        config_overrides={"show_labels": False, "show_gc": False, "show_skew": False},
    )
    svg = canvas.tostring()
    assert svg.count("data-record-row=") == 4
    assert 'data-record-row="0"' in svg
    assert 'data-record-column="1"' in svg


def test_bottom_legend_follows_last_resolved_row() -> None:
    records = _records(*(1000 for _index in range(10)))
    svg = assemble_linear_diagram_from_records(
        records,
        layout=LinearMultiRecordOptions(
            multi_record_positions=tuple(
                f"#{index + 1}@{1 if index < 5 else 2}"
                for index in range(10)
            ),
        ),
        legend="bottom",
        plot_title="Resolved rows",
        plot_title_position="bottom",
        config_overrides={
            "show_labels": False,
            "show_gc": False,
            "show_skew": False,
        },
    ).tostring()
    root = ET.fromstring(svg)
    namespace = {"svg": "http://www.w3.org/2000/svg"}
    groups = {
        group.attrib["id"]: group
        for group in root.findall(".//svg:g", namespace)
        if "id" in group.attrib
    }

    def translate_y(group: ET.Element) -> float:
        match = re.fullmatch(
            r"translate\(([-+0-9.eE]+),([-+0-9.eE]+)\)",
            group.attrib["transform"],
        )
        assert match is not None
        return float(match.group(2))

    last_row_axis = translate_y(groups["record_10_record_10"])
    legend_top = translate_y(groups["legend"])
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
    config_dict["canvas"]["show_labels"] = False
    config_dict["canvas"]["linear"]["track_layout"] = "above"
    config_dict["canvas"]["linear"]["keep_definition_left_aligned"] = True
    config_dict["canvas"]["linear"]["ruler_on_axis"] = ruler_on_axis
    definition_cfg = config_dict["objects"]["definition"]["linear"]
    definition_cfg["show_replicon"] = False
    definition_cfg["show_accession"] = False
    definition_cfg["show_length"] = False

    svg = assemble_linear_diagram_from_records(
        records,
        layout=LinearMultiRecordOptions(
            record_gap_px=24,
            multi_record_positions=("#1@1", "#2@1"),
        ),
        legend="bottom",
        config_dict=config_dict,
    ).tostring()
    root = ET.fromstring(svg)
    namespace = {"svg": "http://www.w3.org/2000/svg"}
    groups = {
        group.attrib["id"]: group
        for group in root.findall(".//svg:g", namespace)
        if "id" in group.attrib
    }

    first_record = groups["record_1_record_1"]
    second_record = groups["record_2_record_2"]
    first_local_definition = groups["record_1_definition_record_1"]
    first_row_definition = groups["record_1_definition_record_1_row"]
    second_local_definition = groups["record_2_definition_record_2"]

    def text_values(group: ET.Element) -> list[str]:
        return [
            "".join(text.itertext())
            for text in group.findall(".//svg:text", namespace)
        ]

    def translate(group: ET.Element) -> tuple[float, float]:
        match = re.fullmatch(
            r"translate\(([-+0-9.eE]+),([-+0-9.eE]+)\)",
            group.attrib["transform"],
        )
        assert match is not None
        return float(match.group(1)), float(match.group(2))

    assert text_values(first_row_definition) == ["TUMSAT-TG-2018", "chromosome 1"]
    assert text_values(first_local_definition) == []
    assert text_values(second_local_definition) == ["chromosome 2"]
    assert 0 <= translate(first_row_definition)[0] < translate(first_record)[0]

    first_feature_y_values = [
        float(y_value)
        for path in first_record.findall(".//svg:path", namespace)
        for _x_value, y_value in re.findall(
            r"[ML]\s*([-+0-9.eE]+)\s*,?\s*([-+0-9.eE]+)",
            path.attrib.get("d", ""),
        )
    ]
    assert first_feature_y_values
    expected_definition_center_y = translate(first_record)[1] + 0.5 * (
        min(first_feature_y_values) + max(first_feature_y_values)
    )
    assert translate(first_row_definition)[1] == pytest.approx(
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
    feature_top = translate(second_record)[1] + min(feature_y_values)
    local_text = second_local_definition.find(".//svg:text", namespace)
    assert local_text is not None
    local_font_size = float(local_text.attrib["font-size"])
    local_bottom = translate(second_local_definition)[1] + (0.5 * local_font_size)
    assert local_bottom <= feature_top
    assert feature_top - local_bottom >= (
        float(config_dict["canvas"]["linear"]["vertical_padding"]) - 0.5
    )


def test_multi_record_layout_rejects_normalize_length() -> None:
    with pytest.raises(ValidationError, match="normalize_length"):
        assemble_linear_diagram_from_records(
            _records(100, 100),
            layout=LinearMultiRecordOptions(
                multi_record_positions=("#1@1", "#2@1"),
            ),
            legend="none",
            config_overrides={"normalize_length": True},
        )


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
