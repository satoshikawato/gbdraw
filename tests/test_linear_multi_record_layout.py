from __future__ import annotations

import pytest
from Bio.Seq import Seq
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
from gbdraw.exceptions import ValidationError
from gbdraw.layout.linear_multi_record import (
    LinearRecordMeasurement,
    RecordKey,
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
