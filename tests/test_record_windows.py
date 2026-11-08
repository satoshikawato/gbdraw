"""Several windows of one record in one Linear row (D-01 of the multi-contig-locus campaign).

A record keeps one identity and lists its windows; each window is drawn as its own unit, and the windows of one
record share a row, ordered by column, between other records' units.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from gbdraw.api.record_planning import record_input_manifest_from_paths, record_input_manifest_from_table
from gbdraw.api.request_render import plan_linear_request, resolve_request_records
from gbdraw.api.requests import (
    CircularDiagramRequest,
    GenBankInputSource,
    LinearDiagramRequest,
    LinearMultiRecordOptions,
    RecordCardinality,
    RecordInput,
    RecordPresentation,
)
from gbdraw.api import render_request
from gbdraw.core.record_metadata import record_identity_key
from gbdraw.exceptions import ValidationError
from gbdraw.io.record_select import RecordSelector
from gbdraw.io.regions import RecordWindow, parse_region_spec
from gbdraw.session_request_codec import (
    CANONICAL_REQUEST_SCHEMA,
    decode_canonical_request,
    encode_canonical_request,
)

REPO = Path(__file__).resolve().parents[1]
MITO = REPO / "tests" / "test_inputs" / "HmmtDNA.gbk"
TWO_RECORDS = REPO / "tests" / "fixtures" / "web_batch_two_records.gb"


def _windowed(*windows: RecordWindow, **presentation: object) -> RecordInput:
    return RecordInput(
        GenBankInputSource(MITO),
        windows=windows,
        record_key="mito",
        presentation=RecordPresentation(**presentation),  # type: ignore[arg-type]
    )


def test_each_window_of_a_record_is_a_drawn_unit_of_that_record() -> None:
    request = LinearDiagramRequest(records=(_windowed(RecordWindow(100, 900), RecordWindow(5000, 5600, True)),))
    collection = resolve_request_records(request)
    assert [len(record) for record in collection.records] == [801, 601]
    assert [item.record_key for item in collection.provenance] == ["mito~w1", "mito~w2"]
    assert {item.record_identity_key for item in collection.provenance} == {"mito"}
    assert [record_identity_key(record) for record in collection.records] == ["mito", "mito"]
    # The second window is drawn reversed: its first base is the window's last.
    assert [item.window.span.reverse_complement for item in collection.provenance if item.window] == [False, True]


def test_a_record_with_one_window_draws_as_its_region_did() -> None:
    windowed = resolve_request_records(LinearDiagramRequest(records=(_windowed(RecordWindow(100, 900)),)))
    region = resolve_request_records(
        LinearDiagramRequest(records=(RecordInput(GenBankInputSource(MITO), region=parse_region_spec("100-900"),
                                                  record_key="mito"),))
    )
    assert [item.record_key for item in windowed.provenance] == ["mito"]
    assert str(windowed.records[0].seq) == str(region.records[0].seq)
    assert "gbdraw_window_of" not in windowed.records[0].annotations


def test_the_windows_of_a_record_share_one_row_without_grid_rows() -> None:
    request = LinearDiagramRequest(
        records=(
            _windowed(RecordWindow(100, 900), RecordWindow(5000, 5600)),
            RecordInput(GenBankInputSource(MITO), record_key="other", region=parse_region_spec("1-500")),
        )
    )
    plan = plan_linear_request(request)
    assert plan.layout is not None
    assert plan.layout.multi_record_positions == ("#1@1", "#2@1", "#3@2")


def test_windows_sit_between_other_records_by_column() -> None:
    request = LinearDiagramRequest(
        records=(
            _windowed(RecordWindow(100, 900, column=1), RecordWindow(5000, 5600, column=3), grid_row=1),
            RecordInput(GenBankInputSource(MITO), record_key="other", region=parse_region_spec("1-500"),
                        presentation=RecordPresentation(grid_row=1, grid_column=2)),
        ),
        layout=LinearMultiRecordOptions(),
    )
    plan = plan_linear_request(request)
    assert plan.layout is not None
    assert plan.layout.multi_record_positions == ("#1@1", "#3@1", "#2@1")


def test_a_region_is_read_as_one_window_with_its_selector() -> None:
    record = RecordInput(GenBankInputSource(TWO_RECORDS), region=parse_region_spec("TESTB:5-40:rc"))
    assert record.windows == (RecordWindow(5, 40, True),)
    assert record.selector == RecordSelector("TESTB", "TESTB", None)
    assert record.region is None


@pytest.mark.parametrize(
    ("make", "message"),
    [
        (lambda: _windowed(RecordWindow(100, 900), RecordWindow(900, 1000)), "overlap"),
        (lambda: _windowed(RecordWindow(100, 900), RecordWindow(1000, 1100), reverse_complement=True),
         "each window"),
        (lambda: _windowed(RecordWindow(100, 900, column=2)), "one window is placed by its presentation"),
        (lambda: _windowed(RecordWindow(100, 900), RecordWindow(1000, 1100), grid_row=1, grid_column=1),
         "own column"),
        (lambda: CircularDiagramRequest(records=(_windowed(RecordWindow(100, 900), RecordWindow(1000, 1100)),)),
         "Circular"),
    ],
)
def test_windows_that_cannot_be_drawn_are_refused(make, message: str) -> None:
    with pytest.raises(ValidationError, match=message):
        make()


def test_repeated_regions_of_one_record_are_its_windows_on_the_command_line() -> None:
    manifest = record_input_manifest_from_paths(
        gbk_paths=[str(MITO)],
        gff_paths=None,
        fasta_paths=None,
        cardinalities=RecordCardinality.EXACTLY_ONE,
        selectors=None,
        reverse_flags=None,
        labels=None,
        subtitles=None,
        regions=["100-900", "5000-5600:rc"],
    )
    request = LinearDiagramRequest(records=manifest.records, record_options=manifest.record_options)
    collection = resolve_request_records(request)
    assert [item.record_key for item in collection.provenance] == ["record-1~w1", "record-1~w2"]
    assert plan_linear_request(request).layout.multi_record_positions == ("#1@1", "#2@1")  # type: ignore[union-attr]


def _table(tmp_path: Path, rows: list[str]) -> Path:
    table = tmp_path / "records.tsv"
    table.write_text(
        "gbk\trecord_label\trow\tcolumn\tregion\treverse_complement\n" + "".join(f"{row}\n" for row in rows),
        encoding="utf-8",
    )
    return table


def test_records_table_rows_of_one_record_in_one_row_are_its_windows(tmp_path: Path) -> None:
    table = _table(tmp_path, [
        f"{MITO}\tM\t1\t1\t100-900\tfalse",
        f"{MITO}\tT\t1\t2\t\tfalse",
        f"{MITO}\tM\t1\t3\t5000-5600\ttrue",
        f"{MITO}\tM\t2\t1\t1-50\tfalse",
    ])
    manifest = record_input_manifest_from_table(str(table))
    assert [record.record_key for record in manifest.records] == ["record-1", "record-2", "record-3"]
    assert manifest.records[0].windows == (RecordWindow(100, 900, False, 1), RecordWindow(5000, 5600, True, 3))
    assert manifest.records[2].windows == (RecordWindow(1, 50),)
    request = LinearDiagramRequest(records=manifest.records, layout=LinearMultiRecordOptions())
    plan = plan_linear_request(request)
    assert plan.layout is not None
    # Record 1's windows are units 1 and 2, record 2 is unit 3, record 3 unit 4.
    assert plan.layout.multi_record_positions == ("#1@1", "#3@1", "#2@1", "#4@2")


def test_records_table_rows_of_one_record_must_agree_on_its_label(tmp_path: Path) -> None:
    table = _table(tmp_path, [f"{MITO}\tM\t1\t1\t100-900\tfalse", f"{MITO}\tN\t1\t2\t1000-1100\tfalse"])
    with pytest.raises(ValidationError, match="record_label"):
        record_input_manifest_from_table(str(table))


def test_request_schema_10_writes_windows_and_schema_9_reads_a_region_as_one(tmp_path: Path) -> None:
    request = LinearDiagramRequest(records=(_windowed(RecordWindow(100, 900), RecordWindow(5000, 5600, True)),))
    encoded = encode_canonical_request(request)
    assert encoded.payload["schema"] == CANONICAL_REQUEST_SCHEMA == 10
    [record] = encoded.payload["records"]
    assert "region" not in record
    assert record["windows"] == [
        {"start": 100, "end": 900, "reverseComplement": False, "column": None},
        {"start": 5000, "end": 5600, "reverseComplement": True, "column": None},
    ]
    resources = {
        resource.resource_id: resource.source_path or tmp_path / resource.name for resource in encoded.resources
    }
    decoded = decode_canonical_request(encoded.payload, resource_paths=resources, output_directory=tmp_path)
    assert decoded.records[0].windows == request.records[0].windows

    legacy = dict(encoded.payload, schema=9)
    legacy["records"] = [dict(record, region={"selector": None, "start": 100, "end": 900, "reverseComplement": True})]
    del legacy["records"][0]["windows"]
    decoded = decode_canonical_request(legacy, resource_paths=resources, output_directory=tmp_path)
    assert decoded.records[0].windows == (RecordWindow(100, 900, True),)


def test_a_row_of_windows_renders(tmp_path: Path) -> None:
    request = LinearDiagramRequest(
        records=(
            _windowed(RecordWindow(100, 900, column=1), RecordWindow(5000, 5600, True, column=3), grid_row=1),
            RecordInput(GenBankInputSource(TWO_RECORDS), record_key="pair",
                        selector=RecordSelector("#1", None, 0),
                        presentation=RecordPresentation(grid_row=1, grid_column=2)),
        ),
        layout=LinearMultiRecordOptions(),
    )
    from dataclasses import replace

    rendered = render_request(replace(request, output=replace(request.output, output_directory=tmp_path)))
    [svg] = [path for path in rendered.output_paths if path.suffix == ".svg"]
    text = svg.read_text(encoding="utf-8")
    assert text.count('data-record-key="mito~w') >= 2 or "mito~w2" in text
