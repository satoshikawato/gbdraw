"""Several windows of one record in one Linear row (D-01 of the multi-contig-locus campaign).

A record keeps one identity and lists its windows; each window is drawn as its own unit, and the windows of one
record share a row, ordered by column, between other records' units.
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pytest

import re

from gbdraw.annotations import AnnotationSet, CoordinateSpan, RegionAnnotation, resolve_annotations
from gbdraw.api.options import LinearDiagramOptions
from gbdraw.features.ids import feature_instances
from gbdraw.features.overrides import FeatureOverride
from gbdraw.features.placement import FeaturePlacementOverride, FeaturePlacementTarget
from gbdraw.api.record_planning import record_input_manifest_from_paths, record_input_manifest_from_table
from gbdraw.api.request_render import plan_linear_request, resolve_request_records
from gbdraw.api.requests import (
    CircularDiagramRequest,
    GenBankInputSource,
    LinearDiagramRequest,
    LinearMultiRecordOptions,
    RecordCardinality,
    RecordDisplayOptions,
    RecordInput,
    RecordPresentation,
)
from gbdraw.api import render_request
from gbdraw.core.record_metadata import _source_feature_index, record_identity_key
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
    # D-20: read back, the region is the one window named by the record's selector.
    assert record.region == parse_region_spec("TESTB:5-40:rc")


def test_region_reads_back_none_without_windows_and_refuses_several() -> None:
    assert RecordInput(GenBankInputSource(MITO)).region is None
    assert RecordInput(GenBankInputSource(MITO), None, parse_region_spec("100-900")).region == parse_region_spec(
        "100-900"
    )
    with pytest.raises(ValidationError, match="windows"):
        _ = _windowed(RecordWindow(100, 900), RecordWindow(1000, 1100)).region


def test_replace_keeps_the_window_of_a_record_built_from_a_region() -> None:
    record = RecordInput(GenBankInputSource(TWO_RECORDS), region=parse_region_spec("TESTB:5-40:rc"), record_key="k")
    relabelled = replace(record, presentation=RecordPresentation(label="B"))
    assert relabelled.windows == record.windows and relabelled.selector == record.selector
    assert relabelled.region == record.region
    assert replace(record, windows=()).region is None
    assert record.oriented(False).windows == (RecordWindow(5, 40, False),)


@pytest.mark.parametrize(
    ("make", "message"),
    [
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


_GROUPED = "gbk\trecord_label\trow\tcolumn\tregion\treverse_complement\trecord_group"


def _table(tmp_path: Path, rows: list[str], header: str = _GROUPED) -> Path:
    table = tmp_path / "records.tsv"
    table.write_text(header + "\n" + "".join(f"{row}\n" for row in rows), encoding="utf-8")
    return table


def test_records_table_rows_sharing_a_record_group_are_the_windows_of_one_record(tmp_path: Path) -> None:
    # D-32: rows are windows of one record only when they share a record_group value; the record's label is the
    # one its rows give.
    table = _table(tmp_path, [
        f"{MITO}\tM\t1\t1\t100-900\tfalse\tmito",
        f"{MITO}\tT\t1\t2\t\tfalse\t",
        f"{MITO}\t\t1\t3\t5000-5600\ttrue\tmito",
        f"{MITO}\tM\t2\t1\t1-50\tfalse\t",
    ])
    manifest = record_input_manifest_from_table(str(table))
    assert [record.record_key for record in manifest.records] == ["record-1", "record-2", "record-3"]
    assert manifest.records[0].windows == (RecordWindow(100, 900, False, 1), RecordWindow(5000, 5600, True, 3))
    assert manifest.records[0].presentation.label == "M"
    assert manifest.records[2].windows == (RecordWindow(1, 50),)
    request = LinearDiagramRequest(records=manifest.records, layout=LinearMultiRecordOptions())
    plan = plan_linear_request(request)
    assert plan.layout is not None
    # Record 1's windows are units 1 and 2, record 2 is unit 3, record 3 unit 4.
    assert plan.layout.multi_record_positions == ("#1@1", "#3@1", "#2@1", "#4@2")


_V0145 = "gbk\trecord_label\ttopology\trow\tcolumn\tregion"


def _v0145_table(tmp_path: Path) -> Path:
    # v0.14.5 drew these rows as two records side by side, each with its own label and topology.
    return _table(
        tmp_path,
        [f"{MITO}\tA\tlinear\t1\t1\t1000-6000", f"{MITO}\tB\tcircular\t1\t2\t4000-9000:rc"],
        _V0145,
    )


def test_a_v0145_records_table_naming_one_record_twice_in_a_row_draws_two_records(tmp_path: Path) -> None:
    manifest = record_input_manifest_from_table(str(_v0145_table(tmp_path)))
    # The inputs v0.14.5 built from these rows: one record per row, each cropped by its region.
    assert manifest.records == tuple(
        RecordInput(
            source=GenBankInputSource(str(MITO)),
            cardinality=RecordCardinality.EXACTLY_ONE,
            selector=None,
            display=RecordDisplayOptions(is_circular=circular),
            region=parse_region_spec(region),
            presentation=RecordPresentation(label=label, grid_row=1, grid_column=column),
            record_key=f"record-{column}",
        )
        for label, circular, column, region in (("A", False, 1, "1000-6000"), ("B", True, 2, "4000-9000:rc"))
    )
    collection = resolve_request_records(
        LinearDiagramRequest(records=manifest.records, layout=LinearMultiRecordOptions())
    )
    assert [(item.record_key, item.record_identity_key) for item in collection.provenance] == [
        ("record-1", "record-1"), ("record-2", "record-2"),
    ]


def test_companion_bindings_of_a_v0145_records_table_name_the_rows(tmp_path: Path) -> None:
    manifest = record_input_manifest_from_table(str(_v0145_table(tmp_path)))
    request = LinearDiagramRequest(
        records=manifest.records,
        layout=LinearMultiRecordOptions(),
        options=LinearDiagramOptions(
            feature_overrides=(FeatureOverride("record-2", _IN_OVERLAP, label_visibility="on", label_text="EDITED"),),
        ),
    )
    plan = plan_linear_request(request)
    assert plan.inputs is not None
    assert [[o.label_text for o in unit.overrides.values()] for unit in plan.inputs.record_features] == [
        [], ["EDITED"],
    ]
    records = resolve_request_records(request).records
    bundle = resolve_annotations(
        (AnnotationSet("s", (RegionAnnotation("a", CoordinateSpan(RecordSelector("#2", None, 1), 4500, 5000)),)),),
        records,
        mode="linear",
    )
    assert [item.record_index for item in bundle.annotations] == [1]


@pytest.mark.parametrize(
    ("rows", "message"),
    [
        ([f"{MITO}\tM\t1\t1\t100-900\tfalse\tg", f"{MITO}\tM\t2\t1\t1000-1100\tfalse\tg"], "same row"),
        ([f"{MITO}\tM\t1\t1\t100-900\tfalse\tg", f"{TWO_RECORDS}\tM\t1\t2\t1-10\tfalse\tg"], "same gbk"),
        ([f"{MITO}\tM\t1\t1\t100-900\tfalse\tg", f"{MITO}\tM\t1\t2\t\tfalse\tg"], "region"),
        ([f"{MITO}\tM\t1\t1\t100-900\tfalse\tg", f"{MITO}\tN\t1\t2\t1000-1100\tfalse\tg"],
         "different record_label"),
    ],
    ids=["row", "file", "region", "label"],
)
def test_rows_of_one_record_group_that_disagree_are_refused(tmp_path: Path, rows: list[str], message: str) -> None:
    with pytest.raises(ValidationError, match=message):
        record_input_manifest_from_table(str(_table(tmp_path, rows)))


def test_rows_of_one_record_group_that_disagree_on_the_record_id_or_topology_are_refused(tmp_path: Path) -> None:
    header = "gbk\trecord_id\ttopology\trow\tcolumn\tregion\trecord_group"
    for rows, message in (
        ([f"{TWO_RECORDS}\tTESTA\t\t1\t1\t1-10\tg", f"{TWO_RECORDS}\tTESTB\t\t1\t2\t1-10\tg"], "same record_id"),
        ([f"{MITO}\t\tlinear\t1\t1\t1-10\tg", f"{MITO}\t\tcircular\t1\t2\t20-30\tg"], "different topology"),
    ):
        with pytest.raises(ValidationError, match=message):
            record_input_manifest_from_table(str(_table(tmp_path, rows, header)))


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
    rendered = render_request(replace(request, output=replace(request.output, output_directory=tmp_path)))
    [svg] = [path for path in rendered.output_paths if path.suffix == ".svg"]
    text = svg.read_text(encoding="utf-8")
    assert text.count('data-record-key="mito~w') >= 2 or "mito~w2" in text


# A feature of HmmtDNA inside 4000-6000 (ND2, 4469..5511) and one across 6500-7000 (COX1, 5903..7445).
_IN_OVERLAP = "fbe3a7c0c"
_ACROSS_GAP = "f48b7bf2f"


def test_overlapping_windows_of_a_record_are_drawn_both() -> None:
    # D-18: windows of one record may overlap; each draws its own span.
    collection = resolve_request_records(
        LinearDiagramRequest(records=(_windowed(RecordWindow(1000, 6000), RecordWindow(4000, 9000)),))
    )
    assert [len(record) for record in collection.records] == [5001, 5001]


def test_a_records_table_with_overlapping_crops_of_one_record_in_one_row_draws_both(tmp_path: Path) -> None:
    # Two rows of one record_group draw as two windows of one record.
    table = _table(tmp_path, [f"{MITO}\tM\t1\t1\t1000-6000\tfalse\tg", f"{MITO}\tM\t1\t2\t4000-9000\tfalse\tg"])
    manifest = record_input_manifest_from_table(str(table))
    request = LinearDiagramRequest(
        records=manifest.records, layout=LinearMultiRecordOptions(),
    )
    rendered = render_request(replace(request, output=replace(request.output, output_directory=tmp_path)))
    [svg] = [path for path in rendered.output_paths if path.suffix == ".svg"]
    text = svg.read_text(encoding="utf-8")
    ids = re.findall(r"\sid=\"([^\"]+)\"", text)
    assert len(ids) == len(set(ids))
    # ND2 lies in both crops: drawn twice, once in each unit.
    assert len(re.findall(rf'data-gbdraw-stable-feature-id="{_IN_OVERLAP}"', text)) >= 2


@pytest.mark.parametrize(
    ("windows", "feature"),
    [
        ((RecordWindow(1000, 6000), RecordWindow(4000, 9000)), _IN_OVERLAP),
        ((RecordWindow(1000, 6500), RecordWindow(7000, 9000)), _ACROSS_GAP),
    ],
    ids=["inside-an-overlap", "across-a-left-out-interval"],
)
def test_a_feature_drawn_in_two_windows_is_one_instance_edited_once(windows, feature: str) -> None:
    request = LinearDiagramRequest(
        records=(_windowed(*windows),),
        options=LinearDiagramOptions(
            feature_overrides=(FeatureOverride("mito", feature, label_visibility="on", label_text="EDITED"),),
            feature_placements=(FeaturePlacementOverride("mito", feature, FeaturePlacementTarget("lane", "above", 1)),),
        ),
    )
    plan = plan_linear_request(request)
    assert plan.inputs is not None
    edited = [
        [override.label_text for override in unit.overrides.values()] for unit in plan.inputs.record_features
    ]
    assert edited == [["EDITED"], ["EDITED"]]
    placed = [[item.biological_feature_id for item in unit.foreground] for unit in plan.inputs.record_features]
    assert placed == [[feature], [feature]]
    assert not plan.inputs.feature_identity_notices
    collection = resolve_request_records(request)
    [source_index] = [
        entry.source_feature_index
        for entry in collection.provenance[0].source_feature_catalog
        if entry.biological_feature_id == feature
    ]
    instances = {
        feature_instances(record)(item)
        for record in collection.records
        for item in record.features
        if _source_feature_index(item) == source_index
    }
    assert instances == {f"{feature}@mito"}


def test_an_annotation_of_a_record_is_drawn_in_each_window_it_falls_in() -> None:
    records = resolve_request_records(
        LinearDiagramRequest(records=(_windowed(RecordWindow(1000, 6000), RecordWindow(4000, 9000)),))
    ).records
    selector = RecordSelector("#1", None, 0)
    bundle = resolve_annotations(
        (AnnotationSet("s", (RegionAnnotation("a", CoordinateSpan(selector, 4500, 5000), label="A"),)),),
        records,
        mode="linear",
    )
    assert [item.record_index for item in bundle.annotations] == [0, 1]
    assert sorted(item.label for item in bundle.annotations) == ["", "A"]
    assert not bundle.warnings
    # A span across the left-out interval is clipped to each window it falls in, as one annotation.
    gapped = resolve_request_records(
        LinearDiagramRequest(records=(_windowed(RecordWindow(1000, 3000), RecordWindow(5000, 9000)),))
    ).records
    bundle = resolve_annotations(
        (AnnotationSet("s", (RegionAnnotation("a", CoordinateSpan(selector, 2500, 5500)),)),),
        gapped,
        mode="linear",
    )
    assert [(item.record_index, item.span_bp) for item in bundle.annotations] == [(0, 501), (1, 501)]
    assert [warning.code for warning in bundle.warnings] == ["out_of_bounds_clipped"]


def test_a_wrap_never_cuts_the_windows_of_a_record() -> None:
    request = LinearDiagramRequest(
        records=(_windowed(RecordWindow(1, 6000), RecordWindow(8000, 16000)),),
        options=LinearDiagramOptions(config_overrides={"canvas.linear.wrap_length": 2000}),
    )
    with pytest.raises(ValidationError, match="several windows of one record"):
        plan_linear_request(request)


def test_each_axis_entry_names_its_record_and_window(tmp_path: Path) -> None:
    request = LinearDiagramRequest(records=(_windowed(RecordWindow(100, 900), RecordWindow(5000, 5600)),))
    geometry = plan_linear_request(request).build().drawing._gbdraw_track_slot_geometry["records"]
    assert [(entry["recordKey"], entry.get("windowIndex")) for entry in geometry] == [("mito", 0), ("mito", 1)]


# D-28, D-34: a window of a circular record may run across its origin (start after end). HmmtDNA (16,569 bp) has its
# D-loop across the origin: complement(join(16024..16569,1..576)).
_ACROSS = RecordWindow(16000, 1000)


def _source_record():
    from Bio import SeqIO

    return SeqIO.read(MITO, "genbank")


@pytest.mark.parametrize("reverse", [False, True], ids=["forward", "reversed"])
def test_a_window_across_the_origin_of_a_circular_record_is_one_unit(reverse: bool) -> None:
    from gbdraw.core.record_metadata import _record_source_start

    assert _ACROSS.crosses_origin
    collection = resolve_request_records(LinearDiagramRequest(records=(_windowed(replace(_ACROSS, reverse_complement=reverse)),)))
    [unit] = collection.records
    source = _source_record()
    arc = source.seq[15999:] + source.seq[:1000]
    assert str(unit.seq) == str(arc.reverse_complement() if reverse else arc)
    # The feature across the origin stays whole, one part, and keeps its source identity.
    [dloop] = [feature for feature in unit.features if feature.type == "D-loop"]
    assert len(dloop.location.parts) == 1
    assert (int(dloop.location.start), int(dloop.location.end)) == ((424, 1546) if reverse else (24, 1146))
    [source_dloop] = [feature for feature in source.features if feature.type == "D-loop"]
    assert getattr(dloop, "_gbdraw_source_feature_location_parts") == tuple(
        (int(part.start), int(part.end), part.strand) for part in source_dloop.location.parts
    )
    # The crop map names source bases: the window's source start is its start, also reversed.
    assert _record_source_start(unit) == 16000
    assert (unit.annotations["gbdraw_coord_base"], unit.annotations["gbdraw_coord_step"]) == (
        (1000, -1) if reverse else (16000, 1)
    )


def test_a_window_across_the_origin_of_a_linear_record_is_refused() -> None:
    record = RecordInput(GenBankInputSource(MITO), windows=(_ACROSS,), display=RecordDisplayOptions(is_circular=False))
    with pytest.raises(ValidationError, match="crosses the origin of a linear record"):
        resolve_request_records(LinearDiagramRequest(records=(record,)))
    with pytest.raises(ValidationError, match="across the origin"):
        _ = record.region


def test_the_ruler_and_caption_of_a_window_across_the_origin_name_source_bases(tmp_path: Path) -> None:
    from gbdraw.layout.record_labels import coordinate_map_ticks

    [unit] = resolve_request_records(LinearDiagramRequest(records=(_windowed(_ACROSS),))).records
    ticks = coordinate_map_ticks(unit, len(unit), 500)
    # 16500 before the origin; the multiples restart at the origin: 500 is 570 + 499 bases in.
    assert [tick.base for tick in ticks] == [16500, 500]
    assert [round(tick.position * 1569) for tick in ticks] == [500, 1069]
    [reversed_unit] = resolve_request_records(
        LinearDiagramRequest(records=(_windowed(replace(_ACROSS, reverse_complement=True)),))
    ).records
    assert [tick.base for tick in coordinate_map_ticks(reversed_unit, len(reversed_unit), 500)] == [500, 16500]
    request = LinearDiagramRequest(records=(_windowed(_ACROSS),))
    rendered = render_request(replace(request, output=replace(request.output, output_directory=tmp_path)))
    [svg] = [path for path in rendered.output_paths if path.suffix == ".svg"]
    text = svg.read_text(encoding="utf-8")
    assert 'data-gbdraw-record-source-start="16000"' in text
    assert 'data-gbdraw-record-source-end="1000"' in text


def test_the_caption_of_a_window_across_the_origin_runs_from_its_start_to_its_end() -> None:
    from gbdraw.layout.record_labels import record_map_span
    from gbdraw.core.record_metadata import _source_period, _wrap_source_base

    [unit] = resolve_request_records(LinearDiagramRequest(records=(_windowed(_ACROSS),))).records
    period = _source_period(unit)
    assert period == 16569
    assert tuple(_wrap_source_base(value, period) for value in record_map_span(unit, len(unit))) == (16000, 1000)


def test_an_annotation_in_a_window_across_the_origin_is_drawn_at_its_source_bases() -> None:
    records = resolve_request_records(LinearDiagramRequest(records=(_windowed(_ACROSS),))).records
    selector = RecordSelector("#1", None, 0)
    bundle = resolve_annotations(
        (AnnotationSet("s", (
            RegionAnnotation("after", CoordinateSpan(selector, 100, 300)),
            RegionAnnotation("before", CoordinateSpan(selector, 16201, 16300)),
            RegionAnnotation("cut", CoordinateSpan(selector, 900, 1100)),
        )),),
        records,
        mode="linear",
    )
    assert {item.id: item.segments for item in bundle.annotations} == {
        "after": ((669, 870),), "before": ((201, 301),), "cut": ((1469, 1570),),
    }
    assert [warning.code for warning in bundle.warnings] == ["out_of_bounds_clipped"]


# I-6 (D-34): `#n` counts records in the given row positions and in a comparisons table, not drawn windows.


def test_given_row_positions_name_records_and_place_all_their_windows() -> None:
    request = LinearDiagramRequest(
        records=(
            _windowed(RecordWindow(100, 900), RecordWindow(5000, 5600)),
            RecordInput(GenBankInputSource(MITO), record_key="other", region=parse_region_spec("1-500")),
        ),
        layout=LinearMultiRecordOptions(multi_record_positions=("#2@1", "#1@2")),
    )
    plan = plan_linear_request(request)
    assert plan.layout is not None
    assert plan.layout.multi_record_positions == ("#3@1", "#1@2", "#2@2")
    with pytest.raises(ValidationError, match="out of range"):
        plan_linear_request(replace(request, layout=LinearMultiRecordOptions(multi_record_positions=("#3@1",))))


def test_a_comparisons_table_names_a_windowed_record_by_its_record_number(tmp_path: Path) -> None:
    from gbdraw.api.record_planning import resolve_linear_options

    records = resolve_request_records(
        LinearDiagramRequest(
            records=(
                _windowed(RecordWindow(100, 900), RecordWindow(5000, 5600)),
                RecordInput(GenBankInputSource(MITO), record_key="other", region=parse_region_spec("1-500")),
            ),
        )
    ).records
    layout = LinearMultiRecordOptions(multi_record_positions=("#1@1", "#2@1", "#3@2"))
    blast = tmp_path / "pair.tsv"
    blast.write_text("q\ts\t90\t100\t0\t0\t1\t100\t1\t100\t1e-30\t200\n", encoding="utf-8")
    table = tmp_path / "comparisons.tsv"
    table.write_text(f"blast\tquery\tsubject\n{blast}\t#2\t#1\n", encoding="utf-8")
    with pytest.raises(ValidationError, match="drawn in several windows"):
        resolve_linear_options(LinearDiagramOptions(comparison_table_file=str(table)), records=records, layout=layout)
    # `#2` is the second record (the third drawn unit), not the first record's second window.
    table.write_text(f"blast\tquery\tsubject\n{blast}\t#2\t#3\n", encoding="utf-8")
    with pytest.raises(ValidationError, match="unresolved record selector '#3'"):
        resolve_linear_options(LinearDiagramOptions(comparison_table_file=str(table)), records=records, layout=layout)


def _saved_request(session: Path) -> dict:
    import json

    document = json.loads(session.read_text(encoding="utf-8"))
    return document["drawings"][0]["renderRequest"]


def _cli(monkeypatch: pytest.MonkeyPatch, arguments: list[str]) -> None:
    import sys

    from gbdraw.cli import main

    monkeypatch.setattr(sys, "argv", ["gbdraw", *arguments])
    main()


@pytest.mark.parametrize(
    "regions",
    [["100-900", "5000-5600:rc"], ["16000-200"]],
    ids=["two-windows", "one-reversed-crop"],
)
def test_a_cli_session_saves_a_record_with_its_windows_and_replays_the_figure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, regions: list[str]
) -> None:
    def main(arguments: list[str]) -> None:
        _cli(monkeypatch, arguments)

    monkeypatch.chdir(tmp_path)
    arguments = ["linear", "--gbk", str(MITO), "-o", "w", "-f", "svg", "--save_session"]
    for region in regions:
        arguments += ["--region", region]
    main(arguments)
    [record] = _saved_request(tmp_path / "w.gbdraw-session.json")["records"]
    assert record["recordKey"] == "record-1"
    assert record["source"]["kind"] == "genbank" and record["selector"] is None
    expected = [parse_region_spec(region) for region in regions]
    assert [(w["start"], w["end"], w["reverseComplement"]) for w in record["windows"]] == [
        (spec.start, spec.end, spec.reverse_complement) for spec in expected
    ]
    main(["render", "--session", "w.gbdraw-session.json", "-o", "replay", "-f", "svg"])
    assert (tmp_path / "replay.svg").read_bytes() == (tmp_path / "w.svg").read_bytes()


def test_a_cli_session_of_a_records_table_group_keeps_the_window_columns(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.chdir(tmp_path)
    table = _table(
        tmp_path,
        [
            f"{MITO}\t\tM\t1\t1\t100-900\tg",
            f"{TWO_RECORDS}\tTESTA\tT\t1\t2\t\t",
            f"{MITO}\t\t\t1\t3\t5000-5600:rc\tg",
            f"{TWO_RECORDS}\tTESTB\tN\t2\t1\t\t",
        ],
        "gbk\trecord_id\trecord_label\trow\tcolumn\tregion\trecord_group",
    )
    _cli(monkeypatch, ["linear", "--records_table", str(table), "-o", "g", "-f", "svg", "--save_session"])
    records = _saved_request(tmp_path / "g.gbdraw-session.json")["records"]
    assert [(record["recordKey"], [w["column"] for w in record["windows"]]) for record in records] == [
        ("record-1", [1, 3]), ("record-2", []), ("record-3", []),
    ]
    _cli(monkeypatch, ["render", "--session", "g.gbdraw-session.json", "-o", "replay", "-f", "svg"])
    assert (tmp_path / "replay.svg").read_bytes() == (tmp_path / "g.svg").read_bytes()
