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
from gbdraw.core.record_metadata import _source_feature_index, _source_feature_location_parts, record_identity_key
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
    assert _source_feature_location_parts(dloop) == tuple(
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
    # Read back, the window is the region that writes it (``:origin``, PD-OI-122).
    assert record.region is not None and record.region.raw == "16000-1000:origin"


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
    from tests.utils.session_files import read_session_payload

    document = read_session_payload(session)
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
    [record] = _saved_request(tmp_path / "w.gbdraw")["records"]
    assert record["recordKey"] == "record-1"
    assert record["source"]["kind"] == "genbank" and record["selector"] is None
    expected = [parse_region_spec(region) for region in regions]
    assert [(w["start"], w["end"], w["reverseComplement"]) for w in record["windows"]] == [
        (spec.start, spec.end, spec.reverse_complement) for spec in expected
    ]
    main(["render", "--session", "w.gbdraw", "-o", "replay", "-f", "svg"])
    assert (tmp_path / "replay.svg").read_bytes() == (tmp_path / "w.svg").read_bytes()
    # A Session the Web saved before it wrote this form names the file's only record; it draws the same figure.
    from tests.utils.session_files import read_session_payload, write_session_file

    document = read_session_payload(tmp_path / "w.gbdraw")
    document["drawings"][0]["renderRequest"]["records"][0]["selector"] = {"kind": "recordId", "value": "NC_012920.1"}
    write_session_file(tmp_path / "named.gbdraw", document)
    main(["render", "--session", "named.gbdraw", "-o", "named", "-f", "svg"])
    assert (tmp_path / "named.svg").read_bytes() == (tmp_path / "w.svg").read_bytes()


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
    records = _saved_request(tmp_path / "g.gbdraw")["records"]
    assert [(record["recordKey"], [w["column"] for w in record["windows"]]) for record in records] == [
        ("record-1", [1, 3]), ("record-2", []), ("record-3", []),
    ]
    _cli(monkeypatch, ["render", "--session", "g.gbdraw", "-o", "replay", "-f", "svg"])
    assert (tmp_path / "replay.svg").read_bytes() == (tmp_path / "g.svg").read_bytes()


def test_a_records_table_names_each_record_of_a_file_the_web_draws_whole_in_a_placed_card(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """OV-3072: a records table writes a file drawn whole one row per record, each with its column; a Session of the
    Web's File card of that file, placed at the first of those columns, draws the same figure. A File card's records
    are keyed ``<card>:<n>``, the table's ``record-<row>``."""
    from tests.utils.session_files import read_session_payload, write_unchecked_session_file

    monkeypatch.chdir(tmp_path)
    table = _table(
        tmp_path,
        [f"{TWO_RECORDS}\t#1\t1\t1", f"{TWO_RECORDS}\t#2\t1\t2", f"{MITO}\t\t1\t3"],
        "gbk\trecord_id\trow\tcolumn",
    )
    _cli(monkeypatch, ["linear", "--records_table", str(table), "-o", "t", "-f", "svg", "--save_session"])
    document = read_session_payload(tmp_path / "t.gbdraw")
    request = document["drawings"][0]["renderRequest"]
    first, _second, last = request["records"]
    request["records"] = [
        {**first, "cardinality": "all", "selector": None},
        {**last, "presentation": {**last["presentation"], "gridColumn": 2}},
    ]
    write_unchecked_session_file(tmp_path / "card.gbdraw", document)
    _cli(monkeypatch, ["render", "--session", "card.gbdraw", "-o", "card", "-f", "svg"])
    card = (tmp_path / "card.svg").read_text(encoding="utf-8")
    assert 'data-record-key="record-1:2"' in card
    assert card.replace("record-1:2", "record-2").replace("record-1:1", "record-1") == (
        tmp_path / "t.svg"
    ).read_text(encoding="utf-8")


def test_a_cli_session_of_a_records_table_drawing_one_file_record_twice_replays_the_figure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """OV-2716: two rows of one file record bind that record twice, each with its own window."""
    monkeypatch.chdir(tmp_path)
    table = _table(
        tmp_path,
        [f"{MITO}\t\tM\t1\t1\t100-900", f"{MITO}\t\tT\t1\t2\t2000-3000:rc", f"{TWO_RECORDS}\tTESTA\tN\t2\t1\t"],
        "gbk\trecord_id\trecord_label\trow\tcolumn\tregion",
    )
    _cli(monkeypatch, ["linear", "--records_table", str(table), "-o", "v", "-f", "svg", "--save_session"])
    request = _saved_request(tmp_path / "v.gbdraw")
    records = request["records"]
    # Both rows read the file itself, each its own window; no in-memory GenBank copy.
    assert [record["source"]["kind"] for record in records] == ["genbank", "genbank", "genbank"]
    assert records[0]["source"] == records[1]["source"]
    assert [
        [(w["start"], w["end"], w["reverseComplement"]) for w in record["windows"]] for record in records
    ] == [[(100, 900, False)], [(2000, 3000, True)], []]
    _cli(monkeypatch, ["render", "--session", "v.gbdraw", "-o", "replay", "-f", "svg"])
    assert (tmp_path / "replay.svg").read_bytes() == (tmp_path / "v.svg").read_bytes()


@pytest.mark.parametrize("second", ["5000-5600", "5000-5600:rc"], ids=["forward", "reversed"])
def test_a_session_of_a_resolved_request_keeps_the_windows_of_a_file_record(tmp_path: Path, second: str) -> None:
    """OV-2923: a request whose records are the planner's drawn units (``resolve_request``) is saved with each file
    record and its windows, so the replay draws the "//" gap between the windows as the live figure does."""
    from gbdraw.api import (
        load_session_document,
        materialize_session,
        render_session,
        resolve_request,
        save_session_document,
    )
    from gbdraw.api.requests import RenderOutputRequest

    output = RenderOutputRequest(output_prefix="live", output_directory=tmp_path, overwrite=True)
    resolved = resolve_request(LinearDiagramRequest(
        records=(
            _windowed(RecordWindow.from_region(parse_region_spec("100-900")),
                      RecordWindow.from_region(parse_region_spec(second)), grid_row=1),
            RecordInput(GenBankInputSource(TWO_RECORDS), selector=RecordSelector(raw="TESTA", record_id="TESTA", record_index=None),
                        record_key="testa", presentation=RecordPresentation(grid_row=2)),
        ),
        options=LinearDiagramOptions(),
        layout=LinearMultiRecordOptions(),
        output=output,
    ))
    render_request(resolved)
    live = (tmp_path / "live.svg").read_text(encoding="utf-8")
    session = tmp_path / "saved.gbdraw"
    save_session_document(session, resolved)
    records = _saved_request(session)["records"]
    assert [(record["recordKey"], len(record["windows"])) for record in records] == [("mito", 2), ("testa", 0)]
    replay_dir = tmp_path / "replay"
    replay_dir.mkdir()
    with materialize_session(load_session_document(session), output_directory=replay_dir) as materialized:
        render_session(materialized)
    replay = next(replay_dir.glob("*.svg")).read_text(encoding="utf-8")
    assert live.count('data-gbdraw-row-mark="gap"') == 1
    assert replay == live


@pytest.mark.parametrize("reverse", [False, True], ids=["forward", "reversed"])
def test_the_transform_of_a_window_across_the_origin_maps_bases_on_both_sides_of_it(reverse: bool) -> None:
    """The transform of a crossing window maps a base after the origin into the window, not before its start."""
    collection = resolve_request_records(
        LinearDiagramRequest(records=(_windowed(replace(_ACROSS, reverse_complement=reverse)),))
    )
    [transform] = collection.transforms
    # Base 577 (tRNA-Phe start, after the origin) and base 16500 (before it), both inside the window.
    expected = {577: 423, 16500: 1069} if reverse else {577: 1146, 16500: 500}
    assert {base: transform.source_base_to_display_index(base) for base in expected} == expected
    assert {transform.display_index_to_source_base(index): index for index in expected.values()} == expected
    assert transform.source_boundary_to_display_offset(1000 if not reverse else 15999) == 1570


def test_anchor_feature_aligns_a_feature_after_the_origin_of_a_window_across_it() -> None:
    from gbdraw.api import FeatureSelector

    request = LinearDiagramRequest(
        records=(_windowed(_ACROSS), RecordInput(GenBankInputSource(MITO), region=parse_region_spec("1-5000"))),
        options=LinearDiagramOptions(anchor_feature=FeatureSelector(value="tRNA-Phe", key="product")),
    )
    plan = plan_linear_request(request)
    assert plan.feature_anchor is not None
    # tRNA-Phe (577..647) is base index 1146 of the window across the origin and 576 of the region 1-5000.
    assert plan.feature_anchor.statuses == ("reference", "aligned")
    assert plan.feature_anchor.starts == (1146.0, 576.0)


@pytest.mark.parametrize("reverse", [False, True], ids=["forward", "reversed"])
def test_a_wrapped_window_across_the_origin_draws_its_local_parts_where_they_lie(reverse: bool, tmp_path: Path) -> None:
    from gbdraw.layout.record_coordinates import SourceInterval

    request = LinearDiagramRequest(records=(_windowed(replace(_ACROSS, reverse_complement=reverse)),))
    [transform] = resolve_request_records(request).transforms
    wrapped = transform.with_wrap(500, 1570)
    # A local part across the source origin (local 570) stays one piece on the second row.
    [fragment] = wrapped.project_local_parts((SourceInterval(520, 700),))
    assert (fragment.row, fragment.display_start, fragment.display_end) == (1, 20, 200)
    wrapped_request = replace(
        request,
        options=LinearDiagramOptions(config_overrides={"canvas.linear.wrap_length": 500}),
        output=replace(request.output, output_directory=tmp_path),
    )
    [svg] = [path for path in render_request(wrapped_request).output_paths if path.suffix == ".svg"]
    assert len(set(re.findall(r'data-gbdraw-wrap-row="(\d+)"', svg.read_text(encoding="utf-8")))) == 4


@pytest.mark.parametrize("reverse", [False, True], ids=["forward", "reversed"])
def test_depth_on_a_window_across_the_origin_reads_both_sides_of_it(reverse: bool) -> None:
    import pandas as pd

    from gbdraw.analysis.depth import depth_df

    [unit] = resolve_request_records(
        LinearDiagramRequest(records=(_windowed(replace(_ACROSS, reverse_complement=reverse)),))
    ).records
    # Depth 1 on bases 1..1000 (after the origin), 2 on 16000..16569 (before it); bases 1001..15999 have none.
    positions = [*range(1, 1001), *range(16000, 16570)]
    table = pd.DataFrame({
        "reference_name": unit.id,
        "position": positions,
        "depth": [1 if position <= 1000 else 2 for position in positions],
    })
    depths = depth_df(unit, table, window=10, step=10, circular=True)["depth"].tolist()
    # 570 bases before the origin, then 1000 after it (read backwards when reversed). Each window is centred on
    # its sample; the first window reads 5 bases of the source before the window (no depth), and the window at
    # the origin reads 5 bases of each side.
    expected = (
        [0.5] + [1.0] * 99 + [1.5] + [2.0] * 56 if reverse else [1.0] + [2.0] * 56 + [1.5] + [1.0] * 99
    )
    assert depths == expected


@pytest.mark.parametrize("reverse", [False, True], ids=["forward", "reversed"])
def test_conservation_hits_on_a_window_across_the_origin_are_drawn_in_its_search_frame(reverse: bool) -> None:
    import pandas as pd

    from gbdraw.analysis.conservation import ConservationSource, _normalize_source_hits_for_record

    [unit] = resolve_request_records(
        LinearDiagramRequest(records=(_windowed(replace(_ACROSS, reverse_complement=reverse)),))
    ).records
    # Frame positions 1..100 (source 16000..16099) and 601..700 (source 31..130, after the origin).
    hits = pd.DataFrame([
        {"source_hit_index": index, "query": unit.id, "subject": "s", "qstart": start, "qend": start + 99,
         "sstart": 1, "send": 100, "identity": 90.0, "alignment_length": 100.0, "mismatches": 0.0,
         "gap_opens": 0.0, "evalue": 1e-30, "bitscore": 200.0}
        for index, start in enumerate((1, 601))
    ])
    rows = _normalize_source_hits_for_record(
        ConservationSource(0, "s", None, hits), track_index=0, reference_side="query", record=unit
    )
    drawn = sorted(zip(rows["draw_start"], rows["draw_end"]))
    assert drawn == ([(870.0, 970.0), (1470.0, 1570.0)] if reverse else [(0.0, 100.0), (600.0, 700.0)])


def _web_catalog_item(request: LinearDiagramRequest, tmp_path: Path) -> dict:
    """The feature catalog item the Web receives for ``request`` (`render_canonical_web_request`)."""
    from gbdraw.api import resolve_request
    from gbdraw.web_support.request_render import render_canonical_web_request

    encoded = encode_canonical_request(resolve_request(request))
    paths = {}
    for index, resource in enumerate(encoded.resources):
        path = tmp_path / f"resource-{index}"
        path.write_bytes(resource.content if resource.content is not None else Path(str(resource.source_path)).read_bytes())
        paths[resource.resource_id] = str(path)
    reply = render_canonical_web_request(encoded.payload, resource_paths=paths, output_directory=tmp_path / "out")
    [item] = reply["featureCatalog"]["items"]
    return item


@pytest.mark.parametrize(
    ("windows", "feature"),
    [
        ((RecordWindow(1000, 6000), RecordWindow(4000, 9000)), _IN_OVERLAP),
        ((RecordWindow(1000, 6500), RecordWindow(7000, 9000)), _ACROSS_GAP),
    ],
    ids=["inside-an-overlap", "across-a-left-out-interval"],
)
def test_the_web_catalog_names_a_windowed_feature_by_its_record_and_python_binds_its_edit(
    windows, feature: str, tmp_path: Path
) -> None:
    """OV-2717: catalog rows name the record (its identity key); a rendered row names its drawn unit by index.

    The Web keys a feature edit by a catalog row's record key and biological feature ID; Python binds that edit to
    the feature in every window that draws it.
    """
    item = _web_catalog_item(LinearDiagramRequest(records=(_windowed(*windows),)), tmp_path)
    assert item["recordKeys"] == ["mito~w1", "mito~w2"], "the drawn units, in drawing order"
    assert {row["recordKey"] for row in item["biologicalFeatures"]} == {"mito"}
    assert len({row["biologicalFeatureId"] for row in item["biologicalFeatures"]}) == len(item["biologicalFeatures"])
    drawn = [row for row in item["features"] if row["biologicalFeatureId"] == feature]
    assert sorted(row["recordIndex"] for row in drawn) == [0, 1], "one rendered row in each window"
    assert {row["recordKey"] for row in item["features"]} == {"mito"}
    [edited_row] = {(row["recordKey"], row["biologicalFeatureId"]) for row in drawn}
    plan = plan_linear_request(
        LinearDiagramRequest(
            records=(_windowed(*windows),),
            options=LinearDiagramOptions(feature_overrides=(FeatureOverride(*edited_row, label_text="EDITED"),)),
        )
    )
    assert plan.inputs is not None
    assert [[override.label_text for override in unit.overrides.values()] for unit in plan.inputs.record_features] == [
        ["EDITED"],
        ["EDITED"],
    ]
    assert not plan.inputs.feature_identity_notices


def test_the_web_catalog_gives_source_coordinates_for_a_window_across_the_origin(tmp_path: Path) -> None:
    """A window across the origin (D-28): every catalog location is a span of the source record, split at its origin."""
    item = _web_catalog_item(LinearDiagramRequest(records=(_windowed(RecordWindow(16000, 700)),)), tmp_path)
    rows = {(row["type"], row.get("start"), row.get("end")): row for row in item["biologicalFeatures"]}
    # tRNA-Phe, 577..647, lies after the origin.
    assert ("tRNA", 576, 647) in rows
    for row in item["biologicalFeatures"]:
        for part in row.get("location_parts") or [{"start": row["start"], "end": row["end"]}]:
            assert 0 <= part["start"] < part["end"] <= 16569, row
    # The D-loop, complement(join(16024..16569,1..576)), is drawn whole and named in both of its source parts.
    [d_loop] = [row for row in item["biologicalFeatures"] if row["type"] == "D-loop"]
    assert sorted((part["start"], part["end"]) for part in d_loop["location_parts"]) == [(0, 576), (16023, 16569)]


@pytest.mark.parametrize("cardinality", [RecordCardinality.ALL, RecordCardinality.EXACTLY_ONE])
def test_windows_on_a_file_of_several_records_without_a_selector_are_refused(cardinality: RecordCardinality) -> None:
    """A card's windows are spans of one record: a file of several records needs the record named (the Web writes
    ``exactly_one`` for a card with windows, and its producers name the record), never every record."""
    source = GenBankInputSource(REPO / "tests" / "test_inputs" / "GCF_000196095.1_ASM19609v1_genomic.gbff")
    windows = (RecordWindow(1, 100), RecordWindow(200, 300))
    with pytest.raises(ValidationError, match=r"selectorless per-input region|requires exactly one record"):
        resolve_request_records(LinearDiagramRequest(records=(RecordInput(source, cardinality=cardinality, windows=windows),)))


def test_window_refusals_of_a_record_input_name_the_record() -> None:
    # PD-OI-122: a record-wide reverse complement together with windows, and several windows in Circular mode, are
    # refused with a reason that names the record.
    two = (RecordWindow(1, 100), RecordWindow(200, 300))
    with pytest.raises(ValidationError, match=r"Record 'mito' \(HmmtDNA\.gbk\)"):
        RecordInput(GenBankInputSource(MITO), windows=two, presentation=RecordPresentation(reverse_complement=True),
                    record_key="mito")
    with pytest.raises(ValidationError, match=r"Record 'mito' \(HmmtDNA\.gbk\)"):
        RecordInput(GenBankInputSource(MITO), windows=(RecordWindow(1, 100, True),),
                    presentation=RecordPresentation(reverse_complement=True), record_key="mito")
    with pytest.raises(ValidationError, match=r"Record 'NC_012920\.1' \(HmmtDNA\.gbk\)"):
        CircularDiagramRequest(records=(RecordInput(
            GenBankInputSource(MITO), selector=RecordSelector("NC_012920.1", "NC_012920.1", None), windows=two
        ),))


def test_a_region_passed_to_replace_replaces_the_windows_as_on_dev() -> None:
    # OV-2769: dataclasses.replace(record, region=...) re-crops and region=None uncrops, as it did before windows.
    cropped = RecordInput(GenBankInputSource(MITO), region=parse_region_spec("NC_012920.1:51-350:rc"), record_key="m")
    assert cropped.windows == (RecordWindow(51, 350, True),)
    whole = replace(cropped, region=None)
    assert whole.windows == () and whole.selector == cropped.selector
    recropped = replace(cropped, region=parse_region_spec("10-20"))
    assert recropped.windows == (RecordWindow(10, 20),) and recropped.selector == cropped.selector
    # Without region, replace keeps the windows (and several windows stay several).
    two = RecordInput(GenBankInputSource(MITO), windows=(RecordWindow(1, 100), RecordWindow(200, 300)), record_key="m")
    assert replace(two, record_key="n").windows == two.windows
    assert replace(cropped, record_key="n").windows == cropped.windows
    # replace passes the record's own windows with the region: the region is the record's one span, also for a
    # record drawn in several windows.
    assert replace(two, region=parse_region_spec("5-9")).windows == (RecordWindow(5, 9),)
    # Windows the caller gives with a region conflict with it, at construction and in replace, and the refusal
    # names the record.
    with pytest.raises(ValidationError, match=r"Record 'm' \(HmmtDNA\.gbk\): give a region or windows, not both"):
        RecordInput(GenBankInputSource(MITO), region=parse_region_spec("1-5"), windows=(RecordWindow(1, 9),),
                    record_key="m")
    with pytest.raises(ValidationError, match="give a region or windows, not both"):
        replace(cropped, region=parse_region_spec("1-5"), windows=[RecordWindow(1, 9)])
    # No windows given is no conflict.
    assert RecordInput(GenBankInputSource(MITO), region=parse_region_spec("1-5"), windows=()).windows == (
        RecordWindow(1, 5),
    )
    # A record with several windows still names windows when region is read.
    with pytest.raises(ValidationError, match="several windows"):
        _ = two.region
