from __future__ import annotations

import ast
from pathlib import Path
from types import SimpleNamespace

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from gbdraw.api.options import LinearDiagramOptions
from gbdraw.api.record_planning import (
    resolve_circular_batch_outputs,
    resolve_linear_options,
    resolve_record_inputs,
)
from gbdraw.api.request_render import resolve_request
from gbdraw.api.requests import (
    CircularBatchOutputPolicy,
    CircularBatchRequest,
    CircularDiagramRequest,
    GenBankInputSource,
    InMemoryRecordSource,
    LinearDiagramRequest,
    RecordCardinality,
    RecordCollectionOptions,
    RecordDisplayOptions,
    RecordInput,
    RecordPresentation,
    RenderOutputRequest,
)
from gbdraw.circular import _circular_cli_record_cardinality
from gbdraw.exceptions import ValidationError
from gbdraw.io.record_select import parse_record_selector
from gbdraw.io.regions import parse_region_spec
from gbdraw.linear import _linear_cli_record_cardinality
from gbdraw.session import (
    build_session_document,
    materialize_session,
    session_to_request,
)
from gbdraw.session_request_codec import encode_canonical_request


def _record(record_id: str, sequence: str = "AAACCG") -> SeqRecord:
    return SeqRecord(
        Seq(sequence),
        id=record_id,
        annotations={"molecule_type": "DNA"},
    )


def _write_genbank(path: Path, *records: SeqRecord) -> None:
    SeqIO.write(list(records), path, "genbank")


def _resolve(
    inputs: tuple[RecordInput, ...],
    loader,
):
    return resolve_record_inputs(
        inputs,
        gff_candidate_features=None,
        gff_keep_all_features=False,
        genbank_loader=loader,
    )


@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("start", [None, 1, 5, 10])
def test_display_context_preserves_source_and_orientation(reverse, start):
    record = _record("circular", "AAACCGTTAA")
    record.annotations["topology"] = "circular"
    resolved = _resolve((RecordInput(
        InMemoryRecordSource(record),
        presentation=RecordPresentation(reverse_complement=reverse),
        display=RecordDisplayOptions(start_coordinate=start),
    ),), None)
    assert resolved.displays[0].start_coordinate == start
    assert resolved.transforms[0].length == 10
    assert resolved.transforms[0].source_base == (10 if reverse else 1)
    assert resolved.transforms[0].source_step == (-1 if reverse else 1)
    assert resolved.provenance[0].display == RecordDisplayOptions(start_coordinate=start)
    assert str(record.seq) == "AAACCGTTAA"
    assert "gbdraw_coord_step" not in record.annotations


def test_display_all_reports_the_invalid_record(tmp_path):
    records = [_record("duplicate", "A" * 10), _record("duplicate", "A" * 4)]
    for record in records:
        record.annotations["topology"] = "circular"
    with pytest.raises(ValidationError, match="source-row:2.*Display start"):
        _resolve((RecordInput(
            GenBankInputSource(tmp_path / "records.gb"),
            cardinality=RecordCardinality.ALL,
            record_key="source-row",
            display=RecordDisplayOptions(start_coordinate=5),
        ),), lambda _paths: records)


def test_record_cardinality_exactly_one_first_and_all(tmp_path: Path) -> None:
    source = GenBankInputSource(tmp_path / "records.gb")
    loaded = [_record("first"), _record("second")]

    with pytest.raises(ValidationError, match="requires exactly one record"):
        _resolve((RecordInput(source=source),), lambda _paths: loaded)

    first = _resolve(
        (
            RecordInput(
                source=source,
                cardinality=RecordCardinality.FIRST,
            ),
        ),
        lambda _paths: loaded,
    )
    all_records = _resolve(
        (
            RecordInput(
                source=source,
                cardinality=RecordCardinality.ALL,
            ),
        ),
        lambda _paths: loaded,
    )

    assert [record.id for record in first.records] == ["first"]
    assert [record.id for record in all_records.records] == ["first", "second"]


def test_record_source_is_loaded_once_for_multiple_selectors(tmp_path: Path) -> None:
    source = GenBankInputSource(tmp_path / "records.gb")
    calls = 0

    def loader(_paths):
        nonlocal calls
        calls += 1
        return [_record("left"), _record("right")]

    resolved = _resolve(
        (
            RecordInput(
                source=source,
                selector=parse_record_selector("#1"),
                record_key="left",
            ),
            RecordInput(
                source=source,
                selector=parse_record_selector("#2"),
                record_key="right",
            ),
        ),
        loader,
    )

    assert calls == 1
    assert [record.id for record in resolved.records] == ["left", "right"]
    assert [item.source_record_index for item in resolved.provenance] == [0, 1]


def test_selection_reverse_region_order_and_provenance(tmp_path: Path) -> None:
    source_path = tmp_path / "records.gb"
    resolved = _resolve(
        (
            RecordInput(
                source=GenBankInputSource(source_path),
                selector=parse_record_selector("chosen"),
                region=parse_region_spec("2-4"),
                presentation=RecordPresentation(
                    label="Shown",
                    subtitle="Detail",
                    reverse_complement=True,
                ),
                record_key="chosen-row",
            ),
        ),
        lambda _paths: [_record("other", "TTTT"), _record("chosen", "AAACCG")],
    )

    record = resolved.records[0]
    provenance = resolved.provenance[0]
    assert str(record.seq) == "GGT"
    assert record.annotations["gbdraw_record_label"] == "Shown"
    assert record.annotations["gbdraw_record_subtitle"] == "Detail"
    assert record.annotations["gbdraw_input_index"] == 0
    assert record.annotations["gbdraw_source_record_index"] == 1
    assert record.annotations["gbdraw_source_file"] == str(source_path)
    assert provenance.source_record_id == "chosen"
    assert provenance.source_record_count == 2
    assert provenance.record_key == "chosen-row"
    assert provenance.selector is not None
    assert provenance.region is not None


def test_region_reverse_crops_once_and_flips_boundary_crossing_feature(
    tmp_path: Path,
) -> None:
    source_record = _record("chosen", "AAACCGTT")
    source_record.features = [
        SeqFeature(FeatureLocation(1, 7, strand=1), type="CDS")
    ]

    resolved = _resolve(
        (
            RecordInput(
                source=GenBankInputSource(tmp_path / "records.gb"),
                selector=parse_record_selector("chosen"),
                region=parse_region_spec("3-6:rc"),
                record_key="chosen-crop",
            ),
        ),
        lambda _paths: [source_record],
    )

    record = resolved.records[0]
    assert str(record.seq) == "CGGT"
    assert len(record.features) == 1
    assert int(record.features[0].location.start) == 0
    assert int(record.features[0].location.end) == 4
    assert record.features[0].location.strand == -1
    assert str(source_record.seq) == "AAACCGTT"
    assert source_record.features[0].location.strand == 1


def test_batch_output_policy_disambiguates_duplicate_record_ids() -> None:
    outputs = resolve_circular_batch_outputs(
        CircularBatchOutputPolicy(),
        (
            _record("duplicate"),
            _record("duplicate_2"),
            _record("duplicate"),
        ),
    )

    assert [output.output_prefix for output in outputs] == [
        "duplicate",
        "duplicate_2",
        "duplicate_3",
    ]


@pytest.mark.parametrize(
    "record_id",
    (
        "nested/record",
        "../record",
        "/absolute/record",
        r"nested\record",
        "record\x00hidden",
        "record\nhidden",
        "NUL",
        "record:stream",
    ),
)
def test_implicit_circular_output_rejects_path_like_record_ids(
    record_id: str,
) -> None:
    with pytest.raises(
        ValidationError,
        match="cannot be used as an implicit output filename prefix",
    ):
        resolve_circular_batch_outputs(
            CircularBatchOutputPolicy(),
            (_record(record_id),),
        )


def test_implicit_circular_grid_output_rejects_path_like_record_id() -> None:
    request = CircularDiagramRequest(
        records=(
            RecordInput(
                source=InMemoryRecordSource(_record("../outside")),
            ),
        ),
        output=RenderOutputRequest(resolve_prefix_from_first_record=True),
    )

    with pytest.raises(
        ValidationError,
        match="cannot be used as an implicit output filename prefix",
    ):
        resolve_request(request)


def test_batch_request_expands_all_records_before_resolving_outputs(
    tmp_path: Path,
) -> None:
    source_path = tmp_path / "duplicates.gb"
    _write_genbank(source_path, _record("duplicate"), _record("duplicate"))
    unresolved = CircularBatchRequest(
        records=(
            RecordInput(
                source=GenBankInputSource(source_path),
                cardinality=RecordCardinality.ALL,
            ),
        ),
        output_policy=CircularBatchOutputPolicy(),
    )

    resolved = resolve_request(unresolved)

    assert isinstance(resolved, CircularBatchRequest)
    assert len(resolved.records) == 2
    assert all(
        record.cardinality is RecordCardinality.EXACTLY_ONE
        for record in resolved.records
    )
    assert [output.output_prefix for output in resolved.outputs] == [
        "duplicate",
        "duplicate_2",
    ]
    assert resolved.output_policy is None


def test_linear_comparison_table_resolves_after_record_expansion(
    tmp_path: Path,
) -> None:
    source_path = tmp_path / "records.gb"
    _write_genbank(source_path, _record("left"), _record("right"))
    blast_path = tmp_path / "pair.tsv"
    blast_path.write_text(
        "left\tright\t99\t6\t0\t0\t1\t6\t1\t6\t1e-20\t100\n",
        encoding="utf-8",
    )
    table_path = tmp_path / "comparisons.tsv"
    table_path.write_text(
        "blast\tquery\tsubject\n"
        f"{blast_path.name}\t#1\t#2\n",
        encoding="utf-8",
    )
    unresolved = LinearDiagramRequest(
        records=(
            RecordInput(
                source=GenBankInputSource(source_path),
                cardinality=RecordCardinality.ALL,
            ),
        ),
        options=LinearDiagramOptions(comparison_table_file=str(table_path)),
    )

    resolved = resolve_request(unresolved)

    assert isinstance(resolved, LinearDiagramRequest)
    assert resolved.options.comparison_table_file is None
    assert resolved.options.linear_comparisons is not None
    comparison = resolved.options.linear_comparisons[0]
    assert comparison.query_record_index == 0
    assert comparison.subject_record_index == 1
    assert len(comparison.matches) == 1


def test_linear_comparison_reader_does_not_hide_unexpected_errors(
    monkeypatch,
) -> None:
    table = SimpleNamespace(
        table_path="comparisons.tsv",
        rows=(
            SimpleNamespace(
                query="#1",
                subject="#2",
                blast="pair.tsv",
                row_number=2,
            ),
        ),
    )
    monkeypatch.setattr(
        "gbdraw.api.record_planning.read_comparisons_table",
        lambda _path: table,
    )

    def fail_reader(*_args, **_kwargs):
        raise RuntimeError("reader implementation bug")

    monkeypatch.setattr(
        "gbdraw.api.record_planning.pd.read_csv",
        fail_reader,
    )

    with pytest.raises(RuntimeError, match="reader implementation bug"):
        resolve_linear_options(
            LinearDiagramOptions(comparison_table_file="comparisons.tsv"),
            records=(_record("left"), _record("right")),
            layout=None,
        )


def test_schema6_round_trips_unresolved_then_materializes_session(
    tmp_path: Path,
) -> None:
    source_path = tmp_path / "records.gb"
    _write_genbank(source_path, _record("one"), _record("two"))
    unresolved = CircularDiagramRequest(
        records=(
            RecordInput(
                source=GenBankInputSource(source_path),
                cardinality=RecordCardinality.ALL,
            ),
        ),
    )

    unresolved_encoded = encode_canonical_request(unresolved)
    assert unresolved_encoded.payload["schema"] == 6
    assert unresolved_encoded.payload["records"][0]["cardinality"] == "all"

    resolved = resolve_request(unresolved)
    encoded = encode_canonical_request(resolved)
    assert len(resolved.records) == 2
    assert all(
        record.cardinality is RecordCardinality.EXACTLY_ONE
        for record in resolved.records
    )
    assert all(row["cardinality"] == "exactly_one" for row in encoded.payload["records"])

    document = build_session_document(unresolved)
    with materialize_session(
        document,
        output_directory=tmp_path,
        temporary_directory=tmp_path / "materialized",
    ) as materialized:
        decoded = session_to_request(materialized)
        assert len(decoded.records) == 2


@pytest.mark.parametrize(
    ("is_gff", "source_count", "comparisons", "expected"),
    (
        (False, 1, False, RecordCardinality.ALL),
        (False, 1, True, RecordCardinality.ALL),
        (False, 2, True, RecordCardinality.FIRST),
        (True, 1, False, RecordCardinality.ALL),
        (True, 1, True, RecordCardinality.FIRST),
    ),
)
def test_cli_legacy_cardinality_is_explicit(
    is_gff: bool,
    source_count: int,
    comparisons: bool,
    expected: RecordCardinality,
) -> None:
    assert _linear_cli_record_cardinality(
        is_gff_source=is_gff,
        source_count=source_count,
        load_comparison=comparisons,
    ) is expected
    assert _circular_cli_record_cardinality() is RecordCardinality.ALL


def test_record_collection_labels_are_strict_strings() -> None:
    with pytest.raises(ValidationError, match="labels must contain strings"):
        RecordCollectionOptions(labels=(None,))  # type: ignore[arg-type]


def test_cli_adapters_do_not_import_or_call_domain_table_readers() -> None:
    forbidden = {
        "apply_region_specs",
        "load_default_colors",
        "load_gbks",
        "load_gff_fasta",
        "read_annotation_table",
        "read_color_table",
        "read_comparisons_table",
        "read_feature_visibility_file",
        "read_filter_list_file",
        "read_label_override_file",
        "read_qualifier_priority_file",
        "read_records_table",
    }
    repo_root = Path(__file__).parents[1]
    for relative_path in ("gbdraw/circular.py", "gbdraw/linear.py"):
        tree = ast.parse((repo_root / relative_path).read_text(encoding="utf-8"))
        referenced = {
            node.id
            for node in ast.walk(tree)
            if isinstance(node, ast.Name)
        }
        imported = {
            alias.name
            for node in ast.walk(tree)
            if isinstance(node, (ast.Import, ast.ImportFrom))
            for alias in node.names
        }
        assert not (forbidden & (referenced | imported))


@pytest.mark.parametrize("topology,override,accepted", [
    ("circular", None, True), ("linear", None, False), (None, None, False),
    ("other", True, True), ("LINEAR", True, True), (" Circular ", None, True),
])
def test_display_resolver_observes_source_topology(topology, override, accepted):
    record = _record("topology", "A" * 10)
    if topology is not None:
        record.annotations["topology"] = topology
    inputs = (RecordInput(InMemoryRecordSource(record),
                          display=RecordDisplayOptions(override, 5)),)
    if not accepted:
        with pytest.raises(ValidationError, match="record-1.*circular record"):
            _resolve(inputs, None)
    else:
        result = _resolve(inputs, None)
        assert result.displays[0].is_circular
        assert result.transforms[0].length == 10


@pytest.mark.parametrize("reverse", [False, True])
def test_display_materialization_and_replanning_keeps_rotation(reverse):
    from gbdraw.api.request_render import plan_request
    record = _record("circular", "AAACCGTTAA")
    record.annotations["topology"] = "circular"
    request = LinearDiagramRequest(records=(RecordInput(
        InMemoryRecordSource(record),
        presentation=RecordPresentation(reverse_complement=reverse),
        display=RecordDisplayOptions(None, 1),
    ),))
    first = plan_request(request)
    second = plan_request(first.request)
    assert first.transforms == second.transforms
    assert first.displays == second.displays
    assert first.records[0].seq == second.records[0].seq
    assert second.request.records[0].display == RecordDisplayOptions(None, 1)
    assert not second.request.records[0].presentation.reverse_complement


def test_display_collection_crop_targets_only_the_selected_record():
    records = [_record("left", "A" * 10), _record("right", "A" * 10)]
    for record in records:
        record.annotations["topology"] = "circular"
    inputs = tuple(RecordInput(InMemoryRecordSource(record),
                   display=RecordDisplayOptions(None, 5 if i == 0 else None))
                   for i, record in enumerate(records))
    result = resolve_record_inputs(
        inputs, record_options=RecordCollectionOptions(regions=(parse_region_spec("#2:2-6:rc"),)),
        gff_candidate_features=None, gff_keep_all_features=False,
    )
    assert [len(record) for record in result.records] == [10, 5]
    assert [item.has_collection_region for item in result.provenance] == [False, True]
    assert [t.length for t in result.transforms] == [10, 10]
    assert result.transforms[1].source_base == 6
    assert result.transforms[1].source_step == -1
    invalid = (inputs[0], RecordInput(inputs[1].source, display=RecordDisplayOptions(None, 5)))
    with pytest.raises(ValidationError, match="record-2.*crop"):
        resolve_record_inputs(invalid,
            record_options=RecordCollectionOptions(regions=(parse_region_spec("#2:2-6"),)),
            gff_candidate_features=None, gff_keep_all_features=False)


def test_display_input_crop_and_known_materialized_crop_validation():
    from gbdraw.api.request_render import plan_request
    record = _record("circular", "A" * 10)
    record.annotations["topology"] = "circular"
    source = InMemoryRecordSource(record)
    with pytest.raises(ValidationError, match="crop"):
        _resolve((RecordInput(source, region=parse_region_spec("2-6"),
                              display=RecordDisplayOptions(None, 5)),), None)
    plan = plan_request(LinearDiagramRequest(records=(RecordInput(source, region=parse_region_spec("2-6:rc")),)))
    again = plan_request(plan.request)
    assert again.transforms == plan.transforms
    assert again.transforms[0].length == 10
    assert again.transforms[0].source_base == 6
    assert len(again.records[0]) == 5
    with pytest.raises(ValidationError, match="crop"):
        _resolve((RecordInput(InMemoryRecordSource(plan.records[0]),
                              display=RecordDisplayOptions(None, 5)),), None)


def test_display_batch_and_duplicate_instances_keep_aligned_context(tmp_path):
    from gbdraw.api.request_render import plan_request
    record = _record("duplicate", "A" * 10)
    record.annotations["topology"] = "circular"
    path = tmp_path / "same.gb"
    _write_genbank(path, record, record)
    source = GenBankInputSource(path)
    inputs = (RecordInput(source, selector=parse_record_selector("#1"), record_key="one",
                          display=RecordDisplayOptions(None, 1)),
              RecordInput(source, selector=parse_record_selector("#2"), record_key="two",
                          display=RecordDisplayOptions(None, 5)),
              RecordInput(source, selector=parse_record_selector("#1"), record_key="again",
                          display=RecordDisplayOptions(None, 10)))
    batch = plan_request(CircularBatchRequest(records=inputs, output_policy=CircularBatchOutputPolicy()))
    items = batch.item_plans()
    assert [i.request.records[0].display.start_coordinate for i in items] == [1, 5, 10]
    assert [i.transforms[0] for i in items] == list(batch.transforms)
    assert [i.displays[0] for i in items] == list(batch.displays)
    assert [i.provenance[0].record_key for i in items] == ["one", "two", "again"]
    for i in items:
        assert len(i.records) == len(i.provenance) == len(i.displays) == len(i.transforms) == 1


def test_display_zero_length_rejected_with_context():
    with pytest.raises(ValidationError, match="record-1.*Complete source length"):
        _resolve((RecordInput(InMemoryRecordSource(_record("empty", ""))),), None)


@pytest.mark.parametrize("crop_rc", [False, True])
@pytest.mark.parametrize("reverse", [False, True])
@pytest.mark.parametrize("recrop", [False, True])
def test_display_unset_preserves_external_crop_source_length(crop_rc, reverse, recrop):
    from copy import deepcopy

    from gbdraw.api.request_render import plan_request
    from gbdraw.io.regions import apply_region_specs

    source = _record("external-crop", "AAACCGTTAA")
    source.annotations["topology"] = "circular"
    original_annotations = deepcopy(source.annotations)
    cropped = apply_region_specs(
        [source], [parse_region_spec("3-7:rc" if crop_rc else "3-7")]
    )[0]
    assert source.annotations == original_annotations
    assert cropped.annotations["gbdraw_source_length"] == 10
    cropped_annotations = deepcopy(cropped.annotations)
    expected_seq = source.seq[2:7]
    if crop_rc:
        expected_seq = expected_seq.reverse_complement()
    base, step = (7, -1) if crop_rc else (3, 1)
    if reverse:
        expected_seq = expected_seq.reverse_complement()
        base, step = base + step * 4, -step
    if recrop:
        recropped = apply_region_specs([cropped], [parse_region_spec("2-4:rc")])[0]
        assert len(recropped) == 3
        assert recropped.annotations["gbdraw_source_length"] == 10
        expected_seq = expected_seq[1:4]
        base += step
        if not reverse:
            expected_seq = expected_seq.reverse_complement()
            base, step = base + step * 2, -step

    plan = plan_request(LinearDiagramRequest(records=(RecordInput(
        InMemoryRecordSource(cropped),
        presentation=RecordPresentation(reverse_complement=reverse),
        region=parse_region_spec("2-4" if reverse else "2-4:rc") if recrop else None,
    ),)))
    again = plan_request(plan.request)
    assert cropped.annotations == cropped_annotations
    assert plan.records[0] is not cropped
    for resolved in (plan, again):
        assert resolved.records[0].seq == expected_seq
        assert resolved.records[0].annotations["gbdraw_source_length"] == 10
        assert resolved.provenance[0].source_length == 10
        assert resolved.transforms[0].length == 10
        assert resolved.transforms[0].source_base == base
        assert resolved.transforms[0].source_step == step
        assert resolved.displays[0].start_coordinate is None
