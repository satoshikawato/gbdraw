from __future__ import annotations

import ast
import dataclasses
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
from gbdraw.api.request_render import build_request_plan_diagram, plan_request, resolve_request
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
from gbdraw.exceptions import GbdrawError, ValidationError
from gbdraw.io.record_select import parse_record_selector
from gbdraw.io.regions import parse_region_spec
from gbdraw.linear import _linear_cli_record_cardinality, linear_main
from gbdraw.session import (
    build_session_document,
    materialize_session,
    session_to_request,
)
from gbdraw.session_request_codec import CANONICAL_REQUEST_SCHEMA, encode_canonical_request


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
    # Cropped to source 2..4 (AAC), then reverse-complemented (OV-1279).
    assert str(record.seq) == "GTT"
    assert provenance.presentation.label == "Shown"
    assert provenance.presentation.subtitle == "Detail"
    assert "gbdraw_record_label" not in record.annotations
    assert "gbdraw_record_subtitle" not in record.annotations
    assert record.annotations["gbdraw_input_index"] == 0
    assert record.annotations["gbdraw_source_record_index"] == 1
    assert record.annotations["gbdraw_source_file"] == str(source_path)
    assert provenance.source_record_id == "chosen"
    assert provenance.source_record_count == 2
    assert provenance.record_key == "chosen-row"
    assert provenance.selector is not None
    assert provenance.region is not None


def test_records_resolved_from_one_source_write_only_their_own_fields(
    tmp_path: Path,
) -> None:
    source_record = _record("chosen", "AAACCGTT")
    source_record.features = [
        SeqFeature(FeatureLocation(1, 7, strand=1), type="CDS", qualifiers={"gene": ["g"]})
    ]
    source = GenBankInputSource(tmp_path / "records.gb")

    resolved = _resolve(
        tuple(
            RecordInput(
                source=source,
                presentation=RecordPresentation(label=label, reverse_complement=reverse),
                record_key=label,
            )
            for label, reverse in (("first", False), ("second", False), ("third", True))
        ),
        lambda _paths: [source_record],
    )

    first, second, third = resolved.records
    assert [item.presentation.label for item in resolved.provenance] == [
        "first",
        "second",
        "third",
    ]
    assert all("gbdraw_record_label" not in record.annotations for record in resolved.records)
    assert "gbdraw_record_label" not in source_record.annotations
    assert first.features is not source_record.features
    # The planned records share the parsed features; nothing writes them.
    assert first.features[0] is source_record.features[0]
    assert third.features[0].location.strand == -1
    assert source_record.features[0].location.strand == 1
    assert source_record.features[0].qualifiers == {"gene": ["g"]}


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
    tmp_path: Path,
) -> None:
    blast = tmp_path / "pair.tsv"
    blast.write_text("left\tright\t90\t100\t0\t0\t1\t100\t1\t100\t1e-20\t200\n", encoding="utf-8")
    table = SimpleNamespace(
        table_path="comparisons.tsv",
        rows=(
            SimpleNamespace(
                query="#1",
                subject="#2",
                blast=str(blast),
                row_number=2,
                source="table",
            ),
        ),
    )
    monkeypatch.setattr(
        "gbdraw.api.record_planning.read_comparisons_table",
        lambda _path: table,
    )

    def fail_reader(*_args, **_kwargs):
        raise RuntimeError("reader implementation bug")

    # The shared comparison reader converts only input errors to ValidationError.
    monkeypatch.setattr(
        "gbdraw.io.comparisons.pd.read_csv",
        fail_reader,
    )

    with pytest.raises(RuntimeError, match="reader implementation bug"):
        resolve_linear_options(
            LinearDiagramOptions(comparison_table_file="comparisons.tsv"),
            records=(_record("left"), _record("right")),
            layout=None,
        )


def test_current_schema_round_trips_unresolved_then_materializes_session(
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
    assert unresolved_encoded.payload["schema"] == CANONICAL_REQUEST_SCHEMA
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


@pytest.mark.linear
def test_cli_rejects_more_blast_files_than_adjacent_record_pairs(tmp_path: Path) -> None:
    # With -b and two files, the CLI loads the first record of each file, so
    # two records have one adjacent pair and the second table has none.
    fixtures = Path(__file__).parent / "fixtures"
    tables = []
    for name, query, subject in (("a", "R2c", "LOCTEST"), ("b", "LOCTEST", "R3c")):
        table = tmp_path / f"{name}.tsv"
        table.write_text(f"{query}\t{subject}\t90.0\t100\t0\t0\t1\t100\t1\t100\t1e-20\t180\n", encoding="utf-8")
        tables.append(str(table))
    prefix = tmp_path / "out"
    session = tmp_path / "out.gbdraw-session.json"

    with pytest.raises(ValidationError, match=r"Too many -b/--blast files \(expected at most 1\)") as raised:
        linear_main(
            [
                "--gbk",
                str(fixtures / "web_comparison_shared_block.gb"),
                str(fixtures / "feature_location_search.gb"),
                "-b",
                *tables,
                "-o",
                str(prefix),
                "-f",
                "svg",
                "--session_output",
                str(session),
            ]
        )

    assert raised.value.diagnostic == {"code": "COMPARISON_INPUT"}
    assert "2 loaded record(s)" in str(raised.value)
    assert sorted(path.name for path in tmp_path.iterdir()) == ["a.tsv", "b.tsv"]


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


@pytest.mark.parametrize("query_reverse,subject_reverse", [(True, False), (False, True), (True, True)])
@pytest.mark.parametrize("cropped", [False, True])
def test_source_bound_comparison_direction_projection_round_trip(tmp_path, query_reverse, subject_reverse, cropped):
    """Saved evidence reprojects without analysis and preserves its source binding.

    The planned request keeps the rows it was given, so a saved Session stays in
    the search frame; the build projects them into the drawn views (OV-399).
    """
    from dataclasses import replace
    import pandas as pd
    from gbdraw.api.record_planning import comparison_record_bindings, project_source_bound_results
    from gbdraw.api.request_render import plan_linear_request
    from gbdraw.features.ids import compute_feature_hash
    from gbdraw.linear_comparison import LinearComparison

    inputs = []
    for key, strand in (("query", 1), ("subject", -1)):
        record = _record(key, "ACGT" * 50)
        record.features = [SeqFeature(FeatureLocation(30, 60, strand=strand), type="CDS")]
        path = tmp_path / (key + ".gbk")
        _write_genbank(path, record)
        inputs.append(RecordInput(GenBankInputSource(path), record_key=key,
            region=parse_region_spec(f"{key}:11-180") if cropped else None))
    baseline = LinearDiagramRequest(records=tuple(inputs))
    before = plan_linear_request(baseline)
    row = {"qstart": 21 if cropped else 31, "qend": 50 if cropped else 60,
           "sstart": 50 if cropped else 60, "send": 21 if cropped else 31,
           "collinearity_orientation": "minus", "bitscore": 123.5}
    for role, record, provenance in zip(("query", "subject"), before.records, before.provenance, strict=True):
        source = provenance.source_feature_catalog[0]
        row.update({f"{role}_feature_index": str(source.source_feature_index),
                    f"{role}_feature_svg_id": source.stable_feature_id,
                    f"{role}_view_feature_svg_id": compute_feature_hash(record.features[0], record_id=record.id)})
    evidence = pd.DataFrame([row])
    options = LinearDiagramOptions(linear_comparisons=(LinearComparison(0, 1, evidence),))
    baseline = replace(baseline, options=options)
    assert plan_linear_request(baseline).request.options.linear_comparisons[0].matches is evidence
    reversed_request = replace(baseline, records=tuple(replace(item,
        region=replace(item.region, reverse_complement=reverse) if item.region else None,
        presentation=RecordPresentation(reverse_complement=reverse if not item.region else False))
        for item, reverse in zip(inputs, (query_reverse, subject_reverse), strict=True)))
    after = plan_linear_request(reversed_request)
    assert after.request.options.linear_comparisons[0].matches is evidence
    projected = project_source_bound_results(
        after.request.options, comparison_record_bindings(after.records, after.provenance, search_frame=False)
    ).options.linear_comparisons[0].matches
    assert evidence.to_dict("records") == [row]
    for role, prefix, record, reverse in zip(("query", "subject"), ("q", "s"), after.records,
                                            (query_reverse, subject_reverse), strict=True):
        assert projected.iloc[0][f"{role}_feature_svg_id"] == row[f"{role}_feature_svg_id"]
        assert projected.iloc[0][f"{role}_view_feature_svg_id"] == compute_feature_hash(record.features[0], record_id=record.id)
        for endpoint in ("start", "end"):
            key = prefix + endpoint
            assert projected.iloc[0][key] == (len(record) + 1 - row[key] if reverse else row[key])
    assert projected.iloc[0].bitscore == row["bitscore"]
    assert projected.iloc[0].collinearity_orientation == ("plus" if query_reverse != subject_reverse else "minus")
    restored = project_source_bound_results(
        replace(options, linear_comparisons=(LinearComparison(0, 1, projected),)),
        comparison_record_bindings(after.records, after.provenance, search_frame=True),
    ).options
    pd.testing.assert_frame_equal(restored.linear_comparisons[0].matches, evidence)
    # Planning checks the bindings, so a Session writer rejects them (review 2).
    corrupted = evidence.copy(deep=True)
    corrupted.at[0, "query_view_feature_svg_id"] = "fnot-source-bound"
    with pytest.raises(ValidationError, match="source/crop binding"):
        plan_linear_request(replace(reversed_request, options=replace(options,
            linear_comparisons=(LinearComparison(0, 1, corrupted),))))
    corrupted.at[0, "query_feature_svg_id"] = "fwrong-source"
    corrupted_request = replace(reversed_request, options=replace(options,
        linear_comparisons=(LinearComparison(0, 1, corrupted),)))
    with pytest.raises(ValidationError, match="source feature ID"):
        plan_linear_request(corrupted_request)
    with pytest.raises(GbdrawError, match="source feature ID"):
        build_session_document(corrupted_request)


_PLANNING_EXAMPLE = Path(__file__).resolve().parents[1] / "examples" / "MellatMJNV.gb"


def test_a_drawn_diagram_ignores_edits_to_the_callers_record() -> None:
    import gbdraw

    record = next(SeqIO.parse(_PLANNING_EXAMPLE, "genbank"))
    diagram = gbdraw.draw_circular(record)
    for feature in record.features:
        if feature.type == "CDS":
            feature.qualifiers["product"] = ["EDITED AFTER DRAWING"]

    assert "EDITED AFTER DRAWING" not in diagram.to_svg(interactive=True)


def test_a_plan_ignores_edits_to_the_callers_record(tmp_path: Path) -> None:
    from Bio.SeqFeature import SimpleLocation

    from gbdraw.api import build_request_plan_diagram, plan_request

    record = next(SeqIO.parse(_PLANNING_EXAMPLE, "genbank"))
    plan = plan_request(
        CircularDiagramRequest(
            records=(RecordInput(source=InMemoryRecordSource(record)),),
            output=RenderOutputRequest(output_directory=tmp_path, overwrite=True),
        )
    )
    before = build_request_plan_diagram(plan).drawing.tostring()
    feature = next(feature for feature in record.features if feature.type == "CDS")
    feature.location = SimpleLocation(
        feature.location.start, feature.location.end, strand=-feature.location.strand
    )

    assert build_request_plan_diagram(plan).drawing.tostring() == before


def _record_snapshot(record: SeqRecord) -> tuple[object, ...]:
    def feature_snapshot(feature: SeqFeature) -> tuple[object, ...]:
        return (
            feature.type,
            feature.id,
            repr(feature.location),
            repr(feature.qualifiers),
            tuple(feature_snapshot(sub) for sub in getattr(feature, "sub_features", ())),
        )

    return (
        record.id,
        record.name,
        record.description,
        str(record.seq),
        repr(record.annotations),
        repr(record.letter_annotations),
        tuple(record.dbxrefs),
        tuple(feature_snapshot(feature) for feature in record.features),
    )


def test_rendering_leaves_the_records_parsed_from_a_source_unchanged(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import pandas as pd

    from gbdraw.api import (
        CircularDiagramOptions,
        build_prepared_interactive_context,
        build_request_plan_diagram,
        plan_request,
    )
    from gbdraw.api import request_render
    from gbdraw.features.overrides import FeatureOverride

    parsed: list[tuple[SeqRecord, tuple[object, ...]]] = []
    load_gbks = request_render.load_gbks

    def recording_loader(paths: list[str]) -> list[SeqRecord]:
        records = load_gbks(paths)
        parsed.extend((record, _record_snapshot(record)) for record in records)
        return records

    monkeypatch.setattr(request_render, "load_gbks", recording_loader)
    source = GenBankInputSource(_PLANNING_EXAMPLE)
    visibility = pd.DataFrame(
        [["*", "CDS", "product", "RCBTB1", "off"]],
        columns=["record_id", "feature_type", "qualifier", "value", "action"],
    )
    inputs = (
        RecordInput(source=source, record_key="plain"),
        RecordInput(
            source=source,
            presentation=RecordPresentation(reverse_complement=True),
            record_key="reversed",
        ),
        RecordInput(source=source, region=parse_region_spec("20000-90000:rc"), record_key="cropped"),
    )
    catalog = plan_request(
        LinearDiagramRequest(records=inputs[:1])
    ).provenance[0].source_feature_catalog
    first_cds = next(entry.biological_feature_id for entry in catalog if entry.feature_type == "CDS")
    parsed.clear()
    overrides = tuple(
        FeatureOverride(key, first_cds, label_visibility="on", label_text="OVERRIDDEN")
        for key in ("plain", "reversed", "cropped")
    )
    requests = (
        LinearDiagramRequest(
            records=inputs,
            options=LinearDiagramOptions(
                feature_visibility_table=visibility, feature_overrides=overrides
            ),
            output=RenderOutputRequest(output_directory=tmp_path, overwrite=True),
        ),
        CircularDiagramRequest(
            records=(inputs[1],),
            options=CircularDiagramOptions(
                feature_visibility_table=visibility, feature_overrides=overrides[1:2]
            ),
            output=RenderOutputRequest(output_directory=tmp_path, overwrite=True),
        ),
    )
    for request in requests:
        prepared = build_request_plan_diagram(plan_request(request))
        assert "OVERRIDDEN" in prepared.drawing.tostring()
        build_prepared_interactive_context(prepared)

    # One parse per request serves all of its inputs.
    assert len(parsed) == 2
    for record, snapshot in parsed:
        assert _record_snapshot(record) == snapshot


MELLAT = Path(__file__).resolve().parents[1] / "examples" / "MellatMJNV.gb"


def _flip_first_cds(record: SeqRecord) -> None:
    feature = next(item for item in record.features if item.type == "CDS")
    feature.location = FeatureLocation(
        feature.location.start, feature.location.end, strand=-feature.location.strand
    )


def _plan_svg(plan: object) -> str:
    return build_request_plan_diagram(plan).drawing.tostring()


def _circular_request(tmp_path: Path) -> CircularDiagramRequest:
    return CircularDiagramRequest(
        records=(RecordInput(source=GenBankInputSource(MELLAT)),),
        output=RenderOutputRequest(output_directory=tmp_path, overwrite=True),
    )


def test_editing_the_record_a_resolved_request_holds_does_not_change_its_plan(tmp_path: Path) -> None:
    resolved = resolve_request(_circular_request(tmp_path))
    plan = plan_request(resolved)
    expected = _plan_svg(plan)
    _flip_first_cds(resolved.records[0].source.record)
    assert _plan_svg(plan) == expected


def test_editing_the_record_swapped_into_a_resolved_request_does_not_change_its_plan(tmp_path: Path) -> None:
    resolved = resolve_request(_circular_request(tmp_path))
    mine = next(SeqIO.parse(MELLAT, "genbank"))
    held = resolved.records[0]
    swapped = dataclasses.replace(
        resolved,
        records=(dataclasses.replace(held, source=dataclasses.replace(held.source, record=mine)),),
    )
    plan = plan_request(swapped)
    expected = _plan_svg(plan)
    _flip_first_cds(mine)
    assert _plan_svg(plan) == expected


def _frame_test_record(length: int = 200) -> SeqRecord:
    from Bio.SeqFeature import CompoundLocation

    record = _record("frame", "ACGT" * (length // 4))
    record.features = [
        SeqFeature(FeatureLocation(30, 60, strand=1), type="CDS"),
        # A joined CDS and one across the origin of a circular record (R1).
        SeqFeature(CompoundLocation([FeatureLocation(10, 20, strand=-1), FeatureLocation(80, 95, strand=-1)]),
                   type="CDS"),
        SeqFeature(CompoundLocation([FeatureLocation(190, 200, strand=1), FeatureLocation(0, 6, strand=1)]),
                   type="CDS"),
    ]
    return record


def _member(index: int, span: tuple[int, int, int | None], protein_id: str = "p"):
    from gbdraw.analysis.protein_colinearity import OrthogroupMember

    return OrthogroupMember(
        orthogroup_id="og1", protein_id=f"{protein_id}{index}", record_index=0, feature_index=index,
        record_id="frame", label=f"cds{index}", start=span[0], end=span[1], strand=span[2],
        feature_svg_id=None, source_protein_id=None,
    )


def _anchor(binding, index: int | None, view: str, start: int, end: int, strand: int | None):
    from gbdraw.analysis.collinearity import CollinearityAnchor

    return CollinearityAnchor(
        query_protein_id="q", subject_protein_id="s", query_record_index=0, subject_record_index=1,
        query_order=1, subject_order=1, query_start=start, query_end=end, subject_start=5, subject_end=10,
        identity=99.0, evalue=0.0, bitscore=50.0, alignment_length=10, query_feature_svg_id="fq",
        subject_feature_svg_id="fs", source="protein_hit", query_unit_id="q", subject_unit_id="s",
        query_unit_kind="cds", subject_unit_kind="cds", query_locus_id=None, subject_locus_id=None,
        query_display_name="q", subject_display_name="s", query_strand=strand, subject_strand=1,
        query_view_feature_svg_id=view, query_feature_index=index,
    )


def test_typed_results_move_between_frames_by_their_tags() -> None:
    # OV-1145: an orthogroup member is tagged by its span, an anchor side by
    # its view ID; an item in the other frame moves, one in the target frame
    # stays, and the round trip is exact (a joined and an origin-spanning CDS
    # included). An untagged side follows its record (Session 39).
    from gbdraw.analysis.collinearity import CollinearityBlock, CollinearityResult
    from gbdraw.analysis.protein_colinearity import OrthogroupResult
    from gbdraw.api.record_planning import (
        ComparisonRecordBinding,
        comparison_record_binding,
        project_source_bound_results,
    )

    record = _frame_test_record()
    from gbdraw.features.source import capture_source_features

    binding = comparison_record_binding(record, capture_source_features(record), False)
    assert binding.opposite_spans[1] == (200 - 95, 200 - 10, 1)
    assert binding.current_spans[2] == binding.opposite_spans[2][:2] + (1,) == (0, 200, 1)
    other = ComparisonRecordBinding({}, {}, {}, 50, {}, {}, False)
    reversed_binding = ComparisonRecordBinding(
        binding.source, binding.opposite, binding.current, binding.length,
        binding.opposite_spans, binding.current_spans, True,
    )
    members = [_member(index, binding.opposite_spans[index]) for index in (0, 1, 2)]
    members.append(_member(3, (1, 2, 1)))
    groups = OrthogroupResult(orthogroups={"og1": members},
                              member_by_protein_id={member.protein_id: member for member in members})
    tagged = _anchor(binding, 0, binding.opposite[0], 170, 141, -1)
    untagged = _anchor(binding, None, "", 3, 20, None)
    blocks = (
        CollinearityBlock("b1", 0, 1, "plus", 1.0, (tagged,)),
        # Its sides disagree: one anchor moves, the other does not.
        CollinearityBlock("b2", 0, 1, "plus", 1.0, (tagged, _anchor(binding, 0, binding.current[0], 31, 60, 1))),
    )
    collinearity = CollinearityResult(blocks=blocks, unblocked_anchors=(untagged,))
    options = LinearDiagramOptions(orthogroups=groups, collinearity_blocks=collinearity)

    projection = project_source_bound_results(options, (binding, other))
    moved = projection.options
    assert [(member.start, member.end, member.strand) for member in moved.orthogroups.orthogroups["og1"]] == [
        binding.current_spans[0], binding.current_spans[1], binding.current_spans[2], (1, 2, 1),
    ]
    assert moved.orthogroups.member_by_protein_id["p0"] is moved.orthogroups.orthogroups["og1"][0]
    (block, kept_block) = moved.collinearity_blocks.blocks
    (anchor,) = block.anchors
    assert (anchor.query_start, anchor.query_end, anchor.query_strand) == (31, 60, 1)
    assert anchor.query_view_feature_svg_id == binding.current[0]
    assert block.orientation == "minus"  # exactly one side moved
    assert kept_block is blocks[1]
    assert moved.collinearity_blocks.unblocked_anchors[0] is untagged  # its record is not reversed
    assert [item.kind for item in projection.unmatched] == ["member", "block"]
    assert projection.unmatched[0].protein_id == "p3"

    back = project_source_bound_results(moved, (reversed_binding, other)).options
    assert repr(back.orthogroups) == repr(groups)
    assert repr(back.collinearity_blocks.blocks[0]) == repr(blocks[0])
    assert project_source_bound_results(moved, (binding, other)).options is moved

    # An untagged side of a reversed record maps x -> L + 1 - x, and a span of
    # mixed strands (strand None) stays descending after the map.
    reversed_untagged = project_source_bound_results(
        LinearDiagramOptions(collinearity_blocks=CollinearityResult(blocks=(), unblocked_anchors=(untagged,))),
        (reversed_binding, other),
    ).options.collinearity_blocks.unblocked_anchors[0]
    assert (reversed_untagged.query_start, reversed_untagged.query_end, reversed_untagged.query_strand) == (
        198, 181, None,
    )


def test_a_projected_typed_result_is_keyed_by_identity_and_rebuilds_on_a_changed_input(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    # OV-1582: a stage that takes a projected orthogroup keys it by an identity
    # from what it was projected from, not by walking every member on each
    # warm render; a changed projection input gives another identity (no stale hit).
    from gbdraw import render_stages
    from gbdraw.analysis.protein_colinearity import OrthogroupResult
    from gbdraw.api.prepared import (
        PreparedBiologicalInputCache,
        PreparedResourceIdentity,
        prepared_value_identity,
        register_prepared_resource_value,
    )
    from gbdraw.api.record_planning import (
        ComparisonRecordBinding,
        comparison_record_binding,
        project_source_bound_results,
    )
    from gbdraw.features.source import capture_source_features
    from gbdraw.render_stages import render_stage
    from gbdraw.session_request_codec import encode_canonical_typed_resource

    record = _frame_test_record()
    binding = comparison_record_binding(record, capture_source_features(record), False)
    other = ComparisonRecordBinding({}, {}, {}, 50, {}, {}, False)
    reversed_binding = ComparisonRecordBinding(
        binding.source, binding.opposite, binding.current, binding.length,
        binding.opposite_spans, binding.current_spans, True,
    )
    members = [_member(index, binding.opposite_spans[index]) for index in (0, 1, 2)]
    groups = OrthogroupResult(orthogroups={"og1": members},
                              member_by_protein_id={member.protein_id: member for member in members})
    resource = tmp_path / "orthogroups.bin"
    resource.write_text("stored orthogroups", encoding="utf-8")
    identities = {resource: PreparedResourceIdentity("orthogroups", "token-orthogroups", 18)}
    cache = PreparedBiologicalInputCache()
    calls = {"key": 0}
    original_key = render_stages._key

    def counting_key(*args, **kwargs):
        calls["key"] += 1
        return original_key(*args, **kwargs)

    monkeypatch.setattr(render_stages, "_key", counting_key)

    def project(bindings: tuple[ComparisonRecordBinding, ...]):
        diagnostics: dict[str, object] = {"metrics": {}}
        with cache.transaction(resource_paths=identities, diagnostics=diagnostics):
            register_prepared_resource_value(groups, resource)
            projected = project_source_bound_results(LinearDiagramOptions(orthogroups=groups), bindings)
            orthogroups = projected.options.orthogroups
            identity = prepared_value_identity(orthogroups)
            calls["key"] = 0
            # A consumer stage, as the derived payload encodes the result.
            render_stage("derived_losat_entries", encode_canonical_typed_resource, "orthogroupResult", orthogroups)
            consumer_keys = calls["key"]
        return projected, identity, consumer_keys, diagnostics["metrics"]

    cold, cold_identity, _, _ = project((binding, other))
    warm, warm_identity, consumer_keys, warm_metrics = project((binding, other))
    expected = project_source_bound_results(LinearDiagramOptions(orthogroups=groups), (binding, other))

    assert repr(cold.options.orthogroups) == repr(warm.options.orthogroups) == repr(expected.options.orthogroups)
    assert cold.options.orthogroups is not groups  # the members moved
    assert warm.unmatched == cold.unmatched == expected.unmatched
    assert warm_identity is not None and warm_identity == cold_identity
    assert consumer_keys <= 2  # the identity, not the members
    assert warm_metrics["renderStageHit:derived_losat_entries"] == 1
    assert "renderStageBuild:derived_losat_entries" not in warm_metrics

    # Another frame: another identity and the projection of that frame.
    flipped, flipped_identity, _, _ = project((reversed_binding, other))
    assert flipped_identity != warm_identity
    assert flipped.options.orthogroups is groups  # already in that frame
    # Another orthogroup resource: rebuilt, never the stored projection.
    identities[resource] = PreparedResourceIdentity("orthogroups", "token-orthogroups-2", 18)
    rebuilt, rebuilt_identity, _, _ = project((binding, other))
    assert rebuilt_identity != warm_identity
    assert repr(rebuilt.options.orthogroups) == repr(expected.options.orthogroups)


@pytest.mark.parametrize("cropped", [False, True])
def test_web_protein_map_binding_equals_the_planned_search_frame_binding(tmp_path: Path, cropped: bool) -> None:
    # OV-1145: the Web builds each record's binding from its protein map with
    # no second plan; for the CDS keys it equals the binding the planner makes
    # for a reverse-complemented (and cropped) record.
    import json

    from gbdraw.api.record_planning import comparison_record_bindings
    from gbdraw.api.request_render import plan_linear_request
    from tests.test_api_request_render import _load_web_helpers

    record = _frame_test_record()
    record.features.append(SeqFeature(FeatureLocation(120, 165, strand=-1), type="CDS"))
    path = tmp_path / "frame.gbk"
    _write_genbank(path, record)
    region = "frame:11-180" if cropped else None
    planned = plan_linear_request(LinearDiagramRequest(records=(RecordInput(
        GenBankInputSource(path), record_key="frame",
        region=dataclasses.replace(parse_region_spec(region), reverse_complement=True) if region else None,
        presentation=RecordPresentation(reverse_complement=not cropped),
    ),)))
    (expected,) = comparison_record_bindings(planned.records, planned.provenance, search_frame=True)
    assert expected.reverse

    helpers = _load_web_helpers()
    extracted = json.loads(str(helpers["extract_cds_protein_fasta"](
        str(path), "genbank", None, region, None, "0", 0, "frame",
    )))
    (actual,) = helpers["_web_protein_map_bindings"]([{
        "record_index": 0, "protein_map": extracted["protein_map"],
        "view_transform": {"length": extracted["record_length"], "reverse": True},
    }])
    keys = sorted(actual.current)
    assert len(keys) == (3 if cropped else 4) and set(keys) <= set(expected.current)
    for field in ("current", "opposite", "current_spans", "opposite_spans"):
        assert {key: getattr(actual, field)[key] for key in keys} == {
            key: getattr(expected, field)[key] for key in keys
        }, field
    assert (actual.length, actual.reverse) == (expected.length, expected.reverse)
