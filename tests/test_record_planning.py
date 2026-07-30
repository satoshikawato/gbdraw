from __future__ import annotations

import ast
from pathlib import Path
from types import SimpleNamespace

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
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
from gbdraw.session_request_codec import (
    CanonicalRequestEncodingError,
    encode_canonical_request,
)


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


def test_schema5_rejects_unresolved_then_round_trips_projection(
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

    with pytest.raises(
        CanonicalRequestEncodingError,
        match="materialized exact-one request",
    ):
        encode_canonical_request(unresolved)

    resolved = resolve_request(unresolved)
    encoded = encode_canonical_request(resolved)
    assert len(resolved.records) == 2
    assert all(
        record.cardinality is RecordCardinality.EXACTLY_ONE
        for record in resolved.records
    )
    assert all("cardinality" not in row for row in encoded.payload["records"])

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
