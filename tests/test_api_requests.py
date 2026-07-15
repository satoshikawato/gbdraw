from __future__ import annotations

import ast
from pathlib import Path

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gbdraw.api.options import CircularMultiRecordOptions, DiagramOptions
from gbdraw.api.requests import (
    CircularDiagramRequest,
    GenBankInputSource,
    GffFastaInputSource,
    InMemoryRecordSource,
    LinearDiagramRequest,
    RecordInput,
    RecordPresentation,
    RenderOutputRequest,
)
from gbdraw.exceptions import ValidationError
from gbdraw.io.record_select import parse_record_selector
from gbdraw.io.regions import parse_region_spec


def _record(
    name: str = "record.gbk",
    *,
    presentation: RecordPresentation | None = None,
) -> RecordInput:
    return RecordInput(
        source=GenBankInputSource(name),
        presentation=presentation or RecordPresentation(),
    )


def test_request_module_does_not_import_cli_or_session_owners() -> None:
    source_path = Path(__file__).parents[1] / "gbdraw" / "api" / "requests.py"
    tree = ast.parse(source_path.read_text(encoding="utf-8"))
    imported_modules = {
        node.module
        for node in ast.walk(tree)
        if isinstance(node, ast.ImportFrom) and node.module is not None
    }
    imported_modules.update(
        alias.name
        for node in ast.walk(tree)
        if isinstance(node, ast.Import)
        for alias in node.names
    )

    forbidden_prefixes = (
        "gbdraw.circular",
        "gbdraw.linear",
        "gbdraw.session_io",
        "gbdraw.cli",
    )
    assert not any(
        module == prefix or module.startswith(f"{prefix}.")
        for module in imported_modules
        for prefix in forbidden_prefixes
    )


def test_record_sources_normalize_materialized_paths() -> None:
    gbk = GenBankInputSource("inputs/record.gbk")
    pair = GffFastaInputSource(Path("record.gff3"), "record.fna")

    assert gbk.path == Path("inputs/record.gbk")
    assert pair.gff_path == Path("record.gff3")
    assert pair.fasta_path == Path("record.fna")


def test_record_source_accepts_an_in_memory_seqrecord() -> None:
    record = SeqRecord(Seq("ATGC"), id="record-a")

    assert InMemoryRecordSource(record).record is record


def test_in_memory_record_source_rejects_other_objects() -> None:
    with pytest.raises(ValidationError, match="SeqRecord"):
        InMemoryRecordSource(object())


@pytest.mark.parametrize(
    "factory",
    [
        lambda: GenBankInputSource(""),
        lambda: GenBankInputSource("."),
        lambda: GffFastaInputSource("record.gff3", ""),
    ],
)
def test_record_sources_reject_missing_paths(factory) -> None:
    with pytest.raises(ValidationError, match="materialized file"):
        factory()


def test_record_input_rejects_duplicate_selector_owners() -> None:
    with pytest.raises(ValidationError, match="either on RecordInput or in its region"):
        RecordInput(
            source=GenBankInputSource("record.gbk"),
            selector=parse_record_selector("record-a"),
            region=parse_region_spec("record-a:1-10"),
        )


def test_record_input_rejects_duplicate_reverse_complement() -> None:
    with pytest.raises(ValidationError, match="both the region and record presentation"):
        RecordInput(
            source=GenBankInputSource("record.gbk"),
            region=parse_region_spec("1-10:rc"),
            presentation=RecordPresentation(reverse_complement=True),
        )


@pytest.mark.parametrize(
    "presentation",
    [
        RecordPresentation(grid_row=1, grid_column=2),
        RecordPresentation(label="record A", subtitle="subtitle", reverse_complement=True),
    ],
)
def test_record_presentation_accepts_typed_values(presentation: RecordPresentation) -> None:
    assert isinstance(presentation, RecordPresentation)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"grid_row": 0},
        {"grid_row": True},
        {"grid_column": 1},
        {"reverse_complement": "yes"},
    ],
)
def test_record_presentation_rejects_invalid_values(kwargs: dict[str, object]) -> None:
    with pytest.raises(ValidationError):
        RecordPresentation(**kwargs)


def test_circular_request_normalizes_records_and_default_multi_layout() -> None:
    request = CircularDiagramRequest(records=[_record("a.gbk"), _record("b.gbk")])

    assert request.records == (_record("a.gbk"), _record("b.gbk"))
    assert request.layout == CircularMultiRecordOptions()


def test_circular_request_accepts_unique_record_grid() -> None:
    request = CircularDiagramRequest(
        records=(
            _record("a.gbk", presentation=RecordPresentation(grid_row=1, grid_column=1)),
            _record("b.gbk", presentation=RecordPresentation(grid_row=1, grid_column=2)),
        )
    )

    assert request.layout is not None


def test_circular_request_rejects_partial_or_duplicate_grid() -> None:
    partial = (
        _record("a.gbk", presentation=RecordPresentation(grid_row=1)),
        _record("b.gbk"),
    )
    with pytest.raises(ValidationError, match="every record"):
        CircularDiagramRequest(records=partial)

    duplicate = (
        _record("a.gbk", presentation=RecordPresentation(grid_row=1, grid_column=1)),
        _record("b.gbk", presentation=RecordPresentation(grid_row=1, grid_column=1)),
    )
    with pytest.raises(ValidationError, match="must be unique"):
        CircularDiagramRequest(records=duplicate)


def test_circular_request_rejects_two_placement_sources() -> None:
    with pytest.raises(ValidationError, match="or layout, not both"):
        CircularDiagramRequest(
            records=(
                _record("a.gbk", presentation=RecordPresentation(grid_row=1)),
                _record("b.gbk", presentation=RecordPresentation(grid_row=2)),
            ),
            layout=CircularMultiRecordOptions(multi_record_positions=("#1@1", "#2@2")),
        )


def test_linear_request_rejects_circular_grid_placement() -> None:
    with pytest.raises(ValidationError, match="only by circular"):
        LinearDiagramRequest(
            records=(_record(presentation=RecordPresentation(grid_row=1)),)
        )


def test_requests_reject_empty_record_sequences() -> None:
    with pytest.raises(ValidationError, match="at least one"):
        CircularDiagramRequest(records=())
    with pytest.raises(ValidationError, match="at least one"):
        LinearDiagramRequest(records=())


def test_requests_reuse_mode_specific_diagram_option_validation() -> None:
    with pytest.raises(ValidationError, match="blast_files"):
        CircularDiagramRequest(
            records=(_record(),),
            options=DiagramOptions(blast_files=("comparison.tsv",)),
        )
    with pytest.raises(ValidationError, match="conservation_blast_files"):
        LinearDiagramRequest(
            records=(_record(),),
            options=DiagramOptions(conservation_blast_files=("comparison.tsv",)),
        )


def test_render_output_request_normalizes_formats_and_paths() -> None:
    request = RenderOutputRequest(
        output_prefix="diagram",
        output_directory="results",
        formats="svg,interactive-svg,svg",
        overwrite=True,
        interactive_metadata_policy="required",
    )

    assert request.output_prefix == "diagram"
    assert request.output_directory == Path("results")
    assert request.formats == ("svg", "interactive_svg")


@pytest.mark.parametrize(
    "kwargs",
    [
        {"output_prefix": ""},
        {"output_prefix": "nested/diagram"},
        {"output_prefix": r"nested\diagram"},
        {"formats": ()},
        {"formats": ("unknown",)},
        {"overwrite": 1},
        {"interactive_metadata_policy": "sometimes"},
        {"formats": ("svg",), "interactive_metadata_policy": "required"},
    ],
)
def test_render_output_request_rejects_invalid_values(kwargs: dict[str, object]) -> None:
    with pytest.raises(ValidationError):
        RenderOutputRequest(**kwargs)
