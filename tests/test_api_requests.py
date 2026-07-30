from __future__ import annotations

import ast
from dataclasses import fields
from pathlib import Path

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gbdraw.api.options import (
    CircularDiagramOptions,
    CircularMultiRecordOptions,
    CircularOutputOptions,
    LinearDiagramOptions,
    LinearOutputOptions,
)
from gbdraw.api.requests import (
    CircularBatchOutputPolicy,
    CircularBatchRequest,
    CircularDiagramRequest,
    GenBankInputSource,
    GffFastaInputSource,
    InMemoryRecordSource,
    LinearDiagramRequest,
    RecordCollectionOptions,
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


def test_new_request_fields_preserve_released_positional_arguments() -> None:
    selector = parse_record_selector("#1")
    presentation = RecordPresentation(label="Shown")
    record = RecordInput(
        GenBankInputSource("record.gbk"),
        selector,
        None,
        presentation,
        "record-key",
    )

    assert record.selector is selector
    assert record.region is None
    assert record.presentation is presentation
    assert record.record_key == "record-key"
    assert record.cardinality.value == "exactly_one"

    output = RenderOutputRequest()
    request = CircularDiagramRequest(
        (record,),
        CircularDiagramOptions(),
        None,
        output,
        "single",
    )

    assert request.grouping == "single"
    assert request.record_options == RecordCollectionOptions()
    assert fields(CircularDiagramOptions)[-1].name == "conservation_table_file"
    assert fields(LinearDiagramOptions)[-1].name == "comparison_table_file"


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
    assert request.grouping == "grid"
    assert request.layout == CircularMultiRecordOptions()


def test_circular_request_allows_explicit_one_record_grid() -> None:
    request = CircularDiagramRequest(
        records=(_record("a.gbk"),),
        grouping="grid",
    )

    assert request.grouping == "grid"
    assert request.layout == CircularMultiRecordOptions()


def test_circular_batch_requires_one_unique_output_per_record() -> None:
    records = (_record("a.gbk"), _record("b.gbk"))
    outputs = (
        RenderOutputRequest(output_prefix="a"),
        RenderOutputRequest(output_prefix="b"),
    )

    request = CircularBatchRequest(records=records, outputs=outputs)

    assert request.grouping == "batch"
    assert request.outputs == outputs
    with pytest.raises(ValidationError, match="one resolved output"):
        CircularBatchRequest(records=records, outputs=outputs[:1])
    with pytest.raises(ValidationError, match="must be unique"):
        CircularBatchRequest(records=records, outputs=(outputs[0], outputs[0]))


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


def test_requests_require_mode_specific_diagram_options() -> None:
    with pytest.raises(ValidationError, match="CircularDiagramOptions"):
        CircularDiagramRequest(
            records=(_record(),),
            options=LinearDiagramOptions(),  # type: ignore[arg-type]
        )
    with pytest.raises(ValidationError, match="LinearDiagramOptions"):
        LinearDiagramRequest(
            records=(_record(),),
            options=CircularDiagramOptions(),  # type: ignore[arg-type]
        )


def test_mode_specific_option_fields_do_not_overlap_other_mode_features() -> None:
    circular_fields = {item.name for item in fields(CircularDiagramOptions)}
    linear_fields = {item.name for item in fields(LinearDiagramOptions)}

    assert {"blast_files", "protein_blastp_mode"}.isdisjoint(circular_fields)
    assert {
        "conservation_blast_files",
        "conservation_fasta_files",
        "conservation_reference",
    }.isdisjoint(linear_fields)


def test_mode_specific_output_options_reject_other_mode_title_positions() -> None:
    with pytest.raises(ValidationError, match="Circular"):
        CircularOutputOptions(plot_title_position="center")
    with pytest.raises(ValidationError, match="Linear"):
        LinearOutputOptions(plot_title_position="none")


@pytest.mark.parametrize(
    ("options_type", "path", "value"),
    [
        (CircularDiagramOptions, "canvas.circular.track_type", "middle"),
        (LinearDiagramOptions, "canvas.linear.track_layout", "above"),
    ],
)
def test_mode_specific_options_accept_owned_config_override_paths(
    options_type,
    path: str,
    value: object,
) -> None:
    options = options_type(config_overrides={path: value})

    assert options.config_overrides == {path: value}


@pytest.mark.parametrize(
    ("options_type", "path", "value", "other_mode"),
    [
        (
            CircularDiagramOptions,
            "canvas.linear.track_layout",
            "above",
            "Linear",
        ),
        (
            CircularDiagramOptions,
            "objects.axis.linear.stroke_color",
            "black",
            "Linear",
        ),
        (
            LinearDiagramOptions,
            "canvas.circular.track_type",
            "middle",
            "Circular",
        ),
    ],
)
def test_mode_specific_options_reject_other_mode_config_override_paths(
    options_type,
    path: str,
    value: object,
    other_mode: str,
) -> None:
    with pytest.raises(
        ValidationError,
        match=rf"cannot target {other_mode} settings: {path}",
    ):
        options_type(config_overrides={path: value})


@pytest.mark.parametrize(
    "kwargs",
    [
        {"conservation_reference": "invalid"},
        {"conservation_ring_width": 0},
        {"conservation_ring_gap": float("inf")},
    ],
)
def test_circular_options_reject_invalid_mode_values(
    kwargs: dict[str, object],
) -> None:
    with pytest.raises(ValidationError):
        CircularDiagramOptions(**kwargs)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"protein_blastp_mode": "invalid"},
        {"pairwise_match_style": "invalid"},
        {"collinearity_unit_mode": "invalid"},
        {"collinearity_anchor_mode": "invalid"},
        {"collinearity_search_scope": "invalid"},
        {"collinearity_color_mode": "invalid"},
        {"orthogroup_membership_mode": "invalid"},
        {"protein_blastp_max_hits": 0},
        {"losatp_threads": 0},
    ],
)
def test_linear_options_reject_invalid_mode_values(
    kwargs: dict[str, object],
) -> None:
    with pytest.raises(ValidationError):
        LinearDiagramOptions(**kwargs)


def test_linear_options_normalize_supported_mode_aliases() -> None:
    options = LinearDiagramOptions(
        collinearity_anchor_mode="top1",
        collinearity_color_mode="identity",
        orthogroup_membership_mode="rbh",
    )

    assert options.collinearity_anchor_mode == "one_to_one"
    assert options.collinearity_color_mode == "average_identity"
    assert options.orthogroup_membership_mode == "anchor_core_v1"


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
        {"output_prefix": "diagram\x00hidden"},
        {"output_prefix": "diagram\nhidden"},
        {"output_prefix": "NUL"},
        {"output_prefix": "CON.txt"},
        {"output_prefix": "COM¹"},
        {"output_prefix": "diagram:stream"},
        {"output_prefix": "diagram?draft"},
        {"output_directory": "results\x00hidden"},
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


def test_circular_batch_output_policy_rejects_control_characters() -> None:
    with pytest.raises(ValidationError):
        CircularBatchOutputPolicy(output_prefix="diagram\x00hidden")
