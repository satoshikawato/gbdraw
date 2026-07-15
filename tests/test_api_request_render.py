from __future__ import annotations

import ast
from pathlib import Path

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from svgwrite import Drawing

import gbdraw.api.request_render as request_render_module
from gbdraw.api.request_render import (
    PreparedDiagramRequest,
    build_request_diagram,
    normalize_request_records,
    render_request,
)
from gbdraw.api.requests import (
    CircularDiagramRequest,
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


def _seqrecord(record_id: str, sequence: str = "AACCGG") -> SeqRecord:
    return SeqRecord(
        Seq(sequence),
        id=record_id,
        annotations={"molecule_type": "DNA"},
    )


def _memory_input(
    record_id: str,
    *,
    presentation: RecordPresentation | None = None,
) -> RecordInput:
    return RecordInput(
        source=InMemoryRecordSource(_seqrecord(record_id)),
        presentation=presentation or RecordPresentation(),
    )


def test_request_render_module_does_not_import_cli_or_session_owners() -> None:
    source_path = Path(request_render_module.__file__)
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


def test_normalize_in_memory_record_applies_region_and_presentation_without_mutation() -> None:
    source_record = _seqrecord("record-a")
    request = LinearDiagramRequest(
        records=(
            RecordInput(
                source=InMemoryRecordSource(source_record),
                region=parse_region_spec("1-4:rc"),
                presentation=RecordPresentation(label=" Label A ", subtitle=" Sub A "),
            ),
        )
    )

    (normalized,) = normalize_request_records(request)

    assert str(normalized.seq) == "GGTT"
    assert normalized.annotations["gbdraw_record_label"] == "Label A"
    assert normalized.annotations["gbdraw_record_subtitle"] == "Sub A"
    assert str(source_record.seq) == "AACCGG"
    assert "gbdraw_record_label" not in source_record.annotations


def test_normalize_in_memory_record_validates_selector() -> None:
    request = LinearDiagramRequest(
        records=(
            RecordInput(
                source=InMemoryRecordSource(_seqrecord("record-a")),
                selector=parse_record_selector("missing"),
            ),
        )
    )

    with pytest.raises(ValidationError, match="did not match"):
        normalize_request_records(request)


def test_normalize_gff_fasta_source_uses_selector() -> None:
    fixture_dir = Path(__file__).parents[1] / "examples" / "gff3_lambda"
    request = LinearDiagramRequest(
        records=(
            RecordInput(
                source=GffFastaInputSource(
                    fixture_dir / "lambda_two_contigs.gff3",
                    fixture_dir / "lambda_two_contigs.fna",
                ),
                selector=parse_record_selector("lambda_left"),
            ),
        ),
    )

    (record,) = normalize_request_records(request)

    assert record.id == "lambda_left"
    assert any(feature.type == "CDS" for feature in record.features)


def test_normalize_record_input_requires_one_resolved_record() -> None:
    fixture_dir = Path(__file__).parents[1] / "examples" / "gff3_lambda"
    request = LinearDiagramRequest(
        records=(
            RecordInput(
                source=GffFastaInputSource(
                    fixture_dir / "lambda_two_contigs.gff3",
                    fixture_dir / "lambda_two_contigs.fna",
                )
            ),
        )
    )

    with pytest.raises(ValidationError, match="exactly one record"):
        normalize_request_records(request)


def test_build_circular_request_derives_row_and_column_order(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records = (_seqrecord("a"), _seqrecord("b"), _seqrecord("c"))
    request = CircularDiagramRequest(
        records=(
            _memory_input("a", presentation=RecordPresentation(grid_row=2, grid_column=1)),
            _memory_input("b", presentation=RecordPresentation(grid_row=1, grid_column=2)),
            _memory_input("c", presentation=RecordPresentation(grid_row=1, grid_column=1)),
        )
    )
    captured: dict[str, object] = {}
    drawing = Drawing("out.svg")

    monkeypatch.setattr(request_render_module, "normalize_request_records", lambda _: records)

    def fake_build(loaded_records, *, options, layout):
        captured["records"] = loaded_records
        captured["options"] = options
        captured["layout"] = layout
        return drawing

    monkeypatch.setattr(request_render_module, "build_circular_multi_diagram", fake_build)

    prepared = build_request_diagram(request)

    assert prepared.drawing is drawing
    assert prepared.mode == "circular"
    assert captured["records"] == records
    assert captured["layout"].multi_record_positions == ("#3@1", "#2@1", "#1@2")


def test_build_linear_request_uses_high_level_builder(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records = (_seqrecord("a"), _seqrecord("b"))
    request = LinearDiagramRequest(records=(_memory_input("a"), _memory_input("b")))
    drawing = Drawing("out.svg")
    captured: dict[str, object] = {}

    monkeypatch.setattr(request_render_module, "normalize_request_records", lambda _: records)

    def fake_build(loaded_records, *, options):
        captured["records"] = loaded_records
        captured["options"] = options
        return drawing

    monkeypatch.setattr(request_render_module, "build_linear_diagram", fake_build)

    prepared = build_request_diagram(request)

    assert prepared.mode == "linear"
    assert prepared.drawing is drawing
    assert captured == {"records": records, "options": request.options}


def test_render_request_passes_output_policy_and_returns_existing_paths(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    request = LinearDiagramRequest(
        records=(_memory_input("a"),),
        output=RenderOutputRequest(
            output_prefix="diagram",
            output_directory=tmp_path,
            formats=("svg", "interactive_svg"),
            overwrite=True,
        ),
    )
    records = (_seqrecord("a"),)
    drawing = Drawing("out.svg")
    prepared = PreparedDiagramRequest(
        mode="linear",
        request=request,
        records=records,
        drawing=drawing,
    )
    expected_paths = [tmp_path / "diagram.svg", tmp_path / "diagram.interactive.svg"]
    captured: dict[str, object] = {}
    context = object()

    monkeypatch.setattr(request_render_module, "build_request_diagram", lambda _: prepared)
    monkeypatch.setattr(request_render_module, "_interactive_context", lambda _: context)

    def fake_save(canvas, formats, **kwargs):
        captured["canvas"] = canvas
        captured["formats"] = formats
        captured.update(kwargs)
        return [str(path) for path in expected_paths]

    monkeypatch.setattr(request_render_module, "save_figure_to", fake_save)

    result = render_request(request)

    assert result.mode == "linear"
    assert result.output_paths == tuple(expected_paths)
    assert captured == {
        "canvas": drawing,
        "formats": ("svg", "interactive_svg"),
        "output_dir": str(tmp_path),
        "output_prefix": "diagram",
        "overwrite": True,
        "interactive_context": context,
    }


@pytest.mark.circular
def test_render_request_circular_smoke_creates_svg(tmp_path: Path) -> None:
    record = _seqrecord("request-smoke", "ATGCGC" * 200)
    record.annotations["topology"] = "circular"
    request = CircularDiagramRequest(
        records=(RecordInput(source=InMemoryRecordSource(record)),),
        output=RenderOutputRequest(
            output_prefix="request-smoke",
            output_directory=tmp_path,
            formats=("svg",),
        ),
    )

    result = render_request(request)

    assert result.mode == "circular"
    assert result.output_paths == (tmp_path / "request-smoke.svg",)
    assert result.output_paths[0].is_file()
