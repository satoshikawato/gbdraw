from __future__ import annotations

import ast
from pathlib import Path
from types import SimpleNamespace

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from pandas import DataFrame
from svgwrite import Drawing

import gbdraw.api.request_render as request_render_module
from gbdraw.analysis.collinearity import LosslessCollinearityParameters
from gbdraw.analysis.protein_colinearity import (
    extract_web_stable_cds_proteins,
)
from gbdraw.api.options import (
    CircularDiagramOptions,
    CircularMultiRecordOptions,
    LinearDiagramOptions,
    LinearMultiRecordOptions,
)
from gbdraw.api.request_render import (
    CircularBatchRenderResult,
    CircularRequestPlan,
    LinearRequestPlan,
    PreparedDiagramRequest,
    build_request_diagram,
    normalize_request_records,
    plan_circular_request,
    plan_linear_request,
    render_request,
)
from gbdraw.api.requests import (
    CircularBatchRequest,
    CircularDiagramRequest,
    GffFastaInputSource,
    InMemoryRecordSource,
    LinearDiagramRequest,
    RecordInput,
    RecordPresentation,
    RenderOutputRequest,
)
from gbdraw.exceptions import ExportError, ValidationError
from gbdraw.io.record_select import parse_record_selector
from gbdraw.io.regions import parse_region_spec
from gbdraw.session_io import validate_current_session_artifacts


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


def _protein_record(record_id: str, protein_id: str) -> SeqRecord:
    record = _seqrecord(record_id, "ATGAAATAG" * 4)
    record.features = [
        SeqFeature(
            FeatureLocation(0, 9, strand=1),
            type="CDS",
            qualifiers={
                "translation": ["MK"],
                "protein_id": [protein_id],
            },
        )
    ]
    return record


def test_rewrite_protein_artifact_references_rejects_key_collisions() -> None:
    with pytest.raises(ValidationError, match="duplicate key"):
        request_render_module.rewrite_protein_artifact_references(
            {"legacy-a": 1, "legacy-b": 2},
            {"legacy-a": "h_same", "legacy-b": "h_same"},
        )


def test_rewrite_protein_artifact_references_updates_compound_legacy_ids() -> None:
    query_id = "p_r_old_0_9_1_deadbeefdead"
    subject_id = "p_r_other_10_19_1_cafebabecafe"

    rewritten = request_render_module.rewrite_protein_artifact_references(
        {
            "query_protein_id": f"{query_id};{subject_id}",
            "supporting_edge": f"{query_id}->{subject_id}:rbh",
        },
        {
            query_id: "h_aaaaaaaaaaaaaaaaaaaaaaaaaa",
            subject_id: "h_bbbbbbbbbbbbbbbbbbbbbbbbbb",
        },
    )

    assert rewritten == {
        "query_protein_id": (
            "h_aaaaaaaaaaaaaaaaaaaaaaaaaa;h_bbbbbbbbbbbbbbbbbbbbbbbbbb"
        ),
        "supporting_edge": (
            "h_aaaaaaaaaaaaaaaaaaaaaaaaaa->h_bbbbbbbbbbbbbbbbbbbbbbbbbb:rbh"
        ),
    }


def _derived_identity_test_context(
    *,
    options: LinearDiagramOptions | None = None,
) -> SimpleNamespace:
    records = (_seqrecord("a"), _seqrecord("b"))
    request = LinearDiagramRequest(
        records=(
            RecordInput(
                source=InMemoryRecordSource(records[0]),
                record_key="record-1",
            ),
            RecordInput(
                source=InMemoryRecordSource(records[1]),
                record_key="record-2",
            ),
        ),
        options=options or LinearDiagramOptions(),
    )
    drawing = Drawing("out.svg")
    drawing._gbdraw_resolved_protein_comparisons = [  # type: ignore[attr-defined]
        DataFrame(columns=request_render_module.COMPARISON_COLUMNS)
    ]
    manifest = SimpleNamespace(
        record_instances={
            "record-1": {
                "runtimeBindingHash": "sha256:runtime-1",
                "displayBindingHash": "sha256:display-1",
            },
            "record-2": {
                "runtimeBindingHash": "sha256:runtime-2",
                "displayBindingHash": "sha256:display-2",
            },
        }
    )
    artifacts = request_render_module._PreparedLinearArtifacts(
        request=request,
        cache=None,
        extraction=SimpleNamespace(identity_manifest=manifest),
        nucleotide_entries=(),
        passthrough_derived_entries=(),
        legacy_candidates=(),
        protein_id_map={},
        source_mode="collinear",
        warnings=(),
    )
    raw_entries = (
        {
            "key": "visible",
            "queryRecordInstanceKey": "record-1",
            "subjectRecordInstanceKey": "record-2",
        },
        {
            "key": "reverse-hidden",
            "queryRecordInstanceKey": "record-2",
            "subjectRecordInstanceKey": "record-1",
        },
        {
            "key": "self-hidden",
            "queryRecordInstanceKey": "record-1",
            "subjectRecordInstanceKey": "record-1",
        },
    )
    return SimpleNamespace(
        records=records,
        request=request,
        drawing=drawing,
        artifacts=artifacts,
        raw_entries=raw_entries,
    )


def _derived_entry_for_params(
    params: LosslessCollinearityParameters | None,
) -> dict[str, object]:
    context = _derived_identity_test_context(
        options=LinearDiagramOptions(collinearity_params=params)
    )
    (entry,) = request_render_module._build_current_derived_entries(
        context.drawing,
        context.request,
        context.records,
        context.artifacts,
        context.raw_entries,
    )
    return dict(entry)


def test_derived_identity_includes_hidden_reverse_and_self_raw_inputs(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    context = _derived_identity_test_context()
    monkeypatch.setattr(
        request_render_module,
        "is_protein_losat_cache_entry",
        lambda _entry: True,
    )

    (entry,) = request_render_module._build_current_derived_entries(
        context.drawing,
        context.request,
        context.records,
        context.artifacts,
        context.raw_entries,
    )

    identity = entry["payload"]["identity"]
    assert identity["rawCacheKeys"] == [
        "reverse-hidden",
        "self-hidden",
        "visible",
    ]
    assert identity["pairs"][0]["rawCacheKeys"] == ["visible"]

    changed_entries = tuple(
        {**raw_entry, "key": "reverse-recomputed"}
        if raw_entry["key"] == "reverse-hidden"
        else raw_entry
        for raw_entry in context.raw_entries
    )
    (changed,) = request_render_module._build_current_derived_entries(
        context.drawing,
        context.request,
        context.records,
        context.artifacts,
        changed_entries,
    )
    assert changed["key"] != entry["key"]


@pytest.mark.parametrize(
    "changed_params",
    (
        LosslessCollinearityParameters(min_anchors=2),
        LosslessCollinearityParameters(max_unit_gap=1),
        LosslessCollinearityParameters(max_diagonal_drift=1),
        LosslessCollinearityParameters(max_conflicts=1),
        LosslessCollinearityParameters(merge_orientation="strand"),
    ),
    ids=(
        "min-anchors",
        "max-unit-gap",
        "max-diagonal-drift",
        "max-conflicts",
        "merge-orientation",
    ),
)
def test_derived_identity_changes_with_each_lossless_collinearity_parameter(
    monkeypatch: pytest.MonkeyPatch,
    changed_params: LosslessCollinearityParameters,
) -> None:
    monkeypatch.setattr(
        request_render_module,
        "is_protein_losat_cache_entry",
        lambda _entry: True,
    )
    baseline = _derived_entry_for_params(LosslessCollinearityParameters())
    unchanged = _derived_entry_for_params(LosslessCollinearityParameters())
    changed = _derived_entry_for_params(changed_params)

    collinear = baseline["payload"]["identity"]["collinear"]
    assert {
        name: collinear[name]
        for name in (
            "minAnchors",
            "maxGeneGap",
            "maxDiagonalDrift",
            "maxConflictsInMergeGap",
            "mergeOrientation",
        )
    } == {
        "minAnchors": 1,
        "maxGeneGap": 0,
        "maxDiagonalDrift": 0,
        "maxConflictsInMergeGap": 0,
        "mergeOrientation": "either",
    }
    assert collinear["parameterIdentity"] == {
        "model": "lossless",
        "parameters": {
            "minAnchors": 1,
            "maxUnitGap": 0,
            "maxDiagonalDrift": 0,
            "maxConflicts": 0,
            "mergeOrientation": "either",
        },
    }
    assert unchanged["key"] == baseline["key"]
    assert changed["key"] != baseline["key"]

@pytest.mark.parametrize("mode", ("orthogroup", "collinear"))
def test_empty_api_derived_result_passes_current_session_validation(
    mode: str,
) -> None:
    records = (
        _protein_record("source-a", "protein-a"),
        _protein_record("source-b", "protein-b"),
    )
    record_keys = ("record-1", "record-2")
    extraction = extract_web_stable_cds_proteins(
        records,
        record_instance_keys=record_keys,
        record_source_ids=tuple(record.id for record in records),
    )
    request = LinearDiagramRequest(
        records=tuple(
            RecordInput(
                source=InMemoryRecordSource(record),
                record_key=record_key,
            )
            for record, record_key in zip(records, record_keys, strict=True)
        ),
        options=LinearDiagramOptions(protein_blastp_mode=mode),
    )
    drawing = Drawing("out.svg")
    drawing._gbdraw_resolved_protein_comparisons = [  # type: ignore[attr-defined]
        DataFrame(columns=request_render_module.COMPARISON_COLUMNS)
    ]
    artifacts = request_render_module._PreparedLinearArtifacts(
        request=request,
        cache=None,
        extraction=extraction,
        nucleotide_entries=(),
        passthrough_derived_entries=(),
        legacy_candidates=(),
        protein_id_map={},
        source_mode=mode,
        warnings=(),
    )

    (entry,) = request_render_module._build_current_derived_entries(
        drawing,
        request,
        records,
        artifacts,
        (),
    )

    assert entry["payload"]["pairs"][0]["rows"] == []
    assert entry["payload"]["pairs"][0]["tsv"] == ""
    assert entry["payload"]["pairs"][0]["hit_count"] == 0
    assert entry["payload"]["orthogroups"] == []
    manifest = extraction.identity_manifest
    assert manifest is not None
    validate_current_session_artifacts(
        {
            "losatCache": {"entries": []},
            "losatDerivedCache": {"entries": [entry]},
            "proteinIdentityManifest": manifest.to_dict(),
        }
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


@pytest.mark.parametrize("module_name", ("circular.py", "linear.py"))
def test_cli_renderers_do_not_call_assemblers_or_export_directly(
    module_name: str,
) -> None:
    source = (Path(request_render_module.__file__).parents[1] / module_name).read_text(
        encoding="utf-8"
    )
    tree = ast.parse(source)
    called_names = {
        node.func.id
        for node in ast.walk(tree)
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
    }

    assert "render_request" in called_names
    assert "save_figure" not in called_names
    assert not any(name.startswith("assemble_") for name in called_names)


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


def test_linear_planner_applies_comparison_selection_after_neutral_gff_load() -> None:
    fixture_dir = Path(__file__).parents[1] / "examples" / "gff3_lambda"
    request = LinearDiagramRequest(
        records=(
            RecordInput(
                source=GffFastaInputSource(
                    fixture_dir / "lambda_two_contigs.gff3",
                    fixture_dir / "lambda_two_contigs.fna",
                )
            ),
        ),
        options=LinearDiagramOptions(blast_files=("comparison.tsv",)),
    )

    assert normalize_request_records(request)[0].id == "lambda_left"


@pytest.mark.parametrize(
    ("topology", "warning"),
    (("circular", None), (None, None), ("linear", "is linear"), ("unknown", "not available")),
)
def test_circular_planner_owns_topology_warnings(
    caplog: pytest.LogCaptureFixture,
    topology: str | None,
    warning: str | None,
) -> None:
    record = _seqrecord("topology")
    if topology is not None:
        record.annotations["topology"] = topology

    plan_circular_request(
        CircularDiagramRequest(
            records=(RecordInput(source=InMemoryRecordSource(record)),),
        )
    )

    assert (warning is not None) == (warning in caplog.text if warning else False)


def test_mode_planners_return_validated_plan_types(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    circular_record = _seqrecord("circular")
    linear_record = _seqrecord("linear")
    circular_request = CircularDiagramRequest(
        records=(_memory_input("circular"),),
        options=CircularDiagramOptions(),
        layout=CircularMultiRecordOptions(
            multi_record_positions=("#1@1",),
        ),
    )
    linear_request = LinearDiagramRequest(
        records=(_memory_input("linear"),),
        options=LinearDiagramOptions(),
        layout=LinearMultiRecordOptions(),
    )

    monkeypatch.setattr(
        request_render_module,
        "normalize_request_records",
        lambda request: (
            circular_record
            if isinstance(request, CircularDiagramRequest)
            else linear_record,
        ),
    )

    circular_plan = plan_circular_request(circular_request)
    linear_plan = plan_linear_request(linear_request)

    assert isinstance(circular_plan, CircularRequestPlan)
    assert circular_plan.mode == "circular"
    assert circular_plan.records == (circular_record,)
    assert circular_plan.layout == circular_request.layout
    assert isinstance(linear_plan, LinearRequestPlan)
    assert linear_plan.mode == "linear"
    assert linear_plan.records == (linear_record,)
    assert linear_plan.layout == linear_request.layout


def test_public_plan_values_reject_inconsistent_construction() -> None:
    circular_request = CircularDiagramRequest(records=(_memory_input("circular"),))
    circular_record = _seqrecord("circular")
    linear_record = _seqrecord("linear")

    with pytest.raises(ValidationError, match="record count"):
        CircularRequestPlan(
            request=circular_request,
            records=(circular_record, _seqrecord("extra")),
            layout=None,
        )
    with pytest.raises(ValidationError, match="layout presence"):
        CircularRequestPlan(
            request=circular_request,
            records=(circular_record,),
            layout=CircularMultiRecordOptions(),
        )
    with pytest.raises(ValidationError, match="LinearDiagramRequest"):
        LinearRequestPlan(
            request=circular_request,  # type: ignore[arg-type]
            records=(linear_record,),
            layout=None,
        )


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
    monkeypatch.setattr(
        request_render_module,
        "_interactive_context",
        lambda _, **__: context,
    )

    def fake_save(canvas, formats, **kwargs):
        captured["canvas"] = canvas
        captured["formats"] = formats
        captured.update(kwargs)
        return [str(path) for path in expected_paths]

    monkeypatch.setattr(request_render_module, "save_figure_to", fake_save)

    result = render_request(request)

    assert result.mode == "linear"
    assert result.output_paths == tuple(expected_paths)
    assert result.interactive_context is context
    assert captured == {
        "canvas": drawing,
        "formats": ("svg", "interactive_svg"),
        "output_dir": str(tmp_path),
        "output_prefix": "diagram",
        "overwrite": True,
        "interactive_context": context,
    }


def test_render_request_fails_when_interactive_metadata_generation_fails(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    request = LinearDiagramRequest(
        records=(_memory_input("a"),),
        output=RenderOutputRequest(
            output_prefix="diagram",
            output_directory=tmp_path,
            formats=("interactive_svg",),
        ),
    )
    prepared = PreparedDiagramRequest(
        mode="linear",
        request=request,
        records=(_seqrecord("a"),),
        drawing=Drawing("out.svg"),
    )
    monkeypatch.setattr(
        request_render_module,
        "build_request_diagram",
        lambda _: prepared,
    )

    def fail_context(*_args, **_kwargs):
        raise RuntimeError("metadata exploded")

    monkeypatch.setattr(
        request_render_module,
        "build_interactive_svg_context",
        fail_context,
    )

    with pytest.raises(
        ExportError,
        match="Interactive SVG metadata generation failed: metadata exploded",
    ):
        render_request(request)

    assert not list(tmp_path.iterdir())


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


@pytest.mark.circular
def test_render_request_preserves_dotted_output_prefix(tmp_path: Path) -> None:
    record = _seqrecord("dotted-prefix", "ATGCGC" * 200)
    record.annotations["topology"] = "circular"
    request = CircularDiagramRequest(
        records=(RecordInput(source=InMemoryRecordSource(record)),),
        output=RenderOutputRequest(
            output_prefix="sample.v1",
            output_directory=tmp_path,
            formats=("svg",),
        ),
    )

    result = render_request(request)

    assert result.output_paths == (tmp_path / "sample.v1.svg",)
    assert result.output_paths[0].is_file()
    assert not (tmp_path / "sample.svg").exists()


@pytest.mark.circular
def test_render_request_circular_batch_creates_one_output_per_record(
    tmp_path: Path,
) -> None:
    records = tuple(
        RecordInput(source=InMemoryRecordSource(_seqrecord(record_id, "ATGCGC" * 50)))
        for record_id in ("batch-a", "batch-b")
    )
    request = CircularBatchRequest(
        records=records,
        outputs=tuple(
            RenderOutputRequest(
                output_prefix=prefix,
                output_directory=tmp_path,
            )
            for prefix in ("batch-a", "batch-b")
        ),
    )

    result = render_request(request)

    assert isinstance(result, CircularBatchRenderResult)
    assert result.output_paths == (
        tmp_path / "batch-a.svg",
        tmp_path / "batch-b.svg",
    )
    assert all(path.is_file() for path in result.output_paths)


@pytest.mark.circular
def test_render_request_circular_batch_loads_comparison_fasta_once(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    comparison_fasta = tmp_path / "comparison.fna"
    comparison_fasta.write_text(">comparison\nAACCGG\n", encoding="utf-8")
    parse_calls: list[object] = []
    parse = request_render_module.SeqIO.parse

    def counting_parse(source, format_name):
        parse_calls.append(source)
        return parse(source, format_name)

    monkeypatch.setattr(request_render_module.SeqIO, "parse", counting_parse)
    request = CircularBatchRequest(
        records=(_memory_input("batch-a"), _memory_input("batch-b")),
        options=CircularDiagramOptions(
            conservation_fasta_files=(str(comparison_fasta),),
        ),
        outputs=tuple(
            RenderOutputRequest(
                output_prefix=prefix,
                output_directory=tmp_path,
                formats=("interactive_svg",),
            )
            for prefix in ("batch-a", "batch-b")
        ),
    )

    result = render_request(request)

    assert isinstance(result, CircularBatchRenderResult)
    assert parse_calls == [str(comparison_fasta)]
    assert all(
        any(
            source["origin"] == "homology-comparison"
            for source in context.sequence_sources
        )
        for context in result.interactive_contexts
        if context is not None
    )


@pytest.mark.parametrize(
    "formats",
    [("svg",), ("interactive_svg",)],
)
def test_render_request_only_requires_comparison_fasta_for_interactive_output(
    tmp_path: Path,
    formats: tuple[str, ...],
) -> None:
    request = CircularDiagramRequest(
        records=(_memory_input("record"),),
        options=CircularDiagramOptions(
            conservation_fasta_files=(str(tmp_path / "missing.fna"),),
        ),
        output=RenderOutputRequest(
            output_prefix="diagram",
            output_directory=tmp_path,
            formats=formats,
        ),
    )

    if "interactive_svg" in formats:
        with pytest.raises(
            ExportError,
            match="Interactive SVG metadata generation failed",
        ):
            render_request(request)
        assert not list(tmp_path.iterdir())
    else:
        result = render_request(request)
        assert result.output_paths


@pytest.mark.circular
def test_render_request_circular_batch_preflights_existing_outputs(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    existing = tmp_path / "batch-b.svg"
    existing.write_text("keep", encoding="utf-8")
    request = CircularBatchRequest(
        records=(_memory_input("batch-a"), _memory_input("batch-b")),
        outputs=tuple(
            RenderOutputRequest(
                output_prefix=prefix,
                output_directory=tmp_path,
            )
            for prefix in ("batch-a", "batch-b")
        ),
    )
    monkeypatch.setattr(
        request_render_module,
        "build_request_diagram",
        lambda _request: pytest.fail("batch outputs must be preflighted before building"),
    )

    with pytest.raises(ValidationError, match="already exist"):
        render_request(request)

    assert not (tmp_path / "batch-a.svg").exists()
    assert existing.read_text(encoding="utf-8") == "keep"
