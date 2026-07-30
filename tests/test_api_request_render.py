from __future__ import annotations

import ast
import inspect
from pathlib import Path
from types import SimpleNamespace

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from pandas import DataFrame
from svgwrite import Drawing

import gbdraw.api.request_render as request_render_module
import gbdraw.api.diagram as diagram_module
import gbdraw.api.io as api_io_module
import gbdraw.api.prepared as prepared_module
import gbdraw.configurators.features as feature_configurator_module
import gbdraw.render.interactive_context as interactive_context_module
from gbdraw.analysis.collinearity import LosslessCollinearityParameters
from gbdraw.analysis.protein_colinearity import (
    OrthogroupMember,
    OrthogroupResult,
    ProteinBlastpResult,
    extract_web_stable_cds_proteins,
)
from gbdraw.api.options import (
    CircularDiagramOptions,
    CircularMultiRecordOptions,
    ColorOptions,
    LinearDiagramOptions,
    LinearMultiRecordOptions,
)
from gbdraw.api.request_render import (
    CircularBatchRenderResult,
    CurrentRequestArtifacts,
    PreparedCircularBatchRequest,
    CircularRequestPlan,
    LinearRequestPlan,
    PreparedDiagramRequest,
    build_request_plan_diagram,
    build_prepared_interactive_context,
    build_request_diagram,
    normalize_request_records,
    plan_request,
    plan_circular_request,
    plan_circular_batch_request,
    plan_linear_request,
    render_request,
)
from gbdraw.api.requests import (
    CircularBatchOutputPolicy,
    CircularBatchRequest,
    CircularDiagramRequest,
    GffFastaInputSource,
    InMemoryRecordSource,
    LinearDiagramRequest,
    RecordCardinality,
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
    metadata = diagram_module.LinearDiagramMetadata(
        protein_comparisons=(
            DataFrame(columns=request_render_module.COMPARISON_COLUMNS),
        )
    )
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
        cache=None,
        extraction=SimpleNamespace(identity_manifest=manifest),
        nucleotide_entries=(),
        passthrough_derived_entries=(),
        source_mode="collinear",
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
        metadata=metadata,
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
        context.metadata,
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
        context.metadata,
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
        context.metadata,
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
    metadata = diagram_module.LinearDiagramMetadata(
        protein_comparisons=(
            DataFrame(columns=request_render_module.COMPARISON_COLUMNS),
        )
    )
    artifacts = request_render_module._PreparedLinearArtifacts(
        cache=None,
        extraction=extraction,
        nucleotide_entries=(),
        passthrough_derived_entries=(),
        source_mode=mode,
    )

    (entry,) = request_render_module._build_current_derived_entries(
        metadata,
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
    artifacts = CurrentRequestArtifacts(
        losat_derived_cache_entries=(entry,),
        protein_identity_manifest=manifest.to_dict(),
    )
    assert artifacts.losat_derived_cache_entries == (entry,)


def test_current_derived_artifacts_require_resolved_manifest_handles() -> None:
    record = _protein_record("source", "protein")
    extraction = extract_web_stable_cds_proteins(
        (record,),
        record_instance_keys=("record-1",),
        record_source_ids=(record.id,),
    )
    manifest = extraction.identity_manifest
    assert manifest is not None
    unresolved_handle = "h_" + ("a" * 26)
    assert all(
        unresolved_handle not in binding["runtimeIds"].values()
        for binding in manifest.record_instances.values()
    )
    entry = {
        "schema": 3,
        "kind": "derived-losatp-payload",
        "idEncoding": "runtime-handle-v1",
        "key": "unresolved-derived",
        "mode": "orthogroup",
        "payload": {
            "orthogroups": [
                {
                    "proteinIds": [unresolved_handle],
                }
            ]
        },
    }

    with pytest.raises(
        ValidationError,
        match="require protein_identity_manifest",
    ):
        CurrentRequestArtifacts(
            losat_derived_cache_entries=(entry,),
        )
    with pytest.raises(
        ValidationError,
        match="unresolved protein references",
    ):
        validate_current_session_artifacts(
            {
                "losatCache": {"entries": []},
                "losatDerivedCache": {"entries": [entry]},
            }
        )
    with pytest.raises(
        ValidationError,
        match="unresolved protein references",
    ):
        CurrentRequestArtifacts(
            losat_derived_cache_entries=(entry,),
            protein_identity_manifest=manifest.to_dict(),
        )
    with pytest.raises(
        ValidationError,
        match="unresolved protein references",
    ):
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


def test_fresh_request_render_has_no_session_compatibility_path(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import gbdraw.api.session_compat as session_compat

    def fail_if_called(*_args, **_kwargs):
        raise AssertionError("fresh rendering entered session compatibility")

    monkeypatch.setattr(
        session_compat,
        "_read_session_artifact_source",
        fail_if_called,
    )
    request = CircularDiagramRequest(records=(_memory_input("fresh"),))

    prepared = build_request_diagram(
        request,
        artifacts=CurrentRequestArtifacts(),
    )

    assert prepared.mode == "circular"
    source = Path(request_render_module.__file__).read_text(encoding="utf-8").lower()
    assert "legacy" not in source
    assert "session_artifacts" not in source
    assert "session_compat" not in source
    assert "protein_id_map" not in source
    assert "session_artifacts" not in inspect.signature(build_request_diagram).parameters
    assert "session_artifacts" not in inspect.signature(render_request).parameters
    fresh_adapter_paths = (
        Path(request_render_module.__file__).parents[1] / "circular.py",
        Path(request_render_module.__file__).parents[1] / "interface.py",
    )
    for adapter_path in fresh_adapter_paths:
        adapter_source = adapter_path.read_text(encoding="utf-8")
        assert "legacy_protein_raw_candidates" not in adapter_source
        assert "legacy_protein_derived_evidence" not in adapter_source
        assert "session_artifacts" not in adapter_source
    linear_source = (
        Path(request_render_module.__file__).parents[1] / "linear.py"
    ).read_text(encoding="utf-8")
    assert 'getattr(render_result, "legacy_protein_raw_candidates"' not in linear_source
    assert (
        'getattr(render_result, "legacy_protein_derived_evidence"'
        not in linear_source
    )
    assert "session_artifacts" not in linear_source


def test_current_request_artifacts_reject_unsupported_schema() -> None:
    from gbdraw.api import CurrentRequestArtifacts as ExportedCurrentRequestArtifacts

    assert ExportedCurrentRequestArtifacts is CurrentRequestArtifacts
    with pytest.raises(ValidationError, match="Unsupported current LOSAT artifact"):
        CurrentRequestArtifacts(
            losat_cache_entries=(
                {
                    "schema": 99,
                    "kind": "raw-losat",
                    "key": "future",
                    "text": "",
                },
            )
        )
    with pytest.raises(
        ValidationError,
        match="Unsupported current derived LOSATP artifact",
    ):
        CurrentRequestArtifacts(
            losat_derived_cache_entries=(
                {
                    "schema": 99,
                    "kind": "derived-losatp-payload",
                    "key": "future",
                    "payload": {},
                },
            )
        )


def test_build_request_diagram_rejects_falsey_untyped_artifacts() -> None:
    request = CircularDiagramRequest(records=(_memory_input("fresh"),))

    with pytest.raises(ValidationError, match="CurrentRequestArtifacts"):
        build_request_diagram(request, artifacts={})  # type: ignore[arg-type]


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


def test_linear_record_first_cardinality_is_explicit_with_comparisons() -> None:
    fixture_dir = Path(__file__).parents[1] / "examples" / "gff3_lambda"
    request = LinearDiagramRequest(
        records=(
            RecordInput(
                source=GffFastaInputSource(
                    fixture_dir / "lambda_two_contigs.gff3",
                    fixture_dir / "lambda_two_contigs.fna",
                ),
                cardinality=RecordCardinality.FIRST,
            ),
        ),
        options=LinearDiagramOptions(blast_files=("comparison.tsv",)),
    )

    assert normalize_request_records(request)[0].id == "lambda_left"


@pytest.mark.linear
def test_prepared_request_resolves_feature_inputs_once_for_gff_and_metadata(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fixture_dir = Path(__file__).parents[1] / "examples" / "gff3_lambda"
    color_file = tmp_path / "colors.tsv"
    color_file.write_text(
        "CDS\tproduct\t.*\t#123456\tCoding sequence\n",
        encoding="utf-8",
    )
    visibility_file = tmp_path / "visibility.tsv"
    visibility_file.write_text(
        "*\tCDS\tproduct\t.*\tshow\n",
        encoding="utf-8",
    )
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
        options=LinearDiagramOptions(
            colors=ColorOptions(color_table_file=str(color_file)),
            feature_visibility_table_file=str(visibility_file),
            selected_features_set=("CDS",),
        ),
    )
    calls = {
        "color_reader": 0,
        "visibility_reader": 0,
        "default_colors": 0,
        "candidate_resolver": 0,
        "visibility_compiler": 0,
        "color_compiler": 0,
    }
    real_color_reader = request_render_module.read_color_table
    real_visibility_reader = request_render_module.read_feature_visibility_file
    real_default_colors = request_render_module.load_default_colors
    real_candidate_resolver = request_render_module.resolve_candidate_feature_types
    real_visibility_compiler = prepared_module.compile_feature_visibility_rules
    real_color_compiler = prepared_module.preprocess_color_tables

    def count(name, function):
        def wrapped(*args, **kwargs):
            calls[name] += 1
            return function(*args, **kwargs)

        return wrapped

    def unexpected(*_args, **_kwargs):
        pytest.fail("resolved feature input was recomputed downstream")

    monkeypatch.setattr(
        request_render_module,
        "read_color_table",
        count("color_reader", real_color_reader),
    )
    monkeypatch.setattr(
        request_render_module,
        "read_feature_visibility_file",
        count("visibility_reader", real_visibility_reader),
    )
    monkeypatch.setattr(
        request_render_module,
        "load_default_colors",
        count("default_colors", real_default_colors),
    )
    monkeypatch.setattr(
        request_render_module,
        "resolve_candidate_feature_types",
        count("candidate_resolver", real_candidate_resolver),
    )
    monkeypatch.setattr(
        prepared_module,
        "compile_feature_visibility_rules",
        count("visibility_compiler", real_visibility_compiler),
    )
    monkeypatch.setattr(
        prepared_module,
        "preprocess_color_tables",
        count("color_compiler", real_color_compiler),
    )
    monkeypatch.setattr(diagram_module, "read_color_table", unexpected)
    monkeypatch.setattr(diagram_module, "read_feature_visibility_file", unexpected)
    monkeypatch.setattr(diagram_module, "load_default_colors", unexpected)
    monkeypatch.setattr(
        feature_configurator_module,
        "compile_feature_visibility_rules",
        unexpected,
    )
    monkeypatch.setattr(
        feature_configurator_module,
        "preprocess_color_tables",
        unexpected,
    )
    monkeypatch.setattr(
        interactive_context_module,
        "compile_feature_visibility_rules",
        unexpected,
    )
    monkeypatch.setattr(
        interactive_context_module,
        "preprocess_color_tables",
        unexpected,
    )
    monkeypatch.setattr(
        api_io_module,
        "resolve_candidate_feature_types",
        unexpected,
    )

    prepared = build_request_diagram(request)
    assert isinstance(prepared, PreparedDiagramRequest)
    context = build_prepared_interactive_context(prepared)

    assert context is not None
    assert calls == {
        "color_reader": 1,
        "visibility_reader": 1,
        "default_colors": 1,
        "candidate_resolver": 1,
        "visibility_compiler": 1,
        "color_compiler": 1,
    }


def test_option_resolution_clears_shadowed_file_inputs_without_reading(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    direct = DataFrame()

    def unexpected(*_args, **_kwargs):
        pytest.fail("an in-memory table must take precedence over its file fallback")

    for name in (
        "read_color_table",
        "load_default_colors",
        "read_feature_visibility_file",
        "read_filter_list_file",
        "read_qualifier_priority_file",
        "read_label_override_file",
    ):
        monkeypatch.setattr(request_render_module, name, unexpected)

    resolved = request_render_module._resolve_request_option_tables(
        LinearDiagramOptions(
            colors=ColorOptions(
                color_table=direct,
                color_table_file="unused-colors.tsv",
                default_colors=direct,
                default_colors_file="unused-defaults.tsv",
            ),
            feature_visibility_table=direct,
            feature_visibility_table_file="unused-visibility.tsv",
        ),
        mode="linear",
        load_comparison_colors=False,
    )

    assert resolved.colors is not None
    assert resolved.colors.color_table is direct
    assert resolved.colors.color_table_file is None
    assert resolved.colors.default_colors is direct
    assert resolved.colors.default_colors_file is None
    assert resolved.feature_visibility_table is direct
    assert resolved.feature_visibility_table_file is None


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
        "_normalize_request_records",
        lambda request, _inputs: (
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

    monkeypatch.setattr(
        request_render_module,
        "_normalize_request_records",
        lambda _request, _inputs: records,
    )

    def fake_build(loaded_records, *, options, layout, **_kwargs):
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


def test_build_linear_request_uses_typed_result_builder(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records = (_seqrecord("a"), _seqrecord("b"))
    request = LinearDiagramRequest(records=(_memory_input("a"), _memory_input("b")))
    drawing = Drawing("out.svg")
    captured: dict[str, object] = {}

    monkeypatch.setattr(
        request_render_module,
        "_normalize_request_records",
        lambda _request, _inputs: records,
    )

    def fake_build(loaded_records, *, options, **_kwargs):
        captured["records"] = loaded_records
        captured["options"] = options
        return drawing

    monkeypatch.setattr(
        request_render_module,
        "build_linear_diagram_result",
        fake_build,
    )

    prepared = build_request_diagram(request)

    assert prepared.mode == "linear"
    assert prepared.drawing is drawing
    assert captured == {"records": records, "options": request.options}


def test_resolved_plan_can_be_built_without_planning_again(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    request = CircularDiagramRequest(records=(_memory_input("planned"),))
    plan = plan_request(request)

    def unexpected_replan(*_args, **_kwargs):
        raise AssertionError("resolved plan was planned again")

    monkeypatch.setattr(
        request_render_module,
        "plan_circular_request",
        unexpected_replan,
    )

    prepared = build_request_plan_diagram(plan)

    assert isinstance(prepared, PreparedDiagramRequest)
    assert prepared.request is plan.request
    assert prepared.records is plan.records


def test_computed_orthogroups_flow_through_typed_metadata_to_interactive_context(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records = (
        _protein_record("source-a", "protein-a"),
        _protein_record("source-b", "protein-b"),
    )
    request = LinearDiagramRequest(
        records=tuple(
            RecordInput(
                source=InMemoryRecordSource(record),
                record_key=f"record-{index + 1}",
            )
            for index, record in enumerate(records)
        ),
        options=LinearDiagramOptions(protein_blastp_mode="orthogroup"),
    )
    computed: dict[str, OrthogroupResult] = {}

    def fake_analysis(
        _records,
        *,
        protein_extraction,
        **_kwargs,
    ) -> ProteinBlastpResult:
        members = [
            OrthogroupMember(
                orthogroup_id="OG0001",
                protein_id=protein.protein_id,
                record_index=protein.record_index,
                feature_index=protein.feature_index,
                record_id=protein.record_id,
                label=protein.label,
                start=protein.start,
                end=protein.end,
                strand=protein.strand,
                feature_svg_id=protein.feature_svg_id,
                source_protein_id=protein.source_protein_id,
            )
            for record_proteins in protein_extraction.proteins_by_record
            for protein in record_proteins
        ]
        orthogroups = OrthogroupResult(
            orthogroups={"OG0001": members},
            member_by_protein_id={
                member.protein_id: member
                for member in members
            },
        )
        computed["value"] = orthogroups
        return ProteinBlastpResult(
            comparisons=[
                DataFrame(columns=request_render_module.COMPARISON_COLUMNS)
            ],
            orthogroups=orthogroups,
        )

    monkeypatch.setattr(
        diagram_module,
        "build_rbh_orthogroup_protein_blastp_comparisons",
        fake_analysis,
    )

    prepared = build_request_diagram(request)
    assert isinstance(prepared, PreparedDiagramRequest)
    context = build_prepared_interactive_context(prepared)

    assert prepared.linear_metadata is not None
    assert prepared.linear_metadata.orthogroups is computed["value"]
    assert [group["id"] for group in context.orthogroups] == ["OG0001"]
    assert {
        member["recordIndex"]
        for member in context.orthogroups[0]["members"]
    } == {0, 1}
    assert not hasattr(prepared.drawing, "_gbdraw_orthogroups")


def test_circular_batch_items_share_prepared_feature_inputs(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    request = CircularBatchRequest(
        records=(_memory_input("batch-a"), _memory_input("batch-b")),
        outputs=(
            RenderOutputRequest(output_prefix="batch-a"),
            RenderOutputRequest(output_prefix="batch-b"),
        ),
    )
    received = []

    def fake_build(_record, **kwargs):
        received.append(kwargs["_resolved_feature_inputs"])
        return Drawing("out.svg")

    monkeypatch.setattr(
        request_render_module,
        "build_circular_diagram",
        fake_build,
    )

    plan = plan_circular_batch_request(request)
    item_plans = plan.item_plans()
    prepared = build_request_diagram(request)

    assert all(item.inputs is plan.inputs for item in item_plans)
    assert isinstance(prepared, PreparedCircularBatchRequest)
    assert all(item.inputs is prepared.inputs for item in prepared.items)
    assert received == [prepared.inputs.features, prepared.inputs.features]


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

    monkeypatch.setattr(
        request_render_module,
        "build_request_plan_diagram",
        lambda _: prepared,
    )
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


def test_render_request_plain_svg_retains_catalog_context_only_when_opted_in(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    request = LinearDiagramRequest(
        records=(_memory_input("catalog-record"),),
        output=RenderOutputRequest(
            output_prefix="catalog-diagram",
            output_directory=tmp_path,
            formats=("svg",),
        ),
    )
    prepared = PreparedDiagramRequest(
        mode="linear",
        request=request,
        records=(_seqrecord("catalog-record"),),
        drawing=Drawing("out.svg"),
    )
    context = object()
    captured: dict[str, object] = {}
    monkeypatch.setattr(
        request_render_module,
        "build_request_plan_diagram",
        lambda _: prepared,
    )
    monkeypatch.setattr(
        request_render_module,
        "build_prepared_interactive_context",
        lambda *_args, **_kwargs: context,
    )

    def fake_save(_canvas, _formats, **kwargs):
        captured.update(kwargs)
        return [str(tmp_path / "catalog-diagram.svg")]

    monkeypatch.setattr(request_render_module, "save_figure_to", fake_save)

    result = render_request(request, include_feature_catalog=True)

    assert result.interactive_context is context
    assert captured["interactive_context"] is None


def test_render_request_plain_svg_skips_catalog_context_by_default(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    request = LinearDiagramRequest(
        records=(_memory_input("default-record"),),
        output=RenderOutputRequest(
            output_prefix="default-diagram",
            output_directory=tmp_path,
            formats=("svg",),
        ),
    )
    prepared = PreparedDiagramRequest(
        mode="linear",
        request=request,
        records=(_seqrecord("default-record"),),
        drawing=Drawing("out.svg"),
    )
    monkeypatch.setattr(
        request_render_module,
        "build_request_plan_diagram",
        lambda _: prepared,
    )
    monkeypatch.setattr(
        request_render_module,
        "build_prepared_interactive_context",
        lambda *_args, **_kwargs: pytest.fail(
            "plain SVG must not build feature context without an explicit opt-in"
        ),
    )
    monkeypatch.setattr(
        request_render_module,
        "save_figure_to",
        lambda *_args, **_kwargs: [str(tmp_path / "default-diagram.svg")],
    )

    result = render_request(request)

    assert result.interactive_context is None


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
        "build_request_plan_diagram",
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
@pytest.mark.parametrize(
    ("formats", "include_feature_catalog"),
    [
        (("interactive_svg",), False),
        (("svg",), True),
    ],
)
def test_render_request_circular_batch_loads_needed_comparison_fasta_once(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    formats: tuple[str, ...],
    include_feature_catalog: bool,
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
                formats=formats,
            )
            for prefix in ("batch-a", "batch-b")
        ),
    )

    result = render_request(
        request,
        include_feature_catalog=include_feature_catalog,
    )

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


def test_prepared_request_memoizes_comparison_fasta_records(
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
    request = CircularDiagramRequest(
        records=(_memory_input("record"),),
        options=CircularDiagramOptions(
            conservation_fasta_files=(str(comparison_fasta),),
        ),
        output=RenderOutputRequest(formats=("interactive_svg",)),
    )
    prepared = build_request_diagram(request)
    assert isinstance(prepared, PreparedDiagramRequest)

    first = request_render_module._comparison_sequence_records(prepared)
    second = request_render_module._comparison_sequence_records(prepared)

    assert first is second
    assert parse_calls == [str(comparison_fasta)]


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
        "build_request_plan_diagram",
        lambda _plan: pytest.fail(
            "batch outputs must be preflighted before building"
        ),
    )

    with pytest.raises(ValidationError, match="already exist"):
        render_request(request)

    assert not (tmp_path / "batch-a.svg").exists()
    assert existing.read_text(encoding="utf-8") == "keep"


@pytest.mark.circular
def test_render_request_circular_batch_rejects_aliased_output_paths_before_build(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    nested = tmp_path / "nested"
    request = CircularBatchRequest(
        records=(_memory_input("batch-a"), _memory_input("batch-b")),
        outputs=(
            RenderOutputRequest(
                output_prefix="diagram",
                output_directory=nested / "..",
            ),
            RenderOutputRequest(
                output_prefix="diagram",
                output_directory=tmp_path,
            ),
        ),
    )
    monkeypatch.setattr(
        request_render_module,
        "build_request_plan_diagram",
        lambda _plan: pytest.fail(
            "aliased batch outputs must be rejected before building"
        ),
    )

    with pytest.raises(ValidationError, match="duplicate file paths"):
        render_request(request)

    assert not (tmp_path / "diagram.svg").exists()


@pytest.mark.circular
def test_render_request_preflights_resolved_batch_policy_before_build(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    occupied = tmp_path / "diagram_2.svg"
    occupied.write_text("keep", encoding="utf-8")
    request = CircularBatchRequest(
        records=(_memory_input("batch-a"), _memory_input("batch-b")),
        output_policy=CircularBatchOutputPolicy(
            output_prefix=tmp_path / "diagram",
        ),
    )
    monkeypatch.setattr(
        request_render_module,
        "build_request_plan_diagram",
        lambda _plan: pytest.fail(
            "resolved batch policy outputs must be preflighted before building"
        ),
    )

    with pytest.raises(ValidationError, match="already exist"):
        render_request(request)

    assert not (tmp_path / "diagram_1.svg").exists()
    assert occupied.read_text(encoding="utf-8") == "keep"


@pytest.mark.circular
@pytest.mark.parametrize("target_kind", ("directory", "file-parent"))
def test_render_request_rejects_nonreplaceable_batch_target_before_build(
    target_kind: str,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    if target_kind == "directory":
        blocked_directory = tmp_path / "blocked.svg"
        blocked_directory.mkdir()
        blocked_output = RenderOutputRequest(
            output_prefix="blocked",
            output_directory=tmp_path,
            overwrite=True,
        )
        expected = "not replaceable files"
    else:
        blocked_parent = tmp_path / "blocked-parent"
        blocked_parent.write_text("not a directory", encoding="utf-8")
        blocked_output = RenderOutputRequest(
            output_prefix="blocked",
            output_directory=blocked_parent,
            overwrite=True,
        )
        expected = "parent path.*not directories"
    request = CircularBatchRequest(
        records=(_memory_input("batch-a"), _memory_input("batch-b")),
        outputs=(
            RenderOutputRequest(
                output_prefix="first",
                output_directory=tmp_path,
                overwrite=True,
            ),
            blocked_output,
        ),
    )
    monkeypatch.setattr(
        request_render_module,
        "build_request_plan_diagram",
        lambda _plan: pytest.fail(
            "nonreplaceable batch targets must be rejected before building"
        ),
    )

    with pytest.raises(ValidationError, match=expected):
        render_request(request)

    assert not (tmp_path / "first.svg").exists()


@pytest.mark.circular
def test_render_request_rejects_implicit_record_path_through_symlink(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    outside = tmp_path / "outside"
    outside.mkdir()
    linked_parent = tmp_path / "linked"
    try:
        linked_parent.symlink_to(outside, target_is_directory=True)
    except OSError as exc:
        pytest.skip(f"directory symlinks are unavailable: {exc}")
    monkeypatch.chdir(tmp_path)
    request = CircularBatchRequest(
        records=(_memory_input("linked/escaped"),),
        output_policy=CircularBatchOutputPolicy(),
    )
    monkeypatch.setattr(
        request_render_module,
        "build_request_plan_diagram",
        lambda _plan: pytest.fail(
            "unsafe implicit output must be rejected before building"
        ),
    )

    with pytest.raises(
        ValidationError,
        match="cannot be used as an implicit output filename prefix",
    ):
        render_request(request)

    assert not (outside / "escaped.svg").exists()


@pytest.mark.parametrize("mode", ("circular", "linear"))
def test_render_request_preflights_every_single_request_format_before_build(
    mode: str,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    svg_path = tmp_path / "diagram.svg"
    svg_path.write_text("keep", encoding="utf-8")
    blocked_png = tmp_path / "diagram.png"
    blocked_png.mkdir()
    output = RenderOutputRequest(
        output_prefix="diagram",
        output_directory=tmp_path,
        formats=("png",),
        overwrite=True,
    )
    request = (
        CircularDiagramRequest(
            records=(_memory_input("record"),),
            output=output,
        )
        if mode == "circular"
        else LinearDiagramRequest(
            records=(_memory_input("record"),),
            output=output,
        )
    )
    monkeypatch.setattr(
        request_render_module,
        "build_request_plan_diagram",
        lambda _plan: pytest.fail(
            "single-request outputs must be preflighted before building"
        ),
    )

    with pytest.raises(ValidationError, match="not replaceable files"):
        render_request(request)

    assert svg_path.read_text(encoding="utf-8") == "keep"
    assert blocked_png.is_dir()
