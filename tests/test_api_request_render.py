from __future__ import annotations

import ast
from dataclasses import replace
from pathlib import Path
from types import SimpleNamespace

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from pandas import DataFrame
from svgwrite import Drawing

import gbdraw.analysis.protein_colinearity as protein_colinearity_module
import gbdraw.api.request_render as request_render_module
from gbdraw.analysis.collinearity import (
    CollinearityParameters,
    LosslessCollinearityParameters,
)
from gbdraw.analysis.protein_colinearity import (
    OrthogroupMember,
    OrthogroupResult,
    build_losat_transport_id,
    extract_web_stable_cds_proteins,
)
from gbdraw.api.options import DiagramOptions
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
from gbdraw.linear_comparison import LinearComparison
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


def _v35_manifest_for_extraction(
    extraction,
    aliases_by_instance: dict[str, str],
) -> tuple[dict[str, object], dict[str, str]]:
    current = extraction.identity_manifest
    assert current is not None
    record_instances: dict[str, object] = {}
    old_ids: dict[str, str] = {}
    for instance_key, alias in aliases_by_instance.items():
        current_binding = current.record_instances[instance_key]
        analysis_id = str(current_binding["recordAnalysisId"])
        analysis = current.record_analyses[analysis_id]
        feature_id = next(iter(current_binding["runtimeIds"]))
        old_id = build_losat_transport_id(
            record_source_id=analysis["recordSourceId"],
            record_instance_id=instance_key,
            display_alias=alias,
            feature_analysis_id=feature_id,
        )
        transport_ids = {feature_id: old_id}
        binding_payload = {
            "recordInstanceKey": instance_key,
            "recordAnalysisId": analysis_id,
            "transportIds": transport_ids,
        }
        record_instances[instance_key] = {
            "schema": 1,
            "recordAnalysisId": analysis_id,
            "bindingHash": protein_colinearity_module._identity_sha256(
                binding_payload
            ),
            "transportIds": transport_ids,
        }
        old_ids[instance_key] = old_id
    return (
        {
            "schema": 1,
            "proteinSets": current.protein_sets,
            "recordAnalyses": current.record_analyses,
            "recordInstances": record_instances,
        },
        old_ids,
    )


def test_rewrite_protein_artifact_references_rejects_key_collisions() -> None:
    with pytest.raises(ValidationError, match="duplicate key"):
        request_render_module.rewrite_protein_artifact_references(
            {"legacy-a": 1, "legacy-b": 2},
            {"legacy-a": "h_same", "legacy-b": "h_same"},
        )


def test_rewrite_protein_artifact_references_updates_compound_legacy_ids() -> None:
    legacy_id = "p_r_old_0_9_1_deadbeefdead"
    feature_id = f"f_{'a' * 64}"
    v35_id = f"record@instance|alias~{feature_id}"

    rewritten = request_render_module.rewrite_protein_artifact_references(
        {
            "query_protein_id": f"{legacy_id};{v35_id}",
            "supporting_edge": f"{legacy_id}->{v35_id}:rbh",
        },
        {
            legacy_id: "h_aaaaaaaaaaaaaaaaaaaaaaaaaa",
            v35_id: "h_bbbbbbbbbbbbbbbbbbbbbbbbbb",
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
    options: DiagramOptions | None = None,
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
        options=options or DiagramOptions(),
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
        v35_candidates=None,
        v35_derived_evidence=(),
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
    params: CollinearityParameters | LosslessCollinearityParameters | None,
) -> dict[str, object]:
    context = _derived_identity_test_context(
        options=DiagramOptions(collinearity_params=params)
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
        "max-gene-gap",
        "max-diagonal-drift",
        "max-conflicts-in-merge-gap",
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


@pytest.mark.parametrize(
    ("field_name", "changed_value"),
    (
        ("min_anchors", 2),
        ("max_gene_gap", 26),
        ("block_merge_gap", 51),
        ("singleton_merge_gap", 26),
        ("max_diagonal_drift", 26),
        ("max_conflicts_in_merge_gap", 2),
        ("max_paralog_links_per_orthogroup", 3),
        ("gap_penalty", 2.0),
        ("nearby_duplicate_window", 1),
        ("score_mode", "bitscore"),
        ("constant_anchor_score", 51.0),
        ("min_block_score", 50.0),
        ("block_evalue", 0.01),
    ),
)
def test_derived_identity_snapshots_every_legacy_collinearity_parameter(
    monkeypatch: pytest.MonkeyPatch,
    field_name: str,
    changed_value: object,
) -> None:
    monkeypatch.setattr(
        request_render_module,
        "is_protein_losat_cache_entry",
        lambda _entry: True,
    )
    baseline_params = CollinearityParameters()
    baseline = _derived_entry_for_params(baseline_params)
    unchanged = _derived_entry_for_params(CollinearityParameters())
    changed = _derived_entry_for_params(
        replace(baseline_params, **{field_name: changed_value})
    )

    parameter_identity = baseline["payload"]["identity"]["collinear"][
        "parameterIdentity"
    ]
    assert parameter_identity["model"] == "legacy"
    assert set(parameter_identity["parameters"]) == {
        "minAnchors",
        "maxGeneGap",
        "blockMergeGap",
        "singletonMergeGap",
        "maxDiagonalDrift",
        "maxConflictsInMergeGap",
        "maxParalogLinksPerOrthogroup",
        "gapPenalty",
        "nearbyDuplicateWindow",
        "scoreMode",
        "constantAnchorScore",
        "minBlockScore",
        "blockEvalue",
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
        options=DiagramOptions(protein_blastp_mode=mode),
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
        v35_candidates=None,
        v35_derived_evidence=(),
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


def test_api_cli_render_migrates_manifest_only_v35_artifact_references(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records = (
        _protein_record("source-left", "current-left"),
        _protein_record("source-right", "current-right"),
    )
    extraction = extract_web_stable_cds_proteins(
        records,
        record_instance_keys=("left", "right"),
        record_source_ids=tuple(record.id for record in records),
    )
    source_manifest, old_ids = _v35_manifest_for_extraction(
        extraction,
        {"left": "old-left", "right": "old-right"},
    )
    comparison = DataFrame(
        [
            {
                "query": old_ids["left"],
                "subject": old_ids["right"],
                "identity": 95.0,
                "alignment_length": 2,
                "mismatches": 0,
                "gap_opens": 0,
                "qstart": 1,
                "qend": 2,
                "sstart": 1,
                "send": 2,
                "evalue": 1e-20,
                "bitscore": 80.0,
                "query_protein_id": old_ids["left"],
                "subject_protein_id": old_ids["right"],
            }
        ]
    )
    left_member = OrthogroupMember(
        orthogroup_id="OG1",
        protein_id=old_ids["left"],
        record_index=0,
        feature_index=0,
        record_id=records[0].id,
        label="Left",
        start=0,
        end=9,
        strand=1,
        feature_svg_id=None,
        source_protein_id=old_ids["left"],
    )
    right_member = OrthogroupMember(
        orthogroup_id="OG1",
        protein_id=old_ids["right"],
        record_index=1,
        feature_index=0,
        record_id=records[1].id,
        label="Right",
        start=0,
        end=9,
        strand=1,
        feature_svg_id=None,
        source_protein_id=old_ids["right"],
    )
    orthogroups = OrthogroupResult(
        orthogroups={"OG1": [left_member, right_member]},
        member_by_protein_id={
            old_ids["left"]: left_member,
            old_ids["right"]: right_member,
        },
    )
    request = LinearDiagramRequest(
        records=(
            RecordInput(
                source=InMemoryRecordSource(records[0]),
                record_key="left",
            ),
            RecordInput(
                source=InMemoryRecordSource(records[1]),
                record_key="right",
            ),
        ),
        options=DiagramOptions(
            linear_comparisons=(LinearComparison(0, 1, comparison),),
            orthogroups=orthogroups,
        ),
    )
    artifacts = {
        "legacyArtifacts": {
            "proteinRawV35Candidates": {
                "schema": 1,
                "sourceManifest": source_manifest,
                "entries": [],
            }
        }
    }
    def fail_losat(*_args, **_kwargs):
        raise AssertionError("LOSAT must not run for precomputed artifacts")

    monkeypatch.setattr(
        protein_colinearity_module,
        "run_losatp_blastp",
        fail_losat,
    )
    monkeypatch.setattr(
        request_render_module,
        "save_figure_to",
        lambda *_args, **_kwargs: [],
    )

    rendered = render_request(request, session_artifacts=artifacts)

    assert rendered.output_paths == ()
    assert rendered.losat_cache_entries == ()
    assert rendered.legacy_protein_raw_v35_candidates is None
    assert set(rendered.protein_id_map or {}) == set(old_ids.values())
    assert all(
        runtime_handle.startswith("h_")
        for runtime_handle in (rendered.protein_id_map or {}).values()
    )
    rewritten_options = rendered.request.options
    rewritten_comparison = rewritten_options.linear_comparisons[0].matches.iloc[0]
    assert rewritten_comparison["query"] == rendered.protein_id_map[old_ids["left"]]
    assert (
        rewritten_comparison["subject"]
        == rendered.protein_id_map[old_ids["right"]]
    )
    rewritten_orthogroups = rewritten_options.orthogroups
    assert set(rewritten_orthogroups.member_by_protein_id) == set(
        rendered.protein_id_map.values()
    )
    assert len(rendered.losat_derived_cache_entries) == 1
    serialized_derived = str(rendered.losat_derived_cache_entries)
    assert old_ids["left"] not in serialized_derived
    assert old_ids["right"] not in serialized_derived


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
    assert result.interactive_context is context
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
