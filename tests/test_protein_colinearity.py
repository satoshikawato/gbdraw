from __future__ import annotations

import copy
from contextlib import ExitStack
import json
import os
from pathlib import Path
import shutil
import subprocess
from types import SimpleNamespace

import pandas as pd
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from svgwrite import Drawing

import gbdraw.api.diagram as api_diagram_module
import gbdraw.api.request_render as request_render_module
import gbdraw.analysis.protein_colinearity as protein_colinearity_module
import gbdraw.linear as linear_cli_module
import gbdraw.losat_setup as losat_setup_module
from gbdraw.api.config import apply_config_overrides
from gbdraw.api.requests import LinearDiagramRequest
from gbdraw.analysis.protein_colinearity import (
    PROTEIN_LOSAT_CACHE_SCHEMA,
    LosatpCacheManager,
    build_legacy_protein_reference_map,
    build_protein_export_id_map,
    build_protein_losat_cache_key,
    build_protein_losat_pair_identity,
    build_protein_runtime_handle,
    build_web_losat_cache_key,
    build_pairwise_protein_blastp_comparisons,
    build_rbh_orthogroup_protein_blastp_comparisons,
    convert_pair_protein_hits_to_genomic_links,
    convert_protein_hits_to_genomic_links,
    extract_cds_proteins,
    extract_protein_identity_manifest,
    extract_web_stable_cds_proteins,
    filter_protein_hits_by_thresholds,
    hydrate_protein_losat_tsv,
    parse_losatp_outfmt6,
    percent_encode_losat_transport_field,
    proteins_to_fasta,
    promote_legacy_protein_raw_cache_entries,
    select_best_hits_per_query,
    select_reciprocal_best_hit_edges,
    select_reciprocal_best_hits,
    select_rbh_orthogroup_edges_from_directional_hits,
    select_top_hits_per_query,
    validate_protein_identity_manifest,
    validate_protein_raw_entry_references,
)
from gbdraw.api.diagram import assemble_linear_diagram_from_records
from gbdraw.api.options import LinearDiagramOptions
from gbdraw.diagrams.linear.orthogroup_alignment import (
    calculate_orthogroup_alignment_canvas_adjustment,
    calculate_orthogroup_alignment_canvas_extents,
    calculate_orthogroup_alignment_offsets,
)
from gbdraw.exceptions import ValidationError
from gbdraw.io.comparisons import COMPARISON_COLUMNS
from gbdraw.io.record_select import reverse_records
from gbdraw.render.groups.linear.pairwise_match import PairWiseMatchGroup


def _record(
    record_id: str,
    sequence: str = "ATGAAATAG" * 20,
    features: list[SeqFeature] | None = None,
) -> SeqRecord:
    record = SeqRecord(Seq(sequence), id=record_id)
    record.features = list(features or [])
    return record


def _orthogroup_alignment_canvas_config() -> object:
    return type(
        "CanvasConfig",
        (),
        {
            "normalize_length": False,
            "align_center": False,
            "alignment_width": 1000.0,
            "longest_genome": 1000,
        },
    )()


def _cds(
    start: int,
    end: int,
    *,
    strand: int = 1,
    qualifiers: dict[str, list[str]] | None = None,
) -> SeqFeature:
    return SeqFeature(
        FeatureLocation(start, end, strand=strand),
        type="CDS",
        qualifiers=qualifiers or {"translation": ["MK*"]},
    )


def _hit_row(
    query: str,
    subject: str,
    *,
    identity: float = 90.0,
    alignment_length: int = 100,
    qstart: object = 1,
    qend: object = 100,
    sstart: object = 1,
    send: object = 100,
    evalue: float = 1e-30,
    bitscore: float = 200.0,
) -> dict[str, object]:
    return {
        "query": query,
        "subject": subject,
        "identity": identity,
        "alignment_length": alignment_length,
        "mismatches": 0,
        "gap_opens": 0,
        "qstart": qstart,
        "qend": qend,
        "sstart": sstart,
        "send": send,
        "evalue": evalue,
        "bitscore": bitscore,
    }


def _web_protein_entry(
    protein_id: str,
    *,
    record_index: int,
    record_id: str,
    feature_index: int = 0,
    start: int = 0,
    end: int = 90,
    strand: int = 1,
    gene: str | None = None,
    product: str | None = None,
    note: str | None = None,
) -> dict[str, object]:
    return {
        "protein_id": protein_id,
        "record_index": record_index,
        "feature_index": feature_index,
        "record_id": record_id,
        "start": start,
        "end": end,
        "strand": strand,
        "label": protein_id,
        "protein_length": 30,
        "feature_svg_id": f"feature_{protein_id}",
        "feature_type": "CDS",
        "feature_hash_start": start,
        "feature_hash_end": end,
        "feature_hash_strand": strand,
        "gene": gene,
        "product": product,
        "note": note,
    }


def _load_web_helper_namespace() -> dict[str, object]:
    helpers_js = Path("gbdraw/web/js/app/python-helpers.js").read_text(encoding="utf-8")
    helper_source = helpers_js.split("`", 1)[1].rsplit("`", 1)[0]
    namespace: dict[str, object] = {}
    exec(helper_source, namespace)
    return namespace


def _current_protein_manifest_payload_for_validation() -> dict[str, object]:
    extraction = extract_protein_identity_manifest(
        [
            _record(
                "record_a",
                features=[
                    _cds(0, 9, qualifiers={"translation": ["MKT"]}),
                    _cds(9, 18, qualifiers={"translation": ["MGG"]}),
                ],
            )
        ],
        record_instance_keys=("row-1",),
    )
    assert extraction.identity_manifest is not None
    return extraction.identity_manifest.to_dict()


@pytest.mark.linear
def test_protein_identity_golden_hash_runtime_handle_and_export_encoding() -> None:
    feature_id = protein_colinearity_module.canonical_feature_analysis_id(
        feature_type="CDS",
        location_operator="join",
        location_parts=[(7438, 8458, 1)],
        strand=1,
        same_location_ordinal=1,
    )

    assert feature_id == "f_c038a3178fe7cfc89c61309ec0bdaf81e01cf84b1a647db670688c7f1eece649"
    assert (
        percent_encode_losat_transport_field("A z@|~%\n雪e\u0301")
        == "A%20z%40%7C%7E%25%0A%E9%9B%AA%C3%A9"
    )
    assert (
        build_protein_runtime_handle(
            feature_analysis_id=feature_id,
            record_instance_key="record-1",
        )
        == "h_xqvkmqu3ozwevaji27xbnjtmau"
    )


@pytest.mark.linear
def test_current_protein_manifest_validator_accepts_extracted_schema2_authority() -> None:
    payload = _current_protein_manifest_payload_for_validation()

    validated = validate_protein_identity_manifest(payload)

    assert validated.to_dict() == payload


@pytest.mark.linear
@pytest.mark.parametrize(
    ("corruption", "message"),
    [
        ("top_level_schema", "must use schema 2"),
        ("protein_set_hash", "Protein set .* hash does not match"),
        ("record_analysis_hash", "Record analysis .* hash does not match"),
        ("runtime_membership", "runtime map does not match its protein set"),
        ("display_membership", "metadata does not match its protein set"),
        ("runtime_binding_hash", "runtime binding hash does not match"),
        ("display_binding_hash", "display binding hash does not match"),
        ("runtime_handle_uniqueness", "deterministic and globally unique"),
    ],
)
def test_current_protein_manifest_validator_rejects_authority_corruption(
    corruption: str,
    message: str,
) -> None:
    payload = copy.deepcopy(_current_protein_manifest_payload_for_validation())
    protein_sets = payload["proteinSets"]
    record_analyses = payload["recordAnalyses"]
    record_instances = payload["recordInstances"]
    assert isinstance(protein_sets, dict)
    assert isinstance(record_analyses, dict)
    assert isinstance(record_instances, dict)
    protein_set = next(iter(protein_sets.values()))
    record_analysis = next(iter(record_analyses.values()))
    binding = record_instances["row-1"]
    assert isinstance(protein_set, dict)
    assert isinstance(record_analysis, dict)
    assert isinstance(binding, dict)
    proteins = protein_set["proteins"]
    runtime_ids = binding["runtimeIds"]
    feature_metadata = binding["featureMetadata"]
    assert isinstance(proteins, list)
    assert isinstance(runtime_ids, dict)
    assert isinstance(feature_metadata, dict)

    if corruption == "top_level_schema":
        payload["schema"] = 1
    elif corruption == "protein_set_hash":
        proteins[0]["aaSha256"] = "0" * 64
    elif corruption == "record_analysis_hash":
        record_analysis["recordSourceId"] = "tampered"
    elif corruption == "runtime_membership":
        runtime_ids.pop(next(iter(runtime_ids)))
    elif corruption == "display_membership":
        feature_metadata.pop(next(iter(feature_metadata)))
    elif corruption == "runtime_binding_hash":
        binding["runtimeBindingHash"] = f"sha256:{'0' * 64}"
    elif corruption == "display_binding_hash":
        binding["displayBindingHash"] = f"sha256:{'0' * 64}"
    elif corruption == "runtime_handle_uniqueness":
        first_feature, second_feature = tuple(runtime_ids)[:2]
        runtime_ids[second_feature] = runtime_ids[first_feature]
    else:
        raise AssertionError(f"Unhandled corruption case: {corruption}")

    with pytest.raises(ValidationError, match=message):
        validate_protein_identity_manifest(payload)


@pytest.mark.linear
@pytest.mark.parametrize(
    "forbidden_key",
    [
        "viewFeatureSvgId",
        "view_feature_svg_id",
        "viewFeatureHashParts",
        "view_feature_hash_parts",
        "renderedFeatureSvgId",
        "rendered_feature_svg_id",
        "renderedSvgId",
        "rendered_svg_id",
        "Rendered-Feature-SVG-ID",
        "queryViewFeatureSvgId",
        "processed_view_feature_svg_id",
        "subjectViewFeatureHashParts",
        "queryRenderedFeatureSvgId",
        "processed_rendered_svg_id",
    ],
)
def test_protein_manifest_rejects_nested_rendered_view_identity(
    forbidden_key: str,
) -> None:
    payload = copy.deepcopy(_current_protein_manifest_payload_for_validation())
    record_analyses = payload["recordAnalyses"]
    record_instances = payload["recordInstances"]
    assert isinstance(record_analyses, dict)
    assert isinstance(record_instances, dict)
    binding = record_instances["row-1"]
    assert isinstance(binding, dict)
    feature_metadata = binding["featureMetadata"]
    assert isinstance(feature_metadata, dict)
    feature_id = next(iter(feature_metadata))
    metadata = feature_metadata[feature_id]
    assert isinstance(metadata, dict)
    metadata["extension"] = {forbidden_key: "presentation-only"}

    analysis_id = str(binding["recordAnalysisId"])
    record_analysis = record_analyses[analysis_id]
    assert isinstance(record_analysis, dict)
    display_payload = protein_colinearity_module._display_binding_payload(
        record_analysis_id=analysis_id,
        record_source_id=str(record_analysis["recordSourceId"]),
        record_instance_key="row-1",
        feature_metadata=feature_metadata,
    )
    binding["displayBindingHash"] = protein_colinearity_module._identity_sha256(
        display_payload
    )

    with pytest.raises(ValidationError, match="rendered-view identity"):
        validate_protein_identity_manifest(payload)


@pytest.mark.linear
def test_protein_manifest_separates_raw_and_derived_invalidation() -> None:
    def make_record(*, protein_id: str, product: str) -> SeqRecord:
        record = _record(
            "record_a",
            features=[
                _cds(
                    0,
                    9,
                    qualifiers={
                        "translation": ["MKT*"],
                        "protein_id": [protein_id],
                        "product": [product],
                    },
                )
            ],
        )
        record.annotations["upload_filename"] = "ignored.gb"
        record.annotations["lastModified"] = 123
        return record

    baseline = extract_protein_identity_manifest(
        [make_record(protein_id="protein-1", product="old product")],
        record_instance_keys=("row-1",),
    )
    annotation_only = extract_protein_identity_manifest(
        [make_record(protein_id="protein-1", product="new product")],
        record_instance_keys=("row-1",),
    )
    alias_changed = extract_protein_identity_manifest(
        [make_record(protein_id="renamed protein", product="old product")],
        record_instance_keys=("row-1",),
    )

    baseline_protein = baseline.proteins_by_record[0][0]
    assert baseline_protein.feature_analysis_id == annotation_only.proteins_by_record[0][0].feature_analysis_id
    assert baseline.protein_set_hashes == annotation_only.protein_set_hashes
    assert baseline.runtime_binding_hashes == annotation_only.runtime_binding_hashes
    assert baseline.display_binding_hashes != annotation_only.display_binding_hashes
    assert baseline_protein.feature_analysis_id == alias_changed.proteins_by_record[0][0].feature_analysis_id
    assert baseline.protein_set_hashes == alias_changed.protein_set_hashes
    assert baseline.runtime_binding_hashes == alias_changed.runtime_binding_hashes
    assert baseline.display_binding_hashes != alias_changed.display_binding_hashes
    assert baseline_protein.protein_id == alias_changed.proteins_by_record[0][0].protein_id


@pytest.mark.linear
def test_protein_raw_identity_is_invariant_to_display_reverse_complement() -> None:
    record = _record(
        "record_a",
        sequence="ATGAAATAG" * 4,
        features=[_cds(3, 18, strand=1, qualifiers={"translation": ["MKT"]})],
    )
    record.annotations["gbdraw_coord_base"] = 1
    record.annotations["gbdraw_coord_step"] = 1
    reversed_record = reverse_records([record], True)[0]

    source = extract_protein_identity_manifest(
        [record],
        record_instance_keys=("row",),
    )
    reversed_view = extract_protein_identity_manifest(
        [reversed_record],
        record_instance_keys=("row",),
    )

    assert (
        source.proteins_by_record[0][0].feature_analysis_id
        == reversed_view.proteins_by_record[0][0].feature_analysis_id
    )
    assert source.protein_set_hashes == reversed_view.protein_set_hashes
    assert source.runtime_binding_hashes == reversed_view.runtime_binding_hashes


@pytest.mark.linear
def test_protein_manifest_compound_location_and_same_location_ordinals_are_stable() -> None:
    location = CompoundLocation(
        [FeatureLocation(0, 6, strand=1), FeatureLocation(12, 18, strand=1)],
        operator="join",
    )
    features = [
        SeqFeature(
            location,
            type="CDS",
            qualifiers={"translation": [sequence], "protein_id": ["duplicate"]},
        )
        for sequence in ("MKK", "MQQ")
    ]
    first = extract_protein_identity_manifest(
        [_record("compound", features=features)],
        record_instance_keys=("row",),
    )
    reordered = extract_protein_identity_manifest(
        [_record("compound", features=list(reversed(features)))],
        record_instance_keys=("row",),
    )

    by_sequence = {protein.sequence: protein for protein in first.proteins_by_record[0]}
    reordered_by_sequence = {
        protein.sequence: protein for protein in reordered.proteins_by_record[0]
    }
    assert {protein.same_location_ordinal for protein in by_sequence.values()} == {1, 2}
    assert all(protein.location_operator == "join" for protein in by_sequence.values())
    assert all(protein.feature_hash_parts == ((0, 6, 1), (12, 18, 1)) for protein in by_sequence.values())
    assert {
        sequence: protein.feature_analysis_id for sequence, protein in by_sequence.items()
    } == {
        sequence: protein.feature_analysis_id
        for sequence, protein in reordered_by_sequence.items()
    }
    assert first.protein_set_hashes == reordered.protein_set_hashes


@pytest.mark.linear
@pytest.mark.parametrize("display_qualifier", ["protein_id", "locus_tag", "ID"])
def test_same_location_protein_identity_is_invariant_to_display_qualifier_changes(
    display_qualifier: str,
) -> None:
    def make_record(first_alias: str) -> SeqRecord:
        return _record(
            "same-location",
            features=[
                _cds(
                    0,
                    9,
                    qualifiers={
                        "translation": ["MKK"],
                        display_qualifier: [first_alias],
                    },
                ),
                _cds(
                    0,
                    9,
                    qualifiers={
                        "translation": ["MQQ"],
                        display_qualifier: ["beta"],
                    },
                ),
            ],
        )

    baseline = extract_protein_identity_manifest(
        [make_record("alpha")],
        record_instance_keys=("row",),
    )
    display_changed = extract_protein_identity_manifest(
        [make_record("zeta")],
        record_instance_keys=("row",),
    )

    baseline_by_sequence = {
        protein.sequence: protein for protein in baseline.proteins_by_record[0]
    }
    changed_by_sequence = {
        protein.sequence: protein for protein in display_changed.proteins_by_record[0]
    }
    assert {
        sequence: (protein.feature_analysis_id, protein.runtime_handle)
        for sequence, protein in baseline_by_sequence.items()
    } == {
        sequence: (protein.feature_analysis_id, protein.runtime_handle)
        for sequence, protein in changed_by_sequence.items()
    }
    assert baseline.protein_set_hashes == display_changed.protein_set_hashes
    assert baseline.runtime_binding_hashes == display_changed.runtime_binding_hashes
    assert baseline.display_binding_hashes != display_changed.display_binding_hashes


@pytest.mark.linear
def test_identical_record_instances_share_content_but_not_runtime_bindings() -> None:
    feature = _cds(
        0,
        9,
        qualifiers={"translation": ["MKT*"], "protein_id": ["same"]},
    )
    records = [
        _record("accession", features=[feature]),
        _record("accession", features=[feature]),
    ]
    extraction = extract_protein_identity_manifest(
        records,
        record_instance_keys=("row-1", "row-2"),
    )

    assert extraction.protein_set_hashes[0] == extraction.protein_set_hashes[1]
    assert extraction.record_analysis_ids[0] == extraction.record_analysis_ids[1]
    assert extraction.runtime_binding_hashes[0] != extraction.runtime_binding_hashes[1]
    assert len(extraction.identity_manifest.protein_sets) == 1
    assert len(extraction.identity_manifest.record_analyses) == 1
    assert len(extraction.protein_map) == 2
    assert extraction.proteins_by_record[0][0].protein_id != extraction.proteins_by_record[1][0].protein_id

    forward = build_protein_losat_pair_identity(
        extraction.identity_manifest,
        query_record_instance_key="row-1",
        subject_record_instance_key="row-2",
    )
    reverse = build_protein_losat_pair_identity(
        extraction.identity_manifest,
        query_record_instance_key="row-2",
        subject_record_instance_key="row-1",
    )
    assert build_protein_losat_cache_key(forward, args=[]) != build_protein_losat_cache_key(reverse, args=[])


@pytest.mark.linear
def test_protein_raw_identity_tracks_candidate_limit_args_exactly() -> None:
    extraction = extract_protein_identity_manifest(
        [
            _record("query", features=[_cds(0, 9)]),
            _record("subject", features=[_cds(3, 12)]),
        ],
        record_instance_keys=("query-row", "subject-row"),
    )
    identity = build_protein_losat_pair_identity(
        extraction.identity_manifest,
        query_record_instance_key="query-row",
        subject_record_instance_key="subject-row",
    )

    omitted = build_protein_losat_cache_key(identity, args=[])
    explicit_none = build_protein_losat_cache_key(identity, args=[])
    finite = build_protein_losat_cache_key(
        identity,
        args=["--max-target-seqs", "7"],
    )

    assert explicit_none == omitted
    assert finite != omitted
    source_scope = build_protein_losat_cache_key(identity, args=[], search_context="a" * 64)
    assert source_scope != omitted
    assert source_scope != build_protein_losat_cache_key(identity, args=[], search_context="b" * 64)
    with pytest.raises(ValidationError, match="search context"):
        build_protein_losat_cache_key(identity, args=[], search_context="invalid")
    assert finite == build_protein_losat_cache_key(
        identity,
        args=["--max-target-seqs", "7"],
    )


@pytest.mark.linear
def test_legacy_protein_artifact_references_resolve_to_runtime_handles() -> None:
    extraction = extract_protein_identity_manifest(
        [
            _record("query", features=[_cds(0, 9), _cds(9, 18)]),
            _record("subject", features=[_cds(3, 12)]),
        ],
        record_instance_keys=("left", "right"),
    )
    legacy_query = protein_colinearity_module._with_stable_web_protein_ids(
        extraction.proteins_by_record[0],
        "r_old_left",
    )
    legacy_subject = protein_colinearity_module._with_stable_web_protein_ids(
        extraction.proteins_by_record[1],
        "r_old_right",
    )

    resolved = build_legacy_protein_reference_map(
        extraction,
        [
            legacy_query[1].protein_id,
            legacy_subject[0].protein_id,
            legacy_query[0].protein_id,
        ],
    )

    assert resolved == {
        legacy_query[0].protein_id: extraction.proteins_by_record[0][0].protein_id,
        legacy_query[1].protein_id: extraction.proteins_by_record[0][1].protein_id,
        legacy_subject[0].protein_id: extraction.proteins_by_record[1][0].protein_id,
    }
    assert all(runtime_handle.startswith("h_") for runtime_handle in resolved.values())


@pytest.mark.linear
def test_legacy_protein_references_resolve_mixed_record_orientations() -> None:
    record = _record(
        "mixed-orientation",
        sequence="ATG" * 12,
        features=[
            _cds(0, 9, qualifiers={"translation": ["MKT"]}),
            _cds(12, 24, strand=-1, qualifiers={"translation": ["MGGG"]}),
        ],
    )
    extraction = extract_protein_identity_manifest(
        [record],
        record_instance_keys=("current-row",),
    )
    forward = protein_colinearity_module._with_stable_web_protein_ids(
        extraction.proteins_by_record[0],
        "r_old",
    )
    reversed_extraction = extract_protein_identity_manifest(
        reverse_records([record], True),
        record_instance_keys=("historical-reversed-row",),
    )
    reversed_legacy = protein_colinearity_module._with_stable_web_protein_ids(
        reversed_extraction.proteins_by_record[0],
        "r_old",
    )
    reversed_reference = next(
        protein.protein_id
        for protein in reversed_legacy
        if protein.sequence == "MGGG"
    )
    current_by_sequence = {
        protein.sequence: protein.protein_id
        for protein in extraction.proteins_by_record[0]
    }

    assert build_legacy_protein_reference_map(
        extraction,
        [forward[0].protein_id, reversed_reference],
    ) == {
        forward[0].protein_id: current_by_sequence["MKT"],
        reversed_reference: current_by_sequence["MGGG"],
    }


@pytest.mark.linear
def test_legacy_protein_artifact_reference_resolution_rejects_ambiguous_records() -> None:
    feature = _cds(0, 9, qualifiers={"translation": ["MKT"]})
    extraction = extract_protein_identity_manifest(
        [
            _record("same", features=[feature]),
            _record("same", features=[feature]),
        ],
        record_instance_keys=("first", "second"),
    )
    legacy = protein_colinearity_module._with_stable_web_protein_ids(
        extraction.proteins_by_record[0],
        "r_old",
    )

    with pytest.raises(ValidationError, match="resolves to 2 current record instances"):
        build_legacy_protein_reference_map(
            extraction,
            [legacy[0].protein_id],
        )


@pytest.mark.linear
def test_web_helper_resolves_legacy_protein_artifact_references_with_python_owner() -> None:
    records = [
        _record("query", features=[_cds(0, 9)]),
        _record("subject", features=[_cds(3, 12)]),
    ]
    extraction = extract_protein_identity_manifest(
        records,
        record_instance_keys=("left", "right"),
    )
    manifest = extraction.identity_manifest
    assert manifest is not None
    legacy_query = protein_colinearity_module._with_stable_web_protein_ids(
        extraction.proteins_by_record[0],
        "r_old_left",
    )[0]
    helpers = _load_web_helper_namespace()
    serialize_protein = helpers["_serialize_cds_protein"]
    raw_records = [
        {
            "proteinMap": {
                protein.protein_id: serialize_protein(protein, records[record_index])
                for protein in proteins
            },
            "fasta": proteins_to_fasta(proteins),
        }
        for record_index, proteins in enumerate(extraction.proteins_by_record)
    ]

    result = json.loads(
        helpers["resolve_legacy_protein_reference_map_json"](
            json.dumps(raw_records),
            json.dumps(manifest.to_dict()),
            json.dumps([legacy_query.protein_id]),
        )
    )

    assert result == {
        "status": "resolved",
        "proteinIdMap": {
            legacy_query.protein_id: extraction.proteins_by_record[0][0].protein_id,
        },
    }


@pytest.mark.linear
def test_protein_raw_tsv_hydration_uses_aliases_and_duplicate_ordinals() -> None:
    query = _record(
        "query",
        features=[
            _cds(
                0,
                9,
                qualifiers={"translation": ["MKT"], "protein_id": ["duplicate"]},
            ),
            _cds(
                9,
                18,
                qualifiers={"translation": ["MGG"], "protein_id": ["duplicate"]},
            ),
        ],
    )
    subject = _record(
        "subject",
        features=[
            _cds(
                0,
                9,
                qualifiers={"translation": ["MQQ"], "protein_id": ["target id"]},
            )
        ],
    )
    extraction = extract_protein_identity_manifest(
        [query, subject],
        record_instance_keys=("query-instance", "subject-instance"),
    )
    manifest = extraction.identity_manifest
    assert manifest is not None
    pair = build_protein_losat_pair_identity(
        manifest,
        query_record_instance_key="query-instance",
        subject_record_instance_key="subject-instance",
    )
    query_handles = [
        protein.protein_id for protein in extraction.proteins_by_record[0]
    ]
    subject_handle = extraction.proteins_by_record[1][0].protein_id
    suffix = "99.1\t3\t0\t0\t1\t3\t1\t3\t1e-5\t20"
    text = (
        "# generated\r\n"
        f"{query_handles[0]}\t{subject_handle}\t{suffix}\r\n"
        "\r\n"
        f"{query_handles[1]}\t{subject_handle}\t{suffix}"
    )
    entry = {
        "schema": PROTEIN_LOSAT_CACHE_SCHEMA,
        "kind": "raw-losat",
        "identityKind": "protein",
        "idEncoding": "runtime-handle-v1",
        "key": build_protein_losat_cache_key(pair, args=[]),
        "text": text,
        "program": "blastp",
        "outfmt": "6",
        "args": [],
        "queryProteinSetHash": pair.query_protein_set_hash,
        "subjectProteinSetHash": pair.subject_protein_set_hash,
        "queryRuntimeBindingHash": pair.query_runtime_binding_hash,
        "subjectRuntimeBindingHash": pair.subject_runtime_binding_hash,
        "queryRecordInstanceKey": pair.query_record_instance_key,
        "subjectRecordInstanceKey": pair.subject_record_instance_key,
    }

    query_exports = build_protein_export_id_map(manifest, "query-instance")
    assert set(query_exports.values()) == {"duplicate~1", "duplicate~2"}
    assert build_protein_export_id_map(
        manifest,
        "subject-instance",
    )[subject_handle] == "target%20id"
    hydrated = hydrate_protein_losat_tsv(entry, manifest)
    assert hydrated.startswith("# generated\r\n")
    assert "\r\n\r\n" in hydrated
    assert not any(handle in hydrated for handle in query_handles)
    assert subject_handle not in hydrated
    raw_rows = [
        line.split("\t")
        for line in text.splitlines()
        if line and not line.startswith("#")
    ]
    hydrated_rows = [
        line.split("\t")
        for line in hydrated.splitlines()
        if line and not line.startswith("#")
    ]
    assert [row[2:] for row in hydrated_rows] == [row[2:] for row in raw_rows]
    assert all(len(row) == 12 for row in hydrated_rows)
    with pytest.raises(ValidationError, match="does not resolve"):
        hydrate_protein_losat_tsv(
            {**entry, "text": text.replace(query_handles[0], "h_aaaaaaaaaaaaaaaaaaaaaaaaaa")},
            manifest,
        )


@pytest.mark.linear
def test_extract_cds_proteins_uses_translation_and_stable_synthetic_ids() -> None:
    records = [
        _record(
            "record_a",
            features=[
                _cds(
                    0,
                    9,
                    qualifiers={
                        "translation": ["MKT*"],
                        "protein_id": ["duplicate"],
                        "locus_tag": ["gene_a"],
                    },
                )
            ],
        ),
        _record(
            "record_b",
            features=[
                _cds(
                    9,
                    18,
                    qualifiers={
                        "translation": ["MGG*"],
                        "protein_id": ["duplicate"],
                    },
                )
            ],
        ),
    ]

    result = extract_cds_proteins(records)

    assert [protein.protein_id for protein in result.protein_map.values()] == [
        "gbd_r0001_cds000001",
        "gbd_r0002_cds000001",
    ]
    first = result.proteins_by_record[0][0]
    assert first.sequence == "MKT"
    assert first.label == "gene_a"
    assert first.start == 0
    assert first.end == 9


@pytest.mark.linear
def test_extract_cds_proteins_prefers_unique_source_protein_id() -> None:
    record = _record(
        "record_a",
        features=[
            _cds(
                0,
                9,
                qualifiers={
                    "translation": ["MKT*"],
                    "protein_id": ["WP_123456789.1"],
                },
            )
        ],
    )

    result = extract_cds_proteins([record])

    protein = result.proteins_by_record[0][0]
    assert protein.protein_id == "WP_123456789.1"
    assert protein.source_protein_id == "WP_123456789.1"
    assert ">WP_123456789.1" in proteins_to_fasta([protein])


@pytest.mark.linear
def test_extract_cds_proteins_carries_annotation_fields() -> None:
    record = _record(
        "record_a",
        features=[
            _cds(
                0,
                9,
                qualifiers={
                    "translation": ["MKT*"],
                    "gene": ["rpoB"],
                    "product": ["DNA-directed RNA polymerase beta subunit"],
                    "note": ["core polymerase subunit"],
                },
            )
        ],
    )

    result = extract_cds_proteins([record])

    protein = result.proteins_by_record[0][0]
    assert protein.gene == "rpoB"
    assert protein.product == "DNA-directed RNA polymerase beta subunit"
    assert protein.note == "core polymerase subunit"


@pytest.mark.linear
def test_extract_cds_proteins_accepts_record_index_offset() -> None:
    record = _record("record_b", features=[_cds(9, 18)])

    result = extract_cds_proteins([record], record_index_offset=1)

    protein = result.proteins_by_record[0][0]
    assert protein.protein_id == "gbd_r0002_cds000001"
    assert protein.record_index == 1


@pytest.mark.linear
def test_extract_cds_proteins_translates_when_translation_missing() -> None:
    record = _record(
        "record_a",
        sequence="ATGAAATAG",
        features=[_cds(0, 9, qualifiers={"locus_tag": ["fallback"]})],
    )

    result = extract_cds_proteins([record])

    assert result.proteins_by_record[0][0].sequence == "MK"


@pytest.mark.linear
def test_extract_cds_proteins_handles_compound_location_span() -> None:
    feature = SeqFeature(
        CompoundLocation(
            [
                FeatureLocation(0, 6, strand=1),
                FeatureLocation(12, 18, strand=1),
            ]
        ),
        type="CDS",
        qualifiers={"translation": ["MKM"]},
    )
    record = _record("record_a", features=[feature])

    result = extract_cds_proteins([record])

    protein = result.proteins_by_record[0][0]
    assert protein.start == 0
    assert protein.end == 18
    assert protein.strand == 1


@pytest.mark.linear
def test_filter_protein_hits_by_thresholds_removes_low_confidence_bridge() -> None:
    raw_hits = pd.DataFrame.from_records(
        [
            _hit_row("a", "b"),
            _hit_row("b", "c", bitscore=20, evalue=1e-1),
        ],
        columns=COMPARISON_COLUMNS,
    )
    filtered_hits = filter_protein_hits_by_thresholds(
        raw_hits,
        bitscore=50,
        evalue=1e-5,
        identity=0,
        alignment_length=0,
    )
    assert filtered_hits[["query", "subject"]].to_records(index=False).tolist() == [
        ("a", "b")
    ]

    records = [
        _record("record_a", features=[_cds(0, 9)]),
        _record("record_b", features=[_cds(9, 18)]),
        _record("record_c", features=[_cds(18, 27)]),
    ]
    calls = 0

    def fake_runner(_query_fasta: str, _subject_fasta: str) -> pd.DataFrame:
        nonlocal calls
        calls += 1
        if calls == 1:
            rows = [_hit_row("gbd_r0001_cds000001", "gbd_r0002_cds000001")]
        else:
            rows = [
                _hit_row(
                    "gbd_r0002_cds000001",
                    "gbd_r0003_cds000001",
                    bitscore=20,
                    evalue=1e-1,
                )
            ]
        return pd.DataFrame.from_records(rows, columns=COMPARISON_COLUMNS)

    result = build_pairwise_protein_blastp_comparisons(
        records,
        runner=fake_runner,
        bitscore=50,
        evalue=1e-5,
        identity=0,
        alignment_length=0,
    )

    assert result.orthogroups is None
    assert result.comparisons[0].iloc[0]["orthogroup_id"] == ""
    assert result.comparisons[1].empty


@pytest.mark.linear
def test_orthogroup_membership_modes_accept_legacy_aliases_for_anchor_core() -> None:
    assert protein_colinearity_module.ORTHOGROUP_MEMBERSHIP_MODES == ("anchor_core_v1",)
    assert protein_colinearity_module.normalize_orthogroup_membership_mode("rbh") == "anchor_core_v1"
    assert protein_colinearity_module.normalize_orthogroup_membership_mode("family_merge") == "anchor_core_v1"
    assert protein_colinearity_module.normalize_orthogroup_membership_mode("local-split") == "anchor_core_v1"
    assert protein_colinearity_module.normalize_orthogroup_membership_mode("density_split") == "anchor_core_v1"
    assert protein_colinearity_module.normalize_orthogroup_membership_mode("outparalog_split") == "anchor_core_v1"


@pytest.mark.linear
def test_select_best_hits_per_query_avoids_secondary_paralog_merge() -> None:
    hits = pd.DataFrame.from_records(
        [
            _hit_row("a", "b1", bitscore=300),
            _hit_row("a", "b2", bitscore=120),
        ],
        columns=COMPARISON_COLUMNS,
    )

    selected = select_best_hits_per_query(hits)

    assert selected[["query", "subject"]].to_records(index=False).tolist() == [
        ("a", "b1")
    ]


@pytest.mark.linear
def test_select_reciprocal_best_hits_rejects_many_to_one_nonreciprocal_hit() -> None:
    hits = pd.DataFrame.from_records(
        [
            _hit_row("a1", "b", bitscore=300),
            _hit_row("a2", "b", bitscore=250),
        ],
        columns=COMPARISON_COLUMNS,
    )

    selected = select_reciprocal_best_hits(hits)

    assert selected[["query", "subject"]].to_records(index=False).tolist() == [
        ("a1", "b")
    ]


@pytest.mark.linear
def test_select_reciprocal_best_hit_edges_requires_directional_reciprocity() -> None:
    forward = pd.DataFrame.from_records(
        [
            _hit_row("a1", "b", bitscore=300),
            _hit_row("a2", "b", bitscore=250),
        ],
        columns=COMPARISON_COLUMNS,
    )
    reverse = pd.DataFrame.from_records(
        [_hit_row("b", "a3", bitscore=400)],
        columns=COMPARISON_COLUMNS,
    )

    selected = select_reciprocal_best_hit_edges(forward, reverse)

    assert selected.empty


@pytest.mark.linear
def test_select_top_hits_per_query_is_deterministic_on_subject_tie() -> None:
    hits = pd.DataFrame.from_records(
        [
            _hit_row("a", "b2", bitscore=300),
            _hit_row("a", "b1", bitscore=300),
        ],
        columns=COMPARISON_COLUMNS,
    )

    selected = select_top_hits_per_query(hits, max_hits=1)

    assert selected.iloc[0]["subject"] == "b1"


@pytest.mark.linear
def test_convert_protein_hits_to_genomic_links_uses_cds_spans_and_strand() -> None:
    records = [
        _record("query_record", features=[_cds(0, 9, strand=1)]),
        _record("subject_record", features=[_cds(20, 50, strand=-1)]),
    ]
    extraction = extract_cds_proteins(records)
    hits = pd.DataFrame.from_records(
        [
            _hit_row(
                "gbd_r0001_cds000001",
                "gbd_r0002_cds000001",
                alignment_length=20,
            )
        ],
        columns=COMPARISON_COLUMNS,
    )

    converted = convert_protein_hits_to_genomic_links(hits, extraction.protein_map)

    row = converted.iloc[0]
    assert row["query"] == "query_record"
    assert row["subject"] == "subject_record"
    assert row["qstart"] == 1
    assert row["qend"] == 9
    assert row["sstart"] == 50
    assert row["send"] == 21
    assert row["alignment_length"] == 20


@pytest.mark.linear
def test_separately_extracted_offset_protein_maps_convert_without_collision() -> None:
    query_record = _record("query_record", features=[_cds(0, 9, strand=1)])
    subject_record = _record("subject_record", features=[_cds(50, 80, strand=1)])
    query_extraction = extract_cds_proteins([query_record], record_index_offset=0)
    subject_extraction = extract_cds_proteins([subject_record], record_index_offset=1)
    protein_map = {
        **query_extraction.protein_map,
        **subject_extraction.protein_map,
    }
    hits = pd.DataFrame.from_records(
        [
            _hit_row(
                "gbd_r0001_cds000001",
                "gbd_r0002_cds000001",
                alignment_length=20,
            )
        ],
        columns=COMPARISON_COLUMNS,
    )

    converted = convert_protein_hits_to_genomic_links(hits, protein_map)

    row = converted.iloc[0]
    assert row["query"] == "query_record"
    assert row["subject"] == "subject_record"
    assert row["qstart"] == 1
    assert row["qend"] == 9
    assert row["sstart"] == 51
    assert row["send"] == 80


@pytest.mark.linear
def test_pair_conversion_disambiguates_same_source_protein_id_between_records() -> None:
    query_record = _record(
        "query_record",
        features=[_cds(0, 9, strand=1, qualifiers={"translation": ["MK*"], "protein_id": ["WP_same.1"]})],
    )
    subject_record = _record(
        "subject_record",
        features=[_cds(50, 80, strand=1, qualifiers={"translation": ["MK*"], "protein_id": ["WP_same.1"]})],
    )
    query_extraction = extract_cds_proteins([query_record], record_index_offset=0)
    subject_extraction = extract_cds_proteins([subject_record], record_index_offset=1)
    hits = pd.DataFrame.from_records(
        [_hit_row("WP_same.1", "WP_same.1", alignment_length=20)],
        columns=COMPARISON_COLUMNS,
    )

    converted = convert_pair_protein_hits_to_genomic_links(
        hits,
        query_extraction.protein_map,
        subject_extraction.protein_map,
    )

    row = converted.iloc[0]
    assert row["query"] == "query_record"
    assert row["subject"] == "subject_record"
    assert row["qstart"] == 1
    assert row["qend"] == 9
    assert row["sstart"] == 51
    assert row["send"] == 80


@pytest.mark.linear
def test_parse_losatp_outfmt6_returns_standard_columns() -> None:
    parsed = parse_losatp_outfmt6(
        "# comment\nq1\ts1\t91.5\t42\t1\t0\t1\t42\t2\t43\t1e-20\t150\n"
    )

    assert parsed.columns.tolist() == list(COMPARISON_COLUMNS)
    assert parsed.iloc[0]["identity"] == pytest.approx(91.5)
    assert parsed.iloc[0]["bitscore"] == pytest.approx(150)


def _protein_map_for_lengths(lengths: dict[str, int]) -> dict[str, protein_colinearity_module.CdsProtein]:
    return {
        protein_id: protein_colinearity_module.CdsProtein(
            protein_id=protein_id,
            record_index=index,
            feature_index=0,
            record_id=f"record_{index}",
            start=0,
            end=int(length) * 3,
            strand=1,
            label=protein_id,
            protein_length=int(length),
            sequence="M" * int(length),
        )
        for index, (protein_id, length) in enumerate(lengths.items())
    }


def _assert_hsp_baseline(hits, protein_map):
    from tests.prototypes.hsp_aggregation import aggregate_hsps_baseline

    before = None if hits is None else hits.copy(deep=True)
    try:
        expected = aggregate_hsps_baseline(hits, protein_map)
    except (ValueError, TypeError, OverflowError, KeyError) as exc:
        with pytest.raises(type(exc)) as caught:
            protein_colinearity_module._aggregate_hsps_by_protein_pair(hits, protein_map)
        assert str(caught.value) == str(exc)
    else:
        actual = protein_colinearity_module._aggregate_hsps_by_protein_pair(hits, protein_map)
        pd.testing.assert_frame_equal(actual, expected, check_exact=True)
    if hits is not None:
        pd.testing.assert_frame_equal(hits, before, check_exact=True)


@pytest.mark.parametrize("seed", range(24))
def test_hsp_accumulation_matches_grouped_oracle(seed):
    import random

    rng = random.Random(seed)
    pm = _protein_map_for_lengths({"q": 100, "s": 80, "t": 200, "zero": 0})
    rows = []
    for i in range(120):
        q, s = rng.choice(list(pm)), rng.choice(list(pm))
        rows.append(_hit_row(q, s, bitscore=rng.choice([0, -1, 50, 100]),
            identity=rng.choice([80, 90]), evalue=rng.choice([0, 1e-10]),
            alignment_length=rng.choice([-1, 0, 20, 40]),
            qstart=rng.randint(-30, 230), qend=rng.randint(-30, 230),
            sstart=rng.randint(-30, 230), send=rng.randint(-30, 230)) | {"source_row": i})
    hits = pd.DataFrame(rows)
    hits.index = [7] * len(hits)  # tie order is positional, not the index label
    _assert_hsp_baseline(hits, pm)


@pytest.mark.parametrize("column", ["bitscore", "evalue", "identity", "alignment_length",
                                    "qstart", "qend", "sstart", "send"])
@pytest.mark.parametrize("value", [float("nan"), float("inf"), -float("inf"), "bad", None, pd.NA])
def test_hsp_private_boundary_and_public_numeric_rejection(column, value):
    from gbdraw.exceptions import ParseError

    pm = _protein_map_for_lengths({"q": 100, "s": 100})
    hits = pd.DataFrame([_hit_row("q", "s"), _hit_row("q", "s") | {column: value}])
    _assert_hsp_baseline(hits, pm)
    for limit in (None, 1):
        with pytest.raises(ParseError, match="non-numeric outfmt 6"):
            select_rbh_orthogroup_edges_from_directional_hits(
                {(0, 1): hits}, pm, orthogroup_member_max_hits=limit)


@pytest.mark.parametrize("dtype", ["object", "string", "category"])
def test_hsp_missing_ids_never_become_valid_strings(dtype):
    pm = _protein_map_for_lengths({"q": 100, "s": 100, "nan": 100, "None": 100, "<NA>": 100})
    rows = [_hit_row(q, s) for q, s in [(None, "s"), ("q", None), (float("nan"), "s"),
            (pd.NA, "s"), ("unknown", "s"), ("nan", "s"), ("q", "s"), ("nan", "s")]]
    hits = pd.DataFrame(rows).astype({"query": dtype, "subject": dtype})
    _assert_hsp_baseline(hits, pm)
    result = protein_colinearity_module._aggregate_hsps_by_protein_pair(hits, pm)
    assert list(zip(result["query"], result["subject"], result["hsp_count"])) == [
        ("nan", "s", 2), ("q", "s", 1)]


@pytest.mark.parametrize("hits", [None, pd.DataFrame(columns=COMPARISON_COLUMNS),
    pd.DataFrame(columns=[*COMPARISON_COLUMNS, "extra"]),
    pd.DataFrame([_hit_row("unknown", "s")]),
    pd.DataFrame([_hit_row("q", "zero")]), pd.DataFrame([_hit_row("q", "negative")])])
def test_hsp_empty_output_columns_and_types(hits):
    _assert_hsp_baseline(hits, _protein_map_for_lengths({"q": 1, "s": 100, "zero": 0, "negative": -1}))


@pytest.mark.parametrize("column,winner", [("bitscore", 101), ("evalue", 0), ("identity", 100),
    ("alignment_length", 101), ("qstart", 0), ("qend", 99), ("sstart", 0), ("send", 99), (None, None)])
def test_hsp_representative_rank_and_complete_tie_order(column, winner):
    pm = _protein_map_for_lengths({"q": 200, "s": 200})
    first = _hit_row("q", "s", bitscore=100, evalue=1e-20, identity=90,
                     alignment_length=100, qstart=1, qend=100, sstart=1, send=100)
    second = first | ({column: winner} if column else {})
    hits = pd.DataFrame([first | {"tag": "first"}, _hit_row("s", "q"), second | {"tag": "second"}])
    _assert_hsp_baseline(hits, pm)
    result = protein_colinearity_module._aggregate_hsps_by_protein_pair(hits, pm)
    assert result["query"].tolist() == ["q", "s"]
    assert result.iloc[0]["tag"] == ("second" if column else "first")


def test_hsp_adjacent_reverse_clamped_union_and_all_selected_hsps(monkeypatch):
    from tests.prototypes.hsp_aggregation import aggregate_hsps_baseline
    from tools.benchmark_protein_comparison import canonical

    pc = protein_colinearity_module
    pm = _protein_map_for_lengths({"q": 100, "s": 100, "t": 100})
    hits = pd.DataFrame([
        _hit_row("q", "s", bitscore=200, alignment_length=30, qstart=-10, qend=20, sstart=20, send=-10),
        _hit_row("q", "t", bitscore=100),
        _hit_row("q", "s", bitscore=190, alignment_length=30, qstart=21, qend=50, sstart=50, send=21),
        _hit_row("q", "s", bitscore=180, alignment_length=40, qstart=150, qend=61, sstart=61, send=150)])
    for limit in (1, None):
        selected = pc._select_member_candidate_hits_per_query(hits, max_hits=limit)
        result = pc._normalize_directional_hit_table(selected, pm, min_coverage=0.0)
        pair = result.loc[result["subject"] == "s"].iloc[0]
        assert pair["hsp_count"] == 3
        assert pair["query_covered_length"] == pair["subject_covered_length"] == 90
        assert pair["total_hsp_alignment_length"] == 100
        assert len(result) == (1 if limit == 1 else 2)
        # String numerics are coerced at the existing normalization boundary.
        string_hits = hits.astype({c: str for c in pc._NUMERIC_COMPARISON_COLUMNS})
        args = ({(0, 1): string_hits}, pm)
        actual = select_rbh_orthogroup_edges_from_directional_hits(*args, orthogroup_member_max_hits=limit)
        with monkeypatch.context() as patch:
            patch.setattr(pc, "_aggregate_hsps_by_protein_pair", aggregate_hsps_baseline)
            expected = select_rbh_orthogroup_edges_from_directional_hits(*args, orthogroup_member_max_hits=limit)
        assert canonical(actual) == canonical(expected)


@pytest.mark.parametrize("first_error,later_error", [(float("nan"), float("inf")), (float("inf"), float("nan"))])
def test_hsp_errors_follow_pair_appearance_not_interleaved_error_row(first_error, later_error):
    pm = _protein_map_for_lengths({"q": 100, "s": 100, "t": 100})
    hits = pd.DataFrame([_hit_row("q", "s"),
        _hit_row("q", "t", alignment_length=later_error),
        _hit_row("q", "s", alignment_length=first_error)])
    _assert_hsp_baseline(hits, pm)


@pytest.mark.parametrize("identifiers", [[1, 2.0], [1, 1.0], [True, 1], ["", "nan"], [None, pd.NaT]])
def test_hsp_group_key_coercion_matches_pandas(identifiers):
    pm = _protein_map_for_lengths({"1": 100, "s": 100, "": 100, "nan": 100})
    hits = pd.DataFrame([_hit_row(q, "s") for q in identifiers])
    hits["query"] = pd.Series(identifiers, dtype=object)
    _assert_hsp_baseline(hits, pm)


def test_hsp_aggregation_iterates_table_once_without_pair_dataframes(monkeypatch):
    from tools.benchmark_protein_comparison import synthetic, SEED

    pc = protein_colinearity_module
    pm, tables = synthetic(pc, "hsp-1000", SEED)
    hits = tables["multi_hsp"]
    original = pd.DataFrame.itertuples
    seen = []
    def iterate(df, *args, **kwargs):
        seen.append(len(df))
        return original(df, *args, **kwargs)
    from pandas.core.groupby.ops import FrameSplitter
    def grouped(*args, **kwargs):
        pytest.fail("HSP aggregation must not construct per-pair DataFrames")
    monkeypatch.setattr(pd.DataFrame, "itertuples", iterate)
    monkeypatch.setattr(FrameSplitter, "_chop", grouped)
    result = pc._aggregate_hsps_by_protein_pair(hits, pm)
    assert seen == [3000]
    assert result["hsp_count"].tolist() == [3] * 1000


@pytest.mark.linear
def test_hsp_union_coverage_uses_merged_intervals() -> None:
    protein_map = _protein_map_for_lengths({"BDT62853.1": 4741, "BDT62565.1": 4468})
    hits = pd.DataFrame.from_records(
        [
            _hit_row("BDT62853.1", "BDT62565.1", identity=58.191, alignment_length=1172, qstart=1200, qend=2318, sstart=1039, send=2075, bitscore=800),
            _hit_row("BDT62853.1", "BDT62565.1", identity=55.136, alignment_length=1032, qstart=33, qend=1034, sstart=48, send=964, bitscore=700),
            _hit_row("BDT62853.1", "BDT62565.1", identity=64.078, alignment_length=824, qstart=2774, qend=3577, sstart=2524, send=3242, bitscore=600),
            _hit_row("BDT62853.1", "BDT62565.1", identity=48.235, alignment_length=340, qstart=2417, qend=2709, sstart=2118, send=2436, bitscore=500),
        ],
        columns=COMPARISON_COLUMNS,
    )

    normalized = protein_colinearity_module._normalize_directional_hit_table(
        hits,
        protein_map,
        min_coverage=0.0,
    )

    assert normalized.shape[0] == 1
    row = normalized.iloc[0]
    assert row["hsp_count"] == 4
    assert row["coverage_source"] == "hsp_union"
    assert row["query_coverage"] > 0.30
    assert row["subject_coverage"] > 0.30
    assert row["min_coverage"] > 0.30
    assert not bool(row["domain_only"])


@pytest.mark.linear
def test_overlapping_hsps_do_not_double_count_coverage() -> None:
    protein_map = _protein_map_for_lengths({"q": 1000, "s": 1000})
    hits = pd.DataFrame.from_records(
        [
            _hit_row("q", "s", alignment_length=300, qstart=1, qend=300, sstart=1, send=300, bitscore=300),
            _hit_row("q", "s", alignment_length=301, qstart=200, qend=500, sstart=200, send=500, bitscore=250),
        ],
        columns=COMPARISON_COLUMNS,
    )

    normalized = protein_colinearity_module._normalize_directional_hit_table(
        hits,
        protein_map,
        min_coverage=0.0,
    )

    row = normalized.iloc[0]
    assert row["hsp_count"] == 2
    assert row["query_covered_length"] == 500
    assert row["subject_covered_length"] == 500
    assert row["total_hsp_alignment_length"] == 601
    assert row["query_coverage"] == pytest.approx(0.5)
    assert row["subject_coverage"] == pytest.approx(0.5)


@pytest.mark.linear
def test_invalid_hsp_coordinates_do_not_inflate_aggregate_coverage() -> None:
    protein_map = _protein_map_for_lengths({"q": 1000, "s": 1000})
    hits = pd.DataFrame.from_records(
        [
            _hit_row("q", "s", alignment_length=500, qstart="bad", qend="bad", sstart="bad", send="bad", bitscore=500),
            _hit_row("q", "s", alignment_length=200, qstart=1, qend=200, sstart=1, send=200, bitscore=300),
        ],
        columns=COMPARISON_COLUMNS,
    )

    aggregated = protein_colinearity_module._aggregate_hsps_by_protein_pair(
        hits,
        protein_map,
    )

    row = aggregated.iloc[0]
    assert row["hsp_count"] == 2
    assert row["bitscore"] == 500
    assert row["query_covered_length"] == 200
    assert row["subject_covered_length"] == 200
    assert row["min_coverage"] == pytest.approx(0.2)


@pytest.mark.linear
def test_run_losatp_blastp_passes_num_threads(monkeypatch: pytest.MonkeyPatch) -> None:
    captured: dict[str, object] = {}

    def fake_run(command, **kwargs):
        captured["command"] = command
        captured["kwargs"] = kwargs
        return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

    monkeypatch.setattr(protein_colinearity_module.subprocess, "run", fake_run)

    protein_colinearity_module.run_losatp_blastp(
        ">query\nM\n",
        ">subject\nM\n",
        losatp_bin="custom-losat",
        max_hits=3,
        threads=4,
    )

    command = captured["command"]
    assert command[0:2] == ["custom-losat", "blastp"]
    assert "-max_hsps" in command
    assert command[command.index("-max_hsps") + 1] == "1"
    assert "-max_target_seqs" in command
    assert command[command.index("-max_target_seqs") + 1] == "3"
    assert "-num_threads" in command
    assert command[command.index("-num_threads") + 1] == "4"


def _write_conda_losat(prefix: Path, *, executable: bool = True) -> Path:
    (prefix / "conda-meta").mkdir(parents=True)
    candidate = prefix / "bin" / "losat"
    candidate.parent.mkdir()
    candidate.write_text("#!/bin/sh\n", encoding="utf-8")
    candidate.chmod(0o755 if executable else 0o644)
    return candidate.absolute()


@pytest.mark.linear
def test_conda_losat_precedes_cache_bundled_and_path_without_side_effects(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    prefix = tmp_path / "conda-b"
    candidate = _write_conda_losat(prefix)
    original_bytes = candidate.read_bytes()
    original_mode = candidate.stat().st_mode

    monkeypatch.setattr(protein_colinearity_module.sys, "prefix", str(prefix))
    monkeypatch.setattr(
        losat_setup_module,
        "managed_losat",
        lambda: (_ for _ in ()).throw(AssertionError("managed cache was read")),
    )
    monkeypatch.setattr(
        protein_colinearity_module,
        "_bundled_losatp_resource",
        lambda: (_ for _ in ()).throw(AssertionError("bundled discovery ran")),
    )
    monkeypatch.setattr(
        protein_colinearity_module,
        "_path_executable",
        lambda _name: (_ for _ in ()).throw(AssertionError("PATH discovery ran")),
    )

    with ExitStack() as stack:
        runtime = protein_colinearity_module._resolve_protein_blastp_runtime(
            "losat", None, stack
        )

    assert runtime.source == "conda"
    assert runtime.executable == str(candidate)
    assert candidate.read_bytes() == original_bytes
    assert candidate.stat().st_mode == original_mode
    assert not (candidate.parent / "LOSAT").exists()


@pytest.mark.linear
def test_conda_losat_identity_uses_sys_prefix_not_environment_or_path(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    prefix_a = tmp_path / "conda-a"
    prefix_b = tmp_path / "conda-b"
    candidate_a = _write_conda_losat(prefix_a)
    candidate_b = _write_conda_losat(prefix_b)
    monkeypatch.setattr(protein_colinearity_module.sys, "prefix", str(prefix_b))
    monkeypatch.setenv("CONDA_PREFIX", str(prefix_a))
    monkeypatch.setenv("PATH", f"{candidate_a.parent}{os.pathsep}{os.environ.get('PATH', '')}")

    with ExitStack() as stack:
        runtime = protein_colinearity_module._resolve_protein_blastp_runtime(
            "losat", None, stack
        )

    assert runtime.source == "conda"
    assert runtime.executable == str(candidate_b)


@pytest.mark.linear
def test_conda_prefix_without_losat_uses_normal_path_fallback(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    prefix_a = tmp_path / "conda-a"
    prefix_b = tmp_path / "conda-b"
    candidate_a = _write_conda_losat(prefix_a)
    (prefix_b / "conda-meta").mkdir(parents=True)
    monkeypatch.setattr(protein_colinearity_module.sys, "prefix", str(prefix_b))
    monkeypatch.setenv("CONDA_PREFIX", str(prefix_a))
    monkeypatch.setattr(losat_setup_module, "managed_losat", lambda: None)
    monkeypatch.setattr(protein_colinearity_module, "_bundled_losatp_resource", lambda: None)
    monkeypatch.setattr(
        protein_colinearity_module,
        "_path_executable",
        lambda name: str(candidate_a) if name == "losat" else None,
    )

    with ExitStack() as stack:
        runtime = protein_colinearity_module._resolve_protein_blastp_runtime(
            "losat", None, stack
        )

    assert runtime.source == "path"
    assert runtime.executable == str(candidate_a)


@pytest.mark.linear
@pytest.mark.parametrize("conda_marker_kind", ["missing", "file"])
def test_non_conda_or_venv_prefix_ignores_base_conda_binary(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    conda_marker_kind: str,
) -> None:
    base_prefix = tmp_path / "base-conda"
    _write_conda_losat(base_prefix)
    venv_prefix = tmp_path / "venv"
    local_candidate = venv_prefix / "bin" / "losat"
    local_candidate.parent.mkdir(parents=True)
    local_candidate.write_text("not selected", encoding="utf-8")
    if conda_marker_kind == "file":
        (venv_prefix / "conda-meta").write_text("not a directory", encoding="utf-8")
    path_candidate = tmp_path / "path" / "losat"
    monkeypatch.setattr(protein_colinearity_module.sys, "prefix", str(venv_prefix))
    monkeypatch.setattr(protein_colinearity_module.sys, "base_prefix", str(base_prefix))
    monkeypatch.setenv("CONDA_PREFIX", str(base_prefix))
    monkeypatch.setattr(losat_setup_module, "managed_losat", lambda: None)
    monkeypatch.setattr(protein_colinearity_module, "_bundled_losatp_resource", lambda: None)
    monkeypatch.setattr(
        protein_colinearity_module,
        "_path_executable",
        lambda name: str(path_candidate) if name == "losat" else None,
    )

    with ExitStack() as stack:
        runtime = protein_colinearity_module._resolve_protein_blastp_runtime(
            "losat", None, stack
        )

    assert runtime.source == "path"
    assert runtime.executable == str(path_candidate)


@pytest.mark.linear
@pytest.mark.parametrize(
    ("candidate_kind", "message"),
    [
        ("broken-link", "broken symbolic link"),
        ("directory", "not a regular file"),
        ("non-executable", "not executable"),
    ],
)
def test_invalid_conda_losat_stops_at_its_path_without_fallback(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    candidate_kind: str,
    message: str,
) -> None:
    prefix = tmp_path / "conda"
    (prefix / "conda-meta").mkdir(parents=True)
    candidate = prefix / "bin" / "losat"
    candidate.parent.mkdir()
    if candidate_kind == "broken-link":
        candidate.symlink_to(prefix / "missing-losat")
    elif candidate_kind == "directory":
        candidate.mkdir()
    else:
        candidate.write_text("#!/bin/sh\n", encoding="utf-8")
        candidate.chmod(0o644)

    monkeypatch.setattr(protein_colinearity_module.sys, "prefix", str(prefix))
    monkeypatch.setattr(
        losat_setup_module,
        "managed_losat",
        lambda: (_ for _ in ()).throw(AssertionError("managed fallback ran")),
    )
    monkeypatch.setattr(
        protein_colinearity_module,
        "_path_executable",
        lambda _name: (_ for _ in ()).throw(AssertionError("PATH fallback ran")),
    )

    with ExitStack() as stack, pytest.raises(ValidationError, match=message) as exc_info:
        protein_colinearity_module._resolve_protein_blastp_runtime(
            "losat", None, stack
        )

    assert str(candidate.absolute()) in str(exc_info.value)


@pytest.mark.linear
def test_conda_losat_execution_failure_does_not_fall_back(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    prefix = tmp_path / "conda"
    candidate = _write_conda_losat(prefix)
    monkeypatch.setattr(protein_colinearity_module.sys, "prefix", str(prefix))
    monkeypatch.setattr(
        protein_colinearity_module,
        "_path_executable",
        lambda _name: (_ for _ in ()).throw(AssertionError("PATH fallback ran")),
    )
    monkeypatch.setattr(
        protein_colinearity_module.subprocess,
        "run",
        lambda *_args, **_kwargs: (_ for _ in ()).throw(OSError("exec format error")),
    )

    with pytest.raises(ValidationError, match="could not be started") as exc_info:
        protein_colinearity_module.run_losatp_blastp(
            ">query\nM\n",
            ">subject\nM\n",
        )

    assert str(candidate) in str(exc_info.value)


@pytest.mark.linear
def test_conda_losat_nonzero_exit_reports_selected_path(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    prefix = tmp_path / "conda"
    candidate = _write_conda_losat(prefix)
    monkeypatch.setattr(protein_colinearity_module.sys, "prefix", str(prefix))
    monkeypatch.setattr(
        protein_colinearity_module.subprocess,
        "run",
        lambda command, **_kwargs: subprocess.CompletedProcess(
            command, 2, stdout="", stderr="bad subject"
        ),
    )

    with pytest.raises(ValidationError, match="bad subject") as exc_info:
        protein_colinearity_module.run_losatp_blastp(
            ">query\nM\n",
            ">subject\nM\n",
        )

    assert str(candidate) in str(exc_info.value)


@pytest.mark.linear
def test_run_losatp_blastp_uses_bundled_binary_by_default(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    bundled_losat = tmp_path / "losat"
    bundled_losat.write_text("#!/bin/sh\n", encoding="utf-8")
    captured: dict[str, object] = {}

    def fake_run(command, **kwargs):
        captured["command"] = command
        captured["kwargs"] = kwargs
        return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

    monkeypatch.setattr(
        protein_colinearity_module,
        "_bundled_losatp_resource",
        lambda: bundled_losat,
    )
    monkeypatch.setattr(protein_colinearity_module, "_conda_losatp_runtime", lambda: None)
    monkeypatch.setattr(losat_setup_module, "managed_losat", lambda: None)
    monkeypatch.setattr(protein_colinearity_module.subprocess, "run", fake_run)

    protein_colinearity_module.run_losatp_blastp(
        ">query\nM\n",
        ">subject\nM\n",
    )

    command = captured["command"]
    assert command[0:2] == [str(bundled_losat), "blastp"]


@pytest.mark.linear
def test_run_losatp_blastp_uses_path_losat_before_blastp(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, object] = {}

    def fake_run(command, **kwargs):
        captured["command"] = command
        captured["kwargs"] = kwargs
        return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

    def fake_path_executable(name: str) -> str | None:
        return {"losat": "/usr/local/bin/losat", "blastp": "/usr/local/bin/blastp"}.get(name)

    monkeypatch.setattr(protein_colinearity_module, "_bundled_losatp_resource", lambda: None)
    monkeypatch.setattr(protein_colinearity_module, "_conda_losatp_runtime", lambda: None)
    monkeypatch.setattr(losat_setup_module, "managed_losat", lambda: None)
    monkeypatch.setattr(protein_colinearity_module, "_path_executable", fake_path_executable)
    monkeypatch.setattr(protein_colinearity_module.subprocess, "run", fake_run)

    protein_colinearity_module.run_losatp_blastp(
        ">query\nM\n",
        ">subject\nM\n",
    )

    command = captured["command"]
    assert command[0:2] == ["/usr/local/bin/losat", "blastp"]


@pytest.mark.linear
def test_run_losatp_blastp_falls_back_to_path_ncbi_blastp(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, object] = {}

    def fake_run(command, **kwargs):
        captured["command"] = command
        captured["kwargs"] = kwargs
        return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

    def fake_path_executable(name: str) -> str | None:
        return "/usr/local/bin/blastp" if name == "blastp" else None

    monkeypatch.setattr(protein_colinearity_module, "_bundled_losatp_resource", lambda: None)
    monkeypatch.setattr(protein_colinearity_module, "_conda_losatp_runtime", lambda: None)
    monkeypatch.setattr(losat_setup_module, "managed_losat", lambda: None)
    monkeypatch.setattr(protein_colinearity_module, "_path_executable", fake_path_executable)
    monkeypatch.setattr(protein_colinearity_module.subprocess, "run", fake_run)

    protein_colinearity_module.run_losatp_blastp(
        ">query\nM\n",
        ">subject\nM\n",
        max_hits=4,
        max_hsps_per_subject=2,
        threads=3,
    )

    command = captured["command"]
    assert command[0] == "/usr/local/bin/blastp"
    assert command[1:6:2] == ["-query", "-subject", "-outfmt"]
    assert "-max_hsps" in command
    assert command[command.index("-max_hsps") + 1] == "2"
    assert "-max_target_seqs" in command
    assert command[command.index("-max_target_seqs") + 1] == "4"
    assert "-num_threads" in command
    assert command[command.index("-num_threads") + 1] == "3"


@pytest.mark.linear
def test_run_losatp_blastp_explicit_ncbi_blastp_bypasses_losat_discovery(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, object] = {}

    def fail_bundled():
        raise AssertionError("bundled LOSAT discovery should not run")

    def fail_path(_name: str):
        raise AssertionError("PATH discovery should not run")

    def fake_run(command, **kwargs):
        captured["command"] = command
        captured["kwargs"] = kwargs
        return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

    monkeypatch.setattr(protein_colinearity_module, "_bundled_losatp_resource", fail_bundled)
    monkeypatch.setattr(
        protein_colinearity_module,
        "_conda_losatp_runtime",
        lambda: (_ for _ in ()).throw(AssertionError("conda discovery should not run")),
    )
    monkeypatch.setattr(protein_colinearity_module, "_path_executable", fail_path)
    monkeypatch.setattr(protein_colinearity_module.subprocess, "run", fake_run)

    protein_colinearity_module.run_losatp_blastp(
        ">query\nM\n",
        ">subject\nM\n",
        ncbi_blastp_bin="/opt/ncbi/bin/blastp",
        max_hits=2,
    )

    command = captured["command"]
    assert command[0] == "/opt/ncbi/bin/blastp"
    assert "blastp" not in command[1:]
    assert "-max_target_seqs" in command
    assert "--max-target-seqs" not in command


@pytest.mark.linear
def test_run_losatp_blastp_rejects_ambiguous_explicit_runtimes() -> None:
    with pytest.raises(ValidationError, match="either --losatp_bin or --ncbi_blastp_bin"):
        protein_colinearity_module.run_losatp_blastp(
            ">query\nM\n",
            ">subject\nM\n",
            losatp_bin="custom-losat",
            ncbi_blastp_bin="/opt/ncbi/bin/blastp",
        )


@pytest.mark.linear
def test_run_losatp_blastp_explicit_losat_bypasses_discovery(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, object] = {}

    def fail_bundled():
        raise AssertionError("bundled LOSAT discovery should not run")

    def fail_path(_name: str):
        raise AssertionError("PATH discovery should not run")

    def fake_run(command, **kwargs):
        captured["command"] = command
        captured["kwargs"] = kwargs
        return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

    monkeypatch.setattr(protein_colinearity_module, "_bundled_losatp_resource", fail_bundled)
    monkeypatch.setattr(
        protein_colinearity_module,
        "_conda_losatp_runtime",
        lambda: (_ for _ in ()).throw(AssertionError("conda discovery should not run")),
    )
    monkeypatch.setattr(protein_colinearity_module, "_path_executable", fail_path)
    monkeypatch.setattr(protein_colinearity_module.subprocess, "run", fake_run)

    protein_colinearity_module.run_losatp_blastp(
        ">query\nM\n",
        ">subject\nM\n",
        losatp_bin="custom-losat",
    )

    assert captured["command"][0:2] == ["custom-losat", "blastp"]


@pytest.mark.linear
def test_run_losatp_blastp_missing_runtime_error_is_actionable(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(protein_colinearity_module, "_bundled_losatp_resource", lambda: None)
    monkeypatch.setattr(protein_colinearity_module, "_conda_losatp_runtime", lambda: None)
    monkeypatch.setattr(losat_setup_module, "managed_losat", lambda: None)
    monkeypatch.setattr(protein_colinearity_module, "_path_executable", lambda _name: None)
    monkeypatch.setattr(
        protein_colinearity_module,
        "_bundled_losatp_platform_dir",
        lambda: "macos-arm64",
    )

    with pytest.raises(ValidationError) as exc_info:
        protein_colinearity_module.run_losatp_blastp(
            ">query\nM\n",
            ">subject\nM\n",
        )

    message = str(exc_info.value)
    assert "Protein BLASTP comparison needs LOSAT or NCBI BLAST+" in message
    assert "macos-arm64" in message
    assert "`losat` was not found on PATH" in message
    assert "`blastp` was not found on PATH" in message
    assert "--ncbi_blastp_bin" in message


@pytest.mark.linear
def test_bundled_losatp_platform_dir_is_platform_specific(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(protein_colinearity_module.sys, "platform", "linux")
    monkeypatch.setattr(protein_colinearity_module.platform, "machine", lambda: "x86_64")
    assert protein_colinearity_module._bundled_losatp_platform_dir() == "linux-x86_64"

    monkeypatch.setattr(protein_colinearity_module.sys, "platform", "darwin")
    monkeypatch.setattr(protein_colinearity_module.platform, "machine", lambda: "arm64")
    assert protein_colinearity_module._bundled_losatp_platform_dir() == "macos-arm64"


@pytest.mark.linear
def test_run_losatp_blastp_omits_hsp_cap_when_requested(monkeypatch: pytest.MonkeyPatch) -> None:
    captured: dict[str, object] = {}

    def fake_run(command, **kwargs):
        captured["command"] = command
        captured["kwargs"] = kwargs
        return subprocess.CompletedProcess(command, 0, stdout="", stderr="")

    monkeypatch.setattr(protein_colinearity_module.subprocess, "run", fake_run)

    protein_colinearity_module.run_losatp_blastp(
        ">query\nM\n",
        ">subject\nM\n",
        losatp_bin="custom-losat",
        max_hsps_per_subject=None,
    )

    assert "-max_hsps" not in captured["command"]


@pytest.mark.linear
def test_build_losat_blastp_command_uses_losat_flags() -> None:
    command = protein_colinearity_module._build_losat_blastp_command(
        executable="losat",
        query_path=Path("query.faa"),
        subject_path=Path("subject.faa"),
        max_hits=5,
        max_hsps_per_subject=2,
        threads=4,
    )

    assert command[:2] == ["losat", "blastp"]
    assert "-max_hsps" in command
    assert "-max_target_seqs" in command
    assert "-num_threads" in command


@pytest.mark.linear
def test_build_ncbi_blastp_command_uses_ncbi_flags() -> None:
    command = protein_colinearity_module._build_ncbi_blastp_command(
        executable="blastp",
        query_path=Path("query.faa"),
        subject_path=Path("subject.faa"),
        max_hits=5,
        max_hsps_per_subject=2,
        threads=4,
    )

    assert command[0] == "blastp"
    assert command[1] == "-query"
    assert "-max_hsps" in command
    assert "-max_target_seqs" in command
    assert "-num_threads" in command


@pytest.mark.linear
def test_build_blastp_commands_omit_hsp_cap_when_requested() -> None:
    losat_command = protein_colinearity_module._build_losat_blastp_command(
        executable="losat",
        query_path=Path("query.faa"),
        subject_path=Path("subject.faa"),
        max_hits=None,
        max_hsps_per_subject=None,
        threads=None,
    )
    ncbi_command = protein_colinearity_module._build_ncbi_blastp_command(
        executable="blastp",
        query_path=Path("query.faa"),
        subject_path=Path("subject.faa"),
        max_hits=None,
        max_hsps_per_subject=None,
        threads=None,
    )

    assert "-max_hsps" not in losat_command
    assert "-max_hsps" not in ncbi_command


@pytest.mark.linear
def test_run_losatp_blastp_calls_raw_output_callback_for_ncbi(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    seen: dict[str, str] = {}

    def fake_run(command, **kwargs):
        return subprocess.CompletedProcess(command, 0, stdout="q\ts\t100\t1\t0\t0\t1\t1\t1\t1\t0\t10\n", stderr="")

    monkeypatch.setattr(protein_colinearity_module.subprocess, "run", fake_run)

    result = protein_colinearity_module.run_losatp_blastp(
        ">query\nM\n",
        ">subject\nM\n",
        ncbi_blastp_bin="/opt/ncbi/bin/blastp",
        raw_output_callback=lambda text: seen.__setitem__("text", text),
    )

    assert seen["text"].startswith("q\ts\t100")
    assert tuple(result.columns) == COMPARISON_COLUMNS
    assert len(result) == 1


@pytest.mark.linear
def test_run_losatp_blastp_ncbi_failure_includes_stderr(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fake_run(command, **kwargs):
        return subprocess.CompletedProcess(command, 2, stdout="", stderr="bad subject")

    monkeypatch.setattr(protein_colinearity_module.subprocess, "run", fake_run)

    with pytest.raises(ValidationError, match="NCBI BLAST\\+ blastp failed.*bad subject"):
        protein_colinearity_module.run_losatp_blastp(
            ">query\nM\n",
            ">subject\nM\n",
            ncbi_blastp_bin="/opt/ncbi/bin/blastp",
        )


@pytest.mark.linear
def test_external_ncbi_blastp_smoke_when_enabled() -> None:
    if os.environ.get("GBDRAW_RUN_EXTERNAL_BLASTP_SMOKE") != "1":
        pytest.skip("Set GBDRAW_RUN_EXTERNAL_BLASTP_SMOKE=1 to run external blastp smoke test.")
    blastp_bin = os.environ.get("GBDRAW_EXTERNAL_BLASTP_BIN") or shutil.which("blastp")
    if not blastp_bin:
        pytest.skip("NCBI BLAST+ blastp executable not found.")

    result = protein_colinearity_module.run_losatp_blastp(
        ">query\nMKTAYIAKQRQISFVKSHFSRQDILD\n",
        ">subject\nMKTAYIAKQRQISFVKSHFSRQDILD\n",
        ncbi_blastp_bin=blastp_bin,
    )

    assert tuple(result.columns) == COMPARISON_COLUMNS


@pytest.mark.linear
def test_build_pairwise_protein_blastp_comparisons_accepts_test_runner_without_orthogroups() -> None:
    records = [
        _record("record_a", features=[_cds(0, 9)]),
        _record("record_b", features=[_cds(9, 18)]),
    ]

    def fake_runner(query_fasta: str, subject_fasta: str) -> pd.DataFrame:
        assert ">gbd_r0001_cds000001" in query_fasta
        assert ">gbd_r0002_cds000001" in subject_fasta
        return pd.DataFrame.from_records(
            [_hit_row("gbd_r0001_cds000001", "gbd_r0002_cds000001")],
            columns=COMPARISON_COLUMNS,
        )

    result = build_pairwise_protein_blastp_comparisons(
        records,
        runner=fake_runner,
    )

    comparisons = result.comparisons
    assert result.orthogroups is None
    assert len(comparisons) == 1
    assert comparisons[0].iloc[0]["query"] == "record_a"
    assert comparisons[0].iloc[0]["subject"] == "record_b"
    assert comparisons[0].iloc[0]["orthogroup_id"] == ""
    assert comparisons[0].iloc[0]["query_protein_id"] == "gbd_r0001_cds000001"
    assert comparisons[0].iloc[0]["subject_protein_id"] == "gbd_r0002_cds000001"


@pytest.mark.linear
def test_pairwise_blastp_uses_web_losat_cache_without_external_run(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records = [
        _record("record_a", features=[_cds(0, 9)]),
        _record("record_b", features=[_cds(9, 18)]),
    ]
    extraction = extract_web_stable_cds_proteins(
        records,
        record_instance_keys=("r_left", "r_right"),
    )
    query_fasta = proteins_to_fasta(extraction.proteins_by_record[0])
    subject_fasta = proteins_to_fasta(extraction.proteins_by_record[1])
    query_id = extraction.proteins_by_record[0][0].protein_id
    subject_id = extraction.proteins_by_record[1][0].protein_id
    pair_identity = build_protein_losat_pair_identity(
        extraction.identity_manifest,
        query_record_instance_key="r_left",
        subject_record_instance_key="r_right",
    )
    cache_key = build_protein_losat_cache_key(
        pair_identity,
        args=["--max-hsps-per-subject", "1"],
    )
    raw_text = (
        f"{query_id}\t{subject_id}\t90\t100\t0\t0\t1\t100\t1\t100\t1e-20\t200\n"
    )
    cache = LosatpCacheManager(
        [
            {
                "schema": PROTEIN_LOSAT_CACHE_SCHEMA,
                "kind": "raw-losat",
                "identityKind": "protein",
                "idEncoding": "runtime-handle-v1",
                "key": cache_key,
                "text": raw_text,
                "program": "blastp",
                "outfmt": "6",
                "args": ["--max-hsps-per-subject", "1"],
                "queryProteinSetHash": pair_identity.query_protein_set_hash,
                "subjectProteinSetHash": pair_identity.subject_protein_set_hash,
                "queryRuntimeBindingHash": pair_identity.query_runtime_binding_hash,
                "subjectRuntimeBindingHash": pair_identity.subject_runtime_binding_hash,
                "queryRecordInstanceKey": pair_identity.query_record_instance_key,
                "subjectRecordInstanceKey": pair_identity.subject_record_instance_key,
            }
        ],
        identity_manifest=extraction.identity_manifest,
    )

    def fail_run(*_args, **_kwargs):
        raise AssertionError("LOSATP should not run on a cache hit")

    monkeypatch.setattr(protein_colinearity_module, "run_losatp_blastp", fail_run)

    result = build_pairwise_protein_blastp_comparisons(
        records,
        losatp_cache=cache,
        protein_extraction=extraction,
        cache_filenames=("record_a.record_b.losatp.tsv",),
    )

    assert result.comparisons[0].iloc[0]["query_protein_id"] == query_id
    entries = cache.session_entries()
    assert entries[0]["key"] == cache_key
    assert entries[0]["display"] is True
    assert entries[0]["filename"] == "record_a.record_b.losatp.tsv"
    assert validate_protein_raw_entry_references(
        entries[0],
        extraction.identity_manifest,
    )
    assert validate_protein_raw_entry_references(
        {**entries[0], "text": ""},
        extraction.identity_manifest,
    )
    assert not validate_protein_raw_entry_references(
        {**entries[0], "text": raw_text.replace(subject_id, "outside-binding")},
        extraction.identity_manifest,
    )
    with pytest.raises(ValidationError, match="IDs and sequences"):
        cache._pair_identity_from_fasta(
            query_fasta.replace("\nMK\n", "\nMA\n"),
            subject_fasta,
        )


@pytest.mark.linear
def test_pairwise_blastp_reuses_only_verified_legacy_web_losat_cache(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records = [
        _record("record_a", features=[_cds(0, 9)]),
        _record("record_b", features=[_cds(9, 18)]),
    ]
    extraction = extract_web_stable_cds_proteins(
        records,
        record_instance_keys=("r_left", "r_right"),
    )
    query_fasta = proteins_to_fasta(extraction.proteins_by_record[0])
    subject_fasta = proteins_to_fasta(extraction.proteins_by_record[1])
    query_id = extraction.proteins_by_record[0][0].protein_id
    subject_id = extraction.proteins_by_record[1][0].protein_id
    legacy_query = protein_colinearity_module._with_stable_web_protein_ids(
        extraction.proteins_by_record[0],
        "r_old_left",
    )
    legacy_subject = protein_colinearity_module._with_stable_web_protein_ids(
        extraction.proteins_by_record[1],
        "r_old_right",
    )
    legacy_query_fasta = proteins_to_fasta(legacy_query)
    legacy_subject_fasta = proteins_to_fasta(legacy_subject)
    legacy_key, query_hash, subject_hash = build_web_losat_cache_key(
        query_fasta=legacy_query_fasta,
        subject_fasta=legacy_subject_fasta,
        args=["--max-hsps-per-subject", "1"],
    )
    raw_text = (
        f"{legacy_query[0].protein_id}\t{legacy_subject[0].protein_id}"
        "\t90\t100\t0\t0\t1\t100\t1\t100\t1e-20\t200\n"
    )
    legacy_entry = {
        "schema": 2,
        "kind": "raw-losat",
        "key": legacy_key,
        "text": raw_text,
        "program": "blastp",
        "outfmt": "6",
        "args": ["--max-hsps-per-subject", "1"],
        "queryCanonicalHash": query_hash,
        "subjectCanonicalHash": subject_hash,
    }
    promotion = promote_legacy_protein_raw_cache_entries(
        [legacy_entry],
        query_proteins=extraction.proteins_by_record[0],
        subject_proteins=extraction.proteins_by_record[1],
        query_fasta=query_fasta,
        subject_fasta=subject_fasta,
        identity_manifest=extraction.identity_manifest,
        expected_args=["--max-hsps-per-subject", "1"],
    )
    assert promotion.promotion is not None
    assert legacy_query[0].protein_id not in promotion.promotion.rewritten_tsv
    cache = LosatpCacheManager([legacy_entry], threads_per_job=32)
    cache.set_protein_extraction(extraction)

    def fail_run(*_args, **_kwargs):
        raise AssertionError("LOSATP should not run on a legacy cache hit")

    monkeypatch.setattr(protein_colinearity_module, "run_losatp_blastp", fail_run)

    result = build_pairwise_protein_blastp_comparisons(
        records,
        losatp_cache=cache,
        protein_extraction=extraction,
        cache_filenames=("record_a.record_b.losatp.tsv",),
    )

    assert result.comparisons[0].iloc[0]["query_protein_id"] == query_id
    assert result.comparisons[0].iloc[0]["subject_protein_id"] == subject_id
    entries = cache.session_entries()
    assert entries[0]["key"] != legacy_key
    assert entries[0]["schema"] == PROTEIN_LOSAT_CACHE_SCHEMA
    assert entries[0]["display"] is True
    assert entries[0]["args"] == ["--max-hsps-per-subject", "1"]
    assert entries[0]["outfmt"] == "6"
    assert cache.has_legacy_candidates is False


@pytest.mark.linear
def test_legacy_protein_cache_promotion_rejects_hash_args_and_ambiguous_empty_output() -> None:
    records = [
        _record("record_a", features=[_cds(0, 9)]),
        _record("record_b", features=[_cds(9, 18)]),
    ]
    extraction = extract_protein_identity_manifest(
        records,
        record_instance_keys=("left", "right"),
    )
    query_proteins, subject_proteins = extraction.proteins_by_record
    query_fasta = proteins_to_fasta(query_proteins)
    subject_fasta = proteins_to_fasta(subject_proteins)
    legacy_query = protein_colinearity_module._with_stable_web_protein_ids(
        query_proteins,
        "legacy_left",
    )
    legacy_subject = protein_colinearity_module._with_stable_web_protein_ids(
        subject_proteins,
        "legacy_right",
    )
    _, query_hash, subject_hash = build_web_losat_cache_key(
        query_fasta=proteins_to_fasta(legacy_query),
        subject_fasta=proteins_to_fasta(legacy_subject),
        args=[],
    )
    valid_text = (
        f"{legacy_query[0].protein_id}\t{legacy_subject[0].protein_id}"
        "\t90\t2\t0\t0\t1\t2\t1\t2\t1e-5\t20\n"
    )
    common = {
        "schema": 2,
        "kind": "raw-losat",
        "program": "blastp",
        "outfmt": "6",
        "queryCanonicalHash": query_hash,
        "subjectCanonicalHash": subject_hash,
    }
    candidates = [
        {**common, "key": "wrong-args", "args": ["--different"], "text": valid_text},
        {**common, "key": "wrong-hash", "args": [], "text": valid_text, "queryCanonicalHash": "0" * 64},
        {
            **common,
            "key": "empty",
            "args": [],
            "text": "",
            "queryCanonicalHash": "unproven-query",
            "subjectCanonicalHash": "unproven-subject",
        },
    ]

    scan = promote_legacy_protein_raw_cache_entries(
        candidates,
        query_proteins=query_proteins,
        subject_proteins=subject_proteins,
        query_fasta=query_fasta,
        subject_fasta=subject_fasta,
        identity_manifest=extraction.identity_manifest,
        expected_args=[],
    )

    assert scan.promotion is None
    assert len(scan.rejections) == 3
    assert "args" in scan.rejections[0].reason
    assert "hash" in scan.rejections[1].reason
    assert "Empty" in scan.rejections[2].reason
    cache = LosatpCacheManager(candidates, identity_manifest=extraction.identity_manifest)
    assert cache.session_entries() == ()
    assert len(cache.legacy_candidate_envelope()["entries"]) == 3

    empty_with_evidence = promote_legacy_protein_raw_cache_entries(
        [
            {
                **common,
                "key": "token-evidence",
                "args": ["--different"],
                "text": valid_text,
            },
            {**common, "key": "verified-empty", "args": [], "text": ""},
        ],
        query_proteins=query_proteins,
        subject_proteins=subject_proteins,
        query_fasta=query_fasta,
        subject_fasta=subject_fasta,
        identity_manifest=extraction.identity_manifest,
        expected_args=[],
    )
    assert empty_with_evidence.promotion is not None
    assert empty_with_evidence.promotion.candidate_index == 1
    assert empty_with_evidence.promotion.rewritten_tsv == ""


@pytest.mark.linear
def test_linear_cli_save_session_writes_web_losat_cache_entries(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    records = [
        _record("record_a", features=[_cds(0, 9)]),
        _record("record_b", features=[_cds(9, 18)]),
    ]
    input_a = tmp_path / "a.gb"
    input_b = tmp_path / "b.gb"
    input_a.write_text("LOCUS       A\n", encoding="utf-8")
    input_b.write_text("LOCUS       B\n", encoding="utf-8")
    output_prefix = tmp_path / "out"
    captured: dict[str, object] = {}

    records_by_path = {
        str(input_a): records[0],
        str(input_b): records[1],
    }

    def fake_load_gbks(paths, **_kwargs):
        return [records_by_path[str(path)] for path in paths]

    monkeypatch.setattr(request_render_module, "load_gbks", fake_load_gbks)
    monkeypatch.setattr(request_render_module, "read_color_table", lambda _path: None)
    monkeypatch.setattr(request_render_module, "read_feature_visibility_file", lambda _path: None)

    real_render = linear_cli_module.render_request

    def capture_render(canonical_request, **kwargs):
        result = real_render(canonical_request, **kwargs)
        captured["canonical_request"] = result.request
        return result

    def fake_losatp(query_fasta: str, subject_fasta: str, **kwargs) -> pd.DataFrame:
        query_id = query_fasta.splitlines()[0][1:].split()[0]
        subject_id = subject_fasta.splitlines()[0][1:].split()[0]
        raw_text = (
            f"{query_id}\t{subject_id}\t90\t100\t0\t0\t1\t100\t1\t100\t1e-20\t200\n"
        )
        callback = kwargs.get("raw_output_callback")
        if callback is not None:
            callback(raw_text)
        return parse_losatp_outfmt6(raw_text)

    monkeypatch.setattr(linear_cli_module, "render_request", capture_render)
    monkeypatch.setattr(protein_colinearity_module, "run_losatp_blastp", fake_losatp)

    linear_cli_module.linear_main(
        [
            "--gbk",
            str(input_a),
            str(input_b),
            "--protein_blastp_mode",
            "pairwise",
            "-o",
            str(output_prefix),
            "-f",
            "svg",
            "--save_session",
        ]
    )

    payload = json.loads(output_prefix.with_suffix(".gbdraw-session.json").read_text(encoding="utf-8"))
    entries = payload["losatCache"]["entries"]
    assert len(entries) == 1
    assert entries[0]["schema"] == PROTEIN_LOSAT_CACHE_SCHEMA
    assert entries[0]["kind"] == "raw-losat"
    assert entries[0]["identityKind"] == "protein"
    assert entries[0]["display"] is True
    assert entries[0]["filename"] == "record_a.record_b.losatp.tsv"
    assert entries[0]["text"].strip()
    assert entries[0]["queryProteinSetHash"].startswith("sha256:")
    assert entries[0]["subjectProteinSetHash"].startswith("sha256:")
    assert payload["proteinIdentityManifest"]["schema"] == 2
    canonical_request = captured["canonical_request"]
    assert isinstance(canonical_request, LinearDiagramRequest)
    assert canonical_request.options.protein_blastp_mode == "pairwise"


@pytest.mark.linear
def test_linear_cli_writes_hydrated_raw_protein_evidence_and_honors_overwrite(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    records = [
        _record(
            "record_a",
            features=[
                _cds(
                    0,
                    9,
                    qualifiers={
                        "translation": ["MKT"],
                        "protein_id": ["query protein"],
                    },
                )
            ],
        ),
        _record(
            "record_b",
            features=[
                _cds(
                    9,
                    18,
                    qualifiers={
                        "translation": ["MKT"],
                        "protein_id": ["subject/protein"],
                    },
                )
            ],
        ),
    ]
    input_a = tmp_path / "a.gb"
    input_b = tmp_path / "b.gb"
    input_a.write_text("LOCUS       A\n", encoding="utf-8")
    input_b.write_text("LOCUS       B\n", encoding="utf-8")
    output_prefix = tmp_path / "out"
    evidence_path = tmp_path / "raw-evidence.tsv"
    records_by_path = {str(input_a): records[0], str(input_b): records[1]}

    monkeypatch.setattr(
        request_render_module,
        "load_gbks",
        lambda paths, **_kwargs: [records_by_path[str(path)] for path in paths],
    )
    monkeypatch.setattr(request_render_module, "read_color_table", lambda _path: None)
    monkeypatch.setattr(
        request_render_module,
        "read_feature_visibility_file",
        lambda _path: None,
    )
    run_count = 0

    def fake_losatp(query_fasta: str, subject_fasta: str, **kwargs) -> pd.DataFrame:
        nonlocal run_count
        run_count += 1
        query_id = query_fasta.splitlines()[0][1:].split()[0]
        subject_id = subject_fasta.splitlines()[0][1:].split()[0]
        raw_text = (
            f"{query_id}\t{subject_id}\t90\t3\t0\t0\t1\t3\t1\t3\t1e-20\t200\n"
        )
        callback = kwargs.get("raw_output_callback")
        if callback is not None:
            callback(raw_text)
        return parse_losatp_outfmt6(raw_text)

    monkeypatch.setattr(protein_colinearity_module, "run_losatp_blastp", fake_losatp)
    args = [
        "--gbk",
        str(input_a),
        str(input_b),
        "--protein_blastp_mode",
        "pairwise",
        "--protein_blastp_output",
        str(evidence_path),
        "-o",
        str(output_prefix),
        "-f",
        "svg",
    ]

    linear_cli_module.linear_main(args)

    evidence = evidence_path.read_text(encoding="utf-8")
    rows = [line.split("\t") for line in evidence.splitlines() if not line.startswith("#")]
    assert "# entry 1: record_a.record_b.losatp.tsv" in evidence
    assert rows == [
        [
            "query%20protein",
            "subject%2Fprotein",
            "90",
            "3",
            "0",
            "0",
            "1",
            "3",
            "1",
            "3",
            "1e-20",
            "200",
        ]
    ]
    assert "h_" not in evidence
    assert run_count == 1

    with pytest.raises(ValidationError, match="already exist"):
        linear_cli_module.linear_main(args)
    assert run_count == 1

    linear_cli_module.linear_main([*args, "--overwrite"])
    assert run_count == 2


@pytest.mark.linear
def test_linear_cli_validates_raw_protein_output_option(
    capsys: pytest.CaptureFixture[str],
) -> None:
    with pytest.raises(SystemExit):
        linear_cli_module._get_args(["--help"])
    assert "--protein_blastp_output TSV" in capsys.readouterr().out

    with pytest.raises(SystemExit, match="2"):
        linear_cli_module._get_args(
            ["--gbk", "a.gb", "b.gb", "--protein_blastp_output", "raw.tsv"]
        )
    with pytest.raises(SystemExit, match="2"):
        linear_cli_module._get_args(
            [
                "--gbk",
                "a.gb",
                "b.gb",
                "--protein_blastp_mode",
                "pairwise",
                "--protein_blastp_output",
                "raw.txt",
            ]
        )
    parsed = linear_cli_module._get_args(
        [
            "--gbk",
            "a.gb",
            "b.gb",
            "--protein_blastp_mode",
            "pairwise",
            "--protein_blastp_output",
            "raw.tsv",
        ]
    )
    assert parsed.protein_blastp_output == "raw.tsv"


@pytest.mark.linear
def test_pairwise_blastp_search_keeps_one_hsp_cap(monkeypatch: pytest.MonkeyPatch) -> None:
    records = [
        _record("record_a", features=[_cds(0, 9)]),
        _record("record_b", features=[_cds(9, 18)]),
    ]
    observed_caps: list[int | None] = []

    def fake_search(
        query_fasta,
        subject_fasta,
        *,
        ncbi_blastp_bin,
        max_hsps_per_subject,
        **_kwargs,
    ):
        assert ncbi_blastp_bin is None
        observed_caps.append(max_hsps_per_subject)
        return pd.DataFrame.from_records(
            [_hit_row("gbd_r0001_cds000001", "gbd_r0002_cds000001")],
            columns=COMPARISON_COLUMNS,
        )

    monkeypatch.setattr(protein_colinearity_module, "run_losatp_blastp", fake_search)

    build_pairwise_protein_blastp_comparisons(records)

    assert observed_caps == [1]


@pytest.mark.linear
def test_orthogroup_blastp_search_omits_one_hsp_cap(monkeypatch: pytest.MonkeyPatch) -> None:
    records = [
        _record("record_a", features=[_cds(0, 9)]),
        _record("record_b", features=[_cds(9, 18)]),
    ]
    observed_caps: list[int | None] = []

    def fake_search(
        query_fasta,
        subject_fasta,
        *,
        ncbi_blastp_bin,
        max_hsps_per_subject,
        **_kwargs,
    ):
        assert ncbi_blastp_bin is None
        observed_caps.append(max_hsps_per_subject)
        query_id = query_fasta.splitlines()[0][1:].split()[0]
        subject_id = subject_fasta.splitlines()[0][1:].split()[0]
        return pd.DataFrame.from_records(
            [_hit_row(query_id, subject_id)],
            columns=COMPARISON_COLUMNS,
        )

    monkeypatch.setattr(protein_colinearity_module, "run_losatp_blastp", fake_search)

    build_rbh_orthogroup_protein_blastp_comparisons(records)

    assert observed_caps == [None, None, None, None]


@pytest.mark.linear
@pytest.mark.parametrize("candidate_limit", (None, 7))
def test_candidate_and_pairwise_display_limits_reach_independent_consumers(
    candidate_limit: int | None,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    records = [
        _record("record_a", features=[_cds(0, 9)]),
        _record(
            "record_b",
            features=[_cds(index * 12, index * 12 + 9) for index in range(3)],
        ),
    ]
    observed_raw_limits: list[int | None] = []

    def fake_search(query_fasta, subject_fasta, *, max_hits, **_kwargs):
        observed_raw_limits.append(max_hits)
        query_id = query_fasta.splitlines()[0][1:].split()[0]
        subject_ids = [
            line[1:].split()[0]
            for line in subject_fasta.splitlines()
            if line.startswith(">")
        ]
        return pd.DataFrame.from_records(
            [
                _hit_row(query_id, subject_id, bitscore=300 - index)
                for index, subject_id in enumerate(subject_ids)
            ],
            columns=COMPARISON_COLUMNS,
        )

    monkeypatch.setattr(protein_colinearity_module, "run_losatp_blastp", fake_search)

    result = build_pairwise_protein_blastp_comparisons(
        records,
        max_hits=2,
        candidate_limit=candidate_limit,
    )

    assert observed_raw_limits == [candidate_limit]
    assert len(result.comparisons[0]) == 2


@pytest.mark.linear
@pytest.mark.parametrize("member_limit", [None, 2, 5])
def test_member_hit_limit_bounds_derived_candidates_and_preserves_their_hsps(
    monkeypatch: pytest.MonkeyPatch,
    member_limit: int | None,
) -> None:
    records = [
        _record("record_a", features=[_cds(0, 9)]),
        _record(
            "record_b",
            features=[_cds(index * 12, index * 12 + 9) for index in range(7)],
        ),
    ]
    extraction = extract_cds_proteins(records)
    query_id = extraction.proteins_by_record[0][0].protein_id
    subject_ids = [protein.protein_id for protein in extraction.proteins_by_record[1]]
    directional_hits = {
        (0, 1): pd.DataFrame.from_records(
            [
                *[_hit_row(query_id, subject_id, bitscore=300 - index * 10)
                  for index, subject_id in enumerate(subject_ids)],
                _hit_row(query_id, subject_ids[1], bitscore=280),
            ],
            columns=COMPARISON_COLUMNS,
        ),
    }
    captured: dict[tuple[int, int], pd.DataFrame] = {}
    real_consumer = (
        protein_colinearity_module._select_anchor_core_orthogroup_edges_from_directional_hits
    )

    def spy_consumer(member_hits, *args, **kwargs):
        captured.update({pair: hits.copy() for pair, hits in member_hits.items()})
        return real_consumer(member_hits, *args, **kwargs)

    monkeypatch.setattr(
        protein_colinearity_module,
        "_select_anchor_core_orthogroup_edges_from_directional_hits",
        spy_consumer,
    )

    select_rbh_orthogroup_edges_from_directional_hits(
        directional_hits,
        extraction.protein_map,
        record_count=2,
        orthogroup_member_max_hits=member_limit,
    )

    observed = captured[(0, 1)]
    assert set(observed["subject"]) == set(subject_ids[:member_limit])
    assert len(observed.loc[observed["subject"] == subject_ids[1]]) == 2


def _long_cds(protein_id: str, length: int = 1000) -> SeqFeature:
    return _cds(
        0,
        int(length) * 3,
        qualifiers={
            "translation": ["M" * int(length) + "*"],
            "protein_id": [protein_id],
            "locus_tag": [protein_id],
        },
    )


def _long_record(record_id: str, protein_id: str, length: int = 1000) -> SeqRecord:
    return _record(
        record_id,
        sequence="ATG" * int(length),
        features=[_long_cds(protein_id, length)],
    )


@pytest.mark.linear
def test_anchor_core_assigns_member_with_multi_hsp_union_coverage() -> None:
    records = [
        _long_record("record_a", "a0"),
        _long_record("record_b", "b0"),
        _long_record("record_c", "c0"),
    ]
    extraction = extract_cds_proteins(records)
    directional_hits = {
        (0, 1): pd.DataFrame.from_records(
            [_hit_row("a0", "b0", alignment_length=900, qstart=1, qend=900, sstart=1, send=900, bitscore=1000)],
            columns=COMPARISON_COLUMNS,
        ),
        (1, 0): pd.DataFrame.from_records(
            [_hit_row("b0", "a0", alignment_length=900, qstart=1, qend=900, sstart=1, send=900, bitscore=1000)],
            columns=COMPARISON_COLUMNS,
        ),
        (2, 1): pd.DataFrame.from_records(
            [
                _hit_row("c0", "b0", alignment_length=180, qstart=1, qend=180, sstart=1, send=180, bitscore=320),
                _hit_row("c0", "b0", alignment_length=191, qstart=220, qend=410, sstart=220, send=410, bitscore=300),
                _hit_row("c0", "b0", alignment_length=201, qstart=500, qend=700, sstart=500, send=700, bitscore=280),
            ],
            columns=COMPARISON_COLUMNS,
        ),
    }

    result = select_rbh_orthogroup_edges_from_directional_hits(
        directional_hits,
        extraction.protein_map,
        record_count=len(records),
        orthogroup_membership_mode="anchor_core_v1",
    ).orthogroups

    assert result.member_by_protein_id["c0"].orthogroup_id == result.member_by_protein_id["a0"].orthogroup_id
    assert result.member_by_protein_id["c0"].role == "coortholog"


@pytest.mark.linear
def test_anchor_core_prefers_membership_support_over_domain_only_top_hit() -> None:
    records = [
        _long_record("record_a", "a0"),
        _long_record("record_b", "b0"),
        _long_record("record_c", "c0"),
    ]
    extraction = extract_cds_proteins(records)
    directional_hits = {
        (0, 1): pd.DataFrame.from_records(
            [_hit_row("a0", "b0", alignment_length=900, qstart=1, qend=900, sstart=1, send=900, bitscore=1000)],
            columns=COMPARISON_COLUMNS,
        ),
        (1, 0): pd.DataFrame.from_records(
            [_hit_row("b0", "a0", alignment_length=900, qstart=1, qend=900, sstart=1, send=900, bitscore=1000)],
            columns=COMPARISON_COLUMNS,
        ),
        (2, 0): pd.DataFrame.from_records(
            [_hit_row("c0", "a0", alignment_length=200, qstart=1, qend=200, sstart=1, send=200, bitscore=400)],
            columns=COMPARISON_COLUMNS,
        ),
        (2, 1): pd.DataFrame.from_records(
            [_hit_row("c0", "b0", alignment_length=400, qstart=1, qend=400, sstart=1, send=400, bitscore=350)],
            columns=COMPARISON_COLUMNS,
        ),
    }

    result = select_rbh_orthogroup_edges_from_directional_hits(
        directional_hits,
        extraction.protein_map,
        record_count=len(records),
        orthogroup_membership_mode="anchor_core_v1",
        max_related_edges_per_orthogroup=2,
    ).orthogroups

    member = result.member_by_protein_id["c0"]
    assert member.orthogroup_id == result.member_by_protein_id["a0"].orthogroup_id
    assert member.confidence == "high"
    assert any(
        edge.edge_kind == "domain_only"
        for edges in result.related_edges_by_orthogroup_id.values()
        for edge in edges
    )


@pytest.mark.linear
def test_anchor_core_keeps_true_domain_only_multi_hsp_hit_unassigned() -> None:
    records = [
        _long_record("record_a", "a0"),
        _long_record("record_b", "b0"),
        _long_record("record_c", "c0"),
    ]
    extraction = extract_cds_proteins(records)
    directional_hits = {
        (0, 1): pd.DataFrame.from_records(
            [_hit_row("a0", "b0", alignment_length=900, qstart=1, qend=900, sstart=1, send=900, bitscore=1000)],
            columns=COMPARISON_COLUMNS,
        ),
        (1, 0): pd.DataFrame.from_records(
            [_hit_row("b0", "a0", alignment_length=900, qstart=1, qend=900, sstart=1, send=900, bitscore=1000)],
            columns=COMPARISON_COLUMNS,
        ),
        (2, 1): pd.DataFrame.from_records(
            [
                _hit_row("c0", "b0", alignment_length=200, qstart=1, qend=200, sstart=1, send=200, bitscore=320),
                _hit_row("c0", "b0", alignment_length=191, qstart=50, qend=240, sstart=50, send=240, bitscore=300),
            ],
            columns=COMPARISON_COLUMNS,
        ),
    }

    result = select_rbh_orthogroup_edges_from_directional_hits(
        directional_hits,
        extraction.protein_map,
        record_count=len(records),
        orthogroup_membership_mode="anchor_core_v1",
        max_related_edges_per_orthogroup=2,
    ).orthogroups

    assert "c0" not in result.member_by_protein_id
    assert any(
        edge.edge_kind == "domain_only"
        for edges in result.related_edges_by_orthogroup_id.values()
        for edge in edges
    )


@pytest.mark.linear
def test_build_rbh_orthogroup_protein_blastp_comparisons_keeps_transitive_all_vs_all_grouping() -> None:
    records = [
        _record("record_a", features=[_cds(0, 9)]),
        _record("record_b", features=[_cds(9, 18)]),
        _record("record_c", features=[_cds(18, 27)]),
    ]
    calls: list[tuple[str, str]] = []
    rows_by_call = [
        [],
        [_hit_row("gbd_r0001_cds000001", "gbd_r0002_cds000001")],
        [_hit_row("gbd_r0002_cds000001", "gbd_r0001_cds000001")],
        [],
        [],
        [],
        [_hit_row("gbd_r0002_cds000001", "gbd_r0003_cds000001")],
        [_hit_row("gbd_r0003_cds000001", "gbd_r0002_cds000001")],
        [],
    ]

    def fake_runner(query_fasta: str, subject_fasta: str) -> pd.DataFrame:
        calls.append((query_fasta, subject_fasta))
        rows = rows_by_call[len(calls) - 1]
        return pd.DataFrame.from_records(rows, columns=COMPARISON_COLUMNS)

    result = build_rbh_orthogroup_protein_blastp_comparisons(
        records,
        runner=fake_runner,
        identity=0,
    )

    comparisons = result.comparisons
    assert len(calls) == 9
    assert result.orthogroups is not None
    assert set(result.orthogroups.member_by_protein_id) == {
        "gbd_r0001_cds000001",
        "gbd_r0002_cds000001",
        "gbd_r0003_cds000001",
    }
    assert comparisons[0].iloc[0]["orthogroup_id"] == "og_1"
    assert comparisons[1].iloc[0]["orthogroup_id"] == "og_1"


@pytest.mark.linear
def test_orthogroup_expanded_display_edges_include_non_rbh_coorthologs() -> None:
    records = [
        _record(
            "record_a",
            features=[
                _cds(0, 9, qualifiers={"translation": ["MK*"], "protein_id": ["a0"]}),
                _cds(12, 21, qualifiers={"translation": ["MK*"], "protein_id": ["a1"]}),
            ],
        ),
        _record(
            "record_b",
            features=[
                _cds(0, 9, qualifiers={"translation": ["MK*"], "protein_id": ["b0"]}),
                _cds(12, 21, qualifiers={"translation": ["MK*"], "protein_id": ["b1"]}),
            ],
        ),
    ]
    extraction = extract_cds_proteins(records)
    directional_hits = {
        (0, 1): pd.DataFrame.from_records(
            [
                _hit_row("a0", "b0", bitscore=300),
                _hit_row("a1", "b0", bitscore=250),
            ],
            columns=COMPARISON_COLUMNS,
        ),
        (1, 0): pd.DataFrame.from_records(
            [
                _hit_row("b0", "a0", bitscore=300),
                _hit_row("b1", "a0", bitscore=260),
            ],
            columns=COMPARISON_COLUMNS,
        ),
    }

    edge_selection = select_rbh_orthogroup_edges_from_directional_hits(
        directional_hits,
        extraction.protein_map,
        record_count=len(records),
        orthogroup_membership_mode="family_merge",
        orthogroup_member_max_hits=2,
        max_related_edges_per_orthogroup=2,
    )

    anchor_edges = edge_selection.adjacent_anchor_edges_by_pair[(0, 1)]
    display_edges = edge_selection.adjacent_display_edges_by_pair[(0, 1)]
    assert {
        (str(row.query), str(row.subject))
        for row in anchor_edges.itertuples(index=False)
    } == {("a0", "b0")}
    assert {
        (str(row.query), str(row.subject))
        for row in display_edges.itertuples(index=False)
    } == {
        ("a0", "b0"),
        ("a0", "b1"),
        ("a1", "b0"),
    }

    converted = convert_protein_hits_to_genomic_links(
        display_edges,
        extraction.protein_map,
        orthogroups=edge_selection.orthogroups,
    )
    edge_kind_by_pair = {
        (str(row.query_protein_id), str(row.subject_protein_id)): str(row.edge_kind)
        for row in converted.itertuples(index=False)
    }
    assert edge_kind_by_pair[("a0", "b0")] == "rbh"
    assert edge_kind_by_pair[("a0", "b1")] == "coortholog"
    assert edge_kind_by_pair[("a1", "b0")] == "coortholog"


@pytest.mark.linear
def test_anchor_core_separates_weakly_bridged_outparalog_families_for_legacy_aliases() -> None:
    records = [
        _record(
            "record_a",
            sequence="ATG" * 260,
            features=[
                _cds(0, 300, qualifiers={"translation": ["M" * 100 + "*"], "protein_id": ["d0"], "gene": ["dnaE"]}),
                _cds(360, 660, qualifiers={"translation": ["M" * 100 + "*"], "protein_id": ["p0"], "gene": ["polC"]}),
            ],
        ),
        _record(
            "record_b",
            sequence="ATG" * 260,
            features=[
                _cds(0, 300, qualifiers={"translation": ["M" * 100 + "*"], "protein_id": ["d1"], "gene": ["dnaE"]}),
                _cds(360, 660, qualifiers={"translation": ["M" * 100 + "*"], "protein_id": ["p1"], "gene": ["polC"]}),
            ],
        ),
        _record(
            "record_c",
            sequence="ATG" * 260,
            features=[
                _cds(0, 300, qualifiers={"translation": ["M" * 100 + "*"], "protein_id": ["d2"], "gene": ["dnaE"]}),
                _cds(360, 660, qualifiers={"translation": ["M" * 100 + "*"], "protein_id": ["p2"], "gene": ["polC"]}),
            ],
        ),
    ]
    extraction = extract_cds_proteins(records)

    directional_hits = {
        (0, 1): pd.DataFrame.from_records(
            [
                _hit_row("d0", "d1", bitscore=320),
                _hit_row("p0", "p1", bitscore=330),
                _hit_row("d0", "p1", bitscore=150),
            ],
            columns=COMPARISON_COLUMNS,
        ),
        (1, 0): pd.DataFrame.from_records(
            [
                _hit_row("d1", "d0", bitscore=320),
                _hit_row("p1", "p0", bitscore=330),
                _hit_row("p1", "d0", bitscore=150),
            ],
            columns=COMPARISON_COLUMNS,
        ),
        (0, 2): pd.DataFrame.from_records(
            [
                _hit_row("d0", "d2", bitscore=315),
                _hit_row("p0", "p2", bitscore=325),
            ],
            columns=COMPARISON_COLUMNS,
        ),
        (2, 0): pd.DataFrame.from_records(
            [
                _hit_row("d2", "d0", bitscore=315),
                _hit_row("p2", "p0", bitscore=325),
            ],
            columns=COMPARISON_COLUMNS,
        ),
        (1, 2): pd.DataFrame.from_records(
            [
                _hit_row("d1", "d2", bitscore=318),
                _hit_row("p1", "p2", bitscore=328),
            ],
            columns=COMPARISON_COLUMNS,
        ),
        (2, 1): pd.DataFrame.from_records(
            [
                _hit_row("d2", "d1", bitscore=318),
                _hit_row("p2", "p1", bitscore=328),
            ],
            columns=COMPARISON_COLUMNS,
        ),
    }

    merged_selection = select_rbh_orthogroup_edges_from_directional_hits(
        directional_hits,
        extraction.protein_map,
        record_count=len(records),
        orthogroup_membership_mode="family_merge",
        orthogroup_member_max_hits=3,
        max_related_edges_per_orthogroup=2,
    )
    merged_groups = {
        orthogroup_id: {member.protein_id for member in members}
        for orthogroup_id, members in merged_selection.orthogroups.orthogroups.items()
    }
    assert merged_groups == {
        "og_1": {"d0", "d1", "d2"},
        "og_2": {"p0", "p1", "p2"},
    }

    split_selection = select_rbh_orthogroup_edges_from_directional_hits(
        directional_hits,
        extraction.protein_map,
        record_count=len(records),
        orthogroup_membership_mode="distribution_split",
        orthogroup_member_max_hits=3,
        max_related_edges_per_orthogroup=2,
    )

    groups = {
        orthogroup_id: {member.protein_id for member in members}
        for orthogroup_id, members in split_selection.orthogroups.orthogroups.items()
    }
    assert groups == {
        "og_1": {"d0", "d1", "d2"},
        "og_2": {"p0", "p1", "p2"},
    }
    assert split_selection.orthogroups.member_by_protein_id["d0"].orthogroup_id != (
        split_selection.orthogroups.member_by_protein_id["p1"].orthogroup_id
    )


@pytest.mark.linear
def test_anchor_core_display_edges_include_direct_adjacent_assigned_same_orthogroup_hits() -> None:
    records = [
        _record(
            "record_a",
            sequence="ATG" * 120,
            features=[
                _cds(0, 300, qualifiers={"translation": ["M" * 100 + "*"], "protein_id": ["a0"]}),
            ],
        ),
        _record(
            "record_b",
            sequence="ATG" * 120,
            features=[
                _cds(0, 300, qualifiers={"translation": ["M" * 100 + "*"], "protein_id": ["b0"]}),
            ],
        ),
        _record(
            "record_c",
            sequence="ATG" * 120,
            features=[
                _cds(0, 300, qualifiers={"translation": ["M" * 100 + "*"], "protein_id": ["c0"]}),
            ],
        ),
    ]
    extraction = extract_cds_proteins(records)
    directional_hits = {
        (0, 1): pd.DataFrame.from_records(
            [
                _hit_row("a0", "b0", bitscore=1000),
            ],
            columns=COMPARISON_COLUMNS,
        ),
        (1, 2): pd.DataFrame.from_records(
            [
                _hit_row("b0", "c0", bitscore=900),
            ],
            columns=COMPARISON_COLUMNS,
        ),
        (2, 1): pd.DataFrame.from_records(
            [
                _hit_row("c0", "b0", bitscore=900),
            ],
            columns=COMPARISON_COLUMNS,
        ),
    }

    edge_selection = select_rbh_orthogroup_edges_from_directional_hits(
        directional_hits,
        extraction.protein_map,
        record_count=len(records),
        orthogroup_membership_mode="distribution_split",
        orthogroup_member_max_hits=1,
        max_related_edges_per_orthogroup=2,
    )

    assert edge_selection.adjacent_anchor_edges_by_pair[(0, 1)].empty
    assert {
        frozenset(member.protein_id for member in members)
        for members in edge_selection.orthogroups.orthogroups.values()
    } == {frozenset({"a0", "b0", "c0"})}
    assert {
        (edge.query_protein_id, edge.subject_protein_id)
        for edges in edge_selection.orthogroups.ortholog_edges_by_orthogroup_id.values()
        for edge in edges
    } == {("a0", "b0"), ("b0", "c0")}

    display_edges = edge_selection.adjacent_display_edges_by_pair[(0, 1)]
    assert {
        (str(row.query), str(row.subject))
        for row in display_edges.itertuples(index=False)
    } == {("a0", "b0")}

    converted = convert_protein_hits_to_genomic_links(
        display_edges,
        extraction.protein_map,
        orthogroups=edge_selection.orthogroups,
    )
    assert converted.iloc[0]["orthogroup_id"] == "og_1"


@pytest.mark.linear
def test_family_merge_display_edges_suppress_already_covered_cross_links() -> None:
    records = [
        _record(
            "record_a",
            features=[
                _cds(0, 9, qualifiers={"translation": ["MK*"], "protein_id": ["a0"]}),
                _cds(12, 21, qualifiers={"translation": ["MK*"], "protein_id": ["a1"]}),
            ],
        ),
        _record(
            "record_b",
            features=[
                _cds(0, 9, qualifiers={"translation": ["MK*"], "protein_id": ["b0"]}),
                _cds(12, 21, qualifiers={"translation": ["MK*"], "protein_id": ["b1"]}),
            ],
        ),
    ]
    extraction = extract_cds_proteins(records)
    directional_hits = {
        (0, 1): pd.DataFrame.from_records(
            [
                _hit_row("a0", "b0", bitscore=300),
                _hit_row("a1", "b1", bitscore=310),
                _hit_row("a0", "b1", bitscore=280),
            ],
            columns=COMPARISON_COLUMNS,
        ),
        (1, 0): pd.DataFrame.from_records(
            [
                _hit_row("b0", "a0", bitscore=300),
                _hit_row("b1", "a1", bitscore=310),
                _hit_row("b1", "a0", bitscore=280),
            ],
            columns=COMPARISON_COLUMNS,
        ),
    }

    edge_selection = select_rbh_orthogroup_edges_from_directional_hits(
        directional_hits,
        extraction.protein_map,
        record_count=len(records),
        orthogroup_membership_mode="family_merge",
        orthogroup_member_max_hits=2,
        max_related_edges_per_orthogroup=2,
    )

    display_edges = edge_selection.adjacent_display_edges_by_pair[(0, 1)]
    assert {
        (str(row.query), str(row.subject))
        for row in display_edges.itertuples(index=False)
    } == {
        ("a0", "b0"),
        ("a0", "b1"),
        ("a1", "b1"),
    }
    assert {
        (edge.query_protein_id, edge.subject_protein_id, edge.edge_kind)
        for edge in edge_selection.orthogroups.ortholog_edges_by_orthogroup_id["og_1"]
    } == {
        ("a0", "b0", "rbh"),
        ("a0", "b1", "coortholog"),
        ("a1", "b1", "rbh"),
    }
    assert edge_selection.orthogroups.related_edges_by_orthogroup_id.get("og_1", ()) == ()


@pytest.mark.linear
def test_family_merge_display_edges_prefer_uncovered_alternative_links() -> None:
    records = [
        _record(
            "record_a",
            features=[
                _cds(0, 9, qualifiers={"translation": ["MK*"], "protein_id": ["a0"]}),
                _cds(12, 21, qualifiers={"translation": ["MK*"], "protein_id": ["a1"]}),
            ],
        ),
        _record(
            "record_b",
            features=[
                _cds(0, 9, qualifiers={"translation": ["MK*"], "protein_id": ["b0"]}),
                _cds(12, 21, qualifiers={"translation": ["MK*"], "protein_id": ["b1"]}),
            ],
        ),
    ]
    extraction = extract_cds_proteins(records)
    directional_hits = {
        (0, 1): pd.DataFrame.from_records(
            [
                _hit_row("a1", "b1", bitscore=400),
                _hit_row("a0", "b1", bitscore=390),
                _hit_row("a0", "b0", bitscore=200),
            ],
            columns=COMPARISON_COLUMNS,
        ),
        (1, 0): pd.DataFrame.from_records(
            [
                _hit_row("b1", "a1", bitscore=400),
                _hit_row("b1", "a0", bitscore=390),
                _hit_row("b0", "a0", bitscore=200),
            ],
            columns=COMPARISON_COLUMNS,
        ),
    }

    edge_selection = select_rbh_orthogroup_edges_from_directional_hits(
        directional_hits,
        extraction.protein_map,
        record_count=len(records),
        orthogroup_membership_mode="family_merge",
        orthogroup_member_max_hits=2,
        max_related_edges_per_orthogroup=2,
    )

    display_edges = edge_selection.adjacent_display_edges_by_pair[(0, 1)]
    assert {
        (str(row.query), str(row.subject))
        for row in display_edges.itertuples(index=False)
    } == {
        ("a1", "b1"),
        ("a0", "b1"),
        ("a0", "b0"),
    }
    assert {
        (edge.query_protein_id, edge.subject_protein_id)
        for edge in edge_selection.orthogroups.ortholog_edges_by_orthogroup_id["og_1"]
    } >= {
        ("a0", "b1"),
        ("a0", "b0"),
    }


@pytest.mark.linear
def test_convert_protein_hits_to_genomic_links_only_sets_matching_orthogroup_id() -> None:
    records = [
        _record("record_a", features=[_cds(0, 9)]),
        _record("record_b", features=[_cds(9, 18)]),
        _record("record_c", features=[_cds(18, 27)]),
    ]
    extraction = extract_cds_proteins(records)
    directional_hits = {
        (0, 1): pd.DataFrame.from_records(
            [_hit_row("gbd_r0001_cds000001", "gbd_r0002_cds000001")],
            columns=COMPARISON_COLUMNS,
        ),
        (1, 0): pd.DataFrame.from_records(
            [_hit_row("gbd_r0002_cds000001", "gbd_r0001_cds000001")],
            columns=COMPARISON_COLUMNS,
        ),
    }
    display_hits = pd.DataFrame.from_records(
        [
            _hit_row("gbd_r0001_cds000001", "gbd_r0002_cds000001"),
            _hit_row("gbd_r0001_cds000001", "gbd_r0003_cds000001"),
        ],
        columns=COMPARISON_COLUMNS,
    )
    orthogroups = select_rbh_orthogroup_edges_from_directional_hits(
        directional_hits,
        extraction.protein_map,
        record_count=len(records),
    ).orthogroups

    converted = convert_protein_hits_to_genomic_links(
        display_hits,
        extraction.protein_map,
        orthogroups=orthogroups,
    )

    assert converted.iloc[0]["orthogroup_id"] == "og_1"
    assert converted.iloc[1]["orthogroup_id"] == ""


@pytest.mark.linear
def test_web_losat_nucleotide_display_transform_keeps_orientation() -> None:
    namespace = _load_web_helper_namespace()
    raw_tsv = "\n".join(
        [
            "\t".join(
                [
                    "query",
                    "subject",
                    "99.0",
                    "40",
                    "0",
                    "0",
                    "10",
                    "20",
                    "5",
                    "15",
                    "1e-20",
                    "120",
                ]
            ),
            "\t".join(
                [
                    "query",
                    "subject",
                    "98.0",
                    "35",
                    "0",
                    "0",
                    "30",
                    "25",
                    "40",
                    "35",
                    "1e-10",
                    "100",
                ]
            ),
        ]
    )

    raw_result = namespace["convert_losat_nucleotide_to_display_tsv"](
        raw_tsv,
        json.dumps({"length": 100, "reverse": True}),
        json.dumps({"length": 80, "reverse": False}),
    )
    result = json.loads(str(raw_result))

    assert "error" not in result
    row = result["rows"][0]
    assert row["qstart"] == 91
    assert row["qend"] == 81
    assert row["sstart"] == 5
    assert row["send"] == 15
    assert len(result["rows"]) == 2


@pytest.mark.linear
def test_web_cds_span_transform_maps_reverse_display_span_and_strand() -> None:
    namespace = _load_web_helper_namespace()

    assert namespace["_web_transform_cds_span"](10, 40, 1, {"length": 100, "reverse": True}) == (60, 90, -1)
    assert namespace["_web_transform_cds_span"](10, 40, -1, {"length": 100, "reverse": True}) == (60, 90, 1)


@pytest.mark.linear
def test_web_extract_cds_protein_fasta_uses_coordinate_stable_ids(tmp_path: Path) -> None:
    namespace = _load_web_helper_namespace()
    record = _record(
        "record_a",
        features=[
            _cds(
                0,
                9,
                qualifiers={
                    "translation": ["MKT*"],
                    "locus_tag": ["gene_a"],
                },
            )
        ],
    )
    record.annotations["molecule_type"] = "DNA"
    gb_path = tmp_path / "input.gb"
    SeqIO.write([record], gb_path, "genbank")

    raw_result = namespace["extract_cds_protein_fasta"](
        str(gb_path),
        "genbank",
        None,
        None,
        None,
        "0",
        0,
        "record_a_region",
    )
    result = json.loads(str(raw_result))

    assert "error" not in result
    protein_id = next(iter(result["protein_map"]))
    assert protein_id.startswith("h_")
    assert len(protein_id) == 28
    assert "gbd_r0001_cds000001" not in result["fasta"]
    assert result["identity_manifest"]["schema"] == 2
    assert result["protein_set_hash"].startswith("sha256:")
    assert result["record_analysis_id"].startswith("sha256:")
    assert result["runtime_binding_hash"].startswith("sha256:")
    assert result["display_binding_hash"].startswith("sha256:")
    boundary_key = json.loads(
        str(
            namespace["build_protein_losat_cache_keys_json"](
                json.dumps(result["identity_manifest"]),
                json.dumps([{
                    "queryRecordInstanceKey": "record_a_region",
                    "subjectRecordInstanceKey": "record_a_region",
                    "expectedOptions": {"program": "blastp", "outfmt": "6", "args": []},
                }]),
            )
        )
    )["keys"][0]
    expected_identity = build_protein_losat_pair_identity(
        result["identity_manifest"],
        query_record_instance_key="record_a_region",
        subject_record_instance_key="record_a_region",
    )
    assert boundary_key == build_protein_losat_cache_key(expected_identity, args=[])


@pytest.mark.linear
def test_web_protein_cache_keys_validate_once_and_preserve_direction_and_options(monkeypatch) -> None:
    namespace = _load_web_helper_namespace()
    extraction = extract_protein_identity_manifest(
        [_record("a", features=[_cds(0, 9, qualifiers={"translation": ["MKT"]})]),
         _record("b", features=[_cds(0, 9, qualifiers={"translation": ["MGG"]})])],
        record_instance_keys=("row-a", "row-b"),
    )
    manifest = extraction.identity_manifest.to_dict()
    pairs = [
        {"queryRecordInstanceKey": query, "subjectRecordInstanceKey": subject,
         "expectedOptions": {"program": "blastp", "outfmt": "6", "args": args}}
        for query, subject, args in [
            ("row-a", "row-b", []), ("row-b", "row-a", []),
            ("row-a", "row-a", []), ("row-a", "row-b", ["--max-target-seqs", "6"]),
            ("row-a", "row-b", []),
        ]
    ]
    pairs.append({**pairs[0], "expectedOptions": {**pairs[0]["expectedOptions"], "searchContext": "a" * 64}})
    expected = [
        build_protein_losat_cache_key(
            build_protein_losat_pair_identity(
                manifest, query_record_instance_key=pair["queryRecordInstanceKey"],
                subject_record_instance_key=pair["subjectRecordInstanceKey"],
            ), args=pair["expectedOptions"]["args"],
            search_context=pair["expectedOptions"].get("searchContext"),
        ) for pair in pairs
    ]
    validation_calls = []

    def validate_once(value):
        validation_calls.append(1)
        return validate_protein_identity_manifest(value)

    monkeypatch.setattr(protein_colinearity_module, "validate_protein_identity_manifest", validate_once)
    result = json.loads(namespace["build_protein_losat_cache_keys_json"](
        json.dumps(manifest), json.dumps(pairs),
    ))
    assert result == {"keys": expected}
    assert len(set(expected)) == 5
    assert validation_calls == [1]
    # A bad later pair must not expose a partial key list.
    pairs[-1]["subjectRecordInstanceKey"] = "unknown"
    rejected = json.loads(namespace["build_protein_losat_cache_keys_json"](
        json.dumps(manifest), json.dumps(pairs),
    ))
    assert "error" in rejected and "keys" not in rejected
    # Validate the full manifest even when a bad record is not in a requested pair.
    pairs = pairs[2:3]
    manifest["recordInstances"]["row-b"]["runtimeBindingHash"] = "invalid"
    rejected = json.loads(namespace["build_protein_losat_cache_keys_json"](
        json.dumps(manifest), json.dumps(pairs),
    ))
    assert "error" in rejected and "keys" not in rejected


@pytest.mark.linear
def test_web_losatp_pairwise_payload_uses_display_view_transform(
    tmp_path: Path,
    stage_web_losatp_transport,
) -> None:
    namespace = _load_web_helper_namespace()
    hits = pd.DataFrame.from_records(
        [_hit_row("qa", "sb")],
        columns=COMPARISON_COLUMNS,
    )
    payload = {
        "records": [
            {
                "recordIndex": 0,
                "recordId": "record_a",
                "proteinMap": {
                    "qa": _web_protein_entry(
                        "qa",
                        record_index=0,
                        record_id="record_a",
                        start=0,
                        end=30,
                        strand=1,
                    )
                },
                "proteinCacheKey": "record-a-cache",
                "viewTransform": {"length": 300, "reverse": True},
            },
            {
                "recordIndex": 1,
                "recordId": "record_b",
                "proteinMap": {
                    "sb": _web_protein_entry(
                        "sb",
                        record_index=1,
                        record_id="record_b",
                        start=100,
                        end=160,
                        strand=-1,
                    )
                },
                "proteinCacheKey": "record-b-cache",
                "viewTransform": {"length": 200, "reverse": False},
            },
        ],
        "pairs": [
            {
                "pairIndex": 0,
                "queryIndex": 0,
                "subjectIndex": 1,
                "cacheKey": "pair-a-b",
                "blastText": hits.to_csv(sep="\t", header=False, index=False, lineterminator="\n"),
            }
        ],
    }

    pairs_path, raw_tsv_path = stage_web_losatp_transport(tmp_path, payload)
    raw_result = namespace["convert_losatp_blastp_pairs_to_genomic_payload"](
        str(pairs_path),
        str(raw_tsv_path),
        "pairwise",
        1,
        50,
        "1e-5",
        0,
        0,
    )
    result = json.loads(str(raw_result))

    assert "error" not in result
    row = result["pairs"][0]["rows"][0]
    assert row["query_protein_id"] == "qa"
    assert row["qstart"] == 300
    assert row["qend"] == 271
    assert row["sstart"] == 160
    assert row["send"] == 101


@pytest.mark.linear
def test_web_losatp_blastp_payload_helper_uses_rbh_edges_for_orthogroups(
    tmp_path: Path,
    stage_web_losatp_transport,
) -> None:
    helpers_js = Path("gbdraw/web/js/app/python-helpers.js").read_text(encoding="utf-8")
    helper_source = helpers_js.split("`", 1)[1].rsplit("`", 1)[0]
    namespace: dict[str, object] = {}
    exec(helper_source, namespace)

    forward_hits = pd.DataFrame.from_records(
        [
            _hit_row("a1", "b", bitscore=300),
            _hit_row("a2", "b", bitscore=250),
        ],
        columns=COMPARISON_COLUMNS,
    )
    reverse_hits = pd.DataFrame.from_records(
        [_hit_row("b", "a1", bitscore=300)],
        columns=COMPARISON_COLUMNS,
    )
    query_map = {
        "a1": _web_protein_entry(
            "a1",
            record_index=0,
            record_id="record_a",
            feature_index=0,
            gene="rpoB",
            product="DNA-directed RNA polymerase beta subunit",
        ),
        "a2": _web_protein_entry(
            "a2",
            record_index=0,
            record_id="record_a",
            feature_index=1,
            start=100,
            end=190,
        ),
    }
    subject_map = {
        "b": _web_protein_entry(
            "b",
            record_index=1,
            record_id="record_b",
            gene="rpoB",
            product="DNA-directed RNA polymerase beta subunit",
        ),
    }
    payload = {
        "records": [
            {
                "recordIndex": 0,
                "recordId": "record_a",
                "proteinMap": query_map,
                "proteinCacheKey": "record-a-cache",
                "viewTransform": {"length": 200, "reverse": False},
            },
            {
                "recordIndex": 1,
                "recordId": "record_b",
                "proteinMap": subject_map,
                "proteinCacheKey": "record-b-cache",
                "viewTransform": {"length": 200, "reverse": True},
            },
        ],
        "pairs": [
            {
                "pairIndex": 0,
                "queryIndex": 0,
                "displayPair": True,
                "subjectIndex": 1,
                "cacheKey": "pair-a-b",
                "blastText": forward_hits.to_csv(
                    sep="\t",
                    header=False,
                    index=False,
                    lineterminator="\n",
                ),
            },
            {
                "pairIndex": 0,
                "queryIndex": 1,
                "subjectIndex": 0,
                "cacheKey": "pair-b-a",
                "blastText": reverse_hits.to_csv(
                    sep="\t",
                    header=False,
                    index=False,
                    lineterminator="\n",
                ),
            },
        ],
    }

    pairs_path, raw_tsv_path = stage_web_losatp_transport(tmp_path, payload)
    raw_result = namespace["convert_losatp_blastp_pairs_to_genomic_payload"](
        str(pairs_path),
        str(raw_tsv_path),
        "orthogroup",
        2,
        50,
        "1e-5",
        0,
        0,
        orthogroup_membership_mode="rbh",
    )
    result = json.loads(str(raw_result))

    assert "error" not in result
    assert result["orthogroupResult"]["schema"] == 3
    assert result["orthogroupResult"]["kind"] == "orthogroupResult"
    assert result["orthogroupResult"]["value"]["type"] == "OrthogroupGraphResult"
    typed_fields = result["orthogroupResult"]["value"]["fields"]
    group_id = next(iter(typed_fields["orthogroups"]))
    group_members = typed_fields["orthogroups"][group_id]
    display_start, display_end, display_strand = namespace["_web_transform_cds_span"](
        subject_map["b"]["start"],
        subject_map["b"]["end"],
        subject_map["b"]["strand"],
        payload["records"][1]["viewTransform"],
    )
    display_feature_svg_id = namespace["_display_feature_svg_id_from_data"](
        subject_map["b"],
        display_start,
        display_end,
        display_strand,
        payload["records"][1]["viewTransform"],
    )
    assert display_feature_svg_id != subject_map["b"]["feature_svg_id"]
    assert len(group_members) == 3
    assert typed_fields["namesByOrthogroupId"][group_id] == "rpoB"
    assert typed_fields["confidenceByOrthogroupId"][group_id] == "high"
    first_candidate = typed_fields["nameCandidatesByOrthogroupId"][group_id][0][
        "fields"
    ]
    assert first_candidate["recordCoverageCount"] == 2
    assert first_candidate["source"] == "gene"
    assert group_members[0]["fields"]["product"] == (
        "DNA-directed RNA polymerase beta subunit"
    )
    subject_member = next(
        member["fields"]
        for member in group_members
        if member["fields"]["proteinId"] == "b"
    )
    assert subject_member["featureSvgId"] == "feature_b"
    assert "orthogroups" not in result
    rows = result["pairs"][0]["rows"]
    assert rows[0]["subject_feature_svg_id"] == "feature_b"
    assert rows[0]["subject_view_feature_svg_id"] == display_feature_svg_id
    assert rows[0]["orthogroup_id"] == "og_1"
    assert rows[0]["edge_kind"] == "rbh"
    assert rows[1]["orthogroup_id"] == "og_1"
    assert rows[1]["edge_kind"] == "coortholog"
    assert len(rows) == 2
    assert result["cache"]["convertedPayloadHit"] is False
    assert result["cache"]["filteredHitCacheMisses"] == 2
    assert result["cache"]["rawTsvEntryCount"] == 2
    assert result["cache"]["rawTsvBytes"] == raw_tsv_path.stat().st_size
    assert 0 < result["cache"]["rawTsvLargestEntryBytes"] <= raw_tsv_path.stat().st_size
    assert result["cache"]["simultaneousParsedTables"] == 2

    repeated_result = json.loads(str(namespace["convert_losatp_blastp_pairs_to_genomic_payload"](
        str(pairs_path),
        str(raw_tsv_path),
        "orthogroup",
        2,
        50,
        "1e-5",
        0,
        0,
        orthogroup_membership_mode="rbh",
    )))
    assert repeated_result["cache"]["convertedPayloadHit"] is True
    assert repeated_result["cache"]["simultaneousParsedTables"] == 0
    assert repeated_result["pairs"] == result["pairs"]

    inactive_pairwise_limit_result = json.loads(str(namespace["convert_losatp_blastp_pairs_to_genomic_payload"](
        str(pairs_path),
        str(raw_tsv_path),
        "orthogroup",
        3,
        50,
        "1e-5",
        0,
        0,
        orthogroup_membership_mode="rbh",
    )))
    assert inactive_pairwise_limit_result["cache"]["convertedPayloadHit"] is True
    assert inactive_pairwise_limit_result["pairs"] == result["pairs"]


@pytest.mark.linear
def test_web_losatp_blastp_payload_helper_rejects_legacy_list_payload(
    tmp_path: Path,
) -> None:
    helpers_js = Path("gbdraw/web/js/app/python-helpers.js").read_text(encoding="utf-8")
    helper_source = helpers_js.split("`", 1)[1].rsplit("`", 1)[0]
    namespace: dict[str, object] = {}
    exec(helper_source, namespace)

    pairs_path = tmp_path / "losatp-pairs.json"
    raw_tsv_path = tmp_path / "losatp-pairs.tsv"
    pairs_path.write_text(json.dumps([]), encoding="utf-8")
    raw_tsv_path.write_bytes(b"")
    raw_result = namespace["convert_losatp_blastp_pairs_to_genomic_payload"](
        str(pairs_path),
        str(raw_tsv_path),
        "pairwise",
        2,
        50,
        "1e-5",
        0,
        0,
    )
    result = json.loads(str(raw_result))

    assert "error" in result
    assert "must be an object" in result["error"]


@pytest.mark.linear
def test_orthogroup_alignment_offsets_align_selected_member_to_representatives() -> None:
    records = [
        _record("record_a", sequence="A" * 1000),
        _record("record_b", sequence="A" * 1000),
    ]
    comparison = pd.DataFrame.from_records(
        [
            {
                **_hit_row("record_a", "record_b", bitscore=200),
                "qstart": 100,
                "qend": 200,
                "sstart": 400,
                "send": 500,
                "query_protein_id": "prot_a",
                "subject_protein_id": "prot_b",
                "query_source_protein_id": "",
                "subject_source_protein_id": "",
                "query_record_index": 0,
                "subject_record_index": 1,
                "query_feature_index": 0,
                "subject_feature_index": 0,
                "query_feature_svg_id": "fanchor",
                "subject_feature_svg_id": "fsubject",
                "orthogroup_id": "og_1",
                "query_orthogroup_representative": True,
                "subject_orthogroup_representative": True,
            }
        ]
    )
    canvas_config = _orthogroup_alignment_canvas_config()

    offsets = calculate_orthogroup_alignment_offsets(
        records,
        [comparison],
        canvas_config,
        "fanchor",
    )

    assert offsets[0] == pytest.approx(0.0)
    assert offsets[1] == pytest.approx(-300.0)


@pytest.mark.linear
def test_orthogroup_alignment_dedup_ignores_public_source_protein_id() -> None:
    records = [
        _record("record_a", sequence="A" * 1000),
        _record("record_b", sequence="A" * 1000),
    ]
    rows = []
    for source_protein_id, bitscore in (("public-a", 200), ("public-b", 150)):
        rows.append(
            {
                **_hit_row("record_a", "record_b", bitscore=bitscore),
                "qstart": 100,
                "qend": 200,
                "sstart": 400,
                "send": 500,
                "query_protein_id": f"runtime-{source_protein_id}",
                "subject_protein_id": "runtime-b",
                "query_source_protein_id": source_protein_id,
                "subject_source_protein_id": "public-subject",
                "query_record_index": 0,
                "subject_record_index": 1,
                "query_feature_index": 4,
                "subject_feature_index": 7,
                "query_feature_svg_id": "fanchor",
                "subject_feature_svg_id": "fsubject",
                "orthogroup_id": "og_1",
                "query_orthogroup_representative": True,
                "subject_orthogroup_representative": True,
            }
        )

    offsets = calculate_orthogroup_alignment_offsets(
        records,
        [pd.DataFrame.from_records(rows)],
        _orthogroup_alignment_canvas_config(),
        "fanchor",
    )

    assert offsets[0] == pytest.approx(0.0)
    assert offsets[1] == pytest.approx(-300.0)


@pytest.mark.linear
def test_orthogroup_alignment_rejects_conflicting_group_for_one_feature() -> None:
    records = [
        _record("record_a", sequence="A" * 1000),
        _record("record_b", sequence="A" * 1000),
    ]
    rows = []
    for orthogroup_id in ("og_1", "og_2"):
        rows.append(
            {
                **_hit_row("record_a", "record_b", bitscore=200),
                "qstart": 100,
                "qend": 200,
                "sstart": 400,
                "send": 500,
                "query_protein_id": "runtime-a",
                "subject_protein_id": f"runtime-{orthogroup_id}",
                "query_source_protein_id": "public-a",
                "subject_source_protein_id": f"public-{orthogroup_id}",
                "query_record_index": 0,
                "subject_record_index": 1,
                "query_feature_index": 4,
                "subject_feature_index": 7 if orthogroup_id == "og_1" else 8,
                "query_feature_svg_id": "fanchor",
                "subject_feature_svg_id": f"fsubject-{orthogroup_id}",
                "orthogroup_id": orthogroup_id,
                "query_orthogroup_representative": True,
                "subject_orthogroup_representative": True,
            }
        )

    with pytest.raises(ValidationError, match="conflicting orthogroups"):
        calculate_orthogroup_alignment_offsets(
            records,
            [pd.DataFrame.from_records(rows)],
            _orthogroup_alignment_canvas_config(),
            "fanchor",
        )


@pytest.mark.linear
def test_orthogroup_alignment_canvas_adjustment_fits_negative_record_offsets() -> None:
    records = [
        _record("record_a", sequence="A" * 1000),
        _record("record_b", sequence="A" * 1000),
    ]
    canvas_config = _orthogroup_alignment_canvas_config()

    shift_x, width_extension = calculate_orthogroup_alignment_canvas_adjustment(
        records,
        canvas_config,
        {1: -300.0},
    )

    assert shift_x == pytest.approx(300.0)
    assert width_extension == pytest.approx(300.0)

    extents = calculate_orthogroup_alignment_canvas_extents(
        records,
        canvas_config,
        {1: -300.0},
    )
    assert extents.ruler_offset_x == pytest.approx(-300.0)
    assert extents.ruler_width == pytest.approx(1300.0)


@pytest.mark.linear
def test_orthogroup_alignment_canvas_adjustment_extends_positive_record_offsets() -> None:
    records = [
        _record("record_a", sequence="A" * 1000),
        _record("record_b", sequence="A" * 1000),
    ]
    canvas_config = _orthogroup_alignment_canvas_config()

    shift_x, width_extension = calculate_orthogroup_alignment_canvas_adjustment(
        records,
        canvas_config,
        {1: 250.0},
    )

    assert shift_x == pytest.approx(0.0)
    assert width_extension == pytest.approx(250.0)

    extents = calculate_orthogroup_alignment_canvas_extents(
        records,
        canvas_config,
        {1: 250.0},
    )
    assert extents.ruler_offset_x == pytest.approx(0.0)
    assert extents.ruler_width == pytest.approx(1250.0)


@pytest.mark.linear
def test_orthogroup_alignment_canvas_extents_use_shifted_record_bounds_for_ruler() -> None:
    records = [
        _record("record_a", sequence="A" * 1000),
        _record("record_b", sequence="A" * 1000),
    ]
    canvas_config = _orthogroup_alignment_canvas_config()

    shift_x, width_extension = calculate_orthogroup_alignment_canvas_adjustment(
        records,
        canvas_config,
        {0: 250.0, 1: 250.0},
    )
    extents = calculate_orthogroup_alignment_canvas_extents(
        records,
        canvas_config,
        {0: 250.0, 1: 250.0},
    )

    assert shift_x == pytest.approx(0.0)
    assert width_extension == pytest.approx(250.0)
    assert extents.ruler_offset_x == pytest.approx(250.0)
    assert extents.ruler_width == pytest.approx(1000.0)

    shift_x, width_extension = calculate_orthogroup_alignment_canvas_adjustment(
        records,
        canvas_config,
        {0: -300.0, 1: -300.0},
    )
    extents = calculate_orthogroup_alignment_canvas_extents(
        records,
        canvas_config,
        {0: -300.0, 1: -300.0},
    )

    assert shift_x == pytest.approx(300.0)
    assert width_extension == pytest.approx(0.0)
    assert extents.ruler_offset_x == pytest.approx(-300.0)
    assert extents.ruler_width == pytest.approx(1000.0)


@pytest.mark.linear
def test_pairwise_match_group_applies_record_specific_alignment_offsets() -> None:
    records = [
        _record("record_a", sequence="A" * 1000),
        _record("record_b", sequence="A" * 1000),
    ]
    comparison = pd.DataFrame.from_records(
        [
            {
                **_hit_row("record_a", "record_b", identity=90),
                "qstart": 100,
                "qend": 200,
                "sstart": 400,
                "send": 500,
            }
        ]
    )
    canvas_config = _orthogroup_alignment_canvas_config()
    blast_config = type(
        "BlastConfig",
        (),
        {
            "fill_color": "#cccccc",
            "identity": 0,
            "min_color": "#eeeeee",
            "max_color": "#111111",
            "fill_opacity": 0.5,
            "stroke_color": "#000000",
            "stroke_width": 0.1,
        },
    )()

    group = PairWiseMatchGroup(
        canvas_config,
        {"record_a": 1000, "record_b": 1000},
        comparison,
        100.0,
        1,
        blast_config,
        records,
        record_offsets_x={1: -300.0},
    ).get_group()

    path = group.elements[0]
    assert 'd="M 100.0,0L200.0,0 L200.0,100.0L100.0,100.0 z"' in path.tostring()


@pytest.mark.linear
def test_assemble_linear_diagram_accepts_precomputed_protein_comparisons(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, object] = {}

    def fake_assemble(*_args, **kwargs):
        captured["comparison_dataframes"] = kwargs.get("comparison_dataframes")
        return Drawing(filename="dummy.svg")

    monkeypatch.setattr(api_diagram_module, "assemble_linear_diagram", fake_assemble)

    records = [_record("record_a"), _record("record_b")]
    comparison = pd.DataFrame.from_records(
        [_hit_row("record_a", "record_b")],
        columns=COMPARISON_COLUMNS,
    )
    canvas = assemble_linear_diagram_from_records(
        records,
        cfg=apply_config_overrides(None, None),
        protein_comparisons=[comparison],
    )

    assert isinstance(canvas, Drawing)
    frames = captured["comparison_dataframes"]
    assert isinstance(frames, list)
    assert len(frames) == 1
    assert frames[0].equals(comparison)


@pytest.mark.linear
def test_build_linear_diagram_forwards_protein_blastp_options(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, object] = {}

    def fake_assemble(*_args, **kwargs):
        captured.update(kwargs)
        return Drawing(filename="dummy.svg")

    monkeypatch.setattr(api_diagram_module, "assemble_linear_diagram_from_records", fake_assemble)

    canvas = api_diagram_module.build_linear_diagram(
        [_record("record_a"), _record("record_b")],
        options=LinearDiagramOptions(
            protein_blastp_mode="orthogroup",
            losatp_bin="custom-losat",
            losatp_threads=8,
            protein_blastp_max_hits=7,
            protein_blastp_candidate_limit=99,
        ),
    )

    assert isinstance(canvas, Drawing)
    assert captured["protein_blastp_mode"] == "orthogroup"
    assert captured["losatp_bin"] == "custom-losat"
    assert captured["losatp_threads"] == 8
    assert captured["protein_blastp_max_hits"] == 7
    assert captured["protein_blastp_candidate_limit"] == 99
    assert captured["orthogroup_membership_mode"] == "anchor_core_v1"
    assert captured["align_orthogroup_feature"] is None


@pytest.mark.linear
def test_build_linear_diagram_forwards_ncbi_blastp_bin(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, object] = {}

    def fake_assemble(*_args, **kwargs):
        captured.update(kwargs)
        return Drawing(filename="dummy.svg")

    monkeypatch.setattr(api_diagram_module, "assemble_linear_diagram_from_records", fake_assemble)

    canvas = api_diagram_module.build_linear_diagram(
        [_record("record_a"), _record("record_b")],
        options=LinearDiagramOptions(
            protein_blastp_mode="pairwise",
            ncbi_blastp_bin="/opt/ncbi/bin/blastp",
        ),
    )

    assert isinstance(canvas, Drawing)
    assert captured["protein_blastp_mode"] == "pairwise"
    assert captured["ncbi_blastp_bin"] == "/opt/ncbi/bin/blastp"


@pytest.mark.linear
def test_build_linear_diagram_forwards_orthogroup_alignment_option(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, object] = {}

    def fake_assemble(*_args, **kwargs):
        captured.update(kwargs)
        return Drawing(filename="dummy.svg")

    monkeypatch.setattr(api_diagram_module, "assemble_linear_diagram_from_records", fake_assemble)

    canvas = api_diagram_module.build_linear_diagram(
        [_record("record_a"), _record("record_b")],
        options=LinearDiagramOptions(
            protein_blastp_mode="orthogroup",
            align_orthogroup_feature="fanchor",
        ),
    )

    assert isinstance(canvas, Drawing)
    assert captured["align_orthogroup_feature"] == "fanchor"


@pytest.mark.linear
def test_linear_cli_rejects_blast_with_protein_blastp_mode() -> None:
    with pytest.raises(SystemExit):
        linear_cli_module.linear_main(
            [
                "--gbk",
                "a.gb",
                "b.gb",
                "-b",
                "a_b.tsv",
                "--protein_blastp_mode",
                "pairwise",
            ]
        )


@pytest.mark.linear
def test_linear_cli_requires_two_records_for_protein_blastp_mode(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(request_render_module, "load_gbks", lambda *_args, **_kwargs: [_record("only")])
    monkeypatch.setattr(request_render_module, "read_color_table", lambda _path: None)
    monkeypatch.setattr(request_render_module, "read_feature_visibility_file", lambda _path: None)

    with pytest.raises(ValidationError, match="requires at least two"):
        linear_cli_module.linear_main(
            [
                "--gbk",
                "dummy.gb",
                "--protein_blastp_mode",
                "pairwise",
                "--format",
                "svg",
            ]
        )


@pytest.mark.linear
def test_linear_cli_forwards_protein_blastp_options(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path,
) -> None:
    records = [_record("record_a"), _record("record_b")]
    captured: dict[str, object] = {}

    monkeypatch.setattr(request_render_module, "load_gbks", lambda *_args, **_kwargs: records)
    monkeypatch.setattr(request_render_module, "read_color_table", lambda _path: None)
    monkeypatch.setattr(request_render_module, "read_feature_visibility_file", lambda _path: None)

    def fake_render(canonical_request, **_kwargs):
        resolved = request_render_module.resolve_request(canonical_request)
        captured["canonical_request"] = resolved
        return SimpleNamespace(
            drawing=Drawing(filename=str(tmp_path / "dummy.svg")),
            interactive_context=None,
            records=tuple(item.source.record for item in resolved.records),
            losat_cache_entries=(),
            losat_derived_cache_entries=(),
            protein_identity_manifest=None,
            request=resolved,
        )

    monkeypatch.setattr(linear_cli_module, "render_request", fake_render)

    linear_cli_module.linear_main(
        [
            "--gbk",
            "a.gb",
            "b.gb",
            "--protein_blastp_mode",
            "orthogroup",
            "--losatp_bin",
            "custom-losat",
            "--losatp_threads",
            "6",
            "--protein_blastp_max_hits",
            "9",
            "--protein_blastp_candidate_limit",
            "123",
            "--format",
            "svg",
            "-o",
            str(tmp_path / "out"),
        ]
    )

    canonical_request = captured["canonical_request"]
    assert isinstance(canonical_request, LinearDiagramRequest)
    options = canonical_request.options
    assert options.protein_blastp_mode == "orthogroup"
    assert options.losatp_bin == "custom-losat"
    assert options.losatp_threads == 6
    assert options.protein_blastp_max_hits == 9
    assert options.protein_blastp_candidate_limit == 123
    assert options.orthogroup_membership_mode == "anchor_core_v1"
    assert options.align_orthogroup_feature is None


@pytest.mark.linear
def test_linear_cli_forwards_ncbi_blastp_bin(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path,
) -> None:
    records = [_record("record_a"), _record("record_b")]
    captured: dict[str, object] = {}

    monkeypatch.setattr(request_render_module, "load_gbks", lambda *_args, **_kwargs: records)
    monkeypatch.setattr(request_render_module, "read_color_table", lambda _path: None)
    monkeypatch.setattr(request_render_module, "read_feature_visibility_file", lambda _path: None)

    def fake_render(canonical_request, **_kwargs):
        resolved = request_render_module.resolve_request(canonical_request)
        captured["canonical_request"] = resolved
        return SimpleNamespace(
            drawing=Drawing(filename=str(tmp_path / "dummy.svg")),
            interactive_context=None,
            records=tuple(item.source.record for item in resolved.records),
            losat_cache_entries=(),
            losat_derived_cache_entries=(),
            protein_identity_manifest=None,
            request=resolved,
        )

    monkeypatch.setattr(linear_cli_module, "render_request", fake_render)

    linear_cli_module.linear_main(
        [
            "--gbk",
            "a.gb",
            "b.gb",
            "--protein_blastp_mode",
            "pairwise",
            "--ncbi_blastp_bin",
            "/opt/ncbi/bin/blastp",
            "--format",
            "svg",
            "-o",
            str(tmp_path / "out"),
        ]
    )

    canonical_request = captured["canonical_request"]
    assert isinstance(canonical_request, LinearDiagramRequest)
    options = canonical_request.options
    assert options.protein_blastp_mode == "pairwise"
    assert options.ncbi_blastp_bin == "/opt/ncbi/bin/blastp"


@pytest.mark.linear
def test_linear_cli_forwards_orthogroup_alignment_option(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path,
) -> None:
    records = [_record("record_a"), _record("record_b")]
    captured: dict[str, object] = {}

    monkeypatch.setattr(request_render_module, "load_gbks", lambda *_args, **_kwargs: records)
    monkeypatch.setattr(request_render_module, "read_color_table", lambda _path: None)
    monkeypatch.setattr(request_render_module, "read_feature_visibility_file", lambda _path: None)

    def fake_render(canonical_request, **_kwargs):
        resolved = request_render_module.resolve_request(canonical_request)
        captured["canonical_request"] = resolved
        return SimpleNamespace(
            drawing=Drawing(filename=str(tmp_path / "dummy.svg")),
            interactive_context=None,
            records=tuple(item.source.record for item in resolved.records),
            losat_cache_entries=(),
            losat_derived_cache_entries=(),
            protein_identity_manifest=None,
            request=resolved,
        )

    monkeypatch.setattr(linear_cli_module, "render_request", fake_render)

    linear_cli_module.linear_main(
        [
            "--gbk",
            "a.gb",
            "b.gb",
            "--protein_blastp_mode",
            "orthogroup",
            "--align_orthogroup_feature",
            "fanchor",
            "--format",
            "svg",
            "-o",
            str(tmp_path / "out"),
        ]
    )

    canonical_request = captured["canonical_request"]
    assert isinstance(canonical_request, LinearDiagramRequest)
    assert canonical_request.options.align_orthogroup_feature == "fanchor"


@pytest.mark.linear
@pytest.mark.parametrize("reverse", [False, True])
def test_orthogroup_display_pairs_do_not_restrict_all_record_membership(reverse: bool) -> None:
    records = [_long_record(f"record_{index}", f"p{index}") for index in range(5)]
    extraction = extract_cds_proteins(records)
    hits = {
        (query, subject): pd.DataFrame.from_records(
            [_hit_row(f"p{query}", f"p{subject}", alignment_length=1000, qend=1000, send=1000, bitscore=1000)],
            columns=COMPARISON_COLUMNS,
        )
        for query in range(5) for subject in range(5)
    }
    expected_pairs = ((0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4))
    if reverse:
        expected_pairs = tuple((subject, query) for query, subject in expected_pairs)
    result = select_rbh_orthogroup_edges_from_directional_hits(
        hits, extraction.protein_map, comparison_pairs=expected_pairs,
    )
    assert set(result.adjacent_display_edges_by_pair) == set(expected_pairs)
    for (query, subject), frame in result.adjacent_display_edges_by_pair.items():
        assert list(zip(frame["query"], frame["subject"])) == [(f"p{query}", f"p{subject}")]
    without_display = select_rbh_orthogroup_edges_from_directional_hits(
        hits, extraction.protein_map, comparison_pairs=(),
    )
    assert without_display.adjacent_display_edges_by_pair == {}
    assert without_display.orthogroups == result.orthogroups
