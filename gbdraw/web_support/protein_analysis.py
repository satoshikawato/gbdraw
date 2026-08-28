"""Thin Worker adapter for canonical browser protein analysis."""

from __future__ import annotations

import copy
from typing import Any, Mapping, Sequence, cast

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]
from pandas import DataFrame  # type: ignore[reportMissingImports]

from gbdraw.analysis.collinearity import CollinearityResult
from gbdraw.analysis.protein_artifacts import validate_current_derived_protein_artifacts
from gbdraw.analysis.protein_colinearity import (
    LosatpCacheManager,
    ProteinExtractionResult,
    build_protein_losat_invocation,
    proteins_to_fasta,
)
from gbdraw.api.diagram import LinearDiagramMetadata, _run_generated_protein_analysis
from gbdraw.api.request_render import _PreparedLinearArtifacts, _build_current_derived_entries
from gbdraw.api.requests import (
    InMemoryRecordSource,
    LinearDiagramRequest,
    RecordInput,
    RecordPresentation,
)
from gbdraw.exceptions import ValidationError
from gbdraw.io.comparisons import COMPARISON_COLUMNS
from gbdraw.session_request_codec import _decode_typed_tree
from gbdraw.session_request_codec import decode_canonical_protein_comparison_intent


def _object(value: object, name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise ValidationError(f"{name} must be an object.")
    return value


def _decode_intent(value: object) -> tuple[Any, dict[str, Any]]:
    intent = _object(value, "proteinIntent")
    if set(intent) != {"generatedProteinComparison", "diagramOptions"}:
        raise ValidationError("proteinIntent has an invalid shape.")
    generated = _object(
        intent["generatedProteinComparison"],
        "proteinIntent.generatedProteinComparison",
    )
    diagram_options = _object(intent["diagramOptions"], "proteinIntent.diagramOptions")
    options = decode_canonical_protein_comparison_intent(
        generated,
        diagram_options=diagram_options,
    )
    settings = _object(generated.get("settings"), "generated protein settings")
    requested = {
        "mode": copy.deepcopy(generated.get("mode")),
        **copy.deepcopy(dict(diagram_options)),
        **copy.deepcopy(dict(settings)),
    }
    if generated.get("pairs"):
        requested["pairs"] = copy.deepcopy(generated["pairs"])
    return options, requested


def _validated_inputs(
    records: Sequence[SeqRecord],
    extraction: ProteinExtractionResult,
    record_payloads: Sequence[Mapping[str, Any]],
    tool_identity: object,
) -> tuple[tuple[SeqRecord, ...], str]:
    resolved_records = tuple(records)
    if not resolved_records or len(resolved_records) != len(record_payloads):
        raise ValidationError("Protein records and transport payloads must align.")
    if len(extraction.proteins_by_record) != len(resolved_records):
        raise ValidationError("Protein extraction and records must align.")
    if extraction.identity_manifest is None:
        raise ValidationError("Protein identity manifest is required.")
    identity = str(tool_identity).strip()
    if not identity:
        raise ValidationError("toolIdentity must be a stable non-empty identity.")
    return resolved_records, identity


def _empty_runner(_query_fasta: str, _subject_fasta: str) -> DataFrame:
    return DataFrame(columns=COMPARISON_COLUMNS)


def _cached_collinearity_result(entry: Mapping[str, Any]) -> CollinearityResult:
    payload = _object(entry["payload"], "derived payload")
    resource = _object(
        payload.get("collinearityResult"),
        "derived payload collinearity result",
    )
    return cast(
        CollinearityResult,
        _decode_typed_tree(
            resource["value"],
            CollinearityResult,
            path="derived payload collinearity result.value",
            resource_schema=2,
        ),
    )


def plan_web_protein_analysis(
    intent_value: object,
    records: Sequence[SeqRecord],
    extraction: ProteinExtractionResult,
    record_payloads: Sequence[Mapping[str, Any]],
    tool_identity: object,
) -> dict[str, Any]:
    """Record raw jobs chosen by the canonical generated-analysis helper."""

    options, requested = _decode_intent(intent_value)
    records, identity = _validated_inputs(
        records,
        extraction,
        record_payloads,
        tool_identity,
    )
    fasta_indexes = {
        proteins_to_fasta(proteins): index
        for index, proteins in enumerate(extraction.proteins_by_record)
    }
    if len(fasta_indexes) != len(records):
        raise ValidationError("Protein records must have distinct runtime bindings.")
    calls: list[tuple[int, int]] = []

    def record_call(query_fasta: str, subject_fasta: str) -> DataFrame:
        try:
            calls.append((fasta_indexes[query_fasta], fasta_indexes[subject_fasta]))
        except KeyError as exc:
            raise ValidationError("Protein helper used an unknown runtime binding.") from exc
        return _empty_runner(query_fasta, subject_fasta)

    _run_generated_protein_analysis(
        records,
        options,
        protein_extraction=extraction,
        runner=record_call,
    )
    display_pairs = (
        tuple(options.protein_comparison_pairs)
        if options.protein_blastp_mode == "pairwise"
        and options.protein_comparison_pairs is not None
        else tuple((index, index + 1) for index in range(len(records) - 1))
    )
    display_ordinals = {pair: index for index, pair in enumerate(display_pairs)}
    jobs = []
    seen_keys: set[str] = set()
    for query_index, subject_index in calls:
        invocation = build_protein_losat_invocation(
            extraction.identity_manifest,
            query_record_instance_key=extraction.record_instance_keys[query_index],
            subject_record_instance_key=extraction.record_instance_keys[subject_index],
            candidate_limit=options.protein_blastp_candidate_limit,
            max_hsps_per_subject=(
                1 if options.protein_blastp_mode == "pairwise" else None
            ),
            tool_identity=identity,
        )
        key = str(invocation["key"])
        if key in seen_keys:
            raise ValidationError("Protein helper planned a duplicate raw identity.")
        seen_keys.add(key)
        pair = (query_index, subject_index)
        display = pair in display_ordinals
        jobs.append(
            {
                "queryIndex": query_index,
                "subjectIndex": subject_index,
                "pairIndex": display_ordinals.get(pair, len(jobs)),
                "display": display,
                "querySequenceKey": str(record_payloads[query_index]["sequenceKey"]),
                "subjectSequenceKey": str(record_payloads[subject_index]["sequenceKey"]),
                "rawEntry": {
                    "schema": 4,
                    "kind": "raw-losat",
                    "identityKind": "protein",
                    "idEncoding": "runtime-handle-v1",
                    **invocation,
                    "filename": (
                        f"{records[query_index].id}.{records[subject_index].id}.losatp.tsv"
                        if display
                        else ""
                    ),
                    "display": display,
                },
            }
        )
    return {"mode": options.protein_blastp_mode, "requested": requested, "jobs": jobs}


def _artifact_request(
    records: Sequence[SeqRecord],
    extraction: ProteinExtractionResult,
    options: Any,
    record_payloads: Sequence[Mapping[str, Any]],
) -> LinearDiagramRequest:
    return LinearDiagramRequest(
        records=tuple(
            RecordInput(
                source=InMemoryRecordSource(record),
                record_key=extraction.record_instance_keys[index],
                presentation=RecordPresentation(
                    reverse_complement=bool(
                        _object(record_payloads[index]["viewTransform"], "viewTransform")[
                            "reverse"
                        ]
                    )
                ),
            )
            for index, record in enumerate(records)
        ),
        options=options,
    )


def _derived_entry(
    metadata: LinearDiagramMetadata,
    request: LinearDiagramRequest,
    records: tuple[SeqRecord, ...],
    extraction: ProteinExtractionResult,
    raw_entries: Sequence[Mapping[str, Any]],
    requested: Mapping[str, Any],
    cached: Mapping[str, Any] | None = None,
) -> Mapping[str, Any]:
    artifacts = _PreparedLinearArtifacts(
        cache=None,
        extraction=extraction,
        nucleotide_entries=(),
        passthrough_derived_entries=((cached,) if cached is not None else ()),
        source_mode=request.options.protein_blastp_mode,
    )
    entries = _build_current_derived_entries(
        metadata,
        request,
        records,
        artifacts,
        raw_entries,
        requested_values=requested,
    )
    if len(entries) != 1:
        raise ValidationError("Protein analysis must produce one derived artifact.")
    return entries[0]


def assemble_web_protein_analysis(
    intent_value: object,
    records: Sequence[SeqRecord],
    extraction: ProteinExtractionResult,
    record_payloads: Sequence[Mapping[str, Any]],
    raw_entries: Sequence[Mapping[str, Any]],
    tool_identity: object,
    cached_derived_entries: Sequence[Mapping[str, Any]] = (),
) -> dict[str, Any]:
    """Admit one derived hit or assemble one result from validated raw entries."""

    options, requested = _decode_intent(intent_value)
    records, identity = _validated_inputs(
        records,
        extraction,
        record_payloads,
        tool_identity,
    )
    planned = plan_web_protein_analysis(
        intent_value,
        records,
        extraction,
        record_payloads,
        identity,
    )
    expected_raw_keys = {str(job["rawEntry"]["key"]) for job in planned["jobs"]}
    admitted_raw = {str(entry.get("key") or ""): entry for entry in raw_entries}
    if (
        len(admitted_raw) != len(raw_entries)
        or set(admitted_raw) != expected_raw_keys
    ):
        raise ValidationError(
            "Protein raw artifacts do not match the canonical planned invocations."
        )
    request = _artifact_request(records, extraction, options, record_payloads)
    expected_metadata = _run_generated_protein_analysis(
        records,
        options,
        protein_extraction=extraction,
        runner=_empty_runner,
    )
    expected = _derived_entry(
        expected_metadata,
        request,
        records,
        extraction,
        raw_entries,
        requested,
    )
    cached_derived = next(
        (
            entry
            for entry in cached_derived_entries
            if entry.get("key") == expected.get("key")
        ),
        None,
    )
    if cached_derived is not None:
        validate_current_derived_protein_artifacts(
            (cached_derived,),
            extraction.identity_manifest.to_dict(),
        )
        cached_metadata = expected_metadata
        if options.protein_blastp_mode == "collinear":
            cached_metadata = _run_generated_protein_analysis(
                records,
                options,
                protein_extraction=extraction,
                collinearity_result=_cached_collinearity_result(cached_derived),
            )
        entry = _derived_entry(
            cached_metadata,
            request,
            records,
            extraction,
            raw_entries,
            requested,
            cached=cached_derived,
        )
        cache_hit = True
    else:
        cache = LosatpCacheManager(
            raw_entries,
            identity_manifest=extraction.identity_manifest,
            tool_identity=identity,
            allow_execution=False,
        )
        cache.set_protein_extraction(extraction)
        metadata = _run_generated_protein_analysis(
            records,
            options,
            losatp_cache=cache,
            protein_extraction=extraction,
        )
        entry = _derived_entry(
            metadata,
            request,
            records,
            extraction,
            cache.session_entries(),
            requested,
        )
        cache_hit = False
    payload = _object(entry["payload"], "derived payload")
    return {
        "cacheHit": cache_hit,
        "derivedEntry": entry,
        "pairs": copy.deepcopy(payload.get("pairs") or []),
        "orthogroupResult": copy.deepcopy(payload.get("orthogroupResult")),
        "collinearityResult": copy.deepcopy(payload.get("collinearityResult")),
        "provenance": copy.deepcopy(payload.get("provenance")),
    }


__all__ = ["assemble_web_protein_analysis", "plan_web_protein_analysis"]
