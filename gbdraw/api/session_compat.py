"""Compatibility adapter from persisted sessions to current typed rendering."""

from __future__ import annotations

import copy
import csv
import io
import json
import re
from dataclasses import dataclass, field, fields, is_dataclass, replace
from pathlib import Path
from typing import Any, Mapping, Sequence

from pandas import DataFrame  # type: ignore[reportMissingImports]

from gbdraw.analysis.protein_colinearity import (
    ProteinExtractionResult,
    build_legacy_protein_reference_map,
    is_legacy_protein_losat_cache_entry,
    is_protein_losat_cache_entry,
    make_legacy_protein_raw_candidate,
    promote_legacy_protein_raw_cache_entries,
    proteins_to_fasta,
    validate_legacy_protein_raw_candidate_envelope,
    validate_protein_raw_entry_references,
)
from gbdraw.exceptions import ValidationError
from gbdraw.features.instance_identity import FEATURE_INSTANCE_HASH_QUALIFIER
from gbdraw.features.semantic_selectors import FEATURE_SEMANTIC_SCOPE_QUALIFIER
from gbdraw.session_io import (
    classify_raw_losat_cache_entry,
    empty_protein_identity_manifest,
    validate_session,
)

from .request_render import (
    CircularBatchRenderResult,
    CircularBatchRequestPlan,
    CurrentRequestArtifacts,
    DiagramRequestPlan,
    LinearRequestPlan,
    PreparedCircularBatchRequest,
    PreparedDiagramRequest,
    RequestRenderResult,
    _extract_linear_request_proteins,
    build_request_plan_diagram,
    plan_request,
    render_prepared_request,
)
from .requests import (
    DiagramRequest,
    LinearDiagramRequest,
)

_LEGACY_LINEAR_TRACK_SLOT_SESSION_VERSION = 32
_LEGACY_MULTILINE_CONSERVATION_LABEL_SESSION_VERSION = 39
_STALE_COLOR_OVERLAY_SESSION_VERSION = 40
_SPECIFIC_COLOR_HEADER = (
    "feature_type",
    "qualifier_key",
    "value",
    "color",
    "caption",
)
_DEFAULT_COLOR_HEADER = ("feature_type", "color")
_COLOR_PATTERN = re.compile(
    r"^(?:none|#(?:[0-9a-f]{3}|[0-9a-f]{4}|[0-9a-f]{6}|[0-9a-f]{8})|[a-z]+)$",
    re.IGNORECASE,
)
_COMPARISON_DEFAULT_COLORS = {
    "pairwise_match": "#d3d3d3",
    "pairwise_match_min": "#ffe7e7",
    "pairwise_match_max": "#ff7272",
    "collinear_block_plus_min": "#f0f1f5",
    "collinear_block_plus": "#8b9cc1",
    "collinear_block_minus_min": "#ffe7e7",
    "collinear_block_minus": "#e15759",
}
_LEGACY_PROTEIN_REFERENCE_RE = re.compile(
    r"p_[A-Za-z0-9._%+-]+?_\d+_\d+_(?:-1|0|1)_[0-9a-f]{12}"
    r"(?:_[2-9][0-9]*)?"
)


@dataclass(frozen=True)
class _V40ColorRecoveryPlan:
    specific_resource_id: str | None
    default_resource_id: str | None
    recovered_palette_name: str | None
_FEATURE_ANALYSIS_REFERENCE_RE = re.compile(r"f_[0-9a-f]{64}")


@dataclass(frozen=True)
class SessionMigrationReport:
    """Compatibility evidence retained when a session enters the current core."""

    protein_raw_candidates: tuple[Mapping[str, Any], ...] = ()
    protein_derived_evidence: tuple[Mapping[str, Any], ...] = ()
    protein_id_map: Mapping[str, str] | None = None
    warnings: tuple[str, ...] = ()


@dataclass(frozen=True)
class SessionCompatibleRequestRenderResult(RequestRenderResult):
    """A current render result paired with persisted migration evidence."""

    migration_report: SessionMigrationReport = field(
        default_factory=SessionMigrationReport
    )

    @property
    def legacy_protein_raw_candidates(self) -> tuple[Mapping[str, Any], ...]:
        return self.migration_report.protein_raw_candidates

    @property
    def legacy_protein_derived_evidence(self) -> tuple[Mapping[str, Any], ...]:
        return self.migration_report.protein_derived_evidence

    @property
    def protein_id_map(self) -> Mapping[str, str]:
        return self.migration_report.protein_id_map or {}

    @property
    def warnings(self) -> tuple[str, ...]:
        return self.migration_report.warnings


@dataclass(frozen=True)
class AdaptedSessionRequest:
    """A persisted request converted to current typed artifacts."""

    request: DiagramRequest
    artifacts: CurrentRequestArtifacts
    migration_report: SessionMigrationReport


@dataclass(frozen=True)
class _SessionArtifactSource:
    current_raw_entries: tuple[Mapping[str, Any], ...]
    current_derived_entries: tuple[Mapping[str, Any], ...]
    protein_identity_manifest: Mapping[str, Any] | None
    protein_source_mode: str | None
    legacy_candidates: tuple[Mapping[str, Any], ...]
    legacy_derived_evidence: tuple[Mapping[str, Any], ...]


def _artifact_entries(
    artifacts: Mapping[str, Any],
    field: str,
) -> tuple[Mapping[str, Any], ...]:
    container = artifacts.get(field)
    if container is None:
        return ()
    if not isinstance(container, Mapping):
        raise ValidationError(f"Session {field} must be an object when present.")
    entries = container.get("entries", [])
    if not isinstance(entries, list) or not all(
        isinstance(entry, Mapping) for entry in entries
    ):
        raise ValidationError(f"Session {field}.entries must be an array of objects.")
    return tuple(copy.deepcopy(dict(entry)) for entry in entries)


def _legacy_candidate_entries(
    artifacts: Mapping[str, Any],
    direct_entries: Sequence[Mapping[str, Any]],
) -> tuple[Mapping[str, Any], ...]:
    candidates: list[Mapping[str, Any]] = []
    legacy_artifacts = artifacts.get("legacyArtifacts")
    if legacy_artifacts is not None and not isinstance(legacy_artifacts, Mapping):
        raise ValidationError("Session legacyArtifacts must be an object.")
    envelope = (
        legacy_artifacts.get("proteinRawCandidates")
        if isinstance(legacy_artifacts, Mapping)
        else None
    )
    if envelope is not None:
        if not isinstance(envelope, Mapping):
            raise ValidationError(
                "legacyArtifacts.proteinRawCandidates must be an object."
            )
        normalized = validate_legacy_protein_raw_candidate_envelope(envelope)
        entries = normalized.get("entries")
        assert isinstance(entries, list)
        candidates.extend(
            candidate for candidate in entries if isinstance(candidate, Mapping)
        )
    candidates.extend(
        make_legacy_protein_raw_candidate(entry) for entry in direct_entries
    )

    result: list[Mapping[str, Any]] = []
    seen: set[str] = set()
    for candidate in candidates:
        original = candidate.get("originalEntry")
        if not isinstance(original, Mapping):
            raise ValidationError("Legacy protein candidate has no original entry.")
        fingerprint = json.dumps(
            original,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
        )
        if fingerprint in seen:
            continue
        seen.add(fingerprint)
        result.append(copy.deepcopy(dict(candidate)))
    return tuple(result)


def _is_derived_entry(entry: Mapping[str, Any], *, schema: int) -> bool:
    return (
        entry.get("schema") == schema
        and entry.get("kind") == "derived-losatp-payload"
        and isinstance(entry.get("key"), str)
        and bool(entry.get("key"))
        and isinstance(entry.get("payload"), Mapping)
    )


def _legacy_derived_evidence_entries(
    artifacts: Mapping[str, Any],
    direct_entries: Sequence[Mapping[str, Any]],
) -> tuple[Mapping[str, Any], ...]:
    evidence = [copy.deepcopy(dict(entry)) for entry in direct_entries]
    legacy_artifacts = artifacts.get("legacyArtifacts")
    envelope = (
        legacy_artifacts.get("proteinDerivedEvidence")
        if isinstance(legacy_artifacts, Mapping)
        else None
    )
    if envelope is not None:
        if (
            not isinstance(envelope, Mapping)
            or envelope.get("schema") != 1
            or not isinstance(envelope.get("entries"), list)
        ):
            raise ValidationError(
                "legacyArtifacts.proteinDerivedEvidence must use schema 1 "
                "with an entries array."
            )
        for index, entry in enumerate(envelope["entries"]):
            if not isinstance(entry, Mapping) or not _is_derived_entry(
                entry,
                schema=1,
            ):
                raise ValidationError(
                    "Invalid legacy protein derived evidence at "
                    f"entries[{index}]."
                )
            evidence.append(copy.deepcopy(dict(entry)))

    result: list[Mapping[str, Any]] = []
    seen: set[str] = set()
    for entry in evidence:
        fingerprint = json.dumps(
            entry,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
        )
        if fingerprint in seen:
            continue
        seen.add(fingerprint)
        result.append(entry)
    return tuple(result)


def _session_protein_mode(artifacts: Mapping[str, Any]) -> str | None:
    config = artifacts.get("config")
    losat = config.get("losat") if isinstance(config, Mapping) else None
    blastp = losat.get("blastp") if isinstance(losat, Mapping) else None
    configured = blastp.get("mode") if isinstance(blastp, Mapping) else None
    if str(configured or "") in {"pairwise", "orthogroup", "collinear"}:
        return str(configured)
    return None


def _read_session_artifact_source(
    artifacts: Mapping[str, Any],
) -> _SessionArtifactSource:
    raw_entries = _artifact_entries(artifacts, "losatCache")
    current_raw: list[Mapping[str, Any]] = []
    direct_legacy: list[Mapping[str, Any]] = []
    for index, entry in enumerate(raw_entries):
        classification = classify_raw_losat_cache_entry(entry)
        if classification in {"protein-current", "nucleotide-current"}:
            current_raw.append(entry)
        elif classification == "protein-legacy":
            direct_legacy.append(entry)
        else:
            raise ValidationError(
                f"Unsupported LOSAT session artifact at losatCache.entries[{index}]."
            )

    current_derived: list[Mapping[str, Any]] = []
    direct_evidence: list[Mapping[str, Any]] = []
    for index, entry in enumerate(
        _artifact_entries(artifacts, "losatDerivedCache")
    ):
        if _is_derived_entry(entry, schema=3) and (
            entry.get("idEncoding") == "runtime-handle-v1"
        ):
            current_derived.append(entry)
        elif _is_derived_entry(entry, schema=1):
            direct_evidence.append(entry)
        else:
            raise ValidationError(
                "Unsupported derived LOSATP session artifact at "
                f"losatDerivedCache.entries[{index}]."
            )

    source_manifest = artifacts.get("proteinIdentityManifest")
    if source_manifest is not None and not isinstance(source_manifest, Mapping):
        raise ValidationError("Session proteinIdentityManifest must be an object.")
    manifest = (
        copy.deepcopy(dict(source_manifest))
        if isinstance(source_manifest, Mapping)
        else None
    )
    return _SessionArtifactSource(
        current_raw_entries=tuple(current_raw),
        current_derived_entries=tuple(current_derived),
        protein_identity_manifest=manifest,
        protein_source_mode=_session_protein_mode(artifacts),
        legacy_candidates=_legacy_candidate_entries(artifacts, direct_legacy),
        legacy_derived_evidence=_legacy_derived_evidence_entries(
            artifacts,
            direct_evidence,
        ),
    )


def _merge_promoted_tsv_id_map(
    target: dict[str, str],
    *,
    legacy_text: str,
    current_text: str,
) -> None:
    def rows(text: str) -> list[list[str]]:
        return [
            line.split("\t")
            for line in str(text).splitlines()
            if line.strip() and not line.lstrip().startswith("#")
        ]

    legacy_rows = rows(legacy_text)
    current_rows = rows(current_text)
    if len(legacy_rows) != len(current_rows):
        raise ValidationError("Promoted LOSATP TSV changed its row count.")
    for old_row, new_row in zip(legacy_rows, current_rows, strict=True):
        if len(old_row) < 2 or len(new_row) < 2:
            raise ValidationError("Promoted LOSATP TSV has an incomplete row.")
        for old_id, new_id in zip(old_row[:2], new_row[:2], strict=True):
            previous = target.get(old_id)
            if previous is not None and previous != new_id:
                raise ValidationError(
                    f"Legacy protein ID {old_id!r} resolves to multiple current IDs."
                )
            target[old_id] = new_id


def _promote_legacy_candidates(
    candidates: Sequence[Mapping[str, Any]],
    extraction: ProteinExtractionResult,
) -> tuple[
    tuple[Mapping[str, Any], ...],
    tuple[Mapping[str, Any], ...],
    Mapping[str, str],
]:
    manifest = extraction.identity_manifest
    if manifest is None:
        raise ValidationError("Protein cache promotion requires an identity manifest.")
    fastas = tuple(proteins_to_fasta(items) for items in extraction.proteins_by_record)
    promoted_entries: list[Mapping[str, Any]] = []
    unresolved_candidates: list[Mapping[str, Any]] = []
    id_map: dict[str, str] = {}

    for candidate in candidates:
        original = candidate.get("originalEntry")
        if not isinstance(original, Mapping) or not is_legacy_protein_losat_cache_entry(
            original
        ):
            raise ValidationError(
                "Legacy protein candidate is not a schema-2 blastp entry."
            )
        raw_args = original.get("args")
        expected_args = (
            tuple(str(arg) for arg in raw_args)
            if isinstance(raw_args, list)
            else ()
        )
        promotion = None
        rejection_reasons: list[str] = []
        for query_index, query_proteins in enumerate(extraction.proteins_by_record):
            for subject_index, subject_proteins in enumerate(
                extraction.proteins_by_record
            ):
                try:
                    scan = promote_legacy_protein_raw_cache_entries(
                        (candidate,),
                        query_proteins=query_proteins,
                        subject_proteins=subject_proteins,
                        query_fasta=fastas[query_index],
                        subject_fasta=fastas[subject_index],
                        identity_manifest=manifest,
                        expected_args=expected_args,
                        expected_program=str(original.get("program") or "blastp"),
                        expected_outfmt=str(original.get("outfmt") or "6"),
                    )
                except ValidationError as exc:
                    rejection_reasons.append(str(exc))
                    continue
                if scan.promotion is not None:
                    promotion = scan.promotion
                    break
                rejection_reasons.extend(item.reason for item in scan.rejections)
            if promotion is not None:
                break

        if promotion is None:
            reason = next(
                (
                    item
                    for item in reversed(rejection_reasons)
                    if item and "hash does not match" not in item
                ),
                "No current record pair matched the legacy cache identity.",
            )
            unresolved_candidates.append(
                make_legacy_protein_raw_candidate(
                    original,
                    state="rejected",
                    rejection_reason=reason,
                )
            )
            continue
        if not validate_protein_raw_entry_references(promotion.entry, manifest):
            raise ValidationError(
                "Promoted protein raw entry failed manifest validation."
            )
        promoted_entries.append(promotion.entry)
        _merge_promoted_tsv_id_map(
            id_map,
            legacy_text=str(original.get("text") or ""),
            current_text=promotion.rewritten_tsv,
        )

    return tuple(promoted_entries), tuple(unresolved_candidates), id_map


def _rewrite_protein_reference_string(
    value: str,
    id_map: Mapping[str, str],
) -> str:
    if value in id_map:
        return str(id_map[value])
    if "p_r_" not in value:
        return value
    return _LEGACY_PROTEIN_REFERENCE_RE.sub(
        lambda match: id_map.get(match.group(0), match.group(0)),
        value,
    )


def rewrite_protein_artifact_references(
    value: Any,
    id_map: Mapping[str, str],
) -> Any:
    """Return a detached compatibility value with verified IDs rewritten."""

    if isinstance(value, str):
        return _rewrite_protein_reference_string(value, id_map)
    if isinstance(value, DataFrame):
        rewritten = value.copy(deep=True)
        for column in rewritten.columns:
            rewritten[column] = rewritten[column].map(
                lambda item: (
                    _rewrite_protein_reference_string(item, id_map)
                    if isinstance(item, str)
                    else item
                )
            )
        return rewritten
    if is_dataclass(value) and not isinstance(value, type):
        return replace(
            value,
            **{
                field.name: rewrite_protein_artifact_references(
                    getattr(value, field.name),
                    id_map,
                )
                for field in fields(value)
            },
        )
    if isinstance(value, Mapping):
        rewritten_mapping: dict[Any, Any] = {}
        for key, item in value.items():
            rewritten_key = rewrite_protein_artifact_references(key, id_map)
            if rewritten_key in rewritten_mapping:
                raise ValidationError(
                    "Protein reference migration produced duplicate key "
                    f"{rewritten_key!r}."
                )
            rewritten_mapping[rewritten_key] = rewrite_protein_artifact_references(
                item,
                id_map,
            )
        return rewritten_mapping
    if isinstance(value, tuple):
        return tuple(
            rewrite_protein_artifact_references(item, id_map) for item in value
        )
    if isinstance(value, list):
        return [
            rewrite_protein_artifact_references(item, id_map) for item in value
        ]
    if isinstance(value, frozenset):
        return frozenset(
            rewrite_protein_artifact_references(item, id_map) for item in value
        )
    if isinstance(value, set):
        return {
            rewrite_protein_artifact_references(item, id_map) for item in value
        }
    return copy.deepcopy(value)


def _contains_legacy_protein_reference(value: Any) -> bool:
    if isinstance(value, str):
        return (
            _LEGACY_PROTEIN_REFERENCE_RE.search(value) is not None
            or _FEATURE_ANALYSIS_REFERENCE_RE.search(value) is not None
        )
    if isinstance(value, DataFrame):
        return any(
            _contains_legacy_protein_reference(item)
            for item in value.to_numpy().ravel()
        )
    if is_dataclass(value) and not isinstance(value, type):
        return any(
            _contains_legacy_protein_reference(getattr(value, field.name))
            for field in fields(value)
        )
    if isinstance(value, Mapping):
        return any(
            _contains_legacy_protein_reference(key)
            or _contains_legacy_protein_reference(item)
            for key, item in value.items()
        )
    if isinstance(value, (tuple, list, set, frozenset)):
        return any(_contains_legacy_protein_reference(item) for item in value)
    return False


def _legacy_protein_reference_ids(value: Any) -> set[str]:
    references: set[str] = set()
    if isinstance(value, str):
        references.update(
            match.group(0)
            for match in _LEGACY_PROTEIN_REFERENCE_RE.finditer(value)
        )
    elif isinstance(value, DataFrame):
        for item in value.to_numpy().ravel():
            references.update(_legacy_protein_reference_ids(item))
    elif is_dataclass(value) and not isinstance(value, type):
        for field in fields(value):
            references.update(
                _legacy_protein_reference_ids(getattr(value, field.name))
            )
    elif isinstance(value, Mapping):
        for key, item in value.items():
            references.update(_legacy_protein_reference_ids(key))
            references.update(_legacy_protein_reference_ids(item))
    elif isinstance(value, (tuple, list, set, frozenset)):
        for item in value:
            references.update(_legacy_protein_reference_ids(item))
    return references


def _request_protein_artifacts(
    request: LinearDiagramRequest,
) -> tuple[Any, ...]:
    options = request.options
    return (
        options.linear_comparisons,
        options.protein_comparisons,
        options.orthogroups,
        options.collinearity_blocks,
        options.align_orthogroup_feature,
    )


def _rewrite_linear_request_protein_references(
    request: LinearDiagramRequest,
    id_map: Mapping[str, str],
) -> LinearDiagramRequest:
    if not id_map:
        rewritten_request = request
    else:
        options = request.options
        rewritten_request = replace(
            request,
            options=replace(
                options,
                linear_comparisons=rewrite_protein_artifact_references(
                    options.linear_comparisons,
                    id_map,
                ),
                protein_comparisons=rewrite_protein_artifact_references(
                    options.protein_comparisons,
                    id_map,
                ),
                orthogroups=rewrite_protein_artifact_references(
                    options.orthogroups,
                    id_map,
                ),
                collinearity_blocks=rewrite_protein_artifact_references(
                    options.collinearity_blocks,
                    id_map,
                ),
                align_orthogroup_feature=rewrite_protein_artifact_references(
                    options.align_orthogroup_feature,
                    id_map,
                ),
            ),
        )
    if any(
        _contains_legacy_protein_reference(value)
        for value in _request_protein_artifacts(rewritten_request)
    ):
        raise ValidationError(
            "Precomputed protein artifacts contain an ID that no verified "
            "session artifact resolved."
        )
    return rewritten_request


def _deduplicate_raw_entries(
    entries: Sequence[Mapping[str, Any]],
) -> tuple[Mapping[str, Any], ...]:
    result: list[Mapping[str, Any]] = []
    by_key: dict[str, Mapping[str, Any]] = {}
    for entry in entries:
        key = str(entry.get("key") or "")
        previous = by_key.get(key)
        if previous is not None:
            if dict(previous) != dict(entry):
                raise ValidationError(
                    f"Session artifact migration produced duplicate key {key!r}."
                )
            continue
        clone = copy.deepcopy(dict(entry))
        by_key[key] = clone
        result.append(clone)
    return tuple(result)


def _adapt_session_plan(
    plan: DiagramRequestPlan,
    session_artifacts: Mapping[str, Any],
) -> AdaptedSessionRequest:
    source = _read_session_artifact_source(session_artifacts)
    request = plan.request
    current_raw = source.current_raw_entries
    manifest = source.protein_identity_manifest
    unresolved = source.legacy_candidates
    id_map: dict[str, str] = {}

    if isinstance(plan, LinearRequestPlan):
        reference_ids: set[str] = set()
        for value in _request_protein_artifacts(plan.request):
            reference_ids.update(_legacy_protein_reference_ids(value))
        current_protein = tuple(
            entry for entry in current_raw if is_protein_losat_cache_entry(entry)
        )
        needs_extraction = bool(
            source.legacy_candidates or reference_ids or current_protein
        )
        promoted: tuple[Mapping[str, Any], ...] = ()
        if needs_extraction:
            if plan.inputs is None:  # pragma: no cover - planner invariant.
                raise ValidationError(
                    "Linear session compatibility plan has no prepared inputs."
                )
            extraction = _extract_linear_request_proteins(
                plan.request,
                plan.records,
                plan.inputs,
            )
            extracted_manifest = extraction.identity_manifest
            assert extracted_manifest is not None
            promoted, unresolved, promoted_id_map = _promote_legacy_candidates(
                source.legacy_candidates,
                extraction,
            )
            id_map.update(promoted_id_map)
            if reference_ids:
                artifact_id_map = build_legacy_protein_reference_map(
                    extraction,
                    sorted(reference_ids),
                )
                for old_id, runtime_handle in artifact_id_map.items():
                    previous = id_map.get(old_id)
                    if previous is not None and previous != runtime_handle:
                        raise ValidationError(
                            f"Protein ID {old_id!r} resolves to multiple runtime handles."
                        )
                    id_map[old_id] = runtime_handle
            reusable_current = tuple(
                entry
                for entry in current_protein
                if validate_protein_raw_entry_references(
                    entry,
                    extracted_manifest,
                )
            )
            nucleotide = tuple(
                entry
                for entry in current_raw
                if not is_protein_losat_cache_entry(entry)
            )
            current_raw = _deduplicate_raw_entries(
                (*reusable_current, *promoted, *nucleotide)
            )
            manifest = extracted_manifest.to_dict()
        request = _rewrite_linear_request_protein_references(
            plan.request,
            id_map,
        )
        if plan.request.options.protein_comparisons is not None and promoted:
            unresolved = ()

    if manifest is None and not any(
        is_protein_losat_cache_entry(entry) for entry in current_raw
    ):
        manifest = empty_protein_identity_manifest()
    artifacts = CurrentRequestArtifacts(
        losat_cache_entries=tuple(current_raw),
        losat_derived_cache_entries=source.current_derived_entries,
        protein_identity_manifest=manifest,
        protein_source_mode=source.protein_source_mode,
    )
    warnings = (
        (
            f"{len(unresolved)} legacy protein raw cache candidate(s) "
            "could not be promoted.",
        )
        if unresolved
        else ()
    )
    return AdaptedSessionRequest(
        request=request,
        artifacts=artifacts,
        migration_report=SessionMigrationReport(
            protein_raw_candidates=tuple(unresolved),
            protein_derived_evidence=source.legacy_derived_evidence,
            protein_id_map=id_map,
            warnings=warnings,
        ),
    )


def adapt_session_request(
    request: DiagramRequest,
    session_artifacts: Mapping[str, Any],
) -> AdaptedSessionRequest:
    """Convert one validated persisted session to the current render contract."""

    if not isinstance(session_artifacts, Mapping):
        raise ValidationError("Session compatibility input must be an object.")
    validate_session(session_artifacts)
    return _adapt_session_plan(plan_request(request), session_artifacts)


def _adjust_migration_report(
    prepared: PreparedDiagramRequest | PreparedCircularBatchRequest,
    report: SessionMigrationReport,
) -> SessionMigrationReport:
    if isinstance(prepared, PreparedCircularBatchRequest):
        return report
    evidence = (
        ()
        if prepared.losat_derived_cache_entries
        else report.protein_derived_evidence
    )
    return replace(
        report,
        protein_derived_evidence=evidence,
    )


def build_session_compatible_request_diagram(
    request: DiagramRequest,
    session_artifacts: Mapping[str, Any],
) -> PreparedDiagramRequest | PreparedCircularBatchRequest:
    """Build a released session after adapting it to current artifacts."""

    if not isinstance(session_artifacts, Mapping):
        raise ValidationError("Session compatibility input must be an object.")
    validate_session(session_artifacts)
    prepared, _report = _build_session_compatible_plan(
        plan_request(request),
        session_artifacts,
    )
    return prepared


def _build_session_compatible_plan(
    plan: DiagramRequestPlan,
    session_artifacts: Mapping[str, Any],
) -> tuple[
    PreparedDiagramRequest | PreparedCircularBatchRequest,
    SessionMigrationReport,
]:
    """Adapt and build one already-resolved session request plan."""

    adapted = _adapt_session_plan(plan, session_artifacts)
    if adapted.request is not plan.request:
        plan = replace(plan, request=adapted.request)
    prepared = build_request_plan_diagram(
        plan,
        artifacts=adapted.artifacts,
    )
    return prepared, _adjust_migration_report(
        prepared,
        adapted.migration_report,
    )


def render_session_compatible_request(
    request: DiagramRequest,
    session_artifacts: Mapping[str, Any],
    *,
    include_feature_catalog: bool = False,
) -> RequestRenderResult | CircularBatchRenderResult:
    """Render a released session through its sole compatibility adapter."""

    if not isinstance(session_artifacts, Mapping):
        raise ValidationError("Session compatibility input must be an object.")
    validate_session(session_artifacts)
    plan = plan_request(request)
    batch_outputs_preflighted = isinstance(plan, CircularBatchRequestPlan)
    plan.preflight_outputs()
    prepared, migration_report = _build_session_compatible_plan(
        plan,
        session_artifacts,
    )
    result = render_prepared_request(
        prepared,
        batch_outputs_preflighted=batch_outputs_preflighted,
        include_feature_catalog=include_feature_catalog,
    )
    if not isinstance(prepared, PreparedDiagramRequest) or not isinstance(
        result,
        RequestRenderResult,
    ):
        return result
    return SessionCompatibleRequestRenderResult(
        **{
            item.name: getattr(result, item.name)
            for item in fields(RequestRenderResult)
        },
        migration_report=migration_report,
    )


def canonical_payload_for_session_decode(
    session_version: int,
    payload: Mapping[str, Any],
) -> dict[str, Any]:
    """Return a detached canonical payload with supported migrations applied."""

    detached = copy.deepcopy(dict(payload))
    if session_version == _LEGACY_MULTILINE_CONSERVATION_LABEL_SESSION_VERSION:
        _migrate_legacy_multiline_conservation_labels(detached)
    if session_version > _LEGACY_LINEAR_TRACK_SLOT_SESSION_VERSION:
        return detached

    diagram_options = detached.get("diagramOptions")
    if not isinstance(diagram_options, Mapping):
        return detached
    detached_diagram_options = dict(diagram_options)
    detached["diagramOptions"] = detached_diagram_options

    tracks = detached_diagram_options.get("tracks")
    if not isinstance(tracks, Mapping):
        return detached
    detached_tracks = dict(tracks)
    detached_diagram_options["tracks"] = detached_tracks

    slots = detached_tracks.get("linearTrackSlots")
    if not isinstance(slots, list):
        return detached
    detached_tracks["linearTrackSlots"] = [
        _migrate_legacy_linear_track_slot(slot) for slot in slots
    ]
    return detached


def canonical_projection_for_session_decode(
    session_version: int,
    session: Mapping[str, Any],
    *,
    resource_paths: Mapping[str, str | Path],
    temp_directory: str | Path,
) -> tuple[dict[str, Any], dict[str, str | Path]]:
    """Project one full session envelope to a detached canonical decode basis.

    Version 40 Web sessions could retain an older canonical color resource after
    accepting editor-owned rules or palette colors.  Recovery is deliberately
    narrow: it requires independent editor sources to agree, prepares both
    overlays before writing either one, and never changes the caller's document
    or materialized-resource mapping.
    """

    raw_payload = session.get("renderRequest")
    if not isinstance(raw_payload, Mapping):
        raise ValidationError("Session renderRequest must be an object.")
    payload = canonical_payload_for_session_decode(session_version, raw_payload)
    projected_paths: dict[str, str | Path] = dict(resource_paths)
    if (
        session_version != _STALE_COLOR_OVERLAY_SESSION_VERSION
        or payload.get("schema") != 5
    ):
        return payload, projected_paths

    diagram_options = payload.get("diagramOptions")
    colors_value = (
        diagram_options.get("colors")
        if isinstance(diagram_options, Mapping)
        else None
    )
    if colors_value is None:
        from gbdraw.session_request_codec import CANONICAL_REQUEST_SCHEMA

        payload["schema"] = CANONICAL_REQUEST_SCHEMA
        return payload, projected_paths

    plan = _plan_v40_color_resource_recovery(
        session,
        payload,
        resource_paths=resource_paths,
    )
    colors = _detached_color_options(payload)
    if plan.specific_resource_id is not None:
        colors["colorTable"] = None
        colors["colorTableFile"] = _recovered_resource_reference(
            session,
            plan.specific_resource_id,
        )
    if plan.default_resource_id is not None:
        colors["defaultColors"] = None
        colors["defaultColorsFile"] = _recovered_resource_reference(
            session,
            plan.default_resource_id,
        )
        colors["defaultColorsPalette"] = plan.recovered_palette_name

    # Keep the temporary-directory argument in this compatibility boundary so
    # callers do not gain a second projection API. Recovery reuses the one
    # deterministic saved resource proved by the full envelope.
    del temp_directory
    from gbdraw.session_request_codec import CANONICAL_REQUEST_SCHEMA

    payload["schema"] = CANONICAL_REQUEST_SCHEMA
    return payload, projected_paths


def _detached_color_options(payload: dict[str, Any]) -> dict[str, Any]:
    diagram_options = payload.get("diagramOptions")
    if not isinstance(diagram_options, Mapping):
        raise ValidationError("Session renderRequest.diagramOptions must be an object.")
    detached_diagram_options = dict(diagram_options)
    payload["diagramOptions"] = detached_diagram_options
    colors = detached_diagram_options.get("colors")
    if not isinstance(colors, Mapping):
        raise ValidationError(
            "Session renderRequest.diagramOptions.colors must be an object."
        )
    detached_colors = dict(colors)
    detached_diagram_options["colors"] = detached_colors
    return detached_colors


def _plan_v40_color_resource_recovery(
    session: Mapping[str, Any],
    payload: dict[str, Any],
    *,
    resource_paths: Mapping[str, str | Path],
) -> _V40ColorRecoveryPlan:
    colors = _detached_color_options(payload)
    canonical_specific_id = _canonical_color_resource_id(
        colors,
        primary_field="colorTableFile",
        fallback_field="colorTable",
        label="specific-color",
    )
    canonical_default_id = _canonical_color_resource_id(
        colors,
        primary_field="defaultColorsFile",
        fallback_field="defaultColors",
        label="default-color",
    )
    canonical_rules = _read_optional_color_rows(
        canonical_specific_id,
        header=_SPECIFIC_COLOR_HEADER,
        resource_paths=resource_paths,
    )
    canonical_defaults = _read_optional_color_rows(
        canonical_default_id,
        header=_DEFAULT_COLOR_HEADER,
        resource_paths=resource_paths,
    )

    config_value = session.get("config")
    config = config_value if isinstance(config_value, Mapping) else {}
    ui_value = session.get("ui")
    ui = ui_value if isinstance(ui_value, Mapping) else {}
    normalized_editor_rules = _normalize_editor_rules(config.get("rules"))
    editor_rules = normalized_editor_rules if normalized_editor_rules else None
    config_defaults = _normalize_editor_defaults(config.get("colors"))
    ui_defaults = _normalize_editor_defaults(ui.get("appliedPaletteColors"))

    effective_rules = editor_rules or canonical_rules or ()
    reserved_qualifier = next(
        (
            rule[1]
            for rule in effective_rules
            if rule[1]
            in {
                FEATURE_INSTANCE_HASH_QUALIFIER,
                FEATURE_SEMANTIC_SCOPE_QUALIFIER,
            }
        ),
        None,
    )
    if reserved_qualifier is not None:
        selector_label = (
            "instance"
            if reserved_qualifier == FEATURE_INSTANCE_HASH_QUALIFIER
            else "semantic"
        )
        raise ValidationError(
            "Session v40 cannot safely promote a schema-5 reserved "
            f"{selector_label} selector; Generate again."
        )
    _validate_v40_caption_colors(
        effective_rules,
        session.get("editorState"),
    )

    specific_mismatch = editor_rules is not None and editor_rules != (canonical_rules or ())
    supported_defaults = (
        config_defaults
        if config_defaults is not None and config_defaults == ui_defaults
        else None
    )
    if (
        config_defaults is not None
        and ui_defaults is not None
        and config_defaults != ui_defaults
        and (
            config_defaults != (canonical_defaults or ())
            or ui_defaults != (canonical_defaults or ())
        )
    ):
        raise ValidationError(
            "Session v40 color recovery found conflicting applied color authorities."
        )
    config_palette_name = _normalized_text(config.get("palette"))
    applied_palette_name = _normalized_text(ui.get("appliedPaletteName"))
    canonical_palette_name = _normalized_text(colors.get("defaultColorsPalette")) or "default"
    supported_palette_name = (
        config_palette_name
        if config_palette_name
        and config_palette_name == applied_palette_name
        else None
    )
    default_mismatch = supported_defaults is not None and (
        supported_defaults != (canonical_defaults or ())
        or (
            supported_palette_name is not None
            and supported_palette_name != canonical_palette_name
        )
    )

    specific_resource_id = None
    if specific_mismatch:
        specific_resource_id = _require_unique_recovery_resource(
            session,
            payload,
            resource_paths=resource_paths,
            kind_pattern=re.compile(r"(?:specific.*color|color.*table)", re.IGNORECASE),
            header=_SPECIFIC_COLOR_HEADER,
            expected=effective_rules,
            label="specific-color",
        )

    default_resource_id = None
    if default_mismatch:
        if supported_palette_name is None:
            raise ValidationError(
                "Session v40 default-color recovery found conflicting palette names."
            )
        assert supported_defaults is not None
        default_resource_id = _require_unique_recovery_resource(
            session,
            payload,
            resource_paths=resource_paths,
            kind_pattern=re.compile(r"(?:default.*color|color.*default)", re.IGNORECASE),
            header=_DEFAULT_COLOR_HEADER,
            expected=supported_defaults,
            label="default-color",
        )

    if specific_mismatch or default_mismatch:
        _validate_v40_saved_color_evidence(
            session,
            expected_rules=effective_rules,
            canonical_defaults=canonical_defaults or (),
            recovered_defaults=supported_defaults or canonical_defaults or (),
        )

    return _V40ColorRecoveryPlan(
        specific_resource_id=specific_resource_id,
        default_resource_id=default_resource_id,
        recovered_palette_name=(supported_palette_name if default_mismatch else None),
    )


def _canonical_color_resource_id(
    colors: Mapping[str, Any],
    *,
    primary_field: str,
    fallback_field: str,
    label: str,
) -> str | None:
    resource_ids: list[str] = []
    for field_name in (primary_field, fallback_field):
        reference = colors.get(field_name)
        if reference is None:
            continue
        if not isinstance(reference, Mapping):
            raise ValidationError(
                f"Session color reference {field_name} is invalid."
            )
        resource_id = _normalized_text(reference.get("resourceId"))
        if not resource_id:
            raise ValidationError(
                f"Session color reference {field_name} has no resourceId."
            )
        resource_ids.append(resource_id)
    if len(set(resource_ids)) > 1:
        raise ValidationError(
            f"Legacy color recovery found ambiguous canonical {label} resources."
        )
    return resource_ids[0] if resource_ids else None


def _read_optional_color_rows(
    resource_id: str | None,
    *,
    header: tuple[str, ...],
    resource_paths: Mapping[str, str | Path],
) -> tuple[tuple[str, ...], ...] | None:
    if resource_id is None:
        return None
    return _read_named_color_rows(
        resource_id,
        header=header,
        resource_paths=resource_paths,
        required=True,
    )


def _read_named_color_rows(
    resource_id: str,
    *,
    header: tuple[str, ...],
    resource_paths: Mapping[str, str | Path],
    required: bool,
) -> tuple[tuple[str, ...], ...] | None:
    raw_path = resource_paths.get(resource_id)
    if raw_path is None:
        if required:
            raise ValidationError(
                f"Session color resource {resource_id!r} is not materialized."
            )
        return None
    try:
        text = Path(raw_path).read_text(encoding="utf-8-sig")
        parsed = [
            tuple(cell.strip() for cell in row)
            for row in csv.reader(io.StringIO(text), delimiter="\t")
            if row and any(cell.strip() for cell in row)
        ]
    except (OSError, UnicodeError, csv.Error) as exc:
        raise ValidationError(
            f"Could not inspect v40 color resource {resource_id!r}."
        ) from exc
    if parsed and parsed[0] == header:
        parsed = parsed[1:]
    normalized: list[tuple[str, ...]] = []
    for index, row in enumerate(parsed, start=1):
        if len(row) != len(header):
            raise ValidationError(
                f"V40 color resource {resource_id!r} row {index} must contain "
                f"{len(header)} columns."
            )
        if header == _SPECIFIC_COLOR_HEADER:
            normalized.append(_normalize_specific_row(row, context=resource_id))
        else:
            normalized.append(_normalize_default_row(row, context=resource_id))
    if header == _DEFAULT_COLOR_HEADER:
        if len({row[0] for row in normalized}) != len(normalized):
            raise ValidationError(
                f"V40 color resource {resource_id!r} contains duplicate feature types."
            )
        normalized_by_type = dict(normalized)
        _normalize_comparison_default_colors(normalized_by_type)
        normalized = list(normalized_by_type.items())
        normalized.sort(key=lambda row: row[0])
    return tuple(normalized)


def _referenced_resource_ids(value: Any, target: set[str] | None = None) -> set[str]:
    result = target if target is not None else set()
    if isinstance(value, list):
        for entry in value:
            _referenced_resource_ids(entry, result)
        return result
    if not isinstance(value, Mapping):
        return result
    resource_id = _normalized_text(value.get("resourceId"))
    if resource_id:
        result.add(resource_id)
    for entry in value.values():
        _referenced_resource_ids(entry, result)
    return result


def _require_unique_recovery_resource(
    session: Mapping[str, Any],
    payload: Mapping[str, Any],
    *,
    resource_paths: Mapping[str, str | Path],
    kind_pattern: re.Pattern[str],
    header: tuple[str, ...],
    expected: tuple[tuple[str, ...], ...],
    label: str,
) -> str:
    resources_value = session.get("resources")
    resources = resources_value if isinstance(resources_value, Mapping) else {}
    referenced = _referenced_resource_ids(payload)
    matches: list[str] = []
    for raw_resource_id, resource in resources.items():
        resource_id = str(raw_resource_id)
        if resource_id in referenced or not isinstance(resource, Mapping):
            continue
        identity = " ".join(
            (
                _normalized_text(resource.get("kind")),
                _normalized_text(resource.get("name")),
            )
        )
        if kind_pattern.search(identity) is None:
            continue
        try:
            rows = _read_named_color_rows(
                resource_id,
                header=header,
                resource_paths=resource_paths,
                required=True,
            )
        except ValidationError:
            continue
        if rows == expected:
            matches.append(resource_id)
    matches.sort()
    if len(matches) != 1:
        evidence = "no deterministic" if not matches else "ambiguous"
        raise ValidationError(
            f"Session v40 {label} mismatch has {evidence} saved resource support."
        )
    return matches[0]


def _recovered_resource_reference(
    session: Mapping[str, Any],
    resource_id: str,
) -> dict[str, str]:
    resources = session.get("resources")
    resource = resources.get(resource_id) if isinstance(resources, Mapping) else None
    if not isinstance(resource, Mapping):
        raise ValidationError(
            f"Session v40 recovery resource {resource_id!r} is missing."
        )
    kind = _normalized_text(resource.get("kind"))
    return {
        "resourceId": resource_id,
        "representation": "canonicalTsv" if "canonical" in kind.lower() else "file",
    }


def _validate_v40_caption_colors(
    rules: tuple[tuple[str, str, str, str, str], ...],
    editor_state_value: Any,
) -> None:
    caption_colors: dict[str, str] = {}
    for _, _, _, color, caption in rules:
        if not caption:
            continue
        previous = caption_colors.setdefault(caption, color)
        if previous != color:
            raise ValidationError(
                f"Legacy color recovery found conflicting colors for caption {caption!r}."
            )

    editor_state = (
        editor_state_value if isinstance(editor_state_value, Mapping) else {}
    )
    legend_value = editor_state.get("legend")
    legend = legend_value if isinstance(legend_value, Mapping) else {}
    entries_value = legend.get("entries")
    entries = entries_value if isinstance(entries_value, list) else []
    saved_colors: dict[str, str] = {}
    for entry in entries:
        if not isinstance(entry, Mapping):
            continue
        caption = _normalized_text(entry.get("caption"))
        if not caption:
            continue
        color = _normalize_color(_normalized_text(entry.get("color")))
        previous = saved_colors.setdefault(caption, color)
        if previous != color or (
            caption in caption_colors and caption_colors[caption] != color
        ):
            raise ValidationError(
                "Legacy color recovery found conflicting saved legend state "
                f"for {caption!r}."
            )


def _validate_v40_saved_color_evidence(
    session: Mapping[str, Any],
    *,
    expected_rules: tuple[tuple[str, str, str, str, str], ...],
    canonical_defaults: tuple[tuple[str, str], ...],
    recovered_defaults: tuple[tuple[str, str], ...],
) -> None:
    editor_state_value = session.get("editorState")
    editor_state = (
        editor_state_value if isinstance(editor_state_value, Mapping) else {}
    )
    catalog = editor_state.get("featureCatalog")
    results = session.get("results")
    if (
        not isinstance(catalog, Mapping)
        or catalog.get("schema") != 3
        or not isinstance(catalog.get("items"), list)
        or not isinstance(results, list)
        or len(catalog["items"]) != len(results)
    ):
        raise ValidationError(
            "Legacy color recovery requires a complete schema-3 Feature catalogue."
        )

    global_record_keys: list[str] = []
    for item in catalog["items"]:
        if not isinstance(item, Mapping) or not isinstance(item.get("recordKeys"), list):
            raise ValidationError(
                "Legacy color recovery found incomplete catalogue record identity."
            )
        for raw_record_key in item["recordKeys"]:
            record_key = _normalized_text(raw_record_key)
            if not record_key:
                raise ValidationError(
                    "Legacy color recovery found incomplete catalogue record identity."
                )
            if record_key not in global_record_keys:
                global_record_keys.append(record_key)

    canonical_default_map = dict(canonical_defaults)
    recovered_default_map = dict(recovered_defaults)
    feature_by_alias: dict[str, set[tuple[str, str]]] = {}
    expected_by_feature: dict[tuple[str, str], tuple[str, str]] = {}
    matched_features: set[tuple[str, str]] = set()
    global_biological_keys: set[tuple[str, str]] = set()
    observed_default_sources: dict[str, set[str]] = {}

    def add_alias(alias: Any, feature_key: tuple[str, str]) -> None:
        normalized = _normalized_text(alias)
        if normalized:
            feature_by_alias.setdefault(normalized, set()).add(feature_key)

    for result_index, (item_value, result_value) in enumerate(
        zip(catalog["items"], results, strict=True)
    ):
        if not isinstance(item_value, Mapping) or not isinstance(result_value, Mapping):
            raise ValidationError(
                "Legacy color recovery found ambiguous catalogue Result ownership."
            )
        if (
            item_value.get("resultIndex") != result_index
            or _normalized_text(item_value.get("resultName"))
            != _normalized_text(result_value.get("name"))
            or not isinstance(item_value.get("biologicalFeatures"), list)
            or not isinstance(item_value.get("features"), list)
        ):
            raise ValidationError(
                "Legacy color recovery found ambiguous catalogue Result ownership."
            )

        biological_by_key: dict[tuple[str, str], Mapping[str, Any]] = {}
        for biological in item_value["biologicalFeatures"]:
            if not isinstance(biological, Mapping):
                raise ValidationError(
                    "Legacy color recovery found ambiguous biological Feature identity."
                )
            feature_key = _biological_feature_key(biological)
            if (
                feature_key is None
                or feature_key in biological_by_key
                or feature_key in global_biological_keys
            ):
                raise ValidationError(
                    "Legacy color recovery found ambiguous biological Feature identity."
                )
            biological_by_key[feature_key] = biological
            global_biological_keys.add(feature_key)

        rendered_keys: set[tuple[str, str]] = set()
        for rendered in item_value["features"]:
            if not isinstance(rendered, Mapping):
                raise ValidationError(
                    "Legacy color recovery found an unresolved rendered Feature."
                )
            feature_key = _biological_feature_key(rendered)
            biological = biological_by_key.get(feature_key) if feature_key else None
            svg_id = _normalized_text(rendered.get("svgId"))
            if biological is None or not svg_id or feature_key in rendered_keys:
                raise ValidationError(
                    "Legacy color recovery found an unresolved rendered Feature."
                )
            assert feature_key is not None
            rendered_keys.add(feature_key)
            matched_rule = _resolve_v40_specific_rule(biological, expected_rules)
            catalog_fill = _normalize_color(_normalized_text(rendered.get("fillColor")))
            result_fill = _saved_svg_fill(result_value.get("content"), svg_id)
            if not catalog_fill or catalog_fill != result_fill:
                raise ValidationError(
                    "Legacy color recovery found conflicting catalogue and Result "
                    "Feature fills."
                )
            if matched_rule is not None:
                expected_color = matched_rule[3]
                if catalog_fill != expected_color:
                    raise ValidationError(
                        "Legacy color recovery found a saved Result that conflicts "
                        "with editor rules."
                    )
                matched_features.add(feature_key)
                expected_by_feature[feature_key] = (expected_color, matched_rule[4])
            else:
                feature_type = _normalized_text(biological.get("type"))
                canonical_color = _feature_default_color(
                    canonical_default_map,
                    feature_type,
                )
                recovered_color = _feature_default_color(
                    recovered_default_map,
                    feature_type,
                )
                if catalog_fill == canonical_color and catalog_fill == recovered_color:
                    source = "both"
                elif catalog_fill == canonical_color:
                    source = "canonical"
                elif catalog_fill == recovered_color:
                    source = "recovered"
                else:
                    raise ValidationError(
                        "Legacy color recovery found an unexplained saved default "
                        "Feature fill."
                    )
                observed_default_sources.setdefault(feature_type, set()).add(source)

            add_alias(f"{feature_key[0]}\0{feature_key[1]}", feature_key)
            add_alias(svg_id, feature_key)
            add_alias(biological.get("stableFeatureId"), feature_key)
            try:
                record_index = global_record_keys.index(feature_key[0])
            except ValueError:
                record_index = -1
            source_feature_index = biological.get("sourceFeatureIndex")
            if (
                record_index >= 0
                and type(source_feature_index) is int
                and source_feature_index >= 0
            ):
                add_alias(
                    f"file{record_index}_f{source_feature_index}",
                    feature_key,
                )

    for sources in observed_default_sources.values():
        if "canonical" in sources and "recovered" in sources:
            raise ValidationError(
                "Legacy color recovery found mixed default-color revisions in "
                "saved Results."
            )

    features_value = session.get("features")
    features = features_value if isinstance(features_value, Mapping) else {}
    overrides = features.get("featureColorOverrides")
    if not isinstance(overrides, Mapping):
        if matched_features:
            raise ValidationError(
                "Legacy color recovery is missing its derived Feature override evidence."
            )
        return

    override_features: set[tuple[str, str]] = set()
    for alias, override in overrides.items():
        targets = feature_by_alias.get(str(alias))
        if targets is None or len(targets) != 1 or not isinstance(override, Mapping):
            raise ValidationError(
                "Legacy color recovery found an ambiguous derived Feature override."
            )
        feature_key = next(iter(targets))
        expected = expected_by_feature.get(feature_key)
        actual = (
            _normalize_color(_normalized_text(override.get("color"))),
            _normalized_text(override.get("caption")),
        )
        if (
            expected is None
            or expected != actual
            or feature_key in override_features
        ):
            raise ValidationError(
                "Legacy color recovery found a conflicting derived Feature override."
            )
        override_features.add(feature_key)
    if override_features != matched_features:
        raise ValidationError(
            "Legacy color recovery found incomplete derived Feature override evidence."
        )


def _biological_feature_key(feature: Mapping[str, Any]) -> tuple[str, str] | None:
    record_key = _normalized_text(feature.get("recordKey"))
    biological_feature_id = _normalized_text(feature.get("biologicalFeatureId"))
    if not record_key or not biological_feature_id:
        return None
    return record_key, biological_feature_id


def _feature_default_color(colors: Mapping[str, str], feature_type: str) -> str:
    return _normalize_color(colors.get(feature_type) or colors.get("default") or "")


def _saved_svg_fill(content: Any, svg_id: str) -> str:
    if not svg_id or any(token in svg_id for token in ('"', "<", ">")):
        raise ValidationError(
            "Legacy color recovery found an invalid rendered Feature ID."
        )
    expression = re.compile(
        rf'<[^>]*\sid="{re.escape(svg_id)}"[^>]*>',
        re.IGNORECASE,
    )
    matches = expression.findall(str(content or ""))
    if len(matches) != 1:
        raise ValidationError(
            f"Legacy color recovery could not resolve saved Feature {svg_id!r}."
        )
    fill_match = re.search(r'\bfill="([^"]+)"', matches[0], re.IGNORECASE)
    if fill_match is None:
        raise ValidationError(
            f"Legacy color recovery found a saved Feature without fill: {svg_id!r}."
        )
    return _normalize_color(fill_match.group(1))


def _resolve_v40_specific_rule(
    feature: Mapping[str, Any],
    rules: tuple[tuple[str, str, str, str, str], ...],
) -> tuple[str, str, str, str, str] | None:
    feature_type = _normalized_text(feature.get("type"))
    for candidate_type in (feature_type, "*"):
        if not candidate_type:
            continue
        for priority in range(5):
            for rule in rules:
                if (
                    rule[0] == candidate_type
                    and _v40_rule_priority(rule[1]) == priority
                    and _v40_rule_matches_feature(feature, rule)
                ):
                    return rule
    return None


def _v40_rule_priority(qualifier: str) -> int:
    if qualifier in {
        FEATURE_INSTANCE_HASH_QUALIFIER,
        FEATURE_SEMANTIC_SCOPE_QUALIFIER,
    }:
        return 0
    normalized = qualifier.lower()
    if normalized == "hash":
        return 1
    if normalized == "record_location":
        return 2
    if normalized == "location":
        return 4
    return 3


def _v40_rule_matches_feature(
    feature: Mapping[str, Any],
    rule: tuple[str, str, str, str, str],
) -> bool:
    feature_type, qualifier, pattern, _, _ = rule
    if feature_type not in {_normalized_text(feature.get("type")), "*"}:
        return False
    if not pattern:
        return False
    if qualifier in {
        FEATURE_INSTANCE_HASH_QUALIFIER,
        FEATURE_SEMANTIC_SCOPE_QUALIFIER,
    }:
        return _normalized_text(feature.get("instanceHash")) == pattern

    normalized_qualifier = qualifier.lower()
    if normalized_qualifier == "hash":
        selector_value = feature.get("selector")
        selector = selector_value if isinstance(selector_value, Mapping) else {}
        candidates = (
            selector.get("hash"),
            feature.get("hash"),
            feature.get("stableFeatureId"),
            feature.get("stable_feature_id"),
            feature.get("renderedFeatureSvgId"),
            feature.get("rendered_feature_svg_id"),
            feature.get("svgId"),
            feature.get("svg_id"),
        )
        return _regex_matches_any(pattern, candidates)
    if normalized_qualifier == "record_location":
        return _regex_matches_any(pattern, (_v40_record_location(feature),))
    if normalized_qualifier == "location":
        return _regex_matches_any(pattern, (_v40_feature_location(feature),))

    qualifiers_value = feature.get("qualifiers")
    qualifiers = qualifiers_value if isinstance(qualifiers_value, Mapping) else {}
    values: Any = None
    for key, candidate_values in qualifiers.items():
        if str(key).lower() == normalized_qualifier:
            values = candidate_values
            break
    if values is None:
        for key, candidate_value in feature.items():
            if str(key).lower() == normalized_qualifier:
                values = candidate_value
                break
    candidates = values if isinstance(values, list) else [values]
    return _regex_matches_any(pattern, candidates)


def _regex_matches_any(pattern: str, candidates: Sequence[Any]) -> bool:
    try:
        expression = re.compile(pattern, re.IGNORECASE)
    except re.error:
        return False
    return any(
        expression.search(str(candidate if candidate is not None else "")) is not None
        for candidate in candidates
    )


def _v40_feature_location(feature: Mapping[str, Any]) -> str:
    selector_value = feature.get("selector")
    selector = selector_value if isinstance(selector_value, Mapping) else {}
    for value in (selector.get("location"), feature.get("location")):
        normalized = _normalized_text(value)
        if normalized:
            return normalized
    start = feature.get("start")
    end = feature.get("end")
    return f"{start}..{end}" if start is not None and end is not None else ""


def _v40_record_location(feature: Mapping[str, Any]) -> str:
    selector_value = feature.get("selector")
    selector = selector_value if isinstance(selector_value, Mapping) else {}
    for value in (
        selector.get("record_location"),
        selector.get("recordLocation"),
        feature.get("record_location"),
        feature.get("recordLocation"),
    ):
        normalized = _normalized_text(value)
        if normalized:
            return normalized
    record_id = next(
        (
            normalized
            for normalized in (
                _normalized_text(feature.get("record_id")),
                _normalized_text(feature.get("recordId")),
                _normalized_text(feature.get("record")),
            )
            if normalized
        ),
        "",
    )
    location = _v40_feature_location(feature)
    strand_value = _normalized_text(feature.get("strand")).lower()
    if strand_value in {"1", "positive", "plus", "forward"}:
        strand = "+"
    elif strand_value in {"-1", "negative", "minus", "reverse"}:
        strand = "-"
    else:
        strand = strand_value
    return f"{record_id}:{location}:{strand}" if record_id and location and strand else ""


def _normalized_text(value: Any) -> str:
    return str(value if value is not None else "").strip()


def _normalize_editor_rules(
    value: Any,
) -> tuple[tuple[str, str, str, str, str], ...] | None:
    if value is None:
        return None
    if not isinstance(value, list):
        raise ValidationError("Session config.rules must be an array when present.")
    result: list[tuple[str, str, str, str, str]] = []
    for index, rule in enumerate(value):
        if not isinstance(rule, Mapping):
            raise ValidationError(f"Session config.rules[{index}] must be an object.")
        row = (
            rule.get("feat"),
            rule.get("qual"),
            rule.get("val"),
            rule.get("color"),
            rule.get("cap", ""),
        )
        result.append(_normalize_specific_row(row, context=f"config.rules[{index}]"))
    return tuple(result)


def _normalize_editor_defaults(
    value: Any,
) -> tuple[tuple[str, str], ...] | None:
    if value is None:
        return None
    if not isinstance(value, Mapping):
        raise ValidationError(
            "Saved applied default colors must be an object when present."
        )
    normalized_by_type = {
        feature_type: color
        for feature_type, color in (
            _normalize_default_row((feature_type, color), context="applied colors")
            for feature_type, color in value.items()
        )
        if re.fullmatch(r"collinear_block_\d+", feature_type) is None
    }
    _normalize_comparison_default_colors(normalized_by_type)
    normalized = list(normalized_by_type.items())
    normalized.sort(key=lambda row: row[0])
    if len({row[0] for row in normalized}) != len(normalized):
        raise ValidationError("Saved applied default colors contain duplicate types.")
    return tuple(normalized)


def _normalize_comparison_default_colors(colors: dict[str, str]) -> None:
    plus_max = colors.get("collinear_block_plus_max")
    if plus_max and "collinear_block_plus" not in colors:
        colors["collinear_block_plus"] = plus_max
    for key, color in _COMPARISON_DEFAULT_COLORS.items():
        colors.setdefault(key, color)


def _normalize_specific_row(
    row: Sequence[Any],
    *,
    context: str,
) -> tuple[str, str, str, str, str]:
    feature_type, qualifier, value, color, caption = (
        _clean_color_cell(cell, context=context) for cell in row
    )
    if not feature_type or not qualifier or not value or not color:
        raise ValidationError(f"{context} contains an incomplete specific-color rule.")
    return feature_type, qualifier, value, _normalize_color(color), caption


def _normalize_default_row(
    row: Sequence[Any],
    *,
    context: str,
) -> tuple[str, str]:
    feature_type, color = (_clean_color_cell(cell, context=context) for cell in row)
    if not feature_type or not color:
        raise ValidationError(f"{context} contains an incomplete default-color row.")
    return feature_type, _normalize_color(color)


def _clean_color_cell(value: Any, *, context: str) -> str:
    if value is None:
        return ""
    if not isinstance(value, str):
        raise ValidationError(f"{context} color-table values must be strings.")
    return value.strip()


def _normalize_color(value: str) -> str:
    normalized = value.strip().lower()
    if normalized and _COLOR_PATTERN.fullmatch(normalized) is None:
        raise ValidationError(
            f"Legacy color recovery found an invalid color value: {value!r}."
        )
    return normalized


def _migrate_legacy_multiline_conservation_labels(
    payload: dict[str, Any],
) -> None:
    """Expand the released v39 Web writer's compact multiline label field."""

    diagram_options = payload.get("diagramOptions")
    if not isinstance(diagram_options, Mapping):
        return
    labels = diagram_options.get("conservationLabels")
    blast_files = diagram_options.get("conservationBlastFiles")
    if (
        not isinstance(labels, list)
        or len(labels) != 1
        or not isinstance(labels[0], str)
        or not isinstance(blast_files, list)
        or len(blast_files) <= 1
    ):
        return
    expanded = [line.strip() for line in labels[0].splitlines() if line.strip()]
    if len(expanded) != len(blast_files):
        return
    migrated_options = dict(diagram_options)
    migrated_options["conservationLabels"] = expanded
    payload["diagramOptions"] = migrated_options


def _migrate_legacy_linear_track_slot(slot: Any) -> Any:
    if isinstance(slot, str):
        return _migrate_legacy_linear_track_slot_spec(slot)
    if not isinstance(slot, Mapping):
        return slot

    migrated = copy.deepcopy(dict(slot))
    if (
        migrated.get("kind") == "linearTrackSlot"
        and str(migrated.get("renderer", "")).strip().lower() == "features"
    ):
        if "height" in migrated:
            migrated["height"] = None
        if "spacing" in migrated:
            migrated["spacing"] = None
    return migrated


def _migrate_legacy_linear_track_slot_spec(spec: str) -> str:
    head, separator, options = spec.partition("@")
    if ":" not in head:
        return spec
    renderer = head.split(":", 1)[1].strip().lower()
    if not separator:
        return spec

    parsed_parts: list[tuple[str, str]] = []
    has_legacy_geometry = False
    for part in options.split(","):
        stripped = part.strip()
        if not stripped:
            parsed_parts.append(("", part))
            continue
        if "=" not in stripped:
            return spec
        raw_key, raw_value = stripped.split("=", 1)
        key = raw_key.strip().lower()
        parsed_parts.append((key, part))
        if key in {"h", "height", "spacing"}:
            has_legacy_geometry = True
        if key in {"renderer", "type"}:
            renderer = raw_value.strip().lower()

    if renderer != "features" or not has_legacy_geometry:
        return spec
    retained = [
        part
        for key, part in parsed_parts
        if key not in {"h", "height", "spacing"}
    ]
    if not any(part.strip() for part in retained):
        return head
    return f"{head}@{','.join(retained)}"


__all__ = [
    "AdaptedSessionRequest",
    "SessionCompatibleRequestRenderResult",
    "SessionMigrationReport",
    "adapt_session_request",
    "build_session_compatible_request_diagram",
    "canonical_payload_for_session_decode",
    "render_session_compatible_request",
    "rewrite_protein_artifact_references",
]
