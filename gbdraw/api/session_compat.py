"""Compatibility adapter from persisted sessions to current typed rendering."""

from __future__ import annotations

import copy
import json
import re
from dataclasses import dataclass, field, fields, is_dataclass, replace
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
_LEGACY_PROTEIN_REFERENCE_RE = re.compile(
    r"p_[A-Za-z0-9._%+-]+?_\d+_\d+_(?:-1|0|1)_[0-9a-f]{12}"
    r"(?:_[2-9][0-9]*)?"
)
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
