"""Compatibility adapter from persisted sessions to current typed rendering."""

from __future__ import annotations

import base64
import copy
import hashlib
import json
import logging
import re
from dataclasses import dataclass, fields, is_dataclass, replace
from functools import lru_cache
from pathlib import Path
from typing import Any, Callable, Collection, Iterable, Mapping, NamedTuple, Sequence, cast

import numpy as np
from pandas import DataFrame

from gbdraw.analysis.protein_colinearity import (
    ProteinBlastpMode,
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
from gbdraw.core.record_metadata import _iter_source_features, _source_feature_index
from gbdraw.exceptions import ValidationError
from gbdraw.features.ids import compute_feature_hash
from gbdraw.layout.similarity_alignment import (
    AlignmentAnchorIdentity,
    AlignmentDecisionStatus,
    AlignmentRecordDecision,
    AlignmentResolutionRationale,
    SimilarityAlignmentPlan,
)
from gbdraw.session_drawings import SessionDrawingArtifacts, drawing_draft_config
from gbdraw.session_io import (
    DRAWN_HASH_RULE_SESSION_MAX_VERSION,
    FEATURE_CATALOG_SCHEMA_BY_SESSION_VERSION,
    DrawnFeature,
    SessionResourceDescriptor,
    classify_raw_losat_cache_entry,
    empty_protein_identity_manifest,
    source_hash_selector_value,
    unmapped_hash_selector_warning,
)

from .migration_report import SessionMigrationReport
from .options import LinearMultiRecordOptions, LinearRecordTranslation, losatp_analysis_mode
from .record_planning import (
    _SOURCE_BINDING_COLUMNS,
    ResolvedRecordCollection,
    ResolvedRecordProvenance,
    UnmatchedFrameItem,
    comparison_record_bindings,
    has_source_bound_results,
    project_similarity_alignment_centers,
    project_source_bound_results,
)

from .request_render import (
    CircularBatchRenderResult,
    CircularBatchRequestPlan,
    CircularRequestPlan,
    CurrentRequestArtifacts,
    DiagramRequestPlan,
    LinearRequestPlan,
    PreparedCircularBatchRequest,
    PreparedDiagramRequest,
    RequestRenderResult,
    _extract_linear_request_proteins,
    _build_request_plan_diagram,
    _request_losatp_job_specs,
    plan_request,
    render_prepared_request,
)
from .requests import (
    CircularBatchRequest,
    CircularDiagramRequest,
    DiagramRequest,
    LinearDiagramRequest,
    _LegacySimilarityAlignment,
    _with_legacy_similarity_alignment,
)
from ..linear_comparison import (
    LinearComparison,
    _endpoint_frame,
    project_search_frame_comparisons,
    reverse_endpoint_table_text,
)

_LEGACY_LINEAR_TRACK_SLOT_SESSION_VERSION = 32
_LEGACY_MULTILINE_CONSERVATION_LABEL_SESSION_VERSION = 39
_LEGACY_PROTEIN_REFERENCE_RE = re.compile(
    r"p_[A-Za-z0-9._%+-]+?_\d+_\d+_(?:-1|0|1)_[0-9a-f]{12}"
    r"(?:_[2-9][0-9]*)?"
)
_FEATURE_ANALYSIS_REFERENCE_RE = re.compile(r"f_[0-9a-f]{64}")

logger = logging.getLogger(__name__)


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
    protein_source_mode: ProteinBlastpMode | None
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


def _session_protein_mode(artifacts: Mapping[str, Any]) -> ProteinBlastpMode | None:
    # The LOSATP settings are Linear settings.
    losat = drawing_draft_config(artifacts, "linear").get("losat")
    blastp = losat.get("blastp") if isinstance(losat, Mapping) else None
    configured = blastp.get("mode") if isinstance(blastp, Mapping) else None
    if str(configured or "") in {"pairwise", "orthogroup", "collinear"}:
        return cast(ProteinBlastpMode, str(configured))
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

    manifest = artifacts.get("proteinIdentityManifest")
    if manifest is not None and not isinstance(manifest, Mapping):
        raise ValidationError("Session proteinIdentityManifest must be an object.")
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


def written_legacy_protein_artifacts(
    candidates: Sequence[Mapping[str, Any]],
    evidence: Sequence[Mapping[str, Any]],
    *,
    promoted: Collection[int],
    has_current_derived: bool,
) -> tuple[tuple[Mapping[str, Any], ...], tuple[Mapping[str, Any], ...]]:
    """The legacy protein raw candidates and derived evidence a written Session keeps.

    ``candidates`` are in their states after this write's attempts and
    ``promoted`` are the indexes it promoted: a promoted candidate leaves the
    envelope (its entry is a current raw entry now); every other candidate
    stays, ``pending`` or ``rejected``. The evidence goes when this write
    promoted a candidate and left none pending, or when the written Session
    holds or built current derived entries (OV-1129); else it stays. The Web
    twin is ``writtenLegacyProteinArtifacts`` (``services/losat-cache.js``);
    ``tests/fixtures/legacy_protein_artifact_vectors.json`` pins both
    (OV-1132, OV-1143).
    """

    kept = tuple(
        candidate for index, candidate in enumerate(candidates) if index not in promoted
    )
    superseded = has_current_derived or (
        bool(promoted) and not any(candidate.get("state") == "pending" for candidate in kept)
    )
    return kept, () if superseded else tuple(evidence)


def _legacy_candidate_attempt_pairs(
    request: LinearDiagramRequest,
    stored_result_mode: ProteinBlastpMode | None,
) -> tuple[tuple[int, int], ...]:
    """The directed record pairs whose legacy candidates this drawing's protein run tries.

    They are the searches of the request's LOSATP run. A request that stores
    its protein results (mode ``"none"``) is computed again in
    ``stored_result_mode``, the Web draft's mode, else the stored result's, as
    the next Web Generate does. A request without a LOSATP search tries none.
    """

    search = request.options.losat_search
    if search is None or search.program != "losatp":
        return ()
    mode = losatp_analysis_mode(search)
    if mode == "none":
        options = request.options
        mode = stored_result_mode or (
            "collinear"
            if options.collinearity_blocks is not None
            else "orthogroup" if options.orthogroups is not None else "pairwise"
        )
    return _request_losatp_job_specs(request, mode, search.pairs)


def _promote_legacy_candidates(
    candidates: Sequence[Mapping[str, Any]],
    extraction: ProteinExtractionResult,
    pairs: Sequence[tuple[int, int]],
) -> tuple[
    tuple[Mapping[str, Any], ...],
    frozenset[int],
    Mapping[str, str],
]:
    """Try each ``pending`` candidate on the searched record ``pairs``.

    Returns the promoted current raw entries, the indexes of the candidates
    they came from and the legacy-to-current protein ID map. An unmatched
    candidate stays as it is; a ``rejected`` one is not tried again.
    """

    manifest = extraction.identity_manifest
    if manifest is None:
        raise ValidationError("Protein cache promotion requires an identity manifest.")
    fastas = tuple(proteins_to_fasta(items) for items in extraction.proteins_by_record)
    promoted_entries: list[Mapping[str, Any]] = []
    promoted_indexes: set[int] = set()
    id_map: dict[str, str] = {}

    for candidate_index, candidate in enumerate(candidates):
        if candidate.get("state") != "pending":
            continue
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
        for query_index, subject_index in pairs:
            try:
                scan = promote_legacy_protein_raw_cache_entries(
                    (candidate,),
                    query_proteins=extraction.proteins_by_record[query_index],
                    subject_proteins=extraction.proteins_by_record[subject_index],
                    query_fasta=fastas[query_index],
                    subject_fasta=fastas[subject_index],
                    identity_manifest=manifest,
                    expected_args=expected_args,
                    expected_program=str(original.get("program") or "blastp"),
                    expected_outfmt=str(original.get("outfmt") or "6"),
                )
            except ValidationError:
                continue
            if scan.promotion is not None:
                promotion = scan.promotion
                break

        if promotion is None:
            continue
        if not validate_protein_raw_entry_references(promotion.entry, manifest):
            raise ValidationError(
                "Promoted protein raw entry failed manifest validation."
            )
        promoted_entries.append(promotion.entry)
        promoted_indexes.add(candidate_index)
        _merge_promoted_tsv_id_map(
            id_map,
            legacy_text=str(original.get("text") or ""),
            current_text=promotion.rewritten_tsv,
        )

    return tuple(promoted_entries), frozenset(promoted_indexes), id_map


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
                if field.init
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


class _ProteinReferences(NamedTuple):
    """The protein references in request artifacts: the legacy protein IDs, and
    whether a feature analysis ID is among them."""

    legacy_ids: set[str]
    feature_analysis: bool


# Exact types that hold no string and no item, so the scan skips them.
_REFERENCE_FREE_TYPES = frozenset(
    {int, float, bool, type(None), np.float64, np.int64, np.bool_}
)


@lru_cache(maxsize=256)
def _init_field_names(cls: type) -> tuple[str, ...]:
    return tuple(field.name for field in fields(cls) if field.init)


def _protein_references(values: Sequence[object]) -> _ProteinReferences:
    """Scan strings, table cells, init fields, mapping keys and items, and
    collection items of ``values`` once."""

    legacy_ids: set[str] = set()
    feature_analysis = False

    def visit(value: object) -> None:
        nonlocal feature_analysis
        value_type = type(value)
        if value_type in _REFERENCE_FREE_TYPES:
            return
        items: Iterable[object]
        if isinstance(value, str):
            # Each pattern starts with its prefix, so a string without it has no match.
            if "p_" in value:
                legacy_ids.update(
                    match.group(0)
                    for match in _LEGACY_PROTEIN_REFERENCE_RE.finditer(value)
                )
            if not feature_analysis and "f_" in value:
                feature_analysis = _FEATURE_ANALYSIS_REFERENCE_RE.search(value) is not None
            return
        if isinstance(value, DataFrame):
            items = value.to_numpy().ravel()
        elif is_dataclass(value) and not isinstance(value, type):
            items = (getattr(value, name) for name in _init_field_names(type(value)))
        elif isinstance(value, Mapping):
            items = (part for pair in value.items() for part in pair)
        elif isinstance(value, (tuple, list, set, frozenset)):
            items = value
        else:
            return
        for item in items:
            visit(item)

    for value in values:
        visit(value)
    return _ProteinReferences(legacy_ids, feature_analysis)


def _request_protein_artifacts(
    request: LinearDiagramRequest,
) -> tuple[Any, ...]:
    options = request.options
    return (
        options.linear_comparisons,
        options.protein_comparisons,
        options.orthogroups,
        options.collinearity_blocks,
    )


@dataclass(frozen=True)
class _LegacyAlignmentMember:
    anchor: AlignmentAnchorIdentity
    representative: bool
    aliases: frozenset[str]


def _mapping_items(value: object) -> tuple[Mapping[str, Any], ...]:
    if not isinstance(value, list):
        return ()
    return tuple(item for item in value if isinstance(item, Mapping))


def _drawn_features(
    record: Any, provenance: ResolvedRecordProvenance, ordinal: int
) -> list[DrawnFeature]:
    """The features of a planned record, the ``ordinal``-th of its Result.

    Release 0.13.0 and Sessions up to 44 hashed the drawn feature, so a cropped
    or reverse-complemented record saved hashes of its drawn coordinates; the
    source feature catalog hashes source coordinates (OV-401). Only these
    legacy readers compute the drawn hash.
    """

    catalog = provenance.source_features.source_hashes
    drawn: list[DrawnFeature] = []
    # A feature's source index, or its ordinal in an untransformed record, as
    # capture_source_features reads it.
    for position, feature in enumerate(_iter_source_features(record.features)):
        index = _source_feature_index(feature)
        index = position if index is None else index
        if index in catalog:
            drawn.append(
                (catalog[index], compute_feature_hash(feature, record_id=record.id), index, ordinal)
            )
    return drawn


def _legacy_display_frame_feature_id(
    plan: LinearRequestPlan | None,
    record_key: str,
    feature_id: str,
) -> str:
    """Bind a release-0.13.0 display-frame feature ID to its source-frame ID.

    An ID that already names a catalog feature is returned as is.
    """

    if plan is None:
        return feature_id
    for record, provenance in zip(plan.records, plan.provenance, strict=True):
        if provenance.record_key != record_key:
            continue
        if feature_id in provenance.source_features.source_hashes.values():
            return feature_id
        matches = {
            source
            for source, drawn_hash, _, _ in _drawn_features(record, provenance, 0)
            if drawn_hash == feature_id
        }
        return next(iter(matches)) if len(matches) == 1 else feature_id
    return feature_id


# A Session 40-44's feature catalog as stored (schema 3) or as the Web reader
# sends it, promoted to a later schema (fields added; ids and groups kept).
_LEGACY_ALIGNMENT_CATALOG_SCHEMAS = frozenset(
    {3, *FEATURE_CATALOG_SCHEMA_BY_SESSION_VERSION.values()}
)


def _legacy_alignment_groups(
    request: LinearDiagramRequest,
    session_artifacts: Mapping[str, Any],
    plan: LinearRequestPlan | None = None,
) -> tuple[tuple[str, tuple[_LegacyAlignmentMember, ...]], ...]:
    record_keys = tuple(record.record_key for record in request.records)
    if any(record_key is None for record_key in record_keys):
        raise ValidationError(
            "Legacy similarity alignment requires stable displayed record keys."
        )
    stable_keys = tuple(str(record_key) for record_key in record_keys)
    editor_state = session_artifacts.get("editorState")
    catalog = editor_state.get("featureCatalog") if isinstance(editor_state, Mapping) else None
    items = (
        _mapping_items(catalog.get("items"))
        if isinstance(catalog, Mapping)
        and catalog.get("schema") in _LEGACY_ALIGNMENT_CATALOG_SCHEMAS
        else ()
    )
    biological = tuple(
        feature
        for item in items
        for feature in _mapping_items(item.get("biologicalFeatures"))
    )
    rendered = tuple(
        feature
        for item in items
        for feature in _mapping_items(item.get("features"))
    )
    catalog_groups = tuple(
        group
        for item in items
        for group in _mapping_items(item.get("orthogroups"))
    )
    state = session_artifacts.get("orthogroupState")
    groups = catalog_groups or (
        _mapping_items(state.get("groups")) if isinstance(state, Mapping) else ()
    )
    use_legacy_groups = not catalog_groups
    result: list[tuple[str, tuple[_LegacyAlignmentMember, ...]]] = []
    for group in groups:
        group_id = group.get("id")
        if not isinstance(group_id, str) or not group_id.strip():
            raise ValidationError("Legacy similarity alignment group has no stable ID.")
        members: list[_LegacyAlignmentMember] = []
        for member in _mapping_items(group.get("members")):
            if use_legacy_groups:
                record_index = member.get("recordIndex")
                if (
                    isinstance(record_index, bool)
                    or not isinstance(record_index, int)
                    or record_index < 0
                    or record_index >= len(stable_keys)
                ):
                    raise ValidationError(
                        "Legacy similarity alignment member has no stable record mapping."
                    )
                record_key = stable_keys[record_index]
                stable_id = next(
                    (
                        value.strip()
                        for name in (
                            "stableFeatureSvgId",
                            "stable_feature_svg_id",
                            "featureSvgId",
                        )
                        if isinstance((value := member.get(name)), str)
                        and value.strip()
                    ),
                    None,
                )
                feature_index = member.get("featureIndex")
                if stable_id is None or not isinstance(feature_index, int) or isinstance(feature_index, bool) or feature_index < 0:
                    raise ValidationError(
                        "Legacy similarity alignment member lacks stable feature metadata."
                    )
                stable_id = _legacy_display_frame_feature_id(plan, record_key, stable_id)
                anchor = AlignmentAnchorIdentity(
                    record_key,
                    stable_id,
                    None,
                    stable_id,
                )
            else:
                # Checked against stable_keys (all str) before use below.
                record_key = cast(str, member.get("recordKey"))
                biological_id = member.get("biologicalFeatureId")
                matches = tuple(
                    feature
                    for feature in biological
                    if feature.get("recordKey") == record_key
                    and feature.get("biologicalFeatureId") == biological_id
                )
                if record_key not in stable_keys or not isinstance(biological_id, str) or len(matches) != 1:
                    raise ValidationError(
                        "Legacy similarity alignment member lacks unique biological identity metadata."
                    )
                source_index = matches[0].get("sourceFeatureIndex")
                stable_feature_id = matches[0].get("stableFeatureId")
                anchor = AlignmentAnchorIdentity(
                    record_key,
                    biological_id,
                    source_index if isinstance(source_index, int) and not isinstance(source_index, bool) and source_index >= 0 else None,
                    # The catalog's stable feature ID, else the biological ID: the
                    # anchor the Web's legacy promotion wrote before Python converted.
                    stable_feature_id
                    if isinstance(stable_feature_id, str) and stable_feature_id.strip()
                    else biological_id,
                )
            rendered_ids = {
                str(feature.get("svgId")).strip()
                for feature in rendered
                if feature.get("recordKey") == anchor.record_key
                and feature.get("biologicalFeatureId") == anchor.biological_feature_id
                and feature.get("svgId")
            }
            aliases = frozenset(
                {
                    str(value).strip()
                    for value in (
                        member.get("biologicalFeatureId"),
                        member.get("stableFeatureSvgId"),
                        member.get("stable_feature_svg_id"),
                        member.get("featureSvgId"),
                        member.get("sourceProteinId"),
                        member.get("proteinId"),
                        member.get("label"),
                    )
                    if value is not None and str(value).strip()
                }
                | rendered_ids
            )
            members.append(
                _LegacyAlignmentMember(
                    anchor,
                    member.get("representative") is True,
                    aliases,
                )
            )
        result.append((group_id.strip(), tuple(members)))
    return tuple(result)


def _legacy_alignment_record_member(
    members: Sequence[_LegacyAlignmentMember],
    record_key: str,
) -> _LegacyAlignmentMember | None:
    matching = [member for member in members if member.anchor.record_key == record_key]
    representatives = [member for member in matching if member.representative]
    if len(representatives) == 1:
        return representatives[0]
    if len(representatives) > 1 or len(matching) > 1:
        raise ValidationError(
            "Legacy similarity alignment has ambiguous members for one displayed record."
        )
    return matching[0] if matching else None


def _legacy_similarity_alignment_plan(
    request: LinearDiagramRequest,
    session_artifacts: Mapping[str, Any],
    *,
    target: str,
    plan: LinearRequestPlan | None = None,
) -> SimilarityAlignmentPlan:
    groups = _legacy_alignment_groups(request, session_artifacts, plan)
    group_matches = [group for group in groups if group[0] == target]
    selected_member: _LegacyAlignmentMember | None = None
    if len(group_matches) > 1:
        raise ValidationError(
            "Legacy similarity alignment group is ambiguous in saved metadata."
        )
    if group_matches:
        group_id, members = group_matches[0]
        selected_member = next(
            (member for member in members if member.representative),
            members[0] if members else None,
        )
    else:
        member_matches = [
            (group_id, members, member)
            for group_id, members in groups
            for member in members
            if target in member.aliases
        ]
        if len(member_matches) != 1:
            raise ValidationError(
                "Legacy similarity alignment target is absent or ambiguous in saved stable metadata."
            )
        group_id, members, selected_member = member_matches[0]
    if selected_member is None:
        raise ValidationError(
            "Legacy similarity alignment has no materializable reference member."
        )
    reference = selected_member.anchor
    decisions = []
    for record in request.records:
        assert record.record_key is not None
        if record.record_key == reference.record_key:
            status = AlignmentDecisionStatus.REFERENCE
            rationale = AlignmentResolutionRationale.REFERENCE
            anchor: AlignmentAnchorIdentity | None = reference
        else:
            member = _legacy_alignment_record_member(members, record.record_key)
            status = (
                AlignmentDecisionStatus.ALIGNED
                if member is not None
                else AlignmentDecisionStatus.SKIPPED
            )
            rationale = (
                AlignmentResolutionRationale.ONLY_USABLE_CANDIDATE
                if member is not None
                else AlignmentResolutionRationale.SKIPPED_NO_CANDIDATE
            )
            anchor = member.anchor if member is not None else None
        decisions.append(
            AlignmentRecordDecision(
                record_key=record.record_key,
                status=status,
                rationale=rationale,
                anchor=anchor,
            )
        )
    return SimilarityAlignmentPlan(
        group_id=group_id,
        reference=reference,
        records=tuple(decisions),
    )


def promote_legacy_session_similarity_alignment_request(
    request: DiagramRequest,
    drawing: SessionDrawingArtifacts,
    *,
    plan: LinearRequestPlan | None = None,
) -> DiagramRequest:
    """Promote one supported old Session alignment into current typed state.

    With the resolved ``plan`` of ``request``, saved feature IDs of displayed
    (reverse-complemented) records also bind to their source features.
    """

    if not isinstance(request, LinearDiagramRequest):
        return request
    session_artifacts = drawing.fields
    orthogroup_state = session_artifacts.get("orthogroupState")
    session_target = (
        orthogroup_state.get("selectedOrthogroupAlignmentFeature")
        if isinstance(orthogroup_state, Mapping)
        else None
    )
    if session_target == "":
        # origin/main's Web Save writes "" when no alignment target is selected.
        session_target = None
    if session_target is not None and (
        not isinstance(session_target, str)
        or not session_target.strip()
        or "\0" in session_target
    ):
        raise ValidationError(
            "Legacy Session similarity alignment target must be non-empty text without NUL."
        )
    session_target = session_target.strip() if isinstance(session_target, str) else None
    legacy = request._legacy_similarity_alignment
    targets = {value for value in (session_target, legacy.target if legacy else None) if value}
    if not targets:
        return request
    if len(targets) != 1:
        raise ValidationError(
            "Legacy Session contains conflicting similarity alignment owners."
        )
    # Session 44 introduced the typed similarity alignment state.
    if drawing.version >= 44:
        raise ValidationError(
            "Current Sessions must store typed similarity alignment state."
        )
    if (
        request.similarity_alignment is not None
        and request.layout is not None
        and request.layout.record_translations
    ):
        if legacy is not None:
            return request
        expected_plan = _legacy_similarity_alignment_plan(
            request,
            session_artifacts,
            target=targets.pop(),
            plan=plan,
        )
        if request.similarity_alignment != expected_plan or any(
            translation.x != 0.0 or translation.y != 0.0
            for translation in request.layout.record_translations
        ):
            raise ValidationError(
                "Legacy Session similarity alignment conflicts with current typed state."
            )
        return _with_legacy_similarity_alignment(
            request,
            _LegacySimilarityAlignment(
                target=session_target or expected_plan.group_id,
                source_schema=int(
                    session_artifacts.get("renderRequest", {}).get("schema", 1)
                ),
            ),
        )
    if request.similarity_alignment is not None or (
        request.layout is not None and request.layout.record_translations
    ):
        raise ValidationError(
            "Legacy Session similarity alignment conflicts with current typed state."
        )
    target = targets.pop()
    alignment = _legacy_similarity_alignment_plan(
        request,
        session_artifacts,
        target=target,
        plan=plan,
    )
    layout = request.layout or LinearMultiRecordOptions()
    promoted = replace(
        request,
        layout=replace(
            layout,
            record_translations=tuple(
                LinearRecordTranslation(record_key=record.record_key)
                for record in request.records
                if record.record_key is not None
            ),
        ),
        similarity_alignment=alignment,
    )
    return _with_legacy_similarity_alignment(
        promoted,
        legacy
        or _LegacySimilarityAlignment(
            target=target,
            source_schema=int(
                session_artifacts.get("renderRequest", {}).get("schema", 1)
            ),
        ),
    )


def project_legacy_similarity_alignment_for_current_write(
    request: DiagramRequest,
    *,
    legacy_source: DiagramRequest | None = None,
) -> DiagramRequest:
    """Remove a materialized runtime string when a typed legacy projection exists."""

    if not isinstance(request, LinearDiagramRequest):
        return request
    source = legacy_source if isinstance(legacy_source, LinearDiagramRequest) else request
    if source._legacy_similarity_alignment is None:
        return request
    if (
        source.similarity_alignment is None
        or source.layout is None
        or not source.layout.record_translations
        or request.similarity_alignment != source.similarity_alignment
    ):
        raise ValidationError(
            "Legacy similarity alignment cannot be projected to a current request."
        )
    return replace(
        request,
        layout=source.layout,
        similarity_alignment=source.similarity_alignment,
    )


def _rewrite_linear_request_protein_references(
    request: LinearDiagramRequest,
    id_map: Mapping[str, str],
    references: _ProteinReferences,
) -> LinearDiagramRequest:
    """Rewrite verified IDs; ``references`` is the scan of ``request``."""

    if not id_map:
        rewritten_request = request
        unresolved = references
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
            ),
        )
        unresolved = _protein_references(_request_protein_artifacts(rewritten_request))
    if unresolved.legacy_ids or unresolved.feature_analysis:
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


def _replace_plan_request(
    plan: DiagramRequestPlan,
    request: DiagramRequest,
) -> DiagramRequestPlan:
    # The replacement request is derived from plan.request, so it has its mode.
    if isinstance(plan, CircularRequestPlan):
        assert isinstance(request, CircularDiagramRequest)
        return replace(plan, request=request)
    if isinstance(plan, CircularBatchRequestPlan):
        assert isinstance(request, CircularBatchRequest)
        return replace(plan, request=request)
    assert isinstance(request, LinearDiagramRequest)
    if plan.request.similarity_alignment == request.similarity_alignment:
        return replace(plan, request=request)
    if plan.request.similarity_alignment is not None:
        raise ValidationError(
            "Session compatibility cannot replace an active similarity alignment plan."
        )
    # Only legacy promotion adds an alignment here, and it adds a resolved plan.
    assert isinstance(request.similarity_alignment, SimilarityAlignmentPlan)
    centers = project_similarity_alignment_centers(
        ResolvedRecordCollection(plan.records, plan.provenance),
        request.similarity_alignment,
    )
    # Promotion adds the record translations of the alignment; keep the
    # placement rows the plan already resolved.
    layout = request.layout
    if layout is not None and plan.layout is not None:
        layout = replace(
            layout, multi_record_positions=plan.layout.multi_record_positions
        )
    return replace(
        plan,
        request=request,
        layout=layout,
        alignment_anchor_centers=centers,
    )


MAIN_DISPLAY_FRAME_SESSION_VERSION = 42
# Sessions up to 47 stored the protein results of a reversed Linear record in
# the drawn frame; Session 48 stores them in the search frame (OV-1145).
DRAWN_FRAME_PROTEIN_RESULT_SESSION_VERSION = 47


def describe_unmatched_frame_item(item: UnmatchedFrameItem) -> str:
    feature = "" if item.feature_index is None else f", source feature {item.feature_index}"
    return (
        f"the stored {item.kind} of protein {item.protein_id!r} (record #{item.record_index + 1}{feature}) "
        "matches neither strand of its record"
    )


def drawn_frame_protein_results_to_search_frame(
    plan: DiagramRequestPlan,
    request: DiagramRequest,
) -> tuple[DiagramRequest, tuple[UnmatchedFrameItem, ...]]:
    """Store the protein results of ``request`` in each record's search frame.

    The core of the Session 47 migrator (OV-1145): every table row and typed
    item tells its frame by its tags, so a tagged item already in the search
    frame stays. An untagged anchor side (Session 39) of a reversed record is
    decided by its record, so only a Session up to 47 may be converted, once.
    An item matching neither frame stays as given and is returned; it draws
    as stored, and a Session writer refuses it.
    """

    if not isinstance(plan, LinearRequestPlan) or not isinstance(request, LinearDiagramRequest):
        return request, ()
    options = request.options
    if not has_source_bound_results(
        options.linear_comparisons, options.protein_comparisons,
        options.orthogroups, options.collinearity_blocks,
    ):
        return request, ()
    projection = project_source_bound_results(
        options,
        comparison_record_bindings(plan.records, plan.provenance, search_frame=True),
    )
    if projection.options is options:
        return request, projection.unmatched
    return replace(request, options=projection.options), projection.unmatched


_STORED_PROTEIN_RESULT_KINDS = frozenset(
    {"precomputedProteinComparison", "orthogroupResult", "collinearityResult"}
)


def stores_protein_results(payload: Mapping[str, Any] | None) -> bool:
    """Whether a stored Linear request payload holds protein results."""

    comparisons = payload.get("comparisons") if isinstance(payload, Mapping) else None
    return isinstance(payload, Mapping) and payload.get("mode") == "linear" and any(
        isinstance(item, Mapping) and item.get("kind") in _STORED_PROTEIN_RESULT_KINDS
        for item in (comparisons if isinstance(comparisons, list) else ())
    )


def search_frame_protein_resources(
    payload: Mapping[str, object],
    request: DiagramRequest,
    converted: DiagramRequest,
) -> dict[str, bytes]:
    """The new bytes of the stored protein resources of ``payload`` that ``converted`` changed.

    ``request`` is ``payload`` decoded and ``converted`` its results in the
    search frame (OV-1145); see :func:`_search_frame_items`. Each changed
    item's bytes go under its stored resource ID, so the request and every
    other resource stay as stored; a resource whose items changed to
    different bytes is refused.
    """

    changed: dict[str, bytes] = {}
    for item in _search_frame_items(payload, request, converted, _STORED_PROTEIN_RESULT_KINDS):
        if item.content is None:
            continue
        if changed.get(item.resource_id, item.content) != item.content:
            raise ValidationError(f"Stored resource {item.resource_id!r} holds two different comparison results.")
        changed[item.resource_id] = item.content
    return changed


_STORED_NUCLEOTIDE_TABLE_KINDS = frozenset({"nucleotideBlast"})


@dataclass(frozen=True)
class SearchFrameNucleotideTables:
    """The resources of a main Session with its nucleotide tables in the search frame.

    ``pair_resources`` names, by the (query, subject) record keys of a
    request item, the resource that holds that item's rows when it is not
    the item's stored resource; the Web comparison card of that record pair
    names it too.
    """

    resources: dict[str, SessionResourceDescriptor]
    pair_resources: dict[tuple[str, str], str]


def search_frame_nucleotide_resources(
    payload: Mapping[str, object],
    request: DiagramRequest,
    converted: DiagramRequest,
    resources: Mapping[str, SessionResourceDescriptor],
) -> SearchFrameNucleotideTables:
    """``resources`` with the stored nucleotide tables of ``payload`` in the search frame.

    ``request`` is ``payload`` decoded and ``converted`` its rows in the
    search frame (a main Session up to 42, PD-OI-073). A stored table keeps
    its resource ID for the items whose rows did not change, or, when every
    item that names it changed, for the first item by record pair, which
    then holds its converted rows. Each other distinct converted table gets
    a new resource as the request encoder allocates one, so the request item
    and the Web comparison card of every record pair read that pair's rows
    (OV-2080, OV-2097).
    """

    from gbdraw.session_resources import SessionResourceTable

    by_resource: dict[str, list[_SearchFrameItem]] = {}
    for item in _search_frame_items(payload, request, converted, _STORED_NUCLEOTIDE_TABLE_KINDS):
        by_resource.setdefault(item.resource_id, []).append(item)
    tables = dict(resources)
    moved: list[tuple[_SearchFrameItem, bytes]] = []
    for resource_id, items in by_resource.items():
        items.sort(key=lambda item: item.record_pair)
        kept = None if any(item.content is None for item in items) else items[0].content
        if kept is not None:
            tables[resource_id] = replace_resource_bytes(resources[resource_id], kept)
        moved.extend((item, item.content) for item in items if item.content not in (None, kept))
    table = SessionResourceTable(tables)
    pair_resources: dict[tuple[str, str], str] = {}
    for item, content in moved:
        stored = resources[item.resource_id]
        resource_id = table.request().add_bytes(
            item.resource_id, kind=stored.get("kind", ""), name=stored.get("name", ""), content=content,
        )
        if resource_id not in tables:
            tables[resource_id] = replace_resource_bytes(
                {**stored, "name": table.descriptors()[resource_id]["name"]}, content,
            )
        query_key, subject_key = (request.records[index].record_key for index in item.record_pair)
        if query_key is not None and subject_key is not None:
            pair_resources[(query_key, subject_key)] = resource_id
    return SearchFrameNucleotideTables(tables, pair_resources)


@dataclass(frozen=True)
class _SearchFrameItem:
    """One stored comparison result: its stored resource and the bytes the conversion gave it.

    ``content`` is ``None`` when the conversion kept the stored bytes;
    ``record_pair`` is the (query, subject) record indices of a
    record-pair result, ``(-1, -1)`` for one of the request.
    """

    resource_id: str
    record_pair: tuple[int, int]
    content: bytes | None


def _search_frame_items(
    payload: Mapping[str, object],
    request: DiagramRequest,
    converted: DiagramRequest,
    kinds: frozenset[str],
) -> list[_SearchFrameItem]:
    """The stored ``kinds`` items of ``payload``, each with its bytes in ``converted``.

    Both requests are encoded alone and their items meet the stored ones by
    identity (:func:`_comparison_items`).
    """

    stored = _comparison_items(payload.get("comparisons"), kinds)
    if not stored:
        return []

    from gbdraw.session_request_codec import encode_canonical_request
    from gbdraw.session_resources import SessionResourceTable
    from gbdraw.session_io import _embedded_resource_bytes

    def items_of(
        value: DiagramRequest,
    ) -> tuple[dict[tuple[object, ...], Mapping[str, object]], dict[str, SessionResourceDescriptor]]:
        table = SessionResourceTable()
        encoded = encode_canonical_request(value, table=table).payload
        return _comparison_items(encoded.get("comparisons"), kinds), table.descriptors()

    # ``request`` as stored: its tables encode to the bytes they were read from.
    as_stored = (_with_drawn_feature_id_marks(request, lambda _comparison: False)
                 if isinstance(request, LinearDiagramRequest) else request)
    (before, before_resources), (after, after_resources) = items_of(as_stored), items_of(converted)
    if not stored.keys() == before.keys() == after.keys():
        raise ValidationError("Stored comparison results do not match their decoded request.")
    items: list[_SearchFrameItem] = []
    for identity, item in stored.items():
        old, new = before[identity], after[identity]
        content = _embedded_resource_bytes(after_resources[_comparison_resource_id(new)])
        query, subject = (identity[1], identity[2]) if len(identity) == 4 else (-1, -1)
        if not isinstance(query, int) or not isinstance(subject, int):
            raise ValidationError("A stored comparison result requires integer record indices.",
                                  diagnostic={"code": "INPUT_INVALID"})
        items.append(_SearchFrameItem(
            _comparison_resource_id(item),
            (query, subject),
            None if content == _embedded_resource_bytes(before_resources[_comparison_resource_id(old)]) else content,
        ))
    return items


_RECORD_PAIR_RESULT_KINDS = frozenset({"nucleotideBlast", "precomputedProteinComparison"})


def _comparison_items(comparisons: object, kinds: frozenset[str]) -> dict[tuple[object, ...], Mapping[str, object]]:
    """The ``kinds`` items of a stored or encoded request's ``comparisons``, by identity.

    An item is its kind and record pair: the stored indices, or under request
    schema 1 (no indices) the k-th item of its kind compares records k and
    k + 1, as the decoder reads it. Orthogroup and collinearity results are
    one per request. A repeated identity is told apart by its occurrence.
    """

    items: dict[tuple[object, ...], Mapping[str, object]] = {}
    per_kind: dict[object, int] = {}
    for item in comparisons if isinstance(comparisons, list) else ():
        if not isinstance(item, Mapping) or item.get("kind") not in kinds:
            continue
        kind = item.get("kind")
        position = per_kind[kind] = per_kind.get(kind, -1) + 1
        pair: tuple[object, ...] = (
            (item.get("queryRecordIndex"), item.get("subjectRecordIndex")) if "queryRecordIndex" in item
            else (position, position + 1) if kind in _RECORD_PAIR_RESULT_KINDS
            else ()
        )
        identity = (kind, *pair)
        occurrence = 0
        while identity + (occurrence,) in items:
            occurrence += 1
        items[identity + (occurrence,)] = item
    return items


def _comparison_resource_id(item: Mapping[str, object]) -> str:
    resource_id = item.get("resourceId")
    if not isinstance(resource_id, str):
        raise ValidationError("A stored comparison result requires a resourceId.", diagnostic={"code": "INPUT_INVALID"})
    return resource_id


def replace_resource_bytes(descriptor: SessionResourceDescriptor, content: bytes) -> SessionResourceDescriptor:
    """``descriptor`` holding ``content``: base64 data, size and any checksum."""

    replaced: SessionResourceDescriptor = {
        **descriptor, "encoding": "base64", "data": base64.b64encode(content).decode("ascii"), "size": len(content),
    }
    checksum = descriptor.get("checksum")
    if checksum not in (None, ""):
        prefix = "sha256:" if str(checksum).lower().startswith("sha256:") else ""
        replaced["checksum"] = prefix + hashlib.sha256(content).hexdigest()
    return replaced


def _main_display_frame_rows_to_search_frame(
    plan: DiagramRequestPlan,
    request: DiagramRequest,
    drawing: SessionDrawingArtifacts,
) -> DiagramRequest:
    """Convert main Session comparison rows to the search frame once (PD-OI-073).

    Sessions up to version 42 (origin/main 4556e04e) stored nucleotide rows
    after the reverse complement of their endpoint records. The x -> L + 1 - x
    map is its own inverse, so the projection that draws search-frame rows also
    converts these rows; the adapted request and any re-saved Session then use
    the search frame. The decoded resource text of a ``nucleotideBlast`` item
    is converted with its rows.
    """

    if (
        not isinstance(plan, LinearRequestPlan)
        or not isinstance(request, LinearDiagramRequest)
        or drawing.version > MAIN_DISPLAY_FRAME_SESSION_VERSION
        or not request.options.linear_comparisons
    ):
        return request
    originals = request.options.linear_comparisons
    projected = project_search_frame_comparisons(originals, plan.records)
    if all(
        before.matches.equals(after.matches)
        for before, after in zip(originals, projected, strict=True)
    ):
        return request
    converted = tuple(
        after
        if before.search_frame_text is None or before.matches.equals(after.matches)
        else replace(
            after,
            search_frame_text=reverse_endpoint_table_text(
                before.search_frame_text,
                *(
                    _endpoint_frame(plan.records[index])
                    for index in (before.query_record_index, before.subject_record_index)
                ),
            ),
        )
        for before, after in zip(originals, projected, strict=True)
    )
    return replace(request, options=replace(request.options, linear_comparisons=converted))


# The (qualifier, value) columns of the request tables whose `hash` rows a
# Web writer of a Session 44 or older filled with drawn hashes: the color
# table (from the specific color rules) and the label override table (from
# the label text edits of features without a unique qualifier).
_HASH_TABLE_COLUMNS = {"color": (1, 2), "label": (2, 3)}


def _main_drawn_hash_selectors_to_source(
    plan: DiagramRequestPlan,
    request: DiagramRequest,
    drawing: SessionDrawingArtifacts,
) -> tuple[DiagramRequest, int]:
    """Name the ``hash`` rows of the color and label override tables of a
    Session 44 or older by source hash (S6).

    The Web writer of such a Session filled those rows with drawn hashes
    (rendered IDs in the label override table), which it matched against the
    drawn feature, so on a cropped or reverse-complemented record they are not
    the source hashes ``hash=`` names now (OV-401, release D-39). Each value
    that names features drawn with its hash
    (:func:`gbdraw.session_io.source_hash_selector_value`) gets their source
    hash; a table file gets a copy beside it. Returns the request and the
    number of values that name features of several source hashes, which stay
    as saved. A request with no record whose drawn hashes differ is returned
    as is.
    """

    if drawing.version > DRAWN_HASH_RULE_SESSION_MAX_VERSION:
        return request, 0
    drawn = [
        feature
        for ordinal, (record, provenance) in enumerate(
            zip(plan.records, plan.provenance, strict=True), start=1
        )
        for feature in _drawn_features(record, provenance, ordinal)
    ]
    if all(source == drawn_hash for source, drawn_hash, _, _ in drawn):
        return request, 0
    unmapped = 0

    def source_value(value: object) -> str:
        nonlocal unmapped
        resolved, source_count = source_hash_selector_value(value, drawn)
        unmapped += source_count > 1
        return str(value) if resolved is None else resolved

    def table_text(text: str, table: str) -> str:
        qualifier_column, value_column = _HASH_TABLE_COLUMNS[table]
        lines = []
        for line in text.splitlines(keepends=True):
            cells = line.rstrip("\r\n").split("\t")
            ending = line[len(line.rstrip("\r\n")):]
            if (
                not line.lstrip().startswith("#")
                and len(cells) > value_column
                and cells[qualifier_column].strip().lower() == "hash"
                and cells[value_column].strip()
            ):
                cells[value_column] = source_value(cells[value_column].strip())
                line = "\t".join(cells) + ending
            lines.append(line)
        return "".join(lines)

    def table_file(path: str | None, table: str) -> str | None:
        if path is None:
            return None
        source = Path(path)
        text = source.read_text(encoding="utf-8")
        migrated = table_text(text, table)
        if migrated == text:
            return path
        copy_path = source.with_name(f"{source.stem}.source-hash{source.suffix}")
        copy_path.write_text(migrated, encoding="utf-8")
        return str(copy_path)

    def table_frame(frame: DataFrame | None, qualifier: str) -> DataFrame | None:
        if frame is None or qualifier not in frame.columns or "value" not in frame.columns:
            return frame
        rows = frame[qualifier].astype(str).str.strip().str.lower() == "hash"
        if not rows.any():
            return frame
        migrated = frame.copy()
        migrated.loc[rows, "value"] = [
            source_value(value) for value in migrated.loc[rows, "value"]
        ]
        return migrated

    options = request.options
    changes: dict[str, Any] = {}
    colors = options.colors
    if colors is not None:
        color_table = table_frame(colors.color_table, "qualifier_key")
        color_table_file = table_file(colors.color_table_file, "color")
        if color_table is not colors.color_table or color_table_file != colors.color_table_file:
            changes["colors"] = replace(
                colors, color_table=color_table, color_table_file=color_table_file
            )
    label_table = table_frame(options.label_override_table, "qualifier")
    if label_table is not options.label_override_table:
        changes["label_override_table"] = label_table
    label_file = table_file(options.label_override_file, "label")
    if label_file != options.label_override_file:
        changes["label_override_file"] = label_file
    if not changes:
        return request, unmapped
    # Each request type takes its own options type, which mypy cannot pair.
    return cast(DiagramRequest, replace(cast(Any, request), options=replace(options, **changes))), unmapped


def _adapt_session_plan(
    request: DiagramRequest,
    plan: DiagramRequestPlan,
    drawing: SessionDrawingArtifacts,
) -> tuple[AdaptedSessionRequest, ProteinExtractionResult | None, DiagramRequestPlan]:
    """Adapt the request and artifacts of ``plan``, the plan of ``request``;
    also return the protein extraction of the plan when the adaptation ran one,
    so its build reuses it, and the plan it adapted (``plan``, or ``request``
    planned again)."""

    source = _read_session_artifact_source(drawing.fields)
    # Tables and annotation targets are prepared with the plan, so a request
    # whose hash selectors change is planned again. ``plan.request`` is a
    # resolved projection (its regions are applied), so the saved request is
    # the one rewritten and planned.
    rewritten, unmapped_hash_selectors = _main_drawn_hash_selectors_to_source(
        plan, request, drawing
    )
    if rewritten is not request:
        plan = plan_request(rewritten)
    request = promote_legacy_session_similarity_alignment_request(
        plan.request,
        drawing,
        plan=plan if isinstance(plan, LinearRequestPlan) else None,
    )
    request = _main_display_frame_rows_to_search_frame(plan, request, drawing)
    unconverted: tuple[UnmatchedFrameItem, ...] = ()
    if drawing.version <= DRAWN_FRAME_PROTEIN_RESULT_SESSION_VERSION:
        request, unconverted = drawn_frame_protein_results_to_search_frame(plan, request)
    if request is not plan.request:
        plan = _replace_plan_request(plan, request)
    current_raw = source.current_raw_entries
    manifest = source.protein_identity_manifest
    promoted_indexes: frozenset[int] = frozenset()
    id_map: dict[str, str] = {}
    extraction: ProteinExtractionResult | None = None

    if isinstance(plan, LinearRequestPlan):
        assert isinstance(request, LinearDiagramRequest)
        references = _protein_references(_request_protein_artifacts(request))
        reference_ids = references.legacy_ids
        current_protein = tuple(
            entry for entry in current_raw if is_protein_losat_cache_entry(entry)
        )
        attempt_pairs = (
            _legacy_candidate_attempt_pairs(request, source.protein_source_mode)
            if any(candidate.get("state") == "pending" for candidate in source.legacy_candidates)
            else ()
        )
        needs_extraction = bool(attempt_pairs or reference_ids or current_protein)
        promoted: tuple[Mapping[str, Any], ...] = ()
        if needs_extraction:
            if plan.inputs is None:  # pragma: no cover - planner invariant.
                raise ValidationError(
                    "Linear session compatibility plan has no prepared inputs."
                )
            extraction = _extract_linear_request_proteins(
                request,
                plan.records,
                plan.inputs,
            )
            extracted_manifest = extraction.identity_manifest
            assert extracted_manifest is not None
            promoted, promoted_indexes, promoted_id_map = _promote_legacy_candidates(
                source.legacy_candidates,
                extraction,
                attempt_pairs,
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
            manifest = extracted_manifest
        request = _rewrite_linear_request_protein_references(
            request,
            id_map,
            references,
        )

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
    unresolved, evidence = written_legacy_protein_artifacts(
        source.legacy_candidates,
        source.legacy_derived_evidence,
        promoted=promoted_indexes,
        has_current_derived=bool(source.current_derived_entries),
    )
    warnings = (
        (
            f"{len(unresolved)} legacy protein raw cache candidate(s) "
            "could not be promoted.",
        )
        if unresolved
        else ()
    ) + (
        (unmapped_hash_selector_warning(unmapped_hash_selectors, drawing.version),)
        if unmapped_hash_selectors
        else ()
    )
    for warning in warnings[len(warnings) - bool(unmapped_hash_selectors):]:
        logger.warning("WARNING: %s", warning)
    return AdaptedSessionRequest(
        request=request,
        artifacts=artifacts,
        migration_report=SessionMigrationReport(
            protein_raw_candidates=unresolved,
            protein_derived_evidence=evidence,
            protein_id_map=id_map,
            warnings=warnings,
            unconverted_protein_results=tuple(map(describe_unmatched_frame_item, unconverted)),
        ),
    ), extraction, plan


def _require_drawing(drawing: object) -> SessionDrawingArtifacts:
    if not isinstance(drawing, SessionDrawingArtifacts):
        raise ValidationError(
            "Session compatibility input must be a drawing of a loaded SessionDocument."
        )
    return drawing


def adapt_session_request(
    request: DiagramRequest,
    drawing: SessionDrawingArtifacts,
) -> AdaptedSessionRequest:
    """Convert one drawing of a validated Session to the current render contract."""

    return _adapt_session_plan(
        request, plan_request(request), _require_drawing(drawing)
    )[0]


def _adjust_migration_report(
    prepared: PreparedDiagramRequest | PreparedCircularBatchRequest,
    adapted: AdaptedSessionRequest,
) -> SessionMigrationReport:
    """The report of ``adapted``; current derived entries supersede the legacy evidence.

    They are the Session's own, which a re-save writes back, or else the ones
    this render built: only a Python API replay (``render_session_drawings``)
    builds them; a CLI render and a batch build none.
    """

    report = adapted.migration_report
    _candidates, evidence = written_legacy_protein_artifacts(
        report.protein_raw_candidates,
        report.protein_derived_evidence,
        promoted=(),
        has_current_derived=bool(
            adapted.artifacts.losat_derived_cache_entries
            or (
                isinstance(prepared, PreparedDiagramRequest)
                and prepared.losat_derived_cache_entries
            )
        ),
    )
    return replace(report, protein_derived_evidence=evidence)


def build_session_compatible_request_diagram(
    request: DiagramRequest,
    drawing: SessionDrawingArtifacts,
) -> PreparedDiagramRequest | PreparedCircularBatchRequest:
    """Build one drawing of a released session after adapting its artifacts."""

    prepared, _report = _build_session_compatible_plan(
        request,
        plan_request(request),
        _require_drawing(drawing),
    )
    return prepared


def _build_session_compatible_plan(
    request: DiagramRequest,
    plan: DiagramRequestPlan,
    drawing: SessionDrawingArtifacts,
) -> tuple[
    PreparedDiagramRequest | PreparedCircularBatchRequest,
    SessionMigrationReport,
]:
    """Adapt and build ``plan``, the already-resolved plan of ``request``."""

    adapted, extraction, plan = _adapt_session_plan(request, plan, drawing)
    if adapted.request is not plan.request:
        # The replacement keeps the records, their inputs and the record keys
        # that the extraction read.
        plan = _replace_plan_request(plan, adapted.request)
    prepared = _build_request_plan_diagram(
        plan,
        artifacts=adapted.artifacts,
        protein_extraction=extraction,
    )
    return prepared, _adjust_migration_report(prepared, adapted)


def render_session_compatible_request(
    request: DiagramRequest,
    drawing: SessionDrawingArtifacts,
    *,
    include_feature_catalog: bool = False,
) -> RequestRenderResult | CircularBatchRenderResult:
    """Render one drawing of a released session through its sole compatibility adapter.

    Both result shapes carry the drawing's ``migration_report``.
    """

    _require_drawing(drawing)
    plan = plan_request(request)
    batch_outputs_preflighted = isinstance(plan, CircularBatchRequestPlan)
    plan.preflight_outputs()
    prepared, migration_report = _build_session_compatible_plan(
        request,
        plan,
        drawing,
    )
    result = render_prepared_request(
        prepared,
        batch_outputs_preflighted=batch_outputs_preflighted,
        include_feature_catalog=include_feature_catalog,
    )
    return replace(result, migration_report=migration_report)


# main's Web writers of Session 39, and of Session 40 before 3bca0e8d9 (#319),
# stored a reversed record's drawn-frame feature ID in a protein table's
# ``*_feature_svg_id`` and wrote no ``*_view_feature_svg_id``; a later Session
# 40 writes both, so the missing column, not the version, tells them apart
# (OV-1183).
DRAWN_FEATURE_ID_SESSION_MAX_VERSION = 40


def _stores_drawn_feature_ids(frame: DataFrame) -> bool:
    # The projection's own predicate: only rows with source bindings are projected.
    if frame.empty or not _SOURCE_BINDING_COLUMNS <= set(frame.columns):
        return False
    return any(f"{role}_view_feature_svg_id" not in frame.columns for role in ("query", "subject"))


def decode_session_request(
    session_version: int,
    payload: Mapping[str, Any],
    *,
    resource_paths: Mapping[str, str | Path],
    output_directory: str | Path,
) -> DiagramRequest:
    """Decode one stored request of a Session of ``session_version``.

    The canonical payload is migrated first, and a stored font size <= 0 reads
    as Auto (``with_stored_font_sizes_as_auto``); a Linear table of a Session up
    to 40 without view ID columns is marked ``drawn_feature_ids``, so the
    planner reads its stored IDs as views (OV-1183).
    """

    from gbdraw.session_request_codec import decode_canonical_request, with_stored_font_sizes_as_auto

    request = decode_canonical_request(
        with_stored_font_sizes_as_auto(canonical_payload_for_session_decode(session_version, payload)),
        resource_paths=resource_paths,
        output_directory=output_directory,
    )
    if session_version > DRAWN_FEATURE_ID_SESSION_MAX_VERSION or not isinstance(request, LinearDiagramRequest):
        return request
    return _with_drawn_feature_id_marks(request, lambda comparison: _stores_drawn_feature_ids(comparison.matches))


def _with_drawn_feature_id_marks(
    request: LinearDiagramRequest,
    marked: Callable[[LinearComparison], bool],
) -> LinearDiagramRequest:
    """``request`` whose comparisons carry ``drawn_feature_ids`` where ``marked`` says so."""

    stored = request.options.linear_comparisons or ()
    comparisons = tuple(
        comparison if comparison.drawn_feature_ids == marked(comparison)
        else replace(comparison, drawn_feature_ids=marked(comparison))
        for comparison in stored
    )
    if all(before is after for before, after in zip(stored, comparisons, strict=True)):
        return request
    adapted = replace(request, options=replace(request.options, linear_comparisons=comparisons))
    # ``replace`` drops the old-schema alignment state (an ``init=False`` field).
    legacy = request._legacy_similarity_alignment
    return adapted if legacy is None else _with_legacy_similarity_alignment(adapted, legacy)


def resolve_session_write_request(request: DiagramRequest) -> DiagramRequest:
    """Resolve ``request`` for a Session writer.

    A table still read with the drawn-frame IDs of a Session up to 40 moves
    to the search frame first, as an upgrade moves it, so the written table
    holds source and view IDs (OV-1183). A row matching neither frame raises.
    """

    plan = plan_request(request)
    resolved = plan.request
    if not isinstance(resolved, LinearDiagramRequest) or not any(
        comparison.drawn_feature_ids for comparison in resolved.options.linear_comparisons or ()
    ):
        return resolved
    converted, unmatched = drawn_frame_protein_results_to_search_frame(plan, resolved)
    if unmatched:
        raise ValidationError(f"Cannot write {describe_unmatched_frame_item(unmatched[0])}.")
    legacy = resolved._legacy_similarity_alignment
    return converted if legacy is None or not isinstance(converted, LinearDiagramRequest) else (
        _with_legacy_similarity_alignment(converted, legacy))


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
    "DRAWN_FEATURE_ID_SESSION_MAX_VERSION",
    "DRAWN_FRAME_PROTEIN_RESULT_SESSION_VERSION",
    "MAIN_DISPLAY_FRAME_SESSION_VERSION",
    "AdaptedSessionRequest",
    "SearchFrameNucleotideTables",
    "SessionMigrationReport",
    "adapt_session_request",
    "build_session_compatible_request_diagram",
    "canonical_payload_for_session_decode",
    "decode_session_request",
    "drawn_frame_protein_results_to_search_frame",
    "replace_resource_bytes",
    "resolve_session_write_request",
    "search_frame_protein_resources",
    "search_frame_nucleotide_resources",
    "stores_protein_results",
    "project_legacy_similarity_alignment_for_current_write",
    "promote_legacy_session_similarity_alignment_request",
    "render_session_compatible_request",
    "rewrite_protein_artifact_references",
    "written_legacy_protein_artifacts",
]
