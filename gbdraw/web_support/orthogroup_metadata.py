from __future__ import annotations

from typing import Any, Callable, Mapping, Sequence

from gbdraw.core.record_metadata import (
    _absolute_display_interval,
    _read_coord_map as _read_record_coord_map,
)


def _text(value: Any) -> str:
    return str(value or "")


def _int_or_none(value: Any) -> int | None:
    if value is None:
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _int_value(value: Any, default: int = 0) -> int:
    parsed = _int_or_none(value)
    return default if parsed is None else parsed


def _float_value(value: Any, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _strand_display(strand: Any) -> str:
    if strand == 1 or str(strand).strip() == "1":
        return "+"
    if strand == -1 or str(strand).strip() == "-1":
        return "-"
    return _text(strand)


def _record_coord_maps(records: Sequence[Any] | None) -> dict[int, tuple[int, int]]:
    if records is None:
        return {}
    return {index: _read_record_coord_map(record) for index, record in enumerate(records)}


def _serialize_orthogroup_name_candidate(candidate: Any) -> dict[str, object]:
    getter = candidate.get if isinstance(candidate, Mapping) else lambda key, default=None: getattr(candidate, key, default)
    return {
        "text": _text(getter("text", "")),
        "source": _text(getter("source", "")),
        "memberCount": _int_value(getter("member_count", 0)),
        "recordCoverageCount": _int_value(getter("record_coverage_count", 0)),
        "representativeCount": _int_value(getter("representative_count", 0)),
        "score": _float_value(getter("score", 0.0)),
    }


def _serialize_ortholog_edge(edge: Any) -> dict[str, object]:
    return {
        "orthogroupId": _text(getattr(edge, "orthogroup_id", "")),
        "sourceRbhOrthogroupId": getattr(edge, "source_rbh_orthogroup_id", None),
        "targetRbhOrthogroupId": getattr(edge, "target_rbh_orthogroup_id", None),
        "queryProteinId": _text(getattr(edge, "query_protein_id", "")),
        "subjectProteinId": _text(getattr(edge, "subject_protein_id", "")),
        "queryRecordIndex": _int_value(getattr(edge, "query_record_index", 0)),
        "subjectRecordIndex": _int_value(getattr(edge, "subject_record_index", 0)),
        "edgeKind": _text(getattr(edge, "edge_kind", "")),
        "renderRole": _text(getattr(edge, "render_role", "")),
        "pathId": getattr(edge, "path_id", None),
        "identity": _float_value(getattr(edge, "identity", 0.0)),
        "evalue": _float_value(getattr(edge, "evalue", 0.0)),
        "bitscore": _float_value(getattr(edge, "bitscore", 0.0)),
        "alignmentLength": _int_value(getattr(edge, "alignment_length", 0)),
    }


def _serialize_ortholog_path(path: Any) -> dict[str, object]:
    return {
        "orthogroupId": _text(getattr(path, "orthogroup_id", "")),
        "pathId": _text(getattr(path, "path_id", "")),
        "proteinIds": [_text(item) for item in (getattr(path, "protein_ids", ()) or ())],
        "edgeIds": [_text(item) for item in (getattr(path, "edge_ids", ()) or ())],
        "sharedProteinIds": [
            _text(item) for item in (getattr(path, "shared_protein_ids", ()) or ())
        ],
    }


FeatureIdMapper = Callable[[int, str], str]


def _serialize_member(
    member: Any,
    feature_id_mapper: FeatureIdMapper | None = None,
    record_coord_maps: Mapping[int, tuple[int, int]] | None = None,
    include_stable_feature_ids: bool = True,
) -> dict[str, object]:
    record_index = int(member.record_index)
    stable_feature_svg_id = str(member.feature_svg_id or "")
    rendered_feature_svg_id = (
        feature_id_mapper(record_index, stable_feature_svg_id)
        if feature_id_mapper is not None and stable_feature_svg_id
        else ""
    )
    start = int(member.start)
    end = int(member.end)
    coord_map = (record_coord_maps or {}).get(record_index)
    if coord_map is not None:
        start, end = _absolute_display_interval(start, end, coord_map[0], coord_map[1])
    payload = {
        "orthogroupId": member.orthogroup_id,
        "proteinId": member.protein_id,
        "sourceProteinId": member.source_protein_id,
        "recordIndex": record_index,
        "recordId": member.record_id,
        "featureIndex": member.feature_index,
        "label": member.label,
        "featureSvgId": stable_feature_svg_id,
        "start": start,
        "end": end,
        "strand": _strand_display(member.strand),
        "representative": member.representative,
        "role": str(getattr(member, "role", "") or ""),
        "confidence": str(getattr(member, "confidence", "") or ""),
        "assignmentReason": str(getattr(member, "assignment_reason", "") or ""),
        "supportingEdges": [
            str(edge_id) for edge_id in (getattr(member, "supporting_edges", ()) or ())
        ],
        "bestCoreSupport": float(getattr(member, "best_core_support", 0.0) or 0.0),
        "secondBestCoreSupport": float(getattr(member, "second_best_core_support", 0.0) or 0.0),
        "gene": getattr(member, "gene", None),
        "product": getattr(member, "product", None),
        "note": getattr(member, "note", None),
        "locusTag": getattr(member, "locus_tag", None),
        "geneId": getattr(member, "gene_id", None),
        "oldLocusTag": getattr(member, "old_locus_tag", None),
    }
    if include_stable_feature_ids:
        payload["stableFeatureSvgId"] = stable_feature_svg_id
        payload["stable_feature_svg_id"] = stable_feature_svg_id
    if rendered_feature_svg_id:
        payload["renderedFeatureSvgId"] = rendered_feature_svg_id
        payload["rendered_feature_svg_id"] = rendered_feature_svg_id
    return payload


def serialize_orthogroups_payload(
    orthogroups: Any,
    *,
    feature_id_mapper: FeatureIdMapper | None = None,
    records: Sequence[Any] | None = None,
    record_coord_maps: Mapping[int, tuple[int, int]] | None = None,
    include_stable_feature_ids: bool = True,
) -> list[dict[str, object]]:
    """Serialize OrthogroupResult to the web interactive SVG payload shape."""

    if orthogroups is None:
        return []
    resolved_coord_maps = record_coord_maps or _record_coord_maps(records)
    groups = getattr(orthogroups, "orthogroups", {}) or {}
    names_by_id = getattr(orthogroups, "names_by_orthogroup_id", {}) or {}
    descriptions_by_id = getattr(orthogroups, "descriptions_by_orthogroup_id", {}) or {}
    candidates_by_id = getattr(orthogroups, "name_candidates_by_orthogroup_id", {}) or {}
    confidence_by_id = getattr(orthogroups, "confidence_by_orthogroup_id", {}) or {}
    rbh_orthogroups = getattr(orthogroups, "rbh_orthogroups", {}) or {}
    ortholog_edges_by_id = getattr(orthogroups, "ortholog_edges_by_orthogroup_id", {}) or {}
    ortholog_paths_by_id = getattr(orthogroups, "ortholog_paths_by_orthogroup_id", {}) or {}
    related_edges_by_id = getattr(orthogroups, "related_edges_by_orthogroup_id", {}) or {}
    scope_by_id = getattr(orthogroups, "scope_by_orthogroup_id", {}) or {}
    source_record_index_by_id = getattr(orthogroups, "source_record_index_by_orthogroup_id", {}) or {}

    payload: list[dict[str, object]] = []
    for orthogroup_id, members in groups.items():
        orthogroup_id = _text(orthogroup_id)
        members = list(members or [])
        member_ids = {_text(getattr(member, "protein_id", "")) for member in members}
        rbh_ids = [
            _text(rbh_id)
            for rbh_id, protein_ids in rbh_orthogroups.items()
            if member_ids.intersection({_text(protein_id) for protein_id in protein_ids})
        ]
        record_coverage_count = len({int(member.record_index) for member in members})
        source_record_index = (
            _int_or_none(source_record_index_by_id[orthogroup_id])
            if orthogroup_id in source_record_index_by_id
            else None
        )
        payload.append(
            {
                "id": orthogroup_id,
                "name": _text(names_by_id.get(orthogroup_id)),
                "description": _text(descriptions_by_id.get(orthogroup_id)),
                "nameConfidence": _text(confidence_by_id.get(orthogroup_id, "none")) or "none",
                "scope": _text(scope_by_id.get(orthogroup_id, "cross_record")) or "cross_record",
                "source_record_index": source_record_index,
                "nameCandidates": [
                    _serialize_orthogroup_name_candidate(candidate)
                    for candidate in (candidates_by_id.get(orthogroup_id, []) or [])
                ],
                "member_count": len(members),
                "record_coverage_count": record_coverage_count,
                "rbhOrthogroupIds": rbh_ids,
                "orthologEdges": [
                    _serialize_ortholog_edge(edge)
                    for edge in (ortholog_edges_by_id.get(orthogroup_id, []) or [])
                ],
                "orthologPaths": [
                    _serialize_ortholog_path(path)
                    for path in (ortholog_paths_by_id.get(orthogroup_id, []) or [])
                ],
                "relatedEdges": [
                    _serialize_ortholog_edge(edge)
                    for edge in (related_edges_by_id.get(orthogroup_id, []) or [])
                ],
                "members": [
                    _serialize_member(
                        member,
                        feature_id_mapper=feature_id_mapper,
                        record_coord_maps=resolved_coord_maps,
                        include_stable_feature_ids=include_stable_feature_ids,
                    )
                    for member in members
                ],
            }
        )
    return payload


def _build_feature_index_entry(
    group: Mapping[str, object],
    member: Mapping[str, object],
) -> dict[str, object]:
    members = group.get("members")
    member_count = _int_value(group.get("member_count"), len(members) if isinstance(members, list) else 0)
    record_coverage = _int_value(group.get("record_coverage_count"))
    return {
        "orthogroup_id": _text(group.get("id")),
        "orthogroup_member_count": member_count,
        "orthogroup_record_coverage": record_coverage,
        "protein_id": _text(member.get("proteinId")),
        "source_protein_id": _text(member.get("sourceProteinId")),
        "orthogroup_representative": bool(member.get("representative")),
        "orthogroup_scope": _text(group.get("scope")),
        "orthogroup_source_record_index": group.get("source_record_index"),
        "orthogroup_member": dict(member),
    }


def enrich_features_with_orthogroups(
    features: Sequence[Mapping[str, object]],
    orthogroups: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    """Return feature payloads with GUI-equivalent orthogroup metadata attached."""

    feature_index: dict[tuple[int, str], dict[str, object]] = {}
    feature_index_by_svg_id: dict[str, dict[str, object]] = {}
    global_feature_owners: dict[str, tuple[int, object]] = {}
    ambiguous_svg_ids: set[str] = set()
    for group in orthogroups:
        if not isinstance(group, Mapping):
            continue
        for member in group.get("members", []) or []:
            if not isinstance(member, Mapping):
                continue
            feature_svg_id = _text(member.get("featureSvgId")).strip()
            stable_feature_svg_id = _text(
                member.get("stableFeatureSvgId") or member.get("stable_feature_svg_id")
            ).strip()
            record_index = _int_or_none(member.get("recordIndex"))
            member_feature_ids = {
                feature_id
                for feature_id in (feature_svg_id, stable_feature_svg_id)
                if feature_id
            }
            if not member_feature_ids or record_index is None:
                continue
            entry = _build_feature_index_entry(group, member)
            owner = (record_index, member.get("featureIndex"))
            for member_feature_id in member_feature_ids:
                feature_index.setdefault((record_index, member_feature_id), entry)
                if member_feature_id in ambiguous_svg_ids:
                    continue
                existing_owner = global_feature_owners.get(member_feature_id)
                if existing_owner is not None and existing_owner != owner:
                    feature_index_by_svg_id.pop(member_feature_id, None)
                    ambiguous_svg_ids.add(member_feature_id)
                    continue
                global_feature_owners[member_feature_id] = owner
                feature_index_by_svg_id[member_feature_id] = entry

    enriched: list[dict[str, object]] = []
    for feature in features:
        next_feature = dict(feature)
        svg_id = _text(next_feature.get("svg_id")).strip()
        stable_svg_id = _text(
            next_feature.get("stable_feature_id")
            or next_feature.get("stable_svg_id")
            or next_feature.get("stableFeatureSvgId")
        ).strip()
        record_index = next(
            (
                parsed
                for value in (
                    next_feature.get("record_idx"),
                    next_feature.get("record_index"),
                    next_feature.get("recordIndex"),
                )
                if (parsed := _int_or_none(value)) is not None
            ),
            None,
        )
        entry = None
        if svg_id and record_index is not None:
            entry = feature_index.get((record_index, svg_id))
        if entry is None and stable_svg_id and record_index is not None:
            entry = feature_index.get((record_index, stable_svg_id))
        if entry is None and svg_id:
            entry = feature_index_by_svg_id.get(svg_id)
        if entry is None and stable_svg_id:
            entry = feature_index_by_svg_id.get(stable_svg_id)
        if entry is not None:
            next_feature.update(entry)
        enriched.append(next_feature)
    return enriched


__all__ = [
    "enrich_features_with_orthogroups",
    "serialize_orthogroups_payload",
]
