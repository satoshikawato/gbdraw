from __future__ import annotations

from typing import Any, Callable, Mapping, Sequence


_COORD_BASE_KEY = "gbdraw_coord_base"
_COORD_STEP_KEY = "gbdraw_coord_step"


def _get_value(source: Any, key: str, default: Any = None) -> Any:
    if isinstance(source, Mapping):
        return source.get(key, default)
    return getattr(source, key, default)


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


def _read_record_coord_map(record: Any) -> tuple[int, int]:
    annotations = getattr(record, "annotations", None) or {}
    try:
        base = int(annotations.get(_COORD_BASE_KEY, 1))
    except (TypeError, ValueError):
        base = 1
    try:
        step = int(annotations.get(_COORD_STEP_KEY, 1))
    except (TypeError, ValueError):
        step = 1
    if step == 0:
        step = 1
    return base, (1 if step > 0 else -1)


def _record_coord_maps(records: Sequence[Any] | None) -> dict[int, tuple[int, int]]:
    if records is None:
        return {}
    return {index: _read_record_coord_map(record) for index, record in enumerate(records)}


def _absolute_display_interval(
    start: int,
    end: int,
    coord_base: int,
    coord_step: int,
) -> tuple[int, int]:
    if end <= start:
        coord = coord_base + (coord_step * start)
        return coord - 1, coord
    first_coord = coord_base + (coord_step * start)
    last_coord = coord_base + (coord_step * (end - 1))
    min_coord = min(first_coord, last_coord)
    max_coord = max(first_coord, last_coord)
    return min_coord - 1, max_coord


def _serialize_orthogroup_name_candidate(candidate: Any) -> dict[str, object]:
    return {
        "text": _text(_get_value(candidate, "text")),
        "source": _text(_get_value(candidate, "source")),
        "memberCount": _int_value(_get_value(candidate, "member_count")),
        "recordCoverageCount": _int_value(_get_value(candidate, "record_coverage_count")),
        "representativeCount": _int_value(_get_value(candidate, "representative_count")),
        "score": _float_value(_get_value(candidate, "score")),
    }


def _serialize_ortholog_edge(edge: Any) -> dict[str, object]:
    return {
        "orthogroupId": _text(_get_value(edge, "orthogroup_id")),
        "sourceRbhOrthogroupId": _get_value(edge, "source_rbh_orthogroup_id"),
        "targetRbhOrthogroupId": _get_value(edge, "target_rbh_orthogroup_id"),
        "queryProteinId": _text(_get_value(edge, "query_protein_id")),
        "subjectProteinId": _text(_get_value(edge, "subject_protein_id")),
        "queryRecordIndex": _int_value(_get_value(edge, "query_record_index")),
        "subjectRecordIndex": _int_value(_get_value(edge, "subject_record_index")),
        "edgeKind": _text(_get_value(edge, "edge_kind")),
        "renderRole": _text(_get_value(edge, "render_role")),
        "pathId": _get_value(edge, "path_id"),
        "identity": _float_value(_get_value(edge, "identity")),
        "evalue": _float_value(_get_value(edge, "evalue")),
        "bitscore": _float_value(_get_value(edge, "bitscore")),
        "alignmentLength": _int_value(_get_value(edge, "alignment_length")),
    }


def _serialize_ortholog_path(path: Any) -> dict[str, object]:
    return {
        "orthogroupId": _text(_get_value(path, "orthogroup_id")),
        "pathId": _text(_get_value(path, "path_id")),
        "proteinIds": [_text(item) for item in (_get_value(path, "protein_ids", ()) or ())],
        "edgeIds": [_text(item) for item in (_get_value(path, "edge_ids", ()) or ())],
        "sharedProteinIds": [
            _text(item) for item in (_get_value(path, "shared_protein_ids", ()) or ())
        ],
    }


FeatureIdMapper = Callable[[int, str], str]


def _serialize_member(
    member: Any,
    feature_id_mapper: FeatureIdMapper | None = None,
    record_coord_maps: Mapping[int, tuple[int, int]] | None = None,
) -> dict[str, object]:
    record_index = _int_value(_get_value(member, "record_index"))
    stable_feature_svg_id = _text(_get_value(member, "feature_svg_id"))
    feature_svg_id = (
        feature_id_mapper(record_index, stable_feature_svg_id)
        if feature_id_mapper is not None and stable_feature_svg_id
        else stable_feature_svg_id
    )
    start = _int_value(_get_value(member, "start"))
    end = _int_value(_get_value(member, "end"))
    coord_map = (record_coord_maps or {}).get(record_index)
    if coord_map is not None:
        start, end = _absolute_display_interval(start, end, coord_map[0], coord_map[1])
    return {
        "orthogroupId": _text(_get_value(member, "orthogroup_id")),
        "proteinId": _text(_get_value(member, "protein_id")),
        "sourceProteinId": _get_value(member, "source_protein_id"),
        "recordIndex": record_index,
        "recordId": _text(_get_value(member, "record_id")),
        "featureIndex": _int_value(_get_value(member, "feature_index")),
        "label": _text(_get_value(member, "label")),
        "featureSvgId": feature_svg_id,
        "stableFeatureSvgId": stable_feature_svg_id,
        "stable_feature_svg_id": stable_feature_svg_id,
        "start": start,
        "end": end,
        "strand": _strand_display(_get_value(member, "strand")),
        "representative": bool(_get_value(member, "representative")),
        "role": _text(_get_value(member, "role")),
        "confidence": _text(_get_value(member, "confidence")),
        "assignmentReason": _text(_get_value(member, "assignment_reason")),
        "supportingEdges": [
            _text(edge_id) for edge_id in (_get_value(member, "supporting_edges", ()) or ())
        ],
        "bestCoreSupport": _float_value(_get_value(member, "best_core_support")),
        "secondBestCoreSupport": _float_value(_get_value(member, "second_best_core_support")),
        "gene": _get_value(member, "gene"),
        "product": _get_value(member, "product"),
        "note": _get_value(member, "note"),
        "locusTag": _get_value(member, "locus_tag"),
        "geneId": _get_value(member, "gene_id"),
        "oldLocusTag": _get_value(member, "old_locus_tag"),
    }


def serialize_orthogroups_payload(
    orthogroups: Any,
    *,
    feature_id_mapper: FeatureIdMapper | None = None,
    records: Sequence[Any] | None = None,
    record_coord_maps: Mapping[int, tuple[int, int]] | None = None,
) -> list[dict[str, object]]:
    """Serialize OrthogroupResult to the web interactive SVG payload shape."""

    if orthogroups is None:
        return []
    resolved_coord_maps = record_coord_maps or _record_coord_maps(records)
    groups = _get_value(orthogroups, "orthogroups", {}) or {}
    names_by_id = _get_value(orthogroups, "names_by_orthogroup_id", {}) or {}
    descriptions_by_id = _get_value(orthogroups, "descriptions_by_orthogroup_id", {}) or {}
    candidates_by_id = _get_value(orthogroups, "name_candidates_by_orthogroup_id", {}) or {}
    confidence_by_id = _get_value(orthogroups, "confidence_by_orthogroup_id", {}) or {}
    rbh_orthogroups = _get_value(orthogroups, "rbh_orthogroups", {}) or {}
    ortholog_edges_by_id = _get_value(orthogroups, "ortholog_edges_by_orthogroup_id", {}) or {}
    ortholog_paths_by_id = _get_value(orthogroups, "ortholog_paths_by_orthogroup_id", {}) or {}
    related_edges_by_id = _get_value(orthogroups, "related_edges_by_orthogroup_id", {}) or {}
    scope_by_id = _get_value(orthogroups, "scope_by_orthogroup_id", {}) or {}
    source_record_index_by_id = (
        _get_value(orthogroups, "source_record_index_by_orthogroup_id", {}) or {}
    )

    payload: list[dict[str, object]] = []
    for orthogroup_id, members in groups.items():
        orthogroup_id = _text(orthogroup_id)
        members = list(members or [])
        member_ids = {_text(_get_value(member, "protein_id")) for member in members}
        rbh_ids = [
            _text(rbh_id)
            for rbh_id, protein_ids in rbh_orthogroups.items()
            if member_ids.intersection({_text(protein_id) for protein_id in protein_ids})
        ]
        record_coverage_count = len(
            {
                int(_get_value(member, "record_index"))
                for member in members
                if _int_or_none(_get_value(member, "record_index")) is not None
            }
        )
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
            if not feature_svg_id or record_index is None:
                continue
            entry = _build_feature_index_entry(group, member)
            feature_index[(record_index, feature_svg_id)] = entry
            feature_index_by_svg_id.setdefault(feature_svg_id, entry)
            if stable_feature_svg_id:
                feature_index.setdefault((record_index, stable_feature_svg_id), entry)
                feature_index_by_svg_id.setdefault(stable_feature_svg_id, entry)

    enriched: list[dict[str, object]] = []
    for feature in features:
        next_feature = dict(feature)
        svg_id = _text(next_feature.get("svg_id")).strip()
        stable_svg_id = _text(
            next_feature.get("stable_svg_id") or next_feature.get("stableFeatureSvgId")
        ).strip()
        record_index = _int_or_none(
            next_feature.get("record_idx", next_feature.get("recordIndex"))
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
