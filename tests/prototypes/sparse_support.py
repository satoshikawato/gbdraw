"""S04 test-only scans frozen from d7695a62 (tree equal to a66237f4)."""
from __future__ import annotations
from typing import Mapping, Sequence
from gbdraw.analysis.protein_colinearity import (
    CdsProtein, OrthogroupResult, _BestCoreEvidence, _CoreSupportCandidate,
    _LocalThreshold, OrthologEdgeKind, _normalized_score_from_row,
    _anchor_core_hit_rank, _row_supports_membership, _row_min_coverage,
    _row_domain_only, _protein_sort_key, _RECORD_LOCAL_CORE_COMPETITION_RATIO,
    _MIN_MEMBERSHIP_MIN_COVERAGE,
)


def _record_local_component_has_competing_core_support(
    member_ids: Sequence[str],
    group_member_ids: Mapping[str, set[str]],
    best_by_direction: Mapping[tuple[str, str], object],
    thresholds: Mapping[str, _LocalThreshold],
    protein_map: Mapping[str, CdsProtein],
    local_support_by_member: Mapping[str, float],
) -> bool:
    if not group_member_ids:
        return False
    credible_core_ids: set[str] = set()
    for member_id in member_ids:
        member_local_support = float(local_support_by_member.get(member_id, 0.0))
        if member_local_support <= 0.0:
            return True
        candidates = [
            candidate
            for group_id, existing_member_ids in group_member_ids.items()
            if (
                candidate := _build_core_support_candidate(
                    member_id,
                    group_id,
                    sorted(existing_member_ids, key=lambda item: _protein_sort_key(protein_map[item])),
                    best_by_direction,
                    thresholds,
                    protein_map,
                )
            )
            is not None
        ]
        if not candidates:
            continue
        strongest = max(
            candidates,
            key=lambda candidate: (
                max(
                    float(candidate.same_record_score),
                    float(candidate.cross_record_score),
                    float(candidate.diagnostic_score),
                ),
                float(candidate.diagnostic_score),
                candidate.group_id,
            ),
        )
        strongest_score = max(
            float(strongest.same_record_score),
            float(strongest.cross_record_score),
            float(strongest.diagnostic_score),
        )
        competition_floor = member_local_support / _RECORD_LOCAL_CORE_COMPETITION_RATIO
        if strongest.domain_only and strongest_score >= competition_floor:
            return True
        for candidate in candidates:
            core_support = max(
                float(candidate.cross_record_score),
                float(candidate.same_record_score),
            )
            if core_support <= 0.0:
                continue
            if candidate.high_confidence_pass:
                return True
            if candidate.domain_only:
                if core_support >= competition_floor:
                    return True
                continue
            if candidate.low_confidence_pass and core_support >= competition_floor:
                credible_core_ids.add(candidate.group_id)
    return len(credible_core_ids) >= 2


def _best_evidence_between_protein_and_members(
    protein_id: str,
    member_ids: Sequence[str],
    best_by_direction: Mapping[tuple[str, str], object],
    protein_map: Mapping[str, CdsProtein],
    *,
    same_record: bool,
) -> _BestCoreEvidence:
    protein = protein_map[protein_id]
    best_support: tuple[float, object | None, str, str] = (0.0, None, "", "")
    best_diagnostic: tuple[float, object | None, str, str] = (0.0, None, "", "")
    member_set = set(member_ids)
    for member_id in member_set:
        if member_id == protein_id or member_id not in protein_map:
            continue
        member = protein_map[member_id]
        is_same_record = int(protein.record_index) == int(member.record_index)
        if is_same_record != bool(same_record):
            continue
        for query_id, subject_id in ((protein_id, member_id), (member_id, protein_id)):
            row = best_by_direction.get((query_id, subject_id))
            if row is None:
                continue
            score = _normalized_score_from_row(row)
            if score <= 0.0:
                continue
            current_diagnostic_row = best_diagnostic[1]
            if (
                current_diagnostic_row is None
                or score > best_diagnostic[0]
                or (
                    score == best_diagnostic[0]
                    and _anchor_core_hit_rank(row, protein_map) < _anchor_core_hit_rank(current_diagnostic_row, protein_map)
                )
            ):
                best_diagnostic = (float(score), row, query_id, subject_id)
            if not _row_supports_membership(row):
                continue
            current_support_row = best_support[1]
            if (
                current_support_row is None
                or score > best_support[0]
                or (
                    score == best_support[0]
                    and _anchor_core_hit_rank(row, protein_map) < _anchor_core_hit_rank(current_support_row, protein_map)
                )
            ):
                best_support = (float(score), row, query_id, subject_id)
    return _BestCoreEvidence(
        support_score=float(best_support[0]),
        support_row=best_support[1],
        support_query_id=best_support[2],
        support_subject_id=best_support[3],
        diagnostic_score=float(best_diagnostic[0]),
        diagnostic_row=best_diagnostic[1],
        diagnostic_query_id=best_diagnostic[2],
        diagnostic_subject_id=best_diagnostic[3],
    )


def _build_core_support_candidate(
    protein_id: str,
    group_id: str,
    member_ids: Sequence[str],
    best_by_direction: Mapping[tuple[str, str], object],
    thresholds: Mapping[str, _LocalThreshold],
    protein_map: Mapping[str, CdsProtein],
) -> _CoreSupportCandidate | None:
    same_evidence = _best_evidence_between_protein_and_members(
        protein_id,
        member_ids,
        best_by_direction,
        protein_map,
        same_record=True,
    )
    cross_evidence = _best_evidence_between_protein_and_members(
        protein_id,
        member_ids,
        best_by_direction,
        protein_map,
        same_record=False,
    )
    if same_evidence.diagnostic_row is None and cross_evidence.diagnostic_row is None:
        return None

    same_score = float(same_evidence.support_score)
    cross_score = float(cross_evidence.support_score)
    support = float(cross_score) + 0.5 * float(same_score)
    diagnostic_score = max(
        float(same_evidence.diagnostic_score),
        float(cross_evidence.diagnostic_score),
    )

    candidate_threshold = thresholds.get(protein_id)
    same_member_id = ""
    if same_evidence.support_row is not None:
        same_member_id = (
            same_evidence.support_subject_id
            if same_evidence.support_query_id == protein_id
            else same_evidence.support_query_id
        )
    same_member_threshold = thresholds.get(same_member_id)
    same_pass = (
        same_evidence.support_row is not None
        and (
            (
                candidate_threshold is not None
                and same_score >= float(candidate_threshold.score)
            )
            or (
                same_member_threshold is not None
                and same_score >= float(same_member_threshold.score)
            )
        )
    )
    cross_pass = (
        cross_evidence.support_row is not None
        and candidate_threshold is not None
        and cross_score >= float(candidate_threshold.score)
    )

    if same_pass:
        evidence_row = same_evidence.support_row
        evidence_query_id = same_evidence.support_query_id
        evidence_subject_id = same_evidence.support_subject_id
        relation_kind: OrthologEdgeKind = "same_record_inparalog"
    elif cross_pass:
        evidence_row = cross_evidence.support_row
        evidence_query_id = cross_evidence.support_query_id
        evidence_subject_id = cross_evidence.support_subject_id
        relation_kind = "coortholog"
    elif same_evidence.support_row is not None and (
        cross_evidence.support_row is None or same_score >= cross_score
    ):
        evidence_row = same_evidence.support_row
        evidence_query_id = same_evidence.support_query_id
        evidence_subject_id = same_evidence.support_subject_id
        relation_kind = "same_record_inparalog"
    elif cross_evidence.support_row is not None:
        evidence_row = cross_evidence.support_row
        evidence_query_id = cross_evidence.support_query_id
        evidence_subject_id = cross_evidence.support_subject_id
        relation_kind = "coortholog"
    elif same_evidence.diagnostic_row is not None and (
        cross_evidence.diagnostic_row is None
        or same_evidence.diagnostic_score >= cross_evidence.diagnostic_score
    ):
        evidence_row = same_evidence.diagnostic_row
        evidence_query_id = same_evidence.diagnostic_query_id
        evidence_subject_id = same_evidence.diagnostic_subject_id
        relation_kind = "same_record_inparalog"
    else:
        evidence_row = cross_evidence.diagnostic_row
        evidence_query_id = cross_evidence.diagnostic_query_id
        evidence_subject_id = cross_evidence.diagnostic_subject_id
        relation_kind = "coortholog"
    if evidence_row is None:
        return None

    min_coverage = _row_min_coverage(evidence_row)
    domain_only = _row_domain_only(evidence_row)
    has_support_row = same_evidence.support_row is not None or cross_evidence.support_row is not None
    high_confidence_pass = bool(has_support_row and (same_pass or cross_pass)) and min_coverage >= _MIN_MEMBERSHIP_MIN_COVERAGE and not domain_only
    low_confidence_pass = bool(has_support_row) and min_coverage >= _MIN_MEMBERSHIP_MIN_COVERAGE and not domain_only
    if same_pass:
        reason = "same-record support passed local anchor threshold"
    elif cross_pass:
        reason = "cross-record support passed local fallback threshold"
    elif domain_only:
        reason = "support is domain-only"
    elif min_coverage < _MIN_MEMBERSHIP_MIN_COVERAGE:
        reason = "support is below membership coverage"
    else:
        reason = "support is separated but below high-confidence threshold"
    return _CoreSupportCandidate(
        group_id=group_id,
        support=support,
        diagnostic_score=diagnostic_score,
        same_record_score=float(same_score),
        cross_record_score=float(cross_score),
        evidence_row=evidence_row,
        evidence_query_id=evidence_query_id,
        evidence_subject_id=evidence_subject_id,
        relation_kind=relation_kind,
        min_coverage=float(min_coverage),
        domain_only=bool(domain_only),
        high_confidence_pass=high_confidence_pass,
        low_confidence_pass=low_confidence_pass,
        reason=reason,
    )


def _edge_metadata_for_protein_pair(
    orthogroups: OrthogroupResult | None,
    orthogroup_id: str,
    query_id: str,
    subject_id: str,
) -> dict[str, object]:
    metadata = {
        "rbh_orthogroup_id": "",
        "ortholog_path_id": "",
        "edge_kind": "",
        "render_role": "",
    }
    if orthogroups is None or not orthogroup_id:
        return metadata
    candidate_edges = [
        *orthogroups.ortholog_edges_by_orthogroup_id.get(orthogroup_id, ()),
        *orthogroups.related_edges_by_orthogroup_id.get(orthogroup_id, ()),
    ]
    for edge in candidate_edges:
        if (
            edge.query_protein_id == query_id
            and edge.subject_protein_id == subject_id
        ) or (
            edge.query_protein_id == subject_id
            and edge.subject_protein_id == query_id
        ):
            source_group = str(edge.source_rbh_orthogroup_id or "")
            target_group = str(edge.target_rbh_orthogroup_id or "")
            if source_group and target_group and source_group != target_group:
                rbh_group = f"{source_group};{target_group}"
            else:
                rbh_group = source_group or target_group
            metadata.update(
                {
                    "rbh_orthogroup_id": rbh_group,
                    "ortholog_path_id": str(edge.path_id or ""),
                    "edge_kind": edge.edge_kind,
                    "render_role": edge.render_role,
                }
            )
            return metadata
    return metadata
