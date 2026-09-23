"""Typed Web adapter for the shared Similarity Group alignment resolver."""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from typing import Any

from gbdraw.exceptions import ValidationError
from gbdraw.layout.similarity_alignment import (
    AlignmentAnchorIdentity,
    AlignmentEvidenceEdge,
    AlignmentRecordChoice,
    AlignmentRecordDecision,
    AmbiguousAlignmentRecord,
    SimilarityAlignmentCandidate,
    SimilarityAlignmentPlan,
    resolve_similarity_alignment,
)


def _object(value: object, keys: set[str], path: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping) or set(value) != keys:
        raise ValidationError(f"{path} has invalid fields.")
    return value


def _list(value: object, path: str) -> Sequence[Any]:
    if isinstance(value, (str, bytes)) or not isinstance(value, Sequence):
        raise ValidationError(f"{path} must be an array.")
    return value


def _text(value: object, path: str) -> str:
    if not isinstance(value, str) or not value.strip() or "\0" in value:
        raise ValidationError(f"{path} must be non-empty text without NUL.")
    return value.strip()


def _boolean(value: object, path: str) -> bool:
    if not isinstance(value, bool):
        raise ValidationError(f"{path} must be a boolean.")
    return value


def _integer(
    value: object,
    path: str,
    *,
    minimum: int = 0,
    nullable: bool = False,
) -> int | None:
    if nullable and value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        suffix = " or null" if nullable else ""
        raise ValidationError(
            f"{path} must be an integer greater than or equal to {minimum}{suffix}."
        )
    return value


def _anchor(value: object, path: str) -> AlignmentAnchorIdentity:
    raw = _object(
        value,
        {
            "recordKey",
            "biologicalFeatureId",
            "sourceFeatureIndex",
            "stableFeatureSvgId",
        },
        path,
    )
    stable_id = raw["stableFeatureSvgId"]
    if stable_id is not None:
        stable_id = _text(stable_id, f"{path}.stableFeatureSvgId")
    return AlignmentAnchorIdentity(
        record_key=_text(raw["recordKey"], f"{path}.recordKey"),
        biological_feature_id=_text(
            raw["biologicalFeatureId"], f"{path}.biologicalFeatureId"
        ),
        source_feature_index=_integer(
            raw["sourceFeatureIndex"],
            f"{path}.sourceFeatureIndex",
            nullable=True,
        ),
        stable_feature_svg_id=stable_id,
    )


def _record_fact(value: object, path: str) -> dict[str, Any]:
    raw = _object(
        value,
        {"recordKey", "recordLength", "region", "presentation"},
        path,
    )
    record_length = _integer(
        raw["recordLength"], f"{path}.recordLength", minimum=1, nullable=True
    )
    presentation = _object(
        raw["presentation"], {"reverseComplement"}, f"{path}.presentation"
    )
    region = raw["region"]
    normalized_region = None
    if region is not None:
        region = _object(
            region,
            {"start", "end", "reverseComplement"},
            f"{path}.region",
        )
        start = _integer(region["start"], f"{path}.region.start", minimum=1)
        end = _integer(region["end"], f"{path}.region.end", minimum=1)
        assert start is not None and end is not None
        if start > end:
            raise ValidationError(f"{path}.region start must not exceed end.")
        normalized_region = {
            "start": start,
            "end": end,
            "reverseComplement": _boolean(
                region["reverseComplement"], f"{path}.region.reverseComplement"
            ),
        }
    return {
        "recordKey": _text(raw["recordKey"], f"{path}.recordKey"),
        "recordLength": record_length,
        "region": normalized_region,
        "presentation": {
            "reverseComplement": _boolean(
                presentation["reverseComplement"],
                f"{path}.presentation.reverseComplement",
            )
        },
    }


def _member_candidate(
    value: object,
    path: str,
    *,
    group_id: str,
    record_facts: Mapping[str, Mapping[str, Any]],
) -> SimilarityAlignmentCandidate:
    raw = _object(
        value,
        {
            "groupId",
            "anchor",
            "sourceStart",
            "sourceEnd",
            "sourceStrand",
            "identityIsUnique",
            "hidden",
            "representative",
            "role",
        },
        path,
    )
    if _text(raw["groupId"], f"{path}.groupId") != group_id:
        raise ValidationError(f"{path} belongs to another Similarity Group.")
    anchor = _anchor(raw["anchor"], f"{path}.anchor")
    record = record_facts.get(anchor.record_key)
    if record is None:
        raise ValidationError(f"{path} belongs to an unknown displayed record.")
    start = _integer(raw["sourceStart"], f"{path}.sourceStart")
    end = _integer(raw["sourceEnd"], f"{path}.sourceEnd")
    assert start is not None and end is not None
    if end < start:
        raise ValidationError(f"{path}.sourceEnd must not precede sourceStart.")
    source_strand = raw["sourceStrand"]
    if source_strand is not None and (
        isinstance(source_strand, bool) or source_strand not in (-1, 1)
    ):
        raise ValidationError(f"{path}.sourceStrand must be -1, 1, or null.")
    role = raw["role"]
    if not isinstance(role, str) or "\0" in role:
        raise ValidationError(f"{path}.role must be text without NUL.")

    region = record["region"]
    effective_reverse = (
        region["reverseComplement"]
        if region is not None
        else record["presentation"]["reverseComplement"]
    )
    source_center = (start + end) / 2
    if region is None:
        center_mappable = True
        record_length = record["recordLength"]
        display_center = (
            float(record_length) - source_center
            if effective_reverse and record_length is not None
            else source_center
        )
    else:
        region_start = region["start"] - 1
        region_end = region["end"]
        center_mappable = region_start <= source_center <= region_end
        display_center = None
        if center_mappable:
            display_center = (
                region_end - source_center
                if effective_reverse
                else source_center - region_start
            )
    displayed_strand = (
        None
        if source_strand is None
        else source_strand * (-1 if effective_reverse else 1)
    )
    return SimilarityAlignmentCandidate(
        group_id=group_id,
        anchor=anchor,
        displayed_strand=displayed_strand,
        center_mappable=center_mappable,
        display_center=display_center,
        identity_is_unique=_boolean(
            raw["identityIsUnique"], f"{path}.identityIsUnique"
        ),
        hidden=_boolean(raw["hidden"], f"{path}.hidden"),
        effective_reverse_complement=effective_reverse,
        representative=_boolean(raw["representative"], f"{path}.representative"),
        role=role.strip(),
    )


def _serialize_anchor(anchor: AlignmentAnchorIdentity) -> dict[str, object]:
    return {
        "recordKey": anchor.record_key,
        "biologicalFeatureId": anchor.biological_feature_id,
        "sourceFeatureIndex": anchor.source_feature_index,
        "stableFeatureSvgId": anchor.stable_feature_svg_id,
    }


def _serialize_decision(decision: AlignmentRecordDecision) -> dict[str, object]:
    return {
        "kind": "decision",
        "recordKey": decision.record_key,
        "status": decision.status.value,
        "rationale": decision.rationale.value,
        "anchor": _serialize_anchor(decision.anchor) if decision.anchor else None,
        "effectiveReverseComplement": decision.effective_reverse_complement,
    }


def _serialize_ambiguity(
    ambiguity: AmbiguousAlignmentRecord,
) -> dict[str, object]:
    return {
        "kind": "ambiguous",
        "recordKey": ambiguity.record_key,
        "candidates": [
            {
                "anchor": _serialize_anchor(candidate.anchor),
                "displayedStrand": candidate.displayed_strand,
                "hidden": candidate.hidden,
                "representative": candidate.representative,
                "role": candidate.role,
            }
            for candidate in ambiguity.candidates
        ],
        "directRbhCandidates": [
            _serialize_anchor(anchor) for anchor in ambiguity.direct_rbh_candidates
        ],
    }


def _serialize_plan(plan: SimilarityAlignmentPlan) -> dict[str, object]:
    return {
        "schema": plan.schema,
        "mode": plan.mode.value,
        "groupId": plan.group_id,
        "reference": _serialize_anchor(plan.reference),
        "records": [
            {
                "recordKey": decision.record_key,
                "status": decision.status.value,
                "rationale": decision.rationale.value,
                "anchor": _serialize_anchor(decision.anchor)
                if decision.anchor
                else None,
                "effectiveReverseComplement": decision.effective_reverse_complement,
            }
            for decision in plan.records
        ],
    }


def resolve_similarity_alignment_payload(payload: object) -> dict[str, object]:
    """Validate Web facts, call the shared resolver, and serialize its result."""

    raw = _object(
        payload,
        {
            "schema",
            "mode",
            "groupId",
            "records",
            "reference",
            "members",
            "directEdges",
            "choices",
        },
        "request",
    )
    if raw["schema"] != 1:
        raise ValidationError("request.schema must be 1.")
    group_id = _text(raw["groupId"], "request.groupId")
    records = [
        _record_fact(value, f"request.records[{index}]")
        for index, value in enumerate(_list(raw["records"], "request.records"))
    ]
    record_keys = [record["recordKey"] for record in records]
    if not record_keys or len(record_keys) != len(set(record_keys)):
        raise ValidationError("request.records must have unique record keys.")
    record_facts = {record["recordKey"]: record for record in records}
    reference = _anchor(raw["reference"], "request.reference")
    candidates = [
        _member_candidate(
            value,
            f"request.members[{index}]",
            group_id=group_id,
            record_facts=record_facts,
        )
        for index, value in enumerate(_list(raw["members"], "request.members"))
    ]
    candidate_keys = {candidate.anchor.canonical_key for candidate in candidates}
    edges = []
    for index, value in enumerate(_list(raw["directEdges"], "request.directEdges")):
        path = f"request.directEdges[{index}]"
        edge = _object(
            value, {"groupId", "query", "subject", "edgeKind"}, path
        )
        edge_group_id = _text(edge["groupId"], f"{path}.groupId")
        query = _anchor(edge["query"], f"{path}.query")
        subject = _anchor(edge["subject"], f"{path}.subject")
        if edge_group_id != group_id or any(
            anchor.canonical_key not in candidate_keys for anchor in (query, subject)
        ):
            raise ValidationError(f"{path} is not direct evidence for current members.")
        edges.append(
            AlignmentEvidenceEdge(
                group_id=edge_group_id,
                query=query,
                subject=subject,
                edge_kind=_text(edge["edgeKind"], f"{path}.edgeKind"),
            )
        )
    choices = []
    for index, value in enumerate(_list(raw["choices"], "request.choices")):
        path = f"request.choices[{index}]"
        choice = _object(value, {"recordKey", "kind", "anchor"}, path)
        choices.append(
            AlignmentRecordChoice(
                record_key=_text(choice["recordKey"], f"{path}.recordKey"),
                kind=_text(choice["kind"], f"{path}.kind"),
                anchor=(
                    _anchor(choice["anchor"], f"{path}.anchor")
                    if choice["anchor"] is not None
                    else None
                ),
            )
        )
    resolution = resolve_similarity_alignment(
        record_keys=record_keys,
        group_id=group_id,
        reference=reference,
        candidates=candidates,
        edges=edges,
        choices=choices,
        mode=_text(raw["mode"], "request.mode"),
    )
    plan = resolution.plan
    return {
        "schema": 1,
        "status": "ambiguous" if resolution.ambiguities else "resolved",
        "mode": resolution.mode.value,
        "groupId": resolution.group_id,
        "reference": _serialize_anchor(resolution.reference),
        "records": [
            _serialize_ambiguity(record)
            if isinstance(record, AmbiguousAlignmentRecord)
            else _serialize_decision(record)
            for record in resolution.records
        ],
        "plan": _serialize_plan(plan) if plan is not None else None,
    }


def resolve_similarity_alignment_json(payload_json: object) -> str:
    """JSON boundary used by the existing lazy diagram Worker."""

    return json.dumps(resolve_similarity_alignment_payload(json.loads(str(payload_json))))
