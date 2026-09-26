"""Typed Web adapter for the shared Similarity Group alignment resolver."""

from __future__ import annotations

import json
import hashlib
from dataclasses import replace
from pathlib import Path
from collections.abc import Mapping, Sequence
from typing import Any

from gbdraw.exceptions import ValidationError
from gbdraw.api.record_planning import (
    ResolvedRecordCollection, project_similarity_alignment_anchor_fact,
)
from gbdraw.api.request_render import plan_linear_request
from gbdraw.api.requests import LinearDiagramRequest
from gbdraw.session_request_codec import decode_canonical_request
from gbdraw.layout.similarity_alignment import (
    AlignmentAnchorIdentity,
    AlignmentEvidenceEdge,
    AlignmentRecordChoice,
    AlignmentRecordDecision,
    AlignmentReviewCandidate,
    AlignmentReviewRow,
    AmbiguousAlignmentRecord,
    SimilarityAlignmentCandidate,
    SimilarityAlignmentPlan,
    resolve_similarity_alignment,
)


def _object(
    value: object,
    keys: set[str],
    path: str,
    *,
    optional: set[str] = frozenset(),
) -> Mapping[str, Any]:
    if (
        not isinstance(value, Mapping)
        or not keys.issubset(value)
        or not set(value).issubset(keys | optional)
    ):
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
    presentation_reverse = _boolean(
        presentation["reverseComplement"], f"{path}.presentation.reverseComplement"
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
        if record_length is not None and end > record_length:
            raise ValidationError(f"{path}.region exceeds recordLength.")
        normalized_region = {
            "start": start,
            "end": end,
            "reverseComplement": _boolean(
                region["reverseComplement"], f"{path}.region.reverseComplement"
            ),
        }
    if normalized_region is not None and normalized_region["reverseComplement"] and presentation_reverse:
        raise ValidationError(
            f"{path} cannot reverse complement both region and presentation."
        )
    return {
        "recordKey": _text(raw["recordKey"], f"{path}.recordKey"),
        "recordLength": record_length,
        "region": normalized_region,
        "presentation": {"reverseComplement": presentation_reverse},
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
        optional={"displayName"},
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
    record_length = record["recordLength"]
    presentation_reverse = record["presentation"]["reverseComplement"]
    region_reverse = region["reverseComplement"] if region is not None else False
    effective_reverse = presentation_reverse ^ region_reverse
    source_center = (start + end) / 2
    if presentation_reverse and record_length is None:
        presented_center = None
    else:
        presented_center = (
            float(record_length) - source_center
            if presentation_reverse else source_center
        )
    if region is None:
        display_center = presented_center
    else:
        region_start = region["start"] - 1
        region_end = region["end"]
        display_center = None
        if presented_center is not None and region_start <= presented_center <= region_end:
            display_center = (
                region_end - presented_center
                if region_reverse else presented_center - region_start
            )
    displayed_strand = (
        None if source_strand is None
        else source_strand * (-1 if effective_reverse else 1)
    )
    display_name = raw.get("displayName", anchor.biological_feature_id)
    return SimilarityAlignmentCandidate(
        group_id=group_id,
        anchor=anchor,
        displayed_strand=displayed_strand,
        center_mappable=display_center is not None,
        display_center=display_center,
        identity_is_unique=_boolean(
            raw["identityIsUnique"], f"{path}.identityIsUnique"
        ),
        hidden=_boolean(raw["hidden"], f"{path}.hidden"),
        representative=_boolean(raw["representative"], f"{path}.representative"),
        role=role.strip(),
        source_start=start,
        source_end=end,
        display_name=_text(display_name, f"{path}.displayName"),
    )


def _serialize_anchor(anchor: AlignmentAnchorIdentity) -> dict[str, object]:
    return {
        "recordKey": anchor.record_key,
        "biologicalFeatureId": anchor.biological_feature_id,
        "sourceFeatureIndex": anchor.source_feature_index,
        "stableFeatureSvgId": anchor.stable_feature_svg_id,
    }


def _serialize_review_candidate(row: AlignmentReviewCandidate) -> dict[str, object]:
    candidate = row.candidate
    return {
        "anchor": _serialize_anchor(candidate.anchor),
        "displayName": candidate.display_name,
        "sourceStart": candidate.source_start,
        "sourceEnd": candidate.source_end,
        "displayCenter": candidate.display_center,
        "displayedStrand": candidate.displayed_strand,
        "hidden": candidate.hidden,
        "representative": candidate.representative,
        "role": candidate.role,
        "usable": row.usable,
        "directEvidence": list(row.direct_evidence),
        "strandRelation": row.strand_relation.value,
    }


def _serialize_decision(
    decision: AlignmentRecordDecision,
    row: AlignmentReviewRow,
) -> dict[str, object]:
    return {
        "kind": "decision",
        "recordKey": decision.record_key,
        "status": decision.status.value,
        "rationale": decision.rationale.value,
        "reviewReason": (
            decision.review_reason.value if decision.review_reason else None
        ),
        "anchor": _serialize_anchor(decision.anchor) if decision.anchor else None,
        "candidates": [_serialize_review_candidate(item) for item in row.candidates],
    }


def _serialize_ambiguity(
    ambiguity: AmbiguousAlignmentRecord,
    row: AlignmentReviewRow,
) -> dict[str, object]:
    return {
        "kind": "ambiguous",
        "recordKey": ambiguity.record_key,
        "candidates": [_serialize_review_candidate(item) for item in row.candidates],
        "directRbhCandidates": [
            _serialize_anchor(anchor) for anchor in ambiguity.direct_rbh_candidates
        ],
        "recommendedAnchor": _serialize_anchor(ambiguity.recommended_anchor),
        "recommendationReason": ambiguity.recommendation_reason.value,
    }


def _serialize_plan(plan: SimilarityAlignmentPlan) -> dict[str, object]:
    return {
        "schema": plan.schema,
        "groupId": plan.group_id,
        "reference": _serialize_anchor(plan.reference),
        "records": [
            {
                "recordKey": decision.record_key,
                "status": decision.status.value,
                "rationale": decision.rationale.value,
                "anchor": _serialize_anchor(decision.anchor)
                if decision.anchor else None,
            }
            for decision in plan.records
        ],
    }



def _projection_context(value: object, resource_paths: Mapping[str, str], output_directory: str):
    raw = _object(value, {"canonicalRequest", "orientations"}, "projection")
    request = decode_canonical_request(
        raw["canonicalRequest"], resource_paths=resource_paths, output_directory=output_directory,
    )
    if not isinstance(request, LinearDiagramRequest):
        raise ValidationError("Alignment projection requires a Linear request.")
    keys = [record.record_key for record in request.records]
    if len(set(keys)) != len(keys) or None in keys:
        raise ValidationError("Alignment projection requires keyed records.")
    # The caller materializes current placements through the composition bridge.
    # The plan is deliberately removed for raw, untranslated placement facts.
    request = replace(request, similarity_alignment=None)
    before = plan_linear_request(request)
    before_collection = ResolvedRecordCollection(before.records, before.provenance)
    orientations = raw["orientations"]
    if orientations is None:
        orientations = {key: item.source_step == 1
                        for key, item in zip(keys, before.transforms, strict=True)}
    else:
        orientations = _object(orientations, set(keys), "projection.orientations")
        orientations = {key: _boolean(value, f"projection.orientations.{key}")
                        for key, value in orientations.items()}
    records = []
    for record in request.records:
        reverse = orientations[record.record_key]
        if record.region is not None:
            records.append(replace(record, region=replace(record.region, reverse_complement=reverse),
                                   presentation=replace(record.presentation, reverse_complement=False)))
        else:
            records.append(replace(record, presentation=replace(record.presentation, reverse_complement=reverse)))
    after = plan_linear_request(replace(request, records=tuple(records)))
    # Counterfactual source transforms are always available for local direction edits.
    flipped_records = []
    for record, transform in zip(request.records, before.transforms, strict=True):
        reverse = transform.source_step == 1
        if record.region is not None:
            flipped_records.append(replace(record, region=replace(record.region, reverse_complement=reverse),
                                           presentation=replace(record.presentation, reverse_complement=False)))
        else:
            flipped_records.append(replace(record, presentation=replace(record.presentation, reverse_complement=reverse)))
    flipped = plan_linear_request(replace(request, records=tuple(flipped_records)))
    return raw, before, before_collection, after, flipped


def _serialize_projection(context, candidates, raw_request):
    raw, before, collection, after, flipped = context
    before_placements = getattr(before.build().drawing, "_gbdraw_alignment_placements")
    after_placements = getattr(after.build().drawing, "_gbdraw_alignment_placements")
    flipped_collection = ResolvedRecordCollection(flipped.records, flipped.provenance)
    fingerprints = {}
    for provenance in collection.provenance:
        for path in provenance.source_paths:
            fingerprints[path] = hashlib.sha256(Path(path).read_bytes()).hexdigest()
    binding_records = []
    records = []
    base = {item.record_key: item for item in (before.layout.record_translations if before.layout else ())}
    for index, provenance in enumerate(collection.provenance):
        key = provenance.record_key
        before_reverse = before.transforms[index].source_step == -1
        after_reverse = after.transforms[index].source_step == -1
        translation = base.get(key)
        record_raw = raw["canonicalRequest"]["records"][index]
        binding_records.append({
            "record": {"recordKey": key, "source": record_raw["source"],
                       "selector": record_raw["selector"], "region": record_raw["region"],
                       "display": {"isCircular": provenance.display.is_circular,
                                   "startCoordinate": provenance.display.start_coordinate},
                       "cardinality": provenance.cardinality.value},
            "reverseComplement": before_reverse,
            "translation": {"x": 0.0 if translation is None else float(translation.x),
                            "y": 0.0 if translation is None else float(translation.y)},
            "sourceRecordIndex": provenance.source_record_index,
            "sourceRecordId": provenance.source_record_id,
            "sourceFingerprints": [fingerprints[path] for path in provenance.source_paths],
        })
        variants = []
        for reverse in (False, True):
            variant_collection = collection if reverse == before_reverse else flipped_collection
            placement = after_placements[index] if reverse == after_reverse else before_placements[index]
            anchors = []
            for candidate in candidates:
                if candidate.anchor.record_key != key:
                    continue
                fact = project_similarity_alignment_anchor_fact(variant_collection, candidate.anchor)
                anchors.append({
                    "anchor": _serialize_anchor(candidate.anchor),
                    "displayedStrand": fact.displayed_strand,
                    "displayCenter": fact.display_center,
                    "centerX": None if fact.display_center is None else placement.x_for_position(fact.display_center),
                })
            variants.append({"reverseComplement": reverse, "axisY": placement.axis_y, "anchors": anchors})
        records.append({
            "recordKey": key, "beforeReverseComplement": before_reverse,
            "base": {"x": 0.0 if translation is None else float(translation.x),
                     "y": 0.0 if translation is None else float(translation.y)},
            "beforeAxisY": before_placements[index].axis_y,
            "beforeAnchors": [
                {"anchor": _serialize_anchor(candidate.anchor),
                 "centerX": before_placements[index].x_for_position(candidate.display_center)}
                for candidate in candidates if candidate.anchor.record_key == key and candidate.center_mappable
            ],
            "variants": variants,
        })
    binding = {"records": binding_records, "groupId": raw_request["groupId"],
               "reference": raw_request["reference"],
               "anchors": [_serialize_anchor(candidate.anchor) for candidate in candidates]}
    return {"binding": hashlib.sha256(json.dumps(binding, sort_keys=True, separators=(",", ":")).encode()).hexdigest(),
            "geometryOrientations": {item.record_key: transform.source_step == -1
                                     for item, transform in zip(after.provenance, after.transforms, strict=True)},
            "records": records}

def resolve_similarity_alignment_payload(
    payload: object, *, projection: object = None,
    resource_paths: Mapping[str, str] | None = None, output_directory: str = "/tmp",
) -> dict[str, object]:
    """Validate Web facts, call the shared resolver, and serialize its result."""

    raw = _object(
        payload,
        {
            "schema",
            "groupId",
            "records",
            "reference",
            "members",
            "directEdges",
            "choices",
        },
        "request",
    )
    if type(raw["schema"]) is not int or raw["schema"] != 2:
        raise ValidationError("request.schema must be 2.")
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
    context = None
    if projection is not None:
        context = _projection_context(projection, resource_paths or {}, output_directory)
        collection = context[2]
        if tuple(item.record_key for item in collection.provenance) != tuple(record_keys):
            raise ValidationError("Alignment projection record binding changed.")
        verified = []
        for candidate, member in zip(candidates, raw["members"], strict=True):
            fact = project_similarity_alignment_anchor_fact(collection, candidate.anchor)
            if (member["sourceStart"], member["sourceEnd"], member["sourceStrand"]) != (
                fact.source_start, fact.source_end, fact.source_strand,
            ):
                raise ValidationError("Alignment projection source anchor facts changed.")
            verified.append(replace(candidate, displayed_strand=fact.displayed_strand,
                                    display_center=fact.display_center, center_mappable=fact.display_center is not None))
        candidates = verified
    candidate_keys = {candidate.anchor.canonical_key for candidate in candidates}
    edges = []
    for index, value in enumerate(_list(raw["directEdges"], "request.directEdges")):
        path = f"request.directEdges[{index}]"
        edge = _object(value, {"groupId", "query", "subject", "edgeKind"}, path)
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
                    if choice["anchor"] is not None else None
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
    )
    plan = resolution.plan
    return {
        "schema": 2,
        "status": "ambiguous" if resolution.ambiguities else "resolved",
        "groupId": resolution.group_id,
        "reference": _serialize_anchor(resolution.reference),
        "referenceDisplayedStrand": resolution.reference_candidate.displayed_strand,
        "referenceDisplayCenter": resolution.reference_candidate.display_center,
        "records": [
            _serialize_ambiguity(record, row)
            if isinstance(record, AmbiguousAlignmentRecord)
            else _serialize_decision(record, row)
            for record, row in zip(
                resolution.records, resolution.review_rows, strict=True
            )
        ],
        "plan": _serialize_plan(plan) if plan is not None else None,
        "projection": _serialize_projection(context, candidates, raw) if context is not None else None,
    }


def resolve_similarity_alignment_json(
    payload_json: object, projection_json: object = "null", resource_paths_json: object = "{}",
    output_directory: str = "/tmp",
) -> str:
    """JSON boundary used by the existing lazy diagram Worker."""

    return json.dumps(resolve_similarity_alignment_payload(
        json.loads(str(payload_json)), projection=json.loads(str(projection_json)),
        resource_paths=json.loads(str(resource_paths_json)), output_directory=output_directory,
    ))
