"""Typed Web adapter for the shared Similarity Group alignment resolver."""

from __future__ import annotations

import json
import hashlib
from dataclasses import dataclass, replace
from pathlib import Path
from collections.abc import Mapping, Sequence, Set
from typing import Any, TypedDict, cast

from gbdraw.exceptions import ValidationError
from gbdraw.api.record_planning import (
    ResolvedRecordCollection, _source_details, alignment_record_keys,
    project_similarity_alignment_anchor_fact, similarity_alignment_evidence_edges,
)
from gbdraw.api.request_render import LinearRequestPlan, plan_linear_request, resolve_request_records
from gbdraw.api.requests import LinearDiagramRequest
from gbdraw.session_request_codec import StoredRecordWindow, decode_canonical_request
from gbdraw.layout.similarity_alignment import (
    SIMILARITY_ALIGNMENT_PLAN_SCHEMA,
    AlignmentAnchorIdentity,
    AlignmentChoiceKind,
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
    optional: Set[str] = frozenset(),
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


class AlignmentHelperPresentation(TypedDict):
    """The presentation fact of a helper request record."""

    reverseComplement: bool


class AlignmentHelperRecord(TypedDict):
    """One displayed record of a helper request, as the Web's ``buildHelperRequest`` sends it.

    ``windows`` and ``presentation`` are the record as request schema 10
    stores them (``StoredRecordWindow``, empty for a whole record); the page
    reads them back to check the reply's projection. Python reads only the
    key: it places each anchor on the committed request of the projection,
    whose records own their windows, orientation and source length (D5-13).
    """

    recordKey: str
    windows: list[StoredRecordWindow]
    presentation: AlignmentHelperPresentation


class _AlignmentHelperCommittedRequest(TypedDict):
    canonicalRequest: dict[str, Any]
    orientations: dict[str, bool] | None


class AlignmentHelperProjection(_AlignmentHelperCommittedRequest, total=False):
    """The committed request a helper call places every candidate on (D5-13), as the Web sends it.

    ``orientations`` are the directions the projection geometry is drawn
    after, ``None`` for every record flipped. ``geometry`` is ``False`` for a
    call that needs only the decisions (the check before a Generate): the
    reply's ``projection`` is then ``None`` and no layout is planned or built.
    """

    geometry: bool


def _record_key(value: object, path: str) -> str:
    """The key of a helper request record."""

    raw = cast(AlignmentHelperRecord, _object(value, {"recordKey", "windows", "presentation"}, path))
    return _text(raw["recordKey"], f"{path}.recordKey")


def _member_candidate(
    value: object,
    path: str,
    *,
    group_id: str,
    collection: ResolvedRecordCollection,
    record_keys: Set[str],
) -> SimilarityAlignmentCandidate:
    """A group member, placed where the committed request draws it (``project_similarity_alignment_anchor_fact``).

    The Web names the member by its identity only: its source coordinates and strand are what the committed request's
    record holds, never a copy the page keeps (an older gbdraw saved a member's strand in the drawn direction).
    """

    raw = _object(
        value,
        {
            "groupId",
            "anchor",
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
    if anchor.record_key not in record_keys:
        raise ValidationError(f"{path} belongs to an unknown displayed record.")
    role = raw["role"]
    if not isinstance(role, str) or "\0" in role:
        raise ValidationError(f"{path}.role must be text without NUL.")

    fact = project_similarity_alignment_anchor_fact(collection, anchor)
    display_name = raw.get("displayName", anchor.biological_feature_id)
    return SimilarityAlignmentCandidate(
        group_id=group_id,
        anchor=anchor,
        displayed_strand=fact.displayed_strand,
        center_mappable=fact.display_center is not None,
        display_center=fact.display_center,
        identity_is_unique=_boolean(
            raw["identityIsUnique"], f"{path}.identityIsUnique"
        ),
        hidden=_boolean(raw["hidden"], f"{path}.hidden"),
        representative=_boolean(raw["representative"], f"{path}.representative"),
        role=role.strip(),
        source_start=fact.source_start,
        source_end=fact.source_end,
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



@dataclass(frozen=True)
class _CommittedPlacement:
    """The committed request a call names (without its plan) and the records it draws, on which every candidate is
    placed, by the renderer's projection, whether or not the call asks for the projection geometry."""

    raw: AlignmentHelperProjection
    request: LinearDiagramRequest
    collection: ResolvedRecordCollection


@dataclass(frozen=True)
class _ProjectionGeometry:
    """The layouts a review draws its previews from: the committed placements and those after ``orientations``."""

    before: LinearRequestPlan
    before_collection: ResolvedRecordCollection
    after: LinearRequestPlan
    flipped_collection: ResolvedRecordCollection


def _committed_placement(value: object, resource_paths: Mapping[str, str], output_directory: str) -> _CommittedPlacement:
    raw = cast(AlignmentHelperProjection, _object(
        value, {"canonicalRequest", "orientations"}, "projection", optional={"geometry"},
    ))
    if "geometry" in raw:
        _boolean(raw["geometry"], "projection.geometry")
    request = decode_canonical_request(
        raw["canonicalRequest"], resource_paths=resource_paths, output_directory=output_directory,
    )
    if not isinstance(request, LinearDiagramRequest):
        raise ValidationError("Alignment projection requires a Linear request.")
    keys = [record.record_key for record in request.records if record.record_key is not None]
    if len(set(keys)) != len(request.records):
        raise ValidationError("Alignment projection requires keyed records.")
    # The caller materializes current placements through the composition bridge.
    # The plan is deliberately removed for raw, untranslated placement facts.
    request = replace(request, similarity_alignment=None)
    # A canonical request carries no planner-only transform (auto_orient, origin_feature: the codec encodes their
    # results), so its resolved records are the records the plan draws.
    collection = resolve_request_records(request)
    if alignment_record_keys(collection) != tuple(keys):
        raise ValidationError("Alignment projection record binding changed.", diagnostic={"code": "INPUT_INVALID"})
    return _CommittedPlacement(raw, request, collection)


def _projection_geometry(committed: _CommittedPlacement) -> _ProjectionGeometry:
    request = committed.request
    keys = [cast(str, record.record_key) for record in request.records]
    before = plan_linear_request(request)
    before_collection = ResolvedRecordCollection(before.records, before.provenance)
    if alignment_record_keys(before_collection) != tuple(keys):
        raise ValidationError("Alignment projection record binding changed.", diagnostic={"code": "INPUT_INVALID"})
    # Counterfactual source transforms are always available for local direction edits.
    flipped_orientations = {key: not _drawn_reverse(before, key) for key in keys}
    raw_orientations = committed.raw["orientations"]
    if raw_orientations is None:
        orientations = flipped_orientations
    else:
        checked = _object(raw_orientations, set(keys), "projection.orientations")
        orientations = {key: _boolean(value, f"projection.orientations.{key}")
                        for key, value in checked.items()}
    after = plan_linear_request(_oriented_request(request, orientations))
    # The flipped records serve only anchor facts, so they are resolved, not
    # planned; without orientations (the stale-target repair) they are the
    # requested ones (OV-1585).
    flipped_collection = (
        ResolvedRecordCollection(after.records, after.provenance)
        if orientations == flipped_orientations
        else resolve_request_records(_oriented_request(request, flipped_orientations))
    )
    return _ProjectionGeometry(before, before_collection, after, flipped_collection)


def _oriented_request(request: LinearDiagramRequest, reverse_by_key: Mapping[str, bool]) -> LinearDiagramRequest:
    """The request with each record drawn in the given absolute direction."""

    records = []
    for record in request.records:
        records.append(record.oriented(reverse_by_key[cast(str, record.record_key)]))
    return replace(request, records=tuple(records))


def _record_units(plan: LinearRequestPlan, record_key: str) -> list[int]:
    """The drawn units of the request record ``record_key``: one, or one per window (D-50)."""

    return [index for index, item in enumerate(plan.provenance) if item.record_identity_key == record_key]


def _drawn_reverse(plan: LinearRequestPlan, record_key: str) -> bool:
    """Whether the request record is drawn reverse-complemented: every unit of it is, as the Web reads a record's
    direction from its windows (``canonicalRecordReverseComplement``)."""

    return all(plan.transforms[index].source_step == -1 for index in _record_units(plan, record_key))


def _serialize_projection(
    committed: _CommittedPlacement,
    geometry: _ProjectionGeometry,
    candidates: Sequence[SimilarityAlignmentCandidate],
    raw_request: Mapping[str, Any],
) -> dict[str, object]:
    raw, request = committed.raw, committed.request
    before, collection, after = geometry.before, geometry.before_collection, geometry.after
    flipped_collection = geometry.flipped_collection
    before_placements = getattr(before.build().drawing, "_gbdraw_alignment_placements")
    after_placements = getattr(after.build().drawing, "_gbdraw_alignment_placements")
    # The content of each source this call staged, hashed once per file. Not
    # ``provenance.source_paths``: records the Worker's cache supplied name the
    # files of the Generate that parsed them.
    fingerprints: dict[int, list[str]] = {}
    hashes: dict[str, str] = {}
    for provenance in collection.provenance:
        _kind, paths = _source_details(request.records[provenance.input_index].source)
        for path in paths:
            if path not in hashes:
                hashes[path] = hashlib.sha256(Path(path).read_bytes()).hexdigest()
        fingerprints[provenance.input_index] = [hashes[path] for path in paths]
    binding_records = []
    records = []
    base = {item.record_key: item for item in (before.layout.record_translations if before.layout else ())}
    # One entry per request record; the windows of a record share its row and move with it, and each anchor is drawn
    # by the window that holds it (D-50).
    # Every record of a committed placement has a key (``_committed_placement``).
    for index, key in enumerate(cast(str, item.record_key) for item in request.records):
        first = _record_units(before, key)[0]
        provenance = collection.provenance[first]
        before_reverse = _drawn_reverse(before, key)
        after_reverse = _drawn_reverse(after, key)
        translation = base.get(key)
        record_raw = raw["canonicalRequest"]["records"][index]
        binding_records.append({
            "record": {"recordKey": key, "source": record_raw["source"],
                       "selector": record_raw["selector"],
                       # Request schema 10 binds a record's windows; an older request its region, as it bound before.
                       **({"windows": record_raw["windows"]} if "windows" in record_raw
                          else {"region": record_raw["region"]}),
                       "display": {"isCircular": provenance.display.is_circular,
                                   "startCoordinate": provenance.display.start_coordinate},
                       "cardinality": provenance.cardinality.value},
            "reverseComplement": before_reverse,
            "translation": {"x": 0.0 if translation is None else float(translation.x),
                            "y": 0.0 if translation is None else float(translation.y)},
            "sourceRecordIndex": provenance.source_record_index,
            "sourceRecordId": provenance.source_record_id,
            "sourceFingerprints": fingerprints[provenance.input_index],
        })
        variants = []
        for reverse in (False, True):
            variant_collection = collection if reverse == before_reverse else flipped_collection
            placements = after_placements if reverse == after_reverse else before_placements
            anchors = []
            for candidate in candidates:
                if candidate.anchor.record_key != key:
                    continue
                fact = project_similarity_alignment_anchor_fact(variant_collection, candidate.anchor)
                anchors.append({
                    "anchor": _serialize_anchor(candidate.anchor),
                    "displayedStrand": fact.displayed_strand,
                    "displayCenter": fact.display_center,
                    "centerX": None if fact.display_center is None
                    else placements[fact.unit_index].x_for_position(fact.display_center),
                })
            variants.append({"reverseComplement": reverse, "axisY": placements[first].axis_y, "anchors": anchors})
        records.append({
            "recordKey": key, "beforeReverseComplement": before_reverse,
            "base": {"x": 0.0 if translation is None else float(translation.x),
                     "y": 0.0 if translation is None else float(translation.y)},
            "beforeAxisY": before_placements[first].axis_y,
            "beforeAnchors": [
                {"anchor": _serialize_anchor(candidate.anchor),
                 "centerX": before_placements[
                     project_similarity_alignment_anchor_fact(collection, candidate.anchor).unit_index
                 ].x_for_position(candidate.display_center)}
                for candidate in candidates if candidate.anchor.record_key == key and candidate.center_mappable
            ],
            "variants": variants,
        })
    binding = {"records": binding_records, "groupId": raw_request["groupId"],
               "reference": raw_request["reference"],
               "anchors": [_serialize_anchor(candidate.anchor) for candidate in candidates]}
    return {"binding": hashlib.sha256(json.dumps(binding, sort_keys=True, separators=(",", ":")).encode()).hexdigest(),
            "geometryOrientations": {item.record_key: _drawn_reverse(after, cast(str, item.record_key))
                                     for item in request.records},
            "records": records}

def _resource_evidence_edges(
    committed: _CommittedPlacement, group_id: str, candidates: Sequence[SimilarityAlignmentCandidate],
) -> list[AlignmentEvidenceEdge]:
    """Bind resource edges to the current members; unbound endpoints are omitted."""

    collection, orthogroups = committed.collection, committed.request.options.orthogroups
    by_key = {candidate.anchor.canonical_key: candidate.anchor for candidate in candidates}

    def current(anchor: AlignmentAnchorIdentity) -> AlignmentAnchorIdentity | None:
        match = by_key.get(anchor.canonical_key)
        if match is None or any(
            left is not None and right is not None and left != right
            for left, right in (
                (anchor.source_feature_index, match.source_feature_index),
                (anchor.stable_feature_svg_id, match.stable_feature_svg_id),
            )
        ):
            return None
        return match

    edges = []
    for edge in similarity_alignment_evidence_edges(collection, orthogroups, group_id):
        query, subject = current(edge.query), current(edge.subject)
        if query is not None and subject is not None:
            edges.append(replace(edge, query=query, subject=subject))
    return edges


def resolve_similarity_alignment_payload(
    payload: object, *, projection: object,
    resource_paths: Mapping[str, str] | None = None, output_directory: str = "/tmp",
) -> dict[str, object]:
    """Validate Web facts, call the shared resolver, and serialize its result.

    ``projection`` names the committed request (``canonicalRequest``) on which every candidate is placed; a payload
    without one is malformed (D5-13). With ``geometry`` false the reply carries the decisions only.
    """

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
    if type(raw["schema"]) is not int or raw["schema"] != SIMILARITY_ALIGNMENT_PLAN_SCHEMA:
        raise ValidationError(f"request.schema must be {SIMILARITY_ALIGNMENT_PLAN_SCHEMA}.")
    group_id = _text(raw["groupId"], "request.groupId")
    committed = _committed_placement(projection, resource_paths or {}, output_directory)
    collection = committed.collection
    record_keys = [
        _record_key(value, f"request.records[{index}]")
        for index, value in enumerate(_list(raw["records"], "request.records"))
    ]
    if not record_keys or len(record_keys) != len(set(record_keys)):
        raise ValidationError("request.records must have unique record keys.")
    if alignment_record_keys(collection) != tuple(record_keys):
        raise ValidationError("Alignment projection record binding changed.", diagnostic={"code": "INPUT_INVALID"})
    reference = _anchor(raw["reference"], "request.reference")
    displayed = set(record_keys)
    candidates = [
        _member_candidate(
            value,
            f"request.members[{index}]",
            group_id=group_id,
            collection=collection,
            record_keys=displayed,
        )
        for index, value in enumerate(_list(raw["members"], "request.members"))
    ]
    # The committed request's orthogroup resource is the only evidence
    # authority; the UI catalog carries no graph.
    if _list(raw["directEdges"], "request.directEdges"):
        raise ValidationError(
            "request.directEdges must be empty; direct evidence comes from "
            "the committed orthogroup resource."
        )
    edges = _resource_evidence_edges(committed, group_id, candidates)
    choices = []
    for index, value in enumerate(_list(raw["choices"], "request.choices")):
        path = f"request.choices[{index}]"
        choice = _object(value, {"recordKey", "kind", "anchor"}, path)
        choices.append(
            AlignmentRecordChoice(
                record_key=_text(choice["recordKey"], f"{path}.recordKey"),
                # AlignmentRecordChoice validates the text as a choice kind.
                kind=cast(AlignmentChoiceKind, _text(choice["kind"], f"{path}.kind")),
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
        "schema": SIMILARITY_ALIGNMENT_PLAN_SCHEMA,
        "status": "ambiguous" if resolution.ambiguities else "resolved",
        "groupId": resolution.group_id,
        "reference": _serialize_anchor(resolution.reference),
        "referenceDisplayedStrand": resolution.reference_candidate.displayed_strand,
        "referenceDisplayCenter": resolution.reference_candidate.display_center,
        "referenceSourceStart": resolution.reference_candidate.source_start,
        "referenceSourceEnd": resolution.reference_candidate.source_end,
        "records": [
            _serialize_ambiguity(record, row)
            if isinstance(record, AmbiguousAlignmentRecord)
            else _serialize_decision(record, row)
            for record, row in zip(
                resolution.records, resolution.review_rows, strict=True
            )
        ],
        "plan": _serialize_plan(plan) if plan is not None else None,
        "projection": (
            _serialize_projection(committed, _projection_geometry(committed), candidates, raw)
            if committed.raw.get("geometry", True) else None
        ),
    }


def resolve_similarity_alignment_json(
    payload_json: object, projection_json: object, resource_paths_json: object = "{}",
    output_directory: str = "/tmp",
) -> str:
    """JSON boundary used by the existing lazy diagram Worker."""

    return json.dumps(resolve_similarity_alignment_payload(
        json.loads(str(payload_json)), projection=json.loads(str(projection_json)),
        resource_paths=json.loads(str(resource_paths_json)), output_directory=output_directory,
    ))
