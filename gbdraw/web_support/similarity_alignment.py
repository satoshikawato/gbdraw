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
    ResolvedRecordCollection, _source_details, alignment_record_keys, project_similarity_alignment_anchor_fact,
    similarity_alignment_evidence_edges,
)
from gbdraw.api.request_render import LinearRequestPlan, plan_linear_request, resolve_request_records
from gbdraw.api.requests import LinearDiagramRequest
from gbdraw.io.regions import RecordWindow, region_slice
from gbdraw.layout.record_coordinates import RecordDisplayTransform
from gbdraw.session_request_codec import StoredRecordWindow, _decode_windows, decode_canonical_request
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

    ``recordLength`` is the complete source length, ``None`` when the Web does
    not know it; ``windows`` are the record's windows as request schema 10
    stores them (``StoredRecordWindow``), empty for a whole record.
    """

    recordKey: str
    recordLength: int | None
    windows: list[StoredRecordWindow]
    presentation: AlignmentHelperPresentation


@dataclass(frozen=True)
class _RecordFact:
    """A validated helper request record: the windows it draws and its orientation."""

    record_key: str
    record_length: int | None
    windows: tuple[RecordWindow, ...]
    presentation_reverse: bool

    def placed(self, source_center: float) -> tuple[float | None, bool]:
        """Where the record draws the source boundary ``source_center``: the display offset on the drawn unit that
        holds it, and whether that unit is drawn reverse-complemented.

        A record drawn in several windows is anchored in the first window, in drawing order, that holds the position,
        as the renderer anchors it (``alignment_anchor_unit``, D-50). The offset is ``None`` when no window holds the
        position, or when the map needs the source length and it is unknown; the orientation is then the first
        window's.
        """

        if not self.windows:
            return self._unit_offset(None, source_center), self.presentation_reverse
        for window in self.windows:
            offset = self._unit_offset(window, source_center)
            if offset is not None:
                return offset, self.presentation_reverse or window.reverse_complement
        return None, self.presentation_reverse or self.windows[0].reverse_complement

    def _unit_offset(self, window: RecordWindow | None, source_center: float) -> float | None:
        """The display offset at which the unit of ``window`` (the whole record for ``None``) draws ``source_center``.

        The crop draws a window as a unit whose coordinate map starts at the
        window's start, or at its end when the unit is reversed (by the window
        or by the presentation, applied after the crop); the unit's
        ``RecordDisplayTransform`` maps the source onto its axis, round the
        origin for a window across it.
        """

        length = self.record_length
        if window is None:
            if not self.presentation_reverse:
                return source_center
            if length is None or source_center > length:
                return None
            transform = RecordDisplayTransform(length, length, -1, None, False)
            return transform.source_position_to_display_offset(source_center)
        reverse = self.presentation_reverse or window.reverse_complement
        if window.crosses_origin and length is None:
            return None
        bound = length if length is not None else window.end
        if not 0.0 <= source_center <= float(bound):
            return None
        span = bound - window.start + 1 + window.end if window.crosses_origin else window.end - window.start + 1
        transform = RecordDisplayTransform(
            bound, window.end if reverse else window.start, -1 if reverse else 1, None,
            window.crosses_origin, display_length=span, crosses_origin=window.crosses_origin,
        )
        offset = transform.source_position_to_display_offset(source_center)
        return offset if 0.0 <= offset <= span else None


def _record_fact(value: object, path: str) -> _RecordFact:
    raw = cast(AlignmentHelperRecord, _object(
        value,
        {"recordKey", "recordLength", "windows", "presentation"},
        path,
    ))
    record_key = _text(raw["recordKey"], f"{path}.recordKey")
    record_length = _integer(
        raw["recordLength"], f"{path}.recordLength", minimum=1, nullable=True
    )
    presentation = _object(
        raw["presentation"], {"reverseComplement"}, f"{path}.presentation"
    )
    presentation_reverse = _boolean(
        presentation["reverseComplement"], f"{path}.presentation.reverseComplement"
    )
    try:
        stored = _decode_windows(raw["windows"], path=f"{path}.windows")
    except ValueError as exc:
        raise ValidationError(f"{path}.windows is invalid: {exc}") from exc
    if any(window.reverse_complement for window in stored) and presentation_reverse:
        raise ValidationError(f"{path} cannot reverse complement both window and presentation.")
    windows = tuple(_drawn_window(window, record_length, record_key) for window in stored)
    return _RecordFact(record_key, record_length, windows, presentation_reverse)


def _drawn_window(window: RecordWindow, record_length: int | None, record_key: str) -> RecordWindow:
    """A window as the crop draws it: its bounds read by the crop's slice rule (``region_slice``).

    The crop clamps an End past the record end and refuses a window that
    starts past it or runs across the origin outside the record, so the helper
    places the record as the drawing places it. Without the record length the
    window is read as stored.
    """

    if record_length is None:
        return window
    _start_0, end_0 = region_slice(record_length, window.start, window.end, record_key)
    return window if window.crosses_origin or end_0 == window.end else replace(window, end=end_0)


def _member_candidate(
    value: object,
    path: str,
    *,
    group_id: str,
    record_facts: Mapping[str, _RecordFact],
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

    display_center, reverse = record.placed((start + end) / 2)
    displayed_strand = (
        None if source_strand is None
        else source_strand * (-1 if reverse else 1)
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
    keys = [record.record_key for record in request.records if record.record_key is not None]
    if len(set(keys)) != len(request.records):
        raise ValidationError("Alignment projection requires keyed records.")
    # The caller materializes current placements through the composition bridge.
    # The plan is deliberately removed for raw, untranslated placement facts.
    request = replace(request, similarity_alignment=None)
    before = plan_linear_request(request)
    before_collection = ResolvedRecordCollection(before.records, before.provenance)
    if alignment_record_keys(before_collection) != tuple(keys):
        raise ValidationError("Alignment projection record binding changed.", diagnostic={"code": "INPUT_INVALID"})
    # Counterfactual source transforms are always available for local direction edits.
    flipped_orientations = {key: not _drawn_reverse(before, key) for key in keys}
    orientations = raw["orientations"]
    if orientations is None:
        orientations = flipped_orientations
    else:
        orientations = _object(orientations, set(keys), "projection.orientations")
        orientations = {key: _boolean(value, f"projection.orientations.{key}")
                        for key, value in orientations.items()}
    after = plan_linear_request(_oriented_request(request, orientations))
    # The flipped records serve only anchor facts, so they are resolved, not
    # planned; without orientations (the stale-target repair) they are the
    # requested ones (OV-1585).
    flipped_collection = (
        ResolvedRecordCollection(after.records, after.provenance)
        if orientations == flipped_orientations
        else resolve_request_records(_oriented_request(request, flipped_orientations))
    )
    return raw, before, before_collection, after, flipped_collection, request.options.orthogroups, request


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


def _serialize_projection(context, candidates, raw_request):
    raw, before, collection, after, flipped_collection, _orthogroups, request = context
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
    for index, key in enumerate(item.record_key for item in request.records):
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

def _resource_evidence_edges(context, group_id, candidates) -> list[AlignmentEvidenceEdge]:
    """Bind resource edges to the current members; unbound endpoints are omitted."""

    collection, orthogroups = context[2], context[5]
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
    if type(raw["schema"]) is not int or raw["schema"] != SIMILARITY_ALIGNMENT_PLAN_SCHEMA:
        raise ValidationError(f"request.schema must be {SIMILARITY_ALIGNMENT_PLAN_SCHEMA}.")
    group_id = _text(raw["groupId"], "request.groupId")
    records = [
        _record_fact(value, f"request.records[{index}]")
        for index, value in enumerate(_list(raw["records"], "request.records"))
    ]
    record_keys = [record.record_key for record in records]
    if not record_keys or len(record_keys) != len(set(record_keys)):
        raise ValidationError("request.records must have unique record keys.")
    record_facts = {record.record_key: record for record in records}
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
        if alignment_record_keys(collection) != tuple(record_keys):
            raise ValidationError("Alignment projection record binding changed.", diagnostic={"code": "INPUT_INVALID"})
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
    direct_edges = _list(raw["directEdges"], "request.directEdges")
    if context is not None:
        # The committed request's orthogroup resource is the only evidence
        # authority; the UI catalog carries no graph.
        if direct_edges:
            raise ValidationError(
                "request.directEdges must be empty; direct evidence comes from "
                "the committed orthogroup resource."
            )
        edges = _resource_evidence_edges(context, group_id, candidates)
    for index, value in enumerate(direct_edges):
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
