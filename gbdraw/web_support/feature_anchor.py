"""Worker helper of the Web's Align at feature (Linear).

Python owns where ``anchor_feature`` places each record: the Web sends its
committed request and the selector, and writes the record translations and
the display starts this helper returns once, as one History step. The helper
runs the planner and the drawing that ``gbdraw linear --anchor_feature`` runs,
so the starts and translations are the ones a CLI Session stores
(``resolved_plan_request``); what each record's status is, and which record is
the reference, is the planner's decision (``resolve_anchor_feature``).
"""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from dataclasses import replace
from typing import TypedDict

from gbdraw.api.options import LinearRecordTranslation
from gbdraw.api.request_render import plan_linear_request, with_anchor_feature_translations
from gbdraw.api.requests import LinearDiagramRequest
from gbdraw.exceptions import ValidationError
from gbdraw.layout.record_placement import FeatureAnchorStatus
from gbdraw.session_request_codec import (
    EncodedRecordTranslation,
    decode_canonical_request,
    encode_record_translations,
)


class FeatureAnchorRecord(TypedDict):
    """One drawn record in the reply: what Align at feature does to it."""

    recordKey: str
    recordId: str
    # rotated: a complete circular record turned so that the feature starts at
    # its left end (``startCoordinate`` is its new display start);
    # reference: when no record is rotated, the first record with the feature,
    # which keeps its place; aligned: moved sideways so that the feature starts
    # below the reference's; unmatched: no feature starts in its drawn region,
    # so it keeps its place and its start.
    status: FeatureAnchorStatus
    # The record's display start (a 1-based source base) in the resolved
    # request: the new start of a rotated record, else the start it had, or
    # null when it has none.
    startCoordinate: int | None
    beforeX: float
    afterX: float


class FeatureAnchorReply(TypedDict):
    """The reply: each record in drawing order, and the layout's translations to write."""

    selector: str
    records: list[FeatureAnchorRecord]
    recordTranslations: list[EncodedRecordTranslation]


def _translation_x(translations: Sequence[LinearRecordTranslation], key: str) -> float:
    return next((item.x for item in translations if item.record_key == key), 0.0)


def resolve_feature_anchor(request: LinearDiagramRequest, selector_text: str) -> FeatureAnchorReply:
    """The record translations and starts ``anchor_feature`` resolves ``request`` to, and what it does to each record.

    The planner's refusals (a selector that no record matches, matches at
    different starts in one record, records that share a row, a wrap that cuts
    a record, a similarity alignment) propagate as they are.
    """

    plan = plan_linear_request(replace(request, options=replace(request.options, anchor_feature=selector_text)))
    anchor = plan.feature_anchor
    assert anchor is not None  # The request sets the anchor, so the planner decided it.
    resolved = with_anchor_feature_translations(plan.request, plan.build().metadata)
    assert resolved.layout is not None  # with_anchor_feature_translations writes one.
    before = request.layout.record_translations if request.layout is not None else ()
    after = resolved.layout.record_translations
    records: list[FeatureAnchorRecord] = [
        {
            "recordKey": translation.record_key,
            "recordId": str(record.id),
            "status": status,
            "startCoordinate": provenance.display.start_coordinate,
            "beforeX": _translation_x(before, translation.record_key),
            "afterX": translation.x,
        }
        for record, provenance, translation, status in zip(
            plan.records, plan.provenance, after, anchor.statuses, strict=True,
        )
    ]
    return {
        "selector": selector_text.strip(),
        "records": records,
        "recordTranslations": encode_record_translations(after),
    }


def resolve_feature_anchor_json(
    request_json: str, selector_json: str, resource_paths_json: str, output_directory: str,
) -> str:
    """JSON boundary of the diagram Worker's ``resolveFeatureAnchor`` helper."""

    payload: object = json.loads(str(request_json))
    selector: object = json.loads(str(selector_json))
    if not isinstance(payload, Mapping) or not isinstance(selector, str):
        raise ValidationError(
            "Align at feature needs a canonical request and a selector text.",
            diagnostic={"code": "HELPER_PROTOCOL"},
        )
    request = decode_canonical_request(
        payload,
        resource_paths=json.loads(str(resource_paths_json)),
        output_directory=output_directory,
    )
    if not isinstance(request, LinearDiagramRequest):
        raise ValidationError(
            "Align at feature needs a Linear request.", diagnostic={"code": "HELPER_PROTOCOL"},
        )
    return json.dumps(resolve_feature_anchor(request, selector))
