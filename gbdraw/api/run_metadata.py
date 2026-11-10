"""What a render reports beside its Results: the run metadata.

This module is the one producer of it for every surface. The Web Worker sends
:func:`build_run_metadata` with its Results (and adds the fields a Session does
not keep); the Session writers of the gbdraw command, for a run and for a
re-save of a Session, store it as it is, and a Web Save writes the values the
Worker sent, so both Sessions hold the same ``runMetadata`` for the same run
(OV-3240, OV-3244). The Web reads each row by ``resultIndex`` and checks ``resultName``
against the Result's name, so each Result is named here as its surface names
it.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING, Literal, NamedTuple

from gbdraw.render.track_slot_metadata import (
    collect_track_slot_geometry_records,
    resolved_track_slot_geometry,
)
from gbdraw.session_contract import SessionRunMetadata

if TYPE_CHECKING:
    from gbdraw.api.request_render import RequestRenderResult


class RenderedResult(NamedTuple):
    """One Result of a render and the name its surface gives it."""

    item: RequestRenderResult
    name: str


def build_run_metadata(
    mode: Literal["circular", "linear"],
    results: Sequence[RenderedResult],
) -> SessionRunMetadata:
    """The run metadata of the rendered Results, in their order: every notice
    list (empty when there is none), and ``trackSlotGeometry`` when a Result
    has any."""

    metadata: SessionRunMetadata = {}
    geometry = resolved_track_slot_geometry(
        mode=mode,
        records=[
            record
            for result_index, (item, name) in enumerate(results)
            for record in collect_track_slot_geometry_records(
                item.drawing, result_index=result_index, result_name=name
            )
        ],
    )
    if geometry is not None:
        metadata["trackSlotGeometry"] = geometry
    metadata.update({
        "annotationWarnings": [
            {
                "code": warning.code,
                "setId": warning.set_id,
                "annotationId": warning.annotation_id,
                "recordId": warning.record_id,
                "recordIndex": warning.record_index,
                "missingCount": warning.missing_count,
                "message": warning.message,
                "resultIndex": result_index,
                "resultName": name,
            }
            for result_index, (item, name) in enumerate(results)
            for warning in item.annotation_warnings
        ],
        "featureIdentityNotices": [
            {
                "recordKey": notice.record_key,
                "biologicalFeatureId": notice.biological_feature_id,
                "status": notice.status,
                "kinds": list(notice.kinds),
                "resultIndex": result_index,
            }
            for result_index, (item, _name) in enumerate(results)
            for notice in item.feature_identity_notices
        ],
        "comparisonWarnings": [
            {
                "code": warning.code,
                "queryRecordIndex": warning.query_record_index,
                "subjectRecordIndex": warning.subject_record_index,
                "queryRecordId": warning.query_record_id,
                "subjectRecordId": warning.subject_record_id,
                "rowCount": warning.row_count,
                "exampleIds": list(warning.example_ids),
                "message": warning.message,
                "resultIndex": result_index,
                "resultName": name,
            }
            for result_index, (item, name) in enumerate(results)
            if item.linear_metadata is not None
            for warning in item.linear_metadata.comparison_warnings
        ],
    })
    return metadata


__all__ = [
    "RenderedResult",
    "build_run_metadata",
]
