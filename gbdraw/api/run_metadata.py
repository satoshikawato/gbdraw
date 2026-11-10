"""What a render reports beside its Results: the run metadata.

This module is the one producer of it for every surface. The Web Worker sends
:func:`build_run_metadata` with its Results (and adds the fields a Session does
not keep); the Session writers of the gbdraw command, for a run and for a
re-save of a Session, store :func:`saved_run_metadata` of it, so a Load of
either Session reads what a Web Save of the same figure holds (OV-3240,
OV-3244). The Web reads each row by ``resultIndex`` and checks ``resultName``
against the Result's name, so each Result is named here as its surface names
it.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING, Literal, NamedTuple, TypedDict

from gbdraw.render.track_slot_metadata import (
    TrackSlotGeometry,
    collect_track_slot_geometry_records,
    resolved_track_slot_geometry,
)

if TYPE_CHECKING:
    from gbdraw.api.request_render import RequestRenderResult


class RunAnnotationWarning(TypedDict):
    """A region annotation that a Result does not draw whole (``ResolutionWarning``)."""

    code: str
    setId: str
    annotationId: str
    recordId: str
    recordIndex: int
    missingCount: int
    message: str
    resultIndex: int
    resultName: str


class RunFeatureIdentityNotice(TypedDict):
    """An identity-addressed feature edit that a Result does not draw (``FeatureIdentityNotice``)."""

    recordKey: str
    biologicalFeatureId: str
    status: Literal["crop_excluded", "absent", "unresolved"]
    kinds: list[str]
    resultIndex: int


class RunComparisonWarning(TypedDict):
    """Comparison table rows whose sequence IDs name no drawn record (``ComparisonRecordIdWarning``)."""

    code: str
    queryRecordIndex: int
    subjectRecordIndex: int
    queryRecordId: str
    subjectRecordId: str
    rowCount: int
    exampleIds: list[str]
    message: str
    resultIndex: int
    resultName: str


class RunMetadata(TypedDict, total=False):
    """The run metadata of one render.

    :func:`build_run_metadata` writes every notice list (empty when there is
    none) and ``trackSlotGeometry`` when a Result has any; the saved form
    (:func:`saved_run_metadata`) keeps only the fields that hold something.
    """

    trackSlotGeometry: TrackSlotGeometry
    annotationWarnings: list[RunAnnotationWarning]
    featureIdentityNotices: list[RunFeatureIdentityNotice]
    comparisonWarnings: list[RunComparisonWarning]


class RenderedResult(NamedTuple):
    """One Result of a render and the name its surface gives it."""

    item: RequestRenderResult
    name: str


def build_run_metadata(
    mode: Literal["circular", "linear"],
    results: Sequence[RenderedResult],
) -> RunMetadata:
    """The run metadata of the rendered Results, in their order."""

    metadata: RunMetadata = {}
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


def saved_run_metadata(metadata: RunMetadata) -> RunMetadata:
    """What a Session written by gbdraw keeps of the run metadata: the fields
    that hold something. A Load reads a missing notice list as an empty one."""

    saved: RunMetadata = {}
    if "trackSlotGeometry" in metadata:
        saved["trackSlotGeometry"] = metadata["trackSlotGeometry"]
    if metadata.get("annotationWarnings"):
        saved["annotationWarnings"] = metadata["annotationWarnings"]
    if metadata.get("featureIdentityNotices"):
        saved["featureIdentityNotices"] = metadata["featureIdentityNotices"]
    if metadata.get("comparisonWarnings"):
        saved["comparisonWarnings"] = metadata["comparisonWarnings"]
    return saved


__all__ = [
    "RenderedResult",
    "RunAnnotationWarning",
    "RunComparisonWarning",
    "RunFeatureIdentityNotice",
    "RunMetadata",
    "build_run_metadata",
    "saved_run_metadata",
]
