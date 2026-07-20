"""Private feature-underlay annotation materialization helpers."""

from __future__ import annotations

from collections.abc import Sequence
from typing import Any

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]

from gbdraw.exceptions import ValidationError
from gbdraw.features.ids import (
    compute_feature_object_hash,
    make_linear_rendered_feature_id,
)
from gbdraw.features.objects import FeatureObject

from .models import (
    RegionAnnotationStyle,
    ResolvedAnnotationBundle,
    ResolvedRegionAnnotation,
)
from .resolve import annotation_midpoint, merge_annotation_segments


AUTO_FEATURE_UNDERLAY_KIND = "gbdraw:auto-feature-underlay"
AUTO_FEATURE_STABLE_ID = "gbdraw:stable-feature-id"
AUTO_FEATURE_RENDERED_ID = "gbdraw:rendered-feature-id"
AUTO_FEATURE_UNDERLAY_SET_ID = "__gbdraw_auto_feature_underlay__"
AUTO_FEATURE_UNDERLAY_SLOT_ID = "__gbdraw_auto_feature_underlay_slot__"
AUTO_FEATURE_UNDERLAY_FILL = "#808080"


def _collision_safe_id(preferred: str, occupied: set[str]) -> str:
    candidate = str(preferred)
    suffix = 2
    while candidate in occupied:
        candidate = f"{preferred}_{suffix}"
        suffix += 1
    return candidate


def _feature_segments(feature: FeatureObject) -> tuple[tuple[int, int], ...]:
    segments: list[tuple[int, int]] = []
    for part in getattr(feature, "coordinates", ()) or ():
        try:
            start, end = int(part.start), int(part.end)
        except (AttributeError, TypeError, ValueError):
            continue
        if end > start:
            segments.append((start, end))
    return merge_annotation_segments(segments)


def merge_feature_underlays(
    bundle: ResolvedAnnotationBundle,
    underlay_features_by_record: Sequence[Sequence[FeatureObject]],
    records: Sequence[SeqRecord],
    *,
    mode: str,
    record_indices: Sequence[int] | None = None,
) -> tuple[ResolvedAnnotationBundle, str | None]:
    """Append private resolved highlights without mutating the public bundle."""

    if mode not in {"circular", "linear"}:
        raise ValidationError("Feature underlay mode must be 'circular' or 'linear'.")
    if len(underlay_features_by_record) != len(records):
        raise ValueError("underlay_features_by_record must match records")
    indices = tuple(record_indices or range(len(records)))
    if len(indices) != len(records):
        raise ValueError("record_indices must match records")
    if not any(underlay_features_by_record):
        return bundle, None

    occupied_set_ids = set(bundle.set_ids)
    occupied_set_ids.update(item.set_id for item in bundle.annotations)
    set_id = _collision_safe_id(AUTO_FEATURE_UNDERLAY_SET_ID, occupied_set_ids)
    resolved: list[ResolvedRegionAnnotation] = []
    record_count = max(indices, default=-1) + 1
    for local_index, (record, features) in enumerate(
        zip(records, underlay_features_by_record)
    ):
        record_index = int(indices[local_index])
        record_length = len(record.seq)
        for feature_index, feature in enumerate(features):
            segments = _feature_segments(feature)
            if not segments:
                continue
            stable_id = compute_feature_object_hash(feature) or str(feature.feature_id)
            rendered_id = (
                make_linear_rendered_feature_id(
                    record_index=record_index,
                    stable_feature_id=stable_id,
                    record_count=record_count,
                )
                if mode == "linear"
                else stable_id
            ) or stable_id
            annotation_id = (
                f"record_{record_index + 1}_{feature.feature_id}_{feature_index + 1}"
            )
            resolved.append(
                ResolvedRegionAnnotation(
                    id=annotation_id,
                    set_id=set_id,
                    record_index=record_index,
                    segments=segments,
                    midpoint_bp=annotation_midpoint(segments, record_length),
                    span_bp=sum(end - start for start, end in segments),
                    label="",
                    mark="highlight",
                    lane=None,
                    style=RegionAnnotationStyle(fill=AUTO_FEATURE_UNDERLAY_FILL),
                    style_is_annotation_override=False,
                    legend_label=None,
                    metadata={
                        AUTO_FEATURE_UNDERLAY_KIND: "true",
                        AUTO_FEATURE_STABLE_ID: stable_id,
                        AUTO_FEATURE_RENDERED_ID: rendered_id,
                    },
                )
            )
    if not resolved:
        return bundle, None
    return (
        ResolvedAnnotationBundle(
            annotations=(*bundle.annotations, *resolved),
            warnings=bundle.warnings,
            set_ids=(*bundle.set_ids, set_id),
        ),
        set_id,
    )


def feature_underlay_anchor_slot_id(slots: Sequence[object]) -> str:
    """Resolve the single enabled feature renderer used as the underlay anchor."""

    feature_slots = [
        slot
        for slot in slots
        if bool(getattr(slot, "enabled", True))
        and str(getattr(slot, "renderer", "")).strip().lower() == "features"
    ]
    if not feature_slots:
        raise ValidationError(
            "Visible feature underlays require an enabled features track slot."
        )
    if len(feature_slots) != 1:
        raise ValidationError(
            "Visible feature underlays require exactly one enabled features track slot."
        )
    return str(getattr(feature_slots[0], "id"))


def feature_underlay_slot_id(slots: Sequence[object]) -> str:
    """Return a deterministic private slot id that cannot collide with user slots."""

    occupied = {str(getattr(slot, "id", "")) for slot in slots}
    return _collision_safe_id(AUTO_FEATURE_UNDERLAY_SLOT_ID, occupied)


def is_auto_feature_underlay(annotation: ResolvedRegionAnnotation) -> bool:
    return annotation.metadata.get(AUTO_FEATURE_UNDERLAY_KIND) == "true"


def apply_feature_underlay_dom_attributes(
    element: Any,
    annotation: ResolvedRegionAnnotation,
    *,
    record_id: str,
    record_index: int,
    part_index: int,
    part_count: int,
    include_stable_id: bool,
) -> None:
    """Mark one private annotation shape as a feature-owned SVG block."""

    rendered_id = annotation.metadata.get(AUTO_FEATURE_RENDERED_ID, "")
    stable_id = annotation.metadata.get(AUTO_FEATURE_STABLE_ID, "")
    if not rendered_id:
        return
    element.attribs["id"] = (
        rendered_id if part_count <= 1 else f"{rendered_id}__part{part_index}"
    )
    element.attribs.update(
        {
            "data-gbdraw-feature-id": rendered_id,
            "data-gbdraw-record-id": str(record_id),
            "data-gbdraw-record-index": str(int(record_index)),
            "data-gbdraw-feature-part": "block",
            "data-gbdraw-auto-feature-underlay": "true",
        }
    )
    if include_stable_id and stable_id:
        element.attribs["data-gbdraw-stable-feature-id"] = stable_id


__all__ = [
    "AUTO_FEATURE_RENDERED_ID",
    "AUTO_FEATURE_STABLE_ID",
    "AUTO_FEATURE_UNDERLAY_KIND",
    "AUTO_FEATURE_UNDERLAY_SET_ID",
    "AUTO_FEATURE_UNDERLAY_SLOT_ID",
    "apply_feature_underlay_dom_attributes",
    "feature_underlay_anchor_slot_id",
    "feature_underlay_slot_id",
    "is_auto_feature_underlay",
    "merge_feature_underlays",
]
