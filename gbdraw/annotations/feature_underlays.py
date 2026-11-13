"""Private feature-underlay annotation materialization helpers."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import replace
from typing import Any

from Bio.SeqRecord import SeqRecord

from gbdraw.exceptions import ValidationError
from gbdraw.features.ids import make_linear_rendered_feature_id
from gbdraw.features.objects import FeatureObject
from gbdraw.layout.record_coordinates import RecordDisplayTransform
from gbdraw.svg.ids import instance_svg_id
from gbdraw.tracks.parsing import SlotView, slots_draw_features

from .models import (
    RegionAnnotationStyle,
    ResolvedAnnotationBundle,
    ResolvedRegionAnnotation,
)
from .resolve import annotation_midpoint, merge_annotation_segments


AUTO_FEATURE_UNDERLAY_KIND = "gbdraw:auto-feature-underlay"
AUTO_FEATURE_STABLE_ID = "gbdraw:stable-feature-id"
AUTO_FEATURE_RENDERED_ID = "gbdraw:rendered-feature-id"
AUTO_FEATURE_DOM_ID = "gbdraw:dom-feature-id"
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
    record_transforms: Sequence[RecordDisplayTransform | None] | None = None,
) -> tuple[ResolvedAnnotationBundle, str | None]:
    """Append private resolved highlights without mutating the public bundle. ``record_transforms``, by record,
    are those that projected the features' display parts; a record on a circle it shares with other records
    (``RecordDisplayTransform.circle``) ends its underlays at its display cut (``mark_parts``)."""

    if mode not in {"circular", "linear"}:
        raise ValidationError("Feature underlay mode must be 'circular' or 'linear'.")
    if len(underlay_features_by_record) != len(records):
        raise ValueError("underlay_features_by_record must match records")
    indices = tuple(record_indices or range(len(records)))
    if len(indices) != len(records):
        raise ValueError("record_indices must match records")
    transforms = tuple(record_transforms or (None,) * len(records))
    if len(transforms) != len(records):
        raise ValueError("record_transforms must match records")
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
        transform = transforms[local_index]
        stable_ids = [
            feature.feature_hash or str(feature.feature_id)
            for feature in features
        ]
        for feature_index, feature in enumerate(features):
            segments = _feature_segments(feature)
            if not segments:
                continue
            stable_id = stable_ids[feature_index]
            rendered_id = (
                make_linear_rendered_feature_id(
                    record_index=record_index,
                    stable_feature_id=stable_id,
                    record_count=record_count,
                )
                if mode == "linear"
                else stable_id
            ) or stable_id
            dom_id = (
                instance_svg_id(rendered_id, feature.rendered_instance_index)
                if feature.rendered_instance_index is not None
                else rendered_id
            )
            annotation_id = (
                f"record_{record_index + 1}_{feature.feature_id}_{feature_index + 1}"
            )
            fill = (
                str(getattr(feature, "color", "") or "").strip()
                or AUTO_FEATURE_UNDERLAY_FILL
            )
            resolved.append(
                ResolvedRegionAnnotation(
                    id=annotation_id,
                    set_id=set_id,
                    record_index=record_index,
                    segments=segments,
                    display_parts=(tuple(part.fragment for part in feature.display_parts if part.kind == "block")
                                   if feature.display_parts is not None else None),
                    shares_circle=transform is not None and transform.circle is not None,
                    midpoint_bp=annotation_midpoint(segments, record_length, crosses_origin=False),
                    span_bp=sum(end - start for start, end in segments),
                    label="",
                    mark="highlight",
                    lane=None,
                    style=RegionAnnotationStyle(fill=fill),
                    style_is_annotation_override=False,
                    legend_label=None,
                    metadata={
                        AUTO_FEATURE_UNDERLAY_KIND: "true",
                        AUTO_FEATURE_STABLE_ID: stable_id,
                        AUTO_FEATURE_RENDERED_ID: rendered_id,
                        AUTO_FEATURE_DOM_ID: dom_id,
                    },
                )
            )
    if not resolved:
        return bundle, None
    return (
        replace(bundle, annotations=(*bundle.annotations, *resolved), set_ids=(*bundle.set_ids, set_id)),
        set_id,
    )


def feature_underlay_anchor_slot_id(
    slots: Sequence[SlotView], underlay_features: Iterable[FeatureObject]
) -> str | None:
    """The enabled features slot the drawn underlays lie under, or ``None`` when the stack draws none.

    A stack draws its feature underlays when it draws its features (:func:`slots_draw_features`, the answer
    the Legend and the input summary read): without an enabled features slot it draws no feature one by one,
    underlays included. With two or more there is no one slot to lie under, and the drawing is refused.
    """

    if not slots_draw_features(slots):
        return None
    feature_slots = [slot for slot in slots if slot.enabled and str(slot.renderer) == "features"]
    if len(feature_slots) > 1:
        types = ", ".join(sorted({str(feature.feature_type) for feature in underlay_features}))
        raise ValidationError(
            f"Feature types drawn as underlay ({types}) need exactly one enabled features track slot; "
            f"{len(feature_slots)} are enabled. Keep one features slot, or draw those types as arrow or "
            "rectangle (--feature_shape TYPE=rectangle).",
            diagnostic={"code": "TRACK_INVALID", "reason": "UNDERLAY_FEATURE_SLOTS"},
        )
    return str(feature_slots[0].id)


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
    dom_id = annotation.metadata.get(AUTO_FEATURE_DOM_ID, "") or rendered_id
    stable_id = annotation.metadata.get(AUTO_FEATURE_STABLE_ID, "")
    if not rendered_id:
        return
    element.attribs["id"] = (
        dom_id if part_count <= 1 else f"{dom_id}__part{part_index}"
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
