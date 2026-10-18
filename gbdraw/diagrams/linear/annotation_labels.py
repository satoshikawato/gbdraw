"""Where a Linear record's annotation labels go along their spans (OV-845).

An annotation label keeps the place ``linear_annotation_labels`` gives it
when its painted letters are clear of the record's feature labels, their
leaders and the annotation labels placed before it. Otherwise it moves along
its mark to the nearest clear place; with none, to the place with the fewest
conflicts, nearest first.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

from gbdraw.annotations import ResolvedAnnotationTrack, annotation_track_params_from_mapping
from gbdraw.labels.linear import calculate_label_bounds, linear_leader_polyline
from gbdraw.layout.spatial import Aabb
from gbdraw.render.drawers.linear.annotations import LinearAnnotationLabel, linear_annotation_labels

from .track_slots import LinearRecordVerticalPlan

AnnotationLabelsBySlot = Mapping[str, tuple[LinearAnnotationLabel, ...]]


def _feature_label_obstacles(
    feature_labels: Sequence[dict[str, object]],
    origin_y: float,
    leader_half_width: float,
) -> list[Aabb]:
    """The record's feature label boxes and leader lines, from the record axis."""

    obstacles: list[Aabb] = []
    for label in feature_labels:
        left, right, top, bottom = calculate_label_bounds(label)
        obstacles.append(Aabb(left, origin_y + top, right, origin_y + bottom))
        leader = linear_leader_polyline(label)
        for (x1, y1), (x2, y2) in zip(leader, leader[1:]):
            obstacles.append(
                Aabb(x1, origin_y + y1, x2, origin_y + y2).expanded(leader_half_width)
            )
    return obstacles


def _place_along_span(
    label: LinearAnnotationLabel, origin_y: float, obstacles: Sequence[Aabb], gap: float
) -> float:
    """The anchor x for ``label`` (in a slot at ``origin_y``): its own when clear, else the nearest clear x over its span."""

    ink = label.ink.translated(0.0, origin_y)
    left, right = ink.min_x - label.x, ink.max_x - label.x
    # Anchor x values that put the letters within ``gap`` of an obstacle in their rows.
    blocked = [
        (box.min_x - gap - right, box.max_x + gap - left)
        for box in obstacles
        if box.min_y < ink.max_y + gap and box.max_y > ink.min_y - gap
    ]

    def conflicts(x: float) -> int:
        return sum(1 for start, end in blocked if start < x < end)

    if conflicts(label.x) == 0:
        return label.x
    low, high = label.span[0] - left, label.span[1] - right
    if high < low:
        return label.x
    candidates = {low, high, *(value for interval in blocked for value in interval if low <= value <= high)}
    return min(sorted(candidates), key=lambda x: (conflicts(x), abs(x - label.x)))


def place_linear_annotation_labels(
    *,
    record_index: int,
    record_length: int,
    sequence_width: float,
    plan: LinearRecordVerticalPlan,
    annotation_track_layouts: Mapping[str, ResolvedAnnotationTrack],
    feature_labels: Sequence[dict[str, object]],
    font_family: str,
    leader_stroke_width: float,
    gap: float,
) -> dict[str, tuple[LinearAnnotationLabel, ...]]:
    """Each annotation slot's labels on one record, placed along their spans (slot-local y)."""

    features_origin_y = next(
        (float(slot.origin_y) for slot in plan.slots if slot.renderer == "features"), 0.0
    )
    obstacles = _feature_label_obstacles(feature_labels, features_origin_y, 0.5 * leader_stroke_width)
    placed: dict[str, tuple[LinearAnnotationLabel, ...]] = {}
    for slot in plan.slots:
        track = annotation_track_layouts.get(slot.id)
        if slot.renderer != "annotations" or track is None:
            continue
        labels = linear_annotation_labels(
            track,
            record_index=record_index,
            record_length=record_length,
            bar_length_px=sequence_width,
            side=slot.side,
            height_px=float(slot.height),
            font_family=font_family,
            params=annotation_track_params_from_mapping(slot.params),
        )
        slot_labels: list[LinearAnnotationLabel] = []
        origin_y = float(slot.origin_y)
        for label in labels:
            moved = label.moved_to(_place_along_span(label, origin_y, obstacles, gap))
            obstacles.append(moved.ink.translated(0.0, origin_y))
            slot_labels.append(moved)
        placed[slot.id] = tuple(slot_labels)
    return placed


__all__ = ["AnnotationLabelsBySlot", "place_linear_annotation_labels"]
