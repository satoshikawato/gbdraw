"""Deterministic radial placement for external circular feature labels."""

from __future__ import annotations

import math
from collections.abc import Sequence
from typing import cast

from ..exceptions import LabelLayoutError, ValidationError
from ..layout.spatial import Aabb, AabbIndex
from ..layout.text_geometry import (
    OrientedTextBox,
    oriented_text_box,
    segment_intersects_polygon,
)
from ..svg.circular_arcs import position_turn_angle, position_turn_rad
from .circular_leaders import leader_crossing_count, leader_segments
from .fit import overlapping_pairs
from .circular_types import CircularLabel, CircularLabelLayoutResult

_TAU = 2.0 * math.pi
_ANGLE_EPSILON = 1e-9
_RADIAL_CLEARANCE_PX = 4.0
# OV-3841: the slanted part of a radial leader rises at least this far from the feature band.
_LEADER_MIN_SLOPE_RAD = math.radians(20.0)
_LEADER_ROOM_TOLERANCE_PX = 1e-3
_CONVERGENCE_LIMIT = 256
# Radius growth at or below this is none: the assembly stops growing the
# canvas, and an inner side waits for growth only above it (OV-694).
RADIAL_GROWTH_TOLERANCE_PX = 1e-3


def _normalize_readable_rotation(rotation_deg: float) -> float:
    normalized = ((float(rotation_deg) + 180.0) % 360.0) - 180.0
    if normalized > 90.0:
        normalized -= 180.0
    elif normalized < -90.0:
        normalized += 180.0
    return normalized


def radial_text_orientation(
    placed_angle_deg: float,
    *,
    is_inner: bool,
) -> tuple[float, str, str]:
    """Return readable SVG rotation, anchor and baseline for a radial label.

    ``placed_angle_deg`` uses genome convention: 0 degrees is twelve o'clock
    and values increase clockwise.
    """
    radial_rotation = (float(placed_angle_deg) % 360.0) - 90.0
    radial_x = math.cos(math.radians(radial_rotation))
    # Treat the exact top/bottom axes as part of the right half-plane. This
    # fixed tie-break prevents tiny input noise from flipping text 180 degrees.
    right_half = radial_x >= -_ANGLE_EPSILON
    if right_half:
        rotation = _normalize_readable_rotation(radial_rotation)
        text_anchor = "end" if is_inner else "start"
    else:
        rotation = _normalize_readable_rotation(radial_rotation + 180.0)
        text_anchor = "start" if is_inner else "end"
    return rotation, text_anchor, "middle"


def _preferred_middle(label: CircularLabel) -> float:
    return float(label.get("preferred_middle", label.get("middle", 0.0)))


def _preferred_angle_rad(label: CircularLabel, total_length: int) -> float:
    return position_turn_rad(_preferred_middle(label), total_length) % _TAU


def _radius_of_point(x_coord: float, y_coord: float) -> float:
    return math.hypot(float(x_coord), float(y_coord))


def _leader_anchor(label: CircularLabel, elbow: tuple[float, float]) -> tuple[float, float]:
    """Where a radial leader ends: the point of its feature's outline edge facing the label nearest ``elbow``."""
    points = label.get("feature_edge_points")
    if not points:
        return float(label.get("feature_anchor_x", 0.0)), float(label.get("feature_anchor_y", 0.0))
    return min(points, key=lambda point: math.dist(point, elbow))


def _leader_room_text_radius(label: CircularLabel, placed_angle_rad: float, *, is_inner: bool) -> float:
    """The text radius whose leader keeps outside the feature band (OV-3841), for a label placed at ``placed_angle_rad``.

    The leader runs from the elbow (``_RADIAL_CLEARANCE_PX`` from the text, at the label's angle) to the point of its
    feature's facing outline nearest the elbow. The elbow keeps ``_RADIAL_CLEARANCE_PX`` off the lane edge facing the
    label, and away from the leader's end by a radial drop that lets the slanted part rise at least
    ``_LEADER_MIN_SLOPE_RAD`` from the band; without it, leaders of a crowded cluster ran side by side along the band
    or crossed. Features are drawn over leaders, so a leader inside the band is hidden.
    """
    radial_angle = placed_angle_rad - (0.5 * math.pi)
    lane_edge = float(label["feature_inner_radius_px"] if is_inner else label["feature_outer_radius_px"])
    sign = -1.0 if is_inner else 1.0
    elbow_radius = lane_edge + sign * _RADIAL_CLEARANCE_PX
    anchor = _leader_anchor(label, (elbow_radius * math.cos(radial_angle), elbow_radius * math.sin(radial_angle)))
    anchor_radius = math.hypot(*anchor)
    turn = abs(((math.atan2(anchor[1], anchor[0]) - radial_angle + math.pi) % _TAU) - math.pi)
    drop = max(
        _RADIAL_CLEARANCE_PX,
        math.tan(_LEADER_MIN_SLOPE_RAD) * anchor_radius * math.sin(min(turn, 0.5 * math.pi)),
    )
    if is_inner:
        elbow_radius = min(elbow_radius, anchor_radius - drop)
    else:
        elbow_radius = max(elbow_radius, anchor_radius + drop)
    return elbow_radius + sign * _RADIAL_CLEARANCE_PX


def _outer_base_anchor(label: CircularLabel, outer_reserved_radius_px: float, total_length: int) -> float:
    elbow_radius = _radius_of_point(label.get("middle_x", 0.0), label.get("middle_y", 0.0))
    configured_minimum = float(label.get("min_outer_start_radius_px", 0.0))
    return max(
        elbow_radius + _RADIAL_CLEARANCE_PX,
        _leader_room_text_radius(label, _preferred_angle_rad(label, total_length), is_inner=False),
        configured_minimum,
        float(outer_reserved_radius_px) + _RADIAL_CLEARANCE_PX,
    )


def _inner_base_anchor(label: CircularLabel, total_length: int) -> float:
    elbow_radius = _radius_of_point(label.get("middle_x", 0.0), label.get("middle_y", 0.0))
    band_clear = _leader_room_text_radius(label, _preferred_angle_rad(label, total_length), is_inner=True)
    configured_maximum = float(label.get("max_inner_start_radius_px", band_clear))
    return max(0.0, min(elbow_radius - _RADIAL_CLEARANCE_PX, band_clear, configured_maximum))


def _seam_order(labels: Sequence[CircularLabel], total_length: int) -> list[CircularLabel]:
    ordered = sorted(
        labels,
        key=lambda label: (
            _preferred_angle_rad(label, total_length),
            str(label.get("stable_id", label.get("feature_id", ""))),
            int(label.get("input_order", 0)),
        ),
    )
    if len(ordered) < 2:
        return ordered
    angles = [_preferred_angle_rad(label, total_length) for label in ordered]
    gaps = [
        ((angles[(index + 1) % len(angles)] - angles[index]) % _TAU)
        for index in range(len(angles))
    ]
    seam_after = max(range(len(gaps)), key=lambda index: (gaps[index], -index))
    start = (seam_after + 1) % len(ordered)
    return ordered[start:] + ordered[:start]


def _pair_clearance_px(first: CircularLabel, second: CircularLabel, spacing_px: float) -> float:
    return (
        0.5 * float(first.get("height_px", 0.0))
        + 0.5 * float(second.get("height_px", 0.0))
        + max(0.0, float(spacing_px))
        + 0.75
    )


def _required_collision_radius(labels: Sequence[CircularLabel], spacing_px: float) -> float:
    if len(labels) < 2:
        return 0.0
    circumference = sum(
        _pair_clearance_px(labels[index - 1], labels[index], spacing_px)
        for index in range(len(labels))
    )
    return circumference / _TAU


def _pack_angles(
    labels: Sequence[CircularLabel],
    total_length: int,
    collision_radius_px: float,
    spacing_px: float,
) -> list[float]:
    """The placed angles of ``labels``: the midpoint of the clockwise and the counterclockwise packing.

    Clockwise packing alone keeps a crowded cluster's first label on its feature and pushes the others clockwise;
    counterclockwise packing does the mirror. Both obey the same cyclic gap constraints, so their midpoint does too
    and keeps the text order, and it spreads a cluster to both sides of its features. It is the one rule for every
    cluster: packing that ran clockwise alone unless its leaders crossed turned a whole cluster by degrees when a
    leader end moved by a pixel (OV-3648).
    """
    if not labels:
        return []
    if len(labels) == 1:
        return [_preferred_angle_rad(labels[0], total_length)]

    radius = max(1.0, float(collision_radius_px))
    preferred: list[float] = []
    wrap_offset = 0.0
    previous_raw: float | None = None
    for label in labels:
        raw_angle = _preferred_angle_rad(label, total_length)
        if previous_raw is not None and raw_angle < previous_raw - _ANGLE_EPSILON:
            wrap_offset += _TAU
        preferred.append(raw_angle + wrap_offset)
        previous_raw = raw_angle
    gaps = [
        _pair_clearance_px(labels[index - 1], labels[index], spacing_px) / radius
        for index in range(len(labels))
    ]

    start_angle = preferred[0]
    placed: list[float] = []
    for _ in range((4 * len(labels)) + 4):
        placed = [start_angle]
        for index in range(1, len(labels)):
            placed.append(max(preferred[index], placed[index - 1] + gaps[index]))
        closing_overflow = (placed[-1] + gaps[0]) - (placed[0] + _TAU)
        if closing_overflow <= _ANGLE_EPSILON:
            break
        start_angle += closing_overflow
    else:  # pragma: no cover - arithmetic guard for non-finite input
        raise LabelLayoutError("radial label cyclic relaxation did not converge")

    end_angle = preferred[-1]
    for _ in range((4 * len(labels)) + 4):
        backward = [end_angle]
        for index in range(len(labels) - 2, -1, -1):
            backward.append(min(preferred[index], backward[-1] - gaps[index + 1]))
        backward.reverse()
        overflow = backward[-1] + gaps[0] - backward[0] - _TAU
        if overflow <= _ANGLE_EPSILON:
            break
        end_angle -= overflow
    else:  # pragma: no cover - arithmetic guard for non-finite input
        raise LabelLayoutError("radial label cyclic relaxation did not converge")

    return [((first + second) / 2.0) % _TAU for first, second in zip(placed, backward, strict=True)]


def _make_placed_label(
    source: CircularLabel,
    *,
    total_length: int,
    placed_angle_rad: float,
    text_anchor_radius_px: float,
) -> CircularLabel:
    label = source.copy()
    is_inner = bool(label.get("is_inner", False))
    placed_angle_deg = math.degrees(placed_angle_rad) % 360.0
    radial_angle = placed_angle_rad - (0.5 * math.pi)
    radial_x = math.cos(radial_angle)
    radial_y = math.sin(radial_angle)
    text_x = float(text_anchor_radius_px) * radial_x
    text_y = float(text_anchor_radius_px) * radial_y
    rotation_deg, text_anchor, dominant_baseline = radial_text_orientation(
        placed_angle_deg,
        is_inner=is_inner,
    )
    box = oriented_text_box(
        text_x,
        text_y,
        float(label.get("width_px", 0.0)),
        float(label.get("height_px", 0.0)),
        rotation_deg,
        text_anchor,
    )

    if is_inner:
        elbow_radius = text_anchor_radius_px + _RADIAL_CLEARANCE_PX
    else:
        elbow_radius = max(0.0, text_anchor_radius_px - _RADIAL_CLEARANCE_PX)
    elbow_x = elbow_radius * radial_x
    elbow_y = elbow_radius * radial_y
    anchor_x, anchor_y = _leader_anchor(source, (elbow_x, elbow_y))

    label.update(
        {
            "placement": "radial",
            "preferred_middle": _preferred_middle(source),
            "preferred_angle_deg": position_turn_angle(_preferred_middle(source), total_length) % 360.0,
            "placed_angle_deg": placed_angle_deg,
            "middle": (placed_angle_deg / 360.0) * float(total_length),
            "start_x": text_x,
            "start_y": text_y,
            "text_x": text_x,
            "text_y": text_y,
            "rotation_deg": rotation_deg,
            "text_anchor": text_anchor,
            "dominant_baseline": dominant_baseline,
            "middle_x": elbow_x,
            "middle_y": elbow_y,
            "feature_anchor_x": anchor_x,
            "feature_anchor_y": anchor_y,
            "leader_start_x": text_x,
            "leader_start_y": text_y,
            "oriented_corners": box.corners,
            "aabb_local": box.aabb,
        }
    )
    return label


def _box_for_label(label: CircularLabel, padding_px: float = 0.0) -> OrientedTextBox:
    return oriented_text_box(
        float(label.get("text_x", label.get("start_x", 0.0))),
        float(label.get("text_y", label.get("start_y", 0.0))),
        float(label.get("width_px", 0.0)),
        float(label.get("height_px", 0.0)),
        float(label.get("rotation_deg", 0.0)),
        str(label.get("text_anchor", "middle")),
        padding_px=padding_px,
    )


def _text_collision_pairs(labels: Sequence[CircularLabel], spacing_px: float) -> list[tuple[int, int]]:
    padding = 0.5 * max(0.0, float(spacing_px))
    return overlapping_pairs(
        [(_box_for_label(label, padding_px=padding).corners,) for label in labels],
        touching=False,
    )


def _leader_text_collision_count(labels: Sequence[CircularLabel]) -> int:
    boxes = [_box_for_label(label) for label in labels]
    box_index = AabbIndex(bucket_size=64.0)
    for index, box in enumerate(boxes):
        box_index.insert(index, Aabb(*box.aabb))
    count = 0
    for label_index, label in enumerate(labels):
        for segment in leader_segments(label):
            for candidate in box_index.query_segment(*segment):
                # Every id inserted into box_index above is a label index (int).
                candidate_index = cast(int, candidate)
                if candidate_index == label_index:
                    continue
                if segment_intersects_polygon(segment, boxes[candidate_index].corners, touching=False):
                    count += 1
    return count


def _order_violation_count(labels: Sequence[CircularLabel]) -> int:
    if len(labels) < 2:
        return 0
    preferred = sorted(range(len(labels)), key=lambda index: float(labels[index]["preferred_angle_deg"]))
    placed = sorted(range(len(labels)), key=lambda index: float(labels[index]["placed_angle_deg"]))
    if not preferred:
        return 0
    doubled = placed + placed
    return 0 if any(doubled[offset : offset + len(preferred)] == preferred for offset in range(len(placed))) else 1


def _content_bounds(labels: Sequence[CircularLabel]) -> tuple[float, float, float, float] | None:
    if not labels:
        return None
    points: list[tuple[float, float]] = []
    for label in labels:
        points.extend(tuple(label.get("oriented_corners", ())))
        for segment in leader_segments(label):
            points.extend(segment)
    if not points:
        return None
    return (
        min(point[0] for point in points),
        min(point[1] for point in points),
        max(point[0] for point in points),
        max(point[1] for point in points),
    )


def _place_side(
    labels: Sequence[CircularLabel],
    *,
    total_length: int,
    spacing_px: float,
    is_inner: bool,
    outer_reserved_radius_px: float,
    inner_reserved_outer_radius_px: float,
) -> tuple[list[CircularLabel], float]:
    copied = [label.copy() for label in labels]
    if not copied:
        return [], 0.0
    ordered = _seam_order(copied, total_length)
    max_width = max(float(label.get("width_px", 0.0)) for label in ordered)
    density_radius = _required_collision_radius(ordered, spacing_px)
    required_growth = 0.0
    minimum_anchor = 0.0
    if is_inner:
        base_anchors = [_inner_base_anchor(label, total_length) for label in ordered]
        minimum_anchor = max(
            float(inner_reserved_outer_radius_px) + max_width + _RADIAL_CLEARANCE_PX,
            density_radius + max_width + 1.0,
        )
        anchor_radius = max(min(base_anchors), minimum_anchor)
        collision_radius = max(1.0, anchor_radius - max_width)
    else:
        base_anchors = [_outer_base_anchor(label, outer_reserved_radius_px, total_length) for label in ordered]
        anchor_radius = max(max(base_anchors), density_radius + 1.0)
        collision_radius = anchor_radius

    pushed_out = False
    for _ in range(_CONVERGENCE_LIMIT):
        angles = _pack_angles(ordered, total_length, collision_radius, spacing_px)
        placed = [
            _make_placed_label(
                label,
                total_length=total_length,
                placed_angle_rad=angle,
                text_anchor_radius_px=anchor_radius,
            )
            for label, angle in zip(ordered, angles, strict=True)
        ]
        # The leaders' room at the packed angles (OV-3841): inner text moves inward to it while it may, else the
        # radius must grow; outer text moves outward to it.
        room = [
            _leader_room_text_radius(label, angle, is_inner=is_inner)
            for label, angle in zip(ordered, angles, strict=True)
        ]
        if is_inner:
            limit = min(min(base_anchors), min(room))
            target = max(limit, minimum_anchor)
            if not pushed_out and anchor_radius > target + _LEADER_ROOM_TOLERANCE_PX:
                collision_radius -= anchor_radius - target
                anchor_radius = target
                continue
            required_growth = max(0.0, anchor_radius - limit)
        elif max(room) > anchor_radius + _LEADER_ROOM_TOLERANCE_PX:
            collision_radius += max(room) - anchor_radius
            anchor_radius = max(room)
            continue
        collisions = _text_collision_pairs(placed, spacing_px)
        leader_text_collisions = _leader_text_collision_count(placed)
        leader_crossings = leader_crossing_count(placed)
        order_violations = _order_violation_count(placed)
        awaiting_inner_preflight = is_inner and required_growth > RADIAL_GROWTH_TOLERANCE_PX
        if not collisions and (
            awaiting_inner_preflight
            or (not leader_text_collisions and not leader_crossings and not order_violations)
        ):
            return placed, required_growth
        radial_deficit = max(
            (
                _pair_clearance_px(placed[left], placed[right], spacing_px)
                for left, right in collisions
            ),
            default=max(
                _RADIAL_CLEARANCE_PX,
                max(float(label.get("height_px", 0.0)) for label in placed)
                + max(0.0, float(spacing_px)),
            ),
        )
        anchor_radius += radial_deficit
        collision_radius += radial_deficit
        pushed_out = True
    raise LabelLayoutError(
        "radial label packing did not converge after deterministic radius expansion; "
        f"side={'inner' if is_inner else 'outer'}, labels={len(ordered)}"
    )


def place_radial_labels(
    outer_labels: Sequence[CircularLabel],
    inner_labels: Sequence[CircularLabel],
    *,
    total_length: int,
    spacing_px: float,
    outer_reserved_radius_px: float = 0.0,
    inner_reserved_outer_radius_px: float = 0.0,
) -> CircularLabelLayoutResult:
    """Place all selected external labels without count or text-length cutoffs."""
    if total_length <= 0:
        raise ValidationError("radial circular labels require a positive genome length")
    outer, outer_growth = _place_side(
        outer_labels,
        total_length=total_length,
        spacing_px=spacing_px,
        is_inner=False,
        outer_reserved_radius_px=outer_reserved_radius_px,
        inner_reserved_outer_radius_px=inner_reserved_outer_radius_px,
    )
    inner, inner_growth = _place_side(
        inner_labels,
        total_length=total_length,
        spacing_px=spacing_px,
        is_inner=True,
        outer_reserved_radius_px=outer_reserved_radius_px,
        inner_reserved_outer_radius_px=inner_reserved_outer_radius_px,
    )
    labels = sorted(
        outer + inner,
        key=lambda label: (float(label["placed_angle_deg"]), bool(label.get("is_inner", False))),
    )
    text_collisions = len(_text_collision_pairs(labels, spacing_px))
    leader_text_collisions = _leader_text_collision_count(labels)
    leader_crossings = leader_crossing_count(labels)
    order_violations = _order_violation_count(outer) + _order_violation_count(inner)
    # The assembly grows the canvas while growth is required and refuses a layout
    # that still collides once it stops (LP-D49 ``RADIAL-JOIN``).
    return CircularLabelLayoutResult(
        labels=tuple(labels),
        content_bounds=_content_bounds(labels),
        required_radius_growth_px=max(outer_growth, inner_growth),
        text_collision_count=text_collisions,
        leader_text_collision_count=leader_text_collisions,
        leader_crossing_count=leader_crossings,
        order_violation_count=order_violations,
    )


__all__ = ["RADIAL_GROWTH_TOLERANCE_PX", "place_radial_labels", "radial_text_orientation"]
