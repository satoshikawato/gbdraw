"""Geometry helpers for text boxes that may be rotated around their SVG anchor."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable, Literal, Sequence, TypeAlias

TextAnchor: TypeAlias = Literal["start", "middle", "end"]
Point: TypeAlias = tuple[float, float]
AabbTuple: TypeAlias = tuple[float, float, float, float]

_EPSILON = 1e-9


@dataclass(frozen=True)
class OrientedTextBox:
    """A measured text rectangle and its axis-aligned broad-phase bounds."""

    anchor_x: float
    anchor_y: float
    width_px: float
    height_px: float
    rotation_deg: float
    text_anchor: TextAnchor
    corners: tuple[Point, Point, Point, Point]
    aabb: AabbTuple


def normalize_text_anchor(value: object) -> TextAnchor:
    anchor = str(value or "middle").strip().lower()
    if anchor not in {"start", "middle", "end"}:
        raise ValueError(f"unsupported text anchor: {value!r}")
    return anchor  # type: ignore[return-value]


def anchor_x_bounds(width_px: float, text_anchor: object) -> tuple[float, float]:
    """Return the unrotated horizontal extent relative to an SVG text anchor."""
    width = max(0.0, float(width_px))
    anchor = normalize_text_anchor(text_anchor)
    if anchor == "start":
        return 0.0, width
    if anchor == "end":
        return -width, 0.0
    return -0.5 * width, 0.5 * width


def rotate_point(point: Point, rotation_deg: float) -> Point:
    """Rotate a point around the origin in SVG's clockwise-positive coordinates."""
    radians = math.radians(float(rotation_deg))
    sin_theta = math.sin(radians)
    cos_theta = math.cos(radians)
    x_coord, y_coord = point
    return (
        (x_coord * cos_theta) - (y_coord * sin_theta),
        (x_coord * sin_theta) + (y_coord * cos_theta),
    )


def text_box_corner_offsets(
    width_px: float,
    height_px: float,
    rotation_deg: float,
    text_anchor: object = "middle",
    *,
    padding_px: float = 0.0,
) -> tuple[Point, Point, Point, Point]:
    """Return clockwise rectangle corners relative to the text anchor.

    Text height is centered on the anchor. Callers should use a matching
    ``dominant-baseline="middle"`` when rendering plain SVG text.
    """
    padding = max(0.0, float(padding_px))
    x_min, x_max = anchor_x_bounds(width_px, text_anchor)
    half_height = 0.5 * max(0.0, float(height_px))
    unrotated = (
        (x_min - padding, -half_height - padding),
        (x_max + padding, -half_height - padding),
        (x_max + padding, half_height + padding),
        (x_min - padding, half_height + padding),
    )
    return tuple(rotate_point(point, rotation_deg) for point in unrotated)  # type: ignore[return-value]


def translate_points(points: Iterable[Point], x_offset: float, y_offset: float) -> tuple[Point, ...]:
    return tuple((x + float(x_offset), y + float(y_offset)) for x, y in points)


def aabb_from_points(points: Sequence[Point]) -> AabbTuple:
    if not points:
        raise ValueError("at least one point is required")
    x_values = [float(point[0]) for point in points]
    y_values = [float(point[1]) for point in points]
    return min(x_values), min(y_values), max(x_values), max(y_values)


def oriented_text_box(
    anchor_x: float,
    anchor_y: float,
    width_px: float,
    height_px: float,
    rotation_deg: float = 0.0,
    text_anchor: object = "middle",
    *,
    padding_px: float = 0.0,
) -> OrientedTextBox:
    """Build an immutable oriented box for measured single-line text."""
    anchor = normalize_text_anchor(text_anchor)
    corners = translate_points(
        text_box_corner_offsets(
            width_px,
            height_px,
            rotation_deg,
            anchor,
            padding_px=padding_px,
        ),
        anchor_x,
        anchor_y,
    )
    return OrientedTextBox(
        anchor_x=float(anchor_x),
        anchor_y=float(anchor_y),
        width_px=max(0.0, float(width_px)),
        height_px=max(0.0, float(height_px)),
        rotation_deg=float(rotation_deg),
        text_anchor=anchor,
        corners=corners,  # type: ignore[arg-type]
        aabb=aabb_from_points(corners),
    )


def _project_polygon(points: Sequence[Point], axis: Point) -> tuple[float, float]:
    axis_x, axis_y = axis
    projections = [(x * axis_x) + (y * axis_y) for x, y in points]
    return min(projections), max(projections)


def convex_polygons_intersect(
    first: Sequence[Point],
    second: Sequence[Point],
    *,
    touching: bool = False,
) -> bool:
    """Return whether two convex polygons intersect using the separating-axis theorem."""
    if len(first) < 3 or len(second) < 3:
        return False
    for polygon in (first, second):
        for index, point in enumerate(polygon):
            next_point = polygon[(index + 1) % len(polygon)]
            edge_x = next_point[0] - point[0]
            edge_y = next_point[1] - point[1]
            axis_len = math.hypot(edge_x, edge_y)
            if axis_len <= _EPSILON:
                continue
            axis = (-edge_y / axis_len, edge_x / axis_len)
            min_first, max_first = _project_polygon(first, axis)
            min_second, max_second = _project_polygon(second, axis)
            if touching:
                separated = max_first < min_second - _EPSILON or max_second < min_first - _EPSILON
            else:
                separated = max_first <= min_second + _EPSILON or max_second <= min_first + _EPSILON
            if separated:
                return False
    return True


def oriented_boxes_intersect(
    first: OrientedTextBox,
    second: OrientedTextBox,
    *,
    touching: bool = False,
) -> bool:
    """Return whether two oriented text boxes overlap."""
    first_left, first_top, first_right, first_bottom = first.aabb
    second_left, second_top, second_right, second_bottom = second.aabb
    if touching:
        broad_phase_separate = (
            first_right < second_left - _EPSILON
            or second_right < first_left - _EPSILON
            or first_bottom < second_top - _EPSILON
            or second_bottom < first_top - _EPSILON
        )
    else:
        broad_phase_separate = (
            first_right <= second_left + _EPSILON
            or second_right <= first_left + _EPSILON
            or first_bottom <= second_top + _EPSILON
            or second_bottom <= first_top + _EPSILON
        )
    if broad_phase_separate:
        return False
    return convex_polygons_intersect(first.corners, second.corners, touching=touching)


def point_in_convex_polygon(point: Point, polygon: Sequence[Point], *, touching: bool = True) -> bool:
    """Return whether a point lies inside a consistently ordered convex polygon."""
    if len(polygon) < 3:
        return False
    signs: list[float] = []
    for index, first in enumerate(polygon):
        second = polygon[(index + 1) % len(polygon)]
        cross = ((second[0] - first[0]) * (point[1] - first[1])) - (
            (second[1] - first[1]) * (point[0] - first[0])
        )
        if abs(cross) <= _EPSILON:
            if not touching:
                return False
            continue
        signs.append(cross)
    return not signs or all(value > 0.0 for value in signs) or all(value < 0.0 for value in signs)


def segments_intersect(
    first: tuple[Point, Point],
    second: tuple[Point, Point],
    *,
    touching: bool = True,
) -> bool:
    """Return whether two closed segments intersect."""
    a, b = first
    c, d = second

    def orientation(p: Point, q: Point, r: Point) -> float:
        return ((q[0] - p[0]) * (r[1] - p[1])) - ((q[1] - p[1]) * (r[0] - p[0]))

    def on_segment(p: Point, q: Point, r: Point) -> bool:
        return (
            min(p[0], r[0]) - _EPSILON <= q[0] <= max(p[0], r[0]) + _EPSILON
            and min(p[1], r[1]) - _EPSILON <= q[1] <= max(p[1], r[1]) + _EPSILON
        )

    o1 = orientation(a, b, c)
    o2 = orientation(a, b, d)
    o3 = orientation(c, d, a)
    o4 = orientation(c, d, b)
    proper = ((o1 > _EPSILON and o2 < -_EPSILON) or (o1 < -_EPSILON and o2 > _EPSILON)) and (
        (o3 > _EPSILON and o4 < -_EPSILON) or (o3 < -_EPSILON and o4 > _EPSILON)
    )
    if proper:
        return True
    if not touching:
        return False
    return (
        (abs(o1) <= _EPSILON and on_segment(a, c, b))
        or (abs(o2) <= _EPSILON and on_segment(a, d, b))
        or (abs(o3) <= _EPSILON and on_segment(c, a, d))
        or (abs(o4) <= _EPSILON and on_segment(c, b, d))
    )


def segment_intersects_polygon(
    segment: tuple[Point, Point],
    polygon: Sequence[Point],
    *,
    touching: bool = True,
) -> bool:
    """Return whether a segment crosses or lies inside a convex polygon."""
    if point_in_convex_polygon(segment[0], polygon, touching=touching):
        return True
    if point_in_convex_polygon(segment[1], polygon, touching=touching):
        return True
    return any(
        segments_intersect(
            segment,
            (polygon[index], polygon[(index + 1) % len(polygon)]),
            touching=touching,
        )
        for index in range(len(polygon))
    )


__all__ = [
    "AabbTuple",
    "OrientedTextBox",
    "Point",
    "TextAnchor",
    "aabb_from_points",
    "anchor_x_bounds",
    "convex_polygons_intersect",
    "normalize_text_anchor",
    "oriented_boxes_intersect",
    "oriented_text_box",
    "point_in_convex_polygon",
    "rotate_point",
    "segment_intersects_polygon",
    "segments_intersect",
    "text_box_corner_offsets",
    "translate_points",
]
