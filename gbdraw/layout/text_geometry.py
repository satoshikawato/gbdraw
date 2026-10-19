"""Geometry helpers for text boxes that may be rotated around their SVG anchor."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Iterable, Literal, Sequence, TypeAlias, cast

TextAnchor: TypeAlias = Literal["start", "middle", "end"]
Point: TypeAlias = tuple[float, float]
AabbTuple: TypeAlias = tuple[float, float, float, float]

# The tolerance of every test here (the separating-axis test and its box broad phase).
GEOMETRY_EPSILON = 1e-9


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
    return cast(TextAnchor, anchor)  # membership checked above


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
    # ``unrotated`` always holds exactly four points.
    return cast(
        tuple[Point, Point, Point, Point],
        tuple(rotate_point(point, rotation_deg) for point in unrotated),
    )


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
        corners=cast(tuple[Point, Point, Point, Point], corners),  # four offsets in, four points out
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
            if axis_len <= GEOMETRY_EPSILON:
                continue
            axis = (-edge_y / axis_len, edge_x / axis_len)
            min_first, max_first = _project_polygon(first, axis)
            min_second, max_second = _project_polygon(second, axis)
            if touching:
                separated = max_first < min_second - GEOMETRY_EPSILON or max_second < min_first - GEOMETRY_EPSILON
            else:
                separated = max_first <= min_second + GEOMETRY_EPSILON or max_second <= min_first + GEOMETRY_EPSILON
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
            first_right < second_left - GEOMETRY_EPSILON
            or second_right < first_left - GEOMETRY_EPSILON
            or first_bottom < second_top - GEOMETRY_EPSILON
            or second_bottom < first_top - GEOMETRY_EPSILON
        )
    else:
        broad_phase_separate = (
            first_right <= second_left + GEOMETRY_EPSILON
            or second_right <= first_left + GEOMETRY_EPSILON
            or first_bottom <= second_top + GEOMETRY_EPSILON
            or second_bottom <= first_top + GEOMETRY_EPSILON
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
        if abs(cross) <= GEOMETRY_EPSILON:
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
    (ax, ay), (bx, by) = first
    (cx, cy), (dx, dy) = second
    # Segments whose boxes lie apart (beyond the touch tolerance) neither cross nor touch.
    if (
        (ax if ax > bx else bx) < (cx if cx < dx else dx) - GEOMETRY_EPSILON
        or (cx if cx > dx else dx) < (ax if ax < bx else bx) - GEOMETRY_EPSILON
        or (ay if ay > by else by) < (cy if cy < dy else dy) - GEOMETRY_EPSILON
        or (cy if cy > dy else dy) < (ay if ay < by else by) - GEOMETRY_EPSILON
    ):
        return False
    o1 = ((bx - ax) * (cy - ay)) - ((by - ay) * (cx - ax))
    o2 = ((bx - ax) * (dy - ay)) - ((by - ay) * (dx - ax))
    o3 = ((dx - cx) * (ay - cy)) - ((dy - cy) * (ax - cx))
    o4 = ((dx - cx) * (by - cy)) - ((dy - cy) * (bx - cx))
    if ((o1 > GEOMETRY_EPSILON and o2 < -GEOMETRY_EPSILON) or (o1 < -GEOMETRY_EPSILON and o2 > GEOMETRY_EPSILON)) and (
        (o3 > GEOMETRY_EPSILON and o4 < -GEOMETRY_EPSILON) or (o3 < -GEOMETRY_EPSILON and o4 > GEOMETRY_EPSILON)
    ):
        return True
    if not touching:
        return False

    def on_segment(p: Point, q: Point, r: Point) -> bool:
        return (
            min(p[0], r[0]) - GEOMETRY_EPSILON <= q[0] <= max(p[0], r[0]) + GEOMETRY_EPSILON
            and min(p[1], r[1]) - GEOMETRY_EPSILON <= q[1] <= max(p[1], r[1]) + GEOMETRY_EPSILON
        )

    a, b, c, d = first[0], first[1], second[0], second[1]
    return (
        (abs(o1) <= GEOMETRY_EPSILON and on_segment(a, c, b))
        or (abs(o2) <= GEOMETRY_EPSILON and on_segment(a, d, b))
        or (abs(o3) <= GEOMETRY_EPSILON and on_segment(c, a, d))
        or (abs(o4) <= GEOMETRY_EPSILON and on_segment(c, b, d))
    )


# A curved text's annular sector is counted in parts this wide or narrower (nearly flat).
SECTOR_PART_DEG = 5.0


def annular_sector_parts(
    inner: float, outer: float, centre_deg: float, half_span_deg: float
) -> tuple[tuple[Point, Point, Point, Point], ...]:
    """The annular sector between ``inner`` and ``outer`` radius around the origin, ``half_span_deg`` each way of
    ``centre_deg`` (clockwise from 3 o'clock, as in SVG), cut into convex parts at most ``SECTOR_PART_DEG`` wide.

    This is the shape of text drawn along a circle (an embedded feature label, a curved region label).
    """

    steps = max(1, math.ceil(2.0 * half_span_deg / SECTOR_PART_DEG))
    angles = [
        math.radians(centre_deg - half_span_deg + 2.0 * half_span_deg * index / steps)
        for index in range(steps + 1)
    ]
    return tuple(
        (
            (inner * math.cos(start), inner * math.sin(start)),
            (outer * math.cos(start), outer * math.sin(start)),
            (outer * math.cos(end), outer * math.sin(end)),
            (inner * math.cos(end), inner * math.sin(end)),
        )
        for start, end in zip(angles, angles[1:])
    )


def grown_convex_polygon(polygon: Sequence[Point], by: float) -> tuple[Point, ...]:
    """A convex polygon with every edge moved ``by`` outward, its corners mitred (a rectangle grows by ``by`` on
    every side; a sector part keeps ``by`` from every edge too).

    A polygon with fewer than three distinct corners grows as its bounding box.
    """

    points = [point for index, point in enumerate(polygon) if point != polygon[index - 1]]
    area = sum(x1 * y2 - x2 * y1 for (x1, y1), (x2, y2) in zip(points, points[1:] + points[:1]))
    if len(points) < 3 or abs(area) <= GEOMETRY_EPSILON:
        left, top, right, bottom = aabb_from_points(polygon)
        return ((left - by, top - by), (right + by, top - by), (right + by, bottom + by), (left - by, bottom + by))
    sign = 1.0 if area > 0.0 else -1.0
    normals: list[Point] = []
    for (x1, y1), (x2, y2) in zip(points, points[1:] + points[:1]):
        length = math.hypot(x2 - x1, y2 - y1)
        normals.append((sign * (y2 - y1) / length, -sign * (x2 - x1) / length))
    grown: list[Point] = []
    for index, (x, y) in enumerate(points):
        (px, py), (nx, ny) = normals[index - 1], normals[index]
        # The corner moves along the bisector of its two edges' normals, so each edge moves ``by``.
        scale = by / (1.0 + px * nx + py * ny)
        grown.append((x + scale * (px + nx), y + scale * (py + ny)))
    return tuple(grown)


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
    "GEOMETRY_EPSILON",
    "OrientedTextBox",
    "Point",
    "SECTOR_PART_DEG",
    "TextAnchor",
    "aabb_from_points",
    "anchor_x_bounds",
    "annular_sector_parts",
    "convex_polygons_intersect",
    "grown_convex_polygon",
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
