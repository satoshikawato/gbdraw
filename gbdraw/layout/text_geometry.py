"""Geometry helpers for text boxes that may be rotated around their SVG anchor."""

from __future__ import annotations

import functools
import math
from dataclasses import dataclass
from typing import TYPE_CHECKING, Iterable, Literal, NamedTuple, Sequence, TypeAlias, cast

import numpy as np
from numpy.typing import NDArray

from ..svg.circular_arcs import circle_point

if TYPE_CHECKING:
    # Read at run time only by its attributes: an import here would join the module cycle through gbdraw.core.text.
    from gbdraw.core.text import GlyphOutlines

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


def reads_upside_down(angle_deg: float) -> bool:
    """Whether text written clockwise along a circle with its middle at ``angle_deg`` (the SVG frame: 0 at 3
    o'clock, positive clockwise) reads upside down: its middle lies on the lower half, between 3 and 9 o'clock.
    There the record band names turn a half turn and the scalar axes' range labels run anticlockwise, so both
    read upright."""

    return 0.0 < angle_deg % 360.0 < 180.0


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
    along_text_px: float = 0.0,
) -> tuple[Point, Point, Point, Point]:
    """Return clockwise rectangle corners relative to the text anchor.

    Text height is centered on the anchor. Callers should use a matching
    ``dominant-baseline="middle"`` when rendering plain SVG text. ``padding_px``
    widens the box on every side, ``along_text_px`` only at both ends of the text
    (in its reading direction, whatever the rotation).
    """
    padding = max(0.0, float(padding_px))
    along = padding + max(0.0, float(along_text_px))
    x_min, x_max = anchor_x_bounds(width_px, text_anchor)
    half_height = 0.5 * max(0.0, float(height_px))
    unrotated = (
        (x_min - along, -half_height - padding),
        (x_max + along, -half_height - padding),
        (x_max + along, half_height + padding),
        (x_min - along, half_height + padding),
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
    edges = [
        (circle_point(inner, angle), circle_point(outer, angle))
        for angle in (centre_deg - half_span_deg + 2.0 * half_span_deg * index / steps for index in range(steps + 1))
    ]
    return tuple(
        (start_in, start_out, end_out, end_in) for (start_in, start_out), (end_in, end_out) in zip(edges, edges[1:])
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


def _cross(a: NDArray[np.float64], b: NDArray[np.float64]) -> NDArray[np.float64]:
    return a[..., 0] * b[..., 1] - a[..., 1] * b[..., 0]


class _PlacedGlyph(NamedTuple):
    """One glyph standing along a circle (``GlyphsAlongCircle``): its outline's edges (``starts`` to ``ends``) and
    the first point of each of its contours, turned about the centre."""

    starts: NDArray[np.float64]
    ends: NDArray[np.float64]
    firsts: NDArray[np.float64]


def _holds(glyph: _PlacedGlyph, points: NDArray[np.float64]) -> bool:
    """Whether the even-odd fill of ``glyph`` holds any of ``points``."""

    starts, ends = glyph.starts, glyph.ends
    x, y = points[:, None, 0], points[:, None, 1]
    spans = (starts[None, :, 1] > y) != (ends[None, :, 1] > y)
    with np.errstate(divide="ignore", invalid="ignore"):
        cross_x = starts[None, :, 0] + (y - starts[None, :, 1]) * (ends[None, :, 0] - starts[None, :, 0]) / (
            ends[None, :, 1] - starts[None, :, 1]
        )
    return bool(((spans & (x < cross_x)).sum(axis=1) % 2 == 1).any())


def _outlines_overlap(first: _PlacedGlyph, second: _PlacedGlyph) -> bool:
    """Whether the even-odd fills of two glyphs overlap: an edge of one crosses an edge of the other, or, when none
    does, one fill holds a contour of the other."""

    a0, a1, b0, b1 = first.starts, first.ends, second.starts, second.ends
    if (
        a0[:, 0].min() > b0[:, 0].max()
        or b0[:, 0].min() > a0[:, 0].max()
        or a0[:, 1].min() > b0[:, 1].max()
        or b0[:, 1].min() > a0[:, 1].max()
    ):
        return False
    a_dir, b_dir = (a1 - a0)[:, None, :], (b1 - b0)[None, :, :]
    a_sides = _cross(a_dir, b0[None, :, :] - a0[:, None, :]) * _cross(a_dir, b1[None, :, :] - a0[:, None, :])
    b_sides = _cross(b_dir, a0[:, None, :] - b0[None, :, :]) * _cross(b_dir, a1[:, None, :] - b0[None, :, :])
    if bool(((a_sides < 0) & (b_sides < 0)).any()):
        return True
    return _holds(second, first.firsts) or _holds(first, second.firsts)


class _GlyphPairs:
    """Every pair of a text's glyphs, each once (``GlyphsAlongCircle.meet``): the pairs of the most glyphs any text has
    had so far, ordered by the later glyph, so the pairs of ``count`` glyphs are their start. One array pair, as large
    as the pairs of the longest text measured, which ``meet`` would otherwise build at every call."""

    def __init__(self) -> None:
        self._later, self._earlier = np.tril_indices(0, -1)

    def of(self, count: int) -> tuple[NDArray[np.intp], NDArray[np.intp]]:
        """The earlier and the later glyph of each pair of ``count`` glyphs."""

        pairs = count * (count - 1) // 2
        if pairs > len(self._later):
            self._later, self._earlier = np.tril_indices(count, -1)
        return self._earlier[:pairs], self._later[:pairs]


_GLYPH_PAIRS = _GlyphPairs()


@dataclass(frozen=True, eq=False)
class GlyphsAlongCircle:
    """One text's glyphs (``text_runs_glyph_outlines``) standing along a circle of ``radius``, its baseline, as an SVG
    ``textPath`` sets it: each glyph upright and unbent at the point of its advance's middle, turned to the circle
    there, its tops toward the centre when ``tops_in``. Every outline point across the ray through its glyph's middle
    (``x``), out along that ray (``out``), and its turn (radians) about the centre from that ray (``own_turns``), which
    whether two glyphs meet and how far the ink turns both read. Built by ``glyphs_along_circle``."""

    outlines: GlyphOutlines
    radius: float
    x: NDArray[np.float64]
    out: NDArray[np.float64]
    own_turns: NDArray[np.float64]

    def meet(self) -> bool:
        """Whether the ink of two of the glyphs overlaps. Letters that stand apart on a line close up at their tops on
        the inner side of a circle."""

        outlines, radius = self.outlines, self.radius
        if not outlines.glyph_count:
            return False
        turn = outlines.middles / radius
        angles = self.own_turns + turn
        first_points = list(outlines.glyph_starts[:-1])
        low, high = np.minimum.reduceat(angles, first_points), np.maximum.reduceat(angles, first_points)
        # Glyphs whose turns from the centre are apart, a whole turn either way included, cannot meet. A whole turn
        # on counts only where the turns span more than one: rounding keeps order, so where the least turn a whole
        # turn on does not come below the greatest, no glyph's does.
        least, most = float(low.min()), float(high.max())
        whole = 2 * math.pi
        shifts = [shift for shift, wraps in ((-whole, least < most - whole), (whole, least + whole < most)) if wraps]
        earlier, later = _GLYPH_PAIRS.of(len(low))
        low_earlier, high_earlier, low_later, high_later = low[earlier], high[earlier], low[later], high[later]
        turns_meet = (low_earlier < high_later) & (low_later < high_earlier)
        for shift in shifts:
            turns_meet |= (low_earlier < high_later + shift) & (low_later + shift < high_earlier)
        pairs = np.flatnonzero(turns_meet)
        if not len(pairs):
            return False
        # Turned clockwise about the centre by its glyph's turn, from the glyph standing at the top of the circle.
        cos, sin = np.cos(turn), np.sin(turn)
        x, out = self.x, self.out
        placed = np.column_stack((cos * x + sin * out, cos * out - sin * x))
        ends = placed[outlines.next_index]
        contour_starts = outlines.contour_starts

        def glyph(index: int) -> _PlacedGlyph:
            start, end = outlines.glyph_starts[index], outlines.glyph_starts[index + 1]
            firsts = contour_starts[(contour_starts >= start) & (contour_starts < end)]
            return _PlacedGlyph(placed[start:end], ends[start:end], placed[firsts])

        return any(_outlines_overlap(glyph(int(earlier[pair])), glyph(int(later[pair]))) for pair in pairs)

    def half_turn(self) -> float | None:
        """How far (radians) the ink reaches about the centre each way from the middle of its ink along the baseline,
        at most: its widest turn, which its letters take where they come nearest the centre, since a glyph does not
        bend to the circle (OV-2851). None when no glyph paints."""

        outlines = self.outlines
        if not outlines.glyph_count:
            return None
        along = outlines.middles + outlines.points[:, 0]
        middle = 0.5 * (float(along.min()) + float(along.max()))
        angles = self.own_turns + (outlines.middles - middle) / self.radius
        return max(-float(angles.min()), float(angles.max()))


def glyphs_along_circle(outlines: GlyphOutlines, radius: float, *, tops_in: bool) -> GlyphsAlongCircle:
    """``GlyphsAlongCircle`` of ``outlines`` on a baseline circle of ``radius``, tops toward the centre when
    ``tops_in``."""

    x = outlines.points[:, 0]
    out = radius + (-1.0 if tops_in else 1.0) * outlines.points[:, 1]
    return GlyphsAlongCircle(outlines, radius, x, out, np.arctan2(x, out))


class _EdgesNearBaseline(NamedTuple):
    """The outline's edges on one side of the baseline (``GlyphsAcrossCircle``), in the order of how near the centre
    they come along the ray: each edge's least ``u`` (``lowest``, ascending) and its index (``order``); and the point
    that comes nearest the centre along the ray (``x``, ``u``)."""

    lowest: NDArray[np.float64]
    order: NDArray[np.intp]
    x: float
    u: float


# Added to the bounds that leave out edges which cannot hold the nearest ink: rounding never leaves out one that can.
_EDGE_BOUND_MARGIN = 1e-9


def _hull(points: NDArray[np.float64]) -> NDArray[np.float64]:
    """The corners of the convex hull of ``points`` (monotone chain)."""

    ordered = sorted(set(map(tuple, points.tolist())))
    if len(ordered) <= 2:
        return np.asarray(ordered, dtype=float).reshape(-1, 2)

    def half(chain_points: Iterable[tuple[float, float]]) -> list[tuple[float, float]]:
        chain: list[tuple[float, float]] = []
        for point in chain_points:
            while len(chain) >= 2 and (chain[-1][0] - chain[-2][0]) * (point[1] - chain[-2][1]) - (
                chain[-1][1] - chain[-2][1]
            ) * (point[0] - chain[-2][0]) <= 0.0:
                chain.pop()
            chain.append(point)
        return chain[:-1]

    return np.asarray(half(ordered) + half(reversed(ordered)), dtype=float)


@dataclass(frozen=True, eq=False)
class GlyphsAcrossCircle:
    """How near the centre and how far from it the ink of one text (``text_runs_glyph_outlines``) comes when it stands
    along a circle, its baseline, as an SVG ``textPath`` sets it (``GlyphsAlongCircle``): each glyph upright and
    unbent at the point of its advance's middle, so a point ``x`` across the ray through that middle and ``y`` up
    from the baseline stands ``hypot(x, radius + u)`` from the centre, whatever its turn, ``u`` being ``y``, or ``-y``
    with its tops toward the centre. A glyph's corners off that ray stand farther from the centre than their height
    along it says (OV-2979). Built by ``glyphs_across_circle``: the outline's edges from ``x``, ``y`` by ``dx``,
    ``dy``, with the terms of ``radii`` that depend on neither the radius nor the side (``along``,
    ``inverse_length``); the corners of the points' convex hull (``hull``), where the ink reaches farthest; and the
    edges near the baseline on each side (``tops_out``, ``tops_in``), where it comes nearest."""

    x: NDArray[np.float64]
    y: NDArray[np.float64]
    dx: NDArray[np.float64]
    dy: NDArray[np.float64]
    along: NDArray[np.float64]
    inverse_length: NDArray[np.float64]
    hull: NDArray[np.float64]
    tops_out: _EdgesNearBaseline
    tops_in: _EdgesNearBaseline

    def radii(self, radius: float, *, tops_in: bool) -> tuple[float, float]:
        """How near the centre and how far from it the ink comes on a baseline circle of ``radius`` that keeps it off
        the centre, its tops toward the centre when ``tops_in``: the nearest point of an edge, the farthest corner.

        No point of an edge that stays farther from the baseline along the ray than the point nearest it is from the
        centre, past the radius, comes nearer the centre than that point: such edges are not measured."""

        side = -1.0 if tops_in else 1.0
        out = radius + side * self.hull[:, 1]
        far = math.sqrt(float((self.hull[:, 0] ** 2 + out * out).max()))
        edges = self.tops_in if tops_in else self.tops_out
        bound = math.hypot(edges.x, radius + edges.u) - radius + _EDGE_BOUND_MARGIN
        near_edges = edges.order[: int(np.searchsorted(edges.lowest, bound, side="right"))]
        dx, dy = self.dx[near_edges], self.dy[near_edges]
        t = np.minimum(np.maximum(-(self.along[near_edges] + side * radius * dy) * self.inverse_length[near_edges], 0.0), 1.0)
        x, out = self.x[near_edges] + t * dx, radius + side * (self.y[near_edges] + t * dy)
        return math.sqrt(float((x * x + out * out).min())), far

    def radius_reaching(self, target: float, *, tops_in: bool, nearest: bool) -> float:
        """The baseline radius at which the ink comes ``target`` from the centre at its nearest (``nearest``), or
        reaches ``target`` at its farthest; -inf when no baseline radius that keeps the ink off the centre is
        needed for it (``nearest``) or works (a corner stands across the ray farther than ``target``).

        A point stands ``target`` out on the baseline radius ``sqrt(target² - x²) - u``. Every corner of the hull
        within ``target`` puts the farthest reach at the least of those; the ink clears ``target`` from the greatest
        over its edges, which along an edge is greatest at an end or where the circle of ``target`` touches it. The
        point that holds the ink off stands at most ``target`` less that radius beyond the baseline, and the radius is
        at least the one the point nearest the baseline asks for: edges beyond are not measured."""

        side = -1.0 if tops_in else 1.0
        if not nearest:
            x, u = self.hull[:, 0], side * self.hull[:, 1]
            if bool((np.abs(x) > target).any()):
                return -math.inf
            return float((np.sqrt(target * target - x * x) - u).min())
        edges = self.tops_in if tops_in else self.tops_out
        near_edges = edges.order
        if abs(edges.x) < target:
            bound = target - (math.sqrt(target * target - edges.x * edges.x) - edges.u) + _EDGE_BOUND_MARGIN
            near_edges = near_edges[: int(np.searchsorted(edges.lowest, bound, side="right"))]
        x, u = self.x[near_edges], side * self.y[near_edges]
        dx, du = self.dx[near_edges], side * self.dy[near_edges]
        with np.errstate(divide="ignore", invalid="ignore"):
            touch = -np.sign(dx) * du * target / np.hypot(dx, du)
            t = (touch - x) / dx
        along = (dx != 0.0) & (t > 0.0) & (t < 1.0)
        x, u = np.concatenate((x, touch[along])), np.concatenate((u, u[along] + t[along] * du[along]))
        within = np.abs(x) < target
        if not within.any():
            return -math.inf
        return float((np.sqrt(target * target - x[within] ** 2) - u[within]).max())


def _edges_near_baseline(
    x: NDArray[np.float64], u: NDArray[np.float64], next_index: NDArray[np.intp]
) -> _EdgesNearBaseline:
    lowest = np.minimum(u, u[next_index])
    order = np.argsort(lowest, kind="stable")
    nearest = int(np.argmin(u))
    return _EdgesNearBaseline(lowest[order], order, float(x[nearest]), float(u[nearest]))


@functools.lru_cache(maxsize=2048)
def glyphs_across_circle(outlines: GlyphOutlines) -> GlyphsAcrossCircle | None:
    """The radial reach of one text's glyphs standing along a circle (``GlyphsAcrossCircle``), one per text like
    ``text_runs_glyph_outlines``. None when no glyph paints."""

    if not outlines.glyph_count:
        return None
    x, y = outlines.points[:, 0].copy(), outlines.points[:, 1].copy()
    following = outlines.next_index
    dx, dy = x[following] - x, y[following] - y
    lengths = dx * dx + dy * dy
    # A point repeated along a contour is an edge of no length: its nearest point is its start.
    inverse_length = np.divide(1.0, lengths, out=np.zeros_like(lengths), where=lengths > 0.0)
    return GlyphsAcrossCircle(
        x,
        y,
        dx,
        dy,
        x * dx + y * dy,
        inverse_length,
        _hull(outlines.points),
        _edges_near_baseline(x, y, following),
        _edges_near_baseline(x, -y, following),
    )


__all__ = [
    "AabbTuple",
    "GEOMETRY_EPSILON",
    "GlyphsAcrossCircle",
    "GlyphsAlongCircle",
    "OrientedTextBox",
    "Point",
    "SECTOR_PART_DEG",
    "TextAnchor",
    "aabb_from_points",
    "anchor_x_bounds",
    "annular_sector_parts",
    "convex_polygons_intersect",
    "glyphs_across_circle",
    "glyphs_along_circle",
    "grown_convex_polygon",
    "normalize_text_anchor",
    "oriented_boxes_intersect",
    "oriented_text_box",
    "point_in_convex_polygon",
    "reads_upside_down",
    "rotate_point",
    "segment_intersects_polygon",
    "segments_intersect",
    "text_box_corner_offsets",
    "translate_points",
]
