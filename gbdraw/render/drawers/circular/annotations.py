"""Circular annotation SVG geometry."""

from __future__ import annotations

import math
import sys
from collections.abc import Callable, Sequence
from dataclasses import dataclass, replace
from functools import partial
from typing import NamedTuple, TypeAlias

from svgwrite.base import BaseElement

from gbdraw.svg.elements import ClipPath, Group, Path, Text, TextPath

from gbdraw.annotations import (
    AnnotationTrackParams,
    ResolvedAnnotationTrack,
    annotation_label_font_size,
    apply_feature_underlay_dom_attributes,
    effective_annotation_style,
    is_auto_feature_underlay,
)
from gbdraw.annotations.label_search import nearest_fewest_conflicts
from gbdraw.annotations.layout import PlacedRegionAnnotation
from gbdraw.annotations.models import RegionAnnotationStyle, ResolvedRegionAnnotation
from gbdraw.layout.record_coordinates import DisplayFragment
from gbdraw.layout.circular import RadialBand
from gbdraw.layout.circular_ink import CircularInk, CircularInkPolygon, CircularInkRing, stroke_pad_px
from gbdraw.layout.spatial import Aabb
from gbdraw.layout.text_geometry import (
    Point,
    annular_sector_parts,
    convex_polygons_intersect,
    grown_convex_polygon,
    segment_intersects_polygon,
)
from gbdraw.render.drawers.annotation_marks import (
    annotation_dom_id,
    cap_size,
    label_ink_offsets,
    mark_reach,
    past_span_edge,
)
from gbdraw.render.group_content import GroupContent
from gbdraw.render.patterns import hatch_pattern, hatch_pattern_id, hatch_pattern_paint
from gbdraw.svg.ids import stable_svg_id
from gbdraw.svg.path_data import path_number

Polygon = tuple[Point, ...]


def _point(position: float, total_length: int, radius: float) -> tuple[float, float]:
    angle = 2.0 * math.pi * (float(position) / max(1, total_length)) - 0.5 * math.pi
    return radius * math.cos(angle), radius * math.sin(angle)


def _arc_path(start: float, end: float, total_length: int, radius: float) -> str:
    x1, y1 = _point(start, total_length, radius)
    x2, y2 = _point(end, total_length, radius)
    span = max(0.0, end - start)
    large = 1 if span > 0.5 * total_length else 0
    return f"M {path_number(x1)},{path_number(y1)} A {path_number(radius)},{path_number(radius)} 0 {large} 1 {path_number(x2)},{path_number(y2)}"


def _band_path(start: float, end: float, total_length: int, inner: float, outer: float) -> str:
    ox1, oy1 = _point(start, total_length, outer)
    ox2, oy2 = _point(end, total_length, outer)
    ix2, iy2 = _point(end, total_length, inner)
    ix1, iy1 = _point(start, total_length, inner)
    span = max(0.0, end - start)
    large = 1 if span > 0.5 * total_length else 0
    return (
        f"M {path_number(ox1)},{path_number(oy1)} A {path_number(outer)},{path_number(outer)} 0 {large} 1 {path_number(ox2)},{path_number(oy2)} "
        f"L {path_number(ix2)},{path_number(iy2)} A {path_number(inner)},{path_number(inner)} 0 {large} 0 {path_number(ix1)},{path_number(iy1)} Z"
    )


def _safe_segments(start: int, end: int, total_length: int):
    if end - start >= total_length:
        middle = start + total_length / 2.0
        return ((float(start), middle), (middle, float(end)))
    return ((float(start), float(end)),)


# The orientations whose labels are drawn along the circle (OV-1210); ``tangent`` keeps straight text on the tangent.
CURVED_LABEL_ORIENTATIONS = frozenset({"auto", "arc"})


@dataclass(frozen=True)
class CircularLabelArc:
    """The circle a curved label's alphabetic baseline follows, and the text path drawn on it.

    The path is centred on ``angle`` (degrees, clockwise from 3 o'clock) and
    spans ``half_span`` degrees each way; ``flipped`` text reads anticlockwise,
    so text in the lower half (past 3 o'clock to 9 o'clock) stays upright.
    """

    radius: float
    angle: float
    half_span: float
    flipped: bool

    @property
    def path(self) -> str:
        """The SVG path the label's ``<textPath>`` follows, in reading order."""

        start, end = self.angle - self.half_span, self.angle + self.half_span
        if self.flipped:
            start, end = end, start
        (x1, y1), (x2, y2) = (
            (self.radius * math.cos(math.radians(value)), self.radius * math.sin(math.radians(value)))
            for value in (start, end)
        )
        large = 1 if 2.0 * self.half_span > 180.0 else 0
        sweep = 0 if self.flipped else 1
        return f"M {path_number(x1)},{path_number(y1)} A {path_number(self.radius)},{path_number(self.radius)} 0 {large} {sweep} {path_number(x2)},{path_number(y2)}"


@dataclass(frozen=True)
class CircularAnnotationLabel:
    """One Circular annotation label as it is drawn, around the circle centre.

    Straight text (``arc`` is None) is drawn ``text-anchor=middle``,
    ``dominant-baseline=central`` at ``x``, ``y``, turned ``rotation`` degrees
    clockwise about that point. Curved text follows ``arc`` on its alphabetic
    baseline, its middle at ``x``, ``y``. ``parts`` are the convex shapes of
    its painted letters (the turned ink box, or an annular sector cut into
    nearly flat parts), which the room, the conflict count, later slots and
    the canvas bounds all read; ``ink_radii`` are how near the centre and how
    far from it those letters come. A ``cut`` label is curved text longer than
    its circle, drawn as the start of its text and an ellipsis (``text``).
    """

    annotation_id: str
    text: str
    x: float
    y: float
    rotation: float
    font_size: float
    color: str
    parts: tuple[Polygon, ...]
    ink_radii: tuple[float, float]
    label_offset: float
    arc: CircularLabelArc | None = None
    cut: bool = False

    @property
    def ink_bounds(self) -> Aabb:
        return _box([point for part in self.parts for point in part])


def _record_lanes(
    track: ResolvedAnnotationTrack,
    *,
    record_index: int,
    inner_radius_px: float,
    outer_radius_px: float,
    params: AnnotationTrackParams,
) -> tuple[list[PlacedRegionAnnotation], float]:
    """The record's placements and the width of one lane."""

    placements = [item for item in track.placements if item.annotation.record_index == record_index]
    lane_count = max((item.lane for item in placements), default=-1) + 1
    usable = max(1.0, outer_radius_px - inner_radius_px - 2.0 * params.padding_px)
    return placements, 14.0 if track.clipped else usable / max(1, lane_count)


def _lane_radii(
    lane: int,
    lane_width: float,
    *,
    side: str,
    inner_radius_px: float,
    outer_radius_px: float,
    params: AnnotationTrackParams,
) -> tuple[float, float]:
    """The inner and outer radius of a lane: ``inside`` slots stack their lanes inward."""

    if side == "inside":
        lane_outer = outer_radius_px - params.padding_px - lane * lane_width
        return lane_outer - lane_width, lane_outer
    lane_inner = inner_radius_px + params.padding_px + lane * lane_width
    return lane_inner, lane_inner + lane_width


def _label_rotation(orientation: str, angle: float) -> float:
    """The turn of a label at ``angle`` (degrees, clockwise from 3 o'clock): along the ring kept upright, along the
    ray, or none."""

    if orientation == "tangent":
        rotation = angle + 90.0
        return rotation + 180.0 if 90.0 < rotation % 360.0 < 270.0 else rotation
    return angle if orientation == "radial" else 0.0


def _ink_extent(
    radius: float, direction: tuple[float, float], offsets: tuple[float, float, float, float]
) -> tuple[float, float]:
    """How near the centre and how far from it a label's ink box comes with its anchor at ``radius``.

    ``direction`` is the unit ray in the text's own frame; the box excludes the centre.
    """

    left, right, top, bottom = offsets
    wx, wy = direction
    # The centre lies at -radius * direction from the anchor.
    ox, oy = -radius * wx, -radius * wy
    near = math.hypot(max(left - ox, 0.0, ox - right), max(top - oy, 0.0, oy - bottom))
    far = max(math.hypot(x - ox, y - oy) for x in (left, right) for y in (top, bottom))
    return near, far


def _label_radius(
    target: float,
    outward: bool,
    direction: tuple[float, float],
    offsets: tuple[float, float, float, float],
) -> float:
    """The anchor radius at which a label's ink comes nearest the centre at ``target`` (``outward``) or reaches
    farthest from it at ``target``."""

    left, right, top, bottom = offsets
    projections = [x * direction[0] + y * direction[1] for x in (left, right) for y in (top, bottom)]
    # The projection on the ray bounds the reach; each step closes the rest, as the reach moves almost one for one.
    radius = target - (min(projections) if outward else max(projections))
    for _ in range(32):
        near, far = _ink_extent(radius, direction, offsets)
        miss = target - (near if outward else far)
        if abs(miss) < 1e-9 or radius + miss <= 0.0:
            break
        radius += miss
    return max(radius, 0.0)


class LeaderSegment(NamedTuple):
    """One straight part of a feature label's leader, with half the width it is drawn at (its label's own)."""

    start: Point
    end: Point
    half_width: float


@dataclass(frozen=True)
class CircularLabelObstacles:
    """What Circular annotation labels keep clear of, around the circle centre (OV-1173).

    ``shapes`` are labels, each as its convex parts (the feature labels as
    the fit count sees them, and annotation labels placed before); a label
    met counts once, whichever of its parts the letters meet. ``segments`` are
    leader lines. Letters keep ``clearance`` from a shape, and ``clearance``
    plus the segment's half width from a segment.
    """

    shapes: tuple[tuple[Polygon, ...], ...] = ()
    segments: tuple[LeaderSegment, ...] = ()
    clearance: float = 0.0

    def with_labels(self, labels: Sequence[CircularAnnotationLabel]) -> CircularLabelObstacles:
        return replace(self, shapes=(*self.shapes, *(label.parts for label in labels)))


def _box(points: Sequence[Point]) -> Aabb:
    return Aabb(min(x for x, _ in points), min(y for _, y in points), max(x for x, _ in points), max(y for _, y in points))


class _ObstacleIndex:
    """Counts what a label's letters meet: the polygons within ``clearance``, the segments within theirs; an obstacle
    met by several parts of a curved label counts once."""

    def __init__(self, obstacles: CircularLabelObstacles) -> None:
        self._obstacles = obstacles
        self._shapes = [
            (_box([point for part in shape for point in part]), [(_box(part), part) for part in shape])
            for shape in obstacles.shapes
        ]
        self._segments = [
            (_box((segment.start, segment.end)), (segment.start, segment.end), obstacles.clearance + segment.half_width)
            for segment in obstacles.segments
        ]

    def conflicts(self, label: CircularAnnotationLabel) -> int:
        near = [(_box(part), part) for part in (grown_convex_polygon(part, self._obstacles.clearance) for part in label.parts)]
        count = sum(
            1
            for box, parts in self._shapes
            if any(
                part_box.intersects(box) and obstacle_box.intersects(part_box) and convex_polygons_intersect(part, obstacle)
                for part_box, part in near
                for obstacle_box, obstacle in parts
            )
        )
        # The letters grown by each distinct segment clearance (one per leader width), made once per call.
        line_near = {
            distance: [(_box(part), part) for part in (grown_convex_polygon(part, distance) for part in label.parts)]
            for distance in {item[2] for item in self._segments}
        }
        return count + sum(
            1
            for box, segment, distance in self._segments
            if any(
                part_box.intersects(box) and segment_intersects_polygon(segment, part, touching=False)
                for part_box, part in line_near[distance]
            )
        )


def _reach(label: CircularAnnotationLabel, inner_radius_px: float, outer_radius_px: float, outward: bool) -> float:
    """How far a label's letters, ``label_offset`` beyond them included, reach past its slot's band on its side."""

    near, far = label.ink_radii
    return far + label.label_offset - outer_radius_px if outward else inner_radius_px - near + label.label_offset


def _along_region(
    annotation: ResolvedRegionAnnotation,
    default: CircularAnnotationLabel,
    label_at: Callable[[float], CircularAnnotationLabel],
    *,
    record_length: int,
    fits: Callable[[CircularAnnotationLabel], bool],
    conflicts: Callable[[CircularAnnotationLabel], int],
) -> CircularAnnotationLabel:
    """The label moved along its own region to the clear place nearest its default (OV-1173, the Linear rule).

    Its centre stays over one of the region's segments, sampled every pixel of
    arc at the label's radius; a place whose letters would reach past the
    slot's label room (``fits``) is not taken. With no clear place the label
    takes the place with the fewest conflicts, nearest its default.
    """

    length = max(1, record_length)
    centre = float(annotation.geometry_midpoint_bp)
    step = max(1e-9, length / (2.0 * math.pi * max(1.0, math.hypot(default.x, default.y))))
    candidates: set[float] = set()
    for start, end in annotation.geometry_segments:
        steps = max(1, math.ceil((float(end) - float(start)) / step))
        for index in range(steps + 1):
            position = float(start) + (float(end) - float(start)) * index / steps
            # Positions are unwrapped around the default so that "nearest" is measured along the circle.
            candidates.add(centre + ((position - centre + 0.5 * length) % length - 0.5 * length))
    labels: dict[float, CircularAnnotationLabel] = {centre: default}

    def count(position: float) -> int:
        if position not in labels:
            labels[position] = label_at(position % length)
        label = labels[position]
        return conflicts(label) if fits(label) else sys.maxsize

    return labels[nearest_fewest_conflicts(centre, candidates, count)]


def circular_annotation_labels(
    track: ResolvedAnnotationTrack,
    *,
    record_index: int,
    record_length: int,
    inner_radius_px: float,
    outer_radius_px: float,
    side: str,
    font_family: str,
    params: AnnotationTrackParams,
    rings: Sequence[RadialBand] = (),
    obstacles: CircularLabelObstacles | None = None,
    covers_slots: bool = False,
) -> tuple[CircularAnnotationLabel, ...]:
    """The labels of one record's annotation track, around the centre of the circle.

    Each label's painted letters keep ``label_offset`` from the painted edge of
    its mark (stroke included), outward or, on an ``inside`` slot, inward, at
    the label's own turn (OV-849, the Linear rule of OV-844). A painted ring of
    ``rings`` (the axis circle) that the letters would meet, or come within
    ``label_offset`` of, is cleared the same way beyond its paint. A span that
    ``covers_slots`` measures the gap from its band's edge, or from the mark's
    paint where it passes that edge (OD-R1174 F1, the Linear span rule), so its
    letters stay off every ring it covers.

    A label sits at its region's midpoint. With ``obstacles``, a label whose
    letters meet one, or a label of this track placed before it, moves along
    its region (``_along_region``) within the room the default places take.
    """

    if not params.show_labels:
        return ()
    placements, lane_width = _record_lanes(
        track, record_index=record_index, inner_radius_px=inner_radius_px, outer_radius_px=outer_radius_px, params=params
    )
    outward = side != "inside"
    entries: list[tuple[ResolvedRegionAnnotation, Callable[[float], CircularAnnotationLabel]]] = []
    for placement in placements:
        annotation = placement.annotation
        if not annotation.label:
            continue
        style = effective_annotation_style(annotation, params)
        lane_inner, lane_outer = _lane_radii(
            placement.lane,
            lane_width,
            side=side,
            inner_radius_px=inner_radius_px,
            outer_radius_px=outer_radius_px,
            params=params,
        )
        # A Circular band fills its whole lane; an arrow cap is drawn as a straight cap.
        reach = mark_reach(annotation, style, lane_width, fill_share=1.0, arrow_heads=False)
        if covers_slots:
            inset = outer_radius_px - lane_outer if outward else lane_inner - inner_radius_px
            beyond = past_span_edge(reach, inset)
            mark_edge = outer_radius_px + beyond if outward else inner_radius_px - beyond
        else:
            centre = 0.5 * (lane_inner + lane_outer)
            mark_edge = centre + reach.from_centre if outward else centre - reach.from_centre
        size = annotation_label_font_size(style, track.label_font_size)
        entries.append(
            (
                annotation,
                partial(
                    _label_at,
                    annotation_id=annotation.id,
                    text=annotation.label,
                    record_length=record_length,
                    mark_edge=mark_edge,
                    outward=outward,
                    gap=float(style.label_offset),
                    font_size=size,
                    color=style.label_color,
                    orientation=style.label_orientation,
                    ink=partial(
                        label_ink_offsets,
                        font_family=font_family,
                        font_size=size,
                        anchor="middle",
                        baseline="alphabetic" if style.label_orientation in CURVED_LABEL_ORIENTATIONS else "central",
                    ),
                    rings=rings,
                ),
            )
        )
    defaults = [label_at(float(annotation.geometry_midpoint_bp)) for annotation, label_at in entries]
    if obstacles is None:
        return tuple(defaults)
    room = max((_reach(label, inner_radius_px, outer_radius_px, outward) for label in defaults), default=0.0)
    placed: list[CircularAnnotationLabel] = []
    for (annotation, label_at), default in zip(entries, defaults):
        index = _ObstacleIndex(obstacles.with_labels(placed))
        placed.append(
            _along_region(
                annotation,
                default,
                label_at,
                record_length=record_length,
                fits=lambda label: _reach(label, inner_radius_px, outer_radius_px, outward) <= room + 1e-9,
                conflicts=index.conflicts,
            )
        )
    return tuple(placed)


# A curved label's text path stays just short of a whole circle (one SVG arc command); the text keeps one em clear of
# each end of it.
_MAX_HALF_PATH_DEG = 179.9
# What ends a curved label cut to fit its circle (OV-1210).
_CUT_MARK = "\u2026"

LabelInk: TypeAlias = Callable[[str], tuple[float, float, float, float]]


def _cleared_radius(
    extent: Callable[[float], tuple[float, float]],
    radius_for: Callable[[float], float],
    *,
    mark_edge: float,
    outward: bool,
    gap: float,
    rings: Sequence[RadialBand],
) -> float:
    """The anchor radius whose letters keep ``gap`` beyond ``mark_edge``, and beyond the paint of any ring of
    ``rings`` they would otherwise meet or come within ``gap`` of."""

    edge = mark_edge
    radius = radius_for(edge + gap if outward else edge - gap)
    for ring in sorted(rings, key=lambda band: band.inner_px if outward else -band.outer_px):
        near, far = extent(radius)
        if outward and edge < ring.outer_px and ring.inner_px < far + gap:
            edge = ring.outer_px
        elif not outward and ring.inner_px < edge and near - gap < ring.outer_px:
            edge = ring.inner_px
        else:
            continue
        radius = radius_for(edge + gap if outward else edge - gap)
    return radius


def _label_at(
    position_bp: float,
    *,
    annotation_id: str,
    text: str,
    record_length: int,
    mark_edge: float,
    outward: bool,
    gap: float,
    font_size: float,
    color: str,
    orientation: str,
    ink: LabelInk,
    rings: Sequence[RadialBand],
) -> CircularAnnotationLabel:
    """One label with its centre on the ray through ``position_bp``, ``gap`` beyond ``mark_edge`` (and any ring).

    ``ink`` measures a text's ink box around the label's anchor: on the
    central baseline for straight text, on the alphabetic baseline for curved
    text (``CURVED_LABEL_ORIENTATIONS``).
    """

    angle = 360.0 * position_bp / max(1, record_length) - 90.0
    clear = partial(_cleared_radius, mark_edge=mark_edge, outward=outward, gap=gap, rings=rings)
    label = partial(
        CircularAnnotationLabel, annotation_id=annotation_id, font_size=font_size, color=color, label_offset=gap
    )
    if orientation in CURVED_LABEL_ORIENTATIONS:
        return _curved_label(label, text, angle, ink, clear, outward=outward, font_size=font_size)
    offsets = ink(text)
    theta = math.radians(angle)
    rotation = _label_rotation(orientation, angle)
    turn = math.radians(rotation)
    # The ray's direction in the text's frame: the screen direction turned back by the label's rotation.
    direction = (math.cos(theta - turn), math.sin(theta - turn))
    radius = clear(
        partial(_ink_extent, direction=direction, offsets=offsets),
        partial(_label_radius, outward=outward, direction=direction, offsets=offsets),
    )
    x, y = radius * math.cos(theta), radius * math.sin(theta)
    left, right, top, bottom = offsets
    cos, sin = math.cos(turn), math.sin(turn)
    return label(
        text=text,
        x=x,
        y=y,
        rotation=rotation,
        ink_radii=_ink_extent(radius, direction, offsets),
        parts=(
            tuple(
                (x + cx * cos - cy * sin, y + cx * sin + cy * cos)
                for cx, cy in ((left, top), (right, top), (right, bottom), (left, bottom))
            ),
        ),
    )


def _curved_label(
    label: Callable[..., CircularAnnotationLabel],
    text: str,
    angle: float,
    ink: LabelInk,
    clear: Callable[[Callable[[float], tuple[float, float]], Callable[[float], float]], float],
    *,
    outward: bool,
    font_size: float,
) -> CircularAnnotationLabel:
    """A label drawn along the circle through its anchor at ``angle``, upright (OV-1210).

    The text keeps one em clear of each end of a path just short of a whole
    circle; a text longer than that is cut: it keeps its start and ends with an
    ellipsis, centred on the label like a whole one, and the label is ``cut``.
    """

    # Text in the lower half reads anticlockwise with its tops toward the centre; elsewhere its tops face out.
    flipped = 0.0 < angle % 360.0 < 180.0
    shown = text
    while True:
        left, right, top, bottom = ink(shown)
        low, high = (top, bottom) if flipped else (-bottom, -top)

        def extent(radius: float, low: float = low, high: float = high) -> tuple[float, float]:
            return radius + low, radius + high

        def radius_for(target: float, low: float = low, high: float = high) -> float:
            # The letters' radial extent moves one for one with the baseline radius; the nearest ink stays off the
            # centre.
            return max(target - (low if outward else high), 1e-6 - low)

        radius = clear(extent, radius_for)
        # How far the ink may reach each way from the label's middle, along its baseline.
        room = math.radians(_MAX_HALF_PATH_DEG) * radius - font_size
        if max(-left, right) <= room or shown == _CUT_MARK:
            break
        shown = _cut_to(text, room, ink)
    near, far = extent(radius)
    # Along the reading direction an offset is an arc length on the baseline; anticlockwise when flipped.
    sign = -1.0 if flipped else 1.0
    centre = angle + sign * math.degrees(0.5 * (left + right) / radius)
    theta = math.radians(angle)
    return label(
        text=shown,
        x=radius * math.cos(theta),
        y=radius * math.sin(theta),
        rotation=0.0,
        ink_radii=(near, far),
        parts=annular_sector_parts(near, far, centre, math.degrees(0.5 * (right - left) / radius)),
        arc=CircularLabelArc(
            radius=radius,
            angle=angle,
            half_span=min(_MAX_HALF_PATH_DEG, math.degrees((max(-left, right) + font_size) / radius)),
            flipped=flipped,
        ),
        cut=shown != text,
    )


def _cut_to(text: str, room: float, ink: LabelInk) -> str:
    """The longest start of ``text``, ended with an ellipsis, whose ink reaches at most ``room`` each way from its
    middle (only the ellipsis when no start fits)."""

    def cut(length: int) -> str:
        return text[:length].rstrip() + _CUT_MARK

    def fits(length: int) -> bool:
        left, right, _top, _bottom = ink(cut(length))
        return max(-left, right) <= room

    low, high = 0, len(text) - 1
    while low < high:
        middle = (low + high + 1) // 2
        if fits(middle):
            low = middle
        else:
            high = middle - 1
    return cut(low)


def circular_annotation_label_room(
    track: ResolvedAnnotationTrack,
    *,
    record_index: int,
    record_length: int,
    band: RadialBand,
    side: str,
    font_family: str,
    params: AnnotationTrackParams,
    rings: Sequence[RadialBand] = (),
    covers_slots: bool = False,
) -> float:
    """How far the record's labels reach beyond the slot's ``band`` on their side (outward, or inward on an
    ``inside`` slot), ``label_offset`` beyond the letters included, so the next ring keeps that gap (OV-1192);
    for a span overlay, ``band`` is the stretched span band and ``covers_slots`` is set.

    Labels moved along their regions stay within this room (``circular_annotation_labels``)."""

    outward = side != "inside"
    labels = circular_annotation_labels(
        track,
        record_index=record_index,
        record_length=record_length,
        inner_radius_px=band.inner_px,
        outer_radius_px=band.outer_px,
        side=side,
        font_family=font_family,
        params=params,
        rings=rings,
        covers_slots=covers_slots,
    )
    return max((_reach(label, band.inner_px, band.outer_px, outward) for label in labels), default=0.0)


def _label_text(
    label: CircularAnnotationLabel, font_family: str, definitions: dict[str, BaseElement], item_id: str
) -> Text:
    """A label's ``<text>``; a curved label's path goes into the track's ``definitions``."""

    paint = {"fill": label.color, "stroke": "none", "font_size": label.font_size, "font_family": font_family}
    if label.arc is None:
        text = Text(
            label.text, insert=(label.x, label.y), text_anchor="middle", dominant_baseline="central", **paint
        )
        if label.rotation:
            text.rotate(label.rotation, center=(label.x, label.y))
        return text
    # The annotation's DOM id carries its record, so records merged onto one canvas keep distinct paths: a merge that
    # renamed a colliding id would not reach the copied <textPath>, which writes its href from its own path object.
    path_id = stable_svg_id("circular_annotation_label_path", item_id)
    path = Path(id=path_id, d=label.arc.path, fill="none", stroke="none")
    definitions[path_id] = path
    text = Text("", text_anchor="middle", **paint)
    text.add(TextPath(path, text=label.text, startOffset="50%"))
    return text


def _band_fill(mark: str, style: RegionAnnotationStyle) -> str:
    """The fill a band or highlight mark is drawn with."""

    return hatch_pattern_paint(style.hatch) if style.hatch else (
        style.fill or ("#94a3b8" if mark == "highlight" else "none")
    )


def _cap_positions(
    start: float, end: float, raw_start: float, raw_end: float, fragment: DisplayFragment | None
) -> tuple[float, ...]:
    """The ends of one drawn part of a bracket that get an end cap: not where the part was split or cut."""

    return tuple(
        position
        for position in (start, end)
        if not (
            fragment is not None
            and (
                (position == start and (start != raw_start or fragment.artificial_start))
                or (position == end and (end != raw_end or fragment.artificial_end))
            )
        )
    )


def circular_annotation_track_ink(
    track: ResolvedAnnotationTrack,
    *,
    record_index: int,
    record_length: int,
    inner_radius_px: float,
    outer_radius_px: float,
    side: str,
    params: AnnotationTrackParams,
    labels: Sequence[CircularAnnotationLabel] = (),
) -> tuple[CircularInk, ...]:
    """What ``circular_annotation_track_content`` paints with the same arguments: each band or highlight (its
    area when filled, else its outline), each line with its end caps, and the labels' letters."""

    if record_length <= 0:
        return ()
    placements, lane_width = _record_lanes(
        track, record_index=record_index, inner_radius_px=inner_radius_px, outer_radius_px=outer_radius_px, params=params
    )
    labels_by_id = {label.annotation_id: label for label in labels}
    inks: list[CircularInk] = []
    for placement in placements:
        annotation = placement.annotation
        auto_feature_underlay = is_auto_feature_underlay(annotation)
        style = effective_annotation_style(annotation, params)
        lane_inner, lane_outer = _lane_radii(
            placement.lane, lane_width, side=side, inner_radius_px=inner_radius_px, outer_radius_px=outer_radius_px,
            params=params,
        )
        radius = 0.5 * (lane_inner + lane_outer)
        pad = 0.0 if auto_feature_underlay else stroke_pad_px(style.stroke, style.stroke_width)
        for segment_index, (raw_start, raw_end) in enumerate(annotation.geometry_segments):
            fragment = annotation.display_parts[segment_index] if annotation.display_parts is not None else None
            for start, end in _safe_segments(raw_start, raw_end, record_length):
                start_turn, end_turn = start / record_length, end / record_length
                if annotation.mark in {"band", "highlight"}:
                    if _band_fill(annotation.mark, style) != "none":
                        inks.append(CircularInkRing(RadialBand(lane_inner - pad, lane_outer + pad), start_turn, end_turn))
                    elif pad > 0:
                        inks.extend(
                            CircularInkRing(RadialBand(edge - pad, edge + pad), start_turn, end_turn)
                            for edge in (lane_inner, lane_outer)
                        )
                        inks.extend(
                            CircularInkPolygon.segment(
                                _point(position, record_length, lane_inner),
                                _point(position, record_length, lane_outer),
                                pad,
                            )
                            for position in (start, end)
                        )
                    continue
                if pad <= 0:
                    continue
                # A round cap reaches half the stroke past each end.
                reach = pad / max(radius, 1e-9) / (2.0 * math.pi)
                inks.append(CircularInkRing(RadialBand(radius - pad, radius + pad), start_turn - reach, end_turn + reach))
                if annotation.mark == "bracket" and style.line_cap != "none":
                    cap = cap_size(lane_width)
                    inks.extend(
                        CircularInkPolygon.segment(
                            _point(position, record_length, radius - cap),
                            _point(position, record_length, radius + cap),
                            pad,
                        )
                        for position in _cap_positions(start, end, raw_start, raw_end, fragment)
                    )
        label = labels_by_id.get(annotation.id)
        if label is not None:
            inks.extend(CircularInkPolygon(tuple(part)) for part in label.parts if len(part) >= 3)
    return tuple(inks)


def circular_annotation_track_content(
    track: ResolvedAnnotationTrack,
    *,
    record_id: str,
    record_index: int,
    record_length: int,
    inner_radius_px: float,
    outer_radius_px: float,
    side: str,
    font_family: str,
    params: AnnotationTrackParams,
    labels: Sequence[CircularAnnotationLabel] = (),
    dom_group_id: str | None = None,
    semantic_slot_id: str | None = None,
) -> GroupContent:
    """Build an annotation track around the circle centre (render stage ``annotation_content``).

    ``labels`` are the track's labels (``circular_annotation_labels``), drawn
    with their annotations. The hatch patterns the marks use, in order of first
    use, and the clip path of a clipped track are the content's definitions.
    """

    group_id = dom_group_id or f"gbdraw-annotation-track-{track.slot_id}-{record_index + 1}"
    group = Group(id=group_id)
    definitions: dict[str, BaseElement] = {}
    if semantic_slot_id:
        group.attribs["data-gbdraw-slot-id"] = str(semantic_slot_id)
        group.attribs["data-gbdraw-slot-renderer"] = "annotations"
    placements, lane_width = _record_lanes(
        track, record_index=record_index, inner_radius_px=inner_radius_px, outer_radius_px=outer_radius_px, params=params
    )
    labels_by_id = {label.annotation_id: label for label in labels}

    for placement in placements:
        annotation = placement.annotation
        auto_feature_underlay = is_auto_feature_underlay(annotation)
        style = effective_annotation_style(annotation, params)
        lane_inner, lane_outer = _lane_radii(
            placement.lane,
            lane_width,
            side=side,
            inner_radius_px=inner_radius_px,
            outer_radius_px=outer_radius_px,
            params=params,
        )
        radius = 0.5 * (lane_inner + lane_outer)
        item_group = Group(
            id=annotation_dom_id(
                record_index=record_index,
                slot_id=track.slot_id,
                set_id=track.set_id,
                annotation_id=annotation.id,
            ),
        )
        if not auto_feature_underlay:
            item_group.attribs.update(
                {
                    "data-gbdraw-annotation-id": annotation.id,
                    "data-gbdraw-annotation-set-id": track.set_id,
                    "data-gbdraw-annotation-track-id": track.slot_id,
                    "data-gbdraw-record-id": record_id,
                    "data-gbdraw-record-index": str(record_index),
                    "data-gbdraw-annotation-mark": annotation.mark,
                    "data-gbdraw-annotation-label": annotation.label,
                }
            )
        dash = ",".join(f"{value:g}" for value in style.stroke_dasharray) or None
        part_count = sum(
            len(_safe_segments(start, end, record_length))
            for start, end in annotation.geometry_segments
        )
        part_index = 0
        for segment_index, (raw_start, raw_end) in enumerate(annotation.geometry_segments):
            fragment = annotation.display_parts[segment_index] if annotation.display_parts is not None else None
            for start, end in _safe_segments(raw_start, raw_end, record_length):
                part_index += 1
                if annotation.mark in {"band", "highlight"}:
                    if style.hatch:
                        definitions.setdefault(hatch_pattern_id(style.hatch), hatch_pattern(style.hatch))
                    fill = _band_fill(annotation.mark, style)
                    path = Path(
                        d=_band_path(start, end, record_length, lane_inner, lane_outer),
                        fill=fill,
                        fill_opacity=style.fill_opacity,
                        stroke="none" if auto_feature_underlay else style.stroke,
                        stroke_width=0 if auto_feature_underlay else style.stroke_width,
                    )
                else:
                    path = Path(
                        d=_arc_path(start, end, record_length, radius),
                        fill="none",
                        stroke=style.stroke,
                        stroke_width=style.stroke_width,
                        stroke_linecap=("butt" if fragment is not None and (fragment.artificial_start or fragment.artificial_end) else "round"),
                    )
                if dash:
                    path.attribs["stroke-dasharray"] = dash
                if auto_feature_underlay:
                    apply_feature_underlay_dom_attributes(
                        path,
                        annotation,
                        record_id=record_id,
                        record_index=record_index,
                        part_index=part_index,
                        part_count=part_count,
                        include_stable_id=False,
                    )
                item_group.add(path)
                if annotation.mark in {"band", "highlight"} and not auto_feature_underlay and fragment is not None and (fragment.artificial_start or fragment.artificial_end):
                    path.attribs["stroke"] = "none"
                    outline = _arc_path(start, end, record_length, lane_inner) + " " + _arc_path(start, end, record_length, lane_outer)
                    for position, capped in ((start, not fragment.artificial_start), (end, not fragment.artificial_end)):
                        if capped:
                            x1, y1 = _point(position, record_length, lane_inner)
                            x2, y2 = _point(position, record_length, lane_outer)
                            outline += f" M {path_number(x1)},{path_number(y1)} L {path_number(x2)},{path_number(y2)}"
                    border = Path(d=outline, fill="none", stroke=style.stroke, stroke_width=style.stroke_width, stroke_linecap="butt")
                    if dash:
                        border.attribs["stroke-dasharray"] = dash
                    item_group.add(border)
                if annotation.mark == "bracket" and style.line_cap != "none":
                    cap = cap_size(lane_width)
                    for position in _cap_positions(start, end, raw_start, raw_end, fragment):
                        x1, y1 = _point(position, record_length, radius - cap)
                        x2, y2 = _point(position, record_length, radius + cap)
                        cap_path = Path(
                            d=f"M {path_number(x1)},{path_number(y1)} L {path_number(x2)},{path_number(y2)}",
                            fill="none",
                            stroke=style.stroke,
                            stroke_width=style.stroke_width,
                        )
                        item_group.add(cap_path)

        label = labels_by_id.get(annotation.id)
        if label is not None:
            item_group.add(_label_text(label, font_family, definitions, str(item_group.attribs["id"])))
        group.add(item_group)
    clip_paths: list[BaseElement] = []
    if track.clipped:
        group.attribs["data-gbdraw-annotation-clipped"] = "true"
        clip_id = (
            stable_svg_id("annotation_clip", "circular", group_id)
            if dom_group_id
            else f"gbdraw-annotation-clip-{track.slot_id}-{record_index + 1}"
        )
        ring = Path(
            d=(
                f"M {path_number(outer_radius_px)},0 A {path_number(outer_radius_px)},{path_number(outer_radius_px)} 0 1 0 "
                f"{path_number(-outer_radius_px)},0 A {path_number(outer_radius_px)},{path_number(outer_radius_px)} 0 1 0 "
                f"{path_number(outer_radius_px)},0 Z M {path_number(inner_radius_px)},0 "
                f"A {path_number(inner_radius_px)},{path_number(inner_radius_px)} 0 1 1 {path_number(-inner_radius_px)},0 "
                f"A {path_number(inner_radius_px)},{path_number(inner_radius_px)} 0 1 1 {path_number(inner_radius_px)},0 Z"
            )
        )
        ring.attribs["clip-rule"] = "evenodd"
        clip_path = ClipPath(id=clip_id)
        clip_path.add(ring)
        clip_paths.append(clip_path)
        group.attribs["clip-path"] = f"url(#{clip_id})"
    return GroupContent.of(group, (*definitions.values(), *clip_paths))


__all__ = [
    "CURVED_LABEL_ORIENTATIONS",
    "CircularAnnotationLabel",
    "CircularLabelArc",
    "CircularLabelObstacles",
    "LeaderSegment",
    "circular_annotation_label_room",
    "circular_annotation_labels",
    "circular_annotation_track_content",
    "circular_annotation_track_ink",
]
