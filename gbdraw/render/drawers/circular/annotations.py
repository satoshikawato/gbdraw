"""Circular annotation SVG geometry."""

from __future__ import annotations

import math
import sys
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field, replace
from functools import partial
from typing import Literal, NamedTuple, TypeAlias

from svgwrite.base import BaseElement

from gbdraw.svg.elements import ClipPath, Group, Path, Text, TextPath

from gbdraw.annotations import (
    AnnotationTrackParams,
    MarkPart,
    ResolvedAnnotationTrack,
    annotation_label_font_size,
    apply_feature_underlay_dom_attributes,
    effective_annotation_style,
    is_auto_feature_underlay,
    mark_parts,
)
from gbdraw.annotations.label_search import nearest_fewest_conflicts
from gbdraw.annotations.layout import (
    LANE_EXTENT_PX,
    LabelExtent,
    PlacedRegionAnnotation,
    lane_count_of,
    lane_edge_px,
    lane_width_px,
)
from gbdraw.annotations.models import RegionAnnotationStyle, ResolvedRegionAnnotation, annotation_solid_fill
from gbdraw.core.text import GlyphOutlines, TextRun, italic_markup_runs, text_runs_glyph_outlines
from gbdraw.layout.circular import RadialBand
from gbdraw.layout.circular_ink import CircularInk, CircularInkPolygon, CircularInkRing, stroke_pad_px
from gbdraw.layout.spatial import Aabb, union_aabbs
from gbdraw.layout.text_geometry import (
    Point,
    annular_sector_parts,
    convex_polygons_intersect,
    glyphs_across_circle,
    glyphs_half_turn_on_circle,
    glyphs_meet_on_circle,
    grown_convex_polygon,
    segment_intersects_polygon,
)
from gbdraw.render.drawers.annotation_marks import (
    ANNOTATION_ID_ATTRIBUTE,
    ANNOTATION_SET_ID_ATTRIBUTE,
    MarkReach,
    annotation_dom_id,
    cap_size,
    gradient_fill_group,
    gradient_steps,
    label_ink_offsets,
    mark_reach,
    past_span_edge,
)
from gbdraw.render.group_content import GroupContent
from gbdraw.render.patterns import hatch_pattern, hatch_pattern_id, hatch_pattern_paint
from gbdraw.svg.ids import stable_svg_id
from gbdraw.svg.circular_arcs import circle_arc, circle_arc_to, circle_point, position_angle
from gbdraw.svg.path_data import path_number, path_point
from gbdraw.svg.text_path import write_text_runs

Polygon = tuple[Point, ...]


def _point(position: float, total_length: int, radius: float) -> tuple[float, float]:
    angle = 2.0 * math.pi * (float(position) / max(1, total_length)) - 0.5 * math.pi
    return radius * math.cos(angle), radius * math.sin(angle)


def _span_deg(start: float, end: float, total_length: int) -> float:
    return 360.0 * max(0.0, end - start) / max(1, total_length)


def _arc_path(start: float, end: float, total_length: int, radius: float) -> str:
    return circle_arc(radius, position_angle(start, max(1, total_length)), _span_deg(start, end, total_length))


def _band_path(start: float, end: float, total_length: int, inner: float, outer: float) -> str:
    start_angle = position_angle(start, max(1, total_length))
    span = _span_deg(start, end, total_length)
    return (
        f"{circle_arc(outer, start_angle, span)} "
        f"L {path_point(*circle_point(inner, start_angle + span))} {circle_arc_to(inner, start_angle + span, -span)} Z"
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

        span = 2.0 * self.half_span
        if self.flipped:
            return circle_arc(self.radius, self.angle + self.half_span, -span)
        return circle_arc(self.radius, self.angle - self.half_span, span)


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
    far from it those letters come. ``runs`` are the styled runs it draws
    (``italic_markup_runs`` of the region's label). A ``cut`` label is curved
    text longer than its circle, drawn as the start of its runs and an ellipsis;
    one cut so far that no letter is left keeps no letter (``keeps_no_letter``).
    """

    annotation_id: str
    runs: tuple[TextRun, ...]
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

    @property
    def keeps_no_letter(self) -> bool:
        """Whether the cut left only the ellipsis (``_cut_to``): such a label is not drawn (RGN-CUT-EMPTY-OMIT)."""

        return self.cut and self.runs == _CUT_ONLY


@dataclass(frozen=True)
class CircularAnnotationLabels:
    """The labels of one record's annotation track (``circular_annotation_labels``): those drawn, and those left
    out because no letter of their text fits their circle (``without_letters``, RGN-CUT-EMPTY-OMIT, OV-2804): the
    labels whose cut keeps no letter, and those the track leaves out because their cut would keep none where their
    text would put their mark (``ResolvedAnnotationTrack.without_letters``, OV-2805), each at its midpoint. A label
    left out is not drawn and takes no place, room, obstacle or canvas bounds; the label fit warning counts it,
    with its font size. ``overlaps`` counts what the drawn labels still overlap where they are drawn: each label or
    leader line of the obstacles, and each label of the track placed before it, that a label's letters meet counts
    once (REGION-WARN, OV-3015); 0 when the labels were placed without obstacles."""

    drawn: tuple[CircularAnnotationLabel, ...]
    without_letters: tuple[CircularAnnotationLabel, ...]
    overlaps: int = 0


def _record_lanes(
    track: ResolvedAnnotationTrack,
    *,
    record_index: int,
    inner_radius_px: float,
    outer_radius_px: float,
    params: AnnotationTrackParams,
) -> tuple[list[PlacedRegionAnnotation], float]:
    """The record's placements and the width of one lane: the lanes share the band less its padding and the slot's
    label room (``lane_width_px``)."""

    placements = [item for item in track.placements if item.annotation.record_index == record_index]
    lane_count = lane_count_of(placements)
    width = lane_width_px(
        outer_radius_px - inner_radius_px, padding_px=params.padding_px, room_px=track.lane_room_px, lane_count=lane_count
    )
    return placements, LANE_EXTENT_PX if track.clipped else width


def _lane_radii(
    track: ResolvedAnnotationTrack,
    placement: PlacedRegionAnnotation,
    lane_width: float,
    *,
    side: str,
    inner_radius_px: float,
    outer_radius_px: float,
    params: AnnotationTrackParams,
) -> tuple[float, float]:
    """The inner and outer radius of a placement's lane: ``inside`` slots stack their lanes inward, each past the
    label room of the gaps before it (``ResolvedAnnotationTrack.room_before``). A lane stacked past the centre (a
    clipped slot of set width) stops at it: a negative radius would draw the mark through the centre, on the far
    side of the circle (OV-2666)."""

    # A circle draws its record on one row.
    room = track.room_before(placement.annotation.record_index, 0, placement.lane)
    if side == "inside":
        lane_outer = lane_edge_px(
            placement.lane, lane_width, padding_px=params.padding_px, room_before_px=room, from_px=outer_radius_px, toward=-1
        )
        return max(0.0, lane_outer - lane_width), max(0.0, lane_outer)
    lane_inner = lane_edge_px(
        placement.lane, lane_width, padding_px=params.padding_px, room_before_px=room, from_px=inner_radius_px
    )
    return lane_inner, lane_inner + lane_width


StrandBands: TypeAlias = Mapping[Literal["+", "-"], RadialBand]


@dataclass(frozen=True)
class AnnotationRing:
    """Where a circular annotation track is drawn: the slot's ring (``CircularResolvedSlot.draw_inner_radius_px`` and
    ``draw_outer_radius_px``), the side its lanes stack from, and each strand's lane of the feature track it covers
    (``CircularResolvedSlot.strand_bands_px``; none when it covers no track that separates the strands)."""

    inner_radius_px: float
    outer_radius_px: float
    side: str
    strand_bands_px: StrandBands = field(default_factory=dict)


def _strand_lanes(
    track: ResolvedAnnotationTrack, *, record_index: int, strand_bands_px: StrandBands, params: AnnotationTrackParams
) -> dict[Literal["+", "-"], tuple[RadialBand, float]]:
    """Each strand lane of the feature track a slot covers (``CircularResolvedSlot.strand_bands_px``), with the width
    of one of the record's lanes in it."""

    return {
        strand: (band, _record_lanes(
            track, record_index=record_index, inner_radius_px=band.inner_px, outer_radius_px=band.outer_px, params=params
        )[1])
        for strand, band in strand_bands_px.items()
    }


def _mark_lane(
    track: ResolvedAnnotationTrack,
    placement: PlacedRegionAnnotation,
    lane_width: float,
    strand_lanes: Mapping[Literal["+", "-"], tuple[RadialBand, float]],
    *,
    ring: AnnotationRing,
    params: AnnotationTrackParams,
) -> tuple[float, float, float]:
    """The inner and outer radius of a placement's lane and the lane's width: in its strand's lane when the slot
    covers a feature track that separates the strands, else in the slot's band."""

    inner_radius_px, outer_radius_px = ring.inner_radius_px, ring.outer_radius_px
    strand = placement.annotation.display_strand
    if strand is not None and strand in strand_lanes:
        band, lane_width = strand_lanes[strand]
        inner_radius_px, outer_radius_px = band.inner_px, band.outer_px
    lane_inner, lane_outer = _lane_radii(
        track, placement, lane_width, side=ring.side, inner_radius_px=inner_radius_px, outer_radius_px=outer_radius_px,
        params=params,
    )
    return lane_inner, lane_outer, lane_width


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
    met by several parts of a curved label counts once. ``add`` takes a label placed since, so a track's labels are
    boxed once each, not once per later label (OV-2664). ``conflicts(label, clearance=0.0)`` counts what the letters
    overlap, with no clearance kept (OV-3015)."""

    def __init__(self, obstacles: CircularLabelObstacles) -> None:
        self._obstacles = obstacles
        self._shapes: list[tuple[Aabb, list[tuple[Aabb, Polygon]]]] = []
        for shape in obstacles.shapes:
            self._add_shape(shape)
        self._segments = [
            (_box((segment.start, segment.end)), (segment.start, segment.end), segment.half_width)
            for segment in obstacles.segments
        ]

    def _add_shape(self, shape: Sequence[Polygon]) -> None:
        self._shapes.append((_box([point for part in shape for point in part]), [(_box(part), part) for part in shape]))

    def add(self, label: CircularAnnotationLabel) -> None:
        """Count ``label`` as an obstacle from now on, as ``CircularLabelObstacles.with_labels`` would."""

        self._add_shape(label.parts)

    def conflicts(self, label: CircularAnnotationLabel, clearance: float | None = None) -> int:
        # An obstacle whose box misses the box around all the grown letters misses each part's box: it is skipped
        # before the parts are paired with it.
        kept = self._obstacles.clearance if clearance is None else clearance
        near = [(_box(part), part) for part in (grown_convex_polygon(part, kept) for part in label.parts)]
        reach = union_aabbs(box for box, _ in near)
        if reach is None:
            return 0
        count = sum(
            1
            for box, parts in self._shapes
            if reach.intersects(box)
            and any(
                part_box.intersects(box) and obstacle_box.intersects(part_box) and convex_polygons_intersect(part, obstacle)
                for part_box, part in near
                for obstacle_box, obstacle in parts
            )
        )
        # The letters grown by each distinct segment clearance (one per leader width), made once per call.
        line_near = {
            distance: [(_box(part), part) for part in (grown_convex_polygon(part, distance) for part in label.parts)]
            for distance in {kept + item[2] for item in self._segments}
        }
        line_reach = {distance: union_aabbs(box for box, _ in parts) for distance, parts in line_near.items()}
        return count + sum(
            1
            for box, segment, half_width in self._segments
            if (reach_at := line_reach[kept + half_width]) is not None
            and reach_at.intersects(box)
            and any(
                part_box.intersects(box) and segment_intersects_polygon(segment, part, touching=False)
                for part_box, part in line_near[kept + half_width]
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
    join_origin: bool,
) -> CircularAnnotationLabel:
    """The label moved along its own region to the clear place nearest its default (OV-1173, the Linear rule).

    Its centre stays over one of the region's segments, sampled every pixel of
    arc at the label's radius; a place whose letters would reach past the
    slot's label room (``fits``) is not taken. With no clear place the label
    takes the place with the fewest conflicts, nearest its default.
    """

    length = max(1, record_length)
    centre = annotation.label_midpoint_bp(join_origin=join_origin)
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
) -> CircularAnnotationLabels:
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
    its region (``_along_region``) within the room the default places take, to
    a place where it keeps a letter. A label whose cut at its midpoint keeps no
    letter, or that the track leaves out for want of a letter, is left out
    (``CircularAnnotationLabels.without_letters``).
    """

    if not params.show_labels:
        return CircularAnnotationLabels((), ())
    placements, lane_width = _record_lanes(
        track, record_index=record_index, inner_radius_px=inner_radius_px, outer_radius_px=outer_radius_px, params=params
    )
    outward = side != "inside"
    entries: list[tuple[ResolvedRegionAnnotation, Callable[[float], CircularAnnotationLabel]]] = []
    for placement in placements:
        annotation = placement.annotation
        # A label the track leaves out (no room, or cut by its clip) is not drawn, placed or given room; one it leaves
        # out for want of a letter is counted below.
        if not (track.draws_label(annotation) or track.leaves_out_letters(annotation)):
            continue
        style = effective_annotation_style(annotation, params)
        lane_inner, lane_outer = _lane_radii(
            track,
            placement,
            lane_width,
            side=side,
            inner_radius_px=inner_radius_px,
            outer_radius_px=outer_radius_px,
            params=params,
        )
        reach = circular_annotation_mark_reach(annotation, style, lane_width)
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
                    runs=italic_markup_runs(annotation.label),
                    record_length=record_length,
                    mark_edge=mark_edge,
                    outward=outward,
                    gap=float(style.label_offset),
                    font_size=size,
                    color=style.label_color,
                    orientation=style.label_orientation,
                    font_family=font_family,
                    rings=rings,
                ),
            )
        )
    defaults = [label_at(annotation.label_midpoint_bp(join_origin=track.join_origin)) for annotation, label_at in entries]
    letterless = [label.keeps_no_letter or track.leaves_out_letters(entry[0]) for entry, label in zip(entries, defaults)]
    without_letters = tuple(label for label, left_out in zip(defaults, letterless) if left_out)
    kept = [(entry, label) for entry, label, left_out in zip(entries, defaults, letterless) if not left_out]
    if obstacles is None:
        return CircularAnnotationLabels(tuple(label for _entry, label in kept), without_letters)
    room = max((_reach(label, inner_radius_px, outer_radius_px, outward) for _entry, label in kept), default=0.0)
    placed: list[CircularAnnotationLabel] = []
    index = _ObstacleIndex(obstacles)
    overlaps = 0
    for (annotation, label_at), default in kept:
        label = _along_region(
            annotation,
            default,
            label_at,
            record_length=record_length,
            fits=lambda label: not label.keeps_no_letter
            and _reach(label, inner_radius_px, outer_radius_px, outward) <= room + 1e-9,
            conflicts=index.conflicts,
            join_origin=track.join_origin,
        )
        # What the placed letters still overlap, without the clearance the places keep (OV-3015).
        overlaps += index.conflicts(label, clearance=0.0)
        placed.append(label)
        index.add(label)
    return CircularAnnotationLabels(tuple(placed), without_letters, overlaps)


def circular_annotation_lanes_meeting(
    track: ResolvedAnnotationTrack,
    params: AnnotationTrackParams,
    extent_px: float,
    *,
    font_family: str,
    record_length: int,
    lanes_from_px: float,
    outward: bool,
    clearance_px: float,
) -> frozenset[tuple[int, int, int]]:
    """The lanes (record index, row, lane) in which ``circular_annotation_labels`` draws two of ``track``'s labels over
    each other, the lanes taking ``extent_px`` from ``lanes_from_px``, the radius they stack from (away from the
    centre, or toward it when not ``outward``), with the track's own labels as the only obstacles, kept
    ``clearance_px`` apart as where the drawing places them (``CircularLabelObstacles.clearance``): each label moved
    along its region where it comes that near one placed before it. The lane packer's measure (``LabelSizing.meet``)."""

    inner, outer = (lanes_from_px, lanes_from_px + extent_px) if outward else (lanes_from_px - extent_px, lanes_from_px)
    met: set[tuple[int, int, int]] = set()
    for record_index in sorted({item.annotation.record_index for item in track.placements}):
        lanes = {item.annotation.id: item.lane for item in track.placements if item.annotation.record_index == record_index}
        drawn = circular_annotation_labels(
            track, record_index=record_index, record_length=record_length, inner_radius_px=max(0.0, inner),
            outer_radius_px=outer, side="outside" if outward else "inside", font_family=font_family, params=params,
            obstacles=CircularLabelObstacles(clearance=clearance_px),
        ).drawn
        for index, first in enumerate(drawn):
            for second in drawn[index + 1 :]:
                lane = lanes[first.annotation_id]
                if lanes[second.annotation_id] == lane and any(
                    convex_polygons_intersect(one, other) for one in first.parts for other in second.parts
                ):
                    met.add((record_index, 0, lane))
    return frozenset(met)


def circular_annotation_labels_across_centre(
    track: ResolvedAnnotationTrack,
    params: AnnotationTrackParams,
    extent_px: float,
    *,
    font_family: str,
    record_length: int,
    lanes_from_px: float,
    clearance_px: float,
) -> frozenset[tuple[int, str]]:
    """The labels (record index, annotation id) of an inside track whose letters ``circular_annotation_labels`` draws
    over the centre of the circle and within their ``label_offset`` of another annotation's painted mark
    (``_mark_ink``), its lanes taking ``extent_px`` toward the centre from ``lanes_from_px`` and its labels placed as
    ``circular_annotation_lanes_meeting`` places them. Letters too wide for the circle their mark leaves them are
    centred on the centre of the circle and come back across it, where no room keeps them clear of the marks they meet
    (OV-2991). The lane packer's measure (``LabelSizing.across_centre``)."""

    ring = AnnotationRing(max(0.0, lanes_from_px - extent_px), lanes_from_px, "inside")
    met: set[tuple[int, str]] = set()
    for record_index in sorted({item.annotation.record_index for item in track.placements}):
        across = [
            label
            for label in circular_annotation_labels(
                track, record_index=record_index, record_length=record_length, inner_radius_px=ring.inner_radius_px,
                outer_radius_px=ring.outer_radius_px, side="inside", font_family=font_family, params=params,
                obstacles=CircularLabelObstacles(clearance=clearance_px),
            ).drawn
            if label.ink_radii[0] <= 0.0
        ]
        if not across:
            continue
        placements, lane_width = _record_lanes(
            track, record_index=record_index, inner_radius_px=ring.inner_radius_px,
            outer_radius_px=ring.outer_radius_px, params=params,
        )
        strand_lanes = _strand_lanes(track, record_index=record_index, strand_bands_px=ring.strand_bands_px, params=params)
        marks = {
            placement.annotation.id: [
                part
                for ink in _mark_ink(
                    track, placement, lane_width, strand_lanes, record_length=record_length, ring=ring, params=params
                )
                for part in _ink_parts(ink)
            ]
            for placement in placements
        }
        met.update(
            (record_index, label.annotation_id)
            for label in across
            if any(
                convex_polygons_intersect(grown_convex_polygon(part, label.label_offset), mark)
                for part in label.parts
                for identity, parts in marks.items()
                if identity != label.annotation_id
                for mark in parts
            )
        )
    return frozenset(met)


def _ink_parts(ink: CircularInk) -> tuple[Sequence[Point], ...]:
    """A mark's paint (``circular_annotation_track_ink``) as convex parts: a ring's arc cut into sector parts
    (``annular_sector_parts``), a polygon as it is."""

    if isinstance(ink, CircularInkPolygon):
        return (ink.points,)
    if not isinstance(ink, CircularInkRing):
        return ()
    span = 1.0 if ink.end_turn - ink.start_turn >= 1.0 else (ink.end_turn - ink.start_turn) % 1.0
    return annular_sector_parts(
        ink.band.inner_px, ink.band.outer_px, 360.0 * (ink.start_turn + 0.5 * span) - 90.0, 180.0 * span
    )


def _label_ink(font_family: str, font_size: float, orientation: str) -> LabelInk:
    """How a label's ink box is measured around its anchor (``_label_at``)."""

    return partial(
        label_ink_offsets,
        font_family=font_family,
        font_size=font_size,
        anchor="middle",
        baseline="alphabetic" if orientation in CURVED_LABEL_ORIENTATIONS else "central",
    )


def circular_annotation_mark_reach(
    annotation: ResolvedRegionAnnotation, style: RegionAnnotationStyle, lane_width: float
) -> MarkReach:
    """How far a Circular mark paints across its lane on its label's side (``MarkReach``), where its label stands
    ``label_offset`` beyond: a band fills its whole lane, and an arrow cap is drawn as a straight cap. The lane
    packer measures the same paint (``LabelSizing.reach``, ``circular_annotation_mark_past_lane``)."""

    return mark_reach(annotation, style, lane_width, fill_share=1.0, arrow_heads=False)


def circular_annotation_mark_past_lane(
    annotation: ResolvedRegionAnnotation, style: RegionAnnotationStyle, lane_width: float
) -> float:
    """How far a Circular mark paints past its lane's edge on its label's side (``circular_annotation_mark_reach``):
    the lane packer's ``MarkReachMeasure``."""

    return circular_annotation_mark_reach(annotation, style, lane_width).past_lane_edge


def circular_annotation_label_extent(
    annotation: ResolvedRegionAnnotation,
    style: RegionAnnotationStyle,
    font_size: float,
    place_bp: float,
    mark_offset_px: float,
    *,
    font_family: str,
    record_length: int,
    lanes_from_px: float,
    outward: bool,
) -> LabelExtent:
    """A region label's letters as ``circular_annotation_labels`` draws them at ``place_bp`` (``_label_at``),
    beside a mark whose edge lies ``mark_offset_px`` from ``lanes_from_px``, the radius the slot's lanes stack
    from: away from the centre, or toward it when not ``outward`` (an ``inside`` slot). The bp the letters reach
    along the circle on each side of the place, how far they reach across the lanes, and whether the label is drawn
    there (not when its cut keeps no letter, ``CircularAnnotationLabel.keeps_no_letter``): the lane packer's measure
    (``LabelSizing.extent``)."""

    edge = lanes_from_px + mark_offset_px if outward else lanes_from_px - mark_offset_px
    label = _label_at(
        place_bp,
        annotation_id="",
        runs=italic_markup_runs(annotation.label),
        record_length=record_length,
        mark_edge=max(1.0, edge),
        outward=outward,
        gap=float(style.label_offset),
        font_size=font_size,
        color="",
        orientation=style.label_orientation,
        font_family=font_family,
        rings=(),
    )
    centre = 360.0 * place_bp / max(1, record_length) - 90.0
    # Each corner's turn from the place, in degrees clockwise.
    turns = [
        (math.degrees(math.atan2(y, x)) - centre + 180.0) % 360.0 - 180.0 for part in label.parts for x, y in part
    ]
    bp_per_degree = record_length / 360.0
    near, far = label.ink_radii
    return LabelExtent(
        min(turns) * bp_per_degree, max(turns) * bp_per_degree, far - near, drawn=not label.keeps_no_letter
    )


# A curved label's text path stays just short of a whole circle (one SVG arc command); the text keeps one em clear of
# each end of it.
_MAX_HALF_PATH_DEG = 179.9
# What ends a curved label cut to fit its circle (OV-1210).
_CUT_MARK = "\u2026"
# A curved label cut to nothing but the ellipsis.
_CUT_ONLY: tuple[TextRun, ...] = ((_CUT_MARK, "normal"),)

LabelInk: TypeAlias = Callable[[Sequence[TextRun]], tuple[float, float, float, float]]


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
    runs: tuple[TextRun, ...],
    record_length: int,
    mark_edge: float,
    outward: bool,
    gap: float,
    font_size: float,
    color: str,
    orientation: str,
    font_family: str,
    rings: Sequence[RadialBand],
) -> CircularAnnotationLabel:
    """One label with its centre on the ray through ``position_bp``, ``gap`` beyond ``mark_edge`` (and any ring).

    Its ink box is measured around its anchor (``_label_ink``): on the central
    baseline for straight text, on the alphabetic baseline for curved text
    (``CURVED_LABEL_ORIENTATIONS``), whose glyphs are also measured one by one
    (``text_runs_glyph_outlines``).
    """

    angle = 360.0 * position_bp / max(1, record_length) - 90.0
    clear = partial(_cleared_radius, mark_edge=mark_edge, outward=outward, gap=gap, rings=rings)
    label = partial(
        CircularAnnotationLabel, annotation_id=annotation_id, font_size=font_size, color=color, label_offset=gap
    )
    ink = _label_ink(font_family, font_size, orientation)
    if orientation in CURVED_LABEL_ORIENTATIONS:
        glyphs = partial(text_runs_glyph_outlines, font_family=font_family, font_size=font_size)
        return _curved_label(label, runs, angle, ink, glyphs, clear, outward=outward, font_size=font_size)
    offsets = ink(runs)
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
        runs=runs,
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
    runs: tuple[TextRun, ...],
    angle: float,
    ink: LabelInk,
    glyphs: Callable[[tuple[TextRun, ...]], GlyphOutlines],
    clear: Callable[[Callable[[float], tuple[float, float]], Callable[[float], float]], float],
    *,
    outward: bool,
    font_size: float,
) -> CircularAnnotationLabel:
    """A label drawn along the circle through its anchor at ``angle``, upright (OV-1210).

    The text keeps one em clear of each end of a path just short of a whole
    circle, measured on the circle its letters come nearest the centre on, and
    no two of its letters' ink overlap where they stand (``glyphs``, one by
    one); a text that does not is cut: it keeps its start and ends with an
    ellipsis, centred on the label like a whole one, and the label is ``cut``.
    """

    # Text in the lower half reads anticlockwise with its tops toward the centre; elsewhere its tops face out.
    flipped = 0.0 < angle % 360.0 < 180.0

    def placed(shown: tuple[TextRun, ...]) -> tuple[float, float, float, float, float]:
        """The ink of ``shown`` along its baseline (left, right), how near the centre and how far from it its
        letters come, and its baseline radius, cleared of its mark (and rings)."""

        left, right, top, bottom = ink(shown)
        low, high = (top, bottom) if flipped else (-bottom, -top)
        # Each glyph stands upright and unbent on the circle, so its corners off the ray through its middle reach
        # farther from the centre than its ink box along the ray (OV-2979). Where no glyph is known, the ink box's
        # extent moves one for one with the baseline radius.
        across = glyphs_across_circle(glyphs(shown))

        def extent(radius: float) -> tuple[float, float]:
            return (radius + low, radius + high) if across is None else across.radii(radius, tops_in=flipped)

        def radius_for(target: float) -> float:
            reach = (
                target - (low if outward else high)
                if across is None
                else across.radius_reaching(target, tops_in=flipped, nearest=outward)
            )
            # The nearest ink stays off the centre.
            return max(reach, 1e-6 - low)

        radius = clear(extent, radius_for)
        return left, right, *extent(radius), radius

    def fits(shown: tuple[TextRun, ...]) -> bool:
        # The ink reaches at most this far each way from the label's middle along the circle its letters come
        # nearest the centre on, where they stand closest together (OV-2806): a circle about as small as the
        # letters, or letters pressed against the centre past their mark, leave no room for one. And no two letters'
        # ink overlaps: letters that stand apart on a line close up at their tops on the inner side of a circle.
        left, right, near, _far, radius = placed(shown)
        return max(-left, right) <= math.radians(_MAX_HALF_PATH_DEG) * near - font_size and not glyphs_meet_on_circle(
            glyphs(shown), radius, tops_in=flipped
        )

    shown = runs if fits(runs) else _cut_to(runs, fits)
    left, right, near, far, radius = placed(shown)
    # The letters' shape turns as far each way as their ink box along the baseline; nearer the centre than the
    # baseline it turns as far as their ink does where that is further: a glyph does not bend to the circle, so there
    # it spans a wider turn than its arc on the baseline (OV-2851). Where no glyph is known, the ink box is measured
    # at the circle the letters come nearest the centre on.
    half_turn = 0.5 * (right - left) / radius
    ink_turn = glyphs_half_turn_on_circle(glyphs(shown), radius, tops_in=flipped)
    inner_turn = 0.5 * (right - left) / near if ink_turn is None else ink_turn
    # Along the reading direction an offset is an arc length on the baseline; anticlockwise when flipped.
    sign = -1.0 if flipped else 1.0
    centre = angle + sign * math.degrees(0.5 * (left + right) / radius)
    theta = math.radians(angle)
    return label(
        runs=shown,
        x=radius * math.cos(theta),
        y=radius * math.sin(theta),
        rotation=0.0,
        ink_radii=(near, far),
        # Letters near the centre cover at most one turn (OV-2665: an inside lane that rooms carry past the centre).
        parts=(
            *annular_sector_parts(near, far, centre, min(180.0, math.degrees(half_turn))),
            *(
                annular_sector_parts(near, min(radius, far), centre, min(180.0, math.degrees(inner_turn)))
                if inner_turn > half_turn and near < radius
                else ()
            ),
        ),
        arc=CircularLabelArc(
            radius=radius,
            angle=angle,
            half_span=min(_MAX_HALF_PATH_DEG, math.degrees((max(-left, right) + font_size) / radius)),
            flipped=flipped,
        ),
        cut=shown != runs,
    )


def _cut_to(runs: tuple[TextRun, ...], fits: Callable[[tuple[TextRun, ...]], bool]) -> tuple[TextRun, ...]:
    """The longest start of ``runs``, ended with an ellipsis in the style of its last run, that ``fits`` its circle
    where its own letters stand (only the ellipsis when no start fits): a longer start reaches further along the
    circle and no less far across it, and holds the letters of a shorter one, so the starts that fit are the shorter
    ones."""

    def cut(length: int) -> tuple[TextRun, ...]:
        kept: list[TextRun] = []
        for text, style in runs:
            if length <= 0:
                break
            kept.append((text[:length], style))
            length -= len(text)
        while kept and not kept[-1][0].rstrip():
            kept.pop()
        if not kept:
            return _CUT_ONLY
        text, style = kept[-1]
        return (*kept[:-1], (text.rstrip() + _CUT_MARK, style))

    low, high = 0, sum(len(text) for text, _style in runs) - 1
    while low < high:
        middle = (low + high + 1) // 2
        if fits(cut(middle)):
            low = middle
        else:
            high = middle - 1
    return cut(low)


class CircularLabelRoom(NamedTuple):
    """What a slot's labels take at a band (``circular_annotation_label_room``): how far they reach beyond it on their
    side (``px``), and whether the band leaves out, for want of a letter, a label the track draws where its lanes were
    packed (``drops_labels``; OV-2848): the slot may not stand there, or moving it would free room by losing labels."""

    px: float
    drops_labels: bool


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
) -> CircularLabelRoom:
    """How far the record's labels reach beyond the slot's ``band`` on their side (outward, or inward on an
    ``inside`` slot), ``label_offset`` beyond the letters included, so the next ring keeps that gap (OV-1192);
    for a span overlay, ``band`` is the stretched span band and ``covers_slots`` is set. Also whether the band leaves
    out a label the track draws (``ResolvedAnnotationTrack.draws_label``) because no letter of it fits its circle
    there (``CircularLabelRoom.drops_labels``).

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
    return CircularLabelRoom(
        max((_reach(label, band.inner_px, band.outer_px, outward) for label in labels.drawn), default=0.0),
        any((record_index, label.annotation_id) not in track.without_letters for label in labels.without_letters),
    )


def _label_text(
    label: CircularAnnotationLabel, font_family: str, definitions: dict[str, BaseElement], item_id: str
) -> Text:
    """A label's ``<text>``, its runs written by ``write_text_runs``; a curved label's path goes into the track's
    ``definitions``."""

    paint = {"fill": label.color, "stroke": "none", "font_size": label.font_size, "font_family": font_family}
    if label.arc is None:
        text = write_text_runs(
            Text("", insert=(label.x, label.y), text_anchor="middle", dominant_baseline="central", **paint), label.runs
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
    text.add(write_text_runs(TextPath(path, text="", startOffset="50%"), label.runs))
    return text


def _band_fill(mark: str, style: RegionAnnotationStyle) -> str:
    """The fill a band or highlight mark is drawn with."""

    return hatch_pattern_paint(style.hatch) if style.hatch else annotation_solid_fill(mark, style)


def _band_paint(
    target: Group,
    definitions: dict[str, BaseElement],
    annotation: ResolvedRegionAnnotation,
    style: RegionAnnotationStyle,
    part: MarkPart,
    start: float,
    end: float,
    record_length: int,
    inner: float,
    outer: float,
) -> str:
    """The fill of one part of a band or highlight mark: its hatch pattern (added to ``definitions``), its solid
    fill, or ``none`` over the steps of its gradient (added to ``target`` before the mark's outline)."""

    if style.hatch:
        definitions.setdefault(hatch_pattern_id(style.hatch), hatch_pattern(style.hatch))
    if style.gradient is None:
        return _band_fill(annotation.mark, style)
    steps = gradient_steps(
        annotation, style.gradient, part, start, end, record_length=record_length,
        length_px=math.radians(_span_deg(start, end, record_length)) * 0.5 * (inner + outer),
    )
    target.add(gradient_fill_group(steps, style.fill_opacity, lambda step, color: Path(
        d=_band_path(step, end, record_length, inner, outer), fill=color, stroke="none",
    )))
    return "none"



def _piece_ends(start: float, end: float, part: MarkPart) -> tuple[tuple[float, bool], tuple[float, bool]]:
    """Each end of one drawn piece of a part (``_safe_segments``) and whether it is drawn: not where the part was
    split into pieces or cut."""

    return (start, start == part.start and part.capped_start), (end, end == part.end and part.capped_end)


def _turn(position: float, record_length: int) -> float:
    """A display position in turns of the circle; a part joined across the origin ends past one turn."""

    turn = position / record_length
    return turn - 1.0 if turn > 1.0 else turn


def _mark_ink(
    track: ResolvedAnnotationTrack,
    placement: PlacedRegionAnnotation,
    slot_lane_width: float,
    strand_lanes: Mapping[Literal["+", "-"], tuple[RadialBand, float]],
    *,
    record_length: int,
    ring: AnnotationRing,
    params: AnnotationTrackParams,
) -> list[CircularInk]:
    """What ``circular_annotation_track_content`` paints for one placement's mark: its band or highlight (its area
    when filled, else its outline), or its line with its end caps."""

    annotation = placement.annotation
    inks: list[CircularInk] = []
    auto_feature_underlay = is_auto_feature_underlay(annotation)
    style = effective_annotation_style(annotation, params)
    lane_inner, lane_outer, lane_width = _mark_lane(track, placement, slot_lane_width, strand_lanes, ring=ring, params=params)
    radius = 0.5 * (lane_inner + lane_outer)
    pad = 0.0 if auto_feature_underlay else stroke_pad_px(style.stroke, style.stroke_width)
    for part in mark_parts(annotation, join_origin=track.join_origin):
        for start, end in _safe_segments(part.start, part.end, record_length):
            start_turn, end_turn = _turn(start, record_length), _turn(end, record_length)
            ends = _piece_ends(start, end, part)
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
                        for position, capped in ends
                        if capped
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
                    for position, capped in ends
                    if capped
                )
    return inks


def circular_annotation_track_ink(
    track: ResolvedAnnotationTrack,
    *,
    record_index: int,
    record_length: int,
    ring: AnnotationRing,
    params: AnnotationTrackParams,
    labels: Sequence[CircularAnnotationLabel] = (),
) -> tuple[CircularInk, ...]:
    """What ``circular_annotation_track_content`` paints with the same arguments: each band or highlight (its
    area when filled, else its outline), each line with its end caps, and the labels' letters."""

    if record_length <= 0:
        return ()
    placements, slot_lane_width = _record_lanes(
        track, record_index=record_index, inner_radius_px=ring.inner_radius_px, outer_radius_px=ring.outer_radius_px,
        params=params,
    )
    strand_lanes = _strand_lanes(track, record_index=record_index, strand_bands_px=ring.strand_bands_px, params=params)
    labels_by_id = {label.annotation_id: label for label in labels}
    inks: list[CircularInk] = []
    for placement in placements:
        annotation = placement.annotation
        inks.extend(
            _mark_ink(track, placement, slot_lane_width, strand_lanes, record_length=record_length, ring=ring, params=params)
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
    ring: AnnotationRing,
    font_family: str,
    params: AnnotationTrackParams,
    labels: Sequence[CircularAnnotationLabel] = (),
    dom_group_id: str | None = None,
    semantic_slot_id: str | None = None,
) -> GroupContent:
    """Build an annotation track around the circle centre (render stage ``annotation_content``).

    ``labels`` are the track's labels (``circular_annotation_labels``), drawn
    with their annotations. A region with a strand is drawn in its strand's
    lane of ``ring.strand_bands_px`` (a slot that covers a feature track that
    separates the strands), any other in the slot's band. The hatch patterns the marks use, in order of first
    use, and the clip path of a clipped track are the content's definitions.
    """

    group_id = dom_group_id or f"gbdraw-annotation-track-{track.slot_id}-{record_index + 1}"
    group = Group(id=group_id)
    definitions: dict[str, BaseElement] = {}
    if semantic_slot_id:
        group.attribs["data-gbdraw-slot-id"] = str(semantic_slot_id)
        group.attribs["data-gbdraw-slot-renderer"] = "annotations"
    placements, slot_lane_width = _record_lanes(
        track, record_index=record_index, inner_radius_px=ring.inner_radius_px, outer_radius_px=ring.outer_radius_px,
        params=params,
    )
    strand_lanes = _strand_lanes(track, record_index=record_index, strand_bands_px=ring.strand_bands_px, params=params)
    labels_by_id = {label.annotation_id: label for label in labels}

    for placement in placements:
        annotation = placement.annotation
        auto_feature_underlay = is_auto_feature_underlay(annotation)
        style = effective_annotation_style(annotation, params)
        lane_inner, lane_outer, lane_width = _mark_lane(track, placement, slot_lane_width, strand_lanes, ring=ring, params=params)
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
                    ANNOTATION_ID_ATTRIBUTE: annotation.id,
                    ANNOTATION_SET_ID_ATTRIBUTE: track.set_id,
                    "data-gbdraw-annotation-track-id": track.slot_id,
                    "data-gbdraw-record-id": record_id,
                    "data-gbdraw-record-index": str(record_index),
                    "data-gbdraw-annotation-mark": annotation.mark,
                    "data-gbdraw-annotation-label": annotation.label,
                }
            )
        dash = ",".join(f"{value:g}" for value in style.stroke_dasharray) or None
        parts = mark_parts(annotation, join_origin=track.join_origin)
        part_count = sum(len(_safe_segments(part.start, part.end, record_length)) for part in parts)
        part_index = 0
        for part in parts:
            for start, end in _safe_segments(part.start, part.end, record_length):
                part_index += 1
                ends = _piece_ends(start, end, part)
                closed = all(capped for _position, capped in ends)
                if annotation.mark in {"band", "highlight"}:
                    fill = _band_paint(item_group, definitions, annotation, style, part, start, end,
                                       record_length, lane_inner, lane_outer)
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
                        stroke_linecap="round" if part.capped_start and part.capped_end else "butt",
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
                if annotation.mark in {"band", "highlight"} and not auto_feature_underlay and not closed:
                    path.attribs["stroke"] = "none"
                    outline = _arc_path(start, end, record_length, lane_inner) + " " + _arc_path(start, end, record_length, lane_outer)
                    for position, capped in ends:
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
                    for position in (position for position, capped in ends if capped):
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
        clip_ring = Path(
            d=f"{circle_arc(ring.outer_radius_px, 0.0, -360.0)} Z {circle_arc(ring.inner_radius_px, 0.0, 360.0)} Z"
        )
        clip_ring.attribs["clip-rule"] = "evenodd"
        clip_path = ClipPath(id=clip_id)
        clip_path.add(clip_ring)
        clip_paths.append(clip_path)
        group.attribs["clip-path"] = f"url(#{clip_id})"
    return GroupContent.of(group, (*definitions.values(), *clip_paths))


__all__ = [
    "CURVED_LABEL_ORIENTATIONS",
    "AnnotationRing",
    "CircularAnnotationLabel",
    "CircularAnnotationLabels",
    "CircularLabelArc",
    "CircularLabelObstacles",
    "CircularLabelRoom",
    "LeaderSegment",
    "circular_annotation_label_extent",
    "circular_annotation_mark_past_lane",
    "circular_annotation_mark_reach",
    "circular_annotation_label_room",
    "circular_annotation_labels",
    "circular_annotation_labels_across_centre",
    "circular_annotation_track_content",
    "circular_annotation_track_ink",
]
