"""Circular annotation SVG geometry."""

from __future__ import annotations

import math
import sys
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field, replace
from functools import partial
from typing import Literal, NamedTuple, TypeAlias, TypeVar

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
from gbdraw.layout.spatial import Aabb, BoxedLetters, BoxedPart, boxed_parts_meet, union_aabbs
from gbdraw.layout.text_geometry import (
    Point,
    annular_sector_parts,
    convex_polygons_intersect,
    GlyphsAlongCircle,
    glyphs_across_circle,
    glyphs_along_circle,
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
from gbdraw.render.drawers.label_backing import (
    BackedRegionLabel,
    SectorBacking,
    TurnedBacking,
    region_label_front,
    sector_backing,
    turned_backing,
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
    far from it those letters come, and ``height`` how tall they stand across
    their baseline (the ink box's top to bottom for straight text, ``ink_radii``
    for curved text). ``runs`` are the styled runs it draws
    (``italic_markup_runs`` of the region's label). A ``cut`` label is curved
    text longer than its circle, drawn as the start of its runs and an ellipsis;
    one cut so far that no letter is left keeps no letter (``keeps_no_letter``).
    A ``backed`` label is drawn in front of the leader lines behind its letters
    on its ``backing_shape`` (``CIR-RGN-FRONT``, as Linear's ``RGN-FRONT``;
    OV-3156): Linear's rounded box turned with straight text, an annular
    sector for curved text (``label_backing``); ``backing`` is its convex parts.
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
    height: float
    arc: CircularLabelArc | None = None
    cut: bool = False
    backing_shape: TurnedBacking | SectorBacking | None = None

    @property
    def backed(self) -> bool:
        return self.backing_shape is not None

    @property
    def backing(self) -> tuple[Polygon, ...]:
        return () if self.backing_shape is None else self.backing_shape.parts

    @property
    def footprint(self) -> tuple[Polygon, ...]:
        """What the drawn label covers, which later labels keep clear of: its letters, and its backing when it has
        one (as Linear's ``LinearAnnotationLabel.footprint``)."""

        return (*self.parts, *self.backing)

    @property
    def ink_bounds(self) -> Aabb:
        return _box([point for part in self.parts for point in part])

    @property
    def keeps_no_letter(self) -> bool:
        """Whether the cut left only the ellipsis (``_cut_to``): such a label is not drawn (RGN-CUT-EMPTY-OMIT)."""

        return self.cut and self.runs == _CUT_ONLY


class RegionLabelPlace(NamedTuple):
    """Where an annotation's label is drawn along its region: the annotation's id and the place in bp, unwrapped
    around the region's label midpoint (``_along_region``)."""

    annotation_id: str
    place_bp: float


@dataclass(frozen=True)
class CircularAnnotationLabels:
    """The labels of one record's annotation track (``circular_annotation_labels``): those drawn, and those left
    out because no letter of their text fits their circle (``without_letters``, RGN-CUT-EMPTY-OMIT, OV-2804): the
    labels whose cut keeps no letter, and those the track leaves out because their cut would keep none where their
    text would put their mark (``ResolvedAnnotationTrack.without_letters``, OV-2805), each at its midpoint. A label
    left out is not drawn and takes no place, room, obstacle or canvas bounds; the label fit warning counts it,
    with its font size. With obstacles, each drawn label also has what it still overlaps where it is drawn
    (``overlapped``): the labels or leader lines of the obstacles, then the labels of the track placed before it, by
    their index, that its letters meet with no clearance kept (REGION-WARN, OV-3015); and its letters and their
    height grown by the obstacles' clearance (``letters``), which the feature labels keep clear when they make way
    (CIR-RGN-YIELD). Both are keyed by annotation id. ``over_feature_letters`` are the places of the drawn labels whose
    letters meet a feature label's letters (OV-3017)."""

    drawn: tuple[CircularAnnotationLabel, ...]
    without_letters: tuple[CircularAnnotationLabel, ...]
    overlapped: Mapping[str, frozenset[int]] = field(default_factory=dict)
    letters: Mapping[str, BoxedLetters] = field(default_factory=dict)
    over_feature_letters: tuple[RegionLabelPlace, ...] = ()

    @property
    def overlaps(self) -> int:
        """What the region label warning counts: each obstacle or earlier label a drawn label overlaps, once."""

        return sum(len(met) for met in self.overlapped.values())


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


class CircularObstacleLabel(NamedTuple):
    """One label that annotation labels keep clear of: its convex ``parts`` and its ``leader`` line's segments."""

    parts: tuple[Polygon, ...]
    leader: tuple[LeaderSegment, ...] = ()
    # A feature label, which makes way for an annotation label over its letters (CIR-RGN-YIELD); annotation labels
    # do not.
    feature: bool = False


@dataclass(frozen=True)
class CircularLabelObstacles:
    """What Circular annotation labels keep clear of, around the circle centre (OV-1173).

    ``labels`` are the feature labels as the fit count sees them, each with
    its leader, and annotation labels placed before. A label met counts once,
    whichever of its parts or leader segments the letters meet, as in Linear
    (OV-3016). Letters keep ``clearance`` from a part, and ``clearance`` plus
    the segment's half width from a leader segment. ``made_way`` are the
    places along their regions that the feature labels made way for
    (CIR-RGN-YIELD, OV-3017): the label of each named annotation is drawn
    there, not moved.
    """

    labels: tuple[CircularObstacleLabel, ...] = ()
    clearance: float = 0.0
    made_way: tuple[RegionLabelPlace, ...] = ()

    def with_labels(self, labels: Sequence[CircularAnnotationLabel]) -> CircularLabelObstacles:
        return replace(self, labels=(*self.labels, *(CircularObstacleLabel(label.footprint) for label in labels)))


def _box(points: Sequence[Point]) -> Aabb:
    return Aabb(min(x for x, _ in points), min(y for _, y in points), max(x for x, _ in points), max(y for _, y in points))


def _grown_label_parts(parts: Sequence[Sequence[Point]], by: float) -> list[BoxedPart]:
    """A label's convex parts grown by ``by``, each with its box (``boxed_parts_meet``)."""

    return [(_box(part), part) for part in (grown_convex_polygon(part, by) for part in parts)]


class _IndexedLabel(NamedTuple):
    """An obstacle label with its boxes: the box around its parts (``None`` without parts), each part with its box,
    and each leader segment with its box and half width."""

    box: Aabb | None
    parts: list[tuple[Aabb, Polygon]]
    segments: list[tuple[Aabb, tuple[Point, Point], float]]
    feature: bool


class _Overlaps(NamedTuple):
    """What a label's letters overlap with no clearance kept (``_ObstacleIndex.overlaps``): the obstacle labels, by
    their index (the warning counts them), whether the letters meet a feature label's letters, not only its leader
    (OV-3017), whether they meet the letters (or backing) of a region label placed before them, and whether they meet
    a feature label's leader (OV-3156).
    """

    met: frozenset[int]
    feature_letters: bool
    region_letters: bool
    feature_leaders: bool

    @property
    def only_leaders(self) -> bool:
        """Whether only feature labels' leader lines run behind the letters: the label may be backed in front of
        them (CIR-RGN-FRONT)."""

        return self.feature_leaders and not (self.feature_letters or self.region_letters)


class _ObstacleIndex:
    """What a label's letters meet: an obstacle label's parts within ``clearance``, or one of its leader segments
    within its own; parts met by several parts of a curved label count once. ``conflicts`` is the placer's cost:
    the labels whose parts the letters meet plus the leader segments they meet, so a place over a label's letters
    and its leader costs more than one over a leader alone. ``overlaps`` counts the labels the letters overlap, with
    no clearance kept, each once, whether its parts, its leader or both are met (OV-3015, OV-3016; one per label,
    as in Linear), and says whether a feature label's letters are among them (OV-3017). ``add`` takes a label placed
    since, so a track's labels are boxed once each, not once per later label (OV-2664)."""

    def __init__(self, obstacles: CircularLabelObstacles) -> None:
        self._obstacles = obstacles
        self._labels: list[_IndexedLabel] = []
        for label in obstacles.labels:
            self._add(label)
        self._half_widths = {segment.half_width for label in obstacles.labels for segment in label.leader}

    def _add(self, label: CircularObstacleLabel) -> None:
        points = [point for part in label.parts for point in part]
        self._labels.append(
            _IndexedLabel(
                _box(points) if points else None,
                [(_box(part), part) for part in label.parts],
                [
                    (_box((segment.start, segment.end)), (segment.start, segment.end), segment.half_width)
                    for segment in label.leader
                ],
                label.feature,
            )
        )

    def add(self, label: CircularAnnotationLabel) -> None:
        """Count ``label`` as an obstacle from now on, as ``CircularLabelObstacles.with_labels`` would."""

        self._add(CircularObstacleLabel(label.footprint))

    def conflicts(self, label: CircularAnnotationLabel) -> int:
        return sum(int(parts) + segments for parts, segments in self._hits(label, self._obstacles.clearance))

    def overlaps(self, label: CircularAnnotationLabel) -> _Overlaps:
        hits = list(zip(self._labels, self._hits(label, 0.0)))
        return _Overlaps(
            frozenset(at for at, (_obstacle, (parts, segments)) in enumerate(hits) if parts or segments),
            any(obstacle.feature and parts for obstacle, (parts, _segments) in hits),
            any(not obstacle.feature and parts for obstacle, (parts, _segments) in hits),
            any(obstacle.feature and segments > 0 for obstacle, (_parts, segments) in hits),
        )

    def _hits(self, label: CircularAnnotationLabel, kept: float) -> list[tuple[bool, int]]:
        """For each obstacle label, whether the letters grown by ``kept`` meet its parts, and how many of its leader
        segments they meet within ``kept`` plus the segment's half width."""

        # An obstacle whose box misses the box around all the grown letters misses each part's box: it is skipped
        # before the parts are paired with it.
        near = _grown_label_parts(label.parts, kept)
        reach = union_aabbs(box for box, _ in near)
        if reach is None:
            return []
        # The letters grown by each distinct segment clearance (one per leader width), made once per call.
        line_near = {
            distance: _grown_label_parts(label.parts, distance)
            for distance in {kept + half_width for half_width in self._half_widths}
        }
        line_reach = {distance: union_aabbs(box for box, _ in parts) for distance, parts in line_near.items()}

        def meets_parts(obstacle: _IndexedLabel) -> bool:
            return obstacle.box is not None and reach.intersects(obstacle.box) and boxed_parts_meet(near, obstacle.parts)

        def segments_met(obstacle: _IndexedLabel) -> int:
            return sum(
                1
                for box, segment, half_width in obstacle.segments
                if (reach_at := line_reach[kept + half_width]) is not None
                and reach_at.intersects(box)
                and any(
                    part_box.intersects(box) and segment_intersects_polygon(segment, part, touching=False)
                    for part_box, part in line_near[kept + half_width]
                )
            )

        return [(meets_parts(obstacle), segments_met(obstacle)) for obstacle in self._labels]


def _origin_distance(polygon: Sequence[Point]) -> float:
    """How near the circle centre a convex polygon that does not hold it comes."""

    def to_segment(a: Point, b: Point) -> float:
        dx, dy = b[0] - a[0], b[1] - a[1]
        length = dx * dx + dy * dy
        t = 0.0 if length == 0.0 else max(0.0, min(1.0, -(a[0] * dx + a[1] * dy) / length))
        return math.hypot(a[0] + t * dx, a[1] + t * dy)

    return min(to_segment(polygon[at], polygon[(at + 1) % len(polygon)]) for at in range(len(polygon)))


def _backed(label: CircularAnnotationLabel, *, outward: bool) -> CircularAnnotationLabel:
    """``label`` on its backing (CIR-RGN-FRONT, ``label_backing``): its letters padded as Linear's, but toward its
    mark no further than ``label_offset`` (the gap the letters keep from the mark), so it never covers its mark."""

    near, far = label.ink_radii
    if label.arc is not None:
        anchor = math.atan2(label.y, label.x)
        turns = [
            anchor + (math.atan2(y, x) - anchor + math.pi) % (2.0 * math.pi) - math.pi for part in label.parts for x, y in part
        ]
        edge = near - label.label_offset if outward else far + label.label_offset
        return replace(label, backing_shape=sector_backing(near, far, min(turns), max(turns), label.font_size, edge=edge))
    # Straight text: its ink box in its own frame (turned back about its centre), and the direction of the centre
    # in that frame (toward the mark when outward).
    turn = math.radians(label.rotation)
    cos, sin = math.cos(turn), math.sin(turn)
    ink = _box([
        (label.x + (x - label.x) * cos + (y - label.y) * sin, label.y - (x - label.x) * sin + (y - label.y) * cos)
        for x, y in label.parts[0]
    ])
    toward = (-label.x * cos - label.y * sin, label.x * sin - label.y * cos)
    if not outward:
        toward = (-toward[0], -toward[1])

    def clear_of_mark(polygon: Sequence[Point]) -> bool:
        if outward:
            return _origin_distance(polygon) >= near - label.label_offset - 1e-6
        return max(math.hypot(x, y) for x, y in polygon) <= far + label.label_offset + 1e-6

    shape = turned_backing(
        ink, label.font_size, rotation=label.rotation, center=(label.x, label.y), toward=toward, clear=clear_of_mark
    )
    return replace(label, backing_shape=shape)


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
) -> tuple[float, CircularAnnotationLabel]:
    """The label moved along its own region to the clear place nearest its default (OV-1173, the Linear rule), and
    that place.

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

    place = nearest_fewest_conflicts(centre, candidates, count)
    return place, labels[place]


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
    a place where it keeps a letter. A label whose letters meet only leader
    lines there is drawn in front of them on a backing (``_backed``,
    CIR-RGN-FRONT) unless the backing would cover another label's letters or
    backing. A label whose letters still meet a feature label's letters there
    is reported
    (``CircularAnnotationLabels.over_feature_letters``); a label the feature
    labels then made way for is drawn at that place
    (``CircularLabelObstacles.made_way``, CIR-RGN-YIELD). A label whose cut
    at its midpoint keeps no letter, or that the track leaves out for want of
    a letter, is left out (``CircularAnnotationLabels.without_letters``).
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
    overlapped: dict[str, frozenset[int]] = {}
    over_letters: list[RegionLabelPlace] = []
    made_way = {place.annotation_id: place.place_bp for place in obstacles.made_way}
    for (annotation, label_at), default in kept:
        # A place the feature labels made way at is kept while its letters still meet neither theirs nor an earlier
        # region label's (a legend may have moved feature labels since); else the label moves along its region.
        pinned = label_at(made_way[annotation.id] % max(1, record_length)) if annotation.id in made_way else None
        # What the placed letters still overlap, without the clearance the places keep (OV-3015).
        met = index.overlaps(pinned) if pinned is not None else None
        if pinned is not None and met is not None and not (met.feature_letters or met.region_letters):
            place, label = made_way[annotation.id], pinned
        else:
            place, label = _along_region(
                annotation,
                default,
                label_at,
                record_length=record_length,
                fits=lambda label: not label.keeps_no_letter
                and _reach(label, inner_radius_px, outer_radius_px, outward) <= room + 1e-9,
                conflicts=index.conflicts,
                join_origin=track.join_origin,
            )
            met = index.overlaps(label)
        if met.only_leaders:
            # Only leader lines run behind the letters: the label is drawn in front of them on a backing that covers
            # no other label's letters or backing (CIR-RGN-FRONT), and a leader behind it is no conflict.
            backed = _backed(label, outward=outward)
            behind = index.overlaps(replace(label, parts=backed.backing))
            if not (behind.feature_letters or behind.region_letters):
                label, met = backed, met._replace(met=frozenset())
        overlapped[annotation.id] = met.met
        if met.feature_letters:
            over_letters.append(RegionLabelPlace(annotation.id, place))
        placed.append(label)
        index.add(label)
    return CircularAnnotationLabels(
        tuple(placed),
        without_letters,
        overlapped,
        {
            label.annotation_id: BoxedLetters(
                tuple(_grown_label_parts(label.parts, obstacles.clearance)), label.height + 2.0 * obstacles.clearance
            )
            for label in placed
        },
        tuple(over_letters),
    )


class CircularLaneLabels(NamedTuple):
    """A track's labels as the lane packer's checks read them (``circular_annotation_lane_labels``): the ring its
    lanes take, the record length they were placed on, and the labels drawn for each record (record index, labels),
    placed among the track's own labels only."""

    ring: AnnotationRing
    record_length: int
    labels: tuple[tuple[int, tuple[CircularAnnotationLabel, ...]], ...]


# ``circular_annotation_lane_labels`` bound to where a slot's lanes stack: the labels of a track, its parameters and
# the extent its lanes take.
LaneLabelPlacer: TypeAlias = Callable[[ResolvedAnnotationTrack, AnnotationTrackParams, float], CircularLaneLabels]


def circular_annotation_lane_labels(
    track: ResolvedAnnotationTrack,
    params: AnnotationTrackParams,
    extent_px: float,
    *,
    font_family: str,
    record_length: int,
    lanes_from_px: float,
    outward: bool,
    clearance_px: float,
) -> CircularLaneLabels:
    """The labels ``circular_annotation_labels`` draws for each record of ``track``, the lanes taking ``extent_px``
    from ``lanes_from_px``, the radius they stack from (away from the centre, or toward it when not ``outward``), with
    the track's own labels as the only obstacles, kept ``clearance_px`` apart as where the drawing places them
    (``CircularLabelObstacles.clearance``): each label moved along its region where it comes that near one placed
    before it. The one placement the lane packer's checks read (``circular_annotation_lanes_meeting``,
    ``circular_annotation_labels_across_centre``)."""

    inner, outer = (lanes_from_px, lanes_from_px + extent_px) if outward else (lanes_from_px - extent_px, lanes_from_px)
    ring = AnnotationRing(max(0.0, inner), outer, "outside" if outward else "inside")
    return CircularLaneLabels(
        ring,
        record_length,
        tuple(
            (
                record_index,
                circular_annotation_labels(
                    track, record_index=record_index, record_length=record_length, inner_radius_px=ring.inner_radius_px,
                    outer_radius_px=ring.outer_radius_px, side=ring.side, font_family=font_family, params=params,
                    obstacles=CircularLabelObstacles(clearance=clearance_px),
                ).drawn,
            )
            for record_index in sorted({item.annotation.record_index for item in track.placements})
        ),
    )


def circular_annotation_lanes_meeting(
    track: ResolvedAnnotationTrack,
    params: AnnotationTrackParams,
    extent_px: float,
    *,
    placed: LaneLabelPlacer,
) -> frozenset[tuple[int, int, int]]:
    """The lanes (record index, row, lane) in which ``circular_annotation_labels`` draws two of ``track``'s labels over
    each other, the lanes taking ``extent_px``, its labels ``placed`` (``circular_annotation_lane_labels``). The lane
    packer's measure (``LabelSizing.meet``)."""

    met: set[tuple[int, int, int]] = set()
    for record_index, drawn in placed(track, params, extent_px).labels:
        lanes = {item.annotation.id: item.lane for item in track.placements if item.annotation.record_index == record_index}
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
    placed: LaneLabelPlacer,
) -> frozenset[tuple[int, str]]:
    """The labels (record index, annotation id) of an inside track whose letters ``circular_annotation_labels`` draws
    over the centre of the circle and within their ``label_offset`` of another annotation's painted mark, or across the
    paint of their own mark (``_mark_ink``), its lanes taking ``extent_px`` toward the centre in the ring, and its labels,
    ``placed`` as ``circular_annotation_lanes_meeting`` reads them (``circular_annotation_lane_labels``).
    Letters too wide for the circle their mark leaves them are centred on the centre of the circle and come back across
    it, where no room keeps them clear of the marks they meet (LBL-ACROSS-CENTRE-OMIT; OV-2991, OV-3035). Letters that
    only come nearer their own mark than their ``label_offset`` are not counted. The lane packer's measure
    (``LabelSizing.across_centre``)."""

    ring, record_length, labels = placed(track, params, extent_px)
    met: set[tuple[int, str]] = set()
    for record_index, drawn in labels:
        across = [label for label in drawn if label.ink_radii[0] <= 0.0]
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
                convex_polygons_intersect(part if identity == label.annotation_id else grown_convex_polygon(part, label.label_offset), mark)
                for part in label.parts
                for identity, parts in marks.items()
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

_Measured = TypeVar("_Measured")


def _once_per_text(measure: Callable[[tuple[TextRun, ...]], _Measured]) -> Callable[[tuple[TextRun, ...]], _Measured]:
    """``measure``, run once for each text one curved label tries (``_curved_label``)."""

    kept: dict[tuple[TextRun, ...], _Measured] = {}

    def measured(shown: tuple[TextRun, ...]) -> _Measured:
        if shown not in kept:
            kept[shown] = measure(shown)
        return kept[shown]

    return measured


def _cleared_radius(
    extent: Callable[[float], tuple[float, float]],
    radius_for: Callable[[float], float],
    *,
    mark_edge: float,
    outward: bool,
    gap: float,
    rings: Sequence[RadialBand],
) -> tuple[float, tuple[float, float]]:
    """The anchor radius whose letters keep ``gap`` beyond ``mark_edge``, and beyond the paint of any ring of
    ``rings`` they would otherwise meet or come within ``gap`` of; and how near the centre and how far from it the
    letters come there (``extent``, measured once for each radius tried)."""

    edge = mark_edge
    radius = radius_for(edge + gap if outward else edge - gap)
    reach: tuple[float, float] | None = None
    for ring in sorted(rings, key=lambda band: band.inner_px if outward else -band.outer_px):
        near, far = reach = extent(radius) if reach is None else reach
        if outward and edge < ring.outer_px and ring.inner_px < far + gap:
            edge = ring.outer_px
        elif not outward and ring.inner_px < edge and near - gap < ring.outer_px:
            edge = ring.inner_px
        else:
            continue
        radius = radius_for(edge + gap if outward else edge - gap)
        reach = None
    return radius, extent(radius) if reach is None else reach


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
    radius, ink_radii = clear(
        partial(_ink_extent, direction=direction, offsets=offsets),
        partial(_label_radius, outward=outward, direction=direction, offsets=offsets),
    )
    x, y = circle_point(radius, angle)
    left, right, top, bottom = offsets
    cos, sin = math.cos(turn), math.sin(turn)
    return label(
        runs=runs,
        x=x,
        y=y,
        rotation=rotation,
        ink_radii=ink_radii,
        height=bottom - top,
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
    clear: Callable[[Callable[[float], tuple[float, float]], Callable[[float], float]], tuple[float, tuple[float, float]]],
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
    glyphs = _once_per_text(glyphs)

    @_once_per_text
    def placed(shown: tuple[TextRun, ...]) -> tuple[float, float, float, float, float]:
        """The ink of ``shown`` along its baseline (left, right), how near the centre and how far from it its
        letters come, and its baseline radius, cleared of its mark (and rings): measured once for each text tried."""

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

        radius, (near, far) = clear(extent, radius_for)
        return left, right, near, far, radius

    @_once_per_text
    def standing(shown: tuple[TextRun, ...]) -> GlyphsAlongCircle:
        """The glyphs of ``shown`` standing along their baseline circle (``placed``)."""

        return glyphs_along_circle(glyphs(shown), placed(shown)[4], tops_in=flipped)

    def fits(shown: tuple[TextRun, ...]) -> bool:
        # The ink reaches at most this far each way from the label's middle along the circle its letters come
        # nearest the centre on, where they stand closest together (OV-2806): a circle about as small as the
        # letters, or letters pressed against the centre past their mark, leave no room for one. And no two letters'
        # ink overlaps: letters that stand apart on a line close up at their tops on the inner side of a circle.
        left, right, near, _far, _radius = placed(shown)
        return max(-left, right) <= math.radians(_MAX_HALF_PATH_DEG) * near - font_size and not standing(shown).meet()

    shown = runs if fits(runs) else _cut_to(runs, fits)
    left, right, near, far, radius = placed(shown)
    # The letters' shape turns as far each way as their ink box along the baseline; nearer the centre than the
    # baseline it turns as far as their ink does where that is further: a glyph does not bend to the circle, so there
    # it spans a wider turn than its arc on the baseline (OV-2851). Where no glyph is known, the ink box is measured
    # at the circle the letters come nearest the centre on.
    half_turn = 0.5 * (right - left) / radius
    ink_turn = standing(shown).half_turn()
    inner_turn = 0.5 * (right - left) / near if ink_turn is None else ink_turn
    # Along the reading direction an offset is an arc length on the baseline; anticlockwise when flipped.
    sign = -1.0 if flipped else 1.0
    centre = angle + sign * math.degrees(0.5 * (left + right) / radius)
    x, y = circle_point(radius, angle)
    return label(
        runs=shown,
        x=x,
        y=y,
        rotation=0.0,
        ink_radii=(near, far),
        height=far - near,
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
    dom_group_id: str,
    semantic_slot_id: str,
) -> GroupContent:
    """Build an annotation track around the circle centre (render stage ``annotation_content``).

    ``labels`` are the track's labels (``circular_annotation_labels``), drawn
    with their annotations. A region with a strand is drawn in its strand's
    lane of ``ring.strand_bands_px`` (a slot that covers a feature track that
    separates the strands), any other in the slot's band. The hatch patterns the marks use, in order of first
    use, and the clip path of a clipped track are the content's definitions.
    """

    group = Group(id=dom_group_id)
    definitions: dict[str, BaseElement] = {}
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
        if label is not None and not label.backed:  # else in front of the leaders (circular_region_label_front_content)
            item_group.add(_label_text(label, font_family, definitions, str(item_group.attribs["id"])))
        group.add(item_group)
    clip_paths: list[BaseElement] = []
    if track.clipped:
        group.attribs["data-gbdraw-annotation-clipped"] = "true"
        clip_id = stable_svg_id("annotation_clip", "circular", dom_group_id)
        clip_ring = Path(
            d=f"{circle_arc(ring.outer_radius_px, 0.0, -360.0)} Z {circle_arc(ring.inner_radius_px, 0.0, 360.0)} Z"
        )
        clip_ring.attribs["clip-rule"] = "evenodd"
        clip_path = ClipPath(id=clip_id)
        clip_path.add(clip_ring)
        clip_paths.append(clip_path)
        group.attribs["clip-path"] = f"url(#{clip_id})"
    return GroupContent.of(group, (*definitions.values(), *clip_paths))


def circular_region_label_front_content(
    track: ResolvedAnnotationTrack,
    *,
    record_index: int,
    font_family: str,
    labels: Sequence[CircularAnnotationLabel],
    dom_group_id: str,
) -> GroupContent | None:
    """The track's backed labels in front of the record's leader lines (CIR-RGN-FRONT, OV-3156), around the circle
    centre; ``None`` when it has none. Each label is drawn on its ``backing_shape``, white at 70 % opacity
    (``region_label_front``, as in Linear); a curved label's path goes into the content's definitions."""

    definitions: dict[str, BaseElement] = {}
    backed: list[BackedRegionLabel] = []
    for label in labels:
        if label.backing_shape is None:
            continue
        item_id = annotation_dom_id(
            record_index=record_index, slot_id=track.slot_id, set_id=track.set_id, annotation_id=label.annotation_id
        )
        backed.append(
            BackedRegionLabel(item_id, label.backing_shape.element(), _label_text(label, font_family, definitions, item_id))
        )
    group = region_label_front(dom_group_id, backed)
    return None if group is None else GroupContent.of(group, tuple(definitions.values()))


__all__ = [
    "CURVED_LABEL_ORIENTATIONS",
    "AnnotationRing",
    "CircularAnnotationLabel",
    "CircularAnnotationLabels",
    "CircularLabelArc",
    "CircularLabelObstacles",
    "RegionLabelPlace",
    "CircularObstacleLabel",
    "CircularLabelRoom",
    "CircularLaneLabels",
    "LaneLabelPlacer",
    "LeaderSegment",
    "circular_annotation_label_extent",
    "circular_annotation_lane_labels",
    "circular_annotation_mark_past_lane",
    "circular_annotation_mark_reach",
    "circular_annotation_label_room",
    "circular_annotation_labels",
    "circular_annotation_labels_across_centre",
    "circular_annotation_track_content",
    "circular_annotation_track_ink",
    "circular_region_label_front_content",
]
