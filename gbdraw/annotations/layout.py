"""Deterministic lane packing shared by circular and linear annotations."""

from __future__ import annotations

import math
from dataclasses import dataclass, field, replace
from typing import Callable, Iterable, Literal, Mapping, NamedTuple, Sequence

from gbdraw.exceptions import ValidationError

from .models import (
    AnnotationTrackParams,
    LabelPlace,
    LabelPosition,
    RegionAnnotationStyle,
    ResolvedRegionAnnotation,
    annotation_label_font_size,
    effective_annotation_style,
)


@dataclass(frozen=True)
class PlacedRegionAnnotation:
    annotation: ResolvedRegionAnnotation
    lane: int


def lane_count_of(placements: Iterable[PlacedRegionAnnotation]) -> int:
    """How many lanes ``placements`` take: one past the deepest. The drawers share a record's band among its lanes by
    it, and the lane packer measures a record's lanes the same way."""

    return max((item.lane for item in placements), default=-1) + 1


class LaneRoom(NamedTuple):
    """Room between lane ``gap`` and lane ``gap + 1`` on one row of one record for the labels that would cross
    the mark of lane ``gap + 1`` (labels toward it) or of lane ``gap`` (labels toward lane 0) (PD-OI-119, lane
    room): their letters and their ``label_offset`` on each side of the letters, beside the painted edge of the
    label's own mark and before the painted edge of the mark it would cross (how far each mark paints past its lane,
    ``LabelSizing.reach``), less the lanes and rooms between a label's mark and that lane."""

    record_index: int
    row: int
    gap: int
    px: float


# The extent of one lane without label room: the lane packer measures lanes by it, and the drawers draw each lane of
# a clipped slot that wide.
LANE_EXTENT_PX = 14.0

# What stops a slot's lanes and rooms (``_LaneBand.stop``): the centre of the circle, where an inside Circular slot
# without a set width runs out of room (lanes past it share the band up to it, as the drawer draws them), or the
# set size of a slot whose overflow is clip (lanes past it are drawn ``lane_extent_px`` apart and the clip cuts them).
RoomStop = Literal["centre", "clip"]


class RoomBounds(NamedTuple):
    """What bounds a track's lane rooms (``layout_annotation_track``): the room each gap keeps at least
    (``at_least``, an earlier packing's), and the most the track may take with its rooms (``limit_px``; None for no
    limit): for an inside Circular slot, the radius of the edge its lanes stack from (the centre of the circle). A
    slot whose size is set and whose overflow is clip is bounded by that size as well, whichever is nearer
    (``_lane_band``)."""

    at_least: tuple[LaneRoom, ...] = ()
    limit_px: float | None = None


def lane_width_px(band_px: float, *, padding_px: float, room_px: float, lane_count: int) -> float:
    """The width of each lane where ``lane_count`` lanes share a band ``band_px`` wide: the band less its two
    paddings and the slot's label room (``ResolvedAnnotationTrack.lane_room_px``), the lane gaps shared out among
    the lanes. The drawers place a slot's lanes by it, and the lane packer measures labels where it puts them."""

    return max(1.0, band_px - 2.0 * padding_px - room_px) / max(1, lane_count)


def lane_edge_px(
    lane: float, width_px: float, *, padding_px: float, room_before_px: float, from_px: float = 0.0, toward: Literal[1, -1] = 1
) -> float:
    """Where a slot's lane ``lane`` starts, from the edge ``from_px`` its lanes stack from (``toward`` 1 when they
    stack toward larger values, -1 toward smaller ones): past the band's padding, the lanes before it (``width_px``
    each, ``lane_width_px``) and the label rooms of the gaps before it (``ResolvedAnnotationTrack.room_before``).
    ``lane + 0.5`` gives the lane's middle. The drawers place lanes by it and the lane packer measures labels there;
    the terms are added in this order so that every caller gets the same pixels."""

    return from_px + toward * padding_px + toward * lane * width_px + toward * room_before_px


class _LaneBand(NamedTuple):
    """The band the lane packer puts a slot's lanes in: ``lane_count`` lanes of ``lane_extent_px`` with
    ``lane_gap_px`` between them and ``padding_px`` on each side, stopped at ``limit_px`` by ``stop`` (``RoomStop``;
    no limit when None): lanes that pass the centre share the band up to it (``lane_width_px``), and lanes that pass
    a clipped set size are drawn ``lane_extent_px`` apart. ``labels_above`` as in ``LabelSizing``. ``set_px`` is
    the extent the drawer gives a slot whose size is set (that size, stopped at the centre for an inside Circular
    slot; None where the slot takes the extent its lanes and rooms need): its lanes share it less the slot's label
    room (``width_px``)."""

    lane_extent_px: float
    lane_gap_px: float
    padding_px: float
    limit_px: float | None
    stop: RoomStop
    labels_above: bool
    set_px: float | None

    def bare_px(self, lane_count: int) -> float:
        """The extent of ``lane_count`` lanes without label room, past any limit."""

        return 2.0 * self.padding_px + lane_count * self.lane_extent_px + max(0, lane_count - 1) * self.lane_gap_px

    def passes(self, lane_count: int) -> bool:
        """Whether ``lane_count`` lanes pass the limit: they share the band up to the centre, or a clipped slot
        draws them ``lane_extent_px`` apart (``width_px``)."""

        return self.limit_px is not None and self.bare_px(lane_count) > self.limit_px

    def drawn_px(self, lane_count: int) -> float:
        """The extent ``lane_count`` lanes take without label room where the drawer draws them: lanes past a clipped
        set size touch, ``lane_extent_px`` apart, between the two paddings."""

        if self.stop == "clip" and self.passes(lane_count):
            return 2.0 * self.padding_px + lane_count * self.lane_extent_px
        return self.bare_px(lane_count)

    def width_px(self, record_lanes: int, slot_lanes: int, *, room_px: float | None = None) -> float:
        """The width of each lane of a record with ``record_lanes`` lanes in a slot whose records take at most
        ``slot_lanes``, where the drawer puts them (``lane_width_px``): the record's lanes share the band the slot's
        lanes take (each record's lanes fill the slot, as the drawers draw them). Lanes that pass the limit share
        the band up to it, as the drawer draws an inside slot clamped at the centre; a clipped slot draws each of
        them ``lane_extent_px`` wide. With the slot's label room (``room_px``, ``ResolvedAnnotationTrack
        .lane_room_px``) the lanes of a slot of set size share that size less the room (``set_px``), as the drawer
        draws them; without it (while the rooms are reckoned) they are measured at the extent they take."""

        lanes = max(1, slot_lanes)
        band = self.bare_px(lanes)
        if self.limit_px is not None and band > self.limit_px:
            if self.stop == "clip":
                return self.lane_extent_px
            band = self.limit_px
        elif self.set_px is not None and room_px is not None:
            return lane_width_px(self.set_px, padding_px=self.padding_px, room_px=room_px, lane_count=record_lanes)
        return lane_width_px(band, padding_px=self.padding_px, room_px=0.0, lane_count=record_lanes)

    def drawn_extent_px(self, lane_count: int, room_px: float) -> float:
        """The extent the drawer gives ``lane_count`` lanes and the slot's label room ``room_px``: the set size
        (``set_px``), or the lanes and the room, stopped at any limit."""

        if self.set_px is not None:
            return self.set_px
        extent = self.bare_px(lane_count) + room_px
        return extent if self.limit_px is None else min(extent, self.limit_px)

    def mark_edge_px(self, lane: int, width: float, *, room_before_px: float, past_lane_edge_px: float) -> float:
        """How far the painted edge of a mark in ``lane`` of lanes ``width`` wide (``width_px``) lies from where the
        lanes stack, on its label's side, past the rooms before its lane: ``past_lane_edge_px`` beyond the lane's
        edge on that side (``LabelSizing.reach``; negative while the paint stays inside the lane)."""

        if self.labels_above:
            return lane_edge_px(lane, width, padding_px=self.padding_px, room_before_px=room_before_px) - past_lane_edge_px
        return lane_edge_px(lane + 1, width, padding_px=self.padding_px, room_before_px=room_before_px) + past_lane_edge_px

    def hides_mark(self, lane: int, width: float, *, room_before_px: float, past_lane_edge_px: float) -> bool:
        """Whether a clipped set size (``limit_px``) hides all of a mark in ``lane`` of lanes ``width`` wide: its
        paint (``past_lane_edge_px`` beyond the lane on each side) starts at or past the set size."""

        start = lane_edge_px(lane, width, padding_px=self.padding_px, room_before_px=room_before_px) - past_lane_edge_px
        return self.stop == "clip" and self.limit_px is not None and start >= self.limit_px


class LabelExtent(NamedTuple):
    """A label's painted letters at its place, as the drawer measures and draws them: the bp they reach along the
    region from the place (``low_bp`` before it, ``high_bp`` after it, both offsets from the place) and the pixels
    they take across the lanes. ``drawn`` is the drawer's verdict: False where it leaves the label out there because
    no letter of it fits its circle (RGN-CUT-EMPTY-OMIT); such a label takes no room, is no obstacle and is not cut
    by a clip."""

    low_bp: float
    high_bp: float
    depth_px: float
    drawn: bool


# A label's letters (``LabelExtent``) from its annotation, its effective style, its font size, the bp it stands at
# and how far its mark's edge lies from the edge the slot's lanes stack from.
LabelMeasure = Callable[[ResolvedRegionAnnotation, RegionAnnotationStyle, float, float, float], LabelExtent]

# How far a mark's paint reaches past the edge of its lane on its label's side (negative while it stays inside the
# lane), from its annotation, its effective style and the width of its lane: the drawer's own measure, where it
# places the label (``label_offset`` beyond that paint).
MarkReachMeasure = Callable[[ResolvedRegionAnnotation, RegionAnnotationStyle, float], float]

# A labelled annotation's letters standing at a bp beside its mark in a lane, as the drawer measures them there
# (``LabelMeasure``); None where the track draws no label for it.
LettersInLane = Callable[[ResolvedRegionAnnotation, float, int], "LabelExtent | None"]

# The lanes (record index, row, lane) in which the drawer draws two labels of a track over each other, from the track,
# its parameters and the extent its lanes and rooms take where the drawer draws them (``_LaneBand.drawn_extent_px``).
LanesMeeting = Callable[["ResolvedAnnotationTrack", AnnotationTrackParams, float], frozenset[tuple[int, int, int]]]

# The labels (record index, annotation id) of an inside Circular track whose letters the drawer draws over the centre of
# the circle onto a mark, from the track, its parameters and the extent its lanes and rooms take where the drawer
# draws them (``_LaneBand.drawn_extent_px``).
LabelsAcrossCentre = Callable[["ResolvedAnnotationTrack", AnnotationTrackParams, float], frozenset[tuple[int, str]]]


@dataclass(frozen=True)
class LabelSizing:
    """How the packer measures a track's labels.

    ``font_size`` is the font size of the drawing's feature labels, which a label whose style sets none takes;
    ``bp_per_px`` the bases one pixel spans on each record (1 when not given). ``extent`` is the drawer's own
    measure of a label's letters at its place (``LabelMeasure``), and ``reach`` its measure of how far a mark paints
    past its lane (``MarkReachMeasure``), from which the label stands. ``labels_above`` says that the labels stand
    toward lane - 1 (Linear labels above their marks); else they stand toward lane + 1 (Linear labels below, and
    Circular, away from the axis on either side). ``meet`` is the drawer's own placement of a track's labels
    (``LanesMeeting``) where a label's letters depend on where its lane is drawn (Circular); None where they do not.
    ``across_centre`` is the drawer's own measure of the labels whose letters it draws over the centre of the circle
    onto a mark (``LabelsAcrossCentre``, an inside Circular slot where it is drawn); None where letters cannot pass it.
    """

    font_size: float
    extent: LabelMeasure
    reach: MarkReachMeasure
    bp_per_px: Mapping[int, float] = field(default_factory=dict)
    labels_above: bool = False
    meet: LanesMeeting | None = None
    across_centre: LabelsAcrossCentre | None = None


@dataclass(frozen=True)
class ResolvedAnnotationTrack:
    slot_id: str
    set_id: str
    placements: tuple[PlacedRegionAnnotation, ...]
    required_extent_px: float
    # The feature label font size of the drawing: the size of a label whose style sets none.
    label_font_size: float
    # Whether the drawing closes into a circle (Circular), so the two sides of a span across the origin always
    # meet (``mark_parts``): the packer placed the labels by it, and the drawers draw the marks and labels by it.
    join_origin: bool
    clipped: bool = False
    lane_gap_px: float = 3.0
    # The room between lanes for labels that would cross another lane's mark (``LaneRoom``).
    lane_rooms: tuple[LaneRoom, ...] = ()
    # Whether the labels stand toward lane 0 (``LabelSizing.labels_above``).
    labels_above: bool = False
    # The labels (record index, annotation id) the track leaves out, and what stopped their room (``RoomStop``; None
    # where nothing bounds the rooms): labels that would cross another lane's mark and were left without their room
    # (an inside Circular slot adds no room past the centre, a clipped slot of set size none past that size), and
    # labels whose letters the edge of a clipped slot would cut (``_cut_by_clip``). Their marks are drawn.
    without_room: tuple[tuple[int, str], ...] = ()
    rooms_stop: RoomStop | None = None
    # The font size of each label the track leaves out (record index, size; ``annotation_label_font_size``), in
    # the order of ``without_room``: what the label fit warning reports for them (``left_out_font_sizes``).
    without_room_font_sizes: tuple[tuple[int, float], ...] = ()
    # The labels (record index, annotation id) the track leaves out because no letter of them fits their circle
    # where their text would put their mark (RGN-CUT-EMPTY-OMIT, OV-2805; ``layout_annotation_track``): their marks
    # stand where they stand without that text. The drawer counts them with the labels whose cut keeps no letter.
    without_letters: tuple[tuple[int, str], ...] = ()

    def labels_without_room(self, record_index: int) -> int:
        """How many of the record's labels the track leaves out (``without_room``)."""

        return sum(1 for record, _identity in self.without_room if record == record_index)

    def left_out_font_sizes(self, record_index: int) -> tuple[float, ...]:
        """The font sizes of the record's labels the track leaves out (``without_room_font_sizes``)."""

        return tuple(size for record, size in self.without_room_font_sizes if record == record_index)

    def draws_label(self, annotation: ResolvedRegionAnnotation) -> bool:
        """Whether the track draws the annotation's label: it has one, and the track does not leave it out
        (``without_room``, ``without_letters``). The drawers, the label placement and the label room a slot reserves
        read it."""

        key = (annotation.record_index, annotation.id)
        return bool(annotation.label) and key not in self.without_room and key not in self.without_letters

    def leaves_out_letters(self, annotation: ResolvedRegionAnnotation) -> bool:
        """Whether the track leaves the annotation's label out because no letter of it fits its circle where its
        text would put its mark (``without_letters``)."""

        return (annotation.record_index, annotation.id) in self.without_letters

    def room_before(self, record_index: int, row: int, lane: int) -> float:
        """How far ``lane`` of the record's ``row`` moves on from where it stands without label room: the rooms of
        the gaps before it there and, where the labels stand toward lane 0, the slot's room the row does not take
        (``lane_room_px``).

        Labels stand away from the axis, so a row that takes less room than the slot keeps its lanes where they
        stand beside the axis, and leaves the rest empty on the labels' side."""

        rooms = [room for room in self.lane_rooms if room.record_index == record_index and room.row == row]
        before = sum(room.px for room in rooms if room.gap < lane)
        if not self.labels_above:
            return before
        return before + (self.lane_room_px - sum(room.px for room in rooms))

    @property
    def lane_room_px(self) -> float:
        """The room the slot keeps for its labels between lanes: the largest a record's row takes."""

        return _largest_room(self.lane_rooms)


def _largest_room(lane_rooms: Sequence[LaneRoom]) -> float:
    totals: dict[tuple[int, int], float] = {}
    for room in lane_rooms:
        totals[(room.record_index, room.row)] = totals.get((room.record_index, room.row), 0.0) + room.px
    return max(totals.values(), default=0.0)


def _letters_at(
    place_bp: float, low_bp: float, high_bp: float, *, record_length: int, join_origin: bool
) -> tuple[tuple[float, float], ...]:
    """The bp a label's letters cover on its row, ``low_bp`` before its place to ``high_bp`` after it.

    On a circle (``join_origin``) the letters go on across the origin; a Linear record cuts them at its ends (on a
    wrapped row, letters past the row's end lie where the row draws nothing).
    """

    start, end = place_bp + low_bp, place_bp + high_bp
    if not join_origin:
        return ((max(0.0, start), min(float(record_length), end)),)
    # A part joined across the origin runs on past the record's end.
    turns = math.floor(start / record_length)
    start, end = start - turns * record_length, end - turns * record_length
    if end > record_length:
        return ((start, float(record_length)), (0.0, end - record_length))
    return ((start, end),)


class _LabelLetters(NamedTuple):
    """Where a label stands, the bp its letters span on its row, and whether they leave the stretch of their mark
    they may cover (``ResolvedRegionAnnotation.label_stretch``)."""

    place: LabelPlace
    intervals: tuple[tuple[float, float], ...]
    leaves_mark: bool


def _label_letters(
    place: LabelPlace,
    stretch: tuple[int, int],
    letters: LabelExtent,
    *,
    record_length: int,
    join_origin: bool,
) -> _LabelLetters:
    """What the lane reservation keeps for a label standing at ``place`` (``label_stretch``): the bp its
    ``letters`` paint there (``join_origin`` as in ``mark_parts``, ``_letters_at``), and whether they leave the
    ``stretch`` of their mark."""

    low, high = stretch
    start, end = place.bp + letters.low_bp, place.bp + letters.high_bp
    intervals = _letters_at(place.bp, letters.low_bp, letters.high_bp, record_length=record_length, join_origin=join_origin)
    return _LabelLetters(place, intervals, start < low or end > high)


def _expanded_intervals(
    annotation: ResolvedRegionAnnotation,
    row: int,
    *,
    record_length: int,
    letters: _LabelLetters | None,
    padding_bp: float,
) -> tuple[tuple[float, float], ...]:
    """What an annotation occupies on one row, with its label's ``letters`` where they leave its mark
    (``_label_letters``)."""

    intervals = [(float(start), float(end)) for start, end in annotation.row_segments(row)]
    if letters is not None and letters.leaves_mark and row == letters.place.part.row:
        intervals.extend(letters.intervals)
    expanded: list[tuple[float, float]] = []
    for start, end in intervals:
        start -= padding_bp
        end += padding_bp
        if start < 0:
            expanded.append((0.0, min(float(record_length), end)))
            expanded.append((max(0.0, record_length + start), float(record_length)))
        elif end > record_length:
            expanded.append((max(0.0, start), float(record_length)))
            expanded.append((0.0, min(float(record_length), end - record_length)))
        else:
            expanded.append((max(0.0, start), min(float(record_length), end)))
    return tuple(sorted((start, end) for start, end in expanded if end > start))


def _overlaps(
    left: Sequence[tuple[float, float]],
    right: Sequence[tuple[float, float]],
) -> bool:
    return any(a_start < b_end and a_end > b_start for a_start, a_end in left for b_start, b_end in right)


def assign_annotation_lanes(
    annotations: Sequence[ResolvedRegionAnnotation],
    *,
    record_lengths: Mapping[int, int],
    letters_in_lane: LettersInLane | None = None,
    padding_bp: Callable[[ResolvedRegionAnnotation], float] | float = 0.0,
    overflow: str = "error",
    join_origin: bool,
    label_position: Callable[[ResolvedRegionAnnotation], LabelPosition] | None = None,
) -> tuple[PlacedRegionAnnotation, ...]:
    """Assign stable lanes using greedy interval coloring, row by row on a wrapped record.

    ``join_origin`` says whether the drawing closes into a circle (``mark_parts``) and ``label_position`` where
    each label stands on it (centred when not given): a lane keeps room for the label's letters there, as the
    drawer measures them beside the annotation's mark in that lane (``letters_in_lane``; no letters when not
    given).

    An annotation keeps one lane on every row it draws on; only annotations on the same row meet.
    Explicit lanes remain fixed. A collision in an explicit lane is rejected by
    the default ``error`` policy and retained for ``compress``/``clip`` so the
    caller can apply its visual policy.
    """

    if overflow not in {"error", "compress", "clip"}:
        raise ValidationError("overflow must be 'error', 'compress', or 'clip'.")

    def occupied_in(annotation: ResolvedRegionAnnotation) -> Callable[[int], dict[int, tuple[tuple[float, float], ...]]]:
        # What the annotation occupies on each of its rows with its mark in a lane: its letters are measured there.
        length = int(record_lengths[annotation.record_index])
        annotation_padding = max(0.0, float(padding_bp(annotation) if callable(padding_bp) else padding_bp))
        stretch = (
            annotation.label_stretch("center" if label_position is None else label_position(annotation), join_origin=join_origin)
            if annotation.label and letters_in_lane is not None
            else None
        )

        def rows_in(lane: int) -> dict[int, tuple[tuple[float, float], ...]]:
            extent = None if stretch is None or letters_in_lane is None else letters_in_lane(annotation, stretch[0].bp, lane)
            letters = (
                None
                if stretch is None or extent is None
                else _label_letters(*stretch, extent, record_length=length, join_origin=join_origin)
            )
            return {
                row: _expanded_intervals(annotation, row, record_length=length, letters=letters, padding_bp=annotation_padding)
                for row in annotation.geometry_rows
            }

        return rows_in

    def first_place(annotation: ResolvedRegionAnnotation) -> tuple[int, int]:
        row = annotation.geometry_rows[0]
        return row, min(start for start, _ in annotation.row_segments(row))

    ordered = sorted(
        annotations,
        key=lambda item: (
            item.record_index,
            first_place(item),
            max(end for _, end in item.row_segments(item.geometry_rows[-1])),
            item.set_id,
            item.id,
        ),
    )
    lane_bounds: dict[tuple[int, int, int], list[tuple[tuple[float, float], ...]]] = {}
    placements: list[PlacedRegionAnnotation] = []
    for annotation in ordered:
        rows_in = occupied_in(annotation)

        def lane_is_taken(lane: int, rows: Mapping[int, tuple[tuple[float, float], ...]]) -> bool:
            return any(
                _overlaps(bounds, occupied)
                for row, bounds in rows.items()
                for occupied in lane_bounds.get((annotation.record_index, row, lane), ())
            )

        if annotation.lane is not None:
            lane = int(annotation.lane)
            rows = rows_in(lane)
            if overflow == "error" and lane_is_taken(lane, rows):
                raise ValidationError(
                    f"Annotation {annotation.set_id}/{annotation.id} conflicts in explicit lane {lane}.",
                    diagnostic={"code": "TRACK_INVALID", "reason": "ANNOTATION_LANE_CONFLICT"},
                )
        else:
            lane = 0
            rows = rows_in(lane)
            while lane_is_taken(lane, rows):
                lane += 1
                rows = rows_in(lane)
        for row, bounds in rows.items():
            lane_bounds.setdefault((annotation.record_index, row, lane), []).append(bounds)
        placements.append(PlacedRegionAnnotation(annotation=annotation, lane=lane))
    placements.sort(
        key=lambda item: (
            item.annotation.record_index,
            item.lane,
            first_place(item.annotation),
            item.annotation.set_id,
            item.annotation.id,
        )
    )
    return tuple(placements)


class _LaneRooms(NamedTuple):
    """The rooms between lanes (``LaneRoom``), and the labels that would cross a mark but were left without the room
    they need (``ResolvedAnnotationTrack.without_room``)."""

    rooms: tuple[LaneRoom, ...]
    without_room: tuple[tuple[int, str], ...] = ()


class _TrackMeasures(NamedTuple):
    """The measures of one slot's annotations as its drawer makes them: each one's effective style (set, slot and
    annotation precedence), its label's place along the region, its label's font size and letters
    (``LabelSizing.extent``), and how far its mark paints past its lane (``LabelSizing.reach``). The labels of
    ``without_letters`` (record index, annotation id) are left out because no letter of them fits their circle
    (``ResolvedAnnotationTrack.without_letters``): the packer places their marks as if they had no text."""

    params: AnnotationTrackParams
    labels: LabelSizing
    join_origin: bool
    without_letters: frozenset[tuple[int, str]]
    # Each label's letters measured so far, by record index, annotation id, place and mark edge: the packing rounds
    # of one slot measure the same labels at the same places again (``_SlotPacking``).
    measured: dict[tuple[int, str, float, float], LabelExtent]

    def has_letters(self, annotation: ResolvedRegionAnnotation) -> bool:
        """Whether the annotation has a label the packer gives letters: one with text, not left out for want of a
        letter (``without_letters``)."""

        return bool(annotation.label) and (annotation.record_index, annotation.id) not in self.without_letters

    def style(self, annotation: ResolvedRegionAnnotation) -> RegionAnnotationStyle:
        return effective_annotation_style(annotation, self.params)

    def label_position(self, annotation: ResolvedRegionAnnotation) -> LabelPosition:
        # A circle draws every label centred (``circular_annotation_labels``).
        return "center" if self.join_origin else self.style(annotation).label_position

    def font_size(self, annotation: ResolvedRegionAnnotation) -> float:
        return annotation_label_font_size(self.style(annotation), self.labels.font_size)

    def extent(self, annotation: ResolvedRegionAnnotation, place_bp: float, mark_offset: float) -> LabelExtent:
        key = (annotation.record_index, annotation.id, place_bp, mark_offset)
        if key not in self.measured:
            self.measured[key] = self.labels.extent(
                annotation, self.style(annotation), self.font_size(annotation), place_bp, mark_offset
            )
        return self.measured[key]

    def row(self, annotation: ResolvedRegionAnnotation) -> int:
        """The row the annotation's label stands on (``label_anchor``)."""

        return annotation.label_anchor(self.label_position(annotation), join_origin=self.join_origin).part.row

    def paint(self, annotation: ResolvedRegionAnnotation, lane_width: float) -> float:
        return self.labels.reach(annotation, self.style(annotation), lane_width)

    def left_out_font_sizes(self, track: ResolvedAnnotationTrack) -> tuple[tuple[int, float], ...]:
        """The font size of each label ``track`` leaves out (``ResolvedAnnotationTrack.without_room_font_sizes``)."""

        by_key = {(item.annotation.record_index, item.annotation.id): item.annotation for item in track.placements}
        return tuple((record, self.font_size(by_key[record, identity])) for record, identity in track.without_room)


class _Letters(NamedTuple):
    """A label's measured letters: where the drawer draws them at any of its record's lane widths (``_hull``), the
    bp they cover on its row (``_letters_at``), and its letters at each of those widths, in their order."""

    extent: LabelExtent
    painted: tuple[tuple[float, float], ...]
    at_widths: tuple[LabelExtent, ...]


def _lane_letters(
    labels: Sequence[tuple[PlacedRegionAnnotation, LabelPlace]],
    lane: int,
    measured: Callable[[PlacedRegionAnnotation, LabelPlace], _Letters],
) -> dict[tuple[int, int], list[tuple[float, float, ResolvedRegionAnnotation]]]:
    """The letters of the labels in ``lane`` where ``measured`` puts them, per record and row, each with its
    annotation, whose mark paints past its lane as ``_lane_rooms`` keeps the marks. They stand past the lane's mark,
    so the room that clears the mark clears them too (OV-2811)."""

    letters: dict[tuple[int, int], list[tuple[float, float, ResolvedRegionAnnotation]]] = {}
    for placement, place in labels:
        if placement.lane == lane:
            label = measured(placement, place)
            if not label.extent.drawn:
                continue
            letters.setdefault((placement.annotation.record_index, place.part.row), []).extend(
                (start, end, placement.annotation) for start, end in label.painted
            )
    return letters


def _label_room(
    label: _Letters,
    annotation: ResolvedRegionAnnotation,
    crossed: Sequence[ResolvedRegionAnnotation],
    *,
    widths: Sequence[float],
    paint_of: Callable[[ResolvedRegionAnnotation, float], float],
    offset: float,
    lanes_between: int,
    rooms_between: float,
) -> float:
    """The room ``annotation``'s label (its ``label`` letters) needs before a lane whose marks or letters it would
    cross (those of ``crossed``), as ``_lane_rooms`` reckons it: at each of its record's lane ``widths`` where its
    letters are drawn, their depth there, ``offset`` on each side of them, its own mark's paint and the furthest
    painting crossed mark's at that width, less ``lanes_between`` lanes of that width and ``rooms_between``. The label
    takes the largest."""

    return max(
        at.depth_px
        + 2.0 * offset
        + paint_of(annotation, width)
        + max(paint_of(other, width) for other in crossed)
        - lanes_between * width
        - rooms_between
        for width, at in zip(widths, label.at_widths)
        if at.drawn
    )


def _lane_rooms(
    placements: Sequence[PlacedRegionAnnotation],
    *,
    record_lengths: Mapping[int, int],
    measures: _TrackMeasures,
    band: _LaneBand,
    join_origin: bool,
    at_least: Sequence[LaneRoom],
    widths: Mapping[int, tuple[float, ...]],
) -> _LaneRooms:
    """The room between lanes for the labels that would cross a mark of another lane on their side (PD-OI-119,
    lane room; LBL-ROOM-BEFORE-MET-LANE), or the letters of its labels (``_lane_letters``), per record and row.

    ``measures`` measures a labelled annotation's letters at its place (``label_anchor``) as the drawer does,
    given how far its mark's painted edge lies from where the lanes start (``_LaneBand.mark_edge_px`` of its lane,
    at each of its record's lane ``widths``, past the rooms of the gaps before its own; labels toward lane 0 stand in Linear, whose letters do not depend on
    it). A label would cross where those letters meet a mark of a lane on its side on the label's row, the next
    lane or one further on. Its room goes into the gap just before that lane, so the lanes between keep their
    spacing beside the label's lane: the letters' depth, the label's ``label_offset`` on each side of them (beside
    its own mark and before the mark it would cross, so the label stands clear of both) and how far the two marks
    paint past their lanes (``_TrackMeasures.paint``, the drawer's ``LabelSizing.reach``: the label's mark, and the
    furthest painting mark it would cross), less what already lies between the label's mark and that lane (a lane
    width for each lane between, and the rooms of the gaps between). The label is measured so at each of its
    record's lane ``widths`` where its letters are drawn, and takes the largest of those rooms. The room of a gap is
    the largest its labels need, and at least its room in ``at_least``. The gaps are
    settled in order from each label's side outward, so the rooms between are known when a gap further on is
    reckoned. Within ``band.limit_px`` a record's row adds each room only if the lanes with it, as the drawer draws
    them (``_LaneBand.drawn_px``), and with the rooms it added before, still end within the limit
    (LBL-ROOM-STOP-AT-CENTRE-OMIT at the centre, LBL-ROOM-CLIP-FIT-OMIT at a clipped set size): a gap takes the
    largest of its labels' rooms that fits, the labels that need more are left without room (they claim no room
    further on), and a later room that fits is added.
    """

    style_of, extent_at, paint_of = measures.style, measures.extent, measures.paint
    lane_count = lane_count_of(placements)
    # Each mark's segments on a row, with its annotation: how far its paint passes its lane (``measures.paint``, at
    # a lane width) toward a label that meets it, as on its own label's side (a mark paints alike on both sides of its
    # lane's centre).
    marks: dict[tuple[int, int, int], list[tuple[float, float, ResolvedRegionAnnotation]]] = {}
    for placement in placements:
        annotation = placement.annotation
        for row in annotation.geometry_rows:
            marks.setdefault((annotation.record_index, row, placement.lane), []).extend(
                (float(start), float(end), annotation) for start, end in annotation.row_segments(row)
            )
    step = -1 if band.labels_above else 1
    every_label = [
        (placement, annotation.label_anchor(measures.label_position(annotation), join_origin=join_origin))
        for placement in placements
        if measures.has_letters(annotation := placement.annotation)
    ]
    labelled = [item for item in every_label if 0 <= item[0].lane + step < lane_count]
    rooms = {(room.record_index, room.row, room.gap): room.px for room in at_least}
    # The gap that parts a lane from the one before it on the labels' side, and each label's letters, measured
    # when its own gap is reckoned (the rooms before it are known then).
    def gap_before(lane: int) -> int:
        return min(lane, lane - step)

    letters: dict[int, _Letters] = {}

    def measured(placement: PlacedRegionAnnotation, place: LabelPlace) -> _Letters:
        """A label's letters where the rooms reckoned so far put its lane."""
        annotation = placement.annotation
        own_gap = gap_before(placement.lane + step)
        before = sum(
            px
            for (record, room_row, room_gap), px in rooms.items()
            if record == annotation.record_index and room_row == place.part.row and room_gap < own_gap
        )
        extents = [
            extent_at(
                annotation, place.bp,
                band.mark_edge_px(placement.lane, width, room_before_px=before, past_lane_edge_px=paint_of(annotation, width)),
            )
            for width in widths[annotation.record_index]
        ]
        # The letters cover the circle as at any of the widths where they are drawn, and take the largest depth they
        # take there.
        hull = _hull(extents)
        extent = extents[0] if hull is None else hull
        length = int(record_lengths[annotation.record_index])
        painted = _letters_at(place.bp, extent.low_bp, extent.high_bp, record_length=length, join_origin=join_origin)
        return _Letters(extent, painted, tuple(extents))

    last_gap = max(lane_count - 2, max((room.gap for room in at_least), default=-1))
    bare = band.drawn_px(lane_count)
    used: dict[tuple[int, int], float] = {}
    left_out: set[tuple[int, str]] = set()
    for gap in range(last_gap + 1) if step > 0 else range(last_gap, -1, -1):
        met_lane = gap + 1 if step > 0 else gap
        met_letters = _lane_letters(every_label, met_lane, measured)
        # The room each label of this gap needs, per record and row.
        needs: dict[tuple[int, int, int], dict[str, float]] = {}
        for index, (placement, place) in enumerate(labelled):
            lanes_on = (met_lane - placement.lane) * step
            annotation = placement.annotation
            # A label left without a room is left out (``without_room``), and claims no room further on.
            if lanes_on < 1 or (annotation.record_index, annotation.id) in left_out:
                continue
            row = place.part.row
            label = letters.get(index) or letters.setdefault(index, measured(*labelled[index]))
            crossed = [
                other
                for start, end, other in (
                    *marks.get((annotation.record_index, row, met_lane), ()),
                    *met_letters.get((annotation.record_index, row), ()),
                )
                if _overlaps(label.painted, ((start, end),))
            ]
            if not (crossed and label.extent.drawn):  # one the drawer leaves out (``LabelExtent.drawn``) claims none
                continue
            own_gap = gap_before(placement.lane + step)
            between = sum(
                rooms.get((annotation.record_index, row, between_gap), 0.0)
                for between_gap in range(min(own_gap, gap), max(own_gap, gap) + 1)
                if between_gap != gap
            )
            room = _label_room(
                label, annotation, crossed, widths=widths[annotation.record_index], paint_of=paint_of,
                offset=float(style_of(annotation).label_offset), lanes_between=lanes_on - 1, rooms_between=between,
            )
            if room > 0.0:
                needs.setdefault((annotation.record_index, row, gap), {})[annotation.id] = room
        limit = band.limit_px
        for key in {*needs, *(key for key in rooms if key[2] == gap)}:
            by_label = needs.get(key, {})
            wanted = [*by_label.values(), *((rooms[key],) if key in rooms else ())]
            if limit is None:
                rooms[key] = max(wanted)
                continue
            # Each room on its own: the gap takes the largest of its rooms that still fits, and the labels that
            # need more are left without room.
            record_row = (key[0], key[1])
            before = used.get(record_row, 0.0)
            fitting = [px for px in wanted if not bare + before + px > limit]
            left_out.update((key[0], identity) for identity, px in by_label.items() if bare + before + px > limit)
            if fitting:
                rooms[key] = max(fitting)
                used[record_row] = before + rooms[key]
            else:
                rooms.pop(key, None)
    return _LaneRooms(
        tuple(LaneRoom(record, row, gap, px) for (record, row, gap), px in sorted(rooms.items())), tuple(sorted(left_out))
    )


def _record_lane_counts(placements: Sequence[PlacedRegionAnnotation]) -> dict[int, int]:
    """How many lanes each record's placements take: each record's lanes fill the slot where the drawers draw them
    (``_LaneBand.width_px``)."""

    return {
        record: lane_count_of(item for item in placements if item.annotation.record_index == record)
        for record in {placement.annotation.record_index for placement in placements}
    }


def _cut_by_clip(
    track: ResolvedAnnotationTrack,
    band: _LaneBand,
    *,
    labelled: bool,
    measures: _TrackMeasures,
) -> ResolvedAnnotationTrack:
    """``track``, leaving out (``ResolvedAnnotationTrack.without_room``) the labels of a clipped track whose letters
    the edge of its set size (``band.limit_px``) would cut, or whose own mark the clip hides on the row the label
    stands on (``_LaneBand.hides_mark``; LBL-ROOM-CLIP-FIT-OMIT); unchanged when it is not clipped or draws no
    labels of its own (not ``labelled``). Measured as the drawer measures them (``_TrackMeasures.extent``) beside their
    mark's paint (``_TrackMeasures.paint``, the drawer's ``LabelSizing.reach``) where its lane lies with the track's rooms
    (``_LaneBand.mark_edge_px``, ``ResolvedAnnotationTrack.room_before``), the letters and their ``label_offset``
    pass the edge the lanes stack from (labels toward lane 0) or the set size (labels away from it)."""

    if not (track.clipped and labelled):
        return track
    limit = float(band.limit_px or 0.0)
    cut: set[tuple[int, str]] = set()
    for placed in _placed_labels(track, band, measures):
        # A label the drawer leaves out for want of a letter (``LabelExtent.drawn``) is counted there, not here.
        if not placed.extent.drawn:
            continue
        placement, edge = placed.placement, placed.edge
        reach = placed.extent.depth_px + float(measures.style(placement.annotation).label_offset)
        hidden = band.hides_mark(placement.lane, placed.width, room_before_px=placed.room, past_lane_edge_px=placed.paint)
        if hidden or ((edge - reach < 0.0) if band.labels_above else (edge + reach > limit)):
            cut.add((placement.annotation.record_index, placement.annotation.id))
    return replace(track, without_room=tuple(sorted({*track.without_room, *cut}))) if cut else track


def _across_centre(
    track: ResolvedAnnotationTrack, band: _LaneBand, labels: LabelSizing, params: AnnotationTrackParams, *, labelled: bool
) -> ResolvedAnnotationTrack:
    """``track``, leaving out (``ResolvedAnnotationTrack.without_room``) the labels whose letters the drawer draws over
    the centre of the circle onto a mark (``LabelSizing.across_centre``; OV-2991): no room keeps such letters clear of
    the marks they meet there, so they are left out as a label whose room would carry its lanes past the centre is
    (LBL-ROOM-STOP-AT-CENTRE-OMIT; in a clipped slot of set size, LBL-ROOM-CLIP-FIT-OMIT), and counted where the
    track's rooms stop (``rooms_stop``). Unchanged where the slot draws no labels of its own (not ``labelled``), its
    letters cannot pass the centre, or its rooms have no stop (``band.limit_px``: a set size with overflow error or
    compress, which the Contract does not name for such a label; OV-2991-B)."""

    if not labelled or labels.across_centre is None or band.limit_px is None:
        return track
    met = labels.across_centre(track, params, band.drawn_extent_px(lane_count_of(track.placements), track.lane_room_px))
    if not met:
        return track
    return replace(track, without_room=tuple(sorted({*track.without_room, *met})), rooms_stop=band.stop)


class _LanePlace(NamedTuple):
    """Where the drawer puts one lane of a record's row: the width of the record's lanes (``_LaneBand.width_px``
    with the track's label room) and the rooms of the gaps before the lane (``ResolvedAnnotationTrack.room_before``)."""

    width: float
    room_before: float


# The places (``_LanePlace``) the packing rounds of a slot drew each lane at, by record index, row and lane
# (``_SlotPacking.settled``).
_LanePlaces = Mapping[tuple[int, int, int], frozenset[_LanePlace]]


def _lane_place(track: ResolvedAnnotationTrack, band: _LaneBand, record_index: int, row: int, lane: int) -> _LanePlace:
    """Where the drawer puts ``lane`` of the record's ``row`` in ``track`` (``_LanePlace``)."""

    width = band.width_px(
        _record_lane_counts(track.placements)[record_index], lane_count_of(track.placements), room_px=track.lane_room_px
    )
    return _LanePlace(width, track.room_before(record_index, row, lane))


class _PlacedLabel(NamedTuple):
    """A label the track draws, measured where the track places it (``_placed_labels``): its letters
    (``LabelExtent``), how far its mark's painted edge lies from where the lanes stack (``edge``), the width of its
    record's lanes, how far its mark paints past its lane and the rooms before its lane on its row."""

    placement: PlacedRegionAnnotation
    extent: LabelExtent
    edge: float
    width: float
    paint: float
    room: float


def _placed_labels(track: ResolvedAnnotationTrack, band: _LaneBand, measures: _TrackMeasures) -> Iterable[_PlacedLabel]:
    """Each label ``track`` draws (``ResolvedAnnotationTrack.draws_label``), measured as the drawer measures it
    (``_TrackMeasures.extent``) beside its mark's paint (``_TrackMeasures.paint``, the drawer's ``LabelSizing.reach``)
    where its lane lies with the track's rooms (``_LaneBand.mark_edge_px``, ``ResolvedAnnotationTrack.room_before``)
    and is as wide as the drawer draws it (``_LaneBand.width_px`` with the track's label room)."""

    for placement in track.placements:
        annotation = placement.annotation
        if not track.draws_label(annotation):
            continue
        place = annotation.label_anchor(measures.label_position(annotation), join_origin=track.join_origin)
        width, room = _lane_place(track, band, annotation.record_index, place.part.row, placement.lane)
        paint = measures.paint(annotation, width)
        edge = band.mark_edge_px(placement.lane, width, room_before_px=room, past_lane_edge_px=paint)
        yield _PlacedLabel(placement, measures.extent(annotation, place.bp, edge), edge, width, paint, room)


def _lettered(track: ResolvedAnnotationTrack, band: _LaneBand, measures: _TrackMeasures) -> frozenset[tuple[int, str]]:
    """The labels ``track`` draws that keep a letter where it places them (``_placed_labels``, ``LabelExtent.drawn``),
    by record index and annotation id."""

    return frozenset(
        (placed.placement.annotation.record_index, placed.placement.annotation.id)
        for placed in _placed_labels(track, band, measures)
        if placed.extent.drawn
    )


def _keeping_no_letter(track: ResolvedAnnotationTrack, band: _LaneBand, measures: _TrackMeasures) -> frozenset[tuple[int, str]]:
    """The labels ``track`` draws that keep no letter where it places them (``_placed_labels``,
    ``LabelExtent.drawn``), by record index and annotation id."""

    return frozenset(
        (placed.placement.annotation.record_index, placed.placement.annotation.id)
        for placed in _placed_labels(track, band, measures)
        if not placed.extent.drawn
    )


def _hull(extents: Iterable[LabelExtent]) -> LabelExtent | None:
    """The letters a label paints at any of its places (``extents``) where the drawer draws them
    (``LabelExtent.drawn``), or None where it draws them at none."""

    drawn = [extent for extent in extents if extent.drawn]
    if not drawn:
        return None
    return LabelExtent(
        min(extent.low_bp for extent in drawn), max(extent.high_bp for extent in drawn),
        max(extent.depth_px for extent in drawn), True,
    )


def _letters_in_lanes(
    band: _LaneBand,
    measures: _TrackMeasures,
    *,
    shows: Callable[[ResolvedRegionAnnotation], bool],
    by_lane: bool,
    drawn_at: _LanePlaces,
) -> LettersInLane:
    """A label's letters beside its mark in the lane the lane packer tries (``LettersInLane``), as the drawer
    measures them (``_TrackMeasures.extent``): where ``band`` puts that lane as the last of the lanes so far, before
    the label rooms between the lanes are known (``_lane_rooms``), or, when not ``by_lane``, beside a mark in lane
    0 whatever the lane; and at every place an earlier packing round of the slot drew that lane (``drawn_at``, with
    its rooms and the width the drawer draws its lanes at, OV-2849). The lane keeps the letters at all those places:
    a room moves the lanes after it further from the axis, and in an inside slot the letters there span more of the
    circle; a slot of set width shares that width among its lanes.

    Only the labels ``shows`` gives letters have them, and only at the places where the drawer draws them
    (``LabelExtent.drawn``). A Linear label's letters do not depend on its lane (``LabelSizing.extent``). An inside
    Circular slot whose lanes, measured by lane, would pass the centre (``_LaneBand.passes``) is measured in lane 0
    before its lanes are drawn: near the centre a label's letters cover most of the circle, so each such label would
    take a lane of its own, and the lanes would leave no room before the centre."""

    def letters(annotation: ResolvedRegionAnnotation, place_bp: float, lane: int) -> LabelExtent | None:
        if not shows(annotation):
            return None
        at = lane if by_lane else 0
        first = band.width_px(at + 1, at + 1)
        edges = {
            band.mark_edge_px(at, first, room_before_px=0.0, past_lane_edge_px=measures.paint(annotation, first)),
            *(
                band.mark_edge_px(
                    lane, place.width, room_before_px=place.room_before,
                    past_lane_edge_px=measures.paint(annotation, place.width),
                )
                for place in drawn_at.get((annotation.record_index, measures.row(annotation), lane), ())
            ),
        }
        return _hull(measures.extent(annotation, place_bp, edge) for edge in sorted(edges))

    return letters


def _lane_band(
    track_params: AnnotationTrackParams,
    *,
    available_extent_px: float | None,
    room_bounds: RoomBounds,
    labels_above: bool,
) -> _LaneBand:
    """The band a slot's lanes and rooms are packed in (``_LaneBand``), stopped at the nearer of its bounds: a slot
    whose size is set and whose overflow is clip at that size (LBL-ROOM-CLIP-FIT-OMIT), and an inside Circular slot
    at the centre of the circle (``room_bounds.limit_px``, LBL-ROOM-STOP-AT-CENTRE-OMIT), which also bounds a clipped
    inside slot whose set width reaches past the centre. A slot whose size is set with overflow error or compress
    has no stop: that setting decides on the room it needs. A span over other slots (``to_slot``) takes its extent
    from them, not from the slot's own size."""

    centre = room_bounds.limit_px
    limit: float | None = centre
    stop: RoomStop = "centre"
    set_px: float | None = None
    if available_extent_px is not None:
        clips = track_params.overflow == "clip" and track_params.to_slot is None
        if clips and (centre is None or float(available_extent_px) <= centre):
            limit, stop = float(available_extent_px), "clip"
        elif not clips:
            limit = None
        if track_params.to_slot is None:
            set_px = float(available_extent_px) if centre is None else min(float(available_extent_px), centre)
    return _LaneBand(
        LANE_EXTENT_PX, float(track_params.lane_gap_px), float(track_params.padding_px), limit, stop, labels_above, set_px
    )


class _SetSizeFit(NamedTuple):
    """How a slot's lanes and label room take its set size (``_set_size_fit``): the extent they take, the gap between
    the lanes, and whether the drawer draws them (False where it raises instead)."""

    required_px: float
    lane_gap_px: float
    draws: bool


def _set_size_fit(params: AnnotationTrackParams, band: _LaneBand, lane_count: int, room_px: float) -> _SetSizeFit:
    """How ``lane_count`` lanes and the label room ``room_px`` take the extent the drawer gives a slot of set size
    (``_LaneBand.set_px``: an inside Circular slot's set size stops at the centre, OV-2990): a slot whose overflow is
    error is not drawn where they pass it, and one whose overflow is compress narrows the gaps between the lanes to fit
    it and is not drawn where they pass it with no gap. A slot without a set size, or a span over other slots
    (``to_slot``), which takes its extent from them, is drawn."""

    lane_gap = float(params.lane_gap_px)
    required = band.bare_px(lane_count) + room_px
    available_px = band.set_px
    if available_px is None or not required > available_px:
        return _SetSizeFit(required, lane_gap, True)
    if params.overflow == "error":
        return _SetSizeFit(required, lane_gap, False)
    if params.overflow == "compress" and lane_count > 1:
        lane_gap = max(
            0.0,
            min(lane_gap, (available_px - 2.0 * params.padding_px - lane_count * band.lane_extent_px - room_px) / (lane_count - 1)),
        )
        required = 2.0 * params.padding_px + lane_count * band.lane_extent_px + (lane_count - 1) * lane_gap + room_px
        return _SetSizeFit(required, lane_gap, not required > available_px)
    return _SetSizeFit(required, lane_gap, True)


class _SlotPacking(NamedTuple):
    """One slot's annotations and what packs them (``layout_annotation_track``): its lanes, their rooms and the
    labels left out, for a set of labels left out for want of a letter (``settled``)."""

    slot_id: str
    set_id: str
    selected: tuple[ResolvedRegionAnnotation, ...]
    record_lengths: Mapping[int, int]
    params: AnnotationTrackParams
    labels: LabelSizing
    join_origin: bool
    band: _LaneBand
    room_bounds: RoomBounds
    # Whether the slot draws labels of its own and gives them room (not a span over other slots).
    labelled: bool
    measured: dict[tuple[int, str, float, float], LabelExtent]

    def padding(self, annotation: ResolvedRegionAnnotation) -> float:
        scale = float(self.labels.bp_per_px.get(annotation.record_index, 1.0))
        style = effective_annotation_style(annotation, self.params)
        cap = style.stroke_width * 2.0 if annotation.mark == "bracket" else 0.0
        return (self.params.padding_px + cap) * scale

    def packed(self, measures: _TrackMeasures, drawn_at: _LanePlaces) -> tuple[PlacedRegionAnnotation, ...]:
        """The lanes, each keeping its labels' letters where the lane packer tries it and where earlier rounds drew it
        (``_letters_in_lanes``)."""

        highlight_track = self.params.marks == ("highlight",)

        def packed_by(by_lane: bool) -> tuple[PlacedRegionAnnotation, ...]:
            return assign_annotation_lanes(
                self.selected,
                record_lengths=self.record_lengths,
                letters_in_lane=_letters_in_lanes(
                    self.band, measures, by_lane=by_lane, drawn_at=drawn_at,
                    shows=lambda annotation: self.params.show_labels and measures.has_letters(annotation),
                ),
                padding_bp=self.padding,
                overflow="compress" if highlight_track else self.params.overflow,
                join_origin=self.join_origin,
                label_position=measures.label_position,
            )

        # A Circular label's letters depend on its lane; an inside slot whose lanes would pass the centre measured
        # so measures them all in its outer lane (``_letters_in_lanes``).
        placements = packed_by(self.join_origin)
        if self.band.stop == "centre" and self.band.passes(lane_count_of(placements)) and not highlight_track:
            placements = packed_by(False)
        if highlight_track:
            placements = tuple(replace(item, lane=0) for item in placements)
        return placements

    def resolved(self, placements: tuple[PlacedRegionAnnotation, ...], measures: _TrackMeasures) -> ResolvedAnnotationTrack:
        """The track of ``placements`` with their rooms (``_lane_rooms``), its extent not yet reckoned."""

        lane_count = lane_count_of(placements)

        def rooms(at_least: Sequence[LaneRoom], drawn_rooms_px: Sequence[float]) -> _LaneRooms:
            # Each record's lanes are measured at the extent they take and at the width the drawer draws them with
            # each label room of ``drawn_rooms_px`` (``_LaneBand.width_px``).
            widths = {
                record: tuple(dict.fromkeys((
                    self.band.width_px(lanes, lane_count),
                    *(self.band.width_px(lanes, lane_count, room_px=room) for room in drawn_rooms_px),
                )))
                for record, lanes in _record_lane_counts(placements).items()
            }
            if not self.labelled:
                return _LaneRooms(())
            return _lane_rooms(
                placements, record_lengths=self.record_lengths, measures=measures, band=self.band,
                join_origin=self.join_origin, at_least=at_least, widths=widths,
            )

        lane_rooms = rooms(self.room_bounds.at_least, ())
        # A slot of set size draws its lanes in that size less the slot's label room, whether they and the rooms fit
        # it or not, not at the extent they take (OV-2849): where they do not fit, it draws them narrower, and a label
        # that keeps no letter at that extent may keep letters where its lane is drawn (OV-2989). The rooms are
        # measured again, each round also where the drawer draws the lanes with the room of every round before, while
        # a round's room is larger than every room before. A round keeps the rooms of the one before, so the room of a
        # slot that does not clip only grows, and its last round measured the lanes where the drawer draws them. A
        # clipped slot keeps only the rooms that fit its set size: a room that grows can leave no space for one further
        # on, which is then dropped, and the slot can end at a room no round measured. The rooms are bounded by the
        # labels' letters, and finitely many numbers lie below that bound, so the rounds end. A slot the drawer does not
        # draw (``_set_size_fit``) is not measured again: its rooms stay where its lanes keep their extent, as at the
        # size it then names.
        drawn_rooms: list[float] = []
        lane_count = lane_count_of(placements)
        while (
            self.band.set_px is not None
            and (room := _largest_room(lane_rooms.rooms)) > max(drawn_rooms, default=-1.0)
            and _set_size_fit(self.params, self.band, lane_count, room).draws
        ):
            drawn_rooms.append(room)
            lane_rooms = rooms((*self.room_bounds.at_least, *lane_rooms.rooms), drawn_rooms)
        return ResolvedAnnotationTrack(
            slot_id=self.slot_id, set_id=self.set_id, placements=placements, required_extent_px=0.0,
            label_font_size=float(self.labels.font_size), join_origin=self.join_origin, lane_rooms=lane_rooms.rooms,
            labels_above=self.labels.labels_above, without_room=lane_rooms.without_room,
            rooms_stop=None if self.band.limit_px is None else self.band.stop,
            without_letters=tuple(sorted(measures.without_letters)),
        )

    def settled(self, without_letters: frozenset[tuple[int, str]]) -> tuple[_TrackMeasures, ResolvedAnnotationTrack]:
        """The lanes and rooms with the labels of ``without_letters`` left out (``_TrackMeasures``), each lane keeping
        its labels' letters where the drawer draws it (OV-2849).

        A first round packs the lanes measured where the lane packer tries them. While the drawer draws two labels
        of one lane over each other (``LabelSizing.meet``), another round measures that lane also where the drawer
        draws it in the round before (``_lane_place``): a lane measured where it is drawn keeps each label's letters
        there, so its labels do not meet at their places. A round whose labels do not meet ends the rounds. A way of
        putting the labels in lanes that comes back has the rounds measure each lane, from then on, at every place a
        round drew it where its labels met, and a round so measured that adds no such place ends them too. The rounds end:
        a round's rooms, and so where it draws its lanes, follow from its lanes alone (``resolved``), so the places
        come from the finitely many ways to put the labels in lanes, and so do the rounds before one comes back.
        Labels in explicit lanes may meet (``assign_annotation_lanes``)."""

        measures = _TrackMeasures(self.params, self.labels, self.join_origin, without_letters, self.measured)
        drawn_at: dict[tuple[int, int, int], frozenset[_LanePlace]] = {}
        # Every place a round drew a lane where its labels met, kept from when a way of putting them in lanes comes
        # back (``every``).
        met_at: dict[tuple[int, int, int], frozenset[_LanePlace]] = {}
        every = False
        # The track of each way of putting the labels in lanes, by record index, annotation id and lane.
        tracks: dict[tuple[tuple[int, str, int], ...], ResolvedAnnotationTrack] = {}
        while True:
            placements = self.packed(measures, drawn_at)
            lanes = tuple((item.annotation.record_index, item.annotation.id, item.lane) for item in placements)
            measured_at_every = every
            every = every or lanes in tracks
            if lanes not in tracks:
                tracks[lanes] = self.resolved(placements, measures)
            track = tracks[lanes]
            if self.labels.meet is None or not self.labelled:
                return measures, track
            lane_count = lane_count_of(placements)
            met = {
                key: _lane_place(track, self.band, *key)
                for key in self.labels.meet(
                    replace(track, clipped=self.band.stop == "clip" and self.band.passes(lane_count)),
                    self.params,
                    self.band.drawn_extent_px(lane_count, track.lane_room_px),
                )
            }
            if not met or (
                measured_at_every and all(place in met_at.get(key, frozenset()) for key, place in met.items())
            ):
                return measures, track
            for key, place in met.items():
                met_at[key] = met_at.get(key, frozenset()) | {place}
            drawn_at = met_at if every else {key: frozenset((place,)) for key, place in met.items()}

    def without_letters(self) -> tuple[_TrackMeasures, ResolvedAnnotationTrack]:
        """The lanes and rooms with the labels left out that keep no letter where the lanes put them
        (RGN-CUT-EMPTY-OMIT): such a label is not drawn, and its text moves no band (OV-2805), unless that costs
        another label its letters.

        Each label that keeps no letter is packed again as if it had no text, one at a time, and kept so if every
        label drawn with letters before still is; else it keeps its text in the packing and is left out where its
        text put its mark, so that what was drawn stays drawn. Each label is tried once, so this ends. Then each label
        packed without text that keeps a letter where the packing puts its mark is given back its text, one at a time,
        and kept so if it is then drawn with letters and every label drawn with letters before still is (OV-2850): a
        label left out while it stood in a lane further on is drawn where the packing that left out others moved it.
        Each label given back shrinks the labels packed without text, so this ends too."""

        measures, track = self.settled(frozenset())
        if not self.labelled:
            return measures, track
        tried: set[tuple[int, str]] = set()
        while more := sorted(_keeping_no_letter(track, self.band, measures) - tried):
            tried.add(more[0])
            again = self.settled(measures.without_letters | {more[0]})
            if _lettered(track, self.band, measures) <= _lettered(again[1], self.band, again[0]):
                measures, track = again
            else:
                track = replace(track, without_letters=tuple(sorted({*track.without_letters, more[0]})))
        given = True
        while given:
            given = False
            for key in sorted(measures.without_letters):
                if not _keeps_a_letter(track, self.band, measures, key):
                    continue
                again = self.settled(measures.without_letters - {key})
                drawn = _lettered(again[1], self.band, again[0])
                if key in drawn and _lettered(track, self.band, measures) <= drawn:
                    (measures, track), given = again, True
        # The labels that keep their text in the packing but no letter where it puts them are left out as well.
        left_out = {*track.without_letters, *_keeping_no_letter(track, self.band, measures)}
        return measures, replace(track, without_letters=tuple(sorted(left_out)))


def _keeps_a_letter(
    track: ResolvedAnnotationTrack, band: _LaneBand, measures: _TrackMeasures, key: tuple[int, str]
) -> bool:
    """Whether the label ``key`` (record index, annotation id), left out for want of a letter, keeps one where
    ``track`` places its mark (``_placed_labels``)."""

    shown = replace(track, without_letters=tuple(item for item in track.without_letters if item != key))
    return any(
        (placed.placement.annotation.record_index, placed.placement.annotation.id) == key and placed.extent.drawn
        for placed in _placed_labels(shown, band, measures)
    )


def layout_annotation_track(
    slot_id: str,
    set_id: str,
    annotations: Sequence[ResolvedRegionAnnotation],
    *,
    record_lengths: Mapping[int, int],
    params: AnnotationTrackParams | None = None,
    available_extent_px: float | None = None,
    labels: LabelSizing,
    join_origin: bool,
    room_bounds: RoomBounds = RoomBounds(),
    checks_set_size: bool = True,
) -> ResolvedAnnotationTrack:
    """Pack one set and calculate its radial/vertical required extent.

    ``labels`` measures the labels (``LabelSizing``); a lane keeps a label's letters as the drawer measures them
    beside the mark in it, where the drawer draws that lane (``_SlotPacking.settled``). ``join_origin`` is True for a
    Circular drawing (``mark_parts``). Lanes lie where ``_LaneBand`` puts them (``lane_edge_px``, as the drawers
    place them) and part by the room of a label that would cross the mark of a lane on its side, within
    ``room_bounds`` (``_lane_rooms``); the track's required extent grows by the largest room a record's row takes
    (``ResolvedAnnotationTrack.lane_room_px``). A label that keeps no letter where its lane is drawn is left out
    (``_SlotPacking.without_letters``). A clipped slot of set size (``available_extent_px``) is bounded by
    that size (``_lane_band``); whether it clips is decided by its lanes alone, before the rooms, so the packer
    measures its lanes where the drawer draws them, and it leaves out the labels its edge would cut
    (``_cut_by_clip``). A span over other slots (``to_slot``, ``cover_anchor``) gets no room: where it covers a slot
    its labels stand outside the slots it covers, and whether it covers one is known only after the lanes are packed.

    A slot whose lanes and rooms do not fit the extent the drawer gives its set size raises (``_set_size_fit``),
    naming the extent they need, unless not ``checks_set_size``: a packing made where the slot is not drawn (the
    Circular planner's at the axis circle, packed again where the slot is drawn) does not raise, so the size named is
    measured where the slot is drawn (OV-3024)."""

    track_params = params or AnnotationTrackParams(set_id=set_id)
    if track_params.set_id != set_id:
        raise ValidationError(f"Annotation slot {slot_id!r} requests set {track_params.set_id!r}, not {set_id!r}.")
    band = _lane_band(
        track_params, available_extent_px=available_extent_px, room_bounds=room_bounds,
        labels_above=labels.labels_above,
    )
    packing = _SlotPacking(
        str(slot_id),
        str(set_id),
        tuple(
            item
            for item in annotations
            if item.set_id == set_id and (track_params.marks is None or item.mark in track_params.marks)
        ),
        record_lengths,
        track_params,
        labels,
        join_origin,
        band,
        room_bounds,
        track_params.show_labels and not (track_params.to_slot is not None or track_params.cover_anchor),
        {},
    )
    measures, track = packing.without_letters()
    placements = track.placements
    room = track.lane_room_px
    lane_count = lane_count_of(placements)
    fit = _set_size_fit(track_params, band, lane_count, room)
    if not fit.draws and checks_set_size and band.set_px is not None:
        # Lanes that need more than the room before the centre of an inside Circular slot fit no set size: that room
        # bounds them, whatever the size. Otherwise the size named (rounded up, so it holds the lanes) draws.
        centre = room_bounds.limit_px
        past_centre = centre is not None and fit.required_px > centre
        bound = f"{centre:.2f}px lies before the centre of the circle" if past_centre else f"{band.set_px:.2f}px is available"
        if track_params.overflow == "error":
            raise ValidationError(
                f"Annotation slot {slot_id!r} requires {math.ceil(fit.required_px * 100.0) / 100.0:.2f}px but only {bound}."
            )
        where = f"the {centre:.2f}px before the centre of the circle" if past_centre else f"{band.set_px:.2f}px"
        raise ValidationError(f"Annotation slot {slot_id!r} cannot fit {lane_count} lanes in {where}.")
    # Its rooms fit within the set size (``_lane_rooms``), so only lanes that pass it alone are clipped.
    clipped = (
        track_params.overflow == "clip"
        and track_params.to_slot is None
        and available_extent_px is not None
        and band.bare_px(lane_count) + room > available_extent_px
        and band.passes(lane_count)
    )
    track = replace(track, required_extent_px=float(fit.required_px), clipped=clipped, lane_gap_px=fit.lane_gap_px)
    track = _cut_by_clip(track, band, labelled=packing.labelled, measures=measures)
    track = _across_centre(track, band, labels, track_params, labelled=packing.labelled)
    return replace(track, without_room_font_sizes=measures.left_out_font_sizes(track))


__all__ = [
    "LANE_EXTENT_PX",
    "LabelExtent",
    "LabelMeasure",
    "LabelSizing",
    "LabelsAcrossCentre",
    "LaneRoom",
    "LettersInLane",
    "MarkReachMeasure",
    "PlacedRegionAnnotation",
    "ResolvedAnnotationTrack",
    "RoomBounds",
    "RoomStop",
    "assign_annotation_lanes",
    "lane_count_of",
    "lane_edge_px",
    "lane_width_px",
    "layout_annotation_track",
]
