"""Deterministic lane packing shared by circular and linear annotations."""

from __future__ import annotations

import math
from dataclasses import dataclass, field, replace
from typing import Callable, Mapping, NamedTuple, Sequence

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
    mark_geometry: object | None = None
    label_geometry: object | None = None
    occupied_bounds: tuple[tuple[float, float], ...] = ()


class LaneRoom(NamedTuple):
    """Room between lane ``gap`` and lane ``gap + 1`` on one row of one record for the labels that would cross
    the mark of lane ``gap + 1`` (labels toward it) or of lane ``gap`` (labels toward lane 0) (PD-OI-119, lane
    room): their letters, their ``label_offset`` and half of each of the two marks' strokes, less the lanes and
    rooms between a label's mark and that lane."""

    record_index: int
    row: int
    gap: int
    px: float


class RoomBounds(NamedTuple):
    """What bounds a track's lane rooms (``layout_annotation_track``): the room each gap keeps at least
    (``at_least``, an earlier packing's), and the most the track may take with its rooms (``limit_px``; None for no
    limit): an inside Circular slot without a set width, the radius of the edge its lanes stack from."""

    at_least: tuple[LaneRoom, ...] = ()
    limit_px: float | None = None


class LabelExtent(NamedTuple):
    """A label's painted letters at its place, as the drawer measures and draws them: the bp they reach along the
    region from the place (``low_bp`` before it, ``high_bp`` after it, both offsets from the place) and the pixels
    they take across the lanes."""

    low_bp: float
    high_bp: float
    depth_px: float


# A label's letters (``LabelExtent``) from its annotation, its effective style, its font size, the bp it stands at
# and how far its mark's edge lies from the edge the slot's lanes stack from.
LabelMeasure = Callable[[ResolvedRegionAnnotation, RegionAnnotationStyle, float, float, float], LabelExtent]

# A labelled annotation's letters standing at a bp beside its mark in a lane, as the drawer measures them there
# (``LabelMeasure``); None where the track draws no label for it.
LettersInLane = Callable[[ResolvedRegionAnnotation, float, int], "LabelExtent | None"]


@dataclass(frozen=True)
class LabelSizing:
    """How the packer measures a track's labels.

    ``font_size`` is the font size of the drawing's feature labels, which a label whose style sets none takes;
    ``bp_per_px`` the bases one pixel spans on each record (1 when not given). ``extent`` is the drawer's own
    measure of a label's letters at its place (``LabelMeasure``). ``labels_above`` says that the labels stand
    toward lane - 1 (Linear labels above their marks); else they stand toward lane + 1 (Linear labels below, and
    Circular, away from the axis on either side).
    """

    font_size: float
    extent: LabelMeasure
    bp_per_px: Mapping[int, float] = field(default_factory=dict)
    labels_above: bool = False


@dataclass(frozen=True)
class ResolvedAnnotationTrack:
    slot_id: str
    set_id: str
    placements: tuple[PlacedRegionAnnotation, ...]
    required_extent_px: float
    occupied_bounds: tuple[tuple[float, float], ...]
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
    # How many labels would cross another lane's mark and were left without their room: an inside Circular slot
    # adds no room past the centre (``RoomBounds.limit_px``).
    labels_without_room: int = 0

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
    drawer measures them beside the annotation's mark in that lane (``letters_in_lane``, LBL-RESERVE-INK; no
    letters when not given).

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
            annotation.label_stretch(
                "center" if label_position is None else label_position(annotation), join_origin=join_origin
            )
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
                row: _expanded_intervals(
                    annotation, row, record_length=length, letters=letters, padding_bp=annotation_padding
                )
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
        placements.append(
            PlacedRegionAnnotation(
                annotation=annotation,
                lane=lane,
                occupied_bounds=tuple(bound for bounds in rows.values() for bound in bounds),
            )
        )
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
    """The rooms between lanes (``LaneRoom``), and which labels need each one: the annotation ids by record, row
    and gap."""

    rooms: tuple[LaneRoom, ...]
    labels: Mapping[tuple[int, int, int], frozenset[str]]


def _lane_rooms(
    placements: Sequence[PlacedRegionAnnotation],
    *,
    record_lengths: Mapping[int, int],
    style_of: Callable[[ResolvedRegionAnnotation], RegionAnnotationStyle],
    extent_at: Callable[[ResolvedRegionAnnotation, float, float], LabelExtent],
    mark_offset_px: Callable[[PlacedRegionAnnotation, float], float],
    lane_pitch_px: float,
    labels_above: bool,
    join_origin: bool,
    label_position: Callable[[ResolvedRegionAnnotation], LabelPosition],
    at_least: Sequence[LaneRoom],
) -> _LaneRooms:
    """The room between lanes for the labels that would cross a mark of another lane on their side (PD-OI-119,
    lane room; LBL-ROOM-BEFORE-MET-LANE), per record and row.

    ``extent_at`` measures a labelled annotation's letters at its place (``label_anchor``) as the drawer does,
    given how far its mark's painted edge lies from where the lanes start (``mark_offset_px`` of its placement and
    the rooms of the gaps before its own; labels toward lane 0 stand in Linear, whose letters do not depend on
    it). A label would cross where those letters meet a mark of a lane on its side on the label's row, the next
    lane or one further on. Its room goes into the gap just before that lane, so the lanes between keep their
    spacing beside the label's lane: the letters' depth, the label's ``label_offset`` and half of each of the two
    strokes (the label's mark and the widest mark it would cross), less what already lies between the label's mark
    and that lane (``lane_pitch_px`` for each lane between, and the rooms of the gaps between). The room of a gap
    is the largest its labels need, and at least its room in ``at_least``. The gaps are settled in order from each
    label's side outward, so the rooms between are known when a gap further on is reckoned.
    """

    marks: dict[tuple[int, int, int], list[tuple[float, float, float]]] = {}
    for placement in placements:
        annotation = placement.annotation
        stroke = float(style_of(annotation).stroke_width)
        for row in annotation.geometry_rows:
            marks.setdefault((annotation.record_index, row, placement.lane), []).extend(
                (float(start), float(end), stroke) for start, end in annotation.row_segments(row)
            )
    lane_count = max((placement.lane for placement in placements), default=-1) + 1
    step = -1 if labels_above else 1
    labelled = [
        (placement, annotation.label_anchor(label_position(annotation), join_origin=join_origin))
        for placement in placements
        if (annotation := placement.annotation).label and 0 <= placement.lane + step < lane_count
    ]
    rooms = {(room.record_index, room.row, room.gap): room.px for room in at_least}
    needed_by: dict[tuple[int, int, int], set[str]] = {}
    # The gap that parts a lane from the one before it on the labels' side, and each label's letters, measured
    # when its own gap is reckoned (the rooms before it are known then).
    def gap_before(lane: int) -> int:
        return min(lane, lane - step)

    letters: dict[int, tuple[LabelExtent, tuple[tuple[float, float], ...]]] = {}

    def letters_of(index: int) -> tuple[LabelExtent, tuple[tuple[float, float], ...]]:
        if index not in letters:
            placement, place = labelled[index]
            annotation = placement.annotation
            own_gap = gap_before(placement.lane + step)
            before = sum(
                px
                for (record, room_row, room_gap), px in rooms.items()
                if record == annotation.record_index and room_row == place.part.row and room_gap < own_gap
            )
            extent = extent_at(annotation, place.bp, mark_offset_px(placement, before))
            letters[index] = extent, _letters_at(
                place.bp,
                extent.low_bp,
                extent.high_bp,
                record_length=int(record_lengths[annotation.record_index]),
                join_origin=join_origin,
            )
        return letters[index]

    for gap in range(lane_count - 1) if step > 0 else range(lane_count - 2, -1, -1):
        met_lane = gap + 1 if step > 0 else gap
        for index, (placement, place) in enumerate(labelled):
            lanes_on = (met_lane - placement.lane) * step
            if lanes_on < 1:
                continue
            annotation = placement.annotation
            row = place.part.row
            extent, painted = letters_of(index)
            crossed = [
                stroke
                for start, end, stroke in marks.get((annotation.record_index, row, met_lane), ())
                if _overlaps(painted, ((start, end),))
            ]
            if not crossed:
                continue
            style = style_of(annotation)
            own_gap = gap_before(placement.lane + step)
            between = sum(
                rooms.get((annotation.record_index, row, between_gap), 0.0)
                for between_gap in range(min(own_gap, gap), max(own_gap, gap) + 1)
                if between_gap != gap
            )
            room = (
                extent.depth_px
                + float(style.label_offset)
                + 0.5 * (float(style.stroke_width) + max(crossed))
                - (lanes_on - 1) * float(lane_pitch_px)
                - between
            )
            if room > 0.0:
                key = (annotation.record_index, row, gap)
                rooms[key] = max(rooms.get(key, 0.0), room)
                needed_by.setdefault(key, set()).add(annotation.id)
    return _LaneRooms(
        tuple(LaneRoom(record, row, gap, px) for (record, row, gap), px in sorted(rooms.items())),
        {key: frozenset(ids) for key, ids in needed_by.items()},
    )


def _rooms_within(
    lane_rooms: _LaneRooms, *, bare_extent_px: float, limit_px: float | None
) -> tuple[tuple[LaneRoom, ...], int]:
    """The rooms a slot keeps when, with them, it may take at most ``limit_px`` (``bare_extent_px`` without them;
    all of them without a limit): gap by gap from the edge the lanes stack from, each record's row adds rooms until
    one would carry it past the limit, and adds none after it (LBL-ROOM-STOP-AT-CENTRE-WARN). With how many labels
    lose a room they need."""

    if limit_px is None:
        return lane_rooms.rooms, 0
    kept: list[LaneRoom] = []
    stopped: set[tuple[int, int]] = set()
    used: dict[tuple[int, int], float] = {}
    left_out: set[tuple[int, str]] = set()
    for room in sorted(lane_rooms.rooms, key=lambda item: (item.record_index, item.row, item.gap)):
        key = (room.record_index, room.row)
        if key not in stopped and bare_extent_px + used.get(key, 0.0) + room.px <= limit_px:
            kept.append(room)
            used[key] = used.get(key, 0.0) + room.px
            continue
        stopped.add(key)
        left_out.update(
            (room.record_index, identity) for identity in lane_rooms.labels.get((room.record_index, room.row, room.gap), ())
        )
    return tuple(kept), len(left_out)


def layout_annotation_track(
    slot_id: str,
    set_id: str,
    annotations: Sequence[ResolvedRegionAnnotation],
    *,
    record_lengths: Mapping[int, int],
    params: AnnotationTrackParams | None = None,
    available_extent_px: float | None = None,
    lane_extent_px: float = 14.0,
    labels: LabelSizing,
    join_origin: bool,
    room_bounds: RoomBounds = RoomBounds(),
) -> ResolvedAnnotationTrack:
    """Pack one set and calculate its radial/vertical required extent.

    ``labels`` measures the labels (``LabelSizing``); a lane keeps a label's letters as the drawer measures them
    beside the mark in it (LBL-RESERVE-INK). ``join_origin`` is True for a Circular drawing (``mark_parts``). Lanes
    part by the room of a label that would cross the mark of a lane on its side (``_lane_rooms``) within
    ``room_bounds`` (``_rooms_within``), and the track's required extent grows by the largest room a record's row
    takes (``ResolvedAnnotationTrack.lane_room_px``). A span over other slots (``to_slot``, ``cover_anchor``) gets
    none: where it covers a slot its labels stand outside the slots it covers, and whether it covers one on a record
    is known only after the lanes are packed."""

    track_params = params or AnnotationTrackParams(set_id=set_id)
    if track_params.set_id != set_id:
        raise ValidationError(f"Annotation slot {slot_id!r} requests set {track_params.set_id!r}, not {set_id!r}.")
    marks = track_params.marks
    selected = tuple(item for item in annotations if item.set_id == set_id and (marks is None or item.mark in marks))

    def style_of(annotation: ResolvedRegionAnnotation) -> RegionAnnotationStyle:
        return effective_annotation_style(annotation, track_params)

    def extent_at(annotation: ResolvedRegionAnnotation, place_bp: float, mark_offset: float) -> LabelExtent:
        style = style_of(annotation)
        size = annotation_label_font_size(style, labels.font_size)
        return labels.extent(annotation, style, size, place_bp, mark_offset)

    def lane_pitch(lane_count: int) -> float:
        # One lane of a slot sized by its ``lane_count`` lanes: the lane gaps shared out among them, as drawn.
        return float(lane_extent_px) + float(track_params.lane_gap_px) * (lane_count - 1) / lane_count

    def mark_offset(annotation: ResolvedRegionAnnotation, lane: int, lane_count: int, room_before: float) -> float:
        # The mark's painted edge on its label's side, from where the lanes stack, past the rooms before its lane.
        half_stroke = 0.5 * float(style_of(annotation).stroke_width)
        if labels.labels_above:
            return float(track_params.padding_px) + lane * lane_pitch(lane_count) + room_before - half_stroke
        return float(track_params.padding_px) + (lane + 1) * lane_pitch(lane_count) + room_before + half_stroke

    def letters_in_lane(annotation: ResolvedRegionAnnotation, place_bp: float, lane: int) -> LabelExtent | None:
        # The letters beside the mark in a lane the packer tries, of a slot with as many lanes as that one needs.
        if not track_params.show_labels or not annotation.label:
            return None
        return extent_at(annotation, place_bp, mark_offset(annotation, lane, lane + 1, 0.0))

    def padding(annotation: ResolvedRegionAnnotation) -> float:
        scale = float(labels.bp_per_px.get(annotation.record_index, 1.0))
        stroke_width = effective_annotation_style(annotation, track_params).stroke_width
        cap = stroke_width * 2.0 if annotation.mark == "bracket" else 0.0
        return (track_params.padding_px + cap) * scale

    def label_position(annotation: ResolvedRegionAnnotation) -> LabelPosition:
        # A circle draws every label centred (``circular_annotation_labels``).
        return "center" if join_origin else effective_annotation_style(annotation, track_params).label_position

    highlight_track = track_params.marks == ("highlight",)
    placements = assign_annotation_lanes(
        selected,
        record_lengths=record_lengths,
        letters_in_lane=letters_in_lane,
        padding_bp=padding,
        overflow="compress" if highlight_track else track_params.overflow,
        join_origin=join_origin,
        label_position=label_position,
    )
    if highlight_track:
        placements = tuple(replace(item, lane=0) for item in placements)
    spans_slots = track_params.to_slot is not None or track_params.cover_anchor
    lane_count = max((item.lane for item in placements), default=-1) + 1
    lane_gap = float(track_params.lane_gap_px)
    bare = 2.0 * float(track_params.padding_px) + lane_count * float(lane_extent_px) + max(0, lane_count - 1) * lane_gap
    lane_rooms, labels_without_room = (
        ((), 0)
        if spans_slots or not track_params.show_labels
        else _rooms_within(
            _lane_rooms(
                placements,
                record_lengths=record_lengths,
                style_of=style_of,
                extent_at=extent_at,
                mark_offset_px=lambda placement, room_before: mark_offset(
                    placement.annotation, placement.lane, max(1, lane_count), room_before
                ),
                lane_pitch_px=lane_pitch(max(1, lane_count)),
                labels_above=labels.labels_above,
                join_origin=join_origin,
                label_position=label_position,
                at_least=room_bounds.at_least,
            ),
            bare_extent_px=bare,
            limit_px=room_bounds.limit_px,
        )
    )
    room = _largest_room(lane_rooms)
    required = bare + room
    clipped = False
    # A to_slot span takes its extent from the two slots, not from the slot's own size.
    if (
        track_params.to_slot is None
        and available_extent_px is not None
        and required > float(available_extent_px)
    ):
        if track_params.overflow == "error":
            raise ValidationError(
                f"Annotation slot {slot_id!r} requires {required:.2f}px but only "
                f"{float(available_extent_px):.2f}px is available."
            )
        if track_params.overflow == "compress" and lane_count > 1:
            lane_gap = max(
                0.0,
                min(
                    lane_gap,
                    (float(available_extent_px) - 2.0 * track_params.padding_px - lane_count * lane_extent_px - room)
                    / (lane_count - 1),
                ),
            )
            required = 2.0 * track_params.padding_px + lane_count * lane_extent_px + (lane_count - 1) * lane_gap + room
            if required > float(available_extent_px):
                raise ValidationError(
                    f"Annotation slot {slot_id!r} cannot fit {lane_count} lanes in {available_extent_px:.2f}px."
                )
        elif track_params.overflow == "clip":
            clipped = True
    occupied = tuple(bound for item in placements for bound in item.occupied_bounds)
    return ResolvedAnnotationTrack(
        slot_id=str(slot_id),
        set_id=str(set_id),
        placements=placements,
        required_extent_px=float(required),
        occupied_bounds=occupied,
        label_font_size=float(labels.font_size),
        join_origin=join_origin,
        clipped=clipped,
        lane_gap_px=lane_gap,
        lane_rooms=lane_rooms,
        labels_above=labels.labels_above,
        labels_without_room=labels_without_room,
    )


__all__ = [
    "LabelExtent",
    "LabelMeasure",
    "LabelSizing",
    "LaneRoom",
    "LettersInLane",
    "PlacedRegionAnnotation",
    "ResolvedAnnotationTrack",
    "RoomBounds",
    "assign_annotation_lanes",
    "layout_annotation_track",
]
