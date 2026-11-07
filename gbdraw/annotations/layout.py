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
    """Room between lane ``gap`` and lane ``gap + 1`` of one record for the labels that would cross the next
    lane's mark (PD-OI-119, lane room): their letters, their ``label_offset`` and their mark's stroke."""

    record_index: int
    gap: int
    px: float


@dataclass(frozen=True)
class LabelSizing:
    """How the packer measures a track's labels.

    ``font_size`` is the font size of the drawing's feature labels, which a label whose style sets none takes;
    ``bp_per_px`` the bases one pixel spans on each record (1 when not given). ``depth_px`` gives the pixels a
    label's letters take across the lanes: from its text, its effective style, its font size and the bp it stands
    at. ``labels_above`` says that the labels stand toward lane - 1 (Linear labels above their marks); else they
    stand toward lane + 1 (Linear labels below, and Circular, away from the axis on either side).
    """

    font_size: float
    depth_px: Callable[[str, RegionAnnotationStyle, float, float], float]
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
    # The room between lanes for labels that would cross the next lane's mark (``LaneRoom``).
    lane_rooms: tuple[LaneRoom, ...] = ()

    def room_before(self, record_index: int, lane: int) -> float:
        """The label room between the record's lanes from lane 0 up to ``lane``: how far ``lane`` moves on."""

        return sum(room.px for room in self.lane_rooms if room.record_index == record_index and room.gap < lane)

    def record_room(self, record_index: int) -> float:
        """All the label room between the record's lanes: the slot's lanes share what is left."""

        return _record_room(self.lane_rooms, record_index)


def _record_room(lane_rooms: Sequence[LaneRoom], record_index: int) -> float:
    return sum(room.px for room in lane_rooms if room.record_index == record_index)


class _LabelLetters(NamedTuple):
    """Where a label stands, the bp its letters span on its row, and whether they leave the stretch of their mark
    they may cover (``ResolvedRegionAnnotation.label_stretch``)."""

    place: LabelPlace
    intervals: tuple[tuple[float, float], ...]
    leaves_mark: bool


def _label_letters(
    annotation: ResolvedRegionAnnotation,
    *,
    record_length: int,
    label_width_bp: float,
    join_origin: bool,
    label_position: LabelPosition,
) -> _LabelLetters:
    """Where an annotation's label stands and what its letters cover (``join_origin`` as in ``mark_parts``).

    A start label runs right from its place, an end label left, a centred one both ways. On a circle the letters
    go on across the origin; a Linear row stops them at its ends.
    """

    place, (low, high) = annotation.label_stretch(label_position, join_origin=join_origin)
    lead = {"start": 0.0, "end": 1.0}.get(label_position, 0.5) * label_width_bp
    start, end = place.bp - lead, place.bp - lead + label_width_bp
    leaves = start < low or end > high
    if not join_origin:
        return _LabelLetters(place, ((max(0.0, start), min(float(record_length), end)),), leaves)
    # A part joined across the origin runs on past the record's end.
    turns = math.floor(start / record_length)
    start, end = start - turns * record_length, end - turns * record_length
    if end > record_length:
        return _LabelLetters(place, ((start, float(record_length)), (0.0, end - record_length)), leaves)
    return _LabelLetters(place, ((start, end),), leaves)


def _expanded_intervals(
    annotation: ResolvedRegionAnnotation,
    row: int,
    *,
    record_length: int,
    label_width_bp: float,
    padding_bp: float,
    join_origin: bool,
    label_position: LabelPosition,
) -> tuple[tuple[float, float], ...]:
    """What an annotation occupies on one row, with its label's letters where they leave its mark
    (``_label_letters``)."""

    intervals = [(float(start), float(end)) for start, end in annotation.row_segments(row)]
    if annotation.label and label_width_bp > 0.0:
        letters = _label_letters(
            annotation,
            record_length=record_length,
            label_width_bp=label_width_bp,
            join_origin=join_origin,
            label_position=label_position,
        )
        if letters.leaves_mark and row == letters.place.part.row:
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
    label_width_bp: Callable[[ResolvedRegionAnnotation], float] | None = None,
    padding_bp: Callable[[ResolvedRegionAnnotation], float] | float = 0.0,
    overflow: str = "error",
    join_origin: bool,
    label_position: Callable[[ResolvedRegionAnnotation], LabelPosition] | None = None,
) -> tuple[PlacedRegionAnnotation, ...]:
    """Assign stable lanes using greedy interval coloring, row by row on a wrapped record.

    ``join_origin`` says whether the drawing closes into a circle (``mark_parts``) and ``label_position`` where
    each label stands on it (centred when not given): the lane keeps room for the label there.

    An annotation keeps one lane on every row it draws on; only annotations on the same row meet.
    Explicit lanes remain fixed. A collision in an explicit lane is rejected by
    the default ``error`` policy and retained for ``compress``/``clip`` so the
    caller can apply its visual policy.
    """

    if overflow not in {"error", "compress", "clip"}:
        raise ValidationError("overflow must be 'error', 'compress', or 'clip'.")
    width_func = label_width_bp or (lambda _annotation: 0.0)
    bounds_by_annotation: dict[tuple[str, str, int], dict[int, tuple[tuple[float, float], ...]]] = {}
    for annotation in annotations:
        length = int(record_lengths[annotation.record_index])
        annotation_padding = padding_bp(annotation) if callable(padding_bp) else padding_bp
        bounds_by_annotation[(annotation.set_id, annotation.id, annotation.record_index)] = {
            row: _expanded_intervals(
                annotation,
                row,
                record_length=length,
                label_width_bp=max(0.0, float(width_func(annotation))),
                padding_bp=max(0.0, float(annotation_padding)),
                join_origin=join_origin,
                label_position="center" if label_position is None else label_position(annotation),
            )
            for row in annotation.geometry_rows
        }

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
        key = (annotation.set_id, annotation.id, annotation.record_index)
        rows = bounds_by_annotation[key]

        def lane_is_taken(lane: int) -> bool:
            return any(
                _overlaps(bounds, occupied)
                for row, bounds in rows.items()
                for occupied in lane_bounds.get((annotation.record_index, row, lane), ())
            )

        if annotation.lane is not None:
            lane = int(annotation.lane)
            if overflow == "error" and lane_is_taken(lane):
                raise ValidationError(
                    f"Annotation {annotation.set_id}/{annotation.id} conflicts in explicit lane {lane}.",
                    diagnostic={"code": "TRACK_INVALID", "reason": "ANNOTATION_LANE_CONFLICT"},
                )
        else:
            lane = 0
            while lane_is_taken(lane):
                lane += 1
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


def _lane_rooms(
    placements: Sequence[PlacedRegionAnnotation],
    *,
    record_lengths: Mapping[int, int],
    label_width_bp: Callable[[ResolvedRegionAnnotation], float],
    label_room_px: Callable[[ResolvedRegionAnnotation, float], float],
    labels_above: bool,
    join_origin: bool,
    label_position: Callable[[ResolvedRegionAnnotation], LabelPosition],
) -> tuple[LaneRoom, ...]:
    """The room between lanes for the labels that would cross a mark of the next lane on their side (PD-OI-119, lane room).

    A label would cross where its letters' bp at its place (``_label_letters``) meet a mark of the next lane on
    the label's row. The room between two lanes is the largest ``label_room_px`` (at the label's place) of those
    labels.
    """

    marks: dict[tuple[int, int, int], list[tuple[float, float]]] = {}
    for placement in placements:
        annotation = placement.annotation
        for row in annotation.geometry_rows:
            marks.setdefault((annotation.record_index, row, placement.lane), []).extend(
                (float(start), float(end)) for start, end in annotation.row_segments(row)
            )
    rooms: dict[tuple[int, int], float] = {}
    for placement in placements:
        annotation = placement.annotation
        width = max(0.0, float(label_width_bp(annotation)))
        next_lane = placement.lane - 1 if labels_above else placement.lane + 1
        if width <= 0.0 or next_lane < 0:
            continue
        letters = _label_letters(
            annotation,
            record_length=int(record_lengths[annotation.record_index]),
            label_width_bp=width,
            join_origin=join_origin,
            label_position=label_position(annotation),
        )
        if _overlaps(letters.intervals, marks.get((annotation.record_index, letters.place.part.row, next_lane), ())):
            key = (annotation.record_index, min(placement.lane, next_lane))
            rooms[key] = max(rooms.get(key, 0.0), float(label_room_px(annotation, letters.place.bp)))
    return tuple(LaneRoom(record, gap, px) for (record, gap), px in sorted(rooms.items()))


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
) -> ResolvedAnnotationTrack:
    """Pack one set and calculate its radial/vertical required extent.

    ``labels`` measures the labels (``LabelSizing``). ``join_origin`` is True for a Circular drawing
    (``mark_parts``). Two lanes part by the room of a label that would cross the next lane's mark
    (``_lane_rooms``), and the track's required extent grows by the largest room of a record; a span over other
    slots (``to_slot``, ``cover_anchor``) puts its labels outside them and gets none."""

    track_params = params or AnnotationTrackParams(set_id=set_id)
    if track_params.set_id != set_id:
        raise ValidationError(
            f"Annotation slot {slot_id!r} requests set {track_params.set_id!r}, not {set_id!r}."
        )
    selected = tuple(
        item
        for item in annotations
        if item.set_id == set_id
        and (track_params.marks is None or item.mark in track_params.marks)
    )

    def label_width(annotation: ResolvedRegionAnnotation) -> float:
        if not track_params.show_labels or not annotation.label:
            return 0.0
        style = effective_annotation_style(annotation, track_params)
        width_px = len(annotation.label) * annotation_label_font_size(style, labels.font_size) * 0.6
        scale = float(labels.bp_per_px.get(annotation.record_index, 1.0))
        return width_px * scale

    def label_room(annotation: ResolvedRegionAnnotation, place_bp: float) -> float:
        # The letters keep label_offset from their mark's paint, half its stroke past the lane; the next mark's
        # paint reaches about as far into the room.
        style = effective_annotation_style(annotation, track_params)
        depth = labels.depth_px(annotation.label, style, annotation_label_font_size(style, labels.font_size), place_bp)
        return float(depth) + float(style.label_offset) + float(style.stroke_width)

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
        label_width_bp=label_width,
        padding_bp=padding,
        overflow="compress" if highlight_track else track_params.overflow,
        join_origin=join_origin,
        label_position=label_position,
    )
    if highlight_track:
        placements = tuple(replace(item, lane=0) for item in placements)
    spans_slots = track_params.to_slot is not None or track_params.cover_anchor
    lane_rooms = (
        ()
        if spans_slots
        else _lane_rooms(
            placements,
            record_lengths=record_lengths,
            label_width_bp=label_width,
            label_room_px=label_room,
            labels_above=labels.labels_above,
            join_origin=join_origin,
            label_position=label_position,
        )
    )
    room = max((_record_room(lane_rooms, index) for index in record_lengths), default=0.0)
    lane_count = max((item.lane for item in placements), default=-1) + 1
    lane_gap = float(track_params.lane_gap_px)
    required = (
        2.0 * float(track_params.padding_px)
        + lane_count * float(lane_extent_px)
        + max(0, lane_count - 1) * lane_gap
        + room
    )
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
    )


__all__ = [
    "LabelSizing",
    "LaneRoom",
    "PlacedRegionAnnotation",
    "ResolvedAnnotationTrack",
    "assign_annotation_lanes",
    "layout_annotation_track",
]
