"""Deterministic lane packing shared by circular and linear annotations."""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Callable, Mapping, Sequence

from gbdraw.exceptions import ValidationError

from .models import AnnotationTrackParams, ResolvedRegionAnnotation


@dataclass(frozen=True)
class PlacedRegionAnnotation:
    annotation: ResolvedRegionAnnotation
    lane: int
    mark_geometry: object | None = None
    label_geometry: object | None = None
    occupied_bounds: tuple[tuple[float, float], ...] = ()


@dataclass(frozen=True)
class ResolvedAnnotationTrack:
    slot_id: str
    set_id: str
    placements: tuple[PlacedRegionAnnotation, ...]
    required_extent_px: float
    occupied_bounds: tuple[tuple[float, float], ...]
    clipped: bool = False
    lane_gap_px: float = 3.0


def _expanded_intervals(
    annotation: ResolvedRegionAnnotation,
    *,
    record_length: int,
    label_width_bp: float,
    padding_bp: float,
) -> tuple[tuple[float, float], ...]:
    intervals = [(float(start), float(end)) for start, end in annotation.segments]
    if annotation.label and label_width_bp > annotation.span_bp:
        half = label_width_bp / 2.0
        center = float(annotation.midpoint_bp)
        start, end = center - half, center + half
        if start < 0:
            intervals.extend(((0.0, end), (record_length + start, float(record_length))))
        elif end > record_length:
            intervals.extend(((start, float(record_length)), (0.0, end - record_length)))
        else:
            intervals.append((start, end))
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
) -> tuple[PlacedRegionAnnotation, ...]:
    """Assign stable lanes using greedy interval coloring.

    Explicit lanes remain fixed. A collision in an explicit lane is rejected by
    the default ``error`` policy and retained for ``compress``/``clip`` so the
    caller can apply its visual policy.
    """

    if overflow not in {"error", "compress", "clip"}:
        raise ValidationError("overflow must be 'error', 'compress', or 'clip'.")
    width_func = label_width_bp or (lambda _annotation: 0.0)
    bounds_by_annotation: dict[tuple[str, str, int], tuple[tuple[float, float], ...]] = {}
    for annotation in annotations:
        length = int(record_lengths[annotation.record_index])
        annotation_padding = padding_bp(annotation) if callable(padding_bp) else padding_bp
        bounds_by_annotation[(annotation.set_id, annotation.id, annotation.record_index)] = _expanded_intervals(
            annotation,
            record_length=length,
            label_width_bp=max(0.0, float(width_func(annotation))),
            padding_bp=max(0.0, float(annotation_padding)),
        )

    ordered = sorted(
        annotations,
        key=lambda item: (
            item.record_index,
            min(start for start, _ in item.segments),
            max(end for _, end in item.segments),
            item.set_id,
            item.id,
        ),
    )
    lane_bounds: dict[tuple[int, int], list[tuple[tuple[float, float], ...]]] = {}
    placements: list[PlacedRegionAnnotation] = []
    for annotation in ordered:
        key = (annotation.set_id, annotation.id, annotation.record_index)
        bounds = bounds_by_annotation[key]
        if annotation.lane is not None:
            lane = int(annotation.lane)
            collisions = lane_bounds.get((annotation.record_index, lane), ())
            if overflow == "error" and any(_overlaps(bounds, occupied) for occupied in collisions):
                raise ValidationError(
                    f"Annotation {annotation.set_id}/{annotation.id} conflicts in explicit lane {lane}."
                )
        else:
            lane = 0
            while any(
                _overlaps(bounds, occupied)
                for occupied in lane_bounds.get((annotation.record_index, lane), ())
            ):
                lane += 1
        lane_bounds.setdefault((annotation.record_index, lane), []).append(bounds)
        placements.append(
            PlacedRegionAnnotation(
                annotation=annotation,
                lane=lane,
                occupied_bounds=bounds,
            )
        )
    placements.sort(
        key=lambda item: (
            item.annotation.record_index,
            item.lane,
            min(start for start, _ in item.annotation.segments),
            item.annotation.set_id,
            item.annotation.id,
        )
    )
    return tuple(placements)


def layout_annotation_track(
    slot_id: str,
    set_id: str,
    annotations: Sequence[ResolvedRegionAnnotation],
    *,
    record_lengths: Mapping[int, int],
    params: AnnotationTrackParams | None = None,
    available_extent_px: float | None = None,
    lane_extent_px: float = 14.0,
    bp_per_px: Mapping[int, float] | None = None,
) -> ResolvedAnnotationTrack:
    """Pack one set and calculate its radial/vertical required extent."""

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
    font_default = 10.0

    def label_width(annotation: ResolvedRegionAnnotation) -> float:
        if not track_params.show_labels or not annotation.label:
            return 0.0
        font_size = annotation.style.label_font_size or font_default
        width_px = len(annotation.label) * float(font_size) * 0.6
        scale = float((bp_per_px or {}).get(annotation.record_index, 1.0))
        return width_px * scale

    def padding(annotation: ResolvedRegionAnnotation) -> float:
        scale = float((bp_per_px or {}).get(annotation.record_index, 1.0))
        cap = annotation.style.stroke_width * 2.0 if annotation.mark == "bracket" else 0.0
        return (track_params.padding_px + cap) * scale

    highlight_track = track_params.marks == ("highlight",)
    placements = assign_annotation_lanes(
        selected,
        record_lengths=record_lengths,
        label_width_bp=label_width,
        padding_bp=padding,
        overflow="compress" if highlight_track else track_params.overflow,
    )
    if highlight_track:
        placements = tuple(replace(item, lane=0) for item in placements)
    lane_count = max((item.lane for item in placements), default=-1) + 1
    lane_gap = float(track_params.lane_gap_px)
    required = (
        2.0 * float(track_params.padding_px)
        + lane_count * float(lane_extent_px)
        + max(0, lane_count - 1) * lane_gap
    )
    clipped = False
    if available_extent_px is not None and required > float(available_extent_px):
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
                    (float(available_extent_px) - 2.0 * track_params.padding_px - lane_count * lane_extent_px)
                    / (lane_count - 1),
                ),
            )
            required = 2.0 * track_params.padding_px + lane_count * lane_extent_px + (lane_count - 1) * lane_gap
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
        clipped=clipped,
        lane_gap_px=lane_gap,
    )


__all__ = [
    "PlacedRegionAnnotation",
    "ResolvedAnnotationTrack",
    "assign_annotation_lanes",
    "layout_annotation_track",
]
