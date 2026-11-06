"""Annotation geometry both diagram modes share: DOM ids, the painted edge of a mark, and label ink."""

from __future__ import annotations

import hashlib
import math
import re
from dataclasses import dataclass, replace
from typing import Literal

from collections.abc import Callable

from svgwrite.base import BaseElement

from gbdraw.annotations import (
    AnnotationGradient,
    RegionAnnotationStyle,
    ResolvedRegionAnnotation,
    is_auto_feature_underlay,
)
from gbdraw.core.color import interpolate_color
from gbdraw.core.text import TEXT_ANCHOR_SHARE, text_ink_offsets
from gbdraw.layout.text_geometry import TextAnchor
from gbdraw.svg.elements import Group


def dom_token(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9_.-]+", "-", str(value)).strip("-") or "annotation"
    return cleaned[:48]


def annotation_dom_id(
    *,
    record_index: int,
    slot_id: str,
    set_id: str,
    annotation_id: str,
) -> str:
    raw = f"{record_index}|{slot_id}|{set_id}|{annotation_id}"
    digest = hashlib.sha1(raw.encode("utf-8")).hexdigest()[:10]
    return f"gbdraw-annotation-{dom_token(slot_id)}-{dom_token(annotation_id)}-{digest}"


def cap_size(lane_width: float) -> float:
    """Half the length of a bracket's end caps across its lane."""

    return min(5.0, 0.35 * lane_width)


@dataclass(frozen=True)
class MarkPart:
    """One drawn part of a mark, in display bp on its ``row``, and whether it draws each of its ends.

    An end it does not draw is a cut, where the region goes on: the display cut, a row boundary, or the origin of
    a span across it. A cut gets no edge stroke, cap or arrow head. A part joined across the origin ends past the
    record's length. ``local_start`` is where the part starts on its record (record-local bp); the part runs on
    from there as it runs on the display, past the record's end when it is joined across the origin.
    """

    start: int
    end: int
    capped_start: bool
    capped_end: bool
    local_start: int
    row: int = 0


def mark_parts(annotation: ResolvedRegionAnnotation, *, join_origin: bool) -> tuple[MarkPart, ...]:
    """The parts a mark draws, in the order of its geometry segments: the one owner of where a mark is cut.

    A span across the origin (``crosses_origin``) is cut there: its part that ends at the record's end and its
    part that starts at 0 leave those ends open, or with ``join_origin`` (a circle drawn whole) become one part.
    """

    fragments = annotation.display_parts
    parts = [
        MarkPart(start, end, True, True, start)
        for start, end in annotation.segments
    ] if fragments is None else [
        MarkPart(
            part.display_start, part.display_end, not part.artificial_start, not part.artificial_end, part.local_start,
            part.row,
        )
        for part in fragments
    ]
    if not annotation.crosses_origin:
        return tuple(parts)
    local_bounds = (
        annotation.segments if fragments is None else [(part.local_start, part.local_end) for part in fragments]
    )
    # The resolver sets crosses_origin only when a segment starts at 0 and the last one ends at the record's end.
    record_end = max(end for _start, end in local_bounds)
    before = next(index for index, (_start, end) in enumerate(local_bounds) if end == record_end)
    after = next(index for index, (start, _end) in enumerate(local_bounds) if start == 0)
    if join_origin:
        joined = parts[after]
        parts[before] = replace(
            parts[before],
            end=parts[before].end + joined.end - joined.start,
            capped_end=joined.capped_end,
        )
        del parts[after]
        return tuple(parts)
    parts[before] = replace(parts[before], capped_end=False)
    parts[after] = replace(parts[after], capped_start=False)
    return tuple(parts)


def _half_stroke(annotation: ResolvedRegionAnnotation, style: RegionAnnotationStyle) -> float:
    return 0.0 if is_auto_feature_underlay(annotation) else 0.5 * float(style.stroke_width)


@dataclass(frozen=True)
class MarkReach:
    """How far a mark's paint reaches across its lane on the label side, stroke included: from the
    lane centre, and past the lane's edge (negative while the paint stays inside the lane)."""

    from_centre: float
    past_lane_edge: float


def mark_reach(
    annotation: ResolvedRegionAnnotation,
    style: RegionAnnotationStyle,
    lane_width: float,
    *,
    fill_share: float,
    arrow_heads: bool,
) -> MarkReach:
    """The painted edge of a mark on its label side (``MarkReach``).

    ``fill_share`` is the part of its lane a ``highlight`` or ``band`` fills;
    ``arrow_heads`` says whether an ``arrow`` cap is drawn as a head of two
    diagonals rather than as a straight cap.
    """

    half_stroke = _half_stroke(annotation, style)
    half_lane = 0.5 * lane_width
    if annotation.mark in {"band", "highlight"}:
        half_fill = 0.5 * fill_share * lane_width
        return MarkReach(half_fill + half_stroke, half_stroke - (half_lane - half_fill))
    reach = half_stroke
    if annotation.mark == "bracket" and style.line_cap != "none":
        cap = cap_size(lane_width)
        parts = mark_parts(annotation, join_origin=False)
        if any(part.capped_start for part in parts):
            reach = max(reach, cap)
        if any(part.capped_end for part in parts):
            # An arrow's butt-capped diagonals reach half the stroke times cos 45 degrees past the cap.
            head = half_stroke * math.sqrt(0.5) if arrow_heads and style.line_cap == "arrow" else 0.0
            reach = max(reach, cap + head)
    return MarkReach(reach, reach - half_lane)


def past_span_edge(reach: MarkReach, inset: float) -> float:
    """How far a span's mark paints past the span's edge on the label side, when its lane lies ``inset`` inside
    that edge; zero while the paint stays within the span. A span's label keeps ``label_offset`` beyond that."""

    return max(0.0, reach.past_lane_edge - inset)


# A gradient fill is drawn as solid steps of at most this many pixels along the mark, and at most this many steps.
GRADIENT_STEP_PX = 2.0
GRADIENT_MAX_STEPS = 48


def _gradient_share(annotation: ResolvedRegionAnnotation, record_length: int) -> Callable[[float], float]:
    """Where a record-local position lies along the region, from 0 at its start to 1 at its end.

    A region across a circular record's origin (``crosses_origin``) runs from the start of its segment that ends
    at the record's end over the origin to the end of its segment that starts at 0.
    """

    segments = annotation.segments
    length = max(1, int(record_length))
    crosses_origin = annotation.crosses_origin
    if crosses_origin:
        low = float(next(start for start, end in segments if end == length))
        high = float(next(end for start, end in segments if start == 0) + length)
    else:
        low = float(min(start for start, _end in segments))
        high = float(max(end for _start, end in segments))
    span = max(1e-9, high - low)

    def share(position: float) -> float:
        if crosses_origin and position < low:
            position += length
        return min(1.0, max(0.0, (position - low) / span))

    return share


def gradient_steps(
    annotation: ResolvedRegionAnnotation,
    gradient: AnnotationGradient,
    part: MarkPart,
    start: float,
    end: float,
    *,
    record_length: int,
    length_px: float,
) -> tuple[tuple[float, str], ...]:
    """The solid steps of a gradient fill over one drawn part ``start`` to ``end`` (its display frame).

    Each step is its start and its colour, the colour at its middle; a step is
    drawn from its start to ``end``, over the steps before it, so that adjacent
    steps leave no seam. The colour follows the region on its record: the
    scale at ``query_start`` where what it stands for starts, by the
    gradient's ``direction`` on the source strand and the record's
    ``source_step``, to ``query_end`` where it ends, the same over every part
    and row of the region.
    """

    share = _gradient_share(annotation, record_length)
    forward = (gradient.direction == "forward") == (annotation.source_step == 1)
    # The region covers query_start..query_end of the scale, from its query start.
    scale_start, span = gradient.query_start, gradient.query_end - gradient.query_start
    count = max(1, min(GRADIENT_MAX_STEPS, math.ceil(max(0.0, float(length_px)) / GRADIENT_STEP_PX)))
    width = (float(end) - float(start)) / count
    steps = []
    for index in range(count):
        middle = float(start) + (index + 0.5) * width
        position = share(part.local_start + (middle - part.start))
        steps.append((
            float(start) + index * width,
            interpolate_color(gradient.start_color, gradient.end_color, scale_start + span * (
                position if forward else 1.0 - position
            )),
        ))
    return tuple(steps)


def gradient_fill_group(
    steps: tuple[tuple[float, str], ...],
    opacity: float,
    shape: Callable[[float, str], BaseElement],
) -> Group:
    """The steps of a gradient fill, each drawn by ``shape(step_start, color)``, in one group at the mark's
    ``fill_opacity``: the group is composited once, so the stacked steps show no overlap."""

    group = Group(opacity=opacity)
    for step_start, color in steps:
        group.add(shape(step_start, color))
    return group


# The em box a label falls back to when its font cannot be read, above and below each baseline it is drawn on (em).
_FALLBACK_EM_BOX = {"central": (-0.5, 0.5), "alphabetic": (-0.8, 0.2)}


def label_ink_offsets(
    text: str,
    font_family: str,
    font_size: float,
    anchor: TextAnchor,
    baseline: Literal["central", "alphabetic"] = "central",
) -> tuple[float, float, float, float]:
    """``text_ink_offsets`` of a label on its ``baseline`` (central for straight text, alphabetic for text along a
    path), or an em box when the font cannot be read."""

    offsets = text_ink_offsets(text, font_family, font_size, dominant_baseline=baseline, text_anchor=anchor)
    if offsets is None:
        # The no-font estimate: 0.6 em per character.
        width = 0.6 * len(text) * font_size
        origin = -TEXT_ANCHOR_SHARE[anchor] * width
        top, bottom = _FALLBACK_EM_BOX[baseline]
        return origin, origin + width, top * font_size, bottom * font_size
    return offsets


__all__ = [
    "GRADIENT_MAX_STEPS",
    "GRADIENT_STEP_PX",
    "MarkPart",
    "MarkReach",
    "annotation_dom_id",
    "cap_size",
    "dom_token",
    "gradient_fill_group",
    "gradient_steps",
    "label_ink_offsets",
    "mark_parts",
    "mark_reach",
    "past_span_edge",
]
