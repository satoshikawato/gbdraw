"""Annotation geometry both diagram modes share: DOM ids, the painted edge of a mark, and label ink."""

from __future__ import annotations

import hashlib
import math
import re
from dataclasses import dataclass, replace
from typing import Literal

from gbdraw.annotations import RegionAnnotationStyle, ResolvedRegionAnnotation, is_auto_feature_underlay
from gbdraw.core.text import TEXT_ANCHOR_SHARE, text_ink_offsets
from gbdraw.layout.text_geometry import TextAnchor


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
    record's length.
    """

    start: int
    end: int
    capped_start: bool
    capped_end: bool
    row: int = 0


def mark_parts(annotation: ResolvedRegionAnnotation, *, join_origin: bool) -> tuple[MarkPart, ...]:
    """The parts a mark draws, in the order of its geometry segments: the one owner of where a mark is cut.

    A span across the origin (``crosses_origin``) is cut there: its part that ends at the record's end and its
    part that starts at 0 leave those ends open, or with ``join_origin`` (a circle drawn whole) become one part.
    """

    fragments = annotation.display_parts
    parts = [
        MarkPart(start, end, True, True)
        for start, end in annotation.segments
    ] if fragments is None else [
        MarkPart(part.display_start, part.display_end, not part.artificial_start, not part.artificial_end, part.row)
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
    "MarkPart",
    "MarkReach",
    "annotation_dom_id",
    "cap_size",
    "dom_token",
    "label_ink_offsets",
    "mark_parts",
    "mark_reach",
    "past_span_edge",
]
