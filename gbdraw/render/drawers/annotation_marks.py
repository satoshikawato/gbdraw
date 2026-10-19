"""Annotation geometry both diagram modes share: DOM ids, the painted edge of a mark, and label ink."""

from __future__ import annotations

import hashlib
import math
import re
from dataclasses import dataclass

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


def capped_ends(annotation: ResolvedRegionAnnotation, part_index: int) -> tuple[bool, bool]:
    """Whether a mark part draws its start and end: a display fragment's cut ends are open."""

    fragment = annotation.display_parts[part_index] if annotation.display_parts is not None else None
    return (
        fragment is None or not fragment.artificial_start,
        fragment is None or not fragment.artificial_end,
    )


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
        ends = [capped_ends(annotation, index) for index in range(len(annotation.geometry_segments))]
        if any(start for start, _end in ends):
            reach = max(reach, cap)
        if any(end for _start, end in ends):
            # An arrow's butt-capped diagonals reach half the stroke times cos 45 degrees past the cap.
            head = half_stroke * math.sqrt(0.5) if arrow_heads and style.line_cap == "arrow" else 0.0
            reach = max(reach, cap + head)
    return MarkReach(reach, reach - half_lane)


def past_span_edge(reach: MarkReach, inset: float) -> float:
    """How far a span's mark paints past the span's edge on the label side, when its lane lies ``inset`` inside
    that edge; zero while the paint stays within the span. A span's label keeps ``label_offset`` beyond that."""

    return max(0.0, reach.past_lane_edge - inset)


def label_ink_offsets(
    text: str, font_family: str, font_size: float, anchor: TextAnchor, *, dominant_baseline: str = "central"
) -> tuple[float, float, float, float]:
    """``text_ink_offsets`` of a label (central baseline unless named), or an em box about its point
    when the font cannot be read."""

    offsets = text_ink_offsets(text, font_family, font_size, dominant_baseline=dominant_baseline, text_anchor=anchor)
    if offsets is None:
        # The fallback width of calculate_bbox_dimensions.
        width = 0.6 * len(text) * font_size
        origin = -TEXT_ANCHOR_SHARE[anchor] * width
        return origin, origin + width, -0.5 * font_size, 0.5 * font_size
    return offsets


__all__ = [
    "MarkReach",
    "annotation_dom_id",
    "cap_size",
    "capped_ends",
    "dom_token",
    "label_ink_offsets",
    "mark_reach",
    "past_span_edge",
]
