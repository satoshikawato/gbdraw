"""Semantic candidate extraction for circular feature labels."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from ..features.coordinates import get_strand
from .circular_types import CircularLabelCandidate
from .filtering import get_label_text


def _candidate_segment(feature: Any, total_length: int) -> tuple[float, float, float, float, str]:
    block_coordinates = [
        coordinate
        for coordinate in getattr(feature, "location", ())
        if getattr(coordinate, "kind", None) == "block"
    ]
    if len(block_coordinates) == 2:
        starts = [int(coordinate.start) for coordinate in block_coordinates]
        ends = [int(coordinate.end) for coordinate in block_coordinates]
        left = [coordinate for coordinate in block_coordinates if int(coordinate.start) <= 1]
        right = [coordinate for coordinate in block_coordinates if int(coordinate.end) >= total_length]
        if len(left) == 1 and len(right) == 1 and left[0] is not right[0]:
            start = float(max(starts))
            end = float(min(ends))
            span = (float(total_length) - start) + end
            middle = (start + (0.5 * span)) % float(total_length)
            return start, end, middle, span, get_strand(block_coordinates[0].strand)

    best: tuple[float, float, float, float, str] | None = None
    for coordinate in getattr(feature, "coordinates", ()):
        start = float(coordinate.start)
        end = float(coordinate.end)
        span = abs(end - start) + 1.0
        candidate = (start, end, 0.5 * (start + end), span, get_strand(coordinate.strand))
        if best is None or span > best[3]:
            best = candidate
    return best or (0.0, 0.0, 0.0, 0.0, "undefined")


def build_circular_label_candidates(
    feature_dict: Mapping[str, Any],
    total_length: int,
    label_filtering: Mapping[str, Any],
    *,
    font_family: str,
    font_size: float,
    measure_text,
) -> tuple[CircularLabelCandidate, ...]:
    """Select and measure each circular feature label exactly once."""
    candidates: list[CircularLabelCandidate] = []
    for input_order, (stable_id, feature) in enumerate(feature_dict.items()):
        text = get_label_text(feature, label_filtering)
        if not text:
            continue
        width_px, height_px = measure_text(text, font_family, font_size)
        start, end, middle, span, strand = _candidate_segment(feature, total_length)
        coordinates = tuple(
            (float(coordinate.start), float(coordinate.end))
            for coordinate in getattr(feature, "coordinates", ())
        )
        candidates.append(
            CircularLabelCandidate(
                stable_id=str(stable_id),
                input_order=input_order,
                text=str(text),
                font_family=str(font_family),
                font_size=float(font_size),
                width_px=float(width_px),
                height_px=float(height_px),
                feature_type=str(getattr(feature, "feature_type", getattr(feature, "type", ""))),
                strand=strand,
                track_id=int(getattr(feature, "feature_track_id", 0)),
                directional=bool(getattr(feature, "is_directional", False)),
                segment_start_bp=start,
                segment_end_bp=end,
                segment_middle_bp=middle,
                segment_span_bp=span,
                feature_coordinates=coordinates,
            )
        )
    return tuple(candidates)


__all__ = ["build_circular_label_candidates"]
