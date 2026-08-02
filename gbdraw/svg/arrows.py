#!/usr/bin/env python
# coding: utf-8

import math
from typing import Literal


ArrowHeadLengthRatio = Literal["auto"] | float


def set_arrow_shoulder(feat_strand: str, arrow_end: float, cds_arrow_length: float) -> float:
    """
    Calculates the shoulder position for an arrow representing a genomic feature.
    """
    if feat_strand == "positive":
        shoulder = float(arrow_end - cds_arrow_length)
    else:
        shoulder = float(arrow_end + cds_arrow_length)
    return shoulder


def calculate_circular_arrow_length(total_length: int) -> float:
    """Return circular arrow head length (bp) derived from genome size."""
    min_arrow_length = 30.0
    max_arrow_length = 700.0
    param_a = 3.0
    param_b = 5.0
    return min_arrow_length + (max_arrow_length - min_arrow_length) * (
        1.0 / (1.0 + math.exp(-param_a * (math.log10(total_length) - param_b)))
    )


def resolve_arrow_head_length_px(
    head_length_ratio: ArrowHeadLengthRatio,
    feature_thickness_px: float,
    automatic_head_length_px: float,
    shaft_width_ratio: float = 1.0,
) -> float:
    """Resolve a head in display space, extending Auto for a narrowed shaft."""
    if head_length_ratio == "auto":
        return float(automatic_head_length_px) + abs(
            float(feature_thickness_px)
        ) * (1.0 - float(shaft_width_ratio))
    return abs(float(feature_thickness_px)) * float(head_length_ratio)


def circular_head_length_bp_to_px(
    head_length_bp: float,
    total_length: int,
    center_radius_px: float,
) -> float:
    """Convert a Circular head span from base pairs to arc length in pixels."""
    return (
        (2.0 * math.pi * abs(float(center_radius_px)))
        * (float(head_length_bp) / float(total_length))
    )


def circular_head_length_px_to_bp(
    head_length_px: float,
    total_length: int,
    center_radius_px: float,
) -> float:
    """Convert a Circular display-space head length to a floating-point bp span."""
    circumference_px = 2.0 * math.pi * abs(float(center_radius_px))
    if circumference_px == 0.0:
        raise ValueError("center_radius_px must be nonzero")
    return float(head_length_px) * float(total_length) / circumference_px


def resolve_circular_arrow_head_length_bp(
    head_length_ratio: ArrowHeadLengthRatio,
    feature_thickness_px: float,
    automatic_head_length_bp: float,
    total_length: int,
    center_radius_px: float,
    shaft_width_ratio: float = 1.0,
) -> float:
    """Resolve a Circular head span, keeping full-width Auto exactly legacy."""
    if head_length_ratio == "auto" and float(shaft_width_ratio) == 1.0:
        return float(automatic_head_length_bp)
    automatic_head_length_px = (
        circular_head_length_bp_to_px(
            automatic_head_length_bp,
            total_length,
            center_radius_px,
        )
        if head_length_ratio == "auto"
        else 0.0
    )
    head_length_px = resolve_arrow_head_length_px(
        head_length_ratio,
        feature_thickness_px,
        automatic_head_length_px,
        shaft_width_ratio,
    )
    return circular_head_length_px_to_bp(
        head_length_px,
        total_length,
        center_radius_px,
    )


def cap_arrow_head_length(feature_length: float, head_length: float) -> float:
    """Cap a nonnegative head length at the available terminal block length."""
    return min(abs(float(feature_length)), max(0.0, float(head_length)))


def is_legacy_arrow_short(feature_length: float, head_length: float) -> bool:
    """Return the legacy five-vertex arrow's strict short-feature decision."""
    return abs(float(feature_length)) < float(head_length)


def has_arrow_shaft(feature_length: float, head_length: float) -> bool:
    """Return whether an arrow has positive shaft length."""
    return abs(float(feature_length)) > float(head_length)


def calculate_arrow_shaft_bounds(
    full_inner_or_top: float,
    center: float,
    full_outer_or_bottom: float,
    shaft_width_ratio: float,
) -> tuple[float, float]:
    """Interpolate reduced shaft edges from the center toward full-width edges."""
    ratio = float(shaft_width_ratio)
    center_value = float(center)
    return (
        center_value + (float(full_inner_or_top) - center_value) * ratio,
        center_value + (float(full_outer_or_bottom) - center_value) * ratio,
    )


__all__ = [
    "ArrowHeadLengthRatio",
    "calculate_arrow_shaft_bounds",
    "calculate_circular_arrow_length",
    "cap_arrow_head_length",
    "circular_head_length_bp_to_px",
    "circular_head_length_px_to_bp",
    "has_arrow_shaft",
    "is_legacy_arrow_short",
    "resolve_arrow_head_length_px",
    "resolve_circular_arrow_head_length_bp",
    "set_arrow_shoulder",
]


