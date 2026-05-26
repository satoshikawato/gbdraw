"""SVG path helpers for circular conservation rings."""

from __future__ import annotations

import math


def _point_on_radius(radius: float, position: float, total_length: int) -> tuple[float, float]:
    angle = math.radians(360.0 * (float(position) / float(total_length)) - 90.0)
    return float(radius) * math.cos(angle), float(radius) * math.sin(angle)


def generate_full_annulus_path_desc(
    *,
    inner_radius_px: float,
    outer_radius_px: float,
    total_length: int,
) -> str:
    """Return a closed annulus path using two half arcs per edge."""

    midpoint = float(total_length) / 2.0
    inner_start = _point_on_radius(inner_radius_px, 0.0, total_length)
    inner_mid = _point_on_radius(inner_radius_px, midpoint, total_length)
    outer_start = _point_on_radius(outer_radius_px, 0.0, total_length)
    outer_mid = _point_on_radius(outer_radius_px, midpoint, total_length)
    return (
        f"M {inner_start[0]},{inner_start[1]} "
        f"A{inner_radius_px},{inner_radius_px} 0 1 1 {inner_mid[0]},{inner_mid[1]} "
        f"A{inner_radius_px},{inner_radius_px} 0 1 1 {inner_start[0]},{inner_start[1]} "
        f"L {outer_start[0]},{outer_start[1]} "
        f"A{outer_radius_px},{outer_radius_px} 0 1 0 {outer_mid[0]},{outer_mid[1]} "
        f"A{outer_radius_px},{outer_radius_px} 0 1 0 {outer_start[0]},{outer_start[1]} z"
    )


def generate_annular_hsp_path_desc(
    *,
    draw_start: float,
    draw_end: float,
    total_length: int,
    inner_radius_px: float,
    outer_radius_px: float,
    full_reference: bool = False,
) -> str:
    """Return a closed annular sector for a normalized HSP span."""

    if full_reference:
        return generate_full_annulus_path_desc(
            inner_radius_px=inner_radius_px,
            outer_radius_px=outer_radius_px,
            total_length=total_length,
        )

    span = max(0.0, min(float(total_length), float(draw_end)) - max(0.0, float(draw_start)))
    start = max(0.0, min(float(total_length), float(draw_start)))
    end = max(0.0, min(float(total_length), float(draw_end)))
    angle_deg = 360.0 * span / float(total_length)

    inner_start = _point_on_radius(inner_radius_px, start, total_length)
    inner_end = _point_on_radius(inner_radius_px, end, total_length)
    outer_start = _point_on_radius(outer_radius_px, start, total_length)
    outer_end = _point_on_radius(outer_radius_px, end, total_length)

    if angle_deg > 20.0:
        midpoint = start + (span / 2.0)
        inner_mid = _point_on_radius(inner_radius_px, midpoint, total_length)
        outer_mid = _point_on_radius(outer_radius_px, midpoint, total_length)
        return (
            f"M {inner_start[0]},{inner_start[1]} "
            f"A{inner_radius_px},{inner_radius_px} 0 0 1 {inner_mid[0]},{inner_mid[1]} "
            f"A{inner_radius_px},{inner_radius_px} 0 0 1 {inner_end[0]},{inner_end[1]} "
            f"L {outer_end[0]},{outer_end[1]} "
            f"A{outer_radius_px},{outer_radius_px} 0 0 0 {outer_mid[0]},{outer_mid[1]} "
            f"A{outer_radius_px},{outer_radius_px} 0 0 0 {outer_start[0]},{outer_start[1]} z"
        )

    return (
        f"M {inner_start[0]},{inner_start[1]}"
        f"A{inner_radius_px},{inner_radius_px} 0 0 1 {inner_end[0]},{inner_end[1]}"
        f" L{outer_end[0]},{outer_end[1]}"
        f"A{outer_radius_px},{outer_radius_px} 0 0 0 {outer_start[0]},{outer_start[1]} z"
    )


__all__ = ["generate_annular_hsp_path_desc", "generate_full_annulus_path_desc"]
