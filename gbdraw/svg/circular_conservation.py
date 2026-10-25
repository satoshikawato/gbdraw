"""SVG path helpers for circular conservation rings."""

from __future__ import annotations

from .circular_arcs import circle_arc_to, circle_point, position_angle
from .path_data import path_point


def _annular_sector(start: float, span: float, total_length: int, inner_radius_px: float, outer_radius_px: float) -> str:
    start_angle = position_angle(start, total_length)
    span_deg = 360.0 * span / float(total_length)
    end_angle = start_angle + span_deg
    return (
        f"M {path_point(*circle_point(inner_radius_px, start_angle))} "
        f"{circle_arc_to(inner_radius_px, start_angle, span_deg)} "
        f"L {path_point(*circle_point(outer_radius_px, end_angle))} "
        f"{circle_arc_to(outer_radius_px, end_angle, -span_deg)} z"
    )


def generate_full_annulus_path_desc(
    *,
    inner_radius_px: float,
    outer_radius_px: float,
    total_length: int,
) -> str:
    """Return a closed annulus path: the inner circle, then the outer circle back."""

    return _annular_sector(0.0, float(total_length), total_length, inner_radius_px, outer_radius_px)


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
    start = max(0.0, min(float(total_length), float(draw_start)))
    span = max(0.0, min(float(total_length), float(draw_end)) - max(0.0, float(draw_start)))
    return _annular_sector(start, span, total_length, inner_radius_px, outer_radius_px)


__all__ = ["generate_annular_hsp_path_desc", "generate_full_annulus_path_desc"]
