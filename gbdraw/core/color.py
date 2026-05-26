#!/usr/bin/env python
# coding: utf-8

import re


_HEX_COLOR_PATTERN = re.compile(r"^#(?P<value>[0-9A-Fa-f]{3}|[0-9A-Fa-f]{6})$")

COLLINEAR_ORIENTATION_COLOR_KEYS = {
    "plus": "collinear_block_plus",
    "minus": "collinear_block_minus",
}
COLLINEAR_ORIENTATION_MIN_COLOR_KEYS = {
    "plus": "collinear_block_plus_min",
    "minus": "collinear_block_minus_min",
}
DEFAULT_COLLINEAR_ORIENTATION_MIN_COLORS = {
    "plus": "#f0f1f5",
    "minus": "#FFE7E7",
}
DEFAULT_COLLINEAR_ORIENTATION_COLORS = {
    "plus": "#8b9cc1",
    "minus": "#E15759",
}


def normalize_hex_color(color: str) -> str:
    text = str(color or "").strip()
    match = _HEX_COLOR_PATTERN.match(text)
    if not match:
        raise ValueError(f"Expected #RGB or #RRGGBB color, got {color!r}")
    value = match.group("value")
    if len(value) == 3:
        value = "".join(char * 2 for char in value)
    return f"#{value.lower()}"


def _hex_to_rgb(color: str) -> tuple[int, int, int]:
    normalized = normalize_hex_color(color)
    return (
        int(normalized[1:3], 16),
        int(normalized[3:5], 16),
        int(normalized[5:7], 16),
    )


def blend_color(color: str, target_color: str, factor: float) -> str:
    """Blend color toward target_color by factor in [0, 1]."""
    r_min, g_min, b_min = _hex_to_rgb(color)
    r_max, g_max, b_max = _hex_to_rgb(target_color)
    clamped_factor = max(0.0, min(1.0, float(factor)))

    r = int(r_min + (r_max - r_min) * clamped_factor)
    g = int(g_min + (g_max - g_min) * clamped_factor)
    b = int(b_min + (b_max - b_min) * clamped_factor)

    return f"#{r:02x}{g:02x}{b:02x}"


def tint_color(color: str, factor: float = 0.82) -> str:
    return blend_color(color, "#ffffff", factor)


def interpolate_color(color_min: str, color_max: str, factor: float) -> str:
    r_min, g_min, b_min = _hex_to_rgb(color_min)
    r_max, g_max, b_max = _hex_to_rgb(color_max)

    r = int(r_min + (r_max - r_min) * factor)
    g = int(g_min + (g_max - g_min) * factor)
    b = int(b_min + (b_max - b_min) * factor)

    return f"#{r:02x}{g:02x}{b:02x}"


__all__ = [
    "COLLINEAR_ORIENTATION_COLOR_KEYS",
    "COLLINEAR_ORIENTATION_MIN_COLOR_KEYS",
    "DEFAULT_COLLINEAR_ORIENTATION_MIN_COLORS",
    "DEFAULT_COLLINEAR_ORIENTATION_COLORS",
    "blend_color",
    "interpolate_color",
    "normalize_hex_color",
    "tint_color",
]


