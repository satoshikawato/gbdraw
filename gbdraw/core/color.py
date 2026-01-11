#!/usr/bin/env python
# coding: utf-8

def interpolate_color(color_min: str, color_max: str, factor: float) -> str:
    r_min, g_min, b_min = (
        int(color_min[1:3], 16),
        int(color_min[3:5], 16),
        int(color_min[5:7], 16),
    )
    r_max, g_max, b_max = (
        int(color_max[1:3], 16),
        int(color_max[3:5], 16),
        int(color_max[5:7], 16),
    )

    r = int(r_min + (r_max - r_min) * factor)
    g = int(g_min + (g_max - g_min) * factor)
    b = int(b_min + (b_max - b_min) * factor)

    return f"#{r:02x}{g:02x}{b:02x}"


__all__ = ["interpolate_color"]


