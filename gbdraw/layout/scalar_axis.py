"""Shared scalar-track axis helpers."""

from __future__ import annotations

import math
from collections.abc import Callable
from typing import Any


LINEAR_SCALAR_AXIS_DEFAULT_FONT_MIN_PX = 5.0
LINEAR_SCALAR_AXIS_DEFAULT_FONT_MAX_PX = 8.0
LINEAR_SCALAR_AXIS_DEFAULT_FONT_HEIGHT_FRACTION = 0.7


def has_tick_value(values: list[float], candidate: float) -> bool:
    return any(math.isclose(candidate, value, rel_tol=1e-9, abs_tol=1e-9) for value in values)


def scalar_axis_tick_values(
    axis_min: float,
    axis_max: float,
    *,
    show_ticks: bool,
    large_tick_interval: float | None,
) -> list[float]:
    if not show_ticks:
        return []
    if axis_max <= axis_min:
        return [axis_min]
    if large_tick_interval is None:
        return [axis_min, axis_max]

    interval = float(large_tick_interval)
    if interval <= 0:
        return [axis_min, axis_max]

    tick_values: list[float] = []
    start = math.ceil((axis_min - 1e-9) / interval) * interval
    if start > axis_min + 1e-9:
        tick_values.append(axis_min)
    tick = start
    while tick <= axis_max + 1e-9 and len(tick_values) < 100:
        if tick >= axis_min - 1e-9:
            tick_values.append(0.0 if abs(tick) < 1e-9 else float(tick))
        tick += interval
    if not tick_values or not math.isclose(tick_values[-1], axis_max, rel_tol=1e-9, abs_tol=1e-9):
        tick_values.append(axis_max)

    deduplicated: list[float] = []
    for value in tick_values:
        if not deduplicated or not math.isclose(value, deduplicated[-1], rel_tol=1e-9, abs_tol=1e-9):
            deduplicated.append(value)
    return deduplicated


def scalar_axis_small_tick_values(
    axis_min: float,
    axis_max: float,
    *,
    show_ticks: bool,
    small_tick_interval: float | None,
    large_tick_values: list[float],
) -> list[float]:
    if not show_ticks or small_tick_interval is None or axis_max <= axis_min:
        return []

    interval = float(small_tick_interval)
    if interval <= 0:
        return []

    tick_values: list[float] = []
    tick = math.ceil((axis_min - 1e-9) / interval) * interval
    while tick <= axis_max + 1e-9 and len(tick_values) < 500:
        normalized_tick = 0.0 if abs(tick) < 1e-9 else float(tick)
        if (
            normalized_tick >= axis_min - 1e-9
            and not has_tick_value(large_tick_values, normalized_tick)
            and not has_tick_value(tick_values, normalized_tick)
        ):
            tick_values.append(normalized_tick)
        tick += interval
    return tick_values


def scaled_scalar_fraction(value: float, axis_min: float, axis_max: float) -> float:
    if axis_max <= axis_min:
        return 0.0
    return max(0.0, min(1.0, (float(value) - float(axis_min)) / (float(axis_max) - float(axis_min))))


def format_percent_tick(value: float) -> str:
    value = 0.0 if not math.isfinite(value) else float(value)
    if abs(value - round(value)) < 0.05 or abs(value) >= 100.0:
        return f"{value:.0f}%"
    if abs(value) >= 10.0:
        return f"{value:.1f}%"
    return f"{value:.2f}".rstrip("0").rstrip(".") + "%"


def linear_scalar_axis_tick_font_size_px(
    axis_config: Any,
    track_height_px: float,
) -> float:
    tick_font_size = getattr(axis_config, "tick_font_size", None)
    if tick_font_size is not None:
        return max(0.0, float(tick_font_size))
    return max(
        LINEAR_SCALAR_AXIS_DEFAULT_FONT_MIN_PX,
        min(
            LINEAR_SCALAR_AXIS_DEFAULT_FONT_MAX_PX,
            float(track_height_px) * LINEAR_SCALAR_AXIS_DEFAULT_FONT_HEIGHT_FRACTION,
        ),
    )


ScalarTickFormatter = Callable[[float], str]


__all__ = [
    "ScalarTickFormatter",
    "format_percent_tick",
    "has_tick_value",
    "linear_scalar_axis_tick_font_size_px",
    "scalar_axis_small_tick_values",
    "scalar_axis_tick_values",
    "scaled_scalar_fraction",
]
