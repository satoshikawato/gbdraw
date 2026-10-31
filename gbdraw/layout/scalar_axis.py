"""Shared scalar-track axis helpers."""

from __future__ import annotations

import math
from collections.abc import Callable
from typing import Any, NamedTuple, cast

import numpy as np
from numpy.typing import ArrayLike, NDArray

from gbdraw.auto_sizes import scalar_axis_tick_font_size

from .record_coordinates import DisplaySeries, RecordDisplayTransform, SeriesPoint


# A Linear left-facing scalar axis (depth, GC%): its line width, its large tick length and the gap from a
# tick to its label, which is drawn left of the tick with text-anchor ``end``. The drawers and the record's
# left room (the figure bound and the definition column) read these.
LINEAR_SCALAR_AXIS_STROKE_PX = 0.8
LINEAR_SCALAR_AXIS_TICK_PX = 3.0
LINEAR_SCALAR_AXIS_LABEL_GAP_PX = 1.5


def project_scalar_samples(
    positions, values, transform: RecordDisplayTransform
) -> tuple[DisplaySeries, ...]:
    """Project one existing continuous track after its domain/value adapters.

    GC windows and Depth bins are record-local: depth_df maps source positions
    through the record coordinate map (missing input bases mean zero coverage).
    Projection never recalculates windows, normalizes values, or flips skew.
    """
    points = tuple(SeriesPoint(float(position), float(value))
                   for position, value in zip(positions, values, strict=True))
    # Every caller reaches this point only when ``transform.projects``,
    # where project_local_series returns DisplaySeries (it returns the input points otherwise).
    return cast(tuple[DisplaySeries, ...], transform.project_local_series(points))


class ProjectedScalarLine(NamedTuple):
    """Projected samples of a continuous track in display order, its drawn segments joined.

    ``segment_starts`` holds the index of the first sample of each drawn segment after the first; a track
    path keeps the ends of every segment."""

    positions: NDArray[np.float64]
    values: NDArray[np.float64]
    segment_starts: tuple[int, ...]


def project_scalar_line(
    positions: ArrayLike, values: ArrayLike, transform: RecordDisplayTransform
) -> ProjectedScalarLine:
    """:func:`project_scalar_samples` as one line: positions, values and where each drawn segment starts."""

    segments = project_scalar_samples(positions, values, transform)
    lengths = [len(segment.points) for segment in segments]
    return ProjectedScalarLine(
        np.array([point.position for segment in segments for point in segment.points], dtype=float),
        np.array([point.value for segment in segments for point in segment.points], dtype=float),
        tuple(int(start) for start in np.cumsum(lengths)[:-1]),
    )


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


def _format_depth_tick(value: float) -> str:
    value = 0.0 if not math.isfinite(value) else float(value)
    if abs(value - round(value)) < 0.05 or abs(value) >= 100.0:
        return f"{value:.0f}x"
    if abs(value) >= 10.0:
        return f"{value:.1f}x"
    return f"{value:.2f}".rstrip("0").rstrip(".") + "x"


def _depth_axis_bounds(depth_df: Any, min_depth: float | None, max_depth: float | None) -> tuple[float, float]:
    axis_min = float(min_depth) if min_depth is not None else 0.0
    if max_depth is not None:
        axis_max = float(max_depth)
    elif "depth" in depth_df.columns and not depth_df.empty:
        axis_max = float(depth_df["depth"].max())
    else:
        axis_max = axis_min
    if not math.isfinite(axis_max):
        axis_max = axis_min
    return axis_min, max(axis_min, axis_max)


def _scale_depth_value(value: float, *, normalize: bool) -> float:
    value = 0.0 if not math.isfinite(value) else float(value)
    if not normalize:
        return value
    return math.log10(value) + 1.0 if value > 0 else 0.0


def _scaled_depth_fraction(
    value: float,
    axis_min: float,
    axis_max: float,
    *,
    normalize: bool,
) -> float:
    return scaled_scalar_fraction(
        _scale_depth_value(value, normalize=normalize),
        _scale_depth_value(axis_min, normalize=normalize),
        _scale_depth_value(axis_max, normalize=normalize),
    )


def _prepare_depth_plot_dataframe(
    depth_df: Any,
    min_depth: float | None,
    max_depth: float | None,
    *,
    normalize: bool,
) -> Any:
    plot_df = depth_df.copy()
    axis_min, axis_max = _depth_axis_bounds(plot_df, min_depth, max_depth)
    if axis_max <= axis_min:
        plot_df["depth_normalized"] = 0.0
        return plot_df
    plot_df["depth_normalized"] = plot_df["depth"].map(
        lambda value: _scaled_depth_fraction(
            float(value), axis_min, axis_max, normalize=normalize
        )
    )
    return plot_df


def linear_scalar_axis_tick_font_size_px(
    axis_config: Any,
    track_height_px: float,
) -> float:
    tick_font_size = getattr(axis_config, "tick_font_size", None)
    if tick_font_size is not None:
        return max(0.0, float(tick_font_size))
    return scalar_axis_tick_font_size("linear", track_height_px)


ScalarTickFormatter = Callable[[float], str]


__all__ = [
    "LINEAR_SCALAR_AXIS_LABEL_GAP_PX",
    "LINEAR_SCALAR_AXIS_STROKE_PX",
    "LINEAR_SCALAR_AXIS_TICK_PX",
    "ProjectedScalarLine",
    "ScalarTickFormatter",
    "format_percent_tick",
    "has_tick_value",
    "linear_scalar_axis_tick_font_size_px",
    "project_scalar_line",
    "scalar_axis_small_tick_values",
    "scalar_axis_tick_values",
    "scaled_scalar_fraction",
]
