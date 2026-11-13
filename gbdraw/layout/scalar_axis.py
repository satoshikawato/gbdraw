"""Shared scalar-track axis helpers."""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from typing import Any, NamedTuple, TypeVar, cast

import numpy as np
from numpy.typing import ArrayLike, NDArray
from pandas import DataFrame

from gbdraw.auto_sizes import scalar_axis_tick_font_size
from gbdraw.config.models.objects import ScalarAxisLabelStyle

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


def held_to_end(
    positions: ArrayLike, values: ArrayLike, end: float
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """A Linear track's samples drawn on to ``end``, the end of its axis: past the last sample the line
    holds that sample's value, the value of the last window or bin, which reaches it (OV-2871). Samples
    that already reach ``end`` are kept as they are."""

    position_array = np.asarray(positions, dtype=float)
    value_array = np.asarray(values, dtype=float)
    if not len(position_array) or position_array[-1] >= end:
        return position_array, value_array
    return np.append(position_array, float(end)), np.append(value_array, value_array[-1])


class ScalarRows:
    """One track's samples on the rows of a wrapped record: projected once, read row by row."""

    def __init__(self, transform: RecordDisplayTransform) -> None:
        self._transform = transform
        self._rows: dict[int, ProjectedScalarLine] | None = None
        self._plot_data: object | None = None

    def line(self, positions: ArrayLike, values: ArrayLike, row: int) -> ProjectedScalarLine:
        if self._rows is None:
            if self._transform.start_coordinate is None:
                # On an affine axis a sample's local position is its display offset (as in
                # ``RecordDisplayTransform.project_local_series``), so the line is held on to the drawn end
                # before it is split into rows; a row past the last sample draws that value. A display
                # start closes the periodic line at the end itself.
                positions, values = held_to_end(positions, values, self._transform.drawn_length)
            self._rows = {
                segment.row: ProjectedScalarLine(
                    np.array([point.position for point in segment.points], dtype=float),
                    np.array([point.value for point in segment.points], dtype=float),
                    (),
                )
                for segment in project_scalar_samples(positions, values, self._transform)
            }
        empty = np.array([], dtype=float)
        return self._rows.get(row, ProjectedScalarLine(empty, empty, ()))


class ScalarRow(NamedTuple):
    """The row of a wrapped record a track path draws."""

    samples: ScalarRows
    row: int


_PlotData = TypeVar("_PlotData")


def plot_data_once(row: ScalarRow | None, build: Callable[[], _PlotData]) -> _PlotData:
    """A track's plot data over the whole record: built once however many rows draw it."""

    if row is None:
        return build()
    samples = row.samples
    if samples._plot_data is None:
        samples._plot_data = build()
    return cast(_PlotData, samples._plot_data)


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


ScalarTickFormatter = Callable[[float], str]


class ScalarAxisLargeTick(NamedTuple):
    """One large tick of a scalar axis (depth, GC content in percent) as drawn: its value, its distance from
    the axis baseline along the axis (``offset_px``: a Circular radius or a Linear height above the baseline)
    and its label, or none when the axis draws no label there."""

    value: float
    offset_px: float
    label: str | None


def scalar_axis_kept_labels(positions_px: Sequence[float], line_px: float) -> tuple[int, ...]:
    """Which labels of an axis are drawn so that no two overlap, as indices into ``positions_px``: the centres
    of the labelled ticks along the axis, lowest value first. The label of the last tick (the upper bound) is
    always drawn; then each label from the bottom up whose centre is at least one label line (``line_px``, the
    label font size) from every label already drawn. Only labels are thinned; every tick is still drawn."""

    if not positions_px:
        return ()
    top = len(positions_px) - 1
    kept = [top]
    for index in range(top):
        if all(abs(positions_px[index] - positions_px[other]) >= line_px - 1e-9 for other in kept):
            kept.append(index)
    return tuple(sorted(kept))


def scalar_axis_large_ticks(
    axis_min: float,
    axis_max: float,
    *,
    show_ticks: bool,
    large_tick_interval: float | None,
    length_px: float,
    fraction: Callable[[float], float],
    format_label: ScalarTickFormatter,
    label_line_px: float,
    label_axis_min: bool = True,
) -> tuple[ScalarAxisLargeTick, ...]:
    """The large ticks of a scalar axis ``length_px`` long, each with the label it draws: the one derivation
    that the drawers, the Circular ink and the Linear left room read. ``fraction`` places a value on the axis
    (0 at the baseline, 1 at the end); ``label_axis_min`` false leaves the lower bound unlabelled."""

    values = scalar_axis_tick_values(
        axis_min, axis_max, show_ticks=show_ticks, large_tick_interval=large_tick_interval
    )
    offsets = [float(length_px) * fraction(value) for value in values]
    labelled = [
        index
        for index, value in enumerate(values)
        if label_axis_min or not math.isclose(value, axis_min, rel_tol=1e-9, abs_tol=1e-9)
    ]
    kept = {
        labelled[index]
        for index in scalar_axis_kept_labels([offsets[index] for index in labelled], label_line_px)
    }
    return tuple(
        ScalarAxisLargeTick(value, offset, format_label(value) if index in kept else None)
        for index, (value, offset) in enumerate(zip(values, offsets))
    )


class ScalarAxisLabels(NamedTuple):
    """What a depth or GC Percent axis draws besides its line: its large ticks, each numbered or not, and its
    range label (``None`` in the ``numbered`` style or without ticks)."""

    ticks: tuple[ScalarAxisLargeTick, ...]
    range_label: str | None


def _axis_labels(
    ticks: tuple[ScalarAxisLargeTick, ...], label_style: ScalarAxisLabelStyle, range_label: str
) -> ScalarAxisLabels:
    if label_style == "numbered":
        return ScalarAxisLabels(ticks, None)
    return ScalarAxisLabels(tuple(tick._replace(label=None) for tick in ticks), range_label if ticks else None)


def percent_axis_labels(
    min_percent: float,
    max_percent: float,
    *,
    label_style: ScalarAxisLabelStyle,
    show_ticks: bool,
    large_tick_interval: float | None,
    length_px: float,
    label_line_px: float,
) -> ScalarAxisLabels:
    """The large ticks and labels of a GC content axis in percent: every large tick may carry a number."""

    ticks = scalar_axis_large_ticks(
        min_percent,
        max_percent,
        show_ticks=show_ticks,
        large_tick_interval=large_tick_interval,
        length_px=length_px,
        fraction=lambda value: scaled_scalar_fraction(value, min_percent, max_percent),
        format_label=format_percent_tick,
        label_line_px=label_line_px,
    )
    return _axis_labels(ticks, label_style, scalar_axis_range_label(min_percent, max_percent, PERCENT_UNIT))


def depth_axis_labels(
    axis_min: float,
    axis_max: float,
    *,
    label_style: ScalarAxisLabelStyle,
    min_depth: float | None,
    normalize: bool,
    show_ticks: bool,
    large_tick_interval: float | None,
    length_px: float,
    label_line_px: float,
) -> ScalarAxisLabels:
    """The large ticks and labels of a depth axis, linear or log-scaled: numbered, the lower bound has no number
    without a set ``min_depth``; the range label starts at the axis minimum (0 when none is set)."""

    ticks = scalar_axis_large_ticks(
        axis_min,
        axis_max,
        show_ticks=show_ticks,
        large_tick_interval=large_tick_interval,
        length_px=length_px,
        fraction=lambda value: _scaled_depth_fraction(value, axis_min, axis_max, normalize=normalize),
        format_label=_format_depth_tick,
        label_line_px=label_line_px,
        label_axis_min=min_depth is not None,
    )
    return _axis_labels(
        ticks, label_style, scalar_axis_range_label(axis_min, axis_max, DEPTH_UNIT, log_scale=normalize)
    )


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


def _format_scalar_number(value: float) -> str:
    """A scalar axis value as its tick labels and range label print it, without the unit."""

    value = 0.0 if not math.isfinite(value) else float(value)
    if abs(value - round(value)) < 0.05 or abs(value) >= 100.0:
        return f"{value:.0f}"
    if abs(value) >= 10.0:
        return f"{value:.1f}"
    return f"{value:.2f}".rstrip("0").rstrip(".")


PERCENT_UNIT = "%"
DEPTH_UNIT = "x"


def format_percent_tick(value: float) -> str:
    return _format_scalar_number(value) + PERCENT_UNIT


def _format_depth_tick(value: float) -> str:
    return _format_scalar_number(value) + DEPTH_UNIT


def scalar_axis_range_label(axis_min: float, axis_max: float, unit: str, *, log_scale: bool = False) -> str:
    """The one label of an axis in the ``range`` style: ``[min–max<unit>]`` with an en dash, ``, log`` added on a
    log-scaled axis (Owner 2026-11-08, D-16)."""

    scale = ", log" if log_scale else ""
    return f"[{_format_scalar_number(axis_min)}\u2013{_format_scalar_number(axis_max)}{unit}{scale}]"


def _depth_axis_bounds(depth_df: DataFrame, min_depth: float | None, max_depth: float | None) -> tuple[float, float]:
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
    """A depth on the track's scale: itself, or on the log scale ``log10(depth) + 1`` (0 at 0x), so 1x, 10x and
    100x take evenly spaced heights."""

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
    depth_df: DataFrame,
    min_depth: float | None,
    max_depth: float | None,
    *,
    normalize: bool,
) -> DataFrame:
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


__all__ = [
    "plot_data_once",
    "LINEAR_SCALAR_AXIS_LABEL_GAP_PX",
    "LINEAR_SCALAR_AXIS_STROKE_PX",
    "LINEAR_SCALAR_AXIS_TICK_PX",
    "ProjectedScalarLine",
    "ScalarAxisLabels",
    "ScalarAxisLargeTick",
    "ScalarTickFormatter",
    "depth_axis_labels",
    "format_percent_tick",
    "has_tick_value",
    "linear_scalar_axis_tick_font_size_px",
    "percent_axis_labels",
    "project_scalar_line",
    "scalar_axis_kept_labels",
    "scalar_axis_large_ticks",
    "scalar_axis_range_label",
    "scalar_axis_small_tick_values",
    "scalar_axis_tick_values",
    "scaled_scalar_fraction",
]
