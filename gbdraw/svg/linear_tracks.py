#!/usr/bin/env python
# coding: utf-8

from typing import NamedTuple

import numpy as np
from numpy.typing import ArrayLike, NDArray

from ..layout.linear_coords import normalize_position_to_linear_track


from gbdraw.layout.record_coordinates import RecordDisplayTransform
from ..layout.scalar_axis import ProjectedScalarLine, ScalarRow, SeriesSamples, held_to_end, project_scalar_line
from .path_data import path_line_to, path_number
from .track_paths import track_line_to, track_polyline

def _samples(
    positions: ArrayLike,
    values: ArrayLike,
    record_len: int,
    record_transform: RecordDisplayTransform | None,
    row: ScalarRow | None,
) -> ProjectedScalarLine:
    """The samples a Linear track path draws, in display order (projected when the record is transformed;
    of one row, local to it, on a wrapped record), on to the end of the ``record_len`` bases drawn."""

    if row is not None:
        return row.samples.line(positions, values, row.row)
    if record_transform is not None and record_transform.projects:
        # A display start closes the periodic line at the end (``RecordDisplayTransform.project_series``).
        return project_scalar_line(positions, values, record_transform)
    return ProjectedScalarLine(*held_to_end(positions, values, record_len), ())


class LinearTrackFrame(NamedTuple):
    """Where a Linear numeric series track draws: its top left corner, the bases and width of what it draws (one
    row of a wrapped record), and its height."""

    start_x: float
    start_y: float
    record_len: int
    alignment_width: float
    track_height: float


def _area_path(start_x: float, base_y: float, xs: NDArray[np.float64], ys: NDArray[np.float64],
               segment_starts: tuple[int, ...]) -> str:
    """The area between a track's line and its baseline at ``base_y``: from the baseline at ``start_x`` along the
    line, down to the baseline and back."""

    return (
        f"M{path_number(start_x)} {path_number(base_y)}"
        + track_line_to(xs, ys, segment_starts)
        + path_line_to([xs[-1], start_x], [base_y, base_y])
        + "z"
    )


def linear_series_baseline_y(samples: SeriesSamples, frame: LinearTrackFrame) -> float:
    """The y of the baseline a numeric series path closes on, and below which the part below the baseline is
    painted: ``start_y`` for a deviation series, else the band's bottom raised by the band fraction
    ``samples.baseline``."""

    if samples.deviation is not None:
        return float(frame.start_y)
    return float(frame.start_y) + float(frame.track_height) - (float(frame.track_height) * samples.baseline)


def _line_desc(xs: NDArray[np.float64], ys: NDArray[np.float64], segment_starts: tuple[int, ...]) -> str:
    """The open line of a ``line`` series through its samples."""

    return track_polyline(xs, ys, segment_starts)


def linear_series_path_desc(
    samples: SeriesSamples,
    frame: LinearTrackFrame,
    record_transform: RecordDisplayTransform | None = None,
    row: ScalarRow | None = None,
) -> str:
    """The path of a numeric series track. A deviation series draws each sample's deviation from the mean of the
    series' windows with a value, closed to the baseline at ``start_y``; a window without a value lies on the
    baseline, so the area breaks there (OV-2799). A band series draws each value as the fraction of the height it
    fills from the bottom of the track, closed to its baseline (``linear_series_baseline_y``). A ``line`` series
    is its line alone, open."""

    start_x, start_y, record_len, alignment_width, track_height = frame
    scale = samples.deviation
    if scale is not None:
        positions, values, segment_starts = _samples(
            samples.positions, samples.values.fillna(scale.mean), record_len, record_transform, row
        )
        if not len(positions):
            return ""
        xs = normalize_position_to_linear_track(positions, record_len, alignment_width)
        max_diff = scale.max_diff
        # A completely flat series draws a flat line (no division by zero).
        ys = np.zeros_like(values) if max_diff == 0 else -(0.5 * track_height * ((values - scale.mean) / max_diff))
        if samples.style == "line":
            return _line_desc(xs, ys, segment_starts)
        return _area_path(start_x, start_y, xs, ys, segment_starts)
    if record_len <= 0 or track_height <= 0:
        return ""
    bottom_y = float(start_y) + float(track_height)
    positions, values, segment_starts = _samples(
        # Unprojected positions are whole bases.
        samples.positions if record_transform is not None and record_transform.projects
        else samples.positions.to_numpy().astype(np.int64),
        samples.values,
        record_len,
        record_transform,
        row,
    )
    if not len(positions):
        return ""
    xs = normalize_position_to_linear_track(positions, record_len, alignment_width)
    ys = bottom_y - (float(track_height) * np.clip(values, 0.0, 1.0))
    if samples.style == "line":
        return _line_desc(xs, ys, segment_starts)
    return _area_path(start_x, linear_series_baseline_y(samples, frame), xs, ys, segment_starts)


__all__ = ["LinearTrackFrame", "linear_series_baseline_y", "linear_series_path_desc"]
