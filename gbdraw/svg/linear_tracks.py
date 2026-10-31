#!/usr/bin/env python
# coding: utf-8

import numpy as np
from numpy.typing import ArrayLike
from pandas import DataFrame, Series

from ..layout.linear_coords import normalize_position_to_linear_track


from gbdraw.layout.record_coordinates import RecordDisplayTransform
from ..layout.scalar_axis import ProjectedScalarLine, project_scalar_line
from .path_data import path_line_to, path_number
from .track_paths import track_line_to

def _samples(
    positions: ArrayLike, values: ArrayLike, record_transform: RecordDisplayTransform | None
) -> ProjectedScalarLine:
    """The samples a Linear track path draws, in display order (projected when the record is transformed)."""

    if record_transform is not None and record_transform.projects:
        return project_scalar_line(positions, values, record_transform)
    return ProjectedScalarLine(np.asarray(positions, dtype=float), np.asarray(values, dtype=float), ())


def _deviation_path_desc(
    start_x: float,
    start_y: float,
    series: Series,
    max_diff: float,
    record_len: int,
    alignment_width: float,
    track_height: float,
    record_transform: RecordDisplayTransform | None,
) -> str:
    """A GC content (deviation) or skew path: the line of each sample's deviation from the mean, closed to
    the baseline at ``start_y``."""

    positions, values, segment_starts = _samples(series.index, series, record_transform)
    xs = normalize_position_to_linear_track(positions, record_len, alignment_width)
    # A completely flat series draws a flat line (no division by zero).
    ys = np.zeros_like(values) if max_diff == 0 else -(0.5 * track_height * ((values - float(series.mean())) / max_diff))
    return (
        f"M{path_number(start_x)} {path_number(start_y)}"
        + track_line_to(xs, ys, segment_starts)
        + path_line_to([xs[-1] if len(xs) else start_x, start_x], [start_y, start_y])
        + "z"
    )


def calculate_gc_content_path_desc(
    start_x: float,
    start_y: float,
    gc_df: DataFrame,
    record_len: int,
    alignment_width: float,
    track_height: float,
    dinucleotide: str,
    record_transform: RecordDisplayTransform | None = None,
) -> str:
    if gc_df.empty:
        return ""
    series = gc_df[f"{dinucleotide} content"]
    max_diff = float((series - float(series.mean())).abs().max())
    return _deviation_path_desc(
        start_x, start_y, series, max_diff, record_len, alignment_width, track_height, record_transform
    )


def calculate_gc_skew_path_desc(
    start_x: float,
    start_y: float,
    skew_df: DataFrame,
    record_len: int,
    alignment_width: float,
    track_height: float,
    record_transform: RecordDisplayTransform | None = None,
) -> str:
    if skew_df.empty:
        return ""
    column = [col for col in skew_df.columns if "skew" in col and "cumulative" not in col.lower()][0]
    series = skew_df[column]
    max_diff = float((series - float(series.mean())).abs().max()) or 1.0
    return _deviation_path_desc(
        start_x, start_y, series, max_diff, record_len, alignment_width, track_height, record_transform
    )


def calculate_linear_scalar_area_path_desc(
    start_x: float,
    start_y: float,
    scalar_df: DataFrame,
    record_len: int,
    alignment_width: float,
    track_height: float,
    *,
    value_column: str = "value_normalized",
    position_column: str = "position",
    record_transform: RecordDisplayTransform | None = None,
) -> str:
    """Return a filled linear area path for normalized scalar values."""

    if scalar_df.empty or record_len <= 0 or track_height <= 0:
        return ""

    baseline_y = float(start_y) + float(track_height)
    positions, values, segment_starts = _samples(
        # Unprojected positions are whole bases.
        scalar_df[position_column] if record_transform is not None and record_transform.projects
        else scalar_df[position_column].to_numpy().astype(np.int64),
        scalar_df[value_column],
        record_transform,
    )
    if not len(positions):
        return ""
    x_values = normalize_position_to_linear_track(positions, record_len, alignment_width)
    y_values = baseline_y - (float(track_height) * np.clip(values, 0.0, 1.0))

    final_x = normalize_position_to_linear_track(
        record_len, record_len, alignment_width
    )
    path_segments = [
        f"M{path_number(start_x)} {path_number(baseline_y)}",
        track_line_to(x_values, y_values, segment_starts),
        path_line_to([final_x, final_x, start_x], [y_values[-1], baseline_y, baseline_y]),
        "z",
    ]
    return "".join(path_segments)


def calculate_depth_path_desc(
    start_x: float,
    start_y: float,
    depth_df: DataFrame,
    record_len: int,
    alignment_width: float,
    track_height: float,
    record_transform: RecordDisplayTransform | None = None,
) -> str:
    """Return a filled linear area path for binned depth coverage."""

    return calculate_linear_scalar_area_path_desc(
        start_x,
        start_y,
        depth_df,
        record_len,
        alignment_width,
        track_height,
        value_column="depth_normalized",
        record_transform=record_transform,
    )


__all__ = [
    "calculate_depth_path_desc",
    "calculate_gc_content_path_desc",
    "calculate_gc_skew_path_desc",
    "calculate_linear_scalar_area_path_desc",
]


