#!/usr/bin/env python
# coding: utf-8

import math
from typing import NamedTuple

import numpy as np
from numpy.typing import NDArray
from pandas import DataFrame
from gbdraw.svg.elements import Circle


from gbdraw.layout.record_coordinates import RecordDisplayTransform
from ..layout.scalar_axis import project_scalar_line
from .circular_arcs import circle_arc, circle_arc_to, circle_point, position_angle
from .path_data import path_point
from .track_paths import track_line_to, track_polyline


class CircularTrackLine(NamedTuple):
    """The line a Circular track path draws: display positions, their radii, and where each drawn segment
    after the first starts."""

    positions: NDArray[np.float64]
    radii: NDArray[np.float64]
    segment_starts: tuple[int, ...] = ()


class CircularAreaLine(NamedTuple):
    """The outer line of a Circular scalar area and the baseline radius the area is filled down to."""

    positions: NDArray[np.float64]
    radii: NDArray[np.float64]
    baseline_radius: float
    segment_starts: tuple[int, ...] = ()

def generate_circle_path_desc(radius: float, norm_factor: float) -> str:
    """The closed circle of ``radius * norm_factor`` from 12 o'clock, clockwise (two arcs)."""
    return circle_arc(radius * norm_factor, -90.0, 360.0) + "z"


def circular_gc_content_line(
    radius: float,
    record_len: int,
    gc_df: DataFrame,
    track_width: float,
    norm_factor: float,
    dinucleotide: str,
    record_transform: RecordDisplayTransform | None = None,
) -> CircularTrackLine | None:
    """The GC content line the path draws: display positions and their radii; the area lies between it and the
    circle at ``radius * norm_factor``."""

    if gc_df.empty or record_len <= 0:
        return None

    norm_radius: float = radius * norm_factor

    column: str = f"{dinucleotide} content"
    mean = gc_df[column].mean()
    max_diff = (gc_df[column] - mean).abs().max()

    diff = gc_df[column] - mean

    if max_diff == 0:
        radius_of_coordinate = np.full(len(gc_df), norm_radius)
    else:
        radius_of_coordinate = (norm_radius + (0.5 * track_width * (diff / max_diff))).to_numpy(dtype=float)

    positions = gc_df.index.to_numpy(dtype=float)
    if record_transform is not None and record_transform.start_coordinate is not None:
        return CircularTrackLine(*project_scalar_line(positions, radius_of_coordinate, record_transform))
    if len(positions) > 0 and np.isclose(positions[0], 0.0):
        positions = np.append(positions, float(record_len))
        radius_of_coordinate = np.append(radius_of_coordinate, radius_of_coordinate[0])
    return CircularTrackLine(positions, np.asarray(radius_of_coordinate, dtype=float))


def generate_circular_gc_content_path_desc(
    radius: float,
    record_len: int,
    gc_df: DataFrame,
    track_width: float,
    norm_factor: float,
    dinucleotide: str,
    record_transform: RecordDisplayTransform | None = None,
) -> str:
    """
    Generates the SVG path description for circular GC content representation using vectorization.
    """
    line = circular_gc_content_line(
        radius, record_len, gc_df, track_width, norm_factor, dinucleotide, record_transform
    )
    if line is None:
        return ""
    return _closed_line_desc(line, float(record_len)) + generate_circle_path_desc(radius, norm_factor)


def _closed_line_desc(line: CircularTrackLine, total_len: float) -> str:
    """The closed line of a GC content (deviation) or skew path: it starts at 12 o'clock on the radius of the
    first sample, runs through every sample, and closes; the baseline circle follows as its own subpath."""

    positions, radii, segment_starts = line
    angles_rad = np.radians(360.0 * (positions / total_len) - 90)
    x_coords = np.concatenate(([radii[0] * math.cos(math.radians(-90))], radii * np.cos(angles_rad)))
    y_coords = np.concatenate(([radii[0] * math.sin(math.radians(-90))], radii * np.sin(angles_rad)))
    return track_polyline(x_coords, y_coords, [start + 1 for start in segment_starts]) + "z"


def circular_gc_skew_line(
    radius: float,
    df: DataFrame,
    total_len: int,
    track_width: float,
    norm_factor: float,
    dinucleotide: str,
    record_transform: RecordDisplayTransform | None = None,
) -> CircularTrackLine | None:
    """The skew line the path draws: display positions and their radii; the areas lie between it and the circle
    at ``radius * norm_factor``."""

    if df.empty or total_len <= 0:
        return None
    norm_radius: float = radius * norm_factor

    column: str = f"{dinucleotide} skew"
    mean = df[column].mean()
    max_diff = (df[column] - mean).abs().max()

    diff = df[column] - mean
    if max_diff == 0:
        radius_of_coordinate = np.full(len(df), norm_radius)
    else:
        radius_of_coordinate = norm_radius + (0.5 * track_width * (diff / max_diff))

    positions = df.index
    if record_transform is not None and record_transform.start_coordinate is not None:
        return CircularTrackLine(*project_scalar_line(positions, radius_of_coordinate, record_transform))
    return CircularTrackLine(positions.to_numpy(dtype=float), np.asarray(radius_of_coordinate, dtype=float))


def generate_circular_gc_skew_path_desc(
    radius: float,
    df: DataFrame,
    total_len: int,
    track_width: float,
    norm_factor: float,
    dinucleotide: str,
    record_transform: RecordDisplayTransform | None = None,
) -> str:
    line = circular_gc_skew_line(radius, df, total_len, track_width, norm_factor, dinucleotide, record_transform)
    if line is None:
        return ""
    return _closed_line_desc(line, float(total_len)) + generate_circle_path_desc(radius, norm_factor)


def circular_scalar_area_line(
    radius: float,
    record_len: int,
    scalar_df: DataFrame,
    track_width: float,
    norm_factor: float,
    *,
    value_column: str = "value_normalized",
    position_column: str = "position",
    close_at_record_len: bool = False,
    record_transform: RecordDisplayTransform | None = None,
) -> CircularAreaLine | None:
    """The outer line of a scalar area as the path draws it: display positions, their radii, and the baseline
    radius the area is filled down to."""

    if scalar_df.empty or record_len <= 0 or track_width <= 0:
        return None

    baseline_radius = max(0.0, (float(radius) * float(norm_factor)) - (0.5 * float(track_width)))
    values = scalar_df[value_column].to_numpy(dtype=float)
    values = np.clip(values, 0.0, 1.0)
    positions = scalar_df[position_column].to_numpy(dtype=float)
    segment_starts: tuple[int, ...] = ()
    if record_transform is not None and record_transform.start_coordinate is not None:
        positions, values, segment_starts = project_scalar_line(positions, values, record_transform)
    elif close_at_record_len and len(positions) > 0 and np.isclose(positions[0], 0.0):
        positions = np.append(positions, float(record_len))
        values = np.append(values, values[0])
    return CircularAreaLine(positions, baseline_radius + (float(track_width) * values), baseline_radius, segment_starts)


def generate_circular_scalar_area_path_desc(
    radius: float,
    record_len: int,
    scalar_df: DataFrame,
    track_width: float,
    norm_factor: float,
    *,
    value_column: str = "value_normalized",
    position_column: str = "position",
    close_at_record_len: bool = False,
    record_transform: RecordDisplayTransform | None = None,
) -> str:
    """Return an annular filled area path for normalized scalar values."""

    line = circular_scalar_area_line(
        radius,
        record_len,
        scalar_df,
        track_width,
        norm_factor,
        value_column=value_column,
        position_column=position_column,
        close_at_record_len=close_at_record_len,
        record_transform=record_transform,
    )
    if line is None:
        return ""
    positions, outer_radii, baseline_radius, segment_starts = line
    if len(positions) == 0:
        return ""
    angles_rad = np.radians(360.0 * (positions / float(record_len)) - 90.0)
    first = position_angle(positions[0], record_len)
    last = position_angle(positions[-1], record_len)
    # The outer line, then the baseline back from the last sample to the first as arcs (Circular wrap kept).
    return (
        f"M {path_point(*circle_point(baseline_radius, first))}"
        + track_line_to(outer_radii * np.cos(angles_rad), outer_radii * np.sin(angles_rad), segment_starts)
        + f"L {path_point(*circle_point(baseline_radius, last))} "
        + circle_arc_to(baseline_radius, last, first - last)
        + "z"
    )


def generate_circular_depth_path_desc(
    radius: float,
    record_len: int,
    depth_df: DataFrame,
    track_width: float,
    norm_factor: float,
    record_transform: RecordDisplayTransform | None = None,
) -> str:
    """Return an annular filled area path for binned circular depth coverage."""

    return generate_circular_scalar_area_path_desc(
        radius,
        record_len,
        depth_df,
        track_width,
        norm_factor,
        value_column="depth_normalized",
        record_transform=record_transform,
    )


def draw_circle_path(radius: float, stroke_color: str, stroke_width: float) -> Circle:
    """
    Draws a circle path for the circular canvas.
    """
    circle_path = Circle(
        center=(0, 0),
        r=radius,
        stroke=stroke_color,
        stroke_width=stroke_width,
        fill="none",
    )
    return circle_path


__all__ = [
    "draw_circle_path",
    "generate_circular_depth_path_desc",
    "generate_circular_scalar_area_path_desc",
    "generate_circle_path_desc",
    "generate_circular_gc_content_path_desc",
    "generate_circular_gc_skew_path_desc",
]


