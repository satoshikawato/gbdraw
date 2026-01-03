#!/usr/bin/env python
# coding: utf-8

import math
import numpy as np
from pandas import DataFrame
from svgwrite.shapes import Circle


def generate_circle_path_desc(radius: float, norm_factor: float) -> str:
    """
    Generates the SVG path description for a circle.
    """
    norm_radius: float = radius * norm_factor
    circle_radius: float = norm_radius
    circle_start_x: float = circle_radius * math.cos(math.radians(360.0 * 0 - 90))
    circle_start_y: float = circle_radius * math.sin(math.radians(360.0 * 0 - 90))
    circle_desc: str = ""
    start_coordinate: str = "M {} {}".format(circle_start_x, circle_start_y)
    circle_desc += start_coordinate
    n_list: list[str] = []
    for n in range(0, 10, 1):
        deg: float = (n + 1) * 0.1
        n_x: float = circle_radius * math.cos(math.radians(360.0 * deg - 90))
        n_y: float = circle_radius * math.sin(math.radians(360.0 * deg - 90))
        n_list.append("A {} {} 0 0 1 {} {}".format(circle_radius, circle_radius, n_x, n_y))
    circle_desc += "".join(n_list)
    circle_desc += "z"
    return circle_desc


def generate_circular_gc_content_path_desc(
    radius: float,
    record_len: int,
    gc_df: DataFrame,
    track_width: float,
    norm_factor: float,
    dinucleotide: str,
) -> str:
    """
    Generates the SVG path description for circular GC content representation using vectorization.
    """
    norm_radius: float = radius * norm_factor

    column: str = f"{dinucleotide} content"
    mean = gc_df[column].mean()
    max_diff = (gc_df[column] - mean).abs().max()

    diff = gc_df[column] - mean

    if max_diff == 0:
        radius_of_coordinate = np.full(len(gc_df), norm_radius)
    else:
        radius_of_coordinate = norm_radius + (0.5 * track_width * (diff / max_diff))

    angles_rad = np.radians(360.0 * (gc_df.index / record_len) - 90)

    x_coords = radius_of_coordinate * np.cos(angles_rad)
    y_coords = radius_of_coordinate * np.sin(angles_rad)

    start_x = radius_of_coordinate[0] * math.cos(math.radians(-90))
    start_y = radius_of_coordinate[0] * math.sin(math.radians(-90))
    path_segments = [f"M{start_x} {start_y}"]

    path_segments.extend([f"L{x} {y}" for x, y in zip(x_coords, y_coords)])

    gc_desc: str = "".join(path_segments)
    gc_desc += "z"
    circle_desc: str = generate_circle_path_desc(radius, norm_factor)
    gc_desc += circle_desc
    return gc_desc


def generate_circular_gc_skew_path_desc(
    radius: float,
    df: DataFrame,
    total_len: int,
    track_width: float,
    norm_factor: float,
    dinucleotide: str,
) -> str:
    norm_radius: float = radius * norm_factor

    column: str = f"{dinucleotide} skew"
    mean = df[column].mean()
    max_diff = (df[column] - mean).abs().max()

    diff = df[column] - mean
    if max_diff == 0:
        radius_of_coordinate = np.full(len(df), norm_radius)
    else:
        radius_of_coordinate = norm_radius + (0.5 * track_width * (diff / max_diff))

    angles_rad = np.radians(360.0 * (df.index / total_len) - 90)

    x_coords = radius_of_coordinate * np.cos(angles_rad)
    y_coords = radius_of_coordinate * np.sin(angles_rad)

    skew_start_x = radius_of_coordinate[0] * math.cos(math.radians(-90))
    skew_start_y = radius_of_coordinate[0] * math.sin(math.radians(-90))

    path_segments = [f"M{skew_start_x} {skew_start_y}"]
    path_segments.extend([f"L{x} {y}" for x, y in zip(x_coords, y_coords)])

    skew_desc: str = "".join(path_segments)
    skew_desc += "z"
    circle_desc: str = generate_circle_path_desc(radius, norm_factor)
    skew_desc += circle_desc
    return skew_desc


def draw_circle_path(radius: float, stroke_color: str, stroke_width: float) -> Circle:
    """
    Draws a circle path for the circular canvas.
    """
    circle_path = Circle(
        center=(0, 0),
        r=radius,  # type: ignore
        stroke=stroke_color,
        stroke_width=stroke_width,
        fill="none",
    )
    return circle_path


__all__ = [
    "draw_circle_path",
    "generate_circle_path_desc",
    "generate_circular_gc_content_path_desc",
    "generate_circular_gc_skew_path_desc",
]


