#!/usr/bin/env python
# coding: utf-8

import math
from typing import Dict, Tuple, Union

from .arrows import set_arrow_shoulder
from ..layout.circular import calculate_feature_position_factors_circular


def generate_circular_intron_path(
    radius: float,
    coord_dict: Dict[str, Union[str, int]],
    total_length: int,
    track_ratio: float,
    cds_ratio: float,
    offset: float,
    track_type: str,
    strandedness: bool,
) -> list[str]:
    """
    Generates the SVG path description for an intron feature on a circular canvas.
    """
    coord_strand: str = str(coord_dict["coord_strand"])

    coord_start: int = int(coord_dict["coord_start"])
    coord_end: int = int(coord_dict["coord_end"])
    intron_strand_dict: dict[str, Tuple[str, str]] = {
        "positive": (" 0 0 1 ", " 0 0 0 "),
        "negative": (" 0 0 0 ", " 0 0 1 "),
    }
    params = intron_strand_dict[coord_strand]
    if coord_start > coord_end:
        param = params[1]
    else:
        param = params[0]

    factors: list[float] = calculate_feature_position_factors_circular(
        total_length, coord_strand, track_ratio, cds_ratio, offset, track_type, strandedness
    )
    start_x_1: float = (radius * factors[1]) * math.cos(
        math.radians(360.0 * ((coord_start) / total_length) - 90)
    )
    start_y_1: float = (radius * factors[1]) * math.sin(
        math.radians(360.0 * ((coord_start) / total_length) - 90)
    )
    end_x_1: float = (radius * factors[1]) * math.cos(
        math.radians(360.0 * ((coord_end) / total_length) - 90)
    )
    end_y_1: float = (radius * factors[1]) * math.sin(
        math.radians(360.0 * ((coord_end) / total_length) - 90)
    )
    feature_path: str = (
        "M "
        + str(start_x_1)
        + ","
        + str(start_y_1)
        + "A"
        + str(radius)
        + ","
        + str(radius)
        + param
        + str(end_x_1)
        + ","
        + str(end_y_1)
    )
    return ["line", feature_path]


def generate_circular_arrowhead_path(
    radius: float,
    coord_dict: Dict[str, Union[str, int]],
    total_length: int,
    cds_arrow_length: float,
    track_ratio: float,
    cds_ratio: float,
    offset: float,
    track_type: str,
    strandedness: bool,
) -> list[str]:
    """
    Generates the SVG path description for an arrowhead feature on a circular canvas.
    """
    coord_strand: str = str(coord_dict["coord_strand"])
    factors: list[float] = calculate_feature_position_factors_circular(
        total_length, coord_strand, track_ratio, cds_ratio, offset, track_type, strandedness
    )
    coord_start = int(coord_dict["coord_start"])
    coord_end = int(coord_dict["coord_end"])
    coord_len: float = abs(int(int(coord_dict["coord_end"]) - int(coord_dict["coord_start"])))
    coord_strand = str(coord_dict["coord_strand"])
    arrow_strand_dict: Dict[str, Tuple[int, int, str, str]] = {
        "positive": (coord_start, coord_end, " 0 0 1 ", " 0 0 0 "),
        "negative": (coord_end, coord_start, " 0 0 0 ", " 0 0 1 "),
    }
    arrow_start, arrow_end, param_1, param_2 = arrow_strand_dict[coord_strand]

    if abs(coord_len) < cds_arrow_length:
        point_x: float = (radius * factors[1]) * math.cos(
            math.radians(360.0 * (arrow_end / total_length) - 90)
        )
        point_y: float = (radius * factors[1]) * math.sin(
            math.radians(360.0 * (arrow_end / total_length) - 90)
        )
        start_x_1: float = (radius * factors[0]) * math.cos(
            math.radians(360.0 * (arrow_start / total_length) - 90)
        )
        start_y_1: float = (radius * factors[0]) * math.sin(
            math.radians(360.0 * (arrow_start / total_length) - 90)
        )
        start_x_2: float = (radius * factors[2]) * math.cos(
            math.radians(360.0 * (arrow_start / total_length) - 90)
        )
        start_y_2: float = (radius * factors[2]) * math.sin(
            math.radians(360.0 * (arrow_start / total_length) - 90)
        )
        feature_path: str = (
            "M "
            + str(start_x_1)
            + ","
            + str(start_y_1)
            + " L"
            + str(point_x)
            + ","
            + str(point_y)
            + " L"
            + str(start_x_2)
            + ","
            + str(start_y_2)
            + " z"
        )
    else:
        point_x = (radius * factors[1]) * math.cos(
            math.radians(360.0 * (arrow_end / total_length) - 90)
        )
        point_y = (radius * factors[1]) * math.sin(
            math.radians(360.0 * (arrow_end / total_length) - 90)
        )
        start_x_1 = (radius * factors[0]) * math.cos(
            math.radians(360.0 * (arrow_start / total_length) - 90)
        )
        start_y_1 = (radius * factors[0]) * math.sin(
            math.radians(360.0 * (arrow_start / total_length) - 90)
        )
        start_x_2 = (radius * factors[2]) * math.cos(
            math.radians(360.0 * (arrow_start / total_length) - 90)
        )
        start_y_2 = (radius * factors[2]) * math.sin(
            math.radians(360.0 * (arrow_start / total_length) - 90)
        )

        shoulder: float = set_arrow_shoulder(coord_strand, arrow_end, cds_arrow_length)

        end_x_1: float = (radius * factors[0]) * math.cos(
            math.radians(360.0 * ((shoulder) / total_length) - 90)
        )
        end_y_1: float = (radius * factors[0]) * math.sin(
            math.radians(360.0 * ((shoulder) / total_length) - 90)
        )
        end_x_2: float = (radius * factors[2]) * math.cos(
            math.radians(360.0 * ((shoulder) / total_length) - 90)
        )
        end_y_2: float = (radius * factors[2]) * math.sin(
            math.radians(360.0 * ((shoulder) / total_length) - 90)
        )

        arc_len_bp = abs(coord_len) - cds_arrow_length
        angle_deg = 360.0 * arc_len_bp / total_length

        if angle_deg > 20:
            if coord_strand == "positive":
                mid_pos = arrow_start + (shoulder - arrow_start) / 2
            else:
                mid_pos = shoulder + (arrow_start - shoulder) / 2

            mid_x_1 = (radius * factors[0]) * math.cos(
                math.radians(360.0 * (mid_pos / total_length) - 90)
            )
            mid_y_1 = (radius * factors[0]) * math.sin(
                math.radians(360.0 * (mid_pos / total_length) - 90)
            )
            mid_x_2 = (radius * factors[2]) * math.cos(
                math.radians(360.0 * (mid_pos / total_length) - 90)
            )
            mid_y_2 = (radius * factors[2]) * math.sin(
                math.radians(360.0 * (mid_pos / total_length) - 90)
            )

            sweep_flag_1 = param_1.strip().split(" ")[2]
            segment_param_1 = f" 0 0 {sweep_flag_1} "
            sweep_flag_2 = param_2.strip().split(" ")[2]
            segment_param_2 = f" 0 0 {sweep_flag_2} "

            outer_arc_path = (
                f"A{radius * factors[0]},{radius * factors[0]}{segment_param_1}{mid_x_1},{mid_y_1} "
                f"A{radius * factors[0]},{radius * factors[0]}{segment_param_1}{end_x_1},{end_y_1}"
            )
            inner_arc_path = (
                f"A{radius * factors[2]},{radius * factors[2]}{segment_param_2}{mid_x_2},{mid_y_2} "
                f"A{radius * factors[2]},{radius * factors[2]}{segment_param_2}{start_x_2},{start_y_2}"
            )
            feature_path = (
                f"M {start_x_1},{start_y_1} {outer_arc_path} L {point_x},{point_y} "
                f"L {end_x_2},{end_y_2} {inner_arc_path} z"
            )

        else:
            feature_path = (
                "M "
                + str(start_x_1)
                + ","
                + str(start_y_1)
                + "A"
                + str(radius)
                + ","
                + str(radius)
                + param_1
                + str(end_x_1)
                + ","
                + str(end_y_1)
                + " L"
                + str(point_x)
                + ","
                + str(point_y)
                + " L"
                + str(end_x_2)
                + ","
                + str(end_y_2)
                + "A"
                + str(radius)
                + ","
                + str(radius)
                + param_2
                + str(start_x_2)
                + ","
                + str(start_y_2)
                + " z"
            )

    return ["block", feature_path]


def generate_circular_rectangle_path(
    radius: float,
    coord_dict: Dict[str, Union[str, int]],
    total_length: int,
    track_ratio: float,
    cds_ratio,
    offset,
    track_type: str,
    strandedness: bool,
) -> list[str]:
    """
    Generates the SVG path description for a rectangular feature on a circular canvas.
    """
    coord_strand: str = str(coord_dict["coord_strand"])
    factors: list[float] = calculate_feature_position_factors_circular(
        total_length, coord_strand, track_ratio, cds_ratio, offset, track_type, strandedness
    )
    coord_start: int = int(coord_dict["coord_start"])
    coord_end: int = int(coord_dict["coord_end"])

    rect_strand_dict = {
        "positive": (coord_start, coord_end, " 0 0 1 ", " 0 0 0 "),
        "negative": (coord_end, coord_start, " 0 0 0 ", " 0 0 1 "),
    }

    rect_start, rect_end, param_1, param_2 = rect_strand_dict[coord_strand]

    coord_len_bp = coord_end - coord_start if coord_end >= coord_start else (coord_end - coord_start) + total_length
    angle_deg = 360.0 * coord_len_bp / total_length

    start_x_1: float = (radius * factors[0]) * math.cos(
        math.radians(360.0 * (rect_start / total_length) - 90)
    )
    start_y_1: float = (radius * factors[0]) * math.sin(
        math.radians(360.0 * (rect_start / total_length) - 90)
    )
    start_x_2: float = (radius * factors[2]) * math.cos(
        math.radians(360.0 * (rect_start / total_length) - 90)
    )
    start_y_2: float = (radius * factors[2]) * math.sin(
        math.radians(360.0 * (rect_start / total_length) - 90)
    )
    end_x_1: float = (radius * factors[0]) * math.cos(
        math.radians(360.0 * ((rect_end) / total_length) - 90)
    )
    end_y_1: float = (radius * factors[0]) * math.sin(
        math.radians(360.0 * ((rect_end) / total_length) - 90)
    )
    end_x_2: float = (radius * factors[2]) * math.cos(
        math.radians(360.0 * ((rect_end) / total_length) - 90)
    )
    end_y_2: float = (radius * factors[2]) * math.sin(
        math.radians(360.0 * ((rect_end) / total_length) - 90)
    )

    outer_arc_path = ""
    inner_arc_path = ""

    if angle_deg > 20:
        mid_pos = (coord_start + coord_len_bp / 2) % total_length
        mid_x_1 = (radius * factors[0]) * math.cos(
            math.radians(360.0 * (mid_pos / total_length) - 90)
        )
        mid_y_1 = (radius * factors[0]) * math.sin(
            math.radians(360.0 * (mid_pos / total_length) - 90)
        )
        mid_x_2 = (radius * factors[2]) * math.cos(
            math.radians(360.0 * (mid_pos / total_length) - 90)
        )
        mid_y_2 = (radius * factors[2]) * math.sin(
            math.radians(360.0 * (mid_pos / total_length) - 90)
        )

        sweep_flag_1 = param_1.strip().split(" ")[2]
        segment_param_1 = f" 0 0 {sweep_flag_1} "
        sweep_flag_2 = param_2.strip().split(" ")[2]
        segment_param_2 = f" 0 0 {sweep_flag_2} "

        outer_arc_path = (
            f"A{radius * factors[0]},{radius * factors[0]}{segment_param_1}{mid_x_1},{mid_y_1} "
            f"A{radius * factors[0]},{radius * factors[0]}{segment_param_1}{end_x_1},{end_y_1}"
        )
        inner_arc_path = (
            f"A{radius * factors[2]},{radius * factors[2]}{segment_param_2}{mid_x_2},{mid_y_2} "
            f"A{radius * factors[2]},{radius * factors[2]}{segment_param_2}{start_x_2},{start_y_2}"
        )
        feature_path: str = (
            f"M {start_x_1},{start_y_1} {outer_arc_path} L {end_x_2},{end_y_2} {inner_arc_path} z"
        )
    else:
        feature_path = (
            "M "
            + str(start_x_1)
            + ","
            + str(start_y_1)
            + "A"
            + str(radius)
            + ","
            + str(radius)
            + param_1
            + str(end_x_1)
            + ","
            + str(end_y_1)
            + " L"
            + str(end_x_2)
            + ","
            + str(end_y_2)
            + "A"
            + str(radius)
            + ","
            + str(radius)
            + param_2
            + str(start_x_2)
            + ","
            + str(start_y_2)
            + " z"
        )

    return ["block", feature_path]


__all__ = [
    "generate_circular_arrowhead_path",
    "generate_circular_intron_path",
    "generate_circular_rectangle_path",
]


