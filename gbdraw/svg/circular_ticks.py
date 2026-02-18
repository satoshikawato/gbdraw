#!/usr/bin/env python
# coding: utf-8

import math
from typing import Literal

from svgwrite.path import Path
from svgwrite.text import Text, TextPath

from ..core.text import calculate_bbox_dimensions


def _tick_path_ratio_table(track_channel: str, track_type: str, strandedness: bool) -> dict[str, tuple[float, float]]:
    """Return radial ratios for tick line paths (small/large)."""
    if strandedness:
        if track_channel == "long":
            if track_type == "middle":
                return {"small": (0.915, 0.93), "large": (0.91, 0.93)}
            if track_type == "spreadout":
                return {"small": (0.985, 1.0), "large": (0.98, 1.0)}
            if track_type == "tuckin":
                return {"small": (0.845, 0.86), "large": (0.84, 0.86)}
            return {"small": (1.06, 1.08), "large": (0.98, 1.0)}
        if track_type == "middle":
            return {"small": (0.885, 0.90), "large": (0.88, 0.90)}
        if track_type == "spreadout":
            return {"small": (0.985, 1.0), "large": (0.98, 1.0)}
        if track_type == "tuckin":
            return {"small": (0.985, 1.0), "large": (0.98, 1.0)}
        return {"small": (1.06, 1.08), "large": (0.98, 1.0)}

    if track_channel == "long":
        if track_type == "middle":
            return {"small": (0.945, 0.96), "large": (0.94, 0.96)}
        if track_type == "spreadout":
            return {"small": (0.985, 1.0), "large": (0.98, 1.0)}
        if track_type == "tuckin":
            return {"small": (0.885, 0.90), "large": (0.88, 0.90)}
        return {"small": (1.06, 1.08), "large": (0.98, 1.0)}
    if track_type == "middle":
        return {"small": (0.925, 0.94), "large": (0.92, 0.94)}
    if track_type == "spreadout":
        return {"small": (0.985, 1.0), "large": (0.98, 1.0)}
    if track_type == "tuckin":
        return {"small": (0.985, 1.0), "large": (0.98, 1.0)}
    return {"small": (1.06, 1.08), "large": (0.98, 1.0)}


def get_circular_tick_path_ratio_bounds(total_len: int, track_type: str, strandedness: bool) -> tuple[float, float]:
    """Return (min_ratio, max_ratio) across small/large tick path ratios."""
    track_channel = "short" if total_len < 50000 else "long"
    ratio_table = _tick_path_ratio_table(track_channel, track_type, strandedness)
    ratio_values = [ratio for pair in ratio_table.values() for ratio in pair]
    return min(ratio_values), max(ratio_values)


def generate_circular_tick_paths(
    radius: float,
    total_len: int,
    size: str,
    ticks: list,
    tick_width: float,
    track_type: str,
    strandedness: bool,
) -> list[Path]:
    """
    Generates SVG path descriptions for tick marks on a circular canvas.
    """
    tick_paths_list: list[Path] = []
    track_channel = "short" if total_len < 50000 else "long"
    ratio = _tick_path_ratio_table(track_channel, track_type, strandedness)

    prox, dist = ratio[size]
    for tick in ticks:
        prox_x: float = (radius * prox) * math.cos(math.radians(360.0 * (tick / total_len) - 90))
        prox_y: float = (radius * prox) * math.sin(math.radians(360.0 * (tick / total_len) - 90))
        dist_x: float = (radius * dist) * math.cos(math.radians(360.0 * (tick / total_len) - 90))
        dist_y: float = (radius * dist) * math.sin(math.radians(360.0 * (tick / total_len) - 90))
        tick_path_desc: str = "M " + str(prox_x) + "," + str(prox_y) + " L" + str(dist_x) + "," + str(dist_y) + " z"
        tick_path = Path(d=tick_path_desc, stroke="gray", stroke_width=tick_width)
        tick_paths_list.append(tick_path)
    return tick_paths_list


def set_tick_label_anchor_value(
    total_len: int, tick: float
) -> tuple[Literal["middle", "start", "end"], Literal["text-after-edge", "middle", "hanging"]]:
    """
    Determines the anchor and baseline values for tick labels based on their position.
    """
    angle: float = 360.0 * (tick / total_len)
    if 0 <= angle < 45:
        anchor_value, baseline_value = "middle", "text-after-edge"
    elif 45 <= angle < 155:
        anchor_value, baseline_value = "start", "middle"
    elif 155 <= angle < 205:
        anchor_value, baseline_value = "middle", "hanging"
    elif 205 <= angle < 315:
        anchor_value, baseline_value = "end", "middle"
    else:
        anchor_value, baseline_value = "middle", "text-after-edge"
    return anchor_value, baseline_value


def generate_circular_tick_labels(
    radius: float,
    total_len: int,
    size: str,
    ticks: list,
    stroke: str,
    fill: str,
    font_size: float,
    font_weight: str,
    font_family: str,
    track_type: str,
    strandedness: bool,
    dpi: int,
) -> list[Text]:
    tick_label_paths_list: list[Text] = []
    if total_len < 50000:
        track_channel = "short"
    else:
        track_channel = "long"
    if strandedness is True:
        if track_channel == "long":
            if track_type == "middle":
                ratio = {"small": [0.88, 1.10], "large": [0.88, 1.13]}
            elif track_type == "spreadout":
                ratio = {"small": [0.94, 1.21], "large": [0.94, 1.24]}
            elif track_type == "tuckin":
                ratio = {"small": [0.882, 1.10], "large": [0.82, 1.13]}
            else:
                ratio = {"small": [0.98, 1.0], "large": [0.98, 1.0]}
        else:
            if track_type == "middle":
                ratio = {"small": [0.85, 1.09], "large": [0.85, 1.12]}
            elif track_type == "spreadout":
                ratio = {"small": [0.95, 1.21], "large": [0.95, 1.24]}
            elif track_type == "tuckin":
                ratio = {"small": [1.03, 1.21], "large": [1.03, 1.24]}
            else:
                ratio = {"small": [0.98, 1.0], "large": [0.98, 1.0]}
    else:
        if track_channel == "long":
            if track_type == "middle":
                ratio = {"small": [0.90, 1.12], "large": [0.90, 1.15]}
            elif track_type == "spreadout":
                ratio = {"small": [0.95, 1.21], "large": [0.95, 1.24]}
            elif track_type == "tuckin":
                ratio = {"small": [0.84, 1.14], "large": [0.84, 1.17]}
            else:
                ratio = {"small": [0.98, 1.0], "large": [0.98, 1.0]}
        else:
            if track_type == "middle":
                ratio = {"small": [0.89, 1.10], "large": [0.89, 1.13]}
            elif track_type == "spreadout":
                ratio = {"small": [0.96, 1.21], "large": [0.96, 1.24]}
            elif track_type == "tuckin":
                ratio = {"small": [1.03, 1.21], "large": [1.03, 1.24]}
            else:
                ratio = {"small": [0.98, 1.0], "large": [0.98, 1.0]}
    prox, dist = ratio[size]
    for tick in ticks:
        anchor_value, baseline_value = set_tick_label_anchor_value(total_len, tick)

        angle = 360.0 * (tick / total_len)
        if total_len < 1000000:
            label_text: str = str(int(tick / 1000)) + " kbp"
        else:
            label_text = str((tick / 1000000)) + " Mbp"

        bbox_width_px, bbox_height_px = calculate_bbox_dimensions(label_text, font_family, font_size, dpi)
        center_offset = bbox_height_px / 4
        label_as_feature_length = total_len * 1.5 * bbox_width_px / (2 * math.pi * radius)
        label_start = tick - (label_as_feature_length / 2)
        label_end = tick + (label_as_feature_length / 2)
        if 0 <= angle < 90:
            param = " 0 0 1 "
            start_x_1: float = (radius * prox - center_offset) * math.cos(
                math.radians(360.0 * (label_start / total_len) - 90)
            )
            start_y_1: float = (radius * prox - center_offset) * math.sin(
                math.radians(360.0 * (label_start / total_len) - 90)
            )
            end_x: float = (radius * prox - center_offset) * math.cos(
                math.radians(360.0 * (label_end / total_len) - 90)
            )
            end_y: float = (radius * prox - center_offset) * math.sin(
                math.radians(360.0 * (label_end / total_len) - 90)
            )
        if 90 <= angle < 270:
            param = " 1 0 0 "
            start_x_1 = (radius * prox + center_offset) * math.cos(
                math.radians(360.0 * (label_end / total_len) - 90)
            )
            start_y_1 = (radius * prox + center_offset) * math.sin(
                math.radians(360.0 * (label_end / total_len) - 90)
            )
            end_x = (radius * prox + center_offset) * math.cos(
                math.radians(360.0 * (label_start / total_len) - 90)
            )
            end_y = (radius * prox + center_offset) * math.sin(
                math.radians(360.0 * (label_start / total_len) - 90)
            )
        elif 270 <= angle <= 360:
            param = " 0 0 1 "
            start_x_1 = (radius * prox - center_offset) * math.cos(
                math.radians(360.0 * (label_start / total_len) - 90)
            )
            start_y_1 = (radius * prox - center_offset) * math.sin(
                math.radians(360.0 * (label_start / total_len) - 90)
            )
            end_x = (radius * prox - center_offset) * math.cos(
                math.radians(360.0 * (label_end / total_len) - 90)
            )
            end_y = (radius * prox - center_offset) * math.sin(
                math.radians(360.0 * (label_end / total_len) - 90)
            )
        label_axis_path_desc: str = (
            "M "
            + str(start_x_1)
            + ","
            + str(start_y_1)
            + "A"
            + str(radius)
            + ","
            + str(radius)
            + param
            + str(end_x)
            + ","
            + str(end_y)
        )
        label_axis_path = Path(d=label_axis_path_desc, stroke="none", fill="none")
        text_path = Text("")
        text_path.add(
            TextPath(
                label_axis_path,
                text=label_text,
                startOffset="50%",
                method="align",
                text_anchor="middle",
                font_size=font_size,
                font_style="normal",
                font_weight="normal",
                font_family=font_family,
            )
        )
        tick_label_paths_list.append(label_axis_path)
        tick_label_paths_list.append(text_path)
    return tick_label_paths_list


__all__ = [
    "get_circular_tick_path_ratio_bounds",
    "generate_circular_tick_labels",
    "generate_circular_tick_paths",
    "set_tick_label_anchor_value",
]


