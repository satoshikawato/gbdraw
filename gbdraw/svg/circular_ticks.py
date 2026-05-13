#!/usr/bin/env python
# coding: utf-8

import math
from dataclasses import dataclass
from typing import Literal

from svgwrite.path import Path
from svgwrite.text import Text, TextPath

from ..core.text import calculate_bbox_dimensions


_TICK_LABEL_MIN_MARGIN_PX = 2.0
_TICK_LABEL_MARGIN_FONT_FACTOR = 0.15


@dataclass(frozen=True)
class CircularTickLabelGeometry:
    path_radius_px: float
    radial_inner_px: float
    radial_outer_px: float
    bbox_width_px: float
    bbox_height_px: float
    text_anchor: str = "middle"
    dominant_baseline: str = "middle"


def _tick_label_margin_px(font_size: float, tick_width: float) -> float:
    return max(_TICK_LABEL_MIN_MARGIN_PX, float(font_size) * _TICK_LABEL_MARGIN_FONT_FACTOR) + (
        max(0.0, float(tick_width)) / 2.0
    )


def get_circular_tick_intervals(total_len: int, manual_interval: int | None = None) -> tuple[int, int]:
    """Return (large, small) tick intervals for circular scales."""
    if manual_interval is not None and manual_interval > 0:
        tick_large = int(manual_interval)
        tick_small = int(manual_interval) // 10
        return tick_large, tick_small

    if total_len <= 30000:
        return 1000, 100
    if 30000 < total_len <= 50000:
        return 5000, 1000
    if 50000 < total_len <= 150000:
        return 10000, 1000
    if 150000 < total_len <= 1000000:
        return 50000, 10000
    if 1000000 < total_len <= 10000000:
        return 500000, 100000
    return 1000000, 200000


def _format_tick_label_text(tick: int, total_len: int) -> str:
    """Format circular tick label text with existing kbp/Mbp rules."""
    if total_len < 1000000:
        return f"{int(tick / 1000)} kbp"
    return f"{tick / 1000000} Mbp"


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


def _resolve_tick_track_channel(
    total_len: int,
    tick_track_channel_override: str | None = None,
) -> Literal["short", "long"]:
    """Resolve tick ratio channel from override or record length."""
    normalized = str(tick_track_channel_override or "").strip().lower()
    if normalized == "short":
        return "short"
    if normalized == "long":
        return "long"
    return "short" if total_len < 50000 else "long"


def get_circular_tick_path_ratio_bounds(
    total_len: int,
    track_type: str,
    strandedness: bool,
    tick_track_channel_override: str | None = None,
) -> tuple[float, float]:
    """Return (min_ratio, max_ratio) across small/large tick path ratios."""
    track_channel = _resolve_tick_track_channel(
        total_len,
        tick_track_channel_override=tick_track_channel_override,
    )
    ratio_table = _tick_path_ratio_table(track_channel, track_type, strandedness)
    ratio_values = [ratio for pair in ratio_table.values() for ratio in pair]
    return min(ratio_values), max(ratio_values)


def _default_tick_length_px(
    radius: float,
    *,
    length_reference_radius_px: float | None = None,
) -> float:
    reference_radius = float(length_reference_radius_px) if length_reference_radius_px is not None else float(radius)
    return max(6.0, 0.025 * reference_radius)


def get_circular_tick_path_radius_bounds(
    center_radius_px: float,
    total_len: int,
    size: str,
    track_type: str,
    strandedness: bool,
    tick_track_channel_override: str | None = None,
    tick_side: str = "legacy",
    tick_length_px: float | None = None,
    length_reference_radius_px: float | None = None,
) -> tuple[float, float]:
    """Return the radial bounds occupied by circular tick marks."""
    normalized_side = str(tick_side or "legacy").strip().lower()
    radius = float(center_radius_px)
    if normalized_side == "legacy":
        track_channel = _resolve_tick_track_channel(
            total_len,
            tick_track_channel_override=tick_track_channel_override,
        )
        ratio = _tick_path_ratio_table(track_channel, track_type, strandedness)
        prox, dist = ratio[size]
        prox_radius = radius * prox
        dist_radius = radius * dist
    else:
        if normalized_side in {"none", ""}:
            return radius, radius
        large_len = (
            float(tick_length_px)
            if tick_length_px is not None
            else _default_tick_length_px(radius, length_reference_radius_px=length_reference_radius_px)
        )
        length_px = large_len if size == "large" else max(3.0, large_len * 0.6)
        if normalized_side == "inside":
            prox_radius = radius
            dist_radius = radius - length_px
        elif normalized_side == "outside":
            prox_radius = radius
            dist_radius = radius + length_px
        elif normalized_side == "both":
            prox_radius = radius - (length_px / 2.0)
            dist_radius = radius + (length_px / 2.0)
        else:
            prox_radius = radius
            dist_radius = radius + length_px
    inner, outer = sorted((float(prox_radius), float(dist_radius)))
    return max(0.0, inner), max(0.0, outer)


def _tick_label_ratio_table(track_channel: str, track_type: str, strandedness: bool) -> dict[str, tuple[float, float]]:
    """Return radial ratios for circular tick labels (small/large)."""
    if strandedness is True:
        if track_channel == "long":
            if track_type == "middle":
                return {"small": (0.88, 1.10), "large": (0.88, 1.13)}
            if track_type == "spreadout":
                return {"small": (0.94, 1.21), "large": (0.94, 1.24)}
            if track_type == "tuckin":
                return {"small": (0.882, 1.10), "large": (0.82, 1.13)}
            return {"small": (0.98, 1.0), "large": (0.98, 1.0)}
        if track_type == "middle":
            return {"small": (0.85, 1.09), "large": (0.85, 1.12)}
        if track_type == "spreadout":
            return {"small": (0.95, 1.21), "large": (0.95, 1.24)}
        if track_type == "tuckin":
            return {"small": (1.03, 1.21), "large": (1.03, 1.24)}
        return {"small": (0.98, 1.0), "large": (0.98, 1.0)}

    if track_channel == "long":
        if track_type == "middle":
            return {"small": (0.90, 1.12), "large": (0.90, 1.15)}
        if track_type == "spreadout":
            return {"small": (0.95, 1.21), "large": (0.95, 1.24)}
        if track_type == "tuckin":
            return {"small": (0.84, 1.14), "large": (0.84, 1.17)}
        return {"small": (0.98, 1.0), "large": (0.98, 1.0)}
    if track_type == "middle":
        return {"small": (0.89, 1.10), "large": (0.89, 1.13)}
    if track_type == "spreadout":
        return {"small": (0.96, 1.21), "large": (0.96, 1.24)}
    if track_type == "tuckin":
        return {"small": (1.03, 1.21), "large": (1.03, 1.24)}
    return {"small": (0.98, 1.0), "large": (0.98, 1.0)}


def get_circular_tick_label_radius_bounds(
    center_radius_px: float,
    total_len: int,
    track_type: str,
    strandedness: bool,
    font_size: float,
    font_family: str,
    dpi: int,
    manual_interval: int | None = None,
    tick_track_channel_override: str | None = None,
    label_side: str = "legacy",
    tick_side: str = "legacy",
    tick_length_px: float | None = None,
    tick_width: float = 0.0,
    length_reference_radius_px: float | None = None,
) -> tuple[float, float] | None:
    """Return (inner, outer) radial bounds occupied by large circular tick labels."""
    if str(label_side or "legacy").strip().lower() in {"none", ""}:
        return None
    tick_large, _ = get_circular_tick_intervals(total_len, manual_interval=manual_interval)
    if tick_large <= 0:
        return None

    ticks_large_nonzero = [tick for tick in range(0, total_len, tick_large) if tick != 0]
    if not ticks_large_nonzero:
        return None

    min_radius = float("inf")
    max_radius = float("-inf")
    for tick in ticks_large_nonzero:
        label_text = _format_tick_label_text(tick, total_len)
        geometry = resolve_circular_tick_label_geometry(
            center_radius_px=center_radius_px,
            total_len=total_len,
            size="large",
            tick=tick,
            label_text=label_text,
            font_size=font_size,
            font_family=font_family,
            track_type=track_type,
            strandedness=strandedness,
            dpi=dpi,
            tick_track_channel_override=tick_track_channel_override,
            label_side=label_side,
            tick_side=tick_side,
            tick_length_px=tick_length_px,
            tick_width=tick_width,
            length_reference_radius_px=length_reference_radius_px,
        )

        min_radius = min(min_radius, geometry.radial_inner_px)
        max_radius = max(max_radius, geometry.radial_outer_px)

    if not math.isfinite(min_radius) or not math.isfinite(max_radius):
        return None
    if max_radius < min_radius:
        min_radius, max_radius = max_radius, min_radius
    return max(0.0, min_radius), max(0.0, max_radius)


def get_circular_tick_label_radius_profile(
    total_len: int,
    track_type: str,
    strandedness: bool,
    font_size: float,
    font_family: str,
    dpi: int,
    manual_interval: int | None = None,
    tick_track_channel_override: str | None = None,
) -> tuple[float, float] | None:
    """Return the legacy large-label radius ratio and radial text extent."""
    tick_large, _ = get_circular_tick_intervals(total_len, manual_interval=manual_interval)
    if tick_large <= 0:
        return None

    ticks_large_nonzero = [tick for tick in range(0, total_len, tick_large) if tick != 0]
    if not ticks_large_nonzero:
        return None

    track_channel = _resolve_tick_track_channel(
        total_len,
        tick_track_channel_override=tick_track_channel_override,
    )
    ratio_table = _tick_label_ratio_table(track_channel, track_type, strandedness)
    label_radius_ratio, _ = ratio_table["large"]

    max_extent_px = 0.0
    for tick in ticks_large_nonzero:
        label_text = _format_tick_label_text(tick, total_len)
        _, bbox_height_px = calculate_bbox_dimensions(label_text, font_family, font_size, dpi)
        max_extent_px = max(max_extent_px, float(bbox_height_px) / 4.0)

    return float(label_radius_ratio), max(0.0, max_extent_px)


def _raw_tick_label_base_radius(
    *,
    center_radius_px: float,
    total_len: int,
    size: str,
    track_type: str,
    strandedness: bool,
    font_size: float,
    tick_track_channel_override: str | None,
    label_side: str,
    tick_length_px: float | None,
    length_reference_radius_px: float | None,
) -> float | None:
    normalized_side = str(label_side or "legacy").strip().lower()
    radius = float(center_radius_px)
    if normalized_side in {"none", ""}:
        return None
    if normalized_side == "legacy":
        track_channel = _resolve_tick_track_channel(
            total_len,
            tick_track_channel_override=tick_track_channel_override,
        )
        ratio = _tick_label_ratio_table(track_channel, track_type, strandedness)
        prox, _ = ratio[size]
        return radius * prox

    large_len = (
        float(tick_length_px)
        if tick_length_px is not None
        else _default_tick_length_px(radius, length_reference_radius_px=length_reference_radius_px)
    )
    offset = large_len + float(font_size) * 0.9
    if normalized_side == "inside":
        return radius - offset
    return radius + offset


def _separate_label_path_radius_from_tick_band(
    *,
    path_radius: float,
    label_inner_offset_px: float,
    label_outer_offset_px: float,
    tick_inner_px: float,
    tick_outer_px: float,
    label_side: str,
    font_size: float,
    tick_width: float,
) -> float:
    normalized_side = str(label_side or "legacy").strip().lower()
    tick_inner, tick_outer = sorted((float(tick_inner_px), float(tick_outer_px)))
    tick_midpoint = (tick_inner + tick_outer) / 2.0
    desired_side = normalized_side
    if desired_side not in {"inside", "outside"}:
        desired_side = "inside" if float(path_radius) <= tick_midpoint else "outside"

    margin_px = _tick_label_margin_px(font_size, tick_width)
    inner = float(path_radius) - float(label_inner_offset_px)
    outer = float(path_radius) + float(label_outer_offset_px)
    adjusted = float(path_radius)
    if desired_side == "inside":
        max_outer = tick_inner - margin_px
        if outer > max_outer:
            adjusted -= outer - max_outer
    else:
        min_inner = tick_outer + margin_px
        if inner < min_inner:
            adjusted += min_inner - inner
    return max(1.0, adjusted)


def _tick_label_radial_offsets(
    *,
    bbox_height_px: float,
    dominant_baseline: str,
    reverse_path: bool,
) -> tuple[float, float]:
    """Return (inner_offset, outer_offset) from the text path radius."""
    bbox_height = max(0.0, float(bbox_height_px))
    baseline = str(dominant_baseline or "middle").strip().lower()
    if baseline == "middle":
        half_height = bbox_height / 2.0
        return half_height, half_height
    if baseline == "hanging":
        return (0.0, bbox_height) if reverse_path else (bbox_height, 0.0)
    if baseline == "text-after-edge":
        return (bbox_height, 0.0) if reverse_path else (0.0, bbox_height)
    half_height = bbox_height / 2.0
    return half_height, half_height


def resolve_circular_tick_label_geometry(
    *,
    center_radius_px: float,
    total_len: int,
    size: str,
    tick: int,
    label_text: str,
    font_size: float,
    font_family: str,
    track_type: str,
    strandedness: bool,
    dpi: int,
    tick_track_channel_override: str | None = None,
    label_side: str = "legacy",
    tick_side: str = "legacy",
    tick_length_px: float | None = None,
    tick_width: float = 0.0,
    length_reference_radius_px: float | None = None,
) -> CircularTickLabelGeometry:
    """Resolve one circular tick label path and its radial text footprint."""
    label_radius_base = _raw_tick_label_base_radius(
        center_radius_px=center_radius_px,
        total_len=total_len,
        size=size,
        track_type=track_type,
        strandedness=strandedness,
        font_size=font_size,
        tick_track_channel_override=tick_track_channel_override,
        label_side=label_side,
        tick_length_px=tick_length_px,
        length_reference_radius_px=length_reference_radius_px,
    )
    if label_radius_base is None:
        return CircularTickLabelGeometry(0.0, 0.0, 0.0, 0.0, 0.0)

    bbox_width_px, bbox_height_px = calculate_bbox_dimensions(label_text, font_family, font_size, dpi)
    _anchor_value, dominant_baseline = set_tick_label_anchor_value(total_len, tick)
    angle = 360.0 * (tick / total_len)
    reverse_path = 90 <= angle < 270
    path_radius = float(label_radius_base)
    label_inner_offset_px, label_outer_offset_px = _tick_label_radial_offsets(
        bbox_height_px=float(bbox_height_px),
        dominant_baseline=dominant_baseline,
        reverse_path=reverse_path,
    )
    tick_inner, tick_outer = get_circular_tick_path_radius_bounds(
        center_radius_px=center_radius_px,
        total_len=total_len,
        size=size,
        track_type=track_type,
        strandedness=strandedness,
        tick_track_channel_override=tick_track_channel_override,
        tick_side=tick_side,
        tick_length_px=tick_length_px,
        length_reference_radius_px=length_reference_radius_px,
    )
    if str(tick_side or "legacy").strip().lower() not in {"none", ""}:
        path_radius = _separate_label_path_radius_from_tick_band(
            path_radius=path_radius,
            label_inner_offset_px=label_inner_offset_px,
            label_outer_offset_px=label_outer_offset_px,
            tick_inner_px=tick_inner,
            tick_outer_px=tick_outer,
            label_side=label_side,
            font_size=font_size,
            tick_width=tick_width,
        )
    return CircularTickLabelGeometry(
        path_radius_px=path_radius,
        radial_inner_px=max(0.0, path_radius - label_inner_offset_px),
        radial_outer_px=max(0.0, path_radius + label_outer_offset_px),
        bbox_width_px=float(bbox_width_px),
        bbox_height_px=float(bbox_height_px),
        text_anchor="middle",
        dominant_baseline=dominant_baseline,
    )


def generate_circular_tick_paths(
    radius: float,
    total_len: int,
    size: str,
    ticks: list,
    tick_width: float,
    track_type: str,
    strandedness: bool,
    tick_track_channel_override: str | None = None,
    tick_side: str = "legacy",
    tick_length_px: float | None = None,
    length_reference_radius_px: float | None = None,
) -> list[Path]:
    """
    Generates SVG path descriptions for tick marks on a circular canvas.
    """
    tick_paths_list: list[Path] = []
    if str(tick_side or "legacy").strip().lower() in {"none", ""}:
        return []
    prox_radius, dist_radius = get_circular_tick_path_radius_bounds(
        center_radius_px=radius,
        total_len=total_len,
        size=size,
        track_type=track_type,
        strandedness=strandedness,
        tick_track_channel_override=tick_track_channel_override,
        tick_side=tick_side,
        tick_length_px=tick_length_px,
        length_reference_radius_px=length_reference_radius_px,
    )
    for tick in ticks:
        prox_x: float = prox_radius * math.cos(math.radians(360.0 * (tick / total_len) - 90))
        prox_y: float = prox_radius * math.sin(math.radians(360.0 * (tick / total_len) - 90))
        dist_x: float = dist_radius * math.cos(math.radians(360.0 * (tick / total_len) - 90))
        dist_y: float = dist_radius * math.sin(math.radians(360.0 * (tick / total_len) - 90))
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


def _tick_label_path_range(
    *,
    tick: float,
    label_as_feature_length: float,
    reverse_path: bool,
    text_anchor: str,
) -> tuple[float, float, str]:
    label_length = float(label_as_feature_length)
    anchor = str(text_anchor or "middle").strip().lower()
    if anchor == "start":
        start_tick = float(tick)
        end_tick = float(tick) - label_length if reverse_path else float(tick) + label_length
        return start_tick, end_tick, "0%"
    if anchor == "end":
        start_tick = float(tick) + label_length if reverse_path else float(tick) - label_length
        end_tick = float(tick)
        return start_tick, end_tick, "100%"
    half_length = label_length / 2.0
    if reverse_path:
        return float(tick) + half_length, float(tick) - half_length, "50%"
    return float(tick) - half_length, float(tick) + half_length, "50%"


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
    tick_track_channel_override: str | None = None,
    label_side: str = "legacy",
    tick_side: str = "legacy",
    tick_length_px: float | None = None,
    tick_width: float = 0.0,
    length_reference_radius_px: float | None = None,
) -> list[Text]:
    tick_label_paths_list: list[Text] = []
    normalized_side = str(label_side or "legacy").strip().lower()
    if normalized_side in {"none", ""}:
        return []
    for tick in ticks:
        angle = 360.0 * (tick / total_len)
        label_text = _format_tick_label_text(tick, total_len)

        geometry = resolve_circular_tick_label_geometry(
            center_radius_px=radius,
            total_len=total_len,
            size=size,
            tick=tick,
            label_text=label_text,
            font_size=font_size,
            font_family=font_family,
            track_type=track_type,
            strandedness=strandedness,
            dpi=dpi,
            tick_track_channel_override=tick_track_channel_override,
            label_side=label_side,
            tick_side=tick_side,
            tick_length_px=tick_length_px,
            tick_width=tick_width,
            length_reference_radius_px=length_reference_radius_px,
        )
        path_radius = geometry.path_radius_px
        label_as_feature_length = total_len * 1.5 * geometry.bbox_width_px / (2 * math.pi * max(1.0, path_radius))
        reverse_path = 90 <= angle < 270
        label_start, label_end, start_offset = _tick_label_path_range(
            tick=tick,
            label_as_feature_length=label_as_feature_length,
            reverse_path=reverse_path,
            text_anchor=geometry.text_anchor,
        )
        if not reverse_path:
            param = " 0 0 1 "
            start_x_1: float = path_radius * math.cos(
                math.radians(360.0 * (label_start / total_len) - 90)
            )
            start_y_1: float = path_radius * math.sin(
                math.radians(360.0 * (label_start / total_len) - 90)
            )
            end_x: float = path_radius * math.cos(
                math.radians(360.0 * (label_end / total_len) - 90)
            )
            end_y: float = path_radius * math.sin(
                math.radians(360.0 * (label_end / total_len) - 90)
            )
        else:
            param = " 1 0 0 "
            start_x_1 = path_radius * math.cos(
                math.radians(360.0 * (label_start / total_len) - 90)
            )
            start_y_1 = path_radius * math.sin(
                math.radians(360.0 * (label_start / total_len) - 90)
            )
            end_x = path_radius * math.cos(
                math.radians(360.0 * (label_end / total_len) - 90)
            )
            end_y = path_radius * math.sin(
                math.radians(360.0 * (label_end / total_len) - 90)
            )
        label_axis_path_desc: str = (
            "M "
            + str(start_x_1)
            + ","
            + str(start_y_1)
            + "A"
            + str(path_radius)
            + ","
            + str(path_radius)
            + param
            + str(end_x)
            + ","
            + str(end_y)
        )
        label_axis_path = Path(d=label_axis_path_desc, stroke="none", fill="none")
        text_path = Text("", dominant_baseline=geometry.dominant_baseline)
        text_path.add(
            TextPath(
                label_axis_path,
                text=label_text,
                startOffset=start_offset,
                method="align",
                text_anchor=geometry.text_anchor,
                font_size=font_size,
                font_style="normal",
                font_weight="normal",
                font_family=font_family,
                dominant_baseline=geometry.dominant_baseline,
            )
        )
        tick_label_paths_list.append(label_axis_path)
        tick_label_paths_list.append(text_path)
    return tick_label_paths_list


__all__ = [
    "CircularTickLabelGeometry",
    "get_circular_tick_intervals",
    "get_circular_tick_label_radius_bounds",
    "get_circular_tick_label_radius_profile",
    "get_circular_tick_path_radius_bounds",
    "get_circular_tick_path_ratio_bounds",
    "generate_circular_tick_labels",
    "generate_circular_tick_paths",
    "resolve_circular_tick_label_geometry",
    "set_tick_label_anchor_value",
]


