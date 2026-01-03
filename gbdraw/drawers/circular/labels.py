#!/usr/bin/env python
# coding: utf-8

import math
from typing import Literal, Tuple

from svgwrite.container import Group  # type: ignore[reportMissingImports]
from svgwrite.path import Path  # type: ignore[reportMissingImports]
from svgwrite.text import Text, TextPath  # type: ignore[reportMissingImports]

from ...config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from ...core.sequence import determine_length_parameter  # type: ignore[reportMissingImports]
from ...layout.circular import calculate_feature_position_factors_circular  # type: ignore[reportMissingImports]
from ...layout.common import calculate_cds_ratio  # type: ignore[reportMissingImports]
from ...svg.text_path import generate_text_path  # type: ignore[reportMissingImports]
from ...text import calculate_bbox_dimensions  # type: ignore[reportMissingImports]


class LabelDrawer:
    def __init__(self, config_dict: dict, cfg: GbdrawConfig | None = None) -> None:
        self.config_dict = config_dict
        cfg = cfg or GbdrawConfig.from_dict(config_dict)
        self._cfg = cfg
        self.strandedness: bool = cfg.canvas.strandedness

    def set_feature_label_anchor_value(
        self,
        total_len: int,
        tick: float,
        start_x: float,
        is_inner: bool = False,
    ) -> Tuple[Literal["middle", "start", "end"], Literal["text-after-edge", "middle", "hanging"]]:
        angle = (360.0 * (tick / total_len)) % 360

        if start_x > 0:
            anchor_value = "end" if is_inner else "start"
        else:
            anchor_value = "start" if is_inner else "end"

        if is_inner:
            if 0 <= angle < 10:
                baseline_value = "hanging"
            elif 10 <= angle < 170:
                baseline_value = "middle"
            elif 170 <= angle < 190:
                baseline_value = "middle"
            elif 190 <= angle < 350:
                baseline_value = "middle"
            else:
                baseline_value = "hanging"
        else:
            if 0 <= angle < 10:
                baseline_value = "text-after-edge"
            elif 10 <= angle < 170:
                baseline_value = "middle"
            elif 170 <= angle < 190:
                baseline_value = "hanging"
            elif 190 <= angle < 350:
                baseline_value = "middle"
            else:
                baseline_value = "text-after-edge"

        return anchor_value, baseline_value

    def embed_label(self, group, label, radius, record_length, track_ratio):
        length_param = determine_length_parameter(record_length, self._cfg.labels.length_threshold.circular)
        track_ratio_factor = self._cfg.canvas.circular.track_ratio_factors[length_param][0]
        cds_ratio, offset = calculate_cds_ratio(track_ratio, length_param, track_ratio_factor)
        factors: list[float] = calculate_feature_position_factors_circular(
            record_length, label["strand"], track_ratio, cds_ratio, offset, self.track_type, self.strandedness
        )
        angle = 360.0 * (label["middle"] / record_length)
        font_px = float(str(self.font_size).rstrip("ptpx"))
        _, bbox_h = calculate_bbox_dimensions(label["label_text"], self.font_family, font_px, dpi=96)
        center_offset = bbox_h / 4

        if 0 <= angle < 90:
            param = " 0 0 1 "
            start_x_1: float = (radius * factors[1] - center_offset) * math.cos(
                math.radians(360.0 * (label["start"] / record_length) - 90)
            )
            start_y_1: float = (radius * factors[1] - center_offset) * math.sin(
                math.radians(360.0 * (label["start"] / record_length) - 90)
            )
            end_x: float = (radius * factors[1] - center_offset) * math.cos(
                math.radians(360.0 * ((label["end"]) / record_length) - 90)
            )
            end_y: float = (radius * factors[1] - center_offset) * math.sin(
                math.radians(360.0 * ((label["end"]) / record_length) - 90)
            )
            label_axis_path_desc: str = (
                "M " + str(start_x_1) + "," + str(start_y_1) + "A" + str(radius - center_offset) + "," + str(radius - center_offset) + param + str(end_x) + "," + str(end_y)
            )
        if 90 <= angle < 270:
            param = " 1 0 0 "
            start_x_1 = (radius * factors[1] + center_offset) * math.cos(
                math.radians(360.0 * ((label["end"]) / record_length) - 90)
            )
            start_y_1 = (radius * factors[1] + center_offset) * math.sin(
                math.radians(360.0 * ((label["end"]) / record_length) - 90)
            )
            end_x = (radius * factors[1] + center_offset) * math.cos(
                math.radians(360.0 * (label["start"] / record_length) - 90)
            )
            end_y = (radius * factors[1] + center_offset) * math.sin(
                math.radians(360.0 * (label["start"] / record_length) - 90)
            )
            label_axis_path_desc = (
                "M " + str(start_x_1) + "," + str(start_y_1) + "A" + str(radius + center_offset) + "," + str(radius + center_offset) + param + str(end_x) + "," + str(end_y)
            )
        elif 270 <= angle <= 360:
            param = " 0 0 1 "
            start_x_1 = (radius * factors[1] - center_offset) * math.cos(
                math.radians(360.0 * (label["start"] / record_length) - 90)
            )
            start_y_1 = (radius * factors[1] - center_offset) * math.sin(
                math.radians(360.0 * (label["start"] / record_length) - 90)
            )
            end_x = (radius * factors[1] - center_offset) * math.cos(
                math.radians(360.0 * ((label["end"]) / record_length) - 90)
            )
            end_y = (radius * factors[1] - center_offset) * math.sin(
                math.radians(360.0 * ((label["end"]) / record_length) - 90)
            )
            label_axis_path_desc = (
                "M " + str(start_x_1) + "," + str(start_y_1) + "A" + str(radius - center_offset) + "," + str(radius - center_offset) + param + str(end_x) + "," + str(end_y)
            )

        label_axis_path = Path(d=label_axis_path_desc, stroke="none", fill="none")
        text_path = Text("")
        text_path.add(
            TextPath(
                label_axis_path,
                text=label["label_text"],
                startOffset="50%",
                method="align",
                text_anchor="middle",
                font_size=self.font_size,
                font_style="normal",
                font_weight="normal",
                font_family=self.font_family,
                dominant_baseline="auto",
            )
        )
        group.add(label_axis_path)
        group.add(text_path)
        return group

    def add_label_on_the_rim(self, group, label, radius, record_length):
        anchor_value, baseline_value = self.set_feature_label_anchor_value(
            record_length, label["middle"], label["start_x"], label.get("is_inner", False)
        )
        start_x_1 = label["start_x"]
        start_y_1 = label["start_y"]
        label_path = generate_text_path(
            label["label_text"],
            start_x_1,
            start_y_1,
            interval=0,
            font_size=self.font_size,
            font_weight="normal",
            font=self.font_family,
            dominant_baseline=baseline_value,
            text_anchor=anchor_value,
        )
        group.add(label_path)
        return group

    def draw(self, label, group, record_length, radius, track_ratio):
        cfg = self._cfg
        length_threshold = cfg.labels.length_threshold.circular
        length_param = determine_length_parameter(record_length, length_threshold)
        self.font_size = cfg.labels.font_size.for_length_param(length_param)
        self.font_family = cfg.objects.text.font_family
        self.track_type = cfg.canvas.circular.track_type
        self.strandedness = cfg.canvas.strandedness
        if label["is_embedded"] is True:
            group = self.embed_label(group, label, radius, record_length, track_ratio)
        else:
            group = self.add_label_on_the_rim(group, label, radius, record_length)
        return group


__all__ = ["LabelDrawer"]


