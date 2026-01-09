#!/usr/bin/env python
# coding: utf-8

from typing import Optional

from svgwrite.container import Group
from svgwrite.shapes import Line
from svgwrite.text import Text

from ....core.text import calculate_bbox_dimensions
from ....config.models import GbdrawConfig  # type: ignore[reportMissingImports]


class LengthBarGroup:
    """
    Handles the creation and display of a length scale in a linear layout.
    Supports two styles: 'bar' (a single scale bar) and 'ruler' (a full-axis ruler).
    """

    def __init__(
        self,
        fig_width: int,
        alignment_width: float,
        longest_genome: int,
        config_dict: dict,
        canvas_config: dict,
        group_id="length_bar",
        cfg: GbdrawConfig | None = None,
    ) -> None:
        """
        Initializes the LengthBarGroup with the given parameters.
        """
        cfg = cfg or GbdrawConfig.from_dict(config_dict)
        self._cfg = cfg

        # --- 1. Load all settings from the 'scale' section of the config ---
        scale_config = cfg.objects.scale
        self.canvas_config = canvas_config
        self.length_param = self.canvas_config.length_param
        self.stroke_color = scale_config.stroke_color
        self.stroke_width = scale_config.stroke_width
        self.font_size = scale_config.font_size.for_length_param(self.length_param)
        self.font_weight = scale_config.font_weight
        self.font_family = cfg.objects.text.font_family
        self.style = scale_config.style
        self.manual_interval = scale_config.interval
        self.scale_group_width: float = 0
        self.scale_group_height: float = 0
        self.dpi = cfg.canvas.dpi
        # --- 2. Set other properties ---
        self.longest_genome: int = longest_genome
        self.alignment_width: float = alignment_width
        self.group_id: str = group_id

        self.scale_group = Group(id=self.group_id)

        # --- 3. Branch to the appropriate setup method based on style ---
        if self.style == "ruler":
            self.setup_scale_ruler()
        else:
            self.setup_scale_bar()

    def setup_scale_bar(self) -> None:
        """
        Sets up the 'bar' style.
        Uses manual_interval if provided, otherwise calculates an automatic interval.
        """
        if self.manual_interval is not None and self.manual_interval > 0:
            # Use the user-provided interval for the bar length
            self.tick = self.manual_interval
            self.label_text = self._format_tick_label(self.manual_interval, for_bar_style=True)
        else:
            # Automatically determine the best interval and label
            self.define_ticks_by_length_for_bar()

        # Configure geometry and add elements to the SVG group
        self.config_bar_geometry()
        self.add_bar_elements_to_group()

    def setup_scale_ruler(self) -> None:
        """
        Sets up the 'ruler' style.
        Uses manual_interval if provided, otherwise falls back to an automatic interval.
        """
        # Determine the tick interval
        if self.manual_interval is not None and self.manual_interval > 0:
            tick_interval = self.manual_interval
        else:
            self.define_ticks_by_length_for_bar()  # Use the same auto-logic to get a sensible interval
            tick_interval = self.tick

        if tick_interval <= 0 or self.longest_genome <= 0:
            return

        # Draw the main axis line
        main_axis = Line(
            start=(0, 0),
            end=(self.alignment_width, 0),
            stroke=self.stroke_color,
            stroke_width=self.stroke_width,
        )
        self.scale_group.add(main_axis)
        scale_ruler_length = self.alignment_width

        # Draw ticks and labels from 0 up to the end of the genome
        position = 0
        first_tick_bbox_width = 0
        last_tick_bbox_width = 0
        max_tick_bbox_height = 0
        while True:
            x_pos = (position / self.longest_genome) * self.alignment_width

            # Draw tick line
            tick_line = Line(
                start=(x_pos, 0 - (0.5 * self.stroke_width)),
                end=(x_pos, 10),
                stroke=self.stroke_color,
                stroke_width=self.stroke_width,
            )
            self.scale_group.add(tick_line)

            # Draw label
            label_text = self._format_tick_label(position)
            text_element = Text(
                label_text,
                insert=(x_pos, 15),
                stroke="none",
                fill="black",
                font_size=self.font_size,
                font_weight=self.font_weight,
                font_family=self.font_family,
                text_anchor="middle",
                dominant_baseline="hanging",
            )
            self.scale_group.add(text_element)
            bbox_width, bbox_height = calculate_bbox_dimensions(
                label_text, self.font_family, self.font_size, self.dpi
            )
            max_tick_bbox_height = max(max_tick_bbox_height, bbox_height)
            if first_tick_bbox_width == 0:
                first_tick_bbox_width = bbox_width
            if position >= self.longest_genome:
                last_tick_bbox_width = bbox_width
                self.scale_group_width = (
                    (0.5 * first_tick_bbox_width) + scale_ruler_length + (0.5 * last_tick_bbox_width)
                )
                self.scale_group_height: float = (0.5 * self.stroke_width) + 15 + max_tick_bbox_height
                break  # Exit after drawing the last tick

            position += tick_interval
            # Ensure the final tick is exactly at the end
            if position > self.longest_genome:
                position = self.longest_genome

    def _format_tick_label(self, position: int, for_bar_style: bool = False) -> str:
        """
        Formats a position number into a readable label (e.g., "0.5M", "10kbp").
        """
        if position == 0:
            return "0"

        unit = "bp"

        if position >= 1_000_000:
            label = f"{position / 1_000_000:.1f} " + ("M" + unit)
        elif not for_bar_style and position >= 100_000:  # For ruler style, use M for values like 500k
            if self.longest_genome >= 1_000_000:
                label = f"{position / 1_000_000:.1f} " + "Mbp"
            else:
                label = f"{position // 1_000} k" + (unit)
        elif position >= 1_000:
            label = f"{position // 1_000} k" + (unit)
        else:
            label = f"{position} {unit}".strip()

        return label

    def define_ticks_by_length_for_bar(self) -> None:
        """
        Automatically determines a sensible tick interval for the 'bar' style
        or as a fallback for the 'ruler' style.
        """
        thresholds = [
            (2000, 100),
            (20000, 1000),
            (50000, 5000),
            (150000, 10000),
            (250000, 50000),
            (1000000, 100000),
            (2000000, 200000),
            (5000000, 500000),
            (float("inf"), 1000000),
        ]
        for threshold, tick_val in thresholds:
            if self.longest_genome < threshold:
                self.tick = tick_val
                self.label_text = self._format_tick_label(tick_val, for_bar_style=True)
                return

    def config_bar_geometry(self) -> None:
        """Configures the geometry (length and position) for the 'bar' style."""
        if self.longest_genome > 0:
            self.bar_length = self.alignment_width * (self.tick / self.longest_genome)
        else:
            self.bar_length = 0
        self.end_x = self.alignment_width
        self.start_x = self.end_x - self.bar_length
        self.start_y = 0
        self.end_y = 0

    def create_scale_path_linear(self) -> Line:
        """Creates the SVG line element for the 'bar' style."""
        return Line(
            start=(self.start_x, self.start_y),
            end=(self.end_x, self.end_y),
            stroke=self.stroke_color,
            stroke_width=self.stroke_width,
            fill="none",
        )

    def create_scale_text_linear(self) -> Text:
        """Creates the SVG text element for the 'bar' style label."""
        self.tick_bbox_width, self.tick_bbox_height = calculate_bbox_dimensions(
            self.label_text, self.font_family, self.font_size, self.dpi
        )
        self.scale_group_width = self.tick_bbox_width + 10 + self.bar_length
        self.scale_group_height = self.tick_bbox_height
        return Text(
            self.label_text,
            insert=(self.start_x - 10, self.start_y),
            stroke="none",
            fill="black",
            font_size=self.font_size,
            font_weight=self.font_weight,
            font_family=self.font_family,
            text_anchor="end",
            dominant_baseline="middle",
        )

    def add_bar_elements_to_group(self) -> None:
        """Adds the line and text elements for the 'bar' style to the SVG group."""
        scale_path = self.create_scale_path_linear()
        scale_text = self.create_scale_text_linear()

        self.scale_group.add(scale_path)
        self.scale_group.add(scale_text)

    def get_group(self) -> Group:
        """Retrieves the final SVG group containing the scale visualization."""
        return self.scale_group


__all__ = ["LengthBarGroup"]


