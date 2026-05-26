#!/usr/bin/env python
# coding: utf-8

import math

from svgwrite.container import Group
from svgwrite.gradients import LinearGradient
from svgwrite.path import Path

from ....core.text import calculate_bbox_dimensions
from ....svg.text_path import generate_text_path


class LegendGroup:
    def __init__(self, canvas_config, legend_config, legend_table):
        self.legend_group = Group(id="legend")
        self.canvas_config = canvas_config
        self.legend_config = legend_config
        self.legend_table = legend_table
        self.font_family = self.legend_config.font_family
        self.font_size = self.legend_config.font_size
        self.dpi = self.canvas_config.dpi

        self.color_rect_size = self.legend_config.color_rect_size
        self.num_of_lines = len(self.legend_table.keys())
        self.add_elements_to_group()

    def create_rectangle_path_for_legend(self) -> str:
        normalized_start: float = 0
        normalized_end: float = self.color_rect_size
        start_y_top, start_y_bottom = -0.5 * self.color_rect_size, 0.5 * self.color_rect_size
        end_y_top, end_y_bottom = -0.5 * self.color_rect_size, 0.5 * self.color_rect_size
        rectangle_path: str = (
            f"M {normalized_start},{start_y_top} L {normalized_end},{end_y_top} "
            f"L {normalized_end},{end_y_bottom} L {normalized_start},{start_y_bottom} z"
        )
        return rectangle_path

    def _gradient_entries(self) -> list[tuple[str, dict]]:
        return [
            (str(key), properties)
            for key, properties in self.legend_table.items()
            if properties.get("type") == "gradient"
        ]

    def _gradient_id(self, key: str, properties: dict) -> str:
        raw = f"{key}_{properties['min_color']}_{properties['max_color']}"
        safe = "".join(char if char.isalnum() else "_" for char in raw).strip("_")
        return f"circular_legend_grad_{safe or 'identity'}"

    def _min_label_text(self, value: object) -> str:
        min_identity = float(value or 0)
        if min_identity == int(min_identity):
            return f"{int(min_identity)}%"
        return f"{min_identity}%"

    def _build_compact_gradient_legend(self, gradient_entries: list[tuple[str, dict]]) -> tuple[Group, float, float]:
        group = Group(id="conservation_identity_legend")
        font = self.font_family
        bar_width = 10 * self.color_rect_size
        label_gap = 0.2 * self.color_rect_size
        label_width = max(
            calculate_bbox_dimensions(str(key), self.font_family, self.font_size, self.dpi)[0]
            for key, _ in gradient_entries
        )
        bar_x = float(label_width) + label_gap
        total_width = bar_x + bar_width
        row_height = (24 / 14) * self.color_rect_size
        path_desc = (
            f"M 0,{-self.color_rect_size / 2} L {bar_width},{-self.color_rect_size / 2} "
            f"L {bar_width},{self.color_rect_size / 2} L 0,{self.color_rect_size / 2} z"
        )

        for idx, (key, properties) in enumerate(gradient_entries):
            row_y = self.color_rect_size / 2 + idx * row_height
            gradient_id = self._gradient_id(key, properties)
            entry_group = Group(debug=False)
            entry_group.attribs["data-legend-key"] = str(key)

            gradient = LinearGradient(start=(0, 0), end=("100%", 0), id=gradient_id)
            gradient.add_stop_color(offset="0%", color=properties["min_color"])
            gradient.add_stop_color(offset="100%", color=properties["max_color"])
            entry_group.add(gradient)

            label_path = generate_text_path(
                key,
                0,
                0,
                0,
                self.font_size,
                "normal",
                font,
                dominant_baseline="central",
                text_anchor="start",
            )
            label_path.translate(0, row_y)
            entry_group.add(label_path)

            grad_rect = Path(
                d=path_desc,
                fill=f"url(#{gradient_id})",
                stroke=properties["stroke"],
                stroke_width=properties["width"],
            )
            grad_rect.translate(bar_x, row_y)
            entry_group.add(grad_rect)
            group.add(entry_group)

        min_label_text = self._min_label_text(gradient_entries[0][1].get("min_value", 0))
        scale_y = self.color_rect_size + (len(gradient_entries) - 1) * row_height + 2
        min_label = generate_text_path(
            min_label_text,
            0,
            0,
            0,
            self.font_size,
            "normal",
            font,
            dominant_baseline="hanging",
            text_anchor="start",
        )
        min_label.translate(bar_x, scale_y)
        group.add(min_label)

        max_label = generate_text_path(
            "100%",
            0,
            0,
            0,
            self.font_size,
            "normal",
            font,
            dominant_baseline="hanging",
            text_anchor="end",
        )
        max_label.translate(bar_x + bar_width, scale_y)
        group.add(max_label)

        _, min_label_height = calculate_bbox_dimensions(
            min_label_text, self.font_family, self.font_size, self.dpi
        )
        _, max_label_height = calculate_bbox_dimensions(
            "100%", self.font_family, self.font_size, self.dpi
        )
        total_height = scale_y + max(float(min_label_height), float(max_label_height))
        return group, total_width, total_height

    def _build_gradient_legend(self, gradient_entries: list[tuple[str, dict]]) -> tuple[Group, float, float]:
        if len(gradient_entries) > 1:
            return self._build_compact_gradient_legend(gradient_entries)

        group = Group(id="conservation_identity_legend")
        font = self.font_family
        bar_width = 10 * self.color_rect_size
        row_height = (24 / 14) * self.color_rect_size
        path_desc = (
            f"M 0,{-self.color_rect_size / 2} L {bar_width},{-self.color_rect_size / 2} "
            f"L {bar_width},{self.color_rect_size / 2} L 0,{self.color_rect_size / 2} z"
        )
        y_offset = 0.0
        total_width = bar_width

        for key, properties in gradient_entries:
            gradient_id = self._gradient_id(key, properties)
            entry_group = Group(debug=False)
            entry_group.attribs["data-legend-key"] = str(key)

            gradient = LinearGradient(start=(0, 0), end=("100%", 0), id=gradient_id)
            gradient.add_stop_color(offset="0%", color=properties["min_color"])
            gradient.add_stop_color(offset="100%", color=properties["max_color"])
            entry_group.add(gradient)

            title_path = generate_text_path(
                key,
                0,
                0,
                0,
                self.font_size,
                "normal",
                font,
                dominant_baseline="hanging",
                text_anchor="middle",
            )
            title_width, title_height = calculate_bbox_dimensions(
                str(key), self.font_family, self.font_size, self.dpi
            )
            title_path.translate(bar_width / 2.0, y_offset)
            entry_group.add(title_path)
            total_width = max(total_width, float(title_width))
            y_offset += float(title_height) + (self.color_rect_size / 2.0)

            grad_rect = Path(
                d=path_desc,
                fill=f"url(#{gradient_id})",
                stroke=properties["stroke"],
                stroke_width=properties["width"],
            )
            grad_rect.translate(0, y_offset)
            entry_group.add(grad_rect)

            min_label_text = self._min_label_text(properties.get("min_value", 0))
            min_label = generate_text_path(
                min_label_text,
                0,
                0,
                0,
                self.font_size,
                "normal",
                font,
                dominant_baseline="hanging",
                text_anchor="start",
            )
            min_label.translate(0, y_offset + self.color_rect_size / 2.0 + 2)
            entry_group.add(min_label)
            max_label = generate_text_path(
                "100%",
                0,
                0,
                0,
                self.font_size,
                "normal",
                font,
                dominant_baseline="hanging",
                text_anchor="end",
            )
            max_label.translate(bar_width, y_offset + self.color_rect_size / 2.0 + 2)
            entry_group.add(max_label)
            _, label_height = calculate_bbox_dimensions(
                "100%", self.font_family, self.font_size, self.dpi
            )
            y_offset += self.color_rect_size / 2.0 + 2 + float(label_height)
            group.add(entry_group)
            y_offset += row_height * 0.35

        return group, total_width, max(0.0, y_offset)

    def add_elements_to_group(self):
        path_desc = (
            f"M {0},{-0.5 * self.color_rect_size} "
            f"L {self.legend_config.legend_width},{-0.5 * self.color_rect_size} "
            f"L {self.legend_config.legend_width},{self.legend_config.legend_height -0.5 * self.color_rect_size} "
            f"L {0},{self.legend_config.legend_height -0.5 * self.color_rect_size} z"
        )
        rect_path = Path(d=path_desc, fill="none", stroke="none", stroke_width=0)
        self.legend_group.add(rect_path)
        path_desc = self.create_rectangle_path_for_legend()
        font = self.font_family
        line_margin = (24 / 14) * self.color_rect_size
        x_margin = (22 / 14) * self.color_rect_size

        horizontal_layout = self.canvas_config.legend_position in {"top", "bottom"}
        gradient_entries = self._gradient_entries()
        gradient_width = float(self.legend_config.pairwise_legend_width) if gradient_entries else 0.0
        if horizontal_layout:
            wrap_width = float(self.legend_config.legend_width)
            if gradient_entries:
                wrap_width = max(self.color_rect_size, wrap_width - gradient_width - x_margin)
            if wrap_width <= 0:
                wrap_width = float("inf")

            solid_entries: list[tuple[str, dict, float]] = []
            for key, properties in self.legend_table.items():
                if properties.get("type") != "solid":
                    continue
                text_width, _ = calculate_bbox_dimensions(
                    str(key), self.font_family, self.font_size, self.dpi
                )
                entry_width = float(text_width) + (2 * x_margin)
                solid_entries.append((str(key), properties, entry_width))

            rows: list[list[tuple[str, dict, float]]] = []
            current_row: list[tuple[str, dict, float]] = []
            current_row_width = 0.0
            for key, properties, entry_width in solid_entries:
                exceeds_wrap = (
                    math.isfinite(wrap_width)
                    and current_row
                    and ((current_row_width + x_margin + entry_width) > wrap_width)
                )
                if exceeds_wrap:
                    rows.append(current_row)
                    current_row = []
                    current_row_width = 0.0
                current_row.append((key, properties, entry_width))
                current_row_width += entry_width
            if current_row:
                rows.append(current_row)

            for row_index, row_entries in enumerate(rows):
                row_y = row_index * line_margin
                row_width = sum(float(entry[2]) for entry in row_entries)
                if math.isfinite(wrap_width):
                    row_start_x = x_margin + max(0.0, (wrap_width - row_width) * 0.5)
                else:
                    row_start_x = x_margin
                current_x = row_start_x

                for key, properties, entry_width in row_entries:
                    # Create entry group with data attribute for identification
                    entry_group = Group(debug=False)
                    entry_group.attribs["data-legend-key"] = str(key)

                    rect_path = Path(
                        d=path_desc,
                        fill=properties["fill"],
                        stroke=properties["stroke"],
                        stroke_width=properties["width"],
                    )
                    rect_path.translate(current_x - x_margin, row_y)
                    entry_group.add(rect_path)
                    legend_path = generate_text_path(
                        key, 0, 0, 0, self.font_size, "normal", font, dominant_baseline="central", text_anchor="start"
                    )
                    legend_path.translate(current_x, row_y)
                    entry_group.add(legend_path)
                    self.legend_group.add(entry_group)
                    current_x += float(entry_width)
            if gradient_entries:
                gradient_group, gradient_group_width, gradient_height = self._build_gradient_legend(gradient_entries)
                gradient_group.translate(
                    max(wrap_width + x_margin, self.legend_config.legend_width - gradient_group_width),
                    max(0.0, (self.legend_config.legend_height - gradient_height) / 2.0),
                )
                self.legend_group.add(gradient_group)
        else:
            count = 0
            for key, properties in self.legend_table.items():
                if properties.get("type") == "solid":
                    # Create entry group with data attribute for identification
                    entry_group = Group(debug=False)
                    entry_group.attribs["data-legend-key"] = str(key)

                    rect_path = Path(
                        d=path_desc,
                        fill=properties["fill"],
                        stroke=properties["stroke"],
                        stroke_width=properties["width"],
                    )
                    rect_path.translate(0, count * line_margin)
                    entry_group.add(rect_path)
                    legend_path = generate_text_path(
                        key, 0, 0, 0, self.font_size, "normal", font, dominant_baseline="central", text_anchor="start"
                    )
                    legend_path.translate(x_margin, count * line_margin)
                    entry_group.add(legend_path)
                    self.legend_group.add(entry_group)
                    count += 1
            if gradient_entries:
                gradient_group, gradient_group_width, _ = self._build_gradient_legend(gradient_entries)
                gradient_group.translate(
                    max(0.0, (float(self.legend_config.legend_width) - float(gradient_group_width)) / 2.0),
                    count * line_margin + (line_margin * 0.5 if count else 0.0),
                )
                self.legend_group.add(gradient_group)
        return self.legend_group

    def get_group(self) -> Group:
        return self.legend_group


__all__ = ["LegendGroup"]
