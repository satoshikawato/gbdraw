#!/usr/bin/env python
# coding: utf-8

from typing import Mapping

from svgwrite.container import Group
from svgwrite.gradients import LinearGradient
from svgwrite.path import Path

from ....legend.circular_layout import CircularGradientLegendLayout, build_circular_legend_layout
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
        canvas_width = float(
            getattr(
                self.canvas_config,
                "total_width",
                getattr(self.legend_config, "legend_width", 0.0),
            )
        )
        configured_width = float(getattr(self.legend_config, "legend_width", 0.0) or 0.0)
        configured_pairwise_width = float(
            getattr(self.legend_config, "pairwise_legend_width", 0.0) or 0.0
        )
        self.layout = build_circular_legend_layout(
            self.legend_table,
            legend_position=str(getattr(self.canvas_config, "legend_position", "right")),
            canvas_width=canvas_width,
            font_family=self.font_family,
            font_size=float(self.font_size),
            dpi=int(self.dpi),
            color_rect_size=float(self.color_rect_size),
            legend_width=configured_width if configured_width > 0 else None,
            pairwise_legend_width=(
                configured_pairwise_width if configured_pairwise_width > 0 else None
            ),
        )
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

    def _gradient_id(self, key: str, properties: Mapping[str, object]) -> str:
        raw = f"{key}_{properties['min_color']}_{properties['max_color']}"
        safe = "".join(char if char.isalnum() else "_" for char in raw).strip("_")
        return f"circular_legend_grad_{safe or 'identity'}"

    def _build_gradient_legend(self, layout: CircularGradientLegendLayout) -> Group:
        group = Group(id="conservation_identity_legend")
        font = self.font_family
        path_desc = (
            f"M 0,{-self.color_rect_size / 2} L {layout.bar_width},{-self.color_rect_size / 2} "
            f"L {layout.bar_width},{self.color_rect_size / 2} L 0,{self.color_rect_size / 2} z"
        )

        if layout.compact:
            for entry in layout.compact_entries:
                properties = entry.properties
                gradient_id = self._gradient_id(entry.key, properties)
                entry_group = Group(debug=False)
                entry_group.attribs["data-legend-key"] = str(entry.key)

                gradient = LinearGradient(start=(0, 0), end=("100%", 0), id=gradient_id)
                gradient.add_stop_color(offset="0%", color=properties["min_color"])
                gradient.add_stop_color(offset="100%", color=properties["max_color"])
                entry_group.add(gradient)

                label_path = generate_text_path(
                    entry.key,
                    0,
                    0,
                    0,
                    self.font_size,
                    "normal",
                    font,
                    dominant_baseline="central",
                    text_anchor="start",
                )
                label_path.translate(0, entry.label_y)
                entry_group.add(label_path)

                grad_rect = Path(
                    d=path_desc,
                    fill=f"url(#{gradient_id})",
                    stroke=properties["stroke"],
                    stroke_width=properties["width"],
                )
                grad_rect.translate(layout.bar_x, entry.bar_y)
                entry_group.add(grad_rect)
                group.add(entry_group)

            min_label = generate_text_path(
                layout.min_label_text,
                0,
                0,
                0,
                self.font_size,
                "normal",
                font,
                dominant_baseline="hanging",
                text_anchor="start",
            )
            min_label.translate(layout.bar_x, layout.scale_y)
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
            max_label.translate(layout.bar_x + layout.bar_width, layout.scale_y)
            group.add(max_label)
            return group

        for entry in layout.single_entries:
            properties = entry.properties
            gradient_id = self._gradient_id(entry.key, properties)
            entry_group = Group(debug=False)
            entry_group.attribs["data-legend-key"] = str(entry.key)

            gradient = LinearGradient(start=(0, 0), end=("100%", 0), id=gradient_id)
            gradient.add_stop_color(offset="0%", color=properties["min_color"])
            gradient.add_stop_color(offset="100%", color=properties["max_color"])
            entry_group.add(gradient)

            title_path = generate_text_path(
                entry.key,
                0,
                0,
                0,
                self.font_size,
                "normal",
                font,
                dominant_baseline="hanging",
                text_anchor="middle",
            )
            title_path.translate(entry.title_x, entry.title_y)
            entry_group.add(title_path)

            grad_rect = Path(
                d=path_desc,
                fill=f"url(#{gradient_id})",
                stroke=properties["stroke"],
                stroke_width=properties["width"],
            )
            grad_rect.translate(entry.bar_x, entry.bar_y)
            entry_group.add(grad_rect)

            min_label = generate_text_path(
                layout.min_label_text,
                0,
                0,
                0,
                self.font_size,
                "normal",
                font,
                dominant_baseline="hanging",
                text_anchor="start",
            )
            min_label.translate(entry.min_label_x, entry.scale_label_y)
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
            max_label.translate(entry.max_label_x, entry.scale_label_y)
            entry_group.add(max_label)
            group.add(entry_group)

        return group

    def add_elements_to_group(self):
        path_desc = (
            f"M {0},{-0.5 * self.color_rect_size} "
            f"L {self.layout.width},{-0.5 * self.color_rect_size} "
            f"L {self.layout.width},{self.layout.height -0.5 * self.color_rect_size} "
            f"L {0},{self.layout.height -0.5 * self.color_rect_size} z"
        )
        rect_path = Path(d=path_desc, fill="none", stroke="none", stroke_width=0)
        self.legend_group.add(rect_path)
        path_desc = self.create_rectangle_path_for_legend()
        font = self.font_family
        feature_group = Group(id="feature_legend") if self.layout.gradient is not None else self.legend_group
        for entry in self.layout.solid_entries:
            entry_group = Group(debug=False)
            entry_group.attribs["data-legend-key"] = str(entry.key)

            rect_path = Path(
                d=path_desc,
                fill=entry.properties["fill"],
                stroke=entry.properties["stroke"],
                stroke_width=entry.properties["width"],
            )
            rect_path.translate(entry.rect_x, entry.rect_y)
            entry_group.add(rect_path)
            legend_path = generate_text_path(
                entry.key,
                0,
                0,
                0,
                self.font_size,
                "normal",
                font,
                dominant_baseline="central",
                text_anchor="start",
            )
            legend_path.translate(entry.text_x, entry.text_y)
            entry_group.add(legend_path)
            feature_group.add(entry_group)

        if self.layout.gradient is not None and self.layout.solid_entries:
            self.legend_group.add(feature_group)
        if self.layout.gradient is not None:
            gradient_group = self._build_gradient_legend(self.layout.gradient)
            gradient_group.translate(self.layout.gradient_x, self.layout.gradient_y)
            self.legend_group.add(gradient_group)
        return self.legend_group

    def get_group(self) -> Group:
        return self.legend_group


__all__ = ["LegendGroup"]
