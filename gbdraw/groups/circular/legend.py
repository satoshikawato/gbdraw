#!/usr/bin/env python
# coding: utf-8

from svgwrite.container import Group
from svgwrite.path import Path

from ...svg.text_path import generate_text_path


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

    def add_elements_to_group(self):
        count = 0
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
        for key, properties in self.legend_table.items():
            if properties["type"] == "solid":
                rect_path = Path(
                    d=path_desc,
                    fill=properties["fill"],
                    stroke=properties["stroke"],
                    stroke_width=properties["width"],
                )
                rect_path.translate(0, count * line_margin)
                self.legend_group.add(rect_path)
                legend_path = generate_text_path(
                    key, 0, 0, 0, self.font_size, "normal", font, dominant_baseline="central", text_anchor="start"
                )
                legend_path.translate(x_margin, count * line_margin)
                self.legend_group.add(legend_path)
                count += 1
        return self.legend_group

    def get_group(self) -> Group:
        return self.legend_group


__all__ = ["LegendGroup"]


