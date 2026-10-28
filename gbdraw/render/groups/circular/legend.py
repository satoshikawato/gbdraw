#!/usr/bin/env python
# coding: utf-8

from typing import Mapping

from gbdraw.svg.elements import Group, LinearGradient, Path

from ....configurators.legend import LegendMeasurement
from ....legend.circular_layout import (
    RING_BAND_LABEL,
    RING_INNER_CAPTION,
    RING_OUTER_CAPTION,
    CircularGradientLegendLayout,
    CircularRingLegendLayout,
)
from ....legend.table import LegendRowKey, LegendTable, write_legend_row_names
from ....render.group_content import GroupContent
from ....svg.text_path import generate_text_path
from ....svg.ids import stable_svg_id
from ....svg.path_data import path_number, path_point


class LegendGroup:
    def __init__(
        self,
        legend_measurement: LegendMeasurement,
        legend_table: LegendTable,
    ) -> None:
        self.legend_group = Group(id="legend")
        self.legend_measurement = legend_measurement
        self.legend_table = legend_table
        self._name_entry = write_legend_row_names(legend_table)
        self.font_family = self.legend_measurement.font_family
        self.font_size = self.legend_measurement.font_size

        self.color_rect_size = self.legend_measurement.color_rect_size
        self.num_of_lines = len(self.legend_table.keys())
        self.layout = self.legend_measurement.circular_layout
        if self.layout is None:
            if self.legend_table:
                raise ValueError(
                    "Circular legend rendering requires a measured circular layout."
                )
        else:
            self.add_elements_to_group()

    def create_rectangle_path_for_legend(self) -> str:
        normalized_start: float = 0
        normalized_end: float = self.color_rect_size
        start_y_top, start_y_bottom = -0.5 * self.color_rect_size, 0.5 * self.color_rect_size
        end_y_top, end_y_bottom = -0.5 * self.color_rect_size, 0.5 * self.color_rect_size
        rectangle_path: str = (
            f"M {path_point(normalized_start, start_y_top)} L {path_point(normalized_end, end_y_top)} "
            f"L {path_point(normalized_end, end_y_bottom)} L {path_point(normalized_start, start_y_bottom)} z"
        )
        return rectangle_path

    def _gradient_id(self, key: LegendRowKey, properties: Mapping[str, object]) -> str:
        return stable_svg_id(
            "circular_legend_grad",
            "circular-conservation-gradient",
            key,
            properties["min_color"],
            properties["max_color"],
        )

    def _text(self, text: str, x: float, y: float, *, baseline: str = "central", anchor: str = "start"):
        path = generate_text_path(
            text, 0, 0, 0, self.font_size, "normal", self.font_family, dominant_baseline=baseline, text_anchor=anchor
        )
        path.translate(x, y)
        return path

    def _build_ring_legend(self, layout: CircularRingLegendLayout) -> Group:
        """The similarity-ring block (PD-OI-110): ``outer``, one row per drawn ring (``data-source-index``,
        ``data-track-index``), ``inner``, and the identity band; the band and captions are not Legend rows."""

        group = self._comparison_legend_group()
        group.add(self._text(RING_OUTER_CAPTION, 0, layout.outer_caption_y))
        swatch_path = self.create_rectangle_path_for_legend()
        for row in layout.rows:
            row_group = Group()
            drawn_name = self._name_entry(row_group.attribs, row.key)
            row_group.attribs["data-source-index"] = str(row.properties["source_index"])
            row_group.attribs["data-track-index"] = str(row.properties["track_index"])
            if layout.show_swatches:
                swatch = Path(d=swatch_path, fill=row.properties["fill"], stroke="none", stroke_width=0)
                swatch.translate(row.swatch_x, row.y)
                row_group.add(swatch)
            row_group.add(self._text(drawn_name, row.text_x, row.y))
            group.add(row_group)
        group.add(self._text(RING_INNER_CAPTION, 0, layout.inner_caption_y))

        band = Group()
        band.attribs["data-gbdraw-role"] = "conservation-identity-band"
        properties = layout.rows[0].properties
        gradient_id = stable_svg_id(
            "circular_legend_grad",
            "circular-conservation-gradient",
            RING_BAND_LABEL,
            properties["min_color"],
            properties["max_color"],
        )
        gradient = LinearGradient(start=(0, 0), end=("100%", 0), id=gradient_id)
        gradient.add_stop_color(offset="0%", color=properties["min_color"])
        gradient.add_stop_color(offset="100%", color=properties["max_color"])
        band.add(gradient)
        half = self.color_rect_size / 2
        bar = Path(
            d=(
                f"M 0,{path_number(-half)} L {path_point(layout.bar_width, -half)} "
                f"L {path_point(layout.bar_width, half)} L 0,{path_number(half)} z"
            ),
            fill=f"url(#{gradient_id})",
            stroke="none",
            stroke_width=0,
        )
        bar.translate(0, layout.band_y)
        band.add(bar)
        band.add(self._text(RING_BAND_LABEL, layout.band_label_x, layout.band_y))
        band.add(self._text(layout.min_label_text, layout.min_label_x, layout.scale_label_y, baseline="hanging"))
        band.add(self._text("100%", layout.max_label_x, layout.scale_label_y, baseline="hanging", anchor="end"))
        group.add(band)
        return group

    def _comparison_legend_group(self) -> Group:
        group = Group(id="conservation_identity_legend")
        group.attribs["data-gbdraw-role"] = "comparison-legend"
        group.attribs["data-gbdraw-orientation"] = "circular"
        return group

    def _build_gradient_legend(self, layout: CircularGradientLegendLayout) -> Group:
        group = self._comparison_legend_group()
        font = self.font_family
        path_desc = (
            f"M 0,{path_number(-self.color_rect_size / 2)} L {path_point(layout.bar_width, -self.color_rect_size / 2)} "
            f"L {path_point(layout.bar_width, self.color_rect_size / 2)} L 0,{path_number(self.color_rect_size / 2)} z"
        )

        if layout.compact:
            for compact_entry in layout.compact_entries:
                properties = compact_entry.properties
                gradient_id = self._gradient_id(compact_entry.key, properties)
                entry_group = Group()
                drawn_name = self._name_entry(entry_group.attribs, compact_entry.key)

                gradient = LinearGradient(start=(0, 0), end=("100%", 0), id=gradient_id)
                gradient.add_stop_color(offset="0%", color=properties["min_color"])
                gradient.add_stop_color(offset="100%", color=properties["max_color"])
                entry_group.add(gradient)

                label_path = generate_text_path(
                    drawn_name,
                    0,
                    0,
                    0,
                    self.font_size,
                    "normal",
                    font,
                    dominant_baseline="central",
                    text_anchor="start",
                )
                label_path.translate(compact_entry.label_x, compact_entry.label_y)
                entry_group.add(label_path)

                grad_rect = Path(
                    d=path_desc,
                    fill=f"url(#{gradient_id})",
                    stroke=properties["stroke"],
                    stroke_width=properties["width"],
                )
                grad_rect.translate(layout.bar_x, compact_entry.bar_y)
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
            min_label.translate(layout.min_label_x, layout.scale_label_y)
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
            max_label.translate(layout.max_label_x, layout.scale_label_y)
            group.add(max_label)
            return group

        for entry in layout.single_entries:
            properties = entry.properties
            gradient_id = self._gradient_id(entry.key, properties)
            entry_group = Group()
            drawn_name = self._name_entry(entry_group.attribs, entry.key)

            gradient = LinearGradient(start=(0, 0), end=("100%", 0), id=gradient_id)
            gradient.add_stop_color(offset="0%", color=properties["min_color"])
            gradient.add_stop_color(offset="100%", color=properties["max_color"])
            entry_group.add(gradient)

            title_path = generate_text_path(
                drawn_name,
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
            min_label.translate(layout.min_label_x, layout.scale_label_y)
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
            max_label.translate(layout.max_label_x, layout.scale_label_y)
            entry_group.add(max_label)
            group.add(entry_group)

        return group

    def add_elements_to_group(self):
        layout = self.layout
        assert layout is not None
        path_desc = (
            f"M {path_point(0, -0.5 * self.color_rect_size)} "
            f"L {path_point(layout.width, -0.5 * self.color_rect_size)} "
            f"L {path_point(layout.width, layout.height -0.5 * self.color_rect_size)} "
            f"L {path_point(0, layout.height -0.5 * self.color_rect_size)} z"
        )
        rect_path = Path(d=path_desc, fill="none", stroke="none", stroke_width=0)
        self.legend_group.add(rect_path)
        path_desc = self.create_rectangle_path_for_legend()
        font = self.font_family
        feature_group = Group(id="feature_legend") if layout.gradient is not None else self.legend_group
        for entry in layout.solid_entries:
            entry_group = Group()
            drawn_name = self._name_entry(entry_group.attribs, entry.key)

            rect_path = Path(
                d=path_desc,
                fill=entry.properties["fill"],
                stroke=entry.properties["stroke"],
                stroke_width=entry.properties["width"],
            )
            rect_path.translate(entry.rect_x, entry.rect_y)
            entry_group.add(rect_path)
            legend_path = generate_text_path(
                drawn_name,
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

        if layout.gradient is not None and layout.solid_entries:
            self.legend_group.add(feature_group)
        if layout.gradient is not None:
            gradient_group = (
                self._build_ring_legend(layout.gradient)
                if isinstance(layout.gradient, CircularRingLegendLayout)
                else self._build_gradient_legend(layout.gradient)
            )
            gradient_group.translate(layout.gradient_x, layout.gradient_y)
            self.legend_group.add(gradient_group)
        return self.legend_group

    def get_group(self) -> Group:
        return self.legend_group


def circular_legend_content(legend_measurement: LegendMeasurement, legend_table: LegendTable) -> GroupContent:
    """The ``legend_content`` stage of a Circular canvas outside the record build."""

    return GroupContent.of(LegendGroup(legend_measurement, legend_table).get_group())


__all__ = ["LegendGroup", "circular_legend_content"]
