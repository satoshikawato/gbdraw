#!/usr/bin/env python
# coding: utf-8

from gbdraw.svg.elements import Group, LinearGradient, Path

from ....configurators.legend import LegendMeasurement
from ....legend.circular_layout import (
    RING_BAND_LABEL,
    RING_INNER_CAPTION,
    RING_OUTER_CAPTION,
    CircularGradientLegendLayout,
    CircularRingLegendLayout,
)
from ....legend.table import as_legend_table, LegendTable, write_legend_row_names
from ....core.text import italic_markup_runs
from ....legend.metrics import write_gradient_scale_names
from ....render.group_content import GroupContent
from ....svg.text_path import generate_runs_text_path
from ....svg.ids import stable_svg_id
from ....svg.result_records import LEGEND_ORDER_RECORD
from ....svg.path_data import path_number, path_point

# The colour of the similarity-ring block's captions (``outer``, ``inner``, ``identity``) and scale labels; the
# ring names keep the Legend's text colour. Drawn only here; the Web's Legend layout reads positions, not colours.
RING_LEGEND_MUTED_TEXT_COLOR = "#5a6273"


class LegendGroup:
    def __init__(
        self,
        legend_measurement: LegendMeasurement,
        legend_table: LegendTable,
    ) -> None:
        self.legend_group = Group(id="legend")
        self.legend_measurement = legend_measurement
        self.legend_table = as_legend_table(legend_table)
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

    def _text(
        self,
        text: str,
        x: float,
        y: float,
        *,
        baseline: str = "central",
        anchor: str = "start",
        fill: str | None = None,
    ):
        # A caption's ``<i>`` runs draw in italics; the fixed texts are one upright run.
        path = generate_runs_text_path(
            italic_markup_runs(text), self.font_size, self.font_family, dominant_baseline=baseline, text_anchor=anchor
        )
        if fill is not None:
            path.attribs["fill"] = fill
        path.translate(x, y)
        return path

    def _add_ring_legend(self, group: Group, layout: CircularRingLegendLayout) -> None:
        """The similarity-ring block (PD-OI-110): ``outer``, one row per drawn ring (``data-source-index``,
        ``data-track-index``), ``inner``, and the identity band; the band and captions are not Legend rows."""

        group.add(self._text(RING_OUTER_CAPTION, 0, layout.outer_caption_y, fill=RING_LEGEND_MUTED_TEXT_COLOR))
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
        group.add(self._text(RING_INNER_CAPTION, 0, layout.inner_caption_y, fill=RING_LEGEND_MUTED_TEXT_COLOR))

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
        band.add(self._text(RING_BAND_LABEL, layout.band_label_x, layout.band_y, fill=RING_LEGEND_MUTED_TEXT_COLOR))
        scale = {"baseline": "hanging", "fill": RING_LEGEND_MUTED_TEXT_COLOR}
        band.add(self._text(layout.min_label_text, layout.min_label_x, layout.scale_label_y, **scale))
        band.add(self._text(layout.max_label_text, layout.max_label_x, layout.scale_label_y, anchor="end", **scale))
        group.add(band)

    def _add_gradient_rows(self, group: Group, layout: CircularGradientLegendLayout) -> None:
        """The gradient rows: each row's caption and bar, then one scale line per group of rows that share scale
        texts (a compact block), or one row's centred title, its bar and its scale line."""

        half = self.color_rect_size / 2
        bar_path = (
            f"M 0,{path_number(-half)} L {path_point(layout.bar_width, -half)} "
            f"L {path_point(layout.bar_width, half)} L 0,{path_number(half)} z"
        )
        rows = [
            (entry.key, entry.properties, entry.label_x, entry.label_y, layout.bar_x, entry.bar_y)
            for entry in layout.compact_entries
        ] + [
            (entry.key, entry.properties, entry.title_x, entry.title_y, entry.bar_x, entry.bar_y)
            for entry in layout.single_entries
        ]
        caption = {"baseline": "central", "anchor": "start"} if layout.compact else {
            "baseline": "hanging",
            "anchor": "middle",
        }
        for key, properties, caption_x, caption_y, bar_x, bar_y in rows:
            gradient_id = stable_svg_id(
                "circular_legend_grad", "circular-gradient-row", key, properties["min_color"], properties["max_color"]
            )
            entry_group = Group()
            drawn_name = self._name_entry(entry_group.attribs, key)
            write_gradient_scale_names(entry_group.attribs, properties)
            gradient = LinearGradient(start=(0, 0), end=("100%", 0), id=gradient_id)
            gradient.add_stop_color(offset="0%", color=properties["min_color"])
            gradient.add_stop_color(offset="100%", color=properties["max_color"])
            entry_group.add(gradient)
            entry_group.add(self._text(drawn_name, caption_x, caption_y, **caption))
            bar = Path(
                d=bar_path,
                fill=f"url(#{gradient_id})",
                stroke=properties["stroke"],
                stroke_width=properties["width"],
            )
            bar.translate(bar_x, bar_y)
            entry_group.add(bar)
            group.add(entry_group)
        for scale in layout.scales:
            group.add(self._text(scale.min_label_text, scale.min_label_x, scale.scale_label_y, baseline="hanging"))
            group.add(
                self._text(
                    scale.max_label_text, scale.max_label_x, scale.scale_label_y, baseline="hanging", anchor="end"
                )
            )

    def _comparison_legend_group(self) -> Group:
        group = Group(id="conservation_identity_legend")
        group.attribs["data-gbdraw-role"] = "comparison-legend"
        group.attribs["data-gbdraw-orientation"] = "circular"
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
        # The side block: the ring block and the gradient rows below it.
        has_side_block = layout.has_gradient
        feature_group = Group(id="feature_legend") if has_side_block else self.legend_group
        order_record = self.legend_table.order_record()
        if order_record is not None:
            feature_group.attribs[LEGEND_ORDER_RECORD] = order_record
        # The shown rows where the layout placed them; a hidden row (a Legend row edit) at a shown row's place.
        feature_rows = self.legend_table.feature_rows({entry.key: entry for entry in layout.solid_entries})
        for key, properties, entry in feature_rows:
            entry_group = Group()
            drawn_name = self._name_entry(entry_group.attribs, key)

            rect_path = Path(
                d=path_desc,
                fill=properties["fill"],
                stroke=properties["stroke"],
                stroke_width=properties["width"],
            )
            legend_path = generate_runs_text_path(
                italic_markup_runs(drawn_name), self.font_size, font, dominant_baseline="central", text_anchor="start"
            )
            if entry is not None:
                rect_path.translate(entry.rect_x, entry.rect_y)
                legend_path.translate(entry.text_x, entry.text_y)
            entry_group.add(rect_path)
            entry_group.add(legend_path)
            feature_group.add(entry_group)

        if has_side_block and feature_rows:
            self.legend_group.add(feature_group)
        if has_side_block:
            # One side block: the ring block, then the gradient rows below it.
            side_group = self._comparison_legend_group()
            if layout.gradient is not None:
                self._add_ring_legend(side_group, layout.gradient)
            if layout.gradient_rows is not None:
                self._add_gradient_rows(side_group, layout.gradient_rows)
            side_group.translate(layout.gradient_x, layout.gradient_y)
            self.legend_group.add(side_group)
        return self.legend_group

    def get_group(self) -> Group:
        return self.legend_group


def circular_legend_content(legend_measurement: LegendMeasurement, legend_table: LegendTable) -> GroupContent:
    """The ``legend_content`` stage of a Circular canvas outside the record build."""

    return GroupContent.of(LegendGroup(legend_measurement, legend_table).get_group())


__all__ = ["LegendGroup", "circular_legend_content"]
