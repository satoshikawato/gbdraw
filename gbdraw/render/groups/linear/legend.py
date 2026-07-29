#!/usr/bin/env python
# coding: utf-8

from svgwrite.container import Group
from svgwrite.gradients import LinearGradient
from svgwrite.path import Path

from ....config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from ....configurators.legend import LegendMeasurement
from ....legend.linear_layout import (
    LinearFeatureLegendLayout,
    LinearGradientLegendLayout,
    LinearOrientationLegendLayout,
)
from ....svg.ids import stable_svg_id
from ....svg.text_path import generate_text_path


class LegendGroup:
    """Render premeasured horizontal and vertical linear legend layouts."""

    def __init__(
        self,
        canvas_config,
        legend_measurement: LegendMeasurement,
        legend_table,
        *,
        cfg: GbdrawConfig,
    ):
        self.legend_group = Group(id="legend")
        self._cfg = cfg
        self.canvas_config = canvas_config
        self.legend_measurement = legend_measurement
        self.legend_table = legend_table
        self.font_family: str = legend_measurement.font_family
        self.font_weight: str = legend_measurement.font_weight
        self.font_size: float = legend_measurement.font_size
        self.rect_size: float = legend_measurement.color_rect_size
        self.legend_position = self.canvas_config.legend_position
        self.dpi: int = cfg.canvas.dpi
        self.layout = legend_measurement.linear_layout
        if self.layout is None:
            if legend_table:
                raise ValueError(
                    "Linear legend rendering requires a measured linear layout."
                )
        else:
            self._build_dual_legends()

    @property
    def horizontal_legend_width(self) -> float:
        return self.layout.horizontal.width if self.layout is not None else 0.0

    @property
    def horizontal_legend_height(self) -> float:
        return self.layout.horizontal.height if self.layout is not None else 0.0

    @property
    def vertical_legend_width(self) -> float:
        return self.layout.vertical.width if self.layout is not None else 0.0

    @property
    def vertical_legend_height(self) -> float:
        return self.layout.vertical.height if self.layout is not None else 0.0

    @property
    def legend_width(self) -> float:
        return self.layout.active.width if self.layout is not None else 0.0

    @property
    def legend_height(self) -> float:
        return self.layout.active.height if self.layout is not None else 0.0

    def create_rectangle_path_for_legend(self) -> str:
        start_y_top = -self.rect_size / 2
        start_y_bottom = self.rect_size / 2
        return (
            f"M 0,{start_y_top} L {self.rect_size},{start_y_top} "
            f"L {self.rect_size},{start_y_bottom} L 0,{start_y_bottom} z"
        )

    @staticmethod
    def _pairwise_group(orientation: str) -> Group:
        group = Group(id=f"pairwise_legend_{orientation}", debug=False)
        group.attribs["data-gbdraw-role"] = "comparison-legend"
        group.attribs["data-gbdraw-orientation"] = orientation
        return group

    @staticmethod
    def _build_pairwise_gradient_id(
        key: str,
        properties: dict,
        orientation: str,
    ) -> str:
        return stable_svg_id(
            "blast_legend_grad",
            "linear-pairwise-gradient",
            key,
            properties["min_color"],
            properties["max_color"],
            namespace=orientation,
        )

    def _build_feature_legend(
        self,
        layout: LinearFeatureLegendLayout,
        *,
        group_id: str,
    ) -> Group:
        group = Group(id=group_id)
        path_desc = self.create_rectangle_path_for_legend()
        for entry in layout.entries:
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
                self.font_family,
                dominant_baseline="central",
                text_anchor="start",
            )
            legend_path.translate(entry.text_x, entry.text_y)
            entry_group.add(legend_path)
            group.add(entry_group)
        return group

    def _add_gradient_scale_labels(
        self,
        group: Group,
        layout: LinearGradientLegendLayout,
    ) -> None:
        label_0 = generate_text_path(
            layout.min_label_text,
            0,
            0,
            0,
            self.font_size,
            "normal",
            self.font_family,
            dominant_baseline="hanging",
            text_anchor="start",
        )
        label_0.translate(layout.min_label_x, layout.scale_label_y)
        group.add(label_0)

        label_100 = generate_text_path(
            "100%",
            0,
            0,
            0,
            self.font_size,
            "normal",
            self.font_family,
            dominant_baseline="hanging",
            text_anchor="end",
        )
        label_100.translate(layout.max_label_x, layout.scale_label_y)
        group.add(label_100)

    def _build_pairwise_legend(
        self,
        layout: LinearGradientLegendLayout,
        orientation: str,
    ) -> Group:
        group = self._pairwise_group(orientation)
        path_desc = (
            f"M 0,{-self.rect_size / 2} "
            f"L {layout.bar_width},{-self.rect_size / 2} "
            f"L {layout.bar_width},{self.rect_size / 2} "
            f"L 0,{self.rect_size / 2} z"
        )
        for entry in layout.entries:
            properties = entry.properties
            gradient_id = self._build_pairwise_gradient_id(
                entry.key,
                properties,
                orientation,
            )
            entry_group = Group(debug=False)
            entry_group.attribs["data-legend-key"] = str(entry.key)

            gradient = LinearGradient(
                start=(0, 0),
                end=("100%", 0),
                id=gradient_id,
            )
            gradient.add_stop_color(
                offset="0%",
                color=properties["min_color"],
            )
            gradient.add_stop_color(
                offset="100%",
                color=properties["max_color"],
            )
            entry_group.add(gradient)

            title_path = generate_text_path(
                entry.key,
                0,
                0,
                0,
                self.font_size,
                "normal",
                self.font_family,
                dominant_baseline=(
                    "central" if layout.compact else "hanging"
                ),
                text_anchor="start" if layout.compact else "middle",
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
            if not layout.compact:
                self._add_gradient_scale_labels(entry_group, layout)
            group.add(entry_group)

        if layout.compact:
            self._add_gradient_scale_labels(group, layout)
        return group

    def _build_orientation_group(
        self,
        *,
        group_id: str,
        orientation: str,
        layout: LinearOrientationLegendLayout,
        hidden: bool,
    ) -> Group:
        group = Group(id=group_id)
        feature_group = self._build_feature_legend(
            layout.feature,
            group_id=(
                "feature_legend_h"
                if orientation == "h"
                else "feature_legend_v"
            ),
        )
        if layout.feature_x or layout.feature_y:
            feature_group.translate(layout.feature_x, layout.feature_y)
        group.add(feature_group)

        if layout.gradient is not None:
            pairwise_group = self._build_pairwise_legend(
                layout.gradient,
                orientation,
            )
            if orientation == "h":
                pairwise_group.translate(layout.gradient_x, 0)
                if layout.gradient_y:
                    pairwise_group.translate(0, layout.gradient_y)
            else:
                pairwise_group.translate(0, layout.gradient_y)
            group.add(pairwise_group)
        if hidden:
            group.attribs["display"] = "none"
        return group

    def _build_dual_legends(self) -> None:
        if self.layout is None:
            return
        is_horizontal = self.layout.active_orientation == "horizontal"
        horizontal_group = self._build_orientation_group(
            group_id="legend_horizontal",
            orientation="h",
            layout=self.layout.horizontal,
            hidden=not is_horizontal,
        )
        vertical_group = self._build_orientation_group(
            group_id="legend_vertical",
            orientation="v",
            layout=self.layout.vertical,
            hidden=is_horizontal,
        )
        self.legend_group.add(horizontal_group)
        self.legend_group.add(vertical_group)

    def get_group(self) -> Group:
        return self.legend_group

    def get_horizontal_dimensions(self) -> tuple[float, float]:
        return self.horizontal_legend_width, self.horizontal_legend_height

    def get_vertical_dimensions(self) -> tuple[float, float]:
        return self.vertical_legend_width, self.vertical_legend_height


__all__ = ["LegendGroup"]
