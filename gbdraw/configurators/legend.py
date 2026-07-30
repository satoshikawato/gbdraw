#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from gbdraw.config.models import (  # type: ignore[reportMissingImports]
    CircularRenderProfile,
    LinearRenderProfile,
    RenderProfile,
)
from gbdraw.legend.circular_layout import (  # type: ignore[reportMissingImports]
    CircularLegendLayout,
    build_circular_legend_layout,
)
from gbdraw.legend.linear_layout import (  # type: ignore[reportMissingImports]
    LinearLegendLayout,
    build_linear_legend_layout,
)

if TYPE_CHECKING:
    from gbdraw.canvas import CircularCanvasConfigurator, LinearCanvasConfigurator


@dataclass(frozen=True, slots=True)
class LegendMeasurement:
    """Immutable legend geometry and rendering inputs for one canvas."""

    font_family: str
    font_weight: str
    font_size: float
    color_rect_size: float
    dpi: int
    legend_width: float
    legend_height: float
    total_feature_legend_width: float
    pairwise_legend_width: float
    num_of_lines: int
    num_of_columns: int
    num_of_items_per_line: int
    has_gradient: bool
    circular_layout: CircularLegendLayout | None = None
    linear_layout: LinearLegendLayout | None = None


class LegendDrawingConfigurator:
    def __init__(
        self,
        color_table,
        default_colors,
        selected_features_set,
        profile: RenderProfile,
        gc_config,
        skew_config,
        feature_config,
        legend_table=None,
        blast_config=None,
        canvas_config=None,
    ) -> None:
        cfg = profile.config
        self.color_table = color_table
        self.default_colors = default_colors
        self.selected_features_set = selected_features_set
        self.canvas_config = canvas_config
        self._profile = profile
        self.gc_config = gc_config
        self.skew_config = skew_config
        self.feature_config = feature_config
        self.blast_config = blast_config
        self.length_param = self.canvas_config.length_param
        self.font_size: float = cfg.objects.legends.font_size.for_length_param(self.length_param)
        self.font_weight: str = cfg.objects.legends.font_weight
        self.font_family: str = cfg.objects.text.font_family
        self.text_anchor: str = cfg.objects.legends.text_anchor
        self.color_rect_size: float = cfg.objects.legends.color_rect_size.for_length_param(self.length_param)
        self.dominant_baseline: str = cfg.objects.legends.dominant_baseline
        self.dpi: int = cfg.canvas.dpi

    def _measurement(
        self,
        *,
        legend_width: float,
        legend_height: float,
        total_feature_legend_width: float,
        pairwise_legend_width: float,
        num_of_lines: int,
        num_of_columns: int,
        num_of_items_per_line: int,
        has_gradient: bool,
        circular_layout: CircularLegendLayout | None = None,
        linear_layout: LinearLegendLayout | None = None,
    ) -> LegendMeasurement:
        return LegendMeasurement(
            font_family=self.font_family,
            font_weight=self.font_weight,
            font_size=float(self.font_size),
            color_rect_size=float(self.color_rect_size),
            dpi=int(self.dpi),
            legend_width=float(legend_width),
            legend_height=float(legend_height),
            total_feature_legend_width=float(total_feature_legend_width),
            pairwise_legend_width=float(pairwise_legend_width),
            num_of_lines=int(num_of_lines),
            num_of_columns=int(num_of_columns),
            num_of_items_per_line=int(num_of_items_per_line),
            has_gradient=bool(has_gradient),
            circular_layout=circular_layout,
            linear_layout=linear_layout,
        )

    def _empty_measurement(self) -> LegendMeasurement:
        return self._measurement(
            legend_width=0,
            legend_height=0,
            total_feature_legend_width=0,
            pairwise_legend_width=0,
            num_of_lines=0,
            num_of_columns=0,
            num_of_items_per_line=0,
            has_gradient=False,
        )

    def measure_legend(
        self,
        legend_table,
        canvas_config: CircularCanvasConfigurator | LinearCanvasConfigurator,
    ) -> LegendMeasurement:
        if not legend_table:
            return self._empty_measurement()
        if isinstance(self._profile, CircularRenderProfile):
            layout = build_circular_legend_layout(
                legend_table,
                legend_position=canvas_config.legend_position,
                canvas_width=float(canvas_config.total_width),
                font_family=self.font_family,
                font_size=float(self.font_size),
                dpi=int(self.dpi),
                color_rect_size=float(self.color_rect_size),
            )
            return self._measurement(
                legend_width=layout.width,
                legend_height=layout.height,
                total_feature_legend_width=layout.feature_width,
                pairwise_legend_width=layout.pairwise_legend_width,
                num_of_lines=layout.num_lines,
                num_of_columns=layout.num_columns,
                num_of_items_per_line=layout.num_items_per_line,
                has_gradient=layout.has_gradient,
                circular_layout=layout,
            )
        if not isinstance(self._profile, LinearRenderProfile):
            raise TypeError(
                "Legend measurement requires a CircularRenderProfile or LinearRenderProfile."
            )
        layout = build_linear_legend_layout(
            legend_table,
            legend_position=canvas_config.legend_position,
            canvas_width=float(canvas_config.total_width),
            font_family=self.font_family,
            font_size=float(self.font_size),
            dpi=int(self.dpi),
            color_rect_size=float(self.color_rect_size),
        )
        active = layout.active
        solid_count = len(active.feature.entries)
        num_lines = active.feature.num_lines
        num_items_per_line = (
            max(1, (solid_count + max(1, num_lines) - 1) // max(1, num_lines))
            if solid_count
            else 0
        )
        return self._measurement(
            legend_width=active.width,
            legend_height=active.height,
            total_feature_legend_width=active.feature.width,
            pairwise_legend_width=(
                active.gradient.width
                if active.gradient is not None
                else 0.0
            ),
            num_of_lines=num_lines,
            num_of_columns=(
                len(legend_table)
                if layout.active_orientation == "horizontal"
                else 1
            ),
            num_of_items_per_line=num_items_per_line,
            has_gradient=layout.has_gradient,
            linear_layout=layout,
        )


__all__ = ["LegendDrawingConfigurator", "LegendMeasurement"]
