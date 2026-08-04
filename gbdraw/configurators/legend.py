#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Mapping

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
from gbdraw.legend.metrics import (  # type: ignore[reportMissingImports]
    legend_line_height,
    legend_text_x_offset,
)
from gbdraw.layout.spatial import Aabb, union_aabbs  # type: ignore[reportMissingImports]


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
    local_bounds: Aabb = Aabb(0.0, 0.0, 0.0, 0.0)

    @property
    def reflow_metrics(self) -> dict[str, float]:
        """Return Python-owned constraints for browser legend editing."""

        return {
            "colorRectSize": float(self.color_rect_size),
            "lineHeight": legend_line_height(self.color_rect_size),
            "textXOffset": legend_text_x_offset(self.color_rect_size),
        }


def _legend_half_stroke_width(
    legend_table: Mapping[object, Mapping[str, object]],
) -> float:
    max_width = 0.0
    for properties in legend_table.values():
        if str(properties.get("stroke", "none")).strip().lower() == "none":
            continue
        raw_width = properties.get("width", 0.0)
        try:
            width = float(raw_width or 0.0)
        except (TypeError, ValueError) as exc:
            raise ValueError("legend stroke widths must be finite non-negative numbers") from exc
        if not math.isfinite(width) or width < 0.0:
            raise ValueError("legend stroke widths must be finite non-negative numbers")
        max_width = max(max_width, width)
    return 0.5 * max_width


def _entry_half_stroke_width(properties: Mapping[str, object]) -> float:
    return _legend_half_stroke_width({"entry": properties})


def _circular_legend_local_bounds(
    layout: CircularLegendLayout,
    *,
    color_rect_size: float,
) -> Aabb:
    rect_size = float(color_rect_size)
    bounds: list[Aabb] = [
        Aabb(
            0.0,
            -0.5 * rect_size,
            float(layout.width),
            float(layout.height) - (0.5 * rect_size),
        )
    ]
    for entry in layout.solid_entries:
        bounds.append(
            Aabb(
                float(entry.rect_x),
                float(entry.rect_y) - (0.5 * rect_size),
                float(entry.rect_x) + rect_size,
                float(entry.rect_y) + (0.5 * rect_size),
            ).expanded(_entry_half_stroke_width(entry.properties))
        )
    gradient = layout.gradient
    if gradient is not None:
        entries = (*gradient.compact_entries, *gradient.single_entries)
        for entry in entries:
            bounds.append(
                Aabb(
                    float(layout.gradient_x) + float(gradient.bar_x),
                    float(layout.gradient_y) + float(entry.bar_y) - (0.5 * rect_size),
                    float(layout.gradient_x) + float(gradient.bar_x) + float(gradient.bar_width),
                    float(layout.gradient_y) + float(entry.bar_y) + (0.5 * rect_size),
                ).expanded(_entry_half_stroke_width(entry.properties))
            )
    combined = union_aabbs(bounds)
    if combined is None:  # pragma: no cover - the layout rectangle is always present
        raise RuntimeError("circular legend bounds unexpectedly empty")
    return combined


def _linear_legend_local_bounds(
    layout: LinearLegendLayout,
    *,
    color_rect_size: float,
) -> Aabb:
    active = layout.active
    rect_size = float(color_rect_size)
    bounds: list[Aabb] = [
        Aabb(0.0, 0.0, float(active.width), float(active.height))
    ]
    for entry in active.feature.entries:
        bounds.append(
            Aabb(
                float(active.feature_x) + float(entry.rect_x),
                float(active.feature_y) + float(entry.rect_y) - (0.5 * rect_size),
                float(active.feature_x) + float(entry.rect_x) + rect_size,
                float(active.feature_y) + float(entry.rect_y) + (0.5 * rect_size),
            ).expanded(_entry_half_stroke_width(entry.properties))
        )
    if active.gradient is not None:
        for entry in active.gradient.entries:
            bounds.append(
                Aabb(
                    float(active.gradient_x) + float(entry.bar_x),
                    float(active.gradient_y) + float(entry.bar_y) - (0.5 * rect_size),
                    float(active.gradient_x) + float(entry.bar_x) + float(active.gradient.bar_width),
                    float(active.gradient_y) + float(entry.bar_y) + (0.5 * rect_size),
                ).expanded(_entry_half_stroke_width(entry.properties))
            )
    combined = union_aabbs(bounds)
    if combined is None:  # pragma: no cover - the layout rectangle is always present
        raise RuntimeError("linear legend bounds unexpectedly empty")
    return combined


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
        local_bounds: Aabb | None = None,
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
            local_bounds=(
                local_bounds
                if local_bounds is not None
                else Aabb(0.0, 0.0, float(legend_width), float(legend_height))
            ),
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
        *,
        placement: str,
        wrap_width: float,
    ) -> LegendMeasurement:
        if not legend_table:
            return self._empty_measurement()
        wrap_width = float(wrap_width)
        if not math.isfinite(wrap_width) or wrap_width < 0.0:
            raise ValueError("legend wrap width must be a finite non-negative number")
        _legend_half_stroke_width(legend_table)
        if isinstance(self._profile, CircularRenderProfile):
            layout = build_circular_legend_layout(
                legend_table,
                legend_position=str(placement),
                canvas_width=wrap_width,
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
                local_bounds=_circular_legend_local_bounds(
                    layout,
                    color_rect_size=float(self.color_rect_size),
                ),
            )
        if not isinstance(self._profile, LinearRenderProfile):
            raise TypeError(
                "Legend measurement requires a CircularRenderProfile or LinearRenderProfile."
            )
        layout = build_linear_legend_layout(
            legend_table,
            legend_position=str(placement),
            canvas_width=wrap_width,
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
            local_bounds=_linear_legend_local_bounds(
                layout,
                color_rect_size=float(self.color_rect_size),
            ),
        )


__all__ = ["LegendDrawingConfigurator", "LegendMeasurement"]
