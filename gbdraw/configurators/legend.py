#!/usr/bin/env python
# coding: utf-8

from gbdraw.config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from gbdraw.core.text import calculate_bbox_dimensions  # type: ignore[reportMissingImports]
from gbdraw.legend.circular_layout import build_circular_legend_layout  # type: ignore[reportMissingImports]


class LegendDrawingConfigurator:
    def __init__(
        self,
        color_table,
        default_colors,
        selected_features_set,
        config_dict,
        gc_config,
        skew_config,
        feature_config,
        legend_table=None,
        blast_config=None,
        canvas_config=None,
        cfg: GbdrawConfig | None = None,
    ) -> None:
        cfg = cfg or GbdrawConfig.from_dict(config_dict)
        self.color_table = color_table
        self.default_colors = default_colors
        self.selected_features_set = selected_features_set
        self.config_dict = config_dict
        self.canvas_config = canvas_config
        self.show_gc = cfg.canvas.show_gc
        self.show_skew = cfg.canvas.show_skew
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
        self.num_of_columns: int = 1
        self.num_of_lines: int = 1
        self.num_of_items_per_line: int = 1
        self.has_gradient: bool = False
        self.dpi: int = cfg.canvas.dpi
        self.legend_width: float = 0
        self.total_feature_legend_width: float = 0
        self.pairwise_legend_width: float = 0

    def calculate_max_bbox_dimensions(self, legend_table):
        longest_key = max(legend_table.keys(), key=len)
        bbox_width_px, bbox_height_px = calculate_bbox_dimensions(
            longest_key, self.font_family, self.font_size, self.dpi
        )
        return bbox_width_px, bbox_height_px

    def recalculate_legend_dimensions(self, legend_table, canvas_config):
        if not legend_table:
            self.has_gradient = False
            self.pairwise_legend_width = 0
            self.legend_width = 0
            self.legend_height = 0
            self.total_feature_legend_width = 0
            self.num_of_lines = 0
            self.num_of_columns = 0
            self.num_of_items_per_line = 0
            return self
        if canvas_config.__class__.__name__ == "CircularCanvasConfigurator":
            layout = build_circular_legend_layout(
                legend_table,
                legend_position=canvas_config.legend_position,
                canvas_width=float(canvas_config.total_width),
                font_family=self.font_family,
                font_size=float(self.font_size),
                dpi=int(self.dpi),
                color_rect_size=float(self.color_rect_size),
            )
            self.has_gradient = layout.has_gradient
            self.pairwise_legend_width = layout.pairwise_legend_width
            self.legend_width = layout.width
            self.legend_height = layout.height
            self.total_feature_legend_width = layout.feature_width
            self.num_of_lines = layout.num_lines
            self.num_of_columns = layout.num_columns
            self.num_of_items_per_line = layout.num_items_per_line
            return self

        line_margin = (24 / 14) * self.color_rect_size  # move to config
        x_margin = (22 / 14) * self.color_rect_size  # move to config
        total_width = canvas_config.total_width
        gradient_count = sum(
            1 for properties in legend_table.values() if properties.get("type") == "gradient"
        )
        if gradient_count:
            self.has_gradient = True
            self.pairwise_legend_width = 10 * self.color_rect_size
            if gradient_count > 1:
                gradient_label_width = max(
                    calculate_bbox_dimensions(
                        str(key),
                        self.font_family,
                        self.font_size,
                        self.dpi,
                    )[0]
                    for key, properties in legend_table.items()
                    if properties.get("type") == "gradient"
                )
                self.pairwise_legend_width += gradient_label_width + x_margin
        if canvas_config.legend_position == "top" or canvas_config.legend_position == "bottom":
            bbox_list = [
                calculate_bbox_dimensions(item, self.font_family, self.font_size, self.dpi)
                for item in legend_table
            ]
            total_feature_legend_width = sum([bbox[0] for bbox in bbox_list]) + 2 * x_margin * (
                len(legend_table) + 1
            )

            if self.has_gradient:
                total_legend_width = total_feature_legend_width + self.pairwise_legend_width
            else:
                total_legend_width = total_feature_legend_width
            if total_legend_width <= total_width:
                self.legend_width = total_legend_width
                feature_legend_height = self.color_rect_size + 2 * line_margin
                if gradient_count > 1:
                    _, scale_label_height = calculate_bbox_dimensions(
                        "100%",
                        self.font_family,
                        self.font_size,
                        self.dpi,
                    )
                    gradient_legend_height = (
                        self.color_rect_size
                        + (gradient_count - 1) * line_margin
                        + 2
                        + scale_label_height
                    )
                else:
                    gradient_legend_height = gradient_count * (
                        2 * self.color_rect_size + 2 * line_margin
                    )
                self.legend_height = max(feature_legend_height, gradient_legend_height)
                self.total_feature_legend_width = self.legend_width - self.pairwise_legend_width
                self.num_of_columns = len(legend_table)
                return self
            else:
                # Split into two rows or more. Add one item at a time until width exceeds total_width
                current_width = x_margin
                num_lines = 1
                per_line_item_count = 0
                for item in bbox_list:
                    item_width = item[0]
                    if current_width + item_width + self.pairwise_legend_width + x_margin <= total_width:
                        per_line_item_count += 1
                        self.num_of_items_per_line = max(per_line_item_count, self.num_of_items_per_line)
                        current_width += item_width + x_margin + self.pairwise_legend_width
                        if self.has_gradient:
                            current_num_columns = per_line_item_count + 1
                        else:
                            current_num_columns = per_line_item_count
                        if current_num_columns > self.num_of_columns:
                            self.num_of_columns = current_num_columns
                        self.legend_width = max(current_width, self.legend_width)
                    else:
                        num_lines += 1
                        self.num_of_lines = max(self.num_of_lines, num_lines)
                        per_line_item_count = 0
                        current_width = x_margin

                self.total_feature_legend_width = self.legend_width - self.pairwise_legend_width
                self.legend_height = num_lines * (self.color_rect_size + line_margin)
                self.num_of_columns = len(legend_table)
                return self

        else:
            bbox_width_px, _ = self.calculate_max_bbox_dimensions(legend_table)
            self.total_feature_legend_width = x_margin + bbox_width_px
            self.legend_width = max(self.total_feature_legend_width, self.pairwise_legend_width)
            num_lines = 0
            for key, properties in legend_table.items():
                if properties.get("type") == "gradient":
                    self.has_gradient = True
                else:
                    num_lines += 1

            if self.has_gradient:
                if gradient_count > 1:
                    num_lines += gradient_count + 1
                else:
                    num_lines += 2 * gradient_count

            if num_lines > 0:
                self.legend_height = self.color_rect_size + (num_lines - 1) * line_margin
            else:
                self.legend_height = 0
            # self.legend_height = (self.color_rect_size + (len(legend_table.keys()) -1) * line_margin)
            return self


__all__ = ["LegendDrawingConfigurator"]


