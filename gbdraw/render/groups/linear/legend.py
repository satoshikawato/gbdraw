#!/usr/bin/env python
# coding: utf-8

from svgwrite.container import Group
from svgwrite.path import Path
from svgwrite.gradients import LinearGradient

from ....core.text import calculate_bbox_dimensions
from ....svg.text_path import generate_text_path
from ....config.models import GbdrawConfig  # type: ignore[reportMissingImports]


class LegendGroup:
    """
    Generates both horizontal and vertical legend layouts for linear diagrams.
    Both layouts are always generated and added to the legend group.
    The appropriate layout is shown/hidden based on legend_position.
    """

    def __init__(
        self,
        config_dict,
        canvas_config,
        legend_config,
        legend_table,
        cfg: GbdrawConfig | None = None,
    ):
        self.legend_group = Group(id="legend")
        self.config_dict = config_dict
        cfg = cfg or GbdrawConfig.from_dict(config_dict)
        self._cfg = cfg
        self.canvas_config = canvas_config
        self.legend_config = legend_config
        self.legend_table = legend_table
        self.font_family: str = self.legend_config.font_family
        self.font_weight: str = self.legend_config.font_weight
        self.font_size: float = self.legend_config.font_size
        self.rect_size: float = self.legend_config.color_rect_size
        self.legend_position = self.canvas_config.legend_position
        self.line_height: float = (24 / 14) * self.rect_size
        self.text_x_offset: float = (22 / 14) * self.rect_size
        self.num_of_columns = self.legend_config.num_of_columns
        self.has_gradient = self.legend_config.has_gradient
        self.total_feature_legend_width = self.legend_config.total_feature_legend_width
        self.dpi: int = cfg.canvas.dpi
        self.pairwise_legend_width = self.legend_config.pairwise_legend_width

        # Dimensions for both layouts
        self.horizontal_legend_width: float = 0
        self.horizontal_legend_height: float = 0
        self.vertical_legend_width: float = 0
        self.vertical_legend_height: float = 0

        # Current active legend dimensions (based on position)
        self.legend_width: float = 0
        self.legend_height: float = 0

        self._build_dual_legends()

    def create_rectangle_path_for_legend(self) -> str:
        start_y_top = -self.rect_size / 2
        start_y_bottom = self.rect_size / 2
        return f"M 0,{start_y_top} L {self.rect_size},{start_y_top} L {self.rect_size},{start_y_bottom} L 0,{start_y_bottom} z"

    def _calculate_entry_widths(self) -> list[float]:
        """Pre-calculate widths for each legend entry."""
        entry_widths = []
        for key, properties in self.legend_table.items():
            if properties["type"] == "solid":
                bbox_width, _ = calculate_bbox_dimensions(
                    str(key), self.font_family, self.font_size, self.dpi
                )
                entry_widths.append(
                    self.rect_size + self.text_x_offset + bbox_width + self.text_x_offset
                )
        return entry_widths

    def _build_horizontal_feature_legend(self) -> tuple[Group, float, float]:
        """Build horizontal layout (entries side by side, for top/bottom positions)."""
        group = Group()
        path_desc = self.create_rectangle_path_for_legend()
        font = self.font_family

        y_offset = self.rect_size / 2
        current_x_offset = 0
        max_height = self.line_height
        current_row_width = 0
        max_row_width = 0

        # For horizontal layout, use canvas width as the wrap limit, not total_feature_legend_width
        # total_feature_legend_width may be calculated for vertical layout (single item width)
        # when legend_position is left/right at generation time
        horizontal_wrap_width = self.canvas_config.total_width if self.canvas_config else 0
        if self.has_gradient and horizontal_wrap_width > 0:
            reserved_width = self.pairwise_legend_width + self.text_x_offset
            min_width = self.rect_size + self.text_x_offset * 2
            horizontal_wrap_width = max(horizontal_wrap_width - reserved_width, min_width)

        for key, properties in self.legend_table.items():
            if properties["type"] != "solid":
                continue

            bbox_width, _ = calculate_bbox_dimensions(
                str(key), self.font_family, self.font_size, self.dpi
            )
            entry_width = self.rect_size + self.text_x_offset + bbox_width + self.text_x_offset

            # Check if we need to wrap to next row (use canvas width for horizontal layout)
            if (
                horizontal_wrap_width > 0
                and current_x_offset + entry_width > horizontal_wrap_width
            ):
                max_row_width = max(max_row_width, current_row_width)
                current_x_offset = 0
                current_row_width = 0
                y_offset += self.line_height
                max_height += self.line_height

            # Create entry group with data attribute for identification
            entry_group = Group(debug=False)
            entry_group.attribs["data-legend-key"] = str(key)

            # Add rectangle
            rect_path = Path(
                d=path_desc,
                fill=properties["fill"],
                stroke=properties["stroke"],
                stroke_width=properties["width"],
            )
            rect_path.translate(current_x_offset, y_offset)
            entry_group.add(rect_path)

            # Add text
            legend_path = generate_text_path(
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
            legend_path.translate(current_x_offset + self.text_x_offset, y_offset)
            entry_group.add(legend_path)
            group.add(entry_group)

            current_x_offset += entry_width
            current_row_width += entry_width

        max_row_width = max(max_row_width, current_row_width)
        width = max_row_width
        height = max_height

        return group, width, height

    def _build_vertical_feature_legend(self) -> tuple[Group, float, float]:
        """Build vertical layout (entries stacked, for left/right positions)."""
        group = Group()
        path_desc = self.create_rectangle_path_for_legend()
        font = self.font_family

        y_offset = self.rect_size / 2
        max_width = 0

        for key, properties in self.legend_table.items():
            if properties["type"] != "solid":
                continue

            bbox_width, _ = calculate_bbox_dimensions(
                str(key), self.font_family, self.font_size, self.dpi
            )

            # Create entry group with data attribute for identification
            entry_group = Group(debug=False)
            entry_group.attribs["data-legend-key"] = str(key)

            # Add rectangle
            rect_path = Path(
                d=path_desc,
                fill=properties["fill"],
                stroke=properties["stroke"],
                stroke_width=properties["width"],
            )
            rect_path.translate(0, y_offset)
            entry_group.add(rect_path)

            # Add text
            legend_path = generate_text_path(
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
            legend_path.translate(self.text_x_offset, y_offset)
            entry_group.add(legend_path)
            group.add(entry_group)

            entry_width = self.text_x_offset + bbox_width
            max_width = max(max_width, entry_width)
            y_offset += self.line_height

        height = y_offset
        width = max_width

        return group, width, height

    def _build_pairwise_legend(self) -> tuple[Group, float, float]:
        """Build pairwise (gradient) legend for BLAST comparisons."""
        group = Group(id="pairwise_legend")
        font = self.font_family
        grad_bar_width = self.pairwise_legend_width

        gradient_y_offset = 0

        for key, properties in self.legend_table.items():
            if properties["type"] != "gradient":
                continue

            gradient_id = f"blast_legend_grad_{abs(hash(properties['min_color'] + properties['max_color']))}"

            # Create gradient
            gradient = LinearGradient(start=(0, 0), end=("100%", 0), id=gradient_id)
            gradient.add_stop_color(offset="0%", color=properties["min_color"])
            gradient.add_stop_color(offset="100%", color=properties["max_color"])
            group.add(gradient)

            # Title
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
            _, title_path_bbox_height = calculate_bbox_dimensions(
                str(key), self.font_family, self.font_size, self.dpi
            )
            title_x_offset = grad_bar_width / 2
            title_path.translate(title_x_offset, gradient_y_offset)
            group.add(title_path)
            gradient_y_offset += title_path_bbox_height + (self.rect_size / 2)

            # Gradient rectangle
            grad_rect_path_desc = f"M 0,{-self.rect_size / 2} L {grad_bar_width},{-self.rect_size / 2} L {grad_bar_width},{self.rect_size / 2} L 0,{self.rect_size / 2} z"
            grad_rect = Path(
                d=grad_rect_path_desc,
                fill=f"url(#{gradient_id})",
                stroke=properties["stroke"],
                stroke_width=properties["width"],
            )
            grad_rect.translate(0, gradient_y_offset)
            group.add(grad_rect)

            # Labels
            min_identity = properties.get("min_value", 0)
            if min_identity == int(min_identity):
                min_label_text = f"{int(min_identity)}%"
            else:
                min_label_text = f"{min_identity}%"

            label_0 = generate_text_path(
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
            label_0.translate(0, gradient_y_offset + self.rect_size / 2 + 2)
            group.add(label_0)

            label_100 = generate_text_path(
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
            label_100.translate(grad_bar_width, gradient_y_offset + self.rect_size / 2 + 2)
            group.add(label_100)

            _, min_label_bbox_height = calculate_bbox_dimensions(
                min_label_text, self.font_family, self.font_size, self.dpi
            )
            _, max_label_bbox_height = calculate_bbox_dimensions(
                "100%", self.font_family, self.font_size, self.dpi
            )
            label_bbox_height = max(min_label_bbox_height, max_label_bbox_height)
            gradient_y_offset += label_bbox_height + self.rect_size / 2 + 2

        return group, grad_bar_width, gradient_y_offset

    def _build_dual_legends(self):
        """Build both horizontal and vertical legend layouts."""
        # Check if we have any legend entries
        has_solid = any(p["type"] == "solid" for p in self.legend_table.values())
        has_gradient = any(p["type"] == "gradient" for p in self.legend_table.values())

        if not has_solid and not has_gradient:
            return

        # Build feature legends for both orientations
        h_feature_group, h_feature_width, h_feature_height = (
            self._build_horizontal_feature_legend()
        )
        v_feature_group, v_feature_width, v_feature_height = (
            self._build_vertical_feature_legend()
        )

        # Build pairwise legend (same for both orientations, positioned differently)
        pairwise_width, pairwise_height = 0, 0
        if has_gradient:
            _, pairwise_width, pairwise_height = self._build_pairwise_legend()

        # Create horizontal legend group
        horizontal_group = Group(id="legend_horizontal")
        horizontal_feature_group = Group(id="feature_legend_h")
        for child in h_feature_group.elements:
            horizontal_feature_group.add(child)

        if has_gradient:
            h_pairwise_group, _, _ = self._build_pairwise_legend()
            # Position pairwise legend to the right of feature legend
            h_pairwise_group.translate(h_feature_width + self.text_x_offset, 0)
            # Vertically center both
            if h_feature_height > pairwise_height:
                h_pairwise_group.translate(0, (h_feature_height - pairwise_height) / 2)
            elif pairwise_height > h_feature_height:
                horizontal_feature_group.translate(0, (pairwise_height - h_feature_height) / 2)
            horizontal_group.add(horizontal_feature_group)
            horizontal_group.add(h_pairwise_group)
            self.horizontal_legend_width = h_feature_width + self.text_x_offset + pairwise_width
            self.horizontal_legend_height = max(h_feature_height, pairwise_height)
        else:
            horizontal_group.add(horizontal_feature_group)
            self.horizontal_legend_width = h_feature_width
            self.horizontal_legend_height = h_feature_height

        # Create vertical legend group
        vertical_group = Group(id="legend_vertical")
        vertical_feature_group = Group(id="feature_legend_v")
        for child in v_feature_group.elements:
            vertical_feature_group.add(child)

        if has_gradient:
            v_pairwise_group, _, _ = self._build_pairwise_legend()
            # Position pairwise legend below feature legend
            v_pairwise_group.translate(0, v_feature_height + self.line_height / 2)
            # Horizontally center
            if pairwise_width > v_feature_width:
                vertical_feature_group.translate((pairwise_width - v_feature_width) / 2, 0)
            elif v_feature_width > pairwise_width:
                v_pairwise_group.translate((v_feature_width - pairwise_width) / 2, 0)
            vertical_group.add(vertical_feature_group)
            vertical_group.add(v_pairwise_group)
            self.vertical_legend_width = max(v_feature_width, pairwise_width)
            self.vertical_legend_height = v_feature_height + self.line_height / 2 + pairwise_height
        else:
            vertical_group.add(vertical_feature_group)
            self.vertical_legend_width = v_feature_width
            self.vertical_legend_height = v_feature_height

        # Set visibility based on current legend position
        is_horizontal = self.legend_position in ("top", "bottom")
        if is_horizontal:
            vertical_group.attribs["display"] = "none"
            self.legend_width = self.horizontal_legend_width
            self.legend_height = self.horizontal_legend_height
        else:
            horizontal_group.attribs["display"] = "none"
            self.legend_width = self.vertical_legend_width
            self.legend_height = self.vertical_legend_height

        # Add both groups to the main legend group
        self.legend_group.add(horizontal_group)
        self.legend_group.add(vertical_group)

    def get_group(self) -> Group:
        """
        Retrieves the SVG group containing the figure legends.

        Returns:
            Group: The SVG group with figure legends.
        """
        return self.legend_group

    def get_horizontal_dimensions(self) -> tuple[float, float]:
        """Get dimensions of the horizontal legend layout."""
        return self.horizontal_legend_width, self.horizontal_legend_height

    def get_vertical_dimensions(self) -> tuple[float, float]:
        """Get dimensions of the vertical legend layout."""
        return self.vertical_legend_width, self.vertical_legend_height


__all__ = ["LegendGroup"]
