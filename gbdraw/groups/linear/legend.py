#!/usr/bin/env python
# coding: utf-8

from svgwrite.container import Group
from svgwrite.path import Path
from svgwrite.gradients import LinearGradient

from ...text import calculate_bbox_dimensions
from ...svg.text_path import generate_text_path
from ...config.models import GbdrawConfig  # type: ignore[reportMissingImports]


class LegendGroup:
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
        self.line_height: float = (24/14) * self.rect_size
        self.text_x_offset: float = (22/14) * self.rect_size
        self.num_of_columns = self.legend_config.num_of_columns
        self.has_gradient = self.legend_config.has_gradient
        self.total_feature_legend_width = self.legend_config.total_feature_legend_width
        self.dpi: int = cfg.canvas.dpi
        self.legend_width: float = 0
        self.feature_legend_width: float = 0
        self.feature_legend_height: float = 0
        self.pairwise_legend_width = self.legend_config.pairwise_legend_width
        self.pairwise_legend_height: float = 0
        self.legend_height: float = 0
        self.add_elements_to_group()

    def create_rectangle_path_for_legend(self) -> str:
        start_y_top = -self.rect_size / 2
        start_y_bottom = self.rect_size / 2
        return f"M 0,{start_y_top} L {self.rect_size},{start_y_top} L {self.rect_size},{start_y_bottom} L 0,{start_y_bottom} z"

    def add_elements_to_group(self):
        self.feature_legend_height = self.rect_size
        y_offset = self.rect_size / 2
        initial_y_offset = y_offset
        path_desc = self.create_rectangle_path_for_legend()
        font = self.font_family
        feature_legend_group = Group(id="feature_legend")
        max_bbox_width = 0
        gradient_y_offset = 0
        grad_bar_width = 0
        current_feature_legend_width = 0
        pairwise_legend_group = Group(id="pairwise_legend")
        if self.num_of_columns > 1:
            if self.has_gradient:
                max_items_per_column = self.num_of_columns - 1
            else:
                max_items_per_column = self.num_of_columns
            current_column = 0
            current_x_offset = 0
            for key, properties in self.legend_table.items():
                if properties['type'] == 'solid':
                    bbox_width, _ = calculate_bbox_dimensions(str(key), self.font_family, self.font_size, self.dpi)
                    max_bbox_width = max(max_bbox_width, bbox_width)
                    if current_feature_legend_width + self.text_x_offset + bbox_width + self.text_x_offset > self.total_feature_legend_width:
                        current_feature_legend_width = 0
                        current_column = 0
                        current_x_offset = 0
                        y_offset += self.line_height
                        self.feature_legend_height += self.line_height
                    rect_path = Path(
                        d=path_desc,
                        fill=properties['fill'],
                        stroke=properties['stroke'],
                        stroke_width=properties['width'])
                    rect_path.translate(current_x_offset, y_offset)
                    feature_legend_group.add(rect_path)
                    current_x_offset +=  self.text_x_offset
                    legend_path = generate_text_path(
                            key, 0, 0, 0, self.font_size, "normal", font, 
                            dominant_baseline='central', text_anchor="start"
                        )
                    legend_path.translate(current_x_offset, y_offset)
                    feature_legend_group.add(legend_path)
                    current_x_offset += (bbox_width + self.text_x_offset)
                    current_feature_legend_width = current_x_offset 
                    self.feature_legend_width = max(current_feature_legend_width, self.feature_legend_width)
                    current_column += 1

                elif properties['type'] == 'gradient':
                    grad_bar_x_start = 0
                    title_x_offset = 0
                    grad_bar_width = self.legend_config.pairwise_legend_width
                    gradient_y_offset = 0
                    gradient_id = f"blast_legend_grad_{abs(hash(properties['min_color'] + properties['max_color']))}"
                    
                    gradient = LinearGradient(start=(0, 0), end=("100%", 0), id=gradient_id)
                    gradient.add_stop_color(offset="0%", color=properties['min_color'])
                    gradient.add_stop_color(offset="100%", color=properties['max_color'])
                    pairwise_legend_group.add(gradient)
                    title_path = generate_text_path(
                        key, 0, 0, 0, self.font_size, "normal", font, 
                        dominant_baseline='hanging', text_anchor="middle"
                    )
                    title_path_bbox_width, title_path_bbox_height = calculate_bbox_dimensions(str(key), self.font_family, self.font_size, self.dpi)
                    title_x_offset = grad_bar_x_start + (grad_bar_width / 2)
                    title_path.translate(title_x_offset, gradient_y_offset) 
                    pairwise_legend_group.add(title_path) 
                    gradient_y_offset += title_path_bbox_height + (self.rect_size / 2)
                    
                    grad_rect_path_desc = f"M 0,{-self.rect_size / 2} L {grad_bar_width},{-self.rect_size / 2} L {grad_bar_width},{self.rect_size / 2} L 0,{self.rect_size / 2} z"
                    grad_rect = Path(
                        d=grad_rect_path_desc,
                        fill=f"url(#{gradient_id})",
                        stroke=properties['stroke'],
                        stroke_width=properties['width']
                    )
                    grad_rect.translate(grad_bar_x_start, gradient_y_offset) 
                    pairwise_legend_group.add(grad_rect)
                    min_identity = properties.get('min_value', 0)
                    if min_identity == int(min_identity):
                        min_label_text = f"{int(min_identity)}%"
                    else:
                        min_label_text = f"{min_identity}%"
                    label_0 = generate_text_path(
                        min_label_text, 0, 0, 0, self.font_size, "normal", font,
                        dominant_baseline='hanging', text_anchor="start"
                     )
                    label_0.translate(grad_bar_x_start, gradient_y_offset + self.rect_size / 2 + 2) 
                    pairwise_legend_group.add(label_0)
                    _, min_label_bbox_height = calculate_bbox_dimensions(min_label_text, self.font_family, self.font_size, self.dpi)

                    label_100 = generate_text_path(
                        "100%", 0, 0, 0, self.font_size, "normal", font,
                        dominant_baseline='hanging', text_anchor="end"
                    )
                    label_100.translate(grad_bar_x_start + grad_bar_width, gradient_y_offset + self.rect_size / 2 + 2) 
                    pairwise_legend_group.add(label_100)
                    _, max_label_bbox_height = calculate_bbox_dimensions("100%", self.font_family, self.font_size, self.dpi)
                    pairwise_legend_group.add(label_100)
                    label_bbox_height = max(min_label_bbox_height, max_label_bbox_height)
                    gradient_y_offset += (label_bbox_height + self.rect_size / 2 + 2)

                    self.pairwise_legend_width = grad_bar_width
                    self.pairwise_legend_height = gradient_y_offset
            if self.pairwise_legend_width >0:
                self.legend_width = self.feature_legend_width + self.pairwise_legend_width + self.text_x_offset
            else:
                self.legend_width = self.feature_legend_width
            self.legend_height = max(self.feature_legend_height, self.pairwise_legend_height)
            
            if self.feature_legend_height > self.pairwise_legend_height:
                pairwise_legend_group.translate(self.feature_legend_width + self.text_x_offset, (self.feature_legend_height - self.pairwise_legend_height)/2)
            else:
                pairwise_legend_group.translate(self.feature_legend_width  + self.text_x_offset, 0)
                feature_legend_group.translate(0, (self.pairwise_legend_height - self.feature_legend_height)/2)

            self.legend_group.add(pairwise_legend_group)
            self.legend_group.add(feature_legend_group)
        else:
            for key, properties in self.legend_table.items():
                if properties['type'] == 'solid':
                    current_feature_legend_width = 0
                    rect_path = Path(
                        d=path_desc,
                        fill=properties['fill'],
                        stroke=properties['stroke'],
                        stroke_width=properties['width'])
                    rect_path.translate(0, y_offset)
                    feature_legend_group.add(rect_path)
                    bbox_width, bbox_height = calculate_bbox_dimensions(key, self.font_family, self.font_size, self.dpi)
                    max_bbox_width = max(max_bbox_width, bbox_width)
                    legend_path = generate_text_path(
                        key, 0, 0, 0, self.font_size, "normal", font, 
                        dominant_baseline='central', text_anchor="start"
                    )
                    legend_path.translate(self.text_x_offset, y_offset)
                    current_feature_legend_width = self.text_x_offset + bbox_width
                    feature_legend_group.add(legend_path)
                    
                    y_offset += self.line_height 
                    self.feature_legend_width = max(current_feature_legend_width, self.feature_legend_width)
                    
                elif properties['type'] == 'gradient':
                    grad_bar_x_start = 0
                    if self.legend_position == 'top' or self.legend_position == 'bottom':
                        gradient_y_offset = 0
                    else:
                        gradient_y_offset = y_offset
                    grad_bar_width = self.pairwise_legend_width
                    gradient_id = f"blast_legend_grad_{abs(hash(properties['min_color'] + properties['max_color']))}"
                    
                    gradient = LinearGradient(start=(0, 0), end=("100%", 0), id=gradient_id)
                    gradient.add_stop_color(offset="0%", color=properties['min_color'])
                    gradient.add_stop_color(offset="100%", color=properties['max_color'])
                    pairwise_legend_group.add(gradient)
                    title_path = generate_text_path(
                        key, 0, 0, 0, self.font_size, "normal", font, 
                        dominant_baseline='hanging', text_anchor="middle"
                    )
                    title_path_bbox_width, title_path_bbox_height = calculate_bbox_dimensions(key, self.font_family, self.font_size, self.dpi)
                    title_x_offset = grad_bar_x_start + (grad_bar_width / 2)

                    title_path.translate(title_x_offset, gradient_y_offset) 
                    pairwise_legend_group.add(title_path) 
                    gradient_y_offset += (title_path_bbox_height + self.rect_size / 2)  # self.line_height 


                    
                    grad_rect_path_desc = f"M 0,{-self.rect_size / 2} L {grad_bar_width},{-self.rect_size / 2} L {grad_bar_width},{self.rect_size / 2} L 0,{self.rect_size / 2} z"
                    grad_rect = Path(
                        d=grad_rect_path_desc,
                        fill=f"url(#{gradient_id})",
                        stroke=properties['stroke'],
                        stroke_width=properties['width']
                    )
                    grad_rect.translate(grad_bar_x_start, gradient_y_offset) 
                    pairwise_legend_group.add(grad_rect)
                    min_identity = properties.get('min_value', 0)
                    if min_identity == int(min_identity):
                        min_label_text = f"{int(min_identity)}%"
                    else:
                        min_label_text = f"{min_identity}%"
                    label_0 = generate_text_path(
                        min_label_text, 0, 0, 0, self.font_size, "normal", font,
                        dominant_baseline='hanging', text_anchor="start"
                    )
                    label_0.translate(grad_bar_x_start, gradient_y_offset + self.rect_size / 2 + 2)
                    _, min_label_bbox_height = calculate_bbox_dimensions(min_label_text, self.font_family, self.font_size, self.dpi)
                    pairwise_legend_group.add(label_0)
                    label_100 = generate_text_path(
                        "100%", 0, 0, 0, self.font_size, "normal", font,
                        dominant_baseline='hanging', text_anchor="end"
                    )
                    label_100.translate(grad_bar_x_start + grad_bar_width, gradient_y_offset + self.rect_size / 2 + 2) 
                    _, max_label_bbox_height = calculate_bbox_dimensions("100%", self.font_family, self.font_size, self.dpi)
                    pairwise_legend_group.add(label_100)
                    label_bbox_height = max(min_label_bbox_height, max_label_bbox_height)
                    gradient_y_offset += label_bbox_height 
            if self.pairwise_legend_width > 0:
                if self.pairwise_legend_width > self.feature_legend_width:
                    feature_legend_group.translate((self.pairwise_legend_width - self.feature_legend_width)/2, 0)
                    self.legend_group.add(pairwise_legend_group)
                else:
                    self.legend_group.add(pairwise_legend_group)
            self.legend_width = max(self.feature_legend_width, self.pairwise_legend_width)
            self.legend_height = max(gradient_y_offset, y_offset) + initial_y_offset
            self.legend_group.add(feature_legend_group)
        return self.legend_group
    def get_group(self) -> Group:
        """
        Retrieves the SVG group containing the figure legends.

        Returns:
            Group: The SVG group with figure legends.
        """
        return self.legend_group


__all__ = ["LegendGroup"]


