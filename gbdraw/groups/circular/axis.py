#!/usr/bin/env python
# coding: utf-8

from svgwrite.container import Group
from svgwrite.shapes import Circle

from ...config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from ...svg.circular_tracks import draw_circle_path


class AxisGroup:
    """
    Creates a group for the axis on a circular canvas.
    """

    def __init__(
        self,
        radius: float,
        config_dict: dict,
        canvas_config: dict,
        cfg: GbdrawConfig | None = None,
    ) -> None:
        self.radius: float = radius
        self.axis_group = Group(id="Axis")
        self.config_dict = config_dict
        self.canvas_config = canvas_config
        self.length_param = self.canvas_config.length_param
        cfg = cfg or GbdrawConfig.from_dict(config_dict)
        self.stroke_color = cfg.objects.axis.circular.stroke_color
        self.stroke_width = cfg.objects.axis.circular.stroke_width.for_length_param(self.length_param)
        self.add_elements_to_group()

    def draw_circular_axis(self) -> None:
        self.circular_axis: Circle = draw_circle_path(self.radius, self.stroke_color, self.stroke_width)

    def add_elements_to_group(self) -> None:
        self.draw_circular_axis()
        self.axis_group.add(self.circular_axis)

    def get_group(self) -> Group:
        return self.axis_group


__all__ = ["AxisGroup"]


