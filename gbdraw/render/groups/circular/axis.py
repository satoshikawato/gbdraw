#!/usr/bin/env python
# coding: utf-8

from gbdraw.svg.elements import Circle, Group, Path

from ....canvas import CircularCanvasGeometry
from ....config.models import GbdrawConfig
from ....svg.circular_arcs import circle_arc
from ....svg.circular_tracks import draw_circle_path


class AxisGroup:
    """
    Creates a group for the axis on a circular canvas.
    """

    def __init__(
        self,
        radius: float,
        geometry: CircularCanvasGeometry,
        *,
        cfg: GbdrawConfig,
        unit_arc: tuple[float, float] | None = None,
    ) -> None:
        self.radius: float = radius
        # The arc (start and span in degrees) of a record that shares its circle with other records.
        self.unit_arc = unit_arc
        self.axis_group = Group(id="Axis")
        self.length_param = geometry.length_param
        self.stroke_color = cfg.objects.axis.circular.stroke_color
        self.stroke_width = cfg.objects.axis.circular.stroke_width.for_length_param(self.length_param)
        self.add_elements_to_group()

    def draw_circular_axis(self) -> None:
        self.circular_axis: Circle | Path = (
            draw_circle_path(self.radius, self.stroke_color, self.stroke_width)
            if self.unit_arc is None
            else Path(
                d=circle_arc(self.radius, *self.unit_arc),
                stroke=self.stroke_color,
                stroke_width=self.stroke_width,
                fill="none",
            )
        )

    def add_elements_to_group(self) -> None:
        self.draw_circular_axis()
        self.axis_group.add(self.circular_axis)

    def get_group(self) -> Group:
        return self.axis_group


__all__ = ["AxisGroup"]


