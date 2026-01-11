#!/usr/bin/env python
# coding: utf-8

"""Positioning helpers for circular diagram assembly."""

from __future__ import annotations

from svgwrite.container import Group  # type: ignore[reportMissingImports]

from ...canvas import CircularCanvasConfigurator  # type: ignore[reportMissingImports]


def center_group_on_canvas(group: Group, canvas_config: CircularCanvasConfigurator) -> Group:
    """
    Centers a given SVG group on the canvas based on the canvas configuration.

    Parameters:
    group (Group): The SVG group to be centered.
    canvas_config (CircularCanvasConfigurator): The configuration of the circular canvas.

    Returns:
    Group: The centered SVG group.
    """
    group.translate(canvas_config.offset_x, canvas_config.offset_y)
    return group


def place_legend_on_canvas(group: Group, canvas_config: CircularCanvasConfigurator) -> Group:
    """
    Places a legend group at the configured legend position on the canvas.

    Parameters:
    group (Group): The SVG group containing the legend.
    canvas_config (CircularCanvasConfigurator): The configuration of the circular canvas.

    Returns:
    Group: The positioned legend group.
    """
    group.translate(canvas_config.legend_offset_x, canvas_config.legend_offset_y)
    return group


__all__ = [
    "center_group_on_canvas",
    "place_legend_on_canvas",
]
