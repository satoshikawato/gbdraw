#!/usr/bin/env python
# coding: utf-8

"""Positioning helpers for linear diagram assembly."""

from __future__ import annotations

from svgwrite.container import Group  # type: ignore[reportMissingImports]

from ...canvas import LinearCanvasConfigurator  # type: ignore[reportMissingImports]


def position_record_group(record_group: Group, offset: float, offset_x: float, canvas_config: LinearCanvasConfigurator) -> Group:
    """Translates the record group to its designated position on the canvas."""
    record_group.translate(offset_x + canvas_config.horizontal_offset, offset)
    return record_group


def position_linear_track_group(
    group: Group,
    offset_y: float,
    offset_x: float,
    canvas_config: LinearCanvasConfigurator,
    track_offset_y: float,
) -> Group:
    """Positions a custom linear track group relative to the record axis."""
    group.translate(
        offset_x + canvas_config.horizontal_offset,
        offset_y + float(track_offset_y),
    )
    return group


def position_record_definition_group(
    record_definition_group: Group, offset: float, offset_x: float, canvas_config: LinearCanvasConfigurator
) -> Group:
    """Places the record definition group in the correct position on the canvas."""
    record_definition_group.translate(canvas_config.horizontal_offset - offset_x, offset)
    return record_definition_group


__all__ = [
    "position_linear_track_group",
    "position_record_group",
    "position_record_definition_group",
]


