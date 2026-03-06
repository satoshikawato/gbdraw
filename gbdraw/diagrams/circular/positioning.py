#!/usr/bin/env python
# coding: utf-8

"""Positioning helpers for circular diagram assembly."""

from __future__ import annotations

from typing import Literal, cast

from svgwrite.container import Group  # type: ignore[reportMissingImports]

from ...canvas import CircularCanvasConfigurator  # type: ignore[reportMissingImports]

DefinitionPosition = Literal["center", "top", "bottom"]
_SUPPORTED_DEFINITION_POSITIONS = {"center", "top", "bottom"}
DEFAULT_DEFINITION_EDGE_MARGIN_PX = 24.0


def normalize_definition_position(
    position: str,
    *,
    argument_name: str = "definition_position",
) -> DefinitionPosition:
    normalized = str(position).strip().lower()
    if normalized not in _SUPPORTED_DEFINITION_POSITIONS:
        raise ValueError(
            f"{argument_name} must be one of: center, top, bottom"
        )
    return cast(DefinitionPosition, normalized)


def _parse_svg_number(value: object, *, default: float = 0.0) -> float:
    if value is None:
        return float(default)
    raw = str(value).strip()
    if raw.endswith("px"):
        raw = raw[:-2]
    try:
        return float(raw)
    except (TypeError, ValueError):
        return float(default)


def _definition_group_vertical_bounds(group: Group) -> tuple[float, float]:
    min_y: float | None = None
    max_y: float | None = None
    for element in getattr(group, "elements", []):
        attribs = getattr(element, "attribs", None)
        if not isinstance(attribs, dict):
            continue
        if "y" not in attribs:
            continue
        y_value = _parse_svg_number(attribs.get("y"), default=0.0)
        font_size = _parse_svg_number(attribs.get("font-size"), default=0.0)
        half_height = 0.5 * font_size if font_size > 0 else 0.0
        top = y_value - half_height
        bottom = y_value + half_height
        min_y = top if min_y is None else min(min_y, top)
        max_y = bottom if max_y is None else max(max_y, bottom)
    if min_y is None or max_y is None:
        return 0.0, 0.0
    return float(min_y), float(max_y)


def _definition_edge_translate_y(
    *,
    position: DefinitionPosition,
    canvas_height: float,
    bounds_min_y: float,
    bounds_max_y: float,
    margin_px: float,
) -> float:
    if position == "top":
        return float(margin_px) - float(bounds_min_y)
    return float(canvas_height) - float(margin_px) - float(bounds_max_y)


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


def place_definition_group_on_canvas(
    group: Group,
    canvas_config: CircularCanvasConfigurator,
    *,
    position: DefinitionPosition | str = "center",
    margin_px: float = DEFAULT_DEFINITION_EDGE_MARGIN_PX,
) -> Group:
    normalized = normalize_definition_position(str(position))
    if normalized == "center":
        return center_group_on_canvas(group, canvas_config)

    bounds_min_y, bounds_max_y = _definition_group_vertical_bounds(group)
    translate_x = 0.5 * float(canvas_config.total_width)
    translate_y = _definition_edge_translate_y(
        position=normalized,
        canvas_height=float(canvas_config.total_height),
        bounds_min_y=bounds_min_y,
        bounds_max_y=bounds_max_y,
        margin_px=float(margin_px),
    )
    group.translate(translate_x, translate_y)
    return group


def place_definition_group_on_size(
    group: Group,
    *,
    canvas_width: float,
    canvas_height: float,
    position: DefinitionPosition | str = "center",
    center_x: float | None = None,
    center_y: float | None = None,
    margin_px: float = DEFAULT_DEFINITION_EDGE_MARGIN_PX,
) -> Group:
    normalized = normalize_definition_position(str(position))
    if normalized == "center":
        group.translate(
            float(center_x) if center_x is not None else (0.5 * float(canvas_width)),
            float(center_y) if center_y is not None else (0.5 * float(canvas_height)),
        )
        return group

    bounds_min_y, bounds_max_y = _definition_group_vertical_bounds(group)
    translate_x = 0.5 * float(canvas_width)
    translate_y = _definition_edge_translate_y(
        position=normalized,
        canvas_height=float(canvas_height),
        bounds_min_y=bounds_min_y,
        bounds_max_y=bounds_max_y,
        margin_px=float(margin_px),
    )
    group.translate(translate_x, translate_y)
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
    "normalize_definition_position",
    "place_definition_group_on_canvas",
    "place_definition_group_on_size",
    "place_legend_on_canvas",
]
