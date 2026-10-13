#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations

from dataclasses import dataclass

from gbdraw.svg.elements import Group, TSpan, Text

from ....core.text import calculate_bbox_dimensions, parse_mixed_content_text
from ....layout.spatial import Aabb
from ....render.group_content import GroupContent


@dataclass(frozen=True)
class PlotTitleContent:
    """The shared plot title of a Linear canvas and its bounds around its anchor."""

    content: GroupContent
    local_bounds: Aabb


def plot_title_content(
    title: str,
    *,
    font_family: str,
    font_size: float,
    dpi: int,
) -> PlotTitleContent:
    """The ``plot_title_content`` stage: one centered, mixed-style title line."""

    parts = parse_mixed_content_text(str(title or "").strip())
    plain = "".join(str(part.get("text") or "") for part in parts)
    width, height = calculate_bbox_dimensions(plain if plain else " ", font_family, float(font_size), dpi)
    text = Text(
        "",
        insert=(0, 0),
        stroke="none",
        fill="black",
        font_size=float(font_size),
        font_weight="normal",
        font_family=font_family,
        text_anchor="middle",
        dominant_baseline="middle",
    )
    for part in parts:
        part_text = part.get("text")
        if part_text is None:
            continue
        if part.get("italic"):
            text.add(TSpan(str(part_text), font_style="italic"))
        else:
            text.add(TSpan(str(part_text)))
    group = Group(id="plot_title")
    group.add(text)
    return PlotTitleContent(
        content=GroupContent.of(group),
        local_bounds=Aabb(-0.5 * float(width), -0.5 * float(height), 0.5 * float(width), 0.5 * float(height)),
    )


__all__ = ["PlotTitleContent", "plot_title_content"]
