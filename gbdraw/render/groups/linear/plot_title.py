#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations

from gbdraw.svg.elements import Group, Text

from ....core.text import measure_mixed_text_line, parse_mixed_content_text
from ....layout.spatial import Aabb
from ....render.group_content import GroupContent
from ....svg.text_path import add_text_runs


def plot_title_content(
    title: str,
    *,
    font_family: str,
    font_size: float,
) -> GroupContent:
    """The ``plot_title_content`` stage: one centered, mixed-style title line."""

    parts = parse_mixed_content_text(str(title or "").strip())
    styled = [
        (str(part_text), "italic" if part.get("italic") else "normal")
        for part in parts
        if (part_text := part.get("text"))
    ]
    # The painted box of the drawn tspans, as Circular titles measure theirs (OV-1608).
    width, height = measure_mixed_text_line(styled or [(" ", "normal")], font_family, float(font_size))
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
    add_text_runs(text, styled)
    group = Group(id="plot_title")
    group.add(text)
    return GroupContent.of(
        group,
        local_bounds=Aabb(-0.5 * float(width), -0.5 * float(height), 0.5 * float(width), 0.5 * float(height)),
    )


__all__ = ["plot_title_content"]
