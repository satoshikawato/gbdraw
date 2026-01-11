#!/usr/bin/env python
# coding: utf-8

from svgwrite.text import Text, TSpan, TextPath


def generate_name_path(
    name_parts: list,
    title_x: float,
    title_y: float,
    interval: float,
    font_size: str,
    font_weight: str,
    font_family: str,
) -> Text:
    """
    Creates a Text element consisting of multiple spans (e.g., italic parts).
    """
    t = Text(
        "",
        insert=(title_x, title_y + interval),
        font_size=font_size,
        font_weight=font_weight,
        font_family=font_family,
        text_anchor="middle",
        dominant_baseline="middle",
    )
    for part in name_parts:
        part_text = part.get("text")
        if part_text is None or not part_text.strip():
            continue
        if part.get("italic"):
            t.add(TSpan(part_text, font_style="italic"))
        else:
            t.add(TSpan(part_text))
    return t


def generate_text_path(
    text: str,
    title_x: float,
    title_y: float,
    interval: float,
    font_size: float,
    font_weight: str,
    font: str,
    dominant_baseline: str = "middle",
    text_anchor: str = "middle",
) -> Text:
    """
    Creates a simple Text element at a position.

    NOTE: Despite the name, this is currently a plain Text factory; kept for compatibility.
    """
    return Text(
        text,
        insert=(title_x, title_y + interval),
        font_size=font_size,
        font_weight=font_weight,
        font_family=font,
        text_anchor=text_anchor,
        dominant_baseline=dominant_baseline,
    )


__all__ = ["generate_name_path", "generate_text_path"]


