#!/usr/bin/env python
# coding: utf-8

from collections.abc import Iterable, Sequence
from typing import TypeVar

from gbdraw.svg.elements import Text, TextPath, TSpan

# A ``text`` or ``textPath`` element: either holds a text node or ``tspan`` runs.
TextElement = TypeVar("TextElement", Text, TextPath)


def generate_name_path(
    name_parts: list,
    title_x: float,
    title_y: float,
    interval: float,
    font_size: float,
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
    font_style: str = "normal",
) -> Text:
    """
    Creates a simple Text element at a position.

    ``font_style`` is written only when it is not the initial value ``normal``.

    NOTE: Despite the name, this is currently a plain Text factory; kept for compatibility.
    """
    style = {} if font_style == "normal" else {"font_style": font_style}
    return Text(
        text,
        insert=(title_x, title_y + interval),
        font_size=font_size,
        font_weight=font_weight,
        font_family=font,
        text_anchor=text_anchor,
        dominant_baseline=dominant_baseline,
        **style,
    )


def add_text_runs(text: TextElement, runs: Iterable[tuple[str, str]]) -> TextElement:
    """Append one ``tspan`` per ``(text, font_style)`` run to ``text``; an italic run carries ``font-style``."""

    for part, font_style in runs:
        text.add(TSpan(part, font_style="italic") if font_style == "italic" else TSpan(part))
    return text


def write_text_runs(element: TextElement, runs: Sequence[tuple[str, str]]) -> TextElement:
    """Write ``runs`` into an empty ``text`` or ``textPath``: one upright run as its plain text node, styled runs as
    ``tspan``s (``add_text_runs``)."""

    if len(runs) == 1 and runs[0][1] == "normal":
        element.text = runs[0][0]
        return element
    return add_text_runs(element, runs)


def generate_runs_text_path(
    runs: Sequence[tuple[str, str]],
    font_size: float,
    font: str,
    *,
    dominant_baseline: str,
    text_anchor: str,
) -> Text:
    """A Text element at the origin drawing ``runs`` (``write_text_runs``)."""

    text = generate_text_path(
        "", 0, 0, 0, font_size, "normal", font, dominant_baseline=dominant_baseline, text_anchor=text_anchor
    )
    return write_text_runs(text, runs)


__all__ = ["add_text_runs", "generate_name_path", "generate_runs_text_path", "generate_text_path", "write_text_runs"]


