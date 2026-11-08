#!/usr/bin/env python
# coding: utf-8

from collections.abc import Iterable, Sequence

from gbdraw.svg.elements import Text, TSpan


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


def add_text_runs(text: Text, runs: Iterable[tuple[str, str]]) -> Text:
    """Append one ``tspan`` per ``(text, font_style)`` run to ``text``; an italic run carries ``font-style``."""

    for part, font_style in runs:
        text.add(TSpan(part, font_style="italic") if font_style == "italic" else TSpan(part))
    return text


def generate_runs_text_path(
    runs: Sequence[tuple[str, str]],
    font_size: float,
    font: str,
    *,
    dominant_baseline: str,
    text_anchor: str,
) -> Text:
    """A Text element at the origin drawing ``runs``: one upright run as a plain text node (``generate_text_path``),
    styled runs as ``tspan``s (``add_text_runs``)."""

    if len(runs) == 1 and runs[0][1] == "normal":
        return generate_text_path(
            runs[0][0], 0, 0, 0, font_size, "normal", font, dominant_baseline=dominant_baseline, text_anchor=text_anchor
        )
    text = generate_text_path(
        "", 0, 0, 0, font_size, "normal", font, dominant_baseline=dominant_baseline, text_anchor=text_anchor
    )
    return add_text_runs(text, runs)


__all__ = ["add_text_runs", "generate_name_path", "generate_runs_text_path", "generate_text_path"]


