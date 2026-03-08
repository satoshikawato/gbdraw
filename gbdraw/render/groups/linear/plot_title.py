#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations

from svgwrite.container import Group
from svgwrite.text import TSpan, Text

from ....config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from ....core.text import calculate_bbox_dimensions, parse_mixed_content_text


class PlotTitleGroup:
    """Build a shared plot-title group for linear diagrams."""

    def __init__(
        self,
        title: str,
        config_dict: dict,
        *,
        font_size: float = 32.0,
        cfg: GbdrawConfig | None = None,
    ) -> None:
        cfg = cfg or GbdrawConfig.from_dict(config_dict)
        self._cfg = cfg
        self.title = str(title or "").strip()
        self.font_family = cfg.objects.text.font_family
        self.font_size = float(font_size)
        self.font_weight = "normal"
        self.text_anchor = "middle"
        self.dominant_baseline = "middle"
        self.group = Group(id="plot_title")

        self._parts = parse_mixed_content_text(self.title)
        plain = "".join(str(part.get("text") or "") for part in self._parts)
        measured_text = plain if plain else " "
        self.text_bbox_width, self.text_bbox_height = calculate_bbox_dimensions(
            measured_text,
            self.font_family,
            self.font_size,
            cfg.canvas.dpi,
        )
        self._build()

    def _build(self) -> None:
        text = Text(
            "",
            insert=(0, 0),
            stroke="none",
            fill="black",
            font_size=self.font_size,
            font_weight=self.font_weight,
            font_family=self.font_family,
            text_anchor=self.text_anchor,
            dominant_baseline=self.dominant_baseline,
        )
        for part in self._parts:
            part_text = part.get("text")
            if part_text is None:
                continue
            if part.get("italic"):
                text.add(TSpan(str(part_text), font_style="italic"))
            else:
                text.add(TSpan(str(part_text)))
        self.group.add(text)

    def get_group(self) -> Group:
        return self.group


__all__ = ["PlotTitleGroup"]

