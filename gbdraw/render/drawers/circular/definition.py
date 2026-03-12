#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations

import math

from svgwrite.container import Group

from ....config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from ....core.text import calculate_bbox_dimensions
from ....svg.text_path import generate_name_path, generate_text_path


class DefinitionDrawer:
    """
    Draws the definition section (species/strain/accession/len/GC) on a circular canvas.
    """

    def __init__(self, config_dict: dict, cfg: GbdrawConfig | None = None) -> None:
        cfg = cfg or GbdrawConfig.from_dict(config_dict)
        self.interval = cfg.objects.definition.circular.interval
        self.font_size = cfg.objects.definition.circular.font_size
        self.font_family = cfg.objects.text.font_family

    @staticmethod
    def build_lines(
        species_parts: list[dict],
        strain_parts: list[dict],
        organelle_parts: list[dict],
        replicon_parts: list[dict],
        accession: str,
        record_length: int,
        gc_percent: float,
        *,
        show_accession: bool = True,
        show_length: bool = True,
        show_gc: bool = True,
    ) -> list[dict]:
        lines: list[dict] = []

        if species_parts and species_parts[0]["text"] is not None:
            lines.append({"kind": "name", "parts": species_parts})
        if strain_parts and strain_parts[0]["text"] is not None:
            lines.append({"kind": "name", "parts": strain_parts})
        if organelle_parts and organelle_parts[0]["text"] is not None:
            lines.append({"kind": "name", "parts": organelle_parts})
        if replicon_parts and replicon_parts[0]["text"] is not None:
            lines.append({"kind": "name", "parts": replicon_parts})
        if show_accession and accession.strip():
            lines.append({"kind": "plain", "text": accession})

        if show_length:
            lines.append({"kind": "plain", "text": "{:,} bp".format(record_length)})

        if show_gc:
            lines.append({"kind": "plain", "text": f"{gc_percent}% GC"})

        return lines

    @staticmethod
    def _line_text(line: dict) -> str:
        if line.get("kind") == "name":
            return "".join(
                str(part.get("text") or "")
                for part in list(line.get("parts") or [])
                if str(part.get("text") or "").strip()
            )
        return str(line.get("text") or "")

    def measure_max_radius(
        self,
        lines: list[dict],
        title_x: float,
        title_y: float,
        *,
        font_size: float,
        dpi: int,
    ) -> float:
        n = len(lines)
        if n == 0:
            return 0.0

        step = self.interval * 2
        start_offset = -step * (n - 2) / 2
        max_radius = 0.0

        for i, line in enumerate(lines):
            line_text = self._line_text(line)
            if not line_text:
                continue

            bbox_width, bbox_height = calculate_bbox_dimensions(
                line_text,
                self.font_family,
                float(font_size),
                int(dpi),
            )
            center_y = float(title_y) + float(start_offset + (i * step))
            half_width = 0.5 * float(bbox_width)
            half_height = 0.5 * float(bbox_height)
            corners = (
                (float(title_x) - half_width, center_y - half_height),
                (float(title_x) - half_width, center_y + half_height),
                (float(title_x) + half_width, center_y - half_height),
                (float(title_x) + half_width, center_y + half_height),
            )
            max_radius = max(max_radius, *(math.hypot(x, y) for x, y in corners))

        return float(max_radius)

    def draw(
        self,
        definition_group: Group,
        title_x: float,
        title_y: float,
        species_parts: list[dict],
        strain_parts: list[dict],
        organelle_parts: list[dict],
        replicon_parts: list[dict],
        gc_percent: float,
        accession: str,
        record_length: int,
        *,
        show_accession: bool = True,
        show_length: bool = True,
        show_gc: bool = True,
        font_size: float | None = None,
        dpi: int = 96,
        name_font_weight: str = "bold",
        lines: list[dict] | None = None,
    ) -> Group:
        active_font_size = float(font_size) if font_size is not None else float(self.font_size)
        resolved_lines = list(lines or self.build_lines(
            species_parts,
            strain_parts,
            organelle_parts,
            replicon_parts,
            accession,
            record_length,
            gc_percent,
            show_accession=show_accession,
            show_length=show_length,
            show_gc=show_gc,
        ))

        n = len(resolved_lines)
        if n == 0:
            return definition_group
        definition_group.attribs["data-definition-max-radius"] = (
            f"{self.measure_max_radius(resolved_lines, title_x, title_y, font_size=active_font_size, dpi=dpi):.6f}"
        )
        step = self.interval * 2
        start_offset = -step * (n - 2) / 2

        for i, line in enumerate(resolved_lines):
            y_shift = start_offset + i * step

            if line["kind"] == "name":
                name_path = generate_name_path(
                    line["parts"],
                    title_x,
                    title_y,
                    y_shift,
                    active_font_size,
                    name_font_weight,
                    self.font_family,
                )
                definition_group.add(name_path)
            else:
                text_path = generate_text_path(
                    line["text"],
                    title_x,
                    title_y,
                    y_shift,
                    active_font_size,
                    "normal",
                    self.font_family,
                )
                definition_group.add(text_path)

        return definition_group


__all__ = ["DefinitionDrawer"]


