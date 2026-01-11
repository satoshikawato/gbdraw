#!/usr/bin/env python
# coding: utf-8

from svgwrite.container import Group

from ....config.models import GbdrawConfig  # type: ignore[reportMissingImports]
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

    def draw(
        self,
        definition_group: Group,
        title_x: float,
        title_y: float,
        species_parts: list,
        strain_parts: list,
        organelle_parts: list,
        repicon_parts: list,
        gc_percent: float,
        accession: str,
        record_length: int,
    ) -> Group:
        lines: list[dict] = []

        if species_parts and species_parts[0]["text"] is not None:
            lines.append({"kind": "name", "parts": species_parts})
        if strain_parts and strain_parts[0]["text"] is not None:
            lines.append({"kind": "name", "parts": strain_parts})
        if organelle_parts and organelle_parts[0]["text"] is not None:
            lines.append({"kind": "name", "parts": organelle_parts})
        if repicon_parts and repicon_parts[0]["text"] is not None:
            lines.append({"kind": "name", "parts": repicon_parts})
        if accession.strip():
            lines.append({"kind": "plain", "text": accession})

        length_text = "{:,} bp".format(record_length)
        lines.append({"kind": "plain", "text": length_text})

        gc_text = f"{gc_percent}% GC"
        lines.append({"kind": "plain", "text": gc_text})

        n = len(lines)
        step = self.interval * 2
        start_offset = -step * (n - 2) / 2

        for i, line in enumerate(lines):
            y_shift = start_offset + i * step

            if line["kind"] == "name":
                name_path = generate_name_path(
                    line["parts"],
                    title_x,
                    title_y,
                    y_shift,
                    self.font_size,
                    "bold",
                    self.font_family,
                )
                definition_group.add(name_path)
            else:
                text_path = generate_text_path(
                    line["text"],
                    title_x,
                    title_y,
                    y_shift,
                    self.font_size,
                    "normal",
                    self.font_family,
                )
                definition_group.add(text_path)

        return definition_group


__all__ = ["DefinitionDrawer"]


