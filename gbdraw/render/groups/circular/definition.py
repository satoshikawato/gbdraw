#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations

import math
from typing import Dict, List, Literal, Optional, cast

from Bio.SeqRecord import SeqRecord
from svgwrite.container import Group

from ....analysis.gc import calculate_gc_percent
from ....canvas import CircularCanvasConfigurator
from ....config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from ...drawers.circular.definition import DefinitionDrawer
from ....core.text import parse_mixed_content_text


CircularDefinitionProfile = Literal["full", "record_summary", "shared_common"]
_SUPPORTED_DEFINITION_PROFILES = {"full", "record_summary", "shared_common"}
_TextPart = Dict[str, str | bool | None]


def _normalize_definition_profile(profile: str) -> CircularDefinitionProfile:
    normalized = str(profile).strip().lower()
    if normalized not in _SUPPORTED_DEFINITION_PROFILES:
        raise ValueError(
            "definition_profile must be one of: full, record_summary, shared_common"
        )
    return cast(CircularDefinitionProfile, normalized)


def _has_visible_text(parts: list[_TextPart]) -> bool:
    return any(
        isinstance(part.get("text"), str) and str(part.get("text")).strip()
        for part in parts
    )


def _clone_nonempty_parts(parts: list[_TextPart]) -> list[_TextPart]:
    cloned: list[_TextPart] = []
    for part in parts:
        text = part.get("text")
        if not isinstance(text, str) or not text:
            continue
        cloned.append({"text": text, "italic": bool(part.get("italic"))})
    return cloned


def _merge_name_parts_with_single_space(
    species_parts: list[_TextPart],
    strain_parts: list[_TextPart],
) -> list[_TextPart]:
    merged: list[_TextPart] = []
    species_clean = _clone_nonempty_parts(species_parts)
    strain_clean = _clone_nonempty_parts(strain_parts)

    if _has_visible_text(species_clean):
        merged.extend(species_clean)
    if _has_visible_text(strain_clean):
        if _has_visible_text(merged):
            first_text = strain_clean[0].get("text")
            if isinstance(first_text, str):
                strain_clean[0]["text"] = f" {first_text.lstrip()}"
        merged.extend(strain_clean)

    if _has_visible_text(merged):
        return merged
    return cast(list[_TextPart], parse_mixed_content_text(""))


class DefinitionGroup:
    """
    Responsible for creating and managing a group for displaying genomic definition information on a circular canvas.
    """

    def __init__(
        self,
        gb_record: SeqRecord,
        canvas_config: CircularCanvasConfigurator,
        config_dict: dict,
        species: Optional[str] = None,
        strain: Optional[str] = None,
        definition_profile: CircularDefinitionProfile | str = "full",
        definition_group_id: str | None = None,
        cfg: GbdrawConfig | None = None,
    ) -> None:
        self.gb_record: SeqRecord = gb_record
        self.canvas_config: CircularCanvasConfigurator = canvas_config
        self.species: str | None = species
        self.strain: str | None = strain
        self.definition_profile: CircularDefinitionProfile = _normalize_definition_profile(
            str(definition_profile)
        )
        self.replicon: str | None = None
        self.organelle: str | None = None
        self.record_name: str = ""
        self.config_dict: dict = config_dict
        cfg = cfg or GbdrawConfig.from_dict(config_dict)
        self._cfg = cfg
        self.interval = cfg.objects.definition.circular.interval
        self.font_size = cfg.objects.definition.circular.font_size
        self.shared_font_size = cfg.objects.definition.circular.shared_font_size
        self.font = cfg.objects.text.font_family
        self.track_id: str = str(self.gb_record.id).replace(" ", "_")
        self.definition_group_id: str = (
            str(definition_group_id) if definition_group_id else f"{self.track_id}_definition"
        )
        self.definition_group = Group(id=self.definition_group_id)
        self.radius: float = self.canvas_config.radius
        self.calculate_coordinates()
        self.find_organism_name()
        self.add_circular_definitions()

    def calculate_coordinates(self) -> None:
        self.end_x_1: float = (self.radius) * math.cos(math.radians(360.0 * 0 - 90))
        self.end_x_2: float = (self.radius) * math.cos(math.radians(360.0 * (0.5) - 90))
        self.end_y_1: float = (self.radius) * math.cos(math.radians(360.0 * (0.25) - 90))
        self.end_y_2: float = (self.radius) * math.cos(math.radians(360.0 * (0.75) - 90))
        self.title_x: float = (self.end_x_1 + self.end_x_2) / 2
        self.title_y: float = (self.end_y_1 + self.end_y_2) / 2

    def find_organism_name(self) -> None:
        strain_name: str = ""
        record_name: str = ""

        for feature in self.gb_record.features:
            if feature.type == "source":
                if "organism" in feature.qualifiers.keys():
                    record_name = feature.qualifiers["organism"][0]
                elif "organism" in self.gb_record.annotations.keys():
                    record_name = str(self.gb_record.annotations["organism"])
                if "isolate" in feature.qualifiers.keys():
                    strain_name = feature.qualifiers["isolate"][0]
                elif "strain" in feature.qualifiers.keys():
                    strain_name = feature.qualifiers["strain"][0]
                if "chromosome" in feature.qualifiers.keys():
                    self.replicon = f"Chromosome {feature.qualifiers['chromosome'][0]}"
                elif "plasmid" in feature.qualifiers.keys():
                    self.replicon = feature.qualifiers["plasmid"][0]
                if "organelle" in feature.qualifiers.keys():
                    self.organelle = feature.qualifiers["organelle"][0]

        if self.species:
            self.species_parts: List[Dict[str, str | bool | None]] = parse_mixed_content_text(self.species)
        else:
            self.species_parts = parse_mixed_content_text(record_name)
        self.record_name = str(record_name).strip() or str(self.gb_record.id)

        if self.strain:
            self.strain_parts: List[Dict[str, str | bool | None]] = parse_mixed_content_text(self.strain)
        else:
            self.strain_parts = parse_mixed_content_text(strain_name)

        if self.replicon:
            self.replicon_parts: List[Dict[str, str | bool | None]] = parse_mixed_content_text(self.replicon)
        else:
            self.replicon_parts = parse_mixed_content_text("")

        if self.organelle:
            self.organelle_parts: List[Dict[str, str | bool | None]] = parse_mixed_content_text(self.organelle)
        else:
            self.organelle_parts = parse_mixed_content_text("")

    def add_circular_definitions(self) -> None:
        record_length: int = len(self.gb_record.seq)
        accession: str = self.gb_record.id
        gc_percent: float = calculate_gc_percent(self.gb_record.seq)
        show_accession = True
        show_length = True
        show_gc = True

        species_parts = self.species_parts
        strain_parts = self.strain_parts
        organelle_parts = self.organelle_parts
        replicon_parts = self.replicon_parts

        if self.definition_profile == "record_summary":
            species_parts = parse_mixed_content_text("")
            strain_parts = parse_mixed_content_text("")
            organelle_parts = parse_mixed_content_text("")
            replicon_parts = parse_mixed_content_text("")
        elif self.definition_profile == "shared_common":
            species_parts = _merge_name_parts_with_single_space(species_parts, strain_parts)
            strain_parts = parse_mixed_content_text("")
            organelle_parts = parse_mixed_content_text("")
            replicon_parts = parse_mixed_content_text("")
            show_accession = False
            show_length = False
            show_gc = False
        active_font_size = (
            self.shared_font_size if self.definition_profile == "shared_common" else self.font_size
        )
        active_name_font_weight = "normal" if self.definition_profile == "shared_common" else "bold"

        self.definition_group: Group = DefinitionDrawer(self.config_dict, cfg=self._cfg).draw(
            self.definition_group,
            self.title_x,
            self.title_y,
            species_parts,
            strain_parts,
            organelle_parts,
            replicon_parts,
            gc_percent,
            accession,
            record_length,
            show_accession=show_accession,
            show_length=show_length,
            show_gc=show_gc,
            font_size=active_font_size,
            name_font_weight=active_name_font_weight,
        )

    def get_group(self) -> Group:
        return self.definition_group


__all__ = ["CircularDefinitionProfile", "DefinitionGroup"]


