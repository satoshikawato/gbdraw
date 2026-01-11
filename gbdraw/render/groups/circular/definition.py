#!/usr/bin/env python
# coding: utf-8

import math
from typing import Optional, List, Dict

from Bio.SeqRecord import SeqRecord
from svgwrite.container import Group

from ....analysis.gc import calculate_gc_percent
from ....canvas import CircularCanvasConfigurator
from ....config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from ...drawers.circular.definition import DefinitionDrawer
from ....core.text import parse_mixed_content_text


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
        cfg: GbdrawConfig | None = None,
    ) -> None:
        self.gb_record: SeqRecord = gb_record
        self.canvas_config: CircularCanvasConfigurator = canvas_config
        self.species: str | None = species
        self.strain: str | None = strain
        self.replicon: str | None = None
        self.organelle: str | None = None
        self.config_dict: dict = config_dict
        cfg = cfg or GbdrawConfig.from_dict(config_dict)
        self._cfg = cfg
        self.interval = cfg.objects.definition.circular.interval
        self.font_size = cfg.objects.definition.circular.font_size
        self.font = cfg.objects.text.font_family
        self.track_id: str = str(self.gb_record.id).replace(" ", "_")
        self.definition_group_id: str = f"{self.track_id}_definition"
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
        self.definition_group: Group = DefinitionDrawer(self.config_dict, cfg=self._cfg).draw(
            self.definition_group,
            self.title_x,
            self.title_y,
            self.species_parts,
            self.strain_parts,
            self.organelle_parts,
            self.replicon_parts,
            gc_percent,
            accession,
            record_length,
        )

    def get_group(self) -> Group:
        return self.definition_group


__all__ = ["DefinitionGroup"]


