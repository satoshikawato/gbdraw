#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING, Dict, List, Literal, NamedTuple, Optional, cast

from Bio.SeqRecord import SeqRecord
from gbdraw.svg.elements import Group

from ....analysis.gc import calculate_gc_percent
from ....canvas import CircularCanvasGeometry
from ....config.models import GbdrawConfig
from ....core.record_metadata import RecordSourceMetadata, infer_record_source_metadata, organism_name_lines
from ....core.text import (
    measure_mixed_text_line,
    parse_mixed_content_text,
)
from ....layout.spatial import Aabb, union_aabbs
from ....render_stages import render_stage
from ....svg.ids import definition_group_svg_id
from ...drawers.circular.definition import DefinitionDrawer



if TYPE_CHECKING:
    from gbdraw.api.requests import RecordPresentation

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


@dataclass(frozen=True)
class CircularDefinitionStyle:
    """What a Circular definition draws with: its line interval, font size, font family and profile.

    The shared-common profile (the plot title) draws in the plot title font
    size, the others in the definition font size. A render stage keys the
    definition by these values, so a setting it does not draw with does not
    rebuild it (OV-1115).
    """

    interval: int
    font_size: float
    font_family: str
    profile: CircularDefinitionProfile

    @classmethod
    def of(cls, cfg: GbdrawConfig, definition_profile: str) -> CircularDefinitionStyle:
        definition = cfg.objects.definition.circular
        profile = _normalize_definition_profile(definition_profile)
        shared_common = profile == "shared_common"
        return cls(
            interval=definition.interval,
            font_size=definition.plot_title_font_size if shared_common else definition.font_size,
            font_family=cfg.objects.text.font_family,
            profile=profile,
        )


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


class DefinitionRecordFacts(NamedTuple):
    """What a circular definition reads of its record.

    A plot title shows neither the length nor the GC content: it takes these
    facts without the GC (``gc_percent`` None) instead of the record, so a
    render stage that builds a title holds no record.
    """

    record_id: str
    metadata: RecordSourceMetadata
    length: int
    gc_percent: float | None


def _record_gc_percent(record: SeqRecord) -> float:
    return calculate_gc_percent(record.seq)


def definition_record_facts(record: SeqRecord, *, gc_percent: bool = True) -> DefinitionRecordFacts:
    """The facts a circular definition reads of ``record`` (the GC only when asked).

    The GC percent reads only the record, so it is a render stage keyed by the
    record's identity: a warm Generate does not count the bases again.
    """

    return DefinitionRecordFacts(
        record_id=str(record.id),
        metadata=infer_record_source_metadata(record),
        length=len(record.seq),
        gc_percent=render_stage("record_gc_percent", _record_gc_percent, record) if gc_percent else None,
    )


class DefinitionGroup:
    """
    Responsible for creating and managing a group for displaying genomic definition information on a circular canvas.
    """

    def __init__(
        self,
        gb_record: SeqRecord | DefinitionRecordFacts,
        geometry: CircularCanvasGeometry,
        *,
        style: CircularDefinitionStyle,
        species: Optional[str] = None,
        strain: Optional[str] = None,
        plot_title: Optional[str] = None,
        definition_group_id: str | None = None,
        record_index: int = 0,
        record_count: int = 1,
        record_presentation: RecordPresentation | None = None,
    ) -> None:
        self.record_presentation = record_presentation
        self.record_facts: DefinitionRecordFacts = (
            gb_record
            if isinstance(gb_record, DefinitionRecordFacts)
            else definition_record_facts(gb_record)
        )
        self.geometry = geometry
        self.species: str | None = species
        self.strain: str | None = strain
        self.plot_title: str | None = plot_title
        self.definition_profile: CircularDefinitionProfile = style.profile
        self.record_index = int(record_index)
        self.record_count = int(record_count)
        self.replicon: str | None = None
        self.organelle: str | None = None
        self.record_name: str = ""
        self.style = style
        self.font_size = style.font_size
        self.font = style.font_family
        self.track_id: str = self.record_facts.record_id
        self.definition_group_id: str = (
            str(definition_group_id)
            if definition_group_id
            else definition_group_svg_id(
                self.track_id,
                mode="circular",
                record_index=self.record_index,
                record_count=self.record_count,
            )
        )
        self.definition_group: Group = Group(id=self.definition_group_id)
        if self.definition_profile == "shared_common" or self.definition_group_id == "plot_title":
            self.definition_group.attribs["data-gbdraw-role"] = "plot-title"
        else:
            self.definition_group.attribs["data-gbdraw-role"] = "record-definition"
            self.definition_group.attribs["data-gbdraw-definition-part"] = "main"
            self.definition_group.attribs["data-gbdraw-record-id"] = self.record_facts.record_id
            self.definition_group.attribs["data-gbdraw-record-index"] = str(self.record_index)
        self.radius: float = geometry.radius
        self.calculate_coordinates()
        self.find_organism_name()
        self.add_circular_definitions()
        # Painted text bounds in group-local coordinates: one box per text line, and their union.
        self.line_bounds = self._measure_line_bounds()
        self.local_bounds = union_aabbs(self.line_bounds) or Aabb(0.0, 0.0, 0.0, 0.0)

    def calculate_coordinates(self) -> None:
        self.end_x_1: float = (self.radius) * math.cos(math.radians(360.0 * 0 - 90))
        self.end_x_2: float = (self.radius) * math.cos(math.radians(360.0 * (0.5) - 90))
        self.end_y_1: float = (self.radius) * math.cos(math.radians(360.0 * (0.25) - 90))
        self.end_y_2: float = (self.radius) * math.cos(math.radians(360.0 * (0.75) - 90))
        self.title_x: float = (self.end_x_1 + self.end_x_2) / 2
        self.title_y: float = (self.end_y_1 + self.end_y_2) / 2

    def find_organism_name(self) -> None:
        metadata = self.record_facts.metadata
        presentation = self.record_presentation
        explicit_label = str((presentation.label if presentation else None) or "").strip()
        explicit_subtitle = str((presentation.subtitle if presentation else None) or "").strip()
        self.replicon = metadata.replicon
        self.organelle = metadata.organelle
        # One owner chooses the two name lines from the record and the typed
        # text (the request's label and subtitle win over Species and Strain;
        # only a typed Species is cut by the Auto rule).
        lines = organism_name_lines(
            metadata,
            typed_taxon=self.species or "",
            typed_strain=explicit_subtitle or self.strain or "",
            record_label=explicit_label,
        )
        self.species_parts: List[Dict[str, str | bool | None]] = parse_mixed_content_text(lines.taxon)
        self.strain_parts: List[Dict[str, str | bool | None]] = parse_mixed_content_text(lines.strain)
        self.record_name = explicit_label or str(metadata.organism).strip() or self.record_facts.record_id

        if self.replicon:
            self.replicon_parts: List[Dict[str, str | bool | None]] = parse_mixed_content_text(self.replicon)
        else:
            self.replicon_parts = parse_mixed_content_text("")

        if self.organelle:
            self.organelle_parts: List[Dict[str, str | bool | None]] = parse_mixed_content_text(self.organelle)
        else:
            self.organelle_parts = parse_mixed_content_text("")

    def add_circular_definitions(self) -> None:
        record_length: int = self.record_facts.length
        accession: str = self.record_facts.record_id
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
        elif self.definition_profile == "shared_common":
            if isinstance(self.plot_title, str) and self.plot_title.strip():
                species_parts = parse_mixed_content_text(self.plot_title)
            else:
                species_parts = _merge_name_parts_with_single_space(species_parts, strain_parts)
            strain_parts = parse_mixed_content_text("")
            organelle_parts = parse_mixed_content_text("")
            replicon_parts = parse_mixed_content_text("")
            show_accession = False
            show_length = False
            show_gc = False
        gc_percent = self.record_facts.gc_percent
        if show_gc and gc_percent is None:
            raise ValueError("A circular definition that shows the GC content needs the record's GC.")
        active_font_size = self.font_size
        active_name_font_weight = "normal" if self.definition_profile == "shared_common" else "bold"
        self.definition_group = DefinitionDrawer(
            interval=self.style.interval, font_size=self.style.font_size, font_family=self.style.font_family
        ).draw(
            self.definition_group,
            self.title_x,
            self.title_y,
            species_parts,
            strain_parts,
            organelle_parts,
            replicon_parts,
            gc_percent if gc_percent is not None else 0.0,
            accession,
            record_length,
            show_accession=show_accession,
            show_length=show_length,
            show_gc=show_gc,
            font_size=active_font_size,
            name_font_weight=active_name_font_weight,
        )

    def _measure_line_bounds(self) -> tuple[Aabb, ...]:
        """Return authoritative painted bounds of each text line in group-local coordinates."""
        line_bounds: list[Aabb] = []
        for element in getattr(self.definition_group, "elements", ()):
            attribs = getattr(element, "attribs", None)
            if not isinstance(attribs, dict):
                continue

            parts: list[tuple[str, str]] = []
            element_text = str(getattr(element, "text", "") or "")
            if element_text:
                parts.append((element_text, str(attribs.get("font-style", "normal"))))
            for child in getattr(element, "elements", ()):
                child_text = str(getattr(child, "text", "") or "")
                if not child_text:
                    continue
                child_attribs = getattr(child, "attribs", None)
                child_style = (
                    str(child_attribs.get("font-style", attribs.get("font-style", "normal")))
                    if isinstance(child_attribs, dict)
                    else str(attribs.get("font-style", "normal"))
                )
                parts.append((child_text, child_style))
            if not parts:
                continue

            x = float(attribs.get("x", 0.0))
            y = float(attribs.get("y", 0.0))
            font_family = str(attribs.get("font-family", self.font))
            font_size = float(attribs.get("font-size", self.font_size))
            font_weight = str(attribs.get("font-weight", "normal"))
            width, height = measure_mixed_text_line(
                parts,
                font_family,
                font_size,
                self.geometry.dpi,
                font_weight=font_weight,
            )
            text_anchor = str(attribs.get("text-anchor", "start")).strip().lower()
            if text_anchor == "middle":
                left = x - (0.5 * width)
                right = x + (0.5 * width)
            elif text_anchor == "end":
                left = x - width
                right = x
            else:
                left = x
                right = x + width

            dominant_baseline = str(
                attribs.get("dominant-baseline", "auto")
            ).strip().lower()
            if dominant_baseline == "middle":
                top = y - (0.5 * height)
                bottom = y + (0.5 * height)
            else:
                top = y - height
                bottom = y

            line_bounds.append(Aabb(left, top, right, bottom))
        return tuple(line_bounds)

    def get_group(self) -> Group:
        return self.definition_group


__all__ = [
    "CircularDefinitionProfile",
    "DefinitionGroup",
    "DefinitionRecordFacts",
    "definition_record_facts",
]
