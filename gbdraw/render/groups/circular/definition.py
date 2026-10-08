#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations

import math
import re
from functools import lru_cache
from typing import Callable, Dict, List, Literal, Optional, Sequence, cast

from Bio.SeqRecord import SeqRecord
from svgwrite.container import Group

from ....analysis.gc import calculate_gc_percent
from ....canvas import CircularCanvasConfigurator
from ....config.models import GbdrawConfig
from ....core.record_metadata import infer_record_source_metadata
from ....core.numeric import compensated_sum
from ....core.text import (
    calculate_bbox_dimensions,
    calculate_svg_bbox_dimensions,
    parse_mixed_content_text,
)
from ....layout.spatial import Aabb
from ....svg.ids import definition_group_svg_id
from ...drawers.circular.definition import DefinitionDrawer


from gbdraw.layout.record_coordinates import RecordDisplayTransform

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


_WHITESPACE_RUN = re.compile(r"(\s+)")


def _merge_adjacent_parts(parts: Sequence[_TextPart]) -> list[_TextPart]:
    merged: list[_TextPart] = []
    for part in parts:
        if merged and bool(merged[-1]["italic"]) == bool(part["italic"]):
            merged[-1] = {"text": f"{merged[-1]['text']}{part['text']}", "italic": bool(part["italic"])}
        else:
            merged.append({"text": str(part["text"]), "italic": bool(part["italic"])})
    return merged


def _split_words(parts: Sequence[_TextPart]) -> tuple[list[list[_TextPart]], list[_TextPart]]:
    """Split mixed italic/plain parts into words and the whitespace between them."""

    words: list[list[_TextPart]] = []
    separators: list[_TextPart] = []
    current: list[_TextPart] = []
    pending_separator: _TextPart | None = None
    for part in parts:
        text = part.get("text")
        if not isinstance(text, str) or not text:
            continue
        italic = bool(part.get("italic"))
        for piece in _WHITESPACE_RUN.split(text):
            if not piece:
                continue
            if piece.isspace():
                if current:
                    words.append(current)
                    current = []
                    pending_separator = {"text": piece, "italic": italic}
                elif pending_separator is not None:
                    pending_separator = {
                        "text": f"{pending_separator['text']}{piece}",
                        "italic": bool(pending_separator["italic"]),
                    }
                continue
            if not current and words:
                separators.append(pending_separator or {"text": " ", "italic": False})
                pending_separator = None
            current.append({"text": piece, "italic": italic})
    if current:
        words.append(current)
    return words, separators


def wrap_definition_line_parts(
    parts: Sequence[_TextPart],
    line_count: int,
    measure_width: Callable[[str], float],
) -> list[list[_TextPart]]:
    """Wrap one definition line at word boundaries into balanced lines.

    The text, its italic runs and the whitespace inside each line are kept; the
    whitespace at a break is consumed. The partition minimizes the widest line.
    With fewer words than ``line_count`` every word gets its own line.
    """

    words, separators = _split_words(parts)
    count = min(max(1, int(line_count)), len(words))

    def line_parts(start: int, end: int) -> list[_TextPart]:
        segments: list[_TextPart] = []
        for index in range(start, end):
            if index > start:
                segments.append(separators[index - 1])
            segments.extend(words[index])
        return _merge_adjacent_parts(segments)

    if count <= 1:
        return [line_parts(0, len(words))] if words else []

    @lru_cache(maxsize=None)
    def width(start: int, end: int) -> float:
        return float(measure_width("".join(str(part["text"]) for part in line_parts(start, end))))

    @lru_cache(maxsize=None)
    def best(start: int, lines: int) -> tuple[float, tuple[int, ...]]:
        if lines == 1:
            return width(start, len(words)), ()
        chosen: tuple[float, tuple[int, ...]] | None = None
        for end in range(start + 1, len(words) - lines + 2):
            tail_width, tail_breaks = best(end, lines - 1)
            candidate = (max(width(start, end), tail_width), (end, *tail_breaks))
            if chosen is None or candidate[0] < chosen[0] - 1e-9:
                chosen = candidate
        assert chosen is not None
        return chosen

    _widest, breaks = best(0, count)
    bounds = (0, *breaks, len(words))
    return [line_parts(start, end) for start, end in zip(bounds, bounds[1:])]


class DefinitionGroup:
    """
    Responsible for creating and managing a group for displaying genomic definition information on a circular canvas.
    """

    def __init__(
        self,
        gb_record: SeqRecord,
        canvas_config: CircularCanvasConfigurator,
        *,
        cfg: GbdrawConfig,
        species: Optional[str] = None,
        strain: Optional[str] = None,
        plot_title: Optional[str] = None,
        definition_profile: CircularDefinitionProfile | str = "full",
        definition_group_id: str | None = None,
        record_index: int = 0,
        record_count: int = 1,
        record_transform: RecordDisplayTransform | None = None,
        species_line_count: int = 1,
    ) -> None:
        self.record_transform = record_transform
        self.species_line_count = max(1, int(species_line_count))
        self.gb_record: SeqRecord = gb_record
        self.canvas_config: CircularCanvasConfigurator = canvas_config
        self.species: str | None = species
        self.strain: str | None = strain
        self.plot_title: str | None = plot_title
        self.definition_profile: CircularDefinitionProfile = _normalize_definition_profile(
            str(definition_profile)
        )
        self.record_index = int(record_index)
        self.record_count = int(record_count)
        self.replicon: str | None = None
        self.organelle: str | None = None
        self.record_name: str = ""
        self._cfg = cfg
        self.interval = cfg.objects.definition.circular.interval
        self.font_size = cfg.objects.definition.circular.font_size
        self.plot_title_font_size = cfg.objects.definition.circular.plot_title_font_size
        self.font = cfg.objects.text.font_family
        self.track_id: str = str(self.gb_record.id)
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
        self.definition_group: Group = Group(id=self.definition_group_id, debug=False)
        if self.definition_profile == "shared_common" or self.definition_group_id == "plot_title":
            self.definition_group.attribs["data-gbdraw-role"] = "plot-title"
        else:
            self.definition_group.attribs["data-gbdraw-role"] = "record-definition"
            self.definition_group.attribs["data-gbdraw-definition-part"] = "main"
            self.definition_group.attribs["data-gbdraw-record-id"] = str(self.gb_record.id)
            self.definition_group.attribs["data-gbdraw-record-index"] = str(self.record_index)
        self.radius: float = self.canvas_config.radius
        self.calculate_coordinates()
        self.find_organism_name()
        self.add_circular_definitions()
        self.local_bounds = self._measure_local_bounds()

    def calculate_coordinates(self) -> None:
        self.end_x_1: float = (self.radius) * math.cos(math.radians(360.0 * 0 - 90))
        self.end_x_2: float = (self.radius) * math.cos(math.radians(360.0 * (0.5) - 90))
        self.end_y_1: float = (self.radius) * math.cos(math.radians(360.0 * (0.25) - 90))
        self.end_y_2: float = (self.radius) * math.cos(math.radians(360.0 * (0.75) - 90))
        self.title_x: float = (self.end_x_1 + self.end_x_2) / 2
        self.title_y: float = (self.end_y_1 + self.end_y_2) / 2

    def find_organism_name(self) -> None:
        metadata = infer_record_source_metadata(self.gb_record)
        strain_name = metadata.strain
        record_name = metadata.organism
        annotations = getattr(self.gb_record, "annotations", {}) or {}
        explicit_label = str(annotations.get("gbdraw_record_label") or "").strip()
        explicit_subtitle = str(
            annotations.get("gbdraw_record_subtitle") or ""
        ).strip()
        self.replicon = metadata.replicon
        self.organelle = metadata.organelle

        if explicit_label:
            self.species_parts: List[Dict[str, str | bool | None]] = parse_mixed_content_text(explicit_label)
        elif self.species:
            self.species_parts = parse_mixed_content_text(self.species)
        else:
            self.species_parts = parse_mixed_content_text(record_name)
        self.record_name = explicit_label or str(record_name).strip() or str(self.gb_record.id)

        if explicit_subtitle:
            self.strain_parts: List[Dict[str, str | bool | None]] = parse_mixed_content_text(explicit_subtitle)
        elif self.strain:
            self.strain_parts = parse_mixed_content_text(self.strain)
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
        active_font_size = (
            self.plot_title_font_size if self.definition_profile == "shared_common" else self.font_size
        )
        active_name_font_weight = "normal" if self.definition_profile == "shared_common" else "bold"
        species_line_parts: list[list[_TextPart]] | None = None
        if self.definition_profile == "full" and self.species_line_count > 1:
            wrapped = wrap_definition_line_parts(
                cast(list[_TextPart], species_parts),
                self.species_line_count,
                lambda text: calculate_bbox_dimensions(
                    text,
                    self.font,
                    active_font_size,
                    int(self.canvas_config.dpi),
                )[0],
            )
            if len(wrapped) > 1:
                species_line_parts = wrapped

        self.definition_group = DefinitionDrawer(cfg=self._cfg).draw(
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
            record_transform=self.record_transform,
            species_line_parts=species_line_parts,
        )

    def _measure_local_bounds(self) -> Aabb:
        """Return authoritative painted text bounds in group-local coordinates."""
        min_x: float | None = None
        min_y: float | None = None
        max_x: float | None = None
        max_y: float | None = None

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
            widths: list[float] = []
            heights: list[float] = []
            for text, font_style in parts:
                width, height = calculate_svg_bbox_dimensions(
                    text,
                    font_family,
                    font_size,
                    self._cfg.canvas.dpi,
                    font_weight=font_weight,
                    font_style=font_style,
                )
                widths.append(float(width))
                heights.append(float(height))

            width = compensated_sum(widths)
            height = max(heights, default=font_size)
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

            min_x = left if min_x is None else min(min_x, left)
            min_y = top if min_y is None else min(min_y, top)
            max_x = right if max_x is None else max(max_x, right)
            max_y = bottom if max_y is None else max(max_y, bottom)

        if min_x is None or min_y is None or max_x is None or max_y is None:
            return Aabb(0.0, 0.0, 0.0, 0.0)
        return Aabb(min_x, min_y, max_x, max_y)

    def get_group(self) -> Group:
        return self.definition_group


__all__ = ["CircularDefinitionProfile", "DefinitionGroup", "wrap_definition_line_parts"]
