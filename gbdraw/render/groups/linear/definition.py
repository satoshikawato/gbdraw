#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations

from dataclasses import dataclass

from Bio.SeqRecord import SeqRecord
from svgwrite.container import Group
from svgwrite.text import Text, TSpan

from ....config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from ....core.record_metadata import infer_record_source_metadata
from ....core.text import calculate_bbox_dimensions, create_text_element, parse_mixed_content_text

_COORD_BASE_KEY = "gbdraw_coord_base"
_COORD_STEP_KEY = "gbdraw_coord_step"
_REGION_APPLIED_KEY = "gbdraw_region_applied"


@dataclass
class _DefinitionLine:
    kind: str
    text: str
    width: float
    height: float
    parts: list[dict[str, str | bool | None]] | None = None
    y: float = 0.0


class DefinitionGroup:
    """Build the stacked linear definition block for one record."""

    def __init__(
        self,
        record: SeqRecord,
        config_dict: dict,
        canvas_config: dict,
        title_start_x: float = 0,
        title_start_y: float = -9,
        length_start_x: float = 0,
        length_start_y: float = 9,
        cfg: GbdrawConfig | None = None,
    ) -> None:
        self.record: SeqRecord = record
        self.title_start_x: float = title_start_x
        self.title_start_y: float = title_start_y
        self.length_start_x: float = length_start_x
        self.length_start_y: float = length_start_y
        self.definition_bounding_box_width: float = 0.0
        self.definition_bounding_box_height: float = 0.0
        self.canvas_config = canvas_config
        self.length_param = self.canvas_config.length_param
        cfg = cfg or GbdrawConfig.from_dict(config_dict)
        self._cfg = cfg
        def_cfg = cfg.objects.definition.linear
        self.linear_definition_stroke: float = def_cfg.stroke
        self.linear_definition_fill: str = def_cfg.fill
        self.linear_definition_font_size: float = def_cfg.font_size.for_length_param(self.length_param)
        self.linear_definition_font_weight: str = def_cfg.font_weight
        self.linear_definition_font_family: str = cfg.objects.text.font_family
        self.interval: int = def_cfg.interval
        self.show_replicon: bool = bool(def_cfg.show_replicon)
        self.show_accession: bool = bool(def_cfg.show_accession)
        self.show_length: bool = bool(def_cfg.show_length)
        self.dpi: int = cfg.canvas.dpi
        self.linear_text_anchor: str = def_cfg.text_anchor
        self.linear_dominant_baseline: str = def_cfg.dominant_baseline
        self.get_definition_details()
        self.definition_lines = self._build_definition_lines()
        self.calculate_start_coordinates()
        self.definition_group = Group(id=self.definition_group_id)
        self.add_elements_to_group()

    def calculate_start_coordinates(self) -> None:
        if not self.definition_lines:
            self.definition_bounding_box_width = 0.0
            self.definition_bounding_box_height = 0.0
            return

        max_width = max(line.width for line in self.definition_lines)
        total_height = sum(line.height for line in self.definition_lines)
        if len(self.definition_lines) > 1:
            total_height += self.interval * (len(self.definition_lines) - 1)

        self.definition_bounding_box_width = max_width
        self.definition_bounding_box_height = total_height

        current_top = -(total_height / 2.0)
        for line in self.definition_lines:
            line.y = current_top + (line.height / 2.0)
            current_top += line.height + self.interval

    def get_definition_details(self) -> None:
        """Resolve all labels that may participate in the stacked definition."""
        self.track_id = str(self.record.id)
        self.definition_group_id = f"{self.track_id.replace(' ', '_')}_definition"

        override = None
        if getattr(self.record, "annotations", None):
            override = self.record.annotations.get("gbdraw_record_label")
        if override is not None:
            override = str(override).strip()

        self.record_name: str = override if override else ""
        self.record_name_parts = parse_mixed_content_text(self.record_name)
        self.record_name_plain = "".join(part.get("text") or "" for part in self.record_name_parts).strip()

        metadata = infer_record_source_metadata(self.record)
        self.replicon_label = str(metadata.replicon or "").strip()
        self.accession_label = self.track_id

        self.record_length: int = len(self.record.seq)
        if self._is_region_applied():
            self.length_label = self._format_coordinate_span()
        else:
            self.length_label = "{:,} bp".format(self.record_length)

    def _build_definition_lines(self) -> list[_DefinitionLine]:
        lines: list[_DefinitionLine] = []

        if self.record_name_plain:
            width, height = calculate_bbox_dimensions(
                self.record_name_plain,
                self.linear_definition_font_family,
                self.linear_definition_font_size,
                self.dpi,
            )
            lines.append(
                _DefinitionLine(
                    kind="mixed",
                    text=self.record_name_plain,
                    width=width,
                    height=height,
                    parts=self.record_name_parts,
                )
            )

        if self.show_replicon and self.replicon_label:
            width, height = calculate_bbox_dimensions(
                self.replicon_label,
                self.linear_definition_font_family,
                self.linear_definition_font_size,
                self.dpi,
            )
            lines.append(
                _DefinitionLine(
                    kind="plain",
                    text=self.replicon_label,
                    width=width,
                    height=height,
                )
            )

        if self.show_accession and self.accession_label.strip():
            width, height = calculate_bbox_dimensions(
                self.accession_label,
                self.linear_definition_font_family,
                self.linear_definition_font_size,
                self.dpi,
            )
            lines.append(
                _DefinitionLine(
                    kind="plain",
                    text=self.accession_label,
                    width=width,
                    height=height,
                )
            )

        if self.show_length and self.length_label.strip():
            width, height = calculate_bbox_dimensions(
                self.length_label,
                self.linear_definition_font_family,
                self.linear_definition_font_size,
                self.dpi,
            )
            lines.append(
                _DefinitionLine(
                    kind="plain",
                    text=self.length_label,
                    width=width,
                    height=height,
                )
            )

        return lines

    def _is_region_applied(self) -> bool:
        annotations = getattr(self.record, "annotations", None) or {}
        raw = annotations.get(_REGION_APPLIED_KEY, False)
        if isinstance(raw, str):
            return raw.strip().lower() in {"1", "true", "yes", "on"}
        return bool(raw)

    def _format_coordinate_span(self) -> str:
        annotations = getattr(self.record, "annotations", None) or {}
        try:
            start_coord = int(annotations.get(_COORD_BASE_KEY, 1))
        except (TypeError, ValueError):
            start_coord = 1
        try:
            step = int(annotations.get(_COORD_STEP_KEY, 1))
        except (TypeError, ValueError):
            step = 1
        if step == 0:
            step = 1
        step = 1 if step > 0 else -1
        end_coord = start_coord + (step * max(0, self.record_length - 1))
        return f"{start_coord}-{end_coord}"

    def add_elements_to_group(self) -> None:
        """Add all visible definition lines to the SVG group."""
        for line in self.definition_lines:
            if line.kind == "mixed":
                self.definition_group.add(self._create_name_text(line.parts or [], line.y))
            else:
                self.definition_group.add(
                    create_text_element(
                        line.text,
                        0,
                        line.y,
                        self.linear_definition_font_size,
                        self.linear_definition_font_weight,
                        self.linear_definition_font_family,
                        text_anchor=self.linear_text_anchor,
                        dominant_baseline=self.linear_dominant_baseline,
                    )
                )

    def _create_name_text(self, parts: list[dict[str, str | bool | None]], y_pos: float) -> Text:
        text_el = Text(
            "",
            insert=(0, y_pos),
            stroke=self.linear_definition_stroke,
            fill=self.linear_definition_fill,
            font_size=self.linear_definition_font_size,
            font_weight=self.linear_definition_font_weight,
            font_family=self.linear_definition_font_family,
            text_anchor=self.linear_text_anchor,
            dominant_baseline=self.linear_dominant_baseline,
        )
        for part in parts:
            part_text = part.get("text")
            if part_text is None:
                continue
            if part.get("italic"):
                text_el.add(TSpan(part_text, font_style="italic"))
            else:
                text_el.add(TSpan(part_text))
        return text_el

    def get_group(self) -> Group:
        return self.definition_group


__all__ = ["DefinitionGroup"]
