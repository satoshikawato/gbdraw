#!/usr/bin/env python
# coding: utf-8

from Bio.SeqRecord import SeqRecord
from svgwrite.container import Group
from svgwrite.text import Text, TSpan

from ....core.text import create_text_element, calculate_bbox_dimensions, parse_mixed_content_text
from ....config.models import GbdrawConfig  # type: ignore[reportMissingImports]


class DefinitionGroup:
    """
    Handles the creation and display of definition details for a SeqRecord in a linear layout.

    Attributes:
        record (SeqRecord): GenBank record containing genomic data.
        title_start_x (float): Horizontal start position for the title text.
        title_start_y (float): Vertical start position for the title text.
        length_start_x (float): Horizontal start position for the length text.
        length_start_y (float): Vertical start position for the length text.
        linear_definition_stroke (float): Stroke width for the text.
        linear_definition_fill (str): Fill color for the text.
        linear_definition_font_size (str): Font size for the text.
        linear_definition_font_weight (str): Font weight for the text.
        linear_definition_font_family (str): Font family for the text.
        linear_text_anchor (str): Text anchor position.
        linear_dominant_baseline (str): Dominant baseline position for the text.
    """

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
        """
        Initializes the DefinitionGroup with the given SeqRecord and configuration.

        Args:
            record (SeqRecord): The genomic record to be used for the definition details.
            config_dict (dict): Configuration dictionary with styling parameters.
            title_start_x (float): Horizontal start position for the title text.
            title_start_y (float): Vertical start position for the title text.
            length_start_x (float): Horizontal start position for the length text.
            length_start_y (float): Vertical start position for the length text.
        """
        self.record: SeqRecord = record
        self.title_start_x: float = title_start_x
        self.title_start_y: float = title_start_y
        self.length_start_x: float = length_start_x
        self.length_start_y: float = length_start_y
        self.definition_bounding_box_width: float = 0
        self.definition_bounding_box_height: float = 0
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
        self.dpi: int = cfg.canvas.dpi
        self.linear_text_anchor: str = def_cfg.text_anchor
        self.linear_dominant_baseline: str = def_cfg.dominant_baseline
        self.get_id_and_length()
        self.name_bounding_box_width, self.name_bounding_box_height = calculate_bbox_dimensions(
            self.record_name_plain, self.linear_definition_font_family, self.linear_definition_font_size, self.dpi
        )
        self.length_bounding_box_width, self.length_bounding_box_height = calculate_bbox_dimensions(
            self.length_label, self.linear_definition_font_family, self.linear_definition_font_size, self.dpi
        )
        self.calculate_start_coordinates()
        self.definition_group = Group(id=self.track_id)
        self.add_elements_to_group()

    def calculate_start_coordinates(self) -> None:
        max_width = max(self.name_bounding_box_width, self.length_bounding_box_width)
        total_height = self.name_bounding_box_height + self.length_bounding_box_height + self.interval
        self.definition_bounding_box_width = max_width
        self.definition_bounding_box_height = total_height
        self.title_start_x = 0
        self.title_start_y = -(total_height / 2) + (self.name_bounding_box_height / 2)
        self.length_start_x = 0
        self.length_start_y = (total_height / 2) - (self.length_bounding_box_height / 2)

    def get_id_and_length(self) -> None:
        """
        Extracts the ID and length of the SeqRecord for display purposes.
        """
        self.track_id = str(self.record.id)
        override = None
        if getattr(self.record, "annotations", None):
            override = self.record.annotations.get("gbdraw_record_label")
        if override is not None:
            override = str(override).strip()
        self.record_name: str = override if override else self.track_id
        parts = parse_mixed_content_text(self.record_name)
        self.record_name_plain = "".join(part.get("text") or "" for part in parts)
        self.record_length: int = len(self.record.seq)
        self.length_label: str = "{:,} bp".format(self.record_length)

    def add_elements_to_group(self) -> None:
        """
        Adds the definition elements (like record name and length) to the group.
        """
        self.name_path: Text = self._create_name_text()
        self.length_path: Text = create_text_element(
            self.length_label,
            self.length_start_x,
            self.length_start_y,
            self.linear_definition_font_size,
            self.linear_definition_font_weight,
            self.linear_definition_font_family,
            text_anchor=self.linear_text_anchor,
            dominant_baseline=self.linear_dominant_baseline,
        )
        self.definition_group.add(self.name_path)
        self.definition_group.add(self.length_path)

    def _create_name_text(self) -> Text:
        parts = parse_mixed_content_text(self.record_name)
        text_el = Text(
            "",
            insert=(self.title_start_x, self.title_start_y),
            stroke="none",
            fill="black",
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
        """
        Retrieves the SVG group containing the definition details.

        Returns:
            Group: The SVG group with definition details.
        """
        return self.definition_group


__all__ = ["DefinitionGroup"]


