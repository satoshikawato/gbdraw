#!/usr/bin/env python
# coding: utf-8

"""External (non-embedded) label layer for circular diagrams."""

from __future__ import annotations

from typing import Optional, Dict

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]
from pandas import DataFrame  # type: ignore[reportMissingImports]
from svgwrite.container import Group  # type: ignore[reportMissingImports]
from svgwrite.shapes import Line  # type: ignore[reportMissingImports]

from ....canvas import CircularCanvasConfigurator  # type: ignore[reportMissingImports]
from ....config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from ....features.objects import FeatureObject  # type: ignore[reportMissingImports]
from ....features.colors import preprocess_color_tables  # type: ignore[reportMissingImports]
from ....features.factory import create_feature_dict  # type: ignore[reportMissingImports]
from ....labels.filtering import preprocess_label_filtering  # type: ignore[reportMissingImports]
from ....labels.placement import prepare_label_list  # type: ignore[reportMissingImports]
from ...drawers.circular.labels import LabelDrawer  # type: ignore[reportMissingImports]
from ....configurators import FeatureDrawingConfigurator  # type: ignore[reportMissingImports]


class LabelsGroup:
    """Draws non-embedded labels + leader lines as a separate SVG group."""

    def __init__(
        self,
        gb_record: SeqRecord,
        canvas_config: CircularCanvasConfigurator,
        feature_config: FeatureDrawingConfigurator,
        config_dict: dict,
        *,
        outer_arena: tuple[float, float] | None = None,
        cfg: GbdrawConfig | None = None,
    ) -> None:
        self.gb_record: SeqRecord = gb_record
        self.canvas_config: CircularCanvasConfigurator = canvas_config
        self.feature_config: FeatureDrawingConfigurator = feature_config
        self.config_dict: dict = config_dict
        self.outer_arena = outer_arena
        cfg = cfg or GbdrawConfig.from_dict(config_dict)
        self._cfg = cfg

        raw_show_labels = cfg.canvas.show_labels
        self.show_labels = (raw_show_labels != "none") if isinstance(raw_show_labels, str) else bool(raw_show_labels)

        self.resolve_overlaps = cfg.canvas.resolve_overlaps
        self.label_stroke_width = cfg.labels.stroke_width.for_length_param(self.canvas_config.length_param)
        self.label_stroke_color = cfg.labels.stroke_color.label_stroke_color
        self.label_filtering = cfg.labels.filtering.as_dict()

        self.labels_group: Group = self.setup_labels_group()

    def setup_labels_group(self) -> Group:
        group = Group(id="labels")
        if not self.show_labels:
            return group

        selected_features_set: str = self.feature_config.selected_features_set
        color_table: Optional[DataFrame] = self.feature_config.color_table
        default_colors: Optional[DataFrame] = self.feature_config.default_colors
        label_filtering = preprocess_label_filtering(self.label_filtering)
        color_table, default_colors = preprocess_color_tables(color_table, default_colors)
        feature_dict: Dict[str, FeatureObject] = create_feature_dict(
            self.gb_record,
            color_table,
            selected_features_set,
            default_colors,
            self.canvas_config.strandedness,
            self.resolve_overlaps,
            label_filtering,
        )

        record_length: int = len(self.gb_record.seq)
        label_list = prepare_label_list(
            feature_dict,
            record_length,
            self.canvas_config.radius,
            self.canvas_config.track_ratio,
            self.config_dict,
            cfg=self._cfg,
            outer_arena=self.outer_arena,
        )

        drawer = LabelDrawer(self.config_dict, cfg=self._cfg)
        for label in label_list:
            if label.get("is_embedded"):
                continue
            # Leader lines first (so they appear behind features)
            line_path = Line(
                start=(label["middle_x"], label["middle_y"]),
                end=(label["start_x"], label["start_y"]),
                stroke=self.label_stroke_color,
                stroke_width=self.label_stroke_width,
            )
            group.add(line_path)
            line_path2 = Line(
                start=(label["middle_x"], label["middle_y"]),
                end=(label["feature_middle_x"], label["feature_middle_y"]),
                stroke=self.label_stroke_color,
                stroke_width=self.label_stroke_width,
            )
            group.add(line_path2)
            # Label text after lines (so it appears on top)
            group = drawer.draw(label, group, record_length, self.canvas_config.radius, self.canvas_config.track_ratio)

        return group

    def get_group(self) -> Group:
        return self.labels_group


__all__ = ["LabelsGroup"]


