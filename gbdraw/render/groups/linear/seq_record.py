#!/usr/bin/env python
# coding: utf-8

from typing import Optional

from pandas import DataFrame
from Bio.SeqRecord import SeqRecord
from svgwrite.container import Group
from svgwrite.shapes import Line
from svgwrite.text import Text

from ....canvas import LinearCanvasConfigurator
from ....config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from ...drawers.linear.features import FeatureDrawer
from ...drawers.linear.labels import LabelDrawer
from ....labels.placement import prepare_label_list_linear
from ....features.factory import create_feature_dict
from ....features.colors import preprocess_color_tables
from ....labels.filtering import preprocess_label_filtering
from ....configurators import FeatureDrawingConfigurator
from .length_bar import (
    RULER_LABEL_OFFSET,
    RULER_TICK_LENGTH,
    auto_linear_tick_interval,
    format_linear_tick_label,
)


class SeqRecordGroup:
    """Manages the visualization of a SeqRecord in a linear layout."""

    def __init__(
        self,
        gb_record: SeqRecord,
        canvas_config: LinearCanvasConfigurator,
        feature_config: FeatureDrawingConfigurator,
        config_dict: dict,
        precalculated_labels: Optional[list] = None,
        cfg: GbdrawConfig | None = None,
    ) -> None:
        self.gb_record = gb_record
        self.canvas_config = canvas_config
        self.length_param = self.canvas_config.length_param
        self.feature_config = feature_config
        self.config_dict = config_dict
        self.precalculated_labels = precalculated_labels
        cfg = cfg or GbdrawConfig.from_dict(config_dict)
        self._cfg = cfg

        raw_show_labels = cfg.canvas.show_labels
        self.show_labels = (raw_show_labels != "none") if isinstance(raw_show_labels, str) else bool(raw_show_labels)

        self.label_stroke_color = cfg.labels.stroke_color.label_stroke_color
        # Keep legacy behavior: linear label leader line width uses the "long" value.
        self.label_stroke_width = cfg.labels.stroke_width.long
        self.label_filtering = cfg.labels.filtering.as_dict()
        self.separate_strands = self.canvas_config.strandedness
        self.resolve_overlaps = self.canvas_config.resolve_overlaps
        self.normalize_length = cfg.canvas.linear.normalize_length
        scale_cfg = cfg.objects.scale
        self.scale_style = str(scale_cfg.style).strip().lower()
        self.scale_stroke_color = scale_cfg.stroke_color
        self.scale_label_color = scale_cfg.label_color
        self.scale_stroke_width = scale_cfg.stroke_width
        self.scale_font_size = scale_cfg.font_size.for_length_param(self.length_param)
        self.ruler_label_font_size = scale_cfg.ruler_label_font_size.for_length_param(self.length_param)
        self.scale_font_weight = scale_cfg.font_weight
        self.scale_font_family = cfg.objects.text.font_family
        self.scale_interval = scale_cfg.interval
        self.auto_scale_interval = auto_linear_tick_interval(max(1, int(self.canvas_config.longest_genome)))
        self.scale_label_context_length = max(1, int(self.canvas_config.longest_genome))
        self.axis_stroke_width = self._cfg.objects.axis.linear.stroke_width.for_length_param(self.length_param)
        self.track_layout = str(self.canvas_config.track_layout).strip().lower()
        self.ruler_on_axis = bool(cfg.canvas.linear.ruler_on_axis)
        self.axis_ruler_enabled = (
            self.ruler_on_axis
            and self.scale_style == "ruler"
            and self.track_layout in {"above", "below"}
        )
        self.record_group: Group = self.setup_record_group()

    def draw_linear_axis(self, alignment_width: float, genome_size_normalization_factor: float) -> Line:
        """
        Draws a linear axis for the genomic record visualization.

        Args:
            alignment_width (float): Width of the alignment area.
            genome_size_normalization_factor (float): Normalization factor based on the genome size.

        Returns:
            Line: An SVG line element representing the linear axis.
        """
        bar_length: float = alignment_width * genome_size_normalization_factor
        linear_axis_stroke_color: str = self._cfg.objects.axis.linear.stroke_color
        start_x: float = 0
        start_y: float = 0
        end_x: float = bar_length
        end_y = 0
        axis_path = Line(
            start=(start_x, start_y),
            end=(end_x, end_y),
            stroke=linear_axis_stroke_color,
            stroke_width=self.axis_stroke_width,
            fill="none",
        )
        return axis_path

    def _resolve_record_coordinate_span(self, record_length: int) -> tuple[int, int]:
        annotations = getattr(self.gb_record, "annotations", None) or {}
        try:
            start_coord = int(annotations.get("gbdraw_coord_base", 1))
        except (TypeError, ValueError):
            start_coord = 1
        try:
            step = int(annotations.get("gbdraw_coord_step", 1))
        except (TypeError, ValueError):
            step = 1
        if step == 0:
            step = 1
        step = 1 if step > 0 else -1
        end_coord = start_coord + (step * max(0, record_length - 1))
        return start_coord, end_coord

    def _axis_tick_coordinates(self, start_coord: int, end_coord: int, tick_interval: int) -> list[int]:
        if tick_interval <= 0:
            return []
        min_coord = min(start_coord, end_coord)
        max_coord = max(start_coord, end_coord)
        first_tick = (min_coord // tick_interval) * tick_interval
        if first_tick < min_coord:
            first_tick += tick_interval
        coordinates: list[int] = []
        tick = first_tick
        while tick < max_coord:
            if tick > min_coord:
                coordinates.append(int(tick))
            tick += tick_interval
        return sorted(coordinates, reverse=(start_coord > end_coord))

    def _draw_axis_ruler(self, group: Group, *, bar_length: float, record_length: int) -> None:
        if not self.axis_ruler_enabled or bar_length <= 0 or record_length <= 0:
            return

        start_coord, end_coord = self._resolve_record_coordinate_span(record_length)
        coord_span = abs(end_coord - start_coord)
        interval = (
            int(self.scale_interval)
            if (self.scale_interval is not None and int(self.scale_interval) > 0)
            else int(self.auto_scale_interval)
        )
        if interval <= 0:
            return

        tick_coords = self._axis_tick_coordinates(start_coord, end_coord, interval)
        if not tick_coords:
            return

        label_y = RULER_LABEL_OFFSET
        tick_start_y = 0 - (0.5 * self.axis_stroke_width)
        tick_end_y = RULER_TICK_LENGTH
        dominant_baseline = "hanging"
        if self.track_layout == "below":
            label_y = -RULER_LABEL_OFFSET
            tick_start_y = 0 + (0.5 * self.axis_stroke_width)
            tick_end_y = -RULER_TICK_LENGTH
            dominant_baseline = "text-after-edge"

        for coord in tick_coords:
            if coord_span <= 0:
                x_pos = 0.0
            else:
                x_pos = bar_length * (abs(coord - start_coord) / float(coord_span))
            if x_pos < 0.0:
                x_pos = 0.0
            elif x_pos > bar_length:
                x_pos = bar_length

            tick_line = Line(
                start=(x_pos, tick_start_y),
                end=(x_pos, tick_end_y),
                stroke=self.scale_stroke_color,
                stroke_width=self.axis_stroke_width,
            )
            group.add(tick_line)

            label_text = format_linear_tick_label(
                int(coord),
                context_length=self.scale_label_context_length,
                tick_interval=interval,
            )
            text_element = Text(
                label_text,
                insert=(x_pos, label_y),
                stroke="none",
                fill=self.scale_label_color,
                font_size=self.ruler_label_font_size,
                font_weight=self.scale_font_weight,
                font_family=self.scale_font_family,
                text_anchor="middle",
                dominant_baseline=dominant_baseline,
            )
            group.add(text_element)

    def draw_record(
        self,
        feature_dict: dict,
        record_length: int,
        cds_height: float,
        alignment_width: float,
        genome_size_normalization_factor: float,
        separate_strands: bool,
        arrow_length: float,
        group: Group,
        label_list: list,
    ) -> Group:
        """Draws the genomic features onto the provided SVG group."""
        bar_length: float = alignment_width * genome_size_normalization_factor
        # Draw the axis
        axis_path = self.draw_linear_axis(alignment_width, genome_size_normalization_factor)
        group.add(axis_path)
        self._draw_axis_ruler(group, bar_length=bar_length, record_length=record_length)

        # Process labels if enabled
        if self.show_labels:
            for label in label_list:
                if not label["is_embedded"]:
                    label_middle_y = float(label["middle_y"])
                    feature_middle_y = float(label["feature_middle_y"])
                    label_height = float(label["height_px"])
                    if label_middle_y >= feature_middle_y:
                        # Label is below the feature track, so connect to its upper edge.
                        label_edge_y = label_middle_y - (0.45 * label_height)
                    else:
                        # Label is above the feature track, so connect to its lower edge.
                        label_edge_y = label_middle_y + (0.45 * label_height)
                    line_path = Line(
                        start=(label["middle"], label["feature_middle_y"]),
                        end=(label["middle"], label_edge_y),
                        stroke=self.label_stroke_color,
                        stroke_width=self.label_stroke_width,
                    )
                    group.add(line_path)

        # Draw features
        for feature_object in feature_dict.values():
            feature_strand = feature_object.strand
            group = FeatureDrawer(self.feature_config).draw(
                feature_object=feature_object,
                group=group,
                genome_length=record_length,
                cds_height=cds_height,
                alignment_width=alignment_width,
                normalization_factor=genome_size_normalization_factor,
                feature_strand=feature_strand,
                separate_strands=separate_strands,
                arrow_length=arrow_length,
                track_layout=self.canvas_config.track_layout,
                track_axis_gap=self.canvas_config.track_axis_gap,
            )

        # Add labels
        if self.show_labels:
            for label in label_list:
                group = LabelDrawer(self.config_dict).draw(label, group)

        return group

    def setup_record_group(self) -> Group:
        """
        Sets up the SVG group for the SeqRecord visualization.

        Returns:
            Group: The SVG group prepared for SeqRecord visualization.
        """
        alignment_width: float = self.canvas_config.alignment_width
        cds_height: float = self.canvas_config.cds_height
        longest_genome: int = self.canvas_config.longest_genome
        arrow_length: float = self.canvas_config.arrow_length
        color_table: DataFrame | None = self.feature_config.color_table

        separate_strands = self.canvas_config.strandedness
        resolve_overlaps = self.canvas_config.resolve_overlaps
        label_filtering = self.label_filtering  # type: ignore
        track_id = self.gb_record.id

        record_group: Group = Group(id=track_id)
        record_group = Group(id=track_id)

        record_length: int = len(self.gb_record.seq)

        if self.normalize_length:
            genome_size_normalization_factor = 1.0
        else:
            genome_size_normalization_factor: float = record_length / longest_genome

        selected_features_set: str = self.feature_config.selected_features_set

        default_colors: DataFrame | None = self.feature_config.default_colors
        label_filtering = preprocess_label_filtering(self.label_filtering)
        color_table, default_colors = preprocess_color_tables(color_table, default_colors)
        feature_dict, _ = create_feature_dict(
            self.gb_record,
            color_table,
            selected_features_set,
            default_colors,
            separate_strands,
            resolve_overlaps,
            label_filtering,
            directional_feature_types=self.feature_config.directional_feature_types,
        )
        label_list = []
        if self.show_labels:
            label_list = prepare_label_list_linear(
                feature_dict,
                record_length,
                alignment_width,
                genome_size_normalization_factor,
                cds_height,
                separate_strands,
                self.canvas_config.track_layout,
                self.canvas_config.track_axis_gap,
                self.config_dict,
                cfg=self._cfg,
            )

        record_group: Group = self.draw_record(
            feature_dict,
            record_length,
            cds_height,
            alignment_width,
            genome_size_normalization_factor,
            separate_strands,
            arrow_length,
            record_group,
            label_list,
        )
        return record_group

    def get_group(self) -> Group:
        """
        Retrieves the SVG group containing the SeqRecord visualization.

        Returns:
            Group: The SVG group with SeqRecord visualization elements.
        """
        return self.record_group


__all__ = ["SeqRecordGroup"]


