#!/usr/bin/env python
# coding: utf-8

from typing import Optional

from pandas import DataFrame
from Bio.SeqRecord import SeqRecord
from svgwrite.container import Group
from svgwrite.shapes import Line

from ....canvas import LinearCanvasConfigurator
from ....config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from ...drawers.linear.features import FeatureDrawer
from ...drawers.linear.labels import LabelDrawer
from ....labels.placement import prepare_label_list_linear
from ....features.factory import create_feature_dict
from ....features.colors import preprocess_color_tables
from ....labels.filtering import preprocess_label_filtering
from ....configurators import FeatureDrawingConfigurator


class SeqRecordGroup:
    """Manages the visualization of a SeqRecord in a linear layout."""

    def __init__(
        self,
        gb_record: SeqRecord,
        canvas_config: LinearCanvasConfigurator,
        feature_config: FeatureDrawingConfigurator,
        config_dict: dict,
        precalculated_labels: Optional[list] = None,
        precalculated_feature_dict: Optional[dict] = None,
        cfg: GbdrawConfig | None = None,
    ) -> None:
        self.gb_record = gb_record
        self.canvas_config = canvas_config
        self.length_param = self.canvas_config.length_param
        self.feature_config = feature_config
        self.config_dict = config_dict
        self.precalculated_labels = precalculated_labels
        self.precalculated_feature_dict = precalculated_feature_dict
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
        linear_axis_stroke_width: float = self._cfg.objects.axis.linear.stroke_width.for_length_param(self.length_param)
        start_x: float = 0
        start_y: float = 0
        end_x: float = bar_length
        end_y = 0
        axis_path = Line(
            start=(start_x, start_y),
            end=(end_x, end_y),
            stroke=linear_axis_stroke_color,
            stroke_width=linear_axis_stroke_width,
            fill="none",
        )
        return axis_path

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
        # Draw the axis
        axis_path = self.draw_linear_axis(alignment_width, genome_size_normalization_factor)
        group.add(axis_path)

        # Process labels if enabled
        if self.show_labels:
            for label in label_list:
                if not label["is_embedded"]:
                    line_path = Line(
                        start=(label["middle"], label["feature_middle_y"]),
                        end=(label["middle"], label["middle_y"] + 0.45 * label["height_px"]),
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
        separate_strands = self.canvas_config.strandedness
        resolve_overlaps = self.canvas_config.resolve_overlaps
        track_id = self.gb_record.id

        record_group: Group = Group(id=track_id)
        record_group = Group(id=track_id)

        record_length: int = len(self.gb_record.seq)

        if self.normalize_length:
            genome_size_normalization_factor = 1.0
        else:
            genome_size_normalization_factor: float = record_length / longest_genome

        selected_features_set: str = self.feature_config.selected_features_set
        feature_dict = self.precalculated_feature_dict
        if feature_dict is None:
            color_table: DataFrame | None = self.feature_config.color_table
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
            )
        label_list = []
        if self.show_labels:
            if self.precalculated_labels is not None:
                label_list = self.precalculated_labels
            else:
                label_list = prepare_label_list_linear(
                    feature_dict,
                    record_length,
                    alignment_width,
                    genome_size_normalization_factor,
                    cds_height,
                    separate_strands,
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


