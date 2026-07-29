#!/usr/bin/env python
# coding: utf-8

from collections import Counter
from typing import Optional

from Bio.SeqRecord import SeqRecord
from svgwrite.container import Group
from svgwrite.shapes import Line
from svgwrite.text import Text

from ....canvas import LinearCanvasConfigurator
from ...drawers.linear.features import FeatureDrawer
from ...drawers.linear.labels import LabelDrawer
from ....labels.linear import prepare_label_list_linear
from ....layout.linear import LinearFeatureLaneGeometry, LinearRecordRenderContext
from ....features.factory import FeatureBuildResult
from ....configurators import FeatureDrawingConfigurator
from ....svg.ids import record_group_svg_id
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
        feature_layers: FeatureBuildResult,
        render_context: LinearRecordRenderContext,
        precalculated_labels: Optional[list] = None,
        draw_features: bool = True,
        label_font_size: float | None = None,
        orthogroup_label_member_ids: set[str] | None = None,
        orthogroup_label_top_member_ids: set[str] | None = None,
        record_index: int = 0,
        record_count: int = 1,
        group_id: str | None = None,
        sequence_width: float | None = None,
        record_local_ruler: bool = False,
        feature_offset_y: float = 0.0,
        feature_lane_geometry: LinearFeatureLaneGeometry | None = None,
    ) -> None:
        self.gb_record = gb_record
        self.canvas_config = canvas_config
        self.length_param = self.canvas_config.length_param
        self.feature_config = feature_config
        self.feature_layers = feature_layers
        self.precalculated_labels = precalculated_labels
        self.draw_features_enabled = bool(draw_features)
        self.orthogroup_label_member_ids = orthogroup_label_member_ids
        self.orthogroup_label_top_member_ids = orthogroup_label_top_member_ids
        self.record_index = int(record_index)
        self.record_count = int(record_count)
        self.record_group_id = (
            str(group_id)
            if group_id
            else record_group_svg_id(
                gb_record.id,
                mode="linear",
                record_index=self.record_index,
                record_count=self.record_count,
            )
        )
        self.sequence_width = float(sequence_width) if sequence_width is not None else None
        self.record_local_ruler = bool(record_local_ruler)
        self.feature_offset_y = float(feature_offset_y)
        self.feature_lane_geometry = feature_lane_geometry
        self.render_context = render_context
        cfg = render_context.profile.config
        self._cfg = cfg

        label_scope = render_context.profile.label_scope
        self.show_labels = (
            label_scope == "all"
            or label_scope == "orthogroup_top"
            or (label_scope == "first" and self.record_index == 0)
        )

        self.label_stroke_color = cfg.labels.stroke_color.label_stroke_color
        # Keep legacy behavior: linear label leader line width uses the "long" value.
        self.label_stroke_width = cfg.labels.stroke_width.long
        self.normalize_length = cfg.canvas.linear.normalize_length
        scale_cfg = cfg.objects.scale
        self.scale_style = str(scale_cfg.style).strip().lower()
        self.scale_stroke_color = scale_cfg.stroke_color
        self.scale_label_color = scale_cfg.label_color
        self.scale_stroke_width = scale_cfg.stroke_width
        self.scale_font_size = scale_cfg.font_size.for_length_param(self.length_param)
        self.ruler_label_font_size = scale_cfg.ruler_label_font_size.for_length_param(self.length_param)
        self.label_font_size = (
            float(label_font_size)
            if label_font_size is not None
            else cfg.labels.font_size.linear.for_length_param(self.length_param)
        )
        self.scale_font_weight = scale_cfg.font_weight
        self.scale_font_family = cfg.objects.text.font_family
        self.scale_interval = scale_cfg.interval
        self.auto_scale_interval = auto_linear_tick_interval(max(1, int(self.canvas_config.longest_genome)))
        self.scale_label_context_length = max(1, int(self.canvas_config.longest_genome))
        self.axis_stroke_width = self._cfg.objects.axis.linear.stroke_width.for_length_param(self.length_param)
        self.record_group: Group = self.setup_record_group()

    @property
    def separate_strands(self) -> bool:
        return self.render_context.profile.strandedness

    @property
    def track_layout(self) -> str:
        return self.render_context.feature_track_layout

    @property
    def axis_ruler_enabled(self) -> bool:
        return self.render_context.axis_ruler_enabled

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
        if self.record_local_ruler:
            return 0, record_length
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

        feature_group = group
        has_feature_shift = abs(self.feature_offset_y) > 1e-9
        if has_feature_shift:
            feature_group = Group(id=f"{self.record_group_id}_features", debug=False)

        # Process labels if enabled
        if self.show_labels and self.draw_features_enabled:
            for label in label_list:
                if label.get("leader_line"):
                    line_path = Line(
                        start=(label["leader_start_x"], label["leader_start_y"]),
                        end=(label["leader_end_x"], label["leader_end_y"]),
                        stroke=self.label_stroke_color,
                        stroke_width=self.label_stroke_width,
                    )
                    feature_group.add(line_path)
                elif not label["is_embedded"]:
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
                    feature_group.add(line_path)

        # Draw features
        if self.draw_features_enabled:
            feature_drawer = FeatureDrawer(self.feature_config)
            feature_id_counts = Counter(
                feature_drawer.get_feature_data_id(feature_object)
                for feature_object in feature_dict.values()
            )
            for feature_object in feature_dict.values():
                feature_strand = feature_object.strand
                stable_feature_id = feature_drawer.get_feature_data_id(feature_object)
                feature_group = feature_drawer.draw(
                    feature_object=feature_object,
                    group=feature_group,
                    genome_length=record_length,
                    cds_height=cds_height,
                    alignment_width=alignment_width,
                    normalization_factor=genome_size_normalization_factor,
                    feature_strand=feature_strand,
                    separate_strands=separate_strands,
                    arrow_length=arrow_length,
                    track_layout=self.track_layout,
                    track_axis_gap=self.canvas_config.track_axis_gap,
                    feature_lane_geometry=self.feature_lane_geometry,
                    record_index=self.record_index,
                    record_count=self.record_count,
                    feature_instance_id=(
                        feature_object.feature_id
                        if stable_feature_id
                        and feature_id_counts[stable_feature_id] > 1
                        else None
                    ),
                )

        # Add labels
        if self.show_labels:
            for label in label_list:
                feature_group = LabelDrawer().draw(label, feature_group)

        if has_feature_shift:
            feature_group.translate(0.0, self.feature_offset_y)
            feature_group.attribs["data-gbdraw-feature-offset-y"] = f"{self.feature_offset_y:g}"
            group.add(feature_group)

        return group

    def setup_record_group(self) -> Group:
        """
        Sets up the SVG group for the SeqRecord visualization.

        Returns:
            Group: The SVG group prepared for SeqRecord visualization.
        """
        alignment_width: float = (
            self.sequence_width
            if self.sequence_width is not None
            else self.canvas_config.alignment_width
        )
        cds_height: float = self.canvas_config.cds_height
        longest_genome: int = self.canvas_config.longest_genome
        arrow_length: float = self.canvas_config.arrow_length

        separate_strands = self.separate_strands
        record_group: Group = Group(id=self.record_group_id, debug=False)
        record_group.attribs["data-gbdraw-record-id"] = str(self.gb_record.id)
        record_group.attribs["data-gbdraw-record-index"] = str(self.record_index)

        record_length: int = len(self.gb_record.seq)

        if self.sequence_width is not None:
            genome_size_normalization_factor = 1.0
        elif self.normalize_length:
            genome_size_normalization_factor = 1.0
        else:
            genome_size_normalization_factor: float = record_length / longest_genome

        feature_dict = self.feature_layers.foreground_features
        label_list = []
        if self.show_labels and self.draw_features_enabled:
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
                    self.track_layout,
                    self.canvas_config.track_axis_gap,
                    cfg=self._cfg,
                    label_font_size=self.label_font_size,
                    orthogroup_label_member_ids=self.orthogroup_label_member_ids,
                    orthogroup_label_top_member_ids=self.orthogroup_label_top_member_ids,
                    feature_lane_geometry=self.feature_lane_geometry,
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
