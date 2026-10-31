#!/usr/bin/env python
# coding: utf-8

from typing import NamedTuple

from Bio.SeqRecord import SeqRecord
from gbdraw.svg.elements import Group, Line, Rect, Text

from ....auto_sizes import linear_tick_interval
from ....canvas.linear import LinearCanvasGeometry
from ....core.record_metadata import _read_coord_map, _record_source_start
from ...drawers.linear.features import FeatureDrawer
from ...drawers.label_backing import feature_label_box, label_halo, underline_casing
from ...drawers.linear.labels import LabelDrawer
from ....labels.label_ink import linear_label_ink
from ....labels.label_style import label_leader_width_px
from ....labels.linear import linear_leader_polyline, linear_record_labels, place_linear_label_row
from ....labels.linear_types import LinearLabel
from ....config.models.render_profiles import LinearRenderProfile
from ....layout.linear import LinearFeatureLaneGeometry, LinearFeatureTrackLayout
from ....features.factory import FeatureBuildResult
from ....configurators import FeatureDrawingConfigurator
from ....svg.ids import record_group_svg_id, stable_svg_id
from ....layout.record_labels import source_ruler_ticks
from ...label_binding import bind_label_part
from .feature_identity import LinearFeatureDomIndex
from .length_bar import (
    RULER_LABEL_OFFSET,
    RULER_TICK_LENGTH,
    format_linear_tick_label,
)


from gbdraw.layout.record_coordinates import RecordDisplayTransform

class _LabelParts(NamedTuple):
    """A feature label's drawn style parts: whether it has a leader, the leader line the leader layer draws
    (two points, or none when the whole leader lies within its casing), the leader's tail in front of its
    underline casing, the casing, the underline and the box (None when not drawn)."""

    has_leader: bool
    leader: tuple[tuple[float, float], ...]
    tail: Line | None = None
    casing: Line | None = None
    underline: Line | None = None
    box: Rect | None = None


class SeqRecordGroup:
    """Manages the visualization of a SeqRecord in a linear layout."""

    def __init__(
        self,
        gb_record: SeqRecord,
        geometry: LinearCanvasGeometry,
        feature_config: FeatureDrawingConfigurator,
        feature_layers: FeatureBuildResult,
        profile: LinearRenderProfile,
        feature_track_layout: LinearFeatureTrackLayout,
        axis_ruler_enabled: bool,
        feature_dom_index: LinearFeatureDomIndex,
        precalculated_labels: list[LinearLabel] | None = None,
        draw_features: bool = True,
        label_font_size: float | None = None,
        orthogroup_label_member_ids: set[str | int] | None = None,
        orthogroup_label_top_member_ids: set[str | int] | None = None,
        record_index: int = 0,
        record_count: int = 1,
        group_id: str | None = None,
        record_local_ruler: bool = False,
        feature_offset_y: float = 0.0,
        feature_lane_geometry: LinearFeatureLaneGeometry | None = None,
        record_transform: RecordDisplayTransform | None = None,
        *,
        sequence_width: float,
    ) -> None:
        self.record_transform = record_transform
        self.gb_record = gb_record
        self.geometry = geometry
        self.length_param = geometry.length_param
        self.feature_config = feature_config
        self.feature_layers = feature_layers
        self.feature_dom_index = feature_dom_index
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
        # The record's drawn axis width, from its placement (``solve_linear_horizontal_layout``).
        self.sequence_width = float(sequence_width)
        self.record_local_ruler = bool(record_local_ruler)
        self.feature_offset_y = float(feature_offset_y)
        self.feature_lane_geometry = feature_lane_geometry
        self.profile = profile
        self._feature_track_layout = feature_track_layout
        self._axis_ruler_enabled = bool(axis_ruler_enabled)
        cfg = profile.config
        self._cfg = cfg

        self.show_labels = profile.record_has_labels(self.record_index)

        self.label_stroke_color = cfg.labels.stroke_color.label_stroke_color
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
        self.auto_scale_interval = linear_tick_interval(max(1, int(self.geometry.longest_genome)))
        self.scale_label_context_length = max(1, int(self.geometry.longest_genome))
        self.axis_stroke_width = self._cfg.objects.axis.linear.stroke_width.for_length_param(self.length_param)
        self.record_group: Group = self.setup_record_group()

    @property
    def separate_strands(self) -> bool:
        return self.profile.strandedness

    @property
    def track_layout(self) -> str:
        return self._feature_track_layout

    @property
    def axis_ruler_enabled(self) -> bool:
        return self._axis_ruler_enabled

    def draw_linear_axis(self, alignment_width: float) -> Line:
        """
        Draws a linear axis for the genomic record visualization.

        Args:
            alignment_width (float): Width of the alignment area.

        Returns:
            Line: An SVG line element representing the linear axis.
        """
        bar_length: float = alignment_width
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

    def _leader_part(self, start: tuple[float, float], end: tuple[float, float], feature_id: str, width: float) -> Line:
        """One straight part of a feature label's leader, or its underline: one style for both (LP-D95)."""

        part = Line(start=start, end=end, stroke=self.label_stroke_color, stroke_width=width)
        bind_label_part(part, feature_id)
        return part

    def _label_style_parts(self, label: LinearLabel, leader: tuple[tuple[float, float], ...]) -> _LabelParts:
        """A label's drawn parts in the style of its mode (``FeatureLabelStyle``)."""

        style = self.profile.label_style
        if not style.draws_from_ink:
            return _LabelParts(bool(leader), leader)
        ink = linear_label_ink(label, leader, style)
        feature_id = label["feature_id"]
        box = feature_label_box(ink, feature_id) if ink.draws_box(style) else None
        if not leader:
            # A turned label on its feature has its box too (OV-1361).
            return _LabelParts(False, leader, box=box)
        width = label_leader_width_px(label)
        underline = casing = tail = None
        line = ink.underline
        if line is not None:
            # The leader ends on its underline, which continues it (LP-D88), or stops above the letters.
            leader = (leader[0], line[0])
            if line[1] != line[2]:
                underline = self._leader_part(line[1], line[2], feature_id, width)
                # The leader layer draws the leader up to its casing; its tail is drawn in front of the casing.
                cased = underline_casing(line[1], line[2], float(label["font_size"]))
                entry = cased.entry(leader[0], line[0])
                casing = cased.element(feature_id)
                tail = self._leader_part(entry, line[0], feature_id, width)
                leader_layer = () if entry == leader[0] else (leader[0], entry)
                return _LabelParts(True, leader_layer, tail, casing, underline, box)
        return _LabelParts(True, leader, None, None, None, box)

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
        projected_ticks = None
        if self.record_transform is not None and self.record_transform.projects:
            projected_ticks = {tick.base: tick.offset for tick in source_ruler_ticks(self.record_transform, interval)}
            tick_coords = list(projected_ticks)
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
            if projected_ticks is not None:
                x_pos = bar_length * projected_ticks[coord] / record_length
            elif coord_span <= 0:
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
        separate_strands: bool,
        arrow_length: float,
        group: Group,
        label_list: list[LinearLabel],
    ) -> Group:
        """Draws the genomic features onto the provided SVG group."""
        bar_length: float = alignment_width
        # Draw the axis
        axis_path = self.draw_linear_axis(alignment_width)
        group.add(axis_path)
        self._draw_axis_ruler(group, bar_length=bar_length, record_length=record_length)

        feature_group = group
        has_feature_shift = abs(self.feature_offset_y) > 1e-9
        if has_feature_shift:
            feature_group = Group(id=f"{self.record_group_id}_features")

        # Process labels if enabled. Each label gets its feature identity on a
        # copy: the prepared labels are this content's arguments.
        styled: list[_LabelParts] = []
        if self.show_labels and self.draw_features_enabled:
            bound_labels: list[LinearLabel] = []
            for label in label_list:
                source_feature_index = label["source_feature_index"]
                feature_id = (
                    None
                    if source_feature_index is None
                    else self.feature_dom_index.by_source_index.get((self.record_index, source_feature_index))
                )
                if not feature_id:
                    raise ValueError(
                        "Prepared linear label has no rendered feature identity."
                    )
                label = {**label, "feature_id": feature_id}
                bound_labels.append(label)
                parts = self._label_style_parts(label, linear_leader_polyline(label))
                styled.append(parts)
                if parts.leader:
                    feature_group.add(
                        self._leader_part(
                            parts.leader[0], parts.leader[1], feature_id, label_leader_width_px(label)
                        )
                    )
            label_list = bound_labels

        # Draw features
        if self.draw_features_enabled:
            feature_drawer = FeatureDrawer(self.feature_config)
            for feature_object in feature_dict.values():
                feature_strand = feature_object.strand
                feature_group = feature_drawer.draw(
                    feature_object=feature_object,
                    group=feature_group,
                    genome_length=record_length,
                    cds_height=cds_height,
                    alignment_width=alignment_width,
                    feature_strand=feature_strand,
                    separate_strands=separate_strands,
                    arrow_length=arrow_length,
                    track_layout=self.track_layout,
                    track_axis_gap=self.geometry.track_axis_gap,
                    feature_lane_geometry=self.feature_lane_geometry,
                    record_index=self.record_index,
                    record_count=self.record_count,
                    feature_instance_id=feature_object.rendered_instance_index,
                )

        # Add labels: in front of every leader and glyph, the boxes, the halos, the underline casings, the
        # leader tails, the underlines, then the letters (a halo or a casing must not erode an underline).
        if self.show_labels:
            letters_at = len(feature_group.elements)
            for label in label_list:
                feature_group = LabelDrawer().draw(label, feature_group)
            texts = feature_group.elements[letters_at:]
            halos = [
                label_halo(text, stable_svg_id("label_text", self.record_group_id, label["feature_id"], index))
                for index, (label, text, parts) in enumerate(zip(label_list, texts, styled))
                if parts.has_leader and self.profile.label_style.protection == "halo"
            ]
            feature_group.elements[letters_at:letters_at] = [
                *(parts.box for parts in styled if parts.box is not None),
                *halos,
                *(parts.casing for parts in styled if parts.casing is not None),
                *(parts.tail for parts in styled if parts.tail is not None),
                *(parts.underline for parts in styled if parts.underline is not None),
            ]

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
        alignment_width: float = self.sequence_width
        cds_height: float = self.geometry.cds_height
        arrow_length: float = self.geometry.arrow_length

        separate_strands = self.separate_strands
        record_group: Group = Group(id=self.record_group_id)
        record_group.attribs["data-gbdraw-record-id"] = str(self.gb_record.id)
        record_group.attribs["data-gbdraw-record-index"] = str(self.record_index)

        record_length: int = len(self.gb_record.seq)
        # Input-file span of a cropped or reverse-complemented record (PD-OI-076):
        # match popups and FASTA headers report source coordinates through it.
        source_base, source_step = _read_coord_map(self.gb_record)
        if record_length and (source_base, source_step) != (1, 1):
            source_low = _record_source_start(self.gb_record)
            record_group.attribs["data-gbdraw-record-source-start"] = str(source_low)
            record_group.attribs["data-gbdraw-record-source-end"] = str(source_low + record_length - 1)
            record_group.attribs["data-gbdraw-record-source-step"] = str(source_step)

        feature_dict = self.feature_layers.foreground_features
        label_list: list[LinearLabel] = []
        if self.show_labels and self.draw_features_enabled:
            if self.precalculated_labels is not None:
                label_list = self.precalculated_labels
            else:
                record_labels = linear_record_labels(
                    feature_dict,
                    record_length,
                    alignment_width,
                    cds_height,
                    separate_strands,
                    self.track_layout,
                    self.geometry.track_axis_gap,
                    profile=self.profile,
                    label_font_size=self.label_font_size,
                    orthogroup_label_member_ids=self.orthogroup_label_member_ids,
                    orthogroup_label_top_member_ids=self.orthogroup_label_top_member_ids,
                    labels_overrides_only=not self.profile.labels_in_scope(
                        self.record_index
                    ),
                    feature_lane_geometry=self.feature_lane_geometry,
                    record_transform=self.record_transform,
                )
                label_list = place_linear_label_row(
                    [(record_labels, 0.0)], profile=self.profile, track_layout=self.track_layout
                )[0]

        record_group = self.draw_record(
            feature_dict,
            record_length,
            cds_height,
            alignment_width,
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
