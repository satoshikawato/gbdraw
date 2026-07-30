#!/usr/bin/env python
# coding: utf-8

from collections import Counter
from typing import Optional, Dict

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]
from svgwrite.container import Group  # type: ignore[reportMissingImports]

from ....canvas import CircularCanvasConfigurator  # type: ignore[reportMissingImports]
from ....core.sequence import determine_length_parameter  # type: ignore[reportMissingImports]
from ....features.factory import FeatureBuildResult  # type: ignore[reportMissingImports]
from ....features.objects import FeatureObject  # type: ignore[reportMissingImports]
from ....labels.circular import prepare_label_list  # type: ignore[reportMissingImports]
from ....layout.circular import CircularRecordRenderContext  # type: ignore[reportMissingImports]
from ....svg.ids import record_group_svg_id
from ...drawers.circular.labels import LabelDrawer  # type: ignore[reportMissingImports]
from ...drawers.circular.features import FeatureDrawer  # type: ignore[reportMissingImports]
from ....configurators import FeatureDrawingConfigurator  # type: ignore[reportMissingImports]


class SeqRecordGroup:
    """
    Manages the creation and visualization of a SeqRecord (genomic data) on a circular canvas.
    """

    def __init__(
        self,
        gb_record: SeqRecord,
        canvas_config: CircularCanvasConfigurator,
        feature_config: FeatureDrawingConfigurator,
        feature_layers: FeatureBuildResult,
        render_context: CircularRecordRenderContext,
        precalculated_labels: Optional[list[dict]] = None,
        feature_track_ratio_factor_override: float | None = None,
        feature_anchor_radius_px: float | None = None,
        record_index: int = 0,
        group_id: str | None = None,
        slot_id: str | None = None,
        feature_dom_namespace: str | None = None,
    ) -> None:
        self.gb_record: SeqRecord = gb_record
        self.canvas_config: CircularCanvasConfigurator = canvas_config
        self.feature_config: FeatureDrawingConfigurator = feature_config
        self.feature_layers = feature_layers
        self.render_context = render_context
        cfg = render_context.profile.config

        self.length_threshold = cfg.labels.length_threshold.circular
        self.length_param = determine_length_parameter(len(gb_record.seq), self.length_threshold)
        self.font_family: str = cfg.objects.text.font_family

        self.show_labels = render_context.profile.labels_enabled

        self.label_stroke_width = cfg.labels.stroke_width.for_length_param(self.length_param)
        self.label_stroke_color = cfg.labels.stroke_color.label_stroke_color
        self.font_size = cfg.labels.font_size.for_length_param(self.length_param)
        self.dpi = cfg.canvas.dpi
        self.strandedness = render_context.profile.strandedness
        self.track_ratio_factors = cfg.canvas.circular.track_ratio_factors[self.length_param]
        self.track_ratio = self.canvas_config.track_ratio
        self.precalculated_labels: Optional[list[dict]] = precalculated_labels
        self.feature_track_ratio_factor_override = feature_track_ratio_factor_override
        self.feature_anchor_radius_px = feature_anchor_radius_px
        self.record_index = int(record_index)
        self.record_group_id = (
            str(group_id)
            if group_id
            else record_group_svg_id(
                gb_record.id,
                mode="circular",
                record_index=self.record_index,
            )
        )
        self.slot_id = str(slot_id) if slot_id else None
        self.feature_dom_namespace = (
            str(feature_dom_namespace) if feature_dom_namespace else None
        )
        self.record_group: Group = self.setup_record_group()

    def draw_record(self, feature_dict: Dict[str, FeatureObject], record_length: int, group: Group) -> Group:
        label_list = []
        feature_anchor_radius = (
            float(self.feature_anchor_radius_px)
            if self.feature_anchor_radius_px is not None
            else float(self.canvas_config.radius)
        )
        if self.show_labels is True:
            # Reuse pre-calculated labels when available to avoid repeating heavy placement work.
            if self.precalculated_labels is not None:
                label_list = self.precalculated_labels
            else:
                label_list = prepare_label_list(
                    feature_dict,
                    record_length,
                    feature_anchor_radius,
                    self.track_ratio,
                    self.render_context.profile,
                    feature_track_ratio_factor_override=self.feature_track_ratio_factor_override,
                    feature_layout=self.render_context.feature_layout,
                    track_preset=self.render_context.track_preset,
                    feature_lane_direction=str(
                        self.render_context.feature_lane_direction
                    ),
                )
        feature_track_ratio_factor = (
            float(self.feature_track_ratio_factor_override)
            if self.feature_track_ratio_factor_override is not None
            else float(self.track_ratio_factors[0])
        )
        feature_drawer = FeatureDrawer(
            self.feature_config,
            self.render_context.feature_layout,
        )
        feature_id_counts = Counter(
            feature_drawer.get_feature_data_id(feature_object)
            for feature_object in feature_dict.values()
        )
        for feature_object in feature_dict.values():
            stable_feature_id = feature_drawer.get_feature_data_id(feature_object)
            group = feature_drawer.draw(
                feature_object,
                group,
                record_length,
                feature_anchor_radius,
                self.canvas_config.track_ratio,
                feature_track_ratio_factor,
                self.render_context.track_preset,
                self.strandedness,
                self.length_param,
                feature_instance_id=(
                    feature_object.source_feature_index
                    if stable_feature_id
                    and feature_id_counts[stable_feature_id] > 1
                    else None
                ),
                feature_dom_namespace=self.feature_dom_namespace,
            )
        if self.show_labels:
            label_drawer = LabelDrawer(profile=self.render_context.profile)
            for label in label_list:
                if label["is_embedded"]:
                    group = label_drawer.draw(
                        label,
                        group,
                        record_length,
                        feature_anchor_radius,
                        self.canvas_config.track_ratio,
                        feature_track_ratio_factor_override=self.feature_track_ratio_factor_override,
                        track_preset=self.render_context.track_preset,
                    )
        return group

    def setup_record_group(self) -> Group:
        feature_dict = self.feature_layers.foreground_features
        record_group = Group(id=self.record_group_id, debug=False)
        record_group.attribs["data-gbdraw-record-id"] = str(self.gb_record.id)
        record_group.attribs["data-gbdraw-record-index"] = str(self.record_index)
        if self.slot_id:
            record_group.attribs["data-gbdraw-slot-id"] = self.slot_id
            record_group.attribs["data-gbdraw-slot-renderer"] = "features"
        record_length: int = len(self.gb_record.seq)

        record_group: Group = self.draw_record(feature_dict, record_length, record_group)
        return record_group

    def get_group(self) -> Group:
        return self.record_group


__all__ = ["SeqRecordGroup"]
