#!/usr/bin/env python
# coding: utf-8

from typing import Dict, Sequence

from Bio.SeqRecord import SeqRecord
from gbdraw.svg.elements import Group

from ....canvas import CircularCanvasGeometry
from ....features.factory import FeatureBuildResult
from ....features.objects import FeatureObject
from ....labels.circular_types import CircularLabel
from ....layout.circular import CircularRecordRenderContext
from ....svg.ids import record_group_svg_id
from ...drawers.circular.labels import LabelDrawer
from ...drawers.circular.features import FeatureDrawer
from ....configurators import FeatureDrawingConfigurator


class SeqRecordGroup:
    """
    Manages the creation and visualization of a SeqRecord (genomic data) on a circular canvas.
    """

    def __init__(
        self,
        gb_record: SeqRecord,
        geometry: CircularCanvasGeometry,
        feature_config: FeatureDrawingConfigurator,
        feature_layers: FeatureBuildResult,
        render_context: CircularRecordRenderContext,
        precalculated_labels: Sequence[CircularLabel] | None = None,
        feature_track_ratio_factor_override: float | None = None,
        feature_anchor_radius_px: float | None = None,
        record_index: int = 0,
        group_id: str | None = None,
        slot_id: str | None = None,
        feature_dom_namespace: str | None = None,
    ) -> None:
        self.gb_record: SeqRecord = gb_record
        self.geometry = geometry
        self.feature_config: FeatureDrawingConfigurator = feature_config
        self.feature_layers = feature_layers
        self.render_context = render_context
        cfg = render_context.profile.config

        self.length_param = geometry.length_param
        self.font_family: str = cfg.objects.text.font_family

        self.show_labels = render_context.profile.labels_enabled

        self.font_size = cfg.labels.font_size.for_length_param(self.length_param)
        self.strandedness = render_context.profile.strandedness
        self.track_ratio_factors = geometry.track_ratio_factors
        self.track_ratio = geometry.track_ratio
        self.precalculated_labels: Sequence[CircularLabel] | None = precalculated_labels
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
        feature_anchor_radius = (
            float(self.feature_anchor_radius_px)
            if self.feature_anchor_radius_px is not None
            else float(self.geometry.radius)
        )
        # The record assembly places every drawn label (and checks their fit) once.
        label_list = (self.precalculated_labels or []) if self.show_labels is True else []
        feature_track_ratio_factor = (
            float(self.feature_track_ratio_factor_override)
            if self.feature_track_ratio_factor_override is not None
            else float(self.track_ratio_factors[0])
        )
        feature_drawer = FeatureDrawer(
            self.feature_config,
            self.render_context.feature_layout,
        )
        for feature_object in feature_dict.values():
            group = feature_drawer.draw(
                feature_object,
                group,
                record_length,
                feature_anchor_radius,
                self.geometry.track_ratio,
                feature_track_ratio_factor,
                self.render_context.track_preset,
                self.strandedness,
                self.length_param,
                feature_instance_id=feature_object.rendered_instance_index,
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
                        self.geometry.track_ratio,
                        feature_track_ratio_factor_override=self.feature_track_ratio_factor_override,
                        track_preset=self.render_context.track_preset,
                    )
        return group

    def setup_record_group(self) -> Group:
        feature_dict = self.feature_layers.foreground_features
        record_group = Group(id=self.record_group_id)
        record_group.attribs["data-gbdraw-record-id"] = str(self.gb_record.id)
        record_group.attribs["data-gbdraw-record-index"] = str(self.record_index)
        if self.slot_id:
            record_group.attribs["data-gbdraw-slot-id"] = self.slot_id
            record_group.attribs["data-gbdraw-slot-renderer"] = "features"
        record_length: int = len(self.gb_record.seq)

        record_group = self.draw_record(feature_dict, record_length, record_group)
        return record_group

    def get_group(self) -> Group:
        return self.record_group


__all__ = ["SeqRecordGroup"]
