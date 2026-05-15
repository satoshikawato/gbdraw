#!/usr/bin/env python
# coding: utf-8

from typing import Optional, Dict

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]
from pandas import DataFrame  # type: ignore[reportMissingImports]
from svgwrite.container import Group  # type: ignore[reportMissingImports]

from ....canvas import CircularCanvasConfigurator  # type: ignore[reportMissingImports]
from ....config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from ....core.sequence import determine_length_parameter  # type: ignore[reportMissingImports]
from ....features.factory import create_feature_dict  # type: ignore[reportMissingImports]
from ....features.colors import preprocess_color_tables  # type: ignore[reportMissingImports]
from ....features.objects import FeatureObject  # type: ignore[reportMissingImports]
from ....labels.filtering import preprocess_label_filtering  # type: ignore[reportMissingImports]
from ....labels.placement import prepare_label_list  # type: ignore[reportMissingImports]
from ...drawers.circular.labels import LabelDrawer  # type: ignore[reportMissingImports]
from ...drawers.circular.features import FeatureDrawer  # type: ignore[reportMissingImports]
from ....configurators import FeatureDrawingConfigurator  # type: ignore[reportMissingImports]
from ....diagrams.circular.radial_layout import CircularFeatureLayout  # type: ignore[reportMissingImports]


class SeqRecordGroup:
    """
    Manages the creation and visualization of a SeqRecord (genomic data) on a circular canvas.
    """

    def __init__(
        self,
        gb_record: SeqRecord,
        canvas_config: CircularCanvasConfigurator,
        feature_config: FeatureDrawingConfigurator,
        config_dict: dict,
        cfg: GbdrawConfig | None = None,
        precomputed_feature_dict: Optional[Dict[str, FeatureObject]] = None,
        precalculated_labels: Optional[list[dict]] = None,
        feature_track_ratio_factor_override: float | None = None,
        feature_anchor_radius_px: float | None = None,
    ) -> None:
        self.gb_record: SeqRecord = gb_record
        self.canvas_config: CircularCanvasConfigurator = canvas_config
        self.feature_config: FeatureDrawingConfigurator = feature_config
        self.config_dict: dict = config_dict
        cfg = cfg or GbdrawConfig.from_dict(config_dict)
        self._cfg = cfg

        self.length_threshold = cfg.labels.length_threshold.circular
        self.length_param = determine_length_parameter(len(gb_record.seq), self.length_threshold)
        self.font_family: str = cfg.objects.text.font_family

        raw_show_labels = cfg.canvas.show_labels
        self.show_labels = (raw_show_labels != "none") if isinstance(raw_show_labels, str) else bool(raw_show_labels)

        self.label_stroke_width = cfg.labels.stroke_width.for_length_param(self.length_param)
        self.label_stroke_color = cfg.labels.stroke_color.label_stroke_color
        self.label_filtering = cfg.labels.filtering.as_dict()
        self.font_size = cfg.labels.font_size.for_length_param(self.length_param)
        self.dpi = cfg.canvas.dpi
        self.track_type = getattr(self.canvas_config, "circular_track_preset", cfg.canvas.circular.track_type)
        self.feature_lane_direction = getattr(self.canvas_config, "circular_feature_lane_direction", None)
        if self.feature_lane_direction is None:
            preset = str(self.track_type).strip().lower()
            if preset == "middle":
                self.feature_lane_direction = "split"
            elif preset == "spreadout":
                self.feature_lane_direction = "outside"
            else:
                self.feature_lane_direction = "inside"
        self.strandedness = cfg.canvas.strandedness
        self.resolve_overlaps = cfg.canvas.resolve_overlaps
        self.split_overlaps_by_strand = (
            bool(self.resolve_overlaps)
            and (not bool(self.strandedness))
            and str(self.feature_lane_direction).strip().lower() == "split"
        )
        self.track_ratio_factors = cfg.canvas.circular.track_ratio_factors[self.length_param]
        self.track_ratio = self.canvas_config.track_ratio
        self.precomputed_feature_dict: Optional[Dict[str, FeatureObject]] = precomputed_feature_dict
        self.precalculated_labels: Optional[list[dict]] = precalculated_labels
        self.feature_track_ratio_factor_override = feature_track_ratio_factor_override
        self.feature_anchor_radius_px = feature_anchor_radius_px
        self.feature_layout: CircularFeatureLayout | None = getattr(
            self.canvas_config,
            "circular_feature_layout",
            None,
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
                    self.config_dict,
                    cfg=self._cfg,
                    feature_track_ratio_factor_override=self.feature_track_ratio_factor_override,
                    feature_layout=self.feature_layout,
                    track_preset=self.track_type,
                    feature_lane_direction=str(self.feature_lane_direction),
                )
        feature_track_ratio_factor = (
            float(self.feature_track_ratio_factor_override)
            if self.feature_track_ratio_factor_override is not None
            else float(self.track_ratio_factors[0])
        )
        for feature_object in feature_dict.values():
            group = FeatureDrawer(self.feature_config, self.feature_layout).draw(
                feature_object,
                group,
                record_length,
                feature_anchor_radius,
                self.canvas_config.track_ratio,
                feature_track_ratio_factor,
                self.track_type,
                self.strandedness,
                self.length_param,
            )
        if self.show_labels:
            label_drawer = LabelDrawer(self.config_dict, cfg=self._cfg)
            for label in label_list:
                if label["is_embedded"]:
                    group = label_drawer.draw(
                        label,
                        group,
                        record_length,
                        feature_anchor_radius,
                        self.canvas_config.track_ratio,
                        feature_track_ratio_factor_override=self.feature_track_ratio_factor_override,
                        track_preset=self.track_type,
                    )
        return group

    def setup_record_group(self) -> Group:
        selected_features_set: str = self.feature_config.selected_features_set
        color_table: Optional[DataFrame] = self.feature_config.color_table
        default_colors: Optional[DataFrame] = self.feature_config.default_colors
        if self.precomputed_feature_dict is not None:
            feature_dict = self.precomputed_feature_dict
        else:
            compute_label_text = self.show_labels and self.precalculated_labels is None
            label_filtering = (
                preprocess_label_filtering(self.label_filtering)
                if compute_label_text
                else {}
            )
            color_table, default_colors = preprocess_color_tables(color_table, default_colors)
            feature_dict, _ = create_feature_dict(
                self.gb_record,
                color_table,
                selected_features_set,
                default_colors,
                self.strandedness,
                self.resolve_overlaps,
                label_filtering,
                split_overlaps_by_strand=self.split_overlaps_by_strand,
                directional_feature_types=self.feature_config.directional_feature_types,
                feature_visibility_rules=self.feature_config.feature_visibility_rules,
                compute_label_text=compute_label_text,
            )
        track_id: str = self.gb_record.id
        record_group = Group(id=track_id)
        record_length: int = len(self.gb_record.seq)

        record_group: Group = self.draw_record(feature_dict, record_length, record_group)
        return record_group

    def get_group(self) -> Group:
        return self.record_group


__all__ = ["SeqRecordGroup"]


