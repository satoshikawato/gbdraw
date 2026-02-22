#!/usr/bin/env python
# coding: utf-8

from typing import Optional, List, Dict

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
        self.track_type = cfg.canvas.circular.track_type
        self.strandedness = cfg.canvas.strandedness
        self.resolve_overlaps = cfg.canvas.resolve_overlaps
        self.split_overlaps_by_strand = (
            bool(self.resolve_overlaps)
            and (not bool(self.strandedness))
            and str(self.track_type).strip().lower() == "middle"
        )
        self.track_ratio_factors = cfg.canvas.circular.track_ratio_factors[self.length_param]
        self.track_ratio = self.canvas_config.track_ratio
        self.precomputed_feature_dict: Optional[Dict[str, FeatureObject]] = precomputed_feature_dict
        self.precalculated_labels: Optional[list[dict]] = precalculated_labels
        self.feature_track_ratio_factor_override = feature_track_ratio_factor_override
        self.record_group: Group = self.setup_record_group()

    def draw_record(self, feature_dict: Dict[str, FeatureObject], record_length: int, group: Group) -> Group:
        label_list = []
        if self.show_labels is True:
            # Reuse pre-calculated labels when available to avoid repeating heavy placement work.
            if self.precalculated_labels is not None:
                label_list = self.precalculated_labels
            else:
                label_list = prepare_label_list(
                    feature_dict,
                    record_length,
                    self.canvas_config.radius,
                    self.track_ratio,
                    self.config_dict,
                    cfg=self._cfg,
                    feature_track_ratio_factor_override=self.feature_track_ratio_factor_override,
                )
        feature_track_ratio_factor = (
            float(self.feature_track_ratio_factor_override)
            if self.feature_track_ratio_factor_override is not None
            else float(self.track_ratio_factors[0])
        )
        for feature_object in feature_dict.values():
            group = FeatureDrawer(self.feature_config).draw(
                feature_object,
                group,
                record_length,
                self.canvas_config.radius,
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
                        self.canvas_config.radius,
                        self.canvas_config.track_ratio,
                        feature_track_ratio_factor_override=self.feature_track_ratio_factor_override,
                    )
        return group

    def setup_record_group(self) -> Group:
        selected_features_set: str = self.feature_config.selected_features_set
        color_table: Optional[DataFrame] = self.feature_config.color_table
        default_colors: Optional[DataFrame] = self.feature_config.default_colors
        if self.precomputed_feature_dict is not None:
            feature_dict = self.precomputed_feature_dict
        else:
            label_filtering = preprocess_label_filtering(self.label_filtering)
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
            )
        track_id: str = self.gb_record.id
        record_group = Group(id=track_id)
        record_length: int = len(self.gb_record.seq)

        record_group: Group = self.draw_record(feature_dict, record_length, record_group)
        return record_group

    def get_group(self) -> Group:
        return self.record_group


__all__ = ["SeqRecordGroup"]


