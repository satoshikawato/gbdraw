#!/usr/bin/env python
# coding: utf-8

"""Pre-calculation helpers for linear diagram assembly.

These functions compute layout constraints (definition widths, label heights) before
final canvas sizing.
"""

from __future__ import annotations

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]

from ...canvas import LinearCanvasConfigurator  # type: ignore[reportMissingImports]
from ...config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from ...configurators import FeatureDrawingConfigurator  # type: ignore[reportMissingImports]
from ...features.colors import preprocess_color_tables  # type: ignore[reportMissingImports]
from ...features.factory import create_feature_dict  # type: ignore[reportMissingImports]
from ...render.groups.linear import DefinitionGroup  # type: ignore[reportMissingImports]
from ...labels.filtering import preprocess_label_filtering  # type: ignore[reportMissingImports]
from ...labels.placement import prepare_label_list_linear  # type: ignore[reportMissingImports]


def _precalculate_definition_widths(
    records: list[SeqRecord],
    config_dict: dict,
    canvas_config,
    cfg: GbdrawConfig | None = None,
) -> float:
    """
    Pre-calculates the maximum definition width among all records.
    """
    max_definition_width = 0
    if not records:
        return 0

    cfg = cfg or GbdrawConfig.from_dict(config_dict)
    for record in records:
        def_group = DefinitionGroup(record, config_dict, canvas_config, cfg=cfg)
        if def_group.definition_bounding_box_width > max_definition_width:
            max_definition_width = def_group.definition_bounding_box_width
    return max_definition_width


def _precalculate_label_dimensions(
    records: list[SeqRecord],
    feature_config: FeatureDrawingConfigurator,
    canvas_config: LinearCanvasConfigurator,
    config_dict: dict,
    cfg: GbdrawConfig | None = None,
    *,
    feature_dicts: dict[str, dict] | None = None,
) -> tuple[float, dict, dict]:
    """Pre-calculates label placements for all records to determine the required canvas height."""

    cfg = cfg or GbdrawConfig.from_dict(config_dict)
    raw_show_labels = cfg.canvas.show_labels
    show_labels_mode = raw_show_labels if isinstance(raw_show_labels, str) else ("all" if raw_show_labels else "none")

    if show_labels_mode == "none":
        return 0, {}, {}

    max_required_height = 0
    all_labels_by_record = {}
    record_label_heights = {}  # Store height required for labels per record
    normalize_length = cfg.canvas.linear.normalize_length

    for i, record in enumerate(records):
        if show_labels_mode == "first" and i > 0:
            all_labels_by_record[record.id] = []
            record_label_heights[record.id] = 0
            continue

        feature_dict = None
        if feature_dicts is not None:
            feature_dict = feature_dicts.get(record.id)
        if feature_dict is None:
            color_table, default_colors = preprocess_color_tables(
                feature_config.color_table, feature_config.default_colors
            )
            label_filtering = preprocess_label_filtering(cfg.labels.filtering.as_dict())
            feature_dict, _ = create_feature_dict(
                record,
                color_table,
                feature_config.selected_features_set,
                default_colors,
                canvas_config.strandedness,
                canvas_config.resolve_overlaps,
                label_filtering,
            )

        record_length = len(record.seq)
        if normalize_length:
            genome_size_normalization_factor = 1.0
        else:
            genome_size_normalization_factor = record_length / canvas_config.longest_genome

        label_list = prepare_label_list_linear(
            feature_dict,
            record_length,
            canvas_config.alignment_width,
            genome_size_normalization_factor,
            canvas_config.cds_height,
            canvas_config.strandedness,
            config_dict,
            cfg=cfg,
        )
        all_labels_by_record[record.id] = label_list

        min_y_coord_record = 0
        for label in label_list:
            if not label["is_embedded"]:
                label_bottom_y = label["middle_y"] - (label["height_px"] / 2)
                if label_bottom_y < min_y_coord_record:
                    min_y_coord_record = label_bottom_y

        record_height = abs(min_y_coord_record)
        record_label_heights[record.id] = record_height

        if record_height > max_required_height:
            max_required_height = record_height

    return max_required_height, all_labels_by_record, record_label_heights


__all__ = [
    "_precalculate_definition_widths",
    "_precalculate_label_dimensions",
]


