#!/usr/bin/env python
# coding: utf-8

"""Pre-calculation helpers for linear diagram assembly.

These functions compute layout constraints (definition widths, label heights) before
final canvas sizing.
"""

from __future__ import annotations

import math

from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]
from ...canvas import LinearCanvasConfigurator  # type: ignore[reportMissingImports]
from ...config.models import GbdrawConfig  # type: ignore[reportMissingImports]
from ...configurators import FeatureDrawingConfigurator  # type: ignore[reportMissingImports]
from ...core.sequence import determine_length_parameter  # type: ignore[reportMissingImports]
from ...features.colors import preprocess_color_tables  # type: ignore[reportMissingImports]
from ...features.factory import create_feature_dict  # type: ignore[reportMissingImports]
from ...features.objects import FeatureObject  # type: ignore[reportMissingImports]
from ...render.groups.linear import DefinitionGroup  # type: ignore[reportMissingImports]
from ...labels.filtering import preprocess_label_filtering  # type: ignore[reportMissingImports]
from ...labels.linear import calculate_label_y_bounds, prepare_label_list_linear  # type: ignore[reportMissingImports]
from .orthogroup_alignment import OrthogroupLabelEligibility, orthogroup_label_sets_for_record


FeatureDict = dict[str, FeatureObject]


def _resolve_linear_diagram_label_font_size(
    records: list[SeqRecord],
    *,
    show_labels_mode: str,
    canvas_config: LinearCanvasConfigurator,
    cfg: GbdrawConfig,
) -> float:
    """Resolve one readable auto label size for all labels drawn in a linear diagram."""

    label_font_sizes: list[float] = []
    threshold = int(cfg.labels.length_threshold.linear)
    for index, record in enumerate(records):
        if show_labels_mode == "first" and index > 0:
            continue
        length_param = determine_length_parameter(len(record.seq), threshold)
        label_font_sizes.append(cfg.labels.font_size.linear.for_length_param(length_param))

    if label_font_sizes:
        return max(label_font_sizes)
    return cfg.labels.font_size.linear.for_length_param(canvas_config.length_param)


def _precalculate_definition_metrics(
    records: list[SeqRecord],
    config_dict: dict,
    canvas_config,
    cfg: GbdrawConfig | None = None,
) -> tuple[float, list[float], list[float]]:
    """
    Pre-calculate definition widths and heights for all records.
    """
    max_definition_width = 0
    definition_heights: list[float] = []
    definition_half_heights: list[float] = []
    if not records:
        return 0.0, definition_heights, definition_half_heights

    cfg = cfg or GbdrawConfig.from_dict(config_dict)
    for record in records:
        def_group = DefinitionGroup(record, config_dict, canvas_config, cfg=cfg)
        if def_group.definition_bounding_box_width > max_definition_width:
            max_definition_width = def_group.definition_bounding_box_width
        definition_height = float(def_group.definition_bounding_box_height)
        definition_heights.append(definition_height)
        definition_half_heights.append(0.5 * definition_height)
    return math.ceil(max_definition_width), definition_heights, definition_half_heights


def _precalculate_definition_widths(
    records: list[SeqRecord],
    config_dict: dict,
    canvas_config,
    cfg: GbdrawConfig | None = None,
) -> float:
    max_definition_width, _definition_heights, _definition_half_heights = _precalculate_definition_metrics(
        records,
        config_dict,
        canvas_config,
        cfg=cfg,
    )
    return max_definition_width


def _precalculate_feature_dicts(
    records: list[SeqRecord],
    feature_config: FeatureDrawingConfigurator,
    canvas_config: LinearCanvasConfigurator,
    config_dict: dict,
    cfg: GbdrawConfig | None = None,
    orthogroup_label_eligibility: OrthogroupLabelEligibility | None = None,
) -> list[FeatureDict]:
    """Build feature objects once per record for the linear assembly pipeline."""

    cfg = cfg or GbdrawConfig.from_dict(config_dict)
    raw_show_labels = cfg.canvas.show_labels
    show_labels_mode = raw_show_labels if isinstance(raw_show_labels, str) else ("all" if raw_show_labels else "none")
    color_table, default_colors = preprocess_color_tables(
        feature_config.color_table,
        feature_config.default_colors,
    )
    label_filtering = preprocess_label_filtering(cfg.labels.filtering.as_dict())

    feature_dicts: list[FeatureDict] = []
    for i, record in enumerate(records):
        compute_label_text = (
            show_labels_mode == "all"
            or (show_labels_mode == "first" and i == 0)
            or show_labels_mode == "orthogroup_top"
        )
        feature_dict, _ = create_feature_dict(
            record,
            color_table,
            feature_config.selected_features_set,
            default_colors,
            canvas_config.strandedness,
            canvas_config.resolve_overlaps,
            label_filtering if compute_label_text else {},
            directional_feature_types=feature_config.directional_feature_types,
            feature_visibility_rules=feature_config.feature_visibility_rules,
            compute_label_text=compute_label_text,
        )
        feature_dicts.append(feature_dict)
    return feature_dicts


def _precalculate_label_dimensions(
    records: list[SeqRecord],
    feature_config: FeatureDrawingConfigurator,
    canvas_config: LinearCanvasConfigurator,
    config_dict: dict,
    cfg: GbdrawConfig | None = None,
    precomputed_feature_dicts: list[FeatureDict] | None = None,
    orthogroup_label_eligibility: OrthogroupLabelEligibility | None = None,
) -> tuple[float, list[list[dict]], list[float]]:
    """Pre-calculates label placements for all records to determine the required canvas height."""

    cfg = cfg or GbdrawConfig.from_dict(config_dict)
    raw_show_labels = cfg.canvas.show_labels
    show_labels_mode = raw_show_labels if isinstance(raw_show_labels, str) else ("all" if raw_show_labels else "none")

    if show_labels_mode == "none":
        return 0, [[] for _ in records], [0.0 for _ in records]

    label_font_size = _resolve_linear_diagram_label_font_size(
        records,
        show_labels_mode=show_labels_mode,
        canvas_config=canvas_config,
        cfg=cfg,
    )
    max_required_height = 0
    all_labels_by_record: list[list[dict]] = []
    record_label_heights: list[float] = []
    normalize_length = cfg.canvas.linear.normalize_length
    if precomputed_feature_dicts is None:
        color_table, default_colors = preprocess_color_tables(
            feature_config.color_table,
            feature_config.default_colors,
        )
        label_filtering = preprocess_label_filtering(cfg.labels.filtering.as_dict())

    for i, record in enumerate(records):
        if show_labels_mode == "first" and i > 0:
            all_labels_by_record.append([])
            record_label_heights.append(0.0)
            continue
        member_ids, top_member_ids = orthogroup_label_sets_for_record(
            orthogroup_label_eligibility if show_labels_mode == "orthogroup_top" else None,
            i,
        )

        if precomputed_feature_dicts is not None:
            feature_dict = precomputed_feature_dicts[i]
        else:
            feature_dict, _ = create_feature_dict(
                record,
                color_table,
                feature_config.selected_features_set,
                default_colors,
                canvas_config.strandedness,
                canvas_config.resolve_overlaps,
                label_filtering,
                directional_feature_types=feature_config.directional_feature_types,
                feature_visibility_rules=feature_config.feature_visibility_rules,
            )

        record_length = len(record.seq)
        if normalize_length:
            genome_size_normalization_factor = 1.0
        else:
            genome_size_normalization_factor = record_length / canvas_config.longest_genome

        # Label x-positions must use the same final axis width as record rendering.
        label_list = prepare_label_list_linear(
            feature_dict,
            record_length,
            canvas_config.alignment_width,
            genome_size_normalization_factor,
            canvas_config.cds_height,
            canvas_config.strandedness,
            canvas_config.track_layout,
            canvas_config.track_axis_gap,
            config_dict,
            cfg=cfg,
            label_font_size=label_font_size,
            orthogroup_label_member_ids=member_ids,
            orthogroup_label_top_member_ids=top_member_ids,
        )
        all_labels_by_record.append(label_list)

        min_y_coord_record = 0.0
        for label in label_list:
            label_top_y, _ = calculate_label_y_bounds(label)
            is_above_feature_embedded = (
                bool(label.get("is_embedded"))
                and label_top_y < float(label.get("feature_top_y", 0.0))
            )
            if bool(label.get("is_embedded")) and not is_above_feature_embedded:
                continue
            if label_top_y < min_y_coord_record:
                min_y_coord_record = label_top_y

        record_height = abs(min_y_coord_record)
        record_label_heights.append(record_height)

        if record_height > max_required_height:
            max_required_height = record_height

    return max_required_height, all_labels_by_record, record_label_heights


__all__ = [
    "FeatureDict",
    "_precalculate_definition_metrics",
    "_precalculate_definition_widths",
    "_precalculate_feature_dicts",
    "_precalculate_label_dimensions",
    "_resolve_linear_diagram_label_font_size",
]


