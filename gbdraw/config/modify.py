#!/usr/bin/env python
# coding: utf-8

from .models import GbdrawConfig  # type: ignore[reportMissingImports]

def suppress_gc_content_and_skew(suppress_gc: bool, suppress_skew: bool) -> tuple[bool, bool]:
    show_gc, show_skew = True, True
    if suppress_gc is True:
        show_gc = False
    if suppress_skew is True:
        show_skew = False
    return show_gc, show_skew


def update_config_value(config_dict, path, value):
    keys = path.split(".")
    for key in keys[:-1]:
        config_dict = config_dict.setdefault(key, {})
    config_dict[keys[-1]] = value


def modify_config_dict(
    config_dict,
    block_stroke_width=None,
    block_stroke_color=None,
    circular_axis_stroke_color=None,
    circular_axis_stroke_width=None,
    linear_axis_stroke_color=None,
    linear_axis_stroke_width=None,
    line_stroke_color=None,
    line_stroke_width=None,
    gc_stroke_color=None,
    linear_definition_font_size=None,
    circular_definition_font_size=None,
    label_font_size=None,
    label_placement=None,
    label_rotation=None,
    show_gc=None,
    show_skew=None,
    show_labels=None,
    align_center=None,
    linear_track_layout=None,
    linear_track_axis_gap=None,
    cicular_width_with_labels=None,
    track_type=None,
    strandedness=None,
    resolve_overlaps=None,
    allow_inner_labels=None,
    label_radius_offset=None,
    label_blacklist=None,
    label_whitelist=None,
    qualifier_priority=None,
    outer_label_x_radius_offset=None,
    outer_label_y_radius_offset=None,
    inner_label_x_radius_offset=None,
    inner_label_y_radius_offset=None,
    comparison_height=None,
    font_family=None,
    default_cds_height=None,
    gc_height=None,
    scale_style=None,
    scale_stroke_color=None,
    scale_stroke_width=None,
    scale_font_size=None,
    scale_interval=None,
    blast_color_min=None,
    blast_color_max=None,
    legend_box_size=None,
    legend_font_size=None,
    normalize_length=None,
) -> dict:

    cfg = GbdrawConfig.from_dict(config_dict)

    label_font_size_circular_long = (
        label_font_size if label_font_size is not None else cfg.labels.font_size.long
    )
    label_font_size_circular_short = (
        label_font_size if label_font_size is not None else cfg.labels.font_size.short
    )
    label_font_size_linear_long = (
        label_font_size if label_font_size is not None else cfg.labels.font_size.linear.long
    )
    label_font_size_linear_short = (
        label_font_size if label_font_size is not None else cfg.labels.font_size.linear.short
    )

    if linear_definition_font_size is not None:
        linear_definition_font_size_short = linear_definition_font_size
        linear_definition_font_size_long = linear_definition_font_size
    else:
        linear_definition_font_size_short = cfg.objects.definition.linear.font_size.short
        linear_definition_font_size_long = cfg.objects.definition.linear.font_size.long

    if default_cds_height is not None:
        default_cds_height_short = default_cds_height
        default_cds_height_long = default_cds_height
    else:
        default_cds_height_short = cfg.canvas.linear.default_cds_height.short
        default_cds_height_long = cfg.canvas.linear.default_cds_height.long

    if circular_definition_font_size is not None:
        circular_definition_font_interval = float(circular_definition_font_size) + 2
    else:
        circular_definition_font_interval = None

    if legend_box_size is not None:
        legend_box_size_short = legend_box_size
        legend_box_size_long = legend_box_size
    else:
        legend_box_size_short = cfg.objects.legends.color_rect_size.short
        legend_box_size_long = cfg.objects.legends.color_rect_size.long

    if legend_font_size is not None:
        legend_font_size_short = legend_font_size
        legend_font_size_long = legend_font_size
    else:
        legend_font_size_short = cfg.objects.legends.font_size.short
        legend_font_size_long = cfg.objects.legends.font_size.long

    if scale_font_size is not None:
        scale_font_size_short = scale_font_size
        scale_font_size_long = scale_font_size
    else:
        scale_font_size_short = cfg.objects.scale.font_size.short
        scale_font_size_long = cfg.objects.scale.font_size.long

    if linear_axis_stroke_width is not None:
        linear_axis_stroke_width_short = linear_axis_stroke_width
        linear_axis_stroke_width_long = linear_axis_stroke_width
    else:
        linear_axis_stroke_width_short = cfg.objects.axis.linear.stroke_width.short
        linear_axis_stroke_width_long = cfg.objects.axis.linear.stroke_width.long

    if line_stroke_width is not None:
        line_stroke_width_short = line_stroke_width
        line_stroke_width_long = line_stroke_width
    else:
        line_stroke_width_short = cfg.objects.features.line_stroke_width.short
        line_stroke_width_long = cfg.objects.features.line_stroke_width.long

    if block_stroke_width is not None:
        block_stroke_width_short = block_stroke_width
        block_stroke_width_long = block_stroke_width
    else:
        block_stroke_width_short = cfg.objects.features.block_stroke_width.short
        block_stroke_width_long = cfg.objects.features.block_stroke_width.long

    if circular_axis_stroke_width is not None:
        circular_axis_stroke_width_short = circular_axis_stroke_width
        circular_axis_stroke_width_long = circular_axis_stroke_width
    else:
        circular_axis_stroke_width_short = cfg.objects.axis.circular.stroke_width.short
        circular_axis_stroke_width_long = cfg.objects.axis.circular.stroke_width.long

    mapping = {
        "block_stroke_width_short": "objects.features.block_stroke_width.short",
        "block_stroke_width_long": "objects.features.block_stroke_width.long",
        "block_stroke_color": "objects.features.block_stroke_color",
        "circular_axis_stroke_color": "objects.axis.circular.stroke_color",
        "circular_axis_stroke_width_short": "objects.axis.circular.stroke_width.short",
        "circular_axis_stroke_width_long": "objects.axis.circular.stroke_width.long",
        "linear_axis_stroke_color": "objects.axis.linear.stroke_color",
        "linear_axis_stroke_width_short": "objects.axis.linear.stroke_width.short",
        "linear_axis_stroke_width_long": "objects.axis.linear.stroke_width.long",
        "line_stroke_color": "objects.features.line_stroke_color",
        "line_stroke_width_short": "objects.features.line_stroke_width.short",
        "line_stroke_width_long": "objects.features.line_stroke_width.long",
        "gc_stroke_color": "objects.gc_content.stroke_color",
        "show_gc": "canvas.show_gc",
        "show_skew": "canvas.show_skew",
        "show_labels": "canvas.show_labels",
        "align_center": "canvas.linear.align_center",
        "linear_track_layout": "canvas.linear.track_layout",
        "linear_track_axis_gap": "canvas.linear.track_axis_gap",
        "cicular_width_with_labels": "canvas.circular.width.with_labels",
        "track_type": "canvas.circular.track_type",
        "strandedness": "canvas.strandedness",
        "resolve_overlaps": "canvas.resolve_overlaps",
        "allow_inner_labels": "canvas.circular.allow_inner_labels",
        "label_radius_offset": "labels.radius_factor",
        "label_blacklist": "labels.filtering.blacklist_keywords",
        "label_whitelist": "labels.filtering.whitelist_df",
        "qualifier_priority": "labels.filtering.qualifier_priority_df",
        "outer_label_x_radius_offset": "labels.unified_adjustment.outer_labels.x_radius_offset",
        "outer_label_y_radius_offset": "labels.unified_adjustment.outer_labels.y_radius_offset",
        "inner_label_x_radius_offset": "labels.unified_adjustment.inner_labels.x_radius_offset",
        "inner_label_y_radius_offset": "labels.unified_adjustment.inner_labels.y_radius_offset",
        "comparison_height": "canvas.linear.comparison_height",
        "font_family": "objects.text.font_family",
        "default_cds_height_short": "canvas.linear.default_cds_height.short",
        "default_cds_height_long": "canvas.linear.default_cds_height.long",
        "gc_height": "canvas.linear.default_gc_height",
        "scale_style": "objects.scale.style",
        "scale_stroke_color": "objects.scale.stroke_color",
        "scale_stroke_width": "objects.scale.stroke_width",
        "scale_font_size_short": "objects.scale.font_size.short",
        "scale_font_size_long": "objects.scale.font_size.long",
        "scale_interval": "objects.scale.interval",
        "blast_color_min": "objects.blast_match.min_color",
        "blast_color_max": "objects.blast_match.max_color",
        "legend_box_size_short": "objects.legends.color_rect_size.short",
        "legend_box_size_long": "objects.legends.color_rect_size.long",
        "legend_font_size_short": "objects.legends.font_size.short",
        "legend_font_size_long": "objects.legends.font_size.long",
        "label_font_size_circular_short": "labels.font_size.short",
        "label_font_size_circular_long": "labels.font_size.long",
        "label_font_size_linear_short": "labels.font_size.linear.short",
        "label_font_size_linear_long": "labels.font_size.linear.long",
        "label_placement": "labels.linear.placement",
        "label_rotation": "labels.linear.rotation",
        "linear_definition_font_size_short": "objects.definition.linear.font_size.short",
        "linear_definition_font_size_long": "objects.definition.linear.font_size.long",
        "circular_definition_font_size": "objects.definition.circular.font_size",
        "circular_definition_font_interval": "objects.definition.circular.interval",
        "normalize_length": "canvas.linear.normalize_length",
    }

    for param, path in mapping.items():
        value = locals()[param]
        if value is not None:
            # Convert label_blacklist string to list if needed
            if param == "label_blacklist" and isinstance(value, str):
                value = [kw.strip() for kw in value.split(",") if kw.strip()]
            # Skip updating whitelist_df and qualifier_priority_df if they're already DataFrames
            # (they were processed from file paths before this function was called)
            # Also skip if the value is an empty string and existing value is None or DataFrame
            if param in ["label_whitelist", "qualifier_priority"]:
                keys = path.split(".")
                target_dict = config_dict
                for key in keys[:-1]:
                    target_dict = target_dict.setdefault(key, {})
                existing_value = target_dict.get(keys[-1])
                # If existing value is a DataFrame, don't overwrite it with a string
                if existing_value is not None and hasattr(existing_value, 'empty'):
                    continue
                # If value is empty string and existing is None, skip (empty string means "no file")
                if value == "" and existing_value is None:
                    continue
            update_config_value(config_dict, path, value)

    return config_dict


__all__ = ["modify_config_dict", "suppress_gc_content_and_skew", "update_config_value"]


