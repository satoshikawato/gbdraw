#!/usr/bin/env python
# coding: utf-8

from typing import List, Optional, Set, Tuple

from pandas import DataFrame


def prepare_legend_table(
    gc_config,
    skew_config,
    feature_config,
    features_present,
    blast_config=None,
    has_blast: bool = False,
    used_color_rules: Optional[Set[Tuple[str, str]]] = None,
):
    """
    Prepare the legend table for the diagram.

    Args:
        used_color_rules: Optional set of (caption, color) tuples that were actually
            matched during feature creation. If provided, only these rules will be
            included in the legend. If None, all rules for present feature types
            will be included (legacy behavior).
    """
    legend_table = dict()
    color_table: Optional[DataFrame] = feature_config.color_table
    default_colors: DataFrame = feature_config.default_colors
    features_present: List[str] = features_present
    block_stroke_color: str = feature_config.block_stroke_color
    block_stroke_width: float = feature_config.block_stroke_width
    show_gc = gc_config.show_gc
    gc_stroke_color: str = gc_config.stroke_color
    gc_stroke_width: float = gc_config.stroke_width
    gc_high_fill_color: str = gc_config.high_fill_color
    gc_low_fill_color: str = gc_config.low_fill_color
    show_skew = skew_config.show_skew
    skew_high_fill_color: str = skew_config.high_fill_color
    skew_low_fill_color: str = skew_config.low_fill_color
    skew_stroke_color: str = skew_config.stroke_color
    skew_stroke_width: float = skew_config.stroke_width
    dinucleotide = gc_config.dinucleotide
    feature_specific_colors = dict()
    if color_table is not None and not color_table.empty:
        for _, row in color_table.iterrows():
            feature_type = row["feature_type"]
            if feature_type not in feature_specific_colors:
                feature_specific_colors[feature_type] = []
            feature_specific_colors[feature_type].append((row["caption"], row["color"]))
    for selected_feature in features_present:
        if selected_feature in feature_specific_colors.keys():
            has_matching_rules = False
            for entry in feature_specific_colors[selected_feature]:
                specific_caption = entry[0]
                specific_fill_color = entry[1]
                # Only add to legend if this rule was actually used (or if used_color_rules not provided)
                if used_color_rules is None or (specific_caption, specific_fill_color) in used_color_rules:
                    has_matching_rules = True
                    legend_table[specific_caption] = {
                        "type": "solid",
                        "fill": specific_fill_color,
                        "stroke": block_stroke_color,
                        "width": block_stroke_width,
                    }
            # Only add "other X" entry if at least one specific rule was used
            if has_matching_rules:
                if selected_feature == "CDS":
                    new_selected_key_name = "other proteins"
                else:
                    new_selected_key_name = f"other {selected_feature}s"
                feature_fill_color = default_colors[default_colors["feature_type"] == selected_feature][
                    "color"
                ].values[0]
                legend_table[new_selected_key_name] = {
                    "type": "solid",
                    "fill": feature_fill_color,
                    "stroke": block_stroke_color,
                    "width": block_stroke_width,
                }
            elif used_color_rules is None:
                # Legacy behavior: add "other X" even if no rules matched
                if selected_feature == "CDS":
                    new_selected_key_name = "other proteins"
                else:
                    new_selected_key_name = f"other {selected_feature}s"
                feature_fill_color = default_colors[default_colors["feature_type"] == selected_feature][
                    "color"
                ].values[0]
                legend_table[new_selected_key_name] = {
                    "type": "solid",
                    "fill": feature_fill_color,
                    "stroke": block_stroke_color,
                    "width": block_stroke_width,
                }
            else:
                # used_color_rules is provided but no specific rules matched
                # Just add the feature type with default color
                matching_rows = default_colors[default_colors["feature_type"] == selected_feature]
                if not matching_rows.empty:
                    feature_fill_color = matching_rows["color"].values[0]
                else:
                    feature_fill_color = default_colors[default_colors["feature_type"] == "default"]["color"].values[0]
                legend_table[selected_feature] = {
                    "type": "solid",
                    "fill": feature_fill_color,
                    "stroke": block_stroke_color,
                    "width": block_stroke_width,
                }
        else:
            matching_rows = default_colors[default_colors["feature_type"] == selected_feature]
            if not matching_rows.empty:
                feature_fill_color = default_colors[default_colors["feature_type"] == selected_feature][
                    "color"
                ].values[0]
            else:
                feature_fill_color = default_colors[default_colors["feature_type"] == "default"][
                    "color"
                ].values[0]
            legend_table[selected_feature] = {
                "type": "solid",
                "fill": feature_fill_color,
                "stroke": block_stroke_color,
                "width": block_stroke_width,
            }
    if show_gc:
        if gc_high_fill_color == gc_low_fill_color:
            legend_table[f"{dinucleotide} content"] = {
                "type": "solid",
                "fill": gc_high_fill_color,
                "stroke": gc_stroke_color,
                "width": gc_stroke_width,
            }
        else:
            legend_table[f"{dinucleotide} content (+)"] = {
                "type": "solid",
                "fill": gc_high_fill_color,
                "stroke": gc_stroke_color,
                "width": gc_stroke_width,
            }
            legend_table[f"{dinucleotide} content (-)"] = {
                "type": "solid",
                "fill": gc_low_fill_color,
                "stroke": gc_stroke_color,
                "width": gc_stroke_width,
            }
    if show_skew:
        if skew_high_fill_color == skew_low_fill_color:
            legend_table[f"{dinucleotide} skew"] = {
                "type": "solid",
                "fill": skew_high_fill_color,
                "stroke": skew_stroke_color,
                "width": skew_stroke_width,
            }
        else:
            legend_table[f"{dinucleotide} skew (+)"] = {
                "type": "solid",
                "fill": skew_high_fill_color,
                "stroke": skew_stroke_color,
                "width": skew_stroke_width,
            }
            legend_table[f"{dinucleotide} skew (-)"] = {
                "type": "solid",
                "fill": skew_low_fill_color,
                "stroke": skew_stroke_color,
                "width": skew_stroke_width,
            }
    if has_blast and blast_config:
        legend_table["Pairwise match identity"] = {
            "type": "gradient",
            "min_color": blast_config.min_color,
            "max_color": blast_config.max_color,
            "stroke": "none",
            "width": 0,
            "min_value": blast_config.identity,
        }
    return legend_table


__all__ = ["prepare_legend_table"]


