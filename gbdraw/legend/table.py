#!/usr/bin/env python
# coding: utf-8

from typing import Iterable, List, Optional, Sequence, Set, Tuple

from pandas import DataFrame

from gbdraw.analysis.collinearity import normalize_collinearity_color_mode
from gbdraw.core.color import (
    DEFAULT_COLLINEAR_ORIENTATION_COLORS,
    DEFAULT_COLLINEAR_ORIENTATION_MIN_COLORS,
    normalize_hex_color,
)
from gbdraw.exceptions import ValidationError


def _unique_legend_key(legend_table: dict, preferred: str) -> str:
    if preferred not in legend_table:
        return preferred
    suffix = 2
    while f"{preferred} ({suffix})" in legend_table:
        suffix += 1
    return f"{preferred} ({suffix})"


def _normalize_collinearity_legend_color_mode(value: object) -> str | None:
    if not str(value or "").strip():
        return None
    try:
        return normalize_collinearity_color_mode(str(value))
    except ValidationError:
        return None


def collect_comparison_collinearity_color_modes(
    comparison_dataframes: Sequence[DataFrame] | None,
) -> set[str]:
    """Collect recognized collinearity color modes from comparison metadata."""

    color_modes: set[str] = set()
    if not comparison_dataframes:
        return color_modes
    for frame in comparison_dataframes:
        if "collinearity_color_mode" not in frame.columns:
            continue
        for value in frame["collinearity_color_mode"].dropna():
            normalized = _normalize_collinearity_legend_color_mode(value)
            if normalized is not None:
                color_modes.add(normalized)
    return color_modes


def configure_pairwise_identity_legend_from_comparisons(
    blast_config,
    comparison_dataframes: Sequence[DataFrame] | None,
    *,
    additional_color_modes: Iterable[object] | None = None,
) -> set[str]:
    """Configure pairwise identity legend labels from collinearity metadata."""

    color_modes = collect_comparison_collinearity_color_modes(comparison_dataframes)
    for value in additional_color_modes or ():
        normalized = _normalize_collinearity_legend_color_mode(value)
        if normalized is not None:
            color_modes.add(normalized)

    if not color_modes:
        return color_modes

    if hasattr(blast_config, "pairwise_identity_legend_entries"):
        delattr(blast_config, "pairwise_identity_legend_entries")
    if hasattr(blast_config, "pairwise_identity_legend_label"):
        delattr(blast_config, "pairwise_identity_legend_label")

    blast_config.hide_pairwise_identity_legend = color_modes == {"orientation"}
    if "orientation_identity" in color_modes:
        orientation_min_colors = getattr(
            blast_config,
            "collinearity_orientation_min_colors",
            DEFAULT_COLLINEAR_ORIENTATION_MIN_COLORS,
        )
        orientation_colors = getattr(
            blast_config,
            "collinearity_orientation_colors",
            DEFAULT_COLLINEAR_ORIENTATION_COLORS,
        )
        blast_config.pairwise_identity_legend_entries = [
            {
                "label": "Collinear",
                "min_color": orientation_min_colors.get(
                    "plus", DEFAULT_COLLINEAR_ORIENTATION_MIN_COLORS["plus"]
                ),
                "max_color": orientation_colors.get(
                    "plus", DEFAULT_COLLINEAR_ORIENTATION_COLORS["plus"]
                ),
            },
            {
                "label": "Inverted",
                "min_color": orientation_min_colors.get(
                    "minus", DEFAULT_COLLINEAR_ORIENTATION_MIN_COLORS["minus"]
                ),
                "max_color": orientation_colors.get(
                    "minus", DEFAULT_COLLINEAR_ORIENTATION_COLORS["minus"]
                ),
            },
        ]
    elif "average_identity" in color_modes:
        blast_config.pairwise_identity_legend_label = "Average identity"

    return color_modes


def _legend_fill_identity(color: object) -> str:
    """Return the color identity used to compare two solid legend rows."""

    from gbdraw.io.colors import resolve_color_to_hex

    text = str(color or "").strip()
    if text.lower() == "none":
        return "none"
    try:
        return normalize_hex_color(resolve_color_to_hex(text))
    except (ValueError, ValidationError):
        return text.lower()


def _other_feature_legend_key(feature_type: str) -> str:
    return "other proteins" if feature_type == "CDS" else f"other {feature_type}s"


def _resolve_legend_feature_color(default_colors: DataFrame, feature_type: str) -> str:
    matching_rows = default_colors[default_colors["feature_type"] == feature_type]
    if not matching_rows.empty:
        return str(matching_rows["color"].values[0])

    default_rows = default_colors[default_colors["feature_type"] == "default"]
    if not default_rows.empty:
        return str(default_rows["color"].values[0])

    return "#d3d3d3"


def _generated_legend_fills(
    feature_specific_colors: dict,
    features_present: Sequence[str],
    default_colors: DataFrame,
    *,
    used_color_rules: Optional[Set[Tuple[str, str]]],
    default_used_features: Set[str],
    gc_config,
    skew_config,
    depth_config,
    show_gc: bool,
    show_skew: bool,
    show_depth: bool,
) -> dict[str, str]:
    """Return the renderer-generated legend rows and their fill identities."""

    fills: dict[str, str] = {}
    for feature_type in features_present:
        feature_color = _legend_fill_identity(
            _resolve_legend_feature_color(default_colors, feature_type)
        )
        if feature_type not in feature_specific_colors:
            fills[feature_type] = feature_color
        elif used_color_rules is None:
            fills[_other_feature_legend_key(feature_type)] = feature_color
        elif any(
            (caption, color) in used_color_rules
            for caption, color in feature_specific_colors[feature_type]
            if str(caption or "").strip()
        ):
            if feature_type in default_used_features:
                fills[_other_feature_legend_key(feature_type)] = feature_color
        else:
            fills[feature_type] = feature_color
    dinucleotide = gc_config.dinucleotide
    if show_depth and depth_config is not None:
        fills["Depth"] = _legend_fill_identity(depth_config.fill_color)
    for show, config, name in ((show_gc, gc_config, "content"), (show_skew, skew_config, "skew")):
        if not show:
            continue
        if config.high_fill_color == config.low_fill_color:
            fills[f"{dinucleotide} {name}"] = _legend_fill_identity(config.high_fill_color)
        else:
            fills[f"{dinucleotide} {name} (+)"] = _legend_fill_identity(config.high_fill_color)
            fills[f"{dinucleotide} {name} (-)"] = _legend_fill_identity(config.low_fill_color)
    return fills


def prepare_legend_table(
    gc_config,
    skew_config,
    feature_config,
    features_present,
    blast_config=None,
    has_blast: bool = False,
    used_color_rules: Optional[Set[Tuple[str, str]]] = None,
    default_used_features: Optional[Set[str]] = None,
    depth_config=None,
    *,
    show_gc: bool,
    show_skew: bool,
    show_depth: bool,
):
    """
    Prepare the legend table for the diagram.

    Args:
        used_color_rules: Optional set of (caption, color) tuples that were actually
            matched during feature creation. If provided, only these rules will be
            included in the legend. If None, all rules for present feature types
            will be included (legacy behavior).
        default_used_features: Optional set of feature types that fell back to
            default colors. Used to decide whether to include "other X" entries
            when specific rules exist.
    """
    legend_table = dict()
    color_table: Optional[DataFrame] = feature_config.color_table
    default_colors: DataFrame = feature_config.default_colors
    features_present: List[str] = features_present
    block_stroke_color: str = feature_config.block_stroke_color
    block_stroke_width: float = feature_config.block_stroke_width
    gc_stroke_color: str = gc_config.stroke_color
    gc_stroke_width: float = gc_config.stroke_width
    gc_high_fill_color: str = gc_config.high_fill_color
    gc_low_fill_color: str = gc_config.low_fill_color
    skew_high_fill_color: str = skew_config.high_fill_color
    skew_low_fill_color: str = skew_config.low_fill_color
    skew_stroke_color: str = skew_config.stroke_color
    skew_stroke_width: float = skew_config.stroke_width
    dinucleotide = gc_config.dinucleotide
    feature_specific_colors = dict()
    default_used_features = default_used_features or set()
    if color_table is not None and not color_table.empty:
        for _, row in color_table.iterrows():
            feature_type = row["feature_type"]
            if feature_type not in feature_specific_colors:
                feature_specific_colors[feature_type] = []
            feature_specific_colors[feature_type].append((row["caption"], row["color"]))
    # PD-OI-042: the renderer's own row names (feature types, "other ..." and
    # numeric tracks) keep their names. A used rule caption that equals one of
    # them with a different color gets the normalized hex suffix, so both
    # colors keep a solid row; the same color shares one row.
    generated_fills = _generated_legend_fills(
        feature_specific_colors,
        features_present,
        default_colors,
        used_color_rules=used_color_rules,
        default_used_features=default_used_features,
        gc_config=gc_config,
        skew_config=skew_config,
        depth_config=depth_config,
        show_gc=show_gc,
        show_skew=show_skew,
        show_depth=show_depth,
    )
    reserved_captions = dict.fromkeys(
        [
            *generated_fills,
            *(
                str(caption)
                for entries in feature_specific_colors.values()
                for caption, _color in entries
                if str(caption or "").strip()
            ),
        ]
    )
    for selected_feature in features_present:
        if selected_feature in feature_specific_colors.keys():
            has_matching_rules = False
            for entry in feature_specific_colors[selected_feature]:
                specific_caption = entry[0]
                if not str(specific_caption or "").strip():
                    continue
                specific_fill_color = entry[1]
                # Only add to legend if this rule was actually used (or if used_color_rules not provided)
                if used_color_rules is None or (specific_caption, specific_fill_color) in used_color_rules:
                    has_matching_rules = True
                    generated_fill = generated_fills.get(specific_caption)
                    if generated_fill is not None and generated_fill != _legend_fill_identity(
                        specific_fill_color
                    ):
                        specific_caption = _unique_legend_key(
                            reserved_captions,
                            f"{specific_caption} [{_legend_fill_identity(specific_fill_color)}]",
                        )
                        reserved_captions[specific_caption] = None
                    if specific_caption in legend_table:
                        continue
                    legend_table[specific_caption] = {
                        "type": "solid",
                        "fill": specific_fill_color,
                        "stroke": block_stroke_color,
                        "width": block_stroke_width,
                    }
            # Only add "other X" entry if at least one specific rule was used
            if has_matching_rules:
                allow_other_entry = False
                if used_color_rules is None:
                    allow_other_entry = True
                elif selected_feature in default_used_features:
                    allow_other_entry = True
                new_selected_key_name = _other_feature_legend_key(selected_feature)
                if allow_other_entry:
                    feature_fill_color = _resolve_legend_feature_color(default_colors, selected_feature)
                    legend_table[new_selected_key_name] = {
                        "type": "solid",
                        "fill": feature_fill_color,
                        "stroke": block_stroke_color,
                        "width": block_stroke_width,
                    }
            elif used_color_rules is None:
                # Legacy behavior: add "other X" even if no rules matched
                new_selected_key_name = _other_feature_legend_key(selected_feature)
                feature_fill_color = _resolve_legend_feature_color(default_colors, selected_feature)
                legend_table[new_selected_key_name] = {
                    "type": "solid",
                    "fill": feature_fill_color,
                    "stroke": block_stroke_color,
                    "width": block_stroke_width,
                }
            else:
                # used_color_rules is provided but no specific rules matched
                # Just add the feature type with default color
                feature_fill_color = _resolve_legend_feature_color(default_colors, selected_feature)
                legend_table[selected_feature] = {
                    "type": "solid",
                    "fill": feature_fill_color,
                    "stroke": block_stroke_color,
                    "width": block_stroke_width,
                }
        else:
            feature_fill_color = _resolve_legend_feature_color(default_colors, selected_feature)
            legend_table[selected_feature] = {
                "type": "solid",
                "fill": feature_fill_color,
                "stroke": block_stroke_color,
                "width": block_stroke_width,
            }
    if show_depth:
        legend_table["Depth"] = {
            "type": "solid",
            "fill": depth_config.fill_color,
            "stroke": depth_config.stroke_color,
            "width": depth_config.stroke_width,
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
    if has_blast and blast_config and not bool(getattr(blast_config, "hide_pairwise_identity_legend", False)):
        identity_legend_entries = getattr(blast_config, "pairwise_identity_legend_entries", None)
        if identity_legend_entries:
            for entry in identity_legend_entries:
                label = str(entry["label"])
                min_color = entry["min_color"]
                if float(blast_config.identity) >= 100:
                    min_color = entry["max_color"]
                legend_table[label] = {
                    "type": "gradient",
                    "min_color": min_color,
                    "max_color": entry["max_color"],
                    "stroke": "none",
                    "width": 0,
                    "min_value": blast_config.identity,
                }
        else:
            identity_legend_label = str(
                getattr(blast_config, "pairwise_identity_legend_label", "Pairwise match identity")
            )
            min_color = blast_config.min_color
            if float(blast_config.identity) >= 100:
                min_color = blast_config.max_color
            legend_table[identity_legend_label] = {
                "type": "gradient",
                "min_color": min_color,
                "max_color": blast_config.max_color,
                "stroke": "none",
                "width": 0,
                "min_value": blast_config.identity,
            }
    return legend_table


__all__ = [
    "collect_comparison_collinearity_color_modes",
    "configure_pairwise_identity_legend_from_comparisons",
    "prepare_legend_table",
]

