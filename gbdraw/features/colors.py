#!/usr/bin/env python
# coding: utf-8

import hashlib
import re
from typing import Optional, Tuple

from pandas import DataFrame
from Bio.SeqFeature import SeqFeature


def compute_feature_hash(feature: SeqFeature) -> str:
    """
    Compute a stable hash for a feature based on type, start, end, strand.
    This matches the hash computed in render/drawers for SVG IDs.
    For multi-exon features (CompoundLocation), uses first part's coordinates.
    """
    # Use first part's coordinates for CompoundLocation (matches drawer code)
    loc = feature.location
    if hasattr(loc, 'parts') and loc.parts:
        first_part = loc.parts[0]
        start = int(first_part.start)
        end = int(first_part.end)
        strand = first_part.strand
    else:
        start = int(loc.start)
        end = int(loc.end)
        strand = loc.strand
    key = f"{feature.type}:{start}:{end}:{strand}"
    return "f" + hashlib.md5(key.encode()).hexdigest()[:8]


def preprocess_color_tables(color_table: DataFrame, default_colors: DataFrame) -> tuple[dict, dict]:
    """
    Preprocesses color tables to create mappings for feature coloring.

    Returns:
        (color_map, default_color_map)
    """
    # Create a mapping for default colors
    default_color_map = default_colors.set_index("feature_type")["color"].to_dict()

    # Create a nested dictionary for specific color rules
    # feature_type -> qualifier_key -> list of (pattern, color, caption)
    color_map: dict = {}
    if isinstance(color_table, DataFrame) and not color_table.empty:
        for row in color_table.itertuples(index=False):
            # Compile regex pattern for case-insensitive matching
            pattern = re.compile(row.value, re.IGNORECASE)

            # Build nested dictionary structure
            qualifier_rules = color_map.setdefault(row.feature_type, {})
            rule_list = qualifier_rules.setdefault(row.qualifier_key, [])

            # Include caption for legend tracking
            caption = getattr(row, 'caption', '') or ''
            rule_list.append((pattern, row.color, caption))

    return color_map, default_color_map


def get_color_with_info(feature: SeqFeature, color_map: dict, default_color_map: dict) -> Tuple[str, Optional[str]]:
    """
    Determines the color for a given feature based on its type and qualifiers.

    Supports special pseudo-qualifiers for targeting features:
    - 'hash': Match by feature hash (most reliable, based on type+position+strand)
    - 'location': Match by position string "start..end"

    Returns:
        Tuple of (color, caption). Caption is None if default color was used.
    """
    # Check for specific rules first
    rules_for_feature_type = color_map.get(feature.type)
    if rules_for_feature_type:
        # First, check hash pseudo-qualifier (highest priority for individual targeting)
        if 'hash' in rules_for_feature_type:
            feature_hash = compute_feature_hash(feature)
            for rule in rules_for_feature_type['hash']:
                if len(rule) >= 3:
                    pattern, color, caption = rule[0], rule[1], rule[2]
                else:
                    pattern, color = rule[0], rule[1]
                    caption = None
                if pattern.search(feature_hash):
                    return color, caption

        # Then, check regular qualifiers
        for qualifier_key, qualifier_values in feature.qualifiers.items():
            if qualifier_key in rules_for_feature_type:
                # Check each pattern for this qualifier
                for rule in rules_for_feature_type[qualifier_key]:
                    # Support both old (pattern, color) and new (pattern, color, caption) formats
                    if len(rule) >= 3:
                        pattern, color, caption = rule[0], rule[1], rule[2]
                    else:
                        pattern, color = rule[0], rule[1]
                        caption = None
                    for value in qualifier_values:
                        if pattern.search(value):
                            return color, caption

        # Then, check location pseudo-qualifier (for features without unique qualifiers)
        if 'location' in rules_for_feature_type:
            # Create location string: "start..end" (GenBank format)
            location_str = f"{int(feature.location.start)}..{int(feature.location.end)}"
            for rule in rules_for_feature_type['location']:
                if len(rule) >= 3:
                    pattern, color, caption = rule[0], rule[1], rule[2]
                else:
                    pattern, color = rule[0], rule[1]
                    caption = None
                if pattern.search(location_str):
                    return color, caption

    # Fallback to default color if no specific rule matched
    return default_color_map.get(feature.type, "#d3d3d3"), None


def get_color(feature: SeqFeature, color_map: dict, default_color_map: dict) -> str:
    """
    Determines the color for a given feature based on its type and qualifiers.

    This is a convenience wrapper that returns only the color.
    Use get_color_with_info() if you also need the matched caption.
    """
    color, _ = get_color_with_info(feature, color_map, default_color_map)
    return color


def precompute_used_color_rules(
    records,
    color_map: dict,
    default_color_map: dict,
    selected_features_set: set,
) -> set:
    """
    Pre-compute which color rules will be used for a set of records.

    This is useful for generating accurate legends before rendering features.

    Args:
        records: SeqRecord or list of SeqRecords
        color_map: Preprocessed color map from preprocess_color_tables
        default_color_map: Preprocessed default color map
        selected_features_set: Set of feature types to consider

    Returns:
        Set of (caption, color) tuples for rules that match features
    """
    from Bio.SeqRecord import SeqRecord
    if isinstance(records, SeqRecord):
        records = [records]

    used_rules: set = set()
    for record in records:
        for feature in record.features:
            if feature.type not in selected_features_set:
                continue
            color, caption = get_color_with_info(feature, color_map, default_color_map)
            if caption:
                used_rules.add((caption, color))
    return used_rules


__all__ = ["get_color", "get_color_with_info", "precompute_used_color_rules", "preprocess_color_tables"]


