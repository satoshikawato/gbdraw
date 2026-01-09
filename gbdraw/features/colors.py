#!/usr/bin/env python
# coding: utf-8

import re
from pandas import DataFrame
from Bio.SeqFeature import SeqFeature


def preprocess_color_tables(color_table: DataFrame, default_colors: DataFrame) -> tuple[dict, dict]:
    """
    Preprocesses color tables to create mappings for feature coloring.

    Returns:
        (color_map, default_color_map)
    """
    # Create a mapping for default colors
    default_color_map = default_colors.set_index("feature_type")["color"].to_dict()

    # Create a nested dictionary for specific color rules
    # feature_type -> qualifier_key -> list of (pattern, color)
    color_map: dict = {}
    if isinstance(color_table, DataFrame) and not color_table.empty:
        for row in color_table.itertuples(index=False):
            # Compile regex pattern for case-insensitive matching
            pattern = re.compile(row.value, re.IGNORECASE)

            # Build nested dictionary structure
            qualifier_rules = color_map.setdefault(row.feature_type, {})
            rule_list = qualifier_rules.setdefault(row.qualifier_key, [])

            rule_list.append((pattern, row.color))

    return color_map, default_color_map


def get_color(feature: SeqFeature, color_map: dict, default_color_map: dict) -> str:
    """
    Determines the color for a given feature based on its type and qualifiers.

    Supports special 'location' pseudo-qualifier for position-based matching.
    This allows targeting individual features without locus_tag by their
    genomic position (e.g., "1000..2000" format).
    """
    # Check for specific rules first
    rules_for_feature_type = color_map.get(feature.type)
    if rules_for_feature_type:
        # First, check regular qualifiers
        for qualifier_key, qualifier_values in feature.qualifiers.items():
            if qualifier_key in rules_for_feature_type:
                # Check each pattern for this qualifier
                for pattern, color in rules_for_feature_type[qualifier_key]:
                    for value in qualifier_values:
                        if pattern.search(value):
                            return color

        # Check 'location' pseudo-qualifier for position-based matching
        # This allows targeting features without locus_tag (misc_feature, repeat_region, etc.)
        if 'location' in rules_for_feature_type:
            # Create location string in GenBank format: "start..end"
            location_str = f"{int(feature.location.start)}..{int(feature.location.end)}"
            for pattern, color in rules_for_feature_type['location']:
                if pattern.search(location_str):
                    return color

    # Fallback to default color if no specific rule matched
    return default_color_map.get(feature.type, "#d3d3d3")  # Fallback color


__all__ = ["get_color", "preprocess_color_tables"]


