#!/usr/bin/env python
# coding: utf-8

import re
from typing import Any, Optional, Tuple

from pandas import DataFrame
from Bio.SeqFeature import SeqFeature

from .ids import compute_feature_hash as compute_feature_hash
from .selector_values import (
    get_feature_hash,
    get_feature_location_str,
    get_feature_qualifiers,
    get_feature_record_location_str,
    get_qualifier_values,
)


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


def _iter_color_rule_sets(color_map: dict, feature_type: str):
    rules_for_feature_type = color_map.get(feature_type)
    if rules_for_feature_type:
        yield rules_for_feature_type
    wildcard_rules = color_map.get("*")
    if wildcard_rules:
        yield wildcard_rules


def _unpack_color_rule(rule) -> tuple[Any, str, Optional[str]]:
    if len(rule) >= 3:
        return rule[0], rule[1], rule[2]
    return rule[0], rule[1], None


def find_specific_color_rule(
    feature: SeqFeature,
    color_map: dict,
    record_id: Optional[str] = None,
) -> Tuple[str, Optional[str]] | None:
    """Return the first matching specific color rule, if any."""

    feature_type = str(
        getattr(feature, "type", None)
        or getattr(feature, "feature_type", "")
        or ""
    )
    for rules_for_feature_type in _iter_color_rule_sets(color_map, feature_type):
        # First, check hash pseudo-qualifier (highest priority for individual targeting)
        if 'hash' in rules_for_feature_type:
            feature_hash = get_feature_hash(feature, record_id=record_id)
            for rule in rules_for_feature_type['hash']:
                pattern, color, caption = _unpack_color_rule(rule)
                if feature_hash and pattern.search(feature_hash):
                    return color, caption

        if 'record_location' in rules_for_feature_type:
            record_location = get_feature_record_location_str(feature, record_id)
            for rule in rules_for_feature_type['record_location']:
                pattern, color, caption = _unpack_color_rule(rule)
                if record_location and pattern.search(record_location):
                    return color, caption

        # Then, check regular qualifiers
        qualifiers = get_feature_qualifiers(feature)
        for qualifier_key in qualifiers:
            if qualifier_key in rules_for_feature_type:
                # Check each pattern for this qualifier
                for rule in rules_for_feature_type[qualifier_key]:
                    # Support both old (pattern, color) and new (pattern, color, caption) formats
                    pattern, color, caption = _unpack_color_rule(rule)
                    for value in get_qualifier_values(qualifiers, qualifier_key):
                        if pattern.search(value):
                            return color, caption

        # Then, check location pseudo-qualifier (for features without unique qualifiers)
        if 'location' in rules_for_feature_type:
            # Create location string: "start..end" (GenBank format)
            location_str = get_feature_location_str(feature)
            for rule in rules_for_feature_type['location']:
                pattern, color, caption = _unpack_color_rule(rule)
                if location_str and pattern.search(location_str):
                    return color, caption
    return None


def feature_matches_specific_color_rule(
    feature: SeqFeature,
    color_map: dict,
    record_id: Optional[str] = None,
) -> bool:
    return find_specific_color_rule(feature, color_map, record_id=record_id) is not None


def get_color_with_info(
    feature: SeqFeature,
    color_map: dict,
    default_color_map: dict,
    record_id: Optional[str] = None,
) -> Tuple[str, Optional[str]]:
    """
    Determines the color for a given feature based on its type and qualifiers.

    Supports special pseudo-qualifiers for targeting features:
    - 'hash': Match by feature hash (most reliable, based on type+position+strand)
    - 'record_location': Match by "record:start..end:strand"
    - 'location': Match by position string "start..end"

    Returns:
        Tuple of (color, caption). Caption is None if default color was used.
        If record_id is provided, hash matching includes the record id.
    """
    matched_rule = find_specific_color_rule(feature, color_map, record_id=record_id)
    if matched_rule is not None:
        return matched_rule

    # Fallback to default color if no specific rule matched
    return default_color_map.get(feature.type, "#d3d3d3"), None


def get_color(feature: SeqFeature, color_map: dict, default_color_map: dict, record_id: Optional[str] = None) -> str:
    """
    Determines the color for a given feature based on its type and qualifiers.

    This is a convenience wrapper that returns only the color.
    Use get_color_with_info() if you also need the matched caption.
    If record_id is provided, hash matching includes the record id.
    """
    color, _ = get_color_with_info(feature, color_map, default_color_map, record_id=record_id)
    return color


def precompute_used_color_rules(
    records,
    color_map: dict,
    default_color_map: dict,
    selected_features_set: set,
    feature_visibility_rules: list[dict] | None = None,
) -> tuple[set, set]:
    """
    Pre-compute which color rules will be used for a set of records.

    This is useful for generating accurate legends before rendering features.

    Args:
        records: SeqRecord or list of SeqRecords
        color_map: Preprocessed color map from preprocess_color_tables
        default_color_map: Preprocessed default color map
        selected_features_set: Set of feature types to consider

    Returns:
        (used_rules, default_used_features)
        used_rules: Set of (caption, color) tuples for rules that match features
        default_used_features: Set of feature types that fell back to default color
    """
    from Bio.SeqRecord import SeqRecord
    from .visibility import should_render_feature
    if isinstance(records, SeqRecord):
        records = [records]

    used_rules: set = set()
    default_used_features: set = set()
    for record in records:
        for feature in record.features:
            if not should_render_feature(
                feature,
                selected_features_set,
                feature_visibility_rules=feature_visibility_rules,
                record_id=record.id,
                specific_color_rules=color_map,
            ):
                continue
            color, caption = get_color_with_info(feature, color_map, default_color_map, record_id=record.id)
            if caption is None:
                default_used_features.add(feature.type)
            elif caption:
                used_rules.add((caption, color))
    return used_rules, default_used_features


__all__ = [
    "feature_matches_specific_color_rule",
    "find_specific_color_rule",
    "get_color",
    "get_color_with_info",
    "precompute_used_color_rules",
    "preprocess_color_tables",
]
