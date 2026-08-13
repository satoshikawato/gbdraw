#!/usr/bin/env python
# coding: utf-8

import re
from typing import Optional, Tuple

from pandas import DataFrame
from Bio.SeqFeature import SeqFeature

from .instance_identity import (
    FEATURE_INSTANCE_HASH_QUALIFIER,
    FeatureInstanceIdentity,
    FeatureInstanceIdentityPlan,
    validate_feature_instance_hash,
)
from .semantic_selectors import (
    FEATURE_SEMANTIC_SCOPE_QUALIFIER,
    FeatureSemanticSelectorContext,
    parse_feature_semantic_selector,
)
from .ids import compute_feature_hash as compute_feature_hash
from .selector_values import feature_matches_specific_color_rule, find_specific_color_rule


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
        rule_schema = int(color_table.attrs.get("gbdraw_specific_rule_schema", 6))
        for row in color_table.itertuples(index=False):
            qualifier_key = str(row.qualifier_key)
            # The reserved instance selector is an opaque literal. Every other
            # qualifier keeps the legacy case-insensitive regex contract.
            if qualifier_key == FEATURE_INSTANCE_HASH_QUALIFIER and rule_schema >= 6:
                pattern = validate_feature_instance_hash(row.value)
            elif qualifier_key == FEATURE_SEMANTIC_SCOPE_QUALIFIER and rule_schema >= 6:
                pattern = parse_feature_semantic_selector(row.value)
            else:
                pattern = re.compile(row.value, re.IGNORECASE)

            # Build nested dictionary structure
            qualifier_rules = color_map.setdefault(row.feature_type, {})
            rule_list = qualifier_rules.setdefault(qualifier_key, [])

            # Include caption for legend tracking
            caption = getattr(row, 'caption', '') or ''
            rule_list.append((pattern, row.color, caption))

    return color_map, default_color_map










def get_color_with_info(
    feature: SeqFeature,
    color_map: dict,
    default_color_map: dict,
    record_id: Optional[str] = None,
    feature_instance_identity: FeatureInstanceIdentity | None = None,
    feature_semantic_selector_context: FeatureSemanticSelectorContext | None = None,
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
    matched_rule = find_specific_color_rule(
        feature,
        color_map,
        record_id=record_id,
        feature_instance_identity=feature_instance_identity,
        feature_semantic_selector_context=feature_semantic_selector_context,
    )
    if matched_rule is not None:
        return matched_rule

    # Fallback to default color if no specific rule matched
    return default_color_map.get(feature.type, "#d3d3d3"), None


def get_color(
    feature: SeqFeature,
    color_map: dict,
    default_color_map: dict,
    record_id: Optional[str] = None,
    feature_instance_identity: FeatureInstanceIdentity | None = None,
    feature_semantic_selector_context: FeatureSemanticSelectorContext | None = None,
) -> str:
    """
    Determines the color for a given feature based on its type and qualifiers.

    This is a convenience wrapper that returns only the color.
    Use get_color_with_info() if you also need the matched caption.
    If record_id is provided, hash matching includes the record id.
    """
    color, _ = get_color_with_info(
        feature,
        color_map,
        default_color_map,
        record_id=record_id,
        feature_instance_identity=feature_instance_identity,
        feature_semantic_selector_context=feature_semantic_selector_context,
    )
    return color


def precompute_used_color_rules(
    records,
    color_map: dict,
    default_color_map: dict,
    selected_features_set: set,
    feature_visibility_rules: list[dict] | None = None,
    feature_instance_identity_plan: FeatureInstanceIdentityPlan | None = None,
    feature_semantic_selector_context: FeatureSemanticSelectorContext | None = None,
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
            feature_instance_identity = (
                feature_instance_identity_plan.identity_for_feature(feature)
                if feature_instance_identity_plan is not None
                else None
            )
            if not should_render_feature(
                feature,
                selected_features_set,
                feature_visibility_rules=feature_visibility_rules,
                record_id=record.id,
                specific_color_rules=color_map,
                feature_instance_identity=feature_instance_identity,
                feature_semantic_selector_context=feature_semantic_selector_context,
            ):
                continue
            color, caption = get_color_with_info(
                feature,
                color_map,
                default_color_map,
                record_id=record.id,
                feature_instance_identity=feature_instance_identity,
                feature_semantic_selector_context=feature_semantic_selector_context,
            )
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
