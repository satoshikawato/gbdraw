#!/usr/bin/env python
# coding: utf-8

import re
from typing import Any, Iterable, Mapping, Optional, Tuple, cast

from pandas import DataFrame
from Bio.SeqFeature import SeqFeature

from .ids import compute_feature_hash as compute_feature_hash
from .selector_values import ColorRuleList, feature_matches_specific_color_rule, find_specific_color_rule


def build_specific_color_rules(color_table: Optional[DataFrame]) -> dict[str, dict[str, ColorRuleList]]:
    """Compile the rows of a color table into feature_type -> qualifier -> rules.

    Each rule is ``(pattern, color, caption)``; patterns match case-insensitively.
    """

    rows: dict[str, dict[str, list[tuple[re.Pattern[str], str, str]]]] = {}
    if isinstance(color_table, DataFrame) and not color_table.empty:
        for row in color_table.itertuples(index=False):
            pattern = re.compile(row.value, re.IGNORECASE)
            qualifier_rules = rows.setdefault(row.feature_type, {})
            rule_list = qualifier_rules.setdefault(row.qualifier_key, [])
            # Include caption for legend tracking
            caption = getattr(row, 'caption', '') or ''
            rule_list.append((pattern, row.color, caption))
    return {
        feature_type: {key: ColorRuleList(rule_list) for key, rule_list in qualifier_rules.items()}
        for feature_type, qualifier_rules in rows.items()
    }


def build_default_color_map(default_colors: DataFrame) -> dict[str, str]:
    """Map each feature type of the default colors to its color."""

    return cast("dict[str, str]", default_colors.set_index("feature_type")["color"].to_dict())


def preprocess_color_tables(color_table: DataFrame, default_colors: DataFrame) -> tuple[dict, dict]:
    """Return ``(color_map, default_color_map)`` for feature coloring."""

    return build_specific_color_rules(color_table), build_default_color_map(default_colors)


def feature_type_color(default_color_map: Mapping[str, str], feature_type: str) -> str:
    """The color of a feature type: its palette entry, else the palette's ``default``.

    The one owner of this color: drawn features without a matching rule and
    the Legend row of their type both read it (OV-3114).
    """

    if feature_type in default_color_map:
        return default_color_map[feature_type]
    return default_color_map.get("default", "#d3d3d3")


def get_color_with_info(
    feature: SeqFeature,
    color_map: Mapping[str, Any],
    default_color_map: Mapping[str, str],
    record_id: Optional[str] = None,
    *,
    instance: Optional[str],
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
        ``instance`` is the feature's instance (``ids.feature_instances``),
        which a `hash` rule may name; None for a feature without one.
    """
    matched_rule = find_specific_color_rule(feature, color_map, record_id=record_id, instance=instance)
    if matched_rule is not None:
        return matched_rule

    # Fallback to default color if no specific rule matched
    return feature_type_color(default_color_map, feature.type), None


def get_color(
    feature: SeqFeature,
    color_map: Mapping[str, Any],
    default_color_map: Mapping[str, str],
    record_id: Optional[str] = None,
    *,
    instance: Optional[str],
) -> str:
    """
    Determines the color for a given feature based on its type and qualifiers.

    This is a convenience wrapper that returns only the color.
    Use get_color_with_info() if you also need the matched caption.
    If record_id is provided, hash matching includes the record id.
    """
    color, _ = get_color_with_info(feature, color_map, default_color_map, record_id=record_id, instance=instance)
    return color


# A color rule as a Legend reads the feature that takes it: its caption ("" for
# none) and its color.
RuleUse = tuple[str, str]


def legend_rule_use(drawn: Iterable[tuple[str, Optional[RuleUse]]]) -> tuple[set[RuleUse], set[str]]:
    """Which color rules a Legend counts as used, from each drawn feature's type
    and the rule it takes (None: no rule, its type's color).

    The one reading of the drawn features for the Legend (Generate's
    ``precompute_used_color_rules``, the Web's Legend order and the row of each
    feature, ``FeatureLegendRows``): the captioned rules features take, as
    (caption, color), and the types of the features no rule colors. An
    uncaptioned rule names no row.
    """

    used: set[RuleUse] = set()
    default_types: set[str] = set()
    for feature_type, winner in drawn:
        if winner is None:
            default_types.add(feature_type)
        elif winner[0]:
            used.add(winner)
    return used, default_types


def precompute_used_color_rules(
    records,
    color_map: Mapping[str, Any],
    default_color_map: Mapping[str, str],
    selected_features_set: set,
    feature_visibility_rules: list[dict] | None = None,
    record_features=(),
    include_hidden: bool = False,
) -> tuple[set, set]:
    """
    Pre-compute which color rules will be used for a set of records.

    This is useful for generating accurate legends before rendering features.

    Args:
        records: SeqRecord or list of SeqRecords
        color_map: Preprocessed color map from preprocess_color_tables
        default_color_map: Preprocessed default color map
        selected_features_set: Set of feature types to consider
        include_hidden: Also match the features the draft hides, so the result names
            every rule row the records can produce (Legend row facts).

    Returns:
        (used_rules, default_used_features)
        used_rules: Set of (caption, color) tuples for rules that match features
        default_used_features: Set of feature types that fell back to default color
    """
    from Bio.SeqRecord import SeqRecord
    from .ids import feature_instances
    from .overrides import feature_override_lookup
    from .visibility import should_render_feature
    if isinstance(records, SeqRecord):
        records = [records]

    drawn: list[tuple[str, RuleUse | None]] = []
    for index, record in enumerate(records):
        override_of = feature_override_lookup(
            record, record_features[index].overrides if record_features else None
        )
        instance_of = feature_instances(record)
        for feature in record.features:
            instance = instance_of(feature)
            if not include_hidden and not should_render_feature(
                feature,
                selected_features_set,
                feature_visibility_rules=feature_visibility_rules,
                record_id=record.id,
                specific_color_rules=color_map,
                feature_override=override_of(feature),
                instance=instance,
            ):
                continue
            color, caption = get_color_with_info(
                feature, color_map, default_color_map, record_id=record.id, instance=instance
            )
            drawn.append((feature.type, None if caption is None else (caption, color)))
    return legend_rule_use(drawn)


__all__ = [
    "RuleUse",
    "feature_matches_specific_color_rule",
    "find_specific_color_rule",
    "get_color",
    "get_color_with_info",
    "legend_rule_use",
    "precompute_used_color_rules",
    "build_default_color_map",
    "build_specific_color_rules",
    "feature_type_color",
    "preprocess_color_tables",
]
