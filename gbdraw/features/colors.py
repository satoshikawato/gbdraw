#!/usr/bin/env python
# coding: utf-8

import re
from typing import Any, Iterable, Mapping, Optional, Tuple, cast

from pandas import DataFrame
from Bio.SeqFeature import SeqFeature

from gbdraw.core.color import normalize_hex_color
from gbdraw.exceptions import ValidationError

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


# The value of a color rule's color that draws the features it matches in
# their type's color (``feature_type_color``: the palette with the user's
# Default colors), whatever that color is when the figure is drawn. Read in any
# letter case, as ``none`` is; Reset fill color writes it. The Web reads this
# constant through tools/generate_web_contracts.py (``TYPE_COLOR_VALUE``).
TYPE_COLOR_VALUE = "palette"


def is_type_color_value(color: object) -> bool:
    """Whether a rule's color is ``TYPE_COLOR_VALUE`` (any letter case, spaces around it)."""

    return isinstance(color, str) and color.strip().lower() == TYPE_COLOR_VALUE


def rule_color(feature_type: str, color: str, palette: Mapping[str, str]) -> str:
    """The color a rule of ``color`` draws a feature of the type in: the type's
    color (``feature_type_color`` of ``palette``) for ``TYPE_COLOR_VALUE``, else
    ``color`` itself."""

    return feature_type_color(palette, feature_type) if is_type_color_value(color) else color


def color_identity(color: object) -> str:
    """The identity two colors compare by (a name equals its hex), as a Legend compares fills."""

    from gbdraw.io.colors import resolve_color_to_hex

    text = str(color or "").strip()
    if text.lower() == "none":
        return "none"
    try:
        return normalize_hex_color(resolve_color_to_hex(text))
    except (ValueError, ValidationError):
        return text.lower()


def other_feature_row_name(feature_type: str) -> str:
    """The name of a type's "other" Legend row, of its features drawn in the
    type's color beside a captioned rule of the type: "other proteins" for CDS,
    else "other <type>s"."""

    return "other proteins" if feature_type == "CDS" else f"other {feature_type}s"


# A color rule as a Legend reads the feature that takes it: its caption ("" for
# none) and its color.
RuleUse = tuple[str, str]
# A captioned rule a drawn feature takes, as a Legend counts it: the feature's
# type, the caption and the color.
TypeRuleUse = tuple[str, str, str]


def rule_use(feature_type: str, caption: Optional[str], color: str, palette: Mapping[str, str]) -> RuleUse | None:
    """The rule of ``caption`` and ``color`` as the Legend reads a feature of
    the type that takes it: its caption and the color it draws the feature in
    (``rule_color``). None for a rule of ``TYPE_COLOR_VALUE`` without a
    caption: the feature is drawn as no rule would draw it, so it counts as a
    feature of no rule (its type's own or "other" row, Legend contract 05A4-04)."""

    if is_type_color_value(color) and not str(caption or "").strip():
        return None
    return str(caption or ""), rule_color(feature_type, color, palette)


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
        Tuple of (color, caption). Caption is None if default color was used,
        also by a rule of ``TYPE_COLOR_VALUE`` without a caption (``rule_use``).
        If record_id is provided, hash matching includes the record id.
        ``instance`` is the feature's instance (``ids.feature_instances``),
        which a `hash` rule may name; None for a feature without one.
    """
    matched_rule = find_specific_color_rule(feature, color_map, record_id=record_id, instance=instance)
    if matched_rule is not None:
        use = rule_use(feature.type, matched_rule[1], matched_rule[0], default_color_map)
        if use is not None:
            return use[1], matched_rule[1]

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


def legend_rule_use(
    drawn: Iterable[tuple[str, Optional[RuleUse]]],
) -> tuple[set[TypeRuleUse], set[str]]:
    """Which color rules a Legend counts as used, per drawn type, from each drawn
    feature's type and the rule it takes as ``rule_use`` reads it (None: no
    rule, its type's color).

    The one reading of the drawn features for the Legend (Generate's
    ``precompute_used_color_rules``, the Web's Legend order and the row of each
    feature, ``FeatureLegendRows``): the captioned rules features take, as
    (type, caption, color), and the types of the features no rule colors. A
    rule is used for a type only when features of that type take it, so a
    type it colors none of keeps its own row (Legend contract 05A4-04 L-C3,
    OV-3465). An uncaptioned rule names no row.
    """

    used: set[TypeRuleUse] = set()
    default_types: set[str] = set()
    for feature_type, winner in drawn:
        if winner is None:
            default_types.add(feature_type)
        elif winner[0]:
            used.add((feature_type, *winner))
    return used, default_types


def precompute_used_color_rules(
    records,
    color_map: Mapping[str, Any],
    default_color_map: Mapping[str, str],
    selected_features_set: set,
    feature_visibility_rules: list[dict] | None = None,
    record_features=(),
    include_hidden: bool = False,
) -> tuple[set[TypeRuleUse], set[str]]:
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
        used_rules: Set of (feature type, caption, color) tuples for the captioned
            rules features of that type take (``legend_rule_use``)
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
    "TYPE_COLOR_VALUE",
    "TypeRuleUse",
    "feature_matches_specific_color_rule",
    "find_specific_color_rule",
    "get_color",
    "get_color_with_info",
    "legend_rule_use",
    "precompute_used_color_rules",
    "build_default_color_map",
    "build_specific_color_rules",
    "color_identity",
    "feature_type_color",
    "is_type_color_value",
    "other_feature_row_name",
    "preprocess_color_tables",
    "rule_color",
    "rule_use",
]
