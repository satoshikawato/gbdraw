"""CLI-independent values shared after request input normalization."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

from pandas import DataFrame  # type: ignore[reportMissingImports]

from gbdraw.features.colors import preprocess_color_tables
from gbdraw.features.visibility import compile_feature_visibility_rules


@dataclass(frozen=True)
class ResolvedFeatureInputs:
    """Feature tables and compiled rules resolved once for one diagram request."""

    color_table: DataFrame | None
    default_colors: DataFrame
    feature_visibility_table: DataFrame | None
    feature_visibility_rules: list[dict[str, Any]]
    specific_color_rules: Mapping[str, Any]
    default_color_map: Mapping[str, str]


def resolve_feature_inputs(
    *,
    color_table: DataFrame | None,
    default_colors: DataFrame,
    feature_visibility_table: DataFrame | None,
) -> ResolvedFeatureInputs:
    """Compile already-loaded feature inputs into one reusable value."""

    specific_color_rules, default_color_map = preprocess_color_tables(
        color_table,
        default_colors,
    )
    return ResolvedFeatureInputs(
        color_table=color_table,
        default_colors=default_colors,
        feature_visibility_table=feature_visibility_table,
        feature_visibility_rules=compile_feature_visibility_rules(
            feature_visibility_table
        ),
        specific_color_rules=specific_color_rules,
        default_color_map=default_color_map,
    )


__all__ = ["ResolvedFeatureInputs", "resolve_feature_inputs"]
