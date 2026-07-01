#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations

import logging
import re
import warnings
from typing import Any, Optional, Sequence

import pandas as pd
from pandas import DataFrame

from ..exceptions import InputFileError, ParseError, ValidationError
from .selector_values import (
    get_feature_hash,
    get_feature_location_str,
    get_feature_qualifiers,
    get_feature_record_id,
    get_feature_record_location_str,
    get_feature_type,
    get_qualifier_values,
)

logger = logging.getLogger(__name__)

_FEATURE_VISIBILITY_REQUIRED_COLS = [
    "record_id",
    "feature_type",
    "qualifier",
    "value",
    "action",
]
_FEATURE_VISIBILITY_HEADER = (
    "record_id",
    "feature_type",
    "qualifier",
    "value",
    "action",
)

_SHOW_ACTION_TOKENS = {"show", "on"}
_HIDE_ACTION_TOKENS = {"hide", "off", "false", "0"}
_EXCLUDE_MATCHING_ACTION_TOKENS = {"exclude_matching"}


def read_feature_visibility_file(filepath: str) -> Optional[DataFrame]:
    if not filepath:
        return None

    required_cols = _FEATURE_VISIBILITY_REQUIRED_COLS

    try:
        with open(filepath, "r", encoding="utf-8") as handle:
            for line_no, raw_line in enumerate(handle, start=1):
                stripped = raw_line.strip()
                if stripped == "" or stripped.startswith("#"):
                    continue
                if raw_line.rstrip("\r\n").count("\t") + 1 > len(required_cols):
                    logger.error(
                        f"ERROR: Malformed line in feature visibility file '{filepath}' at line {line_no}: "
                        f"expected {len(required_cols)} columns."
                    )
                    raise ParseError(
                        f"Malformed line in feature visibility file '{filepath}' at line {line_no}: "
                        f"expected {len(required_cols)} columns."
                    )
    except FileNotFoundError as e:
        logger.error(f"ERROR: Feature visibility file not found: {e}")
        raise InputFileError(f"Feature visibility file not found: {e}") from e

    try:
        df = pd.read_csv(
            filepath,
            sep="\t",
            header=None,
            names=required_cols,
            dtype=str,
            comment="#",
            keep_default_na=False,
            na_filter=False,
            on_bad_lines="error",
            engine="python",
        )
    except pd.errors.ParserError as e:
        logger.error(f"ERROR: Malformed line in feature visibility file '{filepath}': {e}")
        raise ParseError(
            f"Malformed line in feature visibility file '{filepath}': {e}"
        ) from e
    except Exception as e:
        logger.error(f"ERROR: Failed to read '{filepath}': {e}")
        raise ParseError(f"Failed to read '{filepath}': {e}") from e

    null_rows = df[df.isnull().any(axis=1)]
    if not null_rows.empty:
        for idx, row in null_rows.iterrows():
            missing = [c for c in required_cols if pd.isna(row[c])]
            logger.error(
                f"ERROR: Missing values in '{filepath}' at line {idx+1}. "
                f"Missing columns: {missing}. Row data: {row.to_dict()}"
            )
        raise ValidationError(
            f"Missing values in '{filepath}'. See log for details."
        )

    logger.info(f"Successfully loaded feature visibility table from {filepath}")
    return df


def _is_feature_visibility_header_row(
    record_id: str,
    feature_type: str,
    qualifier: str,
    value_pattern: str,
    action: str,
) -> bool:
    normalized = (
        record_id.strip().lower(),
        feature_type.strip().lower(),
        qualifier.strip().lower(),
        value_pattern.strip().lower(),
        action.strip().lower(),
    )
    if normalized == _FEATURE_VISIBILITY_HEADER:
        return True

    return (
        normalized[0] in {"record_id", "record"}
        and normalized[1] == "feature_type"
        and normalized[2] in {"qualifier", "qualifier_key"}
        and normalized[3] in {"value", "qualifier_value_regex"}
        and normalized[4] == "action"
    )


def _normalize_visibility_action(action: str, row_idx: int) -> str:
    normalized = str(action or "").strip().lower()
    if normalized in _SHOW_ACTION_TOKENS:
        return "show"
    if normalized in _HIDE_ACTION_TOKENS:
        return "off"
    if normalized in _EXCLUDE_MATCHING_ACTION_TOKENS:
        return "exclude_matching"
    if normalized == "suppress":
        message = (
            "The feature visibility action 'suppress' is deprecated; "
            "use 'exclude_matching' instead."
        )
        warnings.warn(message, DeprecationWarning, stacklevel=2)
        logger.warning("WARNING: %s", message)
        return "exclude_matching"
    logger.error(
        "ERROR: Invalid action in feature visibility table at row %s: '%s'", row_idx, action
    )
    raise ValidationError(
        f"Invalid action in feature visibility table at row {row_idx}: '{action}'. "
        "Use show, off, or exclude_matching. Accepted aliases are on, hide, "
        "false, and 0; suppress is deprecated and maps to exclude_matching."
    )


def compile_feature_visibility_rules(
    feature_visibility_df: Optional[DataFrame],
) -> Optional[list[dict[str, Any]]]:
    if feature_visibility_df is None or feature_visibility_df.empty:
        return None

    compiled_rules: list[dict[str, Any]] = []
    for idx, row in enumerate(feature_visibility_df.itertuples(index=False), start=1):
        record_id = str(getattr(row, "record_id", getattr(row, "record", "")) or "").strip()
        feature_type = str(getattr(row, "feature_type", "") or "").strip()
        qualifier = str(getattr(row, "qualifier", getattr(row, "qualifier_key", "")) or "").strip()
        value_pattern = str(getattr(row, "value", "") or "").strip()
        action = str(getattr(row, "action", "") or "").strip()

        if _is_feature_visibility_header_row(
            record_id=record_id,
            feature_type=feature_type,
            qualifier=qualifier,
            value_pattern=value_pattern,
            action=action,
        ):
            continue

        if not record_id:
            logger.error(
                f"ERROR: Missing record_id token in feature visibility table at row {idx}."
            )
            raise ValidationError(
                f"Missing record_id token in feature visibility table at row {idx}."
            )
        if not feature_type:
            logger.error(
                f"ERROR: Missing feature_type token in feature visibility table at row {idx}."
            )
            raise ValidationError(
                f"Missing feature_type token in feature visibility table at row {idx}."
            )
        if not qualifier:
            logger.error(
                f"ERROR: Missing qualifier token in feature visibility table at row {idx}."
            )
            raise ValidationError(
                f"Missing qualifier token in feature visibility table at row {idx}."
            )
        if not value_pattern:
            logger.error(
                f"ERROR: Missing value regex token in feature visibility table at row {idx}."
            )
            raise ValidationError(
                f"Missing value regex token in feature visibility table at row {idx}."
            )

        try:
            pattern = re.compile(value_pattern, re.IGNORECASE)
        except re.error as e:
            logger.error(
                f"ERROR: Invalid regex in feature visibility table at row {idx}: '{value_pattern}' ({e})"
            )
            raise ParseError(
                f"Invalid regex in feature visibility table at row {idx}: '{value_pattern}' ({e})"
            ) from e

        normalized_action = _normalize_visibility_action(action, idx)

        compiled_rules.append(
            {
                "record_id": record_id,
                "feature_type": feature_type,
                "qualifier": qualifier,
                "qualifier_normalized": qualifier.lower(),
                "pattern": pattern,
                "action": normalized_action,
            }
        )

    return compiled_rules or None


def _matches_constraint(rule_token: str, actual_value: Optional[str]) -> bool:
    if str(rule_token) == "*":
        return True
    if actual_value is None:
        return False
    return str(rule_token) == str(actual_value)


def _rule_matches_feature(feature: Any, rule: dict[str, Any], record_id: Optional[str]) -> bool:
    qualifier_key = str(rule["qualifier_normalized"])
    pattern = rule["pattern"]

    if qualifier_key == "record_location":
        feature_record_location = get_feature_record_location_str(feature, record_id)
        return bool(feature_record_location and pattern.search(feature_record_location))

    if qualifier_key == "hash":
        feature_hash = get_feature_hash(feature, record_id)
        return bool(feature_hash and pattern.search(feature_hash))

    if qualifier_key == "location":
        feature_location = get_feature_location_str(feature)
        return bool(feature_location and pattern.search(feature_location))

    qualifiers = get_feature_qualifiers(feature)
    values = get_qualifier_values(qualifiers, str(rule["qualifier"]))
    return any(pattern.search(value) for value in values)


def _first_matching_visibility_rule(
    feature: Any,
    feature_visibility_rules: Optional[list[dict[str, Any]]],
    record_id: Optional[str] = None,
) -> Optional[dict[str, Any]]:
    if not feature_visibility_rules:
        return None

    feature_type = get_feature_type(feature)
    resolved_record_id = get_feature_record_id(feature, record_id)

    for rule in feature_visibility_rules:
        if not _matches_constraint(str(rule["record_id"]), resolved_record_id):
            continue
        if not _matches_constraint(str(rule["feature_type"]), feature_type):
            continue
        if not _rule_matches_feature(feature, rule, record_id=resolved_record_id):
            continue
        return rule
    return None


def resolve_candidate_feature_types(
    selected_features_set: Sequence[str] | set[str] | None,
    *,
    color_table: Optional[DataFrame] = None,
    feature_visibility_table: Optional[DataFrame] = None,
) -> tuple[set[str], bool]:
    """Return GFF3 candidate feature types and whether all types are required."""

    candidate_types = {
        str(feature_type).strip()
        for feature_type in (selected_features_set or [])
        if str(feature_type).strip()
    }

    for table in (color_table, feature_visibility_table):
        if table is None or table.empty or "feature_type" not in table.columns:
            continue
        for raw_feature_type in table["feature_type"]:
            feature_type = str(raw_feature_type or "").strip()
            if not feature_type or feature_type.lower() == "feature_type":
                continue
            if feature_type == "*":
                return candidate_types, True
            candidate_types.add(feature_type)

    return candidate_types, False


def should_render_feature(
    feature: Any,
    selected_features_set: Sequence[str] | set[str] | None,
    feature_visibility_rules: Optional[list[dict[str, Any]]] = None,
    record_id: Optional[str] = None,
    specific_color_rules: Optional[dict] = None,
) -> bool:
    feature_type = get_feature_type(feature)
    selected_set = {str(feature_name) for feature_name in (selected_features_set or [])}

    # Keep legacy behavior when no feature filter is provided.
    if not selected_set:
        base_visible = True
    else:
        base_visible = feature_type in selected_set
        if not base_visible and specific_color_rules:
            from .colors import feature_matches_specific_color_rule

            base_visible = feature_matches_specific_color_rule(
                feature,
                specific_color_rules,
                record_id=record_id,
            )

    rule = _first_matching_visibility_rule(
        feature,
        feature_visibility_rules,
        record_id=record_id,
    )
    if rule is not None:
        action = str(rule["action"])
        if action == "show":
            return True
        if action == "off":
            return False

    return base_visible


def should_include_feature_in_analysis(
    feature: Any,
    feature_visibility_rules: Optional[list[dict[str, Any]]] = None,
    record_id: Optional[str] = None,
) -> bool:
    rule = _first_matching_visibility_rule(
        feature,
        feature_visibility_rules,
        record_id=record_id,
    )
    return not (rule is not None and str(rule["action"]) in {"off", "exclude_matching"})


__all__ = [
    "compile_feature_visibility_rules",
    "read_feature_visibility_file",
    "resolve_candidate_feature_types",
    "should_include_feature_in_analysis",
    "should_render_feature",
]
