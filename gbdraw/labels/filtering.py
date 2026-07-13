#!/usr/bin/env python
# coding: utf-8

import logging
import re
from typing import Any, Optional

import pandas as pd
from pandas import DataFrame

from ..exceptions import InputFileError, ParseError, ValidationError
from ..features.selector_values import (
    _matches_constraint,
    get_feature_hash as _get_feature_hash,
    get_feature_location_str as _get_feature_location_str,
    get_feature_qualifiers as _get_feature_qualifiers,
    get_feature_record_id as _get_feature_record_id,
    get_feature_record_location_str as _get_feature_record_location_str,
    get_feature_type as _get_feature_type,
    get_qualifier_values as _get_qualifier_values,
    normalize_qualifier_values as _normalize_qualifier_values,
)

logger = logging.getLogger(__name__)

DEFAULT_LABEL_PRIORITY = ["product", "gene", "locus_tag", "protein_id", "old_locus_tag", "note"]
_LABEL_OVERRIDE_REQUIRED_COLS = [
    "record_id",
    "feature_type",
    "qualifier",
    "value",
    "label_text",
]
_LABEL_OVERRIDE_HEADER = (
    "record_id",
    "feature_type",
    "qualifier",
    "value",
    "label_text",
)


def read_qualifier_priority_file(filepath: str) -> Optional[DataFrame]:
    if not filepath:
        return None

    required_cols = ["feature_type", "priorities"]

    try:
        df = pd.read_csv(
            filepath,
            sep="\t",
            header=None,
            names=required_cols,
            dtype=str,
            on_bad_lines="error",
            engine="python",
        )
    except pd.errors.ParserError as e:
        logger.error(f"ERROR: Malformed line in qualifier priority file '{filepath}': {e}")
        raise ParseError(
            f"Malformed line in qualifier priority file '{filepath}': {e}"
        ) from e
    except FileNotFoundError as e:
        logger.error(f"ERROR: Qualifier priority file not found: {e}")
        raise InputFileError(f"Qualifier priority file not found: {e}") from e
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

    logger.info(f"Successfully loaded qualifier priority from {filepath}")
    return df


def read_filter_list_file(filepath: str) -> Optional[DataFrame]:
    if not filepath:
        return None

    required_cols = ["feature_type", "qualifier", "keyword"]

    try:
        df = pd.read_csv(
            filepath,
            sep="\t",
            header=None,
            names=required_cols,
            dtype=str,
            comment="#",
            on_bad_lines="error",
            engine="python",
        )
    except pd.errors.ParserError as e:
        logger.error(f"ERROR: Malformed line in filter list file '{filepath}': {e}")
        raise ParseError(
            f"Malformed line in filter list file '{filepath}': {e}"
        ) from e
    except FileNotFoundError as e:
        logger.error(f"ERROR: Filter list file not found: {e}")
        raise InputFileError(f"Filter list file not found: {e}") from e
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

    logger.info(f"Successfully loaded filter list from {filepath}")
    return df


def read_label_override_file(filepath: str) -> Optional[DataFrame]:
    if not filepath:
        return None

    required_cols = _LABEL_OVERRIDE_REQUIRED_COLS

    try:
        with open(filepath, "r", encoding="utf-8") as handle:
            for line_no, raw_line in enumerate(handle, start=1):
                stripped = raw_line.strip()
                if stripped == "" or stripped.startswith("#"):
                    continue
                if raw_line.rstrip("\r\n").count("\t") + 1 > len(required_cols):
                    logger.error(
                        f"ERROR: Malformed line in label override file '{filepath}' at line {line_no}: "
                        f"expected {len(required_cols)} columns."
                    )
                    raise ParseError(
                        f"Malformed line in label override file '{filepath}' at line {line_no}: "
                        f"expected {len(required_cols)} columns."
                    )
    except FileNotFoundError as e:
        logger.error(f"ERROR: Label override file not found: {e}")
        raise InputFileError(f"Label override file not found: {e}") from e

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
        logger.error(f"ERROR: Malformed line in label override file '{filepath}': {e}")
        raise ParseError(
            f"Malformed line in label override file '{filepath}': {e}"
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

    logger.info(f"Successfully loaded label override table from {filepath}")
    return df






























def _is_label_override_header_row(
    record_id: str,
    feature_type: str,
    qualifier: str,
    value_pattern: str,
    label_text: str,
) -> bool:
    normalized = (
        record_id.strip().lower(),
        feature_type.strip().lower(),
        qualifier.strip().lower(),
        value_pattern.strip().lower(),
        label_text.strip().lower(),
    )
    if normalized == _LABEL_OVERRIDE_HEADER:
        return True

    return (
        normalized[0] in {"record_id", "record"}
        and normalized[1] == "feature_type"
        and normalized[2] in {"qualifier", "qualifier_key"}
        and normalized[3] in {"value", "qualifier_value_regex"}
        and normalized[4] == "label_text"
    )


def _build_label_override_rules(label_override_df: Optional[DataFrame]) -> Optional[list[dict[str, Any]]]:
    if label_override_df is None or label_override_df.empty:
        return None

    compiled_rules: list[dict[str, Any]] = []
    for idx, row in enumerate(label_override_df.itertuples(index=False), start=1):
        record_id = str(getattr(row, "record_id", getattr(row, "record", "")) or "").strip()
        feature_type = str(getattr(row, "feature_type", "") or "").strip()
        qualifier = str(getattr(row, "qualifier", getattr(row, "qualifier_key", "")) or "").strip()
        value_pattern = str(getattr(row, "value", "") or "").strip()
        label_text = str(getattr(row, "label_text", "") or "")

        if _is_label_override_header_row(
            record_id=record_id,
            feature_type=feature_type,
            qualifier=qualifier,
            value_pattern=value_pattern,
            label_text=label_text,
        ):
            continue

        if not record_id:
            logger.error(f"ERROR: Missing record_id token in label override table at row {idx}.")
            raise ValidationError(
                f"Missing record_id token in label override table at row {idx}."
            )
        if not feature_type:
            logger.error(f"ERROR: Missing feature_type token in label override table at row {idx}.")
            raise ValidationError(
                f"Missing feature_type token in label override table at row {idx}."
            )
        if not qualifier:
            logger.error(f"ERROR: Missing qualifier token in label override table at row {idx}.")
            raise ValidationError(
                f"Missing qualifier token in label override table at row {idx}."
            )
        if not value_pattern:
            logger.error(f"ERROR: Missing value regex token in label override table at row {idx}.")
            raise ValidationError(
                f"Missing value regex token in label override table at row {idx}."
            )

        try:
            pattern = re.compile(value_pattern, re.IGNORECASE)
        except re.error as e:
            logger.error(
                f"ERROR: Invalid regex in label override table at row {idx}: '{value_pattern}' ({e})"
            )
            raise ParseError(
                f"Invalid regex in label override table at row {idx}: '{value_pattern}' ({e})"
            ) from e

        compiled_rules.append(
            {
                "record_id": record_id,
                "feature_type": feature_type,
                "qualifier": qualifier,
                "qualifier_normalized": qualifier.lower(),
                "pattern": pattern,
                "label_text": label_text,
            }
        )

    return compiled_rules or None


def _build_whitelist_map(whitelist_df: Optional[DataFrame]) -> Optional[dict[str, dict[str, list[re.Pattern[str]]]]]:
    if whitelist_df is None or whitelist_df.empty:
        return None

    whitelist_map: dict[str, dict[str, list[re.Pattern[str]]]] = {}
    for idx, row in enumerate(whitelist_df.itertuples(index=False), start=1):
        feature_type = str(getattr(row, "feature_type", "") or "").strip()
        qualifier = str(getattr(row, "qualifier", "") or "").strip().lower()
        value_pattern = str(getattr(row, "keyword", "") or "").strip()
        if not feature_type or not qualifier:
            continue

        try:
            pattern = re.compile(value_pattern, re.IGNORECASE)
        except re.error as e:
            logger.error(
                f"ERROR: Invalid regex in label whitelist table at row {idx}: '{value_pattern}' ({e})"
            )
            raise ParseError(
                f"Invalid regex in label whitelist table at row {idx}: '{value_pattern}' ({e})"
            ) from e

        whitelist_map.setdefault(feature_type, {}).setdefault(qualifier, []).append(pattern)

    return whitelist_map or None


def _matches_any_pattern(candidate: Optional[str], patterns: list[re.Pattern[str]]) -> bool:
    if candidate is None:
        return False
    return any(pattern.search(candidate) for pattern in patterns)




def _resolve_label_override(
    feature: Any,
    feature_type: str,
    qualifiers: dict,
    base_label: str,
    rules: list[dict[str, Any]],
    record_id: Optional[str],
) -> Optional[str]:
    resolved_record_id = _get_feature_record_id(feature, record_id)
    feature_location = None
    feature_hash = None
    record_location = None

    for rule in rules:
        if not _matches_constraint(rule["record_id"], resolved_record_id):
            continue
        if not _matches_constraint(rule["feature_type"], feature_type):
            continue

        qualifier_key = rule["qualifier_normalized"]
        pattern = rule["pattern"]
        replacement = rule["label_text"]

        if qualifier_key == "record_location":
            if record_location is None:
                record_location = _get_feature_record_location_str(feature, record_id)
            if record_location and pattern.search(record_location):
                return replacement
            continue

        if qualifier_key == "hash":
            if feature_hash is None:
                feature_hash = _get_feature_hash(feature, record_id)
            if feature_hash and pattern.search(feature_hash):
                return replacement
            continue

        if qualifier_key == "location":
            if feature_location is None:
                feature_location = _get_feature_location_str(feature)
            if feature_location and pattern.search(feature_location):
                return replacement
            continue

        if qualifier_key == "label":
            if pattern.search(base_label):
                return replacement
            continue

        values = _get_qualifier_values(qualifiers, rule["qualifier"])
        if any(pattern.search(value) for value in values):
            return replacement
    return None


def preprocess_label_filtering(label_filtering: dict):
    # If already processed, return as-is unless label override rules still need compilation.
    override_df = label_filtering.get("label_override_df")
    if isinstance(override_df, str):
        override_df = None
    already_processed = "whitelist_map" in label_filtering and "priority_map" in label_filtering
    if already_processed and ("label_override_rules" in label_filtering or override_df is None):
        return label_filtering

    whitelist_df = label_filtering.get("whitelist_df")
    # Handle case where whitelist_df is a string (should be DataFrame or None)
    # This can happen if modify_config_dict overwrote a DataFrame with a string
    if isinstance(whitelist_df, str):
        whitelist_df = None
    whitelist_map = _build_whitelist_map(whitelist_df)

    priority_df = label_filtering.get("qualifier_priority_df")
    # Handle case where priority_df is a string (should be DataFrame or None)
    if isinstance(priority_df, str):
        priority_df = None
    if priority_df is not None and not priority_df.empty:
        priority_map = {}
        for row in priority_df.itertuples(index=False):
            priority_map[row.feature_type] = [p.strip() for p in row.priorities.split(",")]
    else:
        priority_map = {}

    label_override_rules = _build_label_override_rules(override_df)

    label_filtering["whitelist_map"] = whitelist_map
    label_filtering["priority_map"] = priority_map
    label_filtering["label_override_rules"] = label_override_rules
    return label_filtering


def get_label_text(feature: Any, label_filtering: dict, record_id: Optional[str] = None) -> str:
    feature_type = _get_feature_type(feature)
    qualifiers = _get_feature_qualifiers(feature)
    whitelist_map = label_filtering.get("whitelist_map")
    blacklist = label_filtering.get("blacklist_keywords", [])
    priority_map = label_filtering.get("priority_map", {})
    label_override_rules = label_filtering.get("label_override_rules")
    hash_rules = None
    non_hash_rules = None
    if label_override_rules:
        hash_rules = [rule for rule in label_override_rules if rule.get("qualifier_normalized") == "hash"]
        non_hash_rules = [rule for rule in label_override_rules if rule.get("qualifier_normalized") != "hash"]

    priority_list = priority_map.get(feature_type, DEFAULT_LABEL_PRIORITY)
    final_label = ""
    for key in priority_list:
        values = _get_qualifier_values(qualifiers, key)
        if values:
            final_label = values[0]
            break

    # Hash-based per-feature overrides are highest priority and can bypass whitelist/blacklist filters.
    if hash_rules:
        hash_override_text = _resolve_label_override(
            feature=feature,
            feature_type=feature_type,
            qualifiers=qualifiers,
            base_label=final_label,
            rules=hash_rules,
            record_id=record_id,
        )
        if hash_override_text is not None:
            return hash_override_text

    if whitelist_map:
        is_eligible = False
        rules = whitelist_map.get(feature_type, {})
        if not isinstance(rules, dict):
            rules = {}

        hash_patterns = rules.get("hash", [])
        if hash_patterns:
            feature_hash = _get_feature_hash(feature, record_id)
            if _matches_any_pattern(feature_hash, hash_patterns):
                is_eligible = True

        if not is_eligible:
            record_location_patterns = rules.get("record_location", [])
            if record_location_patterns:
                feature_record_location = _get_feature_record_location_str(feature, record_id)
                if _matches_any_pattern(feature_record_location, record_location_patterns):
                    is_eligible = True

        if not is_eligible:
            location_patterns = rules.get("location", [])
            if location_patterns:
                feature_location = _get_feature_location_str(feature)
                if _matches_any_pattern(feature_location, location_patterns):
                    is_eligible = True

        if not is_eligible:
            for key, values in qualifiers.items():
                normalized_values = _normalize_qualifier_values(values)
                qualifier_key = str(key).lower()
                patterns = rules.get(qualifier_key, [])
                if patterns and any(_matches_any_pattern(value, patterns) for value in normalized_values):
                    is_eligible = True
                    break
        if not is_eligible:
            return ""

    if not whitelist_map and any(str(bl).lower() in final_label.lower() for bl in blacklist):
        return ""

    if non_hash_rules:
        override_text = _resolve_label_override(
            feature=feature,
            feature_type=feature_type,
            qualifiers=qualifiers,
            base_label=final_label,
            rules=non_hash_rules,
            record_id=record_id,
        )
        if override_text is not None:
            return override_text

    return final_label


__all__ = [
    "get_label_text",
    "preprocess_label_filtering",
    "read_filter_list_file",
    "read_label_override_file",
    "read_qualifier_priority_file",
]


