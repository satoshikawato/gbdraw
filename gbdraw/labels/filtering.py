#!/usr/bin/env python
# coding: utf-8

import logging
import hashlib
import re
from typing import Any, Optional

import pandas as pd
from pandas import DataFrame
from Bio.SeqFeature import SeqFeature

from ..exceptions import InputFileError, ParseError, ValidationError
from ..features.colors import compute_feature_hash

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


def _normalize_qualifier_values(raw_values: Any) -> list[str]:
    if raw_values is None:
        return []
    if isinstance(raw_values, (list, tuple, set)):
        return [str(value) for value in raw_values if value is not None]
    return [str(raw_values)]


def _get_feature_type(feature: Any) -> str:
    feature_type = getattr(feature, "type", None)
    if feature_type is None:
        feature_type = getattr(feature, "feature_type", "")
    return str(feature_type or "")


def _get_feature_qualifiers(feature: Any) -> dict:
    qualifiers = getattr(feature, "qualifiers", {})
    if qualifiers is None:
        return {}
    if isinstance(qualifiers, dict):
        return qualifiers
    try:
        return dict(qualifiers)
    except Exception:
        return {}


def _get_feature_record_id(feature: Any, record_id: Optional[str]) -> Optional[str]:
    if record_id is not None and str(record_id).strip() != "":
        return str(record_id)
    feature_record_id = getattr(feature, "record_id", None)
    if feature_record_id is None or str(feature_record_id).strip() == "":
        return None
    return str(feature_record_id)


def _extract_first_location_part(feature: Any) -> Optional[tuple[int, int, Any]]:
    if isinstance(feature, SeqFeature):
        loc = feature.location
        if hasattr(loc, "parts") and loc.parts:
            part = loc.parts[0]
            return int(part.start), int(part.end), part.strand
        return int(loc.start), int(loc.end), loc.strand

    location = getattr(feature, "location", None)
    if not location:
        return None

    selected = None
    for part in location:
        kind = getattr(part, "kind", None)
        if kind is None and isinstance(part, tuple) and len(part) >= 1:
            kind = part[0]
        if kind == "block":
            selected = part
            break
    if selected is None:
        selected = location[0]

    try:
        start = getattr(selected, "start", selected[3])
        end = getattr(selected, "end", selected[4])
        strand = getattr(selected, "strand", selected[2])
    except Exception:
        return None
    try:
        return int(start), int(end), strand
    except Exception:
        return None


def _normalize_strand_for_hash(strand: Any) -> Any:
    if strand in (None, "", "none", "None", "undefined"):
        return None
    if isinstance(strand, str):
        normalized = strand.strip().lower()
        if normalized in {"positive", "plus", "+", "forward", "1"}:
            return 1
        if normalized in {"negative", "minus", "-", "reverse", "-1"}:
            return -1
        if normalized in {"undefined", "none", ""}:
            return None
        try:
            return int(normalized)
        except ValueError:
            return strand
    if isinstance(strand, (int, float)):
        try:
            return int(strand)
        except Exception:
            return strand
    return strand


def _get_feature_hash(feature: Any, record_id: Optional[str]) -> Optional[str]:
    resolved_record_id = _get_feature_record_id(feature, record_id)
    if isinstance(feature, SeqFeature):
        return compute_feature_hash(feature, record_id=resolved_record_id)

    first_part = _extract_first_location_part(feature)
    feature_type = _get_feature_type(feature)
    if not first_part or not feature_type:
        return None
    start, end, strand = first_part
    normalized_strand = _normalize_strand_for_hash(strand)
    if resolved_record_id is not None:
        key = f"{resolved_record_id}:{feature_type}:{start}:{end}:{normalized_strand}"
    else:
        key = f"{feature_type}:{start}:{end}:{normalized_strand}"
    return "f" + hashlib.md5(key.encode()).hexdigest()[:8]


def _get_feature_location_str(feature: Any) -> Optional[str]:
    if isinstance(feature, SeqFeature):
        try:
            return f"{int(feature.location.start)}..{int(feature.location.end)}"
        except Exception:
            return None

    location = getattr(feature, "location", None)
    if not location:
        return None

    starts: list[int] = []
    ends: list[int] = []
    for part in location:
        kind = getattr(part, "kind", None)
        if kind is None and isinstance(part, tuple) and len(part) >= 1:
            kind = part[0]
        if kind not in {"block", None}:
            continue
        try:
            start = int(getattr(part, "start", part[3]))
            end = int(getattr(part, "end", part[4]))
        except Exception:
            continue
        starts.append(start)
        ends.append(end)
    if not starts or not ends:
        return None
    return f"{min(starts)}..{max(ends)}"


def _normalize_strand_token(strand: Any) -> str:
    if strand in (None, "", "none", "None", "undefined"):
        return "undefined"
    if isinstance(strand, str):
        normalized = strand.strip().lower()
        if normalized in {"positive", "plus", "+", "forward", "1"}:
            return "+"
        if normalized in {"negative", "minus", "-", "reverse", "-1"}:
            return "-"
        if normalized in {"undefined", "none", ""}:
            return "undefined"
        return "undefined"
    if isinstance(strand, (int, float)):
        try:
            numeric = int(strand)
        except Exception:
            return "undefined"
        if numeric == 1:
            return "+"
        if numeric == -1:
            return "-"
        return "undefined"
    return "undefined"


def _get_feature_record_location_str(feature: Any, record_id: Optional[str]) -> Optional[str]:
    resolved_record_id = _get_feature_record_id(feature, record_id)
    if not resolved_record_id:
        return None

    position = _get_feature_position_str(feature)
    if not position:
        return None

    return f"{resolved_record_id}:{position}"


def _get_feature_position_str(feature: Any) -> Optional[str]:
    location = _get_feature_location_str(feature)
    if not location:
        return None

    strand_source = None
    first_part = _extract_first_location_part(feature)
    if first_part is not None:
        strand_source = first_part[2]
    strand = _normalize_strand_token(strand_source)
    return f"{location}:{strand}"


def _get_qualifier_values(qualifiers: dict, qualifier_key: str) -> list[str]:
    for key, values in qualifiers.items():
        if str(key).lower() != str(qualifier_key).lower():
            continue
        return _normalize_qualifier_values(values)
    return []


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


def _matches_constraint(rule_token: str, actual_value: Optional[str]) -> bool:
    if str(rule_token) == "*":
        return True
    if actual_value is None:
        return False
    return str(rule_token) == str(actual_value)


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
    if whitelist_df is not None and not whitelist_df.empty:
        whitelist_map = {}
        for row in whitelist_df.itertuples(index=False):
            feature_type = str(getattr(row, "feature_type", "") or "").strip()
            qualifier = str(getattr(row, "qualifier", "") or "").strip().lower()
            keyword = str(getattr(row, "keyword", "") or "").strip()
            if not feature_type or not qualifier:
                continue
            whitelist_map.setdefault(feature_type, {}).setdefault(qualifier, set()).add(keyword)
    else:
        whitelist_map = None

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

        hash_keywords = rules.get("hash", set())
        if hash_keywords:
            feature_hash = _get_feature_hash(feature, record_id)
            if feature_hash and any(feature_hash == str(keyword) for keyword in hash_keywords):
                is_eligible = True

        if not is_eligible:
            record_location_keywords = rules.get("record_location", set())
            if record_location_keywords:
                feature_record_location = _get_feature_record_location_str(feature, record_id)
                if feature_record_location and any(
                    feature_record_location == str(keyword)
                    for keyword in record_location_keywords
                ):
                    is_eligible = True

        if not is_eligible:
            location_keywords = rules.get("location", set())
            if location_keywords:
                feature_location = _get_feature_location_str(feature)
                if feature_location and any(
                    feature_location == str(keyword)
                    for keyword in location_keywords
                ):
                    is_eligible = True

        if not is_eligible:
            for key, values in qualifiers.items():
                normalized_values = _normalize_qualifier_values(values)
                qualifier_key = str(key).lower()
                if qualifier_key in rules and any(v in rules[qualifier_key] for v in normalized_values):
                    is_eligible = True
                    break
        if not is_eligible:
            return ""

    if not whitelist_map and any(bl in final_label.lower() for bl in blacklist):
        return ""

    if final_label and non_hash_rules:
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


