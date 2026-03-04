#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations

import hashlib
import logging
import re
from typing import Any, Optional, Sequence

import pandas as pd
from Bio.SeqFeature import SeqFeature
from pandas import DataFrame

from ..exceptions import InputFileError, ParseError, ValidationError
from .colors import compute_feature_hash

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

_SHOW_ACTION_TOKENS = {"show", "on", "display", "include", "true", "1"}
_HIDE_ACTION_TOKENS = {"hide", "off", "suppress", "exclude", "false", "0"}


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


def _iter_feature_coordinates(feature: Any) -> list[Any]:
    coordinates = getattr(feature, "coordinates", None)
    if coordinates is None:
        return []
    try:
        return list(coordinates)
    except Exception:
        return []


def _extract_first_coordinate_part(feature: Any) -> Optional[tuple[int, int, Any]]:
    for part in _iter_feature_coordinates(feature):
        try:
            start = int(getattr(part, "start"))
            end = int(getattr(part, "end"))
            strand = getattr(part, "strand", None)
        except Exception:
            continue
        return start, end, strand
    return None


def _extract_coordinate_bounds(feature: Any) -> Optional[tuple[int, int]]:
    starts: list[int] = []
    ends: list[int] = []
    for part in _iter_feature_coordinates(feature):
        try:
            starts.append(int(getattr(part, "start")))
            ends.append(int(getattr(part, "end")))
        except Exception:
            continue
    if not starts or not ends:
        return None
    return min(starts), max(ends)


def _extract_first_location_part(feature: Any) -> Optional[tuple[int, int, Any]]:
    if isinstance(feature, SeqFeature):
        loc = feature.location
        if hasattr(loc, "parts") and loc.parts:
            part = loc.parts[0]
            return int(part.start), int(part.end), part.strand
        return int(loc.start), int(loc.end), loc.strand

    coordinate_part = _extract_first_coordinate_part(feature)
    if coordinate_part is not None:
        return coordinate_part

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

    coordinate_bounds = _extract_coordinate_bounds(feature)
    if coordinate_bounds is not None:
        return f"{coordinate_bounds[0]}..{coordinate_bounds[1]}"

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


def _get_feature_record_location_str(feature: Any, record_id: Optional[str]) -> Optional[str]:
    resolved_record_id = _get_feature_record_id(feature, record_id)
    if not resolved_record_id:
        return None

    position = _get_feature_position_str(feature)
    if not position:
        return None

    return f"{resolved_record_id}:{position}"


def _get_qualifier_values(qualifiers: dict, qualifier_key: str) -> list[str]:
    for key, values in qualifiers.items():
        if str(key).lower() != str(qualifier_key).lower():
            continue
        return _normalize_qualifier_values(values)
    return []


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
        return "hide"
    logger.error(
        "ERROR: Invalid action in feature visibility table at row %s: '%s'", row_idx, action
    )
    raise ValidationError(
        f"Invalid action in feature visibility table at row {row_idx}: '{action}'. "
        f"Use show/on/display/include/true/1 or hide/off/suppress/exclude/false/0."
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
        feature_record_location = _get_feature_record_location_str(feature, record_id)
        return bool(feature_record_location and pattern.search(feature_record_location))

    if qualifier_key == "hash":
        feature_hash = _get_feature_hash(feature, record_id)
        return bool(feature_hash and pattern.search(feature_hash))

    if qualifier_key == "location":
        feature_location = _get_feature_location_str(feature)
        return bool(feature_location and pattern.search(feature_location))

    qualifiers = _get_feature_qualifiers(feature)
    values = _get_qualifier_values(qualifiers, str(rule["qualifier"]))
    return any(pattern.search(value) for value in values)


def should_render_feature(
    feature: Any,
    selected_features_set: Sequence[str] | set[str] | None,
    feature_visibility_rules: Optional[list[dict[str, Any]]] = None,
    record_id: Optional[str] = None,
) -> bool:
    feature_type = _get_feature_type(feature)
    selected_set = {str(feature_name) for feature_name in (selected_features_set or [])}

    # Keep legacy behavior when no feature filter is provided.
    if not selected_set:
        base_visible = True
    else:
        base_visible = feature_type in selected_set

    if not feature_visibility_rules:
        return base_visible

    resolved_record_id = _get_feature_record_id(feature, record_id)

    for rule in feature_visibility_rules:
        if not _matches_constraint(str(rule["record_id"]), resolved_record_id):
            continue
        if not _matches_constraint(str(rule["feature_type"]), feature_type):
            continue
        if not _rule_matches_feature(feature, rule, record_id=resolved_record_id):
            continue
        return str(rule["action"]) == "show"

    return base_visible


__all__ = [
    "compile_feature_visibility_rules",
    "read_feature_visibility_file",
    "should_render_feature",
]
