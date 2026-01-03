#!/usr/bin/env python
# coding: utf-8

import logging
import sys
from typing import Optional

import pandas as pd
from pandas import DataFrame
from Bio.SeqFeature import SeqFeature

logger = logging.getLogger(__name__)


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
        sys.exit(1)
    except FileNotFoundError as e:
        logger.error(f"ERROR: Qualifier priority file not found: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"ERROR: Failed to read '{filepath}': {e}")
        sys.exit(1)

    null_rows = df[df.isnull().any(axis=1)]
    if not null_rows.empty:
        for idx, row in null_rows.iterrows():
            missing = [c for c in required_cols if pd.isna(row[c])]
            logger.error(
                f"ERROR: Missing values in '{filepath}' at line {idx+1}. "
                f"Missing columns: {missing}. Row data: {row.to_dict()}"
            )
        sys.exit(1)

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
        sys.exit(1)
    except FileNotFoundError as e:
        logger.error(f"ERROR: Filter list file not found: {e}")
        sys.exit(1)
    except Exception as e:
        logger.error(f"ERROR: Failed to read '{filepath}': {e}")
        sys.exit(1)

    null_rows = df[df.isnull().any(axis=1)]
    if not null_rows.empty:
        for idx, row in null_rows.iterrows():
            missing = [c for c in required_cols if pd.isna(row[c])]
            logger.error(
                f"ERROR: Missing values in '{filepath}' at line {idx+1}. "
                f"Missing columns: {missing}. Row data: {row.to_dict()}"
            )
        sys.exit(1)

    logger.info(f"Successfully loaded filter list from {filepath}")
    return df


def preprocess_label_filtering(label_filtering: dict):
    # #region agent log
    import json
    log_path = "/mnt/c/Users/kawato/Documents/GitHub/gbdraw/.cursor/debug.log"
    try:
        with open(log_path, "a") as f:
            f.write(json.dumps({"id": "log_preprocess_entry", "timestamp": __import__("time").time(), "location": "filtering.py:96", "message": "preprocess_label_filtering entry", "data": {"label_filtering_keys": list(label_filtering.keys()), "whitelist_df_type": str(type(label_filtering.get("whitelist_df"))), "whitelist_df_value": str(label_filtering.get("whitelist_df"))[:100] if isinstance(label_filtering.get("whitelist_df"), str) else "not_string", "qualifier_priority_df_type": str(type(label_filtering.get("qualifier_priority_df"))), "qualifier_priority_df_value": str(label_filtering.get("qualifier_priority_df"))[:100] if isinstance(label_filtering.get("qualifier_priority_df"), str) else "not_string"}, "sessionId": "debug-session", "runId": "run1", "hypothesisId": "A"}) + "\n")
    except: pass
    # #endregion
    # If already processed, return as-is
    if "whitelist_map" in label_filtering and "priority_map" in label_filtering:
        return label_filtering
    
    whitelist_df = label_filtering.get("whitelist_df")
    # #region agent log
    try:
        with open(log_path, "a") as f:
            f.write(json.dumps({"id": "log_whitelist_check", "timestamp": __import__("time").time(), "location": "filtering.py:98", "message": "whitelist_df check", "data": {"whitelist_df_is_none": whitelist_df is None, "whitelist_df_type": str(type(whitelist_df)), "has_empty_attr": hasattr(whitelist_df, "empty") if whitelist_df is not None else False}, "sessionId": "debug-session", "runId": "run1", "hypothesisId": "A"}) + "\n")
    except: pass
    # #endregion
    # Handle case where whitelist_df is a string (should be DataFrame or None)
    # This can happen if modify_config_dict overwrote a DataFrame with a string
    if isinstance(whitelist_df, str):
        whitelist_df = None
    if whitelist_df is not None and not whitelist_df.empty:
        whitelist_map = {}
        for row in whitelist_df.itertuples(index=False):
            whitelist_map.setdefault(row.feature_type, {}).setdefault(row.qualifier, set()).add(row.keyword)
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

    label_filtering["whitelist_map"] = whitelist_map
    label_filtering["priority_map"] = priority_map
    return label_filtering


def get_label_text(feature: SeqFeature, label_filtering: dict) -> str:
    feature_type = feature.type
    qualifiers = feature.qualifiers
    whitelist_map = label_filtering.get("whitelist_map")
    blacklist = label_filtering.get("blacklist_keywords", [])
    priority_map = label_filtering.get("priority_map", {})

    if whitelist_map:
        is_eligible = False
        rules = whitelist_map.get(feature_type, {})
        for key, values in qualifiers.items():
            if key in rules and any(v in rules[key] for v in values):
                is_eligible = True
                break
        if not is_eligible:
            return ""

    priority_list = priority_map.get(
        feature_type, ["product", "gene", "locus_tag", "protein_id", "old_locus_tag", "note"]
    )
    final_label = ""
    for key in priority_list:
        if key in qualifiers and qualifiers[key]:
            final_label = qualifiers[key][0]
            break

    if not whitelist_map and any(bl in final_label.lower() for bl in blacklist):
        return ""

    return final_label


__all__ = [
    "get_label_text",
    "preprocess_label_filtering",
    "read_filter_list_file",
    "read_qualifier_priority_file",
]


