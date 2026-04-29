#!/usr/bin/env python
# coding: utf-8

"""Depth coverage parsing and binning helpers."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]
from pandas import DataFrame  # type: ignore[reportMissingImports]

from gbdraw.exceptions import ParseError, ValidationError  # type: ignore[reportMissingImports]

DEPTH_COLUMNS = ["reference_name", "position", "depth"]


def _is_int_like(value: str) -> bool:
    try:
        int(str(value).strip())
    except ValueError:
        return False
    return True


def _is_float_like(value: str) -> bool:
    try:
        float(str(value).strip())
    except ValueError:
        return False
    return True


def _first_line_has_header(path: Path) -> bool:
    try:
        with path.open("rt", encoding="utf-8") as handle:
            first_line = handle.readline()
    except OSError as exc:
        raise ParseError(f"Unable to read depth file '{path}': {exc}") from exc

    if not first_line:
        return False
    fields = first_line.rstrip("\n\r").split("\t")
    if len(fields) < 3:
        raise ParseError(
            f"Depth file '{path}' must contain at least 3 tab-separated columns."
        )
    return not (_is_int_like(fields[1]) and _is_float_like(fields[2]))


def read_depth_tsv(path: str) -> DataFrame:
    """Read a samtools-depth-like TSV file.

    The expected columns are reference name, 1-based position, and non-negative
    depth. Headerless files are expected; a single header line is tolerated when
    the position/depth fields on the first line are non-numeric.
    """

    depth_path = Path(path)
    if not depth_path.exists():
        raise ParseError(f"Depth file does not exist: {path}")
    if not depth_path.is_file():
        raise ParseError(f"Depth path is not a file: {path}")

    skiprows = 1 if _first_line_has_header(depth_path) else 0
    try:
        depth_table = pd.read_csv(
            depth_path,
            sep="\t",
            header=None,
            names=DEPTH_COLUMNS,
            usecols=[0, 1, 2],
            skiprows=skiprows,
            dtype={
                "reference_name": "string",
                "position": "int64",
                "depth": "float64",
            },
        )
    except ValueError as exc:
        raise ParseError(f"Unable to parse depth file '{path}': {exc}") from exc

    _validate_depth_table(depth_table)
    return depth_table


def _normalize_depth_table(depth_table: DataFrame) -> DataFrame:
    """Return a DataFrame with canonical depth columns."""

    if set(DEPTH_COLUMNS).issubset(depth_table.columns):
        normalized = depth_table.loc[:, DEPTH_COLUMNS].copy()
    elif len(depth_table.columns) >= 3:
        normalized = depth_table.iloc[:, [0, 1, 2]].copy()
        normalized.columns = DEPTH_COLUMNS
    else:
        raise ValidationError("Depth table must contain at least 3 columns.")

    normalized["reference_name"] = normalized["reference_name"].astype("string")
    normalized["position"] = pd.to_numeric(normalized["position"], errors="raise").astype("int64")
    normalized["depth"] = pd.to_numeric(normalized["depth"], errors="raise").astype("float64")
    _validate_depth_table(normalized)
    return normalized


def _validate_depth_table(depth_table: DataFrame) -> None:
    if depth_table.empty:
        return
    if depth_table["reference_name"].isna().any():
        raise ValidationError("Depth table contains missing reference names.")
    if (depth_table["position"] <= 0).any():
        raise ValidationError("Depth positions must be 1-based positive integers.")
    if (depth_table["depth"] < 0).any():
        raise ValidationError("Depth values must be non-negative.")


def _unique_references(depth_table: DataFrame) -> list[str]:
    values: Iterable[object] = depth_table["reference_name"].dropna().unique()
    return [str(value) for value in values]


def _select_record_depth_rows(record: SeqRecord, depth_table: DataFrame) -> DataFrame:
    record_id = str(record.id)
    exact_rows = depth_table[depth_table["reference_name"] == record_id]
    if not exact_rows.empty:
        return exact_rows

    references = _unique_references(depth_table)
    if len(references) == 1:
        return depth_table[depth_table["reference_name"] == references[0]]

    raise ValidationError(
        "Depth file references do not match record "
        f"'{record_id}'. Available references: {', '.join(references)}"
    )


def _validate_depth_range(min_depth: float | None, max_depth: float | None) -> None:
    if min_depth is not None and min_depth < 0:
        raise ValidationError("min_depth must be >= 0.")
    if max_depth is not None and max_depth < 0:
        raise ValidationError("max_depth must be >= 0.")
    if min_depth is not None and max_depth is not None and min_depth > max_depth:
        raise ValidationError("min_depth must be <= max_depth.")


def depth_df(
    record: SeqRecord,
    depth_table: DataFrame,
    window: int,
    step: int,
    *,
    normalize: bool = False,
    min_depth: float | None = None,
    max_depth: float | None = None,
) -> DataFrame:
    """Bin depth values for a record into fixed genomic windows."""

    if window <= 0:
        raise ValidationError("window must be > 0.")
    if step <= 0:
        raise ValidationError("step must be > 0.")
    _validate_depth_range(min_depth, max_depth)

    normalized_table = _normalize_depth_table(depth_table)
    record_len = len(record.seq)
    if record_len <= 0:
        return pd.DataFrame(columns=["position", "depth", "depth_normalized"])

    rows = _select_record_depth_rows(record, normalized_table)
    rows = rows[(rows["position"] >= 1) & (rows["position"] <= record_len)]

    coverage = np.zeros(record_len, dtype=np.float64)
    if not rows.empty:
        zero_based = rows["position"].to_numpy(dtype=np.int64) - 1
        depths = rows["depth"].to_numpy(dtype=np.float64)
        summed_depths = np.bincount(zero_based, weights=depths, minlength=record_len)
        counts = np.bincount(zero_based, minlength=record_len)
        nonzero = counts > 0
        coverage[nonzero] = summed_depths[nonzero] / counts[nonzero]

    prefix = np.concatenate(([0.0], np.cumsum(coverage)))
    starts = np.arange(0, record_len, int(step), dtype=np.int64)
    ends = np.minimum(starts + int(window), record_len)
    bin_lengths = np.maximum(ends - starts, 1)
    means = (prefix[ends] - prefix[starts]) / bin_lengths

    if min_depth is not None or max_depth is not None:
        lower = float(min_depth) if min_depth is not None else None
        upper = float(max_depth) if max_depth is not None else None
        means = np.clip(means, lower, upper)

    if normalize:
        # Log-scale depth follows IGV-style coverage scaling:
        # 1x, 10x, and 100x occupy evenly spaced heights.
        normalization_min = float(min_depth) if min_depth is not None else 0.0
        normalization_max = float(max_depth) if max_depth is not None else float(np.max(means, initial=0.0))

        def transform(values: np.ndarray | float) -> np.ndarray:
            values_arr = np.asarray(values, dtype=np.float64)
            transformed_arr = np.zeros_like(values_arr, dtype=np.float64)
            positive = values_arr > 0
            transformed_arr[positive] = np.log10(values_arr[positive]) + 1.0
            return transformed_arr

        transformed = transform(means)
        transformed_min = float(transform(normalization_min))
        transformed_max = float(transform(normalization_max))
        if transformed_max <= transformed_min:
            normalized = np.zeros_like(means, dtype=np.float64)
        else:
            normalized = (transformed - transformed_min) / (transformed_max - transformed_min)
            normalized = np.clip(normalized, 0.0, 1.0)
    else:
        normalized = means.astype(np.float64, copy=True)

    return pd.DataFrame(
        {
            "position": starts,
            "depth": means,
            "depth_normalized": normalized,
        }
    )


__all__ = ["DEPTH_COLUMNS", "depth_df", "read_depth_tsv"]
