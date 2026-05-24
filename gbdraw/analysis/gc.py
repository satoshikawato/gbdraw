#!/usr/bin/env python
# coding: utf-8

"""GC content plotting helpers."""

from __future__ import annotations

import numpy as np
import pandas as pd
from pandas import DataFrame  # type: ignore[reportMissingImports]

from gbdraw.exceptions import ValidationError  # type: ignore[reportMissingImports]


def calculate_gc_percent(seq: object) -> float:
    """Calculate GC percentage for a sequence-like object."""

    sequence = str(seq)
    if not sequence:
        return 0.0
    gc_count = sequence.count("G") + sequence.count("C")
    return round(100.0 * gc_count / len(sequence), 2)


def gc_content_percent_df(
    gc_df: DataFrame,
    *,
    dinucleotide: str,
    min_percent: float = 0.0,
    max_percent: float = 100.0,
) -> DataFrame:
    """Return scalar plot data for absolute GC percentage rendering.

    The source GC percentage is preserved in ``percent``/``value``. Clipping is
    applied only to ``value_normalized`` so callers can inspect the unclipped
    computed percentage.
    """

    column = f"{dinucleotide} content"
    columns = ["position", "percent", "value", "value_normalized"]
    if gc_df.empty:
        return pd.DataFrame(columns=columns)
    if column not in gc_df.columns:
        raise ValidationError(f"gc_content data frame must include '{column}'.")

    axis_min = float(min_percent)
    axis_max = float(max_percent)
    if axis_min > axis_max:
        raise ValidationError("gc_content_min_percent must be <= gc_content_max_percent")

    percent_values = pd.to_numeric(gc_df[column], errors="raise").astype(float) * 100.0
    if axis_max <= axis_min:
        normalized = np.zeros(len(percent_values), dtype=np.float64)
    else:
        clipped = np.clip(percent_values.to_numpy(dtype=np.float64), axis_min, axis_max)
        normalized = (clipped - axis_min) / (axis_max - axis_min)

    if "position" in gc_df.columns:
        positions = pd.to_numeric(gc_df["position"], errors="raise").astype(float)
    else:
        positions = pd.Series(gc_df.index, index=gc_df.index, dtype="float64")

    return pd.DataFrame(
        {
            "position": positions.to_numpy(dtype=np.float64),
            "percent": percent_values.to_numpy(dtype=np.float64),
            "value": percent_values.to_numpy(dtype=np.float64),
            "value_normalized": normalized,
        },
        index=gc_df.index,
    )


__all__ = ["calculate_gc_percent", "gc_content_percent_df"]
