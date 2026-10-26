#!/usr/bin/env python
# coding: utf-8

"""GC content plotting helpers."""

from __future__ import annotations

import numpy as np
import pandas as pd
from Bio.SeqRecord import SeqRecord
from pandas import DataFrame

from gbdraw.exceptions import ValidationError
from gbdraw.analysis.skew import (
    counted_dinucleotide,
    counted_sequence,
    window_base_counts,
)


def calculate_gc_percent(seq: object) -> float:
    """Calculate GC percentage for a sequence-like object, counted as the GC track counts."""

    sequence = counted_sequence(seq)
    if not sequence:
        return 0.0
    gc_count = sequence.count(b"G") + sequence.count(b"C")
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


def circular_dinucleotide_content_df(
    record: SeqRecord,
    window: int,
    step: int,
    nt: str,
) -> DataFrame:
    """Return centered circular dinucleotide-content data for circular tracks."""

    nt, nt_1, nt_2, seq_bytes = counted_dinucleotide(record, nt)
    seq_length = len(seq_bytes)
    content_legend = f"{nt} content"

    if seq_length == 0 or window <= 0 or step <= 0:
        return pd.DataFrame(columns=[content_legend])

    starts = np.arange(0, seq_length, step, dtype=np.int64)
    circular_starts = (starts - int(window) // 2) % seq_length
    total_counts = window_base_counts(
        seq_bytes, ord(nt_1), starts=circular_starts, window=window
    ) + window_base_counts(seq_bytes, ord(nt_2), starts=circular_starts, window=window)
    content_values = total_counts / float(window)

    return pd.DataFrame({content_legend: content_values}, index=starts)


__all__ = ["calculate_gc_percent", "circular_dinucleotide_content_df", "gc_content_percent_df"]
