#!/usr/bin/env python
# coding: utf-8

"""GC content plotting helpers."""

from __future__ import annotations

import numpy as np
import pandas as pd
from Bio.SeqRecord import SeqRecord  # type: ignore[reportMissingImports]
from pandas import DataFrame  # type: ignore[reportMissingImports]

from gbdraw.exceptions import ValidationError  # type: ignore[reportMissingImports]
from gbdraw.analysis.skew import _build_prefix_counts, _window_count_from_prefix  # type: ignore[reportMissingImports]


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


def circular_dinucleotide_content_df(
    record: SeqRecord,
    window: int,
    step: int,
    nt: str,
) -> DataFrame:
    """Return centered circular dinucleotide-content data for circular tracks."""

    nt_list = list(str(nt).upper())
    nt_1 = nt_list[0]
    nt_2 = nt_list[1]
    seq_str = str(record.seq).upper()
    seq_length = len(seq_str)
    content_legend = f"{nt} content"

    if seq_length == 0 or window <= 0 or step <= 0:
        return pd.DataFrame(columns=[content_legend])

    seq_bytes = seq_str.encode("ascii")
    prefix_nt_1 = _build_prefix_counts(seq_bytes, ord(nt_1))
    prefix_nt_2 = _build_prefix_counts(seq_bytes, ord(nt_2))
    total_nt_1 = prefix_nt_1[-1]
    total_nt_2 = prefix_nt_2[-1]
    half_window = int(window) // 2
    window_length = float(window)

    starts: list[int] = []
    content_values: list[float] = []
    for position in range(0, seq_length, step):
        circular_start = (int(position) - half_window) % seq_length
        nt_1_count = _window_count_from_prefix(
            prefix_nt_1,
            seq_length,
            start=circular_start,
            window=window,
            total_count=total_nt_1,
        )
        nt_2_count = _window_count_from_prefix(
            prefix_nt_2,
            seq_length,
            start=circular_start,
            window=window,
            total_count=total_nt_2,
        )
        starts.append(position)
        content_values.append((nt_1_count + nt_2_count) / window_length)

    return pd.DataFrame({content_legend: content_values}, index=starts)


__all__ = ["calculate_gc_percent", "circular_dinucleotide_content_df", "gc_content_percent_df"]
