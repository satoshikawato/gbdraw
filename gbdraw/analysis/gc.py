#!/usr/bin/env python
# coding: utf-8

"""GC content plotting helpers."""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np
import pandas as pd
from pandas import DataFrame

from gbdraw.exceptions import ValidationError
from Bio.Seq import MutableSeq, Seq

from gbdraw.io.record_sequences import known_base_counts


def calculate_gc_percent(seq: Seq | MutableSeq | bytes) -> float:
    """Calculate GC percentage for a sequence, counted as the GC track counts: over its known bases (A, C, G, T),
    0 when it has none (OV-2799). The counts come from the owner of record sequences
    (:func:`gbdraw.io.record_sequences.known_base_counts`): a record read from a file reads no base again."""

    return gc_percent_of_sequences((seq,))


def gc_percent_of_sequences(seqs: Iterable[Seq | MutableSeq | bytes]) -> float:
    """The GC percentage of several sequences taken together, counted as the GC track counts: over their known bases
    (A, C, G, T), 0 when they have none (OV-2799). Each sequence's counts come from
    :func:`gbdraw.io.record_sequences.known_base_counts`, as for one sequence."""

    gc_count = 0
    known = 0
    for seq in seqs:
        a, c, g, t = (int(count) for count in known_base_counts(seq))
        gc_count += g + c
        known += a + c + g + t
    if not known:
        return 0.0
    return round(100.0 * gc_count / known, 2)


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
    computed percentage. A window without a value (unknown bases only) keeps
    NaN there and draws on the axis minimum (``value_normalized`` 0), so the
    area breaks; a frame with no value at all draws nothing (OV-2799).
    """

    column = f"{dinucleotide} content"
    columns = ["position", "percent", "value", "value_normalized"]
    if column not in gc_df.columns and not gc_df.empty:
        raise ValidationError(f"gc_content data frame must include '{column}'.")
    if gc_df.empty or not gc_df[column].notna().any():
        return pd.DataFrame(columns=columns)

    axis_min = float(min_percent)
    axis_max = float(max_percent)
    if axis_min > axis_max:
        raise ValidationError("gc_content_min_percent must be <= gc_content_max_percent")

    percent_values = pd.to_numeric(gc_df[column], errors="raise").astype(float) * 100.0
    if axis_max <= axis_min:
        normalized = np.zeros(len(percent_values), dtype=np.float64)
    else:
        clipped = np.clip(percent_values.to_numpy(dtype=np.float64), axis_min, axis_max)
        normalized = np.nan_to_num((clipped - axis_min) / (axis_max - axis_min), nan=0.0)

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
