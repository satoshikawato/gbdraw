#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations

from collections.abc import Sequence
from typing import NamedTuple

import numpy as np
import pandas as pd
from pandas import DataFrame, Series
from Bio.Seq import MutableSeq, Seq
from Bio.SeqRecord import SeqRecord

from gbdraw.analysis.windows import SlidingWindows, sliding_windows
from gbdraw.io.record_sequences import known_base_counts_at, sequence_view
from gbdraw.mode_profiles import validate_dinucleotide


_ROWS = {base: row for row, base in enumerate("ACGT")}  # the rows of ``known_base_counts_at``


def _known_counts_at(segments: Sequence[Seq | MutableSeq], positions: np.ndarray) -> np.ndarray:
    """Row ``b`` (A, C, G, T), column ``k``: the count of base ``b`` in the first ``positions[k]`` bases of the
    segments read one after the other (clipped to their joined length), from the owner of record sequences: a
    segment read from a file reads only the rows at its edges and the positions."""

    counts = np.zeros((4, len(positions)), dtype=np.int64)
    offset = 0
    for segment in segments:
        counts += known_base_counts_at(segment, positions - offset)
        offset += len(segment)
    return counts


def _window_pieces(windows: SlidingWindows) -> list[Seq | MutableSeq]:
    """The bases the windows read, one view of the source per source piece, in local order.

    A piece read against the source (a reverse complement) is read as its reverse complement, as the record's own
    bases are (:func:`gbdraw.io.record_sequences.sequence_view`, a chunk at a time).
    """

    sequence = windows.source.sequence
    return [
        sequence_view(sequence, piece.low, piece.low + piece.length, strand=piece.direction)
        for piece in windows.pieces
    ]


def skew_df(record: SeqRecord, window: int, step: int, nt: str, *, circular: bool) -> DataFrame:
    """The dinucleotide content and skew of every GC track, Circular and Linear: one row per sample.

    The samples and their windows are those of :func:`gbdraw.analysis.windows.sliding_windows` (the index holds
    the samples); ``circular`` says whether the molecule is circular. Content is the share of the two bases among
    the known bases (A, C, G, T) the window reads, and skew their difference over their sum; U is counted as T
    (D-26). N and the other IUPAC codes are unknown bases: a window of unknown bases only has no value (NaN) in
    every column, so its tracks break there, and the cumulative skew adds nothing for it (OV-2799).
    """
    nt = validate_dinucleotide(nt)
    counted = nt.replace("U", "T")
    content_legend = f"{nt} content"
    skew_legend = f"{nt} skew"
    cumulative_skew_legend = f"Cumulative {nt} skew, normalized"

    if len(record.seq) == 0 or window <= 0 or step <= 0:
        return pd.DataFrame(
            columns=[content_legend, skew_legend, cumulative_skew_legend],
        )

    windows = sliding_windows(record, window, step, circular=circular)
    segments = _window_pieces(windows)
    bounds = np.concatenate([windows.starts, windows.ends]) - windows.pieces[0].local
    samples = windows.samples
    prefix = _known_counts_at(segments, bounds)
    known_prefix = prefix.sum(axis=0)
    prefix = np.stack((prefix[_ROWS[counted[0]]], prefix[_ROWS[counted[1]]], known_prefix))
    base1_counts, base2_counts, known_counts = prefix[:, len(samples) :] - prefix[:, : len(samples)]
    total_counts = base1_counts + base2_counts
    # A window without a known base has no value (NaN): its tracks break there (OV-2799). A window with known bases
    # but neither of the two has skew 0.
    valued = known_counts > 0
    content_values = np.divide(total_counts, known_counts, out=np.full(len(samples), np.nan), where=valued)
    skew_values = np.divide(
        base1_counts - base2_counts,
        total_counts,
        out=np.where(valued, 0.0, np.nan),
        where=total_counts != 0,
    )
    # The cumulative sum adds nothing for a window without a value, and has none there itself.
    skew_cumulative_values = np.where(valued, np.cumsum(np.where(valued, skew_values, 0.0)), np.nan)

    max_skew_abs = float(np.abs(skew_values[valued]).max()) if valued.any() else 0.0
    max_skew_cumulative_abs = float(np.abs(skew_cumulative_values[valued]).max()) if valued.any() else 0.0
    factor: float = 0.0 if max_skew_cumulative_abs == 0 else (max_skew_abs / max_skew_cumulative_abs)
    normalized_cumulative = skew_cumulative_values * factor
    df = pd.DataFrame(
        {
            content_legend: content_values,
            skew_legend: skew_values,
            cumulative_skew_legend: normalized_cumulative,
        },
        index=samples,
    )
    return df


class DeviationScale(NamedTuple):
    """The baseline and scale of a GC content or skew track drawn as deviation: the mean of the windows with a
    value, and the largest distance of one of them from it."""

    mean: float
    max_diff: float


def deviation_scale(series: Series) -> DeviationScale | None:
    """The deviation baseline and scale of one ``skew_df`` column, read from the windows with a value only (a window
    of unknown bases has none, OV-2799); None when no window has a value. The track draws a window without a value
    on its baseline (the mean)."""

    if not series.notna().any():
        return None
    mean = float(series.mean())
    return DeviationScale(mean, float((series - mean).abs().max()))


__all__ = [
    "DeviationScale",
    "deviation_scale",
    "skew_df",
]


