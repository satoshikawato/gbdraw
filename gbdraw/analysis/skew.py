#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations

from collections.abc import Iterator, Sequence
from typing import NamedTuple

import numpy as np
import pandas as pd
from pandas import DataFrame, Series
from Bio.Seq import MutableSeq, Seq
from Bio.SeqRecord import SeqRecord

from gbdraw.analysis.windows import SlidingWindows, sliding_windows
from gbdraw.core.sequence_spans import sequence_spans
from gbdraw.io.record_sequences import sequence_view
from gbdraw.mode_profiles import validate_dinucleotide


# The sequence is read in chunks of this many bases, so the mask and the
# occurrence offsets of one chunk are the only arrays beside the positions, and
# every offset fits a 32-bit ``intp`` (Pyodide).
_COUNT_CHUNK = 1 << 20


# Multiplying a word of eight 0/1 bytes by this constant sums them into its top byte.
_BYTE_SUM = np.uint64(0x0101010101010101)
_TOP_BYTE = np.uint64(56)
_EIGHT = np.arange(8)


def _prefix_counts_at(
    segments: Sequence[Seq | MutableSeq | bytes], bases: Sequence[int], positions: np.ndarray
) -> np.ndarray:
    """Row ``i``, column ``k``: the count of ``bases[i]`` in the first ``p`` counted bases (:func:`counted_bases`) of
    the segments read one after the other, for each position ``p = positions[k]`` (0 to their joined length).

    The segments are read once for every base together, a chunk at a time and never joined: the positions are
    sorted once, each chunk takes the positions that fall in it, and each base's mask of the chunk is read as
    8-byte words whose set bytes one multiplication sums, so a position's count is the words before it, the bytes
    of its own word before it, and the total of the chunks before; the cost does not depend on how often the base
    occurs (a run of N fills whole chunks). The mask and the word sums are buffers of one chunk, reused for every
    chunk and base, never an array as long as the sequence.
    """
    length = sum(len(segment) for segment in segments)
    order = np.argsort(positions, kind="stable")
    wanted = np.clip(positions[order], 0, length)
    counts = np.empty((len(bases), len(wanted)), dtype=np.int64)
    running = np.zeros(len(bases), dtype=np.int64)
    size = -(-min(_COUNT_CHUNK, max((len(segment) for segment in segments), default=0)) // 8) * 8
    mask = np.zeros(size, dtype=bool)
    words = np.empty(size // 8, dtype=np.uint64)
    before = np.zeros(size // 8 + 1, dtype=np.int64)
    offset = 0
    for segment in segments:
        start = 0
        for counted in counted_bases(segment):
            chunk = np.frombuffer(counted, dtype=np.uint8)
            first = offset + start
            start += len(chunk)
            lo, hi = np.searchsorted(wanted, (first, first + len(chunk)), side="left")
            word, place = np.divmod((wanted[lo:hi] - first).astype(np.intp), 8)
            in_word = _EIGHT < place[:, None]
            count = -(-len(chunk) // 8)
            chunk_mask, chunk_words = mask[: 8 * count], words[:count]
            chunk_mask[len(chunk) :] = False
            for row, base in enumerate(bases):
                np.equal(chunk, base, out=chunk_mask[: len(chunk)])
                np.multiply(chunk_mask.view(np.uint64), _BYTE_SUM, out=chunk_words)
                np.right_shift(chunk_words, _TOP_BYTE, out=chunk_words)
                np.cumsum(chunk_words.view(np.int64), out=before[1 : count + 1])
                partial = np.count_nonzero(chunk_mask.view(np.uint8).reshape(-1, 8)[word] & in_word, axis=1)
                counts[row, lo:hi] = running[row] + before[word] + partial
                running[row] += before[count]
            del chunk, counted  # one chunk is held at a time, also while the next one is read
        offset += len(segment)
    counts[:, np.searchsorted(wanted, length, side="left") :] = running[:, None]
    prefix = np.empty_like(counts)
    prefix[:, order] = counts
    return prefix


def _counting_table() -> bytes:
    """The ``bytes.translate`` table of :func:`counted_bases`: A, C, G and T in either case read as themselves
    in upper case, U as T, and every other byte as N."""
    table = bytearray(b"N" * 256)
    for byte, counted in zip(b"ACGTUacgtu", b"ACGTTACGTT"):
        table[byte] = counted
    return bytes(table)


_COUNTED_BASES = _counting_table()
_UNKNOWN_BASE = ord("N")


def counted_bases(seq: Seq | MutableSeq | bytes) -> Iterator[bytes]:
    """The bases as every GC and skew count reads them, ``_COUNT_CHUNK`` at a time, in order.

    Case is ignored, so soft-masked FASTA bases count, and U is counted as T
    (D-26). The known bases are A, C, G and T; every other byte (N, the other
    IUPAC codes, a gap) is an unknown base and reads as N, which no content or
    skew value counts (OV-2799). Each chunk is a slice of ``seq`` read through
    ``bytes.translate``, so a count holds one chunk, never a copy of the sequence.
    """

    # ``map`` keeps no chunk it has handed on, so a reader holds one chunk and the next one's read only.
    return map(_count_translated, sequence_spans(seq, _COUNT_CHUNK))


def _count_translated(span: bytes) -> bytes:
    return span.translate(_COUNTED_BASES)


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
    prefix = _prefix_counts_at(segments, (ord(counted[0]), ord(counted[1]), _UNKNOWN_BASE), bounds)
    base1_counts, base2_counts, unknown_counts = prefix[:, len(samples) :] - prefix[:, : len(samples)]
    known_counts = windows.bases - unknown_counts
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
    "counted_bases",
    "deviation_scale",
    "skew_df",
]


