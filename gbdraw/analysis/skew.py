#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations

from typing import Any, Generator

import numpy as np
import pandas as pd
from pandas import DataFrame
from Bio.Seq import MutableSeq, Seq
from Bio.SeqRecord import SeqRecord

from gbdraw.mode_profiles import validate_dinucleotide


def calculate_dinucleotide_skew(seq: str, base1: str, base2: str) -> float:
    """
    Calculate the skew of two nucleotides in a DNA sequence.

    Skew = (Count(Base1) - Count(Base2)) / (Count(Base1) + Count(Base2)).
    """
    base1_count: int = seq.count(base1)
    base2_count: int = seq.count(base2)
    total_count = base1_count + base2_count

    if total_count == 0:
        return 0.0

    skew: float = (base1_count - base2_count) / total_count
    return skew


def sliding_window(seq: str, window: int, step: int) -> Generator[tuple[int, str], Any, None]:
    """
    Generate substrings from a given sequence using a circular sliding window.
    """
    for start in range(0, len(seq), step):
        end: int = start + window
        if end > len(seq):
            overhang_length: int = end - len(seq)
            # assuming circular sequence
            out_seq: str = seq[start : len(seq)] + seq[0:overhang_length]
        else:
            out_seq = seq[start:end]
        yield start, out_seq


# The sequence is read in chunks of this many bases, so the mask and the
# occurrence offsets of one chunk are the only arrays beside the positions, and
# every offset fits a 32-bit ``intp`` (Pyodide).
_COUNT_CHUNK = 1 << 20


def _prefix_counts_at(seq_bytes: bytes, base: int, positions: np.ndarray) -> np.ndarray:
    """``seq_bytes[:p].count(base)`` for each position ``p`` (0 to n).

    The sequence is read once, a chunk at a time: the positions are sorted
    once, each chunk takes the positions that fall in it, reads their counts
    from the chunk's own occurrences, and adds the total of the chunks before
    it. The extra memory is the positions plus one chunk, never an array as
    long as the sequence.
    """
    sequence = np.frombuffer(seq_bytes, dtype=np.uint8)
    length = len(sequence)
    order = np.argsort(positions, kind="stable")
    wanted = np.clip(positions[order], 0, length)
    counts = np.empty(len(wanted), dtype=np.int64)
    running = 0
    for start in range(0, length, _COUNT_CHUNK):
        lo, hi = np.searchsorted(wanted, (start, min(start + _COUNT_CHUNK, length)), side="left")
        chunk_hits = np.flatnonzero(sequence[start : start + _COUNT_CHUNK] == base)
        local = (wanted[lo:hi] - start).astype(np.intp)
        counts[lo:hi] = running + np.searchsorted(chunk_hits, local, side="left")
        running += len(chunk_hits)
    counts[np.searchsorted(wanted, length, side="left") :] = running
    prefix = np.empty(len(positions), dtype=np.int64)
    prefix[order] = counts
    return prefix


def window_base_counts(
    seq_bytes: bytes,
    base: int,
    *,
    starts: np.ndarray,
    window: int,
) -> np.ndarray:
    """Count ``base`` in the circular windows of ``window`` bases beginning at
    ``starts``: the one window count of the GC content and skew tracks.

    A window longer than the sequence counts it once per full cycle; the
    counts are read at the window bounds only (``_prefix_counts_at``).
    """
    seq_length = len(seq_bytes)
    full_cycles, remainder = divmod(window, seq_length)
    end = starts + remainder
    wrapped = end > seq_length
    stops = np.where(wrapped, end - seq_length, end)
    prefix = _prefix_counts_at(seq_bytes, base, np.concatenate([starts, stops, [seq_length]]))
    total_count = int(prefix[-1])
    counts = np.full(len(starts), full_cycles * total_count, dtype=np.int64)
    if remainder == 0:
        return counts
    counts += prefix[len(starts) : 2 * len(starts)] - prefix[: len(starts)]
    counts[wrapped] += total_count
    return counts


_COUNTED_BASES = bytes.maketrans(b"abcdefghijklmnopqrstuvwxyzUu", b"ABCDEFGHIJKLMNOPQRSTUVWXYZTT")


def counted_sequence(seq: object) -> bytes:
    """The bases as every GC and skew count reads them, as one bytes copy.

    Case is ignored, so soft-masked FASTA bases count, and U is counted as T
    (D-26). The copy is made once by ``bytes.translate``; callers count it
    directly and do not encode a str of it.
    """

    raw = bytes(seq) if isinstance(seq, (Seq, MutableSeq, bytes, bytearray)) else str(seq).encode("ascii")
    return raw.translate(_COUNTED_BASES)


def counted_dinucleotide(record: SeqRecord, nt: str) -> tuple[str, str, str, bytes]:
    """Return the display pair, the two counted bases, and the counted sequence.

    U is the same base as T (D-26): it is counted as T in both the pair and the
    sequence, while the display name keeps the requested letters.
    """

    pair = validate_dinucleotide(nt)
    counted = pair.replace("U", "T")
    return pair, counted[0], counted[1], counted_sequence(record.seq)


def skew_df(record: SeqRecord, window: int, step: int, nt: str) -> DataFrame:
    """
    Calculates dinucleotide skew and content in a DNA sequence, returning a DataFrame.
    """
    nt, nt_1, nt_2, seq_bytes = counted_dinucleotide(record, nt)
    seq_length = len(seq_bytes)
    content_legend = f"{nt} content"
    skew_legend = f"{nt} skew"
    cumulative_skew_legend = f"Cumulative {nt} skew, normalized"

    if seq_length == 0 or window <= 0 or step <= 0:
        return pd.DataFrame(
            columns=[content_legend, skew_legend, cumulative_skew_legend],
        )

    starts = np.arange(0, seq_length, step, dtype=np.int64)
    base1_counts = window_base_counts(seq_bytes, ord(nt_1), starts=starts, window=window)
    base2_counts = window_base_counts(seq_bytes, ord(nt_2), starts=starts, window=window)
    total_counts = base1_counts + base2_counts
    skew_values = np.divide(
        base1_counts - base2_counts,
        total_counts,
        out=np.zeros(len(starts), dtype=np.float64),
        where=total_counts != 0,
    )
    content_values = total_counts / float(window)
    skew_cumulative_values = np.cumsum(skew_values)

    max_skew_abs = float(np.abs(skew_values).max())
    max_skew_cumulative_abs = float(np.abs(skew_cumulative_values).max())
    factor: float = 0.0 if max_skew_cumulative_abs == 0 else (max_skew_abs / max_skew_cumulative_abs)
    normalized_cumulative = skew_cumulative_values * factor
    df = pd.DataFrame(
        {
            content_legend: content_values,
            skew_legend: skew_values,
            cumulative_skew_legend: normalized_cumulative,
        },
        index=starts,
    )
    return df


__all__ = [
    "calculate_dinucleotide_skew",
    "counted_dinucleotide",
    "counted_sequence",
    "skew_df",
    "sliding_window",
    "window_base_counts",
]


