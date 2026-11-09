#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd
from pandas import DataFrame
from Bio.Seq import MutableSeq, Seq
from Bio.SeqRecord import SeqRecord

from gbdraw.analysis.windows import SlidingWindows, sliding_windows
from gbdraw.mode_profiles import validate_dinucleotide


# The sequence is read in chunks of this many bases, so the mask and the
# occurrence offsets of one chunk are the only arrays beside the positions, and
# every offset fits a 32-bit ``intp`` (Pyodide).
_COUNT_CHUNK = 1 << 20


def _prefix_counts_at(segments: Sequence[bytes], base: int, positions: np.ndarray) -> np.ndarray:
    """``b"".join(segments)[:p].count(base)`` for each position ``p`` (0 to the joined length).

    The segments are read once, a chunk at a time and never joined: the
    positions are sorted once, each chunk takes the positions that fall in it,
    reads their counts from the chunk's own occurrences, and adds the total of
    the chunks before it. The extra memory is the positions plus one chunk,
    never an array as long as the sequence.
    """
    length = sum(len(segment) for segment in segments)
    order = np.argsort(positions, kind="stable")
    wanted = np.clip(positions[order], 0, length)
    counts = np.empty(len(wanted), dtype=np.int64)
    running = 0
    offset = 0
    for segment in segments:
        sequence = np.frombuffer(segment, dtype=np.uint8)
        for start in range(0, len(sequence), _COUNT_CHUNK):
            chunk = sequence[start : start + _COUNT_CHUNK]
            first = offset + start
            lo, hi = np.searchsorted(wanted, (first, first + len(chunk)), side="left")
            chunk_hits = np.flatnonzero(chunk == base)
            local = (wanted[lo:hi] - first).astype(np.intp)
            counts[lo:hi] = running + np.searchsorted(chunk_hits, local, side="left")
            running += len(chunk_hits)
        offset += len(sequence)
    counts[np.searchsorted(wanted, length, side="left") :] = running
    prefix = np.empty(len(positions), dtype=np.int64)
    prefix[order] = counts
    return prefix


_COUNTED_BASES = bytes.maketrans(b"abcdefghijklmnopqrstuvwxyzUu", b"ABCDEFGHIJKLMNOPQRSTUVWXYZTT")


def counted_sequence(seq: object) -> bytes:
    """The bases as every GC and skew count reads them, as one bytes copy.

    Case is ignored, so soft-masked FASTA bases count, and U is counted as T
    (D-26). The copy is made once by ``bytes.translate``; callers count it
    directly and do not encode a str of it.
    """

    raw = bytes(seq) if isinstance(seq, (Seq, MutableSeq, bytes, bytearray)) else str(seq).encode("ascii")
    return raw.translate(_COUNTED_BASES)


def _counted_bases(windows: SlidingWindows) -> list[bytes]:
    """The counted bases the windows read, one segment per source piece, in local order.

    A piece read against the source (a reverse complement) is read as its reverse complement, as the record's own
    bases are. A piece that is the whole source is counted without a slice copy.
    """

    sequence = windows.source.sequence
    segments: list[bytes] = []
    for piece in windows.pieces:
        bases = sequence if piece.length == len(sequence) else sequence[piece.low : piece.low + piece.length]
        segments.append(counted_sequence(bases if piece.direction > 0 else bases.reverse_complement()))
    return segments


def skew_df(record: SeqRecord, window: int, step: int, nt: str, *, circular: bool) -> DataFrame:
    """The dinucleotide content and skew of every GC track, Circular and Linear: one row per sample.

    The samples and their windows are those of :func:`gbdraw.analysis.windows.sliding_windows` (the index holds
    the samples); ``circular`` says whether the molecule is circular. Content is the share of the two bases among
    the bases a window reads (``window``, or fewer at the ends of a linear record). U is counted as T (D-26).
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
    segments = _counted_bases(windows)
    bounds = np.concatenate([windows.starts, windows.ends]) - windows.pieces[0].local
    samples = windows.samples
    base1_prefix = _prefix_counts_at(segments, ord(counted[0]), bounds)
    base2_prefix = _prefix_counts_at(segments, ord(counted[1]), bounds)
    base1_counts = base1_prefix[len(samples) :] - base1_prefix[: len(samples)]
    base2_counts = base2_prefix[len(samples) :] - base2_prefix[: len(samples)]
    total_counts = base1_counts + base2_counts
    skew_values = np.divide(
        base1_counts - base2_counts,
        total_counts,
        out=np.zeros(len(samples), dtype=np.float64),
        where=total_counts != 0,
    )
    content_values = total_counts / windows.bases
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
        index=samples,
    )
    return df


__all__ = [
    "counted_sequence",
    "skew_df",
]


