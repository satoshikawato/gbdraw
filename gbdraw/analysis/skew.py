#!/usr/bin/env python
# coding: utf-8

from __future__ import annotations

from typing import Any, Generator

import pandas as pd
from pandas import DataFrame
from Bio.SeqRecord import SeqRecord


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


def _build_prefix_counts(seq_bytes: bytes, target_base: int) -> list[int]:
    """Return prefix counts for a single nucleotide across the sequence."""
    prefix = [0] * (len(seq_bytes) + 1)
    running_count = 0
    for idx, value in enumerate(seq_bytes, start=1):
        if value == target_base:
            running_count += 1
        prefix[idx] = running_count
    return prefix


def _window_count_from_prefix(
    prefix: list[int],
    seq_length: int,
    *,
    start: int,
    window: int,
    total_count: int,
) -> int:
    """Count a base in a circular window using prefix sums."""
    if seq_length <= 0 or window <= 0:
        return 0

    full_cycles, remainder = divmod(window, seq_length)
    count = full_cycles * total_count
    if remainder == 0:
        return count

    end = start + remainder
    if end <= seq_length:
        count += prefix[end] - prefix[start]
    else:
        count += (prefix[seq_length] - prefix[start]) + prefix[end - seq_length]
    return count


def skew_df(record: SeqRecord, window: int, step: int, nt: str) -> DataFrame:
    """
    Calculates dinucleotide skew and content in a DNA sequence, returning a DataFrame.
    """
    nt_list = list(str(nt).upper())
    nt_1: str = nt_list[0]
    nt_2: str = nt_list[1]
    seq_str = str(record.seq).upper()
    seq_length = len(seq_str)
    content_legend = f"{nt} content"
    skew_legend = f"{nt} skew"
    cumulative_skew_legend = f"Cumulative {nt} skew, normalized"

    if seq_length == 0 or window <= 0 or step <= 0:
        return pd.DataFrame(
            columns=[content_legend, skew_legend, cumulative_skew_legend],
        )

    seq_bytes = seq_str.encode("ascii")
    nt_1_byte = ord(nt_1)
    nt_2_byte = ord(nt_2)
    prefix_nt_1 = _build_prefix_counts(seq_bytes, nt_1_byte)
    prefix_nt_2 = _build_prefix_counts(seq_bytes, nt_2_byte)
    total_nt_1 = prefix_nt_1[-1]
    total_nt_2 = prefix_nt_2[-1]

    starts: list[int] = []
    content_values: list[float] = []
    skew_values: list[float] = []
    skew_cumulative_values: list[float] = []
    skew_sum = 0.0
    window_length = float(window)

    for start in range(0, seq_length, step):
        base1_count = _window_count_from_prefix(
            prefix_nt_1,
            seq_length,
            start=start,
            window=window,
            total_count=total_nt_1,
        )
        base2_count = _window_count_from_prefix(
            prefix_nt_2,
            seq_length,
            start=start,
            window=window,
            total_count=total_nt_2,
        )
        total_count = base1_count + base2_count
        if total_count == 0:
            skew = 0.0
        else:
            skew = (base1_count - base2_count) / total_count

        starts.append(start)
        content_values.append(total_count / window_length)
        skew_values.append(skew)
        skew_sum += skew
        skew_cumulative_values.append(skew_sum)

    max_skew_abs = max((abs(value) for value in skew_values), default=0.0)
    max_skew_cumulative_abs = max((abs(value) for value in skew_cumulative_values), default=0.0)
    factor: float = 0.0 if max_skew_cumulative_abs == 0 else (max_skew_abs / max_skew_cumulative_abs)
    normalized_cumulative = [value * factor for value in skew_cumulative_values]
    df = pd.DataFrame(
        {
            content_legend: content_values,
            skew_legend: skew_values,
            cumulative_skew_legend: normalized_cumulative,
        },
        index=starts,
    )
    return df


__all__ = ["calculate_dinucleotide_skew", "skew_df", "sliding_window"]


