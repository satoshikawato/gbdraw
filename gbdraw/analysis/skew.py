#!/usr/bin/env python
# coding: utf-8

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


def skew_df(record: SeqRecord, window: int, step: int, nt: str) -> DataFrame:
    """
    Calculates dinucleotide skew and content in a DNA sequence, returning a DataFrame.
    """
    nt_list = list(nt)
    nt_1: str = nt_list[0]
    nt_2: str = nt_list[1]
    skew_sum: float = 0
    skew_dict: dict[int, float] = {}
    content_dict: dict[int, float] = {}
    skew_cumulative_dict: dict[int, float] = {}
    seq: str = record.seq.upper()
    for start, seq_part in sliding_window(seq, window, step):
        skew: float = calculate_dinucleotide_skew(seq_part, nt_1, nt_2)
        dinucleotide_content: float = (seq_part.count(nt_1) + seq_part.count(nt_2)) / len(
            seq_part
        )
        content_dict[start] = dinucleotide_content
        skew_dict[start] = skew
        skew_sum = skew_sum + skew
        skew_cumulative_dict[start] = skew_sum
    max_skew_abs: float = abs(max(skew_dict.values(), key=abs))
    max_skew_cumulative_abs: float = abs(max(skew_cumulative_dict.values(), key=abs))
    # If the cumulative series is all zeros, keep it as all zeros (avoid 0/0).
    factor: float = 0.0 if max_skew_cumulative_abs == 0 else (max_skew_abs / max_skew_cumulative_abs)
    skew_cumulative_dict.update((x, (y * factor)) for x, y in skew_cumulative_dict.items())
    content_legend: str = "{} content".format(nt)
    skew_legend: str = "{} skew".format(nt)
    cumulative_skew_legend: str = "Cumulative {} skew, normalized".format(nt)
    df = pd.DataFrame(
        {
            content_legend: pd.Series(content_dict),
            skew_legend: pd.Series(skew_dict),
            cumulative_skew_legend: pd.Series(skew_cumulative_dict),
        }
    )
    return df


__all__ = ["calculate_dinucleotide_skew", "skew_df", "sliding_window"]


