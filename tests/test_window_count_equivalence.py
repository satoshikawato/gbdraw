"""The vectorized window counts must equal the former per-base Python loops."""

from __future__ import annotations

import random

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gbdraw.analysis.gc import circular_dinucleotide_content_df
from gbdraw.analysis.skew import counted_dinucleotide, skew_df


def _old_prefix(seq_bytes: bytes, target_base: int) -> list[int]:
    prefix = [0] * (len(seq_bytes) + 1)
    running_count = 0
    for idx, value in enumerate(seq_bytes, start=1):
        if value == target_base:
            running_count += 1
        prefix[idx] = running_count
    return prefix


def _old_window_count(prefix, seq_length, *, start, window, total_count) -> int:
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


def _old_skew_df(record, window, step, nt) -> pd.DataFrame:
    nt, nt_1, nt_2, seq_str = counted_dinucleotide(record, nt)
    seq_length = len(seq_str)
    content_legend, skew_legend = f"{nt} content", f"{nt} skew"
    cumulative_legend = f"Cumulative {nt} skew, normalized"
    if seq_length == 0 or window <= 0 or step <= 0:
        return pd.DataFrame(columns=[content_legend, skew_legend, cumulative_legend])
    seq_bytes = seq_str.encode("ascii")
    p1, p2 = _old_prefix(seq_bytes, ord(nt_1)), _old_prefix(seq_bytes, ord(nt_2))
    starts, content, skews, cumulative = [], [], [], []
    skew_sum = 0.0
    for start in range(0, seq_length, step):
        c1 = _old_window_count(p1, seq_length, start=start, window=window, total_count=p1[-1])
        c2 = _old_window_count(p2, seq_length, start=start, window=window, total_count=p2[-1])
        total = c1 + c2
        skew = 0.0 if total == 0 else (c1 - c2) / total
        starts.append(start)
        content.append(total / float(window))
        skews.append(skew)
        skew_sum += skew
        cumulative.append(skew_sum)
    max_skew = max((abs(v) for v in skews), default=0.0)
    max_cumulative = max((abs(v) for v in cumulative), default=0.0)
    factor = 0.0 if max_cumulative == 0 else max_skew / max_cumulative
    return pd.DataFrame(
        {
            content_legend: content,
            skew_legend: skews,
            cumulative_legend: [v * factor for v in cumulative],
        },
        index=starts,
    )


def _old_circular_content_df(record, window, step, nt) -> pd.DataFrame:
    nt, nt_1, nt_2, seq_str = counted_dinucleotide(record, nt)
    seq_length = len(seq_str)
    legend = f"{nt} content"
    if seq_length == 0 or window <= 0 or step <= 0:
        return pd.DataFrame(columns=[legend])
    seq_bytes = seq_str.encode("ascii")
    p1, p2 = _old_prefix(seq_bytes, ord(nt_1)), _old_prefix(seq_bytes, ord(nt_2))
    half_window = int(window) // 2
    starts, values = [], []
    for position in range(0, seq_length, step):
        circular_start = (position - half_window) % seq_length
        c1 = _old_window_count(p1, seq_length, start=circular_start, window=window, total_count=p1[-1])
        c2 = _old_window_count(p2, seq_length, start=circular_start, window=window, total_count=p2[-1])
        starts.append(position)
        values.append((c1 + c2) / float(window))
    return pd.DataFrame({legend: values}, index=starts)


def _record(sequence: str) -> SeqRecord:
    return SeqRecord(Seq(sequence), id="rec", annotations={"molecule_type": "DNA"})


def _random_sequence(rng: random.Random, length: int, alphabet: str) -> str:
    return "".join(rng.choice(alphabet) for _ in range(length))


_CASES = [
    # (length, window, step): windows shorter than, equal to and longer than the record,
    # exact multiples (full cycles with no remainder) and steps that skip the end.
    (1, 1, 1), (1, 5, 1), (7, 3, 1), (7, 7, 2), (7, 14, 3), (7, 15, 7), (60, 20, 7),
    (60, 60, 60), (60, 61, 5), (60, 500, 11), (101, 33, 100), (257, 50, 10),
]


@pytest.mark.parametrize("length,window,step", _CASES)
@pytest.mark.parametrize("alphabet", ["ACGT", "ACGTNRYKMacgtn-", "GGCCa"])
@pytest.mark.parametrize("nt", ["GC", "AT", "AU", "CG", "TA"])
def test_vectorized_window_counts_match_the_python_loops(
    length: int, window: int, step: int, alphabet: str, nt: str
) -> None:
    rng = random.Random(f"{length}-{window}-{step}-{alphabet}-{nt}")
    record = _record(_random_sequence(rng, length, alphabet))

    pd.testing.assert_frame_equal(
        skew_df(record, window, step, nt),
        _old_skew_df(record, window, step, nt),
        check_exact=True,
    )
    pd.testing.assert_frame_equal(
        circular_dinucleotide_content_df(record, window, step, nt),
        _old_circular_content_df(record, window, step, nt),
        check_exact=True,
    )


def test_vectorized_window_counts_match_on_a_long_sequence() -> None:
    rng = random.Random(7)
    record = _record(_random_sequence(rng, 200_000, "ACGTACGTN"))

    pd.testing.assert_frame_equal(
        skew_df(record, 10_000, 500, "GC"),
        _old_skew_df(record, 10_000, 500, "GC"),
        check_exact=True,
    )
    pd.testing.assert_frame_equal(
        circular_dinucleotide_content_df(record, 10_000, 500, "AT"),
        _old_circular_content_df(record, 10_000, 500, "AT"),
        check_exact=True,
    )


@pytest.mark.parametrize("window,step", [(0, 1), (5, 0), (-1, 3)])
def test_degenerate_windows_and_empty_records_return_empty_frames(
    window: int, step: int
) -> None:
    for record in (_record("ACGT"), _record("")):
        pd.testing.assert_frame_equal(
            skew_df(record, window, step, "GC"),
            _old_skew_df(record, window, step, "GC"),
        )
        pd.testing.assert_frame_equal(
            circular_dinucleotide_content_df(record, window, step, "GC"),
            _old_circular_content_df(record, window, step, "GC"),
        )

