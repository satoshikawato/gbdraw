"""The counted bytes and window counts read a sequence without per-base arrays."""

from __future__ import annotations

import tracemalloc

import numpy as np
import pytest
from Bio.Seq import MutableSeq, Seq

from gbdraw.analysis.skew import _prefix_counts_at, counted_sequence, window_base_counts


@pytest.mark.parametrize(
    "seq",
    ["acgtuACGTUnN-", Seq("acgtuACGTUnN-"), MutableSeq("acgtuACGTUnN-"), b"acgtuACGTUnN-"],
)
def test_counted_sequence_is_one_bytes_copy_with_case_ignored_and_u_as_t(seq: object) -> None:
    assert counted_sequence(seq) == b"ACGTTACGTTNN-"


def test_counted_sequence_rejects_non_ascii_like_the_encoder_it_replaces() -> None:
    with pytest.raises(UnicodeEncodeError):
        counted_sequence("ACGTé")


def test_prefix_counts_read_unsorted_duplicate_and_out_of_range_positions() -> None:
    seq = b"ACGGTAGGCCGAT" * 7
    positions = np.array([40, 0, 91, 7, 7, 91, 13, 3, 200, -5], dtype=np.int64)
    for base in b"ACGTN":
        expected = [seq[: max(0, min(int(p), len(seq)))].count(bytes([base])) for p in positions]
        assert _prefix_counts_at(seq, base, positions).tolist() == expected


def test_prefix_counts_agree_across_chunk_boundaries() -> None:
    rng = np.random.default_rng(2)
    seq = rng.choice(np.frombuffer(b"ACGT", dtype=np.uint8), size=3_000_000).tobytes()
    edge = 1 << 20
    positions = np.array([0, 10, 5_000, edge - 1, edge, edge + 1, 2 * edge, 2_900_000, 3_000_000, 2_000_000])
    for base in b"ACGT":
        expected = [seq[:p].count(bytes([base])) for p in positions.tolist()]
        assert _prefix_counts_at(seq, base, positions).tolist() == expected


def test_window_base_counts_need_no_array_as_long_as_the_sequence() -> None:
    rng = np.random.default_rng(1)
    seq = rng.choice(np.frombuffer(b"ACGT", dtype=np.uint8), size=8_000_000).tobytes()
    starts = np.arange(0, len(seq), 50_000, dtype=np.int64)
    tracemalloc.start()
    counts = window_base_counts(seq, ord("G"), starts=starts, window=10_000)
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    assert peak < 6 << 20
    assert counts[3] == seq[150_000:160_000].count(b"G")
