"""The counted bytes and prefix counts read a sequence without per-base arrays."""

from __future__ import annotations

import tracemalloc

import numpy as np
import pytest
from Bio.Seq import MutableSeq, Seq

from gbdraw.analysis import skew
from gbdraw.analysis.skew import _prefix_counts_at, counted_bases


@pytest.mark.parametrize(
    "seq",
    [Seq("acgtuACGTUnN-"), MutableSeq("acgtuACGTUnN-"), b"acgtuACGTUnN-"],
)
def test_counted_bases_read_case_ignored_and_u_as_t(seq: object) -> None:
    # Every byte other than A, C, G, T and U is an unknown base and reads as N (OV-2799).
    assert b"".join(counted_bases(seq)) == b"ACGTTACGTTNNN"


def test_counted_bases_are_read_a_chunk_at_a_time(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(skew, "_COUNT_CHUNK", 4)
    assert list(counted_bases(Seq("acgtuACGTUnN-"))) == [b"ACGT", b"TACG", b"TTNN", b"N"]


def test_prefix_counts_read_unsorted_duplicate_and_out_of_range_positions() -> None:
    seq = b"ACGGTAGGCCGAT" * 7
    positions = np.array([40, 0, 91, 7, 7, 91, 13, 3, 200, -5], dtype=np.int64)
    for base in b"ACGTN":
        expected = [seq[: max(0, min(int(p), len(seq)))].count(bytes([base])) for p in positions]
        assert _prefix_counts_at([seq], (base,), positions)[0].tolist() == expected


def test_prefix_counts_agree_across_chunk_boundaries() -> None:
    rng = np.random.default_rng(2)
    seq = rng.choice(np.frombuffer(b"ACGT", dtype=np.uint8), size=3_000_000).tobytes()
    edge = 1 << 20
    positions = np.array([0, 10, 5_000, edge - 1, edge, edge + 1, 2 * edge, 2_900_000, 3_000_000, 2_000_000])
    for base in b"ACGT":
        expected = [seq[:p].count(bytes([base])) for p in positions.tolist()]
        assert _prefix_counts_at([seq], (base,), positions)[0].tolist() == expected


def test_prefix_counts_of_several_bases_come_from_one_pass_and_equal_their_own_counts() -> None:
    rng = np.random.default_rng(4)
    seq = rng.choice(np.frombuffer(b"ACGTN", dtype=np.uint8), size=2_200_000).tobytes()
    positions = np.array([0, 9, 1 << 20, (1 << 20) + 1, 2_000_000, len(seq)])
    bases = tuple(b"GCN")
    counts = _prefix_counts_at([seq], bases, positions)
    for row, base in enumerate(bases):
        assert counts[row].tolist() == _prefix_counts_at([seq], (base,), positions)[0].tolist()


def test_prefix_counts_of_a_base_that_fills_a_chunk_or_is_absent_from_it() -> None:
    # A long N run fills whole chunks (counted from the other bases' offsets) and has no G (counted as none).
    rng = np.random.default_rng(5)
    edge = 1 << 20
    known = rng.choice(np.frombuffer(b"ACT", dtype=np.uint8), size=edge // 2).tobytes()
    seq = known + b"N" * (2 * edge) + b"ACGTN" * 1000 + known
    positions = np.array([0, 5, edge // 2, edge - 3, edge, edge + 7, 2 * edge + 11, 2 * edge + edge // 2 + 2_500, len(seq)])
    for base in b"NGA":
        expected = [seq[:p].count(bytes([base])) for p in positions.tolist()]
        assert _prefix_counts_at([seq], (base,), positions)[0].tolist() == expected


def test_prefix_counts_read_segments_as_one_sequence_across_chunk_boundaries() -> None:
    rng = np.random.default_rng(3)
    seq = rng.choice(np.frombuffer(b"ACGT", dtype=np.uint8), size=2_500_000).tobytes()
    edge = 1 << 20
    segments = [seq[:7], seq[7 : edge + 3], b"", seq[edge + 3 :]]
    positions = np.array([0, 6, 7, 8, edge, edge + 3, edge + 4, 2 * edge + 5, len(seq), len(seq) + 9, -1])
    for base in b"ACGT":
        expected = [seq[: max(0, p)].count(bytes([base])) for p in positions.tolist()]
        assert _prefix_counts_at(segments, (base,), positions)[0].tolist() == expected


def test_prefix_counts_need_no_array_as_long_as_the_sequence() -> None:
    rng = np.random.default_rng(1)
    seq = rng.choice(np.frombuffer(b"ACGT", dtype=np.uint8), size=8_000_000).tobytes()
    positions = np.arange(0, len(seq), 50_000, dtype=np.int64)
    tracemalloc.start()
    counts = _prefix_counts_at([seq], (ord("G"),), positions)[0]
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    assert peak < 6 << 20
    assert counts[3] == seq[:150_000].count(b"G")
