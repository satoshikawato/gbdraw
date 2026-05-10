from __future__ import annotations

import pytest

from gbdraw.layout.spatial import (
    Aabb,
    AabbIndex,
    Interval,
    IntervalIndex,
    candidate_aabb_pairs,
    split_circular_interval,
)


def _intervals_overlap(left: Interval, right: Interval, padding: float = 0.0) -> bool:
    padded_left = Interval(left.start - padding, left.end + padding)
    return not (padded_left.end < right.start or right.end < padded_left.start)


def _bboxes_overlap(left: Aabb, right: Aabb, padding: float = 0.0) -> bool:
    return not (
        left.max_x + padding < right.min_x
        or right.max_x + padding < left.min_x
        or left.max_y + padding < right.min_y
        or right.max_y + padding < left.min_y
    )


def test_interval_and_aabb_normalize_reversed_bounds() -> None:
    assert Interval(10, 2) == Interval(2, 10)
    assert Aabb(10, 20, 1, 3) == Aabb(1, 3, 10, 20)


@pytest.mark.parametrize("bad_bucket_size", [0, -1, float("inf")])
def test_indexes_require_positive_finite_bucket_size(bad_bucket_size: float) -> None:
    with pytest.raises(ValueError):
        IntervalIndex(bucket_size=bad_bucket_size)
    with pytest.raises(ValueError):
        AabbIndex(bucket_size=bad_bucket_size)


def test_interval_index_query_padding_and_duplicate_suppression() -> None:
    index = IntervalIndex(bucket_size=10)
    index.insert("wide", Interval(0, 25))
    index.insert("near", Interval(31, 33))
    index.insert("far", Interval(80, 90))

    assert index.query(Interval(26, 28)) == ["wide"]
    assert index.query(Interval(28, 29), padding=3) == ["wide", "near"]


def test_aabb_index_query_padding_and_duplicate_suppression() -> None:
    index = AabbIndex(bucket_size=10)
    index.insert("wide", Aabb(0, 0, 25, 25))
    index.insert("near", Aabb(31, 0, 35, 10))
    index.insert("far", Aabb(80, 80, 90, 90))

    assert index.query(Aabb(26, 3, 28, 8)) == ["wide"]
    assert index.query(Aabb(28, 3, 29, 8), padding=3) == ["wide", "near"]


def test_split_circular_interval_splits_origin_spanning_ranges() -> None:
    assert split_circular_interval(10, 20, 100) == [Interval(10, 20)]
    assert split_circular_interval(90, 12, 100) == [Interval(90, 100), Interval(1, 12)]


def test_interval_index_candidates_do_not_miss_overlapping_intervals() -> None:
    intervals = [
        Interval(0, 10),
        Interval(8, 12),
        Interval(30, 40),
        Interval(39, 45),
        Interval(-5, -1),
    ]
    index = IntervalIndex(bucket_size=7)
    for item_idx, interval in enumerate(intervals):
        index.insert(item_idx, interval)

    for query_idx, query_interval in enumerate(intervals):
        candidates = set(index.query(query_interval))
        expected = {
            item_idx
            for item_idx, interval in enumerate(intervals)
            if item_idx != query_idx and _intervals_overlap(query_interval, interval)
        }
        assert expected <= candidates - {query_idx}


def test_aabb_candidate_pairs_do_not_miss_overlapping_bboxes() -> None:
    bboxes = [
        Aabb(0, 0, 10, 10),
        Aabb(8, 2, 12, 6),
        Aabb(30, 30, 35, 35),
        Aabb(34, 34, 40, 40),
        Aabb(-5, -5, -1, -1),
    ]

    candidates = candidate_aabb_pairs(bboxes, bucket_size=7)
    expected = {
        (idx, jdx)
        for idx in range(len(bboxes))
        for jdx in range(idx + 1, len(bboxes))
        if _bboxes_overlap(bboxes[idx], bboxes[jdx])
    }

    assert expected <= candidates


def test_candidate_aabb_pairs_respects_padding() -> None:
    bboxes = [
        Aabb(0, 0, 10, 10),
        Aabb(12, 0, 20, 10),
        Aabb(30, 0, 40, 10),
    ]

    assert (0, 1) not in candidate_aabb_pairs(bboxes, bucket_size=1)
    assert (0, 1) in candidate_aabb_pairs(bboxes, padding=2, bucket_size=1)
    assert (1, 2) not in candidate_aabb_pairs(bboxes, padding=2, bucket_size=1)
