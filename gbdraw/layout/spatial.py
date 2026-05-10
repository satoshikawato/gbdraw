#!/usr/bin/env python
# coding: utf-8

"""Small spatial indexes for layout broad-phase collision filtering."""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Sequence


@dataclass(frozen=True)
class Interval:
    """Normalized one-dimensional interval."""

    start: float
    end: float

    def __post_init__(self) -> None:
        start = float(self.start)
        end = float(self.end)
        if not math.isfinite(start) or not math.isfinite(end):
            raise ValueError("interval bounds must be finite")
        if end < start:
            start, end = end, start
        object.__setattr__(self, "start", start)
        object.__setattr__(self, "end", end)


@dataclass(frozen=True)
class Aabb:
    """Normalized axis-aligned bounding box."""

    min_x: float
    min_y: float
    max_x: float
    max_y: float

    def __post_init__(self) -> None:
        min_x = float(self.min_x)
        min_y = float(self.min_y)
        max_x = float(self.max_x)
        max_y = float(self.max_y)
        if not all(math.isfinite(value) for value in (min_x, min_y, max_x, max_y)):
            raise ValueError("bbox bounds must be finite")
        if max_x < min_x:
            min_x, max_x = max_x, min_x
        if max_y < min_y:
            min_y, max_y = max_y, min_y
        object.__setattr__(self, "min_x", min_x)
        object.__setattr__(self, "min_y", min_y)
        object.__setattr__(self, "max_x", max_x)
        object.__setattr__(self, "max_y", max_y)


def _validate_bucket_size(bucket_size: float) -> float:
    bucket_size = float(bucket_size)
    if not math.isfinite(bucket_size) or bucket_size <= 0.0:
        raise ValueError("bucket_size must be a positive finite value")
    return bucket_size


def _validate_padding(padding: float) -> float:
    padding = float(padding)
    if not math.isfinite(padding) or padding < 0.0:
        raise ValueError("padding must be a non-negative finite value")
    return padding


def _bucket_range(start: float, end: float, bucket_size: float) -> range:
    start_bucket = math.floor(start / bucket_size)
    end_bucket = math.floor(end / bucket_size)
    return range(start_bucket, end_bucket + 1)


class IntervalIndex:
    """Uniform-bucket interval index returning candidate item ids."""

    def __init__(self, bucket_size: float = 100.0) -> None:
        self.bucket_size = _validate_bucket_size(bucket_size)
        self._buckets: dict[int, list[object]] = {}
        self._order: dict[object, int] = {}

    def insert(self, item_id: object, interval: Interval) -> None:
        if item_id not in self._order:
            self._order[item_id] = len(self._order)
        for bucket in _bucket_range(interval.start, interval.end, self.bucket_size):
            self._buckets.setdefault(bucket, []).append(item_id)

    def query(self, interval: Interval, padding: float = 0.0) -> list[object]:
        padding = _validate_padding(padding)
        query_interval = Interval(interval.start - padding, interval.end + padding)
        seen: set[object] = set()
        candidates: list[object] = []
        for bucket in _bucket_range(query_interval.start, query_interval.end, self.bucket_size):
            for item_id in self._buckets.get(bucket, []):
                if item_id in seen:
                    continue
                seen.add(item_id)
                candidates.append(item_id)
        candidates.sort(key=self._order.__getitem__)
        return candidates


class AabbIndex:
    """Uniform-grid AABB index returning candidate item ids."""

    def __init__(self, bucket_size: float = 128.0) -> None:
        self.bucket_size = _validate_bucket_size(bucket_size)
        self._buckets: dict[tuple[int, int], list[object]] = {}
        self._order: dict[object, int] = {}

    def insert(self, item_id: object, bbox: Aabb) -> None:
        if item_id not in self._order:
            self._order[item_id] = len(self._order)
        for bucket_x in _bucket_range(bbox.min_x, bbox.max_x, self.bucket_size):
            for bucket_y in _bucket_range(bbox.min_y, bbox.max_y, self.bucket_size):
                self._buckets.setdefault((bucket_x, bucket_y), []).append(item_id)

    def query(self, bbox: Aabb, padding: float = 0.0) -> list[object]:
        padding = _validate_padding(padding)
        query_bbox = Aabb(
            bbox.min_x - padding,
            bbox.min_y - padding,
            bbox.max_x + padding,
            bbox.max_y + padding,
        )
        seen: set[object] = set()
        candidates: list[object] = []
        for bucket_x in _bucket_range(query_bbox.min_x, query_bbox.max_x, self.bucket_size):
            for bucket_y in _bucket_range(query_bbox.min_y, query_bbox.max_y, self.bucket_size):
                for item_id in self._buckets.get((bucket_x, bucket_y), []):
                    if item_id in seen:
                        continue
                    seen.add(item_id)
                    candidates.append(item_id)
        candidates.sort(key=self._order.__getitem__)
        return candidates


def split_circular_interval(start: int, end: int, genome_length: int | None) -> list[Interval]:
    """Split a circular interval into linear pieces for broad-phase indexing."""
    start = int(start)
    end = int(end)
    if start <= end:
        return [Interval(start, end)]
    if genome_length is None or int(genome_length) <= 0:
        return [Interval(start, end)]
    length = int(genome_length)
    return [Interval(start, length), Interval(1, end)]


def candidate_aabb_pairs(
    bboxes: Sequence[Aabb],
    padding: float = 0.0,
    bucket_size: float = 128.0,
) -> set[tuple[int, int]]:
    """Return candidate bbox index pairs using an AABB broad phase."""
    pairs: set[tuple[int, int]] = set()
    index = AabbIndex(bucket_size=bucket_size)
    for item_idx, bbox in enumerate(bboxes):
        for candidate_idx in index.query(bbox, padding=padding):
            left = int(candidate_idx)
            right = int(item_idx)
            pairs.add((left, right) if left < right else (right, left))
        index.insert(item_idx, bbox)
    return pairs


__all__ = [
    "Aabb",
    "AabbIndex",
    "Interval",
    "IntervalIndex",
    "candidate_aabb_pairs",
    "split_circular_interval",
]
