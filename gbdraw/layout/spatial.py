#!/usr/bin/env python
# coding: utf-8

"""Small spatial indexes for layout broad-phase collision filtering."""

from __future__ import annotations

import math
from collections.abc import Iterable
from dataclasses import dataclass
from typing import Sequence, cast


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

    @property
    def width(self) -> float:
        """Horizontal extent of the box."""
        return self.max_x - self.min_x

    @property
    def height(self) -> float:
        """Vertical extent of the box."""
        return self.max_y - self.min_y

    def translated(self, dx: float, dy: float) -> Aabb:
        """Return a copy translated by ``dx`` and ``dy``."""
        dx = _validate_finite_value(dx, name="bbox translation")
        dy = _validate_finite_value(dy, name="bbox translation")
        return Aabb(
            self.min_x + dx,
            self.min_y + dy,
            self.max_x + dx,
            self.max_y + dy,
        )

    def expanded(self, x: float, y: float | None = None) -> Aabb:
        """Return a copy expanded equally on both sides of each axis."""
        x = _validate_padding(x)
        y = x if y is None else _validate_padding(y)
        return Aabb(
            self.min_x - x,
            self.min_y - y,
            self.max_x + x,
            self.max_y + y,
        )

    def intersects(self, other: Aabb, clearance: float = 0.0) -> bool:
        """Return whether the boxes overlap with less than ``clearance`` between them.

        Edge contact is not an intersection. Consequently, two boxes separated
        by exactly ``clearance`` satisfy that clearance.
        """
        if not isinstance(other, Aabb):
            raise TypeError("other must be an Aabb")
        clearance = _validate_padding(clearance)
        return (
            self.min_x < other.max_x + clearance
            and self.max_x > other.min_x - clearance
            and self.min_y < other.max_y + clearance
            and self.max_y > other.min_y - clearance
        )


def _validate_finite_value(value: float, *, name: str) -> float:
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


def union_aabbs(bboxes: Iterable[Aabb]) -> Aabb | None:
    """Return the union of ``bboxes``, or ``None`` for an empty iterable."""
    iterator = iter(bboxes)
    try:
        first = next(iterator)
    except StopIteration:
        return None
    if not isinstance(first, Aabb):
        raise TypeError("all bboxes must be Aabb instances")

    min_x = first.min_x
    min_y = first.min_y
    max_x = first.max_x
    max_y = first.max_y
    for bbox in iterator:
        if not isinstance(bbox, Aabb):
            raise TypeError("all bboxes must be Aabb instances")
        min_x = min(min_x, bbox.min_x)
        min_y = min(min_y, bbox.min_y)
        max_x = max(max_x, bbox.max_x)
        max_y = max(max_y, bbox.max_y)
    return Aabb(min_x, min_y, max_x, max_y)


def bounding_aabb(
    points: Iterable[tuple[float, float]],
    pad_x: float = 0.0,
    pad_y: float | None = None,
) -> Aabb | None:
    """Return the box around ``points`` grown by ``pad_x`` / ``pad_y`` on each side.

    ``None`` means the extent is unknown: there are no points or a coordinate is
    not finite. The touch queries below treat such a box as touching every box.
    """
    xs: list[float] = []
    ys: list[float] = []
    for x, y in points:
        xs.append(float(x))
        ys.append(float(y))
    if not xs or not all(map(math.isfinite, xs)) or not all(map(math.isfinite, ys)):
        return None
    pad_x = _validate_padding(pad_x)
    pad_y = pad_x if pad_y is None else _validate_padding(pad_y)
    return Aabb(min(xs) - pad_x, min(ys) - pad_y, max(xs) + pad_x, max(ys) + pad_y)


def touching_aabb_pairs(bboxes: Sequence[Aabb | None]) -> list[tuple[int, int]]:
    """Return every index pair whose boxes touch; edge contact counts.

    A ``None`` box (unknown extent) touches every box. A predicate that needs two
    shapes to meet is false for every other pair, so a count over these pairs
    equals the count over all pairs.
    """
    finite = sorted(
        ((idx, bbox) for idx, bbox in enumerate(bboxes) if bbox is not None),
        key=lambda item: item[1].min_x,
    )
    pairs: list[tuple[int, int]] = []
    for position, (idx, bbox) in enumerate(finite):
        for jdx, other in finite[position + 1:]:
            if other.min_x > bbox.max_x:
                break
            if other.min_y > bbox.max_y or other.max_y < bbox.min_y:
                continue
            pairs.append((idx, jdx))
    unknown = [idx for idx, bbox in enumerate(bboxes) if bbox is None]
    for idx in unknown:
        for jdx in range(len(bboxes)):
            if jdx != idx and (bboxes[jdx] is not None or jdx > idx):
                pairs.append((idx, jdx))
    return pairs


def touching_aabb_indices(
    bboxes: Sequence[Aabb | None],
    query: Aabb | None,
    *,
    skip: int | None = None,
    candidates: Iterable[int] | None = None,
) -> list[int]:
    """Return the indexes of ``bboxes`` that touch ``query``; edge contact counts.

    ``candidates`` limits and orders the indexes looked at (default: all, in
    order); ``skip`` is left out. A ``None`` box touches every box.
    """
    indexes = range(len(bboxes)) if candidates is None else candidates
    if query is None:
        return [idx for idx in indexes if idx != skip]
    min_x, min_y, max_x, max_y = query.min_x, query.min_y, query.max_x, query.max_y
    touching: list[int] = []
    for idx in indexes:
        if idx == skip:
            continue
        other = bboxes[idx]
        if other is None or not (
            other.max_x < min_x or other.min_x > max_x or other.max_y < min_y or other.min_y > max_y
        ):
            touching.append(idx)
    return touching


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
            left = int(cast(int, candidate_idx))  # ids inserted above are item indexes
            right = int(item_idx)
            pairs.add((left, right) if left < right else (right, left))
        index.insert(item_idx, bbox)
    return pairs


__all__ = [
    "Aabb",
    "AabbIndex",
    "Interval",
    "IntervalIndex",
    "bounding_aabb",
    "candidate_aabb_pairs",
    "split_circular_interval",
    "touching_aabb_indices",
    "touching_aabb_pairs",
    "union_aabbs",
]
