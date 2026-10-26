#!/usr/bin/env python
# coding: utf-8

"""Small spatial indexes for layout broad-phase collision filtering."""

from __future__ import annotations

import bisect
import math
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from itertools import product
from typing import Sequence, TypeAlias, cast

from .text_geometry import Point

# A uniform-grid bucket: (column, row) of ``floor(coordinate / bucket_size)``.
GridCell: TypeAlias = tuple[int, int]
# The rounding of a segment's cut points, covered on each side of every piece (px).
SEGMENT_PIECE_PAD = 1e-6


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


def _ordered_aabb(min_x: float, min_y: float, max_x: float, max_y: float) -> Aabb:
    """An ``Aabb`` from bounds that are finite floats in order already.

    It skips the constructor's checks, which the hot broad phases below (and the
    label repair's) would otherwise pay twice.
    """
    box = object.__new__(Aabb)
    object.__setattr__(box, "min_x", min_x)
    object.__setattr__(box, "min_y", min_y)
    object.__setattr__(box, "max_x", max_x)
    object.__setattr__(box, "max_y", max_y)
    return box


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
    return _ordered_aabb(min_x, min_y, max_x, max_y)


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
    pad_x = _validate_padding(pad_x) if pad_x else 0.0
    pad_y = pad_x if pad_y is None else _validate_padding(pad_y)
    return finite_aabb(min(xs) - pad_x, min(ys) - pad_y, max(xs) + pad_x, max(ys) + pad_y)


def finite_aabb(min_x: float, min_y: float, max_x: float, max_y: float) -> Aabb | None:
    """Return the box with these bounds (swapped bounds are reordered).

    ``None`` means the extent is unknown: a bound is not finite (as for
    :func:`bounding_aabb`).
    """
    if not (math.isfinite(min_x) and math.isfinite(min_y) and math.isfinite(max_x) and math.isfinite(max_y)):
        return None
    if max_x < min_x:
        min_x, max_x = max_x, min_x
    if max_y < min_y:
        min_y, max_y = max_y, min_y
    return _ordered_aabb(min_x, min_y, max_x, max_y)


def aabbs_touch(first: Aabb | None, second: Aabb | None) -> bool:
    """Whether two boxes touch; edge contact counts and a ``None`` box touches every box."""
    if first is None or second is None:
        return True
    return not (
        second.max_x < first.min_x
        or second.min_x > first.max_x
        or second.max_y < first.min_y
        or second.min_y > first.max_y
    )


def touching_aabb_pairs(bboxes: Sequence[Aabb | None]) -> Iterator[tuple[int, int]]:
    """Yield every index pair whose boxes touch; edge contact counts.

    A ``None`` box (unknown extent) touches every box. A predicate that needs two
    shapes to meet is false for every other pair, so a count over these pairs
    equals the count over all pairs. The pairs are yielded one at a time, not held.
    """
    finite = sorted(
        ((idx, bbox) for idx, bbox in enumerate(bboxes) if bbox is not None),
        key=lambda item: item[1].min_x,
    )
    unknown = [idx for idx, bbox in enumerate(bboxes) if bbox is None]
    return _touching_aabb_pairs(bboxes, finite, unknown)


def _touching_aabb_pairs(
    bboxes: Sequence[Aabb | None],
    finite: list[tuple[int, Aabb]],
    unknown: list[int],
) -> Iterator[tuple[int, int]]:
    for position, (idx, bbox) in enumerate(finite):
        for later in range(position + 1, len(finite)):
            jdx, other = finite[later]
            if other.min_x > bbox.max_x:
                break
            if other.min_y > bbox.max_y or other.max_y < bbox.min_y:
                continue
            yield idx, jdx
    for idx in unknown:
        for jdx in range(len(bboxes)):
            if jdx != idx and (bboxes[jdx] is not None or jdx > idx):
                yield idx, jdx


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
    return [idx for idx in indexes if idx != skip and aabbs_touch(query, bboxes[idx])]


def grid_bucket_size(bboxes: Sequence[Aabb], padding: float = 0.0) -> float:
    """A grid bucket for these boxes: the largest box side plus ``padding``, within 16 to 256 (16 + padding for none)."""
    max_extent = max((max(bbox.max_x - bbox.min_x, bbox.max_y - bbox.min_y) for bbox in bboxes), default=16.0)
    return max(16.0, min(256.0, max_extent + float(padding)))


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


def aabb_cells(bbox: Aabb, bucket_size: float, padding: float = 0.0) -> Iterator[GridCell]:
    """The grid buckets that ``bbox`` grown by ``padding`` meets."""
    return product(
        _bucket_range(bbox.min_x - padding, bbox.max_x + padding, bucket_size),
        _bucket_range(bbox.min_y - padding, bbox.max_y + padding, bucket_size),
    )


def segment_cells(start: Point, end: Point, bucket_size: float) -> set[GridCell]:
    """The grid buckets that cover the segment from ``start`` to ``end``.

    The segment is cut into pieces no longer than one bucket, and each piece adds
    the buckets its box meets, padded by :data:`SEGMENT_PIECE_PAD` for the rounding
    of the cut points. Every point of the segment lies in a returned bucket, so two
    segments that meet share one; a long diagonal covers about three buckets per
    bucket length instead of its whole box.
    """
    ax, ay = start
    dx, dy = end[0] - ax, end[1] - ay
    length = math.hypot(dx, dy)
    if not math.isfinite(length):
        raise ValueError("segment points must be finite")
    pieces = max(1, math.ceil(length / bucket_size))
    cells: set[GridCell] = set()
    px, py = ax, ay
    for step in range(1, pieces + 1):
        if step == pieces:
            qx, qy = end
        else:
            t = step / pieces
            qx, qy = ax + dx * t, ay + dy * t
        rows = _bucket_range(min(py, qy) - SEGMENT_PIECE_PAD, max(py, qy) + SEGMENT_PIECE_PAD, bucket_size)
        for column in _bucket_range(min(px, qx) - SEGMENT_PIECE_PAD, max(px, qx) + SEGMENT_PIECE_PAD, bucket_size):
            for row in rows:
                cells.add((column, row))
        px, py = qx, qy
    return cells


def cell_partners(
    inserted: Iterable[Iterable[GridCell]],
    queried: Iterable[Iterable[GridCell]],
) -> Iterator[set[int]]:
    """Yield, for each item in order, the later items that may touch it.

    Item ``j`` is a partner of item ``i < j`` when a bucket of ``queried[i]`` is
    one of ``inserted[j]``. Each item's buckets are a cover of its shape (its box
    as :func:`aabb_cells` gives it, its segments as :func:`segment_cells` gives
    them); the grid holds every inserted bucket, but only one item's partners at a
    time, so the pairs of a crowded set are never all in memory at once.
    """
    grid: dict[GridCell, list[int]] = {}
    for item_idx, cells in enumerate(inserted):
        for cell in cells:
            grid.setdefault(cell, []).append(item_idx)
    empty: list[int] = []
    for item_idx, cells in enumerate(queried):
        partners: set[int] = set()
        for cell in cells:
            bucket = grid.get(cell, empty)
            # Bucket lists hold indexes in increasing order: the later items are at the end.
            partners.update(bucket[bisect.bisect_right(bucket, item_idx):])
        yield partners


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
        self._buckets: dict[GridCell, list[object]] = {}
        self._order: dict[object, int] = {}

    def insert(self, item_id: object, bbox: Aabb) -> None:
        if item_id not in self._order:
            self._order[item_id] = len(self._order)
        for cell in aabb_cells(bbox, self.bucket_size):
            self._buckets.setdefault(cell, []).append(item_id)

    def remove(self, item_id: object, bbox: Aabb) -> None:
        """Take out ``item_id``, inserted with ``bbox``; it keeps its place in the query order."""
        for cell in aabb_cells(bbox, self.bucket_size):
            self._buckets[cell].remove(item_id)

    def query(self, bbox: Aabb, padding: float = 0.0) -> list[object]:
        """The ids in the buckets the box grown by ``padding`` meets, in first-insertion order."""
        return self._items_in(aabb_cells(bbox, self.bucket_size, _validate_padding(padding)))

    def query_segment(self, start: Point, end: Point) -> list[object]:
        """The ids in the buckets that cover the segment (:func:`segment_cells`), in first-insertion order."""
        return self._items_in(segment_cells(start, end, self.bucket_size))

    def _items_in(self, cells: Iterable[GridCell]) -> list[object]:
        empty: list[object] = []
        candidates = list(dict.fromkeys(item for cell in cells for item in self._buckets.get(cell, empty)))
        candidates.sort(key=self._order.__getitem__)
        return candidates


class TouchingAabbIndex:
    """Boxes by index in a grid, kept as boxes change.

    ``touching(query, skip=...)`` equals ``touching_aabb_indices(boxes, query, skip=...)``
    but looks only at the boxes in the query's buckets. A ``None`` box (unknown
    extent) touches every box.
    """

    def __init__(self, bucket_size: float) -> None:
        self.boxes: list[Aabb | None] = []
        self._unknown: set[int] = set()
        self._index = AabbIndex(bucket_size=bucket_size)
        # The grid answers in first-insertion order. That is index order until an unknown box
        # becomes finite: its index then enters the grid after later ones.
        self._grid_in_index_order = True

    @classmethod
    def of(cls, bboxes: Sequence[Aabb | None]) -> TouchingAabbIndex:
        """The boxes in order, in a grid sized by :func:`grid_bucket_size`."""
        index = cls(grid_bucket_size([bbox for bbox in bboxes if bbox is not None]))
        for bbox in bboxes:
            index.append(bbox)
        return index

    def append(self, bbox: Aabb | None) -> None:
        """Add ``bbox`` as the next index."""
        self.boxes.append(None)
        self._place(len(self.boxes) - 1, bbox)

    def replace(self, idx: int, bbox: Aabb | None) -> None:
        """Box ``idx`` is now ``bbox``."""
        old = self.boxes[idx]
        if old is None:
            self._unknown.discard(idx)
            if bbox is not None:
                self._grid_in_index_order = False
        else:
            self._index.remove(idx, old)
        self._place(idx, bbox)

    def touching(self, query: Aabb | None, *, skip: int | None = None) -> list[int]:
        """The indexes of the boxes that touch ``query``, in index order; ``skip`` is left out."""
        if query is None:
            return [idx for idx in range(len(self.boxes)) if idx != skip]
        boxes = self.boxes
        size = self._index.bucket_size
        buckets = (math.floor(query.max_x / size) - math.floor(query.min_x / size) + 1) * (
            math.floor(query.max_y / size) - math.floor(query.min_y / size) + 1
        )
        if buckets >= len(boxes):
            # A query that spans more buckets than there are boxes: the scan is cheaper.
            return touching_aabb_indices(boxes, query, skip=skip)
        candidates = cast("list[int]", self._index.query(query))
        found = [idx for idx in candidates if idx != skip and aabbs_touch(query, boxes[idx])]
        if self._unknown:
            found.extend(idx for idx in self._unknown if idx != skip)
            found.sort()
        elif not self._grid_in_index_order:
            found.sort()
        return found

    def _place(self, idx: int, bbox: Aabb | None) -> None:
        self.boxes[idx] = bbox
        if bbox is None:
            self._unknown.add(idx)
        else:
            self._index.insert(idx, bbox)


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


def candidate_aabb_partners(
    bboxes: Sequence[Aabb],
    padding: float = 0.0,
    bucket_size: float = 128.0,
) -> Iterator[set[int]]:
    """Yield, for each box in order, the later boxes it may touch (:func:`cell_partners` over box cells).

    Box ``j`` is a partner of box ``i < j`` when a grid bucket holds box ``i`` and meets
    box ``j`` grown by ``padding``.
    """
    padding = _validate_padding(padding)
    bucket_size = _validate_bucket_size(bucket_size)
    return cell_partners(
        (aabb_cells(bbox, bucket_size, padding) for bbox in bboxes),
        (aabb_cells(bbox, bucket_size) for bbox in bboxes),
    )


__all__ = [
    "Aabb",
    "AabbIndex",
    "Interval",
    "IntervalIndex",
    "GridCell",
    "aabb_cells",
    "aabbs_touch",
    "bounding_aabb",
    "candidate_aabb_partners",
    "cell_partners",
    "finite_aabb",
    "grid_bucket_size",
    "segment_cells",
    "split_circular_interval",
    "TouchingAabbIndex",
    "touching_aabb_indices",
    "touching_aabb_pairs",
    "union_aabbs",
]
