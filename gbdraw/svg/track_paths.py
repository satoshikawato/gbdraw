"""The one writer of track paths: Depth, GC content (deviation and percent) and GC skew, in both modes.

A track path writes its numbers through :func:`gbdraw.svg.path_data.path_number`, keeps a point of its
drawn line only where dropping it would move the outline by more than ``TRACK_PATH_TOLERANCE_PX``
(Ramer-Douglas-Peucker on the drawn coordinates, per drawn segment). A Circular area closes on its
baseline through :mod:`gbdraw.svg.circular_arcs`. The data behind a track (bins, windows) and its axes are
not touched here.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .path_data import path_line_to, path_number

TRACK_PATH_TOLERANCE_PX = 0.02


def kept_track_points(xs: NDArray[np.float64], ys: NDArray[np.float64], segment_starts: Sequence[int] = ()) -> NDArray[np.bool_]:
    """Which points of a drawn line the path keeps.

    The ends of the line and of each drawn segment (``segment_starts`` are the indices where a segment after
    the first begins) are always kept. Between two kept points, the point farthest from the straight
    segment joining them is kept while that distance exceeds ``TRACK_PATH_TOLERANCE_PX``; every pending
    stretch is split in the same round (level-synchronous, vectorized). A NaN or infinite coordinate raises
    ``ValueError``.
    """

    if not (np.isfinite(xs).all() and np.isfinite(ys).all()):
        raise ValueError("A track path needs finite coordinates; a sample is NaN or infinite.")
    count = len(xs)
    keep = np.zeros(count, dtype=bool)
    if count <= 2:
        keep[:] = True
        return keep
    # Every index array is np.intp: Pyodide's NumPy is 32-bit and ufunc.reduceat refuses int64 indices there.
    breaks = np.asarray(segment_starts, dtype=np.intp)
    breaks = breaks[(breaks > 0) & (breaks < count)]
    fixed = np.unique(np.concatenate((np.array([0, count - 1], dtype=np.intp), breaks - 1, breaks)))
    keep[fixed] = True
    starts, ends = fixed[:-1], fixed[1:]
    while True:
        lengths = ends - starts - 1
        open_ = lengths > 0
        starts, ends, lengths = starts[open_], ends[open_], lengths[open_]
        if not len(starts):
            return keep
        stretch = np.repeat(np.arange(len(starts), dtype=np.intp), lengths)
        first_offsets = np.cumsum(lengths, dtype=np.intp) - lengths
        points = starts[stretch] + 1 + np.arange(int(lengths.sum()), dtype=np.intp) - np.repeat(first_offsets, lengths)
        ax, ay = xs[starts][stretch], ys[starts][stretch]
        dx, dy = xs[ends][stretch] - ax, ys[ends][stretch] - ay
        px, py = xs[points] - ax, ys[points] - ay
        squared = dx * dx + dy * dy
        along = np.clip(np.divide(px * dx + py * dy, squared, out=np.zeros_like(squared), where=squared > 0), 0.0, 1.0)
        distance = np.hypot(px - along * dx, py - along * dy)
        farthest = np.maximum.reduceat(distance, first_offsets)
        at_farthest = np.flatnonzero(distance == farthest[stretch])
        first_of_stretch = np.r_[True, stretch[at_farthest][1:] != stretch[at_farthest][:-1]]
        split_points = points[at_farthest[first_of_stretch]]
        split = farthest > TRACK_PATH_TOLERANCE_PX
        split_points = split_points[split]
        keep[split_points] = True
        starts = np.concatenate((starts[split], split_points))
        ends = np.concatenate((split_points, ends[split]))


def _kept(xs: ArrayLike, ys: ArrayLike, segment_starts: Sequence[int]) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    x_values, y_values = np.asarray(xs, dtype=float), np.asarray(ys, dtype=float)
    keep = kept_track_points(x_values, y_values, segment_starts)
    return x_values[keep], y_values[keep]


def track_line_to(xs: ArrayLike, ys: ArrayLike, segment_starts: Sequence[int] = ()) -> str:
    """``L x y`` for each kept point of a drawn line (its first point included)."""

    return path_line_to(*_kept(xs, ys, segment_starts))


def track_polyline(xs: ArrayLike, ys: ArrayLike, segment_starts: Sequence[int] = ()) -> str:
    """``M x y`` at the first point of a drawn line, then ``L x y`` for each other kept point."""

    x_values, y_values = _kept(xs, ys, segment_starts)
    return f"M{path_number(x_values[0])} {path_number(y_values[0])}{path_line_to(x_values[1:], y_values[1:])}"


__all__ = [
    "TRACK_PATH_TOLERANCE_PX",
    "kept_track_points",
    "track_line_to",
    "track_polyline",
]
