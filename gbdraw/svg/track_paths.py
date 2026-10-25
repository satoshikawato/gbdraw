"""The one writer of track paths: Depth, GC content (deviation and percent) and GC skew, in both modes.

A track path writes its numbers through :func:`gbdraw.svg.path_data.path_number`, keeps a point of its
drawn line only where dropping it would move the outline by more than ``TRACK_PATH_TOLERANCE_PX``
(Ramer-Douglas-Peucker on the drawn coordinates, per drawn segment), and draws a Circular baseline as two
arcs. The data behind a track (bins, windows) and its axes are not touched here.
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
    stretch is split in the same round (level-synchronous, vectorized).
    """

    count = len(xs)
    keep = np.zeros(count, dtype=bool)
    if count <= 2:
        keep[:] = True
        return keep
    breaks = np.asarray(segment_starts, dtype=np.int64)
    breaks = breaks[(breaks > 0) & (breaks < count)]
    fixed = np.unique(np.concatenate(([0, count - 1], breaks - 1, breaks)))
    keep[fixed] = True
    starts, ends = fixed[:-1], fixed[1:]
    while True:
        lengths = ends - starts - 1
        open_ = lengths > 0
        starts, ends, lengths = starts[open_], ends[open_], lengths[open_]
        if not len(starts):
            return keep
        stretch = np.repeat(np.arange(len(starts)), lengths)
        first_offsets = np.cumsum(lengths) - lengths
        points = starts[stretch] + 1 + np.arange(int(lengths.sum())) - np.repeat(first_offsets, lengths)
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


def _arc_end(radius: float, angle: float) -> str:
    return f"{path_number(radius * np.cos(angle))} {path_number(radius * np.sin(angle))}"


def circular_baseline_back(radius: float, start_angle: float, end_angle: float) -> str:
    """Two arcs on the circle of ``radius`` from ``end_angle`` back to ``start_angle`` (radians, the drawing
    direction of a track runs from start to end; the span is at most a full turn)."""

    r = path_number(radius)
    middle = 0.5 * (start_angle + end_angle)
    return f"A{r} {r} 0 0 0 {_arc_end(radius, middle)}A{r} {r} 0 0 0 {_arc_end(radius, start_angle)}"


def circle_subpath(radius: float) -> str:
    """A closed circle of ``radius`` around the origin, from 12 o'clock clockwise, as two arcs."""

    r = path_number(radius)
    return f"M 0 {path_number(-radius)}A {r} {r} 0 0 1 0 {r}A {r} {r} 0 0 1 0 {path_number(-radius)}z"


__all__ = [
    "TRACK_PATH_TOLERANCE_PX",
    "circle_subpath",
    "circular_baseline_back",
    "kept_track_points",
    "track_line_to",
    "track_polyline",
]
