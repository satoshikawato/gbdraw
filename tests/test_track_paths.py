"""The one track-path writer: RDP at 0.02 SVG px per drawn segment, Circular baselines as two arcs
(incremental-render OV-1230 D-144..D-147)."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest

from gbdraw.svg.circular_tracks import generate_circle_path_desc, generate_circular_scalar_area_path_desc
from gbdraw.svg.track_paths import (
    TRACK_PATH_TOLERANCE_PX,
    circular_baseline_back,
    kept_track_points,
    track_line_to,
    track_polyline,
)


def _reference_keep(xs: np.ndarray, ys: np.ndarray, fixed: list[int]) -> np.ndarray:
    """Recursive Ramer-Douglas-Peucker on each stretch between fixed points, distance to the segment."""

    keep = np.zeros(len(xs), dtype=bool)
    keep[fixed] = True
    stack = list(zip(fixed[:-1], fixed[1:]))
    while stack:
        a, b = stack.pop()
        if b - a < 2:
            continue
        dx, dy = xs[b] - xs[a], ys[b] - ys[a]
        px, py = xs[a + 1:b] - xs[a], ys[a + 1:b] - ys[a]
        squared = dx * dx + dy * dy
        along = np.clip((px * dx + py * dy) / squared, 0, 1) if squared else np.zeros_like(px)
        distance = np.hypot(px - along * dx, py - along * dy)
        index = int(distance.argmax())
        if distance[index] > TRACK_PATH_TOLERANCE_PX:
            k = a + 1 + index
            keep[k] = True
            stack += [(a, k), (k, b)]
    return keep


def _distance_to_polyline(px: float, py: float, xs: np.ndarray, ys: np.ndarray) -> float:
    ax, ay, bx, by = xs[:-1], ys[:-1], xs[1:], ys[1:]
    dx, dy = bx - ax, by - ay
    squared = dx * dx + dy * dy
    along = np.clip(np.divide((px - ax) * dx + (py - ay) * dy, squared, out=np.zeros_like(squared), where=squared > 0), 0, 1)
    return float(np.min(np.hypot(px - ax - along * dx, py - ay - along * dy)))


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_kept_points_match_the_recursive_reference_and_stay_within_tolerance(seed: int) -> None:
    rng = np.random.default_rng(seed)
    xs = np.linspace(0.0, 400.0, 4000)
    ys = np.cumsum(rng.normal(0.0, 0.05, len(xs)))
    breaks = [1200, 2900]
    keep = kept_track_points(xs, ys, breaks)
    assert np.array_equal(keep, _reference_keep(xs, ys, [0, 1199, 1200, 2899, 2900, len(xs) - 1]))
    assert keep[[0, 1199, 1200, 2899, 2900, len(xs) - 1]].all()
    assert 0 < keep.sum() < len(xs)
    kept_x, kept_y = xs[keep], ys[keep]
    assert max(_distance_to_polyline(x, y, kept_x, kept_y) for x, y in zip(xs, ys)) <= TRACK_PATH_TOLERANCE_PX


def test_a_closed_line_keeps_its_shape() -> None:
    """A Circular line ends where it starts: the chord between its ends has no length."""

    angles = np.linspace(-math.pi / 2, 1.5 * math.pi, 721)
    radii = 100.0 + 5.0 * np.sin(6 * angles)
    xs, ys = radii * np.cos(angles), radii * np.sin(angles)
    keep = kept_track_points(xs, ys)
    assert 20 < keep.sum() < len(xs)
    assert max(_distance_to_polyline(x, y, xs[keep], ys[keep]) for x, y in zip(xs, ys)) <= TRACK_PATH_TOLERANCE_PX


def test_collinear_points_are_dropped_and_short_lines_kept() -> None:
    assert track_line_to([0, 1, 2, 3], [5, 5, 5, 5]) == "L0 5L3 5"
    assert track_polyline([0, 1, 2], [0, 0.01, 0]) == "M0 0L2 0"
    assert track_polyline([0, 1, 2], [0, 0.03, 0]) == "M0 0L1 0.03L2 0"
    assert track_line_to([1.5], [2.5]) == "L1.5 2.5"
    assert track_line_to([], []) == ""


def test_circular_baselines_are_two_arcs() -> None:
    assert generate_circle_path_desc(50, 2) == "M 0 -100A 100 100 0 0 1 0 100A 100 100 0 0 1 0 -100z"
    # Back from 3 o'clock to 12 o'clock through the middle angle.
    assert circular_baseline_back(10.0, -math.pi / 2, 0.0) == "A10 10 0 0 0 7.071 -7.071A10 10 0 0 0 0 -10"


def test_circular_scalar_area_draws_the_outer_line_and_an_arc_baseline() -> None:
    plot = pd.DataFrame({"position": [0, 25, 50, 75], "value_normalized": [0.5, 0.5, 0.5, 0.5]})
    desc = generate_circular_scalar_area_path_desc(100, 100, plot, 20, 1.0)
    # Baseline radius 90, outer radius 100: the outer line keeps only the points the arc of samples needs,
    # the baseline returns from the last sample (270 degrees on screen) to the first as two arcs.
    assert desc.startswith("M0 -90L0 -100")
    assert desc.endswith("L-90 0A90 90 0 0 0 63.64 63.64A90 90 0 0 0 0 -90z")
    assert desc.count("A") == 2
