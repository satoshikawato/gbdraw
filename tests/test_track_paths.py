"""The one track-path writer: RDP at 0.02 SVG px per drawn segment, Circular baselines as two arcs
(incremental-render OV-1230 D-144..D-147)."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest

import re

from gbdraw.svg.circular_tracks import (
    generate_circle_path_desc,
    generate_circular_gc_content_path_desc,
    generate_circular_gc_skew_path_desc,
    generate_circular_scalar_area_path_desc,
)
from gbdraw.svg.track_paths import (
    TRACK_PATH_TOLERANCE_PX,
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
    # A full circle is 120 + 240 degrees (no drawn arc near 180 or 360; OV-1840).
    assert generate_circle_path_desc(50, 2) == "M 0,-100 A 100,100 0 0 1 86.603,50 A 100,100 0 1 1 0,-100z"


def test_non_finite_samples_are_rejected() -> None:
    with pytest.raises(ValueError, match="finite"):
        kept_track_points(np.array([0.0, 1.0, 2.0]), np.array([0.0, np.nan, 0.0]))
    with pytest.raises(ValueError, match="finite"):
        track_line_to([0.0, np.inf, 2.0], [0.0, 1.0, 0.0])


def test_circular_scalar_area_draws_the_outer_line_and_an_arc_baseline() -> None:
    plot = pd.DataFrame({"position": [0, 25, 50, 75], "value_normalized": [0.5, 0.5, 0.5, 0.5]})
    desc = generate_circular_scalar_area_path_desc(100, 100, plot, 20, 1.0)
    # Baseline radius 90, outer radius 100: the outer line keeps only the points the arc of samples needs,
    # the baseline returns from the last sample (9 o'clock) to the first (12 o'clock) as two 135-degree arcs.
    assert desc.startswith("M 0,-90L0 -100")
    assert desc.endswith("L -90,0 A 90,90 0 0 0 63.64,63.64 A 90,90 0 0 0 0,-90z")
    assert desc.count("A") == 2


def test_index_arrays_are_platform_integers(monkeypatch: pytest.MonkeyPatch) -> None:
    """Pyodide's NumPy is 32-bit (intp is int32) and refuses int64 indices in ``ufunc.reduceat``; every index
    array the simplification builds must be ``np.intp`` so the Web build draws the same path."""

    import gbdraw.svg.track_paths as track_paths

    real_reduceat = np.maximum.reduceat

    class _Maximum:
        @staticmethod
        def reduceat(values, indices, *args, **kwargs):
            assert np.asarray(indices).dtype == np.intp
            return real_reduceat(values, indices.astype(np.int32), *args, **kwargs)

    class _NumPy:
        maximum = _Maximum

        def __getattr__(self, name):
            return getattr(np, name)

    monkeypatch.setattr(track_paths, "np", _NumPy())
    monkeypatch.setattr(np, "intp", np.int32)
    xs = np.linspace(0.0, 50.0, 400)
    ys = np.sin(xs)
    assert track_paths.kept_track_points(xs, ys, (130, 260)).sum() > 4


# OV-2863: a Circular track line joins its samples along the circle. A short record has few windows (670 bp
# at the Auto 100 bp step is 7 samples, 51 degrees apart); straight chords between them drew a polygon.

_PATH_ROUNDING_PX = 0.001
# The added points keep each piece within the tolerance of the line, and the simplification then moves the
# outline by at most the tolerance again.
_LINE_BOUND_PX = 2 * TRACK_PATH_TOLERANCE_PX + _PATH_ROUNDING_PX


def _first_subpath_vertices(desc: str) -> np.ndarray:
    """The M/L vertices of the first subpath of a path (up to its first ``z`` or arc)."""

    first = re.split(r"[zA]", desc, maxsplit=1)[0]
    pairs = re.findall(r"[ML]\s*(-?[0-9.]+)[ ,](-?[0-9.]+)", first)
    return np.array([[float(x), float(y)] for x, y in pairs])


def _polar_curve(positions: np.ndarray, radii: np.ndarray, total: float, per_sample: int = 400) -> np.ndarray:
    """Points of the line a Circular track draws: angle and radius change linearly between samples."""

    t = np.linspace(0.0, 1.0, per_sample, endpoint=False)
    pos = np.concatenate([a + (b - a) * t for a, b in zip(positions[:-1], positions[1:])] + [positions[-1:]])
    rad = np.concatenate([a + (b - a) * t for a, b in zip(radii[:-1], radii[1:])] + [radii[-1:]])
    angles = np.radians(360.0 * pos / total - 90.0)
    return np.column_stack((rad * np.cos(angles), rad * np.sin(angles)))


def _max_deviation(curve: np.ndarray, vertices: np.ndarray) -> float:
    closed = np.vstack((vertices, vertices[:1]))
    return max(_distance_to_polyline(x, y, closed[:, 0], closed[:, 1]) for x, y in curve)


def _short_gc_df() -> pd.DataFrame:
    starts = np.arange(0, 670, 100)
    content = np.array([0.42, 0.47, 0.44, 0.50, 0.46, 0.43, 0.48])
    skew = np.array([0.05, -0.10, 0.02, -0.04, 0.08, -0.06, 0.01])
    return pd.DataFrame({"GC content": content, "GC skew": skew}, index=starts)


def test_short_record_gc_content_follows_the_circle() -> None:
    df = _short_gc_df()
    desc = generate_circular_gc_content_path_desc(300, 670, df, 60, 1.0, "GC")
    diff = df["GC content"] - df["GC content"].mean()
    radii = (300 + 30 * diff / diff.abs().max()).to_numpy()
    positions = np.append(df.index.to_numpy(dtype=float), 670.0)
    curve = _polar_curve(positions, np.append(radii, radii[0]), 670.0)
    assert _max_deviation(curve, _first_subpath_vertices(desc)) <= _LINE_BOUND_PX


def test_short_record_gc_skew_follows_the_circle_through_its_close() -> None:
    df = _short_gc_df()
    desc = generate_circular_gc_skew_path_desc(200, df, 670, 60, 1.0, "GC")
    diff = df["GC skew"] - df["GC skew"].mean()
    radii = (200 + 30 * diff / diff.abs().max()).to_numpy()
    # The line closes from the last window (600) back to 12 o'clock on the radius of the first.
    positions = np.append(df.index.to_numpy(dtype=float), 670.0)
    curve = _polar_curve(positions, np.append(radii, radii[0]), 670.0)
    assert _max_deviation(curve, _first_subpath_vertices(desc)) <= _LINE_BOUND_PX


def test_a_constant_short_line_is_a_circle() -> None:
    df = pd.DataFrame({"GC content": [0.5] * 7, "GC skew": [0.0] * 7}, index=np.arange(0, 670, 100))
    for desc in (
        generate_circular_gc_content_path_desc(300, 670, df, 60, 1.0, "GC"),
        generate_circular_gc_skew_path_desc(300, df, 670, 60, 1.0, "GC"),
    ):
        vertices = _first_subpath_vertices(desc)
        angles = np.linspace(0.0, 2 * math.pi, 3600, endpoint=False)
        circle = np.column_stack((300 * np.cos(angles), 300 * np.sin(angles)))
        assert _max_deviation(circle, vertices) <= _LINE_BOUND_PX


def test_short_record_scalar_area_outer_line_follows_the_circle() -> None:
    plot = pd.DataFrame({"position": np.arange(0, 670, 100), "value_normalized": [0.2, 0.9, 0.4, 0.7, 0.1, 0.6, 0.3]})
    desc = generate_circular_scalar_area_path_desc(300, 670, plot, 60, 1.0, close_at_record_len=True)
    outer = _first_subpath_vertices(desc)[1:-1]  # drop the baseline ends
    positions = np.append(plot["position"].to_numpy(dtype=float), 670.0)
    radii = 270 + 60 * np.append(plot["value_normalized"].to_numpy(), 0.2)
    curve = _polar_curve(positions, radii, 670.0)
    deviation = max(_distance_to_polyline(x, y, outer[:, 0], outer[:, 1]) for x, y in curve)
    assert deviation <= _LINE_BOUND_PX


def test_densely_sampled_lines_gain_no_points() -> None:
    """Samples already closer than the tolerance draw exactly the polyline through them (long records keep
    their paths)."""

    total = 300_000
    positions = np.arange(0, total, 100, dtype=float)
    rng = np.random.default_rng(3)
    df = pd.DataFrame({"GC content": 0.5 + rng.normal(0, 0.02, len(positions))}, index=positions.astype(int))
    desc = generate_circular_gc_content_path_desc(300, total, df, 60, 1.0, "GC")
    diff = df["GC content"] - df["GC content"].mean()
    radii = (300 + 30 * diff / diff.abs().max()).to_numpy()
    line_positions = np.append(positions, float(total))
    line_radii = np.append(radii, radii[0])
    angles = np.radians(360.0 * line_positions / total - 90)
    xs = np.concatenate(([line_radii[0] * math.cos(math.radians(-90))], line_radii * np.cos(angles)))
    ys = np.concatenate(([line_radii[0] * math.sin(math.radians(-90))], line_radii * np.sin(angles)))
    assert desc.startswith(track_polyline(xs, ys) + "z")
