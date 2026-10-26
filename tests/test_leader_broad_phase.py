"""The leader broad phases: crossing leaders and leaders over label text.

Both count exactly what an every-pair test counts, and both hand the exact test
only shapes that come near each other: a radial leader's long box no longer
makes every leader in that quarter of the map a candidate (OV-1946).
"""

from __future__ import annotations

import math
import random
from collections.abc import Callable, Sequence

import pytest

import gbdraw.labels.circular_radial as circular_radial
import gbdraw.labels.fit as fit
from gbdraw.labels.circular_leaders import leader_segments
from gbdraw.labels.circular_types import CircularLabel
from gbdraw.layout.text_geometry import Point, segment_intersects_polygon

BUCKET_PX = 64.0
# Two shapes that share a 64 px bucket of the covers come within two bucket diagonals of each other.
NEAR_PX = 2.0 * math.sqrt(2.0) * BUCKET_PX + 1e-6

Segment = tuple[Point, Point]


def _point_segment_distance(point: Point, segment: Segment) -> float:
    (ax, ay), (bx, by) = segment
    dx, dy = bx - ax, by - ay
    length_sq = dx * dx + dy * dy
    t = 0.0 if length_sq == 0.0 else max(0.0, min(1.0, ((point[0] - ax) * dx + (point[1] - ay) * dy) / length_sq))
    return math.hypot(point[0] - (ax + t * dx), point[1] - (ay + t * dy))


def _segment_distance(first: Segment, second: Segment) -> float:
    if fit.segments_intersect(first, second, touching=True):
        return 0.0
    return min(
        _point_segment_distance(first[0], second),
        _point_segment_distance(first[1], second),
        _point_segment_distance(second[0], first),
        _point_segment_distance(second[1], first),
    )


def _polyline_distance(first: Sequence[Point], second: Sequence[Point]) -> float:
    return min(
        _segment_distance(left, right)
        for left in zip(first, first[1:])
        for right in zip(second, second[1:])
    )


def _radial_fan(rng: random.Random, count: int, ring_px: float, jitter_deg: float) -> list[tuple[Point, ...]]:
    """Leaders from a feature ring at 350 px out to a label ring, as radial placement draws them."""
    polylines = []
    for index in range(count):
        angle = 2.0 * math.pi * index / max(1, count) + rng.uniform(-0.01, 0.01)
        out = angle + math.radians(rng.uniform(-jitter_deg, jitter_deg))
        polylines.append(
            (
                (350.0 * math.cos(angle), 350.0 * math.sin(angle)),
                (380.0 * math.cos(angle), 380.0 * math.sin(angle)),
                (ring_px * math.cos(out), ring_px * math.sin(out)),
            )
        )
    return polylines


def _random_polylines(rng: random.Random) -> list[tuple[Point, ...]]:
    scale = rng.choice([10.0, 100.0, 1000.0, 10000.0])
    polylines: list[tuple[Point, ...]] = []
    for _ in range(rng.randint(0, 40)):
        if rng.random() < 0.3:
            polylines.extend(_radial_fan(rng, 1, scale, 20.0))
            continue
        points = [(rng.uniform(-scale, scale), rng.uniform(-scale, scale)) for _ in range(rng.choice([0, 1, 2, 3, 4]))]
        if rng.random() < 0.3:
            # Points on bucket lines.
            points = [(round(x / BUCKET_PX) * BUCKET_PX, round(y / BUCKET_PX) * BUCKET_PX) for x, y in points]
        polylines.append(tuple(points))
    return polylines


def _every_crossing_pair(polylines: Sequence[Sequence[Point]]) -> list[tuple[int, int]]:
    return [
        (left, right)
        for left in range(len(polylines))
        for right in range(left + 1, len(polylines))
        if len(polylines[left]) >= 2 and len(polylines[right]) >= 2 and fit.polylines_cross(polylines[left], polylines[right])
    ]


def test_crossing_pairs_equal_every_pair_tested() -> None:
    rng = random.Random(1946)
    cases: list[list[tuple[Point, ...]]] = [_random_polylines(rng) for _ in range(250)]
    cases += [_radial_fan(rng, rng.randint(2, 80), rng.choice([500.0, 3000.0, 9000.0]), 15.0) for _ in range(20)]
    # Leaders that meet exactly on bucket lines and corners.
    cases.append([((0.0, -64.0), (0.0, 64.0)), ((-64.0, 0.0), (64.0, 0.0)), ((-128.0, -128.0), (128.0, 128.0))])
    crossings = 0
    for polylines in cases:
        expected = _every_crossing_pair(polylines)
        crossings += len(expected)
        assert fit.crossing_pairs(polylines) == expected
    assert crossings > 1000


def _counting(monkeypatch: pytest.MonkeyPatch, module: object, name: str) -> list[tuple[object, object]]:
    seen: list[tuple[object, object]] = []
    exact: Callable[..., bool] = getattr(module, name)

    def counted(first: object, second: object, **kwargs: bool) -> bool:
        seen.append((first, second))
        return exact(first, second, **kwargs)

    monkeypatch.setattr(module, name, counted)
    return seen


def test_crossing_pairs_keep_a_crossing_the_cut_point_rounding_could_hide() -> None:
    # The cut points of the long segment are rounded; without SEGMENT_PIECE_PAD (on every
    # piece) the cell of the crossing is missing from its cover and the pair is lost.
    long_leader = ((-37213.81971366056, 18715.974944659632), (1314.7739683781765, -661.2402822229997))
    tiny_leader = ((-2.13453691728464e-12, -6.225769573961266e-13), (-100.0, -100.0))

    assert fit.polylines_cross([*long_leader], [*tiny_leader])
    assert fit.crossing_pairs([list(long_leader), list(tiny_leader)]) == [(0, 1)]


def test_crossing_pairs_test_only_leaders_that_come_near(monkeypatch: pytest.MonkeyPatch) -> None:
    polylines = _radial_fan(random.Random(3), 120, 6000.0, 0.0)
    tested = _counting(monkeypatch, fit, "polylines_cross")

    assert fit.crossing_pairs(polylines) == []
    assert tested
    assert max(_polyline_distance(first, second) for first, second in tested) <= NEAR_PX  # type: ignore[arg-type]


def _text_label(index: int, leader: Sequence[Point], rotation_deg: float) -> CircularLabel:
    feature, elbow, contact = leader
    return {
        "stable_id": f"feature-{index}",
        "label_text": f"label {index}",
        "feature_anchor_x": feature[0],
        "feature_anchor_y": feature[1],
        "middle_x": elbow[0],
        "middle_y": elbow[1],
        "leader_start_x": contact[0],
        "leader_start_y": contact[1],
        "text_x": contact[0],
        "text_y": contact[1],
        "width_px": 90.0,
        "height_px": 12.0,
        "rotation_deg": rotation_deg,
        "text_anchor": "start",
    }


def _radial_text_labels(rng: random.Random, count: int, ring_px: float) -> list[CircularLabel]:
    labels = []
    for index, leader in enumerate(_radial_fan(rng, count, ring_px, 25.0)):
        # Text at several radii, so leaders pass through other labels' text.
        scale = rng.uniform(0.5, 1.0)
        contact = (leader[2][0] * scale, leader[2][1] * scale)
        leader = (leader[0], leader[1], contact)
        labels.append(_text_label(index, leader, math.degrees(math.atan2(contact[1], contact[0])) + rng.uniform(-30, 30)))
    return labels


def _every_leader_text_hit(labels: Sequence[CircularLabel]) -> int:
    corners = [circular_radial._box_for_label(label).corners for label in labels]
    return sum(
        segment_intersects_polygon(segment, corners[other], touching=False)
        for index, label in enumerate(labels)
        for segment in leader_segments(label)
        for other in range(len(labels))
        if other != index
    )


def test_leader_text_collision_count_equals_every_pair_tested() -> None:
    rng = random.Random(284)
    hits = 0
    for _ in range(40):
        labels = _radial_text_labels(rng, rng.randint(0, 60), rng.choice([400.0, 1500.0, 6000.0]))
        expected = _every_leader_text_hit(labels)
        hits += expected
        assert circular_radial._leader_text_collision_count(labels) == expected
    assert hits > 100


def test_leader_text_collision_count_tests_only_text_near_the_leader(monkeypatch: pytest.MonkeyPatch) -> None:
    labels = _radial_text_labels(random.Random(5), 120, 6000.0)
    boxes = [circular_radial._box_for_label(label) for label in labels]
    tested = _counting(monkeypatch, circular_radial, "segment_intersects_polygon")

    circular_radial._leader_text_collision_count(labels)

    assert tested
    for segment, corners in tested:
        box = next(box for box in boxes if box.corners == corners)
        min_x, min_y, max_x, max_y = box.aabb
        edges = [((min_x, min_y), (max_x, min_y)), ((max_x, min_y), (max_x, max_y)), ((max_x, max_y), (min_x, max_y)), ((min_x, max_y), (min_x, min_y))]
        inside = min_x <= segment[0][0] <= max_x and min_y <= segment[0][1] <= max_y  # type: ignore[index]
        assert inside or min(_segment_distance(segment, edge) for edge in edges) <= NEAR_PX  # type: ignore[arg-type]
