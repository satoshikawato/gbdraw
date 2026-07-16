from __future__ import annotations

import math

import pytest

from gbdraw.layout.text_geometry import (
    oriented_boxes_intersect,
    oriented_text_box,
    segment_intersects_polygon,
)


@pytest.mark.parametrize(
    ("anchor", "expected"),
    [
        ("start", (10.0, 30.0)),
        ("middle", (0.0, 20.0)),
        ("end", (-10.0, 10.0)),
    ],
)
def test_oriented_text_box_anchor_bounds(anchor: str, expected: tuple[float, float]) -> None:
    box = oriented_text_box(10.0, 20.0, 20.0, 10.0, 0.0, anchor)
    assert box.aabb[0] == pytest.approx(expected[0])
    assert box.aabb[2] == pytest.approx(expected[1])
    assert box.aabb[1:] != ()


@pytest.mark.parametrize(
    ("rotation_deg", "expected"),
    [
        (0.0, (-10.0, -5.0, 10.0, 5.0)),
        (90.0, (-5.0, -10.0, 5.0, 10.0)),
        (180.0, (-10.0, -5.0, 10.0, 5.0)),
        (270.0, (-5.0, -10.0, 5.0, 10.0)),
    ],
)
def test_oriented_text_box_cardinal_rotation_aabb(
    rotation_deg: float,
    expected: tuple[float, float, float, float],
) -> None:
    box = oriented_text_box(0.0, 0.0, 20.0, 10.0, rotation_deg, "middle")
    assert box.aabb == pytest.approx(expected)


def test_oriented_text_box_padding_expands_all_sides() -> None:
    box = oriented_text_box(0.0, 0.0, 20.0, 10.0, 0.0, "middle", padding_px=2.0)
    assert box.aabb == pytest.approx((-12.0, -7.0, 12.0, 7.0))


def test_oriented_box_sat_and_edge_touch_semantics() -> None:
    first = oriented_text_box(0.0, 0.0, 20.0, 10.0, 30.0, "middle")
    same = oriented_text_box(0.0, 0.0, 20.0, 10.0, -20.0, "middle")
    far = oriented_text_box(100.0, 0.0, 20.0, 10.0, 0.0, "middle")
    touching = oriented_text_box(20.0, 0.0, 20.0, 10.0, 0.0, "middle")
    unrotated = oriented_text_box(0.0, 0.0, 20.0, 10.0, 0.0, "middle")

    assert oriented_boxes_intersect(first, same)
    assert not oriented_boxes_intersect(first, far)
    assert not oriented_boxes_intersect(unrotated, touching)
    assert oriented_boxes_intersect(unrotated, touching, touching=True)


def test_segment_intersects_rotated_text_box() -> None:
    box = oriented_text_box(0.0, 0.0, 20.0, 8.0, 45.0, "middle")
    assert segment_intersects_polygon(((-20.0, 0.0), (20.0, 0.0)), box.corners)
    assert not segment_intersects_polygon(((30.0, 30.0), (40.0, 40.0)), box.corners)
    assert all(math.isfinite(value) for point in box.corners for value in point)


def test_segment_corner_touch_semantics() -> None:
    box = oriented_text_box(0.0, 0.0, 20.0, 10.0)
    segment = ((10.0, 5.0), (20.0, 10.0))
    assert segment_intersects_polygon(segment, box.corners, touching=True)
    assert not segment_intersects_polygon(segment, box.corners, touching=False)
