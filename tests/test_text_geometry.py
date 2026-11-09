from __future__ import annotations

import math

import pytest

from gbdraw.core.text import GlyphOutline, glyph_outlines, text_runs_glyph_outlines
from gbdraw.layout.text_geometry import (
    glyphs_across_circle,
    glyphs_along_circle,
    oriented_boxes_intersect,
    oriented_text_box,
    rotate_point,
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


def test_grown_convex_polygon_moves_every_edge_of_a_sector_part_out_by_the_same_distance() -> None:
    from gbdraw.layout.text_geometry import annular_sector_parts, grown_convex_polygon

    # A turned rectangle grows by the distance on every side; a sector part (a trapezoid) keeps it from every edge.
    rectangle = tuple(rotate_point(point, 30.0) for point in ((0.0, 0.0), (10.0, 0.0), (10.0, 4.0), (0.0, 4.0)))
    (part,) = annular_sector_parts(100.0, 112.0, -90.0, 2.0)
    for polygon in (rectangle, part, tuple(reversed(part))):
        grown = grown_convex_polygon(polygon, 3.0)
        for index, (x1, y1) in enumerate(polygon):
            x2, y2 = polygon[(index + 1) % len(polygon)]
            (gx, gy) = grown[index]
            # The grown corner's distance from the line of each edge it touches.
            distance = abs((x2 - x1) * (y1 - gy) - (x1 - gx) * (y2 - y1)) / math.hypot(x2 - x1, y2 - y1)
            assert distance == pytest.approx(3.0)


@pytest.mark.parametrize(
    ("radius", "tops_in", "meet"),
    (
        (1e6, True, False),
        (500.0, True, False),
        (200.0, True, True),
        (60.0, True, True),
        (200.0, False, False),
        (60.0, False, False),
    ),
)
def test_letters_that_stand_apart_on_a_line_close_up_at_their_tops_on_the_inner_side_of_a_circle(
    radius: float, tops_in: bool, meet: bool
) -> None:
    # "xx" at 60 px: the two x's stand a fraction of a px apart on a line; along a circle of 200 px with their tops
    # toward the centre their top corners cross, and with their tops away from it they part.
    glyphs = text_runs_glyph_outlines((("xx", "normal"),), "Liberation Sans", 60.0)

    assert glyphs.glyph_count == 2
    assert glyphs_along_circle(glyphs, radius, tops_in=tops_in).meet() is meet


def test_glyphs_meet_when_one_holds_the_other_or_their_turns_wrap_around() -> None:
    square = ((-10.0, 0.0), (10.0, 0.0), (10.0, 20.0), (-10.0, 20.0))
    small = ((-2.0, 5.0), (2.0, 5.0), (2.0, 9.0), (-2.0, 9.0))
    # A glyph inside another's ink with no edge crossing: they meet.
    assert glyphs_along_circle(glyph_outlines((GlyphOutline(0.0, (square,)), GlyphOutline(0.0, (small,)))), 100.0, tops_in=False).meet()
    # Inside a hole of the other (even-odd): they do not.
    ring = (square, ((-5.0, 2.0), (5.0, 2.0), (5.0, 18.0), (-5.0, 18.0)))
    assert not glyphs_along_circle(glyph_outlines((GlyphOutline(0.0, ring), GlyphOutline(0.0, (small,)))), 100.0, tops_in=False).meet()
    # Nearly a whole turn on, overlapping the first.
    turn = 2 * math.pi * 100.0
    assert glyphs_along_circle(glyph_outlines((GlyphOutline(0.0, (square,)), GlyphOutline(turn - 5.0, (square,)))), 100.0, tops_in=True).meet()
    assert not glyphs_along_circle(glyph_outlines((GlyphOutline(0.0, (square,)), GlyphOutline(turn / 2, (square,)))), 100.0, tops_in=True).meet()


def test_a_glyph_upright_on_a_circle_reaches_past_its_height_at_its_corners() -> None:
    # OV-2979: a 20 px wide, 20 px high square standing on a circle of 100 px is not bent to it. Tops away from the
    # centre: its bottom edge touches the circle at its middle and its top corners stand hypot(10, 120) out. Tops
    # toward the centre: its top edge comes 80 px near and its bottom corners stand hypot(10, 100) out.
    square = ((-10.0, 0.0), (10.0, 0.0), (10.0, 20.0), (-10.0, 20.0))
    glyphs = glyph_outlines((GlyphOutline(0.0, (square,)),))
    across = glyphs_across_circle(glyphs)
    assert across is not None
    assert across.radii(100.0, tops_in=False) == pytest.approx((100.0, math.hypot(10.0, 120.0)))
    assert across.radii(100.0, tops_in=True) == pytest.approx((80.0, math.hypot(10.0, 100.0)))
    # Two feet apart: the nearest ink is a foot's inner corner, off the ray through the glyph's middle.
    feet = glyph_outlines(
        (GlyphOutline(0.0, (((-15.0, 0.0), (-5.0, 0.0), (-5.0, 20.0), (-15.0, 20.0)), ((5.0, 0.0), (15.0, 0.0), (15.0, 20.0), (5.0, 20.0)))),)
    )
    across = glyphs_across_circle(feet)
    assert across is not None and across.radii(100.0, tops_in=False)[0] == pytest.approx(math.hypot(5.0, 100.0))
    assert glyphs_across_circle(glyph_outlines(())) is None


@pytest.mark.parametrize("tops_in", (False, True))
@pytest.mark.parametrize("text", ("R\u2026", "AVW xx", "Region 3 xxx", "gjpqy"))
def test_the_baseline_radius_reaching_a_glyph_radius_is_the_radius_it_was_measured_on(text: str, tops_in: bool) -> None:
    # Slanted and curved edges included: the radius that brings the nearest or farthest ink to a reach is the one
    # whose ink comes that near or reaches that far.
    outlines = text_runs_glyph_outlines(((text, "normal"),), "Liberation Sans", 60.0)
    across = glyphs_across_circle(outlines)
    assert across is not None
    points = outlines.points.tolist()
    side = -1.0 if tops_in else 1.0
    for radius in (60.0, 99.0, 260.0, 2000.0):
        near, far = across.radii(radius, tops_in=tops_in)
        assert near < far
        # Every edge and corner measured: the edges and corners left out never come nearer or reach farther.
        placed = [(x, radius + side * y) for x, y in points]
        nearest = math.inf
        for (x1, y1), (x2, y2) in zip(placed, (placed[index] for index in outlines.next_index.tolist())):
            length = (x2 - x1) ** 2 + (y2 - y1) ** 2
            t = min(1.0, max(0.0, -(x1 * (x2 - x1) + y1 * (y2 - y1)) / length)) if length > 0.0 else 0.0
            nearest = min(nearest, math.hypot(x1 + t * (x2 - x1), y1 + t * (y2 - y1)))
        assert near == pytest.approx(nearest, abs=1e-9)
        assert far == pytest.approx(max(math.hypot(x, y) for x, y in placed), abs=1e-9)
        assert across.radius_reaching(near, tops_in=tops_in, nearest=True) == pytest.approx(radius, abs=1e-6)
        assert across.radius_reaching(far, tops_in=tops_in, nearest=False) == pytest.approx(radius, abs=1e-6)
    # A corner that stands across the ray farther than the reach: no baseline radius brings it within.
    assert across.radius_reaching(1.0, tops_in=tops_in, nearest=False) == -math.inf
