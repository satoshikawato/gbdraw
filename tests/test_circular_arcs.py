"""Every Circular arc gbdraw draws stays on its circle once its numbers are rounded (OV-1840).

An SVG arc is drawn from its rounded endpoints; near 180 or 360 degrees a 0.0005 px endpoint error moved the
drawn arc up to 0.8 px off its circle. ``gbdraw.svg.circular_arcs`` splits every arc so that none is near
those spans; these tests parse what each writer emits and bound the radial deviation of every drawn arc.
"""

from __future__ import annotations

import math
import re

import pandas as pd
import pytest

from gbdraw.render.drawers.circular.annotations import CircularLabelArc, _arc_path, _band_path
from gbdraw.svg.circular_arcs import arc_pieces, circle_arc
from gbdraw.svg.circular_conservation import generate_annular_hsp_path_desc, generate_full_annulus_path_desc
from gbdraw.svg.circular_features import (
    generate_circular_arrow_path_with_radii,
    generate_circular_intron_path_with_radii,
    generate_circular_narrow_arrow_path_with_radii,
    generate_circular_rectangle_path_with_radii,
)
from gbdraw.svg.circular_tracks import (
    circular_scalar_area_line,
    generate_circle_path_desc,
    generate_circular_scalar_area_path_desc,
)
from gbdraw.svg.path_data import PATH_DECIMALS

# The largest radial distance of a drawn arc from its true (unrounded) circle, in SVG px. The drawn arc has the
# rounded radius and endpoints, so it is off by up to 0.0005 px of rounding plus the bow its centre takes from them.
# Measured over the radii of `_RADII` (20..5000 px, most unrounded), every start in 0.81-degree steps and the spans of
# `_SPANS` plus spans next to 20 and 300 degrees, the worst case is 0.0046 px (r = 1845.4195327 from 294.84 degrees
# along 300 degrees, drawn as 90 + 210); it does not grow with the radius. That is inside
# D-146's 0.02 px simplification plus 0.0005 px rounding; the old equal halves of a near-full circle reached 0.82 px.
ARC_DEVIATION_BOUND_PX = 0.005
_TOKEN = re.compile(r"[MLAZz]|-?\d+(?:\.\d+)?")


def _arcs(d: str):
    """Each arc of a path as (start point, drawn radius, large-arc flag, sweep flag, end point)."""

    tokens = _TOKEN.findall(d)
    i, command, current, start = 0, "", (0.0, 0.0), (0.0, 0.0)
    while i < len(tokens):
        if tokens[i] in "MLAZz":
            command = tokens[i]
            i += 1
            if command in "Zz":
                current = start
            continue
        if command == "M":
            current = start = (float(tokens[i]), float(tokens[i + 1]))
            i += 2
            command = "L"
        elif command == "L":
            current = (float(tokens[i]), float(tokens[i + 1]))
            i += 2
        else:
            rx, _ry, _rotation, large, sweep, x, y = tokens[i : i + 7]
            end = (float(x), float(y))
            yield current, float(rx), int(large), int(sweep), end
            current = end
            i += 7


def _radial_deviation(start, drawn_radius, large, sweep, end, radius) -> float:
    """The largest distance from the origin-centred circle of ``radius`` along the arc drawn with ``drawn_radius``
    (the rounded ``rx``) between the rounded endpoints (SVG F.6.5)."""

    (x1, y1), (x2, y2) = start, end
    dx, dy = (x1 - x2) / 2.0, (y1 - y2) / 2.0
    if dx == 0.0 and dy == 0.0:
        return 0.0
    r = drawn_radius
    scale = (dx * dx + dy * dy) / (r * r)
    if scale > 1.0:
        r *= math.sqrt(scale)
    root = math.sqrt(max(0.0, (r * r * r * r - r * r * (dx * dx + dy * dy)) / (r * r * (dx * dx + dy * dy))))
    sign = -1.0 if large == sweep else 1.0
    cx = sign * root * dy + (x1 + x2) / 2.0
    cy = sign * root * -dx + (y1 + y2) / 2.0
    a1 = math.atan2(y1 - cy, x1 - cx)
    turn = math.atan2(y2 - cy, x2 - cx) - a1
    if sweep and turn < 0:
        turn += 2.0 * math.pi
    if not sweep and turn > 0:
        turn -= 2.0 * math.pi
    return max(
        abs(math.hypot(cx + r * math.cos(a1 + turn * k / 400), cy + r * math.sin(a1 + turn * k / 400)) - radius)
        for k in range(401)
    )


def _worst(d: str, radii: tuple[float, ...]) -> float:
    """The worst deviation of the arcs of ``d``, each from the one of ``radii`` (the true radii the writer was given)
    nearest its drawn radius."""

    return max(
        (_radial_deviation(*arc, min(radii, key=lambda radius: abs(radius - arc[1]))) for arc in _arcs(d)),
        default=0.0,
    )


_SPANS = (
    0.5, 10.0, 19.99, 20.0, 20.01, 90.0, 150.0, 179.9, 180.0, 200.0, 299.9, 299.99, 300.0, 300.01, 358.0, 359.8,
    359.97, 360.0,
)
_RADII = (20.0, 203.83452012050364, 280.15, 312.9515104679078, 387.257, 1000.0, 1845.4195327, 5000.0)


def test_no_drawn_piece_is_near_180_or_360_degrees() -> None:
    for span in _SPANS + tuple(k * 0.37 for k in range(973)):
        for signed in (span, -span):
            pieces = arc_pieces(signed)
            assert math.isclose(sum(pieces), signed, abs_tol=1e-9)
            for piece in pieces:
                assert abs(piece) <= 20.0 or abs(abs(piece) - 180.0) >= 30.0 - 1e-9
                assert abs(piece) <= 240.0 + 1e-9


@pytest.mark.parametrize("radius", _RADII)
@pytest.mark.parametrize("start", [-90.0, -37.3, 0.0, 123.4, 269.97])
def test_every_arc_stays_on_its_circle(radius: float, start: float) -> None:
    for span in _SPANS:
        for signed in (span, -span):
            assert _worst(circle_arc(radius, start, signed), (radius,)) <= ARC_DEVIATION_BOUND_PX


def test_the_review_worst_case_stays_on_its_circle() -> None:
    # r = 387.257, the baseline back from 269.97 degrees to -90 (the last Depth bin one step before the end).
    assert _worst(circle_arc(387.257, 269.97, -359.97), (387.257,)) <= ARC_DEVIATION_BOUND_PX
    # The measured worst case of the bound above.
    assert 0.0045 < _worst(circle_arc(1845.4195327, 294.84, 300.0), (1845.4195327,)) <= ARC_DEVIATION_BOUND_PX


@pytest.mark.parametrize("span", [360.000001, -360.5, 720.0, math.inf, math.nan])
def test_arc_pieces_rejects_more_than_a_full_circle(span: float) -> None:
    with pytest.raises(ValueError, match="360"):
        arc_pieces(span)


_LENGTH = 16_569
_FEATURE_SPANS = (10, 1_000, 8_280, 8_285, 13_000, 16_000, 16_550, 16_568)


@pytest.mark.parametrize("strand", ["positive", "negative"])
@pytest.mark.parametrize("length_bp", _FEATURE_SPANS)
def test_feature_writers_draw_arcs_on_their_circles(strand: str, length_bp: int) -> None:
    coords = {"coord_strand": strand, "coord_start": 5, "coord_end": 5 + length_bp}
    inner, center, outer = 376.83452012050364, 383.9515104679078, 391.06849
    paths = [
        generate_circular_rectangle_path_with_radii(coords, _LENGTH, inner, center, outer)[1],
        generate_circular_arrow_path_with_radii(coords, _LENGTH, 60.0, inner, center, outer)[1],
        generate_circular_narrow_arrow_path_with_radii(coords, _LENGTH, 60.0, inner, center, outer, 0.5)[1],
        generate_circular_intron_path_with_radii(coords, _LENGTH, center)[1],
    ]
    # The arrows draw their body between the middles of the band halves.
    radii = (inner, center, outer, (inner + center) / 2.0, (center + outer) / 2.0)
    for d in paths:
        assert _worst(d, radii) <= ARC_DEVIATION_BOUND_PX, d


def test_conservation_annotation_label_and_track_arcs_stay_on_their_circles() -> None:
    inner, outer, label, track = 390.50120364, 394.32104679078, 410.3712345, 387.25683452012
    depth = pd.DataFrame({"position": [0, 6_000, 12_000, 16_568], "value_normalized": [0.2, 0.9, 0.4, 0.6]})
    area_line = circular_scalar_area_line(track, _LENGTH, depth, 30.0, 1.0)
    assert area_line is not None
    paths = [
        (
            generate_annular_hsp_path_desc(
                draw_start=5, draw_end=16_555, total_length=_LENGTH, inner_radius_px=inner, outer_radius_px=outer
            ),
            (inner, outer),
        ),
        (generate_full_annulus_path_desc(inner_radius_px=inner, outer_radius_px=outer, total_length=_LENGTH), (inner, outer)),
        (_arc_path(3.0, 16_560.0, _LENGTH, label), (label,)),
        (_band_path(3.0, 16_560.0, _LENGTH, inner, outer), (inner, outer)),
        (CircularLabelArc(radius=label, angle=12.3, half_span=89.95, flipped=False).path, (label,)),
        (CircularLabelArc(radius=label, angle=171.0, half_span=89.95, flipped=True).path, (label,)),
        (generate_circle_path_desc(track, 1.0), (track,)),
        # The area closes on its baseline circle.
        (generate_circular_scalar_area_path_desc(track, _LENGTH, depth, 30.0, 1.0), (area_line[2],)),
    ]
    for d, radii in paths:
        assert "A" in d
        assert _worst(d, radii) <= ARC_DEVIATION_BOUND_PX, d


def test_arc_numbers_use_the_path_number_form() -> None:
    d = circle_arc(387.257, 269.97, -359.97)
    numbers = re.findall(r"-?\d+(?:\.\d+)?", d)
    assert all(len(number.partition(".")[2]) <= PATH_DECIMALS for number in numbers)
