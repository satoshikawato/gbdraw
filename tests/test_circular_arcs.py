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
from gbdraw.svg.circular_tracks import generate_circle_path_desc, generate_circular_scalar_area_path_desc
from gbdraw.svg.path_data import PATH_DECIMALS

# The largest radial distance of a drawn arc from its circle, in SVG px. Measured over radii 20..5000 px and every
# start in 0.27-degree steps, the worst case is 0.0025 px (two 149.9-degree halves of a 299.9-degree span; it does
# not grow with the radius). That is far inside D-146's 0.02 px simplification plus 0.0005 px rounding; the old
# equal halves of a near-full circle reached 0.82 px.
ARC_DEVIATION_BOUND_PX = 0.003
_TOKEN = re.compile(r"[MLAZz]|-?\d+(?:\.\d+)?")


def _arcs(d: str):
    """Each arc of a path as (start point, radius, large-arc flag, sweep flag, end point)."""

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


def _radial_deviation(start, radius, large, sweep, end) -> float:
    """The largest distance from the origin-centred circle of ``radius`` along the drawn arc (SVG F.6.5)."""

    (x1, y1), (x2, y2) = start, end
    dx, dy = (x1 - x2) / 2.0, (y1 - y2) / 2.0
    if dx == 0.0 and dy == 0.0:
        return 0.0
    r = radius
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


def _worst(d: str) -> float:
    return max((_radial_deviation(*arc) for arc in _arcs(d)), default=0.0)


_SPANS = (0.5, 10.0, 20.0, 20.01, 90.0, 150.0, 179.9, 180.0, 200.0, 299.9, 300.0, 358.0, 359.8, 359.97, 360.0)


def test_no_drawn_piece_is_near_180_or_360_degrees() -> None:
    for span in _SPANS + tuple(k * 0.37 for k in range(973)):
        for signed in (span, -span):
            pieces = arc_pieces(signed)
            assert math.isclose(sum(pieces), signed, abs_tol=1e-9)
            for piece in pieces:
                assert abs(piece) <= 20.0 or abs(abs(piece) - 180.0) >= 30.0 - 1e-9
                assert abs(piece) <= 240.0 + 1e-9


@pytest.mark.parametrize("radius", [20.0, 280.15, 387.257, 1000.0, 5000.0])
@pytest.mark.parametrize("start", [-90.0, -37.3, 0.0, 123.4, 269.97])
def test_every_arc_stays_on_its_circle(radius: float, start: float) -> None:
    for span in _SPANS:
        for signed in (span, -span):
            assert _worst(circle_arc(radius, start, signed)) <= ARC_DEVIATION_BOUND_PX


def test_the_review_worst_case_stays_on_its_circle() -> None:
    # r = 387.257, the baseline back from 269.97 degrees to -90 (the last Depth bin one step before the end).
    assert _worst(circle_arc(387.257, 269.97, -359.97)) <= ARC_DEVIATION_BOUND_PX


_LENGTH = 16_569
_FEATURE_SPANS = (10, 1_000, 8_280, 8_285, 13_000, 16_000, 16_550, 16_568)


@pytest.mark.parametrize("strand", ["positive", "negative"])
@pytest.mark.parametrize("length_bp", _FEATURE_SPANS)
def test_feature_writers_draw_arcs_on_their_circles(strand: str, length_bp: int) -> None:
    coords = {"coord_strand": strand, "coord_start": 5, "coord_end": 5 + length_bp}
    paths = [
        generate_circular_rectangle_path_with_radii(coords, _LENGTH, 380.123, 387.257, 394.391)[1],
        generate_circular_arrow_path_with_radii(coords, _LENGTH, 60.0, 380.123, 387.257, 394.391)[1],
        generate_circular_narrow_arrow_path_with_radii(coords, _LENGTH, 60.0, 380.123, 387.257, 394.391, 0.5)[1],
        generate_circular_intron_path_with_radii(coords, _LENGTH, 387.257)[1],
    ]
    for d in paths:
        assert _worst(d) <= ARC_DEVIATION_BOUND_PX, d


def test_conservation_annotation_label_and_track_arcs_stay_on_their_circles() -> None:
    paths = [
        generate_annular_hsp_path_desc(
            draw_start=5, draw_end=16_555, total_length=_LENGTH, inner_radius_px=390.5, outer_radius_px=394.321
        ),
        generate_full_annulus_path_desc(inner_radius_px=390.5, outer_radius_px=394.321, total_length=_LENGTH),
        _arc_path(3.0, 16_560.0, _LENGTH, 401.777),
        _band_path(3.0, 16_560.0, _LENGTH, 398.1, 405.9),
        CircularLabelArc(radius=410.37, angle=12.3, half_span=89.95, flipped=False).path,
        CircularLabelArc(radius=410.37, angle=171.0, half_span=89.95, flipped=True).path,
        generate_circle_path_desc(387.257, 1.0),
    ]
    depth = pd.DataFrame({"position": [0, 6_000, 12_000, 16_568], "value_normalized": [0.2, 0.9, 0.4, 0.6]})
    paths.append(generate_circular_scalar_area_path_desc(387.257, _LENGTH, depth, 30.0, 1.0))
    for d in paths:
        assert "A" in d
        assert _worst(d) <= ARC_DEVIATION_BOUND_PX, d


def test_arc_numbers_use_the_path_number_form() -> None:
    d = circle_arc(387.257, 269.97, -359.97)
    numbers = re.findall(r"-?\d+(?:\.\d+)?", d)
    assert all(len(number.partition(".")[2]) <= PATH_DECIMALS for number in numbers)
