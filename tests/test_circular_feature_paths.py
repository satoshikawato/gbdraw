from __future__ import annotations

import math
import re
from types import SimpleNamespace

import pytest
from svgwrite.container import Group

from gbdraw.features.objects import FeatureLocationPart, FeatureObject
from gbdraw.layout.circular import calculate_feature_position_factors_circular
from gbdraw.render.drawers.circular.features import FeatureDrawer
from gbdraw.render.drawers.circular.features import FeaturePathGenerator
from gbdraw.svg.circular_features import generate_circular_arrowhead_path
from gbdraw.svg.circular_features import generate_circular_arrow_path_with_radii
from gbdraw.svg.circular_features import (
    generate_circular_seven_vertex_arrowhead_path,
    generate_circular_seven_vertex_arrowhead_path_with_radii,
)
from gbdraw.svg.circular_features import generate_circular_intron_path
from gbdraw.svg.arrows import calculate_circular_arrow_length


TEST_RADIUS = 250.0
TEST_TRACK_RATIO = 0.19
TEST_CDS_RATIO = 0.03
TEST_OFFSET = 0.0
TEST_TRACK_TYPE = "tuckin"
TEST_STRANDEDNESS = True
TEST_TRACK_ID = 0
ARC_COMMAND_RE = re.compile(
    r"A\s*([-+0-9.eE]+),([-+0-9.eE]+)\s+([01])\s+([01])\s+([01])\s+([-+0-9.eE]+),([-+0-9.eE]+)"
)
MOVE_LINE_COMMAND_RE = re.compile(
    r"[ML]\s*([-+0-9.eE]+),([-+0-9.eE]+)"
)


def _make_feature(location: list[FeatureLocationPart], is_directional: bool = False) -> FeatureObject:
    return FeatureObject(
        feature_id="feature_000000001",
        location=location,
        is_directional=is_directional,
        color="#d3d3d3",
        note="",
        label_text="",
        coordinates=location,
        type="D-loop",
        qualifiers={},
    )


def _make_generator(total_length: int) -> FeaturePathGenerator:
    return FeaturePathGenerator(
        radius=TEST_RADIUS,
        total_length=total_length,
        track_ratio=TEST_TRACK_RATIO,
        cds_ratio=TEST_CDS_RATIO,
        offset=TEST_OFFSET,
        track_type=TEST_TRACK_TYPE,
        strandedness=TEST_STRANDEDNESS,
        track_id=TEST_TRACK_ID,
    )


def test_origin_spanning_non_directional_two_blocks_are_drawn_as_single_block() -> None:
    feature = _make_feature(
        [
            FeatureLocationPart("block", "001", "negative", 16023, 16569, False),
            FeatureLocationPart("line", "001", "negative", 577, 16022, False),
            FeatureLocationPart("block", "002", "negative", 0, 576, True),
        ],
        is_directional=False,
    )

    paths = _make_generator(total_length=16569).generate_circular_gene_path(feature)

    assert len(paths) == 1
    assert paths[0][0] == "block"


def test_origin_spanning_directional_two_blocks_are_drawn_as_single_block() -> None:
    feature = _make_feature(
        [
            FeatureLocationPart("block", "001", "negative", 16023, 16569, False),
            FeatureLocationPart("line", "001", "negative", 577, 16022, False),
            FeatureLocationPart("block", "002", "negative", 0, 576, True),
        ],
        is_directional=True,
    )

    paths = _make_generator(total_length=16569).generate_circular_gene_path(feature)

    assert len(paths) == 1
    assert [path[0] for path in paths] == ["block"]


def test_origin_spanning_directional_with_undefined_strand_falls_back_to_rectangle() -> None:
    feature = _make_feature(
        [
            FeatureLocationPart("block", "001", "undefined", 16023, 16569, False),
            FeatureLocationPart("line", "001", "undefined", 577, 16022, False),
            FeatureLocationPart("block", "002", "undefined", 0, 576, True),
        ],
        is_directional=True,
    )

    paths = _make_generator(total_length=16569).generate_circular_gene_path(feature)

    assert len(paths) == 1
    assert paths[0][0] == "block"


def test_non_origin_spanning_multipart_keeps_line_segment() -> None:
    feature = _make_feature(
        [
            FeatureLocationPart("block", "001", "negative", 900, 920, False),
            FeatureLocationPart("line", "001", "negative", 731, 899, False),
            FeatureLocationPart("block", "002", "negative", 700, 730, True),
        ],
        is_directional=False,
    )

    paths = _make_generator(total_length=1000).generate_circular_gene_path(feature)

    assert [path[0] for path in paths] == ["block", "line", "block"]


def _make_intron_path(strand: str, coord_start: int, coord_end: int, total_length: int = 1000) -> str:
    return generate_circular_intron_path(
        radius=TEST_RADIUS,
        coord_dict={
            "coord_strand": strand,
            "coord_start": coord_start,
            "coord_end": coord_end,
        },
        total_length=total_length,
        track_ratio=TEST_TRACK_RATIO,
        cds_ratio=TEST_CDS_RATIO,
        offset=TEST_OFFSET,
        track_type=TEST_TRACK_TYPE,
        strandedness=TEST_STRANDEDNESS,
        track_id=TEST_TRACK_ID,
    )[1]


def _extract_arc_commands(path_data: str) -> list[tuple[float, float, int, int]]:
    arc_commands: list[tuple[float, float, int, int]] = []
    for rx, ry, _rotation, large_arc, sweep, _x, _y in ARC_COMMAND_RE.findall(path_data):
        arc_commands.append((float(rx), float(ry), int(large_arc), int(sweep)))
    return arc_commands


def _extract_arc_commands_with_endpoints(
    path_data: str,
) -> list[tuple[float, float, int, int, float, float]]:
    arc_commands: list[tuple[float, float, int, int, float, float]] = []
    for rx, ry, _rotation, large_arc, sweep, x_end, y_end in ARC_COMMAND_RE.findall(path_data):
        arc_commands.append(
            (float(rx), float(ry), int(large_arc), int(sweep), float(x_end), float(y_end))
        )
    return arc_commands


def _position_from_xy(x: float, y: float, total_length: int) -> float:
    angle_deg = (math.degrees(math.atan2(y, x)) + 90.0) % 360.0
    return total_length * (angle_deg / 360.0)


def _path_command_endpoints(path_data: str) -> list[tuple[float, float]]:
    """Return path endpoints in command order, including arc endpoints."""
    endpoints: list[tuple[int, float, float]] = []
    for match in MOVE_LINE_COMMAND_RE.finditer(path_data):
        endpoints.append((match.start(), float(match.group(1)), float(match.group(2))))
    for match in ARC_COMMAND_RE.finditer(path_data):
        endpoints.append((match.start(), float(match.group(6)), float(match.group(7))))
    return [(x, y) for _offset, x, y in sorted(endpoints)]


def _radius(point: tuple[float, float]) -> float:
    return math.hypot(*point)


def _expected_intron_radius(strand: str, total_length: int = 1000) -> float:
    factors = calculate_feature_position_factors_circular(
        total_length=total_length,
        strand=strand,
        track_ratio=TEST_TRACK_RATIO,
        cds_ratio=TEST_CDS_RATIO,
        offset=TEST_OFFSET,
        track_type=TEST_TRACK_TYPE,
        strandedness=TEST_STRANDEDNESS,
        track_id=TEST_TRACK_ID,
    )
    return TEST_RADIUS * factors[1]


def test_middle_resolve_overlaps_displaces_negative_tracks_inward_symmetrically() -> None:
    positive = calculate_feature_position_factors_circular(
        total_length=1000,
        strand="positive",
        track_ratio=TEST_TRACK_RATIO,
        cds_ratio=TEST_CDS_RATIO,
        offset=TEST_OFFSET,
        track_type="middle",
        strandedness=False,
        track_id=2,
    )
    negative = calculate_feature_position_factors_circular(
        total_length=1000,
        strand="negative",
        track_ratio=TEST_TRACK_RATIO,
        cds_ratio=TEST_CDS_RATIO,
        offset=TEST_OFFSET,
        track_type="middle",
        strandedness=False,
        track_id=2,
    )
    undefined = calculate_feature_position_factors_circular(
        total_length=1000,
        strand="undefined",
        track_ratio=TEST_TRACK_RATIO,
        cds_ratio=TEST_CDS_RATIO,
        offset=TEST_OFFSET,
        track_type="middle",
        strandedness=False,
        track_id=2,
    )

    assert positive[1] > 1.0
    assert negative[1] < 1.0
    assert undefined[1] > 1.0
    assert math.isclose(positive[1] - 1.0, 1.0 - negative[1], rel_tol=1e-9, abs_tol=1e-9)


@pytest.mark.parametrize(
    ("strand", "coord_start", "coord_end", "expected_sweep"),
    [
        ("positive", 100, 800, 1),
        ("negative", 100, 800, 1),
    ],
)
def test_long_introns_are_split_and_stay_on_track_radius(
    strand: str,
    coord_start: int,
    coord_end: int,
    expected_sweep: int,
) -> None:
    path_data = _make_intron_path(strand, coord_start, coord_end)
    arc_commands = _extract_arc_commands(path_data)

    assert len(arc_commands) == 2
    expected_radius = _expected_intron_radius(strand)
    for rx, ry, large_arc, sweep in arc_commands:
        assert math.isclose(rx, expected_radius, rel_tol=1e-9, abs_tol=1e-9)
        assert math.isclose(ry, expected_radius, rel_tol=1e-9, abs_tol=1e-9)
        assert large_arc == 0
        assert sweep == expected_sweep


@pytest.mark.parametrize(
    ("strand", "coord_start", "coord_end", "expected_sweep"),
    [
        ("positive", 900, 100, 1),
        ("negative", 900, 100, 1),
    ],
)
def test_intron_uses_clockwise_short_arc_when_start_exceeds_end(
    strand: str,
    coord_start: int,
    coord_end: int,
    expected_sweep: int,
) -> None:
    path_data = _make_intron_path(strand, coord_start=coord_start, coord_end=coord_end)
    arc_commands = _extract_arc_commands(path_data)

    assert len(arc_commands) == 2
    assert all(command[3] == expected_sweep for command in arc_commands)
    assert all(command[2] == 0 for command in arc_commands)


def test_origin_spanning_negative_arrow_midpoint_wraps_short_way() -> None:
    total_length = 16569
    arrow_length = 40.0
    _kind, path_data = generate_circular_arrowhead_path(
        radius=TEST_RADIUS,
        coord_dict={
            "coord_strand": "negative",
            "coord_start": 16023,
            "coord_end": 576,
        },
        total_length=total_length,
        cds_arrow_length=arrow_length,
        track_ratio=TEST_TRACK_RATIO,
        cds_ratio=TEST_CDS_RATIO,
        offset=TEST_OFFSET,
        track_type=TEST_TRACK_TYPE,
        strandedness=TEST_STRANDEDNESS,
        track_id=TEST_TRACK_ID,
    )

    arc_commands = _extract_arc_commands_with_endpoints(path_data)
    assert len(arc_commands) >= 2

    first_outer_arc = arc_commands[0]
    assert first_outer_arc[2] == 0
    assert first_outer_arc[3] == 0

    midpoint_pos = _position_from_xy(first_outer_arc[4], first_outer_arc[5], total_length)
    arrow_start = 576.0
    arrow_end = 16023.0
    feature_len_bp = (arrow_start - arrow_end) % total_length
    shaft_len_bp = feature_len_bp - arrow_length
    expected_midpoint = (arrow_start - shaft_len_bp / 2.0) % total_length

    assert math.isclose(midpoint_pos, expected_midpoint, rel_tol=0.0, abs_tol=2.0)


@pytest.mark.parametrize("strand", ["positive", "negative"])
@pytest.mark.parametrize("shaft_width_ratio", [0.25, 0.5, 1.0])
def test_seven_vertex_arrowhead_uses_reduced_shaft_radii(
    strand: str,
    shaft_width_ratio: float,
) -> None:
    _kind, path_data = generate_circular_seven_vertex_arrowhead_path_with_radii(
        coord_dict={
            "coord_strand": strand,
            "coord_start": 100,
            "coord_end": 150,
        },
        total_length=1000,
        head_length_bp=10.0,
        inner_radius_px=80.0,
        center_radius_px=100.0,
        outer_radius_px=120.0,
        shaft_width_ratio=shaft_width_ratio,
    )

    endpoints = _path_command_endpoints(path_data)
    assert len(endpoints) == 7
    shaft_inner = 100.0 - (20.0 * shaft_width_ratio)
    shaft_outer = 100.0 + (20.0 * shaft_width_ratio)
    expected_radii = [
        shaft_inner,
        shaft_inner,
        80.0,
        100.0,
        120.0,
        shaft_outer,
        shaft_outer,
    ]
    assert [_radius(point) for point in endpoints] == pytest.approx(expected_radii)

    positions = [_position_from_xy(x, y, 1000) for x, y in endpoints]
    if strand == "positive":
        expected_positions = [100.0, 140.0, 140.0, 150.0, 140.0, 140.0, 100.0]
    else:
        expected_positions = [150.0, 110.0, 110.0, 100.0, 110.0, 110.0, 150.0]
    assert positions == pytest.approx(expected_positions)


@pytest.mark.parametrize("strand", ["positive", "negative"])
def test_factor_based_seven_vertex_arrowhead_matches_explicit_radii_api(
    strand: str,
) -> None:
    coord_dict = {
        "coord_strand": strand,
        "coord_start": 100,
        "coord_end": 150,
    }
    factors = calculate_feature_position_factors_circular(
        total_length=1000,
        strand=strand,
        track_ratio=TEST_TRACK_RATIO,
        cds_ratio=TEST_CDS_RATIO,
        offset=TEST_OFFSET,
        track_type=TEST_TRACK_TYPE,
        strandedness=TEST_STRANDEDNESS,
        track_id=TEST_TRACK_ID,
    )

    factor_based = generate_circular_seven_vertex_arrowhead_path(
        radius=TEST_RADIUS,
        coord_dict=coord_dict,
        total_length=1000,
        head_length_bp=10.0,
        track_ratio=TEST_TRACK_RATIO,
        cds_ratio=TEST_CDS_RATIO,
        offset=TEST_OFFSET,
        track_type=TEST_TRACK_TYPE,
        strandedness=TEST_STRANDEDNESS,
        shaft_width_ratio=0.5,
        track_id=TEST_TRACK_ID,
    )
    explicit_radii = generate_circular_seven_vertex_arrowhead_path_with_radii(
        coord_dict=coord_dict,
        total_length=1000,
        head_length_bp=10.0,
        inner_radius_px=TEST_RADIUS * factors[0],
        center_radius_px=TEST_RADIUS * factors[1],
        outer_radius_px=TEST_RADIUS * factors[2],
        shaft_width_ratio=0.5,
    )

    assert factor_based == explicit_radii
    assert len(_path_command_endpoints(factor_based[1])) == 7


@pytest.mark.parametrize(
    ("strand", "expected_sweeps"),
    [("positive", [1, 1, 0, 0]), ("negative", [0, 0, 1, 1])],
)
def test_long_factor_based_arrowhead_splits_both_shaft_arcs(
    strand: str,
    expected_sweeps: list[int],
) -> None:
    shaft_width_ratio = 0.5
    _kind, path_data = generate_circular_seven_vertex_arrowhead_path(
        radius=TEST_RADIUS,
        coord_dict={
            "coord_strand": strand,
            "coord_start": 100,
            "coord_end": 500,
        },
        total_length=1000,
        head_length_bp=10.0,
        track_ratio=TEST_TRACK_RATIO,
        cds_ratio=TEST_CDS_RATIO,
        offset=TEST_OFFSET,
        track_type=TEST_TRACK_TYPE,
        strandedness=TEST_STRANDEDNESS,
        shaft_width_ratio=shaft_width_ratio,
        track_id=TEST_TRACK_ID,
    )
    factors = calculate_feature_position_factors_circular(
        total_length=1000,
        strand=strand,
        track_ratio=TEST_TRACK_RATIO,
        cds_ratio=TEST_CDS_RATIO,
        offset=TEST_OFFSET,
        track_type=TEST_TRACK_TYPE,
        strandedness=TEST_STRANDEDNESS,
        track_id=TEST_TRACK_ID,
    )
    center_radius = TEST_RADIUS * factors[1]
    shaft_inner = center_radius + (TEST_RADIUS * factors[0] - center_radius) * 0.5
    shaft_outer = center_radius + (TEST_RADIUS * factors[2] - center_radius) * 0.5
    arcs = _extract_arc_commands(path_data)

    assert 360.0 * (400.0 - 10.0) / 1000.0 > 20.0
    assert len(arcs) == 4
    assert [arc[0] for arc in arcs] == pytest.approx(
        [shaft_inner, shaft_inner, shaft_outer, shaft_outer]
    )
    assert [arc[2] for arc in arcs] == [0, 0, 0, 0]
    assert [arc[3] for arc in arcs] == expected_sweeps


@pytest.mark.parametrize(
    ("strand", "expected_tip_position"),
    [("positive", 100.0), ("negative", 900.0)],
)
def test_origin_spanning_arrowhead_is_coalesced_with_reduced_shaft(
    strand: str,
    expected_tip_position: float,
) -> None:
    feature = FeatureObject(
        feature_id="feature_000000001",
        location=[
            FeatureLocationPart("block", "001", strand, 900, 1000, False),
            FeatureLocationPart("line", "001", strand, 101, 899, False),
            FeatureLocationPart("block", "002", strand, 0, 100, True),
        ],
        is_directional=None,
        color="#d3d3d3",
        note="",
        label_text="",
        coordinates=[],
        type="CDS",
        qualifiers={},
        glyph_kind="arrowhead",
    )

    paths = _make_generator(total_length=1000).generate_circular_gene_path(feature)
    endpoints = _path_command_endpoints(paths[0][1])
    factors = calculate_feature_position_factors_circular(
        total_length=1000,
        strand=strand,
        track_ratio=TEST_TRACK_RATIO,
        cds_ratio=TEST_CDS_RATIO,
        offset=TEST_OFFSET,
        track_type=TEST_TRACK_TYPE,
        strandedness=TEST_STRANDEDNESS,
        track_id=TEST_TRACK_ID,
    )
    full_inner = TEST_RADIUS * factors[0]
    center_radius = TEST_RADIUS * factors[1]
    full_outer = TEST_RADIUS * factors[2]
    shaft_inner = center_radius + (full_inner - center_radius) * 0.5
    shaft_outer = center_radius + (full_outer - center_radius) * 0.5

    assert len(paths) == 1
    assert len(_extract_arc_commands(paths[0][1])) == 4
    assert [_radius(point) for point in endpoints] == pytest.approx(
        [
            shaft_inner,
            shaft_inner,
            shaft_inner,
            full_inner,
            center_radius,
            full_outer,
            shaft_outer,
            shaft_outer,
            shaft_outer,
        ]
    )
    assert _position_from_xy(*endpoints[4], 1000) == pytest.approx(
        expected_tip_position
    )


@pytest.mark.parametrize("head_length_bp", [10.0, 20.0])
def test_arrowhead_without_positive_shaft_is_a_triangle(
    head_length_bp: float,
) -> None:
    _kind, path_data = generate_circular_seven_vertex_arrowhead_path_with_radii(
        coord_dict={
            "coord_strand": "positive",
            "coord_start": 100,
            "coord_end": 110,
        },
        total_length=1000,
        head_length_bp=head_length_bp,
        inner_radius_px=80.0,
        center_radius_px=100.0,
        outer_radius_px=120.0,
        shaft_width_ratio=0.5,
    )

    assert "A" not in path_data
    assert len(_path_command_endpoints(path_data)) == 3


def test_legacy_arrow_keeps_equality_boundary_as_five_vertex_path() -> None:
    _kind, path_data = generate_circular_arrow_path_with_radii(
        coord_dict={
            "coord_strand": "positive",
            "coord_start": 100,
            "coord_end": 110,
        },
        total_length=1000,
        cds_arrow_length=10.0,
        inner_radius_px=80.0,
        center_radius_px=100.0,
        outer_radius_px=120.0,
    )

    assert "A" in path_data
    assert len(_path_command_endpoints(path_data)) == 5


def test_multipart_arrowhead_uses_shaft_width_before_terminal_head() -> None:
    feature = FeatureObject(
        feature_id="feature_000000001",
        location=[
            FeatureLocationPart("block", "001", "positive", 100, 130, False),
            FeatureLocationPart("line", "001", "positive", 131, 159, False),
            FeatureLocationPart("block", "002", "positive", 160, 210, True),
        ],
        is_directional=None,
        color="#d3d3d3",
        note="",
        label_text="",
        coordinates=[],
        type="CDS",
        qualifiers={},
        glyph_kind="arrowhead",
    )
    generator = FeaturePathGenerator(
        radius=100.0,
        total_length=1000,
        track_ratio=0.2,
        cds_ratio=0.2,
        offset=0.0,
        track_type="middle",
        strandedness=False,
        head_length_ratio=0.5,
        shaft_width_ratio=0.5,
    )

    paths = generator.generate_circular_gene_path(feature)

    assert [path[0] for path in paths] == ["block", "line", "block"]
    first_block_radii = [_radius(point) for point in _path_command_endpoints(paths[0][1])]
    terminal_radii = [_radius(point) for point in _path_command_endpoints(paths[2][1])]
    assert max(first_block_radii) - min(first_block_radii) == pytest.approx(10.0)
    assert max(terminal_radii) - min(terminal_radii) == pytest.approx(20.0)


def test_undefined_strand_arrowhead_falls_back_to_rectangle() -> None:
    feature = FeatureObject(
        feature_id="feature_000000001",
        location=[FeatureLocationPart("block", "001", "undefined", 100, 150, True)],
        is_directional=None,
        color="#d3d3d3",
        note="",
        label_text="",
        coordinates=[],
        type="CDS",
        qualifiers={},
        glyph_kind="arrowhead",
    )

    path = _make_generator(total_length=1000).generate_circular_gene_path(feature)[0][1]

    assert len(_path_command_endpoints(path)) == 4


def test_circular_numeric_head_ratio_resolves_in_display_space_without_rounding() -> None:
    generator = FeaturePathGenerator(
        radius=100.0,
        total_length=1000,
        track_ratio=0.2,
        cds_ratio=0.2,
        offset=0.0,
        track_type="middle",
        strandedness=False,
        head_length_ratio=0.75,
    )
    inner_radius, center_radius, outer_radius = generator._feature_radii(
        "positive", None
    )

    resolved = generator._resolved_arrow_length("positive", None)

    expected_head_px = abs(outer_radius - inner_radius) * 0.75
    expected_head_bp = expected_head_px * 1000.0 / (2.0 * math.pi * center_radius)
    assert resolved == pytest.approx(expected_head_bp)
    assert not float(resolved).is_integer()


def test_circular_auto_head_length_keeps_legacy_bp_value_exactly() -> None:
    generator = _make_generator(total_length=1000)

    assert generator._resolved_arrow_length(
        "positive", None
    ) == calculate_circular_arrow_length(1000)


@pytest.mark.parametrize(
    ("strand", "coord_start", "coord_end"),
    [
        ("positive", 501, 500),
        ("negative", 501, 500),
    ],
)
def test_intron_one_bp_gap_is_not_rendered_as_near_full_circle(
    strand: str,
    coord_start: int,
    coord_end: int,
) -> None:
    path_data = _make_intron_path(strand, coord_start=coord_start, coord_end=coord_end)
    arc_commands = _extract_arc_commands(path_data)

    assert len(arc_commands) == 1
    rx, ry, large_arc, sweep = arc_commands[0]
    expected_radius = _expected_intron_radius(strand)
    assert math.isclose(rx, expected_radius, rel_tol=1e-9, abs_tol=1e-9)
    assert math.isclose(ry, expected_radius, rel_tol=1e-9, abs_tol=1e-9)
    assert large_arc == 0
    assert sweep == 1


def test_short_introns_remain_single_arc_and_track_aligned() -> None:
    path_data = _make_intron_path("positive", coord_start=100, coord_end=110)
    arc_commands = _extract_arc_commands(path_data)

    assert len(arc_commands) == 1
    rx, ry, large_arc, sweep = arc_commands[0]
    expected_radius = _expected_intron_radius("positive")
    assert math.isclose(rx, expected_radius, rel_tol=1e-9, abs_tol=1e-9)
    assert math.isclose(ry, expected_radius, rel_tol=1e-9, abs_tol=1e-9)
    assert large_arc == 0
    assert sweep == 1


def test_feature_drawer_draws_lines_before_blocks_within_same_feature() -> None:
    feature = _make_feature(
        [
            FeatureLocationPart("block", "001", "positive", 100, 140, False),
            FeatureLocationPart("line", "001", "positive", 141, 199, False),
            FeatureLocationPart("block", "002", "positive", 200, 240, True),
        ],
        is_directional=False,
    )
    feature.color = "#4aa3df"
    cfg = SimpleNamespace(
        block_fill_color="#d3d3d3",
        block_stroke_color="#000000",
        block_stroke_width=1.0,
        line_stroke_color="#8f8f8f",
        line_stroke_width=1.0,
    )
    drawer = FeatureDrawer(cfg)
    group = Group(id="test")

    group = drawer.draw(
        feature,
        group,
        total_length=1000,
        radius=TEST_RADIUS,
        track_ratio=TEST_TRACK_RATIO,
        track_ratio_factor=1.0,
        track_type=TEST_TRACK_TYPE,
        strandedness=TEST_STRANDEDNESS,
        length_param="short",
    )

    paths = [element for element in group.elements if element.elementname == "path"]
    assert len(paths) == 3
    fills = [str(path.attribs.get("fill", "")) for path in paths]
    assert fills[0] == "none"
    assert fills[1] == feature.color
    assert fills[2] == feature.color
