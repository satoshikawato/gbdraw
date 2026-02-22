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
from gbdraw.svg.circular_features import generate_circular_intron_path


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
