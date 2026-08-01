from __future__ import annotations

import re

import pytest

from gbdraw.features.objects import FeatureLocationPart, FeatureObject
from gbdraw.layout.linear import measure_linear_feature_lanes
from gbdraw.render.drawers.linear.features import FeaturePathGenerator
from gbdraw.svg.linear_features import (
    create_arrow_path_linear,
    create_arrowhead_path_linear,
    create_seven_vertex_arrowhead_path_linear,
)


_POINT_RE = re.compile(
    r"[ML]\s*(-?(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?),"
    r"(-?(?:\d+(?:\.\d*)?|\.\d+)(?:e[+-]?\d+)?)",
    re.IGNORECASE,
)


def _points(path_data: str) -> list[tuple[float, float]]:
    return [(float(x), float(y)) for x, y in _POINT_RE.findall(path_data)]


def _linear_path(
    *,
    strand: str = "positive",
    start: int = 100,
    end: int = 300,
    **overrides,
) -> dict:
    values = {
        "coord_dict": {
            "feat_start": start,
            "feat_end": end,
            "feat_strand": strand,
        },
        "arrow_length": 50.0,
        "cds_height": 20.0,
        "feature_strand": strand,
        "genome_length": 1000,
        "alignment_width": 800.0,
        "genome_size_normalization_factor": 1.0,
        "separate_strands": False,
        "feature_y_positions": (0.0, 10.0, 20.0),
    }
    values.update(overrides)
    return values


def _feature(
    *,
    glyph_kind: str,
    strand: str,
    parts: list[FeatureLocationPart],
    track_id: int = 0,
) -> FeatureObject:
    feature = FeatureObject(
        feature_id="feature_000000001",
        location=parts,
        is_directional=None,
        glyph_kind=glyph_kind,
        color="#cccccc",
        note="",
        label_text="",
        coordinates=parts,
        type="CDS",
        qualifiers={},
    )
    feature.feature_track_id = track_id
    assert feature.strand == strand
    return feature


@pytest.mark.parametrize(
    ("strand", "expected"),
    [
        (
            "positive",
            "M 80.0,0.0 L 200.0,0.0 L 240.0,10.0 "
            "L 200.0,20.0 L 80.0,20.0 z",
        ),
        (
            "negative",
            "M 240.0,0.0 L 120.0,0.0 L 80.0,10.0 "
            "L 120.0,20.0 L 240.0,20.0 z",
        ),
    ],
)
def test_legacy_arrow_path_text_is_unchanged(strand: str, expected: str) -> None:
    arguments = _linear_path(strand=strand)

    assert create_arrow_path_linear(**arguments) == ["block", expected]
    assert create_arrowhead_path_linear(**arguments) == ["block", expected]


@pytest.mark.parametrize(
    ("strand", "expected"),
    [
        (
            "positive",
            [
                (100.0, 5.0),
                (260.0, 5.0),
                (260.0, 0.0),
                (300.0, 10.0),
                (260.0, 20.0),
                (260.0, 15.0),
                (100.0, 15.0),
            ],
        ),
        (
            "negative",
            [
                (300.0, 5.0),
                (140.0, 5.0),
                (140.0, 0.0),
                (100.0, 10.0),
                (140.0, 20.0),
                (140.0, 15.0),
                (300.0, 15.0),
            ],
        ),
    ],
)
def test_arrowhead_has_seven_mirrored_vertices(
    strand: str,
    expected: list[tuple[float, float]],
) -> None:
    result = create_seven_vertex_arrowhead_path_linear(
        **_linear_path(
            strand=strand,
            arrow_length=40.0,
            alignment_width=1000.0,
        )
    )

    assert _points(result[1]) == expected
    assert len(set(_points(result[1]))) == 7


@pytest.mark.parametrize(
    ("shaft_width_ratio", "expected_bounds"),
    [(0.25, (7.5, 12.5)), (0.5, (5.0, 15.0)), (1.0, (0.0, 20.0))],
)
def test_arrowhead_shaft_width_ratios(
    shaft_width_ratio: float,
    expected_bounds: tuple[float, float],
) -> None:
    result = create_seven_vertex_arrowhead_path_linear(
        **_linear_path(
            arrow_length=40.0,
            alignment_width=1000.0,
            shaft_width_ratio=shaft_width_ratio,
        )
    )
    points = _points(result[1])

    assert points[0][1] == pytest.approx(expected_bounds[0])
    assert points[1][1] == pytest.approx(expected_bounds[0])
    assert points[5][1] == pytest.approx(expected_bounds[1])
    assert points[6][1] == pytest.approx(expected_bounds[1])
    assert points[3][1] == pytest.approx(10.0)


def test_full_width_arrowhead_has_same_visible_outline_as_arrow() -> None:
    arguments = _linear_path(
        arrow_length=40.0,
        alignment_width=1000.0,
    )
    arrow_points = _points(create_arrow_path_linear(**arguments)[1])
    arrowhead_points = _points(
        create_seven_vertex_arrowhead_path_linear(
            **arguments,
            shaft_width_ratio=1.0,
        )[1]
    )
    deduplicated_arrowhead = [
        point
        for index, point in enumerate(arrowhead_points)
        if index == 0 or point != arrowhead_points[index - 1]
    ]

    assert deduplicated_arrowhead == arrow_points


@pytest.mark.parametrize("normalization_factor", [0.5, 1.0, 1.75])
def test_numeric_head_ratio_is_resolved_after_record_normalization(
    normalization_factor: float,
) -> None:
    result = create_seven_vertex_arrowhead_path_linear(
        **_linear_path(
            start=100,
            end=500,
            genome_size_normalization_factor=normalization_factor,
            head_length_ratio=2.0,
        )
    )
    points = _points(result[1])

    assert abs(points[3][0] - points[2][0]) == pytest.approx(40.0)


@pytest.mark.parametrize(
    ("feature_length", "legacy_vertex_count", "arrowhead_vertex_count"),
    [(40, 3, 3), (50, 5, 3), (100, 5, 7)],
)
def test_short_and_equality_boundaries_are_shape_specific(
    feature_length: int,
    legacy_vertex_count: int,
    arrowhead_vertex_count: int,
) -> None:
    arguments = _linear_path(
        start=100,
        end=100 + feature_length,
        alignment_width=1000.0,
        arrow_length=50.0,
    )

    assert len(_points(create_arrow_path_linear(**arguments)[1])) == legacy_vertex_count
    assert (
        len(_points(create_seven_vertex_arrowhead_path_linear(**arguments)[1]))
        == arrowhead_vertex_count
    )


@pytest.mark.parametrize(
    ("strand", "track_id"),
    [("positive", 2), ("negative", -2)],
)
def test_multipart_arrowhead_uses_lane_center_and_reduced_width_shafts(
    strand: str,
    track_id: int,
) -> None:
    parts = [
        FeatureLocationPart("block", "001", strand, 100, 180, False),
        FeatureLocationPart("line", "001", strand, 181, 219, False),
        FeatureLocationPart("block", "002", strand, 220, 400, True),
    ]
    feature = _feature(
        glyph_kind="arrowhead",
        strand=strand,
        parts=parts,
        track_id=track_id,
    )
    lane_geometry = measure_linear_feature_lanes(
        {"feature": feature},
        cds_height=20.0,
        separate_strands=True,
    )
    lane = lane_geometry.lane_for(
        strand=strand,
        track_id=track_id,
        separate_strands=True,
    )
    generator = FeaturePathGenerator(
        genome_length=1000,
        alignment_width=1000.0,
        cds_height=20.0,
        genome_size_normalization_factor=1.0,
        feature_strand=strand,
        separate_strands=True,
        arrow_length=40.0,
        feature_lane_geometry=lane_geometry,
        shaft_width_ratio=0.5,
    )

    paths = generator.generate_linear_gene_path(feature)
    first_block = _points(paths[0][1])
    connector = _points(paths[1][1])
    terminal_block = _points(paths[2][1])
    expected_shaft_top = lane.middle_y + (lane.top_y - lane.middle_y) * 0.5
    expected_shaft_bottom = lane.middle_y + (lane.bottom_y - lane.middle_y) * 0.5

    assert sorted({point[1] for point in first_block}) == pytest.approx(
        sorted((expected_shaft_top, expected_shaft_bottom))
    )
    assert [point[1] for point in connector] == pytest.approx(
        [lane.middle_y, lane.middle_y]
    )
    assert terminal_block[0][1] == pytest.approx(expected_shaft_top)
    assert terminal_block[2][1] == pytest.approx(lane.top_y)
    assert terminal_block[3][1] == pytest.approx(lane.middle_y)
    assert terminal_block[4][1] == pytest.approx(lane.bottom_y)
    assert terminal_block[6][1] == pytest.approx(expected_shaft_bottom)


def test_undefined_strand_arrowhead_falls_back_to_full_width_rectangle() -> None:
    part = FeatureLocationPart("block", "001", "undefined", 100, 300, True)
    feature = _feature(
        glyph_kind="arrowhead",
        strand="undefined",
        parts=[part],
    )
    generator = FeaturePathGenerator(
        genome_length=1000,
        alignment_width=1000.0,
        cds_height=20.0,
        genome_size_normalization_factor=1.0,
        feature_strand="undefined",
        separate_strands=False,
        arrow_length=40.0,
        shaft_width_ratio=0.25,
    )

    points = _points(generator.generate_linear_gene_path(feature)[0][1])

    assert len(points) == 4
    assert {point[1] for point in points} == {-10.0, 10.0}
