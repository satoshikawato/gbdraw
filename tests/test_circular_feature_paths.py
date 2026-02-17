from __future__ import annotations

from gbdraw.features.objects import FeatureLocationPart, FeatureObject
from gbdraw.render.drawers.circular.features import FeaturePathGenerator


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
        radius=250.0,
        total_length=total_length,
        track_ratio=0.19,
        cds_ratio=0.03,
        offset=0.0,
        track_type="tuckin",
        strandedness=True,
        track_id=0,
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
