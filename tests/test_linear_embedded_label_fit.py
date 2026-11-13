"""A Linear label drawn inside its feature keeps its letters inside the drawn outline (Owner D-161).

The label's ink box, widened by 0.125 em at each end, lies inside the outline of the block it sits on, as
the block's path draws it; it stays at the middle where it fits there, else slides to the nearest x that
fits, and goes outside (or is not drawn, inside-only) when no x fits. A part a later feature draws over the
block is avoided too (OV-3380, OV-3397).
"""

from __future__ import annotations

import re

import pytest

from gbdraw.canvas.linear import LinearCanvasGeometry
from gbdraw.config.models import GbdrawConfig, LinearRenderProfile
from gbdraw.config.modify import modify_config_dict
from gbdraw.config.toml import load_config_toml
from gbdraw.core.text import text_ink_offsets
from gbdraw.features.objects import FeatureLocationPart, FeatureObject
from gbdraw.labels.linear import linear_record_labels
from gbdraw.render.drawers.linear.features import FeaturePathGenerator
from gbdraw.svg.arrows import LinearBlockOutline, linear_block_glyph, linear_block_outline
from gbdraw.svg.linear_features import (
    create_arrow_path_linear,
    create_arrow_shaft_path_linear,
    create_rectangle_path_linear,
)

FONT_SIZE = 16.0
MARGIN_PX = 0.125 * FONT_SIZE
_POINT = re.compile(r"[ML]\s*(-?[\d.]+),(-?[\d.]+)")


def _points(path: str) -> list[tuple[float, float]]:
    return [(float(x), float(y)) for x, y in _POINT.findall(path)]


def _outline_points(outline: LinearBlockOutline, full_width: bool) -> list[tuple[float, float]]:
    """The vertices the outline says the block has, in the order its path draws them."""

    o = outline
    if o.glyph != "arrow":
        return [(o.tail_x, o.shaft_top_y), (o.tip_x, o.shaft_top_y), (o.tip_x, o.shaft_bottom_y), (o.tail_x, o.shaft_bottom_y)]
    if o.triangle:
        return [(o.tail_x, o.top_y), (o.tip_x, o.middle_y), (o.tail_x, o.bottom_y)]
    if full_width:
        return [(o.tail_x, o.top_y), (o.shoulder_x, o.top_y), (o.tip_x, o.middle_y), (o.shoulder_x, o.bottom_y), (o.tail_x, o.bottom_y)]
    return [
        (o.tail_x, o.shaft_top_y), (o.shoulder_x, o.shaft_top_y), (o.shoulder_x, o.top_y), (o.tip_x, o.middle_y),
        (o.shoulder_x, o.bottom_y), (o.shoulder_x, o.shaft_bottom_y), (o.tail_x, o.shaft_bottom_y),
    ]


_Y = (0.0, 10.0, 20.0)


def _path_arguments(start: int, end: int, strand: str) -> dict:
    return {
        "coord_dict": {"feat_start": start, "feat_end": end, "feat_strand": strand},
        "genome_length": 1000, "alignment_width": 800.0, "cds_height": 20.0, "feature_strand": strand,
        "separate_strands": False, "feature_track_id": 0, "feature_y_positions": _Y,
    }


@pytest.mark.parametrize(
    ("glyph", "start", "end", "strand", "head_ratio", "shaft_ratio"),
    [
        ("arrow", 100, 300, "positive", "auto", 1.0),
        ("arrow", 100, 300, "negative", "auto", 1.0),
        ("arrow", 100, 130, "positive", "auto", 1.0),  # shorter than the head: the legacy triangle
        ("arrow", 100, 300, "negative", 0.64, 0.6),
        ("arrow", 100, 108, "positive", 0.64, 0.6),  # no shaft left: a triangle
        ("shaft", 100, 300, "positive", "auto", 0.6),
        ("rectangle", 100, 300, "undefined", "auto", 1.0),
    ],
)
def test_the_label_fit_reads_the_outline_the_block_is_drawn_with(glyph, start, end, strand, head_ratio, shaft_ratio) -> None:
    arguments = _path_arguments(start, end, strand)
    if glyph == "arrow":
        path = create_arrow_path_linear(
            **arguments, arrow_length=50.0, head_length_ratio=head_ratio, shaft_width_ratio=shaft_ratio
        )[1]
    elif glyph == "shaft":
        path = create_arrow_shaft_path_linear(**arguments, shaft_width_ratio=shaft_ratio)[1]
    else:
        path = create_rectangle_path_linear(**arguments)[1]
    outline = linear_block_outline(glyph, start * 0.8, end * 0.8, strand, _Y, 40.0, head_ratio, shaft_ratio)

    drawn, expected = _points(path), _outline_points(outline, shaft_ratio == 1.0)
    assert len(drawn) == len(expected)
    # pytest.approx compares flat sequences only (pytest 9.2 rejects nested tuples).
    assert [c for point in drawn for c in point] == pytest.approx([c for point in expected for c in point], abs=5e-4)


def test_every_block_of_a_multipart_arrow_is_drawn_as_its_glyph_says() -> None:
    parts = [
        FeatureLocationPart("block", "001", "positive", 100, 200, False),
        FeatureLocationPart("line", "001", "positive", 201, 299, False),
        FeatureLocationPart("block", "002", "positive", 300, 400, True),
    ]
    feature = _feature("multi", parts, "")
    generator = FeaturePathGenerator(1000, 800.0, 20.0, "positive", False, 50.0, feature_lane_geometry=None,
                                     head_length_ratio="auto", shaft_width_ratio=0.6)
    paths = generator.generate_linear_gene_path(feature)
    glyphs = [linear_block_glyph("arrow", part.strand, part.is_last, 0.6) for part in parts if part.kind == "block"]
    assert glyphs == ["shaft", "arrow"]
    assert [len(_points(path[1])) for path in paths if path[0] == "block"] == [4, 7]


def _config(overrides: dict[str, object] | None = None) -> GbdrawConfig:
    values: dict[str, object] = {
        "labels.filtering.blacklist_keywords": [],
        "labels.linear.scope": "all",
    }
    values.update(overrides or {})
    return GbdrawConfig.from_dict(modify_config_dict(load_config_toml("gbdraw.data", "config.toml"), values))


def _feature(feature_id: str, parts: list[FeatureLocationPart], label: str) -> FeatureObject:
    return FeatureObject(
        feature_id=feature_id, location=parts, is_directional=True, color="#54bcf8", note="", label_text=label,
        coordinates=[part for part in parts if part.kind == "block"], type="CDS", qualifiers={"product": [label]},
    )


def _arrow(feature_id: str, start: int, end: int, label: str, strand: str = "positive") -> FeatureObject:
    return _feature(feature_id, [FeatureLocationPart("block", "001", strand, start, end, True)], label)


def _labels(features: dict[str, FeatureObject], arrow_length: float, config: GbdrawConfig | None = None) -> dict:
    """The labels by text; one px per bp, a 20 px tall track and heads of ``arrow_length`` px."""

    geometry = LinearCanvasGeometry(
        alignment_width=1000.0, longest_genome=1000, length_param="short", cds_height=20.0,
        arrow_length=arrow_length, track_axis_gap=None, normalize_length=False,
    )
    record = linear_record_labels(
        features, 1000, 1000.0, geometry, False, "middle", LinearRenderProfile(config or _config()),
        label_font_size=FONT_SIZE,
    )
    return {label["label_text"]: label for label in record.labels}


def _ink(text: str) -> tuple[float, float, float, float]:
    ink = text_ink_offsets(text, _config().objects.text.font_family, FONT_SIZE, dominant_baseline="central",
                           text_anchor="middle")
    assert ink is not None
    return ink


def _inside_convex(polygon: list[tuple[float, float]], point: tuple[float, float]) -> bool:
    signs = set()
    for (x1, y1), (x2, y2) in zip(polygon, polygon[1:] + polygon[:1]):
        cross = (x2 - x1) * (point[1] - y1) - (y2 - y1) * (point[0] - x1)
        if abs(cross) > 1e-6:
            signs.add(cross > 0)
    return len(signs) <= 1


def _box_inside_drawn_arrow(label: dict, anchor_x: float, start: int, end: int, strand: str, arrow_length: float) -> bool:
    """Whether the widened ink box at ``anchor_x`` lies inside the arrow the feature's path draws."""

    y_positions = (label["feature_top_y"], label["feature_middle_y"], label["feature_bottom_y"])
    path = create_arrow_path_linear(
        {"feat_start": start, "feat_end": end, "feat_strand": strand}, arrow_length, 20.0, strand, 1000, 1000.0,
        False, feature_y_positions=y_positions,
    )[1]
    polygon = _points(path)
    left, right, top, bottom = _ink(label["label_text"])
    middle_y = label["middle_y"]
    corners = [(anchor_x + x, middle_y + y) for x in (left - MARGIN_PX, right + MARGIN_PX) for y in (top, bottom)]
    return all(_inside_convex(polygon, corner) for corner in corners)


@pytest.mark.linear
@pytest.mark.parametrize("strand", ["positive", "negative"])
def test_a_label_that_touches_the_head_at_the_middle_slides_toward_the_tail(strand: str) -> None:
    text, head = "ea59", 30.0
    left, right, top, bottom = _ink(text)
    share = max(-top, bottom) / 10.0
    width = right - left + 2 * MARGIN_PX
    length = int(width + head * share) + 3
    assert length < 2 * (max(right, -left) + MARGIN_PX + head * share) - 2, "the middle must not fit"
    start = 400
    end = start + length
    label = _labels({"f": _arrow("f", start, end, text, strand)}, head)[text]

    assert label["is_embedded"]
    middle = label["feature_anchor_x"]
    toward_tail = -1.0 if strand == "positive" else 1.0
    assert (label["middle_x"] - middle) * toward_tail > 0.5
    assert not _box_inside_drawn_arrow(label, middle, start, end, strand, head)
    assert _box_inside_drawn_arrow(label, label["middle_x"], start, end, strand, head)
    # The nearest x that fits: half a px nearer the middle touches the outline.
    assert not _box_inside_drawn_arrow(label, label["middle_x"] - 0.5 * toward_tail, start, end, strand, head)


@pytest.mark.linear
def test_a_label_that_fits_at_the_middle_stays_there() -> None:
    label = _labels({"f": _arrow("f", 100, 400, "int")}, 30.0)["int"]

    assert label["is_embedded"]
    assert label["middle_x"] == label["feature_anchor_x"]


@pytest.mark.linear
@pytest.mark.parametrize("rendering", ["auto", "embedded_only"])
def test_a_label_that_fits_nowhere_inside_goes_outside_or_is_not_drawn(rendering: str) -> None:
    text = "ea59"
    left, right, _top, _bottom = _ink(text)
    # Long enough for the advance and its 0.25 em (Owner D-20), but all head: the ink meets the outline anywhere.
    length = int(right - left + 4 * MARGIN_PX) + 4
    labels = _labels({"f": _arrow("f", 100, 100 + length, text)}, float(length + 5), _config({"labels.rendering": rendering}))

    if rendering == "auto":
        assert not labels[text]["is_embedded"]
    else:
        assert text not in labels


@pytest.mark.linear
def test_a_shaft_too_narrow_for_the_letters_and_a_short_head_send_the_label_outside() -> None:
    text = "ggt"
    narrow = _config({
        "objects.features.arrow_geometry.shaft_width_ratio": 0.6,
        "objects.features.arrow_geometry.head_length_ratio": 0.64,
    })
    full = _config({"objects.features.arrow_geometry.head_length_ratio": 0.64})
    features = {"f": _arrow("f", 100, 160, text)}

    assert not _labels(features, 30.0, narrow)[text]["is_embedded"]
    assert _labels(features, 30.0, full)[text]["is_embedded"]


@pytest.mark.linear
def test_a_label_keeps_off_the_outline_of_a_feature_drawn_after_it() -> None:
    text = "C"
    left, right, _top, _bottom = _ink(text)
    long_feature = _arrow("long", 100, 400, text)
    later = _arrow("later", 240, 262, "")
    after = _labels({"long": long_feature, "later": later}, 30.0)[text]
    before = _labels({"later": later, "long": long_feature}, 30.0)[text]

    assert before["is_embedded"] and before["middle_x"] == before["feature_anchor_x"]
    assert after["is_embedded"]
    box = (after["middle_x"] + left - MARGIN_PX, after["middle_x"] + right + MARGIN_PX)
    assert box[1] <= 240.0 + 1e-6 or box[0] >= 262.0 - 1e-6
    assert after["middle_x"] != after["feature_anchor_x"]
