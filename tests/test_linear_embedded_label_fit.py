"""A Linear label drawn inside its feature keeps its letters inside the drawn outline (Owner D-161; PD-OI-097
scenario revision 29, EMB-OUTLINE-SLIDE).

The label's ink box, widened at each end by half the 0.25 em a label needs inside its feature (Owner D-20),
lies inside the outline of the block it sits on, as the block's path draws it; it stays at the middle where
it fits there, else slides to the nearest x that fits, and goes outside (or is not drawn, inside-only) when
no x fits. A part a later feature draws over the block in the same lane and wrap row is avoided too
(OV-3380, OV-3397).
"""

from __future__ import annotations

import re
from xml.etree import ElementTree as ET

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from gbdraw.canvas.linear import LinearCanvasGeometry
from gbdraw.config.models import GbdrawConfig, LinearRenderProfile
from gbdraw.config.modify import modify_config_dict
from gbdraw.config.toml import load_config_toml
from gbdraw.core.text import text_ink_offsets
from gbdraw.api.options import LinearDiagramOptions
from gbdraw.api.request_render import plan_request
from gbdraw.api.requests import InMemoryRecordSource, LinearDiagramRequest, RecordInput
from gbdraw.features.objects import FeatureLocationPart, FeatureObject
from gbdraw.labels import linear as linear_labels
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


def _box_inside_drawn_arrow(
    label: dict, anchor_x: float, start: int, end: int, strand: str, arrow_length: float, margin_px: float = MARGIN_PX
) -> bool:
    """Whether the ink box widened by ``margin_px`` at each end at ``anchor_x`` lies inside the arrow the
    feature's path draws."""

    y_positions = (label["feature_top_y"], label["feature_middle_y"], label["feature_bottom_y"])
    path = create_arrow_path_linear(
        {"feat_start": start, "feat_end": end, "feat_strand": strand}, arrow_length, 20.0, strand, 1000, 1000.0,
        False, feature_y_positions=y_positions,
    )[1]
    polygon = _points(path)
    left, right, top, bottom = _ink(label["label_text"])
    middle_y = label["middle_y"]
    corners = [(anchor_x + x, middle_y + y) for x in (left - margin_px, right + margin_px) for y in (top, bottom)]
    return all(_inside_convex(polygon, corner) for corner in corners)


@pytest.mark.linear
@pytest.mark.parametrize("strand", ["positive", "negative"])
@pytest.mark.parametrize("margin_em", [0.25, 0.5])
def test_a_label_that_touches_the_head_at_the_middle_slides_toward_the_tail(
    strand: str, margin_em: float, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The ink keeps half the room a label needs inside its feature (Owner D-20) at each end: the one margin
    sets both whether a label may go inside and where it fits."""

    monkeypatch.setattr(linear_labels, "_EMBEDDED_LABEL_MARGIN_EM", margin_em)
    margin_px = margin_em / 2 * FONT_SIZE
    text, head = "ea59", 30.0
    left, right, top, bottom = _ink(text)
    share = max(-top, bottom) / 10.0
    width = right - left + 2 * margin_px
    length = int(width + head * share) + 3
    assert length < 2 * (max(right, -left) + margin_px + head * share) - 2, "the middle must not fit"
    start = 400
    end = start + length
    label = _labels({"f": _arrow("f", start, end, text, strand)}, head)[text]

    assert label["is_embedded"]
    middle = label["feature_anchor_x"]
    toward_tail = -1.0 if strand == "positive" else 1.0
    assert (label["middle_x"] - middle) * toward_tail > 0.5
    assert not _box_inside_drawn_arrow(label, middle, start, end, strand, head, margin_px)
    assert _box_inside_drawn_arrow(label, label["middle_x"], start, end, strand, head, margin_px)
    # The nearest x that fits: half a px nearer the middle touches the outline.
    assert not _box_inside_drawn_arrow(label, label["middle_x"] - 0.5 * toward_tail, start, end, strand, head, margin_px)


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


# A 10 kb record, one px per 5 bp at the default width: a long arrow labelled "C" and a short unlabelled one.
_LONG = (1000, 4000, 1, "C")
_OVER_ITS_MIDDLE = (2400, 2620, 1, "")


def _figure(
    features: list[tuple[int, int, int, str]], overrides: dict[str, object] | None = None
) -> tuple[list[str], str, float]:
    """The ids of the features in the order the figure draws them, and the feature id and x of label "C"."""

    record = SeqRecord(Seq("AAAACCGTGG" * 1000), id="r", annotations={"topology": "linear", "molecule_type": "DNA"})
    record.features = [
        SeqFeature(SimpleLocation(start, end, strand=strand), type="CDS",
                   qualifiers={"product": [product], "locus_tag": [f"t{index}"]})
        for index, (start, end, strand, product) in enumerate(features)
    ]
    config_overrides = {"labels.linear.scope": "all", "canvas.show_gc": False, "canvas.show_skew": False,
                        **(overrides or {})}
    request = LinearDiagramRequest(
        records=(RecordInput(InMemoryRecordSource(record)),),
        options=LinearDiagramOptions(config_overrides=config_overrides),
    )
    root = ET.fromstring(plan_request(request).build().drawing.tostring())
    drawn = [node.get("data-gbdraw-feature-id") for node in root.iter()
             if node.tag.endswith("}path") and node.get("data-gbdraw-feature-part") == "block"]
    (label,) = [node for node in root.iter() if node.tag.endswith("}text") and node.get("data-label-feature-id")]
    assert "".join(label.itertext()) == "C"
    return drawn, label.get("data-label-feature-id"), float(label.get("x"))


@pytest.mark.linear
def test_a_label_keeps_off_exactly_the_features_the_figure_draws_after_its_own() -> None:
    """The "later" features a label keeps off are those the figure draws over its feature: the label code and
    the drawer read one order (``feature_dict``). A drawer that changed its order fails here."""

    _drawn, _feature_id, alone_x = _figure([_LONG])
    outcomes = set()
    for features in ([_LONG, _OVER_ITS_MIDDLE], [_OVER_ITS_MIDDLE, _LONG]):
        drawn, feature_id, label_x = _figure(features)
        (other_id,) = [drawn_id for drawn_id in drawn if drawn_id != feature_id]
        drawn_over = drawn.index(other_id) > drawn.index(feature_id)
        assert (label_x != alone_x) == drawn_over, features
        outcomes.add(drawn_over)
    assert outcomes == {True, False}, "both draw orders must be met"


@pytest.mark.linear
@pytest.mark.parametrize(
    ("overrides", "elsewhere"),
    [
        pytest.param({"canvas.strandedness": True}, (2400, 2620, -1, ""), id="another lane"),
        pytest.param({"canvas.linear.wrap_length": 5000}, (7400, 7620, 1, ""), id="another wrap row"),
    ],
)
def test_a_later_part_in_another_lane_or_wrap_row_does_not_move_a_label(
    overrides: dict[str, object], elsewhere: tuple[int, int, int, str]
) -> None:
    """A later feature at the label's x moves it only when it is drawn in the label's lane and wrap row; the
    other strand's lane, or the next row of a wrapped record at the same x, leaves it at the middle."""

    _drawn, _feature_id, alone_x = _figure([_LONG], overrides)
    assert _figure([_LONG, _OVER_ITS_MIDDLE], overrides)[2] != alone_x
    assert _figure([_LONG, elsewhere], overrides)[2] == alone_x
