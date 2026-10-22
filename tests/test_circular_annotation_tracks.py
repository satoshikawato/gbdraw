import json
import math
import random
import re
import warnings
from pathlib import Path
from typing import NamedTuple

import pandas as pd
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
import pytest
from xml.etree import ElementTree

from gbdraw.api import (
    AnnotationOptions,
    AnnotationSet,
    CircularDiagramOptions,
    CircularRequestTrackOptions,
    CoordinateSpan,
    RegionAnnotation,
    RegionAnnotationStyle,
    parse_record_selector,
)
import gbdraw.api.render as api_render
from gbdraw.exceptions import LabelLeadersCrossWarning
from gbdraw.api.diagram import build_circular_diagram, build_circular_multi_diagram
from gbdraw.core.text import dominant_baseline_shift_em, get_font_vertical_metrics, text_ink
from gbdraw.tracks import (
    CircularTrackSlot,
    normalize_circular_track_slots,
    parse_circular_track_slot,
)
from gbdraw.web_support.request_render import render_embedded_canonical_web_request


def test_circular_annotation_renderer_parses_width() -> None:
    slot = parse_circular_track_slot("regions:annotations@set_id=s,w=28px")
    assert slot.renderer == "annotations"
    assert slot.width is not None and slot.width.value == 28


def test_circular_annotation_overlay_rejects_unknown_anchor_in_complete_list() -> None:
    annotation = CircularTrackSlot(
        "a",
        "annotations",
        side="overlay",
        z=1,
        params={"set_id": "s", "anchor_slot": "missing", "layer": "foreground"},
    )

    with pytest.raises(ValueError, match="unknown anchor_slot='missing'"):
        normalize_circular_track_slots([annotation])


def test_circular_annotation_overlay_rejects_annotation_anchor() -> None:
    anchor = CircularTrackSlot(
        "base_notes",
        "annotations",
        side="outside",
        z=0,
        params={"set_id": "base"},
    )
    annotation = CircularTrackSlot(
        "overlay_notes",
        "annotations",
        side="overlay",
        z=1,
        params={
            "set_id": "overlay",
            "anchor_slot": "base_notes",
            "layer": "foreground",
        },
    )

    with pytest.raises(ValueError, match="cannot use annotation slot 'base_notes' as anchor"):
        normalize_circular_track_slots([anchor, annotation])


def test_circular_annotation_overlay_rejects_spacer_anchor() -> None:
    spacer = CircularTrackSlot(
        "space",
        "spacer",
        side="inside",
    )
    annotation = CircularTrackSlot(
        "overlay_notes",
        "annotations",
        side="overlay",
        z=1,
        params={
            "set_id": "overlay",
            "anchor_slot": "space",
            "layer": "foreground",
        },
    )

    with pytest.raises(ValueError, match="anchor 'space' has no drawable band"):
        normalize_circular_track_slots([spacer, annotation])


def test_circular_annotation_overlay_allows_split_feature_anchor() -> None:
    feature = CircularTrackSlot("features", "features", side="overlay", z=0)
    annotation = CircularTrackSlot(
        "overlay_notes",
        "annotations",
        side="overlay",
        z=1,
        params={
            "set_id": "overlay",
            "anchor_slot": "features",
            "layer": "foreground",
        },
    )

    normalized = normalize_circular_track_slots([feature, annotation])

    assert normalized[0].side == "overlay"
    assert normalized[0].reserve is True
    assert normalized[1].annotation is not None


def test_circular_annotation_cover_anchor_uses_feature_band_width() -> None:
    record = SeqRecord(Seq("A" * 1000), id="r1", name="r1")
    annotation_options = AnnotationOptions(
        sets=(
            AnnotationSet(
                "regions",
                (
                    RegionAnnotation(
                        "band",
                        CoordinateSpan(None, 100, 300),
                        mark="band",
                    ),
                ),
            ),
        )
    )

    def resolved_widths(cover_anchor: bool) -> tuple[float, float]:
        drawing = build_circular_diagram(
            record,
            options=CircularDiagramOptions(
                annotations=annotation_options,
                tracks=CircularRequestTrackOptions(
                    circular_track_slots=(
                        "features:features@side=inside,lane_direction=inside,r=0.75",
                        (
                            "regions:annotations@side=overlay,w=14px,set_id=regions,"
                            "anchor_slot=features,layer=foreground,padding_px=0,z=1,"
                            f"cover_anchor={'true' if cover_anchor else 'false'}"
                        ),
                    ),
                ),
            ),
        )
        slots = {
            slot["slotId"]: slot
            for slot in drawing._gbdraw_track_slot_geometry["records"][0]["slots"]
        }
        return float(slots["features"]["widthPx"]), float(slots["regions"]["widthPx"])

    feature_width, uncovered_width = resolved_widths(False)
    covered_feature_width, covered_width = resolved_widths(True)
    assert uncovered_width == pytest.approx(14.0)
    assert covered_feature_width == pytest.approx(feature_width)
    assert covered_width == pytest.approx(feature_width)


def test_circular_overlay_annotation_can_share_an_on_axis_feature_anchor() -> None:
    record = SeqRecord(Seq("A" * 1000), id="r1", name="r1")
    annotation_options = AnnotationOptions(
        sets=(
            AnnotationSet(
                "regions",
                (
                    RegionAnnotation(
                        "band",
                        CoordinateSpan(None, 100, 300),
                        mark="band",
                    ),
                ),
            ),
        )
    )
    drawing = build_circular_diagram(
        record,
        options=CircularDiagramOptions(
            annotations=annotation_options,
            tracks=CircularRequestTrackOptions(
                circular_track_slots=(
                    "features:features@side=overlay,lane_direction=split",
                    (
                        "regions:annotations@side=overlay,w=14px,set_id=regions,"
                        "anchor_slot=features,layer=foreground,padding_px=0,"
                        "cover_anchor=true,z=1"
                    ),
                ),
            ),
        ),
    )
    slots = {
        slot["slotId"]: slot
        for slot in drawing._gbdraw_track_slot_geometry["records"][0]["slots"]
    }
    assert slots["regions"]["radiusFactor"] == pytest.approx(
        slots["features"]["radiusFactor"]
    )
    assert slots["regions"]["widthPx"] == pytest.approx(
        slots["features"]["widthPx"]
    )


def _span_annotations(record_count: int = 1) -> AnnotationOptions:
    return AnnotationOptions(
        sets=(
            AnnotationSet(
                "regions",
                tuple(
                    RegionAnnotation(
                        f"span{index}",
                        CoordinateSpan(parse_record_selector(f"#{index}"), 100, 300),
                        mark="highlight",
                    )
                    for index in range(1, record_count + 1)
                ),
            ),
        )
    )


def _span_slot(spec: str = "anchor_slot=features,to_slot=gc,layer=underlay,z=-1") -> str:
    return f"regions:annotations@side=overlay,set_id=regions,padding_px=0,{spec}"


def _circular_bands(drawing, record_index: int = 0) -> dict[str, tuple[float, float]]:
    """Each slot's resolved band, and the radii of the drawn ``span`` sector."""

    record = drawing._gbdraw_track_slot_geometry["records"][record_index]
    bands = {
        slot["slotId"]: (
            float(slot["radiusFactor"]) * float(record["axisRadiusPx"]) - 0.5 * float(slot["widthPx"]),
            float(slot["radiusFactor"]) * float(record["axisRadiusPx"]) + 0.5 * float(slot["widthPx"]),
        )
        for slot in record["slots"]
    }
    root = ElementTree.fromstring(drawing.tostring())
    sector = next(
        path
        for element in root.iter()
        if element.attrib.get("data-gbdraw-annotation-id") == f"span{record_index + 1}"
        for path in element
    )
    radii = sorted(float(value) for value in re.findall(r"A ([0-9.]+),", sector.attrib["d"]))
    bands["drawn"] = (radii[0], radii[-1])
    return bands


def _build_circular_span(
    *slots: str | CircularTrackSlot,
    features: tuple[tuple[int, int], ...] = (),
    strands: tuple[int, ...] = (1,),
    **overrides: bool,
):
    record = SeqRecord(Seq("ACGT" * 2500), id="r1", name="r1")
    record.features = [
        SeqFeature(FeatureLocation(start, end, strand=strand), type="CDS")
        for strand in strands
        for start, end in features
    ]
    return build_circular_diagram(
        record,
        options=CircularDiagramOptions(
            annotations=_span_annotations(),
            config_overrides=overrides,
            tracks=CircularRequestTrackOptions(circular_track_slots=slots),
        ),
    )


_SPAN_STACK = ("features:features", "gc:dinucleotide_content", "skew:dinucleotide_skew")


def test_circular_annotation_to_slot_spans_the_hull_of_both_slots() -> None:
    bands = _circular_bands(_build_circular_span(*_SPAN_STACK, _span_slot()))

    hull = (bands["gc"][0], bands["features"][1])
    assert bands["gc"][1] < bands["features"][0]
    assert bands["regions"] == pytest.approx(hull)
    assert bands["drawn"] == pytest.approx(hull, abs=1e-3)

    across = _circular_bands(
        _build_circular_span(
            *_SPAN_STACK, _span_slot("anchor_slot=features,to_slot=skew,layer=underlay,z=-1")
        )
    )
    assert across["drawn"] == pytest.approx((across["skew"][0], across["features"][1]), abs=1e-3)


@pytest.mark.parametrize(
    "spec",
    [
        "anchor_slot=gc,to_slot=features,layer=underlay,z=-1",
        "anchor_slot=features,to_slot=gc,layer=underlay,z=-1,cover_anchor=true",
        "anchor_slot=features,to_slot=gc,layer=underlay,z=-1,w=3px",
        "anchor_slot=features,to_slot=gc,layer=underlay,z=-1,w=200px,r=0.3",
        "anchor_slot=features,to_slot=gc,layer=foreground,z=1",
    ],
)
def test_circular_annotation_to_slot_ignores_order_cover_anchor_and_own_width(spec: str) -> None:
    reference = _circular_bands(_build_circular_span(*_SPAN_STACK, _span_slot()))
    bands = _circular_bands(_build_circular_span(*_SPAN_STACK, _span_slot(spec)))

    assert bands["drawn"] == reference["drawn"]
    assert bands["regions"] == reference["regions"]


def test_circular_annotation_to_slot_padding_shrinks_both_edges() -> None:
    bands = _circular_bands(
        _build_circular_span(
            *_SPAN_STACK,
            "regions:annotations@side=overlay,set_id=regions,padding_px=4,"
            "anchor_slot=features,to_slot=gc,layer=underlay,z=-1",
        )
    )

    assert bands["drawn"] == pytest.approx(
        (bands["gc"][0] + 4.0, bands["features"][1] - 4.0), abs=1e-3
    )


@pytest.mark.parametrize(
    ("spec", "message"),
    [
        ("anchor_slot=features,to_slot=missing,layer=underlay,z=-1", "references unknown to_slot='missing'"),
        ("anchor_slot=features,to_slot=off,layer=underlay,z=-1", "references unknown to_slot='off'"),
        ("anchor_slot=features,to_slot=space,layer=underlay,z=-1", "to_slot 'space' has no drawable band"),
        ("anchor_slot=features,to_slot=ticks,layer=underlay,z=-1", "to_slot 'ticks' has no drawable band"),
        ("anchor_slot=features,to_slot=notes,layer=underlay,z=-1", "cannot use annotation slot 'notes' as to_slot"),
        ("anchor_slot=features,to_slot=gc,layer=underlay,z=0", "must have z less than anchor 'features'"),
        ("anchor_slot=features,to_slot=low,layer=underlay,z=-1", "must have z less than to_slot 'low'"),
        ("anchor_slot=low,to_slot=features,layer=underlay,z=-1", "must have z less than anchor 'low'"),
        ("anchor_slot=features,to_slot=high,layer=foreground,z=1", "must have z greater than to_slot 'high'"),
        ("anchor_slot=features,to_slot=,layer=underlay,z=-1", "to_slot cannot be empty"),
        ("anchor_slot=features,to_slot=features,layer=underlay,z=-1", "to_slot must name a slot other than its anchor_slot 'features'"),
    ],
)
def test_circular_annotation_to_slot_reference_errors(spec: str, message: str) -> None:
    slots = [
        parse_circular_track_slot(text)
        for text in (
            "features:features",
            "gc:dinucleotide_content",
            "off:dinucleotide_skew@enabled=false",
            "space:spacer@w=8px",
            "ticks:ticks",
            "notes:annotations@set_id=regions",
            "low:depth@z=-2",
            "high:dinucleotide_skew@nt=AT,z=2",
        )
    ]

    with pytest.raises(ValueError, match=f"annotation .*slot 'regions'.*{message}|^{message}"):
        normalize_circular_track_slots([*slots, parse_circular_track_slot(_span_slot(spec))])


def test_circular_annotation_to_slot_requires_overlay_side() -> None:
    slot = CircularTrackSlot(
        "regions", "annotations", side="outside", params={"set_id": "regions", "to_slot": "gc"}
    )

    with pytest.raises(ValueError, match="annotation slot 'regions' uses to_slot without side=overlay"):
        normalize_circular_track_slots([parse_circular_track_slot("gc:dinucleotide_content"), slot])


@pytest.mark.parametrize("spanning", [False, True])
def test_circular_overlay_annotation_reads_slot_references_in_any_key_case(spanning: bool) -> None:
    def bands(anchor_key: str, to_key: str) -> dict[str, tuple[float, float]]:
        params: dict[str, object] = {
            "set_id": "regions",
            "layer": "underlay",
            "padding_px": 0,
            anchor_key: "gc",
        }
        if spanning:
            params[to_key] = "skew"
        overlay = CircularTrackSlot("regions", "annotations", side="overlay", z=-1, params=params)
        return _circular_bands(
            _build_circular_span(*(parse_circular_track_slot(text) for text in _SPAN_STACK), overlay)
        )

    assert bands("Anchor_Slot", "To_Slot") == bands("anchor_slot", "to_slot")


@pytest.mark.parametrize(
    ("features_slot", "gc_slot", "stranded"),
    [
        pytest.param("features:features@side=overlay,lane_direction=split", "gc:dinucleotide_content", True, id="split"),
        pytest.param("features:features@lane_direction=outside", "gc:dinucleotide_content", False, id="outside"),
        pytest.param("features:features", "gc:dinucleotide_content", False, id="inside"),
        pytest.param(
            "features:features@side=overlay,lane_direction=split",
            "gc:dinucleotide_content@side=outside",
            True,
            id="split-to-an-outer-slot",
        ),
    ],
)
def test_circular_annotation_to_slot_covers_every_feature_lane(
    features_slot: str, gc_slot: str, stranded: bool
) -> None:
    def build(overlay: str):
        return _build_circular_span(
            features_slot,
            gc_slot,
            overlay,
            features=((100, 900), (500, 1500), (700, 1700)),
            strands=(1, -1) if stranded else (1,),
            **{"canvas.resolve_overlaps": True, "canvas.strandedness": stranded},
        )

    drawing = build(_span_slot())
    bands = _circular_bands(drawing)
    lane_radii = [
        float(value)
        for element in ElementTree.fromstring(drawing.tostring()).iter()
        if "data-gbdraw-feature-id" in element.attrib
        for value in re.findall(r"A ?([0-9.]+),", element.attrib.get("d", ""))
    ]
    lanes = (min(lane_radii), max(lane_radii))

    # Overlap resolution gave the slot more lanes than its first one.
    assert lanes[1] - lanes[0] > bands["features"][1] - bands["features"][0] + 1
    assert bands["drawn"] == pytest.approx(
        (min(lanes[0], bands["gc"][0]), max(lanes[1], bands["gc"][1])), abs=1e-3
    )
    # cover_anchor alone keeps the first lane of each strand, not every lane.
    covered = _circular_bands(
        build(_span_slot("anchor_slot=features,cover_anchor=true,layer=underlay,z=-1"))
    )["drawn"]
    assert lanes[0] <= covered[0] + 1e-3 and covered[1] <= lanes[1] + 1e-3
    assert covered[1] - covered[0] < lanes[1] - lanes[0] - 1


def test_circular_annotation_to_slot_spans_a_depth_slot_without_data_for_a_record() -> None:
    records = [
        SeqRecord(Seq("ACGT" * 600), id=f"rec{index}", name=f"rec{index}") for index in (1, 2)
    ]
    depth = pd.DataFrame(
        {
            "reference_name": ["rec1"] * 2400,
            "position": list(range(1, 2401)),
            "depth": [20.0] * 2400,
        }
    )
    drawing = build_circular_multi_diagram(
        records,
        options=CircularDiagramOptions(
            annotations=_span_annotations(record_count=2),
            depth_track_tables=[[depth], []],
            depth_window=100,
            depth_step=100,
            tracks=CircularRequestTrackOptions(
                circular_track_slots=(
                    "features:features",
                    "depth:depth@track_index=0",
                    _span_slot("anchor_slot=features,to_slot=depth,layer=underlay,z=-1"),
                ),
            ),
        ),
    )

    # Circular keeps the Depth slot's band on a record that has no depth data,
    # so every record gets the same span.
    for record_index in (0, 1):
        bands = _circular_bands(drawing, record_index)
        assert bands["drawn"] == pytest.approx((bands["depth"][0], bands["features"][1]), abs=1e-3)
    assert drawing.tostring().count('data-gbdraw-slot-renderer="depth"') == 1


def test_circular_origin_annotation_renders_two_safe_paths() -> None:
    record = SeqRecord(Seq("A" * 1000), id="r1", name="r1")
    record.annotations["topology"] = "circular"
    annotation = RegionAnnotation(
        "origin",
        CoordinateSpan(None, 900, 100, wraps_origin=True),
        label="Origin",
        mark="bracket",
        style=RegionAnnotationStyle(stroke="#225588"),
    )
    drawing = build_circular_diagram(
        record,
        options=CircularDiagramOptions(
            annotations=AnnotationOptions(sets=(AnnotationSet("regions", (annotation,)),))
        ),
    )
    svg = drawing.tostring()
    assert 'data-gbdraw-annotation-id="origin"' in svg
    assert ">Origin<" in svg
    assert svg.count("A ") >= 2


def test_circular_highlight_automatically_renders_behind_features() -> None:
    record = SeqRecord(Seq("A" * 1000), id="r1", name="r1")
    annotation = RegionAnnotation(
        "highlighted",
        CoordinateSpan(None, 100, 300),
        label="Highlighted region",
        mark="highlight",
    )
    drawing = build_circular_diagram(
        record,
        options=CircularDiagramOptions(
            annotations=AnnotationOptions(
                sets=(AnnotationSet("regions", (annotation,)),)
            )
        ),
    )
    svg = drawing.tostring()
    slots = {
        slot["slotId"]: slot
        for slot in drawing._gbdraw_track_slot_geometry["records"][0]["slots"]
    }

    assert slots["annotations_1"]["side"] == "overlay"
    assert slots["annotations_1"]["radiusFactor"] == slots["features"]["radiusFactor"]
    assert slots["annotations_1"]["widthPx"] == pytest.approx(
        slots["features"]["widthPx"]
    )
    assert 'data-gbdraw-annotation-mark="highlight"' in svg
    assert 'fill="#94a3b8"' in svg
    root = ElementTree.fromstring(svg)
    elements = list(root.iter())
    highlight_group = next(
        element
        for element in elements
        if element.attrib.get("data-gbdraw-annotation-id") == "highlighted"
    )
    record_group = next(
        element
        for element in root
        if element.attrib.get("data-gbdraw-record-id") == "r1"
    )
    assert elements.index(highlight_group) < elements.index(record_group)


def test_circular_automatic_annotation_slots_respect_hidden_scale() -> None:
    record = SeqRecord(Seq("A" * 1000), id="r1", name="r1")
    annotation = RegionAnnotation(
        "highlighted",
        CoordinateSpan(None, 100, 300),
        mark="highlight",
    )

    drawing = build_circular_diagram(
        record,
        options=CircularDiagramOptions(
            config_overrides={"objects.scale.show": False},
            annotations=AnnotationOptions(
                sets=(AnnotationSet("regions", (annotation,)),)
            ),
        ),
    )
    svg = drawing.tostring()
    renderers = {
        slot["renderer"]
        for slot in drawing._gbdraw_track_slot_geometry["records"][0]["slots"]
    }

    assert "ticks" not in renderers
    assert 'data-gbdraw-annotation-id="highlighted"' in svg
    assert 'id="Axis"' in svg


def test_circular_multi_record_annotations_bind_by_record_index() -> None:
    records = [
        SeqRecord(Seq("A" * 1000), id="duplicate", name="duplicate"),
        SeqRecord(Seq("A" * 1200), id="duplicate", name="duplicate"),
    ]
    annotations = AnnotationSet(
        "regions",
        (
            RegionAnnotation(
                "first",
                CoordinateSpan(parse_record_selector("#1"), 100, 200),
            ),
            RegionAnnotation(
                "second",
                CoordinateSpan(parse_record_selector("#2"), 300, 500),
            ),
        ),
    )
    options = CircularDiagramOptions(
        annotations=AnnotationOptions(sets=(annotations,))
    )

    multi_svg = build_circular_multi_diagram(records, options=options).tostring()
    single_svg = build_circular_multi_diagram(records[:1], options=CircularDiagramOptions(
        annotations=AnnotationOptions(
            sets=(AnnotationSet("regions", (annotations.annotations[0],)),)
        )
    )).tostring()

    assert multi_svg.count('data-gbdraw-annotation-id="') == 2
    assert 'data-gbdraw-annotation-id="first"' in single_svg


def test_circular_multi_record_curved_labels_follow_their_own_records_paths() -> None:
    records = [
        SeqRecord(Seq("A" * 1000), id="one", name="one"),
        SeqRecord(Seq("A" * 1600), id="two", name="two"),
    ]
    # A labelled region on each record, each in its record's copy of the slot group.
    regions = tuple(
        RegionAnnotation(f"r{record[1]}", CoordinateSpan(parse_record_selector(record), 100, 300), label="Region")
        for record in ("#1", "#2")
    )
    options = CircularDiagramOptions(
        annotations=AnnotationOptions(sets=(AnnotationSet("regions", regions),)),
        tracks=CircularRequestTrackOptions(
            circular_track_slots=("features:features", "regions:annotations@side=outside,set_id=regions,w=20px")
        ),
    )

    root = ElementTree.fromstring(build_circular_multi_diagram(records, options=options).tostring())
    hrefs = [
        element.attrib["{http://www.w3.org/1999/xlink}href"].removeprefix("#")
        for element in root.iter()
        if element.tag.endswith("}textPath") and element.text == "Region"
    ]
    paths = {element.attrib["id"]: element.attrib["d"] for element in root.iter() if element.attrib.get("id") in hrefs}
    # Each record's label follows its own path, at its own region's middle (OV-1210).
    assert len(set(hrefs)) == 2 and len(paths) == 2 and len(set(paths.values())) == 2


_PAIR = re.compile(r"(-?[0-9.]+(?:e-?[0-9]+)?),(-?[0-9.]+(?:e-?[0-9]+)?)")


def _mark_radii(group) -> tuple[float, float]:
    """The radial band an annotation's drawn mark paints, stroke included."""

    low, high = math.inf, -math.inf
    for element in group:
        if element.tag.rsplit("}", 1)[-1] != "path":
            continue
        stroked = element.attrib.get("stroke", "none") != "none"
        half_stroke = 0.5 * float(element.attrib.get("stroke-width", 0.0)) if stroked else 0.0
        d = element.attrib["d"]
        # A path with an arc is half its stroke thick each way; a straight butt-capped radial cap ends at its ends.
        pad = half_stroke if "A" in d else 0.0
        # Skip an arc's radii and flags ("A rx,ry rotation large sweep"); keep its end point.
        points = _PAIR.findall(re.sub(r"A \S+ \S+ \S+ \S+ ", "A ", d))
        radii = [math.hypot(float(x), float(y)) for x, y in points]
        low, high = min(low, min(radii) - pad), max(high, max(radii) + pad)
    return low, high


class _DrawnLabel(NamedTuple):
    """A drawn annotation label's painted letters, in the annotation track's frame: convex ``parts``, how near the
    circle centre and how far from it they come, and the label's middle (``x``, ``y``)."""

    parts: list[list[tuple[float, float]]]
    near: float
    far: float
    x: float
    y: float


_ARC_PATH = re.compile(r"M (\S+),(\S+) A (\S+),(\S+) 0 ([01]) ([01]) (\S+),(\S+)")


def _drawn_label(text, root) -> _DrawnLabel:
    """A drawn label read back from the SVG: straight text by its anchor and turn, curved text by its text path."""

    family = text.attrib["font-family"].split(",")[0].strip("'\" ")
    size = float(text.attrib["font-size"])
    assert text.attrib["text-anchor"] == "middle"
    text_path = next((child for child in text if child.tag.endswith("}textPath")), None)
    if text_path is None:
        corners = _label_ink_corners(text)
        x, y = float(text.attrib["x"]), float(text.attrib["y"])
        return _DrawnLabel([corners], _nearest_distance(corners), max(math.hypot(px, py) for px, py in corners), x, y)
    href = text_path.attrib.get("{http://www.w3.org/1999/xlink}href") or text_path.attrib["href"]
    (path,) = [element for element in root.iter() if element.attrib.get("id") == href.removeprefix("#")]
    match = _ARC_PATH.fullmatch(path.attrib["d"])
    assert match is not None, path.attrib["d"]
    x1, y1, radius, radius_y, large, sweep, x2, y2 = (float(value) for value in match.groups())
    assert radius == radius_y and math.hypot(x1, y1) == pytest.approx(radius, rel=1e-5)
    flipped = sweep == 0.0
    start, end = math.atan2(y1, x1), math.atan2(y2, x2)
    span = ((start - end) if flipped else (end - start)) % (2.0 * math.pi)
    assert (span > math.pi) == (large == 1.0)
    middle = start - 0.5 * span if flipped else start + 0.5 * span
    assert text_path.attrib["startOffset"] == "50%"
    ink = text_ink(text_path.text or "", family)
    assert ink is not None
    # The whole advance lies on the path, so every glyph is drawn.
    assert radius * span >= ink.advance * size
    # The middle of the advance sits at the middle of the path, on the alphabetic baseline; tops face out unless the
    # text reads anticlockwise.
    left, right = (ink.left - 0.5 * ink.advance) * size / radius, (ink.right - 0.5 * ink.advance) * size / radius
    if flipped:
        near, far = radius - ink.top * size, radius - ink.bottom * size
        low, high = middle - right, middle - left
    else:
        near, far = radius + ink.bottom * size, radius + ink.top * size
        low, high = middle + left, middle + right
    steps = max(1, math.ceil(math.degrees(high - low) / 0.5))
    angles = [low + (high - low) * index / steps for index in range(steps + 1)]
    parts = [
        [(near * math.cos(a), near * math.sin(a)), (far * math.cos(a), far * math.sin(a)),
         (far * math.cos(b), far * math.sin(b)), (near * math.cos(b), near * math.sin(b))]
        for a, b in zip(angles, angles[1:])
    ]
    return _DrawnLabel(parts, near, far, radius * math.cos(middle), radius * math.sin(middle))


def _label_ink_corners(text) -> list[tuple[float, float]]:
    """The four corners of a drawn label's painted letters, in the annotation track's frame."""

    family = text.attrib["font-family"].split(",")[0].strip("'\" ")
    size = float(text.attrib["font-size"])
    ink = text_ink(text.text or "", family)
    metrics = get_font_vertical_metrics(family)
    assert ink is not None and metrics is not None
    shift = dominant_baseline_shift_em(text.attrib["dominant-baseline"], metrics)
    assert shift is not None
    x, y = float(text.attrib["x"]), float(text.attrib["y"])
    origin = -0.5 * ink.advance * size
    left, right = origin + ink.left * size, origin + ink.right * size
    top, bottom = (shift - ink.top) * size, (shift - ink.bottom) * size
    match = re.search(r"rotate\(([-0-9.e]+),([-0-9.e]+),([-0-9.e]+)\)", text.attrib.get("transform", ""))
    angle = math.radians(float(match.group(1))) if match else 0.0
    if match:
        assert (float(match.group(2)), float(match.group(3))) == pytest.approx((x, y))
    cos, sin = math.cos(angle), math.sin(angle)
    return [
        (x + cx * cos - cy * sin, y + cx * sin + cy * cos)
        for cx, cy in ((left, top), (right, top), (right, bottom), (left, bottom))
    ]


def _nearest_distance(corners: list[tuple[float, float]]) -> float:
    """The distance from the circle centre to the nearest point of a convex quadrilateral that excludes it."""

    best = math.inf
    for (x1, y1), (x2, y2) in zip(corners, corners[1:] + corners[:1]):
        dx, dy = x2 - x1, y2 - y1
        t = max(0.0, min(1.0, -(x1 * dx + y1 * dy) / (dx * dx + dy * dy)))
        best = min(best, math.hypot(x1 + t * dx, y1 + t * dy))
    return best


def _circular_label_gaps(drawing) -> dict[str, tuple[str, float]]:
    """Each drawn annotation label: the side of its mark it is on and the gap from the mark's paint to its ink."""

    root = ElementTree.fromstring(drawing.tostring())
    gaps: dict[str, tuple[str, float]] = {}
    for group in root.iter():
        if not group.attrib.get("data-gbdraw-annotation-label"):
            continue
        inner, outer = _mark_radii(group)
        for text in (element for element in group if element.tag.rsplit("}", 1)[-1] == "text"):
            drawn = _drawn_label(text, root)
            if math.hypot(drawn.x, drawn.y) > 0.5 * (inner + outer):
                gaps[group.attrib["data-gbdraw-annotation-id"]] = ("outside", drawn.near - outer)
            else:
                gaps[group.attrib["data-gbdraw-annotation-id"]] = ("inside", inner - drawn.far)
    assert gaps, "no annotation labels drawn"
    return gaps


@pytest.mark.parametrize("side", ["outside", "inside"])
@pytest.mark.parametrize("orientation", ["auto", "tangent", "arc", "radial", "horizontal"])
@pytest.mark.parametrize(
    ("mark", "line_cap"),
    [("highlight", "tick"), ("band", "tick"), ("bracket", "tick"), ("bracket", "none"), ("line", "tick")],
)
def test_circular_annotation_label_ink_keeps_label_offset_clear_of_its_mark(
    mark: str, line_cap: str, orientation: str, side: str
) -> None:
    style = RegionAnnotationStyle(
        stroke_width=1.5, label_offset=4.0, label_orientation=orientation, line_cap=line_cap
    )
    record = SeqRecord(Seq("ACGT" * 2500), id="r1", name="r1")
    record.annotations["topology"] = "circular"
    # Regions around the circle, so tangent labels turn over and horizontal ones meet the ring at every slant.
    regions = tuple(
        RegionAnnotation(f"r{index}", CoordinateSpan(None, start, start + 700), label=text, mark=mark, style=style)
        for index, (start, text) in enumerate(
            ((300, "Region gyp"), (1900, "IRb"), (3500, "rRNA genes"), (5600, "Wide region label"), (8000, "LSC"))
        )
    )
    drawing = build_circular_diagram(
        record,
        options=CircularDiagramOptions(
            annotations=AnnotationOptions(sets=(AnnotationSet("regions", regions),)),
            config_overrides={"labels.circular.scope": "none"},
            tracks=CircularRequestTrackOptions(
                circular_track_slots=(
                    "features:features",
                    f"regions:annotations@side={side},set_id=regions,w=20px",
                )
            ),
        ),
    )

    gaps = _circular_label_gaps(drawing)
    assert set(gaps) == {region.id for region in regions}
    for annotation_id, (drawn_side, gap) in gaps.items():
        # The painted letters keep label_offset clear of the mark's paint (OV-849), on the slot's side.
        assert drawn_side == side, annotation_id
        # Mark paths carry six significant digits (``:g``), about 1e-3 px at these radii.
        assert gap == pytest.approx(4.0, abs=1e-3), annotation_id


@pytest.mark.parametrize("side", ["outside", "inside"])
@pytest.mark.parametrize("orientation", ["auto", "arc", "tangent"])
def test_circular_annotation_labels_curve_along_the_circle_and_stay_upright(orientation: str, side: str) -> None:
    style = RegionAnnotationStyle(label_offset=4.0, label_orientation=orientation)
    record = SeqRecord(Seq("ACGT" * 3000), id="r1", name="r1")
    record.annotations["topology"] = "circular"
    # Labels at 12, 6 o'clock and just either side of 3 and 9 o'clock; a long label on a short region.
    regions = tuple(
        RegionAnnotation(f"r{index}", CoordinateSpan(None, centre - 100, centre + 100), label=text, mark="bracket", style=style)
        for index, (centre, text) in enumerate(
            ((150, "A"), (2850, "IRb"), (3150, "LSC"), (6000, "Mitochondrial control region"), (8850, "SSC"), (9150, "IRa"))
        )
    )
    drawing = build_circular_diagram(
        record,
        options=CircularDiagramOptions(
            annotations=AnnotationOptions(sets=(AnnotationSet("regions", regions),)),
            config_overrides={"labels.circular.scope": "none"},
            tracks=CircularRequestTrackOptions(
                circular_track_slots=("features:features", f"regions:annotations@side={side},set_id=regions,w=20px")
            ),
        ),
    )

    root = ElementTree.fromstring(drawing.tostring())
    texts = {
        group.attrib["data-gbdraw-annotation-id"]: next(child for child in group if child.tag.endswith("}text"))
        for group in root.iter()
        if group.attrib.get("data-gbdraw-annotation-label")
    }
    assert set(texts) == {region.id for region in regions}
    for region in regions:
        text = texts[region.id]
        text_path = next((child for child in text if child.tag.endswith("}textPath")), None)
        if orientation == "tangent":
            # Tangent labels stay straight text turned to the tangent (OV-1210).
            assert text_path is None and "rotate(" in text.attrib["transform"]
            continue
        # Auto and arc labels follow the circle (OV-1210), the whole text drawn, centred over the region's middle.
        assert text_path is not None and _text_of(text) == region.label
        drawn = _drawn_label(text, root)
        middle = (region.target.start + region.target.end) / 2.0
        angle = math.degrees(math.atan2(drawn.y, drawn.x))
        assert (angle - (360.0 * middle / 12000 - 90.0) + 180.0) % 360.0 - 180.0 == pytest.approx(0.0, abs=0.05)
        # Text in the lower half reads anticlockwise so that it stays upright, as embedded feature labels do.
        href = text_path.attrib["{http://www.w3.org/1999/xlink}href"].removeprefix("#")
        (path,) = [element for element in root.iter() if element.attrib.get("id") == href]
        sweep = _ARC_PATH.fullmatch(path.attrib["d"]).group(6)
        assert (sweep == "0") == (0.0 < angle % 360.0 < 180.0), (region.label, angle)
    gaps = _circular_label_gaps(drawing)
    assert all(drawn_side == side and gap == pytest.approx(4.0, abs=1e-3) for drawn_side, gap in gaps.values()), gaps


def test_circular_region_label_longer_than_its_circle_is_cut_and_warns() -> None:
    from gbdraw.api.diagram import label_fit_report
    from gbdraw.labels.fit import label_fit_notice
    from gbdraw.web_support.label_fit import label_fit_report_payload

    record = SeqRecord(Seq("ACGT" * 2500), id="r1", name="r1")
    record.annotations["topology"] = "circular"
    text = " ".join(f"segment{index}" for index in range(40))
    regions = (
        RegionAnnotation("long", CoordinateSpan(None, 2000, 2500), label=text, style=RegionAnnotationStyle(label_font_size=14)),
        RegionAnnotation("short", CoordinateSpan(None, 7000, 7500), label="IRb"),
    )
    with pytest.warns(LabelLeadersCrossWarning, match=r"1 region label cut"):
        drawing = build_circular_diagram(
            record,
            options=CircularDiagramOptions(
                annotations=AnnotationOptions(sets=(AnnotationSet("regions", regions),)),
                config_overrides={"labels.circular.scope": "none"},
                tracks=CircularRequestTrackOptions(
                    circular_track_slots=("features:features", "regions:annotations@side=inside,set_id=regions,r=0.5,w=40px")
                ),
            ),
        )

    root = ElementTree.fromstring(drawing.tostring())
    texts = {
        group.attrib["data-gbdraw-annotation-id"]: next(child for child in group if child.tag.endswith("}text"))
        for group in root.iter()
        if group.attrib.get("data-gbdraw-annotation-label")
    }
    # A label longer than its circle keeps its start, ends with an ellipsis, and lies whole on its path (LP-D87).
    shown = _text_of(texts["long"])
    assert shown.endswith("\u2026") and text.startswith(shown[:-1].rstrip()) and len(shown) < len(text)
    _drawn_label(texts["long"], root)
    assert _text_of(texts["short"]) == "IRb"
    # The cut is counted in the label fit warning, the one region label report the command line and the Web read;
    # it offers no fix.
    report = label_fit_report(drawing)
    assert report is not None and report.drawn and report.region_cuts == 1 and report.rows == ()
    assert label_fit_notice(report) == ("Region labels are cut to fit their circle (1).",)
    payload = label_fit_report_payload(report, "circular")
    assert payload["notice"] == ["Region labels are cut to fit their circle (1)."]
    ((record_payload),) = payload["records"]
    assert record_payload["labelKind"] == "circular_region" and record_payload["counts"] == "1 region label cut"
    assert record_payload["labelCount"] == 2 and record_payload["fontSize"] == 14.0


def test_circular_region_labels_count_an_earlier_curved_label_as_one_obstacle(monkeypatch: pytest.MonkeyPatch) -> None:
    import gbdraw.diagrams.circular.assemble as assemble
    from gbdraw.render.drawers.circular.annotations import _ObstacleIndex

    calls: list = []
    original = assemble.circular_annotation_labels

    def record_call(layout, **kwargs):
        labels = original(layout, **kwargs)
        calls.append((layout.set_id, kwargs["obstacles"], labels))
        return labels

    monkeypatch.setattr(assemble, "circular_annotation_labels", record_call)
    record = SeqRecord(Seq("ACGT" * 2500), id="r1", name="r1")
    record.annotations["topology"] = "circular"
    build_circular_diagram(
        record,
        options=CircularDiagramOptions(
            annotations=AnnotationOptions(sets=(
                AnnotationSet("a", (RegionAnnotation("r0", CoordinateSpan(None, 2000, 3000), label="An earlier region"),)),
                AnnotationSet("b", (RegionAnnotation("r1", CoordinateSpan(None, 7000, 8000), label="A later region"),)),
            )),
            config_overrides={"labels.circular.scope": "none"},
            tracks=CircularRequestTrackOptions(
                circular_track_slots=("features:features", "a:annotations@set_id=a", "b:annotations@set_id=b")
            ),
        ),
    )

    first = next(call for call in calls if call[0] == "a")
    later = next(call for call in calls if call[0] == "b")
    (earlier,) = first[2]
    assert len(earlier.parts) > 1
    # The slot hand-off gives the later slot the earlier curved label as one shape (F-4).
    assert _ObstacleIndex(first[1]).conflicts(earlier) == 0
    assert _ObstacleIndex(later[1]).conflicts(earlier) == 1


def test_circular_annotation_labels_of_the_tobacco_gallery_figure_keep_their_label_offset(
    tmp_path, monkeypatch: pytest.MonkeyPatch
) -> None:
    session = json.loads(
        (Path(__file__).parents[1] / "gbdraw" / "web" / "gallery" / "sessions" / "tobacco-chloroplast.gbdraw-session.json")
        .read_text(encoding="utf-8")
    )
    drawings: list = []
    original = api_render.write_svg_document

    def capture(drawing):
        drawings.append(drawing)
        return original(drawing)

    monkeypatch.setattr(api_render, "write_svg_document", capture)
    render_embedded_canonical_web_request(
        session["renderRequest"], resources=session["resources"], workspace=tmp_path / "workspace"
    )

    (drawing,) = drawings
    gaps = _circular_label_gaps(drawing)
    assert set(gaps) == {"lsc", "irb", "ssc", "ira"}
    for annotation_id, (side, gap) in gaps.items():
        # The figure's region lines draw their labels inside with label_offset 7.
        assert side == "inside", annotation_id
        assert gap == pytest.approx(7.0, abs=1e-3), annotation_id


def test_circular_annotation_label_clears_the_axis_circle_on_its_way() -> None:
    record = SeqRecord(Seq("ACGT" * 2500), id="r1", name="r1")
    record.features = [SeqFeature(FeatureLocation(100, 900, strand=1), type="CDS")]
    tint = RegionAnnotation("tint", CoordinateSpan(None, 300, 2000), label="rRNA genes", mark="highlight")
    # The automatic highlight covers the feature band, just inside the axis circle.
    drawing = build_circular_diagram(
        record,
        options=CircularDiagramOptions(annotations=AnnotationOptions(sets=(AnnotationSet("regions", (tint,)),))),
    )

    root = ElementTree.fromstring(drawing.tostring())
    axis = next(element for element in root.iter() if element.attrib.get("id") == "Axis")[0]
    axis_outer = float(axis.attrib["r"]) + 0.5 * float(axis.attrib["stroke-width"])
    text = next(
        element for element in root.iter() if element.tag.rsplit("}", 1)[-1] == "text" and _text_of(element) == "rRNA genes"
    )
    (side, mark_gap), = _circular_label_gaps(drawing).values()
    assert side == "outside"
    # The letters keep label_offset beyond the axis paint, so they keep more than that from the tint. A text path
    # carries six significant digits (``:g``), about 1e-3 px at this radius.
    assert _drawn_label(text, root).near - axis_outer == pytest.approx(4.0, abs=1e-3)
    assert mark_gap > 4.0


def _text_of(text) -> str:
    """The characters a ``<text>`` draws, in a ``<textPath>`` or not."""

    return "".join(text.itertext())


def _ink_radii(text, root) -> tuple[float, float]:
    """How near the circle centre and how far from it a drawn label's painted letters come."""

    drawn = _drawn_label(text, root)
    return drawn.near, drawn.far


@pytest.mark.parametrize("side", ["outside", "inside"])
@pytest.mark.parametrize("orientation", ["auto", "tangent", "radial", "horizontal"])
@pytest.mark.parametrize("pinned", [False, True], ids=["auto", "pinned-annotations"])
def test_circular_annotation_labels_keep_label_offset_from_automatic_rings(
    side: str, orientation: str, pinned: bool
) -> None:
    style = RegionAnnotationStyle(label_font_size=14, label_offset=4.0, label_orientation=orientation)
    record = SeqRecord(Seq("ACGT" * 2500), id="r1", name="r1")
    record.features = [SeqFeature(FeatureLocation(100, 900, strand=1), type="CDS")]
    regions = tuple(
        RegionAnnotation(f"r{index}", CoordinateSpan(None, start, start + 900), label=text, mark="bracket", style=style)
        for index, (start, text) in enumerate(((300, "Region gyp"), (3500, "rRNA genes"), (6500, "Wide region label")))
    )
    radius = ",r=0.8" if side == "inside" and pinned else ",r=1.2" if pinned else ""
    rings_beyond = (f"gc:dinucleotide_content@side={side}", f"skew:dinucleotide_skew@side={side}")
    annotations = f"regions:annotations@side={side},set_id=regions,w=20px{radius}"
    # Slots are listed from the outermost ring inward, so the GC rings lie on the labels' side.
    slots = (*rings_beyond, annotations) if side == "outside" else (annotations, *rings_beyond)
    drawing = build_circular_diagram(
        record,
        options=CircularDiagramOptions(
            annotations=AnnotationOptions(sets=(AnnotationSet("regions", regions),)),
            config_overrides={"labels.circular.scope": "none"},
            tracks=CircularRequestTrackOptions(circular_track_slots=("features:features", *slots)),
        ),
    )

    geometry = drawing._gbdraw_track_slot_geometry["records"][0]
    rings = {
        slot["slotId"]: (
            float(slot["radiusFactor"]) * float(geometry["axisRadiusPx"]) - 0.5 * float(slot["widthPx"]),
            float(slot["radiusFactor"]) * float(geometry["axisRadiusPx"]) + 0.5 * float(slot["widthPx"]),
        )
        for slot in geometry["slots"]
        if slot["slotId"] in {"gc", "skew"}
    }
    root = ElementTree.fromstring(drawing.tostring())
    texts = [element for element in root.iter() if element.tag.rsplit("}", 1)[-1] == "text" and _text_of(element) in {
        region.label for region in regions
    }]
    assert len(texts) == len(regions)
    for text in texts:
        near, far = _ink_radii(text, root)
        for slot_id, (inner, outer) in rings.items():
            # The automatic rings keep out of the labels' room: label_offset beyond the letters (OV-1192); a curved
            # label's text path carries six significant digits.
            gap = max(inner - far, near - outer)
            assert gap >= 4.0 - 1e-3, (_text_of(text), slot_id, gap)


def _feature_label_corners(text) -> list[tuple[float, float]]:
    """The four corners of a drawn (horizontal) feature label's painted letters."""

    family = text.attrib["font-family"].split(",")[0].strip("'\" ")
    size = float(text.attrib["font-size"])
    ink = text_ink(text.text or "", family)
    metrics = get_font_vertical_metrics(family)
    assert ink is not None and metrics is not None
    shift = dominant_baseline_shift_em(text.attrib.get("dominant-baseline", "alphabetic"), metrics)
    assert shift is not None and "transform" not in text.attrib
    x, y = float(text.attrib["x"]), float(text.attrib["y"])
    origin = x - {"start": 0.0, "middle": 0.5, "end": 1.0}[text.attrib.get("text-anchor", "start")] * ink.advance * size
    left, right = origin + ink.left * size, origin + ink.right * size
    top, bottom = y + (shift - ink.top) * size, y + (shift - ink.bottom) * size
    return [(left, top), (right, top), (right, bottom), (left, bottom)]


def _annotation_label_crossings(drawing) -> dict[str, int]:
    """For each annotation label, how many feature leaders and feature labels its painted letters meet."""

    from gbdraw.layout.text_geometry import convex_polygons_intersect, segment_intersects_polygon

    root = ElementTree.fromstring(drawing.tostring())
    groups = {element.attrib.get("id"): element for element in root.iter() if element.tag.endswith("}g")}
    # The plot groups share one frame: the annotation track, the leaders and the feature labels.
    assert len({groups[name].attrib["transform"] for name in ("label_leaders", "label_text")}) == 1
    leaders = [
        ((float(line.attrib["x1"]), float(line.attrib["y1"])), (float(line.attrib["x2"]), float(line.attrib["y2"])))
        for line in groups["label_leaders"]
    ]
    # The letters, not their halos.
    feature_labels = [_feature_label_corners(text) for text in groups["label_text"] if text.tag.endswith("}text")]
    crossings: dict[str, int] = {}
    for element in root.iter():
        label = element.attrib.get("data-gbdraw-annotation-label")
        if not label:
            continue
        (text,) = [child for child in element if child.tag.endswith("}text")]
        parts = _drawn_label(text, root).parts
        crossings[label] = sum(
            any(segment_intersects_polygon(leader, part, touching=False) for part in parts) for leader in leaders
        ) + sum(any(convex_polygons_intersect(part, box) for part in parts) for box in feature_labels)
    return crossings


def _hmmtdna_with_brackets(*regions: tuple[str, int, int]):
    from Bio import SeqIO

    record = SeqIO.read(Path(__file__).parent / "test_inputs" / "HmmtDNA.gbk", "genbank")
    style = RegionAnnotationStyle(label_font_size=14)
    annotations = tuple(
        RegionAnnotation(f"r{index}", CoordinateSpan(None, start, end), label=text, mark="bracket", style=style)
        for index, (text, start, end) in enumerate(regions)
    )
    return build_circular_diagram(
        record,
        options=CircularDiagramOptions(
            annotations=AnnotationOptions(sets=(AnnotationSet("regions", annotations),)),
            config_overrides={"labels.circular.scope": "outer"},
        ),
    )


def test_circular_annotation_labels_move_along_their_region_off_feature_leaders() -> None:
    drawing = _hmmtdna_with_brackets(("rRNA genes", 648, 3229), ("ND1 region", 3307, 4262), ("COX block", 5904, 9207))

    # Feature leaders crossed "ND1 region" and "COX block" at their midpoints (OV-1173); each label now moves
    # along its own region to the clear place nearest its midpoint (the Linear rule of OV-845).
    assert _annotation_label_crossings(drawing) == {"rRNA genes": 0, "ND1 region": 0, "COX block": 0}
    root = ElementTree.fromstring(drawing.tostring())
    for element in root.iter():
        label = element.attrib.get("data-gbdraw-annotation-label")
        if label:
            (text,) = [child for child in element if child.tag.endswith("}text")]
            drawn = _drawn_label(text, root)
            # The label's centre stays over its own region.
            angle = math.degrees(math.atan2(drawn.y, drawn.x)) + 90.0
            position = (angle % 360.0) / 360.0 * 16569
            start, end = {"rRNA genes": (648, 3229), "ND1 region": (3307, 4262), "COX block": (5904, 9207)}[label]
            assert start - 1 <= position <= end + 1, label


def _feature_label_drawings(root) -> list[tuple[str, str, tuple[tuple[float, float], ...]]]:
    """What the label style draws with each feature label's letters, in the frame of the annotation labels: its
    white box (a rect bound to the label) and its underline (a horizontal line bound to it, with its stroke and
    caps)."""

    (text_layer,) = [element for element in root.iter() if element.attrib.get("id") == "label_text"]
    drawn = []
    for element in text_layer:
        feature = element.attrib.get("data-label-feature-id")
        if feature and element.tag.endswith("}rect"):
            assert "transform" not in element.attrib
            x, y, w, h = (float(element.attrib[key]) for key in ("x", "y", "width", "height"))
            drawn.append(("box", feature, ((x, y), (x + w, y), (x + w, y + h), (x, y + h))))
        elif feature and element.tag.endswith("}line") and element.attrib["y1"] == element.attrib["y2"]:
            x1, x2, y = float(element.attrib["x1"]), float(element.attrib["x2"]), float(element.attrib["y1"])
            half = 0.5 * float(element.attrib["stroke-width"])
            cap = half if element.attrib.get("stroke-linecap") in ("round", "square") else 0.0
            left, right = min(x1, x2) - cap, max(x1, x2) + cap
            if x1 != x2:
                drawn.append(("underline", feature, ((left, y - half), (right, y - half), (right, y + half), (left, y + half))))
    return drawn


@pytest.mark.parametrize(("line", "protection"), [("under", "box"), ("none", "box"), ("under", "halo"), ("none", "halo")])
def test_circular_annotation_labels_keep_off_what_the_feature_label_style_draws(
    line: str, protection: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    """OV-1358: what an annotation label keeps off is everything drawn for a feature label, its letters and, by
    its style, its underline and its white box (the Linear rule, OV-1357): each drawn underline and box is one of
    the obstacle parts (a halo adds none)."""

    from Bio import SeqIO

    import gbdraw.diagrams.circular.assemble as circular_assemble

    seen = []
    original = circular_assemble._feature_label_obstacles

    def spy(*args):
        seen.append(original(*args))
        return seen[-1]

    monkeypatch.setattr(circular_assemble, "_feature_label_obstacles", spy)
    record = SeqIO.read(Path(__file__).parent / "test_inputs" / "HmmtDNA.gbk", "genbank")
    region = RegionAnnotation("r0", CoordinateSpan(None, 8900, 9800), label="Control region", mark="bracket")
    drawing = build_circular_diagram(
        record,
        options=CircularDiagramOptions(
            annotations=AnnotationOptions(sets=(AnnotationSet("regions", (region,)),)),
            config_overrides={
                "labels.circular.scope": "outer", "labels.circular.line": line, "labels.circular.protection": protection,
            },
        ),
    )
    drawn = _feature_label_drawings(ElementTree.fromstring(drawing.tostring()))
    expected = ({"box"} if protection == "box" else set()) | ({"underline"} if line == "under" else set())
    assert {kind for kind, _feature, _polygon in drawn} == expected
    (obstacles,) = seen
    parts = [part for shape in obstacles.shapes for part in shape]
    for kind, feature, polygon in drawn:
        assert any(
            len(part) == 4 and all(math.dist(a, b) < 1e-6 for a, b in zip(sorted(part), sorted(polygon))) for part in parts
        ), (kind, feature)


def test_circular_annotation_label_without_a_clear_place_takes_the_fewest_conflicts() -> None:
    from gbdraw.annotations import annotation_track_params_from_mapping, layout_annotation_track, resolve_annotations
    from gbdraw.layout.text_geometry import convex_polygons_intersect
    from gbdraw.layout.text_geometry import grown_convex_polygon
    from gbdraw.render.drawers.circular.annotations import CircularLabelObstacles, circular_annotation_labels

    record = SeqRecord(Seq("ACGT" * 2500), id="r1", name="r1")
    region = RegionAnnotation("r0", CoordinateSpan(None, 2000, 3000), label="Region", mark="bracket")
    bundle = resolve_annotations((AnnotationSet("s", (region,)),), [record], mode="circular")
    params = annotation_track_params_from_mapping({"set_id": "s"})
    track = layout_annotation_track(
        "s", "s", bundle.annotations, record_lengths={0: 10000}, params=params, label_font_size=10.0
    )

    def labels(*polygons):
        return circular_annotation_labels(
            track,
            record_index=0,
            record_length=10000,
            inner_radius_px=300.0,
            outer_radius_px=320.0,
            side="outside",
            font_family="Liberation Sans",
            params=params,
            obstacles=CircularLabelObstacles(shapes=tuple((polygon,) for polygon in polygons), clearance=3.0) if polygons else None,
        )

    (default,) = labels()
    everywhere = ((-1000.0, -1000.0), (1000.0, -1000.0), (1000.0, 1000.0), (-1000.0, 1000.0))
    # Nothing is clear anywhere: the label stays at its midpoint, drawn.
    assert labels(everywhere) == (default,)
    # One more obstacle on the midpoint: the label moves to the nearest place with one conflict, over its region.
    box = default.ink_bounds
    on_midpoint = ((box.min_x, box.min_y), (box.max_x, box.min_y), (box.max_x, box.max_y), (box.min_x, box.max_y))
    (moved,) = labels(everywhere, on_midpoint)
    assert (moved.x, moved.y) != (default.x, default.y)
    assert not any(convex_polygons_intersect(part, on_midpoint) for part in moved.parts)
    position = (math.degrees(math.atan2(moved.y, moved.x)) + 90.0) % 360.0 / 360.0 * 10000
    assert 2000 <= position <= 3000
    # The nearest such place, sampled every pixel of arc: the letters keep the clearance and less than a pixel more.
    assert not any(convex_polygons_intersect(grown_convex_polygon(part, 3.0), on_midpoint) for part in moved.parts)
    assert any(convex_polygons_intersect(grown_convex_polygon(part, 4.0), on_midpoint) for part in moved.parts)


def _draw_span_regions(slots: tuple[str, ...], style: RegionAnnotationStyle, *, mark: str = "bracket"):
    """Two labelled regions over a record whose GC plots reach their band edges (REVIEW-1192 F1): GC 10 % on the
    first half and 90 % on the second."""

    rng = random.Random(7)
    sequence = "".join(
        rng.choice("GC") if rng.random() < (0.1 if (index // 5000) == 0 else 0.9) else rng.choice("AT")
        for index in range(10000)
    )
    record = SeqRecord(Seq(sequence), id="r1", name="r1")
    record.features = [
        SeqFeature(FeatureLocation(start, start + 600, strand=1 if index % 2 else -1), type="CDS")
        for index, start in enumerate(range(100, 9500, 900))
    ]
    regions = (
        RegionAnnotation("a", CoordinateSpan(None, 1000, 3000), label="Region gyp", mark=mark, style=style),
        RegionAnnotation("b", CoordinateSpan(None, 6000, 7500), label="IRb", mark=mark, style=style),
    )
    return build_circular_diagram(
        record,
        options=CircularDiagramOptions(
            annotations=AnnotationOptions(sets=(AnnotationSet("regions", regions),)),
            config_overrides={"labels.circular.scope": "none"},
            tracks=CircularRequestTrackOptions(circular_track_slots=slots),
        ),
    )


def _slot_bands(drawing) -> dict[str, tuple[float, float]]:
    """Each slot's inner and outer radius, from the drawing's track slot geometry."""

    geometry = drawing._gbdraw_track_slot_geometry["records"][0]
    return {
        slot["slotId"]: (
            float(slot["radiusFactor"]) * float(geometry["axisRadiusPx"]) - 0.5 * float(slot["widthPx"]),
            float(slot["radiusFactor"]) * float(geometry["axisRadiusPx"]) + 0.5 * float(slot["widthPx"]),
        )
        for slot in geometry["slots"]
    }


def _annotation_letter_radii(drawing) -> list[tuple[float, float]]:
    """How near the centre and how far from it each annotation label's painted letters come."""

    root = ElementTree.fromstring(drawing.tostring())
    return [
        _ink_radii(next(child for child in element if child.tag.endswith("}text")), root)
        for element in root.iter()
        if element.attrib.get("data-gbdraw-annotation-label")
    ]


@pytest.mark.parametrize(
    "style",
    [RegionAnnotationStyle(label_offset=4.0), RegionAnnotationStyle(label_offset=4.0, label_font_size=40, label_orientation="radial")],
    ids=["14px", "40px-radial"],
)
def test_circular_span_label_room_moves_its_ring_by_the_room_and_keeps_the_features_on_the_axis(
    style: RegionAnnotationStyle,
) -> None:
    # A span on an automatic inside ring, with an automatic ring between it and the features (REVIEW-1192 delta D1).
    def draw(cover_anchor: str):
        return _draw_span_regions(
            (
                "features:features",
                "skew:dinucleotide_skew@side=inside",
                "gc:dinucleotide_content@side=inside",
                f"regions:annotations@side=overlay,set_id=regions,anchor_slot=gc,cover_anchor={cover_anchor},z=5",
            ),
            style,
        )

    plain, span = _slot_bands(draw("false")), draw("true")
    bands = _slot_bands(span)
    # The span's room never moves the features off the axis.
    assert bands["features"] == pytest.approx(plain["features"], abs=1e-6)
    # The letters keep label_offset from the ring beyond them, and gc moves in no further than they need: they keep
    # less than the rings' facing gap (1 % of the axis radius) more.
    facing_gap = 0.01 * float(span._gbdraw_track_slot_geometry["records"][0]["axisRadiusPx"])
    letters = _annotation_letter_radii(span)
    assert len(letters) == 2
    gap = bands["skew"][0] - max(far for _near, far in letters)
    assert 4.0 - 1e-3 <= gap <= 4.0 + facing_gap + 1e-3, gap
    assert min(near for near, _far in letters) - bands["gc"][1] >= 4.0 - 1e-3


def test_circular_span_label_room_keeps_a_pinned_ring_and_moves_the_automatic_rings() -> None:
    # REVIEW-1192 delta D2: a ring pinned outside the axis keeps its radius and width (GX-17); the automatic
    # features ring under it moves in, so the letters keep label_offset from it.
    def draw(cover_anchor: str):
        return _draw_span_regions(
            (
                "outer_gc:dinucleotide_content@side=outside,r=1.08",
                "features:features",
                "gc:dinucleotide_content@side=inside",
                f"regions:annotations@side=overlay,set_id=regions,anchor_slot=features,cover_anchor={cover_anchor},z=5",
            ),
            RegionAnnotationStyle(label_offset=4.0),
            mark="line",
        )

    plain, span = _slot_bands(draw("false")), draw("true")
    bands = _slot_bands(span)
    assert bands["outer_gc"] == pytest.approx(plain["outer_gc"], abs=1e-6)
    letters = _annotation_letter_radii(span)
    assert letters and all(bands["outer_gc"][0] - far >= 4.0 - 1e-3 for _near, far in letters), letters


def test_circular_rings_keep_their_places_when_they_cannot_make_the_label_room() -> None:
    # OV-1207: a ring pinned below an annotation slot, already at its smallest width, cannot make the labels' room;
    # the layout keeps the places it takes without the room (explicit radii win) instead of failing.
    def bands(show_labels: str) -> dict[str, tuple[float, float]]:
        slots = (
            "features:features",
            f"regions:annotations@side=inside,set_id=regions,w=20px,show_labels={show_labels}",
            "gc:dinucleotide_content@side=inside,r=0.8",
        )
        return _slot_bands(_draw_span_regions(slots, RegionAnnotationStyle(label_offset=4.0)))

    assert bands("true") == pytest.approx(bands("false"), abs=1e-6)


def test_circular_span_label_room_on_an_inside_ring_leaves_the_outside_rows_in_place() -> None:
    # REVIEW-1192 delta D3: a span on an inside ring moves that ring by its room; the tick ring outside and the
    # features keep the places they have without the span.
    def draw(cover_anchor: str):
        return _draw_span_regions(
            (
                "ticks:ticks@side=outside",
                "features:features",
                "gc:dinucleotide_content@side=inside",
                f"regions:annotations@side=overlay,set_id=regions,anchor_slot=gc,cover_anchor={cover_anchor},z=5",
            ),
            RegionAnnotationStyle(label_offset=4.0),
        )

    plain, span = _slot_bands(draw("false")), draw("true")
    bands = _slot_bands(span)
    for slot_id in ("ticks", "features"):
        assert bands[slot_id] == pytest.approx(plain[slot_id], abs=1e-6), slot_id
    letters = _annotation_letter_radii(span)
    assert letters and all(bands["features"][0] - far >= 4.0 - 1e-3 for _near, far in letters), letters


@pytest.mark.parametrize(
    "style",
    [RegionAnnotationStyle(label_offset=4.0), RegionAnnotationStyle(label_offset=4.0, label_font_size=40, label_orientation="radial")],
    ids=["14px", "40px-radial"],
)
def test_circular_span_label_room_past_the_axis_moves_the_outside_rows(style: RegionAnnotationStyle) -> None:
    # A span on the features, the inside row next to the axis: its room reaches past the axis, so the features stay
    # on the axis and the automatic outside rings move out of the letters' way (REVIEW-1192 F1, delta D3).
    def draw(cover_anchor: str):
        return _draw_span_regions(
            (
                "o2:dinucleotide_skew@side=outside",
                "o1:dinucleotide_content@side=outside",
                "features:features",
                "gc:dinucleotide_content@side=inside",
                f"regions:annotations@side=overlay,set_id=regions,anchor_slot=features,cover_anchor={cover_anchor},z=5",
            ),
            style,
        )

    plain, span = _slot_bands(draw("false")), draw("true")
    bands = _slot_bands(span)
    assert bands["features"] == pytest.approx(plain["features"], abs=1e-6)
    letters = _annotation_letter_radii(span)
    assert letters and all(bands["o1"][0] - far >= 4.0 - 1e-3 for _near, far in letters), letters


@pytest.mark.parametrize(
    ("slots", "ring"),
    [
        # OV-1360: a span over split feature lanes holds its room past them, so the automatic ring outside moves out.
        (
            (
                "features:features@lane_direction=split",
                "regions:annotations@side=overlay,set_id=regions,anchor_slot=features,cover_anchor=true,z=5",
                "gc:dinucleotide_content@side=outside",
            ),
            "gc",
        ),
        # OV-1360: a span from split feature lanes to an inside ring keys its room to the features, where its labels go.
        (
            (
                "o1:dinucleotide_content@side=outside",
                "features:features@lane_direction=split",
                "gc:dinucleotide_content@side=inside",
                "regions:annotations@side=overlay,set_id=regions,anchor_slot=features,to_slot=gc,z=5",
            ),
            "o1",
        ),
        # OV-1359: an overlay that spans nothing measures its room from its own band and keeps it past its anchor.
        (
            (
                "o2:dinucleotide_skew@side=outside",
                "o1:dinucleotide_content@side=outside",
                "features:features",
                "gc:dinucleotide_content@side=inside",
                "regions:annotations@side=overlay,set_id=regions,anchor_slot=features,cover_anchor=false,z=5",
            ),
            "o1",
        ),
        # OV-1360: a span from split feature lanes to an outside ring keys its room to that ring, which lies beyond the
        # features, so the automatic ring outside it moves out.
        (
            (
                "o2:dinucleotide_skew@side=outside",
                "o1:dinucleotide_content@side=outside",
                "features:features@lane_direction=split",
                "regions:annotations@side=overlay,set_id=regions,anchor_slot=features,to_slot=o1,z=5",
                "gc:dinucleotide_content@side=inside",
            ),
            "o2",
        ),
    ],
    ids=["span-over-split-features", "span-from-split-features", "overlay-spanning-nothing", "span-to-an-outside-ring"],
)
def test_circular_overlay_labels_keep_their_room_from_the_automatic_ring_outside(
    slots: tuple[str, ...], ring: str
) -> None:
    drawing = _draw_span_regions(slots, RegionAnnotationStyle(label_offset=4.0))
    bands = _slot_bands(drawing)
    letters = _annotation_letter_radii(drawing)
    assert len(letters) == 2 and all(bands[ring][0] - far >= 4.0 - 1e-3 for _near, far in letters), (letters, bands)


@pytest.mark.parametrize(
    ("slots", "dev_width", "dev_bands"),
    [
        (
            ("o1:dinucleotide_content@side=outside", "features:features",
             "regions:annotations@side=overlay,set_id=regions,anchor_slot=features,w=60px,show_labels=false,z=5",
             "gc:dinucleotide_content@side=inside"),
            1144.23,
            {"o1": (393.9, 468.0), "features": (349.05, 386.1), "regions": (337.575, 397.575), "gc": (271.05, 345.15)},
        ),
        (
            ("features:features",
             "regions:annotations@side=overlay,set_id=regions,anchor_slot=features,w=80px,show_labels=false,z=5"),
            929.102,
            {"features": (349.05, 386.1), "regions": (327.575, 407.575)},
        ),
        (
            ("o1:dinucleotide_content@side=outside", "features:features", "gc:dinucleotide_content@side=inside",
             "regions:annotations@side=overlay,set_id=regions,anchor_slot=gc,w=120px,show_labels=false,z=5",
             "sk:dinucleotide_skew@side=inside"),
            1153.24,
            {"o1": (393.9, 468.0), "features": (349.05, 386.1), "gc": (271.05, 345.15), "regions": (248.1, 368.1),
             "sk": (193.05, 267.15)},
        ),
    ],
    ids=["under-an-outside-ring", "outermost", "on-an-inside-ring"],
)
def test_circular_label_less_overlay_wider_than_its_anchor_keeps_every_band(
    slots: tuple[str, ...], dev_width: float, dev_bands: dict[str, tuple[float, float]]
) -> None:
    # OV-1359: an overlay that spans nothing and has no labels needs no room, so a tint wider than its anchor leaves
    # every band and the canvas where they were before the labels' room (dev dc1fae14a).
    drawing = _draw_span_regions(slots, RegionAnnotationStyle(label_offset=4.0))
    bands = _slot_bands(drawing)
    assert {slot: (round(inner, 3), round(outer, 3)) for slot, (inner, outer) in bands.items()} == dev_bands
    assert float(str(drawing["width"]).removesuffix("px")) == pytest.approx(dev_width, abs=1e-3)


def test_circular_rings_that_cannot_make_the_room_of_split_lanes_draw_and_warn() -> None:
    # OV-1360: split feature lanes hold the span's room, which reaches outside the axis; the pinned ring leaves the
    # automatic ring between it and the features no space for it, so that room is dropped, the figure is drawn, and
    # the labels on the ring are counted.
    slots = (
        "og:dinucleotide_content@side=outside,r=1.36,w=60px",
        "o1:dinucleotide_content@side=outside",
        "features:features@lane_direction=split",
        "regions:annotations@side=overlay,set_id=regions,anchor_slot=features,cover_anchor=true,z=5",
        "gc:dinucleotide_content@side=inside",
    )
    with pytest.warns(LabelLeadersCrossWarning, match=r"\(2 region labels on other tracks\)\. Drawn anyway\."):
        drawing = _draw_span_regions(slots, RegionAnnotationStyle(label_offset=4.0))
    bands = _slot_bands(drawing)
    axis = float(drawing._gbdraw_track_slot_geometry["records"][0]["axisRadiusPx"])
    assert sum(bands["og"]) / 2.0 == pytest.approx(1.36 * axis, abs=1e-6)
    letters = _annotation_letter_radii(drawing)
    assert len(letters) == 2 and all(far > bands["o1"][0] for _near, far in letters), (letters, bands)


def test_circular_rings_that_cannot_make_one_sides_label_room_keep_the_other_sides_room() -> None:
    # REVIEW-1192 delta D4: the pinned GC ring inside cannot make the inside slot's room (OV-1207); the outside
    # annotation slot still keeps its room, so the automatic ring beyond it keeps label_offset from its letters.
    drawing = _draw_span_regions(
        (
            "ring:dinucleotide_content@side=outside",
            "out:annotations@side=outside,set_id=regions,w=20px",
            "features:features",
            "in:annotations@side=inside,set_id=regions,w=20px",
            "gc:dinucleotide_content@side=inside,r=0.8",
        ),
        RegionAnnotationStyle(label_offset=4.0),
    )
    bands = _slot_bands(drawing)
    axis = float(drawing._gbdraw_track_slot_geometry["records"][0]["axisRadiusPx"])
    outside = [(near, far) for near, far in _annotation_letter_radii(drawing) if near > axis]
    assert outside and all(bands["ring"][0] - far >= 4.0 - 1e-3 for _near, far in outside), outside
    assert sum(bands["gc"]) / 2.0 == pytest.approx(0.8 * axis, abs=1e-6)


def test_circular_region_labels_left_no_room_by_set_radii_warn() -> None:
    # LP-D86 (OV-1208): the GC ring set to r=0.8 leaves the inside slot's labels no room, so they are drawn on it;
    # the label fit warning counts them, the one report the command line and the Web read, and offers no fix.
    from gbdraw.api.diagram import label_fit_report
    from gbdraw.labels.fit import label_fit_notice
    from gbdraw.web_support.label_fit import label_fit_report_payload

    slots = (
        "ring:dinucleotide_content@side=outside",
        "out:annotations@side=outside,set_id=regions,w=20px",
        "features:features",
        "in:annotations@side=inside,set_id=regions,w=20px",
        "gc:dinucleotide_content@side=inside,r=0.8",
    )
    with pytest.warns(LabelLeadersCrossWarning, match=r"\(2 region labels on other tracks\)\. Drawn anyway\."):
        drawing = _draw_span_regions(slots, RegionAnnotationStyle(label_offset=4.0))
    gc_inner, gc_outer = _slot_bands(drawing)["gc"]
    axis = float(drawing._gbdraw_track_slot_geometry["records"][0]["axisRadiusPx"])
    inside = [(near, far) for near, far in _annotation_letter_radii(drawing) if far < axis]
    assert len(inside) == 2 and all(far > gc_inner and near < gc_outer for near, far in inside), inside
    report = label_fit_report(drawing)
    assert report is not None and report.drawn and report.region_track_overlaps == 2 and report.rows == ()
    notice = "Region labels lie on other tracks: the set track radii leave them no room (2)."
    assert label_fit_notice(report) == (notice,)
    payload = label_fit_report_payload(report, "circular")
    assert payload["notice"] == [notice]
    ((record_payload),) = payload["records"]
    assert record_payload["labelKind"] == "circular_region" and record_payload["labelCount"] == 4
    assert record_payload["counts"] == "2 region labels on other tracks"


@pytest.mark.parametrize(
    "slots",
    [
        # A split features slot is an overlay but takes its band: the span's letters lie on its feature lanes.
        ("features:features@lane_direction=split", "gc:dinucleotide_content@side=inside,r=0.91,w=18px",
         "regions:annotations@side=overlay,set_id=regions,anchor_slot=gc,cover_anchor=true,z=5"),
        # A ticks slot takes its tick labels' band too: the letters lie on the scale numbers.
        ("ticks:ticks@side=outside,r=1.14", "out:annotations@side=outside,set_id=regions,w=20px", "features:features"),
        # A span to an overlay plot inside it keeps its room past the skew ring, so the set radius takes it.
        ("features:features", "skew:dinucleotide_skew@side=inside,r=0.86,w=18px",
         "ov:dinucleotide_content@side=overlay,r=0.6,w=10px",
         "regions:annotations@side=overlay,set_id=regions,anchor_slot=skew,to_slot=ov,z=5"),
        # A span over split feature lanes keeps its room past them, so a ring set on the letters takes it (OV-1360).
        ("features:features@lane_direction=split",
         "regions:annotations@side=overlay,set_id=regions,anchor_slot=features,cover_anchor=true,z=5",
         "gc:dinucleotide_content@side=outside,r=1.16"),
        # An overlay that spans nothing keeps its room past its anchor, so a ring set on the letters takes it when the
        # set ring inside holds the features on the axis (OV-1359).
        ("o1:dinucleotide_content@side=outside,r=1.11", "features:features", "gc:dinucleotide_content@side=inside,r=0.85",
         "regions:annotations@side=overlay,set_id=regions,anchor_slot=features,cover_anchor=false,z=5"),
    ],
    ids=["split-features", "tick-labels", "span-to-an-overlay-plot", "span-over-split-features",
         "overlay-spanning-nothing"],
)
def test_circular_region_labels_on_the_band_a_track_takes_warn(slots: tuple[str, ...]) -> None:
    with pytest.warns(LabelLeadersCrossWarning, match=r"\(2 region labels on other tracks\)\. Drawn anyway\."):
        _draw_span_regions(slots, RegionAnnotationStyle(label_offset=4.0))


@pytest.mark.parametrize(
    "slots",
    [
        ("ring:dinucleotide_content@side=outside", "out:annotations@side=outside,set_id=regions,w=20px", "features:features",
         "in:annotations@side=inside,set_id=regions,w=20px", "gc:dinucleotide_content@side=inside"),
        ("features:features@lane_direction=split", "gc:dinucleotide_content@side=inside,w=18px",
         "regions:annotations@side=overlay,set_id=regions,anchor_slot=gc,cover_anchor=true,z=5"),
        ("ticks:ticks@side=outside", "out:annotations@side=outside,set_id=regions,w=20px", "features:features"),
        ("features:features", "skew:dinucleotide_skew@side=inside,w=18px", "ov:dinucleotide_content@side=overlay,r=0.6,w=10px",
         "regions:annotations@side=overlay,set_id=regions,anchor_slot=skew,to_slot=ov,z=5"),
        # An overlay that spans nothing, or covers split feature lanes, keeps a room past its anchor (OV-1359, OV-1360).
        ("o1:dinucleotide_content@side=outside", "features:features", "gc:dinucleotide_content@side=inside",
         "regions:annotations@side=overlay,set_id=regions,anchor_slot=features,cover_anchor=false,z=5"),
        ("features:features@lane_direction=split",
         "regions:annotations@side=overlay,set_id=regions,anchor_slot=features,cover_anchor=true,z=5",
         "gc:dinucleotide_content@side=outside"),
        ("o1:dinucleotide_content@side=outside", "features:features@lane_direction=split",
         "gc:dinucleotide_content@side=inside",
         "regions:annotations@side=overlay,set_id=regions,anchor_slot=features,to_slot=gc,z=5"),
    ],
    ids=["automatic-radii", "split-features-automatic", "tick-labels-automatic", "span-to-an-overlay-plot-automatic",
         "overlay-spanning-nothing", "span-over-split-features", "span-from-split-features"],
)
def test_circular_region_labels_with_their_room_do_not_warn(slots: tuple[str, ...]) -> None:
    with warnings.catch_warnings():
        warnings.simplefilter("error", LabelLeadersCrossWarning)
        _draw_span_regions(slots, RegionAnnotationStyle(label_offset=4.0))


def test_circular_span_label_room_past_the_axis_is_dropped_when_an_outside_row_cannot_make_it() -> None:
    # REVIEW-1192 E1: a span on the features reaches past the axis; when the automatic outside ring under a pinned
    # one cannot make that room, the retry drops it and the figure renders as on dev (explicit radii win).
    drawing = _draw_span_regions(
        (
            "og:dinucleotide_content@side=outside,r=1.3,w=60px",
            "o1:dinucleotide_content@side=outside",
            "features:features",
            "gc:dinucleotide_content@side=inside",
            "regions:annotations@side=overlay,set_id=regions,anchor_slot=features,cover_anchor=true,z=5",
        ),
        RegionAnnotationStyle(label_offset=4.0),
    )
    bands = _slot_bands(drawing)
    axis = float(drawing._gbdraw_track_slot_geometry["records"][0]["axisRadiusPx"])
    assert sum(bands["og"]) / 2.0 == pytest.approx(1.3 * axis, abs=1e-6)
    assert bands["og"][1] - bands["og"][0] == pytest.approx(60.0, abs=1e-6)
    assert len(_annotation_letter_radii(drawing)) == 2


@pytest.mark.parametrize("mark", ["highlight", "band", "bracket", "line"])
@pytest.mark.parametrize(
    "span",
    [
        "anchor_slot=features,cover_anchor=true",
        "anchor_slot=features,to_slot=gc",
        "anchor_slot=features,to_slot=outer_gc",
    ],
    ids=["cover_anchor", "to_slot-inside", "to_slot-outside"],
)
def test_circular_span_labels_keep_label_offset_beyond_every_covered_ring(span: str, mark: str) -> None:
    drawing = _draw_span_regions(
        (
            "outer_gc:dinucleotide_content@side=outside",
            "features:features",
            "gc:dinucleotide_content@side=inside",
            f"regions:annotations@side=overlay,set_id=regions,{span},z=5",
        ),
        RegionAnnotationStyle(),
        mark=mark,
    )

    bands = _slot_bands(drawing)
    span_outer = bands["regions"][1]
    gaps = _circular_label_gaps(drawing)
    letters = _annotation_letter_radii(drawing)
    assert len(letters) == 2
    for near, far in letters:
        # The letters keep label_offset beyond the span's outer edge, off every ring it covers (OD-R1174 F1),
        # and beyond their own mark's paint.
        assert near - span_outer >= 4.0 - 1e-3, near - span_outer
        if "to_slot=outer_gc" not in span:
            # The span reserves that room beyond itself, so the automatic ring beyond it keeps label_offset too.
            assert bands["outer_gc"][0] - far >= 4.0 - 1e-3, bands["outer_gc"][0] - far
    assert gaps and all(gap >= 4.0 - 1e-3 for _side, gap in gaps.values()), gaps
