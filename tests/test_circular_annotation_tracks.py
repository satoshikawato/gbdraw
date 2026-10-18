import json
import math
import re
from pathlib import Path

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
            corners = _label_ink_corners(text)
            anchor = math.hypot(float(text.attrib["x"]), float(text.attrib["y"]))
            if anchor > 0.5 * (inner + outer):
                gaps[group.attrib["data-gbdraw-annotation-id"]] = ("outside", _nearest_distance(corners) - outer)
            else:
                farthest = max(math.hypot(x, y) for x, y in corners)
                gaps[group.attrib["data-gbdraw-annotation-id"]] = ("inside", inner - farthest)
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
        # The figure's brackets draw their labels inside with label_offset 7.
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
    text = next(element for element in root.iter() if element.tag.rsplit("}", 1)[-1] == "text" and element.text == "rRNA genes")
    (side, mark_gap), = _circular_label_gaps(drawing).values()
    assert side == "outside"
    # The letters keep label_offset beyond the axis paint, so they keep more than that from the tint.
    assert _nearest_distance(_label_ink_corners(text)) - axis_outer == pytest.approx(4.0, abs=1e-6)
    assert mark_gap > 4.0


def _ink_radii(text) -> tuple[float, float]:
    """How near the circle centre and how far from it a drawn label's painted letters come."""

    corners = _label_ink_corners(text)
    return _nearest_distance(corners), max(math.hypot(x, y) for x, y in corners)


@pytest.mark.parametrize("side", ["outside", "inside"])
@pytest.mark.parametrize("orientation", ["tangent", "radial", "horizontal"])
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
    texts = [element for element in root.iter() if element.tag.rsplit("}", 1)[-1] == "text" and element.text in {
        region.label for region in regions
    }]
    assert len(texts) == len(regions)
    for text in texts:
        near, far = _ink_radii(text)
        for slot_id, (inner, outer) in rings.items():
            # The automatic rings keep out of the labels' room: label_offset beyond the letters (OV-1192).
            gap = max(inner - far, near - outer)
            assert gap >= 4.0 - 1e-6, (text.text, slot_id, gap)


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
    feature_labels = [_feature_label_corners(text) for text in groups["label_text"]]
    crossings: dict[str, int] = {}
    for element in root.iter():
        label = element.attrib.get("data-gbdraw-annotation-label")
        if not label:
            continue
        (text,) = [child for child in element if child.tag.endswith("}text")]
        corners = _label_ink_corners(text)
        crossings[label] = sum(segment_intersects_polygon(leader, corners, touching=False) for leader in leaders) + sum(
            convex_polygons_intersect(corners, box) for box in feature_labels
        )
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
            # The label's centre stays over its own region.
            angle = math.degrees(math.atan2(float(text.attrib["y"]), float(text.attrib["x"]))) + 90.0
            position = (angle % 360.0) / 360.0 * 16569
            start, end = {"rRNA genes": (648, 3229), "ND1 region": (3307, 4262), "COX block": (5904, 9207)}[label]
            assert start - 1 <= position <= end + 1, label


def test_circular_annotation_label_without_a_clear_place_takes_the_fewest_conflicts() -> None:
    from gbdraw.annotations import annotation_track_params_from_mapping, layout_annotation_track, resolve_annotations
    from gbdraw.layout.text_geometry import convex_polygons_intersect
    from gbdraw.render.drawers.circular.annotations import CircularLabelObstacles, _grown, circular_annotation_labels

    record = SeqRecord(Seq("ACGT" * 2500), id="r1", name="r1")
    region = RegionAnnotation("r0", CoordinateSpan(None, 2000, 3000), label="Region", mark="bracket")
    bundle = resolve_annotations((AnnotationSet("s", (region,)),), [record], mode="circular")
    params = annotation_track_params_from_mapping({"set_id": "s"})
    track = layout_annotation_track("s", "s", bundle.annotations, record_lengths={0: 10000}, params=params)

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
            obstacles=CircularLabelObstacles(polygons=polygons, clearance=3.0) if polygons else None,
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
    assert not convex_polygons_intersect(moved.corners, on_midpoint)
    position = (math.degrees(math.atan2(moved.y, moved.x)) + 90.0) % 360.0 / 360.0 * 10000
    assert 2000 <= position <= 3000
    # The nearest such place, sampled every pixel of arc: the letters keep the clearance and less than a pixel more.
    assert not convex_polygons_intersect(_grown(moved.corners, 3.0), on_midpoint)
    assert convex_polygons_intersect(_grown(moved.corners, 4.0), on_midpoint)


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
    record = SeqRecord(Seq("ACGT" * 2500), id="r1", name="r1")
    record.features = [
        SeqFeature(FeatureLocation(start, start + 600, strand=1 if index % 2 else -1), type="CDS")
        for index, start in enumerate(range(100, 9500, 900))
    ]
    regions = (
        RegionAnnotation("a", CoordinateSpan(None, 1000, 3000), label="Region gyp", mark=mark),
        RegionAnnotation("b", CoordinateSpan(None, 6000, 7500), label="IRb", mark=mark),
    )
    drawing = build_circular_diagram(
        record,
        options=CircularDiagramOptions(
            annotations=AnnotationOptions(sets=(AnnotationSet("regions", regions),)),
            config_overrides={"labels.circular.scope": "none"},
            tracks=CircularRequestTrackOptions(
                circular_track_slots=(
                    "outer_gc:dinucleotide_content@side=outside",
                    "features:features",
                    "gc:dinucleotide_content@side=inside",
                    f"regions:annotations@side=overlay,set_id=regions,{span},z=5",
                )
            ),
        ),
    )

    geometry = drawing._gbdraw_track_slot_geometry["records"][0]
    (slot,) = [slot for slot in geometry["slots"] if slot["slotId"] == "regions"]
    span_outer = float(slot["radiusFactor"]) * float(geometry["axisRadiusPx"]) + 0.5 * float(slot["widthPx"])
    gaps = _circular_label_gaps(drawing)
    root = ElementTree.fromstring(drawing.tostring())
    for element in root.iter():
        if element.attrib.get("data-gbdraw-annotation-label"):
            (text,) = [child for child in element if child.tag.endswith("}text")]
            near, _far = _ink_radii(text)
            # The letters keep label_offset beyond the span's outer edge, off every ring it covers (OD-R1174 F1),
            # and beyond their own mark's paint.
            assert near - span_outer >= 4.0 - 1e-3, (element.attrib["data-gbdraw-annotation-label"], near - span_outer)
    assert gaps and all(gap >= 4.0 - 1e-3 for _side, gap in gaps.values()), gaps
