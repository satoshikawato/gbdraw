import re

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
from gbdraw.api.diagram import build_circular_diagram, build_circular_multi_diagram
from gbdraw.tracks import (
    CircularTrackSlot,
    normalize_circular_track_slots,
    parse_circular_track_slot,
)


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


def _build_circular_span(*slots: str, features: tuple[tuple[int, int], ...] = (), **overrides: bool):
    record = SeqRecord(Seq("ACGT" * 2500), id="r1", name="r1")
    record.features = [
        SeqFeature(FeatureLocation(start, end, strand=1), type="CDS") for start, end in features
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


def test_circular_annotation_to_slot_covers_inward_feature_overflow_lanes() -> None:
    overlapping = ((100, 900), (500, 1500), (700, 1700))
    slots = (*_SPAN_STACK, _span_slot())
    flat = _circular_bands(_build_circular_span(*slots, features=overlapping))
    stacked = _circular_bands(
        _build_circular_span(
            *slots,
            features=overlapping,
            **{"canvas.resolve_overlaps": True, "canvas.strandedness": False},
        )
    )

    # The hull uses the primary feature lane, as cover_anchor does. Overflow lanes
    # grow toward the GC slot, which moves inward, so the span still contains them.
    assert stacked["features"] == pytest.approx(flat["features"])
    assert stacked["gc"][1] < flat["gc"][1]
    assert stacked["drawn"] == pytest.approx((stacked["gc"][0], stacked["features"][1]), abs=1e-3)


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
