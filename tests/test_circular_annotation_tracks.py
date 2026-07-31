from Bio.Seq import Seq
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


def test_circular_origin_annotation_renders_two_safe_paths() -> None:
    record = SeqRecord(Seq("A" * 1000), id="r1", name="r1")
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
