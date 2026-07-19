from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pytest

from gbdraw.api import (
    AnnotationOptions,
    AnnotationSet,
    CoordinateSpan,
    DiagramOptions,
    RegionAnnotation,
    RegionAnnotationStyle,
    build_circular_diagram,
    build_circular_multi_diagram,
    parse_record_selector,
)
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
        options=DiagramOptions(
            annotations=AnnotationOptions(sets=(AnnotationSet("regions", (annotation,)),))
        ),
    )
    svg = drawing.tostring()
    assert 'data-gbdraw-annotation-id="origin"' in svg
    assert ">Origin<" in svg
    assert svg.count("A ") >= 2


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
    options = DiagramOptions(
        annotations=AnnotationOptions(sets=(annotations,))
    )

    multi_svg = build_circular_multi_diagram(records, options=options).tostring()
    single_svg = build_circular_multi_diagram(records[:1], options=DiagramOptions(
        annotations=AnnotationOptions(
            sets=(AnnotationSet("regions", (annotations.annotations[0],)),)
        )
    )).tostring()

    assert multi_svg.count('data-gbdraw-annotation-id="') == 2
    assert 'data-gbdraw-annotation-id="first"' in single_svg
