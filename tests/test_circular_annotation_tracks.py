from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

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
from gbdraw.tracks import parse_circular_track_slot


def test_circular_annotation_renderer_parses_width() -> None:
    slot = parse_circular_track_slot("regions:annotations@set_id=s,w=28px")
    assert slot.renderer == "annotations"
    assert slot.width is not None and slot.width.value == 28


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
