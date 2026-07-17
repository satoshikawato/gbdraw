from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gbdraw.api import (
    AnnotationOptions,
    AnnotationSet,
    CoordinateSpan,
    DiagramOptions,
    HatchStyle,
    RegionAnnotation,
    RegionAnnotationStyle,
    TrackOptions,
    build_linear_diagram,
)
from gbdraw.tracks import LinearTrackSlot, parse_linear_track_slot


def test_linear_annotation_renderer_parses_typed_params() -> None:
    slot = parse_linear_track_slot("regions:annotations@set_id=s,h=24px,overflow=clip")
    assert slot.renderer == "annotations"
    assert slot.height is not None and slot.height.value == 24


def test_linear_annotation_track_renders_metadata_and_hatch() -> None:
    record = SeqRecord(Seq("A" * 1000), id="r1", name="r1")
    annotation = RegionAnnotation(
        "region",
        CoordinateSpan(None, 100, 300),
        label="Region",
        mark="band",
        style=RegionAnnotationStyle(fill="#88aaff", hatch=HatchStyle(cross=True)),
    )
    drawing = build_linear_diagram(
        [record],
        options=DiagramOptions(
            annotations=AnnotationOptions(sets=(AnnotationSet("regions", (annotation,)),))
        ),
    )
    svg = drawing.tostring()
    assert 'data-gbdraw-annotation-id="region"' in svg
    assert "gbdraw-hatch-" in svg
    assert ">Region<" in svg


def test_linear_annotation_overlay_requires_anchor_and_consistent_z() -> None:
    annotation = LinearTrackSlot(
        "a", "annotations", side="overlay", z=1,
        params={"set_id": "s", "anchor_slot": "features", "layer": "foreground"},
    )
    feature = LinearTrackSlot("features", "features", side="overlay", z=0)
    from gbdraw.tracks import normalize_linear_track_slots

    normalized = normalize_linear_track_slots([feature, annotation])
    assert normalized[1].annotation is not None


def test_linear_annotation_clip_policy_creates_clip_path() -> None:
    record = SeqRecord(Seq("A" * 1000), id="r1", name="r1")
    annotations = (
        RegionAnnotation("one", CoordinateSpan(None, 100, 300), lane=0),
        RegionAnnotation("two", CoordinateSpan(None, 400, 600), lane=1),
    )
    drawing = build_linear_diagram(
        [record],
        options=DiagramOptions(
            annotations=AnnotationOptions(
                sets=(AnnotationSet("regions", annotations),)
            ),
            tracks=TrackOptions(
                linear_track_slots=(
                    "regions:annotations@set_id=regions,h=8px,overflow=clip",
                ),
                linear_track_axis_index=1,
            ),
        ),
    )
    svg = drawing.tostring()
    assert 'data-gbdraw-annotation-clipped="true"' in svg
    assert 'clip-path="url(#gbdraw-annotation-clip-' in svg


def test_annotation_hatch_is_used_in_shared_legend_preview() -> None:
    record = SeqRecord(Seq("A" * 1000), id="r1", name="r1")
    annotation = RegionAnnotation(
        "region",
        CoordinateSpan(None, 100, 300),
        mark="band",
        legend_label="Reviewed region",
        style=RegionAnnotationStyle(hatch=HatchStyle(cross=True)),
    )
    drawing = build_linear_diagram(
        [record],
        options=DiagramOptions(
            annotations=AnnotationOptions(
                sets=(AnnotationSet("regions", (annotation,)),)
            )
        ),
    )
    svg = drawing.tostring()
    assert 'data-legend-key="Reviewed region"' in svg
    assert 'fill="url(#gbdraw-hatch-' in svg
