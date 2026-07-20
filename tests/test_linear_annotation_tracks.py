from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
import pytest
from xml.etree import ElementTree

from gbdraw.api import (
    AnnotationOptions,
    AnnotationSet,
    CoordinateSpan,
    DiagramOptions,
    HatchStyle,
    RegionAnnotation,
    RegionAnnotationStyle,
    TrackOptions,
    LinearMultiRecordOptions,
    assemble_linear_diagram_from_records,
    build_linear_diagram,
    parse_record_selector,
)
from gbdraw.tracks import (
    LinearTrackSlot,
    normalize_linear_track_slots,
    parse_linear_track_slot,
)


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


def test_linear_highlight_automatically_renders_behind_features() -> None:
    record = SeqRecord(Seq("A" * 1000), id="r1", name="r1")
    record.features = [
        SeqFeature(FeatureLocation(120, 250, strand=1), type="CDS")
    ]
    annotation = RegionAnnotation(
        "highlighted",
        CoordinateSpan(None, 100, 300),
        label="Highlighted region",
        mark="highlight",
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
    slots = {
        slot["slotId"]: slot
        for slot in drawing._gbdraw_track_slot_geometry["records"][0]["slots"]
    }

    assert slots["annotations_1"]["side"] == "overlay"
    assert 'data-gbdraw-annotation-mark="highlight"' in svg
    assert 'fill="#94a3b8"' in svg
    assert svg.index('data-gbdraw-annotation-id="highlighted"') < svg.index('<g id="r1"')
    root = ElementTree.fromstring(svg)
    highlight_group = next(
        element
        for element in root.iter()
        if element.attrib.get("data-gbdraw-annotation-id") == "highlighted"
    )
    highlight_rect = next(element for element in highlight_group if element.tag.endswith("rect"))
    assert float(highlight_rect.attrib["height"]) >= 14.0


def test_linear_mixed_annotation_set_splits_highlight_from_lane_marks() -> None:
    record = SeqRecord(Seq("A" * 1000), id="r1", name="r1")
    annotations = (
        RegionAnnotation("highlighted", CoordinateSpan(None, 100, 300), mark="highlight"),
        RegionAnnotation("bracketed", CoordinateSpan(None, 400, 600), mark="bracket"),
    )
    drawing = build_linear_diagram(
        [record],
        options=DiagramOptions(
            annotations=AnnotationOptions(
                sets=(AnnotationSet("regions", annotations),)
            )
        ),
    )
    slots = {
        slot["slotId"]: slot
        for slot in drawing._gbdraw_track_slot_geometry["records"][0]["slots"]
    }

    assert slots["annotations_1"]["side"] == "above"
    assert slots["annotations_1_highlight"]["side"] == "overlay"


def test_linear_annotation_overlay_requires_anchor_and_consistent_z() -> None:
    annotation = LinearTrackSlot(
        "a", "annotations", side="overlay", z=1,
        params={"set_id": "s", "anchor_slot": "features", "layer": "foreground"},
    )
    feature = LinearTrackSlot("features", "features", side="overlay", z=0)
    normalized = normalize_linear_track_slots([feature, annotation])

    assert normalized[0].side == "overlay"
    assert normalized[0].reserve is True
    assert normalized[1].annotation is not None


def test_linear_annotation_overlay_rejects_unknown_anchor_in_complete_list() -> None:
    annotation = LinearTrackSlot(
        "a",
        "annotations",
        side="overlay",
        z=1,
        params={"set_id": "s", "anchor_slot": "missing", "layer": "foreground"},
    )

    with pytest.raises(ValueError, match="unknown anchor_slot='missing'"):
        normalize_linear_track_slots([annotation])


def test_linear_annotation_overlay_rejects_annotation_anchor() -> None:
    anchor = LinearTrackSlot(
        "base_notes",
        "annotations",
        side="above",
        z=0,
        params={"set_id": "base"},
    )
    annotation = LinearTrackSlot(
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
        normalize_linear_track_slots([anchor, annotation])


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


@pytest.mark.linear
def test_multi_record_annotation_lanes_use_final_sequence_width() -> None:
    records = [
        SeqRecord(Seq("A" * 1000), id="r1", name="r1"),
        SeqRecord(Seq("A" * 1000), id="r2", name="r2"),
    ]
    annotations = AnnotationOptions(
        sets=(
            AnnotationSet(
                "regions",
                (
                    RegionAnnotation(
                        "left",
                        CoordinateSpan(parse_record_selector("#1"), 100, 140),
                        label="Extremely long annotation label at final record width",
                    ),
                    RegionAnnotation(
                        "right",
                        CoordinateSpan(parse_record_selector("#1"), 300, 340),
                        label="Extremely long annotation label at final record width",
                    ),
                ),
            ),
        )
    )

    drawing = assemble_linear_diagram_from_records(
        records,
        annotation_options=annotations,
        layout=LinearMultiRecordOptions(
            multi_record_positions=("#1@1", "#2@1"),
        ),
        legend="none",
        config_overrides={
            "show_labels": False,
            "show_gc": False,
            "show_skew": False,
        },
    )
    slots = {
        slot["slotId"]: slot
        for slot in drawing._gbdraw_track_slot_geometry["records"][0]["slots"]
    }

    assert slots["annotations_1"]["heightPx"] == pytest.approx(35.0)
