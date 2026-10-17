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
    CoordinateSpan,
    HatchStyle,
    LinearDiagramOptions,
    LinearTrackOptions,
    RegionAnnotation,
    RegionAnnotationStyle,
    LinearMultiRecordOptions,
    parse_record_selector,
)
from gbdraw.api.config import apply_config_overrides
from gbdraw.api.diagram import assemble_linear_diagram_from_records, build_linear_diagram
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
        options=LinearDiagramOptions(
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
        options=LinearDiagramOptions(
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
    root = ElementTree.fromstring(svg)
    elements = list(root.iter())
    highlight_group = next(
        element
        for element in elements
        if element.attrib.get("data-gbdraw-annotation-id") == "highlighted"
    )
    record_group = next(
        element
        for element in elements
        if element.attrib.get("data-gbdraw-record-id") == "r1"
        and element.attrib.get("id", "").startswith("record_group_")
    )
    assert elements.index(highlight_group) < elements.index(record_group)
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
        options=LinearDiagramOptions(
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


def test_linear_annotation_overlay_rejects_spacer_anchor() -> None:
    spacer = parse_linear_track_slot("space:spacer@side=above,h=12px")
    annotation = LinearTrackSlot(
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
        normalize_linear_track_slots([spacer, annotation])


def test_linear_annotation_cover_anchor_expands_overlay_to_feature_band() -> None:
    record = SeqRecord(Seq("A" * 1000), id="r1", name="r1")
    record.features = [
        SeqFeature(FeatureLocation(120, 250, strand=1), type="CDS")
    ]
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

    def band_height(cover_anchor: bool) -> float:
        drawing = assemble_linear_diagram_from_records(
            [record],
            cfg=apply_config_overrides(None, None),
            legend="none",
            selected_features_set=["CDS"],
            linear_track_slots=[
                "features:features@side=overlay",
                (
                    "regions:annotations@side=overlay,h=14px,set_id=regions,"
                    "anchor_slot=features,layer=foreground,padding_px=0,z=1,"
                    f"cover_anchor={'true' if cover_anchor else 'false'}"
                ),
            ],
            linear_track_axis_index=0,
            annotation_options=annotation_options,
        )
        root = ElementTree.fromstring(drawing.tostring())
        group = next(
            element
            for element in root.iter()
            if element.attrib.get("data-gbdraw-annotation-id") == "band"
        )
        rect = next(element for element in group if element.tag.endswith("rect"))
        return float(rect.attrib["height"])

    uncovered_height = band_height(False)
    covered_height = band_height(True)
    assert uncovered_height == pytest.approx(11.76)
    assert covered_height > uncovered_height


def _cover_anchor_drawing(
    track_layout: str, *, highlight: bool = True, slots: tuple[str, ...] = (), label: str | None = None
):
    record = SeqRecord(Seq("ACGT" * 600), id="rec1", name="rec1")
    record.features = [
        SeqFeature(FeatureLocation(120, 250, strand=1), type="CDS"),
        SeqFeature(FeatureLocation(400, 700, strand=-1), type="CDS"),
    ]
    tint = RegionAnnotation("tint", CoordinateSpan(None, 100, 300), label=label, mark="highlight")
    return assemble_linear_diagram_from_records(
        [record],
        cfg=apply_config_overrides(
            None,
            {"labels.linear.scope": "none", "canvas.linear.track_layout": track_layout},
        ),
        legend="none",
        selected_features_set=["CDS"],
        annotation_options=(
            AnnotationOptions(sets=(AnnotationSet("regions", (tint,)),)) if highlight else None
        ),
        **(
            {
                "linear_track_slots": list(slots),
                "linear_track_axis_index": sum("@side=above" in slot for slot in slots),
            }
            if slots
            else {}
        ),
    )


def _cover_anchor_geometry(drawing) -> dict[str, tuple[float, float] | float]:
    """Feature lanes, features origin, record body, and the drawn tint from the axis."""

    record = drawing._gbdraw_track_slot_geometry["records"][0]
    slots = {slot["slotId"]: slot for slot in record["slots"]}
    features = next(slot for slot in slots.values() if slot["renderer"] == "features")
    geometry: dict[str, tuple[float, float] | float] = {
        "lanes": (features["paintBand"]["topPx"], features["paintBand"]["bottomPx"]),
        "origin": float(features["resolvedOriginPx"]),
        "body": (record["recordBodyBand"]["topPx"], record["recordBodyBand"]["bottomPx"]),
    }
    root = ElementTree.fromstring(drawing.tostring())
    marks = [e for e in root.iter() if e.attrib.get("data-gbdraw-annotation-id") == "tint"]
    if not marks:
        return geometry
    overlay = next(slot for slot in slots.values() if slot["renderer"] == "annotations")
    geometry["plan"] = (overlay["paintBand"]["topPx"], overlay["paintBand"]["bottomPx"])
    parents = {child: parent for parent in root.iter() for child in parent}
    rect = next(element for element in marks[0].iter() if element.tag.endswith("rect"))
    top, element = float(rect.attrib["y"]), rect
    while element is not None:
        top += sum(
            float(value)
            for value in re.findall(r"translate\([^,]+,\s*([^)]+)\)", element.attrib.get("transform", ""))
        )
        element = parents.get(element)
    top -= float(record["axisYpx"])
    geometry["drawn"] = (top, top + float(rect.attrib["height"]))
    return geometry


_COVER_ANCHOR_SLOTS = (
    "features:features@side=below",
    "gc:dinucleotide_content@side=below",
    "regions:annotations@side=overlay,set_id=regions,padding_px=0,show_labels=false,"
    "anchor_slot=features,cover_anchor=true,layer=underlay,z=-1",
)


@pytest.mark.parametrize(
    ("track_layout", "slots"),
    [
        pytest.param("below", (), id="highlight-below"),
        pytest.param("above", (), id="highlight-above"),
        pytest.param("middle", (), id="highlight-middle"),
        pytest.param("middle", _COVER_ANCHOR_SLOTS, id="slot-list-features-below"),
        pytest.param(
            "middle",
            (_COVER_ANCHOR_SLOTS[0].replace("below", "above"), *_COVER_ANCHOR_SLOTS[1:]),
            id="slot-list-features-above",
        ),
        pytest.param(
            "middle",
            (_COVER_ANCHOR_SLOTS[0].replace("below", "overlay"), *_COVER_ANCHOR_SLOTS[1:]),
            id="slot-list-features-overlay",
        ),
    ],
)
def test_linear_cover_anchor_covers_the_feature_lanes_and_moves_nothing(
    track_layout: str, slots: tuple[str, ...]
) -> None:
    plain_slots = tuple(slot for slot in slots if not slot.startswith("regions:"))
    plain = _cover_anchor_geometry(
        _cover_anchor_drawing(track_layout, highlight=False, slots=plain_slots)
    )
    covered = _cover_anchor_geometry(_cover_anchor_drawing(track_layout, slots=slots))

    shift = covered["origin"] - plain["origin"]
    lanes = (plain["lanes"][0] + shift, plain["lanes"][1] + shift)
    # The tint covers the lanes as drawn (OV-474), and the plan records that band (OV-476).
    assert covered["drawn"] == pytest.approx(lanes, abs=1e-6)
    assert covered["plan"] == pytest.approx(lanes, abs=1e-6)
    # Its own nominal band reserves nothing: the lanes and the record stay where they were.
    assert shift == pytest.approx(0.0, abs=1e-6)
    assert covered["lanes"] == pytest.approx(plain["lanes"], abs=1e-6)
    assert covered["body"] == pytest.approx(plain["body"], abs=1e-6)


@pytest.mark.parametrize(
    ("show_labels", "label", "room"),
    [
        pytest.param(False, "Region", 0.0, id="labels-off"),
        # label_offset (4) plus half the label font size (10) above an underlay, as for a span.
        pytest.param(True, "Region", 9.0, id="labelled"),
        # No label text, no room (OV-842).
        pytest.param(True, None, 0.0, id="no-label-text"),
    ],
)
def test_linear_cover_anchor_labels_keep_room_beyond_the_lanes(
    show_labels: bool, label: str | None, room: float
) -> None:
    slots = (
        "features:features@side=above",
        "gc:dinucleotide_content@side=below",
        _COVER_ANCHOR_SLOTS[2].replace("show_labels=false", f"show_labels={str(show_labels).lower()}"),
    )
    geometry = _cover_anchor_geometry(_cover_anchor_drawing("middle", slots=slots, label=label))

    assert geometry["drawn"] == pytest.approx(geometry["lanes"], abs=1e-6)
    assert geometry["body"][0] == pytest.approx(geometry["lanes"][0] - room)


@pytest.mark.parametrize("track_layout", ["below", "above"])
@pytest.mark.parametrize("tint", ["repeat_region", "highlight"])
def test_linear_cover_anchor_on_a_record_without_lane_features_stays_off_the_axis(
    track_layout: str, tint: str
) -> None:
    records = []
    for index, features in enumerate(([("CDS", 100, 900)], [])):
        record = SeqRecord(Seq("ACGT" * 600), id=f"r{index + 1}", name=f"r{index + 1}")
        record.annotations["molecule_type"] = "DNA"
        if tint == "repeat_region":
            features = [*features, ("repeat_region", 1000, 1400)]
        record.features = [SeqFeature(FeatureLocation(a, b, strand=1), type=kind) for kind, a, b in features]
        records.append(record)
    regions = tuple(
        RegionAnnotation(f"t{index}", CoordinateSpan(parse_record_selector(f"#{index}"), 50, 1200),
                         label="Reg", mark="highlight")
        for index in (1, 2)
    )
    drawing = assemble_linear_diagram_from_records(
        records,
        cfg=apply_config_overrides(
            None, {"labels.linear.scope": "none", "canvas.linear.track_layout": track_layout}
        ),
        legend="none",
        selected_features_set=["CDS", "repeat_region"],
        annotation_options=(
            AnnotationOptions(sets=(AnnotationSet("regions", regions),)) if tint == "highlight" else None
        ),
    )
    bands = [
        {slot["renderer"]: slot["paintBand"] for slot in record["slots"] if slot["paintBand"] is not None}
        for record in drawing._gbdraw_track_slot_geometry["records"]
    ]
    tints = [(band["annotations"]["topPx"], band["annotations"]["bottomPx"]) for band in bands]
    if tint == "repeat_region":
        # An underlay feature keeps the nominal lanes, so its tint covers them as on the first record.
        assert tints[1] == pytest.approx(tints[0])
        # ... and the tint is drawn there, not only planned there.
        root = ElementTree.fromstring(drawing.tostring())
        parents = {child: parent for parent in root.iter() for child in parent}
        for rect in root.iter():
            if rect.attrib.get("data-gbdraw-auto-feature-underlay") != "true":
                continue
            top, element = float(rect.attrib["y"]), rect
            while element is not None:
                top += sum(
                    float(value)
                    for value in re.findall(r"translate\([^,]+,\s*([^)]+)\)", element.attrib.get("transform", ""))
                )
                element = parents.get(element)
            top -= float(drawing._gbdraw_track_slot_geometry["records"][int(rect.attrib["data-gbdraw-record-index"])]["axisYpx"])
            drawn = (top, top + float(rect.attrib["height"]))
            assert drawn == pytest.approx((9.5, 54.5) if track_layout == "below" else (-54.5, -9.5))
    # Without lanes the tint keeps its own band on its side of the axis, clear of GC.
    if track_layout == "below":
        assert 0.0 < tints[1][0] and tints[1][1] <= bands[1]["dinucleotide_content"]["topPx"]
    else:
        assert tints[1][1] < 0.0


def test_linear_to_slot_span_covers_the_nominal_lanes_of_an_underlay_only_record() -> None:
    records = []
    for index, features in enumerate((("CDS", "repeat_region"), ("repeat_region",), ())):
        record = SeqRecord(Seq("ACGT" * 600), id=f"r{index + 1}", name=f"r{index + 1}")
        record.annotations["molecule_type"] = "DNA"
        record.features = [
            SeqFeature(FeatureLocation(100 + 900 * offset, 900 + 500 * offset, strand=1), type=kind)
            for offset, kind in enumerate(features)
        ]
        records.append(record)
    regions = tuple(
        RegionAnnotation(f"span{index}", CoordinateSpan(parse_record_selector(f"#{index}"), 50, 1200), mark="highlight")
        for index in (1, 2, 3)
    )
    drawing = assemble_linear_diagram_from_records(
        records,
        cfg=apply_config_overrides(None, {"labels.linear.scope": "none"}),
        legend="none",
        selected_features_set=["CDS", "repeat_region"],
        annotation_options=AnnotationOptions(sets=(AnnotationSet("regions", regions),)),
        linear_track_slots=["features:features@side=below", "gc:dinucleotide_content@side=below", _span_slot()],
        linear_track_axis_index=0,
    )
    spans = [
        next((slot["paintBand"]["topPx"], slot["paintBand"]["bottomPx"]) for slot in record["slots"] if slot["slotId"] == "regions")
        for record in drawing._gbdraw_track_slot_geometry["records"]
    ]
    gc_top = next(
        slot["paintBand"]["topPx"]
        for slot in drawing._gbdraw_track_slot_geometry["records"][2]["slots"]
        if slot["slotId"] == "gc"
    )
    assert spans[1] == pytest.approx(spans[0])
    # A record without any lanes: the span keeps to the GC slot.
    assert spans[2][0] > gc_top


def _span_slot(
    spec: str = "anchor_slot=features,to_slot=gc,layer=underlay,z=-1", *, labels: bool = False
) -> str:
    return (
        "regions:annotations@side=overlay,set_id=regions,padding_px=0,"
        f"show_labels={str(labels).lower()},{spec}"
    )


_SPAN_STACK = (
    "features:features@side=overlay",
    "gc:dinucleotide_content@side=below",
    "skew:dinucleotide_skew@side=below",
)


def _build_linear_span(
    *slots: str | LinearTrackSlot,
    record_count: int = 1,
    features: tuple[tuple[int, int], ...] = ((120, 250),),
    overrides: dict[str, object] | None = None,
    **kwargs: object,
):
    records = []
    for index in range(1, record_count + 1):
        record = SeqRecord(Seq("ACGT" * 600), id=f"rec{index}", name=f"rec{index}")
        record.features = [
            SeqFeature(FeatureLocation(start, end, strand=1), type="CDS")
            for start, end in features
        ]
        records.append(record)
    annotations = AnnotationOptions(
        sets=(
            AnnotationSet(
                "regions",
                tuple(
                    RegionAnnotation(
                        f"span{index}",
                        CoordinateSpan(parse_record_selector(f"#{index}"), 100, 300),
                        label=f"Span {index}",
                        mark="highlight",
                    )
                    for index in range(1, record_count + 1)
                ),
            ),
        )
    )
    return assemble_linear_diagram_from_records(
        records,
        cfg=apply_config_overrides(None, {"labels.linear.scope": "none", **(overrides or {})}),
        legend="none",
        selected_features_set=["CDS"],
        linear_track_slots=list(slots),
        linear_track_axis_index=0,
        annotation_options=annotations,
        **kwargs,
    )


def _linear_bands(drawing, record_index: int = 0) -> dict[str, tuple[float, float] | None]:
    """Each slot's band from its record axis, and the drawn ``span`` rectangle.

    The band of a features slot is its lanes; of a GC content or skew slot, the
    plot height around its baseline; of a Depth slot, the plot below its origin.
    """

    record = drawing._gbdraw_track_slot_geometry["records"][record_index]
    slots = {slot["slotId"]: slot for slot in record["slots"]}
    bands: dict[str, tuple[float, float] | None] = {}
    for slot_id, slot in slots.items():
        origin, height = float(slot["resolvedOriginPx"]), float(slot["heightPx"])
        if not slot["dataAvailable"]:
            bands[slot_id] = None
        elif slot["renderer"] == "features":
            bands[slot_id] = (slot["paintBand"]["topPx"], slot["paintBand"]["bottomPx"])
        elif slot["renderer"] in {"depth", "annotations"}:
            bands[slot_id] = (origin, origin + height)
        else:
            bands[slot_id] = (origin - 0.5 * height, origin + 0.5 * height)

    def group_y(element) -> float:
        return float(re.findall(r"translate\([^,]+,([^)]+)\)", element.attrib["transform"])[-1])

    # Every record draws its regions and gc groups, in record order.
    root = ElementTree.fromstring(drawing.tostring())
    span_group, gc_group = (
        [element for element in root.iter() if element.attrib.get("data-gbdraw-slot-id") == slot_id][
            record_index
        ]
        for slot_id in ("regions", "gc")
    )
    rect = next(element for element in span_group.iter() if element.tag.endswith("rect"))
    axis_y = group_y(gc_group) - float(slots["gc"]["resolvedOriginPx"])
    top = group_y(span_group) + float(rect.attrib["y"]) - axis_y
    bands["drawn"] = (top, top + float(rect.attrib["height"]))
    return bands


def test_linear_annotation_to_slot_spans_from_one_slot_band_to_the_other() -> None:
    bands = _linear_bands(_build_linear_span(*_SPAN_STACK, _span_slot()))

    hull = (bands["features"][0], bands["gc"][1])
    assert bands["features"][1] < bands["gc"][0]
    assert bands["regions"] == pytest.approx(hull)
    assert bands["drawn"] == pytest.approx(hull)

    across = _linear_bands(
        _build_linear_span(
            *_SPAN_STACK, _span_slot("anchor_slot=gc,to_slot=skew,layer=underlay,z=-1")
        )
    )
    assert across["drawn"] == pytest.approx((across["gc"][0], across["skew"][1]))


@pytest.mark.parametrize(
    "spec",
    [
        "anchor_slot=gc,to_slot=features,layer=underlay,z=-1",
        "anchor_slot=features,to_slot=gc,layer=underlay,z=-1,cover_anchor=true",
        "anchor_slot=features,to_slot=gc,layer=underlay,z=-1,h=3px",
        "anchor_slot=features,to_slot=gc,layer=underlay,z=-1,h=200px",
        "anchor_slot=features,to_slot=gc,layer=foreground,z=1",
    ],
)
def test_linear_annotation_to_slot_ignores_order_cover_anchor_and_own_height(spec: str) -> None:
    reference = _linear_bands(_build_linear_span(*_SPAN_STACK, _span_slot()))
    bands = _linear_bands(_build_linear_span(*_SPAN_STACK, _span_slot(spec)))

    assert bands["drawn"] == pytest.approx(reference["drawn"])
    assert {key: bands[key] for key in ("features", "gc", "skew")} == {
        key: reference[key] for key in ("features", "gc", "skew")
    }


def test_linear_annotation_to_slot_padding_shrinks_both_edges() -> None:
    bands = _linear_bands(
        _build_linear_span(
            *_SPAN_STACK,
            "regions:annotations@side=overlay,set_id=regions,padding_px=4,show_labels=false,"
            "anchor_slot=features,to_slot=gc,layer=underlay,z=-1",
        )
    )

    assert bands["drawn"] == pytest.approx((bands["features"][0] + 4.0, bands["gc"][1] - 4.0))


@pytest.mark.parametrize(
    ("spec", "message"),
    [
        ("anchor_slot=features,to_slot=missing,layer=underlay,z=-1", "references unknown to_slot='missing'"),
        ("anchor_slot=features,to_slot=off,layer=underlay,z=-1", "references unknown to_slot='off'"),
        ("anchor_slot=features,to_slot=space,layer=underlay,z=-1", "to_slot 'space' has no drawable band"),
        ("anchor_slot=features,to_slot=notes,layer=underlay,z=-1", "cannot use annotation slot 'notes' as to_slot"),
        ("anchor_slot=features,to_slot=gc,layer=underlay,z=0", "must have z less than anchor 'features'"),
        ("anchor_slot=features,to_slot=low,layer=underlay,z=-1", "must have z less than to_slot 'low'"),
        ("anchor_slot=low,to_slot=features,layer=underlay,z=-1", "must have z less than anchor 'low'"),
        ("anchor_slot=features,to_slot=high,layer=foreground,z=1", "must have z greater than to_slot 'high'"),
        ("anchor_slot=features,to_slot=,layer=underlay,z=-1", "to_slot cannot be empty"),
        ("anchor_slot=features,to_slot=features,layer=underlay,z=-1", "to_slot must name a slot other than its anchor_slot 'features'"),
    ],
)
def test_linear_annotation_to_slot_reference_errors(spec: str, message: str) -> None:
    slots = [
        parse_linear_track_slot(text)
        for text in (
            "features:features@side=overlay",
            "gc:dinucleotide_content",
            "off:dinucleotide_skew@enabled=false",
            "space:spacer@h=8px",
            "notes:annotations@set_id=regions",
            "low:depth@z=-2",
            "high:dinucleotide_skew@nt=AT,z=2",
        )
    ]

    with pytest.raises(ValueError, match=f"annotation .*slot 'regions'.*{message}|^{message}"):
        normalize_linear_track_slots([*slots, parse_linear_track_slot(_span_slot(spec))])


def test_linear_annotation_to_slot_requires_overlay_side() -> None:
    slot = LinearTrackSlot(
        "regions", "annotations", side="above", params={"set_id": "regions", "to_slot": "gc"}
    )

    with pytest.raises(ValueError, match="annotation slot 'regions' uses to_slot without side=overlay"):
        normalize_linear_track_slots([parse_linear_track_slot("gc:dinucleotide_content"), slot])


@pytest.mark.parametrize("spanning", [False, True])
def test_linear_overlay_annotation_reads_slot_references_in_any_key_case(spanning: bool) -> None:
    def bands(anchor_key: str, to_key: str) -> dict[str, tuple[float, float] | None]:
        params: dict[str, object] = {
            "set_id": "regions",
            "layer": "underlay",
            "show_labels": False,
            anchor_key: "gc",
        }
        if spanning:
            params[to_key] = "skew"
        overlay = LinearTrackSlot("regions", "annotations", side="overlay", z=-1, params=params)
        return _linear_bands(
            _build_linear_span(*(parse_linear_track_slot(text) for text in _SPAN_STACK), overlay)
        )

    assert bands("Anchor_Slot", "To_Slot") == bands("anchor_slot", "to_slot")


def test_linear_annotation_to_slot_covers_every_feature_lane() -> None:
    overlapping = ((100, 900), (500, 1500), (700, 1700))
    slots = (*_SPAN_STACK, _span_slot())
    flat = _linear_bands(_build_linear_span(*slots, features=overlapping))
    stacked = _linear_bands(
        _build_linear_span(
            *slots,
            features=overlapping,
            overrides={"canvas.resolve_overlaps": True, "canvas.strandedness": False},
        )
    )

    assert stacked["features"][1] - stacked["features"][0] > flat["features"][1] - flat["features"][0]
    assert stacked["drawn"] == pytest.approx((stacked["features"][0], stacked["gc"][1]))


def test_linear_annotation_to_slot_follows_a_features_slot_off_the_axis() -> None:
    stack = ("features:features@side=below", "gc:dinucleotide_content@side=below")
    # A second overlay on the features slot reserves its own band and pushes that slot away from the axis.
    pushing = (
        "pin:annotations@side=overlay,set_id=regions,show_labels=false,h=40px,"
        "anchor_slot=features,layer=underlay,z=-2"
    )

    def features_origin(drawing) -> float:
        slots = drawing._gbdraw_track_slot_geometry["records"][0]["slots"]
        return float(next(slot for slot in slots if slot["slotId"] == "features")["resolvedOriginPx"])

    plain = _build_linear_span(*stack, _span_slot())
    pushed = _build_linear_span(*stack, pushing, _span_slot())
    shift = features_origin(pushed) - features_origin(plain)
    bands = _linear_bands(pushed)

    assert shift > 0
    assert bands["drawn"] == pytest.approx((_linear_bands(plain)["features"][0] + shift, bands["gc"][1]))


@pytest.mark.parametrize(
    ("spec", "grows"),
    [
        pytest.param("anchor_slot=features,to_slot=gc,layer=underlay,z=-1", True, id="ends-at-record-top"),
        pytest.param("anchor_slot=gc,to_slot=skew,layer=underlay,z=-1", False, id="ends-inside-record"),
    ],
)
def test_linear_annotation_to_slot_labels_add_room_only_at_the_record_edge(
    spec: str, grows: bool
) -> None:
    def build(labels: bool):
        drawing = _build_linear_span(*_SPAN_STACK, _span_slot(spec, labels=labels), record_count=2)
        bodies = [
            (record["recordBodyBand"]["topPx"], record["recordBodyBand"]["bottomPx"])
            for record in drawing._gbdraw_track_slot_geometry["records"]
        ]
        return [_linear_bands(drawing, index) for index in range(2)], bodies

    plain_bands, plain_bodies = build(False)
    labelled_bands, labelled_bodies = build(True)

    # Labels do not move a slot or the span.
    assert labelled_bands == plain_bands
    # label_offset (4) plus half the label font size (10) above an underlay span.
    overhang = 9.0 if grows else 0.0
    assert labelled_bodies == pytest.approx([(top - overhang, bottom) for top, bottom in plain_bodies])


def test_linear_annotation_to_slot_uses_the_bands_each_record_has() -> None:
    depth = pd.DataFrame(
        {
            "reference_name": ["rec1"] * 2400,
            "position": list(range(1, 2401)),
            "depth": [20.0] * 2400,
        }
    )
    drawing = _build_linear_span(
        "features:features@side=overlay",
        "depth:depth@side=below,track_index=0",
        "gc:dinucleotide_content@side=below",
        _span_slot("anchor_slot=features,to_slot=depth,layer=underlay,z=-1"),
        record_count=2,
        depth_track_tables=[[depth], [None]],
        depth_window=100,
        depth_step=100,
        layout=LinearMultiRecordOptions(multi_record_positions=("#1@1", "#2@2")),
    )

    with_depth = _linear_bands(drawing, 0)
    without_depth = _linear_bands(drawing, 1)
    assert with_depth["drawn"] == pytest.approx((with_depth["features"][0], with_depth["depth"][1]))
    # A Depth slot without data for a record has no band there: the span keeps to the other slot.
    assert without_depth["depth"] is None
    assert without_depth["drawn"] == pytest.approx(without_depth["features"])


def test_linear_annotation_clip_policy_creates_clip_path() -> None:
    record = SeqRecord(Seq("A" * 1000), id="r1", name="r1")
    annotations = (
        RegionAnnotation("one", CoordinateSpan(None, 100, 300), lane=0),
        RegionAnnotation("two", CoordinateSpan(None, 400, 600), lane=1),
    )
    drawing = build_linear_diagram(
        [record],
        options=LinearDiagramOptions(
            annotations=AnnotationOptions(
                sets=(AnnotationSet("regions", annotations),)
            ),
            tracks=LinearTrackOptions(
                linear_track_slots=(
                    "regions:annotations@set_id=regions,h=8px,overflow=clip",
                ),
                linear_track_axis_index=1,
            ),
        ),
    )
    svg = drawing.tostring()
    assert 'data-gbdraw-annotation-clipped="true"' in svg
    root = ElementTree.fromstring(svg)
    clipped_group = next(
        element
        for element in root.iter()
        if element.attrib.get("data-gbdraw-annotation-clipped") == "true"
    )
    clip_reference = clipped_group.attrib["clip-path"]
    assert clip_reference.startswith("url(#annotation_clip_")
    clip_id = clip_reference.removeprefix("url(#").removesuffix(")")
    assert any(element.attrib.get("id") == clip_id for element in root.iter())


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
        options=LinearDiagramOptions(
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
        cfg=apply_config_overrides(
            None,
            {
                "labels.linear.scope": "none",
                "canvas.show_gc": False,
                "canvas.show_skew": False,
            },
        ),
        annotation_options=annotations,
        layout=LinearMultiRecordOptions(
            multi_record_positions=("#1@1", "#2@1"),
        ),
        legend="none",
    )
    slots = {
        slot["slotId"]: slot
        for slot in drawing._gbdraw_track_slot_geometry["records"][0]["slots"]
    }

    assert slots["annotations_1"]["heightPx"] == pytest.approx(35.0)
