import re
from pathlib import Path

import pandas as pd
from Bio import SeqIO
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
from gbdraw.core.text import dominant_baseline_shift_em, get_font_vertical_metrics, text_ink
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


def _label_room(text: str, *, stroke_width: float = 1.5, label_offset: float = 4.0, font_size: float = 10.0) -> float:
    """The room a highlight's label takes beyond its lane: half the stroke, label_offset, the letters,
    and label_offset again beyond them (OV-844, review 1170 F1)."""

    ink = text_ink(text, "Arial")
    assert ink is not None
    return 0.5 * stroke_width + 2.0 * label_offset + (ink.top - ink.bottom) * font_size


@pytest.mark.parametrize(
    ("show_labels", "label", "room"),
    [
        pytest.param(False, "Region", 0.0, id="labels-off"),
        # Above an underlay, as for a span: the label's letters and their gap to the tint.
        pytest.param(True, "Region", _label_room("Region"), id="labelled"),
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


_MARKS = ("highlight", "band", "bracket", "line")


@pytest.mark.parametrize("layer", ["underlay", "foreground"])
@pytest.mark.parametrize("mark", _MARKS)
def test_linear_cover_anchor_without_lanes_keeps_the_label_rule_of_its_mark(mark: str, layer: str) -> None:
    def build(cover_anchor: bool):
        return _build_linear_span(
            *_FEATURES_BELOW_STACK,
            "regions:annotations@side=overlay,set_id=regions,anchor_slot=features,"
            f"cover_anchor={str(cover_anchor).lower()},layer={layer},z={-1 if layer == 'underlay' else 1}",
            features=(),
            mark=mark,
        )

    covering, plain = build(True), build(False)
    (label,) = _span_label_marks(covering)

    # Without lanes there is nothing to span (OV-1158): the overlay keeps its own band, and its label
    # keeps label_offset from its own painted mark, as on an overlay that covers nothing.
    assert _ink_gap_to_mark(label) == pytest.approx(4.0, abs=1e-6)
    assert label["ink"] == pytest.approx(_span_label_marks(plain)[0]["ink"], abs=1e-9)
    assert _record_body(covering) == pytest.approx(_record_body(plain), abs=1e-9)


@pytest.mark.parametrize("layer", ["underlay", "foreground"])
@pytest.mark.parametrize("mark", _MARKS)
def test_linear_to_slot_span_without_bodies_keeps_the_label_rule_of_its_mark(mark: str, layer: str) -> None:
    depth = pd.DataFrame({"reference_name": ["rec1"] * 2400, "position": list(range(1, 2401)), "depth": [20.0] * 2400})
    drawing = _build_linear_span(
        "features:features@side=overlay",
        "depth:depth@side=below,track_index=0",
        "gc:dinucleotide_content@side=below",
        "regions:annotations@side=overlay,set_id=regions,anchor_slot=features,to_slot=depth,"
        f"layer={layer},z={-1 if layer == 'underlay' else 1}",
        record_count=2,
        features=(),
        mark=mark,
        depth_track_tables=[[depth], [None]],
        depth_window=100,
        depth_step=100,
        layout=LinearMultiRecordOptions(multi_record_positions=("#1@1", "#2@2")),
    )
    labels = {label["text"]: label for label in _span_label_marks(drawing)}
    bands = _linear_bands(drawing, 0)

    # Record 1 has a Depth plot to span: its label keeps label_offset beyond it (OV-1156).
    box = _label_box(drawing, record_index=0)
    assert bands["depth"] is not None
    gap = bands["depth"][0] - box[1] if layer == "underlay" else box[0] - bands["depth"][1]
    assert gap == pytest.approx(4.0, abs=1e-6)
    # Record 2 has neither lanes nor depth: no span (OV-1158), the label follows its own mark.
    assert _linear_bands(drawing, 1)["depth"] is None
    assert _ink_gap_to_mark(labels["Span 2"]) == pytest.approx(4.0, abs=1e-6)


def _span_slot(
    spec: str = "anchor_slot=features,to_slot=gc,layer=underlay,z=-1",
    *,
    labels: bool = False,
    padding_px: float | None = 0,
) -> str:
    """An annotation overlay slot; ``padding_px=None`` keeps the default padding (2)."""

    padding = "" if padding_px is None else f"padding_px={padding_px},"
    return (
        f"regions:annotations@side=overlay,set_id=regions,{padding}"
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
    label_text: bool = True,
    mark: str = "highlight",
    style: RegionAnnotationStyle | None = None,
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
                        label=f"Span {index}" if label_text else None,
                        mark=mark,
                        **({} if style is None else {"style": style}),
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
    rect = next((element for element in span_group.iter() if element.tag.endswith("rect")), None)
    if rect is not None:
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


def _label_box(drawing, slot_id: str = "regions", record_index: int = 0) -> tuple[float, float]:
    """The vertical box of a slot's drawn labels' painted letters from the record axis."""

    record = drawing._gbdraw_track_slot_geometry["records"][record_index]
    root = ElementTree.fromstring(drawing.tostring())
    parents = {child: parent for parent in root.iter() for child in parent}
    axis_y = float(record["axisYpx"])
    boxes = [
        (ink[1] - axis_y, ink[3] - axis_y)
        for group in root.iter()
        if group.attrib.get("data-gbdraw-annotation-track-id") == slot_id
        and group.attrib.get("data-gbdraw-record-index") == str(record_index)
        for text in group.iter()
        if text.tag.endswith("text") and (text.text or "").strip()
        for ink in (_text_ink_box(text, parents),)
    ]
    assert boxes, f"no {slot_id} labels drawn on record {record_index}"
    return min(top for top, _ in boxes), max(bottom for _, bottom in boxes)


def _assert_clear(box: tuple[float, float], band: tuple[float, float], name: str, gap: float = 4.0) -> None:
    """The letters keep at least ``gap`` (label_offset) from ``band``; touching fails."""

    assert box[1] <= band[0] - gap + 1e-6 or box[0] >= band[1] + gap - 1e-6, (
        f"labels {box} within {gap} of {name} {band}"
    )


def _record_body(drawing, record_index: int = 0) -> tuple[float, float]:
    band = drawing._gbdraw_track_slot_geometry["records"][record_index]["recordBodyBand"]
    return band["topPx"], band["bottomPx"]


_FEATURES_BELOW_STACK = (
    "features:features@side=below",
    "gc:dinucleotide_content@side=below",
    "skew:dinucleotide_skew@side=below",
)

_LABELLED_SPANS = [
    # Review F4: labels below the lanes, where GC starts.
    pytest.param(
        _FEATURES_BELOW_STACK,
        "anchor_slot=features,cover_anchor=true,layer=foreground,z=1",
        ("features",),
        id="cover-anchor-labels-before-gc",
    ),
    # Labels above the lanes, where the axis line is.
    pytest.param(
        _FEATURES_BELOW_STACK,
        "anchor_slot=features,cover_anchor=true,layer=underlay,z=-1",
        (),
        id="cover-anchor-labels-before-axis",
    ),
    pytest.param(
        _SPAN_STACK,
        "anchor_slot=features,to_slot=gc,layer=foreground,z=1",
        ("features", "gc"),
        id="to-slot-labels-before-skew",
    ),
    pytest.param(
        _SPAN_STACK,
        "anchor_slot=gc,to_slot=skew,layer=underlay,z=-1",
        ("features",),
        id="to-slot-labels-before-lanes",
    ),
]


@pytest.mark.parametrize(("stack", "spec", "unmoved"), _LABELLED_SPANS)
def test_linear_span_labels_keep_clear_of_the_next_slot_and_the_axis(
    stack: tuple[str, ...], spec: str, unmoved: tuple[str, ...]
) -> None:
    plain = _build_linear_span(*stack, _span_slot(spec))
    labelled = _build_linear_span(*stack, _span_slot(spec, labels=True))
    box = _label_box(labelled)
    bands = _linear_bands(labelled)
    plain_bands = _linear_bands(plain)

    # The next slot, or the axis line, moves out of the labels (OV-843).
    covered = {
        value for key, value in (part.split("=") for part in spec.split(",")) if key in {"anchor_slot", "to_slot"}
    }
    for slot_id, band in bands.items():
        if slot_id not in {"regions", "drawn", *covered} and band is not None:
            _assert_clear(box, band, slot_id)
    if "features:features@side=below" in stack:
        # Nothing lies above the axis line there without labels: the record body starts at the line.
        axis_top = _record_body(plain)[0]
        _assert_clear(box, (axis_top, -axis_top), "axis")
    # Only on the side with label text; the span still covers its slots.
    for slot_id in unmoved:
        assert bands[slot_id] == pytest.approx(plain_bands[slot_id])
    assert bands["drawn"] == pytest.approx(bands["regions"])
    assert _record_body(labelled)[0] <= box[0] and box[1] <= _record_body(labelled)[1]


@pytest.mark.parametrize(("stack", "spec", "unmoved"), _LABELLED_SPANS)
def test_linear_span_without_label_text_keeps_no_label_room(
    stack: tuple[str, ...], spec: str, unmoved: tuple[str, ...]
) -> None:
    del unmoved
    labels_off = _build_linear_span(*stack, _span_slot(spec))
    no_text = _build_linear_span(*stack, _span_slot(spec, labels=True), label_text=False)

    assert _linear_bands(no_text) == _linear_bands(labels_off)
    assert _record_body(no_text) == pytest.approx(_record_body(labels_off))


def test_linear_plain_overlay_labels_keep_their_room_on_the_anchor() -> None:
    spec = "anchor_slot=gc,layer=foreground,z=1"
    plain = _linear_bands(_build_linear_span(*_SPAN_STACK, _span_slot(spec)))
    drawing = _build_linear_span(*_SPAN_STACK, _span_slot(spec, labels=True))
    labelled = _linear_bands(drawing)

    # An overlay that spans nothing reserves its labels on its anchor, as before: GC stays, skew moves.
    assert labelled["gc"] == pytest.approx(plain["gc"])
    assert labelled["skew"][0] > plain["skew"][0]
    _assert_clear(_label_box(drawing), labelled["skew"], "skew")


def test_linear_highlight_labels_keep_off_the_axis_in_the_below_layout() -> None:
    record = SeqIO.read(Path(__file__).parent / "test_inputs" / "AP027078.gb", "genbank")

    def build(label: str | None):
        crispr = RegionAnnotation("crispr", CoordinateSpan(None, 20000, 40000), label=label, mark="highlight")
        return assemble_linear_diagram_from_records(
            [record],
            cfg=apply_config_overrides(
                None, {"labels.linear.scope": "none", "canvas.linear.track_layout": "below"}
            ),
            legend="none",
            selected_features_set=["CDS", "rRNA", "tRNA"],
            annotation_options=AnnotationOptions(sets=(AnnotationSet("regions", (crispr,)),)),
        )

    # Without label text nothing lies above the axis line, so the record body starts at the line.
    axis_top = _record_body(build(None))[0]
    box = _label_box(build("CRISPR-Cas"), slot_id="annotations_1")

    # The label sits below the axis line with label_offset (4) between them, and the lanes move only
    # that far (default CLI path, OV-843; far side of OV-844).
    assert box[0] == pytest.approx(-axis_top + 4.0)


@pytest.mark.parametrize(
    ("spec", "moved"),
    [
        pytest.param("anchor_slot=features,to_slot=gc,layer=underlay,z=-1", (), id="ends-at-record-top"),
        pytest.param("anchor_slot=gc,to_slot=skew,layer=underlay,z=-1", ("gc", "skew"), id="ends-inside-record"),
    ],
)
def test_linear_annotation_to_slot_labels_keep_room_beyond_the_outermost_slot(
    spec: str, moved: tuple[str, ...]
) -> None:
    def build(labels: bool):
        drawing = _build_linear_span(*_SPAN_STACK, _span_slot(spec, labels=labels), record_count=2)
        return drawing, [_linear_bands(drawing, index) for index in range(2)]

    plain, plain_bands = build(False)
    labelled, labelled_bands = build(True)

    for index in range(2):
        box = _label_box(labelled, record_index=index)
        plain_top, plain_bottom = _record_body(plain, index)
        top, bottom = _record_body(labelled, index)
        shift = labelled_bands[index]["skew"][0] - plain_bands[index]["skew"][0]
        if moved:
            # The labels lie between the lanes and GC: the slots below them move out, the lanes stay.
            assert 0.0 < shift <= _label_room("Span 1")
            for slot_id in moved:
                assert labelled_bands[index][slot_id] == pytest.approx(
                    (plain_bands[index][slot_id][0] + shift, plain_bands[index][slot_id][1] + shift)
                )
            assert labelled_bands[index]["features"] == pytest.approx(plain_bands[index]["features"])
            _assert_clear(box, labelled_bands[index]["features"], "features")
            assert labelled_bands[index]["drawn"] == pytest.approx(
                (labelled_bands[index]["gc"][0], labelled_bands[index]["skew"][1])
            )
            assert (top, bottom) == pytest.approx((plain_top, plain_bottom + shift))
        else:
            # At the record top no slot moves; the label room beyond the lane edge (padding_px 0)
            # joins the record body above an underlay span.
            assert labelled_bands[index].keys() == plain_bands[index].keys()
            for slot_id, band in plain_bands[index].items():
                # Float noise only: the record moves by the room, the bands come back to the same place.
                assert labelled_bands[index][slot_id] == (None if band is None else pytest.approx(band, abs=1e-9))
            assert (top, bottom) == pytest.approx((plain_top - _label_room("Span 1"), plain_bottom))


def _translation(element, parents) -> tuple[float, float]:
    x = y = 0.0
    while element is not None:
        for dx, dy in re.findall(r"translate\(([^,)]+),\s*([^)]+)\)", element.attrib.get("transform", "")):
            x, y = x + float(dx), y + float(dy)
        element = parents.get(element)
    return x, y


def _text_ink_box(text, parents) -> tuple[float, float, float, float]:
    """The box of a drawn text's painted letters (min x, min y, max x, max y) in canvas units."""

    family = text.attrib["font-family"].split(",")[0].strip("'\" ")
    size = float(text.attrib.get("font-size", 10.0))
    ink = text_ink(text.text or "", family)
    metrics = get_font_vertical_metrics(family)
    assert ink is not None and metrics is not None
    shift = dominant_baseline_shift_em(text.attrib.get("dominant-baseline", "auto"), metrics)
    assert shift is not None
    dx, dy = _translation(text, parents)
    share = {"start": 0.0, "middle": 0.5, "end": 1.0}[text.attrib.get("text-anchor", "start")]
    origin = dx + float(text.attrib["x"]) - share * ink.advance * size
    baseline = dy + float(text.attrib["y"]) + shift * size
    return (
        origin + ink.left * size,
        baseline - ink.top * size,
        origin + ink.right * size,
        baseline - ink.bottom * size,
    )


def _span_label_marks(drawing) -> list[dict[str, object]]:
    """Each drawn annotation label: its ink box, its mark's painted band and x range, and the feature
    labels' ink boxes and leader lines of the figure, all in canvas units."""

    root = ElementTree.fromstring(drawing.tostring())
    parents = {child: parent for parent in root.iter() for child in parent}

    def tag(element) -> str:
        return element.tag.rsplit("}", 1)[-1]

    feature_inks = [
        (text.text, _text_ink_box(text, parents))
        for text in root.iter()
        if tag(text) == "text" and text.attrib.get("data-label-feature-id")
    ]
    leaders = []
    for line in root.iter():
        if tag(line) == "line" and line.attrib.get("data-label-feature-id"):
            dx, dy = _translation(line, parents)
            leaders.append(
                (
                    dx + float(line.attrib["x1"]),
                    dy + float(line.attrib["y1"]),
                    dx + float(line.attrib["x2"]),
                    dy + float(line.attrib["y2"]),
                )
            )
    labels: list[dict[str, object]] = []
    for group in root.iter():
        if not group.attrib.get("data-gbdraw-annotation-label"):
            continue
        dx, dy = _translation(group, parents)
        tops: list[float] = []
        bottoms: list[float] = []
        spans: list[tuple[float, float]] = []
        for element in group:
            name = tag(element)
            stroked = element.attrib.get("stroke", "none") != "none"
            half_stroke = 0.5 * float(element.attrib.get("stroke-width", 0.0)) if stroked else 0.0
            if name == "rect":
                x, y = float(element.attrib["x"]), float(element.attrib["y"])
                tops.append(dy + y - half_stroke)
                bottoms.append(dy + y + float(element.attrib["height"]) + half_stroke)
                spans.append((dx + x, dx + x + float(element.attrib["width"])))
            elif name == "line":
                x1, x2 = float(element.attrib["x1"]), float(element.attrib["x2"])
                y1, y2 = float(element.attrib["y1"]), float(element.attrib["y2"])
                # A horizontal stroke is half its width thick each way; a vertical butt-capped one ends at its ends.
                pad = half_stroke if y1 == y2 or element.attrib.get("stroke-linecap") in {"round", "square"} else 0.0
                tops.append(dy + min(y1, y2) - pad)
                bottoms.append(dy + max(y1, y2) + pad)
                if y1 == y2:
                    spans.append((dx + min(x1, x2), dx + max(x1, x2)))
            elif name == "path":
                ys = [float(value) for value in re.findall(r"[-0-9.eE]+,([-0-9.eE]+)", element.attrib["d"])]
                tops.append(dy + min(ys) - half_stroke)
                bottoms.append(dy + max(ys) + half_stroke)
        for text in (element for element in group if tag(element) == "text"):
            labels.append(
                {
                    "text": text.text,
                    "ink": _text_ink_box(text, parents),
                    "mark": (min(tops), max(bottoms)),
                    "spans": sorted(spans),
                    "span": (min(start for start, _ in spans), max(end for _, end in spans)),
                    "feature_inks": feature_inks,
                    "leaders": leaders,
                }
            )
    assert labels, "no annotation labels drawn"
    return labels


def _ink_gap_to_mark(label: dict[str, object]) -> float:
    left, top, right, bottom = label["ink"]
    mark_top, mark_bottom = label["mark"]
    return mark_top - bottom if bottom <= 0.5 * (mark_top + mark_bottom) else top - mark_bottom


@pytest.mark.parametrize("text", ["rRNA genes", "IRb"])
@pytest.mark.parametrize(
    ("track_layout", "slots"),
    [
        # Labels above an underlay highlight: the descender side faces the tint.
        pytest.param("above", (), id="highlight-above"),
        # Labels above the tint, toward the axis line.
        pytest.param("below", (), id="highlight-below"),
        # Labels below a foreground span: the cap side faces the tint.
        pytest.param(
            "middle",
            (*_SPAN_STACK, _span_slot("anchor_slot=features,to_slot=gc,layer=foreground,z=1", labels=True)),
            id="to-slot-foreground",
        ),
    ],
)
def test_linear_span_label_ink_keeps_label_offset_clear_of_its_mark(
    track_layout: str, slots: tuple[str, ...], text: str
) -> None:
    style = RegionAnnotationStyle(stroke_width=0.6, label_offset=4.0)
    record = SeqRecord(Seq("ACGT" * 600), id="rec1", name="rec1")
    record.features = [SeqFeature(FeatureLocation(120, 250, strand=1), type="CDS")]
    tint = RegionAnnotation("tint", CoordinateSpan(None, 100, 1300), label=text, mark="highlight", style=style)
    drawing = assemble_linear_diagram_from_records(
        [record],
        cfg=apply_config_overrides(
            None, {"labels.linear.scope": "none", "canvas.linear.track_layout": track_layout}
        ),
        legend="none",
        selected_features_set=["CDS"],
        annotation_options=AnnotationOptions(sets=(AnnotationSet("regions", (tint,)),)),
        **({"linear_track_slots": list(slots), "linear_track_axis_index": 0} if slots else {}),
    )

    (label,) = _span_label_marks(drawing)
    # The painted letters keep label_offset clear of the tint's stroke (OV-844), descenders included.
    assert _ink_gap_to_mark(label) == pytest.approx(4.0, abs=1e-6)


def _leader_zone_drawing(*, feature_names: tuple[str, ...], track_layout: str = "above"):
    """A tint over features whose labels sit in rows above the lanes, with leaders under the tint label."""

    record = SeqRecord(Seq("ACGT" * 1500), id="rec1", name="rec1")
    record.features = [
        SeqFeature(
            FeatureLocation(2800 + 90 * index, 2830 + 90 * index, strand=1),
            type="CDS",
            qualifiers={"product": [name]},
        )
        for index, name in enumerate(feature_names)
    ]
    tint = RegionAnnotation("tint", CoordinateSpan(None, 1000, 5000), label="Region of interest", mark="highlight")
    return assemble_linear_diagram_from_records(
        [record],
        cfg=apply_config_overrides(
            None, {"labels.linear.scope": "all", "canvas.linear.track_layout": track_layout}
        ),
        legend="none",
        selected_features_set=["CDS"],
        annotation_options=AnnotationOptions(sets=(AnnotationSet("regions", (tint,)),)),
    )


def test_linear_span_label_moves_along_its_span_off_feature_leaders_and_labels() -> None:
    (label,) = _span_label_marks(
        _leader_zone_drawing(feature_names=("hypothetical protein", "ribosomal protein S7", "tRNA-Leu"))
    )
    left, top, right, bottom = label["ink"]

    # No leader runs through the label's letters and no feature label overlaps them (OV-845).
    assert not [
        leader
        for leader in label["leaders"]
        if min(leader[0], leader[2]) <= right
        and max(leader[0], leader[2]) >= left
        and min(leader[1], leader[3]) < bottom
        and max(leader[1], leader[3]) > top
    ]
    assert not [
        name
        for name, box in label["feature_inks"]
        if box[0] < right and box[2] > left and box[1] < bottom and box[3] > top
    ]
    # It stays over its own span.
    assert label["span"][0] <= left and right <= label["span"][1]


def test_linear_span_label_without_conflicts_stays_at_its_label_position() -> None:
    (label,) = _span_label_marks(_leader_zone_drawing(feature_names=()))
    left, _top, right, _bottom = label["ink"]

    # Centred on its span, as before (LP-D41: placements without a defect stay).
    assert 0.5 * (left + right) == pytest.approx(0.5 * sum(label["span"]), abs=0.5)


@pytest.mark.parametrize("track_layout", ["above", "below"])
@pytest.mark.parametrize(
    "style",
    [
        pytest.param(RegionAnnotationStyle(), id="bracket-tick"),
        pytest.param(RegionAnnotationStyle(line_cap="none"), id="bracket-no-cap"),
        pytest.param(RegionAnnotationStyle(stroke_width=3.0), id="bracket-thick"),
    ],
)
@pytest.mark.parametrize("mark", ["highlight", "band", "bracket", "line"])
def test_linear_annotation_label_keeps_label_offset_from_the_painted_mark(
    mark: str, style: RegionAnnotationStyle, track_layout: str
) -> None:
    """Every mark kind: label_offset from its painted edge (rect and stroke, caps, line stroke) to the letters."""

    record = SeqRecord(Seq("ACGT" * 600), id="rec1", name="rec1")
    record.features = [SeqFeature(FeatureLocation(120, 250, strand=1), type="CDS")]
    region = RegionAnnotation("region", CoordinateSpan(None, 100, 1300), label="rRNA genes", mark=mark, style=style)
    drawing = assemble_linear_diagram_from_records(
        [record],
        cfg=apply_config_overrides(
            None, {"labels.linear.scope": "none", "canvas.linear.track_layout": track_layout}
        ),
        legend="none",
        selected_features_set=["CDS"],
        annotation_options=AnnotationOptions(sets=(AnnotationSet("regions", (region,)),)),
    )

    (label,) = _span_label_marks(drawing)
    assert _ink_gap_to_mark(label) == pytest.approx(style.label_offset, abs=1e-6)


def _slot_paint(drawing, slot_id: str, record_index: int = 0) -> tuple[float, float]:
    slot = next(
        slot
        for slot in drawing._gbdraw_track_slot_geometry["records"][record_index]["slots"]
        if slot["slotId"] == slot_id
    )
    return slot["paintBand"]["topPx"], slot["paintBand"]["bottomPx"]


def _real_record_tint(file_name: str, start: int, end: int, label: str | None, track_layout: str):
    record = SeqIO.read(Path(__file__).parent / "test_inputs" / file_name, "genbank")
    tint = RegionAnnotation("tint", CoordinateSpan(None, start, end), label=label, mark="highlight")
    return assemble_linear_diagram_from_records(
        [record],
        cfg=apply_config_overrides(
            None, {"labels.linear.scope": "none", "canvas.linear.track_layout": track_layout}
        ),
        legend="none",
        selected_features_set=["CDS", "rRNA", "tRNA"],
        annotation_options=AnnotationOptions(sets=(AnnotationSet("regions", (tint,)),)),
    )


@pytest.mark.parametrize("track_layout", ["above", "below"])
@pytest.mark.parametrize(
    ("file_name", "start", "end", "label"),
    [
        pytest.param("NC_001879.gbk", 86687, 112029, "IRb", id="tobacco-IRb"),
        pytest.param("AP027078.gb", 20000, 40000, "CRISPR-Cas", id="hepatoplasma-CRISPR"),
    ],
)
def test_linear_highlight_label_keeps_label_offset_on_both_sides(
    file_name: str, start: int, end: int, label: str, track_layout: str
) -> None:
    labelled = _real_record_tint(file_name, start, end, label, track_layout)
    (drawn,) = _span_label_marks(labelled)
    box = _label_box(labelled, slot_id="annotations_1")

    # Own side: the tint's stroke (OV-844).
    assert _ink_gap_to_mark(drawn) == pytest.approx(4.0, abs=1e-6)
    # Far side (review 1170 F1): the axis line in `below` (the record body of the unlabelled figure starts
    # at its lower edge, mirrored), the record's top edge in `above`.
    if track_layout == "below":
        axis_bottom = -_record_body(_real_record_tint(file_name, start, end, None, track_layout))[0]
        assert box[0] - axis_bottom == pytest.approx(4.0, abs=1e-6)
    else:
        assert box[0] - _record_body(labelled)[0] == pytest.approx(4.0, abs=1e-6)


@pytest.mark.parametrize(
    ("stack", "spec", "next_slot"),
    [
        pytest.param(
            _SPAN_STACK, "anchor_slot=features,to_slot=gc,layer=foreground,z=1", "skew", id="to-slot-then-skew"
        ),
        pytest.param(
            _FEATURES_BELOW_STACK,
            "anchor_slot=features,cover_anchor=true,layer=foreground,z=1",
            "gc",
            id="cover-anchor-then-gc",
        ),
    ],
)
def test_linear_span_label_keeps_label_offset_from_the_next_slot(
    stack: tuple[str, ...], spec: str, next_slot: str
) -> None:
    drawing = _build_linear_span(*stack, _span_slot(spec, labels=True))
    (drawn,) = _span_label_marks(drawing)
    box = _label_box(drawing)

    assert _ink_gap_to_mark(drawn) == pytest.approx(4.0, abs=1e-6)
    # The next slot's paint starts label_offset beyond the letters (review 1170 F1).
    assert _slot_paint(drawing, next_slot)[0] - box[1] == pytest.approx(4.0, abs=1e-6)


_COVERING_SPANS = [
    pytest.param(_SPAN_STACK, "anchor_slot=features,to_slot=gc,layer=foreground,z=1", "skew", id="to-slot-below"),
    pytest.param(_SPAN_STACK, "anchor_slot=gc,to_slot=skew,layer=underlay,z=-1", "features", id="to-slot-above"),
    pytest.param(
        _FEATURES_BELOW_STACK, "anchor_slot=features,cover_anchor=true,layer=foreground,z=1", "gc", id="cover-below"
    ),
    pytest.param(
        _FEATURES_BELOW_STACK, "anchor_slot=features,cover_anchor=true,layer=underlay,z=-1", None, id="cover-above"
    ),
]


@pytest.mark.parametrize(
    ("padding_px", "stroke_width"),
    [
        # The default padding holds the default stroke: the letters keep label_offset from the span edge.
        pytest.param(None, 1.5, id="default-padding"),
        # Without padding a highlight's stroke passes the span edge by half its width.
        pytest.param(0, 1.5, id="no-padding"),
        # A stroke wider than twice the padding: only a mark that paints out there passes the edge.
        pytest.param(None, 10.0, id="wide-stroke"),
    ],
)
@pytest.mark.parametrize("mark", _MARKS)
@pytest.mark.parametrize(("stack", "spec", "next_slot"), _COVERING_SPANS)
def test_linear_span_label_stays_outside_every_slot_it_covers(
    stack: tuple[str, ...],
    spec: str,
    next_slot: str | None,
    mark: str,
    padding_px: float | None,
    stroke_width: float,
) -> None:
    def build(labels: bool):
        return _build_linear_span(
            *stack,
            _span_slot(spec, labels=labels, padding_px=padding_px),
            mark=mark,
            style=RegionAnnotationStyle(stroke_width=stroke_width),
        )

    drawing = build(True)
    box = _label_box(drawing)
    bands = _linear_bands(drawing)
    covered = {
        value for key, value in (part.split("=") for part in spec.split(",")) if key in {"anchor_slot", "to_slot"}
    }
    (label,) = _span_label_marks(drawing)
    axis_y = float(drawing._gbdraw_track_slot_geometry["records"][0]["axisYpx"])
    mark_top, mark_bottom = (edge - axis_y for edge in label["mark"])
    span_top, span_bottom = bands["regions"]
    above = box[1] <= span_top

    # Whatever the mark, the letters keep label_offset outside every covered slot (OV-1156) ...
    for slot_id in covered:
        _assert_clear(box, bands[slot_id], slot_id)
    # ... from the span edge, or from the mark's paint where it really passes that edge (review 1176):
    # a highlight's stroke beyond the padding; never a bracket or line at the lane centre.
    edge_gap = span_top - box[1] if above else box[0] - span_bottom
    paint_gap = mark_top - box[1] if above else box[0] - mark_bottom
    assert min(edge_gap, paint_gap) == pytest.approx(4.0, abs=1e-6)
    if mark in {"bracket", "line"}:
        assert edge_gap == pytest.approx(4.0, abs=1e-6)
    if mark == "highlight":
        passes = max(0.0, 0.5 * stroke_width - (2.0 if padding_px is None else padding_px))
        assert edge_gap == pytest.approx(4.0 + passes, abs=1e-6)
    # And label_offset from what lies beyond them.
    if next_slot is not None:
        far = _slot_paint(drawing, next_slot)
        assert min(abs(far[0] - box[1]), abs(box[0] - far[1])) == pytest.approx(4.0, abs=1e-6)
    else:
        # Below the axis line: without label text the record body starts at the line.
        axis_top = _record_body(build(False))[0]
        assert box[0] == pytest.approx(-axis_top + 4.0, abs=1e-6)


def test_linear_hmmtdna_highlight_span_labels_keep_label_offset_from_the_covered_gc_plot() -> None:
    """The figure of review 1174: rRNA genes and Control region below a features-to-GC highlight span."""

    record = SeqIO.read(Path(__file__).parent / "test_inputs" / "HmmtDNA.gbk", "genbank")
    style = RegionAnnotationStyle(stroke_width=0.6)
    regions = (
        RegionAnnotation("rrna", CoordinateSpan(None, 648, 3229), label="rRNA genes", mark="highlight", style=style),
        RegionAnnotation(
            "control", CoordinateSpan(None, 16024, 16569), label="Control region", mark="highlight", style=style
        ),
    )
    drawing = assemble_linear_diagram_from_records(
        [record],
        cfg=apply_config_overrides(None, {"labels.linear.scope": "all"}),
        legend="none",
        selected_features_set=["CDS", "rRNA", "tRNA"],
        annotation_options=AnnotationOptions(sets=(AnnotationSet("regions", regions),)),
        linear_track_slots=[
            *_SPAN_STACK,
            _span_slot("anchor_slot=features,to_slot=gc,layer=foreground,z=5", labels=True, padding_px=None),
        ],
        linear_track_axis_index=0,
    )
    gc_bottom = _linear_bands(drawing)["gc"][1]
    axis_y = float(drawing._gbdraw_track_slot_geometry["records"][0]["axisYpx"])

    # Before OV-1156 the letters kept label_offset from the tint inside the span, 2.3 px from the GC plot.
    for label in _span_label_marks(drawing):
        assert label["ink"][1] - axis_y - gc_bottom == pytest.approx(4.0, abs=1e-6), label["text"]


@pytest.mark.parametrize("mark", ["band", "bracket"])
def test_linear_span_label_place_ignores_the_overlays_own_height(mark: str) -> None:
    spec = "anchor_slot=features,to_slot=gc,layer=foreground,z=1"
    plain = _build_linear_span(*_SPAN_STACK, _span_slot(spec, labels=True), mark=mark)
    tall = _build_linear_span(*_SPAN_STACK, _span_slot(f"{spec},h=150", labels=True), mark=mark)

    # With to_slot the overlay's own h is not read: not for the label, not for its room (OV-1156).
    assert _label_box(tall) == pytest.approx(_label_box(plain), abs=1e-6)
    assert _slot_paint(tall, "skew") == pytest.approx(_slot_paint(plain, "skew"), abs=1e-6)


def _six_kb_drawing(region: RegionAnnotation, features: tuple[tuple[int, int, str], ...]):
    record = SeqRecord(Seq("ACGT" * 1500), id="rec1", name="rec1", annotations={"topology": "circular"})
    record.features = [
        SeqFeature(FeatureLocation(start, end, strand=1), type="CDS", qualifiers={"product": [name]})
        for start, end, name in features
    ]
    return assemble_linear_diagram_from_records(
        [record],
        cfg=apply_config_overrides(None, {"labels.linear.scope": "all", "canvas.linear.track_layout": "above"}),
        legend="none",
        selected_features_set=["CDS"],
        annotation_options=AnnotationOptions(sets=(AnnotationSet("regions", (region,)),)),
    )


def _leaders_through(label: dict[str, object]) -> list[tuple[float, float, float, float]]:
    left, top, right, bottom = label["ink"]
    return [
        leader
        for leader in label["leaders"]
        if min(leader[0], leader[2]) <= right
        and max(leader[0], leader[2]) >= left
        and min(leader[1], leader[3]) < bottom
        and max(leader[1], leader[3]) > top
    ]


def test_linear_label_of_an_origin_wrapping_region_stays_over_its_marks() -> None:
    region = RegionAnnotation(
        "origin", CoordinateSpan(None, 5800, 300, wraps_origin=True), label="Origin region", mark="highlight"
    )
    features = (
        (20, 50, "photosystem II protein D1"),
        (60, 90, "ribosomal protein S7"),
        (5850, 5880, "tRNA-Leu"),
        (5900, 5930, "maturase K"),
    )
    (label,) = _span_label_marks(_six_kb_drawing(region, features))
    centre = 0.5 * (label["ink"][0] + label["ink"][2])

    # Its centre stays over one of its two marks, never over the sequence between them (review 1170 F4).
    assert len(label["spans"]) == 2
    assert any(start <= centre <= end for start, end in label["spans"])


def test_linear_label_wider_than_its_region_moves_off_a_leader() -> None:
    region = RegionAnnotation(
        "narrow", CoordinateSpan(None, 2900, 3200), label="A long label for a short region", mark="highlight"
    )
    (label,) = _span_label_marks(_six_kb_drawing(region, ((2905, 2925, "tRNA-Leu"),)))
    centre = 0.5 * (label["ink"][0] + label["ink"][2])

    # Wider than its mark, it still moves off the leader under its default place, centre over the mark (F5).
    assert label["ink"][2] - label["ink"][0] > label["span"][1] - label["span"][0]
    assert not _leaders_through(label)
    assert label["span"][0] <= centre <= label["span"][1]


def test_linear_label_without_a_clear_place_keeps_its_default_place_off_the_name() -> None:
    region = RegionAnnotation(
        "control", CoordinateSpan(None, 1, 300), label="Mitochondrial control region", mark="highlight"
    )
    drawing = _six_kb_drawing(region, ((190, 210, "tRNA-Leu"),))
    (label,) = _span_label_marks(drawing)
    root = ElementTree.fromstring(drawing.tostring())
    parents = {child: parent for parent in root.iter() for child in parent}
    name = _text_ink_box(next(text for text in root.iter() if text.tag.endswith("text") and text.text == "rec1"), parents)
    left, top, right, bottom = label["ink"]

    # No place with its letters on the bar clears the leader, so the label keeps its default place
    # (OV-1155; its letters start before the bar there) and does not move into the record name (OV-1157).
    # (the ink centre is within half a pixel of the text anchor)
    assert 0.5 * (left + right) == pytest.approx(0.5 * sum(label["span"]), abs=0.5)
    assert right <= name[0] or left >= name[2] or bottom <= name[1] or top >= name[3]


@pytest.mark.parametrize(
    ("x", "segment", "hanging", "placed"),
    [
        # A feature label hangs past the bar start over the label's default letters (-10 to 30): the nearest
        # clear place (x = 15) would put the letters before the bar, so they stop at its start (x = 20).
        pytest.param(10.0, (0.0, 60.0), (-40.0, -8.0), 20.0, id="bar-start"),
        # The same at the bar end (1000): the nearest clear place is x = 985, the letters stop at the end.
        pytest.param(990.0, (940.0, 1000.0), (1008.0, 1040.0), 980.0, id="bar-end"),
    ],
)
def test_linear_moved_label_stops_with_its_letters_at_the_bar_end(
    x: float, segment: tuple[float, float], hanging: tuple[float, float], placed: float
) -> None:
    from gbdraw.diagrams.linear.annotation_labels import _place_along_span
    from gbdraw.layout.spatial import Aabb
    from gbdraw.render.drawers.linear.annotations import LinearAnnotationLabel

    label = LinearAnnotationLabel("a", "Label", x, 0.0, "middle", 10.0, "#000", Aabb(x - 20, -4, x + 20, 4), (segment,))
    feature_label = Aabb(hanging[0], -6.0, hanging[1], 6.0)

    moved = label.moved_to(_place_along_span(label, 0.0, [feature_label], 3.0, 1000.0))
    # Moved off the feature label, and clamped: the letters end exactly at the bar's end (OV-1157).
    assert moved.x == pytest.approx(placed)
    assert (moved.ink.min_x, moved.ink.max_x) == pytest.approx((0.0, 40.0) if x < 500 else (960.0, 1000.0))
    assert moved.ink.min_x >= feature_label.max_x + 3.0 or moved.ink.max_x <= feature_label.min_x - 3.0


def test_linear_label_search_matches_a_full_scan() -> None:
    """The sorted sweep gives the place a scan of every candidate gives (review 1170 F6), letters on the bar."""

    import random

    from gbdraw.diagrams.linear.annotation_labels import _place_along_span
    from gbdraw.layout.spatial import Aabb
    from gbdraw.render.drawers.linear.annotations import LinearAnnotationLabel

    generator = random.Random(1170)
    for _ in range(200):
        segments = tuple(
            sorted(
                (start, start + generator.uniform(5.0, 120.0))
                for start in (generator.uniform(0.0, 900.0) for _ in range(generator.randint(1, 3)))
            )
        )
        x = generator.uniform(segments[0][0], segments[0][1])
        label = LinearAnnotationLabel("a", "Label", x, 0.0, "middle", 10.0, "#000", Aabb(x - 20, -4, x + 20, 4), segments)
        obstacles = [
            Aabb(left, generator.uniform(-12, 2), left + generator.uniform(0.0, 40.0), generator.uniform(3, 12))
            for left in (generator.uniform(-50.0, 1000.0) for _ in range(generator.randint(0, 40)))
        ]
        gap = 3.0
        blocked = [
            (box.min_x - gap - 20, box.max_x + gap + 20)
            for box in obstacles
            if box.min_y < 4 + gap and box.max_y > -4 - gap
        ]

        def conflicts(value: float) -> int:
            return sum(1 for start, end in blocked if start < value < end)

        if conflicts(x) == 0:
            expected = x
        else:
            candidates = {x}
            # The letters (20 px each side of x) stay on a 1000 px bar.
            for low, high in ((max(low, 20.0), min(high, 980.0)) for low, high in segments):
                if low <= high:
                    candidates |= {low, high, *(v for interval in blocked for v in interval if low <= v <= high)}
            expected = min(sorted(candidates), key=lambda value: (conflicts(value), abs(value - x)))
        assert _place_along_span(label, 0.0, obstacles, gap, 1000.0) == pytest.approx(expected, abs=1e-9)


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
