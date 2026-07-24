from __future__ import annotations

import xml.etree.ElementTree as ET

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from gbdraw.api.diagram import (
    assemble_circular_diagram_from_record,
    assemble_linear_diagram_from_records,
)
from gbdraw.exceptions import ValidationError
from gbdraw.features.colors import compute_feature_hash
from gbdraw.web_support.feature_metadata import extract_features_from_records_payload


SVG_NS = {"svg": "http://www.w3.org/2000/svg"}
PRIVATE_SLOT_ID = "__gbdraw_auto_feature_underlay_slot__"


def _record(*, compound_repeat: bool = False, include_cds: bool = True) -> SeqRecord:
    record = SeqRecord(Seq("A" * 400), id="record_1", name="record_1")
    repeat_location = (
        CompoundLocation(
            [
                FeatureLocation(20, 80, strand=1),
                FeatureLocation(140, 180, strand=1),
            ]
        )
        if compound_repeat
        else FeatureLocation(20, 180, strand=1)
    )
    record.features = [
        SeqFeature(
            repeat_location,
            type="repeat_region",
            qualifiers={"rpt_family": ["family-a"], "rpt_type": ["direct"]},
        )
    ]
    if include_cds:
        record.features.append(
            SeqFeature(
                FeatureLocation(50, 120, strand=1),
                type="CDS",
                qualifiers={"product": ["overlapping protein"]},
            )
        )
    return record


def _underlay_shapes(svg_text: str) -> list[ET.Element]:
    root = ET.fromstring(svg_text)
    return [
        element
        for element in root.findall(".//*")
        if element.attrib.get("data-gbdraw-auto-feature-underlay") == "true"
    ]


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_default_repeat_underlay_has_feature_identity_and_paints_first(mode: str) -> None:
    record = _record()
    repeat = record.features[0]
    feature_id = compute_feature_hash(repeat, record_id=record.id)
    drawing = (
        assemble_circular_diagram_from_record(
            record,
            selected_features_set=["repeat_region", "CDS"],
            legend="none",
        )
        if mode == "circular"
        else assemble_linear_diagram_from_records(
            [record],
            selected_features_set=["repeat_region", "CDS"],
            legend="none",
        )
    )
    svg = drawing.tostring()
    shapes = _underlay_shapes(svg)

    assert len(shapes) == 1
    assert shapes[0].attrib["id"] == feature_id
    assert shapes[0].attrib["data-gbdraw-feature-id"] == feature_id
    assert shapes[0].attrib["data-gbdraw-feature-part"] == "block"
    assert shapes[0].attrib["data-gbdraw-record-id"] == record.id
    assert shapes[0].attrib["data-gbdraw-record-index"] == "0"
    assert shapes[0].attrib["fill"] == "#808080"
    assert shapes[0].attrib["fill-opacity"] == "0.2"
    assert shapes[0].attrib["stroke"] == "none"
    assert shapes[0].attrib["stroke-width"] == "0"
    assert "data-gbdraw-annotation-id" not in svg
    assert svg.index(f'id="gbdraw-annotation-track-{PRIVATE_SLOT_ID}-1"') < svg.index(
        f'id="{record.id}"'
    )
    if mode == "linear":
        assert shapes[0].attrib["data-gbdraw-stable-feature-id"] == feature_id

    metadata = extract_features_from_records_payload(
        [record], selected_features=["repeat_region", "CDS"]
    )
    repeat_metadata = next(
        feature for feature in metadata["features"] if feature["svg_id"] == feature_id
    )
    assert repeat_metadata["type"] == "repeat_region"
    assert repeat_metadata["qualifiers"]["rpt_family"] == ["family-a"]
    assert repeat_metadata["qualifiers"]["rpt_type"] == ["direct"]


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_compound_underlay_preserves_blocks_without_connector(mode: str) -> None:
    record = _record(compound_repeat=True, include_cds=False)
    repeat_id = compute_feature_hash(record.features[0], record_id=record.id)
    drawing = (
        assemble_circular_diagram_from_record(
            record,
            selected_features_set=["repeat_region"],
            legend="none",
        )
        if mode == "circular"
        else assemble_linear_diagram_from_records(
            [record],
            selected_features_set=["repeat_region"],
            legend="none",
        )
    )
    shapes = _underlay_shapes(drawing.tostring())

    assert [shape.attrib["id"] for shape in shapes] == [
        f"{repeat_id}__part1",
        f"{repeat_id}__part2",
    ]
    assert {shape.attrib["data-gbdraw-feature-part"] for shape in shapes} == {"block"}
    if mode == "linear":
        first_x = float(shapes[0].attrib["x"])
        first_end = first_x + float(shapes[0].attrib["width"])
        second_x = float(shapes[1].attrib["x"])
        assert first_end < second_x


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_underlay_uses_actual_custom_feature_slot(mode: str) -> None:
    record = _record(include_cds=False)
    drawing = (
        assemble_circular_diagram_from_record(
            record,
            selected_features_set=["repeat_region"],
            circular_track_slots=["genes:features@side=overlay,z=5"],
            legend="none",
        )
        if mode == "circular"
        else assemble_linear_diagram_from_records(
            [record],
            selected_features_set=["repeat_region"],
            linear_track_slots=["genes:features@side=overlay,z=5"],
            legend="none",
        )
    )
    geometry = drawing._gbdraw_track_slot_geometry
    slots = geometry["records"][0]["slots"]

    assert [slot["slotId"] for slot in slots] == [PRIVATE_SLOT_ID, "genes"]
    assert len(_underlay_shapes(drawing.tostring())) == 1
    if mode == "circular":
        assert slots[0]["widthPx"] == pytest.approx(slots[1]["widthPx"])


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_visible_underlay_requires_enabled_feature_slot(mode: str) -> None:
    record = _record(include_cds=False)
    with pytest.raises(ValidationError, match="enabled features track slot"):
        if mode == "circular":
            assemble_circular_diagram_from_record(
                record,
                selected_features_set=["repeat_region"],
                circular_track_slots=["axis:ticks"],
                legend="none",
            )
        else:
            assemble_linear_diagram_from_records(
                [record],
                selected_features_set=["repeat_region"],
                linear_track_slots=["gap:spacer@side=above"],
                legend="none",
            )


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_explicit_repeat_rectangle_disables_private_underlay(mode: str) -> None:
    record = _record(include_cds=False)
    drawing = (
        assemble_circular_diagram_from_record(
            record,
            selected_features_set=["repeat_region"],
            feature_shapes={"repeat_region": "rectangle"},
            legend="none",
        )
        if mode == "circular"
        else assemble_linear_diagram_from_records(
            [record],
            selected_features_set=["repeat_region"],
            feature_shapes={"repeat_region": "rectangle"},
            legend="none",
        )
    )
    svg = drawing.tostring()

    assert PRIVATE_SLOT_ID not in svg
    assert not _underlay_shapes(svg)
    assert compute_feature_hash(record.features[0], record_id=record.id) in svg
