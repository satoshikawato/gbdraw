from __future__ import annotations

import xml.etree.ElementTree as ET

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from gbdraw.annotations.feature_underlays import (
    AUTO_FEATURE_UNDERLAY_FILL,
    merge_feature_underlays,
)
from gbdraw.annotations.models import ResolvedAnnotationBundle
from gbdraw.api import LinearMultiRecordOptions
from gbdraw.api.config import apply_config_overrides
from gbdraw.api.diagram import (
    assemble_circular_diagram_from_record,
    assemble_linear_diagram_from_records,
)
from gbdraw.exceptions import ValidationError
from gbdraw.features.colors import compute_feature_hash
from gbdraw.features.objects import FeatureLocationPart, FeatureObject
from gbdraw.io.colors import load_default_colors
from gbdraw.web_support.feature_metadata import extract_features_from_records_payload


SVG_NS = {"svg": "http://www.w3.org/2000/svg"}
PRIVATE_SLOT_ID = "__gbdraw_auto_feature_underlay_slot__"
_DEFAULT_CFG = apply_config_overrides(None, None)
_LINEAR_CFG = apply_config_overrides(
    None,
    {
        "labels.linear.scope": "none",
        "canvas.show_gc": False,
        "canvas.show_skew": False,
    },
)


def _record(
    *,
    record_id: str = "record_1",
    length: int = 400,
    compound_repeat: bool = False,
    include_cds: bool = True,
) -> SeqRecord:
    record = SeqRecord(Seq("A" * length), id=record_id, name=record_id)
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


def _depth_table(record_id: str, record_length: int) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "reference_name": [record_id] * 4,
            "position": [1, record_length // 3, 2 * record_length // 3, record_length],
            "depth": [10, 30, 20, 40],
        }
    )


def _underlay_shapes(svg_text: str) -> list[ET.Element]:
    root = ET.fromstring(svg_text)
    return [
        element
        for element in root.findall(".//*")
        if element.attrib.get("data-gbdraw-auto-feature-underlay") == "true"
    ]


@pytest.mark.linear
@pytest.mark.parametrize(
    ("axis_index", "slot_specs", "expected_depth_side"),
    [
        (
            0,
            [
                "features:features@side=overlay",
                "depth_3:depth@side=below,track_index=0",
            ],
            "below",
        ),
        (
            1,
            [
                "depth_3:depth@side=above,track_index=0",
                "features:features@side=overlay",
            ],
            "above",
        ),
        (
            2,
            [
                "depth_3:depth@side=above,track_index=0",
                "features:features@side=above",
            ],
            "above",
        ),
    ],
    ids=["axis-start", "axis-middle-incident", "axis-end"],
)
def test_linear_underlay_preserves_public_axis_sides(
    axis_index: int,
    slot_specs: list[str],
    expected_depth_side: str,
) -> None:
    record = _record(include_cds=False)
    drawing = assemble_linear_diagram_from_records(
        [record],
        cfg=_LINEAR_CFG,
        selected_features_set=["repeat_region"],
        depth_track_tables=[[_depth_table(record.id, len(record.seq))]],
        linear_track_slots=slot_specs,
        linear_track_axis_index=axis_index,
        legend="none",
    )
    svg = drawing.tostring()
    slots = drawing._gbdraw_track_slot_geometry["records"][0]["slots"]
    depth = next(slot for slot in slots if slot["slotId"] == "depth_3")
    depth_band = depth["paintBand"]

    assert [slot["slotId"] for slot in slots] == [
        PRIVATE_SLOT_ID,
        *(spec.split(":", 1)[0] for spec in slot_specs),
    ]
    assert depth["side"] == expected_depth_side
    assert depth_band is not None
    if expected_depth_side == "above":
        assert depth_band["bottomPx"] < 0
    else:
        assert depth_band["topPx"] > 0
    assert svg.index(f'data-gbdraw-slot-id="{PRIVATE_SLOT_ID}"') < svg.index(
        'data-gbdraw-slot-id="features"'
    )


@pytest.mark.linear
def test_linear_underlay_with_axis_omitted_keeps_explicit_sides() -> None:
    record = _record(include_cds=False)
    drawing = assemble_linear_diagram_from_records(
        [record],
        cfg=_LINEAR_CFG,
        selected_features_set=["repeat_region"],
        depth_track_tables=[[_depth_table(record.id, len(record.seq))]],
        linear_track_slots=[
            "depth_3:depth@side=above,track_index=0",
            "features:features@side=overlay",
        ],
        legend="none",
    )
    slots = drawing._gbdraw_track_slot_geometry["records"][0]["slots"]

    assert [slot["slotId"] for slot in slots] == [
        PRIVATE_SLOT_ID,
        "depth_3",
        "features",
    ]
    assert next(slot for slot in slots if slot["slotId"] == "depth_3")[
        "side"
    ] == "above"


@pytest.mark.linear
def test_linear_multi_record_underlay_keeps_axis_during_final_relayout() -> None:
    records = [
        _record(record_id="record_1", length=400, include_cds=False),
        _record(record_id="record_2", length=700, include_cds=False),
    ]
    drawing = assemble_linear_diagram_from_records(
        records,
        cfg=_LINEAR_CFG,
        layout=LinearMultiRecordOptions(
            multi_record_positions=("#1@1", "#2@1"),
        ),
        selected_features_set=["repeat_region"],
        depth_track_tables=[
            [_depth_table(record.id, len(record.seq))] for record in records
        ],
        linear_track_slots=[
            "depth_3:depth@side=above,track_index=0",
            "features:features@side=overlay",
        ],
        linear_track_axis_index=1,
        legend="none",
    )
    resolved_records = drawing._gbdraw_track_slot_geometry["records"]

    assert len(resolved_records) == 2
    assert resolved_records[0]["axisYpx"] == pytest.approx(
        resolved_records[1]["axisYpx"]
    )
    for resolved_record in resolved_records:
        assert [slot["slotId"] for slot in resolved_record["slots"]] == [
            PRIVATE_SLOT_ID,
            "depth_3",
            "features",
        ]
        depth = next(
            slot for slot in resolved_record["slots"] if slot["slotId"] == "depth_3"
        )
        assert depth["side"] == "above"
        assert depth["paintBand"]["bottomPx"] < 0


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_default_repeat_underlay_has_feature_identity_and_paints_first(mode: str) -> None:
    record = _record()
    repeat = record.features[0]
    feature_id = compute_feature_hash(repeat, record_id=record.id)
    default_colors = load_default_colors("", "default")
    expected_fill = default_colors.loc[
        default_colors["feature_type"] == "repeat_region", "color"
    ].iloc[0]
    drawing = (
        assemble_circular_diagram_from_record(
            record,
            cfg=_DEFAULT_CFG,
            default_colors=default_colors,
            selected_features_set=["repeat_region", "CDS"],
            legend="none",
        )
        if mode == "circular"
        else assemble_linear_diagram_from_records(
            [record],
            cfg=_DEFAULT_CFG,
            default_colors=default_colors,
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
    assert shapes[0].attrib["fill"] == expected_fill
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
def test_underlays_keep_each_specific_color_rule_fill(mode: str) -> None:
    record = _record(include_cds=False)
    record.features.append(
        SeqFeature(
            FeatureLocation(220, 320, strand=1),
            type="repeat_region",
            qualifiers={"rpt_family": ["family-b"], "rpt_type": ["inverted"]},
        )
    )
    color_table = pd.DataFrame(
        [
            ["repeat_region", "rpt_family", "^family-a$", "#ff00aa", "Family A"],
            ["repeat_region", "rpt_family", "^family-b$", "#00aaff", "Family B"],
        ],
        columns=["feature_type", "qualifier_key", "value", "color", "caption"],
    )
    drawing = (
        assemble_circular_diagram_from_record(
            record,
            cfg=_DEFAULT_CFG,
            color_table=color_table,
            selected_features_set=["repeat_region"],
            legend="none",
        )
        if mode == "circular"
        else assemble_linear_diagram_from_records(
            [record],
            cfg=_DEFAULT_CFG,
            color_table=color_table,
            selected_features_set=["repeat_region"],
            legend="none",
        )
    )
    fills_by_feature_id = {
        shape.attrib["id"]: shape.attrib["fill"]
        for shape in _underlay_shapes(drawing.tostring())
    }

    assert fills_by_feature_id == {
        compute_feature_hash(record.features[0], record_id=record.id): "#ff00aa",
        compute_feature_hash(record.features[1], record_id=record.id): "#00aaff",
    }


def test_empty_feature_color_keeps_defensive_underlay_fill() -> None:
    record = _record(include_cds=False)
    coordinates = [
        FeatureLocationPart("block", "001", "positive", 20, 180, True)
    ]
    feature = FeatureObject(
        feature_id="feature_without_color",
        location=coordinates,
        is_directional=False,
        color="",
        note="",
        label_text="",
        coordinates=coordinates,
        type="repeat_region",
        qualifiers={},
        record_id=record.id,
    )

    merged, _set_id = merge_feature_underlays(
        ResolvedAnnotationBundle(()),
        [[feature]],
        [record],
        mode="linear",
    )

    assert merged.annotations[0].style.fill == AUTO_FEATURE_UNDERLAY_FILL


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_compound_underlay_preserves_blocks_without_connector(mode: str) -> None:
    record = _record(compound_repeat=True, include_cds=False)
    repeat_id = compute_feature_hash(record.features[0], record_id=record.id)
    drawing = (
        assemble_circular_diagram_from_record(
            record,
            cfg=_DEFAULT_CFG,
            selected_features_set=["repeat_region"],
            legend="none",
        )
        if mode == "circular"
        else assemble_linear_diagram_from_records(
            [record],
            cfg=_DEFAULT_CFG,
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
            cfg=_DEFAULT_CFG,
            selected_features_set=["repeat_region"],
            circular_track_slots=["genes:features@side=overlay,z=5"],
            legend="none",
        )
        if mode == "circular"
        else assemble_linear_diagram_from_records(
            [record],
            cfg=_DEFAULT_CFG,
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
                cfg=_DEFAULT_CFG,
                selected_features_set=["repeat_region"],
                circular_track_slots=["axis:ticks"],
                legend="none",
            )
        else:
            assemble_linear_diagram_from_records(
                [record],
                cfg=_DEFAULT_CFG,
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
            cfg=_DEFAULT_CFG,
            selected_features_set=["repeat_region"],
            feature_shapes={"repeat_region": "rectangle"},
            legend="none",
        )
        if mode == "circular"
        else assemble_linear_diagram_from_records(
            [record],
            cfg=_DEFAULT_CFG,
            selected_features_set=["repeat_region"],
            feature_shapes={"repeat_region": "rectangle"},
            legend="none",
        )
    )
    svg = drawing.tostring()

    assert PRIVATE_SLOT_ID not in svg
    assert not _underlay_shapes(svg)
    assert compute_feature_hash(record.features[0], record_id=record.id) in svg
