from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from collections import Counter

import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from gbdraw.annotations import (
    AnnotationOptions,
    AnnotationSet,
    CoordinateSpan,
    HatchStyle,
    RegionAnnotation,
    RegionAnnotationStyle,
)
from gbdraw.api.config import apply_config_overrides
from gbdraw.api.diagram import (
    assemble_circular_diagram_from_record,
    assemble_circular_diagram_from_records,
)
from gbdraw.io.comparisons import COMPARISON_COLUMNS
from gbdraw.io.record_select import parse_record_selector
from gbdraw.render.interactive_svg import enrich_svg

_URL_REFERENCE_RE = re.compile(r"url\(#([^)]+)\)")


def _assert_unique_ids_and_resolved_local_references(svg: str) -> None:
    root = ET.fromstring(svg)
    ids = [
        element.attrib["id"]
        for element in root.iter()
        if element.attrib.get("id")
    ]
    assert not {
        identifier: count
        for identifier, count in Counter(ids).items()
        if count > 1
    }

    references: set[str] = set()
    for element in root.iter():
        for name, value in element.attrib.items():
            references.update(_URL_REFERENCE_RE.findall(value))
            if name.rsplit("}", 1)[-1] == "href" and value.startswith("#"):
                references.add(value[1:])
    assert references <= set(ids)


def _comparison_frame(record_id: str) -> pd.DataFrame:
    row = [
        record_id,
        "comparison",
        90.0,
        100,
        0,
        0,
        10,
        100,
        20,
        110,
        1e-20,
        200,
    ]
    return pd.DataFrame([row], columns=COMPARISON_COLUMNS)


def test_conservation_group_ids_distinguish_lossy_labels_and_fixed_groups() -> None:
    record = SeqRecord(Seq("A" * 1000), id="reference", name="reference")
    comparison = _comparison_frame(record.id)
    labels = ("a b", "a_b", "identity_legend")

    def render() -> str:
        return assemble_circular_diagram_from_record(
            record,
            conservation_dataframes=[comparison, comparison, comparison],
            conservation_reference="query",
            conservation_labels=labels,
            legend="right",
            cfg=apply_config_overrides(None, {
                "labels.circular.scope": "none",
                "canvas.show_gc": False,
                "canvas.show_skew": False,
            }),
        ).tostring()

    first = render()
    second = render()
    _assert_unique_ids_and_resolved_local_references(first)

    def conservation_ids(svg: str) -> dict[str, str]:
        return {
            element.attrib["data-track-label"]: element.attrib["id"]
            for element in ET.fromstring(svg).iter()
            if element.tag.rsplit("}", 1)[-1] == "g"
            and element.attrib.get("data-track-label") in labels
        }

    first_ids = conservation_ids(first)
    assert first_ids == conservation_ids(second)
    assert set(first_ids) == set(labels)
    assert len(set(first_ids.values())) == len(labels)
    assert "conservation_identity_legend" not in first_ids.values()
    assert all(identifier.startswith("track_slot_") for identifier in first_ids.values())


def test_circular_multi_merge_preserves_and_namespaces_defs_references() -> None:
    records = [
        SeqRecord(Seq("A" * 1000), id="record-a", name="record-a"),
        SeqRecord(Seq("A" * 1000), id="record-b", name="record-b"),
    ]
    hatch = HatchStyle(cross=True)
    style = RegionAnnotationStyle(hatch=hatch)
    annotations = tuple(
        RegionAnnotation(
            f"record-{record_number}-lane-{lane}",
            CoordinateSpan(
                parse_record_selector(f"#{record_number}"),
                lane * 100,
                lane * 100 + 80,
            ),
            lane=lane - 1,
            mark="band",
            style=style,
        )
        for record_number in (1, 2)
        for lane in (1, 2)
    )
    svg = assemble_circular_diagram_from_records(
        records,
        annotation_options=AnnotationOptions(
            sets=(AnnotationSet("regions", annotations),)
        ),
        circular_track_slots=[
            "regions:annotations@set_id=regions,w=8px,overflow=clip",
            "features:features",
            "ticks:ticks",
        ],
        legend="none",
        cfg=apply_config_overrides(None, {
            "labels.circular.scope": "none",
            "canvas.show_gc": False,
            "canvas.show_skew": False,
        }),
    ).tostring()

    _assert_unique_ids_and_resolved_local_references(svg)
    root = ET.fromstring(svg)
    clip_ids = {
        element.attrib["id"]
        for element in root.iter()
        if element.tag.rsplit("}", 1)[-1] == "clipPath"
    }
    pattern_ids = {
        element.attrib["id"]
        for element in root.iter()
        if element.tag.rsplit("}", 1)[-1] == "pattern"
    }
    assert len(clip_ids) == 2
    assert len(pattern_ids) == 2
    assert any(identifier.endswith("__record_2") for identifier in pattern_ids)

    _assert_unique_ids_and_resolved_local_references(enrich_svg(svg))
