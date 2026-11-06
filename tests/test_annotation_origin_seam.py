"""A span across the origin is one region: no edge, cap or arrow head where it crosses (OV-2573).

Circular draws it as one part; Linear draws the two parts at the record's ends and leaves those ends open.
"""

from __future__ import annotations

import re
from xml.etree import ElementTree

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

import gbdraw
from gbdraw.api import (
    AnnotationOptions,
    AnnotationSet,
    CoordinateSpan,
    FeatureSpan,
    RecordDisplayOptions,
    RegionAnnotation,
    RegionAnnotationStyle,
)

_SVG = "{http://www.w3.org/2000/svg}"
_LENGTH = 10_000
# 1-based starts of the display: unset, and one that puts the origin in the middle of the drawing.
_STARTS = (None, 5001)


def _record() -> SeqRecord:
    record = SeqRecord(
        Seq("ACGT" * (_LENGTH // 4)), id="r1", name="r1", annotations={"topology": "circular", "molecule_type": "DNA"}
    )
    # A directional feature across the origin, the target of the arrow case.
    record.features = [
        SeqFeature(
            CompoundLocation([FeatureLocation(9000, _LENGTH, strand=1), FeatureLocation(0, 1000, strand=1)]),
            type="CDS",
            qualifiers={"locus_tag": ["ORIGIN_CDS"]},
        )
    ]
    return record


def _annotation(mark: str, line_cap: str) -> RegionAnnotation:
    target = (
        FeatureSpan(None, ("locus_tag=ORIGIN_CDS",))
        if line_cap == "arrow"
        else CoordinateSpan(None, 9001, 1000, wraps_origin=True)
    )
    fill = {"fill": "#ffcc00"} if mark in {"band", "highlight"} else {}
    return RegionAnnotation(
        "origin",
        target,
        label="Origin",
        mark=mark,
        style=RegionAnnotationStyle(stroke="#225588", line_cap=line_cap, **fill),
    )


def _item(mode: str, mark: str, line_cap: str, start: int | None) -> ElementTree.Element:
    options_type = gbdraw.CircularOptions if mode == "circular" else gbdraw.LinearOptions
    draw = gbdraw.draw_circular if mode == "circular" else gbdraw.draw_linear
    diagram = draw(
        _record(),
        options=options_type(
            annotations=AnnotationOptions(sets=(AnnotationSet("regions", (_annotation(mark, line_cap),)),))
        ),
        record_displays=[RecordDisplayOptions(start_coordinate=start)],
    )
    root = ElementTree.fromstring(diagram._drawing.tostring())
    (item,) = [
        group for group in root.iter(f"{_SVG}g") if group.get("data-gbdraw-annotation-id") == "origin"
    ]
    return item


_CASES = (
    ("band", "tick"),
    ("highlight", "tick"),
    ("bracket", "tick"),
    ("bracket", "arrow"),
)


@pytest.mark.parametrize("start", _STARTS)
@pytest.mark.parametrize(("mark", "line_cap"), _CASES)
def test_circular_draws_a_span_across_the_origin_as_one_part(mark: str, line_cap: str, start: int | None) -> None:
    item = _item("circular", mark, line_cap, start)
    paths = [path for path in item.iter(f"{_SVG}path") if path.get("stroke") not in {None, "none"}]

    if mark in {"band", "highlight"}:
        # One closed shape with one stroke.
        assert len(paths) == 1
        assert paths[0].get("d", "").rstrip().endswith("Z")
    else:
        arcs = [path for path in paths if " A " in f" {path.get('d', '')} "]
        caps = [path for path in paths if path not in arcs]
        # One arc and a cap at each of its two ends, none at the origin.
        assert (len(arcs), len(caps)) == (1, 2)


def _number(text: str | None) -> float:
    return round(float(text or 0.0), 3)


def _expected_ends(spans: list[tuple[float, float]]) -> set[float]:
    """The x of the span's two outer ends, from its two parts in drawing order.

    Unrotated, the parts lie at the record's ends (from 0, and up to the end of the axis), so the record's ends
    are cuts; rotated, the parts meet at the origin, so the edge they share is the cut.
    """

    (first_start, first_end), (second_start, second_end) = spans
    if first_end == second_start:
        return {first_start, second_end}
    return {first_end, second_start}


@pytest.mark.parametrize("start", _STARTS)
@pytest.mark.parametrize(("mark", "line_cap"), _CASES)
def test_linear_leaves_the_ends_of_a_span_across_the_origin_open(mark: str, line_cap: str, start: int | None) -> None:
    item = _item("linear", mark, line_cap, start)

    if mark in {"band", "highlight"}:
        rects = list(item.iter(f"{_SVG}rect"))
        assert len(rects) == 2
        spans = sorted(
            (_number(rect.get("x")), _number(float(rect.get("x", 0)) + float(rect.get("width", 0)))) for rect in rects
        )
        # The rects carry no stroke of their own; their outlines draw only the span's two outer ends across the lane.
        assert {rect.get("stroke") for rect in rects} == {"none"}
        verticals = {
            _number(x1)
            for path in item.iter(f"{_SVG}path")
            for x1, x2 in re.findall(r"M ([-\d.]+),[-\d.]+ L ([-\d.]+),[-\d.]+", path.get("d", ""))
            if _number(x1) == _number(x2)
        }
        assert verticals == _expected_ends(spans)
        return

    lines = list(item.iter(f"{_SVG}line"))
    bars = [line for line in lines if line.get("y1") == line.get("y2")]
    caps = [line for line in lines if line.get("x1") == line.get("x2")]
    assert len(bars) == 2
    spans = sorted((_number(bar.get("x1")), _number(bar.get("x2"))) for bar in bars)
    expected = _expected_ends(spans)
    heads = [
        _number(match.group(1))
        for path in item.iter(f"{_SVG}path")
        if (match := re.search(r"L ([-\d.]+),[-\d.]+ L", path.get("d", "")))
    ]
    if line_cap == "arrow":
        # The feature runs forward: a cap where its first part starts, a head where its last part ends, nothing at
        # the origin.
        starts = {start for start, _end in spans}
        cap_at, head_at = sorted(expected, key=lambda x: x not in starts)
        assert ({_number(cap.get("x1")) for cap in caps}, set(heads)) == ({cap_at}, {head_at})
    else:
        assert sorted(_number(cap.get("x1")) for cap in caps) == sorted(expected)
        assert heads == []


def test_circular_bracket_over_a_whole_record_is_capped_only_at_its_ends() -> None:
    annotation = RegionAnnotation(
        "origin", CoordinateSpan(None, 1, _LENGTH), label="Origin", mark="bracket", style=RegionAnnotationStyle()
    )
    diagram = gbdraw.draw_circular(
        _record(), options=gbdraw.CircularOptions(annotations=AnnotationOptions(sets=(AnnotationSet("r", (annotation,)),)))
    )
    root = ElementTree.fromstring(diagram._drawing.tostring())
    (item,) = [group for group in root.iter(f"{_SVG}g") if group.get("data-gbdraw-annotation-id") == "origin"]
    caps = [path.get("d", "") for path in item.iter(f"{_SVG}path") if " A " not in f" {path.get('d', '')} "]

    # The drawer splits a whole turn in two pieces; the split is no end (OV-2574): both caps stand at the top.
    assert len(caps) == 2
    assert all(re.match(r"M 0,-[\d.]+ L 0,-[\d.]+$", cap) for cap in caps)
