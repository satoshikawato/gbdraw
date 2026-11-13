"""A span across the origin is one region: no edge, cap or arrow head where it crosses (OV-2573).

Circular draws it as one part; Linear draws the two parts at the record's ends and leaves those ends open, or, from
a display start where the parts meet, one part (OV-2577).
"""

from __future__ import annotations

import math
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
from gbdraw.annotations import ResolvedRegionAnnotation

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


def _stretches(spans: list[tuple[float, float]]) -> list[tuple[float, float]]:
    """The drawn x ranges merged where they touch."""

    merged: list[tuple[float, float]] = []
    for start, end in sorted(spans):
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def _expected_ends(spans: list[tuple[float, float]]) -> set[float]:
    """The x of the span's two own ends.

    Unrotated, the span is drawn as two stretches at the record's ends (from 0, and up to the end of the axis):
    the record's ends are cuts, and the span's own ends are the inner ends. Rotated, its parts meet at the origin
    and draw one stretch, whose two ends are the span's.
    """

    stretches = _stretches(spans)
    if len(stretches) == 1:
        return set(stretches[0])
    (_first_start, first_end), (second_start, _second_end) = stretches
    return {first_end, second_start}


def _vertical_edges(item: ElementTree.Element) -> set[float]:
    """The x of every stroked edge across the lane: the vertical strokes of outline paths, and both sides of a
    rect drawn with its own stroke."""

    edges = {
        _number(x1)
        for path in item.iter(f"{_SVG}path")
        for x1, x2 in re.findall(r"M ([-\d.]+),[-\d.]+ L ([-\d.]+),[-\d.]+", path.get("d", ""))
        if _number(x1) == _number(x2)
    }
    for rect in item.iter(f"{_SVG}rect"):
        if rect.get("stroke") not in {None, "none"}:
            x = float(rect.get("x", 0))
            edges.update({_number(x), _number(x + float(rect.get("width", 0)))})
    return edges


@pytest.mark.parametrize("start", _STARTS)
@pytest.mark.parametrize(("mark", "line_cap"), _CASES)
def test_linear_leaves_the_ends_of_a_span_across_the_origin_open(mark: str, line_cap: str, start: int | None) -> None:
    item = _item("linear", mark, line_cap, start)
    # The span covers 2000 bp of 10000 on an axis of 2000 px: two stretches of 200 px at the record's ends, or,
    # rotated, one stretch of 400 px.
    expected_stretches = [(0.0, 200.0), (1800.0, 2000.0)] if start is None else [(800.0, 1200.0)]

    if mark in {"band", "highlight"}:
        spans = [
            (_number(rect.get("x")), _number(float(rect.get("x", 0)) + float(rect.get("width", 0))))
            for rect in item.iter(f"{_SVG}rect")
        ]
        assert _stretches(spans) == expected_stretches
        # Edges across the lane stand only at the span's own two ends, none at the origin or a record end.
        assert _vertical_edges(item) == _expected_ends(spans)
        return

    lines = list(item.iter(f"{_SVG}line"))
    bars = [line for line in lines if line.get("y1") == line.get("y2")]
    caps = [line for line in lines if line.get("x1") == line.get("x2")]
    spans = [(_number(bar.get("x1")), _number(bar.get("x2"))) for bar in bars]
    assert _stretches(spans) == expected_stretches
    expected = _expected_ends(spans)
    heads = [
        _number(match.group(1))
        for path in item.iter(f"{_SVG}path")
        if (match := re.search(r"L ([-\d.]+),[-\d.]+ L", path.get("d", "")))
    ]
    if line_cap == "arrow":
        # The feature runs forward: a cap where it starts, a head where it ends, nothing at the origin.
        cap_at, head_at = (
            sorted(expected, reverse=True) if start is None else sorted(expected)
        )
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


def _bracket_item(mode: str, record: SeqRecord, target: CoordinateSpan | FeatureSpan, start: int | None = None):
    annotation = RegionAnnotation(
        "origin", target, label="Origin", mark="bracket", style=RegionAnnotationStyle(stroke="#225588")
    )
    options_type = gbdraw.CircularOptions if mode == "circular" else gbdraw.LinearOptions
    draw = gbdraw.draw_circular if mode == "circular" else gbdraw.draw_linear
    diagram = draw(
        record,
        options=options_type(annotations=AnnotationOptions(sets=(AnnotationSet("r", (annotation,)),))),
        record_displays=[RecordDisplayOptions(start_coordinate=start)],
    )
    root = ElementTree.fromstring(diagram._drawing.tostring())
    (item,) = [group for group in root.iter(f"{_SVG}g") if group.get("data-gbdraw-annotation-id") == "origin"]
    return item


def _circular_cap_degrees(item: ElementTree.Element) -> list[float]:
    """The clockwise angle from 12 o'clock of each straight cap (a path with no arc), in degrees."""

    angles = []
    for path in item.iter(f"{_SVG}path"):
        match = re.fullmatch(r"M ([-\d.]+),([-\d.]+) L ([-\d.]+),([-\d.]+)", path.get("d", "").strip())
        if match:
            x1, y1, x2, y2 = map(float, match.groups())
            angles.append(round(math.degrees(math.atan2((x1 + x2) / 2, -(y1 + y2) / 2)) % 360.0, 1))
    return sorted(angles)


@pytest.mark.parametrize(
    ("start", "target", "ends"),
    (
        (501, CoordinateSpan(None, 9001, 1000, wraps_origin=True), (9000, 1000)),
        (9501, CoordinateSpan(None, 9001, 1000, wraps_origin=True), (9000, 1000)),
        (5001, CoordinateSpan(None, 4001, 6000), (4000, 6000)),
    ),
)
def test_circular_display_cut_inside_a_span_joins_its_parts(
    start: int, target: CoordinateSpan, ends: tuple[int, int]
) -> None:
    item = _bracket_item("circular", _record(), target, start)
    arcs = [path for path in item.iter(f"{_SVG}path") if " A " in f" {path.get('d', '')} "]

    # The circle goes on at the display cut (12 o'clock) as at the origin: one arc (OV-2654).
    assert len(arcs) == 1
    # Caps stand only at the span's own ends, not at the display cut or the origin.
    offset = start - 1
    expected = sorted(round(360.0 * ((bp - offset) % _LENGTH) / _LENGTH, 1) for bp in ends)
    assert _circular_cap_degrees(item) == expected


@pytest.mark.parametrize("mark", ("band", "highlight"))
def test_circular_fill_across_the_display_cut_is_one_closed_shape(mark: str) -> None:
    annotation = RegionAnnotation(
        "origin", CoordinateSpan(None, 4001, 6000), label="Origin", mark=mark,
        style=RegionAnnotationStyle(stroke="#225588", fill="#ffcc00"),
    )
    diagram = gbdraw.draw_circular(
        _record(),
        options=gbdraw.CircularOptions(annotations=AnnotationOptions(sets=(AnnotationSet("r", (annotation,)),))),
        record_displays=[RecordDisplayOptions(start_coordinate=5001)],
    )
    root = ElementTree.fromstring(diagram._drawing.tostring())
    (item,) = [group for group in root.iter(f"{_SVG}g") if group.get("data-gbdraw-annotation-id") == "origin"]
    paths = [path for path in item.iter(f"{_SVG}path") if path.get("stroke") not in {None, "none"}]

    # One closed shape with one stroke, no open edges at 12 o'clock (OV-2654).
    assert len(paths) == 1
    assert paths[0].get("d", "").rstrip().endswith("Z")


def _two_end_record(topology: str, *, region_applied: bool) -> SeqRecord:
    record = _record()
    record.annotations["topology"] = topology
    if region_applied:
        record.annotations.update(gbdraw_region_applied=True, gbdraw_source_length=_LENGTH)
    return record


@pytest.mark.parametrize(("topology", "region_applied"), (("linear", False), ("circular", True)))
def test_segments_at_both_ends_of_a_record_without_an_origin_keep_their_caps(topology: str, region_applied: bool) -> None:
    target = FeatureSpan(None, ("locus_tag=ORIGIN_CDS",), envelope="segments")
    item = _bracket_item("circular", _two_end_record(topology, region_applied=region_applied), target)

    # Two separate arcs, each capped at both its ends: no record origin joins them.
    assert len([path for path in item.iter(f"{_SVG}path") if " A " in f" {path.get('d', '')} "]) == 2
    assert len(_circular_cap_degrees(item)) == 4
    linear = _bracket_item("linear", _two_end_record(topology, region_applied=region_applied), target)
    assert len([line for line in linear.iter(f"{_SVG}line") if line.get("x1") == line.get("x2")]) == 4


@pytest.mark.parametrize("mark", ("band", "highlight"))
def test_circular_fill_over_a_whole_record_has_no_radial_edge_halfway_round(mark: str) -> None:
    annotation = RegionAnnotation(
        "origin",
        CoordinateSpan(None, 1, _LENGTH),
        label="Origin",
        mark=mark,
        style=RegionAnnotationStyle(stroke="#225588", fill="#ffcc00"),
    )
    diagram = gbdraw.draw_circular(
        _record(), options=gbdraw.CircularOptions(annotations=AnnotationOptions(sets=(AnnotationSet("r", (annotation,)),)))
    )
    root = ElementTree.fromstring(diagram._drawing.tostring())
    (item,) = [group for group in root.iter(f"{_SVG}g") if group.get("data-gbdraw-annotation-id") == "origin"]
    edges = [
        edge
        for path in item.iter(f"{_SVG}path")
        if path.get("stroke") not in {None, "none"}
        for edge in re.findall(r"M [-\d.]+,[-\d.]+ L [-\d.]+,[-\d.]+", path.get("d", ""))
    ]

    # The drawer splits a whole turn in two pieces; every radial edge stands at the top, where the turn starts and
    # ends, none at the split halfway round (OV-2574).
    assert edges
    assert all(edge.startswith("M 0,-") and " L 0,-" in edge for edge in edges)


def test_resolved_annotation_crossing_the_origin_needs_a_segment_at_zero() -> None:
    def build(segments: tuple[tuple[int, int], ...]) -> ResolvedRegionAnnotation:
        return ResolvedRegionAnnotation(
            id="a", set_id="s", record_index=0, segments=segments, midpoint_bp=0.0, span_bp=1, label="", mark="band",
            lane=None, style=RegionAnnotationStyle(), crosses_origin=True,
        )

    assert build(((9000, _LENGTH), (0, 1000))).crosses_origin
    for segments in (((0, 1000),), ((9000, _LENGTH), (10, 1000))):
        with pytest.raises(ValueError, match="crosses the origin"):
            build(segments)
