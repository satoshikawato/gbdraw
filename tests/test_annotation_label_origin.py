"""A region label stands on the centre of the mark as it is drawn, also across the origin (OV-2577, OV-2578).

One rule decides whether a span crosses the origin (``segments_cross_origin``): the midpoint walks the span's
stretch across the origin (``region_run``). With display parts the label stands halfway along the region where the
display leaves its parts whole, else on the longest part ``mark_parts`` draws.
"""

from __future__ import annotations

import math
from xml.etree import ElementTree

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

import gbdraw
from gbdraw.annotations import AnnotationGradient, MarkPart, ResolvedRegionAnnotation, mark_parts, resolve_annotations
from gbdraw.api import (
    AnnotationOptions,
    AnnotationSet,
    CoordinateSpan,
    FeatureSpan,
    RecordDisplayOptions,
    RegionAnnotation,
    RegionAnnotationStyle,
)
from gbdraw.layout.record_coordinates import DisplayFragment, RecordDisplayTransform, SourceInterval
from gbdraw.render.drawers.annotation_marks import gradient_steps

_SVG = "{http://www.w3.org/2000/svg}"
_LENGTH = 10_000


def _record(topology: str = "circular", *, region_applied: bool = False) -> SeqRecord:
    record = SeqRecord(
        Seq("ACGT" * (_LENGTH // 4)), id="r1", name="r1", annotations={"topology": topology, "molecule_type": "DNA"}
    )
    record.features = [
        SeqFeature(
            CompoundLocation([FeatureLocation(9000, _LENGTH, strand=1), FeatureLocation(0, 1000, strand=1)]),
            type="CDS",
            qualifiers={"locus_tag": ["ENDS"]},
        ),
        SeqFeature(FeatureLocation(5000, 5500, strand=1), type="CDS", qualifiers={"locus_tag": ["MIDDLE"]}),
    ]
    if region_applied:
        record.annotations.update(gbdraw_region_applied=True, gbdraw_source_length=_LENGTH)
    return record


def _resolved(record: SeqRecord, *selectors: str) -> ResolvedRegionAnnotation:
    annotation = RegionAnnotation("a", FeatureSpan(None, selectors, envelope="segments"), label="A")
    (resolved,) = resolve_annotations(
        AnnotationOptions(sets=(AnnotationSet("s", (annotation,)),)), [record], mode="circular"
    ).annotations
    return resolved


@pytest.mark.parametrize(("topology", "region_applied"), (("linear", False), ("circular", True)))
def test_midpoint_of_segments_at_both_ends_of_a_record_without_an_origin_walks_them(
    topology: str, region_applied: bool
) -> None:
    resolved = _resolved(_record(topology, region_applied=region_applied), "locus_tag=ENDS")

    # The record has no origin to cross: halfway along 0-1000 then 9000-10000, not bp 0.
    assert not resolved.crosses_origin
    assert resolved.midpoint_bp == 1000.0


def test_three_segments_across_the_origin_keep_the_ascending_walk_and_gradient_extent() -> None:
    resolved = _resolved(_record(), "locus_tag=ENDS", "locus_tag=MIDDLE")
    gradient = AnnotationGradient("#0000ff", "#ff0000")
    middle = MarkPart(5000, 5500, True, True, 5000)

    # Which way such a span runs is open (OV-2656): its midpoint stays halfway along 0-1000, 5000-5500 and
    # 9000-10000 in ascending order, and its gradient scale spans only its stretch across the origin.
    assert resolved.crosses_origin
    assert resolved.midpoint_bp == 5250.0
    steps = gradient_steps(resolved, gradient, middle, 5000, 5500, record_length=_LENGTH, length_px=10.0)
    assert {color for _start, color in steps} == {"#ff0000"}


def _fragments(start: int | None, wrap: int | None = None, segments=((0, 1000), (9000, _LENGTH))) -> tuple[DisplayFragment, ...]:
    transform = RecordDisplayTransform(
        length=_LENGTH, source_base=1, source_step=1, start_coordinate=start, is_circular=True
    )
    if wrap is not None:
        transform = transform.with_wrap(wrap, _LENGTH)
    parts = transform.project_local_parts(tuple(SourceInterval(low, high) for low, high in segments))
    assert all(isinstance(part, DisplayFragment) for part in parts)
    return tuple(part for part in parts if isinstance(part, DisplayFragment))


def _origin_annotation(fragments: tuple[DisplayFragment, ...]) -> ResolvedRegionAnnotation:
    return ResolvedRegionAnnotation(
        id="a", set_id="s", record_index=0, segments=((0, 1000), (9000, _LENGTH)), midpoint_bp=0.0, span_bp=2000,
        label="A", mark="bracket", lane=None, style=RegionAnnotationStyle(), display_parts=fragments,
        crosses_origin=True,
    )


@pytest.mark.parametrize("join_origin", (True, False))
def test_label_of_a_span_across_the_origin_stands_on_its_joined_part(join_origin: bool) -> None:
    # Display start 5001: the two sides meet at bp 5000 of the display, on one row.
    annotation = _origin_annotation(_fragments(5001))

    parts = mark_parts(annotation, join_origin=join_origin)
    # Drawn as one stretch from 4000 to 6000, capped only at its two outer ends.
    assert [(part.start, part.end) for part in parts if part.capped_start] == [(4000, 6000)]
    assert [part.end for part in parts if part.capped_end] == [6000]
    assert annotation.label_midpoint_bp(join_origin=join_origin) == 5000.0
    assert annotation.label_room_bp(join_origin=join_origin) == 2000


def test_circle_joins_the_sides_where_the_display_cut_is_the_origin() -> None:
    annotation = _origin_annotation(_fragments(1))

    # On the circle the sides meet across the cut: the label stands at the origin, 12 o'clock.
    assert annotation.label_midpoint_bp(join_origin=True) % _LENGTH == 0.0
    # On a Linear row they lie at its two ends, each open where the span goes on over the origin.
    assert {(part.start, part.capped_start, part.end, part.capped_end) for part in mark_parts(annotation, join_origin=False)} == {
        (0, False, 1000, True), (9000, True, _LENGTH, False)
    }


@pytest.mark.parametrize(("start", "wrap"), ((None, 4000), (5001, 5000)))
def test_wrapped_rows_keep_the_sides_of_the_origin_apart(start: int | None, wrap: int) -> None:
    annotation = _origin_annotation(_fragments(start, wrap))

    # The sides lie on two rows, each open where the span goes on, and the label on one of them.
    assert all(
        part.capped_start != part.capped_end for part in mark_parts(annotation, join_origin=False)
    )
    assert annotation.label_room_bp(join_origin=False) == 1000


# Three segments: 1000-3000, 4000-5000 and 6000-6500; halfway along them (1750 bp in) is bp 2750.
_EXONS = ((1000, 3000), (4000, 5000), (6000, 6500))


def _exons(fragments: tuple[DisplayFragment, ...]) -> ResolvedRegionAnnotation:
    return ResolvedRegionAnnotation(
        id="a", set_id="s", record_index=0, segments=_EXONS, midpoint_bp=2750.0, span_bp=3500, label="A",
        mark="band", lane=None, style=RegionAnnotationStyle(), display_parts=fragments,
    )


@pytest.mark.parametrize(("start", "wrap", "expected"), ((9001, None, 3750.0), (None, 4000, 2750.0)))
def test_label_of_a_region_the_display_does_not_cut_stands_at_its_midpoint(
    start: int | None, wrap: int | None, expected: float
) -> None:
    annotation = _exons(_fragments(start, wrap, _EXONS))

    # The display start (bp 9000 at the left end) or the rows of 4000 bp leave every segment whole: the label
    # stands halfway along the region as without them, not on the centre of the longest segment (OV-2653).
    assert annotation.label_midpoint_bp(join_origin=False) == expected
    assert annotation.label_row(join_origin=False) == 0


def test_label_of_a_region_the_display_cut_cuts_stands_on_its_longest_part() -> None:
    # Display start 9501 cuts the side 9000-10000 at the display's ends: the label keeps the centre of the
    # longest part (OV-2654 and the cut case of OV-2653 are open).
    annotation = _origin_annotation(_fragments(9501))

    assert annotation.label_midpoint_bp(join_origin=True) == 750.0


def _label_item(
    mode: str, start: int | None, position: str = "center", *, feature: bool = False
) -> ElementTree.Element:
    target = (
        FeatureSpan(None, ("locus_tag=ENDS",))
        if feature
        else CoordinateSpan(None, 9001, 1000, wraps_origin=True)
    )
    annotation = RegionAnnotation(
        "origin",
        target,
        label="Origin",
        mark="band",
        style=RegionAnnotationStyle(fill="#ffcc00", label_orientation="horizontal", label_position=position),
    )
    options_type = gbdraw.CircularOptions if mode == "circular" else gbdraw.LinearOptions
    draw = gbdraw.draw_circular if mode == "circular" else gbdraw.draw_linear
    diagram = draw(
        _record(),
        options=options_type(annotations=AnnotationOptions(sets=(AnnotationSet("regions", (annotation,)),))),
        record_displays=[RecordDisplayOptions(start_coordinate=start)],
    )
    root = ElementTree.fromstring(diagram._drawing.tostring())
    (item,) = [group for group in root.iter(f"{_SVG}g") if group.get("data-gbdraw-annotation-id") == "origin"]
    return item


def _text_x(item: ElementTree.Element) -> float:
    (text,) = item.iter(f"{_SVG}text")
    return float(text.get("x", "nan"))


def test_circular_label_with_a_display_start_stands_on_the_centre_of_the_arc() -> None:
    (text,) = _label_item("circular", 5001).iter(f"{_SVG}text")
    x, y = float(text.get("x", "nan")), float(text.get("y", "nan"))

    # The arc runs from 0.4 to 0.6 turn (bp 9000 to 1000 with the origin at 6 o'clock): its label at 0.5 turn.
    assert math.degrees(math.atan2(x, -y)) % 360.0 == pytest.approx(180.0, abs=0.5)


def test_linear_label_with_a_display_start_stands_on_the_centre_of_the_mark() -> None:
    # The band runs from bp 9000 to 1000, drawn from x 800 to 1200 of 2000 with the origin at x 1000.
    assert _text_x(_label_item("linear", 5001)) == pytest.approx(1000.0)


# The Linear x of a start or an end label, for each kind of span (OV-2651, OV-2655). Without a display start the
# record's ends cut the region: a feature span's segments run 9000-10000 then 0-1000, so its start stands at x 1800
# and its end at x 200; a coordinate span's segments are ascending, so at x 0 and 2000. From display start 5001 the
# mark is one stretch from x 800 to 1200.
@pytest.mark.parametrize(
    ("start", "feature", "position", "expected"),
    (
        (None, True, "start", 1800.0),
        (None, True, "end", 200.0),
        (None, False, "start", 0.0),
        (None, False, "end", 2000.0),
        (5001, True, "start", 800.0),
        (5001, True, "end", 1200.0),
        (5001, False, "start", 800.0),
        (5001, False, "end", 1200.0),
    ),
)
def test_linear_label_at_the_start_or_end_stands_at_that_end_of_the_mark(
    start: int | None, feature: bool, position: str, expected: float
) -> None:
    assert _text_x(_label_item("linear", start, position, feature=feature)) == pytest.approx(expected)
