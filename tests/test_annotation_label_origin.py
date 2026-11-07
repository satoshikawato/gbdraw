"""A region label stands on the centre of the mark as it is drawn, also across the origin (OV-2577, OV-2578).

One rule decides whether a span crosses the origin (``segments_cross_origin``): the midpoint walks the span in the
order it runs (``region_run``), and with a display start the label stands on the longest part ``mark_parts`` draws.
"""

from __future__ import annotations

import math
from xml.etree import ElementTree

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

import gbdraw
from gbdraw.annotations import ResolvedRegionAnnotation, mark_parts, resolve_annotations
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


def test_midpoint_of_three_segments_across_the_origin_walks_the_joined_span() -> None:
    resolved = _resolved(_record(), "locus_tag=ENDS", "locus_tag=MIDDLE")

    # The span runs from bp 9000 over the origin to bp 1000, then 5000-5500: halfway (1250 bp in) is bp 250.
    assert resolved.crosses_origin
    assert resolved.midpoint_bp == 250.0


def _fragments(start: int | None, wrap: int | None = None) -> tuple[DisplayFragment, ...]:
    transform = RecordDisplayTransform(
        length=_LENGTH, source_base=1, source_step=1, start_coordinate=start, is_circular=True
    )
    if wrap is not None:
        transform = transform.with_wrap(wrap, _LENGTH)
    parts = transform.project_local_parts((SourceInterval(0, 1000), SourceInterval(9000, _LENGTH)))
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

    (part,) = mark_parts(annotation, join_origin=join_origin)
    assert (part.start, part.end, part.capped_start, part.capped_end) == (4000, 6000, True, True)
    assert annotation.label_midpoint_bp(join_origin=join_origin) == 5000.0
    assert annotation.label_room_bp(join_origin=join_origin) == 2000


def test_circle_joins_the_sides_where_the_display_cut_is_the_origin() -> None:
    annotation = _origin_annotation(_fragments(1))

    # On the circle the sides meet across the cut (one arc, 9000 to 11000); on a Linear row they lie at its ends.
    assert annotation.label_midpoint_bp(join_origin=True) == 10_000.0
    assert len(mark_parts(annotation, join_origin=False)) == 2


@pytest.mark.parametrize(("start", "wrap"), ((None, 4000), (5001, 5000)))
def test_wrapped_rows_keep_the_sides_of_the_origin_apart(start: int | None, wrap: int) -> None:
    annotation = _origin_annotation(_fragments(start, wrap))

    parts = mark_parts(annotation, join_origin=False)
    # The sides lie on two rows: two parts, open where the span goes on, and the label on one of them.
    assert len(parts) == 2
    assert len({part.row for part in parts}) == 2
    assert annotation.label_room_bp(join_origin=False) == 1000


def _label_item(mode: str, start: int | None, position: str = "center") -> ElementTree.Element:
    annotation = RegionAnnotation(
        "origin",
        CoordinateSpan(None, 9001, 1000, wraps_origin=True),
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


def test_circular_label_with_a_display_start_stands_on_the_centre_of_the_arc() -> None:
    (text,) = _label_item("circular", 5001).iter(f"{_SVG}text")
    x, y = float(text.get("x", "nan")), float(text.get("y", "nan"))

    # The arc runs from 0.4 to 0.6 turn (bp 9000 to 1000 with the origin at 6 o'clock): its label at 0.5 turn.
    assert math.degrees(math.atan2(x, -y)) % 360.0 == pytest.approx(180.0, abs=0.5)


def test_linear_label_with_a_display_start_stands_on_the_centre_of_the_mark() -> None:
    item = _label_item("linear", 5001)
    (text,) = item.iter(f"{_SVG}text")
    rects = list(item.iter(f"{_SVG}rect"))

    # One band where the two sides meet, and its label over its centre.
    assert len(rects) == 1
    centre = float(rects[0].get("x", "nan")) + 0.5 * float(rects[0].get("width", "nan"))
    assert float(text.get("x", "nan")) == pytest.approx(centre)


@pytest.mark.parametrize("position", ("start", "end"))
def test_linear_label_at_the_start_or_end_stands_at_that_end_of_the_mark(position: str) -> None:
    item = _label_item("linear", 5001, position)
    (text,) = item.iter(f"{_SVG}text")
    (rect,) = item.iter(f"{_SVG}rect")
    left = float(rect.get("x", "nan"))
    right = left + float(rect.get("width", "nan"))

    # The display start puts the region's segment at 0 first in its parts; the label still goes to the mark's end
    # on its side (OV-2651).
    assert float(text.get("x", "nan")) == pytest.approx(left if position == "start" else right)
