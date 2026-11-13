"""A region label stands on the centre of the mark as it is drawn, also across the origin (OV-2577, OV-2578).

One rule decides whether a span crosses the origin (``segments_cross_origin``): the midpoint walks the span's
segments in ``region_run`` order, over the origin from one side into the other. With display parts the label stands halfway along the region where the
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
from gbdraw.layout.record_coordinates import CirclePlace, DisplayFragment, RecordDisplayTransform, SourceInterval
from gbdraw.render.drawers.annotation_marks import gradient_steps

_SVG = "{http://www.w3.org/2000/svg}"
_LENGTH = 10_000


def _half_stroke_reach(_annotation, style, _lane_width: float) -> float:
    # The packer tests' own mark paint (``LabelSizing.reach``): half the stroke past the lane's edge.
    return 0.5 * float(style.stroke_width)



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


def test_three_segments_across_the_origin_run_from_the_stretch_over_the_origin() -> None:
    resolved = _resolved(_record(), "locus_tag=ENDS", "locus_tag=MIDDLE")
    gradient = AnnotationGradient("#0000ff", "#ff0000")
    middle = MarkPart(5000, 5500, True, True, 5000)

    # PD-OI-119: the region runs 9000-10000, over the origin 0-1000, then 5000-5500. Its
    # midpoint is 1250 bp along that run, bp 250, and its gradient spans the run from bp 9000 to bp 5500 (one record
    # length on), so 5000-5500 takes the last part of the scale.
    assert resolved.crosses_origin
    assert resolved.midpoint_bp == 250.0
    steps = gradient_steps(resolved, gradient, middle, 5000, 5500, record_length=_LENGTH, length_px=10.0)
    assert [steps[0][1], steps[-1][1]] == ["#ed0011", "#fd0001"]


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
    assert annotation.label_stretch("center", join_origin=join_origin)[1] == (4000, 6000)


def test_circle_joins_the_sides_where_the_display_cut_is_the_origin() -> None:
    annotation = _origin_annotation(_fragments(1))

    # On the circle the sides meet across the cut: the label stands at the origin, 12 o'clock.
    assert annotation.label_midpoint_bp(join_origin=True) % _LENGTH == 0.0
    # On a Linear row they lie at its two ends, each open where the span goes on over the origin.
    assert {(part.start, part.capped_start, part.end, part.capped_end) for part in mark_parts(annotation, join_origin=False)} == {
        (0, False, 1000, True), (9000, True, _LENGTH, False)
    }


@pytest.mark.parametrize(
    ("start", "wrap", "centre", "first", "last"),
    (
        # Rows of 4000 bp: 9000-10000 at the start of the last row (row 2), 0-1000 at the start of row 0. Only the
        # origin parts them: the label takes the part where the region starts.
        (None, 4000, (2, 1500.0), (2, 1000.0), (0, 1000.0)),
        # From display start 5001 in rows of 5000 bp: 9000-10000 ends row 0 and 0-1000 starts row 1. The row
        # boundary parts them: the label takes the earlier row.
        (5001, 5000, (0, 4500.0), (0, 4000.0), (1, 1000.0)),
    ),
)
def test_wrapped_rows_keep_the_sides_of_the_origin_apart(
    start: int | None, wrap: int, centre: tuple[int, float], first: tuple[int, float], last: tuple[int, float]
) -> None:
    annotation = _origin_annotation(_fragments(start, wrap))

    # The sides lie on two rows, each open where the span goes on, and each label on the row of its place.
    assert all(
        part.capped_start != part.capped_end for part in mark_parts(annotation, join_origin=False)
    )
    places = {position: annotation.label_anchor(position, join_origin=False) for position in ("center", "start", "end")}
    assert [(place.part.row, place.bp) for place in places.values()] == [centre, first, last]
    # The centred label may cover the 1000 bp part it stands on.
    low, high = annotation.label_stretch("center", join_origin=False)[1]
    assert high - low == 1000


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
    assert annotation.label_place(join_origin=False) == (expected, mark_parts(annotation, join_origin=False)[0])


def test_label_of_a_region_the_display_cut_cuts_stands_on_its_longest_part() -> None:
    # PD-OI-119: display start 9501 cuts the side 9000-10000 at the ends of a Linear bar (500 bp at each end); 0-1000
    # follows it from display bp 500 to 1500: the label keeps the centre of that longest part, bp 750 (x 150 of a
    # bar of 2000 px).
    annotation = _origin_annotation(_fragments(9501))

    assert annotation.label_midpoint_bp(join_origin=False) == 750.0


def test_a_circle_joins_a_region_at_its_display_cut() -> None:
    # The circle goes on at the display cut (display bp 0, bp 9500): the region is one part from display bp 9500
    # over the cut to 11500, capped at its own ends, and its label stands halfway along it, at bp 0 of the record,
    # display bp 500 (OV-2654).
    annotation = _origin_annotation(_fragments(9501))

    assert mark_parts(annotation, join_origin=True) == (MarkPart(9500, 11500, True, True, 9000),)
    assert annotation.label_midpoint_bp(join_origin=True) == 500.0


def _on_a_shared_circle(start: int | None, segments: tuple[tuple[int, int], ...]) -> ResolvedRegionAnnotation:
    # The record is drawn from display bp 3000 of a circle of 25000 display bp that it shares with other records.
    transform = RecordDisplayTransform(_LENGTH, 1, 1, start, True).with_circle(
        CirclePlace(offset=3000, length=25_000, scale_length=_LENGTH), _LENGTH
    )
    parts = transform.project_local_parts(tuple(SourceInterval(low, high) for low, high in segments))
    return ResolvedRegionAnnotation(
        id="a", set_id="s", record_index=0, segments=segments, midpoint_bp=0.0, span_bp=2000, label="A",
        mark="band", lane=None, style=RegionAnnotationStyle(),
        display_parts=tuple(part for part in parts if isinstance(part, DisplayFragment)),
        shares_circle=transform.circle is not None, crosses_origin=segments[0][0] == 0 and len(segments) > 1,
    )


def test_a_record_that_shares_one_circle_ends_a_region_at_its_display_cut() -> None:
    # From display start 5001 the record's arc runs from display bp 3000 to 13000 and its display cut (bp 5000) ends
    # it, where the gap to the next record starts: 4000-6000 keeps its two parts, open at the cut, and its label
    # the centre of its longest part, the earlier of two equal ones (PD-OI-119; OV-2654 review).
    annotation = _on_a_shared_circle(5001, ((4000, 6000),))

    assert mark_parts(annotation, join_origin=True) == (
        MarkPart(12_000, 13_000, True, False, 4000), MarkPart(3000, 4000, False, True, 5000)
    )
    assert annotation.label_place(join_origin=True) == (3500.0, MarkPart(3000, 4000, False, True, 5000))


@pytest.mark.parametrize("shared", (False, True))
def test_a_region_knows_whether_its_record_shares_its_circle(shared: bool) -> None:
    from types import SimpleNamespace

    from gbdraw.annotations.planning import prepare_annotation_track_slots

    transform = RecordDisplayTransform(_LENGTH, 1, 1, 5001, True)
    if shared:
        transform = transform.with_circle(CirclePlace(offset=3000, length=25_000, scale_length=_LENGTH), _LENGTH)
    _slots, bundle, _auto = prepare_annotation_track_slots(
        AnnotationOptions(sets=(AnnotationSet("s", (RegionAnnotation("x", CoordinateSpan(None, 4001, 6000)),)),)),
        [_record()], None, mode="circular", default_slots=list,
        slot_factory=lambda **fields: SimpleNamespace(width=None, enabled=True, **fields), record_transforms=[transform],
    )

    # The display parts and the fact that the record's arc ends at its display cut come from one transform.
    (annotation,) = bundle.annotations
    assert annotation.display_parts is not None
    assert annotation.shares_circle is shared


@pytest.mark.parametrize(
    ("start", "expected"),
    (
        # Without a display start the origin lies at the two ends of the record's arc: the sides stay apart there,
        # each open where the span goes on (OV-3519).
        (None, (MarkPart(3000, 4000, False, True, 0), MarkPart(12_000, 13_000, True, False, 9000))),
        # From display start 5001 the sides meet in the middle of the arc, at display bp 8000: one part.
        (5001, (MarkPart(7000, 9000, True, True, 9000),)),
    ),
)
def test_a_record_that_shares_one_circle_joins_the_origin_only_where_its_sides_meet(
    start: int | None, expected: tuple[MarkPart, ...]
) -> None:
    annotation = _on_a_shared_circle(start, ((0, 1000), (9000, _LENGTH)))

    assert mark_parts(annotation, join_origin=True) == expected


@pytest.mark.parametrize(
    ("start", "segments", "first", "last"),
    (
        # 4000-6000 from display start 5001: one part from display bp 9000 over the cut to 11000.
        (5001, ((4000, 6000),), 9000.0, 11_000.0),
        # 9000-10000 and 0-1000 from 501 or 9501: one part over the display cut and the origin.
        (501, ((0, 1000), (9000, _LENGTH)), 8500.0, 10_500.0),
        (9501, ((0, 1000), (9000, _LENGTH)), 9500.0, 11_500.0),
    ),
)
def test_start_and_end_stand_on_the_part_a_circle_joins_at_its_display_cut(
    start: int, segments: tuple[tuple[int, int], ...], first: float, last: float
) -> None:
    annotation = ResolvedRegionAnnotation(
        id="a", set_id="s", record_index=0, segments=segments, midpoint_bp=0.0, span_bp=2000, label="A",
        mark="band", lane=None, style=RegionAnnotationStyle(), display_parts=_fragments(start, segments=segments),
        crosses_origin=len(segments) > 1,
    )
    (part,) = mark_parts(annotation, join_origin=True)

    # The first and last bases stand on the joined part, where it draws them (review LOW 3: no StopIteration).
    assert annotation.label_anchor("start", join_origin=True) == (first, part)
    assert annotation.label_anchor("end", join_origin=True) == (last, part)


def test_a_row_boundary_ties_to_the_earlier_row() -> None:
    # PD-OI-119: rows of 4000 bp cut 3000-5000 into two equal parts, 3000-4000 on row 0 and 0-1000 on row 1: the
    # label takes the earlier row, at display bp 3500 (x 1750).
    segments = ((3000, 5000),)
    annotation = ResolvedRegionAnnotation(
        id="a", set_id="s", record_index=0, segments=segments, midpoint_bp=4000.0, span_bp=2000, label="A",
        mark="band", lane=None, style=RegionAnnotationStyle(), display_parts=_fragments(None, 4000, segments),
    )
    place = annotation.label_place(join_origin=False)

    assert (place.part.row, place.bp) == (0, 3500.0)


_TARGETS = {
    # A feature span whose envelope wraps the origin: segments 9000-10000 then 0-1000.
    "feature": FeatureSpan(None, ("locus_tag=ENDS",)),
    # The same feature's segments: 0-1000 and 9000-10000.
    "segments": FeatureSpan(None, ("locus_tag=ENDS",), envelope="segments"),
    "coordinate": CoordinateSpan(None, 9001, 1000, wraps_origin=True),
    # Unequal sides: 9500-10000 and 0-2000.
    "unequal": CoordinateSpan(None, 9501, 2000, wraps_origin=True),
    # Spans that the display cut of display start 5001 (bp 5000) crosses: 4000-6000 in equal sides, 4500-6000 not.
    "cut": CoordinateSpan(None, 4001, 6000),
    "cut_unequal": CoordinateSpan(None, 4501, 6000),
}


def _label_item(
    mode: str, start: int | None, position: str = "center", *, feature: bool = False, kind: str | None = None
) -> ElementTree.Element:
    target = _TARGETS[kind or ("feature" if feature else "coordinate")]
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


@pytest.mark.parametrize(("kind", "degrees"), (("cut", 0.0), ("cut_unequal", 9.0)))
def test_circular_label_stands_halfway_along_a_span_the_display_cut_crosses(kind: str, degrees: float) -> None:
    # From display start 5001 the display cut (12 o'clock, bp 5000) lies inside the span, and the circle goes on
    # there: the label stands halfway along the span (bp 5000 or 5250), not on the centre of one side (OV-2654).
    (text,) = _label_item("circular", 5001, kind=kind).iter(f"{_SVG}text")
    x, y = float(text.get("x", "nan")), float(text.get("y", "nan"))
    angle = math.degrees(math.atan2(x, -y)) % 360.0

    assert min(abs(angle - degrees), 360.0 - abs(angle - degrees)) == pytest.approx(0.0, abs=0.5)


def test_linear_label_with_a_display_start_stands_on_the_centre_of_the_mark() -> None:
    # The band runs from bp 9000 to 1000, drawn from x 800 to 1200 of 2000 with the origin at x 1000.
    assert _text_x(_label_item("linear", 5001)) == pytest.approx(1000.0)


# PD-OI-119: a start label stands at the region's first base and an end label at its last, for
# every kind of span, with and without a display start. A span across the origin starts at bp 9000 (x 1800 of
# 2000) and ends at bp 1000 (x 200); from display start 5001 the mark is one stretch from x 800 to 1200; from
# display start 9501 the region starts at x 1900 and ends at x 300.
@pytest.mark.parametrize("kind", ("feature", "segments", "coordinate"))
@pytest.mark.parametrize(
    ("start", "position", "expected"),
    (
        (None, "start", 1800.0),
        (None, "end", 200.0),
        (1, "start", 1800.0),
        (1, "end", 200.0),
        (5001, "start", 800.0),
        (5001, "end", 1200.0),
        (9501, "start", 1900.0),
        (9501, "end", 300.0),
    ),
)
def test_linear_label_at_the_start_or_end_stands_at_the_first_or_last_base(
    kind: str, start: int | None, position: str, expected: float
) -> None:
    assert _text_x(_label_item("linear", start, position, kind=kind)) == pytest.approx(expected)


@pytest.mark.parametrize("start", (None, 1))
def test_unequal_sides_start_and_end_at_the_first_and_last_base(start: int | None) -> None:
    # Region 9501-2000: its first base is bp 9500 (x 1900) and its last bp 2000 (x 400).
    assert _text_x(_label_item("linear", start, "start", kind="unequal")) == pytest.approx(1900.0)
    assert _text_x(_label_item("linear", start, "end", kind="unequal")) == pytest.approx(400.0)


def _anchor_record() -> SeqRecord:
    record = _record()
    record.features.append(
        SeqFeature(
            CompoundLocation(
                [FeatureLocation(1000, 3000, strand=1), FeatureLocation(4100, 5600, strand=1), FeatureLocation(5700, 7200, strand=1)]
            ),
            type="CDS",
            qualifiers={"locus_tag": ["EXONS"]},
        )
    )
    return record


def _linear_labels(
    selectors: tuple[str, ...],
    positions: tuple[str, ...],
    *,
    start: int | None = None,
    wrap: int | None = None,
    record: SeqRecord | None = None,
) -> dict[str, tuple[float, float]]:
    """The x and track-local y of the label of a band over ``selectors`` at each label position."""

    places = {}
    for position in positions:
        annotation = RegionAnnotation(
            "x",
            FeatureSpan(None, selectors, envelope="segments"),
            label="X",
            mark="band",
            style=RegionAnnotationStyle(fill="#ffcc00", label_orientation="horizontal", label_position=position),
        )
        diagram = gbdraw.draw_linear(
            record or _anchor_record(),
            options=gbdraw.LinearOptions(
                annotations=AnnotationOptions(sets=(AnnotationSet("regions", (annotation,)),)),
                config_overrides={"canvas.linear.wrap_length": wrap} if wrap is not None else {},
            ),
            record_displays=[RecordDisplayOptions(start_coordinate=start)],
        )
        root = ElementTree.fromstring(diagram._drawing.tostring())
        (item,) = [group for group in root.iter(f"{_SVG}g") if group.get("data-gbdraw-annotation-id") == "x"]
        (text,) = item.iter(f"{_SVG}text")
        places[position] = (float(text.get("x", "nan")), float(text.get("y", "nan")))
    return places


def test_wrapped_start_and_end_labels_stand_at_the_first_and_last_base() -> None:
    # Rows of 4000 bp: 1000-3000 on row 0, 4100-5600 and 5700-7200 on row 1, where the midpoint (bp 4600) puts
    # the centred label. The region's first base is on row 0 (x 500) and its last on row 1 (x 1600; PD-OI-119).
    places = _linear_labels(("locus_tag=EXONS",), ("start", "end", "center"), wrap=4000)

    assert places["start"][0] == pytest.approx(500.0)
    assert places["end"][0] == pytest.approx(1600.0)
    assert places["center"][0] == pytest.approx(300.0)
    assert places["start"][1] < places["end"][1] == places["center"][1]


def test_wrapped_labels_of_a_span_across_the_origin_stand_on_the_rows_of_their_places() -> None:
    # Rows of 4000 bp: the region starts at bp 9000 on row 2 (x 500) and ends at bp 1000 on row 0 (x 500); the
    # centred label takes 9000-10000, where the region starts (x 750 on row 2).
    places = _linear_labels(("locus_tag=ENDS",), ("start", "end", "center"), wrap=4000)

    assert [places[position][0] for position in ("start", "end", "center")] == pytest.approx([500.0, 500.0, 750.0])
    assert places["start"][1] == places["center"][1] > places["end"][1]


def test_wrapped_labels_with_a_display_start_stand_on_the_rows_of_their_places() -> None:
    # Display start 5001, rows of 5000 bp (2000 px): the region runs from display bp 4000 (row 0, x 1600) over the
    # row boundary to display bp 1000 of row 1 (x 400); the row boundary parts it into two equal parts, and the
    # centred label takes the earlier row (x 1800).
    places = _linear_labels(("locus_tag=ENDS",), ("start", "end", "center"), start=5001, wrap=5000)

    assert [places[position][0] for position in ("start", "end", "center")] == pytest.approx([1600.0, 400.0, 1800.0])
    assert places["start"][1] == places["center"][1] < places["end"][1]


def test_centred_label_of_three_parts_the_origin_cuts_stands_on_the_longest_drawn_part() -> None:
    # PD-OI-119, delegated D-A5: a CDS join(9901..10000,1..100,3001..8000) without a display start is drawn as three
    # parts, two of 100 bp at the bar's ends and 3000-8000. The centred label stands on the longest drawn part
    # (bp 5500, x 1100), not on the tie of the two parts at the origin (bp 9950) nor at the run's midpoint (bp 5400).
    record = _record()
    record.features.append(
        SeqFeature(
            CompoundLocation(
                [FeatureLocation(9900, _LENGTH, strand=1), FeatureLocation(0, 100, strand=1), FeatureLocation(3000, 8000, strand=1)]
            ),
            type="CDS",
            qualifiers={"locus_tag": ["THREE"]},
        )
    )

    for start in (None, 1):
        places = _linear_labels(("locus_tag=THREE",), ("center",), start=start, record=record)
        assert places["center"][0] == pytest.approx(1100.0)


@pytest.mark.parametrize(("start", "expected"), ((None, (1800.0, 1100.0, 1900.0)), (5001, (800.0, 100.0, 1050.0))))
def test_three_segments_across_the_origin_start_end_and_centre_follow_the_run(
    start: int | None, expected: tuple[float, float, float]
) -> None:
    # The region runs 9000-10000, 0-1000, 5000-5500 (PD-OI-119): it starts at bp 9000 and ends at bp 5500.
    # Without a display start the record's ends cut it, so the centred label stands on the centre of the longest
    # part; 0-1000 and 9000-10000 tie, and it takes 9000-10000, where the region starts (Q2). From display start 5001 every part is whole and the centred
    # label stands at its midpoint, bp 250 (x 1050).
    places = _linear_labels(("locus_tag=ENDS", "locus_tag=MIDDLE"), ("start", "end", "center"), start=start)

    assert tuple(places[position][0] for position in ("start", "end", "center")) == pytest.approx(expected)


@pytest.mark.parametrize(("kind", "expected"), (("feature", 1900.0), ("coordinate", 1900.0), ("unequal", 200.0)))
@pytest.mark.parametrize("start", (None, 1))
def test_centred_label_of_a_span_cut_at_the_origin_stands_on_its_longer_part(
    kind: str, start: int | None, expected: float
) -> None:
    # PD-OI-119: the bar's ends cut the span at the origin; the label stands on the centre of the
    # longer part, and of two equal parts on the one where the region starts, 9000-10000 at the bar's right end.
    assert _text_x(_label_item("linear", start, kind=kind)) == pytest.approx(expected)


def test_joined_arc_gradient_of_three_segments_across_the_origin_takes_the_first_part_of_the_run() -> None:
    resolved = _resolved(_record(), "locus_tag=ENDS", "locus_tag=MIDDLE")
    gradient = AnnotationGradient("#0000ff", "#ff0000")
    joined = MarkPart(9000, 11_000, True, True, 9000)

    # The scale runs from blue at bp 9000 to red at bp 5500 one record length on (PD-OI-119): the arc over the origin
    # takes its first part.
    steps = gradient_steps(resolved, gradient, joined, 9000, 11_000, record_length=_LENGTH, length_px=400.0)
    assert [steps[0][1], dict(steps)[10_000.0], steps[-1][1]] == ["#0000fe", "#2800d6", "#4d00b1"]


def _letters_6px_each(annotation, style, font_size, _place_bp, _mark_offset):
    """A stand-in for the drawer's measure (``LabelSizing.extent``): letters 0.6 em wide each, at 1 bp per px, running
    right from a start label's place, left from an end label's, both ways from a centred one's; 10 px deep."""

    from gbdraw.annotations import LabelExtent

    width = 0.6 * len(annotation.label) * font_size
    lead = {"start": 0.0, "end": 1.0}.get(style.label_position, 0.5) * width
    return LabelExtent(-lead, width - lead, 10.0, drawn=True)


@pytest.mark.parametrize(("position", "lanes"), (("center", {"a": 0, "b": 0}), ("start", {"a": 0, "b": 1})))
def test_the_lane_packer_keeps_room_where_the_label_stands(position: str, lanes: dict[str, int]) -> None:
    # A's label is 60 bp wide (10 letters of 6 px at 1 bp per px), wider than its 10 bp region. Centred on bp 105
    # it ends at bp 135, before B at bp 150; standing at A's first base (bp 100) it reaches bp 160, so B takes the
    # next lane (OV-2659).
    from gbdraw.annotations import LabelSizing, layout_annotation_track

    style = RegionAnnotationStyle(label_position=position)
    annotations = tuple(
        ResolvedRegionAnnotation(
            id=identity, set_id="s", record_index=0, segments=(segment,), midpoint_bp=(segment[0] + segment[1]) / 2,
            span_bp=segment[1] - segment[0], label=label, mark="band", lane=None, style=style,
        )
        for identity, segment, label in (("a", (100, 110), "X" * 10), ("b", (150, 160), ""))
    )
    track = layout_annotation_track(
        "regions", "s", annotations, record_lengths={0: 1000}, labels=LabelSizing(reach=_half_stroke_reach, font_size=10.0, extent=_letters_6px_each),
        join_origin=False
    )

    assert {placement.annotation.id: placement.lane for placement in track.placements} == lanes


def _packed_lanes(
    labelled: ResolvedRegionAnnotation, neighbour: tuple[int, int], *, length: int, row_length: int | None = None
) -> dict[str, int]:
    """The lanes of ``labelled`` (its label 25 letters, 150 bp at 1 bp per px) and an unlabelled ``neighbour``."""

    from gbdraw.annotations import LabelSizing, layout_annotation_track

    other = ResolvedRegionAnnotation(
        id="n", set_id="s", record_index=0, segments=(neighbour,), midpoint_bp=sum(neighbour) / 2,
        span_bp=neighbour[1] - neighbour[0], label="", mark="band", lane=None, style=RegionAnnotationStyle(),
        display_parts=None if labelled.display_parts is None else _fragments(None, row_length, (neighbour,)),
    )
    track = layout_annotation_track(
        "regions", "s", (labelled, other), record_lengths={0: row_length or length},
        labels=LabelSizing(reach=_half_stroke_reach, font_size=10.0, extent=_letters_6px_each), join_origin=False,
    )
    return {placement.annotation.id: placement.lane for placement in track.placements}


def _labelled(segments: tuple[tuple[int, int], ...], position: str, **fields: object) -> ResolvedRegionAnnotation:
    return ResolvedRegionAnnotation(
        id="a", set_id="s", record_index=0, segments=segments, midpoint_bp=0.0,
        span_bp=sum(end - start for start, end in segments), label="X" * 25, mark="band", lane=None,
        style=RegionAnnotationStyle(label_position=position), **fields,
    )


def test_the_packer_measures_a_centred_label_on_the_part_it_stands_on() -> None:
    # Two 100 bp parts at the bar's ends: the label stands on 900-1000, where the region starts (PD-OI-119), and
    # its 150 bp run from bp 875 past the bar's end, over the neighbour at 860-890: the neighbour moves on.
    labelled = _labelled(((0, 100), (900, 1000)), "center", crosses_origin=True)

    assert _packed_lanes(labelled, (860, 890), length=1000) == {"a": 0, "n": 1}


def test_the_packer_measures_an_end_label_on_its_row_without_wrapping() -> None:
    # The end label stands at bp 50, the region's last base, and runs left to the bar's start; it does not wrap to
    # the bar's other end, so a neighbour at 920-940 keeps lane 0.
    labelled = _labelled(((0, 50), (950, 1000)), "end", crosses_origin=True)

    assert _packed_lanes(labelled, (920, 940), length=1000) == {"a": 0, "n": 0}


def test_the_packer_measures_a_wrapped_start_label_on_the_part_it_stands_on() -> None:
    # Rows of 4000 bp: the start label stands at bp 1000 on row 0 over the 100 bp part 1000-1100 and runs to bp
    # 1150, over the neighbour at 1120-1160 on the same row; the longer part on row 1 does not make room for it.
    segments = ((1000, 1100), (4100, 7000))
    labelled = _labelled(segments, "start", display_parts=_fragments(None, 4000, segments))

    assert _packed_lanes(labelled, (1120, 1160), length=_LENGTH, row_length=4000) == {"a": 0, "n": 1}


def test_wrapped_rows_without_a_display_start_tie_to_the_part_where_the_region_starts() -> None:
    # Rows of 4000 bp: 0-1000 on the first row, 9000-10000 at the end of the last; equal parts, so the label takes
    # the last row's part, where the region starts (PD-OI-119).
    annotation = _origin_annotation(_fragments(None, 4000))
    place = annotation.label_place(join_origin=False)

    assert (place.part.row, place.bp) == (2, 1500.0)


def test_a_display_seam_ties_to_the_earlier_row_then_further_left() -> None:
    # Region 7000-9000 from display start 8001: two equal parts, 9000-10000 and 0-1000 on the display; the
    # display cut makes them, so the label takes the one further left (PD-OI-119, as before).
    fragments = _fragments(8001, None, ((7000, 9000),))
    annotation = ResolvedRegionAnnotation(
        id="a", set_id="s", record_index=0, segments=((7000, 9000),), midpoint_bp=8000.0, span_bp=2000, label="A",
        mark="band", lane=None, style=RegionAnnotationStyle(), display_parts=fragments,
    )

    assert annotation.label_midpoint_bp(join_origin=False) == 500.0


def test_the_run_follows_the_record_not_the_strand() -> None:
    # A minus-strand region across the origin starts at bp 9000 and ends at bp 1000 like a plus-strand one.
    for strand in ("+", "-"):
        annotation = ResolvedRegionAnnotation(
            id="a", set_id="s", record_index=0, segments=((0, 1000), (9000, _LENGTH)), midpoint_bp=0.0, span_bp=2000,
            label="A", mark="band", lane=None, style=RegionAnnotationStyle(), crosses_origin=True, strand=strand,
        )
        assert annotation.label_anchor("start", join_origin=False).bp == 9000.0
        assert annotation.label_anchor("end", join_origin=False).bp == 1000.0


@pytest.mark.parametrize("position", ("start", "end"))
def test_circular_labels_stand_centred_whatever_their_position(position: str) -> None:
    (text,) = _label_item("circular", 5001, position).iter(f"{_SVG}text")
    x, y = float(text.get("x", "nan")), float(text.get("y", "nan"))

    assert math.degrees(math.atan2(x, -y)) % 360.0 == pytest.approx(180.0, abs=0.5)
