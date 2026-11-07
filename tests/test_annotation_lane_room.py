"""Room between two annotation lanes for a label that would cross the next lane's mark (PD-OI-119, OV-2575).

Where a lane's label, at its place along its region, would cross a mark in the next lane on the label's side, the
two lanes are drawn apart by the room the label takes (its letters and its ``label_offset``); the label stays
beside its own mark in that room. Lanes whose labels cross nothing keep their spacing. The track grows by the room;
a slot whose size the user fixed adds it to its required extent, and its ``overflow`` decides (PD-OI-119).

"Would cross" reads the letters as the drawer measures and draws them at their place. The room is counted per record
and row; records and rows that share the slot and need less keep their lanes where they are.
"""

from __future__ import annotations

import math
import random
import re
from pathlib import Path
from xml.etree import ElementTree

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import gbdraw
import gbdraw.annotations.layout as packer
from gbdraw.annotations import (
    AnnotationTrackParams,
    LabelExtent,
    LabelSizing,
    LaneRoom,
    ResolvedRegionAnnotation,
    layout_annotation_track,
)
from gbdraw.api import AnnotationOptions, AnnotationSet, CoordinateSpan, RegionAnnotation, RegionAnnotationStyle
from gbdraw.exceptions import ValidationError
from gbdraw.io.record_select import RecordSelector
from gbdraw.render.drawers.annotation_marks import label_ink_offsets
from gbdraw.layout.text_geometry import SECTOR_PART_DEG
from gbdraw.render.drawers.circular.annotations import _label_at, _label_ink

_SVG = "{http://www.w3.org/2000/svg}"
_OFFSET = RegionAnnotationStyle().label_offset
# The letters keep the label offset from the painted edge of their mark, half its stroke past the band.
_PAINT = RegionAnnotationStyle().stroke_width / 2
_EDGE = _OFFSET + _PAINT
_HMMT = Path(__file__).parent / "test_inputs" / "HmmtDNA.gbk"


def _record(name: str = "r1") -> SeqRecord:
    return SeqRecord(Seq("ACGT" * 2500), id=name, name=name, annotations={"topology": "circular", "molecule_type": "DNA"})


def _region(
    identity: str, start: int, end: int, label: str, orientation: str = "horizontal", record: str | None = None
) -> RegionAnnotation:
    selector = None if record is None else RecordSelector(record, record, None)
    return RegionAnnotation(
        identity,
        CoordinateSpan(selector, start, end),
        label=label,
        mark="band",
        style=RegionAnnotationStyle(fill="#ffcc00", label_orientation=orientation),
    )


def _svg(
    mode: str,
    regions: tuple[RegionAnnotation, ...],
    slot: str | None = None,
    *,
    records: tuple[SeqRecord, ...] | None = None,
    overrides: dict[str, object] | None = None,
) -> str:
    linear = mode == "linear"
    extra: dict[str, object] = {"config_overrides": overrides} if overrides else {}
    if slot is not None:
        slots = ("features:features", slot)
        extra["tracks"] = gbdraw.LinearTrackOptions(slots=slots) if linear else gbdraw.CircularTrackOptions(slots=slots)
    options = (gbdraw.LinearOptions if linear else gbdraw.CircularOptions)(
        annotations=AnnotationOptions(sets=(AnnotationSet("regions", regions),)), **extra
    )
    draw = gbdraw.draw_linear if linear else gbdraw.draw_circular
    return draw(list(records) if records else _record(), options=options)._drawing.tostring()


def _groups(svg: str) -> dict[str, ElementTree.Element]:
    """Each annotation's group, by record index and id (``"1:c"``; the id alone on the first record)."""

    groups: dict[str, ElementTree.Element] = {}
    for group in ElementTree.fromstring(svg).iter(f"{_SVG}g"):
        identity = group.get("data-gbdraw-annotation-id")
        if identity:
            record = group.get("data-gbdraw-record-index")
            groups.setdefault(identity if record in {None, "0"} else f"{record}:{identity}", group)
    return groups


def _items(mode: str, b_span: tuple[int, int], slot: str | None = None, *, labels: bool = True) -> dict[str, ElementTree.Element]:
    """Two bands on a 10 kb record: ``a`` over 1000-5000 (lane 0) and ``b`` over ``b_span`` (lane 1 where they meet)."""

    regions = tuple(
        _region(identity, start, end, label if labels else "")
        for identity, (start, end), label in (("a", (1000, 5000), "Alpha"), ("b", b_span, "Beta"))
    )
    return _groups(_svg(mode, regions, slot))


def _band(item: ElementTree.Element) -> tuple[float, float]:
    """The top and bottom y of a Linear band."""

    (rect,) = item.iter(f"{_SVG}rect")
    top = float(rect.get("y", "nan"))
    return top, top + float(rect.get("height", "nan"))


def _span_x(item: ElementTree.Element) -> tuple[float, float]:
    """The left and right x of a Linear band."""

    (rect,) = item.iter(f"{_SVG}rect")
    left = float(rect.get("x", "nan"))
    return left, left + float(rect.get("width", "nan"))


def _rows(item: ElementTree.Element) -> dict[str, tuple[float, float]]:
    """The top and bottom y of a wrapped Linear band on each of its rows, in that row's frame."""

    return {
        row.get("data-gbdraw-wrap-row", ""): _band(row)
        for row in item.iter(f"{_SVG}g")
        if row.get("data-gbdraw-wrap-row") is not None
    }


def _radii(item: ElementTree.Element) -> tuple[float, float]:
    """The inner and outer radius of a Circular band."""

    (path,) = item.iter(f"{_SVG}path")
    radii = [float(value) for value in re.findall(r"A ([\d.]+),", path.get("d", ""))]
    return min(radii), max(radii)


def _ink(item: ElementTree.Element) -> tuple[float, float, float, float]:
    """The box a horizontal label's letters paint: left, right, top, bottom."""

    (text,) = item.iter(f"{_SVG}text")
    x, y = float(text.get("x", "nan")), float(text.get("y", "nan"))
    left, right, top, bottom = label_ink_offsets(
        "".join(text.itertext()), text.get("font-family", ""), float(text.get("font-size", "nan")), "middle"
    )
    return x + left, x + right, y + top, y + bottom


def _ink_corners(item: ElementTree.Element) -> list[tuple[float, float]]:
    """The corners of a straight label's ink box, turned by its ``rotate`` about its anchor."""

    (text,) = item.iter(f"{_SVG}text")
    x, y = float(text.get("x", "nan")), float(text.get("y", "nan"))
    turn = re.match(r"rotate\(([-\d.e]+)", text.get("transform", ""))
    angle = math.radians(float(turn.group(1))) if turn else 0.0
    left, right, top, bottom = label_ink_offsets(
        "".join(text.itertext()), text.get("font-family", ""), float(text.get("font-size", "nan")), "middle"
    )
    return [
        (x + dx * math.cos(angle) - dy * math.sin(angle), y + dx * math.sin(angle) + dy * math.cos(angle))
        for dx, dy in ((left, top), (right, top), (right, bottom), (left, bottom))
    ]


def _ink_radii(item: ElementTree.Element) -> tuple[float, float]:
    """How near the circle's centre and how far from it a straight label's letters come."""

    corners = _ink_corners(item)

    def nearest(start: tuple[float, float], end: tuple[float, float]) -> float:
        # The distance from the centre to one edge of the box.
        (x1, y1), (x2, y2) = start, end
        t = max(0.0, min(1.0, -(x1 * (x2 - x1) + y1 * (y2 - y1)) / ((x2 - x1) ** 2 + (y2 - y1) ** 2)))
        return math.hypot(x1 + t * (x2 - x1), y1 + t * (y2 - y1))

    near = min(nearest(corners[index], corners[(index + 1) % 4]) for index in range(4))
    return near, max(math.hypot(cx, cy) for cx, cy in corners)


def _ink_polar(item: ElementTree.Element, record_length: int = 10_000) -> list[tuple[float, float]]:
    """Points over a straight label's ink box on a Circular drawing, each as its radius and its bp."""

    (x0, y0), (x1, y1), _corner, (x3, y3) = _ink_corners(item)
    points = []
    for across in range(201):
        for down in range(21):
            u, v = across / 200, down / 20
            x, y = x0 + u * (x1 - x0) + v * (x3 - x0), y0 + u * (y1 - y0) + v * (y3 - y0)
            points.append((math.hypot(x, y), math.degrees(math.atan2(x, -y)) % 360 / 360 * record_length))
    return points


def test_linear_label_above_gets_room_from_the_lane_it_would_cross() -> None:
    # b (lane 1) labels above, at bp 3000, over a's mark in lane 0: the lanes part by the label's room.
    items = _items("linear", (2500, 3500))
    _a_top, a_bottom = _band(items["a"])
    b_top, _b_bottom = _band(items["b"])
    _left, _right, top, bottom = _ink(items["b"])

    assert a_bottom <= top
    assert bottom == pytest.approx(b_top - _EDGE, abs=0.01)


def test_linear_label_below_gets_room_from_the_lane_it_would_cross() -> None:
    # On a slot below the axis a (lane 0) labels below, at bp 3000, over b's mark in lane 1.
    items = _items("linear", (2500, 3500), "regions:annotations@side=below,set_id=regions")
    _a_top, a_bottom = _band(items["a"])
    b_top, _b_bottom = _band(items["b"])
    _left, _right, top, bottom = _ink(items["a"])

    assert bottom <= b_top
    assert top == pytest.approx(a_bottom + _EDGE, abs=0.01)


@pytest.mark.parametrize("slot", (None, "regions:annotations@side=inside,set_id=regions,r=0.6"))
def test_circular_label_gets_room_from_the_lane_it_would_cross(slot: str | None) -> None:
    # a (lane 0) labels away from the axis at bp 3000, over b's mark in lane 1: outward on an outside slot, inward
    # on an inside one.
    items = _items("circular", (2500, 3500), slot)
    a_inner, a_outer = _radii(items["a"])
    b_inner, b_outer = _radii(items["b"])
    near, far = _ink_radii(items["a"])

    if slot is None:
        assert far <= b_inner
        assert near == pytest.approx(a_outer + _EDGE, abs=0.05)
    else:
        assert near >= b_outer
        assert far == pytest.approx(a_inner - _EDGE, abs=0.05)


def _without_room(monkeypatch: pytest.MonkeyPatch) -> None:
    """Draw as if no label asked for room."""

    monkeypatch.setattr(packer, "_lane_rooms", lambda *_args, **_kwargs: ())


@pytest.mark.parametrize("mode", ("linear", "circular"))
def test_lanes_whose_labels_cross_nothing_keep_their_spacing(mode: str, monkeypatch: pytest.MonkeyPatch) -> None:
    # b (4000-9000) shares 4000-5000 with a, so it takes lane 1. Linear: b's label at bp 6500 is past a's mark;
    # Circular: a's label at bp 3000 is short of b's.
    regions = (_region("a", 1000, 5000, "Alpha"), _region("b", 4000, 9000, "Beta"))
    drawn = _svg(mode, regions)
    labelled, bare = _groups(drawn), _items(mode, (4000, 9000), labels=False)
    shape = _band if mode == "linear" else _radii

    assert [shape(labelled[key]) for key in "ab"] == [shape(bare[key]) for key in "ab"]
    _without_room(monkeypatch)
    assert _svg(mode, regions) == drawn


@pytest.mark.parametrize("orientation", ("horizontal", "radial"))
def test_a_circular_label_whose_letters_cross_nothing_parts_no_lanes(
    orientation: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    # a's label stands at bp 2500, at 3 o'clock: a horizontal or radial label spans only its letters' height along
    # the circle there, short of b's mark (2700-3000) in the next lane.
    regions = (_region("a", 2400, 2600, "LONGLABELNAME", orientation), _region("b", 2700, 3000, ""))
    drawn = _svg("circular", regions)
    items = _groups(drawn)

    assert _radii(items["b"])[0] == pytest.approx(_radii(items["a"])[1])
    _without_room(monkeypatch)
    assert _svg("circular", regions) == drawn


def test_an_outside_slot_the_plan_places_away_from_the_axis_is_measured_where_it_is_drawn(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # The plan places "regions" outside a 150 px slot, about 550 px from the centre. There a's label at bp 5000
    # spans about 4967-5032 bp, short of b (from 5040) in the next lane; at the axis circle it would reach past it.
    sets = (
        AnnotationSet("first", (_region("f1", 1, 9999, ""),)),
        AnnotationSet("regions", (_region("a", 4900, 5100, "iiiiiiii"), _region("b", 5040, 6000, ""))),
    )
    slots = (
        "features:features",
        "regions:annotations@side=outside,set_id=regions",
        "first:annotations@side=outside,set_id=first,w=150px",
    )
    options = gbdraw.CircularOptions(
        annotations=AnnotationOptions(sets=sets), tracks=gbdraw.CircularTrackOptions(slots=slots)
    )
    drawn = gbdraw.draw_circular(_record(), options=options)._drawing.tostring()
    items = _groups(drawn)

    assert _radii(items["a"])[0] > 540.0
    assert _radii(items["b"])[0] == pytest.approx(_radii(items["a"])[1])
    _without_room(monkeypatch)
    assert gbdraw.draw_circular(_record(), options=options)._drawing.tostring() == drawn


def test_radial_labels_that_grow_the_circle_measure_an_outside_slot_where_it_is_drawn(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # Radial feature labels grow the circle, which moves the outside slot out by about 20 px. Where it is drawn,
    # a's letters at its middle span bp 7968-8604 (OV-2663), short of b (from 8624) in the next lane; where the slot
    # was before the growth they would reach past it.
    record = gbdraw.read_genbank(str(_HMMT))[0]
    regions = (_region("a", 7456, 9112, "Alpha region label"), _region("b", 8624, 9224, ""))
    options = gbdraw.CircularOptions(
        annotations=AnnotationOptions(sets=(AnnotationSet("regions", regions),)),
        tracks=gbdraw.CircularTrackOptions(slots=("features:features", "regions:annotations@side=outside,set_id=regions")),
        config_overrides={"labels.circular.placement": "radial", "labels.circular.scope": "both"},
    )
    drawn = gbdraw.draw_circular(record, options=options)._drawing.tostring()
    items = _groups(drawn)

    assert _radii(items["b"])[0] == pytest.approx(_radii(items["a"])[1])
    _without_room(monkeypatch)
    assert gbdraw.draw_circular(record, options=options)._drawing.tostring() == drawn


def test_the_definition_scan_measures_each_label_room_once(monkeypatch: pytest.MonkeyPatch) -> None:
    # 150 labelled regions on an inside slot need lane rooms that carry the slot over the centre definition, so the
    # label fit check scans smaller definition fonts, each laying the tracks out again (OV-2664). A slot's label room
    # is measured once for each band, axis and side it is asked for, however often the layouts ask.
    import gbdraw.diagrams.circular.assemble as assemble_module

    record = gbdraw.read_genbank(str(_HMMT))[0]
    length = len(record.seq)
    rng = random.Random(1)
    for _ in range(600):  # Seeded random regions: 150 labelled ones inside the slot.
        rng.randint(1, 2)
    regions = tuple(
        RegionAnnotation(f"q{index}", CoordinateSpan(None, start, start + rng.randint(50, 400)), label=f"q{index}")
        for index in range(150)
        for start in (rng.randint(1, length - 400),)
    )
    measured: list[tuple[object, ...]] = []
    original = assemble_module.circular_annotation_label_room

    def counting(track: object, **kwargs: object) -> float:
        # The track itself is kept, so no other track takes its id while the drawing runs.
        measured.append((track, kwargs["band"], kwargs["side"], tuple(kwargs["rings"]), kwargs["params"], kwargs["covers_slots"]))
        return original(track, **kwargs)

    monkeypatch.setattr(assemble_module, "circular_annotation_label_room", counting)
    options = gbdraw.CircularOptions(
        annotations=AnnotationOptions(sets=(AnnotationSet("i1", regions),)),
        tracks=gbdraw.CircularTrackOptions(slots=("features:features", "i1:annotations@side=inside,set_id=i1")),
    )
    with pytest.warns(UserWarning, match="definition lines do not fit"):
        gbdraw.draw_circular(record, options=options)

    distinct = {(id(track), *rest) for track, *rest in measured}
    assert distinct and len(measured) == len(distinct)


def test_a_curved_label_pushed_to_the_centre_is_measured_within_one_turn() -> None:
    # Lane rooms can carry an inside lane's mark edge past the centre (OV-2665). There a curved label is cut to its
    # ellipsis, whose letters stand on the baseline, so its radius nears zero; its letters still cover at most one
    # turn of the circle, so its sector is cut into at most one turn of parts. (The extent folds angles into one turn
    # either way, so the parts count is what a missing limit shows.)
    style = RegionAnnotationStyle()
    label = _label_at(
        0.0,
        annotation_id="",
        text="ace",
        record_length=10_000,
        mark_edge=1.0,
        outward=False,
        gap=float(style.label_offset),
        font_size=12.0,
        color="",
        orientation=style.label_orientation,
        ink=_label_ink("Arial", 12.0, style.label_orientation),
        rings=(),
    )

    assert len(label.parts) <= math.ceil(360.0 / SECTOR_PART_DEG)


def test_a_wide_linear_label_gets_room_where_its_letters_cross_the_next_lane() -> None:
    # b (lane 1) labels above at bp 3350; its wide letters reach over c (to 3050) and a (from 3720) in lane 0.
    regions = (_region("c", 2000, 3050, ""), _region("b", 3000, 3700, "WWWWWWWW"), _region("a", 3720, 5000, ""))
    items = _groups(_svg("linear", regions))
    left, right, top, bottom = _ink(items["b"])
    b_top, _b_bottom = _band(items["b"])

    assert left < _span_x(items["c"])[1] and right > _span_x(items["a"])[0]
    assert top >= max(_band(items["c"])[1], _band(items["a"])[1])
    assert bottom == pytest.approx(b_top - _EDGE, abs=0.01)


@pytest.mark.parametrize(
    ("orientation", "radius", "b_span"),
    (
        *((orientation, radius, (2500, 3500)) for radius in (",r=0.6", "") for orientation in ("tangent", "horizontal")),
        # A slot centred on its radius widens both ways by the room: its lanes' edge moves out, where the letters
        # would cross nothing, so a room that a gap needed where the slot was drawn before stays.
        ("tangent", ",r=0.6", (3640, 4440)),
        ("tangent", ",r=0.6", (3640, 3680)),
        ("horizontal", ",r=0.6", (3500, 3900)),
    ),
)
def test_an_inside_circular_label_is_measured_on_the_slots_side_and_radius(
    orientation: str, radius: str, b_span: tuple[int, int]
) -> None:
    # On an inside slot at 0.6 of the radius, or where the plan places it (about 330 px from the centre, well inside
    # the axis circle), a (lane 0) labels inward at bp 3000 over b's mark in lane 1.
    regions = (_region("a", 1000, 5000, "A rather long region label text", orientation), _region("b", *b_span, ""))
    items = _groups(_svg("circular", regions, f"regions:annotations@side=inside,set_id=regions{radius}"))
    a_inner, _a_outer = _radii(items["a"])
    b_inner, b_outer = _radii(items["b"])
    _near, far = _ink_radii(items["a"])

    # A Circular band fills its lane: the room leaves the letters right at the next mark's paint, or past it.
    assert not [
        (radius, bp)
        for radius, bp in _ink_polar(items["a"])
        if b_inner - _PAINT < radius < b_outer + _PAINT - 0.05 and b_span[0] <= bp <= b_span[1]
    ]
    assert far == pytest.approx(a_inner - _EDGE, abs=0.05)


def test_a_record_whose_labels_cross_nothing_keeps_its_lanes_beside_one_that_needs_room() -> None:
    # r1: b's label (lane 1, above) crosses a's mark, so r1's lanes part and the slot grows toward the axis by the
    # room (b moves down by it). r2 shares the slot: d's label stands past c's mark, so r2's lanes keep their
    # spacing and their places beside the axis (they move down with the slot's lower edge), the room left empty
    # above them on the labels' side.
    def drawn(beta: str) -> dict[str, ElementTree.Element]:
        regions = (
            _region("a", 1000, 5000, "Alpha", record="r1"),
            _region("b", 2500, 3500, beta, record="r1"),
            _region("c", 1000, 5000, "Gamma", record="r2"),
            _region("d", 4000, 9000, "Delta", record="r2"),
        )
        return _groups(_svg("linear", regions, records=(_record("r1"), _record("r2"))))

    parted, plain = drawn("Beta"), drawn("")
    room = _band(parted["b"])[0] - _band(plain["b"])[0]

    assert room > 0.0
    assert _band(parted["a"]) == _band(plain["a"])
    assert [y for key in ("1:c", "1:d") for y in _band(parted[key])] == pytest.approx(
        [y + room for key in ("1:c", "1:d") for y in _band(plain[key])], abs=1e-9
    )


def test_a_row_whose_labels_cross_nothing_keeps_its_lanes_beside_one_that_needs_room() -> None:
    # Rows of 4 kb. Row 0: b's label (lane 1, above) crosses a's mark, so the lanes part there and the slot grows
    # toward the axis on every row. Row 1: d (lane 1) has no label and c's stands above lane 0, so its lanes keep
    # their spacing and their places beside the axis.
    def drawn(beta: str) -> dict[str, ElementTree.Element]:
        regions = (
            _region("a", 1000, 3500, "Alpha"),
            _region("b", 2000, 3000, beta),
            _region("c", 4500, 7500, "Gamma"),
            _region("d", 6000, 7800, ""),
        )
        return _groups(_svg("linear", regions, overrides={"canvas.linear.wrap_length": 4000}))

    parted, plain = drawn("Beta"), drawn("")
    room = _rows(parted["b"])["0"][0] - _rows(plain["b"])["0"][0]

    assert room > 0.0
    assert [y for key in "cd" for y in _rows(parted[key])["1"]] == pytest.approx(
        [y + room for key in "cd" for y in _rows(plain[key])["1"]], abs=1e-9
    )


def test_lanes_without_room_keep_their_places_to_the_last_bit() -> None:
    # Four unlabelled bands in a slot 83.7 px high: each lane's middle is the padding and (lane + 1/2) lane heights.
    spans = ((1000, 5000), (2000, 6000), (3000, 7000), (3500, 7500))
    regions = tuple(_region(identity, start, end, "") for identity, (start, end) in zip("abcd", spans))
    items = _groups(_svg("linear", regions, "regions:annotations@side=above,set_id=regions,h=83.7px"))
    padding = AnnotationTrackParams(set_id="regions").padding_px
    lane = (83.7 - 2.0 * padding) / 4

    for index, key in enumerate("abcd"):
        centre = padding + (index + 0.5) * lane
        assert _band(items[key])[0] == centre - 0.5 * 0.84 * lane


def _packed(
    *,
    overflow: str = "error",
    available: float | None = None,
    above: bool = False,
    join_origin: bool = False,
    b_stroke: float | None = None,
    **params: object,
):
    """``a`` over 100-500 labelled ``A`` and ``b`` over 200-300 (stroke ``b_stroke``), with letters 10 bp each way
    of their place and 10 px deep."""

    styles = {"a": RegionAnnotationStyle(), "b": RegionAnnotationStyle(stroke_width=b_stroke) if b_stroke else RegionAnnotationStyle()}
    annotations = tuple(
        ResolvedRegionAnnotation(
            id=identity, set_id="s", record_index=0, segments=(segment,), midpoint_bp=(segment[0] + segment[1]) / 2,
            span_bp=segment[1] - segment[0], label=label, mark="band", lane=None, style=styles[identity],
        )
        for identity, segment, label in (("a", (100, 500), "A"), ("b", (200, 300), "B" if above else ""))
    )
    return layout_annotation_track(
        "regions",
        "s",
        annotations,
        record_lengths={0: 1000},
        params=AnnotationTrackParams(set_id="s", overflow=overflow, **params),
        available_extent_px=available,
        labels=LabelSizing(font_size=10.0, extent=lambda *_: LabelExtent(-10.0, 10.0, 10.0), labels_above=above),
        join_origin=join_origin,
    )


# Two lanes of 14 px, a lane gap of 3 px and padding of 2 px on each side. The room holds the letters (10 px), the
# label offset and half of each of the two strokes: the label's own mark past its lane, and the mark it would cross.
_BASE = 2 * 2.0 + 2 * 14.0 + 3.0
_ROOM = 10.0 + _OFFSET + RegionAnnotationStyle().stroke_width


@pytest.mark.parametrize(("above", "join_origin"), ((False, False), (True, False), (False, True)))
def test_the_track_grows_by_the_room_between_the_two_lanes(above: bool, join_origin: bool) -> None:
    track = _packed(above=above, join_origin=join_origin)

    assert {item.annotation.id: item.lane for item in track.placements} == {"a": 0, "b": 1}
    assert track.lane_rooms == (LaneRoom(record_index=0, row=0, gap=0, px=_ROOM),)
    assert track.required_extent_px == pytest.approx(_BASE + _ROOM)


def test_a_span_over_other_slots_gets_no_room() -> None:
    track = _packed(anchor_slot="features", cover_anchor=True)

    assert track.lane_rooms == ()
    assert track.required_extent_px == pytest.approx(_BASE)


def test_a_fixed_slot_with_overflow_error_rejects_the_room() -> None:
    with pytest.raises(ValidationError, match="requires"):
        _packed(available=_BASE + 1.0)


def test_a_fixed_slot_with_overflow_compress_narrows_the_lane_gap() -> None:
    track = _packed(overflow="compress", available=_BASE + _ROOM - 2.0)

    assert track.lane_gap_px == pytest.approx(1.0)
    assert track.required_extent_px == pytest.approx(_BASE + _ROOM - 2.0)
    with pytest.raises(ValidationError, match="cannot fit"):
        _packed(overflow="compress", available=_BASE + _ROOM - 4.0)


def test_a_fixed_slot_with_overflow_clip_clips() -> None:
    track = _packed(overflow="clip", available=_BASE + 1.0)

    assert track.clipped
    assert track.lane_rooms == (LaneRoom(record_index=0, row=0, gap=0, px=_ROOM),)


def test_the_room_clears_a_wider_stroke_of_the_mark_it_would_cross() -> None:
    track = _packed(b_stroke=4.0)
    room = 10.0 + _OFFSET + 0.5 * (RegionAnnotationStyle().stroke_width + 4.0)

    assert track.lane_rooms == (LaneRoom(record_index=0, row=0, gap=0, px=room),)
    assert track.required_extent_px == pytest.approx(_BASE + room)


def test_the_slot_grows_by_the_largest_room_of_a_record_not_their_sum() -> None:
    # Record 0 needs room between its lanes; record 1 has the same lanes with no label.
    annotations = tuple(
        ResolvedRegionAnnotation(
            id=f"{identity}{record}", set_id="s", record_index=record, segments=(segment,),
            midpoint_bp=(segment[0] + segment[1]) / 2, span_bp=segment[1] - segment[0],
            label=label if record == 0 else "", mark="band", lane=None, style=RegionAnnotationStyle(),
        )
        for record in (0, 1)
        for identity, segment, label in (("a", (100, 500), "A"), ("b", (200, 300), ""))
    )
    track = layout_annotation_track(
        "regions",
        "s",
        annotations,
        record_lengths={0: 1000, 1: 1000},
        params=AnnotationTrackParams(set_id="s"),
        labels=LabelSizing(font_size=10.0, extent=lambda *_: LabelExtent(-10.0, 10.0, 10.0)),
        join_origin=False,
    )

    assert track.lane_rooms == (LaneRoom(record_index=0, row=0, gap=0, px=_ROOM),)
    assert track.lane_room_px == pytest.approx(_ROOM)
    assert track.required_extent_px == pytest.approx(_BASE + _ROOM)
    assert (track.room_before(0, 0, 1), track.room_before(1, 0, 1)) == (_ROOM, 0.0)
