"""Room between two annotation lanes for a label that would cross the next lane's mark (PD-OI-119, OV-2575).

Where a lane's label, at its place along its region, would cross a mark in the next lane on the label's side, the
two lanes are drawn apart by the room the label takes (its letters and its ``label_offset``); the label stays
beside its own mark in that room. Lanes whose labels cross nothing keep their spacing. The track grows by the room;
a slot whose size the user fixed adds it to its required extent, and its ``overflow`` decides (PD-OI-119); with
``overflow=clip`` it adds only the rooms that keep its lanes within that size and draws every mark; a label without
room, or one whose letters that size would cut, is left out, and gbdraw warns (LBL-ROOM-CLIP-FIT-OMIT). An inside
Circular slot stops its rooms at the centre the same way (LBL-ROOM-STOP-AT-CENTRE-OMIT).

"Would cross" reads the letters as the drawer measures and draws them at their place. The room is counted per record
and row; records and rows that share the slot and need less keep their lanes where they are.
"""

from __future__ import annotations

import math
import random
import re
import warnings
from dataclasses import replace
from pathlib import Path
from typing import NamedTuple
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
from gbdraw.render.drawers.circular.annotations import CircularAnnotationLabel, _label_at

_SVG = "{http://www.w3.org/2000/svg}"
_OFFSET = RegionAnnotationStyle().label_offset
# The letters keep the label offset from the painted edge of their mark, half its stroke past the band.
_PAINT = RegionAnnotationStyle().stroke_width / 2
_EDGE = _OFFSET + _PAINT
_HMMT = Path(__file__).parent / "test_inputs" / "HmmtDNA.gbk"


def _half_stroke_reach(_annotation, style, _lane_width: float) -> float:
    # The packer tests' own mark paint (``LabelSizing.reach``): half the stroke past the lane's edge.
    return 0.5 * float(style.stroke_width)



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
        (("".join(text.itertext()), "normal"),), text.get("font-family", ""), float(text.get("font-size", "nan")), "middle"
    )
    return x + left, x + right, y + top, y + bottom


def _ink_corners(item: ElementTree.Element) -> list[tuple[float, float]]:
    """The corners of a straight label's ink box, turned by its ``rotate`` about its anchor."""

    (text,) = item.iter(f"{_SVG}text")
    x, y = float(text.get("x", "nan")), float(text.get("y", "nan"))
    turn = re.match(r"rotate\(([-\d.e]+)", text.get("transform", ""))
    angle = math.radians(float(turn.group(1))) if turn else 0.0
    left, right, top, bottom = label_ink_offsets(
        (("".join(text.itertext()), "normal"),), text.get("font-family", ""), float(text.get("font-size", "nan")), "middle"
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
    # b (lane 1) labels above, at bp 3000, over a's mark in lane 0: the lanes part by the label's room, which keeps
    # the label offset on each side of the letters (from b's mark and before a's).
    items = _items("linear", (2500, 3500))
    _a_top, a_bottom = _band(items["a"])
    b_top, _b_bottom = _band(items["b"])
    _left, _right, top, bottom = _ink(items["b"])

    # The band's fill is narrower than its lane, so a's painted edge lies a little further off than the room counts.
    assert top >= a_bottom + _EDGE - 0.01
    assert bottom == pytest.approx(b_top - _EDGE, abs=0.01)


def test_linear_label_below_gets_room_from_the_lane_it_would_cross() -> None:
    # On a slot below the axis a (lane 0) labels below, at bp 3000, over b's mark in lane 1.
    items = _items("linear", (2500, 3500), "regions:annotations@side=below,set_id=regions")
    _a_top, a_bottom = _band(items["a"])
    b_top, _b_bottom = _band(items["b"])
    _left, _right, top, bottom = _ink(items["a"])

    assert bottom <= b_top - _EDGE + 0.01
    assert top == pytest.approx(a_bottom + _EDGE, abs=0.01)


@pytest.mark.parametrize("slot", (None, "regions:annotations@side=inside,set_id=regions,r=0.6"))
def test_circular_label_gets_room_from_the_lane_it_would_cross(slot: str | None) -> None:
    # a (lane 0) labels away from the axis at bp 3000, over b's mark in lane 1: outward on an outside slot, inward
    # on an inside one. The letters keep the label offset from both marks.
    items = _items("circular", (2500, 3500), slot)
    a_inner, a_outer = _radii(items["a"])
    b_inner, b_outer = _radii(items["b"])
    near, far = _ink_radii(items["a"])

    if slot is None:
        assert far == pytest.approx(b_inner - _EDGE, abs=0.05)
        assert near == pytest.approx(a_outer + _EDGE, abs=0.05)
    else:
        assert near == pytest.approx(b_outer + _EDGE, abs=0.05)
        assert far == pytest.approx(a_inner - _EDGE, abs=0.05)


def _without_room(monkeypatch: pytest.MonkeyPatch) -> None:
    """Draw as if no label asked for room."""

    monkeypatch.setattr(packer, "_lane_rooms", lambda *_args, **_kwargs: packer._LaneRooms(()))


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
    # the circle there, short of b's mark (2550-3000) in the next lane.
    regions = (_region("a", 2400, 2600, "LONGLABELNAME", orientation), _region("b", 2550, 3000, ""))
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


def _inside_slot_over_the_definition() -> tuple[SeqRecord, gbdraw.CircularOptions]:
    """200 labelled regions on an inside slot: their lane rooms carry the slot over the centre definition, so the
    label fit check scans smaller definition fonts, each verified with a full assembly (OV-2664)."""

    record = gbdraw.read_genbank(str(_HMMT))[0]
    length = len(record.seq)
    rng = random.Random(1)
    for _ in range(600):  # Seeded random regions: 200 labelled ones inside the slot.
        rng.randint(1, 2)
    regions = tuple(
        RegionAnnotation(f"q{index}", CoordinateSpan(None, start, start + rng.randint(50, 400)), label=f"q{index}")
        for index in range(200)
        for start in (rng.randint(1, length - 400),)
    )
    options = gbdraw.CircularOptions(
        annotations=AnnotationOptions(sets=(AnnotationSet("i1", regions),)),
        tracks=gbdraw.CircularTrackOptions(slots=("features:features", "i1:annotations@side=inside,set_id=i1")),
    )
    return record, options


def test_the_definition_scan_measures_each_label_room_once(monkeypatch: pytest.MonkeyPatch) -> None:
    # The scan lays the tracks out again for each smaller definition font (OV-2664). A slot's label room is measured
    # once for each band, axis and side it is asked for, however often the layouts ask.
    import gbdraw.diagrams.circular.assemble as assemble_module

    record, options = _inside_slot_over_the_definition()
    measured: list[tuple[object, ...]] = []
    original = assemble_module.circular_annotation_label_room

    def counting(track: object, **kwargs: object) -> float:
        # The track itself is kept, so no other track takes its id while the drawing runs.
        measured.append((track, kwargs["band"], kwargs["side"], tuple(kwargs["rings"]), kwargs["params"], kwargs["covers_slots"]))
        return original(track, **kwargs)

    monkeypatch.setattr(assemble_module, "circular_annotation_label_room", counting)
    with pytest.warns(UserWarning, match="definition lines"):
        gbdraw.draw_circular(record, options=options)

    distinct = {(id(track), *rest) for track, *rest in measured}
    assert distinct and len(measured) == len(distinct)


def test_the_definition_scan_packs_and_places_each_slot_once_for_the_same_inputs(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    # The scan's assemblies differ in the definition font only (OV-2664): a slot packed from the same edge with the
    # same rooms is packed once, and its labels drawn at the same radii among the same obstacles are placed once,
    # however many rows the check verifies. Placing a track's labels counts each label against one obstacle index
    # that takes every label placed before it, not an index made again for each label. (The lane packer places a
    # track's labels among no other obstacles where it packs them, ``circular_annotation_lanes_meeting``.)
    import gbdraw.diagrams.circular.assemble as assemble_module
    import gbdraw.render.drawers.circular.annotations as drawer
    from gbdraw.api.diagram import label_fit_report

    record, options = _inside_slot_over_the_definition()
    packed: list[tuple[object, ...]] = []
    placed: list[tuple[object, ...]] = []
    indexes: list[int] = []
    pack, place, index_init = (
        assemble_module.layout_annotation_track,
        assemble_module.circular_annotation_labels,
        drawer._ObstacleIndex.__init__,
    )

    def counting_pack(slot_id: str, set_id: str, annotations: object, **kwargs: object) -> object:
        labels = kwargs["labels"]
        assert isinstance(labels, LabelSizing)
        packed.append((slot_id, tuple(sorted(labels.extent.keywords.items())), kwargs["room_bounds"]))
        return pack(slot_id, set_id, annotations, **kwargs)

    placing: list[bool] = []

    def counting_place(track: object, **kwargs: object) -> object:
        if kwargs.get("obstacles") is not None:
            placed.append((kwargs["inner_radius_px"], kwargs["outer_radius_px"], kwargs["obstacles"]))
        placing.append(True)
        try:
            return place(track, **kwargs)
        finally:
            placing.pop()

    def counting_index(self: object, obstacles: object) -> None:
        if placing:
            indexes.append(1)
        index_init(self, obstacles)

    monkeypatch.setattr(assemble_module, "layout_annotation_track", counting_pack)
    monkeypatch.setattr(assemble_module, "circular_annotation_labels", counting_place)
    monkeypatch.setattr(drawer._ObstacleIndex, "__init__", counting_index)
    with pytest.warns(UserWarning, match="definition lines"):
        drawing = gbdraw.draw_circular(record, options=options)

    report = label_fit_report(drawing._drawing)
    assert report is not None and any(row.kind == "smaller_definition_font" for row in report.rows)
    assert packed and len(packed) == len(set(packed))
    assert placed and len(placed) == len(set(placed))
    assert len(indexes) == len(placed)


def test_a_curved_label_pushed_to_the_centre_is_measured_within_one_turn() -> None:
    # Lane rooms can carry an inside lane's mark edge past the centre (OV-2665). There a curved label is cut to its
    # ellipsis, whose letters stand on the baseline, so its radius nears zero; its letters still cover at most one
    # turn of the circle, so its sector is cut into at most one turn of parts. (The extent folds angles into one turn
    # either way, so the parts count is what a missing limit shows.)
    style = RegionAnnotationStyle()
    label = _label_at(
        0.0,
        annotation_id="",
        runs=(("ace", "normal"),),
        record_length=10_000,
        mark_edge=1.0,
        outward=False,
        gap=float(style.label_offset),
        font_size=12.0,
        color="",
        orientation=style.label_orientation,
        font_family="Arial",
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
        labels=LabelSizing(reach=_half_stroke_reach, font_size=10.0, extent=lambda *_: LabelExtent(-10.0, 10.0, 10.0, drawn=True), labels_above=above),
        join_origin=join_origin,
    )


# Two lanes of 14 px, a lane gap of 3 px and padding of 2 px on each side. The room holds the letters (10 px), the
# label offset on each side of them (beside the label's own mark and before the mark it would cross) and half of
# each of the two strokes: the label's own mark past its lane, and the mark it would cross.
_BASE = 2 * 2.0 + 2 * 14.0 + 3.0
_ROOM = 10.0 + 2.0 * _OFFSET + RegionAnnotationStyle().stroke_width


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


def _nested_regions(count: int, step: int, size: float, orientation: str) -> tuple[RegionAnnotation, ...]:
    """``count`` nested bands from 1501-3500 inward by ``step`` bp, labelled at ``size`` pt in ``orientation``."""

    return tuple(
        RegionAnnotation(
            f"r{index}",
            CoordinateSpan(None, 1501 + step * index, 3500 - step * index),
            label=f"Region {index + 1} label text",
            mark="band",
            style=RegionAnnotationStyle(fill="#ffcc00", label_orientation=orientation, label_font_size=size),
        )
        for index in range(count)
    )


@pytest.mark.parametrize(
    ("mode", "regions", "slot"),
    (
        (
            "linear",
            tuple(_region(f"r{index}", 1501 + 120 * index, 3500 - 120 * index, f"Region {index + 1} label text") for index in range(6)),
            "regions:annotations@side=above,set_id=regions,h={}px,overflow=error",
        ),
        # OV-3024: an inside slot whose curved labels are deeper nearer the centre. The planner's packing at the axis
        # circle, further from the centre than where the slot is drawn, named a size 0.12 px too small.
        ("circular", _nested_regions(7, 100, 24, "arc"), "regions:annotations@side=inside,set_id=regions,w={}px,overflow=error"),
    ),
)
def test_a_slot_that_raises_names_the_smallest_size_that_draws_its_lanes(
    mode: str, regions: tuple[RegionAnnotation, ...], slot: str
) -> None:
    # Nested bands in a slot 10 px wide: the lanes and their rooms do not fit, so the slot raises without drawing. Its
    # rooms are measured where the lanes keep their extent, as at the size it names, not where a narrower size would
    # squeeze them, and where the slot is drawn: the slot draws at that size and not 0.01 px below it.
    with pytest.raises(ValidationError, match="requires") as raised:
        _svg(mode, regions, slot.format(10))
    match = re.search(r"requires ([0-9.]+)px", str(raised.value))
    assert match is not None
    named = float(match.group(1))

    _svg(mode, regions, slot.format(f"{named:.2f}"))
    with pytest.raises(ValidationError, match="requires"):
        _svg(mode, regions, slot.format(f"{named - 0.01:.2f}"))


@pytest.mark.parametrize(
    ("overflow", "refusal"),
    (
        ("error", r"requires [0-9.]+px but only ([0-9.]+)px is available before the centre of the circle"),
        ("compress", r"cannot fit 7 lanes in ([0-9.]+)px before the centre of the circle"),
    ),
)
def test_an_inside_slot_whose_set_width_passes_the_centre_is_checked_against_the_room_before_it(
    overflow: str, refusal: str
) -> None:
    # OV-2990: seven nested bands labelled at 12 pt in an inside slot 640 px wide, past the centre. The drawer gives
    # the slot only the room before the centre; its lanes and their rooms take more, so the slot's overflow refuses
    # it there, naming that room, rather than drawing the bands 0 px wide with their labels stacked at the centre.
    slot = f"regions:annotations@side=inside,set_id=regions,w=640px,overflow={overflow}"
    with pytest.raises(ValidationError) as raised:
        _svg("circular", _nested_regions(7, 100, 12, "horizontal"), slot)
    match = re.search(refusal, str(raised.value))
    assert match is not None, str(raised.value)
    assert float(match.group(1)) < 640.0


def test_a_fixed_slot_with_overflow_compress_narrows_the_lane_gap() -> None:
    track = _packed(overflow="compress", available=_BASE + _ROOM - 2.0)

    assert track.lane_gap_px == pytest.approx(1.0)
    assert track.required_extent_px == pytest.approx(_BASE + _ROOM - 2.0)
    with pytest.raises(ValidationError, match="cannot fit"):
        _packed(overflow="compress", available=_BASE + _ROOM - 4.0)


def test_a_clipped_slot_of_set_size_adds_no_room_past_that_size() -> None:
    # LBL-ROOM-CLIP-FIT-OMIT (OV-2672): the room would carry lane 1 past the set size, which the clip would cut. It is
    # not added: the lanes keep their spacing within the size, and A's label is left without room.
    track = _packed(overflow="clip", available=_BASE + 1.0)

    assert not track.clipped
    assert track.lane_rooms == ()
    assert track.required_extent_px == pytest.approx(_BASE)
    assert (track.without_room, track.rooms_stop, track.labels_without_room(0)) == (((0, "a"),), "clip", 1)


def test_a_clipped_slot_of_set_size_keeps_the_room_that_fits() -> None:
    track = _packed(overflow="clip", available=_BASE + _ROOM + 1.0)

    assert not track.clipped
    assert track.lane_rooms == (LaneRoom(record_index=0, row=0, gap=0, px=_ROOM),)
    assert track.without_room == ()


def _clipped(annotations, letters, available: float):
    return layout_annotation_track(
        "regions", "s", annotations, record_lengths={0: 1000},
        params=AnnotationTrackParams(set_id="s", overflow="clip"), available_extent_px=available,
        labels=LabelSizing(reach=_half_stroke_reach, font_size=10.0, extent=letters), join_origin=False,
    )


def test_a_clipped_slot_keeps_the_room_of_a_label_that_fits_where_another_of_its_gap_does_not() -> None:
    # LBL-ROOM-CLIP-FIT-OMIT, as at the centre: A needs 69.5 px, B 17.5 px; the slot of 55 px has 20 px past its
    # lanes. The gap takes B's room, and only A is left without room.
    annotations = tuple(
        ResolvedRegionAnnotation(
            id=identity, set_id="s", record_index=0, segments=(segment,), midpoint_bp=(segment[0] + segment[1]) / 2,
            span_bp=segment[1] - segment[0], label=label, mark="band", lane=lane, style=RegionAnnotationStyle(),
        )
        for identity, segment, label, lane in (
            ("A", (100, 200), "LONG", 0), ("B", (600, 700), "x", 0), ("C", (100, 200), "", 1), ("D", (600, 700), "", 1)
        )
    )

    def letters(annotation, _style, _font_size, _place_bp, _mark_offset):
        return LabelExtent(-10.0, 10.0, 60.0 if annotation.label == "LONG" else 8.0, drawn=True)

    track = _clipped(annotations, letters, _BASE + 20.0)

    assert track.lane_rooms == (LaneRoom(0, 0, 0, pytest.approx(8.0 + 2.0 * _OFFSET + RegionAnnotationStyle().stroke_width)),)
    assert track.without_room == ((0, "A"),)
    assert not track.clipped


def test_a_clipped_linear_slot_measures_its_lanes_where_the_drawer_draws_them() -> None:
    # Three lanes pass the set size of 51 px (52 px with their gaps), so the slot clips and the drawer draws each lane
    # 14 px high, with no gap. A's letters (lane 0) pass b's mark in lane 1 and meet c's in lane 2: the room before
    # lane 2 counts lane 1 at the 14 px it is drawn (18 px of letters, offsets and strokes less 14 px), and it fits in
    # the 5 px the touching lanes leave. Measured at the lanes' closed form (16 px each) it was 2 px short.
    style = RegionAnnotationStyle()
    annotations = tuple(
        ResolvedRegionAnnotation(
            id=identity, set_id="s", record_index=0, segments=(segment,), midpoint_bp=(segment[0] + segment[1]) / 2,
            span_bp=segment[1] - segment[0], label=label, mark="band", lane=lane, style=style,
        )
        for identity, segment, label, lane in (("a", (100, 500), "A", 0), ("b", (600, 700), "", 1), ("c", (250, 350), "", 2))
    )
    depth = 18.0 - 2.0 * style.label_offset - style.stroke_width
    letters = lambda *_: LabelExtent(-10.0, 10.0, depth, drawn=True)  # noqa: E731

    track = _clipped(annotations, letters, 51.0)

    assert track.clipped
    assert track.lane_rooms == (LaneRoom(0, 0, 1, pytest.approx(4.0)),)
    # One pixel less and the room no longer fits: the lanes touch up to the set size, and A is left without room.
    assert _clipped(annotations, letters, 49.0).lane_rooms == ()


_NESTED = (
    _region("a", 1001, 7000, "Region A with a long name"),
    _region("b", 2501, 9000, "B label here"),
    _region("c", 3001, 6000, "C"),
)


def _marked(item: ElementTree.Element) -> bool:
    """Whether an annotation's group paints its mark (a band's shape, a bracket's lines)."""

    return any(True for tag in ("rect", "path", "line") for _ in item.iter(f"{_SVG}{tag}"))


def _labelled(item: ElementTree.Element) -> bool:
    return any(True for _ in item.iter(f"{_SVG}text"))


def _clip_bounds(svg: str, linear: bool, size: int) -> tuple[float, float] | None:
    """Where a clipped slot's clip lets its marks paint: Linear, from the slot's top to its set size; Circular, between
    the radii of its clip ring. None when the slot draws no clip."""

    clips = list(ElementTree.fromstring(svg).iter(f"{_SVG}clipPath"))
    if not clips:
        return None
    if linear:
        return 0.0, float(size)
    radii = [float(value) for value in re.findall(r"A ([\d.]+),", "".join(item.get("d", "") for item in clips[0]))]
    return min(radii), max(radii)


@pytest.mark.parametrize(
    ("mode", "size", "drawn", "partly_clipped"),
    (("linear", 40, "", "c"), ("linear", 60, "a", ""), ("circular", 40, "", "c"), ("circular", 60, "c", "")),
)
def test_a_clipped_slot_of_set_size_draws_every_mark_and_leaves_out_the_labels_without_room(
    mode: str, size: int, drawn: str, partly_clipped: str
) -> None:
    # LBL-ROOM-CLIP-FIT-OMIT (OV-2672, Q12 X): three nested regions in a slot of set size with overflow=clip. Every
    # mark is drawn and starts within the set size; no room carries one past it. At 40 px the three lanes alone pass
    # the set size, so the clip cuts the outer mark (c) as on dev; at 60 px every mark is whole. A label whose room
    # does not fit, or whose letters the edge of the set size would cut, is not drawn; gbdraw warns with the count of
    # the labels it left out. At 40 px no label fits; at 60 px one does: A's above the outer Linear lane, C's outside
    # the outer Circular lane.
    from gbdraw.api.diagram import label_fit_report
    from gbdraw.labels.fit import label_fit_notice

    linear = mode == "linear"
    side = "side=above,h" if linear else "side=outside,w"
    slot = f"regions:annotations@{side}={size}px,set_id=regions,overflow=clip"
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        drawing = (gbdraw.draw_linear if linear else gbdraw.draw_circular)(
            [_record()] if linear else _record(),
            options=(gbdraw.LinearOptions if linear else gbdraw.CircularOptions)(
                annotations=AnnotationOptions(sets=(AnnotationSet("regions", _NESTED),)),
                tracks=(gbdraw.LinearTrackOptions if linear else gbdraw.CircularTrackOptions)(
                    slots=("features:features", slot)
                ),
            ),
        )._drawing
    svg = drawing.tostring()
    items = _groups(svg)
    paints = {name: (_band if linear else _signed_radii)(items[name]) for name in "abc"}
    bounds = _clip_bounds(svg, linear, size)
    left_out = 3 - len(drawn)

    assert all(_marked(items[name]) for name in "abc"), sorted(items)
    assert "".join(name for name in "abc" if _labelled(items[name])) == drawn
    if bounds is None:
        # Not clipped: every mark is whole within the set size.
        assert partly_clipped == ""
        assert max(end for _start, end in paints.values()) - min(start for start, _end in paints.values()) <= size
    else:
        low, high = bounds
        assert all(low <= start < high for start, _end in paints.values()), (paints, bounds)
        assert "".join(name for name in "abc" if paints[name][1] > high + 1e-6) == partly_clipped, (paints, bounds)
    counted = "region label" if left_out == 1 else "region labels"
    assert any(f"{left_out} {counted} left out in clipped slots" in str(item.message) for item in caught), [
        str(item.message) for item in caught
    ]
    report = label_fit_report(drawing)
    assert report is not None and report.drawn
    assert (
        "Region labels are left out: their clipped slots have no room for them, would cut their letters, "
        f"or hide their marks ({left_out})."
    ) in label_fit_notice(report)


def test_a_clipped_slot_leaves_out_a_label_its_edge_would_cut() -> None:
    # LBL-ROOM-CLIP-FIT-OMIT (Q12 X): three lanes pass the set size of 49 px, so the slot clips. Z's letters, outside
    # the outer lane, would pass the set size and be cut: Z's label is left out. X's letters cross no mark and end
    # within the size: X's label is drawn.
    annotations = tuple(
        ResolvedRegionAnnotation(
            id=identity, set_id="s", record_index=0, segments=(segment,), midpoint_bp=(segment[0] + segment[1]) / 2,
            span_bp=segment[1] - segment[0], label=label, mark="band", lane=lane, style=RegionAnnotationStyle(),
        )
        for identity, segment, label, lane in (("x", (100, 200), "X", 0), ("y", (600, 700), "", 1), ("z", (250, 350), "Z", 2))
    )
    track = _clipped(annotations, lambda *_: LabelExtent(-10.0, 10.0, 8.0, drawn=True), 49.0)
    by_id = {placement.annotation.id: placement.annotation for placement in track.placements}

    assert track.clipped
    assert (track.without_room, track.rooms_stop) == (((0, "z"),), "clip")
    assert track.draws_label(by_id["x"]) and not track.draws_label(by_id["z"])


def _clipped_linear(regions: tuple[RegionAnnotation, ...], height: int) -> tuple[str, list[str]]:
    """A Linear slot of ``height`` px above the features, with overflow clip, and the warnings drawing it gives."""

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        svg = _svg("linear", regions, f"regions:annotations@side=above,h={height}px,set_id=regions,overflow=clip")
    return svg, [str(item.message) for item in caught]


def _marked_region(identity: str, start: int, end: int, label: str, mark: str, size: float | None = None) -> RegionAnnotation:
    return RegionAnnotation(
        identity,
        CoordinateSpan(None, start, end),
        label=label,
        mark=mark,
        style=RegionAnnotationStyle(fill="#ffcc00", label_font_size=size),
    )


def test_a_clipped_slot_measures_the_cut_from_the_paint_of_a_line_mark() -> None:
    # LBL-ROOM-CLIP-FIT-OMIT: a line mark paints its stroke at the middle of its lane, and its label stands
    # label_offset beyond that paint (``linear_annotation_mark_reach``). Three nested line marks pass 30 px, so the
    # slot clips; Alpha's 4 pt letters above the line in lane 0 end within the slot, so they are drawn whole (as on
    # dev), not left out as if the line were painted at its lane's edge.
    regions = (
        _marked_region("a", 1001, 9000, "Alpha", "line", 4.0),
        _marked_region("b", 2001, 3000, "", "line"),
        _marked_region("c", 2201, 2800, "", "line"),
    )
    svg, messages = _clipped_linear(regions, 30)
    items = _groups(svg)

    assert ElementTree.fromstring(svg).find(f".//{_SVG}clipPath") is not None
    assert _labelled(items["a"]), messages
    _left, _right, top, _bottom = _ink(items["a"])
    assert top >= 0.0
    assert not any("left out" in message for message in messages), messages


def test_a_clipped_slot_reports_the_left_out_labels_at_their_own_size() -> None:
    # The warning names the size of the region labels it counts: Alpha's own 4 pt, though every one is left out
    # (its band's letters would pass the slot's top), not the feature labels' size.
    regions = (
        _marked_region("a", 1001, 9000, "Alpha", "band", 4.0),
        _marked_region("b", 2001, 3000, "", "band"),
        _marked_region("c", 2201, 2800, "", "band"),
    )
    svg, messages = _clipped_linear(regions, 30)

    assert not _labelled(_groups(svg)["a"])
    assert "r1: 1 region label does not fit at 4 pt (1 region label left out in clipped slots). The figure is drawn." in (
        messages
    ), messages


@pytest.mark.parametrize(("height", "drawn"), ((30, False), (40, True)))
def test_a_clipped_slot_leaves_out_a_label_whose_mark_the_clip_hides(height: int, drawn: bool) -> None:
    # LBL-ROOM-CLIP-FIT-OMIT (OV-2675, Q16 P): three nested bands pass the set size, so the slot clips. At 30 px the
    # clip hides all of C's band in lane 2 (its paint starts past 30 px), so C's label, which needs no room and whose
    # letters fit, is left out with the labels the clip cuts and counted in the same warning. At 40 px part of C's
    # band shows, and its label is drawn.
    regions = (
        _marked_region("a", 1001, 5000, "", "band"),
        _marked_region("b", 1501, 2500, "", "band"),
        _marked_region("c", 2001, 9000, "C", "band"),
    )
    svg, messages = _clipped_linear(regions, height)
    items = _groups(svg)
    top, _bottom = _band(items["c"])

    assert all(_marked(items[name]) for name in "abc")
    assert (top >= height) is not drawn
    assert _labelled(items["c"]) is drawn
    left_out = [message for message in messages if "1 region label left out in clipped slots" in message]
    assert bool(left_out) is not drawn, messages


def test_a_clipped_linear_slot_clips_only_across_its_lanes() -> None:
    # LBL-ROOM-CLIP-FIT-OMIT: the set size bounds a slot across its lanes, not along the bar. D's label at the
    # bar's start hangs past its left end (the figure grows to hold it); the clip of the slot, clipped by the
    # three nested lanes of a, b and c, keeps all of its letters, so no fragment of it is drawn.
    regions = (
        _marked_region("a", 6001, 9000, "", "band", 8.0),
        _marked_region("b", 6501, 8000, "", "band", 8.0),
        _marked_region("c", 7001, 7500, "", "band", 8.0),
        _marked_region("d", 1, 200, "Dlabel long long text", "band", 8.0),
    )
    with warnings.catch_warnings(record=True):
        warnings.simplefilter("always")
        svg = _svg("linear", regions, "regions:annotations@side=below,h=40px,set_id=regions,overflow=clip")
    left, right, _top, _bottom = _ink(_groups(svg)["d"])
    (rect,) = ElementTree.fromstring(svg).find(f".//{_SVG}clipPath").iter(f"{_SVG}rect")
    clip_left = float(rect.get("x", "nan"))
    clip_right = clip_left + float(rect.get("width", "nan"))

    assert left < 0.0
    assert clip_left <= left and right <= clip_right, (left, right, clip_left, clip_right)


def test_the_room_clears_a_wider_stroke_of_the_mark_it_would_cross() -> None:
    track = _packed(b_stroke=4.0)
    room = 10.0 + 2.0 * _OFFSET + 0.5 * (RegionAnnotationStyle().stroke_width + 4.0)

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
        labels=LabelSizing(reach=_half_stroke_reach, font_size=10.0, extent=lambda *_: LabelExtent(-10.0, 10.0, 10.0, drawn=True)),
        join_origin=False,
    )

    assert track.lane_rooms == (LaneRoom(record_index=0, row=0, gap=0, px=_ROOM),)
    assert track.lane_room_px == pytest.approx(_ROOM)
    assert track.required_extent_px == pytest.approx(_BASE + _ROOM)
    assert (track.room_before(0, 0, 1), track.room_before(1, 0, 1)) == (_ROOM, 0.0)


def _signed_radii(item: ElementTree.Element) -> tuple[float, float]:
    """The inner and outer radius of a Circular band, as its arcs give them (a negative one draws through the
    centre)."""

    (path,) = item.iter(f"{_SVG}path")
    radii = [float(value) for value in re.findall(r"A (-?[\d.]+),", path.get("d", ""))]
    return min(radii), max(radii)


def test_a_linear_label_gets_room_before_the_lane_further_on_whose_mark_it_would_cross() -> None:
    # LBL-ROOM-BEFORE-MET-LANE (OV-2661): c (lane 2) labels above at bp 4000. b (lane 1) ends at 2500, short of the
    # letters, which reach a's mark in lane 0. The room goes between a and b: b keeps its place beside c.
    regions = (_region("a", 1001, 5000, ""), _region("b", 1501, 2500, ""), _region("c", 2001, 6000, "Gamma"))
    items, bare = _groups(_svg("linear", regions)), _groups(_svg("linear", tuple(replace(r, label="") for r in regions)))
    _left, _right, top, bottom = _ink(items["c"])
    c_top, _c_bottom = _band(items["c"])
    moved = c_top - _band(bare["c"])[0]

    assert top >= _band(items["a"])[1] + _PAINT
    assert bottom == pytest.approx(c_top - _EDGE, abs=0.01)
    assert _band(items["a"]) == _band(bare["a"])
    assert moved > 0.0
    assert _band(items["b"]) == pytest.approx(tuple(y + moved for y in _band(bare["b"])), abs=1e-9)


def test_a_circular_label_gets_room_before_the_lane_further_on_whose_mark_it_would_cross() -> None:
    # LBL-ROOM-BEFORE-MET-LANE (OV-2661): c (lane 0) labels outward at bp 2500, where a horizontal label's width runs
    # across the lanes. b (lane 1) spans only 1901-2100 along the circle, clear of the letters, which reach a's mark
    # in lane 2. The room goes between b and a: b keeps its place beside c.
    regions = (
        _region("c", 1801, 3200, "Gamma region"),
        _region("b", 1901, 2100, ""),
        _region("a", 2001, 3000, ""),
    )
    items = _groups(_svg("circular", regions))
    a_inner, a_outer = _radii(items["a"])
    _b_inner, b_outer = _radii(items["b"])

    assert _radii(items["b"])[0] == pytest.approx(_radii(items["c"])[1])
    assert a_inner > b_outer
    assert not [
        (radius, bp)
        for radius, bp in _ink_polar(items["c"])
        if a_inner - _PAINT < radius < a_outer + _PAINT and 2001 <= bp <= 3000
    ]


def _nested_inside(
    count: int,
    labels: tuple[str, ...] | None = None,
    extra: tuple[RegionAnnotation, ...] = (),
    slot: str = "",
    orientation: str = "horizontal",
) -> tuple[dict[str, ElementTree.Element], list[warnings.WarningMessage], object]:
    """``count`` nested bands inside the axis around 3 o'clock, each labelled in ``orientation`` (``labels``, by
    default "Region N label text"), so that each label would cross the next lane's mark
    (``followups2/questions/qrender.py q7-<count>``), and the ``extra`` regions; ``slot`` adds settings to the
    slot."""

    texts = labels or tuple(f"Region {index + 1} label text" for index in range(count))
    regions = tuple(
        _region(f"r{index}", 1501 + 100 * index, 3500 - 100 * index, texts[index], orientation) for index in range(count)
    ) + extra
    options = gbdraw.CircularOptions(
        annotations=AnnotationOptions(sets=(AnnotationSet("regions", regions),)),
        tracks=gbdraw.CircularTrackOptions(slots=("features:features", f"regions:annotations@side=inside,set_id=regions{slot}")),
    )
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        drawing = gbdraw.draw_circular(_record(), options=options)._drawing
    return _groups(drawing.tostring()), caught, drawing


@pytest.mark.parametrize(
    ("count", "without_room", "slot"), ((4, 1, ""), (5, 2, ""), (5, 2, ",w=1000px,overflow=clip"))
)
def test_an_inside_slot_adds_no_room_past_the_centre_and_warns(count: int, without_room: int, slot: str) -> None:
    # LBL-ROOM-STOP-AT-CENTRE-OMIT (OV-2665, OV-2666, Q13 X): the rooms of four or five crossing labels would carry
    # the lanes past the centre. The rooms that fit from the slot's outer edge are kept; the others are not added,
    # every band keeps its width on its own side of the centre, and every mark is drawn. The labels without room are
    # not drawn, and the warning counts them. A clipped slot whose set width reaches past the centre stops there
    # too: the centre is the nearer bound.
    from gbdraw.api.diagram import label_fit_report
    from gbdraw.labels.fit import label_fit_notice
    from gbdraw.web_support.label_fit import label_fit_report_payload

    items, caught, drawing = _nested_inside(count, slot=slot)
    bands = [_signed_radii(items[f"r{index}"]) for index in range(count)]
    widths = [outer - inner for inner, outer in bands]

    assert all(inner > 0.0 for inner, _outer in bands), bands
    # One width (to the SVG's three decimals of each radius).
    assert max(widths) - min(widths) < 0.002 and min(widths) > 10.0, bands
    assert all(_marked(items[f"r{index}"]) for index in range(count))
    assert sum(_labelled(items[f"r{index}"]) for index in range(count)) == count - without_room
    counted = "region label" if without_room == 1 else "region labels"
    assert any(
        f"{without_room} {counted} left out for want of room before the centre" in str(item.message) for item in caught
    ), [str(item.message) for item in caught]
    report = label_fit_report(drawing)
    assert report is not None and report.drawn
    sentence = f"Region labels are left out: there is no room for them before the centre ({without_room})."
    assert sentence in label_fit_notice(report)
    # The Web shows the same sentence over the figure, and the region labels get no row to choose.
    payload = label_fit_report_payload(report, "circular")
    assert sentence in payload["notice"]
    (region,) = [record for record in payload["records"] if record["labelKind"] == "circular_region"]
    assert region["counts"] == f"{without_room} {counted} left out for want of room before the centre"
    assert not [row for row in report.rows if row.kind not in {"smaller_definition_font", "narrow_tracks"}]


def _explicit_lane_inside(
    second: tuple[int, int] = (2505, 2525),
) -> tuple[str, list[warnings.WarningMessage], object]:
    """Two short bands at 3 o'clock, the second over ``second`` (next to the first by default), both in explicit
    lane 0 of an inside slot with ``overflow=compress``, with long level labels."""

    regions = tuple(
        replace(_region(identity, start, end, label), lane=0)
        for identity, (start, end), label in zip(("a", "b"), ((2470, 2490), second), ("A long region label", "Another long region label"))
    )
    options = gbdraw.CircularOptions(
        annotations=AnnotationOptions(sets=(AnnotationSet("regions", regions),)),
        tracks=gbdraw.CircularTrackOptions(
            slots=("features:features", "regions:annotations@side=inside,set_id=regions,overflow=compress")
        ),
    )
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        drawing = gbdraw.draw_circular(_record(), options=options)._drawing
    return drawing.tostring(), caught, drawing


def test_inside_region_labels_over_each_other_are_drawn_and_warned() -> None:
    # OV-3008, OV-3015 (REGION-WARN, as in Linear): explicit lanes under overflow=compress keep two labels in one
    # lane although their letters meet, and the placer finds no clear place along their short regions. The figure
    # is drawn as before, and the warning counts the overlap; the Web shows the Linear sentence over the figure,
    # with no row to choose for the region labels.
    from gbdraw.api.diagram import label_fit_report
    from gbdraw.labels.fit import label_fit_notice
    from gbdraw.web_support.label_fit import label_fit_report_payload

    svg, caught, drawing = _explicit_lane_inside()

    assert all(_labelled(item) for item in _groups(svg).values())
    assert any("1 region label overlap" in str(item.message) for item in caught), [str(item.message) for item in caught]
    report = label_fit_report(drawing)
    assert report is not None and report.drawn and report.region_overlaps == 1
    sentence = "Region labels overlap other labels or lines (1)."
    assert sentence in label_fit_notice(report)
    payload = label_fit_report_payload(report, "circular")
    assert sentence in payload["notice"]
    assert not [row for row in report.rows if row.kind not in {"smaller_definition_font", "narrow_tracks"}]


def test_inside_region_labels_clear_of_each_other_add_no_warning() -> None:
    # OV-3015: the same labels with their regions apart keep clear of each other and add no region label overlap.
    from gbdraw.api.diagram import label_fit_report

    _svg_text, caught, drawing = _explicit_lane_inside((3000, 3020))

    assert not any("region label overlap" in str(item.message) for item in caught)
    report = label_fit_report(drawing)
    assert report is None or report.region_overlaps == 0


def test_an_inside_label_that_keeps_no_letter_is_left_out_and_counted() -> None:
    # RGN-CUT-EMPTY-OMIT (OV-2804, PD-OI-097 scenario revision 22): ten nested bands with labels along the circle.
    # The innermost band's label would be cut so far that no letter is left, only the ellipsis: it is not drawn,
    # and the warning counts it as a label left out, not as a cut one. Its mark is drawn.
    from gbdraw.api.diagram import label_fit_report
    from gbdraw.labels.fit import label_fit_notice
    from gbdraw.web_support.label_fit import label_fit_report_payload

    items, caught, drawing = _nested_inside(10, orientation="auto")

    assert all(_marked(items[f"r{index}"]) for index in range(10))
    assert not _labelled(items["r9"])
    drawn = ["".join(text.itertext()) for item in items.values() for text in item.iter(f"{_SVG}text")]
    assert drawn == [f"Region {index + 1} label text" for index in range(7)], drawn
    counts = "2 region labels left out for want of room before the centre, 1 region label left out: no letter fits its circle"
    assert any(counts in str(item.message) for item in caught), [str(item.message) for item in caught]
    report = label_fit_report(drawing)
    assert report is not None and report.drawn
    assert (report.region_cuts, report.region_cuts_left_out, report.region_rooms_past_centre) == (0, 1, 2)
    sentence = "Region labels are left out: no letter of them fits their circle (1)."
    assert sentence in label_fit_notice(report)
    payload = label_fit_report_payload(report, "circular")
    assert sentence in payload["notice"]
    (region,) = [record for record in payload["records"] if record["labelKind"] == "circular_region"]
    # Every region label of the record, drawn or left out, is counted.
    assert region["counts"] == counts and region["labelCount"] == 10


class _LeftOut(NamedTuple):
    bands: dict[str, tuple[float, float]]
    without_letters: frozenset[str]
    messages: list[str]
    items: dict[str, ElementTree.Element]
    # The labels drawn where the drawer last placed the slot's labels among their obstacles.
    drawn: tuple[CircularAnnotationLabel, ...]


def _bands_and_left_out(
    regions: tuple[tuple[str, int, int, str], ...],
    slot: str,
    orientation: str,
    empty: frozenset[str] = frozenset(),
    size: float = 60,
    sizes: dict[str, float] | None = None,
) -> _LeftOut:
    """The band radii of ``regions`` (id, start, end, label) labelled at ``size`` pt (or their size in ``sizes``) in
    ``orientation`` in an inside slot (``slot`` adds settings), the labels of ``empty`` without text, the labels left
    out because no letter of them
    fits their circle (``CircularAnnotationLabels.without_letters``), as the assembler gets them, the warnings, each
    annotation's group and the labels drawn."""

    import gbdraw.diagrams.circular.assemble as assemble_module

    left_out: list[str] = []
    drawn: list[tuple[CircularAnnotationLabel, ...]] = []
    place = assemble_module.circular_annotation_labels

    def spying(track: object, **kwargs: object) -> object:
        labels = place(track, **kwargs)
        left_out.extend(label.annotation_id for label in labels.without_letters)
        if kwargs.get("obstacles") is not None:
            drawn.append(labels.drawn)
        return labels

    annotations = tuple(
        RegionAnnotation(
            identity,
            CoordinateSpan(None, start, end),
            label="" if identity in empty else label,
            mark="band",
            style=RegionAnnotationStyle(
                fill="#ffcc00", label_orientation=orientation, label_font_size=(sizes or {}).get(identity, size)
            ),
        )
        for identity, start, end, label in regions
    )
    options = gbdraw.CircularOptions(
        annotations=AnnotationOptions(sets=(AnnotationSet("regions", annotations),)),
        tracks=gbdraw.CircularTrackOptions(slots=("features:features", f"regions:annotations@side=inside,set_id=regions{slot}")),
    )
    with pytest.MonkeyPatch.context() as patch, warnings.catch_warnings(record=True) as caught:
        patch.setattr(assemble_module, "circular_annotation_labels", spying)
        warnings.simplefilter("always")
        items = _groups(gbdraw.draw_circular(_record(), options=options)._drawing.tostring())
    bands = {identity: _signed_radii(items[identity]) for identity, *_rest in regions}
    return _LeftOut(bands, frozenset(left_out), [str(item.message) for item in caught], items, drawn[-1] if drawn else ())


_SEVEN_NESTED = tuple((f"r{index}", 1501 + 100 * index, 3500 - 100 * index, f"Region {index + 1} label text") for index in range(7))
# OV-2805 (``ov2804/m1/case.tsv``, with ``w=600px,overflow=error``): r1 meets r2's letters in r2's lane, and the lane
# further in, where its letters meet nothing, lies past the label rooms so near the centre that no letter of it fits
# there (the set width stretches the lanes toward the centre).
_DODGING_INWARD = (
    ("r0", 418, 2243, "Region 0 xxxxxxxxxx"),
    ("r1", 8990, 10000, "Region 1 xxxxxxxxxxxx"),
    ("r2", 2125, 3625, "Region 2 xxxxxxxxxxxxxxxxxxxx"),
    ("r3", 2184, 5656, "Region 3 xxx"),
)


@pytest.mark.parametrize(
    ("regions", "orientation", "slot"),
    (
        (_SEVEN_NESTED, "arc", ""),
        (_SEVEN_NESTED, "arc", ",w=600px,overflow=error"),
        (_DODGING_INWARD, "auto", ",w=600px,overflow=error"),
    ),
)
def test_a_label_that_keeps_no_letter_moves_no_band(
    regions: tuple[tuple[str, int, int, str], ...], orientation: str, slot: str
) -> None:
    # RGN-CUT-EMPTY-OMIT (OV-2804, OV-2805): a label that keeps no letter on its circle is not drawn, so it claims no
    # room between the lanes, is no obstacle there, and does not move its own band inward to keep its letters clear
    # of other labels: every band stands where it stands when the labels left out have no text.
    drawn = _bands_and_left_out(regions, slot, orientation)

    assert drawn.without_letters, drawn.bands
    assert _bands_and_left_out(regions, slot, orientation, drawn.without_letters).bands == drawn.bands


# Reviewer probe figures (``random.Random(7)``, 48 pt, inside slot without a set width): figure 29, where the lane
# packer put r0 and r1 in one lane, measured where that lane stands before the rooms between the lanes are known
# (OV-2849), and figure 72, with ``w=600px,overflow=error``, where r4 was left out in an early packing round and stayed
# out although the last packing gives it letters (OV-2850).
_SHARED_LANE = (
    ("r0", 4802, 5087, "Region 0 xxxxxxxxxxxxxx"), ("r1", 3037, 3782, "Region 1 xxxxxxxx"),
    ("r2", 7305, 7419, "Region 2 xxxxxxxx"), ("r3", 5967, 7414, "Region 3 xxxxxxxxxxxxxxxxx"),
    ("r4", 5301, 6402, "Region 4 x"), ("r5", 5072, 6064, "Region 5 xxxxxxxxxxx"),
    ("r6", 2998, 3102, "Region 6 xxxxxxxxxx"), ("r7", 6253, 6696, "Region 7 xxxxxxxxxxxxxxx"),
)
_LEFT_OUT_EARLY = (
    ("r0", 5750, 8350, "Region 0 xxxxxxxx"), ("r1", 2602, 4028, "Region 1 xxxxxxxxxxxxxxxxxxx"),
    ("r2", 4506, 8313, "Region 2 xxxxxxxxxxxxxx"), ("r3", 2353, 3494, "Region 3 xxxxxxxxxxxxxxxx"),
    ("r4", 7867, 8820, "Region 4 xxxxxxxxxxxxxxxxxx"), ("r5", 4307, 6929, "Region 5 xxxxxxxxxxxxxxxx"),
    ("r6", 3890, 5296, "Region 6 xxxxxxxxxxx"), ("r7", 604, 1518, "Region 7 xxxxx"),
    ("r8", 6611, 7371, "Region 8 xxxxxxxxxxxxxxxxxxxx"),
)


@pytest.mark.parametrize(
    ("regions", "slot", "size"),
    ((_DODGING_INWARD, ",w=600px,overflow=error", 60), (_SHARED_LANE, "", 48)),
)
def test_labels_that_share_a_lane_stay_apart_where_the_rooms_carry_that_lane(
    regions: tuple[tuple[str, int, int, str], ...], slot: str, size: float
) -> None:
    # OV-2849: the rooms between lanes carry an inside slot's inner lanes toward the centre, and a set width shares
    # its extent among the lanes; there a label's letters span more of the circle than where the lane packer first
    # measured them. The packer measures each lane's letters where the drawer draws that lane, so the labels it puts
    # in one lane do not meet.
    from gbdraw.layout.text_geometry import convex_polygons_intersect

    drawn = _bands_and_left_out(regions, slot, "auto", size=size).drawn
    met = [
        (first.annotation_id, second.annotation_id)
        for index, first in enumerate(drawn)
        for second in drawn[index + 1 :]
        if any(convex_polygons_intersect(one, other) for one in first.parts for other in second.parts)
    ]

    assert len(drawn) > 1 and not met, met


def test_a_label_drawn_where_a_set_width_squeezes_the_lanes_gets_its_room() -> None:
    # OV-2989: nine nested bands at 24 pt along their circles, 10 px apart, in an inside slot whose set width reaches
    # past the centre (``w=600px,overflow=compress,lane_gap_px=10``). The lanes and their rooms take more than the
    # centre leaves, so the drawer narrows the gaps and shares the band up to the centre less the rooms among the
    # lanes, narrower than the extent they take. (The figure first reported, ten bands with ``overflow=error``, is
    # refused there now, OV-2990.) r5's label claims no room where the lanes stand at that extent, but is drawn where
    # the drawer draws its lane: it gets its room there, as every label drawn does, and no drawn label's letters come
    # within its label offset of another lane's painted mark (measured only at the extent, r5's came within it of r6's).
    regions = tuple(
        (f"r{index}", 1501 + 90 * index, 3500 - 90 * index, f"Region {index + 1} label text") for index in range(9)
    )
    drawn = _bands_and_left_out(regions, ",w=600px,overflow=compress,lane_gap_px=10", "arc", size=24)

    assert "r5" in {label.annotation_id for label in drawn.drawn}, drawn.messages
    met = [
        (label.annotation_id, identity, label.ink_radii, band)
        for label in drawn.drawn
        for identity, band in drawn.bands.items()
        if identity != label.annotation_id
        # The band radii are written to 0.001 px.
        and label.ink_radii[0] < band[1] + _EDGE - 0.002
        and band[0] - _EDGE + 0.002 < label.ink_radii[1]
    ]
    assert not met, met


def test_a_label_two_lanes_from_the_mark_it_crosses_gets_its_room_where_a_set_width_squeezes_the_lanes() -> None:
    # The set width squeezes the lanes as above (seven nested bands 20 px apart, ``overflow=compress``), so each lane
    # between a label and the mark it crosses gives it less than the extent the lanes take. c's 150 pt label stands in
    # lane 0 beside its mark and crosses the marks of b (lane 1) and a (lane 2): its room counts the lanes between
    # where the drawer draws them, and its letters stay a label offset clear of a's painted mark (counted at the
    # extent, they came 7.5 px into it).
    regions = (
        *((f"r{index}", 1501 + 90 * index, 3500 - 90 * index, f"Region {index + 1} label text") for index in range(7)),
        ("c", 5801, 7200, "G"),
        ("b", 5900, 6100, ""),
        ("a", 6050, 7600, ""),
    )
    drawn = _bands_and_left_out(
        regions, ",w=600px,overflow=compress,lane_gap_px=20", "arc", frozenset({"a", "b"}), 24, {"c": 150}
    )

    (label,) = [label for label in drawn.drawn if label.annotation_id == "c"]
    # The band radii are written to 0.001 px.
    assert label.ink_radii[0] >= drawn.bands["a"][1] + _EDGE - 0.002, (label.ink_radii, drawn.bands["a"])


def test_a_label_left_out_in_an_early_packing_is_drawn_where_the_last_packing_gives_it_letters() -> None:
    # OV-2850: r4 keeps no letter in the lane an early packing round puts it in, so that round packs it without text;
    # the last packing puts it in the outer lane, where its letters fit, and it is drawn there, as before.
    drawn = _bands_and_left_out(_LEFT_OUT_EARLY, ",w=600px,overflow=error", "auto", size=48)

    assert "r4" in {label.annotation_id for label in drawn.drawn}, drawn.messages


def test_an_inside_slot_stands_nowhere_that_leaves_out_a_label_its_lanes_draw() -> None:
    # OV-2848: four nested bands at 56 pt in an inside slot without a set width. Further in, where the letters of the
    # other labels find more room, the innermost label keeps no letter; the slot does not move there, so every label
    # is drawn, and the lanes settle before the assembler's last pass. At 60 pt the innermost label keeps no letter
    # wherever the slot stands: its corners come within its label offset of its mark (OV-2979). It alone is left out.
    import gbdraw.diagrams.circular.assemble as assemble_module

    # Whether each packing where the layout draws the slots changed it, for each layout of the drawing.
    layouts: list[list[bool]] = []
    place, repack = assemble_module._place_annotation_tracks_where_drawn, assemble_module._annotation_tracks_packed_where_drawn

    def placing(*args: object, **kwargs: object) -> object:
        layouts.append([])
        return place(*args, **kwargs)

    def repacking(*args: object, **kwargs: object) -> dict[str, object]:
        result = repack(*args, **kwargs)
        layouts[-1].append(bool(result))
        return result

    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(assemble_module, "_place_annotation_tracks_where_drawn", placing)
        patch.setattr(assemble_module, "_annotation_tracks_packed_where_drawn", repacking)
        drawn = _bands_and_left_out(_SEVEN_NESTED[:4], "", "auto", size=56)

    assert not drawn.without_letters, drawn.messages
    assert all(_labelled(drawn.items[identity]) for identity, *_rest in _SEVEN_NESTED[:4])
    assert layouts and all(changed[-1:] == [False] for changed in layouts), layouts

    larger = _bands_and_left_out(_SEVEN_NESTED[:4], "", "auto")
    assert larger.without_letters == {"r3"}, larger.messages
    assert all(_labelled(larger.items[identity]) for identity, *_rest in _SEVEN_NESTED[:3])


def test_a_label_left_out_for_want_of_a_letter_keeps_its_band_where_it_stands_without_text() -> None:
    # OV-2805: r1's letters meet r2's in r2's lane, and the lane further in lies so near the centre that no letter
    # of r1 fits there. r1's band stays in r2's lane, beside r2's mark, not at the centre.
    drawn = _bands_and_left_out(_DODGING_INWARD, ",w=600px,overflow=error", "auto")

    assert "r1" in drawn.without_letters
    assert drawn.bands["r1"] == pytest.approx(drawn.bands["r2"], abs=0.01), drawn.bands


def test_a_label_without_letters_in_a_clipped_slot_is_counted_once_as_without_letters() -> None:
    # RGN-CUT-EMPTY-OMIT with LBL-ROOM-CLIP-FIT-OMIT: seven nested bands at 60 pt in a clipped inside slot 300 px
    # wide. The clip leaves out four labels; the innermost keeps no letter on its circle and is counted once, as a
    # label without letters, not also as cut by the clip.
    left_out = _bands_and_left_out(_SEVEN_NESTED, ",w=300px,overflow=clip", "arc")
    not_drawn = sorted(identity for identity, item in left_out.items.items() if not _labelled(item))

    assert left_out.without_letters == {"r6"}
    assert not_drawn == ["r2", "r3", "r4", "r5", "r6"]
    counts = "4 region labels left out in clipped slots, 1 region label left out: no letter fits its circle"
    assert any(counts in message for message in left_out.messages), left_out.messages


def test_the_warning_reports_a_label_without_letters_at_its_own_size() -> None:
    # The label fit warning names the largest font size among the labels it counts: the innermost label, at 20 pt,
    # keeps no letter on its circle, and the warning reports 20 pt, not the 14 pt of the labels drawn.
    from gbdraw.api.diagram import label_fit_report
    from gbdraw.web_support.label_fit import label_fit_report_payload

    innermost = RegionAnnotation(
        "z", CoordinateSpan(None, 2401, 2600), label="Region z label text", mark="band",
        style=RegionAnnotationStyle(fill="#ffcc00", label_orientation="auto", label_font_size=20),
    )
    items, caught, drawing = _nested_inside(9, orientation="auto", extra=(innermost,))

    assert not _labelled(items["z"])
    assert any("10 region labels at 20 pt" in str(item.message) for item in caught), [str(item.message) for item in caught]
    payload = label_fit_report_payload(label_fit_report(drawing), "circular")
    (region,) = [record for record in payload["records"] if record["labelKind"] == "circular_region"]
    assert region["fontSize"] == 20.0 and "1 region label left out: no letter fits its circle" in region["counts"]


def test_a_label_without_letters_puts_no_ink_on_the_definition() -> None:
    # A label left out because no letter of it fits its circle paints nothing: the definition line its lone
    # ellipsis met at the centre no longer overlaps a track (``ov2804/scripts/probe.py``, random set 2).
    regions = (
        ("q0", 5916, 6116, "Region label", "bracket"), ("q1", 3477, 3507, "origin of replication", "band"),
        ("q2", 7057, 10000, "origin of replication", "bracket"), ("q3", 8916, 10000, "origin of replication", "bracket"),
        ("q4", 589, 619, "Region label", "line"), ("q5", 5218, 8218, "A long region label text that runs on", "band"),
        ("q6", 9184, 9384, "IRa", "band"), ("q7", 391, 591, "Region label", "band"),
        ("q8", 2240, 3040, "origin of replication", "band"), ("q9", 7302, 10000, "origin of replication", "bracket"),
        ("q10", 5797, 6597, "A long region label text that runs on", "band"),
        ("q11", 6552, 9552, "origin of replication", "band"),
        ("q12", 8029, 8829, "A long region label text that runs on", "bracket"), ("q13", 7450, 10000, "Region label", "line"),
    )
    annotations = tuple(
        RegionAnnotation(
            identity, CoordinateSpan(None, start, end), label=label, mark=mark,
            style=RegionAnnotationStyle(fill="#ffcc00", label_orientation="auto"),
        )
        for identity, start, end, label, mark in regions
    )
    options = gbdraw.CircularOptions(
        annotations=AnnotationOptions(sets=(AnnotationSet("regions", annotations),)),
        tracks=gbdraw.CircularTrackOptions(slots=("features:features", "regions:annotations@side=inside,set_id=regions,r=0.3")),
    )
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        gbdraw.draw_circular(_record(), options=options)
    messages = [str(item.message) for item in caught]

    assert any("1 region label left out: no letter fits its circle" in message for message in messages), messages
    assert not any("overlap" in message and "definition line" in message for message in messages), messages


def test_an_inside_slot_adds_each_room_that_fits_before_the_centre() -> None:
    # LBL-ROOM-STOP-AT-CENTRE-OMIT: a room that would carry the lanes past the centre is not added; each room is
    # judged on its own. Five nested bands: the rooms of r0's and r1's (short) labels fit, r2's (about 68 px) would
    # pass the centre, and r3's short label "x" needs about 15 px, which still fits after it: r4 moves on by it.
    labels = ("Region 1 label text", "Region 2", "Region 3 label text", "x", "")
    items, caught, _drawing = _nested_inside(5, labels)
    bands = {name: _signed_radii(items[name]) for name in ("r2", "r3", "r4")}

    # The lanes touch where no room parts them (each lane takes its share of the lane gaps).
    assert bands["r2"][0] - bands["r3"][1] == pytest.approx(0.0, abs=0.01), bands
    assert bands["r3"][0] - bands["r4"][1] > 10.0, bands
    assert any("1 region label left out for want of room before the centre" in str(item.message) for item in caught), [
        str(item.message) for item in caught
    ]


@pytest.mark.parametrize(
    ("mode", "regions", "shares_lane"),
    (
        # The letters (bp 1711-2390) reach over b's mark; 0.6 em would keep only bp 1834-2266.
        ("linear", (("a", 2001, 2100, "WWWWWW", "horizontal"), ("b", 2301, 2600, "", "horizontal")), False),
        # The letters (bp 1979-2122) end short of b; 0.6 em would keep up to bp 2266.
        ("linear", (("a", 2001, 2100, "iiiiii", "horizontal"), ("b", 2201, 2500, "", "horizontal")), True),
        ("circular", (("a", 2001, 2050, "WWWWWW", "tangent"), ("b", 2151, 2400, "", "tangent")), False),
        # The letters end short of b at bp 2101, where 0.6 em kept b in the next lane.
        ("circular", (("a", 2001, 2050, "iiiiii", "tangent"), ("b", 2101, 2350, "", "tangent")), True),
    ),
    ids=("linear-wide", "linear-narrow", "circular-wide", "circular-narrow"),
)
def test_the_lane_packer_keeps_a_labels_letters_as_the_drawer_measures_them(
    mode: str, regions: tuple[tuple[str, int, int, str, str], ...], shares_lane: bool
) -> None:
    # The lane packer keeps a's letters as the drawer measures them: a's label leaves its short mark, and b shares
    # a's lane only where the letters end short of b's mark.
    items = _groups(_svg(mode, tuple(_region(*region) for region in regions)))
    shape = _band if mode == "linear" else _radii

    assert (shape(items["a"]) == shape(items["b"])) is shares_lane


def test_a_slot_whose_lanes_stay_where_it_is_drawn_lays_the_tracks_out_once(monkeypatch: pytest.MonkeyPatch) -> None:
    # The letters measured where an outside slot is drawn span fewer bp than at the axis circle, but a's lane stays:
    # the packing there equals the one at the axis circle, so the tracks are laid out once. (A packing that kept the
    # bp its letters span at the radius they were measured at would differ, and lay the tracks out twice.)
    from gbdraw.diagrams.circular import assemble

    original = assemble._place_annotation_tracks_where_drawn
    passes: list[int] = []

    def counting(place, slots, tracks):
        calls: list[object] = []

        def counted(*args):
            calls.append(args)
            return place(*args)

        result = original(counted, slots, tracks)
        passes.append(len(calls))
        return result

    monkeypatch.setattr(assemble, "_place_annotation_tracks_where_drawn", counting)
    regions = (_region("a", 2001, 2050, "WWWWWW", "tangent"), _region("b", 6001, 6400, ""))
    _svg("circular", regions, "regions:annotations@side=outside,set_id=regions")

    assert passes and set(passes) == {1}, passes


def _two_stacks(gap_bp: int) -> tuple[ResolvedRegionAnnotation, ...]:
    """Two stacks of three nested labelled regions (a0-a2 from bp 100, b0-b2 ``gap_bp`` further on): each stack takes
    three lanes."""

    return tuple(
        ResolvedRegionAnnotation(
            id=f"{stack}{index}", set_id="s", record_index=0, segments=((start + index, start + 100),),
            midpoint_bp=start + 50.0, span_bp=100 - index, label=f"{stack}{index}", mark="band", lane=None,
            style=RegionAnnotationStyle(),
        )
        for stack, start in (("a", 100), ("b", 100 + gap_bp))
        for index in range(3)
    )


def _inside_lanes(annotations, letters, limit: float) -> dict[str, int]:
    track = layout_annotation_track(
        "regions", "s", annotations, record_lengths={0: 10_000}, params=AnnotationTrackParams(set_id="s"),
        labels=LabelSizing(reach=_half_stroke_reach, font_size=10.0, extent=letters), join_origin=True, room_bounds=packer.RoomBounds(limit_px=limit),
    )
    return {item.annotation.id: item.lane for item in track.placements}


def test_an_inside_slot_measures_each_lane_where_it_is_drawn() -> None:
    # The lane packer measures letters lane by lane: the third lane lies further in, where they span more of the
    # circle (400 bp instead of 50). Measured there, a2's letters reach b2's mark 300 bp on, so b2 takes a lane of
    # its own.
    def letters(_annotation, _style, _font_size, _place_bp, mark_offset):
        reach = 400.0 if mark_offset > 40.0 else 50.0
        return LabelExtent(-reach, reach, 10.0, drawn=True)

    lanes = _inside_lanes(_two_stacks(300), letters, limit=1_000.0)

    assert lanes == {"a0": 0, "a1": 1, "a2": 2, "b0": 0, "b1": 1, "b2": 3}


def test_an_inside_slot_whose_lanes_would_pass_the_centre_measures_its_letters_in_its_outer_lane() -> None:
    # The slot has 40 px before the centre, room for two lanes of 14 px. A label in a lane past that, measured
    # where the band would put it, covers the whole circle, so each would take a lane of its own and the lanes would
    # pass the centre. The slot measures every lane beside a mark in its outer lane instead: the two stacks share
    # their three lanes.
    def letters(_annotation, _style, _font_size, _place_bp, mark_offset):
        reach = 50.0 if mark_offset < 20.0 else 4_900.0
        return LabelExtent(-reach, reach, 10.0, drawn=True)

    lanes = _inside_lanes(_two_stacks(5_000), letters, limit=40.0)

    assert lanes == {"a0": 0, "a1": 1, "a2": 2, "b0": 0, "b1": 1, "b2": 2}


def test_a_clipped_inside_slot_stops_its_lanes_at_the_centre() -> None:
    # OV-2666: an inside slot of set width with overflow=clip, 25 px from the centre, whose five lanes pass the
    # centre. They stop at it: no band is drawn at a negative radius, through the centre.
    regions = tuple(
        _region(f"r{index}", 1501 + 100 * index, 3500 - 100 * index, f"Region {index + 1} label text")
        for index in range(5)
    )
    items = _groups(
        _svg("circular", regions, "regions:annotations@side=inside,set_id=regions,r=25px,w=60px,overflow=clip")
    )
    bands = [_signed_radii(items[f"r{index}"]) for index in range(5)]

    assert all(inner >= 0.0 for inner, _outer in bands), bands
    assert bands[-1] == (0.0, 0.0)
    ((path,),) = [list(items["r4"].iter(f"{_SVG}path"))]
    assert "nan" not in path.get("d", "").lower()


def test_an_inside_slot_keeps_the_room_of_a_label_that_fits_where_another_of_its_gap_does_not() -> None:
    # LBL-ROOM-STOP-AT-CENTRE-OMIT, each room on its own: two lane-0 labels would cross lane 1's marks. A needs
    # 69.5 px, B 17.5 px; the lanes take 35 px of the 55 px before the centre. The gap takes B's room, and only A is
    # left without room.
    annotations = tuple(
        ResolvedRegionAnnotation(
            id=identity, set_id="s", record_index=0, segments=(segment,), midpoint_bp=(segment[0] + segment[1]) / 2,
            span_bp=segment[1] - segment[0], label=label, mark="band", lane=lane, style=RegionAnnotationStyle(),
        )
        for identity, segment, label, lane in (
            ("A", (100, 200), "LONG", 0), ("B", (600, 700), "x", 0), ("C", (100, 200), "", 1), ("D", (600, 700), "", 1)
        )
    )

    def letters(annotation, _style, _font_size, _place_bp, _mark_offset):
        return LabelExtent(-10.0, 10.0, 60.0 if annotation.label == "LONG" else 8.0, drawn=True)

    track = layout_annotation_track(
        "regions", "s", annotations, record_lengths={0: 1000}, params=AnnotationTrackParams(set_id="s"),
        labels=LabelSizing(reach=_half_stroke_reach, font_size=10.0, extent=letters), join_origin=True,
        room_bounds=packer.RoomBounds(limit_px=_BASE + 20.0),
    )

    assert track.lane_rooms == (LaneRoom(0, 0, 0, pytest.approx(8.0 + 2.0 * _OFFSET + RegionAnnotationStyle().stroke_width)),)
    assert track.labels_without_room(0) == 1


def test_an_inside_slot_keeps_a_smaller_room_of_a_gap_whose_largest_room_does_not_fit() -> None:
    # Five nested bands: the rooms of r0's and r1's (short) labels fit, r2's would pass the centre. q2 shares r2's
    # lane elsewhere on the circle, labelled "x", over q3's mark in r3's lane: the gap keeps q2's short room, and only
    # r2's label is left without room.
    labels = ("Region 1 label text", "Region 2", "Region 3 label text", "", "")
    extra = (replace(_region("q2", 6001, 6600, "x"), lane=2), replace(_region("q3", 6001, 6600, ""), lane=3))
    items, caught, _drawing = _nested_inside(5, labels, extra)
    q2, q3 = _signed_radii(items["q2"]), _signed_radii(items["q3"])

    assert q2[0] - q3[1] > 10.0, (q2, q3)
    assert any("1 region label left out for want of room before the centre" in str(item.message) for item in caught), [
        str(item.message) for item in caught
    ]


def test_a_linear_label_is_measured_with_the_lanes_of_its_own_record() -> None:
    # LBL-ROOM-BEFORE-MET-LANE: in one Linear slot, record 0 takes six lanes and record 1 three, so each of record
    # 1's lanes is drawn twice as high (the drawers share the slot among a record's own lanes). c's letters (30 px
    # deep) pass b's short mark in lane 1 and meet a's in lane 2. The room before lane 2 counts lane 1 at the height
    # it is drawn, so the letters keep the label offset before a's mark.
    from gbdraw.annotations.layout import lane_edge_px, lane_width_px

    def annotation(identity: str, record: int, segment: tuple[int, int], label: str, lane: int) -> ResolvedRegionAnnotation:
        return ResolvedRegionAnnotation(
            id=identity, set_id="s", record_index=record, segments=(segment,), midpoint_bp=(segment[0] + segment[1]) / 2,
            span_bp=segment[1] - segment[0], label=label, mark="band", lane=lane, style=RegionAnnotationStyle(),
        )

    annotations = tuple(annotation(f"z{lane}", 0, (100, 900), "", lane) for lane in range(6)) + (
        annotation("c", 1, (100, 900), "Gamma", 0), annotation("b", 1, (100, 150), "", 1), annotation("a", 1, (400, 600), "", 2)
    )
    track = layout_annotation_track(
        "regions", "s", annotations, record_lengths={0: 1000, 1: 1000}, params=AnnotationTrackParams(set_id="s"),
        labels=LabelSizing(reach=_half_stroke_reach, font_size=10.0, extent=lambda *_: LabelExtent(-50.0, 50.0, 30.0, drawn=True)), join_origin=False,
    )
    padding, half_stroke = AnnotationTrackParams(set_id="s").padding_px, 0.5 * RegionAnnotationStyle().stroke_width
    height = lane_width_px(track.required_extent_px, padding_px=padding, room_px=track.lane_room_px, lane_count=3)
    letters_end = lane_edge_px(1, height, padding_px=padding, room_before_px=0.0) + half_stroke + _OFFSET + 30.0
    a_paint = lane_edge_px(2, height, padding_px=padding, room_before_px=track.room_before(1, 0, 2)) - half_stroke

    assert [(room.record_index, room.gap) for room in track.lane_rooms] == [(1, 1)]
    assert a_paint - letters_end == pytest.approx(_OFFSET, abs=1e-9)
