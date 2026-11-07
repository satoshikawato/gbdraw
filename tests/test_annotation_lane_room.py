"""Room between two annotation lanes for a label that would cross the next lane's mark (PD-OI-119, OV-2575).

Where a lane's label, at its place along its region, would cross a mark in the next lane on the label's side, the
two lanes are drawn apart by the room the label takes (its letters and its ``label_offset``); the label stays
beside its own mark in that room. Lanes whose labels cross nothing keep their spacing. The track grows by the room;
a slot whose size the user fixed adds it to its required extent, and its ``overflow`` decides (PD-OI-119).
"""

from __future__ import annotations

import math
import re
from xml.etree import ElementTree

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import gbdraw
from gbdraw.annotations import AnnotationTrackParams, LabelSizing, LaneRoom, ResolvedRegionAnnotation, layout_annotation_track
from gbdraw.api import AnnotationOptions, AnnotationSet, CoordinateSpan, RegionAnnotation, RegionAnnotationStyle
from gbdraw.exceptions import ValidationError
from gbdraw.render.drawers.annotation_marks import label_ink_offsets

_SVG = "{http://www.w3.org/2000/svg}"
_OFFSET = RegionAnnotationStyle().label_offset
# The letters keep the label offset from the painted edge of their mark, half its stroke past the band.
_EDGE = _OFFSET + RegionAnnotationStyle().stroke_width / 2


def _items(mode: str, b_span: tuple[int, int], slot: str | None = None, *, labels: bool = True) -> dict[str, ElementTree.Element]:
    """Two bands on a 10 kb record: ``a`` over 1000-5000 (lane 0) and ``b`` over ``b_span`` (lane 1 where they meet)."""

    record = SeqRecord(
        Seq("ACGT" * 2500), id="r1", name="r1", annotations={"topology": "circular", "molecule_type": "DNA"}
    )
    regions = tuple(
        RegionAnnotation(
            identity,
            CoordinateSpan(None, start, end),
            label=label if labels else "",
            mark="band",
            style=RegionAnnotationStyle(fill="#ffcc00", label_orientation="horizontal"),
        )
        for identity, (start, end), label in (("a", (1000, 5000), "Alpha"), ("b", b_span, "Beta"))
    )
    linear = mode == "linear"
    tracks = {}
    if slot is not None:
        slots = ("features:features", slot)
        tracks = {"tracks": gbdraw.LinearTrackOptions(slots=slots) if linear else gbdraw.CircularTrackOptions(slots=slots)}
    options = (gbdraw.LinearOptions if linear else gbdraw.CircularOptions)(
        annotations=AnnotationOptions(sets=(AnnotationSet("regions", regions),)), **tracks
    )
    diagram = (gbdraw.draw_linear if linear else gbdraw.draw_circular)(record, options=options)
    root = ElementTree.fromstring(diagram._drawing.tostring())
    return {
        group.get("data-gbdraw-annotation-id", ""): group
        for group in root.iter(f"{_SVG}g")
        if group.get("data-gbdraw-annotation-id") in {"a", "b"}
    }


def _band(item: ElementTree.Element) -> tuple[float, float]:
    """The top and bottom y of a Linear band."""

    (rect,) = item.iter(f"{_SVG}rect")
    top = float(rect.get("y", "nan"))
    return top, top + float(rect.get("height", "nan"))


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


def _ink_radii(item: ElementTree.Element) -> tuple[float, float]:
    """How near the circle's centre and how far from it a horizontal label's letters come."""

    left, right, top, bottom = _ink(item)
    near = math.hypot(max(left, 0.0, -right), max(top, 0.0, -bottom))
    far = max(math.hypot(x, y) for x in (left, right) for y in (top, bottom))
    return near, far


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


def test_lanes_whose_labels_cross_nothing_keep_their_spacing() -> None:
    # b (4000-9000) shares 4000-5000 with a, so it takes lane 1; its label at bp 6500 is past a's mark.
    labelled = _items("linear", (4000, 9000))
    bare = _items("linear", (4000, 9000), labels=False)

    assert [_band(labelled[key]) for key in "ab"] == [_band(bare[key]) for key in "ab"]


def _packed(
    *, overflow: str = "error", available: float | None = None, above: bool = False, join_origin: bool = False, **params: object
):
    """``a`` over 100-500 labelled ``A`` and ``b`` over 200-300, with letters 10 px deep."""

    annotations = tuple(
        ResolvedRegionAnnotation(
            id=identity, set_id="s", record_index=0, segments=(segment,), midpoint_bp=(segment[0] + segment[1]) / 2,
            span_bp=segment[1] - segment[0], label=label, mark="band", lane=None, style=RegionAnnotationStyle(),
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
        labels=LabelSizing(font_size=10.0, depth_px=lambda _text, _style, _size, _bp: 10.0, labels_above=above),
        join_origin=join_origin,
    )


# Two lanes of 14 px, a lane gap of 3 px and padding of 2 px on each side. The room holds the letters (10 px), the
# label offset and the stroke: half of it is the label's own mark past its lane, half the mark it would cross.
_BASE = 2 * 2.0 + 2 * 14.0 + 3.0
_ROOM = 10.0 + _OFFSET + RegionAnnotationStyle().stroke_width


@pytest.mark.parametrize(("above", "join_origin"), ((False, False), (True, False), (False, True)))
def test_the_track_grows_by_the_room_between_the_two_lanes(above: bool, join_origin: bool) -> None:
    track = _packed(above=above, join_origin=join_origin)

    assert {item.annotation.id: item.lane for item in track.placements} == {"a": 0, "b": 1}
    assert track.lane_rooms == (LaneRoom(record_index=0, gap=0, px=_ROOM),)
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
    assert track.lane_rooms == (LaneRoom(record_index=0, gap=0, px=_ROOM),)
