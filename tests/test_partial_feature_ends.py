"""Partial feature ends are drawn open: no arrowhead, no stroke on that edge.

A partial end is written in the source location (``<1..100``, ``200..>300``) or
by a gbdraw crop (``--region``, a records-table ``region``). It is drawn like an
artificial display seam: the feature continues past it.
"""
from __future__ import annotations

import math
import re
from typing import NamedTuple
from xml.etree import ElementTree as ET

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import (
    AfterPosition,
    BeforePosition,
    CompoundLocation,
    ExactPosition,
    SeqFeature,
    SimpleLocation,
)
from Bio.SeqRecord import SeqRecord

from gbdraw.api.config import apply_config_overrides
from gbdraw.api.options import CircularDiagramOptions, LinearDiagramOptions
from gbdraw.api.request_render import plan_request
from gbdraw.api.requests import (
    CircularDiagramRequest,
    InMemoryRecordSource,
    LinearDiagramRequest,
    RecordDisplayOptions,
    RecordInput,
    RecordPresentation,
)
from gbdraw.crop_genbank import crop_and_shift_features
from gbdraw.io.regions import parse_region_spec

LENGTH = 1000


class Block(NamedTuple):
    """One drawn block, ordered by display position: open low edge, open high edge, arrowhead."""

    open_low: bool
    open_high: bool
    arrow: bool


def _record(location) -> SeqRecord:
    record = SeqRecord(Seq("ACGT" * (LENGTH // 4)), id="partial", annotations={
        "topology": "circular", "molecule_type": "DNA",
    })
    record.features = [SeqFeature(location, type="CDS", qualifiers={"locus_tag": ["cds1"]})]
    return record


def _svg(mode: str, record: SeqRecord, *, start: int | None = None, reverse: bool = False,
         region: str | None = None) -> ET.Element:
    options_type = CircularDiagramOptions if mode == "circular" else LinearDiagramOptions
    request_type = CircularDiagramRequest if mode == "circular" else LinearDiagramRequest
    cfg = apply_config_overrides(None, {
        "labels.circular.scope": "none", "labels.linear.scope": "none",
        "canvas.show_gc": False, "canvas.show_skew": False,
        "objects.features.arrow_geometry.head_length_ratio": 1.0,
    })
    request = request_type(records=(RecordInput(
        InMemoryRecordSource(record),
        region=parse_region_spec(region) if region is not None else None,
        display=RecordDisplayOptions(start_coordinate=start),
        presentation=RecordPresentation(reverse_complement=reverse),
    ),), options=options_type(config=cfg))
    result = plan_request(request).build()
    drawing = result.drawing if hasattr(result, "drawing") else result
    return ET.fromstring(drawing.tostring())


_COMMAND = re.compile(r"([MLAz])([^MLAz]*)")
_NUMBER = re.compile(r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?")


def _segments(d: str) -> list[tuple[str, tuple[float, float], tuple[float, float]]]:
    segments: list[tuple[str, tuple[float, float], tuple[float, float]]] = []
    current = first = (0.0, 0.0)
    for command, args in _COMMAND.findall(d):
        numbers = [float(value) for value in _NUMBER.findall(args)]
        if command == "M":
            current = first = (numbers[0], numbers[1])
        elif command == "L":
            point = (numbers[0], numbers[1])
            segments.append(("L", current, point))
            current = point
        elif command == "A":
            point = (numbers[5], numbers[6])
            segments.append(("A", current, point))
            current = point
        else:
            segments.append(("L", current, first))
            current = first
    return segments


def _points(d: str) -> list[tuple[float, float]]:
    return [point for segment in _segments(d) for point in segment[1:]]


def _close(a: float, b: float) -> bool:
    return abs(a - b) < 1e-6


def _linear_block(fill: str, stroke: str) -> Block:
    """Edges are vertical stroked segments; the arrowhead is the vertex on the middle line."""
    points = _points(fill)
    ys = sorted({round(y, 6) for _, y in points})
    low, high = min(x for x, _ in points), max(x for x, _ in points)
    head = None
    if len(ys) == 3:
        head = next(x for x, y in points if _close(y, ys[1]))
    stroked = {
        round(a[0], 6) for kind, a, b in _segments(stroke)
        if kind == "L" and _close(a[0], b[0]) and not _close(a[1], b[1])
    }
    expected_closed = {round(x, 6) for x in (low, high) if head is None or not _close(x, head)}
    assert stroked <= expected_closed, (stroked, expected_closed)
    return Block(
        open_low=round(low, 6) not in stroked and (head is None or not _close(low, head)),
        open_high=round(high, 6) not in stroked and (head is None or not _close(high, head)),
        arrow=head is not None,
    )


# Path numbers have three decimals; on a ~350 px radius that is well under 0.05 bp.
_BP_TOLERANCE = 0.05


def _position(point: tuple[float, float], centre: float = LENGTH / 2) -> float:
    """The bp position of a point, unwrapped to within half a record of ``centre``."""
    angle = math.degrees(math.atan2(point[1], point[0]))
    position = (angle + 90.0) / 360.0 * LENGTH
    return centre + (position - centre + LENGTH / 2) % LENGTH - LENGTH / 2


def _centre(points: list[tuple[float, float]]) -> float:
    x = sum(px / math.hypot(px, py) for px, py in points)
    y = sum(py / math.hypot(px, py) for px, py in points)
    return _position((x, y)) % LENGTH


def _circular_block(fill: str, stroke: str) -> Block:
    """Edges are radial stroked segments; the arrowhead is the vertex on the centre radius."""
    points = _points(fill)
    radii: list[float] = []
    for radius in sorted(math.hypot(*p) for p in points):
        if not radii or radius - radii[-1] > 0.05:
            radii.append(radius)
    centre = _centre(points)

    def _at(point: tuple[float, float]) -> float:
        return _position(point, centre)

    positions = [_at(p) for p in points]
    low, high = min(positions), max(positions)
    head = None
    if len(radii) == 3:
        head = next(_at(p) for p in points if abs(math.hypot(*p) - radii[1]) <= 0.05)
    stroked = [
        _at(a) for kind, a, b in _segments(stroke)
        if kind == "L" and abs(_at(a) - _at(b)) < _BP_TOLERANCE
        and abs(math.hypot(*a) - math.hypot(*b)) > 0.05
    ]

    def is_stroked(position: float) -> bool:
        return any(abs(position - value) < _BP_TOLERANCE for value in stroked)

    def is_head(position: float) -> bool:
        return head is not None and abs(position - head) < _BP_TOLERANCE

    assert all(is_stroked(p) for p in stroked)
    assert not any(is_head(p) for p in stroked)
    return Block(
        open_low=not is_stroked(low) and not is_head(low),
        open_high=not is_stroked(high) and not is_head(high),
        arrow=head is not None,
    )


def _blocks(root: ET.Element, mode: str) -> list[Block]:
    nodes = [node for node in root.iter() if node.get("data-gbdraw-feature-part") == "block"]
    outlines = {node.get("id", "").removesuffix("__outline"): node.get("d", "")
                for node in nodes if node.get("id", "").endswith("__outline")}
    fills = [node for node in nodes if not node.get("id", "").endswith("__outline")]
    read = _linear_block if mode == "linear" else _circular_block
    blocks = []
    for node in fills:
        fill = node.get("d", "")
        outline = outlines.get(node.get("id", ""))
        if outline is not None:
            assert node.get("stroke") == "none"
        blocks.append((min(_points(fill)) if mode == "linear" else 0, read(fill, outline or fill), fill))
    if mode == "linear":
        blocks.sort(key=lambda item: item[0])
    else:
        blocks.sort(key=lambda item: _centre(_points(item[2])))
    return [block for _, block, _ in blocks]


def _reversed(blocks: list[Block]) -> list[Block]:
    return [Block(block.open_high, block.open_low, block.arrow) for block in reversed(blocks)]


def _loc(start, end, strand):
    return SimpleLocation(start, end, strand=strand)


CASES = {
    # A partial 5' end keeps the arrowhead and opens the tail.
    "plus_5p_partial": (_loc(BeforePosition(100), 400, 1), [Block(True, False, True)]),
    # A partial 3' end draws no arrowhead and opens that edge.
    "plus_3p_partial": (_loc(100, AfterPosition(400), 1), [Block(False, True, False)]),
    "minus_3p_partial": (_loc(BeforePosition(100), 400, -1), [Block(True, False, False)]),
    "minus_5p_partial": (_loc(100, AfterPosition(400), -1), [Block(False, True, True)]),
    "both_partial": (_loc(BeforePosition(100), AfterPosition(400), 1), [Block(True, True, False)]),
    "no_strand_partial": (_loc(BeforePosition(100), 400, None), [Block(True, False, False)]),
    "join_outer_partial": (
        CompoundLocation([_loc(BeforePosition(100), 200, 1), _loc(300, AfterPosition(400), 1)]),
        [Block(True, False, False), Block(False, True, False)],
    ),
    "minus_join_outer_partial": (
        CompoundLocation([_loc(300, AfterPosition(400), -1), _loc(BeforePosition(100), 200, -1)]),
        [Block(True, False, False), Block(False, True, False)],
    ),
    # Only the outer ends of a join can be partial; markers on inner ends are ignored.
    "join_inner_markers": (
        CompoundLocation([_loc(100, AfterPosition(200), 1), _loc(BeforePosition(300), 400, 1)]),
        [Block(False, False, False), Block(False, False, True)],
    ),
}


@pytest.mark.parametrize("mode", ["linear", "circular"])
@pytest.mark.parametrize("reverse", [False, True], ids=["forward", "reverse"])
@pytest.mark.parametrize("start", [None, 700], ids=["no_display_start", "display_start"])
@pytest.mark.parametrize("case", sorted(CASES))
def test_source_partial_ends_are_drawn_open(case, start, reverse, mode):
    location, expected = CASES[case]
    blocks = _blocks(_svg(mode, _record(location), start=start, reverse=reverse), mode)
    assert blocks == (_reversed(expected) if reverse else expected)


def test_a_complete_feature_keeps_its_arrowhead_and_closed_tail():
    for mode in ("linear", "circular"):
        assert _blocks(_svg(mode, _record(_loc(100, 400, 1))), mode) == [Block(False, False, True)]


@pytest.mark.parametrize("mode", ["linear", "circular"])
def test_a_partial_end_and_a_display_seam_are_both_open(mode):
    # The seam at 250 splits [<100..400]: the head fragment is open at the seam, the
    # tail fragment at both the partial end and the seam.
    blocks = _blocks(_svg(mode, _record(_loc(BeforePosition(100), 400, 1)), start=251), mode)
    assert sorted(blocks) == sorted([Block(True, False, True), Block(True, True, False)])


REGION_CASES = {
    # The crop cuts the 3' end of a plus-strand CDS: no arrowhead at the cut.
    "plus_cut_3p": (_loc(100, 400, 1), "1-250", [Block(False, True, False)]),
    "plus_cut_5p": (_loc(100, 400, 1), "251-1000", [Block(True, False, True)]),
    # On the minus strand the cut at the low end is the 3' end (OV-2273).
    "minus_cut_3p": (_loc(100, 400, -1), "251-1000", [Block(True, False, False)]),
    "minus_cut_5p": (_loc(100, 400, -1), "1-250", [Block(False, True, True)]),
    # A partial end written in the file survives the crop.
    "source_partial_kept": (_loc(BeforePosition(100), AfterPosition(400), 1), "50-800",
                            [Block(True, True, False)]),
}


@pytest.mark.parametrize("mode", ["linear", "circular"])
@pytest.mark.parametrize("rc", [False, True], ids=["forward", "rc"])
@pytest.mark.parametrize("case", sorted(REGION_CASES))
def test_crop_partial_ends_are_drawn_open(case, rc, mode):
    location, region, expected = REGION_CASES[case]
    spec = f"partial:{region}{':rc' if rc else ''}"
    blocks = _blocks(_svg(mode, _record(location), region=spec), mode)
    assert blocks == (_reversed(expected) if rc else expected)


def _crop(location, start_0, end_0):
    feature = SeqFeature(location, type="CDS")
    (cropped,) = crop_and_shift_features([feature], start_0, end_0)
    return cropped.location


def _kinds(location) -> list[tuple[type, type]]:
    return [(type(part.start), type(part.end)) for part in location.parts]


EXACT, BEFORE, AFTER = ExactPosition, BeforePosition, AfterPosition


@pytest.mark.parametrize(
    ("location", "window", "expected"),
    [
        # Overhangs of the part itself, on both strands (OV-2273: the minus strand lost them).
        (_loc(100, 300, 1), (200, 1000), [(BEFORE, EXACT)]),
        (_loc(100, 300, 1), (0, 200), [(EXACT, AFTER)]),
        (_loc(100, 300, -1), (200, 1000), [(BEFORE, EXACT)]),
        (_loc(100, 300, -1), (0, 200), [(EXACT, AFTER)]),
        (_loc(100, 300, None), (200, 1000), [(BEFORE, EXACT)]),
        # Markers written in the file survive the shift (OV-2274).
        (_loc(BeforePosition(100), AfterPosition(300), 1), (50, 1000), [(BEFORE, AFTER)]),
        (_loc(BeforePosition(100), AfterPosition(300), -1), (50, 1000), [(BEFORE, AFTER)]),
        # A dropped neighbouring part makes the kept outer end partial, in biological order.
        (CompoundLocation([_loc(100, 200, 1), _loc(300, 400, 1)]), (250, 1000), [(BEFORE, EXACT)]),
        (CompoundLocation([_loc(100, 200, 1), _loc(300, 400, 1)]), (0, 250), [(EXACT, AFTER)]),
        (CompoundLocation([_loc(300, 400, -1), _loc(100, 200, -1)]), (250, 1000), [(BEFORE, EXACT)]),
        (CompoundLocation([_loc(300, 400, -1), _loc(100, 200, -1)]), (0, 250), [(EXACT, AFTER)]),
        # A neighbour that ends exactly at the window edge is dropped too.
        (CompoundLocation([_loc(100, 200, 1), _loc(300, 400, 1)]), (200, 1000), [(BEFORE, EXACT)]),
    ],
)
def test_crop_marks_every_cut_end(location, window, expected):
    assert _kinds(_crop(location, *window)) == expected


def test_crop_coordinates_are_unchanged():
    cropped = _crop(_loc(BeforePosition(100), 300, -1), 50, 200)
    assert (int(cropped.start), int(cropped.end), cropped.strand) == (50, 150, -1)


def test_an_outline_with_one_open_end_is_one_joined_stroke():
    # OV-2275: the closed tail and the long edges are one subpath, so the tail corners join.
    from gbdraw.svg.arrows import block_path_with_open_caps

    path = "M 0,0 L 10,0 L 10,5 L 0,5 z"
    _, _, outline = block_path_with_open_caps(path, end_cap="L 10,5", tail_cap="L 0,0")
    assert outline.count("M") == 1
    assert [point for segment in _segments(outline) for point in segment[1:]] == [
        (10, 5), (0, 5), (0, 5), (0, 0), (0, 0), (10, 0),
    ]
