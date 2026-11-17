"""Circular features slots outside the axis circle: each slot draws at its own radii, and the outer leaders of an
outward feature stack leave their feature from its outer edge."""

from __future__ import annotations

import math
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from tests.utils.cli_in_process import run_gbdraw

LAMBDA = Path(__file__).parent / "test_inputs" / "NC_001416.gb"
SVG = "{http://www.w3.org/2000/svg}"
_NUMBER = r"-?\d+(?:\.\d+)?(?:e-?\d+)?"


def _draw(tmp_path: Path, *slots: str, options: tuple[str, ...] = ("--labels", "out")) -> ET.Element:
    args = ["circular", "--gbk", str(LAMBDA), *options, "-o", "out", "-f", "svg"]
    for slot in slots:
        args += ["--circular_track_slot", slot]
    result = run_gbdraw(args, cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    return ET.parse(tmp_path / "out.svg").getroot()


def _path_points(d: str) -> list[tuple[float, float]]:
    """The vertices of a feature path: the points of its moves and lines and the end points of its arcs."""

    tokens = re.findall(rf"[MLAZmlaz]|{_NUMBER}", d)
    points: list[tuple[float, float]] = []
    command = ""
    index = 0
    while index < len(tokens):
        token = tokens[index]
        if token.isalpha():
            command = token
            index += 1
            continue
        if command == "A":
            points.append((float(tokens[index + 5]), float(tokens[index + 6])))
            index += 7
        else:
            points.append((float(token), float(tokens[index + 1])))
            index += 2
    return points


def _slot_radii(root: ET.Element, slot_id: str) -> list[float]:
    group = next(g for g in root.iter(f"{SVG}g") if g.get("data-gbdraw-slot-id") == slot_id)
    return [
        math.hypot(x, y)
        for path in group.iter(f"{SVG}path")
        if path.get("data-gbdraw-feature-id")
        for x, y in _path_points(path.get("d", ""))
    ]


def _slot_paths(root: ET.Element, slot_id: str) -> list[str]:
    group = next(g for g in root.iter(f"{SVG}g") if g.get("data-gbdraw-slot-id") == slot_id)
    return [path.get("d", "") for path in group.iter(f"{SVG}path") if path.get("data-gbdraw-feature-id")]


def test_a_second_features_slot_outside_draws_at_its_own_radii(tmp_path: Path) -> None:
    """OV-3626: the outside slot drew its features at the first slot's radii, one copy over the other."""

    root = _draw(tmp_path, "features:features", "outer:features@side=outside,lane_direction=outside")

    assert _slot_paths(root, "outer") != _slot_paths(root, "features")
    axis = float(next(root.iter(f"{SVG}circle")).get("r", "nan"))
    assert max(_slot_radii(root, "features")) <= axis + 1e-6
    assert min(_slot_radii(root, "outer")) >= axis - 1e-6


def _nearest_to_centre(x1: float, y1: float, x2: float, y2: float) -> float:
    """The distance from the circle's centre to the nearest point of the segment."""

    dx, dy = x2 - x1, y2 - y1
    t = max(0.0, min(1.0, -(x1 * dx + y1 * dy) / ((dx * dx + dy * dy) or 1.0)))
    return math.hypot(x1 + t * dx, y1 + t * dy)


@pytest.mark.parametrize(
    "slots",
    [
        pytest.param(("features:features@side=outside,lane_direction=outside",), id="one-slot-outside"),
        pytest.param(
            ("features:features", "outer:features@side=outside,lane_direction=outside"), id="second-slot-outside"
        ),
    ],
)
def test_outer_leaders_of_an_outward_stack_never_enter_its_lane(tmp_path: Path, slots: tuple[str, ...]) -> None:
    """OV-3627: the leaders ran from the feature's centre to an elbow inside the lane and from there at a shallow
    angle under the neighbouring features of the lane. Only the radial stub at the feature itself may lie in the
    lane; every other leader segment stays past the lane's outer edge. With two slots the labels name the outer copy
    (OV-3626), so the leaders never cross it either."""

    root = _draw(tmp_path, *slots)
    lane_outer = max(_slot_radii(root, slots[-1].split(":", 1)[0]))
    segments = [
        tuple(float(line.get(key, "nan")) for key in ("x1", "y1", "x2", "y2"))
        for line in root.iter(f"{SVG}line")
        if line.get("data-label-feature-id")
    ]
    radial = [abs(math.remainder(math.atan2(y1, x1) - math.atan2(y2, x2), math.tau)) < 1e-3 for x1, y1, x2, y2 in segments]

    assert len(segments) > 40 and not all(radial)
    entering = [segment for segment, is_radial in zip(segments, radial) if not is_radial and _nearest_to_centre(*segment) < lane_outer - 0.5]
    assert entering == []


def _outline(d: str) -> list[tuple[float, float]]:
    """A feature path as a polygon: its moves and lines, and each arc sampled along its circle."""

    points: list[tuple[float, float]] = []
    tokens = re.findall(rf"[MLAZmlaz]|{_NUMBER}", d)
    command = ""
    index = 0
    while index < len(tokens):
        if tokens[index].isalpha():
            command = tokens[index]
            index += 1
            continue
        if command == "A":
            x, y = float(tokens[index + 5]), float(tokens[index + 6])
            if points:
                x0, y0 = points[-1]
                r0, r1 = math.hypot(x0, y0), math.hypot(x, y)
                a0 = math.atan2(y0, x0)
                turn = math.remainder(math.atan2(y, x) - a0, math.tau)
                points += [
                    ((r0 + (r1 - r0) * k / 40) * math.cos(a0 + turn * k / 40), (r0 + (r1 - r0) * k / 40) * math.sin(a0 + turn * k / 40))
                    for k in range(1, 40)
                ]
            points.append((x, y))
            index += 7
        else:
            points.append((float(tokens[index]), float(tokens[index + 1])))
            index += 2
    return points


def _inside(point: tuple[float, float], polygon: list[tuple[float, float]]) -> bool:
    """Whether ``point`` lies inside the polygon (even-odd rule)."""

    x, y = point
    edges = zip(polygon, polygon[1:] + polygon[:1])
    return sum(1 for (x1, y1), (x2, y2) in edges if (y1 > y) != (y2 > y) and x < x1 + (y - y1) * (x2 - x1) / (y2 - y1)) % 2 == 1


def _gap_to_outline(point: tuple[float, float], polygon: list[tuple[float, float]]) -> float:
    """The distance from ``point`` to the polygon's outline, inside or outside it: a point in the middle of a feature
    is as far from its outline as one in the white space beside it."""

    x, y = point
    edges = list(zip(polygon, polygon[1:] + polygon[:1]))
    gaps = []
    for (x1, y1), (x2, y2) in edges:
        dx, dy = x2 - x1, y2 - y1
        t = max(0.0, min(1.0, ((x - x1) * dx + (y - y1) * dy) / ((dx * dx + dy * dy) or 1.0)))
        gaps.append(math.hypot(x1 + t * dx - x, y1 + t * dy - y))
    return min(gaps)


@pytest.mark.parametrize(
    ("slots", "options"),
    [
        pytest.param((), ("--labels", "out", "--track_type", "spreadout"), id="spreadout"),
        pytest.param(("features:features@side=outside,lane_direction=outside",), ("--labels", "out"), id="one-slot-outside"),
        pytest.param(
            ("features:features", "outer:features@side=outside,lane_direction=outside"), ("--labels", "out"), id="two-slots"
        ),
        pytest.param((), ("--labels", "out", "--resolve_overlaps"), id="resolve-overlaps"),
        pytest.param((), ("--labels", "both", "--resolve_overlaps", "--track_type", "spreadout"), id="spreadout-both"),
    ],
)
def test_every_leader_starts_on_its_feature(tmp_path: Path, slots: tuple[str, ...], options: tuple[str, ...]) -> None:
    """A leader that starts at its lane's edge ended in the white space beside an arrowhead, where the head narrows
    toward its tip: it seemed to name nothing (OV-3641, and with the lane-clearance rule of an outward lane). Every
    leader starts on the outline its feature is drawn with at the label's angle: not short of it, not inside the
    feature, and on the copy the label names (the outermost features slot here; with two slots both copies share
    their feature ids, so only the named slot's paths count)."""

    root = _draw(tmp_path, *slots, options=options)
    named_slot = (slots[-1] if slots else "features:features").split(":", 1)[0]
    group = next(g for g in root.iter(f"{SVG}g") if g.get("data-gbdraw-slot-id") == named_slot)
    outlines: dict[str, list[list[tuple[float, float]]]] = {}
    for path in group.iter(f"{SVG}path"):
        if path.get("data-gbdraw-feature-id"):
            outlines.setdefault(path.get("data-gbdraw-feature-id", ""), []).append(_outline(path.get("d", "")))
    ends: dict[str, list[tuple[float, float]]] = {}
    for line in root.iter(f"{SVG}line"):
        x1, y1, x2, y2 = (float(line.get(key, "nan")) for key in ("x1", "y1", "x2", "y2"))
        ends.setdefault(line.get("data-label-feature-id", ""), []).extend([(x1, y1), (x2, y2)])
    ends.pop("", None)

    gaps = {
        feature_id: round(min(_gap_to_outline(point, polygon) for point in points for polygon in outlines[feature_id]), 2)
        for feature_id, points in ends.items()
    }
    assert len(gaps) > 40
    assert {feature_id: gap for feature_id, gap in gaps.items() if gap > 0.5} == {}


@pytest.mark.parametrize(
    "slots",
    [
        pytest.param(
            ("outer:features@side=outside,lane_direction=outside", "ticks:ticks", "features:features"), id="outer-first"
        ),
        pytest.param(
            ("features:features", "outer:features@side=outside,lane_direction=outside", "ticks:ticks@side=outside"),
            id="inner-first",
        ),
    ],
)
def test_radial_label_room_reads_the_lanes_its_labels_name(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, slots: tuple[str, ...]
) -> None:
    """With two features slots and a scale between them, the radial inner labels' reserve read the first features
    slot while the inner labels name the innermost one: with the outer slot first the reserve lay past the inner
    lanes and the figure was refused at any font size (OV-3642). The outer labels' scale floor likewise counted a
    scale inside the outer slot when the inner slot came first. Both now read the lanes the labels name."""

    import gbdraw.diagrams.circular.assemble as assemble_module

    placements = []
    labels_of = assemble_module.circular_labels

    def placing(*args, **kwargs):
        placements.append(labels_of(*args, **kwargs))
        return placements[-1]

    monkeypatch.setattr(assemble_module, "circular_labels", placing)
    # 8.5 pt: the size the refusal of the 14 pt figure names; at 14 pt the innermost slot, inside the scale, leaves
    # the radial inner labels too little room, as a single features slot inside the scale does.
    _draw(tmp_path, *slots, options=("--labels", "both", "--label_placement", "radial", "--label_font_size", "8.5"))

    state = placements[-1].state
    outer, inner = state.label_lanes.outer, state.label_lanes.inner
    assert outer is not None and inner is not None and outer is not inner
    assert state.label_radii.inner_reserve_px <= float(inner.all_band_px.inner_px) + 1e-6
    floor = state.label_radii.outer_floor_px
    assert floor == 0.0 or floor >= float(outer.all_band_px.outer_px) - 1e-6


def _segment_meets_box(start: tuple[float, float], end: tuple[float, float], box) -> bool:
    """Whether the segment passes through the box (Liang-Barsky clipping)."""

    (x1, y1), (x2, y2) = start, end
    low, high = 0.0, 1.0
    for p, q in (
        (x1 - x2, x1 - box.min_x), (x2 - x1, box.max_x - x1), (y1 - y2, y1 - box.min_y), (y2 - y1, box.max_y - y1)
    ):
        if p == 0.0:
            if q < 0.0:
                return False
            continue
        t = q / p
        if p < 0.0:
            low = max(low, t)
        else:
            high = min(high, t)
        if low > high:
            return False
    return True


@pytest.mark.parametrize(
    ("slots", "options"),
    [
        pytest.param((), ("--labels", "out", "--track_type", "spreadout"), id="spreadout"),
        pytest.param(("features:features@side=outside,lane_direction=outside",), ("--labels", "out"), id="one-slot-outside"),
        pytest.param(
            ("features:features@side=outside,lane_direction=outside",), ("--labels", "both"), id="one-slot-outside-both"
        ),
        pytest.param(
            ("features:features", "outer:features@side=outside,lane_direction=outside"), ("--labels", "out"), id="two-slots"
        ),
    ],
)
def test_outer_label_text_stays_off_the_features(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, slots: tuple[str, ...], options: tuple[str, ...]
) -> None:
    """With the outer labels' lanes outside the axis circle, some outer label text sat on features near the bottom
    of the map (OV-3627). The painted letters of every outer label meet no feature of the lanes the labels name: no
    feature outline crosses a label's box, and no box lies inside a feature."""

    import gbdraw.diagrams.circular.assemble as assemble_module
    from gbdraw.labels.circular_leaders import label_ink_box

    placements = []
    labels_of = assemble_module.circular_labels

    def placing(*args, **kwargs):
        placements.append(labels_of(*args, **kwargs))
        return placements[-1]

    monkeypatch.setattr(assemble_module, "circular_labels", placing)
    root = _draw(tmp_path, *slots, options=options)
    placed = placements[-1]
    total_length = placed.state.total_length
    named_slot = (slots[-1] if slots else "features:features").split(":", 1)[0]
    group = next(g for g in root.iter(f"{SVG}g") if g.get("data-gbdraw-slot-id") == named_slot)
    outlines = [_outline(path.get("d", "")) for path in group.iter(f"{SVG}path") if path.get("data-gbdraw-feature-id")]
    boxes = [
        (label["label_text"], label_ink_box(label, total_length))
        for label in placed.labels
        if not label.get("is_embedded") and not label.get("is_inner")
    ]

    assert len(boxes) > 20 and all(box is not None for _text, box in boxes)
    on_features = sorted(
        text
        for text, box in boxes
        if any(
            _segment_meets_box(start, end, box)
            for polygon in outlines
            for start, end in zip(polygon, polygon[1:] + polygon[:1])
        )
        or any(
            _inside(corner, polygon)
            for polygon in outlines
            for corner in ((box.min_x, box.min_y), (box.min_x, box.max_y), (box.max_x, box.min_y), (box.max_x, box.max_y))
        )
    )
    assert on_features == []


def test_an_embedded_label_fits_its_feature_at_the_radius_it_is_drawn_at(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
) -> None:
    """OV-3646: the embed test measured a label at the axis radius, so on a lane well inside the axis (the inner copy
    of two features slots) "AAA family ATPase", "integrase" and "ea47" were embedded and ran past their features.
    Each embedded label's span, at the lane centre its text is drawn on, holds the label's width."""

    from gbdraw.diagrams.circular import assemble

    measured: list[list[dict]] = []
    original = assemble.prepare_label_list

    def capture(*args, **kwargs):
        labels = original(*args, **kwargs)
        measured.append(labels)
        return labels

    monkeypatch.setattr(assemble, "prepare_label_list", capture)
    _draw(
        tmp_path,
        "outer:features@side=outside,lane_direction=outside",
        "ticks:ticks",
        "features:features",
        options=("--labels", "both", "--label_placement", "radial", "--label_font_size", "8.5"),
    )

    embedded = [label for label in measured[-1] if label["is_embedded"]]
    assert {label["feature_center_radius_px"] for label in embedded} != set()
    record_length = 48_502
    short = {
        label["label_text"]: round(
            2 * math.pi * float(label["feature_center_radius_px"]) * (float(label["end"]) - float(label["start"]))
            / record_length - float(label["width_px"]), 1,
        )
        for label in embedded
        if 2 * math.pi * float(label["feature_center_radius_px"]) * (float(label["end"]) - float(label["start"]))
        / record_length < float(label["width_px"])
    }
    assert short == {}
