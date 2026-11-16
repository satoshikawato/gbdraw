"""The order-preserving packer of Circular horizontal labels (CP-D04, CP-D05) and the counts it reads."""

from __future__ import annotations

import copy
import functools
import json
import math
import warnings
from pathlib import Path
from typing import cast

import pytest
from Bio import SeqIO

import gbdraw.diagrams.circular.assemble as assemble_module
import gbdraw.labels.circular_packer as packer_module
import gbdraw.labels.fit as fit_module
from gbdraw.api import (
    CircularDiagramOptions,
    CircularDiagramRequest,
    InMemoryRecordSource,
    RecordInput,
    RenderOutputRequest,
)
from gbdraw.circular import circular_main
from gbdraw.exceptions import LabelLeadersCrossWarning, LabelsDoNotFitError
from gbdraw.layout.spatial import Aabb
from gbdraw.session import build_session_document
from gbdraw.session_io import materialize_embedded_file
from tests.utils.session_drawing_view import drawing_view
from tests.utils.session_files import embedded_session
from gbdraw.web_support.request_render import render_canonical_web_request
import gbdraw.labels.circular as circular_module
from gbdraw.labels.circular import _count_label_overlaps, _packer_side_rules, circular_label_fits
from gbdraw.labels.circular_capacity import LabelArc
from gbdraw.labels.circular_packer import PackerLook, SideRules, pack_side
from gbdraw.labels.circular_repair import MAX_DEFECTS, defect_count, repair_circular_labels, too_close_on_one_line
from gbdraw.labels.circular_leaders import leader_crossing_count, leader_polyline
from gbdraw.labels.circular_types import CircularLabel
from gbdraw.labels.fit import Polyline, crossing_pairs, iter_crossing_pairs, polylines_cross

TEST_INPUTS = Path(__file__).parent / "test_inputs"


def _box_label(start_x: float, start_y: float, width: float, **extra: float | bool) -> CircularLabel:
    """An outer label at 3 o'clock (middle baseline) whose text box runs from ``start_x`` to ``start_x + width``."""

    label = {
        "middle": 250.0,
        "start_x": start_x,
        "start_y": start_y,
        "width_px": width,
        "height_px": 10.0,
        "is_inner": False,
        **extra,
    }
    return cast(CircularLabel, label)


# ---- OV-730: the min-gap overlap count pads its broad phase by the largest pair gap ----------------------


def test_min_gap_overlap_count_finds_a_pair_wider_than_any_single_label_gap() -> None:
    """Labels 1 and 3 (not list neighbours) are 8 px apart; their pair gap is 9.2 px (5 px spacing + 4.2 px).

    Each label alone needs at most 5 px, so a broad phase padded by the single-label gap puts the two boxes in
    different grid cells and never tests the pair (OV-730). The pair-correct pad finds it.
    """

    spaced = {"_enforce_spacing_px": True, "_spacing_px": 5.0}
    labels = [
        _box_label(1000.0, 1000.0, 10.0, **spaced),
        _box_label(57.0, 0.0, 50.0, min_outer_start_radius_px=10.0),
        _box_label(2000.0, 1000.0, 10.0, **spaced),
        _box_label(115.0, 0.0, 50.0, **spaced),
    ]
    assert _count_label_overlaps(labels, 1000, use_min_gap=True) == 1


# ---- the fit check's crossing test, stopping at the first pair --------------------------------------------


REPAIR_FIXTURES = ("c_tobacco_both_f14", "o_lambda_f14", "o_n4_default", "o_plasmid37_f14")


def _external_repair_labels(name: str) -> list[CircularLabel]:
    """The horizontal external labels of a repair fixture (a final layout as ``prepare_label_list`` gives it)."""

    data = json.loads((TEST_INPUTS / "label_repair" / f"{name}.json").read_text())
    return [label for label in data["labels"] if not label.get("is_embedded")]


@pytest.mark.parametrize("name", REPAIR_FIXTURES)
def test_iter_crossing_pairs_yields_every_crossing_pair_in_index_order(name: str) -> None:
    polylines = [leader_polyline(label) for label in _external_repair_labels(name)]
    every_pair = [
        (left, right)
        for left in range(len(polylines))
        for right in range(left + 1, len(polylines))
        if polylines_cross(polylines[left], polylines[right])
    ]
    assert every_pair
    assert list(iter_crossing_pairs(polylines)) == every_pair == crossing_pairs(polylines)


@pytest.mark.parametrize("name", REPAIR_FIXTURES)
def test_leader_crossing_count_stops_at_its_limit(
    name: str, monkeypatch: pytest.MonkeyPatch
) -> None:
    labels = _external_repair_labels(name)
    calls = 0
    original = fit_module.polylines_cross

    def counting(first: Polyline, second: Polyline) -> bool:
        nonlocal calls
        calls += 1
        return original(first, second)

    monkeypatch.setattr(fit_module, "polylines_cross", counting)
    total = leader_crossing_count(labels)
    every_call, calls = calls, 0
    assert leader_crossing_count(labels, limit=1) == 1
    assert calls < every_call
    assert leader_crossing_count(labels, limit=total + 1) == total
    assert leader_crossing_count(labels, limit=0) == 0

    # Without one label of each crossing pair, nothing crosses.
    polylines = [leader_polyline(label) for label in labels]
    dropped = {right for _left, right in crossing_pairs(polylines)}
    kept = [label for index, label in enumerate(labels) if index not in dropped]
    assert leader_crossing_count(kept) == leader_crossing_count(kept, limit=1) == 0


# ---- the packer on one side (fixtures: outer labels as `rearrange_labels_fc` receives them) ---------------


PACKER_FIXTURES = TEST_INPUTS / "label_packer"
# The three records need look 1: none packs with every label on its genome half. P1 14 pt's fewest crossings: 7.
PACKINGS = ["p1_f13", "tobacco_f12", "p1_f14"]


def _packer_fixture(name: str, every: int = 1) -> tuple[list[CircularLabel], int, LabelArc, float]:
    """The fixture's labels (every ``every``-th in genome order), genome length, ellipse and spacing."""

    data = json.loads((PACKER_FIXTURES / f"{name}.json").read_text())
    labels = sorted(data["labels"], key=lambda label: float(label["middle"]))[::every]
    return labels, int(data["total_length"]), LabelArc(*data["arc"]), float(data["spacing_px"])


def _outer_rules(total_length: int, spacing: float, max_crossings: int = MAX_DEFECTS) -> SideRules:
    """The placement's packer rules for the outer side, as `prepare_label_list` passes them."""

    return _packer_side_rules(total_length, spacing, is_outer=True, max_crossings=max_crossings)


def _off_their_half(labels: list[CircularLabel], total_length: int) -> list[str]:
    """The labels whose anchor is not on their feature's genome half (first half right of the axis)."""

    return [
        str(label["label_text"])
        for label in labels
        if (float(label["start_x"]) > 0) != (float(label["middle"]) < total_length / 2)
    ]


@functools.cache
def _shared_packing(name: str, look: PackerLook) -> list[CircularLabel] | None:
    labels, total_length, arc, spacing = _packer_fixture(name)
    return pack_side(labels, arc, look, _outer_rules(total_length, spacing))


def _clockwise(label: CircularLabel) -> float:
    """The angle of a label's anchor, clockwise from 12 o'clock."""

    return math.atan2(float(label["start_x"]), -float(label["start_y"])) % (2 * math.pi)


@pytest.mark.circular
@pytest.mark.parametrize("name", PACKINGS)
def test_a_packed_side_keeps_genome_order_and_its_boxes_apart(name: str) -> None:
    """Cyclic genome order kept; no overlap by the fit check (its box margins are the 0.1 px row and 1 px side gaps)."""

    labels, total_length, _arc, _spacing = _packer_fixture(name)
    before = copy.deepcopy(labels)
    packed = _shared_packing(name, "pole")
    assert labels == before  # the input is not changed
    assert packed is not None and len(packed) == len(labels)
    assert [label["middle"] for label in packed] == sorted(label["middle"] for label in labels)
    angles = [_clockwise(label) for label in packed]
    assert sum(1 for k in range(len(angles)) if angles[k] < angles[k - 1]) <= 1  # one turn, no label out of order
    (fit,) = circular_label_fits(packed, total_length).values()
    assert fit.overlap_pairs == 0
    assert fit.crossing_pairs == leader_crossing_count(packed) == (7 if name == "p1_f14" else 0)


@pytest.mark.circular
@pytest.mark.parametrize(
    ("name", "passing"),
    [
        ("p1_f13", ["SimB", "SimA", "ImcB", "ImcA"]),
        ("tobacco_f12", [
            "tRNA-His", "photosystem II protein D1", "photosystem II protein T", "photosystem II protein N",
            "photosystem II protein H",
        ]),
    ],
)
def test_look_1_moves_only_whole_end_runs_across_a_pole(name: str, passing: list[str]) -> None:
    """Look 1 first moves whole runs at the ends of a genome half; the labels right of the axis stay one run."""

    labels, total_length, arc, spacing = _packer_fixture(name)
    assert pack_side(labels, arc, "half", _outer_rules(total_length, spacing)) is None
    packed = _shared_packing(name, "pole")
    assert packed is not None
    right = [float(label["start_x"]) > 0 for label in packed]
    assert sum(1 for k in range(len(right)) if right[k] != right[k - 1]) == 2
    assert _off_their_half(packed, total_length) == passing


@pytest.mark.circular
@pytest.mark.parametrize("name", PACKINGS)
def test_look_2_keeps_every_label_on_its_genome_half(name: str) -> None:
    """Look 2 on every other label of the dense records: packed, nothing crosses, no label passes a pole."""

    labels, total_length, arc, spacing = _packer_fixture(name, every=2)
    packed = pack_side(labels, arc, "half", _outer_rules(total_length, spacing, 0))
    assert packed is not None and len(packed) == len(labels)
    assert _off_their_half(packed, total_length) == []
    (fit,) = circular_label_fits(packed, total_length).values()
    assert (fit.overlap_pairs, fit.crossing_pairs) == (0, 0)


def _same_line_gaps(labels: list[CircularLabel], total_length: int) -> list[tuple[float, int, int]]:
    """Per pair of labels whose boxes overlap vertically (one text line): the side gap of their boxes in px, indices."""

    boxes = [
        (circular_module._label_x_bounds(label, 0.0), circular_module._label_y_bounds(label, total_length, 0.0))
        for label in labels
    ]
    gaps: list[tuple[float, int, int]] = []
    for i, (first_x, first_y) in enumerate(boxes):
        for j in range(i + 1, len(boxes)):
            second_x, second_y = boxes[j]
            if min(first_y[1], second_y[1]) > max(first_y[0], second_y[0]):
                gap = max(second_x[0] - first_x[1], first_x[0] - second_x[1])
                gaps.append((gap, i, j))
    return gaps


@pytest.mark.circular
@pytest.mark.parametrize(("name", "look"), [("p1_f13", "pole"), ("tobacco_f12", "pole"), ("p1_f13", "half")])
def test_packed_labels_on_one_text_line_keep_a_readable_gap(name: str, look: PackerLook) -> None:
    """Two packed labels on one line, on one side or on both sides of 12 or 6 o'clock, stay ``LINE_GAP_EM`` apart.

    With a 1 px side gap, "Bor-like outer membrane protein" and "head protein" read as one phrase (OV-2036).
    """

    labels, total_length, arc, spacing = _packer_fixture(name, every=1 if look == "pole" else 2)
    packed = pack_side(labels, arc, look, _outer_rules(total_length, spacing, 0 if look == "half" else MAX_DEFECTS))
    assert packed is not None
    gaps = _same_line_gaps(packed, total_length)
    font_size = float(labels[0]["font_size"])
    assert {float(label["font_size"]) for label in labels} == {font_size}
    # Labels share lines on one side, and across the vertical axis at a pole.
    assert any((float(packed[i]["start_x"]) > 0) != (float(packed[j]["start_x"]) > 0) for _gap, i, j in gaps)
    gap, i, j = min(gaps)
    assert gap >= packer_module.LINE_GAP_EM * font_size - 1e-6, (gap, packed[i]["label_text"], packed[j]["label_text"])


@pytest.mark.circular
def test_the_run_splits_are_the_splits_whose_both_sides_fit() -> None:
    """The prefix tables of ``_HalfLook.states`` give the splits a check of every split by ``excess`` gives."""

    labels, total_length, arc, spacing = _packer_fixture("p1_f13")
    items = [packer_module._item(label, arc, _outer_rules(total_length, spacing)) for label in labels]
    line_gap = packer_module.line_gap_px(item.label for item in items)
    look = packer_module._HalfLook(items, arc, total_length, line_gap)
    ring = packer_module._Ring(arc, 0.1, line_gap)
    shifts = [shift for shift in range(-packer_module.MAX_RUN_SHIFT, packer_module.MAX_RUN_SHIFT + 1) if look.valid(shift)]
    every_split: list[tuple[int, int]] = []
    for top in shifts:
        bottom = next((bottom for bottom in sorted(shifts, key=abs) if look.excess(ring, top, bottom) <= 0), None)
        if bottom is not None:
            every_split.append((top, bottom))
    every_split.sort(key=lambda state: abs(state[0]) + abs(state[1]))
    assert look.states(ring) == every_split and every_split


def test_the_packer_is_deterministic() -> None:
    labels, total_length, arc, spacing = _packer_fixture("tobacco_f12")
    first = pack_side(labels, arc, "pole", _outer_rules(total_length, spacing, 0))
    second = pack_side(labels, arc, "pole", _outer_rules(total_length, spacing, 0))
    assert first is not None and first == second == _shared_packing("tobacco_f12", "pole")


def test_a_look_without_a_packing_under_the_crossing_limit_gives_none() -> None:
    """P1 14 pt: look 1's fewest crossings are 7 pairs, so with at most 6 allowed the look has no placement."""

    labels, total_length, arc, spacing = _packer_fixture("p1_f14")
    assert pack_side(labels, arc, "pole", _outer_rules(total_length, spacing, 6)) is None


# ---- a record that does not fit today tries the packer (CP-D04, CP-D05) -----------------------------------


DENSE = ["--separate_strands", "--track_type", "tuckin", "--labels", "out", "--legend", "right"]
P1 = ["--gbk", str(TEST_INPUTS / "NC_005856.1.gb"), *DENSE]
TOBACCO = ["--gbk", str(TEST_INPUTS / "NC_001879.gbk"), *DENSE]
T5 = ["--gbk", str(TEST_INPUTS / "NC_005859.1.gb"), *DENSE]


@pytest.fixture
def looks(monkeypatch: pytest.MonkeyPatch) -> list[tuple[PackerLook, bool]]:
    """The looks the packer tried, each with whether it gave a placement."""

    tried: list[tuple[PackerLook, bool]] = []

    def spy(labels: list[CircularLabel], arc: LabelArc, look: PackerLook, rules: SideRules) -> list[CircularLabel] | None:
        packed = pack_side(labels, arc, look, rules)
        if labels:  # an empty side (no inner labels) is packed as it is
            tried.append((look, packed is not None))
        return packed

    monkeypatch.setattr(circular_module, "pack_side", spy)
    return tried


def _draw(tmp_path: Path, args: list[str]) -> Path:
    """Draw with a leader-crossing warning turned into an error (a figure that does not fit stops)."""

    output = tmp_path / "out"
    with warnings.catch_warnings():
        warnings.simplefilter("error", LabelLeadersCrossWarning)
        circular_main([*args, "-o", str(output), "-f", "svg"])
    return output.with_suffix(".svg")


@pytest.mark.circular
@pytest.mark.parametrize(
    "args",
    [
        pytest.param([*P1, "--label_font_size", "14"], id="p1_14pt_stopped_today"),
        pytest.param([*P1, "--label_font_size", "13"], id="p1_13pt_one_crossing_today"),
        pytest.param([*TOBACCO, "--label_font_size", "12"], id="tobacco_12pt_68_crossings_today"),
    ],
)
def test_a_record_that_does_not_fit_today_is_drawn_packed_in_look_1(
    tmp_path: Path, args: list[str], looks: list[tuple[PackerLook, bool]]
) -> None:
    """Today P1 14 pt stops (5 overlaps, 12 crossings) and P1 13 pt and tobacco 12 pt warn of crossings.

    Look 2 is tried first and has no placement (no packing keeps every label on its genome half); look 1
    then packs them without overlaps or crossings.
    """

    assert _draw(tmp_path, args).is_file()
    assert looks == [("half", False), ("pole", True)]


def test_a_packed_pair_on_one_text_line_closer_than_its_line_gap_is_too_close() -> None:
    """The rule every pass after the packer reads from the labels (``line_gap_px``): the repair and the defect count of
    a packed placement, region labels' make-way and the overlay legend's shift (OV-2037). Today's labels carry no gap."""

    packed = {"line_gap_px": 14.0}
    first = _box_label(100.0, 0.0, 40.0, **packed)
    assert too_close_on_one_line(first, _box_label(150.0, 0.0, 40.0, **packed), 1000)  # 10 px apart
    assert not too_close_on_one_line(first, _box_label(154.0, 0.0, 40.0, **packed), 1000)  # 14 px apart
    assert not too_close_on_one_line(first, _box_label(150.0, 10.0, 40.0, **packed), 1000)  # the next text line
    assert not too_close_on_one_line(first, _box_label(150.0, 0.0, 40.0, is_inner=True, **packed), 1000)
    assert not too_close_on_one_line(_box_label(100.0, 0.0, 40.0), _box_label(150.0, 0.0, 40.0), 1000)


@pytest.mark.circular
def test_the_repair_of_a_packed_record_keeps_the_line_gap(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """P1 14 pt is packed with 7 crossing pairs and repaired; two drawn labels on one line stay ``LINE_GAP_EM`` apart.

    With the repair's 1 px side gap, two "anti-repressor" labels at 6 o'clock ended 0.58 em apart (OV-2037).
    """

    drawn: list[tuple[int, list[CircularLabel]]] = []
    prepare = assemble_module.prepare_label_list

    def recording(feature_dict: object, total_length: int, *args: object, **kwargs: object) -> list[CircularLabel]:
        labels = prepare(feature_dict, total_length, *args, **kwargs)  # type: ignore[arg-type]
        if not kwargs.get("probe"):
            drawn.append((int(total_length), labels))
        return labels

    monkeypatch.setattr(assemble_module, "prepare_label_list", recording)
    assert _draw(tmp_path, [*P1, "--label_font_size", "14"]).is_file()
    total_length, labels = drawn[-1]
    external = [label for label in labels if not label.get("is_embedded")]
    font_size = float(external[0]["font_size"])
    assert {float(label["font_size"]) for label in external} == {font_size}
    gap, i, j = min(_same_line_gaps(external, total_length))
    assert gap >= packer_module.LINE_GAP_EM * font_size - 1e-6, (gap, external[i]["label_text"], external[j]["label_text"])


@pytest.mark.slow
@pytest.mark.circular
def test_t5_keeps_todays_placement_and_warning_when_no_look_packs(
    tmp_path: Path, looks: list[tuple[PackerLook, bool]]
) -> None:
    """T5 at the Auto size: no look packs under the crossing limit, so it is drawn as today, with today's counts."""

    with pytest.warns(LabelLeadersCrossWarning, match=r"0 label pairs overlap, 1713 leader pairs cross"):
        circular_main([*T5, "-o", str(tmp_path / "out"), "-f", "svg"])
    assert looks == [("half", False), ("pole", False)]


@pytest.mark.slow
@pytest.mark.circular
def test_t5_at_12pt_needs_more_room_than_its_ellipse_so_the_packer_is_not_tried(
    tmp_path: Path, looks: list[tuple[PackerLook, bool]]
) -> None:
    with pytest.raises(LabelsDoNotFitError):
        _draw(tmp_path, [*T5, "--label_font_size", "12"])
    assert looks == []


def test_the_repair_takes_the_callers_defect_count() -> None:
    data = json.loads((TEST_INPUTS / "label_repair" / "o_plasmid37_f14.json").read_text())
    # Recorded before a label carried its font style (as in test_circular_label_repair.py); drawn upright.
    labels = [{"font_style": "normal", **label} for label in data["labels"]]
    total_length = int(data["total_length"])
    definition = tuple(Aabb(*box) for box in data["definition_text"])
    defects = defect_count(labels, total_length, definition)
    assert defects == 1
    assert repair_circular_labels(labels, total_length, definition_text=definition, defects=defects) == (
        repair_circular_labels(labels, total_length, definition_text=definition)
    )


@pytest.mark.circular
def test_the_web_render_draws_p1_at_13pt_packed_without_a_fit_report(
    tmp_path: Path, looks: list[tuple[PackerLook, bool]]
) -> None:
    """Generate in the Web runs the same placement: P1 13 pt (one crossing today) is drawn with no label fit report."""

    record = SeqIO.read(str(TEST_INPUTS / "NC_005856.1.gb"), "genbank")
    overrides: dict[str, object] = {
        "canvas.strandedness": True,
        "canvas.circular.track_type": "tuckin",
        "labels.circular.scope": "outer",
        "labels.font_size.short": 13.0,
        "labels.font_size.long": 13.0,
    }
    built = build_session_document(
        CircularDiagramRequest(
            records=(RecordInput(source=InMemoryRecordSource(record)),),
            options=CircularDiagramOptions(config_overrides=overrides),
            output=RenderOutputRequest(output_prefix="p1", formats=("svg",)),
        )
    )
    document = drawing_view(embedded_session(built.to_dict(), built.resources.sources))
    resources = tmp_path / "resources"
    resources.mkdir()
    resource_paths = {
        resource_id: materialize_embedded_file(entry, temp_dir=resources, role=resource_id, prefix_role=False)
        for resource_id, entry in document["resources"].items()
    }
    with warnings.catch_warnings():
        warnings.simplefilter("error", LabelLeadersCrossWarning)
        response = render_canonical_web_request(
            document["renderRequest"], resource_paths=resource_paths, output_directory=tmp_path / "out"
        )
    assert response["results"] and not response["metadata"].get("labelFitReports")
    assert looks == [("half", False), ("pole", True)]
