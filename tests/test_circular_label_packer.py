"""The order-preserving packer of Circular horizontal labels (CP-D04, CP-D05) and the counts it reads."""

from __future__ import annotations

import copy
import functools
import json
import math
from pathlib import Path
from typing import cast

import pytest

import gbdraw.labels.fit as fit_module
from gbdraw.labels.circular import _count_label_overlaps, circular_label_fits
from gbdraw.labels.circular_capacity import LabelArc
from gbdraw.labels.circular_packer import PackerLook, pack_side
from gbdraw.labels.circular_repair import MAX_DEFECTS
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
# (fixture, look): the look's packing; P1 14 pt has no crossing-free packing in look 2 (its fewest: 3 pairs).
PACKINGS = [("p1_f13", "half"), ("tobacco_f12", "half"), ("p1_f14", "half"), ("p1_f13", "pole"), ("tobacco_f12", "pole")]


def _packer_fixture(name: str) -> tuple[list[CircularLabel], int, LabelArc, float]:
    data = json.loads((PACKER_FIXTURES / f"{name}.json").read_text())
    return data["labels"], int(data["total_length"]), LabelArc(*data["arc"]), float(data["spacing_px"])


@functools.cache
def _shared_packing(name: str, look: PackerLook) -> list[CircularLabel] | None:
    labels, total_length, arc, spacing = _packer_fixture(name)
    return pack_side(
        labels, arc, look, total_length=total_length, is_outer=True, spacing_px=spacing, max_crossings=MAX_DEFECTS
    )


def _clockwise(label: CircularLabel) -> float:
    """The angle of a label's anchor, clockwise from 12 o'clock."""

    return math.atan2(float(label["start_x"]), -float(label["start_y"])) % (2 * math.pi)


@pytest.mark.circular
@pytest.mark.parametrize(("name", "look"), PACKINGS)
def test_a_packed_side_keeps_genome_order_and_its_boxes_apart(name: str, look: PackerLook) -> None:
    """Cyclic genome order kept; no overlap by the fit check (its box margins are the 0.1 px row and 1 px side gaps)."""

    labels, total_length, _arc, _spacing = _packer_fixture(name)
    before = copy.deepcopy(labels)
    packed = _shared_packing(name, look)
    assert labels == before  # the input is not changed
    assert packed is not None and len(packed) == len(labels)
    assert [label["middle"] for label in packed] == sorted(label["middle"] for label in labels)
    angles = [_clockwise(label) for label in packed]
    assert sum(1 for k in range(len(angles)) if angles[k] < angles[k - 1]) <= 1  # one turn, no label out of order
    (fit,) = circular_label_fits(packed, total_length).values()
    assert fit.overlap_pairs == 0
    assert fit.crossing_pairs == leader_crossing_count(packed) == (3 if name == "p1_f14" else 0)


@pytest.mark.circular
@pytest.mark.parametrize("name", ["p1_f13", "tobacco_f12", "p1_f14"])
def test_look_2_keeps_each_side_one_run_of_genome_order(name: str) -> None:
    """Look 2: the labels right of the axis are one cyclic run of the genome order (end runs pass a pole whole)."""

    packed = _shared_packing(name, "half")
    assert packed is not None
    right = [float(label["start_x"]) > 0 for label in packed]
    assert sum(1 for k in range(len(right)) if right[k] != right[k - 1]) == 2


def test_the_packer_is_deterministic() -> None:
    labels, total_length, arc, spacing = _packer_fixture("tobacco_f12")
    first = pack_side(labels, arc, "pole", total_length=total_length, is_outer=True, spacing_px=spacing, max_crossings=0)
    second = pack_side(labels, arc, "pole", total_length=total_length, is_outer=True, spacing_px=spacing, max_crossings=0)
    assert first is not None and first == second == _shared_packing("tobacco_f12", "pole")


def test_a_look_without_a_packing_under_the_crossing_limit_gives_none() -> None:
    """P1 14 pt: look 2's fewest crossings are 3 pairs, so with no crossing allowed the look has no placement."""

    labels, total_length, arc, spacing = _packer_fixture("p1_f14")
    assert pack_side(labels, arc, "half", total_length=total_length, is_outer=True, spacing_px=spacing, max_crossings=2) is None
