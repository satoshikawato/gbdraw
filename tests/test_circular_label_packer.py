"""The order-preserving packer of Circular horizontal labels (CP-D04, CP-D05) and the counts it reads."""

from __future__ import annotations

import json
from pathlib import Path
from typing import cast

import pytest

import gbdraw.labels.fit as fit_module
from gbdraw.labels.circular import _count_label_overlaps
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
