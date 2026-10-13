"""Local repair of slight Circular label defects (LP-D50, LP-D57).

Before the fit check counts a record, its final horizontal layout is repaired
when it has 1 to 100 defects: labels move a little, never past a neighbour,
and the repair is kept only when the record then fits.
"""

from __future__ import annotations

import copy
import json
import math
import warnings
from pathlib import Path

import pytest

import gbdraw.labels.circular_repair as repair_module
from gbdraw.circular import circular_main
from gbdraw.exceptions import LabelLeadersCrossWarning, LabelsDoNotFitError
from gbdraw.labels.circular import _current_half_from_x, _label_bbox_min_radius, _label_bbox_rect
from gbdraw.labels.circular_leaders import arrival_contact, chosen_contact_end, label_ink_box
from gbdraw.labels.circular_repair import defect_count, repair_circular_labels
from gbdraw.layout.spatial import Aabb

TEST_INPUTS = Path(__file__).parent / "test_inputs"
FIXTURES = TEST_INPUTS / "label_repair"


def _fixture(name: str) -> tuple[list[dict], int, tuple[Aabb, ...]]:
    """A final layout as ``prepare_label_list`` hands it to the repair (written by the engine of this branch)."""

    data = json.loads((FIXTURES / f"{name}.json").read_text())
    return data["labels"], int(data["total_length"]), tuple(Aabb(*box) for box in data["definition_text"])


def _draw(tmp_path: Path, args: list[str]) -> Path:
    output = tmp_path / "out"
    with warnings.catch_warnings():
        warnings.simplefilter("error", LabelLeadersCrossWarning)
        circular_main([*args, "-o", str(output), "-f", "svg"])
    return output.with_suffix(".svg")


def _angle(label: dict) -> float:
    return math.degrees(math.atan2(float(label["start_y"]), float(label["start_x"])))


def _order_by_side(labels: list[dict]) -> dict[bool, list[tuple[str, float]]]:
    """Each side's labels in feature order, with the angle of their text."""

    sides: dict[bool, list[tuple[str, float]]] = {}
    for label in sorted(labels, key=lambda label: float(label["middle"])):
        if not label.get("is_embedded"):
            sides.setdefault(bool(label.get("is_inner")), []).append((label["label_text"], _angle(label)))
    return sides


def _assert_order_kept(before: list[dict], after: list[dict]) -> None:
    """Neighbours keep their clockwise order (a pair out of order before stays as it was)."""

    for side, old in _order_by_side(before).items():
        new = _order_by_side(after)[side]
        assert [name for name, _ in new] == [name for name, _ in old]
        for index in range(len(old)):
            following = (index + 1) % len(old)
            today = (old[following][1] - old[index][1]) % 360.0
            moved = (new[following][1] - new[index][1]) % 360.0
            if today <= 180.0:
                assert moved <= 180.0, (new[index][0], new[following][0])


# ---- failing-first cases (the study's fixtures) ------------------------------


@pytest.mark.circular
def test_lambda_crossings_are_repaired_by_an_elbow_staircase(tmp_path: Path) -> None:
    """Lambda at 14 pt: 2 crossing leader pairs; three elbows move 0.5 to 1.5 px and it draws without a warning."""

    svg = _draw(
        tmp_path,
        [
            "--gbk", str(TEST_INPUTS / "NC_001416.gb"), "--separate_strands", "--track_type", "tuckin",
            "--labels", "out", "--legend", "right", "--label_font_size", "14",
        ],
    )
    assert svg.is_file()


@pytest.mark.circular
def test_tobacco_middle_overlap_is_repaired_by_moving_toward_the_ring(tmp_path: Path) -> None:
    """Tobacco, `middle` tracks, labels both: 1 overlapping inner pair stopped the run; inner labels move toward the ring."""

    svg = _draw(
        tmp_path,
        ["--gbk", str(TEST_INPUTS / "NC_001879.gbk"), "--track_type", "middle", "--labels", "both"],
    )
    assert svg.is_file()


@pytest.mark.circular
def test_inner_labels_over_the_definition_text_are_repaired(tmp_path: Path) -> None:
    """Lambda at 9 pt, labels both: two inner labels covered the definition text (DEF-WARN); they move toward the ring."""

    svg = _draw(
        tmp_path,
        [
            "--gbk", str(TEST_INPUTS / "NC_001416.gb"), "--separate_strands", "--track_type", "tuckin",
            "--labels", "both", "--legend", "right", "--label_font_size", "9",
        ],
    )
    assert svg.is_file()


def test_plasmid_crossing_is_repaired_by_a_contact_change() -> None:
    """Plasmid NZ_AP024092.1 at 14 pt: 1 crossing; the leader meets another point of its label, nothing moves."""

    labels, total_length, definition = _fixture("o_plasmid37_f14")
    assert defect_count(labels, total_length, definition) == 1

    repaired = repair_circular_labels(labels, total_length, definition_text=definition)

    assert repaired is not None
    assert defect_count(repaired, total_length, definition) == 0
    changed = [new for old, new in zip(labels, repaired) if old.get("leader_contact") != new.get("leader_contact")]
    assert len(changed) == 1
    for old, new in zip(labels, repaired):
        if old.get("is_embedded"):
            continue
        assert (new["start_x"], new["start_y"], new["middle_x"], new["middle_y"]) == (
            old["start_x"], old["start_y"], old["middle_x"], old["middle_y"],
        )
    (contacted,) = changed
    box = label_ink_box(contacted, total_length)
    assert box is not None
    assert contacted["leader_contact"] == arrival_contact(contacted, box)
    assert (contacted["leader_start_x"], contacted["leader_start_y"]) == chosen_contact_end(contacted, total_length)


def test_n4_crossings_at_defaults_are_repaired() -> None:
    """Phage N4 at defaults: 55 crossing leader pairs, under the cap of 100; repaired within the evaluation budget."""

    labels, total_length, definition = _fixture("o_n4_default")
    assert defect_count(labels, total_length, definition) == 55

    repaired = repair_circular_labels(labels, total_length, definition_text=definition)

    assert repaired is not None
    assert defect_count(repaired, total_length, definition) == 0
    _assert_order_kept(labels, repaired)


# ---- rules -------------------------------------------------------------------


@pytest.mark.parametrize("name", ["o_n4_default", "o_plasmid37_f14"])
def test_repairs_stay_within_their_bounds(name: str) -> None:
    """Text moves at most 24 px, an elbow 12 px; a label keeps its half and never comes inside its own elbow circle."""

    labels, total_length, definition = _fixture(name)
    repaired = repair_circular_labels(labels, total_length, definition_text=definition)
    assert repaired is not None
    for old, new in zip(labels, repaired):
        if old.get("is_embedded"):
            assert new == old
            continue
        assert math.dist((old["start_x"], old["start_y"]), (new["start_x"], new["start_y"])) <= 24.0 + 1e-6
        assert math.dist((old["middle_x"], old["middle_y"]), (new["middle_x"], new["middle_y"])) <= 12.0 + 1e-6
        if _current_half_from_x(float(old["start_x"])) != 0:
            assert _current_half_from_x(float(new["start_x"])) == _current_half_from_x(float(old["start_x"]))
        if not new.get("is_inner"):
            elbow = math.hypot(new["middle_x"], new["middle_y"])
            assert _label_bbox_min_radius(new, total_length) >= min(
                _label_bbox_min_radius(old, total_length), math.hypot(old["middle_x"], old["middle_y"])
            ) - 1e-6
            assert elbow <= _label_bbox_min_radius(new, total_length) + 1e-6


def test_no_move_of_any_kind_passes_a_neighbour() -> None:
    """Every candidate the search may apply (elbow, toward the ring, away, contact, push, arc, ramp) keeps the order."""

    labels, total_length, definition = _fixture("o_n4_default")
    search = repair_module._Repair(labels, total_length, definition)
    cluster = min(repair_module._clusters(search.defects()), key=len)
    kinds: set[str] = set()
    for _tier, kind, moved in search._candidates(cluster):
        if not moved or not all(search._within_bounds(index, label) for index, label in moved.items()):
            continue
        if not search._order_kept(moved):
            continue
        kinds.add(kind)
        after = [moved.get(index, label) for index, label in enumerate(search.labels)]
        _assert_order_kept(search.labels, after)
    assert kinds >= {"elbow", "ring", "away", "arc"}


def test_a_repair_that_cannot_clear_every_defect_keeps_the_layout(monkeypatch: pytest.MonkeyPatch) -> None:
    """All or nothing: out of budget, the labels come back untouched (byte for byte in the drawing)."""

    labels, total_length, definition = _fixture("o_n4_default")
    before = copy.deepcopy(labels)
    monkeypatch.setattr(repair_module, "MAX_EVALUATIONS", 50)

    assert repair_circular_labels(labels, total_length, definition_text=definition) is None
    assert labels == before


def test_records_with_more_defects_than_the_cap_are_left_alone(monkeypatch: pytest.MonkeyPatch) -> None:
    labels, total_length, definition = _fixture("o_n4_default")
    monkeypatch.setattr(repair_module, "MAX_DEFECTS", 54)
    monkeypatch.setattr(repair_module, "_Repair", None)  # never reached

    assert repair_circular_labels(labels, total_length, definition_text=definition) is None


def test_the_repair_is_deterministic() -> None:
    labels, total_length, definition = _fixture("o_n4_default")
    first = repair_circular_labels(labels, total_length, definition_text=definition)
    second = repair_circular_labels(labels, total_length, definition_text=definition)
    assert first is not None
    assert json.dumps(first, sort_keys=True) == json.dumps(second, sort_keys=True)


@pytest.mark.circular
def test_the_same_run_draws_the_same_bytes(tmp_path: Path) -> None:
    args = [
        "--gbk", str(TEST_INPUTS / "NC_001416.gb"), "--separate_strands", "--track_type", "tuckin",
        "--labels", "out", "--legend", "right", "--label_font_size", "14",
    ]
    first = _draw(tmp_path / "a", args).read_bytes()
    second = _draw(tmp_path / "b", args).read_bytes()
    assert first == second


@pytest.mark.circular
def test_fit_probes_do_not_repair(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Q7: the fit search's probes place labels without the repair; the drawn placement repairs."""

    import gbdraw.labels.circular as circular_labels

    probes: list[bool] = []
    original = circular_labels.prepare_label_list
    calls: list[int] = []

    def counting(*args, probe: bool = False, **kwargs):
        probes.append(probe)
        return original(*args, probe=probe, **kwargs)

    def repair(labels, total_length, *, definition_text=()):
        calls.append(len(labels))
        return None

    monkeypatch.setattr("gbdraw.diagrams.circular.assemble.prepare_label_list", counting)
    monkeypatch.setattr(repair_module, "repair_circular_labels", repair)
    with pytest.raises(LabelsDoNotFitError):
        circular_main(
            [
                "--gbk", str(TEST_INPUTS / "NC_005856.1.gb"), "--separate_strands", "--track_type", "tuckin",
                "--labels", "out", "--legend", "right", "--label_font_size", "14",
                "-o", str(tmp_path / "out"), "-f", "svg",
            ]
        )
    assert True in probes
    assert len(calls) == probes.count(False)


# ---- contact points (LP-D23, LP-D25, LP-D55) -----------------------------------


def _label_at(elbow: tuple[float, float], feature: tuple[float, float]) -> dict:
    return {"middle_x": elbow[0], "middle_y": elbow[1], "feature_anchor_x": feature[0], "feature_anchor_y": feature[1]}


BOX = Aabb(100.0, 50.0, 180.0, 62.0)


def test_a_level_arrival_meets_the_middle_of_the_near_side() -> None:
    """tRNA-Asp, ribosomal protein S2: the last piece within 25 degrees of horizontal."""

    assert arrival_contact(_label_at((90.0, 54.0), (60.0, 40.0)), BOX) == (0.0, 0.5)
    assert arrival_contact(_label_at((190.0, 57.0), (220.0, 40.0)), BOX) == (1.0, 0.5)


def test_an_arrival_from_the_feature_side_below_meets_that_bottom_corner() -> None:
    """tRNA-Ser: the elbow lies below the box's middle band and the feature is below."""

    assert arrival_contact(_label_at((96.0, 60.0), (60.0, 90.0)), BOX) == (0.0, 1.0)


def test_an_elbow_above_the_box_meets_the_nearest_point_of_the_top_side() -> None:
    assert arrival_contact(_label_at((122.0, 40.0), (122.0, 0.0)), BOX) == (0.25, 0.0)


def test_the_ink_box_lies_inside_the_font_box() -> None:
    labels, total_length, _ = _fixture("o_plasmid37_f14")
    label = next(label for label in labels if not label.get("is_embedded"))
    ink = label_ink_box(label, total_length)
    assert ink is not None
    font = _label_bbox_rect(label, total_length, 0.0)
    assert font[0] - 1e-6 <= ink.min_x and ink.max_x <= font[2] + 1e-6
    assert ink.height < font[3] - font[1]
