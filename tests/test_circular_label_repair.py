"""Local repair of slight Circular label defects (LP-D50, LP-D57) and the yield after it (LP-D59).

Before the fit check counts a record, its final horizontal layout is repaired
when it has 1 to 100 defects: labels move a little, never past a neighbour,
and the repair is kept only when the record then fits. A repaired record then
yields room around its repaired defects, under the same rules.
"""

from __future__ import annotations

import copy
import functools
import json
import math
import warnings
from pathlib import Path

import pytest

import gbdraw.labels.circular_repair as repair_module
from gbdraw.circular import circular_main
from gbdraw.exceptions import LabelLeadersCrossWarning, LabelsDoNotFitError
from gbdraw.labels.circular import (
    _current_half_from_x,
    _label_bbox_min_radius,
    _label_bbox_rect,
    circular_label_defects,
    circular_label_fits,
    circular_label_shape,
)
from gbdraw.labels.circular_leaders import arrival_contact, chosen_contact_end, label_ink_box
from gbdraw.labels.circular_repair import defect_count, make_way_for_region_labels, repair_circular_labels
from gbdraw.labels.fit import LabelDefect
from gbdraw.layout.spatial import Aabb, bounding_aabb, boxed_parts_meet

TEST_INPUTS = Path(__file__).parent / "test_inputs"
FIXTURES = TEST_INPUTS / "label_repair"


def _fixture(name: str) -> tuple[list[dict], int, tuple[Aabb, ...]]:
    """A final layout as ``prepare_label_list`` hands it to the repair (written by the engine of this branch)."""

    data = json.loads((FIXTURES / f"{name}.json").read_text())
    # The layouts were recorded before a label carried its font style; each was drawn upright.
    labels = [{"font_style": "normal", **label} for label in data["labels"]]
    return labels, int(data["total_length"]), tuple(Aabb(*box) for box in data["definition_text"])


# Each fixture's repair (and yield) runs once per test process; tests read copies, so none can change a shared
# result. ``_Repair`` is taken at import, before any test replaces it.
_REPAIR = repair_module._Repair


@functools.cache
def _shared_search(name: str) -> repair_module._Repair:
    labels, total_length, definition = _fixture(name)
    search = _REPAIR(labels, total_length, definition)
    assert search.run()
    return search


def _searched(name: str) -> repair_module._Repair:
    """A copy of the fixture's search state right after the repair's own ``run`` (before the yield)."""

    return copy.deepcopy(_shared_search(name))


def _resume_after_repair(monkeypatch: pytest.MonkeyPatch, name: str) -> list[repair_module._Repair]:
    """``repair_circular_labels`` on this fixture starts the yield from a copy of the repaired search state.

    The entry point still counts the defects, checks the repaired layout and calls ``make_room``; only the
    repair's search, already run by ``_shared_search``, is not repeated. Returns the search states handed out.
    """

    handed: list[repair_module._Repair] = []

    def resumed(labels: list[dict], total_length: int, definition: tuple[Aabb, ...]) -> repair_module._Repair:
        search = _searched(name)
        assert list(labels) == search.original and total_length == search.total_length
        search.run = lambda: True
        handed.append(search)
        return search

    monkeypatch.setattr(repair_module, "_Repair", resumed)
    return handed


@functools.cache
def _shared_result(name: str) -> list[dict] | None:
    labels, total_length, definition = _fixture(name)
    with pytest.MonkeyPatch.context() as monkeypatch:
        _resume_after_repair(monkeypatch, name)
        return repair_circular_labels(labels, total_length, definition_text=definition)


def _result(name: str) -> list[dict] | None:
    """A copy of ``repair_circular_labels`` on the fixture: the repair, then the yield."""

    return copy.deepcopy(_shared_result(name))


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


LAMBDA_14 = [
    "--gbk", str(TEST_INPUTS / "NC_001416.gb"), "--separate_strands", "--track_type", "tuckin",
    "--labels", "out", "--legend", "right", "--label_font_size", "14",
]


@pytest.fixture(scope="module")
def lambda_14_svg(tmp_path_factory: pytest.TempPathFactory) -> Path:
    """Lambda at 14 pt, drawn once for the module (``_draw`` turns a leader-crossing warning into an error)."""

    return _draw(tmp_path_factory.mktemp("lambda_14"), LAMBDA_14)


@pytest.mark.circular
def test_lambda_crossings_are_repaired_by_an_elbow_staircase(lambda_14_svg: Path) -> None:
    """Lambda at 14 pt: 2 crossing leader pairs; three elbows move 0.5 to 1.5 px and it draws without a warning."""

    assert lambda_14_svg.is_file()


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
    """Plasmid NZ_AP024092.1 at 14 pt: 1 crossing; the leader meets another point of its label, nothing moves.

    This is the repair's own result, before the yield (LP-D59) moves its neighbours.
    """

    labels, total_length, definition = _fixture("o_plasmid37_f14")
    assert defect_count(labels, total_length, definition) == 1

    search = repair_module._Repair(labels, total_length, definition)
    assert search.run()
    repaired = search.labels

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

    repaired = _result("o_n4_default")

    assert repaired is not None
    assert defect_count(repaired, total_length, definition) == 0
    _assert_order_kept(labels, repaired)


# ---- rules -------------------------------------------------------------------


@pytest.mark.parametrize("name", ["o_n4_default", "o_plasmid37_f14", "o_lambda_f14"])
def test_repairs_stay_within_their_bounds(name: str) -> None:
    """Text moves at most 24 px, an elbow 12 px; a label keeps its half and never comes inside its own elbow circle."""

    labels, total_length, definition = _fixture(name)
    repaired = _result(name)
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


def test_a_repair_that_cannot_clear_every_defect_returns_none_and_leaves_the_labels_unchanged(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """All or nothing: out of budget, the repair returns None and the caller's labels are unchanged."""

    labels, total_length, definition = _fixture("o_n4_default")
    before = copy.deepcopy(labels)
    monkeypatch.setattr(repair_module, "MAX_EVALUATIONS", 50)

    assert repair_circular_labels(labels, total_length, definition_text=definition) is None
    assert labels == before


def test_one_inner_label_over_the_definition_text_is_repaired() -> None:
    """OV-736: one movable label can have a defect (it covers the centre text), so it is repaired."""

    labels, total_length, definition = _fixture("c_tobacco_both_f14")
    one = [next(label for label in labels if label["label_text"].startswith("RNA polymerase beta subunit"))]
    assert one[0]["is_inner"] and defect_count(one, total_length, definition) == 1

    repaired = repair_circular_labels(one, total_length, definition_text=definition)
    assert repaired is not None and defect_count(repaired, total_length, definition) == 0


def test_a_repair_that_cannot_finish_within_the_budget_stops_early() -> None:
    """A record whose moves clear defects too slowly to reach 0 within the budget stops at the first such move.

    Tobacco both sides at 14 pt (54 defects) clears one defect per about 1,500
    evaluations and ends stuck; without the projection it spends 8,895.
    """

    labels, total_length, definition = _fixture("c_tobacco_both_f14")
    search = repair_module._Repair(labels, total_length, definition)

    assert not search.run()
    assert search.evaluations <= 2_000
    assert repair_circular_labels(labels, total_length, definition_text=definition) is None


def test_records_with_more_defects_than_the_cap_are_left_alone(monkeypatch: pytest.MonkeyPatch) -> None:
    labels, total_length, definition = _fixture("o_n4_default")
    monkeypatch.setattr(repair_module, "MAX_DEFECTS", 54)
    monkeypatch.setattr(repair_module, "_Repair", None)  # never reached

    assert repair_circular_labels(labels, total_length, definition_text=definition) is None


def test_the_repair_is_deterministic() -> None:
    """Two independent runs of the repair and its yield give the same labels.

    Plasmid NZ_AP024092.1 at 14 pt is the cheapest fixture whose repair yields (more than 10 labels move).
    """

    labels, total_length, definition = _fixture("o_plasmid37_f14")
    first = repair_circular_labels(labels, total_length, definition_text=definition)
    second = repair_circular_labels(labels, total_length, definition_text=definition)
    assert first is not None and _moved(_repaired_only("o_plasmid37_f14")[1], first)
    assert json.dumps(first, sort_keys=True) == json.dumps(second, sort_keys=True)
    # The other tests read the shared result, which starts the yield from a copy of the repaired search state.
    assert json.dumps(first, sort_keys=True) == json.dumps(_result("o_plasmid37_f14"), sort_keys=True)


@pytest.mark.circular
def test_the_same_run_draws_the_same_bytes(tmp_path: Path, lambda_14_svg: Path) -> None:
    first = lambda_14_svg.read_bytes()
    second = _draw(tmp_path, LAMBDA_14).read_bytes()
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


# ---- the yield after a repair (LP-D59) -------------------------------------------


def _repaired_only(name: str) -> tuple[list[dict], list[dict], int, tuple[Aabb, ...]]:
    """(engine labels, the repair's own result before the yield, record length, definition text)."""

    labels, total_length, definition = _fixture(name)
    return labels, copy.deepcopy(_shared_search(name).labels), total_length, definition


def _through_pairs(labels: list[dict], total_length: int) -> set[tuple[int, int]]:
    """(leader of, text of) by label index: the repair's own font-box test, over the whole record.

    Indices, not texts: a record can repeat a text ("hypothetical protein"), and the repair keeps the label order.
    """

    search = _REPAIR(labels, total_length, ())
    return set(search._current_through(search.external))


def _moved(first: list[dict], second: list[dict]) -> list[int]:
    return [
        index
        for index, (old, new) in enumerate(zip(first, second))
        if not old.get("is_embedded")
        and max(
            math.dist((old["start_x"], old["start_y"]), (new["start_x"], new["start_y"])),
            math.dist((old["middle_x"], old["middle_y"]), (new["middle_x"], new["middle_y"])),
        )
        > 0.05
    ]


@pytest.mark.parametrize("name", ["o_plasmid37_f14", "o_lambda_f14"])
def test_the_labels_around_repaired_defects_yield_room(name: str) -> None:
    """The window (the original defects' labels and up to 12 neighbours on their side) moves; the record still fits."""

    labels, repaired, total_length, definition = _repaired_only(name)
    yielded = _result(name)

    assert yielded is not None
    assert defect_count(yielded, total_length, definition) == 0
    moved = _moved(repaired, yielded)
    assert len(moved) >= 10
    search = _searched(name)
    window = {index for run in search._window_runs(max(repair_module.YIELD_WINDOWS)) for index in run}
    assert set(moved) <= window
    _assert_order_kept(labels, yielded)


@pytest.mark.parametrize("name", ["o_plasmid37_f14", "o_lambda_f14", "o_n4_default"])
def test_the_yield_never_draws_a_new_leader_through_another_labels_text(name: str) -> None:
    _labels, repaired, total_length, _definition = _repaired_only(name)
    yielded = _result(name)

    assert yielded is not None
    assert _through_pairs(yielded, total_length) <= _through_pairs(repaired, total_length)


def test_a_lift_that_trades_one_leader_through_text_for_another_is_refused() -> None:
    """The preview's guard compared counts: a lift that removed two through pairs and added one passed.

    The yield compares pairs, so every such lift of the lambda 14 pt window is refused.
    """

    search = _searched("o_lambda_f14")
    runs = search._window_runs(max(repair_module.YIELD_WINDOWS))
    search.window = frozenset(index for run in runs for index in run)
    traded = 0
    otherwise_valid = 0
    for run in runs:
        for split in range(1, len(run)):
            for yielders, crowd in ((run[max(0, split - 4) : split], run[split : split + 8]), (run[split : split + 4], run[max(0, split - 8) : split][::-1])):
                for distance in repair_module.LIFT_STEPS_PX:
                    moved = search._lift(yielders, crowd, distance)
                    probe = {index: label.copy() for index, label in moved.items()}
                    if not probe or not all(search._takes_rule_contact(label) for label in probe.values()):
                        continue
                    for label in probe.values():
                        search._finish(label)
                    geometry = {index: search._geometry_of(label) for index, label in probe.items()}
                    before = search._current_through(list(moved))
                    after = search._through(probe, geometry)
                    if after - before and len(after) <= len(before):
                        traded += 1
                        otherwise_valid += (
                            all(search._within_bounds(index, label) for index, label in moved.items())
                            and search._order_kept(moved)
                            and search._row_gaps_kept(moved)
                            and not search._local(geometry) - search._current_local(list(moved))
                        )
                        assert search._try_yield(moved) is None
    assert traded > 0
    assert otherwise_valid > 0  # the through-pair test alone refuses these lifts


def test_a_yield_that_ends_with_a_defect_keeps_the_repaired_layout(monkeypatch: pytest.MonkeyPatch) -> None:
    """All or nothing: when the yield's result does not fit, the repair's own result is returned."""

    labels, repaired, total_length, definition = _repaired_only("o_plasmid37_f14")
    make_room = repair_module._Repair.make_room

    def overlapping(self: repair_module._Repair) -> bool:
        make_room(self)
        first, second = self.order["out"][:2]
        self.labels[second] = {**self.labels[second], "start_x": self.labels[first]["start_x"], "start_y": self.labels[first]["start_y"]}
        return True

    monkeypatch.setattr(repair_module._Repair, "make_room", overlapping)
    result = repair_circular_labels(labels, total_length, definition_text=definition)

    assert result is not None
    assert json.dumps(result, sort_keys=True) == json.dumps(repaired, sort_keys=True)


def test_a_yield_that_ends_with_a_new_leader_through_text_keeps_the_repaired_layout(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    labels, repaired, total_length, definition = _repaired_only("o_plasmid37_f14")
    monkeypatch.setattr(repair_module._Repair, "make_room", lambda self: False)

    result = repair_circular_labels(labels, total_length, definition_text=definition)

    assert result is not None
    assert json.dumps(result, sort_keys=True) == json.dumps(repaired, sort_keys=True)


def test_a_yield_move_that_hands_another_label_a_new_leader_through_text_is_undone(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """N4: a moved label can change the leader end the engine gives an unmoved label; that move is undone."""

    labels, repaired, total_length, definition = _repaired_only("o_n4_default")
    apply_yield = repair_module._Repair._apply_yield
    undone: list[bool] = []

    def recording(self: repair_module._Repair, moved: dict, allowed: frozenset) -> bool:
        before = json.dumps(self.labels, sort_keys=True)
        applied = apply_yield(self, moved, allowed)
        if not applied:
            undone.append(json.dumps(self.labels, sort_keys=True) == before)
        return applied

    monkeypatch.setattr(repair_module._Repair, "_apply_yield", recording)
    _resume_after_repair(monkeypatch, "o_n4_default")
    yielded = repair_circular_labels(labels, total_length, definition_text=definition)

    assert undone and all(undone)  # at least one undo, each back to the layout before the move
    assert yielded is not None
    assert defect_count(yielded, total_length, definition) == 0
    assert _through_pairs(yielded, total_length) <= _through_pairs(repaired, total_length)


def test_a_yield_that_keeps_a_new_leader_through_text_is_rejected_at_the_end(monkeypatch: pytest.MonkeyPatch) -> None:
    """``make_room`` itself checks the record's pairs: here the first move the undo refuses is kept, then nothing moves."""

    labels, repaired, total_length, definition = _repaired_only("o_n4_default")
    apply_yield = repair_module._Repair._apply_yield
    kept: list[dict] = []

    def keeping(self: repair_module._Repair, moved: dict, allowed: frozenset) -> bool:
        if kept:
            return False
        if apply_yield(self, moved, allowed):
            return True
        kept.append(moved)
        for index, label in moved.items():
            self.labels[index] = label
        self._refresh()
        return True

    make_room = repair_module._Repair.make_room
    made_room: list[bool] = []
    repaired_search = _searched("o_n4_default")
    allowed = repaired_search._current_through(repaired_search.external)
    monkeypatch.setattr(repair_module._Repair, "_apply_yield", keeping)
    monkeypatch.setattr(repair_module._Repair, "make_room", lambda self: made_room.append(make_room(self)) or made_room[-1])
    handed = _resume_after_repair(monkeypatch, "o_n4_default")

    result = repair_circular_labels(labels, total_length, definition_text=definition)

    (search,) = handed
    assert made_room == [False]
    assert kept and search._current_through(search.external) - allowed
    assert json.dumps(result, sort_keys=True) == json.dumps(repaired, sort_keys=True)


@pytest.mark.parametrize("name", ["fits", "c_tobacco_both_f14"])
def test_records_the_repair_does_not_change_never_yield(name: str, monkeypatch: pytest.MonkeyPatch) -> None:
    """LP-D59: a record that fits, or whose repair fails, gets no yield call (and no added time)."""

    if name == "fits":
        labels, repaired, total_length, definition = _repaired_only("o_plasmid37_f14")
        labels = repaired
    else:
        labels, total_length, definition = _fixture(name)
    calls: list[int] = []
    monkeypatch.setattr(repair_module._Repair, "make_room", lambda self: calls.append(1) or True)

    assert repair_circular_labels(labels, total_length, definition_text=definition) is None
    assert calls == []


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


# ---- feature labels make way for a region label (CIR-RGN-YIELD, OV-3017) -----


def _region_box_over(name: str, text: str) -> tuple[list[dict], int, tuple[Aabb, ...], Aabb]:
    """The fixture's labels and a region label's box over the letters of the label ``text``."""

    labels, total_length, definition = _fixture(name)
    box = label_ink_box(next(label for label in labels if label["label_text"] == text), total_length)
    assert box is not None
    return labels, total_length, definition, box


def _letters(box: Aabb) -> list[tuple[Aabb, tuple[tuple[float, float], ...]]]:
    """A region label's letters as one part with its box (``boxed_parts_meet``)."""

    return [(box, ((box.min_x, box.min_y), (box.max_x, box.min_y), (box.max_x, box.max_y), (box.min_x, box.max_y)))]


def _no_overlaps(_labels: object, _cleared: object) -> tuple[frozenset[int], ...]:
    return ()


def _over(labels: list[dict], total_length: int, box: Aabb) -> list[str]:
    """The labels whose letters meet a region label's letters ``box``, by the one test the placer asks too."""

    region = _letters(box)
    return [
        label["label_text"]
        for label in labels
        if not label.get("is_embedded")
        and boxed_parts_meet(
            [(part_box, part) for part in circular_label_shape(label, total_length).parts if (part_box := bounding_aabb(part))],
            region,
        )
    ]


def _shifts(before: list[dict], after: list[dict]) -> dict[str, float]:
    return {
        old["label_text"]: math.hypot(float(new["start_x"]) - float(old["start_x"]), float(new["start_y"]) - float(old["start_y"]))
        for old, new in zip(before, after)
        if not old.get("is_embedded")
    }


def test_feature_labels_make_way_for_a_region_label_keeping_side_and_order() -> None:
    """The labels over a region label's letters move off them along their side, in their order; the text moves at
    most the repair's 24 px plus the region label's height, and the figure gains no defect."""

    labels, total_length, definition, box = _region_box_over("o_lambda_f14", "Lom outer membrane protein")
    before = copy.deepcopy(labels)
    made_way = make_way_for_region_labels(
        labels, total_length, [_letters(box)], [0], region_overlaps=_no_overlaps, definition_text=definition
    )

    assert labels == before
    assert made_way is not None
    moved, cleared = made_way
    assert cleared == (0,)
    assert _over(before, total_length, box) and not _over(moved, total_length, box)
    _assert_order_kept(before, moved)
    assert [bool(label.get("is_inner")) for label in moved] == [bool(label.get("is_inner")) for label in before]
    assert defect_count(moved, total_length, definition) <= defect_count(before, total_length, definition)
    limit = repair_module.MAX_TEXT_SHIFT_PX + box.height + 1e-6
    assert all(shift <= limit for shift in _shifts(before, moved).values())
    assert any(old != new for old, new in zip(before, moved))


def test_no_label_moves_without_a_region_label_over_feature_letters() -> None:
    """A figure whose region labels meet no feature label's letters asks no label to make way: it is drawn as
    today (the assembler passes only the region labels drawn over feature letters as blocked)."""

    labels, total_length, definition, box = _region_box_over("o_lambda_f14", "Lom outer membrane protein")
    kwargs = {"region_overlaps": _no_overlaps, "definition_text": definition}
    assert make_way_for_region_labels(labels, total_length, [_letters(box)], [], **kwargs) is None
    assert make_way_for_region_labels(labels, total_length, [[]], [0], **kwargs) is None


@pytest.mark.parametrize(
    ("before_defects", "after_defects"),
    [
        ((), (LabelDefect("crossing", 3, 4),)),
        # OV-3080: one overlap traded for another keeps the count, but adds a pair the labels did not have.
        ((LabelDefect("overlap", 1, 2),), (LabelDefect("overlap", 3, 4),)),
        ((LabelDefect("overlap", 1, 2),), (LabelDefect("definition", 5, 5),)),
    ],
)
def test_a_way_that_adds_a_defect_is_not_kept(
    monkeypatch: pytest.MonkeyPatch, before_defects: tuple, after_defects: tuple
) -> None:
    """A move that clears the region label but adds an overlapping pair, a crossing pair or a label over the
    definition text it did not have is not kept, even when it removes another so that the count stays: the labels
    stay as they are and the region label is drawn and warned as today."""

    labels, total_length, definition, box = _region_box_over("o_lambda_f14", "Lom outer membrane protein")
    before = copy.deepcopy(labels)
    found = iter((before_defects, after_defects))
    monkeypatch.setattr(repair_module, "circular_label_defects", lambda *_args, **_kwargs: next(found))
    assert make_way_for_region_labels(
        labels, total_length, [_letters(box)], [0], region_overlaps=_no_overlaps, definition_text=definition
    ) is None
    assert labels == before


def test_the_fit_counts_are_the_defects_counted() -> None:
    """The fit check counts the defects it finds (``label_defects``), one owner for both (OV-3080)."""

    labels, total_length, definition, _box = _region_box_over("o_lambda_f14", "Lom outer membrane protein")
    fits = circular_label_fits(labels, total_length, definition_text=definition).values()
    defects = circular_label_defects(labels, total_length, definition_text=definition)
    assert [sum(getattr(fit, name) for fit in fits) for name in ("overlap_pairs", "crossing_pairs", "definition_overlaps")] == [
        sum(defect.kind == kind for defect in defects) for kind in ("overlap", "crossing", "definition")
    ]
    assert len(set(defects)) == len(defects) == defect_count(labels, total_length, definition)


def test_a_way_that_raises_what_a_region_label_overlaps_is_not_kept() -> None:
    """A move that clears one region label's letters but adds a label to what another region label overlaps, as its
    warning counts it (a leader line moved behind its letters, say), is not kept, even when it takes another label
    off that region label so that its count stays."""

    labels, total_length, definition, box = _region_box_over("o_lambda_f14", "Lom outer membrane protein")
    far = Aabb(box.min_x + 3000.0, box.min_y, box.max_x + 3000.0, box.max_y)
    counts = iter(((frozenset({1}), frozenset({3})), (frozenset(), frozenset({4}))))
    assert make_way_for_region_labels(
        labels, total_length, [_letters(box), _letters(far)], [0], region_overlaps=lambda _labels, _cleared: next(counts),
        definition_text=definition,
    ) is None


def test_the_shift_limit_holds_for_each_label_across_region_labels() -> None:
    """A label that made way for one region label and then makes way for another moves, in all, no further from
    where it was than 24 px plus the region label's height (CIR-RGN-YIELD): the repair for the second measures each
    shift from the place before the first (``origin``), not from where the first left it. Here the second region
    label lies over the middle of the band that "Lom outer membrane protein" crosses to clear the first: measured
    from where the first left it, the label would clear the second too, 49.43 px from where it was, past the
    48.32 px it may move."""

    labels, total_length, definition, box = _region_box_over("o_lambda_f14", "Lom outer membrane protein")
    lom = next(at for at, label in enumerate(labels) if label["label_text"] == "Lom outer membrane protein")
    first = make_way_for_region_labels(labels, total_length, [_letters(box)], [0], region_overlaps=_no_overlaps, definition_text=definition)
    assert first is not None
    moved_box = label_ink_box(first[0][lom], total_length)
    assert moved_box is not None
    band = Aabb(
        min(box.min_x, moved_box.min_x), min(box.min_y, moved_box.min_y), max(box.max_x, moved_box.max_x), max(box.max_y, moved_box.max_y)
    )
    middle = Aabb(band.min_x, band.min_y + 0.125 * band.height, band.max_x, band.max_y - 0.125 * band.height)
    made_way = make_way_for_region_labels(
        labels, total_length, [_letters(box), _letters(middle)], [0, 1], region_overlaps=_no_overlaps, definition_text=definition
    )

    assert made_way is not None and 0 in made_way[1]
    limit = repair_module.MAX_TEXT_SHIFT_PX + max(box.height, middle.height)
    assert limit == pytest.approx(48.32, abs=0.01)
    assert max(_shifts(labels, made_way[0]).values()) <= limit + 1e-6


def test_a_label_moves_by_the_height_of_the_region_labels_it_makes_way_for() -> None:
    """A label may move 24 px plus the height of the region label it makes way for (CIR-RGN-YIELD), the tallest
    of those it made way for itself when there are several: a taller region label cleared before by other labels
    gives it no more room. Each label's limit is the repair's (``_within_bounds``)."""

    labels, total_length, definition, box = _region_box_over("o_lambda_f14", "Lom outer membrane protein")
    lom = next(at for at, label in enumerate(labels) if label["label_text"] == "Lom outer membrane protein")
    # A taller region label over "Lom outer membrane protein", then one over a label on the far side of the circle.
    tall = Aabb(box.min_x, box.min_y - 4.0, box.max_x, box.max_y + 4.0)
    far = next(at for at, label in enumerate(labels) if label["label_text"] == "NinG/ Rap DNA junction specific endonuclease")
    short = label_ink_box(labels[far], total_length)
    assert short is not None and short.height < tall.height
    limits: list[tuple[float, ...]] = []

    class Recorded(_REPAIR):
        def __init__(self, *args, **kwargs) -> None:
            super().__init__(*args, **kwargs)
            limits.append(self.text_shift_px)

    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(repair_module, "_Repair", Recorded)
        made_way = make_way_for_region_labels(
            labels, total_length, [_letters(tall), _letters(short)], [0, 1], region_overlaps=_no_overlaps, definition_text=definition
        )

    assert made_way is not None and made_way[1] == (0, 1)
    first, second = limits
    assert set(first) == {repair_module.MAX_TEXT_SHIFT_PX + tall.height}
    assert second[far] == pytest.approx(repair_module.MAX_TEXT_SHIFT_PX + short.height)
    assert second[lom] == pytest.approx(repair_module.MAX_TEXT_SHIFT_PX + tall.height)


def test_a_label_makes_way_for_a_region_label_over_what_its_style_draws() -> None:
    """A region label over a feature label's underline or white box, not its letters, is over what the label draws
    (the placer's test): the label makes way for it as for its letters, and is not reported clear unmoved."""

    labels, total_length, definition, box = _region_box_over("o_lambda_f14", "Lom outer membrane protein")
    lom = next(label for label in labels if label["label_text"] == "Lom outer membrane protein")
    letters = bounding_aabb(point for part in circular_label_shape(lom, total_length).parts for point in part)
    assert letters is not None
    beside = Aabb(letters.max_x + 1.0, letters.min_y, letters.max_x + 0.5 * letters.width, letters.max_y)

    def drawn_parts(label: dict) -> tuple:
        # A white box reaching past the letters of this label by half their width.
        if label["label_text"] != "Lom outer membrane protein":
            return ()
        text = bounding_aabb(point for part in circular_label_shape(label, total_length).parts for point in part)
        assert text is not None
        right = text.max_x + 0.5 * text.width
        return (((text.min_x, text.min_y), (right, text.min_y), (right, text.max_y), (text.min_x, text.max_y)),)

    assert not _over([lom], total_length, beside)
    made_way = make_way_for_region_labels(
        labels, total_length, [_letters(beside)], [0], region_overlaps=_no_overlaps, definition_text=definition, drawn_parts=drawn_parts
    )
    if made_way is not None:
        moved = next(label for label in made_way[0] if label["label_text"] == "Lom outer membrane protein")
        assert moved != lom


def test_each_region_label_is_cleared_on_its_own_and_the_others_kept() -> None:
    """Region labels are taken in order, each with every other region label's letters kept clear: one the labels
    cannot clear (here letters over the whole ring of labels) is skipped while the next is cleared, and a way that
    would move a label over another region label's letters is not taken: here "tail terminator" clears a region
    label over its letters, but not at the place it takes when no other region label is drawn, where a second
    region label lies."""

    labels, total_length, definition, box = _region_box_over("o_lambda_f14", "Lom outer membrane protein")
    ring = Aabb(-2000.0, -2000.0, 2000.0, 2000.0)
    skipped = make_way_for_region_labels(
        labels, total_length, [_letters(box), _letters(ring)], [1, 0], region_overlaps=_no_overlaps, definition_text=definition
    )
    assert skipped is not None and skipped[1] == (0,)

    labels, total_length, definition, own = _region_box_over("o_lambda_f14", "tail terminator")
    alone = make_way_for_region_labels(labels, total_length, [_letters(own)], [0], region_overlaps=_no_overlaps, definition_text=definition)
    assert alone is not None
    taken = label_ink_box(next(label for label in alone[0] if label["label_text"] == "tail terminator"), total_length)
    assert taken is not None and not _over(labels, total_length, taken)

    made_way = make_way_for_region_labels(
        labels, total_length, [_letters(own), _letters(taken)], [0], region_overlaps=_no_overlaps, definition_text=definition
    )

    assert made_way is not None and made_way[1] == (0,)
    assert not _over(made_way[0], total_length, own)
    assert not _over(made_way[0], total_length, taken)
