"""Circular outer-label fit check: decision point, search, wording, CLI, and LP-D45 pins."""

from __future__ import annotations

import dataclasses
import hashlib
import json
import sys
import warnings
from collections.abc import Callable, Mapping
from pathlib import Path

import pytest
from Bio import SeqIO

import gbdraw.diagrams.circular.assemble as circular_assemble_module
from gbdraw.api.diagram import label_fit_row
from gbdraw.circular import circular_main
from gbdraw.config.models import GbdrawConfig
import gbdraw.labels.circular_fit as circular_fit_module
import gbdraw.labels.circular_repair as circular_repair_module
from gbdraw.exceptions import LabelLeadersCrossWarning, LabelsDoNotFitError
from gbdraw.labels.linear import (
    FeatureLabelRow,
    FeatureLabelRows,
    RecordRegionLabels,
    check_linear_label_fit,
    keeps_feature_label_rows,
    linear_label_fit,
)
from gbdraw.layout.spatial import Aabb
from tests.utils.label_fit_warnings import assert_no_label_warning
from gbdraw.labels.circular_fit import (
    CircularLabelLayout,
    CircularLabelPlacement,
    check_circular_label_fit,
)
from gbdraw.labels.fit import (
    FONT_SIZE_PATHS,
    PLACEMENT_PATH,
    SCOPE_PATH,
    SEPARATE_FEATURES_PATH,
    SPACING_PATH,
    X_OFFSET_PATH,
    Y_OFFSET_PATH,
    LabelFit,
    LabelFitRecord,
    LabelFitReport,
    LabelFitRow,
    LabelFitVerification,
    LabelShape,
    format_label_fit_report,
    label_fit_notice,
    label_fit_outcome,
    label_fits_by_kind,
    overlapping_pairs,
    polyline_crossing_count,
    row_note,
)

TEST_INPUTS = Path(__file__).parent / "test_inputs"
FITS = LabelFit(0, 0)
OVERLAPS = LabelFit(3, 2)
CROSSES = LabelFit(0, 4)


def _settings(values: Mapping[str, float], *, font: float = 14.0) -> tuple[float, float, float, float]:
    return (
        values.get(X_OFFSET_PATH, 1.0),
        values.get(Y_OFFSET_PATH, 1.0),
        values.get(FONT_SIZE_PATHS[0], font),
        values.get(SPACING_PATH, 3.0),
    )


def _layout(
    record_id: str,
    fits: Callable[[float, float, float, float], bool],
    *,
    fit: LabelFit = OVERLAPS,
    probes: list[Mapping[str, float]] | None = None,
) -> CircularLabelLayout:
    def probe(values: Mapping[str, float]) -> LabelFit:
        if probes is not None:
            probes.append(dict(values))
        return FITS if fits(*_settings(values)) else OVERLAPS

    return CircularLabelLayout(
        record_id=record_id,
        record_index=0,
        placement="horizontal",
        font_size=14.0,
        label_counts={"circular_outer": 111},
        kind_fits={"circular_outer": fit},
        x_radius_offset=1.0,
        y_radius_offset=1.0,
        label_spacing=3.0,
        probe=probe,
    )


def _verifier(*layouts_fits: Callable[[float, float, float, float], bool]):
    calls: list[Mapping[str, float]] = []

    def verify(values: Mapping[str, float]) -> LabelFitVerification:
        calls.append(dict(values))
        settings = _settings(values)
        return LabelFitVerification(
            outcome="fits" if all(fits(*settings) for fits in layouts_fits) else "stop",
            canvas_size=(1000.0 * settings[0], 1000.0 * settings[1]),
        )

    return verify, calls


def _p1_like(x: float, y: float, font: float, spacing: float) -> bool:
    """Fits from Y 1.15, at 12 pt, or with 1.05 of Y per 0.5 pt removed; X and spacing never help."""

    return (y - 1.0) / 0.025 + (14.0 - font) / 0.5 * 1.5 >= 6 - 1e-9


def test_outcome_is_one_decision_point() -> None:
    assert label_fit_outcome([FITS, FITS]) == "fits"
    assert label_fit_outcome([FITS, CROSSES]) == "warn"
    assert label_fit_outcome([CROSSES, OVERLAPS]) == "stop"


def test_fit_counts_take_plain_boxes_and_polylines() -> None:
    """Boxes that touch overlap (closed boxes); leaders that only touch do not cross."""

    def box(min_x: float, min_y: float, max_x: float, max_y: float) -> tuple[tuple[float, float], ...]:
        return ((min_x, min_y), (max_x, min_y), (max_x, max_y), (min_x, max_y))

    boxes = [(box(0, 0, 10, 5),), (box(10, 5, 20, 9),), (box(10.01, 9.5, 20, 12),), (box(300, 300, 310, 305),)]
    assert overlapping_pairs(boxes, touching=True) == [(0, 1)]
    assert overlapping_pairs(boxes, touching=False) == []
    crossing = [((0.0, 0.0), (5.0, 5.0), (10.0, 5.0)), ((0.0, 5.0), (5.0, 0.0))]
    touching = [((0.0, 0.0), (5.0, 5.0)), ((5.0, 5.0), (10.0, 0.0))]
    assert polyline_crossing_count(crossing) == 1
    assert polyline_crossing_count(touching) == 0


@pytest.mark.parametrize("limit_first", [True, False])
@pytest.mark.parametrize("downward", [True, False])
def test_every_search_start_finds_the_mildest_fitting_value(limit_first: bool, downward: bool) -> None:
    for size in range(1, 12):
        grid = [float(index) for index in range(size)]
        for threshold in range(size + 1):
            expected = threshold if threshold < size else None
            for bound in (None, *range(size)):
                found = circular_fit_module._mildest_fitting_index(
                    grid,
                    lambda value: value >= threshold,
                    bound,
                    limit_first=limit_first,
                    downward=downward,
                )
                assert found == expected, (size, threshold, bound)


def test_fitting_records_pass_without_a_search() -> None:
    probes: list[Mapping[str, float]] = []
    verify, calls = _verifier(_p1_like)

    assert check_circular_label_fit([_layout("r", _p1_like, fit=FITS, probes=probes)], verify) is None
    assert probes == [] and calls == []


def test_overlaps_raise_with_six_verified_rows() -> None:
    verify, _calls = _verifier(_p1_like)

    with pytest.raises(LabelsDoNotFitError) as raised:
        check_circular_label_fit([_layout("P1", _p1_like)], verify)

    rows = {row.kind: row for row in raised.value.report.rows}
    assert [row.kind for row in raised.value.report.rows] == [
        "outward_both",
        "outward_vertical",
        "outward_horizontal",
        "smaller_font",
        "tighter_rows",
        "combination",
    ]
    assert rows["outward_vertical"].values == {Y_OFFSET_PATH: 1.15}
    assert rows["outward_vertical"].canvas_size == (1000.0, 1150.0)
    assert rows["outward_both"].values == {X_OFFSET_PATH: 1.15, Y_OFFSET_PATH: 1.15}
    assert rows["outward_horizontal"].values is None
    assert rows["outward_horizontal"].reason == "Does not fit up to X offset 2"
    assert rows["smaller_font"].values == {path: 12.0 for path in FONT_SIZE_PATHS}
    assert rows["tighter_rows"].reason == "Does not fit down to 0.1 px label spacing"
    # OD-F7: half the font row's 2 pt drop, then the smallest Y offset at 13 pt.
    assert rows["combination"].values == {
        **{path: 13.0 for path in FONT_SIZE_PATHS},
        Y_OFFSET_PATH: 1.075,
    }
    assert not raised.value.report.drawn


def test_row_steps_stronger_until_the_full_assembly_fits() -> None:
    def full_needs_more(x: float, y: float, font: float, spacing: float) -> bool:
        return y >= 1.2 - 1e-9 or font <= 11.0

    verify, calls = _verifier(full_needs_more)

    with pytest.raises(LabelsDoNotFitError) as raised:
        check_circular_label_fit([_layout("P1", _p1_like)], verify)

    rows = {row.kind: row for row in raised.value.report.rows}
    assert rows["outward_vertical"].values == {Y_OFFSET_PATH: 1.2}
    assert {Y_OFFSET_PATH: 1.15} in calls and {Y_OFFSET_PATH: 1.175} in calls
    assert rows["smaller_font"].values == {path: 11.0 for path in FONT_SIZE_PATHS}


def test_multi_record_rows_take_the_strongest_value_and_fit_every_record() -> None:
    def second(x: float, y: float, font: float, spacing: float) -> bool:
        return y >= 1.25 - 1e-9 or font <= 10.0

    verify, _calls = _verifier(_p1_like, second)
    layouts = [_layout("first", _p1_like), _layout("second", second, fit=CROSSES)]

    with pytest.raises(LabelsDoNotFitError) as raised:
        check_circular_label_fit(layouts, verify)

    report = raised.value.report
    assert [record.record_id for record in report.records] == ["first", "second"]
    rows = {row.kind: row for row in report.rows}
    assert rows["outward_vertical"].values == {Y_OFFSET_PATH: 1.25}
    assert rows["smaller_font"].values == {path: 10.0 for path in FONT_SIZE_PATHS}
    assert str(raised.value).startswith(
        "Outer labels do not fit in 2 records: first (111 labels at 14 pt: 3 label pairs overlap, "
        "2 leader pairs cross); second (111 labels at 14 pt: 0 label pairs overlap, 4 leader pairs cross)."
        " Nothing was drawn.\n"
    )


def test_crossing_leaders_alone_warn_and_return_the_report() -> None:
    verify, _calls = _verifier(_p1_like)

    with pytest.warns(LabelLeadersCrossWarning) as caught:
        checked = check_circular_label_fit([_layout("lambda", _p1_like, fit=CROSSES)], verify)

    assert checked is not None
    report, _search = checked
    assert report.drawn
    assert [row.status for row in report.rows] == [
        "found", "found", "not_found", "found", "not_found", "found"
    ]
    assert caught[0].message.report is report
    assert str(caught[0].message).startswith(
        "lambda: 111 outer labels do not fit at 14 pt (0 label pairs overlap, 4 leader pairs cross). Drawn anyway.\n"
    )


def test_wording_lists_fitting_rows_and_names_the_others() -> None:
    report = LabelFitReport(
        records=(LabelFitRecord("P1", "circular_outer", 111, 14.0, LabelFit(20, 21)),),
        rows=(
            LabelFitRow("outward_both", "found", {X_OFFSET_PATH: 1.15, Y_OFFSET_PATH: 1.15}, (1.0, 1.0)),
            LabelFitRow("outward_vertical", "found", {Y_OFFSET_PATH: 1.15}, (1.0, 1.0)),
            LabelFitRow("outward_horizontal", "not_found", reason="Does not fit up to X offset 2"),
            LabelFitRow("smaller_font", "found", {path: 12.0 for path in FONT_SIZE_PATHS}, (1.0, 1.0)),
            LabelFitRow("tighter_rows", "not_found", reason="Does not fit down to 0.1 px label spacing"),
            LabelFitRow(
                "combination",
                "found",
                {**{path: 13.0 for path in FONT_SIZE_PATHS}, Y_OFFSET_PATH: 1.05},
                (1.0, 1.0),
            ),
        ),
        drawn=False,
    )

    assert format_label_fit_report(report) == (
        "P1: 111 outer labels do not fit at 14 pt (20 label pairs overlap, 21 leader pairs cross). "
        "Nothing was drawn.\n"
        "Run again with one of:\n"
        "  --outer_label_x_radius_offset 1.15 --outer_label_y_radius_offset 1.15   (outward, both directions)\n"
        "  --outer_label_y_radius_offset 1.15                                      (outward, vertical only)\n"
        "  --label_font_size 12                                                    (smaller font)\n"
        "  --outer_label_y_radius_offset 1.05 --label_font_size 13                 (combination)\n"
        "Outward horizontal only and tighter rows do not fit this map."
    )


def test_warn_search_stops_at_its_time_limit_and_computes_a_row_later(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """LP-D48 SEARCH-BOUNDED-5S: rows past the limit are named, and one row is computed on request."""

    clock = [0.0]
    monkeypatch.setattr(circular_fit_module.time, "monotonic", lambda: clock[0])
    verify, _calls = _verifier(_p1_like)

    def slow(x: float, y: float, font: float, spacing: float) -> bool:
        clock[0] += 0.5
        return _p1_like(x, y, font, spacing)

    with pytest.warns(LabelLeadersCrossWarning) as caught:
        checked = check_circular_label_fit([_layout("lambda", slow, fit=CROSSES)], verify)

    assert checked is not None
    report, search = checked
    statuses = {row.kind: row.status for row in report.rows}
    assert statuses["outward_vertical"] == "found"
    assert statuses["combination"] == "not_computed"
    skipped = [row for row in report.rows if row.status == "not_computed"]
    assert all(row.values is None and row.reason is None for row in skipped)
    assert str(caught[0].message).splitlines()[-1].startswith("Not searched within 5 s: ")

    combination = search.row("combination")
    assert combination.status == "found"
    assert combination.values == {**{path: 13.0 for path in FONT_SIZE_PATHS}, Y_OFFSET_PATH: 1.075}


def test_the_warning_names_the_search_bound_the_caller_passed(monkeypatch: pytest.MonkeyPatch) -> None:
    """OV-708: a caller's bound (the Web passes 0 s) is the one the warning names, not the 5 s default."""

    monkeypatch.setattr(circular_fit_module.time, "monotonic", lambda: 0.0)
    verify, _calls = _verifier(_p1_like)

    with pytest.warns(LabelLeadersCrossWarning) as caught:
        check_circular_label_fit([_layout("lambda", _p1_like, fit=CROSSES)], verify, warn_search_seconds=0.0)

    assert str(caught[0].message).splitlines()[-1].startswith("Not searched within 0 s: ")


def test_overlaps_compute_every_row_without_a_time_limit(monkeypatch: pytest.MonkeyPatch) -> None:
    clock = [0.0]
    monkeypatch.setattr(circular_fit_module.time, "monotonic", lambda: clock[0])
    verify, _calls = _verifier(_p1_like)

    def slow(x: float, y: float, font: float, spacing: float) -> bool:
        clock[0] += 0.5
        return _p1_like(x, y, font, spacing)

    with pytest.raises(LabelsDoNotFitError) as raised:
        check_circular_label_fit([_layout("P1", slow)], verify)

    assert "not_computed" not in {row.status for row in raised.value.report.rows}


# --- inner, radial, definition text and Linear labels (LP-D48 SCOPE-ALL-LABELS, LP-D49) ------------


def _inner_layout(*, placement: CircularLabelPlacement = "horizontal", refused: bool = False, inner: int = 54) -> CircularLabelLayout:
    counts = {"circular_outer": 80, "circular_inner": inner} if inner else {"circular_outer": 80}
    kind = "circular_inner" if inner else "circular_outer"
    return CircularLabelLayout(
        record_id="P1",
        record_index=0,
        placement=placement,
        font_size=14.0,
        label_counts=counts,
        kind_fits={} if refused else {kind: LabelFit(3, 2, 1)},
        x_radius_offset=1.0,
        y_radius_offset=1.0,
        label_spacing=3.0,
        refused=refused,
    )


def _assembly(fits: Callable[[Mapping[str, object]], bool]):
    calls: list[dict[str, object]] = []

    def verify(values: Mapping[str, object]) -> LabelFitVerification:
        calls.append(dict(values))
        if fits(values):
            return LabelFitVerification("fits", (1000.0, 900.0))
        return LabelFitVerification("stop", None, LabelFit(1, 0))

    return verify, calls


def _p1_inner(values: Mapping[str, object]) -> bool:
    """P1 at 14 pt: inner labels fit at 12 pt and from 10 pt down, not at 11 pt (not monotone)."""

    font = values.get(FONT_SIZE_PATHS[0])
    return isinstance(font, float) and (font == 12.0 or font <= 10.0)


def test_inner_rows_verify_the_switches_first_then_scan_the_font_down() -> None:
    """INNER-3, SWITCHES-FIRST: the font row is the first font from the top that fits."""

    verify, calls = _assembly(_p1_inner)

    with pytest.raises(LabelsDoNotFitError) as raised:
        check_circular_label_fit([_inner_layout()], verify)

    assert calls == [
        {SCOPE_PATH: "outer"},
        {PLACEMENT_PATH: "radial"},
        *({path: font for path in FONT_SIZE_PATHS} for font in (13.5, 13.0, 12.5, 12.0)),
    ]
    rows = raised.value.report.rows
    assert [row.kind for row in rows] == ["smaller_font", "labels_out", "radial"]
    assert rows[0].values == {path: 12.0 for path in FONT_SIZE_PATHS}
    assert [row.status for row in rows] == ["found", "not_found", "not_found"]
    assert rows[1].reason == "Labels still do not fit (1 label pair overlaps, 0 leader pairs cross)"
    assert str(raised.value) == (
        "P1: 54 inner labels do not fit at 14 pt (3 label pairs overlap, 2 leader pairs cross, "
        "1 label covers the definition text). Nothing was drawn.\n"
        "Run again with one of:\n"
        "  --label_font_size 12   (smaller font)\n"
        "All labels outside and radial labels do not fit this map."
    )


def test_outer_levers_do_not_move_inner_labels() -> None:
    layouts = [
        _layout("A", _p1_like),
        _inner_layout(),
    ]
    verify, _calls = _assembly(lambda values: values.get(SCOPE_PATH) == "outer")

    with pytest.raises(LabelsDoNotFitError) as raised:
        check_circular_label_fit(layouts, verify)

    rows = {row.kind: row for row in raised.value.report.rows}
    # Nothing was drawn, so the definition row is not offered (D-12).
    assert list(rows) == [
        kind for kind in circular_fit_module.ROW_KINDS if kind not in {"separate_features", "smaller_definition_font", "narrow_tracks"}
    ]
    assert rows["labels_out"].status == "found"
    # The outer levers do not move embedded labels either (OV-725).
    assert rows["outward_vertical"].reason == "Does not move inner, embedded, or radial labels"
    assert "--labels out   (all labels outside)" in str(raised.value)


def test_radial_labels_without_room_stop_and_offer_the_font_row_only_without_inner_labels() -> None:
    """RADIAL-JOIN: the record is refused; labels_out applies only when inner labels exist."""

    verify, _calls = _assembly(lambda values: False)

    with pytest.raises(LabelsDoNotFitError) as raised:
        check_circular_label_fit([_inner_layout(placement="radial", refused=True, inner=0)], verify)

    report = raised.value.report
    (record,) = report.records
    assert record.refused and record.label_kind == "circular_outer"
    assert [row.kind for row in report.rows] == ["smaller_font"]
    assert str(raised.value).startswith(
        "P1: 80 radial outer labels do not fit at 14 pt (they need more room than the circular tracks leave). "
        "Nothing was drawn.\n"
    )


def test_the_web_error_of_a_refused_record_states_the_refusal_not_its_counts() -> None:
    """OV-1194: a refused record's counts are of the layout that was refused, so the
    Web error sentence names the refusal the dialog states, not "Crossing leader pairs: 1"."""

    refused = LabelFitRecord("NC_001879.2", "circular_inner", 84, 8.0, LabelFit(0, 1, 0), refused=True)
    error = LabelsDoNotFitError("x", report=LabelFitReport((refused,), (), drawn=False))

    assert error.diagnostic == {
        "code": "LABEL_LAYOUT",
        "labelFitCounts": "Radial labels need more room than the circular tracks leave.",
    }


def test_one_budget_covers_the_switches_then_the_font_scan(monkeypatch: pytest.MonkeyPatch) -> None:
    """OD-F15: 30 s for inner labels when nothing is drawn; a verification starts only
    if one more assembly ends before the deadline; the error computes the rest."""

    clock = [0.0]
    monkeypatch.setattr(circular_fit_module.time, "monotonic", lambda: clock[0])
    verify, calls = _assembly(_p1_inner)

    def slow(values: Mapping[str, object]) -> LabelFitVerification:
        clock[0] += 12.0
        return verify(values)

    with pytest.raises(LabelsDoNotFitError) as raised:
        check_circular_label_fit([_inner_layout()], slow, assembly_seconds=12.0)

    # labels_out starts at 0 s and radial at 12 s; the font scan would end after 36 s.
    assert calls == [{SCOPE_PATH: "outer"}, {PLACEMENT_PATH: "radial"}]
    rows = {row.kind: row.status for row in raised.value.report.rows}
    assert rows == {"smaller_font": "not_computed", "labels_out": "not_found", "radial": "not_found"}
    assert str(raised.value).splitlines()[-1] == "Not searched within 30 s: a smaller font."
    row = label_fit_row(raised.value, "smaller_font")
    assert row is not None and row.values == {path: 12.0 for path in FONT_SIZE_PATHS}


def test_the_outer_error_path_computes_every_row_past_the_budget(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(circular_fit_module, "SWITCH_SEARCH_SECONDS", 0.0)
    verify, _calls = _verifier(_p1_like)

    with pytest.raises(LabelsDoNotFitError) as raised:
        check_circular_label_fit([_layout("P1", _p1_like)], verify, assembly_seconds=60.0)

    assert "not_computed" not in {row.status for row in raised.value.report.rows}


def test_overlapping_embedded_labels_offer_to_separate_the_features() -> None:
    """LP-D52 E-ROW: the switch is verified first and found only when the figure then fits."""

    layout = dataclasses.replace(
        _inner_layout(inner=0),
        label_counts={"circular_outer": 80, "circular_embedded": 3},
        kind_fits={"circular_embedded": LabelFit(1, 0)},
    )
    verify, calls = _assembly(lambda values: values.get(SEPARATE_FEATURES_PATH) is True)

    with pytest.raises(LabelsDoNotFitError) as raised:
        check_circular_label_fit([layout], verify)

    assert calls[0] == {SEPARATE_FEATURES_PATH: True}
    rows = {row.kind: row for row in raised.value.report.rows}
    assert list(rows) == ["smaller_font", "radial", "separate_features"]
    assert rows["separate_features"].values == {SEPARATE_FEATURES_PATH: True}
    assert "  --resolve_overlaps   (separate overlapping features)" in str(raised.value)

    separated = dataclasses.replace(layout, features_separated=True)
    with pytest.raises(LabelsDoNotFitError) as raised:
        check_circular_label_fit([separated], verify)
    assert "separate_features" not in {row.kind for row in raised.value.report.rows}


@pytest.mark.circular
def test_separated_features_that_leave_no_track_room_are_a_row_not_found(tmp_path: Path) -> None:
    """OV-1266: lambda at 12 pt with gene names first overprints two embedded labels; separating the features
    needs six lanes the inside tracks cannot hold, so that row is not found and the labels ask (E-ROW). The
    definition no longer reserves the centre on its own (DF-D02), so the run reserves it with center_reserved_radius."""

    priority = tmp_path / "priority.tsv"
    priority.write_text("CDS\tgene,product,locus_tag\n", encoding="utf-8")
    args = [
        "--gbk", str(TEST_INPUTS / "NC_001416.gb"),
        "--separate_strands", "--track_type", "tuckin", "-l", "left",
        "--labels", "out", "--qualifier_priority", str(priority), "--label_font_size", "12",
        "--center_reserved_radius", "150",
        "-o", str(tmp_path / "out"), "-f", "svg",
    ]

    with pytest.raises(LabelsDoNotFitError) as raised:
        circular_main(args)

    rows = {row.kind: row for row in raised.value.report.rows}
    assert rows["separate_features"].status == "not_found"
    assert rows["separate_features"].reason == "The tracks do not fit with this setting"


def _fit_assembly(fit_for: Callable[[Mapping[str, object]], LabelFit]):
    """A full assembly whose labels' fit depends on the row's settings."""

    calls: list[dict[str, object]] = []

    def verify(values: Mapping[str, object]) -> LabelFitVerification:
        calls.append(dict(values))
        fit = fit_for(values)
        return LabelFitVerification(label_fit_outcome([fit]), (1000.0, 900.0), fit)

    return verify, calls


def _embedded_overlap_layout() -> CircularLabelLayout:
    return dataclasses.replace(
        _inner_layout(inner=0),
        label_counts={"circular_outer": 80, "circular_embedded": 3},
        kind_fits={"circular_embedded": LabelFit(1, 0)},
    )


def test_a_row_that_leaves_only_crossings_is_offered_after_the_rows_that_fit() -> None:
    """LP-D57 (OV-724): a row that removes every overlap but leaves crossing leaders is
    offered, says what remains, and comes after the rows that fit."""

    def fit_for(values: Mapping[str, object]) -> LabelFit:
        if values.get(SEPARATE_FEATURES_PATH) is True:
            return LabelFit(0, 14)
        font = values.get(FONT_SIZE_PATHS[0])
        return FITS if isinstance(font, float) and font <= 12.0 else LabelFit(1, 14)

    verify, _calls = _fit_assembly(fit_for)
    with pytest.raises(LabelsDoNotFitError) as raised:
        check_circular_label_fit([_embedded_overlap_layout()], verify)

    rows = raised.value.report.rows
    assert [(row.kind, row.status) for row in rows] == [
        ("smaller_font", "found"),
        ("radial", "not_found"),
        ("separate_features", "found_with_notice"),
    ]
    separate = rows[2]
    assert separate.values == {SEPARATE_FEATURES_PATH: True}
    assert separate.canvas_size == (1000.0, 900.0)
    assert separate.remaining == LabelFit(0, 14)
    assert row_note(separate) == "14 leader pairs still cross"
    assert row_note(rows[0]) is None
    assert str(raised.value).splitlines()[1:4] == [
        "Run again with one of:",
        "  --label_font_size 12   (smaller font)",
        "  --resolve_overlaps     (separate overlapping features; 14 leader pairs still cross)",
    ]


def test_the_notice_row_words_crossings_and_the_definition_text() -> None:
    def row(remaining: LabelFit) -> LabelFitRow:
        return LabelFitRow("radial", "found_with_notice", {PLACEMENT_PATH: "radial"}, (1.0, 1.0), remaining=remaining)

    assert row_note(row(LabelFit(0, 1))) == "1 leader pair still crosses"
    assert row_note(row(LabelFit(0, 0, 2))) == "2 labels still cover the definition text"
    assert row_note(row(LabelFit(0, 3, 1))) == "3 leader pairs still cross and 1 label still covers the definition text"


def test_each_lever_offers_a_value_that_leaves_only_crossings_when_none_fits() -> None:
    """LP-D57: every row kind, the font scan and the outer levers too; a value that fits
    stays preferred, and a value that leaves overlaps is still not found."""

    def inner_fit(values: Mapping[str, object]) -> LabelFit:
        font = values.get(FONT_SIZE_PATHS[0])
        if values.get(PLACEMENT_PATH) == "radial":
            return LabelFit(2, 0)
        return LabelFit(0, 0, 1) if isinstance(font, float) and font <= 13.0 else LabelFit(1, 0)

    verify, _calls = _fit_assembly(inner_fit)
    with pytest.raises(LabelsDoNotFitError) as raised:
        check_circular_label_fit([_inner_layout()], verify)
    rows = {row.kind: row for row in raised.value.report.rows}
    assert rows["smaller_font"].status == "found_with_notice"
    assert rows["smaller_font"].values == {path: 13.0 for path in FONT_SIZE_PATHS}
    assert rows["smaller_font"].remaining == LabelFit(0, 0, 1)
    assert rows["radial"].status == "not_found"

    def outer_fit(values: Mapping[str, object]) -> LabelFit:
        y_offset = values.get(Y_OFFSET_PATH, 1.0)
        return LabelFit(0, 2) if isinstance(y_offset, float) and y_offset >= 1.3 - 1e-9 else OVERLAPS

    layout = dataclasses.replace(_layout("P1", _p1_like), probe=lambda values: outer_fit(values))
    verify, _calls = _fit_assembly(outer_fit)
    with pytest.raises(LabelsDoNotFitError) as raised:
        check_circular_label_fit([layout], verify)
    rows = {row.kind: row for row in raised.value.report.rows}
    assert rows["outward_vertical"].status == "found_with_notice"
    assert rows["outward_vertical"].values == {Y_OFFSET_PATH: 1.3}
    assert rows["outward_vertical"].remaining == LabelFit(0, 2)
    assert rows["outward_horizontal"].status == "not_found"
    assert rows["outward_horizontal"].reason == "Does not fit up to X offset 2"
    statuses = [row.status for row in raised.value.report.rows]
    assert statuses == sorted(statuses, key=lambda status: status == "found_with_notice")


def test_a_drawn_figure_offers_only_rows_that_fit() -> None:
    """A figure drawn with crossings already has the notice: a row that leaves crossings is no fix."""

    layout = dataclasses.replace(_inner_layout(), kind_fits={"circular_inner": LabelFit(0, 2)})
    verify, _calls = _fit_assembly(lambda values: LabelFit(0, 1))
    with pytest.warns(LabelLeadersCrossWarning):
        checked = check_circular_label_fit([layout], verify, warn_search_seconds=60.0)
    assert checked is not None
    rows = {row.kind: row for row in checked[0].rows}
    assert rows["labels_out"].status == "not_found"
    assert rows["labels_out"].reason == "Labels still do not fit (0 label pairs overlap, 1 leader pair crosses)"
    assert checked[1].row("labels_out").status == "not_found"


def test_labels_count_once_by_kind_and_leaders_count_against_the_definition_text() -> None:
    def box(min_x: float, min_y: float, max_x: float, max_y: float) -> tuple[tuple[float, float], ...]:
        return ((min_x, min_y), (max_x, min_y), (max_x, max_y), (min_x, max_y))

    shapes = [
        LabelShape("circular_outer", (box(0, 0, 10, 5),), ((5.0, 6.0), (5.0, 20.0))),
        LabelShape("circular_inner", (box(8, 0, 18, 5),), ((12.0, 6.0), (12.0, 20.0))),
        LabelShape("circular_inner", (box(100, 100, 110, 105),), ((105.0, 99.0), (105.0, 90.0))),
        LabelShape("circular_embedded", (box(100, 100, 120, 105),)),
    ]
    definition = [Aabb(99.0, 99.0, 130.0, 104.0)]

    fits = label_fits_by_kind(shapes, definition_text=definition)

    assert fits == {
        "circular_inner": LabelFit(2, 0, 1),
        "circular_outer": LabelFit(0, 0, 0),
        "circular_embedded": LabelFit(0, 0, 0),
    }


def test_linear_labels_on_stacked_features_stop_with_the_separate_features_row() -> None:
    def label(start: float, middle_y: float, *, embedded: bool) -> dict:
        return {
            "start": start, "middle": start + 20.0, "middle_x": start + 20.0, "width_px": 40.0,
            "height_px": 10.0, "middle_y": middle_y, "feature_middle_y": 50.0,
            "is_embedded": embedded, "leader_line": False, "font_size": 8.0,
            "rotation_deg": 0.0, "text_anchor": "middle",
        }

    stacked = [label(0.0, 50.0, embedded=True), label(10.0, 50.0, embedded=True), label(80.0, 20.0, embedded=False)]
    label_fit = linear_label_fit([("rec", stacked)], 8.0)
    assert label_fit.embedded_overlap and label_fit.fit == LabelFit(1, 0)
    calls: list[Mapping[str, object]] = []

    def verify(values: Mapping[str, object]) -> LabelFitVerification:
        calls.append(dict(values))
        return LabelFitVerification("fits", (500.0, 200.0))

    with pytest.raises(LabelsDoNotFitError) as raised:
        check_linear_label_fit(label_fit, verify, features_separated=False)

    assert str(raised.value).startswith("rec: 3 Linear labels do not fit at 8 pt (1 label pair overlaps")
    (row,) = raised.value.report.rows
    assert (row.kind, row.status, row.values) == ("separate_features", "found", {SEPARATE_FEATURES_PATH: True})
    assert calls == [{SEPARATE_FEATURES_PATH: True}]

    # LP-D57: separated features whose leaders then cross are offered with the notice.
    def crossing(values: Mapping[str, object]) -> LabelFitVerification:
        return LabelFitVerification("warn", (500.0, 200.0), LabelFit(0, 2))

    with pytest.raises(LabelsDoNotFitError) as raised:
        check_linear_label_fit(label_fit, crossing, features_separated=False)
    (row,) = raised.value.report.rows
    assert (row.status, row.values, row.remaining) == (
        "found_with_notice", {SEPARATE_FEATURES_PATH: True}, LabelFit(0, 2)
    )

    with pytest.raises(LabelsDoNotFitError) as raised:
        check_linear_label_fit(label_fit, verify, features_separated=True)
    assert raised.value.report.rows == ()

    placed = [label(0.0, 20.0, embedded=False), label(10.0, 22.0, embedded=False)]
    with pytest.raises(LabelsDoNotFitError) as raised:
        check_linear_label_fit(linear_label_fit([("rec", placed)], 8.0), verify, features_separated=False)
    assert raised.value.report.rows == ()
    check_linear_label_fit(
        linear_label_fit([("rec", [label(0.0, 20.0, embedded=False)])], 8.0), verify, features_separated=False
    )


def _no_values(values: Mapping[str, object]) -> LabelFitVerification:
    raise AssertionError(f"no config row expected: {dict(values)}")


def test_region_labels_over_other_labels_warn_with_no_fix_and_never_stop() -> None:
    """LP-D68, LP-D82: a figure with region labels over other labels or lines is drawn and warns
    with the count; no setting is offered (the font and the line count never change)."""

    regions = RecordRegionLabels((10.0, 12.0, 14.0), overlaps=7)
    label_fit = linear_label_fit([("rec", [])], 8.0, [regions])
    assert label_fit.outcome == "warn"
    (record,) = label_fit.records
    assert (record.label_kind, record.label_count, record.font_size, record.fit) == (
        "linear_region", 3, 14.0, LabelFit(0, 0, region_overlaps=7)
    )
    with pytest.warns(LabelLeadersCrossWarning) as caught:
        report = check_linear_label_fit(label_fit, _no_values, features_separated=False)
    assert report is not None and report.drawn and report.region_overlaps == 7 and report.rows == ()
    assert label_fit_notice(report) == ("Region labels overlap other labels or lines (7).",)
    assert str(caught[0].message) == "rec: 3 region labels do not fit at 14 pt (7 region label overlaps). Drawn anyway."


def test_feature_label_rows_are_kept_only_when_the_region_label_that_asked_for_them_clears() -> None:
    """LP-D77, LP-D82 RGN-YIELD-ROW: region label 0 asks for a new row; a second figure that clears
    region label 1 instead, with the same total, is not kept. One that clears label 0 is."""

    rows = FeatureLabelRows((FeatureLabelRow(0, 3, 4),), overlaps=1, region_labels=((0, 0),))
    before = linear_label_fit([("rec", [])], 8.0, [RecordRegionLabels((10.0, 10.0), 2, conflicted=(0, 1))], rows)
    still_asking = linear_label_fit([("rec", [])], 8.0, [RecordRegionLabels((10.0, 10.0), 1, conflicted=(0,))])
    cleared = linear_label_fit([("rec", [])], 8.0, [RecordRegionLabels((10.0, 10.0), 1, conflicted=(1,))])
    assert not keeps_feature_label_rows(before, still_asking)
    assert keeps_feature_label_rows(before, cleared)


def test_overlapping_feature_labels_stop_and_name_the_region_count_too() -> None:
    def label(start: float) -> dict:
        return {
            "start": start, "middle": start + 20.0, "middle_x": start + 20.0, "width_px": 40.0, "height_px": 10.0,
            "middle_y": 20.0, "feature_middle_y": 50.0, "is_embedded": False, "leader_line": False,
            "font_size": 8.0, "rotation_deg": 0.0, "text_anchor": "middle",
        }

    label_fit = linear_label_fit([("rec", [label(0.0), label(10.0)])], 8.0, [RecordRegionLabels((10.0,), 2)])
    assert label_fit.outcome == "stop"
    with pytest.raises(LabelsDoNotFitError) as raised:
        check_linear_label_fit(label_fit, _no_values, features_separated=False)
    assert str(raised.value) == (
        "Labels do not fit: rec: 2 Linear labels at 8 pt (1 label pair overlaps, 0 leader pairs cross); "
        "rec: 1 region label at 10 pt (2 region label overlaps). Nothing was drawn."
    )


# --- engine and CLI -------------------------------------------------------------------------------

LAMBDA_ARGS = [
    "--gbk", str(TEST_INPUTS / "NC_001416.gb"),
    "--separate_strands", "--track_type", "tuckin", "--labels", "out", "--legend", "right",
]


@pytest.mark.circular
@pytest.mark.usefixtures("todays_label_placement_only")
def test_cli_draws_crossing_leaders_and_prints_the_warning(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """LP-D46 CROSS-WARN: lambda at 15 pt has 2 crossing leader pairs and no overlaps (with the repair and packer off).

    14 pt had them while labels were laid out in point-mode boxes; in their drawn boxes 14 pt fits (OV-165)."""

    from gbdraw import cli

    monkeypatch.setattr(circular_repair_module, "MAX_DEFECTS", 0)
    # The warning search stops after 5 s of wall clock (LP-D48), so a slow runner prints
    # "Not searched within 5 s" instead of the rows. A clock that stands still searches every row.
    monkeypatch.setattr(circular_fit_module.time, "monotonic", lambda: 0.0)
    output = tmp_path / "lambda"
    monkeypatch.setattr(
        sys, "argv", ["gbdraw", "circular", *LAMBDA_ARGS, "--label_font_size", "15", "-o", str(output), "-f", "svg"]
    )
    cli.main()

    assert output.with_suffix(".svg").is_file()
    stderr = capsys.readouterr().err
    assert stderr.startswith(
        "WARNING: NC_001416.1: 71 outer labels do not fit at 15 pt (0 label pairs overlap, 2 leader pairs cross). "
        "Drawn anyway.\nRun again with one of:\n"
    )
    assert "--outer_label_y_radius_offset 1.025" in stderr


@pytest.mark.circular
@pytest.mark.usefixtures("todays_label_placement_only")
def test_drawing_with_crossing_leaders_carries_the_report_and_computes_rows(monkeypatch: pytest.MonkeyPatch) -> None:
    """The Web reads the report from the success result and asks for one row at a time.

    The local repair (LP-D57) and the packer clear lambda's 2 crossings; both are off here, as for crossings
    they cannot clear.
    Lambda at 15 pt (14 pt while labels were laid out in point-mode boxes; OV-165).
    """

    from gbdraw.api.diagram import build_circular_diagram, label_fit_report, label_fit_row
    from gbdraw.api.options import CircularDiagramOptions
    from gbdraw.config.modify import with_config_overrides
    from gbdraw.config.toml import load_config_toml

    monkeypatch.setattr(circular_repair_module, "MAX_DEFECTS", 0)
    cfg = with_config_overrides(
        GbdrawConfig.from_dict(load_config_toml("gbdraw.data", "config.toml")),
        {
            "labels.filtering.blacklist_keywords": [],
            "labels.circular.scope": "outer",
            "canvas.strandedness": True,
            "canvas.circular.track_type": "tuckin",
            "labels.font_size.short": 15.0,
            "labels.font_size.long": 15.0,
        },
    )
    record = SeqIO.read(str(TEST_INPUTS / "NC_001416.gb"), "genbank")

    with pytest.warns(LabelLeadersCrossWarning):
        drawing = build_circular_diagram(record, options=CircularDiagramOptions(config=cfg))

    report = label_fit_report(drawing)
    assert report is not None and report.drawn
    assert report.records[0].label_kind == "circular_outer"
    assert report.records[0].fit == LabelFit(0, 2)
    vertical = label_fit_row(drawing, "outward_vertical")
    assert vertical is not None and vertical.status == "found"
    assert vertical.values == {Y_OFFSET_PATH: 1.025}


def _drawn_label_digest(monkeypatch: pytest.MonkeyPatch, tmp_path: Path, args: list[str]) -> tuple[int, str]:
    drawn: list[list[dict]] = []
    original = circular_assemble_module.add_labels_group_on_canvas

    def capture(canvas, *capture_args, **kwargs):
        if kwargs.get("phase") == "leaders":
            drawn.append(
                [
                    label
                    for label in kwargs["precalculated_labels"]
                    if not label.get("is_embedded") and not label.get("is_inner")
                ]
            )
        return original(canvas, *capture_args, **kwargs)

    monkeypatch.setattr(circular_assemble_module, "add_labels_group_on_canvas", capture)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", LabelLeadersCrossWarning)
        circular_main([*args, "-o", str(tmp_path / "out"), "-f", "svg"])
    assert_no_label_warning(caught)
    keys = ("start_x", "start_y", "middle_x", "middle_y", "leader_start_x", "leader_start_y")
    rows = [[label["label_text"], *(round(float(label.get(key, 0.0)), 4) for key in keys)] for label in drawn[-1]]
    return len(rows), hashlib.sha256(json.dumps(rows).encode()).hexdigest()


@pytest.mark.circular
def test_dense_map_that_fits_with_legacy_rows_draws_as_before(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """LP-D45: plasmid NZ_CP125877.1 at 8 pt overlaps with 3 px rows; it draws as on 0.13 (0.1 px rows)."""

    records = SeqIO.parse(str(TEST_INPUTS / "GCF_030060435.1_ASM3006043v1_genomic.gbff"), "genbank")
    record = next(record for record in records if record.id == "NZ_CP125877.1")
    gbk = tmp_path / "NZ_CP125877.1.gb"
    SeqIO.write(record, str(gbk), "genbank")

    digest = _drawn_label_digest(
        monkeypatch,
        tmp_path,
        [
            "--gbk", str(gbk), "--separate_strands", "--track_type", "tuckin", "--labels", "out",
            "--legend", "right", "--label_font_size", "8",
        ],
    )

    # Digest of the same run with labels laid out in their drawn boxes (OV-165; the 179 labels keep their order,
    # with no drawn overlap and no leader through another label's letters; OV-1823 moves two labels 1-2 px out).
    # Before, the digest of origin/dev 3d748b78 (before OV-479) was 0572c062...
    assert digest == (179, "beb197dace5ea1dc78a1811febe973c733c2db35661d82f41a0f5ce32b00cd94")


@pytest.mark.circular
def test_dense_map_whose_spaced_rows_cross_draws_as_before(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """LP-D45: E. coli RNA labels at 8 pt fit with 3 px rows but 1 leader pair crosses; drawn as on 0.13."""

    whitelist = tmp_path / "rna_whitelist.tsv"
    whitelist.write_text("tRNA\tproduct\ttRNA\nrRNA\tproduct\tribosomal RNA\n")

    digest = _drawn_label_digest(
        monkeypatch,
        tmp_path,
        [
            "--gbk", str(TEST_INPUTS / "MG1655.gbk"), "--separate_strands", "--track_type", "tuckin",
            "--labels", "out", "--legend", "right", "--label_font_size", "8",
            "--label_whitelist", str(whitelist),
        ],
    )

    # Digest of the same run with labels laid out in their drawn boxes (OV-165; the 108 labels keep their order,
    # with no drawn overlap). OV-1823 steps six labels out, so 1 leader runs through another label's letters
    # (origin/dev f5992e64a: 3). Before, the digest of origin/dev 3d748b78 (before OV-479) was cf67ec84...
    assert digest == (108, "b5b4275737da5e3af8fa8a8a48964c13769ad0950df6d1a17e406bb74d86bbbe")
