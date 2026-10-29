"""Circular label capacity: stop at once when horizontal labels clearly cannot fit, warn when they may not (PD-OI-097)."""

from __future__ import annotations

import sys
import time
import warnings
from pathlib import Path

import pytest

import gbdraw.labels.circular as circular_module
from gbdraw.circular import circular_main
from gbdraw.exceptions import GbdrawError, LabelsDoNotFitError, LabelsMayNotFitWarning
from gbdraw.labels.circular_capacity import (
    CAPACITY_LEVER_RATIO,
    CAPACITY_STOP_RATIO,
    CAPACITY_WARN_RATIO,
    CircularLabelCapacity,
    LabelArc,
    circular_label_capacity,
    label_capacity_unchecked,
)
from gbdraw.labels.circular_fit import CircularLabelLayout, check_circular_label_fit
from gbdraw.labels.circular_types import CircularLabel
from gbdraw.labels.fit import LabelCapacityRecord, LabelFit, format_label_fit_report

TEST_INPUTS = Path(__file__).parent / "test_inputs"
AP027078_ARGS = ["--gbk", str(TEST_INPUTS / "AP027078.gb"), "--labels", "out"]
T5_ARGS = [
    "--gbk", str(TEST_INPUTS / "NC_005859.1.gb"),
    "--separate_strands", "--track_type", "tuckin", "--labels", "out", "--legend", "right",
]


class _Placed(GbdrawError):
    """Placement was reached."""


def _labels(count: int, *, width: float, height: float) -> list[CircularLabel]:
    return [{"width_px": width, "height_px": height} for _ in range(count)]  # type: ignore[typeddict-item]


def _capacity(labels: list[CircularLabel], arc: LabelArc, *, font: float = 8.0) -> CircularLabelCapacity:
    (capacity,) = circular_label_capacity(
        labels, [], outer_arc=arc, inner_arc=arc, outer_offsets=(1.0, 1.0), font_size=font
    )
    return capacity


# --- the model -------------------------------------------------------------------------------------


def test_labels_that_stack_on_the_ellipse_need_at_most_its_room() -> None:
    """Boxes too wide to slide past each other stack along the vertical travel ``4 * y_radius``."""

    # 4 * 100 px of room; 40 boxes of 10 px (+0.1 margin) less the two pole crossings fill 0.96 of it.
    arc = LabelArc(x_radius=1.0, y_radius=100.0, center_x=0.0)
    capacity = _capacity(_labels(40, width=500.0, height=10.0), arc)
    assert 0.9 < capacity.record.ratio <= 1.0
    # A wider ellipse lets boxes slide sideways and a taller one has more room: both lower the ratio.
    assert _capacity(_labels(40, width=500.0, height=10.0), LabelArc(400.0, 100.0, 0.0)).record.ratio < capacity.record.ratio
    assert _capacity(_labels(40, width=500.0, height=10.0), LabelArc(1.0, 200.0, 0.0)).record.ratio < capacity.record.ratio


def test_the_levels_read_the_ratio_and_the_lever() -> None:
    arc = LabelArc(x_radius=1.0, y_radius=100.0, center_x=0.0)
    # 0.8 of the room or more warns before placement.
    assert _capacity(_labels(30, width=500.0, height=10.0), arc).level == "fits"
    assert _capacity(_labels(36, width=500.0, height=10.0), arc).level == "near"
    # Over 1.5 times, but 6 pt or the outward rows make room: still placed (it warns).
    over = _capacity(_labels(70, width=500.0, height=10.0), arc, font=14.0)
    assert over.record.ratio > CAPACITY_STOP_RATIO and over.lever <= CAPACITY_LEVER_RATIO
    assert over.level == "near"
    # Over 1.5 times and no single change makes room: stop.
    crowded = _capacity(_labels(200, width=500.0, height=10.0), arc, font=8.0)
    assert crowded.lever > CAPACITY_LEVER_RATIO and crowded.level == "over"
    assert CAPACITY_WARN_RATIO == 0.8 and CAPACITY_STOP_RATIO == 1.5 and CAPACITY_LEVER_RATIO == 1.0


def test_inner_labels_read_putting_them_outside_as_their_lever() -> None:
    outer_arc = LabelArc(x_radius=1.0, y_radius=1000.0, center_x=0.0)
    inner_arc = LabelArc(x_radius=1.0, y_radius=50.0, center_x=0.0)
    _outer, inner = circular_label_capacity(
        _labels(10, width=500.0, height=10.0),
        _labels(60, width=500.0, height=10.0),
        outer_arc=outer_arc,
        inner_arc=inner_arc,
        outer_offsets=(1.0, 1.0),
        font_size=6.0,
    )
    assert inner.label_kind == "circular_inner" and inner.record.ratio > CAPACITY_STOP_RATIO
    # At the 6 pt floor the font does not help, but all labels outside fit the outer ellipse.
    assert inner.lever < 1.0 and inner.level == "near"


def test_the_font_remedy_is_left_out_at_the_font_row_floor() -> None:
    arc = LabelArc(x_radius=1.0, y_radius=10.0, center_x=0.0)
    assert _capacity(_labels(5, width=5.0, height=5.0), arc, font=8.0).record.remedies == (
        "filter", "label_text", "font", "linear"
    )
    assert _capacity(_labels(5, width=5.0, height=5.0), arc, font=6.0).record.remedies == (
        "filter", "label_text", "linear"
    )


def test_a_capacity_stop_lists_every_failing_record_and_searches_no_row() -> None:
    """The record over capacity and a record whose labels overlap are both named; no row is verified."""

    capacity = LabelCapacityRecord(562, 8.0, 2.49, ("filter", "label_text", "font", "linear"))

    common = {"placement": "horizontal", "font_size": 8.0, "x_radius_offset": 1.0, "y_radius_offset": 1.0,
              "label_spacing": 3.0}
    over = CircularLabelLayout(
        record_id="A.1", record_index=0, label_counts={"circular_outer": 562}, kind_fits={},
        capacity={"circular_outer": capacity}, **common,  # type: ignore[arg-type]
    )
    overlapping = CircularLabelLayout(
        record_id="B.1", record_index=1, label_counts={"circular_outer": 40},
        kind_fits={"circular_outer": LabelFit(3, 1)}, **common,  # type: ignore[arg-type]
    )

    def verify(values: object) -> object:
        raise AssertionError("a capacity stop verifies no row")

    with pytest.raises(LabelsDoNotFitError) as caught:
        check_circular_label_fit([over, overlapping], verify)  # type: ignore[arg-type]

    report = caught.value.report
    assert caught.value.row_search is None and report.rows == () and not report.drawn
    assert [(record.record_id, record.capacity) for record in report.records] == [("A.1", capacity), ("B.1", None)]
    assert str(caught.value).splitlines() == [
        "A.1: 562 outer labels at 8 pt need about 2.5 times the room around the map, so they were not placed.",
        "B.1: 40 outer labels do not fit at 8 pt (3 label pairs overlap, 1 leader pair crosses). Nothing was drawn.",
        "To make room:",
        "  - label fewer features with --label_blacklist or --label_whitelist;",
        "  - shorten the label text with --qualifier_priority (for example gene names only);",
        "  - use a smaller --label_font_size;",
        "  - or draw a Linear map (gbdraw linear).",
    ]
    assert format_label_fit_report(report) == str(caught.value)


# --- the figure ------------------------------------------------------------------------------------


def _placement_raises(monkeypatch: pytest.MonkeyPatch) -> list[float]:
    """Make horizontal placement raise ``_Placed``; the returned list gets the time it was reached."""

    reached: list[float] = []

    def placed(*args: object, **kwargs: object) -> object:
        reached.append(time.monotonic())
        raise _Placed("placement reached")

    monkeypatch.setattr(circular_module, "rearrange_labels_fc", placed)
    return reached


@pytest.mark.circular
def test_ap027078_stops_at_once_without_placing_or_searching(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """CAPACITY-STOP: 562 outer labels at 8 pt need about 2.5 times the room; nothing is placed and no row searched."""

    reached = _placement_raises(monkeypatch)
    started = time.monotonic()
    with pytest.raises(LabelsDoNotFitError) as caught:
        circular_main([*AP027078_ARGS, "-o", str(tmp_path / "ap"), "-f", "svg"])

    assert not reached
    assert time.monotonic() - started < 20.0
    error = caught.value
    assert error.row_search is None and error.report.rows == ()
    (record,) = error.report.records
    assert record.capacity is not None
    assert (record.record_id, record.label_kind, record.capacity.label_count, record.capacity.font_size) == (
        "AP027078.1", "circular_outer", 562, 8.0
    )
    assert 2.4 < record.capacity.ratio < 2.6
    # OV-2179: the room as a value the Web words in the UI language, not an English sentence.
    assert error.diagnostic == {"code": "LABEL_LAYOUT", "reason": "LABEL_CAPACITY", "labelCapacityOuter": 2.5}


def _capacities(monkeypatch: pytest.MonkeyPatch) -> list[tuple[CircularLabelCapacity, ...]]:
    seen: list[tuple[CircularLabelCapacity, ...]] = []
    check = circular_module.check_label_capacity

    def spy(record_id: str, capacities: tuple[CircularLabelCapacity, ...]) -> None:
        seen.append(tuple(capacities))
        check(record_id, capacities)

    monkeypatch.setattr(circular_module, "check_label_capacity", spy)
    return seen


@pytest.mark.circular
@pytest.mark.parametrize(("font", "ratio"), [(14, (1.3, 1.4)), (16, (1.5, 1.65))])
def test_t5_warns_before_placement_and_is_then_placed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, font: int, ratio: tuple[float, float]
) -> None:
    """CAPACITY-WARN: T5 at 14 pt (ratio 1.36) and at 16 pt (ratio 1.58, over 1.5, but 6 pt makes room: the lever
    keeps the rows found today) warn first, then place as before."""

    reached = _placement_raises(monkeypatch)
    seen = _capacities(monkeypatch)
    warned: list[float] = []
    show = warnings.showwarning

    def record(message, category, *args, **kwargs):  # type: ignore[no-untyped-def]
        if issubclass(category, LabelsMayNotFitWarning):
            warned.append(time.monotonic())
        return show(message, category, *args, **kwargs)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", LabelsMayNotFitWarning)
        warnings.showwarning = record
        with pytest.raises(_Placed):
            circular_main([*T5_ARGS, "--label_font_size", str(font), "-o", str(tmp_path / "t5"), "-f", "svg"])

    (capacity,) = seen[0]
    assert ratio[0] < capacity.record.ratio < ratio[1] and capacity.lever <= CAPACITY_LEVER_RATIO
    (warning,) = [item.message for item in caught if isinstance(item.message, LabelsMayNotFitWarning)]
    assert warning.record_id == "NC_005859.1" and warning.capacity == capacity.record
    assert str(warning) == (
        f"NC_005859.1: 195 outer labels at {font} pt need about {capacity.record.ratio:.1f} times the room around the "
        "map, so they may not fit. Placing them can take minutes."
    )
    assert reached


@pytest.mark.circular
def test_the_check_and_the_placement_read_one_ellipse(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """``label_arc`` is the one owner of the ellipse: the check's outer ellipse is the one placement uses, with the
    same outer label arena (a drawn CLI figure always has one)."""

    arcs: list[tuple[tuple[object, ...], LabelArc]] = []
    owner = circular_module.label_arc

    def spy(*args: object) -> LabelArc:
        arc = owner(*args)  # type: ignore[arg-type]
        arcs.append((args, arc))
        if len(arcs) == 3:
            raise _Placed("placement reached")
        return arc

    monkeypatch.setattr(circular_module, "label_arc", spy)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", LabelsMayNotFitWarning)
        with pytest.raises(_Placed):
            circular_main([*T5_ARGS, "--label_font_size", "8", "-o", str(tmp_path / "t5"), "-f", "svg"])

    (check_outer_args, check_outer), (_inner_args, _inner), (placed_args, placed) = arcs
    assert check_outer_args[4] is True and placed_args[4] is True
    assert check_outer_args[5] is not None and placed_args[5] == check_outer_args[5]
    assert placed == check_outer


@pytest.mark.circular
def test_a_row_verification_places_labels_without_the_check(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    reached = _placement_raises(monkeypatch)
    with label_capacity_unchecked(), pytest.raises(_Placed):
        circular_main([*AP027078_ARGS, "-o", str(tmp_path / "ap"), "-f", "svg"])
    assert reached


def _crowded_map(path: Path) -> Path:
    """A 20 kb map with 100 CDS labels: at 16 pt they need about the room around the map (CAPACITY-WARN)."""

    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord

    record = SeqRecord(Seq("ATGC" * 5000), id="SYN1", name="SYN1", description="crowded map")
    record.annotations.update(molecule_type="DNA", topology="circular")
    for index in range(100):
        start = index * 200
        record.features.append(SeqFeature(
            FeatureLocation(start, start + 150, strand=1 if index % 2 else -1), type="CDS",
            qualifiers={"product": [f"conserved hypothetical protein {index:03d}"]},
        ))
    SeqIO.write(record, path, "genbank")
    return path


@pytest.mark.circular
def test_fit_probes_and_row_verification_never_run_the_check(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A crowded map warns once, then stops with a row found: the row search's probes and the verification of the
    found row place the labels again, and none of them reads or warns the capacity."""

    import gbdraw.labels.circular_fit as fit_module

    # One row keeps the search short; it still probes, then verifies the row it finds.
    monkeypatch.setattr(fit_module, "_OUTER_ROWS", frozenset({"smaller_font"}))
    seen = _capacities(monkeypatch)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", LabelsMayNotFitWarning)
        with pytest.raises(LabelsDoNotFitError) as stopped:
            circular_main([
                "--gbk", str(_crowded_map(tmp_path / "crowded.gb")), "--labels", "out", "--label_font_size", "16",
                "-o", str(tmp_path / "crowded"), "-f", "svg",
            ])

    assert [[capacity.label_kind for capacity in capacities] for capacities in seen] == [["circular_outer"]]
    assert CAPACITY_WARN_RATIO <= seen[0][0].record.ratio < CAPACITY_STOP_RATIO
    assert len([item for item in caught if issubclass(item.category, LabelsMayNotFitWarning)]) == 1
    assert [row.kind for row in stopped.value.report.rows if row.status == "found"] == ["smaller_font"]


# --- the CLI ---------------------------------------------------------------------------------------


@pytest.mark.circular
def test_cli_prints_the_stop_and_its_remedies(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    from gbdraw import cli

    monkeypatch.setattr(
        sys, "argv", ["gbdraw", "circular", *AP027078_ARGS, "-o", str(tmp_path / "ap"), "-f", "svg"]
    )
    with pytest.raises(SystemExit) as exited:
        cli.main()

    assert exited.value.code == 1
    assert capsys.readouterr().err.splitlines()[-6:] == [
        "ERROR: AP027078.1: 562 outer labels at 8 pt need about 2.5 times the room around the map, so they were not "
        "placed. Nothing was drawn.",
        "To make room:",
        "  - label fewer features with --label_blacklist or --label_whitelist;",
        "  - shorten the label text with --qualifier_priority (for example gene names only);",
        "  - use a smaller --label_font_size;",
        "  - or draw a Linear map (gbdraw linear).",
    ]


@pytest.mark.circular
def test_cli_prints_the_warning_before_placement(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    from gbdraw import cli

    _placement_raises(monkeypatch)
    monkeypatch.setattr(
        sys,
        "argv",
        ["gbdraw", "circular", *T5_ARGS, "--label_font_size", "14", "-o", str(tmp_path / "t5"), "-f", "svg"],
    )
    with pytest.raises(SystemExit):
        cli.main()

    err = capsys.readouterr().err.splitlines()
    assert err[-2:] == [
        "WARNING: NC_005859.1: 195 outer labels at 14 pt need about 1.4 times the room around the map, so they may "
        "not fit. Placing them can take minutes.",
        "ERROR: placement reached",
    ]

