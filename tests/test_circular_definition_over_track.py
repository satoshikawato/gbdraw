"""D-12 (OV-1606): a Circular centre definition whose text reaches into a track draws as before and warns.

The warning is the label fit warning (``LabelLeadersCrossWarning``): the CLI prints ``WARNING:`` with the rows
that fit, Python warns with the same report, and the Web reads it as the notice with Choose a fix. The rows are a
smaller definition font (``--definition_font_size``), checked by drawing the whole figure again, and, when the
drawn text enters a track's band, narrower inside tracks (``--center_reserved_radius``, DF-D16 N-a). The warning
ends by suggesting a shorter Species or Strain.
"""

from __future__ import annotations

import dataclasses
import sys
import warnings
import xml.etree.ElementTree as ET
from collections.abc import Mapping
from pathlib import Path

import pytest
from Bio import SeqIO

import gbdraw.labels.circular_fit as circular_fit_module
from gbdraw.exceptions import LabelLeadersCrossWarning, LabelsDoNotFitError
from gbdraw.labels.circular_fit import CircularLabelFitSearch, CircularLabelLayout, check_circular_label_fit
from gbdraw.labels.fit import (
    DEFINITION_FONT_PATH,
    Y_OFFSET_PATH,
    LabelFit,
    LabelFitRecord,
    LabelFitReport,
    LabelFitRow,
    LabelFitValue,
    LabelFitVerification,
    format_label_fit_report,
    label_fit_notice,
)

TEST_INPUTS = Path(__file__).parent / "test_inputs"
MERS_COV = TEST_INPUTS / "MERS-CoV.gbk"
DEFINITION = LabelFit(0, 0, definition_track_overlaps=1)


def _definition_layout(
    record_id: str = "rec", *, record_index: int = 0, also: Mapping[str, LabelFit] | None = None
) -> CircularLabelLayout:
    kind_fits = {"circular_definition": DEFINITION, **(also or {})}
    return CircularLabelLayout(
        record_id=record_id,
        record_index=record_index,
        placement="horizontal",
        font_size=14.0,
        label_counts={"circular_definition": 1, **({"circular_outer": 40} if also else {})},
        kind_fits=kind_fits,
        x_radius_offset=1.0,
        y_radius_offset=1.0,
        label_spacing=3.0,
        probe=(lambda values: LabelFit(0, 0) if values.get(Y_OFFSET_PATH, 1.0) >= 1.1 else LabelFit(0, 2))
        if also
        else None,
        definition_font_size=18.0,
    )


def _verifier(definition_fits_below: float, *, crossings_until_y: float | None = None):
    calls: list[dict[str, LabelFitValue]] = []

    def verify(values: Mapping[str, LabelFitValue]) -> LabelFitVerification:
        calls.append(dict(values))
        font = float(values.get(DEFINITION_FONT_PATH, 18.0))
        y_offset = float(values.get(Y_OFFSET_PATH, 1.0))
        crossings = 2 if crossings_until_y is not None and y_offset < crossings_until_y else 0
        fit = LabelFit(0, crossings, definition_track_overlaps=0 if font <= definition_fits_below else 1)
        return LabelFitVerification("fits" if fit.fits else "warn", (1000.0, 1000.0), fit)

    return verify, calls


def test_definition_over_a_track_warns_with_the_mildest_definition_font(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(circular_fit_module.time, "monotonic", lambda: 0.0)
    verify, calls = _verifier(16.5)

    with pytest.warns(LabelLeadersCrossWarning) as caught:
        checked = check_circular_label_fit([_definition_layout(record_index=3)], verify)

    assert checked is not None
    report, _ = checked
    assert report.drawn and report.definition_track_overlaps == 1
    # The record carries the index its drawn definition group has (``data-gbdraw-record-index``).
    assert report.records == (LabelFitRecord("rec", "circular_definition", 1, 18.0, DEFINITION, record_index=3),)
    assert report.rows == (
        LabelFitRow("smaller_definition_font", "found", {DEFINITION_FONT_PATH: 16.5}, (1000.0, 1000.0)),
    )
    # The font is scanned down one step at a time with the full assembly.
    assert [call[DEFINITION_FONT_PATH] for call in calls] == [17.5, 17.0, 16.5]
    assert str(caught[0].message) == (
        "rec: 1 definition line does not fit at 18 pt (1 line overlaps a track). Drawn anyway.\n"
        "Run again with one of:\n"
        "  --definition_font_size 16.5   (smaller definition font)\n"
        "Or shorten the definition text with --species or --strain (for example, an abbreviation)."
    )


def test_each_row_answers_for_what_it_fixes(monkeypatch: pytest.MonkeyPatch) -> None:
    """Crossing leaders and the definition over a track on one record: the Y offset clears the crossings and the
    font clears the definition, so both rows are found although neither clears both."""

    monkeypatch.setattr(circular_fit_module.time, "monotonic", lambda: 0.0)
    layout = _definition_layout(also={"circular_outer": LabelFit(0, 2)})
    verify, _ = _verifier(17.0, crossings_until_y=1.1)

    with pytest.warns(LabelLeadersCrossWarning):
        checked = check_circular_label_fit([layout], verify)

    assert checked is not None
    rows = {row.kind: row for row in checked[0].rows}
    assert rows["outward_vertical"].status == "found"
    assert rows["outward_vertical"].values == {Y_OFFSET_PATH: 1.1}
    assert rows["smaller_definition_font"].values == {DEFINITION_FONT_PATH: 17.0}
    assert {record.label_kind for record in checked[0].records} == {"circular_outer", "circular_definition"}


def test_one_definition_font_fits_every_record(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(circular_fit_module.time, "monotonic", lambda: 0.0)
    verify, _ = _verifier(15.0)

    with pytest.warns(LabelLeadersCrossWarning) as caught:
        checked = check_circular_label_fit([_definition_layout("A"), _definition_layout("B")], verify)

    assert checked is not None
    assert checked[0].rows[0].values == {DEFINITION_FONT_PATH: 15.0}
    assert str(caught[0].message).startswith(
        "Definition lines do not fit in 2 records: A (1 line at 18 pt: 1 line overlaps a track); "
        "B (1 line at 18 pt: 1 line overlaps a track). Drawn anyway.\n"
    )


def test_a_figure_that_was_not_drawn_offers_no_definition_row(monkeypatch: pytest.MonkeyPatch) -> None:
    """Overlapping labels stop the figure; their rows come first, the definition row after the next Generate."""

    monkeypatch.setattr(circular_fit_module.time, "monotonic", lambda: 0.0)
    layout = _definition_layout(also={"circular_outer": LabelFit(3, 0)})
    verify, _ = _verifier(17.0)

    with pytest.raises(LabelsDoNotFitError) as caught:
        check_circular_label_fit([layout], verify)

    assert "smaller_definition_font" not in {row.kind for row in caught.value.report.rows}


def test_the_notice_names_the_definition_text_over_a_track() -> None:
    record = LabelFitRecord("rec", "circular_definition", 2, 18.0, LabelFit(0, 0, definition_track_overlaps=2))
    report = LabelFitReport(records=(record,), rows=(), drawn=True)

    assert label_fit_notice(report) == ("The definition text overlaps a track (2 lines).",)
    assert format_label_fit_report(report).startswith(
        "rec: 2 definition lines do not fit at 18 pt (2 lines overlap a track). Drawn anyway."
    )


@pytest.mark.circular
def test_cli_draws_mers_cov_with_the_warning_and_the_suggested_font_clears_it(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """OV-1606, DF-D16, DF-D25, OV-1836, DF-D36: MERS-CoV's taxon line is drawn on one line (never wrapped) and
    warns. Its band (the farthest corner of the bold line's ink, OV-1854) is not cleared by narrowing the inside
    tracks to the floor, so they keep their widths, and no center_reserved_radius fits either: the warning offers
    the font at which the line clears the tracks (17 pt, where the floor clears the smaller band) and a shorter
    Species or Strain. With ``--no_track_narrowing`` the font must clear the tracks at their widths (13 pt)."""

    from gbdraw import cli

    monkeypatch.setattr(circular_fit_module.time, "monotonic", lambda: 0.0)

    def draw(name: str, *args: str) -> str:
        output = tmp_path / name
        monkeypatch.setattr(
            sys, "argv", ["gbdraw", "circular", "--gbk", str(MERS_COV), *args, "-o", str(output), "-f", "svg"]
        )
        cli.main()
        assert output.with_suffix(".svg").is_file()
        return capsys.readouterr().err

    assert draw("mers") == (
        "WARNING: NC_019843.3: 1 definition line does not fit at 18 pt (1 line overlaps a track). Drawn anyway.\n"
        "Run again with one of:\n"
        "  --definition_font_size 17   (smaller definition font)\n"
        "Narrowing the inside tracks does not fit this map.\n"
        "Or shorten the definition text with --species or --strain (for example, an abbreviation).\n"
    )
    assert draw("mers_font", "--definition_font_size", "17") == ""
    # The size above it does not clear the one line: the row offers the mildest size.
    assert draw("mers_font_17_5", "--definition_font_size", "17.5").startswith(
        "WARNING: NC_019843.3: 1 definition line does not fit at 17.5 pt (1 line overlaps a track). Drawn anyway.\n"
    )

    assert draw("mers_off", "--no_track_narrowing") == (
        "WARNING: NC_019843.3: 1 definition line does not fit at 18 pt (1 line overlaps a track). Drawn anyway.\n"
        "Run again with one of:\n"
        "  --definition_font_size 13   (smaller definition font)\n"
        "Narrowing the inside tracks does not fit this map.\n"
        "Or shorten the definition text with --species or --strain (for example, an abbreviation).\n"
    )
    assert draw("mers_off_font", "--no_track_narrowing", "--definition_font_size", "13") == ""
    assert draw("mers_off_font_13_5", "--no_track_narrowing", "--definition_font_size", "13.5").startswith(
        "WARNING: NC_019843.3: 1 definition line does not fit at 13.5 pt (1 line overlaps a track). Drawn anyway.\n"
    )


@pytest.mark.circular
def test_cli_names_the_record_of_each_file_it_warns_for(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    """OV-1836: records drawn to separate files warn one by one; each warning names its record."""

    import copy

    from gbdraw import cli

    monkeypatch.setattr(circular_fit_module.time, "monotonic", lambda: 0.0)
    first = SeqIO.read(str(MERS_COV), "genbank")
    second = copy.deepcopy(first)
    second.id = second.name = "MERS_COPY"
    gbk = tmp_path / "two.gb"
    SeqIO.write([first, second], str(gbk), "genbank")
    monkeypatch.setattr(sys, "argv", ["gbdraw", "circular", "--gbk", str(gbk), "-o", str(tmp_path / "out"), "-f", "svg"])
    cli.main()

    warnings_printed = [line for line in capsys.readouterr().err.splitlines() if line.startswith("WARNING: ")]
    assert [line.split(": ", 2)[1] for line in warnings_printed] == ["NC_019843.3", "MERS_COPY"]


@pytest.mark.circular
def test_python_warns_with_the_report_and_a_figure_without_the_overlap_does_not() -> None:
    from gbdraw.api.diagram import build_circular_diagram, label_fit_report, label_fit_row
    from gbdraw.api.options import CircularDiagramOptions

    record = SeqIO.read(str(MERS_COV), "genbank")
    with pytest.warns(LabelLeadersCrossWarning):
        drawing = build_circular_diagram(record, options=CircularDiagramOptions())
    report = label_fit_report(drawing)
    assert report is not None and report.drawn
    assert [(item.label_kind, item.fit) for item in report.records] == [("circular_definition", DEFINITION)]
    row = label_fit_row(drawing, "smaller_definition_font")
    assert row is not None and row.values == {DEFINITION_FONT_PATH: 17.0}

    from gbdraw.config.models import GbdrawConfig
    from gbdraw.config.modify import with_config_overrides
    from gbdraw.config.toml import load_config_toml

    cfg = with_config_overrides(
        GbdrawConfig.from_dict(load_config_toml("gbdraw.data", "config.toml")), {DEFINITION_FONT_PATH: 17.0}
    )
    with warnings.catch_warnings():
        warnings.simplefilter("error", LabelLeadersCrossWarning)
        drawing = build_circular_diagram(record, options=CircularDiagramOptions(config=cfg))
    assert label_fit_report(drawing) is None


def test_the_probe_picks_the_first_step_to_verify_and_a_milder_verified_step_wins(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """The probe (tracks before the labels) says 15.5 pt; the full figure fits from 16.5 pt, which is offered."""

    monkeypatch.setattr(circular_fit_module.time, "monotonic", lambda: 0.0)
    probed: list[float] = []

    def probe(value: float) -> bool:
        probed.append(value)
        return value <= 15.5

    layout = dataclasses.replace(_definition_layout(), definition_probe=probe)
    verify, calls = _verifier(16.5)

    with pytest.warns(LabelLeadersCrossWarning):
        checked = check_circular_label_fit([layout], verify)

    assert checked is not None
    assert checked[0].rows[0].values == {DEFINITION_FONT_PATH: 16.5}
    assert probed == [17.5, 17.0, 16.5, 16.0, 15.5]
    assert [call[DEFINITION_FONT_PATH] for call in calls] == [15.5, 16.0, 16.5, 17.0]


def test_a_verified_definition_row_survives_the_deadline_and_a_later_ask_walks_on(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Review F2: the walk to a milder step stops at the deadline; the step already verified is the row.

    The probe says 16.5 pt and each full redraw takes 2 s of a 5 s budget: 16.5 and 17 are verified, 17.5 would
    overrun. The row is 17 pt, not ``not_computed``, and is not kept: asked again without a deadline, it walks on to
    the mildest step (17.5 pt).
    """

    clock = [0.0]
    monkeypatch.setattr(circular_fit_module.time, "monotonic", lambda: clock[0])
    layout = dataclasses.replace(_definition_layout(), definition_probe=lambda value: value <= 16.5)
    fits_below, calls = _verifier(17.5)

    def verify(values: Mapping[str, LabelFitValue]) -> LabelFitVerification:
        clock[0] += 2.0
        return fits_below(values)

    search = CircularLabelFitSearch([layout], verify, offer_notice_rows=False, assembly_seconds=2.0)

    row = search.row("smaller_definition_font", deadline=5.0)
    assert (row.status, row.values) == ("found", {DEFINITION_FONT_PATH: 17.0})
    assert [call[DEFINITION_FONT_PATH] for call in calls] == [16.5, 17.0]

    row = search.row("smaller_definition_font")
    assert (row.status, row.values) == ("found", {DEFINITION_FONT_PATH: 17.5})
    assert [call[DEFINITION_FONT_PATH] for call in calls] == [16.5, 17.0, 17.5]


def test_a_record_failing_only_by_its_definition_leaves_the_label_rows_alone() -> None:
    """Review F5: the label font scan starts from the records whose labels fail, not from a definition-only one."""

    definition_only = dataclasses.replace(_definition_layout("long"), font_size=10.0, label_spacing=1.0)
    labels = CircularLabelLayout(
        record_id="short",
        record_index=1,
        placement="horizontal",
        font_size=14.0,
        label_counts={"circular_outer": 40},
        kind_fits={"circular_outer": LabelFit(0, 2)},
        x_radius_offset=1.0,
        y_radius_offset=1.0,
        label_spacing=3.0,
    )
    calls: list[dict[str, LabelFitValue]] = []

    def verify(values: Mapping[str, LabelFitValue]) -> LabelFitVerification:
        calls.append(dict(values))
        return LabelFitVerification("fits", (1000.0, 1000.0), LabelFit(0, 0))

    search = CircularLabelFitSearch([definition_only, labels], verify, offer_notice_rows=False)
    row = search.row("smaller_font")

    assert row.values == {"labels.font_size.short": 13.5, "labels.font_size.long": 13.5}


@pytest.mark.circular
def test_the_row_keeps_a_set_definition_interval_as_a_run_again_does() -> None:
    """Review F4: with an interval the user set, the row is checked with that interval, so drawing again with the
    row's font and the same interval does not warn."""

    from gbdraw.api.diagram import build_circular_diagram, label_fit_row
    from gbdraw.api.options import CircularDiagramOptions
    from gbdraw.config.models import GbdrawConfig
    from gbdraw.config.modify import with_config_overrides
    from gbdraw.config.toml import load_config_toml

    record = SeqIO.read(str(MERS_COV), "genbank")
    base = GbdrawConfig.from_dict(load_config_toml("gbdraw.data", "config.toml"))
    interval = {"objects.definition.circular.interval": 25}
    with pytest.warns(LabelLeadersCrossWarning):
        drawing = build_circular_diagram(
            record, options=CircularDiagramOptions(config=with_config_overrides(base, interval))
        )
    row = label_fit_row(drawing, "smaller_definition_font")
    assert row is not None and row.status == "found"

    with warnings.catch_warnings():
        warnings.simplefilter("error", LabelLeadersCrossWarning)
        build_circular_diagram(
            record, options=CircularDiagramOptions(config=with_config_overrides(base, {**interval, **row.values}))
        )


@pytest.mark.circular
@pytest.mark.parametrize(
    ("name", "args"),
    [("MjeNMV.gbk", ["--separate_strands", "--legend", "none", "--species", "MjeNMV"]), ("EDL933.gbk", [])],
    ids=["mjenmv-abbreviated", "edl933-default"],
)
def test_a_definition_clear_of_the_drawn_tracks_does_not_warn(
    name: str, args: list[str], tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    """D-18: the text stays clear of the tracks' paint (EDL933 enters the GC skew ring's band): no warning. MjeNMV's
    full name overlaps the GC rings on one line (next test); the abbreviation the docs advise clears them."""

    from gbdraw import cli

    output = tmp_path / "out"
    monkeypatch.setattr(
        sys, "argv", ["gbdraw", "circular", "--gbk", str(TEST_INPUTS / name), *args, "-o", str(output), "-f", "svg"]
    )
    cli.main()

    assert output.with_suffix(".svg").is_file()
    assert "definition" not in capsys.readouterr().err


@pytest.mark.circular
def test_mjenmv_full_name_is_drawn_on_one_line_and_warns_only_without_narrowing(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    """DF-D25: "Marsupenaeus japonicus endogenous nimavirus" stays one line. Rule L narrows the GC rows (to about
    35.5 px, above the 30 px floor) so it clears them; with narrowing off it is drawn over them and the warning offers
    the one-line font."""

    from gbdraw import cli

    monkeypatch.setattr(circular_fit_module.time, "monotonic", lambda: 0.0)
    output = tmp_path / "out"
    argv = ["gbdraw", "circular", "--gbk", str(TEST_INPUTS / "MjeNMV.gbk"), "--separate_strands", "--legend", "none",
            "-o", str(output), "-f", "svg"]
    monkeypatch.setattr(sys, "argv", argv)
    cli.main()
    assert "definition line" not in capsys.readouterr().err

    monkeypatch.setattr(sys, "argv", [*argv, "--no_track_narrowing", "--overwrite"])
    cli.main()
    assert capsys.readouterr().err.endswith(
        "WARNING: LC738868.1: 1 definition line does not fit at 18 pt (1 line overlaps a track). Drawn anyway.\n"
        "Run again with one of:\n"
        "  --definition_font_size 14.5       (smaller definition font)\n"
        "  --center_reserved_radius 233.39   (narrow the inside tracks)\n"
        "Or shorten the definition text with --species or --strain (for example, an abbreviation).\n"
    )
    root = ET.fromstring(output.with_suffix(".svg").read_text())
    group = root.find(".//{http://www.w3.org/2000/svg}g[@data-gbdraw-role='record-definition']")
    assert group is not None
    first = group.find("{http://www.w3.org/2000/svg}text")
    assert first is not None and "".join(first.itertext()) == "Marsupenaeus japonicus endogenous nimavirus"


def test_a_box_in_a_track_band_meets_only_its_paint() -> None:
    """D-18: the band a track takes is not its paint. A GC area drawn outward from the circle at 100 px leaves the
    band's inner part bare: a box reaching to 98 px misses it, one reaching to 101 px meets it, and 0.5 px is the
    threshold."""

    from gbdraw.layout.circular import RadialBand
    from gbdraw.layout.circular_ink import CircularInkProfile, CircularInkRing, box_meets_ink
    from gbdraw.layout.spatial import Aabb

    outward = CircularInkProfile.between([0, 250, 500, 750], [110.0, 104.0, 110.0, 104.0], record_length=1000,
                                         base_px=100.0)
    band = CircularInkRing(RadialBand(90.0, 120.0))

    def box(reach: float) -> Aabb:
        """A narrow box from the centre up to ``reach`` (its corners reach no farther than 0.1 px past it)."""

        return Aabb(-2.0, -reach, 2.0, 0.0)

    assert box_meets_ink(box(95.0), [band], threshold_px=0.5)
    assert not box_meets_ink(box(98.0), [outward], threshold_px=0.5)
    assert not box_meets_ink(box(100.4), [outward], threshold_px=0.5)
    assert box_meets_ink(box(101.0), [outward], threshold_px=0.5)
    # A ring over part of the circle: a box under its gap misses it.
    arc = CircularInkRing(RadialBand(100.0, 110.0), start_turn=0.25, end_turn=0.75)
    assert not box_meets_ink(box(105.0), [arc], threshold_px=0.5)
    assert box_meets_ink(Aabb(-2.0, 0.0, 2.0, 105.0), [arc], threshold_px=0.5)


# --- Every track kind is judged by its paint (D-18; PD-OI-097 scenario revision 17, Must preserve) ---------------

HMMTDNA = TEST_INPUTS / "HmmtDNA.gbk"


def _drawn_track_ink(monkeypatch: pytest.MonkeyPatch, draw) -> tuple[object, tuple[object, ...]]:
    """Draw a figure and return the radial layout and the paint the definition check read for its record."""

    import gbdraw.diagrams.circular.assemble as assemble_module

    captured: list[tuple[object, tuple[object, ...]]] = []
    original = assemble_module._record_track_ink

    def capture(layout, **kwargs):
        inks = original(layout, **kwargs)
        captured.append((layout, inks))
        return inks

    monkeypatch.setattr(assemble_module, "_record_track_ink", capture)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", LabelLeadersCrossWarning)
        draw()
    assert captured, "the definition check did not run"
    return captured[-1]


def _slot(layout, renderer: str):
    return next(slot for slot in layout.slots if str(slot.renderer) == renderer)


def _lines_at(inks, radius: float, turn: float) -> int:
    """How many lines a small definition box at ``radius`` and display ``turn`` counts against the paint."""

    import math

    from gbdraw.diagrams.circular.assemble import _definition_lines_on_tracks
    from gbdraw.layout.spatial import Aabb

    angle = 2.0 * math.pi * turn - 0.5 * math.pi
    x, y = radius * math.cos(angle), radius * math.sin(angle)
    return _definition_lines_on_tracks(inks, [Aabb(x - 0.75, y - 0.75, x + 0.75, y + 0.75)])


def _assemble(**kwargs):
    from gbdraw.api.diagram import assemble_circular_diagram_from_record
    from gbdraw.config.models import GbdrawConfig
    from gbdraw.config.modify import modify_config_dict
    from gbdraw.config.toml import load_config_toml

    overrides = kwargs.pop("overrides", {})
    record = SeqIO.read(str(HMMTDNA), "genbank")
    cfg = GbdrawConfig.from_dict(modify_config_dict(load_config_toml("gbdraw.data", "config.toml"), overrides))
    return lambda: assemble_circular_diagram_from_record(record, cfg=cfg, legend="none", **kwargs)


@pytest.mark.circular
def test_a_depth_track_counts_its_area_and_axis_not_its_band(monkeypatch: pytest.MonkeyPatch) -> None:
    import pandas as pd

    record = SeqIO.read(str(HMMTDNA), "genbank")
    length = len(record.seq)
    positions = list(range(1, length + 1, 10))
    depth = pd.DataFrame(
        {"reference_name": record.id, "position": positions, "depth": [80 if p <= length // 4 else 0 for p in positions]}
    )
    layout, inks = _drawn_track_ink(
        monkeypatch,
        _assemble(
            depth_table=depth, depth_window=100, depth_step=100,
            circular_track_slots=["features:features", "depth:depth"], circular_track_axis_index=1,
        ),
    )
    slot = _slot(layout, "depth")
    band, drawn = slot.own_reserved_band_px, slot.draw_band_px
    middle = 0.5 * (drawn.inner_px + drawn.outer_px)
    # Where there is no coverage the area is bare: inside the band, off the paint.
    assert _lines_at(inks, middle, 0.6) == 0
    # The room the band keeps inward for the axis tick labels holds no paint away from the axis.
    assert band.inner_px < drawn.inner_px - 1.5
    assert _lines_at(inks, 0.5 * (band.inner_px + drawn.inner_px) - 0.2, 0.6) == 0
    # The covered quarter and the baseline circle are paint.
    assert _lines_at(inks, middle, 0.1) == 1
    assert _lines_at(inks, drawn.inner_px, 0.6) == 1


@pytest.mark.circular
def test_a_conservation_ring_counts_its_hits_and_its_background(monkeypatch: pytest.MonkeyPatch) -> None:
    import pandas as pd

    from gbdraw.io.comparisons import COMPARISON_COLUMNS

    record = SeqIO.read(str(HMMTDNA), "genbank")
    hit = ("query1", record.id, 90.0, 4000, 0, 0, 1, 4000, 1, 4000, 1e-20, 100.0)
    slots = ["features:features", "conservation_1:sequence_conservation@source_index=0"]

    def ring(show_background: bool):
        return _drawn_track_ink(
            monkeypatch,
            _assemble(
                overrides={"objects.conservation.show_background": show_background},
                conservation_dataframes=[pd.DataFrame([hit], columns=COMPARISON_COLUMNS)],
                conservation_reference="subject",
                circular_track_slots=slots, circular_track_axis_index=1,
            ),
        )

    layout, inks = ring(False)
    drawn = _slot(layout, "sequence_conservation").draw_band_px
    middle = 0.5 * (drawn.inner_px + drawn.outer_px)
    assert _lines_at(inks, middle, 0.6) == 0
    assert _lines_at(inks, middle, 0.1) == 1
    layout, inks = ring(True)
    assert _lines_at(inks, middle, 0.6) == 1


@pytest.mark.circular
def test_an_annotation_track_counts_its_marks_not_its_band(monkeypatch: pytest.MonkeyPatch) -> None:
    from gbdraw.api import AnnotationOptions, AnnotationSet, CoordinateSpan, RegionAnnotation, RegionAnnotationStyle

    record = SeqIO.read(str(HMMTDNA), "genbank")
    length = len(record.seq)
    region = RegionAnnotation(
        "block", CoordinateSpan(None, 1, length // 4), mark="band", style=RegionAnnotationStyle(fill="#cccccc")
    )
    layout, inks = _drawn_track_ink(
        monkeypatch,
        _assemble(
            annotation_options=AnnotationOptions(sets=(AnnotationSet("regions", (region,)),)),
            circular_track_slots=["features:features", "regions:annotations@side=inside,set_id=regions,w=20px"],
            circular_track_axis_index=1,
        ),
    )
    drawn = _slot(layout, "annotations").draw_band_px
    middle = 0.5 * (drawn.inner_px + drawn.outer_px)
    assert _lines_at(inks, middle, 0.6) == 0
    assert _lines_at(inks, middle, 0.1) == 1


@pytest.mark.circular
def test_a_tick_track_counts_its_marks_and_labels_not_its_band(monkeypatch: pytest.MonkeyPatch) -> None:
    from gbdraw.svg.circular_ticks import get_circular_tick_intervals

    record = SeqIO.read(str(HMMTDNA), "genbank")
    length = len(record.seq)
    layout, inks = _drawn_track_ink(
        monkeypatch,
        _assemble(circular_track_slots=["features:features", "ticks:ticks"], circular_track_axis_index=1),
    )
    tick_layout = _slot(layout, "ticks").payload
    large, _ = get_circular_tick_intervals(length)
    marks = 0.5 * (tick_layout.tick_band_px.inner_px + tick_layout.tick_band_px.outer_px)
    # Between two ticks there is neither a mark nor a label.
    assert _lines_at(inks, marks, 1.5 * large / length) == 0
    assert tick_layout.label_band_px is not None
    labels = 0.5 * (tick_layout.label_band_px.inner_px + tick_layout.label_band_px.outer_px)
    assert _lines_at(inks, labels, 1.5 * large / length) == 0
    # A mark and a label are paint.
    assert _lines_at(inks, marks, large / length) == 1
    assert _lines_at(inks, labels, large / length) == 1


@pytest.mark.circular
def test_the_row_keeps_a_set_interval_equal_to_the_one_the_font_derives() -> None:
    """Review of 1f079028: an interval set in ``config_overrides`` (as the Web field sends it) to the value the
    font derives (20 at 18 pt) is the user's and is kept while the row is checked, so the row fits when drawn
    again with that interval."""

    from gbdraw.api.diagram import build_circular_diagram, label_fit_row
    from gbdraw.api.options import CircularDiagramOptions

    record = SeqIO.read(str(MERS_COV), "genbank")
    interval = {"objects.definition.circular.interval": 20}
    with pytest.warns(LabelLeadersCrossWarning):
        drawing = build_circular_diagram(record, options=CircularDiagramOptions(config_overrides=interval))
    row = label_fit_row(drawing, "smaller_definition_font")
    assert row is not None and row.status == "found"

    with warnings.catch_warnings():
        warnings.simplefilter("error", LabelLeadersCrossWarning)
        build_circular_diagram(record, options=CircularDiagramOptions(config_overrides={**interval, **row.values}))


def _base_config(overrides: Mapping[str, object] | None = None):
    from gbdraw.config.models import GbdrawConfig
    from gbdraw.config.modify import modify_config_dict
    from gbdraw.config.toml import load_config_toml

    return GbdrawConfig.from_dict(modify_config_dict(load_config_toml("gbdraw.data", "config.toml"), overrides or {}))


@pytest.mark.parametrize(
    ("background", "lines"),
    [({}, 1), ({"background_opacity": 0}, 0), ({"background_color": "none"}, 0)],
    ids=["shown", "transparent", "no-color"],
)
def test_a_conservation_background_counts_only_when_it_shows(background: Mapping[str, object], lines: int) -> None:
    """Review of 5d358bd2: a background with no color or no opacity is not paint; the ring then counts its hits."""

    import types

    import pandas as pd

    from gbdraw.render.groups.circular.conservation import conservation_ring_ink

    cfg = _base_config({f"objects.conservation.{key}": value for key, value in background.items()})
    track = types.SimpleNamespace(track_color=None, hits=pd.DataFrame())
    inks = conservation_ring_ink(
        track, total_length=1000, inner_radius_px=90.0, outer_radius_px=110.0, min_identity=70.0, cfg=cfg
    )
    assert _lines_at(inks, 100.0, 0.5) == lines


@pytest.mark.parametrize(
    ("style", "radius", "lines"),
    [
        ({"percent_background_color": "#737373"}, 100.0, 1),
        ({"percent_background_opacity": 0.0}, 100.0, 0),
        ({"percent_border_color": "#4b5563"}, 110.0, 1),
        ({"percent_border_color": "none"}, 110.0, 0),
    ],
    ids=["background-shown", "background-transparent", "border-shown", "border-no-color"],
)
def test_a_percent_background_and_border_count_only_when_they_show(
    style: Mapping[str, object], radius: float, lines: int
) -> None:
    """Review of 5d358bd2: the percent background and border circles are paint only with a color and opacity."""

    import types

    import pandas as pd

    from gbdraw.render.drawers.circular.gc_content import GcContentDrawer

    settings = {
        "high_fill_color": "#000", "low_fill_color": "#000", "fill_color": "#000", "stroke_color": "none",
        "stroke_width": 0.0, "fill_opacity": 1.0, "mode": "percent", "show_axis": False,
        "percent_background_color": "none", "percent_background_opacity": 1.0,
        "percent_border_color": "#4b5563", "percent_border_width": 2.0, **style,
    }
    if "percent_border_color" in style:
        settings["percent_background_color"] = "none"
    drawer = GcContentDrawer(types.SimpleNamespace(**settings))
    gc_df = pd.DataFrame({"GC content": [0.0] * 4}, index=[0, 250, 500, 750])
    inks = drawer.ink(100.0, gc_df, 1000, 20.0, 1.0, "GC", font_family="Liberation Sans")
    assert _lines_at(inks, radius, 0.5) == lines


WEB_DEFAULT_ARGS = ("--multi_record_canvas", "--track_type", "tuckin", "--gc", "--skew", "--separate_strands")


@pytest.mark.circular
@pytest.mark.parametrize(
    ("gbk", "args"),
    [
        (TEST_INPUTS / "BGC0000712.gbk", WEB_DEFAULT_ARGS),
        (TEST_INPUTS / "NC_002745.gb", WEB_DEFAULT_ARGS),
    ],
    ids=["bgc0000712-web", "nc002745-web"],
)
def test_tracks_narrowed_around_the_definition_leave_no_notice(
    gbk: Path,
    args: tuple[str, ...],
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """OV-1854: the band the inside tracks narrow around and the notice read the same painted letters, so a figure
    whose tracks were narrowed (rule L) has no definition line on a track. Both figures were narrowed with the
    notice while the band read the line boxes."""

    import gbdraw.diagrams.circular.assemble as circular_assemble_module
    from gbdraw import cli

    fits: list[circular_assemble_module._CentreDefinitionFit] = []
    fit = circular_assemble_module._resolve_radial_layout_fitting_definition

    def capture(*fit_args, **fit_kwargs):
        result = fit(*fit_args, **fit_kwargs)
        fits.append(result)
        return result

    monkeypatch.setattr(circular_assemble_module, "_resolve_radial_layout_fitting_definition", capture)
    monkeypatch.setattr(circular_fit_module.time, "monotonic", lambda: 0.0)
    output = tmp_path / "out"
    monkeypatch.setattr(sys, "argv", ["gbdraw", "circular", "--gbk", str(gbk), *args, "-o", str(output), "-f", "svg"])
    cli.main()

    narrowed = any(track.compressed for track in fits[0].layout.tracks)
    notice = "definition line does not fit" in capsys.readouterr().err
    assert narrowed and not notice, (narrowed, notice)
