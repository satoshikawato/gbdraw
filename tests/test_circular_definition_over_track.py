"""D-12 (OV-1606): a Circular centre definition whose text reaches into a track draws as before and warns.

The warning is the label fit warning (``LabelLeadersCrossWarning``): the CLI prints ``WARNING:`` with the rows
that fit, Python warns with the same report, and the Web reads it as the notice with Choose a fix. The one row is
a smaller definition font (``--definition_font_size``), checked by drawing the whole figure again.
"""

from __future__ import annotations

import dataclasses
import sys
import warnings
from collections.abc import Mapping
from pathlib import Path

import pytest
from Bio import SeqIO

import gbdraw.labels.circular_fit as circular_fit_module
from gbdraw.exceptions import LabelLeadersCrossWarning, LabelsDoNotFitError
from gbdraw.labels.circular_fit import CircularLabelLayout, check_circular_label_fit
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
AP027131 = TEST_INPUTS / "AP027131.gb"
DEFINITION = LabelFit(0, 0, definition_track_overlaps=1)


def _definition_layout(record_id: str = "rec", *, also: Mapping[str, LabelFit] | None = None) -> CircularLabelLayout:
    kind_fits = {"circular_definition": DEFINITION, **(also or {})}
    return CircularLabelLayout(
        record_id=record_id,
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
        checked = check_circular_label_fit([_definition_layout()], verify)

    assert checked is not None
    report, _ = checked
    assert report.drawn and report.definition_track_overlaps == 1
    assert report.records == (LabelFitRecord("rec", "circular_definition", 1, 18.0, DEFINITION),)
    assert report.rows == (
        LabelFitRow("smaller_definition_font", "found", {DEFINITION_FONT_PATH: 16.5}, (1000.0, 1000.0)),
    )
    # The font is scanned down one step at a time with the full assembly.
    assert [call[DEFINITION_FONT_PATH] for call in calls] == [17.5, 17.0, 16.5]
    assert str(caught[0].message) == (
        "1 definition line does not fit at 18 pt (1 line overlaps a track). Drawn anyway.\n"
        "Run again with one of:\n"
        "  --definition_font_size 16.5   (smaller definition font)"
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
        "2 definition lines do not fit at 18 pt (2 lines overlap a track). Drawn anyway."
    )


@pytest.mark.circular
def test_cli_draws_ap027131_with_the_warning_and_the_suggested_font_clears_it(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    """OV-1606: the bold first line's corners reach into the GC skew ring of the default AP027131 figure."""

    from gbdraw import cli

    monkeypatch.setattr(circular_fit_module.time, "monotonic", lambda: 0.0)
    output = tmp_path / "ap"
    monkeypatch.setattr(sys, "argv", ["gbdraw", "circular", "--gbk", str(AP027131), "-o", str(output), "-f", "svg"])
    cli.main()

    assert output.with_suffix(".svg").is_file()
    assert capsys.readouterr().err == (
        "WARNING: 1 definition line does not fit at 18 pt (1 line overlaps a track). Drawn anyway.\n"
        "Run again with one of:\n"
        "  --definition_font_size 16.5   (smaller definition font)\n"
    )

    fixed = tmp_path / "ap_fixed"
    monkeypatch.setattr(
        sys,
        "argv",
        ["gbdraw", "circular", "--gbk", str(AP027131), "--definition_font_size", "16.5", "-o", str(fixed), "-f", "svg"],
    )
    cli.main()
    assert fixed.with_suffix(".svg").is_file()
    assert capsys.readouterr().err == ""


@pytest.mark.circular
def test_python_warns_with_the_report_and_a_figure_without_the_overlap_does_not() -> None:
    from gbdraw.api.diagram import build_circular_diagram, label_fit_report, label_fit_row
    from gbdraw.api.options import CircularDiagramOptions

    record = SeqIO.read(str(AP027131), "genbank")
    with pytest.warns(LabelLeadersCrossWarning):
        drawing = build_circular_diagram(record, options=CircularDiagramOptions())
    report = label_fit_report(drawing)
    assert report is not None and report.drawn
    assert [(item.label_kind, item.fit) for item in report.records] == [("circular_definition", DEFINITION)]
    row = label_fit_row(drawing, "smaller_definition_font")
    assert row is not None and row.values == {DEFINITION_FONT_PATH: 16.5}

    from gbdraw.config.models import GbdrawConfig
    from gbdraw.config.modify import with_config_overrides
    from gbdraw.config.toml import load_config_toml

    cfg = with_config_overrides(
        GbdrawConfig.from_dict(load_config_toml("gbdraw.data", "config.toml")), {DEFINITION_FONT_PATH: 16.5}
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
