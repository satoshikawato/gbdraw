"""The label fit check over every label kind (LP-D48 SCOPE-ALL-LABELS, LP-D49).

Inner labels, radial labels and labels over the centre definition text follow
the same rules as outer labels: overlapping labels stop the run, crossing
leaders or covered definition text alone draw and warn.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import pytest

import gbdraw.labels.circular_fit as circular_fit_module
from gbdraw.circular import circular_main
from gbdraw.exceptions import LabelLeadersCrossWarning, LabelsDoNotFitError

TEST_INPUTS = Path(__file__).parent / "test_inputs"
FONT_SIZE_PATHS = ("labels.font_size.short", "labels.font_size.long")


def _both_args(gbk: str, *extra: str) -> list[str]:
    return [
        "--gbk", str(TEST_INPUTS / gbk), "--separate_strands", "--track_type", "tuckin",
        "--labels", "both", "--legend", "right", *extra,
    ]


def _draw(tmp_path: Path, args: list[str]) -> Path:
    output = tmp_path / "out"
    circular_main([*args, "-o", str(output), "-f", "svg"])
    return output.with_suffix(".svg")


@pytest.mark.circular
def test_radial_labels_that_collide_after_the_growth_stops_are_not_drawn(tmp_path: Path) -> None:
    """OV-694 / RADIAL-JOIN: T5 radial labels drew 890 crossing leader pairs; now the run stops."""

    with pytest.raises(LabelsDoNotFitError) as raised:
        _draw(tmp_path, _both_args("NC_005859.1.gb", "--label_placement", "radial"))

    report = raised.value.report
    assert not report.drawn
    (record,) = report.records
    assert record.record_id == "NC_005859.1"
    assert record.refused
    assert [row.kind for row in report.rows] == ["smaller_font", "labels_out"]
    labels_out = report.rows[1]
    assert labels_out.status == "found"
    assert labels_out.values == {"labels.circular.scope": "outer"}
    assert "--labels out" in str(raised.value)
    assert not (tmp_path / "out.svg").exists()


@pytest.mark.circular
def test_inner_labels_over_the_definition_text_draw_and_warn(tmp_path: Path) -> None:
    """DEF-WARN (OV-696): lambda at 9 pt prints two inner labels over its definition text."""

    with pytest.warns(LabelLeadersCrossWarning) as caught:
        svg = _draw(tmp_path, _both_args("NC_001416.gb", "--label_font_size", "9"))

    assert svg.is_file()
    (warning,) = [item.message for item in caught if isinstance(item.message, LabelLeadersCrossWarning)]
    report = warning.report
    assert report.drawn
    (record,) = report.records
    assert record.label_kind == "circular_inner"
    assert (record.fit.overlap_pairs, record.fit.crossing_pairs) == (0, 0)
    assert record.fit.definition_overlaps == 2
    assert [row.kind for row in report.rows] == ["smaller_font", "labels_out", "radial"]
    assert str(warning).startswith(
        "24 inner labels do not fit at 9 pt (0 label pairs overlap, 0 leader pairs cross, "
        "2 labels cover the definition text). Drawn anyway."
    )


@pytest.mark.circular
@pytest.mark.slow
def test_overlapping_inner_labels_are_not_drawn(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """SAME-RULES: tobacco at 14 pt with inner labels has 47 overlapping inner label pairs."""

    # The font scan has its own test; here the switches are verified first (SWITCHES-FIRST).
    monkeypatch.setattr(circular_fit_module, "FONT_SCAN_SECONDS", 0.0, raising=False)
    with pytest.raises(LabelsDoNotFitError) as raised:
        _draw(tmp_path, _both_args("NC_001879.gbk", "--label_font_size", "14"))

    report = raised.value.report
    inner = [record for record in report.records if record.label_kind == "circular_inner"]
    assert inner and inner[0].fit.overlap_pairs > 0
    rows = {row.kind: row for row in report.rows}
    assert list(rows) == ["smaller_font", "labels_out", "radial"]
    assert rows["smaller_font"].status == "not_computed"
    assert rows["radial"].status == "found"
    assert rows["radial"].values == {"labels.circular.placement": "radial"}
    assert raised.value.row_search is not None


@pytest.mark.circular
@pytest.mark.slow
def test_inner_font_row_is_the_largest_font_that_fits(tmp_path: Path) -> None:
    """INNER-3: P1 at 14 pt fits at 12 pt, not at 11 pt, and at 10 pt again; the row offers 12 pt."""

    with pytest.raises(LabelsDoNotFitError) as raised:
        _draw(tmp_path, _both_args("NC_005856.1.gb", "--label_font_size", "14"))

    rows = {row.kind: row for row in raised.value.report.rows}
    font = rows["smaller_font"]
    if font.status == "not_computed":
        search = raised.value.row_search
        assert search is not None
        font = search.row("smaller_font")
    assert font.status == "found"
    assert font.values == {path: 12.0 for path in FONT_SIZE_PATHS}


@pytest.mark.circular
@pytest.mark.slow
def test_crossing_inner_leaders_alone_draw_and_warn(tmp_path: Path) -> None:
    """SAME-RULES with CROSS-WARN: lambda's inner leaders cross but no labels overlap."""

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", LabelLeadersCrossWarning)
        svg = _draw(tmp_path, _both_args("NC_001416.gb"))

    assert svg.is_file()
    (warning,) = [item.message for item in caught if isinstance(item.message, LabelLeadersCrossWarning)]
    inner = [record for record in warning.report.records if record.label_kind == "circular_inner"]
    assert inner and inner[0].fit.overlap_pairs == 0 and inner[0].fit.crossing_pairs > 0
