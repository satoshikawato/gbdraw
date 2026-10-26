"""The label fit check over every label kind (LP-D48 SCOPE-ALL-LABELS, LP-D49).

Inner labels, radial labels and labels over the centre definition text follow
the same rules as outer labels: overlapping labels stop the run, crossing
leaders or covered definition text alone draw and warn.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

import gbdraw.diagrams.circular.assemble as circular_assemble_module
import gbdraw.diagrams.circular.builders as circular_builders_module
import gbdraw.labels.circular_fit as circular_fit_module
import gbdraw.labels.circular_repair as circular_repair_module
from gbdraw.api.config import apply_config_overrides
from gbdraw.api.options import CircularDiagramOptions
from gbdraw.api.request_render import plan_request
from gbdraw.api.requests import CircularDiagramRequest, InMemoryRecordSource, RecordInput, RecordPresentation
from gbdraw.canvas import CircularCanvasConfigurator
from gbdraw.circular import circular_main
from gbdraw.exceptions import LabelLeadersCrossWarning, LabelsDoNotFitError
from gbdraw.layout.spatial import Aabb
from gbdraw.render.groups.circular import DefinitionGroup
from gbdraw.svg.elements import Group

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
def test_inner_labels_over_the_definition_text_draw_and_warn(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """DEF-WARN (OV-696): lambda at 10 pt prints two inner labels over its definition text (with the repair off).

    9 pt did while labels were laid out in point-mode boxes; in their drawn boxes it covers it with one (OV-165)."""

    monkeypatch.setattr(circular_repair_module, "MAX_DEFECTS", 0)
    with pytest.warns(LabelLeadersCrossWarning) as caught:
        svg = _draw(tmp_path, _both_args("NC_001416.gb", "--label_font_size", "10"))

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
        "NC_001416.1: 25 inner labels do not fit at 10 pt (0 label pairs overlap, 0 leader pairs cross, "
        "2 labels cover the definition text). Drawn anyway."
    )


@pytest.mark.circular
@pytest.mark.slow
def test_overlapping_inner_labels_are_not_drawn(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """SAME-RULES: tobacco at 14 pt with inner labels has 47 overlapping inner label pairs."""

    # The search budget has its own test; here every row is left for later and one is
    # computed on demand from the error.
    monkeypatch.setattr(circular_fit_module, "SWITCH_SEARCH_SECONDS", 0.0)
    with pytest.raises(LabelsDoNotFitError) as raised:
        _draw(tmp_path, _both_args("NC_001879.gbk", "--label_font_size", "14"))

    report = raised.value.report
    inner = [record for record in report.records if record.label_kind == "circular_inner"]
    assert inner and inner[0].fit.overlap_pairs > 0
    rows = {row.kind: row for row in report.rows}
    assert list(rows) == ["smaller_font", "labels_out", "radial"]
    assert {row.status for row in rows.values()} == {"not_computed"}
    search = raised.value.row_search
    assert search is not None
    radial = search.row("radial")
    assert radial.status == "found"
    assert radial.values == {"labels.circular.placement": "radial"}


@pytest.mark.circular
@pytest.mark.slow
def test_inner_font_row_is_the_largest_font_that_fits(tmp_path: Path) -> None:
    """INNER-3: P1 at 15 pt; with the local repair (LP-D57) the largest font that fits is 14 pt (inner labels stay in cyclic order since OV-1366).

    One point smaller (14 and 13 pt) while labels were laid out in point-mode boxes (OV-165)."""

    with pytest.raises(LabelsDoNotFitError) as raised:
        _draw(tmp_path, _both_args("NC_005856.1.gb", "--label_font_size", "15"))

    rows = {row.kind: row for row in raised.value.report.rows}
    font = rows["smaller_font"]
    if font.status == "not_computed":
        search = raised.value.row_search
        assert search is not None
        font = search.row("smaller_font")
    assert font.status == "found"
    assert font.values == {path: 14.0 for path in FONT_SIZE_PATHS}


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


@pytest.mark.circular
def test_the_fit_check_measures_the_drawn_definition_of_a_record_with_a_label(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """OV-707: the definition text the check and the repair measure carries the record label, as drawn."""

    measured: list[Aabb] = []
    drawn: list[Aabb] = []
    centre_definition_text = circular_assemble_module._centre_definition_text
    place_definition_group = circular_builders_module.place_definition_group_on_canvas

    def measure(definition: DefinitionGroup) -> tuple[Aabb, ...]:
        measured.append(definition.local_bounds)
        return centre_definition_text(definition)

    def place(group: Group, canvas_config: CircularCanvasConfigurator, *, position: str) -> Group:
        drawn.append(getattr(group, "_gbdraw_local_bounds"))
        return place_definition_group(group, canvas_config, position=position)

    monkeypatch.setattr(circular_assemble_module, "_centre_definition_text", measure)
    monkeypatch.setattr(circular_builders_module, "place_definition_group_on_canvas", place)
    record = SeqRecord(
        Seq("AAAACCGTGG" * 1000),
        id="labelled",
        annotations={"topology": "circular", "molecule_type": "DNA"},
    )
    record.features = [
        SeqFeature(SimpleLocation(2000, 3000, strand=1), type="CDS", qualifiers={"product": ["anchor enzyme"]})
    ]
    config = apply_config_overrides(
        None, {"labels.circular.scope": "outer", "canvas.show_gc": False, "canvas.show_skew": False}
    )
    request = CircularDiagramRequest(
        records=(
            RecordInput(
                InMemoryRecordSource(record),
                presentation=RecordPresentation(label="A record label wider than the accession line"),
            ),
        ),
        options=CircularDiagramOptions(config=config),
    )

    plan_request(request).build()

    assert measured and measured == drawn
