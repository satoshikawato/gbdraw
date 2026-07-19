from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from svgwrite import Drawing

import gbdraw
import gbdraw.interface as interface
from gbdraw.exceptions import ValidationError


def _record(record_id: str = "record") -> SeqRecord:
    return SeqRecord(Seq("ATGC" * 25), id=record_id)


def test_root_namespace_is_the_small_beginner_facing_api() -> None:
    assert gbdraw.__all__ == [
        "CircularLayout",
        "CircularOptions",
        "CircularTrackOptions",
        "ConservationOptions",
        "ConservationTrackOptions",
        "DepthTrackOptions",
        "Diagram",
        "FeatureOptions",
        "LabelOptions",
        "LinearComparisonOptions",
        "LinearLayout",
        "LinearOptions",
        "LinearTrackOptions",
        "Thresholds",
        "TitleOptions",
        "__version__",
        "draw_circular",
        "draw_linear",
        "read_genbank",
        "read_gff",
    ]


def test_draw_circular_dispatches_from_record_count(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[tuple[str, object, object]] = []

    def fake_single(record, *, options):
        calls.append(("single", record, options))
        return Drawing("single.svg")

    def fake_multi(records, *, options, layout):
        calls.append(("multi", tuple(records), layout))
        return Drawing("multi.svg")

    monkeypatch.setattr(interface, "_build_circular_diagram", fake_single)
    monkeypatch.setattr(interface, "_build_circular_multi_diagram", fake_multi)
    monkeypatch.setattr(interface, "_interactive_context", lambda *_args, **_kwargs: None)

    one = _record("one")
    assert interface.draw_circular(one).mode == "circular"
    assert calls[0][0:2] == ("single", one)

    records = [one, _record("two")]
    diagram = interface.draw_circular(
        records,
        layout=interface.CircularLayout(size="equal", positions=("#1@1", "#2@1")),
    )
    assert diagram.records == tuple(records)
    assert calls[1][0] == "multi"
    legacy_layout = calls[1][2]
    assert legacy_layout.multi_record_size_mode == "equal"
    assert legacy_layout.multi_record_positions == ("#1@1", "#2@1")


def test_draw_circular_rejects_grid_layout_for_one_record() -> None:
    with pytest.raises(ValidationError, match="at least two records"):
        interface.draw_circular(_record(), layout=interface.CircularLayout())


def test_grouped_options_compile_to_the_existing_render_engine(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, object] = {}
    color_table = pd.DataFrame(
        [["CDS", "product", "polymerase", "#123456"]],
        columns=["feature_type", "qualifier", "value", "color"],
    )

    def fake_single(_record, *, options):
        captured["options"] = options
        return Drawing("out.svg")

    monkeypatch.setattr(interface, "_build_circular_diagram", fake_single)
    monkeypatch.setattr(interface, "_interactive_context", lambda *_args, **_kwargs: None)

    interface.draw_circular(
        _record(),
        options=interface.CircularOptions(
            features=interface.FeatureOptions(
                types=("CDS",),
                color_table=color_table,
                visibility=Path("visibility.tsv"),
            ),
            labels=interface.LabelOptions(whitelist=Path("labels.tsv")),
            title=interface.TitleOptions(text="Genome", position="top", font_size=18),
            depth_tracks=(
                interface.DepthTrackOptions(
                    source=Path("depth.tsv"),
                    label="Coverage",
                    color="#336699",
                ),
            ),
        ),
    )

    options = captured["options"]
    assert options.selected_features_set == ("CDS",)
    assert options.colors.color_table is color_table
    assert options.feature_visibility_table_file == "visibility.tsv"
    assert options.label_whitelist_file == "labels.tsv"
    assert options.plot_title == "Genome"
    assert options.output.plot_title_position == "top"
    assert options.depth_track_files == [["depth.tsv"]]
    assert options.depth_track_labels == ["Coverage"]


def test_circular_companion_sequence_reaches_interactive_context() -> None:
    reference = _record("reference")
    comparison = _record("comparison")
    blast = pd.DataFrame(
        [["reference", "comparison", 99.0, 20, 0, 0, 1, 20, 20, 1, 1e-20, 50]],
        columns=[
            "query", "subject", "identity", "alignment_length", "mismatches",
            "gap_opens", "qstart", "qend", "sstart", "send", "evalue", "bitscore",
        ],
    )
    options = interface.CircularOptions(
        conservation=interface.ConservationOptions(
            tracks=(
                interface.ConservationTrackOptions(
                    source=blast,
                    comparison_sequence_source=(comparison,),
                ),
            )
        )
    )

    context = interface._interactive_context(
        [reference],
        options=options,
        legacy=interface._circular_options(options, record_count=1),
        mode="circular",
    )

    assert [(source["origin"], source["recordId"]) for source in context.sequence_sources] == [
        ("circular-reference", "reference"),
        ("homology-comparison", "comparison"),
    ]


def test_draw_functions_require_mode_specific_options() -> None:
    with pytest.raises(ValidationError, match="CircularOptions"):
        interface.draw_circular(_record(), options=interface.LinearOptions())  # type: ignore[arg-type]
    with pytest.raises(ValidationError, match="LinearOptions"):
        interface.draw_linear(_record(), options=interface.CircularOptions())  # type: ignore[arg-type]


def test_diagram_save_writes_exactly_the_requested_file(tmp_path: Path) -> None:
    drawing = Drawing("internal-name.svg")
    diagram = interface.Diagram(
        drawing,
        mode="circular",
        records=(_record(),),
        interactive_context=None,
    )
    output = tmp_path / "chosen-name.svg"

    assert diagram.save(output) == output
    assert output.read_text(encoding="utf-8").startswith("<svg")
    assert list(tmp_path.iterdir()) == [output]

    with pytest.raises(ValidationError, match="already exists"):
        diagram.save(output)


def test_read_genbank_accepts_one_path(examples_dir: Path) -> None:
    records = interface.read_genbank(examples_dir / "MjeNMV.gb")
    assert [record.id for record in records] == ["LC738868.1"]


@pytest.mark.circular
def test_beginner_circular_workflow_runs(examples_dir: Path, tmp_path: Path) -> None:
    record = interface.read_genbank(examples_dir / "MjeNMV.gb")[0]
    diagram = interface.draw_circular(
        record,
        options=interface.CircularOptions(
            features=interface.FeatureOptions(types=("CDS",)),
            title=interface.TitleOptions(text="Example genome", position="top"),
        ),
    )

    assert diagram.to_svg().startswith("<svg")
    assert diagram.save(tmp_path / "diagram.svg").is_file()
