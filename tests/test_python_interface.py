from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from svgwrite import Drawing

import gbdraw
import gbdraw.api.request_render as request_render_module
import gbdraw.interface as interface
from gbdraw.api.options import (
    CircularDiagramOptions,
    CircularOutputOptions,
    CircularRequestTrackOptions,
    LinearDiagramOptions,
    LinearOutputOptions,
    LinearRequestTrackOptions,
)
from gbdraw.exceptions import ExportError, ValidationError


def _record(record_id: str = "record") -> SeqRecord:
    return SeqRecord(Seq("ATGC" * 25), id=record_id)


def test_root_namespace_is_the_small_beginner_facing_api() -> None:
    assert gbdraw.__all__ == [
        "CircularLayout",
        "CircularOptions",
        "CircularTrackOptions",
        "ComparisonRingOptions",
        "ComparisonRingTrackOptions",
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


def test_root_and_typed_track_option_tiers_have_distinct_canonical_names() -> None:
    assert gbdraw.CircularTrackOptions is interface.CircularTrackOptions
    assert gbdraw.LinearTrackOptions is interface.LinearTrackOptions
    assert gbdraw.CircularTrackOptions is not CircularRequestTrackOptions
    assert gbdraw.LinearTrackOptions is not LinearRequestTrackOptions


def test_comparison_ring_names_are_canonical_with_conservation_aliases() -> None:
    assert gbdraw.ComparisonRingOptions is interface.ComparisonRingOptions
    assert (
        gbdraw.ComparisonRingTrackOptions
        is interface.ComparisonRingTrackOptions
    )
    assert gbdraw.ConservationOptions is gbdraw.ComparisonRingOptions
    assert (
        gbdraw.ConservationTrackOptions
        is gbdraw.ComparisonRingTrackOptions
    )

    legacy = gbdraw.ConservationOptions(
        tracks=(gbdraw.ConservationTrackOptions(source="hits.tsv"),),
    )

    assert type(legacy).__name__ == "ComparisonRingOptions"
    assert type(legacy.tracks[0]).__name__ == "ComparisonRingTrackOptions"
    assert isinstance(legacy, gbdraw.ComparisonRingOptions)
    assert isinstance(legacy.tracks[0], gbdraw.ComparisonRingTrackOptions)


def test_circular_options_names_comparison_ring_type_in_validation() -> None:
    with pytest.raises(ValidationError, match="ComparisonRingOptions"):
        interface.CircularOptions(comparison_rings=object())  # type: ignore[arg-type]


def test_circular_options_exposes_comparison_rings_with_conservation_alias() -> None:
    canonical = interface.ComparisonRingOptions(
        tracks=(interface.ComparisonRingTrackOptions(source="hits.tsv"),),
    )

    current = interface.CircularOptions(comparison_rings=canonical)
    legacy = interface.CircularOptions(conservation=canonical)

    assert current.comparison_rings is canonical
    assert current.conservation is canonical
    assert legacy.comparison_rings is canonical
    assert legacy.conservation is canonical
    with pytest.raises(ValidationError, match="not both"):
        interface.CircularOptions(
            comparison_rings=canonical,
            conservation=canonical,
        )


def test_circular_options_alias_survives_dataclass_replace() -> None:
    original = interface.CircularOptions(
        comparison_rings=interface.ComparisonRingOptions(reference="query"),
    )
    replacement = interface.ComparisonRingOptions(reference="subject")

    renamed = replace(original, species="Example species")
    replaced_rings = replace(original, comparison_rings=replacement)

    assert renamed.species == "Example species"
    assert renamed.comparison_rings is original.comparison_rings
    assert renamed.conservation is original.comparison_rings
    assert replaced_rings.comparison_rings is replacement
    assert replaced_rings.conservation is replacement
    with pytest.raises(ValidationError, match="not both"):
        replace(original, conservation=replacement)


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

    monkeypatch.setattr(request_render_module, "build_circular_diagram", fake_single)
    monkeypatch.setattr(request_render_module, "build_circular_multi_diagram", fake_multi)
    monkeypatch.setattr(interface, "_interactive_context", lambda *_args, **_kwargs: None)

    one = _record("one")
    assert interface.draw_circular(one).mode == "circular"
    assert calls[0][0] == "single"
    assert calls[0][1].id == "one"
    assert isinstance(calls[0][2], CircularDiagramOptions)

    records = [one, _record("two")]
    diagram = interface.draw_circular(records)
    assert diagram.records == tuple(records)
    assert calls[1][0] == "multi"
    default_layout = calls[1][2]
    assert default_layout.multi_record_size_mode == "auto"
    assert default_layout.multi_record_positions is None

    interface.draw_circular(
        records,
        layout=interface.CircularLayout(size="equal", positions=("#1@1", "#2@1")),
    )
    assert calls[2][0] == "multi"
    legacy_layout = calls[2][2]
    assert legacy_layout.multi_record_size_mode == "equal"
    assert legacy_layout.multi_record_positions == ("#1@1", "#2@1")


def test_root_api_builds_metadata_only_for_explicit_interactive_render(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        request_render_module,
        "build_circular_diagram",
        lambda *_args, **_kwargs: Drawing("circular.svg"),
    )
    calls = 0

    def fail_context(*_args, **_kwargs):
        nonlocal calls
        calls += 1
        raise RuntimeError("metadata exploded")

    monkeypatch.setattr(interface, "_interactive_context", fail_context)

    diagram = interface.draw_circular(_record())

    assert calls == 0
    assert diagram.to_svg().startswith("<svg")
    assert calls == 0
    with pytest.raises(
        ExportError,
        match="Interactive SVG metadata generation failed: metadata exploded",
    ):
        diagram.to_svg(interactive=True)
    assert calls == 1


def test_draw_circular_one_record_layout_reaches_grid_builder(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, object] = {}

    def fake_multi(records, *, options, layout):
        captured["records"] = tuple(records)
        captured["options"] = options
        captured["layout"] = layout
        return Drawing("grid.svg")

    monkeypatch.setattr(request_render_module, "build_circular_multi_diagram", fake_multi)
    monkeypatch.setattr(interface, "_interactive_context", lambda *_args, **_kwargs: None)

    diagram = interface.draw_circular(
        _record("one"),
        layout=interface.CircularLayout(size="equal", positions=("#1@1",)),
    )

    assert diagram.mode == "circular"
    assert [record.id for record in captured["records"]] == ["one"]
    assert isinstance(captured["options"], CircularDiagramOptions)
    assert captured["layout"].multi_record_size_mode == "equal"
    assert captured["layout"].multi_record_positions == ("#1@1",)


@pytest.mark.circular
def test_draw_circular_one_record_layout_real_render() -> None:
    record = _record("one")
    record.annotations["molecule_type"] = "DNA"
    record.annotations["topology"] = "circular"

    diagram = interface.draw_circular(
        record,
        layout=interface.CircularLayout(size="equal", positions=("#1@1",)),
    )

    assert diagram.mode == "circular"
    assert diagram.to_svg().startswith("<svg")


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

    monkeypatch.setattr(request_render_module, "build_circular_diagram", fake_single)
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
    assert isinstance(options, CircularDiagramOptions)
    assert isinstance(options.tracks, CircularRequestTrackOptions)
    assert isinstance(options.output, CircularOutputOptions)
    assert options.selected_features_set == ("CDS",)
    assert options.colors.color_table is color_table
    assert options.feature_visibility_table_file == "visibility.tsv"
    assert options.label_whitelist_file == "labels.tsv"
    assert options.plot_title == "Genome"
    assert options.output.plot_title_position == "top"
    assert options.depth_track_files == [["depth.tsv"]]
    assert options.depth_track_labels == ["Coverage"]


def test_draw_linear_routes_through_typed_request_plan(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    captured: dict[str, object] = {}

    def fake_linear(records, *, options, layout, **_kwargs):
        captured["records"] = tuple(records)
        captured["options"] = options
        captured["layout"] = layout
        return Drawing("linear.svg")

    monkeypatch.setattr(request_render_module, "build_linear_diagram", fake_linear)
    monkeypatch.setattr(interface, "_interactive_context", lambda *_args, **_kwargs: None)

    diagram = interface.draw_linear(
        [_record("one"), _record("two")],
        options=interface.LinearOptions(
            title=interface.TitleOptions(text="Linear", position="top"),
        ),
        layout=interface.LinearLayout(record_gap=30, positions=("#1@1", "#2@2")),
    )

    assert diagram.mode == "linear"
    assert [record.id for record in captured["records"]] == ["one", "two"]
    compiled = captured["options"]
    assert isinstance(compiled, LinearDiagramOptions)
    assert isinstance(compiled.tracks, LinearRequestTrackOptions)
    assert isinstance(compiled.output, LinearOutputOptions)
    assert compiled.output.plot_title_position == "top"
    assert captured["layout"].record_gap_px == 30
    assert captured["layout"].multi_record_positions == ("#1@1", "#2@2")


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
        comparison_rings=interface.ComparisonRingOptions(
            tracks=(
                interface.ComparisonRingTrackOptions(
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


def test_root_api_rejects_wrong_mode_config_override_paths() -> None:
    with pytest.raises(
        ValidationError,
        match="Circular config overrides cannot target Linear settings",
    ):
        interface.draw_circular(
            _record(),
            options=interface.CircularOptions(
                config_overrides={"canvas.linear.track_layout": "above"},
            ),
        )
    with pytest.raises(
        ValidationError,
        match="Linear config overrides cannot target Circular settings",
    ):
        interface.draw_linear(
            _record(),
            options=interface.LinearOptions(
                config_overrides={"canvas.circular.track_type": "middle"},
            ),
        )


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
