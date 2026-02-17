import math
from pathlib import Path
from typing import Any

import pytest
from Bio import SeqIO

import gbdraw.diagrams.circular.assemble as circular_assemble_module
import gbdraw.render.groups.circular.seq_record as circular_seq_record_group_module
import gbdraw.circular as circular_cli_module
from gbdraw.api.diagram import assemble_circular_diagram_from_record
from gbdraw.config.models import GbdrawConfig
from gbdraw.config.modify import modify_config_dict
from gbdraw.config.toml import load_config_toml
from gbdraw.features.colors import preprocess_color_tables
from gbdraw.features.factory import create_feature_dict
from gbdraw.io.colors import load_default_colors
from gbdraw.labels.filtering import preprocess_label_filtering
from gbdraw.tracks import parse_track_specs
from svgwrite import Drawing


SELECTED_FEATURES = ["CDS", "rRNA", "tRNA", "tmRNA", "ncRNA", "misc_RNA", "repeat_region"]


def _load_record():
    input_path = Path(__file__).parent / "test_inputs" / "HmmtDNA.gbk"
    return SeqIO.read(str(input_path), "genbank")


def _make_config_dict(*, show_labels: bool) -> dict:
    config_dict = load_config_toml("gbdraw.data", "config.toml")
    return modify_config_dict(
        config_dict,
        show_labels=show_labels,
        allow_inner_labels=False,
        resolve_overlaps=False,
        strandedness=True,
        track_type="tuckin",
    )


def test_feature_width_override_reaches_feature_drawer(monkeypatch: pytest.MonkeyPatch) -> None:
    record = _load_record()
    config_dict = _make_config_dict(show_labels=False)
    cfg = GbdrawConfig.from_dict(config_dict)

    captured_ratio_factors: list[float] = []

    def fake_draw(
        self,
        feature_object,
        group,
        total_length,
        radius,
        track_ratio,
        track_ratio_factor,
        track_type,
        strandedness,
        length_param,
    ):
        captured_ratio_factors.append(float(track_ratio_factor))
        return group

    monkeypatch.setattr(circular_seq_record_group_module.FeatureDrawer, "draw", fake_draw)

    assemble_circular_diagram_from_record(
        record,
        config_dict=config_dict,
        selected_features_set=SELECTED_FEATURES,
        legend="none",
        track_specs=["features@w=96px"],
    )

    expected_ratio_factor = 96.0 / (cfg.canvas.circular.radius * cfg.canvas.circular.track_ratio)
    assert captured_ratio_factors
    assert all(
        math.isclose(ratio_factor, expected_ratio_factor, rel_tol=1e-6, abs_tol=1e-6)
        for ratio_factor in captured_ratio_factors
    )


def test_feature_width_generates_auto_relayout_overrides(monkeypatch: pytest.MonkeyPatch) -> None:
    record = _load_record()
    config_dict = _make_config_dict(show_labels=True)
    cfg = GbdrawConfig.from_dict(config_dict)

    captured: dict[str, Any] = {}

    def fake_add_axis_group_on_canvas(canvas, canvas_config, config_dict, *, radius_override=None, cfg=None):
        captured["axis_radius"] = radius_override
        return canvas

    def fake_add_tick_group_on_canvas(canvas, gb_record, canvas_config, config_dict, *, radius_override=None, cfg=None):
        captured["ticks_radius"] = radius_override
        return canvas

    def fake_add_gc_content_group_on_canvas(
        canvas,
        gb_record,
        gc_df,
        canvas_config,
        gc_config,
        config_dict,
        *,
        track_width_override=None,
        norm_factor_override=None,
        cfg=None,
    ):
        captured["gc_width"] = track_width_override
        captured["gc_norm"] = norm_factor_override
        return canvas

    def fake_add_gc_skew_group_on_canvas(
        canvas,
        gb_record,
        gc_df,
        canvas_config,
        skew_config,
        config_dict,
        *,
        track_width_override=None,
        norm_factor_override=None,
        cfg=None,
    ):
        captured["skew_width"] = track_width_override
        captured["skew_norm"] = norm_factor_override
        return canvas

    def fake_add_labels_group_on_canvas(
        canvas,
        gb_record,
        canvas_config,
        feature_config,
        config_dict,
        *,
        outer_arena=None,
        cfg=None,
        precomputed_feature_dict=None,
        precalculated_labels=None,
        feature_track_ratio_factor_override=None,
    ):
        captured["labels_outer_arena"] = outer_arena
        captured["labels_feature_ratio"] = feature_track_ratio_factor_override
        return canvas

    def fake_add_record_group_on_canvas(
        canvas,
        record,
        canvas_config,
        feature_config,
        config_dict,
        *,
        cfg=None,
        precomputed_feature_dict=None,
        precalculated_labels=None,
        feature_track_ratio_factor_override=None,
    ):
        captured["record_feature_ratio"] = feature_track_ratio_factor_override
        return canvas

    monkeypatch.setattr(circular_assemble_module, "add_axis_group_on_canvas", fake_add_axis_group_on_canvas)
    monkeypatch.setattr(circular_assemble_module, "add_tick_group_on_canvas", fake_add_tick_group_on_canvas)
    monkeypatch.setattr(circular_assemble_module, "add_gc_content_group_on_canvas", fake_add_gc_content_group_on_canvas)
    monkeypatch.setattr(circular_assemble_module, "add_gc_skew_group_on_canvas", fake_add_gc_skew_group_on_canvas)
    monkeypatch.setattr(circular_assemble_module, "add_labels_group_on_canvas", fake_add_labels_group_on_canvas)
    monkeypatch.setattr(circular_assemble_module, "add_record_group_on_canvas", fake_add_record_group_on_canvas)
    monkeypatch.setattr(
        circular_assemble_module,
        "add_record_definition_group_on_canvas",
        lambda canvas, gb_record, canvas_config, species, strain, config_dict, *, cfg=None: canvas,
    )
    monkeypatch.setattr(
        circular_assemble_module,
        "add_legend_group_on_canvas",
        lambda canvas, canvas_config, legend_config, legend_table: canvas,
    )

    assemble_circular_diagram_from_record(
        record,
        config_dict=config_dict,
        selected_features_set=SELECTED_FEATURES,
        legend="none",
        track_specs=["features@w=120px"],
    )

    assert captured["record_feature_ratio"] is not None
    assert captured["labels_feature_ratio"] is not None
    assert captured["axis_radius"] is not None
    assert captured["ticks_radius"] is not None
    assert captured["gc_norm"] is not None
    assert captured["skew_norm"] is not None
    assert captured["labels_outer_arena"] is not None
    assert float(captured["labels_outer_arena"][1]) > float(captured["labels_outer_arena"][0])
    assert not math.isclose(float(captured["axis_radius"]), float(cfg.canvas.circular.radius), rel_tol=1e-6, abs_tol=1e-6)
    assert not math.isclose(float(captured["ticks_radius"]), float(cfg.canvas.circular.radius), rel_tol=1e-6, abs_tol=1e-6)


def test_explicit_track_placement_beats_auto_relayout(monkeypatch: pytest.MonkeyPatch) -> None:
    record = _load_record()
    config_dict = _make_config_dict(show_labels=True)
    cfg = GbdrawConfig.from_dict(config_dict)
    base_radius = float(cfg.canvas.circular.radius)
    captured: dict[str, Any] = {}

    def fake_add_axis_group_on_canvas(canvas, canvas_config, config_dict, *, radius_override=None, cfg=None):
        captured["axis_radius"] = radius_override
        return canvas

    def fake_add_tick_group_on_canvas(canvas, gb_record, canvas_config, config_dict, *, radius_override=None, cfg=None):
        captured["ticks_radius"] = radius_override
        return canvas

    def fake_add_gc_content_group_on_canvas(
        canvas,
        gb_record,
        gc_df,
        canvas_config,
        gc_config,
        config_dict,
        *,
        track_width_override=None,
        norm_factor_override=None,
        cfg=None,
    ):
        captured["gc_width"] = track_width_override
        captured["gc_norm"] = norm_factor_override
        return canvas

    def fake_add_gc_skew_group_on_canvas(
        canvas,
        gb_record,
        gc_df,
        canvas_config,
        skew_config,
        config_dict,
        *,
        track_width_override=None,
        norm_factor_override=None,
        cfg=None,
    ):
        captured["skew_width"] = track_width_override
        captured["skew_norm"] = norm_factor_override
        return canvas

    def fake_add_labels_group_on_canvas(
        canvas,
        gb_record,
        canvas_config,
        feature_config,
        config_dict,
        *,
        outer_arena=None,
        cfg=None,
        precomputed_feature_dict=None,
        precalculated_labels=None,
        feature_track_ratio_factor_override=None,
    ):
        captured["labels_outer_arena"] = outer_arena
        return canvas

    monkeypatch.setattr(circular_assemble_module, "add_axis_group_on_canvas", fake_add_axis_group_on_canvas)
    monkeypatch.setattr(circular_assemble_module, "add_tick_group_on_canvas", fake_add_tick_group_on_canvas)
    monkeypatch.setattr(circular_assemble_module, "add_gc_content_group_on_canvas", fake_add_gc_content_group_on_canvas)
    monkeypatch.setattr(circular_assemble_module, "add_gc_skew_group_on_canvas", fake_add_gc_skew_group_on_canvas)
    monkeypatch.setattr(circular_assemble_module, "add_labels_group_on_canvas", fake_add_labels_group_on_canvas)
    monkeypatch.setattr(circular_assemble_module, "add_record_group_on_canvas", lambda *args, **kwargs: args[0])
    monkeypatch.setattr(
        circular_assemble_module,
        "add_record_definition_group_on_canvas",
        lambda canvas, gb_record, canvas_config, species, strain, config_dict, *, cfg=None: canvas,
    )
    monkeypatch.setattr(
        circular_assemble_module,
        "add_legend_group_on_canvas",
        lambda canvas, canvas_config, legend_config, legend_table: canvas,
    )

    assemble_circular_diagram_from_record(
        record,
        config_dict=config_dict,
        selected_features_set=SELECTED_FEATURES,
        legend="none",
        track_specs=[
            "features@w=120px",
            "axis@r=0.80",
            "ticks@r=0.83",
            "gc_content@r=0.74,w=22px",
            "gc_skew@r=0.68,w=18px",
            "labels@r=1.28,w=44px",
        ],
    )

    assert math.isclose(float(captured["axis_radius"]), 0.80 * base_radius, rel_tol=1e-6, abs_tol=1e-6)
    assert math.isclose(float(captured["ticks_radius"]), 0.83 * base_radius, rel_tol=1e-6, abs_tol=1e-6)
    assert math.isclose(float(captured["gc_norm"]), 0.74, rel_tol=1e-6, abs_tol=1e-6)
    assert math.isclose(float(captured["skew_norm"]), 0.68, rel_tol=1e-6, abs_tol=1e-6)
    assert math.isclose(float(captured["gc_width"]), 22.0, rel_tol=1e-6, abs_tol=1e-6)
    assert math.isclose(float(captured["skew_width"]), 18.0, rel_tol=1e-6, abs_tol=1e-6)

    expected_center = 1.28 * base_radius
    expected_half_width = 22.0
    assert captured["labels_outer_arena"] is not None
    assert math.isclose(float(captured["labels_outer_arena"][0]), expected_center - expected_half_width, rel_tol=1e-6, abs_tol=1e-6)
    assert math.isclose(float(captured["labels_outer_arena"][1]), expected_center + expected_half_width, rel_tol=1e-6, abs_tol=1e-6)


@pytest.mark.parametrize(
    "spec",
    [
        "features@r=0.82,w=48px",
        "features@ri=0.70,w=48px",
        "features@ro=0.94,w=48px",
    ],
)
def test_features_center_placement_warns_and_is_ignored(spec: str, caplog: pytest.LogCaptureFixture) -> None:
    caplog.clear()
    ts = parse_track_specs([spec], mode="circular")[0]
    with caplog.at_level("WARNING"):
        ratio = circular_assemble_module._resolve_feature_track_ratio_factor_override(
            ts,
            base_radius_px=400.0,
            base_track_ratio=0.2,
        )
    assert ratio == pytest.approx(0.6)
    assert any("supports width only" in message for message in caplog.messages)


def test_cli_feature_width_must_be_positive() -> None:
    with pytest.raises(SystemExit):
        circular_cli_module._get_args(["--gbk", "dummy.gb", "--feature_width", "0"])
    with pytest.raises(SystemExit):
        circular_cli_module._get_args(["--gbk", "dummy.gb", "--feature_width", "-10"])


def test_cli_feature_width_forwards_internal_feature_track_spec(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    record = _load_record()
    captured: dict[str, Any] = {}

    monkeypatch.setattr(circular_cli_module, "load_gbks", lambda paths, mode: [record])
    monkeypatch.setattr(circular_cli_module, "read_color_table", lambda _path: None)
    monkeypatch.setattr(circular_cli_module, "load_default_colors", lambda _path, _palette: None)
    monkeypatch.setattr(circular_cli_module, "save_figure", lambda canvas, formats: None)

    def fake_assemble(*args, **kwargs):
        captured["track_specs"] = kwargs.get("track_specs")
        return Drawing(filename=str(tmp_path / "dummy.svg"))

    monkeypatch.setattr(circular_cli_module, "assemble_circular_diagram_from_record", fake_assemble)

    circular_cli_module.circular_main(
        ["--gbk", "dummy.gb", "--feature_width", "42", "--format", "svg", "-o", str(tmp_path / "out")]
    )

    assert captured["track_specs"] == ["features@w=42px"]


def test_radius_mapper_uses_primary_feature_track_for_auto_relayout_with_resolve_overlaps() -> None:
    record = _load_record()
    config_dict = _make_config_dict(show_labels=True)
    config_dict = modify_config_dict(
        config_dict,
        resolve_overlaps=True,
        strandedness=False,
        track_type="middle",
    )
    cfg = GbdrawConfig.from_dict(config_dict)

    color_table, default_colors = preprocess_color_tables(None, load_default_colors("", "default"))
    label_filtering = preprocess_label_filtering(cfg.labels.filtering.as_dict())
    feature_dict, _ = create_feature_dict(
        record,
        color_table,
        SELECTED_FEATURES,
        default_colors,
        cfg.canvas.strandedness,
        cfg.canvas.resolve_overlaps,
        label_filtering,
    )

    class _CanvasStub:
        radius = cfg.canvas.circular.radius
        track_ratio = cfg.canvas.circular.track_ratio
        length_param = "short"

    feature_ratio_override = 75.0 / (float(_CanvasStub.radius) * float(_CanvasStub.track_ratio))
    _, old_band, new_band = circular_assemble_module._build_feature_radius_mapper(
        feature_dict,
        len(record.seq),
        canvas_config=_CanvasStub,
        cfg=cfg,
        feature_track_ratio_factor_override=feature_ratio_override,
    )

    assert old_band is not None
    assert new_band is not None

    default_ratio_factor = float(cfg.canvas.circular.track_ratio_factors["short"][0])
    default_cds_ratio = float(cfg.canvas.circular.track_ratio) * default_ratio_factor
    expected_old_outer = float(cfg.canvas.circular.radius) * (1.0 + 0.5 * default_cds_ratio)
    assert math.isclose(float(old_band[1]), expected_old_outer, rel_tol=1e-6, abs_tol=1e-6)


def test_feature_width_expand_canvas_when_labels_hidden_with_resolve_overlaps() -> None:
    record = _load_record()
    config_dict = _make_config_dict(show_labels=False)
    config_dict = modify_config_dict(
        config_dict,
        resolve_overlaps=True,
        strandedness=False,
        track_type="middle",
    )
    cfg = GbdrawConfig.from_dict(config_dict)
    expected_base_height = float(cfg.canvas.circular.height)

    canvas = assemble_circular_diagram_from_record(
        record,
        config_dict=config_dict,
        selected_features_set=SELECTED_FEATURES,
        legend="left",
        track_specs=["features@w=75px"],
    )

    rendered_height = float(str(canvas.attribs["height"]).rstrip("px"))
    assert rendered_height > expected_base_height


def test_auto_relayout_core_tracks_are_stable_across_show_labels_toggle() -> None:
    record = _load_record()
    captured_by_show: dict[bool, dict[str, Any]] = {}

    for show_labels in [False, True]:
        config_dict = _make_config_dict(show_labels=show_labels)
        config_dict = modify_config_dict(
            config_dict,
            resolve_overlaps=True,
            strandedness=False,
            track_type="middle",
        )
        captured: dict[str, Any] = {}

        def fake_add_axis_group_on_canvas(canvas, canvas_config, config_dict, *, radius_override=None, cfg=None):
            captured["axis"] = radius_override
            return canvas

        def fake_add_tick_group_on_canvas(canvas, gb_record, canvas_config, config_dict, *, radius_override=None, cfg=None):
            captured["ticks"] = radius_override
            return canvas

        def fake_add_gc_content_group_on_canvas(
            canvas,
            gb_record,
            gc_df,
            canvas_config,
            gc_config,
            config_dict,
            *,
            track_width_override=None,
            norm_factor_override=None,
            cfg=None,
        ):
            captured["gc_norm"] = norm_factor_override
            return canvas

        def fake_add_gc_skew_group_on_canvas(
            canvas,
            gb_record,
            gc_df,
            canvas_config,
            skew_config,
            config_dict,
            *,
            track_width_override=None,
            norm_factor_override=None,
            cfg=None,
        ):
            captured["skew_norm"] = norm_factor_override
            return canvas

        def fake_add_labels_group_on_canvas(
            canvas,
            gb_record,
            canvas_config,
            feature_config,
            config_dict,
            *,
            outer_arena=None,
            cfg=None,
            precomputed_feature_dict=None,
            precalculated_labels=None,
            feature_track_ratio_factor_override=None,
        ):
            captured["outer_arena"] = outer_arena
            return canvas

        original_axis = circular_assemble_module.add_axis_group_on_canvas
        original_ticks = circular_assemble_module.add_tick_group_on_canvas
        original_gc = circular_assemble_module.add_gc_content_group_on_canvas
        original_skew = circular_assemble_module.add_gc_skew_group_on_canvas
        original_labels = circular_assemble_module.add_labels_group_on_canvas
        original_record = circular_assemble_module.add_record_group_on_canvas
        original_definition = circular_assemble_module.add_record_definition_group_on_canvas
        original_legend = circular_assemble_module.add_legend_group_on_canvas
        circular_assemble_module.add_axis_group_on_canvas = fake_add_axis_group_on_canvas
        circular_assemble_module.add_tick_group_on_canvas = fake_add_tick_group_on_canvas
        circular_assemble_module.add_gc_content_group_on_canvas = fake_add_gc_content_group_on_canvas
        circular_assemble_module.add_gc_skew_group_on_canvas = fake_add_gc_skew_group_on_canvas
        circular_assemble_module.add_labels_group_on_canvas = fake_add_labels_group_on_canvas
        circular_assemble_module.add_record_group_on_canvas = lambda canvas, *args, **kwargs: canvas
        circular_assemble_module.add_record_definition_group_on_canvas = (
            lambda canvas, gb_record, canvas_config, species, strain, config_dict, *, cfg=None: canvas
        )
        circular_assemble_module.add_legend_group_on_canvas = (
            lambda canvas, canvas_config, legend_config, legend_table: canvas
        )
        try:
            assemble_circular_diagram_from_record(
                record,
                config_dict=config_dict,
                selected_features_set=SELECTED_FEATURES,
                legend="left",
                track_specs=["features@w=75px"],
            )
        finally:
            circular_assemble_module.add_axis_group_on_canvas = original_axis
            circular_assemble_module.add_tick_group_on_canvas = original_ticks
            circular_assemble_module.add_gc_content_group_on_canvas = original_gc
            circular_assemble_module.add_gc_skew_group_on_canvas = original_skew
            circular_assemble_module.add_labels_group_on_canvas = original_labels
            circular_assemble_module.add_record_group_on_canvas = original_record
            circular_assemble_module.add_record_definition_group_on_canvas = original_definition
            circular_assemble_module.add_legend_group_on_canvas = original_legend

        captured_by_show[show_labels] = captured

    off = captured_by_show[False]
    on = captured_by_show[True]
    assert off.get("outer_arena") is None
    assert on.get("outer_arena") is not None
    assert math.isclose(float(off["axis"]), float(on["axis"]), rel_tol=1e-9, abs_tol=1e-9)
    assert math.isclose(float(off["ticks"]), float(on["ticks"]), rel_tol=1e-9, abs_tol=1e-9)
    assert math.isclose(float(off["gc_norm"]), float(on["gc_norm"]), rel_tol=1e-9, abs_tol=1e-9)
    assert math.isclose(float(off["skew_norm"]), float(on["skew_norm"]), rel_tol=1e-9, abs_tol=1e-9)
