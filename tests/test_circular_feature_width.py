import math
import re
from pathlib import Path
from typing import Any

import pytest
from Bio import SeqIO

import gbdraw.diagrams.circular.assemble as circular_assemble_module
import gbdraw.render.groups.circular.seq_record as circular_seq_record_group_module
import gbdraw.circular as circular_cli_module
import gbdraw.labels.circular as circular_labels_module
from gbdraw.api.diagram import assemble_circular_diagram_from_record
from gbdraw.config.models import GbdrawConfig
from gbdraw.config.modify import modify_config_dict
from gbdraw.config.toml import load_config_toml
from gbdraw.features.colors import preprocess_color_tables
from gbdraw.features.factory import create_feature_dict
from gbdraw.io.colors import load_default_colors
from gbdraw.labels.filtering import preprocess_label_filtering
from gbdraw.svg.circular_ticks import get_circular_tick_path_ratio_bounds
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


def _label_genome_intervals_for_clearance(label: dict, total_length: int) -> list[tuple[float, float]]:
    if total_length <= 0:
        return []
    total_length_f = float(total_length)
    start_x = float(label.get("start_x", 0.0))
    start_y = float(label.get("start_y", 0.0))
    radius = math.hypot(start_x, start_y)
    if radius <= 1e-9:
        middle = float(label.get("middle", 0.0)) % total_length_f
        return [(middle, middle)]
    width_px = max(0.0, float(label.get("width_px", 0.0)))
    if width_px <= 1e-9:
        anchor_angle = (math.degrees(math.atan2(start_y, start_x)) + 90.0) % 360.0
        anchor_position = (anchor_angle / 360.0) * total_length_f
        return [(anchor_position, anchor_position)]
    span_bp = total_length_f * (width_px / (2.0 * math.pi * radius))
    if span_bp >= total_length_f:
        return [(0.0, total_length_f)]
    anchor_angle = (math.degrees(math.atan2(start_y, start_x)) + 90.0) % 360.0
    center_position = (anchor_angle / 360.0) * total_length_f
    half_span = 0.5 * span_bp
    start_pos = center_position - half_span
    end_pos = center_position + half_span
    if start_pos < 0.0:
        return [(start_pos + total_length_f, total_length_f), (0.0, end_pos)]
    if end_pos > total_length_f:
        return [(start_pos, total_length_f), (0.0, end_pos - total_length_f)]
    return [(start_pos, end_pos)]


def _max_outer_feature_radius_for_label(
    intervals: list[tuple[float, float, float]],
    label: dict,
    total_length: int,
) -> float:
    label_intervals = _label_genome_intervals_for_clearance(label, total_length)
    max_outer = 0.0
    for label_start, label_end in label_intervals:
        if label_end <= label_start + 1e-9:
            max_outer = max(
                max_outer,
                circular_labels_module._max_outer_feature_radius_at_position(intervals, label_start, total_length),
            )
            continue
        for feature_start, feature_end, feature_outer in intervals:
            if float(feature_start) < float(label_end) and float(feature_end) > float(label_start):
                max_outer = max(max_outer, float(feature_outer))
    return max_outer


def _annulus_overlaps_band(
    annulus: tuple[float, float],
    band: tuple[float, float],
    tol: float = 1e-6,
) -> bool:
    annulus_inner, annulus_outer = sorted((float(annulus[0]), float(annulus[1])))
    band_inner, band_outer = sorted((float(band[0]), float(band[1])))
    return annulus_inner < (band_outer - tol) and annulus_outer > (band_inner + tol)


def _extract_group_translate(svg_text: str, group_id: str) -> tuple[float, float] | None:
    pattern = (
        rf'<g id="{re.escape(group_id)}"[^>]*\btransform="translate\(\s*'
        r'([-+0-9.eE]+)\s*,\s*([-+0-9.eE]+)\s*\)"'
    )
    match = re.search(pattern, svg_text)
    if match is None:
        return None
    return float(match.group(1)), float(match.group(2))


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


def test_feature_width_75_keeps_outer_labels_outside_local_feature_tracks() -> None:
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

    base_radius = float(cfg.canvas.circular.radius)
    base_track_ratio = float(cfg.canvas.circular.track_ratio)
    feature_ratio_override = 75.0 / (base_radius * base_track_ratio)
    canvas_stub = type(
        "_CanvasStub",
        (),
        {
            "radius": base_radius,
            "track_ratio": base_track_ratio,
            "length_param": "short",
        },
    )()
    radius_mapper, _, _ = circular_assemble_module._build_feature_radius_mapper(
        feature_dict,
        len(record.seq),
        canvas_config=canvas_stub,
        cfg=cfg,
        feature_track_ratio_factor_override=feature_ratio_override,
    )
    assert radius_mapper is not None

    default_anchor, default_arc_outer = circular_assemble_module._default_outer_label_arena(
        canvas_config=canvas_stub,
        cfg=cfg,
    )
    outer_arena = (float(radius_mapper(default_anchor)), float(radius_mapper(default_arc_outer)))

    labels = circular_labels_module.prepare_label_list(
        feature_dict,
        len(record.seq),
        base_radius,
        base_track_ratio,
        config_dict,
        cfg=cfg,
        outer_arena=outer_arena,
        feature_track_ratio_factor_override=feature_ratio_override,
    )
    outer_labels = [label for label in labels if not label.get("is_embedded")]
    assert outer_labels

    feature_intervals = circular_labels_module._build_outer_feature_radius_intervals(
        feature_dict,
        len(record.seq),
        base_radius,
        base_track_ratio,
        cfg,
        feature_track_ratio_factor_override=feature_ratio_override,
    )
    assert feature_intervals

    for label in outer_labels:
        local_outer = _max_outer_feature_radius_for_label(feature_intervals, label, len(record.seq))
        required_radius = local_outer + circular_labels_module.MIN_OUTER_LABEL_TEXT_CLEARANCE_PX
        label_bbox_min_radius = circular_labels_module._label_bbox_min_radius(
            label,
            len(record.seq),
            minimum_margin=0.0,
        )
        assert label_bbox_min_radius >= required_radius - 1e-3


def test_feature_width_75_auto_repositions_ticks_outside_feature_band_when_overlapping(monkeypatch: pytest.MonkeyPatch) -> None:
    record = _load_record()
    config_dict = _make_config_dict(show_labels=False)
    config_dict = modify_config_dict(
        config_dict,
        resolve_overlaps=True,
        strandedness=False,
        track_type="middle",
    )
    cfg = GbdrawConfig.from_dict(config_dict)
    base_radius = float(cfg.canvas.circular.radius)
    base_track_ratio = float(cfg.canvas.circular.track_ratio)
    feature_ratio_override = 75.0 / (base_radius * base_track_ratio)
    default_ratio_factor = float(cfg.canvas.circular.track_ratio_factors["short"][0])

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
    default_band = circular_assemble_module._compute_feature_band_bounds_px(
        feature_dict,
        len(record.seq),
        base_radius_px=base_radius,
        track_ratio=base_track_ratio,
        length_param="short",
        track_ratio_factor=default_ratio_factor,
        cfg=cfg,
        track_id_whitelist={0},
    )
    widened_band = circular_assemble_module._compute_feature_band_bounds_px(
        feature_dict,
        len(record.seq),
        base_radius_px=base_radius,
        track_ratio=base_track_ratio,
        length_param="short",
        track_ratio_factor=feature_ratio_override,
        cfg=cfg,
        track_id_whitelist={0},
    )
    assert default_band is not None
    assert widened_band is not None

    captured: dict[str, Any] = {}

    def fake_add_axis_group_on_canvas(canvas, canvas_config, config_dict, *, radius_override=None, cfg=None):
        captured["axis"] = radius_override
        return canvas

    def fake_add_tick_group_on_canvas(canvas, gb_record, canvas_config, config_dict, *, radius_override=None, cfg=None):
        captured["ticks"] = radius_override
        return canvas

    monkeypatch.setattr(circular_assemble_module, "add_axis_group_on_canvas", fake_add_axis_group_on_canvas)
    monkeypatch.setattr(circular_assemble_module, "add_tick_group_on_canvas", fake_add_tick_group_on_canvas)
    monkeypatch.setattr(circular_assemble_module, "add_gc_content_group_on_canvas", lambda canvas, *args, **kwargs: canvas)
    monkeypatch.setattr(circular_assemble_module, "add_gc_skew_group_on_canvas", lambda canvas, *args, **kwargs: canvas)
    monkeypatch.setattr(circular_assemble_module, "add_labels_group_on_canvas", lambda canvas, *args, **kwargs: canvas)
    monkeypatch.setattr(circular_assemble_module, "add_record_group_on_canvas", lambda canvas, *args, **kwargs: canvas)
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
        track_specs=["features@w=75px"],
    )

    ticks_center = float(captured.get("ticks") if captured.get("ticks") is not None else base_radius)
    tick_min_ratio, tick_max_ratio = get_circular_tick_path_ratio_bounds(
        len(record.seq),
        str(cfg.canvas.circular.track_type),
        bool(cfg.canvas.strandedness),
    )
    tick_annulus = (ticks_center * tick_min_ratio, ticks_center * tick_max_ratio)
    assert not _annulus_overlaps_band(tick_annulus, widened_band)


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


@pytest.mark.parametrize("track_type", ["tuckin", "middle", "spreadout"])
def test_feature_width_keeps_axis_concentric_with_rendered_tracks(track_type: str) -> None:
    record = _load_record()
    config_dict = _make_config_dict(show_labels=True)
    config_dict = modify_config_dict(
        config_dict,
        resolve_overlaps=False,
        strandedness=True,
        track_type=track_type,
    )

    canvas = assemble_circular_diagram_from_record(
        record,
        config_dict=config_dict,
        selected_features_set=SELECTED_FEATURES,
        legend="left",
        track_specs=["features@w=75px"],
    )
    svg_text = canvas.tostring()

    axis_transform = _extract_group_translate(svg_text, "Axis")
    record_transform = _extract_group_translate(svg_text, record.id)
    assert axis_transform is not None
    assert record_transform is not None
    assert axis_transform == pytest.approx(record_transform, abs=1e-6)

    for group_id in ["tick", "gc_content", "skew", "labels"]:
        group_transform = _extract_group_translate(svg_text, group_id)
        if group_transform is not None:
            assert group_transform == pytest.approx(axis_transform, abs=1e-6)
