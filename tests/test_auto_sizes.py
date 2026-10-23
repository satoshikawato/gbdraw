from __future__ import annotations

import json
from dataclasses import asdict
from functools import reduce
from operator import getitem
from pathlib import Path
from types import SimpleNamespace
from typing import Any, cast

import pytest

import gbdraw.api.diagram as diagram_api
from gbdraw.api.config import apply_config_overrides
from gbdraw.auto_sizes import (
    DrawnExtent,
    auto_setting_values,
    circular_tick_intervals,
    circular_track_width_px,
    depth_sliding_window,
    determine_length_parameter,
    linear_label_font_size,
    linear_label_size_lengths,
    linear_tick_interval,
    scalar_axis_tick_font_size,
    sliding_window,
    tick_and_window_tiers,
)
from gbdraw.canvas import CircularCanvasConfigurator, LinearCanvasConfigurator
from gbdraw.config.models import CircularRenderProfile, GbdrawConfig, LinearRenderProfile
from gbdraw.config.models.labels import LinearLabelScope
from gbdraw.diagrams.circular.presets import CircularPresetContext, circular_track_slots_for_preset
from gbdraw.diagrams.linear.precalc import _resolve_linear_diagram_label_font_size
from gbdraw.layout.circular_depth_axis import depth_axis_tick_font_size_px
from gbdraw.layout.scalar_axis import linear_scalar_axis_tick_font_size_px
from gbdraw.mode_profiles import DiagramMode, mode_profiles_payload
from gbdraw.svg.circular_ticks import get_circular_tick_intervals


VECTORS = json.loads(
    (Path(__file__).parent / "fixtures" / "auto_size_vectors.json").read_text(encoding="utf-8")
)


@pytest.fixture(scope="module")
def cfg() -> GbdrawConfig:
    return apply_config_overrides(None, None)


def _record(length: int) -> SimpleNamespace:
    return SimpleNamespace(seq="A" * length)


def _config_value(cfg: GbdrawConfig, path: str) -> object:
    return reduce(getattr, path.split("."), cfg)


@pytest.mark.parametrize(
    ("length", "expected"),
    [(1, "short"), (99_999, "short"), (100_000, "long"), (13_000_000, "long")],
)
def test_size_class_threshold(length: int, expected: str) -> None:
    assert determine_length_parameter(length, 100_000) == expected


@pytest.mark.parametrize(("mode", "boundary"), [("circular", 75_000), ("linear", 100_000)])
def test_packaged_size_class_boundary(cfg: GbdrawConfig, mode: DiagramMode, boundary: int) -> None:
    assert DrawnExtent(mode, (boundary - 1,)).size_class(cfg) == "short"
    assert DrawnExtent(mode, (boundary,)).size_class(cfg) == "long"


@pytest.mark.parametrize(
    ("length", "expected"),
    [
        (1, (1_000, 100)),
        (999_999, (1_000, 100)),
        (1_000_000, (10_000, 1_000)),
        (9_999_999, (10_000, 1_000)),
        (10_000_000, (100_000, 10_000)),
    ],
)
def test_sliding_window_tiers(cfg: GbdrawConfig, length: int, expected: tuple[int, int]) -> None:
    assert sliding_window(length, cfg) == expected


@pytest.mark.parametrize(
    ("window", "step", "expected"),
    [(500, 50, (100, 5)), (1_000, 9, (100, 1)), (1_010, 100, (101, 10)), (100_000, 10_000, (10_000, 1_000))],
)
def test_depth_window_follows_the_gc_window(window: int, step: int, expected: tuple[int, int]) -> None:
    assert depth_sliding_window(window, step) == expected


@pytest.mark.parametrize(
    ("length", "expected"),
    [
        (1, (1_000, 100)),
        (30_000, (1_000, 100)),
        (30_001, (5_000, 1_000)),
        (50_000, (5_000, 1_000)),
        (50_001, (10_000, 1_000)),
        (150_000, (10_000, 1_000)),
        (150_001, (50_000, 10_000)),
        (1_000_000, (50_000, 10_000)),
        (1_000_001, (500_000, 100_000)),
        (10_000_000, (500_000, 100_000)),
        (10_000_001, (1_000_000, 200_000)),
    ],
)
def test_circular_tick_tiers(length: int, expected: tuple[int, int]) -> None:
    assert circular_tick_intervals(length) == expected
    assert get_circular_tick_intervals(length) == expected


def test_circular_manual_tick_interval_bypasses_the_tiers() -> None:
    assert get_circular_tick_intervals(5_000_000, manual_interval=2_000) == (2_000, 200)


@pytest.mark.parametrize(
    ("length", "expected"),
    [
        (1, 100),
        (1_999, 100),
        (2_000, 1_000),
        (19_999, 1_000),
        (20_000, 5_000),
        (49_999, 5_000),
        (50_000, 10_000),
        (149_999, 10_000),
        (150_000, 50_000),
        (249_999, 50_000),
        (250_000, 100_000),
        (999_999, 100_000),
        (1_000_000, 200_000),
        (1_999_999, 200_000),
        (2_000_000, 500_000),
        (4_999_999, 500_000),
        (5_000_000, 1_000_000),
    ],
)
def test_linear_tick_tiers(length: int, expected: int) -> None:
    assert linear_tick_interval(length) == expected


@pytest.mark.parametrize(
    ("mode", "thickness", "expected"),
    [
        ("circular", 10.0, 5.0),
        ("circular", 30.0, 6.6),
        ("circular", 74.1, 8.0),
        ("linear", 5.0, 5.0),
        ("linear", 10.0, 7.0),
        ("linear", 20.0, 8.0),
    ],
)
def test_scalar_axis_tick_font_clamps_a_fraction_of_the_thickness(
    mode: DiagramMode, thickness: float, expected: float
) -> None:
    assert scalar_axis_tick_font_size(mode, thickness) == pytest.approx(expected)


def test_explicit_tick_fonts_bypass_the_auto_size() -> None:
    explicit = SimpleNamespace(tick_font_size=11)
    assert depth_axis_tick_font_size_px(explicit, 74.1) == 11.0
    assert linear_scalar_axis_tick_font_size_px(explicit, 20.0) == 11.0


@pytest.mark.parametrize(
    ("renderer", "factor_index", "scale"),
    [
        ("features", 0, 1.0),
        ("sequence_conservation", 0, 1.0),
        ("dinucleotide_content", 1, 1.0),
        ("depth", 1, 0.5),
        ("dinucleotide_skew", 2, 1.0),
    ],
)
def test_circular_track_width_per_renderer(renderer: str, factor_index: int, scale: float) -> None:
    factors = (0.5, 0.75, 0.25)
    assert circular_track_width_px(renderer, radius=400, track_ratio=0.2, factors=factors) == pytest.approx(
        80 * factors[factor_index] * scale
    )


@pytest.mark.parametrize(
    ("scope_lengths", "size"),
    [
        ((200_000, 20_000), "long"),
        ((20_000, 200_000), "long"),
        ((20_000, 99_999), "short"),
        ((200_000, 100_000), "long"),
    ],
)
def test_linear_label_font_follows_the_longest_record_in_the_label_scope(
    cfg: GbdrawConfig, scope_lengths: tuple[int, ...], size: str
) -> None:
    assert linear_label_font_size(scope_lengths, cfg) == getattr(cfg.labels.font_size.linear, size)


# The Web's Auto note reads which records set the label font from the generated
# mode profiles; the renderer draws with the same answer for every scope.
@pytest.mark.parametrize(
    ("scope", "sized_lengths", "size"),
    [
        ("first", (20_000,), "short"),
        ("all", (20_000, 200_000), "long"),
        ("orthogroup_top", (20_000, 200_000), "long"),
        # No record in scope: labels shown by per-feature overrides take the size of all records.
        ("none", (20_000, 200_000), "long"),
    ],
)
def test_linear_label_font_follows_the_records_of_the_label_scope(
    cfg: GbdrawConfig, scope: LinearLabelScope, sized_lengths: tuple[int, ...], size: str
) -> None:
    lengths = (20_000, 200_000)
    assert tuple(linear_label_size_lengths(lengths, scope)) == sized_lengths
    published = cast(Any, mode_profiles_payload()["autoSizes"])["linearLabelSizeRecords"][scope]
    assert published == ("first" if len(sized_lengths) == 1 else "all")

    profile = LinearRenderProfile(apply_config_overrides(cfg, {"labels.linear.scope": scope}))
    canvas = LinearCanvasConfigurator(
        num_of_entries=len(lengths), longest_genome=max(lengths), profile=profile, legend="none"
    )
    assert canvas.length_param == "long"
    assert _resolve_linear_diagram_label_font_size(
        [_record(length) for length in lengths], canvas_config=canvas, profile=profile
    ) == getattr(cfg.labels.font_size.linear, size)


def test_drawn_extent_size_class_follows_the_longest_record(cfg: GbdrawConfig) -> None:
    assert DrawnExtent("linear", (99_999, 3_000)).size_class(cfg) == "short"
    assert DrawnExtent("circular", (3_000, 75_000)).size_class(cfg) == "long"


def test_drawn_extent_rejects_an_empty_or_unknown_extent() -> None:
    with pytest.raises(ValueError):
        DrawnExtent("linear", ())
    with pytest.raises(ValueError):
        DrawnExtent("linear", (0,))
    with pytest.raises(ValueError):
        DrawnExtent(cast(Any, "radial"), (1_000,))


_LENGTHS = (3_000, 30_000, 200_000, 5_000_000)
_PRESET_SLOT_RENDERERS = {"features": "features", "gc_content": "gc_content", "gc_skew": "gc_skew", "depth": "depth"}


@pytest.mark.parametrize("length", _LENGTHS)
def test_circular_auto_values_equal_render_time_resolution(cfg: GbdrawConfig, length: int) -> None:
    values = auto_setting_values(DrawnExtent("circular", (length,)), cfg)
    profile = CircularRenderProfile(cfg)
    record = _record(length)
    canvas = CircularCanvasConfigurator(output_prefix="t", profile=profile, legend="none", gb_record=record)
    size = canvas.length_param
    for path in (
        "labels.font_size",
        "objects.axis.circular.stroke_width",
        "objects.features.block_stroke_width",
        "objects.features.line_stroke_width",
        "objects.legends.font_size",
        "objects.legends.color_rect_size",
    ):
        assert values[path] == getattr(_config_value(cfg, path), size), path
    slots = circular_track_slots_for_preset(
        "tuckin",
        CircularPresetContext(
            cfg=cfg,
            geometry=canvas.geometry(),
            total_length=length,
            strandedness=False,
            show_features=True,
            show_ticks=False,
            show_depth=True,
            show_gc=True,
            show_skew=True,
            track_ids=tuple(canvas.track_ids.items()),
        ),
    )
    widths = {slot.id: slot.width.value for slot in slots if slot.width is not None}
    for slot_id in _PRESET_SLOT_RENDERERS:
        assert values[f"circularSlot.{slot_id}.width"] == widths[slot_id], slot_id
    window_step = diagram_api._resolve_circular_window_step(record, cfg, window=None, step=None)
    assert (values["window"], values["step"]) == window_step
    assert (values["depth_window"], values["depth_step"]) == diagram_api._resolve_depth_window_step(
        window=window_step[0], step=window_step[1], depth_window=None, depth_step=None
    )
    assert values["objects.scale.interval"] == get_circular_tick_intervals(length)[0]
    auto_font = SimpleNamespace(tick_font_size=None)
    assert values["objects.gc_content.tick_font_size"] == depth_axis_tick_font_size_px(auto_font, widths["gc_content"])
    assert values["objects.depth.tick_font_size"] == depth_axis_tick_font_size_px(auto_font, widths["depth"])


@pytest.mark.parametrize("lengths", [(length,) for length in _LENGTHS] + [(200_000, 30_000)])
def test_linear_auto_values_equal_render_time_resolution(cfg: GbdrawConfig, lengths: tuple[int, ...]) -> None:
    values = auto_setting_values(DrawnExtent("linear", lengths), cfg)
    profile = LinearRenderProfile(apply_config_overrides(cfg, {"labels.linear.scope": "all"}))
    canvas = LinearCanvasConfigurator(
        num_of_entries=len(lengths), longest_genome=max(lengths), profile=profile, legend="none"
    )
    size = canvas.length_param
    assert values["canvas.linear.default_cds_height"] == canvas.default_cds_height
    for path in (
        "canvas.linear.arrow_length_parameter",
        "objects.axis.linear.stroke_width",
        "objects.definition.linear.font_size",
        "objects.scale.font_size",
        "objects.scale.ruler_label_font_size",
        "objects.features.block_stroke_width",
        "objects.features.line_stroke_width",
        "objects.legends.font_size",
        "objects.legends.color_rect_size",
    ):
        assert values[path] == getattr(_config_value(cfg, path), size), path
    assert values["labels.font_size.linear"] == _resolve_linear_diagram_label_font_size(
        [_record(length) for length in lengths], canvas_config=canvas, profile=profile
    )
    assert (values["window"], values["step"]) == sliding_window(max(lengths), cfg)
    assert values["objects.scale.interval"] == linear_tick_interval(max(lengths))
    assert values["canvas.linear.default_gc_height"] == canvas.default_gc_height
    assert values["canvas.linear.depth_height"] == canvas.default_depth_height
    auto_font = SimpleNamespace(tick_font_size=None)
    assert values["objects.gc_content.tick_font_size"] == linear_scalar_axis_tick_font_size_px(
        auto_font, canvas.default_gc_height
    )
    assert values["objects.depth.tick_font_size"] == linear_scalar_axis_tick_font_size_px(
        auto_font, canvas.default_depth_height
    )


def test_mixed_circular_canvas_reports_the_harmonized_long_styles(cfg: GbdrawConfig) -> None:
    lengths = (30_000, 200_000)
    values = auto_setting_values(DrawnExtent("circular", lengths), cfg)
    harmonized = diagram_api._harmonize_multi_record_circular_style_cfg(cfg, record_lengths=lengths)
    profile = CircularRenderProfile(harmonized)
    for length in lengths:
        canvas = CircularCanvasConfigurator(output_prefix="t", profile=profile, legend="none", gb_record=_record(length))
        for path in (
            "labels.font_size",
            "objects.axis.circular.stroke_width",
            "objects.features.block_stroke_width",
            "objects.features.line_stroke_width",
        ):
            assert values[path] == getattr(_config_value(harmonized, path), canvas.length_param), (length, path)
        slots = circular_track_slots_for_preset(
            "tuckin",
            CircularPresetContext(
                cfg=harmonized,
                geometry=canvas.geometry(),
                total_length=length,
                strandedness=False,
                show_features=True,
                show_ticks=False,
                show_depth=True,
                show_gc=True,
                show_skew=True,
                track_ids=tuple(canvas.track_ids.items()),
            ),
        )
        for slot in slots:
            assert slot.width is not None
            assert values[f"circularSlot.{slot.id}.width"] == slot.width.value, (length, slot.id)


def _short_long_paths(node: object, path: tuple[str, ...] = ()) -> list[tuple[str, ...]]:
    if not isinstance(node, dict):
        return []
    found = [path] if {"short", "long"} <= node.keys() else []
    return found + [hit for key, value in node.items() for hit in _short_long_paths(value, (*path, key))]


@pytest.mark.parametrize("lengths", [(306_008, 50_466), (30_000, 200_000)])
def test_mixed_circular_canvas_gives_every_record_the_long_class(cfg: GbdrawConfig, lengths: tuple[int, ...]) -> None:
    """OV-492: the canvas config resolves the long value of every short/long node for every record.

    This checks the config, not its readers; the render check of the readers is
    test_multi_record_canvas_draws_a_short_record_beside_a_long_one_as_long.
    """
    canvas_cfg = diagram_api._harmonize_multi_record_circular_style_cfg(cfg, record_lengths=lengths)
    packaged, resolved = asdict(cfg), asdict(canvas_cfg)
    paths = _short_long_paths(packaged)
    assert ("labels", "radius_factor", "tuckin", "single") in paths
    assert ("canvas", "circular", "track_dict") in paths
    for length in lengths:
        size = determine_length_parameter(length, int(canvas_cfg.labels.length_threshold.circular))
        assert size == DrawnExtent("circular", lengths).size_class(cfg)
        for path in paths:
            assert reduce(getitem, path, resolved)[size] == reduce(getitem, path, packaged)["long"], (length, path)


@pytest.mark.parametrize("lengths", [(30_000,), (30_000, 60_000), (200_000, 80_000), (200_000,)])
def test_circular_canvas_of_one_size_class_keeps_its_config(cfg: GbdrawConfig, lengths: tuple[int, ...]) -> None:
    assert diagram_api._harmonize_multi_record_circular_style_cfg(cfg, record_lengths=lengths) is cfg


def test_dependent_tick_fonts_follow_a_given_track_thickness(cfg: GbdrawConfig) -> None:
    linear = auto_setting_values(DrawnExtent("linear", (3_000,)), cfg, track_thickness_px={"depth": 5.0})
    assert linear["objects.depth.tick_font_size"] == 5.0
    assert linear["objects.gc_content.tick_font_size"] == 8.0
    circular = auto_setting_values(DrawnExtent("circular", (3_000,)), cfg, track_thickness_px={"gc_content": 30.0})
    assert circular["objects.gc_content.tick_font_size"] == pytest.approx(6.6)


def test_auto_values_list_only_the_settings_drawn_in_the_mode(cfg: GbdrawConfig) -> None:
    circular = auto_setting_values(DrawnExtent("circular", (3_000,)), cfg)
    linear = auto_setting_values(DrawnExtent("linear", (3_000,)), cfg)
    assert "labels.font_size.linear" not in circular
    assert not any(key.startswith("circularSlot.") for key in linear)
    assert "labels.font_size" not in linear
    assert "objects.axis.circular.stroke_width" not in linear
    assert "objects.axis.linear.stroke_width" not in circular


@pytest.mark.parametrize("case", VECTORS["autoHints"], ids=lambda case: f"{case['mode']}-{case['length']}")
def test_auto_hint_vectors_equal_python_tiers(cfg: GbdrawConfig, case: dict[str, Any]) -> None:
    length = case["length"]
    interval = circular_tick_intervals(length)[0] if case["mode"] == "circular" else linear_tick_interval(length)
    assert interval == case["scaleInterval"]
    assert sliding_window(length, cfg) == (case["window"], case["step"])


def test_payload_tiers_are_the_python_tiers(cfg: GbdrawConfig) -> None:
    sizes = mode_profiles_payload()["autoSizes"]
    assert isinstance(sizes, dict)
    for key, value in tick_and_window_tiers(cfg).items():
        assert sizes[key] == value
