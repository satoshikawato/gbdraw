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
    ADAPTABLE_SETTINGS,
    DrawnDepthRow,
    DrawnExtent,
    SIZE_CLASSES,
    adapt_explicit_settings,
    auto_setting_values,
    circular_tick_intervals,
    circular_track_width_px,
    depth_axis_tick_font_size_px,
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
from gbdraw.layout.scalar_axis import linear_scalar_axis_tick_font_size_px
from gbdraw.layout.circular_depth_axis import CircularScalarAxes, DepthSlotWidth
from gbdraw.analysis.depth_tracks import DepthTrackStyle, depth_track_config
from gbdraw.configurators import DepthConfigurator
from gbdraw.mode_profiles import DiagramMode
from gbdraw.svg.circular_ticks import get_circular_tick_intervals
from gbdraw.web_support.mode_profile_data import mode_profiles_payload


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
        # Below 10 kb the window is a tenth of the length and the step a tenth of the window (OV-2864, Owner B).
        (1, (1, 1)),
        (671, (67, 6)),
        (2_501, (250, 25)),
        (9_000, (900, 90)),
        (9_999, (999, 99)),
        # At 10 kb the tier's own values: no step in the window across the boundary.
        (10_000, (1_000, 100)),
        (10_001, (1_000, 100)),
        (999_999, (1_000, 100)),
        (1_000_000, (10_000, 1_000)),
        (9_999_999, (10_000, 1_000)),
        (10_000_000, (100_000, 10_000)),
        (99_999_999, (100_000, 10_000)),
        (100_000_000, (1_000_000, 100_000)),
        (248_956_422, (1_000_000, 100_000)),
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
        # OV-2747, OV-2748: above 10 Mbp the intervals coarsen where the default labels come closer than a quarter em.
        (34_236_412, (1_000_000, 200_000)),
        (34_236_413, (2_000_000, 400_000)),
        (112_934_049, (5_000_000, 1_000_000)),
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
        # OV-2747, OV-2748: above 5 Mbp the intervals coarsen where the default ruler labels come closer than a quarter em.
        (46_484_179, 1_000_000),
        (46_484_180, 2_000_000),
        (112_934_049, 5_000_000),
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


def _vector_depth_rows(case: dict[str, Any], cfg: GbdrawConfig) -> tuple[DrawnDepthRow, ...]:
    """The depth rows of a gcAxisTickFonts vector: each series' axis as the drawing configures it, on the custom
    stack's Depth rows (``circularDepthSlots``), else on the one Depth row the Depth width sizes."""
    if not case["showDepth"]:
        return ()
    base = DepthConfigurator(10, 10, CircularRenderProfile(cfg))
    axes = [
        depth_track_config(DepthTrackStyle(None, None, None, track["tickFontSize"]), base, 10, 10, None)
        for track in case.get("depthTracks", [{"tickFontSize": None}])
    ]
    if "circularDepthSlots" not in case:
        return (DrawnDepthRow(axes[0], 0, depth_width=True),)
    return tuple(
        DrawnDepthRow(axes[slot["trackIndex"]], slot["trackIndex"], float(slot["widthPx"]))
        for slot in case["circularDepthSlots"]
    )


@pytest.mark.parametrize("case", VECTORS["gcAxisTickFonts"], ids=lambda case: case["name"])
def test_a_percent_gc_axis_takes_the_tick_font_of_the_depth_axis_beside_it(case: dict[str, Any]) -> None:
    # The drawing draws the GC axis in the font of its first depth row (OV-2873, OV-2875); its Auto value is that
    # font, read from the drawing's depth rows: series fonts and stack rows included (OV-2877).
    cfg = apply_config_overrides(None, case["config"])
    values = auto_setting_values(
        DrawnExtent(case["mode"], tuple(case["recordLengths"]), depth_rows=_vector_depth_rows(case, cfg)),
        cfg,
        track_thickness_px=case["trackThicknessPx"],
    )
    for setting, expected in case["expected"].items():
        assert values[setting] == pytest.approx(expected, abs=1e-9), setting


def test_the_gc_axis_takes_the_drawings_own_depth_font_when_its_first_row_has_none(cfg: GbdrawConfig) -> None:
    # A hand-set Depth tick font is the font a series without its own draws in (clone_depth_config), so the Percent
    # GC axis's Auto value follows it; a series' own font wins.
    percent = apply_config_overrides(None, {"objects.gc_content.mode": "percent"})
    base = DepthConfigurator(10, 10, CircularRenderProfile(percent))
    own = depth_track_config(DepthTrackStyle(None, None, None, 12.0), base, 10, 10, None)
    gc_font = "objects.gc_content.tick_font_size"
    for axis, expected in ((base, 7.0), (own, 12.0)):
        extent = DrawnExtent("circular", (16_569,), depth_rows=(DrawnDepthRow(axis, 0, depth_width=True),))
        assert auto_setting_values(extent, percent, depth_tick_font_size=7.0)[gc_font] == expected


def test_a_new_drawing_reads_its_depth_row_at_the_depth_width_it_keeps(cfg: GbdrawConfig) -> None:
    # The Depth row's width is the setting circularSlot.depth.width. A Circular crop keeps a hand-set 20 px (its
    # Auto value does not change), so its first Depth row still draws 5 px and the GC axis's hand-set font is kept;
    # a Linear drawing has no Depth row beside its GC axis.
    percent = apply_config_overrides(None, {"objects.gc_content.mode": "percent"})
    rows = (DrawnDepthRow(DepthConfigurator(10, 10, CircularRenderProfile(percent)), 0, depth_width=True),)
    explicit = {"circularSlot.depth.width": 20.0, "objects.gc_content.tick_font_size": 7}
    by_class = {setting: dict.fromkeys(SIZE_CLASSES, value) for setting, value in explicit.items()}
    source = DrawnExtent("circular", (2_000_000,), depth_rows=rows)
    crop = adapt_explicit_settings(
        by_class, source=source, target=DrawnExtent("circular", (3_000,), depth_rows=rows), cfg=percent
    )
    assert crop.kept == explicit and crop.reset == ()
    linear = adapt_explicit_settings(by_class, source=source, target=DrawnExtent("linear", (3_000,)), cfg=percent)
    assert {item.setting: (item.source_auto, item.target_auto) for item in linear.reset} == {
        "circularSlot.depth.width": (auto_setting_values(source, percent)["circularSlot.depth.width"], None),
        "objects.gc_content.tick_font_size": (5.0, 8.0),
    }


@pytest.mark.parametrize("case", VECTORS["gcAxisTickFonts"], ids=lambda case: case["name"])
def test_the_gc_axis_tick_font_vectors_are_the_fonts_the_circular_drawing_draws(case: dict[str, Any]) -> None:
    # OV-2875/OV-2876: the radial layout and the drawing read each depth slot's own track (its font) at the slot's
    # requested width; the Percent GC axis takes the font of the first depth slot (CircularScalarAxes).
    cfg = apply_config_overrides(None, case["config"])
    thickness = case["trackThicknessPx"]
    extent = DrawnExtent(case["mode"], tuple(case["recordLengths"]))
    auto = auto_setting_values(extent, cfg)
    depth_width = float(thickness.get("depth", auto["circularSlot.depth.width"]))
    slots = [
        DepthSlotWidth(int(slot["trackIndex"]), float(slot["widthPx"]))
        for slot in case.get("circularDepthSlots", [{"trackIndex": 0, "widthPx": depth_width}])
    ]
    profile = CircularRenderProfile(cfg)
    base = DepthConfigurator(10, 10, profile)
    tracks = {
        index: depth_track_config(DepthTrackStyle(None, None, None, track["tickFontSize"]), base, 10, 10, None)
        for index, track in enumerate(case.get("depthTracks", [{"tickFontSize": None}]))
    }
    gc = cfg.objects.gc_content
    axes = CircularScalarAxes.of(base if case["showDepth"] else None, gc, tracks)
    gc_font = axes.gc_tick_font_size(gc, slots) if gc.mode == "percent" else None
    if gc_font is None:
        gc_font = scalar_axis_tick_font_size(
            "circular", float(thickness.get("gc_content", auto["circularSlot.gc_content.width"]))
        )
    expected = case["expected"]
    assert gc_font == pytest.approx(expected["objects.gc_content.tick_font_size"], abs=1e-9)
    # The Auto value of the drawing's own Depth tick font setting: the size of the width of series 0's first row.
    series_0_width = next(width for index, width in slots if index == 0)
    width_font = scalar_axis_tick_font_size("circular", series_0_width)
    assert width_font == pytest.approx(expected["objects.depth.tick_font_size"], abs=1e-9)
    # The Auto value of series 0's own Tick Font Size field (OV-3496): the font a series without its own draws in,
    # the drawing's Depth tick font when it has one, else that size.
    series_auto = base.tick_font_size if base.tick_font_size is not None else width_font
    assert series_auto == pytest.approx(case["seriesAutoTickFontSize"], abs=1e-9)


def test_auto_values_list_only_the_settings_drawn_in_the_mode(cfg: GbdrawConfig) -> None:
    circular = auto_setting_values(DrawnExtent("circular", (3_000,)), cfg)
    linear = auto_setting_values(DrawnExtent("linear", (3_000,)), cfg)
    assert "labels.font_size.linear" not in circular
    assert not any(key.startswith("circularSlot.") for key in linear)
    assert "labels.font_size" not in linear
    assert "objects.axis.circular.stroke_width" not in linear
    assert "objects.axis.linear.stroke_width" not in circular


@pytest.mark.parametrize("mode", ["circular", "linear"])
@pytest.mark.parametrize("length", [3_000, 2_000_000])
def test_adaptable_settings_name_every_auto_value(cfg: GbdrawConfig, mode: DiagramMode, length: int) -> None:
    extent = DrawnExtent(mode, (length,))
    assert extent.size_class(cfg) == ("short" if length == 3_000 else "long")
    assert list(auto_setting_values(extent, cfg)) == list(ADAPTABLE_SETTINGS[mode])


@pytest.mark.parametrize("case", VECTORS["autoHints"], ids=lambda case: f"{case['mode']}-{case['length']}")
def test_auto_hint_vectors_equal_python_tiers(cfg: GbdrawConfig, case: dict[str, Any]) -> None:
    length = case["length"]
    interval = circular_tick_intervals(length)[0] if case["mode"] == "circular" else linear_tick_interval(length)
    assert interval == case["scaleInterval"]
    assert sliding_window(length, cfg) == (case["window"], case["step"])
    assert depth_sliding_window(case["window"], case["step"]) == (case["depthWindow"], case["depthStep"])


def test_payload_tiers_are_the_python_tiers(cfg: GbdrawConfig) -> None:
    sizes = mode_profiles_payload()["autoSizes"]
    assert isinstance(sizes, dict)
    for key, value in tick_and_window_tiers(cfg).items():
        assert sizes[key] == value
