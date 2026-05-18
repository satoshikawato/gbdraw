from __future__ import annotations

from types import SimpleNamespace

import pytest

from gbdraw.canvas import CircularCanvasConfigurator
from gbdraw.config.models import GbdrawConfig
from gbdraw.config.modify import modify_config_dict
from gbdraw.config.toml import load_config_toml
from gbdraw.configurators import DepthConfigurator
from gbdraw.layout.circular_depth_axis import resolve_depth_axis_footprint
from gbdraw.diagrams.circular.radial_layout import build_circular_feature_layout, resolve_circular_radial_layout
from gbdraw.layout.circular import calculate_feature_position_factors_circular
from gbdraw.tracks import CircularTrackSlot, ScalarSpec


class _Feature:
    def __init__(self, track_id: int) -> None:
        self.feature_track_id = track_id


def test_middle_combined_feature_center_is_axis_radius() -> None:
    layout = build_circular_feature_layout(
        {"a": _Feature(0)},
        axis_radius_px=390.0,
        width_px=74.1,
        track_type="middle",
        strandedness=False,
    )

    assert layout is not None
    assert layout.lanes_by_track_id[0].center_px == pytest.approx(390.0)


def test_spreadout_and_tuckin_feature_lanes_stay_on_expected_side() -> None:
    spreadout = build_circular_feature_layout(
        {"a": _Feature(0), "b": _Feature(1)},
        axis_radius_px=390.0,
        width_px=74.1,
        track_type="spreadout",
        strandedness=False,
    )
    tuckin = build_circular_feature_layout(
        {"a": _Feature(0), "b": _Feature(1)},
        axis_radius_px=390.0,
        width_px=74.1,
        track_type="tuckin",
        strandedness=False,
    )

    assert spreadout is not None
    assert tuckin is not None
    assert all(lane.inner_px > 390.0 for lane in spreadout.lanes_by_track_id.values())
    assert all(lane.outer_px < 390.0 for lane in tuckin.lanes_by_track_id.values())


def test_tuckin_combined_feature_lane_sits_close_to_axis() -> None:
    layout = build_circular_feature_layout(
        {"a": _Feature(0)},
        axis_radius_px=390.0,
        width_px=74.1,
        track_type="tuckin",
        strandedness=False,
    )

    assert layout is not None
    assert layout.lanes_by_track_id[0].center_px == pytest.approx(390.0 - (0.75 * 74.1))
    assert layout.lanes_by_track_id[0].outer_px == pytest.approx(390.0 - (0.25 * 74.1))


def test_tuckin_combined_feature_position_factors_sit_close_to_axis() -> None:
    factors = calculate_feature_position_factors_circular(
        total_length=1000,
        strand="positive",
        track_ratio=0.19,
        cds_ratio=0.1,
        offset=0.0,
        track_type="tuckin",
        strandedness=False,
        track_id=0,
    )

    assert factors == pytest.approx([0.88, 0.93, 0.98])


def test_tuckin_combined_overlap_position_factors_still_move_inward() -> None:
    factors = calculate_feature_position_factors_circular(
        total_length=1000,
        strand="positive",
        track_ratio=0.19,
        cds_ratio=0.1,
        offset=0.0,
        track_type="tuckin",
        strandedness=False,
        track_id=2,
    )

    assert factors == pytest.approx([0.64, 0.69, 0.74])


def test_custom_core_slot_order_places_ticks_outside_features() -> None:
    canvas_config, cfg = _small_radial_canvas()
    layout = resolve_circular_radial_layout(
        total_length=5_500_000,
        canvas_config=canvas_config,
        cfg=cfg,
        slots=[
            CircularTrackSlot(id="ticks", renderer="ticks"),
            CircularTrackSlot(id="features", renderer="features"),
        ],
        feature_dict={"a": _Feature(0)},
        show_features=True,
        show_ticks=True,
    )

    assert layout.features is not None
    assert layout.ticks is not None
    assert layout.ticks.reserved_band_px.inner_px > layout.features.all_band_px.outer_px


def test_custom_core_slot_order_keeps_default_feature_then_ticks_order() -> None:
    canvas_config, cfg = _small_radial_canvas()
    layout = resolve_circular_radial_layout(
        total_length=5_500_000,
        canvas_config=canvas_config,
        cfg=cfg,
        slots=[
            CircularTrackSlot(id="features", renderer="features"),
            CircularTrackSlot(id="ticks", renderer="ticks"),
        ],
        feature_dict={"a": _Feature(0)},
        show_features=True,
        show_ticks=True,
    )

    assert layout.features is not None
    assert layout.ticks is not None
    assert layout.features.all_band_px.inner_px > layout.ticks.reserved_band_px.outer_px


def test_inside_auto_features_fit_when_anchor_sits_outside_free_interval() -> None:
    config_dict = modify_config_dict(
        load_config_toml("gbdraw.data", "config.toml"),
        show_labels=False,
        track_type="tuckin",
        strandedness=True,
    )
    cfg = GbdrawConfig.from_dict(config_dict)
    record = SimpleNamespace(seq="N" * 1000)
    canvas_config = CircularCanvasConfigurator("test", config_dict, "none", record, cfg=cfg)
    canvas_config.radius = 100.0

    layout = resolve_circular_radial_layout(
        total_length=1000,
        canvas_config=canvas_config,
        cfg=cfg,
        slots=[
            CircularTrackSlot(
                id="upper_spacer",
                renderer="spacer",
                radius=ScalarSpec(86.0, "px"),
                width=ScalarSpec(8.0, "px"),
            ),
            CircularTrackSlot(
                id="features",
                renderer="features",
                width=ScalarSpec(10.0, "px"),
            ),
        ],
        feature_dict={"plus": _Feature(0), "minus": _Feature(-1)},
        show_features=True,
        show_ticks=False,
        definition_reserved_radius_px=60.0,
    )
    by_id = {slot.id: slot for slot in layout.slots}

    assert by_id["features"].reserved_band_px.inner_px == pytest.approx(60.0)
    assert by_id["features"].reserved_band_px.outer_px <= by_id["upper_spacer"].reserved_band_px.inner_px
    assert by_id["features"].anchor_radius_px > by_id["features"].reserved_band_px.outer_px


def _small_radial_canvas() -> tuple[CircularCanvasConfigurator, GbdrawConfig]:
    config_dict = modify_config_dict(
        load_config_toml("gbdraw.data", "config.toml"),
        show_labels=False,
        show_gc=True,
        show_skew=True,
        track_type="tuckin",
    )
    cfg = GbdrawConfig.from_dict(config_dict)
    record = SimpleNamespace(seq="N" * 1000)
    canvas_config = CircularCanvasConfigurator("test", config_dict, "none", record, cfg=cfg)
    canvas_config.radius = 100.0
    return canvas_config, cfg


def test_radial_numeric_stack_places_default_inside_tracks_when_space_allows() -> None:
    canvas_config, cfg = _small_radial_canvas()
    layout = resolve_circular_radial_layout(
        total_length=1000,
        canvas_config=canvas_config,
        cfg=cfg,
        slots=[
            CircularTrackSlot(id="gc_content", renderer="dinucleotide_content"),
            CircularTrackSlot(id="at_skew", renderer="dinucleotide_skew"),
        ],
        show_features=False,
        show_ticks=False,
        definition_reserved_radius_px=60.0,
    )
    by_id = {track.id: track for track in layout.tracks}

    assert by_id["gc_content"].draw_band_px.width_px == pytest.approx(19.0)
    assert by_id["at_skew"].draw_band_px.width_px == pytest.approx(19.0)
    assert by_id["gc_content"].draw_band_px.inner_px > by_id["at_skew"].draw_band_px.outer_px
    assert min(track.draw_band_px.inner_px for track in layout.tracks) >= 60.0 - 1e-6


def test_radial_numeric_stack_compresses_group_no_lower_than_readable_minimum() -> None:
    canvas_config, cfg = _small_radial_canvas()
    layout = resolve_circular_radial_layout(
        total_length=1000,
        canvas_config=canvas_config,
        cfg=cfg,
        slots=[
            CircularTrackSlot(id="gc_content", renderer="dinucleotide_content", compress=True),
            CircularTrackSlot(id="at_skew", renderer="dinucleotide_skew", compress=True),
        ],
        show_features=False,
        show_ticks=False,
        definition_reserved_radius_px=65.0,
    )
    by_id = {track.id: track for track in layout.tracks}

    assert by_id["gc_content"].draw_band_px.width_px < 19.0
    assert by_id["at_skew"].draw_band_px.width_px < 19.0
    assert by_id["gc_content"].draw_band_px.width_px >= 12.35 - 1e-6
    assert by_id["at_skew"].draw_band_px.width_px >= 12.35 - 1e-6
    assert min(track.draw_band_px.inner_px for track in layout.tracks) >= 65.0 - 1e-6


def test_radial_numeric_stack_raises_when_readable_minimum_cannot_fit() -> None:
    canvas_config, cfg = _small_radial_canvas()

    with pytest.raises(Exception, match="Circular track slot 'at_skew' cannot fit inside"):
        resolve_circular_radial_layout(
            total_length=1000,
            canvas_config=canvas_config,
            cfg=cfg,
            slots=[
                CircularTrackSlot(id="gc_content", renderer="dinucleotide_content", compress=True),
                CircularTrackSlot(id="gc_skew", renderer="dinucleotide_skew", compress=True, strict=True),
                CircularTrackSlot(id="at_skew", renderer="dinucleotide_skew", compress=True, strict=True),
            ],
            show_features=False,
            show_ticks=False,
            definition_reserved_radius_px=70.0,
        )


def test_outside_auto_numeric_width_is_independent_from_inside_compression() -> None:
    canvas_config, cfg = _small_radial_canvas()
    layout = resolve_circular_radial_layout(
        total_length=1000,
        canvas_config=canvas_config,
        cfg=cfg,
        slots=[
            CircularTrackSlot(id="outer_skew", renderer="dinucleotide_skew", side="outside"),
            CircularTrackSlot(id="gc_content", renderer="dinucleotide_content", compress=True),
            CircularTrackSlot(id="gc_skew", renderer="dinucleotide_skew", compress=True),
        ],
        show_features=False,
        show_ticks=False,
        definition_reserved_radius_px=65.0,
    )
    by_id = {track.id: track for track in layout.tracks}

    assert by_id["outer_skew"].side == "outside"
    assert by_id["outer_skew"].draw_band_px.width_px == pytest.approx(19.0)
    assert by_id["gc_skew"].draw_band_px.width_px < by_id["outer_skew"].draw_band_px.width_px
    assert by_id["outer_skew"].draw_band_px.inner_px >= canvas_config.radius


def test_user_preset_generated_param_has_no_layout_effect_on_pinned_numeric_slot() -> None:
    canvas_config, cfg = _small_radial_canvas()
    base_slot = CircularTrackSlot(
        id="gc_content",
        renderer="dinucleotide_content",
        radius=ScalarSpec(0.5, "factor"),
        width=ScalarSpec(20.0, "px"),
    )
    tagged_slot = CircularTrackSlot(
        id="gc_content",
        renderer="dinucleotide_content",
        radius=ScalarSpec(0.5, "factor"),
        width=ScalarSpec(20.0, "px"),
        params={"_preset_generated": True},
    )

    base = resolve_circular_radial_layout(
        total_length=1000,
        canvas_config=canvas_config,
        cfg=cfg,
        slots=[base_slot],
        show_features=False,
        show_ticks=False,
        definition_reserved_radius_px=55.0,
    ).tracks[0]
    tagged = resolve_circular_radial_layout(
        total_length=1000,
        canvas_config=canvas_config,
        cfg=cfg,
        slots=[tagged_slot],
        show_features=False,
        show_ticks=False,
        definition_reserved_radius_px=55.0,
    ).tracks[0]

    assert tagged.center_radius_px == pytest.approx(base.center_radius_px)
    assert tagged.draw_width_px == pytest.approx(base.draw_width_px)


def test_depth_reserved_band_includes_axis_radial_footprint_without_depth_df() -> None:
    config_dict = modify_config_dict(
        load_config_toml("gbdraw.data", "config.toml"),
        show_labels=False,
        show_gc=False,
        show_skew=False,
        show_depth=True,
    )
    cfg = GbdrawConfig.from_dict(config_dict)
    record = SimpleNamespace(seq="N" * 1000)
    canvas_config = CircularCanvasConfigurator("test", config_dict, "none", record, cfg=cfg)
    canvas_config.radius = 100.0
    depth_config = DepthConfigurator(10, 10, config_dict, cfg=cfg, show_axis=True, show_ticks=True)
    width_px = 24.0

    layout = resolve_circular_radial_layout(
        total_length=1000,
        canvas_config=canvas_config,
        cfg=cfg,
        slots=[
            CircularTrackSlot(
                id="depth",
                renderer="depth",
                radius=ScalarSpec(0.7, "factor"),
                width=ScalarSpec(width_px, "px"),
            )
        ],
        show_features=False,
        show_ticks=False,
        depth_config=depth_config,
    )
    depth = layout.tracks[0]
    footprint = resolve_depth_axis_footprint(depth_config, width_px)

    assert depth.draw_band_px is not None
    assert depth.reserved_band_px is not None
    assert depth.draw_band_px.width_px == pytest.approx(width_px)
    assert depth.reserved_band_px.inner_px == pytest.approx(
        depth.draw_band_px.inner_px - footprint.radial_inner_extra_px
    )
    assert depth.reserved_band_px.outer_px == pytest.approx(
        depth.draw_band_px.outer_px + footprint.radial_outer_extra_px
    )

    hidden_axis_config = DepthConfigurator(10, 10, config_dict, cfg=cfg, show_axis=False)
    hidden_axis_layout = resolve_circular_radial_layout(
        total_length=1000,
        canvas_config=canvas_config,
        cfg=cfg,
        slots=[
            CircularTrackSlot(
                id="depth",
                renderer="depth",
                radius=ScalarSpec(0.7, "factor"),
                width=ScalarSpec(width_px, "px"),
            )
        ],
        show_features=False,
        show_ticks=False,
        depth_config=hidden_axis_config,
    )
    hidden_axis_depth = hidden_axis_layout.tracks[0]
    assert hidden_axis_depth.reserved_band_px == hidden_axis_depth.draw_band_px
