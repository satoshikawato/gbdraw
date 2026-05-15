from __future__ import annotations

from types import SimpleNamespace

import pytest

from gbdraw.canvas import CircularCanvasConfigurator
from gbdraw.config.models import GbdrawConfig
from gbdraw.config.modify import modify_config_dict
from gbdraw.config.toml import load_config_toml
from gbdraw.diagrams.circular.radial_layout import build_circular_feature_layout, resolve_circular_radial_layout
from gbdraw.layout.circular import calculate_feature_position_factors_circular
from gbdraw.tracks import CircularTrackSlot
from gbdraw.tracks.spec import ScalarSpec


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
        honor_core_slot_order=True,
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
        honor_core_slot_order=True,
    )

    assert layout.features is not None
    assert layout.ticks is not None
    assert layout.features.all_band_px.inner_px > layout.ticks.reserved_band_px.outer_px


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


def test_radial_numeric_stack_compresses_no_lower_than_readable_minimum() -> None:
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
        definition_reserved_radius_px=70.0,
    )
    by_id = {track.id: track for track in layout.tracks}

    assert by_id["gc_content"].draw_band_px.width_px == pytest.approx(by_id["at_skew"].draw_band_px.width_px)
    assert by_id["gc_content"].draw_band_px.width_px < 19.0
    assert by_id["at_skew"].draw_band_px.width_px >= 12.35 - 1e-6
    assert min(track.draw_band_px.inner_px for track in layout.tracks) >= 70.0 - 1e-6


def test_outside_auto_numeric_width_tracks_inside_auto_compression() -> None:
    canvas_config, cfg = _small_radial_canvas()
    layout = resolve_circular_radial_layout(
        total_length=1000,
        canvas_config=canvas_config,
        cfg=cfg,
        slots=[
            CircularTrackSlot(
                id="outer_skew",
                renderer="dinucleotide_skew",
                params={"side": "outside"},
            ),
            CircularTrackSlot(id="gc_content", renderer="dinucleotide_content"),
            CircularTrackSlot(id="gc_skew", renderer="dinucleotide_skew"),
        ],
        show_features=False,
        show_ticks=False,
        definition_reserved_radius_px=70.0,
    )
    by_id = {track.id: track for track in layout.tracks}

    assert by_id["gc_skew"].draw_band_px.width_px < 19.0
    assert by_id["outer_skew"].draw_band_px.width_px == pytest.approx(
        by_id["gc_skew"].draw_band_px.width_px
    )


def test_explicit_outside_numeric_width_ignores_inside_auto_compression() -> None:
    canvas_config, cfg = _small_radial_canvas()
    layout = resolve_circular_radial_layout(
        total_length=1000,
        canvas_config=canvas_config,
        cfg=cfg,
        slots=[
            CircularTrackSlot(
                id="outer_skew",
                renderer="dinucleotide_skew",
                width=ScalarSpec(19.0, "px"),
                params={"side": "outside"},
            ),
            CircularTrackSlot(id="gc_content", renderer="dinucleotide_content"),
            CircularTrackSlot(id="gc_skew", renderer="dinucleotide_skew"),
        ],
        show_features=False,
        show_ticks=False,
        definition_reserved_radius_px=70.0,
    )
    by_id = {track.id: track for track in layout.tracks}

    assert by_id["gc_skew"].draw_band_px.width_px < 19.0
    assert by_id["outer_skew"].draw_band_px.width_px == pytest.approx(19.0)


def test_explicit_inside_numeric_group_compresses_below_readable_minimum() -> None:
    canvas_config, cfg = _small_radial_canvas()
    layout = resolve_circular_radial_layout(
        total_length=1000,
        canvas_config=canvas_config,
        cfg=cfg,
        slots=[
            CircularTrackSlot(
                id="gc_content",
                renderer="dinucleotide_content",
                params={"side": "inside", "compress": True},
            ),
            CircularTrackSlot(
                id="gc_skew",
                renderer="dinucleotide_skew",
                params={"side": "inside", "compress": True},
            ),
            CircularTrackSlot(
                id="at_skew",
                renderer="dinucleotide_skew",
                params={"side": "inside", "compress": True},
            ),
        ],
        show_features=False,
        show_ticks=False,
        definition_reserved_radius_px=70.0,
    )
    by_id = {track.id: track for track in layout.tracks}

    assert set(by_id) == {"gc_content", "gc_skew", "at_skew"}
    assert all(track.channel == "inside" for track in layout.tracks)
    assert max(track.draw_band_px.width_px for track in layout.tracks) < 12.35 - 1e-6
    assert min(track.draw_band_px.inner_px for track in layout.tracks) >= 70.0 - 1e-6
    assert max(track.draw_band_px.outer_px for track in layout.tracks) <= canvas_config.radius + 1e-6
    assert layout.outer_content_radius_px <= canvas_config.radius + 1e-6


def test_explicit_inside_numeric_group_does_not_fallback_when_not_strict() -> None:
    canvas_config, cfg = _small_radial_canvas()
    layout = resolve_circular_radial_layout(
        total_length=1000,
        canvas_config=canvas_config,
        cfg=cfg,
        slots=[
            CircularTrackSlot(
                id="gc_content",
                renderer="dinucleotide_content",
                params={"side": "inside", "strict": False, "compress": False},
            ),
            CircularTrackSlot(
                id="gc_skew",
                renderer="dinucleotide_skew",
                params={"side": "inside", "strict": False, "compress": False},
            ),
            CircularTrackSlot(
                id="at_skew",
                renderer="dinucleotide_skew",
                params={"side": "inside", "strict": False, "compress": False},
            ),
        ],
        show_features=False,
        show_ticks=False,
        definition_reserved_radius_px=70.0,
    )

    assert max(track.draw_band_px.outer_px for track in layout.tracks) <= canvas_config.radius + 1e-6
    assert min(track.draw_band_px.inner_px for track in layout.tracks) >= 70.0 - 1e-6


def test_explicit_inside_numeric_group_raises_below_hard_minimum() -> None:
    canvas_config, cfg = _small_radial_canvas()

    with pytest.raises(
        ValueError,
        match="Circular track slot 'at_skew' cannot fit inside the feature/tick stack.*Placement to Outside",
    ):
        resolve_circular_radial_layout(
            total_length=1000,
            canvas_config=canvas_config,
            cfg=cfg,
            slots=[
                CircularTrackSlot(
                    id="gc_content",
                    renderer="dinucleotide_content",
                    params={"side": "inside", "compress": True},
                ),
                CircularTrackSlot(
                    id="gc_skew",
                    renderer="dinucleotide_skew",
                    params={"side": "inside", "compress": True},
                ),
                CircularTrackSlot(
                    id="at_skew",
                    renderer="dinucleotide_skew",
                    params={"side": "inside", "compress": True},
                ),
            ],
            show_features=False,
            show_ticks=False,
            definition_reserved_radius_px=89.0,
        )


def test_implicit_numeric_group_keeps_legacy_fallback() -> None:
    canvas_config, cfg = _small_radial_canvas()
    layout = resolve_circular_radial_layout(
        total_length=1000,
        canvas_config=canvas_config,
        cfg=cfg,
        slots=[
            CircularTrackSlot(id="gc_content", renderer="dinucleotide_content"),
            CircularTrackSlot(id="gc_skew", renderer="dinucleotide_skew"),
            CircularTrackSlot(id="at_skew", renderer="dinucleotide_skew"),
        ],
        show_features=False,
        show_ticks=False,
        definition_reserved_radius_px=70.0,
    )

    assert max(track.draw_band_px.outer_px for track in layout.tracks) > canvas_config.radius


def test_strict_radial_inside_numeric_stack_reports_layout_error() -> None:
    canvas_config, cfg = _small_radial_canvas()

    with pytest.raises(ValueError, match="Circular track slot 'at_skew' cannot fit inside the feature/tick stack"):
        resolve_circular_radial_layout(
            total_length=1000,
            canvas_config=canvas_config,
            cfg=cfg,
            slots=[
                CircularTrackSlot(
                    id="gc_content",
                    renderer="dinucleotide_content",
                    params={"side": "inside", "strict": True, "compress": True},
                ),
                CircularTrackSlot(
                    id="at_skew",
                    renderer="dinucleotide_skew",
                    params={"side": "inside", "strict": True, "compress": True},
                ),
            ],
            show_features=False,
            show_ticks=False,
            definition_reserved_radius_px=93.0,
        )
