from __future__ import annotations

import logging
import re
from pathlib import Path

import pandas as pd
import pytest
from Bio import SeqIO
from svgwrite import Drawing

import gbdraw.circular as circular_cli_module
from gbdraw.api.diagram import assemble_circular_diagram_from_record
from gbdraw.config.modify import modify_config_dict
from gbdraw.config.toml import load_config_toml
from gbdraw.io.colors import load_default_colors
from gbdraw.svg.circular_ticks import (
    get_circular_tick_label_radius_bounds,
    get_circular_tick_path_radius_bounds,
    get_circular_tick_path_ratio_bounds,
    generate_circular_tick_labels,
    resolve_circular_tick_label_geometry,
    set_tick_label_anchor_value,
)
from gbdraw.tracks import (
    CircularTrackPlacement,
    CircularTrackLayoutContext,
    CircularTrackSlot,
    ScalarSpec,
    TrackSpec,
    TrackSpecParseError,
    default_circular_track_slots,
    parse_circular_track_slot,
    parse_circular_track_slots,
    resolve_circular_track_slots,
)


SELECTED_FEATURES = ["CDS", "rRNA", "tRNA", "tmRNA", "ncRNA", "misc_RNA", "repeat_region"]


def _load_record():
    input_path = Path(__file__).parent / "test_inputs" / "HmmtDNA.gbk"
    return SeqIO.read(str(input_path), "genbank")


def _load_edl933_record():
    input_path = Path(__file__).parent / "test_inputs" / "EDL933.gbk"
    return SeqIO.read(str(input_path), "genbank")


def _base_config(*, track_type: str = "middle"):
    return modify_config_dict(
        load_config_toml("gbdraw.data", "config.toml"),
        show_labels=False,
        show_gc=True,
        show_skew=True,
        track_type=track_type,
    )


def _depth_table(record_id: str) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "reference_name": [record_id, record_id, record_id, record_id],
            "position": [1, 2, 50, 100],
            "depth": [10, 20, 40, 80],
        }
    )


def _normalize_svg_auto_ids(svg_text: str) -> str:
    return re.sub(r"id\d+", "id_auto", svg_text)


def _axis_circle_radius(svg_text: str) -> float:
    match = re.search(r'<g id="Axis"[^>]*>.*?<circle\b[^>]*\br="([^"]+)"', svg_text, re.S)
    assert match is not None
    return float(match.group(1))


def _capture_circular_core_geometry(
    monkeypatch: pytest.MonkeyPatch,
    *,
    track_type: str,
    strandedness: bool = True,
    show_depth: bool = False,
    use_default_slots: bool = False,
) -> dict[str, tuple[float, float]]:
    import gbdraw.diagrams.circular.assemble as circular_assemble_module

    record = _load_record()
    config_dict = modify_config_dict(
        load_config_toml("gbdraw.data", "config.toml"),
        show_labels=False,
        show_gc=True,
        show_skew=True,
        track_type=track_type,
        strandedness=strandedness,
    )
    default_colors = load_default_colors("", palette="default")
    captured: dict[str, tuple[float, float]] = {}

    def fake_add_record_group_on_canvas(
        canvas,
        gb_record,
        canvas_config,
        feature_config,
        config_dict,
        *,
        cfg=None,
        precomputed_feature_dict=None,
        precalculated_labels=None,
        feature_track_ratio_factor_override=None,
        feature_anchor_radius_px=None,
    ):
        assert cfg is not None
        ratio_factor = (
            float(feature_track_ratio_factor_override)
            if feature_track_ratio_factor_override is not None
            else float(cfg.canvas.circular.track_ratio_factors[str(canvas_config.length_param)][0])
        )
        captured["features"] = (
            float(feature_anchor_radius_px if feature_anchor_radius_px is not None else canvas_config.radius),
            float(canvas_config.radius) * float(canvas_config.track_ratio) * ratio_factor,
        )
        return canvas

    def fake_add_tick_group_on_canvas(
        canvas,
        gb_record,
        canvas_config,
        config_dict,
        *,
        radius_override=None,
        tick_track_channel_override=None,
        label_side="legacy",
        tick_side="legacy",
        tick_length_px=None,
        cfg=None,
    ):
        center = float(radius_override if radius_override is not None else canvas_config.radius)
        if tick_length_px is not None and str(tick_side).strip().lower() != "legacy":
            width_px = float(tick_length_px)
        else:
            tick_min_ratio, tick_max_ratio = get_circular_tick_path_ratio_bounds(
                len(gb_record.seq),
                str((cfg or canvas_config._cfg).canvas.circular.track_type),
                bool((cfg or canvas_config._cfg).canvas.strandedness),
                tick_track_channel_override=tick_track_channel_override,
            )
            width_px = center * (max(float(tick_min_ratio), float(tick_max_ratio)) - min(float(tick_min_ratio), float(tick_max_ratio)))
        captured["ticks"] = (
            center,
            width_px,
        )
        return canvas

    def capture_numeric_slot(
        key: str,
        canvas_config,
        track_width_override,
        norm_factor_override,
    ) -> None:
        assert track_width_override is not None
        assert norm_factor_override is not None
        captured[key] = (
            float(norm_factor_override) * float(canvas_config.radius),
            float(track_width_override),
        )

    def fake_add_depth_group_on_canvas(
        canvas,
        gb_record,
        depth_df,
        canvas_config,
        depth_config,
        config_dict,
        *,
        track_width_override=None,
        norm_factor_override=None,
        group_id=None,
        cfg=None,
    ):
        capture_numeric_slot("depth", canvas_config, track_width_override, norm_factor_override)
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
        group_id=None,
        cfg=None,
    ):
        capture_numeric_slot("gc_content", canvas_config, track_width_override, norm_factor_override)
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
        group_id=None,
        cfg=None,
    ):
        capture_numeric_slot("gc_skew", canvas_config, track_width_override, norm_factor_override)
        return canvas

    monkeypatch.setattr(circular_assemble_module, "add_record_group_on_canvas", fake_add_record_group_on_canvas)
    monkeypatch.setattr(circular_assemble_module, "add_tick_group_on_canvas", fake_add_tick_group_on_canvas)
    monkeypatch.setattr(circular_assemble_module, "add_depth_group_on_canvas", fake_add_depth_group_on_canvas)
    monkeypatch.setattr(circular_assemble_module, "add_gc_content_group_on_canvas", fake_add_gc_content_group_on_canvas)
    monkeypatch.setattr(circular_assemble_module, "add_gc_skew_group_on_canvas", fake_add_gc_skew_group_on_canvas)

    slots = (
        default_circular_track_slots(show_depth=show_depth, show_gc=True, show_skew=True)
        if use_default_slots
        else None
    )
    kwargs = {}
    if use_default_slots:
        kwargs["circular_track_slots"] = slots
    if show_depth:
        kwargs.update({"depth_table": _depth_table(str(record.id)), "window": 100, "step": 100})

    assemble_circular_diagram_from_record(
        record,
        config_dict=config_dict,
        default_colors=default_colors,
        selected_features_set=SELECTED_FEATURES,
        legend="none",
        **kwargs,
    )
    return captured


def _assert_geometry_matches(
    observed: dict[str, tuple[float, float]],
    expected: dict[str, tuple[float, float]],
) -> None:
    assert observed.keys() == expected.keys()
    for key, expected_pair in expected.items():
        assert observed[key] == pytest.approx(expected_pair)


def test_parse_circular_track_slot_with_duplicate_renderer_params() -> None:
    slot = parse_circular_track_slot(
        "at_skew:dinucleotide_skew@nt=AT,w=24px,r=0.42,z=7,legend_label=AT skew"
    )

    assert slot.id == "at_skew"
    assert slot.renderer == "dinucleotide_skew"
    assert slot.params["nt"] == "AT"
    assert slot.params["legend_label"] == "AT skew"
    assert slot.width is not None
    assert slot.width.resolve(390) == 24
    assert slot.placement is not None
    assert slot.placement.radius is not None
    assert slot.placement.radius.resolve(390) == pytest.approx(163.8)
    assert slot.z == 7


def test_parse_circular_track_slots_rejects_duplicate_ids_and_unknown_renderer() -> None:
    with pytest.raises(TrackSpecParseError, match="duplicate circular track slot id"):
        parse_circular_track_slots(["gc_skew:dinucleotide_skew", "gc_skew:dinucleotide_skew@nt=AT"])

    with pytest.raises(TrackSpecParseError, match="unknown circular track renderer"):
        parse_circular_track_slots(["custom:not_a_renderer"])


def test_parse_circular_track_slots_normalizes_object_renderer_aliases() -> None:
    slots = parse_circular_track_slots([CircularTrackSlot(id="custom_skew", renderer="skew")])

    assert slots[0].renderer == "dinucleotide_skew"


def test_default_circular_track_slots_do_not_include_tick_axis_param() -> None:
    slots = default_circular_track_slots(show_depth=False, show_gc=True, show_skew=True)
    ticks = next(slot for slot in slots if slot.renderer == "ticks")

    assert "axis" not in ticks.params


def test_parse_circular_tick_slot_rejects_axis_param() -> None:
    with pytest.raises(TrackSpecParseError, match="ticks slots no longer accept 'axis'"):
        parse_circular_track_slot("ticks:ticks@axis=false")

    with pytest.raises(TrackSpecParseError, match="ticks slots no longer accept 'axis'"):
        parse_circular_track_slots(
            [CircularTrackSlot(id="ticks", renderer="ticks", params={"axis": True})]
        )


def test_no_custom_and_default_custom_slots_match_tuckin_separate_strands_geometry(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    no_custom = _capture_circular_core_geometry(
        monkeypatch,
        track_type="tuckin",
        strandedness=True,
        use_default_slots=False,
    )
    monkeypatch.undo()
    default_custom = _capture_circular_core_geometry(
        monkeypatch,
        track_type="tuckin",
        strandedness=True,
        use_default_slots=True,
    )

    _assert_geometry_matches(default_custom, no_custom)


def test_no_custom_and_default_custom_slots_match_tuckin_combined_strands_geometry(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    no_custom = _capture_circular_core_geometry(
        monkeypatch,
        track_type="tuckin",
        strandedness=False,
        use_default_slots=False,
    )
    monkeypatch.undo()
    default_custom = _capture_circular_core_geometry(
        monkeypatch,
        track_type="tuckin",
        strandedness=False,
        use_default_slots=True,
    )

    _assert_geometry_matches(default_custom, no_custom)


def test_no_custom_and_default_custom_slots_match_middle_geometry(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    no_custom = _capture_circular_core_geometry(
        monkeypatch,
        track_type="middle",
        use_default_slots=False,
    )
    monkeypatch.undo()
    default_custom = _capture_circular_core_geometry(
        monkeypatch,
        track_type="middle",
        use_default_slots=True,
    )

    _assert_geometry_matches(default_custom, no_custom)


def test_no_custom_and_default_custom_slots_match_spreadout_geometry(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    no_custom = _capture_circular_core_geometry(
        monkeypatch,
        track_type="spreadout",
        use_default_slots=False,
    )
    monkeypatch.undo()
    default_custom = _capture_circular_core_geometry(
        monkeypatch,
        track_type="spreadout",
        use_default_slots=True,
    )

    _assert_geometry_matches(default_custom, no_custom)


def test_no_custom_and_default_custom_slots_match_depth_gc_skew_compression(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    no_custom = _capture_circular_core_geometry(
        monkeypatch,
        track_type="middle",
        show_depth=True,
        use_default_slots=False,
    )
    monkeypatch.undo()
    default_custom = _capture_circular_core_geometry(
        monkeypatch,
        track_type="middle",
        show_depth=True,
        use_default_slots=True,
    )

    _assert_geometry_matches(default_custom, no_custom)


def test_resolve_circular_track_slots_preserves_legacy_defaults_in_compatibility_mode() -> None:
    context = CircularTrackLayoutContext(
        base_radius_px=390.0,
        legacy_centers_px={
            "features": 390.0,
            "ticks": 390.0,
            "gc_content": 249.6,
            "gc_skew": 171.6,
        },
        legacy_widths_px={
            "features": 37.05,
            "ticks": 7.8,
            "gc_content": 74.1,
            "gc_skew": 74.1,
        },
    )

    resolved = resolve_circular_track_slots(
        default_circular_track_slots(show_depth=False, show_gc=True, show_skew=True),
        context=context,
        compatibility_mode=True,
    )
    by_id = {slot.id: slot for slot in resolved}

    assert by_id["features"].center_radius_px == 390.0
    assert by_id["gc_content"].center_radius_px == 249.6
    assert by_id["gc_content"].draw_width_px == pytest.approx(74.1)
    assert by_id["gc_skew"].center_radius_px == 171.6
    assert by_id["gc_skew"].draw_width_px == pytest.approx(74.1)


def test_resolve_circular_track_slots_packs_measured_footprints_without_overlap() -> None:
    context = CircularTrackLayoutContext(
        base_radius_px=390.0,
        legacy_centers_px={"features": 390.0, "ticks": 390.0},
        legacy_widths_px={"features": 40.0, "ticks": 12.0, "gc_content": 30.0},
        default_gap_px=4.0,
        feature_band_offsets_px=(-55.0, -5.0),
        tick_path_ratio_bounds=(0.98, 1.0),
        tick_label_offsets_px=(20.0, 80.0),
    )

    resolved = resolve_circular_track_slots(
        [
            CircularTrackSlot(id="features", renderer="features"),
            CircularTrackSlot(
                id="ticks",
                renderer="ticks",
                params={"label_side": "legacy", "tick_side": "legacy"},
            ),
            CircularTrackSlot(id="gc_content", renderer="dinucleotide_content"),
        ],
        context=context,
    )

    annuli = [
        (slot.reserved_inner_radius_px, slot.reserved_outer_radius_px)
        for slot in resolved
    ]
    assert resolved[0].anchor_radius_px == 390.0
    assert resolved[0].draw_inner_radius_px == 335.0
    for idx, current in enumerate(annuli):
        for other in annuli[idx + 1:]:
            assert current[0] >= other[1] or other[0] >= current[1]


def test_measured_packer_places_legacy_ticks_by_reserved_footprint_not_anchor() -> None:
    context = CircularTrackLayoutContext(
        base_radius_px=500.0,
        legacy_centers_px={"ticks": 500.0},
        legacy_widths_px={"ticks": 10.0, "gc_content": 40.0},
        preferred_layouts_px={"gc_content": (300.0, 40.0)},
        default_gap_px=5.0,
        tick_path_ratio_bounds=(0.84, 0.86),
    )

    resolved = resolve_circular_track_slots(
        [
            CircularTrackSlot(id="gc_content", renderer="dinucleotide_content"),
            CircularTrackSlot(
                id="ticks",
                renderer="ticks",
                params={"label_side": "legacy", "tick_side": "legacy"},
            ),
        ],
        context=context,
    )
    by_id = {slot.id: slot for slot in resolved}
    cursor_outer = by_id["gc_content"].reserved_inner_radius_px - context.default_gap_px

    assert by_id["ticks"].anchor_radius_px > cursor_outer
    assert by_id["ticks"].reserved_outer_radius_px == pytest.approx(cursor_outer)
    assert by_id["gc_content"].reserved_inner_radius_px - by_id["ticks"].reserved_outer_radius_px == pytest.approx(5.0)


def test_reordered_builtin_gc_ticks_skew_keeps_gc_widths() -> None:
    context = CircularTrackLayoutContext(
        base_radius_px=500.0,
        legacy_centers_px={"features": 500.0, "ticks": 500.0},
        legacy_widths_px={
            "features": 60.0,
            "ticks": 10.0,
            "gc_content": 50.0,
            "gc_skew": 50.0,
        },
        preferred_layouts_px={
            "gc_content": (320.0, 50.0),
            "gc_skew": (230.0, 50.0),
        },
        default_gap_px=5.0,
        feature_band_offsets_px=(-70.0, -10.0),
        tick_path_ratio_bounds=(0.84, 0.86),
    )

    resolved = resolve_circular_track_slots(
        [
            CircularTrackSlot(id="features", renderer="features"),
            CircularTrackSlot(id="gc_content", renderer="dinucleotide_content"),
            CircularTrackSlot(
                id="ticks",
                renderer="ticks",
                params={"label_side": "legacy", "tick_side": "legacy"},
            ),
            CircularTrackSlot(id="gc_skew", renderer="dinucleotide_skew"),
        ],
        context=context,
    )
    by_id = {slot.id: slot for slot in resolved}

    assert by_id["gc_content"].draw_width_px == pytest.approx(50.0)
    assert by_id["gc_skew"].draw_width_px == pytest.approx(50.0)
    assert by_id["gc_content"].reserved_inner_radius_px - by_id["ticks"].reserved_outer_radius_px == pytest.approx(5.0)


def test_legacy_tick_reserved_label_bounds_are_remeasured_for_moved_anchor() -> None:
    context = CircularTrackLayoutContext(
        base_radius_px=500.0,
        legacy_centers_px={"ticks": 500.0},
        legacy_widths_px={"ticks": 10.0, "gc_content": 40.0},
        preferred_layouts_px={"gc_content": (300.0, 40.0)},
        default_gap_px=5.0,
        tick_path_ratio_bounds=(0.84, 0.86),
        tick_label_radius_ratio=0.82,
        tick_label_extent_px=6.0,
        tick_label_offsets_px=(100.0, 120.0),
    )

    resolved = resolve_circular_track_slots(
        [
            CircularTrackSlot(id="gc_content", renderer="dinucleotide_content"),
            CircularTrackSlot(
                id="ticks",
                renderer="ticks",
                params={"label_side": "legacy", "tick_side": "legacy"},
            ),
        ],
        context=context,
    )
    ticks = {slot.id: slot for slot in resolved}["ticks"]
    label_center = ticks.anchor_radius_px * 0.82

    assert ticks.reserved_inner_radius_px == pytest.approx(min(ticks.draw_inner_radius_px, label_center - 6.0))
    assert ticks.reserved_outer_radius_px == pytest.approx(max(ticks.draw_outer_radius_px, label_center + 6.0))
    assert ticks.reserved_outer_radius_px < ticks.anchor_radius_px + 100.0


def test_tick_slot_footprint_does_not_reserve_axis_padding() -> None:
    context = CircularTrackLayoutContext(
        base_radius_px=100.0,
        legacy_centers_px={"ticks": 100.0},
        legacy_widths_px={"ticks": 0.0},
    )

    resolved = resolve_circular_track_slots(
        [parse_circular_track_slot("ticks:ticks@label_side=none,tick_side=none")],
        context=context,
        compatibility_mode=True,
    )

    assert resolved[0].draw_inner_radius_px == pytest.approx(100.0)
    assert resolved[0].draw_outer_radius_px == pytest.approx(100.0)
    assert resolved[0].reserved_inner_radius_px == pytest.approx(100.0)
    assert resolved[0].reserved_outer_radius_px == pytest.approx(100.0)


def test_tick_label_footprint_repacks_data_track() -> None:
    context = CircularTrackLayoutContext(
        base_radius_px=100.0,
        legacy_widths_px={"gc_content": 10.0},
        default_gap_px=0.0,
        tick_font_size_px=10.0,
    )

    resolved = resolve_circular_track_slots(
        [
            parse_circular_track_slot("ticks:ticks@r=100px,label_side=inside,tick_side=none"),
            CircularTrackSlot(id="gc_content", renderer="dinucleotide_content"),
        ],
        context=context,
    )
    by_id = {slot.id: slot for slot in resolved}

    assert by_id["ticks"].reserved_inner_radius_px == pytest.approx(82.0)
    assert by_id["ticks"].draw_inner_radius_px == pytest.approx(100.0)
    assert by_id["gc_content"].reserved_outer_radius_px <= by_id["ticks"].reserved_inner_radius_px
    assert by_id["gc_content"].center_radius_px == pytest.approx(77.0)


def test_tick_label_hard_context_remains_compatible() -> None:
    context = CircularTrackLayoutContext(
        base_radius_px=100.0,
        legacy_widths_px={"gc_content": 10.0},
        default_gap_px=0.0,
        tick_font_size_px=10.0,
        tick_labels_hard=True,
    )

    resolved = resolve_circular_track_slots(
        [
            parse_circular_track_slot("ticks:ticks@r=100px,label_side=inside,tick_side=none"),
            CircularTrackSlot(id="gc_content", renderer="dinucleotide_content"),
        ],
        context=context,
    )
    by_id = {slot.id: slot for slot in resolved}

    assert by_id["ticks"].reserved_inner_radius_px == pytest.approx(82.0)
    assert by_id["gc_content"].reserved_outer_radius_px <= by_id["ticks"].reserved_inner_radius_px
    assert by_id["gc_content"].center_radius_px == pytest.approx(77.0)


def test_legacy_tick_label_bounds_keep_margin_from_tick_path_when_radius_shrinks() -> None:
    center_radius = 260.0
    total_len = 5_528_445
    font_size = 14.0
    tick_width = 2.0

    tick_inner, _tick_outer = get_circular_tick_path_radius_bounds(
        center_radius,
        total_len,
        "large",
        "tuckin",
        True,
        tick_side="legacy",
    )
    label_bounds = get_circular_tick_label_radius_bounds(
        center_radius_px=center_radius,
        total_len=total_len,
        track_type="tuckin",
        strandedness=True,
        font_size=font_size,
        font_family="Liberation Sans",
        dpi=72,
        tick_width=tick_width,
    )

    assert label_bounds is not None
    expected_margin = max(2.0, font_size * 0.15) + (tick_width / 2.0)
    assert label_bounds[1] <= tick_inner - expected_margin + 1e-6


def test_generated_tick_label_arc_uses_resolved_geometry_radius() -> None:
    center_radius = 260.0
    total_len = 5_528_445
    tick = 1_000_000
    expected = resolve_circular_tick_label_geometry(
        center_radius_px=center_radius,
        total_len=total_len,
        size="large",
        tick=tick,
        label_text="1.0 Mbp",
        font_size=14.0,
        font_family="Liberation Sans",
        track_type="tuckin",
        strandedness=True,
        dpi=72,
        tick_width=2.0,
    )

    elements = generate_circular_tick_labels(
        center_radius,
        total_len,
        "large",
        [tick],
        "none",
        "black",
        14.0,
        "normal",
        "Liberation Sans",
        "tuckin",
        True,
        72,
        tick_width=2.0,
    )

    assert elements
    match = re.search(r"A([0-9.]+),", elements[0].tostring())
    assert match is not None
    assert float(match.group(1)) == pytest.approx(expected.path_radius_px)
    _expected_anchor, expected_baseline = set_tick_label_anchor_value(total_len, tick)
    text_svg = elements[1].tostring()
    assert 'text-anchor="middle"' in text_svg
    assert f'dominant-baseline="{expected_baseline}"' in text_svg
    assert 'startOffset="50%"' in text_svg


def test_generated_tick_labels_apply_position_dependent_textpath_settings() -> None:
    total_len = 5_528_445
    elements = generate_circular_tick_labels(
        260.0,
        total_len,
        "large",
        [5_000_000, 3_000_000],
        "none",
        "black",
        14.0,
        "normal",
        "Liberation Sans",
        "tuckin",
        True,
        72,
        tick_width=2.0,
    )

    five_mbp_svg = elements[1].tostring()
    three_mbp_svg = elements[3].tostring()
    assert 'text-anchor="middle"' in five_mbp_svg
    assert 'dominant-baseline="text-after-edge"' in five_mbp_svg
    assert 'text-anchor="middle"' in three_mbp_svg
    assert 'dominant-baseline="hanging"' in three_mbp_svg


def test_generated_tick_labels_remain_centered_on_left_side_ticks() -> None:
    total_len = 5_528_445
    elements = generate_circular_tick_labels(
        260.0,
        total_len,
        "large",
        [4_500_000],
        "none",
        "black",
        14.0,
        "normal",
        "Liberation Sans",
        "tuckin",
        True,
        72,
        tick_width=2.0,
    )

    text_svg = elements[1].tostring()
    assert 'text-anchor="middle"' in text_svg
    assert 'startOffset="50%"' in text_svg


def test_legacy_tick_label_path_radius_is_position_independent_for_equal_track() -> None:
    center_radius = 260.0
    total_len = 5_528_445
    common_kwargs = {
        "center_radius_px": center_radius,
        "total_len": total_len,
        "size": "large",
        "font_size": 14.0,
        "font_family": "Liberation Sans",
        "track_type": "tuckin",
        "strandedness": True,
        "dpi": 72,
        "tick_width": 2.0,
    }

    five_mbp = resolve_circular_tick_label_geometry(
        tick=5_000_000,
        label_text="5.0 Mbp",
        **common_kwargs,
    )
    three_mbp = resolve_circular_tick_label_geometry(
        tick=3_000_000,
        label_text="3.0 Mbp",
        **common_kwargs,
    )

    assert five_mbp.path_radius_px == pytest.approx(three_mbp.path_radius_px)
    assert five_mbp.radial_inner_px == pytest.approx(three_mbp.radial_inner_px)
    assert five_mbp.radial_outer_px == pytest.approx(three_mbp.radial_outer_px)


def test_resolve_circular_track_slots_auto_slots_avoid_definition_reserved_band() -> None:
    context = CircularTrackLayoutContext(
        base_radius_px=100.0,
        legacy_centers_px={"gc_content": 80.0, "gc_skew": 45.0},
        legacy_widths_px={"gc_content": 20.0, "gc_skew": 20.0},
        default_gap_px=5.0,
        reserved_bands_px=((0.0, 35.0),),
        min_auto_inner_radius_px=35.0,
    )

    resolved = resolve_circular_track_slots(
        [
            CircularTrackSlot(id="gc_content", renderer="dinucleotide_content"),
            CircularTrackSlot(id="gc_skew", renderer="dinucleotide_skew"),
        ],
        context=context,
    )

    assert min(slot.reserved_inner_radius_px for slot in resolved) >= 35.0 - 1e-6


def test_resolve_circular_track_slots_compresses_implicit_numeric_widths_to_fit_definition_guard(
    caplog: pytest.LogCaptureFixture,
) -> None:
    context = CircularTrackLayoutContext(
        base_radius_px=100.0,
        legacy_centers_px={"gc_content": 80.0, "gc_skew": 45.0},
        legacy_widths_px={"gc_content": 30.0, "gc_skew": 30.0},
        default_gap_px=5.0,
        reserved_bands_px=((0.0, 45.0),),
        min_auto_inner_radius_px=45.0,
    )
    caplog.set_level(logging.INFO, logger="gbdraw.diagrams.circular.slot_layout")

    resolved = resolve_circular_track_slots(
        [
            CircularTrackSlot(id="gc_content", renderer="dinucleotide_content"),
            CircularTrackSlot(id="gc_skew", renderer="dinucleotide_skew"),
        ],
        context=context,
    )
    by_id = {slot.id: slot for slot in resolved}

    assert by_id["gc_content"].draw_width_px == pytest.approx(by_id["gc_skew"].draw_width_px)
    assert by_id["gc_content"].draw_width_px < 30.0
    assert min(slot.reserved_inner_radius_px for slot in resolved) >= 45.0 - 1e-6
    assert "Auto-compressed circular numeric track widths" in caplog.text


def test_resolve_circular_track_slots_compresses_pinned_implicit_numeric_width_without_moving_center() -> None:
    context = CircularTrackLayoutContext(
        base_radius_px=100.0,
        legacy_centers_px={"gc_skew": 45.0},
        legacy_widths_px={"gc_content": 30.0, "gc_skew": 30.0},
        default_gap_px=5.0,
        reserved_bands_px=((0.0, 45.0),),
        min_auto_inner_radius_px=45.0,
    )

    resolved = resolve_circular_track_slots(
        [
            parse_circular_track_slot("gc_content:dinucleotide_content@r=80px"),
            CircularTrackSlot(id="gc_skew", renderer="dinucleotide_skew"),
        ],
        context=context,
    )
    by_id = {slot.id: slot for slot in resolved}

    assert by_id["gc_content"].center_radius_px == pytest.approx(80.0)
    assert by_id["gc_content"].draw_width_px == pytest.approx(by_id["gc_skew"].draw_width_px)
    assert by_id["gc_content"].draw_width_px < 30.0
    assert min(slot.reserved_inner_radius_px for slot in resolved) >= 45.0 - 1e-6


def test_resolve_circular_track_slots_warns_for_pinned_definition_overlap(caplog: pytest.LogCaptureFixture) -> None:
    slot = parse_circular_track_slot("gc_skew:dinucleotide_skew@r=30px,w=20px")
    context = CircularTrackLayoutContext(
        base_radius_px=100.0,
        reserved_bands_px=((0.0, 40.0),),
        min_auto_inner_radius_px=40.0,
    )

    resolved = resolve_circular_track_slots([slot], context=context)

    assert resolved[0].center_radius_px == pytest.approx(30.0)
    assert "overlaps a reserved circular layout band" in caplog.text


def test_resolve_circular_track_slots_raises_when_definition_guard_leaves_no_space() -> None:
    context = CircularTrackLayoutContext(
        base_radius_px=100.0,
        legacy_centers_px={"gc_content": 99.0, "gc_skew": 97.0},
        legacy_widths_px={"gc_content": 30.0, "gc_skew": 30.0},
        default_gap_px=5.0,
        reserved_bands_px=((0.0, 98.0),),
        min_auto_inner_radius_px=98.0,
    )

    with pytest.raises(ValueError, match="cannot fit without overlapping the center definition"):
        resolve_circular_track_slots(
            [
                CircularTrackSlot(id="gc_content", renderer="dinucleotide_content"),
                CircularTrackSlot(id="gc_skew", renderer="dinucleotide_skew"),
            ],
            context=context,
        )


def test_resolve_circular_track_slots_does_not_compress_explicit_numeric_widths() -> None:
    context = CircularTrackLayoutContext(
        base_radius_px=100.0,
        legacy_centers_px={"gc_content": 80.0, "gc_skew": 45.0},
        legacy_widths_px={"gc_content": 30.0, "gc_skew": 30.0},
        default_gap_px=5.0,
        reserved_bands_px=((0.0, 45.0),),
        min_auto_inner_radius_px=45.0,
    )

    with pytest.raises(ValueError, match="cannot fit without overlapping the center definition"):
        resolve_circular_track_slots(
            [
                parse_circular_track_slot("gc_content:dinucleotide_content@w=30px"),
                parse_circular_track_slot("gc_skew:dinucleotide_skew@w=30px"),
            ],
            context=context,
        )


def test_resolve_circular_track_slots_does_not_compress_matching_track_spec_widths() -> None:
    context = CircularTrackLayoutContext(
        base_radius_px=100.0,
        legacy_centers_px={"gc_content": 80.0, "gc_skew": 45.0},
        legacy_widths_px={"gc_content": 30.0, "gc_skew": 30.0},
        default_gap_px=5.0,
        reserved_bands_px=((0.0, 45.0),),
        min_auto_inner_radius_px=45.0,
    )
    track_specs = [
        TrackSpec(
            id="gc_content",
            kind="gc_content",
            mode="circular",
            placement=CircularTrackPlacement(width=ScalarSpec(30.0, "px")),
        ),
        TrackSpec(
            id="gc_skew",
            kind="gc_skew",
            mode="circular",
            placement=CircularTrackPlacement(width=ScalarSpec(30.0, "px")),
        ),
    ]

    with pytest.raises(ValueError, match="cannot fit without overlapping the center definition"):
        resolve_circular_track_slots(
            [
                CircularTrackSlot(id="gc_content", renderer="dinucleotide_content"),
                CircularTrackSlot(id="gc_skew", renderer="dinucleotide_skew"),
            ],
            context=context,
            legacy_track_specs=track_specs,
        )


def test_resolve_circular_track_slots_applies_preferred_layouts_by_slot_id_only() -> None:
    context = CircularTrackLayoutContext(
        base_radius_px=100.0,
        legacy_widths_px={"gc_skew": 20.0},
        preferred_layouts_px={"gc_skew": (90.0, 20.0)},
        default_gap_px=5.0,
    )

    resolved = resolve_circular_track_slots(
        [
            CircularTrackSlot(id="gc_skew", renderer="dinucleotide_skew"),
            CircularTrackSlot(id="at_skew", renderer="dinucleotide_skew"),
        ],
        context=context,
    )
    by_id = {slot.id: slot for slot in resolved}

    assert by_id["gc_skew"].center_radius_px == pytest.approx(90.0)
    assert by_id["at_skew"].center_radius_px == pytest.approx(65.0)


def test_resolve_circular_track_slots_auto_numeric_order_ignores_legacy_id_centers() -> None:
    context = CircularTrackLayoutContext(
        base_radius_px=100.0,
        legacy_centers_px={"gc_content": 80.0, "gc_skew": 45.0},
        legacy_widths_px={"gc_content": 20.0, "gc_skew": 20.0},
        default_gap_px=5.0,
        auto_start_radius_px=100.0,
    )

    resolved = resolve_circular_track_slots(
        [
            CircularTrackSlot(id="gc_skew", renderer="dinucleotide_skew"),
            CircularTrackSlot(id="gc_content", renderer="dinucleotide_content"),
        ],
        context=context,
    )

    assert [slot.id for slot in resolved] == ["gc_skew", "gc_content"]
    assert resolved[0].center_radius_px > resolved[1].center_radius_px
    assert resolved[0].center_radius_px == pytest.approx(90.0)
    assert resolved[1].center_radius_px == pytest.approx(65.0)


def test_auto_gap_fallback_does_not_place_later_slot_outside_previous_auto_slot() -> None:
    context = CircularTrackLayoutContext(
        base_radius_px=100.0,
        default_gap_px=0.0,
        feature_band_offsets_px=(-15.0, 0.0),
        reserved_bands_px=((0.0, 20.0),),
        min_auto_inner_radius_px=20.0,
    )

    with pytest.raises(ValueError, match="cannot fit without overlapping the center definition"):
        resolve_circular_track_slots(
            [
                parse_circular_track_slot("features:features@r=100px"),
                parse_circular_track_slot("ticks:ticks@r=60px,w=10px,label_side=none,tick_side=inside"),
                parse_circular_track_slot("gc_content:dinucleotide_content@w=25px"),
                parse_circular_track_slot("gc_skew:dinucleotide_skew@w=25px"),
            ],
            context=context,
        )


def test_ordered_pack_compression_preserves_order() -> None:
    context = CircularTrackLayoutContext(
        base_radius_px=100.0,
        legacy_widths_px={"gc_content": 25.0, "gc_skew": 25.0},
        default_gap_px=0.0,
        feature_band_offsets_px=(-15.0, 0.0),
        reserved_bands_px=((0.0, 20.0),),
        min_auto_inner_radius_px=20.0,
    )

    resolved = resolve_circular_track_slots(
        [
            parse_circular_track_slot("features:features@r=100px"),
            parse_circular_track_slot("ticks:ticks@r=60px,w=10px,label_side=none,tick_side=inside"),
            CircularTrackSlot(id="gc_content", renderer="dinucleotide_content"),
            CircularTrackSlot(id="gc_skew", renderer="dinucleotide_skew"),
        ],
        context=context,
    )
    by_id = {slot.id: slot for slot in resolved}

    assert by_id["gc_content"].center_radius_px > by_id["gc_skew"].center_radius_px
    assert by_id["gc_content"].draw_width_px == pytest.approx(by_id["gc_skew"].draw_width_px)
    assert by_id["gc_content"].draw_width_px < 25.0
    assert by_id["gc_skew"].reserved_inner_radius_px >= 20.0 - 1e-6


def test_explicit_radius_can_override_toolbar_order() -> None:
    context = CircularTrackLayoutContext(
        base_radius_px=100.0,
        legacy_widths_px={"gc_content": 20.0},
        default_gap_px=5.0,
    )

    resolved = resolve_circular_track_slots(
        [
            CircularTrackSlot(id="gc_content", renderer="dinucleotide_content"),
            parse_circular_track_slot("gc_skew:dinucleotide_skew@r=90px,w=10px"),
        ],
        context=context,
    )
    by_id = {slot.id: slot for slot in resolved}

    assert by_id["gc_skew"].center_radius_px > by_id["gc_content"].center_radius_px
    assert by_id["gc_skew"].center_radius_px == pytest.approx(90.0)


@pytest.mark.circular
def test_default_custom_slots_tuckin_preserve_numeric_order_with_feature_footprint(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import gbdraw.diagrams.circular.assemble as circular_assemble_module

    record = _load_edl933_record()
    config_dict = modify_config_dict(
        load_config_toml("gbdraw.data", "config.toml"),
        show_labels=False,
        show_gc=True,
        show_skew=True,
        track_type="tuckin",
        strandedness=True,
    )
    default_colors = load_default_colors("", palette="default")
    captured: dict[str, tuple[float, float]] = {}

    def capture_numeric_slot(
        slot_id: str,
        canvas_config,
        track_width_override,
        norm_factor_override,
    ) -> None:
        assert track_width_override is not None
        assert norm_factor_override is not None
        captured[slot_id] = (
            float(norm_factor_override) * float(canvas_config.radius),
            float(track_width_override),
        )

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
        group_id=None,
        cfg=None,
    ):
        capture_numeric_slot(str(group_id or "gc_content"), canvas_config, track_width_override, norm_factor_override)
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
        group_id=None,
        cfg=None,
    ):
        capture_numeric_slot(str(group_id or "gc_skew"), canvas_config, track_width_override, norm_factor_override)
        return canvas

    monkeypatch.setattr(circular_assemble_module, "add_gc_content_group_on_canvas", fake_add_gc_content_group_on_canvas)
    monkeypatch.setattr(circular_assemble_module, "add_gc_skew_group_on_canvas", fake_add_gc_skew_group_on_canvas)

    assemble_circular_diagram_from_record(
        record,
        config_dict=config_dict,
        default_colors=default_colors,
        selected_features_set=SELECTED_FEATURES,
        legend="none",
        circular_track_slots=default_circular_track_slots(show_depth=False, show_gc=True, show_skew=True),
    )

    assert captured["gc_content"][0] > captured["gc_skew"][0]


@pytest.mark.circular
@pytest.mark.parametrize("track_type", ["tuckin", "middle", "spreadout"])
def test_default_custom_slots_use_ordered_legacy_numeric_lanes(
    monkeypatch: pytest.MonkeyPatch,
    track_type: str,
) -> None:
    import gbdraw.diagrams.circular.assemble as circular_assemble_module

    record = _load_record()
    config_dict = _base_config(track_type=track_type)
    default_colors = load_default_colors("", palette="default")
    captured: dict[str, tuple[float, float]] = {}

    def capture_numeric_slot(
        slot_id: str,
        canvas_config,
        track_width_override,
        norm_factor_override,
    ) -> None:
        assert track_width_override is not None
        assert norm_factor_override is not None
        captured[slot_id] = (
            float(norm_factor_override) * float(canvas_config.radius),
            float(track_width_override),
        )

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
        group_id=None,
        cfg=None,
    ):
        capture_numeric_slot(str(group_id or "gc_content"), canvas_config, track_width_override, norm_factor_override)
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
        group_id=None,
        cfg=None,
    ):
        capture_numeric_slot(str(group_id or "gc_skew"), canvas_config, track_width_override, norm_factor_override)
        return canvas

    monkeypatch.setattr(circular_assemble_module, "add_gc_content_group_on_canvas", fake_add_gc_content_group_on_canvas)
    monkeypatch.setattr(circular_assemble_module, "add_gc_skew_group_on_canvas", fake_add_gc_skew_group_on_canvas)

    assemble_circular_diagram_from_record(
        record,
        config_dict=config_dict,
        default_colors=default_colors,
        selected_features_set=SELECTED_FEATURES,
        legend="none",
        circular_track_slots=default_circular_track_slots(show_depth=False, show_gc=True, show_skew=True),
    )

    assert captured["gc_content"][0] > captured["gc_skew"][0]
    assert captured["gc_content"][1] == pytest.approx(captured["gc_skew"][1])


@pytest.mark.circular
def test_reordered_builtin_numeric_slots_follow_slot_order(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import gbdraw.diagrams.circular.assemble as circular_assemble_module

    record = _load_record()
    config_dict = _base_config(track_type="middle")
    default_colors = load_default_colors("", palette="default")
    captured: dict[str, tuple[float, float]] = {}

    def capture_numeric_slot(
        slot_id: str,
        canvas_config,
        track_width_override,
        norm_factor_override,
    ) -> None:
        assert track_width_override is not None
        assert norm_factor_override is not None
        captured[slot_id] = (
            float(norm_factor_override) * float(canvas_config.radius),
            float(track_width_override),
        )

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
        group_id=None,
        cfg=None,
    ):
        capture_numeric_slot(str(group_id or "gc_content"), canvas_config, track_width_override, norm_factor_override)
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
        group_id=None,
        cfg=None,
    ):
        capture_numeric_slot(str(group_id or "gc_skew"), canvas_config, track_width_override, norm_factor_override)
        return canvas

    monkeypatch.setattr(circular_assemble_module, "add_gc_content_group_on_canvas", fake_add_gc_content_group_on_canvas)
    monkeypatch.setattr(circular_assemble_module, "add_gc_skew_group_on_canvas", fake_add_gc_skew_group_on_canvas)

    assemble_circular_diagram_from_record(
        record,
        config_dict=config_dict,
        default_colors=default_colors,
        selected_features_set=SELECTED_FEATURES,
        legend="none",
        circular_track_slots=[
            CircularTrackSlot(id="features", renderer="features"),
            CircularTrackSlot(
                id="ticks",
                renderer="ticks",
                params={"placement": "legacy_axis", "label_side": "legacy", "tick_side": "legacy"},
            ),
            CircularTrackSlot(id="gc_skew", renderer="dinucleotide_skew"),
            CircularTrackSlot(id="gc_content", renderer="dinucleotide_content"),
        ],
    )

    assert captured["gc_skew"][0] > captured["gc_content"][0]


@pytest.mark.circular
def test_edl933_reordered_gc_ticks_skew_uses_measured_tick_footprint(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import gbdraw.diagrams.circular.assemble as circular_assemble_module

    record = _load_edl933_record()
    config_dict = modify_config_dict(
        load_config_toml("gbdraw.data", "config.toml"),
        show_labels=False,
        show_gc=True,
        show_skew=True,
        track_type="tuckin",
        strandedness=True,
    )
    default_colors = load_default_colors("", palette="default")
    captured: dict[str, float | tuple[float, float]] = {}

    def fake_add_record_group_on_canvas(canvas, *args, **kwargs):
        return canvas

    def capture_numeric_slot(
        slot_id: str,
        canvas_config,
        track_width_override,
        norm_factor_override,
        cfg,
    ) -> None:
        assert track_width_override is not None
        assert norm_factor_override is not None
        assert cfg is not None
        captured[slot_id] = (
            float(norm_factor_override) * float(canvas_config.radius),
            float(track_width_override),
        )
        captured["default_gc_width"] = (
            float(canvas_config.radius)
            * float(canvas_config.track_ratio)
            * float(cfg.canvas.circular.track_ratio_factors[str(canvas_config.length_param)][1])
        )
        captured["default_gap"] = max(1.0, 0.01 * float(canvas_config.radius))

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
        group_id=None,
        cfg=None,
    ):
        capture_numeric_slot(str(group_id or "gc_content"), canvas_config, track_width_override, norm_factor_override, cfg)
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
        group_id=None,
        cfg=None,
    ):
        capture_numeric_slot(str(group_id or "gc_skew"), canvas_config, track_width_override, norm_factor_override, cfg)
        return canvas

    def fake_add_tick_group_on_canvas(
        canvas,
        gb_record,
        canvas_config,
        config_dict,
        *,
        radius_override=None,
        tick_track_channel_override=None,
        label_side="legacy",
        tick_side="legacy",
        tick_length_px=None,
        cfg=None,
    ):
        assert cfg is not None
        center = float(radius_override if radius_override is not None else canvas_config.radius)
        tick_min_ratio, tick_max_ratio = get_circular_tick_path_ratio_bounds(
            len(gb_record.seq),
            str(cfg.canvas.circular.track_type),
            bool(cfg.canvas.strandedness),
            tick_track_channel_override=tick_track_channel_override,
        )
        captured["ticks"] = (center, center * max(float(tick_min_ratio), float(tick_max_ratio)))
        return canvas

    monkeypatch.setattr(circular_assemble_module, "add_record_group_on_canvas", fake_add_record_group_on_canvas)
    monkeypatch.setattr(circular_assemble_module, "add_gc_content_group_on_canvas", fake_add_gc_content_group_on_canvas)
    monkeypatch.setattr(circular_assemble_module, "add_gc_skew_group_on_canvas", fake_add_gc_skew_group_on_canvas)
    monkeypatch.setattr(circular_assemble_module, "add_tick_group_on_canvas", fake_add_tick_group_on_canvas)
    monkeypatch.setattr(circular_assemble_module, "_definition_reserved_radius_px", lambda *args, **kwargs: 0.0)

    assemble_circular_diagram_from_record(
        record,
        config_dict=config_dict,
        default_colors=default_colors,
        selected_features_set=SELECTED_FEATURES,
        legend="none",
        circular_track_slots=[
            CircularTrackSlot(id="features", renderer="features"),
            CircularTrackSlot(id="gc_content", renderer="dinucleotide_content", params={"nt": "GC"}),
            CircularTrackSlot(
                id="ticks",
                renderer="ticks",
                params={"label_side": "legacy", "tick_side": "legacy"},
            ),
            CircularTrackSlot(id="gc_skew", renderer="dinucleotide_skew", params={"nt": "GC"}),
        ],
    )

    gc_center, gc_width = captured["gc_content"]  # type: ignore[misc]
    _, ticks_reserved_outer = captured["ticks"]  # type: ignore[misc]

    assert gc_width == pytest.approx(float(captured["default_gc_width"]))
    assert captured["gc_skew"][1] == pytest.approx(float(captured["default_gc_width"]))  # type: ignore[index]
    assert (float(gc_center) - (0.5 * float(gc_width))) - float(ticks_reserved_outer) == pytest.approx(
        float(captured["default_gap"]),
        abs=0.05,
    )


@pytest.mark.circular
def test_default_custom_slots_with_depth_use_outer_to_inner_numeric_lanes(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import gbdraw.diagrams.circular.assemble as circular_assemble_module

    record = _load_record()
    config_dict = _base_config(track_type="middle")
    default_colors = load_default_colors("", palette="default")
    depth_table = _depth_table(str(record.id))
    captured: dict[str, tuple[float, float]] = {}

    def capture_numeric_slot(
        slot_id: str,
        canvas_config,
        track_width_override,
        norm_factor_override,
    ) -> None:
        assert track_width_override is not None
        assert norm_factor_override is not None
        captured[slot_id] = (
            float(norm_factor_override) * float(canvas_config.radius),
            float(track_width_override),
        )

    def fake_add_depth_group_on_canvas(
        canvas,
        gb_record,
        depth_df,
        canvas_config,
        depth_config,
        config_dict,
        *,
        track_width_override=None,
        norm_factor_override=None,
        group_id=None,
        cfg=None,
    ):
        capture_numeric_slot(str(group_id or "depth"), canvas_config, track_width_override, norm_factor_override)
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
        group_id=None,
        cfg=None,
    ):
        capture_numeric_slot(str(group_id or "gc_content"), canvas_config, track_width_override, norm_factor_override)
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
        group_id=None,
        cfg=None,
    ):
        capture_numeric_slot(str(group_id or "gc_skew"), canvas_config, track_width_override, norm_factor_override)
        return canvas

    monkeypatch.setattr(circular_assemble_module, "add_depth_group_on_canvas", fake_add_depth_group_on_canvas)
    monkeypatch.setattr(circular_assemble_module, "add_gc_content_group_on_canvas", fake_add_gc_content_group_on_canvas)
    monkeypatch.setattr(circular_assemble_module, "add_gc_skew_group_on_canvas", fake_add_gc_skew_group_on_canvas)

    assemble_circular_diagram_from_record(
        record,
        config_dict=config_dict,
        default_colors=default_colors,
        selected_features_set=SELECTED_FEATURES,
        legend="none",
        depth_table=depth_table,
        window=100,
        step=100,
        circular_track_slots=default_circular_track_slots(show_depth=True, show_gc=True, show_skew=True),
    )

    assert captured["depth"][0] > captured["gc_content"][0] > captured["gc_skew"][0]
    assert captured["depth"][1] < captured["gc_content"][1]
    assert captured["gc_content"][1] == pytest.approx(captured["gc_skew"][1])


@pytest.mark.circular
def test_custom_duplicate_skew_with_depth_tuckin_avoids_definition(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import gbdraw.diagrams.circular.assemble as circular_assemble_module

    record = _load_record()
    config_dict = _base_config(track_type="tuckin")
    default_colors = load_default_colors("", palette="default")
    depth_table = _depth_table(str(record.id))
    captured: dict[str, tuple[float, float] | float] = {}

    def capture_numeric_slot(
        slot_id: str,
        canvas_config,
        track_width_override,
        norm_factor_override,
    ) -> None:
        assert track_width_override is not None
        assert norm_factor_override is not None
        captured[slot_id] = (
            float(norm_factor_override) * float(canvas_config.radius),
            float(track_width_override),
        )

    def fake_add_depth_group_on_canvas(
        canvas,
        gb_record,
        depth_df,
        canvas_config,
        depth_config,
        config_dict,
        *,
        track_width_override=None,
        norm_factor_override=None,
        group_id=None,
        cfg=None,
    ):
        capture_numeric_slot(str(group_id or "depth"), canvas_config, track_width_override, norm_factor_override)
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
        group_id=None,
        cfg=None,
    ):
        capture_numeric_slot(str(group_id or "gc_content"), canvas_config, track_width_override, norm_factor_override)
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
        group_id=None,
        cfg=None,
    ):
        capture_numeric_slot(str(group_id or "gc_skew"), canvas_config, track_width_override, norm_factor_override)
        if "definition_reserved" not in captured:
            assert cfg is not None
            captured["definition_reserved"] = circular_assemble_module._definition_reserved_radius_px(
                gb_record,
                canvas_config,
                None,
                None,
                config_dict,
                cfg=cfg,
                plot_title=None,
                definition_profile="full",
            )
        return canvas

    monkeypatch.setattr(circular_assemble_module, "add_depth_group_on_canvas", fake_add_depth_group_on_canvas)
    monkeypatch.setattr(circular_assemble_module, "add_gc_content_group_on_canvas", fake_add_gc_content_group_on_canvas)
    monkeypatch.setattr(circular_assemble_module, "add_gc_skew_group_on_canvas", fake_add_gc_skew_group_on_canvas)

    assemble_circular_diagram_from_record(
        record,
        config_dict=config_dict,
        default_colors=default_colors,
        selected_features_set=SELECTED_FEATURES,
        legend="none",
        depth_table=depth_table,
        window=100,
        step=100,
        circular_track_slots=[
            *default_circular_track_slots(show_depth=True, show_gc=True, show_skew=True),
            "at_skew:dinucleotide_skew@nt=AT,w=20px",
        ],
    )

    definition_reserved = float(captured["definition_reserved"])
    for slot_id in ("depth", "gc_content", "gc_skew", "at_skew"):
        center_px, width_px = captured[slot_id]  # type: ignore[misc]
        assert center_px - (0.5 * width_px) >= definition_reserved - 1e-6


@pytest.mark.circular
def test_api_circular_track_slots_render_duplicate_dinucleotide_skew_slots() -> None:
    record = _load_record()
    config_dict = modify_config_dict(
        load_config_toml("gbdraw.data", "config.toml"),
        show_labels=False,
        show_gc=False,
        show_skew=False,
    )
    default_colors = load_default_colors("", palette="default")

    canvas = assemble_circular_diagram_from_record(
        record,
        config_dict=config_dict,
        default_colors=default_colors,
        selected_features_set=SELECTED_FEATURES,
        legend="right",
        circular_track_slots=[
            CircularTrackSlot(id="features", renderer="features"),
            "ticks:ticks",
            "gc_skew:dinucleotide_skew@nt=GC,w=20px",
            "at_skew:dinucleotide_skew@nt=AT,w=20px",
        ],
    )
    svg_text = canvas.tostring()

    assert 'id="gc_skew"' in svg_text
    assert 'id="at_skew"' in svg_text
    assert "GC skew" in svg_text
    assert "AT skew" in svg_text


@pytest.mark.circular
def test_slot_mode_draws_axis_when_ticks_are_disabled() -> None:
    record = _load_record()
    config_dict = modify_config_dict(
        load_config_toml("gbdraw.data", "config.toml"),
        show_labels=False,
        show_gc=False,
        show_skew=False,
    )
    default_colors = load_default_colors("", palette="default")

    canvas = assemble_circular_diagram_from_record(
        record,
        config_dict=config_dict,
        default_colors=default_colors,
        selected_features_set=SELECTED_FEATURES,
        legend="none",
        circular_track_slots=[CircularTrackSlot(id="ticks", renderer="ticks", enabled=False)],
    )

    assert 'id="Axis"' in canvas.tostring()


@pytest.mark.circular
def test_slot_mode_tick_radius_does_not_move_axis(monkeypatch: pytest.MonkeyPatch) -> None:
    import gbdraw.diagrams.circular.assemble as circular_assemble_module

    record = _load_record()
    config_dict = modify_config_dict(
        load_config_toml("gbdraw.data", "config.toml"),
        show_labels=False,
        show_gc=False,
        show_skew=False,
    )
    default_colors = load_default_colors("", palette="default")
    captured: dict[str, float | None] = {}

    def fake_add_tick_group_on_canvas(
        canvas,
        gb_record,
        canvas_config,
        config_dict,
        *,
        radius_override=None,
        tick_track_channel_override=None,
        label_side="legacy",
        tick_side="legacy",
        tick_length_px=None,
        cfg=None,
    ):
        captured["tick_radius"] = radius_override
        return canvas

    monkeypatch.setattr(circular_assemble_module, "add_tick_group_on_canvas", fake_add_tick_group_on_canvas)

    canvas = assemble_circular_diagram_from_record(
        record,
        config_dict=config_dict,
        default_colors=default_colors,
        selected_features_set=SELECTED_FEATURES,
        legend="none",
        circular_track_slots=["ticks:ticks@r=250px,w=12px,label_side=none,tick_side=inside"],
    )

    assert captured["tick_radius"] == pytest.approx(250.0)
    assert _axis_circle_radius(canvas.tostring()) == pytest.approx(390.0)


def test_cli_circular_track_order_forwards_slots(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    record = _load_record()
    captured: dict[str, object] = {}

    monkeypatch.setattr(circular_cli_module, "load_gbks", lambda paths, mode: [record])
    monkeypatch.setattr(circular_cli_module, "read_color_table", lambda _path: None)
    monkeypatch.setattr(circular_cli_module, "load_default_colors", lambda _path, _palette: None)
    monkeypatch.setattr(circular_cli_module, "save_figure", lambda canvas, formats: None)

    def fake_assemble(*args, **kwargs):
        captured["circular_track_slots"] = kwargs.get("circular_track_slots")
        return Drawing(filename=str(tmp_path / "dummy.svg"))

    monkeypatch.setattr(circular_cli_module, "assemble_circular_diagram_from_record", fake_assemble)

    circular_cli_module.circular_main(
        [
            "--gbk",
            "dummy.gb",
            "--circular_track_order",
            "features,ticks,gc_skew",
            "--nt",
            "AT",
            "--format",
            "svg",
            "-o",
            str(tmp_path / "out"),
        ]
    )

    slots = captured["circular_track_slots"]
    assert slots is not None
    assert [slot.id for slot in slots] == ["features", "ticks", "gc_skew"]
    assert slots[2].renderer == "dinucleotide_skew"
    assert slots[2].params["nt"] == "AT"


def test_cli_circular_track_slot_forwards_raw_specs(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    record = _load_record()
    captured: dict[str, object] = {}

    monkeypatch.setattr(circular_cli_module, "load_gbks", lambda paths, mode: [record])
    monkeypatch.setattr(circular_cli_module, "read_color_table", lambda _path: None)
    monkeypatch.setattr(circular_cli_module, "load_default_colors", lambda _path, _palette: None)
    monkeypatch.setattr(circular_cli_module, "save_figure", lambda canvas, formats: None)

    def fake_assemble(*args, **kwargs):
        captured["circular_track_slots"] = kwargs.get("circular_track_slots")
        return Drawing(filename=str(tmp_path / "dummy.svg"))

    monkeypatch.setattr(circular_cli_module, "assemble_circular_diagram_from_record", fake_assemble)

    circular_cli_module.circular_main(
        [
            "--gbk",
            "dummy.gb",
            "--circular_track_slot",
            "features:features",
            "--circular_track_slot",
            "at_skew:dinucleotide_skew@nt=AT,w=24px",
            "--format",
            "svg",
            "-o",
            str(tmp_path / "out"),
        ]
    )

    assert captured["circular_track_slots"] == [
        "features:features",
        "at_skew:dinucleotide_skew@nt=AT,w=24px",
    ]
