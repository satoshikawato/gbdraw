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
    assert by_id["gc_content"].width_px == 74.1
    assert by_id["gc_skew"].center_radius_px == 171.6
    assert by_id["gc_skew"].width_px == 74.1


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

    assert by_id["gc_content"].width_px == pytest.approx(by_id["gc_skew"].width_px)
    assert by_id["gc_content"].width_px < 30.0
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
    assert by_id["gc_content"].width_px == pytest.approx(by_id["gc_skew"].width_px)
    assert by_id["gc_content"].width_px < 30.0
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
