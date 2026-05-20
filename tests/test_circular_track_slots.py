from __future__ import annotations

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
from gbdraw.svg.circular_ticks import get_circular_tick_path_ratio_bounds
from gbdraw.tracks import (
    CircularTrackSlot,
    CircularTrackSlotParseError,
    ScalarSpec,
    default_circular_track_slots,
    normalize_circular_track_slots,
    parse_circular_track_slot,
    parse_circular_track_slots,
)


SELECTED_FEATURES = ["CDS", "rRNA", "tRNA", "tmRNA", "ncRNA", "misc_RNA", "repeat_region"]


def _load_record():
    input_path = Path(__file__).parent / "test_inputs" / "HmmtDNA.gbk"
    return SeqIO.read(str(input_path), "genbank")


def _load_edl933_record():
    input_path = Path(__file__).parent / "test_inputs" / "EDL933.gbk"
    return SeqIO.read(str(input_path), "genbank")


def _load_mjenmv_record():
    input_path = Path(__file__).parent / "test_inputs" / "MjeNMV.gb"
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
    custom_slots: list[CircularTrackSlot] | None = None,
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
        track_preset=None,
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

    slots = custom_slots if custom_slots is not None else (
        default_circular_track_slots(show_depth=show_depth, show_gc=True, show_skew=True)
        if use_default_slots
        else None
    )
    kwargs = {}
    if use_default_slots or custom_slots is not None:
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


def _capture_circular_radial_layout(
    monkeypatch: pytest.MonkeyPatch,
    *,
    track_type: str,
    circular_track_slots: list[CircularTrackSlot] | None = None,
    input_filename: str = "HmmtDNA.gbk",
):
    import gbdraw.diagrams.circular.assemble as circular_assemble_module

    input_path = Path(__file__).parent / "test_inputs" / input_filename
    record = SeqIO.read(str(input_path), "genbank")
    config_dict = modify_config_dict(
        load_config_toml("gbdraw.data", "config.toml"),
        show_labels=False,
        show_gc=True,
        show_skew=True,
        track_type=track_type,
        strandedness=True,
    )
    default_colors = load_default_colors("", palette="default")
    captured: dict[str, object] = {}

    def capture_axis(canvas, canvas_config, *args, **kwargs):
        captured["radial_layout"] = canvas_config.circular_radial_layout
        return canvas

    monkeypatch.setattr(circular_assemble_module, "add_axis_group_on_canvas", capture_axis)
    kwargs = {}
    if circular_track_slots is not None:
        kwargs["circular_track_slots"] = circular_track_slots

    assemble_circular_diagram_from_record(
        record,
        config_dict=config_dict,
        default_colors=default_colors,
        selected_features_set=SELECTED_FEATURES,
        legend="none",
        **kwargs,
    )
    return captured["radial_layout"]


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
    assert slot.radius is not None
    assert slot.radius.resolve(390) == pytest.approx(163.8)
    assert slot.z == 7


def test_parse_circular_track_slots_rejects_duplicate_ids_and_unknown_renderer() -> None:
    with pytest.raises(CircularTrackSlotParseError, match="duplicate circular track slot id"):
        parse_circular_track_slots(["gc_skew:dinucleotide_skew", "gc_skew:dinucleotide_skew@nt=AT"])

    with pytest.raises(CircularTrackSlotParseError, match="unknown circular track renderer"):
        parse_circular_track_slots(["custom:not_a_renderer"])


def test_parse_circular_track_slots_normalizes_object_renderer_aliases() -> None:
    slots = parse_circular_track_slots([CircularTrackSlot(id="custom_skew", renderer="skew")])

    assert slots[0].renderer == "dinucleotide_skew"


def test_default_circular_track_slots_do_not_include_tick_axis_param() -> None:
    slots = default_circular_track_slots(show_depth=False, show_gc=True, show_skew=True)
    ticks = next(slot for slot in slots if slot.renderer == "ticks")

    assert "axis" not in ticks.params


def test_parse_circular_tick_slot_rejects_axis_param() -> None:
    with pytest.raises(CircularTrackSlotParseError, match="ticks slots no longer accept 'axis'"):
        parse_circular_track_slot("ticks:ticks@axis=false")

    with pytest.raises(CircularTrackSlotParseError, match="ticks slots no longer accept 'axis'"):
        parse_circular_track_slots(
            [CircularTrackSlot(id="ticks", renderer="ticks", params={"axis": True})]
        )


@pytest.mark.parametrize(
    "spec",
    ["axis", "axis@r=0.80", "gc_skew@w=5px"],
)
def test_circular_track_slot_shortcuts_are_rejected(spec: str) -> None:
    with pytest.raises(CircularTrackSlotParseError, match="require '<slot_id>:<renderer>@"):
        parse_circular_track_slot(spec)


@pytest.mark.parametrize(
    "spec",
    [
        "features:features@ri=0.75",
        "features:features@ro=0.82",
        "features:features@innerRadius=0.75",
        "features:features@outerRadius=0.82",
        "gc_skew:dinucleotide_skew@gap=4px",
        "gc_skew:dinucleotide_skew@gap_after=4px",
        "gc_skew:dinucleotide_skew@gapAfter=4px",
    ],
)
def test_circular_track_slot_rejects_obsolete_geometry_keys(spec: str) -> None:
    with pytest.raises(CircularTrackSlotParseError, match="no longer supported"):
        parse_circular_track_slot(spec)


@pytest.mark.parametrize("param_key", ["side", "radius", "width", "spacing", "gap_after", "inner_radius"])
def test_normalize_circular_track_slots_rejects_generic_layout_keys_in_params(param_key: str) -> None:
    with pytest.raises(ValueError, match="generic layout field"):
        normalize_circular_track_slots(
            [
                CircularTrackSlot(
                    id="gc_content",
                    renderer="dinucleotide_content",
                    params={param_key: "inside" if param_key == "side" else "1"},
                )
            ]
        )


def test_custom_slots_ignore_track_type_for_explicit_geometry(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    custom_slots = [
        CircularTrackSlot(
            id="features",
            renderer="features",
            radius=ScalarSpec(0.95, "factor"),
            width=ScalarSpec(36.0, "px"),
            params={"lane_direction": "inside"},
        ),
        CircularTrackSlot(
            id="ticks",
            renderer="ticks",
            side="outside",
            radius=ScalarSpec(1.05, "factor"),
            width=ScalarSpec(8.0, "px"),
            params={"tick_side": "outside", "label_side": "none"},
        ),
        CircularTrackSlot(
            id="gc_content",
            renderer="dinucleotide_content",
            radius=ScalarSpec(0.64, "factor"),
            width=ScalarSpec(30.0, "px"),
        ),
        CircularTrackSlot(
            id="gc_skew",
            renderer="dinucleotide_skew",
            radius=ScalarSpec(0.44, "factor"),
            width=ScalarSpec(24.0, "px"),
        ),
    ]
    tuckin = _capture_circular_core_geometry(
        monkeypatch,
        track_type="tuckin",
        custom_slots=custom_slots,
    )
    monkeypatch.undo()
    spreadout = _capture_circular_core_geometry(
        monkeypatch,
        track_type="spreadout",
        custom_slots=custom_slots,
    )

    _assert_geometry_matches(spreadout, tuckin)


@pytest.mark.parametrize("track_type", ["tuckin", "middle", "spreadout"])
def test_default_custom_slots_auto_place_near_preset_without_radius_inheritance(
    monkeypatch: pytest.MonkeyPatch,
    track_type: str,
) -> None:
    preset_layout = _capture_circular_radial_layout(
        monkeypatch,
        track_type=track_type,
        input_filename="MG1655.gbk",
    )
    monkeypatch.undo()
    custom_layout = _capture_circular_radial_layout(
        monkeypatch,
        track_type=track_type,
        circular_track_slots=default_circular_track_slots(show_depth=False, show_gc=True, show_skew=True),
        input_filename="MG1655.gbk",
    )

    preset_by_id = {slot.id: slot for slot in preset_layout.slots}  # type: ignore[attr-defined]
    custom_by_id = {slot.id: slot for slot in custom_layout.slots}  # type: ignore[attr-defined]
    default_spacing_px = max(1.0, 0.01 * float(custom_layout.axis.radius_px))  # type: ignore[attr-defined]

    assert custom_by_id.keys() == preset_by_id.keys()
    for slot_id, custom_slot in custom_by_id.items():
        assert not custom_slot.explicit_anchor
        assert custom_slot.anchor_radius_px == pytest.approx(
            preset_by_id[slot_id].anchor_radius_px,
            abs=2.0 * default_spacing_px,
        )


def test_custom_slot_order_places_ticks_outside_feature_band_when_ordered_before_features(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import gbdraw.diagrams.circular.assemble as circular_assemble_module

    record = _load_edl933_record()
    config_dict = modify_config_dict(
        load_config_toml("gbdraw.data", "config.toml"),
        show_labels=False,
        show_gc=False,
        show_skew=False,
        track_type="tuckin",
        strandedness=True,
    )
    default_colors = load_default_colors("", palette="default")
    captured: dict[str, object] = {}

    def fake_add_record_group_on_canvas(canvas, *args, **kwargs):
        canvas_config = args[1]
        captured["radial_layout"] = canvas_config.circular_radial_layout
        return canvas

    def fake_add_tick_group_on_canvas(canvas, *args, **kwargs):
        canvas_config = args[1]
        captured["radial_layout"] = canvas_config.circular_radial_layout
        return canvas

    monkeypatch.setattr(circular_assemble_module, "add_record_group_on_canvas", fake_add_record_group_on_canvas)
    monkeypatch.setattr(circular_assemble_module, "add_tick_group_on_canvas", fake_add_tick_group_on_canvas)

    assemble_circular_diagram_from_record(
        record,
        config_dict=config_dict,
        default_colors=default_colors,
        selected_features_set=SELECTED_FEATURES,
        legend="none",
        circular_track_slots=[
            CircularTrackSlot(id="ticks", renderer="ticks"),
            CircularTrackSlot(id="features", renderer="features"),
        ],
    )

    layout = captured["radial_layout"]
    assert layout.features is not None
    assert layout.ticks is not None
    assert layout.ticks.anchor_radius_px > layout.features.anchor_radius_px
    assert layout.ticks.tick_band_px.inner_px >= layout.features.all_band_px.outer_px
    assert layout.ticks.label_band_px is not None
    assert layout.ticks.label_band_px.inner_px >= layout.features.all_band_px.outer_px


def test_parse_circular_track_slot_stores_layout_fields_on_slot() -> None:
    slot = parse_circular_track_slot(
        "at_skew:dinucleotide_skew@nt=AT,w=24px,spacing=4px,side=inside,strict=true,compress=true,reserve=true,z=7"
    )

    assert slot.side == "inside"
    assert slot.strict is True
    assert slot.compress is True
    assert slot.reserve is True
    assert slot.spacing is not None
    assert slot.spacing.resolve(390.0) == pytest.approx(4.0)
    assert "side" not in slot.params
    assert "strict" not in slot.params
    assert "compress" not in slot.params
    assert "reserve" not in slot.params

    normalized = normalize_circular_track_slots([slot])[0]
    assert normalized.side == "inside"
    assert normalized.strict is True
    assert normalized.compress is True
    assert normalized.reserve is True
    assert normalized.params["nt"] == "AT"


def test_normalize_circular_track_slots_derives_feature_side_from_lane_direction() -> None:
    inside, split, outside = normalize_circular_track_slots(
        [
            CircularTrackSlot(id="features_inside", renderer="features", params={"lane_direction": "inside"}),
            CircularTrackSlot(id="features_split", renderer="features", params={"lane_direction": "split"}),
            CircularTrackSlot(id="features_outside", renderer="features", params={"lane_direction": "outside"}),
        ]
    )

    assert inside.side == "inside"
    assert split.side == "overlay"
    assert split.reserve is True
    assert outside.side == "outside"


@pytest.mark.circular
def test_circular_preset_slots_do_not_emit_origin_metadata() -> None:
    from gbdraw.canvas import CircularCanvasConfigurator
    from gbdraw.config.models import GbdrawConfig
    from gbdraw.diagrams.circular.presets import (
        CircularPresetContext,
        circular_radial_plan_for_preset,
        circular_track_slots_for_preset,
    )

    record = _load_record()
    config_dict = _base_config(track_type="tuckin")
    cfg = GbdrawConfig.from_dict(config_dict)
    canvas_config = CircularCanvasConfigurator("test", config_dict, "none", record, cfg=cfg)
    context = CircularPresetContext(
        cfg=cfg,
        canvas_config=canvas_config,
        total_length=len(record.seq),
        strandedness=bool(cfg.canvas.strandedness),
        show_features=True,
        show_ticks=True,
        show_depth=False,
        show_gc=True,
        show_skew=True,
    )

    slots = circular_track_slots_for_preset("tuckin", context)
    assert all("_preset_generated" not in dict(slot.params) for slot in slots)

    plan = circular_radial_plan_for_preset("tuckin", context)
    assert {"gc_content", "gc_skew"} <= set(plan.preferred_anchor_slot_ids)


@pytest.mark.circular
def test_blank_builtin_numeric_slots_use_preset_preferred_anchor_without_explicit_radius() -> None:
    from gbdraw.canvas import CircularCanvasConfigurator
    from gbdraw.config.models import GbdrawConfig
    from gbdraw.diagrams.circular.presets import CircularPresetContext, circular_track_slots_from_preset_order
    from gbdraw.diagrams.circular.radial_layout import resolve_circular_radial_layout

    record = _load_record()
    config_dict = _base_config(track_type="middle")
    cfg = GbdrawConfig.from_dict(config_dict)
    canvas_config = CircularCanvasConfigurator("test", config_dict, "none", record, cfg=cfg)
    context = CircularPresetContext(
        cfg=cfg,
        canvas_config=canvas_config,
        total_length=len(record.seq),
        strandedness=bool(cfg.canvas.strandedness),
        show_features=True,
        show_ticks=True,
        show_depth=False,
        show_gc=True,
        show_skew=True,
    )

    plan = circular_track_slots_from_preset_order(
        default_circular_track_slots(show_depth=False, show_gc=True, show_skew=True),
        "middle",
        context,
    )
    plan_by_id = {slot.id: slot for slot in plan.slots}

    assert plan_by_id["gc_content"].radius is None
    assert plan_by_id["gc_skew"].radius is None
    assert "_preferred_anchor_radius" in plan_by_id["gc_content"].params
    assert "_preferred_anchor_radius" in plan_by_id["gc_skew"].params

    layout = resolve_circular_radial_layout(
        total_length=len(record.seq),
        canvas_config=canvas_config,
        cfg=cfg,
        slots=plan.slots,
        preferred_anchor_slot_ids=plan.preferred_anchor_slot_ids,
    )
    by_id = {slot.id: slot for slot in layout.slots}
    length_param = str(canvas_config.length_param)
    track_dict = cfg.canvas.circular.track_dict[length_param]["middle"]

    assert not by_id["gc_content"].explicit_anchor
    assert not by_id["gc_skew"].explicit_anchor
    assert by_id["gc_content"].anchor_radius_px == pytest.approx(
        float(canvas_config.radius) * float(track_dict["2"])
    )
    assert by_id["gc_skew"].anchor_radius_px == pytest.approx(
        float(canvas_config.radius) * float(track_dict["3"])
    )


@pytest.mark.circular
def test_default_preset_tick_draw_uses_resolved_tick_options(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import gbdraw.diagrams.circular.assemble as circular_assemble_module

    record = _load_record()
    config_dict = modify_config_dict(
        load_config_toml("gbdraw.data", "config.toml"),
        show_labels=False,
        show_gc=False,
        show_skew=False,
        track_type="spreadout",
        strandedness=True,
    )
    default_colors = load_default_colors("", palette="default")
    captured: dict[str, object] = {}

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
        track_preset=None,
        cfg=None,
    ):
        captured["label_side"] = label_side
        captured["tick_side"] = tick_side
        captured["tick_length_px"] = tick_length_px
        captured["track_preset"] = track_preset
        return canvas

    monkeypatch.setattr(circular_assemble_module, "add_tick_group_on_canvas", fake_add_tick_group_on_canvas)

    assemble_circular_diagram_from_record(
        record,
        config_dict=config_dict,
        default_colors=default_colors,
        selected_features_set=SELECTED_FEATURES,
        legend="none",
    )

    assert captured["label_side"] != "legacy"
    assert captured["tick_side"] != "legacy"
    assert captured["tick_length_px"] is not None
    assert captured["track_preset"] == "spreadout"


@pytest.mark.circular
def test_default_preset_slots_compress_to_clear_center_definition(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import gbdraw.diagrams.circular.assemble as circular_assemble_module
    from gbdraw.config.models import GbdrawConfig

    record = _load_mjenmv_record()
    config_dict = _base_config(track_type="tuckin")
    cfg = GbdrawConfig.from_dict(config_dict)
    default_colors = load_default_colors("", palette="default")
    captured: dict[str, object] = {}

    def capture_numeric_slot(
        canvas_config,
    ) -> None:
        captured["radial_layout"] = canvas_config.circular_radial_layout
        captured["definition_reserved"] = circular_assemble_module._definition_reserved_radius_px(
            record,
            canvas_config,
            None,
            None,
            config_dict,
            cfg=cfg,
            plot_title=None,
            definition_profile="full",
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
        capture_numeric_slot(canvas_config)
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
        capture_numeric_slot(canvas_config)
        return canvas

    monkeypatch.setattr(circular_assemble_module, "add_gc_content_group_on_canvas", fake_add_gc_content_group_on_canvas)
    monkeypatch.setattr(circular_assemble_module, "add_gc_skew_group_on_canvas", fake_add_gc_skew_group_on_canvas)

    assemble_circular_diagram_from_record(
        record,
        config_dict=config_dict,
        default_colors=default_colors,
        selected_features_set=SELECTED_FEATURES,
        legend="none",
    )

    layout = captured["radial_layout"]
    by_id = {track.id: track for track in layout.tracks}  # type: ignore[attr-defined]
    definition_reserved = float(captured["definition_reserved"])
    assert by_id["gc_content"].reserved_inner_radius_px >= definition_reserved - 1e-6
    assert by_id["gc_skew"].reserved_inner_radius_px >= definition_reserved - 1e-6
    assert by_id["gc_content"].compressed or by_id["gc_skew"].compressed


@pytest.mark.circular
def test_default_custom_slots_tuckin_inherit_preset_when_inside_numeric_tracks_are_tight(
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

    assert {"gc_content", "gc_skew"} <= set(captured)
    assert captured["gc_content"][0] > captured["gc_skew"][0]


@pytest.mark.circular
def test_explicit_pure_auto_slots_tuckin_raise_when_inside_numeric_tracks_cannot_fit(
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

    monkeypatch.setattr(circular_assemble_module, "add_gc_content_group_on_canvas", lambda *args, **kwargs: args[0])
    monkeypatch.setattr(circular_assemble_module, "add_gc_skew_group_on_canvas", lambda *args, **kwargs: args[0])

    with pytest.raises(Exception, match="cannot fit inside"):
        assemble_circular_diagram_from_record(
            record,
            config_dict=config_dict,
            default_colors=default_colors,
            selected_features_set=SELECTED_FEATURES,
            legend="none",
            circular_track_slots=[
                CircularTrackSlot(id="features", renderer="features"),
                CircularTrackSlot(id="ticks", renderer="ticks"),
                CircularTrackSlot(
                    id="gc_content",
                    renderer="dinucleotide_content",
                    side="inside",
                    width=ScalarSpec(240.0, "px"),
                    compress=True,
                    params={"nt": "GC"},
                ),
                CircularTrackSlot(
                    id="gc_skew",
                    renderer="dinucleotide_skew",
                    side="inside",
                    width=ScalarSpec(240.0, "px"),
                    compress=True,
                    params={"nt": "GC"},
                ),
            ],
        )


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
                params={"label_side": "outside", "tick_side": "inside"},
            ),
            CircularTrackSlot(id="gc_skew", renderer="dinucleotide_skew"),
            CircularTrackSlot(id="gc_content", renderer="dinucleotide_content"),
        ],
    )

    assert captured["gc_skew"][0] > captured["gc_content"][0]


@pytest.mark.circular
def test_edl933_ticks_before_features_use_measured_tick_footprint(
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
    captured: dict[str, object] = {}

    def fake_add_record_group_on_canvas(canvas, *args, **kwargs):
        canvas_config = args[1]
        captured["radial_layout"] = canvas_config.circular_radial_layout
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

    def fake_add_tick_group_on_canvas(canvas, *args, **kwargs):
        canvas_config = args[1]
        captured["radial_layout"] = canvas_config.circular_radial_layout
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
            CircularTrackSlot(
                id="ticks",
                renderer="ticks",
                params={"label_side": "inside", "tick_side": "inside"},
            ),
            CircularTrackSlot(id="features", renderer="features"),
            CircularTrackSlot(id="gc_content", renderer="dinucleotide_content", params={"nt": "GC"}),
            CircularTrackSlot(id="gc_skew", renderer="dinucleotide_skew", params={"nt": "GC"}),
        ],
    )

    gc_center, gc_width = captured["gc_content"]  # type: ignore[misc]
    layout = captured["radial_layout"]

    assert 0.0 < gc_width <= float(captured["default_gc_width"]) + 1e-6
    assert 0.0 < captured["gc_skew"][1] <= float(captured["default_gc_width"]) + 1e-6  # type: ignore[index]
    assert captured["gc_content"][0] > captured["gc_skew"][0]  # type: ignore[index]
    assert layout.features is not None
    assert layout.ticks is not None
    assert layout.ticks.reserved_band_px.inner_px >= layout.features.all_band_px.outer_px
    assert (float(gc_center) - (0.5 * float(gc_width))) > 0.0


@pytest.mark.circular
def test_inside_order_reserves_stranded_feature_stack_between_numeric_slots(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import gbdraw.diagrams.circular.assemble as circular_assemble_module

    record = _load_record()
    config_dict = modify_config_dict(
        load_config_toml("gbdraw.data", "config.toml"),
        show_labels=False,
        show_gc=True,
        show_skew=True,
        track_type="tuckin",
        strandedness=True,
    )
    default_colors = load_default_colors("", palette="default")
    captured: dict[str, object] = {}

    def capture_record_layout(canvas, gb_record, canvas_config, *args, **kwargs):
        captured["radial_layout"] = canvas_config.circular_radial_layout
        return canvas

    monkeypatch.setattr(circular_assemble_module, "add_record_group_on_canvas", capture_record_layout)
    monkeypatch.setattr(circular_assemble_module, "add_tick_group_on_canvas", lambda *args, **kwargs: args[0])
    monkeypatch.setattr(circular_assemble_module, "add_gc_content_group_on_canvas", lambda *args, **kwargs: args[0])
    monkeypatch.setattr(circular_assemble_module, "add_gc_skew_group_on_canvas", lambda *args, **kwargs: args[0])

    assemble_circular_diagram_from_record(
        record,
        config_dict=config_dict,
        default_colors=default_colors,
        selected_features_set=SELECTED_FEATURES,
        legend="none",
        circular_track_slots=[
            CircularTrackSlot(id="ticks", renderer="ticks"),
            CircularTrackSlot(id="gc_content", renderer="dinucleotide_content", params={"nt": "GC"}),
            CircularTrackSlot(id="features", renderer="features"),
            CircularTrackSlot(id="gc_skew", renderer="dinucleotide_skew", params={"nt": "GC"}),
        ],
    )

    layout = captured["radial_layout"]
    by_id = {slot.id: slot for slot in layout.slots}  # type: ignore[attr-defined]

    assert by_id["features"].reserved_width_px > by_id["features"].resolved_width_px
    assert by_id["gc_content"].packing_band_px.center_px > by_id["ticks"].packing_band_px.center_px
    assert by_id["ticks"].packing_band_px.center_px > by_id["features"].packing_band_px.center_px
    assert by_id["features"].packing_band_px.center_px > by_id["gc_skew"].packing_band_px.center_px


@pytest.mark.circular
def test_order_only_gc_content_falls_back_later_auto_tracks_outside_when_needed(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import gbdraw.diagrams.circular.assemble as circular_assemble_module

    record = _load_mjenmv_record()
    config_dict = modify_config_dict(
        load_config_toml("gbdraw.data", "config.toml"),
        show_labels=False,
        show_gc=True,
        show_skew=True,
        track_type="tuckin",
        strandedness=True,
    )
    default_colors = load_default_colors("", palette="default")
    captured: dict[str, object] = {}

    def capture_record_layout(canvas, gb_record, canvas_config, *args, **kwargs):
        captured["radial_layout"] = canvas_config.circular_radial_layout
        return canvas

    monkeypatch.setattr(circular_assemble_module, "add_record_group_on_canvas", capture_record_layout)
    monkeypatch.setattr(circular_assemble_module, "add_tick_group_on_canvas", lambda *args, **kwargs: args[0])
    monkeypatch.setattr(circular_assemble_module, "add_gc_content_group_on_canvas", lambda *args, **kwargs: args[0])
    monkeypatch.setattr(circular_assemble_module, "add_gc_skew_group_on_canvas", lambda *args, **kwargs: args[0])

    assemble_circular_diagram_from_record(
        record,
        config_dict=config_dict,
        default_colors=default_colors,
        selected_features_set=SELECTED_FEATURES,
        legend="none",
        circular_track_slots=[
            CircularTrackSlot(id="gc_content", renderer="dinucleotide_content", params={"nt": "GC"}),
            CircularTrackSlot(id="features", renderer="features"),
            CircularTrackSlot(id="gc_skew", renderer="dinucleotide_skew", params={"nt": "GC"}),
            CircularTrackSlot(id="ticks", renderer="ticks"),
        ],
    )

    layout = captured["radial_layout"]
    by_id = {slot.id: slot for slot in layout.slots}  # type: ignore[attr-defined]
    assert by_id["gc_content"].packing_band_px.center_px > by_id["features"].packing_band_px.center_px
    assert by_id["features"].packing_band_px.center_px > by_id["gc_skew"].packing_band_px.center_px
    assert by_id["gc_skew"].packing_band_px.center_px > by_id["ticks"].packing_band_px.center_px


@pytest.mark.circular
def test_order_only_numeric_before_ticks_reserves_inner_numeric_space(
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
    captured: dict[str, object] = {}

    def capture_layout(canvas, *args, **kwargs):
        canvas_config = args[1]
        captured["radial_layout"] = canvas_config.circular_radial_layout
        return canvas

    def passthrough(canvas, *args, **kwargs):
        return canvas

    monkeypatch.setattr(circular_assemble_module, "add_record_group_on_canvas", capture_layout)
    monkeypatch.setattr(circular_assemble_module, "add_tick_group_on_canvas", capture_layout)
    monkeypatch.setattr(circular_assemble_module, "add_gc_content_group_on_canvas", passthrough)
    monkeypatch.setattr(circular_assemble_module, "add_gc_skew_group_on_canvas", passthrough)

    assemble_circular_diagram_from_record(
        record,
        config_dict=config_dict,
        default_colors=default_colors,
        selected_features_set=SELECTED_FEATURES,
        legend="none",
        circular_track_slots=[
            CircularTrackSlot(id="features", renderer="features"),
            CircularTrackSlot(id="gc_content", renderer="dinucleotide_content", params={"nt": "GC"}),
            CircularTrackSlot(id="ticks", renderer="ticks"),
            CircularTrackSlot(id="gc_skew", renderer="dinucleotide_skew", params={"nt": "GC"}),
        ],
    )

    layout = captured["radial_layout"]
    by_id = {slot.id: slot for slot in layout.slots}  # type: ignore[attr-defined]

    assert not by_id["ticks"].explicit_anchor
    assert by_id["features"].packing_band_px.center_px > by_id["gc_content"].packing_band_px.center_px
    assert by_id["gc_content"].packing_band_px.center_px > by_id["ticks"].packing_band_px.center_px
    assert by_id["ticks"].packing_band_px.center_px > by_id["gc_skew"].packing_band_px.center_px


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
    assert captured["gc_content"][1] > 0
    assert captured["gc_skew"][1] > 0


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
            "at_skew:dinucleotide_skew@nt=AT,w=20px,side=outside",
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
def test_api_circular_track_slots_distribute_extra_dinucleotide_slots_evenly(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import gbdraw.diagrams.circular.assemble as circular_assemble_module

    record = _load_edl933_record()
    config_dict = modify_config_dict(
        load_config_toml("gbdraw.data", "config.toml"),
        show_labels=False,
        show_gc=True,
        show_skew=True,
        track_type="middle",
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
        window=1000,
        step=1000,
        circular_track_slots=[
            CircularTrackSlot(id="features", renderer="features"),
            CircularTrackSlot(
                id="ticks",
                renderer="ticks",
                params={"label_side": "outside", "tick_side": "inside"},
            ),
            "gc_content:dinucleotide_content@nt=GC,side=outside",
            "gc_skew:dinucleotide_skew@nt=GC,side=outside",
            "gc_skew_2:dinucleotide_skew@nt=AT,side=outside",
        ],
    )

    assert set(captured) == {"gc_content", "gc_skew", "gc_skew_2"}
    annuli = {
        slot_id: (center_px - (0.5 * width_px), center_px + (0.5 * width_px))
        for slot_id, (center_px, width_px) in captured.items()
    }
    widths = [captured[slot_id][1] for slot_id in ("gc_content", "gc_skew", "gc_skew_2")]
    ordered_annuli = sorted(annuli.values(), key=lambda item: item[0])

    assert min(widths) > 30.0
    for previous, current in zip(ordered_annuli, ordered_annuli[1:]):
        assert current[0] >= previous[1] - 1e-6


@pytest.mark.circular
def test_api_explicit_inside_duplicate_dinucleotide_skew_places_when_space_is_reserved(
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
    captured: dict[str, object] = {}

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
        captured["radial_layout"] = canvas_config.circular_radial_layout

    def fake_add_record_group_on_canvas(canvas, *args, **kwargs):
        return canvas

    def fake_add_tick_group_on_canvas(canvas, *args, **kwargs):
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

    monkeypatch.setattr(circular_assemble_module, "add_record_group_on_canvas", fake_add_record_group_on_canvas)
    monkeypatch.setattr(circular_assemble_module, "add_tick_group_on_canvas", fake_add_tick_group_on_canvas)
    monkeypatch.setattr(circular_assemble_module, "add_gc_content_group_on_canvas", fake_add_gc_content_group_on_canvas)
    monkeypatch.setattr(circular_assemble_module, "add_gc_skew_group_on_canvas", fake_add_gc_skew_group_on_canvas)

    assemble_circular_diagram_from_record(
        record,
        config_dict=config_dict,
        default_colors=default_colors,
        selected_features_set=SELECTED_FEATURES,
        legend="none",
        window=1000,
        step=1000,
        circular_track_slots=[
            *default_circular_track_slots(show_depth=False, show_gc=True, show_skew=True),
            "gc_skew_2:dinucleotide_skew@nt=AT,side=inside,compress=true,strict=false",
        ],
    )

    assert {"gc_content", "gc_skew", "gc_skew_2"} <= set(captured)
    layout = captured["radial_layout"]
    by_id = {slot.id: slot for slot in layout.slots}  # type: ignore[attr-defined]
    assert by_id["gc_skew_2"].packing_band_px is not None


@pytest.mark.circular
def test_api_circular_track_slots_distribute_repeated_depth_slots_evenly(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import gbdraw.diagrams.circular.assemble as circular_assemble_module

    record = _load_edl933_record()
    config_dict = modify_config_dict(
        load_config_toml("gbdraw.data", "config.toml"),
        show_labels=False,
        show_gc=False,
        show_skew=False,
        track_type="middle",
        strandedness=True,
    )
    default_colors = load_default_colors("", palette="default")
    depth_table = _depth_table(str(record.id))
    captured: dict[str, tuple[float, float]] = {}

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
        assert track_width_override is not None
        assert norm_factor_override is not None
        captured[str(group_id or "depth")] = (
            float(norm_factor_override) * float(canvas_config.radius),
            float(track_width_override),
        )
        return canvas

    monkeypatch.setattr(circular_assemble_module, "add_depth_group_on_canvas", fake_add_depth_group_on_canvas)

    assemble_circular_diagram_from_record(
        record,
        config_dict=config_dict,
        default_colors=default_colors,
        selected_features_set=SELECTED_FEATURES,
        legend="none",
        depth_table=depth_table,
        window=1000,
        step=1000,
        circular_track_slots=[
            "depth:depth",
            "depth_2:depth",
            "depth_3:depth",
        ],
    )

    assert set(captured) == {"depth", "depth_2", "depth_3"}
    annuli = {
        slot_id: (center_px - (0.5 * width_px), center_px + (0.5 * width_px))
        for slot_id, (center_px, width_px) in captured.items()
    }
    widths = [captured[slot_id][1] for slot_id in ("depth", "depth_2", "depth_3")]
    gaps = [
        annuli["depth"][0] - annuli["depth_2"][1],
        annuli["depth_2"][0] - annuli["depth_3"][1],
    ]

    assert min(widths) > 20.0
    assert widths[0] == pytest.approx(widths[1])
    assert widths[1] == pytest.approx(widths[2])
    assert gaps[0] == pytest.approx(gaps[1])


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
        track_preset=None,
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
