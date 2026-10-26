"""Rule L (PD-OI-078): the inside tracks around the Circular centre definition.

The preset lanes start right under the tick labels. When the definition's band reaches the inside tracks, they move
outward at their own widths, then narrow, each Auto width to ``canvas.circular.track_narrowing_min_px`` at the least;
when that cannot clear the band, they keep their widths and the definition is drawn over them. A typed width, a track
narrower than the floor and a pinned radius are kept. An explicit ``center_reserved_radius`` narrows as before
(DF-D05).
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest
from Bio import SeqIO

import gbdraw.diagrams.circular.assemble as circular_assemble_module
from gbdraw.api.diagram import assemble_circular_diagram_from_record
from gbdraw.canvas import CircularCanvasConfigurator
from gbdraw.cli import main as gbdraw_main
from gbdraw.config.models import CircularRenderProfile, GbdrawConfig
from gbdraw.config.modify import modify_config_dict
from gbdraw.config.toml import load_config_toml
from gbdraw.diagrams.circular.presets import (
    CircularPresetContext,
    circular_radial_plan_for_preset,
    circular_track_slots_from_preset_order,
)
from gbdraw.diagrams.circular.radial_layout import resolve_circular_radial_layout
from gbdraw.exceptions import ValidationError
from gbdraw.io.colors import load_default_colors
from gbdraw.layout.circular import CircularRadialLayout, CircularResolvedSlot
from gbdraw.tracks import parse_circular_track_slots

pytestmark = pytest.mark.circular

ROOT = Path(__file__).resolve().parent.parent
HMMT = ROOT / "tests" / "test_inputs" / "HmmtDNA.gbk"
MJENMV = ROOT / "examples" / "MjeNMV.gb"
SELECTED_FEATURES = ["CDS", "rRNA", "tRNA", "tmRNA", "ncRNA", "misc_RNA", "repeat_region"]
EPS = 1e-6


def _config(overrides: dict[str, object] | None = None, *, strandedness: bool = True) -> GbdrawConfig:
    return GbdrawConfig.from_dict(
        modify_config_dict(
            load_config_toml("gbdraw.data", "config.toml"),
            {
                "labels.circular.scope": "none",
                "canvas.show_gc": True,
                "canvas.show_skew": True,
                "canvas.circular.track_type": "tuckin",
                "canvas.strandedness": strandedness,
                **(overrides or {}),
            },
        )
    )


def _layout(
    band_px: float | None,
    *,
    overrides: dict[str, object] | None = None,
    strandedness: bool = True,
    slot_specs: list[str] | None = None,
    center_reserved_radius_explicit: bool = False,
) -> CircularRadialLayout:
    """The HmmtDNA tracks around a centre band of ``band_px`` (the default stack, or ``slot_specs`` on the preset)."""

    record = SeqIO.read(str(HMMT), "genbank")
    cfg = _config(overrides, strandedness=strandedness)
    canvas = CircularCanvasConfigurator("t", CircularRenderProfile(cfg), "none", record)
    context = CircularPresetContext(
        cfg=cfg,
        geometry=canvas.geometry(),
        total_length=len(record.seq),
        strandedness=strandedness,
        show_features=True,
        show_ticks=True,
        show_depth=False,
        show_gc=True,
        show_skew=True,
        track_ids=tuple(canvas.track_ids.items()),
    )
    plan = (
        circular_radial_plan_for_preset("tuckin", context)
        if slot_specs is None
        else circular_track_slots_from_preset_order(parse_circular_track_slots(slot_specs), "tuckin", context)
    )
    return resolve_circular_radial_layout(
        total_length=len(record.seq),
        geometry=canvas.geometry(),
        cfg=canvas.profile.config,
        slots=plan.slots,
        preferred_anchor_slot_ids=plan.preferred_anchor_slot_ids,
        definition_reserved_radius_px=band_px,
        center_reserved_radius_explicit=center_reserved_radius_explicit,
    )


def _by_id(layout: CircularRadialLayout) -> dict[str, CircularResolvedSlot]:
    return {slot.id: slot for slot in layout.slots}


def _inner_edge(layout: CircularRadialLayout) -> float:
    """The innermost reserved edge of the inside tracks: the band clears them when it ends at or inside this."""

    return min(
        float(slot.reserved_band_px.inner_px)
        for slot in layout.slots
        if slot.side == "inside" and slot.reserved_band_px is not None
    )


def _geometry(layout: CircularRadialLayout) -> dict[str, tuple[float, float]]:
    return {slot.id: (float(slot.anchor_radius_px), float(slot.resolved_width_px)) for slot in layout.slots}


GC_ROWS = ("gc_content", "gc_skew")
TYPED_WIDTH_SLOTS = [
    "features:features",
    "ticks:ticks",
    "gc_content:dinucleotide_content@side=inside,w=42px",
    "gc_skew:dinucleotide_skew@side=inside,w=34px",
]


def test_without_separate_strands_the_lanes_start_under_the_tick_labels() -> None:
    # OV-1850: with one feature lane the preset lanes kept a gap under the tick labels.
    gaps = {}
    for strandedness in (False, True):
        by_id = _by_id(_layout(None, strandedness=strandedness))
        gaps[strandedness] = float(by_id["ticks"].reserved_band_px.inner_px) - float(
            by_id["gc_content"].reserved_band_px.outer_px
        )
    assert gaps[False] == pytest.approx(gaps[True], abs=EPS)
    assert gaps[False] == pytest.approx(0.0, abs=EPS)


def test_the_band_moves_rows_outward_before_any_track_is_narrowed() -> None:
    # Typed widths narrower than their lanes leave room between the lanes: the band moves the rows into it.
    natural = _by_id(_layout(None, slot_specs=TYPED_WIDTH_SLOTS))
    band = float(natural["gc_skew"].reserved_band_px.inner_px) + 10.0
    moved = _layout(band, slot_specs=TYPED_WIDTH_SLOTS)
    by_id = _by_id(moved)
    assert by_id["gc_skew"].anchor_radius_px > natural["gc_skew"].anchor_radius_px + 9.0
    for slot_id in GC_ROWS:
        assert by_id[slot_id].resolved_width_px == pytest.approx(natural[slot_id].resolved_width_px)
        assert not by_id[slot_id].compressed, slot_id
    assert _inner_edge(moved) >= band - EPS
    assert moved.definition_reserved_band_px is not None
    assert moved.definition_reserved_band_px.outer_px == pytest.approx(band)


def test_auto_widths_narrow_to_clear_the_band_and_stay_at_or_above_the_floor() -> None:
    natural = _layout(None)
    band = 200.0
    assert _inner_edge(natural) < band
    narrowed = _layout(band)
    by_id = _by_id(narrowed)
    for slot_id in GC_ROWS:
        assert by_id[slot_id].compressed, slot_id
        assert 30.0 - EPS <= by_id[slot_id].resolved_width_px < _by_id(natural)[slot_id].resolved_width_px - 1.0
    assert _inner_edge(narrowed) >= band - EPS


def test_widths_stay_when_narrowing_to_the_floor_cannot_clear_the_band() -> None:
    natural = _layout(None)
    over = _layout(260.0)
    assert _geometry(over) == pytest.approx(_geometry(natural))
    assert not any(slot.compressed for slot in over.slots)
    assert _inner_edge(over) < 260.0
    assert over.definition_reserved_band_px is not None
    assert over.definition_reserved_band_px.outer_px == pytest.approx(260.0)


def test_turning_narrowing_off_keeps_the_widths() -> None:
    natural = _layout(None)
    kept = _layout(200.0, overrides={"canvas.circular.track_narrowing": False})
    assert _geometry(kept) == pytest.approx(_geometry(natural))
    assert not any(slot.compressed for slot in kept.slots)


@pytest.mark.parametrize(("floor_px", "narrows"), [(45.0, True), (60.0, False)])
def test_the_floor_comes_from_the_config(floor_px: float, narrows: bool) -> None:
    # At a 200 px band the GC rows clear it at about 52 px.
    natural = _by_id(_layout(None))
    by_id = _by_id(_layout(200.0, overrides={"canvas.circular.track_narrowing_min_px": floor_px}))
    for slot_id in GC_ROWS:
        assert by_id[slot_id].compressed is narrows, slot_id
        if narrows:
            assert by_id[slot_id].resolved_width_px >= floor_px - EPS
        else:
            assert by_id[slot_id].resolved_width_px == pytest.approx(natural[slot_id].resolved_width_px)


def test_a_track_narrower_than_the_floor_keeps_its_width() -> None:
    # A floor above the 74 px rows: a band that a few px of narrowing would clear leaves them as they are.
    natural = _layout(None)
    kept = _layout(170.0, overrides={"canvas.circular.track_narrowing_min_px": 100.0})
    assert _geometry(kept) == pytest.approx(_geometry(natural))
    assert not any(slot.compressed for slot in kept.slots)


@pytest.mark.parametrize(
    "gc_slots",
    [
        ["gc_content:dinucleotide_content@side=inside,w=74px", "gc_skew:dinucleotide_skew@side=inside,w=74px"],
        ["my_gc:dinucleotide_content@side=inside,w=74px", "my_skew:dinucleotide_skew@side=inside,w=74px"],
    ],
    ids=["preset_lane_ids", "custom_ids"],
)
def test_a_typed_width_is_never_narrowed(gc_slots: list[str]) -> None:
    # OV-1853: a width typed on a preset-lane slot was narrowed (to 57.3 px on MjeNMV); one on a custom slot was not.
    layout = _layout(200.0, slot_specs=["features:features", "ticks:ticks", *gc_slots])
    inside = [slot for slot in layout.slots if slot.renderer.startswith("dinucleotide")]
    assert [slot.resolved_width_px for slot in inside] == pytest.approx([74.0, 74.0])
    assert not any(slot.compressed for slot in inside)


@pytest.mark.parametrize(("band", "dev_width_px"), [(170.0, 68.45), (200.0, 51.8)])
def test_an_explicit_center_reserved_radius_narrows_typed_lane_widths_as_before(band: float, dev_width_px: float) -> None:
    # DF-D38: --gc_content_width 74 --gc_skew_width 74 with --center_reserved_radius 170 (or 200) drew on dev, the GC
    # rows narrowed to 68.45 (51.8) px; the typed-width exemption belongs to rule L's automatic fit only.
    slots = [
        "features:features",
        "ticks:ticks",
        "gc_content:dinucleotide_content@side=inside,w=74px",
        "gc_skew:dinucleotide_skew@side=inside,w=74px",
    ]
    layout = _layout(band, slot_specs=slots, center_reserved_radius_explicit=True)
    by_id = _by_id(layout)
    for slot_id in GC_ROWS:
        assert by_id[slot_id].compressed, slot_id
        assert by_id[slot_id].resolved_width_px == pytest.approx(dev_width_px), slot_id
    assert _inner_edge(layout) >= band - EPS


CUSTOM_STACK_SLOTS = [
    "features:features",
    "ticks:ticks",
    "gc_content:dinucleotide_content",
    "gc_skew:dinucleotide_skew",
    "at:dinucleotide_skew@nt=AT",
]


def test_turning_narrowing_off_keeps_the_auto_widths_of_a_custom_stack() -> None:
    # A custom stack's Auto widths narrow for the definition (all three rows at a 180 px band); off keeps them.
    natural = _layout(None, slot_specs=CUSTOM_STACK_SLOTS)
    narrowed = _by_id(_layout(180.0, slot_specs=CUSTOM_STACK_SLOTS))
    assert all(narrowed[slot_id].compressed for slot_id in ("gc_content", "gc_skew", "at"))
    kept = _layout(180.0, slot_specs=CUSTOM_STACK_SLOTS, overrides={"canvas.circular.track_narrowing": False})
    assert _geometry(kept) == pytest.approx(_geometry(natural))
    assert not any(slot.compressed for slot in kept.slots)


def test_an_explicit_center_reserved_radius_narrows_auto_widths_below_the_floor() -> None:
    # DF-D05: the user's own radius narrows the inside tracks to their readable minimums, below the floor; the
    # definition's band of the same size leaves their widths.
    explicit = _by_id(_layout(245.0, center_reserved_radius_explicit=True))
    definition = _by_id(_layout(245.0))
    for slot_id in GC_ROWS:
        assert explicit[slot_id].compressed, slot_id
        assert explicit[slot_id].resolved_width_px < 30.0
        assert not definition[slot_id].compressed, slot_id


def test_a_pinned_radius_stays_where_the_user_put_it() -> None:
    slots = [
        "features:features",
        "ticks:ticks",
        "gc_content:dinucleotide_content",
        "gc_skew:dinucleotide_skew@side=inside,r=0.45",
    ]
    natural = _by_id(_layout(None, slot_specs=slots))
    band = float(natural["gc_skew"].reserved_band_px.inner_px) + 20.0
    by_id = _by_id(_layout(band, slot_specs=slots))
    assert by_id["gc_skew"].anchor_radius_px == pytest.approx(natural["gc_skew"].anchor_radius_px)
    assert by_id["gc_skew"].explicit_anchor


def _drawn_fit(monkeypatch: pytest.MonkeyPatch, overrides: dict[str, object]) -> CircularRadialLayout:
    """The MjeNMV figure's own layout (the first fit; later ones probe smaller definition fonts)."""

    fits: list[CircularRadialLayout] = []
    fit = circular_assemble_module._resolve_radial_layout_fitting_definition

    def capture(*args, **kwargs):
        result = fit(*args, **kwargs)
        fits.append(result.layout)
        return result

    monkeypatch.setattr(circular_assemble_module, "_resolve_radial_layout_fitting_definition", capture)
    assemble_circular_diagram_from_record(
        SeqIO.read(str(MJENMV), "genbank"),
        cfg=_config(overrides),
        default_colors=load_default_colors("", palette="default"),
        selected_features_set=SELECTED_FEATURES,
        legend="none",
    )
    return fits[0]


@pytest.mark.parametrize(
    ("overrides", "narrows"),
    [
        ({}, True),
        ({"canvas.circular.track_narrowing_min_px": 50.0}, False),
        ({"canvas.circular.track_narrowing": False}, False),
    ],
    ids=["default_floor", "floor_50", "off"],
)
def test_the_mjenmv_name_narrows_the_gc_rows_only_down_to_the_floor(
    monkeypatch: pytest.MonkeyPatch, overrides: dict[str, object], narrows: bool
) -> None:
    # Separate strands (the Web default): the one-line name needs the GC rows at about 47 px. A 30 px floor allows it;
    # a 50 px floor or narrowing off keeps 74.1 px, and the name is drawn over the tracks (PD-OI-097).
    layout = _drawn_fit(monkeypatch, overrides)
    by_id = _by_id(layout)
    band = layout.definition_reserved_band_px
    assert band is not None
    for slot_id in GC_ROWS:
        assert by_id[slot_id].compressed is narrows, slot_id
        if narrows:
            assert 30.0 <= by_id[slot_id].resolved_width_px < 60.0
        else:
            assert by_id[slot_id].resolved_width_px == pytest.approx(74.1)
    assert (_inner_edge(layout) >= float(band.outer_px) - EPS) is narrows


def _cli_fit(monkeypatch: pytest.MonkeyPatch, args: list[str]) -> CircularRadialLayout:
    """The layout of the first definition fit of ``gbdraw circular`` run with ``args``."""

    fits: list[CircularRadialLayout] = []
    fit = circular_assemble_module._resolve_radial_layout_fitting_definition

    def capture(*fit_args, **kwargs):
        result = fit(*fit_args, **kwargs)
        fits.append(result.layout)
        return result

    monkeypatch.setattr(circular_assemble_module, "_resolve_radial_layout_fitting_definition", capture)
    monkeypatch.setattr("sys.argv", ["gbdraw", "circular", *args, "-f", "svg"])
    gbdraw_main()
    return fits[0]


@pytest.mark.parametrize(
    ("width_args", "narrows"),
    [([], True), (["--gc_content_width", "74", "--gc_skew_width", "74"], False)],
    ids=["auto_widths", "typed_shortcut_widths"],
)
def test_the_cli_width_shortcuts_are_typed_widths(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, width_args: list[str], narrows: bool
) -> None:
    # OV-1853: --gc_content_width and --gc_skew_width write the slots' widths, which the layout never narrows.
    by_id = _by_id(
        _cli_fit(monkeypatch, ["--gbk", str(MJENMV), "--separate_strands", *width_args, "-o", str(tmp_path / "mjenmv")])
    )
    for slot_id in GC_ROWS:
        assert by_id[slot_id].compressed is narrows, slot_id
        if not narrows:
            assert by_id[slot_id].resolved_width_px == pytest.approx(74.0)


@pytest.mark.parametrize(("width_args", "narrows"), [([], True), (["--depth_width", "37"], False)], ids=["auto", "typed"])
def test_a_typed_depth_width_is_kept_under_rule_l(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, width_args: list[str], narrows: bool
) -> None:
    # At a 16 pt definition the MjeNMV name narrows the inside rows down to the 30 px floor: the Auto Depth row
    # (37 px) narrows with them, a typed --depth_width stays.
    depth = tmp_path / "depth.tsv"
    depth.write_text("".join(f"LC738868.1\t{position}\t{10 + position % 7}\n" for position in range(1, 306008, 500)))
    layout = _cli_fit(
        monkeypatch,
        [
            "--gbk", str(MJENMV), "--separate_strands", "--depth_track", str(depth), "--definition_font_size", "16",
            *width_args, "-o", str(tmp_path / "depth"),
        ],
    )
    row = _by_id(layout)["depth"]
    assert row.compressed is narrows
    if not narrows:
        assert row.resolved_width_px == pytest.approx(37.0)
    assert _inner_edge(layout) >= float(layout.definition_reserved_band_px.outer_px) - EPS


def test_an_explicit_gc_content_radius_is_kept_under_rule_l(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    # The MjeNMV name's band reaches the GC content row placed at 0.6 of the radius; the row stays there.
    layout = _cli_fit(
        monkeypatch,
        ["--gbk", str(MJENMV), "--separate_strands", "--gc_content_radius", "0.6", "-o", str(tmp_path / "radius")],
    )
    row = _by_id(layout)["gc_content"]
    band = layout.definition_reserved_band_px
    assert band is not None and float(row.reserved_band_px.inner_px) < float(band.outer_px)
    assert row.anchor_radius_px == pytest.approx(0.6 * float(layout.axis.radius_px))


def test_the_cli_writes_the_narrowing_settings_into_the_config(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.chdir(tmp_path)
    session = tmp_path / "narrowing.gbdraw-session.json.gz"
    monkeypatch.setattr(
        "sys.argv",
        [
            "gbdraw", "circular", "--gbk", str(HMMT), "--no_track_narrowing", "--track_narrowing_min_px", "44",
            "-o", "narrowing", "--session_output", str(session),
        ],
    )
    gbdraw_main()
    circular = json.loads(gzip.decompress(session.read_bytes()))["renderRequest"]["diagramOptions"]["config"]["canvas"][
        "circular"
    ]
    assert circular["track_narrowing"] is False
    assert circular["track_narrowing_min_px"] == 44.0


def test_the_cli_refuses_a_negative_floor(capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr("sys.argv", ["gbdraw", "circular", "--gbk", str(HMMT), "--track_narrowing_min_px", "-1", "-o", "x"])
    with pytest.raises(SystemExit):
        gbdraw_main()
    assert "--track_narrowing_min_px must be a finite number >= 0" in capsys.readouterr().err


@pytest.mark.parametrize(("value", "reason"), [(-1, "NONNEGATIVE"), (float("nan"), "FINITE"), (True, "FINITE")])
def test_the_config_refuses_a_floor_that_is_not_a_nonnegative_number(value: object, reason: str) -> None:
    with pytest.raises(ValidationError) as caught:
        _config({"canvas.circular.track_narrowing_min_px": value})
    assert caught.value.diagnostic == {
        "code": "INPUT_INVALID",
        "reason": reason,
        "configPath": "canvas.circular.track_narrowing_min_px",
    }


@pytest.mark.parametrize(
    ("leaf", "value", "reason"),
    [("track_narrowing_min_px", -1, "NONNEGATIVE"), ("track_narrowing", "no", "BOOLEAN")],
)
def test_a_stored_full_config_with_a_bad_narrowing_leaf_names_it(leaf: str, value: object, reason: str) -> None:
    raw = load_config_toml("gbdraw.data", "config.toml")
    raw["canvas"]["circular"][leaf] = value
    with pytest.raises(ValidationError) as caught:
        GbdrawConfig.from_dict(raw)
    assert caught.value.diagnostic == {
        "code": "INPUT_INVALID",
        "reason": reason,
        "configPath": f"canvas.circular.{leaf}",
    }


def test_a_stored_full_config_without_the_narrowing_leaves_takes_the_defaults() -> None:
    # A Session written before rule L stores a full config without them.
    raw = load_config_toml("gbdraw.data", "config.toml")
    del raw["canvas"]["circular"]["track_narrowing"], raw["canvas"]["circular"]["track_narrowing_min_px"]
    circular = GbdrawConfig.from_dict(raw).canvas.circular
    assert (circular.track_narrowing, circular.track_narrowing_min_px) == (True, 30.0)
