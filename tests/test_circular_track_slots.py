from __future__ import annotations

from pathlib import Path

import pytest
from Bio import SeqIO
from svgwrite import Drawing

import gbdraw.circular as circular_cli_module
from gbdraw.api.diagram import assemble_circular_diagram_from_record
from gbdraw.config.modify import modify_config_dict
from gbdraw.config.toml import load_config_toml
from gbdraw.io.colors import load_default_colors
from gbdraw.tracks import (
    CircularTrackLayoutContext,
    CircularTrackSlot,
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
            "ticks:ticks@axis=true",
            "gc_skew:dinucleotide_skew@nt=GC,w=20px",
            "at_skew:dinucleotide_skew@nt=AT,w=20px",
        ],
    )
    svg_text = canvas.tostring()

    assert 'id="gc_skew"' in svg_text
    assert 'id="at_skew"' in svg_text
    assert "GC skew" in svg_text
    assert "AT skew" in svg_text


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
