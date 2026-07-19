from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from svgwrite import Drawing

import gbdraw.linear as linear_cli_module
from gbdraw.api import assemble_linear_diagram_from_records
from gbdraw.canvas import LinearCanvasConfigurator
from gbdraw.config.models import GbdrawConfig
from gbdraw.config.toml import load_config_toml
from gbdraw.exceptions import ValidationError
from gbdraw.io.colors import load_default_colors
from gbdraw.diagrams.linear.track_slots import LinearTrackLayout, resolve_linear_track_layout
from gbdraw.tracks import (
    LinearTrackSlot,
    LinearTrackSlotParseError,
    default_linear_track_slots,
    linear_track_slots_from_order,
    normalize_linear_track_slots,
    normalize_linear_track_slots_with_axis,
    parse_linear_track_slot,
    parse_linear_track_slots,
)


def _record(record_id: str = "rec1") -> SeqRecord:
    record = SeqRecord(Seq("ATGCGC" * 200), id=record_id, description=record_id)
    record.features = [
        SeqFeature(FeatureLocation(10, 120, strand=1), type="CDS"),
        SeqFeature(FeatureLocation(180, 260, strand=-1), type="tRNA"),
    ]
    return record


def _depth_table(reference: str = "rec1") -> pd.DataFrame:
    return pd.DataFrame(
        {
            "reference_name": [reference, reference, reference, reference],
            "position": [1, 2, 10, 30],
            "depth": [10, 20, 40, 80],
        }
    )


def _resolve_layout(slot_specs: list[str]) -> tuple[LinearTrackLayout, LinearCanvasConfigurator]:
    config_dict = load_config_toml("gbdraw.data", "config.toml")
    cfg = GbdrawConfig.from_dict(config_dict)
    canvas_config = LinearCanvasConfigurator(
        num_of_entries=1,
        longest_genome=1200,
        config_dict=config_dict,
        legend="none",
        cfg=cfg,
    )
    slots = normalize_linear_track_slots_with_axis(parse_linear_track_slots(slot_specs), None)
    return (
        resolve_linear_track_layout(
            slots,
            canvas_config=canvas_config,
            cfg=cfg,
        ),
        canvas_config,
    )


def _layout_track(layout: LinearTrackLayout, track_id: str):
    track = next((slot for slot in layout.slots if slot.id == track_id), None)
    assert track is not None
    return track


def _extract_group_y(svg_text: str, group_id: str) -> float:
    match = re.search(
        rf'<g id="{re.escape(group_id)}" transform="translate\([^,]+,([\-0-9.]+)\)">',
        svg_text,
    )
    assert match is not None
    return float(match.group(1))


def _extract_axis_baseline_y(svg_text: str, group_id: str) -> float:
    group_match = re.search(
        rf'<g id="{re.escape(group_id)}">(.*?)</g>',
        svg_text,
        flags=re.DOTALL,
    )
    assert group_match is not None
    for line_match in re.finditer(r"<line\b[^>]*>", group_match.group(1)):
        attrs = dict(re.findall(r'([A-Za-z_:][\w:.-]*)="([^"]*)"', line_match.group(0)))
        if "y1" in attrs and attrs.get("y2") == attrs["y1"]:
            return float(attrs["y1"])
    raise AssertionError(f"axis baseline line not found for {group_id}")


def _extract_group_fragment(svg_text: str, group_id: str) -> str:
    match = re.search(rf'<g id="{re.escape(group_id)}"[^>]*>.*?</g>', svg_text, flags=re.DOTALL)
    assert match is not None
    return match.group(0)


def test_parse_linear_track_slot_with_layout_fields() -> None:
    slot = parse_linear_track_slot("depth_1:depth@track_index=0,h=35px,spacing=8,z=2,side=below")

    assert slot.id == "depth_1"
    assert slot.renderer == "depth"
    assert slot.side == "below"
    assert slot.height is not None
    assert slot.height.resolve(1.0) == pytest.approx(35.0)
    assert slot.spacing is not None
    assert slot.spacing.resolve(1.0) == pytest.approx(8.0)
    assert slot.z == 2
    assert slot.params["track_index"] == 0


@pytest.mark.parametrize("track_index", [True, 1.9, -0.2, "1.9", "not-an-integer"])
def test_normalize_linear_depth_track_index_rejects_non_integer_identity(
    track_index: object,
) -> None:
    with pytest.raises(ValueError, match="invalid track_index"):
        normalize_linear_track_slots(
            [
                LinearTrackSlot(
                    id="depth",
                    renderer="depth",
                    params={"track_index": track_index},
                )
            ]
        )


@pytest.mark.parametrize("track_index", [2, "2", np.int64(2)])
def test_normalize_linear_depth_track_index_accepts_integer_identity(
    track_index: object,
) -> None:
    normalized = normalize_linear_track_slots(
        [
            LinearTrackSlot(
                id="depth",
                renderer="depth",
                params={"track_index": track_index},
            )
        ]
    )

    assert normalized[0].params["track_index"] == 2
    assert isinstance(normalized[0].params["track_index"], int)


def test_linear_resolved_track_retains_spacing_after_px() -> None:
    layout, _canvas_config = _resolve_layout(
        [
            "features:features@side=overlay",
            "depth_1:depth@track_index=0,h=10px,spacing=8px",
        ]
    )

    depth_track = _layout_track(layout, "depth_1")

    assert depth_track.height == pytest.approx(10.0)
    assert depth_track.spacing_after_px == pytest.approx(8.0)


def test_linear_track_slot_geometry_metadata_keeps_duplicate_record_instances() -> None:
    canvas = assemble_linear_diagram_from_records(
        [_record("duplicate"), _record("duplicate")],
        legend="none",
        linear_track_slots=[
            "features:features@side=overlay",
            "gc_content:dinucleotide_content@h=20px,spacing=8px",
        ],
        linear_track_axis_index=0,
    )

    geometry = getattr(canvas, "_gbdraw_track_slot_geometry", None)

    assert geometry["schema"] == 2
    assert geometry["mode"] == "linear"
    assert geometry["source"] == "resolved"
    assert [record["recordIndex"] for record in geometry["records"]] == [0, 1]
    assert [record["recordId"] for record in geometry["records"]] == ["duplicate", "duplicate"]
    gc_slots = [
        next(slot for slot in record["slots"] if slot["slotId"] == "gc_content")
        for record in geometry["records"]
    ]
    assert [slot["slotIndex"] for slot in gc_slots] == [1, 1]
    assert [slot["spacingAfterPx"] for slot in gc_slots] == [pytest.approx(8.0), pytest.approx(8.0)]
    assert all(slot["dataAvailable"] for slot in gc_slots)
    assert all(slot["paintBand"] is not None for slot in gc_slots)
    assert all(slot["reserveBand"] is not None for slot in gc_slots)
    assert all(record["recordBodyBand"] is not None for record in geometry["records"])
    assert all(record["comparisonExclusionBand"] is not None for record in geometry["records"])
    assert all(record["canvasBand"] is not None for record in geometry["records"])


def test_linear_record_plans_expand_only_records_with_extra_feature_lanes() -> None:
    shallow = _record("shallow")
    deep = _record("deep")
    shallow.features = [
        SeqFeature(FeatureLocation(100, 400, strand=-1), type="CDS"),
    ]
    deep.features = [
        SeqFeature(FeatureLocation(100, 400, strand=-1), type="CDS")
        for _ in range(3)
    ]
    canvas = assemble_linear_diagram_from_records(
        [shallow, deep],
        legend="none",
        depth_track_tables=[[_depth_table("shallow")], [_depth_table("deep")]],
        linear_track_slots=[
            "features:features@side=overlay",
            "depth:depth@track_index=0,side=below,h=10px",
        ],
        config_overrides={
            "show_labels": False,
            "show_gc": False,
            "show_skew": False,
            "strandedness": True,
            "resolve_overlaps": True,
        },
    )
    records = canvas._gbdraw_track_slot_geometry["records"]
    features = [
        next(slot for slot in record["slots"] if slot["slotId"] == "features")
        for record in records
    ]
    depth = [
        next(slot for slot in record["slots"] if slot["slotId"] == "depth")
        for record in records
    ]

    assert features[1]["reserveBand"]["bottomPx"] > features[0]["reserveBand"]["bottomPx"]
    assert depth[1]["resolvedOriginPx"] > depth[0]["resolvedOriginPx"]


def test_parse_linear_track_slot_aliases() -> None:
    content = parse_linear_track_slot("gc_content:gc_content@nt=at")
    skew = parse_linear_track_slot("gc_skew:skew@dinucleotide=gc")

    assert content.renderer == "dinucleotide_content"
    assert content.params["nt"] == "AT"
    assert skew.renderer == "dinucleotide_skew"
    assert skew.params["nt"] == "GC"


def test_parse_linear_track_slot_canonicalizes_skew_color_aliases() -> None:
    alias = parse_linear_track_slot(
        "at_skew:dinucleotide_skew@nt=AT,high_color=tomato,low_color=#2a9d8f"
    )
    canonical = parse_linear_track_slot(
        "at_skew:dinucleotide_skew@nt=AT,positive_color=#e76f51,high_color=tomato"
    )

    assert alias.params["positive_color"] == "tomato"
    assert alias.params["negative_color"] == "#2a9d8f"
    assert "high_color" not in alias.params
    assert "low_color" not in alias.params
    assert canonical.params["positive_color"] == "#e76f51"


def test_linear_track_slots_reject_invalid_inputs() -> None:
    with pytest.raises(LinearTrackSlotParseError, match="duplicate"):
        parse_linear_track_slots(["a:features", "a:gc_content"])

    with pytest.raises(LinearTrackSlotParseError, match="unknown"):
        parse_linear_track_slot("bad:not_a_renderer")

    with pytest.raises(LinearTrackSlotParseError, match="generic layout field"):
        parse_linear_track_slots(
            [
                LinearTrackSlot(
                    id="bad",
                    renderer="depth",
                    params={"height": "20px"},
                )
            ]
        )

    with pytest.raises(ValueError, match="overlay"):
        normalize_linear_track_slots_with_axis(
            [LinearTrackSlot(id="bad", renderer="gc_content", side="overlay")],
            None,
        )


def test_normalize_linear_track_slots_with_axis_derives_sides() -> None:
    slots = parse_linear_track_slots(
        [
            "gc_skew:gc_skew@h=28px",
            "features:features",
            "gc_content:gc_content@h=30px",
        ]
    )

    normalized = normalize_linear_track_slots_with_axis(slots, 2)

    assert [slot.side for slot in normalized] == ["above", "above", "below"]
    assert [slot.id for slot in normalized] == ["gc_skew", "features", "gc_content"]


def test_normalize_linear_track_slots_with_axis_rejects_explicit_side_conflicts() -> None:
    slots = parse_linear_track_slots(
        [
            "gc_skew:gc_skew@side=below",
            "features:features",
            "gc_content:gc_content",
        ]
    )

    with pytest.raises(ValueError, match="conflicts with --linear_track_axis_index"):
        normalize_linear_track_slots_with_axis(slots, 1)


def test_normalize_linear_track_slots_with_axis_allows_feature_overlay_at_axis() -> None:
    slots = parse_linear_track_slots(
        [
            "gc_skew:gc_skew",
            "features:features@side=overlay",
            "gc_content:gc_content",
        ]
    )

    normalized = normalize_linear_track_slots_with_axis(slots, 1)

    assert [slot.side for slot in normalized] == ["above", "overlay", "below"]


def test_normalize_linear_track_slots_with_axis_rejects_overlay_off_axis() -> None:
    slots = parse_linear_track_slots(
        [
            "features:features@side=overlay",
            "gc_content:gc_content",
        ]
    )

    with pytest.raises(ValueError, match="conflicts with --linear_track_axis_index"):
        normalize_linear_track_slots_with_axis(slots, 1)


def test_default_and_order_linear_track_slots() -> None:
    defaults = default_linear_track_slots(
        show_depth=True,
        depth_track_count=2,
        show_gc=True,
        show_skew=False,
        track_layout="middle",
    )
    assert [slot.id for slot in defaults] == ["features", "depth_1", "depth_2", "gc_content"]
    assert defaults[0].side == "overlay"
    assert defaults[1].params["track_index"] == 0
    assert defaults[2].params["track_index"] == 1

    ordered = linear_track_slots_from_order(
        "gc_skew,features,gc_content",
        show_gc=True,
        show_skew=True,
        track_layout="above",
    )
    assert [slot.id for slot in ordered] == ["gc_skew", "features", "gc_content"]
    assert ordered[1].side == "above"


def test_resolve_linear_track_layout_clears_depth_above_axis_when_features_are_below() -> None:
    layout, canvas_config = _resolve_layout(
        [
            "depth:depth@side=above",
            "features:features@side=below",
        ]
    )
    depth = _layout_track(layout, "depth")

    assert depth.y_offset + depth.height < 0
    assert 0.0 - (depth.y_offset + depth.height) == pytest.approx(canvas_config.vertical_padding)


def test_resolve_linear_track_layout_clears_depth_below_axis_when_features_are_above() -> None:
    layout, canvas_config = _resolve_layout(
        [
            "features:features@side=above",
            "depth:depth@side=below",
        ]
    )
    depth = _layout_track(layout, "depth")

    assert depth.y_offset > 0
    assert depth.y_offset == pytest.approx(canvas_config.vertical_padding)


def test_resolve_linear_track_layout_does_not_clear_depth_after_same_side_features() -> None:
    layout, _canvas_config = _resolve_layout(
        [
            "features:features@side=below",
            "depth:depth@side=below",
        ]
    )
    features = _layout_track(layout, "features")
    depth = _layout_track(layout, "depth")

    assert features.y_offset == pytest.approx(0.0)
    assert depth.y_offset == pytest.approx(0.0)


def test_resolve_linear_track_layout_clears_numeric_track_below_axis() -> None:
    layout, canvas_config = _resolve_layout(["gc_skew:gc_skew@side=below,h=20px"])
    gc_skew = _layout_track(layout, "gc_skew")

    assert gc_skew.y_offset - gc_skew.top_extent == pytest.approx(canvas_config.vertical_padding)


def test_assemble_linear_custom_slots_places_tracks_above_and_below_axis() -> None:
    svg = assemble_linear_diagram_from_records(
        [_record()],
        legend="none",
        config_overrides={"show_gc": True, "show_skew": True},
        linear_track_slots=[
            "gc_skew:gc_skew@side=above,h=28px",
            "features:features@side=overlay",
            "gc_content:gc_content@side=below,h=30px",
        ],
    ).tostring()

    axis_y = _extract_group_y(svg, "rec1")
    skew_y = _extract_group_y(svg, "gc_skew")
    content_y = _extract_group_y(svg, "gc_content")

    assert skew_y < axis_y
    assert content_y > axis_y


def test_assemble_linear_dinucleotide_skew_slot_uses_custom_colors_and_legend() -> None:
    default_colors = load_default_colors("", palette="default")
    svg = assemble_linear_diagram_from_records(
        [_record()],
        legend="right",
        default_colors=default_colors,
        config_overrides={"show_gc": True, "show_skew": True},
        linear_track_slots=[
            "features:features@side=overlay",
            "gc_skew:dinucleotide_skew@nt=GC,h=20px",
            "at_skew:dinucleotide_skew@nt=AT,h=20px,positive_color=tomato,negative_color=#2a9d8f,legend_label=AT skew",
        ],
    ).tostring()
    at_skew = _extract_group_fragment(svg, "at_skew")
    gc_skew = _extract_group_fragment(svg, "gc_skew")
    default_skew_high = default_colors.loc[default_colors["feature_type"] == "skew_high", "color"].iloc[0]
    default_skew_low = default_colors.loc[default_colors["feature_type"] == "skew_low", "color"].iloc[0]

    assert 'fill="#FF6347"' in at_skew
    assert 'fill="#2a9d8f"' in at_skew
    assert f'fill="{default_skew_high}"' in gc_skew
    assert f'fill="{default_skew_low}"' in gc_skew
    assert "AT skew (+)" in svg
    assert "AT skew (-)" in svg
    assert svg.count('fill="#FF6347"') >= 2
    assert svg.count('fill="#2a9d8f"') >= 2


def test_assemble_linear_dinucleotide_skew_slot_rejects_invalid_custom_color() -> None:
    with pytest.raises(ValidationError, match="Unknown color name"):
        assemble_linear_diagram_from_records(
            [_record()],
            legend="none",
            config_overrides={"show_gc": True, "show_skew": True},
            linear_track_slots=[
                "features:features@side=overlay",
                "at_skew:dinucleotide_skew@nt=AT,positive_color=not-a-color",
            ],
        )


def test_assemble_linear_custom_depth_above_axis_keeps_depth_axis_clear() -> None:
    svg = assemble_linear_diagram_from_records(
        [_record()],
        legend="none",
        depth_table=_depth_table("rec1"),
        linear_track_slots=[
            "depth:depth@side=above,h=10px",
            "features:features@side=below",
        ],
    ).tostring()
    cfg = GbdrawConfig.from_dict(load_config_toml("gbdraw.data", "config.toml"))

    genome_axis_y = _extract_group_y(svg, "rec1")
    depth_group_y = _extract_group_y(svg, "depth")
    depth_axis_y = depth_group_y + _extract_axis_baseline_y(svg, "depth_axis")

    assert depth_axis_y < genome_axis_y
    assert genome_axis_y - depth_axis_y == pytest.approx(cfg.canvas.linear.vertical_padding)


def test_linear_cli_forwards_track_slots_to_api(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    captured: dict[str, object] = {}

    monkeypatch.setattr(linear_cli_module, "load_gbks", lambda *args, **kwargs: [_record()])
    monkeypatch.setattr(linear_cli_module, "read_color_table", lambda _path: None)
    monkeypatch.setattr(linear_cli_module, "load_default_colors", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(linear_cli_module, "read_feature_visibility_file", lambda _path: None)
    monkeypatch.setattr(linear_cli_module, "save_figure", lambda canvas, formats: None)

    def fake_assemble(*args, **kwargs):
        captured.update(kwargs)
        return Drawing(filename=str(tmp_path / "dummy.svg"))

    monkeypatch.setattr(linear_cli_module, "assemble_linear_diagram_from_records", fake_assemble)

    linear_cli_module.linear_main(
        [
            "--gbk",
            "dummy.gb",
            "--show_gc",
            "--show_skew",
            "--linear_track_axis_index",
            "1",
            "--linear_track_slot",
            "gc_skew:gc_skew@h=28px",
            "--linear_track_slot",
            "features:features",
            "--linear_track_slot",
            "gc_content:gc_content@h=30px",
            "--format",
            "svg",
            "-o",
            str(tmp_path / "out"),
        ]
    )

    slots = captured["linear_track_slots"]
    assert [slot.id for slot in slots] == ["gc_skew", "features", "gc_content"]
    assert captured["linear_track_axis_index"] == 1
