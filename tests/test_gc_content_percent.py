from __future__ import annotations

import copy
import re
import xml.etree.ElementTree as ET
from types import SimpleNamespace

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from svgwrite import Drawing

import gbdraw.circular as circular_cli_module
from gbdraw.api.input_summary import InputSummary
import gbdraw.linear as linear_cli_module
import gbdraw.api.request_render as request_render_module
from gbdraw.analysis.gc import (
    calculate_gc_percent,
    gc_content_percent_df,
)
from gbdraw.analysis.skew import skew_df
from gbdraw.api.config import apply_config_overrides
from gbdraw.api.diagram import (
    assemble_circular_diagram_from_record,
    assemble_linear_diagram_from_records,
)
from gbdraw.config.models import GbdrawConfig, LinearRenderProfile
from gbdraw.config.toml import load_config_toml
from gbdraw.exceptions import ValidationError
from gbdraw.render.groups.circular.definition import definition_record_facts
from gbdraw.svg.circular_ticks import circular_tick_label_margin_px
from gbdraw.svg.circular_tracks import (
    generate_circular_depth_path_desc,
    generate_circular_scalar_area_path_desc,
)
from gbdraw.api.assembly_inputs import CircularAssemblyInputs, LinearAssemblyInputs
from gbdraw.layout.record_coordinates import RecordDisplayTransform

SVG_NAMESPACE = "http://www.w3.org/2000/svg"


def _make_record(record_id: str = "rec1", length: int = 160) -> SeqRecord:
    record = SeqRecord(Seq("ATGC" * (length // 4)), id=record_id, name=record_id)
    record.annotations["molecule_type"] = "DNA"
    record.features = [
        SeqFeature(
            FeatureLocation(10, 50, strand=1),
            type="CDS",
            qualifiers={"product": ["test protein"]},
        )
    ]
    return record


def _make_depth_table(record: SeqRecord) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "reference_name": [record.id, record.id, record.id],
            "position": [1, max(1, len(record.seq) // 2), len(record.seq)],
            "depth": [10, 30, 50],
        }
    )


def _axis_font_sizes(svg: str, axis_id: str) -> set[str]:
    match = re.search(rf'<g id="{re.escape(axis_id)}".*?</g>', svg, re.S)
    assert match is not None
    return set(re.findall(r'font-size="([^"]+)"', match.group(0)))


def _circular_slot_parts(svg: str, slot_id: str) -> dict[str, list[ET.Element]]:
    """The direct children of one Circular slot group, found by its slot id (OV-3084).

    A Circular slot group renames its children's ids, so a test finds the
    percent band by its shape: ``background`` is the band circle, ``axis`` the
    group that carries tick labels, and ``border`` the group of edge circles.
    """

    (group,) = [
        element
        for element in ET.fromstring(svg).iter(f"{{{SVG_NAMESPACE}}}g")
        if element.get("data-gbdraw-slot-id") == slot_id
    ]
    children = list(group)
    groups = [child for child in children if child.tag == f"{{{SVG_NAMESPACE}}}g"]
    return {
        "background": [child for child in children if child.tag == f"{{{SVG_NAMESPACE}}}circle"],
        "axis": [child for child in groups if child.find(f".//{{{SVG_NAMESPACE}}}text") is not None],
        "border": [child for child in groups if child.find(f".//{{{SVG_NAMESPACE}}}text") is None],
    }


def _circular_axis_font_sizes(svg: str, slot_id: str) -> set[str]:
    (axis,) = _circular_slot_parts(svg, slot_id)["axis"]
    return {str(element.get("font-size")) for element in axis.iter() if element.get("font-size") is not None}


def _slot_group_id(svg: str, renderer: str) -> str:
    """The id of the one group that draws ``renderer``, found by its slot tag."""

    (group_id,) = [
        str(element.get("id"))
        for element in ET.fromstring(svg).iter()
        if element.get("data-gbdraw-slot-renderer") == renderer
    ]
    return group_id


def _svg_slot_translate_y(svg: str, renderer: str) -> float:
    root = ET.fromstring(svg)
    group = next(
        element
        for element in root.iter()
        if element.get("data-gbdraw-slot-renderer") == renderer
    )
    translations = re.findall(
        r"translate\(\s*[-+0-9.eE]+(?:px)?[\s,]+([-+0-9.eE]+)(?:px)?\s*\)",
        group.get("transform", ""),
    )
    assert translations
    return sum(float(value) for value in translations)


def test_gc_content_config_defaults_preserve_deviation_mode() -> None:
    cfg = GbdrawConfig.from_dict(load_config_toml("gbdraw.data", "config.toml"))

    assert cfg.objects.gc_content.mode == "deviation"
    assert cfg.objects.gc_content.min_percent == pytest.approx(0)
    assert cfg.objects.gc_content.max_percent == pytest.approx(100)
    assert cfg.objects.gc_content.large_tick_interval == pytest.approx(20)
    assert cfg.objects.gc_content.small_tick_interval is None
    assert cfg.objects.gc_content.tick_font_size is None
    assert cfg.objects.gc_content.percent_background_color == "#e5e5e5"
    assert cfg.objects.gc_content.percent_background_opacity == pytest.approx(1.0)
    assert cfg.objects.gc_content.percent_border_color == "#4b5563"
    assert cfg.objects.gc_content.percent_border_width == pytest.approx(0.8)


def test_gc_content_config_rejects_invalid_mode_ranges_and_ticks() -> None:
    base = load_config_toml("gbdraw.data", "config.toml")

    invalid_mode = copy.deepcopy(base)
    invalid_mode["objects"]["gc_content"]["mode"] = "absolute"
    with pytest.raises(ValidationError, match="gc_content_mode"):
        GbdrawConfig.from_dict(invalid_mode)

    invalid_range = copy.deepcopy(base)
    invalid_range["objects"]["gc_content"]["min_percent"] = 70
    invalid_range["objects"]["gc_content"]["max_percent"] = 30
    with pytest.raises(ValidationError, match="gc_content_min_percent"):
        GbdrawConfig.from_dict(invalid_range)

    invalid_tick = copy.deepcopy(base)
    invalid_tick["objects"]["gc_content"]["large_tick_interval"] = 0
    with pytest.raises(ValidationError, match="gc_content"):
        GbdrawConfig.from_dict(invalid_tick)

    invalid_opacity = copy.deepcopy(base)
    invalid_opacity["objects"]["gc_content"]["percent_background_opacity"] = 1.5
    with pytest.raises(ValidationError, match="gc_content_percent_background_opacity"):
        GbdrawConfig.from_dict(invalid_opacity)

    invalid_border = copy.deepcopy(base)
    invalid_border["objects"]["gc_content"]["percent_border_width"] = -0.1
    with pytest.raises(ValidationError, match="gc_content_percent_border_width"):
        GbdrawConfig.from_dict(invalid_border)


@pytest.mark.parametrize(
    ("sequence", "expected"),
    (("gcgc", 100.0), ("GcGc", 100.0), ("atgc", 50.0), ("AUGC", 50.0), ("", 0.0)),
)
def test_calculate_gc_percent_counts_soft_masked_bases(sequence: str, expected: float) -> None:
    # OV-1150: FASTA-sourced sequences keep their case; soft-masked bases count.
    assert calculate_gc_percent(Seq(sequence)) == expected


def test_definition_gc_percent_agrees_with_the_gc_track_for_a_mixed_case_record() -> None:
    record = SeqRecord(Seq("gcgcATATggCCaatt"), id="soft_masked")
    length = len(record.seq)

    track_df = skew_df(record, window=length, step=length, nt="GC", circular=True)

    assert definition_record_facts(record).gc_percent == 50.0
    assert track_df.loc[0, "GC content"] * 100.0 == pytest.approx(50.0)


def test_gc_windows_are_centred_on_their_sample_and_wrap_at_the_origin() -> None:
    record = SeqRecord(Seq("GGGGAAAA"), id="origin")

    frame = skew_df(record, window=4, step=4, nt="GC", circular=True)

    # The window of sample 0 is bases -2..1: AA (wrapped from the end) and GG.
    assert frame.loc[0, "GC content"] == pytest.approx(0.5)
    assert frame.loc[4, "GC content"] == pytest.approx(0.5)


def test_gc_windows_longer_than_the_record_count_it_per_cycle() -> None:
    record = SeqRecord(Seq("GAAA"), id="short")

    content_df = skew_df(record, window=10, step=10, nt="GC", circular=True)

    assert content_df.loc[0, "GC content"] == pytest.approx(0.3)


def test_gc_content_percent_df_preserves_source_percent_and_clips_normalized_values() -> None:
    source = pd.DataFrame({"GC content": [0.2, 0.5, 0.8]}, index=[0, 10, 20])

    plot_df = gc_content_percent_df(source, dinucleotide="GC", min_percent=30, max_percent=70)

    assert plot_df["position"].tolist() == pytest.approx([0, 10, 20])
    assert plot_df["percent"].tolist() == pytest.approx([20, 50, 80])
    assert plot_df["value"].tolist() == pytest.approx([20, 50, 80])
    assert plot_df["value_normalized"].tolist() == pytest.approx([0, 0.5, 1])


@pytest.mark.circular
def test_circular_gc_content_percent_mode_adds_percent_axis_without_depth_axis() -> None:
    svg = assemble_circular_diagram_from_record(
        _make_record(),
        CircularAssemblyInputs(
            legend="none",
            cfg=apply_config_overrides(None, {
                "canvas.show_skew": False,
                "objects.gc_content.mode": "percent",
                "objects.gc_content.min_percent": 0,
                "objects.gc_content.max_percent": 100,
                "objects.gc_content.large_tick_interval": 50,
                "objects.gc_content.axis_labels": "numbered",
            }),
            window=20,
            step=20,
        ),
    ).tostring()

    parts = _circular_slot_parts(svg, "gc_content")
    assert [len(parts[name]) for name in ("background", "axis", "border")] == [1, 1, 1]
    assert 'data-gbdraw-slot-renderer="depth"' not in svg
    assert ">0%<" in svg
    assert ">50%<" in svg
    assert ">100%<" in svg


@pytest.mark.linear
def test_linear_gc_content_percent_mode_adds_percent_axis_without_depth_axis() -> None:
    svg = assemble_linear_diagram_from_records(
        [_make_record()],
        LinearAssemblyInputs(
            legend="none",
            cfg=apply_config_overrides(None, {
                "canvas.show_gc": True,
                "canvas.show_skew": False,
                "objects.gc_content.mode": "percent",
                "objects.gc_content.min_percent": 0,
                "objects.gc_content.max_percent": 100,
                "objects.gc_content.large_tick_interval": 50,
                "objects.gc_content.axis_labels": "numbered",
            }),
            window=20,
            step=20,
        ),
    ).tostring()

    gc_id = _slot_group_id(svg, "dinucleotide_content")
    assert f'id="{gc_id}_axis"' in svg
    assert f'id="{gc_id}_percent_background"' in svg
    assert f'id="{gc_id}_percent_border"' in svg
    assert 'data-gbdraw-slot-renderer="depth"' not in svg
    assert ">0%<" in svg
    assert ">50%<" in svg
    assert ">100%<" in svg


@pytest.mark.linear
def test_linear_gc_content_percent_mode_reserves_full_height_before_gc_skew() -> None:
    from gbdraw.canvas import LinearCanvasConfigurator
    from gbdraw.config.modify import modify_config_dict

    config_dict = modify_config_dict(
        load_config_toml('gbdraw.data', 'config.toml'),
        {
            "canvas.show_gc": True,
            "canvas.show_skew": True,
            "objects.gc_content.mode": 'percent',
        },
    )
    cfg = GbdrawConfig.from_dict(config_dict)
    canvas_config = LinearCanvasConfigurator(
        num_of_entries=1,
        longest_genome=160,
        profile=LinearRenderProfile(cfg),
        legend="none",
    )

    gc_bottom = canvas_config.gc_content_track_offset + canvas_config.gc_height
    skew_top = canvas_config.gc_skew_track_offset - (0.5 * canvas_config.skew_height)

    assert skew_top == pytest.approx(gc_bottom)

    canvas = assemble_linear_diagram_from_records(
        [_make_record()],
        LinearAssemblyInputs(
            legend="none",
            cfg=apply_config_overrides(None, {
                "canvas.show_gc": True,
                "canvas.show_skew": True,
                "objects.gc_content.mode": "percent",
            }),
            window=20,
            step=20,
        ),
    )
    svg = canvas.tostring()
    gc_y = _svg_slot_translate_y(svg, "dinucleotide_content")
    skew_y = _svg_slot_translate_y(svg, "dinucleotide_skew")
    slots = {
        slot["slotId"]: slot
        for slot in canvas._gbdraw_track_slot_geometry["records"][0]["slots"]
    }

    assert slots["gc_skew"]["reserveBand"]["topPx"] == pytest.approx(
        slots["gc_content"]["reserveBand"]["bottomPx"]
    )
    assert skew_y > gc_y + canvas_config.gc_height


def test_default_gc_content_deviation_mode_does_not_emit_percent_axis() -> None:
    circular_svg = assemble_circular_diagram_from_record(
        _make_record(),
        CircularAssemblyInputs(
            legend="none",
            cfg=apply_config_overrides(None, {"canvas.show_skew": False}),
            window=20,
            step=20,
        ),
    ).tostring()
    linear_svg = assemble_linear_diagram_from_records(
        [_make_record()],
        LinearAssemblyInputs(
            legend="none",
            cfg=apply_config_overrides(
                None, {"canvas.show_gc": True, "canvas.show_skew": False}
            ),
            window=20,
            step=20,
        ),
    ).tostring()

    assert _circular_slot_parts(circular_svg, "gc_content") == {"background": [], "axis": [], "border": []}
    assert f'id="{_slot_group_id(linear_svg, "dinucleotide_content")}_axis"' not in linear_svg
    assert "percent_background" not in circular_svg
    assert "percent_border" not in circular_svg
    assert "percent_background" not in linear_svg
    assert "percent_border" not in linear_svg


@pytest.mark.circular
def test_circular_gc_content_percent_container_ids_are_slot_scoped() -> None:
    svg = assemble_circular_diagram_from_record(
        _make_record(),
        CircularAssemblyInputs(
            legend="none",
            circular_track_slots=[
                "gc_content:dinucleotide_content@nt=GC,w=20px",
                "at_content:dinucleotide_content@nt=AT,w=20px",
            ],
            cfg=apply_config_overrides(None, {
                "objects.gc_content.mode": "percent",
                "objects.gc_content.show_axis": False,
            }),
            window=20,
            step=20,
        ),
    ).tostring()

    root = ET.fromstring(svg)
    groups = {
        element.attrib["data-gbdraw-slot-id"]: element
        for element in root.iter()
        if element.attrib.get("data-gbdraw-slot-id")
    }
    child_ids_by_slot = {
        slot_id: {
            child.attrib["id"]
            for child in group.iter()
            if child is not group and child.attrib.get("id")
        }
        for slot_id, group in groups.items()
    }

    assert set(child_ids_by_slot) == {"gc_content", "at_content"}
    assert all(len(child_ids) == 2 for child_ids in child_ids_by_slot.values())
    assert child_ids_by_slot["gc_content"].isdisjoint(
        child_ids_by_slot["at_content"]
    )
    assert all(
        child_id.startswith("track_slot_child_")
        for child_ids in child_ids_by_slot.values()
        for child_id in child_ids
    )


def test_circular_scalar_area_closure_is_opt_in_for_depth_paths() -> None:
    depth_plot_df = pd.DataFrame(
        {
            "position": [0, 50],
            "depth_normalized": [0.25, 0.75],
        }
    )

    default_path = generate_circular_scalar_area_path_desc(
        100,
        100,
        depth_plot_df,
        20,
        1.0,
        value_column="depth_normalized", record_transform=_circular_transform(100)
    )
    depth_path = generate_circular_depth_path_desc(100, 100, depth_plot_df, 20, 1.0, record_transform=_circular_transform(100))
    closed_path = generate_circular_scalar_area_path_desc(
        100,
        100,
        depth_plot_df,
        20,
        1.0,
        value_column="depth_normalized",
        close_at_record_len=True, record_transform=_circular_transform(100)
    )

    assert depth_path == default_path
    assert closed_path != default_path


@pytest.mark.circular
def test_circular_gc_content_percent_axis_inherits_depth_tick_font_size() -> None:
    record = _make_record()

    svg = assemble_circular_diagram_from_record(
        record,
        CircularAssemblyInputs(
            legend="none",
            cfg=apply_config_overrides(None, {
                "canvas.show_depth": True,
                "canvas.show_skew": False,
                "objects.gc_content.mode": "percent",
                "objects.gc_content.large_tick_interval": 50,
                "objects.depth.tick_font_size": 11,
            }),
            depth_table=_make_depth_table(record),
            window=20,
            step=20,
        ),
    ).tostring()

    assert _circular_axis_font_sizes(svg, "depth") == {"11.0"}
    assert _circular_axis_font_sizes(svg, "gc_content") == {"11.0"}


@pytest.mark.linear
def test_linear_gc_content_percent_axis_inherits_depth_tick_font_size() -> None:
    record = _make_record()

    svg = assemble_linear_diagram_from_records(
        [record],
        LinearAssemblyInputs(
            legend="none",
            cfg=apply_config_overrides(None, {
                "canvas.show_depth": True,
                "canvas.show_gc": True,
                "canvas.show_skew": False,
                "objects.gc_content.mode": "percent",
                "objects.gc_content.large_tick_interval": 50,
                "objects.depth.tick_font_size": 11,
            }),
            depth_table=_make_depth_table(record),
            window=20,
            step=20,
        ),
    ).tostring()

    assert _axis_font_sizes(svg, f"{_slot_group_id(svg, 'depth')}_axis") == {"11.0"}
    assert _axis_font_sizes(svg, f"{_slot_group_id(svg, 'dinucleotide_content')}_axis") == {"11.0"}


def test_circular_cli_gc_percent_options_forward_to_api(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path,
    stub_typed_request_export,
) -> None:
    captured: dict[str, object] = {}

    monkeypatch.setattr(request_render_module, "load_gbks", lambda paths, **_kwargs: [_make_record()])
    monkeypatch.setattr(request_render_module, "read_color_table", lambda _path: None)

    def fake_assemble(*args, **kwargs):
        captured["options"] = kwargs["options"]
        return Drawing(filename=str(tmp_path / "dummy.svg"))

    monkeypatch.setattr(request_render_module, "build_circular_diagram", fake_assemble)

    circular_cli_module.circular_main(
        [
            "--gbk",
            "dummy.gb",
            "--gc_content_mode",
            "percent",
            "--gc_content_min_percent",
            "30",
            "--gc_content_max_percent",
            "70",
            "--gc_content_tick_interval",
            "5",
            "--gc_content_small_tick_interval",
            "2.5",
            "--gc_content_tick_font_size",
            "9",
            "--hide_gc_content_axis",
            "--hide_gc_content_ticks",
            "--format",
            "svg",
            "-o",
            str(tmp_path / "out"),
        ]
    )

    cfg = captured["options"].config
    assert isinstance(cfg, GbdrawConfig)
    assert cfg.objects.gc_content.mode == "percent"
    assert cfg.objects.gc_content.min_percent == pytest.approx(30)
    assert cfg.objects.gc_content.max_percent == pytest.approx(70)
    assert cfg.objects.gc_content.large_tick_interval == pytest.approx(5)
    assert cfg.objects.gc_content.small_tick_interval == pytest.approx(2.5)
    assert cfg.objects.gc_content.tick_font_size == pytest.approx(9)
    assert cfg.objects.gc_content.show_axis is False
    assert cfg.objects.gc_content.show_ticks is False


def test_linear_cli_gc_percent_options_forward_to_api(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path,
) -> None:
    captured: dict[str, object] = {}

    monkeypatch.setattr(request_render_module, "load_gbks", lambda *args, **kwargs: [_make_record()])
    monkeypatch.setattr(request_render_module, "read_color_table", lambda _path: None)

    def fake_render_request(request, **_kwargs):
        resolved = request_render_module.resolve_request(request)
        captured["canonical_request"] = resolved
        return SimpleNamespace(
            drawing=Drawing(filename=str(tmp_path / "dummy.svg")),
            interactive_context=None,
            records=tuple(item.source.record for item in resolved.records),
            losat_cache_entries=(),
            losat_derived_cache_entries=(),
            protein_identity_manifest=None,
            request=resolved,
            annotation_warnings=(),
            feature_identity_notices=(),
            input_summary=InputSummary(()),
            linear_metadata=None,
        )

    monkeypatch.setattr(linear_cli_module, "render_request", fake_render_request)

    linear_cli_module.linear_main(
        [
            "--gbk",
            "dummy.gb",
            "--gc",
            "--gc_content_mode",
            "percent",
            "--gc_content_min_percent",
            "35",
            "--gc_content_max_percent",
            "65",
            "--gc_content_large_tick_interval",
            "10",
            "--gc_content_small_tick_interval",
            "5",
            "--gc_content_tick_font_size",
            "8",
            "--hide_gc_content_axis",
            "--hide_gc_content_ticks",
            "--format",
            "svg",
            "-o",
            str(tmp_path / "out"),
        ]
    )

    cfg = captured["canonical_request"].options.config
    assert isinstance(cfg, GbdrawConfig)
    assert cfg.objects.gc_content.mode == "percent"
    assert cfg.objects.gc_content.min_percent == pytest.approx(35)
    assert cfg.objects.gc_content.max_percent == pytest.approx(65)
    assert cfg.objects.gc_content.large_tick_interval == pytest.approx(10)
    assert cfg.objects.gc_content.small_tick_interval == pytest.approx(5)
    assert cfg.objects.gc_content.tick_font_size == pytest.approx(8)
    assert cfg.objects.gc_content.show_axis is False
    assert cfg.objects.gc_content.show_ticks is False


def test_gc_percent_cli_rejects_invalid_ranges() -> None:
    with pytest.raises(SystemExit):
        circular_cli_module._get_args(
            [
                "--gbk",
                "dummy.gb",
                "--gc_content_mode",
                "percent",
                "--gc_content_min_percent",
                "70",
                "--gc_content_max_percent",
                "30",
            ]
        )
    with pytest.raises(SystemExit):
        linear_cli_module._get_args(
            [
                "--gbk",
                "dummy.gb",
                "--gc_content_mode",
                "percent",
                "--gc_content_min_percent",
                "70",
                "--gc_content_max_percent",
                "30",
            ]
        )


def _circular_transform(length: int) -> RecordDisplayTransform:
    """The transform of a whole record drawn as circular: its lines close round the origin."""

    return RecordDisplayTransform(int(length), 1, 1, None, True)


def _tick_label_path_radii(svg: str) -> list[float]:
    """The radii of the circles the scale's tick labels are written along (in the scale's group)."""

    ticks = _slot_group(svg, "ticks")
    return [
        float(match.group(1))
        for path in ticks.iter()
        if path.tag.endswith("path") and str(path.get("id", "")).startswith("circular_tick_label_path_")
        if (match := re.search(r" A ([0-9.]+),", str(path.get("d", "")))) is not None
    ]


def _slot_group(svg: str, slot_id: str) -> ET.Element:
    """The group a track slot draws into, found by its slot id."""

    group = next(
        (element for element in ET.fromstring(svg).iter() if element.get("data-gbdraw-slot-id") == slot_id), None
    )
    assert group is not None
    return group


def _percent_outer_edge_px(svg: str) -> float:
    """Where the Percent GC band's paint ends outward: its background ring's outer edge."""

    # The band's background ring is the first child of the gc_content slot's group.
    background = _slot_group(svg, "gc_content")[0]
    assert background.tag.endswith("circle")
    return float(background.get("r", "nan")) + 0.5 * float(background.get("stroke-width", "nan"))


@pytest.mark.circular
@pytest.mark.parametrize(
    ("paint", "clearance"),
    [
        ({}, True),
        # A band that paints no edge (no background, no border) keeps the tick labels' band as before.
        ({"objects.gc_content.percent_background_color": "none", "objects.gc_content.percent_border_width": 0}, False),
    ],
)
def test_the_tick_labels_keep_their_margin_from_a_percent_gc_band_that_paints_its_edge(paint, clearance) -> None:
    # OV-3116: the tick labels inside the ticks stood on the Percent GC band's painted outer edge, the descenders of
    # "Mbp" on its border line. Inside the ticks the labels stand on their path (above it at the top, below it at
    # the bottom), so the path is their inner edge.
    overrides = {"objects.gc_content.mode": "percent", **paint}
    svg = assemble_circular_diagram_from_record(
        _make_record(length=1_600_000),
        CircularAssemblyInputs(
            legend="none", cfg=apply_config_overrides(None, overrides), window=20_000, step=20_000
        ),
    ).tostring()
    cfg = apply_config_overrides(None, overrides)
    margin = circular_tick_label_margin_px(
        float(cfg.objects.ticks.tick_labels.font_size), float(cfg.objects.ticks.tick_width)
    )
    radii = _tick_label_path_radii(svg)
    assert radii
    assert min(radii) - _percent_outer_edge_px(svg) == pytest.approx(margin if clearance else 0.0, abs=1e-3)


@pytest.mark.parametrize(
    ("words", "drawn"),
    [([], 1), (["--multi_record_canvas"], 3), (["--one_circle"], 3)],
    ids=["one-record", "multi-record-canvas", "one-circle"],
)
def test_a_cli_render_takes_the_gc_percent_of_each_drawn_record_once(
    words: list[str], drawn: int, tmp_path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The definition's facts, GC % included, are built once per record (once per circle on one circle) and handed
    to every reader: the band the tracks keep clear, the fit probes, the centre text and the drawn group (OV-3043:
    the GC % of a record was counted once per reader)."""

    import gbdraw.analysis.gc as gc_module
    from Bio import SeqIO
    from gbdraw.circular import circular_main

    counted: list[int] = []
    real = gc_module.known_base_counts

    def known_base_counts(seq):  # type: ignore[no-untyped-def]
        counted.append(len(seq))
        return real(seq)

    monkeypatch.setattr(gc_module, "known_base_counts", known_base_counts)
    records = [SeqRecord(Seq("GCAT" * (3_000 + 500 * index)), id=f"gc{index}", name=f"gc{index}") for index in range(3)]
    for record in records:
        record.annotations["molecule_type"] = "DNA"
    source = tmp_path / "records.gb"
    SeqIO.write(records if drawn > 1 else records[:1], source, "genbank")
    circular_main(["--gbk", str(source), *words, "-o", str(tmp_path / "out"), "-f", "svg"])

    assert sorted(counted) == sorted(len(record.seq) for record in records[:drawn])
