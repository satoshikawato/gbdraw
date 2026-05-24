from __future__ import annotations

import copy
import re

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from svgwrite import Drawing

import gbdraw.circular as circular_cli_module
import gbdraw.linear as linear_cli_module
from gbdraw.analysis.gc import gc_content_percent_df
from gbdraw.api.diagram import (
    assemble_circular_diagram_from_record,
    assemble_linear_diagram_from_records,
)
from gbdraw.config.models import GbdrawConfig
from gbdraw.config.toml import load_config_toml
from gbdraw.exceptions import ValidationError


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


def test_gc_content_config_defaults_preserve_deviation_mode() -> None:
    cfg = GbdrawConfig.from_dict(load_config_toml("gbdraw.data", "config.toml"))

    assert cfg.objects.gc_content.mode == "deviation"
    assert cfg.objects.gc_content.min_percent == pytest.approx(0)
    assert cfg.objects.gc_content.max_percent == pytest.approx(100)
    assert cfg.objects.gc_content.large_tick_interval == pytest.approx(20)
    assert cfg.objects.gc_content.small_tick_interval is None
    assert cfg.objects.gc_content.tick_font_size is None


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
        legend="none",
        config_overrides={
            "show_skew": False,
            "gc_content_mode": "percent",
            "gc_content_min_percent": 0,
            "gc_content_max_percent": 100,
            "gc_content_large_tick_interval": 50,
        },
        window=20,
        step=20,
    ).tostring()

    assert 'id="gc_content_axis"' in svg
    assert 'id="depth_axis"' not in svg
    assert ">0%<" in svg
    assert ">50%<" in svg
    assert ">100%<" in svg


@pytest.mark.linear
def test_linear_gc_content_percent_mode_adds_percent_axis_without_depth_axis() -> None:
    svg = assemble_linear_diagram_from_records(
        [_make_record()],
        legend="none",
        config_overrides={
            "show_gc": True,
            "show_skew": False,
            "gc_content_mode": "percent",
            "gc_content_min_percent": 0,
            "gc_content_max_percent": 100,
            "gc_content_large_tick_interval": 50,
        },
        window=20,
        step=20,
    ).tostring()

    assert 'id="gc_content_axis"' in svg
    assert 'id="depth_axis"' not in svg
    assert ">0%<" in svg
    assert ">50%<" in svg
    assert ">100%<" in svg


def test_default_gc_content_deviation_mode_does_not_emit_percent_axis() -> None:
    circular_svg = assemble_circular_diagram_from_record(
        _make_record(),
        legend="none",
        config_overrides={"show_skew": False},
        window=20,
        step=20,
    ).tostring()
    linear_svg = assemble_linear_diagram_from_records(
        [_make_record()],
        legend="none",
        config_overrides={"show_gc": True, "show_skew": False},
        window=20,
        step=20,
    ).tostring()

    assert 'id="gc_content_axis"' not in circular_svg
    assert 'id="gc_content_axis"' not in linear_svg


@pytest.mark.circular
def test_circular_gc_content_percent_axis_inherits_depth_tick_font_size() -> None:
    record = _make_record()

    svg = assemble_circular_diagram_from_record(
        record,
        legend="none",
        config_overrides={
            "show_depth": True,
            "show_skew": False,
            "gc_content_mode": "percent",
            "gc_content_large_tick_interval": 50,
            "depth_tick_font_size": 11,
        },
        depth_table=_make_depth_table(record),
        window=20,
        step=20,
    ).tostring()

    assert _axis_font_sizes(svg, "depth_axis") == {"11.0"}
    assert _axis_font_sizes(svg, "gc_content_axis") == {"11.0"}


@pytest.mark.linear
def test_linear_gc_content_percent_axis_inherits_depth_tick_font_size() -> None:
    record = _make_record()

    svg = assemble_linear_diagram_from_records(
        [record],
        legend="none",
        config_overrides={
            "show_depth": True,
            "show_gc": True,
            "show_skew": False,
            "gc_content_mode": "percent",
            "gc_content_large_tick_interval": 50,
            "depth_tick_font_size": 11,
        },
        depth_table=_make_depth_table(record),
        window=20,
        step=20,
    ).tostring()

    assert _axis_font_sizes(svg, "depth_axis") == {"11.0"}
    assert _axis_font_sizes(svg, "gc_content_axis") == {"11.0"}


def test_circular_cli_gc_percent_options_forward_to_api(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path,
) -> None:
    captured: dict[str, object] = {}

    monkeypatch.setattr(circular_cli_module, "load_gbks", lambda paths, mode: [_make_record()])
    monkeypatch.setattr(circular_cli_module, "read_color_table", lambda _path: None)
    monkeypatch.setattr(circular_cli_module, "load_default_colors", lambda _path, _palette: None)
    monkeypatch.setattr(circular_cli_module, "save_figure", lambda canvas, formats: None)

    def fake_assemble(*args, **kwargs):
        captured.update(kwargs)
        return Drawing(filename=str(tmp_path / "dummy.svg"))

    monkeypatch.setattr(circular_cli_module, "assemble_circular_diagram_from_record", fake_assemble)

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

    cfg = captured["cfg"]
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

    monkeypatch.setattr(linear_cli_module, "load_gbks", lambda *args, **kwargs: [_make_record()])
    monkeypatch.setattr(linear_cli_module, "read_color_table", lambda _path: None)
    monkeypatch.setattr(linear_cli_module, "load_default_colors", lambda *args, **kwargs: None)
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

    cfg = captured["cfg"]
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
