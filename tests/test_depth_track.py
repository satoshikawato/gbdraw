from __future__ import annotations

import re
from pathlib import Path

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from svgwrite import Drawing

import gbdraw.circular as circular_cli_module
import gbdraw.linear as linear_cli_module
from gbdraw.analysis.depth import depth_df, read_depth_tsv
from gbdraw.api.diagram import (
    assemble_circular_diagram_from_records,
    assemble_circular_diagram_from_record,
    assemble_linear_diagram_from_records,
)
from gbdraw.exceptions import ValidationError


def _make_record(record_id: str = "rec1", length: int = 120) -> SeqRecord:
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


def _depth_table(reference: str = "rec1") -> pd.DataFrame:
    return pd.DataFrame(
        {
            "reference_name": [reference, reference, reference, reference],
            "position": [1, 2, 10, 30],
            "depth": [10, 20, 40, 80],
        }
    )


def _constant_depth_table(reference: str, depth: float, length: int = 40) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "reference_name": [reference] * length,
            "position": list(range(1, length + 1)),
            "depth": [depth] * length,
        }
    )


def _write_depth_file(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def _svg_group_translate_y(svg: str, group_id: str) -> float:
    match = re.search(rf'<g id="{re.escape(group_id)}" transform="translate\([^,]+,([^)]+)\)"', svg)
    assert match is not None
    return float(match.group(1))


def test_read_depth_tsv_headerless_and_headered(tmp_path: Path) -> None:
    headerless = _write_depth_file(tmp_path / "depth.tsv", "rec1\t1\t10\nrec1\t2\t12\n")
    headered = _write_depth_file(
        tmp_path / "depth_header.tsv",
        "reference\tposition\tdepth\nrec1\t1\t10\nrec1\t2\t12\n",
    )

    assert read_depth_tsv(str(headerless)).to_dict("list") == {
        "reference_name": ["rec1", "rec1"],
        "position": [1, 2],
        "depth": [10.0, 12.0],
    }
    assert read_depth_tsv(str(headered)).to_dict("list") == {
        "reference_name": ["rec1", "rec1"],
        "position": [1, 2],
        "depth": [10.0, 12.0],
    }


def test_depth_df_matching_fallback_and_mismatch() -> None:
    record = _make_record("rec1")
    exact = depth_df(record, _depth_table("rec1"), window=10, step=10)
    fallback = depth_df(record, _depth_table("other"), window=10, step=10)

    assert exact["depth"].iloc[0] == pytest.approx(7.0)
    assert fallback["depth"].iloc[0] == pytest.approx(7.0)

    multi_ref = pd.DataFrame(
        {
            "reference_name": ["other1", "other2"],
            "position": [1, 1],
            "depth": [1, 2],
        }
    )
    with pytest.raises(ValidationError, match="do not match"):
        depth_df(record, multi_ref, window=10, step=10)


def test_depth_df_missing_positions_clipping_and_normalization() -> None:
    record = _make_record("rec1", length=40)
    df = depth_df(
        record,
        _depth_table("rec1"),
        window=10,
        step=10,
        normalize=True,
        min_depth=5,
        max_depth=20,
    )

    assert df["position"].tolist() == [0, 10, 20, 30]
    assert df["depth"].tolist() == pytest.approx([7.0, 5.0, 8.0, 5.0])
    assert df["depth_normalized"].min() >= 0
    assert df["depth_normalized"].max() <= 1


def test_depth_df_log_scale_makes_ten_x_twice_one_x() -> None:
    record = _make_record("rec1", length=20)
    depth_table = pd.DataFrame(
        {
            "reference_name": ["rec1"] * 20,
            "position": list(range(1, 21)),
            "depth": ([1] * 10) + ([10] * 10),
        }
    )

    log_scaled = depth_df(record, depth_table, window=10, step=10, normalize=True)
    linear = depth_df(record, depth_table, window=10, step=10)

    assert log_scaled["depth_normalized"].tolist() == pytest.approx([0.5, 1.0])
    assert linear["depth_normalized"].tolist() == pytest.approx([1.0, 10.0])


def test_depth_config_defaults_to_linear_scale() -> None:
    from gbdraw.config.models import GbdrawConfig
    from gbdraw.config.toml import load_config_toml

    cfg = GbdrawConfig.from_dict(load_config_toml("gbdraw.data", "config.toml"))

    assert cfg.objects.depth.normalize is False
    assert cfg.objects.depth.show_axis is True


@pytest.mark.linear
def test_linear_depth_gc_track_offsets_reserve_visual_gap() -> None:
    from gbdraw.canvas import LinearCanvasConfigurator
    from gbdraw.config.models import GbdrawConfig
    from gbdraw.config.modify import modify_config_dict
    from gbdraw.config.toml import load_config_toml

    config_dict = modify_config_dict(
        load_config_toml("gbdraw.data", "config.toml"),
        show_gc=True,
        show_skew=True,
        show_depth=True,
    )
    cfg = GbdrawConfig.from_dict(config_dict)
    canvas_config = LinearCanvasConfigurator(
        num_of_entries=1,
        longest_genome=120,
        config_dict=config_dict,
        legend="none",
        cfg=cfg,
    )

    depth_bottom = canvas_config.depth_track_offset + canvas_config.depth_height
    gc_top = canvas_config.gc_content_track_offset - (0.5 * canvas_config.gc_height)
    skew_top = canvas_config.gc_skew_track_offset - (0.5 * canvas_config.skew_height)
    gc_bottom = canvas_config.gc_content_track_offset + (0.5 * canvas_config.gc_height)

    assert gc_top - depth_bottom == pytest.approx(cfg.canvas.linear.depth_padding)
    assert skew_top == pytest.approx(gc_bottom)
    assert canvas_config.depth_padding + canvas_config.gc_padding + canvas_config.skew_padding == pytest.approx(
        canvas_config.plot_tracks_height
    )

    svg = assemble_linear_diagram_from_records(
        [_make_record()],
        legend="none",
        depth_table=_depth_table("rec1"),
        config_overrides={"show_gc": True, "show_skew": True},
    ).tostring()
    depth_y = _svg_group_translate_y(svg, "depth")
    gc_y = _svg_group_translate_y(svg, "gc_content")
    skew_y = _svg_group_translate_y(svg, "gc_skew")

    assert (gc_y - (0.5 * canvas_config.gc_height)) - (depth_y + canvas_config.depth_height) == pytest.approx(
        cfg.canvas.linear.depth_padding
    )
    assert skew_y - (0.5 * canvas_config.skew_height) == pytest.approx(gc_y + (0.5 * canvas_config.gc_height))


@pytest.mark.circular
def test_circular_depth_track_is_optional() -> None:
    record = _make_record()

    without_depth = assemble_circular_diagram_from_record(
        record,
        legend="none",
        config_overrides={"show_gc": False, "show_skew": False},
    ).tostring()
    with_depth = assemble_circular_diagram_from_record(
        record,
        legend="none",
        depth_table=_depth_table("rec1"),
        config_overrides={"show_gc": True, "show_skew": True},
    ).tostring()

    assert 'id="depth"' not in without_depth
    assert 'id="depth"' in with_depth
    assert 'id="depth_axis"' in with_depth
    assert ">0x<" not in with_depth
    assert ">1.5x<" in with_depth
    assert 'id="gc_content"' in with_depth


@pytest.mark.circular
def test_circular_depth_origin_label_requires_explicit_min() -> None:
    record = _make_record()
    svg = assemble_circular_diagram_from_record(
        record,
        legend="none",
        depth_table=_depth_table("rec1"),
        config_overrides={
            "show_gc": False,
            "show_skew": False,
            "depth_min": 5,
        },
        window=10,
        step=10,
        depth_window=10,
        depth_step=10,
    ).tostring()

    depth_axis_svg = svg[svg.index('id="depth_axis"') : svg.index('id="depth_axis"') + 900]
    assert ">5x<" in depth_axis_svg
    assert ">0x<" not in depth_axis_svg


@pytest.mark.circular
def test_circular_depth_axis_can_be_hidden_without_hiding_depth_track() -> None:
    record = _make_record()
    svg = assemble_circular_diagram_from_record(
        record,
        legend="none",
        depth_table=_depth_table("rec1"),
        config_overrides={
            "show_gc": False,
            "show_skew": False,
            "depth_show_axis": False,
        },
        window=10,
        step=10,
    ).tostring()

    assert 'id="depth"' in svg
    assert 'id="depth_axis"' not in svg
    assert ">1.5x<" not in svg


@pytest.mark.circular
def test_circular_depth_compresses_gc_skew_to_preserve_definition_space(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import gbdraw.diagrams.circular.assemble as circular_assemble_module
    from gbdraw.config.models import GbdrawConfig
    from gbdraw.config.toml import load_config_toml

    record = _make_record(length=120)
    captured: dict[str, float] = {}

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
        cfg=None,
    ):
        assert track_width_override is not None
        assert norm_factor_override is not None
        captured["gc_width"] = float(track_width_override)
        captured["gc_center"] = float(norm_factor_override) * float(canvas_config.radius)
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
        cfg=None,
    ):
        assert track_width_override is not None
        assert norm_factor_override is not None
        captured["skew_width"] = float(track_width_override)
        captured["skew_center"] = float(norm_factor_override) * float(canvas_config.radius)
        return canvas

    monkeypatch.setattr(
        circular_assemble_module,
        "add_gc_content_group_on_canvas",
        fake_add_gc_content_group_on_canvas,
    )
    monkeypatch.setattr(
        circular_assemble_module,
        "add_gc_skew_group_on_canvas",
        fake_add_gc_skew_group_on_canvas,
    )

    assemble_circular_diagram_from_record(
        record,
        legend="none",
        depth_table=_depth_table("rec1"),
        config_overrides={"show_gc": True, "show_skew": True},
        window=10,
        step=10,
    )

    cfg = GbdrawConfig.from_dict(load_config_toml("gbdraw.data", "config.toml"))
    length_param = "short"
    base_radius = float(cfg.canvas.circular.radius)
    track_ratio = float(cfg.canvas.circular.track_ratio)
    track_type = str(cfg.canvas.circular.track_type)
    track_dict = cfg.canvas.circular.track_dict[length_param][track_type]
    default_gc_width = base_radius * track_ratio * float(cfg.canvas.circular.track_ratio_factors[length_param][1])
    default_skew_width = base_radius * track_ratio * float(cfg.canvas.circular.track_ratio_factors[length_param][2])

    no_depth_skew_inner = (base_radius * float(track_dict["3"])) - (0.5 * default_skew_width)
    actual_skew_inner = captured["skew_center"] - (0.5 * captured["skew_width"])
    assert actual_skew_inner == pytest.approx(no_depth_skew_inner)
    assert captured["gc_width"] < default_gc_width
    assert captured["skew_width"] < default_skew_width

    depth_width = default_gc_width * 0.5
    depth_inner = (base_radius * float(track_dict["2"])) - (0.5 * depth_width)
    actual_gc_outer = captured["gc_center"] + (0.5 * captured["gc_width"])
    assert actual_gc_outer < depth_inner


@pytest.mark.circular
def test_circular_depth_preserves_explicit_gc_skew_track_specs(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import gbdraw.diagrams.circular.assemble as circular_assemble_module

    record = _make_record(length=120)
    captured: dict[str, float] = {}

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
        cfg=None,
    ):
        captured["gc_width"] = float(track_width_override)
        captured["gc_norm"] = float(norm_factor_override)
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
        cfg=None,
    ):
        captured["skew_width"] = float(track_width_override)
        captured["skew_norm"] = float(norm_factor_override)
        return canvas

    monkeypatch.setattr(
        circular_assemble_module,
        "add_gc_content_group_on_canvas",
        fake_add_gc_content_group_on_canvas,
    )
    monkeypatch.setattr(
        circular_assemble_module,
        "add_gc_skew_group_on_canvas",
        fake_add_gc_skew_group_on_canvas,
    )

    assemble_circular_diagram_from_record(
        record,
        legend="none",
        depth_table=_depth_table("rec1"),
        config_overrides={"show_gc": True, "show_skew": True},
        track_specs=["gc_content@r=0.44,w=12px", "gc_skew@r=0.24,w=10px"],
        window=10,
        step=10,
    )

    assert captured["gc_width"] == pytest.approx(12.0)
    assert captured["gc_norm"] == pytest.approx(0.44)
    assert captured["skew_width"] == pytest.approx(10.0)
    assert captured["skew_norm"] == pytest.approx(0.24)


@pytest.mark.circular
def test_circular_depth_axis_uses_record_start_and_tick_options() -> None:
    record = _make_record()
    svg = assemble_circular_diagram_from_record(
        record,
        legend="none",
        depth_table=_depth_table("rec1"),
        config_overrides={
            "show_gc": False,
            "show_skew": False,
            "depth_tick_interval": 4,
            "depth_tick_font_size": 11,
        },
        window=10,
        step=10,
        depth_window=10,
        depth_step=10,
    ).tostring()

    depth_axis_svg = svg[svg.index('id="depth_axis"') : svg.index('id="depth_axis"') + 900]
    assert re.search(r'<line[^>]*x1="0" x2="0" y1="-[^"]+" y2="-[^"]+"', depth_axis_svg)
    assert 'font-size="11.0"' in depth_axis_svg
    assert ">4x<" in depth_axis_svg


@pytest.mark.circular
def test_circular_depth_ticks_use_right_side_and_small_ticks_are_unlabeled() -> None:
    record = _make_record(length=40)
    svg = assemble_circular_diagram_from_record(
        record,
        legend="none",
        depth_table=_constant_depth_table("rec1", 100, length=40),
        config_overrides={
            "show_gc": False,
            "show_skew": False,
            "depth_normalize": False,
            "depth_min": 0,
            "depth_max": 100,
            "depth_large_tick_interval": 50,
            "depth_small_tick_interval": 25,
        },
        window=10,
        step=10,
    ).tostring()

    depth_axis_svg = svg[svg.index('id="depth_axis"') : svg.index('id="depth_axis"') + 1400]
    assert 'x1="-3.0"' not in depth_axis_svg
    assert re.search(r'<line[^>]*x1="0" x2="3.0" y1="-[^"]+" y2="-[^"]+"', depth_axis_svg)
    assert re.search(r'<line[^>]*x1="0" x2="2.0" y1="-[^"]+" y2="-[^"]+"', depth_axis_svg)
    assert ">50x<" in depth_axis_svg
    assert ">25x<" not in depth_axis_svg
    assert ">75x<" not in depth_axis_svg


@pytest.mark.linear
def test_depth_ticks_can_be_hidden_without_hiding_axes() -> None:
    record = _make_record()
    svg = assemble_linear_diagram_from_records(
        [record],
        legend="none",
        depth_table=_depth_table("rec1"),
        config_overrides={
            "show_gc": False,
            "show_skew": False,
            "depth_show_ticks": False,
        },
        window=10,
        step=10,
    ).tostring()

    assert 'id="depth_axis"' in svg
    assert ">0x<" not in svg


@pytest.mark.linear
def test_depth_axis_can_be_hidden_without_hiding_depth_track() -> None:
    record = _make_record()
    svg = assemble_linear_diagram_from_records(
        [record],
        legend="none",
        depth_table=_depth_table("rec1"),
        config_overrides={
            "show_gc": False,
            "show_skew": False,
            "depth_show_axis": False,
        },
        window=10,
        step=10,
    ).tostring()

    assert 'id="depth"' in svg
    assert 'id="depth_axis"' not in svg
    assert ">1.5x<" not in svg


@pytest.mark.linear
def test_linear_shared_depth_axis_uses_common_auto_max() -> None:
    records = [_make_record("rec1", length=40), _make_record("rec2", length=40)]
    svg = assemble_linear_diagram_from_records(
        records,
        legend="none",
        depth_tables=[
            _constant_depth_table("rec1", 10, length=40),
            _constant_depth_table("rec2", 100, length=40),
        ],
        config_overrides={
            "show_gc": False,
            "show_skew": False,
            "depth_share_axis": True,
            "depth_normalize": False,
        },
        window=10,
        step=10,
    ).tostring()

    assert svg.count(">100x<") == 2


@pytest.mark.circular
def test_circular_multi_record_shared_depth_axis_uses_common_auto_max() -> None:
    records = [_make_record("rec1", length=40), _make_record("rec2", length=40)]
    svg = assemble_circular_diagram_from_records(
        records,
        legend="none",
        depth_tables=[
            _constant_depth_table("rec1", 10, length=40),
            _constant_depth_table("rec2", 100, length=40),
        ],
        config_overrides={
            "show_gc": False,
            "show_skew": False,
            "depth_share_axis": True,
            "depth_normalize": False,
        },
        window=10,
        step=10,
    ).tostring()

    assert svg.count(">100x<") == 2


@pytest.mark.linear
def test_linear_depth_track_is_optional() -> None:
    record = _make_record()

    without_depth = assemble_linear_diagram_from_records(
        [record],
        legend="none",
        config_overrides={"show_gc": False, "show_skew": False},
    ).tostring()
    with_depth = assemble_linear_diagram_from_records(
        [record],
        legend="none",
        depth_table=_depth_table("rec1"),
        config_overrides={"show_gc": True, "show_skew": True},
    ).tostring()

    assert 'id="depth"' not in without_depth
    assert 'id="depth"' in with_depth
    assert 'id="depth_axis"' in with_depth
    assert ">0x<" not in with_depth
    assert ">1.5x<" in with_depth
    assert 'id="gc_content"' in with_depth
    assert with_depth.index('id="depth"') < with_depth.index('id="gc_content"')


@pytest.mark.linear
def test_linear_depth_origin_label_requires_explicit_min() -> None:
    record = _make_record()
    svg = assemble_linear_diagram_from_records(
        [record],
        legend="none",
        depth_table=_depth_table("rec1"),
        config_overrides={
            "show_gc": False,
            "show_skew": False,
            "depth_min": 5,
        },
        window=10,
        step=10,
        depth_window=10,
        depth_step=10,
    ).tostring()

    depth_axis_svg = svg[svg.index('id="depth_axis"') : svg.index('id="depth_axis"') + 900]
    assert ">5x<" in depth_axis_svg
    assert ">0x<" not in depth_axis_svg


@pytest.mark.linear
def test_linear_depth_defaults_to_linear_scale_for_plotting() -> None:
    record = _make_record(length=40)
    svg = assemble_linear_diagram_from_records(
        [record],
        legend="none",
        depth_table=_depth_table("rec1"),
        config_overrides={
            "show_gc": False,
            "show_skew": False,
        },
        window=10,
        step=10,
        depth_window=10,
        depth_step=10,
    ).tostring()

    match = re.search(r'id="depth"[^>]*><path d="([^"]+)"', svg)
    assert match is not None
    depth_path = match.group(1)
    assert "L0.0 1.25" in depth_path
    assert "L500.0 10.0" in depth_path
    assert "L1000.0 0.0" in depth_path


@pytest.mark.linear
def test_linear_depth_window_step_can_differ_from_gc(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import gbdraw.diagrams.linear.assemble as linear_assemble_module

    record = _make_record()
    captured: dict[str, int] = {}
    real_depth_df = linear_assemble_module.build_depth_df

    def capturing_depth_df(record_arg, table_arg, window_arg, step_arg, **kwargs):
        captured["window"] = int(window_arg)
        captured["step"] = int(step_arg)
        return real_depth_df(record_arg, table_arg, window_arg, step_arg, **kwargs)

    monkeypatch.setattr(linear_assemble_module, "build_depth_df", capturing_depth_df)

    assemble_linear_diagram_from_records(
        [record],
        legend="none",
        depth_table=_depth_table("rec1"),
        config_overrides={"show_gc": True, "show_skew": True},
        window=20,
        step=20,
        depth_window=10,
        depth_step=5,
    )

    assert captured == {"window": 10, "step": 5}


@pytest.mark.linear
def test_linear_depth_window_step_defaults_to_tenth_of_gc(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import gbdraw.diagrams.linear.assemble as linear_assemble_module

    record = _make_record()
    captured: dict[str, int] = {}
    real_depth_df = linear_assemble_module.build_depth_df

    def capturing_depth_df(record_arg, table_arg, window_arg, step_arg, **kwargs):
        captured["window"] = int(window_arg)
        captured["step"] = int(step_arg)
        return real_depth_df(record_arg, table_arg, window_arg, step_arg, **kwargs)

    monkeypatch.setattr(linear_assemble_module, "build_depth_df", capturing_depth_df)

    assemble_linear_diagram_from_records(
        [record],
        legend="none",
        depth_table=_depth_table("rec1"),
        config_overrides={"show_gc": True, "show_skew": True},
        window=100,
        step=100,
    )

    assert captured == {"window": 100, "step": 10}


@pytest.mark.circular
def test_circular_depth_window_step_can_differ_from_gc(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import gbdraw.api.diagram as diagram_module

    record = _make_record()
    captured: dict[str, int] = {}
    real_depth_df = diagram_module.build_depth_df

    def capturing_depth_df(record_arg, table_arg, window_arg, step_arg, **kwargs):
        captured["window"] = int(window_arg)
        captured["step"] = int(step_arg)
        return real_depth_df(record_arg, table_arg, window_arg, step_arg, **kwargs)

    monkeypatch.setattr(diagram_module, "build_depth_df", capturing_depth_df)

    assemble_circular_diagram_from_record(
        record,
        legend="none",
        depth_table=_depth_table("rec1"),
        config_overrides={"show_gc": True, "show_skew": True},
        window=20,
        step=20,
        depth_window=10,
        depth_step=5,
    )

    assert captured == {"window": 10, "step": 5}


@pytest.mark.circular
def test_circular_depth_window_step_defaults_to_tenth_of_gc(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import gbdraw.api.diagram as diagram_module

    record = _make_record()
    captured: dict[str, int] = {}
    real_depth_df = diagram_module.build_depth_df

    def capturing_depth_df(record_arg, table_arg, window_arg, step_arg, **kwargs):
        captured["window"] = int(window_arg)
        captured["step"] = int(step_arg)
        return real_depth_df(record_arg, table_arg, window_arg, step_arg, **kwargs)

    monkeypatch.setattr(diagram_module, "build_depth_df", capturing_depth_df)

    assemble_circular_diagram_from_record(
        record,
        legend="none",
        depth_table=_depth_table("rec1"),
        config_overrides={"show_gc": True, "show_skew": True},
        window=100,
        step=100,
    )

    assert captured == {"window": 100, "step": 10}


def test_circular_cli_depth_options_forward_to_api(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    record = _make_record()
    captured: dict[str, object] = {}
    depth_file = tmp_path / "depth.tsv"
    depth_file.write_text("rec1\t1\t10\n", encoding="utf-8")

    monkeypatch.setattr(circular_cli_module, "load_gbks", lambda paths, mode: [record])
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
            "--depth",
            str(depth_file),
            "--depth_width",
            "22",
            "--depth_window",
            "10",
            "--depth_step",
            "5",
            "--share_depth_axis",
            "--hide_depth_axis",
            "--depth_tick_interval",
            "4",
            "--depth_small_tick_interval",
            "2",
            "--depth_tick_font_size",
            "9",
            "--hide_depth_ticks",
            "--no_depth_log_scale",
            "--format",
            "svg",
            "-o",
            str(tmp_path / "out"),
        ]
    )

    assert captured["depth_file"] == str(depth_file)
    assert captured["depth_window"] == 10
    assert captured["depth_step"] == 5
    assert captured["track_specs"] == ["depth@w=22px"]
    cfg = captured["cfg"]
    assert cfg.objects.depth.show_axis is False
    assert cfg.objects.depth.show_ticks is False
    assert cfg.objects.depth.tick_interval == pytest.approx(4)
    assert cfg.objects.depth.large_tick_interval == pytest.approx(4)
    assert cfg.objects.depth.small_tick_interval == pytest.approx(2)
    assert cfg.objects.depth.tick_font_size == pytest.approx(9)
    assert cfg.objects.depth.normalize is False
    assert cfg.objects.depth.share_axis is True


def test_linear_cli_depth_options_forward_to_api(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    record = _make_record()
    captured: dict[str, object] = {}
    depth_file = tmp_path / "depth.tsv"
    depth_file.write_text("rec1\t1\t10\n", encoding="utf-8")

    monkeypatch.setattr(linear_cli_module, "load_gbks", lambda *args, **kwargs: [record])
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
            "--depth",
            str(depth_file),
            "--depth_height",
            "24",
            "--depth_window",
            "10",
            "--depth_step",
            "5",
            "--share_depth_axis",
            "--hide_depth_axis",
            "--depth_large_tick_interval",
            "4",
            "--depth_small_tick_interval",
            "2",
            "--depth_tick_font_size",
            "9",
            "--hide_depth_ticks",
            "--no_depth_log_scale",
            "--format",
            "svg",
            "-o",
            str(tmp_path / "out"),
        ]
    )

    assert captured["depth_files"] == [str(depth_file)]
    assert captured["depth_window"] == 10
    assert captured["depth_step"] == 5
    cfg = captured["cfg"]
    assert cfg.objects.depth.show_axis is False
    assert cfg.objects.depth.show_ticks is False
    assert cfg.objects.depth.tick_interval == pytest.approx(4)
    assert cfg.objects.depth.large_tick_interval == pytest.approx(4)
    assert cfg.objects.depth.small_tick_interval == pytest.approx(2)
    assert cfg.objects.depth.tick_font_size == pytest.approx(9)
    assert cfg.objects.depth.normalize is False
    assert cfg.objects.depth.share_axis is True


@pytest.mark.linear
def test_linear_depth_dataframe_is_computed_once_per_record(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import gbdraw.diagrams.linear.assemble as linear_assemble_module

    record = _make_record()
    call_count = 0
    real_depth_df = linear_assemble_module.build_depth_df

    def counting_depth_df(*args, **kwargs):
        nonlocal call_count
        call_count += 1
        return real_depth_df(*args, **kwargs)

    monkeypatch.setattr(linear_assemble_module, "build_depth_df", counting_depth_df)

    canvas = assemble_linear_diagram_from_records(
        [record],
        legend="none",
        depth_table=_depth_table("rec1"),
        config_overrides={"show_gc": True, "show_skew": True},
        window=10,
        step=10,
        depth_window=10,
        depth_step=10,
    )

    assert call_count == 1
    assert canvas.tostring().count("L") < 200
