from __future__ import annotations

import re
from pathlib import Path

import pandas as pd
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord
from svgwrite import Drawing

import gbdraw.circular as circular_cli_module
import gbdraw.linear as linear_cli_module
from gbdraw.analysis.depth import clear_depth_tsv_cache, depth_df, read_depth_tsv
from gbdraw.analysis.depth_tracks import (
    DepthTrackSpec,
    build_depth_track_dataframes,
    depth_track_count,
    depth_track_heights,
    index_depth_track_row,
    normalize_depth_tracks,
    representative_depth_tracks,
)
from gbdraw.api.diagram import (
    assemble_circular_diagram_from_records,
    assemble_circular_diagram_from_record,
    assemble_linear_diagram_from_records,
)
from gbdraw.config.toml import load_config_toml
from gbdraw.configurators import DepthConfigurator
from gbdraw.exceptions import ValidationError
from gbdraw.io.comparisons import COMPARISON_COLUMNS
from gbdraw.linear_comparison import LinearComparison


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


def test_sparse_depth_normalization_preserves_logical_indices_count_and_heights() -> None:
    records = [_make_record("rec1", length=40), _make_record("rec2", length=40)]

    normalized = normalize_depth_tracks(
        records,
        depth_track_tables=[
            [_constant_depth_table("rec1", 10), None],
            [None, _constant_depth_table("rec2", 50)],
        ],
        depth_track_labels=["Sample A", "Sample B"],
        depth_track_heights=[12, 28],
    )

    assert normalized is not None
    assert [[track.track_index for track in row] for row in normalized] == [[0], [1]]
    assert [[track.id for track in row] for row in normalized] == [["depth_1"], ["depth_2"]]
    assert depth_track_count(normalized) == 2
    assert depth_track_heights(normalized) == [12, 28]


@pytest.mark.parametrize(
    ("track_indices", "message"),
    [
        ([-1], "negative track_index=-1"),
        ([0, 0], "duplicate track_index=0"),
    ],
)
def test_depth_track_row_rejects_invalid_logical_indices(
    track_indices: list[int],
    message: str,
) -> None:
    row = [
        DepthTrackSpec(
            id=f"depth_{index}",
            label=f"Depth {index}",
            table=_depth_table(),
            track_index=track_index,
        )
        for index, track_index in enumerate(track_indices)
    ]

    with pytest.raises(ValidationError, match=message):
        index_depth_track_row(row, record_index=0)


def test_sparse_depth_shared_axis_and_representatives_are_per_logical_track() -> None:
    records = [_make_record("rec1", length=40), _make_record("rec2", length=40)]
    normalized = normalize_depth_tracks(
        records,
        depth_track_tables=[
            [_constant_depth_table("rec1", 10), None],
            [None, _constant_depth_table("rec2", 50)],
        ],
        depth_track_labels=["Sample A", "Sample B"],
        depth_track_colors=["#112233", "#445566"],
    )
    assert normalized is not None
    config_dict = load_config_toml("gbdraw.data", "config.toml")
    depth_config = DepthConfigurator(
        10,
        10,
        config_dict,
        normalize=False,
        share_axis=True,
    )

    depth_data = build_depth_track_dataframes(
        records,
        normalized,
        base_config=depth_config,
        depth_df_builder=lambda _record, table, _window, _step, **_kwargs: table,
    )
    representatives = representative_depth_tracks(depth_data)

    assert [track.track_index for track in representatives] == [0, 1]
    assert [track.label for track in representatives] == ["Sample A", "Sample B"]
    assert [track.config.fill_color for track in representatives] == ["#112233", "#445566"]
    assert depth_data[0][0].config.max_depth == pytest.approx(10)
    assert depth_data[1][0].config.max_depth == pytest.approx(50)


def test_depth_normalization_rejects_globally_empty_logical_column() -> None:
    records = [_make_record("rec1", length=40), _make_record("rec2", length=40)]

    with pytest.raises(ValidationError, match="Depth track 2 is empty"):
        normalize_depth_tracks(
            records,
            depth_track_tables=[
                [_constant_depth_table("rec1", 10), None],
                [_constant_depth_table("rec2", 20), None],
            ],
        )


def _write_depth_file(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def _svg_group_translate_y(svg: str, group_id: str) -> float:
    match = re.search(rf'<g id="{re.escape(group_id)}" transform="translate\([^,]+,([^)]+)\)"', svg)
    assert match is not None
    return float(match.group(1))


def _svg_group(svg: str, group_id: str) -> str:
    match = re.search(rf'<g id="{re.escape(group_id)}".*?</g>', svg, flags=re.DOTALL)
    assert match is not None
    return match.group(0)


def _svg_line_y_span(svg_fragment: str) -> float:
    y_values = [
        float(value)
        for value in re.findall(r'\by[12]="([^"]+)"', svg_fragment)
    ]
    assert y_values
    return max(y_values) - min(y_values)


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


def test_read_depth_tsv_uses_cache_until_file_changes(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
    request: pytest.FixtureRequest,
) -> None:
    clear_depth_tsv_cache()
    request.addfinalizer(clear_depth_tsv_cache)
    depth_file = _write_depth_file(tmp_path / "depth.tsv", "rec1\t1\t10\nrec1\t2\t12\n")
    real_read_csv = pd.read_csv
    read_count = 0

    def counting_read_csv(*args, **kwargs):
        nonlocal read_count
        read_count += 1
        return real_read_csv(*args, **kwargs)

    monkeypatch.setattr(pd, "read_csv", counting_read_csv)

    first = read_depth_tsv(str(depth_file))
    first.loc[0, "depth"] = 999.0
    second = read_depth_tsv(str(depth_file))

    assert read_count == 1
    assert second["depth"].tolist() == [10.0, 12.0]

    depth_file.write_text("rec1\t1\t11\nrec1\t2\t12\nrec1\t3\t13\n", encoding="utf-8")
    refreshed = read_depth_tsv(str(depth_file))

    assert read_count == 2
    assert refreshed["depth"].tolist() == [11.0, 12.0, 13.0]


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

    canvas = assemble_linear_diagram_from_records(
        [_make_record()],
        legend="none",
        depth_table=_depth_table("rec1"),
        config_overrides={"show_gc": True, "show_skew": True},
    )
    svg = canvas.tostring()
    depth_y = _svg_group_translate_y(svg, "depth")
    gc_y = _svg_group_translate_y(svg, "gc_content")
    skew_y = _svg_group_translate_y(svg, "gc_skew")
    geometry = canvas._gbdraw_track_slot_geometry["records"][0]
    slots = {slot["slotId"]: slot for slot in geometry["slots"]}

    assert (
        slots["gc_content"]["reserveBand"]["topPx"]
        - slots["depth"]["reserveBand"]["bottomPx"]
    ) == pytest.approx(
        cfg.canvas.linear.depth_padding
    )
    assert slots["gc_skew"]["reserveBand"]["topPx"] == pytest.approx(
        slots["gc_content"]["reserveBand"]["bottomPx"]
    )
    assert depth_y == pytest.approx(geometry["axisYpx"] + slots["depth"]["resolvedOriginPx"])
    assert gc_y == pytest.approx(geometry["axisYpx"] + slots["gc_content"]["resolvedOriginPx"])
    assert skew_y == pytest.approx(geometry["axisYpx"] + slots["gc_skew"]["resolvedOriginPx"])


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
        cfg=None,
    ):
        assert track_width_override is not None
        assert norm_factor_override is not None
        captured["depth_width"] = float(track_width_override)
        captured["depth_center"] = float(norm_factor_override) * float(canvas_config.radius)
        return canvas

    monkeypatch.setattr(
        circular_assemble_module,
        "add_depth_group_on_canvas",
        fake_add_depth_group_on_canvas,
    )
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
    default_gc_width = base_radius * track_ratio * float(cfg.canvas.circular.track_ratio_factors[length_param][1])
    default_skew_width = base_radius * track_ratio * float(cfg.canvas.circular.track_ratio_factors[length_param][2])
    default_depth_width = default_gc_width * 0.5

    assert captured["gc_width"] <= default_gc_width + 1e-9
    assert captured["skew_width"] <= default_skew_width + 1e-9
    assert captured["depth_width"] <= default_depth_width + 1e-9

    depth_inner = captured["depth_center"] - (0.5 * captured["depth_width"])
    depth_outer = captured["depth_center"] + (0.5 * captured["depth_width"])
    gc_inner = captured["gc_center"] - (0.5 * captured["gc_width"])
    actual_gc_outer = captured["gc_center"] + (0.5 * captured["gc_width"])
    skew_inner = captured["skew_center"] - (0.5 * captured["skew_width"])
    skew_outer = captured["skew_center"] + (0.5 * captured["skew_width"])
    assert skew_inner >= 0.0
    assert depth_inner < depth_outer
    assert actual_gc_outer < depth_inner
    assert skew_outer < gc_inner


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
        group_id=None,
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
        group_id=None,
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
        circular_track_slots=[
            "gc_content:dinucleotide_content@r=0.44,w=12px",
            "gc_skew:dinucleotide_skew@r=0.24,w=10px",
        ],
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
def test_linear_multiple_depth_tracks_have_unique_ids() -> None:
    records = [_make_record("rec1", length=40), _make_record("rec2", length=40)]
    canvas = assemble_linear_diagram_from_records(
        records,
        legend="none",
        depth_track_tables=[
            [
                _constant_depth_table("rec1", 10, length=40),
                _constant_depth_table("rec1", 20, length=40),
            ],
            [
                _constant_depth_table("rec2", 30, length=40),
                _constant_depth_table("rec2", 40, length=40),
            ],
        ],
        depth_track_labels=["Sample A", "Sample B"],
        depth_track_colors=["#111111", "#222222"],
        config_overrides={"show_gc": False, "show_skew": False},
        window=10,
        step=10,
        depth_window=10,
        depth_step=10,
    )
    svg = canvas.tostring()

    for group_id in (
        "depth_1_record_1",
        "depth_2_record_1",
        "depth_1_record_2",
        "depth_2_record_2",
    ):
        assert svg.count(f'id="{group_id}"') == 1
        assert svg.count(f'id="{group_id}_axis"') == 1
    assert 'fill="#111111"' in svg
    assert 'fill="#222222"' in svg


@pytest.mark.linear
def test_linear_multiple_depth_track_heights_are_per_track() -> None:
    record = _make_record("rec1", length=40)
    canvas = assemble_linear_diagram_from_records(
        [record],
        legend="none",
        depth_track_tables=[
            [
                _constant_depth_table("rec1", 10, length=40),
                _constant_depth_table("rec1", 40, length=40),
            ]
        ],
        depth_track_heights=[12, 28],
        config_overrides={"show_gc": False, "show_skew": False},
        window=10,
        step=10,
        depth_window=10,
        depth_step=10,
    )
    svg = canvas.tostring()

    assert _svg_line_y_span(_svg_group(svg, "depth_1_axis")) == pytest.approx(12)
    assert _svg_line_y_span(_svg_group(svg, "depth_2_axis")) == pytest.approx(28)


@pytest.mark.linear
def test_linear_custom_depth_slot_uses_track_height_when_slot_height_is_auto() -> None:
    record = _make_record("rec1", length=40)
    svg = assemble_linear_diagram_from_records(
        [record],
        legend="none",
        depth_track_tables=[
            [
                _constant_depth_table("rec1", 10, length=40),
                _constant_depth_table("rec1", 40, length=40),
            ]
        ],
        depth_track_heights=[12, 28],
        linear_track_slots=[
            "features:features@side=overlay",
            "depth_a:depth@track_index=0",
            "depth_b:depth@track_index=1",
        ],
        config_overrides={"show_gc": False, "show_skew": False},
        window=10,
        step=10,
        depth_window=10,
        depth_step=10,
    ).tostring()

    assert _svg_line_y_span(_svg_group(svg, "depth_a_axis")) == pytest.approx(12)
    assert _svg_line_y_span(_svg_group(svg, "depth_b_axis")) == pytest.approx(28)


@pytest.mark.linear
def test_linear_depth_track_tick_options_are_per_track() -> None:
    record = _make_record("rec1", length=40)
    svg = assemble_linear_diagram_from_records(
        [record],
        legend="none",
        depth_track_tables=[
            [
                _constant_depth_table("rec1", 10, length=40),
                _constant_depth_table("rec1", 40, length=40),
            ]
        ],
        depth_track_large_tick_intervals=[5, 20],
        depth_track_small_tick_intervals=["auto", 10],
        depth_track_tick_font_sizes=[8, 14],
        config_overrides={"show_gc": False, "show_skew": False},
        window=10,
        step=10,
        depth_window=10,
        depth_step=10,
    ).tostring()

    depth_1_axis = _svg_group(svg, "depth_1_axis")
    depth_2_axis = _svg_group(svg, "depth_2_axis")

    assert ">5x<" in depth_1_axis
    assert ">20x<" not in depth_1_axis
    assert ">20x<" in depth_2_axis
    assert ">40x<" in depth_2_axis


@pytest.mark.circular
def test_circular_multiple_depth_tracks_render_distinct_rings() -> None:
    record = _make_record("rec1", length=40)
    svg = assemble_circular_diagram_from_record(
        record,
        legend="none",
        depth_track_tables=[
            [
                _constant_depth_table("rec1", 10, length=40),
                _constant_depth_table("rec1", 40, length=40),
            ]
        ],
        depth_track_labels=["Sample A", "Sample B"],
        config_overrides={"show_gc": False, "show_skew": False},
        window=10,
        step=10,
        depth_window=10,
        depth_step=10,
    ).tostring()

    assert svg.count('id="depth_1"') == 1
    assert svg.count('id="depth_1_axis"') == 1
    assert svg.count('id="depth_2"') == 1
    assert svg.count('id="depth_2_axis"') == 1


@pytest.mark.circular
def test_circular_custom_depth_slot_selects_track_index() -> None:
    record = _make_record("rec1", length=40)
    svg = assemble_circular_diagram_from_record(
        record,
        legend="none",
        depth_track_tables=[
            [
                _constant_depth_table("rec1", 10, length=40),
                _constant_depth_table("rec1", 40, length=40),
            ]
        ],
        circular_track_slots=[
            "features:features",
            "ticks:ticks",
            "selected_depth:depth@track_index=1,label=Sample B",
        ],
        config_overrides={"show_gc": False, "show_skew": False},
        window=10,
        step=10,
        depth_window=10,
        depth_step=10,
    ).tostring()

    assert 'id="selected_depth"' in svg
    assert 'id="selected_depth_axis"' in svg
    assert 'id="depth_1"' not in svg


@pytest.mark.linear
def test_depth_share_axis_is_per_logical_track() -> None:
    records = [_make_record("rec1", length=40), _make_record("rec2", length=40)]
    svg = assemble_linear_diagram_from_records(
        records,
        legend="none",
        depth_track_tables=[
            [
                _constant_depth_table("rec1", 10, length=40),
                _constant_depth_table("rec1", 5, length=40),
            ],
            [
                _constant_depth_table("rec2", 100, length=40),
                _constant_depth_table("rec2", 50, length=40),
            ],
        ],
        config_overrides={
            "show_gc": False,
            "show_skew": False,
            "depth_share_axis": True,
        },
        window=10,
        step=10,
        depth_window=10,
        depth_step=10,
    ).tostring()

    assert ">100x<" in svg
    assert ">50x<" in svg


@pytest.mark.linear
def test_linear_depth_track_partial_record_rows_are_not_shared() -> None:
    records = [
        _make_record("rec1", length=40),
        _make_record("rec2", length=40),
        _make_record("rec3", length=40),
    ]
    svg = assemble_linear_diagram_from_records(
        records,
        legend="none",
        depth_track_tables=[
            [_constant_depth_table("rec1", 10, length=40)],
            [_constant_depth_table("rec2", 20, length=40)],
            [],
        ],
        config_overrides={"show_gc": False, "show_skew": False},
        window=10,
        step=10,
        depth_window=10,
        depth_step=10,
    ).tostring()

    assert svg.count('id="depth_record_1"') == 1
    assert svg.count('id="depth_record_2"') == 1
    assert svg.count('id="depth_record_1_axis"') == 1
    assert svg.count('id="depth_record_2_axis"') == 1


@pytest.mark.parametrize(
    ("present_record_index", "depth_side", "feature_side"),
    [
        (0, "above", "below"),
        (1, "below", "above"),
    ],
)
@pytest.mark.linear
def test_linear_custom_depth_slot_skips_leading_or_trailing_missing_record(
    present_record_index: int,
    depth_side: str,
    feature_side: str,
) -> None:
    records = [_make_record("rec1", length=40), _make_record("rec2", length=40)]
    depth_rows: list[list[pd.DataFrame | None]] = [[None], [None]]
    present_record_id = records[present_record_index].id
    depth_rows[present_record_index][0] = _constant_depth_table(present_record_id, 25, length=40)

    canvas = assemble_linear_diagram_from_records(
        records,
        legend="none",
        depth_track_tables=depth_rows,
        linear_track_slots=[
            f"depth:depth@track_index=0,side={depth_side}",
            f"features:features@side={feature_side}",
        ],
        config_overrides={"show_gc": False, "show_skew": False},
        window=10,
        step=10,
        depth_window=10,
        depth_step=10,
    )
    svg = canvas.tostring()

    present_number = present_record_index + 1
    missing_number = 2 if present_number == 1 else 1
    assert svg.count(f'id="depth_record_{present_number}"') == 1
    assert svg.count(f'id="depth_record_{present_number}_axis"') == 1
    assert f'id="depth_record_{missing_number}"' not in svg
    assert f'id="depth_record_{missing_number}_axis"' not in svg
    geometry = canvas._gbdraw_track_slot_geometry["records"]
    depth_slots = [
        next(slot for slot in record_geometry["slots"] if slot["slotId"] == "depth")
        for record_geometry in geometry
    ]
    assert depth_slots[present_record_index]["dataAvailable"] is True
    assert depth_slots[present_record_index]["paintBand"] is not None
    assert depth_slots[1 - present_record_index]["dataAvailable"] is False
    assert depth_slots[1 - present_record_index]["paintBand"] is None
    present_slot = depth_slots[present_record_index]
    missing_slot = depth_slots[1 - present_record_index]
    assert present_slot["reserveBand"]["bottomPx"] > present_slot["reserveBand"]["topPx"]
    assert missing_slot["reserveBand"]["bottomPx"] == pytest.approx(
        missing_slot["reserveBand"]["topPx"]
    )
    present_body = geometry[present_record_index]["recordBodyBand"]
    missing_body = geometry[1 - present_record_index]["recordBodyBand"]
    if depth_side == "above":
        assert missing_body["topPx"] > present_body["topPx"]
    else:
        assert missing_body["bottomPx"] < present_body["bottomPx"]


@pytest.mark.linear
def test_sparse_depth_does_not_separate_pairwise_match_from_missing_record_features() -> None:
    records = [_make_record("rec1", length=40), _make_record("rec2", length=40)]
    comparison = LinearComparison(
        0,
        1,
        pd.DataFrame(
            [["rec1", "rec2", 90.0, 20, 0, 0, 5, 25, 5, 25, 1e-10, 100]],
            columns=COMPARISON_COLUMNS,
        ),
    )
    canvas = assemble_linear_diagram_from_records(
        records,
        legend="none",
        depth_track_tables=[[_constant_depth_table("rec1", 25, length=40)], [None]],
        linear_track_slots=[
            "depth:depth@track_index=0,side=above",
            "features:features@side=above",
        ],
        linear_comparisons=[comparison],
        config_overrides={"show_gc": False, "show_skew": False},
        window=10,
        step=10,
        depth_window=10,
        depth_step=10,
    )

    missing_record = canvas._gbdraw_track_slot_geometry["records"][1]
    slots = {slot["slotId"]: slot for slot in missing_record["slots"]}
    assert slots["depth"]["dataAvailable"] is False
    assert slots["depth"]["reserveBand"]["topPx"] == pytest.approx(
        slots["depth"]["reserveBand"]["bottomPx"]
    )
    assert missing_record["comparisonExclusionBand"]["topPx"] == pytest.approx(
        slots["features"]["paintBand"]["topPx"]
    )


@pytest.mark.linear
def test_linear_diagonal_sparse_depth_binding_survives_slot_order_reversal() -> None:
    records = [_make_record("rec1", length=40), _make_record("rec2", length=40)]
    depth_rows = [
        [_constant_depth_table("rec1", 10, length=40), None],
        [None, _constant_depth_table("rec2", 50, length=40)],
    ]

    def render(depth_slot_ids: list[str]) -> str:
        slots = {
            "depth_a": "depth_a:depth@track_index=0,side=below",
            "depth_b": "depth_b:depth@track_index=1,side=below",
        }
        return assemble_linear_diagram_from_records(
            records,
            legend="none",
            depth_track_tables=depth_rows,
            depth_track_labels=["Sample A", "Sample B"],
            depth_track_colors=["#112233", "#445566"],
            linear_track_slots=[
                *(slots[slot_id] for slot_id in depth_slot_ids),
                "features:features@side=overlay",
            ],
            config_overrides={"show_gc": False, "show_skew": False},
            window=10,
            step=10,
            depth_window=10,
            depth_step=10,
        ).tostring()

    original = render(["depth_a", "depth_b"])
    reversed_order = render(["depth_b", "depth_a"])

    for svg in (original, reversed_order):
        assert 'fill="#112233"' in _svg_group(svg, "depth_a_record_1")
        assert 'fill="#445566"' in _svg_group(svg, "depth_b_record_2")
        assert 'id="depth_a_record_2"' not in svg
        assert 'id="depth_b_record_1"' not in svg
    assert _svg_group_translate_y(original, "depth_a_record_1") == pytest.approx(
        _svg_group_translate_y(reversed_order, "depth_a_record_1")
    )
    assert _svg_group_translate_y(original, "depth_b_record_2") == pytest.approx(
        _svg_group_translate_y(reversed_order, "depth_b_record_2")
    )


@pytest.mark.linear
def test_linear_custom_slot_can_select_only_second_logical_depth_track() -> None:
    record = _make_record("rec1", length=40)

    svg = assemble_linear_diagram_from_records(
        [record],
        legend="none",
        depth_track_tables=[
            [
                _constant_depth_table("rec1", 10, length=40),
                _constant_depth_table("rec1", 50, length=40),
            ]
        ],
        depth_track_colors=["#112233", "#445566"],
        linear_track_slots=[
            "selected_depth:depth@track_index=1,side=below",
            "features:features@side=overlay",
        ],
        config_overrides={"show_gc": False, "show_skew": False},
        window=10,
        step=10,
        depth_window=10,
        depth_step=10,
    ).tostring()

    assert 'fill="#445566"' in _svg_group(svg, "selected_depth")
    assert 'id="selected_depth_axis"' in svg


@pytest.mark.linear
def test_linear_custom_depth_slot_legend_uses_only_selected_sparse_logical_track() -> None:
    records = [_make_record("rec1", length=40), _make_record("rec2", length=40)]
    legend_label = "Selected Sample B"

    svg = assemble_linear_diagram_from_records(
        records,
        legend="right",
        depth_track_tables=[
            [_constant_depth_table("rec1", 10, length=40), None],
            [None, _constant_depth_table("rec2", 50, length=40)],
        ],
        depth_track_labels=["Original Sample A", "Original Sample B"],
        depth_track_colors=["#112233", "#445566"],
        linear_track_slots=[
            f"selected_depth:depth@track_index=1,side=below,legend_label={legend_label}",
            "features:features@side=overlay",
        ],
        config_overrides={"show_gc": False, "show_skew": False},
        window=10,
        step=10,
        depth_window=10,
        depth_step=10,
    ).tostring()

    assert 'id="selected_depth_record_1"' not in svg
    assert 'id="selected_depth_record_1_axis"' not in svg
    assert svg.count('id="selected_depth_record_2"') == 1
    assert svg.count('id="selected_depth_record_2_axis"') == 1
    assert 'fill="#445566"' in _svg_group(svg, "selected_depth_record_2")
    assert "#112233" not in svg
    assert f'data-legend-key="{legend_label}"' in svg
    assert "Original Sample A" not in svg
    assert "Original Sample B" not in svg


@pytest.mark.linear
def test_linear_depth_slot_rejects_globally_out_of_range_logical_index() -> None:
    record = _make_record("rec1", length=40)

    with pytest.raises(
        ValidationError,
        match="Depth slot 'missing_depth' track_index=2.*range 0\\.\\.1",
    ):
        assemble_linear_diagram_from_records(
            [record],
            legend="none",
            depth_track_tables=[
                [
                    _constant_depth_table("rec1", 10, length=40),
                    _constant_depth_table("rec1", 50, length=40),
                ]
            ],
            linear_track_slots=[
                "missing_depth:depth@track_index=2",
                "features:features@side=overlay",
            ],
            config_overrides={"show_gc": False, "show_skew": False},
        )


@pytest.mark.parametrize("present_record_index", [0, 1])
@pytest.mark.circular
def test_circular_multi_record_partial_depth_skips_paint_but_reserves_slot(
    present_record_index: int,
    caplog: pytest.LogCaptureFixture,
) -> None:
    records = [_make_record("rec1", length=40), _make_record("rec2", length=40)]
    depth_rows: list[list[pd.DataFrame | None]] = [[None], [None]]
    present_record_id = records[present_record_index].id
    depth_rows[present_record_index][0] = _constant_depth_table(present_record_id, 25, length=40)

    canvas = assemble_circular_diagram_from_records(
        records,
        legend="none",
        depth_track_tables=depth_rows,
        config_overrides={"show_gc": False, "show_skew": False},
        window=10,
        step=10,
        depth_window=10,
        depth_step=10,
    )
    svg = canvas.tostring()
    geometry = getattr(canvas, "_gbdraw_track_slot_geometry", None)

    assert svg.count('id="depth"') == 1
    present_number = present_record_index + 1
    missing_number = 2 if present_number == 1 else 1
    assert svg.count(f'id="depth_record_{present_number}_axis"') == 1
    assert f'id="depth_record_{missing_number}_axis"' not in svg
    assert len(geometry["records"]) == 2
    assert all(
        any(slot["slotId"] == "depth" for slot in record_geometry["slots"])
        for record_geometry in geometry["records"]
    )
    assert "depth data are unavailable" not in caplog.text


@pytest.mark.circular
def test_circular_multi_record_diagonal_sparse_depth_uses_logical_binding() -> None:
    records = [_make_record("rec1", length=40), _make_record("rec2", length=40)]

    canvas = assemble_circular_diagram_from_records(
        records,
        legend="none",
        depth_track_tables=[
            [_constant_depth_table("rec1", 10, length=40), None],
            [None, _constant_depth_table("rec2", 50, length=40)],
        ],
        depth_track_labels=["Sample A", "Sample B"],
        depth_track_colors=["#112233", "#445566"],
        config_overrides={"show_gc": False, "show_skew": False},
        window=10,
        step=10,
        depth_window=10,
        depth_step=10,
    )
    svg = canvas.tostring()
    geometry = getattr(canvas, "_gbdraw_track_slot_geometry", None)

    assert svg.count('id="depth_1_record_1"') == 1
    assert svg.count('id="depth_1_record_1_axis"') == 1
    assert svg.count('id="depth_2_record_2"') == 1
    assert svg.count('id="depth_2_record_2_axis"') == 1
    assert 'id="depth_1_record_2"' not in svg
    assert 'id="depth_2_record_1"' not in svg
    assert 'fill="#112233"' in _svg_group(svg, "depth_1_record_1")
    assert 'fill="#445566"' in _svg_group(svg, "depth_2_record_2")
    assert all(
        {slot["slotId"] for slot in record_geometry["slots"]} >= {"depth_1", "depth_2"}
        for record_geometry in geometry["records"]
    )


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

    assert captured["depth_table"].to_dict("list") == {
        "reference_name": ["rec1"],
        "position": [1],
        "depth": [10.0],
    }
    assert captured["depth_window"] == 10
    assert captured["depth_step"] == 5
    slots = captured["circular_track_slots"]
    by_id = {slot.id: slot for slot in slots}
    assert by_id["depth"].width.resolve(390.0) == pytest.approx(22.0)
    cfg = captured["cfg"]
    assert cfg.objects.depth.show_axis is False
    assert cfg.objects.depth.show_ticks is False
    assert cfg.objects.depth.tick_interval == pytest.approx(4)
    assert cfg.objects.depth.large_tick_interval == pytest.approx(4)
    assert cfg.objects.depth.small_tick_interval == pytest.approx(2)
    assert cfg.objects.depth.tick_font_size == pytest.approx(9)
    assert cfg.objects.depth.normalize is False
    assert cfg.objects.depth.share_axis is True


def test_circular_cli_reads_depth_file_once_for_multiple_records(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    records = [_make_record("rec1"), _make_record("rec2")]
    depth_file = tmp_path / "depth.tsv"
    depth_file.write_text("rec1\t1\t10\nrec2\t1\t20\n", encoding="utf-8")
    read_count = 0
    captured_tables: list[pd.DataFrame | None] = []
    real_read_depth_tsv = circular_cli_module.read_depth_tsv

    def counting_read_depth_tsv(path: str) -> pd.DataFrame:
        nonlocal read_count
        read_count += 1
        return real_read_depth_tsv(path)

    monkeypatch.setattr(circular_cli_module, "read_depth_tsv", counting_read_depth_tsv)
    monkeypatch.setattr(circular_cli_module, "load_gbks", lambda paths, mode: records)
    monkeypatch.setattr(circular_cli_module, "read_color_table", lambda _path: None)
    monkeypatch.setattr(circular_cli_module, "load_default_colors", lambda _path, _palette: None)
    monkeypatch.setattr(circular_cli_module, "save_figure", lambda canvas, formats: None)

    def fake_assemble(*args, **kwargs):
        captured_tables.append(kwargs.get("depth_table"))
        return Drawing(filename=str(tmp_path / f"dummy_{len(captured_tables)}.svg"))

    monkeypatch.setattr(circular_cli_module, "assemble_circular_diagram_from_record", fake_assemble)

    circular_cli_module.circular_main(
        [
            "--gbk",
            "dummy1.gb",
            "dummy2.gb",
            "--depth",
            str(depth_file),
            "--format",
            "svg",
            "-o",
            str(tmp_path / "out"),
        ]
    )

    assert read_count == 1
    assert len(captured_tables) == 2
    assert captured_tables[0] is captured_tables[1]


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


def test_linear_cli_repeated_depth_track_forwards_record_major_files(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    records = [_make_record("rec1"), _make_record("rec2")]
    captured: dict[str, object] = {}

    monkeypatch.setattr(linear_cli_module, "load_gbks", lambda *args, **kwargs: records)
    monkeypatch.setattr(linear_cli_module, "read_color_table", lambda _path: None)
    monkeypatch.setattr(linear_cli_module, "load_default_colors", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(linear_cli_module, "read_feature_visibility_file", lambda _path: None)
    monkeypatch.setattr(linear_cli_module, "save_figure", lambda canvas, formats: None)

    def fake_assemble(*args, **kwargs):
        captured.update(kwargs)
        return Drawing(filename=str(tmp_path / "dummy.svg"))

    monkeypatch.setattr(linear_cli_module, "assemble_linear_diagram_from_records", fake_assemble)

    paths = [str(tmp_path / name) for name in ("r1_a.tsv", "r2_a.tsv", "r1_b.tsv", "r2_b.tsv")]
    linear_cli_module.linear_main(
        [
            "--gbk",
            "rec1.gb",
            "rec2.gb",
            "--depth_track",
            paths[0],
            paths[1],
            "--depth_track",
            paths[2],
            paths[3],
            "--depth_track_label",
            "A",
            "B",
            "--depth_track_color",
            "#111111",
            "#222222",
            "--depth_track_height",
            "12",
            "28",
            "--depth_track_large_tick_interval",
            "10",
            "20",
            "--depth_track_small_tick_interval",
            "auto",
            "5",
            "--depth_track_tick_font_size",
            "8",
            "12",
            "--format",
            "svg",
            "-o",
            str(tmp_path / "out"),
        ]
    )

    assert captured["depth_files"] is None
    assert captured["depth_track_files"] == [
        [paths[0], paths[2]],
        [paths[1], paths[3]],
    ]
    assert captured["depth_track_labels"] == ["A", "B"]
    assert captured["depth_track_colors"] == ["#111111", "#222222"]
    assert captured["depth_track_heights"] == ["12", "28"]
    assert captured["depth_track_large_tick_intervals"] == ["10", "20"]
    assert captured["depth_track_small_tick_intervals"] == ["auto", "5"]
    assert captured["depth_track_tick_font_sizes"] == ["8", "12"]


def test_linear_cli_depth_track_placeholders_keep_record_slots(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    records = [_make_record("rec1"), _make_record("rec2")]
    captured: dict[str, object] = {}

    monkeypatch.setattr(linear_cli_module, "load_gbks", lambda *args, **kwargs: records)
    monkeypatch.setattr(linear_cli_module, "read_color_table", lambda _path: None)
    monkeypatch.setattr(linear_cli_module, "load_default_colors", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(linear_cli_module, "read_feature_visibility_file", lambda _path: None)
    monkeypatch.setattr(linear_cli_module, "save_figure", lambda canvas, formats: None)

    def fake_assemble(*args, **kwargs):
        captured.update(kwargs)
        return Drawing(filename=str(tmp_path / "dummy.svg"))

    monkeypatch.setattr(linear_cli_module, "assemble_linear_diagram_from_records", fake_assemble)

    paths = [str(tmp_path / name) for name in ("r1_a.tsv", "r2_b.tsv")]
    linear_cli_module.linear_main(
        [
            "--gbk",
            "rec1.gb",
            "rec2.gb",
            "--depth_track",
            paths[0],
            "",
            "--depth_track",
            "",
            paths[1],
            "--format",
            "svg",
            "-o",
            str(tmp_path / "out"),
        ]
    )

    assert captured["depth_track_files"] == [
        [paths[0], None],
        [None, paths[1]],
    ]


@pytest.mark.linear
def test_linear_cli_sparse_depth_tracks_render_end_to_end(tmp_path: Path) -> None:
    records = [_make_record("rec1", length=40), _make_record("rec2", length=40)]
    record_paths = [tmp_path / "rec1.gb", tmp_path / "rec2.gb"]
    for record, path in zip(records, record_paths):
        SeqIO.write(record, path, "genbank")

    depth_paths = [tmp_path / "rec1.tsv", tmp_path / "rec2.tsv"]
    for record, path, depth in zip(records, depth_paths, (10, 50)):
        _write_depth_file(
            path,
            "".join(f"{record.id}\t{position}\t{depth}\n" for position in range(1, 41)),
        )

    output_prefix = tmp_path / "sparse-depth"
    linear_cli_module.linear_main(
        [
            "--gbk",
            *(str(path) for path in record_paths),
            "--depth_track",
            str(depth_paths[0]),
            "",
            "--depth_track",
            "",
            str(depth_paths[1]),
            "--depth_track_label",
            "Sample A",
            "Sample B",
            "--depth_track_color",
            "#112233",
            "#445566",
            "--linear_track_slot",
            "depth_a:depth@track_index=0,side=above",
            "--linear_track_slot",
            "depth_b:depth@track_index=1,side=above",
            "--linear_track_slot",
            "features:features@side=below",
            "--linear_track_axis_index",
            "2",
            "--format",
            "svg",
            "--legend",
            "none",
            "-o",
            str(output_prefix),
        ]
    )

    svg = output_prefix.with_suffix(".svg").read_text(encoding="utf-8")
    assert 'fill="#112233"' in _svg_group(svg, "depth_a_record_1")
    assert 'fill="#445566"' in _svg_group(svg, "depth_b_record_2")
    assert 'id="depth_a_record_2"' not in svg
    assert 'id="depth_b_record_1"' not in svg


def test_circular_cli_repeated_depth_track_forwards_record_major_files(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    records = [_make_record("rec1"), _make_record("rec2")]
    captured: dict[str, object] = {}

    monkeypatch.setattr(circular_cli_module, "load_gbks", lambda paths, mode: records)
    monkeypatch.setattr(circular_cli_module, "read_color_table", lambda _path: None)
    monkeypatch.setattr(circular_cli_module, "load_default_colors", lambda _path, _palette: None)
    monkeypatch.setattr(circular_cli_module, "save_figure", lambda canvas, formats: None)

    def fake_assemble(*args, **kwargs):
        captured.update(kwargs)
        return Drawing(filename=str(tmp_path / "dummy.svg"))

    monkeypatch.setattr(circular_cli_module, "assemble_circular_diagram_from_records", fake_assemble)

    paths = [str(tmp_path / name) for name in ("r1_a.tsv", "r2_a.tsv", "r1_b.tsv", "r2_b.tsv")]
    circular_cli_module.circular_main(
        [
            "--gbk",
            "rec1.gb",
            "rec2.gb",
            "--multi_record_canvas",
            "--depth_track",
            paths[0],
            paths[1],
            "--depth_track",
            paths[2],
            paths[3],
            "--depth_track_large_tick_interval",
            "10",
            "20",
            "--depth_track_small_tick_interval",
            "auto",
            "5",
            "--depth_track_tick_font_size",
            "8",
            "12",
            "--format",
            "svg",
            "-o",
            str(tmp_path / "out"),
        ]
    )

    assert captured["depth_table"] is None
    assert captured["depth_track_files"] == [
        [paths[0], paths[2]],
        [paths[1], paths[3]],
    ]
    assert captured["depth_track_large_tick_intervals"] == ["10", "20"]
    assert captured["depth_track_small_tick_intervals"] == ["auto", "5"]
    assert captured["depth_track_tick_font_sizes"] == ["8", "12"]


@pytest.mark.circular
def test_circular_cli_sparse_depth_tracks_render_separate_outputs(tmp_path: Path) -> None:
    records = [_make_record("rec1", length=40), _make_record("rec2", length=40)]
    record_paths = [tmp_path / "rec1.gb", tmp_path / "rec2.gb"]
    for record, path in zip(records, record_paths):
        SeqIO.write(record, path, "genbank")

    depth_paths = [tmp_path / "rec1.tsv", tmp_path / "rec2.tsv"]
    for record, path, depth in zip(records, depth_paths, (10, 50)):
        _write_depth_file(
            path,
            "".join(f"{record.id}\t{position}\t{depth}\n" for position in range(1, 41)),
        )

    output_prefix = tmp_path / "sparse-circular"
    circular_cli_module.circular_main(
        [
            "--gbk",
            *(str(path) for path in record_paths),
            "--depth_track",
            str(depth_paths[0]),
            "",
            "--depth_track",
            "",
            str(depth_paths[1]),
            "--depth_track_color",
            "#112233",
            "#445566",
            "--format",
            "svg",
            "--legend",
            "none",
            "-o",
            str(output_prefix),
        ]
    )

    first_svg = Path(f"{output_prefix}_1.svg").read_text(encoding="utf-8")
    second_svg = Path(f"{output_prefix}_2.svg").read_text(encoding="utf-8")
    assert 'fill="#112233"' in _svg_group(first_svg, "depth_1")
    assert 'id="depth_2"' not in first_svg
    assert 'fill="#445566"' in _svg_group(second_svg, "depth_2")
    assert 'id="depth_1"' not in second_svg


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
