from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from svgwrite import Drawing

import gbdraw.circular as circular_cli_module
import gbdraw.linear as linear_cli_module
from gbdraw.analysis.depth_tracks import (
    DepthTrackSpec,
    build_depth_track_dataframes,
)
from gbdraw.cli_utils.table_adapters import (
    load_depth_track_table,
    load_input_table_records,
    load_track_table_slots,
)
from gbdraw.cli_utils.tables import read_headered_tsv_table
from gbdraw.config.models import GbdrawConfig
from gbdraw.config.toml import load_config_toml
from gbdraw.configurators import DepthConfigurator
from gbdraw.exceptions import ValidationError


def _record(record_id: str, length: int = 40) -> SeqRecord:
    record = SeqRecord(Seq("A" * length), id=record_id, name=record_id)
    record.annotations["molecule_type"] = "DNA"
    return record


def _depth_table(reference: str, depth: float, length: int = 40) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "reference_name": [reference] * length,
            "position": list(range(1, length + 1)),
            "depth": [depth] * length,
        }
    )


def _write_depth_file(path: Path, reference: str, depth: float) -> str:
    path.write_text(
        "\n".join(f"{reference}\t{position}\t{depth}" for position in range(1, 41)) + "\n",
        encoding="utf-8",
    )
    return str(path)


def test_headered_table_keeps_hash_record_selector(tmp_path: Path) -> None:
    table = tmp_path / "depth_table.tsv"
    table.write_text(
        "# comment\nrecord_id\ttrack_id\tfile\n#1\tdepth_A\tdepth.tsv\n",
        encoding="utf-8",
    )

    rows = read_headered_tsv_table(
        str(table),
        required=("record_id", "track_id", "file"),
        optional=(),
        table_name="depth_track_table",
    )

    assert len(rows) == 1
    assert rows[0].cell("record_id") == "#1"


def test_input_table_row_local_record_selector_and_input_id(tmp_path: Path) -> None:
    table = tmp_path / "inputs.tsv"
    table.write_text(
        "input_id\tinput_type\tgbk\trecord_id\tlabel\n"
        "sample\tgbk\tdummy.gb\t#2\tSample label\n",
        encoding="utf-8",
    )
    source_records = [_record("recA"), _record("recB")]

    def fake_load_gbk(*_args, **_kwargs):
        return source_records

    records = load_input_table_records(
        str(table),
        mode="linear",
        load_gbk_records=fake_load_gbk,
        load_gff_records=lambda *_args, **_kwargs: [],
    )

    assert [record.id for record in records] == ["recB"]
    assert records[0].annotations["gbdraw_input_id"] == "sample"
    assert records[0].annotations["gbdraw_record_label"] == "Sample label"


def test_depth_track_table_sparse_rows_preserve_track_ids(tmp_path: Path) -> None:
    records = [_record("rec1"), _record("rec2")]
    records[0].annotations["gbdraw_input_id"] = "ref"
    records[1].annotations["gbdraw_input_id"] = "sample"
    depth_a = _write_depth_file(tmp_path / "depth_a.tsv", "rec1", 10)
    depth_a_sample = _write_depth_file(tmp_path / "depth_a_sample.tsv", "rec2", 20)
    depth_b_sample = _write_depth_file(tmp_path / "depth_b_sample.tsv", "rec2", 50)
    table = tmp_path / "depth_table.tsv"
    table.write_text(
        "record_id\ttrack_id\tfile\ttrack_label\ttrack_color\torder\n"
        f"*\tdepth_A\t{Path(depth_a).name}\tSample A\t#4A90E2\t1\n"
        f"input:sample\tdepth_A\t{Path(depth_a_sample).name}\tSample A\t#4A90E2\t1\n"
        f"input:sample\tdepth_B\t{Path(depth_b_sample).name}\tSample B\t#E45756\t2\n",
        encoding="utf-8",
    )

    result = load_depth_track_table(str(table), mode="linear", records=records)

    assert result.track_ids == ("depth_A", "depth_B")
    assert [[spec.id for spec in row] for row in result.record_depth_tracks] == [
        ["depth_A"],
        ["depth_A", "depth_B"],
    ]
    assert result.metadata_by_track_id["depth_B"].fill_color == "#E45756"


def test_depth_track_table_rejects_mode_specific_size_columns(tmp_path: Path) -> None:
    record = _record("rec1")
    depth_path = _write_depth_file(tmp_path / "depth.tsv", "rec1", 10)
    table = tmp_path / "depth_table.tsv"
    table.write_text(
        "record_id\ttrack_id\tfile\ttrack_height\n"
        f"#1\tdepth_A\t{Path(depth_path).name}\t12\n",
        encoding="utf-8",
    )

    with pytest.raises(ValidationError, match="track_height is only valid in linear mode"):
        load_depth_track_table(str(table), mode="circular", records=[record])


def test_track_table_depth_track_id_converts_to_slot_param(tmp_path: Path) -> None:
    table = tmp_path / "track_table.tsv"
    table.write_text(
        "slot_id\trenderer\torder\tside\ttrack_id\theight\tenabled\n"
        "features\tfeatures\t1\toverlay\t\t\ttrue\n"
        "depth_B\tdepth\t2\tbelow\tdepth_B\t18px\ttrue\n",
        encoding="utf-8",
    )

    result = load_track_table_slots(
        str(table),
        mode="linear",
        depth_track_ids=("depth_A", "depth_B"),
    )

    depth_slot = result.slots[1]
    assert depth_slot.params["track_id"] == "depth_B"
    assert depth_slot.params["track_index"] == 1
    assert depth_slot.height is not None
    assert depth_slot.height.resolve(1) == pytest.approx(18)


def test_share_depth_axis_groups_sparse_tracks_by_id() -> None:
    records = [_record("rec1"), _record("rec2")]
    record_depth_tracks = [
        [
            DepthTrackSpec(
                id="depth_B",
                label="Depth B",
                table=_depth_table("rec1", 100),
            )
        ],
        [
            DepthTrackSpec(
                id="depth_A",
                label="Depth A",
                table=_depth_table("rec2", 5),
            ),
            DepthTrackSpec(
                id="depth_B",
                label="Depth B",
                table=_depth_table("rec2", 50),
            ),
        ],
    ]
    config_dict = load_config_toml("gbdraw.data", "config.toml")
    cfg = GbdrawConfig.from_dict(config_dict)
    depth_config = DepthConfigurator(
        window=10,
        step=10,
        config_dict=config_dict,
        cfg=cfg,
        share_axis=True,
    )

    depth_data = build_depth_track_dataframes(
        records,
        record_depth_tracks,
        base_config=depth_config,
        window_steps=[(10, 10), (10, 10)],
    )

    assert depth_data[0][0].id == "depth_B"
    assert depth_data[0][0].config.max_depth == pytest.approx(100)
    assert depth_data[1][1].config.max_depth == pytest.approx(100)
    assert depth_data[1][0].id == "depth_A"
    assert depth_data[1][0].config.max_depth == pytest.approx(5)


def test_linear_cli_forwards_table_inputs_to_api(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    input_table = tmp_path / "inputs.tsv"
    input_table.write_text(
        "input_id\tinput_type\tgbk\trecord_id\tlabel\n"
        "ref\tgbk\tref.gb\t#1\tReference\n"
        "sample\tgbk\tsample.gb\t#1\tSample\n",
        encoding="utf-8",
    )
    sample_depth = _write_depth_file(tmp_path / "sample_depth.tsv", "sample", 30)
    depth_table = tmp_path / "depth_table.tsv"
    depth_table.write_text(
        "record_id\ttrack_id\tfile\ttrack_label\ttrack_height\n"
        f"input:sample\tdepth_B\t{Path(sample_depth).name}\tSample B\t18\n",
        encoding="utf-8",
    )
    track_table = tmp_path / "track_table.tsv"
    track_table.write_text(
        "slot_id\trenderer\torder\tside\ttrack_id\theight\n"
        "features\tfeatures\t1\toverlay\t\t\n"
        "sample_depth\tdepth\t2\tbelow\tdepth_B\t\n",
        encoding="utf-8",
    )
    captured: dict[str, object] = {}

    def fake_load_gbk(paths, *_args, **_kwargs):
        stem = Path(paths[0]).stem
        return [_record(stem)]

    monkeypatch.setattr(linear_cli_module, "load_gbks", fake_load_gbk)
    monkeypatch.setattr(linear_cli_module, "read_color_table", lambda _path: None)
    monkeypatch.setattr(linear_cli_module, "load_default_colors", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(linear_cli_module, "read_feature_visibility_file", lambda _path: None)
    monkeypatch.setattr(linear_cli_module, "save_figure", lambda _canvas, _formats: None)

    def fake_assemble(*_args, **kwargs):
        captured.update(kwargs)
        return Drawing(filename=str(tmp_path / "dummy.svg"))

    monkeypatch.setattr(linear_cli_module, "assemble_linear_diagram_from_records", fake_assemble)

    linear_cli_module.linear_main(
        [
            "--input_table",
            str(input_table),
            "--depth_track_table",
            str(depth_table),
            "--track_table",
            str(track_table),
            "--format",
            "svg",
            "-o",
            str(tmp_path / "out"),
        ]
    )

    records = captured["records"]
    assert [record.annotations["gbdraw_input_id"] for record in records] == ["ref", "sample"]
    assert captured["depth_track_files"] is None
    record_depth_tracks = captured["record_depth_tracks"]
    assert [[spec.id for spec in row] for row in record_depth_tracks] == [[], ["depth_B"]]
    slots = captured["linear_track_slots"]
    assert slots[1].params["track_id"] == "depth_B"
    assert slots[1].params["track_index"] == 0


def test_circular_cli_depth_table_width_applies_to_default_slots(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    input_table = tmp_path / "inputs.tsv"
    input_table.write_text(
        "input_id\tinput_type\tgbk\trecord_id\tlabel\n"
        "ref\tgbk\tref.gb\t#1\tReference\n",
        encoding="utf-8",
    )
    depth_path = _write_depth_file(tmp_path / "depth.tsv", "ref", 20)
    depth_table = tmp_path / "depth_table.tsv"
    depth_table.write_text(
        "record_id\ttrack_id\tfile\ttrack_label\ttrack_width\n"
        f"input:ref\tdepth_A\t{Path(depth_path).name}\tSample A\t14\n",
        encoding="utf-8",
    )
    captured: dict[str, object] = {}

    monkeypatch.setattr(circular_cli_module, "load_gbks", lambda paths, *_args, **_kwargs: [_record(Path(paths[0]).stem)])
    monkeypatch.setattr(circular_cli_module, "read_color_table", lambda _path: None)
    monkeypatch.setattr(circular_cli_module, "read_feature_visibility_file", lambda _path: None)
    monkeypatch.setattr(circular_cli_module, "load_default_colors", lambda *_args, **_kwargs: None)
    monkeypatch.setattr(circular_cli_module, "save_figure", lambda _canvas, _formats: None)

    def fake_assemble(*_args, **kwargs):
        captured.update(kwargs)
        return Drawing(filename=str(tmp_path / "dummy.svg"))

    monkeypatch.setattr(circular_cli_module, "assemble_circular_diagram_from_record", fake_assemble)

    circular_cli_module.circular_main(
        [
            "--input_table",
            str(input_table),
            "--depth_track_table",
            str(depth_table),
            "--format",
            "svg",
            "-o",
            str(tmp_path / "out"),
        ]
    )

    assert captured["depth_track_files"] is None
    assert [[spec.id for spec in row] for row in captured["record_depth_tracks"]] == [["depth_A"]]
    slots = captured["circular_track_slots"]
    depth_slot = next(slot for slot in slots if slot.renderer == "depth")
    assert depth_slot.params["track_id"] == "depth_A"
    assert depth_slot.width is not None
    assert depth_slot.width.resolve(1) == pytest.approx(14)


def test_linear_cli_table_examples_run_end_to_end(project_root: Path, tmp_path: Path) -> None:
    examples_dir = project_root / "examples"
    output_prefix = tmp_path / "linear_cli_table"
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "gbdraw.cli",
            "linear",
            "--input_table",
            str(examples_dir / "cli_table_inputs.tsv"),
            "--depth_track_table",
            str(examples_dir / "cli_table_depth_tracks.tsv"),
            "--track_table",
            str(examples_dir / "cli_table_track_slots.tsv"),
            "-o",
            str(output_prefix),
            "-f",
            "svg",
        ],
        cwd=str(project_root),
        capture_output=True,
        text=True,
        timeout=120,
    )

    output_svg = output_prefix.with_suffix(".svg")
    assert result.returncode == 0, result.stdout + result.stderr
    assert output_svg.exists()
    assert "<svg" in output_svg.read_text(encoding="utf-8")


def test_circular_cli_table_examples_run_end_to_end(project_root: Path, tmp_path: Path) -> None:
    examples_dir = project_root / "examples"
    output_prefix = tmp_path / "circular_cli_table"
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "gbdraw.cli",
            "circular",
            "--input_table",
            str(examples_dir / "cli_table_circular_inputs.tsv"),
            "--depth_track_table",
            str(examples_dir / "cli_table_circular_depth_tracks.tsv"),
            "--track_table",
            str(examples_dir / "cli_table_circular_track_slots.tsv"),
            "--suppress_gc",
            "--suppress_skew",
            "-o",
            str(output_prefix),
            "-f",
            "svg",
        ],
        cwd=str(project_root),
        capture_output=True,
        text=True,
        timeout=120,
    )

    output_svg = output_prefix.with_suffix(".svg")
    assert result.returncode == 0, result.stdout + result.stderr
    assert output_svg.exists()
    assert "<svg" in output_svg.read_text(encoding="utf-8")
