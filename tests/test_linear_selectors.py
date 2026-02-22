from __future__ import annotations

import re
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

INPUT_MG1655 = Path(__file__).parent / "test_inputs" / "MG1655.gbk"


def _write_multi_record_gbk(path: Path) -> dict[str, int]:
    rec_a_len = 1234
    rec_b_len = 567

    rec_a = SeqRecord(
        Seq("A" * rec_a_len),
        id="RecA",
        name="RecA",
        description="Record A",
    )
    rec_a.annotations["molecule_type"] = "DNA"
    rec_a.features = [SeqFeature(FeatureLocation(0, 120), type="CDS")]

    rec_b = SeqRecord(
        Seq("T" * rec_b_len),
        id="RecB",
        name="RecB",
        description="Record B",
    )
    rec_b.annotations["molecule_type"] = "DNA"
    rec_b.features = [SeqFeature(FeatureLocation(0, 80), type="CDS")]

    SeqIO.write([rec_a, rec_b], path, "genbank")
    return {"RecA": rec_a_len, "RecB": rec_b_len}


@pytest.mark.linear
def test_linear_record_id_and_label(temp_output_dir: Path, gbdraw_runner) -> None:
    gbk_path = temp_output_dir / "multi_records.gbk"
    lengths = _write_multi_record_gbk(gbk_path)

    returncode, output, svg_path = gbdraw_runner.run_linear(
        [gbk_path],
        "linear_record_id_label",
        temp_output_dir,
        extra_args=["--record_id", "RecB", "--record_label", "CustomLabel", "--legend", "none"],
    )

    assert returncode == 0, f"gbdraw failed: {output}"
    svg_content = svg_path.read_text()
    assert "CustomLabel" in svg_content
    assert f"{lengths['RecB']:,} bp" in svg_content
    assert f"{lengths['RecA']:,} bp" not in svg_content


@pytest.mark.linear
def test_linear_region_crop_with_record_selector(temp_output_dir: Path, gbdraw_runner) -> None:
    gbk_path = temp_output_dir / "multi_records_region.gbk"
    lengths = _write_multi_record_gbk(gbk_path)

    returncode, output, svg_path = gbdraw_runner.run_linear(
        [gbk_path],
        "linear_region_crop",
        temp_output_dir,
        extra_args=["--region", "RecA:1-200", "--legend", "none"],
    )

    assert returncode == 0, f"gbdraw failed: {output}"
    svg_content = svg_path.read_text()
    assert "200 bp" in svg_content
    assert f"{lengths['RecB']:,} bp" in svg_content
    assert f"{lengths['RecA']:,} bp" not in svg_content


@pytest.mark.linear
def test_linear_region_crop_with_file_prefix(temp_output_dir: Path, gbdraw_runner) -> None:
    gbk_path = temp_output_dir / "multi_records_fileprefix.gbk"
    lengths = _write_multi_record_gbk(gbk_path)
    file_selector = gbk_path.name

    returncode, output, svg_path = gbdraw_runner.run_linear(
        [gbk_path],
        "linear_region_fileprefix",
        temp_output_dir,
        extra_args=["--region", f"{file_selector}:RecA:1-200:rc", "--legend", "none"],
    )

    assert returncode == 0, f"gbdraw failed: {output}"
    assert "reverse complement" in output
    svg_content = svg_path.read_text()
    assert "200 bp" in svg_content
    assert f"{lengths['RecB']:,} bp" in svg_content
    assert f"{lengths['RecA']:,} bp" not in svg_content


@pytest.mark.linear
def test_linear_reverse_complement_flag(temp_output_dir: Path, gbdraw_runner) -> None:
    gbk_path = temp_output_dir / "multi_records_rc.gbk"
    _write_multi_record_gbk(gbk_path)

    returncode, output, svg_path = gbdraw_runner.run_linear(
        [gbk_path],
        "linear_reverse_complement",
        temp_output_dir,
        extra_args=["--record_id", "RecA", "--reverse_complement", "true", "--legend", "none"],
    )

    assert returncode == 0, f"gbdraw failed: {output}"
    assert "Reverse complemented 1 record(s)." in output
    assert svg_path.exists()


@pytest.mark.linear
def test_linear_region_ruler_on_axis_uses_absolute_coordinates(
    temp_output_dir: Path,
    gbdraw_runner,
) -> None:
    gbk_path = temp_output_dir / "multi_records_region_axis.gbk"
    _write_multi_record_gbk(gbk_path)

    returncode, output, svg_path = gbdraw_runner.run_linear(
        [gbk_path],
        "linear_region_axis_coords",
        temp_output_dir,
        extra_args=[
            "--record_id",
            "RecA",
            "--region",
            "RecA:101-300",
            "--track_layout",
            "above",
            "--scale_style",
            "ruler",
            "--ruler_on_axis",
            "--scale_interval",
            "100",
            "--legend",
            "none",
        ],
    )

    assert returncode == 0, f"gbdraw failed: {output}"
    svg_content = svg_path.read_text(encoding="utf-8")
    assert 'id="length_bar"' not in svg_content
    assert "101 bp" in svg_content
    assert "300 bp" in svg_content


@pytest.mark.linear
def test_linear_region_ruler_on_axis_rc_labels_descend_left_to_right(
    temp_output_dir: Path,
    gbdraw_runner,
) -> None:
    gbk_path = temp_output_dir / "multi_records_region_axis_rc.gbk"
    _write_multi_record_gbk(gbk_path)

    returncode, output, svg_path = gbdraw_runner.run_linear(
        [gbk_path],
        "linear_region_axis_coords_rc",
        temp_output_dir,
        extra_args=[
            "--record_id",
            "RecA",
            "--region",
            "RecA:101-300:rc",
            "--track_layout",
            "above",
            "--scale_style",
            "ruler",
            "--ruler_on_axis",
            "--scale_interval",
            "100",
            "--legend",
            "none",
        ],
    )

    assert returncode == 0, f"gbdraw failed: {output}"
    svg_content = svg_path.read_text(encoding="utf-8")
    assert 'id="length_bar"' not in svg_content
    assert "300 bp" in svg_content
    assert "101 bp" in svg_content
    assert svg_content.index("300 bp") < svg_content.index("101 bp")


@pytest.mark.linear
def test_linear_region_ruler_on_axis_uses_span_based_units_for_high_coordinates(
    temp_output_dir: Path,
    gbdraw_runner,
) -> None:
    returncode, output, svg_path = gbdraw_runner.run_linear(
        [INPUT_MG1655],
        "linear_region_axis_high_coords",
        temp_output_dir,
        extra_args=[
            "--region",
            "NC_000913.3:1028779-1035047",
            "--track_layout",
            "above",
            "--scale_style",
            "ruler",
            "--ruler_on_axis",
            "--legend",
            "none",
        ],
    )

    assert returncode == 0, f"gbdraw failed: {output}"
    svg_content = svg_path.read_text(encoding="utf-8")
    assert 'id="length_bar"' not in svg_content
    assert re.search(r">\d+(?:\.\d+)? kbp<", svg_content) is not None
    assert "Mbp" not in svg_content
