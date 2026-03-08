from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

INPUT_MG1655 = Path(__file__).parent / "test_inputs" / "MG1655.gbk"
SVG_NS = {"svg": "http://www.w3.org/2000/svg"}


def _extract_group_translate_y(root: ET.Element, group_id: str) -> float:
    group = root.find(f".//svg:g[@id='{group_id}']", SVG_NS)
    assert group is not None
    transform = group.attrib.get("transform", "")
    match = re.search(r"translate\(\s*[-0-9.]+(?:px)?[\s,]+([-0-9.]+)", transform)
    assert match is not None
    return float(match.group(1))


def _extract_group_font_sizes(root: ET.Element, group_id: str) -> list[float]:
    group = root.find(f".//svg:g[@id='{group_id}']", SVG_NS)
    assert group is not None
    sizes: list[float] = []
    for text in group.findall(".//svg:text", SVG_NS):
        raw_size = text.attrib.get("font-size")
        if raw_size is None:
            continue
        try:
            sizes.append(float(raw_size))
        except ValueError:
            continue
    return sizes


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
    assert "1-200" in svg_content
    assert "200 bp" not in svg_content
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
    assert "200-1" in svg_content
    assert "200 bp" not in svg_content
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
def test_linear_reverse_complement_without_region_keeps_length_label(
    temp_output_dir: Path,
    gbdraw_runner,
) -> None:
    gbk_path = temp_output_dir / "multi_records_rc_length.gbk"
    lengths = _write_multi_record_gbk(gbk_path)

    returncode, output, svg_path = gbdraw_runner.run_linear(
        [gbk_path],
        "linear_reverse_complement_keeps_length",
        temp_output_dir,
        extra_args=["--record_id", "RecA", "--reverse_complement", "true", "--legend", "none"],
    )

    assert returncode == 0, f"gbdraw failed: {output}"
    svg_content = svg_path.read_text()
    assert f"{lengths['RecA']:,} bp" in svg_content
    assert "1234-1" not in svg_content


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
    assert "101 bp" not in svg_content
    assert "300 bp" not in svg_content
    assert "200 bp" in svg_content


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
            "50",
            "--legend",
            "none",
        ],
    )

    assert returncode == 0, f"gbdraw failed: {output}"
    svg_content = svg_path.read_text(encoding="utf-8")
    assert 'id="length_bar"' not in svg_content
    assert "300 bp" not in svg_content
    assert "101 bp" not in svg_content
    assert "250 bp" in svg_content
    assert "200 bp" in svg_content
    assert "150 bp" in svg_content
    assert svg_content.index("250 bp") < svg_content.index("200 bp") < svg_content.index("150 bp")


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


@pytest.mark.linear
def test_linear_plot_title_drawn_once_and_coexists_with_record_labels(
    temp_output_dir: Path,
    gbdraw_runner,
) -> None:
    gbk_path = temp_output_dir / "multi_records_plot_title.gbk"
    _write_multi_record_gbk(gbk_path)

    returncode, output, svg_path = gbdraw_runner.run_linear(
        [gbk_path],
        "linear_plot_title_with_labels",
        temp_output_dir,
        extra_args=[
            "--plot_title",
            "Global Plot Title",
            "--record_label",
            "Record A Label",
            "--record_label",
            "Record B Label",
            "--legend",
            "none",
        ],
    )

    assert returncode == 0, f"gbdraw failed: {output}"
    svg_content = svg_path.read_text(encoding="utf-8")
    root = ET.fromstring(svg_content)
    plot_title_groups = root.findall(".//svg:g[@id='plot_title']", SVG_NS)
    assert len(plot_title_groups) == 1
    assert "Global Plot Title" in "".join(plot_title_groups[0].itertext())
    assert "Record A Label" in svg_content
    assert "Record B Label" in svg_content
    assert root.find(".//svg:g[@id='RecA']", SVG_NS) is not None
    assert root.find(".//svg:g[@id='RecB']", SVG_NS) is not None


@pytest.mark.linear
def test_linear_plot_title_position_changes_vertical_location(
    temp_output_dir: Path,
    gbdraw_runner,
) -> None:
    gbk_path = temp_output_dir / "multi_records_plot_title_position.gbk"
    _write_multi_record_gbk(gbk_path)

    y_positions: dict[str, float] = {}
    for position in ("top", "center", "bottom"):
        returncode, output, svg_path = gbdraw_runner.run_linear(
            [gbk_path],
            f"linear_plot_title_{position}",
            temp_output_dir,
            extra_args=[
                "--plot_title",
                "Position Test",
                "--plot_title_position",
                position,
                "--legend",
                "none",
            ],
        )
        assert returncode == 0, f"gbdraw failed: {output}"
        root = ET.fromstring(svg_path.read_text(encoding="utf-8"))
        y_positions[position] = _extract_group_translate_y(root, "plot_title")

    assert y_positions["top"] < y_positions["center"] < y_positions["bottom"]


@pytest.mark.linear
def test_linear_plot_title_font_size_default_and_override(
    temp_output_dir: Path,
    gbdraw_runner,
) -> None:
    gbk_path = temp_output_dir / "multi_records_plot_title_size.gbk"
    _write_multi_record_gbk(gbk_path)

    default_returncode, default_output, default_svg = gbdraw_runner.run_linear(
        [gbk_path],
        "linear_plot_title_default_size",
        temp_output_dir,
        extra_args=["--plot_title", "Default Size", "--legend", "none"],
    )
    assert default_returncode == 0, f"gbdraw failed: {default_output}"
    default_root = ET.fromstring(default_svg.read_text(encoding="utf-8"))
    default_sizes = _extract_group_font_sizes(default_root, "plot_title")
    assert default_sizes == [32.0]

    override_returncode, override_output, override_svg = gbdraw_runner.run_linear(
        [gbk_path],
        "linear_plot_title_override_size",
        temp_output_dir,
        extra_args=[
            "--plot_title",
            "Override Size",
            "--plot_title_font_size",
            "40",
            "--legend",
            "none",
        ],
    )
    assert override_returncode == 0, f"gbdraw failed: {override_output}"
    override_root = ET.fromstring(override_svg.read_text(encoding="utf-8"))
    override_sizes = _extract_group_font_sizes(override_root, "plot_title")
    assert override_sizes == [40.0]
