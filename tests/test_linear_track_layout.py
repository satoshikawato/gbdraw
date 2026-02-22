from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path

import pytest

from gbdraw.layout.linear import calculate_feature_position_factors_linear
from gbdraw.linear import _parse_linear_track_axis_gap, _parse_linear_track_layout


INPUT_GBK = Path(__file__).parent / "test_inputs" / "MjeNMV.gb"
INPUT_MG1655 = Path(__file__).parent / "test_inputs" / "MG1655.gbk"


def _run_linear_with_gbks(
    tmp_path: Path,
    gbk_paths: list[Path],
    extra_args: list[str],
    *,
    output_name: str = "linear_track_layout_test",
) -> tuple[int, str, str, Path]:
    cmd = [
        sys.executable,
        "-m",
        "gbdraw.cli",
        "linear",
        "--gbk",
        *(str(path) for path in gbk_paths),
        "-o",
        output_name,
        "-f",
        "svg",
        "--legend",
        "none",
    ] + list(extra_args)
    result = subprocess.run(cmd, cwd=str(tmp_path), capture_output=True, text=True, timeout=240)
    return result.returncode, result.stdout, result.stderr, (tmp_path / f"{output_name}.svg")


def _run_linear(tmp_path: Path, extra_args: list[str]) -> tuple[int, str, str, Path]:
    return _run_linear_with_gbks(tmp_path, [INPUT_GBK], extra_args)


@pytest.mark.linear
def test_parse_linear_track_layout_aliases() -> None:
    assert _parse_linear_track_layout("above") == "above"
    assert _parse_linear_track_layout("middle") == "middle"
    assert _parse_linear_track_layout("below") == "below"
    assert _parse_linear_track_layout("spreadout") == "above"
    assert _parse_linear_track_layout("tuckin") == "below"
    with pytest.raises(argparse.ArgumentTypeError):
        _parse_linear_track_layout("invalid")


@pytest.mark.linear
def test_parse_linear_track_axis_gap() -> None:
    assert _parse_linear_track_axis_gap("auto") is None
    assert _parse_linear_track_axis_gap("0") == pytest.approx(0.0)
    assert _parse_linear_track_axis_gap("12.5") == pytest.approx(12.5)
    with pytest.raises(argparse.ArgumentTypeError):
        _parse_linear_track_axis_gap("-1")
    with pytest.raises(argparse.ArgumentTypeError):
        _parse_linear_track_axis_gap("invalid")


@pytest.mark.linear
def test_linear_layout_factors_non_stranded_above_below() -> None:
    above_track0 = calculate_feature_position_factors_linear("positive", 0, False, track_layout="above")
    above_track1 = calculate_feature_position_factors_linear("positive", 1, False, track_layout="above")
    above_track2 = calculate_feature_position_factors_linear("positive", 2, False, track_layout="above")
    assert all(v < 0 for v in above_track0)
    assert above_track2[1] < above_track1[1] < above_track0[1]

    below_track0 = calculate_feature_position_factors_linear("positive", 0, False, track_layout="below")
    below_track1 = calculate_feature_position_factors_linear("positive", 1, False, track_layout="below")
    below_track2 = calculate_feature_position_factors_linear("positive", 2, False, track_layout="below")
    assert all(v > 0 for v in below_track0)
    assert below_track2[1] > below_track1[1] > below_track0[1]


@pytest.mark.linear
def test_linear_middle_layout_separate_strands_is_axis_symmetric() -> None:
    pos0 = calculate_feature_position_factors_linear("positive", 0, True, track_layout="middle")
    neg1 = calculate_feature_position_factors_linear("negative", -1, True, track_layout="middle")
    assert abs(pos0[1]) == pytest.approx(abs(neg1[1]))
    assert abs(pos0[2]) == pytest.approx(abs(neg1[0]))

    pos1 = calculate_feature_position_factors_linear("positive", 1, True, track_layout="middle")
    neg2 = calculate_feature_position_factors_linear("negative", -2, True, track_layout="middle")
    assert abs(pos1[1]) == pytest.approx(abs(neg2[1]))
    assert abs(pos1[2]) == pytest.approx(abs(neg2[0]))


@pytest.mark.linear
def test_linear_layout_axis_gap_factor_increases_distance_from_axis() -> None:
    auto_above = calculate_feature_position_factors_linear("positive", 0, False, track_layout="above")
    explicit_above = calculate_feature_position_factors_linear(
        "positive",
        0,
        False,
        track_layout="above",
        axis_gap_factor=0.8,
    )
    assert explicit_above[2] < auto_above[2]

    auto_below = calculate_feature_position_factors_linear("positive", 0, False, track_layout="below")
    explicit_below = calculate_feature_position_factors_linear(
        "positive",
        0,
        False,
        track_layout="below",
        axis_gap_factor=0.8,
    )
    assert explicit_below[0] > auto_below[0]


@pytest.mark.linear
def test_linear_layout_factors_separate_strands_above_below_keep_order() -> None:
    above_pos_0 = calculate_feature_position_factors_linear("positive", 0, True, track_layout="above")
    above_pos_2 = calculate_feature_position_factors_linear("positive", 2, True, track_layout="above")
    above_neg_1 = calculate_feature_position_factors_linear("negative", -1, True, track_layout="above")
    above_neg_3 = calculate_feature_position_factors_linear("negative", -3, True, track_layout="above")
    assert all(v < 0 for v in above_pos_0)
    assert all(v < 0 for v in above_neg_1)
    assert above_pos_0[1] < above_neg_1[1]
    assert above_pos_2[1] < above_pos_0[1]
    assert above_neg_3[1] < above_neg_1[1]

    below_pos_0 = calculate_feature_position_factors_linear("positive", 0, True, track_layout="below")
    below_pos_2 = calculate_feature_position_factors_linear("positive", 2, True, track_layout="below")
    below_neg_1 = calculate_feature_position_factors_linear("negative", -1, True, track_layout="below")
    below_neg_3 = calculate_feature_position_factors_linear("negative", -3, True, track_layout="below")
    assert all(v > 0 for v in below_pos_0)
    assert all(v > 0 for v in below_neg_1)
    assert below_pos_0[1] < below_neg_1[1]
    assert below_pos_2[1] > below_pos_0[1]
    assert below_neg_3[1] > below_neg_1[1]


@pytest.mark.linear
@pytest.mark.parametrize("layout", ["above", "below", "spreadout", "tuckin"])
def test_linear_track_layout_cli_generates_svg(tmp_path: Path, layout: str) -> None:
    returncode, stdout, stderr, output_svg = _run_linear(
        tmp_path,
        ["--track_layout", layout],
    )
    assert returncode == 0, f"stdout={stdout}\nstderr={stderr}"
    assert output_svg.exists()


@pytest.mark.linear
@pytest.mark.parametrize("layout", ["above", "below"])
def test_linear_track_layout_with_separate_and_resolve_overlaps(tmp_path: Path, layout: str) -> None:
    returncode, stdout, stderr, output_svg = _run_linear(
        tmp_path,
        ["--track_layout", layout, "--separate_strands", "--resolve_overlaps"],
    )
    assert returncode == 0, f"stdout={stdout}\nstderr={stderr}"
    assert output_svg.exists()


@pytest.mark.linear
def test_linear_track_axis_gap_cli_generates_svg(tmp_path: Path) -> None:
    returncode, stdout, stderr, output_svg = _run_linear(
        tmp_path,
        ["--track_layout", "below", "--track_axis_gap", "12"],
    )
    assert returncode == 0, f"stdout={stdout}\nstderr={stderr}"
    assert output_svg.exists()


@pytest.mark.linear
@pytest.mark.parametrize("layout", ["above", "below"])
def test_linear_ruler_on_axis_hides_bottom_length_bar(tmp_path: Path, layout: str) -> None:
    returncode, stdout, stderr, output_svg = _run_linear(
        tmp_path,
        ["--track_layout", layout, "--scale_style", "ruler", "--ruler_on_axis", "--scale_interval", "50000"],
    )
    assert returncode == 0, f"stdout={stdout}\nstderr={stderr}"
    svg_content = output_svg.read_text(encoding="utf-8")
    assert 'id="length_bar"' not in svg_content
    assert "kbp" in svg_content


@pytest.mark.linear
@pytest.mark.parametrize(
    ("extra_args", "expect_warning"),
    [
        (["--track_layout", "middle", "--scale_style", "ruler", "--ruler_on_axis"], True),
        (["--track_layout", "above", "--scale_style", "bar", "--ruler_on_axis"], True),
        (["--track_layout", "above", "--scale_style", "ruler"], False),
    ],
)
def test_linear_ruler_on_axis_invalid_conditions_keep_legacy_length_bar(
    tmp_path: Path,
    extra_args: list[str],
    expect_warning: bool,
) -> None:
    returncode, stdout, stderr, output_svg = _run_linear(tmp_path, extra_args)
    assert returncode == 0, f"stdout={stdout}\nstderr={stderr}"
    svg_content = output_svg.read_text(encoding="utf-8")
    assert 'id="length_bar"' in svg_content
    merged_output = f"{stdout}\n{stderr}"
    warning_text = "--ruler_on_axis is ignored unless --scale_style ruler and --track_layout above|below are set."
    if expect_warning:
        assert warning_text in merged_output
    else:
        assert warning_text not in merged_output


@pytest.mark.linear
def test_linear_ruler_on_axis_uses_shared_auto_interval_across_records(tmp_path: Path) -> None:
    returncode, stdout, stderr, output_svg = _run_linear_with_gbks(
        tmp_path,
        [INPUT_MG1655, INPUT_MG1655, INPUT_MG1655],
        [
            "--track_layout",
            "above",
            "--scale_style",
            "ruler",
            "--ruler_on_axis",
            "--region",
            "#1:1028779-1035047",
            "--region",
            "#2:1028779-1035047",
            "--region",
            "#3:1028779-1029277",
        ],
        output_name="linear_track_layout_shared_interval",
    )
    assert returncode == 0, f"stdout={stdout}\nstderr={stderr}"
    svg_content = output_svg.read_text(encoding="utf-8")
    group_matches = re.findall(r'<g id="NC_000913\.3"[^>]*>(.*?)</g>', svg_content)
    axis_groups = [g for g in group_matches if re.search(r'<line[^>]*y2="10\.0"', g) is not None]
    assert len(axis_groups) == 3
    tick_counts = [len(re.findall(r'<line[^>]*y2="10\.0"', g)) for g in axis_groups]
    assert tick_counts[2] <= 1
    assert tick_counts[0] > tick_counts[2]


@pytest.mark.linear
def test_linear_ruler_on_axis_defaults_tick_and_label_colors_to_axis(tmp_path: Path) -> None:
    returncode, stdout, stderr, output_svg = _run_linear(
        tmp_path,
        ["--track_layout", "above", "--scale_style", "ruler", "--ruler_on_axis", "--scale_interval", "50000"],
    )
    assert returncode == 0, f"stdout={stdout}\nstderr={stderr}"
    svg_content = output_svg.read_text(encoding="utf-8")
    axis_match = re.search(
        r'<line fill="none" stroke="([^"]+)" stroke-width="([^"]+)"[^>]*x1="0"[^>]*x2="2000\.0"[^>]*y1="0"[^>]*y2="0"',
        svg_content,
    )
    tick_match = re.search(
        r'<line[^>]*stroke="([^"]+)"[^>]*stroke-width="([^"]+)"[^>]*y2="10\.0"',
        svg_content,
    )
    label_match = re.search(r'<text[^>]*fill="([^"]+)"[^>]*>[^<]*(?:bp|kbp|Mbp)</text>', svg_content)
    assert axis_match is not None
    assert tick_match is not None
    assert label_match is not None
    assert tick_match.group(1) == axis_match.group(1)
    assert float(tick_match.group(2)) == pytest.approx(float(axis_match.group(2)))
    assert label_match.group(1) == axis_match.group(1)


@pytest.mark.linear
def test_linear_ruler_label_options_apply_to_axis_ruler(tmp_path: Path) -> None:
    returncode, stdout, stderr, output_svg = _run_linear(
        tmp_path,
        [
            "--track_layout",
            "above",
            "--scale_style",
            "ruler",
            "--ruler_on_axis",
            "--scale_interval",
            "50000",
            "--ruler_label_font_size",
            "22",
            "--ruler_label_color",
            "tomato",
        ],
    )
    assert returncode == 0, f"stdout={stdout}\nstderr={stderr}"
    svg_content = output_svg.read_text(encoding="utf-8")
    ruler_text_tags = re.findall(r'<text[^>]*>[^<]*(?:bp|kbp|Mbp)</text>', svg_content)
    assert any('fill="tomato"' in tag and 'font-size="22.0"' in tag for tag in ruler_text_tags)


@pytest.mark.linear
def test_linear_ruler_label_options_apply_to_bottom_ruler(tmp_path: Path) -> None:
    returncode, stdout, stderr, output_svg = _run_linear(
        tmp_path,
        [
            "--scale_style",
            "ruler",
            "--scale_interval",
            "50000",
            "--ruler_label_font_size",
            "21",
            "--ruler_label_color",
            "teal",
        ],
    )
    assert returncode == 0, f"stdout={stdout}\nstderr={stderr}"
    svg_content = output_svg.read_text(encoding="utf-8")
    length_bar_match = re.search(r'<g id="length_bar"[^>]*>(.*?)</g>', svg_content)
    assert length_bar_match is not None
    ruler_text_tags = re.findall(r'<text[^>]*>[^<]*(?:bp|kbp|Mbp)</text>', length_bar_match.group(1))
    assert any('fill="teal"' in tag and 'font-size="21.0"' in tag for tag in ruler_text_tags)
