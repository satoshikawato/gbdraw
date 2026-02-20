from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

import pytest

from gbdraw.layout.linear import calculate_feature_position_factors_linear
from gbdraw.linear import _parse_linear_track_layout


INPUT_GBK = Path(__file__).parent / "test_inputs" / "MjeNMV.gb"


def _run_linear(tmp_path: Path, extra_args: list[str]) -> tuple[int, str, str, Path]:
    output_name = "linear_track_layout_test"
    cmd = [
        sys.executable,
        "-m",
        "gbdraw.cli",
        "linear",
        "--gbk",
        str(INPUT_GBK),
        "-o",
        output_name,
        "-f",
        "svg",
        "--legend",
        "none",
    ] + list(extra_args)
    result = subprocess.run(cmd, cwd=str(tmp_path), capture_output=True, text=True, timeout=240)
    return result.returncode, result.stdout, result.stderr, (tmp_path / f"{output_name}.svg")


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

