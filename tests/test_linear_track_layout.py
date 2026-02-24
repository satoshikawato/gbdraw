from __future__ import annotations

import argparse
import re
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

from gbdraw.layout.linear import calculate_feature_position_factors_linear
from gbdraw.linear import _parse_linear_track_axis_gap, _parse_linear_track_layout
from gbdraw.render.groups.linear.length_bar import RULER_TICK_LENGTH


INPUT_GBK = Path(__file__).parent / "test_inputs" / "MjeNMV.gb"
INPUT_MELA_GBK = Path(__file__).parent / "test_inputs" / "MelaMJNV.gb"
INPUT_MG1655 = Path(__file__).parent / "test_inputs" / "MG1655.gbk"
INPUT_HMMTDNA = Path(__file__).parent / "test_inputs" / "HmmtDNA.gbk"
INPUT_MJE_MELA_BLAST = Path(__file__).parent / "test_inputs" / "MjeNMV.MelaMJNV.tblastx.out"


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


def _extract_axis_group_y(svg_content: str, record_id: str) -> float:
    match = re.search(
        rf'<g id="{re.escape(record_id)}" transform="translate\([^,]+,([0-9.]+)\)"><line[^>]*y1="0" y2="0"',
        svg_content,
    )
    assert match is not None
    return float(match.group(1))


def _extract_comparison_group_y(svg_content: str) -> float:
    match = re.search(r'<g id="comparison1" transform="translate\([^,]+,([0-9.]+)\)">', svg_content)
    assert match is not None
    return float(match.group(1))


def _extract_first_comparison_path_ys(svg_content: str) -> tuple[float, float, float, float]:
    match = re.search(r'<g id="comparison1"[^>]*><path d="([^"]+)"', svg_content)
    assert match is not None
    coords = [float(value) for value in re.findall(r"-?\d+(?:\.\d+)?", match.group(1))]
    assert len(coords) >= 8
    return coords[1], coords[3], coords[5], coords[7]


def _is_ruler_tick_y2(y2: float) -> bool:
    return abs(abs(y2) - float(RULER_TICK_LENGTH)) < 1e-6


def _count_ruler_ticks(svg_fragment: str) -> int:
    count = 0
    for match in re.finditer(r'<line[^>]*y2="(-?\d+(?:\.\d+)?)"', svg_fragment):
        if _is_ruler_tick_y2(float(match.group(1))):
            count += 1
    return count


def _extract_ruler_text_tags(svg_fragment: str) -> list[str]:
    return re.findall(r'<text[^>]*>[^<]*(?:bp|kbp|Mbp)</text>', svg_fragment)


def _extract_min_absolute_feature_y(svg_content: str) -> float:
    root = ET.fromstring(svg_content)
    namespace = {"svg": "http://www.w3.org/2000/svg"}
    min_absolute_y: float | None = None

    for group in root.findall("svg:g", namespace):
        group_id = group.attrib.get("id", "")
        if group_id == "length_bar" or group_id == "legend" or group_id.startswith("comparison"):
            continue
        axis_line = group.find('svg:line[@y1="0"][@y2="0"]', namespace)
        if axis_line is None:
            continue
        transform = group.attrib.get("transform", "")
        transform_match = re.search(r'translate\([^,]+,([\-0-9.]+)\)', transform)
        if transform_match is None:
            continue
        axis_y = float(transform_match.group(1))
        for path in group.findall("svg:path", namespace):
            coords = [float(value) for value in re.findall(r"-?\d+(?:\.\d+)?", path.attrib.get("d", ""))]
            y_values = coords[1::2]
            if not y_values:
                continue
            absolute_top_y = axis_y + min(y_values)
            if min_absolute_y is None or absolute_top_y < min_absolute_y:
                min_absolute_y = absolute_top_y

    assert min_absolute_y is not None
    return min_absolute_y


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
@pytest.mark.parametrize("i", [1, 2, 3, 4])
def test_linear_below_displaced_positive_tracks_align_with_matching_negative_tracks(i: int) -> None:
    pos_i = calculate_feature_position_factors_linear("positive", i, True, track_layout="below")
    neg_i = calculate_feature_position_factors_linear("negative", -i, True, track_layout="below")
    assert pos_i[1] == pytest.approx(neg_i[1])


@pytest.mark.linear
@pytest.mark.parametrize("i", [0, 1, 2, 3])
def test_linear_above_positive_tracks_align_with_negative_tracks_offset_by_two(i: int) -> None:
    pos_i = calculate_feature_position_factors_linear("positive", i, True, track_layout="above")
    neg_offset_i = calculate_feature_position_factors_linear("negative", -(i + 2), True, track_layout="above")
    assert pos_i[1] == pytest.approx(neg_offset_i[1])


@pytest.mark.linear
def test_linear_above_displaced_positive_tracks_move_outward_from_axis() -> None:
    pos_0 = calculate_feature_position_factors_linear("positive", 0, True, track_layout="above")
    pos_1 = calculate_feature_position_factors_linear("positive", 1, True, track_layout="above")
    pos_2 = calculate_feature_position_factors_linear("positive", 2, True, track_layout="above")
    assert pos_1[1] < pos_0[1]
    assert pos_2[1] < pos_1[1]


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
def test_linear_pairwise_hits_span_between_axes_for_non_middle_layouts(tmp_path: Path, layout: str) -> None:
    returncode, stdout, stderr, output_svg = _run_linear_with_gbks(
        tmp_path,
        [INPUT_GBK, INPUT_MELA_GBK],
        ["--track_layout", layout, "-b", str(INPUT_MJE_MELA_BLAST)],
        output_name=f"linear_pairwise_axis_{layout}",
    )
    assert returncode == 0, f"stdout={stdout}\nstderr={stderr}"
    svg_content = output_svg.read_text(encoding="utf-8")

    first_axis_y = _extract_axis_group_y(svg_content, "LC738868.1")
    second_axis_y = _extract_axis_group_y(svg_content, "LC738874.1")
    comparison_group_y = _extract_comparison_group_y(svg_content)
    path_y1, path_y2, path_y3, path_y4 = _extract_first_comparison_path_ys(svg_content)
    axis_delta = second_axis_y - first_axis_y

    assert comparison_group_y == pytest.approx(first_axis_y)
    assert path_y1 == pytest.approx(0.0)
    assert path_y2 == pytest.approx(0.0)
    assert path_y3 == pytest.approx(axis_delta)
    assert path_y4 == pytest.approx(axis_delta)


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
def test_linear_middle_layout_keeps_feature_top_inside_canvas_with_resolve_overlaps(tmp_path: Path) -> None:
    returncode, stdout, stderr, output_svg = _run_linear_with_gbks(
        tmp_path,
        [INPUT_HMMTDNA],
        ["--track_layout", "middle", "--separate_strands", "--resolve_overlaps", "--show_labels", "none"],
        output_name="linear_middle_resolve_overlaps_canvas_top",
    )
    assert returncode == 0, f"stdout={stdout}\nstderr={stderr}"
    svg_content = output_svg.read_text(encoding="utf-8")
    assert _extract_min_absolute_feature_y(svg_content) >= -1e-6


@pytest.mark.linear
def test_linear_above_layout_keeps_feature_top_inside_canvas_with_resolve_overlaps(tmp_path: Path) -> None:
    returncode, stdout, stderr, output_svg = _run_linear_with_gbks(
        tmp_path,
        [INPUT_HMMTDNA],
        ["--track_layout", "above", "--separate_strands", "--resolve_overlaps", "--show_labels", "none"],
        output_name="linear_above_resolve_overlaps_canvas_top",
    )
    assert returncode == 0, f"stdout={stdout}\nstderr={stderr}"
    svg_content = output_svg.read_text(encoding="utf-8")
    assert _extract_min_absolute_feature_y(svg_content) >= -1e-6


@pytest.mark.linear
@pytest.mark.parametrize("layout", ["above", "middle"])
def test_linear_layout_preserves_top_margin_when_resolve_overlaps_enabled(
    tmp_path: Path,
    layout: str,
) -> None:
    base_args = [
        "--track_layout",
        layout,
        "--show_labels",
        "none",
        "--scale_style",
        "ruler",
        "--ruler_on_axis",
    ]
    base_returncode, base_stdout, base_stderr, base_output_svg = _run_linear_with_gbks(
        tmp_path,
        [INPUT_HMMTDNA],
        base_args,
        output_name=f"linear_{layout}_top_margin_resolve_baseline",
    )
    resolved_returncode, resolved_stdout, resolved_stderr, resolved_output_svg = _run_linear_with_gbks(
        tmp_path,
        [INPUT_HMMTDNA],
        [*base_args, "--resolve_overlaps"],
        output_name=f"linear_{layout}_top_margin_resolve_enabled",
    )
    assert base_returncode == 0, f"stdout={base_stdout}\nstderr={base_stderr}"
    assert resolved_returncode == 0, f"stdout={resolved_stdout}\nstderr={resolved_stderr}"
    base_svg_content = base_output_svg.read_text(encoding="utf-8")
    resolved_svg_content = resolved_output_svg.read_text(encoding="utf-8")
    base_top_margin = _extract_min_absolute_feature_y(base_svg_content)
    resolved_top_margin = _extract_min_absolute_feature_y(resolved_svg_content)
    assert resolved_top_margin == pytest.approx(base_top_margin, abs=1e-6)


@pytest.mark.linear
@pytest.mark.parametrize("resolve_overlaps", [False, True])
def test_linear_above_separate_strands_uses_middle_top_margin_floor(
    tmp_path: Path,
    resolve_overlaps: bool,
) -> None:
    shared_args = [
        "--separate_strands",
        "--show_labels",
        "none",
        "--scale_style",
        "ruler",
        "--ruler_on_axis",
    ]
    maybe_resolve = ["--resolve_overlaps"] if resolve_overlaps else []

    middle_returncode, middle_stdout, middle_stderr, middle_output_svg = _run_linear_with_gbks(
        tmp_path,
        [INPUT_HMMTDNA],
        ["--track_layout", "middle", *shared_args, *maybe_resolve],
        output_name=f"linear_middle_separate_top_margin_{'resolve' if resolve_overlaps else 'base'}",
    )
    above_returncode, above_stdout, above_stderr, above_output_svg = _run_linear_with_gbks(
        tmp_path,
        [INPUT_HMMTDNA],
        ["--track_layout", "above", *shared_args, *maybe_resolve],
        output_name=f"linear_above_separate_top_margin_{'resolve' if resolve_overlaps else 'base'}",
    )

    assert middle_returncode == 0, f"stdout={middle_stdout}\nstderr={middle_stderr}"
    assert above_returncode == 0, f"stdout={above_stdout}\nstderr={above_stderr}"

    middle_top_margin = _extract_min_absolute_feature_y(middle_output_svg.read_text(encoding="utf-8"))
    above_top_margin = _extract_min_absolute_feature_y(above_output_svg.read_text(encoding="utf-8"))
    assert above_top_margin == pytest.approx(middle_top_margin, abs=1e-6)


@pytest.mark.linear
@pytest.mark.parametrize(
    ("scale_style", "axis_args"),
    [
        ("bar", []),
        ("ruler", ["--ruler_on_axis"]),
    ],
)
@pytest.mark.parametrize("resolve_overlaps", [False, True])
def test_linear_above_non_stranded_uses_middle_top_margin_floor(
    tmp_path: Path,
    scale_style: str,
    axis_args: list[str],
    resolve_overlaps: bool,
) -> None:
    shared_args = [
        "--show_labels",
        "none",
        "--scale_style",
        scale_style,
        *axis_args,
    ]
    maybe_resolve = ["--resolve_overlaps"] if resolve_overlaps else []

    middle_returncode, middle_stdout, middle_stderr, middle_output_svg = _run_linear_with_gbks(
        tmp_path,
        [INPUT_HMMTDNA],
        ["--track_layout", "middle", *shared_args, *maybe_resolve],
        output_name=f"linear_middle_non_stranded_top_margin_{scale_style}_{'resolve' if resolve_overlaps else 'base'}",
    )
    above_returncode, above_stdout, above_stderr, above_output_svg = _run_linear_with_gbks(
        tmp_path,
        [INPUT_HMMTDNA],
        ["--track_layout", "above", *shared_args, *maybe_resolve],
        output_name=f"linear_above_non_stranded_top_margin_{scale_style}_{'resolve' if resolve_overlaps else 'base'}",
    )

    assert middle_returncode == 0, f"stdout={middle_stdout}\nstderr={middle_stderr}"
    assert above_returncode == 0, f"stdout={above_stdout}\nstderr={above_stderr}"

    middle_top_margin = _extract_min_absolute_feature_y(middle_output_svg.read_text(encoding="utf-8"))
    above_top_margin = _extract_min_absolute_feature_y(above_output_svg.read_text(encoding="utf-8"))
    assert above_top_margin == pytest.approx(middle_top_margin, abs=1e-6)


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
def test_linear_default_axis_color_stays_legacy_without_ruler_on_axis(tmp_path: Path) -> None:
    returncode, stdout, stderr, output_svg = _run_linear(
        tmp_path,
        ["--track_layout", "above", "--scale_style", "ruler", "--scale_interval", "50000"],
    )
    assert returncode == 0, f"stdout={stdout}\nstderr={stderr}"
    svg_content = output_svg.read_text(encoding="utf-8")
    axis_match = re.search(
        r'<line fill="none" stroke="([^"]+)" stroke-width="([^"]+)"[^>]*x1="0"[^>]*x2="2000\.0"[^>]*y1="0"[^>]*y2="0"',
        svg_content,
    )
    assert axis_match is not None
    assert axis_match.group(1) == "lightgray"


@pytest.mark.linear
def test_linear_ruler_tick_length_is_two_thirds_of_legacy_default() -> None:
    assert float(RULER_TICK_LENGTH) == pytest.approx(10.0 * (2.0 / 3.0))


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
    axis_groups = [g for g in group_matches if _count_ruler_ticks(g) > 0]
    assert len(axis_groups) == 3
    tick_counts = [_count_ruler_ticks(g) for g in axis_groups]
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
    tick_matches = re.finditer(
        r'<line[^>]*stroke="([^"]+)"[^>]*stroke-width="([^"]+)"[^>]*y2="(-?\d+(?:\.\d+)?)"',
        svg_content,
    )
    tick_match = next((match for match in tick_matches if _is_ruler_tick_y2(float(match.group(3)))), None)
    label_match = re.search(r'<text[^>]*fill="([^"]+)"[^>]*>[^<]*(?:bp|kbp|Mbp)</text>', svg_content)
    assert axis_match is not None
    assert tick_match is not None
    assert label_match is not None
    assert axis_match.group(1) == "dimgray"
    assert tick_match.group(1) == axis_match.group(1)
    assert float(tick_match.group(2)) == pytest.approx(float(axis_match.group(2)))
    assert label_match.group(1) == axis_match.group(1)


@pytest.mark.linear
def test_linear_ruler_on_axis_defaults_label_font_size_to_long_ruler_value(tmp_path: Path) -> None:
    returncode, stdout, stderr, output_svg = _run_linear(
        tmp_path,
        ["--track_layout", "above", "--scale_style", "ruler", "--ruler_on_axis", "--scale_interval", "50000"],
    )
    assert returncode == 0, f"stdout={stdout}\nstderr={stderr}"
    svg_content = output_svg.read_text(encoding="utf-8")
    ruler_text_tags = _extract_ruler_text_tags(svg_content)
    assert any('font-size="12.0"' in tag for tag in ruler_text_tags)


@pytest.mark.linear
def test_linear_ruler_bottom_defaults_label_font_size_to_long_ruler_value(tmp_path: Path) -> None:
    returncode, stdout, stderr, output_svg = _run_linear(
        tmp_path,
        ["--scale_style", "ruler", "--scale_interval", "50000"],
    )
    assert returncode == 0, f"stdout={stdout}\nstderr={stderr}"
    svg_content = output_svg.read_text(encoding="utf-8")
    length_bar_match = re.search(r'<g id="length_bar"[^>]*>(.*?)</g>', svg_content)
    assert length_bar_match is not None
    ruler_text_tags = _extract_ruler_text_tags(length_bar_match.group(1))
    assert any('font-size="12.0"' in tag for tag in ruler_text_tags)


@pytest.mark.linear
def test_linear_ruler_uses_scale_font_size_when_ruler_size_unset(tmp_path: Path) -> None:
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
            "--scale_font_size",
            "19",
        ],
    )
    assert returncode == 0, f"stdout={stdout}\nstderr={stderr}"
    svg_content = output_svg.read_text(encoding="utf-8")
    ruler_text_tags = _extract_ruler_text_tags(svg_content)
    assert any('font-size="19.0"' in tag for tag in ruler_text_tags)


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
    ruler_text_tags = _extract_ruler_text_tags(svg_content)
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
    ruler_text_tags = _extract_ruler_text_tags(length_bar_match.group(1))
    assert any('fill="teal"' in tag and 'font-size="21.0"' in tag for tag in ruler_text_tags)
