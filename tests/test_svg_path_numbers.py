"""Every SVG path writes its numbers with 3 decimals in the shortest form (perf-audit D-13 D2)."""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from gbdraw.svg.path_data import path_line_to, path_number, path_numbers

# A number in path data: no exponent, no trailing zero, no "-0", at most 3 decimals.
_PATH_DATA_NUMBER = re.compile(r"^-?(?:0|[1-9]\d*)(?:\.\d{0,2}[1-9])?$")
_PATH_D = re.compile(r'<path\b[^>]*?\sd="([^"]*)"')
_TOKEN = re.compile(r"[-+0-9.eE]+")


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        (0, "0"),
        (0.0, "0"),
        (-0.0, "0"),
        (-0.0004, "0"),
        (0.0005, "0.001"),
        (12.0, "12"),
        (100.0, "100"),
        (1.5, "1.5"),
        (-1.25, "-1.25"),
        (3.14159265, "3.142"),
        (6.123233995736766e-17, "0"),
        (1e-05, "0"),
        (123456789.98765, "123456789.988"),
        (1e20, "100000000000000000000"),
    ],
)
def test_path_number_writes_three_decimals_in_shortest_form(value: float, expected: str) -> None:
    assert path_number(value) == expected
    assert _PATH_DATA_NUMBER.match(expected)


def test_array_and_polyline_forms_use_the_scalar_formatter() -> None:
    values = [0.1 + 0.2, -0.0, 2.5000001, 7]
    assert path_numbers(values) == [path_number(value) for value in values]
    assert path_line_to([1.23456, 2.0], [-0.00001, 3.9999]) == "L1.235 0L2 4"


def _annotation_table(directory: Path) -> Path:
    table = directory / "annotations.tsv"
    table.write_text(
        "set_id\tid\tmark\tstart\tend\tlabel\n"
        "regions\tline1\tline\t1000\t9000\tLine\n"
        "regions\tband1\tband\t12000\t30000\tBand\n"
        "regions\tbracket1\tbracket\t40000\t52000\tBracket\n",
        encoding="utf-8",
    )
    return table


_CASES = {
    "circular": (
        "circular",
        ["MjeNMV.gb"],
        [],
        ["--labels", "both", "--gc_content_mode", "percent", "--depth_track", "MjeNMV.DRR271272.depth.tsv"],
    ),
    "circular_deviation_skew": ("circular", ["MjeNMV.gb"], [], ["--separate_strands"]),
    "linear": (
        "linear",
        ["MjeNMV.gb", "MelaMJNV.gb"],
        ["MjeNMV.MelaMJNV.tblastx.out"],
        ["--gc", "--skew", "--show_labels", "all", "--depth_track", "MjeNMV.DRR271272.depth.tsv", "MelaMJNV.DRR394934.depth.tsv"],
    ),
    "linear_percent": ("linear", ["MjeNMV.gb"], [], ["--gc", "--gc_content_mode", "percent"]),
}


@pytest.mark.parametrize("case", sorted(_CASES))
def test_every_drawn_path_writes_three_decimals(case, tmp_path, find_test_input, gbdraw_runner) -> None:
    mode, gbks, blasts, args = _CASES[case]
    inputs = [find_test_input(name) for name in gbks + blasts]
    resolved_args = [str(find_test_input(arg)) if arg.endswith(".depth.tsv") else arg for arg in args]
    if not all(inputs) or "None" in resolved_args:
        pytest.skip("test inputs not found")
    if len(gbks) == 1:  # the table names no record
        resolved_args += ["--annotation_table", str(_annotation_table(tmp_path))]
    returncode, output, svg_path = gbdraw_runner.run(
        mode,
        inputs[: len(gbks)],
        case,
        tmp_path,
        blast_files=inputs[len(gbks):],
        extra_args=resolved_args,
    )
    assert returncode == 0, output
    svg = svg_path.read_text(encoding="utf-8")
    paths = _PATH_D.findall(svg)
    assert paths
    offending = sorted(
        {token for d in paths for token in _TOKEN.findall(d) if not _PATH_DATA_NUMBER.match(token)}
    )
    assert not offending, f"{len(offending)} path numbers not in 3-decimal shortest form: {offending[:10]}"
