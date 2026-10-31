"""Length options read a bp count with an optional unit; coordinates stay integers."""

from __future__ import annotations

import argparse

import pytest

from gbdraw.circular import _get_args as circular_args
from gbdraw.core.numeric import argparse_length_bp, argparse_scale_interval, parse_length_bp
from gbdraw.linear import _get_args as linear_args


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("50000", 50_000),
        ("50000bp", 50_000),
        ("50 kb", 50_000),
        ("50k", 50_000),
        ("50kb", 50_000),
        ("50kbp", 50_000),
        ("50KB", 50_000),
        ("0.5Mb", 500_000),
        ("0.5m", 500_000),
        ("1.5mbp", 1_500_000),
        ("1.5kb", 1_500),
        (".5kb", 500),
        (" 2kb ", 2_000),
    ],
)
def test_parse_length_bp_reads_a_count_with_a_unit(text: str, expected: int) -> None:
    assert parse_length_bp(text) == expected


@pytest.mark.parametrize(
    "text",
    [
        "",
        "0",
        "0kb",
        "-5",
        "+5",
        "-1kb",
        "1e3",
        "0.0004kb",  # less than one base
        "1.5",  # a base count is whole
        "1.5bp",
        "5 k b",
        "5gb",
        "kb",
        "1_000",
        "١٠",  # another script's digits
        "٥kb",
        "nan",
        "inf",
        "9007199254740992",  # a JavaScript reader would not keep it
        "1" + "0" * 5000,
    ],
)
def test_parse_length_bp_refuses_anything_else(text: str) -> None:
    assert parse_length_bp(text) is None
    with pytest.raises(argparse.ArgumentTypeError, match="invalid length"):
        argparse_length_bp(text)


@pytest.mark.parametrize("value", [True, False, None, 5, 5.0])
def test_parse_length_bp_reads_text_only(value: object) -> None:
    assert parse_length_bp(value) is None  # type: ignore[arg-type]


def test_scale_interval_keeps_the_integer_spellings_a_legacy_session_replays() -> None:
    assert argparse_scale_interval("5kb") == 5_000
    assert argparse_scale_interval("5000") == 5_000
    assert argparse_scale_interval("+5") == 5
    # The caller decides what a value <= 0 means (D-04 in gbdraw/linear.py).
    assert argparse_scale_interval("0") == 0
    assert argparse_scale_interval("-1") == -1
    with pytest.raises(argparse.ArgumentTypeError, match="invalid length"):
        argparse_scale_interval("5x")


@pytest.mark.parametrize("get_args", [linear_args, circular_args])
def test_both_modes_read_scale_interval_with_a_unit(get_args) -> None:
    assert get_args(["--gbk", "a.gb", "--scale_interval", "5kb"]).scale_interval == 5_000
    assert get_args(["--gbk", "a.gb", "--scale_interval", "5000"]).scale_interval == 5_000


@pytest.mark.parametrize("get_args", [linear_args, circular_args])
def test_scale_interval_at_most_zero_stays_an_error_on_the_command_line(get_args, capsys) -> None:
    with pytest.raises(SystemExit):
        get_args(["--gbk", "a.gb", "--scale_interval", "0"])
    assert "--scale_interval must be > 0" in capsys.readouterr().err


def test_wrap_length_is_auto_or_a_length_and_absent_by_default() -> None:
    assert linear_args(["--gbk", "a.gb"]).wrap_length is None
    assert linear_args(["--gbk", "a.gb", "--wrap_length", "auto"]).wrap_length == "auto"
    assert linear_args(["--gbk", "a.gb", "--wrap_length", "AUTO"]).wrap_length == "auto"
    assert linear_args(["--gbk", "a.gb", "--wrap_length", "50kb"]).wrap_length == 50_000


@pytest.mark.parametrize("text", ["0", "none", "-5kb", "1e3"])
def test_wrap_length_refuses_other_text(text: str, capsys) -> None:
    with pytest.raises(SystemExit):
        linear_args(["--gbk", "a.gb", "--wrap_length", text])
    assert "--wrap_length" in capsys.readouterr().err


def test_region_coordinates_stay_integers() -> None:
    from gbdraw.io.regions import parse_region_spec

    assert parse_region_spec("rec:1000-5000").start == 1000
    with pytest.raises(ValueError, match="Invalid region spec"):
        parse_region_spec("rec:1kb-5kb")
