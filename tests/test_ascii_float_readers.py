"""No Python reader of a float in user or file input accepts other digits or underscores.

``float()`` reads the decimal digits of every script (``float("١.٥") == 1.5``)
and digit-grouping underscores (``float("1_0") == 10.0``). Such text is read
as not a number, through each reader's own error path; everything else
``float()`` accepts still reads as before. The one owner of the grammar is
``gbdraw.core.numeric`` (the float half of the integer rule in
``tests/test_ascii_integer_readers.py``).
"""

from __future__ import annotations

import argparse
import re
from collections.abc import Callable

import pytest

# Arabic-Indic one-point-five, fullwidth three, an ASCII and a Devanagari digit, and
# underscore groups.
REFUSED_TEXT = ("١.٥", "\uff13", "1\u0967", "1_0", "1_000.5")
# What ``float()`` keeps accepting: sign, surrounding whitespace (other spaces too, as
# ``float()`` strips them), exponent, bare fraction or point. The non-finite spellings
# are checked separately.
FULLWIDTH = "".join(chr(0xFF10 + digit) for digit in range(10))
ACCEPTED_TEXT = {
    "2": 2.0, "-2.5": -2.5, "+.5": 0.5, "5.": 5.0, " 1e1\t": 10.0, "\n1E+1 ": 10.0,
    "\x0b3\x0c": 3.0, "\u00a02\u2003": 2.0,
}


def test_ascii_float_keeps_what_float_reads_in_ascii() -> None:
    import math

    from gbdraw.core.numeric import argparse_float, ascii_float, parse_ascii_float

    for text, value in ACCEPTED_TEXT.items():
        assert ascii_float(text) == value == argparse_float(text) == parse_ascii_float(text), text
    for text in ("inf", "-Infinity", "+iNF"):
        assert math.isinf(ascii_float(text))
    assert math.isnan(ascii_float("nan"))
    assert ascii_float(3) == 3.0 and ascii_float(True) == 1.0  # not text: float() unchanged


@pytest.mark.parametrize("text", [*REFUSED_TEXT, "", "x", "0x10", "1.5\x1c"])
def test_ascii_float_refuses_other_characters_with_the_errors_of_float(text: str) -> None:
    from gbdraw.core.numeric import argparse_float, ascii_float, parse_ascii_float

    import re

    with pytest.raises(ValueError, match=f"^{re.escape(f'could not convert string to float: {text!r}')}$"):
        ascii_float(text)
    assert parse_ascii_float(text) is None
    with pytest.raises(argparse.ArgumentTypeError, match=f"^{re.escape(f'invalid float value: {text!r}')}$"):
        argparse_float(text)


def test_ascii_int_is_int_without_other_digits_or_underscores() -> None:
    import re

    from gbdraw.core.numeric import ascii_int

    assert [ascii_int(value) for value in ("12", " -3\t", "\u00a04 ", 7.9, True)] == [12, -3, 4, 7, 1]
    for text in ("١٢", "\uff13", "1_0", "x", ""):
        with pytest.raises(ValueError, match=f"^{re.escape(f'invalid literal for int() with base 10: {text!r}')}$"):
            ascii_int(text)


def test_no_argparse_option_reads_its_float_with_bare_float() -> None:
    import ast
    from pathlib import Path

    import gbdraw

    offenders = []
    for path in sorted(Path(gbdraw.__file__).parent.rglob("*.py")):
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            if isinstance(node, ast.keyword) and node.arg == "type" and getattr(node.value, "id", "") == "float":
                offenders.append(f"{path.name}:{node.value.lineno}")
    assert offenders == []


def _float_options(mode: str) -> list[tuple[str, int]]:
    """The options whose argparse type reads a float, with how many values each takes, from the parser itself."""

    import os
    import sys
    from unittest.mock import patch

    from gbdraw import circular, linear

    class _Captured(Exception):
        pass

    captured: list[argparse.ArgumentParser] = []

    def capture(parser: argparse.ArgumentParser, *_args: object) -> None:
        captured.append(parser)
        raise _Captured

    module = circular if mode == "circular" else linear
    with (
        patch.object(argparse.ArgumentParser, "parse_args", capture),
        patch.object(sys, "argv", ["gbdraw"]),
        patch.dict(os.environ, {"COLUMNS": "80"}),
        pytest.raises(_Captured),
    ):
        module._get_args([])
    return [
        (action.option_strings[-1], action.nargs if isinstance(action.nargs, int) else 1)
        for action in captured[0]._actions
        if getattr(action.type, "__name__", "") in {"float", "argparse_float"}
    ]


def _with_value(option: str, count: int, text: str) -> list[str]:
    """The words that give ``option`` the value ``text``; an option taking more values gets ``text`` last."""

    return [f"{option}={text}"] if count == 1 else [option, *["50"] * (count - 1), text]


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_every_cli_float_option_refuses_non_ascii_text_and_underscores(
    mode: str, capsys: pytest.CaptureFixture[str]
) -> None:
    from gbdraw import circular, linear

    get_args = circular._get_args if mode == "circular" else linear._get_args
    options = _float_options(mode)
    assert len(options) > 15
    for option, count in options:
        for text in ("١.٥", "1_0"):
            with pytest.raises(SystemExit) as raised:
                get_args(["--gbk", "dummy.gb", *_with_value(option, count, text)])
            assert raised.value.code == 2, (option, text)
            assert f"invalid float value: {text!r}" in capsys.readouterr().err, (option, text)


def test_cli_float_options_keep_what_float_read() -> None:
    from gbdraw import circular, linear

    for text, value in (("0.25", 0.25), (" +.5 ", 0.5), ("2.5e-1", 0.25), ("5.", 5.0)):
        assert circular._get_args(["--gbk", "dummy.gb", f"--gc_content_radius={text}"]).gc_content_radius == value
        assert linear._get_args(["--gbk", "dummy.gb", f"--comparison_height={text}"]).comparison_height == value


@pytest.mark.parametrize("text", REFUSED_TEXT)
def test_linear_track_axis_gap_refuses_non_ascii_text(text: str, capsys: pytest.CaptureFixture[str]) -> None:
    from gbdraw import linear

    for ok, value in (("1.5", 1.5), (" 2 ", 2.0), ("auto", None), ("1e1", 10.0)):
        assert linear._get_args(["--gbk", "dummy.gb", f"--track_axis_gap={ok}"]).track_axis_gap == value
    with pytest.raises(SystemExit):
        linear._get_args(["--gbk", "dummy.gb", f"--track_axis_gap={text}"])
    assert "track axis gap must be a non-negative number or 'auto'" in capsys.readouterr().err


@pytest.mark.parametrize("mode", ["circular", "linear"])
@pytest.mark.parametrize("text", ["١", "1_0"])
def test_depth_track_values_refuse_non_ascii_text(
    mode: str, text: str, capsys: pytest.CaptureFixture[str]
) -> None:
    from gbdraw import circular, linear

    get_args = circular._get_args if mode == "circular" else linear._get_args
    option = "--depth_track_tick_font_size"
    assert get_args(["--gbk", "dummy.gb", option, "8", "auto"]).depth_track_tick_font_size == ["8", "auto"]
    with pytest.raises(SystemExit):
        get_args(["--gbk", "dummy.gb", option, text])
    assert f"{option} values must be numbers or auto" in capsys.readouterr().err


def _outcome(read: Callable[[str], object], text: str) -> tuple[str, object]:
    """What a reader does with ``text``: the value it returns or the error type it raises."""

    try:
        return ("value", read(text))
    except Exception as exc:  # noqa: BLE001 - the error type is the outcome
        return ("raises", type(exc))


def _readers() -> dict[str, Callable[[str], object]]:
    from gbdraw.analysis.depth_tracks import _expand_track_float_metadata
    from gbdraw.annotations.io import _number
    from gbdraw.annotations.models import _finite
    from gbdraw.api.record_planning import _optional_positive_float
    from gbdraw.cli_utils.common import parse_arrow_head_length_ratio_arg, parse_arrow_shaft_width_ratio_arg
    from gbdraw.config.models.canvas import _positive_finite_float
    from gbdraw.config.models.objects import (
        _optional_auto_positive_float,
        _optional_depth_bound,
        _optional_gc_content_bound,
        _optional_gc_content_positive,
    )
    from gbdraw.core.fonts import _css_weight
    from gbdraw.definition_line_styles import (
        normalize_definition_line_font_size,
        normalize_definition_line_font_weight,
    )
    from gbdraw.session_io import _is_depth_text
    from gbdraw.session_request_codec import _expanded_depth_numeric_metadata
    from gbdraw.tracks.scalars import ScalarSpec
    from gbdraw.web_support.feature_catalog import _integer_or_none
    from gbdraw.web_support.orthogroup_metadata import _float_value

    return {
        "--arrow_head_length_ratio": parse_arrow_head_length_ratio_arg,
        "--arrow_shaft_width_ratio": parse_arrow_shaft_width_ratio_arg,
        "config depth bound": _optional_depth_bound,
        "config GC content bound": _optional_gc_content_bound,
        "config GC content tick": _optional_gc_content_positive,
        "config conservation ring width": lambda text: _optional_auto_positive_float(text, field_name="w"),
        "config comparison height": lambda text: _positive_finite_float(text, field_name="h"),
        "definition line font size": normalize_definition_line_font_size,
        "definition line font weight": normalize_definition_line_font_weight,
        "font weight": _css_weight,
        "Session depth table value": _is_depth_text,
        "depth track metadata": lambda text: _expand_track_float_metadata([text], track_count=1, field_name="f"),
        "request depth metadata": lambda text: _expanded_depth_numeric_metadata([text], track_count=1, field_name="f"),
        "depth track option": lambda text: _optional_positive_float(text, field_name="f"),
        "annotation table number": lambda text: _number(text, row_number=1, column="c"),
        "annotation number": lambda text: _finite(text, "n"),
        "track slot scalar": ScalarSpec.parse,
        "feature catalog coordinate": _integer_or_none,
        "orthogroup support": _float_value,
    }


@pytest.mark.parametrize("name", sorted(_readers()))
def test_float_readers_read_non_ascii_text_as_not_a_number(name: str) -> None:
    read = _readers()[name]
    number = {"--arrow_shaft_width_ratio": "0.25", "font weight": "700"}.get(name, "400")
    not_a_number = _outcome(read, "x")
    assert _outcome(read, number) != not_a_number, name
    other_digits = (number.translate(str.maketrans("0123456789", script)) for script in ("٠١٢٣٤٥٦٧٨٩", FULLWIDTH))
    for text in (*other_digits, f"{number[:-1]}_{number[-1]}"):
        assert _outcome(read, text) == not_a_number, (name, text)
    assert _outcome(read, f"\u00a0{number} ") == _outcome(read, number), name


@pytest.mark.parametrize("text", ["١", "1_0"])
def test_linear_canvas_track_axis_gap_reads_non_ascii_text_as_invalid(text: str) -> None:
    import copy

    from gbdraw.config.models.canvas import LinearCanvasConfig
    from gbdraw.config.toml import load_config_toml

    raw = load_config_toml("gbdraw.data", "config.toml")["canvas"]["linear"]
    for value, expected in (("2", 2.0), (text, None), ("x", None)):
        linear = copy.deepcopy(raw)
        linear["track_axis_gap"] = value
        assert LinearCanvasConfig.from_dict(linear).track_axis_gap == expected, value


def test_annotation_table_integer_coordinates_refuse_non_ascii_digits() -> None:
    """OV-1006: the annotation table read start and end with a bare ``int()``."""

    from gbdraw.annotations.io import _number
    from gbdraw.exceptions import ValidationError

    assert _number(" +12 ", row_number=1, column="start", integer=True) == 12
    for text in ("١٢", "1_2", "1.0"):
        with pytest.raises(ValidationError, match="expected a integer"):
            _number(text, row_number=1, column="start", integer=True)


ConfigPath = tuple[str | int, ...]


def _numeric_config_leaves(value: object, path: ConfigPath = ()) -> list[ConfigPath]:
    items = value.items() if isinstance(value, dict) else enumerate(value) if isinstance(value, list) else None
    if items is not None:
        return [leaf for key, item in items for leaf in _numeric_config_leaves(item, (*path, key))]
    return [path] if isinstance(value, (int, float)) and not isinstance(value, bool) else []


def _value_at(raw: object, path: ConfigPath) -> object:
    value = raw
    for key in path:
        value = value[key]  # type: ignore[index]
    return value


def _config_with(raw: dict[str, object], path: ConfigPath, value: object) -> tuple[str, object]:
    import copy

    from gbdraw.config.models import GbdrawConfig

    edited = copy.deepcopy(raw)
    _value_at(edited, path[:-1])[path[-1]] = value  # type: ignore[index]
    return _outcome(GbdrawConfig.from_dict, edited)


def test_config_numbers_read_non_ascii_text_as_not_a_number() -> None:
    """Review r1 #1: a Session's or request's ``config`` (``GbdrawConfig.from_dict``)
    read its numbers with a bare ``float()`` or ``int()``, so ``radius: "٣٩٠"``
    drew a radius of 390. Each number of the default config now reads other
    digits and underscores as it reads ``"x"``, and its ASCII text as before
    (a leaf that takes only a number refuses all text).
    """

    from gbdraw.config.toml import load_config_toml

    raw = load_config_toml("gbdraw.data", "config.toml")
    leaves = _numeric_config_leaves(raw)
    assert {("canvas", "circular", "radius"), ("canvas", "dpi"), ("objects", "sliding_window", "default", 0)} <= set(leaves)
    to_arabic_indic = str.maketrans("0123456789", "٠١٢٣٤٥٦٧٨٩")
    failures = []
    for path in leaves:
        number = _value_at(raw, path)
        text = repr(number)
        read = _config_with(raw, path, number)
        not_a_number = _config_with(raw, path, "x")
        if not_a_number == read:
            continue  # GbdrawConfig does not read this leaf (png_output).
        ascii_text = _config_with(raw, path, f" {text} ")
        if ascii_text != read and ascii_text[0] != not_a_number[0]:
            failures.append((path, text, "ASCII text read differently"))
        grouped = f"{text[0]}_{text[1:]}" if len(text) > 1 and text[1].isdigit() else None
        for edited in filter(None, (text.translate(to_arabic_indic), grouped)):
            outcome = _config_with(raw, path, edited)
            if outcome == read or outcome[0] != not_a_number[0]:
                failures.append((path, edited, outcome[0]))
    assert failures == []


NON_FINITE_TEXT = ("nan", "inf", "-inf", "Infinity", "-NaN")


@pytest.mark.parametrize("text", NON_FINITE_TEXT)
def test_argparse_float_refuses_non_finite_values(text: str) -> None:
    """OV-275: ``nan`` and ``inf`` reached the layout and ended in a raw traceback."""

    from gbdraw.core.numeric import argparse_float, parse_ascii_finite_float

    assert parse_ascii_finite_float(text) is None
    assert parse_ascii_finite_float(" 2.5 ") == 2.5
    with pytest.raises(argparse.ArgumentTypeError, match="must be a finite number"):
        argparse_float(text)


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_every_cli_float_option_refuses_non_finite_values(mode: str, capsys: pytest.CaptureFixture[str]) -> None:
    """OV-275: every float option shares one finite-float type, so none lets ``nan`` through."""

    from gbdraw import circular, linear

    get_args = circular._get_args if mode == "circular" else linear._get_args
    options = _float_options(mode)
    assert len(options) > 15
    for option, count in options:
        for text in ("nan", "inf", "-inf"):
            with pytest.raises(SystemExit) as raised:
                get_args(["--gbk", "dummy.gb", *_with_value(option, count, text)])
            assert raised.value.code == 2, (option, text)
            err = capsys.readouterr().err
            refusals = [f"argument {option}: invalid float value: {text!r}"]
            if count > 1 and text.startswith("-"):
                # Before Python 3.13 argparse reads a separate "-inf" word as an option, so one value is missing.
                refusals.append(f"argument {option}: expected {count} arguments")
            assert any(refusal in err for refusal in refusals), (option, text, err)


@pytest.mark.parametrize(
    ("mode", "argv", "message"),
    [
        ("linear", ["--track_axis_gap=nan"], "track axis gap must be a non-negative number or 'auto'"),
        ("linear", ["--track_axis_gap=inf"], "track axis gap must be a non-negative number or 'auto'"),
        ("circular", ["--arrow_head_length_ratio=nan"], "arrow head length ratio must be 'auto'"),
        ("circular", ["--depth_track_tick_font_size", "nan"], "values must be numbers or auto"),
        ("linear", ["--depth_track_height", "inf"], "values must be numbers or auto"),
        # A slot scalar read nan and failed later with an unrelated overlap error.
        ("circular", ["--circular_track_slot", "features:features@r=nan"], "must be a finite number"),
        ("circular", ["--circular_track_slot", "features:features@w=infpx"], "must be a finite number"),
        ("circular", ["--circular_track_slot", "features:features@r=-inf%"], "must be a finite number"),
    ],
)
def test_cli_text_number_options_refuse_non_finite_values(
    mode: str, argv: list[str], message: str, capsys: pytest.CaptureFixture[str]
) -> None:
    from gbdraw import circular, linear

    get_args = circular._get_args if mode == "circular" else linear._get_args
    with pytest.raises(SystemExit):
        get_args(["--gbk", "dummy.gb", *argv])
    assert message in capsys.readouterr().err


@pytest.mark.parametrize("value", [float("nan"), float("inf"), float("-inf")])
def test_config_override_floats_must_be_finite(value: float) -> None:
    """OV-275: ``config_overrides={"canvas.circular.radius": nan}`` ended in
    ``ValueError: bbox bounds must be finite``; the override validator owns the domain."""

    from gbdraw.config.modify import validate_config_overrides
    from gbdraw.exceptions import ValidationError

    assert validate_config_overrides({"canvas.circular.radius": 300}) == (("canvas.circular.radius", 300),)
    for path in ("canvas.circular.radius", "labels.unified_adjustment.outer_labels.x_radius_offset"):
        with pytest.raises(ValidationError, match=re.escape(path)) as raised:
            validate_config_overrides({path: value})
        assert raised.value.diagnostic == {"code": "INPUT_INVALID", "configPath": path, "reason": "FINITE"}


def test_python_draw_refuses_a_non_finite_config_override() -> None:
    from pathlib import Path

    import gbdraw
    from gbdraw.exceptions import ValidationError

    records = gbdraw.read_genbank(Path(__file__).parent / "test_inputs" / "NC_001416.gb")
    with pytest.raises(ValidationError, match="canvas.circular.radius"):
        gbdraw.draw_circular(
            records,
            options=gbdraw.CircularOptions(config_overrides={"canvas.circular.radius": float("nan")}),
        ).to_svg()
