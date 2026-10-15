"""Every Python reader of an integer in user or file input accepts ASCII digits only.

Python's ``\\d`` and ``str.isdigit()`` accept other Unicode digits, and ``int()``
reads the decimal ones (``int("١") == 1``). gbdraw writes ASCII, so a value with
other digits is hand-edited or foreign data and is refused with the reader's own
error. The one owner of the grammar is ``gbdraw.core.numeric``.
"""

from __future__ import annotations

import pytest

from gbdraw.exceptions import ValidationError

# Arabic-Indic one, superscript two (``isdigit()`` is true, ``int()`` raises), fullwidth three.
NON_ASCII_DIGITS = ("١", "²", "３")


def test_ascii_integer_parsers_accept_only_ascii_decimal_text() -> None:
    from gbdraw.core.numeric import parse_ascii_integer, parse_ascii_nonnegative_integer

    assert parse_ascii_nonnegative_integer("0") == 0
    assert parse_ascii_nonnegative_integer("0012") == 12
    assert parse_ascii_integer("+7") == 7
    assert parse_ascii_integer("-7") == -7
    for bad in ("", " 1", "1 ", "1\n", "+1", "1_0", "1.0", "-", "١", "²", "３", "1２"):
        assert parse_ascii_nonnegative_integer(bad) is None, bad
    for bad in ("", " 1", "1\n", "+-1", "1_0", "١", "-²"):
        assert parse_ascii_integer(bad) is None, bad
    assert parse_ascii_nonnegative_integer("9" * 5000) is None


@pytest.mark.parametrize("digit", NON_ASCII_DIGITS)
def test_record_selector_refuses_non_ascii_digits(digit: str) -> None:
    from gbdraw.io.record_select import parse_record_selector

    assert parse_record_selector("#2").record_index == 1  # type: ignore[union-attr]
    with pytest.raises(ValueError, match="Invalid record selector"):
        parse_record_selector(f"#{digit}")


@pytest.mark.parametrize("digit", NON_ASCII_DIGITS)
def test_region_spec_refuses_non_ascii_digits(digit: str) -> None:
    from gbdraw.io.regions import parse_region_spec

    assert parse_region_spec("#2:1-5").record_index == 1
    with pytest.raises(ValueError, match="Invalid record index in region spec"):
        parse_region_spec(f"#{digit}:1-5")
    with pytest.raises(ValueError, match="Invalid region spec"):
        parse_region_spec(f"#1:{digit}-5")
    with pytest.raises(ValueError, match="Invalid region spec"):
        parse_region_spec(f"#1:1-{digit}")


@pytest.mark.parametrize("digit", NON_ASCII_DIGITS)
@pytest.mark.parametrize("compatibility", ("circular", "linear"))
def test_record_row_position_refuses_non_ascii_digits(digit: str, compatibility: str) -> None:
    from gbdraw.layout.record_placement import parse_record_row_position

    assert parse_record_row_position("#2@3", _compatibility=compatibility) == ("#2", 3)  # type: ignore[arg-type]
    with pytest.raises(ValidationError, match="positive integer row"):
        parse_record_row_position(f"#1@{digit}", _compatibility=compatibility)  # type: ignore[arg-type]
    if compatibility == "circular":  # linear resolves the selector against the records later
        with pytest.raises(ValidationError, match="Invalid record selector"):
            parse_record_row_position(f"#{digit}@1", _compatibility=compatibility)  # type: ignore[arg-type]


@pytest.mark.parametrize("digit", NON_ASCII_DIGITS)
def test_linear_record_placement_selector_refuses_non_ascii_digits(digit: str) -> None:
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    from gbdraw.layout.record_placement import resolve_record_row_positions

    records = [SeqRecord(Seq("ACGT"), id="r0")]
    assert resolve_record_row_positions(records, ("#1@1",)) is not None
    with pytest.raises(ValidationError, match="is invalid"):
        resolve_record_row_positions(records, (f"#{digit}@1",))


@pytest.mark.parametrize("digit", NON_ASCII_DIGITS)
def test_linear_comparison_selector_refuses_non_ascii_digits(digit: str) -> None:
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    from gbdraw.api.record_planning import _resolve_linear_comparison_selector

    records = [SeqRecord(Seq("ACGT"), id=f"r{i}") for i in range(2)]
    kwargs = {"table_path": "t.tsv", "row_number": 1, "column": "query"}
    assert _resolve_linear_comparison_selector(records, "#2", **kwargs) == 1
    with pytest.raises(ValidationError, match="invalid record selector"):
        _resolve_linear_comparison_selector(records, f"#{digit}", **kwargs)


@pytest.mark.parametrize("digit", NON_ASCII_DIGITS)
def test_interactive_svg_integers_and_id_suffixes_are_ascii(digit: str) -> None:
    from gbdraw.render import interactive_svg as svg

    assert svg._strict_nonnegative_integer("12") == 12
    assert svg._strict_nonnegative_integer(" 7 ") == 7
    assert svg._strict_nonnegative_integer(digit) is None
    assert svg._FEATURE_RECORD_SUFFIX_RE.search("f_record_2") is not None
    assert svg._FEATURE_RECORD_SUFFIX_RE.search(f"f_record_{digit}") is None
    assert svg._FEATURE_ELEMENT_SUFFIX_RE.search("f__part2") is not None
    assert svg._FEATURE_ELEMENT_SUFFIX_RE.search(f"f__part{digit}") is None
    assert svg._FEATURE_CONNECTOR_SUFFIX_RE.search(f"f__line{digit}") is None


@pytest.mark.parametrize("digit", NON_ASCII_DIGITS)
def test_orthogroup_metadata_alias_refuses_non_ascii_digits(digit: str) -> None:
    from gbdraw.web_support.orthogroup_metadata import _consistent_nonnegative_alias

    assert _consistent_nonnegative_alias({"a": "3"}, ["a"]) == (3, True)
    assert _consistent_nonnegative_alias({"a": digit}, ["a"]) == (None, False)


@pytest.mark.parametrize("digit", NON_ASCII_DIGITS)
def test_session_slot_and_depth_position_refuse_non_ascii_digits(digit: str) -> None:
    from gbdraw import session_io

    assert session_io._parse_slot("files.seqs[1].gb") == ["files", "seqs", 1, "gb"]
    with pytest.raises(ValidationError, match="Invalid session slot"):
        session_io._parse_slot(f"files.seqs[{digit}].gb")
    assert session_io._parse_positive_safe_integer("5") == 5
    assert session_io._parse_positive_safe_integer("+5") == 5
    assert session_io._parse_positive_safe_integer("-5") is None
    assert session_io._parse_positive_safe_integer(digit) is None


@pytest.mark.parametrize("digit", NON_ASCII_DIGITS)
def test_track_index_identity_refuses_non_ascii_digits(digit: str) -> None:
    from gbdraw.tracks.parsing import parse_nonnegative_integer

    assert parse_nonnegative_integer("4", field_name="track_index") == 4
    assert parse_nonnegative_integer("+4", field_name="track_index") == 4
    with pytest.raises(ValueError, match="track_index must be a nonnegative integer"):
        parse_nonnegative_integer(digit, field_name="track_index")


CLI_INTEGER_OPTIONS = {
    "circular": (
        ("--scale_interval", "scale_interval"), ("--window", "window"), ("-s", "step"),
        ("--feature_overlap_tolerance_bp", "feature_overlap_tolerance_bp"),
        ("--display_start_coordinate", "display_start_coordinate"),
        ("--depth_window", "depth_window"), ("--circular_track_axis_index", "circular_track_axis_index"),
    ),
    "linear": (
        ("--scale_interval", "scale_interval"), ("--alignment_length", "alignment_length"),
        ("--collinear_min_anchors", "collinear_min_anchors"), ("--depth_step", "depth_step"),
        ("--linear_track_axis_index", "linear_track_axis_index"),
    ),
}


def _cli_get_args(mode: str):
    from gbdraw import circular, linear

    return circular._get_args if mode == "circular" else linear._get_args


@pytest.mark.parametrize("mode", sorted(CLI_INTEGER_OPTIONS))
@pytest.mark.parametrize("text", ["٣", "²", "1_0", "1.5", "", "9" * 5000])
def test_cli_integer_options_refuse_non_ascii_digits_and_underscores(
    mode: str, text: str, capsys: pytest.CaptureFixture[str]
) -> None:
    for option, _dest in CLI_INTEGER_OPTIONS[mode]:
        with pytest.raises(SystemExit) as raised:
            _cli_get_args(mode)(["--gbk", "dummy.gb", option, text])
        assert raised.value.code == 2
        assert f"invalid int value: {text!r}" in capsys.readouterr().err


@pytest.mark.parametrize("mode", sorted(CLI_INTEGER_OPTIONS))
def test_cli_integer_options_keep_signs_and_surrounding_spaces(mode: str) -> None:
    for option, dest in CLI_INTEGER_OPTIONS[mode]:
        if dest.endswith("_axis_index"):  # needs a track order; the refusal test above covers it
            continue
        for text, value in (("4", 4), ("+4", 4), (" 4 ", 4), ("007", 7)):
            args = _cli_get_args(mode)(["--gbk", "dummy.gb", f"{option}={text}"])
            assert getattr(args, dest) == value, (option, text)


@pytest.mark.parametrize("option", ["--losatp_max_target_seqs", "--losatp_member_max_hits"])
@pytest.mark.parametrize("text", ["٣", "²", "1_0", "0", "-2", "x", "9" * 5000])
def test_linear_losatp_limits_refuse_non_ascii_digits_and_underscores(
    option: str, text: str, capsys: pytest.CaptureFixture[str]
) -> None:
    from gbdraw import linear

    dest = option.lstrip("-")
    for ok, value in (("4", 4), (" +4 ", 4), ("none", None), ("AUTO", None), ("", None)):
        assert getattr(linear._get_args(["--gbk", "dummy.gb", option, ok]), dest) == value
    with pytest.raises(SystemExit):
        linear._get_args(["--gbk", "dummy.gb", option, text])
    assert "must be a positive integer or 'none'" in capsys.readouterr().err


@pytest.mark.parametrize("digit", NON_ASCII_DIGITS)
def test_track_slot_z_refuses_non_ascii_digits_underscores_and_long_digits(digit: str) -> None:
    from gbdraw.tracks.circular import CircularTrackSlotParseError, parse_circular_track_slot
    from gbdraw.tracks.linear import LinearTrackSlotParseError, parse_linear_track_slot

    for parse, error in (
        (parse_circular_track_slot, CircularTrackSlotParseError),
        (parse_linear_track_slot, LinearTrackSlotParseError),
    ):
        assert parse("a:gc_content@z=3").z == 3
        assert parse("a:gc_content@z=-3").z == -3
        for bad in (digit, "1_0", "9" * 5000):
            with pytest.raises(error, match="z must be an integer"):
                parse(f"a:gc_content@z={bad}")


def test_crop_genbank_coordinates_refuse_non_ascii_digits(
    monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    from gbdraw import crop_genbank

    base = ["crop", "-i", "in.gb", "-o", "out.gb"]
    monkeypatch.setattr("sys.argv", [*base, "-s", "10", "-e", "20"])
    args = crop_genbank._get_args()
    assert (args.start, args.end) == (10, 20)
    monkeypatch.setattr("sys.argv", [*base, "-s", "١٠", "-e", "20"])
    with pytest.raises(SystemExit):
        crop_genbank._get_args()
    assert "invalid int value: '١٠'" in capsys.readouterr().err


def test_no_argparse_option_reads_its_integer_with_bare_int() -> None:
    import ast
    from pathlib import Path

    import gbdraw

    offenders = []
    for path in sorted(Path(gbdraw.__file__).parent.rglob("*.py")):
        for node in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            if isinstance(node, ast.keyword) and node.arg == "type" and getattr(node.value, "id", "") == "int":
                offenders.append(f"{path.name}:{node.value.lineno}")
    assert offenders == []


def test_argparse_integer_keeps_the_signs_and_spaces_of_int() -> None:
    import argparse

    from gbdraw.core.numeric import argparse_integer

    assert [argparse_integer(text) for text in ("5", "-5", "+5", " 5 ", "\t5\n", "007")] == [5, -5, 5, 5, 5, 7]
    for bad in ("", " ", "٣", "²", "1_0", "1.0", "0x1", "--1"):
        with pytest.raises(argparse.ArgumentTypeError, match=f"invalid int value: {bad!r}"):
            argparse_integer(bad)


@pytest.mark.parametrize("digit", NON_ASCII_DIGITS)
def test_losat_positive_int_option_refuses_non_ascii_digits(digit: str) -> None:
    import argparse

    from gbdraw.cli_utils.losat_output import parse_positive_int

    assert parse_positive_int(" 4 ") == 4
    with pytest.raises(argparse.ArgumentTypeError, match="positive integer"):
        parse_positive_int(digit)
    with pytest.raises(argparse.ArgumentTypeError, match="positive integer"):
        parse_positive_int("1_0")


@pytest.mark.parametrize("digit", NON_ASCII_DIGITS)
def test_orthogroup_int_or_none_refuses_non_ascii_digits(digit: str) -> None:
    from gbdraw.web_support.orthogroup_metadata import _int_or_none

    assert _int_or_none("3") == 3 and _int_or_none(3.0) == 3 and _int_or_none(None) is None
    assert _int_or_none(digit) is None
    assert _int_or_none("1_0") is None


@pytest.mark.parametrize("digit", NON_ASCII_DIGITS)
def test_comparison_feature_indexes_refuse_non_ascii_digits(digit: str) -> None:
    import pandas as pd

    from gbdraw.api.record_planning import ComparisonRecordBinding, project_comparison_frame

    binding = ComparisonRecordBinding({1: "a"}, {1: "a"}, {1: "b"}, 100)

    def frame(index: str) -> pd.DataFrame:
        return pd.DataFrame(
            {
                "query_feature_index": [index], "query_feature_svg_id": ["a"],
                "subject_feature_index": ["1"], "subject_feature_svg_id": ["a"],
            }
        )

    assert project_comparison_frame(frame("1"), binding, binding) is not None
    with pytest.raises(ValidationError, match="feature indexes must be integers"):
        project_comparison_frame(frame(digit), binding, binding)
    with pytest.raises(ValidationError, match="feature indexes must be integers"):
        project_comparison_frame(frame("1_0"), binding, binding)


def test_over_long_digit_strings_are_refused_by_their_readers_not_raised() -> None:
    """OV-961: a string over int()'s 4300-digit limit used to escape as a raw ValueError."""

    from gbdraw import session_io
    from gbdraw.config.models.root import _is_font_weight
    from gbdraw.io.regions import parse_region_spec
    from gbdraw.web_support.feature_catalog import _match_endpoint_indexes

    long = "9" * 5000
    with pytest.raises(ValidationError, match="Invalid session slot"):
        session_io._parse_slot(f"files.seqs[{long}].gb")
    with pytest.raises(ValueError, match="Invalid region spec"):
        parse_region_spec(f"#1:{long}-5")
    with pytest.raises(ValueError, match="Invalid region spec"):
        parse_region_spec(f"#1:5-{long}")
    assert _is_font_weight("700") and not _is_font_weight("٧٠٠") and not _is_font_weight(long)
    from gbdraw.exceptions import GbdrawError

    with pytest.raises(GbdrawError, match="invalid source feature index"):
        _match_endpoint_indexes({"i": long}, ("i",), description="match")
    with pytest.raises(GbdrawError, match="invalid source feature index"):
        _match_endpoint_indexes({"i": "١"}, ("i",), description="match")
    assert _match_endpoint_indexes({"i": "12"}, ("i",), description="match") == [12]


def test_schema_three_promotion_ignores_a_non_ascii_or_over_long_disambiguated_index() -> None:
    from gbdraw.web_support.feature_catalog import promote_legacy_feature_catalog

    ids = ["s~7", "s~١", "s~²", "s~" + "9" * 5000]
    features = [{"stableFeatureId": "s", "biologicalFeatureId": value} for value in ids]
    promoted = promote_legacy_feature_catalog({"schema": 3, "items": [{"biologicalFeatures": features}]})
    assert [f.get("sourceFeatureIndex") for f in promoted["items"][0]["biologicalFeatures"]] == [7, None, None, None]


@pytest.mark.parametrize("digit", NON_ASCII_DIGITS)
def test_cds_reading_frame_refuses_non_ascii_digits_and_long_digits(digit: str) -> None:
    from gbdraw.core.sequence import _cds_reading_frame

    assert _cds_reading_frame({"codon_start": ["2"]}) == 2
    assert _cds_reading_frame({"phase": ["1"]}) == 2
    for qualifiers, name in (({"codon_start": [digit]}, "codon_start"), ({"phase": [digit]}, "phase"),
                             ({"codon_start": ["9" * 5000]}, "codon_start")):
        with pytest.raises(ValueError, match=f"{name} is invalid"):
            _cds_reading_frame(qualifiers)


