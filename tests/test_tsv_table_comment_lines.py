"""Whole-line comments only in the label and visibility TSV tables (OV-20).

A ``#`` that is not the first non-blank character of a line is part of the
cell value. The Web writes such values (``foo#bar``, ``Gene #1``) to these
tables, so a reader that cuts a line at the first ``#`` rejects the table or
silently shortens the value.
"""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Callable

import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

import gbdraw.circular as circular_cli_module
from gbdraw.exceptions import ParseError
from gbdraw.features.visibility import read_feature_visibility_file
from gbdraw.labels.filtering import read_filter_list_file, read_label_override_file

# (reader, columns, a row whose last cell contains an inline "#")
_READERS: list[Any] = [
    pytest.param(
        read_feature_visibility_file,
        ["record_id", "feature_type", "qualifier", "value", "action"],
        ["*", "CDS", "product", "foo#bar", "off"],
        id="feature-visibility",
    ),
    pytest.param(
        read_label_override_file,
        ["record_id", "feature_type", "qualifier", "value", "label_text"],
        ["*", "CDS", "product", "kinase", "Gene #1"],
        id="label-override",
    ),
    pytest.param(
        read_filter_list_file,
        ["feature_type", "qualifier", "keyword"],
        ["CDS", "product", "protein #1"],
        id="label-filter-list",
    ),
]


def _write(path: Path, lines: list[str]) -> str:
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return str(path)


def _without_hash(row: list[str]) -> list[str]:
    return [cell.replace("#", "") for cell in row]


@pytest.mark.parametrize(("reader", "columns", "row"), _READERS)
def test_inline_hash_is_part_of_every_cell_value(
    tmp_path: Path, reader: Callable[[str], Any], columns: list[str], row: list[str]
) -> None:
    variants = [
        row,
        # A cell that starts with "#" and a "#" in the middle of an earlier cell.
        [*row[:-1], "#" + row[-1]],
        [row[0], row[1], row[2] + "#x", *row[3:]],
    ]
    path = _write(tmp_path / "table.tsv", ["\t".join(variant) for variant in variants])

    frame = reader(path)

    assert list(frame.columns) == columns
    assert frame.values.tolist() == variants


@pytest.mark.parametrize(("reader", "columns", "row"), _READERS)
def test_whole_line_comments_are_still_skipped(
    tmp_path: Path, reader: Callable[[str], Any], columns: list[str], row: list[str]
) -> None:
    row = _without_hash(row)
    path = _write(
        tmp_path / "table.tsv",
        [
            "# " + "\t".join(columns),
            "\t".join(row),
            "   # an indented comment\twith\ttabs",
            "#" + "\t".join(row),
            "",
            "\t".join(row),
        ],
    )

    frame = reader(path)

    assert frame.values.tolist() == [row, row]


@pytest.mark.parametrize(("reader", "columns", "row"), _READERS)
def test_utf8_bom_before_a_first_line_comment_or_row_is_ignored(
    tmp_path: Path, reader: Callable[[str], Any], columns: list[str], row: list[str]
) -> None:
    row = _without_hash(row)
    for first_line in ("# " + "\t".join(columns), "\t".join(row)):
        path = tmp_path / "bom.tsv"
        path.write_bytes(
            b"\xef\xbb\xbf" + (first_line + "\n" + "\t".join(row) + "\n").encode("utf-8")
        )

        frame = reader(str(path))

        assert frame.values.tolist()[-1] == row
        assert frame.values.tolist()[0] == row


@pytest.mark.parametrize(("reader", "columns", "row"), _READERS)
def test_comment_lines_do_not_shift_parse_error_line_numbers(
    tmp_path: Path, reader: Callable[[str], Any], columns: list[str], row: list[str]
) -> None:
    row = _without_hash(row)
    path = _write(
        tmp_path / "table.tsv",
        ["# comment", "\t".join(row), "# comment", "\t".join([*row, "extra"])],
    )

    with pytest.raises(ParseError, match=r"line 4\b"):
        reader(path)


def test_table_readers_do_not_use_the_pandas_comment_option() -> None:
    """Guard: ``comment=`` cuts a value at its first ``#``; the shared helper must own comments."""
    package = Path(__file__).resolve().parents[1] / "gbdraw"
    offenders = [
        f"{path.relative_to(package.parent)}:{number}"
        for path in sorted(package.rglob("*.py"))
        for number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1)
        if re.search(r"\bcomment\s*=\s*['\"]", line)
    ]
    assert offenders == []


def _write_genbank(path: Path) -> None:
    record = SeqRecord(
        Seq("ACGT" * 150),
        id="rec1",
        name="rec1",
        description="inline hash fixture",
        annotations={"molecule_type": "DNA", "topology": "circular"},
    )
    record.features = [
        SeqFeature(FeatureLocation(10, 190, strand=1), type="CDS",
                   qualifiers={"product": ["kinase #1"]}),
        SeqFeature(FeatureLocation(250, 430, strand=-1), type="CDS",
                   qualifiers={"product": ["other protein"]}),
    ]
    SeqIO.write(record, path, "genbank")


def _run_circular(tmp_path: Path, *options: str) -> str:
    gbk = tmp_path / "rec1.gb"
    _write_genbank(gbk)
    circular_cli_module.circular_main(
        ["--gbk", str(gbk), "-f", "svg", "-o", str(tmp_path / "out"), *options]
    )
    return (tmp_path / "out.svg").read_text(encoding="utf-8")


def test_cli_visibility_row_with_inline_hash_hides_the_matching_feature(tmp_path: Path) -> None:
    table = tmp_path / "visibility.tsv"
    table.write_text("*\tCDS\tproduct\tkinase #1\toff\n", encoding="utf-8")

    svg = _run_circular(tmp_path, "--feature_visibility_table", str(table))

    assert len(set(re.findall(r'data-gbdraw-feature-id="([^"]+)"', svg))) == 1


def test_cli_label_table_row_keeps_inline_hash_in_the_label_text(tmp_path: Path) -> None:
    table = tmp_path / "labels.tsv"
    table.write_text("*\tCDS\tproduct\tkinase\tGene #1\n", encoding="utf-8")

    svg = _run_circular(tmp_path, "--labels", "--label_table", str(table))

    root = ET.fromstring(svg)
    texts = ["".join(node.itertext()) for node in root.iter("{http://www.w3.org/2000/svg}text")]
    assert "Gene #1" in texts
