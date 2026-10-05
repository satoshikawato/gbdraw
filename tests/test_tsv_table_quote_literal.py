"""A ``"`` is part of the value in the styling and override TSV tables (OV-22).

The Web writes these tables with each value as typed, so a rule value that
starts with a quote (``"lead``) must read back unchanged and must not break
Generate. The readers therefore turn CSV quoting off.
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
from gbdraw.features.visibility import read_feature_visibility_file
from gbdraw.io.colors import load_default_colors, read_color_table
from gbdraw.labels.filtering import (
    read_filter_list_file,
    read_label_override_file,
    read_qualifier_priority_file,
)

# (reader, a row whose last cell carries the quote variant, which cell index)
_READERS: list[Any] = [
    pytest.param(
        read_feature_visibility_file,
        ["*", "CDS", "product", "{v}", "off"],
        3,
        id="feature-visibility",
    ),
    pytest.param(
        read_label_override_file,
        ["*", "CDS", "product", "kinase", "{v}"],
        4,
        id="label-override",
    ),
    pytest.param(
        read_filter_list_file,
        ["CDS", "product", "{v}"],
        2,
        id="label-filter-list",
    ),
    pytest.param(
        read_qualifier_priority_file,
        ["CDS", "{v}"],
        1,
        id="qualifier-priority",
    ),
    pytest.param(
        read_color_table,
        ["CDS", "product", "{v}", "#ff0000"],
        2,
        id="specific-colors",
    ),
]

_VALUES = ['"quoted"', '"lead', 'trail"', 'say "hi" now', '"a b"']


def _read_values(reader: Callable[[str], Any], path: Path, template: list[str], values: list[str]) -> list[list[str]]:
    path.write_text(
        "".join("\t".join(cell.replace("{v}", value) for cell in template) + "\n" for value in values),
        encoding="utf-8",
    )
    frame = reader(str(path))
    return [[str(cell) for cell in row[: len(template)]] for row in frame.values.tolist()]


@pytest.mark.parametrize(("reader", "template", "index"), _READERS)
def test_a_quote_is_part_of_the_cell_value(
    tmp_path: Path, reader: Callable[[str], Any], template: list[str], index: int
) -> None:
    rows = _read_values(reader, tmp_path / "table.tsv", template, _VALUES)

    assert [row[index] for row in rows] == _VALUES


@pytest.mark.parametrize(("reader", "template", "index"), _READERS)
def test_an_unclosed_leading_quote_does_not_swallow_later_rows(
    tmp_path: Path, reader: Callable[[str], Any], template: list[str], index: int
) -> None:
    values = ['"lead', "after one", "after two"]

    rows = _read_values(reader, tmp_path / "table.tsv", template, values)

    assert [row[index] for row in rows] == values


def test_default_color_table_keeps_a_quote_in_the_feature_type(tmp_path: Path) -> None:
    path = tmp_path / "default_colors.tsv"
    path.write_text('"lead\t#ff0000\n"quoted"\t#00ff00\n', encoding="utf-8")

    frame = load_default_colors(str(path)).set_index("feature_type")

    assert frame.loc['"lead', "color"] == "#ff0000"
    assert frame.loc['"quoted"', "color"] == "#00ff00"


def _write_genbank(path: Path) -> None:
    record = SeqRecord(
        Seq("ACGT" * 150),
        id="rec1",
        name="rec1",
        description="quote fixture",
        annotations={"molecule_type": "DNA", "topology": "circular"},
    )
    record.features = [
        SeqFeature(FeatureLocation(10, 190, strand=1), type="CDS",
                   qualifiers={"product": ['"lead kinase']}),
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


def _label_texts(svg: str) -> list[str]:
    root = ET.fromstring(svg)
    return ["".join(node.itertext()) for node in root.iter("{http://www.w3.org/2000/svg}text")]


def test_generate_accepts_a_rule_value_with_an_unclosed_quote(tmp_path: Path) -> None:
    visibility = tmp_path / "visibility.tsv"
    visibility.write_text('*\tCDS\tproduct\t"lead kinase\toff\n', encoding="utf-8")
    labels = tmp_path / "labels.tsv"
    labels.write_text('*\tCDS\tproduct\tother protein\t"Gene\n', encoding="utf-8")
    whitelist = tmp_path / "whitelist.tsv"
    whitelist.write_text('CDS\tproduct\t"lead\nCDS\tproduct\tother\n', encoding="utf-8")

    svg = _run_circular(
        tmp_path,
        "--labels",
        "--feature_visibility_table", str(visibility),
        "--label_table", str(labels),
        "--label_whitelist", str(whitelist),
    )

    assert len(set(re.findall(r'data-gbdraw-feature-id="([^"]+)"', svg))) == 1
    assert '"Gene' in _label_texts(svg)


def test_generate_keeps_the_quotes_of_a_label_text(tmp_path: Path) -> None:
    labels = tmp_path / "labels.tsv"
    labels.write_text('*\tCDS\tproduct\tother protein\t"Gene"\n', encoding="utf-8")

    svg = _run_circular(tmp_path, "--labels", "--label_table", str(labels))

    assert '"Gene"' in _label_texts(svg)
