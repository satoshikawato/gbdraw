"""A row with more columns than the table has is an error, never a shifted row (OV-25).

pandas takes the extra leading cells of a long row as the index, so
``CDS  product  two  words`` was read as ``product / two / words``. Every reader
of the styling and override tables rejects such a row with the file and line.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Callable

import pytest

from gbdraw.exceptions import ParseError
from gbdraw.features.visibility import read_feature_visibility_file
from gbdraw.io.colors import load_default_colors, read_color_table
from gbdraw.labels.filtering import (
    read_filter_list_file,
    read_label_override_file,
    read_qualifier_priority_file,
)

# (reader, a valid row, the same row plus one more cell)
_READERS: list[Any] = [
    pytest.param(
        read_feature_visibility_file,
        ["*", "CDS", "product", "kinase", "off"],
        id="feature-visibility",
    ),
    pytest.param(
        read_label_override_file,
        ["*", "CDS", "product", "kinase", "Gene"],
        id="label-override",
    ),
    pytest.param(
        read_filter_list_file,
        ["CDS", "product", "kinase"],
        id="label-filter-list",
    ),
    pytest.param(
        read_qualifier_priority_file,
        ["CDS", "gene,product"],
        id="qualifier-priority",
    ),
    pytest.param(
        read_color_table,
        ["CDS", "product", "kinase", "#ff0000", "caption"],
        id="specific-colors",
    ),
    pytest.param(
        lambda path: load_default_colors(path).set_index("feature_type").reset_index(),
        ["CDS", "#ff0000"],
        id="default-colors",
    ),
]


def _write(path: Path, rows: list[list[str]], *, before: tuple[str, ...] = ()) -> str:
    path.write_text(
        "".join(line + "\n" for line in before) + "".join("\t".join(row) + "\n" for row in rows),
        encoding="utf-8",
    )
    return str(path)


@pytest.mark.parametrize("position", ["first", "later"])
@pytest.mark.parametrize(("reader", "row"), _READERS)
def test_a_row_with_an_extra_column_is_rejected_with_its_line(
    tmp_path: Path, reader: Callable[[str], Any], row: list[str], position: str
) -> None:
    long_row = [*row, "extra"]
    rows = [long_row, row, row] if position == "first" else [row, long_row, row]
    path = _write(tmp_path / "table.tsv", rows, before=("",))
    line = 2 if position == "first" else 3

    with pytest.raises(ParseError, match=rf"line {line}\b") as raised:
        reader(path)

    assert "table.tsv" in str(raised.value)


@pytest.mark.parametrize(("reader", "row"), _READERS)
def test_a_row_that_splits_a_value_on_a_tab_is_not_read_shifted(
    tmp_path: Path, reader: Callable[[str], Any], row: list[str]
) -> None:
    split = [*row[:-1], "two", "words"]

    with pytest.raises(ParseError, match=r"line 1\b"):
        reader(_write(tmp_path / "table.tsv", [split]))


@pytest.mark.parametrize(("reader", "row"), _READERS)
def test_valid_rows_are_unchanged(
    tmp_path: Path, reader: Callable[[str], Any], row: list[str]
) -> None:
    path = _write(tmp_path / "table.tsv", [row, row], before=("",))

    frame = reader(path)

    assert row in [list(map(str, values))[: len(row)] for values in frame.values.tolist()]
