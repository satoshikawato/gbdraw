#!/usr/bin/env python
# coding: utf-8

"""Shared reading of the styling, override, and annotation TSV tables.

A whole-line comment is blanked, and a ``"`` is part of the cell value: the Web
writes each value as typed, so these tables are read with CSV quoting off.
"""

from __future__ import annotations

import csv
import io
import logging
import re
from contextlib import contextmanager
from contextvars import ContextVar
from typing import Any, Callable, Iterator, NamedTuple, Sequence

import pandas as pd

from ..exceptions import ParseError

logger = logging.getLogger(__name__)


def read_table_lines(filepath: str) -> list[str]:
    """Return the file's lines; a whole-line comment or a whitespace-only line becomes a blank line.

    A comment is a line whose first non-blank character is ``#``. A ``#`` after
    other text is part of the cell value (``foo#bar``, ``Gene #1``), so the
    readers must not use the pandas ``comment`` option, which cuts a line at
    its first ``#``. A line of only tabs and spaces has no cell, as in the Web
    readers (``services/file-imports.js``). A blanked line keeps its position,
    so line numbers in later error messages still match the file. A leading
    UTF-8 byte order mark is dropped, so it cannot hide a first-line comment.
    """
    with open(filepath, "r", encoding="utf-8-sig") as handle:
        return [
            "\n" if not raw_line.strip() or raw_line.lstrip().startswith("#") else raw_line
            for raw_line in handle
        ]


def table_text_stream(lines: list[str]) -> io.StringIO:
    """Return the lines from :func:`read_table_lines` as one stream for ``pandas.read_csv``."""
    return io.StringIO("".join(lines))


class LegacyRows(NamedTuple):
    """How a Session 31-39 table reads a row: the cells it requires, whether ``[`` lines are
    skipped, and which complete rows hold valid values (``None``: every one)."""

    complete: Callable[[Sequence[str]], bool]
    skips_sections: bool = False
    valid: Callable[[Sequence[str]], bool] | None = None


# Sessions 31-39 stored the Default colors, Label whitelist and Qualifier
# priority tables as their Web writer wrote them, before cell values were
# normalized: a tab in a value made extra cells and a line break a short row.
# While such a Session is read (`legacy_table_rows`), these tables read a row as
# the current writer writes it, as Web Load does (services/file-imports.js
# `tableCells`, OV-40): the extra cells join the last cell with one space, and a
# row without its required cells is dropped. A Default colors row whose color is
# outside the Default colors forms is dropped too and listed as invalid, as Web
# Load does (`parseColorTable`, OV-272). The resource bytes stay as saved.
def _keyed_row(cells: Sequence[str]) -> bool:
    return bool(cells[0] and cells[1])


def _default_colors_row_valid(cells: Sequence[str]) -> bool:
    """The ``feature_type<TAB>color`` header, or a row with a documented user color."""

    from .colors import is_user_color  # colors.py reads its tables through this module

    key, color = cells[0], cells[1]
    return (key.lower() == "feature_type" and color.lower() == "color") or is_user_color(color)


LEGACY_TABLE_ROWS = {
    "label-whitelist": LegacyRows(_keyed_row),
    "qualifier-priority": LegacyRows(_keyed_row),
    "default-colors": LegacyRows(_keyed_row, skips_sections=True, valid=_default_colors_row_valid),
}
_READS_LEGACY_TABLE_ROWS: ContextVar[bool] = ContextVar("gbdraw_legacy_table_rows", default=False)


@contextmanager
def legacy_table_rows(active: bool) -> Iterator[None]:
    """Read the tables of :data:`LEGACY_TABLE_ROWS` as Session 31-39 tables while ``active``."""

    token = _READS_LEGACY_TABLE_ROWS.set(active or _READS_LEGACY_TABLE_ROWS.get())
    try:
        yield
    finally:
        _READS_LEGACY_TABLE_ROWS.reset(token)


def repair_legacy_table_text(
    text: str, rows: LegacyRows, *, columns: int
) -> tuple[str, list[dict[str, Any]]]:
    """Return the table as the current writer writes it and its repaired rows by 1-based line.

    A blank, comment or skipped section line stays; a dropped row becomes a blank
    line, so later line numbers still match the file.
    """

    lines: list[str] = []
    repairs: list[dict[str, Any]] = []
    for row, line in enumerate(text.split("\n"), start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#") or (rows.skips_sections and stripped.startswith("[")):
            lines.append(line)
            continue
        parts = line.rstrip("\r").split("\t")
        cells = [
            re.sub(r"[\t\r\n]+", " ", cell).strip()
            for cell in (*parts[: columns - 1], "\t".join(parts[columns - 1 :]))
        ]
        if len(parts) < columns or not rows.complete(cells):
            repairs.append({"row": row, "repair": "dropped"})
            lines.append("")
            continue
        if len(parts) > columns:
            repairs.append({"row": row, "repair": "joined"})
        if rows.valid is not None and not rows.valid(cells):
            repairs.append({"row": row, "repair": "invalid"})
            lines.append("")
            continue
        lines.append("\t".join(cells))
    return "\n".join(lines), repairs


def read_literal_table(
    source: Any,
    *,
    names: list[str],
    label: str,
    filepath: str | None = None,
    legacy_rows: LegacyRows | None = None,
    **options: Any,
) -> pd.DataFrame:
    """Read a tab-separated styling table; every cell is text and a ``"`` is a plain character.

    ``source`` is a path, whose comments and whitespace-only lines are blanked as in
    :func:`read_table_lines`, or the stream from :func:`table_text_stream`. A blank line
    has no row, and a row's index is its zero-based line in the file. A row with more
    cells than ``names`` raises :class:`ParseError` with the file and line: pandas would
    otherwise take the extra leading cells as an index and read the row shifted. A row
    with fewer cells is left to the caller's missing-value check. ``label`` names the
    table in the message and ``filepath`` names the file when ``source`` is a stream.
    Within :func:`legacy_table_rows`, a table with ``legacy_rows`` is first read as
    :func:`repair_legacy_table_text` returns it, and its repaired rows are logged.
    ``options`` are further ``pandas.read_csv`` options and may replace the defaults (the
    python engine, an error on a row with the wrong number of fields).
    """
    if isinstance(source, io.StringIO):
        text = source.getvalue()
    else:
        filepath = str(source)
        text = "".join(read_table_lines(filepath))
    if legacy_rows is not None and _READS_LEGACY_TABLE_ROWS.get():
        text, repairs = repair_legacy_table_text(text, legacy_rows, columns=len(names))
        for kind, reason in (
            ("joined", "had extra cells, joined into the last column with one space"),
            ("dropped", "lacked a required column and were dropped"),
            ("invalid", "had a value outside the documented forms and were dropped"),
        ):
            lines = [str(repair["row"]) for repair in repairs if repair["repair"] == kind]
            if lines:
                logger.warning(
                    f"WARNING: The {label} '{filepath}' of a Session 31-39 was read as the "
                    f"current version writes it: line(s) {', '.join(lines)} {reason}."
                )
    rows = [
        (line_no, raw_line)
        for line_no, raw_line in enumerate(text.split("\n"), start=1)
        if raw_line.strip()
    ]
    for line_no, raw_line in rows:
        if raw_line.rstrip("\r").count("\t") + 1 > len(names):
            message = (
                f"Malformed line in {label} '{filepath}' at line {line_no}: "
                f"expected {len(names)} columns."
            )
            logger.error(f"ERROR: {message}")
            raise ParseError(message)
    options.setdefault("engine", "python")
    options.setdefault("on_bad_lines", "error")
    frame = pd.read_csv(
        io.StringIO("\n".join(raw_line for _line_no, raw_line in rows)),
        sep="\t",
        header=None,
        names=names,
        dtype=str,
        quoting=csv.QUOTE_NONE,
        **options,
    )
    # One row per non-blank line; its index is the zero-based file line, so
    # ``index + 1`` in a later message is the line in the file.
    frame.index = pd.Index([line_no - 1 for line_no, _raw_line in rows])
    return frame
