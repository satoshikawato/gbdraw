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
from typing import Any

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


def read_literal_table(
    source: Any,
    *,
    names: list[str],
    label: str,
    filepath: str | None = None,
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
    ``options`` are further ``pandas.read_csv`` options and may replace the defaults (the
    python engine, an error on a row with the wrong number of fields).
    """
    if isinstance(source, io.StringIO):
        text = source.getvalue()
    else:
        filepath = str(source)
        text = "".join(read_table_lines(filepath))
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
