#!/usr/bin/env python
# coding: utf-8

"""Shared reading of the styling and override TSV tables.

A whole-line comment is blanked, and a ``"`` is part of the cell value: the Web
writes each value as typed, so these tables are read with CSV quoting off.
"""

from __future__ import annotations

import csv
import io
from typing import Any

import pandas as pd


def read_table_lines(filepath: str) -> list[str]:
    """Return the file's lines; a whole-line comment becomes a blank line.

    A comment is a line whose first non-blank character is ``#``. A ``#`` after
    other text is part of the cell value (``foo#bar``, ``Gene #1``), so the
    readers must not use the pandas ``comment`` option, which cuts a line at
    its first ``#``. A blanked line keeps its position, so line numbers in later
    error messages still match the file. A leading UTF-8 byte order mark is
    dropped, so it cannot hide a first-line comment.
    """
    with open(filepath, "r", encoding="utf-8-sig") as handle:
        return [
            "\n" if raw_line.lstrip().startswith("#") else raw_line
            for raw_line in handle
        ]


def table_text_stream(lines: list[str]) -> io.StringIO:
    """Return the lines from :func:`read_table_lines` as one stream for ``pandas.read_csv``."""
    return io.StringIO("".join(lines))


def read_literal_table(source: Any, *, names: list[str], **options: Any) -> pd.DataFrame:
    """Read a tab-separated styling table; every cell is text and a ``"`` is a plain character.

    ``source`` is a path or the stream from :func:`table_text_stream`. ``options`` are
    further ``pandas.read_csv`` options and may replace the defaults (the python
    engine, an error on a row with the wrong number of fields).
    """
    options.setdefault("engine", "python")
    options.setdefault("on_bad_lines", "error")
    return pd.read_csv(
        source,
        sep="\t",
        header=None,
        names=names,
        dtype=str,
        quoting=csv.QUOTE_NONE,
        **options,
    )
