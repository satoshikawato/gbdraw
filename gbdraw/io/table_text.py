#!/usr/bin/env python
# coding: utf-8

"""Text of the label and visibility TSV tables, with whole-line comments blanked."""

from __future__ import annotations

import io


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
