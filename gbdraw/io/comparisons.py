#!/usr/bin/env python
# coding: utf-8
"""Read BLAST outfmt 6/7 comparison tables through one reader and one normalizer."""

from __future__ import annotations

import codecs
import csv
import io
import logging
import os
from pathlib import Path
from typing import List

import numpy as np
import pandas as pd
from pandas import DataFrame
from pandas.api.types import is_bool_dtype, is_numeric_dtype

from ..configurators import BlastMatchConfigurator
from ..exceptions import ValidationError

logger = logging.getLogger(__name__)

COMPARISON_COLUMNS = (
    "query",
    "subject",
    "identity",
    "alignment_length",
    "mismatches",
    "gap_opens",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
)
_ID_COLUMNS = COMPARISON_COLUMNS[:2]
_NUMERIC_COLUMNS = COMPARISON_COLUMNS[2:]
_INTEGER_COLUMNS = frozenset(
    ("alignment_length", "mismatches", "gap_opens", "qstart", "qend", "sstart", "send")
)


class _InvalidCell(Exception):
    """One invalid value, located by data-row position before it is reported."""

    def __init__(self, row: int, column: str, detail: str) -> None:
        super().__init__(detail)
        self.row = row
        self.column = column
        self.detail = detail

    def location(self) -> str:
        return f"column {COMPARISON_COLUMNS.index(self.column) + 1} ({self.column})"


def _empty_comparison_frame() -> DataFrame:
    return DataFrame({column: pd.Series(dtype=object) for column in COMPARISON_COLUMNS})


def _normalize(frame: DataFrame) -> DataFrame:
    """Type the 12 positional columns of ``frame`` or raise ``_InvalidCell``."""

    frame = frame.reset_index(drop=True)
    first_invalid: tuple[int, int, str, str] | None = None

    def note(rows: np.ndarray, order: int, column: str, detail_for_row) -> None:
        nonlocal first_invalid
        if rows.size == 0:
            return
        row = int(rows[0])
        if first_invalid is None or (row, order) < first_invalid[:2]:
            first_invalid = (row, order, column, detail_for_row(row))

    normalized: dict[str, pd.Series] = {}
    for order, column in enumerate(COMPARISON_COLUMNS):
        series = frame[column]
        if column in _ID_COLUMNS:
            values = series.to_numpy(dtype=object)
            missing = pd.isna(values) | (values == "")
            note(np.flatnonzero(missing), order, column, lambda _row: "the sequence ID is empty")
            normalized[column] = series.astype(str)
            continue
        numeric = series if is_numeric_dtype(series) else pd.to_numeric(series, errors="coerce")
        values = numeric.to_numpy(dtype=float, na_value=np.nan)
        # Booleans are numeric to pandas but are never BLAST values.
        invalid = np.full(len(values), True) if is_bool_dtype(numeric) else ~np.isfinite(values)
        integer = column in _INTEGER_COLUMNS
        if integer:
            invalid |= np.isfinite(values) & (values != np.floor(values))
        kind = "is not an integer" if integer else "is not a finite number"
        note(
            np.flatnonzero(invalid),
            order,
            column,
            lambda row, series=series, kind=kind: f"{str(series.iloc[row])!r} {kind}",
        )
        if integer and numeric.dtype.kind != "i" and not invalid.any():
            numeric = numeric.astype("int64")
        normalized[column] = numeric
    if first_invalid is not None:
        row, _order, column, detail = first_invalid
        raise _InvalidCell(row, column, detail)
    # ``frame`` is already private to this call, so its columns need no copy.
    return DataFrame(normalized, columns=list(COMPARISON_COLUMNS), copy=False)


def normalize_comparison_dataframe(dataframe: DataFrame) -> DataFrame:
    """Return the 12 typed BLAST outfmt 6 columns of a comparison DataFrame.

    Named columns are selected by name; otherwise the first 12 columns are used
    by position. Invalid values raise ValidationError with the 1-based row.
    """

    if not isinstance(dataframe, DataFrame):
        raise ValidationError("A comparison table must be a pandas DataFrame.")
    if set(COMPARISON_COLUMNS).issubset(set(dataframe.columns)):
        frame = dataframe.loc[:, list(COMPARISON_COLUMNS)]
    elif len(dataframe.columns) >= len(COMPARISON_COLUMNS):
        frame = dataframe.iloc[:, : len(COMPARISON_COLUMNS)].copy()
        frame.columns = list(COMPARISON_COLUMNS)
    else:
        raise ValidationError(
            "A comparison DataFrame must contain the 12 BLAST outfmt 6 columns "
            f"({', '.join(COMPARISON_COLUMNS)})."
        )
    try:
        return _normalize(frame)
    except _InvalidCell as cell:
        raise ValidationError(
            f"Comparison DataFrame row {cell.row + 1}, {cell.location()}: {cell.detail}."
        ) from None


_LINE_STARTS_NEEDING_CHECK = frozenset((b"", b"#", b" ", b"\t"))


def _is_data_line(line: bytes) -> bool:
    # outfmt 7 comments start a line; a "#" inside a field is data.
    text = line.strip()
    return bool(text) and not text.startswith(b"#")


def _short_row_error(name: str, line_number: int, field_count: int) -> ValidationError:
    return ValidationError(
        f"{name}: line {line_number}: expected at least {len(COMPARISON_COLUMNS)} "
        f"tab-separated BLAST outfmt 6 columns; found {field_count}."
    )


def _scan_lines(lines: list[bytes], name: str) -> tuple[list[int], int]:
    """Return 0-based non-data line indexes and the widest data row (0 if none)."""

    skipped: list[int] = []
    widest = 0
    for index, line in enumerate(lines):
        width = line.count(b"\t") + 1
        # Fast path: an ordinary data row needs no whitespace or comment check.
        if width >= len(COMPARISON_COLUMNS) and line[:1] not in _LINE_STARTS_NEEDING_CHECK:
            widest = max(widest, width)
            continue
        if not _is_data_line(line):
            skipped.append(index)
            continue
        if width < len(COMPARISON_COLUMNS):
            raise _short_row_error(name, index + 1, width)
        widest = max(widest, width)
    return skipped, widest


def read_comparison_table(
    source: str | os.PathLike[str] | io.StringIO,
    *,
    label: str | None = None,
) -> DataFrame:
    """Read one BLAST outfmt 6/7 table into the 12 typed comparison columns.

    ``source`` is a UTF-8 file path or an ``io.StringIO`` buffer. Lines whose
    first non-blank character is ``#`` are outfmt 7 comments; blank and
    comment-only tables have no rows. Columns after the first 12 are ignored
    with an INFO log. A missing or unreadable file, a row with fewer than 12
    columns, or a value of the wrong type raises ValidationError with the line.
    """

    try:
        if isinstance(source, io.StringIO):
            name = label or "comparison table"
            data = source.getvalue().encode("utf-8")
        else:
            name = label or str(source)
            path = Path(source)
            if not path.is_file():
                raise ValidationError(f"{name}: the comparison file does not exist or is not a file.")
            data = path.read_bytes()
    except (OSError, UnicodeError) as exc:
        raise ValidationError(f"{name}: the comparison file could not be read: {exc}") from exc
    # bytes.splitlines() breaks at \n, \r\n and \r, as the pandas C parser does.
    lines = data.splitlines()
    if lines and lines[0].startswith(codecs.BOM_UTF8):
        lines[0] = lines[0][len(codecs.BOM_UTF8) :]
    skipped, widest = _scan_lines(lines, name)
    if not widest:
        return _empty_comparison_frame()
    if widest > len(COMPARISON_COLUMNS):
        logger.info(
            "INFO: %s: using the first 12 BLAST outfmt 6 columns and ignoring %d extra column(s).",
            name,
            widest - len(COMPARISON_COLUMNS),
        )
    try:
        raw = pd.read_csv(
            io.BytesIO(data),
            sep="\t",
            header=None,
            names=COMPARISON_COLUMNS,
            usecols=range(len(COMPARISON_COLUMNS)),
            dtype={column: str for column in _ID_COLUMNS},
            na_filter=False,
            quoting=csv.QUOTE_NONE,
            skiprows=skipped or None,
            encoding="utf-8-sig",
        )
    except UnicodeError as exc:
        raise ValidationError(f"{name}: the comparison file could not be read: {exc}") from exc
    except pd.errors.ParserError as exc:
        raise ValidationError(f"{name}: the comparison table could not be parsed: {exc}") from exc
    try:
        return _normalize(raw)
    except _InvalidCell as cell:
        data_line_numbers = (
            number for number, line in enumerate(lines, start=1) if _is_data_line(line)
        )
        line_number = next(
            (number for index, number in enumerate(data_line_numbers) if index == cell.row),
            cell.row + 1,
        )
        raise ValidationError(
            f"{name}: line {line_number}, {cell.location()}: {cell.detail}."
        ) from None


def filter_comparison_dataframe(
    df: DataFrame, blast_config: BlastMatchConfigurator
) -> DataFrame:
    """Apply pairwise match thresholds to a comparison DataFrame."""

    evalue_threshold: float = blast_config.evalue
    bitscore_threshold: float = blast_config.bitscore
    identity_threshold: float = blast_config.identity
    alignment_length_threshold: int = blast_config.alignment_length
    return df[
        (df["evalue"] <= evalue_threshold)
        & (df["bitscore"] >= bitscore_threshold)
        & (df["identity"] >= identity_threshold)
        & (df["alignment_length"] >= alignment_length_threshold)
    ]


def load_comparisons(
    comparison_files: List[str], blast_config: BlastMatchConfigurator
) -> List[DataFrame]:
    """Read and filter one comparison table per file, keeping their order.

    Every file must be readable: a skipped file would shift later tables onto
    the wrong record pair, so a missing or malformed file raises ValidationError.
    """

    logger.info(
        "INFO: BLAST output visualization settings: e-value threshold: {}; bitscore threshold: {}; identity threshold: {}; alignment length threshold: {}".format(
            blast_config.evalue,
            blast_config.bitscore,
            blast_config.identity,
            blast_config.alignment_length,
        )
    )
    comparison_list: list[DataFrame] = []
    logger.info("INFO: Loading comparison file(s)...")
    for comparison_file in comparison_files:
        logger.info("INFO: Loading {}".format(comparison_file))
        comparison_list.append(
            filter_comparison_dataframe(read_comparison_table(comparison_file), blast_config)
        )
    logger.info("INFO:             ... finished loading comparison file(s)")
    return comparison_list


__all__ = [
    "COMPARISON_COLUMNS",
    "filter_comparison_dataframe",
    "load_comparisons",
    "normalize_comparison_dataframe",
    "read_comparison_table",
]
