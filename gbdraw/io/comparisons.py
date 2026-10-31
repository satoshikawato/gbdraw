#!/usr/bin/env python
# coding: utf-8
"""Read BLAST outfmt 6/7 comparison tables through one reader and one normalizer."""

from __future__ import annotations

import codecs
import csv
import io
import logging
import os
import sys
from pathlib import Path
from dataclasses import dataclass
from typing import Literal, Protocol

import numpy as np
import pandas as pd
from pandas import DataFrame
from pandas.api.types import is_bool_dtype, is_numeric_dtype

from ..exceptions import ValidationError
from .table_text import EXACT_FLOAT_READ_OPTIONS

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
# Producer diagnostics (gbdraw.web_support.error_adapter): identifiers and
# 1-based line/column locators only, never a file name or cell value.
_UNREADABLE = {"code": "INPUT_UNREADABLE", "field": "comparison"}
# The outfmt 7 "# Fields:" names of the 12 standard columns; BLAST+ names the
# IDs after the -outfmt specifier (qseqid "query id", qaccver "query acc.ver").
_STANDARD_FIELD_NAMES = (
    "query acc.ver", "subject acc.ver", "% identity", "alignment length", "mismatches",
    "gap opens", "q. start", "q. end", "s. start", "s. end", "evalue", "bit score",
)
_FIELD_NAMES = {
    **dict(zip(_STANDARD_FIELD_NAMES, COMPARISON_COLUMNS)),
    **dict.fromkeys(("query id", "query acc.", "query gi"), "query"),
    **dict.fromkeys(("subject id", "subject acc.", "subject gi"), "subject"),
}
_STANDARD_POSITIONS = tuple(range(len(COMPARISON_COLUMNS)))


class _InvalidCell(Exception):
    """One invalid value, located by data-row position before it is reported."""

    def __init__(self, row: int, column: str, reason: str, detail: str) -> None:
        super().__init__(detail)
        self.row = row
        self.column = column
        self.reason = reason
        self.detail = detail
        # The 1-based column in the file; read_comparison_table sets it for a
        # table read in the order of its "# Fields:" line.
        self.column_number = COMPARISON_COLUMNS.index(column) + 1

    def location(self) -> str:
        return f"column {self.column_number} ({self.column})"

    def diagnostic(self, row: int) -> dict[str, object]:
        return {
            "code": "COMPARISON_INPUT",
            "reason": self.reason,
            "row": row,
            "column": self.column_number,
        }


def _empty_comparison_frame() -> DataFrame:
    return DataFrame({column: pd.Series(dtype=object) for column in COMPARISON_COLUMNS})


def _normalize(frame: DataFrame) -> DataFrame:
    """Type the 12 positional columns of ``frame`` or raise ``_InvalidCell``."""

    frame = frame.reset_index(drop=True)
    first_invalid: tuple[int, int, str, str, str] | None = None

    def note(rows: np.ndarray, order: int, column: str, reason: str, detail_for_row) -> None:
        nonlocal first_invalid
        if rows.size == 0:
            return
        row = int(rows[0])
        if first_invalid is None or (row, order) < first_invalid[:2]:
            first_invalid = (row, order, column, reason, detail_for_row(row))

    normalized: dict[str, pd.Series] = {}
    for order, column in enumerate(COMPARISON_COLUMNS):
        series = frame[column]
        if column in _ID_COLUMNS:
            values = series.to_numpy(dtype=object)
            missing = pd.isna(values) | (values == "")
            note(
                np.flatnonzero(missing),
                order,
                column,
                "REQUIRED",
                lambda _row: "the sequence ID is empty",
            )
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
            "INTEGER" if integer else "FINITE",
            lambda row, series=series, kind=kind: f"{str(series.iloc[row])!r} {kind}",
        )
        if integer and numeric.dtype.kind != "i" and not invalid.any():
            numeric = numeric.astype("int64")
        normalized[column] = numeric
    if first_invalid is not None:
        row, _order, column, reason, detail = first_invalid
        raise _InvalidCell(row, column, reason, detail)
    # ``frame`` is already private to this call, so its columns need no copy.
    return DataFrame(normalized, columns=list(COMPARISON_COLUMNS), copy=False)


def normalize_comparison_dataframe(dataframe: DataFrame) -> DataFrame:
    """Return the 12 typed BLAST outfmt 6 columns of a comparison DataFrame.

    Named columns are selected by name; otherwise the first 12 columns are used
    by position. Invalid values raise ValidationError with the 1-based row.
    """

    if not isinstance(dataframe, DataFrame):
        raise ValidationError(
            "A comparison table must be a pandas DataFrame.",
            diagnostic={"code": "COMPARISON_INPUT", "field": "comparison"},
        )
    if set(COMPARISON_COLUMNS).issubset(set(dataframe.columns)):
        frame = dataframe.loc[:, list(COMPARISON_COLUMNS)]
    elif len(dataframe.columns) >= len(COMPARISON_COLUMNS):
        frame = dataframe.iloc[:, : len(COMPARISON_COLUMNS)].copy()
        frame.columns = list(COMPARISON_COLUMNS)
    else:
        raise ValidationError(
            "A comparison DataFrame must contain the 12 BLAST outfmt 6 columns "
            f"({', '.join(COMPARISON_COLUMNS)}).",
            diagnostic={
                "code": "COMPARISON_INPUT",
                "reason": "FIELDS",
                "columnCount": len(COMPARISON_COLUMNS),
            },
        )
    try:
        return _normalize(frame)
    except _InvalidCell as cell:
        raise ValidationError(
            f"Comparison DataFrame row {cell.row + 1}, {cell.location()}: {cell.detail}.",
            diagnostic=cell.diagnostic(cell.row + 1),
        ) from None


_LINE_STARTS_NEEDING_CHECK = frozenset((b"", b"#", b" ", b"\t"))


def _is_data_line(line: bytes) -> bool:
    # outfmt 7 comments start a line; a "#" inside a field is data.
    text = line.strip()
    return bool(text) and not text.startswith(b"#")


def _short_row_error(name: str, line_number: int, field_count: int, required: int) -> ValidationError:
    return ValidationError(
        f"{name}: line {line_number}: expected at least {required} "
        f"tab-separated BLAST outfmt 6 columns; found {field_count}.",
        diagnostic={
            "code": "COMPARISON_INPUT",
            "reason": "FIELDS",
            "row": line_number,
            "columnCount": required,
        },
    )


def _field_positions(line: bytes, line_number: int, name: str) -> tuple[int, ...]:
    """The 0-based file position of each standard column named by a "# Fields:" line.

    A line that names every standard field gives their order (CI-07d). A line
    whose recognized fields all sit at their standard positions keeps the
    positional read; any other line would misread the table, so it fails.
    """

    names = line.strip()[len(b"# Fields:"):].decode("utf-8", "replace").split(",")
    positions: dict[str, int] = {}
    for position, field_name in enumerate(names):
        column = _FIELD_NAMES.get(field_name.strip().lower())
        if column is not None:
            positions.setdefault(column, position)
    if len(positions) == len(COMPARISON_COLUMNS):
        return tuple(positions[column] for column in COMPARISON_COLUMNS)
    if all(position == COMPARISON_COLUMNS.index(column) for column, position in positions.items()):
        return _STANDARD_POSITIONS
    missing = [label for column, label in zip(COMPARISON_COLUMNS, _STANDARD_FIELD_NAMES) if column not in positions]
    raise ValidationError(
        f"{name}: line {line_number}: the outfmt 7 '# Fields:' line lists the standard BLAST "
        f"outfmt 6 fields in another order but lacks {', '.join(missing)}.",
        diagnostic={"code": "COMPARISON_INPUT", "reason": "OUTFMT7_FIELDS", "row": line_number},
    )


def _scan_lines(lines: list[bytes], name: str) -> tuple[list[int], int, tuple[int, ...]]:
    """Return non-data line indexes, the widest data row (0 if none), and column positions.

    Indexes are 0-based; the positions come from the first outfmt 7 "# Fields:" line.
    """

    skipped: list[int] = []
    widest = 0
    positions: tuple[int, ...] | None = None
    narrow_width, narrow_index = sys.maxsize, 0  # the narrowest data row
    for index, line in enumerate(lines):
        width = line.count(b"\t") + 1
        # Fast path: an ordinary data row needs no whitespace or comment check.
        if width >= len(COMPARISON_COLUMNS) and line[:1] not in _LINE_STARTS_NEEDING_CHECK:
            widest = max(widest, width)
            if width < narrow_width:
                narrow_width, narrow_index = width, index
            continue
        if not _is_data_line(line):
            skipped.append(index)
            if positions is None and line.strip().startswith(b"# Fields:"):
                positions = _field_positions(line, index + 1, name)
            continue
        if width < len(COMPARISON_COLUMNS):
            raise _short_row_error(name, index + 1, width, len(COMPARISON_COLUMNS))
        widest = max(widest, width)
        if width < narrow_width:
            narrow_width, narrow_index = width, index
    positions = positions or _STANDARD_POSITIONS
    if widest and narrow_width <= max(positions):
        raise _short_row_error(name, narrow_index + 1, narrow_width, max(positions) + 1)
    return skipped, widest, positions


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
    columns, or a value of the wrong type raises ValidationError with the line;
    its ``diagnostic`` carries the line as ``row`` and, for a value, the column.
    """

    try:
        if isinstance(source, io.StringIO):
            name = label or "comparison table"
            data = source.getvalue().encode("utf-8")
        else:
            name = label or str(source)
            path = Path(source)
            if not path.is_file():
                raise ValidationError(
                    f"{name}: the comparison file does not exist or is not a file.",
                    diagnostic=_UNREADABLE,
                )
            data = path.read_bytes()
    except (OSError, UnicodeError) as exc:
        raise ValidationError(
            f"{name}: the comparison file could not be read: {exc}", diagnostic=_UNREADABLE
        ) from exc
    # bytes.splitlines() breaks at \n, \r\n and \r, as the pandas C parser does.
    lines = data.splitlines()
    if lines and lines[0].startswith(codecs.BOM_UTF8):
        lines[0] = lines[0][len(codecs.BOM_UTF8) :]
    skipped, widest, positions = _scan_lines(lines, name)
    if not widest:
        return _empty_comparison_frame()
    if widest > len(COMPARISON_COLUMNS):
        logger.info(
            "INFO: %s: using the 12 standard BLAST outfmt 6 columns and ignoring %d extra column(s).",
            name,
            widest - len(COMPARISON_COLUMNS),
        )
    file_columns = [f"column_{position}" for position in range(max(positions) + 1)]
    try:
        raw = pd.read_csv(
            io.BytesIO(data),
            sep="\t",
            header=None,
            names=file_columns,
            usecols=[file_columns[position] for position in positions],
            dtype={file_columns[positions[index]]: str for index in range(len(_ID_COLUMNS))},
            na_filter=False,
            quoting=csv.QUOTE_NONE,
            skiprows=skipped or None,
            encoding="utf-8-sig",
            **EXACT_FLOAT_READ_OPTIONS,
        )
        raw = raw.loc[:, [file_columns[position] for position in positions]]
        raw.columns = list(COMPARISON_COLUMNS)
    except UnicodeError as exc:
        raise ValidationError(
            f"{name}: the comparison file could not be read: {exc}", diagnostic=_UNREADABLE
        ) from exc
    except pd.errors.ParserError as exc:
        raise ValidationError(
            f"{name}: the comparison table could not be parsed: {exc}",
            diagnostic={"code": "COMPARISON_INPUT"},
        ) from exc
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
        cell.column_number = positions[COMPARISON_COLUMNS.index(cell.column)] + 1
        raise ValidationError(
            f"{name}: line {line_number}, {cell.location()}: {cell.detail}.",
            diagnostic=cell.diagnostic(line_number),
        ) from None


class ComparisonFilter(Protocol):
    """The thresholds a comparison table is filtered by."""

    @property
    def evalue(self) -> float: ...

    @property
    def bitscore(self) -> float: ...

    @property
    def identity(self) -> float: ...

    @property
    def alignment_length(self) -> int: ...


def read_comparison_resource(source: str | os.PathLike[str]) -> DataFrame:
    """``read_comparison_table`` of one comparison file, parsed once per Worker resource.

    In a Web render the staged resource's identity keys the frame in the
    prepared-input cache, so a table whose resource is unchanged is not parsed
    again; elsewhere the file is read. The frame is shared between renders:
    a reader copies it before changing it.
    """

    # Deferred: gbdraw.api imports this module.
    from gbdraw.api.prepared import get_or_build_decoded_resource

    return get_or_build_decoded_resource(
        Path(source), ("comparison-table-v1",), lambda: read_comparison_table(source)
    )


ComparisonThresholdName = Literal["evalue", "bitscore", "identity", "alignment_length"]
# A row passes a threshold at or below it (``evalue``) or at or above it (the others).
_THRESHOLD_UPPER_BOUNDS: dict[ComparisonThresholdName, bool] = {
    "evalue": True, "bitscore": False, "identity": False, "alignment_length": False,
}


def _threshold_passes(df: DataFrame, blast_config: ComparisonFilter) -> dict[ComparisonThresholdName, pd.Series]:
    """Each threshold's test of every row: the one rule the filter and its tally read."""

    return {
        name: (df[name] <= getattr(blast_config, name)) if upper else (df[name] >= getattr(blast_config, name))
        for name, upper in _THRESHOLD_UPPER_BOUNDS.items()
    }


def filter_comparison_dataframe(
    df: DataFrame, blast_config: ComparisonFilter
) -> DataFrame:
    """Apply pairwise match thresholds to a comparison DataFrame."""

    passes = _threshold_passes(df, blast_config)
    return df[passes["evalue"] & passes["bitscore"] & passes["identity"] & passes["alignment_length"]]


@dataclass(frozen=True)
class ComparisonThresholdTally:
    """The rows of a comparison table one threshold leaves out.

    ``best_value`` is the value nearest the threshold among those rows: the lowest
    e-value, else the highest value.
    """

    name: ComparisonThresholdName
    threshold: float
    rows_left_out: int
    best_value: float


@dataclass(frozen=True)
class ComparisonRowTally:
    """How many rows of a comparison table were read and kept, and why the others were left out."""

    rows_read: int
    rows_kept: int
    # Only the thresholds that leave a row out; a row below two thresholds counts in both.
    thresholds: tuple[ComparisonThresholdTally, ...]


def comparison_row_tally(df: DataFrame, blast_config: ComparisonFilter) -> ComparisonRowTally:
    """The tally of :func:`filter_comparison_dataframe` on ``df``, by the same threshold tests."""

    passes = _threshold_passes(df, blast_config)
    kept = passes["evalue"] & passes["bitscore"] & passes["identity"] & passes["alignment_length"]
    thresholds = []
    for name, upper in _THRESHOLD_UPPER_BOUNDS.items():
        values = df.loc[~passes[name], name]
        if values.empty:
            continue
        thresholds.append(ComparisonThresholdTally(
            name=name,
            threshold=float(getattr(blast_config, name)),
            rows_left_out=int(values.size),
            best_value=float(values.min() if upper else values.max()),
        ))
    return ComparisonRowTally(rows_read=int(len(df)), rows_kept=int(kept.sum()), thresholds=tuple(thresholds))


__all__ = [
    "COMPARISON_COLUMNS",
    "ComparisonFilter",
    "ComparisonRowTally",
    "ComparisonThresholdName",
    "ComparisonThresholdTally",
    "comparison_row_tally",
    "filter_comparison_dataframe",
    "normalize_comparison_dataframe",
    "read_comparison_resource",
    "read_comparison_table",
]
