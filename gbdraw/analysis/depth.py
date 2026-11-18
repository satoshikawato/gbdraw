#!/usr/bin/env python
# coding: utf-8

"""Depth coverage parsing and binning helpers.

A depth file is a table of numbers along a reference. Its format (``DepthFileFormat``) is stated, or ``auto``:
then it is told apart by the column count of its first data line (never by the file's name, which the Web does not
keep): two columns are a 1-based position and a value (``position_value``), three are a reference name, a 1-based
position and a value (``samtools``, ``samtools depth``), and four or more are a reference name, a 0-based start, an
exclusive end and a value (``bedgraph``). A ``bed`` file is never told by its columns: its rows are features, each
adding its score (column 5, else 1) to every base it covers. Values may be negative. Every format is read into one
parsed form, :class:`DepthCoverage`.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Callable, Generator, Iterable, Iterator, Literal, Mapping, TypeAlias, get_args

import numpy as np
import pandas as pd
from Bio.SeqRecord import SeqRecord
from numpy.typing import NDArray
from pandas import DataFrame

from gbdraw.analysis.windows import SlidingWindows, SourcePiece, sliding_windows
from gbdraw.exceptions import ParseError, ValidationError
from gbdraw.io.input_text import GZIP_ERRORS, open_input_text
from gbdraw.io.table_text import EXACT_FLOAT_READ_OPTIONS

DEPTH_COLUMNS = ["reference_name", "position", "depth"]
# How a depth file reads: by its column count (``auto``), or as the format stated.
DepthFileFormat = Literal["auto", "samtools", "bedgraph", "bed", "position_value"]
DEPTH_FILE_FORMATS: tuple[DepthFileFormat, ...] = get_args(DepthFileFormat)
# The value of a window: the mean of the bases it counts, their largest value, or their sum.
DepthAggregation = Literal["mean", "max", "sum"]
DEPTH_AGGREGATIONS: tuple[DepthAggregation, ...] = get_args(DepthAggregation)
# The columns of each form, by the column count of the first data line (four or more: bedGraph).
_POSITION_VALUE_COLUMNS = ("position", "depth")
_INTERVAL_COLUMNS = ("reference_name", "start", "end", "depth")
_FORM_COLUMNS = {2: _POSITION_VALUE_COLUMNS, 3: tuple(DEPTH_COLUMNS), 4: _INTERVAL_COLUMNS}
# A BED row: a feature and its optional name and score. The name is not read.
_BED_COLUMNS = ("reference_name", "start", "end", "name", "depth")
# The columns a stated format reads, and how many a data line needs at least.
_FORMAT_COLUMNS: Mapping[DepthFileFormat, tuple[tuple[str, ...], int]] = {
    "samtools": (tuple(DEPTH_COLUMNS), 3),
    "bedgraph": (_INTERVAL_COLUMNS, 4),
    "bed": (_BED_COLUMNS, 3),
    "position_value": (_POSITION_VALUE_COLUMNS, 2),
}
_TEXT_COLUMNS = frozenset({"reference_name", "name"})
_COLUMN_DTYPES = {
    # A few names repeat on every row: codes, not one string per row.
    "reference_name": "category",
    "position": "int64",
    "start": "int64",
    "end": "int64",
    "depth": "float64",
}
# A line before the rows that is not one: a comment, or a UCSC track or browser line.
_PREAMBLE_LINE = re.compile(r"#|(?:track|browser)(?: |\r?$)")
# The header line ``samtools depth -H`` writes: the rows are samtools depth rows, one value column per sample.
_SAMTOOLS_HEADER = "#CHROM\tPOS\t"
# The reference of a file without reference names (two columns): every record that reads the file takes it.
_UNNAMED_REFERENCE = ""
_DEPTH_TSV_CACHE_MAXSIZE = 8
# Rows per chunk read: a chunk's transient columns stay small beside the stored depths.
_CHUNK_ROWS = 1 << 20
# Bases summed per chunk when binning: the running sums stay small beside the stored depths.
_SUM_CHUNK_BASES = 1 << 20
# Bases of bedGraph rows spread per chunk: a spread base holds about four times the transient bytes of a read
# row (its position, row and value arrays and the store's copies), so a quarter of a read chunk.
_SPREAD_BASES = 1 << 18
# Whole depths take the first of these that holds them all (signed when one is negative); others are float64.
_WHOLE_DEPTH_DTYPES = (np.dtype(np.uint8), np.dtype(np.uint16), np.dtype(np.uint32))
_SIGNED_DEPTH_DTYPES = (np.dtype(np.int8), np.dtype(np.int16), np.dtype(np.int32))
DepthArray: TypeAlias = (
    NDArray[np.uint8] | NDArray[np.uint16] | NDArray[np.uint32]
    | NDArray[np.int8] | NDArray[np.int16] | NDArray[np.int32]
    | NDArray[np.float64]
)
_UNREADABLE = {"code": "INPUT_UNREADABLE"}


def _invalid_depth(reason: str) -> dict[str, str]:
    return {"code": "DEPTH_INVALID", "reason": reason}


def _is_int_like(value: str) -> bool:
    try:
        int(str(value).strip())
    except ValueError:
        return False
    return True


def _is_float_like(value: str) -> bool:
    try:
        float(str(value).strip())
    except ValueError:
        return False
    return True


def checked_depth_file_format(file_format: str) -> DepthFileFormat:
    """``file_format`` as one of ``DEPTH_FILE_FORMATS``; another value is refused."""

    for known in DEPTH_FILE_FORMATS:
        if file_format == known:
            return known
    raise ValidationError(f"Depth file format must be one of {', '.join(DEPTH_FILE_FORMATS)}; got {file_format!r}.")


@dataclass(frozen=True)
class _DepthForm:
    """How the rows of one depth file read: the columns of its form (a ``name`` column is skipped), the lines
    before its first row, and whether its rows stack (``bed``: each adds its value to the bases it covers)."""

    columns: tuple[str, ...]
    skiprows: int
    stacks: bool = False

    @property
    def read_columns(self) -> tuple[str, ...]:
        return tuple(column for column in self.columns if column != "name")


def _is_header(fields: list[str], columns: tuple[str, ...], file_format: DepthFileFormat = "auto") -> bool:
    """Whether the first data line names its columns: a number column of it is not a number (a BED score of ``.``
    is no score, not a name)."""

    return not all(
        (_is_float_like(field) or (file_format == "bed" and field == ".")) if column == "depth" else _is_int_like(field)
        for field, column in zip(fields, columns)
        if column not in _TEXT_COLUMNS
    )


def _fields(line: str) -> list[str]:
    """The tab-separated fields of a line without its line ending: a trailing tab starts an empty field."""

    return line.rstrip("\r\n").split("\t")


def _form_columns(fields: list[str], file_format: DepthFileFormat = "auto") -> tuple[str, ...]:
    """The columns a data line of ``fields`` reads: by its column count (``auto``), else those of the format a line
    of that many fields holds (a BED line of three fields has no name or score)."""

    if file_format == "auto":
        return _FORM_COLUMNS[min(len(fields), 4)]
    columns, _needed = _FORMAT_COLUMNS[file_format]
    return columns if file_format != "bed" else columns[: max(3, min(len(fields), len(columns)))]


def _depth_form(path: Path, file_format: DepthFileFormat = "auto") -> _DepthForm:
    """The form of a depth file from its first data line, after any comment, track and browser lines.

    By its column count (``auto``), that line has two columns, three, or four or
    more; a stated format reads its own columns of each row (the first value
    column of a several-sample ``samtools`` file) and refuses a line of fewer. A
    comment line that reads as a row is one (a reference may be named ``#1``).
    The rows after the header ``samtools depth -H`` writes are samtools depth
    rows, however many samples they hold, and their first sample is read. The
    first data line is a header line, skipped too, when one of its number
    columns is not a number. An empty file reads as three columns, as before (a
    stated format: its own).
    """

    skipped = 0
    samtools = False
    try:
        with open_input_text(path, encoding="utf-8") as handle:
            for line in handle:
                fields = _fields(line)
                if not _PREAMBLE_LINE.match(line) or (
                    len(fields) >= 2 and not _is_header(fields, _form_columns(fields, file_format), file_format)
                ):
                    break
                samtools = samtools or line.startswith(_SAMTOOLS_HEADER)
                skipped += 1
            else:
                line = ""
    except (OSError, *GZIP_ERRORS) as exc:
        raise ParseError(
            f"Unable to read depth file '{path}': {exc}", diagnostic=_UNREADABLE
        ) from exc

    stacks = file_format == "bed"
    if file_format != "auto":
        columns, needed = _FORMAT_COLUMNS[file_format]
        if not line:
            return _DepthForm(columns[:needed], skipped, stacks)
        if len(fields) < needed:
            raise ParseError(
                f"Depth file '{path}' read as {file_format} must have {needed} or more tab-separated columns "
                f"({', '.join(columns[:needed])}).",
                diagnostic={"code": "TABLE_INVALID", "reason": "DEPTH_COLUMNS"},
            )
        columns = _form_columns(fields, file_format)
        return _DepthForm(columns, skipped + int(_is_header(fields, columns, file_format)), stacks)
    if not line:
        return _DepthForm(tuple(DEPTH_COLUMNS), skipped)
    if len(fields) < 2:
        raise ParseError(
            f"Depth file '{path}' must have two tab-separated columns (position, value), three "
            "(reference, position, value) or four (reference, start, end, value).",
            diagnostic={"code": "TABLE_INVALID", "reason": "DEPTH_COLUMNS"},
        )
    columns = tuple(DEPTH_COLUMNS) if samtools else _form_columns(fields)
    return _DepthForm(columns, skipped + int(_is_header(fields, columns)))


def _depth_tsv_cache_key(path: str) -> tuple[str, int, int]:
    depth_path = Path(path)
    if not depth_path.exists():
        raise ParseError(f"Depth file does not exist: {path}", diagnostic=_UNREADABLE)
    if not depth_path.is_file():
        raise ParseError(f"Depth path is not a file: {path}", diagnostic=_UNREADABLE)
    try:
        file_stat = depth_path.stat()
    except OSError as exc:
        raise ParseError(
            f"Unable to read depth file '{path}': {exc}", diagnostic=_UNREADABLE
        ) from exc
    try:
        resolved_path = str(depth_path.resolve())
    except OSError:
        resolved_path = str(depth_path.absolute())
    return resolved_path, int(file_stat.st_mtime_ns), int(file_stat.st_size)


@lru_cache(maxsize=_DEPTH_TSV_CACHE_MAXSIZE)
def _read_depth_tsv_cached(resolved_path: str, mtime_ns: int, size: int, file_format: DepthFileFormat) -> DepthCoverage:
    return _parse_depth_file(resolved_path, file_format)


@dataclass(frozen=True, eq=False)
class DepthCoverage:
    """The parsed form of depth rows: per reference, the depth at each 1-based source position.

    ``depths[name][position - 1]`` is the mean of the rows for that position
    (their sum, in row order, divided by their count; a ``bed`` file's sum) and
    0 where no row names it; a bedGraph or BED row names the positions
    ``start + 1`` to ``end``. Names keep
    the order of their first row; a file without reference names (two columns)
    has one reference named ``""``. Each array is read-only and
    has the smallest dtype that holds every value exactly: ``uint8``,
    ``uint16`` or ``uint32`` for whole depths, ``int8``, ``int16`` or ``int32``
    when one is negative, ``float64`` otherwise.
    ``extents[name]`` is the first and the last position the rows name: the
    part of the reference the file covers.
    """

    depths: Mapping[str, DepthArray]
    extents: Mapping[str, tuple[int, int]]

    @property
    def nbytes(self) -> int:
        return sum(int(array.nbytes) for array in self.depths.values())


def _smallest_exact_dtype(
    values: NDArray[np.float64], dtype: np.dtype, held: tuple[float, float] = (0.0, 0.0)
) -> np.dtype:
    """The first whole dtype of the ladder (signed when a value is negative) that holds ``values`` and the range
    ``held`` of the values already stored in ``dtype``, exactly; ``float64`` once a value is not whole."""

    if dtype == np.float64 or not values.size:
        return dtype
    whole = np.isfinite(values) & (values == np.floor(values))
    if not whole.all():
        return np.dtype(np.float64)
    smallest, largest = min(held[0], float(values.min())), max(held[1], float(values.max()))
    for candidate in _WHOLE_DEPTH_DTYPES if smallest >= 0 else _SIGNED_DEPTH_DTYPES:
        limits = np.iinfo(candidate)
        if limits.min <= smallest and largest <= limits.max:
            return candidate
    return np.dtype(np.float64)


class _SortedDepths:
    """One reference's rows while each comes after the last: depths are stored, not summed."""

    def __init__(self) -> None:
        self.values: DepthArray = np.zeros(0, dtype=np.uint8)
        self.length = 0
        self.first = 0
        # The least and the largest value stored (0 holds where no row is).
        self.held = (0.0, 0.0)

    def add(self, positions: NDArray[np.int64], depths: NDArray[np.float64]) -> bool:
        """Store rows that follow every earlier one; False (nothing stored) when they do not."""

        if positions[0] <= self.length or (positions.size > 1 and bool((np.diff(positions) <= 0).any())):
            return False
        if not self.length:
            self.first = int(positions[0])
        dtype = _smallest_exact_dtype(depths, self.values.dtype, self.held)
        if dtype != self.values.dtype:
            self.values = self.values.astype(dtype)
        self.held = (min(self.held[0], float(depths.min())), max(self.held[1], float(depths.max())))
        last = int(positions[-1])
        if last > self.values.size:
            # Grow in place where the allocator can; the new tail is zero.
            self.values.resize(max(last, self.values.size * 3 // 2), refcheck=False)
        # ``0.0 + depth`` as the row table summed it: -0.0 is stored as 0.0.
        self.values[positions - 1] = depths + 0.0
        self.length = last
        return True

    def finish(self) -> tuple[DepthArray, tuple[int, int]]:
        self.values.resize(self.length, refcheck=False)
        return self.values, (self.first, self.length)


class _SummedDepths:
    """One reference's rows in any order: sums in row order and counts, divided at the end unless the rows stack
    (``bed``: the sum is the value)."""

    def __init__(self, stacks: bool) -> None:
        self.stacks = stacks
        self.sums: NDArray[np.float64] = np.zeros(0, dtype=np.float64)
        self.counts: NDArray[np.uint32] = np.zeros(0, dtype=np.uint32)
        self.first = 0
        self.last = 0

    def add(self, positions: NDArray[np.int64], depths: NDArray[np.float64]) -> None:
        last = int(positions.max())
        first = int(positions.min())
        self.first = first if not self.last else min(self.first, first)
        self.last = max(self.last, last)
        if last > self.sums.size:
            size = max(last, self.sums.size * 3 // 2)
            self.sums.resize(size, refcheck=False)
            self.counts.resize(size, refcheck=False)
        indexes = positions - 1
        np.add.at(self.sums, indexes, depths)
        np.add.at(self.counts, indexes, 1)

    def finish(self) -> tuple[DepthArray, tuple[int, int]]:
        # Up to the last position, as the sorted store keeps: the length never depends on how the rows were chunked.
        self.sums.resize(self.last, refcheck=False)
        self.counts.resize(self.last, refcheck=False)
        means = self.sums
        if not self.stacks:
            present = self.counts > 0
            means[present] /= self.counts[present]
        self.counts = np.zeros(0, dtype=np.uint32)
        dtype = _smallest_exact_dtype(means, np.dtype(np.uint8))
        return (means if dtype == np.float64 else means.astype(dtype)), (self.first, self.last)


def _reference_rows(chunk: DataFrame) -> Iterator[tuple[str, NDArray[np.int64], NDArray[np.float64]]]:
    """The positions and depths of each reference in a chunk, references in the order of their first row."""

    names = chunk["reference_name"]
    positions = chunk["position"].to_numpy(dtype=np.int64)
    depths = chunk["depth"].to_numpy(dtype=np.float64)
    references = [str(name) for name in names.unique()]
    if len(references) == 1:
        yield references[0], positions, depths
        return
    for reference in references:
        mask = (names == reference).to_numpy(dtype=bool)
        yield reference, positions[mask], depths[mask]


def _coverage(references: Mapping[str, tuple[DepthArray, tuple[int, int]]]) -> DepthCoverage:
    for array, _extent in references.values():
        array.flags.writeable = False
    return DepthCoverage(
        {name: array for name, (array, _extent) in references.items()},
        {name: extent for name, (_array, extent) in references.items()},
    )


def _summed_coverage(chunks: Iterable[DataFrame], *, stacks: bool = False) -> DepthCoverage:
    references: dict[str, _SummedDepths] = {}
    for chunk in chunks:
        for reference, positions, depths in _reference_rows(chunk):
            references.setdefault(reference, _SummedDepths(stacks)).add(positions, depths)
    return _coverage({name: summed.finish() for name, summed in references.items()})


def _sorted_coverage(chunks: Iterable[DataFrame]) -> DepthCoverage | None:
    """The coverage of a file whose rows of each reference come in increasing positions, else None."""

    references: dict[str, _SortedDepths] = {}
    for chunk in chunks:
        for reference, positions, depths in _reference_rows(chunk):
            if not references.setdefault(reference, _SortedDepths()).add(positions, depths):
                return None
    return _coverage({name: stored.finish() for name, stored in references.items()})


def _named_rows(chunk: DataFrame) -> DataFrame:
    """A chunk of rows with its reference names: a file without them names its one reference ``""``."""

    if "reference_name" not in chunk:
        codes = np.zeros(len(chunk), dtype=np.int8)
        chunk.insert(0, "reference_name", pd.Categorical.from_codes(codes, categories=[_UNNAMED_REFERENCE]))
    return chunk


def _interval_bases(rows: DataFrame) -> Iterator[DataFrame]:
    """The per-base rows of bedGraph or BED rows: bases ``start + 1`` to ``end`` of a row take its value.

    The bases are spread in row order, ``_SPREAD_BASES`` at a time, so a long row
    spans several chunks and each chunk stays small beside the stored depths.
    """

    names = rows["reference_name"]
    codes = names.cat.codes.to_numpy()
    values = rows["depth"].to_numpy(dtype=np.float64)
    starts = rows["start"].to_numpy(dtype=np.int64)
    lengths = rows["end"].to_numpy(dtype=np.int64) - starts
    # Row i spreads bases spread_ends[i] - lengths[i] .. spread_ends[i] - 1 of all the chunk's bases.
    spread_ends = np.cumsum(lengths)
    spread_starts = spread_ends - lengths
    # Base k of all the bases, spread from row i, is position k + offsets[i].
    offsets = starts - spread_starts + 1
    total = int(spread_ends[-1]) if spread_ends.size else 0
    for first in range(0, total, _SPREAD_BASES):
        last = min(first + _SPREAD_BASES, total)
        low = int(np.searchsorted(spread_ends, first, side="right"))
        high = int(np.searchsorted(spread_starts, last, side="left"))
        # Repeat counts are the platform's index type: int32 in the Web app's Pyodide, where int64 counts raise.
        taken = (np.minimum(spread_ends[low:high], last) - np.maximum(spread_starts[low:high], first)).astype(np.intp)
        row_of_base = np.repeat(np.arange(low, high), taken)
        positions = np.arange(first, last, dtype=np.int64)
        positions += offsets[row_of_base]
        yield DataFrame(
            {
                "reference_name": pd.Categorical.from_codes(codes[row_of_base], dtype=names.dtype),
                "position": positions,
                "depth": values[row_of_base],
            },
            copy=False,
        )


def _intervals_overlap_or_go_back(rows: DataFrame, last_ends: dict[str, int]) -> bool:
    """Whether a reference's bedGraph rows overlap or go back, in these rows or after the earlier ones.

    ``last_ends`` holds each reference's last end so far and is updated. Rows of
    a reference must come in increasing order without overlapping, as
    ``bedtools genomecov`` writes them and ``bedGraphToBigWig`` requires: then
    no base is named twice.
    """

    names = rows["reference_name"]
    order = np.argsort(names.cat.codes.to_numpy(), kind="stable")
    codes = names.cat.codes.to_numpy()[order]
    starts = rows["start"].to_numpy(dtype=np.int64)[order]
    ends = rows["end"].to_numpy(dtype=np.int64)[order]
    same_reference = codes[1:] == codes[:-1]
    failed = bool((same_reference & (starts[1:] < ends[:-1])).any())
    firsts = np.flatnonzero(np.concatenate(([True], ~same_reference)))
    lasts = np.flatnonzero(np.concatenate((~same_reference, [True])))
    for first, last in zip(firsts.tolist(), lasts.tolist()):
        name = str(names.cat.categories[codes[first]])
        failed = failed or int(starts[first]) < last_ends.get(name, 0)
        last_ends[name] = int(ends[last])
    return failed


def _scored_rows(chunk: DataFrame) -> DataFrame:
    """BED rows with their values: the score, else 1 (no score column, or ``.`` for none)."""

    if "depth" not in chunk:
        chunk["depth"] = 1.0
    else:
        chunk["depth"] = chunk["depth"].fillna(1.0)
    return chunk


def _depth_chunks(resolved_path: str, file_format: DepthFileFormat = "auto") -> Generator[DataFrame, None, None]:
    """The per-base rows of one depth file in file order, a chunk at a time, each one valid.

    A file with an invalid row raises once it is read to the end, so a value
    that does not parse anywhere in it is reported first, and then the first
    failed check in the order of ``_DEPTH_CHECKS``, as for the whole table.
    """

    depth_path = Path(resolved_path)
    form = _depth_form(depth_path, file_format)
    columns = form.read_columns
    # A BED score of ``.`` is no score; nothing else reads as missing there but an empty reference name.
    missing = {"keep_default_na": False, "na_values": {"reference_name": [""], "depth": ["."]}} if form.stacks else {}
    failed: set[str] = set()
    last_ends: dict[str, int] = {}
    try:
        with open_input_text(depth_path, encoding="utf-8") as handle:
            reader = pd.read_csv(
                handle,
                sep="\t",
                header=None,
                names=list(columns),
                usecols=[index for index, column in enumerate(form.columns) if column != "name"],
                skiprows=form.skiprows,
                dtype={column: _COLUMN_DTYPES[column] for column in columns},
                chunksize=_CHUNK_ROWS,
                **missing,
                **EXACT_FLOAT_READ_OPTIONS,
            )
            with reader:
                for read in reader:
                    chunk = _named_rows(_scored_rows(read) if form.stacks else read)
                    failed.update(_failed_depth_checks(chunk))
                    if (
                        not failed
                        and "start" in chunk
                        and not form.stacks
                        and _intervals_overlap_or_go_back(chunk, last_ends)
                    ):
                        failed.add("INTERVAL")
                    if failed:
                        continue
                    if "start" in chunk:
                        yield from _interval_bases(chunk)
                    else:
                        yield chunk
    except (ValueError, *GZIP_ERRORS) as exc:
        raise ParseError(
            f"Unable to parse depth file '{resolved_path}': {exc}",
            diagnostic=_invalid_depth("DEPTH_VALUES"),
        ) from exc
    _raise_first_failed_check(failed)


def _parse_depth_file(resolved_path: str, file_format: DepthFileFormat = "auto") -> DepthCoverage:
    """Parse and validate one depth file into its coverage (see ``parse_depth_tsv``).

    Rows in increasing positions per reference (a samtools depth file) are
    stored as they are read; any other file is read again and summed (a
    ``bed`` file's rows stack: the sum is the value).
    """

    chunks = _depth_chunks(resolved_path, file_format)
    try:
        coverage = _sorted_coverage(chunks)
    finally:
        chunks.close()
    if coverage is None:
        coverage = _summed_coverage(_depth_chunks(resolved_path, file_format), stacks=file_format == "bed")
    return coverage


def clear_depth_tsv_cache() -> None:
    """Clear the process-local cache of parsed depth files used by renders."""

    _read_depth_tsv_cached.cache_clear()


def parse_depth_tsv(path: str, file_format: DepthFileFormat = "auto") -> DepthCoverage:
    """Parse one depth file, read in ``file_format``, without the process cache of renders.

    The Worker's decoded-resource layer owns the coverage it returns and shares
    it, read-only, with every record and render.
    """

    file_format = checked_depth_file_format(file_format)
    resolved_path, _mtime_ns, _size = _depth_tsv_cache_key(path)
    return _parse_depth_file(resolved_path, file_format)


def read_cached_depth_coverage(path: str, file_format: DepthFileFormat = "auto") -> DepthCoverage:
    """The process cache's coverage of one file read in ``file_format`` (CLI and Python API renders), shared
    read-only."""

    file_format = checked_depth_file_format(file_format)
    resolved_path, mtime_ns, size = _depth_tsv_cache_key(path)
    return _read_depth_tsv_cached(resolved_path, mtime_ns, size, file_format)


def read_depth_tsv(path: str, file_format: DepthFileFormat = "auto") -> DataFrame:
    """Read a depth file as a table of its per-base rows: reference name, 1-based position, depth.

    The file is in any form of this module, read by its column count or in
    ``file_format``: two columns (position, value),
    three (reference, position, value) or four or more (bedGraph, whose rows
    give one row per base they cover; a BED file's rows give one row per base
    each covers). Headerless files are expected; comment,
    track and browser lines before the rows are skipped, and a single header
    line is tolerated when one of its number fields is non-numeric. Reference
    names are a ``string`` column. Each call reads the file and returns a new
    table.
    """

    file_format = checked_depth_file_format(file_format)
    resolved_path, _mtime_ns, _size = _depth_tsv_cache_key(path)
    chunks = list(_depth_chunks(resolved_path, file_format))
    table = (
        pd.concat(chunks, ignore_index=True)
        if chunks
        else pd.DataFrame({column: pd.Series(dtype=_COLUMN_DTYPES[column]) for column in DEPTH_COLUMNS})
    )
    table["reference_name"] = table["reference_name"].astype("string")
    return table


def _table_coverage(depth_table: DataFrame) -> DepthCoverage:
    """The coverage of an API table: its depth columns, by name or the first three, made numeric."""

    if set(DEPTH_COLUMNS).issubset(depth_table.columns):
        normalized = depth_table.loc[:, DEPTH_COLUMNS].copy()
    elif len(depth_table.columns) >= 3:
        normalized = depth_table.iloc[:, [0, 1, 2]].copy()
        normalized.columns = DEPTH_COLUMNS
    else:
        raise ValidationError(
            "Depth table must contain at least 3 columns.",
            diagnostic={"code": "TABLE_INVALID", "reason": "THREE_COLUMNS"},
        )

    normalized["reference_name"] = normalized["reference_name"].astype("string")
    try:
        normalized["position"] = pd.to_numeric(normalized["position"], errors="raise").astype("int64")
        normalized["depth"] = pd.to_numeric(normalized["depth"], errors="raise").astype("float64")
    except (TypeError, ValueError) as exc:
        raise ValidationError(
            "Depth table positions and values must be numeric.",
            diagnostic=_invalid_depth("DEPTH_VALUES"),
        ) from exc
    _raise_first_failed_check(_failed_depth_checks(normalized))
    return _summed_coverage([normalized])


# The checks of depth rows in the order a failure is reported: reason, message, and whether rows fail it.
_DEPTH_CHECKS: tuple[tuple[str, str, Callable[[DataFrame], bool]], ...] = (
    (
        "REFERENCE_REQUIRED",
        "Depth table contains missing reference names.",
        lambda rows: bool(rows["reference_name"].isna().any()),
    ),
    (
        "POSITIVE_INTEGER",
        "Depth positions must be 1-based positive integers.",
        lambda rows: "position" in rows and bool((rows["position"] <= 0).any()),
    ),
    (
        "INTERVAL",
        "Depth rows of four or more columns are bedGraph rows (reference, 0-based start, end, value): "
        "each start must be 0 or more and each end greater than its start, and the rows of a reference must "
        "come in increasing order without overlapping (sort -k1,1 -k2,2n). A samtools depth file of several "
        "samples has one value column per sample: keep the reference, position and one sample's column, or "
        "write it with samtools depth -H, whose header line marks its rows.",
        lambda rows: "start" in rows and bool(((rows["start"] < 0) | (rows["end"] <= rows["start"])).any()),
    ),
    (
        "FINITE",
        "Depth values must be finite numbers.",
        lambda rows: not bool(np.isfinite(rows["depth"].to_numpy(dtype=np.float64)).all()),
    ),
)


def _failed_depth_checks(rows: DataFrame) -> set[str]:
    return {reason for reason, _message, fails in _DEPTH_CHECKS if fails(rows)}


def _raise_first_failed_check(failed: set[str]) -> None:
    for reason, message, _fails in _DEPTH_CHECKS:
        if reason in failed:
            raise ValidationError(message, diagnostic=_invalid_depth(reason))


def _select_record_reference(record: SeqRecord, coverage: DepthCoverage) -> str:
    """The reference of the depth rows a record draws: its own ID, else the file's only reference."""

    record_id = str(record.id)
    if record_id in coverage.depths:
        return record_id
    references = list(coverage.depths)
    if len(references) == 1:
        return references[0]
    raise ValidationError(
        "Depth file references do not match record "
        f"'{record_id}'. Available references: {', '.join(references)}",
        diagnostic=_invalid_depth("REFERENCE_MISMATCH"),
    )


def _counted_runs(windows: SlidingWindows, record_length: int, extent: tuple[int, int]) -> list[list[int]]:
    """The local positions ``[start, stop)`` a window's value reads (D-64), as merged runs in local order.

    Every base of the molecule the record draws counts; a base past the drawn record or crop counts only within
    ``extent``, the first to the last position the depth rows of the record name. A position inside both that
    the rows leave out adds 0 to the sum and still counts. A file of the drawn region alone then reads no flank,
    and a file of the whole record reads the true flanks.
    """

    source = windows.source
    first, last = extent
    if source.step > 0:
        extent_indexes = (first - source.base, last - source.base)
    else:
        extent_indexes = (source.base - last, source.base - first)

    def local_run(piece: SourcePiece, low: int, high: int) -> tuple[int, int]:
        """The local positions ``[start, stop)`` of ``piece`` whose source indexes lie in ``low .. high``."""
        if piece.direction > 0:
            start, stop = piece.local + low - piece.index, piece.local + high - piece.index + 1
        else:
            start, stop = piece.local + piece.index - high, piece.local + piece.index - low + 1
        return max(start, piece.local), min(stop, piece.local + piece.length)

    # The source indexes of the drawn bases, and those the depth rows cover.
    known = [extent_indexes]
    for piece in windows.pieces:
        start, stop = max(piece.local, 0), min(piece.local + piece.length, record_length)
        if start < stop:
            ends = (piece.index + piece.direction * (start - piece.local), piece.index + piece.direction * (stop - 1 - piece.local))
            known.append((min(ends), max(ends)))
    counted = [(0, record_length)]
    for piece in windows.pieces:
        for low, high in known:
            start, stop = local_run(piece, low, high)
            if start < stop:
                counted.append((start, stop))
    counted.sort()
    merged: list[list[int]] = []
    for start, stop in counted:
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], stop)
        else:
            merged.append([start, stop])
    return merged


def _window_bases_read(windows: SlidingWindows, counted: list[list[int]]) -> NDArray[np.int64]:
    """How many bases each window's mean is over: its bases in the ``counted`` runs (``_counted_runs``)."""

    bases = np.zeros(windows.starts.size, dtype=np.int64)
    for start, stop in counted:
        bases += np.clip(np.minimum(windows.ends, stop) - np.maximum(windows.starts, start), 0, None)
    return bases


def _shown_pieces(windows: SlidingWindows, depths: DepthArray) -> list[tuple[int, DepthArray]]:
    """The per-base depths the windows read, in local order: for each source piece, its first local position and
    the depths it reads there, backwards when the record reads its source backwards (a reverse complement).
    Positions the depth file does not hold are left out (they read 0)."""

    source = windows.source
    available = int(depths.size)
    pieces: list[tuple[int, DepthArray]] = []
    for piece in windows.pieces:
        # Local position piece.local + k reads source coordinate first + run_step * k.
        first = source.base + source.step * piece.index
        run_step = source.step * piece.direction
        if run_step > 0:
            skip, keep = max(0, 1 - first), min(piece.length, available - first + 1)
        else:
            skip, keep = max(0, first - available), min(piece.length, first)
        if skip >= keep:
            continue
        coordinate = first + run_step * skip
        shown = (
            depths[coordinate - 1 : coordinate - 1 + keep - skip]
            if run_step > 0
            else depths[coordinate - (keep - skip) : coordinate][::-1]
        )
        pieces.append((piece.local + skip, shown))
    return pieces


def _window_depth_maxima(windows: SlidingWindows, depths: DepthArray, counted: list[list[int]]) -> NDArray[np.float64]:
    """The largest per-base depth of each window among the bases its mean counts (``counted``): a counted base the
    file leaves out is 0 (D-24).

    The depths are laid out along the local axis in their own dtype, a base no window counts at the dtype's least
    value; each window is a run of the elementary intervals between window bounds, whose maxima are reduced once
    and then combined two overlapping power-of-two runs at a time (a sparse table), so no window reads its bases
    again. Every index array is ``np.intp``: Pyodide's NumPy is 32-bit.
    """

    starts, ends = windows.starts, windows.ends
    # The local axis the pieces map (``sliding_windows``); the counted runs lie on it.
    low = min(int(starts.min()), 0)
    high = max(int(ends.max()), max((stop for _start, stop in counted), default=0))
    values = np.zeros(high - low, dtype=depths.dtype)
    for first, shown in _shown_pieces(windows, depths):
        values[first - low : first - low + int(shown.size)] = shown
    read = np.zeros(high - low, dtype=bool)
    for start, stop in counted:
        read[start - low : stop - low] = True
    values[~read] = np.iinfo(depths.dtype).min if depths.dtype.kind in "iu" else -np.inf
    bounds, bound_of = np.unique(np.concatenate((starts, ends)) - low, return_inverse=True)
    bounds = bounds.astype(np.intp)
    level = [np.maximum.reduceat(values[: int(bounds[-1])], bounds[:-1])]
    first_interval, last_interval = bound_of[: starts.size].astype(np.intp), bound_of[starts.size :].astype(np.intp)
    spans = last_interval - first_interval
    widest = int(spans.max())
    while (1 << len(level)) <= widest:
        half = 1 << (len(level) - 1)
        level.append(np.maximum(level[-1][:-half], level[-1][half:]))
    powers = np.floor(np.log2(np.maximum(spans, 1))).astype(np.intp)
    maxima = np.empty(starts.size, dtype=np.float64)
    for power in np.unique(powers).tolist():
        rows = np.flatnonzero(powers == power)
        table = level[power]
        maxima[rows] = np.maximum(table[first_interval[rows]], table[last_interval[rows] - (1 << power)])
    return maxima


def _window_depth_sums(windows: SlidingWindows, depths: DepthArray) -> NDArray[np.float64]:
    """The depth summed over each window ``windows.starts[i] .. windows.ends[i] - 1`` of the record's local axis.

    Each source piece of the windows (:func:`gbdraw.analysis.windows.sliding_windows`) reads source coordinates
    in a run, backwards when the record reads its source backwards (a reverse complement); positions the depth
    file does not hold add 0. Each sum is ``prefix[end] - prefix[start]`` of the running float64 sums of the
    per-base depths along the local axis. The running sums are kept only at window bounds: the bases are summed
    one chunk of ``_SUM_CHUNK_BASES`` at a time, each chunk carrying on from the last sum of the one before, so
    every sum is the same chain of additions as one pass along the axis.
    """

    starts, ends = windows.starts, windows.ends
    bounds, bound_of = np.unique(np.concatenate((starts, ends)), return_inverse=True)
    prefix = np.zeros(bounds.size, dtype=np.float64)
    pieces = _shown_pieces(windows, depths)
    carry = np.float64(0.0)
    covered = int(bounds[0])
    # The first base starts the chain (adding 0.0 would turn a -0.0 depth into 0.0).
    summed = False
    for low, shown in pieces:
        # Bounds between the last shown base and this piece hold the running sum so far.
        prefix[np.searchsorted(bounds, covered, side="right") : np.searchsorted(bounds, low, side="right")] = carry
        # prefix[low + 1 + i] is the running sum up to shown[i].
        buffer = np.empty(min(_SUM_CHUNK_BASES, int(shown.size)), dtype=np.float64)
        for offset in range(0, int(shown.size), buffer.size):
            sums = buffer[: min(buffer.size, int(shown.size) - offset)]
            sums[:] = shown[offset : offset + sums.size]
            if summed:
                sums[0] += carry
            summed = True
            np.cumsum(sums, out=sums)
            carry = sums[-1]
            chunk_first = low + 1 + offset
            inside = slice(*np.searchsorted(bounds, (chunk_first, chunk_first + sums.size)))
            prefix[inside] = sums[bounds[inside] - chunk_first]
        covered = low + int(shown.size)
    # Past the last shown base the running sum stays at its total.
    if pieces:
        prefix[np.searchsorted(bounds, covered, side="right") :] = carry
    return prefix[bound_of[starts.size :]] - prefix[bound_of[: starts.size]]


def _validate_depth_range(min_depth: float | None, max_depth: float | None) -> None:
    """A value range may be negative; its minimum is at most its maximum."""

    if min_depth is not None and max_depth is not None and min_depth > max_depth:
        raise ValidationError("min_depth must be <= max_depth.")


def _refuse_negative_values(depths: DepthArray) -> None:
    """Refuse depths of which one is negative: a log-scaled axis draws only values of 0 or more."""

    if depths.dtype.kind != "u" and depths.size and float(depths.min()) < 0:
        raise ValidationError("Depth values must be non-negative.", diagnostic=_invalid_depth("NONNEGATIVE"))


def depth_df(
    record: SeqRecord,
    depth_table: DepthCoverage | DataFrame,
    window: int,
    step: int,
    *,
    circular: bool,
    min_depth: float | None = None,
    max_depth: float | None = None,
    agg: DepthAggregation = "mean",
    log_scale: bool = False,
) -> DataFrame:
    """The value of each sliding window of a record (``agg``: the mean depth, the largest or the sum): one row per
    sample.

    ``depth_table`` is a parsed file or a table of rows (an API table, which is
    converted here). Positions are 1-based source coordinates of the named
    record; the record coordinate map maps them to local positions, flipped for
    a reverse complement. The samples (``position``, 0-based local) and their
    windows are those of :func:`gbdraw.analysis.windows.sliding_windows`, centred
    on the sample; ``circular`` says whether the molecule is circular. A crop
    reads its flanking positions of the source where the depth rows cover them;
    a mean is over the bases a window reads (``_counted_runs``), the largest
    value is among the same bases (a base the rows leave out is 0) and the sum
    is over every base the window reads; each value is clipped to
    ``min_depth``..``max_depth``. A ``log_scale`` track refuses a negative
    value of the record's reference. Where a value is drawn on the track,
    linear or log scaled, is the drawers' (``_prepare_depth_plot_dataframe``).
    """

    if window <= 0:
        raise ValidationError("window must be > 0.")
    if step <= 0:
        raise ValidationError("step must be > 0.")
    if agg not in DEPTH_AGGREGATIONS:
        raise ValidationError(f"Depth agg must be one of {', '.join(DEPTH_AGGREGATIONS)}; got {agg!r}.")
    _validate_depth_range(min_depth, max_depth)

    coverage_source = depth_table if isinstance(depth_table, DepthCoverage) else _table_coverage(depth_table)
    record_len = len(record.seq)
    if record_len <= 0:
        return pd.DataFrame(columns=["position", "depth"])

    windows = sliding_windows(record, int(window), int(step), circular=circular)
    reference = _select_record_reference(record, coverage_source)
    depths = coverage_source.depths[reference]
    if log_scale:
        _refuse_negative_values(depths)
    counted = _counted_runs(windows, record_len, coverage_source.extents[reference])
    if agg == "max":
        means = _window_depth_maxima(windows, depths, counted)
    elif agg == "sum":
        means = _window_depth_sums(windows, depths)
    else:
        means = _window_depth_sums(windows, depths) / _window_bases_read(windows, counted)

    if min_depth is not None or max_depth is not None:
        lower = float(min_depth) if min_depth is not None else None
        upper = float(max_depth) if max_depth is not None else None
        means = np.clip(means, lower, upper)

    return pd.DataFrame(
        {
            "position": windows.samples,
            "depth": means,
        }
    )


__all__ = [
    "DEPTH_AGGREGATIONS",
    "DEPTH_COLUMNS",
    "DEPTH_FILE_FORMATS",
    "DepthAggregation",
    "DepthCoverage",
    "DepthFileFormat",
    "checked_depth_file_format",
    "clear_depth_tsv_cache",
    "depth_df",
    "parse_depth_tsv",
    "read_cached_depth_coverage",
    "read_depth_tsv",
]
