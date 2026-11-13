#!/usr/bin/env python
# coding: utf-8

"""Depth coverage parsing and binning helpers."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Callable, Generator, Iterable, Iterator, Mapping, TypeAlias

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
_DEPTH_TSV_CACHE_MAXSIZE = 8
# Rows per chunk read: a chunk's transient columns stay small beside the stored depths.
_CHUNK_ROWS = 1 << 20
# Bases summed per chunk when binning: the running sums stay small beside the stored depths.
_SUM_CHUNK_BASES = 1 << 20
_EMPTY_DEPTH_DTYPES = {"reference_name": "category", "position": "int64", "depth": "float64"}
# Whole depths take the first of these that holds the largest; others are float64.
_WHOLE_DEPTH_DTYPES = (np.dtype(np.uint8), np.dtype(np.uint16), np.dtype(np.uint32))
DepthArray: TypeAlias = NDArray[np.uint8] | NDArray[np.uint16] | NDArray[np.uint32] | NDArray[np.float64]
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


def _first_line_has_header(path: Path) -> bool:
    try:
        with open_input_text(path, encoding="utf-8") as handle:
            first_line = handle.readline()
    except (OSError, *GZIP_ERRORS) as exc:
        raise ParseError(
            f"Unable to read depth file '{path}': {exc}", diagnostic=_UNREADABLE
        ) from exc

    if not first_line:
        return False
    fields = first_line.rstrip("\n\r").split("\t")
    if len(fields) < 3:
        raise ParseError(
            f"Depth file '{path}' must contain at least 3 tab-separated columns.",
            diagnostic={"code": "TABLE_INVALID", "reason": "THREE_COLUMNS"},
        )
    return not (_is_int_like(fields[1]) and _is_float_like(fields[2]))


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
def _read_depth_tsv_cached(resolved_path: str, mtime_ns: int, size: int) -> DepthCoverage:
    return _parse_depth_file(resolved_path)


@dataclass(frozen=True, eq=False)
class DepthCoverage:
    """The parsed form of depth rows: per reference, the depth at each 1-based source position.

    ``depths[name][position - 1]`` is the mean of the rows for that position
    (their sum, in row order, divided by their count) and 0 where no row names
    it. Names keep the order of their first row. Each array is read-only and
    has the smallest dtype that holds every value exactly: ``uint8``,
    ``uint16`` or ``uint32`` for whole depths, ``float64`` otherwise.
    ``extents[name]`` is the first and the last position the rows name: the
    part of the reference the file covers.
    """

    depths: Mapping[str, DepthArray]
    extents: Mapping[str, tuple[int, int]]

    @property
    def nbytes(self) -> int:
        return sum(int(array.nbytes) for array in self.depths.values())


def _smallest_exact_dtype(values: NDArray[np.float64], dtype: np.dtype) -> np.dtype:
    """The first dtype of the ladder from ``dtype`` that holds ``values`` exactly."""

    if dtype == np.float64 or not values.size:
        return dtype
    whole = np.isfinite(values) & (values == np.floor(values))
    if not whole.all():
        return np.dtype(np.float64)
    largest = float(values.max())
    for candidate in _WHOLE_DEPTH_DTYPES[_WHOLE_DEPTH_DTYPES.index(dtype):]:
        if largest <= np.iinfo(candidate).max:
            return candidate
    return np.dtype(np.float64)


class _SortedDepths:
    """One reference's rows while each comes after the last: depths are stored, not summed."""

    def __init__(self) -> None:
        self.values: DepthArray = np.zeros(0, dtype=np.uint8)
        self.length = 0
        self.first = 0

    def add(self, positions: NDArray[np.int64], depths: NDArray[np.float64]) -> bool:
        """Store rows that follow every earlier one; False (nothing stored) when they do not."""

        if positions[0] <= self.length or (positions.size > 1 and bool((np.diff(positions) <= 0).any())):
            return False
        if not self.length:
            self.first = int(positions[0])
        dtype = _smallest_exact_dtype(depths, self.values.dtype)
        if dtype != self.values.dtype:
            self.values = self.values.astype(dtype)
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
    """One reference's rows in any order: sums in row order and counts, divided at the end."""

    def __init__(self) -> None:
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
        present = self.counts > 0
        means = self.sums
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


def _summed_coverage(chunks: Iterable[DataFrame]) -> DepthCoverage:
    references: dict[str, _SummedDepths] = {}
    for chunk in chunks:
        for reference, positions, depths in _reference_rows(chunk):
            references.setdefault(reference, _SummedDepths()).add(positions, depths)
    return _coverage({name: summed.finish() for name, summed in references.items()})


def _sorted_coverage(chunks: Iterable[DataFrame]) -> DepthCoverage | None:
    """The coverage of a file whose rows of each reference come in increasing positions, else None."""

    references: dict[str, _SortedDepths] = {}
    for chunk in chunks:
        for reference, positions, depths in _reference_rows(chunk):
            if not references.setdefault(reference, _SortedDepths()).add(positions, depths):
                return None
    return _coverage({name: stored.finish() for name, stored in references.items()})


def _depth_chunks(resolved_path: str) -> Generator[DataFrame, None, None]:
    """The rows of one depth file in file order, a chunk at a time, each one valid.

    A file with an invalid row raises once it is read to the end, so a value
    that does not parse anywhere in it is reported first, and then the first
    failed check in the order of ``_DEPTH_CHECKS``, as for the whole table.
    """

    depth_path = Path(resolved_path)
    skiprows = 1 if _first_line_has_header(depth_path) else 0
    failed: set[str] = set()
    try:
        with open_input_text(depth_path, encoding="utf-8") as handle:
            reader = pd.read_csv(
                handle,
                sep="\t",
                header=None,
                names=DEPTH_COLUMNS,
                usecols=[0, 1, 2],
                skiprows=skiprows,
                dtype={
                    # A few names repeat on every row: codes, not one string per row.
                    "reference_name": "category",
                    "position": "int64",
                    "depth": "float64",
                },
                chunksize=_CHUNK_ROWS,
                **EXACT_FLOAT_READ_OPTIONS,
            )
            with reader:
                for chunk in reader:
                    failed.update(_failed_depth_checks(chunk))
                    if not failed:
                        yield chunk
    except (ValueError, *GZIP_ERRORS) as exc:
        raise ParseError(
            f"Unable to parse depth file '{resolved_path}': {exc}",
            diagnostic=_invalid_depth("DEPTH_VALUES"),
        ) from exc
    _raise_first_failed_check(failed)


def _parse_depth_file(resolved_path: str) -> DepthCoverage:
    """Parse and validate one depth file into its coverage (see ``parse_depth_tsv``).

    Rows in increasing positions per reference (a samtools depth file) are
    stored as they are read; any other file is read again and summed.
    """

    chunks = _depth_chunks(resolved_path)
    try:
        coverage = _sorted_coverage(chunks)
    finally:
        chunks.close()
    if coverage is None:
        coverage = _summed_coverage(_depth_chunks(resolved_path))
    return coverage


def clear_depth_tsv_cache() -> None:
    """Clear the process-local cache of parsed depth files used by renders."""

    _read_depth_tsv_cached.cache_clear()


def parse_depth_tsv(path: str) -> DepthCoverage:
    """Parse one depth file without the process cache of renders.

    The Worker's decoded-resource layer owns the coverage it returns and shares
    it, read-only, with every record and render.
    """

    resolved_path, _mtime_ns, _size = _depth_tsv_cache_key(path)
    return _parse_depth_file(resolved_path)


def read_cached_depth_coverage(path: str) -> DepthCoverage:
    """The process cache's coverage of one file (CLI and Python API renders), shared read-only."""

    resolved_path, mtime_ns, size = _depth_tsv_cache_key(path)
    return _read_depth_tsv_cached(resolved_path, mtime_ns, size)


def read_depth_tsv(path: str) -> DataFrame:
    """Read a samtools-depth-like TSV file as a table of its rows.

    The expected columns are reference name, 1-based position, and non-negative
    depth. Headerless files are expected; a single header line is tolerated when
    the position/depth fields on the first line are non-numeric. Reference names
    are a ``string`` column. Each call reads the file and returns a new table.
    """

    resolved_path, _mtime_ns, _size = _depth_tsv_cache_key(path)
    chunks = list(_depth_chunks(resolved_path))
    table = (
        pd.concat(chunks, ignore_index=True)
        if chunks
        else pd.DataFrame({column: pd.Series(dtype=dtype) for column, dtype in _EMPTY_DEPTH_DTYPES.items()})
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
        lambda rows: bool((rows["position"] <= 0).any()),
    ),
    (
        "NONNEGATIVE",
        "Depth values must be non-negative.",
        lambda rows: bool((rows["depth"] < 0).any()),
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


def _window_bases_read(windows: SlidingWindows, record_length: int, extent: tuple[int, int]) -> NDArray[np.int64]:
    """How many bases each window's mean is over (D-64).

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
    bases = np.zeros(windows.starts.size, dtype=np.int64)
    for start, stop in merged:
        bases += np.clip(np.minimum(windows.ends, stop) - np.maximum(windows.starts, start), 0, None)
    return bases


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
    if min_depth is not None and min_depth < 0:
        raise ValidationError("min_depth must be >= 0.")
    if max_depth is not None and max_depth < 0:
        raise ValidationError("max_depth must be >= 0.")
    if min_depth is not None and max_depth is not None and min_depth > max_depth:
        raise ValidationError("min_depth must be <= max_depth.")


def depth_df(
    record: SeqRecord,
    depth_table: DepthCoverage | DataFrame,
    window: int,
    step: int,
    *,
    circular: bool,
    min_depth: float | None = None,
    max_depth: float | None = None,
) -> DataFrame:
    """The mean depth of each sliding window of a record: one row per sample.

    ``depth_table`` is a parsed file or a table of rows (an API table, which is
    converted here). Positions are 1-based source coordinates of the named
    record; the record coordinate map maps them to local positions, flipped for
    a reverse complement. The samples (``position``, 0-based local) and their
    windows are those of :func:`gbdraw.analysis.windows.sliding_windows`, centred
    on the sample; ``circular`` says whether the molecule is circular. A crop
    reads its flanking positions of the source where the depth rows cover them;
    a mean is over the bases a window reads (``_window_bases_read``), clipped
    to ``min_depth``..``max_depth``. Where a mean is drawn on the track, linear
    or log scaled, is the drawers' (``_prepare_depth_plot_dataframe``).
    """

    if window <= 0:
        raise ValidationError("window must be > 0.")
    if step <= 0:
        raise ValidationError("step must be > 0.")
    _validate_depth_range(min_depth, max_depth)

    coverage_source = depth_table if isinstance(depth_table, DepthCoverage) else _table_coverage(depth_table)
    record_len = len(record.seq)
    if record_len <= 0:
        return pd.DataFrame(columns=["position", "depth"])

    windows = sliding_windows(record, int(window), int(step), circular=circular)
    reference = _select_record_reference(record, coverage_source)
    depths = coverage_source.depths[reference]
    bases = _window_bases_read(windows, record_len, coverage_source.extents[reference])
    means = _window_depth_sums(windows, depths) / bases

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
    "DEPTH_COLUMNS",
    "DepthCoverage",
    "clear_depth_tsv_cache",
    "depth_df",
    "parse_depth_tsv",
    "read_cached_depth_coverage",
    "read_depth_tsv",
]
