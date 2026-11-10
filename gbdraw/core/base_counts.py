"""How bases are counted: the one counting rule of every GC count, skew count and letter check.

Case is ignored, so soft-masked bases count, and U is counted as T (D-26). The
known bases are A, C, G and T; every other byte (N, the other IUPAC codes, a gap)
is an unknown base, which no content or skew value counts (OV-2799). A sequence
is a nucleotide sequence when it holds IUPAC nucleotide codes only, in either case
(D-17), and a protein sequence when it holds letters and ``*`` only (PD-OI-117).

Bases are counted in rows of ``ROW_BASES`` bases: row ``k`` of a sequence holds
the counts of its bases ``k * ROW_BASES`` up to the next row (:func:`count_rows`,
read straight from the bytes, never translated or copied; a byte that is not a
base, such as a line end, belongs to no row). A count of any span is then whole
rows summed and the bases of its edge rows counted, and a span whose ends are
multiples of ``ROW_BASES`` reads no base: the Auto GC windows of a record of
1 Mbp or more have every bound there (``gbdraw.auto_sizes.sliding_window``).
:mod:`gbdraw.io.record_sequences` keeps the rows of every record sequence it reads.
"""

from __future__ import annotations

import functools
import string
from typing import Final, NamedTuple

import numpy as np
from numpy.typing import NDArray

ROW_BASES: Final = 1000  # the bases of one count row: rows begin at every multiple of it
SLICE_BYTES: Final = 64 << 10  # the bytes one histogram call reads

NUCLEOTIDE_CODES: Final = b"ACGTURYSWKMBDHVNacgturyswkmbdhvn"
PROTEIN_CODES: Final = string.ascii_letters.encode("ascii") + b"*"

# What a sequence holds, from its bytes: nucleotide codes only; else letters and "*" only; else another byte.
NUCLEOTIDE: Final = 0
PROTEIN: Final = 1
OTHER: Final = 2

# The columns of a row's counts: the known bases, then all bases, then the bases outside each alphabet.
A, C, G, T, BASES, BEYOND_NUCLEOTIDE, BEYOND_PROTEIN = range(7)
KNOWN: Final = slice(A, T + 1)
_COLUMNS: Final = 7
# The counting rule of A, C, G and T, in column order: the one owner every count below derives from.
_LETTERS: Final = ((A, b"Aa"), (C, b"Cc"), (G, b"Gg"), (T, b"TtUu"))


def count_matrix(deleted: bytes = b"") -> NDArray[np.int64]:
    """What each byte value adds to each column of a row's counts; a byte in ``deleted`` is not a base and adds nothing."""

    matrix = np.zeros((256, _COLUMNS), dtype=np.int64)
    for column, letters in _LETTERS:
        matrix[list(letters), column] = 1
    kept = np.ones(256, dtype=bool)
    kept[list(deleted)] = False
    matrix[kept, BASES] = 1
    for column, codes in ((BEYOND_NUCLEOTIDE, NUCLEOTIDE_CODES), (BEYOND_PROTEIN, PROTEIN_CODES)):
        matrix[kept, column] = 1
        matrix[list(codes), column] = 0
    matrix[~kept, :] = 0
    return matrix


BASES_MATRIX: Final = count_matrix()  # bytes that are all bases

# The bytes of each column's letters in a row, one after another, and where each column's begin.
_LETTER_BYTES: Final = np.frombuffer(b"".join(letters for _, letters in _LETTERS), dtype=np.uint8).astype(np.intp)
_LETTER_COLUMNS: Final = np.cumsum([0] + [len(letters) for _, letters in _LETTERS[:-1]]).astype(np.intp)


class CountedRows(NamedTuple):
    """What :func:`count_rows` counted in some bytes.

    ``rows`` (A, C, G, T; U as T, either case): row 0 continues the row the bytes begin in (zero when they begin
    one), then one row per row that begins in them, whose first base is at offset ``starts[j]`` of the bytes;
    ``histograms`` how often each byte value occurs in each ``SLICE_BYTES`` of them; ``bases`` the bases they hold.
    """

    rows: NDArray[np.int64]
    starts: NDArray[np.int64]
    histograms: NDArray[np.int64]
    bases: int


@functools.lru_cache(maxsize=8)
def _deletion(deleted: bytes) -> tuple[NDArray[np.bool_], int, NDArray[np.bool_]]:
    """The deletion table, the highest deleted value, and the values up to it that are bases."""

    table = np.zeros(256, dtype=bool)
    table[list(deleted)] = True
    top = max(deleted, default=-1)
    low_bases = ~table
    low_bases[top + 1:] = False
    return table, top, low_bases


def _slice_rows(
    values: NDArray[np.uint8], gone: NDArray[np.intp], phase: int,
) -> tuple[NDArray[np.intp], NDArray[np.int64], NDArray[np.int64]]:
    """The offsets of the rows that begin in a slice whose deleted bytes are at ``gone``, ``phase`` bases into a
    row, the letters of each row the slice touches (the first continues a row), and the slice's histogram."""

    firsts = np.arange(-phase % ROW_BASES, len(values) - len(gone), ROW_BASES)
    # The ``b``-th base of the slice is byte ``b`` plus the deleted bytes before it.
    begins = firsts + np.searchsorted(gone - np.arange(len(gone)), firsts, side="right") if len(gone) else firsts
    cuts = np.concatenate(([0], begins)) if not len(begins) or begins[0] else begins
    key = np.uint16 if len(cuts) <= 256 else np.intp  # a key: a byte's row in the slice, then its value
    sizes = np.empty(len(cuts), dtype=np.intp)
    np.subtract(cuts[1:], cuts[:-1], out=sizes[:-1])
    sizes[-1] = len(values) - cuts[-1]
    keys = np.repeat(np.arange(0, len(cuts) * 256, 256, dtype=key), sizes)
    keys += values
    pieces = np.bincount(keys, minlength=len(cuts) * 256).reshape(len(cuts), 256)
    letters = np.add.reduceat(pieces[:, _LETTER_BYTES], _LETTER_COLUMNS, axis=1)
    if len(cuts) > len(begins):
        return begins, letters, pieces.sum(axis=0)
    return begins, np.concatenate((np.zeros((1, 4), dtype=np.int64), letters)), pieces.sum(axis=0)


def count_rows(data: bytes | memoryview, before: int, deleted: bytes = b"") -> CountedRows:
    """Count the bases of ``data`` in rows of ``ROW_BASES``, ``before`` bases of the sequence before it.

    The one producer of count rows. A byte in ``deleted`` is not a base. Each slice is one histogram call whose
    keys say the row of each byte: a row's first byte is found from the positions of the bytes up to the highest
    deleted value, taken as deleted; when the slice's histogram shows one of them is a base, the slice is counted
    again with the deletion table. A row's counts are read from its letters' columns, so the work does not grow
    with the number of rows.
    """

    table, top, low_bases = _deletion(deleted)
    rows = [np.zeros((1, 4), dtype=np.int64)]
    starts: list[NDArray[np.int64]] = []
    histograms: list[NDArray[np.int64]] = []
    bases = before
    for low in range(0, len(data), SLICE_BYTES):
        values = np.frombuffer(data, dtype=np.uint8, count=min(SLICE_BYTES, len(data) - low), offset=low)
        gone = np.flatnonzero(values <= top) if top >= 0 else np.zeros(0, dtype=np.intp)
        begins, letters, histogram = _slice_rows(values, gone, bases)
        if histogram[low_bases].any():  # a byte taken as deleted is a base
            gone = gone[table[values[gone]]]
            begins, letters, histogram = _slice_rows(values, gone, bases)
        rows[-1][-1] += letters[0]
        if len(letters) > 1:
            rows.append(letters[1:])
        starts.append(begins + low)
        histograms.append(histogram)
        bases += len(values) - len(gone)
    return CountedRows(  # int64 on every platform: numpy's own integers are 32 bits in Pyodide
        (np.concatenate(rows) if len(rows) > 1 else rows[0]).astype(np.int64, copy=False),
        (np.concatenate(starts) if starts else np.zeros(0, dtype=np.int64)).astype(np.int64, copy=False),
        (np.stack(histograms) if histograms else np.zeros((0, 256), dtype=np.int64)).astype(np.int64, copy=False),
        bases - before,
    )


class RowCounts:
    """The count rows of one sequence, added in order (:func:`count_rows`): a row begun by some bytes is ended
    by the next."""

    def __init__(self) -> None:
        self._rows: list[NDArray[np.uint16]] = []

    def add(self, rows: NDArray[np.int64]) -> None:
        if self._rows:
            self._rows[-1][-1] += rows[0].astype(np.uint16)
        else:
            assert not rows[0].any()  # bytes before the first base hold no letter
        if len(rows) > 1:
            self._rows.append(rows[1:].astype(np.uint16))

    def joined(self) -> NDArray[np.uint16]:
        """Row ``k``: the counts of A, C, G and T in bases ``k * ROW_BASES`` up to the next row (``uint16``)."""

        joined = np.concatenate(self._rows) if self._rows else np.zeros((0, 4), dtype=np.uint16)
        self._rows = [joined]
        return joined


def byte_counts(data: bytes | memoryview, matrix: NDArray[np.int64] = BASES_MATRIX) -> NDArray[np.int64]:
    """The counts of all of ``data`` (7 columns: A, C, G, T, all bases, bases outside each alphabet)."""

    return np.bincount(np.frombuffer(data, dtype=np.uint8), minlength=256) @ matrix


def alphabet_of(counts: NDArray[np.int64]) -> int:
    """What the bases of ``counts`` (one row or several) hold: ``NUCLEOTIDE``, ``PROTEIN`` or ``OTHER``."""

    total = counts.reshape(-1, _COLUMNS).sum(axis=0)
    return OTHER if total[BEYOND_PROTEIN] else PROTEIN if total[BEYOND_NUCLEOTIDE] else NUCLEOTIDE


# Multiplying a word of eight 0/1 bytes by this constant sums them into its top byte.
_BYTE_SUM = np.uint64(0x0101010101010101)
_TOP_BYTE = np.uint64(56)
_EIGHT = np.arange(8)
_LOWER = np.uint8(0x20)
# Each letter once in lower case (``| 0x20`` folds a byte's case before it is compared).
_KNOWN_LETTERS: Final = tuple((column, bytes((letter,))) for column, letters in _LETTERS
                              for letter in sorted(set(letters.lower())))


def known_prefix_counts(bases: bytes, positions: NDArray[np.intp]) -> NDArray[np.int64]:
    """Column ``k``: the counts of A, C, G and T (rows) in ``bases[:positions[k]]`` (``0 <= positions[k] <= len``).

    One pass over ``bases`` for every position together: each letter's mask is read as 8-byte words whose set
    bytes one multiplication sums, so a position's count is the words before it and the bytes of its own word
    before it. ``| 0x20`` folds the case of a letter and maps no other byte to a, c, g, t or u.
    """

    counts = np.zeros((4, len(positions)), dtype=np.int64)
    if not len(bases) or not len(positions):
        return counts
    padded = np.zeros(-(-len(bases) // 8) * 8, dtype=np.uint8)
    np.bitwise_or(np.frombuffer(bases, dtype=np.uint8), _LOWER, out=padded[: len(bases)])
    word, place = np.divmod(positions.astype(np.intp), 8)
    in_word = _EIGHT < place[:, None]
    mask = np.empty(len(padded), dtype=bool)
    words = np.empty(len(padded) // 8, dtype=np.uint64)
    before = np.zeros(len(words) + 1, dtype=np.int64)
    for row, letter in _KNOWN_LETTERS:
        np.equal(padded, ord(letter), out=mask)
        mask[len(bases):] = False
        np.multiply(mask.view(np.uint64), _BYTE_SUM, out=words)
        np.right_shift(words, _TOP_BYTE, out=words)
        np.cumsum(words.view(np.int64), out=before[1:])
        safe = np.minimum(word, len(words) - 1)
        partial = np.count_nonzero(mask.reshape(-1, 8)[safe] & in_word, axis=1)
        counts[row] += before[word] + np.where(word < len(words), partial, 0)
    return counts


__all__ = [
    "A",
    "BASES",
    "BASES_MATRIX",
    "BEYOND_NUCLEOTIDE",
    "BEYOND_PROTEIN",
    "C",
    "G",
    "KNOWN",
    "NUCLEOTIDE",
    "NUCLEOTIDE_CODES",
    "OTHER",
    "PROTEIN",
    "PROTEIN_CODES",
    "ROW_BASES",
    "SLICE_BYTES",
    "T",
    "alphabet_of",
    "CountedRows",
    "RowCounts",
    "byte_counts",
    "count_matrix",
    "count_rows",
    "known_prefix_counts",
]
