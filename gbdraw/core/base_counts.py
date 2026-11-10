"""How bases are counted: the one counting rule of every GC count, skew count and letter check.

Case is ignored, so soft-masked bases count, and U is counted as T (D-26). The
known bases are A, C, G and T; every other byte (N, the other IUPAC codes, a gap)
is an unknown base, which no content or skew value counts (OV-2799). A sequence
is a nucleotide sequence when it holds IUPAC nucleotide codes only, in either case
(D-17), and a protein sequence when it holds letters and ``*`` only (PD-OI-117).

Bytes are counted in rows of at most ``ROW_BYTES`` (one histogram of the byte
values per row, read straight from the bytes, never translated or copied), so a
count of any span is whole rows summed and its two edge rows counted.
:mod:`gbdraw.io.record_sequences` keeps the rows of every record sequence it reads.
"""

from __future__ import annotations

import string
from typing import Final

import numpy as np
from numpy.typing import NDArray

ROW_BYTES: Final = 4096  # most bytes of one counted row
SLICE_BYTES: Final = 64 << 10  # the bytes one histogram call reads: whole rows

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


def count_matrix(deleted: bytes = b"") -> NDArray[np.int64]:
    """What each byte value adds to each column of a row's counts; a byte in ``deleted`` is not a base and adds nothing."""

    matrix = np.zeros((256, _COLUMNS), dtype=np.int64)
    for column, letters in ((A, b"Aa"), (C, b"Cc"), (G, b"Gg"), (T, b"TtUu")):
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

# Row ``i`` of a slice adds ``256 * i`` to its bytes' values, so one histogram call counts every row of the slice.
_ROW_KEYS: Final = (np.arange(SLICE_BYTES, dtype=np.uint32) // ROW_BYTES * 256).astype(np.uint16)


def byte_histograms(data: bytes | memoryview) -> NDArray[np.int64]:
    """Row ``i``: how often each byte value occurs in bytes ``ROW_BYTES * i`` onwards of ``data`` (rows x 256)."""

    histograms: list[NDArray[np.int64]] = []
    for start in range(0, len(data), SLICE_BYTES):
        size = min(SLICE_BYTES, len(data) - start)
        values = np.frombuffer(data, dtype=np.uint8, count=size, offset=start)
        rows = -(-size // ROW_BYTES)
        histograms.append(np.bincount(_ROW_KEYS[:size] + values, minlength=rows * 256).reshape(rows, 256))
    if not histograms:
        return np.zeros((0, 256), dtype=np.int64)
    return histograms[0] if len(histograms) == 1 else np.concatenate(histograms)


def row_counts(data: bytes | memoryview, matrix: NDArray[np.int64] = BASES_MATRIX) -> NDArray[np.int64]:
    """The counts of each row of ``data`` (rows x 7: A, C, G, T, all bases, bases outside each alphabet)."""

    return byte_histograms(data) @ matrix


def alphabet_of(counts: NDArray[np.int64]) -> int:
    """What the bases of ``counts`` (one row or several) hold: ``NUCLEOTIDE``, ``PROTEIN`` or ``OTHER``."""

    total = counts.reshape(-1, _COLUMNS).sum(axis=0)
    return OTHER if total[BEYOND_PROTEIN] else PROTEIN if total[BEYOND_NUCLEOTIDE] else NUCLEOTIDE


# Multiplying a word of eight 0/1 bytes by this constant sums them into its top byte.
_BYTE_SUM = np.uint64(0x0101010101010101)
_TOP_BYTE = np.uint64(56)
_EIGHT = np.arange(8)
_LOWER = np.uint8(0x20)
_KNOWN_LETTERS: Final = ((A, b"a"), (C, b"c"), (G, b"g"), (T, b"t"), (T, b"u"))


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
    "ROW_BYTES",
    "SLICE_BYTES",
    "T",
    "alphabet_of",
    "byte_histograms",
    "count_matrix",
    "known_prefix_counts",
    "row_counts",
]
