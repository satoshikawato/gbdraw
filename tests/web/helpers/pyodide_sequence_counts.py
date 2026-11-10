"""Sequence counts in Pyodide (wasm32), where numpy's default integer and ``intp`` are 32 bits (OV-3356, R4-1).

Run by tests/web/pyodide-sequence-counts.test.mjs with the gbdraw sources on ``sys.path``. Each check raises
AssertionError with what differs.
"""

import gzip
import os
import tempfile

import numpy as np
from Bio.Seq import Seq

from gbdraw.core.base_counts import KNOWN, byte_counts, count_rows
from gbdraw.io import record_sequences
from gbdraw.io.genome import read_sequence_file_records
from gbdraw.io.record_sequences import _SequenceFile, _SequenceTable, _index_lines, known_base_counts_at

assert np.dtype(np.intp).itemsize == 4, "not a 32-bit platform: these checks would pass anywhere"
CHUNK = record_sequences._CHUNK_BYTES


def _oracle(bases: bytes, positions: list[int]) -> list[list[int]]:
    """A, C, G, T (U as T, either case) in ``bases[:p]`` for each position, counted by Python."""

    upper = bases.upper().replace(b"U", b"T")
    return [[upper[:p].count(letter) for p in positions] for letter in b"ACGT"]


def chunk_offsets_past_two_gib() -> None:
    """A run's chunks keep their file offsets exact at, across and past 2 GiB."""

    with tempfile.TemporaryDirectory() as folder:
        path = os.path.join(folder, "x.fa")
        with open(path, "wb") as handle:
            handle.write(b">x\nACGT\n")
        source = _SequenceFile(path)
        fasta = (b"ACGTNacgtn" * 7 + b"\n") * 2000
        genbank = b"".join(f"{i * 60 + 1:>9} ".encode() + b" ".join([b"acgtacgtac"] * 6) + b"\n" for i in range(2000))
        for raw, fmt in ((fasta, "fasta"), (genbank, "genbank")):
            for offset in (2**31 - 1000, 2**31 + 5, 5 * 2**30):
                table = _SequenceTable(source, key=source.table_key())
                _index_lines(table, offset, raw, fmt)
                ends = [start + size for start, size in zip(table.offsets, table.sizes)]
                expected = [offset] + ends[:-1]  # each chunk begins where the one before it ends
                assert list(table.offsets) == expected and ends[-1] == offset + len(raw), (
                    fmt, offset, list(table.offsets[:3]), expected[:3])
                assert max(table.sizes) <= CHUNK, (fmt, offset, max(table.sizes))


def counts_equal_a_copy() -> None:
    """Counts at every kind of position equal a whole copy's, for bases in memory and read from files."""

    rng = np.random.default_rng(5)
    letters = np.frombuffer(b"ACGTNacgtuU", dtype=np.uint8)
    bases = letters[rng.integers(0, len(letters), size=250_001)].tobytes()
    totals = count_rows(bases, 0).rows.sum(axis=0)
    assert totals.tolist() == byte_counts(bases)[KNOWN].tolist(), totals
    positions = [0, 1, 999, 1000, 1001, 65_535, 65_536, 123_457, 249_000, 250_000, 250_001]
    with tempfile.TemporaryDirectory() as folder:
        plain = os.path.join(folder, "x.fa")
        text = b">x\n" + b"\n".join(bases[i:i + 61] for i in range(0, len(bases), 61)) + b"\n"
        with open(plain, "wb") as handle:
            handle.write(text)
        packed = plain + ".gz"
        with gzip.open(packed, "wb") as handle:
            handle.write(text)
        sequences = [Seq(bases)] + [read_sequence_file_records(path, "fasta")[0].seq for path in (plain, packed)]
        for sequence in sequences:
            got = known_base_counts_at(sequence, np.array(positions, dtype=np.int64))
            assert got.tolist() == _oracle(bases, positions), (type(sequence), got.tolist())


chunk_offsets_past_two_gib()
counts_equal_a_copy()
