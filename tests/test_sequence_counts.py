"""One producer of base counts: the parse pass counts a record's bases, and every count of a span, a window or the
letter check reads those counts (``gbdraw.io.record_sequences``), never a whole copy of the sequence.

The oracle is the count of a whole copy, as gbdraw counted before: A, C, G and T in either case, U as T.
"""

from __future__ import annotations

import gzip
import tracemalloc
from collections.abc import Iterator
from pathlib import Path

import numpy as np
import pytest
from Bio.Seq import MutableSeq, Seq

from gbdraw.analysis import skew
from gbdraw.core import base_counts
from gbdraw.exceptions import ValidationError
from gbdraw.io import record_sequences
from gbdraw.io.genome import load_query_fasta, read_sequence_file_records, require_nucleotide_records
from gbdraw.io.record_sequences import known_base_counts, known_base_counts_at, sequence_view, two_turns


def _oracle(seq: object, positions: object) -> list[list[int]]:
    """The counts of A, C, G and T before each position, from a whole copy (positions clipped to it)."""

    bases = bytes(seq).upper().replace(b"U", b"T")  # type: ignore[call-overload]
    clipped = [max(0, min(int(p), len(bases))) for p in positions]  # type: ignore[attr-defined]
    return [[bases[:p].count(letter) for p in clipped] for letter in (b"A", b"C", b"G", b"T")]


def _fasta(path: Path, bases: bytes, width: int = 61, *, compress: bool = False) -> Path:
    text = b">bases\n" + b"\n".join(bases[i : i + width] for i in range(0, len(bases), width)) + b"\n"
    path.write_bytes(gzip.compress(text) if compress else text)
    return path


def _genbank(path: Path, bases: bytes) -> Path:
    lines = [f"LOCUS       {'bases':<16}{len(bases):>12} bp    DNA     linear   UNK 01-JAN-1980", "ORIGIN"]
    for start in range(0, len(bases), 60):
        chunk = bases[start : start + 60].decode()
        lines.append(f"{start + 1:>9} " + " ".join(chunk[i : i + 10] for i in range(0, len(chunk), 10)))
    path.write_text("\n".join([*lines, "//", ""]))
    return path


@pytest.fixture
def small_rows(monkeypatch: pytest.MonkeyPatch) -> Iterator[int]:
    """Count rows of 8 raw bytes, so a short sequence has many row edges."""

    rows = 8
    monkeypatch.setattr(base_counts, "ROW_BYTES", rows)
    monkeypatch.setattr(base_counts, "_ROW_KEYS", np.arange(base_counts.SLICE_BYTES, dtype=np.uint32) // rows * 256)
    monkeypatch.setattr(record_sequences, "ROW_BYTES", rows)
    yield rows


def _random_bases(seed: int, size: int, letters: bytes = b"ACGTNacgtuURY-") -> bytes:
    rng = np.random.default_rng(seed)
    return np.frombuffer(letters, dtype=np.uint8)[rng.integers(0, len(letters), size=size)].tobytes()


@pytest.mark.parametrize("held", ["fasta", "fasta.gz", "genbank", "Seq", "MutableSeq", "bytes"])
def test_counts_of_any_span_equal_a_whole_copy_s_count(held: str, tmp_path: Path, small_rows: int) -> None:
    bases = _random_bases(7, 3_000, b"ACGTNacgtuURY") if held == "genbank" else _random_bases(7, 3_000)
    if held.startswith("fasta"):
        (record,) = read_sequence_file_records(str(_fasta(tmp_path / f"b.{held}", bases, compress=held.endswith("gz"))), "fasta")
        seq: object = record.seq
        bases = bytes(record.seq)
    elif held == "genbank":
        (record,) = read_sequence_file_records(str(_genbank(tmp_path / "b.gb", bases)), "genbank")
        seq = record.seq
        bases = bytes(record.seq)
    else:
        seq = {"Seq": Seq, "MutableSeq": MutableSeq, "bytes": bytes}[held](bases)
    positions = np.array([40, 0, 2_999, 7, 7, 8, 9, 16, 3_000, 3_200, -5, 1_234], dtype=np.int64)
    assert known_base_counts_at(seq, positions).tolist() == _oracle(bases, positions)
    assert known_base_counts(seq).tolist() == [row[0] for row in _oracle(bases, [len(bases)])]
    if not isinstance(seq, Seq):
        return
    for start, stop in ((0, 3_000), (5, 13), (8, 16), (1_001, 2_500), (17, 17)):
        for strand in (1, -1):
            view = sequence_view(seq, start, stop, strand=strand)  # a crop, reversed or not
            inside = np.array([0, 1, len(view) // 2, max(0, len(view) - 1), len(view)], dtype=np.int64)
            assert known_base_counts_at(view, inside).tolist() == _oracle(bytes(view), inside), (start, stop, strand)
    across = two_turns(sequence_view(seq, 2_000, 3_000, strand=-1))  # reversed and across the origin
    inside = np.arange(0, len(across) + 1, 37, dtype=np.int64)
    assert known_base_counts_at(across, inside).tolist() == _oracle(bytes(across), inside)


def test_window_counts_at_row_edges_equal_a_copy_at_any_row_size(tmp_path: Path, small_rows: int) -> None:
    bases = _random_bases(3, 5_000)
    (record,) = read_sequence_file_records(str(_fasta(tmp_path / "w.fasta", bases, width=13)), "fasta")
    segments = [sequence_view(record.seq, 4_000, 5_000, strand=-1), Seq(b""), record.seq[:7], sequence_view(record.seq, 7, 4_000)]
    joined = b"".join(bytes(segment) for segment in segments)
    positions = np.array([0, 6, 7, 8, 999, 1_000, 1_001, 1_007, 1_008, 4_999, 5_000, 5_009, -1], dtype=np.int64)
    assert skew._known_counts_at(segments, positions).tolist() == _oracle(joined, positions)


def test_counts_of_bases_in_memory_need_no_array_as_long_as_the_sequence() -> None:
    bases = _random_bases(1, 8_000_000, b"ACGT")
    positions = np.arange(0, len(bases), 50_000, dtype=np.int64)
    tracemalloc.start()
    counts = known_base_counts_at(Seq(bases), positions)
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    assert peak < 6 << 20
    assert counts[2, 3] == bases[:150_000].count(b"G")


def test_a_record_read_from_a_file_is_counted_by_its_parse_pass_and_not_read_again(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    bases = _random_bases(5, 300_000, b"ACGTN")
    (record,) = read_sequence_file_records(str(_fasta(tmp_path / "r.fasta", bases)), "fasta")
    read: list[int] = []
    real = record_sequences._SequenceFile.read_at
    monkeypatch.setattr(record_sequences._SequenceFile, "read_at",
                        lambda self, offset, size: read.append(size) or real(self, offset, size))
    assert known_base_counts(record.seq).tolist() == [row[0] for row in _oracle(bases, [len(bases)])]
    assert read == []  # the whole record's counts are its rows' counts
    positions = np.array([10_000, 150_001, 299_999], dtype=np.int64)
    assert known_base_counts_at(record.seq, positions).tolist() == _oracle(bases, positions)
    assert len(read) == 3 and max(read) <= base_counts.ROW_BYTES  # one edge row each


# The letter check: the same refusals, with the same letter, from the counts of the parse pass.


@pytest.mark.parametrize(
    ("text", "letter"),
    [
        (b"acgtnACGTN" * 500 + b"U", None),
        (b"ACGTRYSWKMBDHVN" * 600 + b"acgturyswkmbdhvn", None),
        (b"ACGT" * 3_000 + b"-ACGT", "-"),
        (b"ACGT" * 3_000 + b"x", "x"),
        (b"ACGT" * 3_000 + b"E", "E"),
        (b"ACGT" * 3_000 + b"1", "1"),
        (b"ACGT" * 3_000 + b"\xc3\xa9", "\xe9"),  # a letter of two bytes is one letter (was an uncaught decode error)
        (b"*ACGT", "*"),
    ],
)
def test_the_letter_check_refuses_the_first_letter_outside_the_nucleotide_codes(
    text: bytes, letter: str | None, tmp_path: Path, small_rows: int
) -> None:
    path = _fasta(tmp_path / "q.fasta", text, width=70)
    records = read_sequence_file_records(str(path), "fasta")
    if letter is None:
        assert require_nucleotide_records(records, str(path)) == records
        assert load_query_fasta(str(path))[0][1] == "nucleotide"
        return
    with pytest.raises(ValidationError) as refused:
        require_nucleotide_records(records, str(path))
    assert refused.value.diagnostic["reason"] == "NOT_NUCLEOTIDE"
    assert refused.value.diagnostic["codepoint"] == ord(letter)
    assert f"holds {letter!r}, which is not a nucleotide code" in str(refused.value)
    if letter.isascii() and (letter.isalpha() or letter == "*"):
        assert load_query_fasta(str(path))[0][1] == "protein"
    else:
        with pytest.raises(ValidationError) as query:
            load_query_fasta(str(path))
        assert query.value.diagnostic["reason"] == "NOT_SEQUENCE"
        assert query.value.diagnostic["codepoint"] == ord(letter)


def test_the_letter_check_of_bases_in_memory_reads_the_same_rule() -> None:
    from gbdraw.io.genome import _first_letter_outside

    assert _first_letter_outside(Seq("ACGTN" * 5_000 + "acgu"), base_counts.NUCLEOTIDE) is None
    assert _first_letter_outside(Seq("ACGTN" * 5_000 + "QZ"), base_counts.NUCLEOTIDE) == "Q"
    assert _first_letter_outside(Seq("ACGTN" * 5_000 + "QZ-"), base_counts.PROTEIN) == "-"
    assert _first_letter_outside(MutableSeq("ACGTN" * 5_000 + "QZ-"), base_counts.PROTEIN) == "-"


def test_every_byte_counts_as_the_counting_rule_says() -> None:
    counts = base_counts.row_counts(bytes(range(256)))[0]
    assert counts[base_counts.KNOWN].tolist() == [2, 2, 2, 4]
    assert counts[base_counts.BASES] == 256
    assert counts[base_counts.BEYOND_NUCLEOTIDE] == 256 - len(base_counts.NUCLEOTIDE_CODES)
    assert counts[base_counts.BEYOND_PROTEIN] == 256 - len(base_counts.PROTEIN_CODES)
    every = np.frombuffer(bytes(range(256)) * 3, dtype=np.uint8)
    positions = np.arange(0, len(every) + 1, dtype=np.intp)
    assert base_counts.known_prefix_counts(every.tobytes(), positions).tolist() == _oracle(every.tobytes(), positions)


def test_a_render_reads_a_file_record_s_bases_once_after_its_parse_pass(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The letter check, the definition's GC % and the GC/skew windows read the parse pass's counts; the one whole
    read left is the record's content digest. Window edges read their own rows, each once, so at most one more read
    of the bases however dense the windows (dev read this record six times whole)."""

    from tests.utils.cli_in_process import run_gbdraw

    bases = _random_bases(11, 1_200_000, b"ACGTNacgt")
    path = _fasta(tmp_path / "g.fasta", bases, width=80)
    whole: list[int] = []
    edges: list[int] = []
    read, row_bases = record_sequences.RecordSequence._read, record_sequences._SequenceTable._row_bases
    monkeypatch.setattr(record_sequences.RecordSequence, "_read",
                        lambda self, start, stop: whole.append(stop - start) or read(self, start, stop))
    monkeypatch.setattr(record_sequences._SequenceTable, "_row_bases",
                        lambda self, row: edges.append(len(data := row_bases(self, row))) or data)
    result = run_gbdraw(["circular", "--fasta", str(path), "-o", str(tmp_path / "out"), "-f", "svg"], cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    assert sum(whole) == len(bases)
    assert sum(edges) <= len(bases)
