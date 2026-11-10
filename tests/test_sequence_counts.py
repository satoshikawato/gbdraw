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
    """Count rows of 8 bases, so a short sequence has many row edges, and slices of 64 bytes, so a row and a line
    cross them."""

    rows = 8
    monkeypatch.setattr(base_counts, "ROW_BASES", rows)
    monkeypatch.setattr(record_sequences, "ROW_BASES", rows)
    monkeypatch.setattr(base_counts, "SLICE_BYTES", 64)
    monkeypatch.setattr(record_sequences, "SLICE_BYTES", 64)
    monkeypatch.setattr(record_sequences, "_CHUNK_BYTES", 64)
    yield rows


def _random_bases(seed: int, size: int, letters: bytes = b"ACGTNacgtuURY-") -> bytes:
    rng = np.random.default_rng(seed)
    return np.frombuffer(letters, dtype=np.uint8)[rng.integers(0, len(letters), size=size)].tobytes()


@pytest.mark.parametrize("held", ["fasta", "fasta.gz", "genbank", "Seq", "MutableSeq"])
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
        seq = {"Seq": Seq, "MutableSeq": MutableSeq}[held](bases)
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


_VARIANTS = Path(__file__).parent / "fixtures" / "input_variants"


def _irregular_genbank(path: Path, bases: bytes) -> Path:
    """Sequence lines in the usual layout, then some that are not (a digit among the letters, which Biopython keeps
    as a base: read line by line), then the usual again."""

    lines = []
    for number, start in enumerate(range(0, len(bases), 60)):
        chunk = bases[start : start + 60].decode()
        if 10 <= number < 20:
            lines.append(f"{start + 1:>9} " + " ".join(chunk[i : i + 7] + "7" for i in range(0, len(chunk), 7)))
        else:
            lines.append(f"{start + 1:>9} " + " ".join(chunk[i : i + 10] for i in range(0, len(chunk), 10)))
    length = sum(len(line[10:].replace(" ", "")) for line in lines)
    head = [f"LOCUS       {'odd':<16}{length:>12} bp    DNA     linear   UNK 01-JAN-1980", "ORIGIN"]
    path.write_text("\n".join([*head, *lines, "//", ""]))
    return path


def _layouts(tmp_path: Path) -> Iterator[tuple[str, str, str]]:
    bases = _random_bases(17, 2_345, b"ACGTNacgtuURY")
    upper = bases.upper()
    yield "irregular genbank", str(_irregular_genbank(tmp_path / "odd.gb", upper)), "genbank"
    widths = b"".join(bases[i : i + w] + (b"\r\n" if w % 2 else b"\n") for i, w in zip(range(0, 2_345, 50), [50] * 47))
    mixed = tmp_path / "mixed.fasta"
    mixed.write_bytes(b">mixed\n" + bases[:7] + b"\n" + bases[7:300] + b"\n" + widths[300:] + b"\n")
    yield "fasta of mixed widths and line ends", str(mixed), "fasta"
    wide = tmp_path / "wide.fasta"
    wide.write_bytes(b">wide\n" + bases[:1_000] + "\u00e9".encode() + bases[1_000:] + b"\n")
    yield "fasta with a letter of two bytes", str(wide), "fasta"
    for name, fmt in (("easyfig72_embl.embl", "genbank"), ("crlf_line_endings.gb", "genbank"),
                      ("lowercase_sequence_n_runs.gb", "genbank"), ("ggn117_crlf.fasta", "fasta"),
                      ("clinker67_single_line_fasta.fasta", "fasta"), ("cgview7_multi_record.gbff", "genbank")):
        yield name, str(_VARIANTS / name), fmt


def test_counts_at_every_position_equal_a_copy_in_every_layout(tmp_path: Path, small_rows: int) -> None:
    """Rows cross lines, chunks, kinds of chunk (a run, lines read line by line) and slices: every count equals the
    count of a whole copy, at every position, of every record."""

    for name, path, fmt in _layouts(tmp_path):
        records = read_sequence_file_records(path, fmt)
        assert records, name
        for record in records:
            bases = bytes(record.seq)
            positions = np.arange(0, len(bases) + 2, dtype=np.int64)
            assert known_base_counts_at(record.seq, positions).tolist() == _oracle(bases, positions), name
            assert known_base_counts_at(Seq(bases), positions).tolist() == _oracle(bases, positions), name


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
    # 10,000 is a row end, which reads nothing; the others read the bases of their row before them, at most a row.
    assert len(read) == 2 and max(read) <= base_counts.ROW_BASES + base_counts.ROW_BASES // 61 + 2  # and its line ends


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
    counts = base_counts.byte_counts(bytes(range(256)))
    assert counts[base_counts.KNOWN].tolist() == [2, 2, 2, 4]
    assert counts[base_counts.BASES] == 256
    assert counts[base_counts.BEYOND_NUCLEOTIDE] == 256 - len(base_counts.NUCLEOTIDE_CODES)
    assert counts[base_counts.BEYOND_PROTEIN] == 256 - len(base_counts.PROTEIN_CODES)
    every = np.frombuffer(bytes(range(256)) * 3, dtype=np.uint8)
    positions = np.arange(0, len(every) + 1, dtype=np.intp)
    assert base_counts.known_prefix_counts(every.tobytes(), positions).tolist() == _oracle(every.tobytes(), positions)
    # The count rows read the same rule (R4-2): their totals over every byte are the histogram's.
    rows = base_counts.count_rows(every.tobytes(), 0).rows
    assert rows.sum(axis=0).tolist() == base_counts.byte_counts(every.tobytes())[base_counts.KNOWN].tolist()


@pytest.mark.parametrize("compress", [False, True], ids=["plain", "gzip"])
def test_a_render_reads_no_whole_file_record_after_its_parse_pass(
    compress: bool, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The letter check, the definition's GC % and the GC/skew windows read the parse pass's counts, and the
    record's content digest is taken only when a reader compares it, which a render of a file does not: no read
    of the whole record is left. Window edges read their own rows, each once, so at most one more read of the
    bases however dense the windows, and of a gzip input at most one more expansion (dev read this record six
    times whole; a gzip input was also written out expanded)."""

    from gbdraw.io import input_text
    from tests.utils.cli_in_process import run_gbdraw

    bases = _random_bases(11, 1_200_000, b"ACGTNacgt")
    path = _fasta(tmp_path / "g.fasta", bases, width=80, compress=compress)
    whole: list[int] = []
    edges: list[int] = []
    expanded: list[int] = []
    read, row_bases = record_sequences.RecordSequence._read, record_sequences._SequenceTable.row_bases
    expand = input_text.GzipExpansion.read
    monkeypatch.setattr(record_sequences.RecordSequence, "_read",
                        lambda self, start, stop: whole.append(stop - start) or read(self, start, stop))
    monkeypatch.setattr(record_sequences._SequenceTable, "row_bases",
                        lambda self, *args: edges.append(len(data := row_bases(self, *args))) or data)
    monkeypatch.setattr(input_text.GzipExpansion, "read",
                        lambda self, size: expanded.append(len(data := expand(self, size))) or data)
    result = run_gbdraw(["circular", "--fasta", str(path), "-o", str(tmp_path / "out"), "-f", "svg"], cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    assert sum(whole) == 0
    assert sum(edges) <= len(bases)
    text_bytes = len(bases) + len(bases) // 80 + 8
    assert sum(expanded) <= (2 * text_bytes + record_sequences._BLOCK_BYTES if compress else 0)


# GD-13: rows of ``ROW_BASES`` bases, so the Auto windows of a record of 1 Mbp or more read no bases.


def test_the_auto_window_bounds_of_records_of_1_mbp_or_more_fall_on_row_ends() -> None:
    """Samples sit at multiples of the step and a window reaches half its width either side
    (``gbdraw.analysis.windows``), so every bound is a multiple of the row size when both are."""

    from gbdraw.auto_sizes import sliding_window
    from gbdraw.config.models import GbdrawConfig
    from gbdraw.config.toml import load_config_toml

    cfg = GbdrawConfig.from_dict(load_config_toml("gbdraw.data", "config.toml"))
    for length in (1_000_000, 1_234_567, 9_999_999, 10_000_000, 99_999_999, 100_000_000, 248_387_328, 3_000_000_000):
        window, step = sliding_window(length, cfg)
        assert step % base_counts.ROW_BASES == 0 and window // 2 % base_counts.ROW_BASES == 0, (length, window, step)


@pytest.mark.parametrize("topology", ["linear", "circular"])
def test_the_auto_windows_of_a_record_of_1_mbp_or_more_read_no_bases(
    topology: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The window counts are the parse pass's rows summed: no edge row is read, except round the origin of a
    circular record whose length is not a multiple of the row size (the windows there start ``length - window // 2``
    into it: two rows at most). The figure is the one dev drew (``tests/reference_outputs`` and the A/B)."""

    from tests.utils.cli_in_process import run_gbdraw

    bases = _random_bases(13, 1_234_567, b"ACGTNacgt")
    path = _fasta(tmp_path / "m.fasta", bases, width=80)
    edges: list[int] = []
    row_bases = record_sequences._SequenceTable.row_bases
    monkeypatch.setattr(record_sequences._SequenceTable, "row_bases",
                        lambda self, *args: edges.append(len(data := row_bases(self, *args))) or data)
    args = ["circular", "--fasta", str(path), "-o", str(tmp_path / "out"), "-f", "svg"]
    if topology == "circular":
        args += ["--record_topology", "circular"]
    result = run_gbdraw(args, cwd=tmp_path)
    assert result.returncode == 0, result.stderr
    assert sum(edges) == 0 if topology == "linear" else len(edges) <= 12 and sum(edges) <= 12 * base_counts.ROW_BASES


# The holders of count rows (review of G1b part 2: S-A, S-B, S-C).


def test_a_file_s_records_let_it_go_without_a_cyclic_collection_after_counts(tmp_path: Path) -> None:
    """Rows hold no reader of their bases, so no cycle keeps a table, its file, its handle or its cache blocks
    alive after the last record goes (a bound method in the rows did until a cyclic collection)."""

    import gc
    import weakref

    gc.collect()
    gc.disable()
    try:
        records = read_sequence_file_records(str(_fasta(tmp_path / "c.fasta", _random_bases(23, 50_000))), "fasta")
        table = records[0].seq._data._bases
        known_base_counts_at(records[0].seq, np.array([12_345, 50_000], dtype=np.int64))
        known_base_counts_at(sequence_view(records[0].seq, 101, 40_000, strand=-1), np.array([7, 999], dtype=np.int64))
        gone = [weakref.ref(table), weakref.ref(table.source)]
        del records, table
        assert [ref() for ref in gone] == [None, None]
        held = len(record_sequences._HELD)
        seq = Seq(_random_bases(29, 5_000))
        known_base_counts_at(sequence_view(seq, 3, 4_000), np.array([1_234], dtype=np.int64))
        assert len(record_sequences._HELD) == held + 1
        del seq
        assert len(record_sequences._HELD) == held
    finally:
        gc.enable()


def test_first_counts_from_many_threads_agree(tmp_path: Path) -> None:
    """A record's rows are built when its parse pass ends, and bases held in memory build theirs once under a
    lock, so threads asking for the first counts at once all get them (lazy rows raised IndexError)."""

    import sys
    import threading

    bases = _random_bases(31, 300_000)
    (record,) = read_sequence_file_records(str(_fasta(tmp_path / "t.fasta", bases)), "fasta")
    assert record.seq._data._bases._rows is not None  # built by the parse pass, before any reader
    expected = [row[0] for row in _oracle(bases, [len(bases)])]
    interval = sys.getswitchinterval()
    sys.setswitchinterval(1e-6)
    try:
        for seq in (record.seq, Seq(bases), Seq(bases)):
            results: list[list[int]] = []
            barrier = threading.Barrier(8)

            def count(seq: Seq = seq) -> None:
                barrier.wait()
                results.append(known_base_counts(sequence_view(seq, 0, len(bases))).tolist())

            threads = [threading.Thread(target=count) for _ in range(8)]
            for thread in threads:
                thread.start()
            for thread in threads:
                thread.join()
            assert results == [expected] * 8
    finally:
        sys.setswitchinterval(interval)


def test_bases_in_memory_asked_from_many_threads_are_counted_once(monkeypatch: pytest.MonkeyPatch) -> None:
    """One holder per bytes object, built under a lock, and its rows counted once under its own lock: threads that
    first count the same new bases (each through its own ``Seq``) count them once together. Slow holders and slow
    counts open the window each lock closes; without either lock the bases are counted several times."""

    import threading
    import time

    bases = _random_bases(33, 20_000)
    counted: list[int] = []
    count_rows = record_sequences.count_rows

    def slow_count_rows(data: bytes, *args: object, **kwargs: object) -> base_counts.CountedRows:
        counted.append(len(data))
        time.sleep(0.01)
        return count_rows(data, *args, **kwargs)  # type: ignore[arg-type]

    class SlowHeldBases(record_sequences._HeldBases):
        __slots__ = ()

        def __init__(self, *args: object, **kwargs: object) -> None:
            time.sleep(0.01)
            super().__init__(*args, **kwargs)  # type: ignore[arg-type]

    monkeypatch.setattr(record_sequences, "count_rows", slow_count_rows)
    monkeypatch.setattr(record_sequences, "_HeldBases", SlowHeldBases)
    seqs = [Seq(bases) for _ in range(8)]
    assert len({id(seq._data) for seq in seqs}) == 1  # one bytes object
    results: list[list[int]] = []
    barrier = threading.Barrier(len(seqs))

    def count(seq: Seq) -> None:
        barrier.wait()
        results.append(known_base_counts(seq).tolist())

    threads = [threading.Thread(target=count, args=(seq,)) for seq in seqs]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    assert results == [[row[0] for row in _oracle(bases, [len(bases)])]] * len(seqs)
    assert sum(counted) == len(bases)


def test_an_api_render_of_a_read_record_counts_its_bases_once(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A record the Python API reads holds its bases in memory and keeps the rows its parse pass counted: the
    render counts no base again (it counted the 2 Mbp eight more times, rebuilding rows per call)."""

    import gbdraw

    path = _fasta(tmp_path / "api.fasta", _random_bases(19, 2_000_000, b"ACGTNacgt"), width=80)
    counted = [0]
    real = record_sequences.count_rows

    def count_rows(data: bytes, *args: object) -> object:
        counted[0] += len(data)
        return real(data, *args)  # type: ignore[arg-type]

    monkeypatch.setattr(record_sequences, "count_rows", count_rows)
    records = gbdraw.read_fasta(str(path))
    parsed = counted[0]
    assert parsed <= path.stat().st_size  # the parse pass, once
    assert "<svg" in gbdraw.draw_circular(records).to_svg()
    assert counted[0] == parsed
