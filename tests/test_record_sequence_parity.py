"""Record sequences read from their file equal the ``Seq`` Biopython builds from the whole text.

The readers of :mod:`gbdraw.io.genome` hand Biopython and BCBio-GFF the lines that are
not sequence and give each record a :class:`gbdraw.io.record_sequences.RecordSequence`.
Here every GenBank, EMBL, FASTA and GFF3 input under ``tests/`` (and gzip copies of
some) is read both ways: the sequence, every feature's extracted sequence, the
parser's warnings and its refusals must be the same.
"""

from __future__ import annotations

import copy
import gzip
import hashlib
import pickle
import warnings
from collections.abc import Callable, Iterator
from itertools import chain
from pathlib import Path

import pytest
from BCBio import GFF
from Bio import BiopythonParserWarning, SeqIO
from Bio.SeqRecord import SeqRecord

from gbdraw.exceptions import ParseError
from gbdraw.io.genome import parse_gff_fasta, read_flat_file_records, read_sequence_file_records
from gbdraw.io.input_text import (
    LineStream,
    flat_file_format,
    flat_file_stream,
    gff3_stream,
    read_sequence_text,
)
from gbdraw.io import record_sequences
from gbdraw.io.record_sequences import RecordSequence

TESTS = Path(__file__).parent
_FLAT = {".gb", ".gbk", ".gbff", ".genbank", ".embl", ".dbj"}
_FASTA = {".fa", ".fasta", ".fna", ".faa", ".fas"}


def _inputs(suffixes: set[str]) -> list[Path]:
    return sorted(
        path for path in TESTS.rglob("*")
        if path.is_file() and path.suffix.lower() in suffixes and "node_modules" not in path.parts
    )


def _biopython_flat(path: Path) -> list[SeqRecord]:
    """Today's reading: the whole text, sequence lines included, through Biopython."""

    def parse(lines: Iterator[str]) -> list[SeqRecord]:
        head: list[str] = []
        for line in lines:
            head.append(line)
            if line.strip():
                break
        fmt = flat_file_format(head[-1] if head else "")
        return list(SeqIO.parse(flat_file_stream(chain(head, lines), fmt), fmt))

    return read_sequence_text(path, parse)


def _biopython_fasta(path: Path) -> list[SeqRecord]:
    def parse(lines: Iterator[str]) -> list[SeqRecord]:
        for line in lines:
            if line.startswith(">"):
                return list(SeqIO.parse(LineStream(chain([line], lines)), "fasta"))
        return []

    return read_sequence_text(path, parse)


def _digest(record: SeqRecord) -> tuple[int, str | None]:
    if not record.seq.defined:
        return len(record.seq), None
    return len(record.seq), hashlib.sha256(bytes(record.seq)).hexdigest()


def _feature_digests(record: SeqRecord) -> list[str | None]:
    if not record.seq.defined:
        return []
    digests: list[str | None] = []
    for feature in record.features:
        try:
            digests.append(hashlib.sha256(bytes(feature.extract(record.seq))).hexdigest())
        except ValueError:  # a location the record cannot hold is refused both ways
            digests.append(None)
    return digests


def _assert_same(expected: list[SeqRecord], actual: list[SeqRecord]) -> None:
    assert [record.id for record in actual] == [record.id for record in expected]
    for want, got in zip(expected, actual):
        assert _digest(got) == _digest(want), want.id
        assert _feature_digests(got) == _feature_digests(want), want.id


def _warned(read: Callable[[], object]) -> tuple[object, set[str]]:
    """What ``read`` returns or the parser error that stops it, and the parser warnings it gives (each once)."""

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", BiopythonParserWarning)
        try:
            result: object = read()
        except ParseError as error:  # gbdraw's reader names the parser's error as its cause
            cause = error.__cause__
            result = ("refused", type(cause).__name__, str(cause))
        except ValueError as error:
            result = ("refused", type(error).__name__, str(error))
    return result, {str(item.message) for item in caught if issubclass(item.category, BiopythonParserWarning)}


@pytest.mark.parametrize("path", _inputs(_FLAT), ids=lambda path: str(path.relative_to(TESTS)))
def test_flat_file_sequences_equal_biopython(path: Path) -> None:
    _assert_same(_biopython_flat(path), read_flat_file_records(str(path)))


@pytest.mark.parametrize("path", _inputs(_FASTA), ids=lambda path: str(path.relative_to(TESTS)))
def test_fasta_sequences_equal_biopython(path: Path) -> None:
    _assert_same(_biopython_fasta(path), read_sequence_file_records(str(path), "fasta"))


def _gff3_pairs() -> list[tuple[Path, Path]]:
    pairs = []
    for gff in _inputs({".gff3", ".gff"}):
        partner = next((gff.with_suffix(suffix) for suffix in (".fasta", ".fa", ".fna") if gff.with_suffix(suffix).is_file()), None)
        if partner is None and "##FASTA" in gff.read_text(encoding="utf-8", errors="replace"):
            partner = gff
        if partner is not None:
            pairs.append((gff, partner))
    return pairs


@pytest.mark.parametrize("pair", _gff3_pairs(), ids=lambda pair: str(pair[0].relative_to(TESTS)))
def test_gff3_records_equal_bcbio(pair: tuple[Path, Path]) -> None:
    gff, fasta = pair
    whole = read_sequence_text(gff, lambda lines: list(GFF.parse(gff3_stream(lines))))
    by_id = {record.id: record for record in _biopython_fasta(fasta)}
    for record in whole:
        record.seq = by_id[record.id].seq
    actual = parse_gff_fasta(str(gff), str(fasta))
    by_actual = {record.id: record for record in actual}
    assert set(by_actual) == {record.id for record in whole}
    for want in whole:
        got = by_actual[want.id]
        assert _digest(got) == _digest(want)
        assert [feature.id for feature in got.features] == [feature.id for feature in want.features]


@pytest.mark.parametrize("name", ["test_inputs/NC_001416.gb", "test_inputs/MjeNMV.gbk"])
def test_gzip_flat_file_expands_once_and_reads_the_same(name: str, tmp_path: Path) -> None:
    source = TESTS / name
    packed = tmp_path / (source.name + ".gz")
    packed.write_bytes(gzip.compress(source.read_bytes()))
    _assert_same(_biopython_flat(source), read_flat_file_records(str(packed)))


def test_gzip_fasta_reads_the_same(tmp_path: Path) -> None:
    source = next(iter(_inputs(_FASTA)))
    packed = tmp_path / "entries.fasta.gz"
    packed.write_bytes(gzip.compress(source.read_bytes()))
    _assert_same(_biopython_fasta(source), read_sequence_file_records(str(packed), "fasta"))


_LOCUS = "LOCUS       SYN                     {length} bp    DNA     linear   UNK 01-JAN-2020\n"
_HEAD = "DEFINITION  synthetic.\nACCESSION   SYN\nVERSION     SYN.1\nFEATURES             Location/Qualifiers\n" \
        "     source          1..{length}\n     CDS             3..14\n"


def _genbank(sequence_lines: str, length: int = 24, end: str = "//\n") -> str:
    return _LOCUS.format(length=length) + _HEAD.format(length=length) + "ORIGIN\n" + sequence_lines + end


_SEQUENCE_CASES = {
    "standard": "        1 acgtacgtac gtacgtacgt acgt\n",
    "upper case": "        1 ACGTACGTAC GTACGTACGT ACGT\n",
    "trailing spaces": "        1 acgtacgtac gtacgtacgt acgt   \n",
    "blank line": "        1 acgtacgtac\n\n       11 gtacgtacgt acgt\n",
    "one space too many": "         1 acgtacgtac gtacgtacgt acgt\n",
    "ten-digit coordinate": "1000000001 acgtacgtac gtacgtacgt acgt\n",
    "digit in the sequence": "        1 acgtacg1ac gtacgtacgt acgt\n",
    "tab in the sequence": "        1 acgtac\tgtac gtacgtacgt acgt\n",
    "no coordinates": "acgtacgtacgtacgtacgtacgt\n",
    "mal-formed": "1 acgtacgtacgtacgtacgtacgt\n",
    "star and dash": "        1 acgt*cgtac gtacg-acgt acgt\n",
    "many lines": "".join(f"{1 + 60 * i:>9} " + " ".join(["acgtacgtac"] * 6) + "\n" for i in range(3000)),
    "length differs": "        1 acgtacgtac gtacgt\n",
    "control character": "        1 acgtacgtac gtacgtacgt acgt\x1c\n",
}


@pytest.mark.parametrize("ending", ["\n", "\r\n", "\r"], ids=["lf", "crlf", "cr"])
@pytest.mark.parametrize("case", sorted(_SEQUENCE_CASES))
def test_genbank_sequence_lines_read_as_biopython_reads_them(case: str, ending: str, tmp_path: Path) -> None:
    lines = _SEQUENCE_CASES[case]
    length = 3000 * 60 if case == "many lines" else 24
    path = tmp_path / "synthetic.gb"
    path.write_bytes(_genbank(lines, length).replace("\n", ending).encode("utf-8"))
    expected, expected_warnings = _warned(lambda: _biopython_flat(path))
    actual, actual_warnings = _warned(lambda: read_flat_file_records(str(path)))
    assert actual_warnings == expected_warnings
    if isinstance(expected, tuple):
        assert actual == expected
    else:
        assert isinstance(actual, list) and isinstance(expected, list)
        _assert_same(expected, actual)


def test_genbank_sequence_without_its_end_line(tmp_path: Path) -> None:
    path = tmp_path / "cut.gb"
    path.write_text(_genbank("        1 acgtacgtac gtacgtacgt acgt\n", end=""), encoding="utf-8")
    expected, expected_warnings = _warned(lambda: _biopython_flat(path))
    actual, actual_warnings = _warned(lambda: read_flat_file_records(str(path)))
    assert actual_warnings == expected_warnings
    assert isinstance(actual, list) and isinstance(expected, list)
    _assert_same(expected, actual)


_FASTA_CASES = {
    "plain": ">a one\nACGTACGT\nacgt\n>b\nNNNN\n",
    "spaces and tabs": ">a\nAC GT\tAC\n  GT \n",
    "text before the first entry": "comment\n\n>a\nACGT\n",
    "empty entry": ">a\n>b\nACGT\n",
    "one line, no newline": ">a\nACGTACGTACGT",
    "not ascii": ">a\nACGTé\n",
}


def test_a_contig_line_after_the_sequence_is_refused(tmp_path: Path) -> None:
    """Biopython refuses a record whose sequence lines a CONTIG line follows; the reader refuses its file."""

    path = tmp_path / "contig.gb"
    path.write_text(
        "LOCUS       X1                        20 bp    DNA     linear   BCT 01-JAN-2000\n"
        "DEFINITION  x.\nACCESSION   X1\nVERSION     X1.1\n"
        "FEATURES             Location/Qualifiers\n     source          1..20\n"
        "ORIGIN\n        1 acgtacgtac gtacgtacgt\nCONTIG      join(A1.1:1..20)\n//\n",
        encoding="utf-8",
    )
    with pytest.raises(AssertionError):
        list(SeqIO.parse(path, "genbank"))
    with pytest.raises(ParseError) as refused:
        read_flat_file_records(str(path))
    assert "CONTIG" in str(refused.value.__cause__)


@pytest.mark.parametrize("ending", ["\n", "\r\n", "\r"], ids=["lf", "crlf", "cr"])
@pytest.mark.parametrize("case", sorted(_FASTA_CASES))
def test_fasta_lines_read_as_biopython_reads_them(case: str, ending: str, tmp_path: Path) -> None:
    path = tmp_path / "entries.fasta"
    path.write_bytes(_FASTA_CASES[case].replace("\n", ending).encode("utf-8"))
    _assert_same(_biopython_fasta(path), read_sequence_file_records(str(path), "fasta"))


def test_embl_sequence_reads_as_biopython_reads_it() -> None:
    path = next(iter(_inputs({".embl"})))
    _assert_same(_biopython_flat(path), read_flat_file_records(str(path)))


@pytest.mark.parametrize("bulk", [None, 1 << 12], ids=["chunks", "bulk"])
def test_record_sequence_reads_slices_and_copies(bulk: int | None, tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A short read reads its chunks through the cache, a long one is decoded in bulk around it."""

    if bulk is not None:
        monkeypatch.setattr(record_sequences, "_BLOCK_BYTES", bulk)
    lines = _SEQUENCE_CASES["many lines"]
    path = tmp_path / "long.gb"
    path.write_text(_genbank(lines, 3000 * 60), encoding="utf-8")
    (record,) = read_flat_file_records(str(path))
    (want,) = _biopython_flat(path)
    data = record.seq._data
    assert isinstance(data, RecordSequence)
    whole = bytes(want.seq)
    for key in (slice(5, 70000), slice(None, None, -1), slice(179990, 10, -7), slice(3, 100, 3), slice(-5, None),
                slice(1000, 41000), slice(0, 0), slice(0, None)):
        assert data[key] == whole[key]
    assert record.seq[-1] == want.seq[-1] and record.seq[0] == want.seq[0]
    assert str(record.seq.reverse_complement()) == str(want.seq.reverse_complement())
    assert copy.deepcopy(record.seq)._data is data and copy.copy(record.seq)._data is data
    assert bytes(pickle.loads(pickle.dumps(record.seq))) == whole
    with pytest.raises(IndexError):
        data[len(whole)]
