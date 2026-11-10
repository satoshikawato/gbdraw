"""Crops, reverse complements and windows read the record sequence through views, and read the same bases as copies.

:func:`gbdraw.io.record_sequences.sequence_view` and :func:`~gbdraw.io.record_sequences.two_turns`
make a crop, a reverse-complemented record and a crop across the origin read the bases
of the record they come from, and the GC counts and the crop digest read a sequence a
span at a time. Here each of them is compared with the whole copies gbdraw made before
(kept below as oracles): the bases, every feature's ``extract``, the crop digest, the
GC and skew windows and the definition's GC percent, for every GenBank, FASTA and GFF3
input under ``tests/`` (and gzip copies), and the reverse complement with Biopython's
over every byte value.
"""

from __future__ import annotations

import copy
import gzip
import hashlib
import pickle
import tracemalloc
from collections.abc import Iterator
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from Bio import SeqIO
from Bio.Seq import MutableSeq, Seq
from Bio.SeqFeature import SeqFeature, SimpleLocation
from Bio.SeqRecord import SeqRecord

from gbdraw.analysis import skew
from gbdraw.analysis.gc import calculate_gc_percent
from gbdraw.analysis.windows import cut_digest
from gbdraw.core import sequence_spans as sequence_spans_module
from gbdraw.core.record_metadata import record_content_digest
from gbdraw.exceptions import InputFileError
from gbdraw.io import record_sequences
from gbdraw.io.genome import parse_gff_fasta, read_flat_file_records, read_sequence_file_records
from gbdraw.io.record_select import reverse_records
from gbdraw.io.record_sequences import RecordSequence, same_bases, sequence_view, two_turns
from gbdraw.io.regions import RegionSpec, crop_record_to_region

TESTS = Path(__file__).parent
_FLAT = {".gb", ".gbk", ".gbff", ".genbank", ".embl", ".dbj"}
_FASTA = {".fa", ".fasta", ".fna", ".fas"}


# Oracles: the whole copies gbdraw made before (dev 18247abd3).


def _old_counted_sequence(seq: object) -> bytes:
    table = bytearray(b"N" * 256)
    for byte, counted in zip(b"ACGTUacgtu", b"ACGTTACGTT"):
        table[byte] = counted
    raw = bytes(seq) if isinstance(seq, (Seq, MutableSeq, bytes, bytearray)) else str(seq).encode("ascii")
    return raw.translate(bytes(table))


def _old_window_pieces(windows) -> list[bytes]:
    sequence = windows.source.sequence
    segments: list[bytes] = []
    for piece in windows.pieces:
        bases = sequence if piece.length == len(sequence) else sequence[piece.low : piece.low + piece.length]
        segments.append(_old_counted_sequence(bases if piece.direction > 0 else bases.reverse_complement()))
    return segments


def _old_gc_percent(seq: object) -> float:
    sequence = _old_counted_sequence(seq)
    known = len(sequence) - sequence.count(b"N")
    if not known:
        return 0.0
    return round(100.0 * (sequence.count(b"G") + sequence.count(b"C")) / known, 2)


def _old_digest(record: SeqRecord) -> str:
    digest = hashlib.blake2b(digest_size=16)
    digest.update(f"{record.id}\0{record.name}\0{record.description}\0".encode())
    digest.update(bytes(record.seq))
    digest.update(
        "".join(f"\0{feature.type}\0{feature.location}\0{feature.qualifiers!r}" for feature in record.features).encode()
    )
    annotations = getattr(record, "annotations", None) or {}
    for key in sorted(annotations):
        if not str(key).startswith("gbdraw_"):
            digest.update(f"\0{key}\0{annotations[key]}".encode())
    return digest.hexdigest()


def _old_skew_df(record: SeqRecord, window: int, step: int, *, circular: bool) -> pd.DataFrame:
    with pytest.MonkeyPatch.context() as patch:
        patch.setattr(skew, "_window_pieces", _old_window_pieces)
        return skew.skew_df(record, window, step, "GC", circular=circular)


# Inputs.


def _inputs(suffixes: set[str]) -> list[Path]:
    return sorted(
        path for path in TESTS.rglob("*")
        if path.is_file() and path.suffix.lower() in suffixes and "node_modules" not in path.parts
    )


def _gff3_pairs() -> list[tuple[Path, Path]]:
    pairs = []
    for gff in _inputs({".gff3", ".gff"}):
        partner = next((gff.with_suffix(s) for s in (".fasta", ".fa", ".fna") if gff.with_suffix(s).is_file()), None)
        if partner is None and "##FASTA" in gff.read_text(encoding="utf-8", errors="replace"):
            partner = gff
        if partner is not None:
            pairs.append((gff, partner))
    return pairs


_CASES = (
    [("flat", str(path.relative_to(TESTS))) for path in _inputs(_FLAT)]
    + [("fasta", str(path.relative_to(TESTS))) for path in _inputs(_FASTA)]
    + [("gff3", str(gff.relative_to(TESTS))) for gff, _ in _gff3_pairs()]
    + [("gzip", "test_inputs/NC_001416.gb"), ("gzip", "test_inputs/MjeNMV.gbk")]
)


def _records(kind: str, name: str, tmp_path: Path) -> list[SeqRecord]:
    path = TESTS / name
    if kind == "flat":
        records = read_flat_file_records(str(path))
    elif kind == "fasta":
        records = read_sequence_file_records(str(path), "fasta")
    elif kind == "gff3":
        fasta = dict(_gff3_pairs())[path]
        records = parse_gff_fasta(str(path), str(fasta))
    else:
        packed = tmp_path / (path.name + ".gz")
        packed.write_bytes(gzip.compress(path.read_bytes()))
        records = read_flat_file_records(str(packed))
    return [record for record in records if record.seq.defined and len(record.seq) >= 4]


def _spec(start: int, end: int, reverse: bool = False) -> RegionSpec:
    return RegionSpec(raw=f"{start}-{end}", file_selector=None, record_id=None, record_index=None,
                      start=start, end=end, reverse_complement=reverse)


def _crops(record: SeqRecord, whole: bytes) -> Iterator[tuple[str, SeqRecord, bytes]]:
    """Each kind of crop of ``record`` with the bases a copy gives it (an RNA record is complemented as RNA)."""

    rna = "RNA" in str(record.annotations.get("molecule_type", ""))

    def reverse_complement(bases: bytes) -> bytes:
        return bytes(Seq(bases).reverse_complement_rna() if rna else Seq(bases).reverse_complement())

    length = len(whole)
    span = min(length // 2, 300_000)
    start, end = length // 4 + 1, length // 4 + span
    crop = crop_record_to_region(record, _spec(start, end), circular=True)
    yield "crop", crop, whole[start - 1 : end]
    reversed_crop = crop_record_to_region(record, _spec(start, end, True), circular=True)
    yield "reverse-complement crop", reversed_crop, reverse_complement(whole[start - 1 : end])
    across_start, across_end = length - span // 2 + 1, max(1, span // 2)
    if across_end < across_start:
        across = crop_record_to_region(record, _spec(across_start, across_end), circular=True)
        yield "origin-crossing crop", across, (whole + whole)[across_start - 1 : length + across_end]
        across_rc = crop_record_to_region(record, _spec(across_start, across_end, True), circular=True)
        expected = reverse_complement((whole + whole)[across_start - 1 : length + across_end])
        yield "origin-crossing reverse-complement crop", across_rc, expected
    inner = len(crop.seq)
    sub_start, sub_end = inner // 3 + 1, max(inner // 3 + 1, 2 * inner // 3)
    yield "crop of a crop", crop_record_to_region(crop, _spec(sub_start, sub_end)), whole[start - 1 : end][
        sub_start - 1 : sub_end]
    sub = crop_record_to_region(reversed_crop, _spec(sub_start, sub_end, True))
    yield "reverse crop of a reverse crop", sub, reverse_complement(
        reverse_complement(whole[start - 1 : end])[sub_start - 1 : sub_end])


def _extracts(record: SeqRecord, seq: Seq) -> list[bytes | None]:
    found: list[bytes | None] = []
    for feature in record.features:
        try:
            found.append(bytes(feature.extract(seq)))
        except ValueError:
            found.append(None)
    return found


@pytest.mark.parametrize("case", _CASES, ids=lambda case: f"{case[0]}:{case[1]}")
def test_crops_read_the_bases_copies_read(case: tuple[str, str], tmp_path: Path) -> None:
    for record in _records(*case, tmp_path):
        whole = bytes(record.seq)
        length = len(whole)
        window, step = max(8, length // 64), max(4, length // 128)
        pd.testing.assert_frame_equal(
            skew.skew_df(record, window, step, "GC", circular=True), _old_skew_df(record, window, step, circular=True))
        assert calculate_gc_percent(record.seq) == _old_gc_percent(record.seq)
        for kind, crop, expected in _crops(record, whole):
            label = f"{record.id} {kind}"
            assert isinstance(crop.seq._data, RecordSequence), label  # read from the record, not copied
            assert bytes(crop.seq) == expected, label
            assert _extracts(crop, crop.seq) == _extracts(crop, Seq(expected)), label
            assert record_content_digest(crop) == _old_digest(crop), label
            assert cut_digest(crop) == _old_digest(crop), label
            assert calculate_gc_percent(crop.seq) == _old_gc_percent(expected), label
            if len(crop.seq):
                small = max(4, len(crop.seq) // 16)
                for circular in (True, False):
                    pd.testing.assert_frame_equal(
                        skew.skew_df(crop, small, max(2, small // 2), "GC", circular=circular),
                        _old_skew_df(crop, small, max(2, small // 2), circular=circular),
                        obj=label,
                    )


# Every byte value, in both cases, with IUPAC codes, U and gaps.
_ALPHABET = bytes(range(256)) + b"acgtuACGTUnNrRyYsSwWkKmMbBdDhHvV-.*" * 3


def _sequence_file(tmp_path: Path, bases: bytes) -> Seq:
    path = tmp_path / "alphabet.fasta"
    printable = bytes(byte for byte in bases if 33 <= byte < 127 and byte != ord(">"))
    path.write_bytes(b">alphabet\n" + b"\n".join(printable[i : i + 60] for i in range(0, len(printable), 60)) + b"\n")
    (record,) = read_sequence_file_records(str(path), "fasta")
    assert isinstance(record.seq._data, RecordSequence)
    assert bytes(record.seq) == printable
    return record.seq


@pytest.mark.parametrize("held", ["in memory", "from its file"])
def test_views_equal_biopython_slices_and_reverse_complements(held: str, tmp_path: Path) -> None:
    seq = Seq(_ALPHABET) if held == "in memory" else _sequence_file(tmp_path, _ALPHABET)
    whole = bytes(seq)
    length = len(whole)
    for start, stop in ((0, length), (3, 200), (17, 18), (40, 40), (-30, None), (5, 10_000), (90, 20)):
        plain = Seq(whole)[start:stop]
        assert bytes(sequence_view(seq, start, stop)) == bytes(plain)
        assert bytes(sequence_view(seq, start, stop, strand=-1)) == bytes(plain.reverse_complement())
        assert bytes(sequence_view(seq, start, stop, strand=-1, rna=True)) == bytes(plain.reverse_complement_rna())
    # A view of a view reads the first sequence's bases, whatever the turns, strands and complements on the way.
    reverse = sequence_view(seq, strand=-1)
    twice = sequence_view(reverse, strand=-1, rna=True)
    expected_twice = Seq(whole).reverse_complement().reverse_complement_rna()
    assert bytes(twice) == bytes(expected_twice)
    turns = two_turns(reverse)
    assert bytes(turns) == bytes(Seq(whole).reverse_complement()) * 2
    across = sequence_view(turns, length - 7, length + 9, strand=-1)
    assert bytes(across) == bytes((Seq(whole).reverse_complement() * 2)[length - 7 : length + 9].reverse_complement())
    for view in (reverse, twice, turns, across):
        data = view._data
        assert isinstance(data, RecordSequence)
        assert data._bases is (seq._data._bases if isinstance(seq._data, RecordSequence) else seq._data)
        expected = bytes(view)
        for key in (slice(2, 50), slice(None, None, -3), slice(-9, None), slice(7, 7), slice(1, 90, 4)):
            assert data[key] == expected[key]
        assert view[0] == chr(expected[0]) and view[-1] == chr(expected[-1])
        assert copy.copy(view)._data is data and copy.deepcopy(view)._data is data
        assert bytes(pickle.loads(pickle.dumps(view))) == expected
    assert sequence_view(seq) is seq


def test_sequences_a_view_cannot_read_take_biopython_s_own_operations() -> None:
    mutable = MutableSeq("ACGTUN")
    assert bytes(sequence_view(mutable, 1, 5, strand=-1)) == bytes(Seq("CGTU").reverse_complement())
    assert bytes(two_turns(mutable)) == b"ACGTUNACGTUN"
    undefined = Seq(None, 30)
    assert len(sequence_view(undefined, 5, 12, strand=-1)) == 7
    assert not sequence_view(undefined, 5, 12).defined
    assert len(two_turns(undefined)) == 60


def test_same_bases_reads_a_span_at_a_time(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(sequence_spans_module, "SPAN_BASES", 3)
    assert list(sequence_spans_module.sequence_spans(b"ACGTACG")) == [b"ACG", b"TAC", b"G"]
    assert same_bases(Seq("ACGTACGTA"), sequence_view(Seq("TACGTACGTAC"), 1, 10))
    assert not same_bases(Seq("ACGTACGTA"), Seq("ACGTACGTT"))
    assert not same_bases(Seq("ACGT"), Seq("ACGTA"))


@pytest.mark.parametrize("molecule_type", ["DNA", "mRNA", "ss-RNA", None])
def test_a_reversed_record_equals_biopython_s_reverse_complement(molecule_type: str | None) -> None:
    record = SeqRecord(Seq("AACGTUuacgnRYK-"), id="r1", name="r1", description="one", dbxrefs=["db:1"],
                       annotations={} if molecule_type is None else {"molecule_type": molecule_type},
                       letter_annotations={"phred_quality": list(range(15))})
    record.features = [SeqFeature(SimpleLocation(2, 9, strand=1), type="CDS")]
    (reversed_record,) = reverse_records((record,), True)
    want = record.reverse_complement(id=True, name=True, description=True, features=True, annotations=True,
                                     letter_annotations=True, dbxrefs=True)
    assert isinstance(reversed_record.seq._data, RecordSequence)
    assert bytes(reversed_record.seq) == bytes(want.seq)
    assert (reversed_record.id, reversed_record.name, reversed_record.description, reversed_record.dbxrefs) == (
        want.id, want.name, want.description, want.dbxrefs)
    assert reversed_record.letter_annotations == want.letter_annotations
    assert {key: value for key, value in reversed_record.annotations.items() if not key.startswith("gbdraw_")} == (
        want.annotations)


def test_a_protein_record_has_no_reverse_complement() -> None:
    record = SeqRecord(Seq("MKV"), id="p", annotations={"molecule_type": "protein"})
    with pytest.raises(ValueError, match="Proteins do not have complements!"):
        reverse_records((record,), True)


def test_a_view_of_a_file_that_changed_is_refused(tmp_path: Path) -> None:
    seq = _sequence_file(tmp_path, b"ACGT" * 5000)
    crop = sequence_view(seq, 100, 9000, strand=-1)
    record_sequences._CACHE.drop_file(seq._data._bases.source.key)
    (tmp_path / "alphabet.fasta").write_bytes(b">alphabet\n" + b"T" * 20000 + b"\n")
    with pytest.raises(InputFileError, match="alphabet.fasta changed"):
        bytes(crop)


def test_gc_windows_of_a_record_past_the_cache_hold_no_copy_of_it(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The GC and skew windows of a crop and of its reverse complement, and the GC percent, read a chunk at a time."""

    monkeypatch.setattr(record_sequences, "_CACHE_BYTES", 1 << 20)
    length = 24_000_000
    rng = np.random.default_rng(7)
    line = np.frombuffer(b"ACGTN", dtype=np.uint8)[rng.integers(0, 5, size=length, dtype=np.uint8)].tobytes()
    path = tmp_path / "long.fasta"
    with path.open("wb") as handle:
        handle.write(b">long\n")
        for start in range(0, length, 80):
            handle.write(line[start : start + 80] + b"\n")
    expected_gc = _old_gc_percent(line)
    del line
    (record,) = read_sequence_file_records(str(path), "fasta")
    record.annotations["topology"] = "circular"
    crop = crop_record_to_region(record, _spec(2_000_001, 20_000_000), circular=True)
    reversed_crop = crop_record_to_region(record, _spec(2_000_001, 20_000_000, True), circular=True)
    tracemalloc.start()
    try:
        windows = [skew.skew_df(drawn, 100_000, 50_000, "GC", circular=True) for drawn in (record, crop, reversed_crop)]
        gc_percent = calculate_gc_percent(record.seq)
        _, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    assert gc_percent == expected_gc
    pd.testing.assert_frame_equal(windows[2], _old_skew_df(reversed_crop, 100_000, 50_000, circular=True))
    assert peak < 12 << 20, peak  # the copies took 2 bytes per base of each: 36 MB and more for the crop alone


@pytest.mark.parametrize("chunk", [1, 3, 8, 13, 64])
def test_the_windows_of_a_reversed_crop_count_as_a_copy_at_any_chunk_size(
    chunk: int, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A reversed piece is read forward and reversed a chunk at a time; no chunk boundary changes a window."""

    rng = np.random.default_rng(chunk)
    bases = np.frombuffer(b"ACGTNacgt", dtype=np.uint8)[rng.integers(0, 9, size=5_000, dtype=np.uint8)].tobytes()
    (record,) = read_sequence_file_records(str(_fasta(tmp_path, bases)), "fasta")
    record.annotations["topology"] = "circular"
    monkeypatch.setattr(skew, "_COUNT_CHUNK", chunk)
    for spec in (_spec(1_001, 4_000, True), _spec(4_001, 1_000, True)):
        drawn = crop_record_to_region(record, spec, circular=True)
        for circular in (True, False):
            pd.testing.assert_frame_equal(skew.skew_df(drawn, 200, 70, "GC", circular=circular),
                                          _old_skew_df(drawn, 200, 70, circular=circular))


def _fasta(tmp_path: Path, bases: bytes) -> Path:
    path = tmp_path / "bases.fasta"
    path.write_bytes(b">bases\n" + b"\n".join(bases[i : i + 70] for i in range(0, len(bases), 70)) + b"\n")
    return path


def _old_crop_genbank(path: Path, out: Path, start: int, end: int) -> None:
    """``crop_genbank``'s run before: its own parse of the first record, and a copy of the crop's bases.

    The feature keys are given their INSDC spelling, as the shared reader gives them (D-G1b-5): the one change in
    what it writes from a file both read.
    """

    from gbdraw.crop_genbank import check_start_end_coords, crop_and_shift_features, write_cropped_genbank
    from gbdraw.io.input_text import LineStream, canonicalize_feature_types, read_sequence_text

    record = read_sequence_text(path, lambda lines: next(SeqIO.parse(LineStream(lines), "genbank")))
    canonicalize_feature_types(record.features)
    start_0, end_0 = check_start_end_coords(len(record.seq), start, end, str(record.id))
    cropped = SeqRecord(record.seq[start_0:end_0], id=record.id, name=record.name, description=record.description,
                        dbxrefs=record.dbxrefs.copy(), annotations=record.annotations)
    cropped.features = crop_and_shift_features(record.features, start_0, end_0)
    write_cropped_genbank(cropped, out, record.annotations["accessions"][0], start, end)


@pytest.mark.parametrize(
    "name", [str(path.relative_to(TESTS)) for path in _inputs({".gb", ".gbk", ".gbff"})], ids=lambda name: name)
def test_crop_genbank_writes_what_it_wrote_with_its_own_parse(
    name: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from gbdraw import crop_genbank

    path = TESTS / name
    try:
        (record, *_) = read_flat_file_records(str(path))
    except Exception as refusal:  # noqa: BLE001 - a file both refuse is not cropped either way
        try:
            _old_crop_genbank(path, tmp_path / "old.gb", 1, 10)
        except Exception:  # noqa: BLE001
            pytest.skip("neither reader reads it")
        pytest.fail(f"crop_genbank read this file with its own parse; the shared reader refuses it: {refusal}")
    if not record.seq.defined or not record.annotations.get("accessions"):
        pytest.skip("no sequence or accession to crop")
    length = len(record.seq)
    start, end = length // 3 + 1, max(length // 3 + 1, 2 * length // 3)
    out = tmp_path / "new.gb"
    monkeypatch.setattr("sys.argv", ["crop_genbank", "-i", str(path), "-o", str(out), "-s", str(start), "-e", str(end)])
    crop_genbank.main()
    try:
        _old_crop_genbank(path, tmp_path / "old.gb", start, end)
    except (AttributeError, TypeError, ValueError):
        # Its own parse refused a file the drawing reads (tabs in the layout, spaces in a LOCUS name); it now reads
        # the file as the drawing does.
        assert bytes(SeqIO.read(out, "genbank").seq) == bytes(record.seq)[start - 1 : end]
        return
    assert out.read_bytes() == (tmp_path / "old.gb").read_bytes()


def test_crop_genbank_writes_feature_keys_in_their_insdc_spelling(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A key written in another case (``cds``, ``Misc_Feature``) is written as the drawing reads it (D-G1b-5)."""

    from gbdraw import crop_genbank

    text = (TESTS / "test_inputs/NC_001416.gb").read_text()
    text = text.replace("     CDS             ", "     cds             ").replace(
        "     misc_feature    ", "     Misc_Feature    ")
    source = tmp_path / "lower.gb"
    source.write_text(text)
    out = tmp_path / "out.gb"
    monkeypatch.setattr("sys.argv", ["crop_genbank", "-i", str(source), "-o", str(out), "-s", "1", "-e", "48502"])
    crop_genbank.main()
    keys = {feature.type for feature in SeqIO.read(out, "genbank").features}
    assert {"CDS", "misc_feature"} <= keys and not {"cds", "Misc_Feature"} & keys


def test_crop_genbank_reads_an_embl_file(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """``crop_genbank`` reads the flat files the drawing reads, EMBL among them; its own parse read only GenBank."""

    from gbdraw import crop_genbank

    out = tmp_path / "out.gb"
    source = TESTS / "fixtures/input_variants/easyfig72_embl.embl"
    monkeypatch.setattr("sys.argv", ["crop_genbank", "-i", str(source), "-o", str(out), "-s", "101", "-e", "700"])
    crop_genbank.main()
    (record,) = read_flat_file_records(str(source))
    assert bytes(SeqIO.read(out, "genbank").seq) == bytes(record.seq)[100:700]


@pytest.mark.parametrize("name", ["fixtures/input_variants/cgview7_multi_record.gbff", "fixtures/b_dup_ids.gb"])
@pytest.mark.parametrize("packed", [False, True], ids=["plain", "gzip"])
def test_the_first_record_alone_reads_as_it_does_among_all(name: str, packed: bool, tmp_path: Path) -> None:
    """``crop_genbank`` reads the first record only (``first_only``), as its own parse did, not every record."""

    path = TESTS / name
    if packed:
        path = tmp_path / (path.name + ".gz")
        path.write_bytes(gzip.compress((TESTS / name).read_bytes()))
    every = read_flat_file_records(str(path))
    (first,) = read_flat_file_records(str(path), first_only=True)
    assert len(every) > 1
    assert (first.id, bytes(first.seq), [str(f.location) for f in first.features]) == (
        every[0].id, bytes(every[0].seq), [str(f.location) for f in every[0].features])
