"""Record sequences and the files they are read from.

A record the render path reads keeps its sequence in its input file
(:mod:`gbdraw.io.record_sequences`). Records handed to a caller hold their
sequence instead, so the caller may rewrite the file. A file that changes under
records that still read it is refused by name, never read as other bases, and
the number of inputs is not bounded by the files a process may hold open.
"""

from __future__ import annotations

import gzip
import os
import random
import shutil
import subprocess
import sys
import tempfile
import textwrap
from pathlib import Path

import pytest
from Bio import SeqIO

import gbdraw
import gbdraw.api
from gbdraw.exceptions import InputFileError
from gbdraw.io import record_sequences
from gbdraw.io.genome import _read_fasta_records, read_flat_file_records
from gbdraw.io.record_sequences import RecordSequence

ROOT = Path(__file__).resolve().parents[1]
INPUTS = ROOT / "tests" / "test_inputs"
GENBANK = INPUTS / "NC_001416.gb"


def _reads_its_file(record) -> bool:
    data = getattr(record.seq, "_data", None)
    return isinstance(data, RecordSequence) and isinstance(data._bases, record_sequences._SequenceTable)


def _sequence(path: Path) -> str:
    return str(next(SeqIO.parse(path, "genbank")).seq)


def _fasta(tmp_path: Path) -> Path:
    path = tmp_path / "a.fasta"
    path.write_text(">a\nACGTACGTAC\nGTACGT\n>b\nTTTTGGGG\n", encoding="utf-8")
    return path


def test_records_given_to_callers_hold_their_sequences(tmp_path: Path) -> None:
    """The public readers detach every record from its file: writing them back over it keeps the bases."""

    path = tmp_path / "a.gb"
    shutil.copy(GENBANK, path)
    expected = _sequence(GENBANK)
    records = gbdraw.read_genbank(path)
    SeqIO.write(records, path, "genbank")
    assert _sequence(path) == expected

    assert not any(_reads_its_file(record) for record in gbdraw.api.load_gbks([str(GENBANK)]))
    assert not any(_reads_its_file(record) for record in gbdraw.read_fasta(_fasta(tmp_path)))
    gff = INPUTS / "NC_013668.gff3"
    assert not any(_reads_its_file(record) for record in gbdraw.read_gff(gff, INPUTS / "NC_013668.fasta"))


def test_the_render_path_reads_sequences_from_the_file() -> None:
    assert all(_reads_its_file(record) for record in read_flat_file_records(str(GENBANK)))


def test_records_read_through_a_link_keep_reading_after_the_link_goes(tmp_path: Path) -> None:
    """The Web Worker links each input of a request into that request's folder and removes the folder afterwards.

    Records Python built from the link (and keeps for the next request by the input's content token) read the
    file the link names, which stays while the token does.
    """

    cached = tmp_path / "cache" / "token.bin"
    cached.parent.mkdir()
    shutil.copy(GENBANK, cached)
    workspace = tmp_path / "request-1"
    workspace.mkdir()
    (workspace / "0001.bin").symlink_to(cached)
    records = read_flat_file_records(str(workspace / "0001.bin"))
    shutil.rmtree(workspace)
    assert str(records[0].seq) == _sequence(GENBANK)


def test_a_holder_of_records_asks_whether_their_files_still_hold_their_bytes(tmp_path: Path) -> None:
    """A file written again with the same bytes is a new file to the records: the holder parses it again."""

    path = tmp_path / "a.gb"
    shutil.copy(GENBANK, path)
    records = read_flat_file_records(str(path))
    assert record_sequences.sequence_files_unchanged(records)
    assert record_sequences.sequence_files_unchanged(gbdraw.api.load_gbks([str(path)]))  # in memory
    data = path.read_bytes()
    path.unlink()
    path.write_bytes(data)
    assert not record_sequences.sequence_files_unchanged(records)
    path.unlink()
    assert not record_sequences.sequence_files_unchanged(records)


@pytest.mark.parametrize("change", ["same size", "shorter", "removed"])
def test_a_file_that_changed_under_its_records_is_refused_by_name(change: str, tmp_path: Path) -> None:
    path = tmp_path / "a.gb"
    shutil.copy(GENBANK, path)
    record = read_flat_file_records(str(path))[0]
    data = path.read_bytes()
    stamp = path.stat().st_mtime_ns
    if change == "same size":  # other bases at the same places
        origin = data.index(b"ORIGIN")
        path.write_bytes(data[:origin] + data[origin:].replace(b"g", b"c"))
        os.utime(path, ns=(stamp + 1_000_000_000, stamp + 1_000_000_000))
    elif change == "shorter":
        path.write_bytes(data[: len(data) // 2])
    else:
        path.unlink()
    with pytest.raises(InputFileError, match="a.gb changed after gbdraw read it"):
        str(record.seq)


def test_more_inputs_than_open_files_allowed(tmp_path: Path) -> None:
    """Records of many inputs read their sequences under a low open-file limit (a file is opened per read)."""

    script = textwrap.dedent(
        f"""
        import os, resource, shutil, sys
        from Bio import SeqIO
        from gbdraw.io.genome import read_flat_file_records
        expected = str(next(SeqIO.parse({str(GENBANK)!r}, "genbank")).seq)
        paths = []
        for index in range(120):
            path = os.path.join({str(tmp_path)!r}, f"a{{index}}.gb")
            shutil.copy({str(GENBANK)!r}, path)
            paths.append(path)
        resource.setrlimit(resource.RLIMIT_NOFILE, (64, resource.getrlimit(resource.RLIMIT_NOFILE)[1]))
        records = [record for path in paths for record in read_flat_file_records(path)]
        assert all(str(record.seq) == expected for record in records)
        assert len(os.listdir("/proc/self/fd")) < 40, os.listdir("/proc/self/fd")
        print("ok", len(records))
        """
    )
    environment = dict(os.environ, PYTHONPATH=str(ROOT))
    result = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True, env=environment, timeout=300)
    assert result.returncode == 0, result.stderr[-2000:]
    assert result.stdout.split() == ["ok", "120"]


def test_records_released_while_the_reader_is_locked_do_not_wait_on_it(tmp_path: Path) -> None:
    """Records may be freed at any allocation, inside the reader's locked sections too (by a cyclic collection).

    Their release (closing the file, dropping its cached blocks) must not wait on
    those locks. The records of two inputs are released
    here while every lock of the reader is held, once by reference count and once
    by a collection run from inside a read. A subprocess with a timeout turns a
    hang into a failure.
    """

    script = textwrap.dedent(
        f"""
        import gc, gzip, os
        from gbdraw.io import record_sequences
        from gbdraw.io.genome import read_flat_file_records
        source = {str(GENBANK)!r}
        expected = None
        def records_of(name, compressed):
            path = os.path.join({str(tmp_path)!r}, name)
            data = open(source, "rb").read()
            open(path, "wb").write(gzip.compress(data) if compressed else data)
            records = read_flat_file_records(path)
            assert str(records[0].seq)  # its file is open, its blocks cached
            return records
        locks = [record_sequences._OPEN._lock, record_sequences._CACHE._lock]
        for compressed in (False, True):
            records = records_of(f"held{{compressed}}.gb", compressed)
            for lock in locks:
                lock.acquire()
            del records  # released by reference count, under every lock
            for lock in locks:
                lock.release()
        held = [records_of("cycle.gb.gz", True)]
        held.append(held)
        real_open = record_sequences.open_input_bytes
        def collecting_open(path):
            held.clear()
            gc.collect()  # the cycle goes inside the read that opens the next file
            return real_open(path)
        record_sequences.open_input_bytes = collecting_open
        records = read_flat_file_records(source)
        record_sequences._OPEN.read(records[0].seq._data._bases.source, 0, 10)
        assert str(records[0].seq) == str(read_flat_file_records(source)[0].seq)
        record_sequences._OPEN.read(records[0].seq._data._bases.source, 0, 10)  # releases what is queued
        assert not record_sequences._RELEASED
        print("ok")
        """
    )
    environment = dict(os.environ, PYTHONPATH=str(ROOT))
    try:
        result = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True, env=environment, timeout=60)
    except subprocess.TimeoutExpired:
        pytest.fail("releasing records waited on a lock the reader held")
    assert result.returncode == 0, result.stderr[-2000:]
    assert result.stdout.split() == ["ok"]


def test_a_gzip_input_is_read_where_it_is_and_nothing_of_it_is_written(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The records of a gzip input read the input itself, at offsets through points of its expansion; neither a
    pass nor a read nor a refusal writes a file (the Web would pay for one in the browser's storage)."""

    from gbdraw.io import input_text

    scratch = tmp_path / "scratch"
    scratch.mkdir()
    monkeypatch.setenv("TMPDIR", str(scratch))
    monkeypatch.setattr(tempfile, "tempdir", None)
    path, _ = _first_of_many(tmp_path)
    records = read_flat_file_records(str(path))
    assert {record.seq._data._bases.source.path for record in records} == {str(path)}
    assert [str(record.seq) for record in records] == [_sequence(GENBANK)] * 13
    monkeypatch.setattr(input_text, "MAX_EXPANDED_INPUT_BYTES", 1 << 20)
    with pytest.raises(InputFileError, match="expands to more than 1 MiB"):
        read_flat_file_records(str(path))
    assert list(scratch.iterdir()) == []


def test_reads_at_offsets_give_the_bases_of_the_plain_file_in_any_order(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A read starts from the point before it, or goes on from the last read; every base equals the plain file's."""

    monkeypatch.setattr(record_sequences, "_POINT_BYTES", 96 << 10)
    monkeypatch.setattr(record_sequences, "_BLOCK_BYTES", 64 << 10)
    monkeypatch.setattr(record_sequences, "_CACHE", record_sequences._BlockCache())
    path, _ = _first_of_many(tmp_path)
    (tmp_path / "many.gb").write_bytes(GENBANK.read_bytes() * 13)
    records = read_flat_file_records(str(path))
    plain = [str(record.seq) for record in read_flat_file_records(str(tmp_path / "many.gb"))]
    assert len(records[0].seq._data._bases.source.points) > 5
    chooser = random.Random(3)
    for _ in range(400):
        index = chooser.randrange(13)
        start = chooser.randrange(len(plain[index]))
        stop = min(len(plain[index]), start + chooser.choice((1, 900, 70_000)))
        assert str(records[index].seq[start:stop]) == plain[index][start:stop]


def test_reads_in_order_expand_each_byte_of_a_gzip_input_once(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Short reads in order (the edges of windows) go on from the last read: the input is expanded once more, not
    once from a point per read."""

    from gbdraw.io import input_text

    monkeypatch.setattr(record_sequences, "_CACHE", record_sequences._BlockCache())
    path, expanded = _first_of_many(tmp_path)
    records = read_flat_file_records(str(path))
    counted = [0]
    real_read = input_text.GzipExpansion.read

    def read(self, size):
        data = real_read(self, size)
        counted[0] += len(data)
        return data

    monkeypatch.setattr(input_text.GzipExpansion, "read", read)
    for record in records:
        table = record.seq._data._bases
        for row in range(0, len(table.rows().counts), 7):
            table.row_bases(row, record_sequences.ROW_BASES - 1)
    assert counted[0] <= expanded + record_sequences._BLOCK_BYTES


def test_the_parse_pass_reads_each_file_once_in_order(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """One sequential read names the format and indexes the file; no read at an offset."""

    opened: list[str] = []
    real_open = record_sequences.open_input_bytes

    def counted(path):
        opened.append(os.fspath(path))
        return real_open(path)

    def no_offset_read(*_args, **_kwargs):
        raise AssertionError("the parse pass read at an offset")

    monkeypatch.setattr(record_sequences, "open_input_bytes", counted)
    monkeypatch.setattr(record_sequences._OPEN, "read", no_offset_read)
    read_flat_file_records(str(GENBANK))
    _read_fasta_records(str(_fasta(tmp_path)))
    assert opened == [str(GENBANK), str(tmp_path / "a.fasta")]


def _first_of_many(tmp_path: Path) -> tuple[Path, int]:
    """A gzip file of lambda (177 kB of text) followed by twelve more records: 2.3 MB expanded."""

    text = GENBANK.read_bytes()
    path = tmp_path / "many.gb.gz"
    path.write_bytes(gzip.compress(text * 13))
    return path, len(text) * 13


def test_a_reader_that_stops_after_the_first_record_expands_only_what_it_read(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The parse pass expands a gzip input as it reads it, so a first-record read expands no more of the file."""

    monkeypatch.setattr(record_sequences, "_BLOCK_BYTES", 64 << 10)
    monkeypatch.setattr(record_sequences, "_POINT_BYTES", 64 << 10)
    path, expanded = _first_of_many(tmp_path)
    (first,) = read_flat_file_records(str(path), first_only=True)
    assert str(first.seq) == _sequence(GENBANK)
    assert first.seq._data._bases.source.points[-1].expanded < expanded // 4


@pytest.mark.parametrize("threads", [True, False], ids=["thread", "no-thread"])
def test_a_gzip_input_expands_ahead_of_its_parse_pass_and_stops_with_it(
    threads: bool, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The expansion runs on a thread of its own beside the pass (in turn where no thread can start, as in
    Pyodide); a pass that stops early or is refused stops it."""

    import threading

    from gbdraw.io import input_text

    if not threads:
        def no_thread(self):
            raise RuntimeError("can't start new thread")

        monkeypatch.setattr(threading.Thread, "start", no_thread)
    monkeypatch.setattr(record_sequences, "_BLOCK_BYTES", 64 << 10)
    path, _ = _first_of_many(tmp_path)
    (first,) = read_flat_file_records(str(path), first_only=True)
    assert [str(record.seq) for record in [first, *read_flat_file_records(str(path))]] == [_sequence(GENBANK)] * 14
    monkeypatch.setattr(input_text, "MAX_EXPANDED_INPUT_BYTES", 1 << 20)
    with pytest.raises(InputFileError, match="expands to more than 1 MiB"):
        read_flat_file_records(str(path))
    assert not [thread for thread in threading.enumerate() if thread.name == "gbdraw-gzip-ahead"]


def test_a_parse_that_fails_stops_the_expansion_ahead_of_it_while_its_error_is_held(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A parser error ends the pass at once: the expansion ahead of it stops and lets the file go, even while the
    error, and the frames its traceback holds, are kept."""

    import threading

    monkeypatch.setattr(record_sequences, "_BLOCK_BYTES", 64 << 10)
    path, _ = _first_of_many(tmp_path)

    def parse(lines, fmt):
        for _ in zip(range(100), lines):
            pass
        raise ValueError("a record the parser refuses")

    with pytest.raises(ValueError, match="parser refuses") as refused:
        record_sequences.read_flat_records(str(path), parse, lambda first: "genbank")
    assert refused.traceback
    assert not [thread for thread in threading.enumerate() if thread.name == "gbdraw-gzip-ahead"]


@pytest.mark.parametrize("slow", ["reader", "maker"])
def test_blocks_made_ahead_come_in_order_at_most_the_bound_ahead_with_an_error_after_them(slow: str) -> None:
    """The thread makes blocks while fewer than ``_AHEAD_BLOCKS`` wait, so a fast maker stops at the bound and a
    fast reader waits for each block; an error comes after every block made before it."""

    import threading
    import time

    made: list[int] = []
    taken: list[int] = []
    ahead: list[int] = []

    def blocks():
        for index in range(40):
            if slow == "maker":
                time.sleep(0.001)
            made.append(index)
            ahead.append(len(made) - len(taken))
            yield bytes([index])
        raise ValueError("refused after 40 blocks")

    with pytest.raises(ValueError, match="after 40 blocks"):
        for block in record_sequences._ahead(blocks()):
            if slow == "reader":
                time.sleep(0.001)
            taken.append(block[0])
    assert taken == list(range(40))
    # One block in the reader's hands, fewer than the bound waiting, and the one the maker is making.
    assert max(ahead) <= record_sequences._AHEAD_BLOCKS + 1
    assert not [thread for thread in threading.enumerate() if thread.name == "gbdraw-gzip-ahead"]


def test_the_cap_counts_what_a_reader_expanded_so_far(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """A file past the cap whose first record lies within it gives that record to a reader that stops there."""

    from gbdraw.io import input_text

    monkeypatch.setattr(record_sequences, "_BLOCK_BYTES", 64 << 10)
    monkeypatch.setattr(input_text, "MAX_EXPANDED_INPUT_BYTES", 1 << 20)
    path, _ = _first_of_many(tmp_path)
    (first,) = read_flat_file_records(str(path), first_only=True)
    assert str(first.seq) == _sequence(GENBANK)
    with pytest.raises(InputFileError, match="expands to more than 1 MiB"):
        read_flat_file_records(str(path))


def test_a_cap_refuses_only_the_bytes_past_it_at_any_block_size(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """The cap is byte-exact: a reader that stops within it reads, whatever the block size and the cap."""

    from gbdraw.io import input_text

    path, expanded = _first_of_many(tmp_path)
    first_record = len(GENBANK.read_bytes())
    for cap in (first_record + 22_801, first_record + 1):
        monkeypatch.setattr(input_text, "MAX_EXPANDED_INPUT_BYTES", cap)
        (first,) = read_flat_file_records(str(path), first_only=True)
        assert str(first.seq) == _sequence(GENBANK)
        with pytest.raises(InputFileError, match="expands to more than"):
            read_flat_file_records(str(path))
    monkeypatch.setattr(input_text, "MAX_EXPANDED_INPUT_BYTES", expanded)
    assert len(read_flat_file_records(str(path))) == 13  # a file of exactly the cap reads whole


def _many_with_a_cp1252_byte(tmp_path: Path) -> tuple[Path, bytes]:
    """Thirteen records, the sixth with a byte that is not UTF-8 (é in cp1252) in its definition."""

    text = GENBANK.read_bytes()
    odd = text.replace(b"DEFINITION  ", b"DEFINITION  caf\xe9 ", 1)
    plain = text * 5 + odd + text * 7
    path = tmp_path / "many.gb.gz"
    path.write_bytes(gzip.compress(plain))
    return path, plain


def test_a_second_pass_expands_the_input_again_and_its_records_read_it(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from gbdraw.io import input_text

    monkeypatch.setattr(record_sequences, "_BLOCK_BYTES", 64 << 10)
    monkeypatch.setattr(record_sequences, "_POINT_BYTES", 64 << 10)
    monkeypatch.setattr(input_text.locale, "getpreferredencoding", lambda _do_setlocale=True: "cp1252")
    path, plain = _many_with_a_cp1252_byte(tmp_path)
    (tmp_path / "plain.gb").write_bytes(plain)
    records = read_flat_file_records(str(path))
    expected = read_flat_file_records(str(tmp_path / "plain.gb"))
    assert [(r.id, r.description, str(r.seq)) for r in records] == [(r.id, r.description, str(r.seq)) for r in expected]
    assert "café" in records[5].description


def test_a_gzip_input_replaced_between_two_passes_is_refused(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    from gbdraw.io import input_text

    monkeypatch.setattr(record_sequences, "_BLOCK_BYTES", 64 << 10)
    path, plain = _many_with_a_cp1252_byte(tmp_path)

    def replaced(_do_setlocale: bool = True) -> str:  # asked between the UTF-8 pass and the second pass
        path.write_bytes(gzip.compress(plain.replace(b"ACCESSION   ", b"ACCESSION   X", 13)))
        return "cp1252"

    monkeypatch.setattr(input_text.locale, "getpreferredencoding", replaced)
    with pytest.raises(InputFileError, match="changed after gbdraw read it"):
        read_flat_file_records(str(path))
