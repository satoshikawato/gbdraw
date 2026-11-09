"""Record sequences and the files they are read from.

A record the render path reads keeps its sequence in its input file
(:mod:`gbdraw.io.record_sequences`). Records handed to a caller hold their
sequence instead, so the caller may rewrite the file. A file that changes under
records that still read it is refused by name, never read as other bases, and
the number of inputs is not bounded by the files a process may hold open.
"""

from __future__ import annotations

import gc
import gzip
import os
import shutil
import subprocess
import sys
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
    return isinstance(getattr(record.seq, "_data", None), RecordSequence)


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


def test_a_gzip_input_is_expanded_once_into_a_file_removed_with_its_records(tmp_path: Path) -> None:
    path = tmp_path / "a.gb.gz"
    path.write_bytes(gzip.compress(GENBANK.read_bytes()))
    records = read_flat_file_records(str(path))
    copies = [Path(table.source.path) for table in (record.seq._data._table for record in records)]
    assert copies[0].parent == Path(record_sequences._expansion_directory())
    assert copies[0].exists() and copies[0] != path
    assert str(records[0].seq) == _sequence(GENBANK)
    del records
    gc.collect()
    assert not copies[0].exists()


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
