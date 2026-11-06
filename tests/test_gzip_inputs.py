"""A gzip-compressed input file reads as the plain file it holds, told by its first two bytes."""

from __future__ import annotations

import base64
import gzip
import json
import sys
from pathlib import Path

import pytest

from gbdraw import cli, read_genbank, session_io
from gbdraw.analysis.depth import parse_depth_tsv
from gbdraw.exceptions import InputFileError, ParseError
from gbdraw.io import input_text
from gbdraw.io.comparisons import read_comparison_table
from gbdraw.io.genome import load_gbks, parse_gff_fasta
from gbdraw.io.input_text import is_gzip_file, open_input_text, read_input_bytes
from gbdraw.session_io import load_session
from gbdraw.web_support.error_adapter import serialize_web_error

INPUTS = Path(__file__).resolve().parent / "test_inputs"
EXAMPLES = Path(__file__).resolve().parents[1] / "examples"
GENBANK = INPUTS / "NC_001416.gb"


def _gzipped(source: Path, target: Path) -> str:
    target.write_bytes(gzip.compress(source.read_bytes()))
    return str(target)


def _facts(records) -> list[tuple[str, int, int]]:
    return [(record.id, len(record.seq), len(record.features)) for record in records]


def test_a_gzip_genbank_file_reads_as_the_plain_one(tmp_path: Path) -> None:
    gzipped = _gzipped(GENBANK, tmp_path / "NC_001416.gb.gz")
    assert _facts(load_gbks([gzipped])) == _facts(load_gbks([str(GENBANK)]))
    assert _facts(read_genbank(gzipped)) == _facts(read_genbank(GENBANK))


def test_a_gzip_gff3_and_fasta_pair_reads_as_the_plain_pair(tmp_path: Path) -> None:
    gff = _gzipped(INPUTS / "NC_013668.gff3", tmp_path / "NC_013668.gff3.gz")
    fasta = _gzipped(INPUTS / "NC_013668.fasta", tmp_path / "NC_013668.fna.gz")
    plain = parse_gff_fasta(str(INPUTS / "NC_013668.gff3"), str(INPUTS / "NC_013668.fasta"))
    assert _facts(parse_gff_fasta(gff, fasta)) == _facts(plain)


def test_gzip_tables_read_as_the_plain_ones(tmp_path: Path) -> None:
    blast = EXAMPLES / "MjeNMV.MelaMJNV.tblastx.out"
    assert read_comparison_table(_gzipped(blast, tmp_path / "blast.tsv.gz")).equals(read_comparison_table(str(blast)))
    depth = INPUTS / "MjeNMV.DRR271272.depth.tsv"
    assert parse_depth_tsv(_gzipped(depth, tmp_path / "depth.tsv.gz")).equals(parse_depth_tsv(str(depth)))


def test_the_bytes_tell_gzip_from_plain_not_the_name(tmp_path: Path) -> None:
    plain_named_gz = tmp_path / "plain.gb.gz"
    plain_named_gz.write_bytes(GENBANK.read_bytes())
    gzip_named_gb = Path(_gzipped(GENBANK, tmp_path / "compressed.gb"))
    assert (is_gzip_file(plain_named_gz), is_gzip_file(gzip_named_gb)) == (False, True)
    for path in (plain_named_gz, gzip_named_gb):
        assert read_input_bytes(path) == GENBANK.read_bytes()
        with open_input_text(path) as handle:
            assert handle.readline().startswith("LOCUS")


@pytest.mark.parametrize("damage", ["truncated", "corrupt"])
def test_a_gzip_genbank_file_cut_off_or_damaged_is_unreadable_input(tmp_path: Path, damage: str) -> None:
    compressed = gzip.compress(GENBANK.read_bytes())
    path = tmp_path / "PRIVATE.gb.gz"
    path.write_bytes(compressed[: len(compressed) // 2] if damage == "truncated" else b"\x1f\x8b\x08\x00" + bytes(range(128, 256)))
    with pytest.raises(ParseError) as caught:
        load_gbks([str(path)])
    assert caught.value.diagnostic == {"code": "INPUT_UNREADABLE", "reason": "GENBANK_PARSE"}


def _base64_resources(value, found=None) -> list[dict]:
    """Every Session resource descriptor (``encoding: base64``) anywhere in the document."""
    found = [] if found is None else found
    if isinstance(value, dict):
        if value.get("encoding") == "base64" and isinstance(value.get("data"), str):
            found.append(value)
        for member in value.values():
            _base64_resources(member, found)
    elif isinstance(value, list):
        for member in value:
            _base64_resources(member, found)
    return found


def test_the_cli_embeds_the_plain_file_of_a_gzip_input_in_the_session(gbdraw_runner, tmp_path: Path) -> None:
    gzipped = Path(_gzipped(GENBANK, tmp_path / "NC_001416.gb.gz"))
    session_path = tmp_path / "gzip.gbdraw-session.json"
    code, output, svg_path = gbdraw_runner.run(
        "circular", [gzipped], "gzip", tmp_path, extra_args=["--session_output", str(session_path)]
    )
    assert code == 0, output
    assert svg_path.is_file()
    genbank = [resource for resource in _base64_resources(load_session(session_path)) if resource.get("kind") == "genbank"]
    assert len(genbank) == 1
    assert base64.b64decode(genbank[0]["data"]) == GENBANK.read_bytes()
    assert genbank[0]["name"] == "NC_001416.gb.gz"


# A gzip input expands to at most MAX_EXPANDED_INPUT_BYTES; past it every reader refuses the file by name
# (PD-OI-116 in docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md).
# The tests lower the cap, so their files stay small.
SMALL_CAP = 4096


@pytest.fixture
def small_cap(monkeypatch: pytest.MonkeyPatch) -> int:
    monkeypatch.setattr(input_text, "MAX_EXPANDED_INPUT_BYTES", SMALL_CAP)
    return SMALL_CAP


def _gzip_of(path: Path, content: bytes) -> Path:
    path.write_bytes(gzip.compress(content))
    return path


def _assert_refused(caught: pytest.ExceptionInfo[BaseException], path: Path) -> None:
    assert isinstance(caught.value, InputFileError)
    message = str(caught.value)
    assert str(path) in message
    assert f"{SMALL_CAP:,} bytes" in message
    assert "Decompress it" in message
    assert caught.value.diagnostic == {"code": "INPUT_UNREADABLE", "reason": "GZIP_EXPANDED_SIZE", "limit": 0}


def test_the_cap_is_512_mib() -> None:
    assert input_text.MAX_EXPANDED_INPUT_BYTES == 536_870_912


def test_a_gzip_input_past_the_cap_is_refused_by_name(small_cap: int, tmp_path: Path) -> None:
    path = _gzip_of(tmp_path / "PRIVATE-big.gb.gz", b"x" * (small_cap + 1))
    with pytest.raises(InputFileError) as caught:
        read_input_bytes(path)
    _assert_refused(caught, path)
    with pytest.raises(InputFileError) as caught, open_input_text(path) as handle:
        handle.read()
    _assert_refused(caught, path)


def test_a_gzip_input_of_exactly_the_cap_and_a_larger_plain_file_read(small_cap: int, tmp_path: Path) -> None:
    content = b"x" * small_cap
    at_cap = _gzip_of(tmp_path / "at-cap.gb.gz", content)
    assert read_input_bytes(at_cap) == content
    with open_input_text(at_cap) as handle:
        assert handle.read() == content.decode()
    plain = tmp_path / "plain.gb"
    plain.write_bytes(b"x" * (small_cap * 4))
    assert len(read_input_bytes(plain)) == small_cap * 4
    with open_input_text(plain) as handle:
        assert len(handle.read()) == small_cap * 4


def test_every_reader_lets_the_refusal_through_unchanged(small_cap: int, tmp_path: Path) -> None:
    """No reader rewords the refusal as a generic unreadable input."""

    readers = [
        (GENBANK, lambda path: load_gbks([str(path)])),
        (INPUTS / "MjeNMV.DRR271272.depth.tsv", lambda path: parse_depth_tsv(str(path))),
        (EXAMPLES / "MjeNMV.MelaMJNV.tblastx.out", lambda path: read_comparison_table(str(path))),
        (INPUTS / "NC_013668.fasta", lambda path: parse_gff_fasta(str(INPUTS / "NC_013668.gff3"), str(path))),
    ]
    for source, read in readers:
        path = _gzip_of(tmp_path / f"{source.name}.gz", source.read_bytes())
        with pytest.raises(InputFileError) as caught:
            read(path)
        _assert_refused(caught, path)
        assert serialize_web_error(caught.value, operation="generate", stage="render")["context"] == {
            "reason": "GZIP_EXPANDED_SIZE", "limit": 0,
        }


def test_a_gzip_session_is_not_capped(small_cap: int, tmp_path: Path) -> None:
    session = tmp_path / "big.gbdraw-session.json.gz"
    session.write_bytes(gzip.compress(json.dumps({"padding": "x" * (small_cap * 4)}).encode()))
    assert len(session_io._read_session_text(session)) > small_cap * 4


def test_the_cli_names_the_refused_file_and_the_cap(
    small_cap: int, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    path = _gzipped(GENBANK, tmp_path / "big.gb.gz")
    monkeypatch.setattr(sys, "argv", ["gbdraw", "circular", "--gbk", path, "-o", str(tmp_path / "out"), "-f", "svg"])
    with pytest.raises(SystemExit) as caught:
        cli.main()
    assert caught.value.code == 1
    error = capsys.readouterr().err
    assert error.startswith(f"ERROR: {path}: ")
    assert f"{small_cap:,} bytes" in error
