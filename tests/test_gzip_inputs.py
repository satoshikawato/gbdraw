"""A gzip-compressed input file reads as the plain file it holds, told by its first two bytes."""

from __future__ import annotations

import base64
import gzip
import json
from pathlib import Path

import pytest

from gbdraw import read_genbank
from gbdraw.analysis.depth import parse_depth_tsv
from gbdraw.exceptions import ParseError
from gbdraw.io.comparisons import read_comparison_table
from gbdraw.io.genome import load_gbks, parse_gff_fasta
from gbdraw.io.input_text import is_gzip_file, open_input_text, read_input_bytes
from gbdraw.session_io import load_session

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
