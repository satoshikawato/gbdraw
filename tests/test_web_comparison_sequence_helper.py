"""The Web worker helper for Circular ring comparison files (design D12, PR-6).

The Web reads a ring's comparison file through the one Python reader and hashes
the same LOSAT query FASTA that the CLI ring search hashes, so FASTA, GenBank
and DDBJ files of one sequence give one raw key in the CLI and the Web.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from gbdraw.circular import circular_main
from gbdraw.comparisons.losat_jobs import nucleotide_losat_cache_key, sha256_text
from gbdraw.exceptions import ValidationError
from gbdraw.io.comparison_sequences import read_comparison_sequence_file
from gbdraw.web_support.comparison_sequences import read_comparison_sequence_json
from gbdraw.web_support.error_adapter import serialize_web_error
from tests.test_comparison_sequences import ddbj_flat_file
from tests.test_losat_linear_nucleotide import _fake_losat, _fresh_probe_caches, _searches  # noqa: F401

REPO_ROOT = Path(__file__).resolve().parents[1]
TUTORIAL = REPO_ROOT / "gbdraw" / "web" / "tutorial-data"
HUMAN = TUTORIAL / "human-mitochondrion" / "HmmtDNA.gbk"
DANIO_FASTA = TUTORIAL / "metazoan-mitochondria-comparison" / "NC_002333.2.fna"
DANIO_GENBANK = TUTORIAL / "metazoan-mitochondria-four" / "NC_002333.2.gb"


def _loose_fasta(path: Path) -> Path:
    """The Danio genome as lowercase 70-column FASTA with a description."""

    record = read_comparison_sequence_file(DANIO_FASTA).records[0]
    sequence = str(record.seq).lower()
    lines = [f">{record.id} Danio rerio mitochondrion, complete genome"]
    lines += [sequence[start:start + 70] for start in range(0, len(sequence), 70)]
    path.write_text("\r\n".join(lines) + "\r\n", encoding="utf-8")
    return path


def _inputs(tmp_path: Path) -> list[Path]:
    record = read_comparison_sequence_file(DANIO_FASTA).records[0]
    ddbj = ddbj_flat_file(
        tmp_path / "danio.ddbj",
        accession="NC_002333.2",
        sequence=str(record.seq),
        definition="Danio rerio mitochondrion, complete genome",
    )
    return [_loose_fasta(tmp_path / "danio-loose.fa"), DANIO_GENBANK, ddbj]


def test_helper_gives_one_query_for_fasta_genbank_and_ddbj(tmp_path: Path) -> None:
    results = [json.loads(read_comparison_sequence_json(str(path))) for path in _inputs(tmp_path)]

    assert [result["format"] for result in results] == ["fasta", "genbank", "genbank"]
    assert len({result["fasta"] for result in results}) == 1
    fasta = results[0]["fasta"]
    assert fasta.startswith(">NC_002333.2\n")
    body = fasta.splitlines()[1:]
    assert all(len(line) <= 60 for line in body) and "".join(body).isupper()
    assert [result["recordIds"] for result in results] == [["NC_002333.2"]] * 3
    # The default label of the shared reader (FASTA: file name; flat file: DEFINITION).
    assert [result["label"] for result in results] == [
        "danio-loose.fa",
        "Danio rerio mitochondrion, complete genome",
        "Danio rerio mitochondrion, complete genome",
    ]


@pytest.mark.circular
def test_helper_query_has_the_cli_ring_raw_key(tmp_path: Path) -> None:
    inputs = _inputs(tmp_path)
    helper_hashes = {
        sha256_text(json.loads(read_comparison_sequence_json(str(path)))["fasta"])
        for path in inputs
    }
    losat, log = _fake_losat(tmp_path)
    prefix = tmp_path / "rings"
    circular_main([
        "--gbk", str(HUMAN), "--losat", "losatn", "--losat_bin", losat,
        "--conservation_sequence", *map(str, inputs),
        "--save_session", "-o", str(prefix), "-f", "svg",
    ])
    # One genome in three layouts is one raw key: one search, one entry.
    assert len(_searches(log)) == 1
    session = json.loads(prefix.with_suffix(".gbdraw-session.json").read_text(encoding="utf-8"))
    entries = session["losatCache"]["entries"]
    assert {entry["queryCanonicalHash"] for entry in entries} == helper_hashes
    for entry in entries:
        assert entry["key"] == nucleotide_losat_cache_key(
            program="blastn",
            args=entry["args"],
            query_hash=entry["queryCanonicalHash"],
            subject_hash=entry["subjectCanonicalHash"],
            flow="circular-conservation",
        )


@pytest.mark.parametrize("content", ["", "LOCUS       EMPTY 0 bp DNA\nORIGIN\n//\n"])
def test_helper_rejects_a_file_without_sequence(tmp_path: Path, content: str) -> None:
    path = tmp_path / "empty.gb"
    path.write_text(content, encoding="utf-8")
    with pytest.raises(ValidationError) as excinfo:
        read_comparison_sequence_json(str(path))
    assert excinfo.value.diagnostic["reason"] == "SEQUENCE_MISSING"
    serialized = serialize_web_error(
        excinfo.value, operation="readComparisonSequence", stage="helper"
    )
    assert serialized["code"] == "INPUT_UNREADABLE"
    assert serialized["operation"] == "readComparisonSequence"
    assert serialized["context"]["reason"] == "SEQUENCE_MISSING"
