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
    # The label the flat file names (DEFINITION); a FASTA file names none and the
    # Web keeps its file-name rule (the helper sees only a staged file name).
    assert [result["recordLabel"] for result in results] == [
        None,
        "Danio rerio mitochondrion, complete genome",
        "Danio rerio mitochondrion, complete genome",
    ]
    assert "label" not in results[0]


def _flat_file_variant(path: Path, *, definition: str, organism: str | None) -> Path:
    text = DANIO_GENBANK.read_text(encoding="utf-8")
    head, _sep, tail = text.partition("DEFINITION")
    rest = tail.split("\nACCESSION", 1)[1]
    text = head + f"DEFINITION  {definition}\nACCESSION" + rest
    if organism is None:
        # Drop the SOURCE block (SOURCE, ORGANISM and its lineage lines).
        kept, in_source = [], False
        for line in text.splitlines():
            if line.startswith("SOURCE"):
                in_source = True
            elif in_source and line[:1] not in (" ", ""):
                in_source = False
            if not in_source:
                kept.append(line)
        text = "\n".join(kept) + "\n"
    elif organism != "Danio rerio":
        text = text.replace("  ORGANISM  Danio rerio", f"  ORGANISM  {organism}", 1)
    path.write_text(text, encoding="utf-8")
    return path


@pytest.mark.parametrize(
    ("definition", "organism", "expected"),
    [
        ("Danio rerio mitochondrion, complete genome.", "Danio rerio",
         "Danio rerio mitochondrion, complete genome"),
        (".", "Danio rerio", "Danio rerio"),
        (".", "Synthetic organism c", "Synthetic organism c"),
        (".", None, None),
    ],
)
def test_helper_record_label_is_the_cli_ring_default(
    tmp_path: Path, definition: str, organism: str | None, expected: str | None
) -> None:
    """D12: the Web ring default is the CLI default the flat file names itself."""

    path = _flat_file_variant(tmp_path / "ring.gb", definition=definition, organism=organism)
    result = json.loads(read_comparison_sequence_json(str(path)))
    assert result["recordLabel"] == expected
    cli_label = read_comparison_sequence_file(path).label
    assert cli_label == (expected if expected is not None else "ring.gb")


def test_helper_record_label_reads_the_first_record(tmp_path: Path) -> None:
    first = _flat_file_variant(tmp_path / "first.gb", definition="First record.", organism="Danio rerio")
    second = _flat_file_variant(tmp_path / "second.gb", definition="Second record.", organism="Danio rerio")
    path = tmp_path / "two-records.gb"
    path.write_text(first.read_text(encoding="utf-8") + second.read_text(encoding="utf-8"), encoding="utf-8")
    result = json.loads(read_comparison_sequence_json(str(path)))
    assert len(result["recordIds"]) == 2
    assert result["recordLabel"] == read_comparison_sequence_file(path).label == "First record"


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
