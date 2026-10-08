"""The comparison-genome sequence reader and --conservation_sequence (design D12, D18)."""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import pytest

from gbdraw.api import CircularDiagramOptions, read_conservation_table
from gbdraw.circular import circular_main
from gbdraw.exceptions import ValidationError
from gbdraw.io.comparison_sequences import read_comparison_sequence_file
from gbdraw.session_io import canonicalize_cli_invocation

REPO_ROOT = Path(__file__).resolve().parents[1]
TUTORIAL = REPO_ROOT / "gbdraw" / "web" / "tutorial-data"
HUMAN = TUTORIAL / "human-mitochondrion" / "HmmtDNA.gbk"
DANIO_FASTA = TUTORIAL / "metazoan-mitochondria-comparison" / "NC_002333.2.fna"
DANIO_GENBANK = TUTORIAL / "metazoan-mitochondria-four" / "NC_002333.2.gb"
DANIO_TSV = TUTORIAL / "metazoan-mitochondria-comparison" / "danio-human.tlosatx.tsv"
V39_SESSION = REPO_ROOT / "tests" / "fixtures" / "sessions" / "conservation-fasta.v39.gbdraw-session.json.gz"


def ddbj_flat_file(path: Path, *, accession: str, sequence: str, definition: str) -> Path:
    """A DDBJ getentry-style flat file (DDBJ LOCUS layout, no BASE COUNT)."""

    lines = [
        f"LOCUS       {accession.split('.')[0]:<16} {len(sequence):>11} bp    DNA     circular VRT 01-JAN-2024",
        f"DEFINITION  {definition}.",
        f"ACCESSION   {accession.split('.')[0]}",
        f"VERSION     {accession}",
        "KEYWORDS    .",
        "SOURCE      mitochondrion Danio rerio (zebrafish)",
        "  ORGANISM  Danio rerio",
        "            Eukaryota; Metazoa; Chordata; Craniata; Vertebrata.",
        "FEATURES             Location/Qualifiers",
        f"     source          1..{len(sequence)}",
        '                     /organism="Danio rerio"',
        '                     /mol_type="genomic DNA"',
        "ORIGIN      ",
    ]
    lowered = sequence.lower()
    for start in range(0, len(lowered), 60):
        chunk = lowered[start:start + 60]
        blocks = " ".join(chunk[offset:offset + 10] for offset in range(0, len(chunk), 10))
        lines.append(f"{start + 1:>9} {blocks}")
    lines.append("//")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _danio_ddbj(tmp_path: Path) -> Path:
    fasta = read_comparison_sequence_file(DANIO_FASTA)
    return ddbj_flat_file(
        tmp_path / "NC_002333.2.ddbj",
        accession="NC_002333.2",
        sequence=str(fasta.records[0].seq),
        definition="Danio rerio mitochondrion, complete genome",
    )


def _diagnostic(excinfo: pytest.ExceptionInfo) -> dict:
    return dict(excinfo.value.diagnostic or {})


# The one reader (D12).


def test_fasta_genbank_and_ddbj_give_the_same_records(tmp_path: Path) -> None:
    fasta = read_comparison_sequence_file(DANIO_FASTA)
    genbank = read_comparison_sequence_file(DANIO_GENBANK)
    ddbj = read_comparison_sequence_file(_danio_ddbj(tmp_path))

    assert (fasta.format, genbank.format, ddbj.format) == ("fasta", "genbank", "genbank")
    for parsed in (genbank, ddbj):
        assert [(r.id, r.name, str(r.seq).upper()) for r in parsed.records] == [
            (r.id, r.name, str(r.seq).upper()) for r in fasta.records
        ]
    assert fasta.records[0].id == "NC_002333.2" and len(fasta.records[0]) == 16_596
    # Default ring labels: FASTA file name without the extension (D-03); flat file DEFINITION.
    assert fasta.label == "NC_002333.2"
    assert genbank.label == ddbj.label == "Danio rerio mitochondrion, complete genome"


def test_flat_file_label_falls_back_to_the_organism(tmp_path: Path) -> None:
    text = DANIO_GENBANK.read_text(encoding="utf-8")
    head, _sep, tail = text.partition("DEFINITION")
    rest = tail.split("\nACCESSION", 1)[1]
    path = tmp_path / "no-definition.gb"
    path.write_text(head + "DEFINITION  .\nACCESSION" + rest, encoding="utf-8")
    assert read_comparison_sequence_file(path).label == "Danio rerio"


_LABEL_CASES = json.loads(
    (REPO_ROOT / "tests" / "fixtures" / "comparison_ring_default_label_cases.json").read_text(encoding="utf-8")
)["cases"]


@pytest.mark.parametrize("case", _LABEL_CASES, ids=[case["fileName"] for case in _LABEL_CASES])
def test_unlabelled_file_takes_the_shared_file_stem(tmp_path: Path, case: dict) -> None:
    # The Web defaultConservationSeriesLabel runs the same vectors.
    path = tmp_path / case["fileName"]
    path.write_text(">s1\nACGT\n", encoding="utf-8")
    assert read_comparison_sequence_file(path).label == case["expected"]
    path.write_text(DANIO_GENBANK.read_text(encoding="utf-8").replace(
        "DEFINITION  Danio rerio mitochondrion, complete genome.", "DEFINITION  ."
    ).replace("ORGANISM  Danio rerio", "ORGANISM  ."), encoding="utf-8")
    flat = read_comparison_sequence_file(path)
    assert (flat.format, flat.record_label, flat.label) == ("genbank", None, case["expected"])


@pytest.mark.parametrize("variant", ["empty-origin", "contig-only"])
def test_flat_file_without_sequence_is_sequence_missing(tmp_path: Path, variant: str) -> None:
    text = DANIO_GENBANK.read_text(encoding="utf-8")
    head = text.split("\nORIGIN", 1)[0]
    body = (
        "\nORIGIN      \n//\n"
        if variant == "empty-origin"
        else "\nCONTIG      join(AB000001.1:1..16596)\n//\n"
    )
    path = tmp_path / f"{variant}.gb"
    path.write_text(head + body, encoding="utf-8")
    with pytest.raises(ValidationError, match="has no sequence") as excinfo:
        read_comparison_sequence_file(path)
    assert _diagnostic(excinfo) == {
        "code": "INPUT_UNREADABLE",
        "field": "comparison_sequence",
        "reason": "SEQUENCE_MISSING",
    }


def test_unknown_format_is_unreadable(tmp_path: Path) -> None:
    path = tmp_path / "table.tsv"
    path.write_text("qseqid\tsseqid\n", encoding="utf-8")
    with pytest.raises(ValidationError, match="not FASTA") as excinfo:
        read_comparison_sequence_file(path)
    assert _diagnostic(excinfo)["code"] == "INPUT_UNREADABLE"


@pytest.mark.circular
def test_genbank_sequence_source_gives_the_same_interactive_svg(tmp_path: Path) -> None:
    outputs = []
    for name, source in (("fasta", DANIO_FASTA), ("genbank", DANIO_GENBANK)):
        prefix = tmp_path / name / "rings"
        prefix.parent.mkdir()
        circular_main([
            "--gbk", str(HUMAN), "--conservation_blast", str(DANIO_TSV),
            "--conservation_sequence", str(source), "--conservation_reference", "subject",
            "--conservation_labels", "Danio", "--identity", "40",
            "-f", "interactive_svg", "-o", str(prefix),
        ])
        outputs.append(sorted(prefix.parent.glob("*.svg")))
    fasta_files, genbank_files = outputs
    assert [p.name for p in fasta_files] == [p.name for p in genbank_files]
    assert any(p.name.endswith(".interactive.svg") for p in fasta_files)
    for left, right in zip(fasta_files, genbank_files):
        assert left.read_bytes() == right.read_bytes()


# --conservation_fasta and comparison_fasta are retired (D18, OD-2).


@pytest.mark.circular
def test_conservation_fasta_is_rejected_naming_the_replacement(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    with pytest.raises(SystemExit) as excinfo:
        circular_main([
            "--gbk", str(HUMAN), "--conservation_blast", str(DANIO_TSV),
            "--conservation_fasta", str(DANIO_FASTA), str(DANIO_FASTA),
            "-o", str(tmp_path / "out"),
        ])
    assert excinfo.value.code == 2
    assert "--conservation_fasta was retired; use --conservation_sequence" in capsys.readouterr().err


def test_comparison_fasta_column_is_rejected_naming_the_replacement(tmp_path: Path) -> None:
    table = tmp_path / "rings.tsv"
    table.write_text(f"blast\tcomparison_fasta\n{DANIO_TSV}\t{DANIO_FASTA}\n", encoding="utf-8")
    with pytest.raises(ValidationError, match="comparison_fasta was retired; use comparison_sequence"):
        read_conservation_table(str(table))

    table.write_text(f"blast\tcomparison_sequence\n{DANIO_TSV}\t{DANIO_GENBANK}\n", encoding="utf-8")
    parsed = read_conservation_table(str(table))
    assert parsed.comparison_sequence_files == [str(DANIO_GENBANK)]
    assert [dep.column for dep in parsed.path_dependencies] == ["blast", "comparison_sequence"]


def test_python_field_is_renamed_without_alias() -> None:
    with pytest.raises(TypeError):
        CircularDiagramOptions(conservation_fasta_files=[str(DANIO_FASTA)])  # type: ignore[call-arg]
    options = CircularDiagramOptions(
        conservation_blast_files=[str(DANIO_TSV)],
        conservation_sequence_files=[str(DANIO_GENBANK)],
    )
    assert tuple(options.conservation_sequence_files) == (str(DANIO_GENBANK),)


@pytest.mark.circular
def test_saved_session_keeps_the_persisted_wire_name(tmp_path: Path) -> None:
    prefix = tmp_path / "saved"
    circular_main([
        "--gbk", str(HUMAN), "--conservation_blast", str(DANIO_TSV),
        "--conservation_sequence", str(DANIO_GENBANK), "--conservation_reference", "subject",
        "-o", str(prefix), "-f", "svg", "--save_session",
    ])
    session = json.loads(prefix.with_suffix(".gbdraw-session.json").read_text(encoding="utf-8"))
    options = session["renderRequest"]["diagramOptions"]
    assert options["conservationFastaFiles"] == [
        {"resourceId": "conservation-fasta-files-1", "representation": "file"}
    ]
    assert "conservationSequenceFiles" not in options
    assert "--conservation_sequence" in session["cliInvocation"]["args"]


@pytest.mark.circular
def test_main_v39_session_with_conservation_fasta_replays(tmp_path: Path) -> None:
    """Main first-parent sessions (v32-v39) may record --conservation_fasta."""

    session_path = tmp_path / "v39.gbdraw-session.json"
    session_path.write_bytes(gzip.decompress(V39_SESSION.read_bytes()))
    session = json.loads(session_path.read_text(encoding="utf-8"))
    replay = tmp_path / "replay"
    circular_main(["--session", str(session_path), "-o", str(replay), "-f", "svg"])
    assert replay.with_suffix(".svg").is_file()

    recorded = session["cliInvocation"]["args"]
    assert recorded[recorded.index("--conservation_fasta") + 1] == "NC_002333.2.fna"
    args, _bindings = canonicalize_cli_invocation(recorded, (), mode="circular")
    assert "--conservation_fasta" not in args
    assert args[args.index("--conservation_sequence") + 1] == "NC_002333.2.fna"
    inline, _ = canonicalize_cli_invocation(["--conservation_fasta=a.fna"], (), mode="circular")
    assert inline == ["--conservation_sequence=a.fna"]
