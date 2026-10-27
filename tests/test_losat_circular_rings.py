"""Circular similarity rings from LOSATN / TLOSATX (design 3.3, 3.4, 3.7, D12, PR-4)."""

from __future__ import annotations

import base64
from dataclasses import replace
import json
from pathlib import Path
import re

import pytest

import gbdraw.comparisons.losat_runtime as runtime_module
from gbdraw import (
    CircularOptions,
    ComparisonRingOptions,
    ComparisonRingTrackOptions,
    draw_circular,
    read_genbank,
)
from gbdraw.api import (
    CircularDiagramOptions,
    CircularDiagramRequest,
    LosatSearchOptions,
    RecordInput,
    RenderOutputRequest,
    render_request,
)
from gbdraw.api.requests import GenBankInputSource
from gbdraw.circular import circular_main
from gbdraw.comparisons.losat_runtime import resolve_losat_runtime
from gbdraw.exceptions import ValidationError
from gbdraw.session_request_codec import CanonicalRequestEncodingError, encode_canonical_request
from tests.test_comparison_sequences import ddbj_flat_file
from tests.test_losat_linear_nucleotide import _fake_losat, _fresh_probe_caches, _searches  # noqa: F401

REPO_ROOT = Path(__file__).resolve().parents[1]
TUTORIAL = REPO_ROOT / "gbdraw" / "web" / "tutorial-data"
HUMAN = TUTORIAL / "human-mitochondrion" / "HmmtDNA.gbk"
COMPARISON = TUTORIAL / "metazoan-mitochondria-comparison"
FOUR = TUTORIAL / "metazoan-mitochondria-four"
FASTAS = [COMPARISON / f"{accession}.fna" for accession in ("NC_002333.2", "NC_024511.2", "NC_001328.1")]
FROZEN_TSVS = [
    COMPARISON / name
    for name in (
        "danio-human.tlosatx.tsv",
        "drosophila-human.tlosatx.tsv",
        "caenorhabditis-human.tlosatx.tsv",
    )
]
T_CLI_09_SVG = REPO_ROOT / "docs" / "images" / "t-cli-09" / "precomputed_circular_rings.svg"
# T-CLI-09 labels, colors, thresholds and canvas options (design 3.2 example 5).
T_CLI_09_STYLE = [
    "--conservation_labels", "Danio rerio (NC_002333.2)",
    "Drosophila melanogaster (NC_024511.2)", "Caenorhabditis elegans (NC_001328.1)",
    "--conservation_colors", "#4E79A7", "#F28E2B", "#59A14F",
    "--bitscore", "50", "--evalue", "1e-5", "--identity", "40", "--alignment_length", "50",
    "--conservation_ring_width", "18", "--conservation_ring_gap", "4",
    "--species", "<i>Homo sapiens</i>",
    "--qualifier_priority", str(TUTORIAL / "shared" / "cds_gene_qualifier_priority.tsv"),
    "--track_type", "middle", "--labels", "out", "--definition_font_size", "18",
    "--plot_title", "Precomputed TLOSATX rings around Homo sapiens mtDNA",
    "--plot_title_position", "bottom", "--legend", "right", "-f", "svg",
]


def _session(prefix: Path) -> dict:
    return json.loads(prefix.with_suffix(".gbdraw-session.json").read_text(encoding="utf-8"))


def _diagnostic(excinfo: pytest.ExceptionInfo) -> dict:
    return dict(excinfo.value.diagnostic or {})


@pytest.fixture(scope="module")
def native_losat() -> None:
    try:
        runtime = resolve_losat_runtime("tlosatx")
    except ValidationError as error:
        pytest.skip(f"No native LOSAT runtime resolves: {error}")
    if runtime.kind != "losat":
        pytest.skip(f"Resolved runtime is {runtime.kind}, not native LOSAT.")


# Acceptance with the resolved native LOSAT (design 4 PR-4).


@pytest.mark.circular
def test_cli_tlosatx_reproduces_the_t_cli_09_figure(native_losat, tmp_path: Path) -> None:
    prefix = tmp_path / "precomputed_circular_rings"
    raw = tmp_path / "raw"
    circular_main([
        "--gbk", str(HUMAN), "--losat", "tlosatx", "--losat_gencode", "2",
        "--conservation_sequence", *map(str, FASTAS),
        "--conservation_losat_gencode", "2", "5", "5",
        *T_CLI_09_STYLE, "--losat_output_dir", str(raw), "--save_session",
        "-o", str(prefix),
    ])

    assert prefix.with_suffix(".svg").read_bytes() == T_CLI_09_SVG.read_bytes()
    names = [f"{path.stem}.circular_conservation.tlosatx.tsv" for path in FASTAS]
    for name, frozen in zip(names, FROZEN_TSVS):
        assert (raw / name).read_bytes() == frozen.read_bytes()
    manifest = (raw / "conservation.tsv").read_text(encoding="utf-8").splitlines()
    assert manifest[0] == "blast\tcomparison_sequence\tlabel\tcolor"
    assert [row.split("\t")[0] for row in manifest[1:]] == names

    # The manifest is a --conservation_table without LOSAT.
    again = tmp_path / "from-table"
    style = [token for token in T_CLI_09_STYLE]
    for flag in ("--conservation_labels", "--conservation_colors"):
        index = style.index(flag)
        del style[index:index + 4]
    circular_main([
        "--gbk", str(HUMAN), "--conservation_table", str(raw / "conservation.tsv"),
        "--conservation_reference", "subject", *style, "-o", str(again),
    ])
    assert again.with_suffix(".svg").read_bytes() == T_CLI_09_SVG.read_bytes()

    session = _session(prefix)
    entries = session["losatCache"]["entries"]
    assert [entry["flow"] for entry in entries] == ["circular-conservation"] * 3
    assert [entry["args"] for entry in entries] == [
        ["--query-gencode", "2", "--db-gencode", "2"],
        ["--query-gencode", "5", "--db-gencode", "2"],
        ["--query-gencode", "5", "--db-gencode", "2"],
    ]
    assert entries[0]["program"] == "tblastx" and entries[0]["runtime"]["program"] == "tblastx"
    assert [entry["filename"] for entry in entries] == names
    options = session["renderRequest"]["diagramOptions"]
    assert len(options["conservationBlastFiles"]) == len(options["conservationFastaFiles"]) == 3
    assert options["conservationReference"] == "subject"


@pytest.mark.circular
def test_losatn_ring_is_identical_for_fasta_genbank_and_ddbj(native_losat, tmp_path: Path) -> None:
    """D12: the comparison genome may be FASTA, GenBank or DDBJ."""

    sequence = "".join(
        line.strip() for line in FASTAS[0].read_text(encoding="utf-8").splitlines()[1:]
    )
    ddbj = ddbj_flat_file(
        tmp_path / "NC_002333.2.ddbj",
        accession="NC_002333.2",
        sequence=sequence,
        definition="Danio rerio mitochondrion, complete genome",
    )
    svgs = []
    for name, source in (("fasta", FASTAS[0]), ("genbank", FOUR / "NC_002333.2.gb"), ("ddbj", ddbj)):
        prefix = tmp_path / name / "ring"
        prefix.parent.mkdir()
        circular_main([
            "--gbk", str(HUMAN), "--losat", "losatn", "--losatn_task", "blastn",
            "--conservation_sequence", str(source), "--conservation_labels", "Danio",
            "--identity", "60", "-f", "svg", "-o", str(prefix),
        ])
        svgs.append(prefix.with_suffix(".svg").read_bytes())
    assert svgs[0].count(b"data-gbdraw-match-id") > 0
    assert svgs[0] == svgs[1] == svgs[2]


@pytest.mark.circular
def test_genbank_without_sequence_stops_the_ring(tmp_path: Path) -> None:
    text = (FOUR / "NC_002333.2.gb").read_text(encoding="utf-8")
    empty = tmp_path / "empty.gb"
    empty.write_text(text.split("\nORIGIN", 1)[0] + "\nORIGIN      \n//\n", encoding="utf-8")
    losat, log = _fake_losat(tmp_path)
    with pytest.raises(ValidationError, match="has no sequence") as excinfo:
        circular_main([
            "--gbk", str(HUMAN), "--losat", "losatn", "--losat_bin", losat,
            "--conservation_sequence", str(empty), "-o", str(tmp_path / "out"),
        ])
    assert _diagnostic(excinfo) == {
        "code": "INPUT_UNREADABLE",
        "field": "comparison_sequence",
        "reason": "SEQUENCE_MISSING",
    }
    assert _searches(log) == []


# Session and replay.


@pytest.mark.circular
def test_saved_session_replays_without_losat(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    losat, log = _fake_losat(tmp_path)
    prefix = tmp_path / "saved"
    circular_main([
        "--gbk", str(HUMAN), "--losat", "tlosatx", "--losat_bin", losat,
        "--conservation_sequence", str(FASTAS[0]), str(FOUR / "NC_024511.2.gb"),
        "--save_session", "-o", str(prefix), "-f", "svg",
    ])
    searches = _searches(log)
    assert len(searches) == 2
    # TLOSATX tables default to 1 and are always passed (Web ring flow, D17).
    assert all("--query-gencode" in call or "-query_gencode" in call for call in searches)
    session = _session(prefix)
    entries = session["losatCache"]["entries"]
    assert [entry["args"] for entry in entries] == [["--query-gencode", "1", "--db-gencode", "1"]] * 2
    assert entries[0]["runtime"]["source"] == "explicit"
    labels = session["renderRequest"]["diagramOptions"]["conservationLabels"]
    # Default labels: FASTA file name without the extension (D-03); GenBank DEFINITION.
    assert labels == ["NC_002333.2", "Drosophila melanogaster mitochondrion, complete genome"]

    def no_runtime(*_args, **_kwargs):
        raise AssertionError("Session replay resolved a LOSAT runtime.")

    monkeypatch.setattr(runtime_module, "resolve_losat_runtime", no_runtime)
    replay = tmp_path / "replay"
    circular_main([
        "--session", str(prefix.with_suffix(".gbdraw-session.json")),
        "-o", str(replay), "-f", "svg",
    ])
    assert replay.with_suffix(".svg").read_bytes() == prefix.with_suffix(".svg").read_bytes()


@pytest.mark.circular
def test_rings_whose_files_share_a_stem_save_and_replay(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    # PR6-B1: X.fna and X.gb gave both ring TSVs one Session resource name.
    losat, _log = _fake_losat(tmp_path)
    prefix, raw = tmp_path / "saved", tmp_path / "raw"
    circular_main([
        "--gbk", str(HUMAN), "--losat", "losatn", "--losat_bin", losat,
        "--conservation_sequence", str(FASTAS[0]), str(FOUR / "NC_002333.2.gb"),
        "--losat_output_dir", str(raw), "--save_session", "-o", str(prefix), "-f", "svg",
    ])
    names = [
        "NC_002333.2.circular_conservation.losatn.tsv",
        "NC_002333.2.circular_conservation.losatn.2.tsv",
    ]
    manifest = (raw / "conservation.tsv").read_text(encoding="utf-8").splitlines()
    assert [row.split("\t")[0] for row in manifest[1:]] == names
    session = _session(prefix)
    rings = session["renderRequest"]["diagramOptions"]["conservationBlastFiles"]
    assert [session["resources"][ring["resourceId"]]["name"] for ring in rings] == names
    assert [(raw / name).read_text(encoding="utf-8") for name in names] == [
        base64.b64decode(session["resources"][ring["resourceId"]]["data"]).decode("utf-8")
        for ring in rings
    ]

    def no_runtime(*_args, **_kwargs):
        raise AssertionError("Session replay resolved a LOSAT runtime.")

    monkeypatch.setattr(runtime_module, "resolve_losat_runtime", no_runtime)
    replay = tmp_path / "replay"
    circular_main([
        "--session", str(prefix.with_suffix(".gbdraw-session.json")),
        "-o", str(replay), "-f", "svg",
    ])
    assert replay.with_suffix(".svg").read_bytes() == prefix.with_suffix(".svg").read_bytes()


@pytest.mark.circular
@pytest.mark.parametrize("losat", [True, False], ids=["losatn", "precomputed"])
def test_ring_files_with_one_basename_in_two_directories_save_and_replay(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, losat: bool
) -> None:
    # FIX-B2-N1: Session resources were named by basename, so a/X.fna and
    # b/X.fna (or a/ring.tsv and b/ring.tsv) stopped --save_session.
    for directory, fasta, frozen in (("a", FASTAS[0], FROZEN_TSVS[0]), ("b", FASTAS[1], FROZEN_TSVS[1])):
        (tmp_path / directory).mkdir()
        (tmp_path / directory / "X.fna").write_bytes(fasta.read_bytes())
        (tmp_path / directory / "ring.tsv").write_bytes(frozen.read_bytes())
    sequences = ["--conservation_sequence", str(tmp_path / "a" / "X.fna"), str(tmp_path / "b" / "X.fna")]
    if losat:
        rings = ["--losat", "losatn", "--losat_bin", _fake_losat(tmp_path)[0]]
        expected = {"X.fna", "X.2.fna"}
    else:
        # No --conservation_labels: each ring is labelled with its file stem (D-03).
        rings = ["--conservation_blast", str(tmp_path / "a" / "ring.tsv"), str(tmp_path / "b" / "ring.tsv")]
        expected = {"X.fna", "X.2.fna", "ring.tsv", "ring.2.tsv"}
    prefix = tmp_path / "saved"
    circular_main([
        "--gbk", str(HUMAN), *rings, *sequences, "--save_session", "-o", str(prefix), "-f", "svg",
    ])
    session = _session(prefix)
    assert expected <= {resource["name"] for resource in session["resources"].values()}
    if not losat:
        # The replay reads ring.2.tsv; the stored labels keep the drawn ones.
        assert session["renderRequest"]["diagramOptions"]["conservationLabels"] == ["ring", "ring"]

    def no_runtime(*_args, **_kwargs):
        raise AssertionError("Session replay resolved a LOSAT runtime.")

    monkeypatch.setattr(runtime_module, "resolve_losat_runtime", no_runtime)
    replay = tmp_path / "replay"
    circular_main([
        "--session", str(prefix.with_suffix(".gbdraw-session.json")),
        "-o", str(replay), "-f", "svg",
    ])
    assert replay.with_suffix(".svg").read_bytes() == prefix.with_suffix(".svg").read_bytes()


@pytest.mark.circular
def test_conservation_table_drives_ring_searches(tmp_path: Path) -> None:
    losat, log = _fake_losat(tmp_path)
    table = tmp_path / "rings.tsv"
    table.write_text(
        "comparison_sequence\tlabel\tlosat_gencode\n"
        f"{FASTAS[1]}\tDrosophila\t5\n{FASTAS[2]}\tCaenorhabditis\t\n",
        encoding="utf-8",
    )
    circular_main([
        "--gbk", str(HUMAN), "--losat", "tlosatx", "--losat_bin", losat,
        "--conservation_table", str(table), "--save_session",
        "-o", str(tmp_path / "out"), "-f", "svg",
    ])
    entries = _session(tmp_path / "out")["losatCache"]["entries"]
    assert [entry["args"][1] for entry in entries] == ["5", "1"]

    table.write_text(f"blast\tcomparison_sequence\n{FROZEN_TSVS[0]}\t{FASTAS[0]}\n", encoding="utf-8")
    with pytest.raises(ValidationError, match="column blast cannot be combined with --losat") as excinfo:
        circular_main([
            "--gbk", str(HUMAN), "--losat", "tlosatx", "--losat_bin", losat,
            "--conservation_table", str(table), "-o", str(tmp_path / "bad"),
        ])
    assert _diagnostic(excinfo)["reason"] == "RING_LOSAT_INPUT"


# Rejections (design 3.3, 3.7).


@pytest.mark.circular
@pytest.mark.parametrize(
    ("argv", "reason", "message"),
    [
        (["--losat", "losatp", "--conservation_sequence", str(FASTAS[0])],
         "RING_LOSAT_PROGRAM", "losatp is not available for rings"),
        (["--losat", "losatn", "--conservation_sequence", str(FASTAS[0]),
          "--conservation_blast", str(FROZEN_TSVS[0])],
         "RING_LOSAT_INPUT", "--conservation_blast cannot be combined with --losat"),
        (["--losat", "losatn", "--conservation_sequence", str(FASTAS[0]),
          "--conservation_reference", "query"],
         "RING_LOSAT_INPUT", "--conservation_reference must be auto or subject"),
        (["--losat", "losatn", "--conservation_sequence", str(FASTAS[0]),
          "--conservation_losat_gencode", "5"],
         "LOSAT_OPTION_PROGRAM", "--conservation_losat_gencode does not apply to --losat losatn"),
        (["--losat", "losatn", "--conservation_sequence", str(FASTAS[0]), "--losat_gencode", "2"],
         "LOSAT_OPTION_PROGRAM", "--losat_gencode does not apply to --losat losatn"),
        (["--losat", "tlosatx"],
         "RING_LOSAT_INPUT", "one comparison sequence file per ring"),
    ],
)
def test_ring_losat_rejections(tmp_path: Path, argv: list[str], reason: str, message: str) -> None:
    with pytest.raises(ValidationError, match=message) as excinfo:
        circular_main(["--gbk", str(HUMAN), *argv, "-o", str(tmp_path / "out")])
    assert _diagnostic(excinfo)["code"] == "COMPARISON_INPUT"
    assert _diagnostic(excinfo)["reason"] == reason


def test_typed_ring_options_reject_linear_only_search_fields() -> None:
    with pytest.raises(ValidationError, match="Linear diagrams only") as excinfo:
        CircularDiagramOptions(
            losat_search=LosatSearchOptions(program="losatn", pairs=((0, 1),)),
            conservation_sequence_files=[str(FASTAS[0])],
        )
    assert _diagnostic(excinfo)["reason"] == "RING_LOSAT_INPUT"


# Typed and introductory APIs.


@pytest.mark.circular
def test_typed_and_introductory_apis_run_ring_tlosatx(tmp_path: Path) -> None:
    losat, log = _fake_losat(tmp_path)
    from gbdraw.api import LosatRuntimeOptions

    typed = render_request(
        CircularDiagramRequest(
            records=(RecordInput(source=GenBankInputSource(path=str(HUMAN))),),
            options=CircularDiagramOptions(
                losat_search=LosatSearchOptions(
                    program="tlosatx",
                    record_gencodes=(2,),
                    runtime=LosatRuntimeOptions(losat_executable=losat),
                ),
                conservation_sequence_files=[str(FASTAS[0])],
                conservation_losat_gencodes=[2],
                conservation_labels=["Danio"],
            ),
            output=RenderOutputRequest(
                output_prefix="typed", output_directory=tmp_path, formats=("svg",)
            ),
        )
    )
    assert typed.request.options.losat_search is None
    assert typed.request.options.conservation_search_results[0].name == (
        "NC_002333.2.circular_conservation.tlosatx.tsv"
    )
    assert [entry["args"] for entry in typed.losat_cache_entries] == [
        ["--query-gencode", "2", "--db-gencode", "2"]
    ]

    diagram = draw_circular(
        read_genbank(HUMAN),
        options=CircularOptions(
            comparison_rings=ComparisonRingOptions(
                losat="tlosatx",
                reference_gencode=2,
                losat_executable=losat,
                tracks=[
                    ComparisonRingTrackOptions(
                        comparison_sequence_source=FASTAS[0], losat_gencode=2, label="Danio"
                    )
                ],
            )
        ),
    )
    intro = Path(diagram.save(tmp_path / "intro.svg")).read_bytes()
    typed_svg = (tmp_path / "typed.svg").read_bytes()
    assert intro.count(b"data-gbdraw-match-id") == typed_svg.count(b"data-gbdraw-match-id") > 0
    assert len(_searches(log)) == 2

    with pytest.raises(ValidationError, match="file path"):
        draw_circular(
            read_genbank(HUMAN),
            options=CircularOptions(
                comparison_rings=ComparisonRingOptions(
                    losat="losatn",
                    tracks=[ComparisonRingTrackOptions(comparison_sequence_source=read_genbank(HUMAN))],
                )
            ),
        )


def test_unresolved_ring_intent_is_never_encoded() -> None:
    request = CircularDiagramRequest(
        records=(RecordInput(source=GenBankInputSource(path=str(HUMAN))),),
        options=CircularDiagramOptions(
            losat_search=LosatSearchOptions(program="losatn"),
            conservation_sequence_files=[str(FASTAS[0])],
        ),
    )
    with pytest.raises(CanonicalRequestEncodingError, match="losatn ring search"):
        encode_canonical_request(request)


# D-03: a precomputed table ring without a label is named by its file name
# without the last extension, as in the Web; Sessions written before that
# store no label for it and keep their full-file-name label on replay.
PRECOMPUTED_TSV = FROZEN_TSVS[2]
UNLABELLED_RING_SESSIONS = {
    # main fe6861f0 CLI (Session 44) and release 0.13.0 CLI (Session 30), with the
    # label their own runs drew (OV-202: not the replay's temporary file name);
    # see tests/fixtures/sessions/precomputed-ring.provenance.json.
    "main-v44": (
        REPO_ROOT / "tests" / "fixtures" / "sessions" / "precomputed-ring-unlabelled-cli.v44.gbdraw-session.json.gz",
        "caenorhabditis-human.tlosatx.tsv",
    ),
    "release-v30": (
        REPO_ROOT / "tests" / "fixtures" / "sessions" / "precomputed-ring-unlabelled-cli.v30.gbdraw-session.json.gz",
        "caenorhabditis-human.tlosatx.tsv",
    ),
}


def _ring_label_texts(svg: bytes) -> set[str]:
    return set(re.findall(r'data-track-label="([^"]*)"', svg.decode("utf-8")))


@pytest.mark.circular
def test_an_unlabelled_precomputed_ring_takes_the_file_stem(tmp_path: Path) -> None:
    prefix = tmp_path / "fresh"
    circular_main([
        "--gbk", str(HUMAN), "--conservation_blast", str(PRECOMPUTED_TSV), "--conservation_reference", "subject",
        "--save_session", "-o", str(prefix), "-f", "svg",
    ])
    assert _ring_label_texts(prefix.with_suffix(".svg").read_bytes()) == {"caenorhabditis-human.tlosatx"}
    # The writer stores the label, so a Session never relies on the reader's default.
    assert _session(prefix)["renderRequest"]["diagramOptions"]["conservationLabels"] == ["caenorhabditis-human.tlosatx"]

    typed_prefix = "typed"
    render_request(
        CircularDiagramRequest(
            records=(RecordInput(source=GenBankInputSource(path=str(HUMAN))),),
            options=CircularDiagramOptions(
                conservation_blast_files=[str(PRECOMPUTED_TSV)], conservation_reference="subject"
            ),
            output=RenderOutputRequest(output_prefix=typed_prefix, output_directory=tmp_path, formats=("svg",)),
        )
    )
    assert _ring_label_texts((tmp_path / f"{typed_prefix}.svg").read_bytes()) == {"caenorhabditis-human.tlosatx"}


@pytest.mark.circular
@pytest.mark.parametrize("name", sorted(UNLABELLED_RING_SESSIONS))
def test_a_session_without_a_precomputed_ring_label_keeps_the_file_name(
    name: str, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import gzip

    monkeypatch.chdir(tmp_path)
    fixture, drawn_label = UNLABELLED_RING_SESSIONS[name]
    session = json.loads(gzip.decompress(fixture.read_bytes()))
    labelled = json.loads(json.dumps(session))
    if "renderRequest" in labelled:
        assert "conservationLabels" not in labelled["renderRequest"]["diagramOptions"]
        labelled["renderRequest"]["diagramOptions"]["conservationLabels"] = [drawn_label]
    else:
        args = labelled["cliInvocation"]["args"]
        assert "--conservation_labels" not in args
        index = args.index("-o")
        labelled["cliInvocation"]["args"] = [*args[:index], "--conservation_labels", drawn_label, *args[index:]]
    svgs = {}
    for variant, document in (("stored", session), ("labelled", labelled)):
        path = tmp_path / f"{variant}.gbdraw-session.json"
        path.write_text(json.dumps(document), encoding="utf-8")
        circular_main(["--session", str(path), "-o", variant, "-f", "svg"])
        svgs[variant] = (tmp_path / f"{variant}.svg").read_bytes()
    assert _ring_label_texts(svgs["stored"]) == {drawn_label}
    assert svgs["stored"] == svgs["labelled"]


# The subject database is the whole subject-side source file (Owner D-25, OV-1682).

_SUBJECT_SEQUENCES = {
    "SUBA": "acgtacgtac" * 7,
    "SUBB": "ggccttaagg" * 9,
    "SUBC": "ttgacatgca" * 8,
    "SUBD": "catgcatgaa" * 6,
}
_CON_RECORD = """LOCUS       SUBCON                    20 bp    DNA     linear   CON 01-JAN-1980
DEFINITION  contig record without a sequence.
ACCESSION   SUBCON
VERSION     SUBCON.1
KEYWORDS    .
SOURCE      .
  ORGANISM  .
            .
FEATURES             Location/Qualifiers
CONTIG      join(SUBA.1:1..10,SUBB.1:1..10)
//
"""


def _genbank_text(record_id: str, sequence: str) -> str:
    from io import StringIO

    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    record = SeqRecord(
        Seq(sequence), id=record_id, name=record_id, description=f"{record_id} record",
        annotations={"molecule_type": "DNA", "topology": "linear"},
    )
    handle = StringIO()
    SeqIO.write([record], handle, "genbank")
    return handle.getvalue()


def _subject_file(tmp_path: Path, *, name: str = "subject.gb", con: bool = True, bom: bool = False) -> Path:
    """SUBA, SUBB, a CON record (no sequence), SUBC."""

    texts = [_genbank_text("SUBA", _SUBJECT_SEQUENCES["SUBA"]), _genbank_text("SUBB", _SUBJECT_SEQUENCES["SUBB"])]
    if con:
        texts.append(_CON_RECORD)
    texts.append(_genbank_text("SUBC", _SUBJECT_SEQUENCES["SUBC"]))
    path = tmp_path / name
    path.write_bytes((b"\xef\xbb\xbf" if bom else b"") + "".join(texts).encode("ascii"))
    return path


def _expected_database(*record_ids: str) -> str:
    from gbdraw.comparisons.losat_jobs import losat_fasta_text

    return "".join(losat_fasta_text(record_id, _SUBJECT_SEQUENCES[record_id]) for record_id in record_ids)


@pytest.fixture
def captured_ring_searches(monkeypatch: pytest.MonkeyPatch) -> list[tuple[str, str]]:
    """The (query, subject) FASTA of every ring search; the search returns no rows."""

    import gbdraw.comparisons.circular_losat as circular_losat_module

    searches: list[tuple[str, str]] = []

    def fake_search(_spec, query_fasta, subject_fasta, **_kwargs):
        searches.append((query_fasta, subject_fasta))
        return ""

    monkeypatch.setattr(circular_losat_module, "run_losat_search", fake_search)
    return searches


def _ring_request(tmp_path: Path, *inputs: RecordInput) -> CircularDiagramRequest:
    query = tmp_path / "query.fasta"
    query.write_text(">QUERY1\nACGTACGTACGTACGTACGT\n", encoding="ascii")
    return CircularDiagramRequest(
        records=tuple(inputs),
        options=CircularDiagramOptions(
            losat_search=LosatSearchOptions(program="losatn"),
            conservation_sequence_files=[str(query)],
        ),
    )


@pytest.mark.circular
@pytest.mark.parametrize(
    "case",
    ["records not drawn", "single presentation", "crop", "reverse complement crop"],
)
def test_ring_database_is_the_whole_reference_file(
    case: str, tmp_path: Path, captured_ring_searches: list[tuple[str, str]]
) -> None:
    from gbdraw.api.request_render import plan_circular_request
    from gbdraw.api.requests import RecordPresentation
    from gbdraw.comparisons.losat_jobs import nucleotide_losat_cache_key, sha256_text
    from gbdraw.io.record_select import parse_record_selector
    from gbdraw.io.regions import parse_region_spec

    source = GenBankInputSource(path=str(_subject_file(tmp_path)))
    inputs = {
        # SUBB is OFF: only SUBA and SUBC are drawn.
        "records not drawn": (
            RecordInput(source=source, selector=parse_record_selector("#1")),
            RecordInput(source=source, selector=parse_record_selector("#4")),
        ),
        "single presentation": (RecordInput(source=source, selector=parse_record_selector("#2")),),
        "crop": (RecordInput(source=source, region=parse_region_spec("SUBB:11-50")),),
        "reverse complement crop": (
            RecordInput(
                source=source, region=parse_region_spec("SUBB:11-50"),
                presentation=RecordPresentation(reverse_complement=True),
            ),
        ),
    }[case]
    plan = plan_circular_request(_ring_request(tmp_path, *inputs))

    expected = _expected_database("SUBA", "SUBB", "SUBC")
    assert [subject for _query, subject in captured_ring_searches] == [expected]
    (entry,) = plan.losat_cache_entries
    assert entry["subjectCanonicalHash"] == sha256_text(expected)
    assert entry["key"] == nucleotide_losat_cache_key(
        program="blastn", args=("--task", "megablast"),
        query_hash=sha256_text(">QUERY1\nACGTACGTACGTACGTACGT\n"),
        subject_hash=sha256_text(expected), flow="circular-conservation",
    )


@pytest.mark.circular
def test_ring_database_reads_the_reference_file_once(
    tmp_path: Path, captured_ring_searches: list[tuple[str, str]], monkeypatch: pytest.MonkeyPatch
) -> None:
    """The planner's parse of the GenBank source is the database: no second parse in one request."""

    import gbdraw.io.genome as genome_module
    from gbdraw.api.request_render import plan_circular_request
    from gbdraw.io.record_select import parse_record_selector

    parses: list[str] = []
    real_parse = genome_module.SeqIO.parse

    def counting_parse(handle, fmt, *args, **kwargs):
        parses.append(fmt)
        return real_parse(handle, fmt, *args, **kwargs)

    monkeypatch.setattr(genome_module.SeqIO, "parse", counting_parse)
    source = GenBankInputSource(path=str(_subject_file(tmp_path)))
    plan_circular_request(_ring_request(tmp_path, RecordInput(source=source, selector=parse_record_selector("#2"))))
    # The comparison FASTA is parsed too; the reference GenBank file once.
    assert [fmt for fmt in parses if fmt == "genbank"] == ["genbank"]
    assert [subject for _query, subject in captured_ring_searches] == [_expected_database("SUBA", "SUBB", "SUBC")]


@pytest.mark.circular
def test_ring_database_of_gff3_is_every_fasta_entry(
    tmp_path: Path, captured_ring_searches: list[tuple[str, str]], monkeypatch: pytest.MonkeyPatch
) -> None:
    """GFF3 mode searches every entry of the FASTA file, also one the GFF3 file does not name."""

    import gbdraw.io.genome as genome_module
    from gbdraw.api.request_render import plan_circular_request
    from gbdraw.api.requests import GffFastaInputSource

    gff = tmp_path / "subject.gff3"
    gff.write_text(
        "##gff-version 3\nSUBB\t.\tgene\t5\t40\t.\t+\t.\tID=gene1\n",
        encoding="ascii",
    )
    fasta = tmp_path / "subject.fasta"
    fasta.write_text(
        "".join(f">{record_id} entry\n{_SUBJECT_SEQUENCES[record_id]}\n" for record_id in ("SUBA", "SUBB")),
        encoding="ascii",
    )
    fasta_reads: list[str] = []
    real_reader = genome_module._read_fasta_records

    def counting_reader(path):
        fasta_reads.append(str(path))
        return real_reader(path)

    monkeypatch.setattr(genome_module, "_read_fasta_records", counting_reader)
    plan_circular_request(_ring_request(tmp_path, RecordInput(source=GffFastaInputSource(gff, fasta))))
    assert [subject for _query, subject in captured_ring_searches] == [_expected_database("SUBA", "SUBB")]
    # The GFF3 parse keeps only the entries the GFF3 names, so the database reads the FASTA file a second time.
    assert fasta_reads == [str(fasta), str(fasta)]


@pytest.mark.circular
def test_ring_database_of_a_file_with_a_byte_order_mark_is_what_the_loader_reads(
    tmp_path: Path, captured_ring_searches: list[tuple[str, str]]
) -> None:
    from gbdraw.api.request_render import plan_circular_request
    from gbdraw.io.genome import load_gbks
    from gbdraw.io.record_select import parse_record_selector

    path = _subject_file(tmp_path, con=False, bom=True)
    loaded = [str(record.id) for record in load_gbks([str(path)])]
    plan_circular_request(_ring_request(
        tmp_path, RecordInput(source=GenBankInputSource(path=str(path)), selector=parse_record_selector("#1")),
    ))
    assert [subject for _query, subject in captured_ring_searches] == [_expected_database(*loaded)]


@pytest.mark.circular
def test_ring_database_of_several_sources_is_each_file_whole_in_request_order(
    tmp_path: Path, captured_ring_searches: list[tuple[str, str]]
) -> None:
    from gbdraw.api.request_render import plan_circular_batch_request
    from gbdraw.api.requests import CircularBatchRequest
    from gbdraw.io.record_select import parse_record_selector

    first = tmp_path / "first.gb"
    first.write_text(_genbank_text("SUBD", _SUBJECT_SEQUENCES["SUBD"]), encoding="ascii")
    second = _subject_file(tmp_path, con=False)
    single = _ring_request(
        tmp_path,
        RecordInput(source=GenBankInputSource(path=str(first))),
        RecordInput(source=GenBankInputSource(path=str(second)), selector=parse_record_selector("#2")),
        RecordInput(source=GenBankInputSource(path=str(first))),
    )
    # The first file is named twice and searched once.
    expected = _expected_database("SUBD", "SUBA", "SUBB", "SUBC")
    from gbdraw.api.request_render import plan_circular_request

    plan_circular_request(single)
    batch = CircularBatchRequest(
        records=single.records, options=single.options,
        outputs=tuple(RenderOutputRequest(output_prefix=f"item{index}") for index in range(3)),
    )
    plan_circular_batch_request(batch)
    assert [subject for _query, subject in captured_ring_searches] == [expected, expected]


@pytest.mark.circular
def test_ring_database_of_in_memory_records_is_the_records_given(
    tmp_path: Path, captured_ring_searches: list[tuple[str, str]]
) -> None:
    """The Python API passes records, not files: each record is its own source, as given."""

    from gbdraw.api.request_render import plan_circular_request
    from gbdraw.api.requests import InMemoryRecordSource

    records = read_genbank(_subject_file(tmp_path, con=False))
    plan_circular_request(_ring_request(
        tmp_path, *(RecordInput(source=InMemoryRecordSource(record)) for record in (records[2], records[0])),
    ))
    assert [subject for _query, subject in captured_ring_searches] == [_expected_database("SUBC", "SUBA")]


def test_ring_database_without_a_sequence_stops_the_search() -> None:
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    from gbdraw.comparisons.circular_losat import ring_subject_fasta

    with pytest.raises(ValidationError, match="no record with a sequence") as excinfo:
        ring_subject_fasta([SeqRecord(Seq(None, length=20), id="SUBCON")])
    assert _diagnostic(excinfo)["reason"] == "SEQUENCE_MISSING"


# Ring tables are in the search frame of the drawn record (PD-OI-073, Owner D-28).

# Whole-file search rows: SUBB 21-40, SUBB 45-60 (crosses the end of a crop to 11-50),
# SUBB 70-61 (outside it), SUBB 15-5 (reverse, crosses its start), and SUBA 1-20.
_WHOLE_FILE_RING_ROWS = "".join(
    f"QUERY1\t{subject}\t99.0\t20\t0\t0\t1\t20\t{sstart}\t{send}\t1e-30\t100\n"
    for subject, sstart, send in (("SUBB", 21, 40), ("SUBB", 45, 60), ("SUBB", 70, 61), ("SUBB", 15, 5), ("SUBA", 1, 20))
)


def _rows_with_subject_spans(rows: str, spans: list[tuple[int, int] | None]) -> str:
    """``rows`` with each row's subject span replaced, or left out where the span is None."""

    kept = []
    for line, span in zip(rows.splitlines(), spans):
        if span is not None:
            columns = line.split("\t")
            columns[8], columns[9] = (str(value) for value in span)
            kept.append("\t".join(columns) + "\n")
    return "".join(kept)


class _AllHits:
    evalue = 1.0
    bitscore = 0.0
    identity = 0.0
    alignment_length = 0


def _ring_hit_spans(plan) -> list[tuple[float, float]]:
    from gbdraw.analysis.conservation import (
        load_conservation_sources,
        normalize_conservation_tracks_for_record,
        read_conservation_inputs,
    )

    loaded = load_conservation_sources(
        blast_config=_AllHits(),
        conservation_files=read_conservation_inputs(plan.request.options.conservation_search_results),
        labels=["ring"],
    )
    (track,) = normalize_conservation_tracks_for_record(
        loaded, displayed_records=plan.records, record=plan.records[0], conservation_reference="subject",
    )
    return sorted((row.draw_start, row.draw_end) for row in track.hits.itertuples())


@pytest.fixture
def whole_file_ring_search(monkeypatch: pytest.MonkeyPatch) -> list[str]:
    """Each ring search returns ``_WHOLE_FILE_RING_ROWS``; the list holds the subject FASTA of each search."""

    import gbdraw.comparisons.circular_losat as circular_losat_module

    subjects: list[str] = []

    def fake_search(_spec, _query_fasta, subject_fasta, **_kwargs):
        subjects.append(subject_fasta)
        return _WHOLE_FILE_RING_ROWS

    monkeypatch.setattr(circular_losat_module, "run_losat_search", fake_search)
    return subjects


@pytest.mark.circular
@pytest.mark.parametrize("reverse", [False, True], ids=["crop", "reverse complement crop"])
def test_ring_rows_of_a_crop_are_a_table_of_the_cropped_record(
    reverse: bool, tmp_path: Path, whole_file_ring_search: list[str]
) -> None:
    """The ring rows are what a search of the cropped record reports; the raw cache keeps the search output."""

    from gbdraw.api.request_render import plan_circular_request
    from gbdraw.api.requests import RecordPresentation
    from gbdraw.io.regions import parse_region_spec

    source = GenBankInputSource(path=str(_subject_file(tmp_path)))
    plan = plan_circular_request(_ring_request(tmp_path, RecordInput(
        source=source, region=parse_region_spec("SUBB:11-50"),
        presentation=RecordPresentation(reverse_complement=reverse),
    )))

    (ring,) = plan.request.options.conservation_search_results
    # SUBB 21-40 is 11-30 of the 40 bp crop; 45-60 is clipped to 35-40; 15-5 to 5-1; 70-61 and SUBA are left out.
    assert ring.text == _rows_with_subject_spans(_WHOLE_FILE_RING_ROWS, [(11, 30), (35, 40), None, (5, 1), None])
    (entry,) = plan.losat_cache_entries
    assert entry["text"] == _WHOLE_FILE_RING_ROWS
    expected = [(0.0, 5.0), (10.0, 30.0), (34.0, 40.0)] if not reverse else [(0.0, 6.0), (10.0, 30.0), (35.0, 40.0)]
    assert _ring_hit_spans(plan) == expected


@pytest.mark.circular
def test_ring_rows_of_an_already_cropped_in_memory_record_planned_again(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, whole_file_ring_search: list[str]
) -> None:
    """A materialized crop is searched as given; its rows land where a search of the file puts them (review N3)."""

    import gbdraw.comparisons.circular_losat as circular_losat_module
    from gbdraw.api.request_render import plan_circular_request
    from gbdraw.io.regions import parse_region_spec

    first = plan_circular_request(_ring_request(tmp_path, RecordInput(
        source=GenBankInputSource(path=str(_subject_file(tmp_path))), region=parse_region_spec("SUBB:11-50"),
    )))
    # The database is the 40 bp crop as given: the search reports SUBB 21-40 of the file as 11-30.
    crop_rows = _rows_with_subject_spans(_WHOLE_FILE_RING_ROWS, [(11, 30), None, None, None, None])
    monkeypatch.setattr(circular_losat_module, "run_losat_search", lambda *_args, **_kwargs: crop_rows)
    ring_options = _ring_request(tmp_path, *first.request.records).options
    again = plan_circular_request(replace(first.request, options=ring_options))

    assert again.request.options.conservation_search_results[0].text == crop_rows
    assert _ring_hit_spans(again) == [(10.0, 30.0)]


def test_ring_rows_conversion_maps_in_memory_database_records_through_their_coordinate_map() -> None:
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    from gbdraw.comparisons.circular_losat import ring_rows_in_search_frame

    drawn = SeqRecord(Seq("A" * 40), id="SUBB", annotations={"gbdraw_coord_base": 50, "gbdraw_coord_step": -1})
    # The database record is the same crop reverse complemented: its 1-20 is source 50-31.
    database = SeqRecord(Seq("A" * 40), id="SUBB", annotations={"gbdraw_coord_base": 50, "gbdraw_coord_step": -1})
    rows = "QUERY1\tSUBB\t99.0\t20\t0\t0\t1\t20\t1\t20\t1e-30\t100\n"

    converted = ring_rows_in_search_frame(rows, drawn_records=[drawn], subject_records=[database])

    # Source 50-31 is 40-21 of the crop on its source strand.
    assert converted == "QUERY1\tSUBB\t99.0\t20\t0\t0\t1\t20\t40\t21\t1e-30\t100\n"


def test_ring_rows_conversion_stops_when_one_id_has_two_frames() -> None:
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    from gbdraw.comparisons.circular_losat import ring_rows_in_search_frame

    crops = [
        SeqRecord(Seq("A" * 40), id="SUBB", annotations={"gbdraw_coord_base": base, "gbdraw_coord_step": 1})
        for base in (11, 31)
    ]
    with pytest.raises(ValidationError, match="two records with this ID") as excinfo:
        ring_rows_in_search_frame(_WHOLE_FILE_RING_ROWS, drawn_records=crops)
    assert _diagnostic(excinfo) == {"code": "COMPARISON_INPUT", "reason": "UNIQUE_IDS", "field": "records"}


@pytest.mark.circular
def test_cli_session_with_a_crop_and_a_search_frame_ring_table_replays_its_hits(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """A CLI Session with crop-frame ring rows draws them on the crop, as before #1314 (PD-OI-073).

    The Session is written by this branch from a precomputed table; 0.14.0 CLI Sessions have the same format
    (SESSION_VERSION 48) and the same crop-frame rows.
    """

    records = tmp_path / "records.tsv"
    records.write_text(f"gbk\tregion\n{HUMAN}\t501-800\n", encoding="utf-8")
    table = tmp_path / "ring.tsv"
    table.write_text("query1\tNC_012920.1\t99.0\t100\t0\t0\t1\t100\t101\t200\t1e-50\t500\n", encoding="utf-8")
    prefix = tmp_path / "saved"
    circular_main([
        "--records_table", str(records), "--conservation_blast", str(table),
        "--conservation_reference", "subject", "--save_session", "-o", str(prefix), "-f", "svg",
    ])
    hit = re.compile(r'<path d="([^"]+)"[^>]*data-sstart="([0-9]+)"')
    (saved,) = hit.findall(prefix.with_suffix(".svg").read_text(encoding="utf-8"))
    assert saved[1] == "101"

    def no_runtime(*_args, **_kwargs):
        raise AssertionError("Session replay resolved a LOSAT runtime.")

    monkeypatch.setattr(runtime_module, "resolve_losat_runtime", no_runtime)
    replay = tmp_path / "replay"
    circular_main(["--session", str(prefix.with_suffix(".gbdraw-session.json")), "-o", str(replay), "-f", "svg"])
    assert hit.findall(replay.with_suffix(".svg").read_text(encoding="utf-8")) == [saved]


@pytest.mark.circular
def test_interactive_ring_popup_reference_span_of_a_crop_slices_the_hit_bases(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """The popup slices the embedded cropped record at the row's span: the rows are crop-frame (review B3).

    The slice models the popup runtime's JavaScript slice in Python.
    """

    import gbdraw.comparisons.circular_losat as circular_losat_module
    from Bio import SeqIO

    sequence = str(next(SeqIO.parse(str(HUMAN), "genbank")).seq).upper()
    query = tmp_path / "query.fasta"
    query.write_text(f">QUERY1\n{sequence[620:700]}\n", encoding="ascii")
    monkeypatch.setattr(
        circular_losat_module, "run_losat_search",
        lambda *_args, **_kwargs: "QUERY1\tNC_012920.1\t100.0\t80\t0\t0\t1\t80\t621\t700\t1e-40\t300\n",
    )
    records = tmp_path / "records.tsv"
    records.write_text(f"gbk\tregion\n{HUMAN}\t501-800\n", encoding="utf-8")
    prefix = tmp_path / "crop"
    circular_main([
        "--records_table", str(records), "--losat", "losatn", "--losat_bin", "/bin/true",
        "--conservation_sequence", str(query), "-o", str(prefix), "-f", "interactive_svg",
    ])
    svg = prefix.with_suffix(".interactive.svg").read_text(encoding="utf-8")

    ((start, end),) = {(int(s), int(e)) for s, e in re.findall(r'"sstart":"(\d+)","send":"(\d+)"', svg)}
    reference = re.findall(r'"sequence":"([ACGTN]+)"', svg)[0]
    assert len(reference) == 300
    # The interactive runtime slices the embedded reference at sstart..send (buildEmbeddedMatchEntry).
    assert reference[start - 1:end] == sequence[620:700]
