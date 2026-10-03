"""Circular similarity rings from LOSATN / TLOSATX (design 3.3, 3.4, 3.7, D12, PR-4)."""

from __future__ import annotations

from contextlib import ExitStack
import json
from pathlib import Path

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
    with ExitStack() as stack:
        try:
            runtime = resolve_losat_runtime("tlosatx", stack=stack)
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
    # Default labels: FASTA file name; GenBank DEFINITION.
    assert labels == ["NC_002333.2.fna", "Drosophila melanogaster mitochondrion, complete genome"]

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
