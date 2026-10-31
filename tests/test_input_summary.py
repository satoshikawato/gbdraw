"""The input summary: what a render read and drew, logged by the CLI and carried to the Web.

A record with no feature to draw and a comparison table whose rows the thresholds
all leave out draw nothing; the run still draws and says why (one owner,
``gbdraw.api.input_summary``). A UTF-8 byte order mark no longer hides a file's
first record (OV-2354).
"""

from __future__ import annotations

import base64
import subprocess
import sys
from pathlib import Path

import pytest

import gbdraw
from gbdraw.api import (
    GenBankInputSource,
    LinearDiagramOptions,
    LinearDiagramRequest,
    RecordInput,
    RenderOutputRequest,
    render_request,
)
from gbdraw.api.input_summary import input_summary_lines, input_warnings
from gbdraw.io.genome import parse_gff_fasta, read_flat_file_records
from gbdraw.session import build_session_document
from gbdraw.web_support.request_render import (
    _render_staged_canonical_web_request_with_prepared_inputs,
    render_embedded_canonical_web_request,
)

ROOT = Path(__file__).resolve().parents[1]
VARIANTS = ROOT / "tests" / "fixtures" / "input_variants"
ZERO_FEATURES = VARIANTS / "zero_features_record.gb"
QUERY = VARIANTS / "prokka_locus_tag_only.gb"
SUBJECT = VARIANTS / "rast_seed_dbxref.gb"
ROWS_BELOW_IDENTITY = VARIANTS / "blast_rows_below_identity.tsv"
BOM = b"\xef\xbb\xbf"

ZERO_FEATURES_WARNING = (
    "WARNING: Record EMPTY001.1 (zero_features_record.gb) has no feature to draw: it has no features."
)
FILTERED_OUT_WARNING = (
    "WARNING: Comparison table blast_rows_below_identity.tsv (record #1 PROKKA01.1 and record #2 RAST0001.1): "
    "the thresholds leave out all 5 rows, so no match is drawn: 5 rows below identity 70% (best 65.2%)."
)


def _gbdraw(*args: str, cwd: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, "-c", "import sys; from gbdraw.cli import main; sys.argv[0] = 'gbdraw'; sys.exit(main())", *args],
        cwd=cwd, capture_output=True, text=True, env={"PYTHONPATH": str(ROOT), "PATH": ""}, check=False,
    )


def _lines(run: subprocess.CompletedProcess[str]) -> list[str]:
    """The CLI's log lines (INFO and WARNING go to stdout, ERROR to stderr)."""
    return (run.stdout + run.stderr).splitlines()


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_cli_warns_about_a_record_with_no_feature_and_still_draws(tmp_path: Path, mode: str) -> None:
    run = _gbdraw(mode, "--gbk", str(ZERO_FEATURES), "-o", "out", "-f", "svg", cwd=tmp_path)

    assert run.returncode == 0, run.stderr
    lines = _lines(run)
    assert "INFO: Record EMPTY001.1 (zero_features_record.gb), record 1 of 1, 1,000 bp: drew 0 features." in lines
    assert ZERO_FEATURES_WARNING in lines
    assert (tmp_path / "out.svg").is_file()


def test_cli_names_the_types_a_record_has_when_none_is_drawn(tmp_path: Path) -> None:
    run = _gbdraw("circular", "--gbk", str(QUERY), "-k", "misc_RNA", "-o", "out", "-f", "svg", cwd=tmp_path)

    assert run.returncode == 0, run.stderr
    assert (
        "WARNING: Record PROKKA01.1 (prokka_locus_tag_only.gb) has no feature to draw: its 2 CDS, 1 rRNA, "
        "1 source, 1 tRNA features are of types not drawn; add a type to the feature types to draw "
        "(-k/--features) to draw it."
    ) in _lines(run)


def test_cli_warns_when_the_thresholds_leave_out_every_comparison_row(tmp_path: Path) -> None:
    run = _gbdraw(
        "linear", "--gbk", str(QUERY), str(SUBJECT), "-b", str(ROWS_BELOW_IDENTITY), "--identity", "70",
        "-o", "out", "-f", "svg", cwd=tmp_path,
    )

    assert run.returncode == 0, run.stderr
    lines = _lines(run)
    assert (
        "INFO: Comparison table blast_rows_below_identity.tsv (record #1 PROKKA01.1 and record #2 RAST0001.1): "
        "kept 0 of 5 rows; 5 rows below identity 70% (best 65.2%)."
    ) in lines
    assert FILTERED_OUT_WARNING in lines
    assert (tmp_path / "out.svg").is_file()


def test_cli_reports_each_threshold_and_keeps_quiet_when_rows_remain(tmp_path: Path) -> None:
    run = _gbdraw(
        "linear", "--gbk", str(QUERY), str(SUBJECT), "-b", str(ROWS_BELOW_IDENTITY), "--identity", "60",
        "--bitscore", "100", "-o", "out", "-f", "svg", cwd=tmp_path,
    )

    assert run.returncode == 0, run.stderr
    assert (
        "INFO: Comparison table blast_rows_below_identity.tsv (record #1 PROKKA01.1 and record #2 RAST0001.1): "
        "kept 1 of 5 rows; 4 rows below bitscore 100 (best 90); 3 rows below identity 60% (best 58.4%)."
    ) in _lines(run)
    assert not any(line.startswith("WARNING: Comparison table") for line in _lines(run))


def test_cli_stops_with_an_error_when_a_file_has_no_record(tmp_path: Path) -> None:
    empty = tmp_path / "empty.gb"
    empty.write_text("")
    run = _gbdraw("circular", "--gbk", str(empty), "-o", "out", "-f", "svg", cwd=tmp_path)

    assert run.returncode == 1
    assert f"ERROR: {empty} is empty: it has no GenBank, DDBJ or EMBL flat file record." in _lines(run)
    assert not (tmp_path / "out.svg").exists()


def test_gff3_summary_counts_the_types_the_type_filter_left_out(tmp_path: Path) -> None:
    run = _gbdraw(
        "circular", "--gff", str(VARIANTS / "ncbi_gene_mrna_cds.gff3"), "--fasta", str(VARIANTS / "ncbi_gene_mrna_cds.fasta"),
        "-o", "out", "-f", "svg", cwd=tmp_path,
    )

    assert run.returncode == 0, run.stderr
    assert (
        "INFO: Record ctg1 (ncbi_gene_mrna_cds.gff3), record 1 of 1, 3,000 bp: drew 2 features (2 CDS); "
        "left out 3 exon, 2 gene, 1 mRNA, 1 region (types not drawn)."
    ) in _lines(run)


def test_cli_summary_names_a_cropped_region(tmp_path: Path) -> None:
    run = _gbdraw("linear", "--gbk", str(QUERY), "--region", "PROKKA01.1:1001-2000", "-o", "out", "-f", "svg", cwd=tmp_path)

    assert run.returncode == 0, run.stderr
    assert (
        "INFO: Record PROKKA01.1 (prokka_locus_tag_only.gb), record 1 of 1, 1,001..2,000 of 3,000 bp: "
        "drew 3 features (1 CDS, 1 rRNA, 1 tRNA); left out 1 source (type not drawn)."
    ) in _lines(run)


def test_web_reply_carries_the_lines_and_warnings_the_cli_logs(tmp_path: Path) -> None:
    """Parity: the Web shows the summary the CLI logs, from the same owner."""

    request = LinearDiagramRequest(
        records=tuple(RecordInput(source=GenBankInputSource(path)) for path in (QUERY, SUBJECT)),
        options=LinearDiagramOptions(blast_files=(str(ROWS_BELOW_IDENTITY),), identity=70.0),
        output=RenderOutputRequest(output_prefix="summary", formats=("svg",), output_directory=tmp_path / "cli"),
    )
    rendered = render_request(request)

    document = build_session_document(request).to_dict()
    reply = render_embedded_canonical_web_request(
        document["renderRequest"], resources=document["resources"], workspace=tmp_path / "web",
    )
    summary = rendered.input_summary
    assert reply["metadata"]["inputSummaries"] == [{
        "resultIndex": 0,
        "resultName": "summary.svg",
        "lines": list(input_summary_lines(summary)),
        "warnings": [{"code": "comparison_rows_filtered_out", "message": FILTERED_OUT_WARNING.removeprefix("WARNING: ")}],
    }]
    assert [item.code for item in input_warnings(summary)] == ["comparison_rows_filtered_out"]

    # The Worker stages each resource under its content token; the summary names the file the page read.
    workspace = tmp_path / "gbdraw-web-render-1"
    (workspace / "resources").mkdir(parents=True)
    (workspace / ".gbdraw-worker-render-workspace").touch()
    cache = tmp_path / "gbdraw-web-render-resource-cache"
    cache.mkdir()
    paths, names = {}, {}
    for index, (resource_id, entry) in enumerate(document["resources"].items(), start=1):
        token = cache / f"{index:064x}.bin"
        token.write_bytes(base64.b64decode(entry["data"], validate=True))
        paths[resource_id] = workspace / "resources" / f"{index:04d}.bin"
        paths[resource_id].symlink_to(token)
        names[resource_id] = entry["name"]
    staged = _render_staged_canonical_web_request_with_prepared_inputs(
        document["renderRequest"], resource_paths=paths, workspace=workspace, _resource_names=names,
    )
    assert staged["metadata"]["inputSummaries"] == reply["metadata"]["inputSummaries"]


def test_web_reply_warns_about_a_record_with_no_feature(tmp_path: Path) -> None:
    from gbdraw.api import CircularDiagramRequest

    document = build_session_document(CircularDiagramRequest(
        records=(RecordInput(source=GenBankInputSource(ZERO_FEATURES)),),
        output=RenderOutputRequest(output_prefix="empty", formats=("svg",)),
    )).to_dict()
    reply = render_embedded_canonical_web_request(
        document["renderRequest"], resources=document["resources"], workspace=tmp_path / "web",
    )

    assert reply["metadata"]["inputSummaries"][0]["warnings"] == [
        {"code": "no_drawable_features", "message": ZERO_FEATURES_WARNING.removeprefix("WARNING: ")}
    ]


# OV-2354: Biopython reads a leading U+FEFF as text before the first record and skips it.
def test_a_byte_order_mark_keeps_the_first_genbank_record(tmp_path: Path) -> None:
    text = (VARIANTS / "cgview7_multi_record.gbff").read_bytes()
    path = tmp_path / "bom.gbff"
    path.write_bytes(BOM + text)

    assert [record.id for record in read_flat_file_records(str(path))] == ["MULTI001.1", "MULTI002.1", "MULTI003.1"]
    assert [record.id for record in gbdraw.read_genbank(str(path))] == ["MULTI001.1", "MULTI002.1", "MULTI003.1"]


def test_a_byte_order_mark_keeps_an_embl_file_embl(tmp_path: Path) -> None:
    path = tmp_path / "bom.embl"
    path.write_bytes(BOM + (VARIANTS / "easyfig72_embl.embl").read_bytes())

    assert [record.id for record in read_flat_file_records(str(path))] == [
        record.id for record in read_flat_file_records(str(VARIANTS / "easyfig72_embl.embl"))
    ]


def test_a_byte_order_mark_keeps_the_first_fasta_entry_and_gff3_line(tmp_path: Path) -> None:
    gff = tmp_path / "bom.gff3"
    fasta = tmp_path / "bom.fasta"
    gff.write_bytes(BOM + (VARIANTS / "ncbi_gene_mrna_cds.gff3").read_bytes())
    fasta.write_bytes(BOM + (VARIANTS / "ncbi_gene_mrna_cds.fasta").read_bytes())

    plain = parse_gff_fasta(str(VARIANTS / "ncbi_gene_mrna_cds.gff3"), str(VARIANTS / "ncbi_gene_mrna_cds.fasta"))
    marked = parse_gff_fasta(str(gff), str(fasta))
    assert [(record.id, len(record), len(record.features)) for record in marked] == [
        (record.id, len(record), len(record.features)) for record in plain
    ]


def test_a_protein_comparison_counts_the_matches_drawn_and_no_threshold(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path,
) -> None:
    """LOSATP hits are filtered before they reach the drawing (OV-2389): the summary says so."""

    from dataclasses import replace

    import gbdraw.analysis.protein_colinearity as protein_colinearity_module
    from tests.test_render_stage_protein_analysis import _fake_losatp, _request

    monkeypatch.setattr(protein_colinearity_module, "run_losatp_blastp", _fake_losatp)
    request, _identities = _request("pairwise", tmp_path)
    request = replace(request, options=replace(request.options, identity=99.0), output=RenderOutputRequest(
        output_prefix="protein", formats=("svg",), output_directory=tmp_path / "out",
    ))
    summary = render_request(request).input_summary

    lines = [line for line in input_summary_lines(summary) if line.startswith("Protein comparison")]
    assert lines and all(line.endswith(": drew 0 protein matches.") for line in lines), lines
    assert {item.code for item in input_warnings(summary)} == {"comparison_table_empty"}
    assert all(item.message.endswith("has no protein match to draw.") for item in input_warnings(summary))


def test_a_cropped_record_with_nothing_to_draw_names_its_region() -> None:
    from gbdraw.api.input_summary import FeatureTypeTally, InputSummary, RecordInputSummary

    record = RecordInputSummary(
        record_index=0, record_id="ctg1", source_name="genes.gff3", source_record_index=0,
        source_record_count=1, drawn_length=101, source_length=3000, region=(2900, 3000),
        features=(FeatureTypeTally("exon", 0, 0, 0, 3),),
    )
    assert [item.message for item in input_warnings(InputSummary((record,)))] == [
        "Record ctg1 (genes.gff3) has no feature to draw: its region 2,900..3,000 holds no feature of a type drawn."
    ]
    # The GFF3 load counts the whole record, before the crop, so the line says so.
    assert input_summary_lines(InputSummary((record,))) == (
        "Record ctg1 (genes.gff3), record 1 of 1, 2,900..3,000 of 3,000 bp: drew 0 features; "
        "the whole record has 3 exon of type not drawn.",
    )


def test_a_byte_order_mark_on_a_line_of_its_own_is_not_text(tmp_path: Path) -> None:
    from gbdraw.io.genome import genbank_absence_diagnostic

    path = tmp_path / "marked.gb"
    path.write_bytes(BOM + b"\n>s1\nACGT\n")
    assert genbank_absence_diagnostic(str(path)) == {"code": "NO_RECORDS", "reason": "FASTA_IN_GENBANK"}
