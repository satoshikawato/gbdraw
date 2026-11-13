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
from gbdraw.session import SessionDocument, build_session_document
from tests.utils.session_files import embedded_session
from gbdraw.web_support.request_render import (
    _render_staged_canonical_web_request_with_prepared_inputs,
    render_embedded_canonical_web_request,
)

ROOT = Path(__file__).resolve().parents[1]
VARIANTS = ROOT / "tests" / "fixtures" / "input_variants"
ZERO_FEATURES = VARIANTS / "zero_features_record.gb"
QUERY = VARIANTS / "prokka_locus_tag_only.gb"
WITH_REPEATS = ROOT / "tests" / "test_inputs" / "NC_001879.gbk"
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
    """The CLI's log lines: every level goes to stderr, and stdout stays empty."""
    assert run.stdout == ""
    return run.stderr.splitlines()


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_cli_warns_about_a_record_with_no_feature_and_still_draws(tmp_path: Path, mode: str) -> None:
    run = _gbdraw(mode, "--gbk", str(ZERO_FEATURES), "-o", "out", "-f", "svg", cwd=tmp_path)

    assert run.returncode == 0, run.stderr
    lines = _lines(run)
    assert "INFO: Record EMPTY001.1 (zero_features_record.gb), record 1 of 1, 1,000 bp: drew 0 features." in lines
    assert ZERO_FEATURES_WARNING in lines
    assert (tmp_path / "out.svg").is_file()


def test_one_circle_summary_says_a_record_left_off_the_circle_is_not_drawn(tmp_path: Path) -> None:
    """A record shorter than the one-circle minimum is not drawn: its line says so, and it gets no
    no-feature warning (OV-3160)."""
    run = _gbdraw(
        "circular", "--gbk", str(QUERY), str(ZERO_FEATURES), "--one_circle", "--one_circle_min_length", "2000",
        "-o", "out", "-f", "svg", cwd=tmp_path,
    )

    assert run.returncode == 0, run.stderr
    lines = _lines(run)
    assert (
        "INFO: Record EMPTY001.1 (zero_features_record.gb), record 1 of 1, 1,000 bp: not drawn, "
        "shorter than 2,000 bp, the shortest record drawn on one circle."
    ) in lines
    assert ZERO_FEATURES_WARNING not in lines
    assert any(line.startswith("INFO: Record PROKKA01.1 (prokka_locus_tag_only.gb), record 1 of 1, 3,000 bp: drew ")
               for line in lines)


def test_cli_names_the_types_a_record_has_when_none_is_drawn(tmp_path: Path) -> None:
    run = _gbdraw("circular", "--gbk", str(QUERY), "-k", "misc_RNA", "-o", "out", "-f", "svg", cwd=tmp_path)

    assert run.returncode == 0, run.stderr
    assert (
        "WARNING: Record PROKKA01.1 (prokka_locus_tag_only.gb) has no feature to draw: its 2 CDS, 1 rRNA, "
        "1 source, 1 tRNA features are of types not drawn; add a type to the feature types to draw "
        "(-k/--features) to draw it."
    ) in _lines(run)


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_cli_says_a_track_stack_without_a_features_slot_draws_no_feature(tmp_path: Path, mode: str) -> None:
    """The stack draws no feature, so the line claims none and no input warning is raised (OV-2797)."""
    slot = f"--{mode}_track_slot"
    run = _gbdraw(
        mode, "--gbk", str(QUERY), slot, "gc:gc_content", slot, "skew:gc_skew", "-o", "out", "-f", "svg", cwd=tmp_path
    )

    assert run.returncode == 0, run.stderr
    lines = _lines(run)
    assert (
        "INFO: Record PROKKA01.1 (prokka_locus_tag_only.gb), record 1 of 1, 3,000 bp: "
        "drew no features (no enabled features track slot)."
    ) in lines
    assert not any(line.startswith("WARNING:") for line in lines), lines


REGION_SET_TABLE = "set_id\tid\tmark\tstart\tend\tlabel\nislands\tisland1\tbracket\t500\t1500\tIsland 1\n"
REGION_SET_LEFT_OUT_WARNING = (
    "WARNING: Region annotation set islands is not drawn: no annotations track slot names it. "
    "Add an annotations track slot for it (set_id=islands) to draw it."
)


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_cli_warns_once_about_a_region_set_the_track_stack_leaves_out(tmp_path: Path, mode: str) -> None:
    """A custom stack draws only the sets its annotations slots name; the run says which set it leaves out and how
    to draw it, once (OV-3007)."""
    table = tmp_path / "regions.tsv"
    table.write_text(REGION_SET_TABLE, encoding="utf-8")
    slot = f"--{mode}_track_slot"
    run = _gbdraw(
        mode, "--gbk", str(QUERY), "--annotation_table", str(table), slot, "features:features",
        "-o", "out", "-f", "svg", cwd=tmp_path,
    )

    assert run.returncode == 0, run.stderr
    assert [line for line in _lines(run) if line.startswith("WARNING:")] == [REGION_SET_LEFT_OUT_WARNING]
    assert (tmp_path / "out.svg").is_file()


@pytest.mark.parametrize("mode", ["circular", "linear"])
@pytest.mark.parametrize("slots", [(), ("features:features", "islands:annotations@set_id=islands")])
def test_cli_keeps_quiet_when_every_region_set_is_drawn(tmp_path: Path, mode: str, slots: tuple[str, ...]) -> None:
    """Without a stack every set gets its automatic slot; a stack that names the set draws it (OV-3007)."""
    table = tmp_path / "regions.tsv"
    table.write_text(REGION_SET_TABLE, encoding="utf-8")
    extra = [arg for spec in slots for arg in (f"--{mode}_track_slot", spec)]
    run = _gbdraw(
        mode, "--gbk", str(QUERY), "--annotation_table", str(table), *extra, "-o", "out", "-f", "svg", cwd=tmp_path,
    )

    assert run.returncode == 0, run.stderr
    assert not any(line.startswith("WARNING:") for line in _lines(run)), run.stderr


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_cli_draws_a_stack_without_a_features_slot_for_a_record_with_underlays(tmp_path: Path, mode: str) -> None:
    """repeat_region is drawn as an underlay by default; a stack without a features slot draws it no more than the
    other features, and the run succeeds (OV-2801)."""
    slot = f"--{mode}_track_slot"
    run = _gbdraw(mode, "--gbk", str(WITH_REPEATS), slot, "gc:gc_content", "-o", "out", "-f", "svg", cwd=tmp_path)

    assert run.returncode == 0, run.stderr
    lines = _lines(run)
    assert (
        "INFO: Record NC_001879.2 (NC_001879.gbk), record 1 of 1, 155,943 bp: "
        "drew no features (no enabled features track slot)."
    ) in lines
    assert not any(line.startswith("WARNING:") for line in lines), lines
    svg = (tmp_path / "out.svg").read_text()
    assert 'data-gbdraw-auto-feature-underlay="true"' not in svg
    assert 'data-legend-key="repeat_region"' not in svg


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


def _embedded(document: SessionDocument) -> dict:
    return embedded_session(document.to_dict(), document.resources.sources)


def test_web_reply_carries_the_lines_and_warnings_the_cli_logs(tmp_path: Path) -> None:
    """Parity: the Web shows the summary the CLI logs, from the same owner."""

    request = LinearDiagramRequest(
        records=tuple(RecordInput(source=GenBankInputSource(path)) for path in (QUERY, SUBJECT)),
        options=LinearDiagramOptions(blast_files=(str(ROWS_BELOW_IDENTITY),), identity=70.0),
        output=RenderOutputRequest(output_prefix="summary", formats=("svg",), output_directory=tmp_path / "cli"),
    )
    rendered = render_request(request)

    document = _embedded(build_session_document(request))
    reply = render_embedded_canonical_web_request(
        document["drawings"][0]["renderRequest"], resources=document["resources"], workspace=tmp_path / "web",
    )
    summary = rendered.input_summary
    assert reply["metadata"]["inputSummaries"] == [{
        "resultIndex": 0,
        "resultName": "summary.svg",
        "lines": list(input_summary_lines(summary)),
        "warnings": [{"code": "comparison_rows_filtered_out", "message": FILTERED_OUT_WARNING.removeprefix("WARNING: ")}],
        "records": [{"sourceLength": record.source_length, "region": None}
                    for record in summary.records],
        "linearWrap": None,
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
        document["drawings"][0]["renderRequest"], resource_paths=paths, workspace=workspace, _resource_names=names,
    )
    assert staged["metadata"]["inputSummaries"] == reply["metadata"]["inputSummaries"]


def test_web_reply_warns_about_a_record_with_no_feature(tmp_path: Path) -> None:
    from gbdraw.api import CircularDiagramRequest

    document = _embedded(build_session_document(CircularDiagramRequest(
        records=(RecordInput(source=GenBankInputSource(ZERO_FEATURES)),),
        output=RenderOutputRequest(output_prefix="empty", formats=("svg",)),
    )))
    reply = render_embedded_canonical_web_request(
        document["drawings"][0]["renderRequest"], resources=document["resources"], workspace=tmp_path / "web",
    )

    assert reply["metadata"]["inputSummaries"][0]["warnings"] == [
        {"code": "no_drawable_features", "message": ZERO_FEATURES_WARNING.removeprefix("WARNING: ")}
    ]



@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_python_api_logs_a_region_set_the_track_stack_leaves_out(
    tmp_path: Path, mode: str, caplog: pytest.LogCaptureFixture
) -> None:
    """``draw_circular`` and ``draw_linear`` log the input summary's warnings, as the CLI does (OV-3007)."""
    from gbdraw import CircularOptions, CircularTrackOptions, LinearOptions, LinearTrackOptions, read_genbank
    from gbdraw.interface import AnnotationOptions

    table = tmp_path / "regions.tsv"
    table.write_text(REGION_SET_TABLE, encoding="utf-8")
    records = read_genbank(str(QUERY))
    annotations = AnnotationOptions(table_file=str(table))
    with caplog.at_level("WARNING", logger="gbdraw.api.input_summary"):
        if mode == "circular":
            gbdraw.draw_circular(records, options=CircularOptions(
                annotations=annotations, tracks=CircularTrackOptions(slots=("features:features",))))
        else:
            gbdraw.draw_linear(records, options=LinearOptions(
                annotations=annotations, tracks=LinearTrackOptions(slots=("features:features",))))

    assert [record.getMessage() for record in caplog.records if record.levelname == "WARNING"] == [
        REGION_SET_LEFT_OUT_WARNING.removeprefix("WARNING: ")
    ]

@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_web_reply_names_a_region_set_the_track_stack_leaves_out(tmp_path: Path, mode: str) -> None:
    """The Web shows the warning the CLI logs, from the same owner (OV-3007)."""
    from gbdraw.api import (
        AnnotationOptions,
        CircularDiagramOptions,
        CircularDiagramRequest,
        CircularRequestTrackOptions,
        LinearRequestTrackOptions,
    )

    table = tmp_path / "regions.tsv"
    table.write_text(REGION_SET_TABLE, encoding="utf-8")
    annotations = AnnotationOptions(table_file=str(table))
    records = (RecordInput(source=GenBankInputSource(QUERY)),)
    output = RenderOutputRequest(output_prefix="stack", formats=("svg",))
    request = (
        CircularDiagramRequest(
            records=records, output=output,
            options=CircularDiagramOptions(
                annotations=annotations, tracks=CircularRequestTrackOptions(circular_track_slots=("features:features",)),
            ),
        )
        if mode == "circular"
        else LinearDiagramRequest(
            records=records, output=output,
            options=LinearDiagramOptions(
                annotations=annotations, tracks=LinearRequestTrackOptions(linear_track_slots=("features:features",)),
            ),
        )
    )
    document = _embedded(build_session_document(request))
    reply = render_embedded_canonical_web_request(
        document["drawings"][0]["renderRequest"], resources=document["resources"], workspace=tmp_path / "web",
    )

    assert reply["metadata"]["inputSummaries"][0]["warnings"] == [
        {"code": "annotation_set_not_drawn", "message": REGION_SET_LEFT_OUT_WARNING.removeprefix("WARNING: ")}
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


@pytest.mark.parametrize("mode", ["pairwise", "similarity_groups"])
def test_web_reply_words_a_stored_protein_pair_without_a_match_as_the_cli(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path, mode: str,
) -> None:
    """OV-3512: after a Generate the Web sends each LOSATP pair as a stored
    protein comparison table; a pair without a match is a protein comparison
    with no protein match to draw, as the CLI's fresh render says."""

    from dataclasses import replace

    import gbdraw.analysis.protein_colinearity as protein_colinearity_module
    from gbdraw.api.request_render import stored_protein_result_request
    from tests.test_render_stage_protein_analysis import _no_hit_losatp, _no_search, _request

    monkeypatch.setattr(protein_colinearity_module, "run_losatp_blastp", _no_hit_losatp)
    request, _identities = _request(mode, tmp_path)
    request = replace(request, output=RenderOutputRequest(
        output_prefix="protein", formats=("svg",), output_directory=tmp_path / "out",
    ))
    rendered = render_request(request)
    fresh = [item.message for item in input_warnings(rendered.input_summary)]
    assert fresh and all(message.endswith("has no protein match to draw.") for message in fresh), fresh

    stored = stored_protein_result_request(rendered.request, rendered.linear_metadata)
    document = _embedded(build_session_document(stored))
    comparisons = document["drawings"][0]["renderRequest"]["comparisons"]
    assert "precomputedProteinComparison" in {comparison["kind"] for comparison in comparisons}
    monkeypatch.setattr(protein_colinearity_module, "run_losatp_blastp", _no_search)
    reply = render_embedded_canonical_web_request(
        document["drawings"][0]["renderRequest"], resources=document["resources"], workspace=tmp_path / "web",
    )
    (summary,) = reply["metadata"]["inputSummaries"]
    assert [warning["message"] for warning in summary["warnings"]] == fresh
    assert summary["lines"] == list(input_summary_lines(rendered.input_summary))


def test_a_cropped_record_with_nothing_to_draw_names_its_region() -> None:
    from gbdraw.api.input_summary import FeatureTypeTally, InputSummary, RecordInputSummary

    record = RecordInputSummary(
        record_index=0, record_id="ctg1", source_name="genes.gff3", source_paths=("genes.gff3",), source_record_index=0,
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


def test_a_window_across_the_origin_names_its_span_and_its_export_check_reads_it() -> None:
    # A window across the origin of a circular record: the summary names the source span the drawing draws, start
    # after end, and the wrapped export check reads the drawn bases of that span (670 bp), not an ordered pair.
    import json

    from gbdraw.api.input_summary import input_summary_payload, summarize_inputs
    from gbdraw.api.request_render import resolve_request_records
    from gbdraw.io.regions import RecordWindow
    from gbdraw.web_support.export_limit import check_web_export_limit_json

    request = LinearDiagramRequest(records=(
        RecordInput(GenBankInputSource(ROOT / "tests/test_inputs/HmmtDNA.gbk"), windows=(RecordWindow(16000, 100),),
                    record_key="mito"),
    ))
    resolved = resolve_request_records(request)
    summary = summarize_inputs(resolved.records, resolved.provenance, selected_features=["CDS"],
                               feature_visibility_rules=None, specific_color_rules=None)
    assert summary.records[0].region == (16000, 100)
    assert input_summary_lines(summary)[0].startswith(
        "Record NC_012920.1 (HmmtDNA.gbk), record 1 of 1, 16,000..100 across the origin of 16,569 bp:"
    )
    records = input_summary_payload(summary, result_index=0, result_name="x")["records"]
    reply = json.loads(check_web_export_limit_json(json.dumps(
        {"format": "png", "width": 1000.0, "height": 500.0, "dpi": 96, "linearWrap": {"requested": 500, "records": records}}
    )))
    assert reply["ok"] is True
