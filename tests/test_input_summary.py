"""The input summary: what a render read and drew, logged by the CLI and carried to the Web.

A record with no feature to draw and a comparison table whose rows the thresholds
all leave out draw nothing; the run still draws and says why (one owner,
``gbdraw.api.input_summary``). A UTF-8 byte order mark no longer hides a file's
first record (OV-2354).
"""

from __future__ import annotations

import base64
import gzip
import json
import subprocess
import sys
from collections.abc import Sequence
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
from gbdraw.api.input_summary import InputSummary, input_summary_lines, input_warnings
from gbdraw.io.comparisons import ComparisonRowTally
from gbdraw.linear_comparison import ComparisonTableTally
from gbdraw.io.genome import parse_gff_fasta, read_flat_file_records
from gbdraw.session import SessionDocument, build_session_document
from tests.utils.session_files import embedded_session, load_session_mapping
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


_PAIR = "record #1 PROKKA01.1 and record #2 RAST0001.1"


def _staged_reply(payload: dict, resources: dict, folder: Path, names: dict[str, str] | None = None) -> dict:
    """The Worker's reply: each resource staged under a content token, named by the page's ``resourceManifest``
    name (``names``; else the Session's name, as after a Load)."""

    workspace = folder / "gbdraw-web-render-1"
    (workspace / "resources").mkdir(parents=True)
    (workspace / ".gbdraw-worker-render-workspace").touch()
    cache = folder / "gbdraw-web-render-resource-cache"
    cache.mkdir()
    paths, manifest = {}, {}
    for index, (resource_id, entry) in enumerate(resources.items(), start=1):
        token = cache / f"{index:064x}.bin"
        token.write_bytes(base64.b64decode(entry["data"], validate=True))
        paths[resource_id] = workspace / "resources" / f"{index:04d}.bin"
        paths[resource_id].symlink_to(token)
        manifest[resource_id] = (names or {}).get(resource_id, entry["name"])
    return _render_staged_canonical_web_request_with_prepared_inputs(
        payload, resource_paths=paths, workspace=workspace, _resource_names=manifest,
    )


def _comparison_names(texts: list[str]) -> list[str | None]:
    """The name each comparison line gives its table, or None for a comparison named by its records."""

    names: list[str | None] = []
    for text in texts:
        if text.startswith("Comparison table "):
            names.append(text.removeprefix("Comparison table ").partition(" (record #")[0])
        elif text.startswith("Comparison of "):
            names.append(None)
    return names


def _blast_request(folder: Path, files: list[tuple[str, bytes]]) -> LinearDiagramRequest:
    """QUERY and SUBJECT compared by each of ``files`` (name, bytes): one as a BLAST file, more as tables given
    without their files' text, as a Python caller gives them."""

    paths = []
    for index, (name, content) in enumerate(files, start=1):
        path = folder / "in" / str(index) / name
        path.parent.mkdir(parents=True)
        path.write_bytes(content)
        paths.append(path)
    return LinearDiagramRequest(
        records=tuple(RecordInput(source=GenBankInputSource(path)) for path in (QUERY, SUBJECT)),
        options=LinearDiagramOptions(
            linear_comparisons=tuple(
                _blast_comparison(path) for path in paths
            ) if len(paths) > 1 else None,
            blast_files=(str(paths[0]),) if len(paths) == 1 else None,
            identity=70.0,
        ),
        output=RenderOutputRequest(output_prefix="summary", formats=("svg",), output_directory=folder / "fresh"),
    )


def _blast_document(folder: Path, files: list[tuple[str, bytes]]) -> dict:
    return _embedded(build_session_document(_blast_request(folder, files)))


def _blast_comparison(path: Path) -> object:
    from gbdraw.io.comparisons import read_comparison_table
    from gbdraw.linear_comparison import LinearComparison

    return LinearComparison(0, 1, read_comparison_table(path), source_name=path.name)


def _older_entry(document: dict, resource_id: str, name: str, *, data: bytes | None = None) -> dict:
    """``document`` as a writer before the stored name wrote it: the comparison resource under ``resource_id``
    and ``name``, the entry without ``sourceName``."""

    payload = document["drawings"][0]["renderRequest"]
    (comparison,) = [item for item in payload["comparisons"] if item["kind"] == "nucleotideBlast"]
    entry = document["resources"].pop(comparison["resourceId"])
    if data is not None:
        entry = {**entry, "data": base64.b64encode(data).decode("ascii"), "size": len(data)}
    document["resources"][resource_id] = {**entry, "name": name}
    comparison["resourceId"] = resource_id
    comparison.pop("sourceName", None)
    return document


def _reencoded(document: dict, folder: Path) -> dict:
    """``document`` decoded by the Python API and written again (OV-3557 review M2)."""

    from gbdraw.session import materialize_session, session_to_request

    with materialize_session(load_session_mapping(document), output_directory=folder / "decoded") as materialized:
        request = session_to_request(materialized)
        return _embedded(build_session_document(request))


def _cli_replay_names(document: dict, folder: Path) -> list[str | None]:
    """The comparison names the CLI replay (`gbdraw render --session`) of ``document`` logs."""

    from gbdraw.session import materialize_session, render_session

    with materialize_session(load_session_mapping(document), output_directory=folder / "replay") as materialized:
        summary = render_session(materialized).input_summary
    return _comparison_names(list(input_summary_lines(summary)))


_ROWS = ROWS_BELOW_IDENTITY.read_bytes()
_HITS = [("hits.tsv", _ROWS)]
# Each case: the Session, the names the page's manifest gives its resources (None: the Session's names, as after a
# Load), and the names the CLI replay and the Web show.
_NAMING_CASES = {
    "plain": (lambda folder: _blast_document(folder, _HITS), None, ["hits.tsv"]),
    "empty": (lambda folder: _blast_document(folder, [("hits.tsv", b"")]), None, ["hits.tsv"]),
    "sanitized": (lambda folder: _blast_document(folder, [("my hits.tsv", _ROWS)]), None, ["my hits.tsv"]),
    "two-of-one-name": (
        lambda folder: _blast_document(folder, [("hits.tsv", _ROWS), ("hits.tsv", _ROWS + _ROWS.splitlines(keepends=True)[-1])]),
        None, ["hits.tsv", "hits.tsv"],
    ),
    "re-encoded": (lambda folder: _reencoded(_blast_document(folder, _HITS), folder), None, ["hits.tsv"]),
    # The Web writes the File's name and names the resource "<id>-<sanitized name>"; Generate sends the File.
    "web-upload": (
        lambda folder: _with_source_name(_older_entry(
            _blast_document(folder, _HITS), "comparison-nucleotide-1", "comparison-nucleotide-1-my_hits.tsv",
        ), "my hits.tsv"),
        {"comparison-nucleotide-1": "my hits.tsv"}, ["my hits.tsv"],
    ),
    "gbdraw-made": (
        lambda folder: _with_source_name(_older_entry(
            _blast_document(folder, _HITS), "comparison-losat-1", "comparison-losat-1.tsv",
        ), None),
        None, [None],
    ),
    # OV-3580: a Web LOSATN Generate stores its result as a table gbdraw made; the page stages it under its resource
    # name, which never names the comparison.
    "web-losat-generate": (
        lambda folder: _with_source_name(_older_entry(
            _blast_document(folder, _HITS), "comparison-resolved-1", "comparison-resolved-1-comparison-resolved-1.tsv",
        ), None),
        {"comparison-resolved-1": "comparison-resolved-1-comparison-resolved-1.tsv"}, [None],
    ),
    # Sessions written before the stored name: the name is read from the resource's name.
    "older-web-upload": (
        lambda folder: _older_entry(_blast_document(folder, _HITS), "comparison-nucleotide-1", "comparison-nucleotide-1-hits.tsv"),
        None, ["hits.tsv"],
    ),
    "older-web-losat": (
        lambda folder: _older_entry(
            _blast_document(folder, _HITS), "comparison-resolved-1", "comparison-resolved-1-comparison-resolved-1.tsv",
        ),
        None, [None],
    ),
    "older-cli-losat": (
        lambda folder: _older_entry(_blast_document(folder, _HITS), "comparison-losat-1", "comparison-losat-1.tsv"),
        None, [None],
    ),
    "older-renamed-losat": (
        lambda folder: _older_entry(
            _blast_document(folder, _HITS), "resource-0001",
            "resource-0001-comparison-resolved-1-comparison-resolved-1.tsv",
        ),
        None, [None],
    ),
    "older-renamed-upload": (
        lambda folder: _older_entry(
            _blast_document(folder, _HITS), "resource-0001", "resource-0001-comparison-nucleotide-1-hits.tsv",
        ),
        None, ["hits.tsv"],
    ),
    "older-gzip-ordinal": (
        lambda folder: _older_entry(
            _blast_document(folder, _HITS), "comparison-losat-1-2", "comparison-losat-1.tsv.2.gz",
            data=gzip.compress(_ROWS, mtime=0),
        ),
        None, [None],
    ),
}


def _with_source_name(document: dict, name: str | None) -> dict:
    for item in document["drawings"][0]["renderRequest"]["comparisons"]:
        if item["kind"] == "nucleotideBlast":
            item["sourceName"] = name
    return document


@pytest.mark.parametrize("case", list(_NAMING_CASES))
def test_a_stored_comparison_is_named_alike_by_the_cli_replay_and_the_web(tmp_path: Path, case: str) -> None:
    """OV-3557 (review H2, M1, M2, L1): one stored owner names a comparison. The Session stores the name of the
    file the user gave (none for a table gbdraw made); the CLI replay and the Web Worker read it through one
    function, and the page's manifest name never overrides it. Older Sessions are read by their resource names."""

    build, manifest, shown = _NAMING_CASES[case]
    document = build(tmp_path / "build")
    payload = document["drawings"][0]["renderRequest"]
    reply = _staged_reply(payload, document["resources"], tmp_path / "web", manifest)
    (summary,) = reply["metadata"]["inputSummaries"]
    web = _comparison_names(summary["lines"])
    cli = _cli_replay_names(document, tmp_path / "cli")
    assert (cli, web) == (shown, shown)


def _record_names(texts: list[str]) -> list[str]:
    """The file name each record line gives its record."""

    return [text.partition(" (")[2].partition("), ")[0] for text in texts if text.startswith("Record ")]


def test_a_web_saved_record_is_named_by_its_file_by_the_cli_replay_and_the_web(tmp_path: Path) -> None:
    """OV-3591: a Web Save names each record resource "<id>-<file>" and stores the file the user chose in
    ``webFiles.resourceOriginalNames``. The CLI replay names each record by that stored name, as its fresh render
    did. The Web half supplies the manifest names itself (the stored names, as the page sends them after a Load), so
    it checks only that the Worker names records by the manifest; the page side, which takes the manifest name
    from the Session resource view, is pinned by tests/web/diagram-resource-staging.test.mjs (OV-3591 case)."""

    request = _blast_request(tmp_path / "build", [("hits.tsv", _ROWS)])
    document = _embedded(build_session_document(request))
    fresh = _record_names(list(input_summary_lines(render_request(request).input_summary)))
    assert fresh == [QUERY.name, SUBJECT.name]
    stored = {}
    for record in document["drawings"][0]["renderRequest"]["records"]:
        resource_id = record["source"]["resourceId"]
        entry = document["resources"][resource_id]
        stored[resource_id] = entry["name"]
        entry["name"] = f"{resource_id}-{entry['name']}"
    document.setdefault("webFiles", {})["resourceOriginalNames"] = stored
    payload = document["drawings"][0]["renderRequest"]
    reply = _staged_reply(payload, document["resources"], tmp_path / "web", stored)
    (summary,) = reply["metadata"]["inputSummaries"]
    from gbdraw.session import materialize_session, render_session

    with materialize_session(load_session_mapping(document), output_directory=tmp_path / "replay") as materialized:
        cli = _record_names(list(input_summary_lines(render_session(materialized).input_summary)))
    assert (cli, _record_names(summary["lines"])) == (fresh, fresh)


_SESSIONS = ROOT / "tests" / "fixtures" / "sessions"
# OV-3660: each fixture and the names Web Load gives its record files (logs of origin/dev 3377d9a7e, through
# config.js importSession and Python's upgrade). The rule that names each: a binding name (Session 41), a CLI file
# binding's slot (Sessions 31-33, CLI 42, and a Circular composite's components of CLI 44), the "<resource ID>-"
# prefix (Session 33), and a stored name (Session 50).
_LOADED_RECORD_NAMES = {
    "single.v41-bindings1.json": ["original.gb"],
    "lambda_basic_linear.v31-schema1.gbdraw-session.json.gz": ["NC_001416.gb"],
    "se08-main-linear-cli.v42.gbdraw-session.json.gz": ["record.gb"],
    "feature-edits-linear-crop-rc.v33.gbdraw-session.json.gz": ["record-0.gb", "record-1.gb"],
    "BGC0000708-BGC0000713.schema-v2.gbdraw-session.json.gz": [
        "BGC0000708.gbk", "BGC0000709.gbk", "BGC0000711.gbk", "BGC0000712.gbk", "BGC0000713.gbk",
    ],
    "HmmtDNA-wrap-5kb.v50-web0145.gbdraw-session.json.gz": ["HmmtDNA.gbk"],
    "composite-circular-three-files.v44-schema8.gbdraw-session.json.gz": [
        "SARS-CoV-1.gbk", "SARS-CoV-2.gbk", "MERS-CoV.gbk",
    ],
}


@pytest.mark.parametrize("fixture", list(_LOADED_RECORD_NAMES))
def test_the_cli_replay_names_an_older_sessions_records_as_web_load_does(tmp_path: Path, fixture: str) -> None:
    """OV-3660: `gbdraw render --session` names each record by the file Web Load shows for it, for a Session that
    stores no name: Python's Session migration names the files once, for both."""

    completed = _gbdraw("render", "--session", str(_SESSIONS / fixture), "-o", "out", "-f", "svg", cwd=tmp_path)
    assert completed.returncode == 0, completed.stderr[-2000:]
    lines = [line.removeprefix("INFO: ") for line in completed.stderr.splitlines() if line.startswith("INFO: Record ")]
    assert _record_names(lines) == _LOADED_RECORD_NAMES[fixture]


_RESOURCE_ID_PREFIX_VECTORS = json.loads(
    (ROOT / "tests" / "fixtures" / "resource-id-prefix-vectors.json").read_text(encoding="utf-8")
)["cases"]


@pytest.mark.parametrize("case", _RESOURCE_ID_PREFIX_VECTORS, ids=lambda case: case["case"])
def test_python_names_a_resource_that_stores_no_name_as_the_shared_vector_says(case: dict[str, str]) -> None:
    """OV-3660: Python's prefix strip gives each case the name the Web's ``withoutResourceIdPrefix`` gives it
    (tests/web/resource-id-prefix.test.mjs reads the same vector; tools/twins.json registers the twin)."""

    from gbdraw.session_migration import _without_resource_id_prefix

    assert _without_resource_id_prefix(case["resourceId"], case["name"]) == case["shown"]


@pytest.mark.parametrize(
    ("binding_name", "stored", "shown"),
    [
        # A binding's name comes first, as in Web Load (``resolveResourceFile``), then the stored name.
        ("a.gb", {"resource-0001": "b.gb"}, "a.gb"),
        (None, {"resource-0001": "b.gb"}, "b.gb"),
        # A schema 1 binding without a name: the file is named as a resource that stores none.
        ("", {}, "x.gb"),
        (None, {}, "x.gb"),
    ],
)
def test_loaded_names_follow_web_loads_order(binding_name: str | None, stored: dict[str, str], shown: str) -> None:
    """OV-3660, review r1 LOW 2: Python names a loaded file in the order Web Load does."""

    from gbdraw.session_migration import loaded_resource_display_names

    binding = {"resourceId": "resource-0001", **({} if binding_name is None else {"name": binding_name})}
    session = {
        "resources": {"resource-0001": {"kind": "web-file", "name": "resource-0001-x.gb"}},
        "webFiles": {"bindings": {"c_gb": binding}, "resourceOriginalNames": stored},
    }
    assert loaded_resource_display_names(session, cli_file_bindings=False) == {"resource-0001": shown}


def test_the_cli_replay_names_two_files_of_one_name_as_its_fresh_render(tmp_path: Path) -> None:
    """OV-3686: a CLI Save renames the second of two files named alike ("record.2.gb") and binds it by the name
    given, which Web Load shows; the replay names it so too, as the render that saved it did."""

    for folder, source in (("a", QUERY), ("b", SUBJECT)):
        (tmp_path / folder).mkdir()
        (tmp_path / folder / "record.gb").write_bytes(source.read_bytes())
    saved = _gbdraw("linear", "--gbk", "a/record.gb", "b/record.gb", "-o", "dup", "-f", "svg", "--save_session",
                    cwd=tmp_path)
    replayed = _gbdraw("render", "--session", "dup.gbdraw", "-o", "replay", "-f", "svg", cwd=tmp_path)
    assert (saved.returncode, replayed.returncode) == (0, 0), saved.stderr[-2000:] + replayed.stderr[-2000:]
    names = [
        _record_names([line.removeprefix("INFO: ") for line in run.stderr.splitlines() if line.startswith("INFO: Record ")])
        for run in (saved, replayed)
    ]
    assert names == [["record.gb", "record.gb"]] * 2


@pytest.mark.parametrize("name", ["hits.tsv", "my hits.tsv"])
def test_a_session_stores_the_name_its_fresh_render_shows(tmp_path: Path, name: str) -> None:
    """The fresh render names a comparison file by its name as given, unsanitized, and the Session stores that name
    on the comparison (review L1)."""

    request = _blast_request(tmp_path, [(name, _ROWS)])
    summary = render_request(request).input_summary
    assert _comparison_names([item.message for item in input_warnings(summary)]) == [name]
    document = _embedded(build_session_document(request))
    (comparison,) = [
        item for item in document["drawings"][0]["renderRequest"]["comparisons"] if item["kind"] == "nucleotideBlast"
    ]
    assert comparison["sourceName"] == name


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


def test_one_file_given_twice_warns_once_in_the_web_reply(tmp_path: Path) -> None:
    """The Web reply carries one warning of one input, as ``gbdraw linear`` logs it."""
    document = _embedded(build_session_document(LinearDiagramRequest(
        records=(RecordInput(source=GenBankInputSource(ZERO_FEATURES)),) * 2,
        output=RenderOutputRequest(output_prefix="twice", formats=("svg",)),
    )))
    reply = render_embedded_canonical_web_request(
        document["drawings"][0]["renderRequest"], resources=document["resources"], workspace=tmp_path / "web",
    )

    assert reply["metadata"]["inputSummaries"][0]["warnings"] == [
        {"code": "no_drawable_features", "message": ZERO_FEATURES_WARNING.removeprefix("WARNING: ")}
    ]


def test_records_without_a_file_each_warn_in_the_python_api(caplog: pytest.LogCaptureFixture) -> None:
    """A record held in memory has no file to name it: each one drawn is its own input and warns."""
    from gbdraw import read_genbank

    records = read_genbank(str(ZERO_FEATURES))
    with caplog.at_level("WARNING", logger="gbdraw.api.input_summary"):
        gbdraw.draw_linear(records + records)

    warning = ZERO_FEATURES_WARNING.removeprefix("WARNING: ").replace("(zero_features_record.gb)", "(in-memory record)")
    assert [record.getMessage() for record in caplog.records if record.levelname == "WARNING"] == [warning, warning]


def test_two_comparison_tables_of_one_name_and_record_pair_each_warn() -> None:
    """Two empty tables of one base name (two folders) for one record pair are two inputs."""
    empty = ComparisonTableTally("hits.tsv", 0, 1, ComparisonRowTally(rows_read=0, rows_kept=0, thresholds=()))

    warnings = input_warnings(InputSummary(records=(), comparisons=(empty, empty)))

    assert [item.message for item in warnings] == [
        "Comparison table hits.tsv (record #1 and record #2) has no rows, so no match is drawn."
    ] * 2


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


def test_a_protein_comparison_of_listed_record_pairs_keeps_its_protein_wording(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path,
) -> None:
    """LOSATP Pairwise with listed pairs: the pair comparisons are protein comparisons, one with no match."""

    from dataclasses import replace

    import gbdraw.analysis.protein_colinearity as protein_colinearity_module
    from tests.test_render_stage_protein_analysis import _fake_losatp, _request

    def search(query_fasta: str, subject_fasta: str, **kwargs) -> object:
        # Record #3 (rec2) shares no protein with record #2: that pair finds nothing.
        if "rec2_" in subject_fasta:
            kwargs["raw_output_callback"]("")
            kwargs["runtime_callback"]({"kind": "losat", "version": "test", "program": "blastp"})
            return protein_colinearity_module.parse_losatp_outfmt6("")
        return _fake_losatp(query_fasta, subject_fasta, **kwargs)

    monkeypatch.setattr(protein_colinearity_module, "run_losatp_blastp", search)
    request, _identities = _request("pairwise", tmp_path)
    options = replace(
        request.options,
        losat_search=replace(request.options.losat_search, pairs=((0, 1), (1, 2))),
    )
    request = replace(request, options=options, output=RenderOutputRequest(
        output_prefix="pairs", formats=("svg",), output_directory=tmp_path / "out",
    ))
    summary = render_request(request).input_summary

    lines = [line for line in input_summary_lines(summary) if "omparison of" in line]
    assert "Protein comparison of record #1 rec0 and record #2 rec1: drew 3 protein matches." in lines, lines
    assert "Protein comparison of record #2 rec1 and record #3 rec2: drew 0 protein matches." in lines, lines
    warnings = [item.message for item in input_warnings(summary)]
    assert len(warnings) == 1 and warnings[0].endswith("has no protein match to draw."), warnings


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


def _genbank_text(name: str, topology: str, locations: Sequence[tuple[str, str]]) -> str:
    """A 100 bp GenBank record whose CDS features have ``locations``, each with a ``locus_tag`` or a ``product``."""

    features = "".join(
        f"     CDS             {location}\n                     /{qualifier}\n" for location, qualifier in locations
    )
    sequence = "acgtgcatgc" * 10
    origin = "".join(
        f"{start + 1:>9} " + " ".join(sequence[at:at + 10] for at in range(start, min(start + 60, 100), 10)) + "\n"
        for start in range(0, 100, 60)
    )
    return (
        f"LOCUS       {name:<16}         100 bp    DNA     {topology} UNK 01-JAN-1980\n"
        f"DEFINITION  .\nACCESSION   {name}\nVERSION     {name}\nKEYWORDS    .\nSOURCE      .\n  ORGANISM  .\n"
        f"            .\nFEATURES             Location/Qualifiers\n     source          1..100\n{features}"
        f"ORIGIN\n{origin}//\n"
    )


# OV-3678: a 100 bp circular record with a CDS that fits and one that starts after the record's end.
OUTSIDE_RECORD = (("21..90", 'product="p"'), ("121..240", 'product="q"'))
OUTSIDE_RECORD_WARNING = (
    "WARNING: Record tiny (tiny.gb): left out 1 feature that does not fit in the record's 100 bp: CDS 121..240 (q)."
)


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_cli_leaves_out_a_feature_outside_its_record_and_names_it(tmp_path: Path, mode: str) -> None:
    """A feature that does not fit in its record is left out with a warning that names it (OV-3678): Circular
    stopped with "A circle arc spans at most 360 degrees", Linear drew it past the end of the record."""
    (tmp_path / "tiny.gb").write_text(_genbank_text("tiny", "circular", OUTSIDE_RECORD), encoding="utf-8")

    run = _gbdraw(mode, "--gbk", "tiny.gb", "-o", "out", "-f", "svg", cwd=tmp_path)

    assert run.returncode == 0, run.stderr
    lines = _lines(run)
    assert OUTSIDE_RECORD_WARNING in lines
    assert "INFO: Record tiny (tiny.gb), record 1 of 1, 100 bp: drew 1 feature (1 CDS); left out 1 source (type not drawn)." in lines


def test_web_reply_names_a_feature_outside_its_record(tmp_path: Path) -> None:
    from gbdraw.api import CircularDiagramRequest

    path = tmp_path / "tiny.gb"
    path.write_text(_genbank_text("tiny", "circular", OUTSIDE_RECORD), encoding="utf-8")
    document = _embedded(build_session_document(CircularDiagramRequest(
        records=(RecordInput(source=GenBankInputSource(path)),),
        output=RenderOutputRequest(output_prefix="tiny", formats=("svg",)),
    )))
    reply = render_embedded_canonical_web_request(
        document["drawings"][0]["renderRequest"], resources=document["resources"], workspace=tmp_path / "web",
    )

    assert reply["metadata"]["inputSummaries"][0]["warnings"] == [
        {"code": "feature_outside_record", "message": OUTSIDE_RECORD_WARNING.removeprefix("WARNING: ")}
    ]


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_python_api_names_a_feature_outside_its_record(
    tmp_path: Path, mode: str, caplog: pytest.LogCaptureFixture
) -> None:
    from gbdraw import read_genbank

    path = tmp_path / "tiny.gb"
    path.write_text(_genbank_text("tiny", "circular", OUTSIDE_RECORD), encoding="utf-8")
    records = read_genbank(str(path))
    with caplog.at_level("WARNING", logger="gbdraw.api.input_summary"):
        diagram = (gbdraw.draw_circular if mode == "circular" else gbdraw.draw_linear)(records)

    assert [record.getMessage() for record in caplog.records if record.levelname == "WARNING"] == [
        OUTSIDE_RECORD_WARNING.removeprefix("WARNING: ").replace("(tiny.gb)", "(in-memory record)")
    ]
    assert diagram is not None


def test_the_warning_names_each_feature_by_its_type_name_and_location(tmp_path: Path) -> None:
    """Up to three features are named, by locus_tag (else gene, label, product), the rest counted."""
    path = tmp_path / "many.gb"
    path.write_text(_genbank_text("many", "linear", (
        ("1..10", 'locus_tag="in1"'),
        ("complement(101..110)", 'locus_tag="out1"'),
        ("join(1..10,150..160)", 'gene="out2"'),
        ("1..150", 'label="out3"'),
        ("200..210", 'note="no name"'),
    )), encoding="utf-8")

    run = _gbdraw("linear", "--gbk", "many.gb", "-o", "out", "-f", "svg", cwd=tmp_path)

    assert run.returncode == 0, run.stderr
    assert (
        "WARNING: Record many (many.gb): left out 4 features that do not fit in the record's 100 bp: "
        "CDS complement(101..110) (out1), CDS join(1..10,150..160) (out2), CDS 1..150 (out3), and 1 more."
    ) in _lines(run)


@pytest.mark.parametrize("mode", ["circular", "linear"])
def test_features_that_fit_their_record_draw_as_before_without_a_warning(tmp_path: Path, mode: str) -> None:
    """Joins over the origin, partial ends, and a span over a circular record's origin written past its end (the
    GFF3 way) all fit: each is drawn and nothing warns."""
    path = tmp_path / "fit.gb"
    path.write_text(_genbank_text("fit", "circular", (
        ("join(91..100,1..10)", 'locus_tag="join"'),
        ("<1..30", 'locus_tag="partial5"'),
        ("70..>100", 'locus_tag="partial3"'),
        ("81..110", 'locus_tag="past_end"'),
        ("1..100", 'locus_tag="whole"'),
    )), encoding="utf-8")

    run = _gbdraw(mode, "--gbk", "fit.gb", "-o", "out", "-f", "svg", cwd=tmp_path)

    assert run.returncode == 0, run.stderr
    lines = _lines(run)
    assert not [line for line in lines if line.startswith("WARNING")]
    assert "INFO: Record fit (fit.gb), record 1 of 1, 100 bp: drew 5 features (5 CDS); left out 1 source (type not drawn)." in lines


def test_features_after_one_left_out_keep_their_source_identity(tmp_path: Path) -> None:
    """Leaving out a feature that does not fit moves no other feature's source index or hash: edits and colors that
    name a later feature still find it."""
    from gbdraw.api.request_render import resolve_request_records
    from gbdraw.core.record_metadata import _feature_source_index_map
    from gbdraw.features.ids import compute_source_feature_hash

    path = tmp_path / "tiny.gb"
    path.write_text(_genbank_text("tiny", "circular", (("121..240", 'locus_tag="out"'), ("21..90", 'locus_tag="in"'))),
                    encoding="utf-8")
    (source,) = read_flat_file_records(str(path))
    resolved = resolve_request_records(LinearDiagramRequest(records=(RecordInput(GenBankInputSource(path)),)))

    (record,) = resolved.records
    indexes = _feature_source_index_map(record.features)
    assert [(feature.qualifiers.get("locus_tag"), indexes[id(feature)]) for feature in record.features] == [
        (None, 0), (["in"], 2),
    ]
    assert compute_source_feature_hash(record.features[1], record) == compute_source_feature_hash(
        source.features[2], source
    )
    assert [item.location for item in resolved.provenance[0].features_outside_record] == ["121..240"]
