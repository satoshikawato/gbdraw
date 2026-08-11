"""Capture whole-record nucleotide-comparison GUI How-to journeys."""

from __future__ import annotations

import csv
import hashlib
import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from decimal import Decimal
from pathlib import Path
from typing import Any, Mapping

from Bio import SeqIO
from playwright.sync_api import BrowserType, Page, expect

from assertions.svg_semantics import (
    assert_plain_gui_losatn_svg,
    assert_static_svg_safety,
    inspect_first_circular_svg,
    inspect_gui_losatn_svg,
)
from config import (
    ACTION_TIMEOUT_MS,
    FIRST_CIRCULAR_FIXTURE_PATH,
    FIRST_CIRCULAR_FIXTURE_SHA256,
    FIRST_CIRCULAR_FIXTURE_SIZE,
    FIRST_LINEAR_FIXTURE_PATH,
    FIRST_LINEAR_FIXTURE_SHA256,
    FIRST_LINEAR_FIXTURE_SIZE,
    FIRST_LINEAR_LABEL_RULE_PATH,
    FIRST_LINEAR_LABEL_RULE_SHA256,
    FIRST_LINEAR_LABEL_RULE_SIZE,
    GUI_CIRCULAR_RINGS_SCREENSHOT_NAMES,
    GUI_CIRCULAR_LAYOUT_FIXTURES,
    GUI_TLOSATX_SCREENSHOT_NAMES,
    GUI_LOSATN_DE3_FIXTURE_PATH,
    GUI_LOSATN_DE3_FIXTURE_SHA256,
    GUI_LOSATN_DE3_FIXTURE_SIZE,
    GUI_UPLOADED_COMPARISON_SCREENSHOT_NAMES,
    WEB_ROOT,
    WORKER_READY_TIMEOUT_MS,
)
from flows.web_capture import (
    assert_fixture_identity,
    assert_output_paths,
    capture_screenshot,
    fit_complete_linear_preview,
    generate_and_inspect,
    linear_pair,
    open_browser_capture,
    open_linear_comparison_disclosure,
    select_linear_losat_mode,
    wait_for_worker,
)


EXPECTED_UPLOADED_COMPARISON_SVG = "uploaded_comparison.svg"
EXPECTED_TLOSATX_SVG = "lambda-de3-tlosatx.svg"
EXPECTED_TLOSATX_TSV = "lambda-de3.tlosatx.tsv"
EXPECTED_CIRCULAR_RING_SVG = "circular_similarity_rings.svg"
EXPECTED_CIRCULAR_SPAN_FASTA = "comparison_spans.fasta"
CIRCULAR_RING_TITLE = "Homo sapiens mtDNA (NC_012920.1) similarity rings"
CIRCULAR_RING_LABELS = (
    "Danio rerio (NC_002333.2)",
    "Drosophila melanogaster (NC_024511.2)",
    "Caenorhabditis elegans (NC_001328.1)",
)

TLOSATX_REFERENCE_PATH = (
    WEB_ROOT
    / "tutorial-data"
    / "lambda-de3-comparison"
    / EXPECTED_TLOSATX_TSV
)
TLOSATX_REFERENCE_SIZE = 28_303
TLOSATX_REFERENCE_SHA256 = (
    "483e98f8b3dce172523cf00c82eb9d47e4faee13e4781b06f3d58ea4fb63532d"
)


@dataclass(frozen=True)
class ComparisonCaptureResult:
    screenshot_bytes: dict[str, int]
    final_svg_semantics: dict[str, Any]
    download: dict[str, Any]
    tsv_download: dict[str, Any] | None = None
    fasta_download: dict[str, Any] | None = None
    popup: dict[str, Any] | None = None
    source_records: tuple[dict[str, Any], ...] = ()


def _read_blast_rows(path: Path) -> tuple[tuple[str, ...], ...]:
    rows = []
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if not row or row[0].startswith("#"):
                continue
            if len(row) < 12:
                raise AssertionError(f"Malformed BLAST row in {path.name}: {row!r}")
            rows.append(tuple(row[:12]))
    return tuple(rows)


def _qualified_tlosatx_rows() -> tuple[tuple[str, ...], ...]:
    assert_fixture_identity(
        TLOSATX_REFERENCE_PATH,
        expected_size=TLOSATX_REFERENCE_SIZE,
        expected_sha256=TLOSATX_REFERENCE_SHA256,
    )
    rows = _read_blast_rows(TLOSATX_REFERENCE_PATH)
    if len(rows) != 397:
        raise AssertionError(f"Expected 397 frozen TLOSATX rows, found {len(rows)}")
    endpoints = {(row[0], row[1]) for row in rows}
    if endpoints != {("NC_001416.1", "NC_042057.1")}:
        raise AssertionError(f"Unexpected frozen TLOSATX endpoints: {endpoints!r}")
    retained = tuple(row for row in rows if int(row[3]) >= 1_000)
    if len(retained) != 7:
        raise AssertionError(
            f"Expected seven TLOSATX rows at 1,000 aa, found {len(retained)}"
        )
    return retained


def _numeric_blast_signature(row: tuple[str, ...]) -> tuple[Any, ...]:
    return (*row[:2], *(Decimal(value) for value in row[2:]))


def _rendered_match_signature(match: dict[str, Any]) -> tuple[Any, ...]:
    row = (
        str(match.get("query", "")),
        str(match.get("subject", "")),
        str(match.get("identity", "")),
        str(match.get("alignmentLength", "")),
        str(match.get("mismatches", "")),
        str(match.get("gapOpens", "")),
        str(match.get("qstart", "")),
        str(match.get("qend", "")),
        str(match.get("sstart", "")),
        str(match.get("send", "")),
        str(match.get("evalue", "")),
        str(match.get("bitscore", "")),
    )
    return _numeric_blast_signature(row)


def _assert_tlosatx_linear_svg(report: dict[str, Any]) -> None:
    assert_plain_gui_losatn_svg({**report, "matches": []})
    actual = sorted(
        _rendered_match_signature(match) for match in report.get("matches", [])
    )
    expected = sorted(
        _numeric_blast_signature(row) for row in _qualified_tlosatx_rows()
    )
    if actual != expected:
        raise AssertionError(
            "Rendered TLOSATX links differ from the seven qualified frozen rows"
        )
    for match in report.get("matches", []):
        qstart = int(match["qstart"])
        qend = int(match["qend"])
        sstart = int(match["sstart"])
        send = int(match["send"])
        if not (1 <= min(qstart, qend) <= max(qstart, qend) <= 48_502):
            raise AssertionError(f"TLOSATX query span is outside complete Lambda: {match}")
        if not (1 <= min(sstart, send) <= max(sstart, send) <= 42_925):
            raise AssertionError(f"TLOSATX subject span is outside complete DE3: {match}")
    assert_static_svg_safety(report)


def _assert_comparison_svg_download(
    path: Path,
    *,
    expected_matches: int,
) -> dict[str, Any]:
    contents = path.read_bytes()
    if len(contents) < 20_000:
        raise AssertionError(f"Downloaded comparison SVG is too small: {len(contents)}")
    root = ET.fromstring(contents)
    if root.tag.rsplit("}", 1)[-1] != "svg":
        raise AssertionError("Downloaded comparison output is not an SVG")
    elements = list(root.iter())
    if any(element.tag.rsplit("}", 1)[-1] == "script" for element in elements):
        raise AssertionError("Static comparison SVG contains a script")
    record_ids = {
        element.attrib.get("data-gbdraw-record-id", "") for element in elements
    }
    if not {"NC_001416.1", "NC_042057.1"} <= record_ids:
        raise AssertionError("Downloaded comparison SVG lost a whole input record")
    matches = [
        element
        for element in elements
        if element.attrib.get("data-gbdraw-pairwise-match-id")
        and element.attrib.get("data-match-kind") == "pairwise"
    ]
    if len(matches) != expected_matches:
        raise AssertionError(
            f"Expected {expected_matches} downloaded links, found {len(matches)}"
        )
    for match in matches:
        if (
            match.attrib.get("data-query-record-id") != "NC_001416.1"
            or match.attrib.get("data-subject-record-id") != "NC_042057.1"
        ):
            raise AssertionError("Downloaded comparison endpoints changed")
    return {
        "filename": path.name,
        "bytes": len(contents),
        "record_ids": sorted(record_ids - {""}),
        "matches": len(matches),
    }


def _assert_tlosatx_download(path: Path) -> dict[str, Any]:
    contents = path.read_bytes()
    expected = TLOSATX_REFERENCE_PATH.read_bytes()
    if contents != expected:
        raise AssertionError("Browser TLOSATX TSV differs from the frozen qualified evidence")
    rows = _read_blast_rows(path)
    if len(rows) != 397:
        raise AssertionError(f"Expected 397 downloaded TLOSATX rows, found {len(rows)}")
    return {
        "filename": path.name,
        "bytes": len(contents),
        "rows": len(rows),
        "sha256": hashlib.sha256(contents).hexdigest(),
        "endpoints": sorted({(row[0], row[1]) for row in rows}),
    }


def _load_complete_linear_inputs(page: Page) -> None:
    linear = page.get_by_role("button", name="Linear", exact=True)
    linear.click()
    expect(linear).to_have_attribute("aria-pressed", "true")
    page.get_by_role("radio", name="GenBank", exact=True).check()
    expect(page.get_by_role("status").filter(has_text="Current:")).to_contain_text(
        "Current: No comparison"
    )
    add_sequence = page.get_by_role(
        "button", name="Add sequence", exact=True
    )
    expect(add_sequence).to_have_count(2)
    add_sequence.first.click()
    page.get_by_test_id("linear-genbank-1").set_input_files(FIRST_LINEAR_FIXTURE_PATH)
    page.get_by_test_id("linear-genbank-2").set_input_files(
        GUI_LOSATN_DE3_FIXTURE_PATH
    )
    selections = page.get_by_role("group", name="GenBank File selection", exact=True)
    expect(selections).to_have_count(2)
    expect(selections.nth(0)).to_contain_text("NC_001416.gb")
    expect(selections.nth(1)).to_contain_text("NC_042057.1.gb")
    for index, record_id in ((1, "NC_001416.1"), (2, "NC_042057.1")):
        selector = page.get_by_label(
            f"Record selector for sequence {index}", exact=True
        )
        expect(selector).to_be_enabled(timeout=WORKER_READY_TIMEOUT_MS)
        expect(selector).to_contain_text(
            record_id, timeout=WORKER_READY_TIMEOUT_MS
        )


def _set_linear_comparison_filter(page: Page, minimum_length: int) -> None:
    settings = open_linear_comparison_disclosure(
        page,
        "settings",
        "Comparison Settings",
    )
    minimum = settings.get_by_label(
        "Linear comparison minimum alignment length", exact=True
    )
    minimum.fill(str(minimum_length))
    minimum.press("Tab")
    expect(minimum).to_have_value(str(minimum_length))


def _download_svg(
    page: Page,
    download_dir: Path,
    expected_name: str,
    *,
    expected_matches: int,
) -> dict[str, Any]:
    with page.expect_download(timeout=ACTION_TIMEOUT_MS) as download_info:
        page.get_by_role("button", name="SVG", exact=True).click()
    download = download_info.value
    if download.failure() is not None:
        raise AssertionError(f"SVG download failed: {download.failure()}")
    if download.suggested_filename != expected_name:
        raise AssertionError(
            f"Expected SVG filename {expected_name}, found {download.suggested_filename}"
        )
    path = download_dir / download.suggested_filename
    download.save_as(path)
    return _assert_comparison_svg_download(path, expected_matches=expected_matches)


def capture_gui_uploaded_comparison(
    browser_type: BrowserType,
    base_url: str,
    output_paths: Mapping[str, Path],
    download_dir: Path,
) -> ComparisonCaptureResult:
    """Run H-GUI-04 with the frozen whole-Lambda versus whole-DE3 table."""

    for path, size, digest in (
        (
            FIRST_LINEAR_FIXTURE_PATH,
            FIRST_LINEAR_FIXTURE_SIZE,
            FIRST_LINEAR_FIXTURE_SHA256,
        ),
        (
            GUI_LOSATN_DE3_FIXTURE_PATH,
            GUI_LOSATN_DE3_FIXTURE_SIZE,
            GUI_LOSATN_DE3_FIXTURE_SHA256,
        ),
    ):
        assert_fixture_identity(path, expected_size=size, expected_sha256=digest)
    _qualified_tlosatx_rows()
    assert_output_paths(
        output_paths,
        GUI_UPLOADED_COMPARISON_SCREENSHOT_NAMES,
        "H-GUI-04",
    )
    download_dir.mkdir(parents=True, exist_ok=True)

    capture = open_browser_capture(browser_type, base_url)
    page = capture.page
    screenshot_bytes: dict[str, int] = {}
    try:
        page.goto(base_url, wait_until="domcontentloaded")
        wait_for_worker(page)
        _load_complete_linear_inputs(page)

        commands = page.get_by_role(
            "group", name="Set all adjacent comparisons", exact=True
        )
        use_upload = commands.get_by_role(
            "button",
            name="Use uploaded BLAST TSV for all adjacent pairs",
            exact=True,
        )
        use_upload.click()
        expect(page.get_by_role("status").filter(has_text="Current:")).to_contain_text(
            "Current: Upload BLAST TSV for all adjacent pairs"
        )
        open_linear_comparison_disclosure(
            page,
            "settings",
            "Comparison Settings",
        )
        selected_pairs = open_linear_comparison_disclosure(
            page,
            "selected-pairs",
            "Selected pairs",
        )
        pair = linear_pair(page, 1, 2)
        edge_upload = pair.get_by_role(
            "radio", name="Upload BLAST TSV", exact=True
        )
        expect(edge_upload).to_be_checked()
        blast_input = pair.get_by_label(
            "BLAST TSV for #1 to #2", exact=True
        )
        blast_input.set_input_files(TLOSATX_REFERENCE_PATH)
        blast_selection = pair.get_by_role(
            "group", name="BLAST TSV for #1 to #2 selection", exact=True
        )
        expect(blast_selection).to_contain_text(EXPECTED_TLOSATX_TSV)
        _set_linear_comparison_filter(page, 1_000)

        prefix = page.get_by_label("Output Prefix", exact=True)
        prefix.fill("uploaded_comparison")
        expect(prefix).to_have_value("uploaded_comparison")
        expect(edge_upload).to_be_checked()
        selected_pairs.scroll_into_view_if_needed()
        blast_selection.scroll_into_view_if_needed()
        blast_selection.hover()
        page.mouse.wheel(0, 240)
        screenshot_bytes["comparison-plan.png"] = capture_screenshot(
            page, output_paths["comparison-plan.png"], "Linear"
        )

        final_report = generate_and_inspect(
            page,
            inspect_gui_losatn_svg,
            _assert_tlosatx_linear_svg,
        )
        fit_complete_linear_preview(page)
        screenshot_bytes["comparison-result.png"] = capture_screenshot(
            page, output_paths["comparison-result.png"], "Linear"
        )
        download_report = _download_svg(
            page,
            download_dir,
            EXPECTED_UPLOADED_COMPARISON_SVG,
            expected_matches=7,
        )
        capture.assert_clean()
    finally:
        capture.close()

    return ComparisonCaptureResult(
        screenshot_bytes=screenshot_bytes,
        final_svg_semantics=final_report,
        download=download_report,
    )


def capture_gui_tlosatx(
    browser_type: BrowserType,
    base_url: str,
    output_paths: Mapping[str, Path],
    download_dir: Path,
) -> ComparisonCaptureResult:
    """Run H-GUI-05 as a serial one-thread browser TLOSATX search."""

    for path, size, digest in (
        (
            FIRST_LINEAR_FIXTURE_PATH,
            FIRST_LINEAR_FIXTURE_SIZE,
            FIRST_LINEAR_FIXTURE_SHA256,
        ),
        (
            GUI_LOSATN_DE3_FIXTURE_PATH,
            GUI_LOSATN_DE3_FIXTURE_SIZE,
            GUI_LOSATN_DE3_FIXTURE_SHA256,
        ),
    ):
        assert_fixture_identity(path, expected_size=size, expected_sha256=digest)
    _qualified_tlosatx_rows()
    assert_output_paths(output_paths, GUI_TLOSATX_SCREENSHOT_NAMES, "H-GUI-05")
    download_dir.mkdir(parents=True, exist_ok=True)

    capture = open_browser_capture(browser_type, base_url)
    page = capture.page
    screenshot_bytes: dict[str, int] = {}
    try:
        page.goto(base_url, wait_until="domcontentloaded")
        wait_for_worker(page)
        _load_complete_linear_inputs(page)

        commands = page.get_by_role(
            "group", name="Set all adjacent comparisons", exact=True
        )
        run_losat = commands.get_by_role(
            "button", name="Run LOSAT for all adjacent pairs", exact=True
        )
        run_losat.click()
        expect(page.get_by_role("status").filter(has_text="Current:")).to_contain_text(
            "Current: Run LOSAT for all adjacent pairs"
        )
        settings = open_linear_comparison_disclosure(
            page,
            "settings",
            "Comparison Settings",
        )
        select_linear_losat_mode(
            settings,
            label="TLOSATX",
            mode_key="tblastx",
        )
        for index in (1, 2):
            record_options = page.get_by_role(
                "button",
                name=f"Record options for sequence {index}",
                exact=True,
            )
            record_options.click()
            gencode = page.get_by_label(
                f"TLOSATX gencode for sequence {index}", exact=True
            )
            gencode.fill("1")
            expect(gencode).to_have_value("1")
            record_options.click()

        advanced = open_linear_comparison_disclosure(
            page,
            "advanced",
            "Advanced comparison and layout",
        )
        execution = advanced.get_by_role(
            "combobox", name="LOSAT execution", exact=True
        )
        execution.select_option("serial")
        expect(execution).to_have_value("serial")
        total_threads = advanced.get_by_role(
            "combobox", name="LOSAT total threads", exact=True
        )
        total_threads.select_option("1")
        expect(total_threads).to_have_value("1")
        parallel_runs = advanced.get_by_role(
            "combobox", name="LOSAT parallel runs", exact=True
        )
        parallel_runs.select_option("1")
        expect(parallel_runs).to_have_value("1")
        threads_per_run = advanced.get_by_role(
            "combobox", name="LOSAT threads per run", exact=True
        )
        expect(threads_per_run).to_be_disabled()

        prefix = page.get_by_label("Output Prefix", exact=True)
        prefix.fill("lambda-de3-tlosatx")
        expect(prefix).to_have_value("lambda-de3-tlosatx")
        raw_filename = advanced.get_by_role(
            "textbox", name="Raw LOSAT filename for #1 to #2", exact=True
        )
        raw_filename.fill(EXPECTED_TLOSATX_TSV)
        raw_filename.press("Tab")
        expect(raw_filename).to_have_value(EXPECTED_TLOSATX_TSV)
        _set_linear_comparison_filter(page, 1_000)

        settings.get_by_role(
            "group", name="LOSAT Mode", exact=True
        ).scroll_into_view_if_needed()
        screenshot_bytes["tlosatx-settings.png"] = capture_screenshot(
            page, output_paths["tlosatx-settings.png"], "Linear"
        )

        final_report = generate_and_inspect(
            page,
            inspect_gui_losatn_svg,
            _assert_tlosatx_linear_svg,
        )
        fit_complete_linear_preview(page)
        screenshot_bytes["tlosatx-result.png"] = capture_screenshot(
            page, output_paths["tlosatx-result.png"], "Linear"
        )

        tsv_button = advanced.get_by_role(
            "button", name="Save Raw LOSAT TSV for #1 to #2", exact=True
        )
        expect(tsv_button).to_be_enabled()
        with page.expect_download(timeout=ACTION_TIMEOUT_MS) as tsv_info:
            tsv_button.click()
        tsv_download = tsv_info.value
        if tsv_download.failure() is not None:
            raise AssertionError(f"TLOSATX TSV download failed: {tsv_download.failure()}")
        if tsv_download.suggested_filename != EXPECTED_TLOSATX_TSV:
            raise AssertionError(
                f"Unexpected TLOSATX TSV filename: {tsv_download.suggested_filename}"
            )
        tsv_path = download_dir / tsv_download.suggested_filename
        tsv_download.save_as(tsv_path)
        tsv_report = _assert_tlosatx_download(tsv_path)

        download_report = _download_svg(
            page,
            download_dir,
            EXPECTED_TLOSATX_SVG,
            expected_matches=7,
        )
        capture.assert_clean()
    finally:
        capture.close()

    return ComparisonCaptureResult(
        screenshot_bytes=screenshot_bytes,
        final_svg_semantics=final_report,
        download=download_report,
        tsv_download=tsv_report,
    )


def _write_complete_fasta(source_path: Path, output_path: Path) -> dict[str, Any]:
    records = list(SeqIO.parse(source_path, "genbank"))
    if len(records) != 1:
        raise AssertionError(f"Expected one complete record in {source_path.name}")
    record = records[0]
    if str(record.annotations.get("topology", "")).lower() != "circular":
        raise AssertionError(f"{record.id} is not naturally circular")
    if "complete" not in record.description.lower():
        raise AssertionError(f"{record.id} is not described as complete")
    sequence = str(record.seq).upper()
    wrapped = "\n".join(
        sequence[index : index + 60] for index in range(0, len(sequence), 60)
    )
    output_path.write_text(
        f">{record.id} {record.description}\n{wrapped}\n",
        encoding="ascii",
    )
    reparsed = list(SeqIO.parse(output_path, "fasta"))
    if (
        len(reparsed) != 1
        or reparsed[0].id != record.id
        or str(reparsed[0].seq).upper() != sequence
    ):
        raise AssertionError(f"Whole-record FASTA conversion changed {record.id}")
    return {
        "record_id": record.id,
        "length": len(record),
        "topology": record.annotations.get("topology"),
        "organism": record.annotations.get("organism"),
        "fasta": output_path,
    }


def _prepare_complete_mitochondrial_inputs(
    download_dir: Path,
) -> tuple[dict[str, Any], ...]:
    reports = []
    for path, size, digest, record_id, length, organism in (
        GUI_CIRCULAR_LAYOUT_FIXTURES
    ):
        assert_fixture_identity(path, expected_size=size, expected_sha256=digest)
        if record_id == "NC_012920.1":
            records = list(SeqIO.parse(path, "genbank"))
            record = records[0]
            reports.append(
                {
                    "record_id": record.id,
                    "length": len(record),
                    "topology": record.annotations.get("topology"),
                    "organism": record.annotations.get("organism"),
                    "genbank": path,
                }
            )
            continue
        fasta_path = download_dir / f"{record_id}.fna"
        report = _write_complete_fasta(path, fasta_path)
        if report["record_id"] != record_id or report["length"] != length:
            raise AssertionError(f"Unexpected mitochondrial source record: {report!r}")
        if report["organism"] != organism:
            raise AssertionError(f"Unexpected mitochondrial organism: {report!r}")
        reports.append(report)
    return tuple(reports)


def _inspect_circular_rings(result_region: Any) -> dict[str, Any]:
    report = inspect_first_circular_svg(result_region)
    ring_report = result_region.evaluate(
        r"""
        root => {
          const svg = root.querySelector('svg');
          const rings = Array.from(svg?.querySelectorAll(
            'g[data-gbdraw-slot-renderer="sequence_conservation"]'
          ) || []);
          const hits = Array.from(svg?.querySelectorAll(
            '[data-gbdraw-match-id][data-match-kind="homology"]'
          ) || []);
          return {
            rings: rings.map((ring) => ({
              label: String(ring.getAttribute('data-track-label') || ''),
              sourceIndex: Number(ring.getAttribute('data-source-index')),
              trackIndex: Number(ring.getAttribute('data-track-index')),
              color: String(ring.getAttribute('data-track-color') || ''),
              slotId: String(ring.getAttribute('data-gbdraw-slot-id') || '')
            })),
            hits: hits.map((hit) => ({
              id: String(hit.getAttribute('data-gbdraw-match-id') || ''),
              label: String(hit.getAttribute('data-track-label') || ''),
              referenceSide: String(hit.getAttribute('data-reference-side') || ''),
              referenceRecordId: String(hit.getAttribute('data-reference-record-id') || ''),
              query: String(hit.getAttribute('data-query-record-id') || ''),
              subject: String(hit.getAttribute('data-subject-record-id') || ''),
              qstart: Number(hit.getAttribute('data-qstart')),
              qend: Number(hit.getAttribute('data-qend')),
              sstart: Number(hit.getAttribute('data-sstart')),
              send: Number(hit.getAttribute('data-send')),
              identity: Number(hit.getAttribute('data-identity')),
              alignmentLength: Number(hit.getAttribute('data-alignment-length'))
            }))
          };
        }
        """
    )
    report.update(ring_report)
    return report


def _interval_union_length(intervals: list[tuple[int, int]]) -> int:
    if not intervals:
        return 0
    merged = []
    for start, end in sorted((min(a, b), max(a, b)) for a, b in intervals):
        if not merged or start > merged[-1][1] + 1:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return sum(end - start + 1 for start, end in merged)


def _assert_circular_ring_svg(report: dict[str, Any]) -> None:
    if report.get("svgCount") != 1:
        raise AssertionError(f"Expected one Circular SVG, found {report.get('svgCount')}")
    if "NC_012920.1" not in set(report.get("recordIds", [])):
        raise AssertionError("The Circular reference is not complete human mtDNA")
    if report.get("featureElementCount") != 37:
        raise AssertionError(
            f"Expected 37 human mitochondrial feature elements, found {report.get('featureElementCount')}"
        )
    texts = set(report.get("texts", []))
    for expected in (
        "NC_012920.1",
        "16,569 bp",
        *CIRCULAR_RING_LABELS,
    ):
        if expected not in texts:
            raise AssertionError(f"Circular comparison SVG is missing {expected!r}")
    if not any("Homo sapiens" in text for text in texts):
        raise AssertionError("Circular comparison SVG is missing the human label")

    rings = report.get("rings", [])
    if [ring.get("label") for ring in rings] != list(CIRCULAR_RING_LABELS):
        raise AssertionError(f"Unexpected comparison-ring order: {rings!r}")
    hits = report.get("hits", [])
    if len(hits) < 3:
        raise AssertionError(f"Expected retained HSPs in three rings, found {len(hits)}")
    hit_labels = {hit.get("label") for hit in hits}
    if hit_labels != set(CIRCULAR_RING_LABELS):
        raise AssertionError(f"One or more complete mtDNA comparisons are empty: {hit_labels}")
    for hit in hits:
        if hit.get("referenceSide") != "subject":
            raise AssertionError(f"The displayed reference side changed: {hit!r}")
        if hit.get("referenceRecordId") != "NC_012920.1":
            raise AssertionError(f"A ring is not mapped to whole human mtDNA: {hit!r}")
        if hit.get("subject") != "NC_012920.1":
            raise AssertionError(f"TLOSATX subject is not human mtDNA: {hit!r}")
        if not (1 <= min(hit["sstart"], hit["send"]) <= 16_569):
            raise AssertionError(f"Human mtDNA HSP starts outside the complete record: {hit!r}")
        if max(hit["sstart"], hit["send"]) > 16_569:
            raise AssertionError(f"Human mtDNA HSP ends outside the complete record: {hit!r}")
    union_coverage = _interval_union_length(
        [(hit["sstart"], hit["send"]) for hit in hits]
    )
    if union_coverage <= 0 or union_coverage > 16_569:
        raise AssertionError(f"Invalid human mtDNA HSP union coverage: {union_coverage}")
    report["unionCoverageBp"] = union_coverage
    assert_static_svg_safety(report)


def _assert_gene_labels_only(report: dict[str, Any]) -> None:
    human = SeqIO.read(FIRST_CIRCULAR_FIXTURE_PATH, "genbank")
    gene_labels = {
        str(feature.qualifiers["gene"][0])
        for feature in human.features
        if feature.type == "CDS" and feature.qualifiers.get("gene")
    }
    product_labels = {
        str(feature.qualifiers["product"][0])
        for feature in human.features
        if feature.type == "CDS" and feature.qualifiers.get("product")
    }
    texts = set(report.get("texts", []))
    missing = gene_labels - texts
    if missing:
        raise AssertionError(f"Human mitochondrial gene labels are missing: {sorted(missing)}")
    leaked_products = product_labels & texts
    if leaked_products:
        raise AssertionError(
            f"CDS product labels were drawn instead of gene labels: {sorted(leaked_products)}"
        )
    report["cdsGeneLabels"] = sorted(gene_labels)
    report["cdsProductLabels"] = []


def _assert_circular_download(path: Path, expected_ring_count: int) -> dict[str, Any]:
    contents = path.read_bytes()
    if len(contents) < 20_000:
        raise AssertionError(f"Downloaded Circular SVG is too small: {len(contents)}")
    root = ET.fromstring(contents)
    elements = list(root.iter())
    if any(element.tag.rsplit("}", 1)[-1] == "script" for element in elements):
        raise AssertionError("Static Circular SVG contains a script")
    rings = [
        element
        for element in elements
        if element.attrib.get("data-gbdraw-slot-renderer")
        == "sequence_conservation"
    ]
    if len(rings) != expected_ring_count:
        raise AssertionError(
            f"Expected {expected_ring_count} downloaded rings, found {len(rings)}"
        )
    hits = [
        element
        for element in elements
        if element.attrib.get("data-match-kind") == "homology"
        and element.attrib.get("data-gbdraw-match-id")
    ]
    if len(hits) < expected_ring_count:
        raise AssertionError("Downloaded Circular SVG lost comparison HSPs")
    if {hit.attrib.get("data-reference-record-id") for hit in hits} != {
        "NC_012920.1"
    }:
        raise AssertionError("Downloaded ring reference metadata changed")
    return {
        "filename": path.name,
        "bytes": len(contents),
        "rings": len(rings),
        "hits": len(hits),
    }


def _assert_span_fasta(path: Path) -> dict[str, Any]:
    records = list(SeqIO.parse(path, "fasta"))
    if len(records) != 2:
        raise AssertionError(f"Expected two HSP span FASTA records, found {len(records)}")
    lengths = [len(record) for record in records]
    if not all(length > 0 for length in lengths):
        raise AssertionError(f"Downloaded HSP FASTA contains an empty span: {lengths}")
    ids = {record.id for record in records}
    if not any("NC_012920.1" in record_id for record_id in ids):
        raise AssertionError(f"HSP FASTA lacks the human reference span: {sorted(ids)}")
    return {
        "filename": path.name,
        "bytes": path.stat().st_size,
        "records": len(records),
        "record_ids": sorted(ids),
        "lengths": lengths,
    }


def _fit_circular_ring_preview(page: Page) -> None:
    reset_zoom = page.get_by_role("button", name="Reset zoom", exact=True)
    zoom_out = page.get_by_role("button", name="Zoom out", exact=True)
    for _ in range(8):
        if "50%" in reset_zoom.inner_text():
            break
        zoom_out.click()
    else:
        raise AssertionError("Could not reach the documented 50% Circular zoom")

    search = page.get_by_role("searchbox", name="Search features", exact=True)
    search_box = search.bounding_box()
    if search_box is None:
        raise AssertionError("Could not resolve the feature-search palette")
    drag_x = search_box["x"] - 4
    drag_y = search_box["y"] + search_box["height"] + 3
    page.mouse.move(drag_x, drag_y)
    page.mouse.down()
    page.mouse.move(drag_x, 2, steps=10)
    page.mouse.up()

    result_region = page.get_by_role("region", name="Result Preview", exact=True)
    preview_box = result_region.bounding_box()
    if preview_box is None:
        raise AssertionError("Could not resolve the Circular result preview")
    y = preview_box["y"] + (preview_box["height"] * 0.72)
    page.mouse.move(preview_box["x"] + (preview_box["width"] * 0.70), y)
    page.mouse.down()
    page.mouse.move(
        preview_box["x"] + (preview_box["width"] * 0.22),
        y,
        steps=10,
    )
    page.mouse.up()
    page.mouse.click(
        preview_box["x"] + (preview_box["width"] * 0.82),
        preview_box["y"] + (preview_box["height"] * 0.90),
    )
    page.wait_for_timeout(250)


def capture_gui_circular_rings(
    browser_type: BrowserType,
    base_url: str,
    output_paths: Mapping[str, Path],
    download_dir: Path,
) -> ComparisonCaptureResult:
    """Run H-GUI-06 with four complete, naturally circular mtDNA records."""

    assert_fixture_identity(
        FIRST_CIRCULAR_FIXTURE_PATH,
        expected_size=FIRST_CIRCULAR_FIXTURE_SIZE,
        expected_sha256=FIRST_CIRCULAR_FIXTURE_SHA256,
    )
    assert_fixture_identity(
        FIRST_LINEAR_LABEL_RULE_PATH,
        expected_size=FIRST_LINEAR_LABEL_RULE_SIZE,
        expected_sha256=FIRST_LINEAR_LABEL_RULE_SHA256,
    )
    assert_output_paths(
        output_paths,
        GUI_CIRCULAR_RINGS_SCREENSHOT_NAMES,
        "H-GUI-06",
    )
    download_dir.mkdir(parents=True, exist_ok=True)
    source_records = _prepare_complete_mitochondrial_inputs(download_dir)

    capture = open_browser_capture(browser_type, base_url)
    page = capture.page
    screenshot_bytes: dict[str, int] = {}
    try:
        page.goto(base_url, wait_until="domcontentloaded")
        wait_for_worker(page)
        circular = page.get_by_role("button", name="Circular", exact=True)
        circular.click()
        expect(circular).to_have_attribute("aria-pressed", "true")
        page.get_by_role("radio", name="GenBank", exact=True).check()
        page.get_by_label("GenBank/DDBJ File", exact=True).set_input_files(
            FIRST_CIRCULAR_FIXTURE_PATH
        )
        expect(
            page.get_by_role(
                "group", name="GenBank/DDBJ File selection", exact=True
            )
        ).to_contain_text("HmmtDNA.gbk")
        prefix = page.get_by_label("Output Prefix", exact=True)
        prefix.fill("circular_similarity_rings")
        expect(prefix).to_have_value("circular_similarity_rings")
        species = page.get_by_label("Species", exact=True)
        species.fill("<i>Homo sapiens</i>")
        expect(species).to_have_value("<i>Homo sapiens</i>")
        separate_strands = page.get_by_label("Separate Strands", exact=True)
        separate_strands.uncheck()
        expect(separate_strands).not_to_be_checked()
        track_preset = page.get_by_label("Track Preset", exact=True)
        track_preset.select_option("middle")
        expect(track_preset).to_have_value("middle")

        page.get_by_label("Pairwise Comparisons", exact=True).click()
        run_losat = page.get_by_role("radio", name="Run LOSAT", exact=True).last
        run_losat.check()
        expect(run_losat).to_be_checked()
        tlosatx = page.get_by_role("radio", name="TLOSATX", exact=True).last
        tlosatx.check()
        expect(tlosatx).to_be_checked()
        reference_gencode = page.get_by_label(
            "Circular reference gencode", exact=True
        )
        reference_gencode.fill("2")
        expect(reference_gencode).to_have_value("2")

        comparison_records = source_records[1:]
        genetic_codes = (2, 5, 5)
        for index, (record, label, gencode) in enumerate(
            zip(
                comparison_records,
                CIRCULAR_RING_LABELS,
                genetic_codes,
                strict=True,
            ),
            start=1,
        ):
            with page.expect_file_chooser() as chooser_info:
                page.get_by_role(
                    "button", name=re.compile(r"Add Seq")
                ).last.click()
            chooser_info.value.set_files(record["fasta"])
            ring_label = page.get_by_label(
                f"Comparison ring label {index}", exact=True
            )
            ring_label.fill(label)
            expect(ring_label).to_have_value(label)
            subject_gencode = page.get_by_label(
                f"Comparison subject gencode {index}", exact=True
            )
            subject_gencode.fill(str(gencode))
            expect(subject_gencode).to_have_value(str(gencode))

        for label, value in (
            ("Circular comparison minimum bitscore", "50"),
            ("Circular comparison maximum e-value", "1e-5"),
            ("Circular comparison minimum identity", "40"),
            ("Circular comparison minimum alignment length", "50"),
        ):
            control = page.get_by_label(label, exact=True)
            control.fill(value)
            control.press("Tab")
            expect(control).to_have_value(value)
        ring_width = page.get_by_label("Circular comparison ring width", exact=True)
        ring_width.fill("18")
        ring_width.press("Tab")
        expect(ring_width).to_have_value("18")
        ring_gap = page.get_by_label("Circular comparison ring gap", exact=True)
        ring_gap.fill("4")
        ring_gap.press("Tab")
        expect(ring_gap).to_have_value("4")

        page.get_by_label("Labels", exact=True).click()
        label_mode = page.get_by_label("Label Mode", exact=True)
        label_mode.select_option("out")
        expect(label_mode).to_have_value("out")
        page.get_by_label("Priority File (TSV)", exact=True).set_input_files(
            FIRST_LINEAR_LABEL_RULE_PATH
        )
        page.get_by_label("Title & Legend", exact=True).click()
        title = page.get_by_label("Plot Title", exact=True)
        title.fill(CIRCULAR_RING_TITLE)
        expect(title).to_have_value(CIRCULAR_RING_TITLE)
        page.get_by_label("Plot Title Position", exact=True).select_option("top")
        page.get_by_label("Legend Position", exact=True).select_option("right")
        page.get_by_label("Definition Font Size", exact=True).fill("18")

        page.get_by_label("Comparison ring label 1", exact=True).scroll_into_view_if_needed()
        screenshot_bytes["ring-settings.png"] = capture_screenshot(
            page, output_paths["ring-settings.png"], "Circular"
        )

        final_report = generate_and_inspect(
            page,
            _inspect_circular_rings,
            _assert_circular_ring_svg,
        )
        _assert_gene_labels_only(final_report)
        _fit_circular_ring_preview(page)
        screenshot_bytes["ring-result.png"] = capture_screenshot(
            page, output_paths["ring-result.png"], "Circular"
        )

        with page.expect_download(timeout=ACTION_TIMEOUT_MS) as svg_info:
            page.get_by_role("button", name="SVG", exact=True).click()
        svg_download = svg_info.value
        if svg_download.failure() is not None:
            raise AssertionError(f"Circular ring SVG download failed: {svg_download.failure()}")
        if svg_download.suggested_filename != EXPECTED_CIRCULAR_RING_SVG:
            raise AssertionError(
                f"Unexpected Circular ring SVG filename: {svg_download.suggested_filename}"
            )
        svg_path = download_dir / svg_download.suggested_filename
        svg_download.save_as(svg_path)
        download_report = _assert_circular_download(svg_path, expected_ring_count=3)

        first_match = page.get_by_role("button", name="Pairwise match 1", exact=True)
        expect(first_match).to_be_visible(timeout=WORKER_READY_TIMEOUT_MS)
        first_match.focus()
        first_match.press("Enter")
        popup = page.get_by_role("dialog", name="Pairwise match details", exact=True)
        expect(popup).to_be_visible()
        popup_text = popup.inner_text()
        for expected in ("Homology ring match", "NC_012920.1", "span"):
            if expected not in popup_text:
                raise AssertionError(f"Circular HSP popup is missing {expected!r}")
        fasta_buttons = popup.get_by_role("button", name=re.compile(r"FASTA"))
        expect(fasta_buttons).to_have_count(3)
        with page.expect_download(timeout=ACTION_TIMEOUT_MS) as fasta_info:
            fasta_buttons.last.click()
        fasta_download = fasta_info.value
        if fasta_download.failure() is not None:
            raise AssertionError(f"HSP FASTA download failed: {fasta_download.failure()}")
        fasta_path = download_dir / EXPECTED_CIRCULAR_SPAN_FASTA
        fasta_download.save_as(fasta_path)
        fasta_report = _assert_span_fasta(fasta_path)

        screenshot_bytes["hsp-popup.png"] = capture_screenshot(
            page, output_paths["hsp-popup.png"], "Circular"
        )
        popup_report = {
            "title": "Homology ring match",
            "text": popup_text,
            "download_suggested_filename": fasta_download.suggested_filename,
        }
        capture.assert_clean()
    finally:
        capture.close()

    return ComparisonCaptureResult(
        screenshot_bytes=screenshot_bytes,
        final_svg_semantics=final_report,
        download=download_report,
        fasta_download=fasta_report,
        popup=popup_report,
        source_records=source_records,
    )


__all__ = [
    "capture_gui_circular_rings",
    "capture_gui_tlosatx",
    "capture_gui_uploaded_comparison",
]
