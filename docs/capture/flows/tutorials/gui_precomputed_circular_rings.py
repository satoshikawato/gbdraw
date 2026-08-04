"""Capture the precomputed circular-comparison-rings project tutorial."""

from __future__ import annotations

import csv
import re
from pathlib import Path
from typing import Any, Mapping

from Bio import SeqIO
from playwright.sync_api import BrowserType, expect

from assertions.svg_semantics import assert_static_svg_safety, inspect_first_circular_svg
from config import (
    ACTION_TIMEOUT_MS,
    FIRST_CIRCULAR_FIXTURE_PATH,
    FIRST_CIRCULAR_FIXTURE_SHA256,
    FIRST_CIRCULAR_FIXTURE_SIZE,
    FIRST_LINEAR_LABEL_RULE_PATH,
    FIRST_LINEAR_LABEL_RULE_SHA256,
    FIRST_LINEAR_LABEL_RULE_SIZE,
    GUI_PRECOMPUTED_CIRCULAR_RINGS_SCREENSHOT_NAMES,
    WEB_ROOT,
    WORKER_READY_TIMEOUT_MS,
)
from flows.how_to.nucleotide_comparisons import (
    CIRCULAR_RING_LABELS,
    ComparisonCaptureResult,
    _assert_circular_download,
    _assert_circular_ring_svg,
    _assert_gene_labels_only,
    _assert_span_fasta,
    _fit_circular_ring_preview,
    _inspect_circular_rings,
)
from flows.web_capture import (
    assert_fixture_identity,
    assert_output_paths,
    capture_screenshot,
    generate_and_inspect,
    open_browser_capture,
    wait_for_worker,
)


SCENARIO_ID = "T-GUI-06"
OUTPUT_SVG = "precomputed_circular_rings.svg"
OUTPUT_FASTA = "circular_hsp_spans.fasta"
TITLE = "Precomputed TLOSATX rings around Homo sapiens mtDNA"
COMPARISON_ROOT = WEB_ROOT / "tutorial-data" / "metazoan-mitochondria-comparison"
COMPARISON_FIXTURES = (
    (
        COMPARISON_ROOT / "danio-human.tlosatx.tsv",
        19_284,
        "4135025ba9fc6f346551757ad3d842b30c2f813bb2a412882d571997ff95c597",
        "NC_002333.2",
        276,
        68,
        COMPARISON_ROOT / "NC_002333.2.fna",
        16_886,
        "2fc1138b90e4ae0e38f7ca9d8e0dfa944842ea1c175670dccfe0c292d873ed93",
        16_596,
    ),
    (
        COMPARISON_ROOT / "drosophila-human.tlosatx.tsv",
        6_587,
        "1d9287499564e24bc52063625be977c71976c2e5f2bad236ec08733584358117",
        "NC_024511.2",
        93,
        24,
        COMPARISON_ROOT / "NC_024511.2.fna",
        19_863,
        "cf44d9db8d78344182e49e12e43593f03bd009c2d4bfb1ba0a37d463cf44731e",
        19_524,
    ),
    (
        COMPARISON_ROOT / "caenorhabditis-human.tlosatx.tsv",
        4_681,
        "b7244a0154965935d60694549a8d2e779c8f0788737b2a1f4f1f96ebc57bd416",
        "NC_001328.1",
        66,
        14,
        COMPARISON_ROOT / "NC_001328.1.fna",
        14_037,
        "6d3258a80aa12f399f801e3eb1296821b1eebb40c7ca41342b4aec2b765133ee",
        13_794,
    ),
)


def _validate_comparison_fixtures() -> tuple[dict[str, Any], ...]:
    reports = []
    for (
        table,
        table_size,
        table_sha,
        query_id,
        raw_count,
        retained_count,
        fasta,
        fasta_size,
        fasta_sha,
        sequence_length,
    ) in COMPARISON_FIXTURES:
        assert_fixture_identity(
            table, expected_size=table_size, expected_sha256=table_sha
        )
        with table.open(encoding="utf-8", newline="") as handle:
            rows = [
                row
                for row in csv.reader(handle, delimiter="\t")
                if row and not row[0].startswith("#")
            ]
        if len(rows) != raw_count or any(len(row) != 12 for row in rows):
            raise AssertionError(f"Unexpected TLOSATX rows in {table.name}")
        if {(row[0], row[1]) for row in rows} != {
            (query_id, "NC_012920.1")
        }:
            raise AssertionError(f"Unexpected TLOSATX endpoints in {table.name}")
        retained = [
            row
            for row in rows
            if float(row[11]) >= 50
            and float(row[10]) <= 1e-5
            and float(row[2]) >= 40
            and int(row[3]) >= 50
        ]
        if len(retained) != retained_count:
            raise AssertionError(f"Unexpected retained rows in {table.name}")

        assert_fixture_identity(
            fasta, expected_size=fasta_size, expected_sha256=fasta_sha
        )
        records = list(SeqIO.parse(fasta, "fasta"))
        if (
            len(records) != 1
            or records[0].id != query_id
            or len(records[0]) != sequence_length
        ):
            raise AssertionError(f"Unexpected comparison FASTA: {fasta.name}")
        reports.append(
            {
                "record_id": query_id,
                "length": sequence_length,
                "table": table,
                "fasta": fasta,
                "raw_rows": raw_count,
                "retained_rows": retained_count,
            }
        )
    return tuple(reports)


def _assert_plain_reference(report: Mapping[str, Any]) -> None:
    if report.get("svgCount") != 1:
        raise AssertionError("Step 2 did not render one SVG")
    if "NC_012920.1" not in report.get("recordIds", []):
        raise AssertionError("Step 2 does not identify the human reference")
    if report.get("featureElementCount") != 37:
        raise AssertionError("Step 2 does not contain all 37 feature elements")
    if report.get("comparisonMatches") or report.get("matches"):
        raise AssertionError("Step 2 already contains comparison matches")
    assert_static_svg_safety(dict(report))


def _assert_precomputed_result(report: dict[str, Any]) -> None:
    _assert_circular_ring_svg(report)
    if len(report.get("hits", [])) != 106:
        raise AssertionError(
            f"Expected 106 retained precomputed HSPs, found {len(report.get('hits', []))}"
        )


def capture_gui_precomputed_circular_rings(
    browser_type: BrowserType,
    base_url: str,
    output_paths: Mapping[str, Path],
    download_dir: Path,
) -> ComparisonCaptureResult:
    """Upload three frozen TLOSATX tables and export an inspected ring map."""

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
    comparison_reports = _validate_comparison_fixtures()
    assert_output_paths(
        output_paths,
        GUI_PRECOMPUTED_CIRCULAR_RINGS_SCREENSHOT_NAMES,
        SCENARIO_ID,
    )
    download_dir.mkdir(parents=True, exist_ok=True)

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
        page.get_by_label("Output Prefix", exact=True).fill(
            "precomputed_circular_rings"
        )
        page.get_by_label("Species", exact=True).fill("<i>Homo sapiens</i>")
        page.get_by_label("Separate Strands", exact=True).uncheck()
        screenshot_bytes[GUI_PRECOMPUTED_CIRCULAR_RINGS_SCREENSHOT_NAMES[0]] = (
            capture_screenshot(
                page,
                output_paths[GUI_PRECOMPUTED_CIRCULAR_RINGS_SCREENSHOT_NAMES[0]],
                "Circular",
            )
        )

        generate_and_inspect(
            page, inspect_first_circular_svg, _assert_plain_reference
        )
        _fit_circular_ring_preview(page)
        screenshot_bytes[GUI_PRECOMPUTED_CIRCULAR_RINGS_SCREENSHOT_NAMES[1]] = (
            capture_screenshot(
                page,
                output_paths[GUI_PRECOMPUTED_CIRCULAR_RINGS_SCREENSHOT_NAMES[1]],
                "Circular",
            )
        )

        page.get_by_label("Pairwise Comparisons", exact=True).click()
        upload = page.get_by_role("radio", name="Upload BLAST", exact=True)
        upload.check()
        expect(upload).to_be_checked()
        page.get_by_label("BLAST outfmt 6/7 files", exact=True).set_input_files(
            [report["table"] for report in comparison_reports]
        )
        reference_side = page.get_by_label(
            "Circular comparison reference side", exact=True
        )
        reference_side.select_option("subject")
        expect(reference_side).to_have_value("subject")

        companion_labels = page.locator("label").filter(
            has_text="Comparison FASTA (optional)"
        )
        expect(companion_labels).to_have_count(3)
        for index, (report, label) in enumerate(
            zip(comparison_reports, CIRCULAR_RING_LABELS, strict=True)
        ):
            companion_labels.nth(index).locator("input[type=file]").set_input_files(
                report["fasta"]
            )
            ring_label = page.get_by_label(
                f"Comparison ring label {index + 1}", exact=True
            )
            ring_label.fill(label)
            expect(ring_label).to_have_value(label)

        for label, value in (
            ("Circular comparison minimum bitscore", "50"),
            ("Circular comparison maximum e-value", "1e-5"),
            ("Circular comparison minimum identity", "40"),
            ("Circular comparison minimum alignment length", "50"),
            ("Circular comparison ring width", "18"),
            ("Circular comparison ring gap", "4"),
        ):
            control = page.get_by_label(label, exact=True)
            control.fill(value)
            control.press("Tab")
            expect(control).to_have_value(value)

        page.get_by_label("Labels", exact=True).click()
        page.get_by_label("Label Mode", exact=True).select_option("out")
        page.get_by_label("Priority File (TSV)", exact=True).set_input_files(
            FIRST_LINEAR_LABEL_RULE_PATH
        )
        page.get_by_label("Title & Legend", exact=True).click()
        page.get_by_label("Plot Title", exact=True).fill(TITLE)
        page.get_by_label("Plot Title Position", exact=True).select_option("top")
        page.get_by_label("Legend Position", exact=True).select_option("right")
        page.get_by_label("Definition Font Size", exact=True).fill("18")

        page.get_by_label(
            "Comparison ring label 1", exact=True
        ).scroll_into_view_if_needed()
        screenshot_bytes[GUI_PRECOMPUTED_CIRCULAR_RINGS_SCREENSHOT_NAMES[2]] = (
            capture_screenshot(
                page,
                output_paths[GUI_PRECOMPUTED_CIRCULAR_RINGS_SCREENSHOT_NAMES[2]],
                "Circular",
            )
        )

        final_report = generate_and_inspect(
            page, _inspect_circular_rings, _assert_precomputed_result
        )
        _assert_gene_labels_only(final_report)
        _fit_circular_ring_preview(page)
        screenshot_bytes[GUI_PRECOMPUTED_CIRCULAR_RINGS_SCREENSHOT_NAMES[3]] = (
            capture_screenshot(
                page,
                output_paths[GUI_PRECOMPUTED_CIRCULAR_RINGS_SCREENSHOT_NAMES[3]],
                "Circular",
            )
        )

        with page.expect_download(timeout=ACTION_TIMEOUT_MS) as svg_info:
            page.get_by_role("button", name="SVG", exact=True).click()
        svg_download = svg_info.value
        if svg_download.failure() is not None:
            raise AssertionError(f"Circular ring SVG download failed: {svg_download.failure()}")
        if svg_download.suggested_filename != OUTPUT_SVG:
            raise AssertionError(
                f"Unexpected Circular ring SVG filename: {svg_download.suggested_filename}"
            )
        svg_path = download_dir / OUTPUT_SVG
        svg_download.save_as(svg_path)
        download_report = _assert_circular_download(svg_path, expected_ring_count=3)

        first_match = page.get_by_role("button", name="Pairwise match 1", exact=True)
        expect(first_match).to_be_visible(timeout=WORKER_READY_TIMEOUT_MS)
        first_match.focus()
        first_match.press("Enter")
        popup = page.get_by_role("dialog", name="Pairwise match details", exact=True)
        expect(popup).to_be_visible()
        popup_text = popup.inner_text()
        for expected in ("Homology ring match", "NC_012920.1", "Reference side", "subject"):
            if expected.casefold() not in popup_text.casefold():
                raise AssertionError(f"Circular HSP popup is missing {expected!r}")
        fasta_buttons = popup.get_by_role("button", name=re.compile(r"FASTA"))
        expect(fasta_buttons).to_have_count(3)
        with page.expect_download(timeout=ACTION_TIMEOUT_MS) as fasta_info:
            fasta_buttons.last.click()
        fasta_download = fasta_info.value
        fasta_path = download_dir / OUTPUT_FASTA
        fasta_download.save_as(fasta_path)
        fasta_report = _assert_span_fasta(fasta_path)
        if not any(
            any(record_id in fasta_id for fasta_id in fasta_report["record_ids"])
            for record_id in ("NC_002333.2", "NC_024511.2", "NC_001328.1")
        ):
            raise AssertionError("HSP FASTA lacks its comparison-record span")

        screenshot_bytes[GUI_PRECOMPUTED_CIRCULAR_RINGS_SCREENSHOT_NAMES[4]] = (
            capture_screenshot(
                page,
                output_paths[GUI_PRECOMPUTED_CIRCULAR_RINGS_SCREENSHOT_NAMES[4]],
                "Circular",
            )
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
        source_records=(
            {
                "record_id": "NC_012920.1",
                "length": 16_569,
                "genbank": FIRST_CIRCULAR_FIXTURE_PATH,
            },
            *comparison_reports,
        ),
    )


__all__ = ["capture_gui_precomputed_circular_rings"]
