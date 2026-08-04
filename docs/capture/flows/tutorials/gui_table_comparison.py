"""Capture the GUI variant of the four-record Majanivirus comparison."""

from __future__ import annotations

import hashlib
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import Path
from typing import Any, Mapping

from Bio import SeqIO
from playwright.sync_api import BrowserType, Page, expect

from assertions.svg_semantics import inspect_first_linear_svg, inspect_svg_file
from config import ACTION_TIMEOUT_MS, WEB_ROOT, WORKER_READY_TIMEOUT_MS
from flows.how_to.nucleotide_comparisons import ComparisonCaptureResult
from flows.web_capture import (
    assert_fixture_identity,
    assert_output_paths,
    capture_screenshot,
    fit_complete_linear_preview,
    generate_and_inspect,
    open_browser_capture,
    set_feature_search_visible,
    wait_for_worker,
)


SCENARIO_ID = "T-GUI-11"
SCREENSHOT_NAMES = ("comparison-settings.png", "comparison-result.png")
OUTPUT_PREFIX = "table_driven_comparison"
OUTPUT_FILENAME = f"{OUTPUT_PREFIX}.svg"
FIXTURE_ROOT = WEB_ROOT / "tutorial-data" / "majanivirus-table-comparison"
RECORDS = (
    (
        "MjeNMV.gb",
        527_653,
        "20b7af3666c66be83002c4f8417b2b07664e7d7d2e8ea6a6c7c589bd5157c0ba",
        "LC738868.1",
        306_008,
        "MjeNMV",
        False,
    ),
    (
        "MelaMJNV.gb",
        497_312,
        "e226b4b6651a5e70a5a402c7553a231a082f2714bee92623d78040ba119649a0",
        "LC738874.1",
        287_061,
        "MelaMJNV",
        False,
    ),
    (
        "PemoMJNVA.gb",
        505_839,
        "1f662eb6ba001a6c396adfe48e67119ecb0ba034b3164df17eae844e13df0374",
        "LC738870.1",
        294_144,
        "PemoMJNVA",
        True,
    ),
    (
        "PeseMJNV.gb",
        500_530,
        "8de24deafddc6275cc39556818606e17ae7fff582e697257f9814fe191d6dadb",
        "LC738873.1",
        291_934,
        "PeseMJNV",
        True,
    ),
)
COMPARISONS = (
    (
        "MjeNMV.MelaMJNV.tblastx.out",
        1_752_746,
        "7d464249ef018f632ab891dd34fe3e3bc0c8420216e5a7062508d6ddd7be0502",
        "LC738868.1",
        "LC738874.1",
        23_709,
        80,
    ),
    (
        "PemoMJNVA.PeseMJNV.tblastx.out",
        2_511_804,
        "fa2c9ed87ed04de0c03295c2aaafbf90b6cec979ec846bc53eb1250296e5a6a4",
        "LC738870.1",
        "LC738873.1",
        34_062,
        2,
    ),
)
EXPECTED_RECORD_IDS = tuple(row[3] for row in RECORDS)
EXPECTED_LABELS = tuple(row[5] for row in RECORDS)
EXPECTED_REVERSE = tuple(row[6] for row in RECORDS)
EXPECTED_PAIRS = {
    ("LC738868.1", "LC738874.1"): 80,
    ("LC738870.1", "LC738873.1"): 2,
}


def _validate_fixtures() -> tuple[dict[str, Any], ...]:
    report = []
    for name, size, digest, record_id, length, label, reverse in RECORDS:
        path = FIXTURE_ROOT / name
        assert_fixture_identity(path, expected_size=size, expected_sha256=digest)
        record = SeqIO.read(path, "genbank")
        if record.id != record_id or len(record) != length:
            raise AssertionError(f"Unexpected complete Majanivirus record in {name}")
        if record.annotations.get("topology") != "linear":
            raise AssertionError(f"{record_id} is not a complete linear source record")
        if "complete genome" not in record.description:
            raise AssertionError(f"{record_id} is not described as a complete genome")
        report.append(
            {
                "filename": name,
                "recordId": record_id,
                "length": length,
                "label": label,
                "reverse": reverse,
            }
        )

    for name, size, digest, query, subject, raw_count, retained_count in COMPARISONS:
        path = FIXTURE_ROOT / name
        assert_fixture_identity(path, expected_size=size, expected_sha256=digest)
        rows = [
            line.split("\t")
            for line in path.read_text(encoding="utf-8").splitlines()
            if line and not line.startswith("#")
        ]
        if len(rows) != raw_count or {len(row) for row in rows} != {12}:
            raise AssertionError(f"Unexpected TBLASTX evidence shape in {name}")
        if {(row[0], row[1]) for row in rows} != {(query, subject)}:
            raise AssertionError(f"Unexpected endpoints in {name}")
        retained = [
            row
            for row in rows
            if float(row[2]) >= 97
            and int(row[3]) >= 500
            and float(row[10]) <= 1e-5
        ]
        if len(retained) != retained_count:
            raise AssertionError(f"Unexpected filtered evidence count in {name}")
    return tuple(report)


def _inspect_comparison_svg(result_region: Any) -> dict[str, Any]:
    report = inspect_first_linear_svg(result_region)
    report["matchStyles"] = result_region.evaluate(
        """
        root => Array.from(root.querySelectorAll('[data-gbdraw-pairwise-match-id]'))
          .map((node) => String(node.getAttribute('data-pairwise-match-style') || ''))
        """
    )
    return report


def _assert_comparison_svg(report: Mapping[str, Any]) -> None:
    if report.get("svgCount") != 1:
        raise AssertionError("Expected one Majanivirus comparison SVG")
    if report.get("scriptCount") != 0 or report.get("eventAttributes"):
        raise AssertionError("The static comparison SVG contains active content")
    if report.get("unsafeHrefs"):
        raise AssertionError("The static comparison SVG contains an unsafe reference")
    placements = report.get("recordPlacements", [])
    actual_order = tuple(item.get("recordId") for item in placements)
    if actual_order != EXPECTED_RECORD_IDS:
        raise AssertionError(f"Unexpected Majanivirus record order: {actual_order!r}")
    if report.get("featureElementCount") != 706:
        raise AssertionError(
            f"Expected 706 feature elements, found {report.get('featureElementCount')}"
        )
    texts = set(report.get("texts", []))
    missing = set(EXPECTED_LABELS) - texts
    if missing:
        raise AssertionError(f"Missing Majanivirus record labels: {sorted(missing)}")
    pairs = Counter(
        (match.get("query"), match.get("subject"))
        for match in report.get("matches", [])
    )
    if dict(pairs) != EXPECTED_PAIRS:
        raise AssertionError(f"Unexpected retained comparison pairs: {pairs!r}")
    styles = report.get("matchStyles", [])
    if len(styles) != 82 or set(styles) != {"curve"}:
        raise AssertionError(f"Unexpected comparison-link styles: {Counter(styles)!r}")


def _upload_gap(page: Page, label: str, evidence_path: Path) -> None:
    gap = page.get_by_text(label, exact=True).locator("xpath=..")
    upload = gap.get_by_role("radio", name="Upload BLAST TSV", exact=True)
    upload.check()
    expect(upload).to_be_checked()
    gap.get_by_label("BLAST TSV", exact=True).set_input_files(evidence_path)
    expect(
        gap.get_by_role("group", name="BLAST TSV selection", exact=True)
    ).to_contain_text(evidence_path.name)


def _capture_state(page: Page) -> dict[str, Any]:
    state = page.evaluate(
        """
        () => {
          const app = window.__GBDRAW_APP__;
          return {
            prefix: String(app?.form?.prefix || ''),
            definitions: (app?.linearSeqs || []).map((seq) => String(seq?.definition || '')),
            reverse: (app?.linearSeqs || []).map((seq) => Boolean(seq?.region_reverse)),
            layoutEnabled: Boolean(app?.linearRecordLayoutEnabled),
            gap: Number(app?.linearRecordGap),
            labels: String(app?.form?.show_labels_linear || ''),
            trackLayout: String(app?.form?.linear_track_layout || ''),
            scaleStyle: String(app?.form?.scale_style || ''),
            rulerOnAxis: Boolean(app?.form?.linear_ruler_on_axis),
            lockDefinition: Boolean(app?.form?.keep_definition_left_aligned),
            legend: String(app?.form?.legend || ''),
            matchStyle: String(app?.adv?.pairwise_match_style || ''),
            comparisonHeight: Number(app?.adv?.comparison_height),
            evalue: String(app?.adv?.evalue || ''),
            identity: Number(app?.adv?.identity),
            alignmentLength: Number(app?.adv?.alignment_length)
          };
        }
        """
    )
    state["rows"] = [
        int(
            page.get_by_label(
                f"Linear record row for sequence {index}", exact=True
            ).input_value()
        )
        for index in range(1, 5)
    ]
    return state


def _assert_state(state: Mapping[str, Any]) -> None:
    expected = {
        "prefix": OUTPUT_PREFIX,
        "definitions": list(EXPECTED_LABELS),
        "reverse": list(EXPECTED_REVERSE),
        "layoutEnabled": True,
        "gap": 28,
        "labels": "none",
        "trackLayout": "above",
        "scaleStyle": "ruler",
        "rulerOnAxis": True,
        "lockDefinition": True,
        "legend": "right",
        "matchStyle": "curve",
        "comparisonHeight": 100,
        "evalue": "1e-5",
        "identity": 97,
        "alignmentLength": 500,
        "rows": [1, 2, 3, 4],
    }
    if dict(state) != expected:
        raise AssertionError(f"Unexpected table-comparison control state: {state!r}")


def _download_svg(page: Page, download_dir: Path) -> dict[str, Any]:
    with page.expect_download(timeout=ACTION_TIMEOUT_MS) as info:
        page.get_by_role("button", name="SVG", exact=True).click()
    download = info.value
    if download.failure() is not None:
        raise AssertionError(f"SVG download failed: {download.failure()}")
    if download.suggested_filename != OUTPUT_FILENAME:
        raise AssertionError(
            f"Expected {OUTPUT_FILENAME}, found {download.suggested_filename}"
        )
    path = download_dir / OUTPUT_FILENAME
    download.save_as(path)
    report = inspect_svg_file(path)
    root = ET.parse(path).getroot()
    report["matchStyles"] = [
        element.attrib.get("data-pairwise-match-style", "")
        for element in root.iter()
        if element.attrib.get("data-gbdraw-pairwise-match-id")
    ]
    _assert_comparison_svg(report)
    return {
        "filename": path.name,
        "bytes": path.stat().st_size,
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "path": str(path),
        "semantics": report,
    }


def capture_gui_table_comparison(
    browser_type: BrowserType,
    base_url: str,
    output_paths: Mapping[str, Path],
    download_dir: Path,
) -> ComparisonCaptureResult:
    """Build the CLI/Python-equivalent four-record comparison in the web app."""

    source_records = _validate_fixtures()
    assert_output_paths(output_paths, SCREENSHOT_NAMES, SCENARIO_ID)
    download_dir.mkdir(parents=True, exist_ok=True)
    capture = open_browser_capture(browser_type, base_url)
    page = capture.page
    screenshot_bytes: dict[str, int] = {}
    try:
        page.goto(base_url, wait_until="domcontentloaded")
        wait_for_worker(page)
        page.get_by_role("button", name="Linear", exact=True).click()
        page.get_by_role("radio", name="GenBank", exact=True).check()
        for _ in range(3):
            page.get_by_role("button", name="Add sequence", exact=True).click()
        for index, record in enumerate(RECORDS, start=1):
            name, _, _, record_id, length, label, reverse = record
            page.get_by_test_id(f"linear-genbank-{index}").set_input_files(
                FIXTURE_ROOT / name
            )
            selector = page.get_by_label(
                f"Record selector for sequence {index}", exact=True
            )
            expect(selector).to_contain_text(
                record_id, timeout=WORKER_READY_TIMEOUT_MS
            )
            definition = page.get_by_label(
                f"Definition for sequence {index}", exact=True
            )
            definition.fill(label)
            page.get_by_label(
                f"Region start for sequence {index}", exact=True
            ).fill("1")
            page.get_by_label(
                f"Region end for sequence {index}", exact=True
            ).fill(str(length))
            reverse_control = page.get_by_label(
                f"Reverse complement for sequence {index}", exact=True
            )
            if reverse:
                reverse_control.check()
            else:
                reverse_control.uncheck()

        arrange = page.get_by_label("Arrange linear records in rows", exact=True)
        arrange.check()
        for index in range(1, 5):
            row = page.get_by_label(
                f"Linear record row for sequence {index}", exact=True
            )
            row.fill(str(index))
            row.press("Tab")
        page.get_by_label("Linear record gap", exact=True).fill("28")

        _upload_gap(
            page,
            "#1 → #2",
            FIXTURE_ROOT / "MjeNMV.MelaMJNV.tblastx.out",
        )
        _upload_gap(
            page,
            "#3 → #4",
            FIXTURE_ROOT / "PemoMJNVA.PeseMJNV.tblastx.out",
        )
        middle_gap = page.get_by_text("#2 → #3", exact=True).locator("xpath=..")
        middle_gap.get_by_role(
            "radio", name="No comparison", exact=True
        ).check()

        pairwise = page.get_by_label("Pairwise Match", exact=True)
        pairwise.click()
        page.get_by_label("Pairwise Match Style", exact=True).select_option("curve")
        page.get_by_label("Pairwise Match Height", exact=True).fill("100")
        page.get_by_label(
            "Linear comparison maximum e-value", exact=True
        ).fill("1e-5")
        page.get_by_label(
            "Linear comparison minimum identity", exact=True
        ).fill("97")
        page.get_by_label(
            "Linear comparison minimum alignment length", exact=True
        ).fill("500")
        pairwise.click()

        page.get_by_label("Track Layout", exact=True).select_option("above")
        lock_definition = page.get_by_role(
            "checkbox", name="Lock Definition Column", exact=True
        )
        lock_definition.check()
        axis = page.get_by_label("Axis & Scale", exact=True)
        axis.click()
        page.get_by_label("Show Coordinate Scale (Linear)", exact=True).check()
        page.get_by_label("Linear scale style", exact=True).select_option("ruler")
        page.get_by_label("Ruler on Axis", exact=True).check()
        axis.click()

        title = page.get_by_label("Title & Legend", exact=True)
        title.click()
        page.get_by_label("Legend Position", exact=True).select_option("right")
        title.click()
        page.get_by_label("Output Prefix", exact=True).fill(OUTPUT_PREFIX)

        state = _capture_state(page)
        _assert_state(state)
        page.get_by_text("#3 → #4", exact=True).scroll_into_view_if_needed()
        screenshot_bytes[SCREENSHOT_NAMES[0]] = capture_screenshot(
            page, output_paths[SCREENSHOT_NAMES[0]], "Linear"
        )

        final_report = generate_and_inspect(
            page, _inspect_comparison_svg, _assert_comparison_svg
        )
        final_report["state"] = state
        fit_complete_linear_preview(page, target_zoom="40%")
        set_feature_search_visible(page, visible=False)
        screenshot_bytes[SCREENSHOT_NAMES[1]] = capture_screenshot(
            page, output_paths[SCREENSHOT_NAMES[1]], "Linear"
        )
        download_report = _download_svg(page, download_dir)
        capture.assert_clean()
    finally:
        capture.close()

    return ComparisonCaptureResult(
        screenshot_bytes=screenshot_bytes,
        final_svg_semantics=final_report,
        download=download_report,
        source_records=source_records,
    )


__all__ = ["capture_gui_table_comparison"]
