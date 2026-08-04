"""Complete T-GUI-03 with a real serial browser LOSATN run."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

from playwright.sync_api import BrowserType, expect

from assertions.downloads import (
    EXPECTED_GUI_LOSATN_SVG,
    EXPECTED_GUI_LOSATN_TSV,
    assert_gui_losatn_svg_download,
    assert_gui_losatn_tsv_download,
)
from assertions.svg_semantics import (
    assert_gui_losatn_svg,
    assert_plain_gui_losatn_svg,
    inspect_gui_losatn_svg,
)
from config import (
    ACTION_TIMEOUT_MS,
    FIRST_LINEAR_FIXTURE_PATH,
    FIRST_LINEAR_FIXTURE_SHA256,
    FIRST_LINEAR_FIXTURE_SIZE,
    GUI_LOSATN_DE3_FIXTURE_PATH,
    GUI_LOSATN_DE3_FIXTURE_SHA256,
    GUI_LOSATN_DE3_FIXTURE_SIZE,
    GUI_LOSATN_REFERENCE_TSV_PATH,
    GUI_LOSATN_REFERENCE_TSV_SHA256,
    GUI_LOSATN_REFERENCE_TSV_SIZE,
    GUI_LOSATN_SCREENSHOT_NAMES,
)
from flows.web_capture import (
    assert_fixture_identity,
    assert_output_paths,
    capture_screenshot,
    fit_complete_linear_preview,
    generate_and_inspect,
    open_browser_capture,
    wait_for_worker,
)


@dataclass(frozen=True)
class CaptureResult:
    screenshot_bytes: dict[str, int]
    first_svg_semantics: dict[str, Any]
    final_svg_semantics: dict[str, Any]
    download: dict[str, Any]
    tsv_download: dict[str, Any]
    popup: dict[str, str]


def _assert_match_popup(page: Any) -> dict[str, str]:
    popup = page.get_by_role("dialog", name="Pairwise match details", exact=True)
    expect(popup).to_be_visible()
    popup_text = popup.inner_text()
    for expected in (
        "Pairwise match",
        "NC_001416.1",
        "NC_042057.1",
        "1..21231",
        "20081..41311",
        "99.981",
        "21232",
        "39185",
    ):
        if expected not in popup_text:
            raise AssertionError(f"LOSATN match popup is missing {expected!r}")
    return {"title": "Pairwise match", "text": popup_text}


def capture_gui_losatn(
    browser_type: BrowserType,
    base_url: str,
    output_paths: Mapping[str, Path],
    download_dir: Path,
) -> CaptureResult:
    """Capture the five approved T-GUI-03 screenshots and both downloads."""

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
        (
            GUI_LOSATN_REFERENCE_TSV_PATH,
            GUI_LOSATN_REFERENCE_TSV_SIZE,
            GUI_LOSATN_REFERENCE_TSV_SHA256,
        ),
    ):
        assert_fixture_identity(path, expected_size=size, expected_sha256=digest)
    assert_output_paths(output_paths, GUI_LOSATN_SCREENSHOT_NAMES, "T-GUI-03")
    download_dir.mkdir(parents=True, exist_ok=True)

    capture = open_browser_capture(browser_type, base_url)
    page = capture.page
    screenshot_bytes: dict[str, int] = {}

    try:
        page.goto(base_url, wait_until="domcontentloaded")
        wait_for_worker(page)

        linear = page.get_by_role("button", name="Linear", exact=True)
        linear.click()
        expect(linear).to_have_attribute("aria-pressed", "true")
        page.get_by_role("radio", name="GenBank", exact=True).check()
        no_comparison = page.get_by_role(
            "radio", name="No comparison", exact=True
        ).first
        no_comparison.check()
        expect(no_comparison).to_be_checked()

        page.get_by_role("button", name="Add sequence", exact=True).click()
        page.get_by_test_id("linear-genbank-1").set_input_files(
            FIRST_LINEAR_FIXTURE_PATH
        )
        page.get_by_test_id("linear-genbank-2").set_input_files(
            GUI_LOSATN_DE3_FIXTURE_PATH
        )
        selected_files = page.get_by_role(
            "group", name="GenBank File selection", exact=True
        )
        expect(selected_files).to_have_count(2)
        expect(selected_files.nth(0)).to_contain_text("NC_001416.gb")
        expect(selected_files.nth(1)).to_contain_text("NC_042057.1.gb")
        screenshot_bytes["01-input-ready.png"] = capture_screenshot(
            page, output_paths["01-input-ready.png"], "Linear"
        )

        first_report = generate_and_inspect(
            page, inspect_gui_losatn_svg, assert_plain_gui_losatn_svg
        )
        fit_complete_linear_preview(page)
        screenshot_bytes["02-first-diagram.png"] = capture_screenshot(
            page, output_paths["02-first-diagram.png"], "Linear"
        )

        run_losat = page.get_by_role("radio", name="Run LOSAT", exact=True).first
        run_losat.check()
        expect(run_losat).to_be_checked()
        losatn = page.get_by_role("radio", name="LOSATN", exact=True).first
        losatn.check()
        expect(losatn).to_be_checked()

        execution = page.get_by_role(
            "combobox", name="LOSAT execution", exact=True
        )
        execution.select_option("serial")
        expect(execution).to_have_value("serial")
        total_threads = page.get_by_role(
            "combobox", name="LOSAT total threads", exact=True
        )
        total_threads.select_option("1")
        expect(total_threads).to_have_value("1")
        parallel_runs = page.get_by_role(
            "combobox", name="LOSAT parallel runs", exact=True
        )
        parallel_runs.select_option("1")
        expect(parallel_runs).to_have_value("1")
        threads_per_run = page.get_by_role(
            "combobox", name="LOSAT threads per run", exact=True
        )
        expect(threads_per_run).to_be_disabled()
        task = page.get_by_role("combobox", name="LOSATN task", exact=True)
        task.select_option("megablast")
        expect(task).to_have_value("megablast")

        output_prefix = page.get_by_role(
            "textbox", name="Output Prefix", exact=True
        )
        output_prefix.fill("lambda-de3-losatn")
        expect(output_prefix).to_have_value("lambda-de3-losatn")
        raw_filename = page.get_by_role(
            "textbox", name="Raw LOSAT filename", exact=True
        )
        raw_filename.fill(EXPECTED_GUI_LOSATN_TSV)
        raw_filename.press("Tab")
        expect(raw_filename).to_have_value(EXPECTED_GUI_LOSATN_TSV)
        pairwise_match = page.get_by_label("Pairwise Match", exact=True)
        pairwise_match.click()
        pairwise_match_height = page.get_by_label(
            "Pairwise Match Height", exact=True
        )
        pairwise_match_height.fill("120")
        pairwise_match_height.press("Tab")
        expect(pairwise_match_height).to_have_value("120")
        pairwise_match.click()
        task.scroll_into_view_if_needed()
        screenshot_bytes["03-losatn-settings.png"] = capture_screenshot(
            page, output_paths["03-losatn-settings.png"], "Linear"
        )

        final_report = generate_and_inspect(
            page, inspect_gui_losatn_svg, assert_gui_losatn_svg
        )
        fit_complete_linear_preview(page)
        screenshot_bytes["04-comparison-result.png"] = capture_screenshot(
            page, output_paths["04-comparison-result.png"], "Linear"
        )

        tsv_button = page.get_by_role(
            "button", name="Save Raw LOSAT TSV", exact=True
        ).first
        expect(tsv_button).to_be_enabled()
        with page.expect_download(timeout=ACTION_TIMEOUT_MS) as tsv_download_info:
            tsv_button.click()
        tsv_download = tsv_download_info.value
        if tsv_download.failure() is not None:
            raise AssertionError(f"LOSATN TSV download failed: {tsv_download.failure()}")
        if tsv_download.suggested_filename != EXPECTED_GUI_LOSATN_TSV:
            raise AssertionError(
                "Unexpected LOSATN TSV download name: "
                f"{tsv_download.suggested_filename}"
            )
        downloaded_tsv = download_dir / tsv_download.suggested_filename
        tsv_download.save_as(downloaded_tsv)
        tsv_report = assert_gui_losatn_tsv_download(
            downloaded_tsv, GUI_LOSATN_REFERENCE_TSV_PATH
        )

        svg_button = page.get_by_role("button", name="SVG", exact=True)
        with page.expect_download(timeout=ACTION_TIMEOUT_MS) as svg_download_info:
            svg_button.click()
        svg_download = svg_download_info.value
        if svg_download.failure() is not None:
            raise AssertionError(f"LOSATN SVG download failed: {svg_download.failure()}")
        if svg_download.suggested_filename != EXPECTED_GUI_LOSATN_SVG:
            raise AssertionError(
                "Unexpected LOSATN SVG download name: "
                f"{svg_download.suggested_filename}"
            )
        downloaded_svg = download_dir / svg_download.suggested_filename
        svg_download.save_as(downloaded_svg)
        svg_report = assert_gui_losatn_svg_download(downloaded_svg)

        first_match = page.get_by_role(
            "button", name="Pairwise match 1", exact=True
        )
        expect(first_match).to_be_visible()
        first_match.focus()
        first_match.press("Enter")
        popup_report = _assert_match_popup(page)
        popup = page.get_by_role("dialog", name="Pairwise match details", exact=True)
        popup.hover()
        page.mouse.wheel(0, 420)
        screenshot_bytes["05-match-popup.png"] = capture_screenshot(
            page, output_paths["05-match-popup.png"], "Linear"
        )

        capture.assert_clean()
    finally:
        capture.close()

    return CaptureResult(
        screenshot_bytes=screenshot_bytes,
        first_svg_semantics=first_report,
        final_svg_semantics=final_report,
        download=svg_report,
        tsv_download=tsv_report,
        popup=popup_report,
    )
