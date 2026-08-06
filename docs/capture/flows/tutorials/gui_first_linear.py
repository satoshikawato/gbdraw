"""Complete T-GUI-02 capture through the public web controls."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

from playwright.sync_api import BrowserType, expect

from assertions.downloads import (
    EXPECTED_FIRST_LINEAR_SVG,
    assert_first_linear_svg_download,
)
from assertions.svg_semantics import (
    assert_finished_linear_svg,
    assert_first_linear_svg,
    inspect_first_linear_svg,
)
from config import (
    ACTION_TIMEOUT_MS,
    FIRST_LINEAR_FIXTURE_PATH,
    FIRST_LINEAR_FIXTURE_SHA256,
    FIRST_LINEAR_FIXTURE_SIZE,
    FIRST_LINEAR_LABEL_RULE_PATH,
    FIRST_LINEAR_LABEL_RULE_SHA256,
    FIRST_LINEAR_LABEL_RULE_SIZE,
    FIRST_LINEAR_SCREENSHOT_NAMES,
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


def capture_first_linear(
    browser_type: BrowserType,
    base_url: str,
    output_paths: Mapping[str, Path],
    download_dir: Path,
) -> CaptureResult:
    """Run T-GUI-02 in one clean context and capture its four images."""

    assert_fixture_identity(
        FIRST_LINEAR_FIXTURE_PATH,
        expected_size=FIRST_LINEAR_FIXTURE_SIZE,
        expected_sha256=FIRST_LINEAR_FIXTURE_SHA256,
    )
    assert_fixture_identity(
        FIRST_LINEAR_LABEL_RULE_PATH,
        expected_size=FIRST_LINEAR_LABEL_RULE_SIZE,
        expected_sha256=FIRST_LINEAR_LABEL_RULE_SHA256,
    )
    assert_output_paths(output_paths, FIRST_LINEAR_SCREENSHOT_NAMES, "T-GUI-02")
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

        genbank = page.get_by_role("radio", name="GenBank", exact=True)
        genbank.check()
        expect(genbank).to_be_checked()
        no_comparison = page.get_by_role("radio", name="No comparison", exact=True)
        no_comparison.check()
        expect(no_comparison).to_be_checked()

        page.get_by_label("GenBank File", exact=True).set_input_files(
            FIRST_LINEAR_FIXTURE_PATH
        )
        selected_file = page.get_by_role(
            "group", name="GenBank File selection", exact=True
        )
        expect(selected_file).to_contain_text("NC_001416.gb")
        screenshot_bytes["01-input-ready.png"] = capture_screenshot(
            page, output_paths["01-input-ready.png"], "Linear"
        )

        first_report = generate_and_inspect(
            page, inspect_first_linear_svg, assert_first_linear_svg
        )
        fit_complete_linear_preview(page, target_zoom="30%", pan_left=True)
        screenshot_bytes["02-first-diagram.png"] = capture_screenshot(
            page, output_paths["02-first-diagram.png"], "Linear"
        )

        output_prefix = page.get_by_label("Output Prefix", exact=True)
        output_prefix.fill("lambda_linear")
        expect(output_prefix).to_have_value("lambda_linear")

        track_layout = page.get_by_label("Track Layout", exact=True)
        track_layout.select_option("middle")
        expect(track_layout).to_have_value("middle")
        separate_strands = page.get_by_label("Separate Strands", exact=True)
        separate_strands.check()
        expect(separate_strands).to_be_checked()

        labels_panel = page.get_by_label("Labels", exact=True)
        labels_panel.click()
        show_labels = page.get_by_label("Show Labels", exact=True)
        show_labels.select_option("all")
        expect(show_labels).to_have_value("all")
        page.get_by_label("Priority File (TSV)", exact=True).set_input_files(
            FIRST_LINEAR_LABEL_RULE_PATH
        )
        selected_priority = page.get_by_role(
            "group", name="Priority File (TSV) selection", exact=True
        )
        expect(selected_priority).to_contain_text("cds_gene_qualifier_priority.tsv")
        labels_panel.click()

        title_legend_panel = page.get_by_label("Title & Legend", exact=True)
        title_legend_panel.click()
        legend_position = page.get_by_label("Legend Position", exact=True)
        legend_position.select_option("left")
        expect(legend_position).to_have_value("left")
        title_legend_panel.click()

        axis_panel = page.get_by_label("Axis & Scale", exact=True)
        axis_panel.click()
        show_scale = page.get_by_label("Show Coordinate Scale (Linear)", exact=True)
        show_scale.check()
        expect(show_scale).to_be_checked()
        scale_style = page.get_by_label("Linear scale style", exact=True)
        scale_style.select_option("ruler")
        expect(scale_style).to_have_value("ruler")
        scale_style.scroll_into_view_if_needed()
        screenshot_bytes["03-layout-settings.png"] = capture_screenshot(
            page, output_paths["03-layout-settings.png"], "Linear"
        )

        final_report = generate_and_inspect(
            page, inspect_first_linear_svg, assert_first_linear_svg
        )
        assert_finished_linear_svg(final_report)
        fit_complete_linear_preview(page, target_zoom="30%", pan_left=True)
        screenshot_bytes["04-finished-diagram.png"] = capture_screenshot(
            page, output_paths["04-finished-diagram.png"], "Linear"
        )

        svg_button = page.get_by_role("button", name="SVG", exact=True)
        expect(svg_button).to_be_visible()
        with page.expect_download(timeout=ACTION_TIMEOUT_MS) as download_info:
            svg_button.click()
        download = download_info.value
        if download.failure() is not None:
            raise AssertionError(f"SVG download failed: {download.failure()}")
        if download.suggested_filename != EXPECTED_FIRST_LINEAR_SVG:
            raise AssertionError(
                "Unexpected SVG download name: "
                f"expected {EXPECTED_FIRST_LINEAR_SVG}, "
                f"found {download.suggested_filename}"
            )
        downloaded_svg = download_dir / download.suggested_filename
        download.save_as(downloaded_svg)
        download_report = assert_first_linear_svg_download(downloaded_svg)

        capture.assert_clean()
    finally:
        capture.close()

    return CaptureResult(
        screenshot_bytes=screenshot_bytes,
        first_svg_semantics=first_report,
        final_svg_semantics=final_report,
        download=download_report,
    )
