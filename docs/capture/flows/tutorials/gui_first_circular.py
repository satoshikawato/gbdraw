"""Complete T-GUI-01 capture through the public web controls."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

from playwright.sync_api import BrowserType, Page, expect

from assertions.downloads import (
    EXPECTED_FIRST_CIRCULAR_SVG,
    assert_first_circular_svg_download,
)
from assertions.svg_semantics import (
    assert_finished_circular_svg,
    assert_first_circular_svg,
    assert_species_markup,
    inspect_first_circular_svg,
)
from config import (
    ACTION_TIMEOUT_MS,
    FIRST_CIRCULAR_SCREENSHOT_NAMES,
    FIRST_CIRCULAR_FIXTURE_PATH,
    FIRST_CIRCULAR_FIXTURE_SHA256,
    FIRST_CIRCULAR_FIXTURE_SIZE,
    FIRST_LINEAR_LABEL_RULE_PATH,
    FIRST_LINEAR_LABEL_RULE_SHA256,
    FIRST_LINEAR_LABEL_RULE_SIZE,
)
from flows.web_capture import (
    assert_fixture_identity,
    assert_output_paths,
    capture_screenshot,
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


def _pan_finished_preview_left(page: Page) -> None:
    """Use the preview's drag interaction to bring the right legend into view."""

    result_region = page.get_by_role("region", name="Result Preview", exact=True)
    box = result_region.bounding_box()
    if box is None:
        raise AssertionError("Could not resolve the Result Preview bounds for panning")
    y = box["y"] + (box["height"] * 0.55)
    page.mouse.move(box["x"] + (box["width"] * 0.43), y)
    page.mouse.down()
    page.mouse.move(box["x"] + (box["width"] * 0.12), y, steps=10)
    page.mouse.up()


def capture_first_circular(
    browser_type: BrowserType,
    base_url: str,
    output_paths: Mapping[str, Path],
    download_dir: Path,
) -> CaptureResult:
    """Run T-GUI-01 Steps 1-5 in one clean context and capture every image."""

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
    assert_output_paths(output_paths, FIRST_CIRCULAR_SCREENSHOT_NAMES, "T-GUI-01")
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

        genbank = page.get_by_role("radio", name="GenBank", exact=True)
        genbank.check()
        expect(genbank).to_be_checked()

        page.get_by_label("GenBank/DDBJ File", exact=True).set_input_files(
            FIRST_CIRCULAR_FIXTURE_PATH
        )
        selected_file = page.get_by_role(
            "group", name="GenBank/DDBJ File selection", exact=True
        )
        expect(selected_file).to_contain_text("HmmtDNA.gbk")
        screenshot_bytes["01-input-ready.png"] = capture_screenshot(
            page, output_paths["01-input-ready.png"], "Circular"
        )

        first_report = generate_and_inspect(
            page, inspect_first_circular_svg, assert_first_circular_svg
        )
        screenshot_bytes["02-first-diagram.png"] = capture_screenshot(
            page, output_paths["02-first-diagram.png"], "Circular"
        )

        output_prefix = page.get_by_label("Output Prefix", exact=True)
        output_prefix.fill("human_mitochondrion")
        expect(output_prefix).to_have_value("human_mitochondrion")
        species = page.get_by_label("Species", exact=True)
        species.fill("<i>Homo sapiens</i>")
        expect(species).to_have_value("<i>Homo sapiens</i>")

        publication_report = generate_and_inspect(
            page, inspect_first_circular_svg, assert_first_circular_svg
        )
        assert_species_markup(publication_report)
        screenshot_bytes["03-publication-label.png"] = capture_screenshot(
            page, output_paths["03-publication-label.png"], "Circular"
        )

        track_preset = page.get_by_label("Track Preset", exact=True)
        track_preset.select_option("middle")
        expect(track_preset).to_have_value("middle")
        separate_strands = page.get_by_label("Separate Strands", exact=True)
        separate_strands.check()
        expect(separate_strands).to_be_checked()
        hide_gc_content = page.get_by_label("Hide GC Content", exact=True)
        hide_gc_content.uncheck()
        expect(hide_gc_content).not_to_be_checked()
        hide_gc_skew = page.get_by_label("Hide GC Skew", exact=True)
        hide_gc_skew.uncheck()
        expect(hide_gc_skew).not_to_be_checked()

        page.get_by_label("Title & Legend", exact=True).click()
        legend_position = page.get_by_label("Legend Position", exact=True)
        legend_position.select_option("right")
        expect(legend_position).to_have_value("right")
        page.get_by_label("Labels", exact=True).click()
        label_mode = page.get_by_label("Label Mode", exact=True)
        label_mode.select_option("out")
        expect(label_mode).to_have_value("out")
        page.get_by_label("Priority File (TSV)", exact=True).set_input_files(
            FIRST_LINEAR_LABEL_RULE_PATH
        )
        expect(
            page.get_by_role(
                "group", name="Priority File (TSV) selection", exact=True
            )
        ).to_contain_text("cds_gene_qualifier_priority.tsv")

        track_preset.scroll_into_view_if_needed()
        expect(track_preset).to_be_visible()
        screenshot_bytes["04-layout-settings.png"] = capture_screenshot(
            page, output_paths["04-layout-settings.png"], "Circular"
        )

        final_report = generate_and_inspect(
            page, inspect_first_circular_svg, assert_first_circular_svg
        )
        assert_finished_circular_svg(final_report)
        zoom_out = page.get_by_role("button", name="Zoom out", exact=True)
        for _ in range(4):
            zoom_out.click()
        expect(page.get_by_role("button", name="Reset zoom", exact=True)).to_contain_text(
            "60%"
        )
        _pan_finished_preview_left(page)
        screenshot_bytes["04-finished-diagram.png"] = capture_screenshot(
            page, output_paths["04-finished-diagram.png"], "Circular"
        )

        svg_button = page.get_by_role("button", name="SVG", exact=True)
        expect(svg_button).to_be_visible()
        svg_button.hover()
        screenshot_bytes["05-export-svg.png"] = capture_screenshot(
            page, output_paths["05-export-svg.png"], "Circular"
        )
        with page.expect_download(timeout=ACTION_TIMEOUT_MS) as download_info:
            svg_button.click()
        download = download_info.value
        if download.failure() is not None:
            raise AssertionError(f"SVG download failed: {download.failure()}")
        if download.suggested_filename != EXPECTED_FIRST_CIRCULAR_SVG:
            raise AssertionError(
                "Unexpected SVG download name: "
                f"expected {EXPECTED_FIRST_CIRCULAR_SVG}, "
                f"found {download.suggested_filename}"
            )
        downloaded_svg = download_dir / download.suggested_filename
        download.save_as(downloaded_svg)
        download_report = assert_first_circular_svg_download(downloaded_svg)
        download_report["sharedFigureSemanticsParity"] = True

        capture.assert_clean()
    finally:
        capture.close()

    return CaptureResult(
        screenshot_bytes=screenshot_bytes,
        first_svg_semantics=first_report,
        final_svg_semantics=final_report,
        download=download_report,
    )
