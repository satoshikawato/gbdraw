"""Complete H-GUI-01 through the public web input controls."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

from playwright.sync_api import BrowserType, expect

from assertions.downloads import (
    EXPECTED_GUI_INPUTS_SVG,
    assert_gui_inputs_svg_download,
)
from assertions.svg_semantics import (
    assert_lambda_input_svg,
    assert_matching_lambda_input_semantics,
    inspect_lambda_input_svg,
)
from config import (
    ACTION_TIMEOUT_MS,
    FIRST_LINEAR_FIXTURE_PATH,
    FIRST_LINEAR_FIXTURE_SHA256,
    FIRST_LINEAR_FIXTURE_SIZE,
    GENERATION_TIMEOUT_MS,
    GUI_INPUTS_FASTA_PATH,
    GUI_INPUTS_FASTA_SHA256,
    GUI_INPUTS_FASTA_SIZE,
    GUI_INPUTS_GFF3_PATH,
    GUI_INPUTS_GFF3_SHA256,
    GUI_INPUTS_GFF3_SIZE,
    GUI_INPUTS_SCREENSHOT_NAMES,
)
from flows.web_capture import (
    assert_fixture_identity,
    assert_output_paths,
    capture_screenshot,
    fit_complete_linear_preview,
    generate_and_inspect,
    generate_and_wait_for_result,
    open_browser_capture,
    wait_for_app_shell,
)


EXPECTED_ID_MISMATCH = (
    "No matching FASTA record found for GFF record NC_001416.1. "
    "Please ensure that all GFF records have corresponding FASTA entries."
)


@dataclass(frozen=True)
class CaptureResult:
    screenshot_bytes: dict[str, int]
    first_svg_semantics: dict[str, Any]
    final_svg_semantics: dict[str, Any]
    download: dict[str, Any]
    validation_message: str


def _write_mismatched_fasta(path: Path) -> None:
    """Write a whole-sequence temporary FASTA with only its record ID changed."""

    source = GUI_INPUTS_FASTA_PATH.read_bytes()
    header, separator, sequence = source.partition(b"\n")
    expected_id = b">NC_001416.1"
    if not separator or not header.startswith(expected_id):
        raise AssertionError("The frozen Lambda FASTA header changed")
    mismatched = b">MISMATCHED_ID" + header[len(expected_id) :] + separator + sequence
    path.write_bytes(mismatched)
    _, written_separator, written_sequence = path.read_bytes().partition(b"\n")
    if written_separator != separator or written_sequence != sequence:
        raise AssertionError("The mismatch fixture did not preserve the whole sequence")


def capture_gui_inputs(
    browser_type: BrowserType,
    base_url: str,
    output_paths: Mapping[str, Path],
    download_dir: Path,
) -> CaptureResult:
    """Run H-GUI-01 in one clean context and capture its three images."""

    for path, size, digest in (
        (
            FIRST_LINEAR_FIXTURE_PATH,
            FIRST_LINEAR_FIXTURE_SIZE,
            FIRST_LINEAR_FIXTURE_SHA256,
        ),
        (GUI_INPUTS_GFF3_PATH, GUI_INPUTS_GFF3_SIZE, GUI_INPUTS_GFF3_SHA256),
        (GUI_INPUTS_FASTA_PATH, GUI_INPUTS_FASTA_SIZE, GUI_INPUTS_FASTA_SHA256),
    ):
        assert_fixture_identity(path, expected_size=size, expected_sha256=digest)
    assert_output_paths(output_paths, GUI_INPUTS_SCREENSHOT_NAMES, "H-GUI-01")
    download_dir.mkdir(parents=True, exist_ok=True)
    mismatched_fasta_path = download_dir / "mismatched_sequence_id.fna"
    _write_mismatched_fasta(mismatched_fasta_path)

    capture = open_browser_capture(browser_type, base_url)
    page = capture.page
    screenshot_bytes: dict[str, int] = {}

    try:
        page.goto(base_url, wait_until="domcontentloaded")
        wait_for_app_shell(page)

        linear = page.get_by_role("button", name="Linear", exact=True)
        linear.click()
        expect(linear).to_have_attribute("aria-pressed", "true")
        expect(page.get_by_role("status").filter(has_text="Current:")).to_contain_text(
            "Current: No comparison"
        )

        genbank = page.get_by_role("radio", name="GenBank", exact=True)
        genbank.check()
        page.get_by_label("GenBank File", exact=True).set_input_files(
            FIRST_LINEAR_FIXTURE_PATH
        )
        expect(
            page.get_by_role("group", name="GenBank File selection", exact=True)
        ).to_contain_text("NC_001416.gb")
        prefix = page.get_by_label("Output Prefix", exact=True)
        prefix.fill("lambda_genbank")
        genbank_report = generate_and_inspect(
            page,
            lambda result: inspect_lambda_input_svg(result, page),
            assert_lambda_input_svg,
        )
        fit_complete_linear_preview(page)
        genbank.scroll_into_view_if_needed()
        screenshot_bytes["genbank-input.png"] = capture_screenshot(
            page, output_paths["genbank-input.png"], "Linear"
        )

        gff3_fasta = page.get_by_role(
            "radio", name="GFF3 + FASTA", exact=True
        )
        gff3_fasta.check()
        expect(gff3_fasta).to_be_checked()
        page.get_by_label("GFF3", exact=True).set_input_files(GUI_INPUTS_GFF3_PATH)
        page.get_by_label("FASTA", exact=True).set_input_files(GUI_INPUTS_FASTA_PATH)
        expect(page.get_by_role("group", name="GFF3 selection", exact=True)).to_contain_text(
            "NC_001416.gff3"
        )
        expect(page.get_by_role("group", name="FASTA selection", exact=True)).to_contain_text(
            "NC_001416.fna"
        )
        prefix.fill("lambda_gff3")
        gff3_report = generate_and_inspect(
            page,
            lambda result: inspect_lambda_input_svg(result, page),
            assert_lambda_input_svg,
        )
        assert_matching_lambda_input_semantics(genbank_report, gff3_report)
        fit_complete_linear_preview(page)
        gff3_fasta.scroll_into_view_if_needed()
        screenshot_bytes["gff3-fasta-input.png"] = capture_screenshot(
            page, output_paths["gff3-fasta-input.png"], "Linear"
        )

        svg_button = page.get_by_role("button", name="SVG", exact=True)
        with page.expect_download(timeout=ACTION_TIMEOUT_MS) as download_info:
            svg_button.click()
        download = download_info.value
        if download.failure() is not None:
            raise AssertionError(f"SVG download failed: {download.failure()}")
        if download.suggested_filename != EXPECTED_GUI_INPUTS_SVG:
            raise AssertionError(
                "Unexpected SVG download name: "
                f"expected {EXPECTED_GUI_INPUTS_SVG}, "
                f"found {download.suggested_filename}"
            )
        downloaded_svg = download_dir / download.suggested_filename
        download.save_as(downloaded_svg)
        download_report = assert_gui_inputs_svg_download(downloaded_svg)

        page.get_by_label("FASTA", exact=True).set_input_files(mismatched_fasta_path)
        expect(page.get_by_role("group", name="FASTA selection", exact=True)).to_contain_text(
            "mismatched_sequence_id.fna"
        )
        generate = page.get_by_role("button", name="Generate Diagram", exact=True)
        generate_and_wait_for_result(page, expected_status="error")
        error_alert = page.get_by_role("alert", name="Generation Error", exact=True)
        expect(error_alert).to_be_visible(timeout=GENERATION_TIMEOUT_MS)
        expect(generate).to_be_enabled(timeout=GENERATION_TIMEOUT_MS)
        expect(error_alert).to_contain_text(EXPECTED_ID_MISMATCH)
        gff3_fasta.scroll_into_view_if_needed()
        screenshot_bytes["id-error.png"] = capture_screenshot(
            page, output_paths["id-error.png"], "Linear"
        )

        capture.assert_clean()
    finally:
        capture.close()

    return CaptureResult(
        screenshot_bytes=screenshot_bytes,
        first_svg_semantics=genbank_report,
        final_svg_semantics=gff3_report,
        download=download_report,
        validation_message=EXPECTED_ID_MISMATCH,
    )
