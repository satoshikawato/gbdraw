"""Capture the annotated-chloroplast project tutorial through visible controls."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Mapping

from playwright.sync_api import BrowserType, expect

from flows.how_to.tracks import (
    EXPECTED_ANNOTATION_SLOT_ORDER,
    GUI_ANNOTATION_GENBANK_PATH,
    GUI_ANNOTATION_GENBANK_SHA256,
    GUI_ANNOTATION_GENBANK_SIZE,
    GUI_ANNOTATION_TABLE_PATH,
    _assert_annotation_svg,
    _assert_safe_svg,
    _assert_slot_snapshot,
    _configure_title,
    _download_svg,
    _fit_circular_preview,
    _inspect_tracks_svg,
    _resize_sidebar,
    _track_slot_snapshot,
    _validate_annotation_fixture,
    _validate_complete_record,
    CaptureResult,
)
from flows.web_capture import (
    assert_output_paths,
    capture_screenshot,
    generate_and_inspect,
    open_browser_capture,
    set_feature_search_visible,
    wait_for_worker,
)


SCENARIO_ID = "T-GUI-05"
SCREENSHOT_NAMES = (
    "01-input-ready.png",
    "02-first-diagram.png",
    "03-annotation-table.png",
    "04-track-settings.png",
    "05-finished-diagram.png",
)
OUTPUT_FILENAME = "annotated_chloroplast_map.svg"


def _assert_plain_plastome(report: Mapping[str, Any]) -> None:
    _assert_safe_svg(report)
    if "NC_001879.2" not in report["recordIds"]:
        raise AssertionError("The first diagram does not identify NC_001879.2")
    if report["featureElementCount"] != 195:
        raise AssertionError("The first diagram does not contain the complete plastome")
    if report.get("logicalFeatureCount") not in (None, 145):
        raise AssertionError("The first diagram has an unexpected logical feature count")
    if report["annotations"]:
        raise AssertionError("The Step 2 diagram already contains region annotations")


def capture_gui_annotated_chloroplast(
    browser_type: BrowserType,
    base_url: str,
    output_paths: Mapping[str, Path],
    download_dir: Path,
) -> CaptureResult:
    """Build and export the complete T-GUI-05 tobacco-plastome figure."""

    genome_report = _validate_complete_record(
        GUI_ANNOTATION_GENBANK_PATH,
        expected_size=GUI_ANNOTATION_GENBANK_SIZE,
        expected_sha256=GUI_ANNOTATION_GENBANK_SHA256,
        expected_id="NC_001879.2",
        expected_length=155_943,
    )
    annotation_report = _validate_annotation_fixture()
    assert_output_paths(output_paths, SCREENSHOT_NAMES, SCENARIO_ID)
    download_dir.mkdir(parents=True, exist_ok=True)

    capture = open_browser_capture(browser_type, base_url)
    page = capture.page
    screenshot_bytes: dict[str, int] = {}
    try:
        page.goto(base_url, wait_until="domcontentloaded")
        wait_for_worker(page)
        _resize_sidebar(page)

        circular = page.get_by_role("button", name="Circular", exact=True)
        circular.click()
        expect(circular).to_have_attribute("aria-pressed", "true")
        page.get_by_role("radio", name="GenBank", exact=True).check()
        page.get_by_label("GenBank/DDBJ File", exact=True).set_input_files(
            GUI_ANNOTATION_GENBANK_PATH
        )
        expect(
            page.get_by_role(
                "group", name="GenBank/DDBJ File selection", exact=True
            )
        ).to_contain_text("NC_001879.gbk")
        page.get_by_label("Output Prefix", exact=True).fill(
            "annotated_chloroplast_map"
        )
        screenshot_bytes[SCREENSHOT_NAMES[0]] = capture_screenshot(
            page, output_paths[SCREENSHOT_NAMES[0]], "Circular"
        )

        generate_and_inspect(page, _inspect_tracks_svg, _assert_plain_plastome)
        _fit_circular_preview(page, target_zoom="60%")
        set_feature_search_visible(page, visible=False)
        screenshot_bytes[SCREENSHOT_NAMES[1]] = capture_screenshot(
            page, output_paths[SCREENSHOT_NAMES[1]], "Circular"
        )

        annotations_section = page.get_by_label("Region Annotations", exact=True)
        annotations_section.click()
        page.get_by_label("Import TSV", exact=True).set_input_files(
            GUI_ANNOTATION_TABLE_PATH
        )
        expect(page.get_by_label("Annotation set id", exact=True)).to_have_value(
            "plastome_regions"
        )
        legend_label = page.get_by_placeholder("Set legend label (optional)")
        legend_label.fill("Plastome structural regions")
        lane_controls = page.get_by_label("Annotation lane", exact=True)
        expect(lane_controls).to_have_count(4)
        for index, lane in enumerate((0, 1, 0, 1)):
            lane_controls.nth(index).fill(str(lane))
        rendered_lanes = tuple(
            page.evaluate(
                """
                () => (window.__GBDRAW_APP__?.annotationSets?.[0]?.annotations || [])
                  .map((item) => Number(item?.lane))
                """
            )
        )
        if rendered_lanes != (0, 1, 0, 1):
            raise AssertionError(
                f"Unexpected visible annotation-lane assignment: {rendered_lanes}"
            )
        annotations_section.scroll_into_view_if_needed()
        screenshot_bytes[SCREENSHOT_NAMES[2]] = capture_screenshot(
            page, output_paths[SCREENSHOT_NAMES[2]], "Circular"
        )

        custom_slots = page.get_by_role(
            "button", name=re.compile(r"Custom Track Slots$"), exact=False
        )
        custom_slots.click()
        use_custom = page.get_by_role(
            "checkbox", name="Use custom stack", exact=True
        )
        use_custom.check()

        skew_slot = page.get_by_role(
            "group", name="Circular track slot gc_skew", exact=True
        )
        skew_slot.get_by_label("Track dinucleotide", exact=True).fill("AT")
        skew_slot.get_by_label("Track legend label", exact=True).fill("AT skew")

        new_renderer = page.get_by_label("New circular track renderer", exact=True)
        new_renderer.select_option("annotations")
        page.get_by_role("button", name=re.compile(r"Add track$"), exact=False).click()
        annotation_slot = page.get_by_role(
            "group", name="Circular track slot annotations", exact=True
        )
        annotation_slot.get_by_label("Annotation set", exact=True).select_option(
            "plastome_regions"
        )
        annotation_slot.get_by_title("Move outside Axis", exact=True).click()
        annotation_slot = page.get_by_role(
            "group", name="Circular track slot annotations", exact=True
        )
        annotation_slot.get_by_label(
            "Show annotation labels", exact=True
        ).check()

        slots = _track_slot_snapshot(page)
        _assert_slot_snapshot(slots, EXPECTED_ANNOTATION_SLOT_ORDER)
        skew_snapshot = next(slot for slot in slots if slot["id"] == "gc_skew")
        annotation_snapshot = next(
            slot for slot in slots if slot["id"] == "annotations"
        )
        if skew_snapshot["params"].get("nt") != "AT":
            raise AssertionError("The custom skew track does not use AT")
        if annotation_snapshot["side"] != "outside":
            raise AssertionError("The annotation track is not outside the axis")

        labels = page.get_by_label("Labels", exact=True)
        labels.click()
        page.get_by_label("Label Mode", exact=True).select_option("none")
        labels.click()
        _configure_title(page, "Complete Nicotiana tabacum plastome regions")

        annotation_slot.scroll_into_view_if_needed()
        screenshot_bytes[SCREENSHOT_NAMES[3]] = capture_screenshot(
            page, output_paths[SCREENSHOT_NAMES[3]], "Circular"
        )

        final_report = generate_and_inspect(
            page, _inspect_tracks_svg, _assert_annotation_svg
        )
        _fit_circular_preview(page, target_zoom="60%")
        feature_popup = page.get_by_role(
            "dialog", name=re.compile(r"^Feature details:")
        )
        if feature_popup.is_visible():
            feature_popup.get_by_role(
                "button", name="Close feature popup", exact=True
            ).click()
        screenshot_bytes[SCREENSHOT_NAMES[4]] = capture_screenshot(
            page, output_paths[SCREENSHOT_NAMES[4]], "Circular"
        )
        download_report = _download_svg(
            page,
            download_dir,
            expected_filename=OUTPUT_FILENAME,
            assert_svg=_assert_annotation_svg,
        )
        capture.assert_clean()
    finally:
        capture.close()

    return CaptureResult(
        screenshot_bytes=screenshot_bytes,
        final_svg_semantics=final_report,
        download=download_report,
        fixture_report={
            **genome_report,
            "annotations": {**annotation_report, "renderedLanes": rendered_lanes},
        },
        track_slots=slots,
    )


__all__ = ["capture_gui_annotated_chloroplast"]
