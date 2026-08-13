"""Capture the Interactive SVG Gallery chloroplast-map tutorial."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Mapping

from playwright.sync_api import BrowserType, Page, expect

from flows.how_to.tracks import (
    GUI_ANNOTATION_GENBANK_PATH,
    GUI_ANNOTATION_GENBANK_SHA256,
    GUI_ANNOTATION_GENBANK_SIZE,
    GUI_ANNOTATION_TABLE_PATH,
    _assert_safe_svg,
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
    assert_fixture_identity,
    assert_output_paths,
    capture_screenshot,
    generate_and_inspect,
    open_browser_capture,
    set_feature_search_visible,
    wait_for_app_shell,
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
SPECIFIC_COLORS_PATH = GUI_ANNOTATION_GENBANK_PATH.with_name(
    "chloroplast_specific_table.tsv"
)
SPECIFIC_COLORS_SIZE = 3_033
SPECIFIC_COLORS_SHA256 = (
    "406f1a042ec072c0026c014fc173950793582ffe4c8c568ec2539e2db98df0ec"
)
QUALIFIER_PRIORITY_PATH = GUI_ANNOTATION_GENBANK_PATH.with_name(
    "qualifier_priority.tsv"
)
GALLERY_FEATURE_TYPES = (
    "CDS",
    "rRNA",
    "tRNA",
    "tmRNA",
    "ncRNA",
    "misc_RNA",
    "rep_origin",
)
GALLERY_SLOT_ORDER = (
    ("features", "features"),
    ("plastome_regions", "annotations"),
    ("gc_content", "dinucleotide_content"),
)


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


def _assert_gallery_chloroplast(report: Mapping[str, Any]) -> None:
    _assert_safe_svg(report)
    if "NC_001879.2" not in report["recordIds"]:
        raise AssertionError("The Gallery-style map does not identify NC_001879.2")
    if report["featureElementCount"] != 197:
        raise AssertionError(
            "Expected 197 rendered Gallery feature elements, found "
            f"{report['featureElementCount']}"
        )
    if report.get("logicalFeatureCount") not in (None, 147):
        raise AssertionError(
            "Expected 147 logical Gallery features, found "
            f"{report['logicalFeatureCount']}"
        )
    slots = tuple(
        (slot["slotId"], slot["renderer"])
        for slot in report["slots"]
        if not str(slot["slotId"]).startswith("__gbdraw_auto_")
    )
    if slots != GALLERY_SLOT_ORDER:
        raise AssertionError(f"Unexpected Gallery chloroplast slots: {slots!r}")
    annotations = {
        (
            item["id"],
            item["setId"],
            item["trackId"],
            item["recordId"],
            item["mark"],
            item["label"],
        )
        for item in report["annotations"]
    }
    expected_annotations = {
        (key, "plastome_regions", "plastome_regions", "NC_001879.2", "bracket", label)
        for key, label in (("lsc", "LSC"), ("irb", "IRb"), ("ssc", "SSC"), ("ira", "IRa"))
    }
    if annotations != expected_annotations:
        raise AssertionError(f"Unexpected plastome region annotations: {annotations!r}")
    texts = set(report["texts"])
    required = {
        "Nicotiana tabacum",
        "NC_001879.2",
        "155,943 bp",
        "LSC",
        "IRb",
        "SSC",
        "IRa",
        "GC content",
        "matK",
        "psaA",
        "photosystem I",
        "photosystem II",
        "RNA polymerase",
        "rep_origin",
    }
    missing = required - texts
    if missing:
        raise AssertionError(f"Gallery chloroplast text is missing: {sorted(missing)}")
    if {"GC skew (+)", "GC skew (-)", "AT skew (+)", "AT skew (-)"} & texts:
        raise AssertionError("The Gallery chloroplast map must not contain a skew track")


def _feature_details(page: Page) -> Any:
    return page.get_by_label("Features", exact=True).locator("xpath=..")


def _configure_feature_types(page: Page) -> None:
    details = _feature_details(page)
    page.get_by_label("Features", exact=True).click()
    current = tuple(
        page.evaluate(
            "() => (window.__GBDRAW_APP__?.adv?.features || []).map(String)"
        )
    )
    for feature_type in current:
        if feature_type in GALLERY_FEATURE_TYPES:
            continue
        feature_name = details.get_by_text(feature_type, exact=True).first
        feature_name.locator("xpath=..").get_by_role("button").click()
    picker = details.locator('select:has(option[value="rep_origin"])')
    expect(picker).to_have_count(1)
    add_button = details.get_by_role("button", name=re.compile(r"Add$"))
    for feature_type in GALLERY_FEATURE_TYPES:
        selected = tuple(
            page.evaluate(
                "() => (window.__GBDRAW_APP__?.adv?.features || []).map(String)"
            )
        )
        if feature_type in selected:
            continue
        picker.select_option(feature_type)
        add_button.click()
    selected = tuple(
        page.evaluate(
            "() => (window.__GBDRAW_APP__?.adv?.features || []).map(String)"
        )
    )
    if selected != GALLERY_FEATURE_TYPES:
        raise AssertionError(f"Unexpected chloroplast feature types: {selected!r}")
    page.get_by_label("Block Stroke Width", exact=True).fill("1")
    page.get_by_label("Line Stroke Width", exact=True).fill("2")
    page.get_by_label("Features", exact=True).click()


def _configure_gallery_presentation(page: Page) -> None:
    page.get_by_label("Output Prefix", exact=True).fill("annotated_chloroplast_map")
    page.get_by_label("Species", exact=True).fill("<i>Nicotiana tabacum</i>")
    page.get_by_label("Track Preset", exact=True).select_option("tuckin")
    page.get_by_label("Separate Strands", exact=True).check()
    page.get_by_label("Hide GC Content", exact=True).uncheck()
    page.get_by_label("Hide GC Skew", exact=True).check()

    _configure_feature_types(page)

    colors = page.get_by_label("Colors", exact=True)
    colors.click()
    page.get_by_label("Specific Table (-t)", exact=True).set_input_files(
        SPECIFIC_COLORS_PATH
    )
    colors.click()

    labels = page.get_by_label("Labels", exact=True)
    labels.click()
    page.get_by_label("Label Mode", exact=True).select_option("both")
    for label, value in (
        ("Outer X Offset", "0.9"),
        ("Outer Y Offset", "0.9"),
        ("Inner X Offset", "0.975"),
        ("Inner Y Offset", "0.975"),
    ):
        page.get_by_text(label, exact=True).locator("xpath=..").locator(
            'input[type="number"]'
        ).fill(value)
    page.get_by_text("Circular Label Placement", exact=True).locator(
        "xpath=.."
    ).locator("select").select_option("radial")
    page.get_by_label("Priority File (TSV)", exact=True).set_input_files(
        QUALIFIER_PRIORITY_PATH
    )
    labels.click()

    axis = page.locator("summary").filter(has_text="Axis & Scale")
    axis.click()
    page.get_by_label("Axis Stroke Width", exact=True).fill("3")
    axis.click()

    title = page.get_by_label("Title & Legend", exact=True)
    title.click()
    page.get_by_label("Plot Title", exact=True).fill("")
    page.get_by_label("Plot Title Position", exact=True).select_option("none")
    page.get_by_label("Definition Font Size", exact=True).fill("28")
    page.get_by_label("Legend Position", exact=True).select_option("upper_left")
    title.click()


def _remove_slot_if_present(page: Page, slot_id: str) -> None:
    group = page.get_by_role(
        "group", name=f"Circular track slot {slot_id}", exact=True
    )
    if group.count():
        group.get_by_title("Remove", exact=True).click()


def _configure_gallery_slots(page: Page) -> tuple[dict[str, Any], ...]:
    custom_slots = page.get_by_role(
        "button", name=re.compile(r"Custom Track Slots$"), exact=False
    )
    custom_slots.click()
    page.get_by_role("checkbox", name="Use custom stack", exact=True).check()
    _remove_slot_if_present(page, "ticks")
    _remove_slot_if_present(page, "gc_skew")

    feature_slot = page.get_by_role(
        "group", name="Circular track slot features", exact=True
    )
    feature_slot.locator("select").last.select_option("split")

    renderer = page.get_by_label("New circular track renderer", exact=True)
    renderer.select_option("annotations")
    page.get_by_role("button", name=re.compile(r"Add track$"), exact=False).click()
    annotation_slot = page.get_by_role(
        "group", name="Circular track slot annotations", exact=True
    )
    annotation_slot.get_by_label("Annotation set", exact=True).select_option(
        "plastome_regions"
    )
    annotation_slot.get_by_label("Annotation placement", exact=True).select_option(
        "inside"
    )
    annotation_slot = page.get_by_role(
        "group", name="Circular track slot annotations", exact=True
    )
    annotation_slot.get_by_label(
        "Circular track slot id annotations", exact=True
    ).fill("plastome_regions")
    annotation_slot = page.get_by_role(
        "group", name="Circular track slot plastome_regions", exact=True
    )
    annotation_slot.get_by_title("Move outward", exact=True).click()
    annotation_slot = page.get_by_role(
        "group", name="Circular track slot plastome_regions", exact=True
    )
    annotation_slot.get_by_title("Width", exact=True).fill("20px")
    annotation_slot.get_by_title("Radius", exact=True).fill("0.65")
    annotation_slot.get_by_title("Inner gap", exact=True).fill("1")
    annotation_slot.get_by_title("Outer gap", exact=True).fill("1")
    annotation_slot.get_by_label("Show annotation labels", exact=True).check()
    annotation_slot.locator(
        'select:has(option[value="compress"])'
    ).select_option("compress")
    annotation_numbers = annotation_slot.locator('input[type="number"]')
    expect(annotation_numbers).to_have_count(2)
    annotation_numbers.last.fill("1")

    gc_slot = page.get_by_role(
        "group", name="Circular track slot gc_content", exact=True
    )
    gc_slot.get_by_title("Width", exact=True).fill("0.08")
    gc_slot.get_by_title("Radius", exact=True).fill("0.56")

    slots = _track_slot_snapshot(page)
    enabled = tuple(
        (slot["id"], slot["renderer"])
        for slot in slots
        if slot["enabled"]
    )
    if enabled != GALLERY_SLOT_ORDER:
        raise AssertionError(f"Unexpected Gallery slot state: {enabled!r}")
    return slots


def capture_gui_annotated_chloroplast(
    browser_type: BrowserType,
    base_url: str,
    output_paths: Mapping[str, Path],
    download_dir: Path,
) -> CaptureResult:
    """Build and export the Gallery-quality T-GUI-05 tobacco plastome."""

    genome_report = _validate_complete_record(
        GUI_ANNOTATION_GENBANK_PATH,
        expected_size=GUI_ANNOTATION_GENBANK_SIZE,
        expected_sha256=GUI_ANNOTATION_GENBANK_SHA256,
        expected_id="NC_001879.2",
        expected_length=155_943,
    )
    annotation_report = _validate_annotation_fixture()
    assert_fixture_identity(
        SPECIFIC_COLORS_PATH,
        expected_size=SPECIFIC_COLORS_SIZE,
        expected_sha256=SPECIFIC_COLORS_SHA256,
    )
    assert_output_paths(output_paths, SCREENSHOT_NAMES, SCENARIO_ID)
    download_dir.mkdir(parents=True, exist_ok=True)

    capture = open_browser_capture(browser_type, base_url)
    page = capture.page
    screenshot_bytes: dict[str, int] = {}
    try:
        page.goto(base_url, wait_until="domcontentloaded")
        wait_for_app_shell(page)
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
        _fit_circular_preview(
            page,
            target_zoom="70%",
            pan_left_ratio=0.0,
        )
        set_feature_search_visible(page, visible=False)
        screenshot_bytes[SCREENSHOT_NAMES[1]] = capture_screenshot(
            page, output_paths[SCREENSHOT_NAMES[1]], "Circular"
        )

        _configure_gallery_presentation(page)
        _fit_circular_preview(
            page,
            target_zoom="70%",
            pan_left_ratio=0.0,
        )
        annotations = page.get_by_label("Region Annotations", exact=True)
        annotations.click()
        page.get_by_label("Import TSV", exact=True).set_input_files(
            GUI_ANNOTATION_TABLE_PATH
        )
        expect(page.get_by_label("Annotation set id", exact=True)).to_have_value(
            "plastome_regions"
        )
        expect(page.get_by_label("Annotation lane", exact=True)).to_have_count(4)
        annotations.scroll_into_view_if_needed()
        screenshot_bytes[SCREENSHOT_NAMES[2]] = capture_screenshot(
            page, output_paths[SCREENSHOT_NAMES[2]], "Circular"
        )

        slots = _configure_gallery_slots(page)
        annotation_slot = page.get_by_role(
            "group", name="Circular track slot plastome_regions", exact=True
        )
        annotation_slot.scroll_into_view_if_needed()
        screenshot_bytes[SCREENSHOT_NAMES[3]] = capture_screenshot(
            page, output_paths[SCREENSHOT_NAMES[3]], "Circular"
        )

        final_report = generate_and_inspect(
            page, _inspect_tracks_svg, _assert_gallery_chloroplast
        )
        _fit_circular_preview(
            page, target_zoom="50%", pan_left_ratio=0.32
        )
        popup = page.get_by_role("dialog", name=re.compile(r"^Feature details:"))
        if popup.is_visible():
            popup.get_by_role(
                "button", name="Close feature popup", exact=True
            ).click()
        screenshot_bytes[SCREENSHOT_NAMES[4]] = capture_screenshot(
            page, output_paths[SCREENSHOT_NAMES[4]], "Circular"
        )
        download_report = _download_svg(
            page,
            download_dir,
            expected_filename=OUTPUT_FILENAME,
            assert_svg=_assert_gallery_chloroplast,
        )
        download_report["galleryFigureSemanticsParity"] = True
        capture.assert_clean()
    finally:
        capture.close()

    return CaptureResult(
        screenshot_bytes=screenshot_bytes,
        final_svg_semantics=final_report,
        download=download_report,
        fixture_report={**genome_report, "annotations": annotation_report},
        track_slots=slots,
    )


__all__ = ["capture_gui_annotated_chloroplast"]
