"""Capture the GUI variant of the mitochondrial feature-presentation project."""

from __future__ import annotations

import hashlib
import re
from pathlib import Path
from typing import Any, Mapping

from Bio import SeqIO
from playwright.sync_api import BrowserType, Page, expect

from assertions.svg_semantics import assert_static_svg_safety
from config import ACTION_TIMEOUT_MS
from flows.how_to.presentation import (
    GENE_PRIORITY_RULE_PATH,
    GENE_PRIORITY_RULE_SHA256,
    GENE_PRIORITY_RULE_SIZE,
    HUMAN_MITOCHONDRION_PATH,
    HUMAN_MITOCHONDRION_SHA256,
    HUMAN_MITOCHONDRION_SIZE,
    _fit_finished_circular_preview,
    _inspect_downloaded_presentation_svg,
    _inspect_presentation_svg,
    _open_human_circular,
    _pan_preview_left,
    _park_feature_search,
    _wait_for_preview_transform,
    CaptureResult,
)
from flows.how_to.tracks import (
    _inspect_tracks_svg,
    _inspect_tracks_svg_file,
    _resize_sidebar,
    _track_slot_snapshot,
)
from flows.web_capture import (
    assert_fixture_identity,
    assert_output_paths,
    capture_screenshot,
    open_browser_capture,
    wait_for_worker,
)


SCENARIO_ID = "T-GUI-10"
SCREENSHOT_NAMES = ("presentation-settings.png", "presentation-result.png")
OUTPUT_PREFIX = "mitochondrial_features_highlighted"
OUTPUT_FILENAME = f"{OUTPUT_PREFIX}.svg"
TITLE = "Human mitochondrial feature presentation"
EXPECTED_SLOTS = (
    ("ticks", "ticks"),
    ("features", "features"),
    ("mitochondrial_regions", "annotations"),
)

COLOR_TABLE = """CDS\tgene\t^ND(4L|[1-6])$\t#3B82F6\tNADH dehydrogenase
CDS\tgene\t^COX[1-3]$\t#EF4444\tCytochrome c oxidase
CDS\tgene\t^ATP[68]$\t#F59E0B\tATP synthase
CDS\tgene\t^CYTB$\t#8B5CF6\tCytochrome b
rRNA\tgene\t^RNR[12]$\t#10B981\tRibosomal RNA
"""
LABEL_WHITELIST = """CDS\tgene\t^(ND1|COX2|ATP8|ATP6|COX3|CYTB)$
rRNA\tgene\t^RNR[12]$
"""
LABEL_OVERRIDES = """record_id\tfeature_type\tqualifier\tvalue\tlabel_text
NC_012920.1\tCDS\tlabel\t^ND1$\tComplex I (ND1)
NC_012920.1\tCDS\tlabel\t^COX2$\tOxidase II
NC_012920.1\trRNA\tlabel\t^s-rRNA$\t12S rRNA
NC_012920.1\trRNA\tlabel\t^l-rRNA$\t16S rRNA
"""
REGION_TABLE = """set_id\tid\tmark\trecord\tstart\tend\tcoordinate_space\twraps_origin\tlabel\tlane\tstroke\tstroke_width\tline_cap\tlabel_color\tlabel_font_size\tlabel_orientation\tlabel_offset
mitochondrial_regions\td_loop\tbracket\tNC_012920.1\t16024\t576\tsource\ttrue\tD-loop\t0\t#202020\t3\ttick\t#202020\t14\ttangent\t7
"""


def _write_tables(download_dir: Path) -> dict[str, Path]:
    tables = {
        "colors": download_dir / "presentation_colors.tsv",
        "whitelist": download_dir / "presentation_labels.tsv",
        "overrides": download_dir / "presentation_label_overrides.tsv",
        "regions": download_dir / "mitochondrial_regions.tsv",
    }
    contents = {
        "colors": COLOR_TABLE,
        "whitelist": LABEL_WHITELIST,
        "overrides": LABEL_OVERRIDES,
        "regions": REGION_TABLE,
    }
    for key, path in tables.items():
        path.write_text(contents[key], encoding="utf-8")
        if path.read_text(encoding="utf-8") != contents[key]:
            raise AssertionError(f"Could not freeze {path.name}")
    return tables


def _validate_source_inputs() -> dict[str, Any]:
    assert_fixture_identity(
        HUMAN_MITOCHONDRION_PATH,
        expected_size=HUMAN_MITOCHONDRION_SIZE,
        expected_sha256=HUMAN_MITOCHONDRION_SHA256,
    )
    assert_fixture_identity(
        GENE_PRIORITY_RULE_PATH,
        expected_size=GENE_PRIORITY_RULE_SIZE,
        expected_sha256=GENE_PRIORITY_RULE_SHA256,
    )
    record = SeqIO.read(HUMAN_MITOCHONDRION_PATH, "genbank")
    if (
        record.id != "NC_012920.1"
        or len(record) != 16_569
        or str(record.annotations.get("topology", "")).lower() != "circular"
        or "complete" not in record.description.lower()
    ):
        raise AssertionError("T-GUI-10 requires complete circular NC_012920.1")
    return {
        "record_id": record.id,
        "length": len(record),
        "topology": record.annotations.get("topology"),
        "description": record.description,
    }


def _remove_slot_if_present(page: Page, slot_id: str) -> None:
    slot = page.get_by_role(
        "group", name=f"Circular track slot {slot_id}", exact=True
    )
    if slot.count():
        slot.get_by_title("Remove", exact=True).click()


def _configure_slots(page: Page) -> tuple[dict[str, Any], ...]:
    page.get_by_role(
        "button", name=re.compile(r"Custom Track Slots$"), exact=False
    ).click()
    page.get_by_role("checkbox", name="Use custom stack", exact=True).check()
    _remove_slot_if_present(page, "gc_content")
    _remove_slot_if_present(page, "gc_skew")

    ticks = page.get_by_role(
        "group", name="Circular track slot ticks", exact=True
    )
    move_outside = ticks.get_by_title("Move outside Axis", exact=True)
    if move_outside.is_enabled():
        move_outside.click()
    ticks = page.get_by_role(
        "group", name="Circular track slot ticks", exact=True
    )
    ticks.locator("select").last.select_option("label_out_tick_in")

    features = page.get_by_role(
        "group", name="Circular track slot features", exact=True
    )
    features.locator("select").last.select_option("split")

    page.get_by_label("New circular track renderer", exact=True).select_option(
        "annotations"
    )
    page.get_by_role("button", name=re.compile(r"Add track$"), exact=False).click()
    annotation = page.get_by_role(
        "group", name="Circular track slot annotations", exact=True
    )
    annotation.get_by_label("Annotation set", exact=True).select_option(
        "mitochondrial_regions"
    )
    annotation.get_by_label("Annotation placement", exact=True).select_option(
        "inside"
    )
    annotation.get_by_label(
        "Circular track slot id annotations", exact=True
    ).fill("mitochondrial_regions")
    annotation = page.get_by_role(
        "group", name="Circular track slot mitochondrial_regions", exact=True
    )
    annotation.get_by_title("Width", exact=True).fill("24px")
    annotation.get_by_label("Show annotation labels", exact=True).check()
    annotation.locator('select:has(option[value="compress"])').select_option(
        "compress"
    )
    annotation.locator('input[type="number"]').last.fill("1")

    slots = _track_slot_snapshot(page)
    actual = tuple(
        (slot["id"], slot["renderer"])
        for slot in slots
        if slot["enabled"]
    )
    if actual != EXPECTED_SLOTS:
        raise AssertionError(f"Unexpected feature-presentation slots: {actual!r}")
    return slots


def _merge_browser_report(result_region: Any) -> dict[str, Any]:
    report = _inspect_presentation_svg(result_region)
    tracks = _inspect_tracks_svg(result_region)
    report["slots"] = tracks["slots"]
    report["annotations"] = tracks["annotations"]
    return report


def _slot_pairs(report: Mapping[str, Any]) -> tuple[tuple[str, str], ...]:
    return tuple(
        (slot["slotId"], slot["renderer"])
        for slot in report.get("slots", [])
        if not str(slot["slotId"]).startswith("__gbdraw_auto_")
    )


def _assert_highlight_svg(report: Mapping[str, Any]) -> None:
    assert_static_svg_safety(dict(report))
    if "NC_012920.1" not in set(report.get("recordIds", [])):
        raise AssertionError("The map does not contain complete NC_012920.1")
    if report.get("featureElementCount") != 37:
        raise AssertionError(
            f"Expected all 37 mitochondrial features, found {report.get('featureElementCount')}"
        )
    texts = set(report.get("texts", []))
    required = {
        "NC_012920.1",
        "16,569 bp",
        TITLE,
        "Complex I (ND1)",
        "Oxidase II",
        "12S rRNA",
        "16S rRNA",
        "D-loop",
    }
    missing = required - texts
    if missing:
        raise AssertionError(f"Missing feature-presentation text: {sorted(missing)}")
    fills = {feature.get("fill") for feature in report.get("features", [])}
    expected_fills = {"#3b82f6", "#ef4444", "#f59e0b", "#8b5cf6", "#10b981"}
    if not expected_fills <= fills:
        raise AssertionError(
            f"Missing feature-presentation fills: {sorted(expected_fills - fills)}"
        )
    if _slot_pairs(report) != EXPECTED_SLOTS:
        raise AssertionError(f"Unexpected rendered slots: {_slot_pairs(report)!r}")
    annotations = {
        (
            row.get("id"),
            row.get("setId"),
            row.get("trackId"),
            row.get("recordId"),
            row.get("mark"),
            row.get("label"),
        )
        for row in report.get("annotations", [])
    }
    expected_annotation = {
        (
            "d_loop",
            "mitochondrial_regions",
            "mitochondrial_regions",
            "NC_012920.1",
            "bracket",
            "D-loop",
        )
    }
    if annotations != expected_annotation:
        raise AssertionError(f"Unexpected D-loop annotation: {annotations!r}")


def _assert_control_state(report: Mapping[str, Any]) -> None:
    state = report.get("state", {})
    expected = {
        "prefix": OUTPUT_PREFIX,
        "title": TITLE,
        "titlePosition": "top",
        "legendPosition": "right",
        "labelsMode": "both",
        "separateStrands": True,
        "trackPreset": "spreadout",
        "palette": "default",
        "filterMode": "Whitelist",
        "priorityFile": GENE_PRIORITY_RULE_PATH.name,
        "arrowHeadLengthRatio": "",
        "arrowShaftWidthRatio": 0.72,
        "blockStrokeColor": "#1f2937",
        "blockStrokeWidth": 1.5,
        "lineStrokeColor": "#9ca3af",
        "lineStrokeWidth": 1.5,
        "resolveOverlaps": False,
    }
    for key, value in expected.items():
        if state.get(key) != value:
            raise AssertionError(
                f"Unexpected {SCENARIO_ID} control state for {key}: {state.get(key)!r}"
            )
    shapes = state.get("featureShapes", {})
    if {key: shapes.get(key) for key in ("CDS", "rRNA", "tRNA")} != {
        "CDS": "arrow",
        "rRNA": "rectangle",
        "tRNA": "arrow",
    }:
        raise AssertionError(f"Unexpected feature shapes: {shapes!r}")


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
    report = _inspect_downloaded_presentation_svg(path)
    tracks = _inspect_tracks_svg_file(path)
    report["slots"] = tracks["slots"]
    report["annotations"] = tracks["annotations"]
    _assert_highlight_svg(report)
    return {
        "filename": path.name,
        "bytes": path.stat().st_size,
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "path": str(path),
        "semantics": report,
    }


def _load_label_overrides(page: Page, path: Path) -> None:
    page.get_by_title("Open editor", exact=True).click()
    drawer = page.locator(".right-drawer")
    expect(drawer).to_have_attribute("aria-hidden", "false")
    drawer.locator("button").filter(has_text="Features").click()
    page.once("dialog", lambda dialog: dialog.accept())
    with page.expect_file_chooser() as chooser:
        drawer.get_by_role("button", name="Load Label TSV", exact=True).click()
    chooser.value.set_files(path)
    page.wait_for_timeout(1_000)
    processing = page.get_by_text("Reflowing label placement...", exact=True)
    if processing.count():
        expect(processing).to_be_hidden(timeout=ACTION_TIMEOUT_MS)
    close = page.get_by_title("Close editor", exact=True)
    if close.is_visible():
        close.click()
    expect(page.get_by_title("Open editor", exact=True)).to_be_visible()


def capture_gui_feature_highlight(
    browser_type: BrowserType,
    base_url: str,
    output_paths: Mapping[str, Path],
    download_dir: Path,
) -> CaptureResult:
    """Build the CLI/Python-equivalent feature presentation in the web app."""

    source_record = _validate_source_inputs()
    assert_output_paths(output_paths, SCREENSHOT_NAMES, SCENARIO_ID)
    download_dir.mkdir(parents=True, exist_ok=True)
    tables = _write_tables(download_dir)

    capture = open_browser_capture(browser_type, base_url)
    page = capture.page
    screenshot_bytes: dict[str, int] = {}
    try:
        page.goto(base_url, wait_until="domcontentloaded")
        wait_for_worker(page)
        _resize_sidebar(page)
        _open_human_circular(page, OUTPUT_PREFIX)

        page.get_by_label("Species", exact=True).fill("<i>Homo sapiens</i>")
        page.get_by_label("Track Preset", exact=True).select_option("spreadout")
        page.get_by_label("Separate Strands", exact=True).check()

        colors = page.get_by_label("Colors", exact=True)
        colors.click()
        expect(page.get_by_label("Palette", exact=True)).to_have_value("default")
        page.get_by_label("Specific Table (-t)", exact=True).set_input_files(
            tables["colors"]
        )
        colors.click()

        features = page.get_by_label("Features", exact=True)
        features.click()
        page.get_by_label("Rendering for CDS", exact=True).select_option("arrow")
        page.get_by_label("Rendering for rRNA", exact=True).select_option(
            "rectangle"
        )
        page.get_by_label("Rendering for tRNA", exact=True).select_option("arrow")
        page.get_by_label("Arrow head length ratio", exact=True).fill("")
        page.get_by_label("Arrow shaft width ratio", exact=True).fill("0.72")
        page.get_by_label("Block stroke color mode", exact=True).select_option(
            "color"
        )
        page.get_by_label("Block stroke color", exact=True).fill("#1F2937")
        page.get_by_label("Block Stroke Width", exact=True).fill("1.5")
        page.get_by_label("Line stroke color mode", exact=True).select_option(
            "color"
        )
        page.get_by_label("Line stroke color", exact=True).fill("#9CA3AF")
        page.get_by_label("Line Stroke Width", exact=True).fill("1.5")
        features.click()

        labels = page.get_by_label("Labels", exact=True)
        labels.click()
        page.get_by_label("Label Mode", exact=True).select_option("both")
        page.get_by_role("button", name="Whitelist", exact=True).click()
        page.get_by_label("Whitelist File", exact=True).set_input_files(
            tables["whitelist"]
        )
        page.get_by_label("Priority File (TSV)", exact=True).set_input_files(
            GENE_PRIORITY_RULE_PATH
        )
        page.get_by_text("Label Rendering", exact=True).locator(
            "xpath=.."
        ).locator("select").select_option("auto")
        labels.click()

        annotations = page.get_by_label("Region Annotations", exact=True)
        annotations.click()
        page.get_by_label("Import TSV", exact=True).set_input_files(
            tables["regions"]
        )
        expect(page.get_by_label("Annotation set id", exact=True)).to_have_value(
            "mitochondrial_regions"
        )
        annotations.click()
        slots = _configure_slots(page)

        axis = page.locator("summary").filter(has_text="Axis & Scale")
        axis.click()
        page.get_by_label("Axis stroke color mode", exact=True).select_option(
            "color"
        )
        page.get_by_label("Axis stroke color", exact=True).fill("#374151")
        page.get_by_label("Axis Stroke Width", exact=True).fill("4")
        axis.click()

        title = page.get_by_label("Title & Legend", exact=True)
        title.click()
        page.get_by_label("Plot Title", exact=True).fill(TITLE)
        page.get_by_label("Plot Title Position", exact=True).select_option("top")
        page.get_by_label("Legend Position", exact=True).select_option("right")
        title.click()

        page.get_by_role(
            "group", name="Circular track slot mitochondrial_regions", exact=True
        ).scroll_into_view_if_needed()
        screenshot_bytes[SCREENSHOT_NAMES[0]] = capture_screenshot(
            page, output_paths[SCREENSHOT_NAMES[0]], "Circular"
        )

        page.get_by_role("button", name="Generate Diagram", exact=True).click()
        result_region = page.get_by_role("region", name="Result Preview", exact=True)
        expect(result_region.locator("svg")).to_have_count(1, timeout=180_000)
        _load_label_overrides(page, tables["overrides"])
        final_report = _merge_browser_report(result_region)
        _assert_control_state(final_report)
        _assert_highlight_svg(final_report)

        _park_feature_search(page)
        _fit_finished_circular_preview(page)
        page.get_by_role("button", name="Zoom out", exact=True).click()
        expect(page.get_by_role("button", name="Reset zoom", exact=True)).to_contain_text(
            "50%"
        )
        _pan_preview_left(page, 0.28)
        _wait_for_preview_transform(page)
        screenshot_bytes[SCREENSHOT_NAMES[1]] = capture_screenshot(
            page, output_paths[SCREENSHOT_NAMES[1]], "Circular"
        )
        download_report = _download_svg(page, download_dir)
        capture.assert_clean()
    finally:
        capture.close()

    return CaptureResult(
        screenshot_bytes=screenshot_bytes,
        final_svg_semantics=final_report,
        download=download_report,
        source_record={**source_record, "trackSlots": slots},
    )


__all__ = ["capture_gui_feature_highlight"]
