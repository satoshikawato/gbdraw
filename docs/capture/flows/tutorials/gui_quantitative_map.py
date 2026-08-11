"""Capture the GUI variant of the quantitative AP027133 genome map."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Mapping

from playwright.sync_api import BrowserType, Page

from flows.how_to.tracks import (
    GUI_QUANTITATIVE_DEPTH_PATH,
    GUI_QUANTITATIVE_GENBANK_PATH,
    GUI_QUANTITATIVE_GENBANK_SHA256,
    GUI_QUANTITATIVE_GENBANK_SIZE,
    _download_svg,
    _fit_circular_preview,
    _inspect_tracks_svg,
    _resize_sidebar,
    _track_slot_snapshot,
    _validate_complete_record,
    _validate_depth_fixture,
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


SCENARIO_ID = "T-GUI-12"
SCREENSHOT_NAMES = ("track-settings.png", "track-result.png")
OUTPUT_PREFIX = "quantitative_genome_map"
OUTPUT_FILENAME = f"{OUTPUT_PREFIX}.svg"
EXPECTED_SLOTS = (
    ("ticks", "ticks"),
    ("features", "features"),
    ("depth_1", "depth"),
    ("gc_content", "dinucleotide_content"),
    ("gc_skew", "dinucleotide_skew"),
)


def _slot_pairs(report: Mapping[str, Any]) -> tuple[tuple[str, str], ...]:
    return tuple(
        (slot["slotId"], slot["renderer"])
        for slot in report.get("slots", [])
        if not str(slot["slotId"]).startswith("__gbdraw_auto_")
    )


def _assert_quantitative_map(report: Mapping[str, Any]) -> None:
    if report.get("svgCount") != 1:
        raise AssertionError("Expected one quantitative-map SVG")
    if report.get("scriptCount") != 0 or report.get("eventAttributes"):
        raise AssertionError("The static quantitative map contains active content")
    if report.get("unsafeHrefs"):
        raise AssertionError("The static quantitative map contains an unsafe reference")
    if "AP027133.1" not in set(report.get("recordIds", [])):
        raise AssertionError("The quantitative map lost complete AP027133.1")
    if report.get("featureElementCount") != 576:
        raise AssertionError(
            f"Expected 576 AP027133.1 features, found {report.get('featureElementCount')}"
        )
    if _slot_pairs(report) != EXPECTED_SLOTS:
        raise AssertionError(f"Unexpected quantitative-map slots: {_slot_pairs(report)!r}")

    by_renderer = {slot["renderer"]: slot for slot in report.get("slots", [])}
    if set(by_renderer["depth"]["tickTexts"]) != {
        "0x",
        "20x",
        "40x",
        "60x",
        "80x",
    }:
        raise AssertionError(
            f"Unexpected depth ticks: {by_renderer['depth']['tickTexts']!r}"
        )
    if set(by_renderer["dinucleotide_content"]["tickTexts"]) != {
        "10%",
        "20%",
        "30%",
        "40%",
        "50%",
        "55%",
    }:
        raise AssertionError(
            "Unexpected GC-content ticks: "
            f"{by_renderer['dinucleotide_content']['tickTexts']!r}"
        )
    required = {
        "AP027133.1",
        "606,194 bp",
        "DRR394922 mean depth",
        "GC content (%)",
        "GC skew (+)",
        "GC skew (-)",
    }
    missing = required - set(report.get("texts", []))
    if missing:
        raise AssertionError(f"Missing quantitative-map text: {sorted(missing)}")


def _number_control(section: Any, label: str) -> Any:
    return section.get_by_text(label, exact=True).locator("xpath=..").locator(
        'input[type="number"]'
    )


def _ensure_side(slot: Any, title: str) -> None:
    button = slot.get_by_title(title, exact=True)
    if button.is_enabled():
        button.click()


def _configure_slots(page: Page) -> tuple[dict[str, Any], ...]:
    page.get_by_role(
        "button", name=re.compile(r"Custom Track Slots$"), exact=False
    ).click()
    page.get_by_role("checkbox", name="Use custom stack", exact=True).check()

    ticks = page.get_by_role(
        "group", name="Circular track slot ticks", exact=True
    )
    _ensure_side(ticks, "Move outside Axis")
    ticks = page.get_by_role(
        "group", name="Circular track slot ticks", exact=True
    )
    ticks.locator("select").last.select_option("label_out_tick_in")

    features = page.get_by_role(
        "group", name="Circular track slot features", exact=True
    )
    _ensure_side(features, "Track on Axis")
    features = page.get_by_role(
        "group", name="Circular track slot features", exact=True
    )
    features.locator("select").last.select_option("split")

    depth = page.get_by_role(
        "group", name="Circular track slot depth", exact=True
    )
    _ensure_side(depth, "Move inside Axis")
    depth = page.get_by_role(
        "group", name="Circular track slot depth", exact=True
    )
    depth.get_by_label("Circular track slot id depth", exact=True).fill("depth_1")
    depth = page.get_by_role(
        "group", name="Circular track slot depth_1", exact=True
    )
    depth.get_by_title("Width", exact=True).fill("52px")

    gc_content = page.get_by_role(
        "group", name="Circular track slot gc_content", exact=True
    )
    _ensure_side(gc_content, "Move inside Axis")
    gc_content = page.get_by_role(
        "group", name="Circular track slot gc_content", exact=True
    )
    gc_content.get_by_title("Width", exact=True).fill("42px")
    gc_content.get_by_label("Track dinucleotide", exact=True).fill("GC")
    gc_content.get_by_label("Track legend label", exact=True).fill(
        "GC content (%)"
    )

    gc_skew = page.get_by_role(
        "group", name="Circular track slot gc_skew", exact=True
    )
    _ensure_side(gc_skew, "Move inside Axis")
    gc_skew = page.get_by_role(
        "group", name="Circular track slot gc_skew", exact=True
    )
    gc_skew.get_by_title("Width", exact=True).fill("34px")
    gc_skew.get_by_label("Track dinucleotide", exact=True).fill("GC")
    gc_skew.get_by_label("Track legend label", exact=True).fill("GC skew")

    slots = _track_slot_snapshot(page)
    actual = tuple(
        (slot["id"], slot["renderer"])
        for slot in slots
        if slot["enabled"]
    )
    if actual != EXPECTED_SLOTS:
        raise AssertionError(f"Unexpected quantitative-map slot state: {actual!r}")
    widths = {
        slot["id"]: str(slot.get("width") or "")
        for slot in page.evaluate(
            """
            () => (window.__GBDRAW_APP__?.adv?.circular_track_slots || [])
              .map((slot) => ({ id: String(slot?.id || ''), width: slot?.width }))
            """
        )
    }
    if (
        widths.get("depth_1") != "52px"
        or widths.get("gc_content") != "42px"
        or widths.get("gc_skew") != "34px"
    ):
        raise AssertionError(f"Unexpected quantitative-map widths: {widths!r}")
    return slots


def _capture_state(page: Page) -> dict[str, Any]:
    return page.evaluate(
        """
        () => {
          const app = window.__GBDRAW_APP__;
          const depth = app?.circularDepthTrackRows?.[0]?.config || app?.depthTracks?.[0] || {};
          return {
            prefix: String(app?.form?.prefix || ''),
            labels: String(app?.form?.labels_mode || ''),
            legend: String(app?.form?.legend || ''),
            separateStrands: Boolean(app?.form?.separate_strands),
            depthWindow: Number(app?.adv?.depth_window_size),
            depthStep: Number(app?.adv?.depth_step_size),
            depthMin: Number(app?.adv?.depth_min),
            depthMax: Number(app?.adv?.depth_max),
            depthLog: Boolean(app?.adv?.depth_normalize),
            depthAxis: Boolean(app?.adv?.depth_show_axis),
            depthTicks: Boolean(app?.adv?.depth_show_ticks),
            depthLabel: String(app?.getDepthTrackLabel?.(0) || depth?.label || ''),
            depthColor: String(app?.getDepthTrackColor?.(0) || depth?.color || '').toLowerCase(),
            depthLargeTick: Number(depth?.large_tick_interval),
            depthSmallTick: Number(depth?.small_tick_interval),
            window: Number(app?.adv?.window_size),
            step: Number(app?.adv?.step_size),
            gcMode: String(app?.adv?.gc_content_mode || ''),
            gcMin: Number(app?.adv?.gc_content_min_percent),
            gcMax: Number(app?.adv?.gc_content_max_percent),
            gcAxis: Boolean(app?.adv?.gc_content_show_axis),
            gcTicks: Boolean(app?.adv?.gc_content_show_ticks),
            gcLargeTick: Number(app?.adv?.gc_content_tick_interval),
            gcSmallTick: Number(app?.adv?.gc_content_small_tick_interval)
          };
        }
        """
    )


def _assert_state(state: Mapping[str, Any]) -> None:
    expected = {
        "prefix": OUTPUT_PREFIX,
        "labels": "none",
        "legend": "right",
        "separateStrands": True,
        "depthWindow": 1,
        "depthStep": 1000,
        "depthMin": 0,
        "depthMax": 80,
        "depthLog": False,
        "depthAxis": True,
        "depthTicks": True,
        "depthLabel": "DRR394922 mean depth",
        "depthColor": "#2563eb",
        "depthLargeTick": 20,
        "depthSmallTick": 10,
        "window": 1000,
        "step": 1000,
        "gcMode": "percent",
        "gcMin": 10,
        "gcMax": 55,
        "gcAxis": True,
        "gcTicks": True,
        "gcLargeTick": 10,
        "gcSmallTick": 5,
    }
    if dict(state) != expected:
        raise AssertionError(f"Unexpected quantitative-map control state: {state!r}")


def capture_gui_quantitative_map(
    browser_type: BrowserType,
    base_url: str,
    output_paths: Mapping[str, Path],
    download_dir: Path,
) -> CaptureResult:
    """Build the CLI/Python-equivalent depth, GC-content, and skew map."""

    genome_report = _validate_complete_record(
        GUI_QUANTITATIVE_GENBANK_PATH,
        expected_size=GUI_QUANTITATIVE_GENBANK_SIZE,
        expected_sha256=GUI_QUANTITATIVE_GENBANK_SHA256,
        expected_id="AP027133.1",
        expected_length=606_194,
    )
    depth_report = _validate_depth_fixture()
    assert_output_paths(output_paths, SCREENSHOT_NAMES, SCENARIO_ID)
    download_dir.mkdir(parents=True, exist_ok=True)

    capture = open_browser_capture(browser_type, base_url)
    page = capture.page
    screenshot_bytes: dict[str, int] = {}
    try:
        page.goto(base_url, wait_until="domcontentloaded")
        wait_for_worker(page)
        _resize_sidebar(page)
        page.get_by_role("button", name="Circular", exact=True).click()
        page.get_by_role("radio", name="GenBank", exact=True).check()
        page.get_by_label("GenBank/DDBJ File", exact=True).set_input_files(
            GUI_QUANTITATIVE_GENBANK_PATH
        )
        page.get_by_label("Output Prefix", exact=True).fill(OUTPUT_PREFIX)
        page.get_by_label("Separate Strands", exact=True).check()

        depth_section = page.get_by_label("Depth TSV tracks", exact=True)
        depth_section.click()
        page.get_by_label("Depth TSV", exact=True).set_input_files(
            GUI_QUANTITATIVE_DEPTH_PATH
        )
        page.get_by_label("Depth Window", exact=True).fill("1")
        page.get_by_label("Depth Step", exact=True).fill("1000")
        page.get_by_label("Depth Min", exact=True).fill("0")
        page.get_by_label("Depth Max", exact=True).fill("80")
        page.get_by_role("checkbox", name="Log Scale", exact=True).uncheck()
        page.get_by_role("checkbox", name="Depth Axis", exact=True).check()
        page.get_by_role("checkbox", name="Depth Ticks", exact=True).check()
        page.get_by_label("Depth legend title", exact=True).fill(
            "DRR394922 mean depth"
        )
        depth_section.locator("xpath=..").locator('input[type="color"]').first.fill(
            "#2563eb"
        )
        page.get_by_label("Depth large tick", exact=True).fill("20")
        page.get_by_label("Depth small tick", exact=True).fill("10")

        dinucleotide_summary = page.get_by_label(
            "Dinucleotide content/skew", exact=True
        )
        dinucleotide_summary.click()
        dinucleotide = dinucleotide_summary.locator("xpath=..")
        _number_control(dinucleotide, "Window").fill("1000")
        _number_control(dinucleotide, "Step").fill("1000")
        dinucleotide.locator('select:has(option[value="percent"])').select_option(
            "percent"
        )
        _number_control(dinucleotide, "GC Min %").fill("10")
        _number_control(dinucleotide, "GC Max %").fill("55")
        dinucleotide.get_by_role(
            "checkbox", name="Percent Axis", exact=True
        ).check()
        dinucleotide.get_by_role(
            "checkbox", name="Percent Ticks", exact=True
        ).check()
        _number_control(dinucleotide, "Large Tick").fill("10")
        _number_control(dinucleotide, "Small Tick").fill("5")

        title = page.get_by_label("Title & Legend", exact=True)
        title.click()
        page.get_by_label("Legend Position", exact=True).select_option("right")
        title.click()

        slots = _configure_slots(page)
        state = _capture_state(page)
        _assert_state(state)
        page.get_by_role(
            "group", name="Circular track slot depth_1", exact=True
        ).scroll_into_view_if_needed()
        screenshot_bytes[SCREENSHOT_NAMES[0]] = capture_screenshot(
            page, output_paths[SCREENSHOT_NAMES[0]], "Circular"
        )

        final_report = generate_and_inspect(
            page, _inspect_tracks_svg, _assert_quantitative_map
        )
        final_report["state"] = state
        _fit_circular_preview(page, target_zoom="70%", pan_left_ratio=0.16)
        set_feature_search_visible(page, visible=False)
        screenshot_bytes[SCREENSHOT_NAMES[1]] = capture_screenshot(
            page, output_paths[SCREENSHOT_NAMES[1]], "Circular"
        )
        download_report = _download_svg(
            page,
            download_dir,
            expected_filename=OUTPUT_FILENAME,
            assert_svg=_assert_quantitative_map,
        )
        capture.assert_clean()
    finally:
        capture.close()

    return CaptureResult(
        screenshot_bytes=screenshot_bytes,
        final_svg_semantics=final_report,
        download=download_report,
        fixture_report={**genome_report, "depth": depth_report, "state": state},
        track_slots=slots,
    )


__all__ = ["capture_gui_quantitative_map"]
