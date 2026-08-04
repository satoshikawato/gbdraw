"""Capture the interactive-SVG and saved-session handoff tutorial."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

from playwright.sync_api import BrowserType, expect

from assertions.svg_semantics import (
    assert_finished_circular_svg,
    inspect_first_circular_svg,
    inspect_svg_file,
)
from config import (
    GUI_INTERACTIVE_HANDOFF_SCREENSHOT_NAMES,
)
from flows.how_to.exports import (
    _download,
    _validate_interactive_behavior,
    _validate_interactive_svg,
)
from flows.how_to.interactive_sessions import (
    _download_svg,
    _frame_finished_preview_with_legend,
    _load_current_session,
    _reset_finished_preview_viewport,
    _save_current_session,
    _search_and_open_feature,
    _stabilize_static_capture_surface,
    _validate_current_session,
)
from flows.human_circular import (
    assert_human_mitochondrion_fixture,
    generate_finished_human_diagram,
    load_raw_human_circular,
)
from flows.web_capture import (
    assert_output_paths,
    capture_screenshot,
    open_browser_capture,
    set_feature_search_visible,
    wait_for_worker,
)


SCENARIO_ID = "T-GUI-09"
INTERACTIVE_PREFIX = "interactive_human_mitochondrion"
INTERACTIVE_SVG = f"{INTERACTIVE_PREFIX}.interactive.svg"
SESSION_TITLE = "interactive_handoff"
SESSION_FILE = f"{SESSION_TITLE}.gbdraw-session.json.gz"
RESTORED_PREFIX = "restored_interactive_figure"
RESTORED_SVG = f"{RESTORED_PREFIX}.svg"


@dataclass(frozen=True)
class InteractiveHandoffResult:
    screenshot_bytes: dict[str, int]
    final_svg_semantics: dict[str, Any]
    download: dict[str, Any]
    downloads: dict[str, dict[str, Any]]
    session_download: dict[str, Any]
    source_record: dict[str, Any]
    feature_popup: dict[str, Any]


def _preview_dimensions(page: Any) -> tuple[float, float]:
    values = page.get_by_role(
        "region", name="Result Preview", exact=True
    ).locator("svg").first.evaluate(
        """
        (svg) => [parseFloat(svg.getAttribute('width')), parseFloat(svg.getAttribute('height'))]
        """
    )
    dimensions = tuple(float(value) for value in values)
    if len(dimensions) != 2 or any(value <= 0 for value in dimensions):
        raise AssertionError(f"Preview has invalid dimensions: {dimensions!r}")
    return dimensions


def _assert_reproduced(
    expected: Mapping[str, Any],
    actual: Mapping[str, Any],
    *,
    compare_label_classification: bool = True,
) -> None:
    """Compare durable figure semantics while tolerating subpixel text metrics."""

    keys = [
        "recordIds",
        "featureIds",
        "texts",
        "groupChildCounts",
        "recordTranslate",
    ]
    if compare_label_classification:
        keys.append("labelTexts")
    for key in keys:
        if expected.get(key) != actual.get(key):
            raise AssertionError(f"Session reproduction changed {key}")
    expected_legend = expected.get("legendTranslate")
    actual_legend = actual.get("legendTranslate")
    if (
        not expected_legend
        or not actual_legend
        or any(
            abs(float(left) - float(right)) >= 1
            for left, right in zip(expected_legend, actual_legend, strict=True)
        )
    ):
        raise AssertionError("Session reproduction moved the legend materially")


def capture_gui_interactive_handoff(
    browser_type: BrowserType,
    base_url: str,
    output_paths: Mapping[str, Path],
    download_dir: Path,
) -> InteractiveHandoffResult:
    """Export an offline interactive SVG and reproduce its session from scratch."""

    source_record = assert_human_mitochondrion_fixture()
    assert_output_paths(
        output_paths,
        GUI_INTERACTIVE_HANDOFF_SCREENSHOT_NAMES,
        SCENARIO_ID,
    )
    download_dir.mkdir(parents=True, exist_ok=True)
    screenshot_bytes: dict[str, int] = {}

    first_capture = open_browser_capture(browser_type, base_url)
    page = first_capture.page
    try:
        page.goto(base_url, wait_until="domcontentloaded")
        wait_for_worker(page)
        _stabilize_static_capture_surface(page)
        load_raw_human_circular(page, output_prefix=INTERACTIVE_PREFIX)
        screenshot_bytes[GUI_INTERACTIVE_HANDOFF_SCREENSHOT_NAMES[0]] = (
            capture_screenshot(
                page,
                output_paths[GUI_INTERACTIVE_HANDOFF_SCREENSHOT_NAMES[0]],
                "Circular",
            )
        )

        source_report = generate_finished_human_diagram(page)
        _reset_finished_preview_viewport(page, target_zoom=60)
        set_feature_search_visible(page, visible=False)
        screenshot_bytes[GUI_INTERACTIVE_HANDOFF_SCREENSHOT_NAMES[1]] = (
            capture_screenshot(
                page,
                output_paths[GUI_INTERACTIVE_HANDOFF_SCREENSHOT_NAMES[1]],
                "Circular",
            )
        )

        dimensions = _preview_dimensions(page)
        interactive_button = page.get_by_role(
            "button", name=re.compile(r"Interactive SVG$")
        )
        interactive_button.scroll_into_view_if_needed()
        screenshot_bytes[GUI_INTERACTIVE_HANDOFF_SCREENSHOT_NAMES[2]] = (
            capture_screenshot(
                page,
                output_paths[GUI_INTERACTIVE_HANDOFF_SCREENSHOT_NAMES[2]],
                "Circular",
            )
        )
        interactive_path = _download(
            page,
            button_name="Interactive SVG",
            expected_name=INTERACTIVE_SVG,
            download_dir=download_dir,
        )
        interactive_report = _validate_interactive_svg(
            interactive_path,
            expected_dimensions=dimensions,
        )
        interactive_report["behavior"] = _validate_interactive_behavior(
            first_capture.browser,
            interactive_path,
        )

        set_feature_search_visible(page, visible=True)
        popup = _search_and_open_feature(
            page,
            query="COX1",
            field="qualifier-value",
            qualifier_key="gene",
            expected_matches=1,
        )
        popup_text = popup.inner_text()
        searchable_text = popup_text + page.get_by_role(
            "region", name="Result Preview", exact=True
        ).inner_text()
        if "COX1" not in searchable_text:
            raise AssertionError("COX1 search did not open the matching feature")
        if "location" not in popup_text.casefold():
            raise AssertionError("COX1 popup is missing its location")
        screenshot_bytes[GUI_INTERACTIVE_HANDOFF_SCREENSHOT_NAMES[3]] = (
            capture_screenshot(
                page,
                output_paths[GUI_INTERACTIVE_HANDOFF_SCREENSHOT_NAMES[3]],
                "Circular",
            )
        )
        popup.get_by_role("button", name="Close feature popup", exact=True).click()
        page.get_by_role("button", name="Clear search", exact=True).click()

        session_path = _save_current_session(
            page,
            download_dir,
            title=SESSION_TITLE,
        )
        if session_path.name != SESSION_FILE:
            raise AssertionError(f"Unexpected saved-session name: {session_path.name}")
        session_report = _validate_current_session(
            session_path,
            expected_title=SESSION_TITLE,
            expected_output_prefix=INTERACTIVE_PREFIX,
        )
        screenshot_bytes[GUI_INTERACTIVE_HANDOFF_SCREENSHOT_NAMES[4]] = (
            capture_screenshot(
                page,
                output_paths[GUI_INTERACTIVE_HANDOFF_SCREENSHOT_NAMES[4]],
                "Circular",
            )
        )
        first_capture.assert_clean()
    finally:
        first_capture.close()

    second_capture = open_browser_capture(browser_type, base_url)
    page = second_capture.page
    try:
        page.goto(base_url, wait_until="domcontentloaded")
        wait_for_worker(page)
        _stabilize_static_capture_surface(page)
        _load_current_session(page, session_path)
        restored_region = page.get_by_role(
            "region", name="Result Preview", exact=True
        )
        expect(restored_region).to_be_visible()
        restored_report = inspect_first_circular_svg(restored_region)
        assert_finished_circular_svg(restored_report)
        _assert_reproduced(source_report, restored_report)

        prefix = page.get_by_label("Output Prefix", exact=True)
        prefix.fill(RESTORED_PREFIX)
        prefix.press("Tab")
        final_report = generate_finished_human_diagram(page)
        _assert_reproduced(source_report, final_report)
        _reset_finished_preview_viewport(page, target_zoom=50)
        final_report["restoredPreviewFrame"] = _frame_finished_preview_with_legend(
            page
        )
        set_feature_search_visible(page, visible=False)
        screenshot_bytes[GUI_INTERACTIVE_HANDOFF_SCREENSHOT_NAMES[5]] = (
            capture_screenshot(
                page,
                output_paths[GUI_INTERACTIVE_HANDOFF_SCREENSHOT_NAMES[5]],
                "Circular",
            )
        )
        restored_path = _download_svg(
            page,
            download_dir=download_dir,
            expected_name=RESTORED_SVG,
        )
        exported = inspect_svg_file(restored_path)
        assert_finished_circular_svg(exported)
        _assert_reproduced(
            source_report,
            exported,
            compare_label_classification=False,
        )
        download_report = {
            "filename": restored_path.name,
            "bytes": restored_path.stat().st_size,
            "semantics": exported,
        }
        second_capture.assert_clean()
    finally:
        second_capture.close()

    return InteractiveHandoffResult(
        screenshot_bytes=screenshot_bytes,
        final_svg_semantics=final_report,
        download=download_report,
        downloads={
            "interactiveSvg": interactive_report,
            "restoredSvg": download_report,
        },
        session_download=session_report,
        source_record=source_record,
        feature_popup={"query": "COX1", "text": popup_text},
    )


__all__ = ["capture_gui_interactive_handoff"]
