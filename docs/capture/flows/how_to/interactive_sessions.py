"""Capture finished-diagram editing and current-session reproduction journeys."""

from __future__ import annotations

import base64
import gzip
import hashlib
import json
import math
import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

from playwright.sync_api import BrowserType, Page, expect

from assertions.svg_semantics import (
    assert_finished_circular_svg,
    assert_first_circular_svg,
    assert_gui_bgc_similarity_groups_svg,
    assert_species_markup,
    inspect_first_circular_svg,
    inspect_gui_bgc_losatp_svg,
    inspect_svg_file,
)
from config import (
    ACTION_TIMEOUT_MS,
    FIRST_CIRCULAR_FIXTURE_SHA256,
    FIRST_CIRCULAR_FIXTURE_SIZE,
    GUI_INTERACTIVE_EDITING_SCREENSHOT_NAMES,
    GUI_SESSION_REPRODUCTION_SCREENSHOT_NAMES,
    VIEWPORT_HEIGHT,
    VIEWPORT_WIDTH,
)
from flows.bgc_losatp import (
    _configure_losatp,
    _set_bgc_inputs,
    _set_gallery_quality_presentation,
    assert_bgc_fixtures,
)
from flows.human_circular import (
    apply_finished_human_settings,
    assert_human_mitochondrion_fixture,
    fit_finished_human_preview,
    generate_finished_human_diagram,
    load_raw_human_circular,
)
from flows.web_capture import (
    assert_output_paths,
    capture_screenshot,
    fit_complete_linear_preview,
    generate_and_inspect,
    open_browser_capture,
    set_feature_search_visible,
    wait_for_worker,
)


EDITED_OUTPUT_PREFIX = "edited_diagram"
EDITED_SVG_NAME = f"{EDITED_OUTPUT_PREFIX}.svg"
EDIT_TARGET_QUERY = "COX"
EDIT_FILL = "#d81b60"
EDIT_STROKE = "#5e35b1"
EDIT_STROKE_WIDTH = 2.5
EDIT_LEGEND = "Oxidative phosphorylation"
SESSION_TITLE = "reproducible_work"
SESSION_FILENAME = f"{SESSION_TITLE}.gbdraw-session.json.gz"
RELOADED_OUTPUT_PREFIX = "reloaded_diagram"
RELOADED_SVG_NAME = f"{RELOADED_OUTPUT_PREFIX}.svg"
CURRENT_SESSION_VERSION = 41
CURRENT_RENDER_REQUEST_SCHEMA = 6
STATIC_CAPTURE_HEADER_STYLE = """
.app-header {
    -webkit-backdrop-filter: none !important;
    backdrop-filter: none !important;
}
"""


@dataclass(frozen=True)
class InteractiveEditingResult:
    """Artifacts and assertions from H-GUI-13."""

    screenshot_bytes: dict[str, int]
    final_svg_semantics: dict[str, Any]
    download: dict[str, Any]
    source_record: dict[str, Any]
    match_popup: dict[str, Any]
    group_popup: dict[str, Any]


@dataclass(frozen=True)
class SessionReproductionResult:
    """Artifacts and assertions from H-GUI-14."""

    screenshot_bytes: dict[str, int]
    final_svg_semantics: dict[str, Any]
    download: dict[str, Any]
    session_download: dict[str, Any]
    source_record: dict[str, Any]


def _stabilize_static_capture_surface(page: Page) -> None:
    """Remove the header compositor blur that makes static pixels nondeterministic."""

    page.add_style_tag(content=STATIC_CAPTURE_HEADER_STYLE)
    expect(page.locator(".app-header")).to_have_css("backdrop-filter", "none")
    page.evaluate(
        """
        async () => {
          if (document.fonts) await document.fonts.ready;
          await new Promise((resolve) => {
            requestAnimationFrame(() => requestAnimationFrame(resolve));
          });
        }
        """
    )


def _reset_finished_preview_viewport(page: Page, *, target_zoom: int) -> None:
    """Reset search-driven pan and reach a pinned zoom through public controls."""

    if target_zoom not in (50, 60):
        raise ValueError(f"Unsupported finished-preview zoom: {target_zoom}%")
    reset_zoom = page.get_by_role("button", name="Reset zoom", exact=True)
    page.get_by_role("button", name="Reset layout", exact=True).click()
    expect(reset_zoom).to_have_text("100%")
    zoom_out = page.get_by_role("button", name="Zoom out", exact=True)
    for expected_zoom in range(90, target_zoom - 1, -10):
        zoom_out.click()
        expect(reset_zoom).to_have_text(f"{expected_zoom}%")
    page.wait_for_timeout(250)
    viewport = page.evaluate(
        """
        () => ({
          panX: Number(window.__GBDRAW_APP__?.canvasPan?.x || 0),
          panY: Number(window.__GBDRAW_APP__?.canvasPan?.y || 0),
          zoom: Number(window.__GBDRAW_APP__?.zoom || 0)
        })
        """
    )
    if (
        not math.isclose(viewport["panX"], 0, abs_tol=0.01)
        or not math.isclose(viewport["panY"], 0, abs_tol=0.01)
        or not math.isclose(viewport["zoom"], target_zoom / 100, abs_tol=0.001)
    ):
        raise AssertionError(f"Finished preview did not reset cleanly: {viewport!r}")


def _feature_search_panel(page: Page) -> Any:
    return page.get_by_role(
        "searchbox", name="Search features", exact=True
    ).locator("xpath=ancestor::div[contains(@class, 'preview-feature-search')]")


def _assert_closed_right_drawer_contract(page: Page) -> None:
    """Check that the closed editor is neither painted nor interactive."""

    for width in (1280, VIEWPORT_WIDTH):
        page.set_viewport_size({"width": width, "height": VIEWPORT_HEIGHT})
        page.wait_for_timeout(350)
        drawer_state = page.locator(".right-drawer").evaluate(
            """
            (drawer) => {
              const style = getComputedStyle(drawer);
              return {
                appOpen: window.__GBDRAW_APP__?.showRightDrawer,
                ariaHidden: drawer.getAttribute('aria-hidden'),
                pointerEvents: style.pointerEvents,
                visibility: style.visibility,
                width: drawer.getBoundingClientRect().width,
                translateX: new DOMMatrixReadOnly(style.transform).m41,
              };
            }
            """
        )
        toggle_state = page.locator(".drawer-toggle").evaluate(
            """
            (toggle) => {
              const rect = toggle.getBoundingClientRect();
              const parentRect = toggle.offsetParent.getBoundingClientRect();
              const style = getComputedStyle(toggle);
              return {
                ariaExpanded: toggle.getAttribute('aria-expanded'),
                className: toggle.className,
                rightGap: parentRect.right - rect.right,
                translateX: new DOMMatrixReadOnly(style.transform).m41,
              };
            }
            """
        )
        if (
            drawer_state.get("appOpen") is not False
            or drawer_state.get("ariaHidden") != "true"
            or drawer_state.get("pointerEvents") != "none"
            or drawer_state.get("visibility") != "hidden"
            or drawer_state.get("translateX", 0) < drawer_state.get("width", 0) - 1
        ):
            raise AssertionError(
                f"The editor drawer is not fully closed at {width}px: "
                f"{drawer_state!r}"
            )
        if (
            toggle_state.get("ariaExpanded") != "false"
            or "drawer-open" in toggle_state.get("className", "")
            or "-translate-x-[" in toggle_state.get("className", "")
            or abs(toggle_state.get("translateX", 1)) > 1
            or abs(toggle_state.get("rightGap", 5)) > 4
        ):
            raise AssertionError(
                f"The closed editor toggle is displaced at {width}px: "
                f"{toggle_state!r}"
            )


def _assert_closed_drawer_does_not_reserve_canvas(page: Page) -> dict[str, float]:
    result_region = page.get_by_role("region", name="Result Preview", exact=True)
    canvas = result_region.locator(
        "div.absolute.inset-0.overflow-auto.flex.p-2"
    ).first
    svg = result_region.locator("svg").first
    legend = result_region.locator("#legend")
    legend_caption = legend.locator("text").filter(has_text=EDIT_LEGEND).first
    expect(legend_caption).to_be_visible()
    boxes = {
        "result": result_region.bounding_box(),
        "canvas": canvas.bounding_box(),
        "svg": svg.bounding_box(),
        "legend": legend.bounding_box(),
        "legendCaption": legend_caption.bounding_box(),
    }
    if any(box is None for box in boxes.values()):
        raise AssertionError(f"Could not measure the finished preview: {boxes!r}")
    result_box = boxes["result"]
    canvas_box = boxes["canvas"]
    svg_box = boxes["svg"]
    legend_box = boxes["legend"]
    legend_caption_box = boxes["legendCaption"]
    assert result_box is not None
    assert canvas_box is not None
    assert svg_box is not None
    assert legend_box is not None
    assert legend_caption_box is not None
    geometry = {
        "canvasRightGap": result_box["x"] + result_box["width"]
        - canvas_box["x"]
        - canvas_box["width"],
        "svgWidth": svg_box["width"],
        "legendRight": legend_box["x"] + legend_box["width"],
        "legendLeft": legend_box["x"],
        "legendCaptionRight": legend_caption_box["x"]
        + legend_caption_box["width"],
        "resultRight": result_box["x"] + result_box["width"],
        "resultLeft": result_box["x"],
    }
    if (
        geometry["canvasRightGap"] > 24
        or geometry["svgWidth"] < 760
        or geometry["legendLeft"] < geometry["resultLeft"] + 4
        or geometry["legendRight"] > geometry["resultRight"] - 4
        or geometry["legendCaptionRight"] > geometry["resultRight"] - 4
    ):
        raise AssertionError(
            "The closed drawer still reserves or clips finished-preview space: "
            f"{geometry!r}"
        )
    return geometry


def _search_and_open_feature(
    page: Page,
    *,
    query: str,
    field: str,
    qualifier_key: str = "",
    use_regex: bool = False,
    expected_matches: int,
) -> Any:
    """Search and open one rendered feature through the preview controls."""

    search = page.get_by_role("searchbox", name="Search features", exact=True)
    search.fill(query)
    field_control = page.get_by_role(
        "combobox", name="Search field", exact=True
    )
    field_control.select_option(field)
    expect(field_control).to_have_value(field)
    qualifier = page.get_by_role(
        "textbox", name="Qualifier key", exact=True
    )
    if qualifier_key:
        expect(qualifier).to_be_enabled()
        qualifier.fill(qualifier_key)
    regex = page.get_by_role("checkbox", name="Regex", exact=True)
    if use_regex:
        regex.check()
        expect(regex).to_be_checked()
    page.get_by_role("button", name="Search features", exact=True).click()
    expect(
        page.get_by_role("status", name="Feature search status", exact=True)
    ).to_have_text(f"1 / {expected_matches} features")
    page.get_by_role("button", name="Open active feature", exact=True).click()
    popup = page.get_by_role("dialog", name=re.compile(r"^Feature details:"))
    expect(popup).to_be_visible()
    page.wait_for_timeout(250)
    return popup


def _active_feature_id(page: Page) -> str:
    active = page.get_by_role("region", name="Result Preview", exact=True)
    feature_id = active.locator(
        ".gbdraw-preview-feature-search-active-match"
    ).first.get_attribute("data-gbdraw-feature-id")
    if not feature_id:
        raise AssertionError("Feature search did not expose an active feature ID")
    return feature_id.split("__part", 1)[0]


def _click_active_feature_with_control(page: Page, feature_id: str) -> None:
    """Ctrl-click an unobscured point on the active rendered feature."""

    result_region = page.get_by_role("region", name="Result Preview", exact=True)
    point = result_region.evaluate(
        r"""
        (root, expectedId) => {
          const normalizeId = (value) => String(value || '').replace(/__part\d+$/, '');
          const matchesExpected = (hit) => {
            const target = hit?.closest?.('[data-gbdraw-feature-id]');
            return normalizeId(target?.getAttribute('data-gbdraw-feature-id')) === expectedId;
          };
          const candidates = Array.from(
            root.querySelectorAll('.gbdraw-preview-feature-search-active-match')
          );
          for (const element of candidates) {
            if (
              typeof element.getTotalLength === 'function' &&
              typeof element.getPointAtLength === 'function'
            ) {
              const length = Number(element.getTotalLength()) || 0;
              const matrix = element.getScreenCTM?.();
              if (length > 0 && matrix) {
                for (let index = 1; index < 40; index += 1) {
                  const local = element.getPointAtLength(length * index / 40);
                  const screen = new DOMPoint(local.x, local.y).matrixTransform(matrix);
                  if (matchesExpected(document.elementFromPoint(screen.x, screen.y))) {
                    return { x: screen.x, y: screen.y };
                  }
                }
              }
            }
            const rect = element.getBoundingClientRect();
            for (const xRatio of [0.2, 0.5, 0.8]) {
              for (const yRatio of [0.2, 0.5, 0.8]) {
                const x = rect.left + rect.width * xRatio;
                const y = rect.top + rect.height * yRatio;
                if (matchesExpected(document.elementFromPoint(x, y))) return { x, y };
              }
            }
          }
          return null;
        }
        """,
        feature_id,
    )
    if not point:
        raise AssertionError(
            f"Could not find an unobscured click point for feature {feature_id}"
        )
    page.mouse.click(float(point["x"]), float(point["y"]))


def _select_all_search_matches(page: Page, expected_matches: int) -> tuple[str, ...]:
    """Select every current search match with visible Ctrl-click interactions."""

    selected_ids: list[str] = []
    page.keyboard.down("Control")
    try:
        for index in range(expected_matches):
            feature_id = _active_feature_id(page)
            selected_ids.append(feature_id)
            _click_active_feature_with_control(page, feature_id)
            expect(
                page.get_by_role(
                    "status", name="Selected feature count", exact=True
                )
            ).to_have_text(f"{index + 1} selected")
            if index + 1 < expected_matches:
                page.get_by_role(
                    "button", name="Next match", exact=True
                ).click()
    finally:
        page.keyboard.up("Control")
    expect(
        page.get_by_role("status", name="Selected feature count", exact=True)
    ).to_have_text(f"{expected_matches} selected")
    return tuple(selected_ids)


def _edit_selected_features(page: Page, feature_ids: tuple[str, ...]) -> None:
    """Apply fill/stroke and stage visibility in the documented bulk editor."""

    result_region = page.get_by_role("region", name="Result Preview", exact=True)
    selected_glyph = result_region.locator(
        f'[data-gbdraw-feature-id="{feature_ids[0]}"]'
    ).first
    fill = page.get_by_label("Selected feature color", exact=True)
    fill.fill(EDIT_FILL)
    expect(fill).to_have_value(EDIT_FILL)
    caption = page.get_by_role(
        "textbox", name="Selected feature legend caption", exact=True
    )
    caption.fill(EDIT_LEGEND)
    page.get_by_role(
        "button",
        name="Apply color and caption to selected features",
        exact=True,
    ).click()
    expect(result_region.locator("#legend")).to_contain_text(
        EDIT_LEGEND, timeout=ACTION_TIMEOUT_MS
    )
    expect(result_region.locator(".gbdraw-feature-selected")).to_have_count(
        3, timeout=ACTION_TIMEOUT_MS
    )
    # Color rules are committed immediately, but the current SVG can retain its
    # pre-rule fill until the next render.  The finished result and downloaded
    # SVG are asserted after a real Generate action below.

    stroke = page.get_by_label("Selected feature stroke color", exact=True)
    stroke.fill(EDIT_STROKE)
    stroke_width = page.get_by_role(
        "spinbutton", name="Selected feature stroke width", exact=True
    )
    stroke_width.fill(str(EDIT_STROKE_WIDTH))
    page.get_by_role(
        "button", name="Apply stroke to selected features", exact=True
    ).click()
    expect(selected_glyph).to_have_attribute(
        "stroke", EDIT_STROKE, timeout=ACTION_TIMEOUT_MS
    )
    expect(selected_glyph).to_have_attribute(
        "stroke-width", str(EDIT_STROKE_WIDTH), timeout=ACTION_TIMEOUT_MS
    )

    visibility = page.get_by_role(
        "combobox", name="Selected feature visibility", exact=True
    )
    visibility.select_option("on")
    expect(visibility).to_have_value("on")
    page.wait_for_timeout(350)


def _inspect_target_style(result_region: Any, feature_id: str) -> dict[str, Any]:
    return result_region.evaluate(
        r"""
        (root, featureId) => {
          const svg = root.getElementsByTagName('svg')[0];
          const elements = Array.from(
            svg?.querySelectorAll('[data-gbdraw-feature-id]') || []
          ).filter((element) => String(
            element.getAttribute('data-gbdraw-feature-id') || ''
          ).replace(/__part\d+$/, '') === featureId);
          const legend = svg?.querySelector('#legend');
          return {
            featureId,
            styles: elements.map((element) => ({
              fill: String(element.getAttribute('fill') || '').toLowerCase(),
              stroke: String(element.getAttribute('stroke') || '').toLowerCase(),
              strokeWidth: Number(element.getAttribute('stroke-width') || 0)
            })),
            legendTransform: String(legend?.getAttribute('transform') || '')
          };
        }
        """,
        feature_id,
    )


def _assert_target_edit(report: dict[str, Any]) -> None:
    styles = report.get("styles", [])
    if not styles:
        raise AssertionError("An edited COX feature is absent from the SVG")
    if not any(style.get("fill") == EDIT_FILL for style in styles):
        raise AssertionError(f"COX fill did not remain {EDIT_FILL}: {styles!r}")
    if not any(
        style.get("stroke") == EDIT_STROKE
        and math.isclose(float(style.get("strokeWidth") or 0), EDIT_STROKE_WIDTH)
        for style in styles
    ):
        raise AssertionError(
            f"COX stroke did not remain {EDIT_STROKE} at {EDIT_STROKE_WIDTH} px"
        )


def _assert_edited_human_svg(report: dict[str, Any]) -> None:
    """Validate the finished state after selected CDS leave the base category."""

    assert_first_circular_svg(report)
    assert_species_markup(report)
    ids = set(report.get("ids", []))
    if not {"label_leaders", "label_text", "legend"}.issubset(ids):
        raise AssertionError("The edited diagram lost labels or its legend")
    counts = report.get("groupChildCounts", {})
    if counts.get("label_leaders", 0) < 10 or counts.get("label_text", 0) < 10:
        raise AssertionError(f"Edited external labels are incomplete: {counts}")
    texts = set(report.get("texts", []))
    expected_legend = {
        "tRNA",
        "rRNA",
        EDIT_LEGEND,
        "other proteins",
        "GC content",
        "GC skew (+)",
        "GC skew (-)",
    }
    if not expected_legend.issubset(texts):
        raise AssertionError(
            f"Edited legend entries are incomplete: {sorted(expected_legend - texts)}"
        )
    record_translate = report.get("recordTranslate")
    legend_translate = report.get("legendTranslate")
    if (
        not record_translate
        or not legend_translate
        or float(legend_translate[0]) <= float(record_translate[0])
    ):
        raise AssertionError("The edited legend is not right of the Circular record")


def _drag_legend(page: Page) -> tuple[str, str]:
    result_region = page.get_by_role("region", name="Result Preview", exact=True)
    before = result_region.locator("#legend").get_attribute("transform") or ""
    toggle = page.get_by_role(
        "button", name="Toggle layout edit mode", exact=True
    )
    toggle.click()
    expect(toggle).to_have_attribute("aria-pressed", "true")
    legend = result_region.locator("#legend")
    expect(legend).to_have_css("cursor", "grab", timeout=ACTION_TIMEOUT_MS)
    handle = legend.evaluate(
        """
        (legend) => {
          for (const element of [legend, ...legend.querySelectorAll('*')]) {
            const rect = element.getBoundingClientRect();
            if (rect.width <= 2 || rect.height <= 2) continue;
            const points = [
              [rect.left + rect.width / 2, rect.top + rect.height / 2],
              [rect.left + 2, rect.top + rect.height / 2],
              [rect.right - 2, rect.top + rect.height / 2],
            ];
            for (const [x, y] of points) {
              const hit = document.elementFromPoint(x, y);
              if (hit && legend.contains(hit)) {
                return { x, y, tag: hit.tagName, text: hit.textContent || '' };
              }
            }
          }
          const svgElement = legend.ownerSVGElement;
          const svg = svgElement?.getBoundingClientRect();
          const rect = legend.getBoundingClientRect();
          const sampleX = rect.left + rect.width / 2;
          const sampleY = rect.top + rect.height / 2;
          const sampleHit = document.elementFromPoint(sampleX, sampleY);
          const serialize = (value) => value ? ({
            left: value.left,
            top: value.top,
            right: value.right,
            bottom: value.bottom,
            width: value.width,
            height: value.height,
          }) : null;
          return {
            error: 'no hit target',
            legend: serialize(rect),
            legendPointerEvents: getComputedStyle(legend).pointerEvents,
            svg: serialize(svg),
            viewBox: svgElement?.getAttribute('viewBox') || '',
            sampleHit: sampleHit ? {
              tag: sampleHit.tagName,
              id: sampleHit.id || '',
              className: String(sampleHit.className || ''),
              ariaHidden: sampleHit.getAttribute('aria-hidden'),
              pointerEvents: getComputedStyle(sampleHit).pointerEvents,
              visibility: getComputedStyle(sampleHit).visibility,
              zIndex: getComputedStyle(sampleHit).zIndex,
            } : null,
          };
        }
        """
    )
    if handle.get("error"):
        raise AssertionError(f"The finished legend is not draggable: {handle!r}")
    start_x = handle["x"]
    start_y = handle["y"]
    page.mouse.move(start_x, start_y)
    page.mouse.down()
    drag_start = page.evaluate(
        """
        ({ x, y }) => {
          const hit = document.elementFromPoint(x, y);
          const legend = hit?.closest?.('#legend');
          return {
            active: legend?.style?.willChange === 'transform',
            hitId: hit?.id || '',
            hitTag: hit?.tagName || '',
            legendAncestor: legend?.id || '',
          };
        }
        """,
        {"x": start_x, "y": start_y},
    )
    if drag_start.get("active") is not True:
        page.mouse.up()
        raise AssertionError(
            f"The public legend drag did not start at {handle!r}: {drag_start!r}"
        )
    page.mouse.move(start_x - 70, start_y - 32, steps=10)
    page.mouse.up()
    page.wait_for_timeout(250)
    after = legend.get_attribute("transform") or ""
    if after == before:
        raise AssertionError("Dragging the legend did not change its SVG transform")
    toggle.click()
    expect(toggle).to_have_attribute("aria-pressed", "false")
    return before, after


def _scroll_finished_preview_to_right(page: Page) -> dict[str, float]:
    """Scroll the public preview viewport until the right-side legend is visible."""

    canvas = page.get_by_role(
        "region", name="Result Preview", exact=True
    ).locator("div.absolute.inset-0.overflow-auto.flex.p-2").first
    metrics = canvas.evaluate(
        """
        (element) => {
          const maxScroll = Math.max(0, element.scrollWidth - element.clientWidth);
          element.scrollLeft = maxScroll * 0.55;
          return {
            clientWidth: element.clientWidth,
            scrollWidth: element.scrollWidth,
            scrollLeft: element.scrollLeft,
          };
        }
        """
    )
    max_scroll = float(metrics["scrollWidth"] - metrics["clientWidth"])
    target_scroll = max_scroll * 0.55
    if max_scroll < 1 or abs(float(metrics["scrollLeft"]) - target_scroll) > 2:
        raise AssertionError(
            f"Could not expose the right-side legend in the preview: {metrics!r}"
        )
    return {key: float(value) for key, value in metrics.items()}


def _frame_finished_preview_with_legend(page: Page) -> dict[str, float]:
    """Scroll a reset 50 percent preview so the record and legend both fit."""

    canvas = page.get_by_role(
        "region", name="Result Preview", exact=True
    ).locator("div.absolute.inset-0.overflow-auto.flex.p-2").first
    geometry = canvas.evaluate(
        """
        (element) => {
          const margin = 16;
          const legend = element.querySelector('#legend');
          const diagramParts = [
            element.querySelector('#Axis'),
            element.querySelector('#label_leaders'),
            element.querySelector('#label_text')
          ].filter(Boolean);
          if (!legend || diagramParts.length !== 3) return null;

          const canvasRect = element.getBoundingClientRect();
          const legendRect = legend.getBoundingClientRect();
          const partRects = diagramParts.map((part) => part.getBoundingClientRect());
          const diagramLeft = Math.min(...partRects.map((rect) => rect.left));
          const diagramRight = Math.max(...partRects.map((rect) => rect.right));
          const safeLeft = canvasRect.left + margin;
          const safeRight = canvasRect.right - margin;
          const minimumDelta = Math.max(0, legendRect.right - safeRight);
          const maximumDelta = Math.max(0, diagramLeft - safeLeft);
          if (minimumDelta > maximumDelta + 1) {
            return {
              error: 'record and legend do not fit',
              minimumDelta,
              maximumDelta
            };
          }
          const delta = (minimumDelta + maximumDelta) / 2;
          const maxScroll = Math.max(0, element.scrollWidth - element.clientWidth);
          element.scrollLeft = Math.min(maxScroll, element.scrollLeft + delta);

          const finalCanvas = element.getBoundingClientRect();
          const finalLegend = legend.getBoundingClientRect();
          const finalParts = diagramParts.map((part) => part.getBoundingClientRect());
          return {
            clientWidth: element.clientWidth,
            scrollWidth: element.scrollWidth,
            scrollLeft: element.scrollLeft,
            canvasLeft: finalCanvas.left,
            canvasRight: finalCanvas.right,
            diagramLeft: Math.min(...finalParts.map((rect) => rect.left)),
            diagramRight: Math.max(...finalParts.map((rect) => rect.right)),
            legendLeft: finalLegend.left,
            legendRight: finalLegend.right,
            margin
          };
        }
        """
    )
    if not geometry or geometry.get("error"):
        raise AssertionError(f"Could not frame the restored legend: {geometry!r}")
    margin = float(geometry["margin"])
    if (
        geometry["diagramLeft"] < geometry["canvasLeft"] + margin - 1
        or geometry["diagramRight"] > geometry["canvasRight"] - margin + 1
        or geometry["legendLeft"] < geometry["canvasLeft"] + margin - 1
        or geometry["legendRight"] > geometry["canvasRight"] - margin + 1
    ):
        raise AssertionError(
            f"The restored record or legend remains clipped: {geometry!r}"
        )
    return {key: float(value) for key, value in geometry.items()}


def _download_svg(
    page: Page,
    *,
    download_dir: Path,
    expected_name: str,
) -> Path:
    with page.expect_download(timeout=ACTION_TIMEOUT_MS) as download_info:
        page.get_by_role("button", name="SVG", exact=True).click()
    download = download_info.value
    if download.failure() is not None:
        raise AssertionError(f"SVG download failed: {download.failure()}")
    if download.suggested_filename != expected_name:
        raise AssertionError(
            f"Expected {expected_name}, downloaded {download.suggested_filename}"
        )
    path = download_dir / expected_name
    download.save_as(path)
    if path.stat().st_size < 10_000:
        raise AssertionError(f"Downloaded SVG is unexpectedly small: {path}")
    return path


def _inspect_downloaded_edit(
    path: Path, feature_ids: tuple[str, ...]
) -> dict[str, Any]:
    report = inspect_svg_file(path)
    _assert_edited_human_svg(report)
    root = ET.parse(path).getroot()
    elements = list(root.iter())
    exported_styles: dict[str, list[dict[str, Any]]] = {}
    for feature_id in feature_ids:
        styles = [
            {
                "fill": element.attrib.get("fill", "").lower(),
                "stroke": element.attrib.get("stroke", "").lower(),
                "strokeWidth": float(
                    element.attrib.get("stroke-width", "0") or 0
                ),
            }
            for element in elements
            if element.attrib.get("data-gbdraw-feature-id", "").split(
                "__part", 1
            )[0]
            == feature_id
        ]
        _assert_target_edit({"styles": styles})
        exported_styles[feature_id] = styles
    report["editedFeatureStyles"] = exported_styles
    if EDIT_LEGEND not in set(report.get("texts", [])):
        raise AssertionError("The edited legend caption is absent from exported SVG")
    return report


def _capture_bgc_popups(
    browser_type: BrowserType,
    base_url: str,
    output_paths: Mapping[str, Path],
) -> tuple[dict[str, int], dict[str, Any], dict[str, Any]]:
    """Generate real LOSATP groups and inspect match and group popups."""

    assert_bgc_fixtures()
    capture = open_browser_capture(browser_type, base_url)
    page = capture.page
    screenshot_bytes: dict[str, int] = {}
    try:
        page.goto(base_url, wait_until="domcontentloaded")
        wait_for_worker(page)
        _stabilize_static_capture_surface(page)
        _set_bgc_inputs(page)
        _set_gallery_quality_presentation(
            page,
            title="LOSATP Similarity groups across five whole BGC records",
            match_gallery_definitions=False,
        )
        _configure_losatp(
            page,
            mode="orthogroup",
            output_prefix="inspect_similarity_groups",
        )
        generate_and_inspect(
            page,
            inspect_gui_bgc_losatp_svg,
            assert_gui_bgc_similarity_groups_svg,
        )
        fit_complete_linear_preview(page, target_zoom="40%")
        page.wait_for_timeout(350)

        match = page.get_by_role(
            "region", name="Result Preview", exact=True
        ).locator(
            '[data-match-kind="orthogroup"][data-orthogroup-id="og_1"]'
        ).first
        expect(match).to_be_visible()
        expect(match).to_have_attribute("role", "button")
        match.click()
        match_popup = page.get_by_role(
            "dialog", name="Pairwise match details", exact=True
        )
        expect(match_popup).to_be_visible()
        page.wait_for_timeout(250)
        match_text = match_popup.inner_text()
        for token in ("Similarity group ID", "Members", "og_"):
            if token.casefold() not in match_text.casefold():
                raise AssertionError(
                    f"Similarity-group match popup is missing {token!r}: "
                    f"{match_text!r}"
                )
        screenshot_bytes["match-popup.png"] = capture_screenshot(
            page, output_paths["match-popup.png"], "Linear"
        )
        page.get_by_role(
            "button", name="Close match popup", exact=True
        ).click()
        group_popup = _search_and_open_feature(
            page,
            query=r"^og_1$",
            field="orthogroup",
            use_regex=True,
            expected_matches=5,
        )
        page.wait_for_timeout(250)
        group_text = group_popup.inner_text()
        for token in ("Similarity group", "og_1", "Members", "Record coverage"):
            if token.casefold() not in group_text.casefold():
                raise AssertionError(
                    f"Similarity-group feature popup is missing {token!r}: "
                    f"{group_text!r}"
                )
        screenshot_bytes["group-popup.png"] = capture_screenshot(
            page, output_paths["group-popup.png"], "Linear"
        )
        page.get_by_title(
            "Open similarity-group editor", exact=True
        ).click()
        orthogroup_editor = page.locator(".right-drawer")
        expect(orthogroup_editor).to_be_visible()
        for token in ("23 similarity groups", "og_1", "5 members", "5 records"):
            expect(orthogroup_editor).to_contain_text(token)
        capture.assert_clean()
    finally:
        capture.close()
    return (
        screenshot_bytes,
        {"text": match_text, "mode": "orthogroup"},
        {"text": group_text, "groupId": "og_1", "members": 5},
    )


def capture_gui_interactive_editing(
    browser_type: BrowserType,
    base_url: str,
    output_paths: Mapping[str, Path],
    download_dir: Path,
) -> InteractiveEditingResult:
    """Run H-GUI-13 from raw HmmtDNA and five raw BGC records."""

    source_record = assert_human_mitochondrion_fixture()
    assert_output_paths(
        output_paths,
        GUI_INTERACTIVE_EDITING_SCREENSHOT_NAMES,
        "H-GUI-13",
    )
    download_dir.mkdir(parents=True, exist_ok=True)
    screenshot_bytes: dict[str, int] = {}

    capture = open_browser_capture(browser_type, base_url)
    page = capture.page
    try:
        page.goto(base_url, wait_until="domcontentloaded")
        wait_for_worker(page)
        _stabilize_static_capture_surface(page)
        load_raw_human_circular(page, output_prefix=EDITED_OUTPUT_PREFIX)
        generate_finished_human_diagram(page)
        fit_finished_human_preview(page)
        page.wait_for_timeout(350)
        popup = _search_and_open_feature(
            page,
            query=EDIT_TARGET_QUERY,
            field="qualifier-value",
            qualifier_key="gene",
            expected_matches=3,
        )
        if "COX1" not in (
            popup.inner_text()
            + _feature_search_panel(page).inner_text()
        ):
            raise AssertionError("COX1 search did not open the matching feature")
        screenshot_bytes["search-popup.png"] = capture_screenshot(
            page, output_paths["search-popup.png"], "Circular"
        )
        page.get_by_role(
            "button", name="Close feature popup", exact=True
        ).click()
        feature_ids = _select_all_search_matches(page, 3)
        if len(set(feature_ids)) != 3:
            raise AssertionError(f"COX search did not select three features: {feature_ids}")
        _edit_selected_features(page, feature_ids)
        _reset_finished_preview_viewport(page, target_zoom=60)
        expect(
            page.get_by_role(
                "status", name="Selected feature count", exact=True
            )
        ).to_have_text("3 selected")
        screenshot_bytes["editor.png"] = capture_screenshot(
            page, output_paths["editor.png"], "Circular"
        )
        selected_status = page.get_by_role(
            "status", name="Selected feature count", exact=True
        )
        page.get_by_role(
            "button", name="Apply visibility to selected features", exact=True
        ).click()
        expect(selected_status).to_have_count(0)
        page.get_by_role("button", name="Clear search", exact=True).click()
        close_toggle = page.get_by_title("Close editor", exact=True)
        if close_toggle.is_visible():
            close_toggle.click()
        expect(page.get_by_title("Open editor", exact=True)).to_be_visible()
        _assert_closed_right_drawer_contract(page)

        generate_and_inspect(
            page,
            inspect_first_circular_svg,
            lambda _report: None,
        )
        result_region = page.get_by_role(
            "region", name="Result Preview", exact=True
        )
        final_report = inspect_first_circular_svg(result_region)
        _assert_edited_human_svg(final_report)
        edited_feature_styles = {}
        for feature_id in feature_ids:
            target_style = _inspect_target_style(result_region, feature_id)
            _assert_target_edit(target_style)
            edited_feature_styles[feature_id] = target_style["styles"]
        if EDIT_LEGEND not in set(final_report.get("texts", [])):
            raise AssertionError("The edited legend caption did not survive regeneration")

        fit_finished_human_preview(page)
        final_report["previewScroll"] = _scroll_finished_preview_to_right(page)
        _, moved_legend_transform = _drag_legend(page)
        legend_transform = _inspect_target_style(
            result_region, feature_ids[0]
        ).get("legendTransform")
        if legend_transform != moved_legend_transform:
            raise AssertionError("The legend drag did not reach the committed SVG")
        final_report = inspect_first_circular_svg(result_region)
        _assert_edited_human_svg(final_report)
        final_report["editedFeatureStyles"] = edited_feature_styles

        expect(page.get_by_title("Open editor", exact=True)).to_be_visible()
        _assert_closed_right_drawer_contract(page)
        set_feature_search_visible(page, visible=False)
        final_report["closedDrawerGeometry"] = (
            _assert_closed_drawer_does_not_reserve_canvas(page)
        )
        screenshot_bytes["edited-result.png"] = capture_screenshot(
            page, output_paths["edited-result.png"], "Circular"
        )
        svg_path = _download_svg(
            page,
            download_dir=download_dir,
            expected_name=EDITED_SVG_NAME,
        )
        exported = _inspect_downloaded_edit(svg_path, feature_ids)
        if exported.get("legendTranslate") != final_report.get("legendTranslate"):
            raise AssertionError("The moved legend changed during SVG export")
        download_report = {
            "filename": svg_path.name,
            "bytes": svg_path.stat().st_size,
            "semantics": exported,
        }
        capture.assert_clean()
    finally:
        capture.close()

    popup_screenshots, match_popup, group_popup = _capture_bgc_popups(
        browser_type,
        base_url,
        output_paths,
    )
    screenshot_bytes.update(popup_screenshots)
    return InteractiveEditingResult(
        screenshot_bytes=screenshot_bytes,
        final_svg_semantics=final_report,
        download=download_report,
        source_record=source_record,
        match_popup=match_popup,
        group_popup=group_popup,
    )


def _semantic_signature(report: dict[str, Any]) -> dict[str, Any]:
    return {
        "recordIds": sorted(set(report.get("recordIds", []))),
        "featureIds": sorted(set(report.get("featureIds", []))),
        "texts": sorted(set(report.get("texts", []))),
        "labelTexts": sorted(set(report.get("labelTexts", []))),
        "groupChildCounts": report.get("groupChildCounts", {}),
        "recordTranslate": report.get("recordTranslate"),
        "legendTranslate": report.get("legendTranslate"),
    }


def _assert_semantically_equivalent(
    expected_report: dict[str, Any], actual_report: dict[str, Any]
) -> None:
    expected = _semantic_signature(expected_report)
    actual = _semantic_signature(actual_report)
    for key in (
        "recordIds",
        "featureIds",
        "texts",
        "labelTexts",
        "groupChildCounts",
    ):
        if actual.get(key) != expected.get(key):
            raise AssertionError(
                "Fresh-context session reload changed the SVG semantics: "
                f"{key}={actual.get(key)!r} != {expected.get(key)!r}"
            )
    for key in ("recordTranslate", "legendTranslate"):
        expected_position = expected.get(key)
        actual_position = actual.get(key)
        if (
            not expected_position
            or not actual_position
            or any(
                abs(float(left) - float(right)) >= 1
                for left, right in zip(
                    expected_position, actual_position, strict=True
                )
            )
        ):
            raise AssertionError(
                "Fresh-context session reload moved the SVG materially: "
                f"{key}={actual_position!r} != {expected_position!r}"
            )


def _validate_current_session(
    path: Path,
    *,
    expected_title: str = SESSION_TITLE,
    expected_output_prefix: str | None = None,
) -> dict[str, Any]:
    contents = path.read_bytes()
    if contents[:2] != b"\x1f\x8b":
        raise AssertionError("Saved session is not gzip-compressed")
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        payload = json.load(handle)
    if payload.get("format") != "gbdraw-session":
        raise AssertionError("Saved session format marker is missing")
    if payload.get("version") != CURRENT_SESSION_VERSION:
        raise AssertionError(
            f"Expected session version {CURRENT_SESSION_VERSION}, "
            f"found {payload.get('version')}"
        )
    render_request = payload.get("renderRequest", {})
    if (
        render_request.get("schema") != CURRENT_RENDER_REQUEST_SCHEMA
        or render_request.get("mode") != "circular"
        or render_request.get("output", {}).get("prefix")
        != (expected_output_prefix or expected_title)
    ):
        raise AssertionError("Saved session does not contain the current Circular request")
    if payload.get("title") != expected_title:
        raise AssertionError("Saved session title does not match its handoff filename")
    if not payload.get("results") or not payload.get("editorState"):
        raise AssertionError("Saved session is missing the committed result or editor state")

    matching_resources = []
    for resource in payload.get("resources", {}).values():
        if resource.get("kind") != "genbank" or resource.get("encoding") != "base64":
            continue
        decoded = base64.b64decode(resource.get("data", ""), validate=True)
        if len(decoded) != FIRST_CIRCULAR_FIXTURE_SIZE:
            continue
        if hashlib.sha256(decoded).hexdigest() == FIRST_CIRCULAR_FIXTURE_SHA256:
            matching_resources.append(resource)
    if len(matching_resources) != 1:
        raise AssertionError("Current session does not embed the exact HmmtDNA input")
    return {
        "filename": path.name,
        "bytes": len(contents),
        "version": payload["version"],
        "renderRequestSchema": render_request["schema"],
        "resources": len(payload.get("resources", {})),
    }


def _save_current_session(
    page: Page,
    download_dir: Path,
    *,
    title: str = SESSION_TITLE,
) -> Path:
    page.once("dialog", lambda dialog: dialog.accept(title))
    with page.expect_download(timeout=ACTION_TIMEOUT_MS) as download_info:
        with page.expect_event("dialog", timeout=ACTION_TIMEOUT_MS) as dialog_info:
            page.get_by_role(
                "button", name="Save Session", exact=True
            ).click()
        dialog = dialog_info.value
        if dialog.type != "prompt" or dialog.message != "Session title":
            raise AssertionError(f"Unexpected session-save dialog: {dialog.message}")
    download = download_info.value
    if download.failure() is not None:
        raise AssertionError(f"Session download failed: {download.failure()}")
    expected_filename = f"{title}.gbdraw-session.json.gz"
    if download.suggested_filename != expected_filename:
        raise AssertionError(
            f"Expected {expected_filename}, downloaded {download.suggested_filename}"
        )
    path = download_dir / expected_filename
    download.save_as(path)
    return path


def _load_current_session(page: Page, path: Path) -> None:
    page.once("dialog", lambda dialog: dialog.accept())
    with page.expect_event("dialog", timeout=ACTION_TIMEOUT_MS) as dialog_info:
        with page.expect_file_chooser(timeout=ACTION_TIMEOUT_MS) as chooser_info:
            page.get_by_role(
                "button", name="Load Session", exact=True
            ).click()
        chooser_info.value.set_files(path)
    dialog = dialog_info.value
    if dialog.message != "Session loaded successfully!":
        raise AssertionError(f"Session reload failed: {dialog.message}")


def _set_history_marker(page: Page, marker: str) -> Any:
    panel = page.get_by_label("Title & Legend", exact=True)
    panel.click()
    title = page.get_by_label("Plot Title", exact=True)
    title.fill(marker)
    title.press("Tab")
    expect(title).to_have_value(marker)
    panel.click()
    return title


def capture_gui_session_reproduction(
    browser_type: BrowserType,
    base_url: str,
    output_paths: Mapping[str, Path],
    download_dir: Path,
) -> SessionReproductionResult:
    """Run H-GUI-14 and reload its real v40 session in a fresh context."""

    source_record = assert_human_mitochondrion_fixture()
    assert_output_paths(
        output_paths,
        GUI_SESSION_REPRODUCTION_SCREENSHOT_NAMES,
        "H-GUI-14",
    )
    download_dir.mkdir(parents=True, exist_ok=True)
    screenshot_bytes: dict[str, int] = {}

    first_capture = open_browser_capture(browser_type, base_url)
    page = first_capture.page
    try:
        page.goto(base_url, wait_until="domcontentloaded")
        wait_for_worker(page)
        _stabilize_static_capture_surface(page)
        load_raw_human_circular(page, output_prefix=SESSION_TITLE)
        generate_finished_human_diagram(page)

        marker = "Undo and redo verification"
        title = _set_history_marker(page, marker)
        undo = page.get_by_role("button", name="Undo", exact=True)
        redo = page.get_by_role("button", name="Redo", exact=True)
        expect(undo).to_be_enabled()
        undo.click()
        expect(title).to_have_value("")
        expect(redo).to_be_enabled()
        redo.click()
        expect(title).to_have_value(marker)
        screenshot_bytes["history-actions.png"] = capture_screenshot(
            page, output_paths["history-actions.png"], "Circular"
        )

        page.once("dialog", lambda dialog: dialog.accept())
        with page.expect_event("dialog", timeout=ACTION_TIMEOUT_MS) as dialog_info:
            page.get_by_role(
                "button", name="Reset Settings", exact=True
            ).click()
        reset_dialog = dialog_info.value
        if reset_dialog.type != "confirm" or not reset_dialog.message.startswith(
            "Reset all settings"
        ):
            raise AssertionError(f"Unexpected reset dialog: {reset_dialog.message}")
        expect(page.get_by_label("Output Prefix", exact=True)).to_have_value("")
        expect(page.get_by_label("Species", exact=True)).to_have_value("")
        expect(page.get_by_label("Track Preset", exact=True)).to_have_value("tuckin")
        expect(
            page.get_by_role(
                "group", name="GenBank/DDBJ File selection", exact=True
            )
        ).to_contain_text("HmmtDNA.gbk")
        expect(
            page.get_by_role("region", name="Result Preview", exact=True)
        ).to_be_visible()

        apply_finished_human_settings(page, output_prefix=SESSION_TITLE)
        source_report = generate_finished_human_diagram(page)
        session_path = _save_current_session(page, download_dir)
        session_report = _validate_current_session(session_path)
        screenshot_bytes["session-download.png"] = capture_screenshot(
            page, output_paths["session-download.png"], "Circular"
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
        _assert_semantically_equivalent(source_report, restored_report)

        prefix = page.get_by_label("Output Prefix", exact=True)
        prefix.fill(RELOADED_OUTPUT_PREFIX)
        prefix.press("Tab")
        expect(prefix).to_have_value(RELOADED_OUTPUT_PREFIX)
        final_report = generate_finished_human_diagram(page)
        _assert_semantically_equivalent(source_report, final_report)
        _reset_finished_preview_viewport(page, target_zoom=60)
        final_report["restoredPreviewFrame"] = (
            _frame_finished_preview_with_legend(page)
        )
        set_feature_search_visible(page, visible=False)
        expect(
            page.get_by_role("button", name="Undo", exact=True)
        ).to_be_disabled()
        _stabilize_static_capture_surface(page)
        screenshot_bytes["reloaded-result.png"] = capture_screenshot(
            page, output_paths["reloaded-result.png"], "Circular"
        )
        svg_path = _download_svg(
            page,
            download_dir=download_dir,
            expected_name=RELOADED_SVG_NAME,
        )
        exported = inspect_svg_file(svg_path)
        assert_finished_circular_svg(exported)
        _assert_semantically_equivalent(source_report, exported)
        download_report = {
            "filename": svg_path.name,
            "bytes": svg_path.stat().st_size,
            "semantics": exported,
        }
        second_capture.assert_clean()
    finally:
        second_capture.close()

    return SessionReproductionResult(
        screenshot_bytes=screenshot_bytes,
        final_svg_semantics=final_report,
        download=download_report,
        session_download=session_report,
        source_record=source_record,
    )


__all__ = [
    "EDITED_SVG_NAME",
    "GUI_INTERACTIVE_EDITING_SCREENSHOT_NAMES",
    "GUI_SESSION_REPRODUCTION_SCREENSHOT_NAMES",
    "RELOADED_SVG_NAME",
    "SESSION_FILENAME",
    "capture_gui_interactive_editing",
    "capture_gui_session_reproduction",
]
