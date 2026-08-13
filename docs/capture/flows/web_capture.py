"""Shared browser mechanics for documentation-owned web journeys."""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from importlib.metadata import version as distribution_version
from pathlib import Path
from typing import Any, Callable, Mapping
from urllib.parse import urlsplit

from PIL import Image
from playwright.sync_api import Browser, BrowserContext, BrowserType, Locator, Page, Route, expect

from config import (
    ACTION_TIMEOUT_MS,
    CHROMIUM_VERSION,
    COLOR_SCHEME,
    DEVICE_SCALE_FACTOR,
    GENERATION_TIMEOUT_MS,
    LOCALE,
    NAVIGATION_TIMEOUT_MS,
    PYTHON_PLAYWRIGHT_VERSION,
    REDUCED_MOTION,
    SCREENSHOT_MAX_BYTES,
    TIMEZONE_ID,
    VIEWPORT_HEIGHT,
    VIEWPORT_WIDTH,
    validate_capture_base_url,
)


@dataclass
class BrowserCapture:
    """One fresh, network-isolated browser context and its findings."""

    browser: Browser
    context: BrowserContext
    page: Page
    external_requests: list[str]
    page_errors: list[str]
    console_errors: list[str]

    def assert_clean(self) -> None:
        """Fail when the journey reached outside the app or logged an error."""

        if self.external_requests:
            raise AssertionError(
                "Capture attempted external requests: "
                f"{sorted(set(self.external_requests))}"
            )
        if self.page_errors:
            raise AssertionError(f"Capture page errors: {self.page_errors}")
        if self.console_errors:
            raise AssertionError(f"Capture console errors: {self.console_errors}")

    def close(self) -> None:
        self.context.close()
        self.browser.close()


def assert_fixture_identity(
    path: Path,
    *,
    expected_size: int,
    expected_sha256: str,
) -> None:
    """Verify that a capture uses the frozen packaged input bytes."""

    if not path.is_file():
        raise FileNotFoundError(f"Missing bundled tutorial fixture: {path}")
    contents = path.read_bytes()
    if len(contents) != expected_size:
        raise AssertionError(
            f"{path.name} size changed: expected {expected_size}, found {len(contents)}"
        )
    digest = hashlib.sha256(contents).hexdigest()
    if digest != expected_sha256:
        raise AssertionError(
            f"{path.name} checksum changed: expected {expected_sha256}, found {digest}"
        )


def assert_output_paths(
    output_paths: Mapping[str, Path],
    expected_names: tuple[str, ...],
    scenario_id: str,
) -> None:
    """Require every manifest-owned image and prepare its directory."""

    if tuple(output_paths) != expected_names:
        raise ValueError(
            f"{scenario_id} requires all {len(expected_names)} manifest-owned screenshot paths"
        )
    for path in output_paths.values():
        path.parent.mkdir(parents=True, exist_ok=True)


def _same_origin(url: str, base_url: str) -> bool:
    candidate = urlsplit(url)
    allowed = urlsplit(base_url)
    return (
        candidate.scheme == allowed.scheme
        and candidate.hostname == allowed.hostname
        and candidate.port == allowed.port
    )


def open_browser_capture(browser_type: BrowserType, base_url: str) -> BrowserCapture:
    """Open the pinned browser with a fresh context and same-origin routing."""

    validate_capture_base_url(base_url)
    installed_playwright = distribution_version("playwright")
    if installed_playwright != PYTHON_PLAYWRIGHT_VERSION:
        raise RuntimeError(
            "Authoritative capture requires Python Playwright "
            f"{PYTHON_PLAYWRIGHT_VERSION}; found {installed_playwright}"
        )

    browser = browser_type.launch(headless=True)
    if browser.version != CHROMIUM_VERSION:
        installed_chromium = browser.version
        browser.close()
        raise RuntimeError(
            f"Authoritative capture requires Chromium {CHROMIUM_VERSION}; "
            f"found {installed_chromium}"
        )

    context = browser.new_context(
        viewport={"width": VIEWPORT_WIDTH, "height": VIEWPORT_HEIGHT},
        device_scale_factor=DEVICE_SCALE_FACTOR,
        locale=LOCALE,
        timezone_id=TIMEZONE_ID,
        color_scheme=COLOR_SCHEME,
        reduced_motion=REDUCED_MOTION,
        service_workers="block",
        accept_downloads=True,
    )
    context.set_default_timeout(ACTION_TIMEOUT_MS)
    context.set_default_navigation_timeout(NAVIGATION_TIMEOUT_MS)

    external_requests: list[str] = []
    page_errors: list[str] = []
    console_errors: list[str] = []

    def enforce_local_requests(route: Route) -> None:
        request_url = route.request.url
        if _same_origin(request_url, base_url):
            route.continue_()
            return
        external_requests.append(request_url)
        route.abort()

    context.route("**/*", enforce_local_requests)
    page = context.new_page()
    page.on("pageerror", lambda error: page_errors.append(str(error)))
    page.on(
        "console",
        lambda message: console_errors.append(message.text)
        if message.type == "error"
        else None,
    )
    return BrowserCapture(
        browser=browser,
        context=context,
        page=page,
        external_requests=external_requests,
        page_errors=page_errors,
        console_errors=console_errors,
    )


def wait_for_app_shell(page: Page, *, timeout_ms: int = GENERATION_TIMEOUT_MS) -> None:
    """Wait for the mounted, palette-backed shell without starting Python."""

    page.wait_for_function(
        """
        () => {
          const app = window.__GBDRAW_APP__;
          return Boolean(
            app && Object.keys(app.paletteDefinitions || {}).length > 0
          );
        }
        """,
        timeout=timeout_ms,
    )
    shell_state = page.evaluate(
        """
        () => {
          const app = window.__GBDRAW_APP__;
          const runtimeFields = [
            'pyodide', 'pyodideReady', 'pyodideLoading', 'pyodideError', 'pyodideStatus'
          ];
          return {
            mounted: Boolean(app),
            paletteDefinitionCount: Object.keys(app?.paletteDefinitions || {}).length,
            mainLoaderPresent: typeof window.loadPyodide === 'function',
            mainRuntimeFields: app
              ? runtimeFields.filter((field) => Object.prototype.hasOwnProperty.call(app, field))
              : [],
            isolated: Boolean(window.crossOriginIsolated)
          };
        }
        """
    )
    if not shell_state["mounted"] or shell_state["paletteDefinitionCount"] == 0:
        raise AssertionError(f"Capture app shell did not become ready: {shell_state!r}")
    if shell_state["mainLoaderPresent"] or shell_state["mainRuntimeFields"]:
        raise AssertionError(
            "Capture app shell exposed a forbidden main-thread Python runtime: "
            f"{shell_state!r}"
        )
    if not shell_state["isolated"]:
        raise AssertionError("Capture page is not cross-origin isolated")


def generate_and_wait_for_result(
    page: Page,
    *,
    timeout_ms: int = GENERATION_TIMEOUT_MS,
    expected_status: str = "ok",
) -> dict[str, Any]:
    """Start Generate and await a committed Result or an explicit app error."""

    generate = page.get_by_role("button", name="Generate Diagram", exact=True)
    expect(generate).to_be_enabled()
    previous_run_marker = page.evaluate(
        "() => String(window.__GBDRAW_APP__?.lastRunInfo?.startedAtIso || '')"
    )
    generate.click()
    page.wait_for_function(
        """
        previous => {
          const app = window.__GBDRAW_APP__;
          if (!app || app.processing) return false;
          const marker = String(app.lastRunInfo?.startedAtIso || '');
          const committed = Boolean(
            marker && marker !== previous && Array.isArray(app.results) && app.results.length > 0
          );
          return committed || Boolean(app.errorLog);
        }
        """,
        arg=previous_run_marker,
        timeout=timeout_ms,
    )
    outcome = page.evaluate(
        """
        () => {
          const app = window.__GBDRAW_APP__;
          return {
            status: app?.errorLog ? 'error' : 'ok',
            errorLog: app?.errorLog || null,
            resultCount: Array.isArray(app?.results) ? app.results.length : 0,
            runMarker: String(app?.lastRunInfo?.startedAtIso || '')
          };
        }
        """
    )
    if outcome["status"] != expected_status:
        raise AssertionError(
            "Generate Diagram settled with an unexpected status: "
            f"expected {expected_status!r}, found {outcome!r}"
        )
    if expected_status == "ok" and outcome["resultCount"] < 1:
        raise AssertionError(f"Generate Diagram did not commit a Result: {outcome!r}")
    expect(generate).to_be_enabled(timeout=timeout_ms)
    return outcome


def linear_pair(page: Page, query_index: int, subject_index: int) -> Locator:
    endpoint_uids = []
    for index in (query_index, subject_index):
        record = page.get_by_role(
            "group", name=f"Linear sequence {index}", exact=True
        )
        uid = record.get_attribute("data-linear-record-uid")
        if not uid:
            raise AssertionError(f"Linear sequence {index} has no stable record UID")
        endpoint_uids.append(uid)
    pair = page.locator(
        f'fieldset[data-edge-key="{endpoint_uids[0]}->{endpoint_uids[1]}"]'
    )
    expect(pair).to_have_count(1)
    return pair


def open_linear_comparison_disclosure(
    page: Page,
    key: str,
    accessible_name: str,
) -> Locator:
    """Open one Linear comparison disclosure through its public summary."""

    details = page.locator(
        f'details[data-linear-comparison-disclosure="{key}"]'
    )
    expect(details).to_have_count(1)
    summary = details.get_by_role(
        "button", name=accessible_name, exact=True
    )
    expect(summary).to_be_visible()
    if details.get_attribute("open") is None:
        summary.click()
    expect(details).to_have_attribute("open", "")
    return details


def select_linear_losat_mode(
    settings: Locator, *, label: str, mode_key: str
) -> Locator:
    """Select one native LOSAT Mode button and verify its projected state."""

    group = settings.get_by_role("group", name="LOSAT Mode", exact=True)
    expect(group).to_be_visible()
    button = group.get_by_role("button", name=label, exact=True)
    expect(button).to_have_attribute(
        "data-linear-comparison-losat-mode-option", mode_key
    )
    button.click()
    expect(button).to_have_attribute("aria-pressed", "true")
    return button


def generate_and_inspect(
    page: Page,
    inspect_svg: Callable[[Any], dict[str, Any]],
    assert_svg: Callable[[dict[str, Any]], None],
    *,
    timeout_ms: int = GENERATION_TIMEOUT_MS,
) -> dict[str, Any]:
    """Click the public Generate action and inspect the resulting SVG."""

    generate_and_wait_for_result(page, timeout_ms=timeout_ms)

    result_region = page.get_by_role("region", name="Result Preview", exact=True)
    expect(result_region).to_be_visible(timeout=timeout_ms)
    report = inspect_svg(result_region)
    assert_svg(report)
    return report


def fit_complete_linear_preview(
    page: Page,
    target_zoom: str = "40%",
) -> None:
    """Use the public preview controls to center a whole Linear diagram in view."""

    reset_zoom = page.get_by_role("button", name="Reset zoom", exact=True)
    zoom_out = page.get_by_role("button", name="Zoom out", exact=True)
    if "100%" not in reset_zoom.inner_text():
        reset_zoom.click()
        expect(reset_zoom).to_contain_text("100%")
    for _ in range(10):
        if target_zoom in reset_zoom.inner_text():
            break
        zoom_out.click()
    else:
        raise AssertionError(
            f"Could not reach the documented Linear preview zoom: {target_zoom}"
        )

    page.wait_for_timeout(250)
    result_region = page.get_by_role("region", name="Result Preview", exact=True)
    canvas = result_region.locator(
        "div.absolute.inset-0.overflow-auto.flex.p-2"
    ).first

    def preview_geometry() -> dict[str, float] | None:
        return canvas.evaluate(
            """
            (element) => {
              const svg = element.querySelector(':scope > div > svg');
              const wrapper = svg?.parentElement;
              if (!wrapper) return null;
              const canvasRect = element.getBoundingClientRect();
              const wrapperRect = wrapper.getBoundingClientRect();
              return {
                canvasLeft: canvasRect.left,
                canvasRight: canvasRect.right,
                canvasWidth: canvasRect.width,
                canvasCenterX: (canvasRect.left + canvasRect.right) / 2,
                wrapperLeft: wrapperRect.left,
                wrapperRight: wrapperRect.right,
                wrapperWidth: wrapperRect.width,
                wrapperCenterX: (wrapperRect.left + wrapperRect.right) / 2,
              };
            }
            """
        )

    geometry = preview_geometry()
    if geometry is None:
        raise AssertionError("Could not resolve the rendered Linear preview")
    if geometry["wrapperWidth"] > geometry["canvasWidth"] + 4:
        raise AssertionError(
            "The documented Linear preview zoom does not fit the complete diagram: "
            f"{geometry!r}"
        )

    delta_x = geometry["canvasCenterX"] - geometry["wrapperCenterX"]
    if abs(delta_x) > 2:
        box = canvas.bounding_box()
        if box is None:
            raise AssertionError("Could not resolve the Result Preview bounds for panning")
        usable_width = box["width"] - 32
        if abs(delta_x) > usable_width:
            raise AssertionError(
                "The Linear preview requires an unsafe centering drag: "
                f"{geometry!r}"
            )
        start_x = box["x"] + (box["width"] - 16 if delta_x < 0 else 16)
        y = box["y"] + (box["height"] * 0.65)
        page.mouse.move(start_x, y)
        page.mouse.down()
        page.mouse.move(start_x + delta_x, y, steps=12)
        page.mouse.up()
        page.wait_for_timeout(250)

    centered = preview_geometry()
    if centered is None:
        raise AssertionError("The rendered Linear preview disappeared while centering")
    center_error = abs(centered["wrapperCenterX"] - centered["canvasCenterX"])
    edge_tolerance = max(
        2,
        ((centered["wrapperWidth"] - centered["canvasWidth"]) / 2) + 2,
    )
    if (
        center_error > 2
        or centered["wrapperLeft"] < centered["canvasLeft"] - edge_tolerance
        or centered["wrapperRight"] > centered["canvasRight"] + edge_tolerance
    ):
        raise AssertionError(
            "The complete Linear diagram is not centered inside Result Preview: "
            f"{centered!r}"
        )
    page.evaluate("() => window.getSelection()?.removeAllRanges()")


def set_feature_search_visible(page: Page, *, visible: bool) -> None:
    """Keep the floating search palette from covering a tutorial figure."""

    palette = page.query_selector(".preview-feature-search")
    if palette is None:
        raise AssertionError("Could not resolve the preview feature-search palette")
    palette.evaluate(
        """
        (element, shouldShow) => {
          if (shouldShow) element.style.removeProperty('display');
          else element.style.setProperty('display', 'none', 'important');
        }
        """,
        visible,
    )
    display = palette.evaluate("element => getComputedStyle(element).display")
    if (display != "none") != visible:
        raise AssertionError(
            f"Feature-search palette visibility did not change: {display!r}"
        )


def capture_screenshot(page: Page, path: Path, mode_name: str) -> int:
    """Write one pinned full-viewport PNG and enforce its size contract."""

    page.get_by_role("button", name=mode_name, exact=True).scroll_into_view_if_needed()
    page.wait_for_function(
        "() => !document.fonts || document.fonts.status === 'loaded'",
        timeout=ACTION_TIMEOUT_MS,
    )
    page.evaluate(
        """
        () => {
          for (const svg of document.querySelectorAll('.result-pane svg')) {
            if (svg.parentElement) svg.parentElement.style.transition = 'none';
          }
        }
        """
    )
    page.evaluate(
        """
        () => new Promise((resolve) => {
          requestAnimationFrame(() => requestAnimationFrame(resolve));
        })
        """
    )
    page.wait_for_timeout(350)
    page.screenshot(
        path=str(path),
        full_page=False,
        animations="disabled",
        caret="hide",
    )
    with Image.open(path) as screenshot:
        if screenshot.size != (VIEWPORT_WIDTH, VIEWPORT_HEIGHT):
            raise AssertionError(
                "Unexpected screenshot dimensions: "
                f"expected {(VIEWPORT_WIDTH, VIEWPORT_HEIGHT)}, found {screenshot.size}"
            )
        if screenshot.mode != "RGB":
            raise AssertionError(
                f"Unexpected screenshot mode: expected RGB, found {screenshot.mode}"
            )
    screenshot_bytes = path.stat().st_size
    if screenshot_bytes > SCREENSHOT_MAX_BYTES:
        raise AssertionError(
            f"Screenshot exceeds {SCREENSHOT_MAX_BYTES} bytes: {path} ({screenshot_bytes})"
        )
    return screenshot_bytes
