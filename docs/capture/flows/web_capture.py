"""Shared browser mechanics for documentation-owned web journeys."""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from importlib.metadata import version as distribution_version
from pathlib import Path
from typing import Any, Callable, Mapping
from urllib.parse import urlsplit

from PIL import Image
from playwright.sync_api import Browser, BrowserContext, BrowserType, Page, Route, expect

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
    WORKER_READY_TIMEOUT_MS,
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


def wait_for_worker(page: Page) -> None:
    """Wait for the packaged diagram worker and verify isolation."""

    page.wait_for_function(
        """
        () => {
          const app = window.__GBDRAW_APP__;
          return Boolean(
            app && (app.diagramGenerationWorkerReady || app.diagramGenerationWorkerError)
          );
        }
        """,
        timeout=WORKER_READY_TIMEOUT_MS,
    )
    worker_state = page.evaluate(
        """
        () => ({
          ready: Boolean(window.__GBDRAW_APP__?.diagramGenerationWorkerReady),
          error: String(window.__GBDRAW_APP__?.diagramGenerationWorkerError || ''),
          status: String(window.__GBDRAW_APP__?.diagramGenerationWorkerStatus || ''),
          isolated: Boolean(window.crossOriginIsolated)
        })
        """
    )
    if not worker_state["ready"]:
        raise AssertionError(
            "Diagram worker did not become ready: "
            f"{worker_state['error'] or worker_state['status']}"
        )
    if not worker_state["isolated"]:
        raise AssertionError("Capture page is not cross-origin isolated")


def generate_and_inspect(
    page: Page,
    inspect_svg: Callable[[Any], dict[str, Any]],
    assert_svg: Callable[[dict[str, Any]], None],
    *,
    timeout_ms: int = GENERATION_TIMEOUT_MS,
) -> dict[str, Any]:
    """Click the public Generate action and inspect the resulting SVG."""

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
          const marker = String(app?.lastRunInfo?.startedAtIso || '');
          return Boolean(app && !app.processing && marker && marker !== previous);
        }
        """,
        arg=previous_run_marker,
        timeout=timeout_ms,
    )
    expect(generate).to_be_enabled(timeout=timeout_ms)

    result_region = page.get_by_role("region", name="Result Preview", exact=True)
    expect(result_region).to_be_visible(timeout=timeout_ms)
    report = inspect_svg(result_region)
    assert_svg(report)
    return report


def fit_complete_linear_preview(
    page: Page,
    target_zoom: str = "40%",
    *,
    pan_left: bool = False,
) -> None:
    """Use the public preview controls to fit a whole Linear record in view."""

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

    if pan_left:
        result_region = page.get_by_role("region", name="Result Preview", exact=True)
        box = result_region.bounding_box()
        if box is None:
            raise AssertionError("Could not resolve the Result Preview bounds for panning")
        y = box["y"] + (box["height"] * 0.82)
        page.mouse.move(box["x"] + (box["width"] * 0.80), y)
        page.mouse.down()
        page.mouse.move(box["x"] + (box["width"] * 0.10), y, steps=12)
        page.mouse.up()
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
