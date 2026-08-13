"""Capture and validate the H-GUI-15 publication export journey."""

from __future__ import annotations

import json
import math
import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping

from PIL import Image
from playwright.sync_api import BrowserType, expect

from assertions.svg_semantics import (
    assert_finished_circular_svg,
    inspect_svg_file,
)
from config import (
    GENERATION_TIMEOUT_MS,
    GUI_EXPORTS_SCREENSHOT_NAMES,
)
from flows.human_circular import (
    assert_human_mitochondrion_fixture,
    fit_finished_human_preview,
    generate_finished_human_diagram,
    load_raw_human_circular,
    pan_finished_human_preview_left,
)
from flows.web_capture import (
    assert_output_paths,
    capture_screenshot,
    open_browser_capture,
    set_feature_search_visible,
    wait_for_app_shell,
)


OUTPUT_PREFIX = "human_mitochondrion"
STATIC_SVG_NAME = f"{OUTPUT_PREFIX}.svg"
INTERACTIVE_SVG_NAME = f"{OUTPUT_PREFIX}.interactive.svg"
PNG_NAME = f"{OUTPUT_PREFIX}.png"
PDF_NAME = f"{OUTPUT_PREFIX}.pdf"
PNG_DPI = 300
PNG_MAGIC = b"\x89PNG\r\n\x1a\n"
INTERACTIVE_ASSET_IDS = {
    "gbdraw-interactive-feature-metadata",
    "gbdraw-interactive-feature-style",
    "gbdraw-interactive-feature-script",
}


@dataclass(frozen=True)
class ExportCaptureResult:
    """Every H-GUI-15 artifact and its parsed contract."""

    screenshot_bytes: dict[str, int]
    final_svg_semantics: dict[str, Any]
    download: dict[str, Any]
    downloads: dict[str, dict[str, Any]]
    source_record: dict[str, Any]


def _local_name(tag: str) -> str:
    return tag.rsplit("}", 1)[-1].lower()


def _svg_dimensions(path: Path) -> tuple[float, float]:
    root = ET.parse(path).getroot()

    def number(value: str | None) -> float | None:
        match = re.match(r"\s*([-+\d.eE]+)", value or "")
        return float(match.group(1)) if match else None

    width = number(root.attrib.get("width"))
    height = number(root.attrib.get("height"))
    if width and height:
        return width, height
    view_box = root.attrib.get("viewBox", "").replace(",", " ").split()
    if len(view_box) == 4:
        return float(view_box[2]), float(view_box[3])
    raise AssertionError(f"SVG has no usable dimensions: {path}")


def _download(
    page: Any,
    *,
    button_name: str,
    expected_name: str,
    download_dir: Path,
) -> Path:
    with page.expect_download(timeout=GENERATION_TIMEOUT_MS) as download_info:
        _export_button(page, button_name).click()
    download = download_info.value
    if download.failure() is not None:
        raise AssertionError(f"{button_name} export failed: {download.failure()}")
    if download.suggested_filename != expected_name:
        raise AssertionError(
            f"Expected {expected_name}, downloaded {download.suggested_filename}"
        )
    path = download_dir / expected_name
    download.save_as(path)
    if path.stat().st_size < 1_000:
        raise AssertionError(f"Exported file is unexpectedly small: {path}")
    return path


def _export_button(page: Any, name: str) -> Any:
    if name == "SVG":
        return page.get_by_role("button", name=name, exact=True)
    return page.get_by_role(
        "button", name=re.compile(rf"{re.escape(name)}$")
    )


def _validate_static_svg(path: Path) -> dict[str, Any]:
    if not path.read_bytes().lstrip().startswith(b"<"):
        raise AssertionError("Static SVG does not begin with XML markup")
    report = inspect_svg_file(path)
    assert_finished_circular_svg(report)
    width, height = _svg_dimensions(path)
    return {
        "filename": path.name,
        "bytes": path.stat().st_size,
        "width": width,
        "height": height,
        "semantics": report,
    }


def _validate_interactive_svg(
    path: Path,
    *,
    expected_dimensions: tuple[float, float],
) -> dict[str, Any]:
    report = inspect_svg_file(path)
    root = ET.parse(path).getroot()
    elements = list(root.iter())
    ids = {element.attrib.get("id", "") for element in elements}
    if root.attrib.get("data-gbdraw-interactive-svg") != "true":
        raise AssertionError("Interactive SVG marker is absent")
    if not INTERACTIVE_ASSET_IDS.issubset(ids):
        raise AssertionError(
            f"Interactive SVG assets are incomplete: {sorted(INTERACTIVE_ASSET_IDS - ids)}"
        )
    if report.get("scriptCount") != 1:
        raise AssertionError("Interactive SVG must contain one bundled runtime script")
    if report.get("eventAttributes") or report.get("unsafeHrefs"):
        raise AssertionError("Interactive SVG contains inline handlers or unsafe links")
    if "NC_012920.1" not in set(report.get("recordIds", [])):
        raise AssertionError("Interactive SVG lost the HmmtDNA record identity")
    if report.get("featureElementCount") != 37:
        raise AssertionError("Interactive SVG lost mitochondrial feature elements")
    metadata = next(
        element
        for element in elements
        if element.attrib.get("id") == "gbdraw-interactive-feature-metadata"
    )
    payload = json.loads(metadata.text or "{}")
    if payload.get("schema") != 3 or len(payload.get("items", [])) != 1:
        raise AssertionError("Interactive SVG does not embed one schema-3 catalog item")
    catalog_item = payload["items"][0]
    if len(catalog_item.get("features", [])) != 37:
        raise AssertionError("Interactive metadata does not contain all 37 features")
    if (
        root.attrib.get("width") != "100vw"
        or root.attrib.get("height") != "100vh"
        or root.attrib.get("preserveAspectRatio") != "xMidYMid meet"
    ):
        raise AssertionError("Interactive SVG does not use its responsive viewport")

    def original_dimension(attribute: str) -> float:
        match = re.match(r"\s*([-+\d.eE]+)", root.attrib.get(attribute, ""))
        if match is None:
            raise AssertionError(f"Interactive SVG has no usable {attribute}")
        return float(match.group(1))

    original_dimensions = (
        original_dimension("data-gbdraw-original-width"),
        original_dimension("data-gbdraw-original-height"),
    )
    if original_dimensions != expected_dimensions:
        raise AssertionError(
            "Interactive SVG did not preserve the static canvas geometry: "
            f"{original_dimensions} != {expected_dimensions}"
        )
    return {
        "filename": path.name,
        "bytes": path.stat().st_size,
        "width": original_dimensions[0],
        "height": original_dimensions[1],
        "viewportWidth": root.attrib["width"],
        "viewportHeight": root.attrib["height"],
        "schema": payload["schema"],
        "features": len(catalog_item["features"]),
        "semantics": report,
    }


def _validate_png(
    path: Path,
    *,
    svg_dimensions: tuple[float, float],
) -> dict[str, Any]:
    if path.read_bytes()[:8] != PNG_MAGIC:
        raise AssertionError("PNG signature is invalid")
    with Image.open(path) as image:
        image.verify()
    with Image.open(path) as image:
        dimensions = image.size
        mode = image.mode
    scale = PNG_DPI / 96
    expected = (int(svg_dimensions[0] * scale), int(svg_dimensions[1] * scale))
    if dimensions != expected:
        raise AssertionError(
            f"Expected {expected} pixels at {PNG_DPI} DPI, found {dimensions}"
        )
    return {
        "filename": path.name,
        "bytes": path.stat().st_size,
        "width": dimensions[0],
        "height": dimensions[1],
        "mode": mode,
        "dpi": PNG_DPI,
    }


def _validate_pdf(
    path: Path,
    *,
    svg_dimensions: tuple[float, float],
) -> dict[str, Any]:
    contents = path.read_bytes()
    if not contents.startswith(b"%PDF-") or b"%%EOF" not in contents[-1024:]:
        raise AssertionError("PDF signature or EOF marker is invalid")
    media_box = re.search(
        rb"/MediaBox\s*\[\s*([-+\d.]+)\s+([-+\d.]+)\s+([-+\d.]+)\s+([-+\d.]+)\s*\]",
        contents,
    )
    if media_box is None:
        raise AssertionError("PDF does not contain a readable MediaBox")
    values = tuple(float(value) for value in media_box.groups())
    page_dimensions = (values[2] - values[0], values[3] - values[1])
    if any(
        not math.isclose(actual, expected, abs_tol=1.0)
        for actual, expected in zip(
            sorted(page_dimensions), sorted(svg_dimensions), strict=True
        )
    ):
        raise AssertionError(
            f"PDF page {page_dimensions} does not match SVG {svg_dimensions}"
        )
    return {
        "filename": path.name,
        "bytes": len(contents),
        "width": page_dimensions[0],
        "height": page_dimensions[1],
    }


def _validate_interactive_behavior(browser: Any, path: Path) -> dict[str, Any]:
    context = browser.new_context()
    page = context.new_page()
    external_requests: list[str] = []
    page_errors: list[str] = []
    console_errors: list[str] = []
    page.on(
        "request",
        lambda request: external_requests.append(request.url)
        if request.url.startswith(("http://", "https://"))
        else None,
    )
    page.on("pageerror", lambda error: page_errors.append(str(error)))
    page.on(
        "console",
        lambda message: console_errors.append(message.text)
        if message.type == "error"
        else None,
    )
    try:
        page.goto(path.resolve().as_uri(), wait_until="load")
        root = page.locator('svg[data-gbdraw-interactive-svg="true"]')
        expect(root).to_be_attached()
        expect(page.locator("#gbdraw-feature-search-controls")).to_be_attached()
        feature = page.locator(
            '[data-gbdraw-interactive-feature="true"]'
        ).first
        expect(feature).to_be_attached()
        feature.click(force=True)
        popup = page.locator("#gbdraw-feature-popup")
        expect(popup).to_be_visible()
        if "Location" not in (popup.text_content() or ""):
            raise AssertionError("Interactive feature popup lacks feature metadata")
        page.locator("[data-close]").click()
        page.get_by_role(
            "button", name="Expand feature search", exact=True
        ).click()
        page.get_by_role(
            "searchbox", name="Search features", exact=True
        ).fill("COX1")
        page.get_by_role(
            "button", name="Search features", exact=True
        ).click()
        matches = page.locator(".gbdraw-interactive-feature--match")
        expect(matches).to_have_count(1)
        if external_requests or page_errors or console_errors:
            raise AssertionError(
                "Standalone interactive SVG was not clean and self-contained: "
                f"requests={external_requests!r}, pageErrors={page_errors!r}, "
                f"consoleErrors={console_errors!r}"
            )
        return {
            "popup": True,
            "searchMatches": matches.count(),
            "externalRequests": 0,
        }
    finally:
        context.close()


def capture_gui_exports(
    browser_type: BrowserType,
    base_url: str,
    output_paths: Mapping[str, Path],
    download_dir: Path,
) -> ExportCaptureResult:
    """Run H-GUI-15 from raw HmmtDNA through all visible export actions."""

    source_record = assert_human_mitochondrion_fixture()
    assert_output_paths(
        output_paths,
        GUI_EXPORTS_SCREENSHOT_NAMES,
        "H-GUI-15",
    )
    download_dir.mkdir(parents=True, exist_ok=True)
    screenshot_bytes: dict[str, int] = {}
    capture = open_browser_capture(browser_type, base_url)
    page = capture.page
    try:
        page.goto(base_url, wait_until="domcontentloaded")
        wait_for_app_shell(page)
        load_raw_human_circular(page, output_prefix=OUTPUT_PREFIX)
        final_report = generate_finished_human_diagram(page)
        page.wait_for_timeout(350)

        dpi = page.get_by_role("combobox", name="PNG DPI", exact=True)
        dpi.select_option(str(PNG_DPI))
        expect(dpi).to_have_value(str(PNG_DPI))
        screenshot_bytes["export-actions.png"] = capture_screenshot(
            page, output_paths["export-actions.png"], "Circular"
        )

        for button_name in ("SVG", "Interactive SVG", "PNG", "PDF"):
            expect(_export_button(page, button_name)).to_be_visible()
        interactive_path = _download(
            page,
            button_name="Interactive SVG",
            expected_name=INTERACTIVE_SVG_NAME,
            download_dir=download_dir,
        )
        static_path = _download(
            page,
            button_name="SVG",
            expected_name=STATIC_SVG_NAME,
            download_dir=download_dir,
        )
        png_path = _download(
            page,
            button_name="PNG",
            expected_name=PNG_NAME,
            download_dir=download_dir,
        )
        pdf_path = _download(
            page,
            button_name="PDF",
            expected_name=PDF_NAME,
            download_dir=download_dir,
        )

        static_report = _validate_static_svg(static_path)
        svg_dimensions = (static_report["width"], static_report["height"])
        interactive_report = _validate_interactive_svg(
            interactive_path,
            expected_dimensions=svg_dimensions,
        )
        interactive_report["behavior"] = _validate_interactive_behavior(
            capture.browser,
            interactive_path,
        )
        png_report = _validate_png(png_path, svg_dimensions=svg_dimensions)
        pdf_report = _validate_pdf(pdf_path, svg_dimensions=svg_dimensions)

        fit_finished_human_preview(page)
        pan_finished_human_preview_left(page)
        set_feature_search_visible(page, visible=False)
        page.wait_for_timeout(250)
        screenshot_bytes["exported-result.png"] = capture_screenshot(
            page, output_paths["exported-result.png"], "Circular"
        )
        capture.assert_clean()
    finally:
        capture.close()

    downloads = {
        "svg": static_report,
        "interactiveSvg": interactive_report,
        "png": png_report,
        "pdf": pdf_report,
    }
    return ExportCaptureResult(
        screenshot_bytes=screenshot_bytes,
        final_svg_semantics=final_report,
        download=static_report,
        downloads=downloads,
        source_record=source_record,
    )


__all__ = [
    "INTERACTIVE_SVG_NAME",
    "PDF_NAME",
    "PNG_NAME",
    "STATIC_SVG_NAME",
    "capture_gui_exports",
]
