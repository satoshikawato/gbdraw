"""Capture real GUI figures into a run-local bundle."""

from __future__ import annotations

import base64
import json
import shutil
from dataclasses import asdict, is_dataclass
from pathlib import Path

from PIL import Image
from playwright.sync_api import Browser, sync_playwright

from config import screenshot_names_for, CHROMIUM_VERSION
from flows.tutorials.gui_first_circular import capture_first_circular
from flows.tutorials.gui_annotated_chloroplast import capture_gui_annotated_chloroplast
from flows.tutorials.gui_losatn import capture_gui_losatn
from flows.tutorials.gui_losatp_groups import capture_gui_losatp_groups
from video.model import sha256
from web_server import CaptureWebServer


INTRO_FLOWS = (
    ("human.circular", "T-GUI-01", capture_first_circular),
    ("tobacco.plastome", "T-GUI-05", capture_gui_annotated_chloroplast),
    ("lambda-de3.comparison", "T-GUI-03", capture_gui_losatn),
    ("bgc.comparison", "T-GUI-04", capture_gui_losatp_groups),
)


def rasterize_svg(browser: Browser, svg: Path, png: Path) -> tuple[int, int]:
    """Render a static SVG in pinned Chromium with a white caption-safe canvas."""

    source = svg.read_bytes()
    url = "data:image/svg+xml;base64," + base64.b64encode(source).decode("ascii")
    context = browser.new_context(viewport={"width": 1920, "height": 1080}, device_scale_factor=1)
    try:
        page = context.new_page()
        page.set_content(
            '<html><head><style>html,body{margin:0;background:white}'
            '#figure{position:absolute;left:80px;top:35px;width:1760px;height:880px;'
            'object-fit:contain}</style></head><body><img id="figure"></body></html>'
        )
        page.locator("#figure").evaluate("(img, url) => { img.src = url; }", url)
        page.wait_for_function("() => document.querySelector('#figure').complete && document.querySelector('#figure').naturalWidth > 0")
        page.evaluate("() => document.fonts.ready")
        png.parent.mkdir(parents=True, exist_ok=True)
        page.screenshot(path=str(png), animations="disabled")
    finally:
        context.close()
    with Image.open(png) as image:
        if image.size != (1920, 1080):
            raise AssertionError(f"Unexpected figure raster size: {image.size}")
        return image.size


def capture_intro_assets(run: Path, assets: dict, *, only: str | None = None) -> None:
    """Execute the four existing semantic GUI journeys in a fresh run."""

    with CaptureWebServer() as server, sync_playwright() as playwright:
        for asset_id, scenario_id, flow in INTRO_FLOWS:
            if only is not None and asset_id != only:
                continue
            evidence_dir = run / "evidence" / asset_id
            screenshot_dir = evidence_dir / "screenshots"
            download_dir = evidence_dir / "downloads"
            paths = {name: screenshot_dir / name for name in screenshot_names_for(scenario_id)}
            result = flow(playwright.chromium, server.base_url, paths, download_dir)
            downloaded = sorted(download_dir.glob("*.svg"))
            if len(downloaded) != 1:
                raise AssertionError(f"{asset_id} yielded {len(downloaded)} SVG downloads")
            svg = run / "figures" / f"{asset_id}.svg"
            svg.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(downloaded[0], svg)
            report_path = evidence_dir / "semantic.json"
            report_path.write_text(
                json.dumps(asdict(result) if is_dataclass(result) else result,
                           indent=2, default=str, ensure_ascii=False),
                encoding="utf-8",
            )
            browser = playwright.chromium.launch(headless=True)
            try:
                if browser.version != CHROMIUM_VERSION:
                    raise RuntimeError(f"Expected Chromium {CHROMIUM_VERSION}, got {browser.version}")
                png = run / "figures" / f"{asset_id}.png"
                width, height = rasterize_svg(browser, svg, png)
            finally:
                browser.close()
            assets[asset_id] = {
                "kind": "image", "path": str(png.relative_to(run)), "sha256": sha256(png),
                "width": width, "height": height,
                "svg": {"path": str(svg.relative_to(run)), "sha256": sha256(svg)},
                "evidence": {"path": str(report_path.relative_to(run)), "sha256": sha256(report_path)},
                "scenario": scenario_id,
            }
            (run / "assets.json").write_text(json.dumps({"schema_version": 1, "assets": assets}, indent=2), encoding="utf-8")
