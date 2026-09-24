"""Record a human-paced Web walkthrough: load, generate, edit, and export.

The journey drives the real local Web app with a visible, eased pointer, typed
input, and wheel scrolling. Chromium's screencast supplies high-resolution
frames only when the page repaints; the pointer, camera, captions, and time
compression are logged as events and composed later by ``walkthrough_render``.
"""

from __future__ import annotations

import base64
import json
import math
import re
import time
from pathlib import Path

from playwright.sync_api import Locator, Page, expect, sync_playwright

from PIL import Image, ImageChops

from config import ACTION_TIMEOUT_MS, REPO_ROOT
from flows.how_to.presentation import HUMAN_MITOCHONDRION_PATH
from flows.web_capture import (
    generate_and_wait_for_result,
    open_browser_capture,
    wait_for_app_shell,
)
from video.walkthrough_render import sha256
from web_server import CaptureWebServer


VIEWPORT = (1920, 1080)
SCALE = 2
FUNCTIONAL_RULES = (
    ("^ND", "#3b82f6", "NADH dehydrogenase"),
    ("^COX", "#ef4444", "Cytochrome c oxidase"),
    ("^ATP", "#f59e0b", "ATP synthase"),
    ("^CYTB$", "#8b5cf6", "Cytochrome b"),
)
# Wider than the default column so short form fields (for example the region
# label) show their whole value.
SETTINGS_RIGHT_EDGE = 450
GALLERY = (
    ("Vnig_TUMSAT-TG-2018", "Multi-replicon genomes"),
    ("hepatoplasmataceae_collinear", "Collinear genome comparisons"),
    ("BGC0000708-BGC0000713", "Biosynthetic gene clusters"),
    ("tobacco-chloroplast", "Plastomes with inverted repeats"),
    ("majanivirus_orthogroup", "Protein-similarity links"),
    ("HmmtDNA_ATskew", "GC and AT skew tracks"),
)
GALLERY_SOURCES = REPO_ROOT / "gbdraw/web/gallery/sources"
FINALE_SIZE = (1600, 900)
FINALE_PORTRAIT_SIZE = (936, 1060)
MENU_ROWS = 8
FINALE_FRAMES = 105
FINALE_ZOOM = 7.0

_TEXT_BOUNDS = """
(pattern) => {
  const region = document.querySelector('[role="region"][aria-label="Result Preview"]');
  const expression = new RegExp(pattern);
  const boxes = Array.from(region.querySelectorAll('svg text'))
    .filter((node) => expression.test(node.textContent || ''))
    .map((node) => node.getBoundingClientRect())
    .filter((box) => box.width > 0 && box.height > 0);
  if (!boxes.length) return null;
  const left = Math.min(...boxes.map((box) => box.left));
  const top = Math.min(...boxes.map((box) => box.top));
  const right = Math.max(...boxes.map((box) => box.right));
  const bottom = Math.max(...boxes.map((box) => box.bottom));
  return {x: left, y: top, width: right - left, height: bottom - top, count: boxes.length};
}
"""

_FEATURE_BOUNDS = """
() => {
  const region = document.querySelector('[role="region"][aria-label="Result Preview"]');
  const boxes = Array.from(region ? region.querySelectorAll('svg [data-gbdraw-feature-id]') : [])
    .map((node) => node.getBoundingClientRect())
    .filter((box) => box.width > 0 || box.height > 0);
  if (!boxes.length) return null;
  const left = Math.min(...boxes.map((box) => box.left));
  const top = Math.min(...boxes.map((box) => box.top));
  const right = Math.max(...boxes.map((box) => box.right));
  const bottom = Math.max(...boxes.map((box) => box.bottom));
  return {x: Math.round(left), y: Math.round(top), width: Math.round(right - left), height: Math.round(bottom - top)};
}
"""

_LARGEST_FEATURE = """
(fill) => {
  const region = document.querySelector('[role="region"][aria-label="Result Preview"]');
  const rows = Array.from(region.querySelectorAll('svg [data-gbdraw-feature-id]'))
    .filter((node) => (node.getAttribute('fill') || '').toLowerCase() === fill)
    .map((node) => ({node, box: node.getBoundingClientRect()}))
    .sort((a, b) => b.box.width * b.box.height - a.box.width * a.box.height);
  if (!rows.length) return null;
  const {box} = rows[0];
  return {x: box.left + box.width / 2, y: box.top + box.height / 2};
}
"""


def _ease(value: float) -> float:
    value = min(1.0, max(0.0, value))
    return value * value * (3 - 2 * value)


def _center(box: dict) -> tuple[float, float]:
    return box["x"] + box["width"] / 2, box["y"] + box["height"] / 2


class Recorder:
    """Drive one page like a person and log everything the editor needs."""

    def __init__(self, page: Page, raw: Path) -> None:
        self.page = page
        self.frames_dir = raw / "frames"
        self.frames_dir.mkdir(parents=True, exist_ok=True)
        self.frames: list[list] = []
        self.events: list[dict] = []
        self.x, self.y = 1240.0, 640.0
        self._cdp = None
        self._origin = 0.0

    # -- log ---------------------------------------------------------------
    def event(self, kind: str, **data) -> None:
        self.events.append({"t": round(time.time() - self._origin, 4), "kind": kind, **data})

    def start(self) -> None:
        cdp = self.page.context.new_cdp_session(self.page)

        def on_frame(params: dict) -> None:
            cdp.send("Page.screencastFrameAck", {"sessionId": params["sessionId"]})
            name = f"{len(self.frames):05d}.jpg"
            (self.frames_dir / name).write_bytes(base64.b64decode(params["data"]))
            stamp = params.get("metadata", {}).get("timestamp") or time.time()
            self.frames.append([round(stamp - self._origin, 4), name])

        cdp.on("Page.screencastFrame", on_frame)
        self._origin = time.time()
        cdp.send("Page.startScreencast", {
            "format": "jpeg", "quality": 92, "everyNthFrame": 1,
            "maxWidth": VIEWPORT[0] * SCALE, "maxHeight": VIEWPORT[1] * SCALE,
        })
        self._cdp = cdp
        self.page.mouse.move(self.x, self.y)
        self.event("move", x=self.x, y=self.y)
        # Force one paint so the timeline starts with a frame.
        self.page.evaluate("() => { document.body.style.outline = '0 solid transparent'; }")
        self.page.wait_for_timeout(400)
        self.page.evaluate("() => { document.body.style.removeProperty('outline'); }")
        self.page.wait_for_timeout(300)
        if not self.frames:
            raise AssertionError("Chromium screencast produced no frames")

    def stop(self) -> None:
        self.event("end")
        self.page.wait_for_timeout(300)
        self._cdp.send("Page.stopScreencast")

    # -- editorial cues ----------------------------------------------------
    def caption(self, step: int | None, title: str = "", subtitle: str = "") -> None:
        self.event("caption", step=step, title=title, subtitle=subtitle)

    def camera(self, rect: dict | None = None, zoom: float = 1.0, duration: float = 0.9,
               *, portrait: dict | None = None) -> None:
        """Frame ``rect`` (CSS px) at ``zoom``; no rect frames the whole page.

        The drawn diagram's bounds are logged as the subject, so tall outputs
        can frame the map where the wide output shows the whole page.
        ``portrait`` overrides the tall framing explicitly.
        """

        if rect is None:
            rect = {"x": 0, "y": 0, "width": VIEWPORT[0], "height": VIEWPORT[1]}
        cx, cy = _center(rect)
        subject = self.page.evaluate(_FEATURE_BOUNDS)
        self.event("camera", cx=round(cx, 1), cy=round(cy, 1), zoom=zoom, duration=duration,
                   rect={key: round(rect[key], 1) for key in ("x", "y", "width", "height")},
                   subject=subject, portrait=portrait)

    def focus(self, locator: Locator, zoom: float = 2.0, duration: float = 0.9) -> None:
        self.camera(locator.bounding_box(), zoom, duration)

    def toast(self, text: str, icon: str) -> None:
        self.event("toast", text=text, icon=icon)

    def mark(self, name: str) -> None:
        """Name a moment so edits can cut the recording by meaning, not by time."""

        self.event("mark", name=name)

    def fast_forward(self, *, factor: float | None = None, target: float | None = None) -> None:
        self.event("speed_start", factor=factor, target=target)

    def normal_speed(self) -> None:
        self.event("speed_end")

    def hold(self, seconds: float) -> None:
        self.page.wait_for_timeout(seconds * 1000)

    # -- pointer -----------------------------------------------------------
    def move_to(self, x: float, y: float, duration: float | None = None) -> None:
        distance = math.hypot(x - self.x, y - self.y)
        if distance < 1:
            return
        if duration is None:
            duration = min(1.05, max(0.38, 0.3 + distance / 1500))
        steps = max(10, int(duration * 50))
        sx, sy = self.x, self.y
        # A small perpendicular bow keeps the path from looking mechanical.
        bow = min(60.0, distance * 0.12)
        nx, ny = -(y - sy) / distance, (x - sx) / distance
        for index in range(1, steps + 1):
            p = _ease(index / steps)
            arc = 4 * p * (1 - p) * bow
            px, py = sx + (x - sx) * p + nx * arc, sy + (y - sy) * p + ny * arc
            self.page.mouse.move(px, py)
            self.event("move", x=round(px, 1), y=round(py, 1))
            self.page.wait_for_timeout(duration * 1000 / steps)
        self.x, self.y = x, y

    def _press(self, x: float, y: float, *, real: bool) -> None:
        self.move_to(x, y)
        self.page.wait_for_timeout(110)
        self.event("down", x=round(x, 1), y=round(y, 1))
        if real:
            self.page.mouse.down()
        self.page.wait_for_timeout(85)
        if real:
            self.page.mouse.up()
        self.event("up", x=round(x, 1), y=round(y, 1))

    def click(self, locator: Locator, *, pause: float = 0.3, left: float | None = None) -> None:
        """Click the center, or ``left`` CSS px from the left edge (for summaries)."""

        expect(locator).to_be_visible()
        box = locator.bounding_box()
        x, y = _center(box)
        if left is not None:
            x = box["x"] + left
        self._press(x, y, real=True)
        self.hold(pause)

    def click_point(self, x: float, y: float, *, pause: float = 0.3) -> None:
        self._press(x, y, real=True)
        self.hold(pause)

    def choose(self, select: Locator, value: str) -> None:
        """Select an option.

        Headless Chromium never paints native popups, so the option list and
        the chosen row are logged for the editor to draw.
        """

        box = select.bounding_box()
        options = select.evaluate(
            "(node, value) => ({labels: Array.from(node.options).map((o) => o.label.trim()),"
            " selected: Array.from(node.options).findIndex((o) => o.value === value)})", value)
        if options["selected"] < 0:
            raise AssertionError(f"{value!r} is not an option of the select")
        labels, selected = options["labels"], options["selected"]
        rows = min(MENU_ROWS, len(labels))
        first = min(max(0, selected - rows // 2), len(labels) - rows)
        x, y = _center(box)
        self._press(x, y, real=False)
        self.event("menu", x=box["x"], y=box["y"], width=box["width"], height=box["height"],
                   options=labels, selected=selected, window=[first, rows], duration=1.05)
        self.hold(0.2)
        # Point at the chosen row of the drawn menu (22 CSS px rows below the box).
        row_y = box["y"] + box["height"] + 6 + (selected - first + 0.5) * 22
        self._press(box["x"] + min(box["width"], 120) * 0.5, row_y, real=False)
        select.select_option(value)
        expect(select).to_have_value(value)
        self.hold(0.45)

    def pick_color(self, swatch: Locator, value: str) -> None:
        box = swatch.bounding_box()
        x, y = _center(box)
        self._press(x, y, real=False)
        self.event("swatch", x=box["x"], y=box["y"], width=box["width"], height=box["height"],
                   value=value, duration=0.8)
        self.hold(0.45)
        swatch.fill(value)
        self.hold(0.45)

    def type(self, field: Locator, text: str, *, delay: float = 0.075) -> None:
        self.click(field, pause=0.12)
        if not field.evaluate("node => node === document.activeElement"):
            raise AssertionError(f"Clicking {field} did not focus it")
        if field.input_value():
            self.page.keyboard.press("Control+A")
        for char in text:
            self.page.keyboard.type(char)
            self.page.wait_for_timeout(delay * 1000)
        expect(field).to_have_value(text)
        self.hold(0.25)

    def reveal(self, locator: Locator, *, anchor: float = 0.38) -> None:
        """Wheel-scroll the settings column until ``locator`` is comfortably visible."""

        scroller = self.page.locator(".settings-scroll")
        area = scroller.bounding_box()
        low, high = area["y"] + 70, area["y"] + area["height"] - 150
        box = locator.bounding_box()
        if low <= box["y"] and box["y"] + box["height"] <= high:
            return
        self.move_to(area["x"] + area["width"] * 0.62, area["y"] + area["height"] * 0.5)
        goal = area["y"] + area["height"] * anchor
        previous, stalled = None, 0
        for _ in range(240):
            box = locator.bounding_box()
            delta = box["y"] - goal
            position = scroller.evaluate("node => node.scrollTop")
            if abs(delta) < 45:
                break
            # Smooth scrolling can lag a wheel event, so only a run of
            # unchanged positions means the end of the column.
            stalled = stalled + 1 if position == previous else 0
            if stalled >= 6:
                break
            previous = position
            # Small wheel steps keep each screencast frame's jump short.
            step = max(-70.0, min(70.0, delta))
            self.page.mouse.wheel(0, step)
            self.page.wait_for_timeout(40)
        self.hold(0.35)
        box = locator.bounding_box()
        if not (low - 40 <= box["y"] and box["y"] + min(box["height"], 80) <= high + 40):
            raise AssertionError(f"Could not scroll {locator} into the settings view: {box}")

    def generate(self, *, target: float = 1.4) -> None:
        def press(button: Locator) -> None:
            self.click(button, pause=0.05)
            # Pull back so the progress overlay and then the result are in view.
            self.camera(duration=0.8)
            self.fast_forward(target=target)

        generate_and_wait_for_result(self.page, click=press)
        self.page.wait_for_function("() => document.fonts.status === 'loaded'")
        self.hold(0.35)
        self.normal_speed()

    # -- page queries ------------------------------------------------------
    def text_bounds(self, pattern: str) -> dict:
        bounds = self.page.evaluate(_TEXT_BOUNDS, pattern)
        if not bounds:
            raise AssertionError(f"No preview text matches {pattern!r}")
        return bounds


def _journey(rec: Recorder, downloads: Path) -> Path:
    page = rec.page
    preview = page.get_by_role("region", name="Result Preview", exact=True)

    # Tall outputs open on the settings column and the empty preview beside it.
    rec.camera(duration=0, portrait={"x": 0, "y": 0, "width": 820, "height": 930})
    rec.caption(1, "Open gbdraw.app", "The whole app runs in your browser")
    rec.hold(1.8)
    upload = page.get_by_role("button", name="Choose GenBank/DDBJ File", exact=True)
    rec.caption(1, "Load a GenBank file", "Human mitochondrial genome, NC_012920.1")
    rec.mark("upload")
    rec.focus(upload, zoom=2.1)
    rec.move_to(620, 380)
    rec.hold(0.5)
    with page.expect_file_chooser() as chooser:
        rec.click(upload, pause=0.1)
    chooser.value.set_files(HUMAN_MITOCHONDRION_PATH)
    rec.toast(HUMAN_MITOCHONDRION_PATH.name, "file")
    expect(page.get_by_role("group", name="GenBank/DDBJ File selection", exact=True)).to_contain_text(
        HUMAN_MITOCHONDRION_PATH.name)
    rec.hold(0.4)
    rec.mark("loaded")
    selection = page.get_by_role("group", name="GenBank/DDBJ File selection", exact=True).bounding_box()
    rec.camera({"x": selection["x"], "y": selection["y"], "width": selection["width"],
                "height": selection["height"] + 260}, zoom=2.1)
    rec.move_to(selection["x"] + 150, selection["y"] + selection["height"] + 70)
    rec.hold(2.4)

    rec.caption(2, "Generate the map", "The first run starts Python inside the browser")
    button = page.get_by_role("button", name="Generate Diagram", exact=True)
    rec.mark("generate")
    rec.focus(button, zoom=1.9)
    rec.hold(0.5)
    rec.generate(target=1.8)
    rec.mark("generated")
    rec.camera(duration=1.1)
    rec.move_to(1500, 700)
    rec.hold(1.2)
    rec.camera(preview.bounding_box(), zoom=1.12, duration=1.6)
    rec.hold(2.2)

    rec.caption(3, "Show feature labels", "Labels › Label Mode: Out")
    labels = page.get_by_label("Labels", exact=True)
    rec.camera(page.locator(".settings-scroll").bounding_box(), zoom=2.0)
    rec.reveal(labels)
    rec.focus(labels, zoom=2.0)
    rec.click(labels, left=70, pause=0.5)
    mode = page.locator("#circular-label-mode")
    rec.reveal(mode)
    rec.mark("labels")
    rec.focus(mode, zoom=2.1)
    rec.choose(mode, "out")
    rec.generate()
    rec.camera(duration=1.0)
    rec.hold(1.4)
    rec.caption(3, "Long product names crowd the map", "")
    right = rec.text_bounds(r"^NADH dehydrogenase subunit [12]$")
    rec.mark("crowded")
    rec.camera(right, zoom=1.8, duration=1.1)
    rec.move_to(1700, 480)
    rec.hold(2.4)

    rec.caption(3, "Prefer short gene symbols", "Qualifier priority: CDS › gene")
    order = page.locator('input[placeholder="product,gene,locus_tag"]')
    rec.camera(page.locator(".settings-scroll").bounding_box(), zoom=2.0)
    rec.reveal(order)
    rec.mark("priority")
    rec.focus(order, zoom=2.2)
    rec.choose(order.locator("xpath=preceding-sibling::select[1]"), "CDS")
    rec.mark("priority-type")
    rec.type(order, "gene")
    rec.click(order.locator("xpath=following-sibling::button[1]"), pause=0.6)
    rec.generate()
    rec.camera(duration=1.0)
    rec.hold(1.0)
    rec.mark("genes")
    rec.camera(right, zoom=1.8, duration=1.1)
    rec.hold(2.2)
    rec.camera(duration=1.0)
    rec.hold(1.0)

    rec.caption(4, "Color genes by function", "Colors › Specific rules match the gene qualifier")
    colors = page.get_by_label("Colors", exact=True)
    rec.camera(page.locator(".settings-scroll").bounding_box(), zoom=2.0)
    rec.reveal(colors)
    rec.focus(colors, zoom=2.0)
    rec.click(colors, left=70, pause=0.5)
    regex = page.locator('input[placeholder="Regex (e.g. hypothetical)"]')
    rule_box = regex.locator('xpath=ancestor::div[contains(@class,"bg-slate-50")][1]')
    for index, (pattern, color, legend) in enumerate(FUNCTIONAL_RULES):
        rec.reveal(rule_box, anchor=0.45)
        rec.mark(f"rule-{index + 1}")
        rec.focus(rule_box, zoom=2.2, duration=0.6)
        if index == 1:
            rec.fast_forward(factor=3.0)
        rec.choose(rule_box.locator("select"), "CDS")
        qualifier = rule_box.get_by_label("Qualifier name", exact=True)
        if qualifier.input_value() != "gene":
            rec.type(qualifier, "gene")
        rec.type(regex, pattern)
        rec.pick_color(rule_box.locator('input[type="color"]'), color)
        rec.type(rule_box.locator('input[placeholder="Legend Caption (Optional)"]'), legend)
        rec.click(rule_box.locator("button").last, pause=0.45)
    rec.normal_speed()
    rec.hold(0.5)
    rec.mark("colors")
    rec.generate()
    rec.camera(duration=1.0)
    rec.hold(1.2)
    rec.mark("legend")
    legend_box = rec.text_bounds(r"^(NADH dehydrogenase|Cytochrome c oxidase|ATP synthase|Cytochrome b)$")
    rec.camera(legend_box, zoom=1.9, duration=1.1)
    rec.hold(2.0)
    rec.camera(duration=1.0)
    rec.hold(0.8)

    rec.caption(5, "Mark a region, even across the origin", "Region Annotations › D-loop, 16,024–576")
    region = page.get_by_label("Region Annotations", exact=True)
    rec.camera(page.locator(".settings-scroll").bounding_box(), zoom=2.0)
    rec.reveal(region)
    rec.focus(region, zoom=2.0)
    rec.click(region, left=70, pause=0.4)
    rec.click(page.get_by_role("button", name="Add set"), pause=0.4)
    add_span = page.get_by_role("button", name="Coordinates")
    rec.reveal(add_span)
    rec.focus(add_span, zoom=2.2)
    rec.click(add_span, pause=0.4)
    start = page.locator('input[placeholder="Start (1-based)"]')
    rec.reveal(start)
    rec.mark("coordinates")
    rec.focus(start, zoom=2.3)
    rec.type(start, "16024")
    rec.type(page.locator('input[placeholder="End (inclusive)"]'), "576")
    rec.type(page.locator('input[placeholder="Label"]'), "D-loop")
    rec.mark("bracket")
    rec.choose(page.locator('select:has(option[value="bracket"])'), "bracket")
    rec.generate()
    rec.camera(duration=1.0)
    rec.hold(1.0)
    dloop = page.locator('[role="region"][aria-label="Result Preview"] svg [data-gbdraw-annotation-id]').first
    rec.mark("dloop")
    rec.focus(dloop, zoom=2.4, duration=1.2)
    rec.move_to(*_center(dloop.bounding_box()))
    rec.hold(2.6)
    rec.camera(duration=1.0)
    rec.hold(0.6)

    rec.caption(6, "Click any feature for details", "")
    point = page.evaluate(_LARGEST_FEATURE, "#ef4444")
    if point is None:
        raise AssertionError("No functionally colored cytochrome c oxidase feature was drawn")
    rec.mark("popup")
    rec.click_point(point["x"], point["y"], pause=0.4)
    popup = page.get_by_role("dialog", name=re.compile(r"^Feature details"))
    expect(popup).to_be_visible()
    rec.mark("popup-open")
    box = popup.bounding_box()
    rec.camera({"x": min(box["x"], point["x"]) - 40, "y": box["y"] - 40,
                "width": max(box["x"] + box["width"], point["x"]) - min(box["x"], point["x"]) + 80,
                "height": box["height"] + 80}, zoom=1.55, duration=1.0)
    rec.hold(3.0)
    rec.click(page.get_by_label("Close feature popup", exact=True), pause=0.3)
    rec.camera(duration=0.9)
    rec.hold(0.6)

    rec.caption(7, "Export a vector SVG", "PNG and PDF are one click away")
    svg_button = page.get_by_role("button", name="SVG", exact=True)
    rec.mark("export")
    rec.focus(svg_button, zoom=2.1)
    rec.hold(0.4)
    with page.expect_download(timeout=ACTION_TIMEOUT_MS) as info:
        rec.click(svg_button, pause=0.1)
    download = info.value
    if download.failure() is not None:
        raise AssertionError(f"SVG download failed: {download.failure()}")
    downloads.mkdir(parents=True, exist_ok=True)
    exported = downloads / download.suggested_filename
    download.save_as(exported)
    rec.toast(f"{exported.name} downloaded ({exported.stat().st_size // 1024} KB)", "download")
    rec.mark("downloaded")
    rec.hold(2.2)
    return exported


def _inspect_export(svg: Path) -> dict:
    text = svg.read_text(encoding="utf-8")
    labels = set(re.findall(r">([^<>]{1,60})</(?:text|textPath)>", text))
    genes = {"ND1", "ND2", "ND3", "ND4", "ND4L", "ND5", "ND6", "COX1", "COX2", "COX3", "ATP6", "ATP8", "CYTB"}
    fills = {fill.lower() for fill in re.findall(r'fill="(#[0-9a-fA-F]{6})"', text)}
    report = {
        "gene_symbols": sorted(genes & labels),
        "product_labels_present": any(label.startswith("NADH dehydrogenase subunit") for label in labels),
        "functional_fills": sorted({color for _, color, _ in FUNCTIONAL_RULES} & fills),
        "legend": sorted({legend for *_, legend in FUNCTIONAL_RULES} & labels),
        "dloop": "D-loop" in labels,
    }
    if (set(report["gene_symbols"]) != genes or report["product_labels_present"]
            or len(report["functional_fills"]) != 4 or len(report["legend"]) != 4 or not report["dloop"]):
        raise AssertionError(f"Exported SVG lacks the recorded edits: {report}")
    return report


def _render_finale(browser, svg: Path, out: Path, size: tuple[int, int]) -> dict:
    """Zoom into the exported file itself; every frame is a fresh vector raster."""

    out.mkdir(parents=True, exist_ok=True)
    context = browser.new_context(viewport={"width": size[0], "height": size[1]}, device_scale_factor=1)
    try:
        page = context.new_page()
        page.set_content('<html><body style="margin:0;background:#fff;overflow:hidden">'
                         f'{svg.read_text(encoding="utf-8")}</body></html>')
        page.evaluate("() => document.fonts.ready")
        geometry = page.evaluate(
            """
            () => {
              const svg = document.querySelector('svg');
              const vb = svg.viewBox.baseVal;
              const base = vb && vb.width ? [vb.x, vb.y, vb.width, vb.height]
                : [0, 0, svg.width.baseVal.value, svg.height.baseVal.value];
              svg.setAttribute('width', innerWidth);
              svg.setAttribute('height', innerHeight);
              svg.setAttribute('preserveAspectRatio', 'xMidYMid meet');
              svg.setAttribute('viewBox', base.join(' '));
              const rect = svg.getBoundingClientRect();
              const scale = Math.min(rect.width / base[2], rect.height / base[3]);
              const offX = (rect.width - base[2] * scale) / 2, offY = (rect.height - base[3] * scale) / 2;
              const target = Array.from(svg.querySelectorAll('text, textPath'))
                .find((node) => node.textContent === 'tRNA-Pro') || svg.querySelector('[data-gbdraw-annotation-id]');
              const box = target.getBoundingClientRect();
              return {base, cx: base[0] + (box.left + box.width / 2 - offX) / scale,
                      cy: base[1] + (box.top + box.height / 2 - offY) / scale};
            }
            """
        )
        x0, y0, w0, h0 = geometry["base"]
        cx, cy = geometry["cx"], geometry["cy"]
        for index in range(FINALE_FRAMES):
            p = _ease(min(1.0, index / (FINALE_FRAMES * 0.8)))
            zoom = math.exp(math.log(FINALE_ZOOM) * p)
            w, h = w0 / zoom, h0 / zoom
            fx = x0 + (cx - x0 - w0 / 2) * (1 - 1 / zoom) + w0 / 2 - w / 2
            fy = y0 + (cy - y0 - h0 / 2) * (1 - 1 / zoom) + h0 / 2 - h / 2
            page.evaluate("(box) => document.querySelector('svg').setAttribute('viewBox', box)",
                          f"{fx:.3f} {fy:.3f} {w:.3f} {h:.3f}")
            page.screenshot(path=str(out / f"{index:04d}.png"))
        return {"frames": FINALE_FRAMES, "zoom": FINALE_ZOOM, "size": list(size)}
    finally:
        context.close()


def _rasterize(browser, svg: Path, png: Path, size: tuple[int, int] = (1400, 1000), scale: float = 2) -> None:
    """Render a static SVG in Chromium, then trim the white margin."""

    context = browser.new_context(viewport={"width": size[0], "height": size[1]}, device_scale_factor=scale)
    try:
        page = context.new_page()
        page.set_content('<html><body style="margin:0;background:#fff">'
                         f'<img id="f" style="width:{size[0]}px;height:{size[1]}px;object-fit:contain"></body></html>')
        url = "data:image/svg+xml;base64," + base64.b64encode(svg.read_bytes()).decode("ascii")
        page.locator("#f").evaluate("(img, url) => { img.src = url; }", url)
        page.wait_for_function("() => document.querySelector('#f').complete && document.querySelector('#f').naturalWidth > 0")
        page.screenshot(path=str(png))
    finally:
        context.close()
    with Image.open(png) as image:
        rgb = image.convert("RGB")
    bbox = ImageChops.difference(rgb, Image.new("RGB", rgb.size, "white")).point(
        lambda value: 255 if value > 12 else 0).getbbox()
    if bbox is None:
        raise AssertionError(f"{svg.name} rendered blank")
    pad = round(24 * scale)
    rgb.crop((max(0, bbox[0] - pad), max(0, bbox[1] - pad),
              min(rgb.width, bbox[2] + pad), min(rgb.height, bbox[3] + pad))).save(png)


def record_walkthrough(run: Path) -> Path:
    """Record the journey into ``run/raw/walkthrough`` and return its manifest."""

    raw = run / "raw" / "walkthrough"
    raw.mkdir(parents=True, exist_ok=True)
    with CaptureWebServer() as server, sync_playwright() as playwright:
        # Screencast frames follow the physical window, so the scale must be
        # forced at launch as well as emulated for the page.
        capture = open_browser_capture(playwright.chromium, server.base_url,
                                       viewport=VIEWPORT, device_scale_factor=SCALE,
                                       launch_args=(f"--force-device-scale-factor={SCALE}",))
        page = capture.page
        try:
            page.goto(server.base_url, wait_until="domcontentloaded")
            wait_for_app_shell(page)
            # The floating feature search would cover the top right of every
            # result; it is the only element hidden for the recording.
            page.add_style_tag(content=".preview-feature-search{display:none !important}")
            circular = page.get_by_role("button", name="Circular", exact=True)
            circular.click()
            expect(circular).to_have_attribute("aria-pressed", "true")
            handle = page.get_by_title("Drag to resize", exact=True).bounding_box()
            page.mouse.move(*_center(handle))
            page.mouse.down()
            page.mouse.move(SETTINGS_RIGHT_EDGE, _center(handle)[1], steps=8)
            page.mouse.up()
            page.wait_for_function("() => document.fonts.status === 'loaded'")
            recorder = Recorder(page, raw)
            recorder.start()
            exported = _journey(recorder, raw / "downloads")
            recorder.stop()
            capture.assert_clean()
        finally:
            capture.close()
        report = _inspect_export(exported)
        browser = playwright.chromium.launch(headless=True)
        try:
            finale = {name: _render_finale(browser, exported, raw / name, size)
                      for name, size in (("finale", FINALE_SIZE), ("finale-portrait", FINALE_PORTRAIT_SIZE))}
            _rasterize(browser, exported, raw / "final-figure.png")
            gallery = []
            for stem, label in GALLERY:
                source = GALLERY_SOURCES / f"{stem}.svg"
                target = raw / "gallery" / f"{stem}.png"
                target.parent.mkdir(exist_ok=True)
                _rasterize(browser, source, target, size=(1600, 1000), scale=1.5)
                gallery.append({"path": str(target.relative_to(raw)), "label": label,
                                "source": str(source.relative_to(REPO_ROOT)), "source_sha256": sha256(source)})
        finally:
            browser.close()
    (raw / "frames.json").write_text(json.dumps(recorder.frames), encoding="utf-8")
    (raw / "events.json").write_text(json.dumps(recorder.events, indent=1), encoding="utf-8")
    manifest = {
        "schema_version": 1,
        "viewport": list(VIEWPORT), "scale": SCALE,
        "source": {"path": HUMAN_MITOCHONDRION_PATH.name, "sha256": sha256(HUMAN_MITOCHONDRION_PATH)},
        "export": {"path": str(exported.relative_to(raw)), "sha256": sha256(exported), "semantics": report},
        "frames": len(recorder.frames), "events": len(recorder.events),
        "duration": recorder.events[-1]["t"], "finale": finale, "gallery": gallery,
    }
    path = raw / "walkthrough.json"
    path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return path
