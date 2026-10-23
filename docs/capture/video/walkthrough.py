"""Record a human-paced Web walkthrough: load, generate, edit, and export.

A journey drives the real local Web app with a visible, eased pointer, typed
input, and wheel scrolling. Chromium's screencast supplies high-resolution
frames only when the page repaints; the pointer, camera, captions, and time
compression are logged as events and composed later by ``walkthrough_render``.
``meet-gbdraw`` is the introduction; ``video.tutorials`` adds one journey per
Gallery tutorial.
"""

from __future__ import annotations

import base64
import json
import math
import re
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

from playwright.sync_api import Locator, Page, expect, sync_playwright

from PIL import Image, ImageChops

from config import ACTION_TIMEOUT_MS, GENERATION_TIMEOUT_MS, REPO_ROOT
from flows.how_to.presentation import HUMAN_MITOCHONDRION_PATH
from flows.web_capture import (
    generate_and_wait_for_result,
    open_browser_capture,
    wait_for_app_shell,
)
from video.walkthrough_render import PORTRAIT, PORTRAIT_PAD, STORIES, _contains, sha256
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
# A portrait take opens the exported file in a tab shaped like the tall video
# window, and fits the preview figure narrow enough for that window to hold.
PORTRAIT_TAB = (round(VIEWPORT[1] * PORTRAIT.aspect), VIEWPORT[1])
PORTRAIT_FIGURE_WIDTH = PORTRAIT_TAB[0] / PORTRAIT_PAD
# The Result Preview region, found by its structure: its aria-label is
# translated with the rest of the UI.
PREVIEW_REGION = '[role="region"]:has(.result-tabs)'

_TEXT_BOUNDS = """
(pattern) => {
  const region = document.querySelector('$REGION');
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

_PREVIEW = """
() => {
  const region = document.querySelector('$REGION');
  const canvas = region?.querySelector('div.absolute.inset-0.overflow-auto');
  const svg = canvas?.querySelector(':scope > div > svg');
  if (!svg) return null;
  const box = (node) => { const b = node.getBoundingClientRect(); return {x: b.left, y: b.top, width: b.width, height: b.height}; };
  // Spots a pan may start from. Spots beside the drawing come first: a press
  // there never reaches the drawing's click handlers, which look for a
  // feature under the release point even after a pan.
  const beside = [], blank = [];
  const c = canvas.getBoundingClientRect();
  for (const fx of [0.04, 0.96, 0.5, 0.25, 0.75]) {
    for (const fy of [0.5, 0.3, 0.7]) {
      const x = c.left + c.width * fx, y = c.top + c.height * fy;
      const hit = document.elementFromPoint(x, y);
      if (hit && canvas.contains(hit) && !hit.closest('[data-gbdraw-feature-id], [data-match-kind], text, button'))
        (svg.contains(hit) ? blank : beside).push([x, y]);
    }
  }
  blank.unshift(...beside);
  // The zoom and layout buttons float over the bottom of the canvas.
  const controls = region.querySelector('.preview-controls')?.getBoundingClientRect();
  const bottom = controls?.height ? Math.min(c.bottom, controls.top - 8) : c.bottom;
  // The preview's own zoom: the scale the app writes on the surface around the drawing.
  const zoom = new DOMMatrixReadOnly(svg.parentElement.style.transform || 'none').a;
  return {canvas: {x: c.left, y: c.top, width: c.width, height: bottom - c.top}, svg: box(svg), blank, zoom};
}
"""

_FEATURE_BOUNDS = """
() => {
  const region = document.querySelector('$REGION');
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
  const region = document.querySelector('$REGION');
  const rows = Array.from(region.querySelectorAll('svg [data-gbdraw-feature-id]'))
    .filter((node) => (node.getAttribute('fill') || '').toLowerCase() === fill)
    .map((node) => ({node, box: node.getBoundingClientRect()}))
    .sort((a, b) => b.box.width * b.box.height - a.box.width * a.box.height);
  if (!rows.length) return null;
  const {box} = rows[0];
  return {x: box.left + box.width / 2, y: box.top + box.height / 2};
}
"""


_TEXT_BOUNDS, _PREVIEW, _FEATURE_BOUNDS, _LARGEST_FEATURE = (
    script.replace("$REGION", PREVIEW_REGION) for script in (_TEXT_BOUNDS, _PREVIEW, _FEATURE_BOUNDS, _LARGEST_FEATURE))


_SET_VALUE = """
(node, value) => {
  node.value = value;
  node.dispatchEvent(new Event('change', {bubbles: true}));
  return node.value;
}
"""


def _ease(value: float) -> float:
    value = min(1.0, max(0.0, value))
    return value * value * (3 - 2 * value)


def _center(box: dict) -> tuple[float, float]:
    return box["x"] + box["width"] / 2, box["y"] + box["height"] / 2


class Recorder:
    """Drive one page like a person and log everything the editor needs."""

    def __init__(self, page: Page, raw: Path, *, portrait: bool = False) -> None:
        self.page = page
        self.app = page                     # the app's tab; ``page`` follows an opened file
        self.portrait = portrait            # a take for the tall video
        self.frames_dir = raw / "frames"
        self.frames_dir.mkdir(parents=True, exist_ok=True)
        self.frames: list[list] = []
        self.events: list[dict] = []
        # Events and the pointer use stage coordinates: the app's page, with a
        # smaller tab placed at ``origin`` on it.
        self.origin = (0.0, 0.0)
        self.offset_frames: dict[str, tuple[float, float]] = {}
        self.x, self.y = 1240.0, 640.0
        self._cdp = None
        self._origin = 0.0

    # -- log ---------------------------------------------------------------
    def event(self, kind: str, **data) -> None:
        self.events.append({"t": round(time.time() - self._origin, 4), "kind": kind, **data})

    def _stage(self, x: float, y: float) -> tuple[float, float]:
        return x + self.origin[0], y + self.origin[1]

    def _stage_rect(self, rect: dict) -> dict:
        x, y = self._stage(rect["x"], rect["y"])
        return {"x": x, "y": y, "width": rect["width"], "height": rect["height"]}

    def _screencast(self) -> None:
        """Stream the current page's repaints into the shared frame log."""

        cdp = self.page.context.new_cdp_session(self.page)
        origin = self.origin

        def on_frame(params: dict) -> None:
            cdp.send("Page.screencastFrameAck", {"sessionId": params["sessionId"]})
            name = f"{len(self.frames):05d}.jpg"
            (self.frames_dir / name).write_bytes(base64.b64decode(params["data"]))
            if origin != (0.0, 0.0):
                self.offset_frames[name] = origin
            stamp = params.get("metadata", {}).get("timestamp") or time.time()
            self.frames.append([round(stamp - self._origin, 4), name])

        cdp.on("Page.screencastFrame", on_frame)
        cdp.send("Page.startScreencast", {
            "format": "jpeg", "quality": 92, "everyNthFrame": 1,
            "maxWidth": VIEWPORT[0] * SCALE, "maxHeight": VIEWPORT[1] * SCALE,
        })
        self._cdp = cdp

    def _first_frame(self, root: str) -> None:
        """Force one paint so the screencast starts with a frame."""

        count = len(self.frames)
        self.page.evaluate(f"() => {{ {root}.style.outline = '0 solid transparent'; }}")
        self.page.wait_for_timeout(400)
        self.page.evaluate(f"() => {{ {root}.style.removeProperty('outline'); }}")
        self.page.wait_for_timeout(300)
        if len(self.frames) == count:
            raise AssertionError("Chromium screencast produced no frames")

    def start(self) -> None:
        self._origin = time.time()
        self._screencast()
        self.page.mouse.move(self.x, self.y)
        self.event("move", x=self.x, y=self.y)
        self._first_frame("document.body")

    def switch_page(self, page: Page, *, origin: tuple[float, float] = (0.0, 0.0)) -> None:
        """Keep recording in ``page``, for example a downloaded file opened in a new tab.

        A tab smaller than the app sits at ``origin`` on the stage; its frames
        are placed there after the recording (``place_offset_frames``).
        """

        self._cdp.send("Page.stopScreencast")
        self.page = page
        self.origin = origin
        page.bring_to_front()
        self._screencast()
        size = page.viewport_size
        # The renderer keeps the camera inside this rect, and a portrait take
        # hides its overview inset while the rect is smaller than the stage.
        self.event("stage", rect={"x": origin[0], "y": origin[1], "width": size["width"], "height": size["height"]})
        x = min(max(self.x, origin[0]), origin[0] + size["width"])
        y = min(max(self.y, origin[1]), origin[1] + size["height"])
        if (x, y) != (self.x, self.y):
            self.x, self.y = x, y
            self.event("move", x=round(x, 1), y=round(y, 1))
        page.mouse.move(self.x - origin[0], self.y - origin[1])
        self._first_frame("document.documentElement")

    def place_offset_frames(self, size: tuple[int, int]) -> None:
        """Put each frame of a smaller tab where that tab sits on the stage."""

        for name, (ox, oy) in self.offset_frames.items():
            path = self.frames_dir / name
            with Image.open(path) as frame:
                stage = Image.new("RGB", size, "white")
                stage.paste(frame.convert("RGB"), (round(ox * SCALE), round(oy * SCALE)))
            stage.save(path, quality=92)

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
            size = self.page.viewport_size
            rect = {"x": 0, "y": 0, "width": size["width"], "height": size["height"]}
        rect = self._stage_rect(rect)
        cx, cy = _center(rect)
        subject = self.page.evaluate(_FEATURE_BOUNDS)
        subject = subject and self._stage_rect(subject)
        portrait = portrait and self._stage_rect(portrait)
        if self.portrait and portrait is None:
            # A portrait take frames any overview on the whole drawn figure.
            figure = self.figure()
            if figure and _contains(rect, figure):
                portrait = figure
        self.event("camera", cx=round(cx, 1), cy=round(cy, 1), zoom=zoom, duration=duration,
                   rect={key: round(rect[key], 1) for key in ("x", "y", "width", "height")},
                   subject=subject, portrait=portrait and {key: round(portrait[key], 1) for key in portrait})

    def figure(self) -> dict | None:
        """The preview figure's visible bounds on the stage, or None before the first result."""

        geometry = self.page.evaluate(_PREVIEW)
        if not geometry:
            return None
        canvas, svg = geometry["canvas"], geometry["svg"]
        left, top = max(svg["x"], canvas["x"]), max(svg["y"], canvas["y"])
        right = min(svg["x"] + svg["width"], canvas["x"] + canvas["width"])
        bottom = min(svg["y"] + svg["height"], canvas["y"] + canvas["height"])
        if right <= left or bottom <= top:
            return None
        return self._stage_rect({"x": left, "y": top, "width": right - left, "height": bottom - top})

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
        """Ease the pointer to page point ``(x, y)``."""

        x, y = self._stage(x, y)
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
        begin = time.monotonic()
        for index in range(1, steps + 1):
            # Each step costs time on a busy machine; skip ahead so the move
            # still takes about ``duration`` (a slow drag would pan slowly).
            share = min(1.0, max(index / steps, (time.monotonic() - begin) / duration))
            p = _ease(share)
            arc = 4 * p * (1 - p) * bow
            px, py = sx + (x - sx) * p + nx * arc, sy + (y - sy) * p + ny * arc
            self.page.mouse.move(px - self.origin[0], py - self.origin[1])
            self.event("move", x=round(px, 1), y=round(py, 1))
            if share >= 1.0:
                break
            self.page.wait_for_timeout(duration * 1000 / steps)
        self.x, self.y = x, y

    def _press(self, x: float, y: float, *, real: bool) -> None:
        self.move_to(x, y)
        self.page.wait_for_timeout(110)
        sx, sy = self._stage(x, y)
        self.event("down", x=round(sx, 1), y=round(sy, 1))
        if real:
            self.page.mouse.down()
        self.page.wait_for_timeout(85)
        if real:
            self.page.mouse.up()
        self.event("up", x=round(sx, 1), y=round(sy, 1))

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

    def drag(self, start: tuple[float, float], end: tuple[float, float], *, duration: float = 0.9) -> None:
        """Press at ``start`` and pull to ``end``, for example to pan the preview."""

        self.move_to(*start)
        self.page.wait_for_timeout(110)
        (sx, sy), (ex, ey) = self._stage(*start), self._stage(*end)
        self.event("down", x=round(sx, 1), y=round(sy, 1))
        self.page.mouse.down()
        self.move_to(*end, duration=duration)
        # Rest before releasing: the preview applies a pan once a frame, and a
        # release before the drawing caught up lands on what passed under it.
        self.page.wait_for_timeout(150)
        self.page.mouse.up()
        self.event("up", x=round(ex, 1), y=round(ey, 1))
        self.hold(0.3)

    def choose(self, select: Locator, value: str, *, in_svg: bool = False, rows: int = MENU_ROWS,
               linger: float = 0.2, keeps: bool = True) -> None:
        """Select an option.

        Headless Chromium never paints native popups, so the option list and
        the chosen row are logged for the editor to draw. ``in_svg`` marks a
        control inside an SVG document, such as an opened Interactive SVG,
        where Playwright's form helpers do not recognize HTML controls.
        ``rows`` is the most rows the drawn list shows, and ``linger`` how long
        it stays open before the row is pointed at. ``keeps`` False marks an
        action option, such as Load font file…, after which the select shows
        its previous value again.
        """

        box = select.bounding_box()
        options = select.evaluate(
            "(node, value) => ({labels: Array.from(node.options).map((o) => o.label.trim()),"
            " selected: Array.from(node.options).findIndex((o) => o.value === value)})", value)
        if options["selected"] < 0:
            raise AssertionError(f"{value!r} is not an option of the select")
        labels, selected = options["labels"], options["selected"]
        rows = min(rows, len(labels))
        first = min(max(0, selected - rows // 2), len(labels) - rows)
        x, y = _center(box)
        self._press(x, y, real=False)
        sx, sy = self._stage(box["x"], box["y"])
        self.event("menu", x=sx, y=sy, width=box["width"], height=box["height"],
                   options=labels, selected=selected, window=[first, rows], duration=0.85 + linger)
        self.hold(linger)
        # Point at the chosen row of the drawn menu (22 CSS px rows below the box).
        row_y = box["y"] + box["height"] + 6 + (selected - first + 0.5) * 22
        self._press(box["x"] + min(box["width"], 120) * 0.5, row_y, real=False)
        if in_svg:
            if select.evaluate(_SET_VALUE, value) != value:
                raise AssertionError(f"{value!r} could not be selected")
        else:
            select.select_option(value)
            if keeps:
                expect(select).to_have_value(value)
        self.hold(0.45)

    def pick_color(self, swatch: Locator, value: str) -> None:
        box = swatch.bounding_box()
        x, y = _center(box)
        self._press(x, y, real=False)
        sx, sy = self._stage(box["x"], box["y"])
        self.event("swatch", x=sx, y=sy, width=box["width"], height=box["height"],
                   value=value, duration=0.8)
        self.hold(0.45)
        swatch.fill(value)
        self.hold(0.45)

    def type(self, field: Locator, text: str, *, delay: float = 0.075, in_svg: bool = False) -> None:
        self.click(field, pause=0.12)
        if not field.evaluate("node => node === document.activeElement"):
            raise AssertionError(f"Clicking {field} did not focus it")
        if field.evaluate("node => node.value"):
            self.page.keyboard.press("Control+A")
        for char in text:
            self.page.keyboard.type(char)
            self.page.wait_for_timeout(delay * 1000)
        if in_svg:
            if field.evaluate("node => node.value") != text:
                raise AssertionError(f"Typing did not set {field} to {text!r}")
        else:
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
        # The last cards cannot scroll up to the anchor; in view is enough.
        if not (area["y"] <= box["y"] and box["y"] + min(box["height"], 80) <= area["y"] + area["height"]):
            raise AssertionError(f"Could not scroll {locator} into the settings view: {box}")

    def open_section(self, summary: Locator, *, mark: str | None = None) -> None:
        """Scroll to a settings card and open it unless it is already open."""

        self.camera(self.page.locator(".settings-scroll").bounding_box(), zoom=2.0)
        self.reveal(summary)
        if mark:
            self.mark(mark)
        self.focus(summary, zoom=2.0)
        if not summary.evaluate("node => node.parentElement.open"):
            self.click(summary, left=70, pause=0.5)

    def upload(self, trigger: Locator, path: Path) -> None:
        """Click a file control and pick ``path`` in the (unpainted) chooser."""

        with self.page.expect_file_chooser() as chooser:
            self.click(trigger, pause=0.1)
        chooser.value.set_files(path)
        self.toast(path.name, "file")

    def set_checked(self, box: Locator, checked: bool, *, pause: float = 0.5) -> None:
        """Point at a checkbox and click it only when it needs to change."""

        if box.is_checked() != checked:
            self.click(box, pause=pause)
        else:
            self.move_to(*_center(box.bounding_box()))
            self.hold(pause)
        expect(box).to_be_checked(checked=checked)

    def generate(self, *, target: float | None = 1.4, timeout_ms: int = GENERATION_TIMEOUT_MS,
                 stages: bool = False, after: float = 0, close: bool = False) -> None:
        """Press Generate Diagram and wait for the result.

        ``target`` fast-forwards the wait to about that many seconds; None
        keeps it at normal speed. ``stages`` frames the stage line that the
        progress overlay shows while Generate runs. ``after`` plays that many
        seconds of the wait at normal speed before the fast-forward starts.
        ``close`` has a tall take frame the stage line alone, large enough to read.
        """

        def press(button: Locator) -> None:
            self.click(button, pause=0.05)
            status = self.page.locator("[data-generation-progress]")
            if stages:
                # The overlay is centered; its stage line changes width as it runs.
                status.wait_for(state="visible", timeout=ACTION_TIMEOUT_MS)
                cx, cy = _center(status.bounding_box())
                self.camera({"x": cx - 400, "y": cy - 190, "width": 800, "height": 300}, zoom=2.4, duration=0.6,
                            portrait={"x": cx - 170, "y": cy - 110, "width": 340, "height": 200} if close else None)
            else:
                # Pull back so the progress overlay and then the result are in view.
                self.camera(duration=0.8)
            if target is not None:
                self.hold(after)
                self.fast_forward(target=target)

        button = self.page.get_by_role("button", name=self.t("Generate Diagram"), exact=True)
        generate_and_wait_for_result(self.page, button=button, click=press, timeout_ms=timeout_ms)
        self.page.wait_for_function("() => document.fonts.status === 'loaded'")
        self.hold(0.35)
        if target is not None:
            self.normal_speed()

    # -- page queries ------------------------------------------------------
    def t(self, text: str, params: dict | None = None) -> str:
        """``text`` as the app shows it in its current language (the app's own ``t``).

        Locators use it for every UI name, since the app translates labels,
        aria-labels, and titles alike. Data such as gene names or file names is
        never translated.
        """

        return self.app.evaluate("([text, params]) => window.__GBDRAW_APP__.t(text, params ?? undefined)",
                                 [text, params])

    def tn(self, count: int, one: str, other: str, params: dict | None = None) -> str:
        """The app's plural text for ``count`` (the app's own ``tn``)."""

        return self.app.evaluate(
            "([count, one, other, params]) => window.__GBDRAW_APP__.tn(count, one, other, params ?? undefined)",
            [count, one, other, params])

    def text_bounds(self, pattern: str) -> dict:
        bounds = self.page.evaluate(_TEXT_BOUNDS, pattern)
        if not bounds:
            raise AssertionError(f"No preview text matches {pattern!r}")
        return bounds


def _meet_gbdraw(rec: Recorder, downloads: Path) -> Path:
    page = rec.page
    preview = page.locator(PREVIEW_REGION)

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
    rec.upload(upload, HUMAN_MITOCHONDRION_PATH)
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

    rec.caption(3, "Choose a track layout", "Layout \u203a Track Preset: Middle")
    rec.open_section(page.get_by_label("Layout", exact=True), mark="layout")
    preset = page.get_by_label("Track Preset", exact=True)
    rec.reveal(preset)
    rec.focus(preset, zoom=2.2)
    rec.choose(preset, "middle")
    rec.caption(3, "Put both strands on one track", "Layout \u203a Separate Strands")
    strands = page.get_by_label("Separate Strands", exact=True)
    rec.focus(strands, zoom=2.2, duration=0.6)
    rec.click(strands, pause=0.5)
    expect(strands).not_to_be_checked()
    rec.generate()
    rec.mark("layout-done")
    rec.camera(duration=1.0)
    rec.move_to(1500, 700)
    rec.hold(2.4)

    rec.caption(4, "Show feature labels", "Labels › Label Mode: Out")
    rec.open_section(page.get_by_label("Labels", exact=True))
    mode = page.locator("#circular-label-mode")
    rec.reveal(mode)
    rec.mark("labels")
    rec.focus(mode, zoom=2.1)
    rec.choose(mode, "out")
    rec.generate()
    rec.camera(duration=1.0)
    rec.hold(1.4)
    rec.caption(4, "Long product names crowd the map", "")
    right = rec.text_bounds(r"^NADH dehydrogenase subunit [12]$")
    rec.mark("crowded")
    rec.camera(right, zoom=1.8, duration=1.1)
    rec.move_to(1700, 480)
    rec.hold(2.4)

    rec.caption(4, "Prefer short gene symbols", "Qualifier priority: CDS › gene")
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

    rec.caption(5, "Color genes by function", "Colors › Specific rules match the gene qualifier")
    rec.open_section(page.get_by_label("Colors", exact=True))
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

    rec.caption(6, "Mark a region, even across the origin", "Region Annotations › D-loop, 16,024–576")
    rec.open_section(page.get_by_label("Region Annotations", exact=True))
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
    dloop = page.locator(f"{PREVIEW_REGION} svg [data-gbdraw-annotation-id]").first
    rec.mark("dloop")
    rec.focus(dloop, zoom=2.4, duration=1.2)
    rec.move_to(*_center(dloop.bounding_box()))
    rec.hold(2.6)
    rec.camera(duration=1.0)
    rec.hold(0.6)

    rec.caption(7, "Click any feature for details", "")
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

    return export_svg(rec, downloads, 8)


def export_svg(rec: Recorder, downloads: Path, step: int) -> Path:
    """Download the result through the SVG button; every journey ends here."""

    page = rec.page
    rec.caption(step, "Export a vector SVG", "PNG and PDF are one click away")
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


def export_interactive_svg(rec: Recorder, downloads: Path, step: int | None) -> Path:
    """Download the result through the Interactive SVG button."""

    page = rec.page
    rec.caption(step, "Export an Interactive SVG", "Every feature's details, in one file")
    # Icon glyphs precede the button text in the accessible name.
    button = page.get_by_role("button", name=re.compile(rf"{re.escape(rec.t('Interactive SVG'))}$"))
    rec.mark("export")
    rec.focus(button, zoom=2.1)
    rec.hold(0.4)
    with page.expect_download(timeout=GENERATION_TIMEOUT_MS) as info:
        rec.click(button, pause=0.1)
    download = info.value
    if download.failure() is not None:
        raise AssertionError(f"Interactive SVG download failed: {download.failure()}")
    downloads.mkdir(parents=True, exist_ok=True)
    exported = downloads / download.suggested_filename
    download.save_as(exported)
    if 'data-gbdraw-interactive-svg="true"' not in exported.read_text(encoding="utf-8")[:20000]:
        raise AssertionError(f"{exported.name} is not an interactive SVG")
    rec.toast(f"{exported.name} downloaded ({exported.stat().st_size // 1024} KB)", "download")
    rec.mark("downloaded")
    rec.hold(1.8)
    return exported


def open_in_browser(rec: Recorder, svg: Path) -> list[str]:
    """Open a downloaded file in a new tab and keep recording there.

    The file is served at a same-origin address so the capture stays
    network-isolated. Returns the tab's page errors and console errors, which
    the journey checks when it is done with the file.
    """

    app = rec.page
    url = re.match(r"^https?://[^/]+", app.url).group(0) + f"/downloads/{svg.name}"
    page = app.context.new_page()
    if rec.portrait:
        # The file fits its drawing to the window, so a tall tab gives a tall drawing.
        page.set_viewport_size({"width": PORTRAIT_TAB[0], "height": PORTRAIT_TAB[1]})
    errors: list[str] = []
    page.on("pageerror", lambda error: errors.append(str(error)))
    page.on("console", lambda message: errors.append(message.text) if message.type == "error" else None)
    page.route(url, lambda route: route.fulfill(path=str(svg), content_type="image/svg+xml"))
    # A large file takes seconds to open; the previous tab stays on screen meanwhile.
    rec.fast_forward(target=1.2)
    page.goto(url, wait_until="load", timeout=GENERATION_TIMEOUT_MS)
    # The embedded script builds its controls once it has read the metadata.
    expect(page.get_by_role("button", name="Expand feature search", exact=True)).to_be_visible(
        timeout=GENERATION_TIMEOUT_MS)
    page.wait_for_function("() => document.fonts.status === 'loaded'")
    if rec.portrait:
        rec.switch_page(page, origin=((VIEWPORT[0] - PORTRAIT_TAB[0]) / 2, 0.0))
        rec.camera(duration=0)   # cut to the tab
    else:
        rec.switch_page(page)
    rec.normal_speed()
    rec.toast(svg.name, "file")
    return errors


def switch_language(rec: Recorder, code: str) -> None:
    """Switch the app's UI language on camera with the header's Language menu."""

    page = rec.page
    menu = page.locator("select[data-ui-language]")
    rec.camera(duration=0)
    rec.caption(None, "Switch the app to Simplified Chinese", "Language, in the header")
    rec.hold(1.6)
    rec.mark("language")
    rec.focus(menu, zoom=2.6)
    rec.hold(0.5)
    rec.choose(menu, code)
    page.wait_for_function("(code) => window.__GBDRAW_APP__.uiLocale.value === code", arg=code)
    page.wait_for_function("() => document.fonts.status === 'loaded'")
    rec.hold(1.4)
    rec.camera(duration=0.9)
    rec.hold(1.6)


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


@dataclass(frozen=True)
class Journey:
    """One recorded story: its inputs, its steps, and what its export must show."""

    name: str
    sources: tuple[Path, ...]
    run: Callable[[Recorder, Path], Path]   # performs the steps and returns the exported SVG
    inspect: Callable[[Path], dict]         # raises when the export lacks the recorded edits
    zoom_text: str                          # exported text that the vector zoom centers on
    mode: str | None = None                 # selected before recording; None shows the choice
    gallery: bool = False                   # render the Gallery tiles for the highlights
    vector_zoom: bool = True                # end on the 7x zoom into the exported file
    portrait: bool = False                  # a take shaped for the tall video
    # How a take in another UI language switches the app to it on camera, before ``run``.
    switch: Callable[[Recorder, str], None] = switch_language


MEET_GBDRAW = Journey("meet-gbdraw", (HUMAN_MITOCHONDRION_PATH,), _meet_gbdraw, _inspect_export,
                      zoom_text="D-loop", mode="Circular", gallery=True)


def journeys() -> dict[str, Journey]:
    from video.promo import PROMOS  # built on this module, so imported late
    from video.tutorials import TUTORIALS

    return {journey.name: journey for journey in (MEET_GBDRAW, *TUTORIALS, *PROMOS)}


def _render_finale(browser, svg: Path, out: Path, size: tuple[int, int], text: str) -> dict:
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
            (text) => {
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
                .find((node) => node.textContent === text);
              if (!target) return null;
              const box = target.getBoundingClientRect();
              return {base, cx: base[0] + (box.left + box.width / 2 - offX) / scale,
                      cy: base[1] + (box.top + box.height / 2 - offY) / scale};
            }
            """,
            text,
        )
        if geometry is None:
            raise AssertionError(f"{svg.name} has no {text!r} text for the vector zoom")
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


def record_walkthrough(run: Path, name: str) -> Path:
    """Record journey ``name`` into ``run/raw/walkthrough`` and return its manifest."""

    journey = journeys()[name]
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
            if journey.mode:
                mode = page.get_by_role("button", name=journey.mode, exact=True)
                mode.click()
                expect(mode).to_have_attribute("aria-pressed", "true")
            handle = page.get_by_title("Drag to resize", exact=True).bounding_box()
            page.mouse.move(*_center(handle))
            page.mouse.down()
            page.mouse.move(SETTINGS_RIGHT_EDGE, _center(handle)[1], steps=8)
            page.mouse.up()
            page.wait_for_function("() => document.fonts.status === 'loaded'")
            recorder = Recorder(page, raw, portrait=journey.portrait)
            recorder.start()
            language = STORIES[name].language
            if language != "en":
                journey.switch(recorder, language)
            exported = journey.run(recorder, raw / "downloads")
            recorder.stop()
            capture.assert_clean()
        finally:
            capture.close()
        recorder.place_offset_frames((VIEWPORT[0] * SCALE, VIEWPORT[1] * SCALE))
        report = journey.inspect(exported)
        browser = playwright.chromium.launch(headless=True)
        try:
            finale = {folder: _render_finale(browser, exported, raw / folder, size, journey.zoom_text)
                      for folder, size in (("finale", FINALE_SIZE), ("finale-portrait", FINALE_PORTRAIT_SIZE))
                      if journey.vector_zoom}
            _rasterize(browser, exported, raw / "final-figure.png")
            gallery = []
            for stem, label in GALLERY if journey.gallery else ():
                source = GALLERY_SOURCES / f"{stem}.svg"
                target = raw / "gallery" / f"{stem}.png"
                target.parent.mkdir(exist_ok=True)
                _rasterize(browser, source, target, size=(1600, 1000), scale=1.5)
                gallery.append({"path": str(target.relative_to(raw)), "label": label,
                                "source": str(source.relative_to(REPO_ROOT)), "source_sha256": sha256(source)})
        finally:
            browser.close()
    # A frame from a tab that stopped streaming can arrive after the next tab's first frame.
    (raw / "frames.json").write_text(json.dumps(sorted(recorder.frames)), encoding="utf-8")
    (raw / "events.json").write_text(json.dumps(recorder.events, indent=1), encoding="utf-8")
    manifest = {
        "schema_version": 1, "journey": journey.name,
        "viewport": list(VIEWPORT), "scale": SCALE,
        "sources": [{"path": source.name, "sha256": sha256(source)} for source in journey.sources],
        "export": {"path": str(exported.relative_to(raw)), "sha256": sha256(exported), "semantics": report},
        "frames": len(recorder.frames), "events": len(recorder.events),
        "duration": recorder.events[-1]["t"], "finale": finale, "gallery": gallery,
    }
    path = raw / "walkthrough.json"
    path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return path
