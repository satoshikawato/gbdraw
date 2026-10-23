"""Record the gbdraw 0.14.0 promotion: its new features, with the app in Simplified Chinese.

One tall take (``walkthrough_render.PROMO_ID``) of the desktop app (1920x1080),
framed for phone feeds. The header Language menu switches the app to 简体中文
on camera; then two Gallery examples, opened with Load example, show the new
features as they are: a loaded font with a Chinese plot title, a record rotated
to start at a feature, a feature moved to a lane while Generate shows its
stages, Save Session with Run info, and a similarity-group alignment. Captions are English source text; the take shows them
only in zh-Hans (``docs/videos/meet-gbdraw/captions.zh-Hans.json``).
"""

from __future__ import annotations

import json
import os
import re
from pathlib import Path

from playwright.sync_api import Locator, expect

from config import ACTION_TIMEOUT_MS, GENERATION_TIMEOUT_MS
from video.tutorials import (_button, _close_feature, _feature_point, _field, _fit_preview, _frame_box, _frame_figure,
                             _label, _look, _open_feature, _pick, _summary, _svg_id)
from video.walkthrough import PREVIEW_REGION, Journey, Recorder, _center
from video.walkthrough_render import FONT_CACHE, PORTRAIT_PAD, PROMO_ID

HMMT_EXAMPLE = "HmmtDNA_basic_circular"
BGC_EXAMPLE = "BGC0000708-BGC0000713"
BGC_RECORDS = ("BGC0000708", "BGC0000709", "BGC0000711", "BGC0000712", "BGC0000713")
TITLE = "人类线粒体基因组"
LANGUAGES = 11          # English and the ten translations
LANE_GENE = "ATP6"      # a gene whose outside label and leader line show the move
LANE_LABEL = "ATP synthase F0 subunit 6"
STAGES_SECONDS = 2.0   # how long Generate's stage line plays at normal speed
MENU_ROW = 22

_FONT_READY = """
() => (document.querySelector('select[data-font-select]')?.selectedOptions[0]?.label || '') === 'Noto Sans SC'
"""
_LOAD_OR_ASK = """
([before, name]) => { const app = window.__GBDRAW_APP__;
  const asked = [...document.querySelectorAll('[role="dialog"] button')]
    .some((button) => button.offsetParent && button.textContent.trim() === name);
  return asked || app.processing || String(app?.results?.[0]?.content || '') !== before; }
"""
_RESULT_CHANGED = """
(before) => { const app = window.__GBDRAW_APP__; const content = String(app?.results?.[0]?.content || '');
  return Boolean(app && !app.processing && !app.similarityAlignmentBusy && content && content !== before); }
"""
_EXAMPLE_LOADED = """
([mode, before]) => { const app = window.__GBDRAW_APP__; const content = String(app?.results?.[0]?.content || '');
  return Boolean(app && app.mode === mode && !app.processing && content && content !== before); }
"""
_BOUNDS = """
([svgId, label]) => {
  const region = document.querySelector('$REGION');
  const boxes = [...region.querySelectorAll(`svg [data-gbdraw-feature-id="${svgId}"]`),
                 ...[...region.querySelectorAll('svg text')].filter((node) => node.textContent === label)]
    .map((node) => node.getBoundingClientRect()).filter((box) => box.width && box.height);
  if (!boxes.length) return null;
  const left = Math.min(...boxes.map((b) => b.left)), top = Math.min(...boxes.map((b) => b.top));
  const right = Math.max(...boxes.map((b) => b.right)), bottom = Math.max(...boxes.map((b) => b.bottom));
  return {x: left, y: top, width: right - left, height: bottom - top};
}
""".replace("$REGION", PREVIEW_REGION)


_GROUP_BOUNDS = """
(group) => {
  const region = document.querySelector('$REGION');
  const ids = (window.__GBDRAW_APP__?.extractedFeatures || []).filter((item) => item.orthogroupId === group)
    .map((item) => item.svg_id);
  const boxes = ids.flatMap((id) => [...region.querySelectorAll(`svg [data-gbdraw-feature-id="${id}"]`)])
    .map((node) => node.getBoundingClientRect()).filter((box) => box.width && box.height);
  if (!boxes.length) return null;
  const top = Math.min(...boxes.map((b) => b.top)) - 40, bottom = Math.max(...boxes.map((b) => b.bottom)) + 40;
  const middle = (Math.min(...boxes.map((b) => b.left)) + Math.max(...boxes.map((b) => b.right))) / 2;
  const width = Math.max(320, (bottom - top) * 0.88);   // a column as wide as the tall window shows it
  return {x: middle - width / 2, y: top, width, height: bottom - top};
}
""".replace("$REGION", PREVIEW_REGION)


def _font_file() -> Path:
    """Noto Sans SC Regular, from the same folders the renderer reads its Chinese captions from."""

    root = Path(os.environ.get("GBDRAW_VIDEO_FONT_DIR") or FONT_CACHE)
    path = root / "NotoSansSC-Regular.otf"
    if not path.is_file():
        raise FileNotFoundError(f"The promotion loads {path} into the app's Font setting")
    return path


def _content(rec: Recorder) -> str:
    return rec.page.evaluate("() => String(window.__GBDRAW_APP__?.results?.[0]?.content || '')")


def _fit(rec: Recorder) -> None:
    """Zoom the new figure out to fit, skimmed: Generate leaves it at 100 %, larger than the tall window."""

    _fit_preview(rec, skim=0.8)


_IDLE = """
() => new Promise((done) => setTimeout(() => requestAnimationFrame(() => requestAnimationFrame(done)), 500))
"""


def _settle(rec: Recorder) -> None:
    """Skim the app's own work after an edit (a live preview render, hover on a large figure).

    The page cannot take pointer moves meanwhile, so the pointer would stand
    still on camera for a second or two.
    """

    rec.fast_forward(target=0.2)
    rec.page.evaluate(_IDLE)
    rec.normal_speed()


def _generate(rec: Recorder, *, target: float, stages: bool = False, after: float = 0) -> None:
    """Generate Diagram, then fit the new figure and show the whole page."""

    _settle(rec)
    button = _button(rec, "Generate Diagram")
    rec.focus(button, zoom=1.9, duration=0.7)
    rec.hold(0.2)
    rec.generate(target=target, stages=stages, after=after, close=stages)
    _fit(rec)
    rec.camera(duration=0.9)


def _switch_to_chinese(rec: Recorder, code: str) -> None:
    """Steps lang and lang2: the header Language menu, open on every language, then 简体中文."""

    page = rec.page
    menu = page.locator("select[data-ui-language]")
    expect(menu.locator("option")).to_have_count(LANGUAGES)
    box = menu.bounding_box()
    # The menu and the list of every language drawn below it.
    area = {"x": box["x"] - 30, "y": box["y"] - 10, "width": box["width"] + 60,
            "height": box["height"] + MENU_ROW * LANGUAGES + 34}
    rec.camera(area, zoom=2.4, duration=0, portrait=area)
    rec.caption(1, "The interface in Simplified Chinese",
                "The Language menu in the header switches at once, without reloading")
    rec.hold(1.7)
    rec.caption(1, "11 interface languages", "Including Simplified and Traditional Chinese")
    rec.choose(menu, code, rows=LANGUAGES, linger=2.0)
    page.wait_for_function("(code) => window.__GBDRAW_APP__.uiLocale.value === code", arg=code)
    page.wait_for_function("() => document.fonts.status === 'loaded'")
    rec.hold(0.6)
    rec.camera(duration=0.8)
    rec.hold(0.6)


def _wheel_into_view(rec: Recorder, area: dict, target: Locator) -> None:
    """Wheel-scroll the list under ``area`` until ``target`` is wholly inside that area."""

    rec.move_to(area["x"] + area["width"] * 0.5, area["y"] + area["height"] * 0.6)
    for _ in range(80):
        box = target.bounding_box()
        if area["y"] <= box["y"] and box["y"] + box["height"] <= area["y"] + area["height"] - 2:
            break
        rec.page.mouse.wheel(0, 50)
        rec.page.wait_for_timeout(40)
    else:
        raise AssertionError(f"Could not scroll {target} into {area}")
    rec.hold(0.4)


def _load_example(rec: Recorder, button: Locator, example_id: str, mode: str) -> None:
    """Load example › ``example_id``; a question about replacing unsaved work is answered with Load example."""

    page = rec.page
    before = _content(rec)
    rec.focus(button, zoom=2.2, duration=0.7)
    rec.hold(0.2)
    rec.click(button, pause=0.3)
    row = page.locator(f'[data-gallery-example="{example_id}"]')
    expect(row).to_be_visible()
    chooser = page.get_by_role("dialog").filter(has=row)
    dialog = chooser.bounding_box()
    rec.camera(dialog, zoom=1.4, duration=0.7, portrait=dialog)
    rec.hold(0.4)
    # The examples scroll inside their list, not the dialog.
    _wheel_into_view(rec, row.locator("xpath=ancestor::ul[1]").bounding_box(), row)
    rec.click(row, pause=0.1)
    rec.fast_forward(target=0.8)
    replace = page.locator('[role="dialog"]:visible').get_by_role("button", name=rec.t("Load example"), exact=True)
    page.wait_for_function(_LOAD_OR_ASK, arg=[before, rec.t("Load example")], timeout=GENERATION_TIMEOUT_MS)
    if replace.count():
        rec.normal_speed()
        question = replace.locator("xpath=ancestor::*[@role='dialog'][1]").bounding_box()
        rec.camera(question, zoom=1.8, duration=0.6, portrait=question)
        rec.hold(0.5)
        rec.click(replace, pause=0.2)
        rec.fast_forward(target=0.8)
    rec.camera(duration=0.8)
    page.wait_for_function(_EXAMPLE_LOADED, arg=[mode, before], timeout=GENERATION_TIMEOUT_MS)
    page.wait_for_function("() => document.fonts.status === 'loaded'")
    rec.hold(0.4)
    rec.normal_speed()
    _fit(rec)


def _font_and_title(rec: Recorder) -> None:
    """Step font: Font › Load font file… (Noto Sans SC), then Plot Title 人类线粒体基因组 at the top."""

    page = rec.page
    rec.caption(3, "Any font for the figure text", "Load a font file; the plot title can be in Chinese too")
    rec.open_section(_summary(rec, "Font"))
    select = page.locator("select[data-font-select]")
    rec.reveal(select)
    rec.focus(select, zoom=2.2, duration=0.6)
    font = _font_file()
    with page.expect_file_chooser() as chooser:
        rec.choose(select, "load-file", keeps=False)
    chooser.value.set_files(font)
    rec.toast(font.name, "file")
    # The first font is read once the diagram runtime has started.
    rec.fast_forward(target=0.8)
    page.wait_for_function(_FONT_READY, timeout=GENERATION_TIMEOUT_MS)
    rec.hold(0.3)
    rec.normal_speed()
    rec.hold(0.2)
    rec.open_section(_summary(rec, "Titles and Record Labels"))
    _field(rec, page.get_by_role("textbox", name=rec.t("Plot Title"), exact=True), TITLE, zoom=2.0)
    _pick(rec, _label(rec, "Plot Title Position"), "top")
    _generate(rec, target=0.8)
    rec.mark("title")
    rec.hold(1.0)
    _look(rec, rec.text_bounds(f"^{TITLE}$"), zoom=2.4, seconds=2.0, pad=50)


def _rotate(rec: Recorder) -> None:
    """Step rotate: COX1 › Rotate record using this feature › Start of the record › Apply now."""

    page = rec.page
    rec.caption(4, "Start a circular record at any base",
                "Rotate record using this feature; the sequence and its coordinates stay the same")
    popup = _open_feature(rec, _feature_point(rec, "(item) => item.gene === arg", "COX1"), "Edit")
    toggle = popup.locator("#feature-record-actions-toggle")
    toggle.evaluate("node => node.scrollIntoView({block: 'center'})")
    page.wait_for_timeout(300)
    _frame_box(rec, popup.bounding_box())
    rec.click(toggle, pause=0.4)
    start = popup.locator("input[name='feature-record-rotation-position'][value='start']")
    start.evaluate("node => node.scrollIntoView({block: 'center'})")
    page.wait_for_timeout(300)
    rec.set_checked(start, True, pause=0.4)
    expect(popup.locator("#feature-record-rotation-preview")).to_contain_text("5,904")
    apply = _button(rec, "Apply now", within=popup)
    apply.evaluate("node => node.scrollIntoView({block: 'end'})")
    page.wait_for_timeout(300)
    _frame_box(rec, popup.bounding_box())
    rec.hold(0.9)
    before = _content(rec)
    rec.click(apply, pause=0.2)
    rec.fast_forward(target=0.8)
    page.wait_for_function(_RESULT_CHANGED, arg=before, timeout=GENERATION_TIMEOUT_MS)
    rec.hold(0.3)
    rec.normal_speed()
    rec.hold(0.6)
    _close_feature(rec)
    rec.mark("rotated")
    _frame_figure(rec)
    rec.hold(2.0)
    _look(rec, rec.text_bounds(r"^cytochrome c oxidase subunit I$"), zoom=2.4, seconds=1.5, pad=110)


def _lane(rec: Recorder) -> None:
    """Steps lane and stages: ATP6 › Feature placement › Outer lane 1, then Generate at normal speed."""

    page = rec.page
    rec.caption(5, "Put a feature on the lane you choose", "Its label and leader line move with it")
    svg_id = _svg_id(rec, "(item) => item.gene === arg", LANE_GENE)
    detail = page.evaluate(_BOUNDS, [svg_id, LANE_LABEL])
    _look(rec, detail, zoom=2.4, seconds=0.7, pad=60)
    popup = _open_feature(rec, _feature_point(rec, "(item) => item.gene === arg", LANE_GENE), "Edit")
    placement = _label(rec, "Feature placement", within=popup)
    placement.evaluate("node => node.scrollIntoView({block: 'center'})")
    page.wait_for_timeout(300)
    rec.focus(placement, zoom=2.4, duration=0.6)
    rec.choose(placement, "outward")
    rec.hold(0.4)
    _close_feature(rec)
    rec.caption(5, "Generate shows each stage as it runs", "Preparing, rendering: progress at a glance")
    # The stage line plays at normal speed; a long render then skims under the badge.
    _generate(rec, target=0.8, stages=True, after=STAGES_SECONDS)
    rec.mark("lane")
    _look(rec, page.evaluate(_BOUNDS, [svg_id, LANE_LABEL]), zoom=2.4, seconds=2.4, pad=60)


def _reset_alignment(rec: Recorder) -> None:
    """The example is already aligned on og_6: Editor › Similarity groups › Reset alignment…, fast-forwarded."""

    page = rec.page
    rec.fast_forward(target=1.6)
    drawer = page.locator("button.drawer-toggle")
    rec.focus(drawer, zoom=2.0, duration=0.5)
    rec.click(drawer, pause=0.6)
    tab = page.locator("button[aria-controls='right-drawer-panel-orthogroups']")
    rec.click(tab, pause=0.5)
    reset = page.locator("[data-similarity-alignment-reset]")
    rec.focus(reset, zoom=2.0, duration=0.5)
    rec.click(reset, pause=0.5)
    dialog = page.locator("[data-similarity-alignment-reset-dialog]")
    expect(dialog).to_be_visible()
    rec.camera(dialog.bounding_box(), zoom=1.6, duration=0.5, portrait=dialog.bounding_box())
    rec.hold(0.6)
    before = _content(rec)
    rec.click(dialog.locator("footer button.btn-primary"), pause=0.2)
    page.wait_for_function(_RESULT_CHANGED, arg=before, timeout=GENERATION_TIMEOUT_MS)
    rec.hold(0.3)
    rec.click(drawer, pause=0.5)
    rec.normal_speed()
    _fit(rec)


def _align(rec: Recorder) -> None:
    """Step align: livU of BGC0000708 › Similarity group og_6 › Align…"""

    page = rec.page
    rec.caption(7, "Align records by a similarity group",
                "Click a gene, choose Align…, and the clusters line up in one step")
    # The header's Load example; the empty preview's button left with the first example.
    _load_example(rec, _button(rec, "Load example").first, BGC_EXAMPLE, "linear")
    _reset_alignment(rec)
    _frame_figure(rec)
    rec.hold(0.7)
    # livU of the first record (BGC0000708), in the group the example was aligned on.
    livu = _feature_point(rec, "(item) => item.gene === arg && item.orthogroupId === 'og_6'", "livU")
    _settle(rec)
    rec.click_point(livu["x"], livu["y"], pause=0.4)
    popup = page.locator(".feature-popup")
    expect(popup).to_be_visible()
    group = popup.locator("[data-feature-popup-group='similarity']")
    group.evaluate("node => node.scrollIntoView({block: 'center'})")
    page.wait_for_timeout(300)
    expect(group).to_contain_text("og_6")
    _frame_box(rec, group.bounding_box(), zoom=2.0)
    rec.hold(0.5)
    align = group.locator("button").filter(has_text=rec.t("Align…")).first
    rec.focus(align, zoom=2.4, duration=0.5)
    before = _content(rec)
    rec.click(align, pause=0.2)
    rec.fast_forward(target=0.8)
    page.wait_for_function(_RESULT_CHANGED, arg=before, timeout=GENERATION_TIMEOUT_MS)
    rec.hold(0.3)
    rec.normal_speed()
    if popup.is_visible():
        _close_feature(rec)
    _fit(rec)
    rec.mark("aligned")
    _frame_figure(rec)
    rec.hold(2.0)
    # The og_6 genes now share one column across the five clusters.
    _look(rec, page.evaluate(_GROUP_BOUNDS, "og_6"), zoom=2.4, seconds=2.2)


_TEXT_BOX = """
(node) => { const range = document.createRange(); range.selectNodeContents(node);
  const box = range.getBoundingClientRect(); return {x: box.left, y: box.top, width: box.width, height: box.height}; }
"""


def _session(rec: Recorder, downloads: Path) -> None:
    """Step session, on the human mitochondrial map after its last Generate: Save Session, then Run info.

    Run info shows the commands of the last Generate; the camera frames the
    Source recipe and Exact replay commands, not the notes below them.
    """

    page = rec.page
    rec.caption(6, "Save a Session and continue any time", "Run info gives the command line that redraws the figure")
    save = page.locator(f"button[aria-label='{rec.t('Save Session')}']")
    rec.focus(save, zoom=2.4, duration=0.8)
    rec.hold(0.2)
    with page.expect_download(timeout=ACTION_TIMEOUT_MS) as info:
        rec.click(save, pause=0.1)
    download = info.value
    if download.failure() is not None:
        raise AssertionError(f"Save Session failed: {download.failure()}")
    downloads.mkdir(parents=True, exist_ok=True)
    session = downloads / download.suggested_filename
    download.save_as(session)
    rec.toast(f"{session.name} downloaded ({session.stat().st_size // 1024} KB)", "download")
    rec.hold(0.9)
    tab = page.get_by_role("button", name=rec.t("Run info"), exact=True)
    _settle(rec)
    rec.focus(tab, zoom=2.4, duration=0.6)
    rec.click(tab, pause=0.6)
    source_heading = page.get_by_text(rec.t("Source recipe"), exact=True)
    replay_heading = page.get_by_text(rec.t("Exact replay"), exact=True)
    # Each heading row's block holds the command, or, for the Source recipe, why it is unavailable.
    source = source_heading.locator("xpath=../../*[last()]")
    replay = replay_heading.locator("xpath=../../pre")
    expect(replay).to_contain_text(" --session ")
    shown = page.evaluate("""() => { const info = window.__GBDRAW_APP__.lastRunInfo || {};
      return {source_available: info.sourceRecipe?.available !== false,
              source: String(info.sourceRecipe?.command || info.sourceRecipe?.unavailableReason || info.command || ''),
              replay: String(info.exactReplay?.command || info.sessionCommand || '')}; }""")
    (downloads / "run-info.json").write_text(json.dumps(shown, ensure_ascii=False, indent=1), encoding="utf-8")
    boxes = [source_heading.bounding_box(), replay_heading.bounding_box(), replay.bounding_box(),
             source.evaluate(_TEXT_BOX) if shown["source_available"] else source.bounding_box(),
             replay.evaluate(_TEXT_BOX)]
    # The tall window shows a column as wide as the Exact replay command, from the Save Session notice
    # down to the Exact replay box; the notes below it fall under the inset or out of the shot. The
    # renderer pads a framed rect by PORTRAIT_PAD, so the rect is that column narrowed by it.
    notice = page.get_by_text(session.name).first.bounding_box()
    left = min(notice["x"], *(box["x"] for box in boxes)) - 12
    width = boxes[4]["x"] + boxes[4]["width"] + 12 - left
    top, bottom = notice["y"] - 8, boxes[2]["y"] + boxes[2]["height"] + 8
    inner = width / PORTRAIT_PAD
    area = {"x": left + (width - inner) / 2, "y": top, "width": inner, "height": bottom - top}
    rec.camera(area, zoom=2.4, duration=0.9, portrait=area)
    rec.move_to(*_center(boxes[4]))
    rec.hold(3.0)
    # Back to the preview for the next example, skimmed.
    rec.fast_forward(target=0.4)
    rec.click(_button(rec, "Preview"), pause=0.2)
    rec.normal_speed()


def _promo(rec: Recorder, downloads: Path) -> Path:
    page = rec.page
    rec.caption(2, "Load an example", "The human mitochondrial genome opens as a finished map")
    _load_example(rec, page.locator("button[aria-describedby='example-session-help']"), HMMT_EXAMPLE, "circular")
    _frame_figure(rec)
    rec.hold(2.0)
    _font_and_title(rec)
    _rotate(rec)
    _lane(rec)
    _session(rec, downloads)
    _align(rec)
    figure = _content(rec)
    # The outro shows the last figure: the Result the app drew, as the SVG it holds.
    exported = downloads / f"{BGC_EXAMPLE}.svg"
    exported.write_text(figure, encoding="utf-8")
    return exported


def _inspect_promo(svg: Path) -> dict:
    text = svg.read_text(encoding="utf-8")
    report = {
        "records": sorted(set(BGC_RECORDS) & set(re.findall(r"BGC\d{7}", text))),
        "similarity_links": text.count('data-match-kind="orthogroup"'),
        "session": sorted(path.name for path in svg.parent.glob("*.gbdraw-session.json*")),
        "run_info": json.loads((svg.parent / "run-info.json").read_text(encoding="utf-8")),
    }
    if len(report["records"]) != len(BGC_RECORDS) or not report["similarity_links"] or not report["session"]:
        raise AssertionError(f"The promotion's last Result or Session is missing: {report}")
    return report


PROMOS = (Journey(PROMO_ID, (), _promo, _inspect_promo, zoom_text="", vector_zoom=False, portrait=True,
                  switch=_switch_to_chinese),)
