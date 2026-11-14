from __future__ import annotations

import io

import pytest
from PIL import Image

import tools.capture_gallery_tutorial_screenshots as capture_tool


pytestmark = pytest.mark.browser


def test_open_select_overlay_keeps_selected_option_of_a_long_list_in_the_crop() -> None:
    playwright_api = pytest.importorskip("playwright.sync_api")
    labels = [f"palette_{index:02d}" for index in range(55)]
    options = "".join(f"<option>{label}</option>" for label in labels)
    html = (
        "<html><body style='margin:0'><div class='card'><div class='card-header'>Colors</div>"
        f"<select>{options}</select></div></body></html>"
    )
    with playwright_api.sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        try:
            page = browser.new_page(viewport={"width": 800, "height": 600})
            page.set_content(html)
            capture_tool.apply_open_select_overlays(
                page, {"selectedLabel": "palette_40", "selector": "select", "placement": "right"}
            )
            clip = capture_tool.capture_open_select_clip(page, capture_tool.normalize_crop_padding(12))
            selected = page.evaluate(
                """() => {
                  const rect = document.querySelector('.gbdraw-capture-select-option.is-selected').getBoundingClientRect();
                  const menu = document.querySelector('.gbdraw-capture-select-menu').getBoundingClientRect();
                  return {top: rect.top + window.scrollY, bottom: rect.bottom + window.scrollY,
                          menuTop: menu.top + window.scrollY, menuBottom: menu.bottom + window.scrollY};
                }"""
            )
        finally:
            browser.close()

    assert selected["menuTop"] <= selected["top"] and selected["bottom"] <= selected["menuBottom"]
    assert clip["y"] <= selected["top"] and selected["bottom"] <= clip["y"] + clip["height"]
    assert clip["height"] < 600


def test_a_capture_hides_what_lies_around_its_subject_and_moves_nothing() -> None:
    """An open-select crop hides the page beside the Settings pane (OV-3395: its edge cut
    figure text and the Layout button), and an open popup hides the Settings pane under it;
    neither moves anything."""

    playwright_api = pytest.importorskip("playwright.sync_api")
    html = (
        "<html><body style='margin:0'><header class='app-header'>gbdraw</header>"
        "<style>footer[data-app-footer] { anchor-name: --app-footer; }"
        ".generate-bar { position: fixed; left: 0; width: 300px; height: 40px; position-anchor: --app-footer;"
        " bottom: anchor(top); }</style>"
        "<main><div class='settings-pane' style='width:300px'><select><option>A</option></select>"
        "<div class='generate-bar'>Generate Diagram</div></div>"
        "<div id='result-pane'>figure</div></main>"
        "<footer data-app-footer style='height:32px'><span>Disclaimer</span></footer></body></html>"
    )
    hidden = """() => Object.fromEntries(['.app-header', '.settings-pane', '#result-pane', '[data-app-footer] span']
      .map((selector) => [selector, getComputedStyle(document.querySelector(selector)).visibility === 'hidden']))"""
    # The pane's Generate bar is anchored to the footer: it must stay drawn, at its place (OV-3395 review).
    bar = """() => { const bar = document.querySelector('.generate-bar'); const rect = bar.getBoundingClientRect();
      return { shown: bar.checkVisibility({ visibilityProperty: true }) && document.elementFromPoint(rect.x + 5, rect.y + 5) === bar,
               top: rect.top }; }"""
    with playwright_api.sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        try:
            page = browser.new_page(viewport={"width": 800, "height": 600})
            page.set_content(html)
            before = page.evaluate("() => document.querySelector('#result-pane').getBoundingClientRect().toJSON()")
            bar_before = page.evaluate(bar)
            capture_tool.hide_settings_pane_under_popup(page)
            no_popup = page.evaluate(hidden)
            capture_tool.hide_page_around_settings_pane(page)
            open_select = page.evaluate(hidden)
            after = page.evaluate("() => document.querySelector('#result-pane').getBoundingClientRect().toJSON()")
            bar_after = page.evaluate(bar)
            page.set_content(html.replace("</main>", "<div class='feature-popup'>dnaA</div></main>"))
            capture_tool.hide_settings_pane_under_popup(page)
            popup = page.evaluate(hidden)
        finally:
            browser.close()

    assert not any(no_popup.values())
    assert open_select == {".app-header": True, ".settings-pane": False, "#result-pane": True, "[data-app-footer] span": True}
    assert before == after
    assert bar_before["shown"] and bar_after == bar_before
    assert popup == {".app-header": False, ".settings-pane": True, "#result-pane": False, "[data-app-footer] span": False}


def test_crop_edges_measure_the_text_a_reader_sees_and_the_pointer_leaves_no_hover() -> None:
    """OV-3410: the text a crop edge must not cut is the text drawn on the page: a line of HTML text,
    an SVG label (not the clickable figure group around it) and a control count; text under a
    popup, hidden text and text scrolled out of its pane do not, and a label half under a popup
    counts only where it shows. OV-3411: a parked pointer leaves no hover from the last click."""

    playwright_api = pytest.importorskip("playwright.sync_api")
    html = (
        "<html><body style='margin:0;font:16px sans-serif'>"
        "<style>button:hover { outline: 3px dashed blue; }</style>"
        "<h3 id='heading' style='position:absolute;left:20px;top:40px;margin:0'>DEFAULT COLORS</h3>"
        "<svg style='position:absolute;left:300px;top:40px' width='200' height='60'>"
        "<g role='button'><text id='label' x='10' y='30' font-size='20'>Homo sapiens</text></g></svg>"
        "<button id='skip' style='position:absolute;left:20px;top:120px'>Skip this record</button>"
        "<p style='position:absolute;left:20px;top:200px;margin:0;visibility:hidden'>hidden text</p>"
        "<p style='position:absolute;left:300px;top:200px;margin:0'>under the popup</p>"
        "<div class='feature-popup' style='position:absolute;left:280px;top:180px;width:200px;height:60px;"
        "background:#fff;z-index:5'></div>"
        "<p style='position:absolute;left:300px;top:400px;margin:0;font-size:20px'>half under the popup</p>"
        "<div style='position:absolute;left:280px;top:380px;width:100px;height:30px;background:#fff;z-index:5'></div>"
        "<div style='position:absolute;left:20px;top:300px;width:200px;height:40px;overflow:hidden'>"
        "<div style='height:200px'></div><p style='margin:0'>scrolled out</p></div>"
        "</body></html>"
    )
    # The hover style a frame later (`matches(':hover')` answers false in headless Chromium).
    hover_outline = """() => new Promise((resolve) => requestAnimationFrame(() => requestAnimationFrame(
      () => resolve(getComputedStyle(document.querySelector('#skip')).outlineStyle === 'dashed'))))"""
    with playwright_api.sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        try:
            page = browser.new_page(viewport={"width": 800, "height": 600})
            page.set_content(html)
            boxes = capture_tool.visible_text_boxes(page, {"x": 0.0, "y": 0.0, "width": 800.0, "height": 600.0})
            page.locator("#skip").click()
            hovered_after_click = page.evaluate(hover_outline)
            capture_tool.park_pointer(page)
            hovered_after_park = page.evaluate(hover_outline)
        finally:
            browser.close()

    texts = sorted({str(box["text"]) for box in boxes})
    assert texts == [
        "DEFAULT COLORS", "Homo sapiens", "Skip this record", "button: Skip this record", "half under the popup",
    ]
    # The half-covered label is an L: an edge through its covered corner cuts nothing a reader sees,
    # one through the part below the popup cuts it.
    half = [box for box in boxes if box["text"] == "half under the popup"]
    page_size = {"width": 800.0, "height": 600.0}
    assert not any(box["x"] <= 330 <= box["x"] + box["width"] and box["y"] <= 403 <= box["y"] + box["height"] for box in half)
    def top_cuts(top: float) -> list[str]:
        clip = {"x": 290.0, "y": top, "width": 80.0, "height": 100.0}
        return [side for side, _ in capture_tool.crop_edge_cuts(clip, half, page_size) if side == "top"]

    assert top_cuts(404.0) == []
    assert top_cuts(416.0)
    heading = next(box for box in boxes if box["text"] == "DEFAULT COLORS")
    assert 36 <= heading["y"] <= 42 and heading["height"] < 30
    assert hovered_after_click and not hovered_after_park


def test_a_highlight_frame_covers_no_text_around_its_element() -> None:
    """The frame around a capture's highlightSelector stays clear of a label 2 px above the
    highlighted field (BGC 11-1: it covered the descenders of Head Length Ratio), and keeps its
    full size around a field with nothing near it, also around a field in a painted pill: the
    field's own control, its value and the pill that draws it are not foreign text (the measure
    leaves out the highlighted element, `exclude`; L2)."""

    playwright_api = pytest.importorskip("playwright.sync_api")
    html = (
        "<html><body style='margin:0;font:16px sans-serif'>"
        "<label style='position:absolute;left:20px;top:20px;line-height:20px'>Head Length Ratio</label>"
        "<input aria-label='near' style='position:absolute;left:20px;top:42px;width:120px;height:24px'>"
        "<input aria-label='alone' value='0.4' style='position:absolute;left:300px;top:200px;width:120px;height:24px'>"
        "<span style='position:absolute;left:300px;top:300px;padding:3px;border:1px solid #cbd5e1;"
        "border-radius:999px;background:#f8fafc'><input aria-label='wrapped' value='0.6' style='width:80px;height:20px'></span>"
        "</body></html>"
    )
    frame = """(label) => [...document.querySelectorAll('input')].map((input) => {
      const style = getComputedStyle(input);
      const rect = input.getBoundingClientRect();
      const reach = parseFloat(style.outlineOffset) + parseFloat(style.outlineWidth);
      return {label: input.getAttribute('aria-label'), top: rect.top - reach, reach,
              offset: parseFloat(style.outlineOffset), glow: style.boxShadow};
    })"""
    with playwright_api.sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        try:
            page = browser.new_page(viewport={"width": 600, "height": 400})
            page.set_content(html)
            capture_tool.draw_highlight_frames(page, "input")
            capture_tool.assert_highlight_frames_cover_no_text(page, "input")
            frames = {item["label"]: item for item in page.evaluate(frame, None)}
            label_bottom = page.evaluate(
                """() => { const range = document.createRange();
                  range.selectNodeContents(document.querySelector('label').firstChild);
                  return range.getBoundingClientRect().bottom; }"""
            )
        finally:
            browser.close()

    assert frames["near"]["top"] >= label_bottom + capture_tool.HIGHLIGHT_FRAME_CLEARANCE - 0.01
    assert frames["alone"]["reach"] == capture_tool.HIGHLIGHT_FRAME_OFFSET + capture_tool.HIGHLIGHT_FRAME_WIDTH
    for label in ("alone", "wrapped"):
        assert frames[label]["offset"] == capture_tool.HIGHLIGHT_FRAME_OFFSET, label
        assert frames[label]["glow"].endswith(" 0px 0px 0px 4px"), label


def _boxes_of(page, region):
    return capture_tool.visible_text_boxes(page, region)


def test_text_counts_as_covered_only_where_something_paints_over_it() -> None:
    """A layer without pointer events that paints (a capture menu) hides the text under it; a clear layer over the
    page hides nothing, with or without pointer events; a figure label without pointer events
    over a painted background is seen (M1)."""

    playwright_api = pytest.importorskip("playwright.sync_api")
    html = (
        "<html><body style='margin:0;font:16px sans-serif'>"
        "<div style='position:absolute;left:0;top:0;width:600px;height:400px;background:#eef'></div>"
        "<span style='position:absolute;left:20px;top:20px'>under the menu</span>"
        "<div class='gbdraw-capture-select-menu' style='position:absolute;left:10px;top:10px;width:200px;"
        "height:40px;background:#fff;pointer-events:none;z-index:5'></div>"
        "<span style='position:absolute;left:20px;top:120px'>under a clear layer</span>"
        "<span style='position:absolute;left:300px;top:120px'>under a clear hit layer</span>"
        "<span style='position:absolute;left:20px;top:220px;pointer-events:none'>figure label</span>"
        "<div style='position:fixed;inset:0;pointer-events:none;z-index:6'></div>"
        "<div style='position:absolute;left:290px;top:110px;width:250px;height:40px;z-index:7'></div>"
        "</body></html>"
    )
    with playwright_api.sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        try:
            page = browser.new_page(viewport={"width": 600, "height": 400})
            page.set_content(html)
            boxes = _boxes_of(page, {"x": 0.0, "y": 0.0, "width": 600.0, "height": 400.0})
            whole = page.evaluate(
                """() => Object.fromEntries([...document.querySelectorAll('span')].map((span) => {
                  const range = document.createRange();
                  range.selectNodeContents(span);
                  const rect = range.getBoundingClientRect();
                  return [span.textContent, rect.width * rect.height];
                }))"""
            )
        finally:
            browser.close()

    seen: dict[str, float] = {}
    for box in boxes:
        seen[box["text"]] = max(seen.get(box["text"], 0.0), box["width"] * box["height"])
    assert "under the menu" not in seen
    for text in ("under a clear layer", "under a clear hit layer", "figure label"):
        assert seen[text] >= whole[text] - 0.5, text


def test_a_control_is_measured_with_the_pill_its_wrapper_paints() -> None:
    """A help button inside a painted pill (a Palette clock note) is measured with the pill, so
    an edge through the pill below the button cuts it (M2)."""

    playwright_api = pytest.importorskip("playwright.sync_api")
    html = (
        "<html><body style='margin:0;font:12px sans-serif'>"
        "<span class='apply-note' style='position:absolute;left:40px;top:40px;display:inline-flex;"
        "padding:3px 6px;border:1px solid #cbd5e1;border-radius:999px;background:#f8fafc'>"
        "<button style='all:unset;width:12px;height:12px'>?</button></span>"
        "</body></html>"
    )
    with playwright_api.sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        try:
            page = browser.new_page(viewport={"width": 300, "height": 200})
            page.set_content(html)
            boxes = _boxes_of(page, {"x": 0.0, "y": 0.0, "width": 300.0, "height": 200.0})
            pill = page.evaluate(
                """() => { const r = document.querySelector('.apply-note').getBoundingClientRect();
                  const b = document.querySelector('button').getBoundingClientRect();
                  return {top: r.top, bottom: r.bottom, buttonBottom: b.bottom}; }"""
            )
        finally:
            browser.close()

    control = [box for box in boxes if box["text"].startswith("button")]
    assert control and min(box["y"] for box in control) <= pill["top"] + 0.01
    assert max(box["y"] + box["height"] for box in control) >= pill["bottom"] - 0.01
    edge = (pill["buttonBottom"] + pill["bottom"]) / 2
    clip = {"x": 0.0, "y": 0.0, "width": 300.0, "height": edge}
    assert capture_tool.crop_edge_cuts(clip, boxes, {"width": 300.0, "height": 200.0})


def test_the_seen_part_of_a_long_line_starts_where_its_cover_ends() -> None:
    """A 1200 px line covered from x 0 to 630: an edge at x 640 crosses its seen part (L1)."""

    playwright_api = pytest.importorskip("playwright.sync_api")
    html = (
        "<html><body style='margin:0;font:20px monospace'>"
        f"<span style='position:absolute;left:0;top:100px;white-space:pre'>{'m' * 100}</span>"
        "<div style='position:absolute;left:0;top:80px;width:630px;height:60px;background:#fff'></div>"
        "</body></html>"
    )
    with playwright_api.sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        try:
            page = browser.new_page(viewport={"width": 1400, "height": 300})
            page.set_content(html)
            width = page.evaluate("() => document.querySelector('span').getBoundingClientRect().width")
            boxes = _boxes_of(page, {"x": 0.0, "y": 0.0, "width": 1400.0, "height": 300.0})
        finally:
            browser.close()

    assert abs(width - 1200.0) < 30
    clip = {"x": 560.0, "y": 60.0, "width": 80.0, "height": 120.0}
    cuts = capture_tool.crop_edge_cuts(clip, boxes, {"width": 1400.0, "height": 300.0})
    assert [side for side, _ in cuts] == ["right"]


def test_text_its_own_box_clips_away_is_not_seen() -> None:
    """A screen-reader line (`.sr-only`: a 1 px box, `overflow:hidden`, `clip: rect(0 0 0 0)`) at the
    end of a card head runs on into the pane's padding but is never drawn: an edge in that padding
    cuts nothing (BGC 7-2, 7-3: the right edge moved out to the panel's border). A line whose
    right half a clip-path hides is seen only where it is drawn."""

    playwright_api = pytest.importorskip("playwright.sync_api")
    html = (
        "<html><body style='margin:0;font:14px sans-serif'>"
        "<style>.sr-only { position:absolute; width:1px; height:1px; padding:0; margin:-1px; overflow:hidden;"
        " clip:rect(0,0,0,0); white-space:nowrap; border-width:0; }</style>"
        "<div class='pane' style='position:absolute;left:20px;top:20px;width:400px;padding:16px;overflow-x:hidden;"
        "background:#fff;border:1px solid #cbd5e1'>"
        "<div style='display:flex;justify-content:space-between'><span>Arrow direction</span>"
        "<span style='display:flex;justify-content:flex-end'>"
        "<span class='sr-only'>Arrow direction for BGC0000713, drawn reversed</span></span></div>"
        "</div><span class='sr-only' style='left:480px;top:120px'>a line that runs out of the viewport</span>"
        "<span id='half' style='position:absolute;left:20px;top:160px;white-space:pre;"
        "clip-path:inset(0 50% 0 0)'>mmmmmmmmmmmmmmmmmmmm</span>"
        "</body></html>"
    )
    with playwright_api.sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        try:
            page = browser.new_page(viewport={"width": 500, "height": 200})
            page.set_content(html)
            boxes = _boxes_of(page, {"x": 0.0, "y": 0.0, "width": 500.0, "height": 200.0})
            half = page.evaluate("() => document.getElementById('half').getBoundingClientRect().toJSON()")
            content_right = page.evaluate(
                """() => { const pane = document.querySelector('.pane').getBoundingClientRect();
                  return pane.right - 1 - 16; }"""
            )
        finally:
            browser.close()

    assert sorted({box["text"] for box in boxes}) == ["Arrow direction", "mmmmmmmmmmmmmmmmmmmm"]
    drawn = [box for box in boxes if box["text"].startswith("m")]
    middle = half["x"] + half["width"] / 2
    assert max(box["x"] + box["width"] for box in drawn) <= middle + 2.5
    boxes = [box for box in boxes if box not in drawn]
    clip = {"x": 10.0, "y": 10.0, "width": content_right + 8.0 - 10.0, "height": 60.0}
    page_size = {"width": 500.0, "height": 200.0}
    padding = capture_tool.normalize_crop_padding(12)
    assert capture_tool.clear_crop_edges(clip, boxes, padding, page_size) == clip


def test_an_opaque_layer_without_pointer_events_covers_the_text_under_it() -> None:
    """A help tip, a hover summary or a readout the page draws without pointer events hides the
    text under it, by what it paints and not by its name (L1); a clear one hides nothing."""

    playwright_api = pytest.importorskip("playwright.sync_api")
    html = (
        "<html><body style='margin:0;font:16px sans-serif'>"
        "<span style='position:absolute;left:20px;top:20px'>under the tip</span>"
        "<div role='tooltip' style='position:fixed;left:10px;top:10px;width:220px;height:40px;"
        "background:#1e293b;color:#fff;pointer-events:none;z-index:50'><span>Help</span></div>"
        "<span style='position:absolute;left:20px;top:120px'>under a clear tip</span>"
        "<div style='position:fixed;left:10px;top:110px;width:220px;height:40px;pointer-events:none;z-index:50'></div>"
        "</body></html>"
    )
    with playwright_api.sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        try:
            page = browser.new_page(viewport={"width": 400, "height": 200})
            page.set_content(html)
            boxes = _boxes_of(page, {"x": 0.0, "y": 0.0, "width": 400.0, "height": 200.0})
            hittable = page.evaluate("() => document.elementFromPoint(200, 45).tagName")
        finally:
            browser.close()

    texts = {box["text"] for box in boxes}
    assert "under the tip" not in texts
    assert {"Help", "under a clear tip"} <= texts
    # The measure leaves the page's pointer events as they were.
    assert hittable == "BODY"


_FIGURE_WITH_LABELS = (
    "<html><body style='margin:0;font:14px sans-serif;background:#fff'>"
    "<svg width='400' height='240' style='position:absolute;left:0;top:0'>"
    "<path id='arc' d='M 40 40 Q 200 160 360 40' fill='none' stroke='#2563eb' stroke-width='12'/>"
    "<text id='inside' x='180' y='70' font-size='14'>dnaB</text>"
    "<text id='across' x='260' y='100' font-size='14'>dnaC</text>"
    "<path id='alone' d='M 40 200 L 120 200' stroke='#2563eb' stroke-width='12'/>"
    "</svg></body></html>"
)
# Moves the label `across` so that its text box ends `by` px below the box the arc's outline is
# drawn around (Chromium draws an SVG element's outline around its stroke box, not its
# getBoundingClientRect).
_PLACE_ACROSS = """([by, edge]) => {
  const label = document.getElementById('across');
  const range = document.createRange();
  range.selectNodeContents(label.firstChild);
  const box = range.getBoundingClientRect();
  label.setAttribute('y', String(parseFloat(label.getAttribute('y')) + edge + by - box.bottom));
}"""


def _arc_outline_box(page) -> dict[str, float]:
    handle = page.query_selector("#arc")
    element = handle.evaluate(capture_tool.HIGHLIGHT_ELEMENT_BOX_SCRIPT)
    return capture_tool.highlight_outline_box(page, handle, element, True)


def _screenshot(page) -> Image.Image:
    with Image.open(io.BytesIO(page.screenshot(type="png"))) as shot:
        return shot.convert("RGB")


def _changed_pixels_in(before: Image.Image, after: Image.Image, boxes) -> int:
    """The pixels in ``boxes`` (CSS px, device scale 1) that differ between two screenshots."""

    count = 0
    for box in boxes:
        for x in range(int(box["x"]), int(box["x"] + box["width"] + 0.999)):
            for y in range(int(box["y"]), int(box["y"] + box["height"] + 0.999)):
                count += before.getpixel((x, y)) != after.getpixel((x, y))
    return count


def test_a_highlight_frame_keeps_off_a_label_inside_or_across_its_figure_path() -> None:
    """OV-3492: a label wholly inside a highlighted figure path's box (ATskew 10-1's l-rRNA) leaves
    the frame its full size; a label 1.5 px across the edge of the box the outline is drawn around
    is passed by the outline at its usual place, and no pixel of the drawn frame lies in either
    label's box. An SVG element draws no glow."""

    playwright_api = pytest.importorskip("playwright.sync_api")
    frames = """() => Object.fromEntries(['arc', 'alone'].map((id) => {
      const style = getComputedStyle(document.getElementById(id));
      return [id, { offset: parseFloat(style.outlineOffset), glow: style.boxShadow }];
    }))"""
    with playwright_api.sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        try:
            page = browser.new_page(viewport={"width": 400, "height": 240})
            page.set_content(_FIGURE_WITH_LABELS)
            outline_box = _arc_outline_box(page)
            page.evaluate(_PLACE_ACROSS, [1.5, outline_box["y"] + outline_box["height"]])
            boxes = _boxes_of(page, {"x": 0.0, "y": 0.0, "width": 400.0, "height": 240.0})
            unframed = _screenshot(page)
            capture_tool.draw_highlight_frames(page, "#arc, #alone")
            drawn = page.evaluate(frames)
            capture_tool.assert_highlight_frames_cover_no_text(page, "#arc, #alone")
            framed = _screenshot(page)
        finally:
            browser.close()

    painted = _changed_pixels_in(unframed, framed, boxes)
    # The outline is drawn: 3 to 6 px outside the left side of the arc's outline box.
    beside_arc = _changed_pixels_in(unframed, framed, [{
        "x": outline_box["x"] - 6.0, "y": outline_box["y"], "width": 3.0, "height": outline_box["height"],
    }])
    assert sorted(box["text"] for box in boxes) == ["dnaB", "dnaC"]
    assert drawn["arc"] == {"offset": capture_tool.HIGHLIGHT_FRAME_OFFSET, "glow": "rgba(245, 158, 11, 0.18) 0px 0px 0px 0px"}
    assert drawn["alone"]["offset"] == capture_tool.HIGHLIGHT_FRAME_OFFSET
    assert beside_arc > 0
    assert painted == 0


def test_no_highlight_frame_is_drawn_across_a_label_it_cannot_pass() -> None:
    """OV-3492: a label 8 px across the edge of a highlighted figure path's outline box lies under
    every outline position within the frame's reach: the capture fails and names the label. A frame
    drawn across a label by other means fails the guard that runs on the captured page."""

    playwright_api = pytest.importorskip("playwright.sync_api")
    with playwright_api.sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        try:
            page = browser.new_page(viewport={"width": 400, "height": 240})
            page.set_content(_FIGURE_WITH_LABELS)
            outline_box = _arc_outline_box(page)
            edge = outline_box["y"] + outline_box["height"]
            page.evaluate(_PLACE_ACROSS, [8.0, edge])
            with pytest.raises(RuntimeError, match="'dnaC'"):
                capture_tool.draw_highlight_frames(page, "#arc")
            page.evaluate(_PLACE_ACROSS, [1.5, edge])
            page.evaluate(
                "(css) => { document.getElementById('arc').style.cssText += `; ${css}`; }",
                capture_tool.highlight_frame_css(-3.0, 0.0),
            )
            with pytest.raises(RuntimeError, match="covers text: 'dnaC'"):
                capture_tool.assert_highlight_frames_cover_no_text(page, "#arc")
        finally:
            browser.close()


def test_text_under_a_translucent_or_border_only_layer_is_seen() -> None:
    """L1: a layer without pointer events hides text only where its paint is opaque. Text under a
    selection marquee (a dashed border and a 12 % tint), under an opaque layer drawn at opacity
    0.5, or inside a border-only frame is seen whole; text under an opaque layer is not."""

    playwright_api = pytest.importorskip("playwright.sync_api")
    layer = "position:absolute;pointer-events:none;z-index:5;"
    html = (
        "<html><body style='margin:0;font:16px sans-serif'>"
        "<span style='position:absolute;left:20px;top:20px'>under the marquee</span>"
        f"<div class='feature-selection-marquee' style='{layer}left:10px;top:10px;width:300px;height:40px;"
        "border:1px dashed #2563eb;background:rgba(37,99,235,0.12)'></div>"
        "<span style='position:absolute;left:20px;top:120px'>under a faded layer</span>"
        f"<div style='{layer}left:10px;top:110px;width:300px;height:40px;background:#fff;opacity:0.5'></div>"
        "<span style='position:absolute;left:20px;top:220px'>inside a border</span>"
        f"<div style='{layer}left:10px;top:210px;width:300px;height:40px;border:2px solid #0f172a'></div>"
        "<span style='position:absolute;left:400px;top:20px'>under a sheet</span>"
        f"<div style='{layer}left:390px;top:10px;width:180px;height:40px;background:#fff'></div>"
        "</body></html>"
    )
    with playwright_api.sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        try:
            page = browser.new_page(viewport={"width": 600, "height": 300})
            page.set_content(html)
            boxes = _boxes_of(page, {"x": 0.0, "y": 0.0, "width": 600.0, "height": 300.0})
            whole = page.evaluate(
                """() => Object.fromEntries([...document.querySelectorAll('span')].map((span) => {
                  const range = document.createRange();
                  range.selectNodeContents(span);
                  const rect = range.getBoundingClientRect();
                  return [span.textContent, rect.width * rect.height];
                }))"""
            )
        finally:
            browser.close()

    seen: dict[str, float] = {}
    for box in boxes:
        seen[box["text"]] = max(seen.get(box["text"], 0.0), box["width"] * box["height"])
    assert "under a sheet" not in seen
    for text in ("under the marquee", "under a faded layer", "inside a border"):
        assert seen.get(text, 0.0) >= whole[text] - 0.5, text


def test_svg_text_in_a_tspan_or_on_a_text_path_is_seen() -> None:
    """Figure text drawn in a `<tspan>` (a definition line) or on a `<textPath>` (a tick label) is
    measured where it is drawn, like any other label, also when the figure sets no pointer events."""

    playwright_api = pytest.importorskip("playwright.sync_api")
    html = (
        "<html><body style='margin:0;font:14px sans-serif'>"
        "<svg width='400' height='300' style='position:absolute;left:0;top:0'>"
        "<defs><path id='ring' d='M 60 200 A 140 140 0 0 1 340 200'/></defs>"
        "<g style='pointer-events:none'>"
        "<text x='120' y='40' font-size='16'><tspan x='120' dy='0'>Homo</tspan>"
        "<tspan x='120' dy='20' font-style='italic'>sapiens</tspan></text>"
        "<text font-size='14'><textPath href='#ring' startOffset='40%'>16 kbp</textPath></text>"
        "</g></svg></body></html>"
    )
    with playwright_api.sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        try:
            page = browser.new_page(viewport={"width": 400, "height": 300})
            page.set_content(html)
            boxes = _boxes_of(page, {"x": 0.0, "y": 0.0, "width": 400.0, "height": 300.0})
            whole = page.evaluate(
                """() => Object.fromEntries([...document.querySelectorAll('tspan, textPath')].map((node) => {
                  const range = document.createRange();
                  range.selectNodeContents(node);
                  return [node.textContent, [...range.getClientRects()].reduce(
                    (area, rect) => area + rect.width * rect.height, 0)];
                }))"""
            )
        finally:
            browser.close()

    seen: dict[str, float] = {}
    for box in boxes:
        seen[box["text"]] = seen.get(box["text"], 0.0) + box["width"] * box["height"]
    assert sorted(seen) == ["16 kbp", "Homo", "sapiens"]
    for text, area in whole.items():
        assert seen[text] >= area - 0.5, text


def test_a_highlight_frame_never_moves_in_over_the_elements_own_caption() -> None:
    """L2: a highlighted row's own caption, 1.5 px inside its top edge, is text a frame must keep off
    like any other: with a neighbour's label 1 px above the row, no outline position is clear, so
    the capture fails and names the caption; a frame forced inside the row fails the guard."""

    playwright_api = pytest.importorskip("playwright.sync_api")
    html = (
        "<html><body style='margin:0;font:14px sans-serif'>"
        "<span id='above' style='position:absolute;left:20px;top:0;line-height:20px'>Records</span>"
        "<div id='row' style='position:absolute;left:20px;top:40px;width:300px;height:40px;background:#f8fafc'>"
        "<span id='caption' style='position:absolute;left:8px;top:6px'>NC_000913 chromosome</span></div>"
        "</body></html>"
    )
    # The label's text box ends 1 px above the row; the caption's text box starts 1.5 px inside it.
    place = """() => {
      const box = (id) => { const range = document.createRange();
        range.selectNodeContents(document.getElementById(id).firstChild); return range.getBoundingClientRect(); };
      const row = document.getElementById('row');
      row.style.top = `${box('above').bottom + 1}px`;
      const caption = document.getElementById('caption');
      caption.style.top = `${parseFloat(caption.style.top) + row.getBoundingClientRect().top + 1.5 - box('caption').top}px`;
    }"""
    with playwright_api.sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        try:
            page = browser.new_page(viewport={"width": 400, "height": 200})
            page.set_content(html)
            page.evaluate(place)
            with pytest.raises(RuntimeError, match="'NC_000913 chromosome'"):
                capture_tool.draw_highlight_frames(page, "#row")
            page.evaluate(
                "(css) => { document.getElementById('row').style.cssText += `; ${css}`; }",
                capture_tool.highlight_frame_css(-3.0, 0.0),
            )
            with pytest.raises(RuntimeError, match="covers text: 'NC_000913 chromosome'"):
                capture_tool.assert_highlight_frames_cover_no_text(page, "#row")
        finally:
            browser.close()


def test_the_frame_guard_sees_a_frame_on_an_element_that_eases_its_style() -> None:
    """Review M1: a field that eases every property (`transition: all`, here after a long delay)
    still shows the guard the frame it paints, because the guard draws it clear and in colour with
    transitions off: a frame forced over the label above it fails the guard. A guard that let the
    field ease would take both shots inside the delay, see no difference and pass it."""

    playwright_api = pytest.importorskip("playwright.sync_api")
    html = (
        "<html><body style='margin:0;font:16px sans-serif'>"
        "<label style='position:absolute;left:20px;top:20px;line-height:20px'>Head Length Ratio</label>"
        "<input aria-label='ratio' value='0.64' style='position:absolute;left:20px;top:42px;width:120px;"
        "height:24px'>"
        "</body></html>"
    )
    with playwright_api.sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        try:
            page = browser.new_page(viewport={"width": 400, "height": 200})
            page.set_content(html)
            page.evaluate(
                "(css) => { document.querySelector('input').style.cssText += `; ${css}`; }",
                capture_tool.highlight_frame_css(3.0, 4.0),
            )
            capture_tool.wait_for_css_transitions(page)
            # The field eases its style from now on, after a delay far longer than two screenshots.
            page.evaluate("() => { document.querySelector('input').style.transition = 'all 10s 10s'; }")
            with pytest.raises(RuntimeError, match="covers text: 'Head Length Ratio'"):
                capture_tool.assert_highlight_frames_cover_no_text(page, "input")
            # The guard leaves the frame drawn and starts no transition.
            settled = page.evaluate(
                """() => { const input = document.querySelector('input');
                  return [getComputedStyle(input).outlineColor, input.getAnimations().length]; }"""
            )
        finally:
            browser.close()

    assert settled == ["rgb(245, 158, 11)", 0]


def test_a_capture_scroll_after_the_record_read_is_not_undone() -> None:
    """HmmtDNA basic 3-1 (OV-3547): a record read that opens the single-record panel puts the
    Settings pane back where it was once the panel has painted, and the app keeps that restore as
    a promise (`waitForCircularRecordScrollAnchor`). The read awaits it, so a scroll the capture
    makes right after it stays. Here the restore waits for the test to release it (or for 2 s, so
    a tool that awaits it goes on): a tool that does not await it scrolls first, and the release
    then undoes the scroll on every run."""

    playwright_api = pytest.importorskip("playwright.sync_api")
    html = (
        "<html><body style='margin:0'>"
        "<div class='settings-scroll' style='height:200px;overflow:auto'>"
        "<div style='height:1000px'></div><div id='preset'>Track Preset</div><div style='height:1000px'></div>"
        "</div><script>"
        "const discovery = { hasInput: true, status: 'pending' };"
        "let anchor = Promise.resolve();"
        "const released = new Promise((resolve) => { window.__releaseRestore = resolve; });"
        "window.__GBDRAW_APP__ = { mode: 'circular', circularRecordDiscoveryState: discovery,"
        "  async inspectCircularSourceRecords() {"
        "    const pane = document.querySelector('.settings-scroll'); const saved = pane.scrollTop;"
        "    setTimeout(() => { discovery.status = 'ready';"
        "      anchor = Promise.race([released, new Promise((resolve) => setTimeout(resolve, 2000))])"
        "        .then(() => { pane.scrollTop = saved; });"
        "    }, 30);"
        "  },"
        "  waitForCircularRecordScrollAnchor: () => anchor };"
        "</script></body></html>"
    )
    with playwright_api.sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        try:
            page = browser.new_page(viewport={"width": 400, "height": 300})
            page.set_content(html)
            capture_tool.inspect_shown_records(page)
            page.evaluate("() => document.getElementById('preset').scrollIntoView({ block: 'center' })")
            page.evaluate("async () => { window.__releaseRestore(); await window.__GBDRAW_APP__.waitForCircularRecordScrollAnchor(); }")
            scroll_top = page.evaluate("() => document.querySelector('.settings-scroll').scrollTop")
        finally:
            browser.close()

    assert scroll_top > 800


def test_a_highlighted_pill_keeps_its_frame_inside_off_a_label_above() -> None:
    """Review D4: a highlighted pill that draws a single field is the shape its frame goes around,
    so the field inside it (whose shape the pill draws) is not text that keeps the frame out: with
    a label 1 px above the pill, the frame moves in, inside the pill and clear of the field's value."""

    playwright_api = pytest.importorskip("playwright.sync_api")
    html = (
        "<html><body style='margin:0;font:14px sans-serif'>"
        "<span id='above' style='position:absolute;left:20px;top:0;line-height:20px'>Records</span>"
        "<div class='pill' style='position:absolute;left:20px;top:40px;padding:4px;border:1px solid #cbd5e1;"
        "background:#f8fafc'><input aria-label='width' value='12' style='display:block;width:120px;height:24px;"
        "margin:0;padding:1px 4px;border:1px solid #94a3b8;font-size:14px'></div>"
        "</body></html>"
    )
    place = """() => { const range = document.createRange();
      range.selectNodeContents(document.getElementById('above').firstChild);
      document.querySelector('.pill').style.top = `${range.getBoundingClientRect().bottom + 1}px`; }"""
    with playwright_api.sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        try:
            page = browser.new_page(viewport={"width": 400, "height": 200})
            page.set_content(html)
            page.evaluate(place)
            capture_tool.draw_highlight_frames(page, "div.pill")
            capture_tool.wait_for_css_transitions(page)
            capture_tool.assert_highlight_frames_cover_no_text(page, "div.pill")
            offset = page.evaluate("() => parseFloat(getComputedStyle(document.querySelector('.pill')).outlineOffset)")
        finally:
            browser.close()

    assert offset < 0


def test_a_highlight_frame_keeps_off_the_value_its_field_shows() -> None:
    """Review D5: the value a highlighted field shows is its own text: its content box stays inside
    the frame and the guard checks it. With a label 1 px above a field whose value sits 2 px inside
    its edge, no outline position is clear, so the capture fails naming the field; a frame forced
    inside the field fails the guard."""

    playwright_api = pytest.importorskip("playwright.sync_api")
    html = (
        "<html><body style='margin:0;font:14px sans-serif'>"
        "<span id='above' style='position:absolute;left:20px;top:0;line-height:20px'>Records</span>"
        "<input aria-label='ratio' value='0.64' style='position:absolute;left:20px;top:40px;width:120px;"
        "height:20px;box-sizing:border-box;margin:0;padding:1px 4px;border:1px solid #94a3b8;"
        "font-size:14px;line-height:16px'>"
        "</body></html>"
    )
    place = """() => { const range = document.createRange();
      range.selectNodeContents(document.getElementById('above').firstChild);
      document.querySelector('input').style.top = `${range.getBoundingClientRect().bottom + 1}px`; }"""
    with playwright_api.sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        try:
            page = browser.new_page(viewport={"width": 400, "height": 200})
            page.set_content(html)
            page.evaluate(place)
            with pytest.raises(RuntimeError, match="'input: ratio'"):
                capture_tool.draw_highlight_frames(page, "input")
            page.evaluate(
                "(css) => { document.querySelector('input').style.cssText += `; ${css}`; }",
                capture_tool.highlight_frame_css(-3.0, 0.0),
            )
            with pytest.raises(RuntimeError, match="covers text: 'input: ratio'"):
                capture_tool.assert_highlight_frames_cover_no_text(page, "input")
        finally:
            browser.close()


def test_text_under_an_opaque_drawer_is_not_seen_past_its_edge() -> None:
    """HmmtDNA_ATskew post 2-1: figure labels run under the Editor drawer (opaque, a 1 px border).
    The hit test snaps a point less than half a pixel left of the drawer onto it, and the drawer
    hides that point: the part a reader sees ends at the drawer's edge, so a crop edge 1 px inside
    the drawer cuts nothing."""

    playwright_api = pytest.importorskip("playwright.sync_api")
    html = (
        "<html><body style='margin:0;font:13px sans-serif'>"
        "<svg width='400' height='200' style='position:absolute;left:0;top:0'>"
        "<text x='40' y='60' font-size='13' style='pointer-events:none'>tRNA-Phe tRNA-Phe tRNA-Phe</text></svg>"
        "<div id='drawer' style='position:absolute;left:0;top:0;width:280px;height:200px;background:#fff;"
        "border-left:1px solid #e2e8f0;z-index:30'></div>"
        "</body></html>"
    )
    # The drawer starts where the measure samples the label 0.5-1 px short of a whole pixel.
    place = """() => {
      const range = document.createRange();
      range.selectNodeContents(document.querySelector('text').firstChild);
      const rect = range.getBoundingClientRect();
      const columns = Math.min(200, Math.max(3, Math.ceil(rect.width / 2)));
      for (let k = Math.floor(columns / 2); k < columns; k += 1) {
        const x = rect.left + (k / columns) * rect.width;
        if (x - Math.floor(x) >= 0.55 && x - Math.floor(x) <= 0.9) {
          document.getElementById('drawer').style.left = `${Math.ceil(x)}px`;
          return Math.ceil(x);
        }
      }
      throw new Error('no sample point to place the drawer at');
    }"""
    with playwright_api.sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        try:
            page = browser.new_page(viewport={"width": 400, "height": 200})
            page.set_content(html)
            left = float(page.evaluate(place))
            boxes = _boxes_of(page, {"x": 0.0, "y": 0.0, "width": 400.0, "height": 200.0})
        finally:
            browser.close()

    assert boxes
    clip = {"x": left + 1.0, "y": 0.0, "width": 200.0, "height": 200.0}
    assert capture_tool.crop_edge_cuts(clip, boxes, {"width": 400.0, "height": 200.0}) == []
