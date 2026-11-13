from __future__ import annotations

import pytest

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
