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


def test_a_highlight_frame_covers_no_text_around_its_element() -> None:
    """The frame around a capture's highlightSelector stays clear of a label 2 px above the
    highlighted field (BGC 11-1: it covered the descenders of Head Length Ratio), and keeps its
    full size around a field with nothing near it."""

    playwright_api = pytest.importorskip("playwright.sync_api")
    html = (
        "<html><body style='margin:0;font:16px sans-serif'>"
        "<label style='position:absolute;left:20px;top:20px;line-height:20px'>Head Length Ratio</label>"
        "<input aria-label='near' style='position:absolute;left:20px;top:42px;width:120px;height:24px'>"
        "<input aria-label='alone' style='position:absolute;left:300px;top:200px;width:120px;height:24px'>"
        "</body></html>"
    )
    frame = """(label) => [...document.querySelectorAll('input')].map((input) => {
      const style = getComputedStyle(input);
      const rect = input.getBoundingClientRect();
      const reach = parseFloat(style.outlineOffset) + parseFloat(style.outlineWidth);
      return {label: input.getAttribute('aria-label'), top: rect.top - reach, reach,
              glow: style.boxShadow};
    })"""
    with playwright_api.sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        try:
            page = browser.new_page(viewport={"width": 600, "height": 400})
            page.set_content(html)
            capture_tool.draw_highlight_frames(page, "input")
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
    assert " 4px" in frames["alone"]["glow"]


def _boxes_of(page, region):
    return capture_tool.visible_text_boxes(page, region)


def test_text_counts_as_covered_only_where_something_paints_over_it() -> None:
    """A capture menu without pointer events hides the text under it; a clear layer over the
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
