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
