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
                page, {"selectedLabel": "palette_40", "cardTitle": "Colors", "placement": "right"}
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
        "<main><div class='settings-pane' style='width:300px'><select><option>A</option></select></div>"
        "<div id='result-pane'>figure</div></main><footer data-app-footer>Disclaimer</footer></body></html>"
    )
    hidden = """() => Object.fromEntries(['.app-header', '.settings-pane', '#result-pane', '[data-app-footer]']
      .map((selector) => [selector, getComputedStyle(document.querySelector(selector)).visibility === 'hidden']))"""
    with playwright_api.sync_playwright() as playwright:
        browser = playwright.chromium.launch()
        try:
            page = browser.new_page(viewport={"width": 800, "height": 600})
            page.set_content(html)
            before = page.evaluate("() => document.querySelector('#result-pane').getBoundingClientRect().toJSON()")
            capture_tool.hide_settings_pane_under_popup(page)
            no_popup = page.evaluate(hidden)
            capture_tool.hide_page_around_settings_pane(page)
            open_select = page.evaluate(hidden)
            after = page.evaluate("() => document.querySelector('#result-pane').getBoundingClientRect().toJSON()")
            page.set_content(html.replace("</main>", "<div class='feature-popup'>dnaA</div></main>"))
            capture_tool.hide_settings_pane_under_popup(page)
            popup = page.evaluate(hidden)
        finally:
            browser.close()

    assert not any(no_popup.values())
    assert open_select == {".app-header": True, ".settings-pane": False, "#result-pane": True, "[data-app-footer]": True}
    assert before == after
    assert popup == {".app-header": False, ".settings-pane": True, "#result-pane": False, "[data-app-footer]": False}
