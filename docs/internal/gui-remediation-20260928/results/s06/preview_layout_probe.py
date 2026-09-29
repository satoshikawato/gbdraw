#!/usr/bin/env python
"""S06 Preview layout probe: search, Editor, toolbar, canvas and Align/Review geometry.

Usage: python preview_layout_probe.py <base-url>
Loads the BGC Gallery Session through the real Load Session button, then measures
each viewport with the Editor closed and open, and the feature popup and the
Similarity drawer actions. 200% zoom is emulated as 720x450 CSS px at device scale 2.
Writes $S00_BASELINE_DIR/evidence/s06-layout/probe.json and screenshots.
"""
import json
import os
import pathlib
import sys

from playwright.sync_api import sync_playwright

WT = pathlib.Path(__file__).resolve().parents[5]
OUT = pathlib.Path(os.environ.get('S00_BASELINE_DIR', '/home/kawato/gbdraw-baselines/gui-remediation-20260928')) / 'evidence/s06-layout'
SESSION = WT / 'gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json'
VIEWPORTS = [('1440x900', 1440, 900, 1), ('1024x740', 1024, 740, 1), ('768x740', 768, 740, 1),
             ('390x844', 390, 844, 1), ('390x740', 390, 740, 1), ('zoom200-1440x900', 720, 450, 2)]

MEASURE_JS = r"""() => {
  const box = (el) => { if (!el) return null; const r = el.getBoundingClientRect();
    const style = getComputedStyle(el);
    return { x: +r.x.toFixed(1), y: +r.y.toFixed(1), w: +r.width.toFixed(1), h: +r.height.toFixed(1),
      visible: style.visibility !== 'hidden' && style.display !== 'none' && r.width > 0 && r.height > 0 }; };
  const q = (s) => document.querySelector(s);
  const layout = q('.preview-editor-layout'), drawer = q('.right-drawer');
  return { viewport: { w: innerWidth, h: innerHeight }, layout: box(layout), search: box(q('.preview-feature-search')),
    drawer: drawer && drawer.getAttribute('aria-hidden') === 'false' ? box(drawer) : null,
    toolbar: box(q('.preview-controls')), canvas: box(q('.preview-viewport')),
    toggle: box(q('.drawer-toggle')), close: box(q('.right-drawer [aria-label="Close editor"], .right-drawer button[title="Close"]')),
    pageScrollWidth: document.documentElement.scrollWidth };
}"""


def overlap(a, b):
    if not a or not b:
        return 0.0
    w = min(a['x'] + a['w'], b['x'] + b['w']) - max(a['x'], b['x'])
    h = min(a['y'] + a['h'], b['y'] + b['h']) - max(a['y'], b['y'])
    return round(max(0.0, w) * max(0.0, h), 1)


def checks(m, compact):
    s, d, t, c, v = m['search'], m['drawer'], m['toolbar'], m['canvas'], m['viewport']
    out = {'searchInsideViewport': bool(s) and s['x'] >= -0.5 and s['x'] + s['w'] <= v['w'] + 0.5,
           'searchWidth': s and s['w'], 'canvasHeight': c and c['h'], 'noHorizontalPageScroll': m['pageScrollWidth'] <= v['w'] + 1}
    if d:
        out.update({'searchDrawerOverlap': overlap(s, d), 'toolbarDrawerOverlap': overlap(t, d),
                    'drawerTopMinusLayoutTop': round(d['y'] - m['layout']['y'], 1),
                    'drawerStartsBelowSearch': d['y'] >= s['y'] + s['h'] - 0.5})
    return out


def main(base_url):
    OUT.mkdir(parents=True, exist_ok=True)
    result = {'url': base_url, 'viewports': {}}
    with sync_playwright() as p:
        browser = p.chromium.launch()
        for name, width, height, scale in VIEWPORTS:
            context = browser.new_context(viewport={'width': width, 'height': height}, device_scale_factor=scale)
            page = context.new_page()
            dialogs = []
            page.on('dialog', lambda d: (dialogs.append(d.message), d.accept()))
            page.goto(f'{base_url}/gbdraw/web/index.html', wait_until='domcontentloaded')
            page.wait_for_function('() => window.__GBDRAW_APP__')
            with page.expect_file_chooser() as chooser:
                page.get_by_role('button', name='Load Session').first.click()
            chooser.value.set_files(str(SESSION))
            page.wait_for_function('() => window.__GBDRAW_APP__.results?.length > 0 && !window.__GBDRAW_APP__.sessionImportPending',
                                   timeout=300000)
            page.wait_for_timeout(1500)
            entry = {}
            compact = page.evaluate("() => getComputedStyle(document.querySelector('[aria-label=\"Result Preview\"]')).getPropertyValue('--alignment-review-compact').trim() === '1'")
            entry['compact'] = compact
            page.locator('.preview-editor-layout').scroll_into_view_if_needed()
            closed = page.evaluate(MEASURE_JS)
            entry['closed'] = {'measure': closed, 'checks': checks(closed, compact)}
            page.locator('.preview-editor-layout').screenshot(path=str(OUT / f'{name}-closed.png'))
            page.evaluate("() => { const app = window.__GBDRAW_APP__; app.openRightDrawerTab('orthogroups'); app.selectedOrthogroupId = 'og_18'; }")
            page.wait_for_timeout(800)
            opened = page.evaluate(MEASURE_JS)
            entry['open'] = {'measure': opened, 'checks': checks(opened, compact)}
            page.locator('.preview-editor-layout').screenshot(path=str(OUT / f'{name}-editor.png'))
            drawer_actions = page.evaluate("""() => {
              const align = document.querySelector('.right-drawer button[title="Align to the selected exact feature"]');
              const review = document.querySelector('[data-similarity-alignment-drawer-review]');
              const r = (e) => e && e.getBoundingClientRect();
              const a = r(align), b = r(review);
              return a && b ? { align: [a.x, a.y, a.width, a.height], review: [b.x, b.y, b.width, b.height],
                sameRow: Math.abs(a.y - b.y) < 2, reviewRightOfAlign: b.x >= a.x + a.width - 1,
                reviewClipped: review.scrollWidth > review.clientWidth + 1 } : null; }""")
            entry['drawerActions'] = drawer_actions
            page.evaluate("() => { window.__GBDRAW_APP__.showRightDrawer = false; }")
            page.wait_for_timeout(500)
            svg_id = page.evaluate("async () => (await import('./js/state.js')).state.extractedFeatures.value.find(f => f.gene === 'livA')?.svg_id")
            feature = page.locator(f'.gbdraw-preview-surface svg [data-gbdraw-feature-id="{svg_id}"]').first
            feature.dispatch_event('click')
            popup = page.locator('.feature-popup[role="dialog"]').first
            popup.wait_for(state='visible', timeout=30000)
            entry['popupActions'] = page.evaluate("""() => {
              const popup = document.querySelector('.feature-popup[role="dialog"]');
              const align = [...popup.querySelectorAll('button')].find(b => b.title === 'Align to this exact feature');
              const review = popup.querySelector('[data-similarity-alignment-popup-review]');
              const help = popup.querySelector('.help-tip:has(#clicked-alignment-action-description) > button');
              const r = (e) => e && e.getBoundingClientRect();
              const a = r(align), b = r(review), h = r(help), pr = r(popup);
              return { sameRow: Math.abs(a.y + a.height / 2 - (b.y + b.height / 2)) < 2,
                reviewRightOfAlign: b.x >= a.x + a.width - 1, reviewClipped: review.scrollWidth > review.clientWidth + 1,
                helpVisible: Boolean(h && h.width > 0), insidePopup: b.x + b.width <= pr.x + pr.width + 0.5,
                alwaysOnParagraph: Boolean([...popup.querySelectorAll('p')].find(p => p.textContent.includes('Align applies resolved targets'))),
                describedBy: align.getAttribute('aria-describedby') }; }""")
            popup.screenshot(path=str(OUT / f'{name}-popup.png'))
            result['viewports'][name] = entry
            print(name, json.dumps({k: entry[k]['checks'] for k in ('closed', 'open')}), json.dumps(drawer_actions),
                  json.dumps(entry['popupActions']), flush=True)
            context.close()
        browser.close()
    (OUT / 'probe.json').write_text(json.dumps(result, indent=1))


if __name__ == '__main__':
    main(sys.argv[1])
