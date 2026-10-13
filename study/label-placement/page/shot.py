"""One look at the built page before publishing: wrap the fragment in a document and screenshot parts of it.

usage: shot.py <out_dir> <shots_dir> [width]
"""
import sys
from pathlib import Path
from playwright.sync_api import sync_playwright

out = Path(sys.argv[1]).resolve(); shots = Path(sys.argv[2]).resolve(); shots.mkdir(parents=True, exist_ok=True)
width = int(sys.argv[3]) if len(sys.argv) > 3 else 1440
body = (out / 'index.html').read_text(encoding='utf-8')
(out / 'preview.html').write_text('<!doctype html><html lang="ja"><head><meta charset="utf-8">'
                                  '<meta name="viewport" content="width=device-width, initial-scale=1"></head><body>'
                                  + body + '</body></html>', encoding='utf-8')
with sync_playwright() as p:
    browser = p.chromium.launch()
    page = browser.new_page(viewport={'width': width, 'height': 1000})
    errors = []
    page.on('console', lambda m: errors.append(m.text) if m.type == 'error' else None)
    page.on('pageerror', lambda e: errors.append(str(e)))
    page.goto((out / 'preview.html').as_uri(), wait_until='load')
    page.wait_for_timeout(800)
    print('scrollWidth', page.evaluate('document.documentElement.scrollWidth'), 'viewport', width,
          'height', page.evaluate('document.documentElement.scrollHeight'))
    page.screenshot(path=str(shots / f'top-{width}.png'), full_page=True, clip={'x': 0, 'y': 0, 'width': width, 'height': 2600})
    for sel, name in [('section.case', 'case1'), ('section.variant', 'variant1'), ('figure.sheet', 'sheet1'), ('table.compare', 'compare')]:
        el = page.query_selector(sel)
        if el:
            el.scroll_into_view_if_needed()
            page.wait_for_timeout(500)
            box = el.bounding_box()
            box['y'] += page.evaluate('window.scrollY')
            page.screenshot(path=str(shots / f'{name}-{width}.png'), full_page=True,
                            clip={'x': 0, 'y': box['y'], 'width': width, 'height': min(box['height'], 1700)})
    print('errors:', errors[:5])
    browser.close()
