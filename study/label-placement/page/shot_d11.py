"""Look at the LP-D11 sections of the built page: usage shot_d11.py <out_dir> <shots_dir>."""
import sys
from pathlib import Path
from playwright.sync_api import sync_playwright

out = Path(sys.argv[1]).resolve(); shots = Path(sys.argv[2]).resolve(); shots.mkdir(parents=True, exist_ok=True)
body = (out / 'index.html').read_text(encoding='utf-8')
(out / 'preview.html').write_text('<!doctype html><html lang="ja"><head><meta charset="utf-8">'
                                  '<meta name="viewport" content="width=device-width, initial-scale=1"></head><body>'
                                  + body + '</body></html>', encoding='utf-8')
with sync_playwright() as p:
    browser = p.chromium.launch()
    for width in (1440, 390):
        page = browser.new_page(viewport={'width': width, 'height': 1000})
        errors = []
        page.on('console', lambda m: errors.append(m.text) if m.type == 'error' else None)
        page.on('pageerror', lambda e: errors.append(str(e)))
        page.goto((out / 'preview.html').as_uri(), wait_until='load')
        page.wait_for_timeout(800)
        print(width, 'scrollWidth', page.evaluate('document.documentElement.scrollWidth'), 'errors', errors[:3])
        ids = page.evaluate("[...document.querySelectorAll('[id^=d11]')].map(e => e.id)")
        print(width, 'ids', ids)
        for sel in ids:
            el = page.query_selector('#' + sel)
            el.scroll_into_view_if_needed(); page.wait_for_timeout(400)
            box = el.bounding_box()
            if box and box['height'] > 4:
                page.screenshot(path=str(shots / f'{sel}-{width}.png'), full_page=True,
                                clip={'x': 0, 'y': box['y'] + page.evaluate('window.scrollY'), 'width': width, 'height': min(box['height'], 2400)})
        broken = page.evaluate("[...document.querySelectorAll('[id^=d11] img')].filter(i => i.complete && !i.naturalWidth).length")
        print(width, 'broken images in d11', broken)
        page.close()
    browser.close()
