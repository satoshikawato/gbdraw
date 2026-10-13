"""Screenshot given element ids of the built page: shot_ids.py <out_dir> <shots_dir> <width> <id> [<id> ...]"""
import sys
from pathlib import Path
from playwright.sync_api import sync_playwright

out = Path(sys.argv[1]).resolve(); shots = Path(sys.argv[2]).resolve(); width = int(sys.argv[3])
with sync_playwright() as p:
    browser = p.chromium.launch()
    page = browser.new_page(viewport={'width': width, 'height': 1000})
    errors = []
    page.on('console', lambda m: errors.append(m.text) if m.type == 'error' else None)
    page.on('pageerror', lambda e: errors.append(str(e)))
    page.goto((out / 'preview.html').as_uri(), wait_until='load')
    for el_id in sys.argv[4:]:
        el = page.query_selector('#' + el_id)
        el.scroll_into_view_if_needed(); page.wait_for_timeout(600)
        box = el.bounding_box(); y = box['y'] + page.evaluate('window.scrollY')
        page.screenshot(path=str(shots / f'{el_id}-{width}.png'), full_page=True, clip={'x': 0, 'y': y, 'width': width, 'height': min(box['height'], 1500)})
    print('scrollWidth', page.evaluate('document.documentElement.scrollWidth'), 'errors:', errors[:5])
    browser.close()
