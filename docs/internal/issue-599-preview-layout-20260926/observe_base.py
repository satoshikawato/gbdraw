"""Archived Issue #599 baseline observer; generates diagnostic JSON.

Use a disposable clone of baseline d457b7189b137185a8dec800819a312c30b969fa.
Copy this script into docs/internal/issue-599-preview-layout-20260926/ there.
Prepare a source-matching browser wheel and run from that repository root:
  python docs/internal/issue-599-preview-layout-20260926/observe_base.py
The output overwrites base-observations.json beside this script in that clone.
Keep the committed baseline evidence unchanged. This observer asserts the bug
exists; it is not an after-fix regression-test expectation.
"""

import functools
import hashlib
import importlib.metadata
import zipfile
import http.server
import json
from pathlib import Path
import subprocess
import threading

from playwright.sync_api import sync_playwright


ROOT = Path(__file__).resolve().parents[3]
OUTPUT = Path(__file__).with_name("base-observations.json")


class QuietHandler(http.server.SimpleHTTPRequestHandler):
    def log_message(self, *_args):
        pass


READ_COMPOSITION = """async () => {
  const app = window.__GBDRAW_APP__;
  const c = await import('./js/app/legend-layout/composition-actions.js');
  const svg = app.svgContainer.querySelector('svg');
  const binding = c.bindCompositionMetadata(svg);
  return {
    deltas: c.compositionUserDeltas(svg),
    cursors: ['legend', 'title'].map(role => ({
      role, cursor: binding[role].targets[0]
        ? getComputedStyle(binding[role].targets[0]).cursor : null
    })),
    primaryIds: binding.primary.targets.map(target => target.id),
    layoutEdit: app.layoutRepositionMode
  };
}"""


def main():
    server = http.server.ThreadingHTTPServer(
        ("127.0.0.1", 0), functools.partial(QuietHandler, directory=str(ROOT))
    )
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    wheel = next((ROOT / "gbdraw/web").glob("gbdraw-*.whl"))
    with zipfile.ZipFile(wheel) as archive:
        mismatches = [name for name in archive.namelist()
            if name.startswith("gbdraw/") and name.endswith((".py", ".toml"))
            and (ROOT / name).exists()
            and archive.read(name) != (ROOT / name).read_bytes()]
    if mismatches:
        server.shutdown()
        server.server_close()
        raise RuntimeError(f"Browser wheel differs from source: {mismatches}")
    observations = {
        "playwrightVersion": importlib.metadata.version("playwright"),
        "wheelSha256": hashlib.sha256(wheel.read_bytes()).hexdigest(),
        "wheelSourceMismatches": mismatches,
        "baseSha": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "scope": "Unchanged runtime; real Generate, actual pointer drag, DOM geometry",
        "generations": [],
        "chrome": [],
    }
    try:
        with sync_playwright() as pw:
            browser = pw.chromium.launch()
            observations["chromiumVersion"] = browser.version
            page = browser.new_page(viewport={"width": 1440, "height": 1000})
            page.set_default_timeout(180_000)
            page.on("dialog", lambda dialog: dialog.accept())
            page.goto(f"http://127.0.0.1:{server.server_port}/gbdraw/web/index.html")
            page.wait_for_function(
                "() => window.__GBDRAW_APP__ && Object.keys(window.__GBDRAW_APP__.paletteDefinitions || {}).length"
            )
            for mode, fixture in [
                ("circular", "HmmtDNA_basic_circular"),
                ("linear", "lambda_basic_linear"),
            ]:
                page.locator(
                    'input[type="file"][accept*="application/json"][accept*="application/gzip"]'
                ).set_input_files(
                    ROOT / f"gbdraw/web/gallery/sessions/{fixture}.gbdraw-session.json"
                )
                page.wait_for_function(
                    "mode => window.__GBDRAW_APP__.mode === mode && window.__GBDRAW_APP__.results.length > 0",
                    arg=mode,
                )
                first = page.evaluate("""async () => {
                  const app = window.__GBDRAW_APP__;
                  app.form.plot_title = 'Issue 599 composition observation';
                  app.adv.plot_title_position = 'top';
                  return await app.runAnalysis();
                }""")
                if first.get("status") != "ok":
                    raise RuntimeError(f"Initial {mode} Generate failed: {first}")
                initial = page.evaluate(READ_COMPOSITION)
                if not initial["layoutEdit"]:
                    page.get_by_role("button", name="Toggle layout edit mode", exact=True).click()
                # Move the camera below the search overlay; no SVG/Result edit.
                page.evaluate("""async () => {
                  const app = window.__GBDRAW_APP__;
                  app.zoom = 0.6;
                  app.canvasPan.y = 150;
                  await window.Vue.nextTick();
                }""")
                page.wait_for_timeout(300)
                selectors = ['#legend text', '#plot_title text']
                if mode == "linear":
                    selectors.append('#length_bar text')
                pointer_hits = []
                for selector in selectors:
                    target = page.locator(f'[aria-label="Result Preview"] {selector}').first
                    box = target.bounding_box()
                    if not box:
                        raise RuntimeError(f"Missing visible drag target: {mode}/{selector}")
                    page.evaluate("""({x,y}) => {
                      const app = window.__GBDRAW_APP__;
                      const viewport = document.querySelector('.preview-viewport').getBoundingClientRect();
                      app.canvasPan.x += viewport.x + viewport.width / 2 - x;
                      app.canvasPan.y += viewport.y + viewport.height / 2 - y;
                    }""", {"x":box["x"]+box["width"]/2,"y":box["y"]+box["height"]/2})
                    page.wait_for_timeout(350)
                    box = target.bounding_box()
                    x, y = box["x"] + box["width"] / 2, box["y"] + box["height"] / 2
                    hit = page.evaluate("""({x,y,selector}) => {
                      const target = document.elementFromPoint(x,y);
                      return {tag:target?.tagName,id:target?.id,
                        intended:Boolean(target?.closest(selector.split(' ')[0]))};
                    }""", {"x": x, "y": y, "selector": selector})
                    pointer_hits.append({"selector": selector, **hit})
                    if not hit["intended"]:
                        raise RuntimeError(f"Drag intercepted: {mode}/{selector}: {hit}")
                    page.mouse.move(x, y)
                    page.mouse.down()
                    page.mouse.move(x + 24, y + 12, steps=6)
                    page.mouse.up()
                    page.wait_for_timeout(300)
                dragged = page.evaluate(READ_COMPOSITION)
                if not all(any(v != 0 for v in dragged["deltas"][role]) for role in ["legend", "title"]):
                    raise RuntimeError(f"Drag did not move {mode} legend/title: {dragged}")
                outcome = page.evaluate("async () => await window.__GBDRAW_APP__.runAnalysis()")
                regenerated = page.evaluate(READ_COMPOSITION)
                observations["generations"].append({
                    "mode": mode, "fixture": fixture, "first": first,
                    "initial": initial, "dragged": dragged, "pointerHits": pointer_hits,
                    "generate": outcome, "regenerated": regenerated,
                })
                print(json.dumps(observations["generations"][-1]), flush=True)
            for width, height in [(1440, 1000), (1024, 768), (900, 768), (768, 768), (390, 844)]:
                page.set_viewport_size({"width": width, "height": height})
                for open_drawer in [False, True]:
                    page.evaluate("""async open => {
                      const app = window.__GBDRAW_APP__;
                      if (open) app.openRightDrawerTab('features');
                      else app.closeRightDrawer();
                      await window.Vue.nextTick();
                    }""", open_drawer)
                    page.wait_for_timeout(350)
                    bounds = page.evaluate("""() => {
                      const rect = selector => {
                        const r = document.querySelector(selector).getBoundingClientRect();
                        return {x:r.x,y:r.y,width:r.width,height:r.height,right:r.right,bottom:r.bottom};
                      };
                      const search = rect('.preview-feature-search');
                      const toolbar = rect('.preview-controls');
                      const preview = rect('.preview-viewport');
                      const intersect = (a,b) => Math.max(0,Math.min(a.right,b.right)-Math.max(a.x,b.x))
                        * Math.max(0,Math.min(a.bottom,b.bottom)-Math.max(a.y,b.y));
                      return {search,toolbar,preview,searchToolbarOverlap:intersect(search,toolbar),
                        searchOutsidePreview:search.x < preview.x-1 || search.right > preview.right+1,
                        searchInlineTransform:document.querySelector('.preview-feature-search').style.transform,
                        searchComputedTransform:getComputedStyle(document.querySelector('.preview-feature-search')).transform};
                    }""")
                    observations["chrome"].append({"viewport": [width, height], "drawerOpen": open_drawer, **bounds})
                    print(json.dumps(observations["chrome"][-1]), flush=True)
            browser.close()
    finally:
        server.shutdown()
        server.server_close()
        OUTPUT.write_text(json.dumps(observations, ensure_ascii=False, indent=2) + "\n")


if __name__ == "__main__":
    main()
