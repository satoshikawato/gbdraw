"""Diagnose late long animation frames after the post-Result comparison switch (not budget evidence)."""
import json
import pathlib
import sys

from playwright.sync_api import sync_playwright

WT = str(pathlib.Path(__file__).resolve().parents[5])
URL = sys.argv[1] if len(sys.argv) > 1 else 'http://127.0.0.1:4305'
REPS = int(sys.argv[2]) if len(sys.argv) > 2 else 60

INIT = r"""(() => {
  const d = window.__DIAG__ = { loaf: [], mutations: [] };
  new PerformanceObserver((list) => { for (const e of list.getEntries()) d.loaf.push({
    start: e.startTime, duration: e.duration, blocking: e.blockingDuration, renderStart: e.renderStart,
    styleAndLayoutStart: e.styleAndLayoutStart, firstUIEventTimestamp: e.firstUIEventTimestamp,
    scripts: (e.scripts || []).map((s) => ({ invoker: s.invoker, invokerType: s.invokerType, fn: s.sourceFunctionName,
      url: (s.sourceURL || '').split('/').slice(-2).join('/'), duration: s.duration, forcedStyleAndLayout: s.forcedStyleAndLayoutDuration })) }); })
    .observe({ type: 'long-animation-frame', buffered: true });
  const describe = (n) => n && n.nodeType === 1 ? (n.tagName.toLowerCase() + (n.id ? '#' + n.id : '') + '.' + String(n.className?.baseVal ?? n.className ?? '').split(' ').slice(0, 3).join('.')) : String(n?.nodeName);
  const start = () => new MutationObserver((records) => { const t = performance.now();
    for (const r of records) d.mutations.push({ t, type: r.type, target: describe(r.target), attr: r.attributeName,
      added: r.addedNodes.length, removed: r.removedNodes.length }); })
    .observe(document.documentElement, { subtree: true, childList: true, attributes: true, characterData: true });
  if (document.documentElement) start(); else addEventListener('DOMContentLoaded', start);
})();"""


def main():
    with sync_playwright() as p:
        browser = p.chromium.launch()
        page = browser.new_context(viewport={'width': 1440, 'height': 900}).new_page()
        page.add_init_script(INIT)
        page.on('dialog', lambda d: d.accept())
        page.goto(f'{URL}/gbdraw/web/index.html', wait_until='domcontentloaded')
        page.wait_for_function('() => Object.keys(window.__GBDRAW_APP__?.paletteDefinitions || {}).length > 0', timeout=180000)
        page.get_by_role('button', name='Linear', exact=True).click()
        for index, name in ((1, 'MG1655.gbk'), (2, 'Sakai.gbk')):
            if page.get_by_test_id(f'linear-genbank-{index}').count() == 0:
                page.get_by_role('button', name='Add sequence', exact=True).last.click()
            page.get_by_test_id(f'linear-genbank-{index}').set_input_files(f'{WT}/tests/test_inputs/{name}')
        page.wait_for_timeout(3000)
        for _ in range(int(sys.argv[3]) if len(sys.argv) > 3 else 1):
            page.get_by_role('button', name='Generate Diagram', exact=True).click()
            page.wait_for_function('() => window.__GBDRAW_APP__.processing === true', timeout=120000)
            page.wait_for_function('() => window.__GBDRAW_APP__.processing === false', timeout=1200000, polling=500)
            page.wait_for_timeout(1500)
        page.evaluate("""async () => {
          const app = window.__GBDRAW_APP__;
          const feature = app.extractedFeatures.find((c) => c?.svg_id && String(c.type || c.feature_type || '').toUpperCase() === 'CDS');
          app.openFeatureEditorFromList(feature, null);
          app.clickedFeature.labelVisibility = 'on';
          const update = app.updateClickedFeatureLabelText();
          for (let i = 0; i < 100 && !app.globalLabelModeDialog?.show; i += 1) await new Promise((r) => setTimeout(r, 50));
          if (app.globalLabelModeDialog?.show) app.handleGlobalLabelModeChoice('whitelist_only');
          await update; app.closeRightDrawer(); }""")
        page.wait_for_timeout(5000)
        late = []
        for rep in range(REPS):
            for target, name in (('losat', 'Run LOSAT for all adjacent pairs'), ('none', 'Set no comparison')):
                page.evaluate('() => { window.__DIAG__.loaf.length = 0; window.__DIAG__.mutations.length = 0; }')
                t0 = page.evaluate('() => performance.now()')
                page.get_by_role('button', name=name, exact=True).click()
                page.wait_for_timeout(1600)
                frames = page.evaluate('() => window.__DIAG__.loaf')
                for f in frames:
                    if f['start'] - t0 > 100 and f['duration'] >= 50:
                        muts = page.evaluate("""([a, b]) => { const m = window.__DIAG__.mutations.filter((x) => x.t >= a && x.t <= b);
                          const c = {}; for (const x of m) { const k = x.type + ' ' + x.target + (x.attr ? '@' + x.attr : ''); c[k] = (c[k] || 0) + 1; }
                          return Object.entries(c).sort((p, q) => q[1] - p[1]).slice(0, 8); }""", [f['start'] - 150, f['start'] + f['duration']])
                        late.append({'rep': rep, 'to': target, 'afterClickMs': round(f['start'] - t0), 'frame': f, 'mutations': muts})
                        print(json.dumps(late[-1])[:1600], flush=True)
        print('late frames', len(late), 'of', REPS * 2, 'switches')
        browser.close()


if __name__ == '__main__':
    main()
