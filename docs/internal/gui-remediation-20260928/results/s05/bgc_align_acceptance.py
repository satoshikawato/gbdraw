#!/usr/bin/env python
"""S05 BGC acceptance on the candidate build, reusing the S00 flow harness.

Usage: python bgc_align_acceptance.py <base-url> [flow ...]
Flows: load-align (Gallery Load, livA Align, livE Align, Save -> fresh Load),
       review-right (livA Review -> All right, livE Align),
       generate-align (Gallery Load, Generate, livA Align), para (Gallery Load, parA Align),
       fresh-losatp (new Linear session from the five BGC files, Generate runs LOSATP, livA Align).
Writes $S00_BASELINE_DIR/evidence/s05-bgc/candidate-<flow>.json and screenshots.
"""
import importlib.util
import json
import os
import pathlib
import re
import sys

from playwright.sync_api import sync_playwright

S00 = pathlib.Path(__file__).resolve().parents[1] / 's00' / 'bgc_align_flows.py'
spec = importlib.util.spec_from_file_location('bgc_align_flows', S00)
flows = importlib.util.module_from_spec(spec)
spec.loader.exec_module(flows)

WT = pathlib.Path(__file__).resolve().parents[5]
flows.OUT = pathlib.Path(os.environ.get('S00_BASELINE_DIR', '/home/kawato/gbdraw-baselines/gui-remediation-20260928')) / 'evidence/s05-bgc'
flows.GENES = ['livA', 'livE', 'parA', 'racM', 'racL']
BGC_FILES = ['BGC0000708.gbk', 'BGC0000709.gbk', 'BGC0000711.gbk', 'BGC0000712.gbk', 'BGC0000713.gbk']
INTERNAL = re.compile(r'\bh_[a-z0-9]{20,}\b|\bf[0-9a-f]{8}\b|\brecord-\d+\b|Internal ID|Feature ID|Record ID')

INSPECTOR_JS = r"""async () => {
  const app = window.__GBDRAW_APP__;
  const previous = { open: app.showRightDrawer, tab: app.rightDrawerTab };
  app.showRightDrawer = true; app.rightDrawerTab = 'orthogroups';
  await window.Vue.nextTick();
  const node = document.querySelector('[data-similarity-alignment-plan-inspector]');
  if (node && !node.open) { node.open = true; await window.Vue.nextTick(); }
  const text = node ? node.innerText.trim() : null;
  const model = JSON.parse(JSON.stringify(app.similarityAlignmentPlanInspector || null));
  app.showRightDrawer = previous.open; app.rightDrawerTab = previous.tab;
  await window.Vue.nextTick();
  return { text, model };
}"""
LOSAT_JS = "() => window.__BGC_WORKER__.constructions.filter(url => url.includes('losat')).length"


CAPTURE_JS = r"""(() => {
  window.__S05_HELPER__ = [];
  const post = Worker.prototype.postMessage;
  Worker.prototype.postMessage = function (message, ...rest) {
    if (message?.type === 'helper' && message.operation === 'resolveSimilarityAlignment') {
      try { window.__S05_HELPER__.push(JSON.parse(JSON.stringify({ request: message.payload?.request,
        projection: message.payload?.projection ?? null }))); } catch (_error) {}
    }
    return post.call(this, message, ...rest);
  };
  window.__S05_HELPER_ERRORS__ = [];
  const Native = window.Worker;
  window.Worker = new Proxy(Native, { construct(target, args) {
    const worker = Reflect.construct(target, args, target);
    worker.addEventListener('message', (event) => {
      const m = event.data || {};
      if (m.type === 'helper' && m.ok === false) {
        try { window.__S05_HELPER_ERRORS__.push(JSON.parse(JSON.stringify(m.error ?? null))); } catch (_error) {}
      }
    });
    return worker;
  } });
})();"""


class Candidate(flows.Run):
    def __init__(self, browser, app, flow):
        super().__init__(browser, app, flow)
        self.ctx.add_init_script(CAPTURE_JS)
        self.worker_console = []

    def step(self, name, extra=None):
        extra = dict(extra or {})
        extra['losatWorkerConstructions'] = self.page.evaluate(LOSAT_JS)
        extra['errorDetail'] = self.page.evaluate(
            "() => JSON.parse(JSON.stringify(window.__GBDRAW_APP__.similarityAlignmentError || window.__GBDRAW_APP__.errorLog || null))")
        extra['inspector'] = self.page.evaluate(INSPECTOR_JS)
        text = extra['inspector']['text'] or ''
        extra['inspectorInternalIds'] = sorted(set(INTERNAL.findall(text)))
        super().step(name, extra)

    def result_shot(self, name):
        has_svg = self.page.evaluate("() => Boolean(document.querySelector('.gbdraw-preview-surface svg'))")
        return super().result_shot(name) if has_svg else None

    def capture_review(self, label):
        review = super().capture_review(label)
        review['displayedInternalIds'] = sorted(set(INTERNAL.findall(review['fullText'])))
        return review

    def generate(self, label):
        key = self.page.evaluate("async () => (await import('./js/state.js')).state.resultGenerationKey.value")
        self.page.get_by_role('button', name='Generate Diagram', exact=True).click()
        self.page.wait_for_function(
            """async (key) => { const { state } = await import('./js/state.js');
               return !state.processing.value && (state.resultGenerationKey.value > key || state.errorLog.value !== null); }""",
            arg=key, timeout=900000, polling=1000)
        self.step(label, {'settleSeconds': flows.settle(self.page),
                          'generationError': self.page.evaluate('() => window.__GBDRAW_APP__.errorLog?.summary || null')})

    def fresh_losatp(self, label):
        """A new Linear session from the five BGC GenBank files; Generate runs LOSATP."""
        if os.environ.get('S05_WORKER_CONSOLE'):
            self.page.on('worker', lambda worker: worker.on('console', lambda m: self.worker_console.append(m.text)))
        self.page.goto(self.url, wait_until='domcontentloaded')
        self.page.wait_for_function('() => Object.keys(window.__GBDRAW_APP__?.paletteDefinitions || {}).length > 0',
                                    timeout=180000)
        self.page.locator('button.app-mode-button[aria-label="Linear"]').click()
        self.page.wait_for_function("() => window.__GBDRAW_APP__.mode === 'linear'")
        for index, name in enumerate(BGC_FILES):
            if index > 0:
                self.page.get_by_role('button', name='Add sequence', exact=True).last.click()
                self.page.wait_for_function(f'() => window.__GBDRAW_APP__.linearSourceGroups.length > {index}')
            with self.page.expect_file_chooser() as chooser:
                self.page.locator('[data-linear-source-card]').nth(index).locator('.upload-zone').first.click()
            chooser.value.set_files(str(WT / 'tests/test_inputs' / name))
            self.page.wait_for_function(
                f"""() => window.__GBDRAW_APP__.linearSourceGroups[{index}]?.sequence?.gb?.name === {name!r}""",
                timeout=300000)
        self.page.get_by_role('button', name='Run LOSAT for all adjacent pairs', exact=True).click()
        self.page.get_by_role('group', name='LOSAT Mode').get_by_role('button', name='LOSATP', exact=True).click()
        self.page.wait_for_function("() => window.__GBDRAW_APP__.losatProgram === 'blastp'")
        # The plain local server is not cross-origin isolated; use serial LOSAT
        # as the existing browser specs do (threaded mode fails here on dev too).
        self.page.evaluate("() => { window.__GBDRAW_APP__.losat.executionMode = 'serial'; }")
        self.log['deviations'].append('LOSAT execution mode set to serial (no cross-origin isolation)')
        flows.settle(self.page)
        self.log['beforeGenerate'] = {'losatWorkerConstructions': self.page.evaluate(LOSAT_JS),
                                      'worker': self.page.evaluate('() => window.__BGC_WORKER__')}
        self.generate(label)

    def save_and_reload(self, label):
        self.close_popup()
        with self.page.expect_download(timeout=300000) as download:
            self.page.evaluate("() => { window.__GBDRAW_APP__.sessionTitle = 's05-bgc'; return window.__GBDRAW_APP__.saveSessionWithTitle(); }")
        saved = flows.OUT / f'{self.tag}-saved.gbdraw-session.json'
        download.value.save_as(str(saved))
        self.page.close()
        self.page = self.ctx.new_page()
        self.page.on('pageerror', lambda e: self.log['pageErrors'].append(str(e)))
        self.page.goto(self.url, wait_until='domcontentloaded')
        self.page.wait_for_function('() => window.__GBDRAW_APP__')
        with self.page.expect_file_chooser() as chooser:
            self.page.get_by_role('button', name='Load Session').click()
        chooser.value.set_files(str(saved))
        self.page.wait_for_function('() => window.__GBDRAW_APP__.results?.length > 0', timeout=300000)
        self.step(label, {'settleSeconds': flows.settle(self.page)})


def main(base_url, names):
    flows.OUT.mkdir(parents=True, exist_ok=True)
    flows.APPS['candidate'] = (f'{base_url}/gbdraw/web/index.html', WT / flows.FIXTURE)
    with sync_playwright() as p:
        browser = p.chromium.launch()
        for name in names:
            run = Candidate(browser, 'candidate', 0)
            run.tag = f'candidate-{name}'
            run.log['flowName'] = name
            try:
                if name == 'fresh-losatp':
                    run.fresh_losatp('01-generate-losatp')
                    with run.page.expect_download(timeout=300000) as download:
                        run.page.evaluate("() => { window.__GBDRAW_APP__.sessionTitle = 's05-fresh'; return window.__GBDRAW_APP__.saveSessionWithTitle(); }")
                    download.value.save_as(str(flows.OUT / f'{run.tag}-generated.gbdraw-session.json'))
                    run.align('livA', '02-livA')
                    (flows.OUT / f'{run.tag}-helper-payloads.json').write_text(
                        json.dumps({'payloads': run.page.evaluate('() => window.__S05_HELPER__'),
                                    'errors': run.page.evaluate('() => window.__S05_HELPER_ERRORS__'),
                                    'workerConsole': [t for t in run.worker_console if 'S05-DEBUG' in t]}, indent=1))
                    continue
                run.load()
                if name == 'load-align':
                    run.align('livA', '01-livA')
                    run.align('livE', '02-livE')
                    run.save_and_reload('03-fresh-load')
                elif name == 'review-right':
                    run.align('livA', '01-livA-review', 'All selected arrows right')
                    run.align('livE', '02-livE')
                elif name == 'generate-align':
                    run.generate('01-generate')
                    run.align('livA', '02-livA')
                elif name == 'para':
                    run.align('parA', '01-parA')
            except Exception as error:
                run.log['fatal'] = repr(error)[:2000]
                print(run.tag, 'FATAL', repr(error)[:500], flush=True)
            finally:
                run.save()
        browser.close()


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2:] or ['load-align', 'review-right', 'generate-align', 'para'])
