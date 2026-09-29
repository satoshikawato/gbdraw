#!/usr/bin/env python
"""S07 integration journey: the fixed S00 inputs in one page, in sequence.

Usage: python integration_journey.py <base-url>
Same page unless noted:
  A1 Load the fixed Vnig v41 Session, A2 Linear, A3 add the two Vibrio GenBank files and Generate (G03).
  A4 After that large Result: comparison commands and LOSATN/LOSATP, one Plot Title keystroke,
     Undo and Redo (G01, History). A5 Generate (G02) and SVG export of the current Result.
  B1 Load the BGC Gallery Session, B2 livA Align (G06), B3 livE Align with Keep (G07), B4 Undo/Redo,
  B5 Editor open/close and feature search (G04/G05), B6 Save -> Load in a fresh page.
Reuses the S00 Vnig inputs, the S04 state snapshot, the S05 BGC Candidate and the S06 layout checks.
Writes $S00_BASELINE_DIR/evidence/s07-integration/candidate-integration.json and screenshots.
Latencies are single journey samples, not budget evidence (the S00 harness owns the budgets).
"""
import hashlib
import importlib.util
import json
import os
import pathlib
import sys

from playwright.sync_api import sync_playwright

RESULTS = pathlib.Path(__file__).resolve().parents[1]
WT = RESULTS.parents[3]


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, RESULTS / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


vnig = load('vnig_linear_evidence', 's00/vnig_linear_evidence.py')
s04 = load('vnig_mode_session_acceptance', 's04/vnig_mode_session_acceptance.py')
s05 = load('bgc_align_acceptance', 's05/bgc_align_acceptance.py')
s06 = load('preview_layout_probe', 's06/preview_layout_probe.py')
OUT = pathlib.Path(os.environ.get('S00_BASELINE_DIR', '/home/kawato/gbdraw-baselines/gui-remediation-20260928')) / 'evidence/s07-integration'
s05.flows.OUT = OUT

# Resolves after the condition holds and two frames have painted; returns ms since the given start.
VISIBLE_JS = """async ([condition, start]) => {
  const test = new Function('app', 'return (' + condition + ');');
  while (!test(window.__GBDRAW_APP__)) await new Promise(r => setTimeout(r, 5));
  await new Promise(r => requestAnimationFrame(() => requestAnimationFrame(r)));
  return +(performance.now() - start).toFixed(1);
}"""


class Journey(s05.Candidate):
    def __init__(self, browser):
        super().__init__(browser, 'candidate', 0)
        self.tag = 'candidate-integration'
        self.log.update({'dialogs': [], 'latenciesMs': {}, 'g03': {}, 'checks': {}})
        self.page.on('dialog', lambda d: (self.log['dialogs'].append(d.message[:300]), d.accept()))

    def timed(self, key, action, condition):
        start = self.page.evaluate('() => performance.now()')
        action()
        self.log['latenciesMs'][key] = self.page.evaluate(VISIBLE_JS, [condition, start])

    def load_file(self, path, label):
        with self.page.expect_file_chooser() as chooser:
            self.page.get_by_role('button', name='Load Session', exact=True).first.click()
        chooser.value.set_files(str(path))
        self.page.wait_for_function('() => !window.__GBDRAW_APP__.sessionImportPending && window.__GBDRAW_APP__.results?.length > 0',
                                    timeout=1_200_000, polling=500)
        self.step(label, {'settleSeconds': s05.flows.settle(self.page)})

    def vnig_linear(self):
        self.page.goto(self.url, wait_until='domcontentloaded')
        self.page.wait_for_function('() => Object.keys(window.__GBDRAW_APP__?.paletteDefinitions || {}).length > 0',
                                    timeout=180_000)
        self.load_file(vnig.SESSIONS['v41'], 'A1-load-vnig-v41')
        self.log['g03']['A1'] = self.page.evaluate(s04.EXTRA)
        self.page.locator('button.app-mode-button[aria-label="Linear"]').click()
        self.page.wait_for_function("() => window.__GBDRAW_APP__.mode === 'linear'")
        self.log['g03']['A2'] = self.page.evaluate(s04.EXTRA)
        for index, genome in enumerate(vnig.GENOMES):
            if index > 0:
                self.page.get_by_role('button', name='Add sequence', exact=True).last.click()
                self.page.wait_for_function(f'() => window.__GBDRAW_APP__.linearSourceGroups.length > {index}')
            with self.page.expect_file_chooser() as chooser:
                self.page.locator('[data-linear-source-card]').nth(index).locator('.upload-zone').first.click()
            chooser.value.set_files(genome)
            self.page.wait_for_function(f"""() => window.__GBDRAW_APP__.linearSourceGroups[{index}]?.records?.length === 2""",
                                        timeout=300_000, polling=500)
        self.generate('A3-generate-vibrio')
        self.log['g03']['A3'] = self.page.evaluate(s04.EXTRA)

    def comparison_and_input(self):
        page, app = self.page, "window.__GBDRAW_APP__"
        pressed = "document.querySelector('[aria-label=\"{0}\"]')?.getAttribute('aria-pressed') === 'true'"
        mode = page.get_by_role('group', name='LOSAT Mode')
        self.timed('noComparisonToLosat', lambda: page.get_by_role('button', name='Run LOSAT for all adjacent pairs', exact=True).click(),
                   pressed.format('Run LOSAT for all adjacent pairs'))
        self.timed('losatnToLosatp', lambda: mode.get_by_role('button', name='LOSATP', exact=True).click(), "app.losatProgram === 'blastp'")
        self.timed('losatpToLosatn', lambda: mode.get_by_role('button', name='LOSATN', exact=True).click(), "app.losatProgram === 'blastn'")
        self.timed('losatToNoComparison', lambda: page.get_by_role('button', name='Set no comparison', exact=True).click(),
                   pressed.format('Set no comparison'))
        title = page.get_by_role('textbox', name='Plot Title')
        if not title.is_visible():
            page.get_by_label('Titles and Record Labels', exact=True).click()
        title.click()
        before = page.evaluate(f'() => {app}.form.plot_title')
        self.timed('plotTitleKeystroke', lambda: page.keyboard.press('X'), f"app.form.plot_title === {json.dumps(before + 'X')}")
        page.keyboard.press('Tab')
        s05.flows.settle(page)
        self.step('A4-comparison-and-input')
        page.get_by_role('button', name='Undo', exact=True).click()
        self.log['checks']['undoPlotTitle'] = page.evaluate(f'() => {app}.form.plot_title')
        page.get_by_role('button', name='Redo', exact=True).click()
        self.log['checks']['redoPlotTitle'] = page.evaluate(f'() => {app}.form.plot_title')
        page.get_by_role('button', name='Undo', exact=True).click()
        self.log['checks']['finalPlotTitle'] = page.evaluate(f'() => {app}.form.plot_title')
        self.log['checks']['comparisonAfterUndo'] = page.evaluate(f'() => {app}.linearComparisonGlobalAction')

    def generate_and_export(self):
        page = self.page
        key = page.evaluate("async () => (await import('./js/state.js')).state.resultGenerationKey.value")
        self.timed('generateAccepted', lambda: page.get_by_role('button', name='Generate Diagram', exact=True).click(),
                   'app.processing === true')
        page.wait_for_function("""async (key) => { const { state } = await import('./js/state.js');
            return !state.processing.value && (state.resultGenerationKey.value > key || state.errorLog.value !== null); }""",
                               arg=key, timeout=900_000, polling=1000)
        self.step('A5-generate', {'settleSeconds': s05.flows.settle(page),
                                  'generationError': page.evaluate('() => window.__GBDRAW_APP__.errorLog?.summary || null')})
        with page.expect_download(timeout=300_000) as download:
            page.get_by_role('button', name='SVG', exact=True).click()
        path = OUT / f'{self.tag}-A5-export.svg'
        download.value.save_as(str(path))
        exported = path.read_text()
        current = page.evaluate("async () => (await import('./js/state.js')).state.results.value[0].content")
        self.log['checks']['export'] = {
            'filename': download.value.suggested_filename, 'bytes': len(exported.encode()),
            'sameAsCurrentResult': hashlib.sha256(exported.encode()).hexdigest() == hashlib.sha256(current.encode()).hexdigest(),
            'currentBytes': len(current.encode()),
            'firstDifference': next((i for i, (a, b) in enumerate(zip(exported, current)) if a != b), None),
            'recordsPresent': {name: name in exported for name in ('NC_022349', 'NC_022359', 'NC_004603', 'NC_004605')}}

    def bgc(self):
        self.load_file(WT / s05.flows.FIXTURE, 'B1-load-bgc')
        self.align('livA', 'B2-livA')
        self.align('livE', 'B3-livE')
        self.page.get_by_role('button', name='Undo', exact=True).click()
        self.step('B4-undo', {'settleSeconds': s05.flows.settle(self.page)})
        self.page.get_by_role('button', name='Redo', exact=True).click()
        self.step('B4-redo', {'settleSeconds': s05.flows.settle(self.page)})

    def editor_and_search(self):
        page = self.page
        page.locator('.preview-editor-layout').scroll_into_view_if_needed()
        self.timed('editorOpen', lambda: page.get_by_role('button', name='Open editor', exact=True).click(), 'app.showRightDrawer === true')
        opened = page.evaluate(s06.MEASURE_JS)
        search = page.get_by_role('searchbox', name='Search features')
        search.fill('livE')
        self.timed('featureSearch', lambda: page.get_by_role('button', name='Search features', exact=True).click(),
                   '(app.previewFeatureSearchMatches || []).length > 0')
        matches = page.evaluate('() => (window.__GBDRAW_APP__.previewFeatureSearchMatches || []).length')
        self.shot('B5-editor-search')
        self.timed('editorClose', lambda: page.locator('.drawer-toggle').click(), 'app.showRightDrawer === false')
        closed = page.evaluate(s06.MEASURE_JS)
        self.log['checks']['editor'] = {'open': s06.checks(opened, False), 'closed': s06.checks(closed, False), 'searchMatches': matches}
        self.step('B5-editor-search')


def main(base_url):
    OUT.mkdir(parents=True, exist_ok=True)
    s05.flows.APPS['candidate'] = (f'{base_url}/gbdraw/web/index.html', WT / s05.flows.FIXTURE)
    with sync_playwright() as p:
        browser = p.chromium.launch()
        run = Journey(browser)
        try:
            run.vnig_linear()
            run.comparison_and_input()
            run.generate_and_export()
            run.bgc()
            run.editor_and_search()
            run.save_and_reload('B6-fresh-load')
        except Exception as error:
            run.log['fatal'] = repr(error)[:2000]
            print(run.tag, 'FATAL', repr(error)[:500], flush=True)
        finally:
            run.save()
            print(json.dumps({k: run.log[k] for k in ('latenciesMs', 'g03', 'checks', 'dialogs', 'pageErrors')}, indent=1)[:6000])
        browser.close()


if __name__ == '__main__':
    main(sys.argv[1])
