#!/usr/bin/env python
"""S04 G03 acceptance: the S00 Vnig journey plus mode round trip, Save -> fresh Load and a rejected Load.

Usage: python vnig_mode_session_acceptance.py <base-url> <v41|v44>
Reuses the S00 evidence harness (results/s00/vnig_linear_evidence.py) for the Worker tracking,
state snapshot, session inputs and genomes. Drives real UI controls except Save, which calls the
same saveSessionWithTitle() action as the Save dialog. Writes <OUT>/cell-candidate-<session>.json.
"""
import gzip
import importlib.util
import json
import os
import pathlib
import sys
import time

from playwright.sync_api import sync_playwright

S00 = pathlib.Path(__file__).resolve().parents[1] / 's00' / 'vnig_linear_evidence.py'
spec = importlib.util.spec_from_file_location('vnig_linear_evidence', S00)
evidence = importlib.util.module_from_spec(spec)
spec.loader.exec_module(evidence)

OUT = pathlib.Path(os.environ.get('S00_BASELINE_DIR', '/home/kawato/gbdraw-baselines/gui-remediation-20260928')) / 'evidence/s04-vnig'
EXTRA = r"""
async () => {
  const { state: s } = await import('./js/state.js');
  return {
    mode: s.mode.value, plot_title: s.form.plot_title, plot_title_font_size: s.adv.plot_title_font_size,
    def_font_size: s.adv.def_font_size, arrange_in_rows: s.linearRecordLayoutEnabled.value,
    linear_show_replicon: s.adv.linear_show_replicon, accession: s.adv.linear_accession_visibility,
    length: s.adv.linear_length_visibility, typography_linked: s.linearTypographyLinked.value,
    results: s.results.value.length, resultGenerationKey: s.resultGenerationKey.value,
    linearSeqs: s.linearSeqs.length
  };
}
"""


def main(base_url, session_key):
    OUT.mkdir(parents=True, exist_ok=True)
    cell = f'candidate-{session_key}'
    record = {'cell': cell, 'url': base_url, 'session': evidence.SESSIONS[session_key], 'stages': {}, 'extra': {},
              'dialogs': [], 'pageErrors': [], 'externalRequests': [], 'timings': {}}
    stage = {'name': 'S0'}
    url = f'{base_url}/gbdraw/web/index.html'

    def write():
        (OUT / f'cell-{cell}.json').write_text(json.dumps(record, indent=1, ensure_ascii=False))

    def snap(page, name):
        page.evaluate('() => window.Vue?.nextTick?.()')
        record['stages'][name] = page.evaluate(evidence.SNAPSHOT, {'stage': name, 'app': 'candidate'})
        record['extra'][name] = page.evaluate(EXTRA)
        write()

    def poll(page, js, timeout_s, interval=1.0):
        end = time.time() + timeout_s
        while time.time() < end:
            if page.evaluate(js):
                return True
            time.sleep(interval)
        return False

    def load_session(page, path, name):
        with page.expect_file_chooser() as chooser:
            page.get_by_role('button', name='Load Session', exact=True).first.click()
        chooser.value.set_files(path)
        end = time.time() + 1200
        while not any(d['stage'] == name for d in record['dialogs']):
            if time.time() > end:
                raise TimeoutError(f'{name}: Load did not report an outcome within 20 minutes')
            page.wait_for_timeout(500)  # a Playwright call dispatches the alert handler; time.sleep does not
        poll(page, '() => !window.__GBDRAW_APP__.sessionImportPending', 120, 0.5)
        time.sleep(2)

    def open_page(context):
        page = context.new_page()
        page.on('pageerror', lambda e: record['pageErrors'].append({'stage': stage['name'], 'message': str(e)}))

        def on_dialog(d):
            record['dialogs'].append({'stage': stage['name'], 'type': d.type, 'message': d.message})
            d.accept()
        page.on('dialog', on_dialog)
        page.goto(url, wait_until='domcontentloaded')
        page.wait_for_function('() => Object.keys(window.__GBDRAW_APP__?.paletteDefinitions || {}).length > 0',
                               timeout=180_000)
        return page

    def switch_mode(page, mode):
        page.locator(f'button.app-mode-button[aria-label="{mode.capitalize()}"]').click()
        poll(page, f"() => window.__GBDRAW_APP__.mode === '{mode}'", 30)
        time.sleep(1)

    with sync_playwright() as p:
        browser = p.chromium.launch()
        context = browser.new_context(viewport={'width': 1440, 'height': 900}, accept_downloads=True)

        def route(r):
            if '127.0.0.1' in r.request.url:
                return r.continue_()
            record['externalRequests'].append(r.request.url)
            return r.abort()
        context.route('**/*', route)
        context.add_init_script(evidence.WORKER_TRACKING)
        page = open_page(context)
        snap(page, 'S0')

        stage['name'] = 'S1'
        load_session(page, evidence.SESSIONS[session_key], 'S1')
        snap(page, 'S1')

        stage['name'] = 'S2'
        switch_mode(page, 'linear')
        snap(page, 'S2')

        stage['name'] = 'S3'
        for index, genome in enumerate(evidence.GENOMES):
            if index > 0:
                page.get_by_role('button', name='Add sequence', exact=True).last.click()
                poll(page, f"() => window.__GBDRAW_APP__.linearSourceGroups.length > {index}", 30, 0.2)
            with page.expect_file_chooser() as chooser:
                page.locator('[data-linear-source-card]').nth(index).locator('.upload-zone').first.click()
            chooser.value.set_files(genome)
            name = pathlib.Path(genome).name
            poll(page, f"""() => {{ const g = window.__GBDRAW_APP__.linearSourceGroups[{index}];
                return Boolean(g && g.sequence?.gb?.name === {json.dumps(name)} && g.records.length === 2); }}""", 300)
        time.sleep(2)
        snap(page, 'S3')

        stage['name'] = 'S4'
        key = page.evaluate("async () => (await import('./js/state.js')).state.resultGenerationKey.value")
        page.get_by_role('button', name='Generate Diagram', exact=True).click()
        record['timings']['S4_settled'] = poll(page, f"""async () => {{ const {{ state }} = await import('./js/state.js');
            return !state.processing.value && (state.resultGenerationKey.value > {key} || state.errorLog.value !== null); }}""",
                                               900, 5)
        time.sleep(2)
        snap(page, 'S4')
        page.screenshot(path=str(OUT / f'{cell}-S4.png'), full_page=True)

        # S5: back to Circular restores the Circular title and fonts.
        stage['name'] = 'S5'
        switch_mode(page, 'circular')
        snap(page, 'S5')

        # S6: Save with the Linear draft edited in this page, then Load into a fresh page.
        stage['name'] = 'S6'
        page.evaluate("() => { window.__GBDRAW_APP__.sessionTitle = 's04-vnig-roundtrip'; }")
        with page.expect_download(timeout=300_000) as download:
            record['saveOutcome'] = page.evaluate('() => window.__GBDRAW_APP__.saveSessionWithTitle()')
        saved = OUT / f'{cell}-roundtrip.gbdraw-session.json'
        download.value.save_as(str(saved))
        raw = saved.read_bytes()
        saved_json = json.loads(gzip.decompress(raw) if raw[:2] == b'\x1f\x8b' else raw)
        record['savedModeProfiles'] = {mode: {k: v for k, v in (profile.get('values') or {}).items()
                                              if k in ('plot_title', 'plot_title_font_size', 'def_font_size')}
                                       for mode, profile in (saved_json['config'].get('modeProfiles') or {}).get('profiles', {}).items()}
        record['savedLinearRecordLayoutEnabled'] = (saved_json['config'].get('linearRecordLayout') or {}).get('enabled')
        page.close()

        stage['name'] = 'S7'
        page = open_page(context)
        load_session(page, str(saved), 'S7')
        snap(page, 'S7')
        stage['name'] = 'S8'
        switch_mode(page, 'linear')
        snap(page, 'S8')

        # S9: a rejected Load (development-only version 43) keeps the current Session.
        stage['name'] = 'S9'
        rejected = dict(saved_json, version=43)
        rejected_path = OUT / f'{cell}-rejected-v43.gbdraw-session.json'
        rejected_path.write_text(json.dumps(rejected))
        before = page.evaluate(EXTRA)
        with page.expect_file_chooser() as chooser:
            page.get_by_role('button', name='Load Session', exact=True).first.click()
        chooser.value.set_files(str(rejected_path))
        # A rejected Load reports through the error log, not an alert.
        rejected_seen = poll(page, """async () => { const { state } = await import('./js/state.js');
            return !window.__GBDRAW_APP__.sessionImportPending && state.errorLog.value !== null; }""", 300, 0.5)
        record['rejectedLoad'] = {'reported': rejected_seen, 'before': before, 'after': page.evaluate(EXTRA),
                                  'error': page.evaluate("""async () => { const { state } = await import('./js/state.js');
                                      const e = state.errorLog.value; return e ? String(e.message || JSON.stringify(e)).slice(0, 300) : null; }""")}
        write()
        browser.close()


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
