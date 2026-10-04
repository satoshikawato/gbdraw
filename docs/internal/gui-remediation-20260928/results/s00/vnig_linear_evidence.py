#!/usr/bin/env python
"""G03 evidence: load the Vnig Circular session, switch to Linear, add two GCF genomes, Generate.

Usage: python vnig_linear_evidence.py <dev|main> <v41|v44>
Writes <OUT>/cell-<app>-<session>.json and <app>-<session>-S2.png / -S4.png (+ dispatched request JSON).
Evidence only: drives real UI controls (Load Session file chooser, Linear mode button,
Linear GenBank upload zones + "Add sequence", "Generate Diagram"); reads state, never writes it.
"""
import json
import os
import pathlib
import sys
import time

from playwright.sync_api import sync_playwright

APPS = {'main': 'http://127.0.0.1:4301/gbdraw/web/index.html',
        'dev': 'http://127.0.0.1:4302/gbdraw/web/index.html'}
WT = str(pathlib.Path(__file__).resolve().parents[5])
OUT = pathlib.Path(os.environ.get('S00_BASELINE_DIR', '/home/kawato/gbdraw-baselines/gui-remediation-20260928')) / 'evidence/vnig-linear'
SESSIONS = {'v41': str(OUT / 'Vnig.v41.gbdraw-session.json.gz'),
            'v44': f'{WT}/gbdraw/web/gallery/sessions/Vnig_TUMSAT-TG-2018.gbdraw-session.json.gz'}
GENOMES = [f'{WT}/tests/test_inputs/GCF_000196095.1_ASM19609v1_genomic.gbff',
           f'{WT}/tests/test_inputs/GCF_000354175.2_ASM35417v2_genomic.gbff']

# Diagram-Worker tracking, adapted from tests/web/helpers/app-lifecycle.cjs (+ keeps the dispatched request).
WORKER_TRACKING = r"""
(() => {
  const activity = { constructions: 0, instances: [] };
  window.__GBDRAW_DIAGRAM_WORKER_ACTIVITY__ = activity;
  const NativeWorker = window.Worker;
  const nativePost = NativeWorker.prototype.postMessage;
  NativeWorker.prototype.postMessage = function (message, transfer) {
    const inst = this.__gbdrawLifecycleActivity;
    if (inst) {
      if (message?.type === 'init') { inst.initializations += 1; inst.events.push('init:request'); }
      else if (message?.type === 'helper') { inst.helpers.push(String(message.operation || '')); inst.events.push(`helper:request:${message.operation}`); }
      else if (message?.type === 'run') {
        let request = null;
        try { request = JSON.parse(JSON.stringify(message.payload?.request ?? null)); } catch (e) { request = { cloneError: String(e) }; }
        inst.runs.push({ requestId: String(message.requestId ?? ''), request });
        inst.events.push('run:request');
      }
    }
    return transfer === undefined ? nativePost.call(this, message) : nativePost.call(this, message, transfer);
  };
  window.Worker = new Proxy(NativeWorker, {
    construct(target, args) {
      const worker = Reflect.construct(target, args, target);
      const url = String(args[0] || '');
      if (!url.includes('diagram-generation-worker.js')) return worker;
      activity.constructions += 1;
      const inst = { id: activity.constructions, url, initializations: 0, helpers: [], runs: [], settlements: [], events: [] };
      activity.instances.push(inst);
      worker.__gbdrawLifecycleActivity = inst;
      worker.addEventListener('message', (event) => {
        const m = event.data || {};
        if (!['init', 'helper', 'run'].includes(m.type)) return;
        inst.settlements.push({ type: m.type, ok: m.ok === true, error: m.ok === false ? String(m.error?.message || m.error || '') : '' });
        inst.events.push(`${m.type}:${m.ok === true ? 'ok' : 'error'}`);
      });
      return worker;
    }
  });
})();
"""

SNAPSHOT = r"""
async ({ stage, app: appName }) => {
  const out = { stage };
  const safe = async (name, fn) => { try { out[name] = await fn(); } catch (e) { out[name] = { error: String(e?.message || e) }; } };
  const { state: s } = await import('./js/state.js');
  const app = window.__GBDRAW_APP__;
  const H = window.__GBDRAW_HISTORY__;
  const unref = (x) => (typeof x === 'function' ? x() : (x && typeof x === 'object' && 'value' in x ? x.value : x));
  const pick = (req) => {
    if (!req || typeof req !== 'object') return req ?? null;
    const co = req.diagramOptions?.configOverrides || {};
    return {
      mode: req.mode, grouping: req.grouping, recordCount: (req.records || []).length,
      records: (req.records || []).map((r) => ({ recordKey: r.recordKey, cardinality: r.cardinality, selector: r.selector, presentation: r.presentation })),
      layout: req.layout ?? null,
      plotTitle: req.diagramOptions?.plotTitle ?? null,
      plotTitlePosition: req.diagramOptions?.output?.plotTitlePosition ?? null,
      configOverrides: Object.fromEntries(Object.entries(co).filter(([k]) => /definition\.linear|replicon|accession|show_length|plot_title/.test(k)))
    };
  };
  await safe('app', () => ({
    mode: s.mode.value, processing: s.processing.value, sessionImportPending: unref(app.sessionImportPending),
    sessionTitle: unref(app.sessionTitle), resultGenerationKey: s.resultGenerationKey.value,
    resultsCount: s.results.value.length,
    errorLog: s.errorLog.value ? JSON.parse(JSON.stringify(s.errorLog.value)) : null
  }));
  await safe('title', () => ({
    stateField: 'form.plot_title', form_plot_title: s.form.plot_title,
    adv_plot_title_position: s.adv.plot_title_position, adv_plot_title_font_size: s.adv.plot_title_font_size,
    dom_PlotTitle_input: document.querySelector('input[aria-label="Plot Title"]')?.value ?? null,
    dom_PlotTitlePosition_select: document.querySelector('select[aria-label="Plot Title Position"]')?.value ?? null
  }));
  await safe('rows', () => {
    const enabled = Boolean(s.linearRecordLayoutEnabled.value);
    const entries = s.linearRecordRows.map((e) => ({ uid: e.uid, row: e.row }));
    const effective = s.linearSeqs.map((q, i) => ({ index: i, row: enabled ? (entries.find((e) => e.uid === q.uid)?.row ?? i + 1) : i + 1 }));
    const rowSet = effective.map((e) => e.row);
    const cb = document.querySelector('input[aria-label="Arrange linear records in rows"]');
    return {
      state_linearRecordLayoutEnabled: enabled, state_linearRecordRows: entries, state_linearRecordGap: s.linearRecordGap.value,
      dom_ArrangeInRows_checkbox: cb ? cb.checked : null,
      dom_rowInputs: [...document.querySelectorAll('input[aria-label^="Linear record row for sequence"]')].map((i) => i.value),
      effectiveRowsReplica: effective, effectiveHasSharedRow: new Set(rowSet).size !== rowSet.length
    };
  });
  await safe('labels', () => {
    const a = s.adv;
    const labelChecked = (text) => {
      const label = [...document.querySelectorAll('label')].find((l) => l.textContent.replace(/\s+/g, ' ').includes(text));
      return label?.querySelector('input[type=checkbox]')?.checked ?? null;
    };
    const sel = (id) => { const e = document.getElementById(id); return e ? { value: e.value, text: e.selectedOptions?.[0]?.textContent?.trim() } : null; };
    return {
      adv_linear_show_replicon: a.linear_show_replicon,
      adv_linear_accession_visibility: a.linear_accession_visibility, adv_linear_length_visibility: a.linear_length_visibility,
      adv_linear_show_accession: a.linear_show_accession, adv_linear_show_length: a.linear_show_length,
      dom_dev_Replicon_checkbox: document.querySelector('input[aria-label="Replicon visibility"]')?.checked ?? null,
      dom_dev_Accession_select: sel('linear-label-visibility-accession'),
      dom_dev_Length_select: sel('linear-label-visibility-length'),
      dom_dev_autoDisclosure: document.querySelector('[data-linear-label-auto-labels]')?.textContent?.trim() ?? null,
      dom_main_ShowReplicon: labelChecked('Show Replicon'), dom_main_ShowAccession: labelChecked('Show Accession'),
      dom_main_ShowLength: labelChecked('Show Length')
    };
  });
  await safe('linear', () => ({
    linearSeqs: s.linearSeqs.length,
    linearSourceGroups: unref(app.linearSourceGroups)?.map((g) => ({ records: g.records.length, gb: g.sequence?.gb?.name ?? null })) ?? null,
    seqs: s.linearSeqs.map((q) => ({ gb: q.gb?.name ?? (q.gb ? '<non-File>' : null), region_record_id: q.region_record_id, definition: q.definition, record_subtitle: q.record_subtitle }))
  }));
  await safe('history', () => ({ undo: H.getUndoCount(), redo: H.getRedoCount(), undoLabel: unref(H.undoLabel) ?? null, redoLabel: unref(H.redoLabel) ?? null }));
  await safe('worker', () => {
    const act = window.__GBDRAW_DIAGRAM_WORKER_ACTIVITY__ || { constructions: 0, instances: [] };
    return {
      constructions: act.constructions,
      initializations: act.instances.reduce((t, i) => t + i.initializations, 0),
      helpers: act.instances.flatMap((i) => i.helpers),
      runs: act.instances.reduce((t, i) => t + i.runs.length, 0),
      settlements: act.instances.flatMap((i) => i.settlements),
      dispatchedRequests: act.instances.flatMap((i) => i.runs.map((r) => pick(r.request)))
    };
  });
  const cfg = await import('./js/services/config.js');
  const sr = await import('./js/services/session-request.js');
  await safe('canonical_committed', () => ({ fn: 'services/config.js getCommittedCanonicalRenderRequest()', request: pick(cfg.getCommittedCanonicalRenderRequest()) }));
  await safe('canonical_draft', () => {
    if (typeof sr.projectGenerationIntent === 'function') {
      const r = sr.projectGenerationIntent({ state: s });
      return { fn: 'services/session-request.js projectGenerationIntent({state}).meaning', status: r.status, error: r.error ?? null, request: pick(r.meaning) };
    }
    const r = sr.buildCanonicalRenderRequest({ state: s, filesData: { ...s.files, linearSeqs: s.linearSeqs, linearComparisons: [] },
      comparisonPlanSnapshot: s.mode.value === 'linear' ? s.linearComparisonResolution.value : null });
    return { fn: 'services/session-request.js buildCanonicalRenderRequest({state, filesData:{...state.files, linearSeqs}, comparisonPlanSnapshot: state.linearComparisonResolution}).renderRequest', request: pick(r.renderRequest) };
  });
  await safe('svg', () => {
    const analyze = (content) => {
      if (!content) return null;
      const doc = new DOMParser().parseFromString(content, 'image/svg+xml');
      const texts = [...doc.querySelectorAll('text')].map((t) => ({ kind: t.getAttribute('data-definition-line-kind'), text: t.textContent.replace(/\s+/g, ' ').trim() })).filter((t) => t.text);
      const kinds = {}; texts.forEach((t) => { if (t.kind) kinds[t.kind] = (kinds[t.kind] || 0) + 1; });
      return {
        chars: content.length,
        recordGroups: [...doc.querySelectorAll('g[id^="record_group"], g[data-record-row], g[data-record-index]')].map((g) => ({ id: g.id, recordId: g.getAttribute('data-gbdraw-record-id'), index: g.getAttribute('data-record-index'), row: g.getAttribute('data-record-row'), column: g.getAttribute('data-record-column'), transform: g.getAttribute('transform') })),
        definitionLineKinds: kinds,
        accessionTexts: texts.filter((t) => /^NC_0(04603|04605|22349|22359)(\.1)?$/.test(t.text) || t.kind === 'accession').map((t) => t.text),
        lengthTexts: texts.filter((t) => /^(3,288,558|1,877,212|3,334,467|1,812,170) bp$/.test(t.text) || t.kind === 'length').map((t) => t.text),
        repliconTexts: texts.filter((t) => t.kind === 'replicon' || /^chromosome/i.test(t.text)).map((t) => t.text),
        titleTexts: texts.filter((t) => /nigripulchritudo|TUMSAT/.test(t.text)).map((t) => t.text),
        nameTexts: texts.filter((t) => t.kind === 'name' || /parahaemolyticus|alginolyticus/.test(t.text)).map((t) => t.text)
      };
    };
    const mounted = s.svgContainer.value?.querySelector('svg');
    const result = s.results.value[s.selectedResultIndex.value];
    return { mountedPreview: mounted ? analyze(mounted.outerHTML) : null, selectedResult: analyze(result?.content || '') };
  });
  return out;
}
"""


def main(app_name, session_key):
    cell = f'{app_name}-{session_key}'
    record = {'cell': cell, 'url': APPS[app_name], 'session': SESSIONS[session_key], 'genomes': GENOMES,
              'stages': {}, 'dialogs': [], 'pageErrors': [], 'consoleErrors': [], 'externalRequests': [], 'timings': {}}
    stage = {'name': 'S0'}

    def snap(page, name):
        page.evaluate('() => window.Vue?.nextTick?.()')
        record['stages'][name] = page.evaluate(SNAPSHOT, {'stage': name, 'app': app_name})
        record['stages'][name]['dialogsSoFar'] = list(record['dialogs'])
        (OUT / f'cell-{cell}.json').write_text(json.dumps(record, indent=1, ensure_ascii=False))

    def shots(page, name):  # full viewport + the two settings disclosures that hold the recorded controls
        page.screenshot(path=str(OUT / f'{cell}-{name}.png'), full_page=True)
        for tag, summary in [('titles', 'summary[aria-label="Titles and Record Labels"], summary[aria-label="Title & Legend"]'),
                             ('layout', 'summary[aria-label="Advanced comparison and layout"]')]:
            loc = page.locator(summary)
            if loc.count():
                loc.first.locator('xpath=..').screenshot(path=str(OUT / f'{cell}-{name}-{tag}.png'))

    def poll(page, js, timeout_s, interval=1.0):  # browser-side predicate polling
        end = time.time() + timeout_s
        while time.time() < end:
            if page.evaluate(js):
                return True
            time.sleep(interval)
        return False

    with sync_playwright() as p:
        browser = p.chromium.launch()
        context = browser.new_context(viewport={'width': 1440, 'height': 900})

        def route(r):
            if '127.0.0.1' in r.request.url:
                return r.continue_()
            record['externalRequests'].append(r.request.url)
            return r.abort()
        context.route('**/*', route)
        context.add_init_script(WORKER_TRACKING)
        page = context.new_page()
        page.on('pageerror', lambda e: record['pageErrors'].append({'stage': stage['name'], 'message': str(e)}))
        page.on('console', lambda m: m.type == 'error' and record['consoleErrors'].append({'stage': stage['name'], 'text': m.text}))

        def on_dialog(d):
            record['dialogs'].append({'stage': stage['name'], 'type': d.type, 'message': d.message})
            d.accept()
        page.on('dialog', on_dialog)

        t0 = time.time()
        page.goto(APPS[app_name], wait_until='domcontentloaded')
        page.wait_for_function('() => Object.keys(window.__GBDRAW_APP__?.paletteDefinitions || {}).length > 0', timeout=180_000)
        snap(page, 'S0')

        # S1: Load Session through the header button -> file chooser.
        stage['name'] = 'S1'
        with page.expect_file_chooser() as chooser:
            page.get_by_role('button', name='Load Session', exact=True).first.click()
        chooser.value.set_files(SESSIONS[session_key])
        t = time.time()
        while time.time() - t < 300 and not any(d['stage'] == 'S1' for d in record['dialogs']):
            if page.evaluate("""async () => { const { state } = await import('./js/state.js');
                    return !window.__GBDRAW_APP__.sessionImportPending && state.errorLog.value !== null; }"""):
                break
            time.sleep(1)
        time.sleep(2)  # let the success/failure alert and post-load watchers settle
        record['timings']['S1_load_s'] = round(time.time() - t, 1)
        snap(page, 'S1')
        loaded = any(d['stage'] == 'S1' and d['message'] == 'Session loaded successfully!' for d in record['dialogs'])
        record['sessionLoaded'] = loaded
        if not loaded:
            record['cellUnavailable'] = 'Session load did not succeed; S2-S4 not run.'
            page.screenshot(path=str(OUT / f'{cell}-S1-load-failed.png'), full_page=True)
            (OUT / f'cell-{cell}.json').write_text(json.dumps(record, indent=1, ensure_ascii=False))
            browser.close()
            return

        # S2: real "Linear" mode button.
        stage['name'] = 'S2'
        page.locator('button.app-mode-button[aria-label="Linear"]').click()
        poll(page, "() => window.__GBDRAW_APP__.mode === 'linear'", 30)
        time.sleep(1)
        snap(page, 'S2')
        # Open the relevant settings disclosures (UI-only) for the screenshots.
        for summary in ['summary[aria-label="Titles and Record Labels"]', 'summary[aria-label="Title & Legend"]',
                        'summary[aria-label="Advanced comparison and layout"]']:
            loc = page.locator(summary)
            if loc.count() and loc.first.evaluate('e => !e.parentElement.open'):
                loc.first.click()
        shots(page, 'S2')

        # S3: two GenBank files through the Linear upload zones (+ "Add sequence" for the second card).
        stage['name'] = 'S3'
        t = time.time()
        for index, genome in enumerate(GENOMES):
            if index > 0:
                page.get_by_role('button', name='Add sequence', exact=True).last.click()
                poll(page, f"() => window.__GBDRAW_APP__.linearSourceGroups.length > {index}", 30, 0.2)
            with page.expect_file_chooser() as chooser:
                page.locator('[data-linear-source-card]').nth(index).locator('.upload-zone').first.click()
            chooser.value.set_files(genome)
            name = pathlib.Path(genome).name
            ok = poll(page, f"""() => {{ const g = window.__GBDRAW_APP__.linearSourceGroups[{index}];
                return Boolean(g && g.sequence?.gb?.name === {json.dumps(name)} && g.records.length === 2); }}""", 300)
            record['timings'][f'S3_file{index + 1}_expanded_to_2_records'] = ok
        time.sleep(2)
        record['timings']['S3_upload_s'] = round(time.time() - t, 1)
        snap(page, 'S3')

        # S4: real "Generate Diagram" button.
        stage['name'] = 'S4'
        key = page.evaluate("async () => (await import('./js/state.js')).state.resultGenerationKey.value")
        t = time.time()
        page.get_by_role('button', name='Generate Diagram', exact=True).click()
        done = poll(page, f"""async () => {{ const {{ state }} = await import('./js/state.js');
            return !state.processing.value && (state.resultGenerationKey.value > {key} || state.errorLog.value !== null); }}""", 900, 5)
        record['timings']['S4_generate_s'] = round(time.time() - t, 1)
        record['timings']['S4_settled_within_15min'] = done
        time.sleep(2)
        snap(page, 'S4')
        shots(page, 'S4')
        svg_text = page.evaluate("async () => { const { state: s } = await import('./js/state.js'); return s.results.value[s.selectedResultIndex.value]?.content || ''; }")
        (OUT / f'{cell}-S4-result.svg').write_text(svg_text)
        dispatched = page.evaluate("""() => (window.__GBDRAW_DIAGRAM_WORKER_ACTIVITY__?.instances || []).flatMap(i => i.runs.map(r => r.request))""")
        (OUT / f'{cell}-dispatched-requests.json').write_text(json.dumps(dispatched, indent=1, ensure_ascii=False))
        record['timings']['total_s'] = round(time.time() - t0, 1)
        (OUT / f'cell-{cell}.json').write_text(json.dumps(record, indent=1, ensure_ascii=False))
        browser.close()


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
