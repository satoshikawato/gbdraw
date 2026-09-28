// Reproduce Issue #619 follow-up boundaries using only public Gallery data.
// Run from the repository root: node docs/internal/issue-619-boundary-followup-20260928/SESSION_RESULTS/probe_s01_browser.cjs /tmp/s01-browser
const { spawn } = require('node:child_process');
const { createHash } = require('node:crypto');
const { mkdir, readFile, writeFile } = require('node:fs/promises');
const { join, resolve } = require('node:path');
const { gunzipSync } = require('node:zlib');
const { chromium } = require('@playwright/test');
const { openApp, getDiagramWorkerActivity } = require('../../../../tests/web/helpers/app-lifecycle.cjs');

const root = process.cwd();
const output = resolve(process.argv[2] || '/tmp/gbdraw-619-s01-browser');
const fixture = resolve('gbdraw/web/gallery/sessions/tobacco-chloroplast.gbdraw-session.json');
const port = 46319;
const hash = bytes => createHash('sha256').update(bytes).digest('hex');
const sleep = ms => new Promise(resolve => setTimeout(resolve, ms));

async function startServer() {
  const server = spawn('python3', ['-m', 'http.server', String(port), '--bind', '127.0.0.1'], {
    cwd: root, stdio: 'ignore'
  });
  for (let i = 0; i < 50; i += 1) {
    if (server.exitCode !== null) throw new Error(`HTTP server exited: ${server.exitCode}`);
    try { if ((await fetch(`http://127.0.0.1:${port}/gbdraw/web/index.html`)).ok) return server; }
    catch (_) { /* startup */ }
    await sleep(100);
  }
  server.kill();
  throw new Error('HTTP server did not start');
}

async function openLoaded(browser, file = fixture) {
  const context = await browser.newContext({ baseURL: `http://127.0.0.1:${port}`, acceptDownloads: true });
  const page = await context.newPage();
  const dialogs = [];
  page.on('dialog', async dialog => { dialogs.push(dialog.message()); await dialog.accept(dialog.type() === 'prompt' ? '' : undefined); });
  await page.addInitScript(() => { window.__GBDRAW_TEST_HOOKS__ = {}; });
  await openApp(page);
  await page.locator('input[type="file"][accept^=".json,"]').setInputFiles(file);
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending && window.__GBDRAW_APP__.results.length === 1, null, { timeout: 180000 });
  return { page, context, dialogs };
}

async function snapshot(page) {
  return page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    const config = await import('./js/services/config.js');
    const slot = state.adv.circular_track_slots.find(row => row.id === 'gc_content');
    const value = slot.width;
    const numberTag = candidate => Number.isNaN(candidate) ? 'NaN'
      : candidate === Infinity ? 'Infinity' : candidate === -Infinity ? '-Infinity' : candidate;
    const rawWidth = value && typeof value === 'object'
      ? { value: numberTag(value.value), unit: value.unit } : numberTag(value);
    const digest = async value => Array.from(new Uint8Array(await crypto.subtle.digest('SHA-256', new TextEncoder().encode(value))))
      .map(byte => byte.toString(16).padStart(2, '0')).join('');
    return {
      mode: state.mode.value, generatedMode: state.generatedMode.value,
      config: config.buildConfigData(), rawWidth,
      request: config.getCommittedCanonicalRenderRequest(),
      resultHashes: await Promise.all(state.results.value.map(item => digest(item.content))),
      history: [window.__GBDRAW_HISTORY__.getUndoCount(), window.__GBDRAW_HISTORY__.getRedoCount()],
      error: state.errorLog.value && { stage: state.errorLog.value.stage, code: state.errorLog.value.code },
      bindingPositions: state.adv.multi_record_positions.map(item => ({ selector: item.selector, row: item.row })),
      annotationTargets: state.annotationSets.map(set => set.annotations.map(item => ({ id: item.id, record: item.target?.record,
        key: item.metadata?._gbdraw_web_target_record_key ?? null })))
    };
  });
}

async function save(page, name) {
  const pending = page.waitForEvent('download', { timeout: 20000 }).catch(() => null);
  await page.getByRole('banner').getByRole('button', { name: 'Save Session', exact: true }).click();
  const download = await pending;
  if (!download) return { saved: false, error: (await snapshot(page)).error };
  const path = join(output, `${name}.gbdraw-session.json.gz`);
  await download.saveAs(path);
  const bytes = await readFile(path);
  return { saved: true, path, sha256: hash(bytes), document: JSON.parse(gunzipSync(bytes)) };
}

const summarize = snapshot => ({
  mode: snapshot.mode, generatedMode: snapshot.generatedMode,
  profileMode: snapshot.config.modeProfiles?.activeMode,
  width: snapshot.rawWidth, firstRuleColor: snapshot.config.rules?.[0]?.color, firstRuleValue: snapshot.config.rules?.[0]?.val, history: snapshot.history,
  requestHash: hash(JSON.stringify(snapshot.request)), resultHashes: snapshot.resultHashes,
  positions: snapshot.bindingPositions, annotationTargets: snapshot.annotationTargets,
  error: snapshot.error
});


async function modeCase(browser) {
  const { page, context } = await openLoaded(browser);
  try {
    await page.getByRole('button', { name: 'Linear', exact: true }).click();
    await page.waitForFunction(() => window.__GBDRAW_APP__.mode === 'linear');
    await page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      await window.__GBDRAW_HISTORY__.runUndoable('S01 pending Linear draft', () => { state.adv.scale_font_size = 19; });
    });
    const before = await snapshot(page);
    const saved = await save(page, 'mode-linear-draft');
    const fresh = await openLoaded(browser, saved.path);
    const restored = await snapshot(fresh.page);
    const resaved = await save(fresh.page, 'mode-linear-draft-resaved');
    await fresh.context.close();
    const summarizeSave = item => ({ sha256: item.sha256, uiMode: item.document.ui?.mode,
      configMode: item.document.config?.modeProfiles?.activeMode, renderMode: item.document.renderRequest?.mode,
      scaleFontSize: item.document.config?.adv?.scale_font_size,
      requestHash: hash(JSON.stringify(item.document.renderRequest)),
      resultHashes: item.document.results?.map(result => hash(result.content)),
      resources: Object.fromEntries(Object.entries(item.document.resources || {})
        .map(([key, value]) => [key, hash(Buffer.from(value.data || '', 'base64'))])) });
    return { before: summarize(before), restored: summarize(restored), saved: summarizeSave(saved), resaved: summarizeSave(resaved) };
  } finally { await context.close(); }
}
async function numericCases(browser) {
  const results = [];
  for (const kind of ['NaN', 'Infinity', 'typedInfinity', 'textInfinity']) {
    const { page, context } = await openLoaded(browser);
    try {
      const before = await snapshot(page);
      await page.evaluate(async kind => {
        const app = window.__GBDRAW_APP__;
        const row = app.adv.circular_track_slots.find(item => item.id === 'gc_content');
        const values = { NaN: NaN, Infinity, typedInfinity: { value: Infinity, unit: 'px' }, textInfinity: { value: 'Infinity', unit: 'px' } };
        await window.__GBDRAW_HISTORY__.runUndoable('S01 injected measure', () => app.updateCircularTrackSlotMeasure(row, 'width', values[kind]));
        await window.Vue.nextTick();
      }, kind);
      const after = await snapshot(page);
      const attempted = await save(page, `numeric-before-undo-${kind}`);
      await page.evaluate(() => window.__GBDRAW_HISTORY__.undo());
      const undone = await snapshot(page);
      await page.evaluate(() => window.__GBDRAW_HISTORY__.redo());
      const redone = await snapshot(page);
      const saved = await save(page, `numeric-${kind}`);
      results.push({ kind, before: summarize(before), after: summarize(after), undo: summarize(undone),
        redo: summarize(redone), beforeUndoSave: { saved: attempted.saved, error: attempted.error,
          width: attempted.document?.config?.adv?.circular_track_slots?.find(row => row.id === 'gc_content')?.width },
        save: { saved: saved.saved, error: saved.error,
          width: saved.document?.config?.adv?.circular_track_slots?.find(row => row.id === 'gc_content')?.width } });
    } finally { await context.close(); }
  }
  return results;
}

async function bindingCase(browser) {
  const { page, context } = await openLoaded(browser);
  try {
    const before = await snapshot(page);
    const beforeSave = await save(page, 'binding-before');
    await page.evaluate(() => {
      window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse = () => { throw new Error('S01 forced post-render failure'); };
    });
    await page.getByRole('button', { name: 'Generate Diagram', exact: true }).click();
    await page.waitForFunction(() => !window.__GBDRAW_APP__.processing, null, { timeout: 180000 });
    const failed = await snapshot(page);
    const failedSave = await save(page, 'binding-failed');
    const loaded = await openLoaded(browser, failedSave.path);
    const restored = await snapshot(loaded.page);
    const resaved = await save(loaded.page, 'binding-reloaded');
    await loaded.context.close();
    await page.evaluate(() => { delete window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse; });
    await page.getByRole('button', { name: 'Generate Diagram', exact: true }).click();
    await page.waitForFunction(() => !window.__GBDRAW_APP__.processing, null, { timeout: 180000 });
    const retry = await snapshot(page);
    const resources = document => Object.fromEntries(Object.entries(document.resources || {})
      .map(([key, value]) => [key, hash(Buffer.from(value.data || '', 'base64'))]));
    const savedSummary = saved => ({ sha256: saved.sha256, version: saved.document?.version,
      requestSchema: saved.document?.renderRequest?.schema,
      positions: saved.document?.config?.adv?.multi_record_positions,
      annotationTargets: saved.document?.config?.annotationSets?.map(set => set.annotations.map(item => ({ id: item.id,
        record: item.target?.record, key: item.metadata?._gbdraw_web_target_record_key ?? null }))),
      resources: resources(saved.document), requestHash: hash(JSON.stringify(saved.document?.renderRequest)),
      resultHashes: saved.document?.results?.map(item => hash(item.content)) });
    return { before: summarize(before), failed: summarize(failed), restored: summarize(restored), retry: summarize(retry),
      beforeSave: savedSummary(beforeSave), failedSave: savedSummary(failedSave), reSave: savedSummary(resaved),
      worker: await getDiagramWorkerActivity(page) };
  } finally { await context.close(); }
}

async function historyCase(browser) {
  const { page, context } = await openLoaded(browser);
  const cdp = await context.newCDPSession(page);
  try {
    await cdp.send('Profiler.enable');
    await cdp.send('Profiler.startPreciseCoverage', { callCount: true, detailed: true });
    const count = async () => {
      const coverage = (await cdp.send('Profiler.takePreciseCoverage')).result;
      return coverage.filter(entry => entry.url.endsWith('/app/rule-matching.js'))
        .flatMap(entry => entry.functions)
        .filter(entry => entry.functionName === 'prepare')
        .reduce((total, entry) => total + Math.max(...entry.ranges.map(range => range.count)), 0);
    };
    const workerBefore = await getDiagramWorkerActivity(page);
    const before = await snapshot(page);
    await page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      await window.__GBDRAW_HISTORY__.runUndoable('S01 config only', () => { state.adv.scale_interval = 12345; });
    });
    const afterEdit = await count();
    await page.evaluate(() => window.__GBDRAW_HISTORY__.undo());
    const afterUndo = await count();
    await page.evaluate(() => window.__GBDRAW_HISTORY__.redo());
    const afterRedo = await count();
    const workerAfter = await getDiagramWorkerActivity(page);
    const after = await snapshot(page);
    const ruleBefore = summarize(after);
    await page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      await window.__GBDRAW_HISTORY__.runUndoable('S01 rule color', () => { state.manualSpecificRules[0].color = '#ff0000'; });
    });
    const ruleEditCalls = await count();
    const ruleEdit = summarize(await snapshot(page));
    await page.evaluate(() => window.__GBDRAW_HISTORY__.undo());
    const ruleUndoCalls = await count();
    const ruleUndo = summarize(await snapshot(page));
    await page.evaluate(() => window.__GBDRAW_HISTORY__.redo());
    const ruleRedoCalls = await count();
    const ruleRedo = summarize(await snapshot(page));
    const workerAfterRules = await getDiagramWorkerActivity(page);
    await page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      await window.__GBDRAW_HISTORY__.runUndoable('S01 rule predicate', () => { state.manualSpecificRules[0].val = 'psaB'; });
    });
    const predicateEditCalls = await count();
    const predicateEdit = summarize(await snapshot(page));
    await page.evaluate(() => window.__GBDRAW_HISTORY__.undo());
    const predicateUndoCalls = await count();
    const predicateUndo = summarize(await snapshot(page));
    await page.evaluate(() => window.__GBDRAW_HISTORY__.redo());
    const predicateRedoCalls = await count();
    const predicateRedo = summarize(await snapshot(page));
    const workerAfterPredicate = await getDiagramWorkerActivity(page);
    await page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      await window.__GBDRAW_HISTORY__.runUndoable('S01 independent setting', () => { state.adv.scale_interval = 23456; });
    });
    const independent = summarize(await snapshot(page));
    return { before: summarize(before), after: summarize(after), prepareCalls: { edit: afterEdit, undo: afterUndo, redo: afterRedo, ruleEdit: ruleEditCalls, ruleUndo: ruleUndoCalls, ruleRedo: ruleRedoCalls, predicateEdit: predicateEditCalls, predicateUndo: predicateUndoCalls, predicateRedo: predicateRedoCalls },
      rule: { before: ruleBefore, edit: ruleEdit, undo: ruleUndo, redo: ruleRedo },
      predicate: { edit: predicateEdit, undo: predicateUndo, redo: predicateRedo, independent },
      workerBefore, workerAfter, workerAfterRules, workerAfterPredicate };
  } finally { await cdp.send('Profiler.stopPreciseCoverage'); await context.close(); }
}

async function main() {
  await mkdir(output, { recursive: true });
  const server = await startServer();
  const browser = await chromium.launch({ headless: true });
  try {
    const result = {
      fixture: { path: fixture, sha256: hash(await readFile(fixture)) },
      versions: { node: process.version, browser: browser.version() }
    };
    if (process.env.GBDRAW_S01_SCOPE === 'binding') {
      result.bindings = await bindingCase(browser);
    } else {
      result.mode = await modeCase(browser);
      result.numeric = await numericCases(browser);
      result.bindings = await bindingCase(browser);
      result.history = await historyCase(browser);
    }
    await writeFile(join(output, 'summary.json'), JSON.stringify(result, null, 2));
    console.log(JSON.stringify({ output: join(output, 'summary.json'),
      mode: result.mode && { saved: result.mode.saved.uiMode, before: result.mode.before.mode, restored: result.mode.restored.mode },
      numeric: result.numeric?.map(item => ({kind:item.kind,after:item.after.width,undo:item.undo.width,redo:item.redo.width,save:item.save})),
      binding: {before: result.bindings.before.positions,failed: result.bindings.failed.positions,restored: result.bindings.restored.positions,
        annotationBefore: result.bindings.before.annotationTargets, annotationFailed: result.bindings.failed.annotationTargets},
      history: result.history?.prepareCalls }, null, 2));
  } finally { await browser.close(); server.kill(); }
}

main().catch(error => { console.error(error); process.exitCode = 1; });
