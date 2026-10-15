// Step 3a (gbdraw/web/CLAUDE.md, Runtime data flow; CW-04 and CW-06 in
// docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md): a Generate whose
// feature catalog build did not change keeps the catalog object the page
// admitted; every other Generate sends and admits it. Each case below is one
// invalidation row of the CW-06 table in the step 3a pull request.
const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { gunzipSync } = require('node:zlib');
const { join } = require('node:path');
const { openApp, evaluateWithRetainedPromise } = require('./helpers/app-lifecycle.cjs');
const { generate, switchMode, download } = require('./helpers/mode-transition.cjs');

const HMMT = readFileSync(join(__dirname, '../test_inputs/HmmtDNA.gbk'), 'utf8');
const SESSIONS = 'gbdraw/web/gallery/sessions';

const start = async (browser) => {
  const context = await browser.newContext({ viewport: { width: 1600, height: 1000 } });
  await context.grantPermissions(['clipboard-read', 'clipboard-write']);
  const page = await context.newPage();
  page.externalRequests = [];
  await context.route('**/*', (route) => {
    if (new URL(route.request().url()).hostname === '127.0.0.1') return route.continue();
    page.externalRequests.push(route.request().url());
    return route.abort();
  });
  await page.addInitScript(() => {
    window.__REUSE__ = { events: [], metrics: [] };
    window.__GBDRAW_TEST_HOOKS__ = {
      onSessionLifecycleEvent: (event) => window.__REUSE__.events.push(event),
      onStructuralMetric: (metric) => window.__REUSE__.metrics.push(metric)
    };
  });
  page.on('dialog', (dialog) => dialog.dismiss());
  await openApp(page);
  return page;
};

const loadSession = async (page, file) => {
  await page.locator('input[accept^=".json,"]').setInputFiles(file);
  await expect.poll(() => page.evaluate(() => !window.__GBDRAW_APP__.sessionImportPending
    && window.__GBDRAW_APP__.extractedFeatures.length > 0), { timeout: 180_000 }).toBe(true);
};

const mark = (page) => page.evaluate(() => ({
  events: window.__REUSE__.events.length,
  metrics: window.__REUSE__.metrics.length,
  undo: window.__GBDRAW_HISTORY__.getUndoCount()
}));

// What one render did: whether it named a held revision, whether Python sent
// the catalog, and whether the page admitted or kept it.
const since = (page, start) => page.evaluate((start) => {
  const events = window.__REUSE__.events.slice(start.events);
  const metrics = {};
  for (const { name, value } of window.__REUSE__.metrics.slice(start.metrics)) {
    metrics[name] = (metrics[name] || 0) + Number(value || 0);
  }
  const transport = events.filter((event) => event.name === 'result-transport-ready');
  return {
    names: events.map((event) => event.name),
    transports: transport.map((event) => ({ held: event.heldFeatureCatalogRevision, catalog: event.featureCatalog })),
    admission: metrics.featureCatalogAdmissionCount || 0,
    reuse: metrics.featureCatalogReuseCount || 0,
    undo: window.__GBDRAW_HISTORY__.getUndoCount() - start.undo
  };
}, start);

const observeGenerate = async (page) => {
  const before = await mark(page);
  await generate(page);
  return since(page, before);
};

const keepShown = (page, key) => page.evaluate(async (key) => {
  const { state } = await import('./js/state.js');
  window[key] = window.Vue.toRaw(state.featureCatalog.value);
  return JSON.stringify(window[key]);
}, key);
const shownIs = (page, key) => page.evaluate(async (key) => {
  const { state } = await import('./js/state.js');
  return Boolean(window[key]) && window.Vue.toRaw(state.featureCatalog.value) === window[key];
}, key);
const shownJson = (page) => page.evaluate(async () => {
  const { state } = await import('./js/state.js');
  return JSON.stringify(window.Vue.toRaw(state.featureCatalog.value));
});
// One History step that sets app fields (`adv.*`, `form.*`) or renames a group.
const undoable = (page, label, { fields = {}, rename = null }) => evaluateWithRetainedPromise(page, async ({ label, fields, rename }) => {
  const app = window.__GBDRAW_APP__;
  await window.__GBDRAW_HISTORY__.runUndoable(label, () => {
    for (const [path, value] of Object.entries(fields)) {
      const [owner, key] = path.split('.');
      app[owner][key] = value;
    }
    if (rename) app.setOrthogroupNameOverride(rename.id, rename.name);
  });
  await window.Vue.nextTick();
}, { label, fields, rename });
const history = (page, action) => evaluateWithRetainedPromise(page, async (action) => {
  await window.__GBDRAW_HISTORY__[action]();
  await window.Vue.nextTick();
}, action);

const SENT = { held: false, catalog: 'sent' };
const HELD_SENT = { held: true, catalog: 'sent' };
const KEPT = { held: true, catalog: 'unchanged' };

test('a Generate keeps the held catalog only while its build is current', async ({ browser }, testInfo) => {
  test.setTimeout(420_000);
  const page = await start(browser);
  try {
    await page.evaluate(async (content) => {
      const app = window.__GBDRAW_APP__;
      app.files.c_gb = new File([content], 'HmmtDNA.gbk', { type: 'text/plain' });
      // The title has its place from the start: moving it there moves the
      // records too, which builds the catalog items again.
      app.adv.plot_title_position = 'top';
      await window.Vue.nextTick();
    }, HMMT);

    console.log('STEP Pw-1', Date.now());
    // Pw-1: the first Generate holds nothing; the catalog is sent and admitted.
    const cold = await observeGenerate(page);
    expect(cold.transports).toEqual([SENT]);
    expect(cold.names).toContain('catalog.admission-started');
    expect(cold.admission).toBe(1);
    await keepShown(page, 'c1');

    console.log('STEP Pw-2', Date.now());
    // Pw-2: no change: Python reuses the build; the page keeps the object and adds no History step.
    const warm = await observeGenerate(page);
    expect(warm.transports).toEqual([KEPT]);
    expect(warm.names).toContain('catalog.reused');
    expect({ admission: warm.admission, reuse: warm.reuse, undo: warm.undo }).toEqual({ admission: 0, reuse: 1, undo: 0 });
    expect(await shownIs(page, 'c1')).toBe(true);

    console.log('STEP Pw-3', Date.now());
    // Pw-3: a render-only change (the plot title text) draws a new figure with the same catalog.
    await undoable(page, 'Set the plot title', { fields: { 'form.plot_title': 'Reuse title' } });
    const retitled = await observeGenerate(page);
    expect(retitled.transports).toEqual([KEPT]);
    expect(retitled.undo).toBe(1);
    expect(await shownIs(page, 'c1')).toBe(true);
    const keptJson = await shownJson(page);

    console.log('STEP Pw-7', Date.now());
    // Pw-7: a new Worker has a new epoch: the same request sends the catalog,
    // and the one the page kept equals a fresh Generate's (D-31).
    await page.evaluate(async () => {
      (await import('./js/services/diagram-generation.js')).disposeDiagramGenerationWorker();
    });
    const fresh = await observeGenerate(page);
    expect(fresh.transports).toEqual([HELD_SENT]);
    expect(fresh.admission).toBe(1);
    expect(await shownJson(page)).toBe(keptJson);
    const c1Json = await keepShown(page, 'c1');

    console.log('STEP Pw-5', Date.now());
    // Pw-5: a drawn-feature change sends a new catalog; Undo and Redo restore
    // each object; only the newest one's builds are current.
    await undoable(page, 'Set the block stroke width', { fields: { 'adv.block_stroke_width': 3 } });
    const changed = await observeGenerate(page);
    expect(changed.transports).toEqual([HELD_SENT]);
    expect(changed.admission).toBe(1);
    await keepShown(page, 'c2');
    await history(page, 'undo');
    expect(await shownIs(page, 'c1')).toBe(true);
    await history(page, 'redo');
    expect(await shownIs(page, 'c2')).toBe(true);
    const again = await observeGenerate(page);
    expect(again.transports).toEqual([KEPT]);
    expect(await shownIs(page, 'c2')).toBe(true);
    await history(page, 'undo');
    await history(page, 'undo');
    expect(await shownIs(page, 'c1')).toBe(true);
    const older = await observeGenerate(page);
    expect(older.transports).toEqual([HELD_SENT]);
    expect(older.admission).toBe(1);
    expect(await shownJson(page)).toBe(c1Json);

    console.log('STEP Pw-6', Date.now());
    // Pw-6: each mode shows its own slot; the Linear Generate evicted the
    // Circular builds, so the Circular catalog is sent again.
    await switchMode(page, 'linear');
    await evaluateWithRetainedPromise(page, async (content) => {
      const app = window.__GBDRAW_APP__;
      app.lInputType = 'gb';
      if (!app.linearSeqs.length) app.addLinearSeq();
      app.setLinearSeqPrimaryFile(0, 'gb', new File([content], 'linear.gbk', { type: 'text/plain', lastModified: 1 }));
      await app.refreshLinearRecordSelectors();
    }, HMMT);
    const linear = await observeGenerate(page);
    expect(linear.transports).toEqual([SENT]);
    await switchMode(page, 'circular');
    const back = await observeGenerate(page);
    expect(back.transports).toEqual([HELD_SENT]);
    expect(back.admission).toBe(1);

    console.log('STEP Pw-9', Date.now());
    // Pw-9: a failed Load rolls back to the same object, still current.
    await keepShown(page, 'c3');
    const rollback = await evaluateWithRetainedPromise(page, async () => {
      const config = await import('./js/services/config.js');
      const file = new File([await (await fetch('/gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json')).arrayBuffer()], 'replacement.json');
      const result = await config.importSession({ target: { files: [file], value: 'selected' } }, {
        beforePreviewMount() { throw new Error('Catalog reuse rollback probe'); }
      });
      return result.status;
    });
    expect(rollback).toBe('error');
    expect(await shownIs(page, 'c3')).toBe(true);
    const afterRollback = await observeGenerate(page);
    expect(afterRollback.transports).toEqual([KEPT]);
    expect(await shownIs(page, 'c3')).toBe(true);

    console.log('STEP Pw-8', Date.now());
    // Pw-8: Save Session writes the kept catalog; a loaded catalog has no revision.
    // An untitled Session asks for a title first, and this page dismisses dialogs.
    await page.evaluate(async () => {
      (await import('./js/state.js')).state.sessionTitle.value = 'reuse';
    });
    const bytes = await download(page, 'Save Session', testInfo.outputPath('reuse.gbdraw-session.json.gz'));
    const saved = JSON.parse((bytes[0] === 0x1f ? gunzipSync(bytes) : bytes).toString('utf8'));
    expect(JSON.stringify(saved.editorState.featureCatalog)).toBe(await shownJson(page));
    await loadSession(page, testInfo.outputPath('reuse.gbdraw-session.json.gz'));
    const loaded = await observeGenerate(page);
    expect(loaded.transports).toEqual([SENT]);
    expect(loaded.admission).toBe(1);
    expect(page.externalRequests).toEqual([]);
  } finally {
    await page.context().close();
  }
});

test('a rerender of a loaded Session names the catalog under the names its Results keep', async ({ browser }) => {
  test.setTimeout(300_000);
  const page = await start(browser);
  try {
    await loadSession(page, `${SESSIONS}/HmmtDNA_basic_circular.gbdraw-session.json`);
    const savedName = await page.evaluate(() => window.__GBDRAW_APP__.results[0].name);
    const reflow = async () => {
      const before = await mark(page);
      await page.evaluate(async () => {
        const { state } = await import('./js/state.js');
        state.labelReflowForceRequestSeq.value += 1;
      });
      await expect.poll(async () => (await since(page, before)).transports.length, { timeout: 180_000 }).toBe(1);
      await expect.poll(() => page.evaluate(async () =>
        (await import('./js/state.js')).state.labelReflowProcessing.value), { timeout: 180_000 }).toBe(false);
      return since(page, before);
    };
    console.log('STEP Pw-4', Date.now());
    // Pw-4: the first reflow after Load holds nothing; Python names the items
    // with the engine's Result names, and the page keeps the saved names (OV-136).
    const first = await reflow();
    expect(first.transports).toEqual([SENT]);
    await keepShown(page, 'reflowed');
    expect(await page.evaluate(() => window.__GBDRAW_APP__.results[0].name)).toBe(savedName);
    const second = await reflow();
    expect(second.transports).toEqual([KEPT]);
    expect(second.names).toContain('catalog.reused');
    expect(await shownIs(page, 'reflowed')).toBe(true);
    expect(await page.evaluate(() => window.__GBDRAW_APP__.results[0].name)).toBe(savedName);
    // A Generate takes the engine's names: the held items do not match, so it is sent.
    const generated = await observeGenerate(page);
    expect(generated.transports).toEqual([SENT]);
    expect(generated.admission).toBe(1);
    const engineName = await page.evaluate(() => window.__GBDRAW_APP__.results[0].name);
    expect(engineName).not.toBe(savedName);
    expect(page.externalRequests).toEqual([]);
  } finally {
    await page.context().close();
  }
});

test('the similarity-group copy feedback clears when a Generate replaces the shown Results, as before', async ({ browser }) => {
  test.setTimeout(300_000);
  const page = await start(browser);
  try {
    await loadSession(page, `${SESSIONS}/majanivirus_orthogroup.gbdraw-session.json.gz`);
    const first = await observeGenerate(page);
    expect(first.transports).toEqual([SENT]);
    await page.evaluate(() => {
      const app = window.__GBDRAW_APP__;
      app.selectOrthogroup(app.orthogroups[0].id);
    });
    const label = () => page.evaluate(() => {
      const app = window.__GBDRAW_APP__;
      return app.orthogroupCopyFeedbackLabel(app.selectedOrthogroup, 'aa', 'Copy aa');
    });
    const copy = () => page.evaluate(async () => {
      const app = window.__GBDRAW_APP__;
      await app.copyOrthogroupSequences(app.selectedOrthogroup, 'aa');
      window.__REUSE__.copiedAt = performance.now();
      return app.orthogroupCopyFeedbackLabel(app.selectedOrthogroup, 'aa', 'Copy aa');
    });
    const elapsed = () => page.evaluate(() => performance.now() - window.__REUSE__.copiedAt);

    console.log('STEP Pw-10', Date.now());
    // Pw-10: copy while the reused Generate is in flight; its activation clears the label.
    await page.evaluate(() => {
      const app = window.__GBDRAW_APP__;
      const hooks = window.__GBDRAW_TEST_HOOKS__;
      const record = hooks.onSessionLifecycleEvent;
      window.__REUSE__.feedback = {};
      hooks.beforeDiagramGenerationResponse = async () => {
        delete hooks.beforeDiagramGenerationResponse;
        await app.copyOrthogroupSequences(app.selectedOrthogroup, 'aa');
        window.__REUSE__.copiedAt = performance.now();
        window.__REUSE__.feedback.before = app.orthogroupCopyFeedbackLabel(app.selectedOrthogroup, 'aa', 'Copy aa');
      };
      hooks.onSessionLifecycleEvent = (event) => {
        record(event);
        if (event.name === 'artifact.activation-completed' && window.__REUSE__.copiedAt) {
          window.__REUSE__.feedback.atActivation = app.orthogroupCopyFeedbackLabel(app.selectedOrthogroup, 'aa', 'Copy aa');
          window.__REUSE__.feedback.elapsed = performance.now() - window.__REUSE__.copiedAt;
        }
      };
    });
    const reused = await observeGenerate(page);
    expect(reused.names).toContain('catalog.reused');
    const feedback = await page.evaluate(() => window.__REUSE__.feedback);
    expect(feedback.elapsed, 'inconclusive: activation came after the label timer').toBeLessThan(1500);
    expect(feedback).toMatchObject({ before: 'Copied!', atActivation: 'Copy aa' });

    // An intent step (a group rename) crosses no Generate: Undo keeps the label.
    const groupId = await page.evaluate(() => window.__GBDRAW_APP__.selectedOrthogroup.id);
    await undoable(page, 'Rename the group', { rename: { id: groupId, name: 'Renamed group' } });
    expect(await copy()).toBe('Copied!');
    await history(page, 'undo');
    expect(await elapsed(), 'inconclusive: Undo came after the label timer').toBeLessThan(1500);
    expect(await label()).toBe('Copied!');

    // Undo of a reused Generate step replaces the shown Results: the label clears.
    await undoable(page, 'Set the plot title', { fields: { 'form.plot_title': 'Reuse title' } });
    const retitled = await observeGenerate(page);
    expect(retitled.transports).toEqual([KEPT]);
    expect(retitled.undo).toBe(1);
    expect(await copy()).toBe('Copied!');
    await history(page, 'undo');
    expect(await elapsed(), 'inconclusive: Undo came after the label timer').toBeLessThan(1500);
    expect(await label()).toBe('Copy aa');
    expect(page.externalRequests).toEqual([]);
  } finally {
    await page.context().close();
  }
});
