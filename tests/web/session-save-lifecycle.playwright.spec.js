const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { join, resolve } = require('node:path');
const { openApp } = require('./helpers/app-lifecycle.cjs');

const repoRoot = resolve(process.env.GBDRAW_REPO || process.cwd());
const sessionFile = {
  name: 'HmmtDNA_basic_circular.gbdraw-session.json',
  mimeType: 'application/json',
  buffer: readFileSync(join(
    repoRoot,
    'gbdraw',
    'web',
    'gallery',
    'sessions',
    'HmmtDNA_basic_circular.gbdraw-session.json'
  ))
};

const installSaveProbe = (page) => page.addInitScript(() => {
  window.__GBDRAW_SAVE_PROBE__ = {
    events: [],
    promptCalls: 0,
    confirmCalls: 0
  };
  window.__GBDRAW_TEST_HOOKS__ = {
    onSessionLifecycleEvent(event) {
      window.__GBDRAW_SAVE_PROBE__.events.push({ ...event });
    }
  };
});

const resetSaveProbe = (page) => page.evaluate(() => {
  const probe = window.__GBDRAW_SAVE_PROBE__;
  probe.events.length = 0;
  probe.promptCalls = 0;
  probe.confirmCalls = 0;
});

const saveProbe = (page) => page.evaluate(() => ({
  ...window.__GBDRAW_SAVE_PROBE__,
  events: window.__GBDRAW_SAVE_PROBE__.events.map((event) => ({ ...event }))
}));

const eventNames = (probe) => probe.events.map(({ name }) => name);

const expectOrderedEvents = (probe, expected) => {
  const names = eventNames(probe);
  let prior = -1;
  for (const name of expected) {
    const index = names.indexOf(name);
    expect(index, `${name} missing from ${names.join(', ')}`).toBeGreaterThan(prior);
    prior = index;
  }
};

const installCompressionGate = (page) => page.evaluate(() => {
  const NativeCompressionStream = window.CompressionStream;
  let release;
  const gate = new Promise((resolveGate) => { release = resolveGate; });
  window.__GBDRAW_SAVE_GATE__ = {
    constructorCalls: 0,
    release,
    restore() {
      window.CompressionStream = NativeCompressionStream;
    }
  };
  window.CompressionStream = class GatedCompressionStream {
    constructor(format) {
      window.__GBDRAW_SAVE_GATE__.constructorCalls += 1;
      const compression = new NativeCompressionStream(format);
      let gated = false;
      const inputGate = new TransformStream({
        async transform(chunk, controller) {
          if (!gated) {
            gated = true;
            await gate;
          }
          controller.enqueue(chunk);
        }
      });
      this.writable = inputGate.writable;
      this.readable = inputGate.readable.pipeThrough(compression);
    }
  };
});

const captureSaveInvariant = (page) => page.evaluate(async () => {
  const { state } = await import('/gbdraw/web/js/state.js');
  const app = window.__GBDRAW_APP__;
  const history = window.__GBDRAW_HISTORY__;
  window.__GBDRAW_SAVE_INVARIANT__ = {
    resultRefs: [...app.results],
    selectedResultIndex: app.selectedResultIndex,
    circularFile: state.files.c_gb,
    losatCache: state.losatCache.value,
    losatEntries: [...state.losatCache.value],
    proteinIdentityManifest: state.proteinIdentityManifest.value,
    undoCount: history.getUndoCount(),
    redoCount: history.getRedoCount()
  };
});

const saveInvariantIsIntact = (page) => page.evaluate(async () => {
  const { state } = await import('/gbdraw/web/js/state.js');
  const app = window.__GBDRAW_APP__;
  const history = window.__GBDRAW_HISTORY__;
  const before = window.__GBDRAW_SAVE_INVARIANT__;
  return {
    resultRefs: before.resultRefs.length === app.results.length
      && before.resultRefs.every((result, index) => result === app.results[index]),
    selectedResultIndex: before.selectedResultIndex === app.selectedResultIndex,
    circularFile: before.circularFile === state.files.c_gb,
    losatCache: before.losatCache === state.losatCache.value,
    losatEntries: JSON.stringify(before.losatEntries) === JSON.stringify([...state.losatCache.value]),
    proteinIdentityManifest: before.proteinIdentityManifest === state.proteinIdentityManifest.value,
    undoCount: before.undoCount === history.getUndoCount(),
    redoCount: before.redoCount === history.getRedoCount()
  };
});

test('Save Session is single-flight, paints pending state, and releases every settlement', async ({
  page
}) => {
  test.setTimeout(240_000);
  await installSaveProbe(page);
  page.on('dialog', (dialog) => dialog.accept());
  await openApp(page);
  await page.locator(
    'input[type="file"][accept*="application/json"][accept*="application/gzip"]'
  ).setInputFiles(sessionFile);
  await page.waitForFunction(() => (
    window.__GBDRAW_APP__?.sessionImportPending === false
    && window.__GBDRAW_APP__?.results?.length === 1
  ));

  const saveButton = page.getByRole('button', { name: 'Save Session', exact: true });
  const saveStatus = page.locator('[data-session-save-status]');
  await captureSaveInvariant(page);
  await resetSaveProbe(page);
  await installCompressionGate(page);

  const downloadPromise = page.waitForEvent('download');
  const start = await page.evaluate(() => {
    const probe = window.__GBDRAW_SAVE_PROBE__;
    window.__GBDRAW_APP__.sessionTitle = '';
    window.prompt = () => {
      probe.promptCalls += 1;
      return 'single-flight-save';
    };
    window.confirm = () => {
      probe.confirmCalls += 1;
      return true;
    };
    const first = window.__GBDRAW_APP__.saveSessionWithTitle();
    const second = window.__GBDRAW_APP__.saveSessionWithTitle();
    window.__GBDRAW_SAVE_PROMISES__ = [first, second];
    return { samePromise: first === second };
  });
  expect(start.samePromise).toBe(true);

  await page.waitForFunction(() => window.__GBDRAW_SAVE_PROBE__.events.some(
    (event) => event.name === 'session-save-compression-start'
  ));
  await expect(saveButton).toBeDisabled();
  await expect(saveButton).toHaveAttribute('aria-busy', 'true');
  await expect(saveStatus).toBeVisible();
  await expect(saveStatus).toHaveText('Saving session…');

  await page.evaluate(() => window.__GBDRAW_SAVE_GATE__.release());
  const settled = await page.evaluate(async () => {
    const [first, second] = await Promise.all(window.__GBDRAW_SAVE_PROMISES__);
    window.__GBDRAW_SAVE_GATE__.restore();
    return {
      statuses: [first?.status, second?.status],
      sameResult: first === second,
      compressionConstructors: window.__GBDRAW_SAVE_GATE__.constructorCalls
    };
  });
  expect(settled).toEqual({
    statuses: ['saved', 'saved'],
    sameResult: true,
    compressionConstructors: 1
  });
  await downloadPromise;

  const successProbe = await saveProbe(page);
  expect(successProbe.promptCalls).toBe(1);
  expect(successProbe.confirmCalls).toBe(0);
  expect(eventNames(successProbe).filter((name) => name === 'session-save-joined')).toHaveLength(1);
  expect(eventNames(successProbe).filter(
    (name) => name === 'session-save-projection-start'
  )).toHaveLength(1);
  expect(eventNames(successProbe).filter(
    (name) => name === 'session-save-compression-start'
  )).toHaveLength(1);
  expect(eventNames(successProbe).filter(
    (name) => name === 'session-save-download-handoff-completed'
  )).toHaveLength(1);
  expectOrderedEvents(successProbe, [
    'session-save-pending-published',
    'session-save-paint-opportunity-completed',
    'session-save-catalog-preparation-start',
    'session-save-projection-start',
    'session-save-projection-end',
    'session-save-compression-start',
    'session-save-compression-end',
    'session-save-download-handoff-completed',
    'session-save-pending-cleared'
  ]);
  await expect(saveButton).toBeEnabled();
  await expect(saveButton).toHaveAttribute('aria-busy', 'false');
  await expect(saveStatus).toHaveCount(0);
  expect(await saveInvariantIsIntact(page)).toEqual({
    resultRefs: true,
    selectedResultIndex: true,
    circularFile: true,
    losatCache: true,
    losatEntries: true,
    proteinIdentityManifest: true,
    undoCount: true,
    redoCount: true
  });

  await resetSaveProbe(page);
  const canceledTitle = await page.evaluate(async () => {
    window.__GBDRAW_APP__.sessionTitle = '';
    window.prompt = () => {
      window.__GBDRAW_SAVE_PROBE__.promptCalls += 1;
      return null;
    };
    const result = await window.__GBDRAW_APP__.saveSessionWithTitle();
    return { result, pending: window.__GBDRAW_APP__.sessionSavePending };
  });
  expect(canceledTitle).toEqual({ result: undefined, pending: false });
  expectOrderedEvents(await saveProbe(page), [
    'session-save-title-canceled',
    'session-save-pending-cleared'
  ]);

  await resetSaveProbe(page);
  const compressionFailure = await page.evaluate(async () => {
    const NativeCompressionStream = window.CompressionStream;
    window.__GBDRAW_APP__.sessionTitle = 'compression-error-save';
    window.CompressionStream = class FailingCompressionStream {
      constructor() {
        throw new Error('controlled compression failure');
      }
    };
    const result = await window.__GBDRAW_APP__.saveSessionWithTitle();
    window.CompressionStream = NativeCompressionStream;
    return {
      result,
      pending: window.__GBDRAW_APP__.sessionSavePending,
      errorSummary: String(window.__GBDRAW_APP__.errorLog?.summary || '')
    };
  });
  expect(compressionFailure.result).toEqual({ status: 'error' });
  expect(compressionFailure.pending).toBe(false);
  expect(compressionFailure.errorSummary).toContain('controlled compression failure');
  expectOrderedEvents(await saveProbe(page), [
    'session-save-compression-start',
    'session-save-error',
    'session-save-pending-cleared'
  ]);
  await expect(saveButton).toBeEnabled();

  await resetSaveProbe(page);
  const repeatCanceled = await page.evaluate(async () => {
    const probe = window.__GBDRAW_SAVE_PROBE__;
    window.__GBDRAW_APP__.sessionTitle = 'single-flight-save';
    window.confirm = () => {
      probe.confirmCalls += 1;
      return false;
    };
    const first = window.__GBDRAW_APP__.saveSessionWithTitle();
    const second = window.__GBDRAW_APP__.saveSessionWithTitle();
    const results = await Promise.all([first, second]);
    return {
      samePromise: first === second,
      sameResult: results[0] === results[1],
      statuses: results.map((result) => result?.status),
      pending: window.__GBDRAW_APP__.sessionSavePending
    };
  });
  expect(repeatCanceled).toEqual({
    samePromise: true,
    sameResult: true,
    statuses: ['canceled', 'canceled'],
    pending: false
  });
  const canceledProbe = await saveProbe(page);
  expect(canceledProbe.confirmCalls).toBe(1);
  expect(eventNames(canceledProbe).filter(
    (name) => name === 'session-save-download-canceled'
  )).toHaveLength(1);
  expect(eventNames(canceledProbe)).not.toContain('session-save-projection-start');
  expectOrderedEvents(canceledProbe, [
    'session-save-pending-published',
    'session-save-paint-opportunity-completed',
    'session-save-download-canceled',
    'session-save-pending-cleared'
  ]);
  await expect(saveButton).toBeEnabled();
});
