const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { join, resolve } = require('node:path');
const { openApp, evaluateWithRetainedPromise } = require('./helpers/app-lifecycle.cjs');

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
    events: []
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
});

// The Session title is asked in the app's dialog, never a native prompt().
const titleDialog = (page, name) => page.locator('[data-choice-prompt]').getByRole('dialog', { name, exact: true });
const titleField = (dialog) => dialog.getByRole('textbox', { name: 'Session title', exact: true });
const recordNativeDialogs = (page) => {
  const seen = [];
  page.on('dialog', (dialog) => { seen.push(dialog.type()); void dialog.dismiss(); });
  return seen;
};

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
  const nativeDialogs = recordNativeDialogs(page);
  await openApp(page);
  await page.locator(
    'input[type="file"][accept*="application/json"][accept*="application/gzip"]'
  ).setInputFiles(sessionFile);
  await page.waitForFunction(() => (
    window.__GBDRAW_APP__?.sessionImportPending === false
    && window.__GBDRAW_APP__?.results?.length === 1
  ));

  const saveButton = page.getByRole('banner').getByRole('button', { name: 'Save Session', exact: true });
  const saveStatus = page.locator('[data-session-save-status]');
  await captureSaveInvariant(page);
  await resetSaveProbe(page);
  await installCompressionGate(page);

  const downloadPromise = page.waitForEvent('download');
  const start = await page.evaluate(() => {
    window.__GBDRAW_APP__.sessionTitle = '';
    const first = window.__GBDRAW_APP__.saveSessionWithTitle();
    const second = window.__GBDRAW_APP__.saveSessionWithTitle();
    window.__GBDRAW_SAVE_PROMISES__ = [first, second];
    return { samePromise: first === second };
  });
  expect(start.samePromise).toBe(true);
  // Save without a title asks for one once; Enter in the field saves.
  const saveTitle = titleDialog(page, 'Save Session');
  await expect(saveTitle).toHaveCount(1);
  await expect(titleField(saveTitle)).toBeFocused();
  await expect(titleField(saveTitle)).toHaveValue('');
  await titleField(saveTitle).fill('single-flight-save');
  await titleField(saveTitle).press('Enter');
  await expect(saveTitle).toHaveCount(0);

  await page.waitForFunction(() => window.__GBDRAW_SAVE_PROBE__.events.some(
    (event) => event.name === 'session-save-compression-start'
  ));
  await expect(saveButton).toBeDisabled();
  await expect(saveButton).toHaveAttribute('aria-busy', 'true');
  await expect(saveStatus).toBeVisible();
  await expect(saveStatus).toHaveText('Saving session…');

  await page.evaluate(() => window.__GBDRAW_SAVE_GATE__.release());
  const settled = await evaluateWithRetainedPromise(page, async () => {
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
  expect(await page.evaluate(() => window.__GBDRAW_APP__.sessionTitle)).toBe('single-flight-save');
  await expect(page.locator('[data-choice-prompt]')).toHaveCount(0);
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
  expect(eventNames(successProbe)).not.toContain('catalog.admission-started');
  expect(successProbe.events.find(
    ({ name }) => name === 'session-save-catalog-preparation-end'
  )?.reusedCommittedSession).toBe(true);
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

  // Escape and Cancel in the title dialog of a Save save nothing, and focus
  // returns to Save Session.
  for (const answer of ['Escape', 'Cancel']) {
    await resetSaveProbe(page);
    await page.evaluate(() => { window.__GBDRAW_APP__.sessionTitle = ''; });
    await saveButton.click();
    await expect(titleField(saveTitle)).toBeFocused();
    if (answer === 'Escape') await titleField(saveTitle).press('Escape');
    else await saveTitle.getByRole('button', { name: 'Cancel', exact: true }).click();
    await expect(saveTitle).toHaveCount(0);
    await expect(saveButton).toBeFocused();
    await expect.poll(async () => eventNames(await saveProbe(page))).toContain('session-save-pending-cleared');
    const canceledProbe = await saveProbe(page);
    expectOrderedEvents(canceledProbe, ['session-save-title-canceled', 'session-save-pending-cleared']);
    expect(eventNames(canceledProbe)).not.toContain('session-save-pending-published');
    expect(await page.evaluate(() => ({
      title: window.__GBDRAW_APP__.sessionTitle, pending: window.__GBDRAW_APP__.sessionSavePending
    }))).toEqual({ title: '', pending: false });
  }

  await resetSaveProbe(page);
  const compressionFailure = await evaluateWithRetainedPromise(page, async () => {
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
  expect(compressionFailure.result).toMatchObject({ status: 'error',
    error: { code: 'UNKNOWN', operation: 'session-save', stage: 'unknown' } });
  expect(compressionFailure.pending).toBe(false);
  expect(compressionFailure.errorSummary).toContain('without recognized diagnostic information');
  expect(compressionFailure.errorSummary).not.toContain('controlled compression failure');
  expectOrderedEvents(await saveProbe(page), [
    'session-save-compression-start',
    'session-save-error',
    'session-save-pending-cleared'
  ]);
  await expect(saveButton).toBeEnabled();

  await resetSaveProbe(page);
  // The second Save of the same file asks once; both calls share that question.
  await page.evaluate(() => {
    window.__GBDRAW_APP__.sessionTitle = 'single-flight-save';
    window.__GBDRAW_SAVE_PROMISES__ = [
      window.__GBDRAW_APP__.saveSessionWithTitle(),
      window.__GBDRAW_APP__.saveSessionWithTitle()
    ];
  });
  const again = page.getByRole('dialog', { name: 'Download the Session again?' });
  await expect(again).toHaveCount(1);
  await again.getByRole('button', { name: 'Cancel', exact: true }).click();
  await expect(again).toHaveCount(0);
  const repeatCanceled = await evaluateWithRetainedPromise(page, async () => {
    const results = await Promise.all(window.__GBDRAW_SAVE_PROMISES__);
    return {
      samePromise: window.__GBDRAW_SAVE_PROMISES__[0] === window.__GBDRAW_SAVE_PROMISES__[1],
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

  // A Session file over 50 MB asks before its download: Cancel
  // downloads nothing, Download saves. Every Blob reports 60 MB here.
  await page.evaluate(() => {
    const { get } = Object.getOwnPropertyDescriptor(Blob.prototype, 'size');
    Object.defineProperty(Blob.prototype, 'size', {
      configurable: true,
      get() { return Math.max(get.call(this), 60 * 1024 * 1024); }
    });
  });
  const large = page.getByRole('dialog', { name: 'Download a large Session?' });
  for (const choice of ['Cancel', 'Download']) {
    await resetSaveProbe(page);
    await page.evaluate(() => { window.__GBDRAW_SAVE_PROMISES__ = [window.__GBDRAW_APP__.saveSessionWithTitle()]; });
    await page.getByRole('dialog', { name: 'Download the Session again?' })
      .getByRole('button', { name: 'Download again', exact: true }).click();
    await expect(large).toContainText('The Session file is 60.0 MB.');
    const downloaded = choice === 'Download' ? page.waitForEvent('download') : null;
    await large.getByRole('button', { name: choice, exact: true }).click();
    await expect(large).toHaveCount(0);
    const status = await evaluateWithRetainedPromise(page, async () => (await window.__GBDRAW_SAVE_PROMISES__[0])?.status);
    expect(status).toBe(choice === 'Download' ? 'saved' : 'canceled');
    if (downloaded) await downloaded;
    const probe = await saveProbe(page);
    const canceled = probe.events.filter(({ name }) => name === 'session-save-download-canceled');
    expect(canceled.map(({ reason }) => reason)).toEqual(choice === 'Download' ? [] : ['large-download']);
  }
  expect(nativeDialogs).toEqual([]);
});

test('Rename Session asks in the app dialog with the title selected; Enter renames, Escape keeps it', async ({ page }) => {
  test.setTimeout(240_000);
  const nativeDialogs = recordNativeDialogs(page);
  await openApp(page);
  await page.locator(
    'input[type="file"][accept*="application/json"][accept*="application/gzip"]'
  ).setInputFiles(sessionFile);
  await page.waitForFunction(() => window.__GBDRAW_APP__?.sessionImportPending === false);
  const before = await page.evaluate(() => window.__GBDRAW_APP__.sessionTitle);
  expect(before.length).toBeGreaterThan(0);
  const opener = page.getByRole('banner').getByRole('button', { name: before, exact: true });
  const rename = titleDialog(page, 'Rename Session');
  const field = titleField(rename);

  await opener.click();
  await expect(field).toBeFocused();
  await expect(field).toHaveValue(before);
  // The whole title is selected, so typing replaces it.
  expect(await field.evaluate((input) => [input.selectionStart, input.selectionEnd])).toEqual([0, before.length]);
  await page.keyboard.press('Escape');
  await expect(rename).toHaveCount(0);
  await expect(opener).toBeFocused();
  expect(await page.evaluate(() => window.__GBDRAW_APP__.sessionTitle)).toBe(before);

  await opener.press('Enter');
  await expect(field).toBeFocused();
  await page.keyboard.type('  Renamed session  ');
  await page.keyboard.press('Enter');
  await expect(rename).toHaveCount(0);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.sessionTitle)).toBe('Renamed session');
  const renamed = page.getByRole('banner').getByRole('button', { name: 'Renamed session', exact: true });
  await expect(renamed).toBeFocused();

  // Rename names its outcome; Cancel keeps the title.
  await renamed.click();
  await field.fill('Not kept');
  await rename.getByRole('button', { name: 'Cancel', exact: true }).click();
  expect(await page.evaluate(() => window.__GBDRAW_APP__.sessionTitle)).toBe('Renamed session');
  await renamed.click();
  await field.fill('Kept');
  await rename.getByRole('button', { name: 'Rename', exact: true }).click();
  expect(await page.evaluate(() => window.__GBDRAW_APP__.sessionTitle)).toBe('Kept');
  expect(nativeDialogs).toEqual([]);
});

test('Rename Session is one History step: Undo and Redo move the title apart from the setting before it', async ({ page }) => {
  test.setTimeout(240_000);
  await openApp(page);
  await page.locator(
    'input[type="file"][accept*="application/json"][accept*="application/gzip"]'
  ).setInputFiles(sessionFile);
  await page.waitForFunction(() => window.__GBDRAW_APP__?.sessionImportPending === false);
  const read = () => page.evaluate(() => ({
    title: window.__GBDRAW_APP__.sessionTitle, gc: window.__GBDRAW_APP__.form.show_gc
  }));
  const history = (method) => evaluateWithRetainedPromise(page, (name) => window.__GBDRAW_HISTORY__[name](), method);
  const start = await read();
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.runUndoable('Change setting', () => {
    const form = window.__GBDRAW_APP__.form;
    form.show_gc = !form.show_gc;
  }));
  const changed = { ...start, gc: !start.gc };
  expect(await read()).toEqual(changed);

  await page.getByRole('banner').getByRole('button', { name: start.title, exact: true }).click();
  const field = titleField(titleDialog(page, 'Rename Session'));
  await expect(field).toBeFocused();
  await field.fill('Renamed');
  await field.press('Enter');
  expect(await read()).toEqual({ ...changed, title: 'Renamed' });
  expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.undoLabel())).toBe('Rename Session');

  await history('undo');
  expect(await read()).toEqual(changed);
  await history('undo');
  expect(await read()).toEqual(start);
  await history('redo');
  expect(await read()).toEqual(changed);
  await history('redo');
  expect(await read()).toEqual({ ...changed, title: 'Renamed' });
});
