const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { join, resolve } = require('node:path');
const { gunzipSync, gzipSync } = require('node:zlib');
const {
  openApp, expectNoLoadRepairDialog, expectSessionLoadedWithoutWarnings, generateAndWaitForResult, answerLoadRepairs
} = require('./helpers/app-lifecycle.cjs');
const { installImportReadGate } = require('./helpers/session-import-gate.cjs');
const { showDrawing } = require('./helpers/drawing-bar.cjs');

const repoRoot = resolve(process.env.GBDRAW_REPO || process.cwd());
const sessionInputSelector =
  'input[type="file"][accept*="application/json"][accept*="application/gzip"]';
const baselineSession = {
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
const replacementSession = {
  name: 'delayed-lambda.gbdraw-session.json.gz',
  mimeType: 'application/gzip',
  buffer: gzipSync(readFileSync(join(
    repoRoot,
    'gbdraw',
    'web',
    'gallery',
    'sessions',
    'lambda_basic_linear.gbdraw-session.json'
  )))
};

// OV-38: the 0.13.0 Gallery Session (version 30, release tag 0.13.0) saved its
// Circular Custom Track Slots rows with `spacing: null` while the slots were off.
const galleryV30Bytes = gunzipSync(readFileSync(join(
  repoRoot, 'tests', 'fixtures', 'sessions', 'BGC0000708-BGC0000713.v30.gbdraw-session.json.gz'
)));
const galleryV30Session = (edit = () => {}) => {
  const document = JSON.parse(galleryV30Bytes);
  edit(document.config.adv);
  return {
    name: 'BGC0000708-BGC0000713.gbdraw-session.json.gz',
    mimeType: 'application/gzip',
    buffer: gzipSync(Buffer.from(JSON.stringify(document)))
  };
};

// `whileLoading` runs while the Load is held before its import work.
const loadBaselineSession = async (page, whileLoading) => {
  if (whileLoading) await installImportReadGate(page, baselineSession.name);
  await page.locator(sessionInputSelector).setInputFiles(baselineSession);
  if (whileLoading) {
    await page.waitForFunction(() => window.__GBDRAW_SESSION_IMPORT_GATE__?.streamInvocations === 1);
    await whileLoading();
    await page.evaluate(() => window.__GBDRAW_SESSION_IMPORT_GATE__.release());
  }
  await expectSessionLoadedWithoutWarnings(page);
  await page.waitForFunction(() => (
    window.__GBDRAW_APP__?.sessionTitle === 'HmmtDNA_basic_circular'
  ));
  await page.waitForFunction(() => (
    window.__GBDRAW_APP__?.sessionImportPending === false
  ));
};

const installLifecycleProbe = async (page) => {
  await page.evaluate(() => {
    window.__GBDRAW_SESSION_LOADING_EVENTS__ = [];
    window.__GBDRAW_TEST_HOOKS__ = {
      onSessionLifecycleEvent(event) {
        window.__GBDRAW_SESSION_LOADING_EVENTS__.push({ ...event });
      }
    };
  });
};

const importSnapshot = (page) => page.evaluate(() => {
  const app = window.__GBDRAW_APP__;
  const selected = app.currentResults()[Number(app.selectedResultIndex) || 0] || null;
  const preview = document.querySelector('.shadow-xl.origin-top > svg');
  return {
    title: String(app.sessionTitle || ''),
    mode: String(app.mode || ''),
    selectedResultName: String(selected?.name || ''),
    selectedResultContent: String(selected?.content || ''),
    previewHtml: String(preview?.outerHTML || '')
  };
});

const expectLifecycleOrder = async (page, names) => {
  const events = await page.evaluate(() => (
    window.__GBDRAW_SESSION_LOADING_EVENTS__.map(({ name }) => name)
  ));
  let previous = -1;
  names.forEach((name) => {
    const current = events.indexOf(name);
    expect(current, `${name} missing from ${events.join(', ')}`).toBeGreaterThan(previous);
    previous = current;
  });
};

test('session loading is painted before import work and prevents duplicate adoption', async ({
  page
}) => {
  test.setTimeout(180_000);
  await openApp(page);
  // The Worker's boot files download after the Load finishes, not at the file
  // choice or while the Load runs.
  const bootUrls = new Set(await page.evaluate(async () => (
    await import('./js/services/pyodide-assets.js')).resolvePyodideBootAssetUrls()));
  /** @type {string[]} */
  const bootRequests = [];
  page.on('request', (request) => { if (bootUrls.has(request.url())) bootRequests.push(request.url()); });
  await loadBaselineSession(page, async () => expect(bootRequests).toEqual([]));
  await expect.poll(() => bootRequests.length).toBeGreaterThan(0);
  const before = await importSnapshot(page);
  expect(before.previewHtml).toContain('<svg');

  await installLifecycleProbe(page);
  await installImportReadGate(page, replacementSession.name);
  const dialogs = [];
  page.on('dialog', async (dialog) => {
    dialogs.push(dialog.message());
    await dialog.accept();
  });
  const notice = page.locator('[data-app-notice]');

  const sessionInput = page.locator(sessionInputSelector);
  const loadButton = page.getByRole('button', { name: 'Load Session', exact: true });
  const loadingStatus = page.locator('[data-session-import-status]');
  await sessionInput.setInputFiles(replacementSession);
  await page.waitForFunction(() => (
    window.__GBDRAW_SESSION_IMPORT_GATE__?.streamInvocations === 1
  ));

  await expect(loadingStatus).toBeVisible();
  await expect(loadingStatus).toHaveText('Loading session…');
  await expect(loadingStatus).toHaveAttribute('role', 'status');
  await expect(loadingStatus).toHaveAttribute('aria-live', 'polite');
  await expect(loadButton).toBeDisabled();
  await expect(loadButton).toHaveAttribute('aria-busy', 'true');
  expect(await importSnapshot(page)).toEqual(before);
  await expectLifecycleOrder(page, [
    'session-import-pending-published',
    'session-import-paint-opportunity-completed',
    'sessionSelection',
    'session-import-worker-start'
  ]);

  await sessionInput.setInputFiles({
    name: 'duplicate-while-pending.gbdraw-session.json',
    mimeType: 'application/json',
    buffer: Buffer.from('{"duplicate":true}')
  });
  expect(await importSnapshot(page)).toEqual(before);
  expect(await page.evaluate(() => (
    window.__GBDRAW_SESSION_LOADING_EVENTS__
      .filter(({ name }) => name === 'sessionSelection').length
  ))).toBe(1);
  expect(await page.evaluate(() => (
    window.__GBDRAW_SESSION_LOADING_EVENTS__
      .filter(({ name }) => name === 'firstCommittedPreview').length
  ))).toBe(0);
  expect(dialogs).toEqual([]);
  await expect(page.locator('[data-app-notice]', { hasText: 'Session loaded' })).toHaveCount(0);

  await page.evaluate(() => window.__GBDRAW_SESSION_IMPORT_GATE__.release());
  await page.waitForFunction(() => window.__GBDRAW_APP__?.sessionImportPending === false);
  await expect(loadingStatus).toBeHidden();
  await expect(loadButton).toBeEnabled();
  await expect(loadButton).toHaveAttribute('aria-busy', 'false');
  await expect(sessionInput).toHaveValue('');
  expect(dialogs).toEqual([]);
  await expect(notice).toContainText(`Session loaded: ${replacementSession.name}`, { timeout: 120_000 });
  await expect(page.locator('[data-app-notice] li')).toHaveCount(0);
  const after = await importSnapshot(page);
  expect(after.title).toBe('lambda_basic_linear');
  expect(after.mode).toBe('linear');
  expect(after.selectedResultName).toBe('lambda_basic_linear');
  expect(after.selectedResultContent).not.toBe(before.selectedResultContent);
  expect(after.previewHtml).toContain('<svg');
  expect(await page.evaluate(() => (
    window.__GBDRAW_SESSION_LOADING_EVENTS__
      .filter(({ name }) => name === 'firstCommittedPreview').length
  ))).toBe(1);

  await sessionInput.setInputFiles(replacementSession);
  await page.waitForFunction(() => (
    window.__GBDRAW_SESSION_LOADING_EVENTS__
      .filter(({ name }) => name === 'sessionSelection').length === 2
    && window.__GBDRAW_APP__?.sessionImportPending === false
  ));
  expect(dialogs).toEqual([]);
  await expect(notice).toContainText(`Session loaded: ${replacementSession.name}`, { timeout: 120_000 });
  await expect(page.locator('[data-app-notice] li')).toHaveCount(0);
  await expect(sessionInput).toHaveValue('');
});

// OV-40: the Session 39 Web writer stored a Label whitelist keyword typed with a
// tab as an extra cell (tests/fixtures/sessions/whitelist-tab-keyword.provenance.json).
// BA-D25: Load lists the row it reads as the current writer writes it in one
// dialog before it changes anything. Cancel leaves the loaded Session, its
// notice, and History as they were; [Load] then loads the file.
test('a lossy Session asks before it loads: Cancel changes nothing, Load loads it', async ({
  page
}) => {
  test.setTimeout(180_000);
  await openApp(page);
  await loadBaselineSession(page);
  const dialogs = [];
  page.on('dialog', async (dialog) => {
    dialogs.push(dialog.message());
    await dialog.accept();
  });
  const before = await importSnapshot(page);
  const visible = () => page.evaluate(() => ({
    notice: document.querySelector('[data-app-notice]')?.textContent?.trim() || '',
    filterMode: window.__GBDRAW_APP__.filterMode,
    canUndo: Boolean(window.__GBDRAW_APP__.canUndoHistory)
  }));
  const shownBefore = await visible();
  const lossy = join(repoRoot, 'tests', 'fixtures', 'sessions', 'whitelist-tab-keyword.v39.gbdraw-session.json.gz');
  // Its whitelist table has a row with extra cells; its feature information
  // is restored (OV-751), so the dialog lists the table row only.
  const repairLines = ['Label whitelist: line 1 had extra cells, joined into the last column with one space.'];
  await page.locator(sessionInputSelector).setInputFiles(lossy);
  await answerLoadRepairs(page, repairLines, 'cancel');
  await page.waitForFunction(() => window.__GBDRAW_APP__?.sessionImportPending === false);
  // The Cancel click leaves the pointer over the preview: a hover leaves an
  // empty style attribute behind, which is no change of the Result.
  await page.mouse.move(0, 0);
  await expect.poll(() => page.locator('[data-gbdraw-hover-opacity]').count()).toBe(0);
  const after = await importSnapshot(page);
  expect({ ...after, previewHtml: after.previewHtml.replaceAll(' style=""', '') }).toEqual(before);
  expect(await visible()).toEqual(shownBefore);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.errorLog)).toBeNull();

  await installLifecycleProbe(page);
  await page.locator(sessionInputSelector).setInputFiles(lossy);
  await answerLoadRepairs(page, repairLines, 'load');
  await expectSessionLoadedWithoutWarnings(page);
  await page.waitForFunction(() => (
    window.__GBDRAW_SESSION_LOADING_EVENTS__.some(({ name }) => name === 'interactiveReady')
    && window.__GBDRAW_APP__?.sessionImportPending === false
  ), null, { timeout: 120_000 });
  expect(await page.evaluate(() => window.__GBDRAW_APP__.errorLog)).toBeNull();
  expect(dialogs).toEqual([]);
  expect(await page.evaluate(() => ({
    filterMode: window.__GBDRAW_APP__.filterMode,
    whitelist: JSON.parse(JSON.stringify(window.__GBDRAW_APP__.manualWhitelist))
  }))).toEqual({
    filterMode: 'Whitelist',
    whitelist: [{ feat: 'CDS', qual: 'product', key: 'cytochrome c oxidase' }]
  });
});

test('failed session loading clears pending state and preserves the prior session', async ({
  page
}) => {
  test.setTimeout(180_000);
  await openApp(page);
  await loadBaselineSession(page);
  const before = await importSnapshot(page);
  const invalidSession = {
    name: 'delayed-invalid.gbdraw-session.json',
    mimeType: 'application/json',
    buffer: Buffer.from('{not valid JSON')
  };

  await installLifecycleProbe(page);
  await installImportReadGate(page, invalidSession.name);
  const dialogs = [];
  page.on('dialog', async (dialog) => {
    dialogs.push(dialog.message());
    await dialog.accept();
  });

  const sessionInput = page.locator(sessionInputSelector);
  const loadButton = page.getByRole('button', { name: 'Load Session', exact: true });
  const loadingStatus = page.locator('[data-session-import-status]');
  await sessionInput.setInputFiles(invalidSession);
  await page.waitForFunction(() => (
    window.__GBDRAW_SESSION_IMPORT_GATE__?.streamInvocations === 1
  ));
  await expect(loadingStatus).toBeVisible();
  await expect(loadButton).toBeDisabled();
  expect(await importSnapshot(page)).toEqual(before);

  await page.evaluate(() => window.__GBDRAW_SESSION_IMPORT_GATE__.release());
  await page.waitForFunction(() => window.__GBDRAW_APP__?.sessionImportPending === false);
  await expect(loadingStatus).toBeHidden();
  await expect(loadButton).toBeEnabled();
  await expect(sessionInput).toHaveValue('');
  expect(dialogs).toEqual([]);
  await expect(page.locator('[data-app-notice]', { hasText: 'Session loaded' })).toHaveCount(0);
  await expect(page.getByRole('alert', { name: 'Session not loaded' })).toContainText(
    'This file is not a gbdraw Session: it is not valid JSON.'
  );
  expect(await page.evaluate(() => window.__GBDRAW_APP__.errorLog)).toMatchObject({
    // The Session import Worker's actual stage (X-01, SE-09).
    code: 'INPUT_INVALID', stage: 'parse', context: { field: 'schema', reason: 'JSON_FORMAT' }
  });
  expect(await importSnapshot(page)).toEqual(before);

  await sessionInput.setInputFiles(invalidSession);
  await page.waitForFunction(() => (
    window.__GBDRAW_SESSION_LOADING_EVENTS__
      .filter(({ name }) => name === 'interactiveReady').length === 2
    && window.__GBDRAW_APP__?.sessionImportPending === false
  ));
  expect(dialogs).toEqual([]);
  await expect(page.locator('[data-app-notice]', { hasText: 'Session loaded' })).toHaveCount(0);
  await expect(page.getByRole('alert', { name: 'Session not loaded' })).toContainText(
    'This file is not a gbdraw Session: it is not valid JSON.'
  );
  expect(await importSnapshot(page)).toEqual(before);
  await expect(sessionInput).toHaveValue('');
});

// OV-372: a Session load error raised in Circular named a file the other mode
// never saw; it stayed on screen after the switch to Linear.
test('a Session load error does not follow a switch to the other diagram mode', async ({ page }) => {
  test.setTimeout(120_000);
  await openApp(page);
  await page.locator(sessionInputSelector).setInputFiles({
    name: 'not-a-session.gbdraw-session.json', mimeType: 'application/json', buffer: Buffer.from('{not valid JSON')
  });
  await page.waitForFunction(() => window.__GBDRAW_APP__?.sessionImportPending === false
    && window.__GBDRAW_APP__?.errorLog !== null);
  await expect(page.locator('.result-pane [role="alert"]')).toHaveCount(1);
  await showDrawing(page, { kind: 'linear' });
  await page.waitForFunction(() => window.__GBDRAW_APP__?.mode === 'linear');
  await expect(page.locator('.result-pane [role="alert"]')).toHaveCount(0);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.errorLog)).toBeNull();
});

// A retired slot field is a name an older writer wrote: Load reads the slots
// without it and with no dialog (receipt 〔出さないもの〕; E-OD03).
test('a 0.13.0 Gallery Session loads, and with its Custom Track Slots turned on it loads without the obsolete field', async ({
  page
}) => {
  test.setTimeout(300_000);
  await openApp(page);
  const dialogs = [];
  page.on('dialog', async (dialog) => {
    dialogs.push(dialog.message());
    await dialog.accept();
  });
  const sessionInput = page.locator(sessionInputSelector);
  const loadAndSettle = async (session) => {
    await installLifecycleProbe(page);
    await sessionInput.setInputFiles(session);
    await page.waitForFunction(() => (
      window.__GBDRAW_SESSION_LOADING_EVENTS__.some(({ name }) => name === 'interactiveReady')
      && window.__GBDRAW_APP__?.sessionImportPending === false
    ), null, { timeout: 240_000 });
  };
  await loadAndSettle(galleryV30Session());
  expect(await page.evaluate(() => window.__GBDRAW_APP__.errorLog)).toBeNull();
  expect(dialogs).toEqual([]);
  await expect(page.locator('[data-app-notice]')).toContainText('Session loaded', { timeout: 120_000 });
  const loaded = await importSnapshot(page);
  expect(loaded).toMatchObject({ title: 'out', mode: 'linear', selectedResultName: 'out' });
  expect(loaded.previewHtml).toContain('<svg');

  await loadAndSettle(galleryV30Session((adv) => {
    adv.circular_track_slots_enabled = true;
  }));
  expect(dialogs).toEqual([]);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.errorLog)).toBeNull();
  await expect(page.locator('[data-app-notice]')).toContainText('Session loaded', { timeout: 120_000 });
  await expect(page.locator('[data-choice-prompt]')).toHaveCount(0);
  expect(await page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    const adv = state.drawingOf('circular').adv;
    return { enabled: adv.circular_track_slots_enabled, spacing: adv.circular_track_slots.some((slot) => Object.hasOwn(slot, 'spacing')) };
  })).toEqual({ enabled: true, spacing: false });
});

// OV-677 (BA-D31): release 0.13.0 kept a Result across a mode switch. A
// Session it saved after Generate and a switch loads with no dialog (the
// feature information is read from the inputs of the Result's mode), opens on
// the mode the Result was drawn in, and keeps the shown mode without a Result.
for (const [name, resultMode, shownMode] of [
  ['mode-switch-circular-result-linear-shown.v30.gbdraw-session.json.gz', 'circular', 'linear'],
  ['mode-switch-linear-result-circular-shown.v30.gbdraw-session.json.gz', 'linear', 'circular']
]) {
  test(`a 0.13.0 Session saved with ${shownMode} shown opens its ${resultMode} Result in ${resultMode} (OV-677)`, async ({ page }) => {
    test.setTimeout(240_000);
    await openApp(page);
    await page.locator(sessionInputSelector).setInputFiles(join(repoRoot, 'tests', 'fixtures', 'sessions', name));
    await expectNoLoadRepairDialog(page);
    await expectSessionLoadedWithoutWarnings(page);
    const shown = () => page.evaluate(() => ({
      mode: window.__GBDRAW_APP__.mode, results: (window.__GBDRAW_APP__.results || []).length
    }));
    expect(await shown()).toEqual({ mode: resultMode, results: 1 });
    expect(await page.evaluate((mode) => window.__GBDRAW_APP__.selectDrawing(mode), shownMode)).toEqual({ status: 'ok' });
    await expect.poll(shown).toEqual({ mode: shownMode, results: 0 });
    expect(await page.evaluate((mode) => window.__GBDRAW_APP__.selectDrawing(mode), resultMode)).toEqual({ status: 'ok' });
    await expect.poll(shown).toEqual({ mode: resultMode, results: 1 });
    expect(await page.evaluate(() => window.__GBDRAW_APP__.errorLog)).toBeNull();
  });
}

// UJ-09 (Owner 2026-10-05): Load Session replaces the work and clears History,
// so it asks first when History changed since the last Save or Load. Cancel
// keeps the work and History; with nothing to lose the file picker opens at once.
test('Load Session asks before it replaces work changed since the last Save or Load', async ({
  page
}) => {
  test.setTimeout(180_000);
  await openApp(page);
  await loadBaselineSession(page);
  const loadButton = page.getByRole('button', { name: 'Load Session', exact: true });
  const confirm = page.getByRole('dialog', { name: 'Replace the current work?', exact: true });
  const undoCount = () => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
  let pickers = 0;
  page.on('filechooser', () => { pickers += 1; });

  await loadButton.click();
  await expect.poll(() => pickers).toBe(1);
  await expect(confirm).toHaveCount(0);

  // An open control transaction asks only when its value changed: focus alone
  // (the transaction a text control begins on focus) is nothing to lose.
  await page.evaluate(async () => {
    window.__focusOnlyTx = await window.__GBDRAW_HISTORY__.begin('Focus only');
    await window.__GBDRAW_APP__.openSessionFilePicker();
  });
  await expect(confirm).toHaveCount(0);
  await expect.poll(() => pickers).toBe(2);
  await page.evaluate(() => {
    window.__GBDRAW_APP__.form.prefix = 'typed-not-committed';
    // The picker call waits for the choice, so it is not awaited here.
    void window.__GBDRAW_APP__.openSessionFilePicker();
  });
  await expect(confirm).toBeVisible();
  await page.keyboard.press('Escape');
  await expect(confirm).toHaveCount(0);
  await page.evaluate(() => {
    window.__GBDRAW_APP__.form.prefix = '';
    window.__GBDRAW_HISTORY__.cancel(window.__focusOnlyTx);
  });
  expect(await undoCount()).toBe(0);

  const prefix = page.locator('#output-prefix');
  await prefix.fill('unsaved-work');
  await prefix.press('Tab');
  await expect.poll(undoCount).toBe(1);
  await loadButton.click();
  await expect(confirm).toBeVisible();
  await expect(confirm).toHaveAttribute('aria-modal', 'true');
  await expect(confirm.getByRole('button', { name: 'Cancel', exact: true })).toBeFocused();
  await page.keyboard.press('Escape');
  await expect(confirm).toHaveCount(0);
  await expect(loadButton).toBeFocused();
  expect(pickers).toBe(2);
  expect(await undoCount()).toBe(1);
  await expect(prefix).toHaveValue('unsaved-work');

  await loadButton.click();
  const chooserPromise = page.waitForEvent('filechooser');
  await confirm.getByRole('button', { name: 'Load Session', exact: true }).click();
  const chooser = await chooserPromise;
  await chooser.setFiles(baselineSession);
  // The page shows the first Load's notice: wait for this Load to start.
  await expectNoLoadRepairDialog(page);
  await expectSessionLoadedWithoutWarnings(page);
  await expect.poll(undoCount).toBe(0);
  await expect(prefix).not.toHaveValue('unsaved-work');

  // A Save is the new point: Load asks again only after a later change.
  await prefix.fill('saved-work');
  await prefix.press('Tab');
  await expect.poll(undoCount).toBe(1);
  const download = page.waitForEvent('download');
  await page.getByRole('button', { name: 'Save Session', exact: true }).click();
  await download;
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.sessionSavePending)).toBe(false);
  const pickersBefore = pickers;
  await loadButton.click();
  await expect.poll(() => pickers).toBe(pickersBefore + 1);
  await expect(confirm).toHaveCount(0);
});

// UJ-06 (Owner 2026-10-08): where this origin serves the Gallery catalog
// (gbdraw.app), the empty state's Load an example lists every catalog entry and
// loads the chosen one's Session through the Load Session path and its UJ-09
// question. Without the catalog (the pip GUI) the button is absent.
const galleryCatalog = JSON.parse(readFileSync(
  join(repoRoot, 'gbdraw', 'web', 'gallery', 'examples.json'), 'utf8'
));
const plainTitle = (title) => title.replace(/<\/?i>/g, '');
const blockExternalRequests = async (page) => {
  const external = [];
  await page.context().route('**/*', (route) => {
    if (new URL(route.request().url()).hostname === '127.0.0.1') return route.continue();
    external.push(route.request().url());
    return route.abort();
  });
  return external;
};
const exampleControls = (page) => ({
  example: page.getByRole('button', { name: 'Load an example', exact: true }),
  header: page.locator('header.app-header').getByRole('button', { name: 'Load example', exact: true }),
  chooser: page.getByRole('dialog', { name: 'Choose an example', exact: true }),
  confirm: page.getByRole('dialog', { name: 'Replace the current work?', exact: true })
});
// Owner 2026-10-08: an example opens fitted. The rendered SVG lies inside the
// preview frame, at a zoom other than the 1.0 its Session stored.
const exampleFit = (page) => page.evaluate(() => {
  const app = window.__GBDRAW_APP__;
  const container = app.canvasContainerRef;
  const frame = container.getBoundingClientRect();
  const box = app.svgContainer.querySelector('svg').getBoundingClientRect();
  const left = frame.left + container.clientLeft;
  const top = frame.top + container.clientTop;
  const inside = box.left >= left - 1 && box.top >= top - 1
    && box.right <= left + container.clientWidth + 1 && box.bottom <= top + container.clientHeight + 1;
  return { inside, zoomChanged: app.zoom !== 1 };
});
const expectExampleFitted = (page) => expect.poll(() => exampleFit(page))
  .toEqual({ inside: true, zoomChanged: true });
const acceptLoadedSession = async (page, title) => {
  return async () => {
    await expectSessionLoadedWithoutWarnings(page);
    await page.waitForFunction((expected) => (
      window.__GBDRAW_APP__?.sessionImportPending === false
      && window.__GBDRAW_APP__.sessionTitle === expected
      && window.__GBDRAW_APP__.results.length === 1
    ), title);
  };
};

test('Load an example lists every Gallery example and loads the chosen Circular or Linear .gz Session', async ({ page }) => {
  test.setTimeout(240_000);
  await page.setViewportSize({ width: 1440, height: 900 });
  const external = await blockExternalRequests(page);
  // The Linear .gz entry answers with the gzip of the small lambda Session, so
  // the file-name route of .gz loading runs without a 6 MB download.
  const gzEntry = galleryCatalog.find((entry) => entry.id === 'majanivirus_orthogroup');
  expect(gzEntry.session.endsWith('.json.gz')).toBe(true);
  expect(gzEntry.tags).toContain('Linear');
  const gzRequests = [];
  await page.route(`**/gbdraw/web/gallery/${gzEntry.session.slice(2)}`, (route) => {
    gzRequests.push(route.request().url());
    return route.fulfill({ status: 200, contentType: 'application/gzip', body: replacementSession.buffer });
  });
  await openApp(page);
  const { example, chooser } = exampleControls(page);
  await example.click();
  await expect(chooser).toBeVisible();
  const items = chooser.locator('[data-gallery-example]');
  await expect(items).toHaveCount(galleryCatalog.length);
  expect(await items.evaluateAll((nodes) => nodes.map((node) => node.dataset.galleryExample)))
    .toEqual(galleryCatalog.map((entry) => entry.id));
  for (const [index, entry] of galleryCatalog.entries()) {
    await expect(items.nth(index)).toContainText(plainTitle(entry.title));
    await expect(items.nth(index)).toContainText(entry.tags.includes('Linear') ? 'Linear' : 'Circular');
  }
  await expect(chooser.locator('[data-gallery-example="tobacco-chloroplast"] i'))
    .toHaveText('Nicotiana tabacum');
  await expect(chooser.locator('[data-gallery-example="lambda_basic_linear"] img'))
    .toHaveAttribute('src', /\/gallery\/thumbnails\/lambda_basic_linear\.webp$/);

  const circularLoaded = await acceptLoadedSession(page, 'HmmtDNA_basic_circular');
  await chooser.locator('[data-gallery-example="HmmtDNA_basic_circular"]').click();
  await circularLoaded();
  expect(await page.evaluate(() => window.__GBDRAW_APP__.mode)).toBe('circular');
  await expectExampleFitted(page);

  await openApp(page);
  await example.click();
  const linearLoaded = await acceptLoadedSession(page, 'lambda_basic_linear');
  await chooser.locator(`[data-gallery-example="${gzEntry.id}"]`).click();
  await linearLoaded();
  expect(await page.evaluate(() => window.__GBDRAW_APP__.mode)).toBe('linear');
  await expectExampleFitted(page);
  expect(gzRequests).toHaveLength(1);
  expect(external).toEqual([]);
});

test('Load an example is absent where this origin serves no Gallery catalog', async ({ page }) => {
  test.setTimeout(180_000);
  let catalogRequested;
  const catalogRead = new Promise((resolve) => { catalogRequested = resolve; });
  await page.route('**/gallery/examples.json', (route) => {
    catalogRequested();
    return route.fulfill({ status: 404, body: 'not found' });
  });
  // The 404 is a network log line, not an application error.
  const diagnostics = await openApp(page, { checkErrors: false });
  await catalogRead;
  await page.evaluate(() => new Promise((resolve) => requestAnimationFrame(() => resolve())));
  expect(diagnostics.pageErrors).toEqual([]);
  await expect(page.getByText('No Circular Result yet. Configure settings and click Generate.')).toBeVisible();
  await expect(exampleControls(page).example).toHaveCount(0);
  await expect(exampleControls(page).header).toHaveCount(0);
  // #1383 (PD-OI-108): the header still opens gbdraw.app's Gallery; the card
  // about the examples' tutorials goes with the examples.
  await expect(page.locator('[data-gallery-link]')).toHaveAttribute('href', 'https://gbdraw.app/gallery/');
  await expect(page.locator('[data-tutorial-card]')).toHaveCount(0);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.errorDisplay)).toBeFalsy();
});

// #1383 (PD-OI-108): the empty Preview, the header, and each chooser row lead to
// the Gallery's tutorials, and the notice of a loaded example links to its own
// tutorial, found by the example's Gallery id. A Session loaded from a file
// offers none, and the link goes with the notice.
const expectNewTab = async (link, href) => {
  await expect(link).toHaveAttribute('href', href);
  await expect(link).toHaveAttribute('target', '_blank');
  await expect(link).toHaveAttribute('rel', 'noopener noreferrer');
};
const tutorialHref = (page, id) => new URL(`/gbdraw/web/gallery/${id}#Tutorial`, page.url()).href;

test('a loaded example links to its Gallery tutorial, and a Session from a file does not', async ({ page }) => {
  test.setTimeout(240_000);
  await page.setViewportSize({ width: 1440, height: 900 });
  const external = await blockExternalRequests(page);
  await openApp(page);
  await expectNewTab(page.locator('[data-gallery-link]'), 'https://gbdraw.app/gallery/');
  await expectNewTab(
    page.locator('[data-tutorial-card]').getByRole('link', { name: 'Open the Interactive Gallery', exact: true }),
    'https://gbdraw.app/gallery/'
  );
  const { example, chooser } = exampleControls(page);
  await example.click();
  const tutorials = chooser.locator('[data-gallery-tutorial]');
  await expect(tutorials).toHaveCount(galleryCatalog.length);
  await expectNewTab(
    chooser.locator('[data-gallery-example="lambda_basic_linear"]').locator('xpath=..').locator('[data-gallery-tutorial]'),
    tutorialHref(page, 'lambda_basic_linear')
  );

  const loaded = await acceptLoadedSession(page, 'HmmtDNA_basic_circular');
  await chooser.locator('[data-gallery-example="HmmtDNA_basic_circular"]').click();
  await loaded();
  const notice = page.locator('[data-app-notice]');
  const offer = notice.locator('[data-app-notice-link]');
  await expect(offer).toContainText('This diagram has a step-by-step tutorial.');
  await expectNewTab(offer.getByRole('link', { name: 'Open the tutorial', exact: true }), tutorialHref(page, 'HmmtDNA_basic_circular'));
  await expect(page.locator('[data-app-notice-announcer]'))
    .toHaveText(/^Session loaded.* This diagram has a step-by-step tutorial\.$/);
  await notice.getByRole('button', { name: 'Dismiss notice', exact: true }).click();
  await expect(notice).toHaveCount(0);

  await page.locator(sessionInputSelector).setInputFiles(baselineSession);
  await expectSessionLoadedWithoutWarnings(page);
  await expect(offer).toHaveCount(0);
  expect(external).toEqual([]);
});

test('an example whose Gallery tutorial is not published offers no tutorial', async ({ page }) => {
  test.setTimeout(240_000);
  await page.route('**/gallery/examples.json', (route) => route.fulfill({
    status: 200,
    contentType: 'application/json',
    body: JSON.stringify(galleryCatalog.map((entry) => (
      entry.id === 'HmmtDNA_basic_circular' ? { ...entry, tutorialStatus: 'draft' } : entry
    )))
  }));
  await openApp(page);
  const { example, chooser } = exampleControls(page);
  await example.click();
  const row = (id) => chooser.locator(`[data-gallery-example="${id}"]`).locator('xpath=..');
  await expect(row('lambda_basic_linear').locator('[data-gallery-tutorial]')).toHaveCount(1);
  await expect(row('HmmtDNA_basic_circular').locator('[data-gallery-tutorial]')).toHaveCount(0);
  const loaded = await acceptLoadedSession(page, 'HmmtDNA_basic_circular');
  await chooser.locator('[data-gallery-example="HmmtDNA_basic_circular"]').click();
  await loaded();
  await expect(page.locator('[data-app-notice-link]')).toHaveCount(0);
});

// #1383 (PD-OI-108): the one-row header gives way without narrowing the
// Session title: where the Gallery link's words do not fit it shows its cap
// alone, and where even the cap would narrow the title the file buttons show
// their icons alone too. Every control keeps its name, and with room the
// words return.
test('the header gives way by icons, never by narrowing the Session title', async ({ page }) => {
  test.setTimeout(180_000);
  await page.setViewportSize({ width: 1440, height: 900 });
  await openApp(page);
  const header = page.locator('header.app-header');
  const words = page.locator('[data-gallery-link] .app-gallery-link-label');
  const titleWidth = () => page.locator('.app-session-button').evaluate((node) => node.getBoundingClientRect().width);
  const fitOf = () => header.evaluate((node) => ['data-file-icons', 'data-gallery-icon', 'data-session-title-alone']
    .map((name) => node.hasAttribute(name)));
  await expect.poll(fitOf).toEqual([false, false, false]);
  await expect(words).toBeVisible();
  const roomyTitle = await titleWidth();
  // Widths are measured on the one-row header with the Circular/Linear buttons
  // in the drawing bar (PD-OI-090): the Gallery cap alone from about 944 px, the
  // file icons from about 808 px, the title alone from about 788 px (English).
  await page.setViewportSize({ width: 880, height: 900 });
  await expect.poll(fitOf).toEqual([false, true, false]);
  await expect(words).toBeHidden();
  await expect(page.getByRole('link', { name: 'Gallery & tutorials', exact: true })).toBeVisible();
  expect(Math.abs(await titleWidth() - roomyTitle)).toBeLessThanOrEqual(1);

  // Where the cap alone still narrows the title, the file buttons show their icons.
  await page.setViewportSize({ width: 799, height: 900 });
  await expect.poll(fitOf).toEqual([true, true, false]);
  await expect(header.locator('.app-file-actions .app-config-button').first()).toHaveCSS('width', '32px');
  await expect(header.locator('.app-toolbar')).toHaveCSS('column-gap', '4px');
  expect(Math.abs(await titleWidth() - roomyTitle)).toBeLessThanOrEqual(1);

  // TL-D06: where the icons and the closed-up gaps still narrow the title, its
  // "Session" caption hides, so the title takes that width.
  await page.setViewportSize({ width: 778, height: 768 });
  await expect.poll(fitOf).toEqual([true, true, true]);
  await expect(header.locator('.app-session-label')).toBeHidden();
  await expect(header.locator('.app-toolbar')).toHaveCSS('column-gap', '2px');

  // The two-row header shows the cap alone in the logo's row, so the toolbar's
  // rows keep their height.
  await page.setViewportSize({ width: 390, height: 844 });
  await expect(words).toBeHidden();
  await expect(page.locator('[data-gallery-link]')).toBeVisible();
  const [logo, cap] = await Promise.all(['.app-logo-title', '[data-gallery-link]']
    .map((selector) => page.locator(selector).evaluate((node) => node.getBoundingClientRect().toJSON())));
  expect(cap.top).toBeGreaterThanOrEqual(logo.top - 1);
  expect(cap.bottom).toBeLessThanOrEqual(logo.bottom + 1);
  expect(cap.left).toBeGreaterThan(logo.right);
  await page.setViewportSize({ width: 1920, height: 900 });
  await expect.poll(fitOf).toEqual([false, false, false]);
  await expect(words).toBeVisible();
  await expect(header.locator('.app-session-label')).toBeVisible();
});

// D-07 (Owner 2026-10-10): the header keeps Load example once a Result exists,
// so a second example loads through the same chooser and UJ-09 question.
test('the header Load example stays after an example loads and loads another', async ({ page }) => {
  test.setTimeout(240_000);
  await page.setViewportSize({ width: 1440, height: 900 });
  await openApp(page);
  const { example, header, chooser, confirm } = exampleControls(page);
  await expect(header).toBeVisible();
  await header.click();
  const circularLoaded = await acceptLoadedSession(page, 'HmmtDNA_basic_circular');
  await chooser.locator('[data-gallery-example="HmmtDNA_basic_circular"]').click();
  await circularLoaded();
  await expect(example).toHaveCount(0);
  await expect(header).toBeVisible();
  await expect(header).toBeEnabled();

  await header.click();
  await expect(chooser).toBeVisible();
  await page.keyboard.press('Escape');
  await expect(chooser).toHaveCount(0);
  await expect(header).toBeFocused();

  const prefix = page.locator('#output-prefix');
  await prefix.fill('unsaved-work');
  await prefix.press('Tab');
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(1);
  await header.click();
  await chooser.locator('[data-gallery-example="lambda_basic_linear"]').click();
  await expect(confirm).toBeVisible();
  const linearLoaded = await acceptLoadedSession(page, 'lambda_basic_linear');
  await confirm.getByRole('button', { name: 'Load example', exact: true }).click();
  await linearLoaded();
  expect(await page.evaluate(() => window.__GBDRAW_APP__.mode)).toBe('linear');
  await expectExampleFitted(page);
  await expect(header).toBeVisible();
});

test('the example chooser cancels with focus back on Load an example and keeps the UJ-09 question', async ({ page }) => {
  test.setTimeout(180_000);
  await openApp(page);
  const { example, chooser, confirm } = exampleControls(page);
  const first = chooser.locator('[data-gallery-example]').first();
  await example.click();
  await expect(chooser).toHaveAttribute('aria-modal', 'true');
  await expect(first).toBeFocused();
  await page.keyboard.press('Shift+Tab');
  await expect(chooser.getByRole('button', { name: 'Cancel', exact: true })).toBeFocused();
  await page.keyboard.press('Escape');
  await expect(chooser).toHaveCount(0);
  await expect(example).toBeFocused();
  await example.click();
  await chooser.getByRole('button', { name: 'Cancel', exact: true }).click();
  await expect(chooser).toHaveCount(0);
  await expect(example).toBeFocused();
  expect(await page.evaluate(() => window.__GBDRAW_APP__.results.length)).toBe(0);

  const prefix = page.locator('#output-prefix');
  await prefix.fill('unsaved-work');
  await prefix.press('Tab');
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(1);
  await example.click();
  await chooser.locator('[data-gallery-example="HmmtDNA_basic_circular"]').click();
  await expect(chooser).toHaveCount(0);
  await expect(confirm).toBeVisible();
  await confirm.getByRole('button', { name: 'Cancel', exact: true }).click();
  await expect(confirm).toHaveCount(0);
  await expect(example).toBeFocused();
  await expect(prefix).toHaveValue('unsaved-work');

  await example.click();
  await chooser.locator('[data-gallery-example="HmmtDNA_basic_circular"]').click();
  const loaded = await acceptLoadedSession(page, 'HmmtDNA_basic_circular');
  await confirm.getByRole('button', { name: 'Load example', exact: true }).click();
  await loaded();
  expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(0);
  await expect(example).toHaveCount(0);
});

// OV-377/OV-378: the first feature popup after a Session load waits for Python
// to match the color rules. That wait is not an edit: Save and Load stay
// available, and the top bar says the popup is being prepared instead of
// "Applying an edit".
test('the first feature popup after a Session load says it is preparing and keeps Save and Load available', async ({ page }) => {
  await openApp(page);
  await page.locator(sessionInputSelector).setInputFiles(join(
    repoRoot, 'gbdraw', 'web', 'gallery', 'sessions', 'BGC0000708-BGC0000713.gbdraw-session.json'
  ));
  await expect(page.locator('[data-app-notice]')).toContainText('Session loaded', { timeout: 120_000 });
  await page.waitForFunction(() => window.__GBDRAW_APP__?.sessionImportPending === false
    && (window.__GBDRAW_APP__?.results || []).length > 0);
  // The livE CDS of record 1; no ribbon lies under its centre.
  const feature = page.locator('.gbdraw-preview-surface svg path[id="fed46a3a6_record_1"]');
  await feature.click({ force: true });
  const preparing = page.locator('[data-feature-details-pending]');
  await expect(preparing).toHaveText('Preparing feature details…');
  await expect(page.locator('[data-session-busy-reason]')).toHaveCount(0);
  await expect(page.getByRole('button', { name: 'Save Session', exact: true })).toBeEnabled();
  await expect(page.getByRole('button', { name: 'Load Session', exact: true })).toBeEnabled();
  await expect(page.locator('.feature-popup')).toBeVisible({ timeout: 180_000 });
  await expect(preparing).toHaveCount(0);
  expect(await page.evaluate(() => Boolean(window.__GBDRAW_APP__.clickedPairwiseMatch))).toBe(false);
});

// The Load notice describes the loaded Session: Generate keeps it; a mode
// switch and Reset Settings, which replace what it describes, clear it. A
// live region that is always rendered announces it.
test('the Load notice stays across Generate and clears on a mode switch and Reset Settings', async ({ page }) => {
  test.setTimeout(240_000);
  await openApp(page);
  await loadBaselineSession(page);
  const notice = page.locator('[data-app-notice]');
  await expect(page.locator('[data-app-notice-announcer]')).toHaveText(/^Session loaded/);
  await generateAndWaitForResult(page);
  await expect(notice).toContainText('Session loaded');
  expect(await page.evaluate(() => window.__GBDRAW_APP__.selectDrawing('linear'))).toEqual({ status: 'ok' });
  await expect(notice).toHaveCount(0);
  await expect(page.locator('[data-app-notice-announcer]')).toHaveCount(1);

  expect(await page.evaluate(() => window.__GBDRAW_APP__.selectDrawing('circular'))).toEqual({ status: 'ok' });
  await loadBaselineSession(page);
  await page.evaluate(() => { window.__GBDRAW_APP__.requestResetSettings(); });
  const reset = page.getByRole('dialog', { name: 'Reset all settings?', exact: true });
  await expect(reset).toContainText('Linear comparisons, which turn off; their uploaded tables stay');
  await reset.getByRole('button', { name: 'Reset all settings', exact: true }).click();
  await expect(notice).toHaveCount(0);
});
