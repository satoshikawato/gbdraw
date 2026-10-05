const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { join, resolve } = require('node:path');
const { gunzipSync, gzipSync } = require('node:zlib');
const { openApp } = require('./helpers/app-lifecycle.cjs');
const { installImportReadGate } = require('./helpers/session-import-gate.cjs');

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

const loadBaselineSession = async (page) => {
  const dialogPromise = page.waitForEvent('dialog');
  await page.locator(sessionInputSelector).setInputFiles(baselineSession);
  const dialog = await dialogPromise;
  expect(dialog.message()).toBe('Session loaded successfully!');
  await dialog.accept();
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
  const selected = app.results?.[Number(app.selectedResultIndex) || 0] || null;
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
  await loadBaselineSession(page);
  const before = await importSnapshot(page);
  expect(before.previewHtml).toContain('<svg');

  await installLifecycleProbe(page);
  await installImportReadGate(page, replacementSession.name);
  const dialogs = [];
  page.on('dialog', async (dialog) => {
    dialogs.push(dialog.message());
    await dialog.accept();
  });

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

  await page.evaluate(() => window.__GBDRAW_SESSION_IMPORT_GATE__.release());
  await page.waitForFunction(() => window.__GBDRAW_APP__?.sessionImportPending === false);
  await expect(loadingStatus).toBeHidden();
  await expect(loadButton).toBeEnabled();
  await expect(loadButton).toHaveAttribute('aria-busy', 'false');
  await expect(sessionInput).toHaveValue('');
  expect(dialogs).toEqual(['Session loaded successfully!']);
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
  expect(dialogs).toEqual([
    'Session loaded successfully!',
    'Session loaded successfully!'
  ]);
  await expect(sessionInput).toHaveValue('');
});

// OV-40: the Session 39 Web writer stored a Label whitelist keyword typed with a
// tab as an extra cell (tests/fixtures/sessions/whitelist-tab-keyword.provenance.json).
// Load reads the row as the current writer writes it and names it in the notice.
test('a Session 39 table row with extra cells loads with a notice that names it', async ({
  page
}) => {
  test.setTimeout(180_000);
  await openApp(page);
  const dialogs = [];
  page.on('dialog', async (dialog) => {
    dialogs.push(dialog.message());
    await dialog.accept();
  });
  await installLifecycleProbe(page);
  await page.locator(sessionInputSelector).setInputFiles(join(
    repoRoot, 'tests', 'fixtures', 'sessions', 'whitelist-tab-keyword.v39.gbdraw-session.json.gz'
  ));
  await page.waitForFunction(() => (
    window.__GBDRAW_SESSION_LOADING_EVENTS__.some(({ name }) => name === 'interactiveReady')
    && window.__GBDRAW_APP__?.sessionImportPending === false
  ), null, { timeout: 120_000 });
  expect(await page.evaluate(() => window.__GBDRAW_APP__.errorLog)).toBeNull();
  expect(dialogs).toEqual(['Session loaded successfully! Some table rows of this older Session were read '
    + 'as the current version writes them. Label whitelist: line 1 had extra cells, joined into the last column '
    + 'with one space.']);
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
  await expect(page.getByRole('alert', { name: 'Operation error' })).toContainText('Use valid JSON.');
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
  await expect(page.getByRole('alert', { name: 'Operation error' })).toContainText('Use valid JSON.');
  expect(await importSnapshot(page)).toEqual(before);
  await expect(sessionInput).toHaveValue('');
});

test('a 0.13.0 Gallery Session loads, and its Custom Track Slots turned on name the obsolete field', async ({
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
  expect(dialogs).toHaveLength(1);
  expect(dialogs[0]).toMatch(/^Session loaded successfully!/);
  const loaded = await importSnapshot(page);
  expect(loaded).toMatchObject({ title: 'out', mode: 'linear', selectedResultName: 'out' });
  expect(loaded.previewHtml).toContain('<svg');

  await loadAndSettle(galleryV30Session((adv) => {
    adv.circular_track_slots_enabled = true;
  }));
  expect(dialogs).toHaveLength(1);
  await expect(page.getByRole('alert', { name: 'Operation error' })).toContainText(
    'The track settings are invalid. Track row 1. Field: spacing. Custom Track Slots no longer read this field.'
  );
  expect(await page.evaluate(() => window.__GBDRAW_APP__.errorLog)).toMatchObject({
    code: 'TRACK_INVALID', context: { field: 'spacing', reason: 'OBSOLETE_TRACK_FIELD', slotIndex: 0 }
  });
  expect(await importSnapshot(page)).toEqual(loaded);
});
