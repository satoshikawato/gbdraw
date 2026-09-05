const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { join, resolve } = require('node:path');
const { gzipSync } = require('node:zlib');
const { openApp } = require('./helpers/app-lifecycle.cjs');

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

const installImportReadGate = async (page, filename) => {
  await page.evaluate((gatedFilename) => {
    const nativeStream = Blob.prototype.stream;
    let releaseGate;
    const gate = new Promise((resolveGate) => {
      releaseGate = resolveGate;
    });
    window.__GBDRAW_SESSION_IMPORT_GATE__ = {
      filename: gatedFilename,
      streamInvocations: 0,
      release: () => releaseGate()
    };
    File.prototype.stream = function gatedSessionStream() {
      if (this.name !== gatedFilename) return nativeStream.call(this);
      window.__GBDRAW_SESSION_IMPORT_GATE__.streamInvocations += 1;
      const source = nativeStream.call(this);
      return new ReadableStream({
        async start(controller) {
          await gate;
          const reader = source.getReader();
          try {
            while (true) {
              const { done, value } = await reader.read();
              if (done) break;
              controller.enqueue(value);
            }
            controller.close();
          } catch (error) {
            controller.error(error);
          } finally {
            reader.releaseLock();
          }
        }
      });
    };
  }, filename);
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
    'gzip-to-text-start'
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
  expect(dialogs).toHaveLength(1);
  expect(dialogs[0]).toMatch(/^Failed to load session:/);
  expect(await importSnapshot(page)).toEqual(before);

  await sessionInput.setInputFiles(invalidSession);
  await page.waitForFunction(() => (
    window.__GBDRAW_SESSION_LOADING_EVENTS__
      .filter(({ name }) => name === 'interactiveReady').length === 2
    && window.__GBDRAW_APP__?.sessionImportPending === false
  ));
  expect(dialogs).toHaveLength(2);
  expect(dialogs[1]).toMatch(/^Failed to load session:/);
  expect(await importSnapshot(page)).toEqual(before);
  await expect(sessionInput).toHaveValue('');
});
