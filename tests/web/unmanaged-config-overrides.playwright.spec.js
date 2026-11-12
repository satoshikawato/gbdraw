const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { join, resolve } = require('node:path');
const { pathToFileURL } = require('node:url');
const { readSessionFile } = require('./helpers/session-file.cjs');
const { shownSavedDrawing } = require('./helpers/saved-session.cjs');
const {
  CURRENT_SESSION_VERSION,
  evaluateWithRetainedPromise,
  generateAndWaitForResult,
  openApp,
  waitForAppShell,
  expectSessionLoadedWithoutWarnings,
} = require('./helpers/app-lifecycle.cjs');

const repoRoot = resolve(process.env.GBDRAW_REPO || process.cwd());
const sessionDocuments = () => import(pathToFileURL(join(__dirname, 'helpers/session-documents.mjs')).href);
const sessionBytesOf = async document => Buffer.from(await (await (await sessionDocuments()).sessionFileOf(document)).arrayBuffer());
let sourceSession;
test.beforeAll(async () => {
  sourceSession = await (await sessionDocuments()).readSavedSession(new Blob([readFileSync(join(
    repoRoot,
    'gbdraw',
    'web',
    'gallery',
    'sessions',
    'HmmtDNA_basic_circular.gbdraw'
  ))]));
});
const unmanagedPath = 'objects.gc_content.percent_background_opacity';

const loadSessionPayload = async (page, payload) => {
  await page.locator('input[accept*=".gbdraw"]').first().setInputFiles({
    name: 'unmanaged-config.gbdraw',
    mimeType: 'application/octet-stream',
    buffer: await sessionBytesOf(payload)
  });
  const notice = page.locator('[data-app-notice]');
  await expect(notice).toContainText('Session loaded', { timeout: 180000 });
  return notice.innerText();
};

const readSessionDownload = async (download) => {
  const path = await download.path();
  expect(path).toBeTruthy();
  return {
    path,
    session: await readSessionFile(path)
  };
};

test('GUI-unmanaged config survives disclosure, Generate, save/reload, reset, and rejects unknown paths', async ({
  page
}) => {
  test.setTimeout(420000);
  await page.addInitScript(() => {
    window.__GBDRAW_RUN_CONFIG_OVERRIDES__ = [];
    const nativePostMessage = Worker.prototype.postMessage;
    Worker.prototype.postMessage = function trackedConfigOverridePostMessage(message, transfer) {
      if (message?.type === 'run') {
        window.__GBDRAW_RUN_CONFIG_OVERRIDES__.push(structuredClone(
          message.payload?.request?.diagramOptions?.configOverrides || {}
        ));
      }
      if (transfer === undefined) return nativePostMessage.call(this, message);
      return nativePostMessage.call(this, message, transfer);
    };
  });
  await openApp(page);

  const preservedSession = structuredClone(sourceSession);
  shownSavedDrawing(preservedSession).renderRequest.diagramOptions.configOverrides[unmanagedPath] = 0.42;
  shownSavedDrawing(preservedSession).renderRequest.diagramOptions.configOverrides[
    'objects.gc_content.mode'
  ] = 'percent';

  expect(await loadSessionPayload(page, preservedSession)).toContain(
    'Session loaded: unmanaged-config.gbdraw'
  );
  await expect(page.locator('[data-app-notice] li')).toHaveCount(0);
  const disclosure = page.locator('[data-unmanaged-config-override]');
  await expect(disclosure).toHaveCount(1);
  await expect(disclosure).toContainText(unmanagedPath);
  await expect(disclosure).toContainText('0.42');

  await page.evaluate(() => {
    window.__GBDRAW_APP__.form.plot_title = 'Preserved overlay journey';
  });
  await generateAndWaitForResult(page);
  expect(await page.evaluate((path) => {
    const runs = window.__GBDRAW_RUN_CONFIG_OVERRIDES__;
    return runs.at(-1)?.[path];
  }, unmanagedPath)).toBe(0.42);

  const downloadPromise = page.waitForEvent('download', { timeout: 120000 });
  await evaluateWithRetainedPromise(page, async () => window.__GBDRAW_APP__.saveSessionWithTitle());
  const saved = await readSessionDownload(await downloadPromise);
  expect(saved.session.version).toBe(CURRENT_SESSION_VERSION);
  const savedDraft = shownSavedDrawing(saved.session).config;
  expect(savedDraft.unmanagedConfigOverrides).toEqual({
    [unmanagedPath]: 0.42
  });
  expect(
    shownSavedDrawing(saved.session).renderRequest.diagramOptions.configOverrides[unmanagedPath]
  ).toBe(0.42);
  expect(savedDraft.form.plot_title).toBe('Preserved overlay journey');

  await page.reload({ waitUntil: 'domcontentloaded' });
  await waitForAppShell(page);
  await page.locator('input[accept*=".gbdraw"]').first().setInputFiles(saved.path);
  await expectSessionLoadedWithoutWarnings(page, { timeout: 180000 });
  await expect(disclosure).toContainText(unmanagedPath);

  await disclosure.getByRole('button', { name: `Reset setting ${unmanagedPath}` })
    .click();
  await expect(disclosure).toHaveCount(0);
  await page.getByRole('button', { name: 'Undo' }).click();
  await expect(disclosure).toContainText(unmanagedPath);
  await page.getByRole('button', { name: 'Redo' }).click();
  await expect(disclosure).toHaveCount(0);

  await generateAndWaitForResult(page);
  expect(await page.evaluate((path) => {
    const runs = window.__GBDRAW_RUN_CONFIG_OVERRIDES__;
    return Object.prototype.hasOwnProperty.call(runs.at(-1) || {}, path);
  }, unmanagedPath)).toBe(false);

  const rejectedSession = structuredClone(sourceSession);
  shownSavedDrawing(rejectedSession).renderRequest.diagramOptions.configOverrides[
    'objects.gc_content.unknown'
  ] = 1;
  const outcome = await evaluateWithRetainedPromise(page, async raw => {
    const file = new File([Uint8Array.from(atob(raw), c => c.charCodeAt(0))], 'invalid-overrides.gbdraw');
    return await window.__GBDRAW_APP__.importSession({ target: { files: [file], value: 'selected' } });
  }, (await sessionBytesOf(rejectedSession)).toString('base64'));
  expect(outcome).toMatchObject({ status: 'error', error: {
    code: 'INPUT_INVALID', operation: 'validateConfigOverrides', stage: 'helper', context: { field: 'configOverrides', reason: 'UNKNOWN_CONFIG_PATH' }
  } });
  await expect(page.getByRole('alert', { name: 'Operation error' }))
    .toContainText('Remove unknown configuration overrides or use a supported setting.');
  expect(JSON.stringify(outcome.error)).not.toContain('objects.gc_content.unknown');
  await expect(disclosure).toHaveCount(0);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.form.plot_title)).toBe(
    'Preserved overlay journey'
  );
});
