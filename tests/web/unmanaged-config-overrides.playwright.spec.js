const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { join, resolve } = require('node:path');
const { gunzipSync } = require('node:zlib');
const {
  generateAndWaitForResult,
  openApp,
  waitForAppShell
} = require('./helpers/app-lifecycle.cjs');

const repoRoot = resolve(process.env.GBDRAW_REPO || process.cwd());
const sourceSession = JSON.parse(readFileSync(join(
  repoRoot,
  'gbdraw',
  'web',
  'gallery',
  'sessions',
  'HmmtDNA_basic_circular.gbdraw-session.json'
), 'utf8'));
const unmanagedPath = 'objects.gc_content.percent_background_opacity';

const loadSessionPayload = async (page, payload) => {
  const dialogPromise = page.waitForEvent('dialog', { timeout: 180000 });
  await page.locator('input[accept^=".json,"]').first().setInputFiles({
    name: 'unmanaged-config.gbdraw-session.json',
    mimeType: 'application/json',
    buffer: Buffer.from(JSON.stringify(payload))
  });
  const dialog = await dialogPromise;
  const message = dialog.message();
  await dialog.accept();
  return message;
};

const readSessionDownload = async (download) => {
  const path = await download.path();
  expect(path).toBeTruthy();
  return {
    path,
    session: JSON.parse(gunzipSync(readFileSync(path)).toString('utf8'))
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
  preservedSession.renderRequest.diagramOptions.configOverrides[unmanagedPath] = 0.42;
  preservedSession.renderRequest.diagramOptions.configOverrides[
    'objects.gc_content.mode'
  ] = 'percent';

  expect(await loadSessionPayload(page, preservedSession)).toBe(
    'Session loaded successfully!'
  );
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
  await page.evaluate(async () => window.__GBDRAW_APP__.saveSessionWithTitle());
  const saved = await readSessionDownload(await downloadPromise);
  expect(saved.session.version).toBe(40);
  expect(saved.session.config.unmanagedConfigOverrides).toEqual({
    [unmanagedPath]: 0.42
  });
  expect(
    saved.session.renderRequest.diagramOptions.configOverrides[unmanagedPath]
  ).toBe(0.42);
  expect(saved.session.config.form.plot_title).toBe('Preserved overlay journey');

  await page.reload({ waitUntil: 'domcontentloaded' });
  await waitForAppShell(page);
  const reloadDialog = page.waitForEvent('dialog', { timeout: 180000 });
  await page.locator('input[accept^=".json,"]').first().setInputFiles(saved.path);
  const loaded = await reloadDialog;
  expect(loaded.message()).toBe('Session loaded successfully!');
  await loaded.accept();
  await expect(disclosure).toContainText(unmanagedPath);

  await disclosure.getByRole('button', { name: `Reset preserved setting ${unmanagedPath}` })
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
  rejectedSession.renderRequest.diagramOptions.configOverrides[
    'objects.gc_content.unknown'
  ] = 1;
  const rejectedMessage = await loadSessionPayload(page, rejectedSession);
  expect(rejectedMessage).toContain('Unknown config override path');
  expect(rejectedMessage).toContain('objects.gc_content.unknown');
  await expect(disclosure).toHaveCount(0);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.form.plot_title)).toBe(
    'Preserved overlay journey'
  );
});
