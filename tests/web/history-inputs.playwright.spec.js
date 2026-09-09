const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { resolve } = require('node:path');
const { gunzipSync } = require('node:zlib');
const { openApp } = require('./helpers/app-lifecycle.cjs');

const seed = resolve('gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json');
const sessionInput = 'input[type="file"][accept^=".json,"]';

const preparePage = async (page) => {
  await page.addInitScript(() => {
    window.__HISTORY_INPUT_EVENTS__ = [];
    window.__GBDRAW_TEST_HOOKS__ = {
      onSessionLifecycleEvent: (event) => window.__HISTORY_INPUT_EVENTS__.push(event)
    };
  });
  page.on('dialog', (dialog) => dialog.accept());
  await openApp(page);
};

const loadSession = async (page, path) => {
  await page.locator(sessionInput).setInputFiles(path);
  await page.waitForFunction(() => (
    window.__HISTORY_INPUT_EVENTS__.some((event) => event.name === 'interactiveReady')
    && !window.__GBDRAW_APP__.sessionImportPending
  ));
};

const generate = async (page) => {
  const before = await page.evaluate(() => window.__HISTORY_INPUT_EVENTS__
    .filter((event) => event.name === 'generate.processing-cleared').length);
  await page.getByRole('button', { name: 'Generate Diagram', exact: true }).click();
  await page.waitForFunction((count) => window.__HISTORY_INPUT_EVENTS__
    .filter((event) => event.name === 'generate.processing-cleared').length > count, before,
  { timeout: 180_000 });
  expect(await page.evaluate(() => window.__GBDRAW_APP__.errorLog?.summary || '')).toBe('');
  return page.evaluate(() => window.__GBDRAW_APP__.results[0].content);
};

const expectHistory = async (page, undo, redo, label) => {
  await expect.poll(() => page.evaluate(() => {
    const history = window.__GBDRAW_HISTORY__;
    return [history.getUndoCount(), history.getRedoCount(), history.undoLabel()];
  })).toEqual([undo, redo, label]);
};

for (const inputMethod of ['keyboard', 'pointer']) {
  test(`Label Mode ${inputMethod} edit has one Undo step and survives generation and fresh Load @pr-smoke`, async ({
    page, context, browser
  }, testInfo) => {
    test.setTimeout(300_000);
    const externalRequests = [];
    const allowLocal = (route) => {
      if (new URL(route.request().url()).origin === 'http://127.0.0.1:4173') {
        return route.continue();
      }
      externalRequests.push(route.request().url());
      return route.abort();
    };
    await context.route('**/*', allowLocal);
    await preparePage(page);
    await loadSession(page, seed);
    const originalSvg = await generate(page);
    const baseline = await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
    await expectHistory(page, baseline, 0, 'Generate diagram');

    const summary = page.locator('summary[aria-label="Labels"]');
    await summary.click();
    const select = page.locator('#circular-label-mode');
    await expect(select).toHaveValue('out');
    if (inputMethod === 'keyboard') {
      await summary.focus();
      await page.keyboard.press('Tab');
      await expect(select).toBeFocused();
      await page.keyboard.press('ArrowDown');
    } else {
      // Open the native popup through its pointer admission path.
      await select.click();
      await page.keyboard.press('ArrowDown');
      await page.keyboard.press('Enter');
    }
    await page.keyboard.press('Tab');
    await expect(select).toHaveValue('both');
    await testInfo.attach('after-select', {
      body: JSON.stringify(await page.evaluate(() => ({
        labels: window.__GBDRAW_APP__.form.labels_mode,
        undo: window.__GBDRAW_HISTORY__.getUndoCount(),
        undoLabel: window.__GBDRAW_HISTORY__.undoLabel()
      }))),
      contentType: 'application/json'
    });
    await expectHistory(page, baseline + 1, 0, 'Change setting');
    expect(await page.evaluate(() => window.__GBDRAW_APP__.results[0].content)).toBe(originalSvg);

    await page.getByRole('button', { name: 'Undo', exact: true }).click();
    await expect(select).toHaveValue('out');
    await expectHistory(page, baseline, 1, 'Generate diagram');
    await page.getByRole('button', { name: 'Redo', exact: true }).click();
    await expect(select).toHaveValue('both');
    await expectHistory(page, baseline + 1, 0, 'Change setting');
    const bothSvg = await generate(page);
    expect(bothSvg).not.toBe(originalSvg);
    await expectHistory(page, baseline + 2, 0, 'Generate diagram');

    // Undo the artifact replacement, then the form edit; neither may consume both.
    await page.getByRole('button', { name: 'Undo', exact: true }).click();
    await expectHistory(page, baseline + 1, 1, 'Change setting');
    await expect(select).toHaveValue('both');
    await page.getByRole('button', { name: 'Undo', exact: true }).click();
    await expect(select).toHaveValue('out');
    await expectHistory(page, baseline, 2, 'Generate diagram');

    const downloadPromise = page.waitForEvent('download');
    await page.getByRole('button', { name: 'Save Session', exact: true }).click();
    const download = await downloadPromise;
    const savedPath = testInfo.outputPath('restored.gbdraw-session.json.gz');
    await download.saveAs(savedPath);
    const saved = JSON.parse(gunzipSync(readFileSync(savedPath)));
    expect(saved.config.form.labels_mode).toBe('out');
    expect(await generate(page)).toBe(originalSvg);

    const freshContext = await browser.newContext({ baseURL: 'http://127.0.0.1:4173' });
    try {
      await freshContext.route('**/*', allowLocal);
      const freshPage = await freshContext.newPage();
      await preparePage(freshPage);
      await loadSession(freshPage, savedPath);
      expect(await freshPage.evaluate(() => window.__GBDRAW_APP__.form.labels_mode)).toBe('out');
      expect(await generate(freshPage)).toBe(originalSvg);
    } finally {
      await freshContext.close();
    }
    expect(externalRequests).toEqual([]);
  });
}
