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
  test(`Label Mode ${inputMethod} edit has one Undo step and survives generation and fresh Load`, async ({
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

test('Linear File removal choices are atomic, undoable, and preserve one slot', async ({ page }) => {
  test.setTimeout(120_000);
  await page.setViewportSize({ width: 390, height: 844 });
  await preparePage(page);
  await page.getByRole('button', { name: 'Linear', exact: true }).click();
  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const first = new File(['LOCUS A\n'], 'first.gbk', { type: 'text/plain', lastModified: 1 });
    const last = new File(['LOCUS B\n'], 'last.gbk', { type: 'text/plain', lastModified: 2 });
    app.linearSeqs[0].gb = first;
    app.linearSeqs[0].region_record_id = 'A1';
    app.linearSeqs[0].definition = 'First one';
    app.linearSeqs[0].depth = [new File(['1\t10\n'], 'first.tsv')];
    app.addLinearSeq();
    app.linearSeqs[1].gb = first;
    app.linearSeqs[1].region_record_id = 'A2';
    app.linearSeqs[1].definition = 'First two';
    app.addLinearSeq();
    app.linearSeqs[2].gb = last;
    app.linearSeqs[2].region_record_id = 'B1';
  });

  const sources = page.locator('[data-linear-source-card]');
  const dialog = page.getByRole('dialog', { name: 'Clear or delete File?' });
  await expect(sources).toHaveCount(2);
  const baseline = await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
  const firstRemove = sources.first().getByRole('button', { name: /Remove$/ });

  await firstRemove.click();
  expect(await page.evaluate(() => ({
    dialog: { ...window.__GBDRAW_APP__.linearSourceRemovalDialog },
    target: window.__GBDRAW_APP__.linearSourceRemovalTarget?.uid || null
  }))).toEqual({
    dialog: { open: true, sourceUid: expect.any(String), origin: 'card' },
    target: expect.any(String)
  });
  await expect(dialog).toBeVisible();
  await expect(dialog).toContainText('first.gbk');
  await expect(dialog).toContainText('2 records');
  const clearFile = dialog.getByRole('button', { name: 'Clear file only', exact: true });
  const cancel = dialog.getByRole('button', { name: 'Cancel', exact: true });
  await expect(clearFile).toBeFocused();
  await page.keyboard.press('Shift+Tab');
  await expect(cancel).toBeFocused();
  await page.keyboard.press('Tab');
  await expect(clearFile).toBeFocused();
  await cancel.click();
  await expect(dialog).toHaveCount(0);
  await expect(firstRemove).toBeFocused();
  expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(baseline);

  await firstRemove.focus();
  await firstRemove.press('Enter');
  await expect(dialog).toBeVisible();
  await page.keyboard.press('Escape');
  await expect(dialog).toHaveCount(0);
  expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(baseline);

  await firstRemove.click();
  await dialog.getByRole('button', { name: 'Clear file only', exact: true }).click();
  await expect(sources).toHaveCount(2);
  await expect(sources.first().getByRole('button', { name: 'Choose GenBank File', exact: true })).toBeFocused();
  expect(await page.evaluate(() => window.__GBDRAW_APP__.linearSeqs.map((sequence) => ({
    uid: sequence.uid,
    file: sequence.gb?.name || null,
    selector: sequence.region_record_id,
    definition: sequence.definition,
    depth: (sequence.depth || []).map((file) => file?.name || null)
  })))).toEqual([
    expect.objectContaining({ file: null, selector: '', definition: '', depth: [null] }),
    expect.objectContaining({ file: 'last.gbk', selector: 'B1' })
  ]);
  await expectHistory(page, baseline + 1, 0, 'Clear Linear File');
  await page.getByRole('button', { name: 'Undo', exact: true }).click();
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.linearSeqs.map(
    (sequence) => sequence.region_record_id
  ))).toEqual(['A1', 'A2', 'B1']);
  await expect.poll(() => page.evaluate(() => [
    window.__GBDRAW_HISTORY__.getUndoCount(),
    window.__GBDRAW_HISTORY__.getRedoCount()
  ])).toEqual([baseline, 1]);
  await page.getByRole('button', { name: 'Redo', exact: true }).click();
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.linearSeqs.map(
    (sequence) => sequence.region_record_id
  ))).toEqual(['', 'B1']);
  await expectHistory(page, baseline + 1, 0, 'Clear Linear File');
  await page.getByRole('button', { name: 'Undo', exact: true }).click();

  const globalRemove = page.getByRole('button', { name: 'Remove last sequence', exact: true });
  await globalRemove.click();
  const globalDialog = page.getByRole('dialog', { name: 'Remove last File?' });
  await expect(globalDialog).toBeVisible();
  await expect(globalDialog).toContainText('last.gbk');
  await globalDialog.getByRole('button', { name: 'Cancel', exact: true }).click();
  expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(baseline);
  await globalRemove.click();
  await globalDialog.getByRole('button', { name: 'Delete card', exact: true }).click();
  await expect(sources).toHaveCount(1);
  await expect(sources.first().getByRole('button', { name: 'Choose GenBank File', exact: true })).toBeFocused();
  await expectHistory(page, baseline + 1, 0, 'Delete Linear File');
  await page.getByRole('button', { name: 'Undo', exact: true }).click();
  await expect(sources).toHaveCount(2);

  await page.getByRole('group', { name: 'Linear File actions' })
    .getByRole('button', { name: 'Add sequence', exact: true }).click();
  await expect(sources).toHaveCount(3);
  const beforeBlankRemoval = await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
  await globalRemove.click();
  await expect(globalDialog).toHaveCount(0);
  await expect(sources).toHaveCount(2);
  await expectHistory(page, beforeBlankRemoval + 1, 0, 'Delete Linear File');
  expect(await page.locator('[data-linear-record-list]').evaluate(
    (element) => element.scrollWidth <= element.clientWidth + 1
  )).toBe(true);
});
