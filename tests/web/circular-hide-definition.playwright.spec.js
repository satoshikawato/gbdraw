// PD-OI-109: Circular Show Definition. Off draws no record's center
// definition, keeps the plot title, hides Keep Full Definition with Plot Title,
// and is saved with the Session and written to Run info as --hide_definition.
const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { gunzipSync } = require('node:zlib');
const { join, resolve } = require('node:path');
const { shownSavedDrawing } = require('./helpers/saved-session.cjs');
const {
  CURRENT_SESSION_VERSION,
  answerSessionTitle,
  expectSessionLoadedWithoutWarnings,
  generateAndWaitForResult,
  openApp,
  waitForAppShell
} = require('./helpers/app-lifecycle.cjs');

const repoRoot = resolve(process.env.GBDRAW_REPO || process.cwd());
const genbankPath = join(repoRoot, 'tests/test_inputs/HmmtDNA.gbk');
const SHOW_PATH = 'objects.definition.circular.show';

const inspectResult = async (page) => page.evaluate(() => {
  const content = window.__GBDRAW_APP__.currentResults()[0]?.content || '';
  const svg = new DOMParser().parseFromString(content, 'image/svg+xml').documentElement;
  return {
    definitions: svg.querySelectorAll('g[data-gbdraw-role="record-definition"]').length,
    titles: [...svg.querySelectorAll('g[data-gbdraw-role="plot-title"]')].map((group) => group.textContent.trim())
  };
});

const readSession = (path) => {
  const raw = readFileSync(path);
  return JSON.parse((raw[0] === 0x1f && raw[1] === 0x8b ? gunzipSync(raw) : raw).toString('utf8'));
};

test('Show Definition off draws no center definition, keeps the plot title, and survives Save and Load', async ({ page }, testInfo) => {
  test.setTimeout(300000);
  await openApp(page);
  // A new Session shows its one Circular drawing; the file is registered as a chosen file is.
  await page.evaluate(async (genbankText) => {
    const app = window.__GBDRAW_APP__;
    const { registerInputFileTree } = await import('./js/services/drawing/inputs.js');
    const { state } = await import('./js/state.js');
    app.cInputType = 'gb';
    app.files.c_gb = registerInputFileTree(state.inputFiles,
      new File([genbankText], 'HmmtDNA.gbk', { type: 'text/plain', lastModified: 1 }));
    app.form.labels_mode = 'none';
    app.form.plot_title = 'Mitochondrial genome';
    app.adv.plot_title_position = 'top';
    await window.Vue.nextTick();
    document.querySelector('[data-setting="titles"] details').open = true;
  }, readFileSync(genbankPath, 'utf8'));

  const showDefinition = page.getByRole('checkbox', { name: 'Show Definition' });
  const keepFull = page.getByText('Keep Full Definition with Plot Title');
  await expect(showDefinition).toBeChecked();
  await expect(keepFull).toBeVisible();
  await generateAndWaitForResult(page);
  expect(await inspectResult(page)).toEqual({ definitions: 1, titles: ['Mitochondrial genome'] });

  await showDefinition.uncheck();
  await expect(keepFull).toBeHidden();
  await generateAndWaitForResult(page);
  expect(await inspectResult(page)).toEqual({ definitions: 0, titles: ['Mitochondrial genome'] });
  const recipe = await page.evaluate(() => window.__GBDRAW_APP__.lastRunInfo?.sourceRecipe?.invocation?.args || []);
  expect(recipe).toContain('--hide_definition');

  const pending = page.waitForEvent('download', { timeout: 180000 });
  await page.getByRole('banner').getByRole('button', { name: 'Save Session', exact: true }).click();
  await answerSessionTitle(page, '');
  const download = await pending;
  const savedPath = testInfo.outputPath(download.suggestedFilename());
  await download.saveAs(savedPath);
  const saved = readSession(savedPath);
  expect(saved.version).toBe(CURRENT_SESSION_VERSION);
  // The shown drawing holds the setting and its request (Session 50: one entry per drawing).
  expect(shownSavedDrawing(saved).config.adv.show_definition).toBe(false);
  expect(shownSavedDrawing(saved).renderRequest.diagramOptions.configOverrides[SHOW_PATH]).toBe(false);

  await page.reload({ waitUntil: 'domcontentloaded' });
  await waitForAppShell(page, { waitForPalette: false });
  await page.locator('input[accept^=".json,"]').first().setInputFiles(savedPath);
  await expectSessionLoadedWithoutWarnings(page, { timeout: 180000 });
  await page.evaluate(() => { document.querySelector('[data-setting="titles"] details').open = true; });
  await expect(showDefinition).not.toBeChecked();
  await expect(keepFull).toBeHidden();
  await generateAndWaitForResult(page);
  expect(await inspectResult(page)).toEqual({ definitions: 0, titles: ['Mitochondrial genome'] });
});
