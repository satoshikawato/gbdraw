const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { join, resolve } = require('node:path');
const { openApp } = require('./helpers/app-lifecycle.cjs');

const repoRoot = resolve(process.env.GBDRAW_REPO || process.cwd());
const genbankPath = join(repoRoot, 'tests/test_inputs/HmmtDNA.gbk');
const sessionInputSelector = 'input[accept^=".json,"]';

const runDiagram = async (page) => page.evaluate(async () => {
  const app = window.__GBDRAW_APP__;
  const result = await app.runAnalysis();
  return {
    result,
    errorSummary: String(app.errorLog?.summary || ''),
    errorDetails: Array.isArray(app.errorLog?.details)
      ? app.errorLog.details.map((detail) => String(detail))
      : []
  };
});

const lengthBarFontSizes = (page) => page.evaluate(() => {
  const content = window.__GBDRAW_APP__.results?.[0]?.content || '';
  const svg = new DOMParser().parseFromString(content, 'image/svg+xml').documentElement;
  const lengthBar = svg.getElementById('length_bar');
  return [...(lengthBar?.querySelectorAll('text') || [])]
    .map((text) => Number(text.getAttribute('font-size')))
    .filter(Number.isFinite);
});

const saveSession = async (page, title) => {
  await page.evaluate((nextTitle) => {
    window.__GBDRAW_APP__.sessionTitle = nextTitle;
  }, title);
  const downloadPromise = page.waitForEvent('download', { timeout: 120000 });
  expect(await page.evaluate(() => window.__GBDRAW_APP__.saveSessionWithTitle()))
    .toMatchObject({ status: 'saved' });
  return (await downloadPromise).path();
};

const importSession = async (page, path) => {
  const dialogPromise = page.waitForEvent('dialog', { timeout: 120000 });
  await page.locator(sessionInputSelector).first().setInputFiles(path);
  const dialog = await dialogPromise;
  expect(dialog.message()).toBe('Session loaded successfully!');
  await dialog.accept();
  await page.waitForFunction(() => window.__GBDRAW_APP__?.mode === 'linear');
};

const focusAfterHistoryCapture = async (page, locator) => {
  await locator.focus();
  await page.evaluate(() => new Promise((resolve) => {
    requestAnimationFrame(() => requestAnimationFrame(resolve));
  }));
};

test('independent Linear typography follows linked, imported, and History journeys @pr-smoke', async ({ page }, testInfo) => {
  test.setTimeout(300000);
  const genbank = readFileSync(genbankPath, 'utf8');

  await openApp(page);
  await page.evaluate(async (genbankText) => {
    const app = window.__GBDRAW_APP__;
    app.mode = 'linear';
    app.lInputType = 'gb';
    app.setLinearSeqPrimaryFile(0, 'gb', new File([genbankText], 'HmmtDNA.gbk', {
      type: 'text/plain',
      lastModified: 1
    }));
    Object.assign(app.form, {
      show_scale: true,
      scale_style: 'bar',
      show_gc: false,
      show_skew: false,
      show_depth: false,
      show_labels_linear: 'none'
    });
    app.setLinearTrackSlotsEnabled(false);
    await window.Vue.nextTick();
  }, genbank);

  const axisCard = page.locator('.card').filter({
    has: page.locator('summary').filter({ hasText: 'Axis & Scale' })
  }).first();
  await axisCard.locator('summary').click();
  const link = page.getByLabel('Link Linear scale and ruler label font sizes');
  const scaleSize = page.getByLabel('Linear scale font size');
  const rulerSize = page.getByLabel('Linear ruler label font size');

  await expect(link).toBeChecked();
  await expect(scaleSize).toHaveValue('');
  await expect(rulerSize).toHaveValue('');
  await expect(rulerSize).toBeDisabled();

  await focusAfterHistoryCapture(page, scaleSize);
  await scaleSize.fill('26');
  await scaleSize.press('Tab');
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_HISTORY__.undoLabel()))
    .toBe('Edit setting');
  expect(await page.evaluate(() => ({
    scale: window.__GBDRAW_APP__.adv.scale_font_size,
    ruler: window.__GBDRAW_APP__.adv.ruler_label_font_size
  }))).toEqual({ scale: 26, ruler: 26 });
  await page.evaluate(() => window.__GBDRAW_APP__.undoHistory());
  await expect(scaleSize).toHaveValue('');
  await expect(rulerSize).toHaveValue('');
  await page.evaluate(() => window.__GBDRAW_APP__.redoHistory());
  await expect(scaleSize).toHaveValue('26');
  await expect(rulerSize).toHaveValue('26');

  await focusAfterHistoryCapture(page, link);
  await link.uncheck();
  await expect(rulerSize).toBeEnabled();
  await expect(scaleSize).toHaveValue('26');
  await expect(rulerSize).toHaveValue('26');
  await focusAfterHistoryCapture(page, rulerSize);
  await rulerSize.fill('13');
  await rulerSize.press('Tab');
  await expect(scaleSize).toHaveValue('26');
  await expect(rulerSize).toHaveValue('13');

  await axisCard.locator('summary').click();
  await axisCard.locator('summary').click();
  await expect(link).not.toBeChecked();
  await expect(scaleSize).toHaveValue('26');
  await expect(rulerSize).toHaveValue('13');

  expect(await runDiagram(page)).toEqual({
    result: { status: 'ok' },
    errorSummary: '',
    errorDetails: []
  });
  expect(await lengthBarFontSizes(page)).toEqual([26]);

  await page.getByLabel('Linear scale style').selectOption('ruler');
  expect(await runDiagram(page)).toEqual({
    result: { status: 'ok' },
    errorSummary: '',
    errorDetails: []
  });
  const rulerFonts = await lengthBarFontSizes(page);
  expect(rulerFonts.length).toBeGreaterThan(1);
  expect(new Set(rulerFonts)).toEqual(new Set([13]));
  await page.screenshot({
    path: testInfo.outputPath('linear-typography-ruler.png'),
    fullPage: true
  });

  const unequalSession = await saveSession(page, 'linear-typography-unequal');
  await page.reload({ waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await importSession(page, unequalSession);
  expect(await page.evaluate(() => ({
    linked: window.__GBDRAW_APP__.linearTypographyLinked,
    scale: window.__GBDRAW_APP__.adv.scale_font_size,
    ruler: window.__GBDRAW_APP__.adv.ruler_label_font_size,
    canUndo: window.__GBDRAW_APP__.canUndoHistory
  }))).toEqual({ linked: false, scale: 26, ruler: 13, canUndo: false });

  await page.locator('.card').filter({
    has: page.locator('summary').filter({ hasText: 'Axis & Scale' })
  }).first().locator('summary').click();
  const importedLink = page.getByLabel('Link Linear scale and ruler label font sizes');
  await focusAfterHistoryCapture(page, importedLink);
  await importedLink.check();
  await expect(page.getByText(
    'Relinking copies the current Scale Font Size to Ruler Label Font Size once.'
  )).toBeVisible();
  expect(await page.evaluate(() => ({
    linked: window.__GBDRAW_APP__.linearTypographyLinked,
    scale: window.__GBDRAW_APP__.adv.scale_font_size,
    ruler: window.__GBDRAW_APP__.adv.ruler_label_font_size
  }))).toEqual({ linked: true, scale: 26, ruler: 26 });
  await page.evaluate(() => window.__GBDRAW_APP__.undoHistory());
  expect(await page.evaluate(() => ({
    linked: window.__GBDRAW_APP__.linearTypographyLinked,
    scale: window.__GBDRAW_APP__.adv.scale_font_size,
    ruler: window.__GBDRAW_APP__.adv.ruler_label_font_size
  }))).toEqual({ linked: false, scale: 26, ruler: 13 });
  await page.evaluate(() => window.__GBDRAW_APP__.redoHistory());

  const linkedScaleSize = page.getByLabel('Linear scale font size');
  await focusAfterHistoryCapture(page, linkedScaleSize);
  await linkedScaleSize.fill('30');
  await linkedScaleSize.press('Tab');
  expect(await page.evaluate(() => ({
    linked: window.__GBDRAW_APP__.linearTypographyLinked,
    scale: window.__GBDRAW_APP__.adv.scale_font_size,
    ruler: window.__GBDRAW_APP__.adv.ruler_label_font_size
  }))).toEqual({ linked: true, scale: 30, ruler: 30 });

  const equalLinkedSession = await saveSession(page, 'linear-typography-linked');
  await page.reload({ waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await importSession(page, equalLinkedSession);
  expect(await page.evaluate(() => ({
    linked: window.__GBDRAW_APP__.linearTypographyLinked,
    scale: window.__GBDRAW_APP__.adv.scale_font_size,
    ruler: window.__GBDRAW_APP__.adv.ruler_label_font_size
  }))).toEqual({ linked: true, scale: 30, ruler: 30 });
});
