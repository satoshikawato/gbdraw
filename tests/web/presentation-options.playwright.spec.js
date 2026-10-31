// Two presentation choices, both off by default. Italic gene names
// (`labels.<mode>.italic_gene_labels`, `--italic_gene_labels`) draws the labels
// the label priority names from the gene qualifier in italics; a live label
// text edit shows what the next Generate draws (R3). Transparent background
// keeps the PNG download transparent, as the CLI's PNG is.
const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { inflateSync } = require('node:zlib');
const { join, resolve } = require('node:path');
const { generateAndWaitForResult } = require('./helpers/app-lifecycle.cjs');
const { openWithGenBank } = require('./helpers/audit-browser.cjs');
const { expectLiveEqualsGenerate, settleLive } = require('./helpers/live-generate-parity.cjs');
const { popupEdit } = require('./helpers/live-generate-parity-steps.cjs');

const repoRoot = resolve(process.env.GBDRAW_REPO || process.cwd());
const HMMT = join(repoRoot, 'tests/test_inputs/HmmtDNA.gbk');

// The font style of each drawn feature label, by its text.
const labelStyles = (page) => page.evaluate(() => {
  const root = window.__GBDRAW_APP__.svgContainer.querySelector('svg');
  return Object.fromEntries([...root.querySelectorAll('text[data-gbdraw-label-binding-schema]')].map((text) => [
    text.textContent.trim(), getComputedStyle(text.querySelector('textPath') || text).fontStyle
  ]));
});

const generate = async (page) => {
  await generateAndWaitForResult(page);
  await settleLive(page);
};

test('Italic gene names draws gene-named labels in italics, and a live label edit equals Generate', async ({ page }) => {
  test.setTimeout(300000);
  await openWithGenBank(page, HMMT, () => {
    const app = window.__GBDRAW_APP__;
    app.form.labels_mode = 'out';
    app.autoLabelReflowEnabled = false;
    app.newPriorityRule.feat = 'CDS';
    app.newPriorityRule.order = 'gene,product';
    app.addPriorityRule();
  });
  const italic = page.getByRole('checkbox', { name: 'Italic gene names' });
  await expect(italic).not.toBeChecked();
  await generate(page);
  const upright = await labelStyles(page);
  expect(upright.ND1).toBe('normal');

  await page.evaluate(() => { window.__GBDRAW_APP__.adv.italic_gene_labels = true; });
  await expect(italic).toBeChecked();
  await generate(page);
  const styles = await labelStyles(page);
  expect(styles.ND1).toBe('italic');
  // tRNA labels keep the default priority, product first.
  expect(Object.entries(styles).filter(([text]) => text.startsWith('tRNA-')).map(([, style]) => style))
    .toEqual(expect.arrayContaining(['normal']));
  expect(Object.entries(styles).filter(([text]) => text.startsWith('tRNA-')).every(([, style]) => style === 'normal')).toBe(true);
  const recipe = await page.evaluate(() => window.__GBDRAW_APP__.lastRunInfo?.sourceRecipe?.invocation?.args || []);
  expect(recipe).toContain('--italic_gene_labels');

  await popupEdit(page, { type: 'CDS', gene: 'ND1' }, { labelText: 'ND1 edited live' });
  const { live } = await expectLiveEqualsGenerate(page, { label: 'label text edit of an italic gene label' });
  expect(Object.values(live.labels).filter(({ text }) => text.includes('ND1 edited live')).map(({ fontStyle }) => fontStyle))
    .toEqual([['italic']]);
});

// The alpha of the PNG's top left pixel: with no row above and no pixel to its
// left, every PNG filter stores that pixel unchanged.
const cornerAlpha = (png) => {
  expect(png.subarray(12, 16).toString('latin1')).toBe('IHDR');
  expect([png[24], png[25]]).toEqual([8, 6]);
  const data = [];
  for (let at = 8; at < png.length;) {
    const length = png.readUInt32BE(at);
    if (png.subarray(at + 4, at + 8).toString('latin1') === 'IDAT') data.push(png.subarray(at + 8, at + 8 + length));
    at += 12 + length;
  }
  return inflateSync(Buffer.concat(data))[4];
};

test('Transparent background keeps the PNG download transparent; off paints it white', async ({ page }) => {
  test.setTimeout(300000);
  await openWithGenBank(page, HMMT, () => { window.__GBDRAW_APP__.form.labels_mode = 'none'; });
  await generate(page);
  const downloadPng = async () => {
    const pending = page.waitForEvent('download', { timeout: 120000 });
    await page.getByRole('button', { name: 'PNG', exact: true }).click();
    return readFileSync(await (await pending).path());
  };
  const transparent = page.getByRole('checkbox', { name: 'Transparent background' });
  await expect(transparent).not.toBeChecked();
  expect(cornerAlpha(await downloadPng())).toBe(255);
  await transparent.check();
  expect(cornerAlpha(await downloadPng())).toBe(0);
});
