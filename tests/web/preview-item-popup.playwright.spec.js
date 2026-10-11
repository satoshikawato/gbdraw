// A click (or Enter) on a Legend row or a record definition in the preview
// opens its popup while Layout edit is off. The popup edits through the Legend
// panel's operations and the sidebar's record settings, one History step each;
// Escape returns focus to the diagram item. With Layout edit on, no popup opens.
const { test, expect } = require('@playwright/test');
const { HMMT_SESSION, loadSessionFile, openFresh, settle } = require('./helpers/audit-browser.cjs');

test.describe.configure({ retries: 0 });

const LINEAR_SESSION = 'gbdraw/web/gallery/sessions/lambda_basic_linear.gbdraw-session.json';
const undoCount = (page) => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
const popup = (page) => page.locator('.preview-item-popup[role="dialog"]');
const legendRow = (page, key) => page.locator(`.gbdraw-preview-surface g[data-legend-key="${key}"]`).first();

test('a Legend row opens its popup; rename, remove, keyboard, and Layout edit', async ({ page }) => {
  test.setTimeout(300_000);
  await openFresh(page);
  await loadSessionFile(page, HMMT_SESSION);

  await legendRow(page, 'tRNA').locator('path').first().click();
  const dialog = page.getByRole('dialog', { name: 'Legend entry: tRNA', exact: true });
  await expect(dialog).toBeVisible();
  await expect(dialog.locator('[data-preview-legend-color-source]')).toHaveText('Default color of tRNA features');

  // Rename: one History step, shown in the Legend panel's rows.
  const before = await undoCount(page);
  await dialog.getByLabel('Caption', { exact: true }).fill('Transfer RNA');
  await dialog.getByLabel('Caption', { exact: true }).press('Enter');
  await settle(page);
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.legendEntries.map((entry) => entry.caption)))
    .toContain('Transfer RNA');
  expect(await undoCount(page)).toBe(before + 1);
  await expect(page.getByRole('dialog', { name: 'Legend entry: Transfer RNA', exact: true })).toBeVisible();

  // Remove from legend: the row moves to the deleted list and the popup closes.
  await popup(page).getByRole('button', { name: 'Remove from legend', exact: true }).click();
  await settle(page);
  await expect(popup(page)).toHaveCount(0);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.deletedLegendEntries.map((entry) => entry.caption)))
    .toContain('Transfer RNA');

  // Keyboard: Enter opens the popup with focus inside; Escape returns focus.
  const key = await page.evaluate(() => window.__GBDRAW_APP__.legendEntries[0].caption);
  const row = legendRow(page, key);
  await row.focus();
  await page.keyboard.press('Enter');
  await expect(popup(page)).toBeVisible();
  expect(await popup(page).evaluate((element) => element.contains(document.activeElement))).toBe(true);
  await page.keyboard.press('Escape');
  await expect(popup(page)).toHaveCount(0);
  expect(await row.evaluate((element) => element === document.activeElement)).toBe(true);

  // Layout edit on: a press drags the legend; no popup opens.
  await page.evaluate(() => { window.__GBDRAW_APP__.layoutRepositionMode = true; });
  await legendRow(page, key).locator('path').first().click();
  await expect(popup(page)).toHaveCount(0);

  // The exported SVG carries none of the preview's popup affordances.
  const exported = await page.evaluate(async () => {
    const { serializeCleanSvg } = await import('./js/services/svg-serialization.js');
    return serializeCleanSvg(window.__GBDRAW_APP__.svgContainer.querySelector('svg'));
  });
  expect(exported).not.toContain('data-gbdraw-preview-item');
});

test('a record definition opens its popup; its edits write the record settings', async ({ page }) => {
  test.setTimeout(300_000);
  await openFresh(page);
  await loadSessionFile(page, LINEAR_SESSION);

  await page.locator('.gbdraw-preview-surface [data-gbdraw-role^="record-definition"] text').first().click();
  const dialog = popup(page);
  await expect(dialog).toBeVisible();
  await expect(dialog).toHaveAttribute('aria-label', /^Record /);
  await expect(dialog.locator('dl')).toContainText('Accession');
  await expect(dialog.locator('dl')).toContainText('Whole record');

  const before = await undoCount(page);
  await dialog.getByRole('checkbox', { name: 'Reverse complement', exact: true }).check();
  await settle(page);
  expect(await undoCount(page)).toBe(before + 1);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.linearSeqs[0].region_reverse)).toBe(true);
  await dialog.getByLabel('Subtitle', { exact: true }).fill('Popup subtitle');
  await dialog.getByLabel('Subtitle', { exact: true }).press('Tab');
  await settle(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.linearSeqs[0].record_subtitle)).toBe('Popup subtitle');
  expect(await undoCount(page)).toBe(before + 2);
  await page.evaluate(() => window.__GBDRAW_HISTORY__.undo());
  await settle(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.linearSeqs[0].record_subtitle || '')).not.toBe('Popup subtitle');
});
