// Sidebar input cards on loaded Gallery Sessions: record readouts before and
// after the source records are inspected.
const { test, expect } = require('@playwright/test');
const { join } = require('node:path');
const { openApp } = require('./helpers/app-lifecycle.cjs');

const loadSession = async (page, name) => {
  page.on('dialog', (dialog) => dialog.accept());
  await page.locator('input[accept^=".json,"]').setInputFiles(
    join(process.cwd(), 'gbdraw/web/gallery/sessions', name));
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending
    && window.__GBDRAW_APP__.results.length > 0, null, { timeout: 180000 });
};

// OV-374: after Inspect, the Display start field shows the start the record
// is drawn from, the same value the readout used to show beside an empty field.
test('OV-374: Display start shows the current start after Inspect', async ({ page }) => {
  test.setTimeout(240000);
  await openApp(page);
  await loadSession(page, 'HmmtDNA_ATskew.gbdraw-session.json');
  await page.locator('[data-circular-inspect]').click();
  const start = page.getByLabel(/^Display start NC_012920\.1/);
  await expect(start).toBeVisible({ timeout: 180000 });
  await expect(start).toHaveValue('1');

  // Reset start leaves with its change; the keyboard focus moves to the field.
  await start.fill('100');
  await start.press('Tab');
  const reset = page.getByRole('button', { name: 'Reset start', exact: true });
  await reset.focus();
  await page.keyboard.press('Enter');
  await expect(start).toBeFocused();
  await expect(start).toHaveValue('1');

  // Clearing the field requests no start, also when none was requested, so
  // it shows the current start again.
  await start.fill('');
  await start.press('Tab');
  await expect(start).toHaveValue('1');

  // With Circular record off the start field is disabled; Reset start then
  // moves the keyboard focus to that checkbox instead of dropping it.
  await start.fill('100');
  await start.press('Tab');
  const topology = page.getByLabel(/^Circular record NC_012920\.1/);
  await topology.uncheck();
  await expect(start).toBeDisabled();
  await reset.focus();
  await page.keyboard.press('Enter');
  await expect(reset).toHaveCount(0);
  await expect(topology).toBeFocused();
});

// OV-373: a single-record Circular Session shows its annotation records as a
// plain record name before Inspect, never "target record unavailable".
test('OV-373: annotation rows name their record before Inspect', async ({ page }) => {
  test.setTimeout(240000);
  await openApp(page);
  await loadSession(page, 'tobacco-chloroplast.gbdraw-session.json');
  const panel = page.locator('details').filter({ has: page.locator('summary[aria-label="Region Annotations"]') });
  await panel.locator('summary').first().click();
  await expect(panel.locator('[data-annotation-record-caption]').first()).toHaveText('Record: NC_001879.2');
  await expect(panel).not.toContainText('unavailable');
});

// A default color edited under "More colors" stays in that grid and keeps
// focus while it is edited; the grouping does not follow each edit.
test('a color edited under More colors stays there and keeps focus', async ({ page }) => {
  test.setTimeout(240000);
  await openApp(page);
  await loadSession(page, 'HmmtDNA_ATskew.gbdraw-session.json');
  await page.locator('summary[aria-label="Colors"]').click();
  const more = page.locator('[data-default-colors-more]');
  await more.locator('summary').click();
  const tmrna = more.getByLabel('tmRNA feature color mode');
  await expect(tmrna).toBeVisible();
  await expect(tmrna).toHaveValue('color');
  await tmrna.focus();
  await page.keyboard.press('ArrowUp');
  await expect(more.getByLabel('tmRNA feature color mode')).toHaveValue('none');
  await expect(more.getByLabel('tmRNA feature color mode')).toBeFocused();
});

// A Default colors file replaces the color table, so a type it colors joins
// the user's colors in the first list, as a color added by hand does.
test('a color from a Default colors file leaves More colors', async ({ page }) => {
  test.setTimeout(240000);
  await openApp(page);
  await loadSession(page, 'HmmtDNA_ATskew.gbdraw-session.json');
  await page.locator('summary[aria-label="Colors"]').click();
  const more = page.locator('[data-default-colors-more]');
  await more.locator('summary').click();
  await expect(more.getByLabel('tmRNA feature color mode')).toBeVisible();
  await page.getByLabel('Default colors file', { exact: true }).setInputFiles({
    name: 'colors.tsv', mimeType: 'text/plain', buffer: Buffer.from('tmRNA\t#ff0000\n')
  });
  await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    await app.waitForAuxiliaryFileImport(app.files.d_color);
  });
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.currentColors.tmRNA)).toBe('#ff0000');
  await expect(more.getByLabel('tmRNA feature color mode')).toHaveCount(0);
  await expect(page.getByLabel('tmRNA feature color mode')).toBeVisible();
});

// OV-379: a loaded Linear Session has not read its records; the record select
// says so instead of a load that never runs.
test('OV-379: Linear record selects say Not inspected after Session load', async ({ page }) => {
  test.setTimeout(240000);
  await openApp(page);
  await page.getByRole('button', { name: 'Linear', exact: true }).click();
  await loadSession(page, 'BGC0000708-BGC0000713.gbdraw-session.json');
  const labels = await page.locator('select[data-record-selector-uid]').evaluateAll((selects) => (
    selects.map((select) => select.options[select.selectedIndex]?.text || '')));
  expect(labels.length).toBeGreaterThan(0);
  expect(new Set(labels)).toEqual(new Set(['Not inspected']));
});
