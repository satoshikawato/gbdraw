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
});

// OV-373: a single-record Circular Session shows its annotation records as a
// plain record name before Inspect, never "target record unavailable".
test('OV-373: annotation rows name their record before Inspect', async ({ page }) => {
  test.setTimeout(240000);
  await openApp(page);
  await loadSession(page, 'tobacco-chloroplast.gbdraw-session.json');
  const panel = page.locator('details').filter({ has: page.locator('summary[aria-label="Region Annotations"]') });
  await panel.locator('summary').first().click();
  await expect(panel).not.toContainText('unavailable');
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
