// @ts-check
// The steps of the drawing menu's "New drawing from a region…" dialog
// (app/region-drawing.js) that several specs take.
const { expect } = require('@playwright/test');

/** @param {import('@playwright/test').Page} page */
const openRegionDialog = async (page) => {
  const menu = page.locator('.drawing-bar-menu');
  await menu.locator('summary').click();
  await menu.getByRole('button', { name: 'New drawing from a region…', exact: true }).click();
  const dialog = page.locator('[data-region-drawing-dialog]');
  // The rows wait for the shown drawing's request, which Python projects (the first call boots it).
  await expect(dialog.locator('[data-region-drawing-row]')).not.toHaveCount(0, { timeout: 180_000 });
  return dialog;
};

/**
 * Check the first record and enter Coordinates `start`..`end`.
 * @param {import('@playwright/test').Locator} dialog @param {number} start @param {number} end
 */
const chooseCoordinates = async (dialog, start, end) => {
  const row = dialog.locator('[data-region-drawing-row]').first();
  await row.locator('[data-region-drawing-check]').check();
  await row.locator('[data-region-drawing-kind]').selectOption('coordinates');
  await row.locator('[data-region-drawing-start]').fill(String(start));
  await row.locator('[data-region-drawing-end]').fill(String(end));
};

module.exports = { chooseCoordinates, openRegionDialog };
