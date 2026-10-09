// D-15 (Owner 2026-10-09): the palette is the base layer and the user default
// colors (Default colors values that differ from the selected palette's color)
// win over it, as `-d` over `-p`. A palette switch or Default colors Reset asks
// before it drops them; a dialog choice is one History step and Cancel records
// none (PD-OI-088). Node cases: tests/web/palette-default-colors.test.mjs.
const { test, expect } = require('@playwright/test');
const { evaluateWithRetainedPromise } = require('./helpers/app-lifecycle.cjs');
const { expectLiveEqualsGenerate, settleLive } = require('./helpers/live-generate-parity.cjs');
const { open } = require('./helpers/live-generate-parity-steps.cjs');

test.describe.configure({ retries: 0 });

const USER_CDS = '#123456';

const openWithUserColor = async (page) => {
  await open(page, { mode: 'circular', results: 'single', reflow: 'off' });
  await evaluateWithRetainedPromise(page, async (color) => {
    const { state } = await import('/gbdraw/web/js/state.js');
    await window.__GBDRAW_HISTORY__.runUndoable('Change color', async () => {
      state.activeDrawing().currentColors.value = { ...state.activeDrawing().currentColors.value, CDS: color };
      await window.Vue.nextTick();
    });
  }, USER_CDS);
  await settleLive(page);
  const colors = page.locator('summary[aria-label="Colors"]');
  if ((await colors.locator('..').getAttribute('open')) === null) await colors.click();
};
const facts = (page) => page.evaluate(() => {
  const app = window.__GBDRAW_APP__;
  return {
    palette: app.selectedPalette,
    colors: { ...app.currentColors },
    undo: window.__GBDRAW_HISTORY__.getUndoCount()
  };
});

test('a palette switch keeps a user default color after asking; Cancel changes nothing', async ({ page }) => {
  test.setTimeout(240_000);
  await openWithUserColor(page);
  // Instant Preview applies the palette live, so the live Result can equal Generate.
  await page.evaluate(async () => {
    (await import('/gbdraw/web/js/state.js')).state.paletteInstantPreviewEnabled.value = true;
  });
  const select = page.getByRole('combobox', { name: 'Palette', exact: true });
  const before = await facts(page);
  const next = await page.evaluate(() => window.__GBDRAW_APP__.paletteNames.find((name) => (
    name !== window.__GBDRAW_APP__.selectedPalette
  )));
  const dialog = page.getByRole('dialog', { name: 'Change palette' });

  await select.selectOption(next);
  await expect(dialog).toBeVisible();
  await expect(select).toHaveValue(before.palette);
  await dialog.getByRole('button', { name: 'Cancel', exact: true }).click();
  await expect(dialog).toHaveCount(0);
  await expect(select).toHaveValue(before.palette);
  expect(await facts(page)).toEqual(before);

  await select.selectOption(next);
  await dialog.getByRole('button', { name: 'Keep my 1 color', exact: true }).click();
  await expect(dialog).toHaveCount(0);
  await expect(select).toHaveValue(next);
  const kept = await facts(page);
  expect(kept.palette).toBe(next);
  expect(kept.colors.CDS).toBe(USER_CDS);
  expect(kept.undo).toBe(before.undo + 1);
  const paletteTrna = await page.evaluate((name) => window.__GBDRAW_APP__.paletteDefinitions[name].tRNA, next);
  expect(kept.colors.tRNA).toBe(paletteTrna);

  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
  await settleLive(page);
  await expect(select).toHaveValue(before.palette);
  expect((await facts(page)).colors).toEqual(before.colors);
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.redo());
  await expect(select).toHaveValue(next);
  // Generate records its own History step, so the parity check comes last.
  await expectLiveEqualsGenerate(page, { label: 'Keep my 1 color' });
});

test('Default colors Reset asks before it discards a user default color', async ({ page }) => {
  test.setTimeout(180_000);
  await openWithUserColor(page);
  const reset = page.locator('h4[aria-label="DEFAULT COLORS"]').locator('..')
    .getByRole('button', { name: 'Reset', exact: true });
  const dialog = page.getByRole('dialog', { name: 'Reset default colors' });
  const before = await facts(page);

  await reset.click();
  await expect(dialog).toBeVisible();
  await dialog.getByRole('button', { name: 'Cancel', exact: true }).click();
  await expect(dialog).toHaveCount(0);
  expect(await facts(page)).toEqual(before);

  await reset.click();
  await dialog.getByRole('button', { name: "Reset to the palette's colors", exact: true }).click();
  await expect(dialog).toHaveCount(0);
  const paletteCds = await page.evaluate(() => (
    window.__GBDRAW_APP__.paletteDefinitions[window.__GBDRAW_APP__.selectedPalette].CDS
  ));
  const after = await facts(page);
  expect(after.colors.CDS).toBe(paletteCds);
  expect(after.undo).toBe(before.undo + 1);

  // Without user default colors, Reset applies at once.
  await reset.click();
  await settleLive(page);
  await expect(dialog).toHaveCount(0);
  expect((await facts(page)).colors.CDS).toBe(paletteCds);
});
