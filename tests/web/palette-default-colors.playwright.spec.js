// PD-OI-089, OIC-029 (Owner 2026-10-09): the palette is the base layer and the user default
// colors (Default colors values that differ from the selected palette's color)
// win over it, as `-d` over `-p`. A palette switch or Default colors Reset asks
// before it drops them; a dialog choice is one History step and Cancel records
// none (PD-OI-088). Node cases: tests/web/palette-default-colors.test.mjs.
const { test, expect } = require('@playwright/test');
const { evaluateWithRetainedPromise } = require('./helpers/app-lifecycle.cjs');
const { expectLiveEqualsGenerate, settleLive } = require('./helpers/live-generate-parity.cjs');
const { generate, open, renameRow } = require('./helpers/live-generate-parity-steps.cjs');

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
// OV-276: the Legend panel row shows the color of the repainted swatch.
const legendRowColor = (page, caption) => page.evaluate((row) => {
  const app = window.__GBDRAW_APP__;
  return String(app.legendEntryColor(app.legendEntries.find((entry) => entry.caption === row)) || '').toLowerCase();
}, caption);
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
  // OV-276: a Default colors edit shows in the Legend panel row.
  await expect.poll(() => legendRowColor(page, 'CDS')).toBe(USER_CDS);
  const next = await page.evaluate(() => window.__GBDRAW_APP__.paletteNames.find((name) => (
    name !== window.__GBDRAW_APP__.selectedPalette
  )));
  // D-24: the dialog names both palettes.
  const dialog = page.getByRole('dialog', { name: `Change palette to "${next}"`, exact: true });

  await select.selectOption(next);
  await expect(dialog).toBeVisible();
  await expect(dialog).toContainText(`from the "${before.palette}" palette`);
  await expect(dialog).not.toContainText('is dropped');
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
  // OV-276: so does a palette switch, on every row of a palette key.
  await expect.poll(() => page.evaluate((name) => {
    const app = window.__GBDRAW_APP__;
    const palette = app.paletteDefinitions[name];
    const rows = app.legendEntries.filter((entry) => entry.caption !== 'CDS' && palette[entry.caption]);
    return rows.length > 0 && rows.every((entry) => (
      String(app.legendEntryColor(entry)).toLowerCase() === palette[entry.caption].toLowerCase()
    ));
  }, next)).toBe(true);
  await expect.poll(() => legendRowColor(page, 'CDS')).toBe(USER_CDS);

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

// OV-281: a Default colors key set to Auto shows the selected palette's color,
// as Generate draws it, and going back to Color starts from that color, so the
// key does not become a user default color.
test('a Default colors key set to Auto shows the palette color, and Color starts from it', async ({ page }) => {
  test.setTimeout(180_000);
  await open(page, { mode: 'circular', results: 'single', reflow: 'off' });
  const colors = page.locator('summary[aria-label="Colors"]');
  if ((await colors.locator('..').getAttribute('open')) === null) await colors.click();
  const panel = page.locator('h4[aria-label="DEFAULT COLORS"]').locator('../..');
  const paletteCds = await page.evaluate(() => String(
    window.__GBDRAW_APP__.paletteDefinitions[window.__GBDRAW_APP__.selectedPalette].CDS
  ).toLowerCase());
  const mode = panel.getByRole('combobox', { name: 'CDS feature color mode', exact: true });

  await mode.selectOption('auto');
  await settleLive(page);
  await expect(panel.locator('input[type="color"][aria-label="CDS feature color"]')).toHaveValue(paletteCds);
  await mode.selectOption('color');
  await settleLive(page);
  expect(String((await facts(page)).colors.CDS).toLowerCase()).toBe(paletteCds);

  // Not a user default color: a palette switch asks nothing.
  const next = await page.evaluate(() => window.__GBDRAW_APP__.paletteNames.find((name) => (
    name !== window.__GBDRAW_APP__.selectedPalette
  )));
  await page.getByRole('combobox', { name: 'Palette', exact: true }).selectOption(next);
  await expect.poll(async () => (await facts(page)).palette).toBe(next);
  await expect(page.getByRole('dialog', { name: `Change palette to "${next}"`, exact: true })).toHaveCount(0);
});

// Q1 B and Q5 A (Owner 2026-10-09). While a palette is queued (Instant Preview
// off), Apply to all on a palette row also shows its color on the Result now;
// switching back to the applied palette asks like any switch and then applies
// live.
test('while a palette is queued, Apply to all shows its color now and switching back asks', async ({ page }) => {
  test.setTimeout(240_000);
  await open(page, { mode: 'circular', results: 'single', reflow: 'off' });
  await page.evaluate(async () => {
    (await import('/gbdraw/web/js/state.js')).state.paletteInstantPreviewEnabled.value = false;
  });
  const colors = page.locator('summary[aria-label="Colors"]');
  if ((await colors.locator('..').getAttribute('open')) === null) await colors.click();
  const select = page.getByRole('combobox', { name: 'Palette', exact: true });
  const applied = await page.evaluate(() => window.__GBDRAW_APP__.selectedPalette);
  const next = await page.evaluate(() => window.__GBDRAW_APP__.paletteNames.find((name) => (
    name !== window.__GBDRAW_APP__.selectedPalette
  )));
  await select.selectOption(next);
  await expect(select).toHaveValue(next);
  const queued = () => page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    return {
      pending: app.pendingPaletteName,
      applied: String(app.appliedPaletteColors.CDS || '').toLowerCase(),
      queuedCds: String(app.pendingPaletteColors.CDS || '').toLowerCase(),
      cds: String(app.currentColors.CDS || '').toLowerCase(),
      rules: app.manualSpecificRules.length
    };
  });
  expect((await queued()).pending).toBe(next);

  // D-24: the scope dialog's default-color line says the color applies now.
  await evaluateWithRetainedPromise(page, async (color) => {
    const app = window.__GBDRAW_APP__;
    await app.openFeatureEditorFromList(app.filteredFeatures.find((item) => item.locus_tag === 'FL1'), null);
    await window.Vue.nextTick();
    await app.updateClickedFeatureColor(color);
  }, USER_CDS);
  const line = page.locator('[data-default-color-scope-line]');
  await expect(line).toContainText(
    'Sets the CDS default color (CDS features without their own color or rule, also hidden ones).'
  );
  await expect(line).toContainText('Applies now, also over the queued palette.');
  await evaluateWithRetainedPromise(page, async () => {
    const app = window.__GBDRAW_APP__;
    await app.handleFeatureStyleScopeChoice('caption');
    app.clickedFeature = null;
  });
  await settleLive(page);
  expect(await queued()).toEqual({ pending: next, applied: USER_CDS, queuedCds: USER_CDS, cds: USER_CDS, rules: 0 });
  await expect.poll(() => legendRowColor(page, 'CDS')).toBe(USER_CDS);

  const dialog = page.getByRole('dialog', { name: `Change palette to "${applied}"`, exact: true });
  await select.selectOption(applied);
  await expect(dialog).toBeVisible();
  // D-24: switching back drops the queued palette, and the colors apply now.
  await expect(dialog).toContainText(`The queued "${next}" palette is dropped and the colors apply now.`);
  await dialog.getByRole('button', { name: 'Keep my 1 color', exact: true }).click();
  await expect(dialog).toHaveCount(0);
  await expect(select).toHaveValue(applied);
  expect(await queued()).toMatchObject({ pending: '', applied: USER_CDS, cds: USER_CDS });
  await expectLiveEqualsGenerate(page, { label: 'switch back to the applied palette, keeping the CDS color' });
});

// OV-282 (D-26): a palette row's color is read from one function, never from
// the color the row recorded at the last Generate. After a Default colors
// edit, a Default colors Reset and a palette switch (Instant Preview), "Use
// existing color" offers the color the CDS row shows, and a Legend editor
// rename of the row keeps its features in that color (OV-311).
test('a palette row offers and keeps the color it shows after palette edits', async ({ page }) => {
  test.setTimeout(300_000);
  await open(page, { mode: 'circular', results: 'single', reflow: 'off' });
  await page.evaluate(async () => {
    (await import('/gbdraw/web/js/state.js')).state.paletteInstantPreviewEnabled.value = true;
  });
  const colors = page.locator('summary[aria-label="Colors"]');
  if ((await colors.locator('..').getAttribute('open')) === null) await colors.click();
  // The "Use existing" color the scope dialog offers for a CDS feature; Cancel changes nothing.
  const offered = () => evaluateWithRetainedPromise(page, async () => {
    const app = window.__GBDRAW_APP__;
    await app.openFeatureEditorFromList(app.filteredFeatures.find((item) => item.locus_tag === 'FL1'), null);
    await window.Vue.nextTick();
    await app.updateClickedFeatureColor('#0f0f0f');
    const color = app.featureStyleScopeDialog.show ? String(app.featureStyleScopeDialog.existingCaptionColor || '').toLowerCase() : 'no dialog';
    if (app.featureStyleScopeDialog.show) await app.handleFeatureStyleScopeChoice('cancel');
    app.clickedFeature = null;
    return color;
  });
  const paletteCds = (name) => page.evaluate((palette) => String(
    window.__GBDRAW_APP__.paletteDefinitions[palette || window.__GBDRAW_APP__.selectedPalette].CDS
  ).toLowerCase(), name);
  const expectShown = async (color, label) => {
    await settleLive(page);
    await expect.poll(() => legendRowColor(page, 'CDS'), { message: label }).toBe(color);
    expect.soft(await offered(), label).toBe(color);
  };
  const original = await paletteCds();
  await expectShown(original, 'Generate');

  await evaluateWithRetainedPromise(page, async (color) => {
    const { state } = await import('/gbdraw/web/js/state.js');
    await window.__GBDRAW_HISTORY__.runUndoable('Change color', async () => {
      state.activeDrawing().currentColors.value = { ...state.activeDrawing().currentColors.value, CDS: color };
      await window.Vue.nextTick();
    });
  }, USER_CDS);
  await expectShown(USER_CDS, 'Default colors edit');
  await generate(page);
  await expectShown(USER_CDS, 'Generate after the edit');

  await page.locator('h4[aria-label="DEFAULT COLORS"]').locator('..').getByRole('button', { name: 'Reset', exact: true }).click();
  await page.getByRole('dialog', { name: 'Reset default colors' })
    .getByRole('button', { name: "Reset to the palette's colors", exact: true }).click();
  await expectShown(original, 'Default colors Reset');

  const next = await page.evaluate(() => window.__GBDRAW_APP__.paletteNames.find((name) => (
    name !== window.__GBDRAW_APP__.selectedPalette
    && window.__GBDRAW_APP__.paletteDefinitions[name].CDS.toLowerCase() !== window.__GBDRAW_APP__.paletteDefinitions[window.__GBDRAW_APP__.selectedPalette].CDS.toLowerCase()
  )));
  await page.getByRole('combobox', { name: 'Palette', exact: true }).selectOption(next);
  await expect.poll(async () => (await facts(page)).palette).toBe(next);
  const nextCds = await paletteCds(next);
  await expectShown(nextCds, 'palette switch');

  await renameRow(page, 'CDS', 'Proteins');
  await settleLive(page);
  expect.soft(await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const { getFeatureFillElements } = await import('/gbdraw/web/js/app/feature-editor/svg-actions.js');
    const svgId = app.extractedFeatures.find((feature) => feature.locus_tag === 'FL1')?.svg_id;
    return String(getFeatureFillElements(document.querySelector('.origin-top svg'), svgId)[0]?.getAttribute('fill') || '').toLowerCase();
  }), 'a renamed CDS row keeps its features in the color it shows').toBe(nextCds);
  await expectLiveEqualsGenerate(page, { label: 'rename of a palette row after a palette switch' });
});
