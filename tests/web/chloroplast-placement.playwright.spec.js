const { test, expect } = require('@playwright/test');
const fs = require('node:fs/promises');
const { openApp } = require('./helpers/app-lifecycle.cjs');

test('chloroplast multipart placement survives strand and label changes and session restore', async ({ page, browser }, testInfo) => {
  test.setTimeout(240000);
  const external = [];
  const blockExternal = (route) => {
    const url = new URL(route.request().url());
    if (['127.0.0.1', 'localhost'].includes(url.hostname)) return route.continue();
    external.push(url.origin);
    return route.abort();
  };
  await page.route('**/*', blockExternal);
  page.on('dialog', (dialog) => dialog.dismiss());
  await page.setViewportSize({ width: 1600, height: 1000 });
  await openApp(page);
  await page.locator('input[accept^=".json,"]').setInputFiles(
    'gbdraw/web/gallery/sessions/tobacco-chloroplast.gbdraw-session.json');
  await expect.poll(() => page.evaluate(() => !window.__GBDRAW_APP__.sessionImportPending
    && window.__GBDRAW_APP__.results.length === 1), { timeout: 180000 }).toBe(true);

  const generate = async (name) => {
    const before = await page.evaluate(async () => (await import('./js/state.js')).state.resultGenerationKey.value);
    await page.getByRole('button', { name: 'Generate Diagram', exact: true }).click();
    await expect.poll(() => page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      return { key: state.resultGenerationKey.value, processing: state.processing.value, error: state.errorLog.value };
    }), { timeout: 180000 }).toMatchObject({ key: before + 1, processing: false, error: null });
    const svg = await page.evaluate(() => window.__GBDRAW_APP__.results[0].content);
    await fs.writeFile(testInfo.outputPath(`${name}.svg`), svg);
    // Check the actual SVG after the Web editor's composition reflow as well.
    const overlap = await page.evaluate(() => {
      const svg = document.querySelector('svg[data-gbdraw-composition]');
      const legend = svg.querySelector('[data-gbdraw-composition-role="legend"]').getBoundingClientRect();
      const axis = svg.querySelector('[id^="Axis"] circle').getBoundingClientRect();
      return legend.left < axis.right && legend.right > axis.left
        && legend.top < axis.bottom && legend.bottom > axis.top;
    });
    expect(overlap).toBe(false);
    for (let step = 0; step < 3; step += 1) {
      await page.getByRole('button', { name: 'Zoom out', exact: true }).click();
    }
    await page.screenshot({ path: testInfo.outputPath(`${name}.png`) });
    return svg;
  };
  const place = async (gene, side) => {
    await page.evaluate((gene) => {
      const app = window.__GBDRAW_APP__;
      const feature = app.extractedFeatures.find((f) => f.type === 'CDS' && f.qualifiers.gene?.[0] === gene);
      if (!feature || feature.location_parts.length < 2) throw new Error(`Multipart feature missing: ${gene}`);
      app.openFeatureEditorFromList(feature);
    }, gene);
    const placement = page.getByRole('combobox', { name: 'Feature placement', exact: true });
    await expect(placement.locator(`option[value=${side}]`)).toHaveJSProperty('disabled', false);
    await placement.selectOption(side);
    await page.getByRole('button', { name: 'Close feature popup', exact: true }).click();
  };
  // The refreshed Gallery session exposes placement immediately on import.
  await place('clpP', 'inward');
  await place('clpP', 'auto');
  await generate('01-loaded-regenerated');
  await place('clpP', 'inward');
  await place('rpl16', 'inward');
  await place('rpoC1', 'inward');
  await place('petB', 'outward');
  await generate('02-separated-placement');
  const strands = page.getByRole('checkbox', { name: 'Separate Strands', exact: true });
  await strands.uncheck();
  await generate('03-combined-both');
  await page.locator('summary[aria-label="Labels"]').click();
  await page.locator('#circular-label-mode').selectOption('out');
  await generate('04-combined-outer');
  await strands.check();
  await generate('05-separated-outer');
  await page.locator('#circular-label-mode').selectOption('both');
  await generate('06-separated-both');
  const expected = await page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    return JSON.parse(JSON.stringify(state.featurePlacementOverrides));
  });
  const pending = page.waitForEvent('download');
  await page.getByRole('button', { name: 'Save Session', exact: true }).click();
  const download = await pending;
  const saved = testInfo.outputPath(download.suggestedFilename());
  await download.saveAs(saved);
  const restoredContext = await browser.newContext({ viewport: { width: 1600, height: 1000 } });
  await restoredContext.route('**/*', blockExternal);
  page = await restoredContext.newPage();
  page.on('dialog', (dialog) => dialog.dismiss());
  await openApp(page);
  await page.locator('input[accept^=".json,"]').setInputFiles(saved);
  await expect.poll(() => page.evaluate(() => !window.__GBDRAW_APP__.sessionImportPending
    && window.__GBDRAW_APP__.results.length === 1), { timeout: 180000 }).toBe(true);
  await place('clpP', 'inward');
  await generate('07-restored');
  expect(await page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    return JSON.parse(JSON.stringify(state.featurePlacementOverrides));
  })).toEqual(expected);
  // The same restored result remains usable at a narrow viewport.
  await page.setViewportSize({ width: 390, height: 844 });
  await expect(page.getByRole('button', { name: 'Generate Diagram', exact: true })).toBeEnabled();
  await page.screenshot({ path: testInfo.outputPath('08-mobile.png') });
  expect(external).toEqual([]);
  await restoredContext.close();
});
