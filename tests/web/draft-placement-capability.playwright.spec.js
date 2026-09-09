const { test, expect } = require('@playwright/test');
const fs = require('node:fs/promises');
const { openApp } = require('./helpers/app-lifecycle.cjs');

for (const mode of ['circular', 'linear']) {
  test(`@pr-smoke ${mode} draft placement capability survives history and dirty session restore`, async ({ browser }, testInfo) => {
    test.setTimeout(240000);
    const external = [];
    let context;
    let page;
    const load = async (file) => {
      if (context) await context.close();
      context = await browser.newContext({ viewport: { width: 1600, height: 1000 } });
      await context.route('**/*', (route) => {
        const url = new URL(route.request().url());
        if (url.origin === 'http://127.0.0.1:4173') return route.continue();
        external.push(url.origin);
        return route.abort();
      });
      page = await context.newPage();
      page.on('dialog', (dialog) => dialog.dismiss());
      await openApp(page);
      await page.locator('input[accept^=".json,"]').setInputFiles(file);
      await expect.poll(() => page.evaluate(() => !window.__GBDRAW_APP__.sessionImportPending
        && window.__GBDRAW_APP__.extractedFeatures.length > 0), { timeout: 180000 }).toBe(true);
    };
    const artifact = () => page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      return { results: window.__GBDRAW_APP__.results.map((result) => result.content),
        geometry: JSON.parse(JSON.stringify(state.trackSlotResolvedGeometry.value)) };
    });
    const generate = async () => {
      const key = await page.evaluate(async () => (await import('./js/state.js')).state.resultGenerationKey.value);
      await page.getByRole('button', { name: 'Generate Diagram', exact: true }).click();
      await expect.poll(() => page.evaluate(async () => {
        const { state } = await import('./js/state.js');
        return { key: state.resultGenerationKey.value, processing: state.processing.value,
          error: state.errorLog.value };
      }), { timeout: 180000 }).toEqual({ key: key + 1, processing: false, error: null });
    };
    const popup = async () => {
      if (!await page.locator('.right-drawer').isVisible()) await page.locator('.drawer-toggle').click();
      await page.locator('.right-drawer').getByRole('button', { name: 'Edit', exact: true }).first().click();
      return page.getByRole('combobox', { name: 'Feature placement', exact: true });
    };
    const closePopup = async () => {
      await page.getByRole('button', { name: 'Close feature popup', exact: true }).click();
      await page.locator('.drawer-toggle').click();
    };
    const sides = mode === 'circular' ? ['outward', 'inward'] : ['above', 'below'];
    const checkChoices = async (enabled) => {
      const placement = await popup();
      for (const side of sides) {
        await expect(placement.locator(`option[value=${side}]`)).toHaveJSProperty('disabled', !enabled);
      }
      for (const value of ['auto', 'main']) {
        await expect(placement.locator(`option[value=${value}]`)).toHaveJSProperty('disabled', false);
      }
      if (!enabled) {
        // An action caller cannot bypass the disabled DOM option.
        const rejected = await page.evaluate((sides) => {
          const app = window.__GBDRAW_APP__;
          const feature = app.extractedFeatures[0];
          const before = app.featurePlacementActions.valueFor(feature);
          const errors = sides.map((side) => {
            try { app.featurePlacementActions.setPlacement([feature], side); return ''; }
            catch (error) { return error.message; }
          });
          return { errors, before, after: app.featurePlacementActions.valueFor(feature) };
        }, sides);
        expect(rejected.errors.every((message) => message.includes('Unavailable'))).toBe(true);
        expect(rejected.after).toBe(rejected.before);
      }
      await closePopup();
    };
    const assertControl = async (changed) => {
      if (mode === 'circular') {
        await expect(page.locator('#circular-track-preset')).toHaveValue(changed ? 'tuckin' : 'middle');
      } else {
        await expect(page.getByRole('checkbox', { name: 'Separate Strands', exact: true })).toBeChecked({ checked: !changed });
      }
    };
    try {
      const seed = mode === 'circular' ? 'HmmtDNA_basic_circular' : 'lambda_basic_linear';
      await load(`gbdraw/web/gallery/sessions/${seed}.gbdraw-session.json`);
      await assertControl(false);
      await generate();
      const baseline = await artifact();
      await checkChoices(mode === 'circular');
      const ordinary = await popup();
      await ordinary.selectOption(mode === 'circular' ? 'outward' : 'main');
      await ordinary.selectOption('auto');
      await closePopup();
      expect(await artifact()).toEqual(baseline);

      if (mode === 'circular') {
        const preset = page.locator('#circular-track-preset');
        await preset.focus();
        await preset.selectOption('tuckin');
        await preset.press('Tab');
      } else {
        await page.getByRole('checkbox', { name: 'Separate Strands', exact: true }).uncheck();
      }
      await assertControl(true);
      await checkChoices(mode === 'linear');
      expect(await artifact()).toEqual(baseline);
      await page.getByRole('button', { name: /^Undo/ }).first().click();
      await assertControl(false);
      await checkChoices(mode === 'circular');
      expect(await artifact()).toEqual(baseline);
      await page.getByRole('button', { name: /^Redo/ }).first().click();
      await assertControl(true);
      await checkChoices(mode === 'linear');
      expect(await artifact()).toEqual(baseline);

      const pending = page.waitForEvent('download');
      await page.getByRole('button', { name: 'Save Session', exact: true }).click();
      const download = await pending;
      const saved = testInfo.outputPath(download.suggestedFilename());
      await download.saveAs(saved);
      await load(saved);
      await assertControl(true);
      expect(await artifact()).toEqual(baseline);
      await checkChoices(mode === 'linear');
      const placement = await popup();
      const valid = mode === 'circular' ? 'main' : 'above';
      await placement.selectOption(valid);
      await expect(placement).toHaveValue(valid);
      await page.screenshot({ path: testInfo.outputPath('dirty-restored-placement.png') });
      await closePopup();
      expect(await artifact()).toEqual(baseline);
      await generate();
      await assertControl(true);
      await checkChoices(mode === 'linear');
      const generated = await artifact();
      expect(generated.results).not.toEqual(baseline.results);
      expect(generated.geometry.records[0].featurePlacementTargets.map((target) => target.side || target.kind))
        .toEqual(mode === 'circular' ? ['main'] : ['main', ...sides]);
      const finalPlacement = await popup();
      await expect(finalPlacement).toHaveValue(valid);
      await closePopup();
      await fs.writeFile(testInfo.outputPath('generated.svg'), generated.results[0]);
      for (let step = 0; step < 3; step += 1) {
        await page.getByRole('button', { name: 'Zoom out', exact: true }).click();
      }
      await page.screenshot({ path: testInfo.outputPath('generated.png') });
      expect(external).toEqual([]);
    } finally {
      if (context) await context.close();
    }
  });
}
