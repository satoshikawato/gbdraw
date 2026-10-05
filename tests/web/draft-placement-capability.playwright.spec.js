const { test, expect } = require('@playwright/test');
const fs = require('node:fs/promises');
const { readFileSync } = require('node:fs');
const { gunzipSync } = require('node:zlib');
const { openApp, reveal } = require('./helpers/app-lifecycle.cjs');
const { load, generate, switchMode, download } = require('./helpers/mode-transition.cjs');

for (const mode of ['circular', 'linear']) {
  test(`${mode} draft placement capability survives history and dirty session restore`, async ({ browser, baseURL }, testInfo) => {
    test.setTimeout(240000);
    const external = [];
    let context;
    let page;
    const load = async (file) => {
      if (context) await context.close();
      context = await browser.newContext({ viewport: { width: 1600, height: 1000 } });
      await context.route('**/*', (route) => {
        const url = new URL(route.request().url());
        if (url.origin === new URL(baseURL).origin) return route.continue();
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
        await expect(await reveal(page.locator('#circular-track-preset'))).toHaveValue(changed ? 'tuckin' : 'middle');
      } else {
        await expect(await reveal(page.getByRole('checkbox', { name: 'Separate Strands', exact: true, includeHidden: true }))).toBeChecked({ checked: !changed });
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
        const preset = await reveal(page.locator('#circular-track-preset'));
        await preset.focus();
        await preset.selectOption('tuckin');
        await preset.press('Tab');
      } else {
        await (await reveal(page.getByRole('checkbox', { name: 'Separate Strands', exact: true, includeHidden: true }))).uncheck();
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

// RC-4: Feature placement rows are keyed by mode-specific record keys, so each
// mode keeps its own rows (R2) and the request carries only its own (OV-08).
// An unsupported lane is a classified failure that names the feature (OV-09).
const HMMT_SESSION = 'gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json';
const placeLane = (page, side) => page.evaluate(async (lane) => {
  const app = window.__GBDRAW_APP__;
  const feature = app.extractedFeatures.find((item) => item.type === 'CDS');
  await app.featurePlacementActions.setPlacement([feature], lane);
  return { key: JSON.stringify([feature.scope, feature.record_key, feature.biological_feature_id]), product: feature.product };
}, side);
const placeOutward = (page) => placeLane(page, 'outward');
const placementState = (page, key) => page.evaluate(async (overrideKey) => {
  const { state } = await import('./js/state.js');
  const { getCommittedCanonicalRenderRequest } = await import('./js/services/config.js');
  const app = window.__GBDRAW_APP__;
  const feature = app.extractedFeatures.find((item) => item.type === 'CDS');
  return { row: state.featurePlacementOverrides[overrideKey]?.placement?.side || null,
    committed: getCommittedCanonicalRenderRequest().diagramOptions.featurePlacements.map((row) => row.placement.side || 'main'),
    pending: app.recordDisplayControls.hasPendingChanges.value,
    popup: feature ? app.featurePlacementActions.valueFor(feature) : null };
}, key);

test('a Circular lane placement waits in Circular while Linear generates (OV-08)', async ({ browser }, testInfo) => {
  test.setTimeout(300000);
  let page = await load(browser, HMMT_SESSION);
  try {
    await generate(page);
    const { key } = await placeOutward(page);
    await generate(page);
    expect(await placementState(page, key)).toEqual({ row: 'outward', committed: ['outward'], pending: false, popup: 'outward' });

    await switchMode(page, 'linear');
    await page.evaluate(async (text) => {
      window.__GBDRAW_APP__.setLinearSeqPrimaryFile(0, 'gb', new File([text], 'HmmtDNA.gbk', { type: 'text/plain', lastModified: 1000 }));
      await window.Vue.nextTick();
    }, readFileSync('tests/test_inputs/HmmtDNA.gbk', 'utf8'));
    await generate(page);
    expect(await placementState(page, key)).toEqual({ row: 'outward', committed: [], pending: false, popup: 'auto' });

    await switchMode(page, 'circular');
    await generate(page);
    expect(await placementState(page, key)).toEqual({ row: 'outward', committed: ['outward'], pending: false, popup: 'outward' });

    // A Session saved in Linear keeps the Circular row for a later return.
    await switchMode(page, 'linear');
    await generate(page);
    const saved = testInfo.outputPath('linear-with-circular-placement.gbdraw-session.json');
    await download(page, 'Save Session', saved);
    await page.context().close();
    page = await (await browser.newContext({ viewport: { width: 1600, height: 1000 } })).newPage();
    await openApp(page);
    await page.locator('input[accept^=".json,"]').setInputFiles(saved);
    const loaded = () => page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      const app = window.__GBDRAW_APP__;
      return { done: !app.sessionImportPending && (app.extractedFeatures.length > 0 || Boolean(state.errorLog.value)),
        error: state.errorLog.value?.summary || null, mode: app.mode };
    });
    await expect.poll(async () => (await loaded()).done, { timeout: 180000 }).toBe(true);
    expect(await loaded()).toEqual({ done: true, error: null, mode: 'linear' });
    expect((await placementState(page, key)).row).toBe('outward');
  } finally {
    await page.context().close();
  }
});

// Q3 (Owner, 2026-10-04): "Track type などを変える時点で "Reset N placements to
// Auto" / 変更を取り消す を選ばせる。" A layout edit that would leave lane
// placements undrawable asks first; Reset is one History step (R10, R11).
const placementDialog = (page) => page.getByRole('dialog', { name: 'Reset Feature placements?', exact: true });
const layoutState = (page, key) => page.evaluate(async (overrideKey) => {
  const { state } = await import('./js/state.js');
  return { trackType: state.form.track_type, separate: state.form.separate_strands,
    row: state.featurePlacementOverrides[overrideKey]?.placement?.side || null,
    undo: window.__GBDRAW_HISTORY__.getUndoCount() };
}, key);
const committedPlacements = (page) => page.evaluate(async () => {
  const { getCommittedCanonicalRenderRequest } = await import('./js/services/config.js');
  return getCommittedCanonicalRenderRequest().diagramOptions.featurePlacements;
});

test('a Track Preset change that drops a lane placement asks first (OV-09)', async ({ browser }, testInfo) => {
  test.setTimeout(240000);
  const page = await load(browser, HMMT_SESSION);
  try {
    const { key } = await placeOutward(page);
    const placed = await layoutState(page, key);
    expect(placed).toMatchObject({ trackType: 'middle', row: 'outward' });
    const preset = await reveal(page.locator('#circular-track-preset'));
    const dialog = placementDialog(page);
    const reset = dialog.getByRole('button', { name: 'Reset 1 placement to Auto', exact: true });
    const cancel = dialog.getByRole('button', { name: 'Cancel change', exact: true });

    // Cancel and Escape keep the preset and the placement and record no step.
    await preset.selectOption('spreadout');
    await expect(dialog).toBeVisible();
    await expect(dialog).toContainText('Changing Track Preset to Spreadout leaves 1 Feature placement without its lane.');
    await expect(reset).toBeFocused();
    await expect(preset).toHaveValue('middle');
    await page.screenshot({ path: testInfo.outputPath('track-preset-dialog.png') });
    await page.keyboard.press('Shift+Tab');
    await expect(cancel).toBeFocused();
    await page.keyboard.press('Tab');
    await expect(reset).toBeFocused();
    await cancel.click();
    await expect(dialog).toBeHidden();
    await expect(preset).toBeFocused();
    await expect(preset).toHaveValue('middle');
    expect(await layoutState(page, key)).toEqual(placed);
    await preset.selectOption('tuckin');
    await expect(reset).toBeFocused();
    await page.keyboard.press('Escape');
    await expect(dialog).toBeHidden();
    await expect(preset).toHaveValue('middle');
    expect(await layoutState(page, key)).toEqual(placed);

    // Reset applies the preset and removes the row as one undoable step.
    await preset.selectOption('spreadout');
    await reset.click();
    await expect(dialog).toBeHidden();
    await expect(preset).toHaveValue('spreadout');
    const applied = { ...placed, trackType: 'spreadout', row: null, undo: placed.undo + 1 };
    expect(await layoutState(page, key)).toEqual(applied);
    await page.getByRole('button', { name: /^Undo/ }).first().click();
    await expect(preset).toHaveValue('middle');
    expect(await layoutState(page, key)).toEqual({ ...placed, undo: placed.undo });
    await page.getByRole('button', { name: /^Redo/ }).first().click();
    await expect(preset).toHaveValue('spreadout');
    expect(await layoutState(page, key)).toEqual(applied);
    await expect(dialog).toHaveCount(0);
    await generate(page);
    expect(await committedPlacements(page)).toEqual([]);
  } finally {
    await page.context().close();
  }
});

test('Separate Strands on with an Above lane placement asks first (OV-09)', async ({ browser }) => {
  test.setTimeout(240000);
  const page = await load(browser, 'gbdraw/web/gallery/sessions/lambda_basic_linear.gbdraw-session.json');
  try {
    const strands = await reveal(page.getByRole('checkbox', { name: 'Separate Strands', exact: true, includeHidden: true }));
    await strands.click();
    await expect(strands).not.toBeChecked();
    const { key } = await placeLane(page, 'above');
    const placed = await layoutState(page, key);
    expect(placed).toMatchObject({ separate: false, row: 'above' });
    const dialog = placementDialog(page);
    await strands.click();
    await expect(dialog).toContainText('Changing Separate Strands to On leaves 1 Feature placement without its lane.');
    await expect(strands).not.toBeChecked();
    await dialog.getByRole('button', { name: 'Cancel change', exact: true }).click();
    await expect(strands).not.toBeChecked();
    expect(await layoutState(page, key)).toEqual(placed);
    await strands.click();
    await dialog.getByRole('button', { name: 'Reset 1 placement to Auto', exact: true }).click();
    await expect(dialog).toBeHidden();
    await expect(strands).toBeChecked();
    expect(await layoutState(page, key)).toEqual({ ...placed, separate: true, row: null, undo: placed.undo + 1 });
    await generate(page);
    expect(await committedPlacements(page)).toEqual([]);
  } finally {
    await page.context().close();
  }
});

// A path that is not a user setting edit, here a Session load, asks nothing;
// Generate keeps the classified diagnostic as the safety net (R6).
test('an unsupported lane from a loaded Session names the feature on Generate (OV-09)', async ({ browser }, testInfo) => {
  test.setTimeout(300000);
  let page = await load(browser, HMMT_SESSION);
  try {
    const { product } = await placeOutward(page);
    const saved = testInfo.outputPath('outward.gbdraw-session.json');
    const session = JSON.parse(gunzipSync(await download(page, 'Save Session', saved)));
    session.config.form.track_type = 'spreadout';
    const edited = testInfo.outputPath('outward-spreadout.gbdraw-session.json');
    await fs.writeFile(edited, JSON.stringify(session));
    await page.context().close();
    page = await load(browser, edited);
    await expect(await reveal(page.locator('#circular-track-preset'))).toHaveValue('spreadout');
    await expect(placementDialog(page)).toHaveCount(0);
    const before = await page.evaluate(async () => (await import('./js/state.js')).state.resultGenerationKey.value);
    await page.getByRole('button', { name: 'Generate Diagram', exact: true }).click();
    await expect.poll(() => page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      return !state.processing.value && Boolean(state.errorLog.value);
    }), { timeout: 180000 }).toBe(true);
    const error = await page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      const { code, context, summary } = state.errorLog.value;
      return { code, reason: context.reason, summary, generation: state.resultGenerationKey.value };
    });
    expect(error).toMatchObject({ code: 'FEATURE_PLACEMENT', reason: 'SPLIT_LANES', generation: before });
    expect(error.summary).toContain(`Feature: ${product}.`);
    expect(error.summary).toContain('Feature placement to Auto or Main');
    await expect(page.getByText(`Feature: ${product}.`, { exact: false }).first()).toBeVisible();
  } finally {
    await page.context().close();
  }
});
