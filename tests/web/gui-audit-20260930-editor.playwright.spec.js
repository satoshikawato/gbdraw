// Web GUI audit 2026-09-30: feature editor, batch Results and legend edits.
// A test marked test.fail(true, '<ID>') asserts the correct behavior of a
// current defect; the PR that fixes the audit ID removes the mark.
const { test, expect } = require('@playwright/test');
const { generateAndWaitForResult } = require('./helpers/app-lifecycle.cjs');
const {
  HMMT_SESSION,
  featurePresentation,
  legendCaptions,
  loadSessionFile,
  openBatch,
  openFresh,
  selectResult,
  settle,
  switchMode
} = require('./helpers/audit-browser.cjs');
const { generate, load } = require('./helpers/mode-transition.cjs');

test.describe.configure({ retries: 0 });

test('batch drawer lists and edits the features of the displayed Result', async ({ page }) => {
  test.fail(true, 'FE-03');
  test.setTimeout(300_000);
  await openBatch(page);
  await selectResult(page, 1);
  await page.locator('.drawer-toggle').click();
  const rows = await page.evaluate(() => window.__GBDRAW_APP__.visibleFeatureRows.map((row) => row.record_id));
  expect(rows.length).toBeGreaterThan(0);
  expect([...new Set(rows)]).toEqual(['TESTB']);
  await page.locator('.right-drawer').getByRole('button', { name: 'Edit', exact: true }).first().click();
  await expect(page.locator('.feature-popup')).toBeVisible();
  expect(await page.evaluate(() => window.__GBDRAW_APP__.clickedFeature?.record_id)).toBe('TESTB');
});

test('batch record-wide color and visibility scopes reach the other Result when it is displayed', async ({ page }) => {
  test.fail(true, 'FE-02');
  test.setTimeout(300_000);
  await openBatch(page);
  await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const colored = app.extractedFeatures.find((feature) => feature.locus_tag === 'TESTA_0001');
    await app.openFeatureEditorFromList(colored, null);
    await app.requestFeatureColorChange(colored, '#ff0000');
    await app.handleFeatureStyleScopeChoice('annotationLabel');
    const hidden = app.extractedFeatures.find((feature) => feature.locus_tag === 'TESTA_0004');
    await app.openFeatureEditorFromList(hidden, null);
    await app.updateClickedFeatureVisibility('off');
    await app.handleFeatureVisibilityScopeChoice('product');
    app.clickedFeature = null;
  });
  await settle(page);
  const edited = await featurePresentation(page, ['TESTA_0001', 'TESTA_0004']);
  expect(edited.TESTA_0001.mounted).toBe('#ff0000|shown');
  expect(edited.TESTA_0004.mounted).toMatch(/\|hidden$/);
  await selectResult(page, 1);
  const displayed = await featurePresentation(page, ['TESTB_0001', 'TESTB_0002', 'TESTB_0004']);
  expect(displayed.TESTB_0001.mounted).toBe('#ff0000|shown');
  expect(displayed.TESTB_0002.mounted).toBe('#ff0000|shown');
  expect(displayed.TESTB_0004.mounted).toMatch(/\|hidden$/);
});

test('batch live legend deletion reaches the other Result when it is displayed', async ({ page }) => {
  test.fail(true, 'PV-09');
  test.setTimeout(300_000);
  await openBatch(page);
  expect(await legendCaptions(page)).toContain('GC content');
  await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    await app.deleteLegendEntry(app.legendEntries.findIndex((entry) => entry.caption === 'GC content'));
  });
  await settle(page);
  expect(await legendCaptions(page)).not.toContain('GC content');
  await selectResult(page, 1);
  expect(await legendCaptions(page)).not.toContain('GC content');
});

test('Selected features annotations generate for a Circular multi-record batch', async ({ page }) => {
  test.fail(true, 'FE-05');
  test.setTimeout(300_000);
  await openBatch(page);
  await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    app.selectedFeatureIds = app.extractedFeatures.filter((feature) => feature.type === 'CDS')
      .slice(0, 2).map((feature) => feature.svg_id);
    await window.Vue.nextTick();
    app.addAnnotationSet();
    await app.addSelectedFeatureAnnotations(app.annotationSets[app.annotationSets.length - 1]);
    await window.Vue.nextTick();
  });
  await generateAndWaitForResult(page);
});

test('Reset fill after a canceled reset dialog uses the reset feature type default', async ({ page }) => {
  test.fail(true, 'FE-10');
  test.setTimeout(300_000);
  await openFresh(page);
  await loadSessionFile(page, HMMT_SESSION);
  const fill = (product) => page.evaluate((name) => {
    const app = window.__GBDRAW_APP__;
    const feature = app.extractedFeatures.find((item) => item.product === name);
    return [...app.svgContainer.querySelector('svg').querySelectorAll(`[data-gbdraw-feature-id="${CSS.escape(feature.svg_id)}"]`)]
      .map((element) => element.getAttribute('fill')).find((value) => value && value !== 'none');
  }, product);
  const openAndReset = (product) => page.evaluate(async (name) => {
    const app = window.__GBDRAW_APP__;
    await app.openFeatureEditorFromList(app.extractedFeatures.find((item) => item.product === name), null);
    await app.resetClickedFeatureFillColor();
  }, product);
  await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    Object.assign(app.newSpecRule, { feat: 'rRNA', qual: 'product', val: '^s-rRNA$', color: '#ff00ff', cap: 'small rRNA' });
    await app.addSpecificRule();
  });
  await settle(page);
  await expect.poll(() => fill('s-rRNA')).toBe('#ff00ff');
  await openAndReset('tRNA-Leu');
  await settle(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.resetColorDialog.show)).toBe(true);
  await page.evaluate(() => window.__GBDRAW_APP__.handleResetColorChoice('cancel'));
  await settle(page);
  await openAndReset('s-rRNA');
  await settle(page);
  const rrnaDefault = await page.evaluate(() => window.__GBDRAW_APP__.appliedPaletteColors.rRNA);
  await expect.poll(() => fill('s-rRNA'), { timeout: 10_000 }).toBe(rrnaDefault);
});

test('Redo of an Exact product hide hides the feature again', async ({ page }) => {
  test.fail(true, 'FE-04');
  test.setTimeout(300_000);
  await openFresh(page);
  await loadSessionFile(page, HMMT_SESSION);
  const id = await page.evaluate(() => window.__GBDRAW_APP__.extractedFeatures
    .find((feature) => feature.product === 'NADH dehydrogenase subunit 1').svg_id);
  const display = () => page.evaluate((featureId) => [...window.__GBDRAW_APP__.svgContainer.querySelector('svg')
    .querySelectorAll(`[data-gbdraw-feature-id="${CSS.escape(featureId)}"]`)]
    .some((element) => element.getAttribute('display') === 'none'), id);
  await page.evaluate(async (featureId) => {
    const app = window.__GBDRAW_APP__;
    await app.openFeatureEditorFromList(app.extractedFeatures.find((feature) => feature.svg_id === featureId), null);
    await app.updateClickedFeatureVisibility('off');
  }, id);
  await page.getByRole('button', { name: /Exact product/ }).click();
  await settle(page);
  expect(await display()).toBe(true);
  await page.evaluate(() => window.__GBDRAW_HISTORY__.undo());
  await settle(page);
  expect(await display()).toBe(false);
  await page.evaluate(() => window.__GBDRAW_HISTORY__.redo());
  await settle(page);
  expect(await display()).toBe(true);
});

test('a rename of a legend entry without features survives Generate', async ({ browser }) => {
  test.fail(true, 'PV-02');
  test.setTimeout(600_000);
  const page = await load(browser, HMMT_SESSION);
  try {
    await generate(page);
    await page.evaluate(async () => {
      const app = window.__GBDRAW_APP__;
      await app.renameLegendEntry(app.legendEntries.findIndex((entry) => entry.caption === 'GC content'), 'GC percent');
    });
    await settle(page);
    expect(await legendCaptions(page)).toContain('GC percent');
    await generate(page);
    expect(await legendCaptions(page, { source: 'result' })).toContain('GC percent');
    expect(await legendCaptions(page, { source: 'result' })).not.toContain('GC content');
  } finally {
    await page.context().close();
  }
});

test('legend order survives Generate', async ({ browser }) => {
  test.fail(true, 'PV-03');
  test.setTimeout(600_000);
  const page = await load(browser, HMMT_SESSION);
  try {
    await generate(page);
    if (!await page.evaluate(() => window.__GBDRAW_APP__.showRightDrawer)) await page.locator('.drawer-toggle').click();
    await page.locator('.right-drawer').getByRole('button', { name: 'Legend' }).click();
    await page.locator('.right-drawer').getByTitle('Sort Z-A', { exact: true }).click();
    await settle(page);
    const sorted = await legendCaptions(page);
    expect(sorted).toEqual([...sorted].sort((left, right) => right.localeCompare(left)));
    await generate(page);
    expect(await legendCaptions(page, { source: 'result' })).toEqual(sorted);
  } finally {
    await page.context().close();
  }
});

test('renaming a feature legend entry to an existing caption offers the conflict dialog', async ({ browser }) => {
  test.fail(true, 'PV-04');
  test.setTimeout(600_000);
  const page = await load(browser, HMMT_SESSION);
  try {
    await generate(page);
    await page.evaluate(async () => {
      const app = window.__GBDRAW_APP__;
      await app.renameLegendEntry(app.legendEntries.findIndex((entry) => entry.caption === 'tRNA'), 'rRNA')
        .catch(() => {});
    });
    await settle(page);
    const outcome = await page.evaluate(() => ({
      errorCode: window.__GBDRAW_APP__.errorLog?.code ?? null,
      dialog: Boolean(window.__GBDRAW_APP__.legendRenameDialog?.show)
    }));
    expect(outcome.errorCode).not.toBe('UNKNOWN');
    expect(outcome.dialog).toBe(true);
  } finally {
    await page.context().close();
  }
});

test('changing the Linear legend position after a legend-free Generate raises no page error', async ({ browser }) => {
  test.fail(true, 'GE-07');
  test.setTimeout(600_000);
  const page = await load(browser, 'gbdraw/web/gallery/sessions/lambda_basic_linear.gbdraw-session.json');
  const pageErrors = [];
  page.on('pageerror', (error) => pageErrors.push(String(error?.message || error)));
  try {
    const position = page.getByLabel('Legend position', { exact: true });
    for (const details of await position.locator('xpath=ancestor::details').all()) {
      if (await details.getAttribute('open') === null) await details.locator(':scope > summary').press('Enter');
    }
    await position.selectOption('none');
    await generate(page);
    await position.selectOption('left');
    await settle(page);
    expect(pageErrors).toEqual([]);
  } finally {
    await page.context().close();
  }
});

test('Escape that closes the Editor returns focus to the Editor toggle', async ({ browser }) => {
  test.fail(true, 'PV-11');
  test.setTimeout(600_000);
  const page = await load(browser, HMMT_SESSION);
  try {
    const toggle = page.locator('.drawer-toggle');
    await toggle.click();
    await expect(page.locator('.right-drawer')).toBeVisible();
    await page.locator('.right-drawer').getByRole('button', { name: 'Legend' }).click();
    await page.keyboard.press('Escape');
    await expect(page.locator('.right-drawer')).toBeHidden();
    expect(await toggle.evaluate((element) => element === document.activeElement)).toBe(true);
  } finally {
    await page.context().close();
  }
});

test('an Undo of a checkpoint edit keeps the feature catalog across a mode round trip', async ({ browser }) => {
  test.fail(true, 'SE-01');
  test.setTimeout(600_000);
  const page = await load(browser, HMMT_SESSION);
  const pageErrors = [];
  page.on('pageerror', (error) => pageErrors.push(String(error?.message || error)));
  try {
    const featureCount = await page.evaluate(() => window.__GBDRAW_APP__.extractedFeatures.length);
    await page.evaluate(async () => {
      const app = window.__GBDRAW_APP__;
      const feature = app.extractedFeatures.find((item) => item.type === 'CDS');
      await app.setFeatureColorValue(feature, '#123456');
    });
    await settle(page);
    await page.evaluate(() => window.__GBDRAW_HISTORY__.undo());
    await settle(page);
    await switchMode(page, 'linear');
    await switchMode(page, 'circular');
    expect(pageErrors).toEqual([]);
    expect(await page.evaluate(() => window.__GBDRAW_APP__.extractedFeatures.length)).toBe(featureCount);
  } finally {
    await page.context().close();
  }
});
