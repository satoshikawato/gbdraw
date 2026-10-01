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
const fs = require('node:fs/promises');
const { download, generate, load } = require('./helpers/mode-transition.cjs');

test.describe.configure({ retries: 0 });

test('batch drawer lists and edits the features of the displayed Result', async ({ page }) => {
  test.setTimeout(300_000);
  await openBatch(page);
  await selectResult(page, 1);
  await page.locator('.drawer-toggle').click();
  const rows = await page.evaluate(() => window.__GBDRAW_APP__.visibleFeatureRows.map((row) => row.record_id));
  expect(rows.length).toBeGreaterThan(0);
  expect([...new Set(rows)]).toEqual(['TESTB']);
  await page.locator('.right-drawer').getByRole('button', { name: 'Edit', exact: true }).first().click();
  await expect(page.locator('.feature-popup')).toBeVisible();
  expect(await page.evaluate(() => window.__GBDRAW_APP__.clickedFeature?.feat?.record_id)).toBe('TESTB');
});

test('batch record-wide color and visibility scopes reach the other Result when it is displayed', async ({ page }) => {
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

// D-07 (PD-OI-062) Must preserve: the displayed Result's export, Save -> Load
// -> Result selection, and Undo reaching every Result.
test('batch edits shown on another Result reach its export, its Session, and Undo', async ({ page, browser }, testInfo) => {
  test.setTimeout(600_000);
  await openBatch(page);
  const undoCountBefore = await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
  const original = await featurePresentation(page, ['TESTA_0001', 'TESTA_0004', 'TESTB_0001', 'TESTB_0004']);
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
    await app.deleteLegendEntry(app.legendEntries.findIndex((entry) => entry.caption === 'GC content'));
  });
  await settle(page);
  await selectResult(page, 1);
  const shown = await featurePresentation(page, ['TESTB_0001', 'TESTB_0004']);
  expect(shown.TESTB_0001.results[1]).toBe('#ff0000|shown');
  expect(shown.TESTB_0004.results[1]).toMatch(/\|hidden$/);
  expect(await legendCaptions(page, { source: 'result' })).not.toContain('GC content');

  const [exportDownload] = await Promise.all([
    page.waitForEvent('download'),
    page.evaluate(() => window.__GBDRAW_APP__.downloadSVG())
  ]);
  const exported = await fs.readFile(await exportDownload.path(), 'utf8');
  expect(exported).not.toContain('data-legend-key="GC content"');
  expect(exported).toContain('#ff0000');

  await page.evaluate(() => { window.__GBDRAW_APP__.sessionTitle = 'batch-edits'; });
  const sessionFile = testInfo.outputPath('batch-edits.gbdraw-session.json.gz');
  await download(page, 'Save Session', sessionFile);
  const fresh = await load(browser, sessionFile);
  try {
    await fresh.locator('h2 select').selectOption({ index: 1 });
    await expect.poll(async () => (await featurePresentation(fresh, ['TESTB_0001'])).TESTB_0001.mounted, {
      timeout: 60_000
    }).toBe('#ff0000|shown');
    await settle(fresh);
    expect((await featurePresentation(fresh, ['TESTB_0004'])).TESTB_0004.mounted).toMatch(/\|hidden$/);
    expect(await legendCaptions(fresh)).not.toContain('GC content');
  } finally {
    await fresh.context().close();
  }

  // Undo back to the Generate also restores the Result selection of that
  // step; each Result then shows its original look when displayed.
  while (await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount()) > undoCountBefore) {
    await page.evaluate(() => window.__GBDRAW_HISTORY__.undo());
    await settle(page);
  }
  const displayed = await page.evaluate(() => window.__GBDRAW_APP__.selectedResultIndex);
  for (const index of [displayed, 1 - displayed]) {
    if (index !== displayed) await selectResult(page, index);
    const tags = index === 0 ? ['TESTA_0001', 'TESTA_0004'] : ['TESTB_0001', 'TESTB_0004'];
    const undone = await featurePresentation(page, tags);
    tags.forEach((tag) => {
      expect(undone[tag].mounted, `${tag} on Result ${index + 1}`).toBe(original[tag].results[index]);
      expect(undone[tag].results[index], `${tag} in Result ${index + 1}`).toBe(original[tag].results[index]);
    });
    expect(await legendCaptions(page)).toContain('GC content');
    expect(await legendCaptions(page, { source: 'result' })).toContain('GC content');
  }
});

test('Selected features annotations generate for a Circular multi-record batch', async ({ page }) => {
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
  // Python binds each explicit target to its own record's output only.
  expect(await page.evaluate(() => window.__GBDRAW_APP__.results.map((result) => (
    new DOMParser().parseFromString(result.content || '', 'image/svg+xml')
      .querySelectorAll('[data-gbdraw-annotation-id]').length > 0
  )))).toEqual([true, false]);
});

test('Reset fill after a canceled reset dialog uses the reset feature type default', async ({ page }) => {
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

// N-16 (D-11, PD-OI-066): a label reflow renders the committed Session with
// the current editor edits. Draft settings that apply on Generate must not
// reach the Result through it.
test('a label reflow leaves Applies on Generate draft settings out of the Result', async ({ browser }) => {
  test.setTimeout(600_000);
  const page = await load(browser, HMMT_SESSION);
  const resultFacts = () => page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    const ingestion = await import('./js/services/svg-result-ingestion.js');
    const { getCommittedCanonicalRenderRequest } = await import('./js/services/config.js');
    const result = state.results.value[state.selectedResultIndex.value];
    const root = new DOMParser().parseFromString(result.content, 'image/svg+xml').documentElement;
    return {
      identity: ingestion.getCommittedSvgResultRuntimeIdentity(result),
      processing: state.labelReflowProcessing.value,
      error: state.labelReflowLastError.value,
      definition: [...root.querySelectorAll('g[data-gbdraw-role="record-definition"] text')]
        .map((node) => `${Number(node.getAttribute('font-size'))}|${node.textContent}`),
      blockStrokeWidths: [...new Set([...root.querySelectorAll('[data-gbdraw-feature-part="block"]')]
        .map((node) => Number(node.getAttribute('stroke-width'))))],
      labels: [...root.querySelectorAll('text')].map((node) => node.textContent),
      request: JSON.stringify(getCommittedCanonicalRenderRequest())
    };
  });
  try {
    await generate(page);
    const before = await resultFacts();
    expect(before.definition.join(' ')).toContain('Homo sapiens');
    await page.evaluate(async () => {
      const app = window.__GBDRAW_APP__;
      Object.assign(app.form, { species: 'Leaked species' });
      Object.assign(app.adv, { def_font_size: 31, block_stroke_width: 7 });
      app.autoLabelReflowEnabled = true;
      await window.Vue.nextTick();
      const feature = app.extractedFeatures.find((item) => item.type === 'CDS');
      await app.openFeatureEditorFromList(feature, null);
      app.clickedFeature.labelText = 'REFLOWED_LABEL';
      await app.updateClickedFeatureLabelText();
      app.clickedFeature = null;
    });
    await expect.poll(async () => {
      const facts = await resultFacts();
      return facts.identity !== before.identity && !facts.processing && facts.labels.includes('REFLOWED_LABEL');
    }, { timeout: 300_000 }).toBe(true);
    const after = await resultFacts();
    expect(after.error).toBeNull();
    expect(after.labels).toContain('REFLOWED_LABEL');
    expect(after.definition).toEqual(before.definition);
    expect(after.blockStrokeWidths).toEqual(before.blockStrokeWidths);
    expect(after.request).toBe(before.request);
  } finally {
    await page.context().close();
  }
});

// N-16 Owner-delegated: Enable Labels keeps its documented live effect. Its
// label selection is committed with the reflowed Result, so the next label
// reflow (from the committed Session) keeps the labels.
test('Enable Labels applies its label selection through the reflow and keeps it', async ({ browser }) => {
  test.setTimeout(600_000);
  const page = await load(browser, HMMT_SESSION);
  const labelState = () => page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    const { getCommittedCanonicalRenderRequest } = await import('./js/services/config.js');
    const result = state.results.value[state.selectedResultIndex.value];
    const root = new DOMParser().parseFromString(result.content, 'image/svg+xml').documentElement;
    return {
      processing: state.labelReflowProcessing.value,
      error: state.labelReflowLastError.value,
      labels: [...root.querySelectorAll('text[data-label-feature-id]')].map((node) => node.textContent),
      scope: getCommittedCanonicalRenderRequest()?.diagramOptions?.configOverrides?.['labels.circular.scope']
    };
  });
  const editLabel = (text) => page.evaluate(async (label) => {
    const app = window.__GBDRAW_APP__;
    const feature = app.extractedFeatures.find((item) => item.type === 'CDS');
    await app.openFeatureEditorFromList(feature, null);
    app.clickedFeature.labelText = label;
    const update = app.updateClickedFeatureLabelText();
    for (let attempt = 0; attempt < 100 && !app.globalLabelModeDialog.show; attempt += 1) {
      await new Promise((resolve) => setTimeout(resolve, 20));
    }
    if (app.globalLabelModeDialog.show) app.handleGlobalLabelModeChoice('show_all');
    await update;
    app.clickedFeature = null;
  }, text);
  try {
    await page.evaluate(() => { window.__GBDRAW_APP__.form.labels_mode = 'none'; });
    await generate(page);
    expect((await labelState()).labels).toEqual([]);
    expect((await labelState()).scope).toBe('none');
    await editLabel('ENABLED_LABEL');
    await expect.poll(async () => {
      const state = await labelState();
      return !state.processing && state.labels.includes('ENABLED_LABEL');
    }, { timeout: 300_000 }).toBe(true);
    expect((await labelState()).scope).toBe('outer');
    const enabledCount = (await labelState()).labels.length;
    expect(enabledCount).toBeGreaterThan(1);
    await page.evaluate(() => { window.__GBDRAW_APP__.autoLabelReflowEnabled = true; });
    await editLabel('SECOND_LABEL');
    await expect.poll(async () => {
      const state = await labelState();
      return !state.processing && state.labels.includes('SECOND_LABEL');
    }, { timeout: 300_000 }).toBe(true);
    const after = await labelState();
    expect(after.error).toBeNull();
    expect(after.labels.length).toBe(enabledCount);
  } finally {
    await page.context().close();
  }
});
