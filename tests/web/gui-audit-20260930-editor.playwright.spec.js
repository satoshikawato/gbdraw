// Web GUI audit 2026-09-30: feature editor, batch Results and legend edits.
// A test marked test.fail(true, '<ID>') asserts the correct behavior of a
// current defect; the PR that fixes the audit ID removes the mark.
const { test, expect } = require('@playwright/test');
const { evaluateWithRetainedPromise, generateAndWaitForResult } = require('./helpers/app-lifecycle.cjs');
const {
  HMMT_SESSION,
  featurePresentation,
  legendCaptions,
  loadSessionFile,
  openBatch,
  openFresh,
  openWithGenBank,
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
  await evaluateWithRetainedPromise(page, async () => {
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
  await evaluateWithRetainedPromise(page, async () => {
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
  await evaluateWithRetainedPromise(page, async () => {
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
    evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.downloadSVG())
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
    await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
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
  await evaluateWithRetainedPromise(page, async () => {
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
  const openAndReset = (product) => evaluateWithRetainedPromise(page, async (name) => {
    const app = window.__GBDRAW_APP__;
    await app.openFeatureEditorFromList(app.extractedFeatures.find((item) => item.product === name), null);
    await app.resetClickedFeatureFillColor();
  }, product);
  await evaluateWithRetainedPromise(page, async () => {
    const app = window.__GBDRAW_APP__;
    Object.assign(app.newSpecRule, { feat: 'rRNA', qual: 'product', val: '^s-rRNA$', color: '#ff00ff', cap: 'small rRNA' });
    await app.addSpecificRule();
  });
  await settle(page);
  await expect.poll(() => fill('s-rRNA')).toBe('#ff00ff');
  await openAndReset('tRNA-Leu');
  await settle(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.resetColorDialog.show)).toBe(true);
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.handleResetColorChoice('cancel'));
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
  await evaluateWithRetainedPromise(page, async (featureId) => {
    const app = window.__GBDRAW_APP__;
    await app.openFeatureEditorFromList(app.extractedFeatures.find((feature) => feature.svg_id === featureId), null);
    await app.updateClickedFeatureVisibility('off');
  }, id);
  await page.getByRole('button', { name: /Exact product/ }).click();
  await settle(page);
  expect(await display()).toBe(true);
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
  await settle(page);
  expect(await display()).toBe(false);
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.redo());
  await settle(page);
  expect(await display()).toBe(true);
});

test('a rename of a legend entry without features survives Generate', async ({ browser }) => {
  test.setTimeout(600_000);
  const page = await load(browser, HMMT_SESSION);
  try {
    await generate(page);
    await evaluateWithRetainedPromise(page, async () => {
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
  test.setTimeout(600_000);
  const page = await load(browser, HMMT_SESSION);
  try {
    await generate(page);
    await evaluateWithRetainedPromise(page, async () => {
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
    // Reset Settings also restores a legend side on the legend-free Result.
    await page.evaluate(() => { window.confirm = () => true; });
    await page.getByRole('button', { name: 'Reset Settings', exact: true }).click();
    await settle(page);
    expect(await page.evaluate(() => window.__GBDRAW_APP__.form.legend)).not.toBe('left');
    expect(pageErrors).toEqual([]);
  } finally {
    await page.context().close();
  }
});

// D-30 (PD-OI-084): the in-place Linear side move could not match the
// renderer (PV-10 measurement in docs/internal/web-gui-audit-20260930), so a
// Linear legend side applies on Generate and the Result is unchanged until then.
test('a Linear legend side change leaves the Result unchanged until Generate applies it', async ({ browser }) => {
  test.setTimeout(600_000);
  const page = await load(browser, 'gbdraw/web/gallery/sessions/lambda_basic_linear.gbdraw-session.json');
  try {
    await generate(page);
    const read = () => page.evaluate(() => {
      const app = window.__GBDRAW_APP__;
      const mounted = app.svgContainer.querySelector('svg');
      return {
        result: app.results[app.selectedResultIndex].content,
        mountedSide: JSON.parse(mounted.getAttribute('data-gbdraw-composition')).legendSide,
        viewBox: mounted.getAttribute('viewBox')
      };
    });
    const before = await read();
    const position = page.getByLabel('Legend position', { exact: true });
    for (const details of await position.locator('xpath=ancestor::details').all()) {
      if (await details.getAttribute('open') === null) await details.locator(':scope > summary').press('Enter');
    }
    await expect(page.getByText('Applies on Generate: legend position, swatch size, and font size.')).toBeVisible();
    await position.selectOption('top');
    await settle(page);
    expect(await read()).toEqual(before);
    await generate(page);
    const after = await read();
    expect(after.mountedSide).toBe('top');
    expect(after.viewBox).not.toBe(before.viewBox);
  } finally {
    await page.context().close();
  }
});

// W1b PV-02/PV-03: a horizontal legend keeps a sort and a featureless rename
// through Generate and a Session round trip.
test('a horizontal legend keeps a sort and a featureless rename through Generate and Load', async ({ browser }, testInfo) => {
  test.setTimeout(600_000);
  const page = await load(browser, HMMT_SESSION);
  try {
    await page.evaluate(() => { window.__GBDRAW_APP__.form.legend = 'top'; });
    await generate(page);
    if (!await page.evaluate(() => window.__GBDRAW_APP__.showRightDrawer)) await page.locator('.drawer-toggle').click();
    await page.locator('.right-drawer').getByRole('button', { name: 'Legend' }).click();
    await page.locator('.right-drawer').getByTitle('Sort Z-A', { exact: true }).click();
    await settle(page);
    await evaluateWithRetainedPromise(page, async () => {
      const app = window.__GBDRAW_APP__;
      await app.renameLegendEntry(app.legendEntries.findIndex((entry) => entry.caption === 'GC content'), 'GC percent');
    });
    await settle(page);
    const edited = await legendCaptions(page);
    expect(edited).toContain('GC percent');
    expect(edited.filter((caption) => caption !== 'GC percent'))
      .toEqual(edited.filter((caption) => caption !== 'GC percent').sort((left, right) => right.localeCompare(left)));
    await generate(page);
    expect(await legendCaptions(page, { source: 'result' })).toEqual(edited);
    const sessionFile = testInfo.outputPath('legend-order.gbdraw-session.json.gz');
    await download(page, 'Save Session', sessionFile);
    const fresh = await load(browser, sessionFile);
    try {
      expect(await legendCaptions(fresh)).toEqual(edited);
      await generate(fresh);
      expect(await legendCaptions(fresh, { source: 'result' })).toEqual(edited);
    } finally {
      await fresh.context().close();
    }
  } finally {
    await page.context().close();
  }
});

// D-09 (PD-OI-064): one canvas padding reaches every batch Result, once.
test('canvas padding reaches every batch Result through Generate and applies once', async ({ page }) => {
  test.setTimeout(600_000);
  await openBatch(page);
  const canvases = () => page.evaluate(() => window.__GBDRAW_APP__.results.map((result) => {
    const svg = new DOMParser().parseFromString(result.content, 'image/svg+xml').documentElement;
    return svg.getAttribute('viewBox').trim().split(/[\s,]+/).map(Number);
  }));
  const base = await canvases();
  await page.getByRole('button', { name: 'Toggle canvas padding controls' }).click();
  const right = page.locator('.preview-canvas-padding input[type=number]').nth(2);
  await right.fill('150');
  await right.press('Tab');
  await settle(page);
  for (let run = 0; run < 2; run += 1) {
    await generateAndWaitForResult(page);
    await settle(page);
    expect((await canvases()).map((box) => box[2])).toEqual(base.map((box) => box[2] + 150));
  }
  await selectResult(page, 1);
  expect(await page.evaluate(() => Number(window.__GBDRAW_APP__.svgContainer.querySelector('svg')
    .getAttribute('viewBox').trim().split(/[\s,]+/)[2]))).toBe(base[1][2] + 150);
  await page.locator('.preview-canvas-padding').getByRole('button', { name: 'Reset', exact: true }).click();
  await settle(page);
  await selectResult(page, 0);
  expect((await canvases()).map((box) => box[2])).toEqual(base.map((box) => box[2]));
});

test('Escape that closes the Editor returns focus to the Editor toggle', async ({ browser }) => {
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
    await evaluateWithRetainedPromise(page, async () => {
      const app = window.__GBDRAW_APP__;
      const feature = app.extractedFeatures.find((item) => item.type === 'CDS');
      await app.setFeatureColorValue(feature, '#123456');
    });
    await settle(page);
    await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
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
    await evaluateWithRetainedPromise(page, async () => {
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

// Owner request (2026-10-04): a Label visibility choice in the feature popup
// applies whatever Show Labels selects. The live reflow and Generate draw the
// same labels (N-16, PD-OI-066).
const labelEditorState = (page, scopePath) => page.evaluate(async (path) => {
  const { state } = await import('./js/state.js');
  const { getCommittedCanonicalRenderRequest } = await import('./js/services/config.js');
  const result = state.results.value[state.selectedResultIndex.value];
  const root = new DOMParser().parseFromString(result.content, 'image/svg+xml').documentElement;
  return {
    processing: state.labelReflowProcessing.value || state.processing.value,
    reflowError: state.labelReflowLastError.value,
    error: state.errorLog.value,
    labels: [...root.querySelectorAll('text[data-label-feature-id]')]
      .filter((node) => node.getAttribute('display') !== 'none')
      .map((node) => [node.getAttribute('data-label-feature-id'), node.textContent]),
    scope: getCommittedCanonicalRenderRequest()?.diagramOptions?.configOverrides?.[path]
  };
}, scopePath);

// `choice` answers the Label Not Shown dialog that a text edit on an
// unlabeled feature opens; `asked` reports whether it opened.
const applyPopupLabel = (page, featureId, text, visibility, choice = null) => evaluateWithRetainedPromise(page, async (edit) => {
  const app = window.__GBDRAW_APP__;
  const feature = app.extractedFeatures.find((item) => item.svg_id === edit.featureId);
  await app.openFeatureEditorFromList(feature, null);
  await window.Vue.nextTick();
  const hint = document.querySelector('.feature-popup')?.textContent
    .includes('This feature has no label in the current Result.');
  app.clickedFeature.labelVisibility = edit.visibility;
  if (edit.text !== null) app.clickedFeature.labelText = edit.text;
  await app.updateClickedFeatureLabelText();
  await window.Vue.nextTick();
  const asked = app.hiddenLabelTextDialog.show
    && Boolean([...document.querySelectorAll('h3')].find((node) => node.textContent === 'Label Not Shown'));
  if (asked && edit.choice) await app.handleHiddenLabelTextChoice(edit.choice);
  app.clickedFeature = null;
  return { hint, asked };
}, { featureId, text, visibility, choice });

test('Label visibility On shows only its label under Show Labels None, live and after Generate', async ({ browser }) => {
  test.setTimeout(600_000);
  const page = await load(browser, HMMT_SESSION);
  const labelState = () => labelEditorState(page, 'labels.circular.scope');
  try {
    await page.evaluate(() => { window.__GBDRAW_APP__.form.labels_mode = 'none'; });
    await generate(page);
    expect((await labelState()).labels).toEqual([]);
    const [forced, textOnly, shown] = await page.evaluate(() => window.__GBDRAW_APP__.extractedFeatures
      .filter((item) => item.type === 'CDS').slice(0, 3).map((item) => item.svg_id));
    expect(await applyPopupLabel(page, forced, 'FORCED_LABEL', 'on')).toEqual({ hint: true, asked: false });
    await expect.poll(async () => {
      const state = await labelState();
      return state.processing ? null : state.labels;
    }, { timeout: 300_000 }).toEqual([[forced, 'FORCED_LABEL']]);
    expect((await labelState()).reflowError).toBeNull();
    expect((await labelState()).scope).toBe('none');
    // Label text alone keeps Default visibility. Keep hidden leaves the feature
    // unlabeled; Show this label sets On. Generate draws the same labels.
    expect(await applyPopupLabel(page, textOnly, 'TEXT_ONLY_LABEL', 'default', 'text_only'))
      .toEqual({ hint: true, asked: true });
    expect(await applyPopupLabel(page, shown, 'SHOWN_LABEL', 'default', 'show'))
      .toEqual({ hint: true, asked: true });
    const expected = [[forced, 'FORCED_LABEL'], [shown, 'SHOWN_LABEL']]
      .sort(([left], [right]) => left.localeCompare(right));
    const sortedLabels = (state) => [...state.labels].sort(([left], [right]) => left.localeCompare(right));
    await expect.poll(async () => {
      const state = await labelState();
      return state.processing ? null : sortedLabels(state);
    }, { timeout: 300_000 }).toEqual(expected);
    await generate(page);
    const after = await labelState();
    expect(sortedLabels(after)).toEqual(expected);
    expect(after.scope).toBe('none');
  } finally {
    await page.context().close();
  }
});

test('Label visibility On and Off apply under First Record Only, live and after Generate', async ({ browser }) => {
  test.setTimeout(600_000);
  const page = await load(browser, 'gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json');
  const labelState = () => labelEditorState(page, 'labels.linear.scope');
  try {
    await generate(page);
    const before = await labelState();
    expect(before.scope).toBe('first');
    const neoR = await page.evaluate(() => window.__GBDRAW_APP__.extractedFeatures
      .find((item) => item.product === 'putative regulator, NeoR')?.svg_id);
    expect(neoR).toMatch(/_record_2$/);
    expect(before.labels.some(([featureId]) => featureId === neoR)).toBe(false);
    const [hidden] = before.labels.find(([featureId]) => /_record_1$/.test(featureId));

    expect(await applyPopupLabel(page, neoR, 'NeoR', 'on')).toEqual({ hint: true, asked: false });
    await expect.poll(async () => {
      const state = await labelState();
      return state.processing ? null : state.labels.filter(([featureId]) => featureId === neoR);
    }, { timeout: 300_000 }).toEqual([[neoR, 'NeoR']]);
    expect(await applyPopupLabel(page, hidden, null, 'off')).toEqual({ hint: false, asked: false });
    await expect.poll(async () => {
      const state = await labelState();
      return state.processing ? null : state.labels.some(([featureId]) => featureId === hidden);
    }, { timeout: 300_000 }).toBe(false);
    expect((await labelState()).reflowError).toBeNull();

    await generate(page);
    const after = await labelState();
    expect(after.scope).toBe('first');
    expect(after.labels.filter(([featureId]) => featureId === neoR)).toEqual([[neoR, 'NeoR']]);
    expect(after.labels.some(([featureId]) => featureId === hidden)).toBe(false);
    expect(after.labels.filter(([featureId]) => !/_record_1$/.test(featureId))).toEqual([[neoR, 'NeoR']]);
  } finally {
    await page.context().close();
  }
});

// FINDINGS OV-05 (color), PD-OI-069: Python matches a `hash` rule by the hash
// without the `__instance_` suffix, so a This feature only fill on a duplicated
// record uses that hash and survives Generate on every copy of the record.
test('a This feature only fill on a duplicated record survives Generate', async ({ page }) => {
  test.setTimeout(600_000);
  await openWithGenBank(page, `${__dirname}/../fixtures/b_dup_ids.gb`, () => {
    Object.assign(window.__GBDRAW_APP__.form, { labels_mode: 'out', multi_record_canvas: true });
  });
  await generateAndWaitForResult(page);
  await settle(page);
  const copies = await page.evaluate(() => window.__GBDRAW_APP__.extractedFeatures
    .filter((feature) => feature.type === 'CDS' && feature.locus_tag === 'TESTA_0002')
    .map((feature) => feature.svg_id).sort());
  expect(copies).toHaveLength(2);
  expect(copies[0]).toContain('__instance_');
  const fills = () => page.evaluate((ids) => {
    const root = window.__GBDRAW_APP__.svgContainer.querySelector('svg');
    const hidden = (element) => {
      for (let node = element; node && node.nodeType === 1; node = node.parentNode) {
        const style = node.getAttribute('style') || '';
        if (node.getAttribute('display') === 'none' || /display\s*:\s*none/.test(style)) return true;
      }
      return false;
    };
    return ids.map((id) => [...root.querySelectorAll(
      `path[data-gbdraw-feature-id="${CSS.escape(id)}"], path[data-gbdraw-rendered-feature-id="${CSS.escape(id)}"]`
    )].filter((element) => !hidden(element))
      .map((element) => element.getAttribute('fill')).find((value) => value && value !== 'none') || null);
  }, copies);
  await evaluateWithRetainedPromise(page, async (id) => {
    const app = window.__GBDRAW_APP__;
    await app.openFeatureEditorFromList(app.extractedFeatures.find((item) => item.svg_id === id), null);
    await app.updateClickedFeatureColor('#c83366');
    if (app.featureStyleScopeDialog.show) await app.handleFeatureStyleScopeChoice('single');
    app.clickedFeature = null;
  }, copies[0]);
  await settle(page);
  const rules = await page.evaluate(() => window.__GBDRAW_APP__.manualSpecificRules.map((rule) => rule.val));
  expect(rules).toEqual([copies[0].replace(/__instance_.*$/, '')]);
  const live = await fills();
  expect(live[0]).toBe('#c83366');
  await generateAndWaitForResult(page);
  await settle(page);
  expect((await fills())[0]).toBe('#c83366');
});

const FORCED_LABEL_FIXTURE = 'tests/fixtures/forced_label_underlay.gb';

const featureIdsByLocator = (page) => page.evaluate(() => {
  const features = window.__GBDRAW_APP__.extractedFeatures;
  const find = (predicate) => features.find(predicate)?.svg_id;
  return {
    fl1: find((item) => item.locus_tag === 'FL1'),
    fl2: find((item) => item.locus_tag === 'FL2'),
    repeat: find((item) => item.type === 'repeat_region'),
    dup: find((item) => item.product === 'dup beta')
  };
});

const waitForLabelReflow = (page) => expect.poll(async () => (
  await labelEditorState(page, 'labels.circular.scope')
).processing, { timeout: 300_000 }).toBe(false);

// OV-05 (design Q4): Label visibility On for one of two features with the same
// type and location is an identity row for that feature, so the live reflow and
// Generate draw its label only. Before identity rows the shared hash row drew
// neither label and Generate failed with LABEL_NOT_DRAWN.
test('a Label visibility On for one of two same-location features draws only its label', async ({ page }) => {
  test.setTimeout(600_000);
  await openWithGenBank(page, FORCED_LABEL_FIXTURE, () => {
    window.__GBDRAW_APP__.form.labels_mode = 'none';
  });
  await generate(page);
  const { dup } = await featureIdsByLocator(page);
  expect(dup).toContain('__instance_');
  const labelState = () => labelEditorState(page, 'labels.circular.scope');

  expect(await applyPopupLabel(page, dup, null, 'on')).toEqual({ hint: true, asked: false });
  await waitForLabelReflow(page);
  const live = await labelState();
  expect(live.reflowError).toBeNull();
  expect(live.labels.map(([featureId]) => featureId)).toEqual([dup]);

  await generate(page);
  expect(await labelState()).toMatchObject({ error: null, reflowError: null, labels: live.labels });
});

// Owner decisions Q1 and Q2 (2026-10-04; FINDINGS OV-06, OV-07, OV-10): an
// Apply of Label visibility On asks when the diagram cannot draw the label.
// Each choice is one History step and Cancel records none; Generate does not
// fail on an On that cannot be drawn, and draws it once it can.
const labelOnFacts = (page) => page.evaluate(async () => {
  const { state } = await import('./js/state.js');
  // Per-feature edits by the rendered ID of the feature in view.
  const edits = (field) => Object.fromEntries(state.extractedFeatures.value.flatMap((feature) => {
    const value = state.featureOverrides[JSON.stringify([feature.scope, feature.record_key, feature.biological_feature_id])]?.[field];
    return value === null || value === undefined ? [] : [[feature.svg_id, value]];
  }));
  return {
    undo: window.__GBDRAW_HISTORY__.getUndoCount(),
    redo: window.__GBDRAW_HISTORY__.getRedoCount(),
    labelVisibility: edits('labelVisibility'),
    labelText: edits('labelText'),
    featureVisibility: edits('featureVisibility')
  };
});

// Applies Label visibility On (and optional text) in the feature popup. The
// Apply settles after the dialog it opens is answered.
const startLabelOn = (page, featureId, text = null) => evaluateWithRetainedPromise(page, async (edit) => {
  const app = window.__GBDRAW_APP__;
  await app.openFeatureEditorFromList(app.extractedFeatures.find((item) => item.svg_id === edit.featureId), null);
  await window.Vue.nextTick();
  app.clickedFeature.labelVisibility = 'on';
  if (edit.text !== null) app.clickedFeature.labelText = edit.text;
  window.__labelOnApply = app.updateClickedFeatureLabelText();
}, { featureId, text });

const answerLabelOn = async (page, title, choice) => {
  const dialog = page.getByRole('dialog', { name: title, exact: true });
  await expect(dialog).toBeVisible({ timeout: 120_000 });
  await dialog.getByRole('button', { name: choice, exact: true }).click();
  await expect(dialog).toHaveCount(0);
  await evaluateWithRetainedPromise(page, async () => {
    await window.__labelOnApply;
    window.__GBDRAW_APP__.clickedFeature = null;
  });
};

// The popup note on why the feature has no label (requirement of Q1, Q2).
const popupLabelHint = (page, featureId) => evaluateWithRetainedPromise(page, async (id) => {
  const app = window.__GBDRAW_APP__;
  await app.openFeatureEditorFromList(app.extractedFeatures.find((item) => item.svg_id === id), null);
  await window.Vue.nextTick();
  const hint = document.querySelector('[data-label-visibility-hint]')?.textContent.trim() || '';
  app.clickedFeature = null;
  return hint;
}, featureId);

const hideFeature = (page, featureId) => evaluateWithRetainedPromise(page, async (id) => {
  const app = window.__GBDRAW_APP__;
  await app.openFeatureEditorFromList(app.extractedFeatures.find((item) => item.svg_id === id), null);
  app.clickedFeature.featureVisibility = 'off';
  await app.updateClickedFeatureVisibility('off');
  if (app.featureVisibilityScopeDialog.show) await app.handleFeatureVisibilityScopeChoice('feature');
  app.clickedFeature = null;
}, featureId);

// Auto Reflow is off, so a hidden feature stays in the Result until a choice
// below queues the label reflow, which draws it no more.
test('Label visibility On for a hidden feature shows the feature and label or keeps the feature hidden', async ({ page }) => {
  test.setTimeout(600_000);
  await openWithGenBank(page, FORCED_LABEL_FIXTURE, () => {
    window.__GBDRAW_APP__.form.labels_mode = 'none';
  });
  await generate(page);
  const { fl1, fl2 } = await featureIdsByLocator(page);
  const labelState = () => labelEditorState(page, 'labels.circular.scope');
  await hideFeature(page, fl1);
  const hidden = await labelOnFacts(page);
  expect(hidden.featureVisibility).toEqual({ [fl1]: 'off' });

  await startLabelOn(page, fl1, 'FL1_SHOWN');
  await answerLabelOn(page, 'Feature Is Hidden', 'Cancel');
  expect(await labelOnFacts(page)).toEqual(hidden);

  await startLabelOn(page, fl1, 'FL1_SHOWN');
  await answerLabelOn(page, 'Feature Is Hidden', 'Show feature and label');
  expect(await labelOnFacts(page)).toEqual({
    ...hidden,
    undo: hidden.undo + 1,
    redo: 0,
    labelVisibility: { [fl1]: 'on' },
    labelText: { [fl1]: 'FL1_SHOWN' },
    featureVisibility: { [fl1]: 'on' }
  });
  await waitForLabelReflow(page);
  expect(await labelState()).toMatchObject({ reflowError: null, labels: [[fl1, 'FL1_SHOWN']] });

  await hideFeature(page, fl2);
  expect(await popupLabelHint(page, fl2)).toBe('This feature has no label in the current Result. The feature is hidden.');
  const shown = await labelOnFacts(page);
  await startLabelOn(page, fl2, 'FL2_KEPT');
  await answerLabelOn(page, 'Feature Is Hidden', 'Keep feature hidden');
  expect(await labelOnFacts(page)).toEqual({
    ...shown,
    undo: shown.undo + 1,
    redo: 0,
    labelVisibility: { [fl1]: 'on', [fl2]: 'on' },
    labelText: { [fl1]: 'FL1_SHOWN', [fl2]: 'FL2_KEPT' },
    featureVisibility: { [fl1]: 'on', [fl2]: 'off' }
  });
  await waitForLabelReflow(page);
  expect(await labelState()).toMatchObject({ reflowError: null, labels: [[fl1, 'FL1_SHOWN']] });

  await generate(page);
  expect(await labelState()).toMatchObject({ error: null, reflowError: null, labels: [[fl1, 'FL1_SHOWN']] });
});

// The given features that the mounted Result draws.
const drawnFeatureIds = (page, featureIds) => page.evaluate((ids) => {
  const root = window.__GBDRAW_APP__.svgContainer.querySelector('svg');
  return ids.filter((id) => [...root.querySelectorAll(`[data-gbdraw-feature-id="${CSS.escape(id)}"]`)]
    .some((element) => !element.closest('[display="none"]')));
}, featureIds);

// R-2 (FINDINGS, remaining issue of #764; Owner Q2, R3, R4): a visibility rule
// that Generate applies hides the feature for the Label On dialog too. The
// rules here, as a loaded visibility TSV can give them, name a record ID with a
// locus_tag regex in another letter case, and the drawn hash in upper case.
// The dialog, the popup note, the live preview, and Generate answer "is this
// feature drawn?" with Python's rule.
const addVisibilityRule = (page, fields) => evaluateWithRetainedPromise(page, async (ruleFields) => {
  const app = window.__GBDRAW_APP__;
  await app.addFeatureVisibilityRule();
  const index = app.featureVisibilityManualRules.length - 1;
  for (const [field, value] of Object.entries(ruleFields)) await app.setFeatureVisibilityRuleField(index, field, value);
}, fields);

// The per-feature edits of one source feature, read by its identity.
const identityEdits = (page, identity) => page.evaluate(async (key) => {
  const { state } = await import('./js/state.js');
  const row = state.featureOverrides[key] || {};
  return { featureVisibility: row.featureVisibility ?? null, labelVisibility: row.labelVisibility ?? null };
}, identity);

test('Label visibility On for a feature that a visibility rule hides asks first, and live equals Generate', async ({ page }) => {
  test.setTimeout(600_000);
  await openWithGenBank(page, FORCED_LABEL_FIXTURE, () => {
    window.__GBDRAW_APP__.form.labels_mode = 'none';
  });
  await generate(page);
  const { fl1, fl2 } = await featureIdsByLocator(page);
  const fl2Identity = await page.evaluate(async (id) => {
    const { featureIdentityKeyOf } = await import('./js/services/feature-placement.js');
    return featureIdentityKeyOf(window.__GBDRAW_APP__.extractedFeatures.find((item) => item.svg_id === id));
  }, fl2);
  const labelState = () => labelEditorState(page, 'labels.circular.scope');
  const shown = async () => ({ labels: (await labelState()).labels, drawn: await drawnFeatureIds(page, [fl1, fl2]) });

  await addVisibilityRule(page, { recordId: 'FORCEDLBL', featureType: 'CDS', qualifier: 'locus_tag', value: '^fl1$', action: 'off' });
  const before = await labelOnFacts(page);
  await startLabelOn(page, fl1, 'FL1_SHOWN');
  await answerLabelOn(page, 'Feature Is Hidden', 'Show feature and label');
  expect(await labelOnFacts(page)).toEqual({
    ...before,
    undo: before.undo + 1,
    redo: 0,
    labelVisibility: { [fl1]: 'on' },
    labelText: { [fl1]: 'FL1_SHOWN' },
    featureVisibility: { [fl1]: 'on' }
  });
  await waitForLabelReflow(page);
  const live = await shown();
  expect(live).toEqual({ labels: [[fl1, 'FL1_SHOWN']], drawn: [fl1, fl2] });
  await generate(page);
  expect((await labelState()).reflowError).toBeNull();
  expect(await shown()).toEqual(live);

  await addVisibilityRule(page, { recordId: '*', featureType: '*', qualifier: 'hash', value: `^${fl2.toUpperCase()}$`, action: 'off' });
  expect(await popupLabelHint(page, fl2)).toBe('This feature has no label in the current Result. The feature is hidden.');
  await startLabelOn(page, fl2, 'FL2_KEPT');
  await answerLabelOn(page, 'Feature Is Hidden', 'Keep feature hidden');
  expect(await identityEdits(page, fl2Identity)).toEqual({ featureVisibility: null, labelVisibility: 'on' });
  await waitForLabelReflow(page);
  const kept = await shown();
  expect(kept).toEqual({ labels: [[fl1, 'FL1_SHOWN']], drawn: [fl1] });
  await generate(page);
  expect((await labelState()).reflowError).toBeNull();
  expect(await shown()).toEqual(kept);
});

// OV-19 (PD-OI-066, R1, R10, R11): a Feature Visibility rule edit in the
// Features panel is live like the other composition edits. Adding, editing,
// moving, and deleting a rule show on the Result at once as Generate draws it,
// and each is one Undo step that Undo and Redo project. A regular expression
// that Generate rejects stays in the draft and reports Generate's error; the
// Result stays as the failed Generate leaves it.
test('a Feature Visibility rule edit shows on the Result at once, as Generate draws it', async ({ page }) => {
  test.setTimeout(600_000);
  await openWithGenBank(page, FORCED_LABEL_FIXTURE);
  await generate(page);
  const { fl1, fl2 } = await featureIdsByLocator(page);
  const drawn = () => drawnFeatureIds(page, [fl1, fl2]);
  const undoCount = () => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
  const edit = (action, ...args) => evaluateWithRetainedPromise(page, ({ name, values }) => window.__GBDRAW_APP__[name](...values),
    { name: action, values: args });
  const reportedError = () => page.evaluate(async () => {
    const error = (await import('./js/state.js')).state.errorLog.value;
    return error && { code: error.code, row: error.context?.row };
  });
  expect(await drawn()).toEqual([fl1, fl2]);

  await addVisibilityRule(page, { recordId: 'FORCEDLBL', featureType: 'CDS', qualifier: 'locus_tag', value: '^fl1$', action: 'off' });
  expect(await drawn()).toEqual([fl2]);
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
  expect(await drawn()).toEqual([fl1, fl2]);
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.redo());
  expect(await drawn()).toEqual([fl2]);

  await edit('setFeatureVisibilityRuleField', 0, 'action', 'show');
  expect(await drawn()).toEqual([fl1, fl2]);
  await addVisibilityRule(page, { recordId: '*', featureType: 'CDS', qualifier: 'locus_tag', value: '^fl', action: 'off' });
  expect(await drawn()).toEqual([fl1]);
  const beforeMove = await undoCount();
  await edit('moveFeatureVisibilityRuleUp', 1);
  expect(await undoCount()).toBe(beforeMove + 1);
  expect(await drawn()).toEqual([]);

  const breakRegex = async () => {
    await edit('setFeatureVisibilityRuleField', 1, 'value', '^fl1(');
    expect(await page.evaluate(() => window.__GBDRAW_APP__.featureVisibilityManualRules[1].value)).toBe('^fl1(');
    expect(await drawn()).toEqual([]);
  };
  const fixRegex = () => edit('setFeatureVisibilityRuleField', 1, 'value', '^fl1$');
  await breakRegex();
  const liveError = await reportedError();
  expect(liveError).toEqual({ code: 'REGEX_SYNTAX', row: 2 });
  await fixRegex();
  expect(await reportedError()).toBeNull();
  expect(await drawn()).toEqual([]);
  await breakRegex();
  expect(await reportedError()).toEqual(liveError);
  await page.getByRole('button', { name: 'Generate Diagram', exact: true }).click();
  await expect.poll(() => page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    return { processing: state.processing.value, operation: state.errorLog.value?.operation };
  }), { timeout: 180_000 }).toEqual({ processing: false, operation: 'generate' });
  expect(await reportedError()).toEqual(liveError);
  expect(await drawn()).toEqual([]);
  await fixRegex();

  const beforeDelete = await undoCount();
  await edit('removeFeatureVisibilityRule', 0);
  expect(await undoCount()).toBe(beforeDelete + 1);
  const live = await drawn();
  expect(live).toEqual([fl1, fl2]);
  await generate(page);
  expect(await drawn()).toEqual(live);
});

// OV-19: the popup note reads Python's matches of the edited rules, also when
// the edit leaves the Result as it was.
test('the feature popup note follows a Feature Visibility rule edit', async ({ page }) => {
  test.setTimeout(600_000);
  await openWithGenBank(page, FORCED_LABEL_FIXTURE, () => {
    window.__GBDRAW_APP__.form.labels_mode = 'none';
  });
  await generate(page);
  const { fl1, fl2 } = await featureIdsByLocator(page);
  await addVisibilityRule(page, { recordId: '*', featureType: 'CDS', qualifier: 'locus_tag', value: '^fl1$', action: 'off' });
  await addVisibilityRule(page, { recordId: '*', featureType: 'CDS', qualifier: 'locus_tag', value: '^fl2$', action: 'off' });
  expect(await drawnFeatureIds(page, [fl1, fl2])).toEqual([]);
  await evaluateWithRetainedPromise(page, async (id) => {
    const app = window.__GBDRAW_APP__;
    await app.openFeatureEditorFromList(app.extractedFeatures.find((item) => item.svg_id === id), null);
  }, fl2);
  const note = page.locator('[data-label-visibility-hint]');
  const hidden = 'This feature has no label in the current Result. The feature is hidden.';
  await expect(note).toHaveText(hidden);
  // The first rule now hides FL2 too; Python has not matched its new regex yet.
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.setFeatureVisibilityRuleField(0, 'value', '^fl[12]$'));
  expect(await drawnFeatureIds(page, [fl1, fl2])).toEqual([]);
  await expect(note).toHaveText(hidden);
});

// R3 (found with OV-19): Load Feature Edits TSV projects visibility as a History
// apply does, so a feature whose edit it clears follows the rules at once.
test('Load Feature Edits TSV applies the visibility rules to a feature whose edit it clears', async ({ page }) => {
  test.setTimeout(600_000);
  await openWithGenBank(page, FORCED_LABEL_FIXTURE);
  await generate(page);
  const { fl1, fl2 } = await featureIdsByLocator(page);
  await addVisibilityRule(page, { recordId: '*', featureType: 'CDS', qualifier: 'locus_tag', value: '^fl1$', action: 'off' });
  await evaluateWithRetainedPromise(page, async (id) => {
    const app = window.__GBDRAW_APP__;
    await app.openFeatureEditorFromList(app.extractedFeatures.find((item) => item.svg_id === id), null);
    app.clickedFeature.featureVisibility = 'on';
    await app.updateClickedFeatureVisibility('on');
    await app.handleFeatureVisibilityScopeChoice('feature');
    app.clickedFeature = null;
  }, fl1);
  await generate(page);
  expect(await drawnFeatureIds(page, [fl1, fl2])).toEqual([fl1, fl2]);
  await evaluateWithRetainedPromise(page, async () => {
    const text = 'record\tfeature_selector\tfeature_visibility\tlabel_visibility\tlabel_text\n';
    const files = [new File([text], 'edits.tsv', { type: 'text/plain' })];
    await window.__GBDRAW_APP__.loadFeatureEditTable({ target: { files, value: '' } });
  });
  expect(await drawnFeatureIds(page, [fl1, fl2])).toEqual([fl2]);
  await generate(page);
  expect(await drawnFeatureIds(page, [fl1, fl2])).toEqual([fl2]);
});

// The given features whose label the mounted Result shows.
const labeledFeatureIds = (page, featureIds) => page.evaluate((ids) => {
  const root = window.__GBDRAW_APP__.svgContainer.querySelector('svg');
  return ids.filter((id) => [...root.querySelectorAll(`text[data-label-feature-id="${CSS.escape(id)}"]`)]
    .some((element) => !element.closest('[display="none"]')));
}, featureIds);
const autoReflowOff = (page) => page.evaluate(() => window.__GBDRAW_APP__.autoLabelReflowEnabled === false);

// OV-35 (Owner 2026-10-05): a feature that a Feature Visibility rule hides
// takes its label with it, as Generate draws it, also with Auto Reflow off and
// before the editor has bound the labels. The rule edit, Undo and Redo (History
// apply), and Generate show the same features and labels.
test('a Feature Visibility rule that hides a feature hides its label too', async ({ page }) => {
  test.setTimeout(600_000);
  await openWithGenBank(page, FORCED_LABEL_FIXTURE, () => {
    window.__GBDRAW_APP__.form.labels_mode = 'out';
  });
  await generate(page);
  expect(await autoReflowOff(page)).toBe(true);
  const { fl1, fl2 } = await featureIdsByLocator(page);
  const shown = async () => ({
    drawn: await drawnFeatureIds(page, [fl1, fl2]),
    labeled: await labeledFeatureIds(page, [fl1, fl2])
  });
  const all = { drawn: [fl1, fl2], labeled: [fl1, fl2] };
  const withoutFl1 = { drawn: [fl2], labeled: [fl2] };
  expect(await shown()).toEqual(all);

  await addVisibilityRule(page, { recordId: 'FORCEDLBL', featureType: 'CDS', qualifier: 'locus_tag', value: '^fl1$', action: 'off' });
  expect(await shown()).toEqual(withoutFl1);
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
  expect(await shown()).toEqual(all);
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.redo());
  expect(await shown()).toEqual(withoutFl1);
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.setFeatureVisibilityRuleField(0, 'action', 'show'));
  expect(await shown()).toEqual(all);
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.setFeatureVisibilityRuleField(0, 'action', 'off'));
  expect(await shown()).toEqual(withoutFl1);
  await generate(page);
  expect(await shown()).toEqual(withoutFl1);
});

// OV-34: with Auto Reflow off, deleting the rule that hid a feature at Generate
// draws the feature and its label live (the automatic rerender), as Generate
// draws them.
test('deleting a Feature Visibility rule draws the feature it hid at Generate', async ({ page }) => {
  test.setTimeout(600_000);
  await openWithGenBank(page, FORCED_LABEL_FIXTURE, () => {
    window.__GBDRAW_APP__.form.labels_mode = 'out';
  });
  await generate(page);
  expect(await autoReflowOff(page)).toBe(true);
  const { fl1, fl2 } = await featureIdsByLocator(page);
  await addVisibilityRule(page, { recordId: 'FORCEDLBL', featureType: 'CDS', qualifier: 'locus_tag', value: '^fl1$', action: 'off' });
  await generate(page);
  const shown = async () => ({
    drawn: await drawnFeatureIds(page, [fl1, fl2]),
    labeled: await labeledFeatureIds(page, [fl1, fl2])
  });
  expect(await shown()).toEqual({ drawn: [fl2], labeled: [fl2] });
  expect(await page.evaluate((id) => Boolean(window.__GBDRAW_APP__.svgContainer
    .querySelector(`[data-gbdraw-feature-id="${CSS.escape(id)}"]`)), fl1)).toBe(false);

  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.removeFeatureVisibilityRule(0));
  await expect.poll(shown, { timeout: 300_000 }).toEqual({ drawn: [fl1, fl2], labeled: [fl1, fl2] });
  await waitForLabelReflow(page);
  expect((await labelEditorState(page, 'labels.circular.scope')).reflowError).toBeNull();
  await generate(page);
  expect(await shown()).toEqual({ drawn: [fl1, fl2], labeled: [fl1, fl2] });
});

// OV-35 on Result display (R3): a rule edited while one batch Result is shown
// hides the feature and its label on the other Result when it is displayed.
test('a Feature Visibility rule hides the label on the other batch Result when it is displayed', async ({ page }) => {
  test.setTimeout(600_000);
  await openBatch(page);
  const ids = await page.evaluate(() => Object.fromEntries(window.__GBDRAW_APP__.extractedFeatures
    .filter((feature) => /_000[14]$/.test(feature.locus_tag || ''))
    .map((feature) => [feature.locus_tag, feature.svg_id])));
  const shownOn = async (tags) => {
    const featureIds = tags.map((tag) => ids[tag]);
    return { drawn: await drawnFeatureIds(page, featureIds), labeled: await labeledFeatureIds(page, featureIds) };
  };
  expect(await shownOn(['TESTA_0001', 'TESTA_0004'])).toEqual({
    drawn: [ids.TESTA_0001, ids.TESTA_0004], labeled: [ids.TESTA_0001, ids.TESTA_0004]
  });
  await addVisibilityRule(page, { recordId: '*', featureType: 'CDS', qualifier: 'locus_tag', value: '_0004$', action: 'off' });
  expect(await shownOn(['TESTA_0001', 'TESTA_0004'])).toEqual({ drawn: [ids.TESTA_0001], labeled: [ids.TESTA_0001] });
  await selectResult(page, 1);
  expect(await shownOn(['TESTB_0001', 'TESTB_0004'])).toEqual({ drawn: [ids.TESTB_0001], labeled: [ids.TESTB_0001] });
});

test('Label visibility On for an underlay feature is kept without a label until the rendering changes', async ({ page }) => {
  test.setTimeout(600_000);
  await openWithGenBank(page, FORCED_LABEL_FIXTURE, () => {
    window.__GBDRAW_APP__.form.labels_mode = 'none';
  });
  await generate(page);
  const { repeat } = await featureIdsByLocator(page);
  const labelState = () => labelEditorState(page, 'labels.circular.scope');
  const before = await labelOnFacts(page);

  await startLabelOn(page, repeat, 'RPT_FORCED');
  await answerLabelOn(page, 'Label Cannot Be Drawn', 'Cancel');
  expect(await labelOnFacts(page)).toEqual(before);

  await startLabelOn(page, repeat, 'RPT_FORCED');
  await answerLabelOn(page, 'Label Cannot Be Drawn', 'Keep without label');
  expect(await labelOnFacts(page)).toEqual({
    ...before,
    undo: before.undo + 1,
    redo: 0,
    labelVisibility: { [repeat]: 'on' },
    labelText: { [repeat]: 'RPT_FORCED' }
  });
  await waitForLabelReflow(page);
  expect(await labelState()).toMatchObject({ reflowError: null, labels: [] });
  expect(await popupLabelHint(page, repeat)).toBe(
    'This feature has no label in the current Result. Labels are not drawn for features drawn as "Underlay". '
    + 'Its Label visibility "On" applies when the label can be drawn.'
  );

  await generate(page);
  expect(await labelState()).toMatchObject({ error: null, labels: [] });
  await page.evaluate(() => window.__GBDRAW_APP__.setFeatureShape('repeat_region', 'rectangle'));
  await generate(page);
  expect(await labelState()).toMatchObject({ error: null, labels: [[repeat, 'RPT_FORCED']] });
});

// Generate keeps Label Rendering = Embedded Only only while Show Labels draws labels.
test('Label visibility On that does not fit with Embedded Only asks after the reflow', async ({ page }) => {
  test.setTimeout(600_000);
  await openWithGenBank(page, FORCED_LABEL_FIXTURE, () => {
    const app = window.__GBDRAW_APP__;
    app.form.labels_mode = 'out';
    app.adv.label_rendering = 'embedded_only';
  });
  await generate(page);
  const { fl1 } = await featureIdsByLocator(page);
  const labelState = () => labelEditorState(page, 'labels.circular.scope');
  const committedRendering = () => page.evaluate(async () => (await import('./js/services/config.js'))
    .getCommittedCanonicalRenderRequest()?.diagramOptions?.configOverrides?.['labels.rendering']);
  expect(await committedRendering()).toBe('embedded_only');
  const initialLabels = (await labelState()).labels;
  expect(initialLabels.some(([featureId]) => featureId === fl1)).toBe(true);
  const longText = 'A_LABEL_TEXT_FAR_TOO_LONG_TO_FIT_INSIDE_ITS_FEATURE_'.repeat(4);
  const before = await labelOnFacts(page);

  // Cancel restores the label intent and redraws the label it showed before.
  await startLabelOn(page, fl1, longText);
  await answerLabelOn(page, 'Label Does Not Fit', 'Cancel');
  expect(await labelOnFacts(page)).toEqual(before);
  await waitForLabelReflow(page);
  expect(await labelState()).toMatchObject({ reflowError: null, labels: initialLabels });

  await startLabelOn(page, fl1, longText);
  await answerLabelOn(page, 'Label Does Not Fit', 'Keep without label');
  expect(await labelOnFacts(page)).toEqual({
    ...before,
    undo: before.undo + 1,
    redo: 0,
    labelVisibility: { [fl1]: 'on' },
    labelText: { [fl1]: longText }
  });
  await waitForLabelReflow(page);
  const kept = initialLabels.filter(([featureId]) => featureId !== fl1);
  expect(await labelState()).toMatchObject({ reflowError: null, labels: kept });
  expect(await popupLabelHint(page, fl1)).toBe(
    'This feature has no label in the current Result. With "Label Rendering" = "Embedded Only", '
    + 'a label is drawn only when it fits inside its feature. '
    + 'Its Label visibility "On" applies when the label can be drawn.'
  );

  await generate(page);
  expect(await labelState()).toMatchObject({ error: null, labels: kept });
  expect(await committedRendering()).toBe('embedded_only');
});

// R-4 (FINDINGS, residual of #764): the Label Not Shown dialog opens when a text
// edit leaves a feature without a label, and its first sentence names the reason
// the popup note names. The two reasons below differ: a feature drawn as
// Underlay, and Show Labels None.
const editTextOnly = (page, locator, text) => evaluateWithRetainedPromise(page, async (edit) => {
  const app = window.__GBDRAW_APP__;
  const feature = app.extractedFeatures.find((item) => item.locus_tag === edit.locator || item.type === edit.locator);
  await app.openFeatureEditorFromList(feature, null);
  await window.Vue.nextTick();
  const note = document.querySelector('[data-label-visibility-hint]')?.textContent.trim() || '';
  app.clickedFeature.labelText = edit.text;
  await app.updateClickedFeatureLabelText();
  await window.Vue.nextTick();
  const heading = [...document.querySelectorAll('h3')].find((node) => node.textContent === 'Label Not Shown');
  const dialog = heading?.parentElement.querySelector('p').textContent.trim() || null;
  if (app.hiddenLabelTextDialog.show) await app.handleHiddenLabelTextChoice('text_only');
  app.clickedFeature = null;
  return { note, dialog };
}, { locator, text });

for (const [name, mode, locator, sentence] of [
  ['a feature drawn as Underlay', 'out', 'repeat_region', 'Labels are not drawn for features drawn as "Underlay".'],
  ['Show Labels None', 'none', 'FL1', 'Labels are set to "None"']
]) {
  test(`the Label Not Shown dialog names the reason: ${name}`, async ({ page }) => {
    test.setTimeout(300_000);
    await openWithGenBank(page, FORCED_LABEL_FIXTURE, () => {});
    await page.evaluate((labelsMode) => { window.__GBDRAW_APP__.form.labels_mode = labelsMode; }, mode);
    await generate(page);
    const { note, dialog } = await editTextOnly(page, locator, 'EDITED_TEXT');
    expect(note).toContain(sentence);
    expect(dialog).toContain(`This feature has no label in the current Result. ${sentence}`);
    expect(dialog).not.toContain('With the current Show Labels and label filter settings');
  });
}
