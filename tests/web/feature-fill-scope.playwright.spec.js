const { test, expect } = require('@playwright/test');
const { evaluateWithRetainedPromise } = require('./helpers/app-lifecycle.cjs');
const { seeds, load, popup: openFeaturePopup } = require('./helpers/mode-transition.cjs');

const loadGallerySession = async (page, filename) => evaluateWithRetainedPromise(page, async (name) => {
  const response = await fetch(`/gbdraw/web/gallery/sessions/${name}`);
  if (!response.ok) throw new Error(`Could not load ${name}: ${response.status}`);
  const bytes = await response.arrayBuffer();
  const file = new File([bytes], name, {
    type: name.endsWith('.gz') ? 'application/gzip' : 'application/json'
  });
  return window.__GBDRAW_APP__.importSession({
    target: { files: [file], value: 'selected' }
  });
}, filename);

const inspectFillState = async (page, target) => page.evaluate(async ({ caption, featureIds }) => {
  const app = window.__GBDRAW_APP__;
  const svg = document.querySelector('.origin-top svg');
  const { getFeatureFillElements } = await import(
    '/gbdraw/web/js/app/feature-editor/svg-actions.js'
  );
  const legendGroup = Array.from(svg.querySelectorAll('g[data-legend-key]'))
    .find((group) => group.getAttribute('data-legend-key') === caption);
  const legendEntry = app.legendEntries.find((entry) => entry.caption === caption);
  return {
    featureFills: Object.fromEntries(featureIds.map((featureId) => [
      featureId,
      getFeatureFillElements(svg, featureId).map((element) => element.getAttribute('fill'))
    ])),
    legendFill: legendGroup
      ?.querySelector('path[fill]:not([fill="none"])')
      ?.getAttribute('fill') || null,
    // OV-282 (D-26): the color the Legend row shows, read through its one reader.
    legendEntryColor: app.legendEntryColor(legendEntry) || null,
    featureColorOverrides: JSON.stringify(app.featureColorOverrides),
    manualSpecificRules: JSON.stringify(app.manualSpecificRules),
    resultContent: app.results[app.selectedResultIndex]?.content || '',
    undoCount: window.__GBDRAW_HISTORY__.getUndoCount(),
    redoCount: window.__GBDRAW_HISTORY__.getRedoCount()
  };
}, target);

const inspectStrokeState = async (page, target) => page.evaluate(async ({ caption, featureIds }) => {
  const app = window.__GBDRAW_APP__;
  const svg = document.querySelector('.origin-top svg');
  const { getFeatureElements } = await import(
    '/gbdraw/web/js/app/feature-editor/svg-actions.js'
  );
  const legendGroup = Array.from(svg.querySelectorAll('g[data-legend-key]'))
    .find((group) => group.getAttribute('data-legend-key') === caption);
  const legendSwatch = legendGroup?.querySelector('path[fill]:not([fill="none"])') || null;
  return {
    featureStrokes: Object.fromEntries(featureIds.map((featureId) => [
      featureId,
      getFeatureElements(svg, featureId).map((element) => {
        const width = element.getAttribute('stroke-width');
        return {
          color: element.getAttribute('stroke'),
          width: width === null ? null : Number(width)
        };
      })
    ])),
    legendStroke: legendSwatch ? {
      color: legendSwatch.getAttribute('stroke'),
      width: legendSwatch.getAttribute('stroke-width') === null
        ? null
        : Number(legendSwatch.getAttribute('stroke-width'))
    } : null,
    featureStrokeOverrides: JSON.stringify(app.featureStrokeOverrides),
    resultContent: app.results[app.selectedResultIndex]?.content || '',
    undoCount: window.__GBDRAW_HISTORY__.getUndoCount(),
    redoCount: window.__GBDRAW_HISTORY__.getRedoCount()
  };
}, target);

const changeColor = async (picker, color) => {
  await picker.evaluate((element, value) => {
    element.value = value;
    element.dispatchEvent(new Event('input', { bubbles: true }));
    element.dispatchEvent(new Event('change', { bubbles: true }));
  }, color);
  // A cold Session prepares Python matches before the synchronous edit phase.
  await picker.page().waitForFunction(() => !window.__GBDRAW_APP__.ruleMatchingPending);
};

test.beforeEach(async ({ page }) => {
  page.on('dialog', (dialog) => dialog.accept());
  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
});

test('inherited Feature fill uses the existing scope dialog with atomic Undo and Redo', async ({
  page
}) => {
  const imported = await loadGallerySession(
    page,
    'HmmtDNA_basic_circular.gbdraw-session.json'
  );
  expect(imported.status).toBe('ok');

  const drawer = page.locator('.right-drawer');
  await page.locator('.drawer-toggle').click();
  await expect(drawer.getByText(/^Features \(\d+\)$/)).toBeVisible();

  const target = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const svg = document.querySelector('.origin-top svg');
    const { formatFeatureLocation } = await import('/gbdraw/web/js/services/feature-utils.js');
    for (const feature of app.visibleFeatureRows) {
      const featureId = String(feature?.svg_id || '').trim();
      if (!featureId || app.getFeatureColorValue(feature) !== null) continue;
      if (!svg.querySelector(`[data-gbdraw-feature-id="${CSS.escape(featureId)}"]`)) continue;
      const legendEntry = app.legendEntries.find((entry) => entry.caption === feature.type);
      if (!legendEntry) continue;
      const featureIds = app.extractedFeatures
        .filter((candidate) => (
          candidate.type === feature.type &&
          app.getFeatureColorValue(candidate) === null &&
          svg.querySelector(
            `[data-gbdraw-feature-id="${CSS.escape(String(candidate.svg_id || ''))}"]`
          )
        ))
        .map((candidate) => candidate.svg_id);
      if (featureIds.length < 2) continue;
      return {
        id: featureId,
        location: formatFeatureLocation(feature),
        color: app.getFeatureColor(feature),
        caption: legendEntry.caption,
        featureIds
      };
    }
    throw new Error('No visible inherited Feature with a shared legend item was found.');
  });

  const featureRow = drawer.locator(`span[title="${target.location}"]`).locator('..');
  const listPicker = featureRow.locator('input[type="color"]');
  await expect(listPicker).toBeEnabled();
  await expect(listPicker).toHaveValue(target.color.toLowerCase());

  await featureRow.getByRole('button', { name: 'Edit' }).click();

  const popup = page.locator('.feature-popup');
  await expect(popup).toBeVisible();
  let popupPicker = popup.locator('input[type="color"][aria-label="Feature fill color"]').first();
  await expect(popupPicker).toBeEnabled();
  await expect(popupPicker).toHaveValue(target.color.toLowerCase());

  const before = await inspectFillState(page, target);
  expect(before.undoCount).toBe(0);
  expect(before.redoCount).toBe(0);

  await changeColor(popupPicker, '#7357c8');
  const scopeDialog = page.getByRole('heading', { name: 'Color Change Scope' }).locator('..');
  await expect(scopeDialog).toBeVisible();
  await expect.poll(
    () => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())
  ).toBe(0);
  await scopeDialog.getByRole('button', { name: 'Cancel' }).click();
  await expect(scopeDialog).toBeHidden();
  expect(await inspectFillState(page, target)).toEqual(before);

  await featureRow.getByRole('button', { name: 'Edit' }).click();
  await expect(popup).toBeVisible();
  popupPicker = popup.locator('input[type="color"][aria-label="Feature fill color"]').first();
  await expect(popupPicker).toBeEnabled();
  await expect(popupPicker).toHaveValue(target.color.toLowerCase());

  const changedColor = '#2f7fc1';
  await changeColor(popupPicker, changedColor);
  await expect(scopeDialog).toBeVisible();
  await expect.poll(
    () => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())
  ).toBe(0);

  const captionScopeButton = scopeDialog.locator('button.bg-green-600');
  await expect(captionScopeButton).toBeVisible();
  await captionScopeButton.click();
  // The choice prepares the rules it adds; on a Session just loaded that
  // starts the Worker runtime (OV-225), so wait for the choice to finish.
  await page.waitForFunction(() => !window.__GBDRAW_APP__.featureStyleScopeDialog.show);
  await expect(scopeDialog).toBeHidden();

  await expect.poll(
    () => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())
  ).toBe(1);
  const changed = await inspectFillState(page, target);
  expect(changed.redoCount).toBe(0);
  Object.values(changed.featureFills).forEach((fills) => {
    expect(fills.length).toBeGreaterThan(0);
    expect(fills.every((fill) => fill === changedColor)).toBe(true);
  });
  expect(changed.legendFill).toBe(changedColor);
  expect(changed.legendEntryColor).toBe(changedColor);

  await page.getByRole('button', { name: 'Undo', exact: true }).click();
  await expect.poll(
    () => page.evaluate(() => window.__GBDRAW_HISTORY__.getRedoCount())
  ).toBe(1);
  const undone = await inspectFillState(page, target);
  expect(undone.featureFills).toEqual(before.featureFills);
  expect(undone.legendEntryColor).toBe(before.legendEntryColor);
  expect(undone.legendFill).toBe(before.legendFill);

  await page.getByRole('button', { name: 'Redo', exact: true }).click();
  await expect.poll(
    () => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())
  ).toBe(1);
  const redone = await inspectFillState(page, target);
  expect(redone.featureFills).toEqual(changed.featureFills);
  expect(redone.legendFill).toBe(changedColor);
  expect(redone.legendEntryColor).toBe(changedColor);
});

test('Feature stroke width steppers defer scope selection and stroke changes keep atomic Undo and Redo', async ({
  page
}) => {
  const imported = await loadGallerySession(
    page,
    'HmmtDNA_basic_circular.gbdraw-session.json'
  );
  expect(imported.status).toBe('ok');

  const drawer = page.locator('.right-drawer');
  await page.locator('.drawer-toggle').click();
  await expect(drawer.getByText(/^Features \(\d+\)$/)).toBeVisible();

  const target = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const svg = document.querySelector('.origin-top svg');
    const { getFeatureElements } = await import(
      '/gbdraw/web/js/app/feature-editor/svg-actions.js'
    );
    const { formatFeatureLocation } = await import('/gbdraw/web/js/services/feature-utils.js');
    for (const feature of app.visibleFeatureRows) {
      const featureId = String(feature?.svg_id || '').trim();
      if (!featureId || getFeatureElements(svg, featureId).length === 0) continue;
      const entry = app.legendEntries.find((candidate) => candidate.caption === feature.type);
      if (!entry) continue;
      const featureIds = app.extractedFeatures
        .filter((candidate) => (
          candidate.type === feature.type &&
          getFeatureElements(svg, String(candidate.svg_id || '')).length > 0
        ))
        .map((candidate) => candidate.svg_id);
      if (featureIds.length < 2) continue;
      return {
        id: featureId,
        location: formatFeatureLocation(feature),
        caption: entry.caption,
        featureIds
      };
    }
    throw new Error('No visible Feature with a shared legend item was found.');
  });

  const featureRow = drawer.locator(`span[title="${target.location}"]`).locator('..');
  await featureRow.getByRole('button', { name: 'Edit' }).click();

  const popup = page.locator('.feature-popup');
  await expect(popup).toBeVisible();
  const strokeMode = popup.getByLabel('Feature stroke color mode', { exact: true });
  if (await strokeMode.inputValue() === 'auto') {
    await strokeMode.selectOption('color');
    await page.waitForFunction(() => !window.__GBDRAW_APP__.ruleMatchingPending);
  }
  let strokePicker = popup.getByLabel('Feature stroke color', { exact: true });
  await expect(strokePicker).toBeEnabled();
  await page.evaluate(() => window.__GBDRAW_HISTORY__.captureBaseline('Stroke color mode'));
  const before = await inspectStrokeState(page, target);
  expect(before.undoCount).toBe(0);
  expect(before.redoCount).toBe(0);

  const scopeDialog = page.getByRole('heading', { name: 'Stroke Change Scope' }).locator('..');
  const strokeWidthInput = popup.getByLabel('Feature stroke width', { exact: true });
  const initialStrokeWidth = Number(await strokeWidthInput.inputValue());
  await strokeWidthInput.press('ArrowUp');
  await expect.poll(async () => Number(await strokeWidthInput.inputValue())).toBeCloseTo(
    initialStrokeWidth + 0.1,
    5
  );
  await expect(scopeDialog).toBeHidden();
  expect(await inspectStrokeState(page, target)).toEqual(before);

  await strokeWidthInput.press('Enter');
  await expect(scopeDialog).toBeVisible();
  expect(await page.evaluate(() => window.__GBDRAW_APP__.featureStyleScopeDialog.strokeWidth))
    .toBeCloseTo(initialStrokeWidth + 0.1, 5);
  await expect.poll(
    () => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())
  ).toBe(0);
  await scopeDialog.getByRole('button', { name: 'Cancel' }).click();
  await expect(scopeDialog).toBeHidden();
  expect(await inspectStrokeState(page, target)).toEqual(before);

  await featureRow.getByRole('button', { name: 'Edit' }).click();
  await expect(popup).toBeVisible();
  strokePicker = popup.getByLabel('Feature stroke color', { exact: true });

  await changeColor(strokePicker, '#8a5a44');
  await expect(scopeDialog).toBeVisible();
  await expect.poll(
    () => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())
  ).toBe(0);
  await scopeDialog.getByRole('button', { name: 'Cancel' }).click();
  await expect(scopeDialog).toBeHidden();
  expect(await inspectStrokeState(page, target)).toEqual(before);

  await featureRow.getByRole('button', { name: 'Edit' }).click();
  await expect(popup).toBeVisible();
  strokePicker = popup.getByLabel('Feature stroke color', { exact: true });
  const changedColor = '#365f8d';
  await changeColor(strokePicker, changedColor);
  await expect(scopeDialog).toBeVisible();
  await expect.poll(
    () => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())
  ).toBe(0);

  const captionScopeButton = scopeDialog.locator('button.bg-green-600');
  await expect(captionScopeButton).toBeVisible();
  await captionScopeButton.click();
  await expect(scopeDialog).toBeHidden();

  await expect.poll(
    () => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())
  ).toBe(1);
  const changed = await inspectStrokeState(page, target);
  expect(changed.redoCount).toBe(0);
  Object.values(changed.featureStrokes).forEach((strokes) => {
    expect(strokes.length).toBeGreaterThan(0);
    expect(strokes.every((stroke) => stroke.color === changedColor)).toBe(true);
  });
  expect(changed.legendStroke?.color).toBe(changedColor);
  expect(changed.resultContent).toContain(`stroke="${changedColor}"`);

  await page.getByRole('button', { name: 'Undo', exact: true }).click();
  await expect.poll(
    () => page.evaluate(() => window.__GBDRAW_HISTORY__.getRedoCount())
  ).toBe(1);
  const undone = await inspectStrokeState(page, target);
  expect(undone.featureStrokes).toEqual(before.featureStrokes);
  expect(undone.legendStroke).toEqual(before.legendStroke);
  expect(undone.featureStrokeOverrides).toBe(before.featureStrokeOverrides);

  await page.getByRole('button', { name: 'Redo', exact: true }).click();
  await expect.poll(
    () => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())
  ).toBe(1);
  const redone = await inspectStrokeState(page, target);
  expect(redone.featureStrokes).toEqual(changed.featureStrokes);
  expect(redone.legendStroke).toEqual(changed.legendStroke);
  expect(redone.featureStrokeOverrides).toBe(changed.featureStrokeOverrides);
});

// OV-161 (gui-fix GX-16): right after a Session load, Cancel on the Color
// Change Scope dialog (Escape or the Cancel button) closes it at once and
// records no History step. Cancel used to run inside a History step, so it
// first waited for that step's intent capture, which took seconds on a
// freshly loaded Session. The dialog is opened from its state, as the
// accessibility spec does, so no earlier History step has captured the
// intent yet.
test('Cancel closes the Color Change Scope dialog at once after a Session load and records no step (OV-161)', async ({
  browser
}) => {
  test.setTimeout(240_000);
  const page = await load(browser, seeds.circular);
  try {
    await openFeaturePopup(page);
    const undoCount = await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
    const dialog = page.getByRole('dialog', { name: 'Color Change Scope' });
    for (const answer of ['Escape', 'Cancel']) {
      await page.evaluate(() => Object.assign(window.__GBDRAW_APP__.featureStyleScopeDialog, { kind: 'fill', show: true }));
      await expect(dialog).toBeVisible();
      // The dialog takes focus (UI-02), so Escape reaches it.
      await expect.poll(() => dialog.evaluate((element) => element.contains(document.activeElement))).toBe(true);
      const started = Date.now();
      if (answer === 'Escape') await page.keyboard.press('Escape');
      else await dialog.getByRole('button', { name: 'Cancel', exact: true }).click();
      await expect(dialog).toHaveCount(0, { timeout: 60_000 });
      expect(Date.now() - started, `milliseconds from ${answer} to a closed dialog`).toBeLessThan(1000);
    }
    expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(undoCount);
  } finally {
    await page.context().close();
  }
});

// PD-OI-088 (OIC-028): a popup fill pick opens the scope dialog ready, with
// focus on its first choice and no busy status. A choice keeps focus inside
// the dialog while it applies (seconds on a Session just loaded, as the choice
// starts the Worker runtime), so Escape reaches the dialog, which ignores it,
// and does not clear the feature selection behind it.
test('a popup scope dialog opens ready and keeps focus while its choice applies', async ({ browser }) => {
  test.setTimeout(180_000);
  const page = await load(browser, 'gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json');
  try {
    await page.evaluate(() => {
      const app = window.__GBDRAW_APP__;
      app.openFeatureEditorFromList(app.extractedFeatures.find((feature) => feature.type === 'tRNA'), { clientX: 220, clientY: 220 });
    });
    const undoCount = await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
    await page.getByRole('dialog', { name: /Feature details:/ }).getByLabel('Feature fill color').first().fill('#c026d3');
    const scope = page.getByRole('dialog', { name: 'Color Change Scope' });
    await expect(scope).toBeVisible();
    await page.waitForFunction(() => !window.__GBDRAW_HISTORY__.mutationPending());
    const focusIn = () => scope.evaluate((element) => ({
      inDialog: element.contains(document.activeElement),
      firstChoice: document.activeElement === element.querySelector('button:not(:disabled)')
    }));
    const atOpen = { focus: await focusIn(), status: await scope.getByRole('status').count() };
    const selected = await page.evaluate(() => {
      const app = window.__GBDRAW_APP__;
      app.selectedFeatureIds = new Set(app.extractedFeatures.filter((feature) => feature.type === 'CDS').map((feature) => feature.svg_id));
      return app.selectedFeatureIds.size;
    });
    expect(selected).toBeGreaterThan(0);
    await scope.getByRole('button').filter({ hasText: 'Apply to all "tRNA"' }).last().click();
    await expect(scope.getByRole('status')).toHaveText('Applying an edit…');
    const whileBusy = await focusIn();
    // W1: no control is enabled, so Tab and Shift+Tab keep focus on the panel.
    await page.keyboard.press('Tab');
    await page.keyboard.press('Shift+Tab');
    const afterTab = await scope.evaluate((element) => document.activeElement === element);
    await page.keyboard.press('Escape');
    const afterEscape = await page.evaluate(() => ({
      selected: window.__GBDRAW_APP__.selectedFeatureIds.size,
      dialogOpen: window.__GBDRAW_APP__.featureStyleScopeDialog.show
    }));
    expect({ atOpen, whileBusy, afterTab, afterEscape }).toEqual({
      atOpen: { focus: { inDialog: true, firstChoice: true }, status: 0 },
      whileBusy: { inDialog: true, firstChoice: false },
      afterTab: true,
      afterEscape: { selected, dialogOpen: true }
    });
    await expect(scope).toHaveCount(0, { timeout: 60_000 });
    expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(undoCount + 1);
  } finally {
    await page.context().close();
  }
});
