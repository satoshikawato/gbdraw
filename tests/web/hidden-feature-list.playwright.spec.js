// R-5 of the override-precedence audit (Owner decision 2026-10-05): a hidden
// feature stays in the Features list and in Search features, so it can be
// shown again. The row's Visibility checkbox shows whether the feature is
// drawn and sets its Feature visibility; the popup of a hidden feature says
// so. Live equals Generate (R1, R3), and a toggle is one History step (R11).
const { test, expect } = require('@playwright/test');
const { openWithGenBank, settle } = require('./helpers/audit-browser.cjs');
const { generate } = require('./helpers/mode-transition.cjs');

test.describe.configure({ retries: 0 });

const FIXTURE = 'tests/fixtures/forced_label_underlay.gb';

const openDrawer = async (page) => {
  if (!await page.locator('.right-drawer').isVisible()) await page.locator('.drawer-toggle').click();
  await expect(page.locator('.right-drawer')).toBeVisible();
};

const visibilityBox = (page, name) => page.locator('.right-drawer')
  .getByRole('checkbox', { name: `Visibility of ${name}`, exact: true });

// The features the mounted Result draws, by locus_tag, product, or note.
const drawnFeatures = (page) => page.evaluate(() => {
  const app = window.__GBDRAW_APP__;
  const root = app.svgContainer.querySelector('svg');
  const drawn = new Set([...root.querySelectorAll('[data-gbdraw-feature-id]')]
    .filter((element) => !element.closest('[display="none"]'))
    .map((element) => element.getAttribute('data-gbdraw-feature-id')));
  return app.extractedFeatures.filter((feature) => drawn.has(feature.svg_id))
    .map((feature) => feature.locus_tag || feature.product || feature.note).sort();
});

// The Features list as the drawer shows it: name and checkbox state.
const listedRows = (page) => page.evaluate(() => {
  const app = window.__GBDRAW_APP__;
  return app.filteredFeatures.map((feature) => [
    feature.locus_tag || feature.product || feature.note,
    app.featureListState(feature).drawn
  ]);
});

const undoCount = (page) => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());

// The forced rerender after a toggle that needs geometry.
const waitForRerender = (page) => expect.poll(() => page.evaluate(async () => {
  const { state } = await import('./js/state.js');
  return state.labelReflowProcessing.value || state.processing.value;
}), { timeout: 300_000 }).toBe(false);

const labels = (page) => page.evaluate(() => {
  const root = window.__GBDRAW_APP__.svgContainer.querySelector('svg');
  return [...root.querySelectorAll('text[data-label-feature-id]')]
    .filter((node) => !node.closest('[display="none"]'))
    .map((node) => node.textContent);
});

test('a hidden feature stays listed; its checkbox shows it again, live and after Generate, and Undo reverts', async ({ page }) => {
  test.setTimeout(600_000);
  await openWithGenBank(page, FIXTURE, () => { window.__GBDRAW_APP__.form.labels_mode = 'out'; });
  await generate(page);
  await settle(page);
  await openDrawer(page);
  expect(await listedRows(page)).toEqual([
    ['FL1', true], ['RPT_ONE', true], ['FL2', true], ['dup alpha', true], ['dup beta', true]
  ]);
  // Every row has the virtual-scroll row height (FEATURE_ROW_HEIGHT_PX).
  const rowHeights = await page.locator('[data-feature-list-row]')
    .evaluateAll((rows) => rows.map((row) => row.getBoundingClientRect().height));
  expect(new Set(rowHeights)).toEqual(new Set([83]));
  const shown = await drawnFeatures(page);

  // Uncheck: one History step; FL1 is hidden live and stays listed.
  const before = await undoCount(page);
  await visibilityBox(page, 'alpha protein').uncheck();
  await settle(page);
  expect(await undoCount(page)).toBe(before + 1);
  expect(await drawnFeatures(page)).toEqual(shown.filter((name) => name !== 'FL1'));
  await generate(page);
  await settle(page);
  // Generate does not draw FL1; the list keeps it with an unchecked box.
  expect(await listedRows(page)).toContainEqual(['FL1', false]);
  await expect(visibilityBox(page, 'alpha protein')).not.toBeChecked();
  const hidden = await drawnFeatures(page);
  expect(hidden).not.toContain('FL1');

  // Check: one History step; the rerender draws FL1 live, as Generate does.
  const beforeShow = await undoCount(page);
  await visibilityBox(page, 'alpha protein').check();
  await settle(page);
  expect(await undoCount(page)).toBe(beforeShow + 1);
  await waitForRerender(page);
  await settle(page);
  const live = await drawnFeatures(page);
  expect(live).toContain('FL1');
  await expect(visibilityBox(page, 'alpha protein')).toBeChecked();

  // Undo reverts the toggle: FL1 is hidden again and its box unchecked.
  await page.evaluate(() => window.__GBDRAW_HISTORY__.undo());
  await settle(page);
  await waitForRerender(page);
  expect(await drawnFeatures(page)).toEqual(hidden);
  await expect(visibilityBox(page, 'alpha protein')).not.toBeChecked();
  await page.evaluate(() => window.__GBDRAW_HISTORY__.redo());
  await settle(page);
  await waitForRerender(page);
  expect(await drawnFeatures(page)).toEqual(live);

  await generate(page);
  await settle(page);
  expect(await drawnFeatures(page)).toEqual(live);
  await expect(visibilityBox(page, 'alpha protein')).toBeChecked();
});

test('Search features finds a hidden feature; Open shows its popup with the hidden note, and On shows it', async ({ page }) => {
  test.setTimeout(600_000);
  await openWithGenBank(page, FIXTURE, () => { window.__GBDRAW_APP__.form.labels_mode = 'out'; });
  await generate(page);
  await settle(page);
  await openDrawer(page);
  await visibilityBox(page, 'beta protein').uncheck();
  await settle(page);
  await generate(page);
  await settle(page);
  expect(await drawnFeatures(page)).not.toContain('FL2');

  const search = page.getByRole('searchbox', { name: 'Search features', exact: true });
  await search.fill('FL2');
  await search.press('Enter');
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.previewFeatureSearchMatches.length)).toBe(1);
  await page.getByRole('button', { name: 'Open active feature', exact: true }).click();
  const popup = page.locator('.feature-popup');
  await expect(popup).toBeVisible();
  expect(await page.evaluate(() => window.__GBDRAW_APP__.clickedFeature.feat.locus_tag)).toBe('FL2');
  await expect(popup.locator('[data-feature-hidden-note]')).toHaveText('This feature is hidden. Choose On to show it.');
  // The popup of a feature the Result does not draw has no point to open at;
  // it opens inside the viewport, also on a phone.
  const insideViewport = async () => {
    const box = await popup.boundingBox();
    const viewport = page.viewportSize();
    return box.x >= 0 && box.y >= 0 && box.x + box.width <= viewport.width && box.y < viewport.height;
  };
  expect(await insideViewport()).toBe(true);
  const desktop = page.viewportSize();
  await page.evaluate(() => { window.__GBDRAW_APP__.clickedFeature = null; });
  await page.setViewportSize({ width: 390, height: 844 });
  await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    await app.openFeatureEditorFromList(app.filteredFeatures.find((feature) => feature.locus_tag === 'FL2'), null);
  });
  await expect(popup).toBeVisible();
  expect(await insideViewport()).toBe(true);
  await page.setViewportSize(desktop);

  const before = await undoCount(page);
  await popup.getByRole('combobox', { name: 'Feature visibility', exact: true }).selectOption('on');
  // FL2 has a product, so the scope dialog asks first.
  const thisFeature = page.getByRole('button', { name: /This feature One rule/ });
  await expect(thisFeature).toBeVisible();
  await thisFeature.click();
  await settle(page);
  expect(await undoCount(page)).toBe(before + 1);
  await expect(popup.locator('[data-feature-hidden-note]')).toHaveCount(0);
  await waitForRerender(page);
  await settle(page);
  const live = await drawnFeatures(page);
  expect(live).toContain('FL2');
  await page.evaluate(() => { window.__GBDRAW_APP__.clickedFeature = null; });
  await generate(page);
  await settle(page);
  expect(await drawnFeatures(page)).toEqual(live);
});

test('a Label On kept with Keep feature hidden is drawn when the checkbox shows the feature', async ({ page }) => {
  test.setTimeout(600_000);
  await openWithGenBank(page, FIXTURE, () => { window.__GBDRAW_APP__.form.labels_mode = 'none'; });
  await generate(page);
  await settle(page);
  await openDrawer(page);
  await visibilityBox(page, 'beta protein').uncheck();
  await settle(page);
  await generate(page);
  await settle(page);
  expect(await labels(page)).toEqual([]);

  // Label On for the hidden FL2 from its popup; keep the feature hidden.
  await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const fl2 = app.filteredFeatures.find((feature) => feature.locus_tag === 'FL2');
    await app.openFeatureEditorFromList(fl2, null);
    await window.Vue.nextTick();
    app.clickedFeature.labelVisibility = 'on';
    app.clickedFeature.labelText = 'FL2_KEPT';
    window.__labelOnApply = app.updateClickedFeatureLabelText();
  });
  const dialog = page.getByRole('dialog', { name: 'Feature Is Hidden', exact: true });
  await expect(dialog).toBeVisible({ timeout: 120_000 });
  await dialog.getByRole('button', { name: 'Keep feature hidden', exact: true }).click();
  await page.evaluate(async () => {
    await window.__labelOnApply;
    window.__GBDRAW_APP__.clickedFeature = null;
  });
  await settle(page);
  await waitForRerender(page);
  expect(await labels(page)).toEqual([]);

  await visibilityBox(page, 'beta protein').check();
  await settle(page);
  await waitForRerender(page);
  await settle(page);
  expect(await labels(page)).toEqual(['FL2_KEPT']);
  expect(await drawnFeatures(page)).toContain('FL2');
  await generate(page);
  await settle(page);
  expect(await labels(page)).toEqual(['FL2_KEPT']);
});
