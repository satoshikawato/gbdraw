// Find in other files (PD-OI-123): the Region popup of a Linear drawing opens the dialog with the picked span as the
// reference; Search runs tblastn (the vendored LOSAT Wasm) for each other File under the LOSAT overlay; the proposal
// comes from Python's plan; Create adds a new Linear drawing right after the source as one project step, which
// Generate draws and Undo removes. The fixture is six MIBiG cluster records of one Linear drawing.
const { test, expect } = require('@playwright/test');
const { generateAndWaitForResult } = require('./helpers/app-lifecycle.cjs');
const { loadSessionFile, openFresh } = require('./helpers/audit-browser.cjs');
const { drawingIds, drawingTab, shownDrawing } = require('./helpers/drawing-bar.cjs');

test.describe.configure({ retries: 0 });

const SESSION = 'tests/test_inputs/BGC0000708-BGC0000713.gbdraw-session.json';

/** A point of a record's axis on the screen (as region-select.playwright.spec.js reads it). */
const axisPoint = (page, recordIndex, position) => page.evaluate(async ({ recordIndex, position }) => {
  const { state } = await import('./js/state.js');
  const { readRegionAxisRecords, regionAxisFrame } = await import('./js/services/region-axis.js');
  const record = readRegionAxisRecords(state.trackSlotResolvedGeometry.value, state.selectedResultIndex.value)
    .find((entry) => entry.recordIndex === recordIndex);
  const frame = regionAxisFrame(state.svgContainer.value.querySelector('svg'), recordIndex);
  const shown = record.axis.rows[0];
  const x = shown.xStartPx + (shown.xEndPx - shown.xStartPx) * (position - shown.displayStart) / (shown.displayEnd - shown.displayStart);
  const point = new DOMPoint(x, (shown.bandTopPx + shown.bandBottomPx) / 2).matrixTransform(frame.getScreenCTM());
  return { x: point.x, y: point.y };
}, { recordIndex, position });

/** @param {import('@playwright/test').Page} page */
const status = (page) => page.evaluate(() => window.__GBDRAW_APP__.findInOtherFiles.dialog.status);
/** @param {import('@playwright/test').Page} page */
const undoLabels = (page) => page.evaluate(() => window.__GBDRAW_HISTORY__.timeline().undo.map((step) => step.label));

test('Linear: Find in other files proposes windows from a tblastn search and Create adds a drawing that Generate draws and Undo removes', async ({ page }) => {
  test.setTimeout(900_000);
  await openFresh(page);
  await loadSessionFile(page, SESSION);
  await generateAndWaitForResult(page);
  const [source] = await drawingIds(page);

  // Select region on the first record, then the popup's entry.
  await page.getByRole('button', { name: 'Fit to preview', exact: true }).click();
  const toggle = page.locator('[data-region-select-toggle]');
  await toggle.click();
  await expect(toggle).toHaveAttribute('aria-pressed', 'true');
  const from = await axisPoint(page, 0, 2000.5);
  const to = await axisPoint(page, 0, 20000.5);
  await page.mouse.move(from.x, from.y);
  await page.mouse.down();
  await page.mouse.move(to.x, to.y, { steps: 8 });
  await page.mouse.up();
  const popup = page.locator('[data-preview-item-popup="region"]');
  const entry = popup.locator('[data-region-find-in-files]');
  await expect(entry).toBeEnabled();
  // After "Use as search query".
  expect(await popup.locator('button').evaluateAll((buttons) => buttons.map((button) => button.hasAttribute('data-region-search-query')
    ? 'query' : button.hasAttribute('data-region-find-in-files') ? 'find' : '').filter(Boolean))).toEqual(['query', 'find']);
  await entry.click();

  const dialog = page.locator('[data-find-files-dialog]');
  await expect(dialog.locator('[data-find-files-genes]')).toBeVisible({ timeout: 180_000 });
  const row = dialog.locator('[data-find-files-row]');
  await expect(row).toHaveCount(1);
  // The picked span, in source coordinates (this record is drawn from a later display start).
  await expect(row.locator('[data-find-files-record]')).not.toHaveValue('');
  expect(Number(await row.locator('[data-find-files-end]').inputValue()) - Number(await row.locator('[data-find-files-start]').inputValue())).toBe(18_000);
  // Every other File is offered, ticked at first.
  await expect(dialog.locator('[data-find-files-target]')).toHaveCount(5);
  await expect(dialog.locator('[data-find-files-target] input:not(:checked)')).toHaveCount(0);
  await expect(dialog.locator('[data-find-files-program]')).toHaveValue('tlosatn');

  const steps = (await undoLabels(page)).length;
  await dialog.locator('[data-find-files-search]').click();
  await expect.poll(() => status(page), { timeout: 600_000 }).toBe('proposal');
  await expect(dialog).toHaveAttribute('data-find-files-screen', 'proposal');
  await expect(dialog.locator('[data-find-files-file]')).toHaveCount(5);
  const windows = await dialog.locator('[data-find-files-table] [data-find-files-window]').count();
  expect(windows).toBeGreaterThan(0);
  await expect(dialog.locator('[data-find-files-count]')).toContainText(' in ');
  expect((await undoLabels(page)).length, 'the search and the proposal make no History step').toBe(steps);

  // Back and Search again, then Create.
  await dialog.locator('[data-find-files-back]').click();
  await expect(dialog).toHaveAttribute('data-find-files-screen', 'search');
  await dialog.locator('[data-find-files-search]').click();
  await expect.poll(() => status(page), { timeout: 600_000 }).toBe('proposal');
  await dialog.locator('[data-find-files-create]').click();
  await expect(dialog).toHaveCount(0);
  const ids = await drawingIds(page);
  expect(ids).toHaveLength(2);
  const added = ids[1];
  expect(ids[0]).toBe(source);
  expect(await shownDrawing(page)).toEqual({ id: added, kind: 'linear' });
  expect((await undoLabels(page))[0], 'Create is one project step').toBe('Find in other files');

  // Generate draws it: the reference row first, labelled with its File's name.
  await generateAndWaitForResult(page);
  const rows = await page.evaluate(() => window.__GBDRAW_APP__.currentResults().length);
  expect(rows).toBeGreaterThan(0);

  // Generate is a step of its own: the first Undo takes the figure back, the second the new drawing.
  const undo = page.getByRole('button', { name: 'Undo', exact: true });
  await undo.click();
  await expect.poll(async () => (await undoLabels(page))[0]).toBe('Find in other files');
  await undo.click();
  await expect.poll(() => drawingIds(page)).toEqual([source]);
  await expect(drawingTab(page, source)).toHaveAttribute('aria-selected', 'true');
});
