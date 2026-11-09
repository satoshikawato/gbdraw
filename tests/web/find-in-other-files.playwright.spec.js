// Find in other files (PD-OI-123): the Region popup of a Linear drawing opens the dialog with the picked span as the
// reference; Search runs under the full-window overlay from Python's first step to its last, with tblastn (the
// vendored LOSAT Wasm) for each other File; the proposal comes from Python's plan, as a table on a wide screen and as
// cards on a phone; Create adds a new Linear drawing right after the source as one project step, which
// Generate draws and Undo removes. The fixture is five MIBiG cluster records (five Files) of one Linear drawing.
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
  // The fit animates the preview transform: read screen points once it has settled (as region-select's spec does).
  await page.getByRole('button', { name: 'Fit to preview', exact: true }).click();
  let last = '';
  await expect.poll(async () => {
    const now = await page.evaluate(() => JSON.stringify(document.querySelector('.gbdraw-preview-surface').getBoundingClientRect()));
    const settled = now === last;
    last = now;
    return settled;
  }, { intervals: [250] }).toBe(true);
  const toggle = page.locator('[data-region-select-toggle]');
  await toggle.click();
  await expect(toggle).toHaveAttribute('aria-pressed', 'true');
  const from = await axisPoint(page, 0, 2000.5);
  const to = await axisPoint(page, 0, 20000.5);
  await page.mouse.move(from.x, from.y);
  await page.mouse.down();
  await page.mouse.move(to.x, to.y, { steps: 8 });
  await expect(page.locator('[data-region-select-readout]')).toBeVisible();
  await page.mouse.up();
  const popup = page.locator('[data-preview-item-popup="region"]');
  const entry = popup.locator('[data-region-find-in-files]');
  await expect(entry).toBeEnabled();
  // After "Use as search query".
  expect(await popup.locator('button').evaluateAll((buttons) => buttons.map((button) => button.hasAttribute('data-region-search-query')
    ? 'query' : button.hasAttribute('data-region-find-in-files') ? 'find' : '').filter(Boolean))).toEqual(['query', 'find']);
  // The popup's span ("2,001-20,001 (18,001 bp)"): the dialog's reference row is that span.
  const [spanStart, spanEnd] = (await page.locator('[data-region-span]').innerText()).split(' ')[0].split('-').map((text) => text.replace(/,/g, ''));
  await entry.click();

  const dialog = page.locator('[data-find-files-dialog]');
  await expect(dialog.locator('[data-find-files-genes]')).toBeVisible({ timeout: 180_000 });
  const row = dialog.locator('[data-find-files-row]');
  await expect(row).toHaveCount(1);
  await expect(row.locator('[data-find-files-record]')).not.toHaveValue('');
  await expect(row.locator('[data-find-files-start]')).toHaveValue(spanStart);
  await expect(row.locator('[data-find-files-end]')).toHaveValue(spanEnd);
  // Every other File is offered, ticked at first.
  await expect(dialog.locator('[data-find-files-target]')).toHaveCount(4);
  await expect(dialog.locator('[data-find-files-target] input:not(:checked)')).toHaveCount(0);
  await expect(dialog.locator('[data-find-files-program]')).toHaveValue('tlosatn');

  const referenceName = await page.evaluate(() => {
    const view = window.__GBDRAW_APP__.findInOtherFiles.dialog;
    return view.files.find((file) => file.records.some((record) => record.recordKey === view.rows[0].recordKey)).name;
  });
  const steps = (await undoLabels(page)).length;
  await dialog.locator('[data-find-files-search]').click();
  // The overlay covers the app for the whole search, its Python steps included, and the dialog steps aside.
  const overlay = page.locator('[data-generation-progress]');
  await expect(overlay).toBeVisible();
  await expect(dialog).toBeHidden();
  const shots = process.env.FIND_FILES_SHOTS;
  if (shots) {
    // The status line of a Python step, when the poll sees it before LOSAT starts.
    await expect.poll(() => overlay.innerText(), { timeout: 60_000, intervals: [50] }).toMatch(/Reading the files|Running LOSAT/);
    await page.screenshot({ path: `${shots}/overlay-search.png` });
  }
  await expect.poll(() => status(page), { timeout: 600_000 }).toBe('proposal');
  await expect(overlay).toBeHidden();
  await expect(dialog).toHaveAttribute('data-find-files-screen', 'proposal');
  await expect(dialog.locator('[data-find-files-file]')).toHaveCount(4);
  // A screen 768 px wide or wider shows the proposal as a table, not as cards.
  await expect(dialog.locator('[data-find-files-table]').first()).toBeVisible();
  await expect(dialog.locator('[data-find-files-cards]').first()).toBeHidden();
  const windows = await dialog.locator('[data-find-files-table] [data-find-files-window]').count();
  expect(windows).toBeGreaterThan(0);
  await expect(dialog.locator('[data-find-files-count]')).toContainText(' in ');
  expect((await undoLabels(page)).length, 'the search and the proposal make no History step').toBe(steps);
  // A narrower screen, such as a phone, shows each file as a card with one block per window.
  const wide = page.viewportSize();
  await page.setViewportSize({ width: 390, height: 844 });
  await expect(dialog.locator('[data-find-files-cards]').first()).toBeVisible();
  await expect(dialog.locator('[data-find-files-table]').first()).toBeHidden();
  expect(await dialog.locator('[data-find-files-cards] [data-find-files-window]').count()).toBe(windows);
  if (shots) await page.screenshot({ path: `${shots}/proposal-phone.png` });
  await page.setViewportSize(wide);
  await expect(dialog.locator('[data-find-files-table]').first()).toBeVisible();

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

  // Generate draws it, the reference row labelled with its File's name.
  await generateAndWaitForResult(page);
  const rows = await page.evaluate(() => window.__GBDRAW_APP__.currentResults().length);
  expect(rows).toBeGreaterThan(0);
  expect(await page.locator('.gbdraw-preview-surface svg text').allTextContents()).toContain(referenceName);

  // Generate is a step of its own: the first Undo takes the figure back, the second the new drawing.
  const undo = page.getByRole('button', { name: 'Undo', exact: true });
  await undo.click();
  await expect.poll(async () => (await undoLabels(page))[0]).toBe('Find in other files');
  await undo.click();
  await expect.poll(() => drawingIds(page)).toEqual([source]);
  await expect(drawingTab(page, source)).toHaveAttribute('aria-selected', 'true');
  // Redo restores it.
  await page.getByRole('button', { name: 'Redo', exact: true }).click();
  await expect.poll(() => drawingIds(page)).toEqual([source, added]);
});

// OV-2995, OV-3003: a reference span long enough that a File gets two windows of one record. The new drawing's LOSAT
// pairs join drawn units (one per window) that share a reference gene; Generate runs LOSAT in the browser for exactly
// those pairs and draws their ribbons, without a choice in the Comparison panel.
test('Linear: the drawing Create makes with several windows of one record draws with Generate', async ({ page }) => {
  test.setTimeout(900_000);
  await openFresh(page);
  await loadSessionFile(page, SESSION);
  await generateAndWaitForResult(page);
  await page.getByRole('button', { name: 'Fit to preview', exact: true }).click();
  let last = '';
  await expect.poll(async () => {
    const now = await page.evaluate(() => JSON.stringify(document.querySelector('.gbdraw-preview-surface').getBoundingClientRect()));
    const settled = now === last;
    last = now;
    return settled;
  }, { intervals: [250] }).toBe(true);
  const toggle = page.locator('[data-region-select-toggle]');
  await toggle.click();
  const from = await axisPoint(page, 0, 1000.5);
  const to = await axisPoint(page, 0, 40000.5);
  await page.mouse.move(from.x, from.y);
  await page.mouse.down();
  await page.mouse.move(to.x, to.y, { steps: 8 });
  await page.mouse.up();
  await page.locator('[data-preview-item-popup="region"] [data-region-find-in-files]').click();
  const dialog = page.locator('[data-find-files-dialog]');
  await expect(dialog.locator('[data-find-files-genes]')).toBeVisible({ timeout: 180_000 });
  await dialog.locator('[data-find-files-search]').click();
  await expect.poll(() => status(page), { timeout: 600_000 }).toBe('proposal');
  await dialog.locator('[data-find-files-create]').click();
  await expect(dialog).toHaveCount(0);
  // A File of the new drawing draws one record in two or more windows.
  const windowCounts = await page.evaluate(async () => (await import('./js/state.js')).state.activeDrawing().linearSeqs
    .map((seq) => (Array.isArray(seq.windows) ? seq.windows.length : 0)));
  expect(Math.max(...windowCounts)).toBeGreaterThan(1);
  // The window pairs Python kept: the saved pipeline's, as unit indexes.
  const savedPairs = await page.evaluate(async () => (await import('./js/state.js')).state.activeDrawing()
    .importedComparisonIntent.source.comparisons.find((comparison) => comparison.kind === 'generatedProteinComparison')
    .pairs.map((pair) => [pair.queryRecordIndex, pair.subjectRecordIndex]));
  expect(savedPairs.length).toBeGreaterThan(0);
  await generateAndWaitForResult(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.currentResults().length)).toBeGreaterThan(0);
  // The request draws one LOSATP table per saved pair, and no other pair; the figure draws their ribbons.
  const drawn = await page.evaluate(async () => {
    const { getCommittedCanonicalSession } = await import('./js/services/config.js');
    const { state } = await import('./js/state.js');
    return {
      pairs: getCommittedCanonicalSession().renderRequest.comparisons
        .filter((comparison) => comparison.kind === 'precomputedProteinComparison')
        .map((comparison) => [comparison.queryRecordIndex, comparison.subjectRecordIndex]),
      ribbons: state.svgContainer.value.querySelectorAll('svg [data-gbdraw-match-id]').length
    };
  });
  const order = (pairs) => pairs.map((pair) => pair.join(':')).sort();
  expect(order(drawn.pairs)).toEqual(order(savedPairs));
  expect(drawn.ribbons).toBeGreaterThan(0);
});
