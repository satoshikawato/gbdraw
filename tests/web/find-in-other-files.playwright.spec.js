// Find in other files (PD-OI-123): the Region popup of a Linear drawing opens the dialog with the picked span as the
// reference; Search runs under the full-window overlay from Python's first step to its last, with tblastn (the
// vendored LOSAT Wasm) for each other File; the proposal comes from Python's plan, as a table on a wide screen and as
// cards on a phone; Create adds a new Linear drawing right after the source as one project step, which
// Generate draws and Undo removes. The fixture is five MIBiG cluster records (five Files) of one Linear drawing.
const { test, expect } = require('@playwright/test');
const { generateAndWaitForResult, evaluateWithRetainedPromise, routeCrossOriginIsolated } = require('./helpers/app-lifecycle.cjs');
const { loadSessionFile, openFresh } = require('./helpers/audit-browser.cjs');
const { drawingIds, drawingTab, showDrawing, shownDrawing } = require('./helpers/drawing-bar.cjs');

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

// The pVA1 plasmid of two Vibrio parahaemolyticus strains (NZ_CP065372.1, 20-082A3; NZ_CP046414.1, 19-021-D1, stored
// from another start) and three contigs of a draft assembly (CAVP11, GCF_019780985.1: NZ_JABKAY010000079.1, 80 and 86)
// that hold the same plasmid in pieces, some of them reverse-complemented.
const VIBRIO = {
  '20-082A3': { file: 'tests/test_inputs/NZ_CP065372.1.gbk', length: 67093 },
  '19-021-D1': { file: 'tests/test_inputs/NZ_CP046414.1.gbk', length: 68848 },
  contigs: 'tests/test_inputs/CAVP11_contigs_79_80_86.gbk'
};

/**
 * One Linear drawing of a reference plasmid, the other plasmid and the contigs; Find in other files with the whole
 * reference plasmid (cut at `start` when the span crosses the origin), Search, Create, then Generate of the new drawing.
 * Generate succeeds; returns the new drawing's cards, the records its request commits, and whether a change is pending.
 */
const findWholePlasmid = async (page, { reference, other, start, end }) => {
  await routeCrossOriginIsolated(page.context());
  await openFresh(page);
  await showDrawing(page, { kind: 'linear' });
  await page.getByTestId('linear-genbank-1').setInputFiles(VIBRIO[reference].file);
  await generateAndWaitForResult(page);
  for (const [index, file] of [[2, VIBRIO[other].file], [3, VIBRIO.contigs]]) {
    await page.evaluate(() => window.__GBDRAW_APP__.addLinearSeq());
    await page.getByTestId(`linear-genbank-${index}`).setInputFiles(file);
  }
  await page.getByRole('button', { name: 'Fit to preview', exact: true }).click();
  let last = '';
  await expect.poll(async () => {
    const now = await page.evaluate(() => JSON.stringify(document.querySelector('.gbdraw-preview-surface').getBoundingClientRect()));
    const settled = now === last;
    last = now;
    return settled;
  }, { intervals: [250] }).toBe(true);
  await page.locator('[data-region-select-toggle]').click();
  const from = await axisPoint(page, 0, 1);
  const to = await axisPoint(page, 0, VIBRIO[reference].length);
  await page.mouse.move(from.x + 1, from.y);
  await page.mouse.down();
  await page.mouse.move(to.x - 1, to.y, { steps: 8 });
  await page.mouse.up();
  await page.locator('[data-preview-item-popup="region"] [data-region-find-in-files]').click();
  const dialog = page.locator('[data-find-files-dialog]');
  await expect(dialog.locator('[data-find-files-genes]')).toBeVisible({ timeout: 300_000 });
  const row = dialog.locator('[data-find-files-row]');
  await row.locator('[data-find-files-start]').fill(String(start));
  await row.locator('[data-find-files-end]').fill(String(end));
  const origin = row.locator('[data-find-files-origin]');
  if ((await origin.isChecked()) !== (start > end)) await origin.click();
  await expect(dialog.locator('[data-find-files-search]')).toBeEnabled({ timeout: 120_000 });
  await dialog.locator('[data-find-files-search]').click();
  await expect.poll(() => status(page), { timeout: 1_200_000 }).toBe('proposal');
  await dialog.locator('[data-find-files-create]').click();
  await expect(dialog).toHaveCount(0);
  const outcome = await evaluateWithRetainedPromise(page, async () => {
    const result = await window.__GBDRAW_APP__.runAnalysis();
    return JSON.parse(JSON.stringify({ status: result?.status, code: result?.error?.code ?? null }));
  });
  expect(outcome).toEqual({ status: 'ok', code: null });
  await expect.poll(() => page.evaluate(async () => Boolean((await import('./js/services/config.js'))
    .getCommittedCanonicalRenderRequest()))).toBe(true);
  const drawn = await page.evaluate(async () => {
    const { getCommittedCanonicalRenderRequest } = await import('./js/services/config.js');
    const { state } = await import('./js/state.js');
    return {
      cards: state.activeDrawing().linearSeqs.map((card) => ({ id: card.region_record_id, reverse: card.region_reverse,
        windows: (card.windows || []).map((window) => [window.start, window.end, window.reverseComplement]) })),
      records: getCommittedCanonicalRenderRequest().records.map((record) => ({ id: record.selector.value,
        reverse: record.presentation.reverseComplement,
        windows: record.windows.map((window) => [window.start, window.end, window.reverseComplement]) })),
      pending: window.__GBDRAW_APP__.recordDisplayControls.hasPendingChanges.value
    };
  });
  return drawn;
};

// OV-3365, OV-3375: a File of the new drawing draws contig 80 in two (-) windows; the first Generate draws the drawing
// (the saved comparison's records and the cards are one list), and no change is pending after it.
test('Linear: the Vibrio whole-plasmid drawing Create makes draws with its first Generate and leaves nothing pending', async ({ page }) => {
  test.setTimeout(1_800_000);
  const found = await findWholePlasmid(page, { reference: '20-082A3', other: '19-021-D1', start: 1, end: 67093 });
  const split = found.cards.find((card) => card.id === 'NZ_JABKAY010000080.1');
  expect(split.windows.length).toBeGreaterThan(1);
  // The request draws the record's windows in the order and orientation of its card.
  expect(found.records.find((record) => record.id === split.id).windows).toEqual(split.windows);
  expect(found.pending).toBe(false);
  await expect(page.getByText('changes pending Generate')).toHaveCount(0);
});

// OV-3376: the proposal draws contig 79 whole and (-) here; Generate draws it reverse-complemented, as its card says.
test('Linear: a (-) window that covers a whole record is drawn reverse-complemented', async ({ page }) => {
  test.setTimeout(1_800_000);
  const found = await findWholePlasmid(page, { reference: '19-021-D1', other: '20-082A3', start: 63090, end: 63089 });
  const card = found.cards.find((entry) => entry.id === 'NZ_JABKAY010000079.1');
  expect(card).toEqual({ id: 'NZ_JABKAY010000079.1', reverse: true, windows: [] });
  expect(found.records.find((record) => record.id === 'NZ_JABKAY010000079.1').reverse).toBe(true);
  expect(found.pending).toBe(false);
  await expect(page.getByText('changes pending Generate')).toHaveCount(0);
});
