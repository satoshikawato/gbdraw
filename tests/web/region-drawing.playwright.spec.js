// New drawing from a region (PD-OI-090 (4)): the drawing menu's "New drawing
// from a region…" opens the dialog of app/region-drawing.js, which derives the
// new drawing from the shown one in Python and adds it as one project step.
// The cases find the new drawing by its drawing ID, never by its tab text.
const { execFileSync } = require('node:child_process');
const { gunzipSync } = require('node:zlib');
const { test, expect } = require('@playwright/test');
const { drawingIds, drawingTab, shownDrawing } = require('./helpers/drawing-bar.cjs');
const { download, generate, load, seeds } = require('./helpers/mode-transition.cjs');
const { shownSavedDrawing } = require('./helpers/saved-session.cjs');
const { expectSessionLoadedWithoutWarnings, openApp } = require('./helpers/app-lifecycle.cjs');
const { chooseCoordinates, openRegionDialog } = require('./helpers/region-drawing.cjs');

test.describe.configure({ retries: 0 });

/** @param {import('@playwright/test').Page} page */
const regionStatus = (page) => page.evaluate(() => window.__GBDRAW_APP__.regionDrawings.dialog.status);

/** @param {import('@playwright/test').Page} page */
const undoLabels = (page) => page.evaluate(() => window.__GBDRAW_HISTORY__.timeline().undo.map((step) => step.label));

test('a region of a Circular record becomes a Linear drawing that Generate draws at the region length and Undo removes after the figure', async ({ browser }) => {
  test.setTimeout(420_000);
  const page = await load(browser, seeds.circular);
  try {
    await generate(page);
    const [source] = await drawingIds(page);
    const start = 1001;
    // 1,000 bp: the Linear Auto scale interval (100 bp) differs from the Circular one of the 16,569 bp record (1,000 bp).
    const end = 2000;

    // Without a size set by hand, a shorter region keeps every setting.
    let dialog = await openRegionDialog(page);
    await chooseCoordinates(dialog, start, end);
    await expect.poll(() => regionStatus(page), { timeout: 60_000 }).toBe('derived');
    const returned = dialog.locator('[data-region-drawing-returned] li');
    const before = await returned.count();
    await dialog.locator('[data-region-drawing-cancel]').click();
    await expect(dialog).toHaveCount(0);
    expect(await drawingIds(page), 'Cancel adds no drawing').toEqual([source]);

    // A scale interval set by hand returns to Auto for the shorter region.
    await page.evaluate(() => { window.__GBDRAW_APP__.adv.scale_interval = 777; });
    dialog = await openRegionDialog(page);
    await chooseCoordinates(dialog, start, end);
    await expect.poll(() => regionStatus(page), { timeout: 60_000 }).toBe('derived');
    await expect(dialog.locator('[data-region-drawing-returned]')).toBeVisible();
    await expect(returned).toHaveCount(before + 1);
    await expect(returned.filter({ hasText: '777' })).toHaveCount(1);
    const notCarried = dialog.locator('[data-region-drawing-not-carried] li');
    await expect(notCarried.or(dialog.getByText('Every setting is carried over.', { exact: true })).first()).toBeVisible();

    const create = dialog.locator('[data-region-drawing-create]');
    await expect(create).toBeEnabled();
    await create.click();
    await expect(dialog).toHaveCount(0);
    await expect.poll(async () => (await drawingIds(page)).length).toBe(2);
    const [added] = (await drawingIds(page)).filter((id) => id !== source);
    await expect(drawingTab(page, added)).toBeVisible();
    await expect(drawingTab(page, added)).toHaveAttribute('aria-selected', 'true');
    expect(await shownDrawing(page)).toEqual({ id: added, kind: 'linear' });
    expect((await undoLabels(page))[0], 'Create is one project step').toBe('New drawing from a region');

    await generate(page);
    const drawn = await page.evaluate(() => window.__GBDRAW_APP__.currentResults().map((result) => String(result.content || '')));
    expect(drawn).toHaveLength(1);
    // The record group carries the drawn span of its source record.
    expect(drawn[0]).toContain(`data-gbdraw-record-source-start="${start}"`);
    expect(drawn[0]).toContain(`data-gbdraw-record-source-end="${end}"`);

    // Generate is a step of its own: the first Undo takes the figure back, the second the new drawing.
    const undo = page.getByRole('button', { name: 'Undo', exact: true });
    await undo.click();
    await expect.poll(async () => (await undoLabels(page))[0]).toBe('New drawing from a region');
    expect(await drawingIds(page), 'Undo of Generate keeps the new drawing').toEqual([source, added]);
    await undo.click();
    await expect.poll(() => drawingIds(page)).toEqual([source]);
    await expect(drawingTab(page, added)).toHaveCount(0);
  } finally {
    await page.context().close();
  }
});

// OV-2489: after a Load the Circular record list stays deferred until Generate; the dialog reads it.
test('OV-2489: after a Session Load without Generate the dialog offers the drawing\'s records', async ({ browser }) => {
  test.setTimeout(300_000);
  const page = await load(browser, seeds.circular);
  try {
    const dialog = await openRegionDialog(page);
    await expect(dialog.locator('[data-region-drawing-row]')).toHaveCount(1);
  } finally {
    await page.context().close();
  }
});

test('Coordinates that start after their end are refused on a linear record with the reason, Create stays disabled, and no drawing is added', async ({ browser }) => {
  test.setTimeout(300_000);
  // Lambda (NC_001416) is a linear record: a start after the end names no span.
  const page = await load(browser, seeds.linear);
  try {
    const ids = await drawingIds(page);
    const steps = await undoLabels(page);
    const dialog = await openRegionDialog(page);
    await chooseCoordinates(dialog, 5000, 100);
    await expect.poll(() => regionStatus(page), { timeout: 60_000 }).toBe('refused');
    const message = dialog.locator('[data-region-drawing-message]');
    await expect(message).toHaveText('Choose a start before the end.');
    await expect(dialog.locator('[data-region-drawing-create]')).toBeDisabled();
    await dialog.locator('[data-region-drawing-cancel]').click();
    await expect(dialog).toHaveCount(0);
    expect(await drawingIds(page)).toEqual(ids);
    expect(await undoLabels(page)).toEqual(steps);
  } finally {
    await page.context().close();
  }
});

test('Coordinates across the origin of a circular record make one window across it, which Generate draws', async ({ browser }) => {
  test.setTimeout(420_000);
  const page = await load(browser, seeds.circular);
  try {
    const [source] = await drawingIds(page);
    const dialog = await openRegionDialog(page);
    await chooseCoordinates(dialog, 16000, 100);
    await expect.poll(() => regionStatus(page), { timeout: 60_000 }).toBe('derived');
    await dialog.locator('[data-region-drawing-create]').click();
    await expect(dialog).toHaveCount(0);
    const [added] = (await drawingIds(page)).filter((id) => id !== source);
    expect(await shownDrawing(page)).toEqual({ id: added, kind: 'linear' });
    // The card's Region block holds the window with "Across the origin" ticked (PD-OI-124): Start 16,000, End 100.
    const region = page.locator('[data-record-region-fields]');
    // The File card may be folded: read the fields by their labels, shown or not.
    await expect(region.locator('input[type="checkbox"][aria-label="Across the origin for sequence 1"]')).toBeChecked();
    await expect(region.getByLabel('Region start for sequence 1')).toHaveValue('16000');
    await expect(region.getByLabel('Region end for sequence 1')).toHaveValue('100');
    await generate(page);
    // Once the record is read, the block gives the bases it draws across the origin of the 16,569 bp record.
    await expect(region.locator('[data-record-window-origin]'))
      .toHaveText('16,000 \u2192 16,569 \u2192 1 \u2192 100: 670 bp across the origin.');
    const drawn = await page.evaluate(() => window.__GBDRAW_APP__.currentResults().map((result) => String(result.content || '')));
    expect(drawn).toHaveLength(1);
    // The record group carries the drawn span: it ends before it starts, across the origin.
    expect(drawn[0]).toContain('data-gbdraw-record-source-start="16000"');
    expect(drawn[0]).toContain('data-gbdraw-record-source-end="100"');
  } finally {
    await page.context().close();
  }
});

// Review S3 (PD-OI-122): a window of a circular record may cross the origin. A Circular Session the Python API saved
// with RecordWindow(16000, 100) loads that window into the record's read-only list (the Start/End fields would read
// 16,000-100 as an error), Generate draws the 670 bp window across the origin, and Save writes the window back.
test('A Circular Session whose record has one window across the origin loads it, Generate draws it, and Save keeps it', async ({ page }, testInfo) => {
  test.setTimeout(420_000);
  const sessionPath = testInfo.outputPath('circular-origin-window.gbdraw-session.json');
  execFileSync('python', ['-c', [
    'import sys',
    'from gbdraw.api import save_session_document',
    'from gbdraw.api.requests import CircularDiagramRequest, GenBankInputSource, RecordInput',
    'from gbdraw.io.regions import RecordWindow',
    'source = GenBankInputSource("tests/test_inputs/HmmtDNA.gbk")',
    'request = CircularDiagramRequest(records=(RecordInput(source, windows=(RecordWindow(16000, 100),), record_key="mito"),))',
    'save_session_document(sys.argv[1], request, title="circular-origin-window", overwrite=True)'
  ].join('\n'), sessionPath], { cwd: process.cwd(), stdio: 'inherit' });
  // The Session holds no Result: Load ends with its notice, without the features a Result lists.
  await openApp(page);
  await page.locator('input[accept^=".json,"]').setInputFiles(sessionPath);
  await expectSessionLoadedWithoutWarnings(page);
  expect(await shownDrawing(page)).toMatchObject({ kind: 'circular' });
  // Generate draws the window across the origin of the 16,569 bp record: 16,000 to the end, then 1 to 100.
  await generate(page);
  const drawn = await page.evaluate(() => window.__GBDRAW_APP__.currentResults().map((result) => String(result.content || '')));
  expect(drawn).toHaveLength(1);
  // A Circular drawing prints its drawn length in the center: 670 bp, not the record's 16,569 bp.
  expect(drawn[0]).toMatch(/>670 bp<\/text>/);
  expect(drawn[0]).not.toContain('16,569 bp');
  // The drawing holds the window in its read-only list, not in the Start/End fields, which would refuse 16,000-100.
  const originWindow = { start: 16000, end: 100, reverseComplement: false, column: null };
  // The form's window has its "Across the origin" box ticked (PD-OI-124); the request writes the four fields above.
  const formWindow = { ...originWindow, acrossOrigin: true };
  expect(await page.evaluate(() => {
    const { circular_region_start, circular_region_end, circular_windows } = window.__GBDRAW_APP__.form;
    return { circular_region_start, circular_region_end, circular_windows };
  })).toEqual({ circular_region_start: null, circular_region_end: null, circular_windows: [formWindow] });
  await page.locator('[data-circular-record-presentation] summary').click();
  await expect(page.getByLabel('Circular region start')).toHaveCount(0);
  await expect(page.locator('[data-circular-record-presentation] [data-record-window]')).toHaveText([/^16,000\u2013100\s*670 bp$/]);
  await page.locator('main svg').last().screenshot({ path: testInfo.outputPath('circular-origin-window.png') });
  await page.screenshot({ path: testInfo.outputPath('circular-origin-window-page.png') });
  // Save writes the window back unchanged.
  const bytes = await download(page, 'Save Session', testInfo.outputPath('saved.gbdraw-session.json'));
  const saved = JSON.parse((bytes[0] === 0x1f ? gunzipSync(bytes) : bytes).toString('utf8'));
  expect(shownSavedDrawing(saved).renderRequest.records[0].windows).toEqual([originWindow]);

  // Review X1 of #1750: the listed window holds the record's orientation. The record popup's Reverse complement
  // reads and writes that window, and Generate draws the window reversed.
  await page.locator('.gbdraw-preview-surface [data-gbdraw-role^="record-definition"] text').first().click();
  const reverse = page.locator('.preview-item-popup[role="dialog"]').getByRole('checkbox', { name: 'Reverse complement', exact: true });
  await expect(reverse).not.toBeChecked();
  await expect(reverse).toBeEnabled();
  await reverse.check();
  expect(await page.evaluate(() => {
    const { circular_reverse, circular_windows } = window.__GBDRAW_APP__.form;
    return { circular_reverse, circular_windows };
  })).toEqual({ circular_reverse: false, circular_windows: [{ ...formWindow, reverseComplement: true }] });
  await generate(page);
  const committed = await page.evaluate(async () => (await import('/gbdraw/web/js/services/config.js')).getCommittedCanonicalRenderRequest());
  expect([committed.records[0].windows, committed.records[0].presentation.reverseComplement])
    .toEqual([[{ ...originWindow, reverseComplement: true }], false]);
  await expect(page.locator('[data-circular-record-presentation] [data-record-window]'))
    .toHaveText([/^16,000\u2013100 \(Reverse complement\)\s*670 bp$/]);
});
