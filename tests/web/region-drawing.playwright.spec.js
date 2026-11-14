// Copy to another drawing (PD-OI-090 (4), PD-OI-115): the drawing menu, a
// record's popup and the records a Shift+drag selects open the dialog of
// app/region-drawing.js, which derives the new drawing from the shown one in
// Python and adds it as one project step. The cases find the new drawing by
// its drawing ID, never by its tab text.
const { execFileSync } = require('node:child_process');
const { writeFileSync } = require('node:fs');
const { generateAndWaitForResult } = require('./helpers/app-lifecycle.cjs');
const { loadSessionFile, openWithGenBank } = require('./helpers/audit-browser.cjs');
const { test, expect } = require('@playwright/test');
const { drawingIds, drawingTab, shownDrawing } = require('./helpers/drawing-bar.cjs');
const { download, generate, load, seeds } = require('./helpers/mode-transition.cjs');
const { shownSavedDrawing } = require('./helpers/saved-session.cjs');
const { readSessionFile } = require('./helpers/session-file.cjs');
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
    expect((await undoLabels(page))[0], 'Create is one project step').toBe('Copy to another drawing');

    await generate(page);
    const drawn = await page.evaluate(() => window.__GBDRAW_APP__.currentResults().map((result) => String(result.content || '')));
    expect(drawn).toHaveLength(1);
    // The record group carries the drawn span of its source record.
    expect(drawn[0]).toContain(`data-gbdraw-record-source-start="${start}"`);
    expect(drawn[0]).toContain(`data-gbdraw-record-source-end="${end}"`);

    // Generate is a step of its own: the first Undo takes the figure back, the second the new drawing.
    const undo = page.getByRole('button', { name: 'Undo', exact: true });
    await undo.click();
    await expect.poll(async () => (await undoLabels(page))[0]).toBe('Copy to another drawing');
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
  const sessionPath = testInfo.outputPath('circular-origin-window.gbdraw');
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
  await page.locator('input[accept*=".gbdraw"]').setInputFiles(sessionPath);
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
  const savedPath = testInfo.outputPath('saved.gbdraw');
  await download(page, 'Save Session', savedPath);
  const saved = await readSessionFile(savedPath);
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

/** A 4 kb circular record with two CDS. @param {string} id */
const smallGenbank = (id) => {
  const sequence = 'ACGTTGCA'.repeat(500);
  const origin = sequence.match(/.{1,60}/g).map((chunk, index) => (
    `${String(index * 60 + 1).padStart(9)} ${chunk.match(/.{1,10}/g).join(' ')}`
  )).join('\n');
  return `LOCUS       ${id}                    4000 bp    DNA     circular UNA 01-JAN-2000
DEFINITION  Two CDS.
ACCESSION   ${id}
VERSION     ${id}
FEATURES             Location/Qualifiers
     CDS             101..600
                     /gene="${id.toLowerCase()}A"
     CDS             complement(1001..1600)
                     /gene="${id.toLowerCase()}B"
ORIGIN
${origin}
//
`;
};

/** The screen box of each drawn record name (its definition group), in drawn order. @param {import('@playwright/test').Page} page */
const recordNameBoxes = (page) => page.evaluate(() => [...document.querySelectorAll(
  '.gbdraw-preview-surface [data-gbdraw-role="record-definition"][data-gbdraw-record-id]'
)].map((group) => {
  const { left, top, right, bottom } = group.getBoundingClientRect();
  return { id: group.getAttribute('data-gbdraw-record-id'), left, top, right, bottom };
}));

/** Rows of the dialog as [record key, checked, kind]. @param {import('@playwright/test').Page} page */
const dialogRows = (page) => page.evaluate(() => window.__GBDRAW_APP__.regionDrawings.dialog.rows
  .map((row) => [row.caption, row.checked, row.kind]));

test('a record popup and the records of a Shift+drag open Copy to another drawing with them whole; two whole records make one Multi-Record Canvas that Undo removes', async ({ page }, testInfo) => {
  test.setTimeout(420_000);
  const input = testInfo.outputPath('two-records.gbk');
  writeFileSync(input, smallGenbank('RECA') + smallGenbank('RECB'));
  await openWithGenBank(page, input, () => { window.__GBDRAW_APP__.form.multi_record_canvas = true; });
  await generateAndWaitForResult(page);
  await page.getByRole('button', { name: 'Fit to preview', exact: true }).click();
  await expect.poll(async () => (await recordNameBoxes(page)).map(({ id }) => id)).toEqual(['RECA', 'RECB']);
  const source = (await shownDrawing(page)).id;
  const dialog = page.locator('[data-region-drawing-dialog]');

  // The record popup: that record chosen whole; one whole record keeps the drawing's mode.
  await page.locator('.gbdraw-preview-surface [data-gbdraw-role="record-definition"][data-gbdraw-record-id="RECA"] text').first().click();
  const copyRecord = page.locator('[data-record-copy-to-drawing]');
  await expect(copyRecord).toBeVisible();
  await copyRecord.click();
  await expect(copyRecord, 'the popup closes as the dialog opens').toHaveCount(0);
  await expect.poll(() => dialogRows(page), { timeout: 180_000 })
    .toEqual([['RECA', true, 'whole'], ['RECB', false, 'whole']]);
  await expect.poll(() => regionStatus(page), { timeout: 180_000 }).toBe('derived');
  await expect(dialog.locator('[data-region-drawing-mode]')).toHaveValue('circular');
  await dialog.locator('[data-region-drawing-cancel]').click();
  await expect(dialog).toHaveCount(0);
  await page.getByRole('button', { name: 'Fit to preview', exact: true }).click();

  // Shift+drag over both names selects the records (and the features in the box); Clear clears both. The box
  // starts and ends on the blank canvas inside the circles: a press on a name opens its popup.
  const shiftDragOverNames = async () => {
    // The fit animates the preview transform: the names are read once they stand still.
    let last = '';
    await expect.poll(async () => {
      const now = JSON.stringify(await recordNameBoxes(page));
      const settled = now === last;
      last = now;
      return settled;
    }, { intervals: [250] }).toBe(true);
    const boxes = JSON.parse(last);
    const left = Math.min(...boxes.map((box) => box.left)) - 24;
    const top = Math.min(...boxes.map((box) => box.top)) - 24;
    const right = Math.max(...boxes.map((box) => box.right)) + 24;
    const bottom = Math.max(...boxes.map((box) => box.bottom)) + 24;
    await page.keyboard.down('Shift');
    await page.mouse.move(left, top);
    await page.mouse.down();
    await page.mouse.move(right, bottom, { steps: 8 });
    await page.mouse.up();
    await page.keyboard.up('Shift');
  };
  await shiftDragOverNames();
  const records = page.locator('[data-selected-records]');
  await expect(records).toContainText('2 records: RECA, RECB');
  await page.getByRole('button', { name: 'Clear selected features', exact: true }).click();
  await expect(records).toHaveCount(0);
  expect(await undoLabels(page), 'a selection is no History step').not.toContain('Copy to another drawing');

  await shiftDragOverNames();
  await records.locator('[data-selected-records-copy]').click();
  await expect.poll(() => dialogRows(page), { timeout: 180_000 })
    .toEqual([['RECA', true, 'whole'], ['RECB', true, 'whole']]);
  await expect.poll(() => regionStatus(page), { timeout: 180_000 }).toBe('derived');
  await expect(dialog.locator('[data-region-drawing-mode]')).toHaveValue('circular');
  await dialog.locator('[data-region-drawing-create]').click();
  await expect(dialog).toHaveCount(0);
  await expect.poll(async () => (await drawingIds(page)).length).toBe(2);
  const [added] = (await drawingIds(page)).filter((id) => id !== source);
  expect(await shownDrawing(page)).toEqual({ id: added, kind: 'circular' });
  expect((await undoLabels(page))[0], 'Copy is one project step').toBe('Copy to another drawing');
  expect(await page.evaluate(() => window.__GBDRAW_APP__.form.multi_record_canvas), 'one Multi-Record Canvas').toBe(true);
  await generateAndWaitForResult(page);
  const drawn = await page.evaluate(() => window.__GBDRAW_APP__.currentResults().map((result) => String(result.content || '')));
  expect(drawn, 'one figure of both records').toHaveLength(1);
  for (const id of ['RECA', 'RECB']) expect(drawn[0]).toContain(`data-gbdraw-record-id="${id}"`);

  const undo = page.getByRole('button', { name: 'Undo', exact: true });
  await undo.click();
  await expect.poll(async () => (await undoLabels(page))[0]).toBe('Copy to another drawing');
  await undo.click();
  await expect.poll(() => drawingIds(page)).toEqual([source]);
});

/** The shown drawing's Linear rows (record captions are their files' record IDs). @param {import('@playwright/test').Page} page */
const linearRowCount = (page) => page.evaluate(async () => (await import('./js/state.js')).state.activeDrawing().linearSeqs.length);
/** The record IDs the shown drawing's Result draws, in drawn order. @param {import('@playwright/test').Page} page */
const drawnRecordIds = (page) => page.evaluate(() => window.__GBDRAW_APP__.currentResults().flatMap((result) => [
  ...new DOMParser().parseFromString(String(result.content || ''), 'image/svg+xml')
    .querySelectorAll('[data-gbdraw-role="record-definition"][data-gbdraw-record-id]')
].map((group) => group.getAttribute('data-gbdraw-record-id'))));
/** Checks the dialog's rows of `captions` whole and waits for the derivation. @param {import('@playwright/test').Page} page */
const checkRecords = async (page, dialog, captions) => {
  for (const caption of captions) {
    await dialog.locator('[data-region-drawing-row]').filter({ hasText: caption }).locator('[data-region-drawing-check]').check();
  }
  await expect.poll(() => regionStatus(page), { timeout: 180_000 }).toBe('derived');
};

// OV-3533 (D-30): the whole records of one file copied to a Circular drawing are one File card, so each keeps the row
// the source drew it on: #1 and #2 of a three-record canvas share its first row.
test('a Circular copy of whole records #1 and #2 of a three-record file draws exactly those two, on one row as the source drew them', async ({ page }, testInfo) => {
  test.setTimeout(420_000);
  const input = testInfo.outputPath('three-records.gbk');
  writeFileSync(input, smallGenbank('RECA') + smallGenbank('RECB') + smallGenbank('RECD'));
  await openWithGenBank(page, input, () => { window.__GBDRAW_APP__.form.multi_record_canvas = true; });
  await generateAndWaitForResult(page);
  await page.getByRole('button', { name: 'Fit to preview', exact: true }).click();
  await expect.poll(async () => (await recordNameBoxes(page)).map(({ id }) => id)).toEqual(['RECA', 'RECB', 'RECD']);
  const rowsOf = async (ids) => {
    const boxes = await recordNameBoxes(page);
    return ids.map((id) => Math.round(boxes.find((box) => box.id === id).top));
  };
  const [topA, topB] = await rowsOf(['RECA', 'RECB']);
  expect(topA, 'the source draws #1 and #2 on one row').toBe(topB);
  const dialog = await openRegionDialog(page);
  await checkRecords(page, dialog, ['RECA', 'RECB']);
  await expect(dialog.locator('[data-region-drawing-mode]')).toHaveValue('circular');
  await dialog.locator('[data-region-drawing-create]').click();
  await expect(dialog).toHaveCount(0);
  await generateAndWaitForResult(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.currentResults().length)).toBe(1);
  expect(await drawnRecordIds(page), 'exactly the two records').toEqual(['RECA', 'RECB']);
  await page.getByRole('button', { name: 'Fit to preview', exact: true }).click();
  await expect.poll(async () => (await recordNameBoxes(page)).map(({ id }) => id)).toEqual(['RECA', 'RECB']);
  const [copyA, copyB] = await rowsOf(['RECA', 'RECB']);
  expect(copyA, 'one row, as the source drew them').toBe(copyB);
});

test('Copy to an existing Linear drawing appends the records as one Undo step of that drawing, which the next Generate draws; Skip it and Add it again; the last target is kept until a Load', async ({ page }, testInfo) => {
  test.setTimeout(600_000);
  const input = testInfo.outputPath('two-records.gbk');
  writeFileSync(input, smallGenbank('RECA') + smallGenbank('RECB'));
  await openWithGenBank(page, input, () => { window.__GBDRAW_APP__.form.multi_record_canvas = true; });
  await generateAndWaitForResult(page);
  const source = (await shownDrawing(page)).id;

  // A new Linear drawing of RECA: the dialog opens on a new drawing at first, then on the drawing that copy made (D-29).
  let dialog = await openRegionDialog(page);
  await expect(dialog.locator('[data-region-drawing-target-new]')).toBeChecked();
  await expect(dialog.locator('[data-region-drawing-target-existing]'), 'no other drawing yet').toBeDisabled();
  await checkRecords(page, dialog, ['RECA']);
  await dialog.locator('[data-region-drawing-mode]').selectOption('linear');
  await expect.poll(() => regionStatus(page), { timeout: 180_000 }).toBe('derived');
  await dialog.locator('[data-region-drawing-create]').click();
  await expect(dialog).toHaveCount(0);
  const target = (await shownDrawing(page)).id;
  expect(target).not.toBe(source);
  expect(await linearRowCount(page)).toBe(1);

  await drawingTab(page, source).click();
  await expect.poll(async () => (await shownDrawing(page)).id).toBe(source);
  dialog = await openRegionDialog(page);
  await expect(dialog.locator('[data-region-drawing-target-existing]'), 'the last target').toBeChecked();
  await expect(dialog.locator('[data-region-drawing-existing]')).toHaveValue(target);
  await checkRecords(page, dialog, ['RECA', 'RECB']);
  const already = dialog.locator('[data-region-drawing-already]');
  await expect(already).toContainText('RECA (whole record) is already in');
  await expect(dialog.locator('[data-region-drawing-skip]'), 'Skip it first').toBeChecked();
  await dialog.locator('[data-region-drawing-create]').click();
  await expect(dialog).toHaveCount(0);
  expect(await shownDrawing(page)).toEqual({ id: target, kind: 'linear' });
  expect(await linearRowCount(page), 'RECB appended, RECA skipped').toBe(2);
  expect((await undoLabels(page))[0], 'one step of the drawing').toBe('Copy to another drawing');
  const undo = page.getByRole('button', { name: 'Undo', exact: true });
  await undo.click();
  await expect.poll(() => linearRowCount(page)).toBe(1);
  expect(await shownDrawing(page)).toEqual({ id: target, kind: 'linear' });
  await page.getByRole('button', { name: 'Redo', exact: true }).click();
  await expect.poll(() => linearRowCount(page)).toBe(2);
  await generateAndWaitForResult(page);
  expect(await drawnRecordIds(page), 'Generate draws the pasted record').toEqual(['RECA', 'RECB']);

  // Add it again adds the record that is there already as another record.
  await drawingTab(page, source).click();
  await expect.poll(async () => (await shownDrawing(page)).id).toBe(source);
  dialog = await openRegionDialog(page);
  await expect(dialog.locator('[data-region-drawing-existing]')).toHaveValue(target);
  await checkRecords(page, dialog, ['RECA']);
  await expect(already).toBeVisible();
  await dialog.locator('[data-region-drawing-add-again]').check();
  await expect.poll(() => regionStatus(page), { timeout: 180_000 }).toBe('derived');
  await expect(already, 'the line stays, with Add it again chosen').toBeVisible();
  await dialog.locator('[data-region-drawing-create]').click();
  await expect(dialog).toHaveCount(0);
  expect(await linearRowCount(page)).toBe(3);
  await generateAndWaitForResult(page);
  expect(await drawnRecordIds(page)).toEqual(['RECA', 'RECB', 'RECA']);

  // A Load replaces the drawings: the dialog opens on a new drawing again.
  const saved = testInfo.outputPath('pasted.gbdraw');
  await download(page, 'Save Session', saved, { title: 'pasted' });
  await loadSessionFile(page, saved);
  dialog = await openRegionDialog(page);
  await expect(dialog.locator('[data-region-drawing-target-new]')).toBeChecked();
});
