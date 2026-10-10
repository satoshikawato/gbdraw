// Record selection (PD-OI-091, DESIGN 9): each record of a multi-record file
// is drawn or left out per drawing. The cases build their GenBank files here.
const { test, expect } = require('@playwright/test');
const {
  evaluateWithRetainedPromise,
  generateAndWaitForResult,
  openApp
} = require('./helpers/app-lifecycle.cjs');

test.describe.configure({ retries: 0 });

const record = (id, { repeats = 100, trna = false } = {}) => {
  const sequence = 'atg'.repeat(repeats);
  const origin = sequence.match(/.{1,60}/g).map((chunk, index) => (
    `${String(index * 60 + 1).padStart(9)} ${chunk.match(/.{1,10}/g).join(' ')}`
  )).join('\n');
  const features = [
    '     CDS             1..90',
    `                     /locus_tag="${id}_1"`,
    ...(trna ? ['     tRNA            121..180', `                     /locus_tag="${id}_t"`] : [])
  ].join('\n');
  return `LOCUS       ${id.padEnd(24)} ${sequence.length} bp    DNA     linear   UNA 01-JAN-2000
DEFINITION  record selection browser test.
ACCESSION   ${id}
VERSION     ${id}
KEYWORDS    .
SOURCE      synthetic construct
  ORGANISM  synthetic construct
            .
FEATURES             Location/Qualifiers
${features}
ORIGIN
${origin}
//
`;
};
const genbank = (records) => records.map((entry) => (typeof entry === 'string' ? record(entry) : record(entry.id, entry))).join('');

const observeRequests = (page) => page.addInitScript(() => {
  window.__GBDRAW_REQUESTS__ = [];
  const nativePostMessage = Worker.prototype.postMessage;
  Worker.prototype.postMessage = function trackedRequest(message, transfer) {
    if (message?.type === 'run' && message.payload?.request) {
      window.__GBDRAW_REQUESTS__.push(structuredClone(message.payload.request));
    }
    return transfer === undefined ? nativePostMessage.call(this, message) : nativePostMessage.call(this, message, transfer);
  };
});
const lastRequest = (page) => page.evaluate(() => structuredClone(window.__GBDRAW_REQUESTS__.at(-1) || null));
const drawnRecordIds = (page) => page.evaluate(() => window.__GBDRAW_APP__.results.map((result) => {
  const svg = new DOMParser().parseFromString(String(result.content || ''), 'image/svg+xml');
  return [...new Set([...svg.querySelectorAll('g[data-gbdraw-record-id]')]
    .map((element) => element.getAttribute('data-gbdraw-record-id')))];
}));

const openLinear = async (page) => {
  await observeRequests(page);
  await openApp(page);
  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    app.setDiagramMode('linear');
    app.lInputType = 'gb';
    Object.assign(app.form, { show_gc: false, show_skew: false, show_labels_linear: 'none' });
  });
};
// Uploads a File into card `index` (or a new File card) as the file input does.
const uploadLinear = (page, name, text, { newFile = false, index = 0 } = {}) => evaluateWithRetainedPromise(page, async (args) => {
  const app = window.__GBDRAW_APP__;
  if (args.newFile) app.addLinearSeq();
  const target = args.newFile ? app.linearSeqs.length - 1 : args.index;
  app.setLinearSeqPrimaryFile(target, 'gb', new File([args.text], args.name, { type: 'text/plain', lastModified: 1 }));
  await app.refreshLinearRecordSelectors();
}, { name, text, newFile, index });
const linearState = (page) => page.evaluate(() => {
  const app = window.__GBDRAW_APP__;
  return app.linearSeqs.map((seq) => ({
    id: seq.region_record_id, drawn: app.recordSelection.isDrawn('linear', seq.uid)
  }));
});
const card = (page, recordId) => page.locator(`[data-linear-record-card="${recordId}"]`);
const cardOf = async (page, recordId) => card(page, await page.evaluate(
  (id) => window.__GBDRAW_APP__.linearSeqs.find((seq) => seq.region_record_id === id).uid, recordId
));
const openRecords = async (page, fileIndex = 0) => {
  const records = page.locator('[data-linear-source-card]').nth(fileIndex).locator('[data-linear-source-records]');
  if (await records.getAttribute('open') === null) await records.locator(':scope > summary').press('Enter');
  return records;
};
const drawCheckbox = async (page, recordId) => {
  const options = (await cardOf(page, recordId)).locator('[data-linear-record-options]');
  if (await options.getAttribute('open') === null) await options.locator(':scope > summary').click();
  return options.getByRole('checkbox', { name: 'Draw this record', exact: true });
};
const undo = (page) => evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
const saveSession = async (page, title) => {
  const download = page.waitForEvent('download');
  await evaluateWithRetainedPromise(page, async (value) => {
    window.__GBDRAW_APP__.sessionTitle = value;
    await window.__GBDRAW_APP__.saveSessionWithTitle();
  }, title);
  return (await download).path();
};
const loadSession = async (page, path) => {
  const loaded = page.waitForEvent('dialog');
  await page.locator('input[accept^=".json,"]').first().setInputFiles(path);
  const dialog = await loaded;
  expect(dialog.message()).toBe('Session loaded successfully!');
  await dialog.accept();
};

test('a Linear record turned OFF leaves the request, its pairs, and its row, and returns to its row when ON', async ({ page }) => {
  test.setTimeout(300_000);
  await openLinear(page);
  await uploadLinear(page, 'three.gb', genbank(['R1', 'R2', 'R3']));
  await expect.poll(() => linearState(page)).toEqual([
    { id: 'R1', drawn: true }, { id: 'R2', drawn: true }, { id: 'R3', drawn: true }
  ]);
  await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    Object.assign(app.adv, { min_bitscore: 0, evalue: 1, identity: 0, alignment_length: 0 });
    app.losatExecution.executionMode = 'serial';
    await app.setLinearRecordLayoutEnabled(true);
    app.linearSeqs.forEach((seq, index) => app.setLinearRecordRow(seq.uid, index + 1));
    await app.setLinearComparisonGlobalAction('losat');
  });
  const records = await openRecords(page);
  await expect(records.locator(':scope > summary')).toHaveText(/Number of records: 3 · 3 drawn/);

  const draw = await drawCheckbox(page, 'R2');
  await draw.uncheck();
  await expect(records.locator(':scope > summary')).toHaveText(/Number of records: 3 · 2 drawn/);
  const r2Options = (await cardOf(page, 'R2')).locator('[data-linear-record-options]');
  await r2Options.locator(':scope > summary').click();
  await expect(r2Options).not.toHaveAttribute('open', '');
  await expect(r2Options.locator('[data-record-off-badge]')).toHaveText('OFF');
  await expect(r2Options.locator(':scope > summary')).toHaveClass(/opacity-60/);

  await generateAndWaitForResult(page);
  const off = await lastRequest(page);
  expect(off.records.map((entry) => entry.presentation.gridRow)).toEqual([1, 3]);
  expect(off.comparisons.map((entry) => [entry.queryRecordIndex, entry.subjectRecordIndex])).toEqual([[0, 1]]);
  expect(await drawnRecordIds(page)).toEqual([['R1', 'R3']]);

  await (await drawCheckbox(page, 'R2')).check();
  await generateAndWaitForResult(page);
  const on = await lastRequest(page);
  expect(on.records.map((entry) => entry.presentation.gridRow)).toEqual([1, 2, 3]);
  expect(on.comparisons.map((entry) => [entry.queryRecordIndex, entry.subjectRecordIndex])).toEqual([[0, 1], [1, 2]]);
  expect(await drawnRecordIds(page)).toEqual([['R1', 'R2', 'R3']]);
});

test('an upload with 25 records opens its list; search, sort, bulk OFF, and the D-06 question', async ({ page, browser }) => {
  test.setTimeout(300_000);
  await openLinear(page);
  const ids = Array.from({ length: 25 }, (_, index) => `contig_${index + 1}`);
  await uploadLinear(page, 'draft.gb', genbank(ids.map((id, index) => ({ id, repeats: 50 + index }))));
  const dialog = page.locator('[data-record-list-dialog]');
  await expect(dialog).toBeVisible();
  await expect(dialog.locator('[data-record-list-count]')).toHaveText('25 of 25 records drawn');
  const rowIds = () => dialog.locator('[data-record-list-row]').evaluateAll((inputs) => inputs.map(
    (input) => input.closest('li').querySelector('span[title]').getAttribute('title')
  ));

  await dialog.getByRole('combobox', { name: 'Sort', exact: true }).selectOption('id-asc');
  expect((await rowIds()).slice(0, 3)).toEqual(['contig_1', 'contig_2', 'contig_3']);
  await dialog.getByRole('combobox', { name: 'Sort', exact: true }).selectOption('length-desc');
  expect((await rowIds())[0]).toBe('contig_25');
  await dialog.getByLabel('Search record IDs', { exact: true }).fill('contig_1');
  expect((await rowIds()).sort()).toEqual(['contig_1', ...ids.slice(9, 19)].sort());
  await dialog.getByRole('button', { name: 'Select none', exact: true }).click();
  await expect(dialog.locator('[data-record-list-count]')).toHaveText('14 of 25 records drawn');
  // The drawing order is the file order whatever the list shows.
  expect((await linearState(page)).map((entry) => entry.id)).toEqual(ids);

  // Select none on every row would empty the File: the question, and Cancel changes nothing.
  await dialog.getByLabel('Search record IDs', { exact: true }).fill('');
  await dialog.getByRole('button', { name: 'Select none', exact: true }).click();
  const question = page.locator('[data-record-remove-file-dialog]');
  await expect(question).toBeVisible();
  await expect(question.getByRole('heading')).toHaveText('Remove the whole file?');
  await expect(question).toContainText('draft.gb');
  await question.getByRole('button', { name: 'Cancel', exact: true }).click();
  await expect(question).toBeHidden();
  await expect(dialog.locator('[data-record-list-count]')).toHaveText('14 of 25 records drawn');
  await dialog.getByRole('button', { name: 'Close', exact: true }).click();
  await expect(dialog).toBeHidden();

  // A Session load reads the File again but opens no list.
  const saved = await saveSession(page, 'record-list');
  const fresh = await browser.newPage();
  try {
    await openApp(fresh);
    await loadSession(fresh, saved);
    await expect.poll(() => fresh.evaluate(() => window.__GBDRAW_APP__.linearSeqs
      .filter((seq) => window.__GBDRAW_APP__.recordSelection.isDrawn('linear', seq.uid)).length)).toBe(14);
    await expect(fresh.locator('[data-record-list-dialog]')).toHaveCount(0);
  } finally {
    await fresh.close();
  }
});

test('D-06 Remove File deletes one of two Files, and clears the only File', async ({ page }) => {
  test.setTimeout(180_000);
  await openLinear(page);
  await uploadLinear(page, 'first.gb', genbank(['A1', 'A2']));
  await uploadLinear(page, 'second.gb', genbank(['B1', 'B2']), { newFile: true });
  await expect.poll(() => linearState(page)).toHaveLength(4);
  await openRecords(page, 0);
  await (await drawCheckbox(page, 'A1')).uncheck();
  // The last drawn record of the File asks first and stays checked.
  await (await drawCheckbox(page, 'A2')).click();
  const question = page.locator('[data-record-remove-file-dialog]');
  await expect(question).toContainText('first.gb');
  await question.getByRole('button', { name: 'Remove File', exact: true }).click();
  await expect(question).toBeHidden();
  await expect.poll(() => linearState(page)).toEqual([{ id: 'B1', drawn: true }, { id: 'B2', drawn: true }]);
  // A1 stays OFF only as long as its File: the removal took the whole File.
  await undo(page);
  await expect.poll(() => linearState(page)).toEqual([
    { id: 'A1', drawn: false }, { id: 'A2', drawn: true }, { id: 'B1', drawn: true }, { id: 'B2', drawn: true }
  ]);

  await page.locator('[data-linear-source-card]').nth(1).locator('[data-record-list-open]').click();
  const list = page.locator('[data-record-list-dialog]');
  await list.getByRole('button', { name: 'Select none', exact: true }).click();
  await expect(question).toContainText('second.gb');
  await question.getByRole('button', { name: 'Remove File', exact: true }).click();
  await expect(list).toBeHidden();
  await expect.poll(() => linearState(page)).toEqual([{ id: 'A1', drawn: false }, { id: 'A2', drawn: true }]);
  await openRecords(page, 0);
  await (await drawCheckbox(page, 'A2')).click();
  await question.getByRole('button', { name: 'Remove File', exact: true }).click();
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.linearSeqs.map((seq) => Boolean(seq.gb)))).toEqual([false]);
});

test('an OFF card keeps its Definition through Save and Load; Delete settings clears it and Undo restores it', async ({ page, browser }) => {
  test.setTimeout(180_000);
  await openLinear(page);
  await uploadLinear(page, 'pair.gb', genbank(['D1', 'D2']));
  await expect.poll(() => linearState(page)).toHaveLength(2);
  await openRecords(page);
  await (await drawCheckbox(page, 'D2')).uncheck();
  const definition = (await cardOf(page, 'D2')).getByLabel('Definition for sequence 2', { exact: true });
  await definition.fill('Kept while OFF');
  const saved = await saveSession(page, 'off-definition');

  const fresh = await browser.newPage();
  try {
    await openApp(fresh);
    await loadSession(fresh, saved);
    await expect.poll(() => linearState(fresh)).toEqual([{ id: 'D1', drawn: true }, { id: 'D2', drawn: false }]);
    expect(await fresh.evaluate(() => window.__GBDRAW_APP__.linearSeqs[1].definition)).toBe('Kept while OFF');
    await openRecords(fresh);
    const d2 = await cardOf(fresh, 'D2');
    const options = d2.locator('[data-linear-record-options]');
    if (await options.getAttribute('open') === null) await options.locator(':scope > summary').click();
    await d2.getByRole('button', { name: 'Delete settings', exact: true }).click();
    await expect.poll(() => fresh.evaluate(() => window.__GBDRAW_APP__.linearSeqs[1].definition)).toBe('');
    await undo(fresh);
    await expect.poll(() => fresh.evaluate(() => window.__GBDRAW_APP__.linearSeqs[1].definition)).toBe('Kept while OFF');
  } finally {
    await fresh.close();
  }
});

test('Circular draws the ON records in the canvas and as separate diagrams', async ({ page }) => {
  test.setTimeout(300_000);
  await observeRequests(page);
  await openApp(page);
  await page.evaluate((text) => {
    const app = window.__GBDRAW_APP__;
    app.setDiagramMode('circular');
    app.setCircularSourceFile('c_gb', new File([text], 'four.gb', { type: 'text/plain', lastModified: 1 }));
    Object.assign(app.form, { multi_record_canvas: true, suppress_gc: true, suppress_skew: true, labels_mode: 'none', legend: 'none' });
  }, genbank(['C1', 'C2', 'C3', 'C4']));
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.circularRecordList.length)).toBe(4);
  await page.locator('[data-record-list-open="circular"]').click();
  const list = page.locator('[data-record-list-dialog]');
  await list.locator('[data-record-list-row="#1"]').uncheck();
  await list.locator('[data-record-list-row="#3"]').uncheck();
  await expect(list.locator('[data-record-list-count]')).toHaveText('2 of 4 records drawn');
  await list.getByRole('button', { name: 'Close', exact: true }).click();
  await expect(page.locator('[data-circular-records-drawn]')).toHaveText('2 of 4 records drawn');

  await generateAndWaitForResult(page);
  expect((await lastRequest(page)).records.map((entry) => entry.recordKey)).toEqual(['record-2', 'record-4']);
  expect(await drawnRecordIds(page)).toEqual([['C2', 'C4']]);

  await page.locator('[data-record-list-open="circular"]').click();
  await list.locator('[data-record-list-row="#2"]').uncheck();
  await list.getByRole('button', { name: 'Close', exact: true }).click();
  await page.evaluate(() => { window.__GBDRAW_APP__.form.multi_record_canvas = false; });
  await generateAndWaitForResult(page);
  // All records stays separate diagrams with one ON record (PD-OI-044).
  expect(await drawnRecordIds(page)).toEqual([['C4']]);
});

// Review F1: a record's OFF state belongs to the file it was read from.
test('a replaced Linear File and another Circular input draw every record', async ({ page }) => {
  test.setTimeout(180_000);
  await openLinear(page);
  await uploadLinear(page, 'old.gb', genbank(['O1', 'O2', 'O3']));
  await expect.poll(() => linearState(page)).toHaveLength(3);
  await openRecords(page);
  await (await drawCheckbox(page, 'O1')).uncheck();
  await uploadLinear(page, 'new.gb', genbank(['N1', 'N2']));
  await expect.poll(() => linearState(page)).toEqual([{ id: 'N1', drawn: true }, { id: 'N2', drawn: true }]);

  await page.evaluate((text) => {
    const app = window.__GBDRAW_APP__;
    app.setDiagramMode('circular');
    app.setCircularSourceFile('c_gb', new File([text], 'three.gb', { type: 'text/plain', lastModified: 1 }));
  }, genbank(['G1', 'G2', 'G3']));
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.circularRecordList.length)).toBe(3);
  await page.locator('[data-record-list-open="circular"]').click();
  await page.locator('[data-record-list-row="#2"]').uncheck();
  await page.locator('[data-record-list-close]').click();
  await page.getByRole('radio', { name: 'GFF3 + FASTA', exact: true }).first().check();
  await page.getByRole('radio', { name: 'GenBank', exact: true }).first().check();
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.circularRecordList.length)).toBe(3);
  expect(await page.evaluate(() => ['#1', '#2', '#3'].map((key) => window.__GBDRAW_APP__.recordSelection.isDrawn('circular', key))))
    .toEqual([true, true, true]);
});

// Review F4: a Legend rename of a row only an OFF record draws survives a
// Generate that replaced another File (the request carries the rename, so the
// row returns renamed when its record is ON again).
test('a Legend rename waits while its record is OFF, also across a source replacement', async ({ page }) => {
  test.setTimeout(300_000);
  await openLinear(page);
  await uploadLinear(page, 'with-trna.gb', genbank(['T1', { id: 'T2', trna: true }]));
  await uploadLinear(page, 'other.gb', genbank(['P1']), { newFile: true });
  await expect.poll(() => linearState(page)).toHaveLength(3);
  await generateAndWaitForResult(page);
  await evaluateWithRetainedPromise(page, async () => {
    const app = window.__GBDRAW_APP__;
    const index = app.legendEntries.findIndex((entry) => entry.caption === 'tRNA');
    if (index < 0) throw new Error(`no tRNA row: ${app.legendEntries.map((entry) => entry.caption).join(', ')}`);
    await app.renameLegendEntry(index, 'transfer RNA');
  });
  await openRecords(page);
  await (await drawCheckbox(page, 'T2')).uncheck();
  await generateAndWaitForResult(page);
  await uploadLinear(page, 'other-v2.gb', genbank(['P2']), { index: 2 });
  await generateAndWaitForResult(page);
  await (await drawCheckbox(page, 'T2')).check();
  await generateAndWaitForResult(page);
  const captions = await page.evaluate(() => window.__GBDRAW_APP__.legendEntries.map((entry) => entry.caption));
  expect(captions).toContain('transfer RNA');
  expect(captions).not.toContain('tRNA');
});
