const { test, expect } = require('@playwright/test');
const { evaluateWithRetainedPromise, openApp, generateAndWaitForResult } = require('./helpers/app-lifecycle.cjs');
const { showDrawing } = require('./helpers/drawing-bar.cjs');

const genbank = (topology, length = 120) => `LOCUS       same                     ${length} bp    DNA     ${topology.padEnd(8)} UNA 01-JAN-2000
DEFINITION  Source discovery regression fixture.
ACCESSION   same
VERSION     same
KEYWORDS    .
SOURCE      synthetic construct
  ORGANISM  synthetic construct
FEATURES             Location/Qualifiers
     CDS             11..35
                     /product="source feature"
ORIGIN
        1 ${'acgt'.repeat(length / 4)}
//
`;

for (const mode of ['circular', 'linear']) {
  test(`record topology discovery uses uploaded source and packaged Worker in ${mode}`, async ({ page }) => {
    test.setTimeout(180000);
    const external = [];
    await page.route('**/*', route => {
      if (new URL(route.request().url()).hostname === '127.0.0.1') return route.continue();
      external.push(route.request().url());
      return route.abort();
    });
    await openApp(page);
    if (mode === 'linear') await showDrawing(page, { kind: 'linear' });
    const upload = mode === 'linear' ? page.getByTestId('linear-genbank-1') : circularFileInput(page, 'gb');
    await upload.setInputFiles({ name: 'same.gbk', mimeType: 'text/plain', buffer: Buffer.from(genbank('circular') + genbank('linear')) });
    if (mode === 'linear') {
      await page.getByRole('button', { name: 'Records for file 1', exact: true }).click();
      await page.getByRole('button', { name: 'Record options for sequence 1', exact: true }).click();
      const selector = page.getByRole('combobox', { name: 'Record selector for sequence 1', exact: true });
      await expect(selector).toBeEnabled();
      await selector.selectOption('#1');
    }
    const discovery = await page.evaluate(async mode => {
      const app = window.__GBDRAW_APP__;
      const { state } = await import('./js/state.js');
      const file = state.inputFiles.fileOf(mode === 'linear' ? app.linearSeqs[0].gb : state.activeDrawing().circularSources[0].gb);
      const { discoverSequenceRecords, normalizeSequenceRecords } = await import('./js/services/record-discovery.js');
      const { runDiagramHelperOperation, DIAGRAM_HELPER_OPERATIONS } = await import('./js/services/diagram-generation.js');
      const fast = await discoverSequenceRecords({ file, format: 'genbank' });
      const response = await runDiagramHelperOperation(DIAGRAM_HELPER_OPERATIONS.LIST_SEQUENCE_RECORDS, {
        format: 'genbank', files: [{ role: 'source', blob: file }]
      });
      return { fast, worker: normalizeSequenceRecords(response.result) };
    }, mode);
    expect(discovery.fast).toEqual(discovery.worker);
    expect(discovery.fast.map(r => [r.selector, r.recordId, r.detectedTopology])).toEqual([
      ['#1', 'same', 'circular'], ['#2', 'same', 'linear']
    ]);
    await generateAndWaitForResult(page);
    expect(await page.evaluate(() => window.__GBDRAW_APP__.results.length)).toBeGreaterThan(0);
    // Each drawn record has its Display start; Reset start shows only beside a requested start.
    await expect(page.getByRole('spinbutton', { name: /^Display start/ })).toHaveCount(mode === 'circular' ? 2 : 1);
    await expect(page.getByRole('button', { name: 'Reset start', exact: true })).toHaveCount(0);
    // OV-212 (PD-OI-087): the feature-based start is only in the feature popup.
    await expect(page.getByRole('button', { name: /^Use selected feature/ })).toHaveCount(0);
    await expect(page.getByRole('spinbutton', { name: 'Display start same #1', exact: true })).toBeEnabled();
    expect(external).toEqual([]);
  });
}

const { readFileSync } = require('node:fs');
const { readSessionFile } = require('./helpers/session-file.cjs');
const { join } = require('node:path');
const {
  assertDiagramWorkerIdle, assertSessionLoadLeftWorkerIdle, getDiagramWorkerActivity
} = require('./helpers/app-lifecycle.cjs');
const { circularFileInput } = require('./helpers/circular-inputs.cjs');
const fixture = (path) => join(process.cwd(), path);
const snapshot = (page) => page.evaluate(async () => {
  const { state } = await import('./js/state.js');
  const circular = state.drawingList().find((drawing) => drawing.mode === 'circular');
  const fileName = (ref) => (ref ? state.inputFiles.fileOf(ref).name : null);
  return {
    status: circular.circularRecordDiscovery.status,
    viewStatus: window.__GBDRAW_APP__.circularRecordDiscoveryState.status,
    rows: circular.circularRecordList.value.map(row => [row.selector, row.record_id]),
    primary: fileName(circular.circularSources[0]?.gb) || fileName(circular.circularSources[0]?.gff) || null,
    result: window.__GBDRAW_APP__.currentResults()[0]?.content || null,
    selector: state.activeDrawing().form.circular_record_selector,
    grouping: state.activeDrawing().adv.circular_grouping_intent,
    grid: state.activeDrawing().form.multi_record_canvas
  };
});
const waitStatus = (page, status) => expect.poll(() => page.evaluate(async () => (
  (await import('./js/state.js')).state.drawingList().find((drawing) => drawing.mode === 'circular').circularRecordDiscovery.status
)), { timeout: 180000 }).toBe(status);
const uploadNative = (page, text, name = 'native.gb') => circularFileInput(page, 'gb')
  .setInputFiles({ name, mimeType: 'text/plain', buffer: Buffer.from(text) });
// The GenBank file of the Circular drawing's one File card, set as its card's
// slot sets it: a new card when there is none; null removes the card.
const defineSetCardGenbank = (page) => page.evaluate(() => {
  window.setCardGenbank = async (file) => {
    const { state } = await import('./js/state.js');
    const app = window.__GBDRAW_APP__;
    const card = state.drawingList().find((drawing) => drawing.mode === 'circular').circularSources[0];
    return card ? app.setCircularSourceFile(card.uid, 'gb', file) : app.addCircularSourceFiles([file]);
  };
});
const loadPreview = async (page) => {
  page.on('dialog', dialog => dialog.accept());
  await page.locator('input[accept*=".gbdraw"]').setInputFiles(fixture(
    'gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw'));
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending
    && window.__GBDRAW_APP__.results.length === 1);
};

test('native upload automatically inspects one, two, duplicate and DDBJ records before Generate', async ({ page }) => {
  await openApp(page);
  for (const [name, text, expected] of [
    ['one.gb', genbank('circular'), [['#1', 'same']]],
    ['two.gb', genbank('circular') + genbank('linear').replaceAll('same', 'second'), [['#1', 'same'], ['#2', 'second']]],
    ['duplicate.gb', genbank('circular') + genbank('linear'), [['#1', 'same'], ['#2', 'same']]],
    ['MjeNMV.ddbj', readFileSync(fixture('tests/test_inputs/MjeNMV.gbk'), 'utf8'), [['#1', 'LC738868.1']]]
  ]) {
    await uploadNative(page, text, name);
    await waitStatus(page, 'ready');
    expect((await snapshot(page)).rows).toEqual(expected);
    await expect(page.getByRole('button', { name: 'Inspect records' })).toHaveCount(0);
    await expect(page.getByRole('spinbutton', { name: /^Display start/ })).toHaveCount(expected.length);
    // The Gallery capture tool awaits this wait after every record read (tools/capture_gallery_tutorial_screenshots.py).
    expect(await page.evaluate(() => window.__GBDRAW_APP__.waitForCircularRecordScrollAnchor())).toBeUndefined();
    await assertDiagramWorkerIdle(page);
  }
});

test('incomplete GFF pair stays idle; complete pair automatically discovers without Python', async ({ page }) => {
  await openApp(page);
  await page.locator('input[type=radio][value=gff]').first().check();
  await circularFileInput(page, 'gff').setInputFiles(fixture('tests/test_inputs/NC_013668.gff3'));
  await waitStatus(page, 'idle');
  await expect(page.locator('[data-circular-discovery-status]')).toHaveText('');
  await circularFileInput(page, 'fasta').setInputFiles(fixture('tests/test_inputs/NC_013668.fasta'));
  await waitStatus(page, 'ready');
  expect((await snapshot(page)).rows).toEqual([['#1', 'NC_013668.3']]);
  await assertDiagramWorkerIdle(page);
});

test('saved preview remains deferred across disclosure and mode changes; Inspect settles without catalog reuse', async ({ page }) => {
  await openApp(page);
  await loadPreview(page);
  const original = await snapshot(page);
  expect(original.status).toBe('deferred');
  expect(original.rows).toEqual([]);
  await expect(page.locator('[data-circular-discovery-status]')).toHaveText('Records not inspected');
  const summary = page.getByLabel('Circular record presentation', { exact: true });
  await summary.click();
  await summary.click();
  await showDrawing(page, { kind: 'linear' });
  await showDrawing(page, { kind: 'circular' });
  expect((await snapshot(page)).status).toBe('deferred');
  await assertSessionLoadLeftWorkerIdle(page);
  const inspectionHistory = await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
  const inspect = page.locator('[data-circular-inspect]');
  await inspect.click();
  await waitStatus(page, 'ready');
  await expect(inspect).toBeFocused();
  expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(inspectionHistory);
  expect((await snapshot(page)).rows).toEqual([['#1', 'NC_012920.1']]);
  expect((await snapshot(page)).result).toBe(original.result);
  await assertDiagramWorkerIdle(page);
});

test('saved preview Generate inspects the source before rendering', async ({ page }) => {
  test.setTimeout(180000);
  await openApp(page);
  await loadPreview(page);
  await assertSessionLoadLeftWorkerIdle(page);
  await generateAndWaitForResult(page);
  expect((await snapshot(page)).status).toBe('ready');
  expect((await snapshot(page)).rows).toEqual([['#1', 'NC_012920.1']]);
  expect((await getDiagramWorkerActivity(page)).runs).toBe(1);
});

test('replacement draft does not inherit saved catalog; invalid, Retry, Replace and Remove keep the Result', async ({ page }) => {
  test.setTimeout(180000);
  await openApp(page);
  await loadPreview(page);
  const original = await snapshot(page);
  await uploadNative(page, genbank('circular') + genbank('linear'), 'replacement.gb');
  await waitStatus(page, 'ready');
  expect((await snapshot(page)).rows).toEqual([['#1', 'same'], ['#2', 'same']]);
  expect((await snapshot(page)).result).toBe(original.result);
  await assertDiagramWorkerIdle(page);
  await uploadNative(page, 'invalid source', 'invalid.gb');
  await waitStatus(page, 'error');
  expect((await snapshot(page)).primary).toBe('invalid.gb');
  expect((await snapshot(page)).result).toBe(original.result);
  await expect(page.locator('[data-circular-discovery-status]')).toHaveText(
    'No records were found. Choose input containing records. The file is not a GenBank/DDBJ flat file: it has no record header line.'
  );
  await expect(page.locator('[data-circular-discovery-status]')).not.toContainText('invalid.gb');
  await page.getByRole('button', { name: 'Retry source inspection' }).click();
  await waitStatus(page, 'error');
  expect((await snapshot(page)).result).toBe(original.result);
  await uploadNative(page, genbank('circular'));
  await waitStatus(page, 'ready');
  // A file leaves with its File card, through the removal dialog.
  await page.getByRole('button', { name: 'Remove File 1', exact: true }).click();
  await page.locator('[data-circular-source-removal-dialog] [role=dialog]').getByRole('button', { name: 'Remove file' }).click();
  await waitStatus(page, 'idle');
  expect((await snapshot(page)).rows).toEqual([]);
  expect((await snapshot(page)).result).toBe(original.result);
});

test('rapid replace/remove, mode and input-type changes reject delayed native completion', async ({ page }) => {
  await openApp(page);
  await defineSetCardGenbank(page);
  const delay = async () => {
    await page.evaluate(async text => {
      const source = new File([text], 'slow.gb');
      const pending = new Promise(resolve => { window.__releaseSlow = () => resolve(new TextEncoder().encode(text).buffer); });
      source.arrayBuffer = () => pending;
      await setCardGenbank(source);
    }, genbank('circular'));
    await waitStatus(page, 'loading');
    await expect(page.getByLabel('Circular record', { exact: true })).toBeDisabled();
  };
  await delay();
  await uploadNative(page, genbank('linear').replaceAll('same', 'new'));
  await waitStatus(page, 'ready');
  await page.evaluate(() => window.__releaseSlow());
  expect((await snapshot(page)).rows).toEqual([['#1', 'new']]);
  await delay();
  await page.evaluate(async () => { await setCardGenbank(null); window.__releaseSlow(); });
  await waitStatus(page, 'idle');
  expect((await snapshot(page)).rows).toEqual([]);
  await delay();
  // Each drawing owns its discovery and a drawing switch runs none: the
  // Circular source is still current while Linear is shown, so its reply
  // settles the Circular drawing and is shown on return.
  await showDrawing(page, { kind: 'linear' });
  await page.evaluate(() => window.__releaseSlow());
  await waitStatus(page, 'ready');
  expect((await snapshot(page)).rows).toEqual([['#1', 'same']]);
  await showDrawing(page, { kind: 'circular' });
  await waitStatus(page, 'ready');
  await delay();
  await page.locator('input[type=radio][value=gff]').first().check();
  await page.evaluate(() => window.__releaseSlow());
  await waitStatus(page, 'idle');
  expect((await snapshot(page)).rows).toEqual([]);
  await assertDiagramWorkerIdle(page);
});

test('existing helper stays loading until settlement and rejects superseded source replies', async ({ page }) => {
  test.setTimeout(180000);
  await page.addInitScript(() => {
    const add = Worker.prototype.addEventListener;
    Worker.prototype.addEventListener = function (type, listener, options) {
      return add.call(this, type, type !== 'message' ? listener : function (event) {
        if (event.data?.type === 'helper' && window.__holdHelper) {
          window.__heldHelpers.push(() => listener.call(this, event));
        } else listener.call(this, event);
      }, options);
    };
    window.__heldHelpers = [];
  });
  await openApp(page);
  await defineSetCardGenbank(page);
  const helperSource = async () => {
    await page.evaluate(async text => {
      const source = new File([text], 'helper.gb');
      const read = source.arrayBuffer.bind(source);
      let first = true;
      source.arrayBuffer = () => {
        if (first) { first = false; return Promise.reject(new Error('forced lightweight reader miss')); }
        return read();
      };
      window.__holdHelper = true;
      await setCardGenbank(source);
    }, genbank('circular'));
    await waitStatus(page, 'loading');
    await expect(page.locator('[data-circular-inspect]')).toHaveAttribute('aria-disabled', 'true');
    expect(await page.evaluate(() => window.__GBDRAW_APP__.inspectCircularSourceRecords())).toMatchObject({ status: 'unavailable' });
    await page.waitForFunction(() => window.__heldHelpers.length > 0, null, { timeout: 180000 });
    expect((await snapshot(page)).status).toBe('loading');
  };
  const settle = () => page.evaluate(() => {
    window.__holdHelper = false;
    window.__heldHelpers.splice(0).forEach(release => release());
  });
  await helperSource();
  await settle();
  await waitStatus(page, 'ready');
  expect((await snapshot(page)).rows).toEqual([['#1', 'same']]);
  await helperSource();
  await uploadNative(page, genbank('linear').replaceAll('same', 'replacement'));
  await waitStatus(page, 'ready');
  await settle();
  expect((await snapshot(page)).rows).toEqual([['#1', 'replacement']]);
  const worker = await getDiagramWorkerActivity(page);
  expect(worker.constructions).toBe(1);
  expect(worker.helpers).toBe(2);
  expect(worker.settledHelpers).toBe(2);
});

// A FASTA given without GFF3 (PD-OI-068 rev 2): Python's FASTA-only loader lists its records in the Worker,
// each with unknown topology, and Generate draws them without features.
const twoRecordFasta = `>alpha first\n${'acgt'.repeat(60)}\n>beta\n${'ggcc'.repeat(45)}\n`;
const drawnText = (page) => page.evaluate(() => window.__GBDRAW_APP__.currentResults()
  .map((result) => String(result.content || '')).join('\n'));
for (const mode of ['circular', 'linear']) {
  test(`a FASTA alone lists its records and Generate draws them in ${mode}`, async ({ page }) => {
    test.setTimeout(240000);
    await openApp(page);
    if (mode === 'linear') await showDrawing(page, { kind: 'linear' });
    await page.locator('input[type=radio][value=fasta]').first().check();
    await expect(page.getByText('Records are drawn without features.')).toBeVisible();
    await expect(page.locator('input[type=radio][value=gff]').first()).toBeVisible();
    const upload = mode === 'linear' ? page.getByTestId('linear-fasta-1') : circularFileInput(page, 'fasta');
    await upload.setInputFiles({ name: 'two.fasta', mimeType: 'text/plain', buffer: Buffer.from(twoRecordFasta) });
    if (mode === 'circular') {
      await waitStatus(page, 'ready');
      expect((await snapshot(page)).rows).toEqual([['#1', 'alpha'], ['#2', 'beta']]);
      // Each record's Circular record row starts off: the file makes no topology claim.
      const circularRows = page.getByRole('checkbox', { name: /^Circular record/ });
      await expect(circularRows).toHaveCount(2);
      for (const row of await circularRows.all()) await expect(row).not.toBeChecked();
    } else {
      await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.linearSeqs.length), { timeout: 180000 }).toBe(2);
      await expect(page.getByLabel('Default definition for file 1', { exact: true })).toBeVisible();
    }
    await generateAndWaitForResult(page);
    const drawn = await drawnText(page);
    expect(drawn).toContain('alpha');
    expect(drawn).toContain('beta');
  });
}

test('a protein FASTA alone is refused before Generate with the record named', async ({ page }) => {
  test.setTimeout(240000);
  await openApp(page);
  await page.locator('input[type=radio][value=fasta]').first().check();
  await circularFileInput(page, 'fasta').setInputFiles({
    name: 'protein.fasta', mimeType: 'text/plain', buffer: Buffer.from('>LvMJNV_ORF1 hypothetical protein\nMKLVEQRSTWYF\n')
  });
  await waitStatus(page, 'error');
  const status = page.locator('[data-circular-discovery-status]');
  await expect(status).toContainText('LvMJNV_ORF1');
  await expect(status).toContainText('"L"');
  await expect(page.locator('.upload-zone.failed')).toHaveCount(1);
  await generateAndWaitForResult(page, { expectedStatus: 'error' });
  expect(await page.evaluate(() => window.__GBDRAW_APP__.results.length)).toBe(0);
});

// In Linear the refusal is in the card's Record options. A refused FASTA file leaves Generate as a refused GenBank
// file of the same card does (D-21): the same protein file, given as each input type.
for (const [inputType, testId, named] of [
  ['gb', 'linear-genbank-1', /looks like FASTA/],
  ['fasta', 'linear-fasta-1', /LvMJNV_ORF1[\s\S]*"L"/]
]) {
  test(`a protein file refused as ${inputType} in linear names the failure in Record options and Generate reports it`, async ({ page }) => {
    test.setTimeout(240000);
    await openApp(page);
    await showDrawing(page, { kind: 'linear' });
    await page.locator(`input[type=radio][value=${inputType}]`).first().check();
    await page.getByTestId(testId).setInputFiles({
      name: 'protein.fasta', mimeType: 'text/plain', buffer: Buffer.from('>LvMJNV_ORF1 hypothetical protein\nMKLVEQRSTWYF\n')
    });
    await expect(page.locator('p.text-red-600', { hasText: named })).toHaveCount(1, { timeout: 180000 });
    await expect(page.locator('.upload-zone.failed')).toHaveCount(1);
    await generateAndWaitForResult(page, { expectedStatus: 'error' });
    expect(await page.evaluate(() => window.__GBDRAW_APP__.results.length)).toBe(0);
  });
}

test('a saved Session with a FASTA source loads its Result without the Worker', async ({ page, context }) => {
  test.setTimeout(300000);
  await openApp(page);
  await page.locator('input[type=radio][value=fasta]').first().check();
  await circularFileInput(page, 'fasta').setInputFiles(fixture('tests/test_inputs/NC_013668.fasta'));
  await waitStatus(page, 'ready');
  await generateAndWaitForResult(page);
  const download = page.waitForEvent('download');
  await evaluateWithRetainedPromise(page, async () => {
    window.__GBDRAW_APP__.sessionTitle = 'fasta-source';
    await window.__GBDRAW_APP__.saveSessionWithTitle();
  });
  const saved = await (await download).path();
  const session = await readSessionFile(saved);
  const drawing = session.drawings.find((entry) => entry.id === session.ui.activeDrawingId);
  expect(drawing.renderRequest.records[0].source.kind).toBe('fasta');

  const fresh = await context.newPage();
  await openApp(fresh);
  fresh.on('dialog', dialog => dialog.accept());
  await fresh.locator('input[accept*=".gbdraw"]').setInputFiles(saved);
  await fresh.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending
    && window.__GBDRAW_APP__.results.length === 1, null, { timeout: 180000 });
  await assertSessionLoadLeftWorkerIdle(fresh);
  await expect(fresh.locator('input[type=radio][value=fasta]').first()).toBeChecked();
  expect((await snapshot(fresh)).status).toBe('deferred');
  expect(await drawnText(fresh)).toContain('NC_013668.3');
});


// A Linear card that held a GFF3 + FASTA pair and was switched to FASTA holds nothing once its FASTA file is
// removed: the call site asks the input type which files a card holds, so the GFF3 file left in the other
// type's slot keeps no card.
test('removing the FASTA file of a Linear card switched from GFF3 + FASTA to FASTA empties the card', async ({ page }) => {
  test.setTimeout(180000);
  await openApp(page);
  await showDrawing(page, { kind: 'linear' });
  await page.locator('input[type=radio][value=gff]').first().check();
  const card = page.locator('[data-linear-source-card]').first();
  await card.locator('input[type=file][accept^=".gff"]').setInputFiles(fixture('tests/test_inputs/NC_013668.gff3'));
  await card.locator('input[type=file][accept^=".fa"]').setInputFiles(fixture('tests/test_inputs/NC_013668.fasta'));
  await expect.poll(() => page.evaluate(() => {
    const [seq] = window.__GBDRAW_APP__.linearSeqs;
    return Boolean(seq.gff && seq.fasta);
  })).toBe(true);
  await page.locator('input[type=radio][value=fasta]').first().check();
  await page.evaluate(() => window.__GBDRAW_APP__.setLinearSeqPrimaryFile(0, 'fasta', null));
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.linearSeqs.map(
    (seq) => ({ gff: seq.gff?.name || null, fasta: seq.fasta?.name || null })
  ))).toEqual([{ gff: null, fasta: null }]);
});
