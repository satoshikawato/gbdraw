const { test, expect } = require('@playwright/test');
const { openApp, generateAndWaitForResult } = require('./helpers/app-lifecycle.cjs');

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
    if (mode === 'linear') await page.getByRole('button', { name: 'Linear', exact: true }).click();
    const upload = mode === 'linear' ? page.getByTestId('linear-genbank-1') : page.getByLabel('GenBank/DDBJ File', { exact: true });
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
      const file = mode === 'linear' ? app.linearSeqs[0].gb : app.files.c_gb;
      const { discoverSequenceRecords, normalizeSequenceRecords } = await import('./js/app/record-discovery.js');
      const { runDiagramHelperOperation, DIAGRAM_HELPER_OPERATIONS } = await import('./js/services/diagram-generation.js');
      const fast = await discoverSequenceRecords({ file, format: 'genbank' });
      const response = await runDiagramHelperOperation(DIAGRAM_HELPER_OPERATIONS.LIST_SEQUENCE_RECORDS, {
        format: 'genbank', files: [{ role: 'source', bytes: await file.arrayBuffer() }]
      });
      return { fast, worker: normalizeSequenceRecords(response.result) };
    }, mode);
    expect(discovery.fast).toEqual(discovery.worker);
    expect(discovery.fast.map(r => [r.selector, r.recordId, r.detectedTopology])).toEqual([
      ['#1', 'same', 'circular'], ['#2', 'same', 'linear']
    ]);
    await generateAndWaitForResult(page);
    expect(await page.evaluate(() => window.__GBDRAW_APP__.results.length)).toBeGreaterThan(0);
    const shortcuts = page.getByRole('button', { name: 'Use selected feature midpoint', exact: true });
    await expect(shortcuts).toHaveCount(mode === 'circular' ? 2 : 1);
    for (const shortcut of await shortcuts.all()) await expect(shortcut).toBeDisabled();
    await expect(page.getByRole('spinbutton', { name: 'Display start same #1', exact: true })).toBeEnabled();
    expect(external).toEqual([]);
  });
}

const { readFileSync } = require('node:fs');
const { join } = require('node:path');
const {
  assertDiagramWorkerIdle, assertSessionLoadLeftWorkerIdle, getDiagramWorkerActivity
} = require('./helpers/app-lifecycle.cjs');
const fixture = (path) => join(process.cwd(), path);
const snapshot = (page) => page.evaluate(async () => {
  const { state } = await import('./js/state.js');
  return {
    status: state.circularRecordDiscovery.status,
    viewStatus: window.__GBDRAW_APP__.circularRecordDiscoveryState.status,
    rows: state.circularRecordList.value.map(row => [row.selector, row.record_id]),
    primary: state.files.c_gb?.name || state.files.c_gff?.name || null,
    result: state.results.value[0]?.content || null,
    selector: state.form.circular_record_selector,
    grouping: state.adv.circular_grouping_intent,
    grid: state.form.multi_record_canvas
  };
});
const waitStatus = (page, status) => expect.poll(() => page.evaluate(async () => (
  (await import('./js/state.js')).state.circularRecordDiscovery.status
)), { timeout: 180000 }).toBe(status);
const uploadNative = (page, text, name = 'native.gb') => page.getByLabel('GenBank/DDBJ File', { exact: true })
  .setInputFiles({ name, mimeType: 'text/plain', buffer: Buffer.from(text) });
const loadPreview = async (page) => {
  page.on('dialog', dialog => dialog.accept());
  await page.locator('input[accept^=".json,"]').setInputFiles(fixture(
    'gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json'));
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
    await expect(page.getByRole('button', { name: 'Load record rotation controls' })).toHaveCount(0);
    await expect(page.getByRole('spinbutton', { name: /^Display start/ })).toHaveCount(expected.length);
    await assertDiagramWorkerIdle(page);
  }
});

test('incomplete GFF pair stays idle; complete pair automatically discovers without Python', async ({ page }) => {
  await openApp(page);
  await page.locator('input[type=radio][value=gff]').first().check();
  await page.getByLabel('GFF3 File', { exact: true }).setInputFiles(fixture('tests/test_inputs/NC_013668.gff3'));
  await waitStatus(page, 'idle');
  await expect(page.locator('[data-circular-discovery-status]')).toContainText('Upload both');
  await page.getByLabel('FASTA File', { exact: true }).setInputFiles(fixture('tests/test_inputs/NC_013668.fasta'));
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
  await page.getByRole('button', { name: 'Linear', exact: true }).click();
  await page.getByRole('button', { name: 'Circular', exact: true }).click();
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
  await page.getByRole('group', { name: 'GenBank/DDBJ File selection', exact: true }).locator('button').click();
  await waitStatus(page, 'idle');
  expect((await snapshot(page)).rows).toEqual([]);
  expect((await snapshot(page)).result).toBe(original.result);
});

test('rapid replace/remove, mode and input-type changes reject delayed native completion', async ({ page }) => {
  await openApp(page);
  const delay = async () => {
    await page.evaluate(async text => {
      const { state } = await import('./js/state.js');
      const source = new File([text], 'slow.gb');
      const pending = new Promise(resolve => { window.__releaseSlow = () => resolve(new TextEncoder().encode(text).buffer); });
      source.arrayBuffer = () => pending;
      state.files.c_gb = source;
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
  await page.evaluate(() => { window.__GBDRAW_APP__.files.c_gb = null; window.__releaseSlow(); });
  await waitStatus(page, 'idle');
  expect((await snapshot(page)).rows).toEqual([]);
  await delay();
  await page.getByRole('button', { name: 'Linear', exact: true }).click();
  await page.evaluate(() => window.__releaseSlow());
  await waitStatus(page, 'idle');
  await page.getByRole('button', { name: 'Circular', exact: true }).click();
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
  const helperSource = async () => {
    await page.evaluate(async text => {
      const { state } = await import('./js/state.js');
      const source = new File([text], 'helper.gb');
      const read = source.arrayBuffer.bind(source);
      let first = true;
      source.arrayBuffer = () => {
        if (first) { first = false; return Promise.reject(new Error('forced lightweight reader miss')); }
        return read();
      };
      window.__holdHelper = true;
      state.files.c_gb = source;
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
