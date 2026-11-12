const { test, expect } = require('@playwright/test');
const { execFile } = require('node:child_process');
const { promisify } = require('node:util');
const path = require('node:path');
const { openApp, generateAndWaitForResult } = require('./helpers/app-lifecycle.cjs');

// A Circular grid or batch Session that the command line saved from a records
// table keeps each crop as its record's range on the unchanged input file. The
// Web loads the range per record, draws it while the drawing shows every
// record, and puts it in the crop fields when Record select chooses that record
// alone. The pure parts are in tests/web/circular-record-ranges.test.mjs.
const FIXTURES = path.resolve('tests/fixtures');
const twoRecords = path.join(FIXTURES, 'web_batch_two_records.gb');
const mito = path.join(FIXTURES, 'sessions', 'cli-web-mito.gb');

const cliSession = async (testInfo, rows, extra = []) => {
  const table = testInfo.outputPath('records.tsv');
  require('node:fs').writeFileSync(table, `gbk\trecord_id\tregion\n${rows.map((row) => row.join('\t')).join('\n')}\n`);
  const output = testInfo.outputPath('cli.gbdraw');
  await promisify(execFile)('python', ['-m', 'gbdraw.cli', 'circular', '--records_table', table, ...extra,
    '-o', testInfo.outputPath('cli'), '-f', 'svg', '--session_output', output],
  { cwd: testInfo.outputDir, env: { ...process.env, PYTHONPATH: process.cwd() }, timeout: 300_000 });
  return output;
};
const loadSession = async (page, file) => {
  page.on('dialog', (dialog) => dialog.accept());
  await openApp(page);
  await page.locator('input[accept*=".gbdraw"]').setInputFiles(file);
  await expect.poll(() => page.evaluate(() => !window.__GBDRAW_APP__.sessionImportPending), { timeout: 180_000 }).toBe(true);
};
const committedWindows = (page) => page.evaluate(async () => {
  const request = (await import('/gbdraw/web/js/services/config.js')).getCommittedCanonicalRenderRequest();
  return request.records.map(({ windows }) => windows.map(({ start, end, reverseComplement }) => ({ start, end, reverseComplement })));
});
const cropFields = (page) => page.evaluate(() => {
  const { circular_region_start: start, circular_region_end: end, circular_reverse: reverse } = window.__GBDRAW_APP__.form;
  return { start, end, reverse };
});
const chooseRecord = async (page, recordId) => {
  const select = page.getByLabel('Circular record', { exact: true });
  const value = () => select.locator('option').evaluateAll((options, key) => {
    const found = options.find((option) => option.value === '' ? key === '' : option.textContent.includes(key));
    return found ? found.value : null;
  }, recordId);
  await expect.poll(value).not.toBeNull();
  await select.selectOption(await value());
};

test('a CLI batch of one crop loads as its range on the input file and Generate draws the crop', async ({ page }, testInfo) => {
  test.setTimeout(600_000);
  await page.setViewportSize({ width: 1280, height: 1000 });
  await loadSession(page, await cliSession(testInfo, [[mito, '', '1001-3000']]));
  await generateAndWaitForResult(page);
  expect(await committedWindows(page)).toEqual([[{ start: 1001, end: 3000, reverseComplement: false }]]);
  const drawn = await page.evaluate(() => window.__GBDRAW_APP__.currentResults().map((result) => String(result.content || '')));
  expect(drawn).toHaveLength(1);
  expect(drawn[0]).toMatch(/>2,000 bp<\/text>/);
  await page.screenshot({ path: testInfo.outputPath('batch-of-one-after.png') });
});

test('Record select puts a batch record\'s loaded range in the crop fields, and all records draw the loaded ranges', async ({ page }, testInfo) => {
  test.setTimeout(600_000);
  await loadSession(page, await cliSession(testInfo, [
    [twoRecords, 'TESTA', '101-500:rc'], [twoRecords, 'TESTB', ''], [mito, '', '1001-3000']
  ]));
  const loaded = [[{ start: 101, end: 500, reverseComplement: true }], [], [{ start: 1001, end: 3000, reverseComplement: false }]];
  // The ranges are the drawing's, not the canvas's: Multi-Record Canvas on draws them on one canvas, and off again
  // one figure per record.
  for (const canvas of [true, false]) {
    await page.evaluate((on) => { window.__GBDRAW_APP__.form.multi_record_canvas = on; }, canvas);
    await generateAndWaitForResult(page);
    expect(await committedWindows(page)).toEqual(loaded);
  }
  await chooseRecord(page, 'NC_012920');
  expect(await cropFields(page)).toEqual({ start: 1001, end: 3000, reverse: false });
  await chooseRecord(page, 'TESTA');
  expect(await cropFields(page)).toEqual({ start: 101, end: 500, reverse: true });
  await chooseRecord(page, 'TESTB');
  expect(await cropFields(page)).toEqual({ start: null, end: null, reverse: false });
  await chooseRecord(page, '');
  await generateAndWaitForResult(page);
  expect(await committedWindows(page)).toEqual(loaded);
});
