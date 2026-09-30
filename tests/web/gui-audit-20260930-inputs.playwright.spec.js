// Web GUI audit 2026-09-30: inputs, record topology, Session Save/Load and
// Reset. A test marked test.fail(true, '<ID>') asserts the correct behavior of
// a current defect; the PR that fixes the audit ID removes the mark.
const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { generateAndWaitForResult } = require('./helpers/app-lifecycle.cjs');
const {
  BATCH_FIXTURE,
  HMMT,
  openFresh,
  openWithGenBank,
  settle
} = require('./helpers/audit-browser.cjs');

test.describe.configure({ retries: 0 });

const BATCH_TEXT = readFileSync(BATCH_FIXTURE, 'utf8');
const FIRST_RECORD_TEXT = `${BATCH_TEXT.split(/^\/\/\s*$/m)[0].trimEnd()}\n//\n`;

const openLinear = async (page) => {
  await openFresh(page);
  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    app.setMode('linear');
    Object.assign(app.form, { legend: 'none', show_gc: false, show_skew: false });
  });
  await settle(page);
};

const linearCard = (page, index) => page.locator('[data-linear-source-card]').nth(index);

const uploadLinear = async (page, index, name, text) => {
  await linearCard(page, index).locator('input[type=file]').first()
    .setInputFiles({ name, mimeType: 'text/plain', buffer: Buffer.from(text) });
  await settle(page);
};

const recordDefinitions = (page) => page.evaluate(() => {
  const content = window.__GBDRAW_APP__.results[0]?.content || '';
  const document = new DOMParser().parseFromString(content, 'image/svg+xml');
  return Object.fromEntries([...document.querySelectorAll('g[data-gbdraw-role="record-definition"]')]
    .map((group) => [
      group.getAttribute('data-gbdraw-record-id'),
      [...group.querySelectorAll('text')].map((text) => text.textContent).join(' / ')
    ]));
});

const attemptSave = async (page) => {
  const downloads = [];
  const onDownload = (download) => downloads.push(download.suggestedFilename());
  page.on('download', onDownload);
  await page.evaluate(() => { window.__GBDRAW_APP__.sessionTitle = 'audit'; });
  await page.getByRole('banner').getByRole('button', { name: 'Save Session', exact: true }).click();
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionSavePending, null, { timeout: 120_000 });
  await settle(page);
  page.off('download', onDownload);
  return {
    downloads,
    error: await page.evaluate(() => {
      const error = window.__GBDRAW_APP__.errorLog;
      return error ? { code: error.code, summary: error.summary, actions: [...(error.actions || [])] } : null;
    })
  };
};

test('a Circular definition edit after Generate leaves the committed Result unchanged', async ({ page }) => {
  test.fail(true, 'IN-01');
  test.setTimeout(300_000);
  await openWithGenBank(page, HMMT, () => {
    Object.assign(window.__GBDRAW_APP__.form, {
      multi_record_canvas: false, circular_region_start: 1000, circular_region_end: 9000
    });
  });
  await generateAndWaitForResult(page);
  const committed = await page.evaluate(() => window.__GBDRAW_APP__.results[0].content);
  expect(committed).toContain('8,001 bp');
  await page.evaluate(() => { window.__GBDRAW_APP__.form.species = 'Homo sapiens'; });
  const rewritten = await page.waitForFunction(
    (original) => window.__GBDRAW_APP__.results[0].content !== original,
    committed,
    { timeout: 15_000 }
  ).then(() => true, () => false);
  expect(rewritten).toBe(false);
});

test('replacing a cropped Linear file with a multi-record file expands it and clears record options', async ({ page }) => {
  test.fail(true, 'IN-04');
  test.setTimeout(300_000);
  await openLinear(page);
  await uploadLinear(page, 0, 'first.gb', FIRST_RECORD_TEXT);
  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const sequence = app.linearSeqs[0];
    app.setLinearRecordCrop(sequence, 'region_start', 100);
    app.setLinearRecordCrop(sequence, 'region_end', 3000);
    sequence.definition = 'Genome A custom';
  });
  await generateAndWaitForResult(page);
  await uploadLinear(page, 0, 'two.gb', BATCH_TEXT);
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.linearSeqs.length), { timeout: 30_000 }).toBe(2);
  const rows = await page.evaluate(() => window.__GBDRAW_APP__.linearSeqs.map((sequence) => ({
    record: sequence.region_record_id,
    start: sequence.region_start ?? null,
    end: sequence.region_end ?? null,
    definition: sequence.definition || ''
  })));
  expect(rows.map(({ start, end, definition }) => ({ start, end, definition })))
    .toEqual([{ start: null, end: null, definition: '' }, { start: null, end: null, definition: '' }]);
  await generateAndWaitForResult(page);
});

// The shared record-discovery vector (G-D) supplies the loader's per-record
// inferred definitions for one file whose records name different organisms.
const TWO_ORGANISMS = JSON.parse(readFileSync('tests/fixtures/record_metadata_inference_cases.json', 'utf8'))
  .discovery.find(({ name }) => name === 'records of one file with different organisms');

test('each record of a multi-record Linear file keeps its own organism', async ({ page }) => {
  test.fail(true, 'IN-06');
  test.setTimeout(300_000);
  await openLinear(page);
  const source = TWO_ORGANISMS.files.source;
  await uploadLinear(page, 0, source.name, `${source.lines.join('\n')}\n`);
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.linearSeqs.length), { timeout: 30_000 }).toBe(2);
  await page.evaluate(() => window.__GBDRAW_APP__.setLinearRecordLayoutEnabled(false));
  await generateAndWaitForResult(page);
  const definitions = await recordDefinitions(page);
  for (const record of TWO_ORGANISMS.expected.records) {
    expect(definitions[record.recordId]).toContain(record.inferredDefinition.replace(/<\/?i>/g, ''));
  }
});

test('Save before the first Generate with an empty Linear card explains the missing input', async ({ page }) => {
  test.fail(true, 'IN-08');
  test.setTimeout(300_000);
  await openLinear(page);
  await uploadLinear(page, 0, 'first.gb', FIRST_RECORD_TEXT);
  await page.evaluate(() => window.__GBDRAW_APP__.addLinearSeq());
  await settle(page);
  const outcome = await attemptSave(page);
  expect(outcome.downloads).toEqual([]);
  expect(outcome.error?.code).toBeTruthy();
  expect(outcome.error.code).not.toBe('UNKNOWN');
});

test('a failed Save does not offer Save Session again', async ({ page }) => {
  test.fail(true, 'N-12');
  test.setTimeout(300_000);
  await openLinear(page);
  await uploadLinear(page, 0, 'first.gb', FIRST_RECORD_TEXT);
  await page.evaluate(() => window.__GBDRAW_APP__.addLinearSeq());
  await settle(page);
  const outcome = await attemptSave(page);
  expect(outcome.downloads).toEqual([]);
  const alert = page.locator('.border-l-red-500');
  await expect(alert).toBeVisible();
  await expect(alert.getByRole('button', { name: 'Save Session', exact: true })).toHaveCount(0);
});

test('Save after loading a v39 Session asks for Generate instead of failing with UNKNOWN', async ({ page }) => {
  test.fail(true, 'SE-05');
  test.setTimeout(300_000);
  await openFresh(page);
  await page.locator('input[accept^=".json,"]')
    .setInputFiles('tests/fixtures/sessions/BGC0000708-BGC0000713.v39.gbdraw-session.json.gz');
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending, null, { timeout: 180_000 });
  await page.waitForFunction(() => window.__GBDRAW_APP__.sessionSaveAvailable, null, { timeout: 180_000 });
  await settle(page);
  const outcome = await attemptSave(page);
  expect(outcome.downloads).toEqual([]);
  expect(outcome.error?.code).toBeTruthy();
  expect(outcome.error.code).not.toBe('UNKNOWN');
  expect(outcome.error.summary).toMatch(/Generate/);
});

test('loading a JSON file that is not a Session reports a recognized diagnostic', async ({ page }) => {
  test.fail(true, 'X-01');
  test.setTimeout(300_000);
  await openFresh(page);
  await page.locator('input[accept^=".json,"]').setInputFiles({
    name: 'not-a-session.json', mimeType: 'application/json', buffer: Buffer.from('[1, 2, 3]\n')
  });
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending, null, { timeout: 60_000 });
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.errorLog?.code ?? null)).not.toBeNull();
  expect(await page.evaluate(() => window.__GBDRAW_APP__.errorLog.code)).not.toBe('UNKNOWN');
});

test('Reset Settings clears Linear per-record display text and the alignment plan', async ({ page }) => {
  test.fail(true, 'SE-10');
  test.setTimeout(300_000);
  await openFresh(page);
  await page.locator('input[accept^=".json,"]')
    .setInputFiles('gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json');
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending
    && window.__GBDRAW_APP__.sessionSaveAvailable, null, { timeout: 300_000 });
  await settle(page);
  const before = await page.evaluate(async () => {
    const { state } = await import('/gbdraw/web/js/state.js');
    return state.linearSeqs.map((sequence) => sequence.record_subtitle || '');
  });
  expect(before.some(Boolean)).toBe(true);
  await page.getByRole('button', { name: 'Reset Settings', exact: true }).click();
  await settle(page);
  const after = await page.evaluate(async () => {
    const { state } = await import('/gbdraw/web/js/state.js');
    return {
      definitions: state.linearSeqs.map((sequence) => sequence.definition || ''),
      subtitles: state.linearSeqs.map((sequence) => sequence.record_subtitle || ''),
      alignmentPlan: Boolean(state.similarityAlignmentPlan?.value)
    };
  });
  expect(after.definitions.every((value) => value === '')).toBe(true);
  expect(after.subtitles.every((value) => value === '')).toBe(true);
  expect(after.alignmentPlan).toBe(false);
});
