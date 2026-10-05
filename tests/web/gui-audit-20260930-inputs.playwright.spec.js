// Web GUI audit 2026-09-30: inputs, record topology, Session Save/Load and
// Reset. A test marked test.fail(true, '<ID>') asserts the correct behavior of
// a current defect; the PR that fixes the audit ID removes the mark.
const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { generateAndWaitForResult, evaluateWithRetainedPromise } = require('./helpers/app-lifecycle.cjs');
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
const SECOND_RECORD_TEXT = `${BATCH_TEXT.split(/^\/\/\s*$/m)[1].trim()}\n//\n`;

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

const RECORD_VECTORS = JSON.parse(readFileSync('tests/fixtures/record_metadata_inference_cases.json', 'utf8'))
  .discovery;
const vectorText = ({ lines }) => `${lines.join('\n')}\n`;
const splitRecords = (file) => readFileSync(file, 'utf8').split(/^\/\/\s*$/m)
  .map((chunk) => chunk.trim()).filter((chunk) => chunk.startsWith('LOCUS')).map((chunk) => `${chunk}\n//\n`);
// Prokka writes empty ACCESSION and VERSION lines; Biopython then names the record by LOCUS.
const prokkaStyle = (record) => record.replace(/^ACCESSION .*$/m, 'ACCESSION   ').replace(/^VERSION .*$/m, 'VERSION');
// R2c 2001..3000 and R3c 1..1000 are one shared 1 kb block.
const [PROKKA_R2C, PROKKA_R3C] = splitRecords('tests/fixtures/web_comparison_shared_block.gb').map(prokkaStyle);

const resultRecordIds = (page) => page.evaluate(() => {
  const content = window.__GBDRAW_APP__.results[0]?.content || '';
  const document = new DOMParser().parseFromString(content, 'image/svg+xml');
  return [...new Set([...document.querySelectorAll('[data-gbdraw-record-id]')]
    .map((element) => element.getAttribute('data-gbdraw-record-id')))].sort();
});

test('Prokka-style GenBank records keep their LOCUS names through Circular Generate and Linear LOSAT', async ({ page }) => {
  test.setTimeout(600_000);
  await openFresh(page);
  await page.getByLabel('GenBank/DDBJ File', { exact: true })
    .setInputFiles({ name: 'R2c.gbk', mimeType: 'text/plain', buffer: Buffer.from(PROKKA_R2C) });
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.circularRecordList
    .map((record) => record.record_id)), { timeout: 60_000 }).toEqual(['R2c']);
  await settle(page);
  await generateAndWaitForResult(page);
  expect(await resultRecordIds(page)).toEqual(['R2c']);

  await page.getByRole('button', { name: 'Linear', exact: true }).click();
  await page.waitForFunction(() => window.__GBDRAW_APP__?.mode === 'linear');
  await page.evaluate(async (texts) => {
    const app = window.__GBDRAW_APP__;
    Object.assign(app.form, { legend: 'none', show_gc: false, show_skew: false });
    while (app.linearSeqs.length < texts.length) app.addLinearSeq();
    texts.forEach((text, index) => app.setLinearSeqPrimaryFile(index, 'gb', new File(
      [text], `prokka-${index + 1}.gbk`, { type: 'text/plain', lastModified: 1000 + index }
    )));
    await window.Vue.nextTick();
    await app.setLinearComparisonGlobalAction('losat');
    app.setLinearComparisonLosatMode('blastn');
    // Serial LOSAT: the Playwright server does not send COOP/COEP.
    app.losat.executionMode = 'serial';
  }, [PROKKA_R2C, PROKKA_R3C]);
  await settle(page);
  await generateAndWaitForResult(page);
  expect(await resultRecordIds(page)).toEqual(['R2c', 'R3c']);
  expect(await page.evaluate(() => {
    const content = window.__GBDRAW_APP__.results[0]?.content || '';
    return new DOMParser().parseFromString(content, 'image/svg+xml')
      .querySelectorAll('[data-gbdraw-pairwise-match-id]').length;
  })).toBeGreaterThan(0);
});

test('GFF3+FASTA records without GFF3 lines are not offered and both modes generate the annotated records', async ({ page }) => {
  test.setTimeout(300_000);
  const pair = RECORD_VECTORS.find(({ name }) => name === 'GFF3+FASTA whose FASTA holds a sequence without GFF3 lines');
  const annotated = pair.expected.records.map(({ recordId }) => recordId);
  const pairFiles = (gff, fasta) => [
    { name: gff.name, mimeType: 'text/plain', buffer: Buffer.from(vectorText(gff)) },
    { name: fasta.name, mimeType: 'text/plain', buffer: Buffer.from(vectorText(fasta)) }
  ];
  const [gffFile, fastaFile] = pairFiles(pair.files.gff, pair.files.fasta);
  await openFresh(page);
  await page.locator('input[type=radio][value=gff]').first().check();
  await page.getByLabel('GFF3 File', { exact: true }).setInputFiles(gffFile);
  await page.getByLabel('FASTA File', { exact: true }).setInputFiles(fastaFile);
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.circularRecordList
    .map((record) => record.record_id)), { timeout: 60_000 }).toEqual(annotated);
  await settle(page);
  await generateAndWaitForResult(page);
  expect(await resultRecordIds(page)).toEqual(annotated);

  await page.getByRole('button', { name: 'Linear', exact: true }).click();
  await page.waitForFunction(() => window.__GBDRAW_APP__?.mode === 'linear');
  await page.evaluate(async (files) => {
    const app = window.__GBDRAW_APP__;
    Object.assign(app.form, { legend: 'none', show_gc: false, show_skew: false });
    app.setLinearInputType('gff');
    app.setLinearSeqPrimaryFile(0, 'gff', new File([files.gff], 'pair.gff3', { type: 'text/plain' }));
    app.setLinearSeqPrimaryFile(0, 'fasta', new File([files.fasta], 'pair.fasta', { type: 'text/plain' }));
    await window.Vue.nextTick();
  }, { gff: vectorText(pair.files.gff), fasta: vectorText(pair.files.fasta) });
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.linearSeqs
    .map((sequence) => sequence.region_record_id)), { timeout: 60_000 }).toEqual(annotated);
  await settle(page);
  await generateAndWaitForResult(page);
  expect(await resultRecordIds(page)).toEqual(annotated);
});

test('a Circular definition edit after Generate leaves the committed Result unchanged', async ({ page }) => {
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
  const rows = () => page.evaluate(() => window.__GBDRAW_APP__.linearSeqs.map((sequence) => ({
    start: sequence.region_start ?? null,
    end: sequence.region_end ?? null,
    definition: sequence.definition || ''
  })));
  // D-32: replacing one record with one record keeps its crop and Definition.
  await uploadLinear(page, 0, 'second.gb', SECOND_RECORD_TEXT);
  expect(await rows()).toEqual([{ start: 100, end: 3000, definition: 'Genome A custom' }]);
  await uploadLinear(page, 0, 'two.gb', BATCH_TEXT);
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.linearSeqs.length), { timeout: 30_000 }).toBe(2);
  expect(await rows())
    .toEqual([{ start: null, end: null, definition: '' }, { start: null, end: null, definition: '' }]);
  await generateAndWaitForResult(page);
});

// The shared record-discovery vector (G-D) supplies the loader's per-record
// inferred definitions for one file whose records name different organisms.
const TWO_ORGANISMS = JSON.parse(readFileSync('tests/fixtures/record_metadata_inference_cases.json', 'utf8'))
  .discovery.find(({ name }) => name === 'records of one file with different organisms');

test('each record of a multi-record Linear file keeps its own organism', async ({ page }) => {
  test.setTimeout(300_000);
  await openLinear(page);
  const source = TWO_ORGANISMS.files.source;
  await uploadLinear(page, 0, source.name, `${source.lines.join('\n')}\n`);
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.linearSeqs.length), { timeout: 30_000 }).toBe(2);
  const expectOwnOrganisms = async () => {
    const definitions = await recordDefinitions(page);
    for (const record of TWO_ORGANISMS.expected.records) {
      expect(definitions[record.recordId]).toContain(record.inferredDefinition.replace(/<\/?i>/g, ''));
    }
  };

  // Reset returns a record to its own inferred definition.
  await page.evaluate(() => { window.__GBDRAW_APP__.linearSeqs[1].definition = 'Record override'; });
  await settle(page);
  await page.getByRole('button', { name: 'Reset Settings', exact: true }).click();
  await settle(page);
  await page.evaluate(() => window.__GBDRAW_APP__.setLinearRecordLayoutEnabled(false));
  await settle(page);
  await generateAndWaitForResult(page);
  await expectOwnOrganisms();

  // A File default the user enters applies to every record of the File.
  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    app.setLinearSourceDefaultDefinition(app.linearSourceGroups[0], 'Shared organism');
  });
  await settle(page);
  await generateAndWaitForResult(page);
  const shared = await recordDefinitions(page);
  expect(Object.keys(shared).sort()).toEqual(TWO_ORGANISMS.expected.records.map(({ recordId }) => recordId));
  for (const text of Object.values(shared)) expect(text.split(' / ')[0]).toBe('Shared organism');

  // A saved Session restores the same inferred definitions.
  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    app.setLinearSourceDefaultDefinition(app.linearSourceGroups[0], '');
  });
  await settle(page);
  await generateAndWaitForResult(page);
  await expectOwnOrganisms();
  const linearState = () => page.evaluate(() => window.__GBDRAW_APP__.linearSeqs.map((sequence) => ({
    definition: sequence.definition, file: sequence.file_definition, inferred: sequence.inferred_definition
  })));
  const before = await linearState();
  expect(before.map(({ inferred }) => inferred))
    .toEqual(TWO_ORGANISMS.expected.records.map(({ inferredDefinition }) => inferredDefinition));
  const downloaded = page.waitForEvent('download');
  await page.evaluate(() => { window.__GBDRAW_APP__.sessionTitle = 'organisms'; });
  await page.getByRole('banner').getByRole('button', { name: 'Save Session', exact: true }).click();
  const sessionPath = await (await downloaded).path();
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionSavePending, null, { timeout: 120_000 });
  await page.locator('input[accept^=".json,"]').setInputFiles(sessionPath);
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending
    && window.__GBDRAW_APP__.sessionSaveAvailable, null, { timeout: 180_000 });
  await settle(page);
  expect(await linearState()).toEqual(before);
  await generateAndWaitForResult(page);
  await expectOwnOrganisms();
});

test('Save before the first Generate with an empty Linear card explains the missing input', async ({ page }) => {
  test.setTimeout(300_000);
  await openLinear(page);
  await uploadLinear(page, 0, 'first.gb', FIRST_RECORD_TEXT);
  await page.evaluate(() => window.__GBDRAW_APP__.addLinearSeq());
  await settle(page);
  const outcome = await attemptSave(page);
  expect(outcome.downloads).toEqual([]);
  expect(outcome.error?.code).toBeTruthy();
  expect(outcome.error.code).not.toBe('UNKNOWN');
  // Save and Generate share one input check and explain the same missing input.
  expect(outcome.error.code).toBe('INPUT_REQUIRED');
  expect(outcome.error.summary).toMatch(/Sequence 2\./);
  await generateAndWaitForResult(page, { expectedStatus: 'error' });
  expect(await page.evaluate(() => {
    const { code, summary } = window.__GBDRAW_APP__.errorLog;
    return { code, summary };
  })).toEqual({ code: outcome.error.code, summary: outcome.error.summary });
});

test('a failed Save does not offer Save Session again', async ({ page }) => {
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

for (const [version, fixture] of [
  ['v39', 'tests/fixtures/sessions/BGC0000708-BGC0000713.v39.gbdraw-session.json.gz'],
  ['v33 (schema-v2)', 'tests/fixtures/sessions/BGC0000708-BGC0000713.schema-v2.gbdraw-session.json.gz']
]) test(`Save after loading a ${version} Session asks for Generate instead of failing with UNKNOWN`, async ({ page }) => {
  test.setTimeout(300_000);
  await openFresh(page);
  await page.locator('input[accept^=".json,"]')
    .setInputFiles(fixture);
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending, null, { timeout: 180_000 });
  await page.waitForFunction(() => window.__GBDRAW_APP__.sessionSaveAvailable, null, { timeout: 180_000 });
  await settle(page);
  const outcome = await attemptSave(page);
  expect(outcome.downloads).toEqual([]);
  expect(outcome.error?.code).toBeTruthy();
  expect(outcome.error.code).not.toBe('UNKNOWN');
  expect(outcome.error.summary).toMatch(/Generate/);
  // D-25 Must preserve: the panel runs the needed Generate, then Save succeeds.
  expect(outcome.error.code).toBe('SESSION_SAVE_REQUIRES_GENERATE');
  await expect(page.locator('[data-session-save-needs-generate]')).toBeVisible();
  const alert = page.locator('.border-l-red-500');
  await expect(alert.getByRole('button', { name: 'Save Session', exact: true })).toHaveCount(0);
  await alert.getByRole('button', { name: 'Generate', exact: true }).click();
  await page.waitForFunction(() => !window.__GBDRAW_APP__.processing, null, { timeout: 300_000 });
  await settle(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.errorLog?.code ?? null)).toBeNull();
  await expect(page.locator('[data-session-save-needs-generate]')).toHaveCount(0);
  const saved = await attemptSave(page);
  expect(saved.error).toBeNull();
  expect(saved.downloads).toHaveLength(1);
});

test('loading a JSON file that is not a Session reports a recognized diagnostic', async ({ page }) => {
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
  const planBefore = await page.evaluate(async () => {
    const { state } = await import('/gbdraw/web/js/state.js');
    return Boolean(state.similarityAlignmentPlan?.value);
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
  // D-15: Undo restores the record display text and the alignment plan.
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
  await settle(page);
  const restored = await page.evaluate(async () => {
    const { state } = await import('/gbdraw/web/js/state.js');
    return {
      subtitles: state.linearSeqs.map((sequence) => sequence.record_subtitle || ''),
      alignmentPlan: Boolean(state.similarityAlignmentPlan?.value)
    };
  });
  expect(restored.subtitles).toEqual(before);
  expect(restored.alignmentPlan).toBe(planBefore);
});
