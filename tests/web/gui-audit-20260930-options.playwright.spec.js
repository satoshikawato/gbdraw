// Web GUI audit 2026-09-30: option values, tracks, History boundaries,
// comparisons and CLI Sessions. A test marked test.fail(true, '<ID>') asserts
// the correct behavior of a current defect; the fixing PR removes the mark.
const { test, expect } = require('@playwright/test');
const { execFile } = require('node:child_process');
const { readFileSync } = require('node:fs');
const path = require('node:path');
const { promisify } = require('node:util');
const { generateAndWaitForResult, reveal } = require('./helpers/app-lifecycle.cjs');
const {
  BATCH_FIXTURE,
  HMMT,
  HMMT_SESSION,
  loadSessionFile,
  openFresh,
  openWithGenBank,
  settle
} = require('./helpers/audit-browser.cjs');

test.describe.configure({ retries: 0 });

const root = process.cwd();
const BATCH_RECORDS = readFileSync(BATCH_FIXTURE, 'utf8').split(/^\/\/\s*$/m)
  .map((chunk) => chunk.trim()).filter((chunk) => chunk.startsWith('LOCUS')).map((chunk) => `${chunk}\n//\n`);

const dinucleotideSection = async (page) => {
  const summary = page.locator('summary[aria-label="Dinucleotide content/skew"]');
  await reveal(summary);
  await summary.evaluate((element) => { element.parentElement.open = true; });
  return summary.locator('xpath=..');
};

const withBoundedWait = (page, predicate, argument, timeout = 10_000) => page
  .waitForFunction(predicate, argument, { timeout }).then(() => true, () => false);

test('a GC window of 0 is rejected instead of silently becoming Auto', async ({ page }) => {
  test.fail(true, 'TR-04');
  test.setTimeout(300_000);
  await openWithGenBank(page, HMMT);
  const window = (await dinucleotideSection(page)).locator('input[type="number"]').first();
  await window.fill('0');
  await window.press('Tab');
  await generateAndWaitForResult(page, { expectedStatus: 'error' });
  expect(await window.inputValue()).toBe('0');
});

test('a decimal GC step is rejected instead of silently becoming Auto', async ({ page }) => {
  test.fail(true, 'GE-08');
  test.setTimeout(300_000);
  await openWithGenBank(page, HMMT);
  const step = (await dinucleotideSection(page)).locator('input[type="number"]').nth(1);
  await step.fill('10.5');
  await step.press('Tab');
  await generateAndWaitForResult(page, { expectedStatus: 'error' });
  expect(await step.inputValue()).toBe('10.5');
});

test('a non-numeric comparison e-value is rejected instead of using the default filter', async ({ page }) => {
  test.fail(true, 'GE-04');
  test.setTimeout(300_000);
  await openFresh(page);
  const examples = path.join(root, 'examples');
  await page.evaluate(async ({ first, second, table }) => {
    const app = window.__GBDRAW_APP__;
    app.mode = 'linear';
    await window.Vue.nextTick();
    if (app.linearSeqs.length < 2) app.addLinearSeq();
    app.setLinearSeqPrimaryFile(0, 'gb', new File([first], 'MjeNMV.gb', { type: 'text/plain' }));
    app.setLinearSeqPrimaryFile(1, 'gb', new File([second], 'MelaMJNV.gb', { type: 'text/plain' }));
    app.linearComparisonPlan.mode = 'adjacent';
    app.linearComparisonPlan.defaultSource = 'upload';
    app.linearComparisonPlan.edges.splice(0, app.linearComparisonPlan.edges.length, {
      id: 'audit-edge', queryUid: app.linearSeqs[0].uid, subjectUid: app.linearSeqs[1].uid,
      included: true, fileActive: true, losatFilenameActive: false, source: 'upload',
      file: new File([table], 'MjeNMV.MelaMJNV.tblastx.out', { type: 'text/plain' }), losatFilename: ''
    });
    await window.Vue.nextTick();
  }, {
    first: readFileSync(path.join(examples, 'MjeNMV.gb'), 'utf8'),
    second: readFileSync(path.join(examples, 'MelaMJNV.gb'), 'utf8'),
    table: readFileSync(path.join(examples, 'MjeNMV.MelaMJNV.tblastx.out'), 'utf8')
  });
  await settle(page);
  const evalue = page.getByLabel('Linear comparison maximum e-value');
  await reveal(evalue);
  await evalue.fill('1e-50x');
  await evalue.blur();
  await generateAndWaitForResult(page, { expectedStatus: 'error' });
});

test('a live stroke edit leaves the committed Result unchanged until Generate', async ({ page }) => {
  test.setTimeout(300_000);
  await openWithGenBank(page, HMMT);
  await generateAndWaitForResult(page);
  const committed = await page.evaluate(() => window.__GBDRAW_APP__.results[0].content);
  const input = await reveal(page.getByLabel('Block Stroke Width', { exact: true }).first());
  await input.fill('5');
  await input.fill('');
  const rewritten = await withBoundedWait(
    page,
    (original) => window.__GBDRAW_APP__.results[0].content !== original,
    committed
  );
  expect(rewritten).toBe(false);
});

const depthTsv = () => {
  const lines = ['reference_name\tposition\tdepth'];
  for (let position = 1; position <= 16569; position += 100) {
    lines.push(`NC_012920.1\t${position}\t${(10 + (position % 1000) / 50).toFixed(3)}`);
  }
  return `${lines.join('\n')}\n`;
};

const openCustomStackPanel = async (page) => {
  const button = page.locator('button[aria-controls="circular-custom-track-slots-panel"]');
  await reveal(button);
  if (await button.getAttribute('aria-expanded') !== 'true') await button.click();
  return page.getByText('Use custom stack', { exact: true }).locator('input');
};

const slots = (page) => page.evaluate(() => window.__GBDRAW_APP__.adv.circular_track_slots
  .map((slot) => ({ id: slot.id, renderer: slot.renderer, enabled: Boolean(slot.enabled) })));

test('removing the Depth file with the uploader Remove leaves a Generate-ready stack', async ({ page }) => {
  test.fail(true, 'TR-03');
  test.setTimeout(300_000);
  await openWithGenBank(page, HMMT);
  await page.evaluate((text) => window.__GBDRAW_APP__.setCircularDepthFile(
    0, new File([text], 'sampleA.depth.tsv', { type: 'text/tab-separated-values' })
  ), depthTsv());
  await (await openCustomStackPanel(page)).check();
  await settle(page);
  const summary = page.locator('summary[aria-label="Depth TSV tracks"]');
  await reveal(summary);
  await summary.evaluate((element) => { element.parentElement.open = true; });
  await summary.locator('xpath=..').locator('[role="group"][aria-label$="selection"] button', { hasText: 'Remove' })
    .first().click();
  await settle(page);
  expect((await slots(page)).filter(({ renderer, enabled }) => renderer === 'depth' && enabled)).toEqual([]);
  await generateAndWaitForResult(page);
});

test('un-hiding GC while the custom stack is inactive restores the GC content row', async ({ page }) => {
  test.fail(true, 'TR-06');
  test.setTimeout(300_000);
  await openWithGenBank(page, HMMT);
  const useStack = await openCustomStackPanel(page);
  const hideGc = await reveal(page.getByRole('checkbox', { name: 'Hide GC Content', exact: true }));
  await useStack.check();
  await hideGc.check();
  await useStack.uncheck();
  await hideGc.uncheck();
  await useStack.check();
  await settle(page);
  expect((await slots(page)).find(({ id }) => id === 'gc_content')?.enabled).toBe(true);
});

test('a preset reset honors Show Coordinate Scale off', async ({ page }) => {
  test.fail(true, 'TR-09');
  test.setTimeout(300_000);
  await openWithGenBank(page, HMMT, () => { window.__GBDRAW_APP__.form.show_scale = false; });
  await (await openCustomStackPanel(page)).check();
  await page.getByRole('button', { name: 'Reset to Middle', exact: true }).click();
  await settle(page);
  expect((await slots(page)).filter(({ renderer, enabled }) => renderer === 'ticks' && enabled)).toEqual([]);
});

test('clicking the text of a checkbox label records an Undo step', async ({ page }) => {
  test.setTimeout(300_000);
  await openFresh(page);
  await loadSessionFile(page, HMMT_SESSION);
  const label = page.locator('label.option-label', { hasText: 'Rich Feature Popup' }).first();
  await reveal(label);
  const before = await page.evaluate(() => ({
    value: window.__GBDRAW_APP__.adv.rich_feature_popup, undo: window.__GBDRAW_HISTORY__.getUndoCount()
  }));
  await label.locator('span').first().click();
  await settle(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.adv.rich_feature_popup)).toBe(!before.value);
  expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(before.undo + 1);
  await page.evaluate(() => window.__GBDRAW_HISTORY__.undo());
  await settle(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.adv.rich_feature_popup)).toBe(before.value);
});

test('a checkbox click while a text field has focus records its own Undo step', async ({ page }) => {
  test.setTimeout(300_000);
  await openFresh(page);
  await loadSessionFile(page, HMMT_SESSION);
  const prefix = await reveal(page.locator('#output-prefix'));
  const checkbox = page.locator('label.option-label', { hasText: 'Rich Feature Popup' }).first()
    .locator('input[type=checkbox]');
  await reveal(checkbox);
  const before = await page.evaluate(() => ({
    value: window.__GBDRAW_APP__.adv.rich_feature_popup, undo: window.__GBDRAW_HISTORY__.getUndoCount()
  }));
  await prefix.fill('audit');
  await checkbox.click();
  await settle(page);
  expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(before.undo + 2);
  await page.evaluate(() => window.__GBDRAW_HISTORY__.undo());
  await settle(page);
  expect(await page.evaluate(() => ({
    value: window.__GBDRAW_APP__.adv.rich_feature_popup, prefix: window.__GBDRAW_APP__.form.prefix
  }))).toEqual({ value: before.value, prefix: 'audit' });
});

test('Undo while a Generate is running is rejected as busy', async ({ page }) => {
  test.setTimeout(300_000);
  await openWithGenBank(page, HMMT, () => { window.__GBDRAW_APP__.adv.scale_interval = 2000; });
  await generateAndWaitForResult(page);
  await page.evaluate(() => { window.__GBDRAW_APP__.adv.scale_interval = 5000; });
  await generateAndWaitForResult(page);
  const committed = () => page.evaluate(async () => {
    const { getCommittedCanonicalRenderRequest } = await import('/gbdraw/web/js/services/config.js');
    return {
      request: JSON.stringify(getCommittedCanonicalRenderRequest()),
      undo: window.__GBDRAW_HISTORY__.getUndoCount(),
      redo: window.__GBDRAW_HISTORY__.getRedoCount()
    };
  });
  const before = await committed();
  await page.evaluate(() => {
    window.__GBDRAW_APP__.adv.scale_interval = 1000;
    window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse = () => new Promise((resolve) => {
      window.__AUDIT_RELEASE_GENERATE__ = resolve;
    });
  });
  const afterEdit = await committed();
  await page.getByRole('button', { name: 'Generate Diagram', exact: true }).click();
  await page.waitForFunction(() => Boolean(window.__AUDIT_RELEASE_GENERATE__), null, { timeout: 120_000 });
  await page.locator('body').press('Control+z');
  await page.evaluate(() => new Promise((resolve) => requestAnimationFrame(() => requestAnimationFrame(resolve))));
  const during = await committed();
  await page.evaluate(() => {
    window.__AUDIT_RELEASE_GENERATE__();
    delete window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse;
  });
  await page.waitForFunction(() => !window.__GBDRAW_APP__.processing, null, { timeout: 120_000 });
  expect(during.request).toBe(before.request);
  expect({ undo: during.undo, redo: during.redo }).toEqual({ undo: afterEdit.undo, redo: afterEdit.redo });
  expect(await page.evaluate(() => window.__GBDRAW_APP__.adv.scale_interval)).toBe(1000);
});

test('default threaded LOSAT without cross-origin isolation reports a recognized diagnostic', async ({ page }) => {
  test.fail(true, 'CO-01');
  test.setTimeout(300_000);
  await openFresh(page);
  await page.getByRole('button', { name: 'Linear', exact: true }).click();
  await page.waitForFunction(() => window.__GBDRAW_APP__?.mode === 'linear');
  expect(await page.evaluate(() => globalThis.crossOriginIsolated)).toBe(false);
  await page.evaluate(async (records) => {
    const app = window.__GBDRAW_APP__;
    while (app.linearSeqs.length < records.length) app.addLinearSeq();
    records.forEach((text, index) => app.setLinearSeqPrimaryFile(
      index, 'gb', new File([text], `record-${index + 1}.gb`, { type: 'text/plain', lastModified: 1000 + index })
    ));
    await window.Vue.nextTick();
    await app.setLinearComparisonGlobalAction('losat');
    app.setLinearComparisonLosatMode('blastp');
  }, BATCH_RECORDS);
  await settle(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.losat.executionMode)).toBe('threaded');
  const outcome = await generateAndWaitForResult(page, { expectedStatus: 'error', requireCommittedResult: false });
  expect(outcome.health.errorCode).not.toBe('UNKNOWN');
  expect(await page.evaluate(() => window.__GBDRAW_APP__.errorLog?.stage)).not.toBe('request-validation');
});

test('a CLI Session keeps its legend position through load and the first Generate', async ({ page }, testInfo) => {
  test.fail(true, 'SE-07');
  test.setTimeout(600_000);
  const prefix = testInfo.outputPath('cli-legend');
  const session = `${prefix}.gbdraw-session.json.gz`;
  await promisify(execFile)('python', [
    '-m', 'gbdraw.cli', 'circular', '--gbk', path.join(root, 'tests/fixtures/sessions/cli-web-mito.gb'),
    '--legend', 'upper_left', '-o', prefix, '--session_output', session
  ], { cwd: testInfo.outputDir, env: { ...process.env, PYTHONPATH: root }, timeout: 600_000, maxBuffer: 1_000_000 });
  await openFresh(page);
  await page.locator('input[accept^=".json,"]').setInputFiles(session);
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending
    && window.__GBDRAW_APP__.results.length > 0, null, { timeout: 300_000 });
  await settle(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.form.legend)).toBe('upper_left');
  await generateAndWaitForResult(page);
  expect(await page.evaluate(async () => (
    (await import('/gbdraw/web/js/state.js')).state.generatedLegendPosition.value
  ))).toBe('upper_left');
});
