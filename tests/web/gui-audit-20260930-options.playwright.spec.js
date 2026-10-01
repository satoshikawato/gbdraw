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
  settle,
  switchMode
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
  test.setTimeout(300_000);
  await openWithGenBank(page, HMMT);
  const window = (await dinucleotideSection(page)).locator('input[type="number"]').first();
  await window.fill('0');
  await window.press('Tab');
  await generateAndWaitForResult(page, { expectedStatus: 'error' });
  expect(await window.inputValue()).toBe('0');
});

test('a decimal GC step is rejected instead of silently becoming Auto', async ({ page }) => {
  test.setTimeout(300_000);
  await openWithGenBank(page, HMMT);
  const step = (await dinucleotideSection(page)).locator('input[type="number"]').nth(1);
  await step.fill('10.5');
  await step.press('Tab');
  await generateAndWaitForResult(page, { expectedStatus: 'error' });
  expect(await step.inputValue()).toBe('10.5');
});

test('a non-numeric comparison e-value is rejected instead of using the default filter', async ({ page }) => {
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

const depthRows = (page, key = 'circular_track_slots') => page.evaluate((slotsKey) => (
  window.__GBDRAW_APP__.adv[slotsKey].filter((slot) => slot.renderer === 'depth').map((slot) => ({
    id: slot.id, enabled: slot.enabled !== false, trackIndex: slot.params?.track_index ?? null
  }))
), key);

const uploaderRemove = (container) => container.locator('[role="group"][aria-label$="selection"] button', { hasText: 'Remove' }).first();

// TR-03 / PD-OI-058: the uploader Remove is a Depth source change, so the
// managed row goes with the file in the same History step.
test('removing the Depth file with the uploader Remove leaves a Generate-ready stack', async ({ page }) => {
  test.setTimeout(300_000);
  await openWithGenBank(page, HMMT);
  await page.evaluate((text) => window.__GBDRAW_APP__.setCircularDepthFile(
    0, new File([text], 'sampleA.depth.tsv', { type: 'text/tab-separated-values' })
  ), depthTsv());
  await (await openCustomStackPanel(page)).check();
  await settle(page);
  expect(await depthRows(page)).toEqual([{ id: 'depth', enabled: true, trackIndex: 0 }]);
  const summary = page.locator('summary[aria-label="Depth TSV tracks"]');
  await reveal(summary);
  await summary.evaluate((element) => { element.parentElement.open = true; });
  await uploaderRemove(summary.locator('xpath=..')).click();
  await settle(page);
  expect((await slots(page)).filter(({ renderer, enabled }) => renderer === 'depth' && enabled)).toEqual([]);
  await page.evaluate(() => window.__GBDRAW_HISTORY__.undo());
  await settle(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.getCircularDepthFile(0)?.name)).toBe('sampleA.depth.tsv');
  expect(await depthRows(page)).toEqual([{ id: 'depth', enabled: true, trackIndex: 0 }]);
  await page.evaluate(() => window.__GBDRAW_HISTORY__.redo());
  await settle(page);
  expect(await depthRows(page)).toEqual([]);
  await generateAndWaitForResult(page);
});

const openLinearStack = async (page) => {
  await openFresh(page);
  await switchMode(page, 'linear');
  await page.getByTestId('linear-genbank-1').setInputFiles(HMMT);
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.linearSeqs[0]?.gb?.name)).toBe('HmmtDNA.gbk');
  const button = page.locator('button[aria-controls="linear-custom-track-slots-panel"]');
  await reveal(button);
  if (await button.getAttribute('aria-expanded') !== 'true') await button.click();
  await page.getByText('Use custom stack', { exact: true }).locator('input').check();
  await settle(page);
  const source = page.locator('[data-linear-source-depth]').first();
  await reveal(source.locator(':scope > summary'));
  await source.evaluate((element) => { element.open = true; });
  return source;
};

const setLinearDepthFromFile = async (page) => {
  await page.getByTestId('linear-source-depth-1-1').setInputFiles({
    name: 'sampleA.depth.tsv', mimeType: 'text/tab-separated-values', buffer: Buffer.from(depthTsv())
  });
  await settle(page);
};

// TR-02 / PD-OI-058: both modes add a managed row when a series gets its first
// source, whether or not the custom stack is in use.
test('a Linear Depth file adds one managed row to the custom stack', async ({ page }) => {
  test.setTimeout(300_000);
  await openLinearStack(page);
  expect(await depthRows(page, 'linear_track_slots')).toEqual([]);
  await setLinearDepthFromFile(page);
  expect(await depthRows(page, 'linear_track_slots')).toEqual([{ id: 'depth', enabled: true, trackIndex: 0 }]);
});

test('deleting the Linear File that held a series source removes its managed row', async ({ page }) => {
  test.setTimeout(300_000);
  await openLinearStack(page);
  await page.evaluate(() => window.__GBDRAW_APP__.addLinearSeq());
  await page.getByTestId('linear-genbank-2').setInputFiles(HMMT);
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.linearSeqs[1]?.gb?.name)).toBe('HmmtDNA.gbk');
  await page.evaluate((text) => {
    const app = window.__GBDRAW_APP__;
    app.setLinearDepthFile(app.linearSeqs[1], 0, new File([text], 'sampleB.depth.tsv', { type: 'text/tab-separated-values' }));
  }, depthTsv());
  await settle(page);
  expect(await depthRows(page, 'linear_track_slots')).toEqual([{ id: 'depth', enabled: true, trackIndex: 0 }]);
  await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    app.requestLinearSourceRemoval(app.linearSourceGroups[1]);
    await app.applyLinearSourceRemoval('delete');
  });
  await settle(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.linearSeqs.length)).toBe(1);
  expect(await depthRows(page, 'linear_track_slots')).toEqual([]);
});

const openLinearDepthStack = async (page) => {
  const source = await openLinearStack(page);
  await setLinearDepthFromFile(page);
  await page.evaluate(() => window.__GBDRAW_APP__.resetLinearTrackSlotsFromSimpleControls());
  await settle(page);
  return source;
};

// TR-03 Linear / PD-OI-083: a File clear keeps the logical series (PD-OI-025);
// an enabled manual row on it shows the row issue and stops Generate.
test('a Linear File clear leaves a manual Depth row whose row issue stops Generate', async ({ page }) => {
  test.setTimeout(300_000);
  const source = await openLinearDepthStack(page);
  expect(await depthRows(page, 'linear_track_slots')).toEqual([{ id: 'depth', enabled: true, trackIndex: 0 }]);
  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    app.setLinearTrackSlotHeight(app.adv.linear_track_slots.find((slot) => slot.renderer === 'depth'), '30');
  });
  await uploaderRemove(source).click();
  await settle(page);
  const state = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const index = app.adv.linear_track_slots.findIndex((slot) => slot.renderer === 'depth');
    return {
      index,
      width: app.linearSeqs[0].depth.length,
      issue: app.linearTrackSlotIssue(app.adv.linear_track_slots[index], index)
    };
  });
  expect(state).toMatchObject({ width: 1, issue: "Linear Depth track 'depth' has no logical Depth source." });
  expect(await depthRows(page, 'linear_track_slots')).toEqual([{ id: 'depth', enabled: true, trackIndex: 0 }]);
  await generateAndWaitForResult(page, { expectedStatus: 'error' });
  expect(await page.evaluate(() => {
    const { code, context } = window.__GBDRAW_APP__.errorLog;
    return { code, reason: context.reason, slotIndex: context.slotIndex };
  })).toEqual({ code: 'TRACK_INVALID', reason: 'REQUIRED', slotIndex: state.index });
  await page.evaluate(() => {
    window.__GBDRAW_APP__.adv.linear_track_slots.find((slot) => slot.renderer === 'depth').enabled = false;
  });
  await settle(page);
  await generateAndWaitForResult(page);
});

// TR-02 Linear / PD-OI-058 May retire: Add Depth TSV series adds a series
// without a source, so it neither adds nor re-enables a row.
test('Add Depth TSV series keeps a disabled Linear Depth row disabled', async ({ page }) => {
  test.setTimeout(300_000);
  await openLinearDepthStack(page);
  await page.evaluate(() => {
    const slot = window.__GBDRAW_APP__.adv.linear_track_slots.find((entry) => entry.renderer === 'depth');
    slot.enabled = false;
    slot.params = { ...slot.params, legend_label: 'Coverage A' };
  });
  const before = await page.evaluate(() => JSON.parse(JSON.stringify(window.__GBDRAW_APP__.adv.linear_track_slots)));
  await page.getByRole('button', { name: 'Add Depth TSV series from file 1', exact: true }).click();
  await settle(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.linearSeqs[0].depth.length)).toBe(2);
  expect(await page.evaluate(() => JSON.parse(JSON.stringify(window.__GBDRAW_APP__.adv.linear_track_slots)))).toEqual(before);
});

test('un-hiding GC while the custom stack is inactive restores the GC content row', async ({ page }) => {
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
