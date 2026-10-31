// The input summary (gbdraw/api/input_summary.py): a record with no feature to draw still draws, and the status
// area shows the warning the command line prints; Run Info lists the summary lines. Python owns every word.
const { readFileSync } = require('node:fs');
const { join } = require('node:path');
const { test, expect } = require('@playwright/test');
const { evaluateWithRetainedPromise, generateAndWaitForResult, openApp } = require('./helpers/app-lifecycle.cjs');

const ZERO_FEATURES = () => readFileSync(join(process.cwd(), 'tests/fixtures/input_variants/zero_features_record.gb'));

test('a record with no feature to draw: the status area says so and Run Info lists the summary', async ({ page }) => {
  test.setTimeout(300000);
  await openApp(page);
  await page.getByLabel('GenBank/DDBJ File', { exact: true })
    .setInputFiles({ name: 'zero_features_record.gb', mimeType: 'text/plain', buffer: ZERO_FEATURES() });
  await page.waitForFunction(() => window.__GBDRAW_APP__.circularRecordList.length === 1);
  await generateAndWaitForResult(page);

  const notice = page.getByTestId('input-summary-notice');
  await expect(notice).toBeVisible();
  await expect(notice.locator('p').first()).toHaveText('Part of the input was not drawn.');
  await expect(notice.locator('[data-input-warning="no_drawable_features"]')).toHaveText(
    'Record EMPTY001.1 (zero_features_record.gb) has no feature to draw: it has no features.'
  );

  await page.evaluate(() => { window.__GBDRAW_APP__.resultPanelTab = 'run-info'; });
  await expect(page.locator('[data-run-info-input-summary] li')).toHaveText([
    'Record EMPTY001.1 (zero_features_record.gb), record 1 of 1, 1,000 bp: drew 0 features.'
  ]);
});

// A record of `bp` bases with a CDS per locus tag, in GenBank text.
const genbankRecord = (name, tags) => {
  const bp = 600;
  const sequence = 'acgtacgtac'.repeat(bp / 10);
  const lines = [];
  for (let start = 0; start < bp; start += 60) {
    const blocks = sequence.slice(start, start + 60).match(/.{1,10}/g).join(' ');
    lines.push(`${String(start + 1).padStart(9)} ${blocks}`);
  }
  const cds = tags.flatMap((tag, index) => [
    `     CDS             ${index * 200 + 10}..${index * 200 + 190}`,
    `                     /locus_tag="${tag}"`,
    `                     /product="protein ${tag}"`
  ]);
  return [
    `LOCUS       ${name.padEnd(16)} ${String(bp).padStart(11)} bp    DNA     linear   UNK 01-JAN-1980`,
    `DEFINITION  ${name}.`, `ACCESSION   ${name}`, `VERSION     ${name}.1`,
    'FEATURES             Location/Qualifiers', `     source          1..${bp}`, ...cds, 'ORIGIN', ...lines, '//', ''
  ].join('\n');
};
const genbankFile = (name, records) => ({ name, mimeType: 'text/plain', buffer: Buffer.from(records.join('')) });

const settled = (page) => page.waitForFunction(() => {
  const app = window.__GBDRAW_APP__;
  const history = window.__GBDRAW_HISTORY__;
  return !app.processing && !app.labelReflowProcessing && !app.ruleMatchingPending
    && !history?.capturing?.value && !history?.restoring?.value;
}, null, { timeout: 300_000 });

// What the page says of the shown Result's input summary.
const summaryShown = async (page) => {
  await settled(page);
  await page.evaluate(() => { window.__GBDRAW_APP__.resultPanelTab = 'run-info'; });
  const runInfo = page.locator('[data-run-info-input-summary]');
  return {
    notice: await page.getByTestId('input-summary-notice').isVisible()
      ? await page.locator('[data-testid="input-summary-notice"] [data-input-warning]').allTextContents() : null,
    staleLine: await page.getByTestId('input-summary-stale').isVisible(),
    runInfo: await runInfo.locator('li').allTextContents(),
    runInfoStale: await runInfo.locator('[data-run-info-input-summary-stale]').isVisible()
  };
};

const addVisibilityRule = (page, fields) => evaluateWithRetainedPromise(page, async (ruleFields) => {
  const app = window.__GBDRAW_APP__;
  await app.addFeatureVisibilityRule();
  const index = app.featureVisibilityManualRules.length - 1;
  for (const [field, value] of Object.entries(ruleFields)) await app.setFeatureVisibilityRuleField(index, field, value);
}, fields);
const setRuleValue = (page, value) => evaluateWithRetainedPromise(
  page, (text) => window.__GBDRAW_APP__.setFeatureVisibilityRuleField(0, 'value', text), value
);
const history = (page, step) => evaluateWithRetainedPromise(page, (name) => window.__GBDRAW_HISTORY__[name](), step);
const openRecords = async (page, file) => {
  page.on('dialog', (dialog) => dialog.accept());
  await openApp(page);
  await page.getByLabel('GenBank/DDBJ File', { exact: true }).setInputFiles(file);
  await page.waitForFunction(() => window.__GBDRAW_APP__.circularRecordList.length > 0);
};

const TWO_CDS = genbankFile('two_cds.gb', [genbankRecord('VIS1', ['VA', 'VB'])]);
const HIDE_VB = 'Record VIS1.1 (two_cds.gb), record 1 of 1, 600 bp: drew 1 feature (1 CDS); left out 1 source (type not drawn); hid 1 CDS (visibility rule or feature edit).';

// Python counts the Feature visibility of the render that drew the Result. A
// live edit the rerender draws (a feature the Result does not draw, a Legend
// row that goes) brings the rerender's summary, so live equals Generate.
test('a live Show and a live hide of every feature that the rerender draws: live equals Generate', async ({ page }) => {
  test.setTimeout(600_000);
  await openRecords(page, TWO_CDS);
  await addVisibilityRule(page, { recordId: '*', featureType: 'CDS', qualifier: 'locus_tag', value: '^V', action: 'off' });
  await generateAndWaitForResult(page);
  const allHidden = await summaryShown(page);
  expect(allHidden.notice).toEqual([expect.stringMatching(/^Record VIS1\.1 \(two_cds\.gb\) has no feature to draw/)]);
  expect(allHidden.runInfoStale).toBe(false);

  // Show one live: the rerender draws VA, so its summary replaces the Generate's.
  await setRuleValue(page, '^VB$');
  const shownOne = await summaryShown(page);
  expect(shownOne).toEqual({ notice: null, staleLine: false, runInfo: [HIDE_VB], runInfoStale: false });
  await generateAndWaitForResult(page);
  expect(await summaryShown(page)).toEqual(shownOne);

  // Hide every feature live: the CDS Legend row goes, so the rerender counts the hide.
  await setRuleValue(page, '^V');
  expect(await summaryShown(page)).toEqual(allHidden);
  await history(page, 'undo');
  expect(await summaryShown(page)).toEqual(shownOne);
  await history(page, 'redo');
  expect(await summaryShown(page)).toEqual(allHidden);
  await generateAndWaitForResult(page);
  expect(await summaryShown(page)).toEqual(allHidden);
});

// A live hide that no render draws (the other CDS keeps the Legend row) leaves
// the counts of the last render: the page marks them so and holds back their
// warning, Undo back to the rendered visibility takes the mark away, Redo
// brings it back, and Generate counts the hide.
test('a live hide no render draws marks the summary as the last render\'s until Undo or Generate', async ({ page }) => {
  test.setTimeout(600_000);
  await openRecords(page, genbankFile('batch.gb', [genbankRecord('EMPTY1', []), genbankRecord('VIS1', ['VA', 'VB'])]));
  await page.evaluate(() => { window.__GBDRAW_APP__.form.multi_record_canvas = true; });
  await generateAndWaitForResult(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.results.length)).toBe(1);
  const rendered = await summaryShown(page);
  expect(rendered.notice).toEqual([expect.stringMatching(/^Record EMPTY1\.1 \(batch\.gb\) has no feature to draw/)]);
  expect(rendered.runInfo[1]).toMatch(/^Record VIS1\.1 \(batch\.gb\), record 2 of 2, 600 bp: drew 2 features/);

  const undoCount = () => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
  const stepsBefore = await undoCount();
  await addVisibilityRule(page, { recordId: '*', featureType: 'CDS', qualifier: 'locus_tag', value: '^VA$', action: 'off' });
  const hidden = await summaryShown(page);
  // Adding the rule and setting its fields are steps of their own.
  const steps = await undoCount() - stepsBefore;
  expect(steps).toBeGreaterThan(0);
  expect(hidden).toEqual({ notice: null, staleLine: true, runInfo: rendered.runInfo, runInfoStale: true });
  await expect(page.getByTestId('input-summary-stale')).toHaveText(
    'The input summary counts the last render. Feature visibility changed since; Generate to count it again.'
  );
  for (let step = 0; step < steps; step += 1) await history(page, 'undo');
  expect(await summaryShown(page)).toEqual(rendered);
  for (let step = 0; step < steps; step += 1) await history(page, 'redo');
  expect(await summaryShown(page)).toEqual(hidden);
  await generateAndWaitForResult(page);
  const counted = await summaryShown(page);
  expect(counted).toEqual({ ...rendered, runInfo: [rendered.runInfo[0], expect.stringMatching(/drew 1 feature \(1 CDS\).*hid 1 CDS/)] });
});

// Each Result of a batch shows its own summary, through a Result switch.
test('a batch Result switch shows each Result\'s own summary', async ({ page }) => {
  test.setTimeout(600_000);
  await openRecords(page, genbankFile('batch.gb', [genbankRecord('EMPTY1', []), genbankRecord('VIS1', ['VA', 'VB'])]));
  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    app.form.multi_record_canvas = false;
    app.adv.circular_grouping_intent = 'batch';
  });
  await generateAndWaitForResult(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.results.length)).toBe(2);
  const select = async (index) => {
    await page.locator('h2 select').selectOption({ index });
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.selectedResultIndex)).toBe(index);
    return summaryShown(page);
  };
  const first = await summaryShown(page);
  expect(first.notice).toEqual([expect.stringMatching(/^Record EMPTY1\.1 \(batch\.gb\) has no feature to draw/)]);
  const second = await select(1);
  expect(second.notice).toBeNull();
  expect(second.runInfo).toEqual([expect.stringMatching(/^Record VIS1\.1 \(batch\.gb\), record 2 of 2, 600 bp: drew 2 features/)]);
  expect(await select(0)).toEqual(first);
});
