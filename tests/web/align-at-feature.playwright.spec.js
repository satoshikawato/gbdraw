// Align at feature (Linear) asks Python for the record translations that
// `gbdraw linear --anchor_feature` draws with, previews which records move, and
// writes the translations once, as one History step, on Apply.
const { test, expect } = require('@playwright/test');
const { execFileSync } = require('node:child_process');
const { copyFileSync, readFileSync, writeFileSync } = require('node:fs');
const path = require('node:path');
const { evaluateWithRetainedPromise, generateAndWaitForResult } = require('./helpers/app-lifecycle.cjs');
const { showDrawing } = require('./helpers/drawing-bar.cjs');
const { readSessionFile } = require('./helpers/session-file.cjs');
const { loadSessionFile, openFresh, settle } = require('./helpers/audit-browser.cjs');

const root = path.resolve(process.env.GBDRAW_REPO || process.cwd());
const data = path.join(root, 'gbdraw', 'web', 'tutorial-data');
const INPUTS = [
  path.join(data, 'human-mitochondrion', 'HmmtDNA.gbk'),
  path.join(data, 'metazoan-mitochondria-four', 'NC_002333.2.gb')
];

const openTwoMitogenomes = async (page) => {
  await openFresh(page);
  await showDrawing(page, { kind: 'linear' });
  await page.evaluate(async (files) => {
    const app = window.__GBDRAW_APP__;
    if (app.linearSeqs.length < 2) app.addLinearSeq();
    files.forEach(({ name, text }, index) => app.setLinearSeqPrimaryFile(index, 'gb', new File([text], name, { type: 'text/plain' })));
    await window.Vue.nextTick();
  }, INPUTS.map((file) => ({ name: path.basename(file), text: readFileSync(file, 'utf8') })));
  await settle(page);
};

// The committed Linear request's record translations and display starts, and each drawn record's translation.
const placement = (page) => page.evaluate(async () => {
  const { getCommittedCanonicalRenderRequest } = await import('/gbdraw/web/js/services/config.js');
  const { state } = await import('/gbdraw/web/js/state.js');
  const committed = getCommittedCanonicalRenderRequest();
  const svg = window.__GBDRAW_APP__.svgContainer.querySelector('svg');
  return {
    committed: committed?.layout?.recordTranslations ?? null,
    records: committed?.records.map(({ recordKey, cardinality, display }) => ({ recordKey, cardinality,
      startCoordinate: display?.startCoordinate ?? null })) ?? null,
    plan: JSON.stringify(state.similarityAlignmentPlan.value),
    drawn: [...svg.querySelectorAll('[data-record-translation-x]')].map((group) => Number(group.getAttribute('data-record-translation-x')))
  };
});
// The drawn figure, for a live = Generate comparison.
const drawnSvg = (page) => page.evaluate(() => window.__GBDRAW_APP__.svgContainer.querySelector('svg').outerHTML);

// The x each record is drawn at; a record without a stored translation is at 0.
const xs = (translations, count) => Array.from({ length: count }, (_, index) => translations?.[index]?.x ?? 0);

// The CLI run of the Web's Source recipe with --anchor_feature, and the request its Session stores.
const cliRequest = async (recipe, testInfo) => {
  const directory = testInfo.outputPath('cli');
  execFileSync('mkdir', ['-p', directory]);
  for (const file of INPUTS) copyFileSync(file, path.join(directory, path.basename(file)));
  const session = path.join(directory, 'anchored.gbdraw');
  const args = recipe[0] === 'linear' ? recipe.slice(1) : recipe;
  execFileSync(process.env.PYTHON || 'python', [
    '-m', 'gbdraw.cli', 'linear', ...args, '--anchor_feature', 'gene=COX1', '--save_session', '--session_output', session
  ], { cwd: directory, stdio: 'pipe', env: { ...process.env, PYTHONPATH: root } });
  return (await readSessionFile(session)).drawings[0].renderRequest;
};

// Open the Record Layout controls, type the selector, and resolve it.
const resolveSelector = async (page, text) => {
  const selector = page.getByLabel('Align at feature', { exact: true });
  await selector.evaluate((element) => {
    for (let details = element.closest('details'); details; details = details.parentElement?.closest('details')) details.open = true;
  });
  await selector.fill(text);
  await page.getByTestId('linear-anchor-feature-resolve').click();
  const preview = page.getByTestId('linear-anchor-feature-preview');
  await expect(preview).toBeVisible({ timeout: 120_000 });
  return preview;
};

// The Display start field of a record, its File card's record options opened.
const displayStart = async (page, name) => {
  for (const disclosure of await page.locator('[data-linear-source-records], [data-linear-record-options]').all()) {
    await disclosure.evaluate((element) => {
      for (let details = element; details; details = details.parentElement?.closest('details')) details.open = true;
    });
  }
  return page.getByRole('spinbutton', { name: `Display start ${name}`, exact: true });
};

// PD-OI-113 scenario revision 3: Apply rotates complete circular records so that the feature starts at their left
// end; their Display start shows the new base, nothing moves sideways, and the next Generate draws the same figure.
test('Align at feature rotates both mitogenomes in one step, Display start shows the new bases, one Undo restores, and Generate draws the same figure', async ({ page }, testInfo) => {
  test.setTimeout(480_000);
  await openTwoMitogenomes(page);
  await generateAndWaitForResult(page);
  const before = await placement(page);
  const beforeSvg = await drawnSvg(page);
  const recipe = await page.evaluate(() => window.__GBDRAW_APP__.lastRunInfo?.sourceRecipe?.invocation?.args || []);
  const human = await displayStart(page, 'NC_012920.1 #1');
  const zebrafish = await displayStart(page, 'NC_002333.2 #1');
  const startsBefore = [await human.inputValue(), await zebrafish.inputValue()];

  const undoCount = () => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
  const stepsBefore = await undoCount();
  const preview = await resolveSelector(page, 'gene=COX1');
  await expect(preview).toContainText('#1 NC_012920.1 starts at base 5904 (rotated)');
  await expect(preview).toContainText('#2 NC_002333.2 starts at base 6425 (rotated)');
  // Nothing changes before Apply.
  expect(await placement(page)).toEqual(before);
  expect(await undoCount()).toBe(stepsBefore);

  await page.getByTestId('linear-anchor-feature-apply').click();
  await expect(page.getByTestId('linear-anchor-feature-notice')).toHaveText('Records aligned at gene=COX1.', { timeout: 180_000 });
  const after = await placement(page);
  const afterSvg = await drawnSvg(page);
  expect(await undoCount()).toBe(stepsBefore + 1);
  await expect(preview).toHaveCount(0);
  expect(after.records.map(({ startCoordinate }) => startCoordinate)).toEqual([5904, 6425]);
  // Rotated, no record moves sideways.
  expect(xs(after.committed, 2)).toEqual([0, 0]);
  expect(after.drawn).toEqual([0, 0]);
  await expect(human).toHaveValue('5904');
  await expect(zebrafish).toHaveValue('6425');
  // The CLI draws the same recipe with --anchor_feature with the same starts and places.
  expect(recipe.length, 'the Source recipe is available').toBeGreaterThan(0);
  const cli = await cliRequest(recipe, testInfo);
  expect(cli.records.map(({ display }) => display.startCoordinate)).toEqual([5904, 6425]);
  expect(xs(cli.layout.recordTranslations, 2)).toEqual([0, 0]);

  // One Undo puts back the starts, the positions and the Display start fields; the notice goes with them.
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
  await settle(page);
  expect(await placement(page)).toEqual(before);
  expect(await drawnSvg(page)).toBe(beforeSvg);
  await expect(human).toHaveValue(startsBefore[0]);
  await expect(zebrafish).toHaveValue(startsBefore[1]);
  await expect(page.getByTestId('linear-anchor-feature-notice')).toHaveCount(0);
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.redo());
  await settle(page);
  expect(await placement(page)).toEqual(after);
  await expect(page.getByTestId('linear-anchor-feature-notice')).toHaveText('Records aligned at gene=COX1.');

  // Generate draws what Apply drew: the starts are in the drafts, and the anchor is not applied again.
  await generateAndWaitForResult(page);
  expect(await placement(page)).toEqual(after);
  expect(await drawnSvg(page)).toBe(afterSvg);
});

// OV-2474: a loaded Session of one file with two records holds one File card for both records. Apply rotates each
// record by itself, one Undo restores, and the next Generate draws what Apply drew.
test('Align at feature on a loaded Session of a two-record file rotates each record, and Generate draws what Apply drew', async ({ page }, testInfo) => {
  test.setTimeout(480_000);
  const directory = testInfo.outputPath('cli');
  execFileSync('mkdir', ['-p', directory]);
  const file = path.join(directory, 'two-mitogenomes.gb');
  writeFileSync(file, INPUTS.map((input) => readFileSync(input, 'utf8')).join(''));
  const session = path.join(directory, 'two-mitogenomes.gbdraw');
  execFileSync(process.env.PYTHON || 'python', ['-m', 'gbdraw.cli', 'linear', '--gbk', file, '-o', 'two',
    '--save_session', '--session_output', session], { cwd: directory, stdio: 'pipe', env: { ...process.env, PYTHONPATH: root } });
  await openFresh(page);
  await loadSessionFile(page, session);
  const loaded = await placement(page);
  const loadedSvg = await drawnSvg(page);
  expect(loaded.records.map(({ startCoordinate }) => startCoordinate)).toEqual([null, null]);

  const preview = await resolveSelector(page, 'gene=COX1');
  await expect(preview).toContainText('#1 NC_012920.1 starts at base 5904 (rotated)');
  await expect(preview).toContainText('#2 NC_002333.2 starts at base 6425 (rotated)');
  await page.getByTestId('linear-anchor-feature-apply').click();
  await expect(page.getByTestId('linear-anchor-feature-notice')).toHaveText('Records aligned at gene=COX1.', { timeout: 180_000 });
  const after = await placement(page);
  const afterSvg = await drawnSvg(page);
  expect(after.records.map(({ startCoordinate }) => startCoordinate)).toEqual([5904, 6425]);
  expect(after.drawn).toEqual([0, 0]);

  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
  await settle(page);
  expect(await placement(page)).toEqual(loaded);
  expect(await drawnSvg(page)).toBe(loadedSvg);
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.redo());
  await settle(page);
  expect(await placement(page)).toEqual(after);

  await generateAndWaitForResult(page);
  const generated = await placement(page);
  expect(generated.records.map(({ startCoordinate }) => startCoordinate)).toEqual([5904, 6425]);
  expect(generated.drawn).toEqual(after.drawn);
  expect(await drawnSvg(page)).toBe(afterSvg);
});

// A Session with an active Similarity alignment: Apply clears the plan in its own step, keeps the
// aligned places as record translations, and one Undo brings back the plan and the translations.
test('Align at feature over an active Similarity alignment clears it in the same step, and one Undo restores both', async ({ page }) => {
  test.setTimeout(480_000);
  await openFresh(page);
  await loadSessionFile(page, path.join(root, 'gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw'));
  const before = await placement(page);
  expect(before.plan).not.toBe('null');
  const undoCount = () => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
  const stepsBefore = await undoCount();

  // Only BGC0000709 has this protein: it is the reference, every other record keeps its aligned place.
  const preview = await resolveSelector(page, 'protein_id=CAF33300.1');
  await expect(preview).toContainText('Apply clears the active similarity alignment.');
  await expect(preview).toContainText('#2 BGC0000709 keeps its place (first record with the feature)');
  await expect(preview).toContainText('#1 BGC0000708 keeps its place (no matching feature)');
  await page.getByTestId('linear-anchor-feature-apply').click();
  await expect(page.getByTestId('linear-anchor-feature-notice')).toHaveText('Records aligned at protein_id=CAF33300.1.', { timeout: 180_000 });
  await expect(page.locator('[data-similarity-alignment-notice]')).toHaveText('Alignment cleared: records aligned at a feature.');
  const after = await placement(page);
  expect(after.plan).toBe('null');
  expect(await undoCount()).toBe(stepsBefore + 1);
  // The records stay where the alignment drew them, now as their own translations.
  expect(after.drawn).toEqual(before.drawn);
  expect(xs(after.committed, after.drawn.length)).toEqual(after.drawn);

  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
  await settle(page);
  expect(await placement(page)).toEqual(before);
  await expect(page.locator('[data-similarity-alignment-notice]')).toHaveCount(0);
  await expect(page.getByTestId('linear-anchor-feature-notice')).toHaveCount(0);
});
