// Align at feature (Linear) asks Python for the record translations that
// `gbdraw linear --anchor_feature` draws with, previews which records move, and
// writes the translations once, as one History step, on Apply.
const { test, expect } = require('@playwright/test');
const { execFileSync } = require('node:child_process');
const { copyFileSync, readFileSync } = require('node:fs');
const path = require('node:path');
const { evaluateWithRetainedPromise, generateAndWaitForResult } = require('./helpers/app-lifecycle.cjs');
const { openFresh, settle } = require('./helpers/audit-browser.cjs');

const root = path.resolve(process.env.GBDRAW_REPO || process.cwd());
const data = path.join(root, 'gbdraw', 'web', 'tutorial-data');
const INPUTS = [
  path.join(data, 'human-mitochondrion', 'HmmtDNA.gbk'),
  path.join(data, 'metazoan-mitochondria-four', 'NC_002333.2.gb')
];

const openTwoMitogenomes = async (page) => {
  await openFresh(page);
  await page.evaluate(async (files) => {
    const app = window.__GBDRAW_APP__;
    app.selectDrawing('linear');
    await window.Vue.nextTick();
    if (app.linearSeqs.length < 2) app.addLinearSeq();
    files.forEach(({ name, text }, index) => app.setLinearSeqPrimaryFile(index, 'gb', new File([text], name, { type: 'text/plain' })));
    await window.Vue.nextTick();
  }, INPUTS.map((file) => ({ name: path.basename(file), text: readFileSync(file, 'utf8') })));
  await settle(page);
};

// The committed Linear request's record translations, and each drawn record's translation.
const placement = (page) => page.evaluate(async () => {
  const { getCommittedCanonicalRenderRequest } = await import('/gbdraw/web/js/services/config.js');
  const committed = getCommittedCanonicalRenderRequest();
  const svg = window.__GBDRAW_APP__.svgContainer.querySelector('svg');
  return {
    committed: committed?.layout?.recordTranslations ?? null,
    drawn: [...svg.querySelectorAll('[data-record-translation-x]')].map((group) => Number(group.getAttribute('data-record-translation-x')))
  };
});

// The x each record is drawn at; a record without a stored translation is at 0.
const xs = (translations, count) => Array.from({ length: count }, (_, index) => translations?.[index]?.x ?? 0);

// The CLI run of the Web's Source recipe with --anchor_feature, and the translations its Session stores.
const cliTranslations = (recipe, testInfo) => {
  const directory = testInfo.outputPath('cli');
  execFileSync('mkdir', ['-p', directory]);
  for (const file of INPUTS) copyFileSync(file, path.join(directory, path.basename(file)));
  const session = path.join(directory, 'anchored.gbdraw-session.json');
  const args = recipe[0] === 'linear' ? recipe.slice(1) : recipe;
  execFileSync(process.env.PYTHON || 'python', [
    '-m', 'gbdraw.cli', 'linear', ...args, '--anchor_feature', 'gene=COX1', '--save_session', '--session_output', session
  ], { cwd: directory, stdio: 'pipe', env: { ...process.env, PYTHONPATH: root } });
  return JSON.parse(readFileSync(session, 'utf8')).renderRequest.layout.recordTranslations;
};

test('Align at feature previews the moved record, writes its translation as one step, and Undo restores it', async ({ page }, testInfo) => {
  test.setTimeout(480_000);
  await openTwoMitogenomes(page);
  await generateAndWaitForResult(page);
  const before = await placement(page);
  const recipe = await page.evaluate(() => window.__GBDRAW_APP__.lastRunInfo?.sourceRecipe?.invocation?.args || []);

  const selector = page.getByLabel('Align at feature', { exact: true });
  await selector.evaluate((element) => {
    for (let details = element.closest('details'); details; details = details.parentElement?.closest('details')) details.open = true;
  });
  await selector.fill('gene=COX1');
  const undoCount = () => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
  const stepsBefore = await undoCount();
  await page.getByTestId('linear-anchor-feature-resolve').click();
  const preview = page.getByTestId('linear-anchor-feature-preview');
  await expect(preview).toContainText('#1 NC_012920.1: keeps its place (first record with the feature)', { timeout: 120_000 });
  await expect(preview).toContainText(/#2 NC_002333\.2: moves -\d/);
  // Nothing changes before Apply.
  expect(await placement(page)).toEqual(before);
  expect(await undoCount()).toBe(stepsBefore);

  await page.getByTestId('linear-anchor-feature-apply').click();
  await expect(page.getByTestId('linear-anchor-feature-notice')).toHaveText('Records aligned at gene=COX1.', { timeout: 180_000 });
  const after = await placement(page);
  expect(await undoCount()).toBe(stepsBefore + 1);
  expect(after.committed.map(({ recordKey }) => recordKey)).toHaveLength(2);
  const [human, zebrafish] = xs(after.committed, 2);
  const [humanBefore, zebrafishBefore] = xs(before.committed, 2);
  expect(human).toBe(humanBefore);
  await expect(preview).toHaveCount(0);
  expect(zebrafish).toBeLessThan(zebrafishBefore);
  expect(after.drawn).toEqual(xs(after.committed, 2));
  // The CLI draws the same recipe with --anchor_feature at the same places.
  expect(recipe.length, 'the Source recipe is available').toBeGreaterThan(0);
  expect(xs(cliTranslations(recipe, testInfo), 2)).toEqual([human, zebrafish]);

  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
  await settle(page);
  expect(await placement(page)).toEqual(before);
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.redo());
  await settle(page);
  expect(await placement(page)).toEqual(after);

  // Generate keeps the written translations: the anchor is not applied again.
  await generateAndWaitForResult(page);
  expect(await placement(page)).toEqual(after);
});
