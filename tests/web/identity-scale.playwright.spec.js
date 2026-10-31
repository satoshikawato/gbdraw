// A fixed identity colour scale (`--identity_scale MIN MAX`, the mode's
// `objects.<blast_match|conservation>.identity_scale` leaf) keeps the match
// colours and the Legend labels of a Linear comparison when the identity
// threshold changes, as the CLI does; a live palette recolour reads the same
// scale as Generate (R3).
const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const path = require('node:path');
const { generateAndWaitForResult } = require('./helpers/app-lifecycle.cjs');
const { openFresh, settle } = require('./helpers/audit-browser.cjs');
const { expectLiveEqualsGenerate, settleLive } = require('./helpers/live-generate-parity.cjs');
const { switchPalette } = require('./helpers/live-generate-parity-steps.cjs');

const root = path.resolve(process.env.GBDRAW_REPO || process.cwd());
const examples = path.join(root, 'examples');

const openLinearComparison = async (page) => {
  await openFresh(page);
  await page.evaluate(async ({ first, second, table }) => {
    const app = window.__GBDRAW_APP__;
    app.selectDrawing('linear');
    await window.Vue.nextTick();
    if (app.linearSeqs.length < 2) app.addLinearSeq();
    app.setLinearSeqPrimaryFile(0, 'gb', new File([first], 'MjeNMV.gb', { type: 'text/plain' }));
    app.setLinearSeqPrimaryFile(1, 'gb', new File([second], 'MelaMJNV.gb', { type: 'text/plain' }));
    app.linearComparisonPlan.mode = 'adjacent';
    app.linearComparisonPlan.defaultSource = 'upload';
    app.linearComparisonPlan.edges.splice(0, app.linearComparisonPlan.edges.length, {
      id: 'identity-scale-edge', queryUid: app.linearSeqs[0].uid, subjectUid: app.linearSeqs[1].uid,
      included: true, fileActive: true, losatFilenameActive: false, source: 'upload',
      file: new File([table], 'MjeNMV.MelaMJNV.tblastx.out', { type: 'text/plain' }), losatFilename: ''
    });
    app.form.show_labels_linear = 'none';
    await window.Vue.nextTick();
  }, {
    first: readFileSync(path.join(examples, 'MjeNMV.gb'), 'utf8'),
    second: readFileSync(path.join(examples, 'MelaMJNV.gb'), 'utf8'),
    table: readFileSync(path.join(examples, 'MjeNMV.MelaMJNV.tblastx.out'), 'utf8')
  });
  await settle(page);
};

// The fill of each drawn match by its identity, and the distinct percent labels in drawn order.
const comparisonView = (page) => page.evaluate(() => {
  const svg = window.__GBDRAW_APP__.svgContainer.querySelector('svg');
  const fills = {};
  for (const match of svg.querySelectorAll('path[data-identity]')) {
    fills[match.getAttribute('data-identity')] = match.getAttribute('fill');
  }
  const labels = [...svg.querySelectorAll('text')].map((text) => text.textContent.trim())
    .filter((text) => /^\d+(?:\.\d+)?%$/.test(text));
  return { fills, labels: [...new Set(labels)] };
});

const setThreshold = async (page, value) => {
  const input = page.getByLabel('Linear comparison minimum identity');
  await input.fill(String(value));
  await input.blur();
};

test('a fixed identity scale keeps match colours and Legend labels when the threshold changes', async ({ page }) => {
  test.setTimeout(420_000);
  await openLinearComparison(page);

  const minimum = page.getByLabel('Linear comparison identity scale minimum');
  await minimum.scrollIntoViewIfNeeded();
  await minimum.fill('30');
  await minimum.blur();
  // The blank maximum takes the top of the identity domain.
  expect(await page.evaluate(() => window.__GBDRAW_APP__.adv.identity_scale)).toEqual([30, 100]);

  await setThreshold(page, 30);
  await generateAndWaitForResult(page);
  await settleLive(page);
  const low = await comparisonView(page);
  const recipe = await page.evaluate(() => window.__GBDRAW_APP__.lastRunInfo?.sourceRecipe?.invocation?.args || []);
  expect(recipe.slice(recipe.indexOf('--identity_scale'), recipe.indexOf('--identity_scale') + 3))
    .toEqual(['--identity_scale', '30', '100']);

  await setThreshold(page, 60);
  await generateAndWaitForResult(page);
  await settleLive(page);
  const high = await comparisonView(page);

  expect(low.labels).toEqual(['30%', '100%']);
  expect(high.labels).toEqual(low.labels);
  const kept = Object.keys(high.fills);
  expect(kept.length).toBeGreaterThan(0);
  expect(Object.keys(low.fills).length).toBeGreaterThan(kept.length);
  for (const identity of kept) expect(high.fills[identity], `identity ${identity}`).toBe(low.fills[identity]);

  // Without the scale the Legend starts at the threshold again.
  await minimum.fill('');
  await minimum.blur();
  expect(await page.evaluate(() => window.__GBDRAW_APP__.adv.identity_scale)).toBeNull();
  await generateAndWaitForResult(page);
  await settleLive(page);
  expect((await comparisonView(page)).labels).toEqual(['60%', '100%']);
});

test('a live palette recolour of a fixed-scale comparison equals Generate', async ({ page }) => {
  test.setTimeout(420_000);
  await openLinearComparison(page);
  await page.evaluate(() => { window.__GBDRAW_APP__.adv.identity_scale = [20, 90]; });
  await generateAndWaitForResult(page);
  await settleLive(page);
  await switchPalette(page, 'forest');
  await settleLive(page);
  const live = await comparisonView(page);
  await expectLiveEqualsGenerate(page, { label: 'palette recolour with a fixed identity scale' });
  const generated = await comparisonView(page);
  expect(live).toEqual(generated);
  expect(generated.labels).toEqual(['20%', '90%']);
});
