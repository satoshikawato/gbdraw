// A fixed identity colour scale (`--identity_scale MIN MAX`, the mode's
// `objects.<blast_match|conservation>.identity_scale` leaf) keeps the match
// colours (Linear) and ring colours (Circular) and the Legend labels when the
// identity threshold changes, as the CLI does; a live palette recolour reads
// the same scale as Generate (R3).
const { test, expect } = require('@playwright/test');
const { execFileSync } = require('node:child_process');
const { readFileSync, writeFileSync } = require('node:fs');
const path = require('node:path');
const { generateAndWaitForResult } = require('./helpers/app-lifecycle.cjs');
const { openFresh, openWithGenBank, settle } = require('./helpers/audit-browser.cjs');
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

// The comparison inputs sit in closed <details> sections; open them, then type and commit as a user does.
const setInput = async (page, label, value) => {
  const input = page.getByLabel(label, { exact: true });
  await input.evaluate((element) => {
    for (let details = element.closest('details'); details; details = details.parentElement?.closest('details')) {
      details.open = true;
    }
  });
  await input.fill(String(value));
  await input.dispatchEvent('change');
};

const setThreshold = (page, value) => setInput(page, 'Linear comparison minimum identity', value);
const setScaleMinimum = (page, value) => setInput(page, 'Linear comparison identity scale minimum', value);

// Each channel of a #rrggbb colour, for a comparison that allows the 1-unit rounding difference
// between the live recolour (rounds) and Python (truncates), OV-2393.
const channels = (color) => [1, 3, 5].map((start) => Number.parseInt(color.slice(start, start + 2), 16));

test('a fixed identity scale keeps match colours and Legend labels when the threshold changes', async ({ page }) => {
  test.setTimeout(420_000);
  await openLinearComparison(page);

  // Without a minimum there is no scale, so the maximum waits for one.
  const maximum = page.getByLabel('Linear comparison identity scale maximum', { exact: true });
  await expect(maximum).toBeDisabled();
  await setScaleMinimum(page, 30);
  await expect(maximum).toBeEnabled();
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
  await setScaleMinimum(page, '');
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
  expect(generated.labels).toEqual(['20%', '90%']);
  expect(live.labels).toEqual(generated.labels);
  expect(Object.keys(live.fills).sort()).toEqual(Object.keys(generated.fills).sort());
  // The live recolour reads the factor Python wrote for the fixed scale; only channel rounding may differ.
  for (const [identity, fill] of Object.entries(generated.fills)) {
    const differences = channels(live.fills[identity]).map((value, index) => Math.abs(value - channels(fill)[index]));
    expect(Math.max(...differences), `identity ${identity}: live ${live.fills[identity]}, Generate ${fill}`)
      .toBeLessThanOrEqual(1);
  }
});

// Three ring hits on the single-record fixture (subject FORCEDLBL), at 40, 60, and 95% identity.
const RING_FIXTURE = 'tests/fixtures/forced_label_underlay.gb';
const RING_BLAST = [40, 60, 95].map((identity, index) => [
  'RING', 'FORCEDLBL', identity.toFixed(1), '300', '0', '0', String(1 + 300 * index), String(300 * (index + 1)),
  String(1 + 300 * index), String(300 * (index + 1)), '1e-50', '200'
].join('\t')).join('\n') + '\n';

const drawRingsWithPython = (testInfo, threshold) => {
  const table = testInfo.outputPath('ring.blast');
  writeFileSync(table, RING_BLAST);
  execFileSync(process.env.PYTHON || 'python', [
    '-m', 'gbdraw.cli', 'circular', '--gbk', path.resolve(root, RING_FIXTURE), '--conservation_blast', table,
    '--identity', String(threshold), '--identity_scale', '20', '90', '-o', testInfo.outputPath(`rings-${threshold}`), '-f', 'svg'
  ], { cwd: root, stdio: 'pipe' });
  return readFileSync(testInfo.outputPath(`rings-${threshold}.svg`), 'utf8');
};

// What comparisonView reads, from an SVG text drawn by Python.
const svgView = (page, text) => page.evaluate((svgText) => {
  const svg = new DOMParser().parseFromString(svgText, 'image/svg+xml').documentElement;
  const fills = {};
  for (const match of svg.querySelectorAll('path[data-identity]')) fills[match.getAttribute('data-identity')] = match.getAttribute('fill');
  const labels = [...svg.querySelectorAll('text')].map((node) => node.textContent.trim())
    .filter((label) => /^\d+(?:\.\d+)?%$/.test(label));
  return { fills, labels: [...new Set(labels)] };
}, text);

test('a fixed identity scale keeps Circular ring colours and Legend labels, with the CLI Legend', async ({ page }, testInfo) => {
  test.setTimeout(420_000);
  await openWithGenBank(page, path.join(root, RING_FIXTURE));
  await page.evaluate(async (table) => {
    const app = window.__GBDRAW_APP__;
    app.circularConservation.source = 'upload';
    app.setCircularConservationUploadFiles([new File([table], 'ring.blast', { type: 'text/plain' })]);
    await window.Vue.nextTick();
  }, RING_BLAST);
  await settle(page);
  const maximum = page.getByLabel('Circular comparison identity scale maximum', { exact: true });
  await expect(maximum).toBeDisabled();
  await setInput(page, 'Circular comparison identity scale minimum', 20);
  await setInput(page, 'Circular comparison identity scale maximum', 90);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.adv.identity_scale)).toEqual([20, 90]);

  const views = {};
  for (const threshold of [30, 50]) {
    await setInput(page, 'Circular comparison minimum identity', threshold);
    await generateAndWaitForResult(page);
    await settleLive(page);
    views[threshold] = await comparisonView(page);
    // The Web gives each ring its own series colour, so the CLI comparison reads the
    // Legend labels and the drawn hits; tests/test_identity_scale.py pins the CLI colours.
    const cli = await svgView(page, drawRingsWithPython(testInfo, threshold));
    expect(views[threshold].labels, `Web Legend equals CLI at --identity ${threshold}`).toEqual(cli.labels);
    expect(Object.keys(views[threshold].fills).sort(), `Web hits equal CLI at --identity ${threshold}`)
      .toEqual(Object.keys(cli.fills).sort());
  }
  expect(views[30].labels).toEqual(['20%', '90%']);
  expect(views[50].labels).toEqual(views[30].labels);
  expect(Object.keys(views[30].fills).sort()).toEqual(['40.0', '60.0', '95.0']);
  expect(Object.keys(views[50].fills).sort()).toEqual(['60.0', '95.0']);
  for (const identity of ['60.0', '95.0']) expect(views[50].fills[identity], `identity ${identity}`).toBe(views[30].fills[identity]);
  await expectLiveEqualsGenerate(page, { label: 'Circular rings with a fixed identity scale' });
});
