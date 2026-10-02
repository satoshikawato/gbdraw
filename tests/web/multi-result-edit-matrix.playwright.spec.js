// G-B (Web GUI audit 2026-09-30; R2, D-07/PD-OI-062): Result topology matrix.
// In a Circular batch (one Result per record), editor edits made while Result 2
// is displayed reach each Result's mounted SVG and its committed content (the
// export and Session source): after a Result 1 -> 2 round trip, after Generate,
// and after Save, a fresh Load and Generate. The drawer lists the displayed
// Result's record (FE-05). Scoped color and visibility reach Result 1; label
// edits stay with their feature (FE-01 to FE-03, PV-09). A Label TSV imported on
// Result 1 reaches Result 2 the same way, and Undo and Redo reach both Results
// as one step (B6, R3). The grid and Linear topologies are in
// live-edit-generate-equivalence (G-A).
const { test, expect } = require('@playwright/test');
const { generateAndWaitForResult } = require('./helpers/app-lifecycle.cjs');
const { editFeature, openBatch, settle } = require('./helpers/audit-browser.cjs');
const { download, load } = require('./helpers/mode-transition.cjs');

test.describe.configure({ retries: 0 });

const RED = '#ff0000';
const TEAL = '#2a9d8f';
const EDITS = [
  // "duplicate protein" annotates TESTA_0001/0002 and TESTB_0001/0002.
  ['annotation-label color scope', (page) => editFeature(page, 'TESTB_0001', { fill: RED, scope: 'annotationLabel' })],
  // "Exact product: spliced minus" covers TESTA_0004 and TESTB_0004.
  ['product visibility scope', (page) => editFeature(page, 'TESTB_0004', { visibility: 'off', visibilityScope: 'product' })],
  ['legend entry color', (page) => page.evaluate((color) => {
    const app = window.__GBDRAW_APP__;
    return app.updateLegendEntryColor(app.legendEntries.findIndex((entry) => entry.caption === 'tRNA'), color);
  }, TEAL)],
  ['label text', (page) => editFeature(page, 'TESTB_0003', { labelText: 'MATRIX_LABEL' })],
  ['label hidden', (page) => editFeature(page, 'TESTB_0005', { labelVisibility: 'off' })]
];

// Per Result: [record, edited (each must fail before the edits), kept]; a check
// is [name, read(state), holds(value)].
const shown = (color) => `${color}|shown`;
// A live hide keeps the glyph with display none; Generate leaves the feature out
// of the drawing and of the feature catalog (null: no element).
const notDrawn = (value) => value === null || /\|hidden$/.test(value || '');
const feature = (tag) => (s) => s.features[tag];
const label = (tag) => (s) => s.labels[tag];
const is = (expected) => (value) => value === expected;
const EXPECTED = [
  ['TESTA', [
    ['TESTA_0001 red', feature('TESTA_0001'), is(shown(RED))],
    ['TESTA_0002 red', feature('TESTA_0002'), is(shown(RED))],
    ['TESTA_0004 hidden', feature('TESTA_0004'), notDrawn],
    ['tRNA swatch teal', (s) => s.legend.tRNA, is(TEAL)]
  ], [
    ['TESTA_0003 label kept', label('TESTA_0003'), is('origin spanning protein')],
    ['TESTA_0005 label kept', label('TESTA_0005'), is('codon start two')]
  ]],
  ['TESTB', [
    ['TESTB_0001 red', feature('TESTB_0001'), is(shown(RED))],
    ['TESTB_0002 red', feature('TESTB_0002'), is(shown(RED))],
    ['TESTB_0004 hidden', feature('TESTB_0004'), notDrawn],
    ['tRNA swatch teal', (s) => s.legend.tRNA, is(TEAL)],
    ['TESTB_0003 label text', label('TESTB_0003'), is('MATRIX_LABEL')],
    ['TESTB_0005 label hidden', label('TESTB_0005'), is('')]
  ], [
    ['TESTB_0006 label kept', label('TESTB_0006'), is('gtg start')]
  ]]
];

const openDrawer = async (page) => {
  if (!await page.locator('.right-drawer').isVisible()) await page.locator('.drawer-toggle').click();
};

const readyCount = (page) => page.evaluate(() => [
  ...(window.__AUDIT_LIFECYCLE__ || []), ...(window.__MODE_EVENTS__ || []).map((event) => event.name)
].filter((name) => name === 'preview.ready-receipt-accepted').length);

const show = async (page, index) => {
  if (await page.evaluate(() => window.__GBDRAW_APP__.selectedResultIndex) !== index) {
    const before = await readyCount(page);
    await page.locator('h2 select').selectOption({ index });
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.selectedResultIndex)).toBe(index);
    await expect.poll(() => readyCount(page), { timeout: 60_000 }).toBeGreaterThan(before);
  }
  await settle(page);
};

// The displayed Result in the mounted SVG and in its committed content.
const resultState = (page, ids) => page.evaluate((known) => {
  const app = window.__GBDRAW_APP__;
  const index = app.selectedResultIndex;
  const hidden = (element) => {
    for (let node = element; node && node.nodeType === 1; node = node.parentNode) {
      if (node.getAttribute('display') === 'none' || /display\s*:\s*none/.test(node.getAttribute('style') || '')) return true;
    }
    return false;
  };
  const describe = (root) => {
    const features = {};
    const labels = {};
    for (const [tag, svgId] of Object.entries(known)) {
      const id = CSS.escape(svgId);
      const parts = [...root.querySelectorAll(`[data-gbdraw-feature-id="${id}"]`)];
      const fill = parts.map((element) => element.getAttribute('fill')).find((value) => value && value !== 'none');
      features[tag] = parts.length ? `${fill}|${parts.some(hidden) ? 'hidden' : 'shown'}` : null;
      labels[tag] = [...root.querySelectorAll(`[data-label-feature-id="${id}"] text, text[data-label-feature-id="${id}"]`)]
        .filter((element) => !hidden(element)).map((element) => element.textContent.trim()).filter(Boolean).join(' ');
    }
    const legend = Object.fromEntries([...root.querySelectorAll('g[data-legend-key]')].filter((entry) => !hidden(entry))
      .map((entry) => [entry.querySelector('text')?.textContent || entry.getAttribute('data-legend-key'),
        [...entry.querySelectorAll('path,rect')].map((element) => element.getAttribute('fill')).find((value) => value && value !== 'none')]));
    return { features, labels, legend };
  };
  return {
    index,
    drawer: [...new Set(app.visibleFeatureRows.map((row) => row.record_id))],
    mounted: describe(app.svgContainer.querySelector('svg')),
    content: describe(new DOMParser().parseFromString(app.results[index].content, 'image/svg+xml').documentElement)
  };
}, ids);

// Feature ids of the displayed Result, read before the edits (a hidden feature
// leaves the catalog after Generate).
const featureIds = (page) => page.evaluate(() => Object.fromEntries(window.__GBDRAW_APP__.extractedFeatures
  .filter((feature) => feature.locus_tag).map((feature) => [feature.locus_tag, feature.svg_id])));

// Failed checks of every Result, displaying each in turn.
const IDS = [];
const violations = async (page, stage, { edited = true } = {}) => {
  const failed = [];
  for (const [index, [record, edits, kept]] of EXPECTED.entries()) {
    await show(page, index);
    if (!edited) IDS[index] = await featureIds(page);
    const state = await resultState(page, IDS[index]);
    if (state.drawer.join() !== record) failed.push(`${stage} Result ${index + 1}: drawer lists ${state.drawer.join()}`);
    for (const source of ['mounted', 'content']) {
      for (const [name, read, holds] of [...(edited ? edits : []), ...kept]) {
        const value = read(state[source]);
        if (!holds(value)) failed.push(`${stage} Result ${index + 1} ${source}: ${name} (got ${JSON.stringify(value)})`);
      }
      if (!edited) {
        for (const [name, read, holds] of edits) {
          if (holds(read(state[source]))) failed.push(`${stage} Result ${index + 1} ${source}: ${name} holds before the edits`);
        }
      }
    }
  }
  return failed;
};

test('Circular batch: edits made on Result 2 reach both Results through selection, Generate, Save and Load', async ({ page, browser }, testInfo) => {
  test.setTimeout(600_000);
  await openBatch(page);
  await openDrawer(page);
  // Negative: every edited check fails on the unedited Results.
  expect(await violations(page, 'before edits', { edited: false })).toEqual([]);

  await show(page, 1);
  for (const [name, apply] of EDITS) {
    await test.step(name, async () => {
      await apply(page);
      await settle(page);
    });
  }
  await show(page, 0);
  expect.soft(await violations(page, 'Result 2 -> 1 -> 2')).toEqual([]);

  await generateAndWaitForResult(page);
  await settle(page);
  expect.soft(await violations(page, 'Generate')).toEqual([]);

  const saved = testInfo.outputPath('batch-matrix.gbdraw-session.json');
  await download(page, 'Save Session', saved);
  const reloaded = await load(browser, saved);
  try {
    await settle(reloaded);
    await generateAndWaitForResult(reloaded);
    await settle(reloaded);
    await openDrawer(reloaded);
    expect.soft(await violations(reloaded, 'Save, Load, Generate')).toEqual([]);
  } finally {
    await reloaded.context().close();
  }
});

// B6 (R2, R3): one Label TSV import while Result 1 is displayed. The rows select
// a Result 2 feature, a feature label on both Results, and (global `label` row)
// one source text on both Results.
const LABEL_TSV = [
  '*\tCDS\tlocus_tag\t^TESTB_0003$\tIMPORTED_B3',
  '*\tCDS\tlabel\t^codon start two$\tIMPORTED_CST',
  '*\t*\tlabel\t^gtg start$\tIMPORTED_GTG'
].join('\n') + '\n';
const ORIGINAL_LABELS = [
  [['TESTA_0003', 'origin spanning protein'], ['TESTA_0005', 'codon start two'], ['TESTA_0006', 'gtg start']],
  [['TESTB_0003', 'origin spanning protein'], ['TESTB_0005', 'codon start two'], ['TESTB_0006', 'gtg start']]
];
const IMPORTED_LABELS = [
  [['TESTA_0003', 'origin spanning protein'], ['TESTA_0005', 'IMPORTED_CST'], ['TESTA_0006', 'IMPORTED_GTG']],
  [['TESTB_0003', 'IMPORTED_B3'], ['TESTB_0005', 'IMPORTED_CST'], ['TESTB_0006', 'IMPORTED_GTG']]
];

// Labels that differ from the expected ones, displaying each Result in turn.
const labelViolations = async (page, stage, ids, expected) => {
  const failed = [];
  for (const [index, labels] of expected.entries()) {
    await show(page, index);
    const state = await resultState(page, ids);
    for (const source of ['mounted', 'content']) {
      for (const [tag, text] of labels) {
        const value = state[source].labels[tag];
        if (value !== text) failed.push(`${stage} Result ${index + 1} ${source}: ${tag} (got ${JSON.stringify(value)})`);
      }
    }
  }
  return failed;
};

test('Circular batch: a Label TSV imported on Result 1 reaches Result 2 live, through Undo and Redo, and after Generate', async ({ page }) => {
  test.setTimeout(600_000);
  await openBatch(page);
  const ids = await featureIds(page);
  expect(await labelViolations(page, 'before import', ids, ORIGINAL_LABELS)).toEqual([]);

  await show(page, 0);
  const undoCount = () => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
  const before = await undoCount();
  await page.evaluate((tsv) => window.__GBDRAW_APP__.loadLabelOverrideTable({
    target: { files: [new File([tsv], 'labels.tsv')], value: 'labels.tsv' }
  }), LABEL_TSV);
  await settle(page);
  expect(await undoCount()).toBe(before + 1);
  expect.soft(await labelViolations(page, 'import', ids, IMPORTED_LABELS)).toEqual([]);

  // Undo and Redo run while Result 2 is displayed; Result 1 follows when shown.
  await page.evaluate(() => window.__GBDRAW_HISTORY__.undo());
  await settle(page);
  expect(await undoCount()).toBe(before);
  expect.soft(await labelViolations(page, 'Undo', ids, ORIGINAL_LABELS)).toEqual([]);

  await page.evaluate(() => window.__GBDRAW_HISTORY__.redo());
  await settle(page);
  expect.soft(await labelViolations(page, 'Redo', ids, IMPORTED_LABELS)).toEqual([]);

  await generateAndWaitForResult(page);
  await settle(page);
  expect.soft(await labelViolations(page, 'Generate', ids, IMPORTED_LABELS)).toEqual([]);
});
