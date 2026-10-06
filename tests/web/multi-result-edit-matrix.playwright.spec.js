// G-B (Web GUI audit 2026-09-30; R2, D-07/PD-OI-062): Result topology matrix.
// In a Circular batch (one Result per record), editor edits made while Result 2
// is displayed reach each Result's mounted SVG and its committed content (the
// export and Session source): after a Result 1 -> 2 round trip, after Generate,
// and after Save, a fresh Load and Generate. The drawer lists the displayed
// Result's record (FE-05). Scoped color and visibility reach Result 1; label
// edits stay with their feature (FE-01 to FE-03, PV-09). A Label TSV imported on
// Result 1 reaches Result 2 the same way, and Undo and Redo reach both Results
// as one step (B6, R3). A Legend sort keeps the entries only one Result draws
// in place through a Result switch, which records no Undo step (B18); Undo and
// Redo of that sort and Sort by default reach each Result without copying
// those entries to another Result (B19, B20). Every
// editor edit kind commits the displayed Result at once, so none waits for a
// Result switch (R1), and Undo and Redo of a drag restore the Result it was
// made on while another Result is shown (B17). The grid and Linear topologies
// are in live-edit-generate-equivalence (G-A).
const { test, expect } = require('@playwright/test');
const { evaluateWithRetainedPromise, generateAndWaitForResult } = require('./helpers/app-lifecycle.cjs');
const { editFeature, legendCaptions, openBatch, settle } = require('./helpers/audit-browser.cjs');
const { download, load } = require('./helpers/mode-transition.cjs');

test.describe.configure({ retries: 0 });

const RED = '#ff0000';
const TEAL = '#2a9d8f';
const EDITS = [
  // "duplicate protein" annotates TESTA_0001/0002 and TESTB_0001/0002.
  ['annotation-label color scope', (page) => editFeature(page, 'TESTB_0001', { fill: RED, scope: 'annotationLabel' })],
  // "Exact product: spliced minus" covers TESTA_0004 and TESTB_0004.
  ['product visibility scope', (page) => editFeature(page, 'TESTB_0004', { visibility: 'off', visibilityScope: 'product' })],
  ['legend entry color', (page) => evaluateWithRetainedPromise(page, (color) => {
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
  await evaluateWithRetainedPromise(page, (tsv) => window.__GBDRAW_APP__.loadLabelOverrideTable({
    target: { files: [new File([tsv], 'labels.tsv')], value: 'labels.tsv' }
  }), LABEL_TSV);
  await settle(page);
  expect(await undoCount()).toBe(before + 1);
  expect.soft(await labelViolations(page, 'import', ids, IMPORTED_LABELS)).toEqual([]);

  // Undo and Redo run while Result 2 is displayed; Result 1 follows when shown.
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
  await settle(page);
  expect(await undoCount()).toBe(before);
  expect.soft(await labelViolations(page, 'Undo', ids, ORIGINAL_LABELS)).toEqual([]);

  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.redo());
  await settle(page);
  expect.soft(await labelViolations(page, 'Redo', ids, IMPORTED_LABELS)).toEqual([]);

  await generateAndWaitForResult(page);
  await settle(page);
  expect.soft(await labelViolations(page, 'Generate', ids, IMPORTED_LABELS)).toEqual([]);
});

// The entries only Result 2 draws once the "This feature only" color of
// TESTB_0001 is a rule: the rule row and the "other proteins" row that Python
// captions the remaining CDS row with (the rerender draws them, OV-43).
const RESULT_2_ONLY = ['duplicate protein', 'other proteins'];

// B18 (D-07, D-08, R11): a Legend Sort Z-A made on Result 2 keeps its order on
// Result 2 through a switch to Result 1 and back, also for the entry only
// Result 2 draws (the "This feature only" color entry of TESTB_0001), and
// Generate gives the same order. The Result picker is navigation, so a switch
// records no Undo step.
test('Circular batch: a Legend sort on Result 2 keeps its order through a Result switch and Generate', async ({ page }) => {
  test.setTimeout(300_000);
  await openBatch(page);
  await show(page, 1);
  await editFeature(page, 'TESTB_0001', { fill: RED });
  await settle(page);
  const undoCount = () => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
  const edited = await undoCount();
  const unsorted = await legendCaptions(page);
  await openDrawer(page);
  await page.locator('.right-drawer').getByRole('button', { name: 'Legend' }).click();
  await page.locator('.right-drawer').getByTitle('Sort Z-A', { exact: true }).click();
  await settle(page);
  const sorted = await legendCaptions(page);
  expect(sorted).toContain('duplicate protein');
  expect(sorted).not.toEqual(unsorted);
  expect(await legendCaptions(page, { source: 'content' })).toEqual(sorted);
  expect(await undoCount()).toBe(edited + 1);

  await show(page, 0);
  expect.soft(await legendCaptions(page), 'Result 1 shows the sort without the Result 2 entries, its own CDS entry after them')
    .toEqual(sorted.filter((caption) => !RESULT_2_ONLY.includes(caption)).concat('CDS'));
  expect.soft(await undoCount(), 'showing Result 1 records no Undo step').toBe(edited + 1);
  await show(page, 1);
  expect.soft(await legendCaptions(page), 'Result 2 keeps its sorted order').toEqual(sorted);
  expect.soft(await legendCaptions(page, { source: 'content' }), 'Result 2 content keeps its sorted order').toEqual(sorted);
  expect.soft(await undoCount(), 'showing Result 2 records no Undo step').toBe(edited + 1);

  // The rerender drew the "This feature only" color as a rule and captioned the
  // remaining CDS row "other proteins" (OV-43), as Generate does, so Generate
  // keeps the order.
  await generateAndWaitForResult(page);
  await settle(page);
  await show(page, 1);
  expect.soft(await legendCaptions(page), 'Generate keeps the sorted order on Result 2').toEqual(sorted);
});

// The displayed Result's Legend: mounted, committed (export and Session
// source), and listed in the Legend panel.
const legendView = async (page) => ({
  mounted: await legendCaptions(page),
  content: await legendCaptions(page, { source: 'content' }),
  listed: await page.evaluate(() => window.__GBDRAW_APP__.legendEntries.map((entry) => entry.caption))
});
const sameView = (captions) => ({ mounted: captions, content: captions, listed: captions });

// Result 2 gets an entry Result 1 does not draw (the "This feature only" color
// entry of TESTB_0001) and a Legend Sort Z-A.
const sortResult2Descending = async (page) => {
  await openBatch(page);
  await show(page, 0);
  const result1 = await legendCaptions(page);
  await show(page, 1);
  await editFeature(page, 'TESTB_0001', { fill: RED });
  await settle(page);
  const result2 = await legendCaptions(page);
  expect(result2.filter((caption) => !result1.includes(caption))).toEqual(RESULT_2_ONLY);
  await openDrawer(page);
  await page.locator('.right-drawer').getByRole('button', { name: 'Legend' }).click();
  await page.locator('.right-drawer').getByTitle('Sort Z-A', { exact: true }).click();
  await settle(page);
  const sorted = await legendCaptions(page);
  expect(sorted).not.toEqual(result2);
  const shared = (captions) => captions.filter((caption) => !RESULT_2_ONLY.includes(caption));
  return {
    result1,
    result2,
    sorted,
    shared,
    onResult1: (captions) => shared(captions).filter((caption) => caption !== 'CDS').concat('CDS')
  };
};

// OV-47: the default Legend order (`originalLegendOrder`) is one inventory,
// taken from the Result displayed when Generate or the automatic rerender
// draws. A rerender while Result 2 is shown drops Result 1's own rows from it,
// so Undo and Sort by default on Result 1 do not return Python's order. A
// per-Result inventory fixes it.
const OV_47 = 'OV-47: the default Legend order is one inventory, taken from the displayed Result; Undo and Sort by default do not return Python\'s order of the other Result';

// B19 (D-07, D-08, R3, R11): Undo and Redo of a Legend sort made on Result 2
// while Result 1 is displayed restore the order on Result 1 and never give it
// the entry only Result 2 draws; Result 2 shows the restored order once it is
// displayed, its own entry following the shared entries.
test('Circular batch: Undo and Redo of a Legend sort on Result 2 keep each Result\'s own entries while Result 1 is shown', async ({ page }) => {
  test.fail(true, OV_47);
  test.setTimeout(300_000);
  const { result1, result2, sorted, shared, onResult1 } = await sortResult2Descending(page);
  await show(page, 0);
  expect(await legendCaptions(page)).toEqual(onResult1(sorted));

  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
  await settle(page);
  expect.soft(await legendView(page), 'Undo on Result 1 gives its order before the sort, without the Result 2 entries')
    .toEqual(sameView(result1));
  await show(page, 1);
  expect.soft(await legendView(page), 'after Undo, Result 2 shows its order before the sort').toEqual(sameView(result2));

  await show(page, 0);
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.redo());
  await settle(page);
  expect.soft(await legendView(page), 'Redo on Result 1 sorts it again, without the Result 2 entry')
    .toEqual(sameView(onResult1(sorted)));
  await show(page, 1);
  const redone = await legendView(page);
  expect.soft(shared(redone.mounted), 'after Redo, Result 2 shows the sort').toEqual(shared(sorted));
  expect.soft(redone.mounted.filter((caption) => caption === 'duplicate protein'), 'Result 2 keeps its own entry')
    .toEqual(['duplicate protein']);
  expect.soft(redone.content, 'Result 2 content matches its preview').toEqual(redone.mounted);
  expect.soft(redone.listed, 'the Legend panel lists Result 2').toEqual(redone.mounted);
});

// B20 (D-07, D-08): Sort by default made on Result 1 reaches Result 2, which
// shows the Sort Z-A made on it, once Result 2 is displayed.
test('Circular batch: Sort by default on Result 1 reaches Result 2 shown in a sorted order', async ({ page }) => {
  test.fail(true, OV_47);
  test.setTimeout(300_000);
  const { result1, result2 } = await sortResult2Descending(page);
  await show(page, 0);
  await openDrawer(page);
  await page.locator('.right-drawer').getByTitle('Sort by default', { exact: true }).click();
  await settle(page);
  expect(await legendCaptions(page)).toEqual(result1);
  await show(page, 1);
  expect.soft(await legendView(page), 'Result 2 shows the default order').toEqual(sameView(result2));
});

// R1 (A1): each editor edit kind writes the displayed Result when it is made,
// so nothing stays pending for a Result switch to persist. Right after the
// edit, the committed content changed and equals the mounted SVG; after a
// switch to Result 1 and back, Result 2 is not reverted and still matches its
// mounted SVG. (Displaying a Result projects the editor intent onto it, D-07,
// so its bytes may change on the way back.)
const committed = (page) => page.evaluate(async () => {
  const { serializeCleanSvg } = await import('/gbdraw/web/js/services/svg-serialization.js');
  const app = window.__GBDRAW_APP__;
  const index = app.selectedResultIndex;
  const content = app.results[index].content;
  const mounted = serializeCleanSvg(app.svgContainer.querySelector('svg'));
  return { index, content, mounted, matchesMounted: mounted === content };
});
// Composition offsets of every Result's committed content.
const offsets = (page) => page.evaluate(async () => {
  const { compositionUserDeltas } = await import('/gbdraw/web/js/app/legend-layout/composition-actions.js');
  return window.__GBDRAW_APP__.results.map((result) => compositionUserDeltas(
    new DOMParser().parseFromString(result.content, 'image/svg+xml').documentElement
  ));
});
const mountedOffsets = (page) => page.evaluate(async () => {
  const { compositionUserDeltas } = await import('/gbdraw/web/js/app/legend-layout/composition-actions.js');
  return compositionUserDeltas(window.__GBDRAW_APP__.svgContainer.querySelector('svg'));
});
// Where the committed content and the mounted SVG first differ.
const firstDiff = ({ content, mounted }) => {
  let at = 0;
  while (at < content.length && content[at] === mounted[at]) at += 1;
  return `at ${at}: ${JSON.stringify(content.slice(at, at + 120))} vs ${JSON.stringify(mounted.slice(at, at + 120))}`;
};

const dragRole = (role, dx, dy) => (page) => page.evaluate(async ({ targetRole, deltaX, deltaY }) => {
  const app = window.__GBDRAW_APP__;
  app.layoutRepositionMode = true;
  await window.Vue.nextTick();
  const svg = app.svgContainer.querySelector('svg');
  const candidates = [...svg.querySelectorAll(`[data-gbdraw-composition-role="${targetRole}"]`)];
  const target = candidates.find((element) => element.id !== 'length_bar') || candidates[0];
  if (!target) throw new Error(`no ${targetRole} drag target`);
  const bounds = target.getBoundingClientRect();
  const x = bounds.left + Math.min(bounds.width / 2, 8);
  const y = bounds.top + Math.min(bounds.height / 2, 8);
  const event = (type, clientX, clientY, buttons) => new MouseEvent(type, {
    bubbles: true, cancelable: true, clientX, clientY, buttons, view: window
  });
  const frame = () => new Promise((resolve) => requestAnimationFrame(resolve));
  target.dispatchEvent(event('mousedown', x, y, 1));
  await frame();
  const moveTarget = targetRole === 'legend' ? svg : document;
  moveTarget.dispatchEvent(event('mousemove', x + deltaX, y + deltaY, 1));
  await frame();
  moveTarget.dispatchEvent(event('mouseup', x + deltaX, y + deltaY, 0));
  await new Promise((resolve) => setTimeout(resolve, 100));
}, { targetRole: role, deltaX: dx, deltaY: dy });

const onApp = (action, arg = null) => (page) => evaluateWithRetainedPromise(page, action, arg);
// A fourth element marks an edit that changes a Legend source, so the
// automatic rerender draws the Result again (OV-43): its content is Python's
// drawing, as after Generate, and shows the edit's text instead of equaling
// the mounted SVG that the preview binder completes.
const SWITCH_EDITS = [
  ['feature fill', (page) => editFeature(page, 'TESTB_0001', { fill: RED }), null, RED],
  ['feature stroke', onApp(async () => {
    const app = window.__GBDRAW_APP__;
    await app.openFeatureEditorFromList(app.extractedFeatures.find((item) => item.locus_tag === 'TESTB_0002'), null);
    await app.updateClickedFeatureStroke('#1d3557', 2.5);
    app.clickedFeature = null;
  })],
  ['feature visibility', (page) => editFeature(page, 'TESTB_0004', { visibility: 'off' })],
  ['label text', (page) => editFeature(page, 'TESTB_0003', { labelText: 'SWITCH_LABEL' })],
  ['legend entry color', onApp((color) => {
    const app = window.__GBDRAW_APP__;
    return app.updateLegendEntryColor(app.legendEntries.findIndex((entry) => entry.caption === 'tRNA'), color);
  }, TEAL)],
  ['legend stroke', onApp(() => {
    const app = window.__GBDRAW_APP__;
    return app.updateLegendEntryStrokeWidth(app.legendEntries.findIndex((entry) => entry.caption === 'tRNA'), 2);
  })],
  ['legend sort', onApp(() => window.__GBDRAW_APP__.sortLegendEntries('desc'))],
  ['legend drag', dragRole('legend', 30, 20)],
  ['diagram drag', dragRole('primary', 25, 15)],
  // A caption change resizes the Legend, which repositions it in place.
  ['legend rename', onApp(async () => {
    const app = window.__GBDRAW_APP__;
    await app.renameLegendEntry(app.legendEntries.findIndex((entry) => entry.caption === 'tRNA'), 'transfer RNA');
  }), null, 'transfer RNA'],
  ['reset positions', onApp(() => window.__GBDRAW_APP__.resetAllPositions())],
  ['canvas padding', onApp(() => { window.__GBDRAW_APP__.canvasPadding.right = 40; })],
  ['palette color', onApp(() => {
    const app = window.__GBDRAW_APP__;
    app.paletteInstantPreviewEnabled = true;
    app.currentColors.CDS = '#457b9d';
  })],
  ['track visibility', onApp(() => { window.__GBDRAW_APP__.form.suppress_gc = true; })]
];

test('Circular batch: each editor edit kind is committed when made, not by the Result switch', async ({ page }) => {
  test.setTimeout(600_000);
  await openBatch(page);
  await show(page, 1);
  const failed = [];
  for (const [name, apply, , rerendered] of SWITCH_EDITS) {
    await test.step(name, async () => {
      const before = await committed(page);
      await apply(page);
      await settle(page);
      const edited = await committed(page);
      if (edited.index !== 1) failed.push(`${name}: displays Result ${edited.index + 1}`);
      if (edited.content === before.content) failed.push(`${name}: Result 2 content unchanged after the edit`);
      if (rerendered && !edited.content.includes(rerendered)) failed.push(`${name}: Result 2 content lacks ${rerendered}`);
      if (!rerendered && !edited.matchesMounted) failed.push(`${name}: Result 2 content differs from the mounted SVG ${firstDiff(edited)}`);
      await show(page, 0);
      await show(page, 1);
      const back = await committed(page);
      if (back.content === before.content) failed.push(`${name}: Result 2 reverted by the switch`);
      if (rerendered && !back.content.includes(rerendered)) failed.push(`${name}: after the switch, Result 2 content lacks ${rerendered}`);
      if (!rerendered && !back.matchesMounted) failed.push(`${name}: after the switch, Result 2 differs from the mounted SVG ${firstDiff(back)}`);
    });
  }
  expect(failed).toEqual([]);

  // History reconciles the displayed Result in place: an Undo of a Result 2
  // drag made while Result 1 is shown leaves each Result matching its mounted SVG.
  await dragRole('primary', -20, 10)(page);
  await settle(page);
  await show(page, 0);
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
  await settle(page);
  for (const index of [0, 1]) {
    await show(page, index);
    const state = await committed(page);
    expect.soft(state.matchesMounted, `after Undo, Result ${index + 1}: ${firstDiff(state)}`).toBe(true);
  }
});

// B17 (D-07, R11): History restores a drag on the Result it was made on. A
// Result switch records no step. Undo made while Result 1 is shown keeps
// Result 1 shown and unmoved and restores Result 2's committed content (the
// Session and export source) at once; Redo drags Result 2 again. Each Result
// keeps matching its mounted SVG.
test('Circular batch: Undo and Redo of drags on Result 2 restore Result 2 while Result 1 is shown', async ({ page }) => {
  test.setTimeout(300_000);
  await openBatch(page);
  await show(page, 1);
  const undoCount = () => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
  const undragged = await offsets(page);
  await dragRole('primary', -20, 10)(page);
  await settle(page);
  await dragRole('legend', 30, 20)(page);
  await settle(page);
  const dragged = await committed(page);
  const draggedOffsets = await offsets(page);
  expect(draggedOffsets[0]).toEqual(undragged[0]);
  expect(draggedOffsets[1].primary).not.toEqual(undragged[1].primary);
  expect(draggedOffsets[1].legend).not.toEqual(undragged[1].legend);
  const steps = await undoCount();
  expect(steps).toBeLessThan(30);
  await show(page, 0);
  expect.soft(await undoCount(), 'a Result switch records no Undo step').toBe(steps);

  for (const _drag of ['legend', 'primary']) await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
  await settle(page);
  expect.soft(await page.evaluate(() => window.__GBDRAW_APP__.selectedResultIndex), 'Undo keeps Result 1 shown').toBe(0);
  expect.soft(await offsets(page), 'Undo restores Result 2 and leaves Result 1').toEqual(undragged);
  await show(page, 1);
  const undone = await committed(page);
  expect.soft(undone.content, 'Undo reverts the Result 2 drag').not.toBe(dragged.content);
  expect.soft(undone.matchesMounted, `after Undo, Result 2: ${firstDiff(undone)}`).toBe(true);
  // Result 1 has no edit, so its content keeps the admitted bytes (zero fast
  // path); its mounted root shows the offsets of that content.
  await show(page, 0);
  expect.soft(await mountedOffsets(page), 'after Undo, Result 1 is shown unmoved').toEqual(undragged[0]);

  for (const _drag of ['primary', 'legend']) await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.redo());
  await settle(page);
  expect.soft(await offsets(page), 'Redo drags Result 2 again and leaves Result 1').toEqual(draggedOffsets);
  await show(page, 1);
  const redone = await committed(page);
  expect.soft(redone.matchesMounted, `after Redo, Result 2: ${firstDiff(redone)}`).toBe(true);
});
