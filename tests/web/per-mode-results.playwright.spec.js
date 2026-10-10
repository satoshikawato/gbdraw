// E1: each diagram mode keeps its own Result (Owner, 2026-10-07). Circular and
// Linear each hold one Result slot. A Generate replaces only its own mode's
// Result, a switch shows the arriving mode's Result (or the empty Preview), and
// Save writes every Result: the shown mode's set at the top level and the other
// mode's set in `otherModeResult` (Session 46). PD-OI-045, PD-OI-052,
// PD-OI-066, OIPC-C06, OIPC-C07.
const { test, expect } = require('@playwright/test');
const { execFileSync } = require('node:child_process');
const { mkdirSync, readFileSync } = require('node:fs');
const { gunzipSync } = require('node:zlib');
const { evaluateWithRetainedPromise, generateAndWaitForResult, openApp } = require('./helpers/app-lifecycle.cjs');
const { openWithGenBank } = require('./helpers/audit-browser.cjs');
const { semanticSnapshot, settleLive } = require('./helpers/live-generate-parity.cjs');
const { popup, closeEditor, loadEditorLegendRows } = require('./helpers/mode-transition.cjs');

test.describe.configure({ retries: 0 });

// History is last in, first out, so an artifact handle or checkpoint restore
// finds the mode it was captured in. A mismatch is counted
// (`artifactRestoreModeMismatchCount`) and must stay 0 on every page.
const trackModeMismatches = (target) => target.addInitScript(() => {
  window.__E1_MODE_MISMATCHES__ = 0;
  window.__E1_READINESS__ = [];
  window.__GBDRAW_TEST_HOOKS__ = {
    ...(window.__GBDRAW_TEST_HOOKS__ || {}),
    onStructuralMetric: (metric) => {
      if (metric.name === 'artifactRestoreModeMismatchCount') window.__E1_MODE_MISMATCHES__ += metric.value;
    },
    // The phases of the Results presented for display.
    onSessionLifecycleEvent: (event) => {
      if (event.name === 'preview.readiness-expectation-registered') window.__E1_READINESS__.push(event.phase);
    }
  };
});
const expectNoModeMismatch = async (target) => {
  expect(await target.evaluate(() => window.__E1_MODE_MISMATCHES__), 'artifactRestoreModeMismatchCount').toBe(0);
};
test.beforeEach(({ page }) => trackModeMismatches(page));
test.afterEach(async ({ page }) => {
  if (!page.isClosed()) await expectNoModeMismatch(page);
});
// A page in a fresh browser context, as a reader opening the saved file.
const freshPage = async (browser) => {
  const context = await browser.newContext({
    baseURL: `http://127.0.0.1:${process.env.GBDRAW_WEB_TEST_PORT || 4173}`, acceptDownloads: true
  });
  await trackModeMismatches(context);
  return context.newPage();
};

const FIXTURE = 'tests/fixtures/forced_label_underlay.gb';
const MODE_NAMES = { circular: 'Circular', linear: 'Linear' };
const EMPTY_TEXT = (mode) => `No ${MODE_NAMES[mode]} Result yet. Configure settings and click Generate.`;

const generate = async (page, options) => {
  const outcome = await generateAndWaitForResult(page, options);
  await settleLive(page);
  return outcome;
};

const showMode = async (page, mode) => {
  await page.getByRole('button', { name: MODE_NAMES[mode], exact: true }).click();
  await page.waitForFunction((value) => window.__GBDRAW_APP__?.mode === value, mode);
  await settleLive(page);
};

const setLinearSource = (page) => page.evaluate(async (text) => {
  window.__GBDRAW_APP__.setLinearSeqPrimaryFile(0, 'gb', new File([text], 'forced_label_underlay.gb', {
    type: 'text/plain', lastModified: 1000
  }));
  await window.Vue.nextTick();
}, readFileSync(FIXTURE, 'utf8'));

// Circular source loaded and a Linear source set, Circular shown.
const openBothSources = async (page) => {
  await openWithGenBank(page, FIXTURE, () => { window.__GBDRAW_APP__.autoLabelReflowEnabled = false; });
  await showMode(page, 'linear');
  await setLinearSource(page);
  await settleLive(page);
  await showMode(page, 'circular');
};

// The displayed Result: its committed identity, mode, and whether it is mounted.
const shown = (page) => page.evaluate(async () => {
  const ingestion = await import('/gbdraw/web/js/services/svg-result-ingestion.js');
  const { state } = await import('/gbdraw/web/js/state.js');
  const result = state.results.value[state.selectedResultIndex.value] || null;
  return {
    mode: state.mode.value,
    generatedMode: state.generatedMode.value,
    count: state.results.value.length,
    identity: result ? ingestion.getCommittedSvgResultRuntimeIdentity(result) : null,
    mounted: result ? ingestion.isCommittedSvgResultMounted(result) : false,
    content: result?.content || '',
    hasSvg: Boolean(state.svgContainer.value?.querySelector('svg'))
  };
});

const expectEmptyPreview = async (page, mode) => {
  const state = await shown(page);
  expect(state).toMatchObject({ mode, count: 0, identity: null, hasSvg: false });
  await expect(page.getByText(EMPTY_TEXT(mode), { exact: true })).toBeVisible();
};

// A Legend or plot title drag in layout mode, as the pointer does it.
const moveDecoration = async (page, role) => {
  await page.evaluate(async (targetRole) => {
    const app = window.__GBDRAW_APP__;
    app.layoutRepositionMode = true;
    await window.Vue.nextTick();
    const svg = app.svgContainer.querySelector('svg');
    const target = svg.querySelector(`[data-gbdraw-composition-role="${targetRole}"]`);
    if (!target) throw new Error(`No ${targetRole} to move`);
    const bounds = target.getBoundingClientRect();
    const x = bounds.left + Math.min(bounds.width / 2, 8);
    const y = bounds.top + Math.min(bounds.height / 2, 8);
    const mouse = (type, dx, dy, buttons) => new MouseEvent(type, {
      bubbles: true, cancelable: true, clientX: x + dx, clientY: y + dy, buttons, view: window
    });
    const frame = () => new Promise((resolve) => requestAnimationFrame(resolve));
    target.dispatchEvent(mouse('mousedown', 0, 0, 1));
    await frame();
    const moveTarget = targetRole === 'legend' ? svg : document;
    moveTarget.dispatchEvent(mouse('mousemove', 12, -9, 1));
    await frame();
    moveTarget.dispatchEvent(mouse('mouseup', 12, -9, 0));
  }, role);
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_HISTORY__?.undoLabel?.() || ''))
    .toBe(role === 'title' ? 'Move plot title' : 'Move legend');
  await page.evaluate(() => { window.__GBDRAW_APP__.layoutRepositionMode = false; });
  await settleLive(page);
};

const displayedMove = (page, role) => page.evaluate(async (targetRole) => {
  const { compositionUserDeltas } = await import('/gbdraw/web/js/app/legend-layout/composition-actions.js');
  const svg = window.__GBDRAW_APP__.svgContainer.querySelector('svg');
  return compositionUserDeltas(svg)[targetRole];
}, role);

const expectSameMove = (actual, expected) => {
  expect(actual).toHaveLength(2);
  actual.forEach((value, index) => expect(value).toBeCloseTo(expected[index], 5));
};

const history = async (page, step) => {
  await page.evaluate((name) => window.__GBDRAW_HISTORY__[name](), step);
  await settleLive(page);
};

const save = async (page, path) => {
  const download = page.waitForEvent('download', { timeout: 180_000 });
  await page.getByRole('banner').getByRole('button', { name: 'Save Session', exact: true }).click();
  await (await download).saveAs(path);
  return JSON.parse(gunzipSync(readFileSync(path)));
};

const dialogPages = new WeakSet();
const load = async (page, path) => {
  if (!dialogPages.has(page)) {
    page.on('dialog', (dialog) => dialog.accept());
    dialogPages.add(page);
  }
  await openApp(page);
  await page.locator('input[accept^=".json,"]').setInputFiles(path);
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending, null, { timeout: 180_000 });
  expect(await page.evaluate(() => window.__GBDRAW_APP__.errorLog)).toBeNull();
  await settleLive(page);
};

// The drawn features and legend rows of an SVG text, positions included.
const drawing = (page, content) => page.evaluate((text) => {
  const root = new DOMParser().parseFromString(text, 'image/svg+xml').documentElement;
  return {
    features: [...root.querySelectorAll('[data-gbdraw-feature-id]')].map((element) => [
      element.getAttribute('data-gbdraw-feature-id'), element.getAttribute('d'), element.getAttribute('fill')
    ]),
    legend: root.querySelector('[data-gbdraw-composition-role="legend"]')?.getAttribute('transform') || null
  };
}, content);

test('each mode keeps its own Result, with its Legend drag, across switches', async ({ page }) => {
  test.setTimeout(240_000);
  await openBothSources(page);
  await generate(page);
  await moveDecoration(page, 'legend');
  const circular = await shown(page);
  const circularMove = await displayedMove(page, 'legend');
  expect(circular).toMatchObject({ mode: 'circular', generatedMode: 'circular', count: 1, mounted: true });

  await showMode(page, 'linear');
  await expectEmptyPreview(page, 'linear');
  await generate(page);
  const linear = await shown(page);
  expect(linear).toMatchObject({ mode: 'linear', generatedMode: 'linear', count: 1, mounted: true });
  expect(linear.identity).not.toBe(circular.identity);
  expect(await displayedMove(page, 'legend')).toEqual([0, 0]);

  await showMode(page, 'circular');
  const returned = await shown(page);
  expect(returned).toMatchObject({ mode: 'circular', generatedMode: 'circular', count: 1, mounted: true });
  expect(returned.identity).toBe(circular.identity);
  expectSameMove(await displayedMove(page, 'legend'), circularMove);

  await showMode(page, 'linear');
  expect((await shown(page)).identity).toBe(linear.identity);
});

// PR-1 (PD-OI-086): each mode's drawing keeps its own Legend edits. A rename,
// an added row, and a moved row made on the Circular Result stay in Circular:
// Linear's first Generate draws its own rows, and the Circular Result keeps the
// edit when it is shown again, with Sort by default back to its own generated
// order. (Until PR-1 these edits reached the other mode: REVIEW-1 M1.)
const legendState = async (page) => ({
  rows: await page.evaluate(() => window.__GBDRAW_APP__.legendEntries.map((entry) => entry.caption)),
  // The drawn rows, sorted: a duplicate row shows twice.
  drawn: (await semanticSnapshot(page)).legend.map(({ caption }) => caption).sort()
});
const sorted = (captions) => [...captions].sort();
const CIRCULAR_ROWS = ['CDS', 'repeat_region', 'GC content', 'GC skew (+)', 'GC skew (-)'];
const LINEAR_ROWS = ['CDS', 'repeat_region'];
const LEGEND_EDITS = [
  {
    name: 'a rename',
    edit: (page) => evaluateWithRetainedPromise(page, async () => {
      const app = window.__GBDRAW_APP__;
      await app.renameLegendEntry(app.legendEntries.findIndex((entry) => entry.caption === 'CDS'), 'Coding');
    }),
    // The Circular rows Sort by default shows.
    byDefault: ['Coding', ...CIRCULAR_ROWS.slice(1)]
  },
  {
    // An editor-added row comes from a Session (R15-2).
    name: 'an added row',
    edit: (page) => loadEditorLegendRows(page, [['Extra', '#7b2cbf']]),
    byDefault: [...CIRCULAR_ROWS, 'Extra']
  },
  {
    name: 'a moved row',
    edit: (page) => page.evaluate(() => {
      const app = window.__GBDRAW_APP__;
      app.moveLegendEntryUp(app.legendEntries.findIndex((entry) => entry.caption === 'repeat_region'));
    }),
    byDefault: CIRCULAR_ROWS
  }
];
for (const { name, edit, byDefault } of LEGEND_EDITS) {
  test(`${name} on the Circular Legend stays in Circular and on return`, async ({ page }) => {
    test.setTimeout(240_000);
    await openBothSources(page);
    await generate(page);
    await edit(page);
    await expect.poll(async () => sorted((await legendState(page)).rows)).toEqual(sorted(byDefault));
    await settleLive(page);
    const circular = (await legendState(page)).rows;
    expect(circular, 'the edit').not.toEqual(CIRCULAR_ROWS);

    await showMode(page, 'linear');
    await generate(page);
    const generated = await legendState(page);
    expect(generated.rows, 'Linear Legend rows').toEqual(LINEAR_ROWS);
    expect(generated.drawn, 'Linear Result').toEqual(sorted(LINEAR_ROWS));

    await showMode(page, 'circular');
    const returned = await legendState(page);
    expect(returned.rows, 'Circular Legend rows after the return').toEqual(circular);
    expect(returned.drawn, 'Circular Result after the return').toEqual(sorted(circular));
    await generate(page);
    expect((await legendState(page)).rows, 'the next Circular Generate').toEqual(circular);
    // Sort by default shows the Circular Result's own generated order.
    await page.evaluate(() => window.__GBDRAW_APP__.sortLegendEntriesByDefault());
    await settleLive(page);
    expect((await legendState(page)).rows, 'Sort by default').toEqual(byDefault);
  });
}

// With both Results, a Result shown again by a switch shows what its next
// Generate draws, in its own order: a rename, an added row, or a delete made on
// the other mode's Result leaves it, and a switch alone reorders nothing.
const bothResults = async (page) => {
  await openBothSources(page);
  await generate(page);
  await showMode(page, 'linear');
  await generate(page);
};
// An editor-added row comes from a Session (R15-2).
const addExtraRow = async (page) => {
  await loadEditorLegendRows(page, [['Extra', '#7b2cbf']]);
  await settleLive(page);
};
const renameRow = async (page, from, to) => {
  await evaluateWithRetainedPromise(page, async ({ source, target }) => {
    const app = window.__GBDRAW_APP__;
    await app.renameLegendEntry(app.legendEntries.findIndex((entry) => entry.caption === source), target);
  }, { source: from, target: to });
  await settleLive(page);
};
// The rows on screen equal the drawn Legend and what the next Generate draws.
const expectShownAsGenerated = async (page, rows, label) => {
  const shownState = await legendState(page);
  expect(shownState.rows, `${label}: Legend rows`).toEqual(rows);
  expect(shownState.drawn, `${label}: drawn Legend`).toEqual(sorted(rows));
};
// Each Legend row's anchor, in reading order: the next Generate lays the rows
// out in the same slots, with no gap.
const legendGeometry = (page) => page.evaluate(() => window.__GBDRAW_APP__.legendEntries
  .map((entry) => [entry.caption, Math.round(Number(entry.xPos) || 0), Math.round(Number(entry.yPos) || 0)]));
const expectGeneratedGeometry = async (page, label) => {
  const shownGeometry = await legendGeometry(page);
  await generate(page);
  expect(await legendGeometry(page), `${label}: the next Generate's Legend slots`).toEqual(shownGeometry);
};
const DEFAULT_ROWS = { circular: CIRCULAR_ROWS, linear: LINEAR_ROWS };

for (const [first, second] of [['circular', 'linear'], ['linear', 'circular']]) {
  test(`a Legend rename made in ${second} with both Results leaves the ${first} Result`, async ({ page }) => {
    test.setTimeout(300_000);
    await bothResults(page);
    await showMode(page, second);
    await renameRow(page, 'CDS', 'Coding');
    const renamed = (await legendState(page)).rows;
    expect(sorted(renamed)).toEqual(sorted(['Coding', ...DEFAULT_ROWS[second].slice(1)]));
    await showMode(page, first);
    await expectShownAsGenerated(page, DEFAULT_ROWS[first], `${first} shown again`);
    await generate(page);
    await expectShownAsGenerated(page, DEFAULT_ROWS[first], `the next ${first} Generate`);
    await showMode(page, second);
    await expectShownAsGenerated(page, renamed, `${second} after a round trip`);
  });
}

test('an added row stays in its mode and moves neither Legend in a round trip', async ({ page }) => {
  test.setTimeout(300_000);
  await bothResults(page);
  await addExtraRow(page);
  const linear = [...LINEAR_ROWS, 'Extra'];
  await expectShownAsGenerated(page, linear, 'Linear after the addition');
  await showMode(page, 'circular');
  await expectShownAsGenerated(page, CIRCULAR_ROWS, 'Circular shown again');
  await showMode(page, 'linear');
  await expectShownAsGenerated(page, linear, 'Linear after a round trip');
  await showMode(page, 'circular');
  await expectShownAsGenerated(page, CIRCULAR_ROWS, 'Circular after a round trip');
  await expectGeneratedGeometry(page, 'Circular');
  await expectShownAsGenerated(page, CIRCULAR_ROWS, 'the next Circular Generate');
  await showMode(page, 'linear');
  await expectGeneratedGeometry(page, 'Linear with the added row');
  await expectShownAsGenerated(page, linear, 'the next Linear Generate');
});

test('an added row deleted in its mode stays deleted and never reaches the other mode', async ({ page }) => {
  test.setTimeout(300_000);
  await bothResults(page);
  await addExtraRow(page);
  await showMode(page, 'circular');
  await expectShownAsGenerated(page, CIRCULAR_ROWS, 'Circular without the added row');
  await showMode(page, 'linear');
  await expectShownAsGenerated(page, [...LINEAR_ROWS, 'Extra'], 'Linear with the added row');
  await evaluateWithRetainedPromise(page, async () => {
    const app = window.__GBDRAW_APP__;
    await app.deleteLegendEntry(app.legendEntries.findIndex((entry) => entry.caption === 'Extra'));
  });
  await settleLive(page);
  await expectShownAsGenerated(page, LINEAR_ROWS, 'Linear after the delete');
  await showMode(page, 'circular');
  await expectShownAsGenerated(page, CIRCULAR_ROWS, 'Circular shown again');
  await showMode(page, 'linear');
  await expectShownAsGenerated(page, LINEAR_ROWS, 'Linear after a round trip');
  await expectGeneratedGeometry(page, 'Linear after the delete');
  await expectShownAsGenerated(page, LINEAR_ROWS, 'the next Linear Generate');
  await showMode(page, 'circular');
  await expectGeneratedGeometry(page, 'Circular');
  await expectShownAsGenerated(page, CIRCULAR_ROWS, 'the next Circular Generate');
});

test('Undo and Redo walk a Generate in each mode across the switch between them', async ({ page }) => {
  test.setTimeout(240_000);
  await openBothSources(page);
  await generate(page);
  const g1 = await shown(page);

  await showMode(page, 'linear');
  await expectEmptyPreview(page, 'linear');
  await generate(page);
  const g2 = await shown(page);

  await history(page, 'undo');
  await expectEmptyPreview(page, 'linear');
  await history(page, 'undo');
  expect(await shown(page)).toMatchObject({ mode: 'circular', identity: g1.identity, mounted: true });
  await history(page, 'redo');
  await expectEmptyPreview(page, 'linear');
  await history(page, 'redo');
  expect(await shown(page)).toMatchObject({ mode: 'linear', identity: g2.identity, mounted: true });
  await showMode(page, 'circular');
  expect(await shown(page)).toMatchObject({ mode: 'circular', identity: g1.identity, mounted: true });
});

test('a failed Linear Generate keeps both Results', async ({ page }) => {
  test.setTimeout(240_000);
  await openBothSources(page);
  await generate(page);
  const circular = await shown(page);
  await showMode(page, 'linear');
  await generate(page);
  const linear = await shown(page);

  await page.evaluate(() => { window.__GBDRAW_APP__.adv.scale_font_size = '1e-50x'; });
  await generate(page, { expectedStatus: 'error', requireCommittedResult: false });
  expect(await shown(page)).toMatchObject({ mode: 'linear', identity: linear.identity, mounted: true });
  await page.evaluate(() => { window.__GBDRAW_APP__.adv.scale_font_size = null; });
  await showMode(page, 'circular');
  expect(await shown(page)).toMatchObject({ mode: 'circular', identity: circular.identity, mounted: true });
});

test('the mode switch waits for Generate (OIPC-C07)', async ({ page }) => {
  test.setTimeout(240_000);
  await openBothSources(page);
  const refused = await evaluateWithRetainedPromise(page, async () => {
    const app = window.__GBDRAW_APP__;
    const run = app.runAnalysis();
    await window.Vue.nextTick();
    const button = [...document.querySelectorAll('button')].find((element) => element.textContent.trim() === 'Linear');
    const outcome = { processing: app.processing, disabled: button.disabled, result: app.setDiagramMode('linear') };
    await run;
    return outcome;
  });
  expect(refused).toMatchObject({ processing: true, disabled: true, result: { status: 'busy' } });
  expect(await shown(page)).toMatchObject({ mode: 'circular', generatedMode: 'circular', count: 1 });
});

test('Save writes every Result and a fresh Load shows each mode\'s own Result', async ({ page, browser }, info) => {
  test.setTimeout(360_000);
  await openBothSources(page);
  await generate(page);
  await moveDecoration(page, 'legend');
  const circular = await shown(page);
  await showMode(page, 'linear');
  await generate(page);
  const linear = await shown(page);

  const saved = await save(page, info.outputPath('both.gbdraw-session.json.gz'));
  expect(saved.ui.mode).toBe('linear');
  expect(saved.renderRequest.mode).toBe('linear');
  expect(saved.otherModeResult.renderRequest.mode).toBe('circular');
  expect(saved.otherModeResult.results).toHaveLength(1);

  const fresh = await freshPage(browser);
  try {
    await load(fresh, info.outputPath('both.gbdraw-session.json.gz'));
    const loadedLinear = await shown(fresh);
    expect(loadedLinear).toMatchObject({ mode: 'linear', generatedMode: 'linear', count: 1, mounted: true });
    expect(await drawing(fresh, loadedLinear.content)).toEqual(await drawing(page, linear.content));
    await showMode(fresh, 'circular');
    const loadedCircular = await shown(fresh);
    expect(loadedCircular).toMatchObject({ mode: 'circular', generatedMode: 'circular', count: 1, mounted: true });
    expect(await drawing(fresh, loadedCircular.content)).toEqual(await drawing(page, circular.content));

    const resaved = await save(fresh, info.outputPath('both-resaved.gbdraw-session.json.gz'));
    expect(resaved.ui.mode).toBe('circular');
    expect(resaved.renderRequest.mode).toBe('circular');
    expect(resaved.otherModeResult.renderRequest.mode).toBe('linear');
    // Load admission may re-serialize a Result's markup; its drawing stays.
    expect(await drawing(fresh, resaved.results[0].content))
      .toEqual(await drawing(page, saved.otherModeResult.results[0].content));
    expect(await drawing(fresh, resaved.otherModeResult.results[0].content))
      .toEqual(await drawing(page, saved.results[0].content));
    await expectNoModeMismatch(fresh);
  } finally {
    await fresh.context().close();
  }
});

test('a Session with one Result opens on that Result\'s mode', async ({ page, browser }, info) => {
  test.setTimeout(240_000);
  await openBothSources(page);
  await generate(page);
  const circular = await shown(page);
  await showMode(page, 'linear');
  await expectEmptyPreview(page, 'linear');

  const saved = await save(page, info.outputPath('one.gbdraw-session.json.gz'));
  expect(saved.ui.mode).toBe('linear');
  expect(saved.renderRequest.mode).toBe('circular');
  expect(Object.hasOwn(saved, 'otherModeResult')).toBe(false);

  const fresh = await freshPage(browser);
  try {
    await load(fresh, info.outputPath('one.gbdraw-session.json.gz'));
    const loaded = await shown(fresh);
    expect(loaded).toMatchObject({ mode: 'circular', generatedMode: 'circular', count: 1, mounted: true });
    expect(await drawing(fresh, loaded.content)).toEqual(await drawing(page, circular.content));
    await showMode(fresh, 'linear');
    await expectEmptyPreview(fresh, 'linear');
    await expectNoModeMismatch(fresh);
  } finally {
    await fresh.context().close();
  }
});

// §7.7 risk 2: the LOSAT caches and protein manifest stay one shared union, so
// a Linear BLASTP Result saved while Circular is shown keeps its evidence.
test('a Linear BLASTP Result saved while Circular is shown keeps its comparison cache', async ({ page, browser }, info) => {
  test.setTimeout(360_000);
  const BGC = 'gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json';
  const cache = (target) => target.evaluate(async () => {
    const { state } = await import('/gbdraw/web/js/state.js');
    return { entries: state.losatCache.value.size, manifest: Object.keys(state.proteinIdentityManifest.value?.proteinSets || {}).length };
  });
  await load(page, BGC);
  const linear = await shown(page);
  expect(linear).toMatchObject({ mode: 'linear', count: 1, mounted: true });
  const loadedCache = await cache(page);
  expect(loadedCache.entries).toBeGreaterThan(0);

  await showMode(page, 'circular');
  await page.getByLabel('GenBank/DDBJ File', { exact: true }).setInputFiles(FIXTURE);
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.circularRecordList.length)).toBeGreaterThan(0);
  await settleLive(page);
  await generate(page);
  expect(await cache(page)).toEqual(loadedCache);

  const saved = await save(page, info.outputPath('blastp.gbdraw-session.json.gz'));
  expect(saved.renderRequest.mode).toBe('circular');
  expect(saved.otherModeResult.renderRequest.mode).toBe('linear');
  expect(saved.losatCache.entries).toHaveLength(loadedCache.entries);

  const fresh = await freshPage(browser);
  try {
    await load(fresh, info.outputPath('blastp.gbdraw-session.json.gz'));
    expect(await shown(fresh)).toMatchObject({ mode: 'circular', count: 1, mounted: true });
    expect(await cache(fresh)).toEqual(loadedCache);
    await showMode(fresh, 'linear');
    const restored = await shown(fresh);
    expect(restored).toMatchObject({ mode: 'linear', generatedMode: 'linear', count: 1, mounted: true });
    expect(await drawing(fresh, restored.content)).toEqual(await drawing(page, linear.content));
    await expectNoModeMismatch(fresh);
  } finally {
    await fresh.context().close();
  }
});

// REVIEW-1 M2 (PD-OI-079): Save keeps every Result, so a Result an older gbdraw
// saved, kept while the other mode is shown, still asks for one Generate. The
// error names its mode, and the error's Generate action shows that mode and
// runs Generate there; Save then writes both Results.
test('Save names the hidden mode whose older Result needs Generate, and Generate runs there', async ({ page }, info) => {
  test.setTimeout(480_000);
  await load(page, 'tests/fixtures/sessions/BGC0000708-BGC0000713.v39.gbdraw-session.json.gz');
  expect(await shown(page)).toMatchObject({ mode: 'linear', count: 1, mounted: true });
  await showMode(page, 'circular');
  await page.getByLabel('GenBank/DDBJ File', { exact: true }).setInputFiles(FIXTURE);
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.circularRecordList.length)).toBeGreaterThan(0);
  await settleLive(page);
  await generate(page);
  const circular = await shown(page);

  const downloads = [];
  const onDownload = (download) => downloads.push(download.suggestedFilename());
  page.on('download', onDownload);
  await page.getByRole('banner').getByRole('button', { name: 'Save Session', exact: true }).click();
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionSavePending, null, { timeout: 120_000 });
  await settleLive(page);
  page.off('download', onDownload);
  expect(downloads).toEqual([]);
  const error = await page.evaluate(() => {
    const { code, summary } = window.__GBDRAW_APP__.errorLog || {};
    return { code, summary };
  });
  expect(error.code).toBe('SESSION_SAVE_REQUIRES_GENERATE');
  expect(error.summary).toContain('Diagram: Linear.');

  await page.locator('.border-l-red-500').getByRole('button', { name: 'Generate', exact: true }).click();
  await page.waitForFunction(() => window.__GBDRAW_APP__.mode === 'linear', null, { timeout: 60_000 });
  await page.waitForFunction(() => !window.__GBDRAW_APP__.processing, null, { timeout: 300_000 });
  await settleLive(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.errorLog)).toBeNull();
  expect(await shown(page)).toMatchObject({ mode: 'linear', generatedMode: 'linear', count: 1, mounted: true });
  const saved = await save(page, info.outputPath('older-linear.gbdraw-session.json.gz'));
  expect(saved.renderRequest.mode).toBe('linear');
  expect(saved.otherModeResult.renderRequest.mode).toBe('circular');
  expect(await drawing(page, saved.otherModeResult.results[0].content)).toEqual(await drawing(page, circular.content));
});

// REVIEW-1 M3: Load admits both Result sets through one function. The Linear
// draft's similarity alignment plan and record translations come from the
// Linear set's committed request wherever the set sits, so a Linear Generate
// after a Load that opened on Circular keeps the alignment.
test('an aligned Linear Result saved while Circular is shown keeps its alignment after Load', async ({ page, browser }, info) => {
  test.setTimeout(480_000);
  const alignment = (target) => target.evaluate(async () => {
    const { state } = await import('/gbdraw/web/js/state.js');
    const { getCommittedCanonicalRenderRequest } = await import('/gbdraw/web/js/services/config.js');
    const committed = getCommittedCanonicalRenderRequest();
    return {
      plan: JSON.stringify(state.similarityAlignmentPlan.value),
      translations: JSON.stringify(state.linearRecordTranslations.value),
      committed: committed?.mode === 'linear' ? JSON.stringify([
        committed.layout?.similarityAlignment ?? null, committed.layout?.recordTranslations ?? null
      ]) : null
    };
  });
  await load(page, 'gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json');
  const aligned = await alignment(page);
  expect(aligned.plan).not.toBe('null');
  expect(aligned.committed).not.toBeNull();
  await showMode(page, 'circular');
  await page.getByLabel('GenBank/DDBJ File', { exact: true }).setInputFiles(FIXTURE);
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.circularRecordList.length)).toBeGreaterThan(0);
  await settleLive(page);
  await generate(page);
  const saved = await save(page, info.outputPath('aligned.gbdraw-session.json.gz'));
  expect(saved.renderRequest.mode).toBe('circular');
  expect(saved.otherModeResult.renderRequest.layout.similarityAlignment).not.toBeNull();

  const fresh = await freshPage(browser);
  try {
    await load(fresh, info.outputPath('aligned.gbdraw-session.json.gz'));
    expect(await shown(fresh)).toMatchObject({ mode: 'circular', count: 1, mounted: true });
    expect(await alignment(fresh)).toMatchObject({ plan: aligned.plan, translations: aligned.translations });
    await showMode(fresh, 'linear');
    expect(await alignment(fresh)).toEqual(aligned);
    await generate(fresh);
    expect(await alignment(fresh), 'Linear Generate keeps the alignment').toEqual(aligned);
    await expectNoModeMismatch(fresh);
  } finally {
    await fresh.context().close();
  }
});

// A two-Result Session with the Linear set at the top level, the Circular set
// in `otherModeResult`, and `ui.mode` Circular. The fixture's drawings are laid
// out Linear first by the Python layout owner (`write_session_drawings`), then
// `gbdraw linear --session ... --session_output` re-saves it in place (D-15-36).
const linearTopResavedSession = (info) => {
  const directory = info.outputPath('cli-resave');
  mkdirSync(directory, { recursive: true });
  execFileSync('python', ['-c', [
    'import sys',
    'from pathlib import Path',
    'from gbdraw.cli_utils import session as cli_session',
    'from gbdraw.session import load_session_document, session_drawing_artifacts',
    'from gbdraw.session_drawings import write_session_drawings',
    'from gbdraw.session_io import load_session, write_session_json',
    'from tests.utils.two_mode_session import two_mode_session',
    'out = Path(sys.argv[1])',
    'two = load_session_document(two_mode_session())',
    "views = {d.mode: dict(session_drawing_artifacts(two, d.id).fields) for d in two.drawings}",
    "write_session_json(out / 'two.gbdraw-session.json', write_session_drawings([views['linear'], views['circular']], active=1))",
    "assert cli_session.render_canonical_session_if_present(load_session(str(out / 'two.gbdraw-session.json')),"
      + " mode='linear', output_override=str(out / 'replayed'), format_override='svg', save_session=False,"
      + " session_output=str(out / 'resaved.gbdraw-session.json'))"
  ].join('\n'), directory], { cwd: process.cwd(), stdio: 'pipe' });
  const path = `${directory}/resaved.gbdraw-session.json`;
  const bytes = readFileSync(path);
  const session = JSON.parse(bytes[0] === 0x1f && bytes[1] === 0x8b ? gunzipSync(bytes) : bytes);
  expect(session.ui.mode).toBe('circular');
  expect(session.renderRequest.mode).toBe('linear');
  expect(session.otherModeResult.renderRequest.mode).toBe('circular');
  return { path, session };
};

// REVIEW-1 P2: the Load of such a Session shows the Circular set, and the
// top-level Linear set waits in its mode: it is never displayed, bound, or
// projected while Circular is shown.
test('a CLI re-saved Session opens on its shown mode with each Result in its own mode', async ({ page }, info) => {
  test.setTimeout(480_000);
  const { path, session } = linearTopResavedSession(info);
  const holds = async (label) => {
    const state = await shown(page);
    expect(state.generatedMode === state.mode || state.count === 0, label).toBe(true);
    return state;
  };
  await load(page, path);
  expect(await page.evaluate(() => window.__E1_READINESS__), 'the Load shows one Result, Circular').toEqual(['session-load']);
  const circular = await holds('Load');
  expect(circular).toMatchObject({ mode: 'circular', count: 1, mounted: true });
  expect(await drawing(page, circular.content)).toEqual(await drawing(page, session.otherModeResult.results[0].content));
  // The Circular Result's own Legend rows, none of the Linear set's.
  const ownRows = await page.evaluate((content) => [...new Set([
    ...new DOMParser().parseFromString(content, 'image/svg+xml').querySelectorAll('g[data-legend-key]')
  ].map((entry) => entry.getAttribute('data-legend-key')))], session.otherModeResult.results[0].content);
  expect((await page.evaluate(() => window.__GBDRAW_APP__.legendEntries.map((entry) => entry.caption))).sort())
    .toEqual(ownRows.sort());
  await showMode(page, 'linear');
  const linear = await holds('switch to Linear');
  expect(linear).toMatchObject({ mode: 'linear', count: 1, mounted: true });
  expect(await drawing(page, linear.content)).toEqual(await drawing(page, session.results[0].content));
  await history(page, 'undo');
  expect((await holds('Undo of the switch')).identity).toBe(circular.identity);
  await history(page, 'redo');
  expect((await holds('Redo of the switch')).identity).toBe(linear.identity);
  await showMode(page, 'circular');
  await holds('switch back');
  const resaved = await save(page, info.outputPath('web-resaved.gbdraw-session.json.gz'));
  expect(resaved.ui.mode).toBe('circular');
  expect(resaved.renderRequest.mode).toBe('circular');
  expect(resaved.otherModeResult.renderRequest.mode).toBe('linear');
});

// A Session opened on the mode of its other Result set (a Web Session that
// `gbdraw linear --session_output` re-saved) shows that Result's own Legend
// rows, and its next Generate draws the same rows.
test('a Session opened on its other Result set generates the Legend it shows', async ({ page, browser }, info) => {
  test.setTimeout(480_000);
  await bothResults(page);
  await showMode(page, 'circular');
  const webPath = info.outputPath('web.gbdraw-session.json.gz');
  await save(page, webPath);
  const resavedPath = info.outputPath('cli-resaved.gbdraw-session.json');
  execFileSync('python', ['-c', [
    'import sys',
    'from gbdraw.cli_utils import session as cli_session',
    'from gbdraw.session_io import load_session',
    "assert cli_session.render_canonical_session_if_present(load_session(sys.argv[1]), mode='linear',"
      + " output_override=sys.argv[3], format_override='svg', save_session=False, session_output=sys.argv[2])"
  ].join('\n'), webPath, resavedPath, info.outputPath('replayed')], { cwd: process.cwd(), stdio: 'pipe' });
  const fresh = await freshPage(browser);
  try {
    await load(fresh, resavedPath);
    expect(await shown(fresh)).toMatchObject({ mode: 'circular', generatedMode: 'circular', count: 1, mounted: true });
    const rows = (await legendState(fresh)).rows;
    expect(sorted(rows)).toEqual(sorted(CIRCULAR_ROWS));
    await expectShownAsGenerated(fresh, rows, 'the opened Result');
    await generate(fresh);
    await expectShownAsGenerated(fresh, rows, 'its next Generate');
    await expectNoModeMismatch(fresh);
  } finally {
    await fresh.context().close();
  }
});

// REVIEW-1 m9: a Load canceled within its History baseline, after it showed
// the Result set of the mode it opens on, restores the displayed artifact (its
// committed Session and match-sequence owner) and the other mode's Result.
test('a Load canceled within its History baseline keeps both modes\' Results', async ({ page }, info) => {
  test.setTimeout(480_000);
  const { path } = linearTopResavedSession(info);
  await openBothSources(page);
  await generate(page);
  await showMode(page, 'linear');
  await generate(page);
  const linear = await shown(page);
  await showMode(page, 'circular');
  const circular = await shown(page);
  const counts = () => page.evaluate(() => [window.__GBDRAW_HISTORY__.getUndoCount(), window.__GBDRAW_HISTORY__.getRedoCount()]);
  const historyBefore = await counts();
  await page.evaluate(async () => {
    const config = await import('/gbdraw/web/js/services/config.js');
    const { state } = await import('/gbdraw/web/js/state.js');
    window.__E1_BEFORE_LOAD__ = {
      committed: config.getCommittedCanonicalSession(),
      owner: state.matchSequenceRegistry.captureTrustedOwner()
    };
    window.__GBDRAW_TEST_HOOKS__.onSessionLifecycleEvent = (event) => {
      if (event.name !== 'history-baseline-start') return;
      window.__E1_CANCELED_IN_BASELINE__ = true;
      config.disposeSessionOperations();
    };
  });
  await page.locator('input[accept^=".json,"]').setInputFiles(path);
  await page.waitForFunction(() => window.__E1_CANCELED_IN_BASELINE__ && !window.__GBDRAW_APP__.sessionImportPending,
    null, { timeout: 180_000 });
  await settleLive(page);
  expect(await page.evaluate(async () => {
    const config = await import('/gbdraw/web/js/services/config.js');
    const { state } = await import('/gbdraw/web/js/state.js');
    return {
      committed: config.getCommittedCanonicalSession() === window.__E1_BEFORE_LOAD__.committed,
      owner: state.matchSequenceRegistry.captureTrustedOwner() === window.__E1_BEFORE_LOAD__.owner
    };
  })).toEqual({ committed: true, owner: true });
  expect(await shown(page)).toMatchObject({ mode: 'circular', generatedMode: 'circular', identity: circular.identity, mounted: true });
  expect(await counts()).toEqual(historyBefore);
  await showMode(page, 'linear');
  expect(await shown(page)).toMatchObject({ mode: 'linear', generatedMode: 'linear', identity: linear.identity, mounted: true });
});

// REVIEW-1 m8: the switch also waits for an Undo or Redo, so a restore always
// finds the mode it was captured in. The mode buttons read the same state.
test('the mode switch waits for Undo and Redo', async ({ page }) => {
  test.setTimeout(240_000);
  await openBothSources(page);
  await generate(page);
  const refused = await evaluateWithRetainedPromise(page, async () => {
    const app = window.__GBDRAW_APP__;
    const history = window.__GBDRAW_HISTORY__;
    const run = history.undo();
    await window.Vue.nextTick();
    const button = [...document.querySelectorAll('button')].find((element) => element.textContent.trim() === 'Linear');
    const outcome = {
      pending: history.restoring.value || Boolean(history.traversalPending?.()),
      disabled: button.disabled,
      result: app.setDiagramMode('linear')
    };
    await run;
    return outcome;
  });
  expect(refused).toMatchObject({ disabled: true, result: { status: 'busy' }, pending: true });
  await settleLive(page);
  // The Undo of the Generate finished in the mode it was made in.
  await expectEmptyPreview(page, 'circular');
});

// REVIEW-2 P-a: an Undo or Redo of a mode switch waits for a label rerender,
// as the mode buttons do.
test('an Undo of a mode switch waits for a label rerender', async ({ page }) => {
  test.setTimeout(300_000);
  await bothResults(page);
  await showMode(page, 'circular');
  await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    await app.setFeatureColorValue(app.extractedFeatures.find((feature) => feature.type === 'CDS'), '#123456');
  });
  await settleLive(page);
  const attempt = await evaluateWithRetainedPromise(page, async () => {
    const { state } = await import('/gbdraw/web/js/state.js');
    const owner = window.__GBDRAW_HISTORY__;
    // The Undo of the color restores the rule table and asks for the rerender.
    await owner.undo();
    for (let wait = 0; wait < 200 && !state.labelReflowProcessing.value; wait += 1) {
      await new Promise((resolve) => setTimeout(resolve, 10));
    }
    const rerendering = state.labelReflowProcessing.value;
    const undo = await owner.undo();
    return { rerendering, undo, mode: state.mode.value };
  });
  expect(attempt).toMatchObject({ rerendering: true, undo: { status: 'busy' }, mode: 'circular' });
  await settleLive(page);
  await history(page, 'undo');
  expect(await shown(page)).toMatchObject({ mode: 'linear', generatedMode: 'linear', count: 1, mounted: true });
});

// REVIEW-2 R: a Load that fails inside its commit rolls back and leaves the
// same Result under the same root; an edit after it commits into that Result.
test('an edit after a rolled-back Load commits into the restored Result', async ({ page }) => {
  test.setTimeout(300_000);
  const LAMBDA = 'gbdraw/web/gallery/sessions/lambda_basic_linear.gbdraw-session.json';
  await load(page, LAMBDA);
  await page.evaluate(async () => {
    window.__E1_ROOT__ = window.__GBDRAW_APP__.svgContainer.querySelector('svg');
    // One failure inside the synchronous commit of the next Load.
    const { state } = await import('/gbdraw/web/js/state.js');
    const target = state.linearRecordTranslations;
    const descriptor = Object.getOwnPropertyDescriptor(Object.getPrototypeOf(target), 'value');
    window.__E1_INJECTED__ = 0;
    Object.defineProperty(target, 'value', {
      configurable: true,
      get() { return descriptor.get.call(this); },
      set(next) {
        if (window.__E1_INJECTED__ === 0 && window.__GBDRAW_APP__.sessionImportPending) {
          window.__E1_INJECTED__ = 1;
          delete target.value;
          throw new Error('injected Load failure');
        }
        descriptor.set.call(this, next);
      }
    });
  });
  await page.locator('input[accept^=".json,"]').setInputFiles(LAMBDA);
  await page.waitForFunction(() => window.__E1_INJECTED__ === 1 && !window.__GBDRAW_APP__.sessionImportPending,
    null, { timeout: 180_000 });
  await settleLive(page);
  const edit = await page.evaluate(async () => {
    const { state } = await import('/gbdraw/web/js/state.js');
    const app = window.__GBDRAW_APP__;
    const root = app.svgContainer.querySelector('svg');
    const outcome = await app.updateLegendEntryColor(0, '#ff00ff');
    await window.Vue.nextTick();
    const result = state.results.value[state.selectedResultIndex.value];
    return {
      sameRoot: root === window.__E1_ROOT__,
      outcome,
      mounted: root.outerHTML.toLowerCase().includes('#ff00ff'),
      committed: String(result?.content || '').toLowerCase().includes('#ff00ff')
    };
  });
  expect(edit).toEqual({ sameRoot: true, outcome: true, mounted: true, committed: true });
});

// UJ-02 (gbdraw-41 journey): a mode without a Result of its own shows no
// Result of the other mode. Its Preview says so, the record display reports no
// change pending Generate, and SVG export (the Result Preview card, shown only
// with a Result) is absent, so no other mode's Result can be exported. Shown
// again, a mode's own Result has its features, opens a feature popup, and
// exports itself.
const UJ02 = [
  { session: 'gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json', mode: 'circular', other: 'linear', features: 37 },
  { session: 'gbdraw/web/gallery/sessions/lambda_basic_linear.gbdraw-session.json', mode: 'linear', other: 'circular', features: 73 }
];
for (const { session, mode, other, features } of UJ02) {
  test(`UJ-02: a ${mode} Result is neither shown nor exported in ${other} and comes back whole`, async ({ page }) => {
    test.setTimeout(300_000);
    const svgButton = () => page.getByRole('button', { name: 'SVG', exact: true });
    const exportSvg = async () => {
      const [download] = await Promise.all([page.waitForEvent('download'), svgButton().click()]);
      const path = await download.path();
      return { name: download.suggestedFilename(), content: readFileSync(path, 'utf8') };
    };
    const featureIds = (content) => [...new Set(content.match(/data-gbdraw-feature-id="[^"]+"/g) || [])].sort();
    const pending = () => page.evaluate(() => {
      const value = window.__GBDRAW_APP__.recordDisplayControls.hasPendingChanges;
      return Boolean(value && typeof value === 'object' ? value.value : value);
    });
    await load(page, session);
    await showMode(page, mode);
    await generate(page);
    const own = await shown(page);
    expect(own).toMatchObject({ mode, generatedMode: mode, count: 1, mounted: true });
    expect(await page.evaluate(() => window.__GBDRAW_APP__.extractedFeatures.length)).toBe(features);
    const exported = await exportSvg();

    await showMode(page, other);
    // (a) the empty Preview, and no Result of the other mode mounted.
    await expectEmptyPreview(page, other);
    // (b) no record display change pending Generate.
    expect(await pending(), 'hasPendingChanges').toBe(false);
    await expect(page.getByText(/has changes pending Generate/)).toHaveCount(0);
    // (c) no export of the other mode's Result.
    await expect(svgButton()).toHaveCount(0);

    await showMode(page, mode);
    // (d) its own Result, with its features.
    expect(await shown(page)).toMatchObject({ mode, generatedMode: mode, identity: own.identity, mounted: true });
    expect(await page.evaluate(() => window.__GBDRAW_APP__.extractedFeatures.length)).toBe(features);
    expect(await pending(), 'hasPendingChanges').toBe(false);
    // (e) a feature popup opens.
    expect((await popup(page)).featureId).toBeTruthy();
    await closeEditor(page);
    // (f) SVG export downloads this Result.
    const again = await exportSvg();
    expect(again.name).toBe(exported.name);
    expect(featureIds(again.content)).toEqual(featureIds(exported.content));
    expect(featureIds(again.content).length).toBeGreaterThan(0);
  });
}

// The displayed Result always belongs to the shown mode (Phase E relies on it):
// after every path that installs a Result or a mode, `generatedMode === mode`
// or no Result is shown. Paths: a switch to a mode with and without a Result,
// a failed Generate, Reset Settings, Undo and Redo through intent (switch),
// artifact-handle (Generate) and checkpoint (Reset Settings) steps, a fresh
// Load of a two-Result Session and of a Session whose shown mode has no Result,
// and the rollback of a failed Load. The fallback UI restore of the History
// service is not reachable from the app (it always passes `applyUiStateData`).
test('the displayed Result belongs to the shown mode after every mode path', async ({ page, browser }, info) => {
  test.setTimeout(480_000);
  const holds = (target) => target.evaluate(async () => {
    const { state } = await import('/gbdraw/web/js/state.js');
    return {
      mode: state.mode.value,
      holds: state.generatedMode.value === state.mode.value || state.results.value.length === 0
    };
  });
  const check = async (target, label) => expect((await holds(target)).holds, label).toBe(true);

  await openBothSources(page);
  await generate(page);
  await check(page, 'Circular Generate');
  await showMode(page, 'linear');
  await check(page, 'switch to a mode without a Result');
  const oneResultFile = info.outputPath('linear-draft.gbdraw-session.json.gz');
  await save(page, oneResultFile);
  await generate(page);
  await check(page, 'Linear Generate');
  await showMode(page, 'circular');
  await check(page, 'switch to a mode with a Result');

  await page.evaluate(() => { window.__GBDRAW_APP__.adv.def_font_size = '1e-50x'; });
  await generate(page, { expectedStatus: 'error', requireCommittedResult: false });
  await check(page, 'failed Generate');
  await page.evaluate(() => { window.__GBDRAW_APP__.adv.def_font_size = null; });
  await settleLive(page);

  await page.evaluate(() => window.__GBDRAW_APP__.resetSettings());
  await settleLive(page);
  await check(page, 'Reset Settings');

  const counts = () => page.evaluate(() => [window.__GBDRAW_HISTORY__.getUndoCount(), window.__GBDRAW_HISTORY__.getRedoCount()]);
  const [undoCount] = await counts();
  for (let step = 0; step < undoCount; step += 1) {
    await history(page, 'undo');
    await check(page, `Undo ${step + 1}`);
  }
  expect((await counts())[0]).toBe(0);
  for (let step = 0; step < undoCount; step += 1) {
    await history(page, 'redo');
    await check(page, `Redo ${step + 1}`);
  }
  const bothFile = info.outputPath('both.gbdraw-session.json.gz');
  await save(page, bothFile);

  const fresh = await freshPage(browser);
  try {
    await load(fresh, bothFile);
    await check(fresh, 'Load of a two-Result Session');
    await showMode(fresh, 'linear');
    await check(fresh, 'switch after Load');
    await load(fresh, oneResultFile);
    expect((await holds(fresh)).mode, 'Load opens on the mode that has a Result').toBe('circular');
    await check(fresh, 'Load of a Session whose shown mode has no Result');
    await load(fresh, bothFile);
    const before = await shown(fresh);
    const failure = await evaluateWithRetainedPromise(fresh, async () => {
      const name = 'tobacco-chloroplast.gbdraw-session.json';
      const file = new File([await (await fetch(`/gbdraw/web/gallery/sessions/${name}`)).arrayBuffer()], name);
      const { importSession } = await import('/gbdraw/web/js/services/config.js');
      const result = await importSession({ target: { files: [file], value: 'selected' } }, {
        beforePreviewMount: () => { throw new Error('forced rollback'); }
      });
      return result.status;
    });
    expect(failure).toBe('error');
    await settleLive(fresh);
    await check(fresh, 'rollback of a failed Load');
    expect(await shown(fresh)).toMatchObject({ mode: before.mode, identity: before.identity });
    await expectNoModeMismatch(fresh);
  } finally {
    await fresh.context().close();
  }
});
