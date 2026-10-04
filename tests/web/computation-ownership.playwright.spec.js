// CW-01..04 (docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md#computation-ownership):
// comparison switching, ordinary input, Editor placement, help-tips, and
// Generate acceptance build no status-only tables, and each Generate builds
// its label table once.
const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { join } = require('node:path');
const { openApp, getDiagramWorkerActivity } = require('./helpers/app-lifecycle.cjs');

const source = readFileSync(join(__dirname, '../test_inputs/HmmtDNA.gbk'), 'utf8');
const SELECTOR_FUNCTIONS = ['buildFeatureMetadataMap', 'buildFeatureSelectorUniquenessIndexFromMetadata'];

const install = (page) => page.addInitScript(() => {
  const cw = window.__CW__ = { metrics: [], events: [], losatConstructions: 0 };
  window.__GBDRAW_TEST_HOOKS__ = {
    onStructuralMetric: (metric) => cw.metrics.push(metric),
    onSessionLifecycleEvent: (event) => cw.events.push(event.name)
  };
  const NativeWorker = window.Worker;
  window.Worker = new Proxy(NativeWorker, {
    construct(target, args) {
      if (/losat/i.test(String(args[0] || ''))) cw.losatConstructions += 1;
      return Reflect.construct(target, args, target);
    }
  });
});

// V8 precise coverage counts real calls; takePreciseCoverage resets counters.
const selectorCalls = async (cdp) => {
  const { result } = await cdp.send('Profiler.takePreciseCoverage');
  const counts = Object.fromEntries(SELECTOR_FUNCTIONS.map((name) => [name, 0]));
  for (const script of result) {
    if (!script.url.includes('/js/app/feature-selector.js')) continue;
    for (const fn of script.functions) {
      if (fn.functionName in counts) counts[fn.functionName] += fn.ranges[0]?.count || 0;
    }
  }
  return counts;
};

const settle = (page) => page.waitForFunction(() => {
  const app = window.__GBDRAW_APP__;
  return !app.processing && !app.labelReflowProcessing && !window.__GBDRAW_HISTORY__?.capturing?.value;
}, null, { timeout: 180_000 }).then(() => page.evaluate(() => new Promise((resolve) => {
  requestAnimationFrame(() => requestAnimationFrame(() => setTimeout(resolve, 0)));
})));

// One operation window: metric deltas, selector calls, Worker dispatch, History.
const observe = async (page, cdp, action) => {
  await settle(page);
  await selectorCalls(cdp);
  const start = await page.evaluate(() => ({
    metrics: window.__CW__.metrics.length, events: window.__CW__.events.length,
    losat: window.__CW__.losatConstructions, undo: window.__GBDRAW_HISTORY__.getUndoCount()
  }));
  const workersBefore = await getDiagramWorkerActivity(page);
  await action();
  await settle(page);
  const selector = await selectorCalls(cdp);
  const workersAfter = await getDiagramWorkerActivity(page);
  const end = await page.evaluate((start) => ({
    metrics: window.__CW__.metrics.slice(start.metrics),
    events: window.__CW__.events.slice(start.events),
    losat: window.__CW__.losatConstructions - start.losat,
    undo: window.__GBDRAW_HISTORY__.getUndoCount() - start.undo
  }), start);
  const counts = {};
  for (const { name, value } of end.metrics) counts[name] = (counts[name] || 0) + value;
  return {
    counts, events: end.events, selector, losat: end.losat, undo: end.undo,
    helpers: workersAfter.helpers - workersBefore.helpers,
    runs: workersAfter.runs - workersBefore.runs
  };
};

// A missing metric is a failure, never an implicit zero.
const count = (op, name) => {
  expect(Object.hasOwn(op.counts, name), `probe missing: ${name}`).toBe(true);
  return op.counts[name];
};
const expectNoStatusWork = (op, label) => {
  for (const name of ['canonicalRequestProjectionCount', 'generatedTableBuildCount', 'labelOverrideTableBuildCount']) {
    expect(op.counts[name] || 0, `${label}: ${name}`).toBe(0);
  }
  expect(op.selector, `${label}: selector metadata/index`).toEqual({
    buildFeatureMetadataMap: 0, buildFeatureSelectorUniquenessIndexFromMetadata: 0
  });
  expect({ helpers: op.helpers, runs: op.runs, losat: op.losat }, label)
    .toEqual({ helpers: 0, runs: 0, losat: 0 });
};

const clickGenerate = async (page) => {
  await page.getByRole('button', { name: 'Generate Diagram', exact: true }).click();
  await expect.poll(() => page.evaluate(() => window.__CW__.events.includes('generate.completed')),
    { timeout: 180_000 }).toBe(true);
  await page.evaluate(() => { window.__CW__.events = window.__CW__.events.filter((name) => name !== 'generate.completed'); });
};

const applyVisibilityOverride = (page) => page.evaluate(async () => {
  const app = window.__GBDRAW_APP__;
  const feature = app.extractedFeatures.find((candidate) => candidate?.svg_id);
  app.openFeatureEditorFromList(feature, null);
  app.clickedFeature.labelVisibility = 'on';
  await app.updateClickedFeatureLabelText();
  app.closeRightDrawer();
  return Object.keys(app.labelVisibilityOverrides).length;
});

test('status, selection, and input build no tables; each Generate builds its label table once', async ({ page }) => {
  test.setTimeout(300_000);
  await install(page);
  const cdp = await page.context().newCDPSession(page);
  await cdp.send('Profiler.enable');
  await cdp.send('Profiler.startPreciseCoverage', { callCount: true, detailed: false });
  await openApp(page);
  await page.evaluate(async (content) => {
    const app = window.__GBDRAW_APP__;
    app.mode = 'linear';
    await window.Vue.nextTick();
    while (app.linearSeqs.length < 2) app.addLinearSeq();
    for (let index = 0; index < 2; index += 1) {
      app.setLinearSeqPrimaryFile(index, 'gb', new File([content], `record-${index}.gbk`, { type: 'text/plain' }));
    }
    await window.Vue.nextTick();
  }, source);

  // Positive controls: Generate builds the request and tables, so the probes are live.
  const empty = await observe(page, cdp, () => clickGenerate(page));
  expect(count(empty, 'canonicalRequestProjectionCount')).toBeGreaterThan(0);
  expect(count(empty, 'generatedTableBuildCount')).toBeGreaterThan(0);
  // CW-03: empty overrides build no override metadata or index.
  expect(empty.counts.labelOverrideTableBuildCount || 0).toBe(0);
  expect(empty.runs).toBe(1);

  expect(await applyVisibilityOverride(page)).toBeGreaterThan(0);
  // CW-02: one label table per Generate, shared by the staged copy and the request.
  const override = await observe(page, cdp, () => clickGenerate(page));
  expect(count(override, 'labelOverrideTableBuildCount')).toBe(1);
  expect(override.selector.buildFeatureMetadataMap).toBeGreaterThan(0);
  expect(await page.evaluate(async () => Boolean((await import('./js/services/config.js'))
    .getCommittedCanonicalRenderRequest()?.diagramOptions?.labelOverrideFile))).toBe(true);
  // CW-04: every Generate renders and admits its catalog, even for unchanged input.
  for (const op of [override, await observe(page, cdp, () => clickGenerate(page))]) {
    expect(op.runs).toBe(1);
    expect(count(op, 'featureCatalogAdmissionCount')).toBe(1);
    expect(count(op, 'svgSanitizationCount')).toBeGreaterThan(0);
  }

  // CW-01: operation acceptance and selection updates start no heavy work.
  const compareButton = (target) => page.getByRole('button', {
    name: target === 'losat' ? 'Run LOSAT for all adjacent pairs' : 'Set no comparison', exact: true
  });
  for (const target of ['losat', 'none']) {
    const op = await observe(page, cdp, () => compareButton(target).click());
    expect(await page.evaluate(() => window.__GBDRAW_APP__.linearComparisonGlobalAction)).toBe(target);
    expect(op.undo).toBe(1);
    expectNoStatusWork(op, `comparison ${target}`);
  }
  await compareButton('losat').click();
  for (const [label, program] of [['LOSATP', 'blastp'], ['LOSATN', 'blastn']]) {
    const op = await observe(page, cdp, () => page.getByRole('group', { name: 'LOSAT Mode' })
      .getByRole('button', { name: label, exact: true }).click());
    expect(await page.evaluate(() => window.__GBDRAW_APP__.losatProgram)).toBe(program);
    expect(op.undo).toBe(1);
    expectNoStatusWork(op, `program ${label}`);
  }
  await compareButton('none').click();
  const prefix = page.locator('#output-prefix');
  await prefix.scrollIntoViewIfNeeded();
  await prefix.click();
  const key = await observe(page, cdp, () => page.keyboard.press('x'));
  await expect(prefix).toHaveValue('x');
  expectNoStatusWork(key, 'keystroke');
  const blur = await observe(page, cdp, () => page.keyboard.press('Tab'));
  expect(blur.undo).toBe(1);
  expectNoStatusWork(blur, 'blur');
  await expect(page.locator('[data-generation-application-feedback]')).toHaveCount(0);

  // CW-01 for presentation-only changes: Editor placement and help-tips.
  const toggle = page.locator('.drawer-toggle');
  for (const expanded of ['true', 'false']) {
    const op = await observe(page, cdp, () => toggle.click());
    await expect(toggle).toHaveAttribute('aria-expanded', expanded);
    expect(op.undo).toBe(0);
    expectNoStatusWork(op, `Editor aria-expanded=${expanded}`);
  }
  const help = page.locator('.help-tip > button[aria-describedby="linear-definition-lock-help"]');
  const tip = await observe(page, cdp, () => help.click());
  await expect(page.locator('[role="tooltip"]')).toContainText('Changes apply on Generate');
  expect(tip.undo).toBe(0);
  expectNoStatusWork(tip, 'help-tip');
});

test('a double press on Generate starts one operation', async ({ page }) => {
  test.setTimeout(180_000);
  await install(page);
  await openApp(page);
  await page.evaluate(async (content) => {
    const app = window.__GBDRAW_APP__;
    app.files.c_gb = new File([content], 'HmmtDNA.gbk', { type: 'text/plain' });
    await window.Vue.nextTick();
  }, source);
  const before = await getDiagramWorkerActivity(page);
  const undoBefore = await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
  await page.getByRole('button', { name: 'Generate Diagram', exact: true }).dblclick();
  await expect.poll(() => page.evaluate(() => window.__CW__.events.includes('generate.completed')),
    { timeout: 180_000 }).toBe(true);
  await settle(page);
  const after = await getDiagramWorkerActivity(page);
  expect(after.runs - before.runs).toBe(1);
  expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount()) - undoBefore).toBe(1);
  expect(await page.evaluate(() => window.__CW__.events
    .filter((name) => name === 'generate.processing-published').length)).toBe(1);
});
