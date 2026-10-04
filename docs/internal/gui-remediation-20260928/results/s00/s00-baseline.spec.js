// S00 baseline harness: evidence only, not a product test. Run through
// playwright.s00.config.js against one source snapshot served on its own
// origin. Timing runs and call-count runs are separate because V8 precise
// coverage changes execution cost.
'use strict';

const { test, expect } = require('@playwright/test');
const { appendFileSync, mkdirSync } = require('node:fs');
const { join, resolve } = require('node:path');

const repo = resolve(__dirname, '../../../../..');
const { openApp, getDiagramWorkerActivity, reveal } = require(join(repo, 'tests/web/helpers/app-lifecycle.cjs'));

const env = process.env;
const TARGET = env.S00_TARGET || 'unknown';
const OUT = env.S00_OUT_DIR || join(__dirname, 'out', TARGET);
const PAGES = Number(env.S00_PAGES || 3);
const WARMUP = Number(env.S00_WARMUP || 2);
const REPS = Number(env.S00_REPS || 8);
const GENERATE_WARM = Number(env.S00_GENERATE_WARM || 7);
const VNIG_PAGES = Number(env.S00_VNIG_PAGES || 3);
const COVERAGE = env.S00_COVERAGE === '1';
// Diagnosis only: frame attribution, heap per frame, mutation targets. Not used for budgets.
const DIAG = env.S00_DIAG === '1';
const QUIET_MS = 750;
const ACCEPT_TIMEOUT_MS = 120_000;
const SETTLE_TIMEOUT_MS = 60_000;
const GENERATE_TIMEOUT_MS = 20 * 60_000;
const INPUTS = {
  first: join(repo, 'tests/test_inputs/MG1655.gbk'),
  second: join(repo, 'tests/test_inputs/Sakai.gbk'),
  vnig: env.S00_VNIG_SESSION
};
const PREFIX_INPUT = '#output-prefix';
const TARGET_FUNCTIONS = [
  'getGenerationApplicationStatus', 'describeGenerationApplication', 'projectGenerationIntent',
  'compareGenerationIntent', 'liveGenerationTables', 'projectCanonicalRenderInput',
  'buildCanonicalRenderRequest', 'addGeneratedTableResources', 'buildLabelOverrideTsv',
  'buildLabelOverrideRows', 'buildFeatureMetadataMap', 'buildFeatureSelectorUniquenessIndexFromMetadata',
  'buildEditableLabelByFeatureId', 'buildFeatureIdsBySourceText', 'validateCurrentWriterActiveConfig',
  'generationComparisonPlan', 'generationMeaning', 'buildHistoryIntent', 'buildArtifactCheckpoint',
  'snapshotSignature', 'serializeFeatureVisibilityRules', 'buildDefaultColorOverrideTsv',
  'projectLinearComparisonUi', 'buildLinearComparisonTimeline'
];

mkdirSync(OUT, { recursive: true });
const write = (record) => appendFileSync(join(OUT, 'samples.jsonl'), `${JSON.stringify({ target: TARGET, ...record })}\n`);

// Page probe. Installed before app-lifecycle's Worker tracking, which wraps it.
const installProbe = (page) => page.addInitScript((diag) => {
  const s = window.__S00__ = {
    lastMutationAt: 0, mutationCount: 0, loaf: [], longtasks: [], events: [], inputs: [],
    workers: {}, historyRevisionAt: 0, lastRevision: undefined, processing: [], arm: null,
    diag, heap: [], mutationBatches: []
  };
  const capped = (list, entry) => { list.push(entry); if (list.length > 4000) list.splice(0, 2000); };
  const now = () => performance.now();
  const observe = (type, map, extra = {}) => {
    try {
      new PerformanceObserver((list) => { for (const entry of list.getEntries()) map(entry); })
        .observe({ type, buffered: true, ...extra });
    } catch { /* unsupported entry type is reported as missing data */ }
  };
  observe('long-animation-frame', (e) => s.loaf.push({ start: e.startTime, duration: e.duration,
    blocking: e.blockingDuration, end: e.startTime + e.duration,
    ...(diag ? { renderStart: e.renderStart, styleAndLayoutStart: e.styleAndLayoutStart,
      scripts: (e.scripts || []).map((x) => ({ invoker: x.invoker, invokerType: x.invokerType,
        fn: x.sourceFunctionName, url: String(x.sourceURL || '').split('/').slice(-2).join('/'),
        start: x.startTime, duration: x.duration, forced: x.forcedStyleAndLayoutDuration })) } : {}) }));
  observe('longtask', (e) => s.longtasks.push({ start: e.startTime, duration: e.duration, end: e.startTime + e.duration }));
  observe('event', (e) => s.events.push({ name: e.name, start: e.startTime, duration: e.duration,
    processingStart: e.processingStart, processingEnd: e.processingEnd, interactionId: e.interactionId || 0 }),
  { durationThreshold: 16 });
  for (const type of ['pointerdown', 'keydown', 'click', 'focusout']) {
    window.addEventListener(type, (event) => {
      if (event.isTrusted) s.inputs.push({ type, ts: event.timeStamp, key: event.key || '' });
    }, { capture: true });
  }
  const describe = (node) => (node?.nodeType === 1
    ? `${node.tagName.toLowerCase()}.${String(node.className?.baseVal ?? node.className ?? '').split(' ').slice(0, 2).join('.')}`
    : String(node?.nodeName));
  const startMutations = () => new MutationObserver((records) => {
    s.mutationCount += records.length;
    s.lastMutationAt = now();
    if (diag) capped(s.mutationBatches, { t: s.lastMutationAt, n: records.length,
      targets: records.slice(0, 3).map((r) => `${r.type}:${describe(r.target)}${r.attributeName ? `@${r.attributeName}` : ''}`) });
  }).observe(document.documentElement, { subtree: true, childList: true, attributes: true, characterData: true });
  if (document.documentElement) startMutations();
  else document.addEventListener('readystatechange', startMutations, { once: true });

  const NativeWorker = window.Worker;
  window.Worker = new Proxy(NativeWorker, {
    construct(target, args) {
      const worker = Reflect.construct(target, args, target);
      const name = String(args[0] || '').split('?')[0].split('/').pop();
      const entry = s.workers[name] ||= { constructed: 0, terminated: 0, posts: {} };
      entry.constructed += 1;
      worker.postMessage = function countedPostMessage(message, transfer) {
        const type = String(message?.type || typeof message);
        entry.posts[type] = (entry.posts[type] || 0) + 1;
        const post = Object.getPrototypeOf(worker).postMessage;
        return transfer === undefined ? post.call(this, message) : post.call(this, message, transfer);
      };
      worker.addEventListener('error', () => { entry.posts['worker-error'] = (entry.posts['worker-error'] || 0) + 1; });
      const terminate = worker.terminate.bind(worker);
      worker.terminate = () => { entry.terminated += 1; return terminate(); };
      return worker;
    }
  });

  const button = (label) => document.querySelector(`button[aria-label="${label}"]`);
  s.predicates = {
    compare: (target) => Boolean(button(target === 'losat' ? 'Run LOSAT for all adjacent pairs' : 'Set no comparison')
      ?.classList.contains('border-blue-500')),
    program: (label) => [...document.querySelectorAll('[data-linear-comparison-losat-mode] button')]
      .find((element) => element.textContent.trim() === label)?.getAttribute('aria-pressed') === 'true',
    inputValue: ({ selector, value }) => document.querySelector(selector)?.value === value,
    blurred: (selector) => document.activeElement !== document.querySelector(selector),
    generateAccepted: () => button('Generate Diagram')?.disabled === true
  };
  // Visible reflection: first frame whose DOM satisfies the predicate, then a
  // task posted from that frame's rAF, which runs after the frame is produced.
  s.armVisible = (kind, arg) => {
    const arm = s.arm = { kind, arg, armedAt: now(), frameAt: null, paintedAt: null };
    const tick = () => {
      if (s.arm !== arm || arm.frameAt !== null) return;
      if (s.predicates[kind](arg)) {
        arm.frameAt = now();
        const channel = new MessageChannel();
        channel.port1.onmessage = () => { arm.paintedAt = now(); };
        channel.port2.postMessage(0);
        return;
      }
      requestAnimationFrame(tick);
    };
    requestAnimationFrame(tick);
  };
  const frameLoop = () => {
    const app = window.__GBDRAW_APP__;
    const revision = window.__GBDRAW_HISTORY__?.revision?.value;
    if (revision !== s.lastRevision) { s.lastRevision = revision; s.historyRevisionAt = now(); }
    if (diag) capped(s.heap, [now(), performance.memory?.usedJSHeapSize ?? null]);
    const processing = Boolean(app?.processing);
    if (processing !== (s.processing.at(-1)?.value ?? false)) s.processing.push({ value: processing, at: now() });
    requestAnimationFrame(frameLoop);
  };
  requestAnimationFrame(frameLoop);
  s.settleState = () => {
    const app = window.__GBDRAW_APP__;
    const history = window.__GBDRAW_HISTORY__;
    const activity = window.__GBDRAW_DIAGRAM_WORKER_ACTIVITY__ || { instances: [] };
    const pendingWorkerRequests = activity.instances.reduce((total, instance) => total
      + instance.helpers.length + instance.runs.length
      - instance.settlements.filter(({ type }) => type !== 'init').length, 0);
    // An armed operation is not settled before its visible reflection, and the
    // quiet window starts no earlier than that reflection.
    const visiblePending = Boolean(s.arm && s.arm.paintedAt == null);
    const lastActivity = Math.max(s.lastMutationAt, s.historyRevisionAt, s.arm?.paintedAt ?? 0,
      ...s.loaf.map(({ end }) => end), ...s.longtasks.map(({ end }) => end), 0);
    return {
      busy: Boolean(app?.processing || app?.labelReflowProcessing || app?.sessionImportPending
        || history?.capturing?.value || history?.restoring?.value || pendingWorkerRequests > 0
        || visiblePending),
      pendingWorkerRequests, visiblePending, lastActivity, now: now()
    };
  };
}, DIAG);

const waitSettled = async (page, timeoutMs = SETTLE_TIMEOUT_MS) => {
  const started = Date.now();
  for (;;) {
    const state = await page.evaluate(() => window.__S00__.settleState());
    if (!state.busy && state.now - state.lastActivity >= QUIET_MS) return { ...state, settled: true };
    if (Date.now() - started > timeoutMs) return { ...state, settled: false };
    await page.waitForTimeout(50);
  }
};

const waitProcessingCycle = async (page, inputIndex, timeoutMs) => {
  const started = Date.now();
  for (;;) {
    const cycle = await page.evaluate((index) => {
      const s = window.__S00__;
      const start = s.inputs.slice(index).find(({ type }) => type === 'pointerdown' || type === 'keydown')?.ts;
      const published = s.processing.find(({ value, at }) => value && at > start);
      const cleared = published && s.processing.find(({ value, at }) => !value && at > published.at);
      return { publishedMs: published ? published.at - start : null, clearedMs: cleared ? cleared.at - start : null };
    }, inputIndex);
    if (cycle.clearedMs !== null) return { ...cycle, outcome: 'cleared' };
    const waited = Date.now() - started;
    if (cycle.publishedMs === null && waited > ACCEPT_TIMEOUT_MS) return { ...cycle, outcome: 'never-published' };
    if (waited > timeoutMs) return { ...cycle, outcome: 'timeout' };
    await page.waitForTimeout(100);
  }
};

const appState = (page) => page.evaluate((selector) => {
  const app = window.__GBDRAW_APP__;
  const history = window.__GBDRAW_HISTORY__;
  return {
    mode: app.mode,
    comparison: app.linearComparisonGlobalAction,
    losatProgram: app.losatProgram,
    prefix: document.querySelector(selector)?.value ?? null,
    processing: Boolean(app.processing),
    resultCount: app.results?.length || 0,
    resultRef: app.results?.[0]?.content?.length || 0,
    labelOverrideCount: Object.keys(app.labelTextFeatureOverrides || {}).length,
    visibilityOverrideCount: Object.keys(app.labelVisibilityOverrides || {}).length,
    undo: history?.getUndoCount?.() ?? null,
    redo: history?.getRedoCount?.() ?? null,
    error: String(app.errorLog?.summary || '')
  };
}, PREFIX_INPUT);

const coverageSession = async (page) => {
  if (!COVERAGE) return null;
  const cdp = await page.context().newCDPSession(page);
  await cdp.send('Profiler.enable');
  await cdp.send('Profiler.startPreciseCoverage', { callCount: true, detailed: false });
  return cdp;
};
// takePreciseCoverage resets counters, so each call returns counts since the previous call.
const takeCounts = async (cdp) => {
  if (!cdp) return null;
  const { result } = await cdp.send('Profiler.takePreciseCoverage');
  const counts = {};
  for (const script of result) {
    if (!script.url.includes('/gbdraw/web/js/')) continue;
    const file = script.url.split('/gbdraw/web/js/')[1].split('?')[0];
    for (const fn of script.functions) {
      const count = fn.ranges[0]?.count || 0;
      if (!count || !TARGET_FUNCTIONS.includes(fn.functionName)) continue;
      const key = `${file}#${fn.functionName}`;
      counts[key] = (counts[key] || 0) + count;
    }
  }
  return counts;
};

// One measured interaction. `trigger` must use real Playwright pointer or
// keyboard input; `expectState` asserts that the selection really changed.
const measure = async (page, cdp, meta, { predicate, arg, trigger, expectState, timeoutMs, terminalProcessing = false }) => {
  const pre = await waitSettled(page);
  expect(pre.settled, `pre-settle ${JSON.stringify(meta)}`).toBe(true);
  const before = await appState(page);
  const beforeWorkers = await getDiagramWorkerActivity(page);
  const marks = await page.evaluate(() => ({
    inputs: window.__S00__.inputs.length, workers: JSON.parse(JSON.stringify(window.__S00__.workers))
  }));
  if (cdp) await takeCounts(cdp);
  await page.evaluate(({ predicate, arg }) => window.__S00__.armVisible(predicate, arg), { predicate, arg });
  await trigger();
  // Generate is terminal only after processing was published and cleared; the
  // pre-processing before publication can be quiet for longer than QUIET_MS.
  const terminal = terminalProcessing ? await waitProcessingCycle(page, marks.inputs, timeoutMs) : null;
  const settle = await waitSettled(page, timeoutMs);
  const counts = await takeCounts(cdp);
  const after = await appState(page);
  const afterWorkers = await getDiagramWorkerActivity(page);
  const observed = await page.evaluate(({ inputIndex, settleAt }) => {
    const s = window.__S00__;
    const input = s.inputs.slice(inputIndex).find(({ type }) => type === 'pointerdown' || type === 'keydown');
    const start = input?.ts ?? null;
    const within = (entry) => start !== null && entry.end >= start && entry.start <= settleAt;
    const loaf = s.loaf.filter(within);
    const events = s.events.filter((entry) => start !== null && entry.start >= start - 1 && entry.start <= settleAt);
    const processingFalse = s.processing.find(({ value, at }) => !value && at > start);
    return {
      startInput: input?.type || null,
      visibleMs: s.arm?.paintedAt != null && start !== null ? s.arm.paintedAt - start : null,
      settleMs: start !== null ? settleAt - start : null,
      maxLoafMs: Math.max(0, ...loaf.map(({ duration }) => duration)),
      loafCount: loaf.length,
      loafAtLeast100: loaf.filter(({ duration }) => duration >= 100).length,
      loafTotalBlockingMs: loaf.reduce((total, { blocking }) => total + (blocking || 0), 0),
      maxEventDurationMs: Math.max(0, ...events.map(({ duration }) => duration)),
      processingEndMs: processingFalse && start !== null ? processingFalse.at - start : null,
      workers: s.workers
    };
  }, { inputIndex: marks.inputs, settleAt: settle.lastActivity });
  if (DIAG) {
    observed.lateFrames = await page.evaluate(({ inputIndex, settleAt }) => {
      const s = window.__S00__;
      const start = s.inputs.slice(inputIndex).find(({ type }) => type === 'pointerdown' || type === 'keydown')?.ts;
      const painted = s.arm?.paintedAt ?? start;
      if (start == null) return [];
      return s.loaf.filter((entry) => entry.start > painted && entry.start <= settleAt).map((entry) => ({
        afterInputMs: entry.start - start, duration: entry.duration, blocking: entry.blocking,
        renderStart: entry.renderStart ? entry.renderStart - entry.start : null,
        styleAndLayoutStart: entry.styleAndLayoutStart ? entry.styleAndLayoutStart - entry.start : null,
        scripts: entry.scripts,
        heap: s.heap.filter(([t]) => t >= entry.start - 120 && t <= entry.end + 60)
          .map(([t, bytes]) => [Math.round(t - entry.start), bytes]),
        mutations: s.mutationBatches.filter(({ t }) => t >= entry.start - 150 && t <= entry.end)
      }));
    }, { inputIndex: marks.inputs, settleAt: settle.lastActivity });
  }
  const sample = {
    ...meta, settled: settle.settled, terminal, ...observed, before, after, counts,
    diagramWorker: {
      constructionsDelta: afterWorkers.constructions - beforeWorkers.constructions,
      initializationsDelta: afterWorkers.initializations - beforeWorkers.initializations,
      helpersDelta: afterWorkers.helpers - beforeWorkers.helpers,
      runsDelta: afterWorkers.runs - beforeWorkers.runs,
      constructionsBefore: beforeWorkers.constructions
    },
    workerPostsBefore: marks.workers
  };
  write({ kind: 'sample', ...sample });
  expect(settle.settled, `settle ${JSON.stringify(meta)}`).toBe(true);
  if (terminal) expect(terminal.outcome, `terminal ${JSON.stringify(meta)}`).toBe('cleared');
  expect(observed.visibleMs, `visible ${JSON.stringify(meta)}`).not.toBeNull();
  expectState(before, after);
  return sample;
};

const openLinearWithGenomes = async (page) => {
  await page.getByRole('button', { name: 'Linear', exact: true }).click();
  for (let index = 1; index <= 2; index += 1) {
    const input = page.getByTestId(`linear-genbank-${index}`);
    if (await input.count() === 0) await page.getByRole('button', { name: 'Add sequence', exact: true }).last().click();
    await page.getByTestId(`linear-genbank-${index}`).setInputFiles(index === 1 ? INPUTS.first : INPUTS.second);
  }
  const settled = await waitSettled(page, 5 * 60_000);
  expect(settled.settled).toBe(true);
  expect((await appState(page)).mode).toBe('linear');
};

const clickCompare = (page, target) => page.getByRole('button', {
  name: target === 'losat' ? 'Run LOSAT for all adjacent pairs' : 'Set no comparison', exact: true
}).click();

const runInteractionSet = async (page, cdp, state, pageIndex) => {
  // Comparison switch: No comparison <-> LOSAT, both directions.
  for (let rep = 0; rep < WARMUP + REPS; rep += 1) {
    for (const target of ['losat', 'none']) {
      await measure(page, cdp, { op: 'compare', direction: `to-${target}`, state, pageIndex, rep, warmup: rep < WARMUP }, {
        predicate: 'compare', arg: target, trigger: () => clickCompare(page, target),
        expectState: (before, after) => {
          expect(before.comparison).not.toBe(target);
          expect(after.comparison).toBe(target);
        }
      });
    }
  }
  // LOSATN <-> LOSATP with LOSAT selected (setup click is not measured).
  await clickCompare(page, 'losat');
  for (let rep = 0; rep < WARMUP + REPS; rep += 1) {
    for (const [label, program] of [['LOSATP', 'blastp'], ['LOSATN', 'blastn']]) {
      await measure(page, cdp, { op: 'program', direction: `to-${label}`, state, pageIndex, rep, warmup: rep < WARMUP }, {
        predicate: 'program', arg: label,
        trigger: () => page.getByRole('group', { name: 'LOSAT Mode' }).getByRole('button', { name: label, exact: true }).click(),
        expectState: (before, after) => {
          expect(before.losatProgram).not.toBe(program);
          expect(after.losatProgram).toBe(program);
        }
      });
    }
  }
  await clickCompare(page, 'none');
  await runTextInputSet(page, cdp, state, pageIndex);
};

// Ordinary input: one keystroke (visible echo and its settlement), then blur
// by Tab, whose settlement includes the History commit.
const runTextInputSet = async (page, cdp, state, pageIndex) => {
  const input = page.locator(PREFIX_INPUT);
  await reveal(input);
  for (let rep = 0; rep < WARMUP + REPS; rep += 1) {
    await input.scrollIntoViewIfNeeded();
    await input.click();
    const current = await input.inputValue();
    const insert = rep % 2 === 0;
    const next = insert ? `${current}x` : current.slice(0, -1);
    const meta = { state, pageIndex, rep, warmup: rep < WARMUP };
    await measure(page, cdp, { ...meta, op: 'input-key', direction: insert ? 'insert' : 'delete' }, {
      predicate: 'inputValue', arg: { selector: PREFIX_INPUT, value: next },
      trigger: () => page.keyboard.press(insert ? 'x' : 'Backspace'),
      expectState: (_before, after) => expect(after.prefix).toBe(next)
    });
    await measure(page, cdp, { ...meta, op: 'input-blur', direction: insert ? 'insert' : 'delete' }, {
      predicate: 'blurred', arg: PREFIX_INPUT,
      trigger: () => page.keyboard.press('Tab'),
      expectState: (_before, after) => expect(after.prefix).toBe(next)
    });
  }
};

const generate = (page, cdp, meta) => measure(page, cdp, { op: 'generate', ...meta }, {
  predicate: 'generateAccepted', arg: null, timeoutMs: GENERATE_TIMEOUT_MS, terminalProcessing: true,
  trigger: () => page.getByRole('button', { name: 'Generate Diagram', exact: true }).click(),
  expectState: (before, after) => {
    expect(after.error, 'Generate error').toBe('');
    expect(after.processing).toBe(false);
    expect(after.resultCount).toBeGreaterThan(0);
  }
});

// Whole-genome Linear output has no editable labels, so the override case uses
// the feature editor's label-visibility action (a real label override owner),
// answering its global-label dialog with "this feature only".
const applyLabelOverride = async (page) => {
  const target = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const feature = app.extractedFeatures.find((candidate) => candidate?.svg_id
      && String(candidate.type || candidate.feature_type || '').toUpperCase() === 'CDS')
      || app.extractedFeatures.find((candidate) => candidate?.svg_id);
    if (!feature) throw new Error('No rendered feature was found.');
    app.openFeatureEditorFromList(feature, null);
    app.clickedFeature.labelVisibility = 'on';
    // Labels are globally off here, so the owner asks how to enable them.
    const update = app.updateClickedFeatureLabelText();
    for (let attempt = 0; attempt < 100 && !app.globalLabelModeDialog?.show; attempt += 1) {
      await new Promise((resolve) => setTimeout(resolve, 50));
    }
    if (app.globalLabelModeDialog?.show) app.handleGlobalLabelModeChoice('whitelist_only');
    await update;
    app.closeRightDrawer();
    return feature.svg_id;
  });
  expect((await waitSettled(page, 5 * 60_000)).settled).toBe(true);
  const state = await appState(page);
  expect(state.labelOverrideCount + state.visibilityOverrideCount).toBeGreaterThan(0);
  write({ kind: 'setup', step: 'label-visibility-override', featureId: target, state });
};

const openInstrumentedApp = async (page) => {
  page.on('dialog', (dialog) => dialog.accept());
  await installProbe(page);
  const cdp = await coverageSession(page);
  const lifecycle = await openApp(page, { checkErrors: false });
  write({ kind: 'environment', lifecycleShell: lifecycle.shell, browser: page.context().browser()?.version(),
    page: await page.evaluate(() => ({ crossOriginIsolated: self.crossOriginIsolated, userAgent: navigator.userAgent,
      hardwareConcurrency: navigator.hardwareConcurrency, devicePixelRatio, viewport: [innerWidth, innerHeight] })) });
  return cdp;
};

test.describe.configure({ mode: 'serial' });

for (let pageIndex = 0; pageIndex < (COVERAGE ? 1 : PAGES); pageIndex += 1) {
  test(`${TARGET} linear MG1655/Sakai page ${pageIndex}`, async ({ page }) => {
    test.setTimeout(90 * 60_000);
    const cdp = await openInstrumentedApp(page);
    await openLinearWithGenomes(page);
    write({ kind: 'setup', step: 'genomes-loaded', pageIndex, workers: await getDiagramWorkerActivity(page) });
    await runInteractionSet(page, cdp, 'pre-generate', pageIndex);
    await generate(page, cdp, { state: 'pre-generate', pageIndex, rep: 0, warmup: false, sequence: 'first' });
    await runInteractionSet(page, cdp, 'post-generate', pageIndex);
    for (let rep = 0; rep < GENERATE_WARM; rep += 1) {
      await generate(page, cdp, { state: 'post-generate', pageIndex, rep, warmup: false, sequence: 'repeat' });
    }
    await applyLabelOverride(page);
    await runInteractionSet(page, cdp, 'post-generate-override', pageIndex);
    await generate(page, cdp, { state: 'post-generate-override', pageIndex, rep: 0, warmup: false, sequence: 'repeat' });
  });
}

for (let pageIndex = 0; pageIndex < (COVERAGE ? 1 : VNIG_PAGES); pageIndex += 1) {
  test(`${TARGET} Vnig saved preview page ${pageIndex}`, async ({ page }) => {
    test.skip(!INPUTS.vnig, 'S00_VNIG_SESSION is not set.');
    test.setTimeout(90 * 60_000);
    const cdp = await openInstrumentedApp(page);
    await page.locator('input[type="file"][accept*="application/json"][accept*="application/gzip"]')
      .setInputFiles(INPUTS.vnig);
    await page.waitForFunction(() => window.__GBDRAW_APP__?.sessionImportPending === false
      && window.__GBDRAW_APP__?.results?.length > 0, null, { timeout: 10 * 60_000 });
    expect((await waitSettled(page, 5 * 60_000)).settled).toBe(true);
    write({ kind: 'setup', step: 'vnig-loaded', pageIndex, state: await appState(page),
      workers: await getDiagramWorkerActivity(page) });
    await runTextInputSet(page, cdp, 'vnig-saved-preview', pageIndex);
    await generate(page, cdp, { state: 'vnig-saved-preview', pageIndex, rep: 0, warmup: false, sequence: 'first' });
    await generate(page, cdp, { state: 'vnig-generated', pageIndex, rep: 0, warmup: false, sequence: 'repeat' });
  });
}

test.afterAll(() => {
  write({ kind: 'run-end', inputs: Object.fromEntries(Object.entries(INPUTS).map(([key, path]) => [key, path || null])),
    settings: { PAGES, WARMUP, REPS, GENERATE_WARM, VNIG_PAGES, COVERAGE, QUIET_MS } });
});
