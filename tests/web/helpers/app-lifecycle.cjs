const { readFileSync } = require('node:fs');
const { join } = require('node:path');
const { expect, test: baseTest } = require('@playwright/test');
const { PYODIDE_SNAPSHOT_TEST_PATH } = require('./pyodide-snapshot-test-server.cjs');

// The current canonical request schema, Session version and feature catalog
// schema, from the constants generated from their Python owners.
const generatedSchema = readFileSync(
  join(__dirname, '..', '..', '..', 'gbdraw', 'web', 'js', 'contracts', 'schema.generated.js'), 'utf8'
);
const generatedConstant = name => Number(generatedSchema.match(new RegExp(`^export const ${name} = (\\d+);$`, 'm'))[1]);
const CURRENT_REQUEST_SCHEMA = generatedConstant('CANONICAL_REQUEST_SCHEMA');
const CURRENT_SESSION_VERSION = generatedConstant('SESSION_VERSION');
const CURRENT_FEATURE_CATALOG_SCHEMA = generatedConstant('FEATURE_CATALOG_SCHEMA');

const DEFAULT_APP_TIMEOUT_MS = 180_000;
const DEFAULT_APP_PATH = '/gbdraw/web/index.html';
const pageDiagnostics = new WeakMap();
const workerTrackingPages = new WeakSet();
// The pages of `warmAppTest` (below): page -> { contextOptions, baseline, ready, savedSessions }.
const warmPages = new WeakMap();

const compactJson = (value, maxLength = 6_000) => {
  let rendered;
  try {
    rendered = JSON.stringify(value, null, 2);
  } catch (error) {
    rendered = JSON.stringify({ diagnosticSerializationError: String(error?.message || error) });
  }
  if (rendered.length <= maxLength) return rendered;
  return `${rendered.slice(0, maxLength)}\n... diagnostics truncated ...`;
};

const installPageErrorCollection = (page) => {
  if (pageDiagnostics.has(page)) return pageDiagnostics.get(page);
  const diagnostics = { pageErrors: [], consoleErrors: [] };
  page.on('pageerror', (error) => diagnostics.pageErrors.push(String(error?.message || error)));
  page.on('console', (message) => {
    if (message.type() === 'error') diagnostics.consoleErrors.push(message.text());
  });
  pageDiagnostics.set(page, diagnostics);
  return diagnostics;
};

const installDiagramWorkerTracking = async (page) => {
  if (workerTrackingPages.has(page)) return;
  await page.addInitScript(() => {
    const activity = {
      constructions: 0,
      instances: []
    };
    window.__GBDRAW_DIAGRAM_WORKER_ACTIVITY__ = activity;

    const NativeWorker = window.Worker;
    const nativePostMessage = NativeWorker.prototype.postMessage;
    NativeWorker.prototype.postMessage = function lifecycleTrackedPostMessage(message, transfer) {
      const instance = this.__gbdrawLifecycleActivity;
      if (instance) {
        const transferList = Array.isArray(transfer) ? transfer : [];
        if (message?.type === 'init') {
          instance.initializations += 1;
          instance.events.push('init:request');
        } else if (message?.type === 'helper') {
          instance.helpers.push({
            requestId: String(message.requestId ?? ''),
            operation: String(message.operation || ''),
            transferCount: transferList.length,
            transferredBytes: transferList.reduce(
              (total, item) => total + Number(item?.byteLength || 0),
              0
            )
          });
          instance.events.push(`helper:request:${String(message.operation || '')}`);
        } else if (message?.type === 'run') {
          const resourceManifest = Array.isArray(message.payload?.resourceManifest)
            ? message.payload.resourceManifest
            : [];
          const stagedResources = Array.isArray(message.payload?.stagedResources)
            ? message.payload.stagedResources
            : [];
          instance.runs.push({
            requestId: String(message.requestId ?? ''),
            referencedResourceCount: resourceManifest.length,
            referencedDeclaredBytes: resourceManifest.reduce(
              (total, resource) => total + Number(resource?.size || 0),
              0
            ),
            stagedResourceCount: stagedResources.length,
            stagedResourceBytes: stagedResources.reduce(
              (total, resource) => total + Number(resource?.bytes?.byteLength || 0),
              0
            ),
            transferCount: transferList.length,
            transferredBytes: transferList.reduce(
              (total, item) => total + Number(item?.byteLength || 0),
              0
            ),
            hasBase64ResourceTable: Boolean(message.payload?.resources)
          });
          instance.events.push('run:request');
        }
      }
      if (transfer === undefined) return nativePostMessage.call(this, message);
      return nativePostMessage.call(this, message, transfer);
    };

    const nativeTerminate = NativeWorker.prototype.terminate;
    NativeWorker.prototype.terminate = function lifecycleTrackedTerminate() {
      const instance = this.__gbdrawLifecycleActivity;
      if (instance) {
        instance.terminated = true;
        instance.events.push('terminate');
      }
      return nativeTerminate.call(this);
    };

    window.Worker = new Proxy(NativeWorker, {
      construct(target, args) {
        const worker = Reflect.construct(target, args, target);
        const url = String(args[0] || '');
        if (!url.includes('diagram-generation-worker.js')) return worker;

        activity.constructions += 1;
        const instance = {
          id: activity.constructions,
          url,
          initializations: 0,
          helpers: [],
          runs: [],
          settlements: [],
          errors: [],
          events: [],
          terminated: false
        };
        activity.instances.push(instance);
        worker.__gbdrawLifecycleActivity = instance;

        worker.addEventListener('message', (event) => {
          const message = event.data || {};
          if (!['init', 'helper', 'run'].includes(message.type)) return;
          // A bounded reply streams parts before its one final settlement.
          if (message.status === 'part') {
            instance.events.push(`${message.type}:part`);
            return;
          }
          const identifier = message.type === 'init' ? message.id : message.requestId;
          instance.settlements.push({
            type: message.type,
            id: String(identifier ?? ''),
            ok: message.ok === true,
            error: message.ok === false
              ? String(message.error?.message || message.error || '')
              : ''
          });
          instance.events.push(`${message.type}:${message.ok === true ? 'ok' : 'error'}`);
        });
        worker.addEventListener('error', (event) => {
          instance.errors.push({
            type: 'error',
            message: String(event?.message || 'Diagram Worker error')
          });
          instance.events.push('worker:error');
        });
        worker.addEventListener('messageerror', () => {
          instance.errors.push({
            type: 'messageerror',
            message: 'Diagram Worker message could not be decoded'
          });
          instance.events.push('worker:messageerror');
        });

        return worker;
      }
    });
  });
  workerTrackingPages.add(page);
};

const summarizeDiagramWorkerActivity = (activity = {}) => {
  const instances = Array.isArray(activity.instances) ? activity.instances : [];
  const settlements = instances.flatMap((instance) => instance.settlements || []);
  return {
    constructions: Number(activity.constructions || 0),
    initializations: instances.reduce(
      (total, instance) => total + Number(instance.initializations || 0),
      0
    ),
    helpers: instances.reduce(
      (total, instance) => total + (Array.isArray(instance.helpers) ? instance.helpers.length : 0),
      0
    ),
    runs: instances.reduce(
      (total, instance) => total + (Array.isArray(instance.runs) ? instance.runs.length : 0),
      0
    ),
    settledInitializations: settlements.filter(({ type }) => type === 'init').length,
    settledHelpers: settlements.filter(({ type }) => type === 'helper').length,
    settledRuns: settlements.filter(({ type }) => type === 'run').length,
    instances
  };
};

const getDiagramWorkerActivity = async (page) => summarizeDiagramWorkerActivity(
  await page.evaluate(() => window.__GBDRAW_DIAGRAM_WORKER_ACTIVITY__ || {
    constructions: 0,
    instances: []
  })
);

const getAppShellSnapshot = (page) => page.evaluate(() => {
  const app = window.__GBDRAW_APP__;
  const mainRuntimeFields = [
    'pyodide',
    'pyodideReady',
    'pyodideLoading',
    'pyodideError',
    'pyodideStatus'
  ];
  return {
    appMounted: Boolean(app),
    paletteDefinitionCount: app && typeof app.paletteDefinitions === 'object'
      ? Object.keys(app.paletteDefinitions).length
      : 0,
    mainLoaderPresent: typeof window.loadPyodide === 'function',
    mainRuntimeFields: app
      ? mainRuntimeFields.filter((field) => Object.prototype.hasOwnProperty.call(app, field))
      : []
  };
});

const getLifecycleDiagnostics = async (page) => {
  const collected = installPageErrorCollection(page);
  return {
    shell: await getAppShellSnapshot(page),
    worker: await getDiagramWorkerActivity(page),
    pageErrors: [...collected.pageErrors],
    consoleErrors: [...collected.consoleErrors]
  };
};

const waitForAppShell = async (
  page,
  { waitForPalette = true, timeout = DEFAULT_APP_TIMEOUT_MS } = {}
) => {
  await page.waitForFunction(() => Boolean(window.__GBDRAW_APP__), null, { timeout });
  if (waitForPalette) {
    await page.waitForFunction(
      () => Object.keys(window.__GBDRAW_APP__?.paletteDefinitions || {}).length > 0,
      null,
      { timeout }
    );
  }
};

const assertAppShellReady = async (
  page,
  { waitForPalette = true, checkErrors = true } = {}
) => {
  const diagnostics = await getLifecycleDiagnostics(page);
  const rendered = compactJson(diagnostics);
  expect(diagnostics.shell.appMounted, rendered).toBe(true);
  if (waitForPalette) {
    expect(diagnostics.shell.paletteDefinitionCount, rendered).toBeGreaterThan(0);
  }
  expect(diagnostics.shell.mainLoaderPresent, rendered).toBe(false);
  expect(diagnostics.shell.mainRuntimeFields, rendered).toEqual([]);
  if (checkErrors) {
    expect(diagnostics.pageErrors, rendered).toEqual([]);
    expect(diagnostics.consoleErrors, rendered).toEqual([]);
  }
  return diagnostics;
};

const openApp = async (
  page,
  {
    path = DEFAULT_APP_PATH,
    waitForPalette = true,
    timeout = DEFAULT_APP_TIMEOUT_MS,
    checkErrors = true
  } = {}
) => {
  installPageErrorCollection(page);
  await installDiagramWorkerTracking(page);
  // A warm page that `warmAppTest` restored is the app as it booted.
  const warm = warmPages.get(page);
  const freshBootPath = path === DEFAULT_APP_PATH && waitForPalette;
  const restored = Boolean(warm?.ready) && freshBootPath;
  if (warm) warm.ready = false;
  if (restored) return assertAppShellReady(page, { waitForPalette, checkErrors });
  await page.goto(path, { waitUntil: 'domcontentloaded' });
  await waitForAppShell(page, { waitForPalette, timeout });
  const diagnostics = await assertAppShellReady(page, { waitForPalette, checkErrors });
  if (warm) warm.baseline = freshBootPath ? await readWarmBaseline(page) : null;
  return diagnostics;
};

// A loaded Session shows the Linear File cards that hold a file folded to
// their header line. A case that acts inside the cards after a load opens
// them first; `reveal` opens the card around one locator.
const openLinearFileCards = (page) => page.evaluate(() => {
  document.querySelectorAll('[data-linear-source-fold]').forEach((card) => { card.open = true; });
});

// Enter toggles a summary wherever its help-tip button sits; a click at the
// summary center can land on that button and leave the section closed.
const reveal = async (locator) => {
  for (const details of await locator.locator('xpath=ancestor::details').all()) {
    if (await details.getAttribute('open') === null) {
      await details.locator(':scope > summary').press('Enter');
    }
  }
  return locator;
};

const assertDiagramWorkerIdle = async (page, label = 'Expected the diagram Worker to remain idle') => {
  const diagnostics = await getLifecycleDiagnostics(page);
  expect(
    {
      constructions: diagnostics.worker.constructions,
      initializations: diagnostics.worker.initializations,
      helpers: diagnostics.worker.helpers,
      runs: diagnostics.worker.runs
    },
    `${label}:\n${compactJson(diagnostics)}`
  ).toEqual({ constructions: 0, initializations: 0, helpers: 0, runs: 0 });
  return diagnostics.worker;
};

const assertSessionLoadLeftWorkerIdle = (page) => assertDiagramWorkerIdle(
  page,
  'Loading a saved preview must not initialize the diagram Worker'
);

// Chromium can collect a protocol promise that page.evaluate awaits for a long
// time (Playwright then reports "Execution context was destroyed"). Start the
// operation once, retain its outcome in the page, poll for settlement, and
// rethrow a rejection. tests/web/playwright-long-app-promises.test.mjs keeps
// long app operations (Generate, Session import and save) on this path.
let retainedEvaluationIndex = 0;
const evaluateWithRetainedPromise = async (page, callback, argument) => {
  const key = `__GBDRAW_TEST_EVALUATION_${++retainedEvaluationIndex}`;
  try {
    await page.evaluate(`(() => {
      const entry = window[${JSON.stringify(key)}] = { settled: false };
      entry.promise = Promise.resolve((${callback.toString()})(${JSON.stringify(argument) ?? 'undefined'}));
      entry.promise.then(value => {
        entry.value = value;
        entry.settled = true;
      }, error => {
        entry.error = error;
        entry.failed = true;
        entry.settled = true;
      });
    })()`);
    // Like page.evaluate, the owning test deadline bounds this operation.
    await page.waitForFunction(key => window[key]?.settled, key, { timeout: 0 });
    return await page.evaluate(key => {
      const entry = window[key];
      if (entry.failed) throw entry.error;
      return entry.value;
    }, key);
  } finally {
    // Best-effort cleanup must not replace the operation's own failure.
    if (!page.isClosed()) await page.evaluate(key => { delete window[key]; }, key).catch(() => {});
  }
};

// No Generate, label reflow, automatic rerender, rule matching, Session import,
// or History step is running, over two frames in a row.
const settleLive = async (page, { timeout = 120_000 } = {}) => {
  const idle = () => {
    const app = window.__GBDRAW_APP__;
    const history = window.__GBDRAW_HISTORY__;
    return !app.processing && !app.labelReflowProcessing && !app.sessionImportPending
      && !app.ruleMatchingPending && !history?.capturing?.value && !history?.restoring?.value
      && !history?.mutationPending?.();
  };
  for (let stableFrames = 0; stableFrames < 2;) {
    await page.waitForFunction(idle, null, { timeout });
    await page.evaluate(() => new Promise((resolve) => (
      requestAnimationFrame(() => requestAnimationFrame(resolve))
    )));
    stableFrames = await page.evaluate(idle) ? stableFrames + 1 : 0;
  }
};

// G-G(3), Web GUI audit 2026-09-30: the shared Generate and Save helpers
// reject three outcomes that previously reached users without failing a test.
// A caller may opt out only by naming the audit ID of the known defect, so the
// PR that fixes that defect finds and removes the opt-out.
const AUDIT_ID = /^(?:[A-Z]{1,3}-\d{2}|N-\d{2})$/;
const knownDefectOptOut = (value, option) => {
  if (value === false || value === undefined || value === null) return '';
  if (typeof value === 'string' && AUDIT_ID.test(value)) return value;
  throw new TypeError(`${option} must name the audit ID of the known defect, for example 'IN-08'.`);
};

const INVALID_RUN_INFO_TEXT = /\b(?:NaN|undefined)\b/;
const findInvalidRunInfoValues = (value, path = 'lastRunInfo', found = []) => {
  if (typeof value === 'number') {
    if (!Number.isFinite(value)) found.push(`${path}=${value}`);
  } else if (typeof value === 'string') {
    if (INVALID_RUN_INFO_TEXT.test(value)) found.push(`${path}=${JSON.stringify(value.slice(0, 200))}`);
  } else if (Array.isArray(value)) {
    value.forEach((item, index) => findInvalidRunInfoValues(item, `${path}[${index}]`, found));
  } else if (value && typeof value === 'object') {
    Object.entries(value).forEach(([key, item]) => findInvalidRunInfoValues(item, `${path}.${key}`, found));
  }
  return found;
};

// An operation that does not clear the alert must not inherit an earlier one.
const readErrorSignature = (page) => page.evaluate(() => {
  const error = window.__GBDRAW_APP__?.errorLog;
  return error ? JSON.stringify([error.code, error.operation, error.stage, error.summary]) : null;
});

const assertOperationHealth = async (page, {
  operation = 'operation',
  allowUnknown = false,
  allowPageErrors = false,
  allowInvalidRunInfo = false,
  errorSignatureBefore = undefined
} = {}) => {
  const unknownDefect = knownDefectOptOut(allowUnknown, 'allowUnknown');
  const pageErrorDefect = knownDefectOptOut(allowPageErrors, 'allowPageErrors');
  const runInfoDefect = knownDefectOptOut(allowInvalidRunInfo, 'allowInvalidRunInfo');
  const observed = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const error = app?.errorLog;
    return {
      errorCode: error?.code ?? null,
      errorSummary: String(error?.summary || ''),
      errorSignature: error ? JSON.stringify([error.code, error.operation, error.stage, error.summary]) : null,
      lastRunInfo: app?.lastRunInfo ?? null
    };
  });
  const inheritedError = errorSignatureBefore !== undefined
    && observed.errorSignature !== null
    && observed.errorSignature === errorSignatureBefore;
  const collected = installPageErrorCollection(page);
  const health = {
    errorCode: inheritedError ? null : observed.errorCode,
    pageErrors: [...collected.pageErrors],
    invalidRunInfo: findInvalidRunInfoValues(observed.lastRunInfo)
  };
  const rendered = compactJson({ operation, ...health, errorSummary: observed.errorSummary });
  if (!unknownDefect) {
    expect(health.errorCode, `${operation} must not report an UNKNOWN diagnostic:\n${rendered}`)
      .not.toBe('UNKNOWN');
  }
  if (!pageErrorDefect) {
    expect(health.pageErrors, `${operation} must not raise an uncaught page error:\n${rendered}`)
      .toEqual([]);
  }
  if (!runInfoDefect) {
    expect(health.invalidRunInfo, `${operation} must not record NaN or undefined in Run Info:\n${rendered}`)
      .toEqual([]);
  }
  return health;
};

const generateAndWaitForResult = async (
  page,
  {
    expectedStatus = 'ok',
    requireCommittedResult = expectedStatus === 'ok',
    allowUnknown = false,
    allowPageErrors = false,
    allowInvalidRunInfo = false
  } = {}
) => {
  const outcome = await evaluateWithRetainedPromise(page, async () => {
    const app = window.__GBDRAW_APP__;
    const result = await app.runAnalysis();
    return {
      result,
      errorSummary: String(app.errorLog?.summary || ''),
      errorDetails: Array.isArray(app.errorLog?.details)
        ? app.errorLog.details.map((detail) => String(detail))
        : [],
      resultCount: Array.isArray(app.results) ? app.results.length : 0,
      committedResult: Boolean(
        Array.isArray(app.results)
        && app.results.length > 0
        && String(app.results[0]?.content || '').includes('<svg')
      )
    };
  });
  const diagnostics = await getLifecycleDiagnostics(page);
  const rendered = compactJson({ outcome, diagnostics });
  expect(['ok', 'error', 'canceled', 'stale'], rendered).toContain(outcome.result?.status);
  if (expectedStatus !== null) expect(outcome.result?.status, rendered).toBe(expectedStatus);
  if (requireCommittedResult) expect(outcome.committedResult, rendered).toBe(true);
  if (outcome.result?.status === 'error') {
    expect(Boolean(outcome.errorSummary || outcome.errorDetails.length), rendered).toBe(true);
  }
  outcome.health = await assertOperationHealth(page, {
    operation: `Generate (${outcome.result?.status})`,
    allowUnknown,
    allowPageErrors,
    allowInvalidRunInfo
  });
  return outcome;
};

// G-C (Web GUI audit 2026-09-30): an operation that is not an edit (Result
// selection, a no-change Generate, a mode round trip, an Undo+Redo pair, an
// unrelated toggle) must leave user-owned state unchanged. The snapshot follows
// the audit harness: config, UI and editor state plus the user-owned feature,
// label, and group maps; the feature catalog, bulk feature tables, and preview
// navigation are excluded. `mode` reads that mode's drawing instead of the
// shown one.
const snapshotUserOwnedState = (page, { mode = null } = {}) => page.evaluate(async (requestedMode) => {
  const { state } = await import('/gbdraw/web/js/state.js');
  const config = await import('/gbdraw/web/js/services/config.js');
  const drawingMode = requestedMode || state.mode.value;
  const drawing = state.drawingOf(drawingMode);
  const editor = config.buildEditorStateData(drawing);
  delete editor.featureCatalog;
  const features = config.buildFeatureStateData(drawing);
  for (const key of ['extractedFeatures', 'biologicalFeatures']) {
    delete features[key];
  }
  const orthogroups = config.buildOrthogroupStateData(drawing);
  orthogroups.groupCount = orthogroups.groups.length;
  delete orthogroups.groups;
  const ui = config.buildUiStateData(drawing, { includePreviewNavigation: false });
  // A drawing reads and saves only its own mode's layout slot (state.js,
  // `modes.<mode>.ui.layoutPreferences`); the other mode's slot is never read.
  ui.layoutPreferences = ui.layoutPreferences?.[drawingMode];
  const history = window.__GBDRAW_HISTORY__;
  return JSON.parse(JSON.stringify({
    state: {
      mode: state.mode.value,
      config: config.buildConfigData(drawing),
      ui,
      editor,
      features,
      orthogroups
    },
    history: { undo: history?.getUndoCount?.() ?? null, redo: history?.getRedoCount?.() ?? null }
  }));
}, mode);

const diffUserOwnedState = (before, after, path = '', changes = []) => {
  if (JSON.stringify(before) === JSON.stringify(after)) return changes;
  const isObject = (value) => value && typeof value === 'object' && !Array.isArray(value);
  if (isObject(before) && isObject(after)) {
    for (const key of [...new Set([...Object.keys(before), ...Object.keys(after)])].sort()) {
      diffUserOwnedState(
        Object.hasOwn(before, key) ? before[key] : '<missing>',
        Object.hasOwn(after, key) ? after[key] : '<missing>',
        path ? `${path}.${key}` : key,
        changes
      );
    }
    return changes;
  }
  if (Array.isArray(before) && Array.isArray(after) && before.length === after.length) {
    before.forEach((item, index) => diffUserOwnedState(item, after[index], `${path}[${index}]`, changes));
    return changes;
  }
  changes.push({ path, before, after });
  return changes;
};

const observeNonEditOperation = async (page, operation) => {
  const before = await snapshotUserOwnedState(page);
  await operation();
  const after = await snapshotUserOwnedState(page);
  return {
    before,
    after,
    changes: diffUserOwnedState(before.state, after.state),
    undoDelta: after.history.undo - before.history.undo
  };
};

// allowedPaths name the settings an operation legitimately owns. A change that
// the operation records as a History step is accepted only with allowRecorded.
const expectNoSilentStateChange = (observation, { label, allowedPaths = [], allowRecorded = false }) => {
  const allowed = (path) => allowedPaths.some((prefix) => (
    path === prefix || path.startsWith(`${prefix}.`) || path.startsWith(`${prefix}[`)
  ));
  const unexpected = observation.changes.filter(({ path }) => !allowed(path));
  if (allowRecorded && observation.undoDelta > 0) return unexpected;
  expect(unexpected, `${label} silently changed user-owned state:\n${compactJson(unexpected)}`).toEqual([]);
  return unexpected;
};

// The warm app page (CD-05, CI cost 2026-10-14). Booting the app and the
// diagram Worker's Python runtime is most of a short browser case. A spec that
// takes its `test` from `warmAppTest` keeps one app page per Playwright
// worker. Before each case the `page` fixture returns that page to the state a
// fresh boot shows and compares it with what the page held right after it
// booted; a difference, or a restore that fails, gives the case a fresh
// context instead. On a restored page `openApp` does not navigate. Each case
// stays its own test, and its `warm-app-page` annotation says whether it
// reused the page or why it did not. The restore keeps the diagram Worker with
// its Python runtime and resource caches: every upload is a new resource
// identity, so the next Generate evicts the earlier entries, as it does when a
// user replaces a file.
const WARM_APP_PAGE = 'warm-app-page';
const WARM_RESTORE_TIMEOUT_MS = 60_000;
// History keeps at most 30 steps (services/history.js); a restore that undoes
// more found a step that does not leave.
const WARM_RESTORE_MAX_UNDO = 40;
// A Generate writes these globals on every run (gbdraw/web/js/app/run-analysis.js);
// a page that has not generated has none. A global missing here costs speed,
// not correctness: every case after a Generate falls back, naming the global.
const GENERATE_GLOBALS = ['__GBDRAW_LAST_LOSAT_TELEMETRY__', '_origPairwiseMin', '_origPairwiseMax'];
// A Session file name (services/config.js `buildSessionFilename`). The app
// remembers the Session file it downloaded last and asks before it downloads
// one of that name again; the restore cannot reach that memory, so a page
// that saved a Session after it booted serves no later case (OV-1001).
const SESSION_FILE = /session\.json\.gz$/;

// What the restore compares: both drawings' user-owned state, the inputs and
// Results, the Session files saved, and what a case can leave in the page
// outside the app.
const readWarmFacts = async (page) => ({
  url: page.url(),
  viewport: page.viewportSize(),
  savedSessions: [...warmPages.get(page).savedSessions],
  circular: (await snapshotUserOwnedState(page, { mode: 'circular' })),
  linear: (await snapshotUserOwnedState(page, { mode: 'linear' })),
  page: await page.evaluate(async () => {
    const { state } = await import('/gbdraw/web/js/state.js');
    const app = window.__GBDRAW_APP__;
    const describe = (value) => {
      if (value instanceof Blob) return `${value.name || 'Blob'} (${value.size} bytes)`;
      if (Array.isArray(value)) return value.map(describe);
      if (value && typeof value === 'object') {
        return Object.fromEntries(Object.entries(value).map(([key, item]) => [key, describe(item)]));
      }
      return value ?? null;
    };
    return {
      results: app.results.length,
      extractedFeatures: app.extractedFeatures.length,
      error: app.errorLog?.code ?? null,
      notice: document.querySelector('[data-app-notice]')?.textContent.trim() || '',
      choicePrompt: Boolean(document.querySelector('[data-choice-prompt]')),
      zoom: app.zoom,
      pan: { ...app.canvasPan },
      inputs: describe({ ...state.cardFiles, ...state.drawingOfMode('circular').files }),
      linearSources: state.drawingOfMode('linear').linearSeqs.map(({ gb, gff, fasta, depth }) => describe({ gb, gff, fasta, depth })),
      openSections: [...document.querySelectorAll('details[open] > summary')].map((summary) => summary.textContent.trim()),
      globals: Object.fromEntries(Object.keys(window).sort().map((key) => [key, typeof window[key]])),
      testHooks: Object.entries(window.__GBDRAW_TEST_HOOKS__ || {}).map(([name, hook]) => `${name}: ${hook}`).sort(),
      storage: [Object.entries(localStorage).sort(), Object.entries(sessionStorage).sort()]
    };
  })
});

// The fresh boot: the facts above, and the values the restore writes back.
const readWarmBaseline = async (page) => ({
  facts: await readWarmFacts(page),
  restorePoint: await page.evaluate(async () => {
    const { state } = await import('/gbdraw/web/js/state.js');
    const config = await import('/gbdraw/web/js/services/config.js');
    return JSON.parse(JSON.stringify({
      shownDrawing: state.activeDrawing().id,
      drawings: Object.fromEntries(['circular', 'linear'].map((mode) => [mode, {
        config: config.buildConfigData(state.drawingOf(mode)),
        ui: config.buildUiStateData(state.drawingOf(mode), { includePreviewNavigation: false })
      }])),
      // The test-owned lists an init script starts (`__AUDIT_LIFECYCLE__`).
      testLists: Object.fromEntries(Object.keys(window)
        .filter((key) => key.startsWith('__') && Array.isArray(window[key]))
        .map((key) => [key, window[key]])),
      storage: [{ ...localStorage }, { ...sessionStorage }]
    }));
  })
});

const returnToFreshBoot = async (page, { facts, restorePoint }) => {
  if (page.url() !== facts.url) return [`the case left the app for ${page.url()}`];
  // A focused field holds an open History step; leaving it commits the step.
  await page.evaluate(() => document.activeElement?.blur?.());
  await settleLive(page, { timeout: WARM_RESTORE_TIMEOUT_MS });
  // History takes back the inputs, Results and edits its steps recorded.
  // Each drawing has its own History and a drawing switch is no step
  // (PD-OI-090), so the restore walks each drawing's History, the drawing
  // shown at boot last, which stays shown.
  const drawingIds = await page.evaluate(async () => (await import('/gbdraw/web/js/state.js')).state
    .drawingList().map((drawing) => drawing.id));
  const shownAtBoot = restorePoint.shownDrawing;
  for (const drawingId of [...drawingIds.filter((id) => id !== shownAtBoot), shownAtBoot]) {
    const switched = await page.evaluate((id) => window.__GBDRAW_APP__.selectDrawing(id), drawingId);
    if (switched && switched.status !== 'ok') return [`the restore could not show the drawing ${drawingId}`];
    await settleLive(page);
    for (let undone = 0; undone < WARM_RESTORE_MAX_UNDO; undone += 1) {
      const steps = await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
      if (steps === 0) break;
      await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
      await settleLive(page);
      if (await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount()) >= steps) break;
    }
  }
  // What a case set outside History: Linear sources set through their setter,
  // drafts and settings (Reset Settings), and the settings of both drawings as
  // the boot read them. The UI part follows the Multi-Record Canvas watcher,
  // which writes layout defaults when the reset turns the canvas back on.
  await evaluateWithRetainedPromise(page, async ({ drawings, shownAtBoot }) => {
    const { state } = await import('/gbdraw/web/js/state.js');
    const config = await import('/gbdraw/web/js/services/config.js');
    const app = window.__GBDRAW_APP__;
    // The source setter acts on the shown drawing, so the Linear drawing is
    // shown while its sources are removed.
    const linear = state.drawingOfMode('linear');
    const linearSeqs = linear.linearSeqs;
    if (linearSeqs.some((seq) => seq.gb || seq.gff || seq.fasta)) {
      app.selectDrawing(linear.id);
      for (let index = linearSeqs.length - 1; index >= 0; index -= 1) {
        for (const field of ['gb', 'gff', 'fasta']) {
          if (linearSeqs[index]?.[field]) await app.setLinearSeqPrimaryFile(index, field, null);
        }
      }
      app.selectDrawing(shownAtBoot);
      await window.Vue.nextTick();
    }
    await app.resetSettings();
    for (const [mode, { config: settings }] of Object.entries(drawings)) {
      config.applyCapturedConfigData(state.drawingOf(mode), settings, { resolveTrackPlacements: false });
    }
    await window.Vue.nextTick();
    for (const [mode, { ui }] of Object.entries(drawings)) {
      config.applyUiStateData(state.drawingOf(mode), ui, { restorePreviewNavigation: false });
    }
  }, { drawings: restorePoint.drawings, shownAtBoot });
  await settleLive(page);
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.captureBaseline('Initial state'));
  await settleLive(page);
  // What a case can leave in the page outside the app.
  await page.evaluate(({ testLists, storage, generateGlobals, openSections }) => {
    generateGlobals.forEach((key) => { delete window[key]; });
    Object.entries(testLists).forEach(([key, items]) => {
      if (Array.isArray(window[key])) window[key].splice(0, window[key].length, ...items);
    });
    [localStorage, sessionStorage].forEach((area, index) => {
      area.clear();
      Object.entries(storage[index]).forEach(([key, value]) => area.setItem(key, value));
    });
    document.querySelectorAll('details[open]').forEach((section) => {
      const summary = section.querySelector(':scope > summary');
      if (!openSections.includes(summary?.textContent.trim())) section.open = false;
    });
  }, { ...restorePoint, generateGlobals: GENERATE_GLOBALS, openSections: facts.page.openSections });
  await page.removeAllListeners('dialog');
  await page.unrouteAll({ behavior: 'ignoreErrors' });
  const diagnostics = installPageErrorCollection(page);
  diagnostics.pageErrors.length = 0;
  diagnostics.consoleErrors.length = 0;
  const brief = (value) => {
    const text = JSON.stringify(value) ?? String(value);
    return text.length > 120 ? `${text.slice(0, 117)}...` : text;
  };
  return diffUserOwnedState(facts, await readWarmFacts(page)).map(({ path, before, after }) => (
    `${path}: ${brief(before)} -> ${brief(after)}`
  ));
};

// The reasons the page cannot serve the next case; none when it can.
const restoreWarmPage = async (page, contextOptions) => {
  const warm = warmPages.get(page);
  if (page.isClosed()) return ['the page closed'];
  if (!warm?.baseline) return ['the case did not open the app'];
  if (warm.contextOptions !== contextOptions) return ['the context options changed'];
  let timer = null;
  try {
    return await Promise.race([
      returnToFreshBoot(page, warm.baseline),
      new Promise((resolve) => {
        timer = setTimeout(() => resolve([`the restore took over ${WARM_RESTORE_TIMEOUT_MS / 1000} s`]), WARM_RESTORE_TIMEOUT_MS);
      })
    ]);
  } catch (error) {
    return [`the restore failed: ${String(error?.message || error).split('\n')[0]}`];
  } finally {
    clearTimeout(timer);
  }
};

const warmAppTest = baseTest.extend({
  warmAppSlot: [async ({}, use) => {
    const slot = { page: null };
    await use(slot);
    if (slot.page && !slot.page.isClosed()) await slot.page.context().close();
  }, { scope: 'worker' }],
  // The project's context options, as the built-in `context` fixture combines them.
  page: async ({
    warmAppSlot, browser, baseURL, viewport, userAgent, deviceScaleFactor, isMobile, hasTouch
  }, use, testInfo) => {
    const options = { baseURL, viewport, userAgent, deviceScaleFactor, isMobile, hasTouch };
    const contextOptions = JSON.stringify(options);
    let page = warmAppSlot.page;
    const reasons = page ? await restoreWarmPage(page, contextOptions) : ['the worker had no app page yet'];
    if (reasons.length) {
      if (page && !page.isClosed()) await page.context().close().catch(() => {});
      page = await (await browser.newContext(options)).newPage();
      const savedSessions = [];
      page.on('download', (file) => {
        if (SESSION_FILE.test(file.suggestedFilename())) savedSessions.push(file.suggestedFilename());
      });
      warmPages.set(page, { contextOptions, baseline: null, ready: false, savedSessions });
    } else {
      warmPages.get(page).ready = true;
    }
    warmAppSlot.page = page;
    testInfo.annotations.push({
      type: WARM_APP_PAGE,
      description: reasons.length ? `fresh: ${reasons.join('; ')}` : 'reused'
    });
    await use(page);
  }
});

const assertWorkerReuseAcrossHelperAndRender = async (page) => {
  const diagnostics = await getLifecycleDiagnostics(page);
  const activity = diagnostics.worker;
  const rendered = compactJson(diagnostics);
  expect(activity.constructions, rendered).toBe(1);
  expect(activity.initializations, rendered).toBe(1);
  expect(activity.helpers, rendered).toBeGreaterThan(0);
  expect(activity.runs, rendered).toBeGreaterThan(0);
  expect(activity.settledInitializations, rendered).toBe(1);
  expect(activity.settledHelpers, rendered).toBeGreaterThan(0);
  expect(activity.settledRuns, rendered).toBeGreaterThan(0);
  expect(activity.instances, rendered).toHaveLength(1);
  expect(activity.instances[0].terminated, rendered).toBe(false);
  const events = activity.instances[0].events || [];
  const firstHelper = events.findIndex((event) => event.startsWith('helper:request:'));
  const firstRun = events.indexOf('run:request');
  expect(firstHelper, rendered).toBeGreaterThan(-1);
  expect(firstRun, rendered).toBeGreaterThan(firstHelper);
  return activity;
};

const assertSingleWorkerRun = async (page) => {
  const diagnostics = await getLifecycleDiagnostics(page);
  const activity = diagnostics.worker;
  const rendered = compactJson(diagnostics);
  expect(activity.constructions, rendered).toBe(1);
  expect(activity.initializations, rendered).toBe(1);
  expect(activity.runs, rendered).toBe(1);
  expect(activity.settledInitializations, rendered).toBe(1);
  expect(activity.settledRuns, rendered).toBe(1);
  return activity;
};

// What a Load that ended or asked shows: no dialog lines, no error.
const expectNoRepairsOrError = async (page) => {
  const prompt = page.locator('[data-choice-prompt]');
  expect(await prompt.locator('li').allTextContents(), 'the Load listed repairs').toEqual([]);
  await expect(prompt).toHaveCount(0);
  expect(await page.evaluate(() => {
    const error = window.__GBDRAW_APP__.errorLog;
    return error ? { code: error.code, summary: error.summary } : null;
  }), 'the Load failed').toBeNull();
};

// A Load of a file that loses nothing ends without the Load repair dialog
// (BA-D25). Called right after the file is chosen, it waits until the new
// Load has started, so a notice of an earlier Load does not end the wait,
// then until the Load ends or the dialog opens. It fails with the
// dialog's lines when the dialog opened, and with the error when the Load
// failed.
const expectNoLoadRepairDialog = async (page, { timeout = 180_000 } = {}) => {
  await page.waitForFunction(() => window.__GBDRAW_APP__?.sessionImportPending, null, { timeout });
  await page.waitForFunction(() => document.querySelector('[data-choice-prompt]') || !window.__GBDRAW_APP__.sessionImportPending,
    null, { timeout });
  await expectNoRepairsOrError(page);
};

// A Session Load succeeded with no warning: the notice says "Session loaded"
// and lists nothing (each warning of a Load is one list item of the notice).
// Its callers may have waited for the Load's end or answered its dialog, so it
// does not wait for the start; a notice or an error counts once no Load runs.
const expectSessionLoadedWithoutWarnings = async (page, { timeout = 120_000 } = {}) => {
  await page.waitForFunction(() => document.querySelector('[data-choice-prompt]') || (!window.__GBDRAW_APP__?.sessionImportPending
    && (window.__GBDRAW_APP__?.errorLog || /Session loaded/.test(document.querySelector('[data-app-notice]')?.textContent || ''))),
  null, { timeout });
  await expectNoRepairsOrError(page);
  await expect(page.locator('[data-app-notice]')).toContainText('Session loaded', { timeout });
  await expect(page.locator('[data-app-notice] li')).toHaveCount(0);
};

// The app asks in its own dialog, which no browser `dialog` handler answers:
// the spec names the question it expects and the action it takes (the
// action's id; the dialog's keep button is not an action). The question of a
// Save comes after the Save's preparation, hence the long wait.
const answerChoicePrompt = async (page, title, action, { timeout = 120_000 } = {}) => {
  const prompt = page.locator('[data-choice-prompt]');
  await expect(prompt.getByRole('heading', { name: title, exact: true })).toBeVisible({ timeout });
  await prompt.locator(`[data-choice-prompt-action="${action}"]`).click();
};

// Save asks for a Session title in the app's dialog when the Session has none.
// The spec answers with `title` ('' saves an untitled file) and Enter.
const answerSessionTitle = async (page, title = '', { timeout = 120_000 } = {}) => {
  const dialog = page.locator('[data-choice-prompt]').getByRole('dialog', { name: 'Save Session', exact: true });
  const field = dialog.getByRole('textbox', { name: 'Session title', exact: true });
  await expect(field).toBeFocused({ timeout });
  await field.fill(title);
  await field.press('Enter');
  await expect(dialog).toHaveCount(0);
};

// BA-D25: a Load that cannot restore everything lists what it cannot in one
// dialog before it changes anything. The spec names the lines it expects and
// answers `load` or `cancel`.
const answerLoadRepairs = async (page, items, answer = 'load', { timeout = 120_000 } = {}) => {
  const prompt = page.locator('[data-choice-prompt]');
  await expect(prompt.getByRole('heading', { name: 'Load this file?', exact: true })).toBeVisible({ timeout });
  await expect(prompt.locator('li')).toHaveText(items);
  if (answer === 'load') await prompt.locator('[data-choice-prompt-action="load"]').click();
  else await prompt.getByRole('button', { name: 'Cancel', exact: true }).click();
  await expect(prompt).toHaveCount(0);
};

// The local static server can reset a connection under CI load (OV-713).
// Browser tests fetch a routed request through this helper, which retries once
// on a connection error and rethrows any other error unchanged.
const CONNECTION_ERROR = /\b(?:ECONNRESET|ECONNREFUSED|ECONNABORTED|EPIPE)\b|socket hang up/i;
const fetchRoute = async (route, options) => {
  try {
    return await route.fetch(options);
  } catch (error) {
    if (!CONNECTION_ERROR.test(String(error?.message || error))) throw error;
    return route.fetch(options);
  }
};

// Serves the app cross-origin isolated (threaded LOSAT needs
// SharedArrayBuffer): every 127.0.0.1 request gets COOP/COEP headers and other
// hosts are aborted. `stubs` maps a path to the fulfill options served instead.
// A served Pyodide snapshot passes through unchanged: the Worker fetches it
// same-origin, which COEP allows without headers, and fulfilling its 100 MB
// body from the test process closes the page (OV-1020).
const CROSS_ORIGIN_ISOLATION_HEADERS = {
  'Cross-Origin-Opener-Policy': 'same-origin',
  'Cross-Origin-Embedder-Policy': 'require-corp'
};
const routeCrossOriginIsolated = (context, stubs = {}) => context.route('**/*', async (route) => {
  const url = new URL(route.request().url());
  if (url.hostname !== '127.0.0.1') return route.abort();
  if (url.pathname.startsWith(PYODIDE_SNAPSHOT_TEST_PATH)) return route.fallback();
  const stub = stubs[url.pathname];
  if (stub) return route.fulfill({ ...stub, headers: CROSS_ORIGIN_ISOLATION_HEADERS });
  const response = await fetchRoute(route);
  return route.fulfill({ response, headers: { ...response.headers(), ...CROSS_ORIGIN_ISOLATION_HEADERS } });
});

module.exports = {
  answerChoicePrompt,
  answerSessionTitle,
  expectNoLoadRepairDialog,
  expectSessionLoadedWithoutWarnings,
  answerLoadRepairs,
  CURRENT_FEATURE_CATALOG_SCHEMA,
  CURRENT_REQUEST_SCHEMA,
  CURRENT_SESSION_VERSION,
  assertDiagramWorkerIdle,
  assertOperationHealth,
  readErrorSignature,
  assertSessionLoadLeftWorkerIdle,
  assertSingleWorkerRun,
  assertWorkerReuseAcrossHelperAndRender,
  diffUserOwnedState,
  evaluateWithRetainedPromise,
  expectNoSilentStateChange,
  fetchRoute,
  generateAndWaitForResult,
  getDiagramWorkerActivity,
  observeNonEditOperation,
  openLinearFileCards,
  openApp,
  reveal,
  routeCrossOriginIsolated,
  settleLive,
  snapshotUserOwnedState,
  waitForAppShell,
  WARM_APP_PAGE,
  warmAppTest
};
