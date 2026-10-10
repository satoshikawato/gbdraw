const { readFileSync } = require('node:fs');
const { join } = require('node:path');
const { expect } = require('@playwright/test');

// The Web writer's current canonical request schema, read from its owner.
const CURRENT_REQUEST_SCHEMA = Number(readFileSync(
  join(__dirname, '..', '..', '..', 'gbdraw', 'web', 'js', 'services', 'session-request.js'), 'utf8'
).match(/^export const CANONICAL_REQUEST_SCHEMA = (\d+);$/m)[1]);
// The current Session writer version and feature catalog schema, from their owners.
const CURRENT_SESSION_VERSION = Number(readFileSync(
  join(__dirname, '..', '..', '..', 'gbdraw', 'session_io.py'), 'utf8'
).match(/^CURRENT_SESSION_VERSION = (\d+)$/m)[1]);
const CURRENT_FEATURE_CATALOG_SCHEMA = Number(readFileSync(
  join(__dirname, '..', '..', '..', 'gbdraw', 'web', 'js', 'services', 'feature-catalog.js'), 'utf8'
).match(/^export const FEATURE_CATALOG_SCHEMA = (\d+);$/m)[1]);

const DEFAULT_APP_TIMEOUT_MS = 180_000;
const pageDiagnostics = new WeakMap();
const workerTrackingPages = new WeakSet();

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
    path = '/gbdraw/web/index.html',
    waitForPalette = true,
    timeout = DEFAULT_APP_TIMEOUT_MS,
    checkErrors = true
  } = {}
) => {
  installPageErrorCollection(page);
  await installDiagramWorkerTracking(page);
  await page.goto(path, { waitUntil: 'domcontentloaded' });
  await waitForAppShell(page, { waitForPalette, timeout });
  return assertAppShellReady(page, { waitForPalette, checkErrors });
};

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
// navigation are excluded.
const snapshotUserOwnedState = (page) => page.evaluate(async () => {
  const { state } = await import('/gbdraw/web/js/state.js');
  const config = await import('/gbdraw/web/js/services/config.js');
  const drawing = state.activeDrawing();
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
  ui.layoutPreferences = ui.layoutPreferences?.[state.mode.value];
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
});

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

module.exports = {
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
  generateAndWaitForResult,
  getDiagramWorkerActivity,
  observeNonEditOperation,
  openApp,
  reveal,
  snapshotUserOwnedState,
  waitForAppShell
};
