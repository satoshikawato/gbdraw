const { expect } = require('@playwright/test');

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

const generateAndWaitForResult = async (
  page,
  { expectedStatus = 'ok', requireCommittedResult = expectedStatus === 'ok' } = {}
) => {
  const outcome = await page.evaluate(async () => {
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
  return outcome;
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
  assertDiagramWorkerIdle,
  assertSessionLoadLeftWorkerIdle,
  assertSingleWorkerRun,
  assertWorkerReuseAcrossHelperAndRender,
  generateAndWaitForResult,
  getDiagramWorkerActivity,
  openApp,
  waitForAppShell
};
