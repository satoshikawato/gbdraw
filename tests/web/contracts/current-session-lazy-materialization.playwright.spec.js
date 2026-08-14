const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { join, resolve } = require('node:path');
const { gunzipSync } = require('node:zlib');
const {
  getDiagramWorkerActivity,
  openApp
} = require('../helpers/app-lifecycle.cjs');

const repoRoot = resolve(process.env.GBDRAW_REPO || process.cwd());
const webSessionRoot = join(repoRoot, 'gbdraw', 'web', 'gallery', 'sessions');
const syntheticSourceName = 'HmmtDNA_basic_circular.gbdraw-session.json';
const syntheticSource = JSON.parse(readFileSync(
  join(webSessionRoot, syntheticSourceName),
  'utf8'
));
const vibrioSession = join(
  webSessionRoot,
  'vibrio-harveyi-group-collinear.gbdraw-session.json.gz'
);
const frozenV39Session = join(
  repoRoot,
  'tests',
  'fixtures',
  'sessions',
  'BGC0000708-BGC0000713.v39.gbdraw-session.json.gz'
);

const STRUCTURAL_METRICS = [
  'base64DecodeCount',
  'decodedByteCount',
  'fileConstructionCount',
  'blobConstructionCount',
  'resourceTextReadCount',
  'resourceByteReadCount',
  'sourceRecoveryCount',
  'workerConstructionCount',
  'workerInitializationCount'
];

const ZERO_PREVIEW_METRICS = Object.fromEntries(
  STRUCTURAL_METRICS.map((name) => [name, 0])
);

const installLazySessionProbe = async (page) => page.addInitScript((metricNames) => {
  const metricMap = () => Object.fromEntries(metricNames.map((name) => [name, 0]));
  const hookMetrics = metricMap();
  const nativeMetrics = metricMap();
  const details = [];
  const lifecycle = [];
  const constructedFiles = new WeakSet();
  const ignoredFiles = new WeakSet();

  const probe = {
    hookMetrics,
    nativeMetrics,
    details,
    lifecycle,
    historyLoaded: false,
    ignoreFile(file) {
      if (file && typeof file === 'object') ignoredFiles.add(file);
    },
    reset() {
      Object.keys(hookMetrics).forEach((name) => {
        hookMetrics[name] = 0;
      });
      metricNames.forEach((name) => {
        nativeMetrics[name] = 0;
      });
      details.length = 0;
      lifecycle.length = 0;
      this.historyLoaded = false;
    },
    snapshot() {
      return {
        structural: Object.fromEntries(metricNames.map((name) => [
          name,
          Math.max(Number(hookMetrics[name] || 0), Number(nativeMetrics[name] || 0))
        ])),
        hookMetrics: { ...hookMetrics },
        nativeMetrics: { ...nativeMetrics },
        details: details.map((detail) => ({ ...detail })),
        lifecycle: lifecycle.map((event) => ({ ...event })),
        historyLoaded: this.historyLoaded
      };
    }
  };
  window.__GBDRAW_LAZY_SESSION_PROBE__ = probe;
  window.__GBDRAW_TEST_HOOKS__ = {
    onStructuralMetric(metric) {
      const name = String(metric?.name || '');
      if (!Object.hasOwn(hookMetrics, name)) hookMetrics[name] = 0;
      hookMetrics[name] += Number(metric.value || 0);
      details.push({ ...metric, timestamp: performance.now() });
    },
    onSessionLifecycleEvent(event) {
      lifecycle.push({ ...event });
    }
  };

  const nativeAtob = window.atob;
  window.atob = function lazySessionTrackedAtob(value) {
    const decoded = nativeAtob.call(this, value);
    nativeMetrics.base64DecodeCount += 1;
    nativeMetrics.decodedByteCount += decoded.length;
    return decoded;
  };

  const NativeBlob = window.Blob;
  window.Blob = new Proxy(NativeBlob, {
    construct(target, args) {
      nativeMetrics.blobConstructionCount += 1;
      return Reflect.construct(target, args, target);
    }
  });
  const NativeFile = window.File;
  window.File = new Proxy(NativeFile, {
    construct(target, args) {
      nativeMetrics.fileConstructionCount += 1;
      const file = Reflect.construct(target, args, target);
      constructedFiles.add(file);
      return file;
    }
  });
  const nativeArrayBuffer = NativeBlob.prototype.arrayBuffer;
  NativeBlob.prototype.arrayBuffer = function lazySessionTrackedArrayBuffer(...args) {
    if (constructedFiles.has(this) && !ignoredFiles.has(this)) {
      nativeMetrics.resourceByteReadCount += 1;
    }
    return nativeArrayBuffer.apply(this, args);
  };
  const nativeText = NativeBlob.prototype.text;
  NativeBlob.prototype.text = function lazySessionTrackedText(...args) {
    if (constructedFiles.has(this) && !ignoredFiles.has(this)) {
      nativeMetrics.resourceTextReadCount += 1;
    }
    return nativeText.apply(this, args);
  };
}, STRUCTURAL_METRICS);

const armHistoryCompletion = (page) => page.evaluate(() => {
  const history = window.__GBDRAW_HISTORY__;
  const probe = window.__GBDRAW_LAZY_SESSION_PROBE__;
  const original = history.captureBaseline;
  probe.historyLoaded = false;
  history.captureBaseline = async (label, ...args) => {
    const result = await original(label, ...args);
    if (label === 'Loaded session') {
      probe.historyLoaded = true;
      history.captureBaseline = original;
    }
    return result;
  };
});

const probeSnapshot = (page) => page.evaluate(() => ({
  ...window.__GBDRAW_LAZY_SESSION_PROBE__.snapshot(),
  savedPreviewVisible: Boolean(document.querySelector('.shadow-xl.origin-top > svg')),
  selectedResultIndex: window.__GBDRAW_APP__?.selectedResultIndex ?? null,
  resultCount: window.__GBDRAW_APP__?.results?.length || 0
}));

const phaseDuration = (events, startName, endName) => {
  const start = events.find(({ name }) => name === startName);
  const end = events.find(({ name }) => name === endName);
  return start && end ? end.timestamp - start.timestamp : null;
};

const PREFLIGHT_SUBPHASES = Object.freeze({
  sessionAuthorityValidationMs: 'session-authority-validation',
  featureCatalogValidationMs: 'feature-catalog-validation',
  editorStateNormalizationMs: 'editor-state-normalization',
  resourceTableAdoptionMs: 'resource-table-adoption',
  losatArtifactValidationMs: 'losat-artifact-validation',
  canonicalRequestProjectionMs: 'canonical-request-projection',
  currentDraftValidationMs: 'current-draft-validation',
  artifactProjectionMs: 'artifact-projection'
});

const loadSyntheticSession = (page, variant = 'normal') => page.evaluate(
  async ({ filename, requestedVariant }) => {
    const response = await fetch(`/gbdraw/web/gallery/sessions/${filename}`);
    if (!response.ok) throw new Error(`Could not load ${filename}: ${response.status}`);
    const session = await response.json();
    session.title = `lazy-${requestedVariant}`;
    session.resources['unused-lazy-contract'] = {
      kind: 'web-file',
      name: 'unused-lazy-contract.txt',
      type: 'text/plain',
      size: 7,
      lastModified: 0,
      encoding: 'base64',
      data: btoa('unused\n')
    };
    const recordResourceId = session.renderRequest.records[0].source.resourceId;
    if (requestedVariant === 'invalid-size') {
      session.resources['unused-lazy-contract'].size = -1;
    } else if (requestedVariant === 'invalid-base64') {
      session.resources[recordResourceId].data = '%%%';
    }
    const file = new File(
      [JSON.stringify(session)],
      `lazy-${requestedVariant}.gbdraw-session.json`,
      { type: 'application/json' }
    );
    const probe = window.__GBDRAW_LAZY_SESSION_PROBE__;
    probe.ignoreFile(file);
    probe.reset();
    const result = await window.__GBDRAW_APP__.importSession({
      target: { files: [file], value: 'selected' }
    });
    return {
      status: result?.status,
      degradedRecovery: Boolean(result?.degradedRecovery),
      message: String(result?.error?.message || '')
    };
  },
  { filename: syntheticSourceName, requestedVariant: variant }
);

const openInstrumentedApp = async (page) => {
  page.on('dialog', (dialog) => dialog.accept());
  await installLazySessionProbe(page);
  await openApp(page);
};

const lifecycleTiming = (snapshot) => {
  const selected = snapshot.lifecycle.find(({ name }) => name === 'sessionSelection');
  const preview = snapshot.lifecycle.find(({ name }) => name === 'firstCommittedPreview');
  const ready = snapshot.lifecycle.find(({ name }) => name === 'interactiveReady');
  return {
    currentSessionPreflightMs: phaseDuration(
      snapshot.lifecycle,
      'current-session-preflight-start',
      'current-session-preflight-end'
    ),
    preflightSubphases: Object.fromEntries(
      Object.entries(PREFLIGHT_SUBPHASES).map(([metric, event]) => [
        metric,
        phaseDuration(snapshot.lifecycle, `${event}-start`, `${event}-end`)
      ])
    ),
    firstPreviewMs: preview && selected ? preview.timestamp - selected.timestamp : null,
    interactiveReadyMs: ready && selected ? ready.timestamp - selected.timestamp : null,
    ready
  };
};

test.describe.configure({ mode: 'serial' });

test('synthetic current session restores and exports without materializing resources', async ({
  page
}) => {
  test.setTimeout(180_000);
  await openInstrumentedApp(page);

  const imported = await loadSyntheticSession(page);
  expect(imported).toEqual({ status: 'ok', degradedRecovery: false, message: '' });
  const preview = await probeSnapshot(page);
  expect(preview.savedPreviewVisible).toBe(true);
  expect(preview.resultCount).toBe(1);
  expect(preview.structural).toEqual(ZERO_PREVIEW_METRICS);

  await page.evaluate(() => {
    window.__GBDRAW_APP__.sessionTitle = 'lazy-unchanged-export';
  });
  const downloadPromise = page.waitForEvent('download', { timeout: 120_000 });
  const saved = await page.evaluate(() => window.__GBDRAW_APP__.saveSessionWithTitle());
  expect(saved.status).toBe('saved');
  const download = await downloadPromise;
  const exported = JSON.parse(gunzipSync(readFileSync(await download.path())).toString('utf8'));
  const afterExport = await probeSnapshot(page);

  expect(afterExport.structural.base64DecodeCount).toBe(0);
  expect(afterExport.hookMetrics.base64EncodeCount || 0).toBe(0);
  for (const resourceId of [
    'record-1-genbank',
    'colors-default-colors'
  ]) {
    expect(exported.resources[resourceId]).toEqual(syntheticSource.resources[resourceId]);
  }
  expect(exported.resources['unused-lazy-contract']).toEqual({
    kind: 'web-file',
    name: 'unused-lazy-contract.txt',
    type: 'text/plain',
    size: 7,
    lastModified: 0,
    encoding: 'base64',
    data: Buffer.from('unused\n').toString('base64')
  });
});

test('real Vibrio preview is lazy and leaves the Worker idle', async ({ page }, testInfo) => {
  test.setTimeout(180_000);
  await openInstrumentedApp(page);
  await page.evaluate(() => window.__GBDRAW_LAZY_SESSION_PROBE__.reset());
  await armHistoryCompletion(page);

  const input = page.locator(
    'input[type="file"][accept*="application/json"][accept*="application/gzip"]'
  );
  await input.setInputFiles(vibrioSession);
  await page.waitForFunction(
    () => window.__GBDRAW_LAZY_SESSION_PROBE__?.historyLoaded === true,
    null,
    { timeout: 180_000 }
  );
  await page.evaluate(() => new Promise((resolve) => (
    requestAnimationFrame(() => requestAnimationFrame(resolve))
  )));

  const snapshot = await probeSnapshot(page);
  const timing = lifecycleTiming(snapshot);
  const preflightStructural = {
    proteinManifestFullValidationCount: Number(
      snapshot.hookMetrics.currentSessionPreflightProteinManifestValidationCount || 0
    ),
    proteinRawTextValidationCount: Number(
      snapshot.hookMetrics.currentSessionPreflightProteinRawTextValidationCount || 0
    )
  };
  expect(snapshot.savedPreviewVisible).toBe(true);
  expect(snapshot.resultCount).toBe(1);
  expect(snapshot.structural).toEqual(ZERO_PREVIEW_METRICS);
  expect(timing.ready).toMatchObject({
    status: 'success',
    degradedRecovery: false
  });
  expect(timing.firstPreviewMs).toBeGreaterThan(0);
  expect(timing.interactiveReadyMs).toBeGreaterThanOrEqual(timing.firstPreviewMs);
  expect(timing.currentSessionPreflightMs).toBeLessThan(30_000);
  expect(preflightStructural.proteinManifestFullValidationCount).toBe(1);
  expect(preflightStructural.proteinRawTextValidationCount).toBeGreaterThan(0);

  await testInfo.attach('vibrio-lazy-session-metrics.json', {
    body: Buffer.from(JSON.stringify({
      structural: snapshot.structural,
      preflightStructural,
      timing
    }, null, 2)),
    contentType: 'application/json'
  });
  console.log(`Vibrio lazy-session evidence: ${JSON.stringify({
    structural: snapshot.structural,
    preflightStructural,
    timing
  })}`);
});

test('Generate materializes only required resources and reuses one Worker', async ({ page }) => {
  test.setTimeout(300_000);
  await openInstrumentedApp(page);
  expect(await loadSyntheticSession(page)).toMatchObject({
    status: 'ok',
    degradedRecovery: false
  });
  const preview = await probeSnapshot(page);
  expect(preview.structural).toEqual(ZERO_PREVIEW_METRICS);

  const generated = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const first = await app.runAnalysis();
    const {
      DIAGRAM_HELPER_OPERATIONS,
      runDiagramHelperOperation
    } = await import('/gbdraw/web/js/services/diagram-generation.js');
    const helper = await runDiagramHelperOperation(
      DIAGRAM_HELPER_OPERATIONS.MEASURE_LEGEND_TEXT,
      { caption: 'lazy worker reuse', fontFamily: 'Arial', fontSize: 14 }
    );
    const second = await app.runAnalysis();
    return {
      first,
      second,
      helperWidth: Number(helper?.result?.width || 0),
      errorSummary: String(app.errorLog?.summary || '')
    };
  });
  expect(generated).toMatchObject({
    first: { status: 'ok' },
    second: { status: 'ok' },
    errorSummary: ''
  });
  expect(generated.helperWidth).toBeGreaterThan(0);

  const snapshot = await probeSnapshot(page);
  const materializedIds = new Set(
    snapshot.details
      .filter(({ name }) => name === 'resourceByteReadCount')
      .map(({ resourceId }) => resourceId)
  );
  expect(materializedIds.has('record-1-genbank')).toBe(true);
  expect(materializedIds.has('unused-lazy-contract')).toBe(false);
  expect([...materializedIds].every((resourceId) => (
    ['record-1-genbank', 'colors-default-colors'].includes(resourceId)
  ))).toBe(true);
  expect(snapshot.structural.base64DecodeCount).toBe(materializedIds.size);

  const previewEvent = snapshot.lifecycle.find(
    ({ name }) => name === 'firstCommittedPreview'
  );
  const workerInitEvent = snapshot.details.find(
    ({ name }) => name === 'workerInitializationCount'
  );
  expect(previewEvent?.timestamp).toBeLessThan(workerInitEvent?.timestamp);
  expect(snapshot.structural.workerConstructionCount).toBe(1);
  expect(snapshot.structural.workerInitializationCount).toBe(1);

  const worker = await getDiagramWorkerActivity(page);
  expect(worker.constructions).toBe(1);
  expect(worker.initializations).toBe(1);
  expect(worker.helpers).toBeGreaterThanOrEqual(1);
  expect(worker.runs).toBe(2);
  expect(worker.instances).toHaveLength(1);
  expect(worker.instances[0].terminated).toBe(false);
});

test('preflight and lazy-access failures preserve the committed preview', async ({ page }) => {
  test.setTimeout(180_000);
  await openInstrumentedApp(page);
  expect(await loadSyntheticSession(page)).toMatchObject({ status: 'ok' });
  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    window.__GBDRAW_LAZY_ROLLBACK_BASELINE__ = {
      file: app.files.c_gb,
      content: app.results[app.selectedResultIndex].content,
      selectedResultIndex: app.selectedResultIndex,
      undoCount: window.__GBDRAW_HISTORY__.getUndoCount()
    };
  });

  const rejected = await loadSyntheticSession(page, 'invalid-size');
  expect(rejected.status).toBe('error');
  expect(rejected.message).toMatch(/unused-lazy-contract.*invalid declared byte size/);
  expect(await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const baseline = window.__GBDRAW_LAZY_ROLLBACK_BASELINE__;
    return {
      sameFile: app.files.c_gb === baseline.file,
      sameContent: app.results[app.selectedResultIndex].content === baseline.content,
      selectedResultIndex: app.selectedResultIndex,
      undoCount: window.__GBDRAW_HISTORY__.getUndoCount(),
      previewVisible: Boolean(document.querySelector('.shadow-xl.origin-top > svg'))
    };
  })).toEqual({
    sameFile: true,
    sameContent: true,
    selectedResultIndex: 0,
    undoCount: 0,
    previewVisible: true
  });

  const corrupt = await loadSyntheticSession(page, 'invalid-base64');
  expect(corrupt).toEqual({ status: 'ok', degradedRecovery: false, message: '' });
  const beforeAccess = await probeSnapshot(page);
  expect(beforeAccess.structural).toEqual(ZERO_PREVIEW_METRICS);
  const access = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const content = app.results[app.selectedResultIndex].content;
    const { readFileText } = await import('/gbdraw/web/js/services/file-content-cache.js');
    let first = '';
    let second = '';
    try {
      await readFileText(app.files.c_gb);
    } catch (error) {
      first = String(error?.message || error);
    }
    try {
      await readFileText(app.files.c_gb);
    } catch (error) {
      second = String(error?.message || error);
    }
    return {
      first,
      second,
      sameContent: app.results[app.selectedResultIndex].content === content,
      previewVisible: Boolean(document.querySelector('.shadow-xl.origin-top > svg'))
    };
  });
  expect(access.first).toMatch(
    /record-1-genbank \(record-1\.gbk\) contains invalid encoded data/
  );
  expect(access.second).toBe(access.first);
  expect(access.sameContent).toBe(true);
  expect(access.previewVisible).toBe(true);
  const afterAccess = await probeSnapshot(page);
  expect(afterAccess.structural.base64DecodeCount).toBe(1);
  expect(afterAccess.structural.resourceByteReadCount).toBe(1);
  expect(afterAccess.structural.workerConstructionCount).toBe(0);
  expect(afterAccess.structural.workerInitializationCount).toBe(0);
});

test('a frozen v39 session round-trips through the legacy migration path', async ({ page }) => {
  test.setTimeout(240_000);
  await openInstrumentedApp(page);
  await armHistoryCompletion(page);
  const input = page.locator(
    'input[type="file"][accept*="application/json"][accept*="application/gzip"]'
  );
  await input.setInputFiles(frozenV39Session);
  await page.waitForFunction(
    () => window.__GBDRAW_LAZY_SESSION_PROBE__?.historyLoaded === true,
    null,
    { timeout: 180_000 }
  );
  expect((await probeSnapshot(page)).savedPreviewVisible).toBe(true);

  await page.evaluate(() => {
    window.__GBDRAW_APP__.sessionTitle = 'legacy-lazy-round-trip';
  });
  const regenerated = await page.evaluate(async () => ({
    result: await window.__GBDRAW_APP__.runAnalysis(),
    errorLog: window.__GBDRAW_APP__.errorLog
  }));
  expect(regenerated.result, JSON.stringify(regenerated.errorLog)).toEqual({ status: 'ok' });
  const downloadPromise = page.waitForEvent('download', { timeout: 120_000 });
  const saveOutcome = await page.evaluate(async () => {
    const result = await window.__GBDRAW_APP__.saveSessionWithTitle();
    return {
      result,
      errorLog: window.__GBDRAW_APP__.errorLog
    };
  });
  expect(saveOutcome.result.status, JSON.stringify(saveOutcome.errorLog)).toBe('saved');
  const roundTripPath = await (await downloadPromise).path();
  const roundTrip = JSON.parse(gunzipSync(readFileSync(roundTripPath)).toString('utf8'));
  expect(roundTrip.version).toBe(40);
  expect(roundTrip.results).toHaveLength(1);
  expect(roundTrip.editorState.featureCatalog).toBeTruthy();

  await page.reload({ waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => Boolean(window.__GBDRAW_APP__), null, {
    timeout: 180_000
  });
  await page.evaluate(() => window.__GBDRAW_LAZY_SESSION_PROBE__.reset());
  await armHistoryCompletion(page);
  await page.locator(
    'input[type="file"][accept*="application/json"][accept*="application/gzip"]'
  ).setInputFiles(roundTripPath);
  await page.waitForFunction(
    () => window.__GBDRAW_LAZY_SESSION_PROBE__?.historyLoaded === true,
    null,
    { timeout: 180_000 }
  );
  const restored = await probeSnapshot(page);
  expect(restored.savedPreviewVisible).toBe(true);
  expect(restored.resultCount).toBe(1);
  expect(restored.structural.workerConstructionCount).toBe(0);
  expect(restored.structural.workerInitializationCount).toBe(0);
});
