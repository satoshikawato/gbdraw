const { test, expect } = require('@playwright/test');
const { spawnSync } = require('node:child_process');
const { readFileSync, writeFileSync } = require('node:fs');
const { join, resolve } = require('node:path');
const { gunzipSync } = require('node:zlib');
const {
  CURRENT_SESSION_VERSION,
  evaluateWithRetainedPromise,
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
const frozenV39Session = join(
  repoRoot,
  'tests',
  'fixtures',
  'sessions',
  'BGC0000708-BGC0000713.v39.gbdraw-session.json.gz'
);
const neutralConservationSession = join(
  repoRoot,
  'tests',
  'fixtures',
  'sessions',
  'synthetic_conservation.gbdraw-session.json.gz'
);
const neutralConservationSessionText = gunzipSync(
  readFileSync(neutralConservationSession)
).toString('utf8');

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

const ZERO_ARTIFACT_HISTORY_BASELINE = Object.freeze({
  artifactCheckpointBuilds: 0,
  artifactCheckpointSignatureComputations: 0,
  artifactSvgBytes: 0,
  intentBuilds: 1,
  intentSignatureComputations: 1,
  undoCount: 0,
  redoCount: 0,
  currentCheckpointAbsent: true
});

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
    historyBaseline: null,
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
      this.historyBaseline = null;
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
        historyLoaded: this.historyLoaded,
        historyBaseline: this.historyBaseline ? { ...this.historyBaseline } : null
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
  if (!history?.initializeIntentBaseline) {
    throw new Error('The lightweight session-import History boundary is unavailable.');
  }
  const original = history.initializeIntentBaseline;
  probe.historyLoaded = false;
  probe.historyBaseline = null;
  history.initializeIntentBaseline = async (label, ...args) => {
    const before = history.getDiagnostics();
    try {
      const result = await original(label, ...args);
      if (label === 'Loaded session') {
        const after = history.getDiagnostics();
        probe.historyBaseline = {
          artifactCheckpointBuilds:
            after.artifactCheckpointBuilds - before.artifactCheckpointBuilds,
          artifactCheckpointSignatureComputations:
            after.artifactCheckpointSignatureComputations
              - before.artifactCheckpointSignatureComputations,
          artifactSvgBytes: after.historySvgBytes - before.historySvgBytes,
          intentBuilds: after.intentBuilds - before.intentBuilds,
          intentSignatureComputations:
            after.intentSignatureComputations - before.intentSignatureComputations,
          undoCount: history.getUndoCount(),
          redoCount: history.getRedoCount(),
          currentCheckpointAbsent: history.getCurrentCheckpoint() === null
        };
        probe.historyLoaded = true;
        history.initializeIntentBaseline = original;
      }
      return result;
    } catch (error) {
      history.initializeIntentBaseline = original;
      throw error;
    }
  };
});

const probeSnapshot = (page) => page.evaluate(() => ({
  ...window.__GBDRAW_LAZY_SESSION_PROBE__.snapshot(),
  savedPreviewVisible: Boolean(document.querySelector('.shadow-xl.origin-top > svg')),
  selectedResultIndex: window.__GBDRAW_APP__?.selectedResultIndex ?? null,
  resultCount: window.__GBDRAW_APP__?.results?.length || 0
}));

const loadSyntheticSession = (page, variant = 'normal') => evaluateWithRetainedPromise(page,
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
      message: String(result?.error?.summary || '')
    };
  },
  { filename: syntheticSourceName, requestedVariant: variant }
);

const openInstrumentedApp = async (page) => {
  page.on('dialog', (dialog) => dialog.accept());
  await installLazySessionProbe(page);
  await openApp(page);
};

test('synthetic current session restores and exports without materializing resources', async ({
  page
}) => {
  test.setTimeout(180_000);
  await openInstrumentedApp(page);
  await armHistoryCompletion(page);

  const imported = await loadSyntheticSession(page);
  expect(imported).toEqual({ status: 'ok', degradedRecovery: false, message: '' });
  const preview = await probeSnapshot(page);
  expect(preview.savedPreviewVisible).toBe(true);
  expect(preview.resultCount).toBe(1);
  expect(preview.structural).toEqual(ZERO_PREVIEW_METRICS);
  expect(preview.historyBaseline).toEqual(ZERO_ARTIFACT_HISTORY_BASELINE);

  const intentHistory = await evaluateWithRetainedPromise(page, async () => {
    const app = window.__GBDRAW_APP__;
    const history = window.__GBDRAW_HISTORY__;
    const loadedValue = Boolean(app.form.show_scale);
    const loadedContent = app.results[app.selectedResultIndex]?.content || '';
    await history.runUndoable('Change coordinate scale', () => {
      app.form.show_scale = !loadedValue;
    });
    const editedValue = Boolean(app.form.show_scale);
    const undoResult = await history.undo();
    const undoValue = Boolean(app.form.show_scale);
    const redoResult = await history.redo();
    return {
      loadedValue,
      editedValue,
      undoResult,
      undoValue,
      redoResult,
      redoValue: Boolean(app.form.show_scale),
      undoCount: history.getUndoCount(),
      redoCount: history.getRedoCount(),
      previewUnchanged: app.results[app.selectedResultIndex]?.content === loadedContent
    };
  });
  expect(intentHistory).toMatchObject({
    editedValue: !intentHistory.loadedValue,
    undoResult: true,
    undoValue: intentHistory.loadedValue,
    redoResult: true,
    redoValue: !intentHistory.loadedValue,
    undoCount: 1,
    redoCount: 0,
    previewUnchanged: true
  });

  await page.evaluate(() => {
    window.__GBDRAW_APP__.sessionTitle = 'lazy-unchanged-export';
  });
  const downloadPromise = page.waitForEvent('download', { timeout: 120_000 });
  const saved = await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.saveSessionWithTitle());
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
  // Save writes only the resources a request or a binding names, as the
  // Python writer does (E1, REVIEW-1 m4); the unused one is dropped unread.
  expect(exported.resources).not.toHaveProperty('unused-lazy-contract');
});

test('neutral cached conservation session regenerates three ordered rings offline', async ({
  browser
}, testInfo) => {
  test.setTimeout(180_000);
  const context = await browser.newContext();
  const page = await context.newPage();
  const externalRequests = [];
  page.on('request', (request) => {
    const url = new URL(request.url());
    if (
      ['http:', 'https:'].includes(url.protocol)
      && !['127.0.0.1', 'localhost'].includes(url.hostname)
    ) externalRequests.push(request.url());
  });

  try {
    await page.addInitScript(() => {
      window.__GBDRAW_NEUTRAL_CONSERVATION_PROBE__ = { losatCalls: 0 };
      window.__GBDRAW_LOSAT_EXECUTOR__ = async () => {
        window.__GBDRAW_NEUTRAL_CONSERVATION_PROBE__.losatCalls += 1;
        throw new Error('Neutral cached replay must not execute LOSAT.');
      };
    });
    await openInstrumentedApp(page);

    const imported = await evaluateWithRetainedPromise(page, async (sessionText) => {
      const file = new File(
        [sessionText],
        'synthetic_conservation.gbdraw-session.json',
        { type: 'application/json', lastModified: 0 }
      );
      const probe = window.__GBDRAW_LAZY_SESSION_PROBE__;
      probe.ignoreFile(file);
      probe.reset();
      const result = await window.__GBDRAW_APP__.importSession({
        target: { files: [file], value: '' }
      });
      const { state } = await import('/gbdraw/web/js/state.js');
      const selected = window.__GBDRAW_APP__.results[
        window.__GBDRAW_APP__.selectedResultIndex
      ];
      const lazyFields = ['text', 'arrayBuffer', 'data', 'resourceId'];
      const resources = [
        state.files.c_gb,
        ...state.files.c_conservation_fastas,
        ...state.files.c_conservation_blasts
      ];
      return {
        status: result?.status,
        content: String(selected?.content || ''),
        referenceName: state.files.c_gb?.name,
        comparisonNames: state.files.c_conservation_fastas.map((fileValue) => fileValue.name),
        cachedResultNames: state.files.c_conservation_blasts.map((fileValue) => fileValue.name),
        cachedResultSource: state.files.c_conservation_blasts_source,
        cacheEntries: state.losatCache.value.size,
        metadataOnly: resources.map((fileValue) => ({
          frozen: Object.isFrozen(fileValue),
          ownFields: lazyFields.filter((field) => Object.hasOwn(fileValue, field))
        })),
        metrics: probe.snapshot()
      };
    }, neutralConservationSessionText);

    expect(imported.status).toBe('ok');
    expect(imported.referenceName).toBe('reference-a.gb');
    expect(imported.comparisonNames).toEqual([
      'comparison-b.fasta',
      'comparison-c.fasta',
      'comparison-d.fasta'
    ]);
    expect(imported.cachedResultNames).toEqual([
      'comparison-b.circular_conservation.losatn.tsv',
      'comparison-c.circular_conservation.losatn.tsv',
      'comparison-d.circular_conservation.losatn.tsv'
    ]);
    expect(imported.cachedResultSource).toBe('losat-cache');
    expect(imported.cacheEntries).toBe(3);
    expect(imported.metadataOnly).toEqual(
      Array.from({ length: 7 }, () => ({ frozen: true, ownFields: [] }))
    );
    expect(imported.metrics.structural).toEqual(ZERO_PREVIEW_METRICS);
    expect(await getDiagramWorkerActivity(page)).toMatchObject({
      constructions: 0,
      initializations: 0,
      helpers: 0,
      runs: 0
    });

    const generated = await evaluateWithRetainedPromise(page, async () => {
      const app = window.__GBDRAW_APP__;
      const result = await app.runAnalysis();
      const selected = app.results[app.selectedResultIndex];
      const content = String(selected?.content || '');
      const documentRoot = new DOMParser().parseFromString(content, 'image/svg+xml');
      const slots = [...documentRoot.querySelectorAll(
        '[data-gbdraw-slot-renderer="sequence_conservation"]'
      )];
      return {
        result,
        errorLog: app.errorLog,
        content,
        slots: slots.map((slot) => ({
          id: slot.getAttribute('data-gbdraw-slot-id'),
          sourceIndex: Number(slot.getAttribute('data-source-index')),
          label: slot.getAttribute('data-track-label'),
          color: slot.getAttribute('data-track-color')
        })),
        renderedMatchCount: documentRoot.querySelectorAll('[data-gbdraw-match-id]').length,
        metrics: window.__GBDRAW_LAZY_SESSION_PROBE__.snapshot(),
        losatCalls: window.__GBDRAW_NEUTRAL_CONSERVATION_PROBE__.losatCalls
      };
    });

    expect(generated.result).toEqual({ status: 'ok' });
    expect(generated.errorLog).toBeNull();
    expect(generated.slots).toEqual([
      { id: 'conservation_1', sourceIndex: 0, label: 'comparison-b', color: '#4e79a7' },
      { id: 'conservation_2', sourceIndex: 1, label: 'comparison-c', color: '#e15759' },
      { id: 'conservation_3', sourceIndex: 2, label: 'comparison-d', color: '#59a14f' }
    ]);
    expect(generated.renderedMatchCount).toBeGreaterThan(0);
    expect(generated.losatCalls).toBe(0);
    expect(externalRequests).toEqual([]);
    expect(
      generated.metrics.details
        .filter(({ name }) => name === 'resourceTextReadCount')
        .map(({ resourceId }) => resourceId)
    ).toEqual(['resource-0001']);
    // Ring files reach the one Python sequence reader as bytes (D12).
    expect(
      generated.metrics.details
        .filter(({ name, resourceId }) => (
          name === 'resourceByteReadCount' && resourceId.startsWith('conservation-losat-fasta-files-')
        ))
        .map(({ resourceId }) => resourceId)
    ).toEqual([
      'conservation-losat-fasta-files-1',
      'conservation-losat-fasta-files-2',
      'conservation-losat-fasta-files-3'
    ]);

    const initialPath = testInfo.outputPath('synthetic-conservation-loaded.svg');
    const generatedPath = testInfo.outputPath('synthetic-conservation-regenerated.svg');
    writeFileSync(initialPath, imported.content, 'utf8');
    writeFileSync(generatedPath, generated.content, 'utf8');
    const comparisonCommand = [
      'import sys',
      'from tests.utils.svg_compare import compare_svgs',
      'result = compare_svgs(sys.argv[1], sys.argv[2])',
      'print(result.message)',
      'print("\\n".join(result.differences))',
      'raise SystemExit(0 if result.equal else 1)'
    ].join(';');
    const semanticComparison = spawnSync(
      process.env.GBDRAW_PYTHON || 'python',
      ['-c', comparisonCommand, initialPath, generatedPath],
      { cwd: repoRoot, encoding: 'utf8' }
    );
    expect(
      semanticComparison.status,
      `${semanticComparison.stdout}\n${semanticComparison.stderr}`
    ).toBe(0);

    const worker = await getDiagramWorkerActivity(page);
    expect(worker.constructions).toBe(1);
    expect(worker.initializations).toBe(1);
    expect(worker.runs).toBe(1);
    expect(worker.settledInitializations).toBe(1);
    expect(worker.settledRuns).toBe(1);
  } finally {
    await context.close();
  }
});

// Two rings of one sequence share one raw LOSAT cache key; the saved Session
// keeps both rings and one entry per raw key (first row's filename, as the CLI
// cache writes it), restores without LOSAT, and the CLI replays it.
test('two rings of one sequence save and restore with one cached LOSAT table', async ({
  browser
}, testInfo) => {
  test.setTimeout(240_000);
  const context = await browser.newContext();
  const page = await context.newPage();
  try {
    await page.addInitScript(() => {
      const probe = { losatCalls: 0 };
      probe.read = async () => {
        const { state } = await import('/gbdraw/web/js/state.js');
        const app = window.__GBDRAW_APP__;
        const content = String(app.results[app.selectedResultIndex]?.content || '');
        const documentRoot = new DOMParser().parseFromString(content, 'image/svg+xml');
        return {
          errorLog: app.errorLog ? JSON.parse(JSON.stringify(app.errorLog)) : null,
          content,
          comparisonNames: state.files.c_conservation_fastas.map(({ name }) => name),
          labels: state.activeDrawing().circularConservation.series.map(({ label }) => label),
          cacheKeys: Array.from(state.losatCache.value.keys()),
          cacheInfo: state.losatCacheInfo.value.map(({ key, filename }) => ({ key, filename })),
          slots: [...documentRoot.querySelectorAll(
            '[data-gbdraw-slot-renderer="sequence_conservation"]'
          )].map((slot) => slot.getAttribute('data-track-label')),
          losatCalls: probe.losatCalls
        };
      };
      window.__GBDRAW_SHARED_RING_PROBE__ = probe;
      window.__GBDRAW_LOSAT_EXECUTOR__ = async () => {
        probe.losatCalls += 1;
        throw new Error('The cached ring replay must not execute LOSAT.');
      };
    });
    await openInstrumentedApp(page);
    const generated = await evaluateWithRetainedPromise(page, async (sessionText) => {
      const file = new File([sessionText], 'synthetic_conservation.gbdraw-session.json', {
        type: 'application/json', lastModified: 0
      });
      window.__GBDRAW_LAZY_SESSION_PROBE__.ignoreFile(file);
      const app = window.__GBDRAW_APP__;
      const imported = await app.importSession({ target: { files: [file], value: '' } });
      const { state } = await import('/gbdraw/web/js/state.js');
      const { readFileText } = await import('/gbdraw/web/js/services/file-content-cache.js');
      const copy = new File(
        [await readFileText(state.files.c_conservation_fastas[1])],
        'comparison-c-copy.fasta',
        { type: 'text/plain', lastModified: 0 }
      );
      app.addCircularConservationComparisonFile({ target: { files: [copy], value: '' } });
      const result = await app.runAnalysis();
      return {
        imported: imported?.status,
        result,
        ...(await window.__GBDRAW_SHARED_RING_PROBE__.read())
      };
    }, neutralConservationSessionText);
    expect(generated.imported).toBe('ok');
    expect(generated.result, JSON.stringify(generated.errorLog)).toEqual({ status: 'ok' });
    expect(generated.losatCalls).toBe(0);
    expect(generated.labels).toEqual([
      'comparison-b', 'comparison-c', 'comparison-d', 'comparison-c-copy'
    ]);
    expect(generated.slots).toEqual(generated.labels);
    expect(generated.cacheKeys).toHaveLength(3);
    expect(generated.cacheInfo.map(({ filename }) => filename)).toEqual([
      'comparison-b.circular_conservation.losatn.tsv',
      'comparison-c.circular_conservation.losatn.tsv',
      'comparison-d.circular_conservation.losatn.tsv',
      'comparison-c-copy.circular_conservation.losatn.tsv'
    ]);
    expect(generated.cacheInfo[3].key).toBe(generated.cacheInfo[1].key);

    const downloadPromise = page.waitForEvent('download', { timeout: 120_000 });
    const saved = await evaluateWithRetainedPromise(page, async () => {
      const app = window.__GBDRAW_APP__;
      app.sessionTitle = 'shared-ring-cache';
      return { result: await app.saveSessionWithTitle(), errorLog: app.errorLog };
    });
    expect(saved.result.status, JSON.stringify(saved.errorLog)).toBe('saved');
    const savedBytes = readFileSync(await (await downloadPromise).path());
    const savedText = savedBytes[0] === 0x1f
      ? gunzipSync(savedBytes).toString('utf8')
      : savedBytes.toString('utf8');
    const savedPath = testInfo.outputPath('shared-ring-cache.gbdraw-session.json');
    writeFileSync(savedPath, savedText, 'utf8');
    const savedRows = generated.cacheInfo.slice(0, 3);
    expect(JSON.parse(savedText).losatCache.entries.map(({ key, filename }) => ({ key, filename })))
      .toEqual(savedRows);

    const restored = await evaluateWithRetainedPromise(page, async (sessionText) => {
      const file = new File([sessionText], 'shared-ring-cache.gbdraw-session.json', {
        type: 'application/json', lastModified: 0
      });
      window.__GBDRAW_LAZY_SESSION_PROBE__.ignoreFile(file);
      const app = window.__GBDRAW_APP__;
      const imported = await app.importSession({ target: { files: [file], value: '' } });
      const afterImport = await window.__GBDRAW_SHARED_RING_PROBE__.read();
      const result = await app.runAnalysis();
      return {
        imported: imported?.status,
        afterImport,
        result,
        ...(await window.__GBDRAW_SHARED_RING_PROBE__.read())
      };
    }, savedText);
    expect(restored.imported, JSON.stringify(restored.afterImport.errorLog)).toBe('ok');
    expect(restored.afterImport.comparisonNames).toEqual(generated.comparisonNames);
    expect(restored.afterImport.cacheKeys).toEqual(generated.cacheKeys);
    expect(restored.afterImport.cacheInfo).toEqual(savedRows);
    expect(restored.result, JSON.stringify(restored.errorLog)).toEqual({ status: 'ok' });
    expect(restored.losatCalls).toBe(0);
    expect(restored.labels).toEqual(generated.labels);
    expect(restored.slots).toEqual(generated.labels);
    expect(restored.cacheInfo).toEqual(generated.cacheInfo);

    const generatedPath = testInfo.outputPath('shared-ring-generated.svg');
    const restoredPath = testInfo.outputPath('shared-ring-restored.svg');
    writeFileSync(generatedPath, generated.content, 'utf8');
    writeFileSync(restoredPath, restored.content, 'utf8');
    const comparison = spawnSync(process.env.GBDRAW_PYTHON || 'python', ['-c', [
      'import sys',
      'from tests.utils.svg_compare import compare_svgs',
      'result = compare_svgs(sys.argv[1], sys.argv[2])',
      'print(result.message)',
      'print("\\n".join(result.differences))',
      'raise SystemExit(0 if result.equal else 1)'
    ].join(';'), generatedPath, restoredPath], { cwd: repoRoot, encoding: 'utf8' });
    expect(comparison.status, `${comparison.stdout}\n${comparison.stderr}`).toBe(0);

    const cliDir = testInfo.outputPath('cli');
    require('node:fs').mkdirSync(cliDir, { recursive: true });
    const cli = spawnSync(process.env.GBDRAW_PYTHON || 'python', [
      '-m', 'gbdraw.cli', 'circular', '--session', savedPath, '-o', 'shared-ring', '-f', 'svg'
    ], {
      cwd: cliDir,
      encoding: 'utf8',
      timeout: 120_000,
      env: { ...process.env, PYTHONPATH: repoRoot }
    });
    expect(cli.status, `${cli.stdout}\n${cli.stderr}`).toBe(0);
    const cliSvg = readFileSync(join(cliDir, 'shared-ring.svg'), 'utf8');
    expect([...cliSvg.matchAll(
      /data-gbdraw-slot-renderer="sequence_conservation"[^>]*data-track-label="([^"]*)"/g
    )].map((match) => match[1])).toEqual(generated.labels);
  } finally {
    await context.close();
  }
});

test('a GenBank ring file reuses the FASTA ring search, records the Web runtime, and gives a runnable recipe', async ({
  browser
}, testInfo) => {
  test.setTimeout(240_000);
  const context = await browser.newContext();
  const page = await context.newPage();
  try {
    await page.addInitScript(() => {
      window.__GBDRAW_RING_PROBE__ = { losatCalls: 0, rows: null };
      window.__GBDRAW_LOSAT_EXECUTOR__ = async (jobs) => {
        window.__GBDRAW_RING_PROBE__.losatCalls += 1;
        const rows = window.__GBDRAW_RING_PROBE__.rows;
        if (!rows) throw new Error('The cached ring replay must not execute LOSAT.');
        return jobs.map((job) => ({ cacheKey: job.cacheKey, text: rows[job.cacheKey] }));
      };
    });
    await openInstrumentedApp(page);
    const accept = await page.locator('input[type="file"][accept*=".ddbj"]').evaluateAll(
      (inputs) => inputs.map((input) => input.getAttribute('accept'))
    );
    const loaded = await evaluateWithRetainedPromise(page, async (sessionText) => {
      const file = new File([sessionText], 'synthetic_conservation.gbdraw-session.json', {
        type: 'application/json', lastModified: 0
      });
      window.__GBDRAW_LAZY_SESSION_PROBE__.ignoreFile(file);
      const result = await window.__GBDRAW_APP__.importSession({ target: { files: [file], value: '' } });
      const app = window.__GBDRAW_APP__;
      const fastaRun = await app.runAnalysis();
      return {
        status: result?.status,
        fastaRun,
        content: String(app.results[app.selectedResultIndex]?.content || '')
      };
    }, neutralConservationSessionText);
    expect(loaded.status).toBe('ok');
    expect(loaded.fastaRun).toEqual({ status: 'ok' });

    const genbank = await evaluateWithRetainedPromise(page, async () => {
      const { state } = await import('/gbdraw/web/js/state.js');
      const { readFileText } = await import('/gbdraw/web/js/services/file-content-cache.js');
      const app = window.__GBDRAW_APP__;
      const fasta = await readFileText(state.files.c_conservation_fastas[1]);
      const [header, ...body] = fasta.trim().split(/\r?\n/);
      const recordId = header.slice(1).split(/\s+/)[0];
      const sequence = body.join('').toLowerCase();
      const origin = [];
      for (let start = 0; start < sequence.length; start += 60) {
        const chunk = sequence.slice(start, start + 60).match(/.{1,10}/g).join(' ');
        origin.push(`${String(start + 1).padStart(9)} ${chunk}`);
      }
      const flatFile = [
        `LOCUS       ${recordId.padEnd(16)} ${String(sequence.length).padStart(11)} bp    DNA     linear   UNK 01-JAN-1980`,
        'DEFINITION  synthetic comparison c.',
        `ACCESSION   ${recordId}`,
        `VERSION     ${recordId}`,
        'FEATURES             Location/Qualifiers',
        'ORIGIN',
        ...origin,
        '//',
        ''
      ].join('\n');
      const files = [...state.files.c_conservation_fastas];
      files[1] = new File([flatFile], 'comparison-c.gbk', { type: 'text/plain', lastModified: 0 });
      state.files.c_conservation_fastas = files;
      const result = await app.runAnalysis();
      const content = String(app.results[app.selectedResultIndex]?.content || '');
      const entries = Array.from(state.losatCache.value.entries());
      return {
        result,
        errorLog: app.errorLog,
        content,
        losatCalls: window.__GBDRAW_RING_PROBE__.losatCalls,
        labels: state.activeDrawing().circularConservation.series.map(({ label }) => label),
        rows: Object.fromEntries(entries.map(([key, entry]) => [key, entry.text])),
        cacheKeys: entries.map(([key]) => key)
      };
    });
    expect(genbank.result, JSON.stringify(genbank.errorLog)).toEqual({ status: 'ok' });
    expect(genbank.losatCalls).toBe(0);
    expect(genbank.labels).toEqual(['comparison-b', 'comparison-c', 'comparison-d']);
    expect(genbank.cacheKeys).toHaveLength(3);

    const compareSvgs = (leftName, left, rightName, right) => {
      const leftPath = testInfo.outputPath(leftName);
      const rightPath = testInfo.outputPath(rightName);
      writeFileSync(leftPath, left, 'utf8');
      writeFileSync(rightPath, right, 'utf8');
      return spawnSync(process.env.GBDRAW_PYTHON || 'python', ['-c', [
        'import sys',
        'from tests.utils.svg_compare import compare_svgs',
        'result = compare_svgs(sys.argv[1], sys.argv[2])',
        'print(result.message)',
        'print("\\n".join(result.differences))',
        'raise SystemExit(0 if result.equal else 1)'
      ].join(';'), leftPath, rightPath], { cwd: repoRoot, encoding: 'utf8' });
    };
    const genbankComparison = compareSvgs('fasta-ring.svg', loaded.content, 'genbank-ring.svg', genbank.content);
    expect(genbankComparison.status, `${genbankComparison.stdout}\n${genbankComparison.stderr}`).toBe(0);

    // A fresh search with the same rows records the Web runtime (D9/D10).
    await page.evaluate(async (rows) => {
      const { state } = await import('/gbdraw/web/js/state.js');
      window.__GBDRAW_RING_PROBE__.rows = rows;
      state.activeDrawing().losat.executionMode = 'serial';
      state.losatCache.value = new Map();
      state.losatCacheInfo.value = [];
      window.__GBDRAW_RING_PROBE__.searchRun = null;
      window.__GBDRAW_APP__.runAnalysis().then(
        (result) => { window.__GBDRAW_RING_PROBE__.searchRun = result; },
        (error) => { window.__GBDRAW_RING_PROBE__.searchRun = { status: 'threw', message: String(error) }; }
      );
    }, genbank.rows);
    await expect.poll(
      () => page.evaluate(() => window.__GBDRAW_RING_PROBE__.searchRun),
      { timeout: 120_000 }
    ).not.toBeNull();
    const searched = await evaluateWithRetainedPromise(page, async () => {
      const { state } = await import('/gbdraw/web/js/state.js');
      const { readFileText } = await import('/gbdraw/web/js/services/file-content-cache.js');
      const app = window.__GBDRAW_APP__;
      return {
        result: window.__GBDRAW_RING_PROBE__.searchRun,
        errorLog: app.errorLog,
        content: String(app.results[app.selectedResultIndex]?.content || ''),
        losatCalls: window.__GBDRAW_RING_PROBE__.losatCalls,
        runtimes: Array.from(state.losatCache.value.values()).map(({ runtime }) => runtime),
        runInfo: JSON.parse(JSON.stringify(app.lastRunInfo)),
        files: {
          'reference-a.gb': await readFileText(state.files.c_gb),
          'comparison-b.fasta': await readFileText(state.files.c_conservation_fastas[0]),
          'comparison-c.gbk': await readFileText(state.files.c_conservation_fastas[1]),
          'comparison-d.fasta': await readFileText(state.files.c_conservation_fastas[2])
        },
        tables: Object.fromEntries(state.losatCacheInfo.value.map(({ key, filename }) => (
          [filename, state.losatCache.value.get(key)?.text || '']
        )))
      };
    });
    expect(searched.result, JSON.stringify(searched.errorLog)).toEqual({ status: 'ok' });
    expect(searched.losatCalls).toBe(1);
    expect(searched.runtimes).toEqual(Array.from({ length: 3 }, () => (
      { kind: 'losat', source: 'wasm', version: null, program: 'blastn' }
    )));
    expect(searched.runInfo.losatRuntimes.map(({ text }) => text)).toEqual([
      'blastn: LOSAT, version not recorded (wasm)'
    ]);
    await page.getByRole('button', { name: /Run info$/ }).click({ timeout: 15_000 });
    await expect(page.locator('[data-run-info-search-runtimes] li')).toHaveText([
      'blastn: LOSAT, version not recorded (wasm)'
    ]);
    const searchedComparison = compareSvgs('fasta-ring-2.svg', loaded.content, 'searched-ring.svg', searched.content);
    expect(searchedComparison.status, `${searchedComparison.stdout}\n${searchedComparison.stderr}`).toBe(0);

    // The Source recipe keeps the table form (D14) and names the GenBank file
    // with --conservation_sequence; the CLI runs it.
    const args = searched.runInfo.sourceRecipe.commandArgs;
    const sequenceIndex = args.indexOf('--conservation_sequence');
    expect(sequenceIndex).toBeGreaterThan(-1);
    expect(args.slice(sequenceIndex + 1, sequenceIndex + 4)).toEqual([
      'comparison-b.fasta', 'comparison-c.gbk', 'comparison-d.fasta'
    ]);
    const workDir = testInfo.outputPath('recipe');
    require('node:fs').mkdirSync(workDir, { recursive: true });
    Object.entries({ ...searched.files, ...searched.tables }).forEach(([name, text]) => {
      writeFileSync(join(workDir, name), text, 'utf8');
    });
    const recipeArgs = args[0] === 'gbdraw' ? args.slice(1) : args;
    writeFileSync(join(workDir, 'recipe-args.json'), JSON.stringify(recipeArgs), 'utf8');
    const recipe = spawnSync(process.env.GBDRAW_PYTHON || 'python', [
      '-m', 'gbdraw.cli', ...recipeArgs.map(String)
    ], {
      cwd: workDir,
      encoding: 'utf8',
      timeout: 120_000,
      env: { ...process.env, PYTHONPATH: repoRoot }
    });
    expect(recipe.status, `${JSON.stringify(recipeArgs)}\n${recipe.stdout}\n${recipe.stderr}`).toBe(0);
  } finally {
    await context.close();
  }
});

// D12: a ring row added from a GenBank or DDBJ file without a typed label is
// named like the CLI names it (first record's DEFINITION, then organism); a
// FASTA row keeps the Web file-name default, a typed label wins, and the
// labels survive Session save and restore.
test('GenBank and DDBJ ring rows added without a label take the CLI default label', async ({
  browser
}, testInfo) => {
  test.setTimeout(240_000);
  const context = await browser.newContext();
  const page = await context.newPage();
  try {
    await page.addInitScript(() => {
      window.__GBDRAW_LOSAT_EXECUTOR__ = async () => {
        throw new Error('The cached ring replay must not execute LOSAT.');
      };
    });
    await openInstrumentedApp(page);
    const loaded = await evaluateWithRetainedPromise(page, async (sessionText) => {
      const file = new File([sessionText], 'synthetic_conservation.gbdraw-session.json', {
        type: 'application/json', lastModified: 0
      });
      window.__GBDRAW_LAZY_SESSION_PROBE__.ignoreFile(file);
      const result = await window.__GBDRAW_APP__.importSession({ target: { files: [file], value: '' } });
      const { state } = await import('/gbdraw/web/js/state.js');
      const { readFileText } = await import('/gbdraw/web/js/services/file-content-cache.js');
      return { status: result?.status, fasta: await readFileText(state.files.c_conservation_fastas[1]) };
    }, neutralConservationSessionText);
    expect(loaded.status).toBe('ok');

    // Flat files of comparison c's sequence: one raw key, so Generate reuses the
    // Session's cached rows and LOSAT does not run.
    const [header, ...body] = loaded.fasta.trim().split(/\r?\n/);
    const recordId = header.slice(1).split(/\s+/)[0];
    const sequence = body.join('').toLowerCase();
    const origin = [];
    for (let start = 0; start < sequence.length; start += 60) {
      const chunk = sequence.slice(start, start + 60).match(/.{1,10}/g).join(' ');
      origin.push(`${String(start + 1).padStart(9)} ${chunk}`);
    }
    const flatFile = ({ definition, division }) => [
      `LOCUS       ${recordId.padEnd(16)} ${String(sequence.length).padStart(11)} bp    DNA     linear   ${division} 01-JAN-1980`,
      `DEFINITION  ${definition}`,
      `ACCESSION   ${recordId}`,
      `VERSION     ${recordId}`,
      'KEYWORDS    .',
      'SOURCE      Synthetic organism c',
      '  ORGANISM  Synthetic organism c',
      '            Unclassified.',
      'FEATURES             Location/Qualifiers',
      'ORIGIN',
      ...origin,
      '//',
      ''
    ].join('\n');
    const ringFiles = [
      { name: 'ring-genbank.gbk', text: flatFile({ definition: 'synthetic comparison c.', division: 'UNK' }) },
      { name: 'ring-ddbj.ddbj', text: flatFile({ definition: '.', division: 'SYN' }) },
      { name: 'ring-typed.gbk', text: flatFile({ definition: 'synthetic comparison c.', division: 'UNK' }), typed: 'Typed ring' },
      { name: 'ring-fasta.fa', text: loaded.fasta }
    ];

    // The CLI default for the same files (gbdraw.io comparison reader).
    const cliDir = testInfo.outputPath('ring-files');
    require('node:fs').mkdirSync(cliDir, { recursive: true });
    ringFiles.forEach(({ name, text }) => writeFileSync(join(cliDir, name), text, 'utf8'));
    const cli = spawnSync(process.env.GBDRAW_PYTHON || 'python', ['-c', [
      'import json, sys',
      'from gbdraw.io.comparison_sequences import read_comparison_sequence_file',
      'print(json.dumps([read_comparison_sequence_file(path).label for path in sys.argv[1:]]))'
    ].join(';'), ...ringFiles.map(({ name }) => join(cliDir, name))], {
      cwd: repoRoot, encoding: 'utf8', env: { ...process.env, PYTHONPATH: repoRoot }
    });
    expect(cli.status, cli.stderr).toBe(0);
    const cliLabels = JSON.parse(cli.stdout);
    expect(cliLabels).toEqual([
      'synthetic comparison c', 'Synthetic organism c', 'synthetic comparison c', 'ring-fasta.fa'
    ]);

    const run = await evaluateWithRetainedPromise(page, async (inputs) => {
      const { state } = await import('/gbdraw/web/js/state.js');
      const app = window.__GBDRAW_APP__;
      for (const { name, text, typed } of inputs) {
        const file = new File([text], name, { type: 'text/plain', lastModified: 0 });
        app.addCircularConservationComparisonFile({ target: { files: [file], value: '' } });
        // Typed before the reader answers: the typed label wins.
        if (typed) state.activeDrawing().circularConservation.series[state.activeDrawing().circularConservation.series.length - 1].label = typed;
      }
      // Generate waits for the pending reads, then draws the rows' labels.
      const result = await app.runAnalysis();
      return {
        result,
        errorLog: app.errorLog,
        labels: state.activeDrawing().circularConservation.series.map(({ label }) => label)
      };
    }, ringFiles);
    expect(run.result, JSON.stringify(run.errorLog)).toEqual({ status: 'ok' });
    const expectedLabels = [
      'comparison-b', 'comparison-c', 'comparison-d',
      cliLabels[0], cliLabels[1], 'Typed ring', 'ring-fasta'
    ];
    expect(run.labels).toEqual(expectedLabels);

    // Session save and restore keep the labels, also for the five rings that
    // share comparison c's sequence and raw cache key.
    const downloadPromise = page.waitForEvent('download', { timeout: 120_000 });
    const saved = await evaluateWithRetainedPromise(page, async () => {
      const app = window.__GBDRAW_APP__;
      app.sessionTitle = 'ring-default-labels';
      return {
        result: await app.saveSessionWithTitle(),
        errorLog: app.errorLog
      };
    });
    const savedLabels = expectedLabels;
    expect(saved.result.status, JSON.stringify(saved.errorLog)).toBe('saved');
    const savedPath = await (await downloadPromise).path();
    const savedBytes = readFileSync(savedPath);
    const savedText = savedBytes[0] === 0x1f ? gunzipSync(savedBytes).toString('utf8') : savedBytes.toString('utf8');
    const findLabels = (value) => {
      if (!value || typeof value !== 'object') return null;
      if (Array.isArray(value.conservationLabels)) return value.conservationLabels;
      for (const child of Object.values(value)) {
        const found = findLabels(child);
        if (found) return found;
      }
      return null;
    };
    expect(findLabels(JSON.parse(savedText))).toEqual(savedLabels);
    const restored = await evaluateWithRetainedPromise(page, async (sessionText) => {
      const file = new File([sessionText], 'ring-default-labels.gbdraw-session.json', {
        type: 'application/json', lastModified: 0
      });
      window.__GBDRAW_LAZY_SESSION_PROBE__.ignoreFile(file);
      const result = await window.__GBDRAW_APP__.importSession({ target: { files: [file], value: '' } });
      const { state } = await import('/gbdraw/web/js/state.js');
      return {
        status: result?.status,
        message: result?.message,
        errorLog: window.__GBDRAW_APP__.errorLog,
        labels: state.activeDrawing().circularConservation.series.map(({ label }) => label)
      };
    }, savedText);
    expect(restored.status, `${restored.message} ${JSON.stringify(restored.errorLog)}`).toBe('ok');
    expect(restored.labels).toEqual(savedLabels);
  } finally {
    await context.close();
  }
});

test('Generate materializes only required resources and reuses one Worker', async ({ page }) => {
  test.setTimeout(300_000);
  await openInstrumentedApp(page);
  await armHistoryCompletion(page);
  expect(await loadSyntheticSession(page)).toMatchObject({
    status: 'ok',
    degradedRecovery: false
  });
  const preview = await probeSnapshot(page);
  expect(preview.structural).toEqual(ZERO_PREVIEW_METRICS);
  expect(preview.historyBaseline).toEqual(ZERO_ARTIFACT_HISTORY_BASELINE);

  const generated = await evaluateWithRetainedPromise(page, async () => {
    const app = window.__GBDRAW_APP__;
    const history = window.__GBDRAW_HISTORY__;
    const svgIdentity = (content) => {
      const documentElement = new DOMParser()
        .parseFromString(String(content || ''), 'image/svg+xml')
        .documentElement;
      const normalized = new XMLSerializer().serializeToString(documentElement);
      let hash = 2166136261;
      for (let index = 0; index < normalized.length; index += 1) {
        hash ^= normalized.charCodeAt(index);
        hash = Math.imul(hash, 16777619);
      }
      return { length: normalized.length, hash: hash >>> 0 };
    };
    const snapshotResultState = () => ({
      mode: app.mode,
      generatedMode: app.generatedMode,
      names: app.results.map((result) => String(result?.name || '')),
      svgIdentities: app.results.map((result) => svgIdentity(result?.content)),
      selectedResultIndex: app.selectedResultIndex,
      resultCount: app.results.length,
      featureCount: app.extractedFeatures.length,
      orthogroupCount: app.orthogroups.length
    });
    const historyDelta = (before, after) => Object.fromEntries([
      'artifactCheckpointBuilds',
      'artifactCheckpointSignatureComputations',
      'historySvgBytes',
      'checkpointEstimatedBytes',
      'generatedArtifactFullCloneCount',
      'generatedArtifactFullSerializationCount',
      'manualCancelFullArtifactSnapshotBuildCount',
      'artifactHandleBeforeBuildCount',
      'artifactHandleAfterBuildCount',
      'artifactFingerprintComparisonCount',
      'artifactReplacementHistoryEntryCount'
    ].map((name) => [name, Number(after[name] || 0) - Number(before[name] || 0)]));
    const loaded = snapshotResultState();
    const beforeFirst = history.getDiagnostics();
    const first = await app.runAnalysis();
    const afterFirst = history.getDiagnostics();
    const generatedState = snapshotResultState();
    const firstUndoCount = history.getUndoCount();
    const firstRedoCount = history.getRedoCount();
    const firstUndoLabel = history.undoLabel();
    const undone = await history.undo();
    const restoredLoadedState = snapshotResultState();
    const redoLabelAfterUndo = history.redoLabel();
    const redone = await history.redo();
    const restoredGeneratedState = snapshotResultState();
    const {
      DIAGRAM_HELPER_OPERATIONS,
      runDiagramHelperOperation
    } = await import('/gbdraw/web/js/services/diagram-generation.js');
    const helper = await runDiagramHelperOperation(
      DIAGRAM_HELPER_OPERATIONS.MEASURE_LEGEND_TEXT,
      { caption: 'lazy worker reuse', fontFamily: 'Arial', fontSize: 14 }
    );
    const beforeSecond = history.getDiagnostics();
    const undoCountBeforeSecond = history.getUndoCount();
    const redoCountBeforeSecond = history.getRedoCount();
    const second = await app.runAnalysis();
    const afterSecond = history.getDiagnostics();
    const undoCountAfterSecond = history.getUndoCount();
    const redoCountAfterSecond = history.getRedoCount();
    const generateProbe = window.__GBDRAW_LAZY_SESSION_PROBE__.snapshot();
    const editableFeature = app.extractedFeatures.find((feature) => (
      app.canEditFeatureColor(feature)
    ));
    if (!editableFeature) {
      throw new Error('The synthetic session has no editable feature for mutation isolation.');
    }
    const currentColor = String(app.getFeatureColorValue(editableFeature) || '').toLowerCase();
    const mutationColor = currentColor === '#123456' ? '#654321' : '#123456';
    const editApplied = await app.setFeatureColorValue(editableFeature, mutationColor);
    const editedState = snapshotResultState();
    const undoEdit = await history.undo();
    const afterUndoEdit = snapshotResultState();
    const undoGenerate = await history.undo();
    const afterUndoGenerate = snapshotResultState();
    const redoGenerate = await history.redo();
    const afterRedoGenerate = snapshotResultState();
    const redoEdit = await history.redo();
    const afterRedoEdit = snapshotResultState();
    return {
      first,
      second,
      loaded,
      generatedState,
      firstUndoCount,
      firstRedoCount,
      firstUndoLabel,
      redoLabelAfterUndo,
      undone,
      restoredLoadedState,
      redone,
      restoredGeneratedState,
      firstHistory: historyDelta(beforeFirst, afterFirst),
      secondHistory: historyDelta(beforeSecond, afterSecond),
      undoCountBeforeSecond,
      redoCountBeforeSecond,
      undoCountAfterSecond,
      redoCountAfterSecond,
      generateProbe,
      mutationIsolation: {
        editApplied,
        editedState,
        undoEdit,
        afterUndoEdit,
        undoGenerate,
        afterUndoGenerate,
        redoGenerate,
        afterRedoGenerate,
        redoEdit,
        afterRedoEdit
      },
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
  expect(generated.firstUndoCount).toBe(1);
  expect(generated.firstRedoCount).toBe(0);
  expect(generated.firstUndoLabel).toBe('Generate diagram');
  expect(generated.undone).toBe(true);
  expect(generated.restoredLoadedState, JSON.stringify(generated, null, 2)).toEqual(
    generated.loaded
  );
  expect(generated.redone).toBe(true);
  expect(generated.redoLabelAfterUndo).toBe('Generate diagram');
  expect(generated.restoredGeneratedState).toEqual(generated.generatedState);
  expect(generated.restoredGeneratedState.selectedResultIndex).toBeGreaterThanOrEqual(0);
  expect(generated.restoredGeneratedState.selectedResultIndex).toBeLessThan(
    generated.restoredGeneratedState.resultCount
  );
  expect(generated.firstHistory).toEqual({
    artifactCheckpointBuilds: 0,
    artifactCheckpointSignatureComputations: 0,
    historySvgBytes: 0,
    checkpointEstimatedBytes: 0,
    generatedArtifactFullCloneCount: 0,
    generatedArtifactFullSerializationCount: 0,
    manualCancelFullArtifactSnapshotBuildCount: 0,
    artifactHandleBeforeBuildCount: 1,
    artifactHandleAfterBuildCount: 1,
    artifactFingerprintComparisonCount: 1,
    artifactReplacementHistoryEntryCount: 1
  });
  expect(generated.secondHistory).toEqual({
    ...generated.firstHistory,
    artifactReplacementHistoryEntryCount: 0
  });
  expect(generated.undoCountAfterSecond).toBe(generated.undoCountBeforeSecond);
  expect(generated.redoCountAfterSecond).toBe(generated.redoCountBeforeSecond);
  expect(generated.mutationIsolation.editApplied).toBe(true);
  expect(generated.mutationIsolation.editedState.svgIdentities).not.toEqual(
    generated.generatedState.svgIdentities
  );
  expect(generated.mutationIsolation.undoEdit).toBe(true);
  expect(generated.mutationIsolation.afterUndoEdit.resultCount).toBe(
    generated.generatedState.resultCount
  );
  expect(generated.mutationIsolation.undoGenerate).toBe(true);
  expect(generated.mutationIsolation.afterUndoGenerate).toEqual(generated.loaded);
  expect(generated.mutationIsolation.redoGenerate).toBe(true);
  expect(generated.mutationIsolation.afterRedoGenerate).toEqual(generated.generatedState);
  expect(generated.mutationIsolation.redoEdit).toBe(true);
  expect(generated.mutationIsolation.afterRedoEdit.resultCount).toBe(
    generated.mutationIsolation.editedState.resultCount
  );

  const snapshot = generated.generateProbe;
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
  expect(snapshot.structural.base64DecodeCount).toBe(
    snapshot.structural.resourceByteReadCount
  );
  expect(snapshot.structural.base64DecodeCount).toBeGreaterThanOrEqual(
    materializedIds.size
  );

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
  // Two Generates, and the automatic rerender of the color edit, which adds a
  // Legend row (OV-42, OV-43, #857).
  expect(worker.runs).toBe(3);
  expect(worker.instances).toHaveLength(1);
  expect(worker.instances[0].terminated).toBe(false);
});

test('render-only Generate reuses preparation and remains undoable', async ({ page }) => {
  test.setTimeout(300_000);
  await openInstrumentedApp(page);
  await armHistoryCompletion(page);
  expect(await loadSyntheticSession(page)).toMatchObject({
    status: 'ok',
    degradedRecovery: false
  });

  const first = await evaluateWithRetainedPromise(page, async () => {
    const app = window.__GBDRAW_APP__;
    const history = window.__GBDRAW_HISTORY__;
    window.__GBDRAW_LAZY_SESSION_PROBE__.reset();
    const before = history.getDiagnostics();
    const result = await app.runAnalysis();
    const content = String(app.results?.[app.selectedResultIndex]?.content || '');
    const after = history.getDiagnostics();
    window.__GBDRAW_PREPARED_CACHE_A__ = content;
    window.__GBDRAW_PREPARED_CACHE_SCALE_A__ = Boolean(app.form.show_scale);
    return {
      result,
      contentLength: content.length,
      showScale: Boolean(app.form.show_scale),
      undoCount: history.getUndoCount(),
      historyEntries: Number(after.artifactReplacementHistoryEntryCount || 0)
        - Number(before.artifactReplacementHistoryEntryCount || 0),
      errorSummary: String(app.errorLog?.summary || '')
    };
  });
  const firstProbe = await probeSnapshot(page);
  expect(first).toMatchObject({
    result: { status: 'ok' },
    undoCount: 1,
    historyEntries: 1,
    errorSummary: ''
  });
  expect(first.contentLength).toBeGreaterThan(0);

  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    app.form.show_scale = !Boolean(app.form.show_scale);
    window.__GBDRAW_LAZY_SESSION_PROBE__.reset();
  });
  const second = await evaluateWithRetainedPromise(page, async () => {
    const app = window.__GBDRAW_APP__;
    const history = window.__GBDRAW_HISTORY__;
    const before = history.getDiagnostics();
    const result = await app.runAnalysis();
    const content = String(app.results?.[app.selectedResultIndex]?.content || '');
    const after = history.getDiagnostics();
    const undoOk = await history.undo();
    const afterUndo = String(app.results?.[app.selectedResultIndex]?.content || '');
    const redoOk = await history.redo();
    const afterRedo = String(app.results?.[app.selectedResultIndex]?.content || '');
    return {
      result,
      outputChanged: content !== window.__GBDRAW_PREPARED_CACHE_A__,
      renderOnlyValueChanged:
        Boolean(app.form.show_scale) !== window.__GBDRAW_PREPARED_CACHE_SCALE_A__,
      historyEntries: Number(after.artifactReplacementHistoryEntryCount || 0)
        - Number(before.artifactReplacementHistoryEntryCount || 0),
      undoOk,
      undoRestoredA: afterUndo === window.__GBDRAW_PREPARED_CACHE_A__,
      redoOk,
      redoRestoredB: afterRedo === content,
      undoCount: history.getUndoCount(),
      redoCount: history.getRedoCount(),
      errorSummary: String(app.errorLog?.summary || '')
    };
  });
  const secondProbe = await probeSnapshot(page);
  expect(second).toEqual({
    result: { status: 'ok' },
    outputChanged: true,
    renderOnlyValueChanged: true,
    historyEntries: 1,
    undoOk: true,
    undoRestoredA: true,
    redoOk: true,
    redoRestoredB: true,
    undoCount: 2,
    redoCount: 0,
    errorSummary: ''
  });

  const firstPython = firstProbe.lifecycle.find(
    ({ name }) => name === 'python-diagnostics'
  );
  const secondPython = secondProbe.lifecycle.find(
    ({ name }) => name === 'python-diagnostics'
  );
  expect(firstPython?.metrics).toMatchObject({
    parsedSourceCacheMissCount: 1,
    parsedSourceParseCount: 1,
    resolvedRecordCacheMissCount: 1,
    resolvedRecordBuildCount: 1,
    interactiveContextCacheMissCount: 1,
    interactiveContextBuildCount: 1,
    interactiveFeatureTraversalCount: 1,
    preparedInputCacheMutationViolationCount: 0
  });
  expect(secondPython?.metrics).toMatchObject({
    parsedSourceCacheHitCount: 1,
    parsedSourceParseCount: 0,
    resolvedRecordCacheHitCount: 1,
    resolvedRecordBuildCount: 0,
    interactiveContextCacheHitCount: 1,
    interactiveContextBuildCount: 0,
    interactiveFeatureTraversalCount: 0,
    preparedInputCacheMutationViolationCount: 0,
    featureCatalogSvgParseCount: 1,
    featureCatalogFullDomTraversalCount: 1
  });
  expect(Number(secondPython?.metrics?.preparedInputCacheRetainedBytes || 0))
    .toBeGreaterThan(0);
  expect(Number(secondPython?.timingsMs?.drawing || 0)).toBeGreaterThan(0);
  expect(Number(secondPython?.timingsMs?.svgWrite || 0)).toBeGreaterThan(0);
  expect(Number(secondProbe.hookMetrics.resourceMaterializationCount || 0)).toBe(0);
  expect(secondProbe.lifecycle.find(
    ({ name }) => name === 'worker-resource-linking-end'
  )).toMatchObject({ newlyStagedResourceBytes: 0 });

  const worker = await getDiagramWorkerActivity(page);
  expect(worker.constructions).toBe(1);
  expect(worker.initializations).toBe(1);
  expect(worker.runs).toBe(2);
  expect(worker.instances).toHaveLength(1);
  expect(worker.instances[0].terminated).toBe(false);
});

test('the real Cancel control restores the committed artifact and leaves no History entry', async ({
  page
}) => {
  test.setTimeout(300_000);
  await openInstrumentedApp(page);
  await armHistoryCompletion(page);
  expect(await loadSyntheticSession(page)).toMatchObject({
    status: 'ok',
    degradedRecovery: false
  });

  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const history = window.__GBDRAW_HISTORY__;
    window.__GBDRAW_CANCEL_BASELINE__ = {
      content: String(app.results?.[app.selectedResultIndex]?.content || ''),
      resultCount: app.results.length,
      selectedResultIndex: app.selectedResultIndex,
      featureCount: app.extractedFeatures.length,
      orthogroupCount: app.orthogroups.length,
      undoCount: history.getUndoCount(),
      redoCount: history.getRedoCount(),
      diagnostics: history.getDiagnostics()
    };
    window.__GBDRAW_CANCEL_BASELINE_REFS__ = {
      result: app.results?.[app.selectedResultIndex] || null,
      extractedFeatures: app.extractedFeatures,
      orthogroups: app.orthogroups
    };
    let releaseResponse;
    const responseGate = new Promise((resolve) => {
      releaseResponse = resolve;
    });
    window.__GBDRAW_CANCEL_RESPONSE_STARTED__ = false;
    window.__GBDRAW_RELEASE_CANCEL_RESPONSE__ = releaseResponse;
    window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse = () => {
      window.__GBDRAW_CANCEL_RESPONSE_STARTED__ = true;
      return responseGate;
    };
    window.__GBDRAW_CANCEL_RUN__ = app.runAnalysis();
  });

  await page.waitForFunction(
    () => window.__GBDRAW_CANCEL_RESPONSE_STARTED__ === true,
    null,
    { timeout: 240_000 }
  );
  await page.getByRole('button', { name: /Cancel$/ }).click();
  const canceled = await evaluateWithRetainedPromise(page, async () => {
    const result = await window.__GBDRAW_CANCEL_RUN__;
    window.__GBDRAW_RELEASE_CANCEL_RESPONSE__?.();
    delete window.__GBDRAW_TEST_HOOKS__.beforeDiagramGenerationResponse;
    const app = window.__GBDRAW_APP__;
    const history = window.__GBDRAW_HISTORY__;
    const baseline = window.__GBDRAW_CANCEL_BASELINE__;
    const baselineRefs = window.__GBDRAW_CANCEL_BASELINE_REFS__;
    const after = history.getDiagnostics();
    const delta = Object.fromEntries([
      'artifactCheckpointBuilds',
      'artifactCheckpointSignatureComputations',
      'historySvgBytes',
      'checkpointEstimatedBytes',
      'generatedArtifactFullCloneCount',
      'generatedArtifactFullSerializationCount',
      'manualCancelFullArtifactSnapshotBuildCount',
      'artifactHandleBeforeBuildCount',
      'artifactHandleAfterBuildCount',
      'artifactFingerprintComparisonCount',
      'artifactReplacementHistoryEntryCount'
    ].map((name) => [
      name,
      Number(after[name] || 0) - Number(baseline.diagnostics[name] || 0)
    ]));
    return {
      result,
      contentRestored:
        String(app.results?.[app.selectedResultIndex]?.content || '') === baseline.content,
      resultCount: app.results.length,
      selectedResultIndex: app.selectedResultIndex,
      featureCount: app.extractedFeatures.length,
      orthogroupCount: app.orthogroups.length,
      baselineFeatureCount: baseline.featureCount,
      baselineOrthogroupCount: baseline.orthogroupCount,
      sameResultObject:
        app.results?.[app.selectedResultIndex] === baselineRefs.result,
      sameExtractedFeatureOwner: app.extractedFeatures === baselineRefs.extractedFeatures,
      sameOrthogroupOwner: app.orthogroups === baselineRefs.orthogroups,
      undoCount: history.getUndoCount(),
      redoCount: history.getRedoCount(),
      processing: app.processing,
      diagnostics: delta
    };
  });
  expect(canceled).toMatchObject({
    result: { status: 'canceled' },
    contentRestored: true,
    resultCount: 1,
    selectedResultIndex: 0,
    sameResultObject: true,
    sameExtractedFeatureOwner: true,
    sameOrthogroupOwner: true,
    undoCount: 0,
    redoCount: 0,
    processing: false,
    diagnostics: {
      artifactCheckpointBuilds: 0,
      artifactCheckpointSignatureComputations: 0,
      historySvgBytes: 0,
      checkpointEstimatedBytes: 0,
      generatedArtifactFullCloneCount: 0,
      generatedArtifactFullSerializationCount: 0,
      manualCancelFullArtifactSnapshotBuildCount: 0,
      artifactHandleBeforeBuildCount: 1,
      artifactHandleAfterBuildCount: 0,
      artifactFingerprintComparisonCount: 0,
      artifactReplacementHistoryEntryCount: 0
    }
  });
  expect(canceled.featureCount).toBe(canceled.baselineFeatureCount);
  expect(canceled.featureCount).toBeGreaterThan(0);
  expect(canceled.orthogroupCount).toBe(canceled.baselineOrthogroupCount);

  const canceledWorker = await getDiagramWorkerActivity(page);
  expect(canceledWorker.runs).toBe(1);
  expect(canceledWorker.instances[0].terminated).toBe(true);

  await page.evaluate(() => window.__GBDRAW_LAZY_SESSION_PROBE__.reset());
  const retry = await evaluateWithRetainedPromise(page, async () => ({
    result: await window.__GBDRAW_APP__.runAnalysis(),
    undoCount: window.__GBDRAW_HISTORY__.getUndoCount(),
    redoCount: window.__GBDRAW_HISTORY__.getRedoCount(),
    previewVisible: Boolean(document.querySelector('.shadow-xl.origin-top > svg')),
    errorSummary: String(window.__GBDRAW_APP__.errorLog?.summary || '')
  }));
  expect(retry).toMatchObject({
    result: { status: 'ok' },
    undoCount: 1,
    redoCount: 0,
    previewVisible: true,
    errorSummary: ''
  });
  const retryProbe = await probeSnapshot(page);
  const retryPython = retryProbe.lifecycle.find(
    ({ name }) => name === 'python-diagnostics'
  );
  expect(retryPython?.metrics).toMatchObject({
    parsedSourceCacheHitCount: 0,
    parsedSourceCacheMissCount: 1,
    parsedSourceParseCount: 1,
    resolvedRecordCacheHitCount: 0,
    resolvedRecordCacheMissCount: 1,
    resolvedRecordBuildCount: 1,
    interactiveContextCacheHitCount: 0,
    interactiveContextCacheMissCount: 1,
    interactiveContextBuildCount: 1
  });
  const afterRetryWorker = await getDiagramWorkerActivity(page);
  expect(afterRetryWorker.constructions).toBe(2);
  expect(afterRetryWorker.instances[1].terminated).toBe(false);
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
  expect(rejected.message).toContain('Load a supported Session file or recreate it with the current writer.');
  expect(await page.evaluate(() => window.__GBDRAW_APP__.errorLog)).toMatchObject({
    code: 'INPUT_INVALID', stage: 'request-validation', context: { field: 'schema', reason: 'RESOURCE_SIZE' }
  });
  expect(rejected.message).not.toContain('unused-lazy-contract');
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
    /record-1-genbank \(record-1-genbank-HmmtDNA\.gbk\) contains invalid encoded data/
  );
  expect(access.second).toBe(access.first);
  expect(access.sameContent).toBe(true);
  expect(access.previewVisible).toBe(true);
  const afterAccess = await probeSnapshot(page);
  expect(afterAccess.structural.base64DecodeCount).toBe(1);
  expect(afterAccess.structural.resourceByteReadCount).toBe(1);
  expect(afterAccess.structural.workerConstructionCount).toBe(0);
  expect(afterAccess.structural.workerInitializationCount).toBe(0);

  const failedGenerate = await evaluateWithRetainedPromise(page, async () => {
    const app = window.__GBDRAW_APP__;
    const history = window.__GBDRAW_HISTORY__;
    const loadedContent = app.results[app.selectedResultIndex]?.content || '';
    const result = await app.runAnalysis();
    return {
      status: result?.status,
      previewCommitted: app.results[app.selectedResultIndex]?.content === loadedContent,
      previewVisible: Boolean(document.querySelector('.shadow-xl.origin-top > svg')),
      undoCount: history.getUndoCount(),
      redoCount: history.getRedoCount(),
      currentCheckpointAbsent: history.getCurrentCheckpoint() === null
    };
  });
  expect(failedGenerate).toEqual({
    status: 'error',
    previewCommitted: true,
    previewVisible: true,
    undoCount: 0,
    redoCount: 0,
    currentCheckpointAbsent: true
  });
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
  const regenerated = await evaluateWithRetainedPromise(page, async () => ({
    result: await window.__GBDRAW_APP__.runAnalysis(),
    errorLog: window.__GBDRAW_APP__.errorLog
  }));
  expect(regenerated.result, JSON.stringify(regenerated.errorLog)).toEqual({ status: 'ok' });
  const downloadPromise = page.waitForEvent('download', { timeout: 120_000 });
  const saveOutcome = await evaluateWithRetainedPromise(page, async () => {
    const result = await window.__GBDRAW_APP__.saveSessionWithTitle();
    return {
      result,
      errorLog: window.__GBDRAW_APP__.errorLog
    };
  });
  expect(saveOutcome.result.status, JSON.stringify(saveOutcome.errorLog)).toBe('saved');
  const roundTripPath = await (await downloadPromise).path();
  const roundTrip = JSON.parse(gunzipSync(readFileSync(roundTripPath)).toString('utf8'));
  expect(roundTrip.version).toBe(CURRENT_SESSION_VERSION);
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
