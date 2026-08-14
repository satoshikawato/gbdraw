const { test, expect } = require('@playwright/test');
const { join, resolve } = require('node:path');
const {
  getDiagramWorkerActivity,
  openApp
} = require('../helpers/app-lifecycle.cjs');

const repoRoot = resolve(process.env.GBDRAW_REPO || process.cwd());
const vibrioSession = join(
  repoRoot,
  'gbdraw',
  'web',
  'gallery',
  'sessions',
  'vibrio-harveyi-group-collinear.gbdraw-session.json.gz'
);

const PREVIEW_ZERO_METRICS = [
  'base64DecodeCount',
  'decodedByteCount',
  'sourceRecoveryCount',
  'workerConstructionCount',
  'workerInitializationCount'
];

const installStructuralProbe = async (page, runnerEvidence) => {
  await page.exposeFunction('__gbdrawCaptureVibrioEvidence', (entry) => {
    runnerEvidence.push(entry);
  });
  await page.addInitScript(() => {
    const metrics = {};
    const details = [];
    const lifecycle = [];
    const probe = {
      metrics,
      details,
      lifecycle,
      reset() {
        Object.keys(metrics).forEach((name) => delete metrics[name]);
        details.length = 0;
        lifecycle.length = 0;
      },
      snapshot() {
        return {
          metrics: { ...metrics },
          details: details.map((detail) => ({ ...detail })),
          lifecycle: lifecycle.map((event) => ({ ...event }))
        };
      }
    };
    window.__GBDRAW_VIBRIO_GENERATE_PROBE__ = probe;
    window.__GBDRAW_TEST_HOOKS__ = {
      onStructuralMetric(metric) {
        const name = String(metric?.name || '');
        metrics[name] = Number(metrics[name] || 0) + Number(metric?.value || 0);
        details.push({ ...metric, timestamp: performance.now() });
        void window.__gbdrawCaptureVibrioEvidence({ type: 'metric', value: metric });
      },
      onSessionLifecycleEvent(event) {
        lifecycle.push({ ...event });
        void window.__gbdrawCaptureVibrioEvidence({ type: 'lifecycle', value: event });
      }
    };
  });
};

const probeSnapshot = (page) => page.evaluate(() => (
  window.__GBDRAW_VIBRIO_GENERATE_PROBE__.snapshot()
));

const phaseDuration = (events, startName, endName) => {
  const start = events.find(({ name }) => name === startName);
  const end = events.find(({ name }) => name === endName);
  return start && end ? end.timestamp - start.timestamp : null;
};

const loadTimings = (events) => {
  const selection = events.find(({ name }) => name === 'sessionSelection');
  const preview = events.find(({ name }) => name === 'firstCommittedPreview');
  const ready = events.find(({ name }) => name === 'interactiveReady');
  return {
    gzipToTextMs: phaseDuration(events, 'gzip-to-text-start', 'gzip-to-text-end'),
    jsonParseMs: phaseDuration(events, 'json-parse-start', 'json-parse-end'),
    currentSessionPreflightMs: phaseDuration(
      events,
      'current-session-preflight-start',
      'current-session-preflight-end'
    ),
    svgAdmissionMs: phaseDuration(events, 'svg-admission-start', 'svg-admission-end'),
    previewMountMs: phaseDuration(events, 'preview-mount-start', 'preview-mount-end'),
    historyBaselineMs: phaseDuration(
      events,
      'history-baseline-start',
      'history-baseline-end'
    ),
    firstPreviewMs: selection && preview ? preview.timestamp - selection.timestamp : null,
    interactiveReadyMs: selection && ready ? ready.timestamp - selection.timestamp : null
  };
};

const workerErrors = (activity) => (
  (activity?.instances || []).flatMap((instance) => instance.errors || [])
);

const workerFailureSettlements = (activity) => (
  (activity?.instances || []).flatMap((instance) => (
    (instance.settlements || []).filter(({ type, ok }) => type === 'run' && !ok)
  ))
);

const errorText = (diagnostic) => [
  diagnostic?.evaluationError,
  diagnostic?.outcome?.thrown,
  diagnostic?.outcome?.errorSummary,
  ...(diagnostic?.outcome?.errorDetails || []),
  ...workerFailureSettlements(diagnostic?.worker).map(({ error }) => error),
  ...workerErrors(diagnostic?.worker).map(({ message }) => message)
].filter(Boolean).join('\n');

const classifyGeneration = (diagnostic) => {
  const message = errorText(diagnostic);
  if (diagnostic.terminal.contextClosed) return 'browser-context-closure';
  if (diagnostic.terminal.pageCrashed || /target (?:page )?crashed/i.test(message)) {
    return 'page-crash';
  }
  if (/RangeError|invalid string length|string length overflow/i.test(message)) {
    return 'range-or-invalid-string-length';
  }
  if (/out of memory|memory access out of bounds|Cannot enlarge memory|WebAssembly\.Memory/i.test(message)) {
    return 'wasm-or-pyodide-out-of-memory';
  }
  if (
    workerErrors(diagnostic.worker).length > 0
    || workerFailureSettlements(diagnostic.worker).length > 0
  ) return 'worker-error';
  if (diagnostic.terminal.pageClosed) return 'page-closure';
  if (diagnostic.evaluationError || diagnostic.outcome?.thrown) return 'evaluation-error';
  if (diagnostic.outcome?.status === 'error') return 'returned-application-error';
  if (diagnostic.outcome?.status === 'ok') return 'success';
  return 'unknown';
};

const invokeGeneration = async (page, terminal, runnerEvidence, terminalSignal) => {
  let outcome = null;
  let evaluationError = '';
  const evaluation = page.evaluate(async () => {
      const app = window.__GBDRAW_APP__;
      const started = performance.now();
      try {
        const result = await app.runAnalysis();
        const selected = app.results?.[app.selectedResultIndex];
        return {
          status: String(result?.status || ''),
          elapsedMs: performance.now() - started,
          thrown: '',
          errorSummary: String(app.errorLog?.summary || ''),
          errorDetails: Array.isArray(app.errorLog?.details)
            ? app.errorLog.details.map((detail) => String(detail))
            : [],
          resultCount: Array.isArray(app.results) ? app.results.length : 0,
          generatedSvgCharacters: String(selected?.content || '').length,
          resultCommitted: String(selected?.content || '').includes('<svg'),
          previewVisible: Boolean(document.querySelector('.shadow-xl.origin-top > svg'))
        };
      } catch (error) {
        return {
          status: 'thrown',
          elapsedMs: performance.now() - started,
          thrown: `${String(error?.name || 'Error')}: ${String(error?.message || error)}`,
          errorSummary: String(app?.errorLog?.summary || ''),
          errorDetails: [],
          resultCount: Array.isArray(app?.results) ? app.results.length : 0,
          generatedSvgCharacters: 0,
          resultCommitted: false,
          previewVisible: Boolean(document.querySelector('.shadow-xl.origin-top > svg'))
        };
      }
    })
    .then((value) => ({ type: 'outcome', value }))
    .catch((error) => ({
      type: 'evaluation-error',
      value: `${String(error?.name || 'Error')}: ${String(error?.message || error)}`
    }));
  const completion = await Promise.race([
    evaluation,
    terminalSignal.then((value) => ({ type: 'terminal', value }))
  ]);
  if (completion.type === 'outcome') {
    outcome = completion.value;
  } else if (completion.type === 'evaluation-error') {
    evaluationError = completion.value;
  } else {
    evaluationError = `Browser terminal event: ${completion.value}`;
  }

  let worker = null;
  if (!page.isClosed() && !terminal.pageCrashed && !terminal.contextClosed) {
    try {
      worker = await getDiagramWorkerActivity(page);
    } catch (error) {
      evaluationError ||= `${String(error?.name || 'Error')}: ${String(error?.message || error)}`;
    }
  }
  const diagnostic = {
    outcome,
    evaluationError,
    terminal: { ...terminal },
    worker,
    lastLifecycle: runnerEvidence.findLast?.(({ type }) => type === 'lifecycle')?.value
      || [...runnerEvidence].reverse().find(({ type }) => type === 'lifecycle')?.value
      || null
  };
  return { ...diagnostic, classification: classifyGeneration(diagnostic) };
};

const assertRecoverableErrorState = async (page, diagnostic, originalPreview) => {
  if (diagnostic.classification !== 'returned-application-error') return null;
  expect(page.isClosed(), JSON.stringify(diagnostic, null, 2)).toBe(false);
  const recovery = await page.evaluate(async ({ expectedPreview }) => {
    const app = window.__GBDRAW_APP__;
    const selected = app.results?.[app.selectedResultIndex];
    const {
      DIAGRAM_HELPER_OPERATIONS,
      runDiagramHelperOperation
    } = await import('/gbdraw/web/js/services/diagram-generation.js');
    const helper = await runDiagramHelperOperation(
      DIAGRAM_HELPER_OPERATIONS.MEASURE_LEGEND_TEXT,
      { caption: 'Vibrio recovery probe', fontFamily: 'Arial', fontSize: 14 }
    );
    return {
      appAlive: Boolean(app),
      previewVisible: Boolean(document.querySelector('.shadow-xl.origin-top > svg')),
      previewPreserved: String(selected?.content || '') === expectedPreview,
      helperWidth: Number(helper?.result?.width || 0)
    };
  }, { expectedPreview: originalPreview });
  expect(recovery, JSON.stringify(diagnostic, null, 2)).toMatchObject({
    appAlive: true,
    previewVisible: true,
    previewPreserved: true
  });
  expect(recovery.helperWidth).toBeGreaterThan(0);
  return recovery;
};

test.describe.configure({ mode: 'serial' });

test('real Vibrio preview regenerates twice through staged binary resources', async ({
  page,
  context
}, testInfo) => {
  const terminal = {
    pageCrashed: false,
    pageClosed: false,
    contextClosed: false
  };
  let signalTerminal;
  const terminalSignal = new Promise((resolveSignal) => {
    signalTerminal = resolveSignal;
  });
  page.on('crash', () => {
    terminal.pageCrashed = true;
    signalTerminal('page-crash');
  });
  page.on('close', () => {
    terminal.pageClosed = true;
    signalTerminal('page-close');
  });
  context.on('close', () => {
    terminal.contextClosed = true;
    signalTerminal('context-close');
  });
  page.on('dialog', (dialog) => dialog.accept());
  const runnerEvidence = [];
  await installStructuralProbe(page, runnerEvidence);
  await openApp(page);
  await page.evaluate(() => window.__GBDRAW_VIBRIO_GENERATE_PROBE__.reset());

  await page.locator(
    'input[type="file"][accept*="application/json"][accept*="application/gzip"]'
  ).setInputFiles(vibrioSession);
  await page.waitForFunction(
    () => window.__GBDRAW_VIBRIO_GENERATE_PROBE__.lifecycle.some(
      ({ name }) => name === 'history-baseline-end'
    ),
    null,
    { timeout: 300_000 }
  );

  const loaded = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const selected = app.results?.[app.selectedResultIndex];
    return {
      originalPreview: String(selected?.content || ''),
      previewVisible: Boolean(document.querySelector('.shadow-xl.origin-top > svg')),
      resultCount: Array.isArray(app.results) ? app.results.length : 0
    };
  });
  const loadProbe = await probeSnapshot(page);
  const ready = loadProbe.lifecycle.find(({ name }) => name === 'interactiveReady');
  const previewMetrics = Object.fromEntries(
    PREVIEW_ZERO_METRICS.map((name) => [name, Number(loadProbe.metrics[name] || 0)])
  );
  expect(loaded.previewVisible).toBe(true);
  expect(loaded.resultCount).toBe(1);
  expect(ready).toMatchObject({ status: 'success', degradedRecovery: false });
  expect(previewMetrics).toEqual(Object.fromEntries(
    PREVIEW_ZERO_METRICS.map((name) => [name, 0])
  ));

  await page.evaluate(() => window.__GBDRAW_VIBRIO_GENERATE_PROBE__.reset());
  const first = await invokeGeneration(page, terminal, runnerEvidence, terminalSignal);
  const firstProbe = !page.isClosed() && !terminal.pageCrashed
    ? await probeSnapshot(page)
    : { metrics: {}, details: [], lifecycle: [] };
  const firstRecovery = await assertRecoverableErrorState(page, first, loaded.originalPreview);
  expect(first.classification, JSON.stringify(first, null, 2)).toBe('success');
  expect(first.outcome).toMatchObject({
    status: 'ok',
    resultCommitted: true,
    previewVisible: true
  });

  const second = await invokeGeneration(page, terminal, runnerEvidence, terminalSignal);
  const cumulativeProbe = !page.isClosed() && !terminal.pageCrashed
    ? await probeSnapshot(page)
    : { metrics: {}, details: [], lifecycle: [] };
  const secondRecovery = await assertRecoverableErrorState(page, second, loaded.originalPreview);
  expect(second.classification, JSON.stringify(second, null, 2)).toBe('success');
  expect(second.outcome).toMatchObject({
    status: 'ok',
    resultCommitted: true,
    previewVisible: true
  });

  const activity = second.worker;
  const runs = activity?.instances?.[0]?.runs || [];
  const lifecycleNames = cumulativeProbe.lifecycle.map(({ name }) => name);
  const resultSvgCharacters = cumulativeProbe.lifecycle
    .filter(({ name }) => name === 'result-svg-characters')
    .map(({ value }) => Number(value || 0));
  const structural = {
    canonicalReplayFullSerializationCount: Number(
      cumulativeProbe.metrics.canonicalReplayFullSerializationCount || 0
    ),
    workerBase64ResourceCloneCount: runs.filter(
      ({ hasBase64ResourceTable }) => hasBase64ResourceTable
    ).length,
    workerBase64ResourceJsonStringifyCount: lifecycleNames.filter(
      (name) => name.startsWith('worker-resources-json-')
    ).length,
    resourceMaterializationCount: Number(
      cumulativeProbe.metrics.resourceMaterializationCount || 0
    ),
    resourceReencodeCount: Number(cumulativeProbe.metrics.resourceReencodeCount || 0),
    resultBinaryDecodeCount: Number(cumulativeProbe.metrics.resultBinaryDecodeCount || 0),
    resultBinaryDecodedBytes: Number(cumulativeProbe.metrics.resultBinaryDecodedBytes || 0),
    resultMetadataBinaryDecodeCount: Number(
      cumulativeProbe.metrics.resultMetadataBinaryDecodeCount || 0
    ),
    resultMetadataBinaryDecodedBytes: Number(
      cumulativeProbe.metrics.resultMetadataBinaryDecodedBytes || 0
    )
  };
  const resultTransportEvents = cumulativeProbe.lifecycle.filter(
    ({ name }) => name === 'result-transport-ready'
  );

  expect(activity.constructions).toBe(1);
  expect(activity.initializations).toBe(1);
  expect(activity.runs).toBe(2);
  expect(activity.instances).toHaveLength(1);
  expect(activity.instances[0].terminated).toBe(false);
  expect(workerErrors(activity)).toEqual([]);
  expect(runs).toHaveLength(2);
  expect(runs[0].referencedResourceCount).toBeGreaterThan(0);
  expect(runs[0].referencedDeclaredBytes).toBeGreaterThan(0);
  expect(runs[0].stagedResourceCount).toBe(runs[0].referencedResourceCount);
  expect(runs[0].stagedResourceBytes).toBe(runs[0].referencedDeclaredBytes);
  expect(runs[0].transferredBytes).toBe(runs[0].stagedResourceBytes);
  expect(runs[1]).toMatchObject({
    referencedResourceCount: runs[0].referencedResourceCount,
    referencedDeclaredBytes: runs[0].referencedDeclaredBytes,
    stagedResourceCount: 0,
    stagedResourceBytes: 0,
    transferredBytes: 0,
    hasBase64ResourceTable: false
  });
  expect(structural).toMatchObject({
    canonicalReplayFullSerializationCount: 0,
    workerBase64ResourceCloneCount: 0,
    workerBase64ResourceJsonStringifyCount: 0,
    resourceMaterializationCount: runs[0].referencedResourceCount,
    resourceReencodeCount: 0,
    resultBinaryDecodeCount: 2,
    resultMetadataBinaryDecodeCount: 2
  });
  expect(structural.resultBinaryDecodedBytes).toBeGreaterThan(0);
  expect(structural.resultMetadataBinaryDecodedBytes).toBeGreaterThan(0);
  expect(resultTransportEvents).toHaveLength(2);
  expect(resultTransportEvents.every(
    ({ transport, bytes }) => transport === 'transferable-binary' && Number(bytes) > 0
  )).toBe(true);
  expect(resultSvgCharacters).toHaveLength(2);
  expect(resultSvgCharacters.every((characters) => characters > 0)).toBe(true);

  const workerPing = await page.evaluate(async () => {
    const {
      DIAGRAM_HELPER_OPERATIONS,
      runDiagramHelperOperation
    } = await import('/gbdraw/web/js/services/diagram-generation.js');
    const response = await runDiagramHelperOperation(
      DIAGRAM_HELPER_OPERATIONS.MEASURE_LEGEND_TEXT,
      { caption: 'Vibrio worker alive', fontFamily: 'Arial', fontSize: 14 }
    );
    return Number(response?.result?.width || 0);
  });
  expect(workerPing).toBeGreaterThan(0);
  const afterPing = await getDiagramWorkerActivity(page);
  expect(afterPing.constructions).toBe(1);
  expect(afterPing.initializations).toBe(1);
  expect(afterPing.instances[0].terminated).toBe(false);

  const report = {
    load: {
      previewMetrics,
      timings: loadTimings(loadProbe.lifecycle)
    },
    generations: {
      first,
      second,
      firstRecovery,
      secondRecovery,
      referencedResourceCount: runs[0].referencedResourceCount,
      referencedDeclaredBytes: runs[0].referencedDeclaredBytes,
      firstTransferredBytes: runs[0].transferredBytes,
      secondTransferredBytes: runs[1].transferredBytes,
      resourceMaterializations: structural.resourceMaterializationCount,
      resourceReencodes: structural.resourceReencodeCount,
      generatedSvgCharacters: resultSvgCharacters,
      firstElapsedMs: first.outcome.elapsedMs,
      secondElapsedMs: second.outcome.elapsedMs,
      firstProbe
    },
    structural,
    worker: afterPing
  };
  await testInfo.attach('vibrio-full-generation-metrics.json', {
    body: Buffer.from(JSON.stringify(report, null, 2)),
    contentType: 'application/json'
  });
  console.log(`Vibrio full-generation evidence: ${JSON.stringify({
    loadTimings: report.load.timings,
    firstGenerateMs: report.generations.firstElapsedMs,
    secondGenerateMs: report.generations.secondElapsedMs,
    referencedResourceCount: report.generations.referencedResourceCount,
    referencedDeclaredBytes: report.generations.referencedDeclaredBytes,
    firstTransferredBytes: report.generations.firstTransferredBytes,
    secondTransferredBytes: report.generations.secondTransferredBytes,
    structural: report.structural
  })}`);
});
