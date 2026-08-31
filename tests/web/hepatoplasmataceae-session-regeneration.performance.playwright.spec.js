const { test, expect } = require('@playwright/test');
const { createHash } = require('node:crypto');
const { existsSync, readFileSync, statSync } = require('node:fs');
const os = require('node:os');
const { join, resolve } = require('node:path');
const {
  getDiagramWorkerActivity,
  openApp
} = require('./helpers/app-lifecycle.cjs');
const {
  evaluateSessionRegenerationContract
} = require('./helpers/session-regeneration-contract.cjs');

const repoRoot = resolve(process.env.GBDRAW_REPO || process.cwd());
const fixturePath = join(
  repoRoot,
  'gbdraw',
  'web',
  'gallery',
  'sessions',
  'hepatoplasmataceae_collinear.gbdraw-session.json.gz'
);
const EXPECTED_FIXTURE_SHA = '2afba7111520b9dd7b00dffd351a1f4a21d4e0d9bfbaebc25ba7b5a937288577';
const EXPECTED_RECORD_IDS = [
  'AP027078.1',
  'AP027131.1',
  'AP027133.1',
  'AP027132.1',
  'NZ_CP006932.1'
];
const EXPECTED_RECORD_KEYS = ['record-1', 'record-2', 'record-3', 'record-4', 'record-5'];
const EXPECTED_RELATIONSHIPS = [
  'AP027078.1->AP027131.1',
  'AP027131.1->AP027133.1',
  'AP027133.1->AP027132.1',
  'AP027132.1->NZ_CP006932.1'
];
const OPERATION_TIMEOUT_MS = 900_000;
const INTERACTION_TIMEOUT_MS = 60_000;

const fixtureSha = createHash('sha256').update(readFileSync(fixturePath)).digest('hex');
const readRepositorySha = () => {
  const marker = join(repoRoot, '.git');
  const gitRoot = statSync(marker).isFile()
    ? resolve(repoRoot, readFileSync(marker, 'utf8').trim().slice('gitdir: '.length))
    : marker;
  const commonRoot = existsSync(join(gitRoot, 'commondir'))
    ? resolve(gitRoot, readFileSync(join(gitRoot, 'commondir'), 'utf8').trim())
    : gitRoot;
  const head = readFileSync(join(gitRoot, 'HEAD'), 'utf8').trim();
  if (!head.startsWith('ref: ')) return head;
  const ref = head.slice('ref: '.length);
  const looseRef = [gitRoot, commonRoot]
    .map((root) => join(root, ...ref.split('/')))
    .find((path) => existsSync(path));
  if (looseRef) return readFileSync(looseRef, 'utf8').trim();
  const packed = readFileSync(join(commonRoot, 'packed-refs'), 'utf8');
  const row = packed.split(/\r?\n/).find((line) => line.endsWith(` ${ref}`));
  if (!row) throw new Error(`Could not resolve repository ref ${ref}.`);
  return row.split(' ')[0];
};
const repositorySha = readRepositorySha();

const installRuntimeProbe = async (page) => page.addInitScript(() => {
  const metrics = Object.create(null);
  const details = [];
  const lifecycle = [];
  const responsiveness = {
    heartbeatCount: 0,
    maximumHeartbeatGapMs: 0,
    longTaskSupported: false,
    longTaskCount: 0,
    longTaskTotalDurationMs: 0,
    longestTaskMs: 0,
    memorySupported: Boolean(performance.memory),
    memoryHighWaterBytes: Number(performance.memory?.usedJSHeapSize || 0)
  };
  let lastHeartbeat = performance.now();

  const resetResponsiveness = () => {
    responsiveness.heartbeatCount = 0;
    responsiveness.maximumHeartbeatGapMs = 0;
    responsiveness.longTaskCount = 0;
    responsiveness.longTaskTotalDurationMs = 0;
    responsiveness.longestTaskMs = 0;
    responsiveness.memoryHighWaterBytes = Number(performance.memory?.usedJSHeapSize || 0);
    lastHeartbeat = performance.now();
  };

  window.__GBDRAW_EXACT_RUNTIME_PROBE__ = {
    phase: 'bootstrap',
    metrics,
    details,
    lifecycle,
    responsiveness,
    reset(phase) {
      this.phase = String(phase || 'phase');
      Object.keys(metrics).forEach((name) => delete metrics[name]);
      details.length = 0;
      lifecycle.length = 0;
      resetResponsiveness();
    },
    snapshot() {
      return {
        phase: this.phase,
        metrics: { ...metrics },
        details: details.map((entry) => ({ ...entry })),
        lifecycle: lifecycle.map((entry) => ({ ...entry })),
        responsiveness: { ...responsiveness }
      };
    }
  };

  window.__GBDRAW_TEST_HOOKS__ = {
    onStructuralMetric(metric) {
      const name = String(metric?.name || '');
      metrics[name] = Number(metrics[name] || 0) + Number(metric?.value || 0);
      details.push({ ...metric, timestamp: performance.now() });
    },
    onSessionLifecycleEvent(event) {
      lifecycle.push({ ...event });
    }
  };

  setInterval(() => {
    const now = performance.now();
    responsiveness.heartbeatCount += 1;
    responsiveness.maximumHeartbeatGapMs = Math.max(
      responsiveness.maximumHeartbeatGapMs,
      now - lastHeartbeat
    );
    lastHeartbeat = now;
    if (performance.memory) {
      responsiveness.memoryHighWaterBytes = Math.max(
        responsiveness.memoryHighWaterBytes,
        Number(performance.memory.usedJSHeapSize || 0)
      );
    }
  }, 100);

  try {
    const observer = new PerformanceObserver((list) => {
      for (const entry of list.getEntries()) {
        responsiveness.longTaskCount += 1;
        responsiveness.longTaskTotalDurationMs += Number(entry.duration || 0);
        responsiveness.longestTaskMs = Math.max(
          responsiveness.longestTaskMs,
          Number(entry.duration || 0)
        );
      }
    });
    observer.observe({ type: 'longtask', buffered: true });
    responsiveness.longTaskSupported = true;
  } catch (_error) {
    responsiveness.longTaskSupported = false;
  }
});

const probeSnapshot = (page) => page.evaluate(() => (
  window.__GBDRAW_EXACT_RUNTIME_PROBE__.snapshot()
));

const resetProbe = (page, phase) => page.evaluate((nextPhase) => {
  window.__GBDRAW_EXACT_RUNTIME_PROBE__.reset(nextPhase);
}, phase);

const workerTotals = (activity) => ({
  constructions: Number(activity.constructions || 0),
  initializations: Number(activity.initializations || 0),
  runs: Number(activity.runs || 0),
  stagedResourceBytes: (activity.instances || []).reduce(
    (total, instance) => total + (instance.runs || []).reduce(
      (runTotal, run) => runTotal + Number(run.stagedResourceBytes || 0),
      0
    ),
    0
  )
});

const sumMetricDetails = (details, name, predicate = () => true) => (
  details
    .filter((detail) => detail.name === name && predicate(detail))
    .reduce((total, detail) => total + Number(detail.value || 0), 0)
);

const lifecycleCount = (events, name) => events.filter((event) => event.name === name).length;

const eventTimestamp = (events, name) => (
  events.find((event) => event.name === name)?.timestamp ?? null
);

const durationBetween = (events, startName, endName) => {
  const start = eventTimestamp(events, startName);
  const end = eventTimestamp(events, endName);
  return start === null || end === null ? null : end - start;
};

const phaseDiagnostics = (snapshot, wallDurationMs) => ({
  totalDurationMs: wallDurationMs,
  workerCompleteToAdmissionCompleteMs: durationBetween(
    snapshot.lifecycle,
    'worker-artifact-identity-end',
    'svg.admission-completed'
  ),
  admissionCompleteToMountObservedMs: durationBetween(
    snapshot.lifecycle,
    'svg.admission-completed',
    'preview.mount-observed'
  ),
  mountObservedToBindCompleteMs: durationBetween(
    snapshot.lifecycle,
    'preview.mount-observed',
    'preview.bind-completed'
  ),
  bindCompleteToPostBindFrameMs: durationBetween(
    snapshot.lifecycle,
    'preview.bind-completed',
    'preview.post-bind-frame-completed'
  ),
  responsiveness: snapshot.responsiveness,
  lastLifecycleEvent: snapshot.lifecycle.at(-1)?.name || ''
});

const deriveContractMetrics = (snapshot, beforeWorker, afterWorker) => {
  const before = workerTotals(beforeWorker);
  const after = workerTotals(afterWorker);
  const readyTimestamp = eventTimestamp(snapshot.lifecycle, 'preview.ready-receipt-accepted');
  const beforeReady = (detail) => (
    readyTimestamp === null || Number(detail.timestamp || 0) < readyTimestamp
  );
  return {
    ...snapshot.metrics,
    rendererExecutionCount: lifecycleCount(snapshot.lifecycle, 'python-wrapper-start'),
    admissionFeatureDomFullScanCount: sumMetricDetails(
      snapshot.details,
      'featureDomFullScanCount',
      (detail) => detail.phase === 'current-worker'
    ),
    admissionLegendDomFullScanCount: sumMetricDetails(
      snapshot.details,
      'legendDomFullScanCount',
      (detail) => detail.phase === 'current-worker'
    ),
    preReadyFeatureDomFullScanCount: sumMetricDetails(
      snapshot.details,
      'featureDomFullScanCount',
      beforeReady
    ),
    preReadyLegendDomFullScanCount: sumMetricDetails(
      snapshot.details,
      'legendDomFullScanCount',
      beforeReady
    ),
    preReadyFeatureSearchIndexBuildCount: sumMetricDetails(
      snapshot.details,
      'featureSearchIndexBuildCount',
      beforeReady
    ),
    workerConstructionCount: after.constructions - before.constructions,
    workerInitializationCount: after.initializations - before.initializations,
    resourceTransferredByteCount: after.stagedResourceBytes - before.stagedResourceBytes
  };
};

const loadSessionThroughUi = async (page, dialogMessages) => {
  await resetProbe(page, 'session-load');
  const chooserPromise = page.waitForEvent('filechooser', { timeout: OPERATION_TIMEOUT_MS });
  await page.getByRole('button', { name: 'Load Session', exact: true }).click();
  const chooser = await chooserPromise;
  await chooser.setFiles(fixturePath);
  await page.waitForFunction(() => (
    window.__GBDRAW_EXACT_RUNTIME_PROBE__?.lifecycle?.some(
      (event) => event.name === 'interactiveReady'
    )
  ), null, { timeout: OPERATION_TIMEOUT_MS });
  await page.evaluate(() => new Promise((resolveFrame) => (
    requestAnimationFrame(() => requestAnimationFrame(resolveFrame))
  )));
  expect(dialogMessages).toContain('Session loaded successfully!');
};

const inspectSavedPreview = (page) => page.evaluate(async () => {
  const { state } = await import('/gbdraw/web/js/state.js');
  const { getCommittedCanonicalSession } = await import(
    '/gbdraw/web/js/services/config.js'
  );
  const app = window.__GBDRAW_APP__;
  const catalog = state.featureCatalog.value;
  const item = catalog?.items?.[0];
  const committed = getCommittedCanonicalSession();
  const precomputed = (committed?.renderRequest?.comparisons || []).filter(
    (comparison) => comparison.kind === 'precomputedProteinComparison'
  );
  return {
    resultCount: app.results.length,
    selectedResultIndex: app.selectedResultIndex,
    rootMounted: Boolean(state.svgContainer.value?.querySelector?.('svg')),
    catalogSchema: catalog?.schema,
    catalogItemCount: catalog?.items?.length || 0,
    catalogResultIndex: item?.resultIndex,
    activePlan: {
      mode: app.linearComparisonPlan.mode,
      defaultSource: app.linearComparisonPlan.defaultSource,
      edgeCount: app.linearComparisonPlan.edges.length
    },
    committedRequest: {
      schema: committed?.renderRequest?.schema,
      precomputedComparisonCount: precomputed.length,
      pairs: precomputed.map((comparison) => [
        comparison.queryRecordIndex,
        comparison.subjectRecordIndex
      ])
    }
  };
});

const exerciseFeatureInteraction = async (page) => {
  const feature = page.locator(
    '.shadow-xl.origin-top > svg '
      + '[data-gbdraw-feature-id]:not([data-gbdraw-auto-feature-underlay="true"])'
  ).first();
  await expect(feature).toBeVisible({ timeout: INTERACTION_TIMEOUT_MS });
  await feature.click();
  const dialog = page.getByRole('dialog', { name: /Feature details:/ });
  await expect(dialog).toBeVisible({ timeout: INTERACTION_TIMEOUT_MS });
  await dialog.getByRole('button', { name: 'Close feature popup' }).click();
  await expect(dialog).toHaveCount(0);
};

const semanticSnapshot = (page) => page.evaluate(async () => {
  const { state } = await import('/gbdraw/web/js/state.js');
  const app = window.__GBDRAW_APP__;
  const catalog = state.featureCatalog.value;
  const item = catalog?.items?.[app.selectedResultIndex];
  const svg = state.svgContainer.value?.querySelector?.('svg');
  if (!svg || !item) throw new Error('The selected generated Result is not mounted with a catalog.');

  const ids = [...svg.querySelectorAll('[id]')].map((element) => element.id);
  const unsafeElements = svg.querySelectorAll('script, foreignObject, iframe, object, embed').length;
  let unsafeAttributes = 0;
  for (const element of svg.querySelectorAll('*')) {
    for (const attribute of element.attributes) {
      const name = attribute.name.toLowerCase();
      const value = attribute.value.trim().toLowerCase();
      if (name.startsWith('on') || value.startsWith('javascript:')) unsafeAttributes += 1;
    }
  }
  const relationships = [...new Set((item.comparisonMatches || []).map((match) => (
    `${match.query_record_id}->${match.subject_record_id}`
  )))];
  const comparisonGroups = [...svg.querySelectorAll(
    'g[data-query-record-index][data-subject-record-index]'
  )].map((group) => (
    `${group.getAttribute('data-query-record-index')}->${group.getAttribute('data-subject-record-index')}`
  ));
  return {
    resultCount: app.results.length,
    selectedResultIndex: app.selectedResultIndex,
    catalogSchema: catalog.schema,
    catalogItemCount: catalog.items.length,
    catalogResultIndex: item.resultIndex,
    recordKeys: item.recordKeys,
    recordIds: item.sequenceSources.map((source) => source.recordId),
    renderedFeatureCount: item.features.length,
    biologicalFeatureCount: item.biologicalFeatures.length,
    orthogroupCount: item.orthogroups.length,
    comparisonMatchCount: item.comparisonMatches.length,
    relationships,
    svgIdCount: ids.length,
    uniqueSvgIdCount: new Set(ids).size,
    featureHookCount: svg.querySelectorAll('[data-gbdraw-feature-id]').length,
    stableFeatureHookCount: svg.querySelectorAll('[data-gbdraw-stable-feature-id]').length,
    matchHookCount: svg.querySelectorAll('[data-gbdraw-match-id]').length,
    pairwiseMatchHookCount: svg.querySelectorAll('[data-gbdraw-pairwise-match-id]').length,
    recordDefinitionCount: svg.querySelectorAll('[data-gbdraw-role="record-definition"]').length,
    comparisonGroups,
    comparisonLegendCount: svg.querySelectorAll('[data-gbdraw-role="comparison-legend"]').length,
    legendGroupCount: svg.querySelectorAll('[data-gbdraw-composition-role="legend"]').length,
    unsafeElements,
    unsafeAttributes
  };
});

const assertSemanticSnapshot = (snapshot) => {
  expect(snapshot).toMatchObject({
    resultCount: 1,
    selectedResultIndex: 0,
    catalogSchema: 3,
    catalogItemCount: 1,
    catalogResultIndex: 0,
    recordKeys: EXPECTED_RECORD_KEYS,
    recordIds: EXPECTED_RECORD_IDS,
    renderedFeatureCount: 2_994,
    biologicalFeatureCount: 5_982,
    orthogroupCount: 554,
    comparisonMatchCount: 500,
    relationships: EXPECTED_RELATIONSHIPS,
    featureHookCount: 2_994,
    stableFeatureHookCount: 2_994,
    matchHookCount: 500,
    pairwiseMatchHookCount: 500,
    recordDefinitionCount: 5,
    comparisonGroups: ['0->1', '1->2', '2->3', '3->4'],
    comparisonLegendCount: 2,
    legendGroupCount: 1,
    unsafeElements: 0,
    unsafeAttributes: 0
  });
  expect(snapshot.svgIdCount).toBeGreaterThan(0);
  expect(snapshot.uniqueSvgIdCount).toBe(snapshot.svgIdCount);
};

const exerciseReadyInteractions = async (page) => {
  const beforeSearch = await probeSnapshot(page);
  expect(Number(beforeSearch.metrics.featureSearchIndexBuildCount || 0)).toBe(0);

  await page.getByRole('searchbox', { name: 'Search features' }).fill('dnaA');
  await page.getByRole('button', { name: 'Search features', exact: true }).click();
  await expect(page.getByRole('button', { name: 'Open active feature' })).toBeEnabled({
    timeout: INTERACTION_TIMEOUT_MS
  });
  await page.getByRole('button', { name: 'Search features', exact: true }).click();
  const afterRepeatedSearch = await probeSnapshot(page);
  expect(Number(afterRepeatedSearch.metrics.featureSearchIndexBuildCount || 0)).toBe(1);

  await page.getByRole('button', { name: 'Open active feature' }).click();
  const featureDialog = page.getByRole('dialog', { name: /Feature details:/ });
  await expect(featureDialog).toBeVisible();
  await featureDialog.getByRole('button', { name: 'Close feature popup' }).click();

  const match = page.locator(
    '.shadow-xl.origin-top > svg [data-gbdraw-pairwise-match-id]'
  ).first();
  await expect(match).toHaveAttribute('role', 'button');
  await match.focus();
  await match.press('Enter');
  const matchDialog = page.getByRole('dialog', { name: 'Pairwise match details' });
  await expect(matchDialog).toBeVisible();
  await matchDialog.getByRole('button', { name: 'Close match popup' }).click();

  await expect(page.getByRole('button', { name: 'SVG', exact: true })).toBeEnabled();
  await expect(page.getByRole('button', { name: 'Save Session', exact: true })).toBeEnabled();
  await page.getByRole('button', { name: 'Clear search', exact: true }).click();
};

const generateThroughUi = async ({
  page,
  profile,
  priorWorker,
  priorBindSequence
}) => {
  await resetProbe(page, profile);
  const startedAt = Date.now();
  await page.getByRole('button', { name: 'Generate Diagram', exact: true }).click();
  await page.waitForFunction(() => (
    window.__GBDRAW_EXACT_RUNTIME_PROBE__?.lifecycle?.some(
      (event) => event.name === 'generate.processing-cleared'
    )
  ), null, { timeout: OPERATION_TIMEOUT_MS });
  const wallDurationMs = Date.now() - startedAt;
  const snapshot = await probeSnapshot(page);
  const worker = await getDiagramWorkerActivity(page);
  const cleared = snapshot.lifecycle.find((event) => event.name === 'generate.processing-cleared');
  expect(cleared?.status).toBe('ok');
  expect(lifecycleCount(snapshot.lifecycle, 'generate.completed')).toBe(1);
  expect(lifecycleCount(snapshot.lifecycle, 'python-wrapper-end')).toBe(1);

  const readyReceipt = snapshot.lifecycle.find(
    (event) => event.name === 'preview.ready-receipt-accepted'
  );
  const rootGeneration = readyReceipt?.rootGeneration;
  const bindSequence = readyReceipt?.bindSequence;
  expect(Number.isSafeInteger(rootGeneration)).toBe(true);
  expect(Number.isSafeInteger(bindSequence)).toBe(true);
  if (priorBindSequence !== null) expect(bindSequence).not.toBe(priorBindSequence);

  const semantic = await semanticSnapshot(page);
  assertSemanticSnapshot(semantic);
  await exerciseReadyInteractions(page);

  const metrics = deriveContractMetrics(snapshot, priorWorker, worker);
  const failures = evaluateSessionRegenerationContract({
    profile,
    metrics,
    events: snapshot.lifecycle,
    resultCount: semantic.resultCount,
    state: { terminalOutcome: 'success', interactiveProbePassed: true }
  });
  expect(failures, JSON.stringify({ failures, metrics, lifecycle: snapshot.lifecycle }, null, 2))
    .toEqual([]);

  const workerErrors = (worker.instances || []).flatMap((instance) => instance.errors || []);
  expect(workerErrors).toEqual([]);
  return {
    snapshot,
    worker,
    rootGeneration,
    bindSequence,
    semantic,
    metrics,
    diagnostics: phaseDiagnostics(snapshot, wallDurationMs)
  };
};

test.describe.configure({ mode: 'serial' });

test('exact saved Session regenerates twice with bounded work and fresh readiness', async ({
  page,
  context,
  browser
}, testInfo) => {
  test.setTimeout(1_800_000);
  page.setDefaultTimeout(INTERACTION_TIMEOUT_MS);
  expect(fixtureSha).toBe(EXPECTED_FIXTURE_SHA);

  const terminal = {
    pageCrashes: [],
    pageClosures: 0,
    contextClosures: 0,
    pageErrors: []
  };
  const dialogMessages = [];
  const evidence = {
    repositorySha,
    fixtureSha,
    browser: {
      project: testInfo.project.name,
      version: browser.version(),
      platform: process.platform,
      node: process.version,
      osRelease: os.release()
    },
    terminal
  };

  page.on('crash', () => terminal.pageCrashes.push('page crash'));
  page.on('close', () => { terminal.pageClosures += 1; });
  context.on('close', () => { terminal.contextClosures += 1; });
  page.on('pageerror', (error) => terminal.pageErrors.push(String(error?.message || error)));
  page.on('dialog', async (dialog) => {
    dialogMessages.push(dialog.message());
    await dialog.accept();
  });

  await installRuntimeProbe(page);
  try {
    await openApp(page);
    await loadSessionThroughUi(page, dialogMessages);
    const savedProbe = await probeSnapshot(page);
    const savedWorker = await getDiagramWorkerActivity(page);
    const saved = await inspectSavedPreview(page);
    expect(saved).toEqual({
      resultCount: 1,
      selectedResultIndex: 0,
      rootMounted: true,
      catalogSchema: 3,
      catalogItemCount: 1,
      catalogResultIndex: 0,
      activePlan: { mode: 'adjacent', defaultSource: 'losat', edgeCount: 0 },
      committedRequest: {
        schema: 5,
        precomputedComparisonCount: 4,
        pairs: [[0, 1], [1, 2], [2, 3], [3, 4]]
      }
    });
    await exerciseFeatureInteraction(page);
    const savedFailures = evaluateSessionRegenerationContract({
      profile: 'saved-preview',
      metrics: {
        ...savedProbe.metrics,
        workerConstructionCount: savedWorker.constructions,
        workerInitializationCount: savedWorker.initializations
      },
      state: {
        selectedResultCount: saved.resultCount,
        currentCatalog: saved.catalogSchema === 3 && saved.catalogResultIndex === 0,
        activeDraftCommittedDistinct: saved.activePlan.edgeCount === 0
          && saved.committedRequest.precomputedComparisonCount === 4,
        interactiveProbePassed: true
      }
    });
    expect(savedFailures).toEqual([]);

    const first = await generateThroughUi({
      page,
      profile: 'first-generate',
      priorWorker: savedWorker,
      priorBindSequence: null
    });
    const second = await generateThroughUi({
      page,
      profile: 'second-generate',
      priorWorker: first.worker,
      priorBindSequence: first.bindSequence
    });

    evidence.savedPreview = {
      lifecycleLastEvent: savedProbe.lifecycle.at(-1)?.name || '',
      worker: workerTotals(savedWorker),
      activeDraftCommittedDistinct: true
    };
    evidence.firstGenerate = {
      readiness: {
        rootGeneration: first.rootGeneration,
        bindSequence: first.bindSequence
      },
      metrics: first.metrics,
      semantic: first.semantic,
      diagnostics: first.diagnostics
    };
    evidence.secondGenerate = {
      readiness: {
        rootGeneration: second.rootGeneration,
        bindSequence: second.bindSequence
      },
      metrics: second.metrics,
      semantic: second.semantic,
      diagnostics: second.diagnostics
    };
    evidence.workerErrors = (second.worker.instances || []).flatMap(
      (instance) => instance.errors || []
    );
    console.log(`GBDRAW_EXACT_EVIDENCE ${JSON.stringify(evidence)}`);

    expect(terminal).toEqual({
      pageCrashes: [],
      pageClosures: 0,
      contextClosures: 0,
      pageErrors: []
    });
    expect(evidence.workerErrors).toEqual([]);
  } finally {
    if (!page.isClosed()) {
      try {
        evidence.lastProbe = await probeSnapshot(page);
        evidence.lastWorker = workerTotals(await getDiagramWorkerActivity(page));
      } catch (error) {
        evidence.diagnosticCollectionError = String(error?.message || error);
      }
    }
    await testInfo.attach('hepatoplasmataceae-session-regeneration-evidence.json', {
      body: Buffer.from(JSON.stringify(evidence, null, 2)),
      contentType: 'application/json'
    });
  }
});
