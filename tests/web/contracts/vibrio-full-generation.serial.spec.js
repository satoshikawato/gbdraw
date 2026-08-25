const { test, expect } = require('@playwright/test');
const { spawnSync } = require('node:child_process');
const { createHash } = require('node:crypto');
const { writeFileSync } = require('node:fs');
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

const featureCatalogDigest = (page) => page.evaluate(async () => {
  const { state } = await import('/gbdraw/web/js/state.js');
  const catalog = state.featureCatalog?.value ?? null;
  const encoded = new TextEncoder().encode(JSON.stringify(catalog));
  const digest = new Uint8Array(await crypto.subtle.digest('SHA-256', encoded));
  return {
    sha256: [...digest].map((value) => value.toString(16).padStart(2, '0')).join(''),
    jsonBytes: encoded.byteLength,
    schema: Number(catalog?.schema || 0),
    itemCount: Array.isArray(catalog?.items) ? catalog.items.length : 0
  };
});

const featureCatalogSummary = (page) => page.evaluate(async () => {
  const { state } = await import('/gbdraw/web/js/state.js');
  const catalog = state.featureCatalog?.value ?? null;
  return {
    schema: Number(catalog?.schema || 0),
    itemCount: Array.isArray(catalog?.items) ? catalog.items.length : 0
  };
});

const svgSourceIdentity = (content) => ({
  characters: String(content || '').length,
  sha256: createHash('sha256').update(String(content || '')).digest('hex')
});

const svgComparisonSummary = (content) => {
  const source = String(content || '');
  return {
    comparisonGroups: (source.match(/<g[^>]*\bid="comparison\d+"/g) || []).length,
    pairwiseMatches: (source.match(/data-gbdraw-pairwise-match-id=/g) || []).length,
    comparisonLegends: (source.match(/data-gbdraw-role="comparison-legend"/g) || []).length
  };
};

const svgLinearGeometrySummary = (content) => {
  const groupTags = String(content || '').match(/<g\b[^>]*>/g) || [];
  const attribute = (tag, name) => (
    tag.match(new RegExp(`\\b${name}="([^"]*)"`))?.[1] ?? null
  );
  const placements = groupTags
    .filter((tag) => (
      attribute(tag, 'data-record-index') !== null
      && attribute(tag, 'data-record-row') !== null
      && attribute(tag, 'data-record-column') !== null
    ))
    .map((tag) => ({
      record: Number(attribute(tag, 'data-record-index')),
      row: Number(attribute(tag, 'data-record-row')),
      column: Number(attribute(tag, 'data-record-column'))
    }))
    .sort((left, right) => left.record - right.record);
  const rowNumbers = [...new Set(placements.map(({ row }) => row))]
    .sort((left, right) => left - right);
  const recordsByRow = rowNumbers.map((row) => placements
    .filter((placement) => placement.row === row)
    .sort((left, right) => left.column - right.column)
    .map(({ record }) => record));
  const definitionGroups = groupTags.filter((tag) => (
    ['record-definition', 'record-definition-row'].includes(
      attribute(tag, 'data-gbdraw-role')
    )
  ));
  return {
    biologicalGroups: placements.length,
    rows: rowNumbers.length,
    rowSizes: recordsByRow.map((records) => records.length),
    placements,
    recordsByRow,
    definitionPseudoPlacements: definitionGroups.filter((tag) => (
      attribute(tag, 'data-record-row') !== null
      || attribute(tag, 'data-record-column') !== null
    )).length
  };
};

const EXPECTED_VIBRIO_LINEAR_GEOMETRY = {
  biologicalGroups: 11,
  rows: 5,
  rowSizes: [3, 2, 2, 2, 2],
  placements: [
    { record: 0, row: 0, column: 0 },
    { record: 1, row: 0, column: 1 },
    { record: 2, row: 0, column: 2 },
    { record: 3, row: 1, column: 0 },
    { record: 4, row: 1, column: 1 },
    { record: 5, row: 2, column: 0 },
    { record: 6, row: 2, column: 1 },
    { record: 7, row: 3, column: 0 },
    { record: 8, row: 3, column: 1 },
    { record: 9, row: 4, column: 0 },
    { record: 10, row: 4, column: 1 }
  ],
  recordsByRow: [
    [0, 1, 2],
    [3, 4],
    [5, 6],
    [7, 8],
    [9, 10]
  ],
  definitionPseudoPlacements: 0
};

const installCanonicalRequestCapture = (page) => page.evaluate(() => {
  const previousPostMessage = Worker.prototype.postMessage;
  window.__GBDRAW_VIBRIO_CANONICAL_REQUESTS__ = [];
  Worker.prototype.postMessage = function captureCanonicalRequest(message, transfer) {
    if (message?.type === 'run' && message.payload?.request) {
      window.__GBDRAW_VIBRIO_CANONICAL_REQUESTS__.push(
        {
          request: structuredClone(message.payload.request),
          resourceManifest: structuredClone(message.payload.resourceManifest || [])
        }
      );
    }
    if (transfer === undefined) return previousPostMessage.call(this, message);
    return previousPostMessage.call(this, message, transfer);
  };
});

const canonicalRequestEvidence = async (page) => {
  const capture = await page.evaluate(() => {
    const requests = [...(window.__GBDRAW_VIBRIO_CANONICAL_REQUESTS__ || [])];
    return {
      requests: requests.map(({ request }) => JSON.stringify(request)),
      resourceManifests: requests.map(
        ({ resourceManifest }) => JSON.stringify(resourceManifest)
      )
    };
  });
  const requestDigests = capture.requests.map(
    (serialized) => createHash('sha256').update(serialized).digest('hex')
  );
  const resourceManifestDigests = capture.resourceManifests.map(
    (serialized) => createHash('sha256').update(serialized).digest('hex')
  );
  return {
    requests: capture.requests.map((serialized, index) => {
      const request = JSON.parse(serialized);
      const resourceManifest = JSON.parse(capture.resourceManifests[index]);
      return {
        sha256: requestDigests[index],
        schema: Number(request.schema || 0),
        mode: String(request.mode || ''),
        recordCount: Array.isArray(request.records) ? request.records.length : 0,
        comparisonCount: Array.isArray(request.comparisons) ? request.comparisons.length : 0,
        comparisonKinds: Array.isArray(request.comparisons)
          ? request.comparisons.map(({ kind }) => String(kind || ''))
          : [],
        layout: request.layout || null,
        colors: request.diagramOptions?.colors || null,
        resourceManifestSha256: resourceManifestDigests[index],
        referencedResourceCount: resourceManifest.length,
        referencedResourceBytes: resourceManifest.reduce(
          (total, { size }) => total + Number(size || 0),
          0
        )
      };
    }),
    repeatComparison: capture.requests.length === 2
      ? {
          equivalent: requestDigests[0] === requestDigests[1]
            && resourceManifestDigests[0] === resourceManifestDigests[1]
        }
      : null
  };
};

const activeIntentSummary = (page) => page.evaluate(async () => {
  const { state } = await import('/gbdraw/web/js/state.js');
  const stable = (value) => {
    if (Array.isArray(value)) return value.map((item) => stable(item));
    if (!value || typeof value !== 'object') return value;
    return Object.fromEntries(
      Object.keys(value).sort().map((key) => [key, stable(value[key])])
    );
  };
  const plain = (value) => JSON.parse(JSON.stringify(value ?? null));
  const overrides = {
    fills: plain(state.featureColorOverrides),
    strokes: plain(state.featureStrokeOverrides),
    visibility: plain(state.featureVisibilityOverrides),
    labelText: plain(state.labelTextFeatureOverrides),
    labelVisibility: plain(state.labelVisibilityOverrides),
    legendColors: plain(state.legendColorOverrides),
    legendStrokes: plain(state.legendStrokeOverrides)
  };
  const digest = async (value) => {
    const bytes = new TextEncoder().encode(JSON.stringify(stable(value)));
    const hash = new Uint8Array(await crypto.subtle.digest('SHA-256', bytes));
    return [...hash].map((item) => item.toString(16).padStart(2, '0')).join('');
  };
  const summary = {
    mode: state.mode.value,
    inputType: state.mode.value === 'circular'
      ? state.cInputType.value
      : state.lInputType.value,
    palette: {
      selected: state.selectedPalette.value,
      currentColorsSha256: await digest(state.currentColors.value),
      currentColorCount: Object.keys(state.currentColors.value || {}).length,
      instantPreview: state.paletteInstantPreviewEnabled.value,
      applied: state.appliedPaletteName.value,
      pending: state.pendingPaletteName.value
    },
    rules: {
      specificCount: state.manualSpecificRules.length,
      specificSha256: await digest(state.manualSpecificRules),
      priorityCount: state.manualPriorityRules.length,
      prioritySha256: await digest(state.manualPriorityRules)
    },
    filters: {
      mode: state.filterMode.value,
      whitelistCount: state.manualWhitelist.length,
      blacklistText: state.manualBlacklist.value
    },
    form: {
      plotTitle: state.form.plot_title,
      labelsMode: state.form.labels_mode,
      showScale: state.form.show_scale,
      legend: state.form.legend
    },
    adv: {
      axisStrokeWidth: state.adv.axis_stroke_width,
      labelFontSize: state.adv.label_font_size,
      featureWidthCircular: state.adv.feature_width_circular,
      plotTitlePosition: state.adv.plot_title_position
    },
    annotationSetIds: state.annotationSets.map(({ id }) => id),
    circularTracks: {
      enabled: state.adv.circular_track_slots_enabled,
      slots: state.adv.circular_track_slots.map(({ id, renderer, enabled }) => ({
        id,
        renderer,
        enabled
      }))
    },
    linearComparisonPlan: plain(state.linearComparisonPlan),
    linearRecordLayout: {
      enabled: state.linearRecordLayoutEnabled.value,
      recordGap: state.linearRecordGap.value,
      rows: plain(state.linearRecordRows)
    },
    layoutPreferences: plain(state.layoutPreferences),
    editorOverrides: {
      counts: Object.fromEntries(
        Object.entries(overrides).map(([name, values]) => [name, Object.keys(values || {}).length])
      ),
      sha256: await digest(overrides)
    }
  };
  return { ...summary, sha256: await digest(summary) };
});

const phaseDuration = (events, startName, endName) => {
  const start = events.find(({ name }) => name === startName);
  const end = events.find(({ name }) => name === endName);
  return start && end ? end.timestamp - start.timestamp : null;
};

const durationOrZero = (events, startName, endName) => (
  phaseDuration(events, startName, endName) ?? 0
);

const HISTORY_DIAGNOSTIC_FIELDS = [
  'artifactCheckpointBuilds',
  'artifactCheckpointSignatureComputations',
  'intentBuilds',
  'intentSignatureComputations',
  'historySvgBytes',
  'checkpointEstimatedBytes',
  'generatedArtifactFullCloneCount',
  'generatedArtifactFullSerializationCount',
  'manualCancelFullArtifactSnapshotBuildCount',
  'artifactHandleBeforeBuildCount',
  'artifactHandleAfterBuildCount',
  'artifactFingerprintComparisonCount',
  'artifactReplacementHistoryEntryCount'
];

const historyStructuralDelta = (events, startName, endName) => {
  const start = events.find(({ name }) => name === startName)?.diagnostics || {};
  const end = events.find(({ name }) => name === endName)?.diagnostics || {};
  return Object.fromEntries(HISTORY_DIAGNOSTIC_FIELDS.map((name) => [
    name,
    Number(end[name] || 0) - Number(start[name] || 0)
  ]));
};

const generationPhaseAttribution = (outcome, probe) => {
  const events = probe.lifecycle || [];
  const python = events.find(({ name }) => name === 'python-diagnostics') || {};
  const pythonTimings = python.timingsMs || {};
  const pythonMetrics = python.metrics || {};
  const initialization = events.find(
    ({ name }) => name === 'worker-initialization-resolved'
  ) || {};
  const resourceLinking = events.find(
    ({ name }) => name === 'worker-resource-linking-end'
  ) || {};
  const losatCache = events.find(
    ({ name }) => name === 'losat-cache-preparation-end'
  ) || {};
  const resultTransportDecodeMs = (
    durationOrZero(events, 'result-binary-decode-start', 'result-binary-decode-end')
    + durationOrZero(
      events,
      'result-metadata-binary-decode-start',
      'result-metadata-binary-decode-end'
    )
  );
  const inputResolutionMs = durationOrZero(
    events,
    'generation-input-resolution-start',
    'generation-input-resolution-end'
  );
  const canonicalFileSerializationMs = durationOrZero(
    events,
    'serialize-canonical-files-start',
    'serialize-canonical-files-end'
  );
  const historyBeforeStructural = historyStructuralDelta(
    events,
    'generate-history-before-start',
    'generate-history-before-end'
  );
  const historyAfterStructural = historyStructuralDelta(
    events,
    'generate-history-after-start',
    'generate-history-after-end'
  );
  return {
    warmGenerateMs: Number(outcome?.elapsedMs || 0),
    historyBeforeMs: durationOrZero(
      events,
      'generate-history-before-start',
      'generate-history-before-end'
    ),
    jsPreparationMs: inputResolutionMs + canonicalFileSerializationMs,
    inputResolutionMs,
    losatCachePreparationMs: durationOrZero(
      events,
      'losat-cache-preparation-start',
      'losat-cache-preparation-end'
    ),
    losatCacheHits: Number(losatCache.cacheHits || 0),
    losatCacheMisses: Number(losatCache.cacheMisses || 0),
    canonicalFileSerializationMs,
    canonicalRequestMs: durationOrZero(
      events,
      'canonical-request-construction-start',
      'canonical-request-construction-end'
    ),
    resourcePrepareMs: durationOrZero(
      events,
      'resource-stage-start',
      'resource-stage-end'
    ),
    workerInitializationMs: Number(initialization.durationMs || 0),
    workerInitializationReused: initialization.reused === true,
    workerWorkspaceMs: durationOrZero(
      events,
      'worker-workspace-preparation-start',
      'worker-workspace-preparation-end'
    ),
    workerResourceLinkingMs: durationOrZero(
      events,
      'worker-resource-linking-start',
      'worker-resource-linking-end'
    ),
    newlyStagedResourceBytes: Number(resourceLinking.newlyStagedResourceBytes || 0),
    pythonInvocationMs: durationOrZero(
      events,
      'python-wrapper-start',
      'python-wrapper-end'
    ),
    pythonDecodeMs: Number(pythonTimings.decode || 0),
    pythonArtifactCopyMs: Number(pythonTimings.artifactCopy || 0),
    pythonArtifactValidationMs: Number(pythonTimings.artifactValidation || 0),
    pythonRecordLoadMs: Number(pythonTimings.recordLoad || 0),
    pythonPreparationMs: Number(pythonTimings.preparation || 0),
    pythonComparisonPreparationMs: Number(
      pythonTimings.comparisonPreparation || 0
    ),
    pythonDrawingMs: Number(pythonTimings.drawing || 0),
    pythonInteractivePreparationMs: Number(
      pythonTimings.interactivePreparation || 0
    ),
    pythonSvgWriteMs: Number(pythonTimings.svgWrite || 0),
    pythonSvgReadbackMs: Number(pythonTimings.svgReadback || 0),
    pythonFeatureCatalogMs: Number(pythonTimings.featureCatalog || 0),
    pythonGeometryMetadataMs: Number(pythonTimings.geometryMetadata || 0),
    pythonArtifactMetrics: { ...pythonMetrics },
    preparedInputCacheStructural: {
      decodedResourceCacheHitCount: Number(
        pythonMetrics.decodedResourceCacheHitCount || 0
      ),
      decodedResourceCacheMissCount: Number(
        pythonMetrics.decodedResourceCacheMissCount || 0
      ),
      decodedResourceBuildCount: Number(
        pythonMetrics.decodedResourceBuildCount || 0
      ),
      parsedSourceCacheHitCount: Number(
        pythonMetrics.parsedSourceCacheHitCount || 0
      ),
      parsedSourceCacheMissCount: Number(
        pythonMetrics.parsedSourceCacheMissCount || 0
      ),
      parsedSourceParseCount: Number(pythonMetrics.parsedSourceParseCount || 0),
      resolvedRecordCacheHitCount: Number(
        pythonMetrics.resolvedRecordCacheHitCount || 0
      ),
      resolvedRecordCacheMissCount: Number(
        pythonMetrics.resolvedRecordCacheMissCount || 0
      ),
      resolvedRecordBuildCount: Number(
        pythonMetrics.resolvedRecordBuildCount || 0
      ),
      interactiveContextCacheHitCount: Number(
        pythonMetrics.interactiveContextCacheHitCount || 0
      ),
      interactiveContextCacheMissCount: Number(
        pythonMetrics.interactiveContextCacheMissCount || 0
      ),
      interactiveContextBuildCount: Number(
        pythonMetrics.interactiveContextBuildCount || 0
      ),
      interactiveFeatureTraversalCount: Number(
        pythonMetrics.interactiveFeatureTraversalCount || 0
      ),
      selectorSafetyScopeBuildCount: Number(
        pythonMetrics.selectorSafetyScopeBuildCount || 0
      ),
      preparedInputCacheEvictionCount: Number(
        pythonMetrics.preparedInputCacheEvictionCount || 0
      ),
      preparedInputCacheRetainedBytes: Number(
        pythonMetrics.preparedInputCacheRetainedBytes || 0
      ),
      preparedInputCacheMutationViolationCount: Number(
        pythonMetrics.preparedInputCacheMutationViolationCount || 0
      )
    },
    featureCatalogStructural: {
      featureCatalogSvgParseCount: Number(
        pythonMetrics.featureCatalogSvgParseCount || 0
      ),
      featureCatalogFullDomTraversalCount: Number(
        pythonMetrics.featureCatalogFullDomTraversalCount || 0
      ),
      featureCatalogRenderedFeatureCollectionCount: Number(
        pythonMetrics.featureCatalogRenderedFeatureCollectionCount || 0
      ),
      featureCatalogCatalogOnlyDomMutationCount: Number(
        pythonMetrics.featureCatalogCatalogOnlyDomMutationCount || 0
      ),
      featureCatalogDomElementCount: Number(
        pythonMetrics.featureCatalogDomElementCount || 0
      ),
      featureCatalogFeatureCandidateCount: Number(
        pythonMetrics.featureCatalogFeatureCandidateCount || 0
      ),
      featureCatalogUniqueRenderedFeatureCount: Number(
        pythonMetrics.featureCatalogUniqueRenderedFeatureCount || 0
      ),
      featureCatalogMatchCandidateCount: Number(
        pythonMetrics.featureCatalogMatchCandidateCount || 0
      ),
      featureCatalogAnnotationCandidateCount: Number(
        pythonMetrics.featureCatalogAnnotationCandidateCount || 0
      )
    },
    workerResultPreparationMs: durationOrZero(
      events,
      'result-object-conversion-start',
      'result-transport-ready'
    ),
    workerCleanupMs: durationOrZero(
      events,
      'worker-cleanup-start',
      'worker-cleanup-end'
    ),
    workerArtifactIdentityMs: durationOrZero(
      events,
      'worker-artifact-identity-start',
      'worker-artifact-identity-end'
    ),
    resultTransportDecodeMs,
    candidateResultValidationMs: durationOrZero(
      events,
      'candidate-result-validation-start',
      'candidate-result-validation-end'
    ),
    featureCatalogAdoptionMs: durationOrZero(
      events,
      'feature-catalog-adoption-start',
      'feature-catalog-adoption-end'
    ),
    legacyFeatureOverrideMigrationMs: durationOrZero(
      events,
      'legacy-feature-override-migration-start',
      'legacy-feature-override-migration-end'
    ),
    ruleDerivedFillOverrideMs: durationOrZero(
      events,
      'rule-derived-fill-override-start',
      'rule-derived-fill-override-end'
    ),
    svgAdmissionMs: durationOrZero(
      events,
      'svg-admission-start',
      'svg-admission-end'
    ),
    resultAdmissionMs: durationOrZero(
      events,
      'result-admission-start',
      'result-admission-end'
    ),
    previewResultCommitMs: durationOrZero(
      events,
      'preview-result-commit-start',
      'preview-result-commit-end'
    ),
    historyAfterMs: durationOrZero(
      events,
      'generate-history-after-start',
      'generate-history-after-end'
    ),
    historyBeforeStructural,
    historyAfterStructural,
    historyStructural: Object.fromEntries(HISTORY_DIAGNOSTIC_FIELDS.map((name) => [
      name,
      Number(historyBeforeStructural[name] || 0) + Number(historyAfterStructural[name] || 0)
    ]))
  };
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
      const history = window.__GBDRAW_HISTORY__;
      const diagnosticFields = [
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
      ];
      const historyBefore = history.getDiagnostics();
      const undoCountBefore = history.getUndoCount();
      const redoCountBefore = history.getRedoCount();
      const started = performance.now();
      try {
        const result = await app.runAnalysis();
        const selected = app.results?.[app.selectedResultIndex];
        const historyAfter = history.getDiagnostics();
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
          previewVisible: Boolean(document.querySelector('.shadow-xl.origin-top > svg')),
          undoCountBefore,
          undoCountAfter: history.getUndoCount(),
          redoCountBefore,
          redoCountAfter: history.getRedoCount(),
          historyStructural: Object.fromEntries(diagnosticFields.map((name) => [
            name,
            Number(historyAfter[name] || 0) - Number(historyBefore[name] || 0)
          ]))
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
  const loadedPreviewIdentity = svgSourceIdentity(loaded.originalPreview);
  const preFirstGenerateActiveIntent = await activeIntentSummary(page);
  const loadProbe = await probeSnapshot(page);
  const ready = loadProbe.lifecycle.find(({ name }) => name === 'interactiveReady');
  const previewMetrics = Object.fromEntries(
    PREVIEW_ZERO_METRICS.map((name) => [name, Number(loadProbe.metrics[name] || 0)])
  );
  const preflightStructural = {
    proteinManifestFullValidationCount: Number(
      loadProbe.metrics.currentSessionPreflightProteinManifestValidationCount || 0
    ),
    proteinRawTextValidationCount: Number(
      loadProbe.metrics.currentSessionPreflightProteinRawTextValidationCount || 0
    )
  };
  expect(loaded.previewVisible).toBe(true);
  expect(loaded.resultCount).toBe(1);
  expect(ready).toMatchObject({ status: 'success', degradedRecovery: false });
  expect(previewMetrics).toEqual(Object.fromEntries(
    PREVIEW_ZERO_METRICS.map((name) => [name, 0])
  ));
  expect(preflightStructural.proteinManifestFullValidationCount).toBe(1);
  expect(preflightStructural.proteinRawTextValidationCount).toBeGreaterThan(0);
  expect(loadProbe.metrics.currentWriterActiveConfigRestoreCount).toBe(1);
  expect(loadProbe.metrics.activeConfigCanonicalOverwriteCount || 0).toBe(0);
  expect(await featureCatalogSummary(page)).toEqual({ schema: 3, itemCount: 1 });
  expect(preFirstGenerateActiveIntent.linearComparisonPlan).toEqual({
    mode: 'adjacent',
    defaultSource: 'losat',
    edges: []
  });
  expect(preFirstGenerateActiveIntent.linearRecordLayout).toMatchObject({
    enabled: true,
    recordGap: 48
  });
  expect(preFirstGenerateActiveIntent.linearRecordLayout.rows).toHaveLength(11);
  await installCanonicalRequestCapture(page);

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
  const firstGeneratedSvg = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    return String(app.results?.[app.selectedResultIndex]?.content || '');
  });
  const firstGeneratedIdentity = svgSourceIdentity(firstGeneratedSvg);

  await page.evaluate(() => window.__GBDRAW_VIBRIO_GENERATE_PROBE__.reset());
  const second = await invokeGeneration(page, terminal, runnerEvidence, terminalSignal);
  const secondProbe = !page.isClosed() && !terminal.pageCrashed
    ? await probeSnapshot(page)
    : { metrics: {}, details: [], lifecycle: [] };
  const secondRecovery = await assertRecoverableErrorState(page, second, loaded.originalPreview);
  expect(second.classification, JSON.stringify(second, null, 2)).toBe('success');
  expect(second.outcome).toMatchObject({
    status: 'ok',
    resultCommitted: true,
    previewVisible: true
  });
  const secondGeneratedSvg = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    return String(app.results?.[app.selectedResultIndex]?.content || '');
  });
  const canonicalRequestCapture = await canonicalRequestEvidence(page);
  const capturedCanonicalRequests = canonicalRequestCapture.requests;
  expect(capturedCanonicalRequests).toHaveLength(2);
  expect(capturedCanonicalRequests[0]).toMatchObject({
    schema: 5,
    mode: 'linear',
    recordCount: 11,
    comparisonCount: 2,
    comparisonKinds: ['collinearityResult', 'generatedProteinComparison'],
    layout: {
      recordGapPx: 48,
      multiRecordPositions: [
        '#1@1', '#2@1', '#3@1', '#4@2', '#5@2', '#6@3',
        '#7@3', '#8@4', '#9@4', '#10@5', '#11@5'
      ]
    }
  });
  expect(canonicalRequestCapture.repeatComparison).toMatchObject({ equivalent: true });
  const loadedComparisonSummary = svgComparisonSummary(loaded.originalPreview);
  const firstComparisonSummary = svgComparisonSummary(firstGeneratedSvg);
  const secondComparisonSummary = svgComparisonSummary(secondGeneratedSvg);
  expect(loadedComparisonSummary).toEqual({
    comparisonGroups: 16,
    pairwiseMatches: 633,
    comparisonLegends: 2
  });
  expect(firstComparisonSummary).toEqual(loadedComparisonSummary);
  expect(secondComparisonSummary).toEqual(firstComparisonSummary);
  for (const [stage, svg] of [
    ['loaded preview', loaded.originalPreview],
    ['first generation', firstGeneratedSvg],
    ['second generation', secondGeneratedSvg]
  ]) {
    expect(
      svgLinearGeometrySummary(svg),
      `Unexpected semantic Linear geometry in ${stage}.`
    ).toEqual(EXPECTED_VIBRIO_LINEAR_GEOMETRY);
  }
  const generatedFeatureCatalogDigest = await featureCatalogDigest(page);

  const activity = second.worker;
  const runs = activity?.instances?.[0]?.runs || [];
  const lifecycleNames = secondProbe.lifecycle.map(({ name }) => name);
  const resultSvgCharacters = [firstProbe, secondProbe].map((probe) => Number(
    probe.lifecycle.find(({ name }) => name === 'result-svg-characters')?.value || 0
  ));
  const artifactFingerprints = [firstProbe, secondProbe].map((probe) => String(
    probe.lifecycle.find(({ name }) => name === 'worker-artifact-identity-end')
      ?.fingerprint || ''
  ));
  const firstStructural = {
    canonicalReplayFullSerializationCount: Number(
      firstProbe.metrics.canonicalReplayFullSerializationCount || 0
    ),
    resourceMaterializationCount: Number(
      firstProbe.metrics.resourceMaterializationCount || 0
    ),
    resourceReencodeCount: Number(firstProbe.metrics.resourceReencodeCount || 0),
    resultBinaryDecodeCount: Number(firstProbe.metrics.resultBinaryDecodeCount || 0),
    resultBinaryDecodedBytes: Number(firstProbe.metrics.resultBinaryDecodedBytes || 0),
    resultMetadataBinaryDecodeCount: Number(
      firstProbe.metrics.resultMetadataBinaryDecodeCount || 0
    ),
    resultMetadataBinaryDecodedBytes: Number(
      firstProbe.metrics.resultMetadataBinaryDecodedBytes || 0
    ),
    legacyFeatureOverrideFullDescriptorComparisonCount: Number(
      firstProbe.metrics.legacyFeatureOverrideFullDescriptorComparisonCount || 0
    ),
    legacyFeatureOverrideIndexedDescriptorComparisonCount: Number(
      firstProbe.metrics.legacyFeatureOverrideIndexedDescriptorComparisonCount || 0
    ),
    legacyFeatureOverrideLegacyFeaturesVisited: Number(
      firstProbe.metrics.legacyFeatureOverrideLegacyFeaturesVisited || 0
    ),
    legacyFeatureOverrideMigrationCallCount: Number(
      firstProbe.metrics.legacyFeatureOverrideMigrationCallCount || 0
    ),
    legacyFeatureOverrideScanSkipCount: Number(
      firstProbe.metrics.legacyFeatureOverrideScanSkipCount || 0
    )
  };
  const secondStructural = {
    canonicalReplayFullSerializationCount: Number(
      secondProbe.metrics.canonicalReplayFullSerializationCount || 0
    ),
    resourceMaterializationCount: Number(
      secondProbe.metrics.resourceMaterializationCount || 0
    ),
    resourceReencodeCount: Number(secondProbe.metrics.resourceReencodeCount || 0),
    resultBinaryDecodeCount: Number(secondProbe.metrics.resultBinaryDecodeCount || 0),
    resultBinaryDecodedBytes: Number(secondProbe.metrics.resultBinaryDecodedBytes || 0),
    resultMetadataBinaryDecodeCount: Number(
      secondProbe.metrics.resultMetadataBinaryDecodeCount || 0
    ),
    resultMetadataBinaryDecodedBytes: Number(
      secondProbe.metrics.resultMetadataBinaryDecodedBytes || 0
    ),
    legacyFeatureOverrideFullDescriptorComparisonCount: Number(
      secondProbe.metrics.legacyFeatureOverrideFullDescriptorComparisonCount || 0
    ),
    legacyFeatureOverrideIndexedDescriptorComparisonCount: Number(
      secondProbe.metrics.legacyFeatureOverrideIndexedDescriptorComparisonCount || 0
    ),
    legacyFeatureOverrideLegacyFeaturesVisited: Number(
      secondProbe.metrics.legacyFeatureOverrideLegacyFeaturesVisited || 0
    ),
    legacyFeatureOverrideMigrationCallCount: Number(
      secondProbe.metrics.legacyFeatureOverrideMigrationCallCount || 0
    ),
    legacyFeatureOverrideScanSkipCount: Number(
      secondProbe.metrics.legacyFeatureOverrideScanSkipCount || 0
    )
  };
  const resultTransportEvents = [firstProbe, secondProbe].flatMap((probe) => (
    probe.lifecycle.filter(({ name }) => name === 'result-transport-ready')
  ));
  const structural = {
    canonicalReplayFullSerializationCount:
      firstStructural.canonicalReplayFullSerializationCount
      + secondStructural.canonicalReplayFullSerializationCount,
    workerBase64ResourceCloneCount: runs.filter(
      ({ hasBase64ResourceTable }) => hasBase64ResourceTable
    ).length,
    workerBase64ResourceJsonStringifyCount: lifecycleNames.filter(
      (name) => name.startsWith('worker-resources-json-')
    ).length,
    resourceMaterializationCount:
      firstStructural.resourceMaterializationCount
      + secondStructural.resourceMaterializationCount,
    resourceReencodeCount:
      firstStructural.resourceReencodeCount + secondStructural.resourceReencodeCount,
    resultBinaryDecodeCount:
      firstStructural.resultBinaryDecodeCount + secondStructural.resultBinaryDecodeCount,
    resultBinaryDecodedBytes:
      firstStructural.resultBinaryDecodedBytes + secondStructural.resultBinaryDecodedBytes,
    resultMetadataBinaryDecodeCount:
      firstStructural.resultMetadataBinaryDecodeCount
      + secondStructural.resultMetadataBinaryDecodeCount,
    resultMetadataBinaryDecodedBytes:
      firstStructural.resultMetadataBinaryDecodedBytes
      + secondStructural.resultMetadataBinaryDecodedBytes
  };
  const firstPhaseAttribution = generationPhaseAttribution(first.outcome, firstProbe);
  const secondPhaseAttribution = generationPhaseAttribution(second.outcome, secondProbe);
  const expectedWarmGenerateHistory = {
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
    artifactReplacementHistoryEntryCount: 0
  };

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
  expect(secondPhaseAttribution).toMatchObject({
    workerInitializationMs: 0,
    workerInitializationReused: true,
    newlyStagedResourceBytes: 0,
    losatCacheHits: 0,
    losatCacheMisses: 0
  });
  expect(firstPhaseAttribution).toMatchObject({
    losatCacheHits: 47,
    losatCacheMisses: 0
  });
  expect(second.outcome.historyStructural).toEqual(expectedWarmGenerateHistory);
  expect(secondPhaseAttribution.historyStructural).toMatchObject(
    expectedWarmGenerateHistory
  );
  expect(second.outcome.undoCountAfter).toBe(second.outcome.undoCountBefore);
  expect(second.outcome.redoCountAfter).toBe(second.outcome.redoCountBefore);
  expect(first.outcome.historyStructural).toMatchObject({
    ...expectedWarmGenerateHistory,
    artifactReplacementHistoryEntryCount: 1
  });
  expect(secondStructural.resourceMaterializationCount).toBe(0);
  expect(secondStructural).toMatchObject({
    legacyFeatureOverrideFullDescriptorComparisonCount: 0,
    legacyFeatureOverrideIndexedDescriptorComparisonCount: 0,
    legacyFeatureOverrideLegacyFeaturesVisited: 0,
    legacyFeatureOverrideMigrationCallCount: 2,
    legacyFeatureOverrideScanSkipCount: 2
  });
  expect(secondPhaseAttribution.pythonInvocationMs).toBeGreaterThan(0);
  expect(firstPhaseAttribution.preparedInputCacheStructural).toMatchObject({
    decodedResourceCacheHitCount: 0,
    decodedResourceCacheMissCount: 1,
    decodedResourceBuildCount: 1,
    parsedSourceCacheHitCount: 0,
    resolvedRecordCacheHitCount: 0,
    resolvedRecordCacheMissCount: 1,
    resolvedRecordBuildCount: 1,
    interactiveContextCacheHitCount: 0,
    interactiveContextCacheMissCount: 1,
    interactiveContextBuildCount: 1,
    interactiveFeatureTraversalCount: 1,
    selectorSafetyScopeBuildCount: 0,
    preparedInputCacheMutationViolationCount: 0
  });
  expect(firstPhaseAttribution.preparedInputCacheStructural
    .parsedSourceCacheMissCount).toBeGreaterThan(0);
  expect(firstPhaseAttribution.preparedInputCacheStructural
    .parsedSourceParseCount).toBe(
      firstPhaseAttribution.preparedInputCacheStructural.parsedSourceCacheMissCount
    );
  expect(firstPhaseAttribution.preparedInputCacheStructural
    .preparedInputCacheRetainedBytes).toBeGreaterThan(0);
  expect(secondPhaseAttribution.preparedInputCacheStructural).toMatchObject({
    decodedResourceCacheHitCount: 1,
    decodedResourceCacheMissCount: 0,
    decodedResourceBuildCount: 0,
    parsedSourceCacheMissCount: 0,
    parsedSourceParseCount: 0,
    resolvedRecordCacheHitCount: 1,
    resolvedRecordCacheMissCount: 0,
    resolvedRecordBuildCount: 0,
    interactiveContextCacheHitCount: 1,
    interactiveContextCacheMissCount: 0,
    interactiveContextBuildCount: 0,
    interactiveFeatureTraversalCount: 0,
    selectorSafetyScopeBuildCount: 0,
    preparedInputCacheMutationViolationCount: 0
  });
  expect(secondPhaseAttribution.preparedInputCacheStructural
    .parsedSourceCacheHitCount).toBeGreaterThan(0);
  expect(secondPhaseAttribution.preparedInputCacheStructural
    .preparedInputCacheRetainedBytes).toBeGreaterThan(0);
  expect(firstPhaseAttribution.featureCatalogStructural).toMatchObject({
    featureCatalogSvgParseCount: 1,
    featureCatalogFullDomTraversalCount: 1,
    featureCatalogRenderedFeatureCollectionCount: 1,
    featureCatalogCatalogOnlyDomMutationCount: 0
  });
  expect(secondPhaseAttribution.featureCatalogStructural).toMatchObject({
    featureCatalogSvgParseCount: 1,
    featureCatalogFullDomTraversalCount: 1,
    featureCatalogRenderedFeatureCollectionCount: 1,
    featureCatalogCatalogOnlyDomMutationCount: 0
  });
  for (const phase of [firstPhaseAttribution, secondPhaseAttribution]) {
    const catalog = phase.featureCatalogStructural;
    expect(catalog.featureCatalogDomElementCount).toBeGreaterThan(0);
    expect(catalog.featureCatalogFeatureCandidateCount).toBeGreaterThan(0);
    expect(catalog.featureCatalogUniqueRenderedFeatureCount).toBeGreaterThan(0);
    expect(catalog.featureCatalogFeatureCandidateCount).toBeGreaterThanOrEqual(
      catalog.featureCatalogUniqueRenderedFeatureCount
    );
    expect(catalog.featureCatalogMatchCandidateCount).toBeGreaterThanOrEqual(0);
    expect(catalog.featureCatalogAnnotationCandidateCount).toBeGreaterThanOrEqual(0);
  }
  expect(secondProbe.lifecycle.filter(({ name }) => name === 'python-diagnostics')).toHaveLength(1);
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
  expect(artifactFingerprints).toHaveLength(2);
  expect(artifactFingerprints.every((fingerprint) => /^[0-9a-f]{64}$/.test(fingerprint)))
    .toBe(true);
  expect(artifactFingerprints[1]).toBe(artifactFingerprints[0]);
  expect(generatedFeatureCatalogDigest).toMatchObject({
    schema: 3,
    itemCount: 1
  });
  expect(generatedFeatureCatalogDigest.sha256).toMatch(/^[0-9a-f]{64}$/);
  expect(generatedFeatureCatalogDigest.jsonBytes).toBeGreaterThan(0);

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

  const loadedSvgPath = testInfo.outputPath('vibrio-loaded-preview.svg');
  const firstSvgPath = testInfo.outputPath('vibrio-generate-1.svg');
  const secondSvgPath = testInfo.outputPath('vibrio-generate-2.svg');
  writeFileSync(loadedSvgPath, loaded.originalPreview, 'utf8');
  writeFileSync(firstSvgPath, firstGeneratedSvg, 'utf8');
  writeFileSync(secondSvgPath, secondGeneratedSvg, 'utf8');
  const svgComparisonCommand = [
    'import sys',
    'from tests.utils.svg_compare import compare_svgs',
    'result = compare_svgs(sys.argv[1], sys.argv[2])',
    'print(result.message)',
    'print("\\n".join(result.differences))',
    'raise SystemExit(0 if result.equal else 1)'
  ].join(';');
  const loadedToFirstComparison = spawnSync(
    process.env.GBDRAW_PYTHON || 'python',
    ['-c', svgComparisonCommand, loadedSvgPath, firstSvgPath],
    { cwd: repoRoot, encoding: 'utf8' }
  );
  expect(
    loadedToFirstComparison.status,
    `${loadedToFirstComparison.stdout}\n${loadedToFirstComparison.stderr}`
  ).toBe(0);
  const repeatComparison = spawnSync(
    process.env.GBDRAW_PYTHON || 'python',
    ['-c', svgComparisonCommand, firstSvgPath, secondSvgPath],
    { cwd: repoRoot, encoding: 'utf8' }
  );
  expect(
    repeatComparison.status,
    `${repeatComparison.stdout}\n${repeatComparison.stderr}`
  ).toBe(0);

  const outputsExactlyEqual = firstGeneratedSvg === secondGeneratedSvg;
  expect(outputsExactlyEqual).toBe(true);
  let historyRoundTrip;
  if (outputsExactlyEqual) {
    historyRoundTrip = await page.evaluate((expectedCharacters) => {
      const history = window.__GBDRAW_HISTORY__;
      const app = window.__GBDRAW_APP__;
      const currentSvg = String(app.results?.[app.selectedResultIndex]?.content || '');
      return {
        undoCountBefore: history.getUndoCount(),
        historyEntryCreated: false,
        undoRedoAttempted: false,
        noOpGeneratePreserved: currentSvg.length === expectedCharacters,
        undoCountAfter: history.getUndoCount(),
        redoCountAfter: history.getRedoCount()
      };
    }, secondGeneratedSvg.length);
    expect(historyRoundTrip).toMatchObject({
      historyEntryCreated: false,
      undoRedoAttempted: false,
      noOpGeneratePreserved: true,
      undoCountAfter: second.outcome.undoCountAfter,
      redoCountAfter: second.outcome.redoCountAfter
    });
  } else {
    historyRoundTrip = await page.evaluate(async ({ firstSvg, secondSvg }) => {
      const history = window.__GBDRAW_HISTORY__;
      const app = window.__GBDRAW_APP__;
      const undoCountBefore = history.getUndoCount();
      const undoOk = await history.undo();
      const afterUndo = String(app.results?.[app.selectedResultIndex]?.content || '');
      const redoOk = await history.redo();
      const afterRedo = String(app.results?.[app.selectedResultIndex]?.content || '');
      return {
        undoCountBefore,
        historyEntryCreated: true,
        undoRedoAttempted: true,
        undoOk,
        redoOk,
        undoRestoredFirstGenerate: afterUndo === firstSvg,
        redoRestoredSecondGenerate: afterRedo === secondSvg,
        undoCountAfter: history.getUndoCount(),
        redoCountAfter: history.getRedoCount()
      };
    }, {
      firstSvg: firstGeneratedSvg,
      secondSvg: secondGeneratedSvg
    });
    expect(historyRoundTrip).toMatchObject({
      historyEntryCreated: true,
      undoRedoAttempted: true,
      undoOk: true,
      redoOk: true,
      undoRestoredFirstGenerate: true,
      redoRestoredSecondGenerate: true,
      redoCountAfter: 0
    });
  }

  const report = {
    load: {
      previewMetrics,
      preflightStructural,
      loadedPreviewIdentity,
      preFirstGenerateActiveIntent,
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
      firstPhaseAttribution,
      secondPhaseAttribution
    },
    perGenerationStructural: {
      first: firstStructural,
      second: secondStructural
    },
    deterministicOutput: {
      comparison: 'tests.utils.svg_compare.compare_svgs',
      loadedToFirstRelationship: 'publication-parity',
      loadedToFirstExactBytesEqual: loaded.originalPreview === firstGeneratedSvg,
      loadedPreviewIdentity,
      firstGeneratedIdentity,
      loadedComparisonSummary,
      firstComparisonSummary,
      secondComparisonSummary,
      capturedCanonicalRequests,
      loadedToFirstSemanticallyEquivalent: loadedToFirstComparison.status === 0,
      repeatedGenerateSemanticallyEquivalent: repeatComparison.status === 0,
      exactBytesEqual: firstGeneratedSvg === secondGeneratedSvg,
      firstCharacters: firstGeneratedSvg.length,
      secondCharacters: secondGeneratedSvg.length,
      featureCatalog: generatedFeatureCatalogDigest,
      artifactFingerprints,
      generatedMetadataFingerprintEqual:
        artifactFingerprints[0] === artifactFingerprints[1]
    },
    historyRoundTrip,
    structural,
    worker: afterPing
  };
  await testInfo.attach('vibrio-full-generation-metrics.json', {
    body: Buffer.from(JSON.stringify(report, null, 2)),
    contentType: 'application/json'
  });
  console.log(`Vibrio full-generation evidence: ${JSON.stringify({
    loadTimings: report.load.timings,
    loadedPreviewIdentity: report.load.loadedPreviewIdentity,
    preFirstGenerateActiveIntent: report.load.preFirstGenerateActiveIntent,
    firstGenerateMs: report.generations.firstElapsedMs,
    secondGenerateMs: report.generations.secondElapsedMs,
    referencedResourceCount: report.generations.referencedResourceCount,
    referencedDeclaredBytes: report.generations.referencedDeclaredBytes,
    firstTransferredBytes: report.generations.firstTransferredBytes,
    secondTransferredBytes: report.generations.secondTransferredBytes,
    secondPhaseAttribution: report.generations.secondPhaseAttribution,
    structural: report.structural,
    deterministicOutput: report.deterministicOutput,
    historyRoundTrip: report.historyRoundTrip
  })}`);
});
