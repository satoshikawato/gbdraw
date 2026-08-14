const { test, expect } = require('@playwright/test');
const { createHash } = require('node:crypto');
const {
  existsSync,
  mkdtempSync,
  readFileSync,
  rmSync
} = require('node:fs');
const { tmpdir } = require('node:os');
const { join, resolve } = require('node:path');
const { gunzipSync } = require('node:zlib');

const repoRoot = resolve(process.env.GBDRAW_REPO || process.cwd());
const fixtureDir = join(repoRoot, 'tests', 'fixtures', 'sessions');
const fixturePath = join(
  fixtureDir,
  'BGC0000708-BGC0000713.schema-v2.gbdraw-session.json.gz'
);
const expected = JSON.parse(readFileSync(join(
  fixtureDir,
  'BGC0000708-BGC0000713.schema-v2.expected.json'
), 'utf8'));
const CURRENT_SESSION_VERSION = 40;
const CURRENT_RENDER_REQUEST_SCHEMA = 5;
const CURRENT_PROTEIN_RAW_SCHEMA = 4;
const CURRENT_PROTEIN_DERIVED_SCHEMA = 3;

let sourceSession;
let savedSessionCount = 0;
let savedSessionDir;

const expandedSessionBytes = (path) => {
  const bytes = readFileSync(path);
  return bytes[0] === 0x1f && bytes[1] === 0x8b ? gunzipSync(bytes) : bytes;
};

const readSession = (path) => JSON.parse(expandedSessionBytes(path).toString('utf8'));

const importSession = async (page, path) => {
  const dialogPromise = page.waitForEvent('dialog', { timeout: 120000 });
  await page.locator('input[accept^=".json,"]').first().setInputFiles(path);
  const dialog = await dialogPromise;
  expect(dialog.message()).toBe('Session loaded successfully!');
  await dialog.accept();
  await page.waitForFunction(() => {
    const app = window.__GBDRAW_APP__;
    return app?.mode === 'linear' && app?.linearSeqs?.length === 5;
  }, null, { timeout: 120000 });
};

const saveSession = async (page) => {
  const downloadPromise = page.waitForEvent('download', { timeout: 120000 });
  const outcome = await page.evaluate(async () => ({
    result: await window.__GBDRAW_APP__.saveSessionWithTitle(),
    errorLog: window.__GBDRAW_APP__.errorLog
  }));
  if (outcome.result?.status !== 'saved') {
    downloadPromise.catch(() => {});
    throw new Error(`Session Save failed: ${JSON.stringify(outcome)}`);
  }
  const download = await downloadPromise;
  savedSessionCount += 1;
  const destination = join(
    savedSessionDir,
    `saved-${savedSessionCount}.gbdraw-session.json.gz`
  );
  await download.saveAs(destination);
  return destination;
};

const downloadSvg = async (page) => {
  const downloadPromise = page.waitForEvent('download', { timeout: 120000 });
  await page.evaluate(() => window.__GBDRAW_APP__.downloadSVG());
  const download = await downloadPromise;
  const path = await download.path();
  expect(path).toBeTruthy();
  return path;
};

const downloadLosatPair = async (page, pairIndex) => {
  const downloadPromise = page.waitForEvent('download', { timeout: 120000 });
  await page.evaluate(
    async (index) => window.__GBDRAW_APP__.downloadLosatPair(index),
    pairIndex
  );
  const download = await downloadPromise;
  const path = await download.path();
  expect(path).toBeTruthy();
  return path;
};

const downloadAllLosatPairs = async (page, expectedCount) => {
  const downloads = [];
  const capture = (download) => downloads.push(download);
  page.on('download', capture);
  try {
    await page.evaluate(async () => window.__GBDRAW_APP__.downloadLosatCache());
    await expect.poll(
      () => downloads.length,
      { timeout: 120000 }
    ).toBe(expectedCount);
    const paths = await Promise.all(downloads.map((download) => download.path()));
    expect(paths.every(Boolean)).toBe(true);
    return paths;
  } finally {
    page.off('download', capture);
  }
};

const dataRows = (text) => String(text || '')
  .split(/\r?\n/)
  .filter((line) => line.trim() && !line.trimStart().startsWith('#'))
  .map((line) => line.split('\t'));

const assertHydratedDownload = (text) => {
  const rows = dataRows(text);
  expect(rows.length).toBeGreaterThan(0);
  expect(rows.every((columns) => columns.length === 12)).toBe(true);
  for (const pattern of expected.downloads.forbiddenPatterns) {
    expect(text).not.toMatch(new RegExp(pattern, 'i'));
  }
  for (const columns of rows) {
    expect(columns[0]).toMatch(/^[A-Za-z0-9._%~-]+$/);
    expect(columns[1]).toMatch(/^[A-Za-z0-9._%~-]+$/);
  }
};

const assertHydratedDownloads = async (page) => {
  const pairText = readFileSync(
    await downloadLosatPair(page, expected.downloads.pairIndex),
    'utf8'
  );
  assertHydratedDownload(pairText);
  expect(dataRows(pairText)[0].slice(0, 2)).toEqual(expected.downloads.firstDataIds);

  const bulkPaths = await downloadAllLosatPairs(
    page,
    expected.downloads.bulkEntryCount
  );
  expect(bulkPaths).toHaveLength(expected.downloads.bulkEntryCount);
  for (const path of bulkPaths) {
    assertHydratedDownload(readFileSync(path, 'utf8'));
  }
};

const assertHydratedUtf8PromptBoundary = async (page) => {
  const probe = await page.evaluate(async () => {
    const {
      LOSAT_EXPORT_CONFIRM_THRESHOLD_BYTES,
      confirmHydratedLosatExport,
      totalHydratedLosatExportBytes
    } = await import('/gbdraw/web/js/app/run-analysis.js');
    const hydrated = { text: '界界' };
    const utf8Bytes = totalHydratedLosatExportBytes([hydrated]);
    let prompts = 0;
    const accepted = confirmHydratedLosatExport(
      utf8Bytes,
      () => {
        prompts += 1;
        return true;
      },
      5
    );
    return {
      threshold: LOSAT_EXPORT_CONFIRM_THRESHOLD_BYTES,
      codeUnits: hydrated.text.length,
      utf8Bytes,
      prompts,
      accepted
    };
  });
  expect(probe).toEqual({
    threshold: expected.downloads.hydratedPromptThresholdBytes,
    codeUnits: 2,
    utf8Bytes: 6,
    prompts: 1,
    accepted: true
  });
};

const installUserUploadedTsv = async (page) => page.evaluate(
  async ({ filename, text }) => {
    const bytes = new TextEncoder().encode(text);
    const app = window.__GBDRAW_APP__;
    app.losatProgram = 'blastn';
    const file = new File(
      [bytes],
      filename,
      { type: 'text/tab-separated-values', lastModified: 0 }
    );
    app.linearComparisonPlan.mode = 'adjacent';
    app.linearComparisonPlan.defaultSource = 'upload';
    app.linearComparisonPlan.edges.splice(0, app.linearComparisonPlan.edges.length, {
      id: 'acceptance-upload-edge',
      queryUid: app.linearSeqs[0].uid,
      subjectUid: app.linearSeqs[1].uid,
      included: true,
      fileActive: true,
      losatFilenameActive: false,
      source: 'upload',
      file,
      losatFilename: ''
    });
    return Array.from(
      new Uint8Array(await file.arrayBuffer())
    );
  },
  {
    filename: expected.downloads.userUploadedFilename,
    text: expected.downloads.userUploadedTsv
  }
);

const readFirstUploadedTsv = async (page) => page.evaluate(async () => {
  const file = window.__GBDRAW_APP__.linearComparisonPlan.edges[0]?.file;
  const { readFileBytes } = await import('/gbdraw/web/js/services/file-content-cache.js');
  return {
    name: file?.name || '',
    bytes: file ? Array.from(await readFileBytes(file)) : []
  };
});

const assertCurrentSessionBoundary = (session) => {
  const rawEntries = session.losatCache?.entries || [];
  const derivedEntries = session.losatDerivedCache?.entries || [];
  expect(session.format).toBe('gbdraw-session');
  expect(session.version).toBe(CURRENT_SESSION_VERSION);
  expect(session.renderRequest?.schema).toBe(CURRENT_RENDER_REQUEST_SCHEMA);
  expect(rawEntries.every((entry) => (
    entry?.identityKind !== 'protein' ||
    entry?.schema === CURRENT_PROTEIN_RAW_SCHEMA
  ))).toBe(true);
  expect(derivedEntries.every(
    (entry) => (
      entry?.schema === CURRENT_PROTEIN_DERIVED_SCHEMA &&
      entry?.idEncoding === 'runtime-handle-v1'
    )
  )).toBe(true);
};

const assertRenamedZeroMtimeResources = (session) => {
  const genbankResources = Object.values(session.resources || {})
    .filter((resource) => resource?.kind === 'genbank');
  expect(genbankResources).toHaveLength(5);
  expect(genbankResources.every((resource) => resource.lastModified === 0)).toBe(true);
  expect(genbankResources.map((resource) => resource.name)).toEqual(
    expect.arrayContaining(Array.from({ length: 5 }, (_, index) => (
      expect.stringContaining(`renamed-cache-input-${index + 1}.gbk`)
    )))
  );
  const originalNames = Object.values(session.webFiles?.resourceOriginalNames || {});
  expect(originalNames).toEqual(expect.arrayContaining(
    Array.from({ length: 5 }, (_, index) => `renamed-cache-input-${index + 1}.gbk`)
  ));
};

const assertCurrentProteinArtifacts = (session) => {
  assertCurrentSessionBoundary(session);
  const manifest = session.proteinIdentityManifest;
  const entries = session.losatCache?.entries || [];
  const derivedEntries = session.losatDerivedCache?.entries || [];
  expect(manifest?.schema).toBe(2);
  expect(entries).toHaveLength(expected.totalPairs);
  expect(entries.every((entry) => (
    entry?.schema === 4 &&
    entry?.kind === 'raw-losat' &&
    entry?.identityKind === 'protein' &&
    entry?.idEncoding === 'runtime-handle-v1' &&
    entry?.program === 'blastp'
  ))).toBe(true);
  expect(derivedEntries).toEqual([]);
  const legacyCandidates = session.legacyArtifacts?.proteinRawCandidates?.entries || [];
  expect(legacyCandidates)
    .toHaveLength(expected.storedRawEntries - expected.totalPairs);
  expect(legacyCandidates.every((candidate) => (
    candidate.state === 'pending' &&
    (sourceSession.losatCache?.entries || []).some(
      (entry) => JSON.stringify(entry) === JSON.stringify(candidate.originalEntry)
    )
  ))).toBe(true);
  expect(session.legacyArtifacts?.proteinDerivedEvidence?.entries || [])
    .toEqual(sourceSession.losatDerivedCache?.entries || []);

  const runtimeOwners = new Map();
  for (const [instanceKey, instance] of Object.entries(manifest.recordInstances || {})) {
    expect(manifest.recordAnalyses?.[instance.recordAnalysisId]).toBeTruthy();
    expect(instance.runtimeBindingHash).toBeTruthy();
    expect(instance.displayBindingHash).toBeTruthy();
    expect(Object.keys(instance.runtimeIds || {}).sort())
      .toEqual(Object.keys(instance.featureMetadata || {}).sort());
    for (const [featureId, runtimeHandle] of Object.entries(instance.runtimeIds || {})) {
      expect(featureId).toMatch(/^f_/);
      expect(runtimeHandle).toMatch(/^h_[a-z2-7]{26}$/);
      expect(runtimeOwners.has(runtimeHandle)).toBe(false);
      runtimeOwners.set(runtimeHandle, { instanceKey, featureId });
    }
  }

  let resolvedReferenceCount = 0;
  for (const entry of entries) {
    const query = manifest.recordInstances?.[entry.queryRecordInstanceKey];
    const subject = manifest.recordInstances?.[entry.subjectRecordInstanceKey];
    expect(query?.runtimeBindingHash).toBe(entry.queryRuntimeBindingHash);
    expect(subject?.runtimeBindingHash).toBe(entry.subjectRuntimeBindingHash);
    expect(manifest.recordAnalyses?.[query.recordAnalysisId]?.proteinSetHash)
      .toBe(entry.queryProteinSetHash);
    expect(manifest.recordAnalyses?.[subject.recordAnalysisId]?.proteinSetHash)
      .toBe(entry.subjectProteinSetHash);
    const queryIds = new Set(Object.values(query.runtimeIds || {}));
    const subjectIds = new Set(Object.values(subject.runtimeIds || {}));
    for (const rawLine of String(entry.text || '').split(/\r?\n/)) {
      const line = rawLine.trim();
      if (!line || line.startsWith('#')) continue;
      const [queryId, subjectId] = rawLine.split('\t');
      expect(queryId).toMatch(/^h_[a-z2-7]{26}$/);
      expect(subjectId).toMatch(/^h_[a-z2-7]{26}$/);
      expect(queryIds.has(queryId)).toBe(true);
      expect(subjectIds.has(subjectId)).toBe(true);
      expect(runtimeOwners.get(queryId)?.instanceKey).toBe(entry.queryRecordInstanceKey);
      expect(runtimeOwners.get(subjectId)?.instanceKey).toBe(entry.subjectRecordInstanceKey);
      resolvedReferenceCount += 2;
    }
  }
  expect(resolvedReferenceCount).toBeGreaterThan(0);

};

const generateWithTelemetry = async (page) => page.evaluate(async () => {
  const app = window.__GBDRAW_APP__;
  const result = await app.runAnalysis();
  return {
    result,
    errorSummary: String(app.errorLog?.summary || ''),
    errorDetails: Array.isArray(app.errorLog?.details)
      ? app.errorLog.details.map((detail) => String(detail))
      : [],
    telemetry: app.lastRunInfo?.losatTelemetry || window.__GBDRAW_LAST_LOSAT_TELEMETRY__ || null,
    executorCalls: Number(window.__GBDRAW_LOSAT_EXECUTOR_CALLS__ || 0)
  };
});

const settleAppRender = async (page) => page.evaluate(async () => {
  await window.Vue.nextTick();
  await new Promise((resolveFrame) => requestAnimationFrame(resolveFrame));
  await window.Vue.nextTick();
});

const migrationUiSnapshot = async (page) => page.evaluate(async () => {
  const { state } = await import('/gbdraw/web/js/state.js');
  return JSON.parse(JSON.stringify({
    orthogroups: state.orthogroups.value,
    selectedOrthogroupAlignmentFeature: state.selectedOrthogroupAlignmentFeature.value,
    selectedOrthogroupId: state.selectedOrthogroupId.value,
    extractedFeatures: state.extractedFeatures.value,
    biologicalFeatures: state.biologicalFeatures.value
  }));
});

const cancelDuringRender = async (page) => page.evaluate(async () => {
  const app = window.__GBDRAW_APP__;
  const { state } = await import('/gbdraw/web/js/state.js');
  const {
    CANONICAL_REQUEST_SCHEMA
  } = await import('/gbdraw/web/js/services/session-request.js');
  const before = {
    proteinIdentityManifest: state.proteinIdentityManifest.value,
    legacyProteinRawCandidates: state.legacyProteinRawCandidates.value,
    legacyProteinDerivedEvidence: state.legacyProteinDerivedEvidence.value,
    losatCache: Array.from(state.losatCache.value.entries()),
    losatDerivedCache: Array.from(state.losatDerivedCache.value.entries()),
    losatCacheInfo: state.losatCacheInfo.value
  };
  const originalWorkerPostMessage = Worker.prototype.postMessage;
  let releaseRunRequest;
  const runRequestIssued = new Promise((resolve) => {
    releaseRunRequest = resolve;
  });
  let heldCanonicalRun = null;
  let cancelInvoked = false;
  Worker.prototype.postMessage = function holdFirstCanonicalRun(message, ...args) {
    if (
      !heldCanonicalRun &&
      message?.type === 'run' &&
      message?.payload?.request?.schema === CANONICAL_REQUEST_SCHEMA
    ) {
      heldCanonicalRun = { worker: this, message, args };
      releaseRunRequest();
      return;
    }
    return originalWorkerPostMessage.call(this, message, ...args);
  };
  try {
    const runPromise = app.runAnalysis();
    await runRequestIssued;
    app.cancelGeneration();
    cancelInvoked = true;
    const result = await runPromise;
    const sameMapEntries = (entries, current) => (
      entries.length === current.size &&
      entries.every(([key, value]) => current.get(key) === value)
    );
    return {
      result,
      runRequestIssued: Boolean(heldCanonicalRun),
      cancelInvoked,
      errorSummary: String(app.errorLog?.summary || ''),
      executorCalls: Number(window.__GBDRAW_LOSAT_EXECUTOR_CALLS__ || 0),
      authorityRestored: (
        state.proteinIdentityManifest.value === before.proteinIdentityManifest &&
        state.legacyProteinRawCandidates.value === before.legacyProteinRawCandidates &&
        state.legacyProteinDerivedEvidence.value ===
          before.legacyProteinDerivedEvidence &&
        sameMapEntries(before.losatCache, state.losatCache.value) &&
        sameMapEntries(before.losatDerivedCache, state.losatDerivedCache.value) &&
        state.losatCacheInfo.value === before.losatCacheInfo
      )
    };
  } finally {
    Worker.prototype.postMessage = originalWorkerPostMessage;
  }
});

const failRendererAfterMigration = async (page) => page.evaluate(async () => {
  const app = window.__GBDRAW_APP__;
  const { state } = await import('/gbdraw/web/js/state.js');
  const {
    CANONICAL_REQUEST_SCHEMA
  } = await import('/gbdraw/web/js/services/session-request.js');
  const before = {
    proteinIdentityManifest: state.proteinIdentityManifest.value,
    legacyProteinRawCandidates: state.legacyProteinRawCandidates.value,
    legacyProteinDerivedEvidence: state.legacyProteinDerivedEvidence.value,
    losatCache: Array.from(state.losatCache.value.entries()),
    losatDerivedCache: Array.from(state.losatDerivedCache.value.entries()),
    losatCacheInfo: state.losatCacheInfo.value
  };
  const originalWorkerPostMessage = Worker.prototype.postMessage;
  let rendererFailureInjected = false;
  Worker.prototype.postMessage = function (...args) {
    const message = args[0];
    if (
      !rendererFailureInjected &&
      message?.type === 'run' &&
      message?.payload?.request?.schema === CANONICAL_REQUEST_SCHEMA &&
      Array.isArray(message?.payload?.resourceManifest) &&
      Array.isArray(message?.payload?.stagedResources) &&
      !Object.hasOwn(message.payload, 'resources')
    ) {
      rendererFailureInjected = true;
      message.payload.request = null;
    }
    return originalWorkerPostMessage.apply(this, args);
  };
  try {
    const result = await app.runAnalysis();
    const sameMapEntries = (entries, current) => (
      entries.length === current.size &&
      entries.every(([key, value]) => current.get(key) === value)
    );
    return {
      result,
      errorSummary: String(app.errorLog?.summary || ''),
      executorCalls: Number(window.__GBDRAW_LOSAT_EXECUTOR_CALLS__ || 0),
      rendererFailureInjected,
      authorityRestored: (
        state.proteinIdentityManifest.value === before.proteinIdentityManifest &&
        state.legacyProteinRawCandidates.value === before.legacyProteinRawCandidates &&
        state.legacyProteinDerivedEvidence.value ===
          before.legacyProteinDerivedEvidence &&
        sameMapEntries(before.losatCache, state.losatCache.value) &&
        sameMapEntries(before.losatDerivedCache, state.losatDerivedCache.value) &&
        state.losatCacheInfo.value === before.losatCacheInfo
      )
    };
  } finally {
    Worker.prototype.postMessage = originalWorkerPostMessage;
  }
});

const assertExpectedTelemetry = (run, migrationExpected = expected) => {
  const diagnostics = `Generation diagnostics:\n${JSON.stringify(run, null, 2)}`;
  expect(run.result, diagnostics).toEqual({ status: 'ok' });
  expect(run.errorSummary, diagnostics).toBe('');
  expect(run.errorDetails, diagnostics).toEqual([]);
  expect(run.telemetry, diagnostics).toMatchObject({
    totalPairs: migrationExpected.totalPairs,
    cacheHits: migrationExpected.cacheHits,
    cacheMisses: migrationExpected.cacheMisses,
    uniqueJobs: migrationExpected.uniqueJobs,
    workerCalls: migrationExpected.workerCalls
  });
  expect(run.executorCalls, diagnostics).toBe(migrationExpected.workerCalls);
};

const inspectLayout = async (page) => {
  return page.evaluate(async () => {
    const { state } = await import('/gbdraw/web/js/state.js');
    const geometry = JSON.parse(JSON.stringify(state.trackSlotResolvedGeometry.value || null));
    const svg = state.svgContainer.value?.querySelector('svg');
    if (!geometry || !svg) {
      return { error: 'Resolved Linear geometry or the preview SVG is unavailable.' };
    }

    const epsilon = 1e-6;
    const records = (Array.isArray(geometry.records) ? geometry.records : [])
      .slice()
      .sort((left, right) => Number(left.recordIndex) - Number(right.recordIndex));
    const directGroups = Array.from(svg.children).filter(
      (element) => element.tagName?.toLowerCase() === 'g'
    );
    const comparisonGroups = directGroups
      .filter((element) => (
        element.hasAttribute('data-query-row') &&
        element.hasAttribute('data-subject-row')
      ))
      .sort((left, right) => (
        Number(left.getAttribute('data-query-row')) -
        Number(right.getAttribute('data-query-row'))
      ));
    const activeBoundaries = new Set(comparisonGroups.map((element) => (
      `${element.getAttribute('data-query-row')}:${element.getAttribute('data-subject-row')}`
    )));

    const xOverlaps = (left, right) => (
      Math.min(Number(left.xEndPx), Number(right.xEndPx)) -
      Math.max(Number(left.xStartPx), Number(right.xStartPx))
    ) > epsilon;
    const pairIsEligible = (left, right, boundaryActive) => (
      (left.kind === 'body' && right.kind === 'body') ||
      (boundaryActive && left.kind === 'comparison' && right.kind === 'comparison') ||
      (left.kind === 'definition' && right.kind === 'definition') ||
      (left.kind === 'definition' && right.kind === 'body') ||
      (left.kind === 'body' && right.kind === 'definition')
    );

    const collisionChecks = [];
    const axisGaps = [];
    for (let boundary = 0; boundary < records.length - 1; boundary += 1) {
      const nextRow = boundary + 1;
      const current = records[boundary];
      const following = records[nextRow];
      const boundaryActive = activeBoundaries.has(`${boundary}:${nextRow}`);
      axisGaps.push(Number(following?.axisYpx) - Number(current?.axisYpx));
      for (const currentBand of current?.collisionBands || []) {
        for (const nextBand of following?.collisionBands || []) {
          if (
            !xOverlaps(currentBand, nextBand) ||
            !pairIsEligible(currentBand, nextBand, boundaryActive)
          ) continue;
          collisionChecks.push({
            boundaryRow: boundary,
            kindPair: `${currentBand.kind}:${nextBand.kind}`,
            gapPx: Number(nextBand.absoluteTopPx) - Number(currentBand.absoluteBottomPx)
          });
        }
      }
    }

    const indexedSemanticGroups = (predicate) => directGroups
      .filter(predicate)
      .map((element) => {
        const rawIndex = element.getAttribute('data-gbdraw-record-index');
        const index = rawIndex === null || String(rawIndex).trim() === ''
          ? Number.NaN
          : Number(rawIndex);
        return Number.isInteger(index) && index >= 0 ? { element, index } : null;
      })
      .filter(Boolean)
      .sort((left, right) => left.index - right.index);
    const definitionRoles = new Set(['record-definition', 'record-definition-row']);
    const recordGroups = indexedSemanticGroups((element) => (
      element.hasAttribute('data-gbdraw-record-id') &&
      !definitionRoles.has(String(element.getAttribute('data-gbdraw-role') || '').trim())
    ));
    const definitionGroups = indexedSemanticGroups((element) => (
      String(element.getAttribute('data-gbdraw-role') || '').trim() === 'record-definition'
    ));
    const screenBox = (element) => {
      const box = element.getBoundingClientRect();
      return {
        id: String(element.id || ''),
        left: Number(box.left),
        right: Number(box.right),
        top: Number(box.top),
        bottom: Number(box.bottom)
      };
    };
    const overlap = (left, right) => ({
      overlapX: Math.min(left.right, right.right) - Math.max(left.left, right.left),
      overlapY: Math.min(left.bottom, right.bottom) - Math.max(left.top, right.top)
    });
    const domChecks = [];
    const addDomCheck = (label, leftElement, rightElement) => {
      if (!leftElement || !rightElement) {
        domChecks.push({ label, missing: true });
        return;
      }
      const left = screenBox(leftElement);
      const right = screenBox(rightElement);
      domChecks.push({ label, missing: false, ...overlap(left, right) });
    };
    for (let index = 0; index < records.length; index += 1) {
      addDomCheck(
        `definition-body:${index}`,
        definitionGroups[index]?.element,
        recordGroups[index]?.element
      );
      if (index >= records.length - 1) continue;
      addDomCheck(
        `body-body:${index}`,
        recordGroups[index]?.element,
        recordGroups[index + 1]?.element
      );
      addDomCheck(
        `definition-definition:${index}`,
        definitionGroups[index]?.element,
        definitionGroups[index + 1]?.element
      );
      addDomCheck(
        `comparison-upper-body:${index}`,
        comparisonGroups[index],
        recordGroups[index]?.element
      );
      addDomCheck(
        `comparison-lower-body:${index}`,
        comparisonGroups[index],
        recordGroups[index + 1]?.element
      );
    }

    return {
      error: '',
      schema: geometry.schema,
      mode: geometry.mode,
      recordCount: records.length,
      comparisonCount: comparisonGroups.length,
      comparisonGroupDiagnostics: Array.from(
        svg.querySelectorAll('g[id^="comparison"]')
      ).map((element) => ({
        id: String(element.id || ''),
        attributes: element.getAttributeNames(),
        queryRecordIndex: element.getAttribute('data-query-record-index'),
        subjectRecordIndex: element.getAttribute('data-subject-record-index'),
        queryRow: element.getAttribute('data-query-row'),
        subjectRow: element.getAttribute('data-subject-row'),
        childCount: element.childElementCount
      })),
      resultComparisonGroup: String(
        state.results.value[state.selectedResultIndex.value]?.content || ''
      ).match(/<g[^>]*id="comparison1"[^>]*>/)?.[0] || '',
      recordGroupCount: recordGroups.length,
      definitionGroupCount: definitionGroups.length,
      axisGaps,
      collisionChecks,
      collisionDomains: Array.from(new Set(
        collisionChecks.map((check) => check.kindPair)
      )).sort(),
      domChecks,
      geometrySignature: JSON.stringify(geometry)
    };
  });
};

const assertLayout = (layout) => {
  expect(layout.error || '').toBe('');
  // Browser run metadata publishes collision bands in the schema-1 wrapper.
  // The canvas-private schema-2 axisGapConstraints are covered by Python unit tests.
  expect(layout.schema).toBe(1);
  expect(layout.mode).toBe('linear');
  expect(layout.recordCount).toBe(5);
  expect(
    layout.comparisonCount,
    JSON.stringify({
      groups: layout.comparisonGroupDiagnostics,
      resultGroup: layout.resultComparisonGroup
    }, null, 2)
  ).toBe(4);
  expect(layout.recordGroupCount).toBe(5);
  expect(layout.definitionGroupCount).toBe(5);
  expect(layout.axisGaps).toHaveLength(4);
  expect(layout.axisGaps.every((gap) => gap > 0)).toBe(true);

  expect(layout.collisionChecks.length).toBeGreaterThan(0);
  expect(layout.collisionChecks.filter((check) => check.gapPx < -1e-6)).toEqual([]);
  expect(layout.collisionDomains).toEqual(expect.arrayContaining([
    'body:body',
    'comparison:comparison',
    'definition:definition'
  ]));

  expect(layout.domChecks.length).toBeGreaterThan(0);
  expect(layout.domChecks.filter((check) => (
    check.missing ||
    (check.overlapX > 1e-6 && check.overlapY > 1e-6)
  ))).toEqual([]);
};

const assertSvgGeometryParity = async (page, exportedPath) => {
  const exportedSvg = readFileSync(exportedPath, 'utf8');
  const parity = await page.evaluate(async (svgText) => {
    const { state } = await import('/gbdraw/web/js/state.js');
    const preview = state.svgContainer.value?.querySelector('svg');
    const parsed = new DOMParser().parseFromString(String(svgText || ''), 'image/svg+xml');
    const exported = parsed.documentElement;
    if (
      !preview ||
      !exported ||
      exported.tagName?.toLowerCase() !== 'svg' ||
      parsed.querySelector('parsererror')
    ) {
      return { error: 'The preview or exported SVG could not be parsed.' };
    }
    const geometryAttributes = [
      'viewBox', 'width', 'height',
      'transform', 'd', 'points', 'x', 'y', 'x1', 'x2', 'y1', 'y2',
      'cx', 'cy', 'r', 'rx', 'ry', 'dx', 'dy', 'font-family', 'font-size',
      'font-weight', 'font-style', 'text-anchor', 'dominant-baseline',
      'stroke-width', 'display', 'href', 'xlink:href'
    ];
    const signature = (root) => [root, ...root.querySelectorAll('*')].map((element) => {
      const tag = element.tagName.toLowerCase();
      const attributes = geometryAttributes
        .filter((name) => element.hasAttribute(name))
        .map((name) => [name, element.getAttribute(name)]);
      const text = tag === 'text' || tag === 'tspan' ? element.textContent : '';
      return [tag, String(element.id || ''), attributes, text];
    });
    const previewSignature = signature(preview);
    const exportedSignature = signature(exported);
    const length = Math.max(previewSignature.length, exportedSignature.length);
    let mismatchIndex = -1;
    for (let index = 0; index < length; index += 1) {
      if (
        JSON.stringify(previewSignature[index] ?? null) !==
        JSON.stringify(exportedSignature[index] ?? null)
      ) {
        mismatchIndex = index;
        break;
      }
    }
    return {
      error: '',
      equal: mismatchIndex === -1,
      previewElementCount: previewSignature.length,
      exportedElementCount: exportedSignature.length,
      mismatchIndex,
      previewMismatch: mismatchIndex < 0 ? null : previewSignature[mismatchIndex],
      exportedMismatch: mismatchIndex < 0 ? null : exportedSignature[mismatchIndex]
    };
  }, exportedSvg);
  expect(parity.error || '').toBe('');
  expect(parity.previewElementCount).toBeGreaterThan(0);
  expect(parity.previewElementCount).toBe(parity.exportedElementCount);
  expect(parity, JSON.stringify(parity, null, 2)).toMatchObject({
    equal: true,
    mismatchIndex: -1
  });
};

test.beforeAll(async () => {
  savedSessionDir = mkdtempSync(join(tmpdir(), 'gbdraw-losat-migration-'));
  expect(existsSync(fixturePath)).toBe(true);
  const fixtureBytes = expandedSessionBytes(fixturePath);
  const fixture = JSON.parse(fixtureBytes.toString('utf8'));
  sourceSession = fixture;
  expect(createHash('sha256').update(fixtureBytes).digest('hex')).toBe(expected.sourceSha256);
  expect(fixture.version).toBe(expected.sessionVersion);
  expect(fixture.renderRequest?.schema).toBe(expected.renderRequestSchema);
  expect(fixture.losatCache?.entries).toHaveLength(expected.storedRawEntries);
});

test.afterAll(async () => {
  rmSync(savedSessionDir, { recursive: true, force: true });
});

test('SVG sanitizer preserves linear comparison row hooks', async ({ page }) => {
  await page.goto('/gbdraw/web/index.html', {
    waitUntil: 'domcontentloaded'
  });
  const sanitized = await page.evaluate(async () => {
    const { sanitizeSvgContent } = await import(
      '/gbdraw/web/js/services/svg-sanitization.js'
    );
    const once = sanitizeSvgContent(
      '<svg xmlns="http://www.w3.org/2000/svg">'
        + '<g id="comparison1" data-query-row="0" data-subject-row="1">'
        + '<path d="M0,0 L1,1" /></g></svg>'
    );
    return sanitizeSvgContent(once);
  });
  expect(sanitized).toContain('data-query-row="0"');
  expect(sanitized).toContain('data-subject-row="1"');
});

test('legacy protein caches migrate, export readable TSV, and preserve uploads', async ({ page }) => {
  test.setTimeout(900000);
  await page.addInitScript(() => {
    window.__GBDRAW_LOSAT_EXECUTOR_CALLS__ = 0;
    window.__GBDRAW_LOSAT_EXECUTOR__ = async () => {
      window.__GBDRAW_LOSAT_EXECUTOR_CALLS__ += 1;
      throw new Error('LOSAT executor must not run during verified cache migration.');
    };
  });

  await page.goto('/gbdraw/web/index.html', { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await importSession(page, fixturePath);

  const renamed = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    for (let index = 0; index < app.linearSeqs.length; index += 1) {
      const source = app.linearSeqs[index].gb;
      app.linearSeqs[index].gb = new File(
        [await source.arrayBuffer()],
        `renamed-cache-input-${index + 1}.gbk`,
        { type: source.type || 'text/plain', lastModified: 0 }
      );
    }
    return app.linearSeqs.map((record) => ({
      name: record.gb.name,
      lastModified: record.gb.lastModified
    }));
  });
  expect(renamed).toEqual(Array.from({ length: 5 }, (_, index) => ({
    name: `renamed-cache-input-${index + 1}.gbk`,
    lastModified: 0
  })));
  const legacyUiBeforeCancel = await migrationUiSnapshot(page);
  const canceledLegacyRun = await cancelDuringRender(page);
  await settleAppRender(page);
  expect(canceledLegacyRun).toMatchObject({
    result: { status: 'canceled' },
    runRequestIssued: true,
    cancelInvoked: true,
    errorSummary: '',
    executorCalls: 0,
    authorityRestored: true
  });
  expect(await migrationUiSnapshot(page)).toEqual(legacyUiBeforeCancel);
  const legacyUiBeforeRenderError = await migrationUiSnapshot(page);
  const failedLegacyRun = await failRendererAfterMigration(page);
  await settleAppRender(page);
  expect(failedLegacyRun).toMatchObject({
    result: { status: 'error' },
    authorityRestored: true,
    executorCalls: 0,
    rendererFailureInjected: true
  });
  expect(failedLegacyRun.errorSummary).not.toBe('');
  expect(await migrationUiSnapshot(page)).toEqual(legacyUiBeforeRenderError);
  const migratedRun = await generateWithTelemetry(page);
  await settleAppRender(page);
  assertExpectedTelemetry(migratedRun);
  const firstLayout = await inspectLayout(page);
  assertLayout(firstLayout);
  await assertSvgGeometryParity(page, await downloadSvg(page));
  await assertHydratedDownloads(page);
  await assertHydratedUtf8PromptBoundary(page);

  const migratedPath = await saveSession(page);
  const migratedSession = readSession(migratedPath);
  assertCurrentProteinArtifacts(migratedSession);
  assertRenamedZeroMtimeResources(migratedSession);

  await page.reload({ waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await importSession(page, migratedPath);
  const regeneratedRun = await generateWithTelemetry(page);
  await settleAppRender(page);
  assertExpectedTelemetry(regeneratedRun);
  expect(regeneratedRun.telemetry?.proteinDerivedPayloadCacheMisses).toBeGreaterThan(0);
  expect(await page.evaluate(async () => {
    const { state } = await import('/gbdraw/web/js/state.js');
    return state.losatDerivedCache.value.size;
  })).toBeGreaterThan(0);
  const secondLayout = await inspectLayout(page);
  assertLayout(secondLayout);
  expect(secondLayout.geometrySignature).toBe(firstLayout.geometrySignature);

  const expectedUploadBytes = await installUserUploadedTsv(page);
  const uploadedSessionPath = await saveSession(page);
  expect(readSession(uploadedSessionPath).losatDerivedCache?.entries || []).toEqual([]);
  await page.reload({ waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await importSession(page, uploadedSessionPath);
  const restoredUpload = await readFirstUploadedTsv(page);
  expect(restoredUpload).toEqual({
    name: expected.downloads.userUploadedFilename,
    bytes: expectedUploadBytes
  });
});
