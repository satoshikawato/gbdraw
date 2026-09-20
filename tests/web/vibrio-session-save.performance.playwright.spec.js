const { test, expect } = require('@playwright/test');
const { createHash } = require('node:crypto');
const { spawnSync } = require('node:child_process');
const { existsSync, mkdirSync, readFileSync } = require('node:fs');
const os = require('node:os');
const { join, resolve } = require('node:path');
const {
  getDiagramWorkerActivity,
  openApp
} = require('./helpers/app-lifecycle.cjs');

const repoRoot = resolve(process.env.GBDRAW_REPO || process.cwd());
const fixturePath = join(
  repoRoot,
  'gbdraw',
  'web',
  'gallery',
  'sessions',
  'vibrio-harveyi-group-collinear.gbdraw-session.json.gz'
);
const fixtureSha = createHash('sha256').update(readFileSync(fixturePath)).digest('hex');
const OPERATION_TIMEOUT_MS = 900_000;
const BASE_SAVE_WALL_MS = 22_341;
const BASE_HEAP_DELTA_BYTES = 3_237_405_402 - 1_339_279_560;
const BASE_HEARTBEAT_GAP_MS = 14_379;

test.use({
  launchOptions: {
    args: ['--enable-precise-memory-info']
  }
});

const installSaveAcceptanceProbe = (page) => page.addInitScript(() => {
  const lifecycle = [];
  const responsiveness = {
    heartbeatCount: 0,
    maximumHeartbeatGapMs: 0,
    memorySupported: Boolean(performance.memory),
    memoryHighWaterBytes: Number(performance.memory?.usedJSHeapSize || 0)
  };
  let lastHeartbeat = performance.now();

  window.__GBDRAW_SAVE_ACCEPTANCE__ = {
    lifecycle,
    responsiveness,
    reset() {
      lifecycle.length = 0;
      responsiveness.heartbeatCount = 0;
      responsiveness.maximumHeartbeatGapMs = 0;
      responsiveness.memoryHighWaterBytes = Number(performance.memory?.usedJSHeapSize || 0);
      lastHeartbeat = performance.now();
    },
    snapshot() {
      return {
        lifecycle: lifecycle.map((event) => ({ ...event })),
        responsiveness: { ...responsiveness },
        usedJsHeapBytes: Number(performance.memory?.usedJSHeapSize || 0)
      };
    }
  };
  window.__GBDRAW_TEST_HOOKS__ = {
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
});

const sessionSummaryScript = String.raw`
const { createHash } = require('node:crypto');
const { readFileSync } = require('node:fs');
const { gunzipSync } = require('node:zlib');
const path = process.argv[1];
const bytes = readFileSync(path);
const text = bytes[0] === 0x1f && bytes[1] === 0x8b
  ? gunzipSync(bytes).toString('utf8')
  : bytes.toString('utf8');
const document = JSON.parse(text);
const digest = (value) => createHash('sha256').update(JSON.stringify(value)).digest('hex');
const catalog = document.editorState?.featureCatalog;
const losatEntries = document.losatCache?.entries || [];
process.stdout.write(JSON.stringify({
  format: document.format,
  version: document.version,
  requestSchema: document.renderRequest?.schema,
  topLevelKeys: Object.keys(document).sort(),
  resourceCount: Object.keys(document.resources || {}).length,
  resultCount: (document.results || []).length,
  catalogItems: (catalog?.items || []).length,
  losatEntries: losatEntries.length,
  hashes: {
    renderRequest: digest(document.renderRequest),
    resources: digest(document.resources),
    webFiles: digest(document.webFiles),
    results: digest(document.results),
    featureCatalog: digest(catalog),
    losatRawTextAuthority: digest(losatEntries
      .map((entry) => [String(entry.key || ''), String(entry.text || '')])
      .sort(([left], [right]) => left.localeCompare(right))),
    proteinIdentityManifest: digest(document.proteinIdentityManifest)
  }
}));
`;

const extractResultSvgScript = String.raw`
const { readFileSync, writeFileSync } = require('node:fs');
const { gunzipSync } = require('node:zlib');
const [sessionPath, outputPath] = process.argv.slice(1);
const bytes = readFileSync(sessionPath);
const text = bytes[0] === 0x1f && bytes[1] === 0x8b
  ? gunzipSync(bytes).toString('utf8')
  : bytes.toString('utf8');
const document = JSON.parse(text);
const content = document.results?.[0]?.content;
if (typeof content !== 'string' || !/<svg[\s>]/.test(content)) {
  throw new Error('Session does not contain an SVG result.');
}
writeFileSync(outputPath, content, 'utf8');
`;

const inspectSession = (path) => {
  const inspected = spawnSync(
    process.execPath,
    ['--max-old-space-size=3072', '-e', sessionSummaryScript, path],
    {
      cwd: repoRoot,
      encoding: 'utf8',
      maxBuffer: 4 * 1024 * 1024,
      timeout: OPERATION_TIMEOUT_MS
    }
  );
  expect(inspected.status, `${inspected.stdout}\n${inspected.stderr}`).toBe(0);
  return JSON.parse(inspected.stdout);
};

const extractResultSvg = (sessionPath, outputPath) => {
  const extracted = spawnSync(
    process.execPath,
    ['--max-old-space-size=3072', '-e', extractResultSvgScript, sessionPath, outputPath],
    {
      cwd: repoRoot,
      encoding: 'utf8',
      timeout: OPERATION_TIMEOUT_MS
    }
  );
  expect(extracted.status, `${extracted.stdout}\n${extracted.stderr}`).toBe(0);
};

const compareResultSvgs = (expectedPath, actualPath) => {
  const comparisonScript = [
    'import sys',
    'from tests.utils.svg_compare import compare_svgs',
    'result = compare_svgs(sys.argv[1], sys.argv[2])',
    'print(result.message)',
    'print("\\n".join(result.differences))',
    'raise SystemExit(0 if result.equal else 1)'
  ].join(';');
  const compared = spawnSync(
    process.env.GBDRAW_PYTHON || 'python',
    ['-c', comparisonScript, expectedPath, actualPath],
    { cwd: repoRoot, encoding: 'utf8', timeout: OPERATION_TIMEOUT_MS }
  );
  expect(compared.status, `${compared.stdout}\n${compared.stderr}`).toBe(0);
};

const crossSurfaceAcceptance = (sessionPath, outputDirectory) => {
  const python = process.env.GBDRAW_PYTHON || 'python';
  mkdirSync(outputDirectory, { recursive: true });
  const readerScript = [
    'import json, sys',
    'from pathlib import Path',
    'from gbdraw.api import load_session_document, materialize_session, session_to_request',
    'document = load_session_document(Path(sys.argv[1]))',
    'with materialize_session(document, output_directory=Path(sys.argv[2])) as materialized:',
    '    request = session_to_request(materialized)',
    "print(json.dumps({'version': document.version, 'mode': document.mode, 'records': len(request.records)}))"
  ].join('\n');
  const reader = spawnSync(python, ['-c', readerScript, sessionPath, outputDirectory], {
    cwd: repoRoot,
    encoding: 'utf8',
    env: { ...process.env, PYTHONPATH: repoRoot },
    timeout: OPERATION_TIMEOUT_MS
  });
  expect(reader.status, `${reader.stdout}\n${reader.stderr}`).toBe(0);

  const outputPrefix = join(outputDirectory, 'vibrio-cli-replay');
  const cli = spawnSync(python, [
    '-m',
    'gbdraw.cli',
    'linear',
    '--session',
    sessionPath,
    '-o',
    outputPrefix,
    '-f',
    'svg',
    '--overwrite'
  ], {
    cwd: outputDirectory,
    encoding: 'utf8',
    env: { ...process.env, PYTHONPATH: repoRoot },
    timeout: OPERATION_TIMEOUT_MS
  });
  expect(cli.status, `${cli.stdout}\n${cli.stderr}`).toBe(0);
  expect(existsSync(`${outputPrefix}.svg`)).toBe(true);
  return {
    reader: JSON.parse(reader.stdout.trim()),
    cliExitCode: cli.status,
    cliSvgBytes: readFileSync(`${outputPrefix}.svg`).byteLength
  };
};

const snapshotLiveState = (page) => page.evaluate(async () => {
  const { state } = await import('/gbdraw/web/js/state.js');
  const app = window.__GBDRAW_APP__;
  const history = window.__GBDRAW_HISTORY__;
  window.__GBDRAW_SAVE_LIVE_STATE__ = {
    results: [...app.results],
    selectedResultIndex: app.selectedResultIndex,
    linearFiles: state.linearSeqs.map((sequence) => sequence.gb),
    losatCache: state.losatCache.value,
    losatEntries: [...state.losatCache.value],
    manifest: state.proteinIdentityManifest.value,
    undoCount: history.getUndoCount(),
    redoCount: history.getRedoCount()
  };
});

const liveStateIntact = (page) => page.evaluate(async () => {
  const { state } = await import('/gbdraw/web/js/state.js');
  const app = window.__GBDRAW_APP__;
  const history = window.__GBDRAW_HISTORY__;
  const before = window.__GBDRAW_SAVE_LIVE_STATE__;
  return {
    results: before.results.length === app.results.length
      && before.results.every((result, index) => result === app.results[index]),
    selectedResultIndex: before.selectedResultIndex === app.selectedResultIndex,
    linearFiles: before.linearFiles.length === state.linearSeqs.length
      && before.linearFiles.every((file, index) => file === state.linearSeqs[index].gb),
    losatCache: before.losatCache === state.losatCache.value,
    losatEntries: before.losatEntries.length === state.losatCache.value.size
      && before.losatEntries.every(([key, value]) => state.losatCache.value.get(key) === value),
    manifest: before.manifest === state.proteinIdentityManifest.value,
    undoCount: before.undoCount === history.getUndoCount(),
    redoCount: before.redoCount === history.getRedoCount()
  };
});

test.describe.configure({ mode: 'serial' });

test('Vibrio Session saves once within memory, responsiveness, and compatibility budgets', async ({
  page,
  context,
  browser
}, testInfo) => {
  test.setTimeout(1_800_000);
  const terminal = { pageErrors: [], crashes: 0 };
  const requests = [];
  const dialogs = [];
  let downloadCount = 0;
  page.on('pageerror', (error) => terminal.pageErrors.push(String(error?.message || error)));
  page.on('crash', () => { terminal.crashes += 1; });
  page.on('request', (request) => requests.push(request.url()));
  page.on('download', () => { downloadCount += 1; });
  page.on('dialog', async (dialog) => {
    dialogs.push({ type: dialog.type(), message: dialog.message() });
    await dialog.accept();
  });

  await installSaveAcceptanceProbe(page);
  await openApp(page);
  await page.locator(
    'input[type="file"][accept*="application/json"][accept*="application/gzip"]'
  ).setInputFiles(fixturePath);
  await page.waitForFunction(() => (
    window.__GBDRAW_APP__?.sessionImportPending === false
    && window.__GBDRAW_APP__?.results?.length === 1
    && document.querySelector('.shadow-xl.origin-top > svg')
  ), null, { timeout: OPERATION_TIMEOUT_MS });
  await snapshotLiveState(page);

  const workerBeforeSave = await getDiagramWorkerActivity(page);
  expect(workerBeforeSave.constructions).toBe(0);
  expect(workerBeforeSave.initializations).toBe(0);
  await page.evaluate(() => {
    window.__GBDRAW_APP__.sessionTitle = 'issue-544-vibrio-save';
    window.__GBDRAW_SAVE_ACCEPTANCE__.reset();
  });
  const before = await page.evaluate(() => window.__GBDRAW_SAVE_ACCEPTANCE__.snapshot());
  const downloadPromise = page.waitForEvent('download', { timeout: OPERATION_TIMEOUT_MS });
  const startedAt = Date.now();
  await page.evaluate(() => {
    window.__GBDRAW_VIBRIO_SAVE__ = window.__GBDRAW_APP__.saveSessionWithTitle();
  });
  await page.waitForFunction(() => window.__GBDRAW_SAVE_ACCEPTANCE__.lifecycle.some(
    (event) => event.name === 'session-save-paint-opportunity-completed'
  ), null, { timeout: OPERATION_TIMEOUT_MS });
  const saveButton = page.getByRole('button', { name: 'Save Session', exact: true });
  await expect(saveButton).toBeDisabled();
  await expect(saveButton).toHaveAttribute('aria-busy', 'true');
  await expect(page.locator('[data-session-save-status]')).toHaveText('Saving session…');

  const outcome = await page.evaluate(() => window.__GBDRAW_VIBRIO_SAVE__);
  const saveWallMs = Date.now() - startedAt;
  expect(outcome.status).toBe('saved');
  const download = await downloadPromise;
  const savedPath = testInfo.outputPath('issue-544-vibrio-save.gbdraw-session.json.gz');
  await download.saveAs(savedPath);
  const after = await page.evaluate(() => window.__GBDRAW_SAVE_ACCEPTANCE__.snapshot());

  await expect(saveButton).toBeEnabled();
  await expect(saveButton).toHaveAttribute('aria-busy', 'false');
  await expect(page.locator('[data-session-save-status]')).toHaveCount(0);
  expect(await liveStateIntact(page)).toEqual({
    results: true,
    selectedResultIndex: true,
    linearFiles: true,
    losatCache: true,
    losatEntries: true,
    manifest: true,
    undoCount: true,
    redoCount: true
  });

  const lifecycleNames = after.lifecycle.map(({ name }) => name);
  for (const name of [
    'session-save-pending-published',
    'session-save-paint-opportunity-completed',
    'session-save-projection-start',
    'session-save-projection-end',
    'session-save-compression-start',
    'session-save-compression-end',
    'session-save-download-handoff-completed',
    'session-save-pending-cleared'
  ]) {
    expect(lifecycleNames.filter((eventName) => eventName === name), name).toHaveLength(1);
  }
  const orderedLifecycle = [
    'session-save-pending-published',
    'session-save-paint-opportunity-completed',
    'session-save-projection-start',
    'session-save-projection-end',
    'session-save-compression-start',
    'session-save-compression-end',
    'session-save-download-handoff-completed',
    'session-save-pending-cleared'
  ].map((name) => lifecycleNames.indexOf(name));
  expect(orderedLifecycle).toEqual([...orderedLifecycle].sort((left, right) => left - right));
  expect(downloadCount).toBe(1);
  expect(dialogs.filter(({ message }) => message.startsWith('Compressed session size is ')))
    .toHaveLength(1);
  expect(after.lifecycle.find(
    ({ name }) => name === 'session-save-catalog-preparation-end'
  )?.reusedCommittedSession).toBe(true);
  expect(lifecycleNames).not.toContain('catalog.admission-started');

  const heapHighWaterBytes = Math.max(
    after.usedJsHeapBytes,
    after.responsiveness.memoryHighWaterBytes
  );
  const heapDeltaBytes = heapHighWaterBytes - before.usedJsHeapBytes;
  const performanceEvidence = {
    singleObservation: true,
    saveWallMs,
    usedJsHeapBeforeBytes: before.usedJsHeapBytes,
    heapHighWaterBytes,
    heapDeltaBytes,
    maximumHeartbeatGapMs: after.responsiveness.maximumHeartbeatGapMs,
    compressedBytes: outcome.blob?.size || readFileSync(savedPath).byteLength,
    lifecycle: after.lifecycle
  };
  console.log(`GBDRAW_ISSUE_544_PERFORMANCE ${JSON.stringify(performanceEvidence)}`);
  await testInfo.attach('issue-544-vibrio-save-performance.json', {
    body: Buffer.from(JSON.stringify(performanceEvidence, null, 2)),
    contentType: 'application/json'
  });
  if (after.responsiveness.memorySupported) {
    expect(heapDeltaBytes).toBeLessThan(BASE_HEAP_DELTA_BYTES);
  }
  expect(
    after.responsiveness.maximumHeartbeatGapMs < 1_000
      || after.responsiveness.maximumHeartbeatGapMs <= BASE_HEARTBEAT_GAP_MS * 0.2
  ).toBe(true);
  expect(saveWallMs).toBeLessThan(BASE_SAVE_WALL_MS);

  const externalRequests = requests.filter((url) => {
    const parsed = new URL(url);
    return !['http://127.0.0.1:4173', 'blob:', 'data:'].includes(parsed.origin)
      && !['blob:', 'data:'].includes(parsed.protocol);
  });
  expect(externalRequests).toEqual([]);
  expect(terminal).toEqual({ pageErrors: [], crashes: 0 });
  expect((await getDiagramWorkerActivity(page)).constructions).toBe(0);

  await context.close();
  const sourceSummary = inspectSession(fixturePath);
  const savedSummary = inspectSession(savedPath);
  expect(sourceSummary).toMatchObject({
    format: 'gbdraw-session',
    version: 41,
    requestSchema: 7,
    resourceCount: 12,
    resultCount: 1,
    catalogItems: 1,
    losatEntries: 59
  });
  expect(savedSummary).toMatchObject({
    format: 'gbdraw-session',
    version: 42,
    requestSchema: 7,
    resourceCount: 12,
    resultCount: 1,
    catalogItems: 1,
    losatEntries: 59
  });
  for (const authority of [
    'renderRequest',
    'resources',
    'featureCatalog',
    'losatRawTextAuthority',
    'proteinIdentityManifest'
  ]) {
    expect(savedSummary.hashes[authority], authority).toBe(sourceSummary.hashes[authority]);
  }
  const sourceSvgPath = testInfo.outputPath('source-result.svg');
  const savedSvgPath = testInfo.outputPath('saved-result.svg');
  extractResultSvg(fixturePath, sourceSvgPath);
  extractResultSvg(savedPath, savedSvgPath);
  compareResultSvgs(sourceSvgPath, savedSvgPath);

  const freshContext = await browser.newContext();
  const freshPage = await freshContext.newPage();
  const freshRequests = [];
  const freshErrors = [];
  freshPage.on('request', (request) => freshRequests.push(request.url()));
  freshPage.on('pageerror', (error) => freshErrors.push(String(error?.message || error)));
  freshPage.on('dialog', (dialog) => dialog.accept());
  await openApp(freshPage);
  await freshPage.locator(
    'input[type="file"][accept*="application/json"][accept*="application/gzip"]'
  ).setInputFiles(savedPath);
  await freshPage.waitForFunction(() => (
    window.__GBDRAW_APP__?.sessionImportPending === false
    && window.__GBDRAW_APP__?.results?.length === 1
    && document.querySelector('.shadow-xl.origin-top > svg')
  ), null, { timeout: OPERATION_TIMEOUT_MS });
  expect((await getDiagramWorkerActivity(freshPage)).constructions).toBe(0);
  await expect(freshPage.getByRole('button', { name: 'Save Session', exact: true })).toBeEnabled();
  await freshPage.getByRole('searchbox', { name: 'Search features' }).fill('dnaA');
  await freshPage.getByRole('button', { name: 'Search features', exact: true }).click();
  await expect(freshPage.getByRole('button', { name: 'Open active feature' })).toBeEnabled();
  const pairwiseMatch = freshPage.locator(
    '.shadow-xl.origin-top > svg [data-gbdraw-pairwise-match-id]'
  ).first();
  await expect(pairwiseMatch).toHaveAttribute('role', 'button');
  await pairwiseMatch.click();
  await expect(freshPage.getByRole('dialog', { name: 'Pairwise match details' })).toBeVisible();
  expect((await getDiagramWorkerActivity(freshPage)).constructions).toBe(0);
  expect(freshErrors).toEqual([]);
  expect(freshRequests.filter((url) => {
    const parsed = new URL(url);
    return !['http://127.0.0.1:4173', 'blob:', 'data:'].includes(parsed.origin)
      && !['blob:', 'data:'].includes(parsed.protocol);
  })).toEqual([]);
  await freshContext.close();

  const crossSurface = crossSurfaceAcceptance(savedPath, testInfo.outputPath('cross-surface'));
  expect(crossSurface.reader).toEqual({ version: 42, mode: 'linear', records: 5 });
  expect(crossSurface.cliExitCode).toBe(0);
  expect(crossSurface.cliSvgBytes).toBeGreaterThan(0);

  const evidence = {
    repositorySha: spawnSync('git', ['rev-parse', 'HEAD'], {
      cwd: repoRoot,
      encoding: 'utf8'
    }).stdout.trim(),
    fixtureSha,
    browser: {
      project: testInfo.project.name,
      version: browser.version(),
      platform: process.platform,
      osRelease: os.release(),
      viewport: testInfo.project.use.viewport,
      launchArgs: ['--enable-precise-memory-info']
    },
    performance: performanceEvidence,
    lifecycle: after.lifecycle,
    sourceSummary,
    savedSummary,
    crossSurface,
    downloadCount,
    terminal,
    externalRequestCount: externalRequests.length
  };
  console.log(`GBDRAW_ISSUE_544_EVIDENCE ${JSON.stringify(evidence)}`);
  await testInfo.attach('issue-544-vibrio-save-evidence.json', {
    body: Buffer.from(JSON.stringify(evidence, null, 2)),
    contentType: 'application/json'
  });
});
