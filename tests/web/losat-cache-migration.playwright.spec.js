const { test, expect } = require('@playwright/test');
const { createHash } = require('node:crypto');
const { createReadStream, existsSync, readFileSync } = require('node:fs');
const { createServer } = require('node:http');
const { extname, join, normalize, resolve, sep } = require('node:path');
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

const contentTypes = {
  '.css': 'text/css; charset=utf-8',
  '.data': 'application/octet-stream',
  '.html': 'text/html; charset=utf-8',
  '.js': 'text/javascript; charset=utf-8',
  '.json': 'application/json; charset=utf-8',
  '.mjs': 'text/javascript; charset=utf-8',
  '.svg': 'image/svg+xml',
  '.wasm': 'application/wasm',
  '.whl': 'application/octet-stream',
  '.woff2': 'font/woff2'
};

let server;
let baseUrl;

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
  await page.evaluate(async () => window.__GBDRAW_APP__.saveSessionWithTitle());
  const download = await downloadPromise;
  const path = await download.path();
  expect(path).toBeTruthy();
  return path;
};

const assertLegacyCandidatesPreserved = (session) => {
  const candidates = session.legacyArtifacts?.proteinRawCandidates?.entries;
  expect(session.format).toBe('gbdraw-session');
  expect(session.losatCache?.entries || []).toHaveLength(0);
  expect(candidates).toHaveLength(expected.storedRawEntries);
  expect(candidates.every((candidate) => (
    candidate.state === 'pending' &&
    candidate.originalEntry?.schema === 2 &&
    candidate.originalEntry?.program === 'blastp'
  ))).toBe(true);
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
  const manifest = session.proteinIdentityManifest;
  const entries = session.losatCache?.entries || [];
  expect(manifest?.schema).toBe(1);
  expect(entries).toHaveLength(expected.totalPairs);
  expect(entries.every((entry) => (
    entry?.schema === 3 &&
    entry?.kind === 'raw-losat' &&
    entry?.identityKind === 'protein' &&
    entry?.program === 'blastp'
  ))).toBe(true);

  const transportOwners = new Map();
  for (const [instanceKey, instance] of Object.entries(manifest.recordInstances || {})) {
    expect(manifest.recordAnalyses?.[instance.recordAnalysisId]).toBeTruthy();
    expect(instance.bindingHash).toBeTruthy();
    for (const [featureId, transportId] of Object.entries(instance.transportIds || {})) {
      expect(featureId).toMatch(/^f_/);
      expect(transportId).not.toMatch(/^p_r_/);
      expect(transportId).toMatch(/@.+\|.+~f_/);
      expect(transportId).not.toMatch(/\s/);
      expect(transportOwners.has(transportId)).toBe(false);
      transportOwners.set(transportId, { instanceKey, featureId });
    }
  }

  let resolvedReferenceCount = 0;
  for (const entry of entries) {
    const query = manifest.recordInstances?.[entry.queryRecordInstanceKey];
    const subject = manifest.recordInstances?.[entry.subjectRecordInstanceKey];
    expect(query?.bindingHash).toBe(entry.queryBindingHash);
    expect(subject?.bindingHash).toBe(entry.subjectBindingHash);
    expect(manifest.recordAnalyses?.[query.recordAnalysisId]?.proteinSetHash)
      .toBe(entry.queryProteinSetHash);
    expect(manifest.recordAnalyses?.[subject.recordAnalysisId]?.proteinSetHash)
      .toBe(entry.subjectProteinSetHash);
    const queryIds = new Set(Object.values(query.transportIds || {}));
    const subjectIds = new Set(Object.values(subject.transportIds || {}));
    for (const rawLine of String(entry.text || '').split(/\r?\n/)) {
      const line = rawLine.trim();
      if (!line || line.startsWith('#')) continue;
      const [queryId, subjectId] = rawLine.split('\t');
      expect(queryId).not.toMatch(/^p_r_/);
      expect(subjectId).not.toMatch(/^p_r_/);
      expect(queryIds.has(queryId)).toBe(true);
      expect(subjectIds.has(subjectId)).toBe(true);
      expect(transportOwners.get(queryId)?.instanceKey).toBe(entry.queryRecordInstanceKey);
      expect(transportOwners.get(subjectId)?.instanceKey).toBe(entry.subjectRecordInstanceKey);
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

const assertExpectedTelemetry = (run) => {
  const diagnostics = `Generation diagnostics:\n${JSON.stringify(run, null, 2)}`;
  expect(run.result, diagnostics).toEqual({ status: 'ok' });
  expect(run.errorSummary, diagnostics).toBe('');
  expect(run.errorDetails, diagnostics).toEqual([]);
  expect(run.telemetry, diagnostics).toMatchObject({
    totalPairs: expected.totalPairs,
    cacheHits: expected.cacheHits,
    cacheMisses: expected.cacheMisses,
    uniqueJobs: expected.uniqueJobs,
    workerCalls: expected.workerCalls
  });
  expect(run.executorCalls, diagnostics).toBe(expected.workerCalls);
};

test.beforeAll(async () => {
  expect(existsSync(fixturePath)).toBe(true);
  const fixtureBytes = expandedSessionBytes(fixturePath);
  const fixture = JSON.parse(fixtureBytes.toString('utf8'));
  expect(createHash('sha256').update(fixtureBytes).digest('hex')).toBe(expected.sourceSha256);
  expect(fixture.version).toBe(expected.sessionVersion);
  expect(fixture.renderRequest?.schema).toBe(expected.renderRequestSchema);
  expect(fixture.losatCache?.entries).toHaveLength(expected.storedRawEntries);
  await new Promise((resolveServer, rejectServer) => {
    server = createServer((request, response) => {
      const url = new URL(request.url || '/', 'http://127.0.0.1');
      const requested = normalize(decodeURIComponent(url.pathname))
        .replace(/^(\.\.(?:\/|\\|$))+/, '');
      const filePath = resolve(repoRoot, requested.replace(/^[/\\]+/, ''));
      if ((!filePath.startsWith(`${repoRoot}${sep}`) && filePath !== repoRoot) || !existsSync(filePath)) {
        response.writeHead(404);
        response.end('Not found');
        return;
      }
      response.writeHead(200, {
        'Content-Type': contentTypes[extname(filePath)] || 'application/octet-stream'
      });
      createReadStream(filePath).pipe(response);
    });
    server.once('error', rejectServer);
    server.listen(0, '127.0.0.1', () => {
      baseUrl = `http://127.0.0.1:${server.address().port}`;
      resolveServer();
    });
  });
});

test.afterAll(async () => {
  await new Promise((resolveClose) => server.close(resolveClose));
});

test('real schema-2 protein cache survives save/load and migrates without LOSAT work', async ({ page }) => {
  test.setTimeout(600000);
  await page.addInitScript(() => {
    window.__GBDRAW_LOSAT_EXECUTOR_CALLS__ = 0;
    window.__GBDRAW_LOSAT_EXECUTOR__ = async () => {
      window.__GBDRAW_LOSAT_EXECUTOR_CALLS__ += 1;
      throw new Error('LOSAT executor must not run during verified cache migration.');
    };
  });

  await page.goto(`${baseUrl}/gbdraw/web/index.html`, { waitUntil: 'domcontentloaded' });
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

  const preGeneratePath = await saveSession(page);
  const preGenerateSession = readSession(preGeneratePath);
  assertLegacyCandidatesPreserved(preGenerateSession);
  assertRenamedZeroMtimeResources(preGenerateSession);

  await page.reload({ waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await importSession(page, preGeneratePath);
  await page.waitForFunction(
    () => window.__GBDRAW_APP__?.pyodideReady === true,
    null,
    { timeout: 240000 }
  );
  assertExpectedTelemetry(await generateWithTelemetry(page));

  const migratedPath = await saveSession(page);
  const migratedSession = readSession(migratedPath);
  assertCurrentProteinArtifacts(migratedSession);
  assertRenamedZeroMtimeResources(migratedSession);

  await page.reload({ waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await importSession(page, migratedPath);
  await page.waitForFunction(
    () => window.__GBDRAW_APP__?.pyodideReady === true,
    null,
    { timeout: 240000 }
  );
  assertExpectedTelemetry(await generateWithTelemetry(page));
});
