const { test, expect } = require('@playwright/test');
const { createReadStream, existsSync, readFileSync } = require('node:fs');
const { createServer } = require('node:http');
const { extname, join, normalize, resolve, sep } = require('node:path');

const repoRoot = resolve(process.env.GBDRAW_REPO || process.cwd());
const contentTypes = {
  '.html': 'text/html; charset=utf-8',
  '.js': 'text/javascript; charset=utf-8',
  '.mjs': 'text/javascript; charset=utf-8',
  '.css': 'text/css; charset=utf-8',
  '.json': 'application/json; charset=utf-8',
  '.svg': 'image/svg+xml',
  '.wasm': 'application/wasm',
  '.whl': 'application/octet-stream'
};

let server;
let baseUrl;

const makeComparisonGenbank = (recordId, base = 'atg') => {
  const sequence = base.repeat(100);
  const origin = sequence.match(/.{1,60}/g).map((chunk, index) => {
    const groups = chunk.match(/.{1,10}/g).join(' ');
    return `${String(index * 60 + 1).padStart(9)} ${groups}`;
  }).join('\n');
  return `LOCUS       ${recordId.padEnd(24)} 300 bp    DNA     linear   UNA 01-JAN-2000
DEFINITION  linear comparison browser test.
ACCESSION   ${recordId}
VERSION     ${recordId}
KEYWORDS    .
SOURCE      synthetic construct
  ORGANISM  synthetic construct
            .
FEATURES             Location/Qualifiers
     CDS             1..90
                     /product="comparison test protein"
ORIGIN
${origin}
//
`;
};

const installDiagramRequestObserver = async (page) => {
  await page.addInitScript(() => {
    window.__GBDRAW_DIAGRAM_RUNS__ = [];
    window.__GBDRAW_RUNTIME_URLS__ = { fetches: [], workers: [] };
    const nativeFetch = window.fetch.bind(window);
    window.fetch = (...args) => {
      window.__GBDRAW_RUNTIME_URLS__.fetches.push(String(args[0]?.url || args[0] || ''));
      return nativeFetch(...args);
    };
    const NativeWorker = window.Worker;
    window.Worker = new Proxy(NativeWorker, {
      construct(target, args) {
        const worker = Reflect.construct(target, args, target);
        const workerUrl = String(args[0] || '');
        window.__GBDRAW_RUNTIME_URLS__.workers.push(workerUrl);
        if (workerUrl.includes('diagram-generation-worker.js')) {
          const nativePostMessage = worker.postMessage.bind(worker);
          worker.postMessage = (message, transfer) => {
            if (message?.type === 'run' && message?.payload?.request) {
              window.__GBDRAW_DIAGRAM_RUNS__.push(
                JSON.parse(JSON.stringify(message.payload.request))
              );
            }
            if (transfer === undefined) return nativePostMessage(message);
            return nativePostMessage(message, transfer);
          };
        }
        return worker;
      }
    });
  });
};

test.beforeAll(async () => {
  await new Promise((resolveServer, rejectServer) => {
    server = createServer((request, response) => {
      const url = new URL(request.url || '/', 'http://127.0.0.1');
      const requestedPath = normalize(decodeURIComponent(url.pathname)).replace(/^(\.\.(?:\/|\\|$))+/, '');
      const filePath = resolve(repoRoot, requestedPath.replace(/^[/\\]+/, ''));
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

test('Linear record rows and N-to-M comparison batches remain keyed by sequence uid', async ({ page }) => {
  await page.goto(`${baseUrl}/gbdraw/web/index.html`, { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);

  const state = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    app.mode = 'linear';
    app.addLinearSeq();
    app.addLinearSeq();
    app.addLinearSeq();
    app.files.linearCanonicalComparisons = [{ id: 'stale' }];
    app.losatCacheInfo = [
      { edgeKey: 'stale->edge', key: 'linear-cache' },
      { key: 'circular-cache' }
    ];
    app.setLinearRecordLayoutEnabled(true);
    app.setLinearRecordRow(app.linearSeqs[0].uid, 1);
    app.setLinearRecordRow(app.linearSeqs[1].uid, 1);
    app.setLinearRecordRow(app.linearSeqs[2].uid, 2);
    app.setLinearRecordRow(app.linearSeqs[3].uid, 2);
    app.addLinearComparisonBatch(true);
    return {
      tokens: app.linearLayoutTokens,
      comparisonMode: app.linearComparisonPlan.mode,
      comparisonCount: app.linearComparisonPlan.edges.filter((item) => item.included).length,
      endpoints: app.linearComparisonPlan.edges
        .filter((item) => item.included)
        .map((item) => [item.queryUid, item.subjectUid]),
      uniqueUids: new Set(app.linearSeqs.map((item) => item.uid)).size,
      canonicalComparisonCount: app.files.linearCanonicalComparisons.length,
      cacheInfo: app.losatCacheInfo.map((entry) => entry.key)
    };
  });

  expect(state.tokens).toEqual(['#1@1', '#2@1', '#3@2', '#4@2']);
  expect(state.comparisonMode).toBe('selected');
  expect(state.comparisonCount).toBe(4);
  expect(state.uniqueUids).toBe(4);
  expect(state.canonicalComparisonCount).toBe(0);
  expect(state.cacheInfo).toEqual(['circular-cache']);
  expect(new Set(state.endpoints.map((pair) => pair.join('->'))).size).toBe(4);
  await expect(page.locator('input[aria-label="Linear record row"]')).toHaveCount(4);
  await expect(page.getByText('All adjacent-row pairs (cross-product)', { exact: true })).toBeVisible();
  const globalComparisonControls = page.locator('[data-capture="linear-blast-source"]');
  await expect(globalComparisonControls.getByRole('radio', { name: 'No comparison' })).toBeVisible();
  await expect(globalComparisonControls.getByRole('radio', { name: 'Run LOSAT' })).toBeVisible();
  await expect(globalComparisonControls.getByRole('radio', { name: 'Upload BLAST TSV' })).toBeVisible();
  await globalComparisonControls.getByRole('radio', { name: 'No comparison' }).check();
  const optedOut = await page.evaluate(() => ({
    mode: window.__GBDRAW_APP__.linearComparisonPlan.mode,
    retainedDrafts: window.__GBDRAW_APP__.linearComparisonPlan.edges.length,
    resolvedEdges: window.__GBDRAW_APP__.linearComparisonResolution.edges.length
  }));
  expect(optedOut).toEqual({ mode: 'none', retainedDrafts: 4, resolvedEdges: 0 });
  await expect(page.locator('[data-capture="linear-losat-settings"]')).toHaveCount(0);
});

test('No comparison completes a real render without touching dormant comparison work', async ({ page }) => {
  test.setTimeout(300000);
  await installDiagramRequestObserver(page);
  await page.addInitScript(() => {
    window.__GBDRAW_LOSAT_EXECUTOR_CALLS__ = 0;
    window.__GBDRAW_LOSAT_EXECUTOR__ = async () => {
      window.__GBDRAW_LOSAT_EXECUTOR_CALLS__ += 1;
      throw new Error('No comparison must not execute LOSAT.');
    };
  });
  await page.goto(`${baseUrl}/gbdraw/web/index.html`, { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await page.evaluate(async (records) => {
    const app = window.__GBDRAW_APP__;
    const { state } = await import('./js/state.js');
    app.mode = 'linear';
    app.lInputType = 'gb';
    app.addLinearSeq();
    app.addLinearSeq();
    records.forEach((content, index) => app.setLinearSeqPrimaryFile(index, 'gb', new File(
      [content], `none-record-${index + 1}.gbk`, { type: 'text/plain', lastModified: index + 1 }
    )));
    Object.assign(app.form, {
      legend: 'bottom', show_gc: false, show_skew: false,
      show_depth: false, show_labels_linear: 'none'
    });
    const dormant = new File([
      'DormantA\tDormantB\t100\t30\t0\t0\t1\t30\t1\t30\t1e-20\t100\n'
    ], 'dormant-comparison.tsv', { type: 'text/tab-separated-values', lastModified: 99 });
    const nativeArrayBuffer = dormant.arrayBuffer.bind(dormant);
    window.__GBDRAW_DORMANT_READS__ = 0;
    dormant.arrayBuffer = () => {
      window.__GBDRAW_DORMANT_READS__ += 1;
      return nativeArrayBuffer();
    };
    const [first, second] = app.linearSeqs;
    app.linearComparisonPlan.edges.splice(0, app.linearComparisonPlan.edges.length, {
      id: 'dormant-none-edge', queryUid: first.uid, subjectUid: second.uid,
      included: true, fileActive: true, losatFilenameActive: true,
      source: 'upload', file: dormant, losatFilename: 'dormant-losat-name.tsv'
    });
    app.setLinearComparisonGlobalAction('none');
    const rawCache = new Map([['dormant-cache', { text: 'must remain unread' }]]);
    const nativeGet = rawCache.get.bind(rawCache);
    window.__GBDRAW_CACHE_LOOKUPS__ = 0;
    rawCache.get = (key) => {
      window.__GBDRAW_CACHE_LOOKUPS__ += 1;
      return nativeGet(key);
    };
    state.losatCache.value = window.Vue.markRaw(rawCache);
  }, [
    makeComparisonGenbank('NoneRecA', 'atg'),
    makeComparisonGenbank('NoneRecB', 'gct'),
    makeComparisonGenbank('NoneRecC', 'tta')
  ]);

  const outcome = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const { state } = await import('./js/state.js');
    const result = await app.runAnalysis();
    const request = window.__GBDRAW_DIAGRAM_RUNS__.at(-1);
    const svg = new DOMParser().parseFromString(app.results[0].content, 'image/svg+xml');
    return {
      result,
      error: app.errorLog,
      request: [request.records.length, request.comparisons.length],
      resolution: [Object.isFrozen(app.linearComparisonResolution), app.linearComparisonResolution.edges.length],
      draft: app.linearComparisonPlan.edges.map((edge) => [
        edge.id, edge.included, edge.fileActive, edge.losatFilenameActive,
        edge.file?.name, edge.losatFilename
      ]),
      skipped: [
        window.__GBDRAW_DORMANT_READS__, window.__GBDRAW_CACHE_LOOKUPS__,
        window.__GBDRAW_LOSAT_EXECUTOR_CALLS__
      ],
      cacheSize: state.losatCache.value.size,
      losatRuntimeUrls: [
        ...window.__GBDRAW_RUNTIME_URLS__.fetches,
        ...window.__GBDRAW_RUNTIME_URLS__.workers
      ].filter((url) => url.includes('losat')),
      comparisonSvgNodes: svg.querySelectorAll(
        '[data-query-row], [data-subject-row], [data-gbdraw-role="comparison-legend"], #pairwise_legend'
      ).length
    };
  });
  expect(outcome).toEqual({
    result: { status: 'ok' }, error: null, request: [3, 0], resolution: [true, 0],
    draft: [['dormant-none-edge', true, true, true, 'dormant-comparison.tsv', 'dormant-losat-name.tsv']],
    skipped: [0, 0, 0], cacheSize: 1, losatRuntimeUrls: [], comparisonSvgNodes: 0
  });
  await expect(page.getByText('Raw LOSAT results', { exact: true })).toHaveCount(0);
});

test('Sparse upload and mixed selected renders keep snapshots and raw cache identity stable', async ({ page }) => {
  test.setTimeout(420000);
  await installDiagramRequestObserver(page);
  await page.addInitScript(() => {
    window.__GBDRAW_LOSAT_EXECUTOR_CALLS__ = 0;
    window.__GBDRAW_LOSAT_EXECUTOR_JOBS__ = [];
    window.__GBDRAW_LOSAT_EXECUTOR__ = async (jobs) => {
      window.__GBDRAW_LOSAT_EXECUTOR_CALLS__ += 1;
      window.__GBDRAW_LOSAT_EXECUTOR_JOBS__.push(jobs.map((job) => [
        job.edgeKey, job.ordinal, job.queryIndex, job.subjectIndex, job.cacheKey
      ]));
      if (window.__GBDRAW_LOSAT_EXECUTOR_CALLS__ > 1) {
        throw new Error('A content-addressed cache hit must not execute another LOSAT job.');
      }
      await new Promise((resolve) => { window.__GBDRAW_RELEASE_LOSAT__ = resolve; });
      return jobs.map((job) => ({
        cacheKey: job.cacheKey,
        text: 'MixedRecB\tMixedRecC\t100\t60\t0\t0\t1\t60\t1\t60\t1e-30\t150\n'
      }));
    };
  });
  await page.goto(`${baseUrl}/gbdraw/web/index.html`, { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  const [uidA, uidB, uidC] = await page.evaluate((records) => {
    const app = window.__GBDRAW_APP__;
    app.mode = 'linear';
    app.lInputType = 'gb';
    app.addLinearSeq();
    app.addLinearSeq();
    records.forEach((content, index) => app.setLinearSeqPrimaryFile(index, 'gb', new File(
      [content], `mixed-record-${index + 1}.gbk`, { type: 'text/plain', lastModified: index + 10 }
    )));
    Object.assign(app.form, {
      legend: 'bottom', show_gc: false, show_skew: false,
      show_depth: false, show_labels_linear: 'none'
    });
    app.setLinearLosatProgram('blastn');
    app.losat.executionMode = 'serial';
    const [first, second, third] = app.linearSeqs;
    const blastRow = 'MixedRecA\tMixedRecB\t100\t60\t0\t0\t1\t60\t1\t60\t1e-30\t150\n';
    app.linearComparisonPlan.mode = 'adjacent';
    app.linearComparisonPlan.defaultSource = 'upload';
    app.linearComparisonPlan.edges.splice(0, app.linearComparisonPlan.edges.length,
      {
        id: 'upload-a-b', queryUid: first.uid, subjectUid: second.uid,
        included: false, fileActive: true, losatFilenameActive: false, source: 'upload',
        file: new File([blastRow], 'mixed-a-to-b.tsv'), losatFilename: ''
      },
      {
        id: 'retained-b-c', queryUid: second.uid, subjectUid: third.uid,
        included: false, fileActive: false, losatFilenameActive: false, source: 'upload',
        file: new File([blastRow], 'retained-b-to-c.tsv'), losatFilename: ''
      }
    );
    return [first.uid, second.uid, third.uid];
  }, [
    makeComparisonGenbank('MixedRecA', 'atg'),
    makeComparisonGenbank('MixedRecB', 'gct'),
    makeComparisonGenbank('MixedRecC', 'tta')
  ]);
  const requestPairs = (entries) => entries.map((entry) => [
    entry.kind, entry.queryRecordIndex, entry.subjectRecordIndex
  ]);
  const sparseUpload = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const result = await app.runAnalysis();
    return {
      result,
      comparisons: window.__GBDRAW_DIAGRAM_RUNS__.at(-1).comparisons,
      executorCalls: window.__GBDRAW_LOSAT_EXECUTOR_CALLS__
    };
  });
  expect(sparseUpload.result).toEqual({ status: 'ok' });
  expect(requestPairs(sparseUpload.comparisons)).toEqual([['nucleotideBlast', 0, 1]]);
  expect(sparseUpload.executorCalls).toBe(0);

  const independentReuse = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    app.setLinearComparisonLosatFilename('retained-b-c', 'retained-name.tsv');
    app.deactivateLinearComparisonLosatFilename('retained-b-c');
    app.reuseLinearComparisonFile('retained-b-c');
    const afterFile = app.linearComparisonPlan.edges.find((edge) => edge.id === 'retained-b-c');
    const fileOnly = [afterFile.fileActive, afterFile.losatFilenameActive, afterFile.source];
    app.deactivateLinearComparisonFile('retained-b-c');
    app.reuseLinearComparisonLosatFilename('retained-b-c');
    const afterName = app.linearComparisonPlan.edges.find((edge) => edge.id === 'retained-b-c');
    return [fileOnly, [afterName.fileActive, afterName.losatFilenameActive, afterName.source]];
  });
  expect(independentReuse).toEqual([[true, false, 'upload'], [false, true, 'losat']]);

  await page.waitForFunction(() => window.__GBDRAW_APP__?.pyodideReady === true, null, {
    timeout: 240000
  });
  const selectedResolution = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const { buildPairwiseLosatJobSpecs, resolveLinearComparisonPlan } = await import(
      './js/app/linear-comparisons.js'
    );
    const [first, second, third] = app.linearSeqs;
    const upload = app.linearComparisonPlan.edges.find((edge) => edge.id === 'upload-a-b').file;
    app.linearComparisonPlan.mode = 'selected';
    app.linearComparisonPlan.defaultSource = 'losat';
    app.linearComparisonPlan.edges.splice(0, app.linearComparisonPlan.edges.length,
      {
        id: 'selected-upload-a-b', queryUid: first.uid, subjectUid: second.uid,
        included: true, fileActive: true, losatFilenameActive: false,
        source: 'upload', file: upload, losatFilename: ''
      },
      {
        id: 'selected-losat-b-c', queryUid: second.uid, subjectUid: third.uid,
        included: true, fileActive: false, losatFilenameActive: true,
        source: 'losat', file: null, losatFilename: 'selected-b-to-c.tsv'
      },
      {
        id: 'omitted-a-c', queryUid: first.uid, subjectUid: third.uid,
        included: false, fileActive: false, losatFilenameActive: false,
        source: 'upload', file: null, losatFilename: ''
      }
    );
    const programJobs = ['blastn', 'tblastx', 'blastp'].map((program) => {
      const resolution = resolveLinearComparisonPlan({
        plan: app.linearComparisonPlan,
        sequences: app.linearSeqs,
        losatProgram: program,
        blastpMode: 'pairwise'
      });
      return buildPairwiseLosatJobSpecs({ resolution, program, blastpMode: 'pairwise' })
        .map((job) => [job.program, job.edgeKey, job.queryIndex, job.subjectIndex]);
    });
    return [app.linearComparisonResolution.valid, Object.isFrozen(app.linearComparisonResolution),
      app.linearRecordLayoutEnabled,
      app.linearComparisonResolution.edges.map((edge) => [edge.edgeKey, edge.ordinal, edge.source]),
      programJobs];
  });
  expect(selectedResolution).toEqual([true, true, false, [
    [`${uidA}->${uidB}`, 0, 'upload'], [`${uidB}->${uidC}`, 1, 'losat']
  ], [
    [['blastn', `${uidB}->${uidC}`, 1, 2]],
    [['tblastx', `${uidB}->${uidC}`, 1, 2]],
    [['blastp', `${uidB}->${uidC}`, 1, 2]]
  ]]);

  await page.evaluate(() => { window.__GBDRAW_MIXED_RUN__ = window.__GBDRAW_APP__.runAnalysis(); });
  await page.waitForFunction(() => window.__GBDRAW_LOSAT_EXECUTOR_CALLS__ === 1);
  await page.evaluate(() => {
    window.__GBDRAW_APP__.setLinearComparisonGlobalAction('none');
    window.__GBDRAW_RELEASE_LOSAT__();
  });
  const firstMixed = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const result = await window.__GBDRAW_MIXED_RUN__;
    return {
      result,
      modeAndDrafts: [app.linearComparisonPlan.mode, app.linearComparisonPlan.edges.length],
      comparisons: window.__GBDRAW_DIAGRAM_RUNS__.at(-1).comparisons,
      cacheInfo: app.losatCacheInfo.map((entry) => [
        entry.edgeKey, entry.ordinal, entry.queryIndex, entry.subjectIndex, entry.filename
      ]),
      telemetry: app.lastRunInfo.losatTelemetry,
      jobs: window.__GBDRAW_LOSAT_EXECUTOR_JOBS__
    };
  });
  expect(firstMixed.result).toEqual({ status: 'ok' });
  expect(firstMixed.modeAndDrafts).toEqual(['none', 3]);
  expect(requestPairs(firstMixed.comparisons)).toEqual([
    ['nucleotideBlast', 0, 1], ['nucleotideBlast', 1, 2]
  ]);
  expect(firstMixed.cacheInfo).toEqual([[
    `${uidB}->${uidC}`, 1, 1, 2, 'selected-b-to-c.tsv'
  ]]);
  expect(firstMixed.telemetry).toMatchObject({ cacheHits: 0, cacheMisses: 1, uniqueJobs: 1 });
  expect(firstMixed.jobs[0][0].slice(0, 4)).toEqual([`${uidB}->${uidC}`, 1, 1, 2]);
  expect(firstMixed.jobs[0][0][4]).not.toBe('');

  const reordered = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    app.linearComparisonPlan.mode = 'selected';
    app.setLinearRecordLayoutEnabled(true);
    app.moveLinearSeqUp(2);
    app.moveLinearSeqUp(1);
    return [app.linearSeqs.map((record) => record.uid), app.linearComparisonResolution.edges.map(
      (edge) => [edge.edgeKey, edge.queryIndex, edge.subjectIndex]
    )];
  });
  expect(reordered).toEqual([[uidC, uidA, uidB], [
    [`${uidA}->${uidB}`, 1, 2], [`${uidB}->${uidC}`, 2, 0]
  ]]);
  const cachedRun = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const result = await app.runAnalysis();
    return {
      result,
      comparisons: window.__GBDRAW_DIAGRAM_RUNS__.at(-1).comparisons,
      executorCalls: window.__GBDRAW_LOSAT_EXECUTOR_CALLS__,
      telemetry: app.lastRunInfo.losatTelemetry,
      cacheInfo: app.losatCacheInfo.map((entry) => [
        entry.edgeKey, entry.queryIndex, entry.subjectIndex, entry.filename
      ])
    };
  });
  expect(cachedRun.result).toEqual({ status: 'ok' });
  expect(cachedRun.executorCalls).toBe(1);
  expect(cachedRun.telemetry).toMatchObject({ cacheHits: 1, cacheMisses: 0, uniqueJobs: 0 });
  expect(requestPairs(cachedRun.comparisons)).toEqual([
    ['nucleotideBlast', 1, 2], ['nucleotideBlast', 2, 0]
  ]);
  expect(cachedRun.cacheInfo).toEqual([[
    `${uidB}->${uidC}`, 2, 0, 'selected-b-to-c.tsv'
  ]]);

  await page.evaluate(() => { window.__GBDRAW_APP__.sessionTitle = 'linear-comparison-browser-matrix'; });
  const sessionDownloadPromise = page.waitForEvent('download', { timeout: 120000 });
  expect((await page.evaluate(() => window.__GBDRAW_APP__.saveSessionWithTitle())).status).toBe('saved');
  const sessionPath = await (await sessionDownloadPromise).path();
  await page.reload({ waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  const dialogPromise = page.waitForEvent('dialog', { timeout: 120000 });
  await page.locator('input[accept^=".json,"]').first().setInputFiles(sessionPath);
  const dialog = await dialogPromise;
  expect(dialog.message()).toBe('Session loaded successfully!');
  await dialog.accept();
  await page.waitForFunction(() => window.__GBDRAW_APP__?.losatCacheInfo?.length === 1);
  const restored = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const { state } = await import('./js/state.js');
    return {
      mode: app.linearComparisonPlan.mode,
      resolution: app.linearComparisonResolution.edges.map(
        (edge) => [edge.edgeKey, edge.queryIndex, edge.subjectIndex]
      ),
      cache: app.losatCacheInfo.map((entry) => [
        entry.edgeKey, entry.ordinal, entry.queryUid, entry.subjectUid,
        entry.queryIndex, entry.subjectIndex, entry.filename
      ]),
      cacheSize: state.losatCache.value.size
    };
  });
  expect(restored).toEqual({
    mode: 'selected',
    resolution: [[`${uidA}->${uidB}`, 1, 2], [`${uidB}->${uidC}`, 2, 0]],
    cache: [[`${uidB}->${uidC}`, 1, uidB, uidC, 2, 0, 'selected-b-to-c.tsv']],
    cacheSize: 1
  });
  const rawDownloadPromise = page.waitForEvent('download', { timeout: 120000 });
  await page.evaluate((edgeKey) => window.__GBDRAW_APP__.downloadLosatPair(edgeKey, ''), `${uidB}->${uidC}`);
  const rawDownload = await rawDownloadPromise;
  expect(rawDownload.suggestedFilename()).toBe('selected-b-to-c.tsv');
  expect(readFileSync(await rawDownload.path(), 'utf8')).toBe(
    'MixedRecB\tMixedRecC\t100\t60\t0\t0\t1\t60\t1\t60\t1e-30\t150\n'
  );
});

test('Candidate render post-processing sanitizes and reapplies stable styles before commit', async ({ page }) => {
  await page.goto(`${baseUrl}/gbdraw/web/index.html`, { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);

  const outcome = await page.evaluate(async () => {
    const { prepareCandidateRenderCommit } = await import('./js/app/candidate-render.js');
    const stableKey = `record-1\u0000feature-1`;
    const sourceResult = {
      name: 'candidate.svg',
      format: 'svg',
      content: [
        '<svg xmlns="http://www.w3.org/2000/svg">',
        '<script>throw new Error("unsafe")</script>',
        '<path id="rendered-1" data-gbdraw-feature-id="shared-feature"',
        ' data-gbdraw-rendered-feature-id="rendered-1"',
        ' data-gbdraw-feature-part="block" fill="#000000" stroke="#000000" stroke-width="1"/>',
        '<path id="rendered-2" data-gbdraw-feature-id="shared-feature"',
        ' data-gbdraw-rendered-feature-id="rendered-2"',
        ' data-gbdraw-feature-part="block" fill="#000000" stroke="#000000" stroke-width="1"/>',
        '</svg>'
      ].join('')
    };
    const catalog = {
      schema: 3,
      items: [{
        resultIndex: 0,
        resultName: 'candidate.svg',
        recordKeys: ['record-1'],
        features: [{
          svgId: 'rendered-1',
          recordKey: 'record-1',
          biologicalFeatureId: 'feature-1'
        }, {
          svgId: 'rendered-2',
          recordKey: 'record-1',
          biologicalFeatureId: 'feature-2'
        }],
        biologicalFeatures: [{
          recordKey: 'record-1',
          biologicalFeatureId: 'feature-1',
          stableFeatureId: 'stable-feature-1',
          type: 'CDS',
          record_idx: 0,
          record_id: 'record-1',
          start: 0,
          end: 10,
          strand: 1,
          qualifiers: {}
        }, {
          recordKey: 'record-1',
          biologicalFeatureId: 'feature-2',
          stableFeatureId: 'stable-feature-2',
          type: 'CDS',
          record_idx: 0,
          record_id: 'record-1',
          start: 20,
          end: 30,
          strand: 1,
          qualifiers: {}
        }],
        orthogroups: [],
        annotations: [],
        comparisonMatches: []
      }]
    };
    const prepared = prepareCandidateRenderCommit({
      results: [sourceResult],
      catalog,
      featureColorOverrides: {
        [stableKey]: { color: '#ff00ff', caption: 'Candidate style' }
      },
      featureStrokeOverrides: {
        [stableKey]: { strokeColor: '#ff00ff', strokeWidth: 5 }
      }
    });
    const svg = new DOMParser().parseFromString(
      prepared.results[0].content,
      'image/svg+xml'
    );
    const feature = svg.getElementById('rendered-1');
    const sibling = svg.getElementById('rendered-2');
    return {
      sourceUnchanged:
        sourceResult.content.includes('<script>')
        && sourceResult.content.includes('fill="#000000"'),
      scriptRemoved: !prepared.results[0].content.includes('<script>'),
      fill: feature?.getAttribute('fill'),
      stroke: feature?.getAttribute('stroke'),
      strokeWidth: feature?.getAttribute('stroke-width'),
      siblingFill: sibling?.getAttribute('fill'),
      renderedBindingPreserved:
        feature?.getAttribute('data-gbdraw-rendered-feature-id') === 'rendered-1',
      featureCount: prepared.featureState.extractedFeatures.length
    };
  });

  expect(outcome).toEqual({
    sourceUnchanged: true,
    scriptRemoved: true,
    fill: '#ff00ff',
    stroke: '#ff00ff',
    strokeWidth: '5',
    siblingFill: '#000000',
    renderedBindingPreserved: true,
    featureCount: 2
  });
});

test('Label-scoped feature colors and legends survive linear regeneration', async ({ page }) => {
  test.setTimeout(240000);
  const makeGenbank = (recordId) => {
    const sequence = 'atg'.repeat(100);
    const origin = sequence.match(/.{1,60}/g).map((chunk, index) => {
      const groups = chunk.match(/.{1,10}/g).join(' ');
      return `${String(index * 60 + 1).padStart(9)} ${groups}`;
    }).join('\n');
    return `LOCUS       ${recordId.padEnd(24)} 300 bp    DNA     linear   UNA 01-JAN-2000
DEFINITION  feature color regeneration test.
ACCESSION   ${recordId}
VERSION     ${recordId}
KEYWORDS    .
SOURCE      synthetic construct
  ORGANISM  synthetic construct
            .
FEATURES             Location/Qualifiers
     CDS             1..90
                     /product="wsv360-like protein"
ORIGIN
${origin}
//
`;
  };

  await page.goto(`${baseUrl}/gbdraw/web/index.html`, { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);

  await page.evaluate(({ firstRecord, secondRecord }) => {
    const app = window.__GBDRAW_APP__;
    app.mode = 'linear';
    app.lInputType = 'gb';
    app.addLinearSeq();
    app.setLinearSeqPrimaryFile(0, 'gb', new File([firstRecord], 'record-a.gbk', {
      type: 'text/plain', lastModified: 1
    }));
    app.setLinearSeqPrimaryFile(1, 'gb', new File([secondRecord], 'record-b.gbk', {
      type: 'text/plain', lastModified: 2
    }));
    Object.assign(app.form, {
      legend: 'bottom',
      show_gc: false,
      show_skew: false,
      show_depth: false,
      show_labels_linear: 'none'
    });
  }, {
    firstRecord: makeGenbank('ColorRecA'),
    secondRecord: makeGenbank('ColorRecB')
  });
  await page.waitForFunction(
    () => window.__GBDRAW_APP__?.pyodideReady === true,
    null,
    { timeout: 180000 }
  );

  const firstRun = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const result = await app.runAnalysis();
    return { result, errorLog: app.errorLog };
  });
  expect(firstRun).toEqual({ result: { status: 'ok' }, errorLog: null });

  const edited = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const target = app.extractedFeatures.find((feature) => feature.product === 'wsv360-like protein');
    if (!target) throw new Error('Expected product feature was not extracted.');
    await app.requestFeatureColorChange(target, '#ff00ff');
    const dialog = {
      displayLabel: app.colorScopeDialog.displayLabel,
      displayLabelSiblingCount: app.colorScopeDialog.displayLabelSiblingCount
    };
    await app.handleColorScopeChoice('displayLabel');
    app.extractedFeatures
      .filter((feature) => feature.product === 'wsv360-like protein')
      .forEach((feature) => {
        const key = String(feature.stable_override_key || feature.id || '');
        if (!key) throw new Error('Expected a stable feature override key.');
        app.featureStrokeOverrides[key] = {
          strokeColor: '#ff00ff',
          strokeWidth: 5
        };
      });
    return {
      dialog,
      rules: app.manualSpecificRules.map(({ feat, qual, val, color, cap }) => ({
        feat, qual, val, color, cap
      }))
    };
  });
  expect(edited.dialog).toEqual({
    displayLabel: 'wsv360-like protein',
    displayLabelSiblingCount: 1
  });
  expect(edited.rules).toEqual([{
    feat: 'CDS',
    qual: 'product',
    val: '^wsv360-like protein$',
    color: '#ff00ff',
    cap: 'wsv360-like protein'
  }]);

  const secondRun = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    app.editableLabels = [{
      key: 'stale-label',
      idx: 1,
      text: 'stale',
      sourceText: 'stale',
      featureId: 'stale-feature',
      kind: 'regular',
      draftText: 'stale'
    }];
    const result = await app.runAnalysis();
    await window.Vue.nextTick();
    await window.Vue.nextTick();
    return {
      result,
      errorLog: app.errorLog,
      hasStaleEditableLabel: app.editableLabels.some(
        (entry) => entry.key === 'stale-label'
      )
    };
  });
  expect(secondRun).toEqual({
    result: { status: 'ok' },
    errorLog: null,
    hasStaleEditableLabel: false
  });

  const regenerated = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const svgText = app.results?.[0]?.content || '';
    const svg = new DOMParser().parseFromString(svgText, 'image/svg+xml').documentElement;
    const featureIds = app.extractedFeatures
      .filter((feature) => feature.product === 'wsv360-like protein')
      .map((feature) => feature.svg_id);
    const featureStyles = featureIds.map((featureId) => {
      const roots = [...svg.querySelectorAll('[data-gbdraw-feature-id]')]
        .filter((element) => element.getAttribute('data-gbdraw-feature-id') === featureId);
      const elements = [...roots, ...roots.flatMap((root) => [...root.querySelectorAll('*')])];
      return {
        fills: elements
          .map((element) => String(element.getAttribute('fill') || '').toLowerCase())
          .filter((fill) => fill && fill !== 'none'),
        strokeColors: elements
          .map((element) => String(element.getAttribute('stroke') || '').toLowerCase())
          .filter(Boolean),
        strokeWidths: elements
          .map((element) => String(element.getAttribute('stroke-width') || ''))
          .filter(Boolean)
      };
    });
    const legendEntries = [...svg.querySelectorAll('[data-legend-key]')]
      .filter((entry) => entry.getAttribute('data-legend-key') === 'wsv360-like protein');
    const legendFills = legendEntries.flatMap((entry) => [...entry.querySelectorAll('[fill]')])
      .map((element) => String(element.getAttribute('fill') || '').toLowerCase())
      .filter((fill) => fill && fill !== 'none');
    return {
      featureIds,
      featureStyles,
      legendEntryCount: legendEntries.length,
      legendFills,
      rules: app.manualSpecificRules.map(({ feat, qual, val, color, cap }) => ({
        feat, qual, val, color, cap
      }))
    };
  });

  expect(regenerated.featureIds).toHaveLength(2);
  regenerated.featureStyles.forEach(({ fills, strokeColors, strokeWidths }) => {
    expect(fills).toContain('#ff00ff');
    expect(strokeColors).toContain('#ff00ff');
    expect(strokeWidths).toContain('5');
  });
  expect(regenerated.legendEntryCount).toBeGreaterThan(0);
  expect(regenerated.legendFills).toContain('#ff00ff');
  expect(regenerated.rules).toEqual(edited.rules);

  const postprocessingFailure = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const targetId = String(app.extractedFeatures?.[0]?.svg_id || '');
    app.selectedFeatureIds = new Set(targetId ? [targetId] : []);
    app.selectedFeatureAnchorId = targetId;
    app.editableLabels = [{
      key: 'postprocess-rollback-label',
      idx: 1,
      text: 'Keep after post-processing failure',
      sourceText: 'Original rollback label',
      featureId: targetId,
      kind: 'regular',
      draftText: 'Keep after post-processing failure'
    }];
    const snapshot = () => JSON.stringify({
      results: app.results,
      resultGenerationKey: app.resultGenerationKey,
      selectedResultIndex: app.selectedResultIndex,
      selectedFeatureIds: [...app.selectedFeatureIds].sort(),
      selectedFeatureAnchorId: app.selectedFeatureAnchorId,
      featureColorOverrides: app.featureColorOverrides,
      featureStrokeOverrides: app.featureStrokeOverrides,
      extractedFeatures: app.extractedFeatures,
      featureCatalog: app.featureCatalog,
      legendEntries: app.legendEntries,
      lastRunInfo: app.lastRunInfo,
      appliedPaletteName: app.appliedPaletteName,
      appliedPaletteColors: app.appliedPaletteColors,
      pendingPaletteName: app.pendingPaletteName,
      pendingPaletteColors: app.pendingPaletteColors,
      editableLabels: app.editableLabels
    });
    const before = snapshot();
    const originalSanitize = window.DOMPurify.sanitize;
    let result;
    try {
      window.DOMPurify.sanitize = () => {
        throw new Error('Forced candidate post-processing failure.');
      };
      result = await app.runAnalysis();
    } finally {
      window.DOMPurify.sanitize = originalSanitize;
    }
    await window.Vue.nextTick();
    const beforeState = JSON.parse(before);
    const afterState = JSON.parse(snapshot());
    return {
      result,
      errorSummary: String(app.errorLog?.summary || ''),
      snapshotPreserved: JSON.stringify(afterState) === before,
      changedFields: Object.keys(beforeState).filter(
        (key) => JSON.stringify(beforeState[key]) !== JSON.stringify(afterState[key])
      )
    };
  });
  expect(postprocessingFailure).toEqual({
    result: { status: 'error' },
    errorSummary: 'Forced candidate post-processing failure.',
    snapshotPreserved: true,
    changedFields: []
  });

  const staleResponse = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    app.errorLog = null;
    const snapshot = () => JSON.stringify({
      results: app.results,
      resultGenerationKey: app.resultGenerationKey,
      selectedResultIndex: app.selectedResultIndex,
      selectedFeatureIds: [...app.selectedFeatureIds].sort(),
      selectedFeatureAnchorId: app.selectedFeatureAnchorId,
      featureColorOverrides: app.featureColorOverrides,
      featureStrokeOverrides: app.featureStrokeOverrides,
      extractedFeatures: app.extractedFeatures,
      featureCatalog: app.featureCatalog,
      legendEntries: app.legendEntries,
      lastRunInfo: app.lastRunInfo,
      appliedPaletteName: app.appliedPaletteName,
      appliedPaletteColors: app.appliedPaletteColors,
      pendingPaletteName: app.pendingPaletteName,
      pendingPaletteColors: app.pendingPaletteColors,
      editableLabels: app.editableLabels
    });
    const before = snapshot();
    const originalPostMessage = Worker.prototype.postMessage;
    let delayedRunObserved = false;
    Worker.prototype.postMessage = function delayedFirstRun(message, transfer) {
      if (!delayedRunObserved && message?.type === 'run') {
        delayedRunObserved = true;
        window.setTimeout(() => {
          originalPostMessage.call(this, message, transfer);
        }, 1500);
        return;
      }
      originalPostMessage.call(this, message, transfer);
    };

    let firstResult;
    let secondResult;
    try {
      const firstRun = app.runAnalysis();
      for (let attempt = 0; attempt < 500 && !delayedRunObserved; attempt += 1) {
        await new Promise((resolve) => window.setTimeout(resolve, 2));
      }
      if (!delayedRunObserved) throw new Error('The first worker request was not observed.');
      const secondRun = app.runAnalysis();
      [firstResult, secondResult] = await Promise.all([firstRun, secondRun]);
    } finally {
      Worker.prototype.postMessage = originalPostMessage;
    }
    await window.Vue.nextTick();
    const beforeState = JSON.parse(before);
    const afterState = JSON.parse(snapshot());
    return {
      firstResult,
      secondResult,
      errorLog: app.errorLog,
      snapshotPreserved: JSON.stringify(afterState) === before,
      changedFields: Object.keys(beforeState).filter(
        (key) => JSON.stringify(beforeState[key]) !== JSON.stringify(afterState[key])
      )
    };
  });
  expect(staleResponse).toEqual({
    firstResult: { status: 'stale' },
    secondResult: { status: 'error' },
    errorLog: null,
    snapshotPreserved: true,
    changedFields: []
  });

  const reset = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const target = app.extractedFeatures.find((feature) => feature.product === 'wsv360-like protein');
    if (!target) throw new Error('Expected a feature to reset.');
    const defaultColor = app.appliedPaletteColors.CDS;
    app.clickedFeature = {
      svg_id: target.svg_id,
      feat: target,
      color: '#ff00ff'
    };
    app.resetColorDialog.defaultColor = defaultColor;
    app.resetColorDialog.caption = 'wsv360-like protein';
    await app.handleResetColorChoice('this');
    return {
      defaultColor,
      targetId: target.svg_id,
      rules: app.manualSpecificRules.map(({ feat, qual, val, color, cap }) => ({
        feat, qual, val, color, cap
      }))
    };
  });
  expect(reset.rules).toContainEqual({
    feat: 'CDS',
    qual: 'hash',
    val: reset.targetId.replace(/_record_\d+$/i, ''),
    color: reset.defaultColor,
    cap: 'other proteins'
  });
  expect(reset.rules).toContainEqual(edited.rules[0]);

  const thirdRun = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const result = await app.runAnalysis();
    const svg = new DOMParser().parseFromString(app.results?.[0]?.content || '', 'image/svg+xml').documentElement;
    const fillsById = Object.fromEntries(
      app.extractedFeatures
        .filter((feature) => feature.product === 'wsv360-like protein')
        .map((feature) => {
          const roots = [...svg.querySelectorAll('[data-gbdraw-feature-id]')]
            .filter((element) => element.getAttribute('data-gbdraw-feature-id') === feature.svg_id);
          const fills = [...roots, ...roots.flatMap((root) => [...root.querySelectorAll('[fill]')])]
            .map((element) => String(element.getAttribute('fill') || '').toLowerCase())
            .filter((fill) => fill && fill !== 'none');
          return [feature.svg_id, fills];
        })
    );
    return { result, errorLog: app.errorLog, fillsById };
  });
  expect(thirdRun.result).toEqual({ status: 'ok' });
  expect(thirdRun.errorLog).toBeNull();
  expect(thirdRun.fillsById[reset.targetId]).toContain(reset.defaultColor.toLowerCase());
  const siblingFills = Object.entries(thirdRun.fillsById)
    .filter(([featureId]) => featureId !== reset.targetId)
    .flatMap(([, fills]) => fills);
  expect(siblingFills).toContain('#ff00ff');
});

test('Region annotations expose and persist an explicit target-record selection', async ({ page }) => {
  await page.goto(`${baseUrl}/gbdraw/web/index.html`, { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);

  const genbank = `LOCUS       RecA                      10 bp    DNA     linear   UNA 01-JAN-2000
DEFINITION  first.
ACCESSION   RecA
VERSION     RecA
KEYWORDS    .
SOURCE      synthetic construct
  ORGANISM  synthetic construct
            .
FEATURES             Location/Qualifiers
ORIGIN
        1 aaaaaaaaaa
//
LOCUS       RecB                      12 bp    DNA     linear   UNA 01-JAN-2000
DEFINITION  second.
ACCESSION   RecB
VERSION     RecB
KEYWORDS    .
SOURCE      synthetic construct
  ORGANISM  synthetic construct
            .
FEATURES             Location/Qualifiers
ORIGIN
        1 cccccccccccc
//
`;
  await page.evaluate((content) => {
    const app = window.__GBDRAW_APP__;
    app.mode = 'linear';
    app.lInputType = 'gb';
    app.setLinearSeqPrimaryFile(0, 'gb', new File([content], 'two-records.gb', {
      type: 'text/plain',
      lastModified: 1
    }));
    const set = app.addAnnotationSet('review');
    app.addCoordinateAnnotation(set, { start: 1, end: 10 });
  }, genbank);

  await page.getByText('Region Annotations', { exact: false }).click();
  const selector = page.getByLabel('Annotation target record');
  await expect(selector).toHaveCount(1);
  await expect(selector).toHaveValue('');
  await expect(selector.locator('option')).toHaveText([
    'Select target record',
    '#1 · RecA · 10 bp',
    '#2 · RecB · 12 bp'
  ], { timeout: 60000 });
  await expect(page.getByText('Choose the record that this annotation targets.')).toBeVisible();

  const rejected = await page.evaluate(() => window.__GBDRAW_APP__.runAnalysis());
  expect(rejected).toEqual({ status: 'error' });
  await expect(page.getByText('Choose a target record for region annotation review/region_1.')).toBeVisible();

  await selector.selectOption({ label: '#2 · RecB · 12 bp' });
  await expect(page.getByText('Choose the record that this annotation targets.')).toHaveCount(0);
  const target = await page.evaluate(() => (
    window.__GBDRAW_APP__.annotationSets[0].annotations[0].target.record
  ));
  expect(target).toEqual({ kind: 'recordId', value: 'RecB' });
});

test('Region annotation IDs accept continuous typing without losing focus', async ({ page }) => {
  await page.goto(`${baseUrl}/gbdraw/web/index.html`, { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);

  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const set = app.addAnnotationSet('review');
    app.addCoordinateAnnotation(set);
  });

  await page.getByText('Region Annotations', { exact: false }).click();
  const idInput = page.getByPlaceholder('annotation_id');
  await idInput.selectText();
  await idInput.pressSequentially('Repeat');

  await expect(idInput).toBeFocused();
  await expect(idInput).toHaveValue('Repeat');
  await idInput.press('Tab');
  await expect.poll(() => page.evaluate(() => (
    window.__GBDRAW_APP__.annotationSets[0].annotations[0].id
  ))).toBe('Repeat');
});

test('GFF annotation targets follow FASTA record order', async ({ page }) => {
  await page.goto(`${baseUrl}/gbdraw/web/index.html`, { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);

  const gff = `##gff-version 3
##sequence-region RecB 1 12
RecB\ttest\tgene\t1\t3\t.\t+\t.\tID=gene_b
##sequence-region RecA 1 10
RecA\ttest\tgene\t2\t4\t.\t+\t.\tID=gene_a
`;
  const fasta = `>RecB
CCCCCCCCCCCC
>RecA
AAAAAAAAAA
`;
  await page.evaluate(({ gffText, fastaText }) => {
    const app = window.__GBDRAW_APP__;
    app.mode = 'linear';
    app.lInputType = 'gff';
    app.setLinearSeqPrimaryFile(0, 'gff', new File([gffText], 'records.gff3', {
      type: 'text/plain',
      lastModified: 2
    }));
    app.setLinearSeqPrimaryFile(0, 'fasta', new File([fastaText], 'records.fasta', {
      type: 'text/plain',
      lastModified: 3
    }));
    const set = app.addAnnotationSet('gff-review');
    app.addCoordinateAnnotation(set, { start: 1, end: 3 });
  }, { gffText: gff, fastaText: fasta });

  await page.getByText('Region Annotations', { exact: false }).click();
  await expect(page.getByLabel('Annotation target record').locator('option')).toHaveText([
    'Select target record',
    '#1 · RecB · 12 bp',
    '#2 · RecA · 10 bp'
  ], { timeout: 60000 });
});
