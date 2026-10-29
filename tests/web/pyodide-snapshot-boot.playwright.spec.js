// The diagram Worker's boot. A chosen file makes the page download the Worker's
// boot files (services/pyodide-assets.js), and the first Generate's Worker
// reads them from the HTTP cache.
//
// The Pyodide snapshot boot (docs/internal/SELECTIVE_CI.md,
// "Pyodide snapshot"): with GBDRAW_TEST_PYODIDE_SNAPSHOT_DIR set, the test
// server serves a snapshot and the Worker restores it; on any restore error the
// page discards that Worker and boots one cold: when the snapshot files fail
// to download, when they were built from other inputs, and when a side module
// is not linked after the restore. A snapshot that stops restoring, after a
// Pyodide upgrade for example, fails this test instead of silently making
// every other spec boot cold.
const { test, expect } = require('@playwright/test');
const { spawn } = require('node:child_process');
const { rmSync } = require('node:fs');
const { join } = require('node:path');
const { fetchRoute, openApp } = require('./helpers/app-lifecycle.cjs');
const { generate } = require('./helpers/mode-transition.cjs');
const { PYODIDE_SNAPSHOT_TEST_PATH, forceColdPyodideBoot } = require('./helpers/pyodide-snapshot-test-server.cjs');

const SESSION = join(process.cwd(), 'gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json');
const GENBANK = join(process.cwd(), 'tests/fixtures/sessions/cli-web-mito.gb');

// Serves the snapshot's state.json changed by `change`.
const changeSnapshotState = (page, change) => page.route(
  (url) => url.pathname === `${PYODIDE_SNAPSHOT_TEST_PATH}state.json`,
  async (route) => {
    const response = await fetchRoute(route);
    const state = await response.json();
    change(state);
    await route.fulfill({ response, json: state });
  }
);

// Drops numpy's side modules from the restore: the snapshot memory still
// holds numpy, imported, with its function table slots empty.
const unlinkNumpy = (state) => {
  const numpy = new Set(state.libs.filter(({ name }) => name.includes('/numpy/')).map(({ name }) => name));
  if (!numpy.size) throw new Error('The snapshot links no numpy side module.');
  state.libs = state.libs.filter(({ name }) => !numpy.has(name));
  state.extraSlots = state.extraSlots.filter(([, library]) => !numpy.has(library));
};

// Loads a saved Circular Session (no Worker yet), Generates, and returns the
// Worker boots the page reported and the generated SVG. `prepare` runs before
// the app opens.
const bootAndGenerate = async (browser, prepare = async () => {}) => {
  const context = await browser.newContext({ viewport: { width: 1600, height: 1000 } });
  const page = await context.newPage();
  try {
    await page.addInitScript(() => {
      window.__WORKER_BOOTS__ = [];
      window.__GBDRAW_TEST_HOOKS__ = {
        onStructuralMetric: (metric) => {
          if (metric.name === 'workerBoot') window.__WORKER_BOOTS__.push(metric.boot);
        }
      };
    });
    await prepare(page);
    page.on('dialog', (dialog) => dialog.dismiss());
    await openApp(page);
    await page.locator('input[accept^=".json,"]').setInputFiles(SESSION);
    await expect.poll(() => page.evaluate(() => !window.__GBDRAW_APP__.sessionImportPending
      && window.__GBDRAW_APP__.extractedFeatures.length > 0), { timeout: 180_000 }).toBe(true);
    await generate(page);
    return await page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      return { boots: window.__WORKER_BOOTS__, svg: state.svgContainer.value?.querySelector('svg')?.outerHTML || '' };
    });
  } finally {
    await context.close();
  }
};

test('the diagram Worker restores the served Pyodide snapshot and boots cold when it cannot', async ({ browser }) => {
  test.skip(!process.env.GBDRAW_TEST_PYODIDE_SNAPSHOT_DIR, 'The test server serves no Pyodide snapshot.');
  test.setTimeout(600_000);
  const restored = await bootAndGenerate(browser);
  expect(restored.boots).toEqual(['snapshot']);
  expect(restored.svg).toContain('<svg');

  for (const [name, prepare, reason] of [
    ['snapshot files that fail to download', forceColdPyodideBoot, /^cold after snapshot error: /],
    [
      'a snapshot built from other inputs',
      (page) => changeSnapshotState(page, (state) => { state.key = '0'.repeat(64); }),
      /^cold after snapshot error: Pyodide snapshot was built from another recipe/
    ],
    ['a side module not linked after the restore', (page) => changeSnapshotState(page, unlinkNumpy), /^cold after snapshot error: /]
  ]) {
    await test.step(name, async () => {
      const fallback = await bootAndGenerate(browser, prepare);
      expect(fallback.boots).toHaveLength(1);
      expect(fallback.boots[0]).toMatch(reason);
      expect(fallback.svg).toBe(restored.svg);
    });
  }
});

// The repository on a plain `python3 -m http.server` at a free port, every
// file fresh for an hour so that no response depends on a file's modification
// time: its origin, and a function that stops it. Rejects when the server
// exits or stays silent for 10 s.
const PLAIN_SERVER = [
  'import http.server',
  'class Handler(http.server.SimpleHTTPRequestHandler):',
  '    def end_headers(self):',
  '        self.send_header("Cache-Control", "max-age=3600")',
  '        super().end_headers()',
  'http.server.test(HandlerClass=Handler, port=0, bind="127.0.0.1")'
].join('\n');

const startPlainServer = () => new Promise((resolveServer, rejectServer) => {
  const server = spawn('python3', ['-u', '-c', PLAIN_SERVER], { cwd: process.cwd(), stdio: ['ignore', 'pipe', 'ignore'] });
  const fail = (message) => { server.kill(); rejectServer(new Error(message)); };
  const timer = setTimeout(() => fail('The test server printed no port.'), 10_000);
  server.once('error', (error) => { clearTimeout(timer); rejectServer(error); });
  server.once('exit', (code) => { clearTimeout(timer); rejectServer(new Error(`The test server exited (${code}) before printing its port.`)); });
  server.stdout.on('data', (chunk) => {
    const port = /port (\d+)/.exec(String(chunk))?.[1];
    if (!port) return;
    clearTimeout(timer);
    resolveServer({ origin: `http://127.0.0.1:${port}`, stop: () => server.kill() });
  });
});

// As a page on gbdraw.app loads: a browser profile with its disk cache (the
// off-the-record cache of a test context keeps no entry as large as
// pyodide.asm.wasm), no routing (a routed request bypasses the HTTP cache) and
// no snapshot.
test('the first Worker boots from the files the page prefetched when a file was chosen', async ({ playwright }, testInfo) => {
  test.setTimeout(300_000);
  const server = await startPlainServer();
  const profile = testInfo.outputPath('profile');
  let context;
  try {
    context = await playwright.chromium.launchPersistentContext(profile, {
      baseURL: server.origin, viewport: { width: 1280, height: 900 }
    });
    const page = context.pages()[0] || await context.newPage();
    await page.addInitScript(() => performance.setResourceTimingBufferSize(10_000));
    await openApp(page);
    const bootUrls = await page.evaluate(async () => (
      await import('./js/services/pyodide-assets.js')).resolvePyodideBootAssetUrls());
    await page.getByLabel('GenBank/DDBJ File', { exact: true }).setInputFiles(GENBANK);
    // The prefetch has finished once the page has a resource entry (written at
    // the end of the body) for every boot file.
    await expect.poll(() => page.evaluate((urls) => {
      const loaded = new Set(performance.getEntriesByType('resource').map(({ name }) => name));
      return urls.filter((url) => !loaded.has(url));
    }, bootUrls), { timeout: 180_000 }).toEqual([]);
    await generate(page);
    const worker = page.workers().find((candidate) => candidate.url().includes('/diagram-generation-worker.js'));
    expect(worker).toBeDefined();
    const requests = await /** @type {import('@playwright/test').Worker} */ (worker).evaluate(() => (
      performance.getEntriesByType('resource').map((entry) => {
        const { name, transferSize, encodedBodySize } = /** @type {PerformanceResourceTiming} */ (entry);
        return { name, transferSize, encodedBodySize };
      })));
    // Everything the Worker requested besides its own modules is a boot file,
    // each requested once.
    const boot = requests.filter(({ name }) => !new URL(name).pathname.startsWith('/gbdraw/web/js/'));
    expect(boot.map(({ name }) => name).sort()).toEqual([...bootUrls].sort());
    // From the HTTP cache: no body crossed the network (a cache hit transfers
    // nothing, a revalidation only headers).
    expect(boot.filter(({ transferSize, encodedBodySize }) => transferSize !== 0 && transferSize >= encodedBodySize)).toEqual([]);
  } finally {
    try {
      await context?.close();
    } finally {
      server.stop();
      rmSync(profile, { recursive: true, force: true });
    }
  }
});
