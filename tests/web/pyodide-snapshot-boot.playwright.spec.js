// The diagram Worker's Pyodide snapshot boot (docs/internal/SELECTIVE_CI.md,
// "Pyodide snapshot"). With GBDRAW_TEST_PYODIDE_SNAPSHOT_DIR set, the test
// server serves a snapshot and the Worker restores it; on any restore error the
// page discards that Worker and boots one cold: when the snapshot files fail
// to download, when they were built from other inputs, and when a side module
// is not linked after the restore. A snapshot that stops restoring, after a
// Pyodide upgrade for example, fails this test instead of silently making
// every other spec boot cold.
const { test, expect } = require('@playwright/test');
const { join } = require('node:path');
const { openApp } = require('./helpers/app-lifecycle.cjs');
const { generate } = require('./helpers/mode-transition.cjs');
const { PYODIDE_SNAPSHOT_TEST_PATH, forceColdPyodideBoot } = require('./helpers/pyodide-snapshot-test-server.cjs');

const SESSION = join(process.cwd(), 'gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json');

// Serves the snapshot's state.json changed by `change`.
const changeSnapshotState = (page, change) => page.route(
  (url) => url.pathname === `${PYODIDE_SNAPSHOT_TEST_PATH}state.json`,
  async (route) => {
    const response = await route.fetch();
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
