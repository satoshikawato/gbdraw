// The diagram Worker's Pyodide snapshot boot (docs/internal/SELECTIVE_CI.md,
// "Pyodide snapshot"). With GBDRAW_TEST_PYODIDE_SNAPSHOT_DIR set, the test
// server serves a snapshot and the Worker restores it; on any restore error the
// page discards that Worker and boots one cold. A snapshot that stops
// restoring, after a Pyodide upgrade for example, fails this test instead of
// silently making every other spec boot cold.
const { test, expect } = require('@playwright/test');
const { join } = require('node:path');
const { openApp } = require('./helpers/app-lifecycle.cjs');
const { generate } = require('./helpers/mode-transition.cjs');
const { forceColdPyodideBoot } = require('./helpers/pyodide-snapshot-test-server.cjs');

const SESSION = join(process.cwd(), 'gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json');

// Loads a saved Circular Session (no Worker yet), Generates, and returns the
// Worker boots the page reported and the generated SVG.
const bootAndGenerate = async (browser, { cold }) => {
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
    if (cold) await forceColdPyodideBoot(page);
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
  test.setTimeout(360_000);
  const restored = await bootAndGenerate(browser, { cold: false });
  expect(restored.boots).toEqual(['snapshot']);
  expect(restored.svg).toContain('<svg');

  const fallback = await bootAndGenerate(browser, { cold: true });
  expect(fallback.boots).toHaveLength(1);
  expect(fallback.boots[0]).toMatch(/^cold after snapshot error: /);
  expect(fallback.svg).toBe(restored.svg);
});
