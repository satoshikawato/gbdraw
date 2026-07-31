const { test, expect } = require('@playwright/test');
const { createReadStream, existsSync, readFileSync } = require('node:fs');
const { createServer } = require('node:http');
const { extname, join, normalize, resolve, sep } = require('node:path');

const repoRoot = resolve(process.env.GBDRAW_REPO || process.cwd());
const genbankPath = join(repoRoot, 'tests/test_inputs/HmmtDNA.gbk');
const contentTypes = {
  '.html': 'text/html; charset=utf-8',
  '.js': 'text/javascript; charset=utf-8',
  '.mjs': 'text/javascript; charset=utf-8',
  '.css': 'text/css; charset=utf-8',
  '.json': 'application/json; charset=utf-8',
  '.svg': 'image/svg+xml',
  '.wasm': 'application/wasm',
  '.whl': 'application/octet-stream',
  '.data': 'application/octet-stream'
};

let server;
let baseUrl;

test.beforeAll(async () => {
  await new Promise((resolveServer, rejectServer) => {
    server = createServer((request, response) => {
      const url = new URL(request.url || '/', 'http://127.0.0.1');
      const requestedPath = normalize(decodeURIComponent(url.pathname))
        .replace(/^(\.\.(?:\/|\\|$))+/, '');
      const filePath = resolve(repoRoot, requestedPath.replace(/^[/\\]+/, ''));
      if (!filePath.startsWith(`${repoRoot}${sep}`) && filePath !== repoRoot) {
        response.writeHead(403);
        response.end('Forbidden');
        return;
      }
      const finalPath = filePath === repoRoot
        ? join(repoRoot, 'gbdraw/web/index.html')
        : filePath;
      if (!existsSync(finalPath)) {
        response.writeHead(404);
        response.end('Not found');
        return;
      }
      response.writeHead(200, {
        'Content-Type': contentTypes[extname(finalPath)] || 'application/octet-stream'
      });
      createReadStream(finalPath).pipe(response);
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

const runDiagram = async (page) => page.evaluate(async () => {
  const app = window.__GBDRAW_APP__;
  const result = await app.runAnalysis();
  return {
    result,
    errorSummary: String(app.errorLog?.summary || ''),
    errorDetails: Array.isArray(app.errorLog?.details)
      ? app.errorLog.details.map((detail) => String(detail))
      : []
  };
});

const inspectResultSvg = async (page) => page.evaluate(() => {
  const content = window.__GBDRAW_APP__.results?.[0]?.content || '';
  const svg = new DOMParser().parseFromString(content, 'image/svg+xml').documentElement;
  const recordGroups = [...svg.querySelectorAll('g[data-gbdraw-record-index]')]
    .filter((group) => String(group.id || '').startsWith('record_group_'));
  return {
    axis: Boolean(svg.getElementById('Axis')),
    ticks: Boolean(
      svg.getElementById('tick') ||
      svg.getElementById('ticks') ||
      svg.querySelector('[data-gbdraw-slot-renderer="ticks"]')
    ),
    lengthBar: Boolean(svg.getElementById('length_bar')),
    recordAxisCount: recordGroups.filter(
      (group) => Boolean(group.querySelector(':scope > line'))
    ).length
  };
});

test('coordinate scale visibility follows simple controls and explicit Circular Ticks slots', async ({ page }) => {
  test.setTimeout(300000);
  const genbank = readFileSync(genbankPath, 'utf8');

  await page.goto(`${baseUrl}/gbdraw/web/index.html`, { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await page.waitForFunction(
    () => window.__GBDRAW_APP__?.diagramGenerationWorkerReady === true,
    null,
    { timeout: 180000 }
  );

  await page.evaluate(async (genbankText) => {
    const app = window.__GBDRAW_APP__;
    app.mode = 'circular';
    app.cInputType = 'gb';
    app.files.c_gb = new File([genbankText], 'HmmtDNA.gbk', {
      type: 'text/plain',
      lastModified: 1
    });
    app.form.show_scale = true;
    app.form.labels_mode = 'none';
    app.setCircularTrackSlotsEnabled(false);
    await window.Vue.nextTick();
  }, genbank);

  expect(await runDiagram(page)).toEqual({
    result: { status: 'ok' },
    errorSummary: '',
    errorDetails: []
  });
  expect(await inspectResultSvg(page)).toMatchObject({
    axis: true,
    ticks: true
  });

  await page.evaluate(() => {
    window.__GBDRAW_APP__.form.show_scale = false;
  });
  expect(await runDiagram(page)).toEqual({
    result: { status: 'ok' },
    errorSummary: '',
    errorDetails: []
  });
  expect(await inspectResultSvg(page)).toMatchObject({
    axis: true,
    ticks: false
  });

  await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    app.setCircularTrackSlotsEnabled(true);
    await window.Vue.nextTick();
  });
  const circularAxisCard = page.locator('summary').filter({ hasText: 'Axis & Scale' });
  await circularAxisCard.click();
  await expect(page.getByLabel('Show Coordinate Scale (Circular)')).toBeDisabled();
  await expect(page.locator('[data-scale-visibility-note]')).toContainText(
    'Use an enabled Ticks slot'
  );
  await expect(page.getByLabel('Circular scale interval')).toBeEnabled();

  expect(await runDiagram(page)).toEqual({
    result: { status: 'ok' },
    errorSummary: '',
    errorDetails: []
  });
  expect(await inspectResultSvg(page)).toMatchObject({
    axis: true,
    ticks: true
  });

  await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const tick = app.adv.circular_track_slots.find((slot) => slot.renderer === 'ticks');
    app.setCircularTrackSlotEnabled(tick, false);
    await window.Vue.nextTick();
  });
  await expect(page.getByLabel('Circular scale interval')).toBeDisabled();

  await page.evaluate(async (genbankText) => {
    const app = window.__GBDRAW_APP__;
    app.mode = 'linear';
    app.lInputType = 'gb';
    app.setLinearSeqPrimaryFile(0, 'gb', new File([genbankText], 'HmmtDNA.gbk', {
      type: 'text/plain',
      lastModified: 2
    }));
    Object.assign(app.form, {
      show_scale: true,
      scale_style: 'bar',
      show_gc: false,
      show_skew: false,
      show_depth: false,
      show_labels_linear: 'none'
    });
    app.setLinearTrackSlotsEnabled(false);
    await window.Vue.nextTick();
  }, genbank);

  expect(await runDiagram(page)).toEqual({
    result: { status: 'ok' },
    errorSummary: '',
    errorDetails: []
  });
  const visibleLinearScale = await inspectResultSvg(page);
  expect(visibleLinearScale.lengthBar).toBe(true);
  expect(visibleLinearScale.recordAxisCount).toBeGreaterThan(0);

  await page.evaluate(() => {
    window.__GBDRAW_APP__.form.show_scale = false;
  });
  const linearAxisCard = page.locator('summary').filter({ hasText: 'Axis & Scale' });
  await linearAxisCard.click();
  await expect(page.getByLabel('Linear scale style')).toBeDisabled();
  await expect(page.getByLabel('Axis stroke color mode')).toBeEnabled();

  expect(await runDiagram(page)).toEqual({
    result: { status: 'ok' },
    errorSummary: '',
    errorDetails: []
  });
  const hiddenLinearScale = await inspectResultSvg(page);
  expect(hiddenLinearScale.lengthBar).toBe(false);
  expect(hiddenLinearScale.recordAxisCount).toBeGreaterThan(0);
});
