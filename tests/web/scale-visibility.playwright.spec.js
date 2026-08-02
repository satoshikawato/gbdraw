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

test('Arrow controls render in both modes and survive a session round trip', async ({ page }) => {
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
      lastModified: 31
    });
    app.adv.features.splice(0, app.adv.features.length, 'CDS');
    app.form.labels_mode = 'none';
    app.form.suppress_gc = true;
    app.form.suppress_skew = true;
    await window.Vue.nextTick();
  }, genbank);

  const featuresCard = page.locator('.card').filter({
    has: page.locator('summary').filter({ hasText: 'Features' })
  }).first();
  await featuresCard.locator('summary').click();
  const rendering = page.getByLabel('Rendering for CDS');
  const headRatio = page.getByLabel('Arrow head length ratio');
  const shaftRatio = page.getByLabel('Arrow shaft width ratio');
  expect(await rendering.locator('option').evaluateAll((options) => (
    options.map((option) => ({ value: option.value, text: option.textContent.trim() }))
  ))).toEqual([
    { value: 'arrow', text: 'Arrow' },
    { value: 'rectangle', text: 'Rectangle' },
    { value: 'underlay', text: 'Underlay' }
  ]);
  await expect(rendering).toHaveValue('arrow');
  await expect(headRatio).toHaveAttribute('min', '0.05');
  await expect(headRatio).toHaveAttribute('step', '0.05');
  expect(await headRatio.getAttribute('max')).toBeNull();
  await expect(shaftRatio).toHaveAttribute('min', '0.05');
  await expect(shaftRatio).toHaveAttribute('max', '1');
  await expect(shaftRatio).toHaveAttribute('step', '0.05');
  await expect(shaftRatio).toHaveValue('1');

  await headRatio.fill('1.25');
  await headRatio.press('ArrowUp');
  await expect(headRatio).toHaveValue('1.3');
  await headRatio.press('ArrowDown');
  await expect(headRatio).toHaveValue('1.25');
  await shaftRatio.fill('0.25');
  await shaftRatio.press('ArrowUp');
  await expect(shaftRatio).toHaveValue('0.3');
  await shaftRatio.press('ArrowDown');
  await expect(shaftRatio).toHaveValue('0.25');

  await headRatio.fill('0');
  const invalidRun = await runDiagram(page);
  expect(invalidRun.result).toEqual({ status: 'error' });
  expect(invalidRun.errorSummary).toContain(
    'Arrow head length ratio must be Auto or a positive finite number.'
  );
  await headRatio.fill('1.25');
  await shaftRatio.fill('0.25');
  await expect(rendering).toHaveValue('arrow');
  await expect(headRatio).toHaveValue('1.25');
  await expect(shaftRatio).toHaveValue('0.25');

  expect(await runDiagram(page)).toEqual({
    result: { status: 'ok' },
    errorSummary: '',
    errorDetails: []
  });
  expect(await page.evaluate(() => {
    const content = window.__GBDRAW_APP__.results?.[0]?.content || '';
    const svg = new DOMParser().parseFromString(content, 'image/svg+xml').documentElement;
    return svg.querySelectorAll('path[data-gbdraw-feature-part="block"]').length;
  })).toBeGreaterThan(0);

  await page.evaluate(async (genbankText) => {
    const app = window.__GBDRAW_APP__;
    app.mode = 'linear';
    app.lInputType = 'gb';
    app.setLinearSeqPrimaryFile(0, 'gb', new File([genbankText], 'HmmtDNA.gbk', {
      type: 'text/plain',
      lastModified: 32
    }));
    Object.assign(app.form, {
      show_gc: false,
      show_skew: false,
      show_depth: false,
      show_labels_linear: 'none'
    });
    await window.Vue.nextTick();
  }, genbank);

  expect(await page.evaluate(() => ({
    shape: window.__GBDRAW_APP__.adv.feature_shapes.CDS,
    head: window.__GBDRAW_APP__.adv.arrow_head_length_ratio,
    shaft: window.__GBDRAW_APP__.adv.arrow_shaft_width_ratio
  }))).toEqual({ shape: 'arrow', head: 1.25, shaft: 0.25 });
  expect(await runDiagram(page)).toEqual({
    result: { status: 'ok' },
    errorSummary: '',
    errorDetails: []
  });

  const linearGeometry = await page.evaluate(() => {
    const content = window.__GBDRAW_APP__.results?.[0]?.content || '';
    const svg = new DOMParser().parseFromString(content, 'image/svg+xml').documentElement;
    const numberPattern = '[+-]?(?:\\d+\\.?\\d*|\\.\\d+)(?:e[+-]?\\d+)?';
    const vertexPattern = new RegExp(`[ML]\\s*(${numberPattern})\\s*,\\s*(${numberPattern})`, 'gi');
    for (const path of svg.querySelectorAll('path[data-gbdraw-feature-part="block"]')) {
      const vertices = [...String(path.getAttribute('d') || '').matchAll(vertexPattern)]
        .map((match) => ({ x: Number(match[1]), y: Number(match[2]) }));
      if (vertices.length !== 7) continue;
      const fullWidth = Math.abs(vertices[4].y - vertices[2].y);
      return {
        d: path.getAttribute('d'),
        vertexCount: vertices.length,
        headRatio: Math.abs(vertices[3].x - vertices[2].x) / fullWidth,
        shaftRatio: Math.abs(vertices[6].y - vertices[0].y) / fullWidth
      };
    }
    return null;
  });
  expect(linearGeometry).not.toBeNull();
  expect(linearGeometry.vertexCount).toBe(7);
  expect(linearGeometry.headRatio).toBeCloseTo(1.25, 6);
  expect(linearGeometry.shaftRatio).toBeCloseTo(0.25, 6);

  await page.evaluate(() => {
    window.__GBDRAW_APP__.sessionTitle = 'arrow-browser-round-trip';
  });
  const sessionDownloadPromise = page.waitForEvent('download', { timeout: 120000 });
  expect((await page.evaluate(() => window.__GBDRAW_APP__.saveSessionWithTitle())).status)
    .toBe('saved');
  const sessionPath = await (await sessionDownloadPromise).path();

  await page.reload({ waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  const dialogPromise = page.waitForEvent('dialog', { timeout: 120000 });
  await page.locator('input[accept^=".json,"]').first().setInputFiles(sessionPath);
  const dialog = await dialogPromise;
  expect(dialog.message()).toBe('Session loaded successfully!');
  await dialog.accept();
  await page.waitForFunction(
    () => (
      window.__GBDRAW_APP__?.mode === 'linear'
      && window.__GBDRAW_APP__?.adv?.arrow_head_length_ratio === 1.25
      && window.__GBDRAW_APP__?.adv?.arrow_shaft_width_ratio === 0.25
    )
  );
  expect(await page.evaluate(() => ({
    mode: window.__GBDRAW_APP__.mode,
    shape: window.__GBDRAW_APP__.adv.feature_shapes.CDS,
    head: window.__GBDRAW_APP__.adv.arrow_head_length_ratio,
    shaft: window.__GBDRAW_APP__.adv.arrow_shaft_width_ratio,
    content: window.__GBDRAW_APP__.results?.[0]?.content || ''
  }))).toMatchObject({
    mode: 'linear',
    shape: 'arrow',
    head: 1.25,
    shaft: 0.25
  });
  expect(await page.evaluate((expectedPath) => (
    String(window.__GBDRAW_APP__.results?.[0]?.content || '').includes(expectedPath)
  ), linearGeometry.d)).toBe(true);

  expect(await page.evaluate(() => {
    const originalConfirm = window.confirm;
    window.confirm = () => true;
    try {
      window.__GBDRAW_APP__.resetSettings();
      return {
        shape: window.__GBDRAW_APP__.adv.feature_shapes.CDS,
        head: window.__GBDRAW_APP__.adv.arrow_head_length_ratio,
        shaft: window.__GBDRAW_APP__.adv.arrow_shaft_width_ratio
      };
    } finally {
      window.confirm = originalConfirm;
    }
  })).toEqual({ shape: 'arrow', head: null, shaft: 1.0 });
});
