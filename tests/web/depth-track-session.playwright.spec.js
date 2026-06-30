const { test, expect } = require('@playwright/test');
const { createReadStream, existsSync } = require('node:fs');
const { createServer } = require('node:http');
const { extname, join, normalize, resolve, sep } = require('node:path');

const repoRoot = resolve(process.env.GBDRAW_REPO || process.cwd());
const sessionPath = join(repoRoot, 'tests/test_inputs/2026-06-16_wssv.gbdraw-session.json');
const bgcSessionPath = join(repoRoot, 'tests/test_inputs/BGC0000708-BGC0000713.gbdraw-session.json');

const contentTypes = {
  '.html': 'text/html; charset=utf-8',
  '.js': 'text/javascript; charset=utf-8',
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
      const requestedPath = normalize(decodeURIComponent(url.pathname)).replace(/^(\.\.(?:\/|\\|$))+/, '');
      const filePath = resolve(repoRoot, requestedPath.replace(/^[/\\]+/, ''));
      if (!filePath.startsWith(`${repoRoot}${sep}`) && filePath !== repoRoot) {
        response.writeHead(403);
        response.end('Forbidden');
        return;
      }
      const finalPath = filePath === repoRoot ? join(repoRoot, 'gbdraw/web/index.html') : filePath;
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
      const address = server.address();
      baseUrl = `http://127.0.0.1:${address.port}`;
      resolveServer();
    });
  });
});

test.afterAll(async () => {
  await new Promise((resolveClose) => server.close(resolveClose));
});

test('WSSV depth session removes stale circular depth metadata and slots', async ({ page }) => {
  page.on('dialog', (dialog) => dialog.accept());
  await page.goto(`${baseUrl}/gbdraw/web/index.html`, { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);

  await page.locator('input[accept=".json"]').first().setInputFiles(sessionPath);
  await page.waitForFunction(() => {
    const app = window.__GBDRAW_APP__;
    return Array.isArray(app?.files?.c_depth) && app.files.c_depth.length === 2;
  });

  const loaded = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    return {
      names: app.files.c_depth.map((file) => file?.name || ''),
      labels: app.adv.depth_tracks.map((track) => track?.label || ''),
      slots: app.adv.circular_track_slots
        .filter((slot) => slot?.renderer === 'depth')
        .map((slot) => ({
          id: slot.id,
          enabled: slot.enabled !== false,
          trackIndex: slot.params?.track_index,
          legendLabel: slot.params?.legend_label
        }))
    };
  });

  expect(loaded.names).toHaveLength(2);
  expect(loaded.labels).toHaveLength(2);
  expect(loaded.labels[1]).toContain('SRR14027721');
  expect(loaded.slots.filter((slot) => slot.enabled).map((slot) => slot.trackIndex)).toEqual([0, 1]);
  expect(loaded.slots.filter((slot) => slot.enabled).map((slot) => slot.legendLabel)).toEqual(loaded.labels);

  await page.evaluate(() => window.__GBDRAW_APP__.removeCircularDepthTrack(1));
  await page.waitForFunction(() => {
    const app = window.__GBDRAW_APP__;
    return Array.isArray(app.files.c_depth) &&
      app.files.c_depth.length === 1 &&
      app.adv.circular_track_slots
        .filter((slot) => slot?.renderer === 'depth' && slot.enabled !== false)
        .every((slot) => Number(slot.params?.track_index) < 1);
  });

  const afterRemoval = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    return {
      names: app.files.c_depth.map((file) => file?.name || ''),
      labels: app.adv.depth_tracks.map((track) => track?.label || ''),
      slots: app.adv.circular_track_slots
        .filter((slot) => slot?.renderer === 'depth')
        .map((slot) => ({
          id: slot.id,
          enabled: slot.enabled !== false,
          trackIndex: slot.params?.track_index,
          legendLabel: slot.params?.legend_label
        }))
    };
  });

  expect(afterRemoval.names).toHaveLength(1);
  expect(afterRemoval.labels[0]).toContain('SRR14027705');
  expect(afterRemoval.slots.filter((slot) => slot.enabled).map((slot) => slot.trackIndex)).toEqual([0]);
  expect(afterRemoval.slots.some((slot) => Number(slot.trackIndex) >= 1)).toBe(false);
  expect(afterRemoval.slots.some((slot) => String(slot.legendLabel || '').includes('SRR14027712'))).toBe(false);
});

test('BGC session keeps restored feature metadata selectable in the preview', async ({ page }) => {
  page.on('dialog', (dialog) => dialog.accept());
  await page.goto(`${baseUrl}/gbdraw/web/index.html`, { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await page.evaluate(() => {
    window.__GBDRAW_APP__.pyodideReady = true;
  });

  await page.locator('input[accept=".json"]').first().setInputFiles(bgcSessionPath);
  await page.waitForFunction(() => window.__GBDRAW_APP__?.results?.length > 0);
  await page.waitForTimeout(250);

  const summary = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const svg =
      document.querySelector('[data-gbdraw-feature-id]')?.ownerSVGElement ||
      document.querySelector('svg');
    const renderedIds = Array.from(
      svg?.querySelectorAll('[data-gbdraw-feature-id], path[id^="f"], polygon[id^="f"], rect[id^="f"]') || []
    )
      .map((el) => el.getAttribute('data-gbdraw-feature-id') || el.id || '')
      .filter(Boolean);
    const uniqueRenderedIds = Array.from(new Set(renderedIds));
    const featureIds = new Set(
      (Array.isArray(app.extractedFeatures) ? app.extractedFeatures : [])
        .map((feature) => String(feature?.svg_id || '').trim())
        .filter(Boolean)
    );
    return {
      extractedCount: Array.isArray(app.extractedFeatures) ? app.extractedFeatures.length : 0,
      uniqueRenderedCount: uniqueRenderedIds.length,
      matchedCount: uniqueRenderedIds.filter((id) => featureIds.has(id)).length
    };
  });

  expect(summary.extractedCount).toBeGreaterThan(0);
  expect(summary.uniqueRenderedCount).toBeGreaterThan(0);
  expect(summary.matchedCount).toBe(summary.uniqueRenderedCount);

  const target = await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const svg =
      document.querySelector('[data-gbdraw-feature-id]')?.ownerSVGElement ||
      document.querySelector('svg');
    const featureIds = new Set(
      (Array.isArray(app.extractedFeatures) ? app.extractedFeatures : [])
        .map((feature) => String(feature?.svg_id || '').trim())
        .filter(Boolean)
    );
    const element = Array.from(
      svg.querySelectorAll('[data-gbdraw-feature-id], path[id^="f"], polygon[id^="f"], rect[id^="f"]')
    ).find((candidate) => featureIds.has(candidate.getAttribute('data-gbdraw-feature-id') || candidate.id || ''));
    element.scrollIntoView({ block: 'center', inline: 'center' });
    await new Promise((resolve) => requestAnimationFrame(() => requestAnimationFrame(resolve)));
    const rect = element.getBoundingClientRect();
    return {
      id: element.getAttribute('data-gbdraw-feature-id') || element.id || '',
      x: rect.left + rect.width / 2,
      y: rect.top + rect.height / 2
    };
  });

  await page.mouse.click(target.x, target.y);
  await page.waitForFunction(
    (featureId) => window.__GBDRAW_APP__?.clickedFeature?.svg_id === featureId,
    target.id
  );
  await expect(page.locator('.feature-popup')).toBeVisible();
});
