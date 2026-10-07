const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { join } = require('node:path');
const { openApp, generateAndWaitForResult } = require('./helpers/app-lifecycle.cjs');

// UI-11 (Web GUI re-audit 2026-10-05, Owner decision 2026-10-07): the preview
// toolbar's Fit button sets the zoom so the whole diagram lies inside the
// preview frame, centred. Generate keeps its first view.
const repoRoot = process.cwd();
const hmmt = readFileSync(join(repoRoot, 'tests/test_inputs/HmmtDNA.gbk'), 'utf8');
const linearSession = join(repoRoot, 'gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json');

const fitButton = (page) => page.getByRole('button', { name: 'Fit to preview', exact: true });

// The computed transform equals the state and stays still for two frames.
const settle = async (page) => {
  await expect.poll(() => page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const transform = new DOMMatrix(getComputedStyle(app.svgContainer).transform);
    return Math.max(Math.abs(transform.a - app.zoom),
      Math.abs(transform.e - app.canvasPan.x), Math.abs(transform.f - app.canvasPan.y));
  })).toBeLessThan(0.001);
  await expect.poll(() => page.evaluate(() => new Promise((resolve) => {
    const sample = () => {
      const box = window.__GBDRAW_APP__.svgContainer.querySelector('svg').getBoundingClientRect();
      return [box.x, box.y, box.width, box.height];
    };
    requestAnimationFrame(() => {
      const before = sample();
      requestAnimationFrame(() => resolve(before.every((value, index) => value === sample()[index])));
    });
  }))).toBe(true);
};

// Margins of the rendered SVG inside the preview frame's visible client area.
const fitGeometry = (page) => page.evaluate(() => {
  const app = window.__GBDRAW_APP__;
  const container = app.canvasContainerRef;
  const frame = container.getBoundingClientRect();
  const box = app.svgContainer.querySelector('svg').getBoundingClientRect();
  const left = frame.left + container.clientLeft;
  const top = frame.top + container.clientTop;
  return {
    zoom: app.zoom,
    pan: { ...app.canvasPan },
    frame: { width: container.clientWidth, height: container.clientHeight },
    unscaled: { width: box.width / app.zoom, height: box.height / app.zoom },
    margins: {
      left: box.left - left,
      right: left + container.clientWidth - box.right,
      top: box.top - top,
      bottom: top + container.clientHeight - box.bottom
    }
  };
});

const expectFitted = (geometry, label) => {
  const detail = `${label}: ${JSON.stringify(geometry)}`;
  for (const margin of Object.values(geometry.margins)) expect(margin, detail).toBeGreaterThanOrEqual(-1);
  expect(Math.abs(geometry.margins.left - geometry.margins.right), detail).toBeLessThanOrEqual(2);
  expect(Math.abs(geometry.margins.top - geometry.margins.bottom), detail).toBeLessThanOrEqual(2);
  // The next zoom step would no longer fit on the limiting axis.
  const largest = Math.min(geometry.frame.width / geometry.unscaled.width,
    geometry.frame.height / geometry.unscaled.height);
  expect(geometry.zoom + 0.1, detail).toBeGreaterThan(largest);
};

const fit = async (page) => {
  await expect(fitButton(page)).toBeVisible();
  await fitButton(page).click();
  await settle(page);
  return fitGeometry(page);
};

// A point that pans the canvas instead of starting a feature or label edit.
const backgroundPoint = (page) => page.evaluate(() => {
  const app = window.__GBDRAW_APP__;
  const rect = app.canvasContainerRef.getBoundingClientRect();
  const svg = app.svgContainer.querySelector('svg');
  for (let y = rect.bottom - 20; y > rect.y + 20; y -= 20) {
    for (let x = rect.x + 20; x < rect.right - 20; x += 40) {
      const target = document.elementFromPoint(x, y);
      if (target === app.canvasContainerRef || target === app.svgContainer || target === svg) return { x, y };
    }
  }
  throw new Error('No visible non-editing background point.');
});

test('Fit shows the whole Circular diagram centred at 1280x800 and 1920x1080, also after zoom and pan', async ({ page }) => {
  test.setTimeout(240000);
  await page.setViewportSize({ width: 1280, height: 800 });
  await openApp(page);
  const chooser = page.waitForEvent('filechooser');
  await page.getByRole('button', { name: 'Choose GenBank/DDBJ File', exact: true }).click();
  await (await chooser).setFiles({ name: 'HmmtDNA.gbk', mimeType: 'text/plain', buffer: Buffer.from(hmmt) });
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.circularRecordList.length),
    { timeout: 60000 }).toBeGreaterThan(0);
  await generateAndWaitForResult(page);
  await settle(page);
  // Generate keeps its first view; only the button fits.
  const firstView = await fitGeometry(page);
  expect(firstView.zoom).toBe(1);
  expect(firstView.margins.bottom).toBeLessThan(-1);
  const undoBefore = await page.evaluate(() => window.__GBDRAW_APP__.canUndoHistory);

  const fitted = await fit(page);
  expectFitted(fitted, '1280x800');
  // A view change: no History step.
  expect(await page.evaluate(() => window.__GBDRAW_APP__.canUndoHistory)).toBe(undoBefore);
  await expect(page.getByRole('button', { name: 'Reset zoom', exact: true }))
    .toHaveText(`${Math.round(fitted.zoom * 100)}%`);

  for (let step = 0; step < 3; step += 1) await page.getByRole('button', { name: 'Zoom in', exact: true }).click();
  const point = await backgroundPoint(page);
  await page.mouse.move(point.x, point.y);
  await page.mouse.down();
  await page.mouse.move(point.x - 120, point.y - 90, { steps: 8 });
  await page.mouse.up();
  await settle(page);
  const moved = await fitGeometry(page);
  expect(moved.zoom).toBeCloseTo(fitted.zoom + 0.3, 5);
  expect(moved.pan).not.toEqual(fitted.pan);
  const refitted = await fit(page);
  expectFitted(refitted, '1280x800 after zoom and pan');
  expect(refitted.zoom).toBe(fitted.zoom);
  expect(refitted.pan.x).toBeCloseTo(fitted.pan.x, 1);
  expect(refitted.pan.y).toBeCloseTo(fitted.pan.y, 1);

  await page.setViewportSize({ width: 1920, height: 1080 });
  expectFitted(await fit(page), '1920x1080');
});

test('Fit shows the whole Linear comparison centred', async ({ page }) => {
  test.setTimeout(180000);
  await page.setViewportSize({ width: 1280, height: 800 });
  await openApp(page);
  page.on('dialog', (dialog) => dialog.dismiss());
  await page.locator('input[type=file][accept*="application/json"][accept*="application/gzip"]')
    .setInputFiles(linearSession);
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.results.length)).toBe(1);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.mode)).toBe('linear');
  await settle(page);
  expectFitted(await fit(page), 'Linear 1280x800');
});
