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

// Margins of the rendered SVG inside the visible preview frame: the canvas
// client area left of an open Editor drawer that lies over it.
const fitGeometry = (page) => page.evaluate(() => {
  const app = window.__GBDRAW_APP__;
  const container = app.canvasContainerRef;
  const frame = container.getBoundingClientRect();
  const box = app.svgContainer.querySelector('svg').getBoundingClientRect();
  const left = frame.left + container.clientLeft;
  const top = frame.top + container.clientTop;
  const bottom = top + container.clientHeight;
  let right = left + container.clientWidth;
  const drawer = document.querySelector('.right-drawer[aria-hidden="false"]')?.getBoundingClientRect();
  const covered = drawer && drawer.top < bottom && drawer.bottom > top && drawer.left < right
    ? right - drawer.left : 0;
  right -= covered;
  return {
    zoom: app.zoom,
    pan: { ...app.canvasPan },
    covered,
    frame: { width: right - left, height: container.clientHeight },
    unscaled: { width: box.width / app.zoom, height: box.height / app.zoom },
    margins: { left: box.left - left, right: right - box.right, top: box.top - top, bottom: bottom - box.bottom }
  };
});

const expectFitted = (geometry, label) => {
  const detail = `${label}: ${JSON.stringify(geometry)}`;
  for (const margin of Object.values(geometry.margins)) expect(margin, detail).toBeGreaterThanOrEqual(-1);
  expect(Math.abs(geometry.margins.left - geometry.margins.right), detail).toBeLessThanOrEqual(2);
  expect(Math.abs(geometry.margins.top - geometry.margins.bottom), detail).toBeLessThanOrEqual(2);
  // Fit fills the frame: a whole percent, and on the limiting axis only the
  // 8 px margins and less than one 1% step stay free.
  expect(Math.round(geometry.zoom * 100) / 100, detail).toBe(geometry.zoom);
  const widthLimited = geometry.frame.width / geometry.unscaled.width
    < geometry.frame.height / geometry.unscaled.height;
  const free = widthLimited ? geometry.margins.left + geometry.margins.right
    : geometry.margins.top + geometry.margins.bottom;
  const size = widthLimited ? geometry.unscaled.width : geometry.unscaled.height;
  expect(free, detail).toBeLessThanOrEqual(16 + 0.01 * size + 2);
};

const appZoom = (page) => page.evaluate(() => window.__GBDRAW_APP__.zoom);

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

  const zoomIn = page.getByRole('button', { name: 'Zoom in', exact: true });
  // From a whole-percent Fit, the buttons step back onto the 0.1 grid.
  let expected = fitted.zoom;
  for (let step = 0; step < 3; step += 1) {
    await zoomIn.click();
    expected = Math.round((expected + 0.1) * 10) / 10;
  }
  const point = await backgroundPoint(page);
  await page.mouse.move(point.x, point.y);
  await page.mouse.down();
  await page.mouse.move(point.x - 120, point.y - 90, { steps: 8 });
  await page.mouse.up();
  await settle(page);
  const moved = await fitGeometry(page);
  expect(moved.zoom).toBe(expected);
  expect(moved.pan).not.toEqual(fitted.pan);
  const refitted = await fit(page);
  expectFitted(refitted, '1280x800 after zoom and pan');
  expect(refitted.zoom).toBe(fitted.zoom);
  expect(refitted.pan.x).toBeCloseTo(fitted.pan.x, 1);
  expect(refitted.pan.y).toBeCloseTo(fitted.pan.y, 1);

  await page.setViewportSize({ width: 1920, height: 1080 });
  expectFitted(await fit(page), '1920x1080');
});

const openLinearSession = async (page) => {
  await page.setViewportSize({ width: 1280, height: 800 });
  await openApp(page);
  page.on('dialog', (dialog) => dialog.dismiss());
  await page.locator('input[type=file][accept*="application/json"][accept*="application/gzip"]')
    .setInputFiles(linearSession);
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.results.length)).toBe(1);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.mode)).toBe('linear');
  await settle(page);
};

test('Fit shows the whole Linear comparison centred, left of an open Editor', async ({ page }) => {
  test.setTimeout(180000);
  await openLinearSession(page);
  const plain = await fit(page);
  expectFitted(plain, 'Linear 1280x800');
  expect(plain.covered).toBe(0);

  // The open Editor drawer lies over the right of the canvas; Fit uses the rest.
  await page.getByRole('button', { name: 'Editor', exact: true }).click();
  await expect(page.locator('.right-drawer')).toHaveAttribute('aria-hidden', 'false');
  await page.locator('.right-drawer').evaluate((drawer) => Promise.all(
    drawer.getAnimations().map((animation) => animation.finished.catch(() => {}))
  ));
  const beside = await fit(page);
  expect(beside.covered).toBeGreaterThan(300);
  expectFitted(beside, 'Linear 1280x800 with the Editor open');
  expect(beside.zoom).toBeLessThanOrEqual(plain.zoom);
});

// GX-02: the buttons keep the wheel's range and 0.1 steps.
test('Zoom in stops at 500% and Zoom out at 10%, on exact 0.1 steps', async ({ page }) => {
  test.setTimeout(180000);
  await openLinearSession(page);
  const zoomIn = page.getByRole('button', { name: 'Zoom in', exact: true });
  const zoomOut = page.getByRole('button', { name: 'Zoom out', exact: true });
  const resetZoom = page.getByRole('button', { name: 'Reset zoom', exact: true });
  await resetZoom.click();
  const steps = [];
  for (let step = 0; step < 45; step += 1) {
    await zoomIn.click();
    steps.push(await appZoom(page));
  }
  expect(steps.filter((value) => value !== Math.round(value * 10) / 10)).toEqual([]);
  expect(steps.slice(0, 3)).toEqual([1.1, 1.2, 1.3]);
  expect(steps.at(-1)).toBe(5);
  await expect(resetZoom).toHaveText('500%');
  for (let step = 0; step < 3; step += 1) await zoomOut.click();
  expect(await appZoom(page)).toBe(4.7);
  for (let step = 0; step < 50; step += 1) await zoomOut.click();
  expect(await appZoom(page)).toBe(0.1);
  await expect(resetZoom).toHaveText('10%');
});
