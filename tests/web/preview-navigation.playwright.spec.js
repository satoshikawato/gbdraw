const { test, expect } = require('@playwright/test');
const { join } = require('node:path');
const { writeFile } = require('node:fs/promises');
const { openApp, getDiagramWorkerActivity } = require('./helpers/app-lifecycle.cjs');

const geometry = (page) => page.evaluate(() => {
  const app = window.__GBDRAW_APP__;
  const wrapper = app.svgContainer;
  const svg = wrapper.querySelector('svg');
  const container = app.canvasContainerRef;
  const matrix = svg.getScreenCTM();
  const point = (x, y) => {
    const p = new DOMPoint(x, y).matrixTransform(matrix);
    return { x: p.x, y: p.y };
  };
  const bounds = wrapper.getBoundingClientRect();
  return {
    zoom: app.zoom, pan: { ...app.canvasPan },
    scroll: { x: container.scrollLeft, y: container.scrollTop },
    points: [point(0, 0), point(100, 100)],
    anchor: { x: bounds.x + bounds.width / 2, y: bounds.y },
    panning: app.isPanning
  };
});

const settle = async (page) => {
  await expect.poll(() => page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const transform = new DOMMatrix(getComputedStyle(app.svgContainer).transform);
    return Math.max(Math.abs(transform.a - app.zoom),
      Math.abs(transform.e - app.canvasPan.x), Math.abs(transform.f - app.canvasPan.y));
  })).toBeLessThan(0.001);
};

const backgroundPoint = (page) => page.evaluate(() => {
  const app = window.__GBDRAW_APP__;
  const container = app.canvasContainerRef;
  const rect = container.getBoundingClientRect();
  const svg = app.svgContainer.querySelector('svg');
  for (let y = rect.bottom - 160; y > rect.y + 30; y -= 30) {
    for (let x = rect.x + 50; x < rect.right - 40; x += 60) {
      const target = document.elementFromPoint(x, y);
      if (target === container || target === app.svgContainer || target === svg) return { x, y };
    }
  }
  throw new Error('No visible non-editing background point.');
});

const move = async (page, point, dx, dy) => {
  await page.mouse.move(point.x, point.y);
  await page.mouse.down();
  await page.mouse.move(point.x + dx, point.y + dy, { steps: 10 });
  await page.mouse.up();
  await settle(page);
};

const expectTranslation = (before, after, dx, dy) => {
  for (let i = 0; i < before.points.length; i += 1) {
    expect(after.points[i].x - before.points[i].x).toBeCloseTo(dx, 1);
    expect(after.points[i].y - before.points[i].y).toBeCloseTo(dy, 1);
  }
  expect(after.panning).toBe(false);
};

test('@pr-smoke background pan preserves screen displacement after text selection and zoom', async ({ page }, testInfo) => {
  test.setTimeout(180000);
  await openApp(page);
  page.on('dialog', (dialog) => dialog.dismiss());
  await page.locator('input[type=file][accept*="application/json"][accept*="application/gzip"]')
    .setInputFiles(join(process.cwd(), 'gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json'));
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.results.length)).toBe(1);
  await settle(page);
  const committed = await page.evaluate(async () => {
    const { getCommittedCanonicalSession } = await import('./js/services/config.js');
    return JSON.stringify(getCommittedCanonicalSession());
  });
  const evidence = [];
  for (const zoom of [0.3, 1, 2.8]) {
    await page.getByRole('button', { name: 'Reset layout', exact: true }).click();
    const direction = zoom < 1 ? 'Zoom out' : 'Zoom in';
    for (let i = 0; i < Math.round(Math.abs(zoom - 1) * 10); i += 1) {
      await page.getByRole('button', { name: direction, exact: true }).click();
    }
    await settle(page);
    const before = await geometry(page);
    const point = await backgroundPoint(page);
    await move(page, point, 60, 80);
    const after = await geometry(page);
    evidence.push({ zoom, point, before, after });
    await writeFile(testInfo.outputPath('navigation.json'), JSON.stringify(evidence, null, 2));
    expectTranslation(before, after, 60, 80);
    // Preserve the existing top-center zoom anchor, including after panning.
    await page.mouse.move(point.x, point.y);
    await page.mouse.wheel(0, -100);
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.zoom)).toBeCloseTo(zoom + 0.1, 5);
    await settle(page);
    const zoomed = await geometry(page);
    expect(zoomed.anchor.x).toBeCloseTo(after.anchor.x, 1);
    expect(zoomed.anchor.y).toBeCloseTo(after.anchor.y, 1);
    await page.mouse.wheel(0, 100);
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.zoom)).toBeCloseTo(zoom, 5);
    await settle(page);
    expectTranslation(after, await geometry(page), 0, 0);
  }
  await page.getByRole('button', { name: 'Reset layout', exact: true }).click();
  await settle(page);
  expect(await geometry(page)).toMatchObject({ zoom: 1, pan: { x: 0, y: 0 }, panning: false });
  expect(await page.evaluate(async () => {
    const { getCommittedCanonicalSession } = await import('./js/services/config.js');
    return JSON.stringify(getCommittedCanonicalSession());
  })).toBe(committed);
  expect(await getDiagramWorkerActivity(page)).toMatchObject({ constructions: 0 });
});

test('@pr-smoke preview pan leaves feature and match gestures available', async ({ page }) => {
  test.setTimeout(180000);
  // Keep the wide comparison's editing targets visible at its existing zoom anchor.
  await page.setViewportSize({ width: 2520, height: 1327 });
  await openApp(page);
  page.on('dialog', (dialog) => dialog.dismiss());
  await page.locator('input[type=file][accept*="application/json"][accept*="application/gzip"]')
    .setInputFiles(join(process.cwd(), 'gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json'));
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.results.length)).toBe(1);
  for (let i = 0; i < 7; i += 1) await page.getByRole('button', { name: 'Zoom out', exact: true }).click();
  await settle(page);

  for (const selector of ['[data-gbdraw-feature-id]', '[data-gbdraw-pairwise-match-id]']) {
    const point = await page.evaluate((selector) => {
      for (const element of window.__GBDRAW_APP__.svgContainer.querySelectorAll(selector)) {
        const bounds = element.getBoundingClientRect();
        const x = bounds.x + bounds.width / 2;
        const y = bounds.y + bounds.height / 2;
        if (document.elementFromPoint(x, y)?.closest(selector)) return { x, y };
      }
      throw new Error(`No visible gesture target: ${selector}`);
    }, selector);
    const before = await geometry(page);
    await move(page, point, 10, 10);
    expectTranslation(before, await geometry(page), 0, 0);
    await page.mouse.click(point.x, point.y);
    const selectionKey = selector.includes('pairwise') ? 'clickedPairwiseMatch' : 'clickedFeature';
    await expect.poll(() => page.evaluate((key) => Boolean(window.__GBDRAW_APP__[key]), selectionKey)).toBe(true);
    await page.keyboard.press('Escape');
  }

  await page.locator('.drawer-toggle').click();
  await expect(page.locator('.right-drawer')).toBeVisible();
  const before = await geometry(page);
  await move(page, await backgroundPoint(page), 30, 20);
  expectTranslation(before, await geometry(page), 30, 20);
  await page.getByRole('button', { name: 'Reset zoom', exact: true }).focus();
  await page.keyboard.press('Enter');
  await settle(page);
  expect(await geometry(page)).toMatchObject({ zoom: 1, pan: { x: 30, y: 20 } });
  await page.getByRole('button', { name: 'Reset layout', exact: true }).click();
  await settle(page);
  expect(await geometry(page)).toMatchObject({ zoom: 1, pan: { x: 0, y: 0 }, panning: false });
  expect(await getDiagramWorkerActivity(page)).toMatchObject({ constructions: 0 });
});
