const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { join } = require('node:path');
const { openApp, generateAndWaitForResult, waitForSessionLoadAtRest } = require('./helpers/app-lifecycle.cjs');

// UI-11 (Web GUI re-audit 2026-10-05, Owner decision 2026-10-07): the preview
// toolbar's Fit button sets the zoom so the whole diagram lies inside the
// preview frame, centred. Generate keeps its first view.
const repoRoot = process.cwd();
const hmmt = readFileSync(join(repoRoot, 'tests/test_inputs/HmmtDNA.gbk'), 'utf8');
const linearSession = join(repoRoot, 'gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw');

const fitButton = (page) => page.getByRole('button', { name: 'Fit to preview', exact: true });

// The computed transform equals the state and stays still for two frames.
// Chromium serializes the computed matrix to six significant digits (a pan of
// -1045.2734375 px reads -1045.27) while the state and the layout keep the
// double, so each value is compared within half a unit of its own sixth
// significant digit: 0.005 px for a pan of 1000 to 9999 px, 5e-7 for a zoom
// below 1. A value on a rounding tie reads exactly half a unit off, and the
// subtraction adds float noise (1.0000000000331966 for the BGC Gallery pan,
// OV-2983), so the bound keeps a 1e-6 margin. A transition frame or a stale
// pan is off by far more.
const SERIALIZED_TIE_MARGIN = 1 + 1e-6;
const settle = async (page) => {
  await expect.poll(() => page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    const transform = new DOMMatrix(getComputedStyle(app.svgContainer).transform);
    const serializedPrecision = (value) => 0.5 * 10 ** (Math.floor(Math.log10(Math.abs(value) || 1)) - 5);
    return Math.max(...[[transform.a, app.zoom], [transform.e, app.canvasPan.x], [transform.f, app.canvasPan.y]]
      .map(([computed, state]) => Math.abs(computed - state) / serializedPrecision(state)));
  })).toBeLessThanOrEqual(SERIALIZED_TIE_MARGIN);
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
  await page.getByRole('button', { name: 'Add GenBank/DDBJ files', exact: true }).click();
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
  await page.locator('input[type=file][accept*=".gbdraw"]')
    .setInputFiles(linearSession);
  await waitForSessionLoadAtRest(page, { timeout: 60_000 });
  expect(await page.evaluate(() => window.__GBDRAW_APP__.results.length)).toBe(1);
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

// OV-2702: a file Session Load shows its figure fitted from the first frame, with no
// transition to the Fit. Every animation frame from the start of the Load until the Load
// has ended and the figure lies still records the surface's computed transform while the
// figure can be seen (the canvas has a size, at a phone width only once the Preview pane
// shows, and the surface is not hidden), whether the Load was pending, and whether the
// figure spilled over the canvas by more than a pixel. The sampling ends on the app's own
// state, not after a delay.
const sampleLoadFrames = (page) => page.evaluate(() => {
  const frames = [];
  let still = 0;
  window.__loadFrames = new Promise((resolve) => {
    const sample = () => {
      const app = window.__GBDRAW_APP__;
      const surface = app.svgContainer;
      const canvas = app.canvasContainerRef;
      if (surface && canvas?.offsetWidth && getComputedStyle(surface).visibility !== 'hidden') {
        const box = surface.querySelector('svg').getBoundingClientRect();
        const frame = canvas.getBoundingClientRect();
        const left = frame.left + canvas.clientLeft;
        const top = frame.top + canvas.clientTop;
        const spill = Math.max(left - box.left, box.right - left - canvas.clientWidth,
          top - box.top, box.bottom - top - canvas.clientHeight);
        frames.push({ transform: getComputedStyle(surface).transform, animations: surface.getAnimations().length,
          pending: Boolean(app.sessionImportPending), spills: spill > 1 });
      }
      const done = !app.sessionImportPending && app.results.length === 1 && surface && !surface.getAnimations().length;
      still = done ? still + 1 : 0;
      if (still >= 3) resolve(frames);
      else requestAnimationFrame(sample);
    };
    requestAnimationFrame(sample);
  });
});
const expectFittedFromFirstFrame = async (page, label) => {
  const frames = await page.evaluate(() => window.__loadFrames);
  const shown = await page.evaluate(() => getComputedStyle(window.__GBDRAW_APP__.svgContainer).transform);
  const detail = `${label}: ${JSON.stringify({ shown, frames })}`;
  expect(frames.length, detail).toBeGreaterThan(0);
  for (const { transform, animations } of frames) expect({ transform, animations }, detail).toEqual({ transform: shown, animations: 0 });
};

// R17-25 (OV-1809): a Session file Load with a Result fits the figure to the
// shown Preview at every width; Generate keeps 100 %.
const circularSession = join(repoRoot, 'gbdraw/web/gallery/sessions/HmmtDNA_ATskew.gbdraw');
for (const viewport of [{ width: 1280, height: 800 }, { width: 390, height: 844 }]) {
  test(`a Session file Load fits the figure to the Preview at ${viewport.width}x${viewport.height}`, async ({ page }) => {
    test.setTimeout(180000);
    await page.setViewportSize(viewport);
    await openApp(page);
    page.on('dialog', (dialog) => dialog.dismiss());
    await sampleLoadFrames(page);
    await page.locator('input[type=file][accept*=".gbdraw"]')
      .setInputFiles(circularSession);
    await waitForSessionLoadAtRest(page, { timeout: 60_000 });
    expect(await page.evaluate(() => window.__GBDRAW_APP__.results.length)).toBe(1);
    await settle(page);
    const loaded = await fitGeometry(page);
    expect(loaded.zoom, JSON.stringify(loaded)).toBeLessThan(1);
    expectFitted(loaded, `Load ${viewport.width}x${viewport.height}`);
    await expectFittedFromFirstFrame(page, `Load ${viewport.width}x${viewport.height}`);
  });
}

// OV-2703: a Session 27-44 file shows its stored figure at once while Python converts the
// file (DESIGN 4.3, PD-OI-044). That figure is fitted from its first frame too, and so is the
// Result the Load ends with: no frame shows a figure spilling over the canvas or animating.
// (The end of the Load may move the fitted figure a few pixels, with the notice that replaces
// the older-Session notice above the canvas.) At a phone width the Settings pane shows during
// the Load and the canvas has no room: a tap on the Preview tab while the Load is pending shows
// the figure fitted from its first frame, the stored figure of a Session 27-44 as well as the
// Result of a current Session, which stays hidden until the Load ends.
const pendingTapCases = [
  { session: 'tests/fixtures/sessions/HmmtDNA_ATskew.v40-schema5.json', viewport: { width: 1280, height: 800 } },
  { session: 'tests/fixtures/sessions/HmmtDNA_ATskew.v40-schema5.json', viewport: { width: 390, height: 844 }, tap: true },
  { session: 'gbdraw/web/gallery/sessions/HmmtDNA_ATskew.gbdraw', viewport: { width: 390, height: 844 }, tap: true }
];
for (const { session, viewport, tap } of pendingTapCases) {
  const older = session.includes('.v40-');
  const title = `${older ? 'a Session 27-44' : 'a current Session'} file Load shows its ${older ? 'stored figure' : 'figure'} fitted`
    + ` from the first frame at ${viewport.width}x${viewport.height}${tap ? ' after an early Preview tap' : ''}`;
  test(title, async ({ page }) => {
    test.setTimeout(240000);
    await page.setViewportSize(viewport);
    await openApp(page);
    page.on('dialog', (dialog) => dialog.dismiss());
    await sampleLoadFrames(page);
    await page.locator('input[type=file][accept*=".gbdraw"]')
      .setInputFiles(join(repoRoot, session));
    if (tap) {
      // The tap comes while the Load is pending: once the stored figure is drawn, else at once.
      await expect.poll(() => page.evaluate((waitForFigure) => {
        const app = window.__GBDRAW_APP__;
        return app.sessionImportPending && (!waitForFigure || app.results.length === 1);
      }, older), { intervals: [20] }).toBe(true);
      await page.locator('[data-phone-pane="preview"]').click();
    }
    const frames = await page.evaluate(() => window.__loadFrames);
    const detail = JSON.stringify(frames);
    // A Session 27-44's stored figure was shown while the Load was pending.
    if (older) expect(frames.some(({ pending }) => pending), detail).toBe(true);
    for (const { animations, spills } of frames) expect({ animations, spills }, detail).toEqual({ animations: 0, spills: false });
  });
}

// R17-25 covers a file Session Load with a Result: a legacy configuration Load
// brings no Result, so the generated figure the user zoomed keeps its zoom.
// (After a Session Load the held view refits to the room the notice changes,
// so this case starts from a Generate.) OV-2201: the R17-25 Fit that ran for a
// legacy Load too refitted it (zoom 1.3 -> 0.43).
test('a legacy configuration Load keeps the zoom of the figure already shown', async ({ page }) => {
  test.setTimeout(240000);
  await page.setViewportSize({ width: 1280, height: 800 });
  await openApp(page);
  page.on('dialog', (dialog) => dialog.dismiss());
  const chooser = page.waitForEvent('filechooser');
  await page.getByRole('button', { name: 'Add GenBank/DDBJ files', exact: true }).click();
  await (await chooser).setFiles({ name: 'HmmtDNA.gbk', mimeType: 'text/plain', buffer: Buffer.from(hmmt) });
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.circularRecordList.length),
    { timeout: 60000 }).toBeGreaterThan(0);
  await generateAndWaitForResult(page);
  await settle(page);
  const zoomIn = page.getByRole('button', { name: 'Zoom in', exact: true });
  for (let step = 0; step < 3; step += 1) await zoomIn.click();
  const zoomed = await appZoom(page);
  await page.locator('input[accept*=".gbdraw"]').setInputFiles({
    name: 'legacy-config.json', mimeType: 'application/json', buffer: Buffer.from(JSON.stringify({ form: {} }))
  });
  await expect(page.getByText('Configuration loaded from an older gbdraw.', { exact: false }).first()).toBeVisible();
  await settle(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.results.length)).toBe(1);
  expect(await appZoom(page)).toBe(zoomed);
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
