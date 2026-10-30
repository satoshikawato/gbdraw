// OV-1036, OV-1038 (R1, PD-OI-037 rev 2 / OIC-024): the timing notes that
// `SETTING_APPLICATION` renders. tests/web/setting-application.test.mjs checks
// that every settings control sits in a row's scope; this checks the rendered
// notes the user sees.
const { test, expect } = require('@playwright/test');
const { openApp } = require('./helpers/app-lifecycle.cjs');
const { showDrawing } = require('./helpers/drawing-bar.cjs');
const { load, seeds } = require('./helpers/mode-transition.cjs');


const openLayout = async (page) => {
  const details = page.locator('[data-setting="layout"] > details');
  if (!await details.evaluate((element) => element.open)) await details.locator(':scope > summary').click();
  return page.locator('[data-setting="layout"]');
};

test('the Layout card says its settings apply on Generate, the track toggles are live, and the stack names its live switches (OV-1036, PD-OI-104 rev 5)', async ({ page }) => {
  await openApp(page);
  for (const [mode, toggles] of [
    ['linear', ['Show GC Content', 'Show GC Skew', 'Show Depth']],
    ['circular', ['Hide GC Content', 'Hide GC Skew', 'Show Depth']]
  ]) {
    await showDrawing(page, { kind: mode });
    const layout = await openLayout(page);
    await expect(layout.locator(':scope > details > [data-applies-on-generate="layout"]')).toHaveText('Applies on Generate');
    const visibility = layout.locator('[data-setting="track-visibility"]');
    await expect(visibility).toHaveCount(1);
    await expect(visibility.locator('.apply-note').first()).toHaveText('Live edit');
    await expect(visibility.locator('[data-applies-on-generate="track-visibility"]')).toHaveCount(0);
    await expect(visibility.locator('label')).toHaveText(toggles.map((name) => new RegExp(name)));
    await expect(layout.locator('[data-setting="track-slots"] [data-applies-on-generate="track-slots"]').first())
      .toHaveText(/Applies on Generate, except turning a GC content, GC skew or Depth row on or off/);
  }
});

const selectPaintedFeature = async (page, body) => {
  await body.scrollIntoViewIfNeeded();
  const point = await body.evaluate((element) => {
    const box = element.getBBox();
    for (let y = 1; y < 30; y += 1) for (let x = 1; x < 30; x += 1) {
      const local = new DOMPoint(box.x + box.width * x / 30, box.y + box.height * y / 30);
      if (!element.isPointInFill(local)) continue;
      const screen = local.matrixTransform(element.getScreenCTM());
      if (document.elementFromPoint(screen.x, screen.y)?.getAttribute('data-gbdraw-feature-id') === element.getAttribute('data-gbdraw-feature-id')) return { x: screen.x, y: screen.y };
    }
    throw new Error('No visible painted feature point.');
  });
  await page.keyboard.down('Control');
  await page.mouse.click(point.x, point.y);
  await page.keyboard.up('Control');
};

test('the selection toolbar placement select has the note of the popup Feature placement (OV-1038)', async ({ browser }) => {
  const page = await load(browser, seeds.linear);
  try {
    const ids = await page.evaluate(() => window.__GBDRAW_APP__.extractedFeatures.map((feature) => feature.svg_id).slice(0, 2));
    for (const id of ids) await selectPaintedFeature(page, page.locator(`[data-gbdraw-feature-id="${id}"]`).first());
    const toolbar = page.locator('.feature-selection-toolbar');
    await expect(toolbar.locator('[aria-label="Selected feature count"]')).toContainText('2 selected');
    const placement = toolbar.locator('[data-setting="feature-placement"]');
    await expect(placement.locator('select[aria-label="Selected feature placements"]')).toBeVisible();
    await expect(placement.locator('[data-applies-on-generate="feature-placement"]')).toHaveText('Applies on Generate');
  } finally {
    await page.context().close();
  }
});
