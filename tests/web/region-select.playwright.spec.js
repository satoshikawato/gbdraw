// Select region (PD-OI-115): the preview toolbar toggle, a drag along a
// record's axis, the band and the Region popup, and its actions. Points are
// placed from the axis Python wrote (services/region-axis.js), so a drag
// lands where Python says a base is drawn.
const { test, expect } = require('@playwright/test');
const { gunzipSync } = require('node:zlib');
const { generateAndWaitForResult } = require('./helpers/app-lifecycle.cjs');
const { HMMT, loadSessionFile, openWithGenBank } = require('./helpers/audit-browser.cjs');
const { download, generate, load, seeds } = require('./helpers/mode-transition.cjs');

test.describe.configure({ retries: 0 });

/**
 * The screen point of display position `position` (bases from display index 0) of the shown Result's record
 * `recordIndex`, on the middle of its feature tracks (Linear: on row `row`).
 * @param {import('@playwright/test').Page} page
 */
const axisPoint = (page, recordIndex, position, row = 0) => page.evaluate(async ({ recordIndex, position, row }) => {
  const { state } = await import('./js/state.js');
  const { readRegionAxisRecords, regionAxisFrame } = await import('./js/services/region-axis.js');
  const record = readRegionAxisRecords(state.trackSlotResolvedGeometry.value, state.selectedResultIndex.value)
    .find((entry) => entry.recordIndex === recordIndex);
  const frame = regionAxisFrame(state.svgContainer.value.querySelector('svg'), recordIndex);
  const { axis } = record;
  let x;
  let y;
  if (axis.kind === 'circular') {
    const radius = (axis.bandInnerPx + axis.bandOuterPx) / 2;
    const angle = (axis.startAngleDeg + axis.degreesPerBp * position) * Math.PI / 180;
    x = radius * Math.cos(angle);
    y = radius * Math.sin(angle);
  } else {
    const shown = axis.rows[row];
    x = shown.xStartPx + (shown.xEndPx - shown.xStartPx) * (position - shown.displayStart) / (shown.displayEnd - shown.displayStart);
    y = (shown.bandTopPx + shown.bandBottomPx) / 2;
  }
  const point = new DOMPoint(x, y).matrixTransform(frame.getScreenCTM());
  return { x: point.x, y: point.y };
}, { recordIndex, position, row });

/** A drag through the given display positions of one record. @param {import('@playwright/test').Page} page */
const dragAlong = async (page, recordIndex, positions) => {
  const [first, ...rest] = positions;
  const start = await axisPoint(page, recordIndex, first);
  await page.mouse.move(start.x, start.y);
  await page.mouse.down();
  for (const position of rest) {
    const point = await axisPoint(page, recordIndex, position);
    await page.mouse.move(point.x, point.y, { steps: 4 });
  }
  await expect(page.locator('[data-region-select-readout]')).toBeVisible();
  await page.mouse.up();
};

const regionPopup = (page) => page.locator('[data-preview-item-popup="region"]');

/**
 * Toggles Select region on after a Generate; Layout edit turns off. The figure is fitted to the preview first, so
 * every point of the drags is on the screen.
 * @param {import('@playwright/test').Page} page
 */
const turnOn = async (page) => {
  await page.getByRole('button', { name: 'Fit to preview', exact: true }).click();
  // The fit animates the preview transform: the points are read once it has settled.
  let last = '';
  await expect.poll(async () => {
    const now = await page.evaluate(() => JSON.stringify(document.querySelector('.gbdraw-preview-surface').getBoundingClientRect()));
    const settled = now === last;
    last = now;
    return settled;
  }, { intervals: [250] }).toBe(true);
  const toggle = page.locator('[data-region-select-toggle]');
  await expect(toggle).toBeEnabled();
  await toggle.click();
  await expect(toggle).toHaveAttribute('aria-pressed', 'true');
  return toggle;
};

/** @param {import('@playwright/test').Page} page */
const regionRows = (page) => page.evaluate(() => JSON.parse(JSON.stringify(
  window.__GBDRAW_APP__.annotationSets.flatMap((set) => set.annotations
    .filter((item) => /^Region \d+$/.test(item.label))
    .map((item) => ({ set: set.id, id: item.id, label: item.label, target: item.target })))
)));

/** @param {import('@playwright/test').Page} page */
const undoLabels = (page) => page.evaluate(() => window.__GBDRAW_HISTORY__.timeline().undo.map((step) => step.label));

test('Circular: a loaded Result without axis disables Select region; after Generate a drag selects a span that Add to annotation set draws', async ({ browser }) => {
  test.setTimeout(420_000);
  const page = await load(browser, seeds.circular);
  try {
    const toggle = page.locator('[data-region-select-toggle]');
    // The Gallery Session was saved before Python wrote the axis: no Python starts at Load.
    await expect(toggle).toBeDisabled();
    await expect(toggle).toHaveAttribute('title', 'Generate to enable region selection');

    await generate(page);
    // Layout edit and Select region turn each other off.
    const layout = page.getByRole('button', { name: 'Toggle layout edit mode', exact: true });
    await layout.click();
    await expect(layout).toHaveAttribute('aria-pressed', 'true');
    await turnOn(page);
    await expect(layout).toHaveAttribute('aria-pressed', 'false');

    // Shift+drag keeps the feature marquee and selects no text (OV-2478).
    const a = await axisPoint(page, 0, 1000);
    const b = await axisPoint(page, 0, 2000);
    await page.keyboard.down('Shift');
    await page.mouse.move(a.x, a.y);
    await page.mouse.down();
    await page.mouse.move(b.x, b.y, { steps: 5 });
    expect(await page.evaluate(() => String(window.getSelection())), 'no native text selection').toBe('');
    await page.mouse.up();
    await page.keyboard.up('Shift');
    await expect(page.locator('[data-region-select-readout]')).toHaveCount(0);
    await page.keyboard.press('Escape');

    const record = await page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      return state.trackSlotResolvedGeometry.value.records[0].axis.source.drawnLength;
    });
    await dragAlong(page, 0, [4000.5, 5000.5, 6000.5]);
    const popup = regionPopup(page);
    await expect(popup).toBeVisible();
    await expect(page.locator('[data-region-select-band]')).toHaveCount(1);
    await expect(page.locator('[data-match-popup-title]')).toHaveText('NC_012920.1');
    await expect(page.locator('[data-preview-item-subtitle]')).toHaveText('4,001-6,001 (2,001 bp)');
    await expect(popup.locator('[data-region-feature-count]')).toHaveText(/^\d+ features?$/);
    await expect(popup.locator('[data-region-new-drawing]')).toBeEnabled();
    expect(record).toBe(16_569);

    await popup.locator('[data-region-add-annotation]').click();
    expect(await regionRows(page)).toEqual([expect.objectContaining({
      label: 'Region 1',
      target: expect.objectContaining({ kind: 'coordinateSpan', start: 4001, end: 6001, wrapsOrigin: false })
    })]);
    expect((await undoLabels(page)).length, 'Add to annotation set is one History step').toBeGreaterThan(0);

    await generate(page);
    const [rowId] = (await regionRows(page)).map((row) => row.id);
    const drawn = await page.evaluate(() => window.__GBDRAW_APP__.currentResults().map((result) => String(result.content || '')).join(''));
    expect(drawn).toContain(`data-gbdraw-annotation-id="${rowId}"`);
    // The band is not part of the Result.
    expect(drawn).not.toContain('region-select-band');
  } finally {
    await page.context().close();
  }
});

test('Circular: a span across the origin leaves New drawing from region disabled with its reason', async ({ browser }) => {
  test.setTimeout(420_000);
  const page = await load(browser, seeds.circular);
  try {
    await generate(page);
    await turnOn(page);
    // From 16,000 through the origin (display index 0 is base 1) to 500.
    await dragAlong(page, 0, [16_000.5, 16_400.5, 16_568.5, 200.5, 500.5]);
    const popup = regionPopup(page);
    await expect(popup).toBeVisible();
    await expect(page.locator('[data-preview-item-subtitle]')).toHaveText('16,001-501 (1,070 bp)');
    await expect(popup.locator('[data-region-new-drawing]')).toBeDisabled();
    await expect(popup.locator('[data-region-new-drawing-reason]')).toHaveText('A region that crosses the origin cannot be drawn on its own.');
    await popup.locator('[data-region-add-annotation]').click();
    expect(await regionRows(page)).toEqual([expect.objectContaining({
      target: expect.objectContaining({ start: 16_001, end: 501, wrapsOrigin: true })
    })]);
    // A click elsewhere clears the selection and its band.
    await page.mouse.click(5, 5);
    await expect(page.locator('[data-region-select-band]')).toHaveCount(0);
  } finally {
    await page.context().close();
  }
});

test('Linear: a drag along a row selects a span that Add to annotation set draws on the next Generate', async ({ browser }) => {
  test.setTimeout(420_000);
  const page = await load(browser, seeds.linear);
  try {
    await generate(page);
    await turnOn(page);
    await dragAlong(page, 0, [10_000.5, 12_000.5, 15_000.5]);
    const popup = regionPopup(page);
    await expect(popup).toBeVisible();
    await expect(page.locator('[data-preview-item-subtitle]')).toHaveText('10,001-15,001 (5,001 bp)');
    await popup.locator('[data-region-add-annotation]').click();
    const [row] = await regionRows(page);
    expect(row).toEqual(expect.objectContaining({ label: 'Region 1', target: expect.objectContaining({ start: 10_001, end: 15_001 }) }));
    await generate(page);
    const drawn = await page.evaluate(() => window.__GBDRAW_APP__.currentResults().map((result) => String(result.content || '')).join(''));
    expect(drawn).toContain(`data-gbdraw-annotation-id="${row.id}"`);
  } finally {
    await page.context().close();
  }
});

test('Save Session keeps the axis of the Result, and its Load enables Select region without Generate', async ({ page }, testInfo) => {
  test.setTimeout(300_000);
  await openWithGenBank(page, HMMT);
  await generateAndWaitForResult(page);
  // The Session is the positive fixture tests/fixtures/sessions/HmmtDNA-region-axis.v51-web.gbdraw-session.json.gz.
  const saved = testInfo.outputPath('region-axis.gbdraw-session.json.gz');
  const session = JSON.parse(gunzipSync(await download(page, 'Save Session', saved, { title: '' })).toString('utf8'));
  const [record] = session.drawings[0].runMetadata.trackSlotGeometry.records;
  expect(record).toMatchObject({ recordKey: expect.any(String), axis: { kind: 'circular', source: { drawnLength: 16_569 } } });
  await loadSessionFile(page, saved);
  await expect(page.locator('[data-region-select-toggle]')).toBeEnabled();
});
