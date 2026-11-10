// Select region (PD-OI-115): the preview toolbar toggle, a drag along a
// record's axis, the band and the Region popup, and its actions. Points are
// placed from the axis Python wrote (services/region-axis.js), so a drag
// lands where Python says a base is drawn.
const { test, expect } = require('@playwright/test');
const { execFileSync } = require('node:child_process');
const { readFileSync, writeFileSync } = require('node:fs');
const { gunzipSync } = require('node:zlib');
const { evaluateWithRetainedPromise, generateAndWaitForResult } = require('./helpers/app-lifecycle.cjs');
const { HMMT, loadSessionFile, openFresh, openWithGenBank } = require('./helpers/audit-browser.cjs');
const { download, generate, load, seeds } = require('./helpers/mode-transition.cjs');
const { showPhonePane } = require('./helpers/drawing-bar.cjs');
const { circularFileInput } = require('./helpers/circular-inputs.cjs');

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
 * The Region popup frame covers neither the band, nor the Legend, nor any feature label or record definition (D-28),
 * nor the preview toolbar where the tool is turned off, the export bar, or the search bar (D-30). The popup is placed
 * again once rendered, so the boxes are read when its place has settled.
 * @param {import('@playwright/test').Page} page
 */
const expectPopupClear = async (page) => {
  const boxes = () => page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    const box = (element) => {
      const { left, right, top, bottom } = element.getBoundingClientRect();
      return { left, right, top, bottom };
    };
    const svg = state.svgContainer.value.querySelector('svg');
    return {
      popup: box(document.querySelector('[data-preview-item-popup="region"]').closest('.pairwise-match-popup')),
      band: box(document.querySelector('[data-region-select-band]')),
      legend: box(svg.querySelector('[data-gbdraw-composition-role="legend"]')),
      toolbar: box(document.querySelector('.preview-controls')),
      bars: Array.from(document.querySelectorAll('.preview-toolbar, .preview-feature-search'), box).filter((bar) => bar.right > bar.left),
      labels: Array.from(svg.querySelectorAll('text[data-label-feature-id], [data-gbdraw-role^="record-definition"]'), box)
        .filter((label) => label.right > label.left)
    };
  });
  let last = '';
  await expect.poll(async () => {
    const now = JSON.stringify((await boxes()).popup);
    const settled = now === last;
    last = now;
    return settled;
  }, { intervals: [250] }).toBe(true);
  const { popup, band, legend, labels, toolbar, bars } = await boxes();
  const meets = (a, b) => a.left < b.right && b.left < a.right && a.top < b.bottom && b.top < a.bottom;
  expect(labels.length, 'the Result draws feature labels').toBeGreaterThan(0);
  expect(meets(popup, band), 'popup over the band').toBe(false);
  expect(meets(popup, legend), 'popup over the Legend').toBe(false);
  expect(labels.filter((label) => meets(popup, label)), 'labels under the popup').toEqual([]);
  expect(toolbar.right > toolbar.left, 'the preview toolbar is shown').toBe(true);
  expect(meets(popup, toolbar), 'popup over the preview toolbar').toBe(false);
  expect(bars.filter((bar) => meets(popup, bar)), 'export or search bar under the popup').toEqual([]);
};

/**
 * Fits the figure to the preview. The fit animates the preview transform, so screen points are read once it has
 * settled.
 * @param {import('@playwright/test').Page} page
 */
const fitToPreview = async (page) => {
  await page.getByRole('button', { name: 'Fit to preview', exact: true }).click();
  let last = '';
  await expect.poll(async () => {
    const now = await page.evaluate(() => JSON.stringify(document.querySelector('.gbdraw-preview-surface').getBoundingClientRect()));
    const settled = now === last;
    last = now;
    return settled;
  }, { intervals: [250] }).toBe(true);
};

/**
 * Toggles Select region on after a Generate; Layout edit turns off. The figure is fitted to the preview first, so
 * every point of the drags is on the screen.
 * @param {import('@playwright/test').Page} page
 */
const turnOn = async (page) => {
  await fitToPreview(page);
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

/**
 * The Region popup's heading: "Region", the record, the span, and the feature count.
 * @param {import('@playwright/test').Page} page
 */
const expectHeading = async (page, record, span) => {
  await expect(page.locator('[data-match-popup-title]')).toHaveText('Region');
  await expect(page.locator('[data-preview-item-subtitle]')).toHaveText(record);
  await expect(page.locator('[data-region-span]')).toHaveText(span);
};

/**
 * The Gallery Session as a released version (0.14.x) wrote it: its Result's track-slot geometry records carry no
 * `axis` and no `recordKey` (#1620 added both). Such a Session still loads its Result as saved until its next
 * Generate (D-36); the Gallery's own Sessions now carry the axis, so the test removes it.
 * @param {import('@playwright/test').TestInfo} testInfo
 */
const sessionWithoutAxis = (testInfo) => {
  const session = JSON.parse(readFileSync(seeds.circular, 'utf8'));
  const records = session.drawings[0].runMetadata.trackSlotGeometry.records;
  expect(records.length).toBeGreaterThan(0);
  for (const record of records) {
    expect(record.axis, 'the Gallery Session carries the axis').toBeDefined();
    delete record.axis;
    delete record.recordKey;
  }
  const path = testInfo.outputPath('HmmtDNA_basic_circular-no-axis.gbdraw-session.json');
  writeFileSync(path, JSON.stringify(session));
  return path;
};

test('Circular: a loaded Result without axis disables Select region; after Generate a drag selects a span that Add to annotation set draws', async ({ browser }, testInfo) => {
  test.setTimeout(420_000);
  const page = await load(browser, sessionWithoutAxis(testInfo));
  try {
    const toggle = page.locator('[data-region-select-toggle]');
    // The Result was saved without the axis: no Python starts at Load.
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
    const stepsBeforeDrag = (await undoLabels(page)).length;
    await dragAlong(page, 0, [4000.5, 5000.5, 6000.5]);
    const popup = regionPopup(page);
    await expect(popup).toBeVisible();
    await expect(page.locator('[data-region-select-band]')).toHaveCount(1);
    await expectHeading(page, 'NC_012920.1', '4,001-6,001 (2,001 bp)');
    await expectPopupClear(page);
    await expect(page.locator('[data-region-feature-count]')).toHaveText(/^\d+ features?$/);
    // The actions, one button each, in the approved order (PD-OI-115).
    expect(await popup.locator('button').evaluateAll((buttons) => buttons.map((button) => button.textContent.trim())))
      .toEqual(['Copy FASTA', expect.stringMatching(/^Select \d+ features?$/), 'Add to annotation set', 'New drawing from region', 'Use as search query']);
    await expect(popup.locator('[data-region-new-drawing]')).toBeEnabled();
    expect(record).toBe(16_569);
    // The selection is transient: the drag and the popup add no History step (Q-A4).
    expect((await undoLabels(page)).length, 'the drag is no History step').toBe(stepsBeforeDrag);

    await popup.locator('[data-region-add-annotation]').click();
    expect(await regionRows(page)).toEqual([expect.objectContaining({
      label: 'Region 1',
      target: expect.objectContaining({ kind: 'coordinateSpan', start: 4001, end: 6001, wrapsOrigin: false })
    })]);
    // The click's step commits after its handlers (app/history-inputs.js).
    await expect.poll(async () => (await undoLabels(page)).length, { message: 'Add to annotation set is one History step' })
      .toBe(stepsBeforeDrag + 1);
    await page.waitForTimeout(500);
    expect((await undoLabels(page)).length, 'and no second step follows').toBe(stepsBeforeDrag + 1);

    await generate(page);
    const [rowId] = (await regionRows(page)).map((row) => row.id);
    const drawn = await page.evaluate(() => window.__GBDRAW_APP__.currentResults().map((result) => String(result.content || '')).join(''));
    expect(drawn).toContain(`data-gbdraw-annotation-id="${rowId}"`);
    // The band is not part of the Result.
    expect(drawn).not.toContain('region-select-band');

    // OV-2482: the drawn region annotation (a group of the record) leaves the record's frame found: a drag still
    // selects.
    await page.keyboard.press('Escape');
    await fitToPreview(page);
    if (await toggle.getAttribute('aria-pressed') !== 'true') await toggle.click();
    await dragAlong(page, 0, [8000.5, 8500.5, 9000.5]);
    await expect(page.locator('[data-region-span]')).toHaveText('8,001-9,001 (1,001 bp)');
  } finally {
    await page.context().close();
  }
});

test('Circular: a span across the origin says so and leaves New drawing from region enabled', async ({ browser }) => {
  test.setTimeout(420_000);
  const page = await load(browser, seeds.circular);
  try {
    await generate(page);
    await turnOn(page);
    // From 16,000 through the origin (display index 0 is base 1) to 500.
    await dragAlong(page, 0, [16_000.5, 16_400.5, 16_568.5, 200.5, 500.5]);
    const popup = regionPopup(page);
    await expect(popup).toBeVisible();
    await expectHeading(page, 'NC_012920.1', '16,001-501 (1,070 bp)');
    await expect(page.locator('[data-region-feature-count]')).toHaveText(/^\d+ features? · crosses the origin$/);
    // Python draws such a span as one window across the origin.
    await expect(popup.locator('[data-region-new-drawing]')).toBeEnabled();
    await expectPopupClear(page);

    // Copy FASTA: the popup started no Python; the click asks Python for the bases on the source strand.
    await page.context().grantPermissions(['clipboard-read', 'clipboard-write'], { origin: new URL(page.url()).origin });
    await page.evaluate(() => {
      const hooks = window.__GBDRAW_TEST_HOOKS__;
      const metric = hooks.onStructuralMetric;
      window.__SPAN_READS__ = 0;
      hooks.onStructuralMetric = (event) => {
        metric?.(event);
        if (event.name === 'diagramHelperRequest' && event.operation === 'readRecordSpanSequences') window.__SPAN_READS__ += 1;
      };
    });
    await popup.locator('[data-region-copy-fasta]').click();
    await expect(popup.locator('[data-region-copy-status]')).toHaveText('Copied FASTA');
    expect(await page.evaluate(() => window.__SPAN_READS__)).toBe(1);
    const fasta = await page.evaluate(() => navigator.clipboard.readText());
    const [header, ...lines] = fasta.trim().split('\n');
    expect(header).toBe('>NC_012920.1:16001-501');
    expect(lines.join('')).toMatch(/^[ACGTN]{1070}$/);
    expect(lines.slice(0, -1).every((line) => line.length === 60)).toBe(true);
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

// PD-OI-115 and PD-OI-117: Use as search query adds the span's bases on the source strand, read by Python, as a query
// named "<record>:<start>-<end>", and opens the Sequence search section; the popup closes, and no History step is made.
test('Circular: Use as search query adds a span across the origin as a query and opens Sequence search', async ({ browser }) => {
  test.setTimeout(420_000);
  const page = await load(browser, seeds.circular);
  try {
    await generate(page);
    await turnOn(page);
    await dragAlong(page, 0, [16_000.5, 16_400.5, 16_568.5, 200.5, 500.5]);
    const popup = regionPopup(page);
    await expect(popup).toBeVisible();
    const steps = (await undoLabels(page)).length;
    await popup.locator('[data-region-search-query]').click();
    await expect(popup).toHaveCount(0);
    const queries = () => page.evaluate(() => window.__GBDRAW_APP__.sequenceSearch.panel.queries.map(({ name, length }) => ({ name, length })));
    await expect.poll(queries).toEqual([{ name: 'NC_012920.1:16001-501', length: 1070 }]);
    const section = page.locator('[data-settings-section="sequence-search"]');
    await expect(section.locator('details')).toHaveAttribute('open', '');
    await expect(section.locator('[data-sequence-search-query]')).toHaveText([/NC_012920\.1:16001-501/]);
    await expect.poll(() => page.evaluate(() => document.activeElement?.closest('[data-settings-section]')?.dataset.settingsSection ?? null))
      .toBe('sequence-search');
    expect((await undoLabels(page)).length, 'no History step').toBe(steps);
    // The same span again adds nothing twice.
    await dragAlong(page, 0, [16_000.5, 16_400.5, 16_568.5, 200.5, 500.5]);
    await regionPopup(page).locator('[data-region-search-query]').click();
    await expect(regionPopup(page)).toHaveCount(0);
    expect(await queries()).toEqual([{ name: 'NC_012920.1:16001-501', length: 1070 }]);
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
    await expectHeading(page, 'NC_001416.1', '10,001-15,001 (5,001 bp)');
    await expectPopupClear(page);
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

test('Circular: on a record that is not circular a drag through 12 o\'clock stops at its last base', async ({ page }) => {
  test.setTimeout(300_000);
  // Lambda's LOCUS says linear: drawn in Circular mode, its ends meet at 12 o'clock but the record has no origin to cross.
  await openWithGenBank(page, 'tests/test_inputs/NC_001416.gb');
  await generateAndWaitForResult(page);
  await turnOn(page);
  await dragAlong(page, 0, [45_500.5, 47_000.5, 48_490.5, 1_000.5, 3_000.5]);
  const popup = regionPopup(page);
  await expect(popup).toBeVisible();
  await expectHeading(page, 'NC_001416.1', '45,501-48,502 (3,002 bp)');
  await expect(page.locator('[data-region-feature-count]')).not.toContainText('crosses the origin');
  await expect(popup.locator('[data-region-new-drawing]')).toBeEnabled();
  await popup.locator('[data-region-add-annotation]').click();
  const [row] = await regionRows(page);
  expect(row.target).toEqual(expect.objectContaining({ start: 45_501, end: 48_502, wrapsOrigin: false }));
  await generateAndWaitForResult(page);
  const drawn = await page.evaluate(() => window.__GBDRAW_APP__.currentResults().map((result) => String(result.content || '')).join(''));
  expect(drawn).toContain(`data-gbdraw-annotation-id="${row.id}"`);
});

test('Circular at phone width: the Region popup stays above the preview toolbar', async ({ browser }) => {
  test.setTimeout(300_000);
  const page = await load(browser, seeds.circular, { width: 390, height: 844 });
  try {
    await generate(page);
    await showPhonePane(page, 'preview');
    await turnOn(page);
    await dragAlong(page, 0, [4000.5, 5000.5, 6000.5]);
    await expect(regionPopup(page)).toBeVisible();
    // The popup is placed again once rendered: read its box when it has settled.
    let last = '';
    await expect.poll(async () => {
      const now = await page.evaluate(() => JSON.stringify(document.querySelector('.pairwise-match-popup').getBoundingClientRect()));
      const settled = now === last;
      last = now;
      return settled;
    }, { intervals: [250] }).toBe(true);
    const [popup, toolbar] = await page.evaluate(() => ['.pairwise-match-popup', '.preview-controls']
      .map((selector) => document.querySelector(selector).getBoundingClientRect().toJSON()));
    expect(toolbar.top, 'the toolbar is in the lower half').toBeGreaterThan(422);
    expect(popup.bottom, 'popup over the preview toolbar').toBeLessThanOrEqual(toolbar.top);
    await expect(page.locator('[data-region-select-toggle]')).toBeVisible();
  } finally {
    await page.context().close();
  }
});

test.describe('touch', () => {
  test.use({ hasTouch: true });

  test('Circular: a touch drag on the tracks selects a span, one elsewhere pans; the canvas takes no touch scroll while the tool is on', async ({ page }) => {
    test.setTimeout(300_000);
    await openWithGenBank(page, HMMT);
    await generateAndWaitForResult(page);
    const canvas = page.locator('.preview-canvas');
    const touchAction = () => canvas.evaluate((element) => getComputedStyle(element).touchAction);
    expect(await touchAction()).not.toBe('none');
    const toggle = await turnOn(page);
    expect(await touchAction()).toBe('none');

    const points = [];
    for (const position of [4_000.5, 4_500.5, 5_000.5, 5_500.5, 6_000.5]) points.push(await axisPoint(page, 0, position));
    const client = await page.context().newCDPSession(page);
    const touch = (type, { x, y }) => client.send('Input.dispatchTouchEvent', { type, touchPoints: type === 'touchEnd' ? [] : [{ x, y }] });
    await touch('touchStart', points[0]);
    for (const point of points.slice(1)) await touch('touchMove', point);
    await expect(page.locator('[data-region-select-readout]')).toBeVisible();
    await touch('touchEnd', points.at(-1));
    await expect(regionPopup(page)).toBeVisible();
    await expect(page.locator('[data-region-span]')).toHaveText('4,001-6,001 (2,001 bp)');

    await page.keyboard.press('Escape');
    // A touch outside the tracks pans the preview through the pan owner, as a mouse drag does. The point lies inside
    // the circle, halfway to the feature band, clear of the centre's record definition (a control, which no drag pans).
    const centre = await page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      const { readRegionAxisRecords, regionAxisFrame } = await import('./js/services/region-axis.js');
      const { axis } = readRegionAxisRecords(state.trackSlotResolvedGeometry.value, state.selectedResultIndex.value)[0];
      const frame = regionAxisFrame(state.svgContainer.value.querySelector('svg'), 0);
      const point = new DOMPoint(0, -axis.bandInnerPx / 2).matrixTransform(frame.getScreenCTM());
      const target = document.elementFromPoint(point.x, point.y);
      return { x: point.x, y: point.y, control: Boolean(target?.closest('button, input, textarea, select, a, [role="button"]')) };
    });
    expect(centre.control).toBe(false);
    const surface = () => page.evaluate(() => {
      const { left, top } = document.querySelector('.gbdraw-preview-surface').getBoundingClientRect();
      return { left, top };
    });
    const before = await surface();
    await touch('touchStart', centre);
    for (const step of [1, 2, 3]) await touch('touchMove', { x: centre.x + 20 * step, y: centre.y + 10 * step });
    await touch('touchEnd', { x: centre.x + 60, y: centre.y + 30 });
    await expect.poll(async () => {
      const now = await surface();
      return [Math.round(now.left - before.left), Math.round(now.top - before.top)];
    }).toEqual([60, 30]);
    await expect(regionPopup(page)).toHaveCount(0);
    await toggle.click();
    await expect(toggle).toHaveAttribute('aria-pressed', 'false');
    expect(await touchAction(), 'the canvas scrolls again once the tool is off').not.toBe('none');
  });
});

test('Linear Wrap Rows: a drag over a drawn feature on the second row selects that feature\'s span', async ({ browser }) => {
  test.setTimeout(420_000);
  const page = await load(browser, seeds.linear);
  try {
    const wrap = page.locator('#linear-wrap-rows');
    await wrap.scrollIntoViewIfNeeded();
    await wrap.check();
    await generate(page);
    await turnOn(page);
    // A feature drawn on row 2, found by the rows Python wrote and its drawn shape.
    const target = await page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      const { readRegionAxisRecords } = await import('./js/services/region-axis.js');
      const [record] = readRegionAxisRecords(state.trackSlotResolvedGeometry.value, state.selectedResultIndex.value);
      const row = record.axis.rows[1];
      const svg = state.svgContainer.value.querySelector('svg');
      const feature = window.__GBDRAW_APP__.extractedFeatures.find((entry) => {
        const start = Number(entry.start);
        const end = Number(entry.end);
        return entry.type === 'CDS' && start > row.displayStart && end < row.displayEnd && end - start > 900
          && svg.querySelector(`[id="${entry.svg_id}"]`);
      });
      const box = svg.querySelector(`[id="${feature.svg_id}"]`).getBoundingClientRect();
      return {
        start: Number(feature.start) + 1,
        end: Number(feature.end),
        rows: record.axis.rows.length,
        left: box.left,
        right: box.right,
        y: box.top + box.height / 2,
        bpPerPx: (row.displayEnd - row.displayStart) / (row.xEndPx - row.xStartPx)
          / svg.querySelector('g[data-gbdraw-record-index]').getScreenCTM().a
      };
    });
    expect(target.rows).toBeGreaterThan(1);
    await page.mouse.move(target.left + 0.5, target.y);
    await page.mouse.down();
    await page.mouse.move((target.left + target.right) / 2, target.y, { steps: 4 });
    await page.mouse.move(target.right - 0.5, target.y, { steps: 4 });
    await expect(page.locator('[data-region-select-readout]')).toBeVisible();
    await page.mouse.up();
    await expect(regionPopup(page)).toBeVisible();
    const [, start, end] = (await page.locator('[data-region-span]').textContent()).replaceAll(',', '').match(/^(\d+)-(\d+) /);
    // The drag starts and ends inside the shape: each end lies within a pixel's worth of bases of the feature's.
    const tolerance = Math.ceil(2 * target.bpPerPx) + 1;
    expect(Math.abs(Number(start) - target.start), `start ${start} vs ${target.start}`).toBeLessThanOrEqual(tolerance);
    expect(Math.abs(Number(end) - target.end), `end ${end} vs ${target.end}`).toBeLessThanOrEqual(tolerance);
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

// OV-3240: a Session the gbdraw command saved carries the axis as a Web Save does, so a drag selects on its loaded
// Result before any Generate.
test('a Session saved by the gbdraw command enables Select region at Load: a drag selects a span without Generate', async ({ browser }, testInfo) => {
  test.setTimeout(300_000);
  const saved = testInfo.outputPath('cli.gbdraw-session.json');
  execFileSync(process.env.PYTHON || 'python', [
    '-m', 'gbdraw.cli', 'circular', '--gbk', HMMT, '-o', testInfo.outputPath('cli'), '-f', 'svg', '--session_output', saved
  ], { cwd: process.cwd(), env: { ...process.env, PYTHONPATH: process.cwd() }, stdio: 'pipe' });
  const page = await load(browser, saved);
  try {
    await turnOn(page);
    await dragAlong(page, 0, [4000.5, 5000.5, 6000.5]);
    await expectHeading(page, 'NC_012920.1', '4,001-6,001 (2,001 bp)');
  } finally {
    await page.context().close();
  }
});

/** A 10 kb circular record whose four CDS lie in its first 1,200 bases, on both strands. */
const arcGenbank = (id = 'ARC1') => {
  const sequence = 'ACGTTGCA'.repeat(1250);
  const origin = sequence.match(/.{1,60}/g).map((chunk, index) => (
    `${String(index * 60 + 1).padStart(9)} ${chunk.match(/.{1,10}/g).join(' ')}`
  )).join('\n');
  const cds = [['f1', '101..300'], ['f2', 'complement(401..600)'], ['f3', '701..900'], ['f4', 'complement(1001..1200)']]
    .map(([gene, location]) => `     CDS             ${location}\n                     /gene="${gene}"\n                     /product="${gene} protein"`).join('\n');
  return `LOCUS       ${id}                   10000 bp    DNA     circular UNA 01-JAN-2000
DEFINITION  Four CDS on one arc.
ACCESSION   ${id}
VERSION     ${id}
KEYWORDS    .
SOURCE      synthetic construct
  ORGANISM  synthetic construct
            .
FEATURES             Location/Qualifiers
${cds}
ORIGIN
${origin}
//
`;
};

test('Circular: Shift+click selects the features between two features along the record, also when they cover one arc', async ({ page }) => {
  test.setTimeout(300_000);
  await openFresh(page);
  await circularFileInput(page, 'gb').setInputFiles({ name: 'ARC1.gbk', mimeType: 'text/plain', buffer: Buffer.from(arcGenbank()) });
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.circularRecordList.length)).toBe(1);
  await generateAndWaitForResult(page);
  await fitToPreview(page);
  // A feature by its record start (0-based): f1 100, f2 400, f3 700, f4 1,000.
  const centre = (start) => page.evaluate((from) => {
    const feature = window.__GBDRAW_APP__.extractedFeatures.find((entry) => Number(entry.start) === from);
    const element = document.querySelector('.preview-canvas svg').querySelector(`[id="${feature.svg_id}"]`);
    const box = element.getBoundingClientRect();
    return { id: feature.svg_id, x: box.left + box.width / 2, y: box.top + box.height / 2 };
  }, start);
  const [f1, f2, f3, f4] = [await centre(100), await centre(400), await centre(700), await centre(1000)];
  await page.keyboard.down('Control');
  await page.mouse.click(f1.x, f1.y);
  await page.keyboard.up('Control');
  await page.keyboard.down('Shift');
  await page.mouse.click(f3.x, f3.y);
  await page.keyboard.up('Shift');
  const selected = await page.evaluate(() => [...window.__GBDRAW_APP__.selectedFeatureIds].sort());
  expect(selected).toEqual([f1.id, f2.id, f3.id].sort());
  expect(selected).not.toContain(f4.id);
});

test('Add to annotation set right after a Load reads the records first and binds the row to its record', async ({ page }, testInfo) => {
  test.setTimeout(420_000);
  // Two records on one canvas: a row needs its record to be drawn. A Load starts no Python, so the record catalog
  // is not read when the Region popup opens (OV-373).
  const input = testInfo.outputPath('two-arcs.gbk');
  writeFileSync(input, arcGenbank('ARC1') + arcGenbank('ARC2'));
  await openWithGenBank(page, input, () => { window.__GBDRAW_APP__.form.multi_record_canvas = true; });
  await generateAndWaitForResult(page);
  const saved = testInfo.outputPath('two-arcs.gbdraw-session.json.gz');
  await download(page, 'Save Session', saved, { title: '' });
  await loadSessionFile(page, saved);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.circularRecordDiscoveryState.status)).not.toBe('ready');
  await turnOn(page);
  await dragAlong(page, 1, [100.5, 600.5, 1_200.5]);
  const popup = regionPopup(page);
  await expectHeading(page, 'ARC2', '101-1,201 (1,101 bp)');
  await popup.locator('[data-region-add-annotation]').click();
  await expect.poll(() => regionRows(page)).toHaveLength(1);
  await expect(popup.locator('[data-region-add-status]')).toHaveCount(0);
  const [row] = await regionRows(page);
  expect(JSON.stringify(row.target.record)).toContain('ARC2');
  await generateAndWaitForResult(page);
  const drawn = await page.evaluate(() => window.__GBDRAW_APP__.currentResults().map((result) => String(result.content || '')).join(''));
  expect(drawn).toContain(`data-gbdraw-annotation-id="${row.id}"`);
});

test('Linear: Add to annotation set whose record read fails adds nothing, takes no History step, and says what the record list says', async ({ browser }, testInfo) => {
  test.setTimeout(420_000);
  const page = await load(browser, seeds.linear);
  try {
    await generate(page);
    // The Gallery Session has a title, so Save asks for none.
    const saved = testInfo.outputPath('lambda-axis.gbdraw-session.json.gz');
    await download(page, 'Save Session', saved);
    await loadSessionFile(page, saved);
    // A Load reads no records (OV-373). From now on the card's file cannot be read, as a file gone from the disk.
    // The card's record select, as the template reads it: disabled until its records are read, and its error.
    const card = () => page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      const [seq] = state.activeDrawing().linearSeqs;
      const app = window.__GBDRAW_APP__;
      return { disabled: app.linearRecordSelectorDisabled(seq), error: app.linearRecordSelectorError(seq) };
    });
    expect(await card()).toEqual({ disabled: true, error: '' });
    await page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      const unreadable = new File(['LOCUS'], 'gone.gb');
      const refuse = () => Promise.reject(new DOMException('The file could not be read.', 'NotReadableError'));
      Object.defineProperties(unreadable, { arrayBuffer: { value: refuse }, text: { value: refuse } });
      state.inputFiles.fileOf = (ref) => (ref ? unreadable : null);
    });
    await turnOn(page);
    await dragAlong(page, 0, [10_000.5, 12_000.5, 15_000.5]);
    const popup = regionPopup(page);
    const before = await undoLabels(page);
    await popup.locator('[data-region-add-annotation]').click();
    const status = popup.locator('[data-region-add-status]');
    await expect(status).toBeVisible();
    const { error: listed } = await card();
    // A read the browser refuses is the input-unreadable failure, not the unknown one (OV-2563).
    expect(listed).toContain('An input could not be read. Replace or reselect it.');
    await expect(status).toHaveText(listed);
    await expect(popup.locator('[data-region-add-annotation]')).toBeEnabled();
    expect(await regionRows(page)).toEqual([]);
    expect(await undoLabels(page)).toEqual(before);
  } finally {
    await page.context().close();
  }
});

// OV-2632, D-58: the lambda card and an HmmtDNA card; after Generate the HmmtDNA card is deleted, so one record is
// left in the draft. A Region added on the HmmtDNA record of the Result is held as "target record unavailable", as a
// search hit is (D-57): Generate refuses it instead of drawing it on lambda.
test('Linear: Add to annotation set on a record whose card was deleted after Generate is refused at Generate, not moved', async ({ browser }) => {
  test.setTimeout(420_000);
  const page = await load(browser, seeds.linear);
  try {
    await page.evaluate(async (text) => {
      const app = window.__GBDRAW_APP__;
      app.addLinearSeq();
      app.setLinearSeqPrimaryFile(1, 'gb', new File([text], 'HmmtDNA.gbk', { type: 'text/plain', lastModified: 1000 }));
      await window.Vue.nextTick();
    }, readFileSync(HMMT, 'utf8'));
    await generate(page);
    await page.evaluate(() => { const app = window.__GBDRAW_APP__; app.requestLinearSourceRemoval(app.linearSourceGroups[1]); });
    await page.locator('[data-linear-source-removal-dialog]').getByRole('button', { name: 'Remove file', exact: true }).click();
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.linearSeqs.length)).toBe(1);
    await turnOn(page);
    await dragAlong(page, 1, [1_000.5, 2_000.5, 3_000.5]);
    const popup = regionPopup(page);
    await expect(popup).toBeVisible();
    await popup.locator('[data-region-add-annotation]').click();
    await expect.poll(() => regionRows(page)).toHaveLength(1);
    const rejected = await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.runAnalysis());
    expect(rejected.status).toBe('error');
    expect(rejected.error).toMatchObject({ code: 'ANNOTATION_TARGET', context: { reason: 'TARGET_RECORD' } });
    await expect(page.getByRole('alert', { name: 'Generation Error' })).toContainText('Choose an available target record.');
  } finally {
    await page.context().close();
  }
});
