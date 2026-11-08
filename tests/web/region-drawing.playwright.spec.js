// New drawing from a region (PD-OI-090 (4)): the drawing menu's "New drawing
// from a region…" opens the dialog of app/region-drawing.js, which derives the
// new drawing from the shown one in Python and adds it as one project step.
// The cases find the new drawing by its drawing ID, never by its tab text.
const { test, expect } = require('@playwright/test');
const { drawingIds, drawingTab, shownDrawing } = require('./helpers/drawing-bar.cjs');
const { generate, load, seeds } = require('./helpers/mode-transition.cjs');
const { chooseCoordinates, openRegionDialog } = require('./helpers/region-drawing.cjs');

test.describe.configure({ retries: 0 });

/** @param {import('@playwright/test').Page} page */
const regionStatus = (page) => page.evaluate(() => window.__GBDRAW_APP__.regionDrawings.dialog.status);

/** @param {import('@playwright/test').Page} page */
const undoLabels = (page) => page.evaluate(() => window.__GBDRAW_HISTORY__.timeline().undo.map((step) => step.label));

test('a region of a Circular record becomes a Linear drawing that Generate draws at the region length and Undo removes after the figure', async ({ browser }) => {
  test.setTimeout(420_000);
  const page = await load(browser, seeds.circular);
  try {
    await generate(page);
    const [source] = await drawingIds(page);
    const start = 1001;
    // 1,000 bp: the Linear Auto scale interval (100 bp) differs from the Circular one of the 16,569 bp record (1,000 bp).
    const end = 2000;

    // Without a size set by hand, a shorter region keeps every setting.
    let dialog = await openRegionDialog(page);
    await chooseCoordinates(dialog, start, end);
    await expect.poll(() => regionStatus(page), { timeout: 60_000 }).toBe('derived');
    const returned = dialog.locator('[data-region-drawing-returned] li');
    const before = await returned.count();
    await dialog.locator('[data-region-drawing-cancel]').click();
    await expect(dialog).toHaveCount(0);
    expect(await drawingIds(page), 'Cancel adds no drawing').toEqual([source]);

    // A scale interval set by hand returns to Auto for the shorter region.
    await page.evaluate(() => { window.__GBDRAW_APP__.adv.scale_interval = 777; });
    dialog = await openRegionDialog(page);
    await chooseCoordinates(dialog, start, end);
    await expect.poll(() => regionStatus(page), { timeout: 60_000 }).toBe('derived');
    await expect(dialog.locator('[data-region-drawing-returned]')).toBeVisible();
    await expect(returned).toHaveCount(before + 1);
    await expect(returned.filter({ hasText: '777' })).toHaveCount(1);
    const notCarried = dialog.locator('[data-region-drawing-not-carried] li');
    await expect(notCarried.or(dialog.getByText('Every setting is carried over.', { exact: true })).first()).toBeVisible();

    const create = dialog.locator('[data-region-drawing-create]');
    await expect(create).toBeEnabled();
    await create.click();
    await expect(dialog).toHaveCount(0);
    await expect.poll(async () => (await drawingIds(page)).length).toBe(2);
    const [added] = (await drawingIds(page)).filter((id) => id !== source);
    await expect(drawingTab(page, added)).toBeVisible();
    await expect(drawingTab(page, added)).toHaveAttribute('aria-selected', 'true');
    expect(await shownDrawing(page)).toEqual({ id: added, kind: 'linear' });
    expect((await undoLabels(page))[0], 'Create is one project step').toBe('New drawing from a region');

    await generate(page);
    const drawn = await page.evaluate(() => window.__GBDRAW_APP__.currentResults().map((result) => String(result.content || '')));
    expect(drawn).toHaveLength(1);
    // The record group carries the drawn span of its source record.
    expect(drawn[0]).toContain(`data-gbdraw-record-source-start="${start}"`);
    expect(drawn[0]).toContain(`data-gbdraw-record-source-end="${end}"`);

    // Generate is a step of its own: the first Undo takes the figure back, the second the new drawing.
    const undo = page.getByRole('button', { name: 'Undo', exact: true });
    await undo.click();
    await expect.poll(async () => (await undoLabels(page))[0]).toBe('New drawing from a region');
    expect(await drawingIds(page), 'Undo of Generate keeps the new drawing').toEqual([source, added]);
    await undo.click();
    await expect.poll(() => drawingIds(page)).toEqual([source]);
    await expect(drawingTab(page, added)).toHaveCount(0);
  } finally {
    await page.context().close();
  }
});

// OV-2489: after a Load the Circular record list stays deferred until Generate; the dialog reads it.
test('OV-2489: after a Session Load without Generate the dialog offers the drawing\'s records', async ({ browser }) => {
  test.setTimeout(300_000);
  const page = await load(browser, seeds.circular);
  try {
    const dialog = await openRegionDialog(page);
    await expect(dialog.locator('[data-region-drawing-row]')).toHaveCount(1);
  } finally {
    await page.context().close();
  }
});

test('Coordinates that start after their end are refused on a linear record with the reason, Create stays disabled, and no drawing is added', async ({ browser }) => {
  test.setTimeout(300_000);
  // Lambda (NC_001416) is a linear record: a start after the end names no span.
  const page = await load(browser, seeds.linear);
  try {
    const ids = await drawingIds(page);
    const steps = await undoLabels(page);
    const dialog = await openRegionDialog(page);
    await chooseCoordinates(dialog, 5000, 100);
    await expect.poll(() => regionStatus(page), { timeout: 60_000 }).toBe('refused');
    const message = dialog.locator('[data-region-drawing-message]');
    await expect(message).toHaveText('Choose a start before the end.');
    await expect(dialog.locator('[data-region-drawing-create]')).toBeDisabled();
    await dialog.locator('[data-region-drawing-cancel]').click();
    await expect(dialog).toHaveCount(0);
    expect(await drawingIds(page)).toEqual(ids);
    expect(await undoLabels(page)).toEqual(steps);
  } finally {
    await page.context().close();
  }
});

test('Coordinates across the origin of a circular record make one window across it, which Generate draws', async ({ browser }) => {
  test.setTimeout(420_000);
  const page = await load(browser, seeds.circular);
  try {
    const [source] = await drawingIds(page);
    const dialog = await openRegionDialog(page);
    await chooseCoordinates(dialog, 16000, 100);
    await expect.poll(() => regionStatus(page), { timeout: 60_000 }).toBe('derived');
    await dialog.locator('[data-region-drawing-create]').click();
    await expect(dialog).toHaveCount(0);
    const [added] = (await drawingIds(page)).filter((id) => id !== source);
    expect(await shownDrawing(page)).toEqual({ id: added, kind: 'linear' });
    // The card holds the window in its read-only list (its Start/End fields would read 16,000-100 as reversed), with
    // the bases it draws across the origin of the 16,569 bp record.
    await expect(page.locator('[data-record-window]')).toHaveText([/^1\. 16,000\u2013100\s*670 bp$/]);
    await generate(page);
    const drawn = await page.evaluate(() => window.__GBDRAW_APP__.currentResults().map((result) => String(result.content || '')));
    expect(drawn).toHaveLength(1);
    // The record group carries the drawn span: it ends before it starts, across the origin.
    expect(drawn[0]).toContain('data-gbdraw-record-source-start="16000"');
    expect(drawn[0]).toContain('data-gbdraw-record-source-end="100"');
  } finally {
    await page.context().close();
  }
});
