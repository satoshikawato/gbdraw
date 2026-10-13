// Labels that do not fit (LP-D43, LP-D45, LP-D48, OD-F14): overlapping labels
// open the dialog after a Generate that drew nothing; crossing leaders show a
// notice bar over the drawn figure, and its "Choose a fix" computes the rows
// Generate did not search. Apply is one History step and Generates; Cancel
// writes nothing.
const { readFileSync } = require('node:fs');
const { join } = require('node:path');
const { test, expect } = require('@playwright/test');
const { evaluateWithRetainedPromise, generateAndWaitForResult, openApp, reveal } = require('./helpers/app-lifecycle.cjs');

const ROW_KINDS = ['outward_both', 'outward_vertical', 'outward_horizontal', 'smaller_font', 'tighter_rows', 'combination'];
const SCREENSHOTS = process.env.GBDRAW_LABEL_FIT_SCREENSHOTS || '';

const SARS = () => readFileSync(join(process.cwd(), 'tests/test_inputs/SARS-CoV-2.gbk'), 'utf8');

const openWithOuterLabels = async (page, { name, text, fontSize = null, rendering = 'auto' }) => {
  await page.setViewportSize({ width: 1440, height: 1000 });
  await openApp(page);
  await page.getByLabel('GenBank/DDBJ File', { exact: true })
    .setInputFiles({ name, mimeType: 'text/plain', buffer: Buffer.from(text) });
  await page.waitForFunction(() => window.__GBDRAW_APP__.circularRecordList.length === 1);
  await page.evaluate(({ size, labelRendering }) => {
    const app = window.__GBDRAW_APP__;
    app.form.labels_mode = 'out';
    if (size !== null) app.adv.label_font_size = size;
    app.adv.label_rendering = labelRendering;
  }, { size: fontSize, labelRendering: rendering });
};
const undoCount = (page) => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
const labelSettings = (page) => page.evaluate(() => {
  const { adv } = window.__GBDRAW_APP__;
  return { x: adv.outer_label_x_offset, y: adv.outer_label_y_offset, font: adv.label_font_size,
    spacing: adv.circular_label_spacing, resolveOverlaps: adv.resolve_overlaps };
});
const rowStatuses = (dialog) => dialog.locator('[data-label-fit-row]').evaluateAll((rows) => (
  rows.map((row) => [row.getAttribute('data-label-fit-row'), row.getAttribute('data-label-fit-status')])
));
const shot = async (locator, name) => {
  if (SCREENSHOTS) await locator.screenshot({ path: join(SCREENSHOTS, name) });
};

// LP-D52 (E-ROW): SARS-CoV-2 with labels Out overprints two embedded labels;
// with Separate Strands off only "separate overlapping features" fits it. (With
// the Web's Separate Strands on, its outer leaders still cross after that fix,
// so no row fits: OV-724.)
test('overlapping labels: the dialog shows every row; Cancel writes nothing; Apply is one step that draws; Undo restores', async ({ page }) => {
  test.setTimeout(600000);
  await openWithOuterLabels(page, { name: 'SARS-CoV-2.gbk', text: SARS() });
  await reveal(page.locator('#circular-track-preset'));
  await page.getByRole('checkbox', { name: 'Separate Strands', exact: true }).uncheck();
  await generateAndWaitForResult(page, { expectedStatus: 'error', requireCommittedResult: false });
  const dialog = page.locator('[data-label-fit-dialog]');
  await expect(dialog).toBeVisible();
  await expect(dialog.getByRole('heading')).toHaveText('Labels do not fit');
  await expect(dialog.locator('[data-label-fit-records]')).toContainText('3 embedded labels at 14 pt (1 label pair overlaps');
  // Python searched every row before the dialog opened, in its order.
  const statuses = await rowStatuses(dialog);
  expect(statuses.map(([kind]) => kind)).toEqual([...ROW_KINDS, 'radial', 'separate_features']);
  expect(statuses.filter(([, status]) => status === 'found').map(([kind]) => kind)).toEqual(['separate_features']);
  expect(statuses.every(([, status]) => status === 'found' || status === 'not_found')).toBe(true);
  for (const [kind, status] of statuses) {
    await expect(dialog.locator(`[data-label-fit-row="${kind}"] input[type="radio"]`))
      .toBeEnabled({ enabled: status === 'found' });
  }
  const separate = dialog.locator('[data-label-fit-row="separate_features"]');
  await expect(separate).toContainText('separate overlapping features');
  await expect(separate).toContainText('Resolve Overlaps On');
  await expect(separate).toContainText(/\d+ × \d+ px/);
  await shot(dialog.getByRole('dialog'), 'label-fit-dialog-separate-features.png');

  const before = { undo: await undoCount(page), settings: await labelSettings(page) };
  expect(before.settings.resolveOverlaps).toBe(false);
  await dialog.getByRole('button', { name: 'Cancel', exact: true }).click();
  await expect(dialog).toBeHidden();
  expect(await undoCount(page)).toBe(before.undo);
  expect(await labelSettings(page)).toEqual(before.settings);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.results.length)).toBe(0);

  await generateAndWaitForResult(page, { expectedStatus: 'error', requireCommittedResult: false });
  await expect(dialog).toBeVisible();
  // The only row with a fix is preselected.
  await expect(separate.locator('input[type="radio"]')).toBeChecked();
  const key = await page.evaluate(async () => (await import('./js/state.js')).state.resultGenerationKey.value);
  await dialog.getByRole('button', { name: 'Apply and Generate', exact: true }).click();
  await expect(dialog).toBeHidden();
  await expect.poll(() => page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    return { key: state.resultGenerationKey.value, processing: state.processing.value, error: state.errorLog.value };
  }), { timeout: 300000 }).toEqual({ key: key + 1, processing: false, error: null });
  expect(await page.evaluate(() => window.__GBDRAW_APP__.results.length)).toBe(1);
  // It fits: no overlap, no crossing, so no notice.
  await expect(page.getByTestId('label-fit-notice')).toBeHidden();
  // The row's settings are one step (PD-OI-088); the Generate that follows is its own.
  expect(await undoCount(page)).toBe(before.undo + 2);
  expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.undoLabel())).toBe('Generate diagram');
  expect(await labelSettings(page)).toEqual({ ...before.settings, resolveOverlaps: true });
  await expect(page.getByLabel('Resolve Overlaps', { exact: true })).toBeChecked();

  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
  expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.undoLabel())).toBe('Fit labels: separate overlapping features');
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
  expect(await labelSettings(page)).toEqual(before.settings);
});

test('crossing leaders: the notice bar opens the dialog, which computes the rows; Cancel keeps the figure', async ({ page }) => {
  test.setTimeout(600000);
  // External labels only: the embedded labels would overlap (test above); outside, leaders cross.
  await openWithOuterLabels(page, { name: 'SARS-CoV-2.gbk', text: SARS(), fontSize: 14, rendering: 'external_only' });
  await generateAndWaitForResult(page);
  const notice = page.getByTestId('label-fit-notice');
  await expect(notice).toBeVisible();
  await expect(notice).toContainText(/Leader lines cross \(\d+ pairs\)\./);
  // OD-F14: Generate searched no row.
  expect(await page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    return state.labelFitReports.value[0].report.rows.map((row) => row.status);
  })).toEqual(ROW_KINDS.map(() => 'not_computed'));
  await shot(notice, 'label-fit-notice.png');

  const before = { undo: await undoCount(page), settings: await labelSettings(page),
    result: await page.evaluate(() => window.__GBDRAW_APP__.results[0].content) };
  await notice.getByRole('button', { name: 'Choose a fix', exact: true }).click();
  const dialog = page.locator('[data-label-fit-dialog]');
  await expect(dialog).toBeVisible();
  await expect(dialog.getByRole('heading')).toHaveText('Labels do not fit');
  await expect(dialog).toContainText('The figure is drawn.');
  // The rows are computed one after another and show as they arrive.
  await expect.poll(async () => (await rowStatuses(dialog))
    .every(([, status]) => status === 'found' || status === 'not_found'), { timeout: 420000 }).toBe(true);
  expect((await rowStatuses(dialog)).map(([kind]) => kind)).toEqual(ROW_KINDS);
  await shot(dialog.getByRole('dialog'), 'label-fit-dialog-crossing.png');

  await dialog.getByRole('button', { name: 'Cancel', exact: true }).click();
  await expect(dialog).toBeHidden();
  await expect(notice).toBeVisible();
  expect(await undoCount(page)).toBe(before.undo);
  expect(await labelSettings(page)).toEqual(before.settings);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.results[0].content)).toBe(before.result);
});
