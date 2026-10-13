// Labels that do not fit (LP-D43, LP-D45, LP-D48, OD-F14): overlapping labels
// open the dialog after a Generate that drew nothing; crossing leaders show a
// notice bar over the drawn figure, and its "Choose a fix" computes the rows
// Generate did not search. Apply is one History step and Generates; Cancel
// writes nothing.
const { readFileSync } = require('node:fs');
const { join } = require('node:path');
const { test, expect } = require('@playwright/test');
const { evaluateWithRetainedPromise, generateAndWaitForResult, openApp } = require('./helpers/app-lifecycle.cjs');

const ROW_KINDS = ['outward_both', 'outward_vertical', 'outward_horizontal', 'smaller_font', 'tighter_rows', 'combination'];
const SCREENSHOTS = process.env.GBDRAW_LABEL_FIT_SCREENSHOTS || '';

// 80 CDS with long products on 20 kb: at 18 pt one label pair overlaps.
const denseGenbank = (count, length) => {
  let seed = 1;
  const bases = Array.from({ length }, () => {
    seed = (seed * 1103515245 + 12345) % 2147483648;
    return 'acgt'[(seed >> 16) & 3];
  }).join('');
  const step = Math.floor(length / count);
  const features = Array.from({ length: count }, (_, index) => {
    const start = index * step + 11;
    const span = `${start}..${start + step - 41}`;
    return [`     CDS             ${index % 2 ? span : `complement(${span})`}`,
      `                     /product="hypothetical protein ${index + 1}"`,
      `                     /locus_tag="D_${String(index + 1).padStart(3, '0')}"`].join('\n');
  });
  const origin = [];
  for (let offset = 0; offset < length; offset += 60) {
    origin.push(`${String(offset + 1).padStart(9, ' ')} ${bases.slice(offset, offset + 60).match(/.{1,10}/g).join(' ')}`);
  }
  return [`LOCUS       DENSE1 ${String(length).padStart(21, ' ')} bp    DNA     circular UNK 01-JAN-1980`,
    'DEFINITION  Synthetic dense labels.', 'ACCESSION   DENSE1', 'VERSION     DENSE1', 'KEYWORDS    .',
    'SOURCE      .', '  ORGANISM  .', '            .', 'FEATURES             Location/Qualifiers',
    ...features, 'ORIGIN', ...origin, '//', ''].join('\n');
};

const openWithOuterLabels = async (page, { name, text, fontSize }) => {
  await page.setViewportSize({ width: 1440, height: 1000 });
  await openApp(page);
  await page.getByLabel('GenBank/DDBJ File', { exact: true })
    .setInputFiles({ name, mimeType: 'text/plain', buffer: Buffer.from(text) });
  await page.waitForFunction(() => window.__GBDRAW_APP__.circularRecordList.length === 1);
  await page.evaluate((size) => {
    const app = window.__GBDRAW_APP__;
    app.form.labels_mode = 'out';
    app.adv.label_font_size = size;
  }, fontSize);
};
const undoCount = (page) => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
const labelSettings = (page) => page.evaluate(() => {
  const { adv } = window.__GBDRAW_APP__;
  return { x: adv.outer_label_x_offset, y: adv.outer_label_y_offset, font: adv.label_font_size,
    spacing: adv.circular_label_spacing };
});
const rowStatuses = (dialog) => dialog.locator('[data-label-fit-row]').evaluateAll((rows) => (
  rows.map((row) => [row.getAttribute('data-label-fit-row'), row.getAttribute('data-label-fit-status')])
));
const shot = async (locator, name) => {
  if (SCREENSHOTS) await locator.screenshot({ path: join(SCREENSHOTS, name) });
};

test('overlapping labels: the dialog shows every row; Cancel writes nothing; Apply is one step that draws; Undo restores', async ({ page }) => {
  test.setTimeout(600000);
  await openWithOuterLabels(page, { name: 'dense.gbk', text: denseGenbank(80, 20000), fontSize: 18 });
  await generateAndWaitForResult(page, { expectedStatus: 'error', requireCommittedResult: false });
  const dialog = page.locator('[data-label-fit-dialog]');
  await expect(dialog).toBeVisible();
  await expect(dialog.getByRole('heading')).toHaveText('Labels do not fit');
  // Python searched every row before the dialog opened, in its order.
  const statuses = await rowStatuses(dialog);
  expect(statuses.map(([kind]) => kind)).toEqual(ROW_KINDS);
  expect(statuses.every(([, status]) => status === 'found' || status === 'not_found')).toBe(true);
  for (const [kind, status] of statuses) {
    await expect(dialog.locator(`[data-label-fit-row="${kind}"] input[type="radio"]`))
      .toBeEnabled({ enabled: status === 'found' });
  }
  const vertical = dialog.locator('[data-label-fit-row="outward_vertical"]');
  await expect(vertical).toHaveAttribute('data-label-fit-status', 'found');
  await expect(vertical).toContainText('Outer Y Offset');
  await expect(vertical).toContainText(/\d+ × \d+ px/);
  await shot(dialog.getByRole('dialog'), 'label-fit-dialog-overlap.png');

  const before = { undo: await undoCount(page), settings: await labelSettings(page) };
  await dialog.getByRole('button', { name: 'Cancel', exact: true }).click();
  await expect(dialog).toBeHidden();
  expect(await undoCount(page)).toBe(before.undo);
  expect(await labelSettings(page)).toEqual(before.settings);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.results.length)).toBe(0);

  await generateAndWaitForResult(page, { expectedStatus: 'error', requireCommittedResult: false });
  await expect(dialog).toBeVisible();
  await vertical.locator('input[type="radio"]').check();
  const offset = Number((await vertical.innerText()).match(/Outer Y Offset ([\d.]+)/)[1]);
  const key = await page.evaluate(async () => (await import('./js/state.js')).state.resultGenerationKey.value);
  await dialog.getByRole('button', { name: 'Apply and Generate', exact: true }).click();
  await expect(dialog).toBeHidden();
  await expect.poll(() => page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    return { key: state.resultGenerationKey.value, processing: state.processing.value, error: state.errorLog.value };
  }), { timeout: 300000 }).toEqual({ key: key + 1, processing: false, error: null });
  expect(await page.evaluate(() => window.__GBDRAW_APP__.results.length)).toBe(1);
  expect(await undoCount(page)).toBe(before.undo + 1);
  expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.undoLabel())).toBe('Fit labels: outward, vertical only');
  expect((await labelSettings(page)).y).toBe(offset);

  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
  expect(await labelSettings(page)).toEqual(before.settings);
});

test('crossing leaders: the notice bar opens the dialog, which computes the rows; Cancel keeps the figure', async ({ page }) => {
  test.setTimeout(600000);
  const text = readFileSync(join(process.cwd(), 'tests/test_inputs/SARS-CoV-2.gbk'), 'utf8');
  await openWithOuterLabels(page, { name: 'SARS-CoV-2.gbk', text, fontSize: 14 });
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
  await expect(dialog.getByRole('heading')).toHaveText('Leader lines cross');
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
