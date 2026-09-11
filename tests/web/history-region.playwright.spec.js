const { test, expect } = require('@playwright/test');
const { load, generate } = require('./helpers/mode-transition.cjs');
const { generateAndWaitForResult } = require('./helpers/app-lifecycle.cjs');

const region = page => page.evaluate(() => {
  const form = window.__GBDRAW_APP__.form;
  return [form.circular_region_start, form.circular_region_end];
});
const edit = async (page, label, value) => {
  const before = await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
  const input = page.getByLabel(label, { exact: true });
  await input.fill(value);
  await input.press('Tab');
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(before + 1);
};

test('Circular region restores number to Auto @pr-smoke', async ({ browser }) => {
  test.setTimeout(300000);
  const page = await load(browser);
  try {
    await generate(page);
    await page.locator('[data-circular-record-presentation] summary').click();
    expect(await region(page)).toEqual([null, null]);
    await edit(page, 'Circular region start', '1000');
    await page.evaluate(() => window.__GBDRAW_HISTORY__.undo());
    expect(await region(page)).toEqual([null, null]);
    await expect(page.getByLabel('Circular region start', { exact: true })).toHaveValue('');
    await page.evaluate(() => window.__GBDRAW_HISTORY__.redo());
    expect(await region(page)).toEqual([1000, null]);
    expect(page.externalRequests).toEqual([]);
  } finally { await page.context().close(); }
});

test('rejected Circular region can Undo both edits back to valid Auto @pr-smoke', async ({ browser }) => {
  test.setTimeout(300000);
  const page = await load(browser);
  try {
    await generate(page);
    const before = await page.evaluate(() => ({
      undo: window.__GBDRAW_HISTORY__.getUndoCount(), result: window.__GBDRAW_APP__.results[0].content
    }));
    await page.locator('[data-circular-record-presentation] summary').click();
    await edit(page, 'Circular region start', '1000');
    await edit(page, 'Circular region end', '500');
    await generateAndWaitForResult(page, { expectedStatus: 'error' });
    expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(before.undo + 2);
    await page.evaluate(() => window.__GBDRAW_HISTORY__.undo());
    await page.evaluate(() => window.__GBDRAW_HISTORY__.undo());
    expect(await region(page)).toEqual([null, null]);
    await expect(page.getByLabel('Circular region start', { exact: true })).toHaveValue('');
    await expect(page.getByLabel('Circular region end', { exact: true })).toHaveValue('');
    await generate(page);
    expect(await page.evaluate(() => window.__GBDRAW_APP__.results[0].content)).toBe(before.result);
    expect(page.externalRequests).toEqual([]);
  } finally { await page.context().close(); }
});
