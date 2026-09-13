const { test, expect } = require('@playwright/test');
const { load, generate } = require('./helpers/mode-transition.cjs');

test('comparison off has its own Undo step after a record definition edit', async ({ browser }) => {
  test.setTimeout(300000);
  const page = await load(browser, 'tests/test_inputs/BGC0000708-BGC0000713.gbdraw-session.json');
  try {
    await generate(page);
    const baseline = await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
    await page.getByRole('button', { name: 'Record options for sequence 2', exact: true }).click();
    const definition = page.getByLabel('Definition for sequence 2', { exact: true });
    await definition.fill('HISTORY_RECORD_TWO');
    await definition.press('Tab');
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(baseline + 1);
    const before = await page.evaluate(async () => ({
      plan: (await import('./js/services/config.js')).buildConfigData().linearComparisonPlan,
      result: window.__GBDRAW_APP__.results[0].content
    }));
    await page.getByRole('button', { name: 'Set no comparison', exact: true }).click();
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.linearComparisonPlan.mode)).toBe('none');
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(baseline + 2);
    await page.getByRole('button', { name: 'Undo', exact: true }).click();
    await expect.poll(() => page.evaluate(async () =>
      (await import('./js/services/config.js')).buildConfigData().linearComparisonPlan)).toEqual(before.plan);
    await expect(definition).toHaveValue('HISTORY_RECORD_TWO');
    expect(await page.evaluate(() => window.__GBDRAW_APP__.results[0].content)).toBe(before.result);
    await page.getByRole('button', { name: 'Redo', exact: true }).click();
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.linearComparisonPlan.mode)).toBe('none');
    await expect(definition).toHaveValue('HISTORY_RECORD_TWO');
    await page.getByRole('button', { name: 'Set no comparison', exact: true }).click();
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(baseline + 2);
    const losat = page.getByRole('button', { name: 'Run LOSAT for all adjacent pairs', exact: true });
    await losat.focus();
    await losat.press('Enter');
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(baseline + 3);
    await page.evaluate(() => window.__GBDRAW_HISTORY__.undo());
    expect(await page.evaluate(() => window.__GBDRAW_APP__.linearComparisonPlan.mode)).toBe('none');
    await expect(definition).toHaveValue('HISTORY_RECORD_TWO');
    expect(page.externalRequests).toEqual([]);
  } finally {
    await page.context().close();
  }
});
