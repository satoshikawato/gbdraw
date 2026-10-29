const { test, expect } = require('@playwright/test');
const { load, generate } = require('./helpers/mode-transition.cjs');
const { evaluateWithRetainedPromise, openLinearFileCards, settleLive } = require('./helpers/app-lifecycle.cjs');

test('comparison off has its own Undo step after a record definition edit', async ({ browser }) => {
  test.setTimeout(300000);
  const page = await load(browser, 'tests/test_inputs/BGC0000708-BGC0000713.gbdraw-session.json');
  try {
    await openLinearFileCards(page);
    await generate(page);
    const baseline = await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
    await page.getByRole('button', { name: 'Record options for sequence 2', exact: true }).click();
    const definition = page.getByLabel('Definition for sequence 2', { exact: true });
    await definition.fill('HISTORY_RECORD_TWO');
    await definition.press('Tab');
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(baseline + 1);
    // A record name is a live edit (PD-OI-103 rev 2): the Result draws it
    // without Generate, and the comparison step starts from that Result.
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.currentResults()[0].content.includes('HISTORY_RECORD_TWO')),
      { timeout: 120000 }).toBe(true);
    await settleLive(page);
    const before = await page.evaluate(async () => ({
      plan: (await import('./js/services/config.js')).buildConfigData((await import('./js/state.js')).state.activeDrawing()).linearComparisonPlan,
      result: window.__GBDRAW_APP__.currentResults()[0].content
    }));
    await page.getByRole('button', { name: 'Set no comparison', exact: true }).click();
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.linearComparisonPlan.mode)).toBe('none');
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(baseline + 2);
    await page.getByRole('button', { name: 'Undo', exact: true }).click();
    await expect.poll(() => page.evaluate(async () =>
      (await import('./js/services/config.js')).buildConfigData((await import('./js/state.js')).state.activeDrawing()).linearComparisonPlan))
      .toEqual(before.plan);
    await expect(definition).toHaveValue('HISTORY_RECORD_TWO');
    expect(await page.evaluate(() => window.__GBDRAW_APP__.currentResults()[0].content)).toBe(before.result);
    await page.getByRole('button', { name: 'Redo', exact: true }).click();
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.linearComparisonPlan.mode)).toBe('none');
    await expect(definition).toHaveValue('HISTORY_RECORD_TWO');
    await page.getByRole('button', { name: 'Set no comparison', exact: true }).click();
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(baseline + 2);
    const losat = page.getByRole('button', { name: 'Run LOSAT for all adjacent pairs', exact: true });
    await losat.focus();
    await losat.press('Enter');
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(baseline + 3);
    await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
    expect(await page.evaluate(() => window.__GBDRAW_APP__.linearComparisonPlan.mode)).toBe('none');
    await expect(definition).toHaveValue('HISTORY_RECORD_TWO');
    expect(page.externalRequests).toEqual([]);
  } finally {
    await page.context().close();
  }
});
