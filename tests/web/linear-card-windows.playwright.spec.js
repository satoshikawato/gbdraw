// The windows of a Linear record card (PD-OI-124): an End typed before the Start shows in order with Reverse
// complement ticked (OV-2884); "+ Add window" turns the Region block into the Windows list, the card counts the bases
// left out between two windows, Generate draws each window, and Remove returns to the Region block. Each action is
// one History step. The Record Layout and the request rules are pinned in linear-card-windows.test.mjs.
const { test, expect } = require('@playwright/test');
const { generate, load, seeds } = require('./helpers/mode-transition.cjs');

test.describe.configure({ retries: 0 });

/** @param {import('@playwright/test').Page} page */
const undoLabels = (page) => page.evaluate(() => window.__GBDRAW_HISTORY__.timeline().undo.map((step) => step.label));

/** @param {import('@playwright/test').Locator} input @param {string} value */
const type = async (input, value) => {
  await input.fill(value);
  await input.press('Tab');
};

/** The source starts and ends of the records the shown Result draws. @param {import('@playwright/test').Page} page */
const drawnSpans = (page) => page.evaluate(() => {
  const svg = new DOMParser().parseFromString(String(window.__GBDRAW_APP__.currentResults()[0].content || ''), 'image/svg+xml');
  return [...svg.querySelectorAll('[data-gbdraw-record-source-start]')].map((node) => [
    Number(node.getAttribute('data-gbdraw-record-source-start')), Number(node.getAttribute('data-gbdraw-record-source-end'))
  ]);
});

test('a Linear card edits its windows: End before Start, Add window, the bases between, Remove, one step each', async ({ browser }) => {
  test.setTimeout(420_000);
  // Lambda (NC_001416, 48,502 bp) is a linear record: its card has no "Across the origin".
  const page = await load(browser, seeds.linear);
  try {
    await generate(page);
    const uid = await page.evaluate(() => window.__GBDRAW_APP__.linearSeqs[0].uid);
    await page.locator(`[data-record-selector-uid="${uid}"]`).evaluate((select) => {
      for (let node = select; node; node = node.parentElement) if (node.tagName === 'DETAILS') node.open = true;
    });
    const region = page.locator('[data-record-region-fields]');
    await expect(region.getByRole('checkbox', { name: 'Across the origin for sequence 1' })).toHaveCount(0);
    const start = region.getByLabel('Region start for sequence 1', { exact: true });
    const end = region.getByLabel('Region end for sequence 1', { exact: true });
    const reverse = region.getByRole('checkbox', { name: 'Reverse complement for sequence 1' });

    // OV-2884: an End before the Start is the reversed span, shown in order with Reverse complement ticked.
    await type(start, '5000');
    await type(end, '1001');
    await expect(start).toHaveValue('1001');
    await expect(end).toHaveValue('5000');
    await expect(reverse).toBeChecked();
    await reverse.uncheck();

    const steps = (await undoLabels(page)).length;
    await page.getByRole('button', { name: 'Add window to sequence 1' }).click();
    const windows = page.locator('[data-record-windows]');
    await expect(windows.locator('[data-record-window-count]')).toHaveText('2 windows of this record');
    await expect(region).toHaveCount(0);
    await expect(windows.locator('[data-record-window]')).toHaveCount(2);
    await expect(windows.getByText('The whole record until Start and End are set.')).toBeVisible();
    const added = await undoLabels(page);
    expect([added.length - steps, added[0]]).toEqual([1, 'Add window']);

    await type(windows.getByLabel('Window 2 start for sequence 1', { exact: true }), '9001');
    await type(windows.getByLabel('Window 2 end for sequence 1', { exact: true }), '12000');
    await expect(windows.locator('[data-record-window-gap]'))
      .toHaveText('4,000 bp of the record lie between windows 1 and 2 and are left out.');
    await generate(page);
    expect(await drawnSpans(page)).toEqual([[1001, 5000], [9001, 12000]]);

    await windows.getByRole('button', { name: 'Remove window 1 of sequence 1' }).click();
    expect((await undoLabels(page))[0]).toBe('Remove window');
    await expect(windows).toHaveCount(0);
    await expect(start).toHaveValue('9001');
    await expect(end).toHaveValue('12000');
    // Undo brings the window back.
    await page.evaluate(async () => { await window.__GBDRAW_HISTORY__.undo(); });
    await expect(page.locator('[data-record-windows] [data-record-window]')).toHaveCount(2);
  } finally {
    await page.context().close();
  }
});
