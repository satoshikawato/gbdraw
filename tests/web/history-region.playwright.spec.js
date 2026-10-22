const { test, expect } = require('@playwright/test');
const { load, generate } = require('./helpers/mode-transition.cjs');
const { generateAndWaitForResult, evaluateWithRetainedPromise } = require('./helpers/app-lifecycle.cjs');

const region = page => page.evaluate(() => {
  const form = window.__GBDRAW_APP__.form;
  return [form.circular_region_start, form.circular_region_end];
});
// Generate may already open the presentation panel; open it only when closed.
const openPresentation = async (page) => {
  const details = page.locator('[data-circular-record-presentation]');
  if (!await details.evaluate((element) => element.open)) await details.locator('summary').click();
};
const edit = async (page, label, value) => {
  const before = await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
  const input = page.getByLabel(label, { exact: true });
  await input.fill(value);
  await input.press('Tab');
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(before + 1);
};

test('Circular region restores number to Auto', async ({ browser }) => {
  test.setTimeout(300000);
  const page = await load(browser);
  try {
    await generate(page);
    await openPresentation(page);
    expect(await region(page)).toEqual([null, null]);
    await edit(page, 'Circular region start', '1000');
    await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
    expect(await region(page)).toEqual([null, null]);
    await expect(page.getByLabel('Circular region start', { exact: true })).toHaveValue('');
    await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.redo());
    expect(await region(page)).toEqual([1000, null]);
    expect(page.externalRequests).toEqual([]);
  } finally { await page.context().close(); }
});

test('rejected Circular region can Undo both edits back to valid Auto', async ({ browser }) => {
  test.setTimeout(300000);
  const page = await load(browser);
  try {
    await generate(page);
    const before = await page.evaluate(() => ({
      undo: window.__GBDRAW_HISTORY__.getUndoCount(), result: window.__GBDRAW_APP__.results[0].content
    }));
    await openPresentation(page);
    await edit(page, 'Circular region start', '1000');
    await edit(page, 'Circular region end', '500');
    await generateAndWaitForResult(page, { expectedStatus: 'error' });
    expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(before.undo + 2);
    await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
    await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
    expect(await region(page)).toEqual([null, null]);
    await expect(page.getByLabel('Circular region start', { exact: true })).toHaveValue('');
    await expect(page.getByLabel('Circular region end', { exact: true })).toHaveValue('');
    await generate(page);
    expect(await page.evaluate(() => window.__GBDRAW_APP__.results[0].content)).toBe(before.result);
    expect(page.externalRequests).toEqual([]);
  } finally { await page.context().close(); }
});

// OV-1097 (Owner decision 2026-10-18): a Region End past the record end is
// drawn to the record end, as Python draws it, and the note beside the Region
// fields says so; a Start past the record end stays an error. The fields keep
// the typed values.
test('a Region End past the record end is drawn to the record end, and the note says so', async ({ browser }) => {
  test.setTimeout(300000);
  const page = await load(browser);
  try {
    await generate(page);
    await openPresentation(page);
    const length = await page.evaluate(() => Number(window.__GBDRAW_APP__.circularRecordList[0].record_length));
    const text = (value) => value.toLocaleString('en-US');
    const note = page.locator('#circular-region-note');
    await edit(page, 'Circular region start', '10000');
    await expect(note).toHaveCount(0);
    await edit(page, 'Circular region end', String(length + 1000));
    await expect(note).toHaveText(
      `Region reaches ${text(length + 1000)}, past the record end (${text(length)} bp); Generate draws to ${text(length)}.`
    );
    await expect(page.getByLabel('Circular region end', { exact: true })).toHaveAttribute('aria-describedby', /circular-region-note/);
    await generateAndWaitForResult(page);
    expect(await region(page)).toEqual([10000, length + 1000]);

    await edit(page, 'Circular region start', String(length + 1));
    await expect(note).toHaveCount(0);
    await generateAndWaitForResult(page, { expectedStatus: 'error', requireCommittedResult: false });
    expect(await page.evaluate(() => {
      const { code, context } = window.__GBDRAW_APP__.errorLog || {};
      return { code, reason: context?.reason };
    })).toEqual({ code: 'REGION_INVALID', reason: 'RECORD_BOUNDS' });

    // Linear: the note follows the record the entry draws.
    const linear = await load(browser, 'gbdraw/web/gallery/sessions/lambda_basic_linear.gbdraw-session.json');
    try {
      // Generate reads the entry's records, so their length is known.
      await generateAndWaitForResult(linear);
      const seq = await linear.evaluate(() => window.__GBDRAW_APP__.linearSeqs[0].uid);
      const end = linear.getByLabel('Region end for sequence 1', { exact: true });
      const start = linear.getByLabel('Region start for sequence 1', { exact: true });
      for (const input of [start, end]) {
        if (!await input.isVisible()) await linear.locator(`[data-record-selector-uid="${seq}"]`).evaluate((select) => {
          for (let node = select; node; node = node.parentElement) if (node.tagName === 'DETAILS') node.open = true;
        });
      }
      await start.fill('40000');
      await start.press('Tab');
      await end.fill('60000');
      await end.press('Tab');
      await expect(linear.locator(`#linear-region-note-${seq}`)).toHaveText(
        /^Region reaches 60,000, past the record end \(48,502 bp\); Generate draws to 48,502\.$/ // lambda: 48,502 bp
      );
      await generateAndWaitForResult(linear);
    } finally { await linear.context().close(); }
  } finally { await page.context().close(); }
});
