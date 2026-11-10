// The records a drawing draws and where a cropped record starts, as the sidebar shows them before Generate.
// OV-3340: the status line under the Circular file card counts the records the next Generate draws: the record the
// Record select names alone, else the ON records. OV-3341: while Display start is unavailable for a cropped record,
// its field shows where the drawn crop starts (the Circular 12 o'clock, the Linear left end), not an end of the whole
// record.
const { test, expect } = require('@playwright/test');
const { evaluateWithRetainedPromise, openApp } = require('./helpers/app-lifecycle.cjs');
const { showDrawing } = require('./helpers/drawing-bar.cjs');

test.describe.configure({ retries: 0 });

const LENGTH = 4000;
const record = (id) => {
  const sequence = 'acgt'.repeat(LENGTH / 4);
  const origin = sequence.match(/.{1,60}/g).map((chunk, index) => (
    `${String(index * 60 + 1).padStart(9)} ${chunk.match(/.{1,10}/g).join(' ')}`
  )).join('\n');
  return `LOCUS       ${id.padEnd(16)} ${String(LENGTH).padStart(8)} bp    DNA     circular UNA 01-JAN-2000
DEFINITION  drawn record surface test.
ACCESSION   ${id}
VERSION     ${id}
KEYWORDS    .
SOURCE      synthetic construct
  ORGANISM  synthetic construct
            .
FEATURES             Location/Qualifiers
     CDS             1..90
                     /locus_tag="${id}_1"
ORIGIN
${origin}
//
`;
};

const openCircularTwoRecords = async (page) => {
  await openApp(page);
  await page.evaluate((text) => {
    const app = window.__GBDRAW_APP__;
    app.form.multi_record_canvas = false;
    void app.addCircularSourceFiles([new File([text], 'two.gb', { type: 'text/plain', lastModified: 1 })]);
  }, record('TESTA') + record('TESTB'));
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.circularRecordDiscoveryState.status)).toBe('ready');
};
const selectTestB = async (page) => {
  const select = page.getByRole('combobox', { name: 'Circular record', exact: true });
  const value = await select.evaluate((element) => [...element.options].find((option) => option.text.startsWith('TESTB')).value);
  await select.selectOption(value);
};

test('the status line counts the record the Record select draws alone (OV-3340)', async ({ page }) => {
  test.setTimeout(180_000);
  await openCircularTwoRecords(page);
  const drawn = page.locator('[data-circular-records-drawn]');
  await expect(drawn).toHaveText('2 of 2 records drawn');
  await selectTestB(page);
  await expect(drawn).toHaveText('1 of 2 records drawn');
  await expect(page.locator('[data-circular-discovery-status]')).toHaveText('2 source records inspected · 1 of 2 records drawn');
  // The automatic option still describes the all-records choice.
  const select = page.getByRole('combobox', { name: 'Circular record', exact: true });
  await expect(select.locator('option').first()).toHaveText('All records (separate diagrams)');
  await expect(page.getByRole('spinbutton', { name: /^Display start/ })).toHaveCount(1);
});

for (const [reverse, start] of [[false, '500'], [true, '3000']]) {
  test(`a cropped Circular record shows its crop start ${start} in Display start (OV-3341)`, async ({ page }) => {
    test.setTimeout(180_000);
    await openCircularTwoRecords(page);
    await selectTestB(page);
    await page.evaluate((reverse) => Object.assign(window.__GBDRAW_APP__.form, {
      circular_region_start: 500, circular_region_end: 3000, circular_reverse: reverse
    }), reverse);
    const field = page.getByRole('spinbutton', { name: 'Display start TESTB #2', exact: true });
    await expect(field).toBeDisabled();
    await expect(field).toHaveValue(start);
    await expect(page.getByText('Display start is unavailable for a cropped record.')).toBeVisible();
  });
}

test('a Linear card drawn in a window shows the window start in Display start (OV-3341)', async ({ page }) => {
  test.setTimeout(180_000);
  await openApp(page);
  await showDrawing(page, { kind: 'linear' });
  await evaluateWithRetainedPromise(page, async (text) => {
    const app = window.__GBDRAW_APP__;
    app.setLinearInputType('gb');
    app.setLinearSeqPrimaryFile(0, 'gb', new File([text], 'one.gb', { type: 'text/plain', lastModified: 1 }));
    await app.refreshLinearRecordSelectors();
  }, record('TESTL'));
  await page.getByRole('button', { name: 'Record options for sequence 1', exact: true }).click();
  const field = page.getByRole('spinbutton', { name: 'Display start TESTL #1', exact: true });
  await expect(field).toBeEnabled();
  await expect(field).toHaveValue('1');
  for (const [reverse, start] of [[false, '500'], [true, '3000']]) {
    await page.evaluate((reverse) => Object.assign(window.__GBDRAW_APP__.linearSeqs[0], {
      region_start: 500, region_end: 3000, region_reverse: reverse
    }), reverse);
    await expect(field).toBeDisabled();
    await expect(field).toHaveValue(start);
  }
});
