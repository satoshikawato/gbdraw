const { test, expect } = require('@playwright/test');
const fs = require('node:fs');
const path = require('node:path');
const { gunzipSync, gzipSync } = require('node:zlib');
const { evaluateWithRetainedPromise } = require('./helpers/app-lifecycle.cjs');
const { openFresh } = require('./helpers/audit-browser.cjs');
const { generate } = require('./helpers/mode-transition.cjs');

// Every Session resource descriptor (`encoding: base64`) anywhere in the document.
const base64Resources = (value, found = []) => {
  if (Array.isArray(value)) value.forEach((member) => base64Resources(member, found));
  else if (value && typeof value === 'object') {
    if (value.encoding === 'base64' && typeof value.data === 'string') found.push(value);
    Object.values(value).forEach((member) => base64Resources(member, found));
  }
  return found;
};

// A gzip GenBank upload (NCBI's .gbff.gz) is inspected and drawn like the plain
// file, and the saved Session holds the plain file, not the gzip bytes.
test('a gzip GenBank upload lists its record, generates, and saves the plain file', async ({ page }, info) => {
  test.setTimeout(240000);
  const gzipped = info.outputPath('NC_001416.gb.gz');
  const plain = fs.readFileSync(path.resolve('tests/test_inputs/NC_001416.gb'));
  fs.writeFileSync(gzipped, gzipSync(plain));
  await openFresh(page);
  await page.getByLabel('GenBank/DDBJ File', { exact: true }).setInputFiles(gzipped);
  await expect.poll(() => page.evaluate(() => (
    window.__GBDRAW_APP__.circularRecordList.map((record) => record.record_id)
  )), { timeout: 60_000 }).toEqual(['NC_001416.1']);
  await generate(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.results.length)).toBe(1);
  // A titled Session saves without the title dialog.
  const pending = page.waitForEvent('download', { timeout: 60000 });
  await evaluateWithRetainedPromise(page, async () => {
    const app = window.__GBDRAW_APP__;
    app.sessionTitle = 'gzip-input';
    await app.saveSessionWithTitle();
  });
  const saved = await (await pending).path();
  const genbank = base64Resources(JSON.parse(gunzipSync(fs.readFileSync(saved)).toString('utf8')))
    .filter((resource) => resource.kind === 'genbank');
  expect(genbank).toHaveLength(1);
  expect(Buffer.from(genbank[0].data, 'base64').equals(plain)).toBe(true);
  // The Web prefixes an embedded resource's name with its resource id.
  expect(genbank[0].name).toMatch(/NC_001416\.gb\.gz$/);
});

// A gzip upload that expands past the cap (Python's MAX_EXPANDED_INPUT_BYTES; a
// test hook stands in for 512 MiB) is refused with its name and the cap, and
// no record is listed (PD-OI-116 in docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md).
test('a gzip upload past the expanded cap is refused with its name and the cap', async ({ page }, info) => {
  test.setTimeout(240000);
  const cap = 1024 * 1024;
  const gzipped = info.outputPath('big genome.gb.gz');
  fs.writeFileSync(gzipped, gzipSync(Buffer.alloc(cap + 1, 0x41)));
  await page.addInitScript((maxExpandedInputBytes) => {
    window.__GBDRAW_TEST_HOOKS__ = { ...(window.__GBDRAW_TEST_HOOKS__ || {}), maxExpandedInputBytes };
  }, cap);
  await openFresh(page);
  await page.getByLabel('GenBank/DDBJ File', { exact: true }).setInputFiles(gzipped);
  await expect(page.getByText('An input could not be read. Replace or reselect it. File: big genome.gb.gz. '
    + 'gbdraw reads gzip files that expand to at most 1 MiB. Decompress the gzip file and choose the plain file.'))
    .toBeVisible({ timeout: 60_000 });
  expect(await page.evaluate(() => window.__GBDRAW_APP__.circularRecordList.length)).toBe(0);
});
