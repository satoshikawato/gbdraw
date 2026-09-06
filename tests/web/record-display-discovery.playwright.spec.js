const { test, expect } = require('@playwright/test');
const { openApp, generateAndWaitForResult } = require('./helpers/app-lifecycle.cjs');

const genbank = (topology, length = 120) => `LOCUS       same                     ${length} bp    DNA     ${topology.padEnd(8)} UNA 01-JAN-2000
DEFINITION  Source discovery regression fixture.
ACCESSION   same
VERSION     same
KEYWORDS    .
SOURCE      synthetic construct
  ORGANISM  synthetic construct
FEATURES             Location/Qualifiers
     CDS             11..35
                     /product="source feature"
ORIGIN
        1 ${'acgt'.repeat(length / 4)}
//
`;

for (const mode of ['circular', 'linear']) {
  test(`record topology discovery uses uploaded source and packaged Worker in ${mode}`, async ({ page }) => {
    test.setTimeout(180000);
    const external = [];
    await page.route('**/*', route => {
      if (new URL(route.request().url()).hostname === '127.0.0.1') return route.continue();
      external.push(route.request().url());
      return route.abort();
    });
    await openApp(page);
    if (mode === 'linear') await page.getByRole('button', { name: 'Linear', exact: true }).click();
    const upload = mode === 'linear' ? page.getByTestId('linear-genbank-1') : page.getByLabel('GenBank/DDBJ File', { exact: true });
    await upload.setInputFiles({ name: 'same.gbk', mimeType: 'text/plain', buffer: Buffer.from(genbank('circular') + genbank('linear')) });
    if (mode === 'linear') {
      await page.getByRole('button', { name: 'Record options for sequence 1', exact: true }).click();
      const selector = page.getByRole('combobox', { name: 'Record selector for sequence 1', exact: true });
      await expect(selector).toBeEnabled();
      await selector.selectOption('#1');
    }
    const discovery = await page.evaluate(async mode => {
      const app = window.__GBDRAW_APP__;
      const file = mode === 'linear' ? app.linearSeqs[0].gb : app.files.c_gb;
      const { discoverSequenceRecords, normalizeSequenceRecords } = await import('./js/app/record-discovery.js');
      const { runDiagramHelperOperation, DIAGRAM_HELPER_OPERATIONS } = await import('./js/services/diagram-generation.js');
      const fast = await discoverSequenceRecords({ file, format: 'genbank' });
      const response = await runDiagramHelperOperation(DIAGRAM_HELPER_OPERATIONS.LIST_SEQUENCE_RECORDS, {
        format: 'genbank', files: [{ role: 'source', bytes: await file.arrayBuffer() }]
      });
      return { fast, worker: normalizeSequenceRecords(response.result) };
    }, mode);
    expect(discovery.fast).toEqual(discovery.worker);
    expect(discovery.fast.map(r => [r.selector, r.recordId, r.detectedTopology])).toEqual([
      ['#1', 'same', 'circular'], ['#2', 'same', 'linear']
    ]);
    await generateAndWaitForResult(page);
    expect(await page.evaluate(() => window.__GBDRAW_APP__.results.length)).toBeGreaterThan(0);
    // New controls cannot be activated until the shared writer is integrated.
    await expect(page.getByRole('button', { name: 'Use selected feature midpoint', exact: true })).toHaveCount(0);
    expect(external).toEqual([]);
  });
}
