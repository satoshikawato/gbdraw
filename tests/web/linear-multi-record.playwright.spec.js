const { test, expect } = require('@playwright/test');
const { createReadStream, existsSync } = require('node:fs');
const { createServer } = require('node:http');
const { extname, join, normalize, resolve, sep } = require('node:path');

const repoRoot = resolve(process.env.GBDRAW_REPO || process.cwd());
const contentTypes = {
  '.html': 'text/html; charset=utf-8',
  '.js': 'text/javascript; charset=utf-8',
  '.css': 'text/css; charset=utf-8',
  '.json': 'application/json; charset=utf-8',
  '.svg': 'image/svg+xml',
  '.wasm': 'application/wasm',
  '.whl': 'application/octet-stream'
};

let server;
let baseUrl;

test.beforeAll(async () => {
  await new Promise((resolveServer, rejectServer) => {
    server = createServer((request, response) => {
      const url = new URL(request.url || '/', 'http://127.0.0.1');
      const requestedPath = normalize(decodeURIComponent(url.pathname)).replace(/^(\.\.(?:\/|\\|$))+/, '');
      const filePath = resolve(repoRoot, requestedPath.replace(/^[/\\]+/, ''));
      if ((!filePath.startsWith(`${repoRoot}${sep}`) && filePath !== repoRoot) || !existsSync(filePath)) {
        response.writeHead(404);
        response.end('Not found');
        return;
      }
      response.writeHead(200, {
        'Content-Type': contentTypes[extname(filePath)] || 'application/octet-stream'
      });
      createReadStream(filePath).pipe(response);
    });
    server.once('error', rejectServer);
    server.listen(0, '127.0.0.1', () => {
      baseUrl = `http://127.0.0.1:${server.address().port}`;
      resolveServer();
    });
  });
});

test.afterAll(async () => {
  await new Promise((resolveClose) => server.close(resolveClose));
});

test('Linear record rows and N-to-M comparison batches remain keyed by sequence uid', async ({ page }) => {
  await page.goto(`${baseUrl}/gbdraw/web/index.html`, { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);

  const state = await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    app.mode = 'linear';
    app.addLinearSeq();
    app.addLinearSeq();
    app.addLinearSeq();
    app.linearRecordLayoutEnabled = true;
    app.syncLinearRecordLayout();
    app.setLinearRecordRow(app.linearSeqs[0].uid, 1);
    app.setLinearRecordRow(app.linearSeqs[1].uid, 1);
    app.setLinearRecordRow(app.linearSeqs[2].uid, 2);
    app.setLinearRecordRow(app.linearSeqs[3].uid, 2);
    app.addLinearComparisonBatch(true);
    return {
      tokens: app.linearLayoutTokens,
      comparisonCount: app.linearComparisons.length,
      endpoints: app.linearComparisons.map((item) => [item.queryUid, item.subjectUid]),
      uniqueUids: new Set(app.linearSeqs.map((item) => item.uid)).size
    };
  });

  expect(state.tokens).toEqual(['#1@1', '#2@1', '#3@2', '#4@2']);
  expect(state.comparisonCount).toBe(4);
  expect(state.uniqueUids).toBe(4);
  expect(new Set(state.endpoints.map((pair) => pair.join('->'))).size).toBe(4);
  await expect(page.locator('input[aria-label="Linear record row"]')).toHaveCount(4);
  await expect(page.getByText('All adjacent-row pairs', { exact: true })).toBeVisible();
});

test('Region annotations expose and persist an explicit target-record selection', async ({ page }) => {
  await page.goto(`${baseUrl}/gbdraw/web/index.html`, { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);

  const genbank = `LOCUS       RecA                      10 bp    DNA     linear   UNA 01-JAN-2000
DEFINITION  first.
ACCESSION   RecA
VERSION     RecA
KEYWORDS    .
SOURCE      synthetic construct
  ORGANISM  synthetic construct
            .
FEATURES             Location/Qualifiers
ORIGIN
        1 aaaaaaaaaa
//
LOCUS       RecB                      12 bp    DNA     linear   UNA 01-JAN-2000
DEFINITION  second.
ACCESSION   RecB
VERSION     RecB
KEYWORDS    .
SOURCE      synthetic construct
  ORGANISM  synthetic construct
            .
FEATURES             Location/Qualifiers
ORIGIN
        1 cccccccccccc
//
`;
  await page.evaluate((content) => {
    const app = window.__GBDRAW_APP__;
    app.mode = 'linear';
    app.lInputType = 'gb';
    app.setLinearSeqPrimaryFile(0, 'gb', new File([content], 'two-records.gb', {
      type: 'text/plain',
      lastModified: 1
    }));
    const set = app.addAnnotationSet('review');
    app.addCoordinateAnnotation(set, { start: 1, end: 10 });
  }, genbank);

  await page.getByText('Region Annotations', { exact: false }).click();
  const selector = page.getByLabel('Annotation target record');
  await expect(selector).toHaveCount(1);
  await expect(selector).toHaveValue('');
  await expect(selector.locator('option')).toHaveText([
    'Select target record',
    '#1 · RecA · 10 bp',
    '#2 · RecB · 12 bp'
  ], { timeout: 60000 });
  await expect(page.getByText('Choose the record that this annotation targets.')).toBeVisible();

  const rejected = await page.evaluate(() => window.__GBDRAW_APP__.runAnalysis());
  expect(rejected).toEqual({ status: 'error' });
  await expect(page.getByText('Choose a target record for region annotation review/region_1.')).toBeVisible();

  await selector.selectOption({ label: '#2 · RecB · 12 bp' });
  await expect(page.getByText('Choose the record that this annotation targets.')).toHaveCount(0);
  const target = await page.evaluate(() => (
    window.__GBDRAW_APP__.annotationSets[0].annotations[0].target.record
  ));
  expect(target).toEqual({ kind: 'recordId', value: 'RecB' });
});

test('GFF annotation targets follow FASTA record order', async ({ page }) => {
  await page.goto(`${baseUrl}/gbdraw/web/index.html`, { waitUntil: 'domcontentloaded' });
  await page.waitForFunction(() => window.__GBDRAW_APP__);

  const gff = `##gff-version 3
##sequence-region RecB 1 12
RecB\ttest\tgene\t1\t3\t.\t+\t.\tID=gene_b
##sequence-region RecA 1 10
RecA\ttest\tgene\t2\t4\t.\t+\t.\tID=gene_a
`;
  const fasta = `>RecB
CCCCCCCCCCCC
>RecA
AAAAAAAAAA
`;
  await page.evaluate(({ gffText, fastaText }) => {
    const app = window.__GBDRAW_APP__;
    app.mode = 'linear';
    app.lInputType = 'gff';
    app.setLinearSeqPrimaryFile(0, 'gff', new File([gffText], 'records.gff3', {
      type: 'text/plain',
      lastModified: 2
    }));
    app.setLinearSeqPrimaryFile(0, 'fasta', new File([fastaText], 'records.fasta', {
      type: 'text/plain',
      lastModified: 3
    }));
    const set = app.addAnnotationSet('gff-review');
    app.addCoordinateAnnotation(set, { start: 1, end: 3 });
  }, { gffText: gff, fastaText: fasta });

  await page.getByText('Region Annotations', { exact: false }).click();
  await expect(page.getByLabel('Annotation target record').locator('option')).toHaveText([
    'Select target record',
    '#1 · RecB · 12 bp',
    '#2 · RecA · 10 bp'
  ], { timeout: 60000 });
});
