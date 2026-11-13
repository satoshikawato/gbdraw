const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { join } = require('node:path');
const { CURRENT_FEATURE_CATALOG_SCHEMA } = require('./helpers/app-lifecycle.cjs');

const EXPORT_MODULE_PATH = '/gbdraw/web/js/services/export.js';
const STANDALONE_MODULE_PATH = '/gbdraw/web/js/services/standalone-interactivity.js';
const STANDALONE_ASSET_PATH = '/gbdraw/web/js/services/standalone-interactivity-assets.js';
const JSPDF_PATH = '/gbdraw/web/vendor/jspdf/jspdf.umd.min.js';
const SVG2PDF_PATH = '/gbdraw/web/vendor/svg2pdf.js/svg2pdf.umd.min.js';

const trackedExportAsset = (url) => [
  EXPORT_MODULE_PATH,
  STANDALONE_MODULE_PATH,
  STANDALONE_ASSET_PATH,
  JSPDF_PATH,
  SVG2PDF_PATH
].find((path) => new URL(url).pathname.endsWith(path)) || '';

const installRequestCounter = (page) => {
  const counts = new Map();
  page.on('request', (request) => {
    const asset = trackedExportAsset(request.url());
    if (asset) counts.set(asset, (counts.get(asset) || 0) + 1);
  });
  return (path) => counts.get(path) || 0;
};

// An Interactive SVG holds every sequence, which Python reads from the request the shown Result was made from
// (PD-OI-128), so the cases that write one load a Gallery Session: its Result keeps its request.
const loadGallerySession = async (page) => {
  await page.locator('input[accept*=".gbdraw"]').setInputFiles(
    join(process.cwd(), 'gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw')
  );
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending && window.__GBDRAW_APP__.extractedFeatures.length);
  return page.evaluate(() => `${window.__GBDRAW_APP__.results[0].name.replace(/\.svg$/, '')}.interactive.svg`);
};

const displayedErrorOperation = (page) => page.evaluate(
  () => window.__GBDRAW_APP__.errorDisplay?.operation || null
);

// The fixture is a legacy imported SVG without composition metadata. Mounting
// it must leave the displayed error unchanged.
const mountExportFixture = async (page, { interactive = false, size = { width: 120, height: 80 } } = {}) => {
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  const errorBeforeMount = await displayedErrorOperation(page);
  await page.evaluate(async ({ interactive, catalogSchema, size }) => {
    const { state } = await import('/gbdraw/web/js/state.js');
    const { admitLegacyImportedResults, createLegacyImportResultSource } = await import(
      '/gbdraw/web/js/services/svg-result-ingestion.js'
    );
    const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
    svg.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
    svg.setAttribute('width', String(size.width));
    svg.setAttribute('height', String(size.height));
    svg.setAttribute('viewBox', `0 0 ${size.width} ${size.height}`);

    const rect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
    rect.setAttribute('id', 'feature-a');
    rect.setAttribute('data-gbdraw-feature-id', 'biological-a');
    rect.setAttribute('data-gbdraw-rendered-feature-id', 'feature-a');
    rect.setAttribute('x', '10');
    rect.setAttribute('y', '10');
    rect.setAttribute('width', '60');
    rect.setAttribute('height', '30');
    rect.setAttribute('fill', '#54bcf8');
    svg.appendChild(rect);

    const committedResult = admitLegacyImportedResults(
      createLegacyImportResultSource([{
        name: 'lazy-export.svg',
        content: new XMLSerializer().serializeToString(svg)
      }])
    )[0];
    state.results.value = [committedResult];
    state.selectedResultIndex.value = 0;
    state.downloadDpi.value = 96;
    state.featureCatalog.value = interactive ? {
      schema: catalogSchema,
      items: [{
        resultIndex: 0,
        resultName: 'lazy-export.svg',
        recordKeys: ['record-a'],
        features: [{
          svgId: 'feature-a',
          recordKey: 'record-a',
          biologicalFeatureId: 'biological-a'
        }],
        biologicalFeatures: [{
          recordKey: 'record-a',
          biologicalFeatureId: 'biological-a',
          type: 'CDS',
          start: 0,
          end: 9,
          anchorProfile: {
            precision: 'exact', operator: 'single', partOrder: 'source-forward', strand: 'unstranded'
          },
          product: 'Lazy export feature'
        }],
        orthogroups: [],
        annotations: [],
        comparisonMatches: []
      }]
    } : { schema: catalogSchema, items: [] };
  }, { interactive, catalogSchema: CURRENT_FEATURE_CATALOG_SCHEMA, size });
  await expect(page.locator('.origin-top svg')).toBeAttached();
  // The preview binds the mounted root in microtasks; two frames cover it.
  await page.evaluate(() => new Promise((resolve) => requestAnimationFrame(() => requestAnimationFrame(resolve))));
  expect(await displayedErrorOperation(page)).toBe(errorBeforeMount);
};

const expectNoPdfOrInteractiveRequests = (requestCount) => {
  expect(requestCount(JSPDF_PATH)).toBe(0);
  expect(requestCount(SVG2PDF_PATH)).toBe(0);
  expect(requestCount(STANDALONE_MODULE_PATH)).toBe(0);
  expect(requestCount(STANDALONE_ASSET_PATH)).toBe(0);
};

for (const failedAsset of [JSPDF_PATH, SVG2PDF_PATH]) {
  test(`PDF retries after a transient ${failedAsset.split('/').at(-1)} failure`, async ({ page }) => {
    await page.goto('/gbdraw/web/index.html');
    await mountExportFixture(page);
    await page.route(`**${failedAsset}`, (route) => route.abort());
    expect(await page.evaluate(() => window.__GBDRAW_APP__.downloadPDF())).toMatchObject({ status: 'error', error: {code:'PDF_LIBRARY',operation:'export-pdf'} });
    await page.unroute(`**${failedAsset}`);
    const pending = page.waitForEvent('download');
    await page.evaluate(() => window.__GBDRAW_APP__.downloadPDF());
    const download = await pending;
    expect(download.suggestedFilename()).toBe('lazy-export.pdf');
    expect(await download.failure()).toBeNull();
  });
}

test('startup and non-PDF exports do not load PDF or interactive payloads', { tag: '@pr-smoke' }, async ({
  page
}) => {
  const requestCount = installRequestCounter(page);
  await page.goto('/gbdraw/web/index.html');
  await page.waitForFunction(() => window.__GBDRAW_APP__);

  expect(requestCount(EXPORT_MODULE_PATH)).toBe(0);
  expectNoPdfOrInteractiveRequests(requestCount);

  await mountExportFixture(page);

  const svgDownload = page.waitForEvent('download');
  await page.evaluate(async () => window.__GBDRAW_APP__.downloadSVG());
  expect((await svgDownload).suggestedFilename()).toBe('lazy-export.svg');
  expect(requestCount(EXPORT_MODULE_PATH)).toBe(1);
  expectNoPdfOrInteractiveRequests(requestCount);

  const pngDownload = page.waitForEvent('download');
  await page.evaluate(async () => window.__GBDRAW_APP__.downloadPNG());
  expect((await pngDownload).suggestedFilename()).toBe('lazy-export.png');
  expect(requestCount(EXPORT_MODULE_PATH)).toBe(1);
  expectNoPdfOrInteractiveRequests(requestCount);
});

// OV-2682: a PNG is the size the command line draws for the same figure and DPI, the size Python's check returns
// (gbdraw/render/export_limits.py, rounded half to even as CairoSVG rounds). At 300 dpi, 1,237 x 100 SVG px is
// 3,865.625 x 312.5: a canvas that truncated drew 3,865 px wide, and a JavaScript Math.round would draw 313 px tall.
test('a PNG download has the pixel size of the command line PNG of the same figure', async ({ page }) => {
  await page.goto('/gbdraw/web/index.html');
  await mountExportFixture(page, { size: { width: 1237, height: 100 } });
  await page.evaluate(async () => {
    const { state } = await import('/gbdraw/web/js/state.js');
    state.downloadDpi.value = 300;
  });
  const pngDownload = page.waitForEvent('download', { timeout: 120000 });
  await page.evaluate(async () => window.__GBDRAW_APP__.downloadPNG());
  const png = readFileSync(await (await pngDownload).path());
  expect(png.subarray(12, 16).toString('latin1')).toBe('IHDR');
  expect([png.readUInt32BE(16), png.readUInt32BE(20)]).toEqual([3866, 312]);
});

test('concurrent first PDF exports load each library once and repeated exports reuse them', async ({
  page
}) => {
  const requestCount = installRequestCounter(page);
  await page.goto('/gbdraw/web/index.html');
  await mountExportFixture(page);

  const firstDownloads = [];
  page.on('download', (download) => firstDownloads.push(download));
  await page.evaluate(async () => Promise.all([
    window.__GBDRAW_APP__.downloadPDF(),
    window.__GBDRAW_APP__.downloadPDF()
  ]));
  await expect.poll(() => firstDownloads.length).toBe(2);
  expect(firstDownloads.map((download) => download.suggestedFilename()))
    .toEqual(['lazy-export.pdf', 'lazy-export.pdf']);
  expect(requestCount(EXPORT_MODULE_PATH)).toBe(1);
  expect(requestCount(JSPDF_PATH)).toBe(1);
  expect(requestCount(SVG2PDF_PATH)).toBe(1);
  expect(requestCount(STANDALONE_MODULE_PATH)).toBe(0);
  expect(requestCount(STANDALONE_ASSET_PATH)).toBe(0);

  const secondDownload = page.waitForEvent('download');
  await page.evaluate(async () => window.__GBDRAW_APP__.downloadPDF());
  expect((await secondDownload).suggestedFilename()).toBe('lazy-export.pdf');
  expect(requestCount(EXPORT_MODULE_PATH)).toBe(1);
  expect(requestCount(JSPDF_PATH)).toBe(1);
  expect(requestCount(SVG2PDF_PATH)).toBe(1);
});

test('interactive assets load only for interactive SVG export and the file is standalone', async ({
  page
}) => {
  test.setTimeout(180000);
  const requestCount = installRequestCounter(page);
  await page.goto('/gbdraw/web/index.html');
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  const filename = await loadGallerySession(page);

  expect(requestCount(EXPORT_MODULE_PATH)).toBe(0);
  expectNoPdfOrInteractiveRequests(requestCount);

  const downloadPromise = page.waitForEvent('download');
  await page.evaluate(async () => window.__GBDRAW_APP__.downloadInteractiveSVG());
  const download = await downloadPromise;
  expect(download.suggestedFilename()).toBe(filename);
  const source = await download.createReadStream().then(async (stream) => {
    const chunks = [];
    for await (const chunk of stream) chunks.push(chunk);
    return Buffer.concat(chunks).toString('utf8');
  });

  expect(requestCount(EXPORT_MODULE_PATH)).toBe(1);
  expect(requestCount(STANDALONE_MODULE_PATH)).toBe(1);
  expect(requestCount(STANDALONE_ASSET_PATH)).toBe(1);
  expect(requestCount(JSPDF_PATH)).toBe(0);
  expect(requestCount(SVG2PDF_PATH)).toBe(0);
  expect(source).toContain('id="gbdraw-interactive-feature-metadata"');
  expect(source).toContain('id="gbdraw-interactive-feature-style"');
  expect(source).toContain('id="gbdraw-interactive-feature-script"');
  expect(source).toContain('data-gbdraw-interactive-svg="true"');
});

// The legacy imported SVG has no request: there is no input to read its sequences from until a Generate, so the
// export writes no file and says why (PD-OI-128).
test('Interactive SVG export of a Result without its request writes no file and says why', async ({ page }) => {
  const requestCount = installRequestCounter(page);
  let downloads = 0;
  page.on('download', () => { downloads += 1; });
  await page.goto('/gbdraw/web/index.html');
  await mountExportFixture(page, { interactive: true });
  const outcome = await page.evaluate(async () => {
    const result = await window.__GBDRAW_APP__.downloadInteractiveSVG();
    return { result, summary: window.__GBDRAW_APP__.errorDisplay?.summary || '' };
  });
  expect(outcome.result).toMatchObject({ status: 'error', error: { code: 'SEQUENCE_REQUIRES_GENERATE' } });
  expect(outcome.summary).toContain('its sequences are read after the next Generate');
  await expect(page.getByRole('alert', { name: 'SVG export error' }))
    .toContainText('its sequences are read after the next Generate');
  expect(requestCount(STANDALONE_MODULE_PATH)).toBe(0);
  expect(downloads).toBe(0);
});

test('PDF library load failures are shown to the user', async ({ page }) => {
  await page.route(`**${JSPDF_PATH}`, (route) => route.abort('failed'));
  await page.goto('/gbdraw/web/index.html');
  await mountExportFixture(page);

  const outcome = await page.evaluate(async () => {
    const result = await window.__GBDRAW_APP__.downloadPDF();
    return {
      result,
      summary: window.__GBDRAW_APP__.errorDisplay?.summary || ''
    };
  });

  expect(outcome.result.status).toBe('error');
  expect(outcome.result.error).toMatchObject({code:'PDF_LIBRARY',operation:'export-pdf',stage:'initialization'});
  expect(outcome.summary).toContain('PDF libraries could not start');
  await expect(page.getByRole('alert', { name: 'PDF export error' }))
    .toContainText('PDF libraries could not start');
});

test('export failures are shown to the user', async ({ page }) => {
  await page.goto('/gbdraw/web/index.html');
  await page.waitForFunction(() => window.__GBDRAW_APP__);

  const outcome = await page.evaluate(async () => {
    const result = await window.__GBDRAW_APP__.downloadSVG();
    return {
      result,
      summary: window.__GBDRAW_APP__.errorDisplay?.summary || ''
    };
  });

  expect(outcome.result.status).toBe('error');
  expect(outcome.result.error).toMatchObject({code:'EXPORT_INPUT',operation:'export-svg',stage:'export-capture'});
  expect(outcome.summary).toContain('current diagram could not be prepared for export');
  await expect(page.getByRole('alert', { name: 'SVG export error' }))
    .toContainText('current diagram could not be prepared for export');
});


test('Interactive SVG retry keeps the actual standalone export operation',async({page})=>{
  test.setTimeout(180000);
  await page.goto('/gbdraw/web/index.html');
  await page.waitForFunction(()=>window.__GBDRAW_APP__);
  const filename=await loadGallerySession(page);
  // The first export finds no diagram: the loaded Result is taken off the page; the Session's request stays.
  await page.evaluate(async()=>{
    const {state}=await import('/gbdraw/web/js/state.js');
    window.__HELD_EXPORT_RESULTS__=state.results.value;
    state.results.value=[];
  });
  const failure=await page.evaluate(()=>window.__GBDRAW_APP__.downloadInteractiveSVG());
  expect(failure).toMatchObject({status:'error',error:{code:'EXPORT_INPUT',operation:'export-svg',stage:'export-capture'}});
  const alert=page.getByRole('alert',{name:'SVG export error'});
  await expect(alert.getByRole('button',{name:'Retry Interactive SVG export',exact:true})).toBeVisible();
  await expect(alert.getByRole('button',{name:'Retry SVG export',exact:true})).toHaveCount(0);
  await page.evaluate(async()=>{
    const {state}=await import('/gbdraw/web/js/state.js');
    state.results.value=window.__HELD_EXPORT_RESULTS__;
    state.selectedResultIndex.value=0;
  });
  await expect(page.locator('.origin-top svg')).toBeAttached();
  const pending=page.waitForEvent('download');
  await alert.getByRole('button',{name:'Retry Interactive SVG export',exact:true}).click();
  const downloaded=await pending;
  expect(downloaded.suggestedFilename()).toBe(filename);
  const stream=await downloaded.createReadStream();const chunks=[];
  for await(const chunk of stream)chunks.push(chunk);
  expect(Buffer.concat(chunks).toString()).toContain('id="gbdraw-interactive-feature-script"');
  await expect(alert).toHaveCount(0);
});
