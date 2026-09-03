const { test, expect } = require('@playwright/test');

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

const mountExportFixture = async (page, { interactive = false } = {}) => {
  await page.waitForFunction(() => window.__GBDRAW_APP__);
  await page.evaluate(async ({ interactive }) => {
    const { state } = await import('/gbdraw/web/js/state.js');
    const { admitLegacyImportedResults, createLegacyImportResultSource } = await import(
      '/gbdraw/web/js/services/svg-result-ingestion.js'
    );
    const svg = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
    svg.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
    svg.setAttribute('width', '120');
    svg.setAttribute('height', '80');
    svg.setAttribute('viewBox', '0 0 120 80');

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
      schema: 3,
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
          product: 'Lazy export feature'
        }],
        orthogroups: [],
        annotations: [],
        comparisonMatches: []
      }]
    } : { schema: 3, items: [] };
  }, { interactive });
  await expect(page.locator('.origin-top svg')).toBeAttached();
};

const expectNoPdfOrInteractiveRequests = (requestCount) => {
  expect(requestCount(JSPDF_PATH)).toBe(0);
  expect(requestCount(SVG2PDF_PATH)).toBe(0);
  expect(requestCount(STANDALONE_MODULE_PATH)).toBe(0);
  expect(requestCount(STANDALONE_ASSET_PATH)).toBe(0);
};

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
  const requestCount = installRequestCounter(page);
  await page.goto('/gbdraw/web/index.html');
  await mountExportFixture(page, { interactive: true });

  expect(requestCount(EXPORT_MODULE_PATH)).toBe(0);
  expectNoPdfOrInteractiveRequests(requestCount);

  const downloadPromise = page.waitForEvent('download');
  await page.evaluate(async () => window.__GBDRAW_APP__.downloadInteractiveSVG());
  const download = await downloadPromise;
  expect(download.suggestedFilename()).toBe('lazy-export.interactive.svg');
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

  expect(outcome.result).toEqual({ status: 'error' });
  expect(outcome.summary).toContain('PDF export failed');
  expect(outcome.summary).toContain('Failed to load the vendored jsPDF library');
  await expect(page.getByRole('alert', { name: 'Generation Error' }))
    .toContainText('PDF export failed');
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

  expect(outcome.result).toEqual({ status: 'error' });
  expect(outcome.summary).toContain('SVG export failed');
  expect(outcome.summary).toContain('No SVG result is available for export');
  await expect(page.getByRole('alert', { name: 'Generation Error' }))
    .toContainText('No SVG result is available for export');
});
