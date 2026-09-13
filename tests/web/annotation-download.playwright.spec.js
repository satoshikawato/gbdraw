const { test, expect } = require('@playwright/test');
const { readFile, writeFile } = require('node:fs/promises');
const { execFileSync } = require('node:child_process');
const { join } = require('node:path');
const { openApp, assertDiagramWorkerIdle, getDiagramWorkerActivity } = require('./helpers/app-lifecycle.cjs');

const panel = (page) => page.locator('details').filter({ has: page.locator('summary[aria-label="Region Annotations"]') });
const effectiveRows = (page) => page.evaluate(() => window.__GBDRAW_APP__.annotationSets.flatMap((set) =>
  set.annotations.map((item) => ({
    setId: set.id, id: item.id, target: item.target, mark: item.mark, label: item.label,
    lane: item.lane, legendLabel: item.legendLabel ?? set.legendLabel,
    style: item.style || set.defaultStyle
  }))
));

const snapshot = (page) => page.evaluate(async () => {
  const app = window.__GBDRAW_APP__;
  const history = window.__GBDRAW_HISTORY__;
  const { getCommittedCanonicalSession } = await import('./js/services/config.js');
  return JSON.stringify({
    draft: app.annotationSets, results: app.results, resultIndex: app.selectedResultIndex,
    annotation: app.selectedAnnotation, features: app.selectedFeatures,
    zoom: app.zoom, pan: app.canvasPan, svg: app.svgContainer?.innerHTML,
    scroll: [app.canvasContainerRef?.scrollLeft, app.canvasContainerRef?.scrollTop],
    committed: getCommittedCanonicalSession(), runInfo: app.lastRunInfo,
    history: [history.getUndoCount(), history.getRedoCount(), history.revision.value,
      history.getCurrentIntentSignature(), history.getCurrentCheckpointSignature(), history.getDiagnostics()]
  });
});

for (const width of [1440, 390]) {
  test(`annotation TSV download/re-import and Python round-trip offline (${width}px)`, async ({ page, browser }, testInfo) => {
    test.setTimeout(180000);
    await page.setViewportSize({ width, height: 960 });
    const external = [];
    const requests = [];
    page.on('request', (request) => requests.push(request.url()));
    await page.context().route('**/*', (route) => {
      const url = new URL(route.request().url());
      if (url.hostname === '127.0.0.1') return route.continue();
      external.push(url.href);
      return route.abort();
    });
    await openApp(page);
    const downloads = [];
    page.on('download', (download) => downloads.push(download));
    await panel(page).locator('summary').click();
    const button = panel(page).getByRole('button', { name: 'Download TSV', exact: true });
    await expect(button).toBeDisabled();
    await expect(button).toHaveAccessibleDescription('Add an annotation row to download TSV.');
    await page.evaluate(() => window.__GBDRAW_APP__.downloadAnnotationTable());
    await panel(page).getByRole('button', { name: /Add set/ }).click();
    await expect(button).toBeDisabled();
    await page.evaluate(() => window.__GBDRAW_APP__.downloadAnnotationTable());
    expect(downloads).toHaveLength(0);

    // A saved Result remains present while its annotation draft is edited.
    await page.locator('input[type=file][accept*="application/json"][accept*="application/gzip"]')
      .setInputFiles(join(process.cwd(), 'gbdraw/web/gallery/sessions/tobacco-chloroplast.gbdraw-session.json'));
    await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending, null, { timeout: 180000 });
    expect(await page.evaluate(() => window.__GBDRAW_APP__.results.length)).toBe(1);
    await page.evaluate(async () => {
      const { normalizeAnnotationSets } = await import('./js/app/annotations/state.js');
      const { coordinateTarget, featureTarget } = await import('./js/app/annotations/target-actions.js');
      const app = window.__GBDRAW_APP__;
      app.annotationSets.splice(0, app.annotationSets.length, ...normalizeAnnotationSets([
        { id: 'empty', annotations: [] },
        { id: 'draft', legendLabel: 'Regions', defaultStyle: { fill: null, stroke: '#123456' }, annotations: [
          { id: 'unique', target: coordinateTarget({ recordId: 'unique', start: 7, end: 21 }), mark: 'band', label: 'Before edit', lane: 1 },
          { id: 'duplicate', target: coordinateTarget({ recordIndex: 1, start: 90, end: 10 }), mark: 'highlight', label: 'Origin span', style: {
            fill: null, strokeWidth: 2, strokeDasharray: [3, 2], lineCap: 'arrow', fillOpacity: 0,
            hatch: { angle: 30, spacing: 6, color: '#654321', width: 2, cross: true },
            labelColor: '#345678', labelFontSize: 14, labelOrientation: 'tangent', labelPosition: 'end', labelOffset: 0
          } }
        ] },
        { id: 'features', annotations: [
          { id: 'gene', target: featureTarget({ recordIndex: 0, selector: 'locus_tag=ABC_1;gene=abc', extent: 'segments', circularPath: 'forward' }), mark: 'bracket', label: 'Gene α', legendLabel: 'Genes' },
          { id: 'local', target: { ...coordinateTarget({ start: 2, end: 8, coordinateSpace: 'local' }), outOfBounds: 'skip' }, mark: 'line' }
        ] }
      ]));
    });
    await panel(page).getByPlaceholder('Label', { exact: true }).first().fill('Edited draft α');
    await panel(page).locator('summary').focus();
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_HISTORY__.capturing.value)).toBe(false);
    const expected = await effectiveRows(page);
    await expect(button).toBeEnabled();
    await button.scrollIntoViewIfNeeded();
    await button.locator('..').screenshot({ path: testInfo.outputPath('annotation-download.png') });
    const bounds = await button.boundingBox();
    expect(bounds.x).toBeGreaterThanOrEqual(0);
    expect(bounds.x + bounds.width).toBeLessThanOrEqual(width);
    // Instrument the shared service's browser URL lifecycle after session loading.
    await page.evaluate(() => {
      window.__annotationUrls = { created: [], revoked: [] };
      const create = URL.createObjectURL.bind(URL);
      const revoke = URL.revokeObjectURL.bind(URL);
      URL.createObjectURL = (blob) => {
        const url = create(blob);
        window.__annotationUrls.created.push({ url, type: blob.type });
        return url;
      };
      URL.revokeObjectURL = (url) => { window.__annotationUrls.revoked.push(url); revoke(url); };
    });
    const before = await snapshot(page);
    const workerBefore = await getDiagramWorkerActivity(page);
    const requestCount = requests.length;
    await page.context().setOffline(true);
    for (let attempt = 0; attempt < 2; attempt += 1) {
      const pending = page.waitForEvent('download');
      await button.click();
      const downloaded = await pending;
      expect(downloaded.suggestedFilename()).toBe('annotations.tsv');
      const path = testInfo.outputPath(`annotations-${attempt}.tsv`);
      await downloaded.saveAs(path);
      expect(await downloaded.failure()).toBeNull();
      expect(await snapshot(page)).toBe(before);
      if (attempt === 0) {
        const expectedPath = testInfo.outputPath('expected-rows.json');
        await writeFile(expectedPath, JSON.stringify(expected, null, 2));
        const pythonResult = execFileSync('python', ['tests/web/helpers/annotation-roundtrip.py', path, expectedPath], { encoding: 'utf8' });
        await writeFile(testInfo.outputPath('python-reader.txt'), pythonResult);
        const context = await browser.newContext({ baseURL: 'http://127.0.0.1:4173', viewport: { width, height: 960 } });
        try {
          await context.route('**/*', (route) => {
            if (new URL(route.request().url()).hostname === '127.0.0.1') return route.continue();
            external.push(route.request().url());
            return route.abort();
          });
          const fresh = await context.newPage();
          await openApp(fresh);
          await context.setOffline(true);
          await panel(fresh).locator('summary').click();
          await panel(fresh).locator('input[type=file]').setInputFiles(path);
          await expect.poll(() => effectiveRows(fresh)).toEqual(expected);
          await assertDiagramWorkerIdle(fresh);
        } finally {
          await context.close();
        }
      } else {
        expect(await readFile(path, 'utf8')).toBe(await readFile(testInfo.outputPath('annotations-0.tsv'), 'utf8'));
      }
    }
    expect(downloads).toHaveLength(2);
    expect(requests).toHaveLength(requestCount);
    expect(external).toEqual([]);
    const urls = await page.evaluate(() => window.__annotationUrls);
    expect(urls.created.map(({ type }) => type)).toEqual(Array(2).fill('text/tab-separated-values;charset=utf-8'));
    expect(urls.revoked).toEqual(urls.created.map(({ url }) => url));
    expect(await getDiagramWorkerActivity(page)).toEqual(workerBefore);
    await writeFile(testInfo.outputPath('browser-evidence.json'), JSON.stringify({ width, external, requests, urls, downloads: downloads.length, unchangedState: true }, null, 2));
  });
}
