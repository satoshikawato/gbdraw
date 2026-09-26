const { test, expect } = require('@playwright/test');
const { readFile, writeFile } = require('node:fs/promises');
const { execFileSync } = require('node:child_process');
const { gunzipSync } = require('node:zlib');
const { join } = require('node:path');
const { openApp, assertDiagramWorkerIdle, getDiagramWorkerActivity, generateAndWaitForResult } = require('./helpers/app-lifecycle.cjs');

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
        const context = await browser.newContext({ baseURL: new URL(page.url()).origin, viewport: { width, height: 960 } });
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

const importSnapshot = async (page) => {
  const value = JSON.parse(await snapshot(page));
  // No-op transactions update diagnostics/revision; entries and signatures must stay unchanged.
  value.history.pop();
  value.history.splice(2, 1);
  return value;
};

for (const width of [1440, 390]) {
  test(`annotation auxiliary file import is atomic, accessible and round-trips (${width}px)`, async ({ page, browser }, testInfo) => {
    test.setTimeout(180000);
    await page.setViewportSize({ width, height: 960 });
    const consoleMessages = [];
    page.on('console', (message) => consoleMessages.push(message.text()));
    const dialogs = [];
    page.on('dialog', async (dialog) => {
      dialogs.push(dialog.message());
      await dialog.accept(dialog.type() === 'prompt' ? 'S01 annotations' : undefined);
    });
    await openApp(page);
    await page.locator('input[type=file][accept*="application/json"][accept*="application/gzip"]')
      .setInputFiles(join(process.cwd(), 'gbdraw/web/gallery/sessions/tobacco-chloroplast.gbdraw-session.json'));
    await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending);
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.results.length)).toBe(1);
    const summary = panel(page).locator('summary');
    await summary.focus();
    await page.keyboard.press('Enter');
    await expect(panel(page)).toHaveAttribute('open', '');
    const fileInput = panel(page).locator('input[type=file]');
    const status = panel(page).getByRole('status');
    await expect(status).toHaveAttribute('aria-live', 'polite');
    await expect(status).toHaveAttribute('aria-atomic', 'true');
    const cases = JSON.parse(await readFile(join(process.cwd(), 'tests/fixtures/annotations/tsv-import-cases.json'), 'utf8'));
    const control = await page.evaluate(async () => {
      const { encodeAnnotationTable } = await import('./js/app/annotations/table-codec.js');
      return encodeAnnotationTable(window.__GBDRAW_APP__.annotationSets);
    });
    const ignored = cases[0].ignored;
    const good = { control, ignored, table: control.trimEnd().split('\n').map((line, index) =>
      `${line}\t${index === 0 ? ignored.join('\t') : 'PRIVATE-CELL\tPRIVATE-GENE\t12345\tred'}`
    ).join('\n') + '\n' };
    const upload = async (text) => {
      await fileInput.setInputFiles({ name: 'annotations.tsv', mimeType: 'text/tab-separated-values', buffer: Buffer.from(text) });
      await expect.poll(() => fileInput.inputValue()).toBe('');
      await expect.poll(() => page.evaluate(() => window.__GBDRAW_HISTORY__.capturing.value)).toBe(false);
    };
    const workerBeforeControl = await getDiagramWorkerActivity(page);
    const pendingChooser = page.waitForEvent('filechooser');
    await panel(page).getByRole('button', { name: 'Import TSV', exact: true }).focus();
    await page.keyboard.press('Enter');
    await (await pendingChooser).setFiles({ name: 'annotations.tsv', mimeType: 'text/tab-separated-values', buffer: Buffer.from(good.control) });
    await expect.poll(() => fileInput.inputValue()).toBe('');
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_HISTORY__.capturing.value)).toBe(false);
    expect(await getDiagramWorkerActivity(page)).toEqual(workerBeforeControl);
    const rows = await effectiveRows(page);
    await generateAndWaitForResult(page);
    const controlSvg = await page.evaluate(() => window.__GBDRAW_APP__.results[0].content);
    const workerBefore = await getDiagramWorkerActivity(page);
    await upload(good.table);
    await expect(status).toContainText('Ignored annotation table columns: notes, gene_desc, pmid, fill_colour.');
    await expect(status).toContainText('not saved in Sessions or TSV re-export');
    expect(await effectiveRows(page)).toEqual(rows);
    expect(await page.evaluate(() => window.__GBDRAW_APP__.results[0].content)).toBe(controlSvg);
    expect(await getDiagramWorkerActivity(page)).toEqual(workerBefore);
    await status.evaluate((element) => element.scrollIntoView({ block: 'center' }));
    const bounds = await status.boundingBox();
    expect(bounds.x).toBeGreaterThanOrEqual(0);
    expect(bounds.x + bounds.width).toBeLessThanOrEqual(width);
    expect(await status.evaluate((element) => element.scrollWidth <= element.clientWidth)).toBe(true);
    await page.screenshot({ path: testInfo.outputPath('annotation-import-notice.png') });
    expect(consoleMessages.join('\n')).not.toContain('PRIVATE-CELL');
    expect(consoleMessages.join('\n')).not.toContain('PRIVATE-GENE');
    const pendingDownload = page.waitForEvent('download');
    await panel(page).getByRole('button', { name: 'Download TSV', exact: true }).click();
    const download = await pendingDownload;
    const tsvPath = testInfo.outputPath('projected-annotations.tsv');
    await download.saveAs(tsvPath);
    const exported = await readFile(tsvPath, 'utf8');
    for (const name of good.ignored) expect(exported.split('\n')[0].split('\t')).not.toContain(name);
    expect(exported).not.toContain('PRIVATE-CELL');
    const expectedPath = testInfo.outputPath('expected-rows.json');
    await writeFile(expectedPath, JSON.stringify(rows));
    const native = execFileSync('python', ['tests/web/helpers/annotation-roundtrip.py', tsvPath, expectedPath], { encoding: 'utf8' });
    await writeFile(testInfo.outputPath('python-reader.txt'), native);
    await generateAndWaitForResult(page);
    expect(await page.evaluate(() => window.__GBDRAW_APP__.results[0].content)).toBe(controlSvg);
    // All malformed inputs must leave the complete prior artifact and draft intact.
    for (const entry of cases.filter((entry) => !entry.valid)) {
      const before = await importSnapshot(page);
      await upload(entry.table);
      expect(await importSnapshot(page), entry.name).toEqual(before);
      await expect(status).toBeEmpty();
    }
    // Exercise actual file-reader rejection and an edit during a delayed read.
    const beforeFailure = await importSnapshot(page);
    await page.evaluate(() => {
      window.__annotationReadOriginal = File.prototype.arrayBuffer;
      File.prototype.arrayBuffer = function () { return Promise.reject(new Error('S01 file-read failure')); };
    });
    await upload(good.table);
    expect(await importSnapshot(page)).toEqual(beforeFailure);
    expect(dialogs.at(-1)).toContain('S01 file-read failure');
    await page.evaluate(() => {
      File.prototype.arrayBuffer = function () { return new Promise((resolve) => { window.__finishAnnotationRead = () => resolve(window.__annotationReadOriginal.call(this)); }); };
    });
    await fileInput.setInputFiles({ name: 'delayed.tsv', mimeType: 'text/tab-separated-values', buffer: Buffer.from(good.table) });
    await page.waitForFunction(() => typeof window.__finishAnnotationRead === 'function');
    await page.evaluate(() => { window.__GBDRAW_APP__.annotationSets[0].annotations[0].label = 'Newer edit'; });
    const resultBeforeStale = await page.evaluate(() => window.__GBDRAW_APP__.results[0].content);
    await page.evaluate(() => window.__finishAnnotationRead());
    await expect.poll(() => fileInput.inputValue()).toBe('');
    expect(await page.evaluate(() => window.__GBDRAW_APP__.annotationSets[0].annotations[0].label)).toBe('Newer edit');
    expect(await page.evaluate(() => window.__GBDRAW_APP__.results[0].content)).toBe(resultBeforeStale);
    await expect(status).toBeEmpty();
    await page.evaluate(() => { File.prototype.arrayBuffer = window.__annotationReadOriginal; });
    await upload(good.table);
    await upload(good.control);
    await expect(status).toBeEmpty();
    const pendingSession = page.waitForEvent('download');
    await page.getByRole('button', { name: 'Save Session', exact: true }).click();
    const sessionDownload = await pendingSession;
    const sessionPath = testInfo.outputPath('annotations.gbdraw-session.json.gz');
    await sessionDownload.saveAs(sessionPath);
    const session = JSON.parse(gunzipSync(await readFile(sessionPath)));
    expect(session.config.annotationSets).toEqual(await page.evaluate(() => JSON.parse(JSON.stringify(window.__GBDRAW_APP__.annotationSets))));
    expect(JSON.stringify(session)).not.toContain('PRIVATE-CELL');
    for (const name of good.ignored) expect(JSON.stringify(session.config.annotationSets)).not.toContain(name);
    expect(JSON.stringify(session)).not.toContain('Ignored annotation table columns:');
    const context = await browser.newContext({ baseURL: new URL(page.url()).origin, viewport: { width, height: 960 } });
    try {
      const fresh = await context.newPage();
      fresh.on('dialog', (dialog) => dialog.accept());
      await openApp(fresh);
      const freshInput = panel(fresh).locator('input[type=file]');
      await freshInput.setInputFiles({ name: 'auxiliary.tsv', mimeType: 'text/tab-separated-values', buffer: Buffer.from(good.table) });
      await expect.poll(() => freshInput.inputValue()).toBe('');
      await expect.poll(() => fresh.evaluate(() => window.__GBDRAW_APP__.annotationImportNotice)).toContain('Ignored annotation table columns:');
      await assertDiagramWorkerIdle(fresh);
      await fresh.locator('input[type=file][accept*="application/json"][accept*="application/gzip"]').setInputFiles(sessionPath);
      await fresh.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending);
      await expect.poll(() => effectiveRows(fresh)).toEqual(rows);
      expect(await fresh.evaluate(() => window.__GBDRAW_APP__.results[0].content)).toBe(controlSvg);
      expect(await fresh.evaluate(() => window.__GBDRAW_APP__.annotationImportNotice)).toBe('');
      expect((await getDiagramWorkerActivity(fresh)).instances.flatMap((instance) => instance.runs)).toEqual([]);
      await generateAndWaitForResult(fresh);
      expect(await fresh.evaluate(() => window.__GBDRAW_APP__.results[0].content)).toBe(controlSvg);
    } finally {
      await context.close();
    }
    await writeFile(testInfo.outputPath('import-evidence.json'), JSON.stringify({ width, cases: cases.length, dialogs, noticesContainNoCells: true, identicalSvg: true, sessionRoundtrip: true }, null, 2));
  });
}
