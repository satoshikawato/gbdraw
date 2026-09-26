const { test, expect } = require('@playwright/test');
const { readFile, writeFile } = require('node:fs/promises');
const { gunzipSync } = require('node:zlib');
const { join } = require('node:path');
const { openApp, generateAndWaitForResult, getDiagramWorkerActivity } = require('./helpers/app-lifecycle.cjs');

const panel = (page) => page.locator('details').filter({ has: page.locator('summary[aria-label="Region Annotations"]') });
const artifact = (page) => page.evaluate(async () => {
  const { state } = await import("./js/state.js");
  return JSON.parse(JSON.stringify({
  results: window.__GBDRAW_APP__.results,
  warnings: window.__GBDRAW_APP__.annotationWarnings,
  geometry: state.trackSlotResolvedGeometry.value
}));
});

for (const width of [1440, 390]) {
  test(`selector warnings belong to successful Gallery Results and survive History/Session (${width}px)`, async ({ page, browser }, testInfo) => {
    test.setTimeout(240000);
    await page.setViewportSize({ width, height: 960 });
    const consoleText = [];
    page.on('console', (message) => consoleText.push(message.text()));
    page.on('dialog', (dialog) => dialog.accept(dialog.type() === 'prompt' ? 'S02 selectors' : undefined));
    await openApp(page);
    await page.locator('input[type=file][accept*="application/json"][accept*="application/gzip"]')
      .setInputFiles(join(process.cwd(), 'gbdraw/web/gallery/sessions/tobacco-chloroplast.gbdraw-session.json'));
    await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending);
    await generateAndWaitForResult(page);
    const control = await artifact(page);
    expect(control.warnings).toEqual([]);
    const tables = await page.evaluate(async () => {
      const { encodeAnnotationTable } = await import('./js/app/annotations/table-codec.js');
      const { featureTarget } = await import('./js/app/annotations/target-actions.js');
      const sets = JSON.parse(JSON.stringify(window.__GBDRAW_APP__.annotationSets));
      const record = sets[0].annotations[0].target.record;
      sets[0].annotations.push({ ...sets[0].annotations[0], id: 'partial-miss', label: 'SKIPPED-PRIVATE-MARK', legendLabel: 'SKIPPED-PRIVATE-LEGEND',
        target: { ...featureTarget({ selector: 'gene=rbcL;gene=PRIVATE-MISSING;gene=PRIVATE-MISSING' }), record } });
      sets[0].annotations.push({ ...sets[0].annotations[0], id: 'all-miss', label: 'SKIPPED-PRIVATE-MARK', legendLabel: 'SKIPPED-PRIVATE-LEGEND',
        target: { ...featureTarget({ selector: 'gene=PRIVATE-MISSING' }), record } });
      const mixed = encodeAnnotationTable(sets);
      const invalidSets = JSON.parse(JSON.stringify(sets));
      invalidSets[0].annotations[0].target.record = { kind: "recordId", value: "INVALID-RECORD" };
      const invalid = encodeAnnotationTable(invalidSets);
      sets[0].annotations.forEach((row) => { row.target = { ...featureTarget({ selector: 'gene=PRIVATE-MISSING' }), record }; });
      return { mixed, invalid, all: encodeAnnotationTable(sets), count: sets[0].annotations.length };
    });
    await panel(page).locator('summary').focus();
    await page.keyboard.press('Enter');
    const chooser = page.waitForEvent('filechooser');
    await panel(page).getByRole('button', { name: 'Import TSV', exact: true }).focus();
    await page.keyboard.press('Enter');
    await (await chooser).setFiles({ name: 'selectors.tsv', mimeType: 'text/tab-separated-values', buffer: Buffer.from(tables.mixed) });
    await expect.poll(() => panel(page).locator('input[type=file]').inputValue()).toBe('');
    // Exercise the real keyboard Generate command, then wait for its explicit lifecycle settlement.
    await page.getByRole('button', { name: 'Generate Diagram', exact: true }).focus();
    await page.keyboard.press('Enter');
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.processing)).toBe(false);
    const mixed = await artifact(page);
    expect(mixed.warnings).toHaveLength(2);
    expect(mixed.warnings.map((warning) => warning.missingCount)).toEqual([1, 1]);
    expect(mixed.results[0].content).toBe(control.results[0].content);
    expect(mixed.geometry).toEqual(control.geometry);
    const status = page.getByTestId('annotation-resolution-notice');
    await expect(status).toHaveAttribute('role', 'status');
    await expect(status).toHaveAttribute('aria-live', 'polite');
    await expect(status).toHaveAttribute('aria-atomic', 'true');
    await expect(status).toContainText('Skipped 2 annotation row(s).');
    await expect(status).toContainText('plastome_regions/partial-miss');
    await expect(status).toContainText('record #1');
    await status.evaluate((element) => element.scrollIntoView({ block: 'center' }));
    expect(await status.evaluate((element) => element.scrollWidth <= element.clientWidth)).toBe(true);
    const box = await status.boundingBox();
    expect(box.x).toBeGreaterThanOrEqual(0);
    expect(box.x + box.width).toBeLessThanOrEqual(width);
    await page.screenshot({ path: testInfo.outputPath('selector-warning.png') });
    await page.evaluate(() => window.__GBDRAW_HISTORY__.undo());
    expect(await artifact(page)).toEqual(control);
    await expect(status).toHaveCount(0);
    await page.evaluate(() => window.__GBDRAW_HISTORY__.redo());
    expect(await artifact(page)).toEqual(mixed);
    await expect(status).toContainText('Skipped 2 annotation row(s).');

    // A binding failure keeps this exact successful artifact and its notices.
    await panel(page).locator('input[type=file]').setInputFiles({ name: 'invalid-binding.tsv', mimeType: 'text/tab-separated-values',
      buffer: Buffer.from(tables.invalid) });
    await expect.poll(() => panel(page).locator('input[type=file]').inputValue()).toBe('');
    await generateAndWaitForResult(page, { expectedStatus: 'error' });
    expect(await artifact(page)).toEqual(mixed);
    await expect(status).toContainText('Skipped 2 annotation row(s).');

    await panel(page).locator('input[type=file]').setInputFiles({ name: 'all-missing.tsv', mimeType: 'text/tab-separated-values', buffer: Buffer.from(tables.all) });
    await expect.poll(() => panel(page).locator('input[type=file]').inputValue()).toBe('');
    await generateAndWaitForResult(page);
    const all = await artifact(page);
    expect(all.warnings).toHaveLength(tables.count);
    expect(all.results[0].content).not.toContain('data-gbdraw-annotation-id');
    expect(all.results[0].content).not.toContain('SKIPPED-PRIVATE-LEGEND');
    expect(all.geometry).toEqual(mixed.geometry);
    await expect(status).toContainText(`Skipped ${tables.count} annotation row(s).`);
    const pendingSvg = page.waitForEvent('download');
    await page.getByRole('button', { name: 'SVG', exact: true }).click();
    const svgPath = testInfo.outputPath('all-missing.svg');
    await (await pendingSvg).saveAs(svgPath);
    expect(await readFile(svgPath, 'utf8')).not.toContain('data-gbdraw-annotation-id');
    const pendingTsv = page.waitForEvent('download');
    await panel(page).getByRole('button', { name: 'Download TSV', exact: true }).click();
    const tsvPath = testInfo.outputPath('all-missing.tsv');
    await (await pendingTsv).saveAs(tsvPath);
    expect(await readFile(tsvPath, 'utf8')).toContain('PRIVATE-MISSING');
    const pendingSession = page.waitForEvent('download');
    await page.getByRole('button', { name: 'Save Session', exact: true }).click();
    const sessionPath = testInfo.outputPath('selectors.gbdraw-session.json.gz');
    await (await pendingSession).saveAs(sessionPath);
    const session = JSON.parse(gunzipSync(await readFile(sessionPath)));
    expect(session.runMetadata.annotationWarnings).toEqual(all.warnings);
    expect(session.config.annotationSets[0].annotations).toHaveLength(tables.count);
    const context = await browser.newContext({ baseURL: new URL(page.url()).origin, viewport: { width, height: 960 } });
    try {
      const fresh = await context.newPage();
      fresh.on('dialog', (dialog) => dialog.accept());
      await openApp(fresh);
      await fresh.locator('input[type=file][accept*="application/json"][accept*="application/gzip"]').setInputFiles(sessionPath);
      await fresh.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending);
      expect(await artifact(fresh)).toEqual(all);
      expect((await getDiagramWorkerActivity(fresh)).instances.flatMap((instance) => instance.runs)).toEqual([]);
      await expect(fresh.getByTestId('annotation-resolution-notice')).toContainText(`Skipped ${tables.count} annotation row(s).`);
      await generateAndWaitForResult(fresh);
      expect(await artifact(fresh)).toEqual(all);
    } finally { await context.close(); }
    // Web projection still supplies an automatic slot for the unchanged rows,
    // even when Python skips every row. Native automatic slots have a separate contract.
    await page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      state.adv.circular_track_slots_enabled = false;
      state.adv.feature_width_circular = 16;
    });
    await generateAndWaitForResult(page);
    const automatic = await page.evaluate(async () => {
      const { getCommittedCanonicalRenderRequest } = await import('./js/services/config.js');
      const { state } = await import('./js/state.js');
      return { request: getCommittedCanonicalRenderRequest(), geometry: state.trackSlotResolvedGeometry.value,
        warnings: JSON.parse(JSON.stringify(state.annotationWarnings.value)) };
    });
    const slot = automatic.request.diagramOptions.tracks.circularTrackSlots.find((item) => item.id === 'annotations_1');
    expect(slot.params.set_id).toBe('plastome_regions');
    expect(slot.side).toBe('outside');
    expect(automatic.geometry.records[0].slots.some((item) => item.slotId === 'annotations_1')).toBe(true);
    expect(automatic.warnings).toHaveLength(tables.count);
    expect(consoleText.join('\n')).not.toContain('PRIVATE-MISSING');
    await writeFile(testInfo.outputPath('selector-evidence.json'), JSON.stringify({ width, count: tables.count, mixed, all }, null, 2));
  });
}
