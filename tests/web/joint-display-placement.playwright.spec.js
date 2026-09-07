const { test, expect } = require('@playwright/test');
const fs = require('node:fs/promises');
const zlib = require('node:zlib');
const path = require('node:path');
const { spawnSync } = require('node:child_process');
const { openApp } = require('./helpers/app-lifecycle.cjs');

const replayEnv = { ...process.env };
delete replayEnv.PYTHONPATH;
delete replayEnv.PYTHONHOME;

test.beforeEach(async ({ page }) => {
  page.on('dialog', async (dialog) => { console.log('Session dialog:', dialog.message()); await dialog.dismiss(); });
  page.on('console', (message) => { if (message.type() === 'error') console.log('Browser error:', message.text()); });
});

const generateFromControl = async (page) => {
  const before = await page.evaluate(async () => (await import('./js/state.js')).state.resultGenerationKey.value);
  await page.getByRole('button', { name: 'Generate Diagram', exact: true }).click();
  await expect.poll(() => page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    return { generation: state.resultGenerationKey.value, processing: state.processing.value, error: state.errorLog.value?.summary || '' };
  }), { timeout: 180000 }).toMatchObject({ generation: before + 1, processing: false, error: '' });
};

const svgStructure = (page, content) => page.evaluate((svg) => {
  const tree = new DOMParser().parseFromString(svg, 'image/svg+xml');
  const visit = (node) => node.nodeType === Node.ELEMENT_NODE
    ? [node.tagName, Object.fromEntries([...node.attributes].map((attr) => [attr.name, attr.value]).sort()), [...node.childNodes].map(visit)]
    : node.textContent;
  return visit(tree.documentElement);
}, content);

const selectPaintedFeature = async (page, body) => {
  await body.scrollIntoViewIfNeeded();
  const point = await body.evaluate((element) => {
    const box = element.getBBox();
    for (let y = 1; y < 30; y += 1) for (let x = 1; x < 30; x += 1) {
      const local = new DOMPoint(box.x + box.width * x / 30, box.y + box.height * y / 30);
      if (!element.isPointInFill(local)) continue;
      const screen = local.matrixTransform(element.getScreenCTM());
      if (document.elementFromPoint(screen.x, screen.y)?.getAttribute('data-gbdraw-feature-id') === element.getAttribute('data-gbdraw-feature-id')) return { x: screen.x, y: screen.y };
    }
    throw new Error('No visible painted feature point.');
  });
  await page.keyboard.down('Control');
  await page.mouse.click(point.x, point.y);
  await page.keyboard.up('Control');
};

const source = `LOCUS       shared                   360 bp    DNA     circular UNA 01-JAN-2000
DEFINITION  Joint surface regression.
ACCESSION   shared
VERSION     shared
KEYWORDS    .
SOURCE      synthetic construct
  ORGANISM  synthetic construct
FEATURES             Location/Qualifiers
     CDS             21..105
                     /locus_tag="chosen"
                     /product="chosen protein"
     CDS             81..150
                     /locus_tag="neighbor"
                     /product="neighbor protein"
ORIGIN
        1 ${'acgt'.repeat(90)}
//
`;

for (const mode of ['circular', 'linear']) {
  test(`@pr-smoke joint rotation placement and saved drafts in ${mode}`, async ({ page, browser }, testInfo) => {
    test.setTimeout(240000);
    await page.addInitScript(() => {
      window.__JOINT_LIFECYCLE__ = [];
      window.__GBDRAW_TEST_HOOKS__ = { onSessionLifecycleEvent: (event) => window.__JOINT_LIFECYCLE__.push(event) };
    });
    await openApp(page);
    if (mode === 'linear') await page.getByRole('button', { name: 'Linear', exact: true }).click();
    const upload = mode === 'linear' ? page.getByTestId('linear-genbank-1')
      : page.getByLabel('GenBank/DDBJ File', { exact: true });
    await upload.setInputFiles({ name: 'shared.gbk', mimeType: 'text/plain', buffer: Buffer.from(source) });
    if (mode === 'linear') await page.getByRole('button', { name: 'Record options for sequence 1', exact: true }).click();
    const start = page.getByRole('spinbutton', { name: 'Display start shared #1', exact: true });
    await expect(start).toBeEnabled();
    await generateFromControl(page);
    const initial = await page.evaluate(() => ({
      result: window.__GBDRAW_APP__.results[0].content,
      feature: window.__GBDRAW_APP__.extractedFeatures[0]
    }));
    expect(initial.result).toContain('<svg');
    await start.fill('71');
    await start.press('Tab');
    await page.locator('.drawer-toggle').click();
    await page.locator('.right-drawer').getByRole('button', { name: 'Edit', exact: true }).first().click();
    const placement = page.getByRole('combobox', { name: 'Feature placement', exact: true });
    await expect(placement).toBeVisible();
    await expect(placement.locator('option[value=main]')).toHaveJSProperty('disabled', false);
    await placement.selectOption('main');
    await expect.poll(() => page.evaluate(() => { const app = window.__GBDRAW_APP__; return app.featurePlacementActions.valueFor(app.extractedFeatures[0]); })).toBe('main');
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.results[0].content)).toBe(initial.result);
    await page.getByRole('button', { name: 'Close feature popup', exact: true }).click();
    await page.locator('.drawer-toggle').click();
    const tolerance = page.getByRole('spinbutton', { name: 'Feature overlap tolerance (bp)', exact: true });
    await tolerance.fill('1');
    await tolerance.press('Tab');
    await generateFromControl(page);
    const reuse = await page.evaluate(() => ({
      python: window.__JOINT_LIFECYCLE__.filter((event) => event.name === 'python-diagnostics').at(-1),
      worker: window.__GBDRAW_DIAGRAM_WORKER_ACTIVITY__
    }));
    expect(reuse.python.metrics).toMatchObject({ parsedSourceCacheHitCount: 1, parsedSourceParseCount: 0,
      preparedInputCacheMutationViolationCount: 0 });
    expect(reuse.worker.constructions).toBe(1);
    expect(reuse.worker.instances[0].runs).toHaveLength(2);
    await fs.writeFile(testInfo.outputPath('joint-cache-reuse.json'), JSON.stringify(reuse, null, 2));
    await expect.poll(() => page.evaluate(() => { const app = window.__GBDRAW_APP__; return app.featurePlacementActions.valueFor(app.extractedFeatures[0]); })).toBe('main');
    const runInfo = await page.evaluate(() => window.__GBDRAW_APP__.lastRunInfo);
    expect(runInfo.sourceRecipe.available, runInfo.sourceRecipe.unavailableReason).toBe(true);
    await page.getByRole('button', { name: /Run info/i }).click();
    const bundlePromise = page.waitForEvent('download');
    await page.getByRole('button', { name: /Download reproducibility files/ }).click();
    const bundle = await bundlePromise;
    const bundlePath = testInfo.outputPath(bundle.suggestedFilename());
    await bundle.saveAs(bundlePath);
    const sourcePath = testInfo.outputPath('shared.gbk');
    const infoPath = testInfo.outputPath('run-info.json');
    await fs.writeFile(sourcePath, source);
    await fs.writeFile(infoPath, JSON.stringify(runInfo));
    const replay = spawnSync('python', [path.resolve('tests/web/helpers/joint-replay.py'), bundlePath,
      sourcePath, infoPath, testInfo.outputPath('replay')], { encoding: 'utf8', env: replayEnv });
    expect(replay.status, replay.stdout + replay.stderr).toBe(0);
    await start.fill('91');
    await start.press('Tab');
    await page.getByRole('button', { name: /^Undo/ }).first().click();
    await expect(start).toHaveValue('71');
    await page.getByRole('button', { name: /^Redo/ }).first().click();
    await expect(start).toHaveValue('91');
    const savedResult = await page.evaluate(async () => {
      const { serializeResults } = await import('./js/services/config.js');
      return serializeResults()[0].content;
    });
    await page.evaluate(() => { window.__GBDRAW_APP__.sessionTitle = 'Joint surfaces'; });
    const downloadPromise = page.waitForEvent('download');
    await page.getByRole('button', { name: 'Save Session', exact: true }).click();
    const download = await downloadPromise;
    const savedPath = testInfo.outputPath(download.suggestedFilename());
    await download.saveAs(savedPath);
    const bytes = await fs.readFile(savedPath);
    const session = JSON.parse((bytes[0] === 0x1f ? zlib.gunzipSync(bytes) : bytes).toString());
    expect(session.version).toBe(41);
    expect(session.renderRequest.schema).toBe(7);
    expect(session.renderRequest.records[0].display.startCoordinate).toBe(71);
    expect(session.renderRequest.diagramOptions.featurePlacements).toHaveLength(1);
    expect(session.renderRequest.diagramOptions.configOverrides['canvas.feature_overlap_tolerance_bp']).toBe(1);
    expect(session.config.recordDisplayDrafts.some((draft) => draft.startCoordinate === 91)).toBe(true);
    expect(session.results[0].content).toBe(savedResult);
    await page.screenshot({ path: testInfo.outputPath('joint-controls.png'), fullPage: true });
    const context = await browser.newContext();
    const loaded = await context.newPage();
    await openApp(loaded);
    await loaded.locator('input[accept^=".json,"]').setInputFiles(savedPath);
    await expect.poll(() => loaded.evaluate(() => window.__GBDRAW_APP__.results.length)).toBe(1);
    const loadedResult = await loaded.evaluate(() => window.__GBDRAW_APP__.results[0].content);
    expect(await svgStructure(loaded, loadedResult)).toEqual(await svgStructure(loaded, savedResult));
    await loaded.locator('.drawer-toggle').click();
    await loaded.locator('.right-drawer').getByRole('button', { name: 'Edit', exact: true }).first().click();
    const loadedPlacement = loaded.getByRole('combobox', { name: 'Feature placement', exact: true });
    // Resolved geometry is regenerated; Auto removal only needs the saved source binding.
    await expect(loadedPlacement.locator('option[value=main]')).toHaveJSProperty('disabled', true);
    await expect(loadedPlacement.locator('option[value=auto]')).toHaveJSProperty('disabled', false);
    await loadedPlacement.selectOption('auto');
    await loaded.getByRole('button', { name: 'Close feature popup', exact: true }).click();
    await loaded.locator('.drawer-toggle').click();
    await loaded.getByRole('button', { name: /^Undo/ }).first().click();
    await expect.poll(() => loaded.evaluate(() => { const app = window.__GBDRAW_APP__; return app.featurePlacementActions.valueFor(app.extractedFeatures[0]); })).toBe('main');
    await generateFromControl(loaded);
    expect(await loaded.evaluate(() => window.__GBDRAW_APP__.results[0].content)).not.toBe(savedResult);
    await context.close();
  });
}

test('historical v40/schema6 saves as the joint format without Generate', async ({ page }, testInfo) => {
  test.setTimeout(180000);
  await openApp(page);
  await page.locator('input[accept^=".json,"]').setInputFiles(
    'tests/fixtures/sessions/rendered-v27.v40-schema6.json.gz');
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.sessionImportPending), { timeout: 180000 }).toBe(false);
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.results.length)).toBe(1);
  await page.evaluate(() => { window.__GBDRAW_APP__.sessionTitle = 'Historical promotion'; });
  const pending = page.waitForEvent('download');
  await page.getByRole('button', { name: 'Save Session', exact: true }).click();
  const download = await pending;
  const savedPath = testInfo.outputPath(download.suggestedFilename());
  await download.saveAs(savedPath);
  const bytes = await fs.readFile(savedPath);
  const saved = JSON.parse(zlib.gunzipSync(bytes));
  expect(saved.version).toBe(41);
  expect(saved.renderRequest.schema).toBe(7);
  expect(saved.renderRequest.records.every((record) => record.display.isCircular === null
    && record.display.startCoordinate === null)).toBe(true);
  expect(saved.renderRequest.diagramOptions.featurePlacements).toEqual([]);
  // This CLI fixture validates preserved config through the existing helper;
  // promotion must never submit a diagram generation request.
  expect(await page.evaluate(() => window.__GBDRAW_DIAGRAM_WORKER_ACTIVITY__.instances
    .reduce((count, instance) => count + instance.runs.length, 0))).toBe(0);
});

for (const mode of ['circular', 'linear']) {
  test(`${mode} duplicate source records keep independent starts through replacement and Undo`, async ({ page }, testInfo) => {
    test.setTimeout(180000);
    await openApp(page);
    if (mode === 'linear') await page.getByRole('button', { name: 'Linear', exact: true }).click();
    const upload = mode === 'linear' ? page.getByTestId('linear-genbank-1')
      : page.getByLabel('GenBank/DDBJ File', { exact: true });
    const bytes = Buffer.from(source + source);
    await upload.setInputFiles({ name: 'same.gbk', mimeType: 'text/plain', buffer: bytes });
    if (mode === 'linear') await page.getByRole('button', { name: 'Record options for sequence 1', exact: true }).click();
    const first = page.getByRole('spinbutton', { name: 'Display start shared #1', exact: true });
    const second = page.getByRole('spinbutton', { name: 'Display start shared #2', exact: true });
    await expect(second).toBeEnabled();
    await first.fill('71'); await first.press('Tab');
    await second.fill('91'); await second.press('Tab');
    await generateFromControl(page);
    const committed = await page.evaluate(async () => (await import('./js/services/config.js')).getCommittedCanonicalRenderRequest());
    expect(committed.records.map((record) => record.display.startCoordinate)).toEqual([71, 91]);
    expect(new Set(committed.records.map((record) => record.recordKey)).size).toBe(2);
    expect(new Set(committed.records.map((record) => record.source.resourceId)).size).toBe(1);
    if (mode === 'linear') expect(committed.records.map((record) => record.presentation.gridRow)).toEqual([1, 1]);
    await upload.setInputFiles({ name: 'same.gbk', mimeType: 'text/plain', buffer: Buffer.from((source + source).replace('chosen protein', 'replacement protein')) });
    await expect(first).toHaveValue('');
    await expect(second).toHaveValue('');
    await page.getByRole('button', { name: /^Undo/ }).first().click();
    await expect(first).toHaveValue('71');
    await expect(second).toHaveValue('91');
    await page.evaluate(() => { window.__GBDRAW_APP__.sessionTitle = 'Duplicate sources'; });
    const pendingSave = page.waitForEvent('download');
    await page.getByRole('button', { name: 'Save Session', exact: true }).click();
    const saved = await pendingSave;
    const savedPath = testInfo.outputPath(saved.suggestedFilename());
    await saved.saveAs(savedPath);
    await page.locator('input[accept^=".json,"]').setInputFiles(savedPath);
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.sessionImportPending), { timeout: 180000 }).toBe(false);
    if (mode === 'linear') {
      for (const button of await page.getByRole('button', { name: /Record options for sequence/ }).all()) {
        if (!await button.evaluate((element) => element.parentElement.open)) await button.click();
      }
    }
    await page.getByRole('button', { name: 'Load record rotation controls', exact: true }).first().click();
    await expect(page.getByRole('spinbutton', { name: 'Display start shared #1', exact: true }).first()).toHaveValue('71');
    await expect(page.getByRole('spinbutton', { name: 'Display start shared #2', exact: true }).last()).toHaveValue('91');
    await generateFromControl(page);
  });
}

for (const mode of ['circular', 'linear']) {
  test(`${mode} real shortcuts topology reset bulk placement and Auto history`, async ({ page }, testInfo) => {
    test.setTimeout(180000);
    await openApp(page);
    if (mode === 'linear') await page.getByRole('button', { name: 'Linear', exact: true }).click();
    const upload = mode === 'linear' ? page.getByTestId('linear-genbank-1') : page.getByLabel('GenBank/DDBJ File', { exact: true });
    await upload.setInputFiles({ name: 'shared.gbk', mimeType: 'text/plain', buffer: Buffer.from(source) });
    if (mode === 'linear') await page.getByRole('button', { name: 'Record options for sequence 1', exact: true }).click();
    const start = page.getByRole('spinbutton', { name: 'Display start shared #1', exact: true });
    const topology = page.getByRole('checkbox', { name: 'Circular record shared #1', exact: true });
    await expect(start).toBeEnabled();
    await start.fill('71'); await start.press('Tab');
    await topology.uncheck();
    await expect(start).toBeDisabled();
    await expect(start).toHaveValue('71');
    await page.evaluate(() => { window.__GBDRAW_APP__.sessionTitle = 'Inactive draft'; });
    const inactiveSave = page.waitForEvent('download');
    await page.getByRole('button', { name: 'Save Session', exact: true }).click();
    const inactiveDownload = await inactiveSave;
    const inactivePath = testInfo.outputPath(inactiveDownload.suggestedFilename());
    await inactiveDownload.saveAs(inactivePath);
    const inactive = JSON.parse(zlib.gunzipSync(await fs.readFile(inactivePath)));
    expect(inactive.results).toEqual([]);
    expect(inactive.config.recordDisplayDrafts[0]).toMatchObject({ topologyOverride: false, startCoordinate: 71 });
    await page.getByRole('button', { name: 'Reset to detected', exact: true }).click();
    await expect(start).toBeEnabled();
    await expect(start).toHaveValue('71');
    await page.getByRole('button', { name: 'Reset start', exact: true }).click();
    await expect(start).toHaveValue('');
    if (mode === 'circular') await page.locator('#circular-track-preset').selectOption('middle');
    await page.getByRole('checkbox', { name: 'Separate Strands', exact: true }).uncheck();
    await generateFromControl(page);
    const features = await page.evaluate(() => window.__GBDRAW_APP__.extractedFeatures.map((f) => ({ id: f.svg_id, parts: f.location_parts, strand: f.strand })));
    const bodies = features.map((feature) => page.locator(`[data-gbdraw-feature-id="${feature.id}"]`).first());
    await selectPaintedFeature(page, bodies[0]);
    await expect(page.getByLabel('Selected feature count')).toContainText('1 selected');
    await page.getByRole('button', { name: 'Use selected feature 5′ end', exact: true }).click();
    await expect(start).toHaveValue('21');
    await page.getByRole('button', { name: 'Use selected feature midpoint', exact: true }).click();
    await expect(start).toHaveValue('63');
    await selectPaintedFeature(page, bodies[1]);
    await expect(page.getByLabel('Selected feature count')).toContainText('2 selected');
    await expect(page.getByRole('button', { name: 'Use selected feature midpoint', exact: true })).toBeDisabled();
    const bulk = page.getByRole('combobox', { name: 'Selected feature placements', exact: true });
    await bulk.selectOption(mode === 'circular' ? 'outward' : 'above');
    await expect.poll(() => page.evaluate(async () => Object.keys((await import('./js/state.js')).state.featurePlacementOverrides).length)).toBe(2);
    await bulk.selectOption('auto');
    await expect.poll(() => page.evaluate(async () => Object.keys((await import('./js/state.js')).state.featurePlacementOverrides).length)).toBe(0);
    await page.getByRole('button', { name: /^Undo/ }).first().click();
    await expect.poll(() => page.evaluate(async () => Object.keys((await import('./js/state.js')).state.featurePlacementOverrides).length)).toBe(2);
    await page.getByRole('button', { name: /^Redo/ }).first().click();
    await expect.poll(() => page.evaluate(async () => Object.keys((await import('./js/state.js')).state.featurePlacementOverrides).length)).toBe(0);
    await generateFromControl(page);
    expect(await page.evaluate(() => window.__GBDRAW_APP__.extractedFeatures.map((f) => ({ id: f.svg_id, parts: f.location_parts, strand: f.strand })))).toEqual(features);
    await page.screenshot({ path: testInfo.outputPath('shortcuts-bulk.png'), fullPage: true });
  });
}

test('paired GFF/FASTA replacement purges source-bound rotation and placement and Undo restores both', async ({ page }) => {
  test.setTimeout(180000);
  await openApp(page);
  await page.getByRole('radio', { name: 'GFF3 + FASTA', exact: true }).check();
  const gff = '##gff-version 3\n##sequence-region shared 1 360\nshared\t.\tCDS\t21\t105\t.\t+\t0\tID=chosen.full.source;product=chosen protein\n';
  const fasta = '>shared\n' + 'ACGT'.repeat(90) + '\n';
  await page.getByLabel('GFF3 File', { exact: true }).setInputFiles({ name: 'same.gff', mimeType: 'text/plain', buffer: Buffer.from(gff) });
  const upload = page.getByLabel('FASTA File', { exact: true });
  await upload.setInputFiles({ name: 'same.fa', mimeType: 'text/plain', buffer: Buffer.from(fasta) });
  const start = page.getByRole('spinbutton', { name: 'Display start shared #1', exact: true });
  await page.getByRole('checkbox', { name: 'Circular record shared #1', exact: true }).check();
  await start.fill('71'); await start.press('Tab');
  await generateFromControl(page);
  await page.locator('.drawer-toggle').click();
  await page.locator('.right-drawer').getByRole('button', { name: 'Edit', exact: true }).first().click();
  await page.getByRole('combobox', { name: 'Feature placement', exact: true }).selectOption('main');
  await page.getByRole('button', { name: 'Close feature popup', exact: true }).click();
  await page.locator('.drawer-toggle').click();
  await generateFromControl(page);
  const before = await page.evaluate(() => window.__GBDRAW_APP__.results[0].content);
  await upload.setInputFiles({ name: 'same.fa', mimeType: 'text/plain', buffer: Buffer.from(fasta.replace('ACGT', 'TCGT')) });
  await expect(start).toHaveValue('');
  await expect.poll(() => page.evaluate(async () => Object.keys((await import('./js/state.js')).state.featurePlacementOverrides).length)).toBe(0);
  await page.getByRole('button', { name: /^Undo/ }).first().click();
  await expect(start).toHaveValue('71');
  await expect.poll(() => page.evaluate(async () => Object.keys((await import('./js/state.js')).state.featurePlacementOverrides).length)).toBe(1);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.results[0].content)).toBe(before);
});

for (const mode of ['circular', 'linear']) for (const intent of ['rotation', 'placement']) {
  test(`${mode} ${intent} alone exports working Source recipe and Exact replay`, async ({ page }, testInfo) => {
    test.setTimeout(180000);
    await openApp(page);
    if (mode === 'linear') await page.getByRole('button', { name: 'Linear', exact: true }).click();
    const upload = mode === 'linear' ? page.getByTestId('linear-genbank-1') : page.getByLabel('GenBank/DDBJ File', { exact: true });
    await upload.setInputFiles({ name: 'shared.gbk', mimeType: 'text/plain', buffer: Buffer.from(source) });
    if (intent === 'rotation') {
      if (mode === 'linear') await page.getByRole('button', { name: 'Record options for sequence 1', exact: true }).click();
      const start = page.getByRole('spinbutton', { name: 'Display start shared #1', exact: true });
      await expect(start).toBeEnabled(); await start.fill('71'); await start.press('Tab');
    } else {
      await generateFromControl(page);
      await page.locator('.drawer-toggle').click();
      await page.locator('.right-drawer').getByRole('button', { name: 'Edit', exact: true }).first().click();
      await page.getByRole('combobox', { name: 'Feature placement', exact: true }).selectOption('main');
      await page.getByRole('button', { name: 'Close feature popup', exact: true }).click();
      await page.locator('.drawer-toggle').click();
      const tolerance = page.getByRole('spinbutton', { name: 'Feature overlap tolerance (bp)', exact: true });
      await tolerance.fill('1'); await tolerance.press('Tab');
    }
    await generateFromControl(page);
    const runInfo = await page.evaluate(() => window.__GBDRAW_APP__.lastRunInfo);
    expect(runInfo.sourceRecipe.available, runInfo.sourceRecipe.unavailableReason).toBe(true);
    await page.getByRole('button', { name: /Run info/i }).click();
    const pending = page.waitForEvent('download');
    await page.getByRole('button', { name: /Download reproducibility files/ }).click();
    const download = await pending; const bundlePath = testInfo.outputPath(download.suggestedFilename());
    await download.saveAs(bundlePath);
    const sourcePath = testInfo.outputPath('shared.gbk'), infoPath = testInfo.outputPath('run-info.json');
    await fs.writeFile(sourcePath, source); await fs.writeFile(infoPath, JSON.stringify(runInfo));
    const expected = intent === 'rotation' ? { start: 71, placements: 0, tolerance: 0 } : { start: null, placements: 1, tolerance: 1 };
    const result = spawnSync('python', [path.resolve('tests/web/helpers/joint-replay.py'), bundlePath, sourcePath, infoPath,
      testInfo.outputPath('replay'), JSON.stringify(expected)], { encoding: 'utf8', env: replayEnv });
    expect(result.status, result.stdout + result.stderr).toBe(0);
  });
}

test('loaded replacement draft cannot place or take shortcuts from the prior Result source', async ({ page }, testInfo) => {
  test.setTimeout(180000);
  await openApp(page);
  const upload = page.getByLabel('GenBank/DDBJ File', { exact: true });
  await upload.setInputFiles({ name: 'shared.gbk', mimeType: 'text/plain', buffer: Buffer.from(source) });
  await generateFromControl(page);
  await upload.setInputFiles({ name: 'shared.gbk', mimeType: 'text/plain', buffer: Buffer.from(source.replace('chosen protein', 'replacement protein')) });
  await page.evaluate(() => { window.__GBDRAW_APP__.sessionTitle = 'Replacement draft'; });
  const pending = page.waitForEvent('download');
  await page.getByRole('button', { name: 'Save Session', exact: true }).click();
  const download = await pending, savedPath = testInfo.outputPath(download.suggestedFilename());
  await download.saveAs(savedPath);
  await page.locator('input[accept^=".json,"]').setInputFiles(savedPath);
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.sessionImportPending), { timeout: 180000 }).toBe(false);
  await page.getByRole('button', { name: 'Load record rotation controls', exact: true }).click();
  await expect(page.getByRole('spinbutton', { name: 'Display start shared #1', exact: true })).toBeEnabled();
  await page.locator('.drawer-toggle').click();
  await page.locator('.right-drawer').getByRole('button', { name: 'Edit', exact: true }).first().click();
  const placement = page.getByRole('combobox', { name: 'Feature placement', exact: true });
  await expect(placement.locator('option[value=main]')).toHaveJSProperty('disabled', true);
  await page.getByRole('button', { name: 'Close feature popup', exact: true }).click();
  await page.locator('.drawer-toggle').click();
  await generateFromControl(page);
  await page.locator('.drawer-toggle').click();
  await page.locator('.right-drawer').getByRole('button', { name: 'Edit', exact: true }).first().click();
  await expect(placement.locator('option[value=main]')).toHaveJSProperty('disabled', false);
});
