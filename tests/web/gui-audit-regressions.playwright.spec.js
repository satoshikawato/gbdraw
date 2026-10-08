const { test, expect } = require('@playwright/test');
const { join } = require('node:path');
const { readFileSync } = require('node:fs');
const { readPdfText } = require('./helpers/pdf-text.cjs');
const { openApp, generateAndWaitForResult, reveal, evaluateWithRetainedPromise } = require('./helpers/app-lifecycle.cjs');

const session = async (page, name = 'HmmtDNA_basic_circular.gbdraw-session.json') => {
  page.on('dialog', (dialog) => dialog.dismiss());
  await openApp(page);
  await page.locator('input[accept^=".json,"]').setInputFiles(join(process.cwd(), 'gbdraw/web/gallery/sessions', name));
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending && window.__GBDRAW_APP__.extractedFeatures.length > 0);
};
// Inspection may already open the presentation panel; open it only when closed.
const openPresentation = async (page) => {
  const details = page.locator('[data-circular-record-presentation]');
  if (!await details.evaluate((element) => element.open)) await details.locator('summary').click();
};
const annotations = (page) => page.locator('details').filter({ has: page.locator('summary[aria-label="Region Annotations"]') });

test('opening saved Circular record controls enables editing before another Generate', async ({ page }) => {
  test.setTimeout(180000);
  await session(page);
  await (await reveal(page.getByLabel('Multi-Record Canvas', { exact: true }))).uncheck();
  await page.locator('[data-circular-inspect]').click();
  await openPresentation(page);
  await expect(page.getByLabel('Circular record label', { exact: true })).toBeEnabled({ timeout: 180000 });
  await expect(page.getByLabel('Circular record', { exact: true })).toContainText('NC_012920.1');
  await page.getByLabel('Circular record label', { exact: true }).fill('Edited before Generate');
  await generateAndWaitForResult(page);
  await expect(page.locator('.origin-top svg')).toContainText('Edited before Generate');
});

test('new annotation IDs still select their editor row after deletion and regeneration', async ({ page }) => {
  test.setTimeout(180000);
  await session(page);
  const panel = annotations(page);
  await panel.locator('summary').press('Enter');
  await panel.getByRole('button', { name: /Add set/ }).click();
  const add = panel.getByRole('button', { name: /Coordinates/ });
  for (let i = 0; i < 3; i += 1) await add.click();
  await panel.getByTitle('Delete annotation', { exact: true }).first().click();
  await add.click();
  for (let i = 0; i < 3; i += 1) {
    await panel.getByPlaceholder('Start (1-based)', { exact: true }).nth(i).fill(String(1000 + i * 3000));
    await panel.getByPlaceholder('End (inclusive)', { exact: true }).nth(i).fill(String(2000 + i * 3000));
  }
  await panel.getByTitle('Fill color', { exact: true }).first().fill('#ff0000');
  await expect(panel.getByTitle('Fill color', { exact: true }).nth(1)).toHaveValue('#94a3b8');
  await generateAndWaitForResult(page);
  const ids = await page.evaluate(() => window.__GBDRAW_APP__.annotationSets[0].annotations.map((item) => item.id));
  expect(new Set(ids).size).toBe(3);
  await expect(page.locator(`.origin-top svg [data-gbdraw-annotation-id="${ids[0]}"] path`).first()).toHaveAttribute('fill', '#ff0000');
  await expect(page.locator(`.origin-top svg [data-gbdraw-annotation-id="${ids[1]}"] path`).first()).toHaveAttribute('fill', '#94a3b8');
  const target = page.locator(`.origin-top svg [data-gbdraw-annotation-id="${ids[2]}"]`).first();
  await target.dispatchEvent('click');
  await expect(panel.locator('.border-blue-400')).toHaveCount(1);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.selectedAnnotation.id)).toBe(ids[2]);
});

test('imported annotation colors edit the rendered row style', async ({ page }) => {
  test.setTimeout(180000);
  await session(page);
  const panel = annotations(page);
  await panel.locator('summary').press('Enter');
  await panel.locator('input[type=file]').setInputFiles({ name: 'regions.tsv', mimeType: 'text/plain', buffer: Buffer.from('set_id\tid\tmark\tstart\tend\tstroke\tfill\nregions\ta\tband\t1000\t3000\t#123456\t#94a3b8\n') });
  await panel.getByTitle('Fill color', { exact: true }).fill('#ff0000');
  await panel.getByTitle('Stroke color', { exact: true }).fill('#00ff00');
  await generateAndWaitForResult(page);
  const path = page.locator('.origin-top svg [data-gbdraw-annotation-id="a"] path').first();
  await expect(path).toHaveAttribute('fill', '#ff0000');
  await expect(path).toHaveAttribute('stroke', '#00ff00');
});

test('selected D-loop without gene or locus_tag generates a feature annotation', async ({ page }) => {
  test.setTimeout(180000);
  await session(page);
  const features = page.locator('details').filter({ has: page.locator('summary[aria-label="Features"]') });
  await features.locator('summary').press('Enter');
  await features.locator('select').filter({ has: page.locator('option[value="D-loop"]') }).selectOption('D-loop');
  await features.getByRole('button', { name: /Add$/ }).click();
  await generateAndWaitForResult(page);
  const id = await page.evaluate(() => window.__GBDRAW_APP__.extractedFeatures.find((feature) => feature.type === 'D-loop').svg_id);
  await page.locator(`.origin-top svg [data-gbdraw-feature-id="${id}"]`).first().dispatchEvent('click', { ctrlKey: true });
  const panel = annotations(page);
  await panel.locator('summary').press('Enter');
  await panel.getByRole('button', { name: /Add set/ }).click();
  await panel.getByRole('button', { name: /Selected features/ }).click();
  await generateAndWaitForResult(page);
  await expect(page.locator('.origin-top svg [data-gbdraw-annotation-id="feature_1"]')).toHaveCount(1);
});

test('label search includes live edits and reports the initial rendered feature count', async ({ page }) => {
  test.setTimeout(180000);
  await session(page);
  await page.getByRole('button', { name: 'Toggle layout edit mode' }).click();
  await expect(page.getByRole('button', { name: 'Toggle layout edit mode' })).toHaveAttribute('aria-pressed', 'true');
  const status = page.getByRole('status', { name: 'Feature search status' });
  await expect(status).toHaveText('0 / 37 features');
  const search = page.getByRole('searchbox', { name: 'Search features', exact: true });
  await search.fill('tRNA');
  await search.press('Enter');
  await page.getByRole('button', { name: 'Open active feature', exact: true }).click();
  await page.locator('.feature-popup input[placeholder="Edit label text"]').fill('AUDIT_RENAMED_FEATURE');
  await page.getByRole('button', { name: 'Apply Label', exact: true }).click();
  await page.waitForFunction(() => !window.__GBDRAW_APP__.labelReflowProcessing && !window.__GBDRAW_APP__.processing);
  await page.getByRole('button', { name: 'Close feature popup', exact: true }).click();
  await page.getByRole('combobox', { name: 'Search field', exact: true }).selectOption('label');
  await search.fill('AUDIT_RENAMED_FEATURE');
  await search.press('Enter');
  await expect(status).toHaveText('1 / 1 features');
  await generateAndWaitForResult(page);
  await expect(status).toHaveText('1 / 1 features');
});

test('similarity group names and descriptions survive temporary mode changes', async ({ page }) => {
  test.setTimeout(180000);
  await session(page, 'majanivirus_orthogroup.gbdraw-session.json.gz');
  await page.locator('.drawer-toggle').click();
  const drawer = page.locator('.right-drawer');
  await drawer.getByRole('button', { name: 'Similarity groups' }).click();
  await page.getByPlaceholder('Rename this similarity group').fill('AUDIT_GROUP');
  await drawer.locator('textarea:visible').fill('AUDIT_DESCRIPTION');
  await drawer.locator('textarea:visible').press('Tab');
  await page.locator('.drawer-toggle').click();
  await page.getByRole('combobox', { name: 'Search field', exact: true }).selectOption('orthogroup');
  const search = page.getByRole('searchbox', { name: 'Search features', exact: true });
  for (const query of ['AUDIT_GROUP', 'AUDIT_DESCRIPTION']) {
    await search.fill(query);
    await search.press('Enter');
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.previewFeatureSearchMatches.length)).toBeGreaterThan(0);
  }
  const before = await page.evaluate(() => ({ names: { ...window.__GBDRAW_APP__.orthogroupNameOverrides }, descriptions: { ...window.__GBDRAW_APP__.orthogroupDescriptionOverrides }, count: window.__GBDRAW_APP__.orthogroupCount }));
  await page.getByRole('button', { name: 'Circular', exact: true }).click();
  await page.getByRole('button', { name: 'Linear', exact: true }).click();
  const inspect = () => page.evaluate(() => ({ names: { ...window.__GBDRAW_APP__.orthogroupNameOverrides }, descriptions: { ...window.__GBDRAW_APP__.orthogroupDescriptionOverrides }, count: window.__GBDRAW_APP__.orthogroupCount }));
  await expect.poll(inspect).toEqual(before);
  await generateAndWaitForResult(page);
  await expect.poll(inspect).toEqual(before);
});

for (const { method, delayedAsset, extension } of [
  { method: 'downloadPDF', delayedAsset: 'vendor/jspdf/jspdf.umd.min.js', extension: 'pdf' },
  { method: 'downloadPDF', delayedAsset: 'js/services/export.js', extension: 'pdf' },
  { method: 'downloadInteractiveSVG', delayedAsset: 'js/services/standalone-interactivity.js', extension: 'interactive.svg' }
]) test(`${method} keeps the selected batch output while loading ${delayedAsset}`, async ({ page }) => {
  test.setTimeout(180000);
  await openApp(page);
  const source = ['BGC0000708', 'BGC0000709'].map((id) => readFileSync(join(process.cwd(), 'tests/test_inputs', `${id}.gbk`), 'utf8')).join('\n');
  await page.getByLabel('GenBank/DDBJ File', { exact: true }).setInputFiles({ name: 'two-records.gbk', mimeType: 'text/plain', buffer: Buffer.from(source) });
  await page.waitForFunction(() => window.__GBDRAW_APP__.circularRecordList.length === 2);
  await (await reveal(page.getByLabel('Multi-Record Canvas', { exact: true }))).uncheck();
  await generateAndWaitForResult(page);
  let release;
  let requested;
  const delayed = new Promise((resolve) => { release = resolve; });
  const requestStarted = new Promise((resolve) => { requested = resolve; });
  await page.route(`**/${delayedAsset}`, async (route) => { requested(); await delayed; await route.continue(); });
  const filename = await page.evaluate((extension) => window.__GBDRAW_APP__.results[0].name.replace(/\.svg$/, `.${extension}`), extension);
  await page.evaluate((method) => { window.__AUDIT_EXPORT__ = window.__GBDRAW_APP__[method](); }, method);
  await requestStarted;
  await page.locator('[aria-label="Result Preview"] h2 select').selectOption('1');
  const pending = page.waitForEvent('download');
  release();
  await page.evaluate(() => window.__AUDIT_EXPORT__);
  const download = await pending;
  expect(download.suggestedFilename()).toBe(filename);
  const bytes = readFileSync(await download.path());
  const content = extension === 'pdf' ? readPdfText(bytes) : bytes.toString('utf8');
  expect(content).toContain('BGC0000708');
  expect(content).not.toContain('BGC0000709');
});

test('malformed annotations preserve the draft and a valid import is undoable', async ({ page }) => {
  await openApp(page);
  const errors = [];
  const dialogs = [];
  page.on('pageerror', (error) => errors.push(error.message));
  page.on('dialog', (dialog) => { dialogs.push(dialog.message()); dialog.dismiss(); });
  const panel = annotations(page);
  await panel.locator('summary').press('Enter');
  await panel.getByRole('button', { name: /Add set/ }).click();
  await panel.getByRole('button', { name: /Coordinates/ }).click();
  const before = await page.evaluate(() => JSON.stringify(window.__GBDRAW_APP__.annotationSets));
  const upload = panel.locator('input[type=file]');
  await upload.setInputFiles({ name: 'bad.tsv', mimeType: 'text/plain', buffer: Buffer.from('set_id\tid\tmark\tfeature_selector\ns\ta\thighlight\t;\n') });
  // FL-14: the panel notice names the error; no browser alert.
  await expect(panel.getByRole('status')).toContainText('feature_selector requires');
  expect(dialogs).toEqual([]);
  expect(await page.evaluate(() => JSON.stringify(window.__GBDRAW_APP__.annotationSets))).toBe(before);
  await upload.setInputFiles({ name: 'good.tsv', mimeType: 'text/plain', buffer: Buffer.from('set_id\tid\tmark\tstart\tend\nloaded\ta\tband\t1\t8\n') });
  await expect(panel.getByLabel('Annotation set id')).toHaveValue('loaded');
  await page.getByRole('button', { name: 'Undo', exact: true }).click();
  await expect.poll(() => page.evaluate(() => JSON.stringify(window.__GBDRAW_APP__.annotationSets))).toBe(before);
  await page.getByRole('button', { name: 'Redo', exact: true }).click();
  await expect(panel.getByLabel('Annotation set id')).toHaveValue('loaded');
  expect(errors).toEqual([]);
});

test('PDF preserves Greek and mathematical text from the existing bundled fonts', async ({ page }) => {
  test.setTimeout(180000);
  await session(page);
  const title = 'β-lactamase α ≥ 95%';
  await page.locator('#circular-species').fill(title);
  await page.locator('#circular-species').press('Tab');
  await generateAndWaitForResult(page);
  const pending = page.waitForEvent('download');
  await page.evaluate(() => window.__GBDRAW_APP__.downloadPDF());
  const bytes = readFileSync(await (await pending).path());
  expect(readPdfText(bytes)).toContain(title);
  expect(bytes.toString('latin1')).toContain('/FontFile2');
});

test('PDF pages convert CSS px to pt and keep spaces in curved and tick labels', async ({ page }) => {
  test.setTimeout(180000);
  await session(page);
  const [width, height] = await page.evaluate(() => {
    const svg = window.__GBDRAW_APP__.svgContainer.querySelector('svg');
    return [parseFloat(svg.getAttribute('width')), parseFloat(svg.getAttribute('height'))];
  });
  const pending = page.waitForEvent('download');
  await page.evaluate(() => window.__GBDRAW_APP__.downloadPDF());
  const bytes = readFileSync(await (await pending).path());
  const mediaBoxes = [...bytes.toString('latin1').matchAll(/MediaBox\s*\[([^\]]*)\]/g)];
  expect(mediaBoxes).toHaveLength(1);
  const [, , boxWidth, boxHeight] = mediaBoxes[0][1].trim().split(/\s+/).map(Number);
  expect(boxWidth).toBeCloseTo(width * 0.75, 1);
  expect(boxHeight).toBeCloseTo(height * 0.75, 1);
  const text = readPdfText(bytes);
  expect(text).toContain('cytochrome c oxidase subunit I');
  expect(text).toContain('1 kbp');
});

test('PDF text flattening places each curved-label character by its UTF-16 index', async ({ page }) => {
  await page.goto('/');
  const outcome = await page.evaluate(async () => {
    const markup = '<svg xmlns="http://www.w3.org/2000/svg" width="400px" height="200px" viewBox="0 0 400 200">'
      + '<defs><path id="arc" d="M 10 150 A 190 190 0 0 1 390 150"/></defs>'
      + '<text font-size="20"><textPath href="#arc">a\u{1F600} b</textPath></text></svg>';
    const live = new DOMParser().parseFromString(markup, 'image/svg+xml').documentElement;
    document.body.appendChild(live);
    const text = live.querySelector('text');
    const expected = [[0, 'a'], [1, '\u{1F600}'], [3, ' '], [4, 'b']].map(([index, character]) => {
      const point = text.getStartPositionOfChar(index);
      return [character, Math.round(point.x * 100) / 100, Math.round(point.y * 100) / 100, !character.trim()];
    });
    const flattened = [];
    const observer = new MutationObserver((records) => records.forEach((record) => {
      record.addedNodes.forEach((node) => {
        if (node.localName !== 'g') return;
        node.querySelectorAll('text').forEach((charText) => flattened.push([
          charText.textContent,
          Math.round(Number(charText.getAttribute('x')) * 100) / 100,
          Math.round(Number(charText.getAttribute('y')) * 100) / 100,
          charText.getAttributeNS('http://www.w3.org/XML/1998/namespace', 'space') === 'preserve'
        ]));
      });
    }));
    observer.observe(document.body, { childList: true, subtree: true });
    const { downloadPDF } = await import('/gbdraw/web/js/services/export.js');
    try {
      await downloadPDF(
        { svg: live.cloneNode(true), name: 'flatten.svg' },
        { loadPdfFont: () => Promise.reject(new Error('font not needed')) }
      );
    } catch {
      // The flattened text is captured before fonts are loaded.
    }
    await new Promise((resolve) => setTimeout(resolve, 0));
    observer.disconnect();
    return { expected, flattened };
  });
  expect(outcome.flattened).toEqual(outcome.expected);
});

test('a failed PDF font request can be retried without reloading the diagram', async ({ page }) => {
  test.setTimeout(180000);
  await session(page);
  await page.route('**/gbdraw-*-py3-none-any.whl*', (route) => route.abort());
  expect(await page.evaluate(() => window.__GBDRAW_APP__.downloadPDF())).toMatchObject({
    status: 'error', error: { code: 'WORKER_INIT', operation: 'export-pdf', stage: 'initialization' }
  });
  await expect(page.getByRole('alert', { name: 'PDF export error' })).toBeVisible();
  expect(await page.evaluate(() => window.__GBDRAW_APP__.errorLog)).toMatchObject({
    operation: 'export-pdf', stage: 'initialization'
  });
  await page.unroute('**/gbdraw-*-py3-none-any.whl*');
  const pending = page.waitForEvent('download');
  await page.evaluate(() => window.__GBDRAW_APP__.downloadPDF());
  expect(await (await pending).failure()).toBeNull();
  expect(await page.evaluate(() => window.__GBDRAW_APP__.errorLog)).toBeNull();
});

test('a delayed label TSV cannot replace the labels from a newer upload', async ({ page }) => {
  test.setTimeout(180000);
  await session(page);
  const outcome = await evaluateWithRetainedPromise(page, async () => {
    const app = window.__GBDRAW_APP__;
    const older = new File(['*\t*\tproduct\t.*\tOLDER\n'], 'older.tsv');
    const newer = new File(['*\t*\tproduct\t.*\tNEWER\n'], 'newer.tsv');
    const bytes = await older.arrayBuffer();
    let release;
    let signalStarted;
    const started = new Promise((resolve) => { signalStarted = resolve; });
    const held = new Promise((resolve) => { release = resolve; });
    older.arrayBuffer = async () => { signalStarted(); await held; return bytes; };
    const input = { files: [older], value: 'older.tsv' };
    const pending = app.loadLabelOverrideTable({ target: input });
    await started;
    input.files = [newer];
    input.value = 'newer.tsv';
    await app.loadLabelOverrideTable({ target: input });
    release();
    await pending;
    return Object.values(app.featureOverrides).map((row) => row.labelText).filter((text) => text !== null);
  });
  expect(outcome.length).toBeGreaterThan(0);
  expect(new Set(outcome)).toEqual(new Set(['NEWER']));
});

// OV-131: a label TSV replaces the label edits only when it applies to a label.
// One that matches no label keeps them, records no Undo step, and says so.
test('a label TSV that matches no label keeps the label edits and records no Undo step', async ({ page }) => {
  test.setTimeout(180000);
  await session(page);
  const alerts = [];
  page.on('dialog', (dialog) => alerts.push(dialog.message()));
  const importTsv = (tsv) => evaluateWithRetainedPromise(page, async (text) => {
    const app = window.__GBDRAW_APP__;
    await app.loadLabelOverrideTable({ target: { files: [new File([text], 'labels.tsv')], value: 'labels.tsv' } });
    await new Promise((resolve) => requestAnimationFrame(() => requestAnimationFrame(resolve)));
    const history = window.__GBDRAW_HISTORY__;
    return { bulk: { ...app.labelTextBulkOverrides }, undoCount: history.getUndoCount(), undoLabel: history.undoLabel() };
  }, tsv);
  const edited = await importTsv('*\t*\tlabel\t^tRNA-Phe$\tMyLabel\n');
  expect(edited).toMatchObject({ bulk: { 'tRNA-Phe': 'MyLabel' }, undoLabel: 'Load label edits' });
  expect(await importTsv('*\t*\tproduct\t^no such product$\tZZ\n')).toEqual(edited);
  const replaced = await importTsv('*\t*\tlabel\t^tRNA-Val$\tOther\n');
  expect(replaced).toEqual({ bulk: { 'tRNA-Val': 'Other' }, undoCount: edited.undoCount + 1, undoLabel: 'Load label edits' });
  expect(alerts).toEqual([
    'Loaded 1 row(s). Applied to 1 label(s).',
    'Loaded 1 row(s). Not applied: no row matched a label of the diagram. The existing label edits were kept.',
    'Loaded 1 row(s). Applied to 1 label(s).'
  ]);
});

test('invalid annotation coordinates preserve the draft and the last successful diagram', async ({ page }) => {
  test.setTimeout(180000);
  await session(page);
  const panel = annotations(page);
  await panel.locator('summary').press('Enter');
  await panel.getByRole('button', { name: /Add set/ }).click();
  await panel.getByRole('button', { name: /Coordinates/ }).click();
  const before = await page.evaluate(() => JSON.stringify(window.__GBDRAW_APP__.annotationSets));
  const dialogs = [];
  page.on('dialog', (dialog) => dialogs.push(dialog.message()));
  for (const [index, start] of ['abc', '0', '-10', '1.5'].entries()) {
    await panel.locator('input[type=file]').setInputFiles({ name: 'bad.tsv', mimeType: 'text/plain', buffer: Buffer.from(`set_id\tid\tmark\tstart\tend\ns\ta\thighlight\t${start}\t3000\n`) });
    await expect.poll(() => dialogs.length).toBe(index + 1);
    expect(dialogs.at(-1)).toContain('positive integers');
    expect(await page.evaluate(() => JSON.stringify(window.__GBDRAW_APP__.annotationSets))).toBe(before);
  }
  const oldResult = await page.evaluate(() => window.__GBDRAW_APP__.results[0].content);
  await panel.getByPlaceholder('Start (1-based)', { exact: true }).fill('1.5');
  await panel.getByPlaceholder('Start (1-based)', { exact: true }).press('Tab');
  expect((await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.runAnalysis())).status).toBe('error');
  expect(await page.evaluate(() => window.__GBDRAW_APP__.errorLog)).toMatchObject({ code: 'ANNOTATION_TARGET', context: { reason: 'POSITIVE_INTEGER' },
    summary: expect.stringContaining('integer greater than zero') });
  expect(await page.evaluate(() => window.__GBDRAW_APP__.results[0].content)).toBe(oldResult);
  await panel.getByPlaceholder('Start (1-based)', { exact: true }).fill('1');
  await generateAndWaitForResult(page);
});
