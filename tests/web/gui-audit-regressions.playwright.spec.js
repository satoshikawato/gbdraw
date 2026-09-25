const { test, expect } = require('@playwright/test');
const { join } = require('node:path');
const { readFileSync } = require('node:fs');
const { readPdfText } = require('./helpers/pdf-text.cjs');
const { openApp, generateAndWaitForResult, reveal } = require('./helpers/app-lifecycle.cjs');

const session = async (page, name = 'HmmtDNA_basic_circular.gbdraw-session.json') => {
  page.on('dialog', (dialog) => dialog.dismiss());
  await openApp(page);
  await page.locator('input[accept^=".json,"]').setInputFiles(join(process.cwd(), 'gbdraw/web/gallery/sessions', name));
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending && window.__GBDRAW_APP__.extractedFeatures.length > 0);
};
const annotations = (page) => page.locator('details').filter({ has: page.locator('summary[aria-label="Region Annotations"]') });

test('opening saved Circular record controls enables editing before another Generate', async ({ page }) => {
  test.setTimeout(180000);
  await session(page);
  await (await reveal(page.getByLabel('Multi-Record Canvas', { exact: true }))).uncheck();
  await page.locator('summary[aria-label="Circular record presentation"]').click();
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
  await panel.locator('summary').click();
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
  await panel.locator('summary').click();
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
  await features.locator('summary').click();
  await features.locator('select').filter({ has: page.locator('option[value="D-loop"]') }).selectOption('D-loop');
  await features.getByRole('button', { name: /Add$/ }).click();
  await generateAndWaitForResult(page);
  const id = await page.evaluate(() => window.__GBDRAW_APP__.extractedFeatures.find((feature) => feature.type === 'D-loop').svg_id);
  await page.locator(`.origin-top svg [data-gbdraw-feature-id="${id}"]`).first().dispatchEvent('click', { ctrlKey: true });
  const panel = annotations(page);
  await panel.locator('summary').click();
  await panel.getByRole('button', { name: /Add set/ }).click();
  await panel.getByRole('button', { name: /Selected features/ }).click();
  await generateAndWaitForResult(page);
  await expect(page.locator('.origin-top svg [data-gbdraw-annotation-id="feature_1"]')).toHaveCount(1);
});

test('label search includes live edits and reports the initial rendered feature count', async ({ page }) => {
  test.setTimeout(180000);
  await session(page);
  const status = page.getByRole('status', { name: 'Feature search status' });
  await expect(status).toHaveText('0 / 37 features');
  const search = page.getByRole('searchbox', { name: 'Search features', exact: true });
  await search.fill('tRNA');
  await search.press('Enter');
  await page.getByRole('button', { name: 'Open active feature', exact: true }).click();
  await page.locator('.feature-popup input[placeholder="Edit label text"]').fill('AUDIT_RENAMED_FEATURE');
  await page.getByRole('button', { name: 'Apply Label', exact: true }).click();
  if (await page.getByRole('heading', { name: 'Enable Labels', exact: true }).isVisible()) {
    await page.getByRole('button', { name: /Show all labels/ }).click();
  }
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
  await page.locator('select[class*="border-green-300"]').selectOption('1');
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
  await panel.locator('summary').click();
  await panel.getByRole('button', { name: /Add set/ }).click();
  await panel.getByRole('button', { name: /Coordinates/ }).click();
  const before = await page.evaluate(() => JSON.stringify(window.__GBDRAW_APP__.annotationSets));
  const upload = panel.locator('input[type=file]');
  await upload.setInputFiles({ name: 'bad.tsv', mimeType: 'text/plain', buffer: Buffer.from('set_id\tid\tmark\tfeature_selector\ns\ta\thighlight\t;\n') });
  await expect.poll(() => dialogs.length).toBe(1);
  expect(dialogs[0]).toContain('feature_selector requires');
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

test('a failed PDF font request can be retried without reloading the diagram', async ({ page }) => {
  test.setTimeout(180000);
  await session(page);
  await page.route('**/gbdraw-*-py3-none-any.whl*', (route) => route.abort());
  expect(await page.evaluate(() => window.__GBDRAW_APP__.downloadPDF())).toEqual({ status: 'error' });
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.errorLog?.type)).toBe('Export error');
  await page.unroute('**/gbdraw-*-py3-none-any.whl*');
  const pending = page.waitForEvent('download');
  await page.evaluate(() => window.__GBDRAW_APP__.downloadPDF());
  expect(await (await pending).failure()).toBeNull();
  expect(await page.evaluate(() => window.__GBDRAW_APP__.errorLog)).toBeNull();
});

test('a delayed label TSV cannot replace the labels from a newer upload', async ({ page }) => {
  test.setTimeout(180000);
  await session(page);
  const outcome = await page.evaluate(async () => {
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
    return Object.values(app.labelTextFeatureOverrides);
  });
  expect(outcome.length).toBeGreaterThan(0);
  expect(new Set(outcome)).toEqual(new Set(['NEWER']));
});

test('invalid annotation coordinates preserve the draft and the last successful diagram', async ({ page }) => {
  test.setTimeout(180000);
  await session(page);
  const panel = annotations(page);
  await panel.locator('summary').click();
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
  expect((await page.evaluate(() => window.__GBDRAW_APP__.runAnalysis())).status).toBe('error');
  expect(await page.evaluate(() => window.__GBDRAW_APP__.errorLog)).toMatchObject({ summary: expect.stringContaining('positive integers') });
  expect(await page.evaluate(() => window.__GBDRAW_APP__.results[0].content)).toBe(oldResult);
  await panel.getByPlaceholder('Start (1-based)', { exact: true }).fill('1');
  await generateAndWaitForResult(page);
});
