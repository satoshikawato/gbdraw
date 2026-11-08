const { test, expect } = require('@playwright/test');
const { execFile } = require('node:child_process');
const { promisify } = require('node:util');
const path = require('node:path');
const { openApp, evaluateWithRetainedPromise } = require('./helpers/app-lifecycle.cjs');

// OV-2838: a CLI Session of `--gbk a.gb b.gb`, where a.gb holds two records,
// keys its records `record-1:1`, `record-1:2`, `record-2`. A feature edit on
// b.gb's record stays on that record through a Web Generate: the Web reads one
// File card per file, keyed as the CLI keys it.
const twoRecords = path.resolve('tests/fixtures/web_batch_two_records.gb');
const mito = path.resolve('tests/fixtures/sessions/cli-web-mito.gb');

test('a feature edit on the second file of a CLI Session stays on its record through Generate', async ({ page }, testInfo) => {
  test.setTimeout(900_000);
  const output = testInfo.outputPath('cli.gbdraw-session.json.gz');
  await promisify(execFile)('python', ['-m', 'gbdraw.cli', 'circular', '--gbk', twoRecords, mito,
    '--multi_record_canvas', '-o', testInfo.outputPath('cli'), '--session_output', output],
  { cwd: testInfo.outputDir, env: { ...process.env, PYTHONPATH: process.cwd() }, timeout: 600_000 });
  page.on('dialog', (dialog) => dialog.accept());
  await openApp(page);
  await page.locator('input[accept^=".json,"]').setInputFiles(output);
  await expect.poll(() => page.evaluate(() => !window.__GBDRAW_APP__.sessionImportPending
    && window.__GBDRAW_APP__.extractedFeatures.length > 0), { timeout: 180_000 }).toBe(true);
  const mitoRecordId = await page.evaluate(() => {
    const keys = new Map(window.__GBDRAW_APP__.extractedFeatures.map((feature) => [feature.record_key, feature.record_id]));
    return keys.get('record-2');
  });
  expect(mitoRecordId).toBeTruthy();
  const target = await page.evaluate((recordId) => window.__GBDRAW_APP__.extractedFeatures
    .find((feature) => feature.record_id === recordId && feature.type === 'CDS'), mitoRecordId);
  await evaluateWithRetainedPromise(page, (feature) => window.__GBDRAW_APP__.openFeatureEditorFromList(feature), target);
  await page.locator('.feature-popup input[placeholder="Edit label text"]').fill('OV2838_LABEL');
  await page.getByRole('button', { name: 'Apply Label', exact: true }).click();
  // The default Result draws no label for it: show this one.
  await page.getByRole('button', { name: 'Show this label', exact: true }).click();
  await page.getByRole('button', { name: 'Close feature popup', exact: true }).click();
  expect(await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.runAnalysis())).toEqual({ status: 'ok' });
  const after = await page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    const config = await import('./js/services/config.js');
    const drawing = state.activeDrawing();
    return {
      overrides: Object.keys(drawing.featureOverrides || {}),
      requestKeys: config.getCommittedCanonicalSession().renderRequest.records.map((record) => record.recordKey),
      featureKeys: [...new Set(window.__GBDRAW_APP__.extractedFeatures.map((feature) => `${feature.record_key}=${feature.record_id}`))]
    };
  });
  await require('node:fs/promises').writeFile(testInfo.outputPath('after.json'), JSON.stringify({ mitoRecordId, ...after }, null, 2));
  // The edit is keyed by the record that holds the mito sequence, and the
  // Result draws the edited label.
  const mitoKey = after.featureKeys.find((entry) => entry.endsWith(`=${mitoRecordId}`)).split('=')[0];
  expect(after.overrides.every((key) => JSON.parse(key)[0] === mitoKey)).toBe(true);
  expect(after.overrides).toHaveLength(1);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.currentResults()[0].content)).toContain('OV2838_LABEL');
});
