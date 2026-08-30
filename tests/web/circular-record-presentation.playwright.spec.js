const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { join, resolve } = require('node:path');
const {
  generateAndWaitForResult,
  openApp,
  waitForAppShell
} = require('./helpers/app-lifecycle.cjs');

const repoRoot = resolve(process.env.GBDRAW_REPO || process.cwd());
const firstGenbank = readFileSync(
  join(repoRoot, 'tests/test_inputs/BGC0000708.gbk'),
  'utf8'
);
const secondGenbank = readFileSync(
  join(repoRoot, 'tests/test_inputs/BGC0000709.gbk'),
  'utf8'
);

const inspectCircularResult = (page) => page.evaluate(() => {
  const app = window.__GBDRAW_APP__;
  const content = String(app.results?.[0]?.content || '');
  const svg = new DOMParser().parseFromString(content, 'image/svg+xml').documentElement;
  return {
    content,
    text: String(svg.textContent || '').replace(/\s+/g, ' ').trim(),
    recordIds: [...svg.querySelectorAll('[data-gbdraw-record-id]')]
      .map((element) => element.getAttribute('data-gbdraw-record-id'))
      .filter(Boolean),
    canonicalRecord: structuredClone(
      window.__GBDRAW_CIRCULAR_REQUESTS__?.at(-1)?.records?.[0] || null
    )
  };
});

test('Circular single-record presentation selects, transforms, titles, and round-trips one record', async ({ page }, testInfo) => {
  test.setTimeout(420000);
  await page.addInitScript(() => {
    window.__GBDRAW_CIRCULAR_REQUESTS__ = [];
    const nativePostMessage = Worker.prototype.postMessage;
    Worker.prototype.postMessage = function trackedCircularRequest(message, transfer) {
      if (message?.type === 'run' && message.payload?.request?.mode === 'circular') {
        window.__GBDRAW_CIRCULAR_REQUESTS__.push(
          structuredClone(message.payload.request)
        );
      }
      if (transfer === undefined) return nativePostMessage.call(this, message);
      return nativePostMessage.call(this, message, transfer);
    };
  });
  await openApp(page, { waitForPalette: false });

  await page.evaluate(({ first, second }) => {
    const app = window.__GBDRAW_APP__;
    app.mode = 'circular';
    app.cInputType = 'gb';
    app.files.c_gb = new File([`${first}\n${second}`], 'two-records.gbk', {
      type: 'text/plain',
      lastModified: 1
    });
    Object.assign(app.form, {
      multi_record_canvas: false,
      suppress_gc: true,
      suppress_skew: true,
      labels_mode: 'none',
      legend: 'none'
    });
    app.sessionTitle = 'circular-record-presentation';
  }, { first: firstGenbank, second: secondGenbank });

  await page.waitForFunction(() => (
    window.__GBDRAW_APP__?.circularRecordList?.length === 2
  ), null, { timeout: 180000 });
  const presentationDetails = page.locator('[data-circular-record-presentation]');
  await presentationDetails.locator('summary').click();
  const recordSelect = page.getByLabel('Circular record', { exact: true });
  await expect(recordSelect).toContainText('BGC0000708');
  await expect(recordSelect).toContainText('BGC0000709');
  await recordSelect.selectOption('BGC0000709');
  await page.getByLabel('Circular record label').fill('長い 環状ゲノム ラベル');
  await page.getByLabel('Circular record subtitle').fill('Selected second record');

  await page.getByLabel('Circular reverse complement').check();
  await generateAndWaitForResult(page);
  const reverseOnly = await inspectCircularResult(page);
  expect(new Set(reverseOnly.recordIds)).toEqual(new Set(['BGC0000709']));
  expect(reverseOnly.text).toContain('長い 環状ゲノム ラベル');
  expect(reverseOnly.text).toContain('Selected second record');
  expect(reverseOnly.canonicalRecord.selector).toEqual({
    kind: 'recordId',
    value: 'BGC0000709'
  });
  expect(reverseOnly.canonicalRecord.region).toBeNull();
  expect(reverseOnly.canonicalRecord.presentation.reverseComplement).toBe(true);

  await page.getByLabel('Circular region start').fill('1000');
  await page.getByLabel('Circular region end').fill('12000');
  await generateAndWaitForResult(page);
  const cropReverse = await inspectCircularResult(page);
  expect(cropReverse.text).toContain('11,001 bp');
  expect(cropReverse.canonicalRecord.selector).toBeNull();
  expect(cropReverse.canonicalRecord.region).toEqual({
    selector: { kind: 'recordId', value: 'BGC0000709' },
    start: 1000,
    end: 12000,
    reverseComplement: true
  });
  expect(cropReverse.canonicalRecord.presentation.reverseComplement).toBe(false);
  const renderedSvg = page.locator('main svg').last();
  await expect(renderedSvg).toBeVisible();
  await renderedSvg.screenshot({
    path: testInfo.outputPath('circular-record-presentation.png')
  });

  const downloadPromise = page.waitForEvent('download', { timeout: 180000 });
  await page.evaluate(() => window.__GBDRAW_APP__.saveSessionWithTitle());
  const savedSessionPath = await (await downloadPromise).path();
  expect(savedSessionPath).toBeTruthy();

  await page.getByLabel('Circular region end').fill('60000');
  const failed = await generateAndWaitForResult(page, {
    expectedStatus: 'error',
    requireCommittedResult: true
  });
  expect(failed.errorSummary).toContain(
    'Circular region End (60000) exceeds the selected record length (50466).'
  );
  expect((await inspectCircularResult(page)).content).toBe(cropReverse.content);

  await page.reload({ waitUntil: 'domcontentloaded' });
  await waitForAppShell(page, { waitForPalette: false });
  const dialogPromise = page.waitForEvent('dialog', { timeout: 180000 });
  await page.locator('input[accept^=".json,"]').first().setInputFiles(savedSessionPath);
  const dialog = await dialogPromise;
  expect(dialog.message()).toBe('Session loaded successfully!');
  await dialog.accept();
  await expect(page.getByLabel('Circular record', { exact: true })).toHaveValue('BGC0000709');
  await expect(page.getByLabel('Circular region start')).toHaveValue('1000');
  await expect(page.getByLabel('Circular region end')).toHaveValue('12000');
  await expect(page.getByLabel('Circular reverse complement')).toBeChecked();

  await generateAndWaitForResult(page);
  const regenerated = await inspectCircularResult(page);
  expect(regenerated.content).toBe(cropReverse.content);
  expect(regenerated.canonicalRecord.region.reverseComplement).toBe(true);
  expect(regenerated.canonicalRecord.presentation.reverseComplement).toBe(false);
  expect(new Set(regenerated.recordIds)).toEqual(new Set(['BGC0000709']));

  await presentationDetails.locator('summary').click();
  await page.getByLabel('Circular record label').fill('');
  await page.getByLabel('Circular record subtitle').fill('');
  await generateAndWaitForResult(page);
  const inferred = await inspectCircularResult(page);
  expect(inferred.text).toContain('Streptomyces fradiae');
  expect(inferred.text).not.toContain('Selected second record');
});
