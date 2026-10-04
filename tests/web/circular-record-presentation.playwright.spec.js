const { test, expect } = require('@playwright/test');
const { readFileSync, writeFileSync } = require('node:fs');
const { join, resolve } = require('node:path');
const {
  evaluateWithRetainedPromise,
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
  writeFileSync(testInfo.outputPath('circular-record-presentation.svg'), cropReverse.content);
  const renderedSvg = page.locator('main svg').last();
  await expect(renderedSvg).toBeVisible();
  await renderedSvg.screenshot({
    path: testInfo.outputPath('circular-record-presentation.png')
  });

  const downloadPromise = page.waitForEvent('download', { timeout: 180000 });
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.saveSessionWithTitle());
  const savedSessionPath = await (await downloadPromise).path();
  expect(savedSessionPath).toBeTruthy();

  await page.getByLabel('Circular region end').fill('60000');
  const failed = await generateAndWaitForResult(page, {
    expectedStatus: 'error',
    requireCommittedResult: true
  });
  expect(failed.errorSummary).toContain(
    'Keep the region within the record length.'
  );
  expect(await page.evaluate(() => window.__GBDRAW_APP__.errorLog)).toMatchObject({
    code: 'REGION_INVALID', operation: 'generate', stage: 'request-validation',
    context: { field: 'region', reason: 'RECORD_BOUNDS' }
  });
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

  await expect(presentationDetails).toHaveAttribute('open', '');
  await page.getByLabel('Circular record label').fill('');
  await page.getByLabel('Circular record subtitle').fill('');
  await generateAndWaitForResult(page);
  const inferred = await inspectCircularResult(page);
  expect(inferred.text).toContain('Streptomyces fradiae');
  expect(inferred.text).not.toContain('Selected second record');
});

const singleSource = readFileSync(join(repoRoot, 'tests/fixtures/sessions/cli-web-mito.gb'), 'utf8');
const panel = (page) => page.locator('[data-circular-record-presentation]');
const ready = (page) => expect.poll(() => page.evaluate(async () => (
  (await import('./js/state.js')).state.circularRecordDiscovery.status
))).toBe('ready');
const nativeUpload = (page, content, name = 'presentation.gb') => page.getByLabel('GenBank/DDBJ File', { exact: true })
  .setInputFiles({ name, mimeType: 'text/plain', buffer: Buffer.from(content) });

for (const width of [1280, 390]) {
  test(`applicable disclosure keeps manual close, focus, scroll and keyboard at ${width}px`, async ({ page }, testInfo) => {
    await page.setViewportSize({ width, height: 900 });
    await openApp(page);
    await page.evaluate(() => { window.__GBDRAW_APP__.form.multi_record_canvas = false; });
    const upload = page.getByRole('button', { name: 'Choose GenBank/DDBJ File', exact: true });
    await upload.focus();
    const uploadAnchor = await upload.boundingBox();
    await nativeUpload(page, singleSource);
    await ready(page);
    await expect(upload).toBeFocused();
    await expect.poll(async () => Math.abs((await upload.boundingBox()).y - uploadAnchor.y)).toBeLessThanOrEqual(2);
    await expect(panel(page)).toHaveAttribute('open', '');
    await expect(page.getByLabel('Circular region start')).toBeEnabled();
    const baseline = await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
    const summary = panel(page).locator('summary');
    await summary.focus();
    await page.keyboard.press('Enter');
    await expect(panel(page)).not.toHaveAttribute('open', '');
    await expect(summary).toBeFocused();
    await page.evaluate(() => { window.__GBDRAW_APP__.form.plot_title = 'Unrelated update'; });
    await expect(panel(page)).not.toHaveAttribute('open', '');
    expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(baseline);
    await page.keyboard.press('Space');
    await expect(panel(page)).toHaveAttribute('open', '');
    await page.keyboard.press('Enter');
    await nativeUpload(page, `${firstGenbank}\n${secondGenbank}`, 'two.gb');
    await ready(page);
    await expect(panel(page)).not.toHaveAttribute('open', '');
    const selector = page.getByLabel('Circular record', { exact: true });
    await selector.focus();
    const before = await selector.boundingBox();
    await selector.selectOption('BGC0000709');
    await expect(panel(page)).toHaveAttribute('open', '');
    await expect(selector).toBeFocused();
    await expect.poll(async () => Math.abs((await selector.boundingBox()).y - before.y))
      .toBeLessThanOrEqual(2);
    expect(await page.evaluate(() => window.__GBDRAW_APP__.form.circular_record_selector)).toBe('BGC0000709');
    expect(await page.evaluate(() => window.__GBDRAW_APP__.form.multi_record_canvas)).toBe(false);
    await summary.focus();
    await page.keyboard.press('Enter');
    await page.evaluate(() => { window.__GBDRAW_APP__.adv.def_font_size += 1; });
    await expect(panel(page)).not.toHaveAttribute('open', '');
    await selector.selectOption('BGC0000708');
    await expect(panel(page)).toHaveAttribute('open', '');
    if (width === 390) {
      expect(await page.evaluate(() => document.documentElement.scrollWidth)).toBeLessThanOrEqual(width);
    }
    await page.locator('.settings-pane').screenshot({ path: testInfo.outputPath(`disclosure-${width}.png`) });
  });
}

test('grid and explicit batch explain applicability outside disclosure and require an explicit one-record choice', async ({ page }) => {
  await openApp(page);
  await page.evaluate(() => { window.__GBDRAW_APP__.form.multi_record_canvas = false; });
  await nativeUpload(page, `${firstGenbank}\n${secondGenbank}`);
  await ready(page);
  await expect(panel(page)).not.toHaveAttribute('open', '');
  await expect(page.getByText('Select one record to edit crop, orientation, and title lines.', { exact: true })).toBeVisible();
  await page.evaluate(() => { window.__GBDRAW_APP__.form.multi_record_canvas = true; });
  await expect(page.getByLabel('Circular record', { exact: true })).toBeDisabled();
  const action = await page.evaluate(() => window.__GBDRAW_APP__.setCircularRecordPresentationSelector('BGC0000708'));
  expect(action.status).toBe('unavailable');
  expect(await page.evaluate(() => window.__GBDRAW_APP__.form.circular_record_selector)).toBe('');
  await expect(page.getByText('Multi-Record Canvas uses a grid. Turn it off to select one record.', { exact: true })).toBeVisible();
  await page.getByRole('button', { name: 'Show Multi-Record Canvas setting' }).click();
  const grid = page.locator('[data-circular-canvas-setting]');
  await expect(grid).toBeFocused();
  await expect(grid).toBeChecked();
  await page.keyboard.press('Space');
  await page.getByLabel('Circular record', { exact: true }).selectOption('BGC0000709');
  await expect(panel(page)).toHaveAttribute('open', '');
  await expect(page.getByLabel('Circular region start')).toBeEnabled();
  // A saved one-record batch remains a batch until the user selects that record.
  await page.evaluate(() => {
    window.__GBDRAW_APP__.form.circular_record_selector = '';
    window.__GBDRAW_APP__.adv.circular_grouping_intent = 'batch';
  });
  await nativeUpload(page, singleSource);
  await ready(page);
  await expect(page.getByLabel('Circular region start')).toBeDisabled();
  await expect(page.getByText('Select one record to edit crop, orientation, and title lines.', { exact: true })).toBeVisible();
  expect(await page.evaluate(() => window.__GBDRAW_APP__.adv.circular_grouping_intent)).toBe('batch');
  await page.getByLabel('Circular record', { exact: true }).selectOption('NC_012920.1');
  await expect(page.getByLabel('Circular region start')).toBeEnabled();
});

test('History restore re-inspects sources and reconciles missing or ambiguous selectors without changing crop intent', async ({ page }) => {
  await openApp(page);
  await page.evaluate(() => { window.__GBDRAW_APP__.form.multi_record_canvas = false; });
  await nativeUpload(page, singleSource);
  await ready(page);
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(1);
  await page.getByRole('button', { name: 'Undo', exact: true }).click();
  await expect.poll(() => page.evaluate(async () => (await import('./js/state.js')).state.circularRecordDiscovery.status)).toBe('idle');
  await page.getByRole('button', { name: 'Redo', exact: true }).click();
  await ready(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.circularRecordList.length)).toBe(1);
  await page.getByLabel('Circular region start').fill('25');
  await page.getByLabel('Circular region end').fill('100');
  await page.getByLabel('Circular region end').blur();
  await page.evaluate(() => { window.__GBDRAW_APP__.form.circular_record_selector = 'missing'; });
  await expect(page.getByLabel('Circular region start')).toBeDisabled();
  await expect(page.getByText("Record selector 'missing' was not found in the current input.", { exact: true })).toBeVisible();
  await nativeUpload(page, `${singleSource}\n${singleSource}`);
  await ready(page);
  await page.evaluate(() => { window.__GBDRAW_APP__.form.circular_record_selector = 'NC_012920.1'; });
  await expect(page.getByLabel('Circular region start')).toBeDisabled();
  await expect(page.getByText("Record selector 'NC_012920.1' is ambiguous in the current input.", { exact: true })).toBeVisible();
  await page.getByLabel('Circular record', { exact: true }).selectOption('#2');
  await expect(page.getByLabel('Circular region start')).toBeEnabled();
  await expect(page.getByLabel('Circular region start')).toHaveValue('25');
  await expect(page.getByLabel('Circular region end')).toHaveValue('100');
});
