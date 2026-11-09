const { test, expect } = require('@playwright/test');
const { readFileSync, writeFileSync } = require('node:fs');
const { join, resolve } = require('node:path');
const {
  evaluateWithRetainedPromise,
  generateAndWaitForResult,
  openApp,
  waitForAppShell,
  expectSessionLoadedWithoutWarnings,
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
  const content = String(app.currentResults()[0]?.content || '');
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

  await page.evaluate(async ({ first, second }) => {
    const app = window.__GBDRAW_APP__;
    app.selectDrawing('circular');
    app.setCircularInputType('gb');
    // The File transition runs the record discovery the record list waits for (R10).
    void app.addCircularSourceFiles('gb', [new File([`${first}\n${second}`], 'two-records.gbk', {
      type: 'text/plain',
      lastModified: 1
    })]);
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
  const chosenIdentity = await chooseRecord(recordSelect, 'BGC0000709');
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
  expect(reverseOnly.canonicalRecord.windows).toEqual([]);
  expect(reverseOnly.canonicalRecord.presentation.reverseComplement).toBe(true);

  await page.getByLabel('Circular region start').fill('1000');
  await page.getByLabel('Circular region end').fill('12000');
  await generateAndWaitForResult(page);
  const cropReverse = await inspectCircularResult(page);
  expect(cropReverse.text).toContain('11,001 bp');
  // Request schema 10: the record's selector names the record, and its one window is the crop.
  expect(cropReverse.canonicalRecord.selector).toEqual({ kind: 'recordId', value: 'BGC0000709' });
  expect(cropReverse.canonicalRecord.windows).toEqual([{
    start: 1000,
    end: 12000,
    reverseComplement: true,
    column: null
  }]);
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

  // An End past the record end is drawn to the record end; a Start past it is refused.
  await page.getByLabel('Circular region end').fill('61000');
  await page.getByLabel('Circular region start').fill('60000');
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
  await page.locator('input[accept^=".json,"]').first().setInputFiles(savedSessionPath);
  await expectSessionLoadedWithoutWarnings(page, { timeout: 180000 });
  await expect(page.getByLabel('Circular record', { exact: true })).toHaveValue(chosenIdentity);
  await expect(page.getByLabel('Circular region start')).toHaveValue('1000');
  await expect(page.getByLabel('Circular region end')).toHaveValue('12000');
  await expect(page.getByLabel('Circular reverse complement')).toBeChecked();

  await generateAndWaitForResult(page);
  const regenerated = await inspectCircularResult(page);
  expect(regenerated.content).toBe(cropReverse.content);
  expect(regenerated.canonicalRecord.windows.map((window) => window.reverseComplement)).toEqual([true]);
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
  (await import('./js/state.js')).state.activeDrawing().circularRecordDiscovery.status
))).toBe('ready');
// The Circular record select names a record by its identity (`<card uid>`,
// or `<card uid>:<k>` in a file of several records): the option of record ID
// `record`, or of the k-th record of the card for `#k`.
const chooseRecord = async (select, record) => {
  const value = () => select.locator('option').evaluateAll((options, key) => {
    const index = /^#(\d+)$/.exec(key);
    const found = options.filter((option) => option.value !== '').find((option) => (index
      ? option.value.endsWith(`:${index[1]}`) : option.textContent.includes(key)));
    return found ? found.value : null;
  }, record);
  await expect.poll(value).not.toBeNull();
  await select.selectOption(await value());
  return value();
};
const nativeUpload = (page, content, name = 'presentation.gb') => page.getByLabel('GenBank/DDBJ File', { exact: true })
  .setInputFiles({ name, mimeType: 'text/plain', buffer: Buffer.from(content) });

for (const width of [1280, 390]) {
  test(`applicable disclosure keeps manual close, focus, scroll and keyboard at ${width}px`, async ({ page }, testInfo) => {
    await page.setViewportSize({ width, height: 900 });
    await openApp(page);
    await page.evaluate(() => { window.__GBDRAW_APP__.form.multi_record_canvas = false; });
    const upload = page.getByRole('button', { name: 'Choose GenBank/DDBJ File', exact: true });
    await upload.focus();
    await nativeUpload(page, singleSource);
    await ready(page);
    // The first file makes File card 1, which takes the place of the empty input's uploader: the focus moves to
    // the card's file button, opened to show it (OV-2905), as a Linear card keeps it on its own.
    const cardFile = page.locator('[data-circular-source-card]').first().getByRole('button', { name: 'Choose GenBank/DDBJ File', exact: true });
    await expect(cardFile).toBeFocused();
    await expect(cardFile).toBeInViewport();
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
    await chooseRecord(selector, 'BGC0000709');
    await expect(panel(page)).toHaveAttribute('open', '');
    await expect(selector).toBeFocused();
    await expect.poll(async () => Math.abs((await selector.boundingBox()).y - before.y))
      .toBeLessThanOrEqual(2);
    // The Record select writes the record's identity.
    expect(await page.evaluate(() => window.__GBDRAW_APP__.form.circular_record_selector)).toBe(await page.evaluate(() => (
      window.__GBDRAW_APP__.circularRecordList.find((record) => record.record_id === 'BGC0000709').identity)));
    expect(await page.evaluate(() => window.__GBDRAW_APP__.form.multi_record_canvas)).toBe(false);
    await summary.focus();
    await page.keyboard.press('Enter');
    await page.evaluate(() => { window.__GBDRAW_APP__.adv.def_font_size += 1; });
    await expect(panel(page)).not.toHaveAttribute('open', '');
    await chooseRecord(selector, 'BGC0000708');
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
  await page.getByRole('button', { name: 'Show Multi-Record Canvas setting' }).click();
  const grid = page.locator('[data-circular-canvas-setting]');
  await expect(grid).toBeFocused();
  await expect(grid).toBeChecked();
  await page.keyboard.press('Space');
  await chooseRecord(page.getByLabel('Circular record', { exact: true }), 'BGC0000709');
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
  await chooseRecord(page.getByLabel('Circular record', { exact: true }), 'NC_012920.1');
  await expect(page.getByLabel('Circular region start')).toBeEnabled();
});

test('History restore re-inspects sources and reconciles missing or ambiguous selectors without changing crop intent', async ({ page }) => {
  await openApp(page);
  await page.evaluate(() => { window.__GBDRAW_APP__.form.multi_record_canvas = false; });
  await nativeUpload(page, singleSource);
  await ready(page);
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(1);
  await page.getByRole('button', { name: 'Undo', exact: true }).click();
  await expect.poll(() => page.evaluate(async () => (await import('./js/state.js')).state.activeDrawing().circularRecordDiscovery.status)).toBe('idle');
  await page.getByRole('button', { name: 'Redo', exact: true }).click();
  await ready(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.circularRecordList.length)).toBe(1);
  await page.getByLabel('Circular region start').fill('25');
  await page.getByLabel('Circular region end').fill('100');
  await page.getByLabel('Circular region end').blur();
  await page.evaluate(() => { window.__GBDRAW_APP__.form.circular_record_selector = 'missing'; });
  await expect(page.getByLabel('Circular region start')).toBeDisabled();
  await expect(page.getByText("No record in the current input matches the record selector 'missing'.", { exact: true })).toBeVisible();
  await nativeUpload(page, `${singleSource}\n${singleSource}`);
  await ready(page);
  await page.evaluate(() => { window.__GBDRAW_APP__.form.circular_record_selector = 'NC_012920.1'; });
  await expect(page.getByLabel('Circular region start')).toBeDisabled();
  await expect(page.getByText("Record selector 'NC_012920.1' is ambiguous in the current input.", { exact: true })).toBeVisible();
  await chooseRecord(page.getByLabel('Circular record', { exact: true }), '#2');
  await expect(page.getByLabel('Circular region start')).toBeEnabled();
  await expect(page.getByLabel('Circular region start')).toHaveValue('25');
  await expect(page.getByLabel('Circular region end')).toHaveValue('100');
});

test('Loaded Session shows the saved record selector as uninspected, not as missing, until inspection', async ({ page }) => {
  test.setTimeout(420000);
  await openApp(page);
  await page.evaluate(() => {
    window.__GBDRAW_APP__.form.multi_record_canvas = false;
    window.__GBDRAW_APP__.sessionTitle = 'ov41-record-selector';
  });
  await nativeUpload(page, singleSource);
  await ready(page);
  const chosen = await chooseRecord(page.getByLabel('Circular record', { exact: true }), 'NC_012920.1');
  await generateAndWaitForResult(page);
  const downloadPromise = page.waitForEvent('download', { timeout: 180000 });
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.saveSessionWithTitle());
  const savedSessionPath = await (await downloadPromise).path();

  await page.reload({ waitUntil: 'domcontentloaded' });
  await waitForAppShell(page, { waitForPalette: false });
  await page.locator('input[accept^=".json,"]').first().setInputFiles(savedSessionPath);
  await expect(page.locator('[data-app-notice]')).toContainText('Session loaded', { timeout: 180000 });

  const select = page.getByLabel('Circular record', { exact: true });
  await expect(select).toHaveValue(chosen);
  const shot = async (name) => {
    if (!process.env.OV41_SHOT_DIR) return;
    const section = page.locator('section[aria-label="Source records"]');
    await section.scrollIntoViewIfNeeded();
    await section.screenshot({ path: join(process.env.OV41_SHOT_DIR, name) });
  };
  // A Session load does not inspect the source; the saved selector is unverified, not missing.
  await expect(page.locator('[data-circular-discovery-status]')).toHaveText('Records not inspected');
  await shot('record-control-loaded.png');
  // The record of a file of one record has no record ID in the request; the record ID the loaded Result drew it
  // with names it, as before the File cards (OV-2959, D-68).
  await expect(select.locator('option:checked')).toHaveText('NC_012920.1 (not inspected)');

  // Inspection resolves the selector to the real record entry.
  await page.getByRole('button', { name: 'Inspect source records', exact: true }).click();
  await ready(page);
  await expect(select).toHaveValue(chosen);
  await expect(select.locator('option:checked')).toHaveText('NC_012920.1 (16,569 bp)');
  await shot('record-control-inspected.png');
});

// CI-03: a replaced or removed multi-record source retires the record selector,
// crop, orientation, and title lines chosen for one of its records, in the
// same History step; one record replacing one record keeps them (D-32).
test('replacing or removing a multi-record source resets its one-record presentation in one History step', async ({ page }) => {
  test.setTimeout(300000);
  await openApp(page);
  await page.evaluate(() => { window.__GBDRAW_APP__.form.multi_record_canvas = false; });
  const presentation = () => page.evaluate(() => {
    const { form, adv } = window.__GBDRAW_APP__;
    return [form.circular_record_selector, form.circular_region_start, form.circular_region_end, form.circular_reverse,
      form.circular_record_label, form.circular_record_subtitle, adv.circular_grouping_intent];
  });
  const clean = ['', null, null, false, '', '', 'auto'];
  // The Record select writes the chosen record's identity.
  let chosen = '';
  const choose = async () => {
    await nativeUpload(page, `${firstGenbank}\n${secondGenbank}`, 'two.gb');
    await ready(page);
    chosen = await chooseRecord(page.getByLabel('Circular record', { exact: true }), 'BGC0000709');
    await page.evaluate(() => Object.assign(window.__GBDRAW_APP__.form, {
      circular_region_start: 500, circular_region_end: 3000, circular_reverse: true,
      circular_record_label: 'Label X', circular_record_subtitle: 'Sub Y'
    }));
    expect(await presentation()).toEqual([chosen, 500, 3000, true, 'Label X', 'Sub Y', 'single']);
  };
  await choose();
  const undoBefore = await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
  await nativeUpload(page, singleSource);
  await ready(page);
  expect(await presentation()).toEqual(clean);
  expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(undoBefore + 1);
  await page.getByRole('button', { name: 'Undo', exact: true }).click();
  await expect.poll(presentation).toEqual([chosen, 500, 3000, true, 'Label X', 'Sub Y', 'single']);
  await page.getByRole('button', { name: 'Redo', exact: true }).click();
  await expect.poll(presentation).toEqual(clean);
  await ready(page);
  await generateAndWaitForResult(page);

  await choose();
  await page.getByRole('button', { name: 'Remove File 1', exact: true }).click();
  await page.locator('[data-circular-source-removal-dialog] [role=dialog]').getByRole('button', { name: 'Remove file', exact: true }).click();
  await expect.poll(presentation).toEqual(clean);
});

// CI-03 / D-32: one record replacing one record with a different record ID
// keeps the crop and Record label but retires the selector naming the old
// record, so Generate draws the new record instead of failing RECORD_SELECTION.
test('one record replacing one record with another ID keeps the crop and still generates', async ({ page }) => {
  test.setTimeout(300000);
  await openApp(page);
  await page.evaluate(() => { window.__GBDRAW_APP__.form.multi_record_canvas = false; });
  await nativeUpload(page, firstGenbank, 'one.gb');
  await ready(page);
  await chooseRecord(page.getByLabel('Circular record', { exact: true }), 'BGC0000708');
  await page.evaluate(() => Object.assign(window.__GBDRAW_APP__.form, {
    circular_region_start: 500, circular_region_end: 3000, circular_record_label: 'Label X'
  }));
  await nativeUpload(page, secondGenbank, 'other.gb');
  await ready(page);
  const outcome = await evaluateWithRetainedPromise(page, async () => {
    const app = window.__GBDRAW_APP__;
    const result = await app.runAnalysis();
    return { status: result?.status, code: app.errorLog?.code || null,
      selector: app.form.circular_record_selector, start: app.form.circular_region_start,
      end: app.form.circular_region_end, label: app.form.circular_record_label };
  });
  expect(outcome).toEqual({ status: 'ok', code: null, selector: '', start: 500, end: 3000, label: 'Label X' });
  expect((await inspectCircularResult(page)).recordIds).toContain('BGC0000709');
});

// OV-3021: a Record Order row names its record after its File card's file name. Every word of
// the row stays whole and in view, the record ID and the OFF badge too; the label wraps between
// words, and only a word longer than the row (a long file name) may end in an ellipsis.
const recordOrderRows = (page) => page.locator('[data-capture="multi-record-order"]').evaluate((root) => (
  [...root.querySelectorAll('select[aria-label^="Row for "]')].map((select) => {
    const label = select.parentElement.firstElementChild;
    const box = label.getBoundingClientRect();
    const words = [];
    const walker = document.createTreeWalker(label, NodeFilter.SHOW_TEXT);
    for (let node = walker.nextNode(); node; node = walker.nextNode()) {
      for (const match of node.data.matchAll(/\S+/g)) {
        const range = document.createRange();
        range.setStart(node, match.index);
        range.setEnd(node, match.index + match[0].length);
        const rects = [...range.getClientRects()].filter((rect) => rect.width > 0);
        words.push({
          word: match[0],
          lines: new Set(rects.map((rect) => Math.round(rect.top))).size,
          inside: rects.every((rect) => rect.left >= box.left - 0.5 && rect.right <= box.right + 0.5)
        });
      }
    }
    return {
      text: label.textContent.replace(/\s+/g, ' ').trim(),
      title: label.getAttribute('title'),
      name: select.getAttribute('aria-label'),
      clipped: label.scrollWidth > label.clientWidth,
      words
    };
  })
));

for (const width of [1280, 390]) {
  test(`Record Order rows of several File cards show each record ID and OFF badge whole at ${width}px`, async ({ page }) => {
    test.setTimeout(300000);
    await page.setViewportSize({ width, height: 900 });
    await openApp(page);
    const longName = 'Streptomyces_fradiae_tylosin_cluster_v2.gbk';
    await page.evaluate(({ first, second, longName: name }) => {
      const app = window.__GBDRAW_APP__;
      const file = (text, fileName) => new File([text], fileName, { type: 'text/plain', lastModified: 1 });
      void app.addCircularSourceFiles('gb', [
        file(`${first}\n${first}`, 'dup.gbk'), file(second, 'b.gbk'), file(`${second}\n${first}`, name)
      ]);
      app.form.multi_record_canvas = true;
    }, { first: firstGenbank, second: secondGenbank, longName });
    await ready(page);
    await page.locator('[data-record-list-open="circular"]').click();
    const list = page.locator('[data-record-list-dialog]');
    await list.locator('[data-record-list-row]').nth(1).uncheck();
    await list.locator('[data-record-list-row]').nth(4).uncheck();
    await list.getByRole('button', { name: 'Close', exact: true }).click();
    await page.getByRole('button', { name: 'Show Multi-Record Canvas setting', exact: true }).click();
    await expect(page.locator('[data-capture="multi-record-order"] select')).toHaveCount(5);

    const rows = await recordOrderRows(page);
    for (const row of rows) {
      // No word wraps inside itself.
      expect(row.words.filter((word) => word.lines !== 1)).toEqual([]);
      // Every word but a file name longer than the row is in view: the record ID and the badge.
      expect(row.words.filter((word) => !word.inside && word.word !== `${longName}:`)).toEqual([]);
      if (!row.text.startsWith(longName)) expect(row.clipped).toBe(false);
      // The whole label stays in the row's accessible name and in its title.
      const label = row.text.replace(/ OFF$/, '');
      expect(row.name).toBe(`Row for ${label}`);
      expect(row.title).toBe(label);
    }
    // The badge is a word of its own, also for a screen reader.
    expect(rows.map((row) => row.text)).toEqual([
      'dup.gbk: BGC0000708 [#1]', 'dup.gbk: BGC0000708 [#2] OFF', 'b.gbk: BGC0000709',
      `${longName}: BGC0000709`, `${longName}: BGC0000708 OFF`
    ]);
    if (width === 390) {
      expect(await page.evaluate(() => document.documentElement.scrollWidth)).toBeLessThanOrEqual(width);
    }
  });
}
