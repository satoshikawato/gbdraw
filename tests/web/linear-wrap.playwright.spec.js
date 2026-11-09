// Wrap Rows (Linear Layout, PD-OI-113): turning it on switches the settings rows
// decide in one History step, Auto shows the row length and rows Python chose,
// comparisons and a record moved sideways keep it from being turned on (also in
// a loaded Session, whose Generate then refuses and keeps the figure), a feature
// cut over rows stays one feature, Save and Load keep the length, and a PNG or a
// wrapped PDF over its size limit is refused before a file is made, also for a
// loaded Result, a cropped record's included.
const { test, expect } = require('@playwright/test');
const { execFileSync } = require('node:child_process');
const { copyFileSync, mkdirSync, readFileSync, writeFileSync } = require('node:fs');
const path = require('node:path');
const { evaluateWithRetainedPromise, generateAndWaitForResult } = require('./helpers/app-lifecycle.cjs');
const { showDrawing } = require('./helpers/drawing-bar.cjs');
const { loadSessionFile, openFresh, settle } = require('./helpers/audit-browser.cjs');
const { download } = require('./helpers/mode-transition.cjs');

const root = path.resolve(process.env.GBDRAW_REPO || process.cwd());
const data = path.join(root, 'gbdraw', 'web', 'tutorial-data');
const LAMBDA = path.join(data, 'lambda', 'NC_001416.gb');
const MITOGENOMES = [
  path.join(data, 'human-mitochondrion', 'HmmtDNA.gbk'),
  path.join(data, 'metazoan-mitochondria-four', 'NC_002333.2.gb')
];
// The texts of the refusals the field's hint shares (REASONS in gbdraw/web/js/utils/error-normalization.js).
const COMPARISONS_TEXT = 'A wrapped figure draws no comparisons.';
const MOVED_TEXT = 'A wrapped record cannot keep a sideways move (from Align at feature, Similarity groups, or a loaded Session). Click Reset record positions under Wrap Rows, or turn Wrap Rows off.';

const openLinear = async (page, files) => {
  await openFresh(page);
  await showDrawing(page, { kind: 'linear' });
  await page.evaluate(async (inputs) => {
    const app = window.__GBDRAW_APP__;
    while (app.linearSeqs.length < inputs.length) app.addLinearSeq();
    inputs.forEach(({ name, text }, index) => app.setLinearSeqPrimaryFile(index, 'gb', new File([text], name, { type: 'text/plain' })));
    await window.Vue.nextTick();
  }, files.map((file) => ({ name: path.basename(file), text: readFileSync(file, 'utf8') })));
  await settle(page);
};

// The five fields one Wrap Rows edit writes, and what the field shows.
const wrapState = (page) => page.evaluate(() => {
  const { form } = window.__GBDRAW_APP__;
  const text = (selector) => document.querySelector(selector)?.textContent?.trim() ?? null;
  return {
    wrap: form.linear_wrap_length,
    layout: form.linear_track_layout,
    scale: form.scale_style,
    normalize: form.normalize_length,
    rulerOnAxis: form.linear_ruler_on_axis,
    middleDisabled: document.querySelector('#linear-track-layout option[value="middle"]')?.disabled ?? null,
    amount: document.querySelector('#linear-wrap-length')?.value ?? null,
    unit: document.querySelector('#linear-wrap-unit')?.value ?? null,
    blocked: text('#linear-wrap-blocked')
  };
});

const wrapBox = (page) => page.locator('#linear-wrap-rows');
// The record translations the page holds (`state.linearRecordTranslations`), the ones the next request carries.
const translations = (page) => page.evaluate(async () => {
  const { state } = await import('/gbdraw/web/js/state.js');
  return state.linearRecordTranslations.value.map(({ x, y }) => ({ x, y }));
});
const undo = async (page) => { await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo()); await settle(page); };
const setWrap = async (page, amount, unit) => {
  if (!(await wrapBox(page).isChecked())) await wrapBox(page).check();
  await page.locator('#linear-wrap-unit').selectOption(unit);
  await page.locator('#linear-wrap-length').fill(String(amount));
  await page.locator('#linear-wrap-length').dispatchEvent('change');
  await settle(page);
};

test('turning Wrap Rows on switches the settings rows decide, one Undo restores them, and Auto shows what Python drew', async ({ page }) => {
  test.setTimeout(300_000);
  await openLinear(page, [LAMBDA]);
  await page.evaluate(() => {
    const { form } = window.__GBDRAW_APP__;
    form.linear_track_layout = 'middle';
    form.scale_style = 'bar';
    form.normalize_length = true;
    form.linear_ruler_on_axis = false;
  });
  await settle(page);
  const before = await wrapState(page);
  expect(before).toMatchObject({ wrap: null, layout: 'middle', scale: 'bar', normalize: true, rulerOnAxis: false, middleDisabled: false });

  await wrapBox(page).scrollIntoViewIfNeeded();
  await wrapBox(page).check();
  await settle(page);
  expect(await wrapState(page)).toMatchObject({
    wrap: 'auto', layout: 'above', scale: 'ruler', normalize: false, rulerOnAxis: true, middleDisabled: true, amount: '', blocked: null
  });
  await expect(page.locator('#linear-wrap-fixed-track-layout')).toHaveText('Features on axis is unavailable while Wrap Rows is on.');
  await expect(page.locator('#linear-wrap-fixed-scale-style')).toHaveText('Bar (Simple) is unavailable while Wrap Rows is on.');
  await expect(page.locator('#linear-wrap-fixed-normalize')).toHaveText('Fixed while Wrap Rows is on.');
  await expect(page.locator('#linear-wrap-fixed-ruler-on-axis')).toHaveText('Fixed while Wrap Rows is on.');

  // One Undo puts back every setting the edit switched.
  await undo(page);
  expect(await wrapState(page)).toEqual(before);
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.redo());
  await settle(page);
  expect((await wrapState(page)).wrap).toBe('auto');

  // After Generate, Auto names the row length and rows Python chose (lambda, 48.5 kb), in full (OV-2467).
  await generateAndWaitForResult(page);
  await settle(page);
  const placeholder = page.locator('[data-linear-wrap-field] .auto-value-placeholder');
  await expect(placeholder).toHaveText('10 kbp per row, 5 rows (auto)');
  const clipped = await placeholder.evaluate((element) => element.scrollWidth > element.clientWidth);
  expect(clipped, 'the Auto text is not cut').toBe(false);
  expect(await page.evaluate(() => new Set([...document.querySelectorAll('[data-gbdraw-wrap-row]')]
    .map((row) => row.getAttribute('data-gbdraw-wrap-row'))).size)).toBe(5);

  // Turning it off writes only the length.
  await wrapBox(page).uncheck();
  await settle(page);
  expect(await wrapState(page)).toMatchObject({ wrap: null, layout: 'above', scale: 'ruler', normalize: false, rulerOnAxis: true, middleDisabled: false });
});

test('comparisons and a record moved by Align at feature keep Wrap Rows from being turned on and name what to remove', async ({ page }) => {
  test.setTimeout(480_000);
  await openLinear(page, MITOGENOMES);
  // NC_002333.2 drawn linear: Align at feature moves it sideways (a circular record would rotate).
  await page.evaluate(async () => {
    const controls = window.__GBDRAW_APP__.recordDisplayControls;
    const rows = controls.allRows.value ?? controls.allRows;
    await controls.setTopology(rows.find((row) => row.recordId === 'NC_002333.2'), false);
  });
  await settle(page);
  await generateAndWaitForResult(page);
  await settle(page);
  await expect(wrapBox(page)).toBeEnabled();

  // Align at feature rotates NC_012920.1 and moves the linear NC_002333.2 sideways; while a record is moved, rows
  // cannot draw it.
  const selector = page.getByLabel('Align at feature', { exact: true });
  await selector.evaluate((element) => {
    for (let details = element.closest('details'); details; details = details.parentElement?.closest('details')) details.open = true;
  });
  await selector.fill('gene=COX1');
  await page.getByTestId('linear-anchor-feature-resolve').click();
  await expect(page.getByTestId('linear-anchor-feature-preview')).toBeVisible({ timeout: 120_000 });
  await page.getByTestId('linear-anchor-feature-apply').click();
  await expect(page.getByTestId('linear-anchor-feature-notice')).toHaveText('Records aligned at gene=COX1.', { timeout: 180_000 });
  await settle(page);
  await expect(wrapBox(page)).toBeDisabled();
  await expect(page.locator('#linear-wrap-blocked')).toHaveText(MOVED_TEXT);
  // Reset record positions moves the records off Align at feature's places, so its notice goes; one Undo puts the
  // positions back and, with them, the block and the notice (Q5 A: Align at feature counts through moved records).
  const reset = page.getByTestId('linear-wrap-reset-positions');
  await reset.click();
  await expect(page.locator('#linear-wrap-blocked')).toHaveCount(0, { timeout: 180_000 });
  await settle(page);
  await expect(wrapBox(page)).toBeEnabled();
  await expect(page.getByTestId('linear-anchor-feature-notice')).toHaveCount(0);
  await undo(page);
  await expect(page.getByTestId('linear-anchor-feature-notice')).toHaveText('Records aligned at gene=COX1.');
  await expect(wrapBox(page)).toBeDisabled();
  await expect(page.locator('#linear-wrap-blocked')).toHaveText(MOVED_TEXT);
  // The Undo of Apply ends Align at feature.
  await undo(page);
  await expect(wrapBox(page)).toBeEnabled();
  await expect(page.locator('#linear-wrap-blocked')).toHaveCount(0);
  await expect(page.getByTestId('linear-anchor-feature-notice')).toHaveCount(0);

  // A comparison of the two records.
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_APP__.setLinearComparisonGlobalAction('losat'));
  await settle(page);
  await expect(wrapBox(page)).toBeDisabled();
  await expect(page.locator('#linear-wrap-blocked')).toContainText(COMPARISONS_TEXT);
});

// Q5 A (1), PD-OI-113 scenario revision 3: Align at feature counts only through the records its Apply moved. Both
// mitogenomes are complete circular records, so Apply rotates them and moves neither: Wrap Rows stays enabled, and
// Generate wraps the rotated records.
test('an Align at feature Apply that only rotated leaves Wrap Rows enabled, and Generate wraps', async ({ page }) => {
  test.setTimeout(480_000);
  await openLinear(page, MITOGENOMES);
  await generateAndWaitForResult(page);
  await settle(page);
  const selector = page.getByLabel('Align at feature', { exact: true });
  await selector.evaluate((element) => {
    for (let details = element.closest('details'); details; details = details.parentElement?.closest('details')) details.open = true;
  });
  await selector.fill('gene=COX1');
  await page.getByTestId('linear-anchor-feature-resolve').click();
  await expect(page.getByTestId('linear-anchor-feature-preview')).toBeVisible({ timeout: 120_000 });
  await page.getByTestId('linear-anchor-feature-apply').click();
  await expect(page.getByTestId('linear-anchor-feature-notice')).toHaveText('Records aligned at gene=COX1.', { timeout: 180_000 });
  await settle(page);
  expect((await translations(page)).map(({ x }) => x)).toEqual([0, 0]);
  await expect(wrapBox(page)).toBeEnabled();
  await expect(page.locator('#linear-wrap-blocked')).toHaveCount(0);

  await setWrap(page, 5, 'kbp');
  await generateAndWaitForResult(page);
  await settle(page);
  const committed = await page.evaluate(async () => {
    const { getCommittedCanonicalRenderRequest } = await import('/gbdraw/web/js/services/config.js');
    return getCommittedCanonicalRenderRequest().records.map(({ display }) => display?.startCoordinate ?? null);
  });
  expect(committed).toEqual([5904, 6425]);
  // Each 16.6 kb record is cut into four 5 kbp rows.
  expect(await page.evaluate(() => new Set([...document.querySelectorAll('[data-gbdraw-wrap-row]')]
    .map((row) => row.getAttribute('data-gbdraw-wrap-row'))).size)).toBeGreaterThan(1);
});

// Q5 A (2): a Similarity alignment counts as moved records, also beside comparisons, so its hint names Reset record
// positions and the button shows. One press puts every record at x 0 and ends the alignment in one History step,
// and one Undo brings back the plan and the positions.
test('a Similarity alignment blocks Wrap Rows beside comparisons, and Reset record positions ends it in one step', async ({ page }) => {
  test.setTimeout(480_000);
  await openFresh(page);
  await loadSessionFile(page, path.join(root, 'gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json'));
  await settle(page);
  const placement = () => page.evaluate(async () => {
    const { state } = await import('/gbdraw/web/js/state.js');
    const svg = window.__GBDRAW_APP__.svgContainer.querySelector('svg');
    return {
      plan: JSON.stringify(state.similarityAlignmentPlan.value),
      translations: state.linearRecordTranslations.value.map(({ recordKey, x, y }) => ({ recordKey, x, y })),
      drawn: [...svg.querySelectorAll('[data-record-translation-x]')].map((group) => Number(group.getAttribute('data-record-translation-x')))
    };
  });
  const before = await placement();
  expect(before.plan).not.toBe('null');
  expect(before.drawn.some((x) => x !== 0), 'the alignment draws records moved sideways').toBe(true);
  await expect(wrapBox(page)).toBeDisabled();
  await expect(page.locator('#linear-wrap-blocked')).toContainText(COMPARISONS_TEXT);
  await expect(page.locator('#linear-wrap-blocked')).toContainText(MOVED_TEXT);
  const reset = page.getByTestId('linear-wrap-reset-positions');
  await expect(reset).toBeVisible();

  const steps = await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
  await reset.click();
  await expect(page.locator('[data-similarity-alignment-notice]')).toHaveText('Alignment cleared: record positions reset.', { timeout: 180_000 });
  await settle(page);
  const after = await placement();
  expect(after.plan).toBe('null');
  expect(after.translations.map(({ x }) => x)).toEqual(before.drawn.map(() => 0));
  expect(after.drawn).toEqual(before.drawn.map(() => 0));
  expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(steps + 1);
  // The comparisons still block; the moved-record reason and its button are gone.
  await expect(page.locator('#linear-wrap-blocked')).toContainText(COMPARISONS_TEXT);
  await expect(page.locator('#linear-wrap-blocked')).not.toContainText(MOVED_TEXT);
  await expect(reset).toHaveCount(0);

  await undo(page);
  expect(await placement()).toEqual(before);
  await expect(page.locator('#linear-wrap-blocked')).toContainText(MOVED_TEXT);
  await expect(reset).toBeVisible();
});

// PD-OI-113 scenario revision 2: a record's place is read alike whatever wrote it. A Session holding a wrap length
// and a record moved sideways (as a Python API Session can; made here from a CLI Session with one translation
// written in) keeps its figure on Load. The field cannot be changed but can be turned off, and a Generate that would
// cut the moved record is refused with the hint's words and keeps the figure. Reset record positions beside the
// field puts the record at x 0 in one step, and its Undo blocks the field again.
test('a loaded Session with a record moved sideways blocks Wrap Rows until Reset record positions, and its Generate refuses and keeps the figure', async ({ page }, testInfo) => {
  test.setTimeout(300_000);
  const directory = testInfo.outputPath('moved');
  mkdirSync(directory, { recursive: true });
  copyFileSync(LAMBDA, path.join(directory, 'NC_001416.gb'));
  const session = path.join(directory, 'moved.gbdraw-session.json');
  execFileSync(process.env.PYTHON || 'python', [
    '-m', 'gbdraw.cli', 'linear', '--gbk', 'NC_001416.gb', '--wrap_length', '5000', '-o', 'lambda',
    '--save_session', '--session_output', session
  ], { cwd: directory, stdio: 'pipe', env: { ...process.env, PYTHONPATH: root } });
  const document = JSON.parse(readFileSync(session, 'utf8'));
  const request = document.drawings[0].renderRequest;
  // The record 40 px to the right, as `LinearRecordTranslation(record_key, x=40)` stores it.
  request.layout = { recordGapPx: 24, multiRecordPositions: null, similarityAlignment: null,
    recordTranslations: [{ recordKey: request.records[0].recordKey, x: 40, y: 0 }] };
  writeFileSync(session, JSON.stringify(document));

  await openFresh(page);
  await loadSessionFile(page, session);
  expect(await wrapState(page)).toMatchObject({ wrap: 5000, amount: '5', unit: 'kbp', blocked: MOVED_TEXT });
  await expect(page.locator('#linear-wrap-length')).toBeDisabled();
  const shown = await page.evaluate(() => window.__GBDRAW_APP__.currentResults()[0].content);
  expect(shown).toContain('data-gbdraw-wrap-row');

  const outcome = await generateAndWaitForResult(page, { expectedStatus: 'error' });
  expect(outcome.errorSummary).toContain(MOVED_TEXT);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.currentResults()[0].content), 'the figure is kept').toBe(shown);

  // The button the hint names (Q4): one press puts every record at x 0 as one History step and draws it, so Wrap
  // Rows can be changed again; one Undo puts the move back, and the field is blocked again.
  const reset = page.getByTestId('linear-wrap-reset-positions');
  await expect(reset).toBeVisible();
  const steps = await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
  await reset.click();
  await expect(page.locator('#linear-wrap-blocked')).toHaveCount(0, { timeout: 180_000 });
  await settle(page);
  expect(await translations(page)).toEqual([{ x: 0, y: 0 }]);
  expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount())).toBe(steps + 1);
  await expect(reset).toHaveCount(0);
  await expect(page.locator('#linear-wrap-length')).toBeEnabled();
  expect(await page.evaluate(() => window.__GBDRAW_APP__.currentResults()[0].content)).toContain('data-gbdraw-wrap-row');
  await undo(page);
  expect(await translations(page)).toEqual([{ x: 40, y: 0 }]);
  expect(await wrapState(page)).toMatchObject({ wrap: 5000, blocked: MOVED_TEXT });
  await expect(page.locator('#linear-wrap-length')).toBeDisabled();
  await expect(reset).toBeVisible();

  // The other way out: Wrap Rows turns off while blocked, and then cannot be turned on again.
  await expect(wrapBox(page)).toBeEnabled();
  await wrapBox(page).uncheck();
  await settle(page);
  expect(await wrapState(page)).toMatchObject({ wrap: null, blocked: MOVED_TEXT });
  await expect(wrapBox(page)).toBeDisabled();
});

test('a feature cut over rows opens one popup from any fragment and highlights every fragment', async ({ page }) => {
  test.setTimeout(300_000);
  await openLinear(page, [LAMBDA]);
  await wrapBox(page).check();
  await generateAndWaitForResult(page);
  await settle(page);
  // A feature drawn on two rows, and the centre of its fragment on the later row where nothing covers it.
  const cut = await page.evaluate(() => {
    const svg = window.__GBDRAW_APP__.svgContainer.querySelector('svg');
    const byFeature = new Map();
    svg.querySelectorAll('[data-gbdraw-feature-part="block"][data-gbdraw-feature-id]').forEach((part) => {
      const row = part.closest('[data-gbdraw-wrap-row]')?.getAttribute('data-gbdraw-wrap-row');
      const id = part.getAttribute('data-gbdraw-feature-id');
      if (!byFeature.has(id)) byFeature.set(id, []);
      byFeature.get(id).push({ part, row: Number(row) });
    });
    for (const [id, parts] of byFeature) {
      if (new Set(parts.map(({ row }) => row)).size < 2) continue;
      const later = parts.reduce((last, entry) => (entry.row > last.row ? entry : last));
      const box = later.part.getBoundingClientRect();
      const x = box.left + box.width / 2;
      const y = box.top + box.height / 2;
      if (document.elementFromPoint(x, y) === later.part) return { id, x, y, rows: new Set(parts.map(({ row }) => row)).size };
    }
    return null;
  });
  expect(cut, 'a cut feature whose later fragment is visible').not.toBeNull();
  expect(cut.rows).toBeGreaterThan(1);
  const fragments = (attribute) => page.evaluate(({ id, name }) => {
    const parts = [...document.querySelectorAll(`[data-gbdraw-feature-part="block"][data-gbdraw-feature-id="${id}"]`)];
    return { count: parts.length, marked: parts.filter((part) => (name.startsWith('.') ? part.classList.contains(name.slice(1)) : part.hasAttribute(name))).length };
  }, { id: cut.id, name: attribute });

  await page.mouse.move(cut.x, cut.y);
  await expect.poll(async () => {
    const { count, marked } = await fragments('data-gbdraw-hover-opacity');
    return marked === count;
  }).toBe(true);
  const hovered = await fragments('data-gbdraw-hover-opacity');
  expect(hovered.count).toBeGreaterThan(1);
  expect(hovered.marked).toBe(hovered.count);

  await page.mouse.click(cut.x, cut.y);
  await expect(page.locator('[role="dialog"].feature-popup')).toHaveCount(1);
  const svgId = await page.evaluate(() => window.__GBDRAW_APP__.clickedFeature?.svg_id ?? null);
  expect(svgId).not.toBeNull();

  // A Shift-drag over the later fragment selects the feature: the selection marks every fragment.
  await page.keyboard.press('Escape');
  await page.keyboard.down('Shift');
  await page.mouse.move(cut.x - 3, cut.y - 3);
  await page.mouse.down();
  await page.mouse.move(cut.x + 3, cut.y + 3, { steps: 4 });
  await page.mouse.up();
  await page.keyboard.up('Shift');
  await expect.poll(() => page.evaluate((id) => window.__GBDRAW_APP__.selectedFeatureIds.has(id), svgId)).toBe(true);
  const selected = await fragments('.gbdraw-feature-selected');
  expect(selected.marked).toBe(selected.count);
});

test('Save then Load keeps the length, a loaded wrapped PDF over its limit and a PNG over its limit make no file', async ({ page }, testInfo) => {
  test.setTimeout(480_000);
  const downloads = [];
  page.on('download', (file) => downloads.push(file.suggestedFilename()));
  await openLinear(page, [LAMBDA]);
  await wrapBox(page).scrollIntoViewIfNeeded();
  // 300 bp rows: 162 rows of lambda, a PDF page of about 16,400 pt, which the CLI refuses.
  await setWrap(page, 300, 'bp');
  expect(await wrapState(page)).toMatchObject({ wrap: 300, amount: '300', unit: 'bp' });
  await generateAndWaitForResult(page);
  const saved = testInfo.outputPath('wrapped.gbdraw-session.json.gz');
  await download(page, 'Save Session', saved, { title: '' });
  await wrapBox(page).uncheck();
  await settle(page);
  expect((await wrapState(page)).wrap).toBe(null);
  await loadSessionFile(page, saved);
  expect(await wrapState(page)).toMatchObject({ wrap: 300, amount: '300', unit: 'bp' });

  // Review 3 F1: the loaded Result has no input summary and its record list waits; Export PDF without Generate
  // reads the record first and refuses the page as the CLI does.
  downloads.length = 0;
  await page.getByRole('button', { name: 'PDF', exact: true }).click();
  await expect(page.getByText('PDF is limited to 14,400 pt per side', { exact: false }).first()).toBeVisible({ timeout: 60_000 });
  await page.waitForTimeout(1_000);
  expect(downloads).toEqual([]);

  // A Session the CLI wrote with --wrap_length shows its length in the field.
  const directory = testInfo.outputPath('cli');
  mkdirSync(directory, { recursive: true });
  copyFileSync(LAMBDA, path.join(directory, 'NC_001416.gb'));
  const cliSession = path.join(directory, 'cli.gbdraw-session.json');
  execFileSync(process.env.PYTHON || 'python', [
    '-m', 'gbdraw.cli', 'linear', '--gbk', 'NC_001416.gb', '--wrap_length', '5000', '-o', 'lambda',
    '--save_session', '--session_output', cliSession
  ], { cwd: directory, stdio: 'pipe', env: { ...process.env, PYTHONPATH: root } });
  await loadSessionFile(page, cliSession);
  expect(await wrapState(page)).toMatchObject({ wrap: 5000, amount: '5', unit: 'kbp' });

  // 1 kbp rows at 600 DPI pass PNG's 32,767 px per side: the message, and no download.
  await setWrap(page, 1, 'kbp');
  await generateAndWaitForResult(page);
  await settle(page);
  await page.evaluate(() => { window.__GBDRAW_APP__.downloadDpi = '600'; });
  downloads.length = 0;
  await page.getByRole('button', { name: 'PNG', exact: true }).click();
  const message = page.getByText('PNG is limited to 32,767 px per side', { exact: false }).first();
  await expect(message).toBeVisible({ timeout: 60_000 });
  await expect(page.getByText('Use SVG, which has no size limit.', { exact: false }).first()).toBeVisible();
  await page.waitForTimeout(1_000);
  expect(downloads).toEqual([]);
});

// D-59 residual (D-70, D-73): a loaded Result whose record is cropped. The Web sends the region the Result drew, from
// its committed request, with the source length of the record the card reads, and Python crops them by the
// drawing's rule: 30 kbp of lambda in 150 bp rows (200 rows, about 20,300 pt) is refused as the CLI refuses it, also
// after the card's Region is edited to 100 bp with no Generate; in 300 bp rows (about 10,200 pt) the PDF is written.
test('a loaded Result with a cropped wrapped record: a PDF over the limit makes no file, one under it downloads', async ({ page }, testInfo) => {
  test.setTimeout(300_000);
  const downloads = [];
  page.on('download', (file) => downloads.push(file.suggestedFilename()));
  const directory = testInfo.outputPath('cropped');
  mkdirSync(directory, { recursive: true });
  copyFileSync(LAMBDA, path.join(directory, 'NC_001416.gb'));
  const cliSession = (rowLength) => {
    const file = path.join(directory, `crop${rowLength}.gbdraw-session.json`);
    execFileSync(process.env.PYTHON || 'python', [
      '-m', 'gbdraw.cli', 'linear', '--gbk', 'NC_001416.gb', '--region', 'NC_001416.1:5001-35000',
      '--wrap_length', String(rowLength), '-o', `crop${rowLength}`, '--save_session', '--session_output', file
    ], { cwd: directory, stdio: 'pipe', env: { ...process.env, PYTHONPATH: root } });
    return file;
  };
  const cardRegion = () => page.evaluate(() => {
    const [seq] = window.__GBDRAW_APP__.linearSeqs;
    return [Number(seq.region_start), Number(seq.region_end)];
  });
  await openFresh(page);
  await loadSessionFile(page, cliSession(150));
  expect(await wrapState(page)).toMatchObject({ wrap: 150, amount: '150', unit: 'bp' });
  expect(await cardRegion()).toEqual([5001, 35000]);
  // Review 1 M2: the edited card would cut nothing, but the Result shown drew 5001..35000.
  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    app.setLinearCardWindow(app.linearSeqs[0], 0, 'end', 5100);
  });
  await settle(page);
  expect(await cardRegion()).toEqual([5001, 5100]);
  downloads.length = 0;
  await page.getByRole('button', { name: 'PDF', exact: true }).click();
  await expect(page.getByText('PDF is limited to 14,400 pt per side', { exact: false }).first()).toBeVisible({ timeout: 60_000 });
  await page.waitForTimeout(1_000);
  expect(downloads).toEqual([]);

  await loadSessionFile(page, cliSession(300));
  expect(await wrapState(page)).toMatchObject({ wrap: 300, amount: '300', unit: 'bp' });
  expect(await cardRegion()).toEqual([5001, 35000]);
  const [pdf] = await Promise.all([
    page.waitForEvent('download', { timeout: 120_000 }),
    page.getByRole('button', { name: 'PDF', exact: true }).click()
  ]);
  expect(pdf.suggestedFilename()).toMatch(/\.pdf$/);
});
