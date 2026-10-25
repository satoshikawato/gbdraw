// Labels that do not fit (LP-D43, LP-D45, LP-D48, OD-F14): overlapping labels
// open the dialog after a Generate that drew nothing; crossing leaders show a
// notice bar over the drawn figure, and its "Choose a fix" computes the rows
// Generate did not search. Apply is one History step and Generates; Cancel
// ("Keep as drawn" over a drawn figure) writes nothing. Linear region labels over other labels (LP-D68) and
// Circular region labels cut to fit their circle (LP-D87) show the notice with
// no fix (LP-D82).
const { readFileSync } = require('node:fs');
const { join } = require('node:path');
const { test, expect } = require('@playwright/test');
const { evaluateWithRetainedPromise, generateAndWaitForResult, openApp } = require('./helpers/app-lifecycle.cjs');

const ROW_KINDS = ['outward_both', 'outward_vertical', 'outward_horizontal', 'smaller_font', 'tighter_rows', 'combination'];
const SCREENSHOTS = process.env.GBDRAW_LABEL_FIT_SCREENSHOTS || '';

const SARS = () => readFileSync(join(process.cwd(), 'tests/test_inputs/SARS-CoV-2.gbk'), 'utf8');
const LAMBDA = () => readFileSync(join(process.cwd(), 'tests/test_inputs/NC_001416.gb'), 'utf8');

const openWithOuterLabels = async (page, { name, text, fontSize = null, rendering = 'auto' }) => {
  await page.setViewportSize({ width: 1440, height: 1000 });
  await openApp(page);
  await page.getByLabel('GenBank/DDBJ File', { exact: true })
    .setInputFiles({ name, mimeType: 'text/plain', buffer: Buffer.from(text) });
  await page.waitForFunction(() => window.__GBDRAW_APP__.circularRecordList.length === 1);
  await page.evaluate(({ size, labelRendering }) => {
    const app = window.__GBDRAW_APP__;
    app.form.labels_mode = 'out';
    if (size !== null) app.adv.label_font_size = size;
    app.adv.label_rendering = labelRendering;
  }, { size: fontSize, labelRendering: rendering });
};
const undoCount = (page) => page.evaluate(() => window.__GBDRAW_HISTORY__.getUndoCount());
const labelSettings = (page) => page.evaluate(() => {
  const { adv } = window.__GBDRAW_APP__;
  return { x: adv.outer_label_x_offset, y: adv.outer_label_y_offset, font: adv.label_font_size,
    spacing: adv.circular_label_spacing, resolveOverlaps: adv.resolve_overlaps };
});
const rowStatuses = (dialog) => dialog.locator('[data-label-fit-row]').evaluateAll((rows) => (
  rows.map((row) => [row.getAttribute('data-label-fit-row'), row.getAttribute('data-label-fit-status')])
));
const shot = async (locator, name) => {
  if (SCREENSHOTS) await locator.screenshot({ path: join(SCREENSHOTS, name) });
};
// The page at desktop and phone width, for the PR's Before | After images.
const pageShots = async (page, name) => {
  if (!SCREENSHOTS) return;
  for (const width of [1440, 390]) {
    await page.setViewportSize({ width, height: 1000 });
    await page.getByText('Result Preview', { exact: true }).scrollIntoViewIfNeeded();
    await page.screenshot({ path: join(SCREENSHOTS, `${name}-${width}.png`) });
  }
  await page.setViewportSize({ width: 1440, height: 1000 });
};

// LP-D52 (E-ROW): SARS-CoV-2 with labels Out overprints two embedded labels.
// With the Web's default Separate Strands, "separate overlapping features"
// removes the overlap but its outer leaders still cross: the row is offered,
// says so, and draws with the notice bar (LP-D57, OV-724).
test('overlapping labels: the dialog shows every row; Cancel writes nothing; Apply is one step that draws; Undo restores', async ({ page }) => {
  test.setTimeout(600000);
  await openWithOuterLabels(page, { name: 'SARS-CoV-2.gbk', text: SARS() });
  await generateAndWaitForResult(page, { expectedStatus: 'error', requireCommittedResult: false });
  const dialog = page.locator('[data-label-fit-dialog]');
  await expect(dialog).toBeVisible();
  await expect(dialog.getByRole('heading')).toHaveText('Labels do not fit');
  await expect(dialog.locator('[data-label-fit-records]')).toContainText('3 embedded labels at 14 pt (1 label pair overlaps');
  // Python searched every row before the dialog opened, in its order.
  const statuses = await rowStatuses(dialog);
  expect(statuses.map(([kind]) => kind)).toEqual([...ROW_KINDS, 'radial', 'separate_features']);
  expect(statuses.filter(([, status]) => status !== 'not_found')).toEqual([['separate_features', 'found_with_notice']]);
  for (const [kind, status] of statuses) {
    await expect(dialog.locator(`[data-label-fit-row="${kind}"] input[type="radio"]`))
      .toBeEnabled({ enabled: status === 'found_with_notice' });
  }
  const separate = dialog.locator('[data-label-fit-row="separate_features"]');
  await expect(separate).toContainText('separate overlapping features');
  await expect(separate).toContainText('Resolve Overlaps On');
  await expect(separate).toContainText(/\d+ × \d+ px/);
  await expect(separate.locator('[data-label-fit-note]')).toHaveText('15 leader pairs still cross');
  await shot(dialog.getByRole('dialog'), 'label-fit-dialog-crossing-row.png');

  const before = { undo: await undoCount(page), settings: await labelSettings(page) };
  expect(before.settings.resolveOverlaps).toBe(false);
  await dialog.getByRole('button', { name: 'Cancel', exact: true }).click();
  await expect(dialog).toBeHidden();
  expect(await undoCount(page)).toBe(before.undo);
  expect(await labelSettings(page)).toEqual(before.settings);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.results.length)).toBe(0);

  await generateAndWaitForResult(page, { expectedStatus: 'error', requireCommittedResult: false });
  await expect(dialog).toBeVisible();
  // The only row with a fix is preselected.
  await expect(separate.locator('input[type="radio"]')).toBeChecked();
  const key = await page.evaluate(async () => (await import('./js/state.js')).state.resultGenerationKey.value);
  await dialog.getByRole('button', { name: 'Apply and Generate', exact: true }).click();
  await expect(dialog).toBeHidden();
  await expect.poll(() => page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    return { key: state.resultGenerationKey.value, processing: state.processing.value, error: state.errorLog.value };
  }), { timeout: 300000 }).toEqual({ key: key + 1, processing: false, error: null });
  expect(await page.evaluate(() => window.__GBDRAW_APP__.results.length)).toBe(1);
  // No overlap; the leaders the row said still cross show the notice bar.
  const notice = page.getByTestId('label-fit-notice');
  await expect(notice).toBeVisible();
  await expect(notice).toContainText('Leader lines cross (15 pairs).');
  await shot(notice, 'label-fit-notice-after-crossing-row.png');
  // The row's settings are one step (PD-OI-088); the Generate that follows is its own.
  expect(await undoCount(page)).toBe(before.undo + 2);
  expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.undoLabel())).toBe('Generate diagram');
  expect(await labelSettings(page)).toEqual({ ...before.settings, resolveOverlaps: true });
  await expect(page.getByLabel('Resolve Overlaps', { exact: true })).toBeChecked();

  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
  expect(await page.evaluate(() => window.__GBDRAW_HISTORY__.undoLabel())).toBe('Fit labels: separate overlapping features');
  await evaluateWithRetainedPromise(page, () => window.__GBDRAW_HISTORY__.undo());
  expect(await labelSettings(page)).toEqual(before.settings);
});

test('crossing leaders: the notice bar opens the dialog, which computes the rows; Keep as drawn keeps the figure', async ({ page }) => {
  test.setTimeout(600000);
  // Lambda with labels Out at 17 pt and no blacklist: no labels overlap and 10 leader pairs cross, the recipe of
  // the Python twin (tests/test_web_label_fit.py). 16 pt crossed while labels were laid out in point-mode boxes
  // (OV-165). SARS-CoV-2 with external labels only at 14 pt crossed only because labels next to the origin could
  // wrap past their neighbours (OV-1366); it fits now.
  await openWithOuterLabels(page, { name: 'NC_001416.gb', text: LAMBDA(), fontSize: 17 });
  await generateAndWaitForResult(page);
  const notice = page.getByTestId('label-fit-notice');
  await expect(notice).toBeVisible();
  await expect(notice).toContainText(/Leader lines cross \(\d+ pairs\)\./);
  // OD-F14: Generate searched no row.
  expect(await page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    return state.labelFitReports.value[0].report.rows.map((row) => row.status);
  })).toEqual(ROW_KINDS.map(() => 'not_computed'));
  await shot(notice, 'label-fit-notice.png');

  const before = { undo: await undoCount(page), settings: await labelSettings(page),
    result: await page.evaluate(() => window.__GBDRAW_APP__.currentResults()[0].content) };
  await notice.getByRole('button', { name: 'Choose a fix', exact: true }).click();
  const dialog = page.locator('[data-label-fit-dialog]');
  await expect(dialog).toBeVisible();
  await expect(dialog.getByRole('heading')).toHaveText('Labels do not fit');
  await expect(dialog).toContainText('The figure is drawn.');
  // The rows are computed one after another and show as they arrive.
  await expect.poll(async () => (await rowStatuses(dialog))
    .every(([, status]) => status === 'found' || status === 'not_found'), { timeout: 420000 }).toBe(true);
  expect((await rowStatuses(dialog)).map(([kind]) => kind)).toEqual(ROW_KINDS);
  await shot(dialog.getByRole('dialog'), 'label-fit-dialog-crossing.png');

  // K-a (PD-OI-097): over a drawn figure the dialog's Cancel reads "Keep as drawn".
  await expect(dialog.getByRole('button', { name: 'Cancel', exact: true })).toHaveCount(0);
  await dialog.getByRole('button', { name: 'Keep as drawn', exact: true }).click();
  await expect(dialog).toBeHidden();
  await expect(notice).toBeVisible();
  expect(await undoCount(page)).toBe(before.undo);
  expect(await labelSettings(page)).toEqual(before.settings);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.currentResults()[0].content)).toBe(before.result);
});

// D-12, D-18 (OV-1606), DF-D16, DF-D25, DF-D36: MERS-CoV's organism name is drawn on one line (never wrapped)
// over the tracks at the Web defaults, with the notice. The inside tracks are narrowed around it already (rule L),
// so the dialog offers the font at which the one line clears the tracks (no narrow row; its data with narrowing
// off is pinned in tests/test_web_label_fit.py), and Edit Species / Strain… (E-b) closes it and shows the
// record's name field without a History step. Apply writes Record Labels > Default font size and Generates a
// figure without the notice.
test('the definition text over a track: the notice offers a smaller font and the name field', async ({ page }) => {
  test.setTimeout(600000);
  await page.setViewportSize({ width: 1440, height: 1000 });
  await openApp(page);
  await page.getByLabel('GenBank/DDBJ File', { exact: true }).setInputFiles({
    name: 'MERS-CoV.gbk', mimeType: 'text/plain',
    buffer: readFileSync(join(process.cwd(), 'tests/test_inputs/MERS-CoV.gbk'))
  });
  await page.waitForFunction(() => window.__GBDRAW_APP__.circularRecordList.length === 1);
  await generateAndWaitForResult(page);
  const notice = page.getByTestId('label-fit-notice');
  await expect(notice).toBeVisible();
  await expect(notice).toContainText('The definition text overlaps a track (1 line).');
  await shot(notice, 'definition-notice.png');
  await pageShots(page, 'definition-notice-page');

  const openDialog = async () => {
    await notice.getByRole('button', { name: 'Choose a fix', exact: true }).click();
    const opened = page.locator('[data-label-fit-dialog]');
    await expect(opened).toBeVisible();
    return opened;
  };
  let dialog = await openDialog();
  await expect(dialog.locator('[data-label-fit-records]')).toContainText('1 definition line at 18 pt (1 line overlaps a track)');
  const row = dialog.locator('[data-label-fit-row="smaller_definition_font"]');
  await expect(row).toHaveAttribute('data-label-fit-status', 'found', { timeout: 300000 });
  expect((await rowStatuses(dialog)).map(([kind]) => kind)).toEqual(['smaller_definition_font']);
  await expect(row).toContainText('Default font size 17 pt');
  await shot(dialog.getByRole('dialog'), 'definition-dialog.png');

  // E-b: the name field of that record, found by its drawn record index; nothing is applied.
  const undo = await undoCount(page);
  const edit = dialog.locator('[data-label-fit-edit-record="0"]');
  await expect(edit).toHaveText('Edit Species / Strain…');
  await edit.click();
  await expect(dialog).toBeHidden();
  await expect(notice).toBeVisible();
  await expect.poll(() => page.evaluate(() => {
    const active = document.activeElement;
    return Boolean(active && (active.id === 'circular-species' || active.closest('[data-circular-record-presentation]')));
  })).toBe(true);
  expect(await undoCount(page)).toBe(undo);

  dialog = await openDialog();
  // The only row that clears the text is preselected.
  await expect(row.locator('input[type="radio"]')).toBeChecked();
  const key = await page.evaluate(async () => (await import('./js/state.js')).state.resultGenerationKey.value);
  await dialog.getByRole('button', { name: 'Apply and Generate', exact: true }).click();
  await expect(dialog).toBeHidden();
  await expect.poll(() => page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    return { key: state.resultGenerationKey.value, processing: state.processing.value, error: state.errorLog.value };
  }), { timeout: 300000 }).toEqual({ key: key + 1, processing: false, error: null });
  await expect(notice).toHaveCount(0);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.adv.def_font_size)).toBe(17);
  // The row's setting is one step; the Generate that follows is its own.
  expect(await undoCount(page)).toBe(undo + 2);
  await pageShots(page, 'definition-fixed-page');
});

// A 6 kb linear record with the given CDS features ([location, product]).
const genbank = (name, features) => {
  const sequence = 'acgt'.repeat(1500);
  const origin = [];
  for (let start = 0; start < sequence.length; start += 60) {
    const groups = sequence.slice(start, start + 60).match(/.{1,10}/g).join(' ');
    origin.push(`${String(start + 1).padStart(9)} ${groups}`);
  }
  const feature = ([location, product]) => [`     CDS             ${location}`,
    `                     /product="${product}"`, `                     /gene="${product}"`];
  return [`LOCUS       ${name.padEnd(16)}        6000 bp    DNA     linear   UNK 01-JAN-1980`,
    `DEFINITION  ${name}.`, `ACCESSION   ${name}`, `VERSION     ${name}`, 'FEATURES             Location/Qualifiers',
    ...features.flatMap(feature), 'ORIGIN', ...origin, '//', ''].join('\n');
};
// The labels are pinned at 10 pt, the size the OV-1155 cases were drawn at (a region label without
// its own font takes the feature label size since #1202).
const regionsTsv = (rows) => ['set_id\tid\tmark\tstart\tend\tlabel\tlabel_font_size',
  ...rows.map(([id, start, end, label]) => `regions\t${id}\thighlight\t${start}\t${end}\t${label}\t10`), ''].join('\n');
const annotationsPanel = (page) => page.locator('details').filter({ has: page.locator('summary[aria-label="Region Annotations"]') });
const generateLinear = async (page, name, features, regions, { separateStrands = true } = {}) => {
  test.setTimeout(600000);
  await page.setViewportSize({ width: 1440, height: 1000 });
  await openApp(page);
  await page.getByRole('button', { name: 'Linear', exact: true }).click();
  await page.waitForFunction(() => window.__GBDRAW_APP__?.mode === 'linear');
  await page.evaluate(async ({ text, file, separate }) => {
    const app = window.__GBDRAW_APP__;
    app.setLinearSeqPrimaryFile(0, 'gb', new File([text], file, { type: 'text/plain', lastModified: 1000 }));
    app.form.show_labels_linear = 'all';
    app.form.linear_track_layout = 'above';
    app.form.separate_strands = separate;
    await window.Vue.nextTick();
  }, { text: genbank(name, features), file: `${name}.gb`, separate: separateStrands });
  await annotationsPanel(page).locator('input[type=file]')
    .setInputFiles({ name: 'regions.tsv', mimeType: 'text/tab-separated-values', buffer: Buffer.from(regionsTsv(regions)) });
  await page.waitForFunction(() => window.__GBDRAW_APP__.annotationSets.length === 1);
  await generateAndWaitForResult(page);
};

// LP-D68, LP-D77, LP-D82 (OV-1155): two close regions whose long labels overlap
// each other. Nothing may move, so the figure is drawn with the overlap; the
// notice counts it and offers no "Choose a fix", as the report has no rows.
test('Linear region labels over other labels: the notice counts them and offers no fix', async ({ page }) => {
  await generateLinear(page, 'adjacent', [['3000..3300', 'photosystem II protein D1']],
    [['r0', 1000, 1100, 'First region with a long label'], ['r1', 1120, 1220, 'Second region with a long label']]);
  const notice = page.getByTestId('label-fit-notice');
  await expect(notice).toBeVisible();
  await expect(notice).toContainText('Region labels overlap other labels or lines (1).');
  await expect(notice.getByRole('button')).toHaveCount(0);
  expect(await page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    return state.labelFitReports.value[0].report.rows.length;
  })).toBe(0);
  await pageShots(page, 'linear-region-notice');
});

// LP-D87 (OV-1210): a Circular region label longer than its circle is drawn
// curved with what fits, its start and an ellipsis; the notice counts the cut
// label and offers no "Choose a fix".
test('a Circular region label longer than its circle: the notice counts it as cut and offers no fix', async ({ page }) => {
  test.setTimeout(600000);
  await page.setViewportSize({ width: 1440, height: 1000 });
  await openApp(page);
  await page.getByLabel('GenBank/DDBJ File', { exact: true }).setInputFiles({ name: 'cut.gb', mimeType: 'text/plain',
    buffer: Buffer.from(genbank('cut', [['3000..3300', 'photosystem II protein D1']])) });
  await page.waitForFunction(() => window.__GBDRAW_APP__.circularRecordList.length === 1);
  await annotationsPanel(page).locator('input[type=file]').setInputFiles({ name: 'regions.tsv',
    mimeType: 'text/tab-separated-values', buffer: Buffer.from(regionsTsv([['r0', 1000, 1100, 'a long region label '.repeat(40).trim()]])) });
  await page.waitForFunction(() => window.__GBDRAW_APP__.annotationSets.length === 1);
  await generateAndWaitForResult(page);
  const notice = page.getByTestId('label-fit-notice');
  await expect(notice).toBeVisible();
  await expect(notice).toContainText('Region labels are cut to fit their circle (1).');
  await expect(notice.getByRole('button')).toHaveCount(0);
  expect(await page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    return state.labelFitReports.value[0].report.rows.length;
  })).toBe(0);
  const shown = await page.evaluate(() => {
    const svg = new DOMParser().parseFromString(window.__GBDRAW_APP__.currentResults()[0].content, 'image/svg+xml');
    return svg.querySelector('[data-gbdraw-annotation-id] textPath')?.textContent ?? '';
  });
  expect(shown.startsWith('a long region label')).toBe(true);
  expect(shown.endsWith('\u2026')).toBe(true);
  await pageShots(page, 'circular-region-cut-notice');
});

// LP-D86 (OV-1208): in the tobacco Gallery figure, the GC ring set to 0.58 x R
// leaves the region labels drawn inside their ring no room; the figure keeps
// the set radii and the notice counts the four labels on the GC ring, with no
// "Choose a fix".
test('Circular region labels left no room by set radii: the notice counts them and offers no fix', async ({ page }) => {
  test.setTimeout(600000);
  await page.setViewportSize({ width: 1440, height: 1000 });
  await openApp(page);
  await page.locator('input[type=file][accept*="application/json"][accept*="application/gzip"]')
    .setInputFiles(join(process.cwd(), 'gbdraw/web/gallery/sessions/tobacco-chloroplast.gbdraw-session.json'));
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending, null, { timeout: 180000 });
  await page.evaluate(() => {
    window.__GBDRAW_APP__.adv.circular_track_slots.find((row) => row.id === 'gc_content').radius = '0.58';
  });
  await generateAndWaitForResult(page);
  const notice = page.getByTestId('label-fit-notice');
  await expect(notice).toBeVisible();
  await expect(notice).toContainText('Region labels lie on other tracks: the set track radii leave them no room (4).');
  await expect(notice.getByRole('button')).toHaveCount(0);
  expect(await page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    return state.labelFitReports.value[0].report.rows.length;
  })).toBe(0);
  await pageShots(page, 'circular-region-no-room-notice');
});

// LP-D82 (e'): on one strand the control region's label meets only the tRNA-Leu leader, so it
// is drawn in front of the leader lines on a white backing, outside its
// annotation's group. A click on it still selects that annotation's row.
test('a region label drawn in front of the leader lines selects its annotation', async ({ page }) => {
  await generateLinear(page, 'left_end',
    [['190..210', 'tRNA-Leu'], ['complement(3000..3300)', 'photosystem II protein D1']],
    [['r0', 1, 300, 'Mitochondrial control region']], { separateStrands: false });
  await expect(page.getByTestId('label-fit-notice')).toHaveCount(0);
  const front = page.locator('.origin-top svg [data-gbdraw-region-label-for]');
  await expect(front).toHaveCount(1);
  await expect(front.locator('rect')).toHaveAttribute('fill-opacity', '0.7');
  await annotationsPanel(page).evaluate((details) => { details.open = true; });
  await front.locator('text').click();
  expect(await page.evaluate(() => window.__GBDRAW_APP__.selectedAnnotation?.id)).toBe('r0');
  await expect(annotationsPanel(page).locator('.border-blue-400')).toHaveCount(1);
});
