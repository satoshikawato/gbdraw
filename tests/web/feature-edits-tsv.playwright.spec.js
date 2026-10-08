// Design Q4 6.4 (Owner Q2 = A): Export Feature Edits TSV writes the per-feature
// Feature visibility, Label visibility, and label text edits as
// --feature_override_table rows, and Load Feature Edits TSV reads them back
// through Python's reader (R4) as one Undo step. Rows that name no record or
// feature of the current diagram are counted, not applied (Owner Q3 = A); a
// malformed table is a classified table error (R6).
const { test, expect } = require('@playwright/test');
const { readFileSync, writeFileSync } = require('node:fs');
const { spawnSync } = require('node:child_process');
const path = require('node:path');
const { generateAndWaitForResult, evaluateWithRetainedPromise } = require('./helpers/app-lifecycle.cjs');
const { BATCH_FIXTURE, openFresh, settle } = require('./helpers/audit-browser.cjs');

test.describe.configure({ retries: 0 });

const [TESTA, TESTB] = readFileSync(BATCH_FIXTURE, 'utf8').split(/^\/\/\s*$/m)
  .map((chunk) => chunk.trim()).filter((chunk) => chunk.startsWith('LOCUS')).map((chunk) => `${chunk}\n//\n`);
const HEADER = 'record\tfeature_selector\tfeature_visibility\tlabel_visibility\tlabel_text';
const replayEnv = { ...process.env };
delete replayEnv.PYTHONPATH;
delete replayEnv.PYTHONHOME;

// TESTA cropped, TESTB reverse-complemented, and a second copy of TESTA.
const INPUTS = [TESTA, TESTB, TESTA];
const openLinear = async (page) => {
  await openFresh(page);
  await page.getByRole('button', { name: 'Linear', exact: true }).click();
  await page.waitForFunction(() => window.__GBDRAW_APP__?.mode === 'linear');
  await page.evaluate(async (inputs) => {
    const app = window.__GBDRAW_APP__;
    while (app.linearSeqs.length < inputs.length) app.addLinearSeq();
    inputs.forEach((text, index) => app.setLinearSeqPrimaryFile(index, 'gb', new File(
      [text], `record-${index}.gb`, { type: 'text/plain', lastModified: 1000 + index }
    )));
    await window.Vue.nextTick();
  }, INPUTS);
  await settle(page);
  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    app.setLinearRecordCrop(app.linearSeqs[0], 'region_start', 201);
    app.setLinearRecordCrop(app.linearSeqs[0], 'region_end', 3800);
    app.linearSeqs[1].region_reverse = true;
    app.form.show_labels_linear = 'all';
    if (!app.adv.features.includes('misc_feature')) app.adv.features.push('misc_feature');
  });
  await settle(page);
  await generateAndWaitForResult(page);
  await settle(page);
};

// The displayed features by record position (#1-based) and type.
const features = (page) => page.evaluate(async () => {
  const { getCommittedCanonicalRenderRequest } = await import('/gbdraw/web/js/services/config.js');
  const keys = getCommittedCanonicalRenderRequest().records.map((record) => record.recordKey);
  return window.__GBDRAW_APP__.extractedFeatures.map((feature) => ({
    svgId: feature.svg_id,
    record: keys.indexOf(feature.record_key) + 1,
    identity: feature.biological_feature_id,
    type: feature.type
  }));
});
const pick = (list, record, type) => list.find((feature) => feature.record === record && feature.type === type);

// One popup edit through the editor actions the controls call.
const edit = (page, svgId, change) => evaluateWithRetainedPromise(page, async ({ id, change: requested }) => {
  const app = window.__GBDRAW_APP__;
  const feature = app.extractedFeatures.find((item) => item.svg_id === id);
  await app.openFeatureEditorFromList(feature, null);
  if (requested.labelText !== undefined || requested.labelVisibility) {
    if (requested.labelText !== undefined) app.clickedFeature.labelText = requested.labelText;
    if (requested.labelVisibility) app.clickedFeature.labelVisibility = requested.labelVisibility;
    // Label Not Shown keeps the Apply open until its choice (Keep hidden here).
    let applying = true;
    const applied = Promise.resolve(app.updateClickedFeatureLabelText()).finally(() => { applying = false; });
    while (applying && !app.hiddenLabelTextDialog?.show) await new Promise((resolve) => setTimeout(resolve, 20));
    if (app.hiddenLabelTextDialog?.show) await app.handleHiddenLabelTextChoice('text_only');
    await applied;
  }
  if (requested.visibility) {
    app.clickedFeature.featureVisibility = requested.visibility;
    await app.updateClickedFeatureVisibility(requested.visibility);
    if (app.featureVisibilityScopeDialog.show) await app.handleFeatureVisibilityScopeChoice('feature');
  }
  app.clickedFeature = null;
}, { id: svgId, change });

// What a Result draws: visible feature paths and visible label texts by rendered ID.
const drawn = (page, ids, from = 'result') => page.evaluate(({ renderedIds, source }) => {
  const app = window.__GBDRAW_APP__;
  const root = source === 'mounted' ? app.svgContainer.querySelector('svg')
    : new DOMParser().parseFromString(app.results[app.selectedResultIndex].content, 'image/svg+xml').documentElement;
  const hidden = (element) => Boolean(element.closest('[display="none"]'));
  return Object.fromEntries(renderedIds.map((id) => {
    const paths = [...root.querySelectorAll(`path[data-gbdraw-feature-id="${CSS.escape(id)}"]`)]
      .filter((element) => !hidden(element));
    const labels = [...root.querySelectorAll(
      `[data-label-feature-id="${CSS.escape(id)}"] text, text[data-label-feature-id="${CSS.escape(id)}"]`
    )].filter((element) => !hidden(element)).map((element) => element.textContent.trim());
    return [id, { drawn: paths.length > 0, labels }];
  }));
}, { renderedIds: ids, source: from });

const draftRows = (page) => page.evaluate(async () => {
  const { getCommittedCanonicalRenderRequest } = await import('/gbdraw/web/js/services/config.js');
  const { state } = await import('/gbdraw/web/js/state.js');
  const keys = getCommittedCanonicalRenderRequest().records.map((record) => record.recordKey);
  return Object.values(state.featureOverrides).map((row) => [keys.indexOf(row.recordKey) + 1,
    row.biologicalFeatureId, row.featureVisibility, row.labelVisibility, row.labelText]).sort();
});

const exportTable = async (page, file) => {
  await page.evaluate(() => window.__GBDRAW_APP__.openRightDrawerTab('features'));
  const [download] = await Promise.all([
    page.waitForEvent('download', { timeout: 60_000 }),
    page.getByRole('button', { name: 'Export Feature Edits TSV', exact: true }).click()
  ]);
  await download.saveAs(file);
  return readFileSync(file, 'utf8');
};

const loadTable = async (page, file, alerts) => {
  const expected = alerts.length + 1;
  await page.evaluate(() => window.__GBDRAW_APP__.openRightDrawerTab('features'));
  const [chooser] = await Promise.all([
    page.waitForEvent('filechooser'),
    page.getByRole('button', { name: 'Load Feature Edits TSV', exact: true }).click()
  ]);
  await chooser.setFiles(file);
  await expect.poll(() => alerts.length, { timeout: 120_000 }).toBe(expected);
  await settle(page);
  return alerts.at(-1);
};

// The Source recipe of a Result drawn without feature edits reads only the
// input files, so the CLI run needs no reproducibility bundle.
const sourceRecipe = async (page) => {
  const recipe = await page.evaluate(() => {
    const source = window.__GBDRAW_APP__.lastRunInfo?.sourceRecipe || {};
    return { available: source.available, command: source.command || '', reason: source.unavailableReason || '',
      generated: (source.generatedFiles || []).length };
  });
  expect(recipe.available, recipe.reason).toBe(true);
  expect(recipe.generated, recipe.command).toBe(0);
  return recipe.command;
};

const resultSvg = async (page, file) => {
  writeFileSync(file, await page.evaluate(() => window.__GBDRAW_APP__.results[0].content));
  return file;
};

test('Export and Load Feature Edits TSV carry per-feature edits by identity, also to the CLI', async ({ page }, testInfo) => {
  test.setTimeout(600_000);
  const alerts = [];
  page.on('dialog', (dialog) => alerts.push(dialog.message()));
  await openLinear(page);
  const baseline = await sourceRecipe(page);
  const baselineSvg = await resultSvg(page, testInfo.outputPath('baseline.gui.svg'));

  let list = await features(page);
  const hidden = pick(list, 3, 'misc_feature');
  const unlabeled = pick(list, 2, 'tRNA');
  const renamed = pick(list, 1, 'tRNA');
  await edit(page, hidden.svgId, { visibility: 'off' });
  await settle(page);
  await edit(page, unlabeled.svgId, { labelVisibility: 'off' });
  await settle(page);
  await edit(page, renamed.svgId, { labelText: 'EDITED_TRNA' });
  await settle(page);
  const edits = await draftRows(page);
  expect(edits).toEqual([
    [1, renamed.identity, null, null, 'EDITED_TRNA'],
    [2, unlabeled.identity, null, 'off', null],
    [3, hidden.identity, 'off', null, null]
  ]);

  // The rows of --feature_override_table: #index of the request record and the
  // feature's biological ID (a duplicated record is told apart by its #index).
  const table = testInfo.outputPath('edits.feature_override_table.tsv');
  const exported = (await exportTable(page, table)).split('\n');
  expect(exported[0]).toBe(HEADER);
  expect(exported.at(-1)).toBe('');
  expect(exported.slice(1, -1).sort()).toEqual([
    `#1\thash=${renamed.identity}\t\t\tEDITED_TRNA`,
    `#2\thash=${unlabeled.identity}\t\toff\t`,
    `#3\thash=${hidden.identity}\toff\t\t`
  ].sort());

  await generateAndWaitForResult(page);
  await settle(page);
  const editedSvg = await resultSvg(page, testInfo.outputPath('edited.gui.svg'));
  // Generate does not draw the hidden feature. Loading the table changes the
  // Legend source (the first misc_feature is the hidden one), so the automatic
  // rerender draws the preview as Generate does, without the feature (OV-42,
  // #857).
  const expectDrawsEdits = async (target, from) => {
    await settle(target);
    list = await features(target);
    const hiddenCopy = pick(list, 3, 'misc_feature');
    expect(hiddenCopy).toBeUndefined();
    const ids = [pick(list, 1, 'misc_feature'), pick(list, 2, 'tRNA'), pick(list, 1, 'tRNA'), pick(list, 3, 'tRNA')]
      .map((feature) => feature.svgId);
    const shown = await drawn(target, ids, from);
    expect(ids.map((id) => shown[id])).toEqual([
      { drawn: true, labels: expect.any(Array) },
      { drawn: true, labels: [] },
      { drawn: true, labels: ['EDITED_TRNA'] },
      { drawn: true, labels: ['tRNA-Leu'] }
    ]);
  };
  await expectDrawsEdits(page, 'result');

  // The baseline Source recipe with the exported table draws the edited Result.
  const replay = spawnSync(process.env.GBDRAW_PYTHON || 'python', [
    path.resolve('tests/web/helpers/source-recipe-replay.py'), testInfo.outputPath('cli'),
    `${baseline} --feature_override_table ${path.basename(table)}`, editedSvg, baselineSvg,
    '-', ...INPUTS.map((text, index) => {
      const file = testInfo.outputPath(`record-${index}.gb`);
      writeFileSync(file, text);
      return file;
    }), table
  ], { encoding: 'utf8', env: replayEnv });
  expect(replay.status, `${baseline}\n${replay.stdout}${replay.stderr}`).toBe(0);

  // A new page has new record keys; the table names records by #index.
  const fresh = await page.context().newPage();
  fresh.on('dialog', (dialog) => alerts.push(dialog.message()));
  await openLinear(fresh);
  expect(await draftRows(fresh)).toEqual([]);
  alerts.length = 0;
  expect(await loadTable(fresh, table, alerts)).toBe('Loaded 3 row(s). Applied 3 feature edit(s).');
  expect(await draftRows(fresh)).toEqual(edits);
  await expectDrawsEdits(fresh, 'mounted');

  await evaluateWithRetainedPromise(fresh, () => window.__GBDRAW_HISTORY__.undo());
  await settle(fresh);
  expect(await draftRows(fresh)).toEqual([]);
  list = await features(fresh);
  const trna = pick(list, 1, 'tRNA').svgId;
  expect((await drawn(fresh, [trna], 'mounted'))[trna].labels).toEqual(['tRNA-Leu']);
  await evaluateWithRetainedPromise(fresh, () => window.__GBDRAW_HISTORY__.redo());
  await settle(fresh);
  expect(await draftRows(fresh)).toEqual(edits);

  await generateAndWaitForResult(fresh);
  await settle(fresh);
  await expectDrawsEdits(fresh, 'result');
});

test('Load Feature Edits TSV counts rows that name no record or feature and rejects a malformed table', async ({ page }, testInfo) => {
  test.setTimeout(300_000);
  const alerts = [];
  page.on('dialog', (dialog) => alerts.push(dialog.message()));
  await openLinear(page);
  await page.evaluate(() => window.__GBDRAW_APP__.openRightDrawerTab('features'));
  await page.getByRole('button', { name: 'Export Feature Edits TSV', exact: true }).click();
  await expect.poll(() => alerts.at(-1)).toBe('No feature edits to export.');

  const list = await features(page);
  const target = pick(list, 2, 'tRNA');
  const table = testInfo.outputPath('unmatched.tsv');
  writeFileSync(table, [
    HEADER,
    `#2\thash=${target.identity}\toff\t\t`,
    `#4\thash=${target.identity}\toff\t\t`,
    '#1\thash=f00000000\toff\t\t',
    ''
  ].join('\n'));
  alerts.length = 0;
  expect(await loadTable(page, table, alerts)).toBe('Loaded 3 row(s). Applied 1 feature edit(s). '
    + '2 row(s) name no record or feature of the current diagram and were not applied.');
  expect(await draftRows(page)).toEqual([[2, target.identity, 'off', null, null]]);

  const malformed = testInfo.outputPath('malformed.tsv');
  writeFileSync(malformed, `${HEADER}\n#2\thash=${target.identity}\thidden\t\t\n`);
  await page.evaluate(() => window.__GBDRAW_APP__.openRightDrawerTab('features'));
  const [chooser] = await Promise.all([
    page.waitForEvent('filechooser'),
    page.getByRole('button', { name: 'Load Feature Edits TSV', exact: true }).click()
  ]);
  await chooser.setFiles(malformed);
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.errorLog?.code || ''), { timeout: 120_000 })
    .toBe('TABLE_INVALID');
  expect(await page.evaluate(() => {
    const { operation, context } = window.__GBDRAW_APP__.errorLog;
    return { operation, row: context?.row };
  })).toEqual({ operation: 'readFeatureOverrideTable', row: 2 });
  await expect(page.getByRole('button', { name: 'Reselect Feature Edits TSV', exact: true })).toBeVisible();
  expect(await draftRows(page)).toEqual([[2, target.identity, 'off', null, null]]);
});
