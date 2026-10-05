// Design Q4 (docs/internal/web-gui-audit-20261004-override-precedence/
// DESIGN-Q4-FEATURE-IDENTITY.md, PR-Q4-4): popup Feature visibility, Label
// visibility, and label text edits name the original-source feature, so they
// keep naming it after crop, reverse complement, reordering, and duplication,
// and the Session replays the same diagram on the CLI. FINDINGS OV-01, OV-02,
// OV-04, OV-05, OV-11, OV-12.
const { test, expect } = require('@playwright/test');
const { readFileSync, writeFileSync } = require('node:fs');
const { gunzipSync } = require('node:zlib');
const { spawnSync } = require('node:child_process');
const path = require('node:path');
const { generateAndWaitForResult, evaluateWithRetainedPromise } = require('./helpers/app-lifecycle.cjs');
const { BATCH_FIXTURE, openFresh, settle } = require('./helpers/audit-browser.cjs');
const { download } = require('./helpers/mode-transition.cjs');

const records = (file) => readFileSync(file, 'utf8').split(/^\/\/\s*$/m)
  .map((chunk) => chunk.trim()).filter((chunk) => chunk.startsWith('LOCUS')).map((chunk) => `${chunk}\n//\n`);
const [TESTA, TESTB] = records(BATCH_FIXTURE);
// OV-01 repro B: the drawn hash of Y after the crop equals X's source hash.
const [COLLIDE] = records('tests/fixtures/b_collide.gb');
// OV-05: two CDS at the same coordinates (biological IDs <hash>~<source index>).
const SAME_COORDINATES = TESTA.replace(
  '     CDS             301..600\n',
  '     CDS             301..600\n                     /locus_tag="TWIN_A"\n'
    + '                     /product="twin product"\n     CDS             301..600\n'
);

const openLinear = async (page, texts, configure) => {
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
  }, texts);
  await settle(page);
  if (configure) await page.evaluate(configure);
  await settle(page);
};

const catalog = (page) => page.evaluate(() => window.__GBDRAW_APP__.extractedFeatures.map((feature) => ({
  svgId: feature.svg_id,
  record: feature.record_id,
  recordKey: feature.record_key,
  identity: feature.biological_feature_id,
  type: feature.type,
  start: feature.start,
  locusTag: feature.locus_tag || ''
})));

// One popup edit through the editor actions the controls call.
const edit = (page, svgId, change) => evaluateWithRetainedPromise(page, async ({ id, change: requested }) => {
  const app = window.__GBDRAW_APP__;
  const feature = app.extractedFeatures.find((item) => item.svg_id === id);
  if (!feature) throw new Error(`no feature ${id}`);
  await app.openFeatureEditorFromList(feature, null);
  if (requested.labelText !== undefined || requested.labelVisibility) {
    if (requested.labelText !== undefined) app.clickedFeature.labelText = requested.labelText;
    if (requested.labelVisibility) app.clickedFeature.labelVisibility = requested.labelVisibility;
    await app.updateClickedFeatureLabelText();
    if (app.hiddenLabelTextDialog?.show) await app.handleHiddenLabelTextChoice('keep');
  }
  if (requested.visibility) {
    app.clickedFeature.featureVisibility = requested.visibility;
    await app.updateClickedFeatureVisibility(requested.visibility);
    if (app.featureVisibilityScopeDialog.show) await app.handleFeatureVisibilityScopeChoice('feature');
  }
  app.clickedFeature = null;
}, { id: svgId, change });

// What the committed Result draws: feature paths and label texts by rendered ID.
const drawn = (page, ids) => page.evaluate((renderedIds) => {
  const app = window.__GBDRAW_APP__;
  const root = new DOMParser().parseFromString(app.results[app.selectedResultIndex].content, 'image/svg+xml')
    .documentElement;
  return Object.fromEntries(renderedIds.map((id) => {
    const paths = [...root.querySelectorAll(`[data-gbdraw-feature-id="${CSS.escape(id)}"]`)]
      .filter((element) => element.localName === 'path' && element.getAttribute('display') !== 'none');
    const labels = [...root.querySelectorAll(
      `[data-label-feature-id="${CSS.escape(id)}"] text, text[data-label-feature-id="${CSS.escape(id)}"]`
    )].filter((element) => !element.closest('[display="none"]'))
      .map((element) => element.textContent.trim());
    return [id, { drawn: paths.length > 0, labels }];
  }));
}, ids);

const labelTextsInResult = (page) => page.evaluate(() => {
  const app = window.__GBDRAW_APP__;
  const root = new DOMParser().parseFromString(app.results[app.selectedResultIndex].content, 'image/svg+xml')
    .documentElement;
  return [...root.querySelectorAll('[data-label-feature-id] text, text[data-label-feature-id]')]
    .map((element) => element.textContent.trim());
});

const replayOnCli = (sessionPath, outputPrefix) => {
  const result = spawnSync(process.env.GBDRAW_PYTHON || 'python', [
    '-c', 'import sys; from gbdraw.cli import main; sys.argv = ["gbdraw", *sys.argv[1:]]; main()',
    'linear', '--session', sessionPath, '-o', outputPrefix, '-f', 'svg'
  ], { encoding: 'utf8', env: { ...process.env } });
  expect(result.status, result.stdout + result.stderr).toBe(0);
  return readFileSync(`${outputPrefix}.svg`, 'utf8');
};

const cropTesta = () => {
  const app = window.__GBDRAW_APP__;
  app.setLinearRecordCrop(app.linearSeqs[0], 'region_start', 201);
  app.setLinearRecordCrop(app.linearSeqs[0], 'region_end', 3800);
  app.form.show_labels_linear = 'all';
  if (!app.adv.features.includes('misc_feature')) app.adv.features.push('misc_feature');
};

test('OV-01: edits on a cropped record reach their feature, also on CLI replay', async ({ page }, testInfo) => {
  test.setTimeout(300_000);
  await openLinear(page, [COLLIDE], cropTesta);
  await generateAndWaitForResult(page);
  await settle(page);
  const features = await catalog(page);
  const at = (type, start) => features.find((feature) => feature.type === type && feature.start === start).svgId;
  const [X, Y, TX, TY] = [at('misc_feature', 2400), at('misc_feature', 2600), at('tRNA', 2549), at('tRNA', 2749)];

  await edit(page, X, { visibility: 'off' });
  await settle(page);
  await edit(page, TX, { labelText: 'EDITED_X' });
  await settle(page);
  await generateAndWaitForResult(page);
  await settle(page);

  const result = await drawn(page, [X, Y, TX, TY]);
  expect(result[X].drawn).toBe(false);
  expect(result[Y].drawn).toBe(true);
  expect(result[TX].labels).toEqual(['EDITED_X']);
  expect(result[TY].labels).not.toContain('EDITED_X');

  const sessionPath = testInfo.outputPath('ov01.gbdraw-session.json');
  await download(page, 'Save Session', sessionPath);
  const replay = replayOnCli(sessionPath, testInfo.outputPath('ov01-replay'));
  expect(replay.match(/EDITED_X/g) || []).toHaveLength(1);
  expect(replay).toContain(`data-gbdraw-feature-id="${Y}"`);
  expect(replay).not.toContain(`data-gbdraw-feature-id="${X}"`);
});

test('OV-02: a record_location color rule on a cropped record paints the same feature live and on Generate', async ({ page }) => {
  test.setTimeout(300_000);
  await openLinear(page, [COLLIDE], cropTesta);
  await generateAndWaitForResult(page);
  await settle(page);
  const features = await catalog(page);
  const X = features.find((feature) => feature.type === 'misc_feature' && feature.start === 2400).svgId;
  const Y = features.find((feature) => feature.type === 'misc_feature' && feature.start === 2600).svgId;
  // The CLI's record_location is the drawn record's: X at 2401..2500 is drawn
  // at 2201..2300 after the crop that starts at 201.
  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    Object.assign(app.newSpecRule, {
      feat: 'misc_feature', qual: 'record_location', val: '^TESTA:2200\\.\\.2300:\\+$', color: '#00aa00', cap: 'drawn X'
    });
    return app.addSpecificRule();
  });
  await settle(page);
  const fill = (id, source) => page.evaluate(({ renderedId, from }) => {
    const app = window.__GBDRAW_APP__;
    const root = from === 'mounted' ? app.svgContainer.querySelector('svg')
      : new DOMParser().parseFromString(app.results[app.selectedResultIndex].content, 'image/svg+xml').documentElement;
    return [...root.querySelectorAll(`path[data-gbdraw-feature-id="${CSS.escape(renderedId)}"]`)]
      .map((element) => element.getAttribute('fill')).find((value) => value && value !== 'none') || null;
  }, { renderedId: id, from: source });
  expect((await fill(X, 'mounted'))?.toLowerCase()).toBe('#00aa00');
  expect((await fill(Y, 'mounted'))?.toLowerCase()).not.toBe('#00aa00');
  await generateAndWaitForResult(page);
  await settle(page);
  expect((await fill(X, 'result'))?.toLowerCase()).toBe('#00aa00');
  expect((await fill(Y, 'result'))?.toLowerCase()).not.toBe('#00aa00');
});

test('OV-04, OV-05: an edit of one copy or one same-coordinate feature reaches only it', async ({ page }) => {
  test.setTimeout(300_000);
  await openLinear(page, [SAME_COORDINATES, TESTA], () => {
    const app = window.__GBDRAW_APP__;
    app.form.show_labels_linear = 'none';
  });
  await generateAndWaitForResult(page);
  await settle(page);
  const features = await catalog(page);
  const [firstKey, secondKey] = [...new Set(features.map((feature) => feature.recordKey))];
  const copy = (recordKey, tag) => features.find((feature) => feature.recordKey === recordKey && feature.locusTag === tag);
  const twins = features.filter((feature) => feature.recordKey === firstKey && feature.type === 'CDS' && feature.start === 300);
  expect(twins).toHaveLength(2);
  expect(new Set(twins.map((feature) => feature.identity.split('~')[0])).size).toBe(1);
  expect(twins.every((feature) => /~\d+$/.test(feature.identity))).toBe(true);

  await edit(page, copy(firstKey, 'TESTA_0002').svgId, { visibility: 'off' });
  await settle(page);
  await edit(page, copy(firstKey, 'TESTA_0006').svgId, { labelVisibility: 'on', labelText: 'COPY1_ONLY' });
  await settle(page);
  const twin = twins[1];
  await edit(page, twin.svgId, { labelVisibility: 'on', labelText: 'TWIN_ONE' });
  await settle(page);
  await generateAndWaitForResult(page);
  await settle(page);

  const ids = [copy(firstKey, 'TESTA_0002').svgId, copy(secondKey, 'TESTA_0002').svgId];
  const result = await drawn(page, ids);
  expect(result[ids[0]].drawn).toBe(false);
  expect(result[ids[1]].drawn).toBe(true);
  const texts = await labelTextsInResult(page);
  expect(texts.filter((text) => text === 'COPY1_ONLY')).toHaveLength(1);
  expect(texts.filter((text) => text === 'TWIN_ONE')).toHaveLength(1);
  expect((await drawn(page, [twin.svgId]))[twin.svgId].labels).toEqual(['TWIN_ONE']);
});

test('OV-11: edits stay on their feature after the crop and reverse complement change', async ({ page }) => {
  test.setTimeout(300_000);
  await openLinear(page, [TESTA, TESTB], () => {
    const app = window.__GBDRAW_APP__;
    app.setLinearRecordCrop(app.linearSeqs[0], 'region_start', 201);
    app.setLinearRecordCrop(app.linearSeqs[0], 'region_end', 3800);
    app.linearSeqs[1].region_reverse = true;
    app.form.show_labels_linear = 'all';
    if (!app.adv.features.includes('misc_feature')) app.adv.features.push('misc_feature');
  });
  await generateAndWaitForResult(page);
  await settle(page);
  const pick = (features, record, type) => features.find((feature) => feature.record === record && feature.type === type);
  let features = await catalog(page);
  const identities = {
    hidden: pick(features, 'TESTA', 'misc_feature').identity,
    text: pick(features, 'TESTA', 'tRNA').identity,
    labelOff: pick(features, 'TESTB', 'tRNA').identity
  };
  await edit(page, pick(features, 'TESTA', 'misc_feature').svgId, { visibility: 'off' });
  await settle(page);
  await edit(page, pick(features, 'TESTA', 'tRNA').svgId, { labelText: 'PROBE_A_TRNA' });
  await settle(page);
  await edit(page, pick(features, 'TESTB', 'tRNA').svgId, { labelVisibility: 'off' });
  await settle(page);
  await generateAndWaitForResult(page);
  await settle(page);

  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    app.setLinearRecordCrop(app.linearSeqs[0], 'region_start', 101);
    app.linearSeqs[1].region_reverse = false;
  });
  await settle(page);
  await generateAndWaitForResult(page);
  await settle(page);
  features = await catalog(page);
  const byIdentity = (identity) => features.find((feature) => feature.identity === identity);
  expect(byIdentity(identities.hidden)).toBeUndefined();
  const textFeature = byIdentity(identities.text);
  const labelOffFeature = byIdentity(identities.labelOff);
  const result = await drawn(page, [textFeature.svgId, labelOffFeature.svgId]);
  expect(result[textFeature.svgId].labels).toEqual(['PROBE_A_TRNA']);
  expect(result[labelOffFeature.svgId].drawn).toBe(true);
  expect(result[labelOffFeature.svgId].labels).toEqual([]);
});

test('OV-12: after reordering records the popup shows the edit of the same feature', async ({ page }) => {
  test.setTimeout(300_000);
  await openLinear(page, [TESTA, TESTB], () => {
    window.__GBDRAW_APP__.form.show_labels_linear = 'all';
  });
  await generateAndWaitForResult(page);
  await settle(page);
  let features = await catalog(page);
  const target = features.find((feature) => feature.record === 'TESTB' && feature.type === 'tRNA');
  await edit(page, target.svgId, { labelVisibility: 'off' });
  await settle(page);
  await page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    await app.setLinearRecordLayoutEnabled(true);
    await window.Vue.nextTick();
    const uid = app.linearSeqs[0].uid;
    app.linearSeqs.forEach((seq) => app.setLinearRecordRow(seq.uid, 1));
    await window.Vue.nextTick();
    await app.moveLinearRecordWithinRow(uid, 1);
  });
  await settle(page);
  await generateAndWaitForResult(page);
  await settle(page);
  features = await catalog(page);
  const moved = features.find((feature) => feature.identity === target.identity && feature.recordKey === target.recordKey);
  expect(moved.svgId).not.toBe(target.svgId);
  const popup = await evaluateWithRetainedPromise(page, async (id) => {
    const app = window.__GBDRAW_APP__;
    const feature = app.extractedFeatures.find((item) => item.svg_id === id);
    await app.openFeatureEditorFromList(feature, null);
    const value = app.clickedFeature.labelVisibility;
    app.clickedFeature = null;
    return value;
  }, moved.svgId);
  expect(popup).toBe('off');
  expect((await drawn(page, [moved.svgId]))[moved.svgId].labels).toEqual([]);
});

// A Session before feature catalogs may restore no feature metadata until
// Generate, so the load waits for the committed Result only.
const loadSession = async (page, file) => {
  await page.locator('input[accept^=".json,"]').setInputFiles(file);
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending
    && window.__GBDRAW_APP__.results.length > 0, null, { timeout: 180_000 });
  await settle(page);
};

// Design Q4 4.3: Sessions 44 and 33 written by first-parent main keep their
// per-feature edits; Generate draws them on the same features.
test('a Session 44 with crop and reverse complement keeps its feature edits', async ({ page }) => {
  test.setTimeout(300_000);
  await openFresh(page);
  await loadSession(page, 'tests/fixtures/sessions/feature-edits-crop-rc.v44.gbdraw-session.json.gz');
  const rows = await page.evaluate(async () => {
    const { state } = await import('/gbdraw/web/js/state.js');
    return Object.values(state.featureOverrides).map((row) => [row.biologicalFeatureId, row.featureVisibility,
      row.labelVisibility, row.labelText]).sort();
  });
  expect(rows).toEqual([
    ['f571eb983', null, null, 'PROBE_A_TRNA'],
    ['f88047061', 'off', 'off', null],
    ['fb5977f81', 'off', null, null]
  ]);
  await generateAndWaitForResult(page);
  await settle(page);
  const features = await catalog(page);
  expect(features.some((feature) => feature.identity === 'fb5977f81')).toBe(false);
  expect(features.some((feature) => feature.identity === 'f88047061')).toBe(false);
  const trna = features.find((feature) => feature.identity === 'f571eb983');
  expect((await drawn(page, [trna.svgId]))[trna.svgId].labels).toEqual(['PROBE_A_TRNA']);
});

test('a Session 33 without a feature catalog keeps its feature edits', async ({ page }) => {
  test.setTimeout(300_000);
  await openFresh(page);
  await loadSession(page, 'tests/fixtures/sessions/feature-edits-circular.v33.gbdraw-session.json.gz');
  await generateAndWaitForResult(page);
  await settle(page);
  const rows = await page.evaluate(async () => {
    const { getCommittedCanonicalRenderRequest } = await import('/gbdraw/web/js/services/config.js');
    return getCommittedCanonicalRenderRequest().diagramOptions.featureOverrides;
  });
  const recordKey = rows[0]?.recordKey;
  expect(rows).toEqual([
    { recordKey, biologicalFeatureId: 'f406d90f1', featureVisibility: 'off', labelVisibility: null, labelText: null },
    { recordKey, biologicalFeatureId: 'f48b7bf2f', featureVisibility: null, labelVisibility: 'off', labelText: null },
    { recordKey, biologicalFeatureId: 'fbe3a7c0c', featureVisibility: null, labelVisibility: null, labelText: 'V33_LABEL' }
  ]);
  const features = await catalog(page);
  expect(features.some((feature) => feature.identity === 'f406d90f1')).toBe(false);
  expect(features.every((feature) => feature.recordKey === recordKey)).toBe(true);
  // The Session shows no labels (labels_mode none); text alone draws none.
  expect(await labelTextsInResult(page)).not.toContain('V33_LABEL');
});

const LINEAR_V33 = 'tests/fixtures/sessions/feature-edits-linear-crop-rc.v33.gbdraw-session.json.gz';

const featureOverrideRows = (page) => page.evaluate(async () => {
  const { state } = await import('/gbdraw/web/js/state.js');
  return Object.values(state.featureOverrides).map((row) => [row.recordKey, row.biologicalFeatureId,
    row.featureVisibility, row.labelVisibility, row.labelText]).sort();
});

const committedRecordKeys = (page) => page.evaluate(async () => {
  const { getCommittedCanonicalRenderRequest } = await import('/gbdraw/web/js/services/config.js');
  return getCommittedCanonicalRenderRequest().records.map((record) => record.recordKey);
});

// A Linear Session 33 keys each edit by `<drawn hash>_record_<n>`, and every
// input's features have record_idx 0. The edits keep their records, also on a
// cropped and a reverse-complemented record whose hidden features the Session's
// metadata no longer lists.
test('a Linear Session 33 with crop and reverse complement keeps its feature edits', async ({ page }) => {
  test.setTimeout(300_000);
  const alerts = [];
  page.on('dialog', (dialog) => alerts.push(dialog.message()));
  await openFresh(page);
  await loadSession(page, LINEAR_V33);
  const [testa, testb] = await committedRecordKeys(page);
  expect(alerts).toEqual(['Session loaded successfully!']);
  expect(await featureOverrideRows(page)).toEqual([
    [testa, 'fef810304', null, 'off', null],
    [testa, 'ffa1f4c4a', 'off', null, null],
    [testb, 'f3f7207f4', 'off', null, null],
    [testb, 'f88047061', null, null, 'V33_LINEAR_TRNA']
  ].sort());
  await generateAndWaitForResult(page);
  await settle(page);
  const features = await catalog(page);
  expect(features.some((feature) => feature.identity === 'ffa1f4c4a')).toBe(false);
  expect(features.some((feature) => feature.identity === 'f3f7207f4')).toBe(false);
  const trna = features.find((feature) => feature.recordKey === testb && feature.identity === 'f88047061');
  const unlabeled = features.find((feature) => feature.recordKey === testa && feature.identity === 'fef810304');
  const shown = await drawn(page, [trna.svgId, unlabeled.svgId]);
  expect(shown[trna.svgId].labels).toEqual(['V33_LINEAR_TRNA']);
  expect(shown[unlabeled.svgId]).toEqual({ drawn: true, labels: [] });
});

// The Linear Session 33 with one more edit that names no feature.
const writeUnmatchedEditSession = (testInfo) => {
  const session = JSON.parse(gunzipSync(readFileSync(LINEAR_V33)));
  session.features.featureVisibilityOverrides.f00000000_record_2 = 'off';
  const file = testInfo.outputPath('linear-v33-unmatched-edit.gbdraw-session.json');
  writeFileSync(file, JSON.stringify(session));
  return file;
};

// Design 4.3 rule 3: the Load notice counts exactly the edits it drops.
test('loading an older Session counts only the feature edits it drops', async ({ page }, testInfo) => {
  test.setTimeout(300_000);
  const file = writeUnmatchedEditSession(testInfo);
  const alerts = [];
  page.on('dialog', (dialog) => alerts.push(dialog.message()));
  await openFresh(page);
  await loadSession(page, file);
  expect(alerts).toEqual([
    'Session loaded successfully! 1 feature edit(s) from an older Session could not be matched to a feature of its saved diagram and were dropped.'
  ]);
  expect((await featureOverrideRows(page)).length).toBe(4);
});

// OV-39: an older Session without a feature catalog reads its sources again
// through the diagram Worker. When the Worker cannot start, the Load fails with
// the runtime diagnostic and keeps the previous Session; it does not drop the
// older Session's edits as unmatched.
test('an older Session whose sources cannot be read again fails to load and keeps the previous Session', async ({ page }, testInfo) => {
  test.setTimeout(300_000);
  const file = writeUnmatchedEditSession(testInfo);
  const alerts = [];
  page.on('dialog', (dialog) => alerts.push(dialog.message()));
  await openFresh(page);
  const runtimeAssets = /pyodide\.asm\.(?:js|wasm)(?:\?|$)/;
  await page.context().route(runtimeAssets, (route) => route.abort());
  // A current Session loads without the diagram Worker.
  await loadSession(page, 'tests/fixtures/sessions/feature-edits-crop-rc.v44.gbdraw-session.json.gz');
  const loaded = () => page.evaluate(async () => {
    const app = window.__GBDRAW_APP__;
    const { state } = await import('/gbdraw/web/js/state.js');
    return { mode: app.mode, results: app.results.map((result) => result.content),
      rows: Object.values(state.featureOverrides).map((row) => [row.biologicalFeatureId, row.featureVisibility,
        row.labelVisibility, row.labelText]).sort(), records: app.linearSeqs.map((seq) => seq.gb?.name || null) };
  });
  const previous = await loaded();
  expect(previous.rows.length).toBe(3);
  alerts.length = 0;
  await page.locator('input[accept^=".json,"]').setInputFiles(file);
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending, null, { timeout: 180_000 });
  await settle(page);
  expect(alerts).toEqual([]);
  expect(await page.evaluate(async () => {
    const { state } = await import('/gbdraw/web/js/state.js');
    const { code, stage, actions } = state.errorLog.value || {};
    return { code, stage, actions };
  })).toEqual({ code: 'WORKER_INIT', stage: 'initialization', actions: ['save-session', 'reload', 'retry'] });
  await expect(page.getByRole('alert')).toContainText('The diagram runtime could not start.');
  expect(await loaded()).toEqual(previous);
  expect(await page.evaluate(() => [window.__GBDRAW_APP__.sessionSaveAvailable,
    window.__GBDRAW_APP__.sessionLoadAvailable])).toEqual([true, true]);
  // Retry once the runtime can start: the Load drops only the unmatched edit.
  await page.context().unroute(runtimeAssets);
  await loadSession(page, file);
  expect(alerts).toEqual([
    'Session loaded successfully! 1 feature edit(s) from an older Session could not be matched to a feature of its saved diagram and were dropped.'
  ]);
  expect((await featureOverrideRows(page)).length).toBe(4);
  await expect(page.getByRole('alert')).toHaveCount(0);
});

const importLabelTsv = (page, text) => evaluateWithRetainedPromise(page, async (tsv) => {
  await window.__GBDRAW_APP__.loadLabelOverrideTable({
    target: { files: [new File([tsv], 'labels.tsv')], value: 'labels.tsv' }
  });
}, text);

const previewLabelCount = (page, text) => page.locator('.origin-top svg text').filter({ hasText: text }).count();

// A Result of feature catalog 3 or 4 (Session 44, the released Gallery Session
// fixture) has
// no drawn selector values. A rendered ID that carries the source hash means
// the record was drawn with its source coordinates, so a `location` or
// `record_location` row matches the source values there.
test('a Label TSV record_location row applies to a catalog 4 Result drawn with source coordinates', async ({ page }) => {
  test.setTimeout(300_000);
  const alerts = [];
  page.on('dialog', (dialog) => alerts.push(dialog.message()));
  await openFresh(page);
  await loadSession(page, 'tests/fixtures/sessions/lambda_basic_linear.v44-schema8.gbdraw-session.json.gz');
  alerts.length = 0;
  await importLabelTsv(page, 'NC_001416.1\tCDS\trecord_location\t^NC_001416\\.1:190\\.\\.736:\\+$\tTSV_RL\n');
  await settle(page);
  expect(alerts).toEqual(['Loaded 1 row(s). Applied to 1 label(s).']);
  expect(await previewLabelCount(page, /^TSV_RL$/)).toBe(1);
});

// Where the drawn values are unknown (a cropped or reverse-complemented record
// of a catalog 4 Result), the import is declined with its reason and keeps the
// label edits it would have replaced; after Generate it applies.
test('a Label TSV location row waits for Generate on a catalog 4 Result of cropped records', async ({ page }) => {
  test.setTimeout(300_000);
  const alerts = [];
  page.on('dialog', (dialog) => alerts.push(dialog.message()));
  await openFresh(page);
  await loadSession(page, 'tests/fixtures/sessions/feature-edits-crop-rc.v44.gbdraw-session.json.gz');
  const before = await featureOverrideRows(page);
  alerts.length = 0;
  const row = '*\ttRNA\trecord_location\t^TESTA:\tTSV_TRNA\n';
  await importLabelTsv(page, row);
  await settle(page);
  expect(alerts).toHaveLength(1);
  expect(alerts[0]).toMatch(/^Loaded 1 row\(s\)\. Not applied: .*Generate/);
  expect(await featureOverrideRows(page)).toEqual(before);
  expect(await previewLabelCount(page, /^PROBE_A_TRNA$/)).toBe(1);
  await generateAndWaitForResult(page);
  await settle(page);
  alerts.length = 0;
  await importLabelTsv(page, row);
  await settle(page);
  expect(alerts).toEqual(['Loaded 1 row(s). Applied to 1 label(s).']);
  expect(await previewLabelCount(page, /^TSV_TRNA$/)).toBe(1);
});

// Owner decision Q1 = A: a Session 44 Feature visibility edit was a `hash`
// row that hid the feature in every copy of a duplicated record; it now names
// the copy that was edited, and Load says the next Generate draws the others.
test('loading a Session 44 says which Feature visibility edits now apply to one copy', async ({ page }) => {
  test.setTimeout(300_000);
  const alerts = [];
  page.on('dialog', (dialog) => alerts.push(dialog.message()));
  await openFresh(page);
  await loadSession(page, 'tests/fixtures/sessions/feature-edits-circular-copies.v44.gbdraw-session.json.gz');
  expect(alerts).toEqual([
    'Session loaded successfully! 1 Feature visibility edit(s) from an older Session hid every feature with '
      + 'the same hash, such as each copy of a duplicated record. Each now applies only to the feature that '
      + 'was edited, so the next Generate draws the others.'
  ]);
});
