// Design Q4 PR-Q4-5 (docs/internal/web-gui-audit-20261004-override-precedence/
// DESIGN-Q4-FEATURE-IDENTITY.md, section 7): an annotation made from selected
// features names each feature by its source identity, so it stays on that
// feature after crop changes, in one copy of a duplicated record, through Save
// and Load, and on CLI replay of the Session. FINDINGS OV-03.
const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { spawnSync } = require('node:child_process');
const { generateAndWaitForResult } = require('./helpers/app-lifecycle.cjs');
const { BATCH_FIXTURE, loadSessionFile, openFresh, openWithGenBank, settle } = require('./helpers/audit-browser.cjs');
const { download } = require('./helpers/mode-transition.cjs');

const records = (file) => readFileSync(file, 'utf8').split(/^\/\/\s*$/m)
  .map((chunk) => chunk.trim()).filter((chunk) => chunk.startsWith('LOCUS')).map((chunk) => `${chunk}\n//\n`);
const [TESTA] = records(BATCH_FIXTURE);
// OV-01 repro B: after the crop, Y (misc_feature 2601..2700) is drawn with the
// hash that X (misc_feature 2401..2500) has in the source.
const [COLLIDE] = records('tests/fixtures/b_collide.gb');
// Two CDS at the same coordinates (biological IDs <hash>~<source index>).
const TWINS = TESTA.replace(
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

const cropTesta = () => {
  const app = window.__GBDRAW_APP__;
  app.setLinearRecordCrop(app.linearSeqs[0], 'region_start', 201);
  app.setLinearRecordCrop(app.linearSeqs[0], 'region_end', 3800);
  if (!app.adv.features.includes('misc_feature')) app.adv.features.push('misc_feature');
};

const generate = async (page) => {
  await generateAndWaitForResult(page);
  await settle(page);
};

const catalog = (page) => page.evaluate(() => window.__GBDRAW_APP__.extractedFeatures.map((feature) => ({
  svgId: feature.svg_id,
  recordKey: feature.record_key,
  identity: feature.biological_feature_id,
  drawnHash: feature.drawnSelector?.hash || null,
  type: feature.type,
  start: feature.start
})));

const panel = (page) => page.locator('details')
  .filter({ has: page.locator('summary[aria-label="Region Annotations"]') });

// Ctrl-click selects each feature in the preview; "Selected features" makes
// one annotation per selected feature in a new set.
const annotateSelection = async (page, svgIds) => {
  for (const id of svgIds) {
    await page.locator(`.origin-top svg [data-gbdraw-rendered-feature-id="${id}"], `
      + `.origin-top svg [data-gbdraw-feature-id="${id}"]`).first().dispatchEvent('click', { ctrlKey: true });
  }
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.selectedFeatureCount)).toBe(svgIds.length);
  const annotations = panel(page);
  if (await annotations.getAttribute('open') === null) await annotations.locator('summary').press('Enter');
  await annotations.getByRole('button', { name: /Add set/ }).click();
  await annotations.getByRole('button', { name: /Selected features/ }).click();
  await settle(page);
};

const annotationTargets = (page) => page.evaluate(() => JSON.parse(JSON.stringify(
  window.__GBDRAW_APP__.annotationSets.flatMap((set) => set.annotations.map((item) => item.target))
)));

// The annotations the committed Result draws, from its feature catalog:
// record position and drawn 0-based half-open segments.
const drawnAnnotations = (page) => page.evaluate(async () => {
  const { state } = await import('/gbdraw/web/js/state.js');
  const { admitFeatureCatalog } = await import('/gbdraw/web/js/services/feature-catalog.js');
  const raw = (value) => window.Vue.toRaw(value);
  return admitFeatureCatalog(raw(state.featureCatalog.value), raw(state.results.value), {
    mode: state.generatedMode.value
  }).annotations.map((item) => ({ id: item.id, recordIndex: item.record_index, segments: item.segments }));
});

const annotationWarnings = (page) => page.evaluate(() => JSON.parse(JSON.stringify(
  window.__GBDRAW_APP__.annotationWarnings
)));

test('OV-03: an annotation of a selected feature on a cropped record stays on it after the crop changes', async ({ page }) => {
  test.setTimeout(300_000);
  await openLinear(page, [COLLIDE], cropTesta);
  await generate(page);
  const features = await catalog(page);
  const X = features.find((feature) => feature.type === 'misc_feature' && feature.start === 2400);
  await annotateSelection(page, [X.svgId]);
  await generate(page);

  // X is drawn at 2201..2300 after the crop that starts at 201 (Y at 2401..2500).
  expect(await drawnAnnotations(page)).toEqual([{ id: 'feature_1', recordIndex: 0, segments: [[2200, 2300]] }]);
  expect(await annotationWarnings(page)).toEqual([]);
  expect(await annotationTargets(page)).toEqual([{
    kind: 'featureIdentity', recordKey: X.recordKey, biologicalFeatureId: X.identity,
    envelope: 'outer_bounds', circularPath: 'shortest'
  }]);
  await expect(panel(page).locator('[data-annotation-feature-identity]'))
    .toHaveText('Selected feature: misc_feature at 2401..2500');

  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    app.setLinearRecordCrop(app.linearSeqs[0], 'region_start', 101);
  });
  await settle(page);
  await generate(page);
  expect(await drawnAnnotations(page)).toEqual([{ id: 'feature_1', recordIndex: 0, segments: [[2300, 2400]] }]);
  expect(await annotationWarnings(page)).toEqual([]);
});

// Design 6.4: the annotation table has no source-identity selector, so the TSV
// names each selected feature by its current record position and drawn hash.
test('the annotation TSV writes a selected feature by its current record position and drawn hash', async ({ page }, testInfo) => {
  test.setTimeout(300_000);
  const alerts = [];
  page.on('dialog', (dialog) => alerts.push(dialog.message()));
  await openLinear(page, [COLLIDE], cropTesta);
  await generate(page);
  const features = await catalog(page);
  const X = features.find((feature) => feature.type === 'misc_feature' && feature.start === 2400);
  const Y = features.find((feature) => feature.type === 'misc_feature' && feature.start === 2600);
  // The drawn hash of Y is the source hash of X (OV-01 repro B).
  expect(Y.drawnHash).toBe(X.identity);
  await annotateSelection(page, [X.svgId]);
  await generate(page);

  const tsvPath = testInfo.outputPath('annotations.tsv');
  alerts.length = 0;
  await download(page, 'Download TSV', tsvPath);
  const [header, row, ...rest] = readFileSync(tsvPath, 'utf8').trim().split('\n').map((line) => line.split('\t'));
  expect(rest).toEqual([]);
  const cell = (name) => row[header.indexOf(name)];
  expect([cell('record'), cell('feature_selector'), cell('envelope'), cell('circular_path')])
    .toEqual(['#1', `hash=${X.drawnHash}`, 'outer_bounds', 'shortest']);
  expect(alerts).toEqual([
    '1 annotation(s) of selected features were written as record=#<position> and feature_selector=hash=<drawn hash>. '
      + 'They name the same features only while the crop, orientation, and record order stay as drawn.'
  ]);

  // The table imported back draws the same annotation in the same placement.
  await panel(page).locator('input[type=file]').setInputFiles(tsvPath);
  await expect.poll(() => annotationTargets(page)).toMatchObject([{
    kind: 'featureSpan', selectors: [{ key: 'hash', value: X.drawnHash }]
  }]);
  await generate(page);
  expect(await drawnAnnotations(page)).toEqual([{ id: 'feature_1', recordIndex: 0, segments: [[2200, 2300]] }]);
});

const replayOnCli = (sessionPath, outputPrefix) => {
  const result = spawnSync(process.env.GBDRAW_PYTHON || 'python', [
    '-c', 'import sys; from gbdraw.cli import main; sys.argv = ["gbdraw", *sys.argv[1:]]; main()',
    'linear', '--session', sessionPath, '-o', outputPrefix, '-f', 'svg'
  ], { encoding: 'utf8', env: { ...process.env } });
  expect(result.status, result.stdout + result.stderr).toBe(0);
  return readFileSync(`${outputPrefix}.svg`, 'utf8');
};

// The drawn shapes of one annotation in an SVG document.
const annotationShapes = (page, svg, id) => page.evaluate(({ text, annotationId }) => {
  const root = new DOMParser().parseFromString(text, 'image/svg+xml').documentElement;
  return [...root.querySelectorAll(`[data-gbdraw-annotation-id="${annotationId}"] rect`)]
    .map((element) => ['x', 'y', 'width', 'height'].map((name) => element.getAttribute(name)));
}, { text: svg, annotationId: id });

test('an annotation of one same-coordinate feature in one copy of a record stays on it through Save, Load, and CLI replay', async ({ page, browser }, testInfo) => {
  test.setTimeout(300_000);
  await openLinear(page, [TWINS, TWINS]);
  await generate(page);
  const features = await catalog(page);
  const secondKey = [...new Set(features.map((feature) => feature.recordKey))][1];
  const twins = features.filter((feature) => feature.recordKey === secondKey
    && feature.type === 'CDS' && feature.start === 300);
  expect(twins).toHaveLength(2);
  const twin = twins.find((feature) => /~1$/.test(feature.identity));
  await annotateSelection(page, [twin.svgId]);
  await generate(page);
  const target = {
    kind: 'featureIdentity', recordKey: secondKey, biologicalFeatureId: twin.identity,
    envelope: 'outer_bounds', circularPath: 'shortest'
  };
  expect(await annotationTargets(page)).toEqual([target]);
  const drawn = [{ id: 'feature_1', recordIndex: 1, segments: [[300, 600]] }];
  expect(await drawnAnnotations(page)).toEqual(drawn);
  const webSvg = await page.evaluate(() => window.__GBDRAW_APP__.results[0].content);
  const shapes = await annotationShapes(page, webSvg, 'feature_1');
  expect(shapes).toHaveLength(1);

  const sessionPath = testInfo.outputPath('twins.gbdraw-session.json');
  await download(page, 'Save Session', sessionPath);
  expect(await annotationShapes(page, replayOnCli(sessionPath, testInfo.outputPath('twins-replay')), 'feature_1'))
    .toEqual(shapes);

  const context = await browser.newContext({ baseURL: new URL(page.url()).origin });
  try {
    const fresh = await context.newPage();
    await openFresh(fresh);
    await fresh.locator('input[accept^=".json,"]').setInputFiles(sessionPath);
    await fresh.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending
      && window.__GBDRAW_APP__.results.length > 0, null, { timeout: 180_000 });
    await settle(fresh);
    expect(await annotationTargets(fresh)).toEqual([target]);
    const caption = await fresh.evaluate(async (identity) => {
      const { getFeatureCaption } = await import('/gbdraw/web/js/services/feature-utils.js');
      return getFeatureCaption(window.__GBDRAW_APP__.extractedFeatures
        .find((feature) => feature.biological_feature_id === identity.id && feature.record_key === identity.key));
    }, { id: twin.identity, key: secondKey });
    const annotations = panel(fresh);
    await annotations.locator('summary').press('Enter');
    await expect(annotations.locator('[data-annotation-feature-identity]')).toHaveText(`Selected feature: ${caption}`);
    await generate(fresh);
    expect(await drawnAnnotations(fresh)).toEqual(drawn);
  } finally {
    await context.close();
  }
});

// A selected feature names its record by key. A request carries only the
// annotations of its own records, so drawing another record keeps the
// annotation in the draft until that record is drawn again (R2).
test('a Circular annotation of a selected feature waits in the draft while another record is drawn', async ({ page }) => {
  test.setTimeout(300_000);
  await openWithGenBank(page, BATCH_FIXTURE, () => {
    const app = window.__GBDRAW_APP__;
    app.form.multi_record_canvas = false;
    app.adv.circular_grouping_intent = 'single';
    app.form.circular_record_selector = '#1';
  });
  await generate(page);
  const cds = (await catalog(page)).find((feature) => feature.type === 'CDS');
  await annotateSelection(page, [cds.svgId]);
  await generate(page);
  const [first] = await drawnAnnotations(page);
  expect(first).toMatchObject({ id: 'feature_1', recordIndex: 0 });

  await page.evaluate(() => { window.__GBDRAW_APP__.form.circular_record_selector = '#2'; });
  await settle(page);
  await generate(page);
  expect(await drawnAnnotations(page)).toEqual([]);
  expect(await annotationWarnings(page)).toEqual([]);
  expect(await annotationTargets(page)).toHaveLength(1);

  await page.evaluate(() => { window.__GBDRAW_APP__.form.circular_record_selector = '#1'; });
  await settle(page);
  await generate(page);
  expect(await drawnAnnotations(page)).toEqual([first]);
});

// R2, OV-21 of the override-precedence audit: a Gallery Session's Linear
// record and a Circular grid of the same file both use the record key
// `record-1`, so a feature has the same [recordKey, feature ID] in both modes.
// An annotation of a feature selected in Circular belongs to the Circular
// drawing: the Linear request does not carry it, and it waits there for Circular.
test('OV-21: an annotation of a selected Circular feature stays out of Linear requests with the same record key', async ({ page }) => {
  test.setTimeout(360_000);
  const switchMode = async (mode) => {
    await page.getByRole('button', { name: mode === 'circular' ? 'Circular' : 'Linear', exact: true }).click();
    await page.waitForFunction((expected) => window.__GBDRAW_APP__?.mode === expected, mode);
    await settle(page);
  };
  const requestTargets = () => page.evaluate(async () => {
    const { getCommittedCanonicalRenderRequest } = await import('/gbdraw/web/js/services/config.js');
    const request = getCommittedCanonicalRenderRequest();
    return {
      mode: request.mode,
      recordKeys: request.records.map((record) => record.recordKey),
      targets: (request.diagramOptions.annotations?.sets || []).flatMap((set) => set.annotations.map((item) => item.target))
    };
  });
  await openFresh(page);
  await loadSessionFile(page, 'gbdraw/web/gallery/sessions/lambda_basic_linear.gbdraw-session.json');
  await switchMode('circular');
  await page.getByLabel('GenBank/DDBJ File', { exact: true }).setInputFiles('tests/test_inputs/NC_001416.gb');
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.circularRecordList.length)).toBeGreaterThan(0);
  // A one-record grid names its record `record-1`, as the Linear Session does.
  await page.evaluate(() => { window.__GBDRAW_APP__.form.multi_record_canvas = true; });
  await settle(page);
  await generate(page);
  const cds = (await catalog(page)).find((feature) => feature.type === 'CDS');
  expect(cds.recordKey).toBe('record-1');
  await annotateSelection(page, [cds.svgId]);
  await generate(page);
  const target = { kind: 'featureIdentity', recordKey: 'record-1', biologicalFeatureId: cds.identity,
    envelope: 'outer_bounds', circularPath: 'shortest' };
  expect(await requestTargets()).toEqual({ mode: 'circular', recordKeys: ['record-1'], targets: [target] });
  const circularDrawn = await drawnAnnotations(page);
  expect(circularDrawn).toMatchObject([{ id: 'feature_1', recordIndex: 0 }]);

  await switchMode('linear');
  await generate(page);
  expect(await requestTargets()).toEqual({ mode: 'linear', recordKeys: ['record-1'], targets: [] });
  expect(await drawnAnnotations(page)).toEqual([]);
  expect(await annotationWarnings(page)).toEqual([]);
  // A mode change and the Linear Generate keep the Circular target in the
  // Circular drawing (R2, PD-OI-086); the Linear drawing has none.
  expect(await annotationTargets(page)).toEqual([]);
  expect(await page.evaluate(async () => {
    const { state } = await import('/gbdraw/web/js/state.js');
    return JSON.parse(JSON.stringify(state.drawings.circular.annotationSets
      .flatMap((set) => set.annotations.map((item) => item.target))));
  })).toEqual([target]);

  await switchMode('circular');
  await generate(page);
  expect(await requestTargets()).toEqual({ mode: 'circular', recordKeys: ['record-1'], targets: [target] });
  expect(await drawnAnnotations(page)).toEqual(circularDrawn);
});

// R-7 (Owner decision 2026-10-05): a Session 44 named a selected feature by
// `hash=` (selected-feature-annotations.provenance.json). Load moves the target
// to the feature's source identity only when its record is drawn untransformed
// and the hash names one feature (feature_1); a hash of two CDS at the same
// coordinates (feature_2) and a target on a reverse-complemented record
// (feature_3) stay. The saved figure, the next Generate, and the CLI replay of
// the same Session draw the same annotations, and feature_1 stays on its
// feature after a later crop.
test('R-7: a Session 44 hash annotation of one untransformed feature moves to its identity without changing the figure', async ({ page }, testInfo) => {
  test.setTimeout(360_000);
  const SESSION = 'tests/fixtures/sessions/selected-feature-annotations.v44.gbdraw-session.json.gz';
  const dialogs = [];
  page.on('dialog', (dialog) => dialogs.push(dialog.message()));
  await openFresh(page);
  await page.locator('input[accept^=".json,"]').setInputFiles(SESSION);
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending
    && window.__GBDRAW_APP__.results.length > 0, null, { timeout: 180_000 });
  await settle(page);
  expect(dialogs).toEqual(['Session loaded successfully! 1 annotation(s) from an older Session named a feature by '
    + 'hash=. Each now names that feature by its source, so it stays on the feature when the crop or orientation changes.']);
  const testa = await page.evaluate(async () => {
    const { getCommittedCanonicalRenderRequest } = await import('/gbdraw/web/js/services/config.js');
    return getCommittedCanonicalRenderRequest().records[0].recordKey;
  });
  const targets = await annotationTargets(page);
  expect(targets[0]).toEqual({ kind: 'featureIdentity', recordKey: testa,
    biologicalFeatureId: 'fef810304', envelope: 'outer_bounds', circularPath: 'shortest' });
  expect(targets.slice(1).map((target) => [target.kind, target.selectors])).toEqual([
    ['featureSpan', [{ key: 'hash', value: 'f3ccacda4' }]],
    ['featureSpan', [{ key: 'hash', value: 'fbc76b20f' }]]
  ]);

  // Every drawn annotation shape, by annotation ID.
  const geometry = (svg) => page.evaluate((text) => {
    const root = new DOMParser().parseFromString(text, 'image/svg+xml').documentElement;
    return [...root.querySelectorAll('[data-gbdraw-annotation-id]')].map((group) => [
      group.getAttribute('data-gbdraw-annotation-id'),
      [group, ...group.querySelectorAll('*')].map((element) => [element.localName,
        ...['x', 'y', 'width', 'height', 'd', 'points', 'transform'].map((name) => element.getAttribute(name))])
    ]).sort(([a], [b]) => a.localeCompare(b));
  }, svg);
  const saved = await geometry(await page.evaluate(() => window.__GBDRAW_APP__.results[0].content));
  expect(saved.map(([id]) => id)).toEqual(['feature_1', 'feature_2']);
  await generate(page);
  expect(await geometry(await page.evaluate(() => window.__GBDRAW_APP__.results[0].content))).toEqual(saved);
  expect(await geometry(replayOnCli(SESSION, testInfo.outputPath('v44-replay')))).toEqual(saved);
  expect(await drawnAnnotations(page)).toEqual([
    { id: 'feature_2', recordIndex: 0, segments: [[300, 600]] },
    { id: 'feature_1', recordIndex: 0, segments: [[1600, 1901]] }
  ]);

  await page.evaluate(() => {
    const app = window.__GBDRAW_APP__;
    app.setLinearRecordCrop(app.linearSeqs[0], 'region_start', 701);
    app.setLinearRecordCrop(app.linearSeqs[0], 'region_end', 3800);
  });
  await settle(page);
  await generate(page);
  expect((await drawnAnnotations(page)).find((item) => item.id === 'feature_1'))
    .toEqual({ id: 'feature_1', recordIndex: 0, segments: [[900, 1201]] });
});
