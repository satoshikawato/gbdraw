// R2, residual R-1 of the override-precedence audit: a Gallery Session's
// Linear record and a Circular grid of the same file both use the record key
// `record-1`, so the same feature has the same [recordKey, feature ID] in both
// modes. A Feature placement (Main) and a Feature visibility edit made in
// Circular apply only to Circular requests; the Linear Generate draws the
// feature with neither, and the edits wait in the draft for Circular.
const { test, expect } = require('@playwright/test');
const { readFileSync } = require('node:fs');
const { gunzipSync } = require('node:zlib');
const { basename } = require('node:path');
const { generateAndWaitForResult, evaluateWithRetainedPromise } = require('./helpers/app-lifecycle.cjs');
const { openFresh, loadSessionFile, settle } = require('./helpers/audit-browser.cjs');
const { download } = require('./helpers/mode-transition.cjs');

const LINEAR_SESSION = 'gbdraw/web/gallery/sessions/lambda_basic_linear.gbdraw-session.json';
const LAMBDA = 'tests/test_inputs/NC_001416.gb';

const committed = (page, id) => page.evaluate(async (featureId) => {
  const { getCommittedCanonicalRenderRequest } = await import('./js/services/config.js');
  const { state } = await import('./js/state.js');
  const request = getCommittedCanonicalRenderRequest();
  return {
    mode: request.mode,
    recordKeys: request.records.map((record) => record.recordKey),
    featurePlacements: request.diagramOptions.featurePlacements,
    featureOverrides: request.diagramOptions.featureOverrides,
    drawn: state.extractedFeatures.value.filter((feature) => feature.biological_feature_id === featureId).length,
    drafts: [Object.keys(state.activeDrawing().featurePlacementOverrides), Object.keys(state.activeDrawing().featureOverrides)]
  };
}, id);

const switchMode = async (page, mode) => {
  await page.getByRole('button', { name: mode === 'circular' ? 'Circular' : 'Linear', exact: true }).click();
  await page.waitForFunction((expected) => window.__GBDRAW_APP__?.mode === expected, mode);
  await settle(page);
};

test('a Circular Main placement and Feature visibility edit stay out of Linear requests with the same record key', async ({ page }) => {
  test.setTimeout(360_000);
  await openFresh(page);
  await loadSessionFile(page, LINEAR_SESSION);
  await switchMode(page, 'circular');
  await page.getByLabel('GenBank/DDBJ File', { exact: true }).setInputFiles(LAMBDA);
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.circularRecordList.length)).toBeGreaterThan(0);
  // A one-record grid names its record `record-1`, as the Linear Session does.
  await page.evaluate(() => { window.__GBDRAW_APP__.form.multi_record_canvas = true; });
  await settle(page);
  await generateAndWaitForResult(page);

  const target = await evaluateWithRetainedPromise(page, async () => {
    const app = window.__GBDRAW_APP__;
    const feature = app.extractedFeatures.find((item) => item.type === 'CDS');
    await app.featurePlacementActions.setPlacement([feature], 'main');
    await app.openFeatureEditorFromList(feature, null);
    app.clickedFeature.featureVisibility = 'off';
    await app.updateClickedFeatureVisibility('off');
    if (app.featureVisibilityScopeDialog.show) await app.handleFeatureVisibilityScopeChoice('feature');
    app.clickedFeature = null;
    return { recordKey: feature.record_key, id: feature.biological_feature_id };
  });
  await settle(page);
  expect(target.recordKey).toBe('record-1');
  const key = JSON.stringify(['circular', target.recordKey, target.id]);
  const circularRows = {
    featurePlacements: [{ recordKey: target.recordKey, biologicalFeatureId: target.id, placement: { kind: 'main' } }],
    featureOverrides: [{ recordKey: target.recordKey, biologicalFeatureId: target.id,
      featureVisibility: 'off', labelVisibility: null, labelText: null }]
  };

  await switchMode(page, 'linear');
  await generateAndWaitForResult(page);
  const linear = await committed(page, target.id);
  expect(linear).toMatchObject({ mode: 'linear', recordKeys: ['record-1'], featurePlacements: [], featureOverrides: [] });
  // The Linear Result draws the feature; the feature reads no edit of its own,
  // through the owners the popup and the Features list read.
  expect(linear.drawn).toBeGreaterThan(0);
  expect(await page.evaluate(async (id) => {
    const app = window.__GBDRAW_APP__;
    const { state } = await import('./js/state.js');
    const { getFeatureVisibilityOverride } = await import('./js/services/feature-visibility.js');
    const feature = app.extractedFeatures.find((item) => item.biological_feature_id === id);
    return [
      feature.scope,
      app.featurePlacementActions.valueFor(feature),
      getFeatureVisibilityOverride(state.activeDrawing().featureOverrides, feature),
      app.featureListState(feature).drawn
    ];
  }, target.id)).toEqual(['linear', 'auto', 'default', true]);
  // A mode change and the Linear Generate keep the Circular rows (R2).
  expect(linear.drafts).toEqual([[key], [key]]);

  await switchMode(page, 'circular');
  await generateAndWaitForResult(page);
  const circular = await committed(page, target.id);
  expect(circular).toMatchObject({ mode: 'circular', recordKeys: ['record-1'], ...circularRows, drawn: 0 });
});

// OV-84 (R2, OIPC-C06): a per-feature stroke is keyed `recordKey\0featureId`
// in the drawing of its mode (PR-1). A Generate of the other mode neither draws
// nor removes it; a saved Session keeps each mode's strokes in its slice, and
// the next Generate of the stroke's mode draws it again.
const STROKE_FIXTURE = 'tests/fixtures/forced_label_underlay.gb';
const STROKE = { strokeColor: '#ff0000', strokeWidth: 3 };

const loadStrokeInput = async (page, mode, path = STROKE_FIXTURE) => {
  if (mode === 'circular') {
    await page.getByLabel('GenBank/DDBJ File', { exact: true }).setInputFiles(path);
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.circularRecordList.length)).toBeGreaterThan(0);
  } else {
    await page.evaluate(async ({ text, name }) => {
      window.__GBDRAW_APP__.setLinearSeqPrimaryFile(0, 'gb', new File([text], name, { type: 'text/plain', lastModified: 1000 }));
      await window.Vue.nextTick();
    }, { text: readFileSync(path, 'utf8'), name: basename(path) });
  }
  await settle(page);
};

// Sets the stroke of the first CDS from its popup ("This feature only").
const strokeFirstCds = (page) => evaluateWithRetainedPromise(page, async ({ strokeColor, strokeWidth }) => {
  const app = window.__GBDRAW_APP__;
  const feature = app.extractedFeatures.find((item) => item.type === 'CDS');
  await app.openFeatureEditorFromList(feature, null);
  await window.Vue.nextTick();
  if (await app.updateClickedFeatureStroke(strokeColor, strokeWidth) !== true) throw new Error('stroke not applied');
  app.clickedFeature = null;
  return feature.stable_override_key;
}, STROKE);

// The draft's stroke on `key`, and whether the displayed Result draws it.
const strokeOf = (page, key) => page.evaluate(async ({ strokeKey, strokeColor, strokeWidth }) => {
  const { getFeatureElements } = await import('./js/services/feature-dom.js');
  const app = window.__GBDRAW_APP__;
  const feature = app.extractedFeatures.find((item) => item.stable_override_key === strokeKey);
  const content = String(app.results[app.selectedResultIndex]?.content || '');
  const svg = new DOMParser().parseFromString(content, 'image/svg+xml').documentElement;
  const elements = feature ? getFeatureElements(svg, feature.svg_id) : [];
  const override = app.featureStrokeOverrides[strokeKey];
  return {
    override: override ? { strokeColor: override.strokeColor, strokeWidth: override.strokeWidth } : null,
    drawn: elements.length > 0 && elements.every((element) => element.getAttribute('stroke') === strokeColor
      && element.getAttribute('stroke-width') === String(strokeWidth))
  };
}, { strokeKey: key, ...STROKE });

// Saves the Session and loads it into a fresh page; returns the loaded strokes.
const strokesAfterSessionRoundTrip = async (page, browser, testInfo, keys, name) => {
  const saved = testInfo.outputPath(`${name}.gbdraw-session.json`);
  await download(page, 'Save Session', saved);
  const context = await browser.newContext({ baseURL: new URL(page.url()).origin });
  try {
    const fresh = await context.newPage();
    await openFresh(fresh);
    await fresh.locator('input[accept^=".json,"]').setInputFiles(saved);
    await fresh.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending
      && window.__GBDRAW_APP__.results.length > 0, null, { timeout: 180_000 });
    await settle(fresh);
    return Object.fromEntries(await Promise.all(keys.map(async (key) => [key, (await strokeOf(fresh, key)).override])));
  } finally {
    await context.close();
  }
};

for (const [mode, otherMode] of [['circular', 'linear'], ['linear', 'circular']]) {
  test(`a ${mode} feature stroke survives a ${otherMode} Generate and a Session saved after it (OV-84)`,
    async ({ page, browser }, testInfo) => {
      test.setTimeout(360_000);
      await openFresh(page);
      if (mode === 'linear') await switchMode(page, 'linear');
      await loadStrokeInput(page, mode);
      await generateAndWaitForResult(page);
      const key = await strokeFirstCds(page);
      await settle(page);
      expect(await strokeOf(page, key)).toEqual({ override: STROKE, drawn: true });

      await switchMode(page, otherMode);
      await loadStrokeInput(page, otherMode);
      await generateAndWaitForResult(page);
      // The other mode draws the same file under its own record key, and its
      // drawing holds no stroke of the first mode.
      expect(await page.evaluate((strokeKey) => window.__GBDRAW_APP__.extractedFeatures
        .some((item) => item.stable_override_key === strokeKey), key)).toBe(false);
      expect(await strokeOf(page, key)).toEqual({ override: null, drawn: false });
      // A stroke of each mode: the Session keeps each in its mode's slice.
      const otherKey = await strokeFirstCds(page);
      await settle(page);
      expect(otherKey).not.toBe(key);
      const bytes = await download(page, 'Save Session', testInfo.outputPath(`${mode}-stroke.gbdraw-session.json`));
      const saved = JSON.parse((bytes[0] === 0x1f ? gunzipSync(bytes) : bytes).toString('utf8'));
      expect(saved.modes[mode].editorState.featureStrokes.overrides[key]).toMatchObject(STROKE);
      expect(saved.modes[otherMode].editorState.featureStrokes.overrides[otherKey]).toMatchObject(STROKE);
      expect(saved.modes[otherMode].editorState.featureStrokes.overrides[key]).toBeUndefined();
      expect(await strokesAfterSessionRoundTrip(page, browser, testInfo, [otherKey], `${mode}-stroke-loaded`))
        .toEqual({ [otherKey]: STROKE });

      await switchMode(page, mode);
      await generateAndWaitForResult(page);
      expect(await strokeOf(page, key)).toEqual({ override: STROKE, drawn: true });
      expect(await strokeOf(page, otherKey)).toEqual({ override: null, drawn: false });
    });
}

// OV-84 (Owner-delegated 2026-10-07): a source replacement retires the stroke
// of a feature the new file does not have, as it does the feature's label
// edits; Generate and Save Session do not fail, and the stroke does not return
// with the old file.
test('a feature stroke on a replaced Circular file is retired without failing Generate or Save Session (OV-84)',
  async ({ page, browser }, testInfo) => {
    test.setTimeout(360_000);
    await openFresh(page);
    await loadStrokeInput(page, 'circular');
    await generateAndWaitForResult(page);
    const key = await strokeFirstCds(page);
    await settle(page);

    await loadStrokeInput(page, 'circular', LAMBDA);
    await generateAndWaitForResult(page);
    expect(await strokeOf(page, key)).toEqual({ override: null, drawn: false });
    expect(await page.evaluate(() => window.__GBDRAW_APP__.results[0].content.includes('stroke="#ff0000"'))).toBe(false);
    expect(await strokesAfterSessionRoundTrip(page, browser, testInfo, [key], 'replaced-stroke'))
      .toEqual({ [key]: null });

    await loadStrokeInput(page, 'circular');
    await generateAndWaitForResult(page);
    expect(await strokeOf(page, key)).toEqual({ override: null, drawn: false });
  });
