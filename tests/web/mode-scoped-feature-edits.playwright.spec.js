// R2, residual R-1 of the override-precedence audit: a Gallery Session's
// Linear record and a Circular grid of the same file both use the record key
// `record-1`, so the same feature has the same [recordKey, feature ID] in both
// modes. A Feature placement (Main) and a Feature visibility edit made in
// Circular apply only to Circular requests; the Linear Generate draws the
// feature with neither, and the edits wait in the draft for Circular.
const { test, expect } = require('@playwright/test');
const { generateAndWaitForResult } = require('./helpers/app-lifecycle.cjs');
const { openFresh, loadSessionFile, settle } = require('./helpers/audit-browser.cjs');

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
    drafts: [Object.keys(state.featurePlacementOverrides), Object.keys(state.featureOverrides)]
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

  const target = await page.evaluate(async () => {
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
  // The Linear Result draws the feature; the popup reads no edit of its own.
  expect(linear.drawn).toBeGreaterThan(0);
  expect(await page.evaluate((id) => {
    const app = window.__GBDRAW_APP__;
    const feature = app.extractedFeatures.find((item) => item.biological_feature_id === id);
    return [feature.scope, app.featurePlacementActions.valueFor(feature), app.getFeatureVisibility(feature)];
  }, target.id)).toEqual(['linear', 'auto', 'default']);
  // A mode change and the Linear Generate keep the Circular rows (R2).
  expect(linear.drafts).toEqual([[key], [key]]);

  await switchMode(page, 'circular');
  await generateAndWaitForResult(page);
  const circular = await committed(page, target.id);
  expect(circular).toMatchObject({ mode: 'circular', recordKeys: ['record-1'], ...circularRows, drawn: 0 });
});
