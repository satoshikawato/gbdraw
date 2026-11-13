const { test, expect } = require('@playwright/test');
const { readFileSync, writeFileSync } = require('node:fs');
const { resolve } = require('node:path');
const { gunzipSync } = require('node:zlib');
const { openApp, getDiagramWorkerActivity, evaluateWithRetainedPromise, generateAndWaitForResult } = require('./helpers/app-lifecycle.cjs');
const { expectLiveEqualsGenerate, settleLive } = require('./helpers/live-generate-parity.cjs');
const { openBatch } = require('./helpers/audit-browser.cjs');
const { appAction, addColorRule, history } = require('./helpers/live-generate-parity-steps.cjs');

const fixture = resolve('gbdraw/web/gallery/sessions/tobacco-chloroplast.gbdraw');

const visibleRuleColor = (page, color = '#ff0000') => page.evaluate(color => {
  const svg = (window.__GBDRAW_APP__.currentResults()[0]?.content || '').toLowerCase();
  const mounted = [...document.querySelectorAll('svg [fill], svg [style]')]
    .filter(node => `${node.getAttribute('fill') || ''} ${node.getAttribute('style') || ''}`.toLowerCase().includes(color));
  return { result: svg.split(color).length - 1, mounted: mounted.length };
}, color);

// The rule the color and predicate edits change: the tRNA-by-gene rule. Its
// features fall to the tRNA-by-product rule, which draws the same Legend row
// in the same color, so a predicate restore returns no Legend row that
// Python has to draw again (a rerender would prepare the rules once more).
const editedRule = (rules) => rules.findIndex((rule) => rule.feat === 'tRNA' && rule.qual === 'gene');

const edit = (page, label, action) => page.evaluate(async ({ label, action, editedRuleSource }) => {
  const { state } = await import('./js/state.js');
  const editedRule = new Function(`return (${editedRuleSource})`)();
  await window.__GBDRAW_HISTORY__.runUndoable(label, () => {
    const rules = state.activeDrawing().manualSpecificRules;
    if (action === 'scale') state.activeDrawing().adv.scale_interval = 12345;
    if (action === 'scaleAgain') state.activeDrawing().adv.scale_interval = 23456;
    if (action === 'color') rules[editedRule(rules)].color = '#ff0000';
    if (action === 'predicate') rules[editedRule(rules)].val = 'trn_NO_MATCH';
  });
}, { label, action, editedRuleSource: editedRule.toString() });

const restore = (page, direction) => page.evaluate(direction => window.__GBDRAW_HISTORY__[direction](), direction);

test('History prepares cold rules once, skips warm config preparation, and rechecks changed predicates', async ({ page }) => {
  test.setTimeout(180_000);
  page.on('dialog', dialog => dialog.accept());
  await openApp(page);
  await page.locator('input[accept*=".gbdraw"]').setInputFiles(fixture);
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending);
  const context = page.context();
  const cdp = await context.newCDPSession(page);
  await cdp.send('Profiler.enable');
  await cdp.send('Profiler.startPreciseCoverage', { callCount: true, detailed: true });
  const prepareCalls = async () => (await cdp.send('Profiler.takePreciseCoverage')).result
    .filter(entry => entry.url.endsWith('/app/rule-matching.js'))
    .flatMap(entry => entry.functions)
    .filter(entry => entry.functionName === 'prepare')
    .reduce((total, entry) => total + Math.max(...entry.ranges.map(range => range.count)), 0);
  try {
    const before = await getDiagramWorkerActivity(page);
    expect(await visibleRuleColor(page)).toEqual({ result: 0, mounted: 0 });
    const baseGreen = await visibleRuleColor(page, '#2c6d3e');
    expect(baseGreen.result).toBeGreaterThan(0);
    expect(baseGreen.mounted).toBeGreaterThan(0);
    await edit(page, 'Unrelated scale', 'scale');
    expect(await prepareCalls()).toBe(0);
    await restore(page, 'undo');
    expect(await prepareCalls()).toBe(1);
    expect(await visibleRuleColor(page)).toEqual({ result: 0, mounted: 0 });
    expect(await visibleRuleColor(page, '#2c6d3e')).toEqual(baseGreen);
    await restore(page, 'redo');
    expect(await prepareCalls()).toBe(0);
    expect(await visibleRuleColor(page, '#2c6d3e')).toEqual(baseGreen);
    const warm = await getDiagramWorkerActivity(page);
    // The Session has only managed `configOverrides`, so `checkUnmanagedConfigOverrides`
    // (`config.js`) returns before starting the Worker: the Load builds none.
    expect(before.constructions).toBe(0);
    expect(warm.constructions).toBe(1);
    expect(warm.helpers - before.helpers).toBe(1);
    expect(warm.runs).toBe(before.runs);
    await edit(page, 'Rule color', 'color');
    expect(await prepareCalls()).toBe(0);
    await restore(page, 'undo');
    expect(await prepareCalls()).toBe(1);
    expect(await visibleRuleColor(page)).toEqual({ result: 0, mounted: 0 });
    expect(await visibleRuleColor(page, '#2c6d3e')).toEqual(baseGreen);
    await restore(page, 'redo');
    expect(await prepareCalls()).toBe(1);
    const colored = await visibleRuleColor(page);
    expect(colored.result).toBeGreaterThan(0);
    expect(colored.mounted).toBeGreaterThan(0);
    await edit(page, 'Rule predicate', 'predicate');
    expect(await prepareCalls()).toBe(0);
    await restore(page, 'undo');
    expect(await prepareCalls()).toBe(1);
    expect(await visibleRuleColor(page)).toEqual(colored);
    await restore(page, 'redo');
    expect(await prepareCalls()).toBe(1);
    expect(await visibleRuleColor(page)).toEqual({ result: 0, mounted: 0 });
    await edit(page, 'Independent scale', 'scaleAgain');
    expect(await prepareCalls()).toBe(0);
    expect(await page.evaluate((editedRuleSource) => {
      const rules = window.__GBDRAW_APP__.manualSpecificRules;
      return rules[new Function(`return (${editedRuleSource})`)()(rules)].val;
    }, editedRule.toString())).toBe('trn_NO_MATCH');
    expect(await visibleRuleColor(page)).toEqual({ result: 0, mounted: 0 });
    await restore(page, 'undo');
    expect(await prepareCalls()).toBe(0);
    await restore(page, 'redo');
    expect(await prepareCalls()).toBe(0);
    expect(await visibleRuleColor(page)).toEqual({ result: 0, mounted: 0 });
    const after = await getDiagramWorkerActivity(page);
    expect(after.constructions).toBe(warm.constructions);
    expect(after.helpers - warm.helpers).toBe(1);
    expect(after.runs).toBe(before.runs);
  } finally {
    await cdp.send('Profiler.stopPreciseCoverage');
  }
});

// D-15-6 (3): Undo and Redo of a rule color and a rule predicate change, made
// through the rule owner, show the Legend rows and fills Generate draws for
// the restored rules. This case guards the rule-owner path of the restore; the
// direct table edit of OV-347 is pinned by the case above and the
// palette-default-colors node case, and its live caption differs from Generate's
// (OV-357, accepted in D-15-38).
test('Redo of a rule predicate change shows what Generate draws', async ({ page }) => {
  test.setTimeout(240_000);
  page.on('dialog', dialog => dialog.accept());
  await openApp(page);
  await page.locator('input[accept*=".gbdraw"]').setInputFiles(fixture);
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending);
  await appAction(page, 'setSpecificRuleField', 0, 'color', '#ff0000');
  await appAction(page, 'setSpecificRuleField', 0, 'val', 'psaA_NO_MATCH');
  await history(page, 'undo');
  await history(page, 'undo');
  await history(page, 'redo');
  await history(page, 'redo');
  await expectLiveEqualsGenerate(page, { label: 'Redo of a rule predicate change' });
});

// OV-3285: a feature color edit and a product-scope hide, both undone, end
// in the restore of a History checkpoint. A rule added after it shows live
// what Generate draws: the restore installs the catalog's own features, on
// which the rule matches are recorded, not copies of them.
const undoneFeatureEdits = async (page) => {
  await openBatch(page);
  await evaluateWithRetainedPromise(page, async () => {
    const app = window.__GBDRAW_APP__;
    const colored = app.extractedFeatures.find((item) => item.locus_tag === 'TESTA_0001');
    await app.openFeatureEditorFromList(colored, null);
    await app.requestFeatureColorChange(colored, '#ff0000');
    await app.handleFeatureStyleScopeChoice('annotationLabel');
    const hidden = app.extractedFeatures.find((item) => item.locus_tag === 'TESTA_0004');
    await app.openFeatureEditorFromList(hidden, null);
    await app.updateClickedFeatureVisibility('off');
    await app.handleFeatureVisibilityScopeChoice('product');
    app.clickedFeature = null;
  });
  await settleLive(page);
  await history(page, 'undo');
  await history(page, 'undo');
};

test('A visibility rule added after the Undo of feature edits shows what Generate draws', async ({ page }) => {
  test.setTimeout(240_000);
  await undoneFeatureEdits(page);
  await evaluateWithRetainedPromise(page, async () => {
    const app = window.__GBDRAW_APP__;
    await app.addFeatureVisibilityRule();
    const index = app.featureVisibilityManualRules.length - 1;
    const rule = { recordId: '*', featureType: 'CDS', qualifier: 'locus_tag', value: '_0004$', action: 'off' };
    for (const [field, value] of Object.entries(rule)) await app.setFeatureVisibilityRuleField(index, field, value);
  });
  await expectLiveEqualsGenerate(page, { label: 'Visibility rule after the Undo of feature edits' });
});

test('A color rule added after the Undo of feature edits shows what Generate draws', async ({ page }) => {
  test.setTimeout(240_000);
  await undoneFeatureEdits(page);
  await addColorRule(page, { feat: 'CDS', qual: 'locus_tag', val: '_0004$', color: '#ff00ff', cap: 'fourth' });
  await expectLiveEqualsGenerate(page, { label: 'Color rule after the Undo of feature edits' });
});

// What a Generate that promotes legacy protein cache entries installs: the
// catalog's own feature objects, its orthogroups and the selected alignment feature.
const migratedFeatureState = (page) => page.evaluate(async () => {
  const { state } = await import('./js/state.js');
  const { featureStateFromCatalog } = await import('./js/services/feature-catalog.js');
  const raw = (value) => window.Vue.toRaw(value);
  const catalog = featureStateFromCatalog(raw(state.featureCatalog.value), { mode: 'linear' });
  return {
    promoted: state.legacyProteinRawCandidates.value.entries.filter((entry) => entry.state === 'promoted').length > 0,
    catalogFeaturesInstalled: raw(state.extractedFeatures.value).every((feature, index) => feature === catalog.extractedFeatures[index])
      && (raw(state.biologicalFeatures.value) || []).every((feature, index) => feature === catalog.biologicalFeatures[index]),
    orthogroups: JSON.stringify(state.orthogroups.value),
    selectedAlignmentFeature: state.selectedOrthogroupAlignmentFeature.value,
    // The Generate draws from current protein handles only: no legacy `p_r_` protein ID is left to rewrite.
    legacyProteinIds: JSON.stringify([
      state.orthogroups.value, state.selectedOrthogroupAlignmentFeature.value, state.extractedFeatures.value
    ]).match(/p_r_[0-9a-f]+_\d+_\d+_/g)?.length ?? 0
  };
});

// OV-3352: a Generate that promotes legacy protein cache entries installs the
// catalog's own features, as every other Generate does, so a rule added after
// it shows live what Generate draws. The v33 Session's Web draft says Pairwise
// while its stored result is an Orthogroup run, as when a user changed the
// mode after Generate and saved: the Load promotes only the cache entries of
// the adjacent pairs, and the Orthogroup Generate promotes the others.
test('A visibility rule added after a Generate that migrates legacy protein caches shows what Generate draws', async ({ page }, testInfo) => {
  test.setTimeout(300_000);
  const session = JSON.parse(gunzipSync(readFileSync(
    resolve('tests/fixtures/sessions/BGC0000708-BGC0000713.schema-v2.gbdraw-session.json.gz')
  )).toString('utf8'));
  session.config.losat.blastp.mode = 'pairwise';
  const pairwiseDraft = testInfo.outputPath('BGC0000708-BGC0000713.v33-pairwise-draft.gbdraw-session.json');
  writeFileSync(pairwiseDraft, JSON.stringify(session));
  page.on('dialog', dialog => dialog.accept());
  await openApp(page);
  await page.locator('input[accept*=".gbdraw"]').setInputFiles(pairwiseDraft);
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending, null, { timeout: 180_000 });
  await page.evaluate(() => { window.__GBDRAW_APP__.losat.blastp.mode = 'orthogroup'; });
  await generateAndWaitForResult(page, { timeout: 300_000 });
  const generated = await migratedFeatureState(page);
  expect(generated).toMatchObject({ promoted: true, catalogFeaturesInstalled: true, legacyProteinIds: 0 });
  // The History step of this Generate: the newest step before the rule edits (`timeline().undo` lists the steps
  // newest first).
  const generateStep = await page.evaluate(() => window.__GBDRAW_HISTORY__.timeline().undo.at(0).sequence);
  await evaluateWithRetainedPromise(page, async () => {
    const app = window.__GBDRAW_APP__;
    await app.addFeatureVisibilityRule();
    const index = app.featureVisibilityManualRules.length - 1;
    const rule = { recordId: '*', featureType: 'CDS', qualifier: 'protein_id', value: '^CAG38690\\.1$', action: 'off' };
    for (const [field, value] of Object.entries(rule)) await app.setFeatureVisibilityRuleField(index, field, value);
  });
  await expectLiveEqualsGenerate(page, { label: 'Visibility rule after a migrating Generate' });
  await test.step('OV-3352: going back in History to the migrating Generate keeps its features', async () => {
    await page.evaluate((sequence) => window.__GBDRAW_HISTORY__.travelTo(sequence), generateStep);
    await settleLive(page);
    // The restored step holds the catalog's own feature objects, not copies of them.
    expect(await migratedFeatureState(page)).toEqual({ ...generated, catalogFeaturesInstalled: true });
  });
});
