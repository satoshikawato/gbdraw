// Shared browser steps for the Web GUI audit 2026-09-30 regression specs.
const { expect } = require('@playwright/test');
const { evaluateWithRetainedPromise, generateAndWaitForResult, openApp } = require('./app-lifecycle.cjs');

const BATCH_FIXTURE = 'tests/fixtures/web_batch_two_records.gb';
const HMMT = 'tests/test_inputs/HmmtDNA.gbk';
const HMMT_SESSION = 'gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json';

const trackPreviewReadiness = (page) => page.addInitScript(() => {
  window.__AUDIT_LIFECYCLE__ = [];
  window.__GBDRAW_TEST_HOOKS__ = {
    ...(window.__GBDRAW_TEST_HOOKS__ || {}),
    onSessionLifecycleEvent: (event) => window.__AUDIT_LIFECYCLE__.push(event.name)
  };
});

const readyCount = (page) => page.evaluate(() => (
  window.__AUDIT_LIFECYCLE__.filter((name) => name === 'preview.ready-receipt-accepted').length
));

const settle = async (page) => {
  await page.waitForFunction(() => {
    const app = window.__GBDRAW_APP__;
    const history = window.__GBDRAW_HISTORY__;
    return !app.processing && !app.labelReflowProcessing && !app.sessionImportPending && !app.ruleMatchingPending
      && !history?.capturing?.value && !history?.restoring?.value;
  });
  await page.evaluate(() => new Promise((resolve) => (
    requestAnimationFrame(() => requestAnimationFrame(resolve))
  )));
};

const openFresh = async (page) => {
  page.on('dialog', (dialog) => dialog.accept());
  await trackPreviewReadiness(page);
  await openApp(page);
};

const uploadCircularGenBank = async (page, path, configure = null) => {
  await page.getByLabel('GenBank/DDBJ File', { exact: true }).setInputFiles(path);
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.circularRecordList.length), {
    timeout: 60_000
  }).toBeGreaterThan(0);
  if (configure) await page.evaluate(configure);
  await settle(page);
};

const openWithGenBank = async (page, path, configure = null) => {
  await openFresh(page);
  await uploadCircularGenBank(page, path, configure);
};

const loadSessionFile = async (page, path) => {
  await page.locator('input[accept^=".json,"]').setInputFiles(path);
  await page.waitForFunction(() => !window.__GBDRAW_APP__.sessionImportPending
    && window.__GBDRAW_APP__.extractedFeatures.length > 0, null, { timeout: 180_000 });
  await settle(page);
};

const openBatch = async (page) => {
  await openWithGenBank(page, BATCH_FIXTURE, () => {
    const app = window.__GBDRAW_APP__;
    app.form.labels_mode = 'out';
    app.form.multi_record_canvas = false;
    app.adv.circular_grouping_intent = 'batch';
    app.autoLabelReflowEnabled = false;
  });
  await generateAndWaitForResult(page);
  expect(await page.evaluate(() => window.__GBDRAW_APP__.results.length)).toBe(2);
  await settle(page);
};

const selectResult = async (page, index) => {
  const before = await readyCount(page);
  await page.locator('h2 select').selectOption({ index });
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.selectedResultIndex)).toBe(index);
  await expect.poll(() => readyCount(page)).toBeGreaterThan(before);
  await settle(page);
};

const switchMode = async (page, mode) => {
  await page.getByRole('button', { name: mode === 'circular' ? 'Circular' : 'Linear', exact: true }).click();
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.mode)).toBe(mode);
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.circularRecordDiscovery?.status))
    .not.toBe('loading');
  await settle(page);
};

// Feature presentation in the mounted SVG and in each committed Result.
const featurePresentation = (page, locusTags) => page.evaluate((tags) => {
  const app = window.__GBDRAW_APP__;
  const mounted = app.svgContainer.querySelector('svg');
  const documents = app.results.map((result) => new DOMParser().parseFromString(result.content || '', 'image/svg+xml'));
  const describe = (root, id) => {
    const elements = [...root.querySelectorAll(`[data-gbdraw-feature-id="${CSS.escape(id)}"]`)];
    if (!elements.length) return null;
    const fill = elements.map((element) => element.getAttribute('fill')).find((value) => value && value !== 'none');
    return `${fill}|${elements.some((element) => element.getAttribute('display') === 'none') ? 'hidden' : 'shown'}`;
  };
  return Object.fromEntries(tags.map((tag) => {
    const feature = app.extractedFeatures.find((item) => item.locus_tag === tag);
    return [tag, {
      mounted: feature ? describe(mounted, feature.svg_id) : null,
      results: feature ? documents.map((document) => describe(document, feature.svg_id)) : []
    }];
  }));
}, locusTags);

const legendCaptions = (page, { source = 'mounted' } = {}) => page.evaluate(async (from) => {
  const { getVisibleFeatureLegendGroup } = await import('/gbdraw/web/js/services/legend-svg.js');
  const app = window.__GBDRAW_APP__;
  const root = from === 'mounted'
    ? app.svgContainer.querySelector('svg')
    : new DOMParser().parseFromString(app.results[app.selectedResultIndex].content, 'image/svg+xml').documentElement;
  const group = getVisibleFeatureLegendGroup(root);
  const entries = [...(group?.querySelectorAll('g[data-legend-key]') || [])].map((entry) => {
    const text = entry.querySelector('text');
    const match = String(text?.getAttribute('transform') || '').match(/translate\(\s*([-\d.e]+)[ ,]+([-\d.e]+)/);
    return { caption: text?.textContent || entry.getAttribute('data-legend-key'), x: match ? Number(match[1]) : 0, y: match ? Number(match[2]) : 0 };
  });
  return entries
    .sort((left, right) => (Math.abs(left.y - right.y) < 1 ? left.x - right.x : left.y - right.y))
    .map(({ caption }) => caption);
}, source);

// One feature edit through the editor actions the drawer controls call. A
// color or visibility edit with a scope dialog takes `scope`/`visibilityScope`
// (default: This feature only / the first scope, This feature). The edit
// awaits Python's rule matches and the projection, so it runs on the retained
// path: Chromium can collect the promise a long page.evaluate awaits.
const editFeature = (page, locusTag, edit) => evaluateWithRetainedPromise(page, async ({ tag, change }) => {
  const app = window.__GBDRAW_APP__;
  const feature = app.extractedFeatures.find((item) => item.locus_tag === tag);
  if (!feature) throw new Error(`no feature ${tag}`);
  await app.openFeatureEditorFromList(feature, null);
  if (change.fill) {
    await app.updateClickedFeatureColor(change.fill);
    if (app.featureStyleScopeDialog.show) await app.handleFeatureStyleScopeChoice(change.scope || 'single');
  }
  if (change.labelText !== undefined || change.labelVisibility) {
    if (change.labelText !== undefined) app.clickedFeature.labelText = change.labelText;
    if (change.labelVisibility) app.clickedFeature.labelVisibility = change.labelVisibility;
    // Label Not Shown keeps the Apply open until its choice (Keep hidden here).
    let applying = true;
    const applied = Promise.resolve(app.updateClickedFeatureLabelText()).finally(() => { applying = false; });
    while (applying && !app.hiddenLabelTextDialog?.show) await new Promise((resolve) => setTimeout(resolve, 20));
    if (app.hiddenLabelTextDialog?.show) await app.handleHiddenLabelTextChoice('text_only');
    await applied;
  }
  if (change.visibility) {
    app.clickedFeature.featureVisibility = change.visibility;
    await app.updateClickedFeatureVisibility(change.visibility);
    if (app.featureVisibilityScopeDialog.show) {
      await app.handleFeatureVisibilityScopeChoice(change.visibilityScope || app.featureVisibilityScopeDialog.scopes[0].id);
    }
  }
  app.clickedFeature = null;
}, { tag: locusTag, change: edit });

module.exports = {
  BATCH_FIXTURE,
  HMMT,
  HMMT_SESSION,
  editFeature,
  featurePresentation,
  legendCaptions,
  loadSessionFile,
  openBatch,
  openFresh,
  openWithGenBank,
  readyCount,
  selectResult,
  settle,
  switchMode,
  trackPreviewReadiness,
  uploadCircularGenBank
};
