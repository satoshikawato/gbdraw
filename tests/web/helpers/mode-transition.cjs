const { expect } = require('@playwright/test');
const fs = require('node:fs/promises');
const { openApp } = require('./app-lifecycle.cjs');

const seeds = {
  circular: 'gbdraw/web/gallery/sessions/HmmtDNA_basic_circular.gbdraw-session.json',
  linear: 'gbdraw/web/gallery/sessions/lambda_basic_linear.gbdraw-session.json'
};

const load = async (browser, file = seeds.circular, viewport = { width: 1600, height: 1000 }) => {
  const context = await browser.newContext({ viewport });
  const page = await context.newPage();
  page.externalRequests = [];
  await context.route('**/*', route => {
    if (new URL(route.request().url()).hostname === '127.0.0.1') return route.continue();
    page.externalRequests.push(route.request().url());
    return route.abort();
  });
  await page.addInitScript(() => {
    window.__MODE_EVENTS__ = [];
    window.__GBDRAW_TEST_HOOKS__ = {
      onSessionLifecycleEvent: event => window.__MODE_EVENTS__.push(event)
    };
  });
  page.on('dialog', dialog => dialog.type() === 'confirm' && dialog.message().startsWith('Download ')
    ? dialog.accept() : dialog.dismiss());
  await openApp(page);
  await page.locator('input[accept^=".json,"]').setInputFiles(file);
  await expect.poll(() => page.evaluate(() => !window.__GBDRAW_APP__.sessionImportPending
    && window.__GBDRAW_APP__.extractedFeatures.length > 0), { timeout: 180000 }).toBe(true);
  return page;
};

const generate = async page => {
  const key = await page.evaluate(async () => (await import('./js/state.js')).state.resultGenerationKey.value);
  await page.getByRole('button', { name: 'Generate Diagram', exact: true }).click();
  await expect.poll(() => page.evaluate(async () => {
    const { state } = await import('./js/state.js');
    return { key: state.resultGenerationKey.value, processing: state.processing.value, error: state.errorLog.value };
  }), { timeout: 180000 }).toEqual({ key: key + 1, processing: false, error: null });
};

const switchMode = async (page, mode) => {
  await page.getByRole('button', { name: mode === 'circular' ? 'Circular' : 'Linear', exact: true }).click();
  await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.mode)).toBe(mode);
  await expect.poll(() => page.evaluate(async () =>
    (await import('./js/state.js')).state.circularRecordDiscovery.status)).not.toBe('loading');
  // Settle the Vue render and its asynchronous preview binder, without a time delay.
  await page.evaluate(async () => { await window.Vue.nextTick(); });
};

const popup = async (page, index = 0) => {
  if (!await page.locator('.right-drawer').isVisible()) await page.locator('.drawer-toggle').click();
  await page.locator('.right-drawer').getByRole('button', { name: 'Edit', exact: true }).nth(index).click();
  await expect(page.locator('.feature-popup')).toBeVisible();
  return page.evaluate(() => {
    const feature = window.__GBDRAW_APP__.clickedFeature;
    return { id: feature.id, featureId: feature.featureId, sourceText: feature.labelSourceText };
  });
};

const closeEditor = async page => {
  if (await page.locator('.feature-popup').isVisible()) {
    await page.getByRole('button', { name: 'Close feature popup', exact: true }).click();
  }
  if (await page.locator('.right-drawer').isVisible()) await page.locator('.drawer-toggle').click();
};

const download = async (page, button, path) => {
  const pending = page.waitForEvent('download');
  await page.getByRole('button', { name: button, exact: true }).click();
  await (await pending).saveAs(path);
  return fs.readFile(path);
};

const snapshot = page => page.evaluate(async () => {
  const { state: s } = await import('./js/state.js');
  const ingestion = await import('./js/services/svg-result-ingestion.js');
  const { getCommittedCanonicalRenderRequest } = await import('./js/services/config.js');
  const root = s.svgContainer.value.querySelector('svg');
  const result = s.results.value[s.selectedResultIndex.value];
  return {
    mode: s.mode.value, generation: s.resultGenerationKey.value,
    labels: { ...s.labelTextFeatureOverrides }, bulkLabels: { ...s.labelTextBulkOverrides },
    labelSources: { ...s.labelTextFeatureOverrideSources }, visibility: { ...s.labelVisibilityOverrides },
    featureVisibility: { ...s.featureVisibilityOverrides }, visibilityRules: [...s.featureVisibilityManualRules],
    colors: { ...s.featureColorOverrides }, rules: s.manualSpecificRules.map(rule => ({ ...rule, fromFile: Boolean(rule.fromFile) })),
    context: s.labelOverrideContextKey.value, featureCount: s.extractedFeatures.value.length,
    resultIdentity: ingestion.getCommittedSvgResultRuntimeIdentity(result),
    markedMounted: ingestion.isCommittedSvgResultMounted(result),
    result: result.content, payload: s.svgContent.value, mounted: root.outerHTML,
    sameRoot: root === window.__MODE_EDITED_ROOT__,
    mountEvents: window.__MODE_EVENTS__.filter(e => e.name === 'preview.mount-observed'),
    request: getCommittedCanonicalRenderRequest()
  };
});

const semantics = (page, content) => page.evaluate(content => {
  const root = new DOMParser().parseFromString(content, 'image/svg+xml').documentElement;
  return [...root.querySelectorAll('[data-gbdraw-feature-id]')].map(element => ({
    id: element.getAttribute('data-gbdraw-feature-id'),
    part: element.getAttribute('data-gbdraw-feature-part'),
    d: element.getAttribute('d'), fill: element.getAttribute('fill'),
    stroke: element.getAttribute('stroke'), display: element.getAttribute('display')
  }));
}, content);

module.exports = { seeds, load, generate, switchMode, popup, closeEditor, download, snapshot, semantics };
