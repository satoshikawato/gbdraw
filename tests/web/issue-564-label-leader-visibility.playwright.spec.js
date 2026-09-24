const { test, expect } = require('@playwright/test');
const { openApp } = require('./helpers/app-lifecycle.cjs');
const {
  seeds,
  load,
  generate,
  download
} = require('./helpers/mode-transition.cjs');

test.describe.configure({ retries: 0 });

const inspectUnit = (page, featureId) => page.evaluate((targetId) => {
  const collect = (root) => Array.from(root.querySelectorAll('[data-label-feature-id]'))
    .filter((element) => element.getAttribute('data-label-feature-id') === targetId)
    .map((element) => ({
      tag: element.localName,
      display: element.getAttribute('display'),
      preview: element.getAttribute('data-gbdraw-label-visibility-preview'),
      schema: element.getAttribute('data-gbdraw-label-binding-schema')
    }));
  const app = window.__GBDRAW_APP__;
  const mounted = document.querySelector('.origin-top svg');
  const result = new DOMParser().parseFromString(
    app.results[app.selectedResultIndex]?.content || '',
    'image/svg+xml'
  );
  return {
    mounted: collect(mounted),
    result: collect(result),
    override: app.labelVisibilityOverrides[targetId]
  };
}, featureId);

const findCurrentLabelTarget = (page, { leaderCount = null, embedded = false } = {}) => (
  page.evaluate(({ requiredLeaders, requireEmbedded }) => {
    const app = window.__GBDRAW_APP__;
    const svg = document.querySelector('.origin-top svg');
    const labels = Array.from(svg.querySelectorAll(
      'text[data-gbdraw-label-binding-schema="1"][data-label-feature-id]'
    ));
    for (const text of labels) {
      const featureId = text.getAttribute('data-label-feature-id');
      const lines = Array.from(svg.querySelectorAll('line[data-label-feature-id]'))
        .filter((line) => line.getAttribute('data-label-feature-id') === featureId);
      const isEmbedded = Boolean(text.querySelector('textPath'));
      if (isEmbedded !== requireEmbedded) continue;
      if (requiredLeaders !== null && lines.length !== requiredLeaders) continue;
      if (!app.getEditableLabelByFeatureId(featureId)) continue;
      if (!app.extractedFeatures.some((feature) => feature.svg_id === featureId)) continue;
      return { featureId, leaderCount: lines.length, text: text.textContent || '' };
    }
    return null;
  }, { requiredLeaders: leaderCount, requireEmbedded: embedded })
);

const setVisibility = (page, featureId, mode) => page.evaluate(async ({ targetId, nextMode }) => {
  const app = window.__GBDRAW_APP__;
  const feature = app.extractedFeatures.find((candidate) => candidate.svg_id === targetId);
  if (!feature) throw new Error(`Feature ${targetId} is unavailable.`);
  await app.openFeatureEditorFromList(feature, null);
  app.clickedFeature.labelVisibility = nextMode;
  await app.updateClickedFeatureLabelText();
}, { targetId: featureId, nextMode: mode });

const expectUnitHidden = (unit, leaderCount) => {
  expect(unit.mounted).toHaveLength(leaderCount + 1);
  expect(unit.result).toHaveLength(leaderCount + 1);
  for (const parts of [unit.mounted, unit.result]) {
    expect(parts.filter(({ tag }) => tag === 'line')).toHaveLength(leaderCount);
    expect(parts.every(({ display }) => display === 'none')).toBe(true);
    expect(parts.every(({ preview }) => preview === 'off')).toBe(true);
  }
};

const expectUnitVisible = (unit, leaderCount) => {
  expect(unit.mounted).toHaveLength(leaderCount + 1);
  expect(unit.result).toHaveLength(leaderCount + 1);
  for (const parts of [unit.mounted, unit.result]) {
    expect(parts.filter(({ tag }) => tag === 'line')).toHaveLength(leaderCount);
    expect(parts.every(({ display }) => display !== 'none')).toBe(true);
    expect(parts.every(({ preview }) => preview === null)).toBe(true);
  }
};

const labelGeometry = (page) => page.evaluate(() => Array.from(
  document.querySelectorAll('.origin-top svg [data-label-feature-id]')
).map((element) => ({
  tag: element.localName,
  id: element.getAttribute('data-label-feature-id'),
  d: element.getAttribute('d'),
  x: element.getAttribute('x'),
  y: element.getAttribute('y'),
  x1: element.getAttribute('x1'),
  y1: element.getAttribute('y1'),
  x2: element.getAttribute('x2'),
  y2: element.getAttribute('y2'),
  transform: element.getAttribute('transform')
})));

test('fresh Linear visibility is atomic through reflow, Save/Load, regeneration, and export', async ({
  browser
}, testInfo) => {
  test.setTimeout(600000);
  const page = await load(browser, seeds.linear);
  let fresh;
  try {
    await generate(page);
    await page.evaluate(() => window.__GBDRAW_APP__.syncLabelEditor());
    const target = await findCurrentLabelTarget(page, { leaderCount: 1 });
    expect(target).not.toBeNull();
    await page.evaluate(() => { window.__GBDRAW_APP__.autoLabelReflowEnabled = false; });

    const geometryBefore = await labelGeometry(page);
    const countersBefore = await page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      return {
        ordinary: state.labelReflowRequestSeq.value,
        forced: state.labelReflowForceRequestSeq.value
      };
    });
    await setVisibility(page, target.featureId, 'off');
    expectUnitHidden(await inspectUnit(page, target.featureId), 1);
    expect(await labelGeometry(page)).toEqual(geometryBefore);
    expect(await page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      return {
        ordinary: state.labelReflowRequestSeq.value,
        forced: state.labelReflowForceRequestSeq.value
      };
    })).toEqual(countersBefore);

    const saved = testInfo.outputPath('issue-564-hidden-linear.gbdraw-session.json.gz');
    await download(page, 'Save Session', saved);
    fresh = await load(browser, saved);
    await fresh.evaluate(() => window.__GBDRAW_APP__.syncLabelEditor());
    expectUnitHidden(await inspectUnit(fresh, target.featureId), 1);
    await fresh.evaluate(() => { window.__GBDRAW_APP__.autoLabelReflowEnabled = false; });
    await setVisibility(fresh, target.featureId, 'default');
    expectUnitVisible(await inspectUnit(fresh, target.featureId), 1);
    await setVisibility(fresh, target.featureId, 'off');
    await generate(fresh);
    const regenerated = await inspectUnit(fresh, target.featureId);
    expect([0, 2]).toContain(regenerated.mounted.length);
    expect(regenerated.mounted.every(({ display }) => display === 'none')).toBe(true);
    expect(regenerated.result.every(({ display }) => display === 'none')).toBe(true);

    const exported = (await download(
      fresh,
      'SVG',
      testInfo.outputPath('issue-564-hidden-linear.svg')
    )).toString();
    const exportUnit = await fresh.evaluate(({ content, targetId }) => {
      const svg = new DOMParser().parseFromString(content, 'image/svg+xml');
      return Array.from(svg.querySelectorAll('[data-label-feature-id]'))
        .filter((element) => element.getAttribute('data-label-feature-id') === targetId)
        .map((element) => ({
          tag: element.localName,
          display: element.getAttribute('display'),
          preview: element.getAttribute('data-gbdraw-label-visibility-preview'),
          schema: element.getAttribute('data-gbdraw-label-binding-schema')
        }));
    }, { content: exported, targetId: target.featureId });
    expect([0, 2]).toContain(exportUnit.length);
    expect(exportUnit.every(({ display }) => display === 'none')).toBe(true);
    expect(exportUnit.every(({ preview }) => preview === 'off')).toBe(true);
    await fresh.screenshot({
      path: testInfo.outputPath('issue-564-hidden-linear.png'),
      fullPage: true
    });

    await setVisibility(page, target.featureId, 'default');
    const restored = await inspectUnit(page, target.featureId);
    expectUnitVisible(restored, 1);
    expect(restored.override).toBeUndefined();

    const autoCounters = await page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      window.__GBDRAW_APP__.autoLabelReflowEnabled = true;
      return {
        ordinary: state.labelReflowRequestSeq.value,
        forced: state.labelReflowForceRequestSeq.value
      };
    });
    await setVisibility(page, target.featureId, 'off');
    expectUnitHidden(await inspectUnit(page, target.featureId), 1);
    expect(await page.evaluate(() => window.__GBDRAW_APP__.labelReflowProcessing)).toBe(true);
    await expect.poll(() => page.evaluate(() => ({
      processing: window.__GBDRAW_APP__.labelReflowProcessing,
      error: window.__GBDRAW_APP__.labelReflowLastError
    })), { timeout: 180000 }).toEqual({ processing: false, error: null });
    const postReflow = await inspectUnit(page, target.featureId);
    expect([0, 2]).toContain(postReflow.mounted.length);
    expect(postReflow.mounted.every(({ display }) => display === 'none')).toBe(true);
    expect(postReflow.result.every(({ display }) => display === 'none')).toBe(true);
    expect(postReflow.override).toBe('off');
    expect(await page.evaluate(async () => {
      const { state } = await import('./js/state.js');
      return {
        ordinary: state.labelReflowRequestSeq.value,
        forced: state.labelReflowForceRequestSeq.value
      };
    })).toEqual({ ordinary: autoCounters.ordinary + 1, forced: autoCounters.forced });
    expect(page.externalRequests).toEqual([]);
    expect(fresh.externalRequests).toEqual([]);
  } finally {
    await page.context().close();
    if (fresh) await fresh.context().close();
  }
});

test('fresh Circular two-segment leader visibility is atomic without forced regeneration', async ({
  browser
}, testInfo) => {
  test.setTimeout(300000);
  const page = await load(browser, seeds.circular);
  try {
    await generate(page);
    await page.evaluate(() => window.__GBDRAW_APP__.syncLabelEditor());
    await page.evaluate(() => { window.__GBDRAW_APP__.autoLabelReflowEnabled = false; });
    const target = await findCurrentLabelTarget(page, { leaderCount: 2 });
    expect(target).not.toBeNull();
    const before = await page.evaluate(async () => (await import('./js/state.js'))
      .state.labelReflowForceRequestSeq.value);
    await setVisibility(page, target.featureId, 'off');
    expectUnitHidden(await inspectUnit(page, target.featureId), 2);
    await setVisibility(page, target.featureId, 'on');
    expectUnitVisible(await inspectUnit(page, target.featureId), 2);
    expect(await page.evaluate(async () => (await import('./js/state.js'))
      .state.labelReflowForceRequestSeq.value)).toBe(before);
    await page.screenshot({
      path: testInfo.outputPath('issue-564-restored-circular.png'),
      fullPage: true
    });
    expect(page.externalRequests).toEqual([]);
  } finally {
    await page.context().close();
  }
});

const prepareLegacyTarget = (page) => page.evaluate(async () => {
  const app = window.__GBDRAW_APP__;
  app.syncLabelEditor();
  const entry = app.editableLabels.find((candidate) => candidate.kind === 'regular'
    && candidate.featureId
    && app.extractedFeatures.some((feature) => feature.svg_id === candidate.featureId));
  if (!entry) throw new Error('No legacy external label target was found.');
  const feature = app.extractedFeatures.find((candidate) => candidate.svg_id === entry.featureId);
  await app.openFeatureEditorFromList(feature, null);
  const { state } = await import('./js/state.js');
  return {
    featureId: entry.featureId,
    mounted: document.querySelector('.origin-top svg').outerHTML,
    result: app.results[app.selectedResultIndex].content,
    forceSeq: state.labelReflowForceRequestSeq.value
  };
});

test('tracked metadata-free Session refreshes once and never exposes a partial visual unit', async ({
  browser
}) => {
  test.setTimeout(300000);
  const page = await load(browser, seeds.linear);
  try {
    const legacy = await prepareLegacyTarget(page);
    expect(legacy.mounted).not.toContain('data-gbdraw-label-binding-schema="1"');
    await page.evaluate(() => {
      window.__GBDRAW_APP__.autoLabelReflowEnabled = false;
      window.__GBDRAW_APP__.clickedFeature.labelVisibility = 'off';
      window.__ISSUE564_LEGACY_APPLY__ = window.__GBDRAW_APP__.updateClickedFeatureLabelText();
    });
    await expect.poll(() => page.evaluate(() => window.__GBDRAW_APP__.labelReflowProcessing))
      .toBe(true);
    const during = await page.evaluate(async () => {
      const app = window.__GBDRAW_APP__;
      const { state } = await import('./js/state.js');
      return {
        mounted: document.querySelector('.origin-top svg').outerHTML,
        result: app.results[app.selectedResultIndex].content,
        override: app.labelVisibilityOverrides[app.clickedFeature.svg_id],
        forceSeq: state.labelReflowForceRequestSeq.value
      };
    });
    expect(during.mounted).toBe(legacy.mounted);
    expect(during.result).toBe(legacy.result);
    expect(during.override).toBe('off');
    expect(during.forceSeq).toBe(legacy.forceSeq + 1);

    await expect.poll(() => page.evaluate(() => ({
      processing: window.__GBDRAW_APP__.labelReflowProcessing,
      error: window.__GBDRAW_APP__.labelReflowLastError
    })), { timeout: 180000 }).toEqual({ processing: false, error: null });
    await page.evaluate(() => window.__ISSUE564_LEGACY_APPLY__);
    const postRefresh = await page.evaluate((targetId) => {
      const app = window.__GBDRAW_APP__;
      const parts = Array.from(document.querySelectorAll(
        '.origin-top [data-label-feature-id]'
      )).filter((element) => element.getAttribute('data-label-feature-id') === targetId);
      return {
        markerCount: document.querySelectorAll(
          '.origin-top text[data-gbdraw-label-binding-schema="1"]'
        ).length,
        override: app.labelVisibilityOverrides[targetId],
        parts: parts.map((element) => ({
          tag: element.localName,
          display: element.getAttribute('display')
        }))
      };
    }, legacy.featureId);
    expect(postRefresh.markerCount).toBeGreaterThan(0);
    expect(postRefresh.override).toBe('off');
    expect(postRefresh.parts.length === 0 || postRefresh.parts.length > 1).toBe(true);
    expect(postRefresh.parts.every(({ display }) => display === 'none')).toBe(true);
    const refreshed = await inspectUnit(page, legacy.featureId);
    expect(refreshed.mounted.length === 0 || refreshed.mounted.length > 1).toBe(true);
    expect(refreshed.mounted.every(({ display }) => display === 'none')).toBe(true);
    expect(refreshed.result.every(({ display }) => display === 'none')).toBe(true);
    expect(page.externalRequests).toEqual([]);
  } finally {
    await page.context().close();
  }
});

test('metadata-free refresh failure retains the old visual and canonical override', async ({
  browser
}) => {
  test.setTimeout(300000);
  const page = await load(browser, seeds.linear);
  try {
    const legacy = await prepareLegacyTarget(page);
    await page.evaluate(() => {
      const hooks = window.__GBDRAW_TEST_HOOKS__ || (window.__GBDRAW_TEST_HOOKS__ = {});
      hooks.beforeDiagramGenerationResponse = async () => {
        throw new Error('Forced Issue 564 label refresh failure.');
      };
      window.__GBDRAW_APP__.autoLabelReflowEnabled = false;
      window.__GBDRAW_APP__.clickedFeature.labelVisibility = 'off';
    });
    await page.evaluate(() => window.__GBDRAW_APP__.updateClickedFeatureLabelText());
    await expect.poll(() => page.evaluate(() => ({
      processing: window.__GBDRAW_APP__.labelReflowProcessing,
      error: window.__GBDRAW_APP__.labelReflowLastError
    })), { timeout: 180000 }).toEqual({
      processing: false,
      error: 'Forced Issue 564 label refresh failure.'
    });
    const failed = await page.evaluate(async (targetId) => {
      const app = window.__GBDRAW_APP__;
      const { state } = await import('./js/state.js');
      return {
        mounted: document.querySelector('.origin-top svg').outerHTML,
        result: app.results[app.selectedResultIndex].content,
        override: app.labelVisibilityOverrides[targetId],
        forceSeq: state.labelReflowForceRequestSeq.value
      };
    }, legacy.featureId);
    expect(failed.mounted).toBe(legacy.mounted);
    expect(failed.result).toBe(legacy.result);
    expect(failed.override).toBe('off');
    expect(failed.forceSeq).toBe(legacy.forceSeq + 1);
    expect(page.externalRequests).toEqual([]);
  } finally {
    await page.context().close();
  }
});

const sameLabelGenbank = (recordId, base) => {
  const sequence = base.repeat(160);
  const origin = sequence.match(/.{1,60}/g).map((chunk, index) => (
    `${String(index * 60 + 1).padStart(9)} ${chunk.match(/.{1,10}/g).join(' ')}`
  )).join('\n');
  return `LOCUS       ${recordId.padEnd(24)} ${sequence.length} bp    DNA     linear   UNA 01-JAN-2000
DEFINITION  Issue 564 same-label record.
ACCESSION   ${recordId}
VERSION     ${recordId}
KEYWORDS    .
SOURCE      synthetic construct
  ORGANISM  synthetic construct
            .
FEATURES             Location/Qualifiers
     CDS             40..50
                     /gene="shared"
                     /product="same label test"
     CDS             160..450
                     /gene="inside"
                     /product="embedded label test"
ORIGIN
${origin}
//
`;
};

test('same label text in two Linear records changes only the selected exact identity', async ({
  browser
}) => {
  test.setTimeout(300000);
  const context = await browser.newContext({ viewport: { width: 1600, height: 1000 } });
  const page = await context.newPage();
  const external = [];
  await context.route('**/*', (route) => {
    if (new URL(route.request().url()).hostname === '127.0.0.1') return route.continue();
    external.push(route.request().url());
    return route.abort();
  });
  page.on('dialog', (dialog) => dialog.dismiss());
  try {
    await openApp(page);
    await page.evaluate(async ({ first, second }) => {
      const app = window.__GBDRAW_APP__;
      app.mode = 'linear';
      app.lInputType = 'gb';
      app.addLinearSeq();
      [first, second].forEach((content, index) => app.setLinearSeqPrimaryFile(
        index,
        'gb',
        new File([content], `same-${index + 1}.gbk`, {
          type: 'text/plain',
          lastModified: index + 1
        })
      ));
      Object.assign(app.form, {
        legend: 'none',
        show_gc: false,
        show_skew: false,
        show_depth: false,
        show_labels_linear: 'all'
      });
      app.adv.label_rendering = 'auto';
      app.autoLabelReflowEnabled = false;
      await app.setLinearComparisonGlobalAction('none');
    }, {
      first: sameLabelGenbank('Issue564A', 'atg'),
      second: sameLabelGenbank('Issue564B', 'gct')
    });
    expect(await page.evaluate(() => window.__GBDRAW_APP__.runAnalysis())).toEqual({ status: 'ok' });
    await page.evaluate(() => window.__GBDRAW_APP__.syncLabelEditor());
    const targets = await page.evaluate(() => {
      const app = window.__GBDRAW_APP__;
      const svg = document.querySelector('.origin-top svg');
      const featureIds = new Set(app.extractedFeatures
        .filter((feature) => feature.qualifiers?.gene?.[0] === 'shared')
        .map((feature) => feature.svg_id));
      return Array.from(svg.querySelectorAll(
        'text[data-gbdraw-label-binding-schema="1"][data-label-feature-id]'
      )).filter((text) => featureIds.has(text.getAttribute('data-label-feature-id'))).map((text) => ({
        featureId: text.getAttribute('data-label-feature-id'),
        text: text.textContent,
        leaderCount: Array.from(svg.querySelectorAll('line[data-label-feature-id]'))
          .filter((line) => line.getAttribute('data-label-feature-id')
            === text.getAttribute('data-label-feature-id')).length
      }));
    });
    expect(targets).toHaveLength(2);
    expect(new Set(targets.map(({ featureId }) => featureId)).size).toBe(2);
    expect(new Set(targets.map(({ text }) => text)).size).toBe(1);
    expect(targets.map(({ leaderCount }) => leaderCount)).toEqual([1, 1]);
    await setVisibility(page, targets[0].featureId, 'off');
    expectUnitHidden(await inspectUnit(page, targets[0].featureId), 1);
    expectUnitVisible(await inspectUnit(page, targets[1].featureId), 1);
    const embedded = await findCurrentLabelTarget(page, { leaderCount: 0 });
    expect(embedded).not.toBeNull();
    await setVisibility(page, embedded.featureId, 'off');
    expectUnitHidden(await inspectUnit(page, embedded.featureId), 0);
    await setVisibility(page, embedded.featureId, 'default');
    expectUnitVisible(await inspectUnit(page, embedded.featureId), 0);
    expect(external).toEqual([]);
  } finally {
    await context.close();
  }
});
