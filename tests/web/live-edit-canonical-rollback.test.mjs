import assert from 'node:assert/strict';
import test from 'node:test';

import { createPreviewRuntime } from '../../gbdraw/web/js/app/preview-runtime.js';
import { createFeatureColorActions } from '../../gbdraw/web/js/app/feature-editor/color-actions.js';
import { createFeatureRuleActions } from '../../gbdraw/web/js/app/feature-editor/rule-actions.js';
import { createFeatureLabelActions } from '../../gbdraw/web/js/app/feature-editor/label-actions.js';
import { createFeatureVisibilityActions } from '../../gbdraw/web/js/app/feature-editor/visibility-actions.js';
import { createLegendEntryActions } from '../../gbdraw/web/js/app/legend/entry-actions.js';
import { createLegendStrokeActions } from '../../gbdraw/web/js/app/legend/stroke-actions.js';
import { createHistoryManager } from '../../gbdraw/web/js/services/history.js';
import { createHistoryFileStore } from '../../gbdraw/web/js/services/history-files.js';
import { createHistorySnapshotService } from '../../gbdraw/web/js/services/history-snapshot.js';
import {
  FakeSvgElement,
  appendEditableLabel,
  appendFeature,
  appendFeatureLegend,
  serializeFakeSvg
} from './helpers/editor-svg-fixture.mjs';

globalThis.CSS = { escape: (value) => String(value) };

const ref = (value) => ({ value });
const clone = (value) => JSON.parse(JSON.stringify(value));

const createHistoryHarness = async (buildIntent = () => ({ stable: true })) => {
  const history = createHistoryManager({
    buildIntent,
    applyIntent: () => {},
    buildCheckpoint: () => ({ results: [] }),
    applyCheckpoint: () => {}
  });
  await history.captureBaseline('transaction baseline');
  return history;
};

const assertFailedHistoryAction = async (label, action, buildIntent) => {
  const history = await createHistoryHarness(buildIntent);
  await assert.rejects(history.runUndoable(label, action), /serialization failed/);
  assert.equal(history.getUndoCount(), 0);
  assert.equal(history.getRedoCount(), 0);
};

const createRuntimeFixture = ({ svg, serializeSvg = null } = {}) => {
  const result = { name: 'result.svg', content: '<svg data-owner="before"></svg>' };
  const state = {
    results: ref([result]),
    selectedResultIndex: ref(0),
    resultGenerationKey: ref('generation-1'),
    skipCaptureBaseConfig: ref(false),
    svgContainer: ref({ querySelector: (selector) => selector === 'svg' ? svg : null })
  };
  const metrics = [];
  globalThis.__GBDRAW_TEST_HOOKS__ = {
    onStructuralMetric: (metric) => metrics.push(metric)
  };
  const runtime = createPreviewRuntime({
    state,
    serializeSvg: serializeSvg || (() => {
      throw new Error('serialization failed');
    })
  });
  runtime.mountResultSvg(0, svg);
  const retainedFeatureIndex = new Map([['retained', []]]);
  const retainedLegendIndex = new Map([['retained', []]]);
  runtime.getActiveRuntime().indexes.features = retainedFeatureIndex;
  runtime.getActiveRuntime().indexes.legend = retainedLegendIndex;
  return {
    metrics,
    result,
    retainedFeatureIndex,
    retainedLegendIndex,
    runtime,
    state
  };
};

const assertRuntimeRestored = (fixture) => {
  const active = fixture.runtime.getActiveRuntime();
  assert.equal(fixture.state.results.value[0], fixture.result);
  assert.equal(fixture.result.content, '<svg data-owner="before"></svg>');
  assert.equal(active.dirty, false);
  assert.deepEqual([...active.dirtyReasons], []);
  assert.equal(active.indexes.features, fixture.retainedFeatureIndex);
  assert.equal(active.indexes.legend, fixture.retainedLegendIndex);
  assert.equal(fixture.state.skipCaptureBaseConfig.value, false);
};

const createColorHarness = ({
  addLegendEntry = null,
  applySpecificRulesToSvg = null,
  onLegendGeometryChanged = null
} = {}) => {
  const svg = new FakeSvgElement('svg');
  const featureElement = appendFeature(svg, 'feature-a');
  const feature = {
    id: 'feature-a',
    svg_id: 'feature-a',
    type: 'CDS',
    product: 'Feature A',
    start: 1,
    end: 10
  };
  const fixture = createRuntimeFixture({ svg });
  const manualSpecificRules = [];
  const featureColorOverrides = {};
  const featureStrokeOverrides = {};
  const legendColorOverrides = {};
  const legendStrokeOverrides = {};
  const legendEntries = ref([]);
  const clickedOwner = {
    svg_id: feature.svg_id,
    feat: feature,
    color: '#111111',
    legendName: 'Feature A',
    strokeColor: '#222222',
    strokeWidth: 1,
    originalStrokeColor: '#222222',
    originalStrokeWidth: 1
  };
  const clickedFeature = ref(clickedOwner);
  const state = {
    ...fixture.state,
    appliedPaletteColors: ref({ CDS: '#cccccc' }),
    manualSpecificRules,
    extractedFeatures: ref([feature]),
    biologicalFeatures: ref([feature]),
    featureColorOverrides,
    svgContainer: fixture.state.svgContainer,
    clickedFeature,
    featureStyleScopeDialog: {},
    resetColorDialog: {},
    legendRenameDialog: {},
    legendEntries,
    legendStrokeOverrides,
    legendColorOverrides,
    originalLegendOrder: ref([]),
    originalLegendColors: ref({}),
    originalSvgStroke: ref({ color: '#222222', width: 1 }),
    featureStrokeOverrides,
    skipCaptureBaseConfig: fixture.state.skipCaptureBaseConfig,
    skipExtractOnSvgChange: ref(false),
    addedLegendCaptions: ref(new Set()),
    semanticFileWatchersSuppressed: ref(false)
  };
  const actions = createFeatureColorActions({
    state,
    nextTick: async () => {},
    legendActions: {
      addLegendEntry: addLegendEntry || (async (caption) => caption),
      removeLegendEntry: () => false,
      updateLegendEntryColorByCaption: () => false,
      compactLegendEntries: () => {},
      extractLegendEntries: () => {},
      onLegendGeometryChanged: onLegendGeometryChanged || (() => {}),
      getAllFeatureLegendGroups: () => []
    },
    svgActions: { applySpecificRulesToSvg: applySpecificRulesToSvg || (() => {}) },
    ruleActions: {
      countFeaturesMatchingRule: () => 0,
      findExistingColorForCaption: () => null,
      findFeaturesWithSameDisplayedLabel: () => [],
      findFeaturesWithSameIndividualLabel: () => [],
      findFeaturesWithSameLegendItem: () => [],
      findMatchingRegexRule: () => null,
      getDisplayedFeatureLabel: () => feature.product,
      getEffectiveLegendCaption: () => feature.product,
      getIndividualFeatureLabel: () => feature.product,
      getFeatureQualifier: () => ({ qual: 'hash', val: feature.svg_id }),
      getLabelSpecificRule: () => null
    },
    featureSvgActions: {
      getFeatureElements: (_svg, featureId) => featureId === feature.svg_id ? [featureElement] : [],
      getFeatureFillElements: (_svg, featureId) => featureId === feature.svg_id ? [featureElement] : []
    },
    previewRuntime: fixture.runtime
  });
  return {
    actions,
    addedLegendCaptions: state.addedLegendCaptions,
    clickedFeature,
    clickedOwner,
    feature,
    featureColorOverrides,
    featureElement,
    featureStrokeOverrides,
    fixture,
    legendColorOverrides,
    legendEntries,
    legendStrokeOverrides,
    manualSpecificRules,
    state,
    svg
  };
};

test('Feature fill production action rolls back canonical state and mounted artifact on serialization failure', async () => {
  const harness = createColorHarness();
  const beforeDom = serializeFakeSvg(harness.svg);
  const beforeClicked = clone(harness.clickedOwner);

  await assertFailedHistoryAction(
    'Feature fill',
    () => harness.actions.setFeatureColorValue(harness.feature, 'none', 'No fill'),
    () => ({ overrides: clone(harness.featureColorOverrides) })
  );

  assert.deepEqual(harness.featureColorOverrides, {});
  assert.deepEqual(harness.manualSpecificRules, []);
  assert.equal(harness.clickedFeature.value, harness.clickedOwner);
  assert.deepEqual(harness.clickedOwner, beforeClicked);
  assert.equal(serializeFakeSvg(harness.svg), beforeDom);
  assertRuntimeRestored(harness.fixture);
  assert.equal(
    harness.fixture.metrics.filter((metric) => metric.name === 'liveEditCanonicalRollbackCount').length,
    1
  );
});

test('Async Feature fill add branch rolls back Legend, rules, DOM, and Result as one transaction', async () => {
  let harness = null;
  let addedLegendNode = null;
  harness = createColorHarness({
    addLegendEntry: async (caption, _color, { lease } = {}) => {
      assert.ok(lease, 'the async Legend helper receives the action lease');
      await Promise.resolve();
      assert.equal(
        harness.state.semanticFileWatchersSuppressed.value,
        true,
        'reactive semantic watchers stay suppressed across the helper await'
      );
      addedLegendNode = new FakeSvgElement('g', { 'data-legend-key': caption });
      lease.mutate(({ mutation }) => mutation.appendChild(harness.svg, addedLegendNode));
      return caption;
    },
    applySpecificRulesToSvg: ({ lease } = {}) => harness.fixture.runtime.applyFeatureFillChanges(
      [{ featureId: harness.feature.svg_id, color: '#abcdef' }],
      { reason: 'async-feature-fill-rule', lease }
    )
  });
  const beforeDom = serializeFakeSvg(harness.svg);
  const beforeClicked = clone(harness.clickedOwner);

  await assertFailedHistoryAction(
    'Async Feature fill',
    () => harness.actions.applyColorToSelectedFeatures(
      [harness.feature],
      '#abcdef',
      'Async Legend'
    ),
    () => ({
      rules: clone(harness.manualSpecificRules),
      overrides: clone(harness.featureColorOverrides)
    })
  );

  assert.deepEqual(harness.manualSpecificRules, []);
  assert.deepEqual(harness.featureColorOverrides, {});
  assert.deepEqual([...harness.addedLegendCaptions.value], []);
  assert.deepEqual(harness.clickedOwner, beforeClicked);
  assert.equal(harness.state.semanticFileWatchersSuppressed.value, false);
  assert.equal(addedLegendNode?.parentNode, null);
  assert.equal(serializeFakeSvg(harness.svg), beforeDom);
  assertRuntimeRestored(harness.fixture);
});

test('Specific color rule add owns async Legend extraction and mounted color in one transaction', async () => {
  const svg = new FakeSvgElement('svg');
  const featureElement = appendFeature(svg, 'feature-a');
  const fixture = createRuntimeFixture({ svg });
  const feature = {
    id: 'feature-a',
    svg_id: 'feature-a',
    type: 'CDS',
    product: 'Feature A',
    start: 1,
    end: 10
  };
  const manualSpecificRules = [];
  const featureColorOverrides = {};
  const legendEntriesOwner = [];
  const originalLegendOrderOwner = [];
  const originalLegendColorsOwner = {};
  const state = {
    ...fixture.state,
    currentColors: ref({}),
    appliedPaletteColors: ref({ CDS: '#cccccc' }),
    newColorFeat: ref(''),
    newColorVal: ref('#000000'),
    manualSpecificRules,
    newSpecRule: {
      feat: 'CDS',
      qual: 'product',
      val: '^Feature A$',
      color: '#abcdef',
      cap: 'Rule Legend'
    },
    specificRulePresets: [],
    selectedSpecificPreset: ref(''),
    specificRulePresetLoading: ref(false),
    manualPriorityRules: [],
    newPriorityRule: {},
    adv: { features: [], feature_shapes: {}, legend_box_size: 14, legend_font_size: 14 },
    newFeatureToAdd: ref(''),
    extractedFeatures: ref([feature]),
    featureColorOverrides,
    editableLabels: ref([]),
    labelTextFeatureOverrides: {},
    labelTextBulkOverrides: {},
    addedLegendCaptions: ref(new Set()),
    fileLegendCaptions: ref(new Set()),
    semanticFileWatchersSuppressed: ref(false),
    legendEntries: ref(legendEntriesOwner),
    originalLegendOrder: ref(originalLegendOrderOwner),
    originalLegendColors: ref(originalLegendColorsOwner)
  };
  let addedLegendNode = null;
  const actions = createFeatureRuleActions({
    state,
    nextTick: async () => {},
    legendActions: {
      addLegendEntry: async (caption, _color, { lease } = {}) => {
        assert.ok(lease);
        await Promise.resolve();
        assert.equal(state.semanticFileWatchersSuppressed.value, true);
        addedLegendNode = new FakeSvgElement('g', { 'data-legend-key': caption });
        lease.mutate(({ mutation }) => mutation.appendChild(svg, addedLegendNode));
        return caption;
      },
      removeLegendEntry: () => false,
      extractLegendEntries: () => {
        state.legendEntries.value = [{ caption: 'Rule Legend', color: '#abcdef' }];
        state.originalLegendOrder.value = ['Rule Legend'];
        state.originalLegendColors.value['Rule Legend'] = '#abcdef';
      },
      onLegendGeometryChanged: () => {}
    },
    svgActions: {
      applyPaletteToSvg: () => false,
      applySpecificRulesToSvg: ({ lease } = {}) => fixture.runtime.applyFeatureFillChanges(
        [{ featureId: feature.svg_id, color: '#abcdef' }],
        { reason: 'specific-rule-fill', lease }
      )
    },
    previewRuntime: fixture.runtime
  });
  const beforeDom = serializeFakeSvg(svg);

  await assertFailedHistoryAction(
    'Add specific color rule',
    () => actions.addSpecificRule(),
    () => ({ rules: clone(manualSpecificRules) })
  );

  assert.deepEqual(manualSpecificRules, []);
  assert.equal(state.legendEntries.value, legendEntriesOwner);
  assert.equal(state.originalLegendOrder.value, originalLegendOrderOwner);
  assert.equal(state.originalLegendColors.value, originalLegendColorsOwner);
  assert.deepEqual(legendEntriesOwner, []);
  assert.deepEqual(originalLegendOrderOwner, []);
  assert.deepEqual(originalLegendColorsOwner, {});
  assert.deepEqual([...state.addedLegendCaptions.value], []);
  assert.deepEqual([...state.fileLegendCaptions.value], []);
  assert.equal(state.newSpecRule.val, '^Feature A$');
  assert.equal(state.semanticFileWatchersSuppressed.value, false);
  assert.equal(addedLegendNode?.parentNode, null);
  assert.equal(serializeFakeSvg(svg), beforeDom);
  assertRuntimeRestored(fixture);
});

test('Specific-rule preset aborts when its leased Result becomes stale during fetch', async () => {
  const svg = new FakeSvgElement('svg');
  const { entry } = appendFeatureLegend(svg, 'Existing', '#112233');
  const fixture = createRuntimeFixture({
    svg,
    serializeSvg: () => serializeFakeSvg(svg)
  });
  const existingRule = {
    feat: 'CDS',
    qual: 'product',
    val: '^Existing$',
    color: '#112233',
    cap: 'Existing'
  };
  const manualSpecificRules = [existingRule];
  const legendEntriesOwner = [{ caption: 'Existing', color: '#112233' }];
  const originalLegendOrderOwner = ['Existing'];
  const originalLegendColorsOwner = { Existing: '#112233' };
  const state = {
    ...fixture.state,
    currentColors: ref({ CDS: '#cccccc' }),
    appliedPaletteColors: ref({ CDS: '#cccccc' }),
    newColorFeat: ref(''),
    newColorVal: ref('#000000'),
    manualSpecificRules,
    newSpecRule: { feat: 'CDS', qual: 'product', val: '', color: '#000000', cap: '' },
    specificRulePresets: [{ id: 'preset', path: '/preset.tsv' }],
    selectedSpecificPreset: ref('preset'),
    specificRulePresetLoading: ref(false),
    manualPriorityRules: [],
    newPriorityRule: {},
    adv: { features: [], feature_shapes: {}, legend_box_size: 14, legend_font_size: 14 },
    newFeatureToAdd: ref(''),
    extractedFeatures: ref([]),
    featureColorOverrides: {},
    editableLabels: ref([]),
    labelTextFeatureOverrides: {},
    labelTextBulkOverrides: {},
    addedLegendCaptions: ref(new Set(['Existing'])),
    fileLegendCaptions: ref(new Set(['Existing'])),
    semanticFileWatchersSuppressed: ref(false),
    legendEntries: ref(legendEntriesOwner),
    originalLegendOrder: ref(originalLegendOrderOwner),
    originalLegendColors: ref(originalLegendColorsOwner)
  };
  const actions = createFeatureRuleActions({
    state,
    nextTick: async () => {},
    legendActions: {
      addLegendEntry: async (caption) => caption,
      removeLegendEntry: (_caption, { lease } = {}) => {
        assert.ok(lease);
        return lease.mutate(({ mutation }) => mutation.removeNode(entry));
      },
      extractLegendEntries: () => {
        state.legendEntries.value = [];
        state.originalLegendOrder.value = [];
        state.originalLegendColors.value = {};
      },
      onLegendGeometryChanged: () => {}
    },
    svgActions: {
      applyPaletteToSvg: () => false,
      applySpecificRulesToSvg: () => false
    },
    previewRuntime: fixture.runtime
  });
  let releaseFetch;
  let markFetchStarted;
  const fetchStarted = new Promise((resolve) => { markFetchStarted = resolve; });
  const originalFetch = globalThis.fetch;
  const originalAlert = globalThis.alert;
  const originalConsoleError = console.error;
  globalThis.fetch = () => new Promise((resolve) => {
    releaseFetch = () => resolve({
      ok: true,
      text: async () => 'CDS\tproduct\t^Preset$\t#abcdef\tPreset\n'
    });
    markFetchStarted();
  });
  globalThis.alert = () => {};
  console.error = () => {};
  const beforeDom = serializeFakeSvg(svg);
  const history = await createHistoryHarness(() => ({ rules: clone(manualSpecificRules) }));
  try {
    const pending = history.runUndoable(
      'Apply specific color preset',
      () => actions.applySpecificRulePreset()
    );
    await fetchStarted;
    fixture.result.content = '<svg data-owner="newer"></svg>';
    releaseFetch();
    await assert.rejects(pending, (error) => error?.name === 'StaleDomEditError');
  } finally {
    globalThis.fetch = originalFetch;
    globalThis.alert = originalAlert;
    console.error = originalConsoleError;
  }

  assert.equal(manualSpecificRules[0], existingRule);
  assert.deepEqual(manualSpecificRules, [existingRule]);
  assert.equal(state.legendEntries.value, legendEntriesOwner);
  assert.equal(state.originalLegendOrder.value, originalLegendOrderOwner);
  assert.equal(state.originalLegendColors.value, originalLegendColorsOwner);
  assert.deepEqual([...state.addedLegendCaptions.value], ['Existing']);
  assert.deepEqual([...state.fileLegendCaptions.value], ['Existing']);
  assert.equal(state.specificRulePresetLoading.value, false);
  assert.equal(state.semanticFileWatchersSuppressed.value, false);
  assert.equal(serializeFakeSvg(svg), beforeDom);
  assert.equal(fixture.state.results.value[0], fixture.result);
  assert.equal(fixture.result.content, '<svg data-owner="newer"></svg>');
  assert.equal(fixture.runtime.getActiveRuntime().dirty, false);
  assert.equal(fixture.runtime.getActiveRuntime().indexes.features, fixture.retainedFeatureIndex);
  assert.equal(fixture.runtime.getActiveRuntime().indexes.legend, fixture.retainedLegendIndex);
  assert.equal(history.getUndoCount(), 0);
  assert.equal(history.getRedoCount(), 0);
});

test('Feature stroke production action rolls back nested override and clicked editor state', async () => {
  const harness = createColorHarness();
  const beforeDom = serializeFakeSvg(harness.svg);
  const beforeClicked = clone(harness.clickedOwner);

  await assertFailedHistoryAction(
    'Feature stroke',
    () => harness.actions.updateClickedFeatureStroke('#abcdef', 3),
    () => ({ overrides: clone(harness.featureStrokeOverrides) })
  );

  assert.deepEqual(harness.featureStrokeOverrides, {});
  assert.equal(harness.clickedFeature.value, harness.clickedOwner);
  assert.deepEqual(harness.clickedOwner, beforeClicked);
  assert.equal(serializeFakeSvg(harness.svg), beforeDom);
  assertRuntimeRestored(harness.fixture);
});

test('Feature visibility production action rolls back direct overrides, rules, and clicked state', async () => {
  const svg = new FakeSvgElement('svg');
  const featureElement = appendFeature(svg, 'feature-a');
  const fixture = createRuntimeFixture({ svg });
  const feature = { svg_id: 'feature-a', type: 'CDS', product: '' };
  const featureVisibilityManualRules = [];
  const featureVisibilityOverrides = {};
  const clickedOwner = { svg_id: feature.svg_id, feat: feature, featureVisibility: 'default' };
  const clickedFeature = ref(clickedOwner);
  const actions = createFeatureVisibilityActions({
    state: {
      ...fixture.state,
      clickedFeature,
      extractedFeatures: ref([feature]),
      orthogroups: ref([]),
      featureVisibilityManualRules,
      featureVisibilityRules: ref([]),
      featureVisibilityOverrides,
      featureVisibilitySelectorCache: {},
      featureVisibilityScopeDialog: {},
      labelLayoutDirtyReason: ref(''),
      results: fixture.state.results,
      selectedResultIndex: fixture.state.selectedResultIndex,
      resultGenerationKey: fixture.state.resultGenerationKey,
      svgContainer: fixture.state.svgContainer
    },
    featureSvgActions: {
      applyVisibilityPreviewChanges: (changes, options) => (
        fixture.runtime.applyFeatureVisibilityChanges(changes, options)
      )
    },
    previewRuntime: fixture.runtime
  });
  const beforeDom = serializeFakeSvg(svg);
  const beforeClicked = clone(clickedOwner);

  await assertFailedHistoryAction(
    'Feature visibility',
    () => actions.updateClickedFeatureVisibility('off'),
    () => ({ overrides: clone(featureVisibilityOverrides) })
  );

  assert.deepEqual(featureVisibilityOverrides, {});
  assert.deepEqual(featureVisibilityManualRules, []);
  assert.equal(clickedFeature.value, clickedOwner);
  assert.deepEqual(clickedOwner, beforeClicked);
  assert.equal(featureElement.getAttribute('display'), null);
  assert.equal(serializeFakeSvg(svg), beforeDom);
  assertRuntimeRestored(fixture);

  feature.product = 'Feature A';
  const beforeManualClicked = clone(clickedOwner);
  await assertFailedHistoryAction(
    'manual Feature visibility rule',
    () => actions.setFeatureVisibility(feature, 'off', {
      triggerReflow: false,
      scope: {
        id: 'product',
        featureType: 'CDS',
        qualifier: 'product',
        value: 'Feature A',
        label: 'Exact product: Feature A'
      }
    }),
    () => ({ rules: clone(featureVisibilityManualRules) })
  );
  assert.deepEqual(featureVisibilityOverrides, {});
  assert.deepEqual(featureVisibilityManualRules, []);
  assert.deepEqual(clickedOwner, beforeManualClicked);
  assert.equal(serializeFakeSvg(svg), beforeDom);
  assertRuntimeRestored(fixture);
});

test('Label text production action journals canonical overrides with its DOM edit', async () => {
  const svg = new FakeSvgElement('svg');
  const label = appendEditableLabel(svg, 'feature-a', 'before label');
  label.setAttribute('data-label-key', 'label-a');
  const defaultQuerySelector = svg.querySelector.bind(svg);
  svg.querySelector = (selector) => (
    selector === 'text[data-label-key="label-a"]' ? label : defaultQuerySelector(selector)
  );
  const fixture = createRuntimeFixture({ svg });
  const labelTextFeatureOverrides = {};
  const labelTextBulkOverrides = {};
  const labelTextFeatureOverrideSources = {};
  const labelVisibilityOverrides = {};
  const labelTextScopeDialog = {
    show: true,
    labelKey: 'label-a',
    sourceText: 'before label',
    newText: 'after label',
    featureId: 'feature-a'
  };
  const actions = createFeatureLabelActions({
    state: {
      ...fixture.state,
      mode: ref('linear'),
      form: { show_labels_linear: 'out' },
      filterMode: ref('none'),
      manualWhitelist: [],
      editableLabels: ref([]),
      extractedFeatures: ref([]),
      clickedFeature: ref(null),
      labelTextScopeDialog,
      labelTextFeatureOverrides,
      labelTextBulkOverrides,
      labelTextFeatureOverrideSources,
      labelVisibilityOverrides,
      labelOverrideContextKey: ref(''),
      labelOverrideBuildWarning: ref(''),
      globalLabelModeDialog: {},
      autoLabelReflowEnabled: ref(false),
      labelReflowRequestSeq: ref(0),
      labelReflowRequestReason: ref(''),
      labelReflowForceRequestSeq: ref(0),
      labelReflowForceRequestReason: ref(''),
      labelReflowLastError: ref(null)
    },
    previewRuntime: fixture.runtime
  });
  const beforeDialog = clone(labelTextScopeDialog);
  const beforeDom = serializeFakeSvg(svg);

  await assertFailedHistoryAction(
    'Label text',
    () => actions.handleLabelTextScopeChoice('single'),
    () => ({ overrides: clone(labelTextFeatureOverrides) })
  );

  assert.deepEqual(labelTextFeatureOverrides, {});
  assert.deepEqual(labelTextBulkOverrides, {});
  assert.deepEqual(labelTextFeatureOverrideSources, {});
  assert.deepEqual(labelVisibilityOverrides, {});
  assert.deepEqual(labelTextScopeDialog, beforeDialog);
  assert.equal(serializeFakeSvg(svg), beforeDom);
  assertRuntimeRestored(fixture);
});

const createLegendHarness = () => {
  const svg = new FakeSvgElement('svg');
  const legend = appendFeatureLegend(svg, 'Legend A', '#111111');
  const fixture = createRuntimeFixture({ svg });
  const legendEntries = ref([{
    caption: 'Legend A',
    originalCaption: 'Legend A',
    color: '#111111',
    featureIds: []
  }]);
  const legendColorOverrides = {};
  const legendStrokeOverrides = {};
  const state = {
    ...fixture.state,
    adv: {},
    legendEntries,
    deletedLegendEntries: ref([]),
    originalLegendOrder: ref(['Legend A']),
    originalLegendColors: ref({ 'Legend A': '#111111' }),
    newLegendCaption: ref(''),
    newLegendColor: ref('#808080'),
    legendStrokeOverrides,
    legendColorOverrides,
    manualSpecificRules: [],
    extractedFeatures: ref([]),
    featureStrokeOverrides: {},
    originalSvgStroke: ref({ color: '#222222', width: 1 })
  };
  const entryActions = createLegendEntryActions({
    state,
    layoutActions: {
      updatePairwiseLegendPositions: () => {},
      reflowDualLegendLayout: () => {},
      compactLegendEntries: () => {}
    },
    previewRuntime: fixture.runtime,
    legendHelperOperations: {
      measureText: async () => ({ width: 10 }),
      generateEntrySvg: async () => '<g></g>'
    }
  });
  const strokeActions = createLegendStrokeActions({ state, previewRuntime: fixture.runtime });
  return {
    entryActions,
    fixture,
    legend,
    legendColorOverrides,
    legendEntries,
    legendStrokeOverrides,
    state,
    strokeActions,
    svg
  };
};

test('Legend fill production action leaves canonical entry state untouched when serialization fails', async () => {
  const harness = createLegendHarness();
  const beforeEntriesOwner = harness.legendEntries.value;
  const beforeEntries = clone(beforeEntriesOwner);
  const beforeDom = serializeFakeSvg(harness.svg);

  await assertFailedHistoryAction(
    'Legend fill',
    () => harness.entryActions.updateLegendEntryColor(0, '#abcdef'),
    () => ({ entries: clone(harness.legendEntries.value) })
  );

  assert.equal(harness.legendEntries.value, beforeEntriesOwner);
  assert.deepEqual(beforeEntriesOwner, beforeEntries);
  assert.deepEqual(harness.legendColorOverrides, {});
  assert.equal(serializeFakeSvg(harness.svg), beforeDom);
  assertRuntimeRestored(harness.fixture);
});

test('Legend stroke production action rolls back its nested override and mounted swatch', async () => {
  const harness = createLegendHarness();
  const beforeEntriesOwner = harness.legendEntries.value;
  const beforeEntries = clone(beforeEntriesOwner);
  const beforeDom = serializeFakeSvg(harness.svg);

  await assertFailedHistoryAction(
    'Legend stroke',
    () => harness.strokeActions.updateLegendEntryStrokeColor(0, '#abcdef'),
    () => ({ overrides: clone(harness.legendStrokeOverrides) })
  );

  assert.equal(harness.legendEntries.value, beforeEntriesOwner);
  assert.deepEqual(beforeEntriesOwner, beforeEntries);
  assert.deepEqual(harness.legendStrokeOverrides, {});
  assert.equal(serializeFakeSvg(harness.svg), beforeDom);
  assertRuntimeRestored(harness.fixture);
});

test('file Legend sync failure preserves canonical array and entry owners', async () => {
  const harness = createLegendHarness();
  harness.legend.entry.setAttribute('data-legend-owner', 'specific-color-file');
  harness.legend.swatch.setAttribute('transform', 'translate(0,7)');
  harness.legend.entry.querySelector('text').setAttribute('transform', 'translate(22,7)');
  const entriesOwner = harness.legendEntries.value;
  const entryOwner = entriesOwner[0];
  const orderOwner = harness.state.originalLegendOrder.value;
  const colorsOwner = harness.state.originalLegendColors.value;
  const beforeEntries = clone(entriesOwner);
  const beforeOrder = clone(orderOwner);
  const beforeColors = clone(colorsOwner);
  const beforeDom = serializeFakeSvg(harness.svg);

  await assert.rejects(
    harness.entryActions.syncFileLegendEntries([
      { caption: 'Legend A', color: '#abcdef' }
    ], {
      previousFileIntents: [{ caption: 'Legend A', color: '#111111' }]
    }),
    /serialization failed/
  );

  assert.equal(harness.legendEntries.value, entriesOwner);
  assert.equal(harness.legendEntries.value[0], entryOwner);
  assert.equal(harness.state.originalLegendOrder.value, orderOwner);
  assert.equal(harness.state.originalLegendColors.value, colorsOwner);
  assert.deepEqual(entriesOwner, beforeEntries);
  assert.deepEqual(orderOwner, beforeOrder);
  assert.deepEqual(colorsOwner, beforeColors);
  assert.equal(serializeFakeSvg(harness.svg), beforeDom);
  assertRuntimeRestored(harness.fixture);
});

test('Legend add keeps its helper failure primary when journal rollback also fails', async () => {
  const svg = new FakeSvgElement('svg', {
    'data-gbdraw-composition-schema': '1',
    'data-gbdraw-composition': JSON.stringify({
      legend: {
        automaticTranslation: [0, 0],
        localBounds: { x: 0, y: 0, width: 100, height: 20 },
        role: 'legend',
        selector: '[data-gbdraw-composition-role="legend"]'
      },
      legendReflow: { colorRectSize: 14, lineHeight: 24, textXOffset: 22 },
      legendSide: 'bottom',
      overlayObstacles: [],
      overlayPolicy: {
        candidateScoreOrder: [
          'totalAnchorDistance',
          'xAnchorDistance',
          'yAnchorDistance',
          'nearEdgeX',
          'nearEdgeY'
        ],
        canvasGrowthCandidateOrder: ['horizontal', 'vertical'],
        canvasGrowthScoreOrder: ['addedArea', 'addedExtent', 'candidateOrder'],
        quadrantBoundaryRatio: 0.5
      },
      primary: {
        automaticTranslation: [0, 0],
        finalBounds: { x: 0, y: 0, width: 100, height: 100 },
        role: 'primary',
        selector: '[data-gbdraw-composition-role="primary"]'
      },
      spacing: {
        dockGapPx: 24,
        edgePaddingPx: 16,
        overlayClearancePx: 8,
        stackGapPx: 20,
        titleGapPx: 20
      },
      title: null,
      titleSide: 'none'
    })
  });
  appendFeatureLegend(svg, 'Legend A', '#111111');
  const fixture = createRuntimeFixture({ svg, serializeSvg: () => '<svg/>' });
  const originalDocument = globalThis.document;
  const originalDomParser = globalThis.DOMParser;
  globalThis.document = {
    createElementNS: () => {
      const entry = new FakeSvgElement('g');
      entry.remove = () => { throw new Error('Legend rollback removal failed'); };
      return entry;
    },
    importNode: (node) => node
  };
  globalThis.DOMParser = class {
    parseFromString() {
      return { querySelector: () => null };
    }
  };
  try {
    const state = {
      ...fixture.state,
      adv: {},
      legendEntries: ref([{ caption: 'Legend A', color: '#111111' }]),
      deletedLegendEntries: ref([]),
      originalLegendOrder: ref(['Legend A']),
      originalLegendColors: ref({ 'Legend A': '#111111' }),
      newLegendCaption: ref(''),
      newLegendColor: ref('#808080'),
      legendStrokeOverrides: {},
      legendColorOverrides: {},
      manualSpecificRules: []
    };
    const actions = createLegendEntryActions({
      state,
      layoutActions: {
        compactLegendEntries: () => {},
        reflowDualLegendLayout: () => {},
        updatePairwiseLegendPositions: () => {
          throw new Error('Legend layout helper failed');
        }
      },
      previewRuntime: fixture.runtime,
      legendHelperOperations: {
        measureText: async () => ({ result: { width: 10 } }),
        generateEntrySvg: async () => ({ result: { rect: '', text: '' } })
      }
    });

    let error = null;
    try {
      await actions.addLegendEntry('Legend B', '#abcdef', { throwOnError: true });
    } catch (caught) {
      error = caught;
    }
    assert.ok(error);
    assert.equal(error.message, 'Legend layout helper failed');
    assert.ok(Array.isArray(error.rollbackErrors));
    assert.match(
      error.rollbackErrors.map((rollbackError) => rollbackError.message).join('\n'),
      /rollback/i
    );
    assertRuntimeRestored(fixture);
  } finally {
    globalThis.document = originalDocument;
    globalThis.DOMParser = originalDomParser;
  }
});

test('stale same-object content settlement rolls back canonical visibility and DOM state', async () => {
  const svg = new FakeSvgElement('svg');
  const featureElement = appendFeature(svg, 'feature-a');
  let fixture = null;
  fixture = createRuntimeFixture({
    svg,
    serializeSvg: () => {
      fixture.result.content = '<svg data-owner="newer-same-object"></svg>';
      return '<svg data-owner="stale-settlement"></svg>';
    }
  });
  const feature = { svg_id: 'feature-a', type: 'CDS', product: 'Feature A' };
  const featureVisibilityOverrides = {};
  const clickedOwner = { svg_id: feature.svg_id, feat: feature, featureVisibility: 'default' };
  const actions = createFeatureVisibilityActions({
    state: {
      ...fixture.state,
      clickedFeature: ref(clickedOwner),
      extractedFeatures: ref([feature]),
      orthogroups: ref([]),
      featureVisibilityManualRules: [],
      featureVisibilityRules: ref([]),
      featureVisibilityOverrides,
      featureVisibilitySelectorCache: {},
      featureVisibilityScopeDialog: {},
      labelLayoutDirtyReason: ref(''),
      results: fixture.state.results,
      selectedResultIndex: fixture.state.selectedResultIndex,
      resultGenerationKey: fixture.state.resultGenerationKey,
      svgContainer: fixture.state.svgContainer
    },
    featureSvgActions: {
      applyVisibilityPreviewChanges: (changes, options) => (
        fixture.runtime.applyFeatureVisibilityChanges(changes, options)
      )
    },
    previewRuntime: fixture.runtime
  });

  const history = await createHistoryHarness(() => ({
    overrides: clone(featureVisibilityOverrides),
    clickedVisibility: clickedOwner.featureVisibility
  }));
  await history.runUndoable('stale Feature visibility', () => actions.setFeatureVisibility(
    feature,
    'off',
    { triggerReflow: false, scope: { id: 'feature' } }
  ));

  assert.deepEqual(featureVisibilityOverrides, {});
  assert.equal(clickedOwner.featureVisibility, 'default');
  assert.equal(featureElement.getAttribute('display'), null);
  assert.equal(fixture.state.results.value[0], fixture.result);
  assert.equal(fixture.result.content, '<svg data-owner="newer-same-object"></svg>');
  assert.equal(fixture.runtime.getActiveRuntime().dirty, false);
  assert.equal(history.getUndoCount(), 0);
  assert.equal(history.getRedoCount(), 0);
  assert.equal(
    fixture.metrics.filter((metric) => metric.name === 'liveEditStaleContentTokenCount').length,
    1
  );
});

test('History Undo serialization failure restores snapshot owners, DOM, Result, and stacks', async () => {
  const svg = new FakeSvgElement('svg');
  const featureElement = appendFeature(svg, 'feature-a');
  const legend = appendFeatureLegend(svg, 'Before', '#111111');
  const fixture = createRuntimeFixture({ svg });
  const beforeLegendOwner = [{ caption: 'Before', color: '#111111' }];
  const featureVisibilityOverrides = {};
  const legendEntries = ref(beforeLegendOwner);
  const clickedOwner = { featureVisibility: 'off', svg_id: 'feature-a' };
  const state = {
    ...fixture.state,
    form: {},
    adv: {},
    files: {},
    linearSeqs: [],
    linearComparisonPlan: { mode: 'adjacent', defaultSource: 'losat', edges: [] },
    semanticFileWatchersSuppressed: ref(false),
    sessionResourceDiscoveryDeferred: ref(false),
    selectedFeatureRecordIdx: ref(0),
    featureColorOverrides: {},
    featureVisibilityManualRules: [],
    featureVisibilityOverrides,
    labelTextFeatureOverrides: {},
    labelTextBulkOverrides: {},
    labelTextFeatureOverrideSources: {},
    labelVisibilityOverrides: {},
    labelOverrideContextKey: ref(''),
    legendEntries,
    deletedLegendEntries: ref([]),
    legendColorOverrides: {},
    legendStrokeOverrides: {},
    addedLegendCaptions: ref(new Set()),
    featureStrokeOverrides: {},
    selectedOrthogroupId: ref(''),
    selectedOrthogroupAlignmentFeature: ref(''),
    orthogroupNameOverrides: {},
    orthogroupDescriptionOverrides: {},
    clickedFeature: ref(clickedOwner),
    clickedPairwiseMatch: ref(null),
    clickedLabel: ref(null),
    featureStyleScopeDialog: { show: false },
    resetColorDialog: { show: false },
    legendRenameDialog: { show: false },
    labelTextScopeDialog: { show: false },
    featureVisibilityScopeDialog: { show: false },
    globalLabelModeDialog: { show: false }
  };
  const snapshots = createHistorySnapshotService({
    state,
    fileStore: createHistoryFileStore(),
    nextTick: async () => {}
  });
  snapshots.setAfterApplyHistoryIntent(async (_intent, { lease }) => {
    assert.ok(lease, 'History replay must reuse the outer mounted transaction');
    lease.mutate(({ mutation }) => {
      if (featureVisibilityOverrides['feature-a'] === 'off') {
        mutation.setAttribute(featureElement, 'display', 'none');
      } else {
        mutation.removeAttribute(featureElement, 'display');
      }
      mutation.setTextContent(
        legend.entry.querySelector('text'),
        legendEntries.value[0]?.caption || ''
      );
      mutation.setAttribute(
        legend.entry,
        'data-legend-key',
        legendEntries.value[0]?.caption || ''
      );
      return true;
    });
  });
  const canonicalState = () => [
    { target: featureVisibilityOverrides },
    { target: legendEntries, key: 'value', deep: true },
    { target: state.deletedLegendEntries, key: 'value', deep: true },
    { target: state.legendColorOverrides },
    { target: state.legendStrokeOverrides },
    { target: state.addedLegendCaptions, key: 'value', deep: true },
    { target: state.clickedFeature, key: 'value', deep: true },
    { target: state.semanticFileWatchersSuppressed, key: 'value' }
  ];
  let ownsMountedTransaction = false;
  const history = createHistoryManager({
    buildIntent: snapshots.buildHistoryIntent,
    applyIntent: (intent, context = {}) => {
      if (context.rollback === true) return snapshots.applyHistoryIntent(intent, context);
      return fixture.runtime.runDomEditTransaction({
        reason: 'history-intent-apply-test',
        canonicalState: canonicalState(),
        invalidateIndexes: ['features', 'legend'],
        action: ({ lease }) => {
          ownsMountedTransaction = Boolean(lease);
          return snapshots.applyHistoryIntent(intent, { ...context, lease });
        }
      }).catch((error) => {
        if (ownsMountedTransaction) {
          Object.defineProperty(error, 'historyRestoredDomains', {
            configurable: true,
            value: Object.freeze(['features', 'editorState'])
          });
        }
        throw error;
      });
    },
    buildCheckpoint: () => ({}),
    applyCheckpoint: async () => {}
  });

  await history.initializeIntentBaseline('History rollback baseline');
  let afterLegendOwner = null;
  await history.runUndoable('Apply editor state', () => {
    featureVisibilityOverrides['feature-a'] = 'off';
    afterLegendOwner = [{ caption: 'After', color: '#222222' }];
    legendEntries.value = afterLegendOwner;
    featureElement.setAttribute('display', 'none');
    legend.entry.setAttribute('data-legend-key', 'After');
    legend.entry.querySelector('text').textContent = 'After';
  });
  const afterEntryOwner = afterLegendOwner[0];
  const beforeDom = serializeFakeSvg(svg);

  await assert.rejects(history.undo(), /serialization failed/);

  assert.equal(featureVisibilityOverrides['feature-a'], 'off');
  assert.equal(legendEntries.value, afterLegendOwner);
  assert.equal(legendEntries.value[0], afterEntryOwner);
  assert.equal(state.clickedFeature.value, clickedOwner);
  assert.equal(serializeFakeSvg(svg), beforeDom);
  assertRuntimeRestored(fixture);
  assert.equal(history.getUndoCount(), 1);
  assert.equal(history.getRedoCount(), 0);
});
