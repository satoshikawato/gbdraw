import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import test from 'node:test';

import { createFeatureBulkStyleActions } from '../../gbdraw/web/js/app/feature-editor/bulk-style-actions.js';
import { createFeatureRuleActions } from '../../gbdraw/web/js/app/feature-editor/rule-actions.js';
import { FEATURE_INSTANCE_HASH_QUALIFIER } from '../../gbdraw/web/js/services/feature-instance-identity.js';
import {
  FEATURE_SEMANTIC_SCOPE_QUALIFIER,
  encodeFeatureSemanticSelector
} from '../../gbdraw/web/js/app/feature-editor/semantic-fill-selectors.js';
import { createResultsManager } from '../../gbdraw/web/js/app/results.js';
import { setupWatchers } from '../../gbdraw/web/js/app/watchers.js';

const ref = (value) => ({ value });
const clone = (value) => structuredClone(value);

const emptyBulkState = () => ({
  manualSpecificRules: [],
  appliedPaletteName: ref('default'),
  appliedPaletteColors: ref({ CDS: '#cccccc' }),
  featureColorOverrides: {},
  results: ref([]),
  legendEntries: ref([]),
  manualLegendEntries: ref([]),
  validatedStyleFingerprintByResultKey: ref(Object.freeze({})),
  semanticStyleRevision: ref(0),
  semanticStyleFingerprint: ref(''),
  documentEpoch: ref(1),
  resultGenerationKey: ref(0),
  featureCatalog: ref({ schema: 3, items: [] }),
  selectedResultIndex: ref(0),
  svgContainer: ref(null),
  clickedFeature: ref(null),
  selectedFeatures: ref([]),
  adv: { features: [] }
});

test('bulk style facade accepts one complete semantic request before first Generate', async () => {
  const state = emptyBulkState();
  const historyCalls = [];
  const history = {
    async runUndoableCommand(label, build) {
      historyCalls.push(label);
      const command = await build();
      return command.apply();
    }
  };
  const actions = createFeatureBulkStyleActions({ state, history });
  const changed = await actions.requestFeatureBulkStyleChange({
    writerKind: 'rule-add',
    label: 'Add specific color rule',
    replacement: {
      rules: [{
        feat: 'CDS', qual: 'product', val: '^Kinase$', color: '#abcdef', cap: 'Core'
      }]
    }
  });
  assert.equal(changed, true);
  assert.equal(historyCalls.length, 1);
  assert.equal(state.manualSpecificRules[0].color, '#abcdef');
  assert.equal(state.semanticStyleRevision.value, 1);
  assert.match(state.semanticStyleFingerprint.value, /^sf1_/);
  assert.deepEqual(state.validatedStyleFingerprintByResultKey.value, {});

  const noOp = await actions.requestFeatureBulkStyleChange({
    replacement: { rules: clone(state.manualSpecificRules) }
  });
  assert.equal(noOp, false);
  assert.equal(state.semanticStyleRevision.value, 1);
});

test('Result-local legend extraction failure leaves bulk semantic and presentation state untouched', async () => {
  const state = emptyBulkState();
  const result = { name: 'unadmitted.svg', content: '<svg />' };
  const results = [result];
  const legend = [{ caption: 'Selected only', color: '#123456', strokeColor: '#654321' }];
  state.featureCatalog.value = { schema: 3, items: [{
    resultKey: 'result-a', resultIndex: 0, resultName: 'unadmitted.svg',
    biologicalFeatures: [], features: []
  }] };
  state.results.value = results;
  state.legendEntries.value = legend;
  state.semanticStyleFingerprint.value = 'sf1_before';
  state.validatedStyleFingerprintByResultKey.value = Object.freeze({ 'result-a': 'sf1_before' });
  let historyBuilds = 0;
  const actions = createFeatureBulkStyleActions({
    state,
    history: {
      async runUndoableCommand(_label, build) {
        historyBuilds += 1;
        return build();
      }
    },
    onError() {}
  });

  assert.equal(await actions.requestFeatureBulkStyleChange({
    replacement: {
      rules: [{ feat: 'CDS', qual: 'product', val: 'Kinase', color: '#abcdef', cap: 'Core' }]
    }
  }), false);
  assert.equal(historyBuilds, 1);
  assert.equal(state.results.value, results);
  assert.equal(state.results.value[0], result);
  assert.equal(state.legendEntries.value, legend);
  assert.deepEqual(state.manualSpecificRules, []);
  assert.deepEqual(state.featureColorOverrides, {});
  assert.equal(state.semanticStyleRevision.value, 0);
  assert.equal(state.semanticStyleFingerprint.value, 'sf1_before');
  assert.deepEqual(state.validatedStyleFingerprintByResultKey.value, { 'result-a': 'sf1_before' });
});

const ruleState = () => ({
  currentColors: ref({ CDS: '#cccccc' }),
  appliedPaletteColors: ref({ CDS: '#cccccc' }),
  newColorFeat: ref('CDS'),
  newColorVal: ref('#abcdef'),
  manualSpecificRules: [{
    feat: 'CDS', qual: 'product', val: '^Kinase$', color: '#111111', cap: 'Core', fromFile: true
  }],
  newSpecRule: { feat: 'tRNA', qual: 'product', val: '^tRNA', color: '#222222', cap: 'Transfer' },
  specificRulePresets: [{ id: 'bakta', label: 'Bakta', path: '/bakta.tsv' }],
  selectedSpecificPreset: ref('bakta'),
  specificRulePresetLoading: ref(false),
  manualPriorityRules: [],
  newPriorityRule: {},
  adv: { features: [], feature_shapes: {}, legend_box_size: 10, legend_font_size: 10 },
  newFeatureToAdd: ref(''),
  extractedFeatures: ref([]),
  featureColorOverrides: {},
  editableLabels: ref([]),
  labelTextFeatureOverrides: {},
  labelTextBulkOverrides: {},
  fileLegendCaptions: ref(new Set(['Core']))
});

const committedRuleHarness = (state, { firstRuleColor = '#111111' } = {}) => {
  const requests = [];
  return {
    requests,
    actions: {
      async requestFeatureBulkStyleChange(request) {
        requests.push(clone(request));
        assert.equal(
          state.manualSpecificRules[0]?.color,
          requests.length === 1 ? firstRuleColor : state.manualSpecificRules[0]?.color,
          'caller must not mutate rules before command acceptance'
        );
        if (request.replacement.rules) {
          state.manualSpecificRules.splice(
            0,
            state.manualSpecificRules.length,
            ...clone(request.replacement.rules)
          );
        }
        if (request.replacement.appliedPaletteColors) {
          state.appliedPaletteColors.value = clone(request.replacement.appliedPaletteColors);
        }
        if (request.presentationPatch?.legendBoxSize !== undefined) {
          state.adv.legend_box_size = request.presentationPatch.legendBoxSize;
        }
        if (request.presentationPatch?.legendFontSize !== undefined) {
          state.adv.legend_font_size = request.presentationPatch.legendFontSize;
        }
        return true;
      }
    }
  };
};

test('specific rule field/add/reorder/remove/clear handlers only submit complete snapshots', async () => {
  const state = ruleState();
  const bulk = committedRuleHarness(state);
  const actions = createFeatureRuleActions({
    state,
    nextTick: async () => {},
    legendActions: {},
    bulkStyleActions: bulk.actions
  });

  assert.equal(await actions.setSpecificRuleField(0, 'color', '#abcdef'), true);
  assert.equal(bulk.requests.at(-1).writerKind, 'rule-field');
  assert.equal(state.manualSpecificRules[0].fromFile, undefined);
  assert.deepEqual([...state.fileLegendCaptions.value], []);

  assert.equal(await actions.addSpecificRule(), true);
  assert.equal(bulk.requests.at(-1).writerKind, 'rule-add');
  assert.equal(state.newSpecRule.val, '');
  assert.equal(state.manualSpecificRules.length, 2);

  assert.equal(await actions.moveSpecificRuleUp(1), true);
  assert.equal(bulk.requests.at(-1).writerKind, 'rule-reorder');
  assert.equal(state.manualSpecificRules[0].feat, 'tRNA');

  assert.equal(await actions.removeSpecificRule(0), true);
  assert.equal(bulk.requests.at(-1).writerKind, 'rule-remove');
  assert.equal(state.manualSpecificRules.length, 1);

  assert.equal(await actions.clearAllSpecificRules(), true);
  assert.equal(bulk.requests.at(-1).writerKind, 'rule-clear');
  assert.deepEqual(state.manualSpecificRules, []);
});

test('generic rule handlers cannot edit, delete, clear, or reorder managed Feature rows', async () => {
  const state = ruleState();
  const exact = {
    feat: 'CDS', qual: FEATURE_INSTANCE_HASH_QUALIFIER,
    val: 'fi1_aaaaaaaaaaaaaaaaaaaaaaaaaa', color: '#abcdef', cap: 'Exact', match: 'literal'
  };
  const semantic = {
    feat: 'CDS', qual: FEATURE_SEMANTIC_SCOPE_QUALIFIER,
    val: encodeFeatureSemanticSelector('feature-type', 'CDS'),
    color: '#123456', cap: 'All CDS', match: 'literal'
  };
  const legacyReservedQualifier = {
    feat: 'CDS', qual: FEATURE_INSTANCE_HASH_QUALIFIER,
    val: '^legacy$', color: '#654321', cap: 'Biological', match: 'regex'
  };
  const legacySemanticQualifier = {
    feat: 'CDS', qual: FEATURE_SEMANTIC_SCOPE_QUALIFIER,
    val: '^bio.*$', color: '#456789', cap: 'Biological semantic spelling', match: 'regex'
  };
  state.manualSpecificRules.splice(
    0,
    state.manualSpecificRules.length,
    exact,
    semantic,
    legacyReservedQualifier,
    legacySemanticQualifier
  );
  const bulk = committedRuleHarness(state, { firstRuleColor: '#abcdef' });
  const actions = createFeatureRuleActions({
    state,
    nextTick: async () => {},
    legendActions: {},
    bulkStyleActions: bulk.actions
  });

  assert.equal(actions.isOpaqueSpecificRule(exact), true);
  assert.equal(actions.isOpaqueSpecificRule(semantic), true);
  assert.equal(actions.isOpaqueSpecificRule(legacyReservedQualifier), false);
  assert.equal(actions.isOpaqueSpecificRule(legacySemanticQualifier), false);
  assert.equal(await actions.setSpecificRuleField(0, 'color', '#000000'), false);
  assert.equal(await actions.removeSpecificRule(1), false);
  assert.equal(await actions.moveSpecificRuleDown(0), false);
  assert.equal(await actions.moveSpecificRuleUp(2), false);
  assert.equal(bulk.requests.length, 0);

  assert.equal(await actions.setSpecificRuleField(2, 'color', '#fedcba'), true);
  assert.equal(state.manualSpecificRules[2].color, '#fedcba');
  assert.equal(await actions.setSpecificRuleField(3, 'color', '#abcdef'), true);
  assert.equal(state.manualSpecificRules[3].color, '#abcdef');
  assert.equal(await actions.clearAllSpecificRules(), true);
  assert.deepEqual(state.manualSpecificRules, [exact, semantic]);
});

test('preset parsing preflights rules and default colors in one bulk request', async () => {
  const state = ruleState();
  const bulk = committedRuleHarness(state);
  const actions = createFeatureRuleActions({
    state,
    nextTick: async () => {},
    legendActions: {},
    bulkStyleActions: bulk.actions
  });
  const originalFetch = globalThis.fetch;
  globalThis.fetch = async () => ({
    ok: true,
    text: async () => 'CDS\tproduct\t^Enzyme$\t#123456\tEnzyme\n'
  });
  try {
    assert.equal(await actions.applySpecificRulePreset(), true);
  } finally {
    globalThis.fetch = originalFetch;
  }
  const request = bulk.requests.at(-1);
  assert.equal(request.writerKind, 'preset');
  assert.equal(request.replacement.rules[0].fromFile, true);
  assert.equal(request.replacement.appliedPaletteColors.CDS, '#cccccc');
  assert.deepEqual(request.presentationPatch, { legendBoxSize: 12, legendFontSize: 12 });
  assert.equal(state.manualSpecificRules[0].cap, 'Enzyme');
  assert.deepEqual([...state.fileLegendCaptions.value], ['Enzyme']);
  assert.equal(state.adv.legend_box_size, 12);
  assert.equal(state.adv.legend_font_size, 12);
});

test('rule mutation fails explicitly when the bulk owner is not injected', async () => {
  const state = ruleState();
  const actions = createFeatureRuleActions({ state, nextTick: async () => {}, legendActions: {} });
  await assert.rejects(
    actions.setSpecificRuleField(0, 'color', '#abcdef'),
    /Bulk Feature style actions are unavailable/
  );
  assert.equal(state.manualSpecificRules[0].color, '#111111');
});

const paletteState = () => ({
  pyodideReady: ref(false),
  svgContent: ref(''),
  mode: ref('linear'),
  shouldDeferCircularPreviewUpdates: ref(false),
  svgContainer: ref(null),
  cInputType: ref('gb'),
  files: {},
  linearSeqs: [],
  form: {},
  adv: {},
  selectedResultIndex: ref(0),
  results: ref([]),
  skipCaptureBaseConfig: ref(false),
  paletteDefinitions: ref({
    default: { CDS: '#cccccc' },
    ocean: { CDS: '#116699', tRNA: '#44aacc' }
  }),
  selectedPalette: ref('default'),
  currentColors: ref({ CDS: '#cccccc' }),
  paletteInstantPreviewEnabled: ref(true),
  appliedPaletteName: ref('default'),
  appliedPaletteColors: ref({ CDS: '#cccccc' }),
  pendingPaletteName: ref(''),
  pendingPaletteColors: ref({}),
  normalizePaletteColors: (colors) => ({ ...colors }),
  normalizePaletteDefinitions: (palettes) => ({ ...palettes })
});

test('palette acceptance uses bulk snapshots while disabled instant preview stays a draft', async () => {
  const state = paletteState();
  const requests = [];
  const bulkStyleActions = {
    async requestFeatureBulkStyleChange(request) {
      requests.push(clone(request));
      state.appliedPaletteName.value = request.replacement.appliedPaletteName;
      state.appliedPaletteColors.value = clone(request.replacement.appliedPaletteColors);
      state.selectedPalette.value = request.replacement.appliedPaletteName;
      state.currentColors.value = clone(request.replacement.appliedPaletteColors);
      state.pendingPaletteName.value = '';
      state.pendingPaletteColors.value = {};
      return true;
    }
  };
  const manager = createResultsManager({
    state,
    getPyodide: () => null,
    legendLayout: { refreshCompositionGeometry() {} },
    bulkStyleActions
  });

  state.selectedPalette.value = 'ocean';
  assert.equal(await manager.updatePalette(), true);
  assert.equal(requests.at(-1).writerKind, 'palette-default-acceptance');
  assert.deepEqual(requests.at(-1).replacement.appliedPaletteColors, {
    CDS: '#116699', tRNA: '#44aacc'
  });
  const appliedRequestCount = requests.length;
  assert.equal(await manager.syncPaletteDraftState(), false);
  assert.equal(requests.length, appliedRequestCount, 'commit-side draft sync must be a no-op');

  state.pendingPaletteName.value = 'ocean';
  state.pendingPaletteColors.value = clone(state.appliedPaletteColors.value);
  assert.equal(await manager.applyPaletteDraftToPreview(), false);
  assert.equal(requests.length, appliedRequestCount);
  assert.equal(state.pendingPaletteName.value, '');
  assert.deepEqual(state.pendingPaletteColors.value, {});

  state.paletteInstantPreviewEnabled.value = false;
  state.selectedPalette.value = 'default';
  const requestCount = requests.length;
  assert.equal(await manager.updatePalette(), false);
  assert.equal(requests.length, requestCount);
  assert.equal(state.pendingPaletteName.value, 'default');

  assert.equal(await manager.applyPaletteDraftToPreview(), true);
  assert.equal(requests.length, requestCount + 1);
  assert.equal(state.pendingPaletteName.value, '');

  state.currentColors.value = { CDS: '#eeeeee' };
  assert.equal(await manager.syncPaletteDraftState(), true);
  assert.equal(requests.at(-1).label, 'Change default feature color');
});

test('watcher source has no rule/palette mutation projection owner and imports call bulk actions', async () => {
  const source = await readFile(
    new URL('../../gbdraw/web/js/app/watchers.js', import.meta.url),
    'utf8'
  );
  assert.doesNotMatch(source, /\(\) => \[\.\.\.manualSpecificRules\]/);
  assert.doesNotMatch(source, /manualSpecificRules\.splice\(0, manualSpecificRules\.length/);
  assert.doesNotMatch(source, /Object\.entries\(colors\)\.forEach/);
  assert.doesNotMatch(source, /prepared\.fileLegendCaptions[\s\S]{0,300}extractLegendEntries\(\)/);
  assert.match(source, /writerKind: 'specific-color-tsv-import'/);
  assert.match(source, /writerKind: 'default-color-tsv-import'/);
  const setupSource = await readFile(
    new URL('../../gbdraw/web/js/app/app-setup.js', import.meta.url),
    'utf8'
  );
  const afterApplyOwner = setupSource.match(
    /setAfterApplyHistoryIntent\([\s\S]*?\n\s*}\);/
  )?.[0] || '';
  assert.doesNotMatch(afterApplyOwner, /applyPaletteToSvg|applySpecificRulesToSvg/);
});

const fileLike = (content) => {
  const bytes = new TextEncoder().encode(content);
  return { arrayBuffer: async () => bytes.slice().buffer };
};

const watcherHarness = () => {
  const required = {
    manualSpecificRules: [
      { feat: 'CDS', qual: 'product', val: '^Manual$', color: '#111111', cap: 'Manual' },
      { feat: 'CDS', qual: 'product', val: '^Old$', color: '#222222', cap: 'Old', fromFile: true }
    ],
    files: {},
    currentColors: ref({ CDS: '#cccccc' }),
    appliedPaletteColors: ref({ CDS: '#cccccc' }),
    semanticFileWatchersSuppressed: ref(false),
    sessionImportRollbackInProgress: ref(false),
    fileLegendCaptions: ref(new Set(['Old'])),
    addedLegendCaptions: ref(new Set(['Old'])),
    manualPriorityRules: [],
    manualWhitelist: [],
    manualBlacklist: ref(''),
    linearSeqs: [],
    form: {},
    adv: {},
    layoutPreferences: { circular: { single: {}, multi: {} } },
    featureVisibilityManualRules: [],
    featureVisibilityOverrides: {},
    featureVisibilitySelectorCache: {},
    labelTextBulkOverrides: {},
    labelTextFeatureOverrides: {},
    labelTextFeatureOverrideSources: {},
    labelVisibilityOverrides: {},
    orthogroupNameOverrides: {},
    orthogroupDescriptionOverrides: {},
    labelTextScopeDialog: {},
    globalLabelModeDialog: {}
  };
  const state = new Proxy(required, {
    get(target, key) {
      if (!(key in target)) target[key] = ref(null);
      return target[key];
    }
  });
  const registrations = [];
  const requests = [];
  let legendExtractions = 0;
  setupWatchers({
    state,
    watch(source, callback) {
      registrations.push({ source, callback });
    },
    nextTick() {},
    onMounted() {},
    legendActions: {
      extractLegendEntries() { legendExtractions += 1; },
      setupLegendDrag() {},
      refreshLegendDragAffordances() {}
    },
    featureActions: {
      attachSvgFeatureHandlers() {}, syncLabelEditor() {}
    },
    legendLayout: {
      applyCanvasPadding() {}, captureBaseConfig() {}, captureOriginalStroke() {},
      repositionForLegendChange() {}, refreshDiagramDragAffordances() {}, setupDiagramDrag() {}
    },
    resultsManager: {
      applyPaletteDraftToPreview() {}, clearPendingPaletteDraft() {},
      scheduleDefinitionUpdate() {}, cancelDefinitionUpdate() {}, syncPaletteDraftState() {}
    },
    resetRightDrawer() {},
    bulkStyleActions: {
      async requestFeatureBulkStyleChange(request) {
        requests.push(clone(request));
        if (request.replacement.rules) {
          state.manualSpecificRules.splice(
            0,
            state.manualSpecificRules.length,
            ...clone(request.replacement.rules)
          );
        }
        if (request.replacement.appliedPaletteColors) {
          state.appliedPaletteColors.value = clone(request.replacement.appliedPaletteColors);
          state.currentColors.value = clone(request.replacement.appliedPaletteColors);
          state.pendingPaletteName.value = '';
          state.pendingPaletteColors.value = {};
        }
        return true;
      }
    }
  });
  const callbackFor = (needle) => registrations.find(({ source }) => (
    typeof source === 'function' && source.toString().includes(needle)
  ))?.callback;
  return { state, requests, callbackFor, legendExtractions: () => legendExtractions };
};

test('specific/default TSV watcher acceptance submits one immutable bulk snapshot', async () => {
  const harness = watcherHarness();
  const beforeRules = clone(harness.state.manualSpecificRules);
  await harness.callbackFor('files.t_color')(
    fileLike('CDS\tproduct\t^Imported$\t#abcdef\tImported\n')
  );
  assert.equal(harness.requests.length, 1);
  assert.equal(harness.requests[0].writerKind, 'specific-color-tsv-import');
  assert.deepEqual(beforeRules, [
    { feat: 'CDS', qual: 'product', val: '^Manual$', color: '#111111', cap: 'Manual' },
    { feat: 'CDS', qual: 'product', val: '^Old$', color: '#222222', cap: 'Old', fromFile: true }
  ]);
  assert.deepEqual(harness.state.manualSpecificRules.map((rule) => rule.cap), ['Manual', 'Imported']);
  assert.deepEqual([...harness.state.fileLegendCaptions.value], ['Imported']);
  assert.equal(harness.legendExtractions(), 0, 'bulk command owns the only legend reconciliation');

  await harness.callbackFor('files.d_color')(fileLike('tRNA\t#123456\n'));
  assert.equal(harness.requests.length, 2);
  assert.equal(harness.requests[1].writerKind, 'default-color-tsv-import');
  assert.deepEqual(harness.requests[1].replacement.appliedPaletteColors, {
    CDS: '#cccccc', tRNA: '#123456'
  });
  assert.deepEqual(harness.state.currentColors.value, { CDS: '#cccccc', tRNA: '#123456' });
});
