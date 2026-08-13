import assert from 'node:assert/strict';
import test from 'node:test';

import {
  buildFeatureLabelCommand,
  buildFeatureLabelSemanticSnapshots
} from '../../gbdraw/web/js/app/feature-editor/label-command.js';
import {
  planFeatureLabelChange,
  resolveFeatureLabelPlan,
  stableFeatureTargetKey
} from '../../gbdraw/web/js/app/feature-editor/label-scope-plan.js';
import { createHistoryManager } from '../../gbdraw/web/js/services/history.js';
import { styleFingerprint } from '../../gbdraw/web/js/services/style-revision.js';

const clone = (value) => structuredClone(value);
const ref = (value) => ({ value });

class CountingRef {
  constructor(value) { this.current = value; this.writes = 0; }
  get value() { return this.current; }
  set value(value) { this.writes += 1; this.current = value; }
}

const hashes = {
  a: 'fi1_aaaaaaaaaaaaaaaaaaaaaaaaaa',
  b: 'fi1_bbbbbbbbbbbbbbbbbbbbbbbbbb',
  c: 'fi1_cccccccccccccccccccccccccc'
};

const biological = (recordKey, id, instanceHash, product) => ({
  recordKey, biologicalFeatureId: id, instanceHash, type: 'CDS', product,
  qualifiers: { product: [product] }
});

const catalog = {
  schema: 3,
  items: [{
    resultKey: 'result-a', resultIndex: 0, resultName: 'a.svg', recordKeys: ['record-a'],
    biologicalFeatures: [biological('record-a', 'feature-a', hashes.a, 'Kinase')],
    features: [{ recordKey: 'record-a', biologicalFeatureId: 'feature-a', svgId: 'svg-a' }]
  }, {
    resultKey: 'result-b', resultIndex: 1, resultName: 'b.svg', recordKeys: ['record-b'],
    biologicalFeatures: [
      biological('record-b', 'feature-b', hashes.b, 'Kinase'),
      biological('record-b', 'feature-c', hashes.c, 'Transporter')
    ],
    features: [
      { recordKey: 'record-b', biologicalFeatureId: 'feature-b', svgId: 'svg-b' },
      { recordKey: 'record-b', biologicalFeatureId: 'feature-c', svgId: 'svg-c' }
    ]
  }]
};

const feature = (resultIndex, featureIndex) => {
  const item = catalog.items[resultIndex];
  return {
    ...item.biologicalFeatures[featureIndex],
    resultKey: item.resultKey,
    resultIndex,
    resultName: item.resultName,
    renderedFeatureSvgId: item.features[featureIndex].svgId
  };
};

const targetA = stableFeatureTargetKey(feature(0, 0));
const targetC = stableFeatureTargetKey(feature(1, 1));

const fillFingerprint = styleFingerprint({
  rules: [], appliedPaletteName: 'default', appliedPaletteColors: { CDS: '#cccccc' }
});

const createState = ({ featureOverrides = {}, bulkOverrides = {}, sources = {} } = {}) => ({
  manualSpecificRules: [],
  appliedPaletteName: ref('default'),
  appliedPaletteColors: ref({ CDS: '#cccccc' }),
  labelTextFeatureOverrides: clone(featureOverrides),
  labelTextBulkOverrides: clone(bulkOverrides),
  labelTextFeatureOverrideSources: clone(sources),
  results: new CountingRef([
    { name: 'a.svg', content: 'result-a:Shared rendered' },
    { name: 'b.svg', content: 'result-b:Shared rendered,Shared rendered' }
  ]),
  documentEpoch: ref(3),
  resultGenerationKey: ref(7),
  semanticStyleRevision: ref(4),
  semanticStyleFingerprint: ref(fillFingerprint),
  validatedStyleFingerprintByResultKey: ref(Object.freeze({
    'result-a': fillFingerprint,
    'result-b': fillFingerprint
  })),
  featureCatalog: ref(catalog),
  clickedFeature: ref({ id: 'clicked' }),
  editableLabels: ref([{ key: 'label-a', featureId: 'svg-a', text: 'Shared rendered' }]),
  selectedFeatures: ref([{ id: 'selected' }]),
  selectedResultIndex: ref(0)
});

const readyPlan = ({ state, scope = 'source-annotation-label', newText = 'Renamed' }) => {
  const plan = planFeatureLabelChange({
    catalog,
    presentationLabelsBySvgId: {
      'svg-a': { text: 'Shared rendered', sourceText: 'Kinase' },
      'svg-b': { text: 'Shared rendered', sourceText: 'Kinase' },
      'svg-c': { text: 'Shared rendered', sourceText: 'Transporter' }
    },
    labelTextFeatureOverrides: state.labelTextFeatureOverrides,
    labelTextBulkOverrides: state.labelTextBulkOverrides,
    labelTextFeatureOverrideSources: state.labelTextFeatureOverrideSources,
    selectedFeatureKeys: [targetA, targetC],
    exactScopeAvailable: true,
    resultGenerationKey: 7,
    documentEpoch: 3,
    styleRevision: state.semanticStyleRevision.value,
    styleFingerprint: fillFingerprint
  }, {
    targetFeatureKey: targetA,
    targetFeatureKeys: [targetA, targetC],
    newText,
    source: scope === 'selected-features' ? 'selection-toolbar' : 'popup',
    semanticScope: scope,
    originResultKey: 'result-a',
    originResultIndex: 0,
    resultGenerationKey: 7,
    documentEpoch: 3,
    styleRevision: state.semanticStyleRevision.value,
    styleFingerprint: fillFingerprint
  });
  assert.equal(plan.status, 'ready', JSON.stringify(plan.diagnostics));
  return resolveFeatureLabelPlan(plan);
};

const projectionFactory = ({ fail = false, calls = [] } = {}) => async ({
  direction,
  results,
  labelsByTargetKey,
  affectedResultKeys,
  targetFeatureKeysByResult,
  semanticScope,
  sourceText,
  mounted
}) => {
  calls.push({ direction, affectedResultKeys: [...affectedResultKeys], semanticScope, sourceText, mounted });
  if (fail) throw new Error('label projection failed');
  const indexByKey = new Map([['result-a', 0], ['result-b', 1]]);
  const affected = new Set(affectedResultKeys);
  const nextResults = results.map((result, index) => {
    const resultKey = index === 0 ? 'result-a' : 'result-b';
    if (!affected.has(resultKey)) return result;
    const values = targetFeatureKeysByResult.get(resultKey)
      .map((key) => labelsByTargetKey[key].text).join(',');
    return { ...result, content: `${resultKey}:${values}` };
  });
  return {
    nextResults,
    admissionMetadataByResultKey: Object.fromEntries(
      affectedResultKeys.map((key) => [key, { before: null, after: null }])
    ),
    preparedMountedSvg: mounted?.resultKey && affected.has(mounted.resultKey)
      ? { resultKey: mounted.resultKey }
      : null,
    counters: {
      affectedResults: affected.size,
      mountedSerializations: mounted?.resultKey && affected.has(mounted.resultKey) ? 1 : 0,
      detachedPasses: affected.size,
      changedResults: affected.size,
      changedLabels: affected.size
    }
  };
};

test('source-label scope stages every Result and commits one semantic/map transaction', async () => {
  const state = createState({
    featureOverrides: { 'svg-a': 'Old exact' },
    sources: { 'svg-a': 'Kinase' }
  });
  const originalResults = state.results.value;
  const plan = readyPlan({ state });
  const semantic = buildFeatureLabelSemanticSnapshots({
    plan,
    catalog,
    labelTextFeatureOverrides: state.labelTextFeatureOverrides,
    labelTextBulkOverrides: state.labelTextBulkOverrides,
    labelTextFeatureOverrideSources: state.labelTextFeatureOverrideSources
  });
  assert.deepEqual(semantic.affectedResultKeys, ['result-a', 'result-b']);
  assert.equal(semantic.after.featureOverrides['svg-a'], undefined);
  assert.equal(semantic.after.bulkOverrides.Kinase, 'Renamed');

  const calls = [];
  let reflowCalls = 0;
  const command = await buildFeatureLabelCommand({
    plan,
    state,
    catalog,
    prepareResultProjection: projectionFactory({ calls }),
    prepareLabelReflow: async () => { reflowCalls += 1; return null; }
  });
  assert.equal(reflowCalls, 1, 'reflow prerequisite belongs to mutation-free preflight');
  assert.equal(state.results.value, originalResults);
  assert.equal(state.results.writes, 0);
  assert.equal(await command.apply(), true);
  assert.equal(state.labelTextBulkOverrides.Kinase, 'Renamed');
  assert.equal(state.labelTextFeatureOverrides['svg-a'], undefined);
  assert.deepEqual(state.results.value.map((result) => result.content), [
    'result-a:Renamed', 'result-b:Renamed'
  ]);
  assert.equal(state.results.writes, 1);
  assert.equal(state.semanticStyleRevision.value, 5);
  assert.equal(state.semanticStyleFingerprint.value, fillFingerprint);
  assert.deepEqual(calls[0].affectedResultKeys, ['result-a', 'result-b']);
  assert.ok(command.estimateBytes() > 0);
});

test('single and selected scopes persist exact rendered identities only', async () => {
  const single = createState();
  const singleCommand = await buildFeatureLabelCommand({
    plan: readyPlan({ state: single, scope: 'single', newText: 'Only A' }),
    state: single,
    catalog,
    prepareResultProjection: projectionFactory()
  });
  const untouched = single.results.value[1];
  assert.equal(await singleCommand.apply(), true);
  assert.equal(single.labelTextFeatureOverrides['svg-a'], 'Only A');
  assert.equal(single.labelTextFeatureOverrideSources['svg-a'], 'Kinase');
  assert.equal(single.labelTextBulkOverrides.Kinase, undefined);
  assert.equal(single.results.value[1], untouched);

  const selected = createState();
  const selectedCommand = await buildFeatureLabelCommand({
    plan: readyPlan({ state: selected, scope: 'selected-features', newText: 'Chosen' }),
    state: selected,
    catalog,
    prepareResultProjection: projectionFactory()
  });
  assert.equal(await selectedCommand.apply(), true);
  assert.deepEqual(selected.labelTextFeatureOverrides, {
    'svg-a': 'Chosen',
    'svg-c': 'Chosen'
  });
  assert.equal(selected.labelTextFeatureOverrides['svg-b'], undefined);
});

test('preflight and synchronous commit failures preserve maps, Results, ledger, and presentation', async () => {
  const preflight = createState();
  const beforePreflight = preflight.results.value;
  await assert.rejects(() => buildFeatureLabelCommand({
    plan: readyPlan({ state: preflight }),
    state: preflight,
    catalog,
    prepareResultProjection: projectionFactory({ fail: true })
  }), /label projection failed/);
  assert.equal(preflight.results.value, beforePreflight);
  assert.deepEqual(preflight.labelTextBulkOverrides, {});
  assert.equal(preflight.semanticStyleRevision.value, 4);

  const commit = createState();
  const beforeResults = commit.results.value;
  const beforeClicked = commit.clickedFeature.value;
  const beforeLedger = commit.validatedStyleFingerprintByResultKey.value;
  const command = await buildFeatureLabelCommand({
    plan: readyPlan({ state: commit }),
    state: commit,
    catalog,
    prepareResultProjection: projectionFactory(),
    refreshPresentation: () => {
      commit.clickedFeature.value = null;
      throw new Error('label presentation failed');
    }
  });
  await assert.rejects(() => command.apply(), /label presentation failed/);
  assert.equal(commit.results.value, beforeResults);
  assert.equal(commit.clickedFeature.value, beforeClicked);
  assert.equal(commit.validatedStyleFingerprintByResultKey.value, beforeLedger);
  assert.deepEqual(commit.labelTextBulkOverrides, {});
  assert.equal(commit.semanticStyleRevision.value, 4);
});

test('Apply, Undo, and Redo use the same all-Result atomic projection', async () => {
  const state = createState();
  const history = createHistoryManager({
    buildIntent: () => ({
      feature: clone(state.labelTextFeatureOverrides),
      bulk: clone(state.labelTextBulkOverrides),
      results: clone(state.results.value),
      revision: state.semanticStyleRevision.value
    }),
    applyIntent: () => {},
    buildCheckpoint: () => ({}),
    applyCheckpoint: () => {}
  });
  await history.captureBaseline();
  const calls = [];
  const command = await buildFeatureLabelCommand({
    plan: readyPlan({ state }),
    state,
    catalog,
    prepareResultProjection: projectionFactory({ calls })
  });
  assert.equal(await history.runUndoableCommand('Change feature label', () => command), true);
  assert.equal(history.getUndoCount(), 1);
  assert.equal(await history.undo(), true);
  assert.deepEqual(state.results.value.map((result) => result.content), [
    'result-a:Shared rendered', 'result-b:Shared rendered'
  ]);
  assert.deepEqual(state.labelTextBulkOverrides, {});
  assert.equal(await history.redo(), true);
  assert.deepEqual(state.results.value.map((result) => result.content), [
    'result-a:Renamed', 'result-b:Renamed'
  ]);
  assert.deepEqual(calls.map((entry) => entry.direction), ['apply', 'undo', 'redo']);
  assert.equal(state.semanticStyleRevision.value, 7);
});

test('History post-apply failure compensates label maps, Results, and presentation exactly', async () => {
  const state = createState();
  let failCapture = false;
  const history = createHistoryManager({
    buildIntent: () => {
      if (failCapture) {
        failCapture = false;
        state.selectedResultIndex.value = 1;
        throw new Error('label intent capture failed');
      }
      return {
        feature: clone(state.labelTextFeatureOverrides),
        bulk: clone(state.labelTextBulkOverrides),
        results: clone(state.results.value),
        revision: state.semanticStyleRevision.value
      };
    },
    applyIntent: () => {},
    buildCheckpoint: () => ({}),
    applyCheckpoint: () => {}
  });
  await history.captureBaseline();
  const command = await buildFeatureLabelCommand({
    plan: readyPlan({ state }),
    state,
    catalog,
    prepareResultProjection: projectionFactory(),
    refreshPresentation: () => {
      state.clickedFeature.value = null;
      return true;
    }
  });
  const beforeResults = state.results.value;
  const beforeClicked = state.clickedFeature.value;
  const beforeEditable = state.editableLabels.value;
  const beforeLedger = state.validatedStyleFingerprintByResultKey.value;
  failCapture = true;
  await assert.rejects(
    () => history.runUndoableCommand('Change feature label', () => command),
    /label intent capture failed/
  );
  assert.equal(state.results.value, beforeResults);
  assert.equal(state.clickedFeature.value, beforeClicked);
  assert.equal(state.editableLabels.value, beforeEditable);
  assert.equal(state.validatedStyleFingerprintByResultKey.value, beforeLedger);
  assert.equal(state.selectedResultIndex.value, 0);
  assert.deepEqual(state.labelTextFeatureOverrides, {});
  assert.deepEqual(state.labelTextBulkOverrides, {});
  assert.equal(state.semanticStyleRevision.value, 4);
  assert.equal(history.getUndoCount(), 0);
  assert.equal(history.getRedoCount(), 0);
});

test('revision, Result bytes, and ledger replacement invalidate a prepared command', async () => {
  for (const mutate of [
    (state) => { state.semanticStyleRevision.value += 1; },
    (state) => { state.results.value[1].content = 'external'; },
    (state) => {
      state.validatedStyleFingerprintByResultKey.value = {
        'result-a': fillFingerprint, 'result-b': fillFingerprint
      };
    }
  ]) {
    const state = createState();
    const command = await buildFeatureLabelCommand({
      plan: readyPlan({ state }),
      state,
      catalog,
      prepareResultProjection: projectionFactory()
    });
    mutate(state);
    assert.equal(await command.apply(), false);
    assert.deepEqual(state.labelTextBulkOverrides, {});
  }
});

test('identical exact assignment is a no-op with no Result preparation', async () => {
  const state = createState({
    featureOverrides: { 'svg-a': 'Same' },
    sources: { 'svg-a': 'Kinase' }
  });
  let preparations = 0;
  const command = await buildFeatureLabelCommand({
    plan: readyPlan({ state, scope: 'single', newText: 'Same' }),
    state,
    catalog,
    prepareResultProjection: async (input) => {
      preparations += 1;
      return projectionFactory()(input);
    }
  });
  assert.equal(command.noop, true);
  assert.equal(command.estimateBytes(), 0);
  assert.equal(await command.apply(), false);
  assert.equal(preparations, 0);
  assert.equal(state.results.writes, 0);
});
