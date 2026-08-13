import assert from 'node:assert/strict';
import test from 'node:test';

import {
  buildFeatureStrokeCommand,
  buildFeatureStrokeSemanticSnapshots,
  deriveFeatureStrokeTargets
} from '../../gbdraw/web/js/app/feature-editor/stroke-command.js';
import {
  planFeatureStrokeChange,
  resolveFeatureStrokePlan,
  stableFeatureTargetKey
} from '../../gbdraw/web/js/app/feature-editor/stroke-scope-plan.js';
import { createHistoryManager } from '../../gbdraw/web/js/services/history.js';
import { styleFingerprint } from '../../gbdraw/web/js/services/style-revision.js';

const clone = (value) => structuredClone(value);
const ref = (value) => ({ value });

class CountingRef {
  constructor(value) {
    this.current = value;
    this.writes = 0;
  }

  get value() { return this.current; }
  set value(value) {
    this.writes += 1;
    this.current = value;
  }
}

const hashes = {
  a: 'fi1_aaaaaaaaaaaaaaaaaaaaaaaaaa',
  b: 'fi1_bbbbbbbbbbbbbbbbbbbbbbbbbb',
  c: 'fi1_cccccccccccccccccccccccccc'
};

const biological = (recordKey, id, instanceHash, legendCaption) => ({
  recordKey,
  biologicalFeatureId: id,
  instanceHash,
  type: 'CDS',
  legendCaption
});

const rendered = (recordKey, id, svgId) => ({
  recordKey,
  biologicalFeatureId: id,
  svgId
});

const catalog = {
  schema: 3,
  items: [{
    resultKey: 'result-a',
    resultIndex: 0,
    resultName: 'a.svg',
    recordKeys: ['record-a'],
    biologicalFeatures: [biological('record-a', 'feature-a', hashes.a, 'Core')],
    features: [rendered('record-a', 'feature-a', 'svg-a')]
  }, {
    resultKey: 'result-b',
    resultIndex: 1,
    resultName: 'b.svg',
    recordKeys: ['record-b'],
    biologicalFeatures: [
      biological('record-b', 'feature-b', hashes.b, 'Core'),
      biological('record-b', 'feature-c', hashes.c, 'Accessory')
    ],
    features: [
      rendered('record-b', 'feature-b', 'svg-b'),
      rendered('record-b', 'feature-c', 'svg-c')
    ]
  }]
};

const catalogFeature = (resultIndex, featureIndex) => {
  const item = catalog.items[resultIndex];
  const feature = item.biologicalFeatures[featureIndex];
  return {
    ...feature,
    resultKey: item.resultKey,
    resultIndex,
    resultName: item.resultName,
    renderedFeatureSvgId: item.features[featureIndex].svgId
  };
};

const targetA = stableFeatureTargetKey(catalogFeature(0, 0));
const targetB = stableFeatureTargetKey(catalogFeature(1, 0));
const targetC = stableFeatureTargetKey(catalogFeature(1, 1));
const pairA = 'record-a\u0000feature-a';
const pairB = 'record-b\u0000feature-b';
const pairC = 'record-b\u0000feature-c';

const fillFingerprint = styleFingerprint({
  rules: [],
  appliedPaletteName: 'default',
  appliedPaletteColors: { CDS: '#cccccc' }
});

const createState = ({ featureOverrides = {}, legendOverrides = {} } = {}) => ({
  manualSpecificRules: [],
  appliedPaletteName: ref('default'),
  appliedPaletteColors: ref({ CDS: '#cccccc' }),
  featureStrokeOverrides: clone(featureOverrides),
  legendStrokeOverrides: clone(legendOverrides),
  originalSvgStroke: ref({ color: '#000000', width: 0.5 }),
  results: new CountingRef([
    { name: 'a.svg', content: 'a:#000000/0.5' },
    { name: 'b.svg', content: 'b:#000000/0.5' }
  ]),
  documentEpoch: ref(3),
  resultGenerationKey: ref(7),
  semanticStyleRevision: ref(4),
  semanticStyleFingerprint: ref(fillFingerprint),
  clickedFeature: ref({ id: 'clicked' }),
  selectedFeatures: ref([{ id: 'selected' }]),
  legendEntries: ref([{ caption: 'Core', featureIds: ['svg-a'] }]),
  selectedResultIndex: ref(0)
});

const readyPlan = ({
  state,
  scope = 'legend-caption',
  value = { kind: 'color', color: '#abcdef', strokeWidth: 2 },
  selected = [targetA, targetC]
}) => {
  const target = catalogFeature(0, 0);
  const planned = planFeatureStrokeChange({
    catalog,
    featureStrokeOverrides: state.featureStrokeOverrides,
    legendStrokeOverrides: state.legendStrokeOverrides,
    originalSvgStroke: state.originalSvgStroke.value,
    selectedFeatureKeys: selected,
    resultGenerationKey: state.resultGenerationKey.value,
    documentEpoch: state.documentEpoch.value,
    styleRevision: state.semanticStyleRevision.value,
    styleFingerprint: fillFingerprint,
    exactScopeAvailable: true
  }, {
    targetFeatureKey: stableFeatureTargetKey(target),
    targetFeatureKeys: selected,
    value,
    source: scope === 'selected-features' ? 'selection-toolbar' : 'popup',
    semanticScope: scope,
    originResultKey: 'result-a',
    originResultIndex: 0,
    resultGenerationKey: state.resultGenerationKey.value,
    documentEpoch: state.documentEpoch.value,
    styleRevision: state.semanticStyleRevision.value,
    styleFingerprint: fillFingerprint
  });
  assert.equal(planned.status, 'ready', JSON.stringify(planned.diagnostics));
  return resolveFeatureStrokePlan(planned);
};

const projectionFactory = ({ fail = false, calls = [] } = {}) => async ({
  direction,
  results,
  strokesByTargetKey,
  affectedResultKeys,
  targetFeatureKeysByResult,
  legendCaption,
  legendStroke,
  mounted
}) => {
  calls.push({ direction, affectedResultKeys: [...affectedResultKeys], mounted, legendCaption, legendStroke });
  if (fail) throw new Error('stroke projection failed');
  const resultIndexByKey = new Map([['result-a', 0], ['result-b', 1]]);
  const affected = new Set(affectedResultKeys);
  const nextResults = results.map((result, index) => {
    const resultKey = index === 0 ? 'result-a' : 'result-b';
    if (!affected.has(resultKey)) return result;
    const keys = targetFeatureKeysByResult.get(resultKey);
    const signature = keys.map((key) => {
      const stroke = strokesByTargetKey[key];
      return `${stroke.strokeColor}/${stroke.strokeWidth}`;
    }).join(',');
    return { ...result, content: `${resultKey}:${signature}` };
  });
  const mountedAffected = Boolean(mounted?.resultKey && affected.has(mounted.resultKey));
  return {
    nextResults,
    admissionMetadataByResultKey: Object.fromEntries(
      affectedResultKeys.map((key) => [key, { before: null, after: null }])
    ),
    preparedMountedSvg: mountedAffected ? { resultKey: mounted.resultKey } : null,
    counters: {
      affectedResults: affected.size,
      mountedSerializations: mountedAffected ? 1 : 0,
      detachedPasses: affected.size - (mountedAffected ? 1 : 0),
      changedResults: affected.size
    }
  };
};

test('legend-caption batch applies all Results through one map commit and one Result swap', async () => {
  const state = createState({
    featureOverrides: { [pairA]: { strokeColor: '#111111', strokeWidth: 1 } }
  });
  const originalResults = state.results.value;
  const calls = [];
  const plan = readyPlan({ state });
  const semantic = buildFeatureStrokeSemanticSnapshots({
    plan,
    catalog,
    featureStrokeOverrides: state.featureStrokeOverrides,
    legendStrokeOverrides: state.legendStrokeOverrides,
    originalSvgStroke: state.originalSvgStroke.value
  });
  assert.deepEqual(semantic.affectedResultKeys, ['result-a', 'result-b']);
  assert.equal(semantic.after.featureStrokeOverrides[pairA], undefined);
  assert.deepEqual(semantic.after.legendStrokeOverrides.Core, {
    strokeColor: '#abcdef', strokeWidth: 2
  });

  const command = await buildFeatureStrokeCommand({
    plan,
    state,
    catalog,
    prepareResultProjection: projectionFactory({ calls })
  });
  assert.equal(state.results.value, originalResults, 'preflight must not publish Results');
  assert.equal(state.results.writes, 0);
  assert.equal(await command.apply(), true);
  assert.deepEqual(state.legendStrokeOverrides.Core, { strokeColor: '#abcdef', strokeWidth: 2 });
  assert.equal(state.featureStrokeOverrides[pairA], undefined);
  assert.deepEqual(state.results.value.map((result) => result.content), [
    'result-a:#abcdef/2', 'result-b:#abcdef/2'
  ]);
  assert.equal(state.results.writes, 1);
  assert.equal(state.semanticStyleRevision.value, 5);
  assert.equal(state.semanticStyleFingerprint.value, fillFingerprint);
  assert.deepEqual(calls[0].affectedResultKeys, ['result-a', 'result-b']);
  assert.equal(command.counters.projectedResults, 2);
  assert.equal(command.counters.detachedPasses, 2);
  assert.equal(command.counters.resultArraySwaps, 1);
  assert.ok(command.estimateBytes() > 0);
});

test('single and selected scopes use exact feature maps without adding caption ownership', async () => {
  const singleState = createState();
  const singlePlan = readyPlan({ state: singleState, scope: 'single' });
  const singleCommand = await buildFeatureStrokeCommand({
    plan: singlePlan,
    state: singleState,
    catalog,
    prepareResultProjection: projectionFactory()
  });
  const untouchedB = singleState.results.value[1];
  assert.equal(await singleCommand.apply(), true);
  assert.deepEqual(singleState.featureStrokeOverrides[pairA], {
    strokeColor: '#abcdef', strokeWidth: 2
  });
  assert.deepEqual(singleState.legendStrokeOverrides, {});
  assert.equal(singleState.results.value[1], untouchedB);

  const selectedState = createState();
  const selectedPlan = readyPlan({ state: selectedState, scope: 'selected-features' });
  const selectedCommand = await buildFeatureStrokeCommand({
    plan: selectedPlan,
    state: selectedState,
    catalog,
    prepareResultProjection: projectionFactory()
  });
  assert.equal(await selectedCommand.apply(), true);
  assert.deepEqual(Object.keys(selectedState.featureStrokeOverrides).sort(), [pairA, pairC].sort());
  assert.equal(selectedState.featureStrokeOverrides[pairB], undefined);
  assert.deepEqual(selectedState.results.value.map((result) => result.content), [
    'result-a:#abcdef/2', 'result-b:#abcdef/2'
  ]);
});

test('inherit removes exact ownership and projects the caption fallback', async () => {
  const state = createState({
    featureOverrides: { [pairA]: { strokeColor: '#111111', strokeWidth: 3 } },
    legendOverrides: { Core: { strokeColor: '#445566', strokeWidth: 1.5 } }
  });
  const plan = readyPlan({ state, scope: 'single', value: { kind: 'inherit' } });
  const effective = deriveFeatureStrokeTargets({
    plan,
    catalog,
    featureStrokeOverrides: {},
    legendStrokeOverrides: state.legendStrokeOverrides,
    originalSvgStroke: state.originalSvgStroke.value
  });
  assert.deepEqual(effective[targetA], {
    strokeColor: '#445566', strokeWidth: 1.5, caption: 'Core'
  });
  const command = await buildFeatureStrokeCommand({
    plan,
    state,
    catalog,
    prepareResultProjection: projectionFactory()
  });
  assert.equal(await command.apply(), true);
  assert.equal(state.featureStrokeOverrides[pairA], undefined);
  assert.equal(state.results.value[0].content, 'result-a:#445566/1.5');
});

test('projection failure and synchronous commit failure restore exact references and revision', async () => {
  const preflightState = createState();
  const preflightReferences = {
    results: preflightState.results.value,
    clicked: preflightState.clickedFeature.value
  };
  await assert.rejects(() => buildFeatureStrokeCommand({
    plan: readyPlan({ state: preflightState }),
    state: preflightState,
    catalog,
    prepareResultProjection: projectionFactory({ fail: true })
  }), /stroke projection failed/);
  assert.equal(preflightState.results.value, preflightReferences.results);
  assert.equal(preflightState.clickedFeature.value, preflightReferences.clicked);
  assert.deepEqual(preflightState.featureStrokeOverrides, {});
  assert.deepEqual(preflightState.legendStrokeOverrides, {});
  assert.equal(preflightState.semanticStyleRevision.value, 4);

  const commitState = createState();
  const beforeResults = commitState.results.value;
  const beforeClicked = commitState.clickedFeature.value;
  const command = await buildFeatureStrokeCommand({
    plan: readyPlan({ state: commitState }),
    state: commitState,
    catalog,
    prepareResultProjection: projectionFactory(),
    refreshPresentation: () => {
      commitState.clickedFeature.value = null;
      throw new Error('presentation failed');
    }
  });
  await assert.rejects(() => command.apply(), /presentation failed/);
  assert.equal(commitState.results.value, beforeResults);
  assert.equal(commitState.clickedFeature.value, beforeClicked);
  assert.deepEqual(commitState.featureStrokeOverrides, {});
  assert.deepEqual(commitState.legendStrokeOverrides, {});
  assert.equal(commitState.semanticStyleRevision.value, 4);
});

test('apply, Undo, and Redo use the same atomic batch projection', async () => {
  const state = createState();
  const history = createHistoryManager({
    buildIntent: () => ({
      featureStrokeOverrides: clone(state.featureStrokeOverrides),
      legendStrokeOverrides: clone(state.legendStrokeOverrides),
      results: clone(state.results.value),
      revision: state.semanticStyleRevision.value
    }),
    applyIntent: () => {},
    buildCheckpoint: () => ({}),
    applyCheckpoint: () => {}
  });
  await history.captureBaseline();
  const calls = [];
  const command = await buildFeatureStrokeCommand({
    plan: readyPlan({ state }),
    state,
    catalog,
    prepareResultProjection: projectionFactory({ calls })
  });
  assert.equal(await history.runUndoableCommand('Change feature stroke', () => command), true);
  assert.equal(history.getUndoCount(), 1);
  assert.equal(state.semanticStyleRevision.value, 5);
  assert.equal(await history.undo(), true);
  assert.deepEqual(state.results.value.map((result) => result.content), [
    'result-a:#000000/0.5', 'result-b:#000000/0.5'
  ]);
  assert.deepEqual(state.legendStrokeOverrides, {});
  assert.equal(state.semanticStyleRevision.value, 6);
  assert.equal(await history.redo(), true);
  assert.deepEqual(state.results.value.map((result) => result.content), [
    'result-a:#abcdef/2', 'result-b:#abcdef/2'
  ]);
  assert.equal(state.semanticStyleRevision.value, 7);
  assert.deepEqual(calls.map((call) => call.direction), ['apply', 'undo', 'redo']);
});

test('History post-apply failure compensates maps, Results, presentation, and revision exactly', async () => {
  const state = createState();
  let failCapture = false;
  const history = createHistoryManager({
    buildIntent: () => {
      if (failCapture) {
        failCapture = false;
        state.selectedResultIndex.value = 1;
        throw new Error('intent capture failed');
      }
      return {
        featureStrokeOverrides: clone(state.featureStrokeOverrides),
        legendStrokeOverrides: clone(state.legendStrokeOverrides),
        results: clone(state.results.value),
        revision: state.semanticStyleRevision.value
      };
    },
    applyIntent: () => {},
    buildCheckpoint: () => ({}),
    applyCheckpoint: () => {}
  });
  await history.captureBaseline();
  const command = await buildFeatureStrokeCommand({
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
  const beforeSelected = state.selectedFeatures.value;
  const beforeLegend = state.legendEntries.value;
  failCapture = true;
  await assert.rejects(
    () => history.runUndoableCommand('Change feature stroke', () => command),
    /intent capture failed/
  );
  assert.equal(state.results.value, beforeResults);
  assert.equal(state.clickedFeature.value, beforeClicked);
  assert.equal(state.selectedFeatures.value, beforeSelected);
  assert.equal(state.legendEntries.value, beforeLegend);
  assert.equal(state.selectedResultIndex.value, 0);
  assert.deepEqual(state.featureStrokeOverrides, {});
  assert.deepEqual(state.legendStrokeOverrides, {});
  assert.equal(state.semanticStyleRevision.value, 4);
  assert.equal(history.getUndoCount(), 0);
  assert.equal(history.getRedoCount(), 0);
});

test('epoch, generation, revision, Result bytes, and mount ownership are stale guards', async () => {
  for (const [field, nextValue] of [
    ['documentEpoch', 4],
    ['resultGenerationKey', 8],
    ['semanticStyleRevision', 5]
  ]) {
    const state = createState();
    const command = await buildFeatureStrokeCommand({
      plan: readyPlan({ state }),
      state,
      catalog,
      prepareResultProjection: projectionFactory()
    });
    state[field].value = nextValue;
    assert.equal(await command.apply(), false, `${field} must invalidate the command`);
    assert.deepEqual(state.legendStrokeOverrides, {});
  }

  const resultState = createState();
  const resultCommand = await buildFeatureStrokeCommand({
    plan: readyPlan({ state: resultState }),
    state: resultState,
    catalog,
    prepareResultProjection: projectionFactory()
  });
  resultState.results.value[1].content = 'externally changed';
  assert.equal(await resultCommand.apply(), false);
  assert.deepEqual(resultState.legendStrokeOverrides, {});

  const mountState = createState();
  const mountedA = { resultIndex: 0, resultKey: 'result-a', svg: { id: 'mounted-a' } };
  const mountedB = { resultIndex: 1, resultKey: 'result-b', svg: { id: 'mounted-b' } };
  let mounted = mountedA;
  let release;
  let started;
  const waiting = new Promise((resolve) => { release = resolve; });
  const entered = new Promise((resolve) => { started = resolve; });
  const building = buildFeatureStrokeCommand({
    plan: readyPlan({ state: mountState }),
    state: mountState,
    catalog,
    getMountedContext: () => mounted,
    commitMountedProjection: () => true,
    restoreMountedProjection: () => true,
    prepareResultProjection: async (input) => {
      started();
      await waiting;
      return projectionFactory()(input);
    }
  });
  await entered;
  mounted = mountedB;
  release();
  await assert.rejects(() => building, /became stale/);
  assert.deepEqual(mountState.legendStrokeOverrides, {});
  assert.equal(mountState.results.writes, 0);
});

test('mounted commit is compensated when a later synchronous seam fails', async () => {
  const state = createState();
  let liveStroke = 'before';
  const mounted = { resultIndex: 0, resultKey: 'result-a', svg: { id: 'mounted-a' } };
  const command = await buildFeatureStrokeCommand({
    plan: readyPlan({ state }),
    state,
    catalog,
    getMountedContext: () => mounted,
    prepareResultProjection: projectionFactory(),
    captureRuntimeState: () => liveStroke,
    commitMountedProjection: () => {
      liveStroke = 'after';
      return true;
    },
    restoreMountedProjection: ({ snapshot }) => {
      liveStroke = snapshot.runtime;
      return true;
    },
    reconcile: () => { throw new Error('reconcile failed'); }
  });
  assert.equal(command.counters.mountedSerializations, 1);
  assert.equal(command.counters.detachedPasses, 1);
  const beforeResults = state.results.value;
  await assert.rejects(() => command.apply(), /reconcile failed/);
  assert.equal(liveStroke, 'before');
  assert.equal(state.results.value, beforeResults);
  assert.deepEqual(state.legendStrokeOverrides, {});
  assert.equal(state.semanticStyleRevision.value, 4);
});

test('an identical exact assignment is a no-op with no projection or History work', async () => {
  const state = createState({
    featureOverrides: { [pairA]: { strokeColor: '#abcdef', strokeWidth: 2 } }
  });
  let preparations = 0;
  const command = await buildFeatureStrokeCommand({
    plan: readyPlan({ state, scope: 'single' }),
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
