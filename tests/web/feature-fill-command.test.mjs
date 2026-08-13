import assert from 'node:assert/strict';
import test from 'node:test';

import {
  buildFeatureFillCommand,
  buildFeatureFillSemanticSnapshots,
  deriveFeatureFillOverrides
} from '../../gbdraw/web/js/app/feature-editor/fill-command.js';
import {
  planFeatureFillChange,
  resolveFeatureFillPlan,
  stableFeatureTargetKey
} from '../../gbdraw/web/js/app/feature-editor/fill-scope-plan.js';
import { styleFingerprint } from '../../gbdraw/web/js/services/style-revision.js';
import { createHistoryManager } from '../../gbdraw/web/js/services/history.js';
import { exactJsonByteLength } from '../../gbdraw/web/js/app/feature-editor/style-snapshot-command.js';

const ref = (value) => ({ value });
const clone = (value) => structuredClone(value);

class PrototypeRef {
  constructor(value) {
    this.current = value;
  }

  get value() {
    return this.current;
  }

  set value(value) {
    this.current = value;
  }
}

class ProxyingArrayRef extends PrototypeRef {
  get value() {
    return this.current;
  }

  set value(value) {
    this.current = Array.isArray(value) ? new Proxy(value, {}) : value;
  }
}

const hashes = {
  a: 'fi1_aaaaaaaaaaaaaaaaaaaaaaaaaa',
  b: 'fi1_bbbbbbbbbbbbbbbbbbbbbbbbbb'
};

const rule = {
  feat: 'CDS', qual: 'product', val: '^Kinase$', color: '#111111', cap: 'Core', fromFile: true
};

const catalog = {
  schema: 3,
  items: [{
    resultKey: 'result-a', resultIndex: 0, resultName: 'a.svg', recordKeys: ['record-a'],
    biologicalFeatures: [{
      recordKey: 'record-a', biologicalFeatureId: 'feature-a', instanceHash: hashes.a,
      type: 'CDS', product: 'Kinase', qualifiers: { product: ['Kinase'] }, renderedLabel: 'Kinase'
    }],
    features: [{ recordKey: 'record-a', biologicalFeatureId: 'feature-a', svgId: 'svg-a' }]
  }, {
    resultKey: 'result-b', resultIndex: 1, resultName: 'b.svg', recordKeys: ['record-b'],
    biologicalFeatures: [{
      recordKey: 'record-b', biologicalFeatureId: 'feature-b', instanceHash: hashes.b,
      type: 'CDS', product: 'Kinase', qualifiers: { product: ['Kinase'] }, renderedLabel: 'Kinase'
    }],
    features: [{ recordKey: 'record-b', biologicalFeatureId: 'feature-b', svgId: 'svg-b' }]
  }]
};

const resultFeature = (item) => ({
  ...item.biologicalFeatures[0],
  resultKey: item.resultKey,
  resultIndex: item.resultIndex,
  resultName: item.resultName,
  renderedFeatureSvgId: item.features[0].svgId,
  legendCaption: 'Core'
});

const beforeFingerprint = styleFingerprint({
  rules: [rule], appliedPaletteName: 'default', appliedPaletteColors: { CDS: '#cccccc' }
});

const readyPlan = ({ value = { kind: 'color', color: '#abcdef' }, scope = 'matching-rule', caption = 'Core' } = {}) => {
  const target = resultFeature(catalog.items[0]);
  const planned = planFeatureFillChange({
    features: catalog.items.map(resultFeature),
    specificRules: [rule],
    paletteColors: { CDS: '#cccccc' },
    resultGenerationKey: 5,
    documentEpoch: 3,
    styleFingerprint: beforeFingerprint,
    exactScopeAvailable: true
  }, {
    targetFeatureKey: stableFeatureTargetKey(target),
    value,
    requestedCaption: caption,
    source: 'popup',
    originResultKey: 'result-a',
    originResultIndex: 0,
    resultGenerationKey: 5,
    documentEpoch: 3,
    styleFingerprint: beforeFingerprint,
    semanticScope: scope
  });
  const candidate = planned.candidates.find((entry) => entry.semanticScope === scope);
  assert.ok(candidate);
  return resolveFeatureFillPlan(planned, candidate.id, {
    conflictChoiceId: planned.status === 'conflict' ? 'separate-caption' : null
  });
};

const createState = () => ({
  manualSpecificRules: clone([rule]),
  appliedPaletteName: ref('default'),
  appliedPaletteColors: ref({ CDS: '#cccccc' }),
  featureColorOverrides: {},
  results: ref([{ name: 'a.svg', content: 'a:#111111' }, { name: 'b.svg', content: 'b:#111111' }]),
  legendEntries: ref([{ caption: 'Core', color: '#111111', featureIds: ['svg-a'] }]),
  validatedStyleFingerprintByResultKey: ref({ 'result-a': beforeFingerprint, 'result-b': beforeFingerprint }),
  semanticStyleRevision: ref(7),
  semanticStyleFingerprint: ref(beforeFingerprint),
  documentEpoch: ref(3),
  resultGenerationKey: ref(5),
  selectedResultIndex: ref(0),
  clickedFeature: ref({ id: 'popup-feature' }),
  selectedFeatures: ref([{ id: 'selected-feature' }])
});

const projectionFactory = ({ failDirection = '', mutations = [] } = {}) => async ({
  direction, results, fillsByTargetKey, affectedResultKeys
}) => {
  mutations.push({ phase: 'prepare', direction, stateResult0: results[0].content });
  if (direction === failDirection) throw new Error('inactive preparation failed');
  assert.deepEqual(affectedResultKeys, ['result-a', 'result-b']);
  const nextResults = results.map((result, index) => {
    const item = catalog.items[index];
    const key = stableFeatureTargetKey(resultFeature(item));
    const color = fillsByTargetKey[key].color;
    return { ...result, content: `${index ? 'b' : 'a'}:${color}` };
  });
  return {
    nextResults,
    selectedLegendEntries: [{
      caption: 'Core',
      color: fillsByTargetKey[stableFeatureTargetKey(resultFeature(catalog.items[0]))].color,
      featureIds: ['svg-a']
    }]
  };
};

test('semantic snapshots preserve rule provenance and derive canonical target overrides', () => {
  const plan = readyPlan();
  const semantic = buildFeatureFillSemanticSnapshots({
    plan,
    catalog,
    rules: [rule],
    appliedPaletteName: 'default',
    appliedPaletteColors: { CDS: '#cccccc' },
    resultFingerprintByKey: { 'result-a': beforeFingerprint, 'result-b': beforeFingerprint }
  });
  assert.equal(semantic.before.rules[0].fromFile, true);
  assert.equal(semantic.after.rules[0].fromFile, true);
  assert.equal(semantic.after.rules[0].color, '#abcdef');
  assert.equal(Object.isFrozen(semantic), true);
  assert.equal(Object.isFrozen(semantic.after.rules), true);
  assert.deepEqual(Object.values(semantic.fillsAfterByTargetKey).map((entry) => entry.color), [
    '#abcdef', '#abcdef'
  ]);

  const derived = deriveFeatureFillOverrides({
    catalog,
    rules: semantic.after.rules,
    paletteColors: semantic.after.appliedPaletteColors
  });
  assert.equal(derived.counters.features, 2);
  assert.equal(derived.counters.ruleResolutionUpperBound, 2);
});

test('feature-type scope persists one future-matching semantic selector, never mass exact rows', () => {
  const semantic = buildFeatureFillSemanticSnapshots({
    plan: readyPlan({ scope: 'feature-type', caption: 'All CDS' }),
    catalog,
    rules: [rule],
    appliedPaletteName: 'default',
    appliedPaletteColors: { CDS: '#cccccc' }
  });
  const semanticRules = semantic.after.rules.filter(
    (entry) => entry.qual === '__gbdraw_semantic_scope__'
  );
  assert.deepEqual(semanticRules, [{
    feat: 'CDS',
    qual: '__gbdraw_semantic_scope__',
    val: 'fs1:feature-type:CDS',
    color: '#abcdef',
    cap: 'All CDS',
    match: 'literal'
  }]);
  assert.equal(
    semantic.after.rules.some((entry) => entry.qual === '__gbdraw_instance_hash__'),
    false
  );
  assert.deepEqual(Object.values(semantic.fillsAfterByTargetKey).map((entry) => entry.color), [
    '#abcdef', '#abcdef'
  ]);

  const futureCatalog = clone(catalog);
  futureCatalog.items[0].biologicalFeatures.push({
    recordKey: 'record-a',
    biologicalFeatureId: 'feature-future',
    instanceHash: 'fi1_cccccccccccccccccccccccccc',
    type: 'CDS',
    product: 'Future protein',
    qualifiers: { product: ['Future protein'] },
    renderedLabel: 'Future protein'
  });
  futureCatalog.items[0].features.push({
    recordKey: 'record-a', biologicalFeatureId: 'feature-future', svgId: 'svg-future'
  });
  const future = deriveFeatureFillOverrides({
    catalog: futureCatalog,
    rules: semantic.after.rules,
    paletteColors: { CDS: '#cccccc' }
  });
  const futureKey = stableFeatureTargetKey({
    resultKey: 'result-a', recordKey: 'record-a', biologicalFeatureId: 'feature-future'
  });
  assert.equal(future.effectiveByTargetKey[futureKey].color, '#abcdef');
});

test('repeated legend-caption edits retain the original semantic selector anchor', () => {
  const first = buildFeatureFillSemanticSnapshots({
    plan: readyPlan({ scope: 'feature-type', caption: 'All CDS' }),
    catalog,
    rules: [rule],
    appliedPaletteName: 'default',
    appliedPaletteColors: { CDS: '#cccccc' }
  });
  const firstFingerprint = styleFingerprint({
    rules: first.after.rules,
    appliedPaletteName: 'default',
    appliedPaletteColors: { CDS: '#cccccc' }
  });
  const captionFeatures = catalog.items.map((item) => ({
    ...resultFeature(item),
    legendCaption: 'All CDS'
  }));
  const planned = planFeatureFillChange({
    features: captionFeatures,
    specificRules: first.after.rules,
    paletteColors: { CDS: '#cccccc' },
    resultGenerationKey: 5,
    documentEpoch: 3,
    styleFingerprint: firstFingerprint,
    exactScopeAvailable: true
  }, {
    targetFeatureKey: stableFeatureTargetKey(captionFeatures[0]),
    value: { kind: 'color', color: '#fedcba' },
    requestedCaption: 'Renamed CDS',
    source: 'legend-editor',
    originResultKey: 'result-a',
    originResultIndex: 0,
    resultGenerationKey: 5,
    documentEpoch: 3,
    styleFingerprint: firstFingerprint,
    semanticScope: 'legend-caption'
  });
  const candidate = planned.candidates.find(
    (entry) => entry.semanticScope === 'legend-caption'
  );
  assert.deepEqual(candidate.durableRuleIntent, {
    kind: 'semantic-selector', selector: 'fs1:feature-type:CDS'
  });
  const resolved = resolveFeatureFillPlan(planned, candidate.id, {
    conflictChoiceId: planned.status === 'conflict' ? 'separate-caption' : null
  });
  const second = buildFeatureFillSemanticSnapshots({
    plan: resolved,
    catalog,
    rules: first.after.rules,
    appliedPaletteName: 'default',
    appliedPaletteColors: { CDS: '#cccccc' }
  });
  const semanticRows = second.after.rules.filter(
    (entry) => entry.qual === '__gbdraw_semantic_scope__'
  );
  assert.equal(semanticRows.length, 1);
  assert.equal(semanticRows[0].val, 'fs1:feature-type:CDS');
  assert.equal(semanticRows[0].color, '#fedcba');
  assert.equal(semanticRows[0].cap, 'Renamed CDS');
  assert.equal(
    second.after.rules.some(
      (entry) => entry.val === 'fs1:base-legend-caption:All%20CDS'
    ),
    false
  );
});

test('retained byte sizing is exact UTF-8 over stable JSON', () => {
  const value = { z: 'é', a: 'genome 🧬' };
  assert.equal(
    exactJsonByteLength(value),
    Buffer.byteLength(JSON.stringify({ a: 'genome 🧬', z: 'é' }), 'utf8')
  );
});

test('multi-Result apply prepares all Results before one synchronous state commit', async () => {
  const state = createState();
  const mutations = [];
  const originalResults = state.results.value;
  const command = await buildFeatureFillCommand({
    plan: readyPlan(),
    state,
    catalog,
    prepareResultProjection: projectionFactory({ mutations }),
    reconcile: () => mutations.push({ phase: 'reconcile', result0: state.results.value[0].content }),
    refreshPresentation: () => mutations.push({ phase: 'refresh', revision: state.semanticStyleRevision.value })
  });

  assert.equal(state.results.value, originalResults, 'preflight must not publish Results');
  assert.equal(state.manualSpecificRules[0].color, '#111111');
  assert.equal(await command.apply(), true);
  assert.deepEqual(state.results.value.map((result) => result.content), ['a:#abcdef', 'b:#abcdef']);
  assert.equal(state.manualSpecificRules[0].color, '#abcdef');
  assert.equal(state.featureColorOverrides['record-a\u0000feature-a'].color, '#abcdef');
  assert.equal(state.semanticStyleRevision.value, 8);
  assert.equal(state.semanticStyleFingerprint.value, styleFingerprint({
    rules: state.manualSpecificRules,
    appliedPaletteName: 'default',
    appliedPaletteColors: { CDS: '#cccccc' }
  }));
  assert.deepEqual(mutations.map((entry) => entry.phase), ['prepare', 'reconcile', 'refresh']);
  assert.equal(command.counters.affectedResults, 2);
  assert.equal(command.counters.projectedResults, 2);
  assert.equal(command.counters.resultArraySwaps, 1);
  assert.ok(command.estimateBytes() > 0);
});

test('Vue-like refs with prototype value accessors are read and written atomically', async () => {
  const state = createState();
  [
    'appliedPaletteName',
    'appliedPaletteColors',
    'results',
    'legendEntries',
    'validatedStyleFingerprintByResultKey',
    'semanticStyleRevision',
    'semanticStyleFingerprint',
    'documentEpoch',
    'resultGenerationKey',
    'clickedFeature',
    'selectedFeatures'
  ].forEach((key) => {
    state[key] = new PrototypeRef(state[key].value);
    assert.equal(Object.prototype.hasOwnProperty.call(state[key], 'value'), false);
  });

  const command = await buildFeatureFillCommand({
    plan: readyPlan(),
    state,
    catalog,
    prepareResultProjection: projectionFactory()
  });
  assert.equal(await command.apply(), true);
  assert.deepEqual(state.results.value.map((result) => result.content), ['a:#abcdef', 'b:#abcdef']);
  assert.equal(state.semanticStyleFingerprint.value, styleFingerprint({
    rules: state.manualSpecificRules,
    appliedPaletteName: 'default',
    appliedPaletteColors: { CDS: '#cccccc' }
  }));
  assert.equal(await command.revert(), true);
  assert.deepEqual(state.results.value.map((result) => result.content), ['a:#111111', 'b:#111111']);
});

test('Vue-like Result array proxying preserves post-commit invariants', async () => {
  const state = createState();
  state.results = new ProxyingArrayRef(state.results.value);
  const command = await buildFeatureFillCommand({
    plan: readyPlan(),
    state,
    catalog,
    prepareResultProjection: projectionFactory()
  });
  assert.equal(await command.apply(), true);
  assert.deepEqual(state.results.value.map((result) => result.content), ['a:#abcdef', 'b:#abcdef']);
  assert.equal(await command.revert(), true);
  assert.deepEqual(state.results.value.map((result) => result.content), ['a:#111111', 'b:#111111']);
});

test('inactive preparation failure leaves every domain byte-for-byte unchanged', async () => {
  const state = createState();
  const references = {
    rules: state.manualSpecificRules,
    overrides: state.featureColorOverrides,
    results: state.results.value,
    legend: state.legendEntries.value,
    ledger: state.validatedStyleFingerprintByResultKey.value
  };
  const before = clone({
    rules: state.manualSpecificRules,
    overrides: state.featureColorOverrides,
    results: state.results.value,
    legend: state.legendEntries.value,
    ledger: state.validatedStyleFingerprintByResultKey.value,
    revision: state.semanticStyleRevision.value,
    fingerprint: state.semanticStyleFingerprint.value
  });
  await assert.rejects(() => buildFeatureFillCommand({
    plan: readyPlan(),
    state,
    catalog,
    prepareResultProjection: projectionFactory({ failDirection: 'apply' })
  }), /inactive preparation failed/);
  assert.deepEqual({
    rules: state.manualSpecificRules,
    overrides: state.featureColorOverrides,
    results: state.results.value,
    legend: state.legendEntries.value,
    ledger: state.validatedStyleFingerprintByResultKey.value,
    revision: state.semanticStyleRevision.value,
    fingerprint: state.semanticStyleFingerprint.value
  }, before);
  assert.equal(state.manualSpecificRules, references.rules);
  assert.equal(state.featureColorOverrides, references.overrides);
  assert.equal(state.results.value, references.results);
  assert.equal(state.legendEntries.value, references.legend);
  assert.equal(state.validatedStyleFingerprintByResultKey.value, references.ledger);
});

test('inherit removes only selected exact rule and exposes the broader inherited rule', () => {
  const exactPlan = readyPlan({ value: { kind: 'color', color: '#222222' }, scope: 'single', caption: 'One' });
  const created = buildFeatureFillSemanticSnapshots({
    plan: exactPlan,
    catalog,
    rules: [rule],
    appliedPaletteName: 'default',
    appliedPaletteColors: { CDS: '#cccccc' }
  });
  assert.equal(created.after.rules.length, 2);
  assert.equal(created.after.rules[0].qual, '__gbdraw_instance_hash__');
  assert.equal(created.after.rules[1].color, '#111111');

  const inheritedPlan = clone(exactPlan);
  inheritedPlan.intent.value = { kind: 'inherit' };
  inheritedPlan.semanticAfter.value = { kind: 'inherit' };
  const inherited = buildFeatureFillSemanticSnapshots({
    plan: inheritedPlan,
    catalog,
    rules: created.after.rules,
    appliedPaletteName: 'default',
    appliedPaletteColors: { CDS: '#cccccc' }
  });
  assert.deepEqual(inherited.after.rules, [{
    feat: 'CDS', qual: 'product', val: '^Kinase$', color: '#111111', cap: 'Core', fromFile: true
  }]);
  assert.equal(inherited.fillsAfterByTargetKey[stableFeatureTargetKey(resultFeature(catalog.items[0]))].color, '#111111');

  const duplicateExact = { ...created.after.rules[0], color: '#333333', cap: 'Shadow', fromFile: true };
  const inheritedWithDuplicate = buildFeatureFillSemanticSnapshots({
    plan: inheritedPlan,
    catalog,
    rules: [created.after.rules[0], duplicateExact, rule],
    appliedPaletteName: 'default',
    appliedPaletteColors: { CDS: '#cccccc' }
  });
  assert.deepEqual(inheritedWithDuplicate.after.rules, [duplicateExact, rule]);
});

test('an exact edit without a replacement caption splits a conflicting inherited legend caption', () => {
  const semantic = buildFeatureFillSemanticSnapshots({
    plan: readyPlan({ value: { kind: 'color', color: '#222222' }, scope: 'single', caption: '' }),
    catalog,
    rules: [rule],
    appliedPaletteName: 'default',
    appliedPaletteColors: { CDS: '#cccccc' }
  });
  assert.equal(semantic.after.rules[0].qual, '__gbdraw_instance_hash__');
  assert.equal(semantic.after.rules[0].cap, 'Core (2)');
});

test('palette fallback captions ignore stale catalogue projections', () => {
  const staleCatalog = clone(catalog);
  staleCatalog.items[0].biologicalFeatures[0].effectiveLegendCaption = 'Retired exact caption';
  const derived = deriveFeatureFillOverrides({
    catalog: staleCatalog,
    rules: [],
    paletteColors: { CDS: '#cccccc' }
  });
  assert.equal(
    derived.effectiveByTargetKey[stableFeatureTargetKey(resultFeature(staleCatalog.items[0]))].caption,
    'CDS'
  );
});

test('explicit no-fill remains distinct from inherit in rules and target projections', () => {
  const semantic = buildFeatureFillSemanticSnapshots({
    plan: readyPlan({ value: { kind: 'none' }, scope: 'matching-rule', caption: 'Core' }),
    catalog,
    rules: [rule],
    appliedPaletteName: 'default',
    appliedPaletteColors: { CDS: '#cccccc' }
  });
  assert.equal(semantic.after.rules.length, 1);
  assert.equal(semantic.after.rules[0].color, 'none');
  assert.deepEqual(
    Object.values(semantic.fillsAfterByTargetKey).map((entry) => entry.color),
    ['none', 'none']
  );
});

test('no-op and stale commands do not mutate or create History entries', async () => {
  const noOpState = createState();
  const noOpPlan = readyPlan({ value: { kind: 'color', color: '#111111' }, scope: 'matching-rule' });
  const noOp = await buildFeatureFillCommand({
    plan: noOpPlan,
    state: noOpState,
    catalog,
    prepareResultProjection: projectionFactory()
  });
  assert.equal(noOp.noop, true);
  assert.equal(noOp.estimateBytes(), 0);
  assert.equal(await noOp.apply(), false);

  const staleState = createState();
  const command = await buildFeatureFillCommand({
    plan: readyPlan(), state: staleState, catalog,
    prepareResultProjection: projectionFactory()
  });
  staleState.resultGenerationKey.value += 1;
  assert.equal(await command.apply(), false);
  assert.equal(staleState.manualSpecificRules[0].color, '#111111');

  const staleTargetState = createState();
  const mutableCatalog = clone(catalog);
  const staleTargetCommand = await buildFeatureFillCommand({
    plan: readyPlan(), state: staleTargetState, catalog: mutableCatalog,
    prepareResultProjection: projectionFactory()
  });
  mutableCatalog.items[1].biologicalFeatures.length = 0;
  assert.equal(await staleTargetCommand.apply(), false);
  assert.equal(staleTargetState.manualSpecificRules[0].color, '#111111');
  assert.deepEqual(staleTargetState.results.value.map((result) => result.content), [
    'a:#111111', 'b:#111111'
  ]);

  const staleArtifactState = createState();
  const staleArtifactCommand = await buildFeatureFillCommand({
    plan: readyPlan(), state: staleArtifactState, catalog,
    prepareResultProjection: projectionFactory()
  });
  staleArtifactState.results.value[1].content = 'b:externally-changed';
  assert.equal(await staleArtifactCommand.apply(), false);
  assert.equal(staleArtifactState.manualSpecificRules[0].color, '#111111');
  assert.equal(staleArtifactState.results.value[1].content, 'b:externally-changed');
});

test('History apply/Undo/Redo uses the same atomic projection and compensation restores exact state', async () => {
  const state = createState();
  const buildIntent = () => ({
    rules: clone(state.manualSpecificRules),
    results: clone(state.results.value),
    revision: state.semanticStyleRevision.value
  });
  const history = createHistoryManager({
    buildIntent,
    applyIntent: () => {},
    buildCheckpoint: () => ({}),
    applyCheckpoint: () => {}
  });
  await history.captureBaseline();
  const command = await buildFeatureFillCommand({
    plan: readyPlan(), state, catalog,
    prepareResultProjection: projectionFactory()
  });
  assert.equal(await history.runUndoableCommand('Change feature fill', () => command), true);
  assert.equal(history.getUndoCount(), 1);
  assert.deepEqual(state.results.value.map((result) => result.content), ['a:#abcdef', 'b:#abcdef']);
  assert.equal(await history.undo(), true);
  assert.deepEqual(state.results.value.map((result) => result.content), ['a:#111111', 'b:#111111']);
  assert.equal(state.semanticStyleRevision.value, 9, 'Undo is a new accepted transition');
  assert.equal(await history.redo(), true);
  assert.deepEqual(state.results.value.map((result) => result.content), ['a:#abcdef', 'b:#abcdef']);
  assert.equal(state.semanticStyleRevision.value, 10);
});

test('History rejects an oversized fill command before any semantic or Result mutation', async () => {
  const state = createState();
  const command = await buildFeatureFillCommand({
    plan: readyPlan(), state, catalog,
    prepareResultProjection: projectionFactory()
  });
  const results = state.results.value;
  const history = createHistoryManager({
    buildIntent: () => ({ revision: state.semanticStyleRevision.value }),
    applyIntent: () => {},
    buildCheckpoint: () => ({}),
    applyCheckpoint: () => {},
    maxBytes: command.estimateBytes() - 1
  });
  await history.captureBaseline();
  const warn = console.warn;
  console.warn = () => {};
  try {
    assert.equal(await history.runUndoableCommand('Change feature fill', () => command), false);
  } finally {
    console.warn = warn;
  }
  assert.equal(state.results.value, results);
  assert.equal(state.manualSpecificRules[0].color, '#111111');
  assert.equal(state.semanticStyleRevision.value, 7);
  assert.equal(history.getUndoCount(), 0);
});

test('a synchronous commit exception self-rolls back structural state and counters remain honest', async () => {
  const state = createState();
  const beforeResults = state.results.value;
  const command = await buildFeatureFillCommand({
    plan: readyPlan(), state, catalog,
    prepareResultProjection: projectionFactory(),
    reconcile: () => { throw new Error('reconcile failed'); }
  });
  await assert.rejects(() => command.apply(), /reconcile failed/);
  assert.equal(state.results.value, beforeResults);
  assert.equal(state.manualSpecificRules[0].color, '#111111');
  assert.equal(state.semanticStyleRevision.value, 7);
  assert.equal(command.counters.targetFeatures, 2);
  assert.equal(command.counters.affectedResults, 2);
});

test('plan fingerprint, Result ledger, and mounted ownership are revalidated around async preflight', async () => {
  const staleStyleState = createState();
  staleStyleState.manualSpecificRules[0].color = '#222222';
  staleStyleState.semanticStyleFingerprint.value = styleFingerprint({
    rules: staleStyleState.manualSpecificRules,
    appliedPaletteName: 'default',
    appliedPaletteColors: { CDS: '#cccccc' }
  });
  let preparations = 0;
  await assert.rejects(() => buildFeatureFillCommand({
    plan: readyPlan(),
    state: staleStyleState,
    catalog,
    prepareResultProjection: async () => {
      preparations += 1;
      return { nextResults: staleStyleState.results.value };
    }
  }), /became stale/);
  assert.equal(preparations, 0, 'a stale semantic plan must fail before projection');

  const staleLedgerState = createState();
  staleLedgerState.validatedStyleFingerprintByResultKey.value = { 'result-a': beforeFingerprint };
  await assert.rejects(() => buildFeatureFillCommand({
    plan: readyPlan(),
    state: staleLedgerState,
    catalog,
    prepareResultProjection: projectionFactory()
  }), /became stale/);

  const mountState = createState();
  const mountedA = { resultIndex: 0, resultKey: 'result-a', svg: { id: 'svg-a' } };
  const mountedB = { resultIndex: 1, resultKey: 'result-b', svg: { id: 'svg-b' } };
  let mounted = mountedA;
  let releasePreparation;
  let preparationStarted;
  const started = new Promise((resolve) => { preparationStarted = resolve; });
  const release = new Promise((resolve) => { releasePreparation = resolve; });
  const building = buildFeatureFillCommand({
    plan: readyPlan(),
    state: mountState,
    catalog,
    getMountedContext: () => mounted,
    prepareResultProjection: async (input) => {
      preparationStarted();
      await release;
      return projectionFactory()(input);
    }
  });
  await started;
  mounted = mountedB;
  releasePreparation();
  await assert.rejects(() => building, /became stale/);
  assert.equal(mountState.results.value[0].content, 'a:#111111');
  assert.equal(mountState.manualSpecificRules[0].color, '#111111');
});

test('History post-operation failures compensate apply, Undo, and Redo to exact pre-attempt state', async () => {
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
        rules: clone(state.manualSpecificRules),
        results: clone(state.results.value),
        revision: state.semanticStyleRevision.value
      };
    },
    applyIntent: () => {},
    buildCheckpoint: () => ({}),
    applyCheckpoint: () => {}
  });
  await history.captureBaseline();

  const firstCommand = await buildFeatureFillCommand({
    plan: readyPlan(), state, catalog,
    prepareResultProjection: projectionFactory(),
    refreshPresentation: () => {
      state.clickedFeature.value = null;
      return true;
    }
  });
  const beforeResults = state.results.value;
  const beforeRule = state.manualSpecificRules[0];
  const beforePaletteColors = state.appliedPaletteColors.value;
  const beforeLegend = state.legendEntries.value;
  const beforeLedger = state.validatedStyleFingerprintByResultKey.value;
  const beforeClickedFeature = state.clickedFeature.value;
  failCapture = true;
  await assert.rejects(
    () => history.runUndoableCommand('Change feature fill', () => firstCommand),
    /intent capture failed/
  );
  assert.equal(state.results.value, beforeResults);
  assert.equal(state.manualSpecificRules[0], beforeRule);
  assert.equal(state.appliedPaletteColors.value, beforePaletteColors);
  assert.equal(state.legendEntries.value, beforeLegend);
  assert.equal(state.validatedStyleFingerprintByResultKey.value, beforeLedger);
  assert.equal(state.clickedFeature.value, beforeClickedFeature);
  assert.equal(state.selectedResultIndex.value, 0);
  assert.equal(state.manualSpecificRules[0].color, '#111111');
  assert.equal(state.semanticStyleRevision.value, 7);
  assert.equal(history.getUndoCount(), 0);
  assert.equal(history.getRedoCount(), 0);

  const command = await buildFeatureFillCommand({
    plan: readyPlan(), state, catalog,
    prepareResultProjection: projectionFactory(),
    refreshPresentation: () => {
      state.clickedFeature.value = null;
      return true;
    }
  });
  assert.equal(await history.runUndoableCommand('Change feature fill', () => command), true);
  const afterApplyResults = state.results.value;
  failCapture = true;
  await assert.rejects(() => history.undo(), /intent capture failed/);
  assert.equal(state.results.value, afterApplyResults);
  assert.equal(state.manualSpecificRules[0].color, '#abcdef');
  assert.equal(state.semanticStyleRevision.value, 8);
  assert.equal(history.getUndoCount(), 1);
  assert.equal(history.getRedoCount(), 0);
  assert.equal(state.selectedResultIndex.value, 0);

  assert.equal(await history.undo(), true);
  const afterUndoResults = state.results.value;
  failCapture = true;
  await assert.rejects(() => history.redo(), /intent capture failed/);
  assert.equal(state.results.value, afterUndoResults);
  assert.equal(state.manualSpecificRules[0].color, '#111111');
  assert.equal(state.semanticStyleRevision.value, 9);
  assert.equal(history.getUndoCount(), 0);
  assert.equal(history.getRedoCount(), 1);
  assert.equal(state.selectedResultIndex.value, 0);
});
