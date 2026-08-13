import assert from 'node:assert/strict';
import test from 'node:test';

import {
  buildFeatureBulkStyleCommand,
  deriveFeatureBulkStyleChange,
  normalizeFeatureBulkStyleSnapshot,
  replaceFeatureBulkStyleSnapshot
} from '../../gbdraw/web/js/app/feature-editor/bulk-style-command.js';
import { deriveFeatureFillOverrides } from '../../gbdraw/web/js/app/feature-editor/fill-command.js';
import { stableFeatureTargetKey } from '../../gbdraw/web/js/app/feature-editor/fill-scope-plan.js';
import { deriveUsedFeatureFillGroupsByResult } from '../../gbdraw/web/js/app/legend/feature-fill-projection.js';
import { createHistoryManager } from '../../gbdraw/web/js/services/history.js';
import { styleFingerprint } from '../../gbdraw/web/js/services/style-revision.js';

const ref = (value) => ({ value });
const clone = (value) => structuredClone(value);

const kinaseRule = {
  feat: 'CDS', qual: 'product', val: '^Kinase$', color: '#111111', cap: 'Core', fromFile: true
};

const catalog = {
  schema: 3,
  items: [{
    resultKey: 'result-a', resultIndex: 0, resultName: 'a.svg',
    biologicalFeatures: [{
      recordKey: 'record-a', biologicalFeatureId: 'feature-a', type: 'CDS',
      qualifiers: { product: ['Kinase'] }
    }],
    features: [{ recordKey: 'record-a', biologicalFeatureId: 'feature-a', svgId: 'svg-a' }]
  }, {
    resultKey: 'result-b', resultIndex: 1, resultName: 'b.svg',
    biologicalFeatures: [{
      recordKey: 'record-b', biologicalFeatureId: 'feature-b', type: 'CDS',
      qualifiers: { product: ['Kinase'] }
    }],
    features: [{ recordKey: 'record-b', biologicalFeatureId: 'feature-b', svgId: 'svg-b' }]
  }, {
    resultKey: 'result-c', resultIndex: 2, resultName: 'c.svg',
    biologicalFeatures: [{
      recordKey: 'record-c', biologicalFeatureId: 'feature-c', type: 'tRNA',
      qualifiers: { product: ['tRNA-Lys'] }
    }],
    features: [{ recordKey: 'record-c', biologicalFeatureId: 'feature-c', svgId: 'svg-c' }]
  }]
};

const baseStyle = normalizeFeatureBulkStyleSnapshot({
  rules: [kinaseRule],
  appliedPaletteName: 'default',
  appliedPaletteColors: { CDS: '#cccccc', tRNA: '#22aa66' }
});

const fakeRoot = (resultKey, payload = '') => ({
  resultKey,
  payload,
  cloneNode() { return fakeRoot(this.resultKey, this.payload); }
});

const sideForStyle = (style) => deriveFeatureBulkStyleChange({
  catalog,
  before: style,
  after: style
}).before;

const resultContent = (side, resultKey) => {
  const item = catalog.items.find((entry) => entry.resultKey === resultKey);
  const fills = item.biologicalFeatures.map((feature) => {
    const targetKey = stableFeatureTargetKey({ ...feature, resultKey });
    const value = side.effectiveByTargetKey[targetKey];
    return `${feature.biologicalFeatureId}:${value.color}:${value.caption}`;
  });
  const legend = side.legendProjections
    .find((projection) => projection.resultKey === resultKey)?.entries
    .map((entry) => `${entry.caption}:${entry.color}:${entry.featureIds.join(',')}`) || [];
  return `${resultKey}|${fills.join('|')}|legend=${legend.join(';')}`;
};

const createState = (style = baseStyle) => {
  const side = sideForStyle(style);
  const fingerprint = styleFingerprint(style);
  return {
    manualSpecificRules: clone(style.rules),
    appliedPaletteName: ref(style.appliedPaletteName),
    appliedPaletteColors: ref(clone(style.appliedPaletteColors)),
    selectedPalette: ref(style.appliedPaletteName),
    currentColors: ref(clone(style.appliedPaletteColors)),
    pendingPaletteName: ref(''),
    pendingPaletteColors: ref({}),
    adv: { legend_box_size: 10, legend_font_size: 10 },
    featureColorOverrides: clone(side.overrides),
    results: ref(catalog.items.map((item) => ({
      name: item.resultName,
      content: resultContent(side, item.resultKey)
    }))),
    legendEntries: ref(clone(side.legendProjections[0].entries)),
    manualLegendEntries: ref([]),
    validatedStyleFingerprintByResultKey: ref(Object.fromEntries(
      catalog.items.map((item) => [item.resultKey, fingerprint])
    )),
    semanticStyleRevision: ref(7),
    semanticStyleFingerprint: ref(fingerprint),
    documentEpoch: ref(3),
    resultGenerationKey: ref(5),
    featureCatalog: ref(catalog),
    selectedResultIndex: ref(0),
    clickedFeature: ref({ id: 'clicked' }),
    selectedFeatures: ref([{ id: 'selected' }])
  };
};

const createHarness = ({
  state,
  failLegendResult = '',
  failResult = '',
  failCommit = false,
  mountedIndex = 0
} = {}) => {
  const mounted = {
    resultIndex: mountedIndex,
    resultKey: catalog.items[mountedIndex]?.resultKey || '',
    svg: mountedIndex === null ? null : fakeRoot(catalog.items[mountedIndex].resultKey, 'mounted-before')
  };
  const calls = { legend: [], result: [], commit: 0, reconcile: 0, restore: 0 };
  const prepareLegendProjection = async ({ projections, affectedResultKeys, mounted: owner }) => {
    calls.legend.push([...affectedResultKeys]);
    const staged = new Map();
    projections.forEach((projection) => {
      if (projection.resultKey === failLegendResult) throw new Error(`legend failed: ${projection.resultKey}`);
      staged.set(projection.resultKey, {
        svg: fakeRoot(projection.resultKey, `legend:${projection.entries.length}`)
      });
    });
    return {
      staged,
      selectedEntries: clone(
        projections.find((projection) => projection.resultKey === owner.resultKey)?.entries || []
      ),
      counters: { measurements: 0 },
      retainedForHistory: { resultKeys: [...affectedResultKeys] }
    };
  };
  const prepareResultProjection = async ({
    to,
    results,
    affectedResultKeys,
    preparedSvgByResultKey,
    targetFeatureKeysByResult,
    mounted: owner
  }) => {
    calls.result.push([...affectedResultKeys]);
    const affected = new Set(affectedResultKeys);
    const admissionMetadataByResultKey = {};
    const nextResults = results.map((result, resultIndex) => {
      const resultKey = catalog.items[resultIndex].resultKey;
      if (!affected.has(resultKey)) return result;
      if (resultKey === failResult) throw new Error(`Result failed: ${resultKey}`);
      assert.ok(preparedSvgByResultKey.get(resultKey)?.cloneNode);
      assert.ok(targetFeatureKeysByResult[resultKey]?.length > 0);
      admissionMetadataByResultKey[resultKey] = {
        before: { token: `before:${resultKey}` },
        after: { token: `after:${resultKey}` }
      };
      return { ...result, content: resultContent(to, resultKey) };
    });
    return {
      previousResults: results,
      nextResults,
      admissionMetadataByResultKey,
      preparedMountedSvg: affected.has(owner.resultKey)
        ? preparedSvgByResultKey.get(owner.resultKey)
        : null,
      counters: {
        affectedResults: affected.size,
        mountedSerializations: affected.has(owner.resultKey) ? 1 : 0,
        detachedPasses: affected.size - (affected.has(owner.resultKey) ? 1 : 0),
        changedResults: affected.size
      },
      retainedForHistory: { affected: [...affected] }
    };
  };
  const captureRuntimeState = () => ({ payload: mounted.svg?.payload || null });
  const restoreRuntimeState = (snapshot) => {
    calls.restore += 1;
    if (mounted.svg) mounted.svg.payload = snapshot?.payload || '';
    return true;
  };
  const commitMountedProjection = ({ prepared }) => {
    calls.commit += 1;
    if (mounted.svg) mounted.svg.payload = prepared.preparedMountedSvg?.payload || 'committed';
    if (failCommit) throw new Error('mounted commit failed');
    return true;
  };
  const options = {
    prepareLegendProjection,
    prepareResultProjection,
    getMountedContext: () => ({ ...mounted }),
    commitMountedProjection,
    captureRuntimeState,
    restoreRuntimeState,
    restoreMountedProjection: () => true,
    reconcile: () => {
      calls.reconcile += 1;
      return true;
    }
  };
  return { mounted, calls, options };
};

const buildCommand = async ({
  state = createState(),
  before = baseStyle,
  after,
  writerKind = 'rule-edit',
  harness = null,
  extra = {}
} = {}) => {
  const activeHarness = harness || createHarness({ state });
  const command = await buildFeatureBulkStyleCommand({
    state,
    catalog,
    before,
    after,
    writerKind,
    ...activeHarness.options,
    ...extra
  });
  return { state, harness: activeHarness, command };
};

test('pure bulk derivation is immutable and identifies all and only affected Results', () => {
  const before = clone(baseStyle);
  const after = replaceFeatureBulkStyleSnapshot(before, {
    rules: [{ ...kinaseRule, color: '#abcdef' }]
  });
  Object.freeze(before.rules[0]);
  Object.freeze(before.rules);
  Object.freeze(before);
  const change = deriveFeatureBulkStyleChange({ catalog, before, after });
  assert.deepEqual(change.affectedResultKeys, ['result-a', 'result-b']);
  assert.deepEqual(change.unaffectedResultKeys, ['result-c']);
  assert.deepEqual(Object.keys(change.targetFeatureKeysByResult), ['result-a', 'result-b']);
  assert.equal(change.counters.changedRenderedFeatures, 2);
  assert.equal(Object.isFrozen(change), true);
  assert.equal(Object.isFrozen(change.after.rules), true);
  assert.equal(change.after.rules[0].fromFile, true);
  assert.equal(before.rules[0].color, '#111111');
});

test('non-Feature palette ownership can extend one atomic snapshot to every affected Result', () => {
  const before = replaceFeatureBulkStyleSnapshot(baseStyle, {
    appliedPaletteColors: {
      ...baseStyle.appliedPaletteColors,
      gc_content: '#777777'
    }
  });
  const after = replaceFeatureBulkStyleSnapshot(before, {
    appliedPaletteColors: {
      ...before.appliedPaletteColors,
      gc_content: '#23ab45'
    }
  });
  const resultKeys = catalog.items.map((item) => item.resultKey);
  const change = deriveFeatureBulkStyleChange({
    catalog,
    before,
    after,
    nonFeatureResultKeys: resultKeys
  });
  assert.deepEqual(change.affectedResultKeys, resultKeys);
  assert.deepEqual(change.featureAffectedResultKeys, []);
  assert.deepEqual(change.nonFeatureAffectedResultKeys, resultKeys);
  assert.deepEqual(change.unaffectedResultKeys, []);
  assert.deepEqual(change.targetFeatureKeysByResult, {});
  assert.equal(change.counters.nonFeatureAffectedResults, 3);
});

test('multi-Result rule edit preflights every legend and Result before one commit', async () => {
  const after = replaceFeatureBulkStyleSnapshot(baseStyle, {
    rules: [{ ...kinaseRule, color: '#abcdef' }]
  });
  const state = createState();
  const beforeResults = state.results.value;
  const beforeUnaffected = beforeResults[2];
  const { command, harness } = await buildCommand({ state, after });
  assert.equal(state.results.value, beforeResults, 'preflight cannot publish state');
  assert.deepEqual(harness.calls.legend, [['result-a', 'result-b']]);
  assert.deepEqual(harness.calls.result, [['result-a', 'result-b']]);
  assert.equal(await command.apply(), true);
  assert.equal(state.manualSpecificRules[0].color, '#abcdef');
  assert.notEqual(state.results.value, beforeResults);
  assert.equal(state.results.value[2], beforeUnaffected);
  assert.equal(state.results.value[0].content.includes('#abcdef'), true);
  assert.equal(state.results.value[1].content.includes('#abcdef'), true);
  assert.equal(state.results.value[2].content.includes('#22aa66'), true);
  assert.equal(state.semanticStyleRevision.value, 8);
  assert.equal(state.semanticStyleFingerprint.value, styleFingerprint(after));
  assert.deepEqual(Object.values(state.validatedStyleFingerprintByResultKey.value), [
    styleFingerprint(after), styleFingerprint(after), styleFingerprint(after)
  ]);
  assert.equal(command.metadata.writerKind, 'rule-edit');
  assert.equal(command.counters.affectedResults, 2);
  assert.equal(command.counters.mountedSerializations, 1);
  assert.equal(command.counters.detachedPasses, 1);
  assert.equal(command.counters.resultArraySwaps, 1);
  assert.ok(command.estimateBytes() > state.results.value[0].content.length);
});

test('complete snapshots support every bulk writer shape without writer-specific mutation branches', async (t) => {
  const cases = [{
    name: 'rule field edit',
    writerKind: 'rule-edit',
    before: baseStyle,
    after: replaceFeatureBulkStyleSnapshot(baseStyle, {
      rules: [{ ...kinaseRule, color: '#222222' }]
    }),
    expected: ['result-a', 'result-b']
  }, {
    name: 'rule add',
    writerKind: 'rule-add',
    before: baseStyle,
    after: replaceFeatureBulkStyleSnapshot(baseStyle, {
      rules: [kinaseRule, {
        feat: 'tRNA', qual: 'product', val: '^tRNA-Lys$', color: '#8844cc', cap: 'Transfer'
      }]
    }),
    expected: ['result-c']
  }, {
    name: 'rule remove',
    writerKind: 'rule-remove',
    before: baseStyle,
    after: replaceFeatureBulkStyleSnapshot(baseStyle, { rules: [] }),
    expected: ['result-a', 'result-b']
  }, {
    name: 'rule reorder',
    writerKind: 'rule-reorder',
    before: normalizeFeatureBulkStyleSnapshot({
      ...baseStyle,
      rules: [
        kinaseRule,
        { feat: 'CDS', qual: 'product', val: 'Kinase', color: '#3377aa', cap: 'Second' }
      ]
    }),
    after: normalizeFeatureBulkStyleSnapshot({
      ...baseStyle,
      rules: [
        { feat: 'CDS', qual: 'product', val: 'Kinase', color: '#3377aa', cap: 'Second' },
        kinaseRule
      ]
    }),
    expected: ['result-a', 'result-b']
  }, {
    name: 'preset replacement',
    writerKind: 'preset',
    before: baseStyle,
    after: replaceFeatureBulkStyleSnapshot(baseStyle, {
      rules: [{ feat: '*', qual: 'product', val: '.+', color: '#dd7733', cap: 'Preset' }]
    }),
    expected: ['result-a', 'result-b', 'result-c']
  }, {
    name: 'TSV import replacement with provenance',
    writerKind: 'tsv-import',
    before: baseStyle,
    after: replaceFeatureBulkStyleSnapshot(baseStyle, {
      rules: [{
        feat: '*', qual: 'product', val: '.+', color: '#55aaee', cap: 'Imported', fromFile: true
      }]
    }),
    expected: ['result-a', 'result-b', 'result-c']
  }, {
    name: 'palette/default acceptance',
    writerKind: 'palette',
    before: normalizeFeatureBulkStyleSnapshot({
      rules: [], appliedPaletteName: 'default',
      appliedPaletteColors: { CDS: '#cccccc', tRNA: '#22aa66' }
    }),
    after: normalizeFeatureBulkStyleSnapshot({
      rules: [], appliedPaletteName: 'ocean',
      appliedPaletteColors: { CDS: '#116699', tRNA: '#44aacc' }
    }),
    expected: ['result-a', 'result-b', 'result-c']
  }];

  for (const scenario of cases) {
    await t.test(scenario.name, async () => {
      const state = createState(scenario.before);
      const change = deriveFeatureBulkStyleChange({
        catalog, before: scenario.before, after: scenario.after
      });
      assert.deepEqual(change.affectedResultKeys, scenario.expected);
      const harness = createHarness({ state });
      const { command } = await buildCommand({
        state,
        before: scenario.before,
        after: scenario.after,
        writerKind: scenario.writerKind,
        harness
      });
      assert.equal(await command.apply(), true);
      assert.equal(command.metadata.writerKind, scenario.writerKind);
      assert.deepEqual(state.manualSpecificRules, clone(scenario.after.rules));
      assert.deepEqual(state.appliedPaletteColors.value, clone(scenario.after.appliedPaletteColors));
      assert.equal(state.semanticStyleFingerprint.value, styleFingerprint(scenario.after));
    });
  }
});

test('provenance-only snapshot is undoable while retaining unchanged Result objects', async () => {
  const before = normalizeFeatureBulkStyleSnapshot({
    ...baseStyle,
    rules: [{ ...kinaseRule, fromFile: false }]
  });
  const after = replaceFeatureBulkStyleSnapshot(before, {
    rules: [{ ...kinaseRule, fromFile: true }]
  });
  assert.equal(styleFingerprint(before), styleFingerprint(after));
  const state = createState(before);
  const results = state.results.value;
  const harness = createHarness({ state });
  const { command } = await buildCommand({ state, before, after, harness, writerKind: 'tsv-import' });
  assert.equal(command.noop, undefined);
  assert.deepEqual(harness.calls.legend, []);
  assert.deepEqual(harness.calls.result, []);
  assert.equal(await command.apply(), true);
  assert.equal(state.manualSpecificRules[0].fromFile, true);
  assert.notEqual(state.results.value, results, 'one shallow Result-array transaction is still published');
  assert.deepEqual(state.results.value, results);
  assert.equal(state.results.value[0], results[0]);
  assert.equal(state.semanticStyleRevision.value, 8);
});

test('exact no-op skips preflight, Result publication, revision, and History work', async () => {
  const state = createState();
  const results = state.results.value;
  const currentColors = state.currentColors.value;
  const pendingPaletteColors = state.pendingPaletteColors.value;
  const harness = createHarness({ state });
  const { command } = await buildCommand({ state, after: baseStyle, harness });
  assert.equal(command.noop, true);
  assert.equal(command.estimateBytes(), 0);
  assert.deepEqual(harness.calls.legend, []);
  assert.deepEqual(harness.calls.result, []);
  assert.equal(await command.apply(), false);
  assert.equal(state.results.value, results);
  assert.equal(state.currentColors.value, currentColors);
  assert.equal(state.pendingPaletteColors.value, pendingPaletteColors);
  assert.equal(state.semanticStyleRevision.value, 7);
});

test('failure in an inactive legend or Result leaves every state reference unchanged', async (t) => {
  const after = replaceFeatureBulkStyleSnapshot(baseStyle, {
    rules: [{ ...kinaseRule, color: '#abcdef' }]
  });
  for (const scenario of [
    { label: 'legend', options: { failLegendResult: 'result-b' }, pattern: /legend failed/ },
    { label: 'Result', options: { failResult: 'result-b' }, pattern: /Result failed/ }
  ]) {
    await t.test(scenario.label, async () => {
      const state = createState();
      const references = {
        rules: state.manualSpecificRules,
        overrides: state.featureColorOverrides,
        results: state.results.value,
        legend: state.legendEntries.value,
        ledger: state.validatedStyleFingerprintByResultKey.value
      };
      const snapshot = clone({
        rules: state.manualSpecificRules,
        overrides: state.featureColorOverrides,
        results: state.results.value,
        legend: state.legendEntries.value,
        ledger: state.validatedStyleFingerprintByResultKey.value,
        revision: state.semanticStyleRevision.value,
        fingerprint: state.semanticStyleFingerprint.value
      });
      const harness = createHarness({ state, ...scenario.options });
      await assert.rejects(() => buildCommand({ state, after, harness }), scenario.pattern);
      assert.deepEqual({
        rules: state.manualSpecificRules,
        overrides: state.featureColorOverrides,
        results: state.results.value,
        legend: state.legendEntries.value,
        ledger: state.validatedStyleFingerprintByResultKey.value,
        revision: state.semanticStyleRevision.value,
        fingerprint: state.semanticStyleFingerprint.value
      }, snapshot);
      assert.equal(state.manualSpecificRules, references.rules);
      assert.equal(state.featureColorOverrides, references.overrides);
      assert.equal(state.results.value, references.results);
      assert.equal(state.legendEntries.value, references.legend);
      assert.equal(state.validatedStyleFingerprintByResultKey.value, references.ledger);
    });
  }
});

test('synchronous commit failure restores semantic state, Results, mounted state, and revision exactly', async () => {
  const after = replaceFeatureBulkStyleSnapshot(baseStyle, {
    rules: [{ ...kinaseRule, color: '#abcdef' }],
    appliedPaletteName: 'ocean',
    appliedPaletteColors: { CDS: '#116699', tRNA: '#44aacc' }
  });
  const state = createState();
  state.selectedPalette.value = 'ocean';
  const references = {
    rules: [...state.manualSpecificRules],
    palette: state.appliedPaletteColors.value,
    results: state.results.value,
    legend: state.legendEntries.value,
    ledger: state.validatedStyleFingerprintByResultKey.value,
    selectedPalette: state.selectedPalette.value,
    currentColors: state.currentColors.value,
    pendingPaletteColors: state.pendingPaletteColors.value,
    legendBoxSize: state.adv.legend_box_size,
    legendFontSize: state.adv.legend_font_size
  };
  const harness = createHarness({ state });
  const initialMountedPayload = harness.mounted.svg.payload;
  const { command } = await buildCommand({
    state,
    after,
    harness,
    extra: {
      presentationPatch: { legendBoxSize: 12, legendFontSize: 12 },
      reconcile: () => { throw new Error('reconcile failed'); }
    }
  });
  await assert.rejects(() => command.apply(), /reconcile failed/);
  assert.deepEqual(state.manualSpecificRules, references.rules);
  assert.equal(state.appliedPaletteColors.value, references.palette);
  assert.equal(state.results.value, references.results);
  assert.equal(state.legendEntries.value, references.legend);
  assert.equal(state.validatedStyleFingerprintByResultKey.value, references.ledger);
  assert.equal(state.selectedPalette.value, references.selectedPalette);
  assert.equal(state.currentColors.value, references.currentColors);
  assert.equal(state.pendingPaletteColors.value, references.pendingPaletteColors);
  assert.equal(state.adv.legend_box_size, references.legendBoxSize);
  assert.equal(state.adv.legend_font_size, references.legendFontSize);
  assert.equal(state.semanticStyleRevision.value, 7);
  assert.equal(state.semanticStyleFingerprint.value, styleFingerprint(baseStyle));
  assert.equal(harness.mounted.svg.payload, initialMountedPayload);
  assert.equal(harness.calls.restore, 1);
});

test('epoch, generation, revision, fingerprint, catalogue, Result, and mount guards reject stale apply', async (t) => {
  const after = replaceFeatureBulkStyleSnapshot(baseStyle, {
    rules: [{ ...kinaseRule, color: '#abcdef' }]
  });
  const cases = [{
    name: 'epoch', mutate: ({ state }) => { state.documentEpoch.value += 1; }
  }, {
    name: 'generation', mutate: ({ state }) => { state.resultGenerationKey.value += 1; }
  }, {
    name: 'revision', mutate: ({ state }) => { state.semanticStyleRevision.value += 1; }
  }, {
    name: 'fingerprint', mutate: ({ state }) => { state.semanticStyleFingerprint.value = 'sf1_stale'; }
  }, {
    name: 'catalogue', mutate: () => { catalog.items[0].biologicalFeatures[0].qualifiers.product = ['Changed']; },
    restore: () => { catalog.items[0].biologicalFeatures[0].qualifiers.product = ['Kinase']; }
  }, {
    name: 'Result', mutate: ({ state }) => { state.results.value[1].content += ':external'; }
  }, {
    name: 'mount', mutate: ({ harness }) => {
      harness.mounted.resultIndex = 2;
      harness.mounted.resultKey = 'result-c';
      harness.mounted.svg = fakeRoot('result-c');
    }
  }, {
    name: 'palette UI draft', mutate: ({ state }) => {
      state.currentColors.value = { ...state.currentColors.value, CDS: '#123456' };
    }
  }];
  for (const scenario of cases) {
    await t.test(scenario.name, async () => {
      const state = createState();
      const harness = createHarness({ state });
      const { command } = await buildCommand({ state, after, harness });
      scenario.mutate({ state, harness });
      try {
        assert.equal(await command.apply(), false);
        assert.equal(state.manualSpecificRules[0].color, '#111111');
        assert.equal(state.semanticStyleRevision.value, scenario.name === 'revision' ? 8 : 7);
      } finally {
        scenario.restore?.();
      }
    });
  }
});

test('History apply, switched-mount Undo, and Redo rerun the same all-Result projection', async () => {
  const after = replaceFeatureBulkStyleSnapshot(baseStyle, {
    rules: [{ ...kinaseRule, color: '#abcdef' }]
  });
  const state = createState();
  const harness = createHarness({ state });
  const history = createHistoryManager({
    buildIntent: () => ({
      rules: clone(state.manualSpecificRules),
      results: clone(state.results.value),
      revision: state.semanticStyleRevision.value
    }),
    applyIntent: () => {},
    buildCheckpoint: () => ({}),
    applyCheckpoint: () => {}
  });
  await history.captureBaseline('bulk style baseline');
  const { command } = await buildCommand({ state, after, harness });
  assert.equal(await history.runUndoableCommand('Bulk style', () => command), true);
  assert.equal(history.getUndoCount(), 1);
  harness.mounted.resultIndex = 1;
  harness.mounted.resultKey = 'result-b';
  harness.mounted.svg = fakeRoot('result-b', 'switched');
  state.selectedResultIndex.value = 1;
  assert.equal(await history.undo(), true);
  assert.equal(state.manualSpecificRules[0].color, '#111111');
  assert.equal(state.selectedResultIndex.value, 1, 'command preserves the switched selection');
  assert.equal(state.semanticStyleRevision.value, 9);
  assert.equal(await history.redo(), true);
  assert.equal(state.manualSpecificRules[0].color, '#abcdef');
  assert.equal(state.semanticStyleRevision.value, 10);
  assert.deepEqual(harness.calls.result, [
    ['result-a', 'result-b'],
    ['result-a', 'result-b'],
    ['result-a', 'result-b']
  ]);
});

test('palette UI refs and Bakta legend settings commit, undo, and redo with the semantic snapshot', async () => {
  const after = replaceFeatureBulkStyleSnapshot(baseStyle, {
    appliedPaletteName: 'ocean',
    appliedPaletteColors: { CDS: '#116699', tRNA: '#44aacc' }
  });
  const state = createState();
  state.selectedPalette.value = 'ocean';
  const originalCurrentColors = state.currentColors.value;
  const originalPendingColors = state.pendingPaletteColors.value;
  const { command } = await buildCommand({
    state,
    after,
    extra: {
      presentationPatch: { legendBoxSize: 12, legendFontSize: 12 }
    }
  });

  assert.equal(await command.apply(), true);
  assert.equal(state.appliedPaletteName.value, 'ocean');
  assert.equal(state.selectedPalette.value, 'ocean');
  assert.deepEqual(state.currentColors.value, after.appliedPaletteColors);
  assert.equal(state.pendingPaletteName.value, '');
  assert.deepEqual(state.pendingPaletteColors.value, {});
  assert.equal(state.adv.legend_box_size, 12);
  assert.equal(state.adv.legend_font_size, 12);

  assert.equal(await command.revert(), true);
  assert.equal(state.appliedPaletteName.value, 'default');
  assert.equal(state.selectedPalette.value, 'default');
  assert.deepEqual(state.currentColors.value, baseStyle.appliedPaletteColors);
  assert.equal(state.pendingPaletteName.value, '');
  assert.deepEqual(state.pendingPaletteColors.value, {});
  assert.equal(state.adv.legend_box_size, 10);
  assert.equal(state.adv.legend_font_size, 10);

  assert.equal(await command.apply({ direction: 'redo' }), true);
  assert.equal(state.selectedPalette.value, 'ocean');
  assert.deepEqual(state.currentColors.value, after.appliedPaletteColors);
  assert.equal(state.adv.legend_box_size, 12);
  assert.notEqual(state.currentColors.value, originalCurrentColors);
  assert.notEqual(state.pendingPaletteColors.value, originalPendingColors);
});

test('undo restores a queued palette draft that existed before explicit acceptance', async () => {
  const after = replaceFeatureBulkStyleSnapshot(baseStyle, {
    appliedPaletteName: 'ocean',
    appliedPaletteColors: { CDS: '#116699', tRNA: '#44aacc' }
  });
  const state = createState();
  state.selectedPalette.value = 'ocean';
  state.currentColors.value = clone(after.appliedPaletteColors);
  state.pendingPaletteName.value = 'ocean';
  state.pendingPaletteColors.value = clone(after.appliedPaletteColors);
  const pendingColors = state.pendingPaletteColors.value;
  const { command } = await buildCommand({ state, after });

  assert.equal(await command.apply(), true);
  assert.equal(state.pendingPaletteName.value, '');
  assert.deepEqual(state.pendingPaletteColors.value, {});
  assert.equal(await command.revert(), true);
  assert.equal(state.appliedPaletteName.value, 'default');
  assert.equal(state.selectedPalette.value, 'ocean');
  assert.deepEqual(state.currentColors.value, after.appliedPaletteColors);
  assert.equal(state.pendingPaletteName.value, 'ocean');
  assert.deepEqual(state.pendingPaletteColors.value, after.appliedPaletteColors);
  assert.notEqual(state.pendingPaletteColors.value, pendingColors);
});

test('History rejects an oversized bulk command before mutation', async () => {
  const after = replaceFeatureBulkStyleSnapshot(baseStyle, {
    rules: [{ ...kinaseRule, color: '#abcdef' }]
  });
  const state = createState();
  const results = state.results.value;
  const { command } = await buildCommand({ state, after });
  const history = createHistoryManager({
    buildIntent: () => ({ revision: state.semanticStyleRevision.value }),
    applyIntent: () => {},
    buildCheckpoint: () => ({}),
    applyCheckpoint: () => {},
    maxBytes: command.estimateBytes() - 1
  });
  await history.captureBaseline('oversize baseline');
  const warn = console.warn;
  console.warn = () => {};
  try {
    assert.equal(await history.runUndoableCommand('Oversized bulk style', () => command), false);
  } finally {
    console.warn = warn;
  }
  assert.equal(state.results.value, results);
  assert.equal(state.manualSpecificRules[0].color, '#111111');
  assert.equal(state.semanticStyleRevision.value, 7);
  assert.equal(history.getUndoCount(), 0);
});

test('History post-apply capture failure uses command compensation and preserves stacks', async () => {
  const after = replaceFeatureBulkStyleSnapshot(baseStyle, {
    rules: [{ ...kinaseRule, color: '#abcdef' }]
  });
  const state = createState();
  const before = {
    rules: clone(state.manualSpecificRules),
    results: state.results.value,
    ledger: state.validatedStyleFingerprintByResultKey.value,
    revision: state.semanticStyleRevision.value,
    fingerprint: state.semanticStyleFingerprint.value
  };
  let failCapture = false;
  const history = createHistoryManager({
    buildIntent: () => {
      if (failCapture) {
        failCapture = false;
        throw new Error('forced bulk post-apply capture failure');
      }
      return { revision: state.semanticStyleRevision.value };
    },
    applyIntent: () => {},
    buildCheckpoint: () => ({}),
    applyCheckpoint: () => {}
  });
  await history.captureBaseline('bulk compensation baseline');
  const { command } = await buildCommand({ state, after });
  failCapture = true;
  await assert.rejects(
    history.runUndoableCommand('Bulk style', () => command),
    /forced bulk post-apply capture failure/
  );
  assert.deepEqual(state.manualSpecificRules, before.rules);
  assert.equal(state.results.value, before.results);
  assert.equal(state.validatedStyleFingerprintByResultKey.value, before.ledger);
  assert.equal(state.semanticStyleRevision.value, before.revision);
  assert.equal(state.semanticStyleFingerprint.value, before.fingerprint);
  assert.equal(history.getUndoCount(), 0);
  assert.equal(history.getRedoCount(), 0);
});

test('bulk command has no rendering-runtime dependency', async () => {
  const source = await import('node:fs/promises').then((fs) => fs.readFile(
    new URL('../../gbdraw/web/js/app/feature-editor/bulk-style-command.js', import.meta.url),
    'utf8'
  ));
  assert.doesNotMatch(source, /pyodide|diagram-generation-worker/i);
  assert.match(source, /buildStyleSnapshotCommand/);
});
