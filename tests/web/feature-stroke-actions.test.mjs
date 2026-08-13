import assert from 'node:assert/strict';
import test from 'node:test';

import { createFeatureStrokeActions } from '../../gbdraw/web/js/app/feature-editor/stroke-actions.js';
import { stableFeatureTargetKey } from '../../gbdraw/web/js/app/feature-editor/stroke-scope-plan.js';

const ref = (value) => ({ value });
const clone = (value) => structuredClone(value);

const hashes = {
  a: 'fi1_aaaaaaaaaaaaaaaaaaaaaaaaaa',
  b: 'fi1_bbbbbbbbbbbbbbbbbbbbbbbbbb',
  c: 'fi1_cccccccccccccccccccccccccc'
};

const feature = ({
  resultKey,
  resultIndex,
  resultName,
  recordKey,
  biologicalFeatureId,
  instanceHash,
  svgId,
  caption
}) => ({
  resultKey,
  resultIndex,
  resultName,
  recordKey,
  biologicalFeatureId,
  instanceHash,
  type: 'CDS',
  svg_id: svgId,
  rendered_feature_svg_id: svgId,
  caption
});

const features = [
  feature({
    resultKey: 'result-a', resultIndex: 0, resultName: 'a.svg',
    recordKey: 'record-a', biologicalFeatureId: 'feature-a', instanceHash: hashes.a,
    svgId: 'svg-a', caption: 'Core'
  }),
  feature({
    resultKey: 'result-b', resultIndex: 1, resultName: 'b.svg',
    recordKey: 'record-b', biologicalFeatureId: 'feature-b', instanceHash: hashes.b,
    svgId: 'svg-b', caption: 'Core'
  }),
  feature({
    resultKey: 'result-b', resultIndex: 1, resultName: 'b.svg',
    recordKey: 'record-b', biologicalFeatureId: 'feature-c', instanceHash: hashes.c,
    svgId: 'svg-c', caption: 'Accessory'
  })
];

const catalog = {
  schema: 3,
  items: [{
    resultKey: 'result-a', resultIndex: 0, resultName: 'a.svg', recordKeys: ['record-a'],
    biologicalFeatures: [{
      recordKey: 'record-a', biologicalFeatureId: 'feature-a', instanceHash: hashes.a, type: 'CDS'
    }],
    features: [{ recordKey: 'record-a', biologicalFeatureId: 'feature-a', svgId: 'svg-a' }]
  }, {
    resultKey: 'result-b', resultIndex: 1, resultName: 'b.svg', recordKeys: ['record-b'],
    biologicalFeatures: [{
      recordKey: 'record-b', biologicalFeatureId: 'feature-b', instanceHash: hashes.b, type: 'CDS'
    }, {
      recordKey: 'record-b', biologicalFeatureId: 'feature-c', instanceHash: hashes.c, type: 'CDS'
    }],
    features: [{ recordKey: 'record-b', biologicalFeatureId: 'feature-b', svgId: 'svg-b' }, {
      recordKey: 'record-b', biologicalFeatureId: 'feature-c', svgId: 'svg-c'
    }]
  }]
};

class FakeNode {
  constructor(value) { this.value = value; }
  cloneNode() { return new FakeNode(this.value); }
}

class FakeSvg {
  constructor(marker) {
    this.attributeMap = new Map([['data-marker', marker]]);
    this.childNodes = [new FakeNode(marker)];
  }

  get attributes() {
    return [...this.attributeMap].map(([name, value]) => ({ name, value }));
  }

  getAttribute(name) { return this.attributeMap.get(name) ?? null; }
  setAttribute(name, value) { this.attributeMap.set(name, String(value)); }
  removeAttribute(name) { this.attributeMap.delete(name); }
  replaceChildren(...children) { this.childNodes = children; }
  cloneNode(deep = false) {
    const copy = new FakeSvg(this.getAttribute('data-marker'));
    copy.attributeMap = new Map(this.attributeMap);
    copy.childNodes = deep ? this.childNodes.map((node) => node.cloneNode(true)) : [];
    return copy;
  }
}

const createHistory = () => {
  const calls = [];
  return {
    calls,
    async runUndoableCommand(label, build) {
      const command = await build();
      calls.push({ label, command });
      if (!command || command.noop) return false;
      return (await command.apply()) !== false;
    }
  };
};

const createHarness = ({
  commandBuilder = async () => ({
    apply: () => true,
    revert: () => true,
    estimateBytes: () => 1
  }),
  resultProjection = () => ({ nextResults: [] }),
  mountedSvg = null,
  exactScopeAvailable = true,
  previewRuntime = null
} = {}) => {
  const state = {
    featureCatalog: ref(catalog),
    extractedFeatures: ref(clone(features)),
    selectedFeatures: ref([]),
    featureStrokeOverrides: {},
    legendStrokeOverrides: {},
    originalSvgStroke: ref({ color: '#000000', width: 0.5 }),
    featureExactScopeAvailable: ref(exactScopeAvailable),
    manualSpecificRules: [],
    appliedPaletteName: ref('default'),
    appliedPaletteColors: ref({ CDS: '#cccccc' }),
    resultGenerationKey: ref(7),
    documentEpoch: ref(3),
    semanticStyleRevision: ref(4),
    semanticStyleFingerprint: ref(''),
    results: ref([{ name: 'a.svg', content: 'a' }, { name: 'b.svg', content: 'b' }]),
    selectedResultIndex: ref(0),
    svgContainer: ref(mountedSvg ? { querySelector: () => mountedSvg } : null),
    clickedFeature: ref(null),
    legendEntries: ref([])
  };
  const history = createHistory();
  const pendingFeatureStrokePlan = ref(null);
  const featureStrokePlanStatus = ref('idle');
  const featureStrokePlanProgress = ref(null);
  const actions = createFeatureStrokeActions({
    state,
    history,
    previewRuntime,
    pendingFeatureStrokePlan,
    featureStrokePlanStatus,
    featureStrokePlanProgress,
    getEffectiveLegendCaption: (entry) => entry.caption,
    commandBuilder,
    resultProjection
  });
  return {
    actions,
    state,
    history,
    pendingFeatureStrokePlan,
    featureStrokePlanStatus,
    featureStrokePlanProgress
  };
};

test('public facade contains only view-model and request/resolve/cancel operations', () => {
  const { actions } = createHarness();
  assert.deepEqual(Object.keys(actions).sort(), [
    'cancelFeatureStrokeScope',
    'getFeatureStrokeViewModel',
    'requestFeatureStrokeChange',
    'requestSelectedFeatureStrokeChange',
    'resolveFeatureStrokeScope'
  ]);
});

test('view model resolves exact, caption, and inherited SVG stroke without locking the picker', () => {
  const { actions, state } = createHarness();
  state.legendStrokeOverrides.Core = { strokeColor: '#334455', strokeWidth: 1.25 };
  let model = actions.getFeatureStrokeViewModel(features[0]);
  assert.equal(model.effectiveColor, '#334455');
  assert.equal(model.effectiveWidth, 1.25);
  assert.equal(model.origin, 'legend-caption');
  assert.equal(model.pickerEnabled, true);
  assert.equal(model.allowNone, false);

  state.featureStrokeOverrides['record-a\u0000feature-a'] = {
    strokeColor: '#abcdef', strokeWidth: 2
  };
  model = actions.getFeatureStrokeViewModel(features[0]);
  assert.equal(model.effectiveColor, '#abcdef');
  assert.equal(model.origin, 'local');
  assert.equal(model.canReset, true);
});

test('group request remains mutation-free until resolving its planned candidate', async () => {
  let commandArguments = null;
  const resultProjection = () => ({ nextResults: [] });
  const harness = createHarness({
    resultProjection,
    commandBuilder: async (input) => {
      commandArguments = input;
      return { apply: () => true, revert: () => true, estimateBytes: () => 1 };
    }
  });
  const requested = await harness.actions.requestFeatureStrokeChange(
    features[0],
    { kind: 'color', color: '#abcdef', strokeWidth: 2 }
  );
  assert.equal(requested, false);
  assert.equal(harness.featureStrokePlanStatus.value, 'needs-scope');
  assert.equal(harness.history.calls.length, 0);
  const pending = harness.pendingFeatureStrokePlan.value;
  assert.deepEqual(pending.candidates.map((candidate) => candidate.semanticScope), [
    'legend-caption', 'single'
  ]);
  const group = pending.candidates[0];
  assert.equal(group.affectedResultCount, 2);

  assert.equal(
    await harness.actions.resolveFeatureStrokeScope(pending.token, group.id),
    true
  );
  assert.equal(harness.history.calls.length, 1);
  assert.equal(commandArguments.plan.semanticScope, 'legend-caption');
  assert.deepEqual(commandArguments.plan.affectedResultKeys, ['result-a', 'result-b']);
  assert.equal(commandArguments.catalog, catalog);
  assert.equal(commandArguments.prepareResultProjection, resultProjection);
  assert.equal(typeof commandArguments.getMountedContext, 'function');
  assert.equal(typeof commandArguments.captureRuntimeState, 'function');
  assert.equal(typeof commandArguments.restoreRuntimeState, 'function');
  assert.equal(typeof commandArguments.commitMountedProjection, 'function');
  assert.equal(typeof commandArguments.restoreMountedProjection, 'function');
  assert.equal(harness.featureStrokePlanStatus.value, 'committed');
  assert.equal(harness.pendingFeatureStrokePlan.value, null);
});

test('selection requests one exact multi-Result command and one-target selection uses single', async () => {
  const plans = [];
  const harness = createHarness({
    commandBuilder: async ({ plan }) => {
      plans.push(plan);
      return { apply: () => true, revert: () => true, estimateBytes: () => 1 };
    }
  });
  assert.equal(await harness.actions.requestSelectedFeatureStrokeChange(
    [features[0], features[2]],
    { kind: 'stroke', strokeColor: '#123456', strokeWidth: 1.5 }
  ), true);
  assert.equal(plans[0].semanticScope, 'selected-features');
  assert.deepEqual(plans[0].affectedResultKeys, ['result-a', 'result-b']);
  assert.deepEqual(
    plans[0].targetsByResult.flatMap((group) => group.featureKeys).sort(),
    [stableFeatureTargetKey(features[0]), stableFeatureTargetKey(features[2])].sort()
  );

  assert.equal(await harness.actions.requestSelectedFeatureStrokeChange(
    [features[2]],
    { kind: 'inherit' }
  ), true);
  assert.equal(plans[1].semanticScope, 'single');
  assert.deepEqual(plans[1].affectedResultKeys, ['result-b']);
  assert.equal(harness.history.calls.length, 2);
});

test('cancel validates the token and supersedes an asynchronous command preparation', async () => {
  let release;
  let entered;
  let applied = 0;
  const waiting = new Promise((resolve) => { release = resolve; });
  const started = new Promise((resolve) => { entered = resolve; });
  const harness = createHarness({
    commandBuilder: async () => {
      entered();
      await waiting;
      return {
        apply: () => { applied += 1; return true; },
        revert: () => true,
        estimateBytes: () => 1
      };
    }
  });
  const request = harness.actions.requestFeatureStrokeChange(
    features[2],
    { kind: 'color', color: '#abcdef' },
    { semanticScope: 'single' }
  );
  await started;
  const token = harness.pendingFeatureStrokePlan.value.token;
  assert.equal(harness.actions.cancelFeatureStrokeScope('wrong-token'), false);
  assert.equal(harness.actions.cancelFeatureStrokeScope(token), true);
  release();
  assert.equal(await request, false);
  assert.equal(applied, 0);
  assert.equal(harness.featureStrokePlanStatus.value, 'idle');
  assert.equal(harness.pendingFeatureStrokePlan.value, null);
});

test('mounted adapters capture, commit, remount, and restore the exact root/runtime state', async () => {
  const mounted = new FakeSvg('before');
  const runtime = { dirty: true, dirtyReasons: new Set(['existing']) };
  const events = [];
  const previewRuntime = {
    getActiveRuntime: () => runtime,
    clearActiveRuntime: () => events.push('clear'),
    mountResultSvg: (index, svg) => {
      events.push(`mount:${index}:${svg.getAttribute('data-marker')}`);
      return runtime;
    },
    invalidatePreviewIndexes: (reason) => events.push(`invalidate:${reason}`)
  };
  let adapterAssertions = 0;
  const harness = createHarness({
    mountedSvg: mounted,
    previewRuntime,
    commandBuilder: async (input) => {
      const snapshot = input.captureRuntimeState();
      assert.equal(snapshot.resultKey, 'result-a');
      assert.equal(snapshot.dirty, true);
      input.commitMountedProjection({
        prepared: { preparedMountedSvg: new FakeSvg('after') }
      });
      assert.equal(mounted.getAttribute('data-marker'), 'after');
      input.reconcile({
        prepared: { projection: { preparedMountedSvg: new FakeSvg('after') } }
      });
      input.restoreMountedProjection();
      assert.equal(input.restoreRuntimeState(snapshot), true);
      assert.equal(mounted.getAttribute('data-marker'), 'before');
      assert.equal(runtime.dirty, true);
      assert.deepEqual([...runtime.dirtyReasons], ['existing']);
      adapterAssertions += 1;
      return { apply: () => true, revert: () => true, estimateBytes: () => 1 };
    }
  });
  assert.equal(await harness.actions.requestFeatureStrokeChange(
    features[2],
    { kind: 'color', color: '#abcdef' },
    { semanticScope: 'single' }
  ), true);
  assert.equal(adapterAssertions, 1);
  assert.ok(events.includes('invalidate:feature-stroke-commit'));
  assert.ok(events.includes('invalidate:feature-stroke-rollback'));
});

test('unsafe legacy exact request fails closed without entering History', async () => {
  const harness = createHarness({ exactScopeAvailable: false });
  const changed = await harness.actions.requestFeatureStrokeChange(
    features[2],
    { kind: 'color', color: '#abcdef' },
    { semanticScope: 'single' }
  );
  assert.equal(changed, false);
  assert.equal(harness.featureStrokePlanStatus.value, 'invalid');
  assert.match(harness.featureStrokePlanProgress.value.message, /Generate to enable exact/);
  assert.equal(harness.history.calls.length, 0);
});
