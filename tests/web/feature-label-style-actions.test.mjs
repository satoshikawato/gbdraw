import assert from 'node:assert/strict';
import test from 'node:test';

import { createFeatureLabelStyleActions } from '../../gbdraw/web/js/app/feature-editor/label-style-actions.js';

const ref = (value) => ({ value });
const hashes = {
  a: 'fi1_aaaaaaaaaaaaaaaaaaaaaaaaaa',
  b: 'fi1_bbbbbbbbbbbbbbbbbbbbbbbbbb'
};

const catalog = {
  schema: 3,
  items: [{
    resultKey: 'result-a', resultIndex: 0, resultName: 'a.svg', recordKeys: ['record-a'],
    biologicalFeatures: [{
      recordKey: 'record-a', biologicalFeatureId: 'feature-a', instanceHash: hashes.a,
      type: 'CDS', product: 'Kinase', qualifiers: { product: ['Kinase'] }
    }],
    features: [{ recordKey: 'record-a', biologicalFeatureId: 'feature-a', svgId: 'svg-a' }]
  }, {
    resultKey: 'result-b', resultIndex: 1, resultName: 'b.svg', recordKeys: ['record-b'],
    biologicalFeatures: [{
      recordKey: 'record-b', biologicalFeatureId: 'feature-b', instanceHash: hashes.b,
      type: 'CDS', product: 'Kinase', qualifiers: { product: ['Kinase'] }
    }],
    features: [{ recordKey: 'record-b', biologicalFeatureId: 'feature-b', svgId: 'svg-b' }]
  }]
};

const target = {
  ...catalog.items[0].biologicalFeatures[0],
  resultKey: 'result-a', resultIndex: 0, resultName: 'a.svg',
  renderedFeatureSvgId: 'svg-a', svg_id: 'svg-a'
};

const createState = () => ({
  manualSpecificRules: [],
  appliedPaletteName: ref('default'),
  appliedPaletteColors: ref({ CDS: '#cccccc' }),
  featureCatalog: ref(catalog),
  biologicalFeatures: ref(catalog.items.flatMap((item, resultIndex) => (
    item.biologicalFeatures.map((feature) => ({
      ...feature,
      resultKey: item.resultKey,
      resultIndex,
      resultName: item.resultName,
      renderedFeatureSvgId: item.features[0].svgId
    }))
  ))),
  extractedFeatures: ref([target, {
    ...catalog.items[1].biologicalFeatures[0],
    resultKey: 'result-b', resultIndex: 1, resultName: 'b.svg',
    renderedFeatureSvgId: 'svg-b', svg_id: 'svg-b'
  }]),
  editableLabels: ref([{
    key: 'label-a', featureId: 'svg-a', text: 'Kinase', sourceText: 'Kinase'
  }]),
  labelTextFeatureOverrides: {},
  labelTextBulkOverrides: {},
  labelTextFeatureOverrideSources: {},
  selectedFeatures: ref([]),
  featureExactScopeAvailable: ref(true),
  resultGenerationKey: ref(7),
  documentEpoch: ref(3),
  semanticStyleRevision: ref(4),
  selectedResultIndex: ref(0),
  results: ref([{ name: 'a.svg', content: 'a' }, { name: 'b.svg', content: 'b' }]),
  svgContainer: ref(null),
  clickedFeature: ref({ feat: target, labelText: 'Kinase' })
});

test('facade keeps scope choice mutation-free and sends one resolved command to History', async () => {
  const state = createState();
  const pending = ref(null);
  const status = ref('idle');
  const progress = ref(null);
  let historyCalls = 0;
  let commandBuilds = 0;
  let applies = 0;
  const history = {
    runUndoableCommand: async (label, build) => {
      historyCalls += 1;
      assert.equal(label, 'Change feature label');
      const command = await build();
      return command.apply();
    }
  };
  const actions = createFeatureLabelStyleActions({
    state,
    history,
    pendingFeatureLabelPlan: pending,
    featureLabelPlanStatus: status,
    featureLabelPlanProgress: progress,
    commandBuilder: async ({ plan }) => {
      commandBuilds += 1;
      assert.equal(plan.semanticScope, 'source-annotation-label');
      assert.deepEqual(plan.affectedResultKeys, ['result-a', 'result-b']);
      return {
        apply: async () => { applies += 1; return true; },
        revert: async () => true,
        snapshot: () => null,
        compensate: () => true,
        estimateBytes: () => 1
      };
    }
  });

  assert.equal(await actions.requestFeatureLabelTextChange(target, 'Renamed'), false);
  assert.equal(status.value, 'needs-scope');
  assert.equal(historyCalls, 0);
  assert.deepEqual(state.labelTextBulkOverrides, {});
  const sourceCandidate = pending.value.candidates.find((candidate) => (
    candidate.semanticScope === 'source-annotation-label'
  ));
  assert.ok(sourceCandidate);
  assert.equal(await actions.resolveFeatureLabelScope(pending.value.token, sourceCandidate.id), true);
  assert.equal(historyCalls, 1);
  assert.equal(commandBuilds, 1);
  assert.equal(applies, 1);
  assert.equal(status.value, 'committed');
  assert.equal(pending.value, null);
});

test('cancel invalidates a pending token without History or semantic mutation', async () => {
  const state = createState();
  const pending = ref(null);
  const status = ref('idle');
  const progress = ref(null);
  let historyCalls = 0;
  const actions = createFeatureLabelStyleActions({
    state,
    history: { runUndoableCommand: async () => { historyCalls += 1; return true; } },
    pendingFeatureLabelPlan: pending,
    featureLabelPlanStatus: status,
    featureLabelPlanProgress: progress
  });
  await actions.requestFeatureLabelTextChange(target, 'Cancelled');
  const token = pending.value.token;
  assert.equal(actions.cancelFeatureLabelScope(token), true);
  assert.equal(pending.value, null);
  assert.equal(status.value, 'idle');
  assert.equal(historyCalls, 0);
  assert.deepEqual(state.labelTextFeatureOverrides, {});
  assert.deepEqual(state.labelTextBulkOverrides, {});
  assert.equal(actions.resolveFeatureLabelScope(token, 'missing'), false);
});

test('view-model owner exposes effective inherited and exact label state', () => {
  const state = createState();
  const refs = { pending: ref(null), status: ref('idle'), progress: ref(null) };
  const actions = createFeatureLabelStyleActions({
    state,
    history: { runUndoableCommand: async () => false },
    pendingFeatureLabelPlan: refs.pending,
    featureLabelPlanStatus: refs.status,
    featureLabelPlanProgress: refs.progress
  });
  const inherited = actions.getFeatureLabelViewModel(target);
  assert.equal(inherited.effectiveText, 'Kinase');
  assert.equal(inherited.origin, 'source-annotation');
  state.labelTextFeatureOverrides['svg-a'] = 'Exact';
  state.labelTextFeatureOverrideSources['svg-a'] = 'Kinase';
  const exact = actions.getFeatureLabelViewModel(target);
  assert.equal(exact.effectiveText, 'Exact');
  assert.equal(exact.origin, 'feature-override');
});
