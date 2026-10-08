import assert from 'node:assert/strict';
import { createFeatureVisibilityActions } from '../../gbdraw/web/js/app/feature-editor/visibility-actions.js';
import {
  requestFeatureVisibilityRules,
  setFeatureVisibilityOverride
} from '../../gbdraw/web/js/services/feature-visibility.js';
import { createRulePreparation } from '../../gbdraw/web/js/app/rule-matching.js';
import { resultCatalogFeatures } from '../../gbdraw/web/js/services/feature-catalog.js';
import { evaluatePythonRules } from './helpers/python-rule-evaluator.mjs';
import { withDrawings } from './helpers/drawing-state.mjs';

const ref = (value) => ({ value });
// The live projection reads Python's rule matches (R4) and the committed
// request's feature types.
const rulePreparationFor = (state) => createRulePreparation({
  state: withDrawings(state),
  evaluate: evaluatePythonRules,
  visibilityRules: () => requestFeatureVisibilityRules(state.featureVisibilityManualRules)
});
const committedRequest = (selectedFeaturesSet) => () => ({ diagramOptions: { selectedFeaturesSet } });
// Per-feature visibility is the identity row of the feature (design Q4).
const identity = (id) => ({ scope: 'circular', record_key: 'record-1', biological_feature_id: `bio-${id}` });
const modes = (overrides) => Object.fromEntries(Object.values(overrides)
  .filter((row) => row.featureVisibility !== null)
  .map((row) => [row.biologicalFeatureId.replace(/^bio-/, ''), row.featureVisibility]));
const clear = (overrides) => Object.keys(overrides).forEach((key) => delete overrides[key]);
const featureA = { svg_id: 'feature-a', type: 'CDS', label: 'A', ...identity('feature-a') };
const featureB = { svg_id: 'feature-b', type: 'CDS', label: 'B', ...identity('feature-b') };
const extractedFeatures = ref([featureA, featureB]);
const orthogroups = ref([]);
const featureVisibilityOverrides = {};
const clickedFeature = ref({ svg_id: 'feature-a', featureVisibility: 'default' });
const featureVisibilityScopeDialog = {};
const selectedResultIndex = ref(0);
const resultGenerationKey = ref('generation-1');
const appliedPreviewChanges = [];
// F-3: each visibility edit hands the feature's label to the label owner,
// through the port the composition root registers (R13).
const labelVisibilityCalls = [];

const actionState = {
  clickedFeature,
  extractedFeatures,
  orthogroups,
  featureVisibilityManualRules: [],
  featureVisibilityRules: ref([]),
  featureOverrides: featureVisibilityOverrides,
  featureVisibilityScopeDialog,
  resultGenerationKey,
  results: ref([{ name: 'one.svg', content: '<svg></svg>' }]),
  selectedResultIndex,
  svgContainer: ref({
    querySelector: (selector) => (selector === 'svg' ? {} : null)
  })
};
const actions = createFeatureVisibilityActions({
  state: withDrawings(actionState),
  rulePreparation: rulePreparationFor(actionState),
  getCommittedRequest: committedRequest(['CDS']),
  applyVisibilityPreviewChanges: (changes, options = {}) => {
    appliedPreviewChanges.push({ changes, reason: options.reason });
    return true;
  },
  ports: {
    applyFeatureVisibilityToLabels: (options = {}) => {
      labelVisibilityCalls.push(options.reflow !== false);
      return true;
    }
  },
  selectResult: (index) => {
    selectedResultIndex.value = index;
    return true;
  }
});

const command = actions.buildSelectedFeaturesVisibilityCommand([featureA, featureB], 'off');
assert.ok(command);
assert.equal(await command.apply(), true);
assert.deepEqual(modes(featureVisibilityOverrides), {
  'feature-a': 'off',
  'feature-b': 'off'
});
assert.equal(appliedPreviewChanges.length, 1);
assert.deepEqual(
  appliedPreviewChanges[0].changes.map((change) => [change.featureId, change.mode]),
  [['feature-a', 'off'], ['feature-b', 'off']]
);

assert.equal(await command.revert(), true);
assert.deepEqual(featureVisibilityOverrides, {});
assert.equal(appliedPreviewChanges.length, 2);
assert.deepEqual(labelVisibilityCalls, [
  true,
  true
]);
assert.deepEqual(
  appliedPreviewChanges[1].changes.map((change) => [change.featureId, change.mode]),
  [['feature-a', 'on'], ['feature-b', 'on']]
);

assert.equal(actions.setFeatureVisibility(featureA, 'off', {
  triggerReflow: false,
  scope: { id: 'feature' }
}), true);
assert.equal(modes(featureVisibilityOverrides)['feature-a'], 'off');
assert.equal(appliedPreviewChanges.length, 3);
assert.deepEqual(labelVisibilityCalls.at(-1), false,
  'the label follows the feature even when the caller declines the reflow');
assert.deepEqual(
  appliedPreviewChanges[2].changes.map((change) => [change.featureId, change.mode]),
  [['feature-a', 'off']]
);
clear(featureVisibilityOverrides);

const sourceIdFeature = {
  svg_id: 'feature-source-id',
  type: 'CDS',
  proteinId: 'h_aaaaaaaaaaaaaaaaaaaaaaaaaa',
  sourceProteinId: 'WP_012345678.1'
};
clickedFeature.value = {
  svg_id: sourceIdFeature.svg_id,
  featureVisibility: 'default',
  feat: sourceIdFeature
};
await actions.updateClickedFeatureVisibility('off');
assert.equal(featureVisibilityScopeDialog.show, true);
assert.ok(featureVisibilityScopeDialog.scopes.some((scope) => (
  scope.id === 'protein_id' &&
  scope.value === 'WP_012345678.1' &&
  scope.label === 'Exact protein ID: WP_012345678.1'
)));
assert.ok(featureVisibilityScopeDialog.scopes.every((scope) => !scope.label.includes('h_')));

const runtimeOnlyFeature = {
  svg_id: 'feature-runtime-only',
  type: 'CDS',
  proteinId: 'h_bbbbbbbbbbbbbbbbbbbbbbbbbb'
};
clickedFeature.value = {
  svg_id: runtimeOnlyFeature.svg_id,
  featureVisibility: 'default',
  feat: runtimeOnlyFeature
};
featureVisibilityScopeDialog.show = false;
await actions.updateClickedFeatureVisibility('off');
assert.equal(featureVisibilityScopeDialog.show, false);
clear(featureVisibilityOverrides);

const reversedFeature = {
  svg_id: 'display-feature_record_3',
  stable_feature_id: 'source-feature',
  record_idx: 2,
  type: 'CDS',
  orthogroupId: 'og_reverse'
};
const groupedFeature = {
  svg_id: 'grouped-feature_record_1',
  stable_feature_id: 'grouped-source-feature',
  record_idx: 0,
  type: 'CDS',
  orthogroupId: 'og_reverse'
};
extractedFeatures.value.push(reversedFeature, groupedFeature);
orthogroups.value = [{
  id: 'og_reverse',
  members: [
    {
      recordIndex: 2,
      featureSvgId: 'source-feature',
      stableFeatureSvgId: 'source-feature',
      renderedFeatureSvgId: 'display-feature_record_3'
    },
    {
      recordIndex: 0,
      featureSvgId: 'grouped-source-feature',
      stableFeatureSvgId: 'grouped-source-feature',
      renderedFeatureSvgId: 'grouped-feature_record_1'
    }
  ]
}];
clickedFeature.value = {
  svg_id: reversedFeature.svg_id,
  featureVisibility: 'default',
  feat: reversedFeature
};
await actions.updateClickedFeatureVisibility('off');
const reversedGroupScope = featureVisibilityScopeDialog.scopes.find((scope) => scope.id === 'orthogroup');
assert.ok(reversedGroupScope);
assert.deepEqual(
  reversedGroupScope.features.map((feature) => feature.svg_id),
  ['display-feature_record_3', 'grouped-feature_record_1']
);

const strictTrigger = {
  svg_id: 'strict-trigger-rendered',
  stable_feature_id: 'strict-trigger-source',
  record_idx: 2,
  type: 'CDS',
  orthogroupId: 'og_strict'
};
const wrongRecordFeature = {
  svg_id: 'shared-rendered',
  stable_feature_id: 'shared-source',
  record_idx: 0,
  type: 'CDS',
  orthogroupId: 'og_strict'
};
extractedFeatures.value = [strictTrigger, wrongRecordFeature];
orthogroups.value = [{
  id: 'og_strict',
  members: [
    {
      recordIndex: 2,
      featureSvgId: 'strict-trigger-source',
      stableFeatureSvgId: 'strict-trigger-source',
      renderedFeatureSvgId: 'strict-trigger-rendered'
    },
    {
      recordIndex: 1,
      featureSvgId: 'shared-source',
      stableFeatureSvgId: 'shared-source',
      renderedFeatureSvgId: 'shared-rendered'
    }
  ]
}];
clickedFeature.value = {
  svg_id: strictTrigger.svg_id,
  featureVisibility: 'default',
  feat: strictTrigger
};
featureVisibilityScopeDialog.show = false;
featureVisibilityScopeDialog.scopes = [];
await actions.updateClickedFeatureVisibility('off');
assert.equal(featureVisibilityScopeDialog.show, false);
clear(featureVisibilityOverrides);

const duplicateA = {
  svg_id: 'duplicate-rendered',
  stable_feature_id: 'duplicate-source',
  record_idx: 1,
  type: 'CDS',
  orthogroupId: 'og_strict'
};
const duplicateB = { ...duplicateA, label: 'duplicate metadata row' };
extractedFeatures.value = [strictTrigger, duplicateA, duplicateB];
orthogroups.value[0].members[1] = {
  recordIndex: 1,
  featureSvgId: 'duplicate-source',
  stableFeatureSvgId: 'duplicate-source',
  renderedFeatureSvgId: 'duplicate-rendered'
};
featureVisibilityScopeDialog.show = false;
featureVisibilityScopeDialog.scopes = [];
await actions.updateClickedFeatureVisibility('off');
assert.equal(featureVisibilityScopeDialog.show, false);
clear(featureVisibilityOverrides);

extractedFeatures.value = [strictTrigger, duplicateA];
orthogroups.value = [orthogroups.value[0], { ...orthogroups.value[0] }];
featureVisibilityScopeDialog.show = false;
featureVisibilityScopeDialog.scopes = [];
await actions.updateClickedFeatureVisibility('off');
assert.equal(featureVisibilityScopeDialog.show, false);
clear(featureVisibilityOverrides);

const previewChangeCountBeforeStaleApply = appliedPreviewChanges.length;
resultGenerationKey.value = 'generation-2';
assert.equal(await command.apply(), false);
assert.deepEqual(featureVisibilityOverrides, {});
assert.equal(appliedPreviewChanges.length, previewChangeCountBeforeStaleApply);

// FE-04: an Exact product hide is committed as an editor qualifier rule. The
// projection that History runs after Undo/Redo resolves it with the action's
// resolver (Python's matches), so the feature stays hidden. R-2: every rule,
// also one typed in the Features panel or loaded from a TSV, hides live as
// Generate hides.
{
  const nd1 = { svg_id: 'nd1', type: 'CDS', qualifiers: { product: ['NADH dehydrogenase subunit 1'] }, ...identity('nd1') };
  const nd1Case = {
    svg_id: 'nd1-case', type: 'CDS', qualifiers: { product: ['nadh DEHYDROGENASE subunit 1'] }, ...identity('nd1-case')
  };
  const nd2 = { svg_id: 'nd2', type: 'CDS', qualifiers: { product: ['NADH dehydrogenase subunit 2'] }, ...identity('nd2') };
  const geneRna = {
    svg_id: 'nd1-rna', type: 'tRNA', qualifiers: { product: ['NADH dehydrogenase subunit 1'] }, ...identity('nd1-rna')
  };
  const manualRules = [];
  const overrides = {};
  const reconciled = [];
  const scopeDialog = {};
  const productState = {
    clickedFeature: ref({ svg_id: nd1.svg_id, featureVisibility: 'default', feat: nd1 }),
    extractedFeatures: ref([nd1, nd1Case, nd2, geneRna]),
    orthogroups: ref([]),
    featureVisibilityManualRules: manualRules,
    featureVisibilityRules: ref([]),
    featureOverrides: overrides,
    featureVisibilityScopeDialog: scopeDialog,
    resultGenerationKey: ref('generation-1'),
    results: ref([{ name: 'one.svg', content: '<svg></svg>' }]),
    selectedResultIndex: ref(0),
    svgContainer: ref({ querySelector: (selector) => (selector === 'svg' ? {} : null) })
  };
  const productPreparation = rulePreparationFor(productState);
  const productActions = createFeatureVisibilityActions({
    state: withDrawings(productState),
    rulePreparation: productPreparation,
    getCommittedRequest: committedRequest(['CDS', 'tRNA']),
    applyVisibilityPreviewChanges: (changes) => {
      reconciled.push(Object.fromEntries(changes.map((change) => [change.featureId, change.mode])));
      return true;
    },
    ports: { applyFeatureVisibilityToLabels: () => true },
    selectResult: () => true
  });
  await productActions.updateClickedFeatureVisibility('off');
  assert.equal(scopeDialog.show, true);
  assert.equal(await productActions.handleFeatureVisibilityScopeChoice('product'), true);
  const hidden = { nd1: 'off', 'nd1-case': 'off', nd2: 'on', 'nd1-rna': 'on' };
  assert.deepEqual(reconciled.at(-1), hidden);
  assert.equal(manualRules.length, 1);
  assert.equal(await productActions.projectFeatureVisibility(), true);
  assert.deepEqual(reconciled.at(-1), hidden);
  setFeatureVisibilityOverride(overrides, nd1Case, 'on');
  await productActions.projectFeatureVisibility();
  assert.equal(reconciled.at(-1)['nd1-case'], 'on', 'a per-feature override takes precedence');
  clear(overrides);
  manualRules.splice(0, manualRules.length, {
    ...manualRules[0], source: 'file', recordId: '*', featureType: '*', qualifier: 'Product', value: 'subunit 2$'
  });
  assert.equal(await productPreparation.prepareDrawn(), true);
  await productActions.projectFeatureVisibility();
  assert.deepEqual(reconciled.at(-1), { nd1: 'on', 'nd1-case': 'on', nd2: 'off', 'nd1-rna': 'on' },
    'a loaded rule hides live as Generate does');
}

// OV-19 (PD-OI-066, R10): every Features panel rule edit is one transition
// that writes the rules and projects them onto the displayed Result as
// Generate draws them. A table Generate rejects stays in the draft, changes
// nothing, and reports Generate's error until an edit leaves a table Generate
// accepts.
{
  const fl1 = { svg_id: 'fl1', type: 'CDS', qualifiers: { locus_tag: ['FL1'] }, ...identity('fl1') };
  const fl2 = { svg_id: 'fl2', type: 'CDS', qualifiers: { locus_tag: ['FL2'] }, ...identity('fl2') };
  const rules = [];
  const shown = { fl1: 'on', fl2: 'on' };
  const labelProjections = [];
  const panelState = {
    clickedFeature: ref(null),
    extractedFeatures: ref([fl1, fl2]),
    orthogroups: ref([]),
    featureVisibilityManualRules: rules,
    featureVisibilityRules: ref([]),
    featureOverrides: {},
    featureVisibilityScopeDialog: {},
    resultGenerationKey: ref('generation-1'),
    results: ref([{ name: 'one.svg', content: '<svg></svg>' }]),
    selectedResultIndex: ref(0),
    svgContainer: ref({ querySelector: (selector) => (selector === 'svg' ? {} : null) }),
    errorLog: ref(null)
  };
  const panel = createFeatureVisibilityActions({
    state: withDrawings(panelState),
    rulePreparation: rulePreparationFor(panelState),
    getCommittedRequest: committedRequest(['CDS']),
    applyVisibilityPreviewChanges: (changes) => changes.reduce((changed, { featureId, mode }) => {
      if (shown[featureId] === mode) return changed;
      shown[featureId] = mode;
      return true;
    }, false),
    ports: { applyFeatureVisibilityToLabels: (options) => labelProjections.push(options) },
    selectResult: () => true
  });
  const field = (index, name, value) => panel.setFeatureVisibilityRuleField(index, name, value);

  await panel.addFeatureVisibilityRule();
  await field(0, 'qualifier', 'locus_tag');
  assert.deepEqual(shown, { fl1: 'on', fl2: 'on' }, 'a rule without a value is not in the request');
  await field(0, 'value', '^fl1$');
  assert.deepEqual(shown, { fl1: 'off', fl2: 'on' }, 'a rule edit hides live');
  assert.deepEqual(labelProjections, [{ rerender: false }], 'the labels follow the features the edit hides');
  await field(0, 'action', 'show');
  assert.deepEqual(shown, { fl1: 'on', fl2: 'on' });
  await panel.addFeatureVisibilityRule();
  await field(1, 'qualifier', 'locus_tag');
  await field(1, 'value', '^fl');
  assert.deepEqual(shown, { fl1: 'on', fl2: 'off' }, 'the first matching rule decides');
  await panel.moveFeatureVisibilityRuleUp(1);
  assert.deepEqual(rules.map((rule) => rule.value), ['^fl', '^fl1$']);
  assert.deepEqual(shown, { fl1: 'off', fl2: 'off' });

  await field(1, 'value', '^fl1(');
  assert.equal(rules[1].value, '^fl1(', 'the draft keeps a regex Generate rejects');
  assert.deepEqual(shown, { fl1: 'off', fl2: 'off' });
  assert.equal(panelState.errorLog.value?.operation, 'evaluateRules');
  await panel.removeFeatureVisibilityRule(0);
  assert.deepEqual(shown, { fl1: 'off', fl2: 'off' }, 'no edit shows a table Generate rejects');
  assert.notEqual(panelState.errorLog.value, null);
  await field(0, 'value', '^fl1$');
  assert.equal(panelState.errorLog.value, null, 'an edit that Generate accepts clears the report');
  assert.deepEqual(shown, { fl1: 'on', fl2: 'on' });

  assert.equal(await panel.removeFeatureVisibilityRule(5), false, 'an edit of no rule changes nothing');
  assert.equal(await panel.moveFeatureVisibilityRuleDown(0), false);
  assert.equal(rules.length, 1);
}

// R13, R3: the owner reaches the label owner through the one port the
// composition root registers once the label owner exists. Its projection, which
// History apply, the display of a Result, and Load Feature Edits TSV call
// through `projectMountedEditorIntent`, hands the label owner every feature it
// hides or shows (OV-35), asks for the rerender only when a History step or a
// loaded table (`rerender`) must draw a feature the Result does not draw, and
// after a loaded table (`reflow`) also places the labels. A later projection
// supersedes one that still waits for its matches.
{
  const recordFeature = (id, start) => ({
    recordKey: 'REC1', biologicalFeatureId: id, record_id: 'REC1', type: 'CDS', start, end: start + 30, strand: 1,
    anchorProfile: { precision: 'exact', operator: 'single', partOrder: 'biological', strand: '+' },
    qualifiers: { locus_tag: [id] }
  });
  // The Result draws A; B is not drawn.
  const catalog = {
    schema: 5,
    items: [{
      resultIndex: 0,
      resultName: 'result-0.svg',
      recordKeys: ['REC1'],
      biologicalFeatures: [recordFeature('A', 0), recordFeature('B', 100)],
      features: [{
        svgId: 'svg-A', recordKey: 'REC1', biologicalFeatureId: 'A', fillColor: '#000000',
        drawnSelector: { hash: 'svg-A', location: null, recordLocation: null }
      }],
      orthogroups: [],
      annotations: [],
      comparisonMatches: []
    }]
  };
  const overrides = {};
  const mounted = { 'svg-A': 'on' };
  let projections = 0;
  const follows = [];
  const portState = {
    clickedFeature: ref(null),
    extractedFeatures: ref([]),
    orthogroups: ref([]),
    featureVisibilityManualRules: [],
    featureVisibilityRules: ref([]),
    featureOverrides: overrides,
    featureVisibilityScopeDialog: {},
    featureCatalog: ref(catalog),
    generatedMode: ref('circular'),
    resultGenerationKey: ref('generation-1'),
    results: ref([{ name: 'result-0.svg', content: '<svg />' }]),
    selectedResultIndex: ref(0),
    svgContainer: ref({ querySelector: (selector) => (selector === 'svg' ? {} : null) }),
    errorLog: ref(null)
  };
  const ports = {};
  const owner = createFeatureVisibilityActions({
    state: withDrawings(portState),
    rulePreparation: rulePreparationFor(portState),
    getCommittedRequest: committedRequest(['CDS']),
    applyVisibilityPreviewChanges: (changes) => {
      projections += 1;
      return changes.reduce((changed, { featureId, mode }) => {
        if (mounted[featureId] === mode) return changed;
        mounted[featureId] = mode;
        return true;
      }, false);
    },
    ports,
    selectResult: () => true
  });
  // The root registers the port after both owners exist.
  ports.applyFeatureVisibilityToLabels = (options) => follows.push(options);
  const featureOf = (id) => resultCatalogFeatures(portState).biological
    .find((feature) => feature.biological_feature_id === id);
  setFeatureVisibilityOverride(overrides, featureOf('B'), 'off');

  assert.equal(await owner.projectFeatureVisibility({ rerender: true }), false);
  assert.deepEqual(follows, [], 'nothing to follow when the Result draws what Generate draws');

  // OV-35: a Result display (or a History step) that hides A hides its label.
  setFeatureVisibilityOverride(overrides, featureOf('A'), 'off');
  assert.equal(await owner.projectFeatureVisibility(), true);
  assert.deepEqual(mounted, { 'svg-A': 'off' });
  assert.deepEqual(follows, [{ reflow: false, rerender: false }], 'the label follows the feature');

  // B is shown, which the Result does not draw.
  setFeatureVisibilityOverride(overrides, featureOf('A'), 'default');
  setFeatureVisibilityOverride(overrides, featureOf('B'), 'on');
  await owner.projectFeatureVisibility();
  assert.deepEqual(follows.at(-1), { reflow: false, rerender: false }, 'a Result display does not rerender');
  await owner.projectFeatureVisibility({ rerender: true });
  assert.deepEqual(follows.at(-1), { reflow: false, rerender: true }, 'a History step rerenders to draw B');
  await owner.projectFeatureVisibility({ rerender: true, reflow: true });
  assert.deepEqual(follows.at(-1), { reflow: true, rerender: true }, 'a loaded table also places the labels');
  setFeatureVisibilityOverride(overrides, featureOf('B'), 'off');
  await owner.projectFeatureVisibility({ rerender: true, reflow: true });
  assert.deepEqual(follows.at(-1), { reflow: true, rerender: false });

  setFeatureVisibilityOverride(overrides, featureOf('A'), 'off');
  const projectionsBefore = projections;
  const followCount = follows.length;
  const [waiting, latest] = await Promise.all([
    owner.projectFeatureVisibility({ rerender: true }),
    owner.projectFeatureVisibility({ rerender: true })
  ]);
  assert.deepEqual([waiting, latest], [false, true], 'the later projection supersedes the waiting one');
  assert.equal(projections, projectionsBefore + 1);
  assert.equal(follows.length, followCount + 1);
}

console.log('feature visibility action tests passed');
