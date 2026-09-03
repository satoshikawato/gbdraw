import assert from 'node:assert/strict';
import test from 'node:test';

import { compileDirectEditorMutationPlan } from '../../gbdraw/web/js/app/candidate-render.js';
import { biologicalFeatureKey } from '../../gbdraw/web/js/services/feature-catalog.js';

const stableKey = biologicalFeatureKey('record-a', 'feature-a');

const admission = () => ({
  resultNames: ['diagram.svg'],
  renderedTargetsByOverrideKey: new Map([[
    stableKey,
    [{ resultIndex: 0, renderedId: 'f0001' }]
  ]]),
  resultIndexesByRenderedId: new Map([['f0001', new Set([0])]])
});

test('empty, default, stale, renderer, manual-rule, and legacy inputs compile to EMPTY', () => {
  const currentAdmission = admission();
  const plan = compileDirectEditorMutationPlan({
    catalogAdmission: currentAdmission,
    featureColorOverrides: {
      stale: '#112233',
      [stableKey]: ''
    },
    featureStrokeOverrides: {
      stale: { strokeColor: '#223344' },
      [stableKey]: {}
    },
    featureVisibilityOverrides: {
      stale: 'off',
      f0001: 'default'
    },
    labelTextFeatureOverrides: { stale: 'ignored' },
    labelVisibilityOverrides: { f0001: 'default' },
    manualSpecificRules: [{ cap: 'manual', color: '#abcdef' }],
    legendEntries: [{ caption: 'file-derived', color: '#abcdef' }],
    originalLegendOrder: ['CDS'],
    addedLegendCaptions: new Set(['file-derived']),
    paletteDefinitions: { default: { CDS: '#000000' } },
    selectedPalette: 'default',
    comparisonSettings: { showLegend: false },
    legacyFeatures: [{ svg_id: 'f0001' }],
    suppressPairwiseIdentityLegend: true
  });

  assert.equal(plan.kind, 'EMPTY');
  assert.equal(Object.isFrozen(plan), true);
  assert.equal(Object.isFrozen(plan.operationsByResult), true);
  assert.equal(Object.values(plan.operationsByResult[0]).flat().length, 0);
});

test('rule-derived fill residue is excluded without matching rules against Features', () => {
  const plan = compileDirectEditorMutationPlan({
    catalogAdmission: admission(),
    featureColorOverrides: {
      [stableKey]: { color: '#abcdef', caption: 'manual' }
    },
    manualSpecificRules: [{ feat: 'CDS', cap: 'manual', color: '#abcdef' }]
  });
  assert.equal(plan.kind, 'EMPTY');
});

test('each direct editor domain and their combination compile to MUTATING', () => {
  const cases = {
    fill: { featureColorOverrides: { [stableKey]: '#112233' } },
    stroke: {
      featureStrokeOverrides: {
        [stableKey]: { strokeColor: '#223344', strokeWidth: 2 }
      }
    },
    visibility: { featureVisibilityOverrides: { f0001: 'off' } },
    labelText: { labelTextFeatureOverrides: { f0001: 'renamed' } },
    labelVisibility: { labelVisibilityOverrides: { f0001: 'off' } },
    legendFill: {
      legendEntries: [{ caption: 'CDS', originalCaption: 'CDS', color: '#334455' }],
      originalLegendOrder: ['CDS'],
      legendColorOverrides: { CDS: '#334455' }
    },
    legendStroke: {
      legendEntries: [{ caption: 'CDS', originalCaption: 'CDS', color: '#aaaaaa' }],
      originalLegendOrder: ['CDS'],
      legendStrokeOverrides: { CDS: { strokeColor: '#445566', strokeWidth: 3 } }
    },
    legendRename: {
      legendEntries: [{
        caption: 'Genes', originalCaption: 'CDS', color: '#aaaaaa', xPos: 20, yPos: 30
      }],
      originalLegendOrder: ['CDS']
    },
    legendDelete: {
      legendEntries: [],
      deletedLegendEntries: [{ caption: 'CDS', originalCaption: 'CDS' }],
      originalLegendOrder: ['CDS']
    },
    legendAdd: {
      legendEntries: [{
        caption: 'New', originalCaption: 'New', color: '#556677', xPos: 40, yPos: 50
      }],
      originalLegendOrder: ['CDS']
    },
    callerTransform: { transformSvg() {} }
  };

  Object.entries(cases).forEach(([name, options]) => {
    const plan = compileDirectEditorMutationPlan({
      catalogAdmission: admission(),
      ...options
    });
    assert.equal(plan.kind, 'MUTATING', name);
  });

  const combined = compileDirectEditorMutationPlan({
    catalogAdmission: admission(),
    ...cases.fill,
    ...cases.stroke,
    ...cases.visibility,
    ...cases.labelText,
    ...cases.labelVisibility,
    ...cases.legendFill,
    ...cases.legendStroke,
    transformSvg() {}
  });
  assert.equal(combined.kind, 'MUTATING');
  assert.equal(Object.isFrozen(combined.operationsByResult[0]), true);
  Object.values(combined.operationsByResult[0]).forEach((entries) => {
    assert.equal(Object.isFrozen(entries), true);
  });
  assert.deepEqual(combined.operationsByResult[0].labelText, [{
    renderedId: 'f0001',
    value: 'renamed'
  }]);
  assert.deepEqual(combined.operationsByResult[0].labelVisibility, [{
    renderedId: 'f0001',
    mode: 'off'
  }]);
});

test('plan construction uses admission indexes and never enumerates catalog Features', () => {
  const currentAdmission = admission();
  currentAdmission.catalog = {
    items: new Proxy([], {
      get() {
        throw new Error('Feature catalog was enumerated');
      }
    })
  };
  const plan = compileDirectEditorMutationPlan({
    catalogAdmission: currentAdmission,
    featureColorOverrides: { [stableKey]: '#112233' },
    manualSpecificRules: [{ feat: 'CDS', cap: 'manual', color: '#abcdef' }]
  });
  assert.equal(plan.kind, 'MUTATING');
  assert.deepEqual(plan.operationsByResult[0].featureFills, [{
    renderedId: 'f0001',
    color: '#112233'
  }]);
});
