import assert from 'node:assert/strict';
import test from 'node:test';

import { compileDirectEditorMutationPlan } from '../../gbdraw/web/js/app/candidate-render.js';
import { biologicalFeatureKey } from '../../gbdraw/web/js/services/feature-catalog.js';

const stableKey = biologicalFeatureKey('record-a', 'feature-a');
// Per-feature visibility and label edits are identity rows (design Q4); each
// Result reaches them through the rendered targets of their identity.
const identityRow = (fields, featureId = 'feature-a') => ({
  [JSON.stringify(['record-a', featureId])]: {
    recordKey: 'record-a',
    biologicalFeatureId: featureId,
    featureVisibility: null,
    labelVisibility: null,
    labelText: null,
    labelSourceText: null,
    ...fields
  }
});

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
    featureOverrides: {
      ...identityRow({ featureVisibility: 'off', labelText: 'ignored' }, 'stale'),
      ...identityRow({ featureVisibility: 'exclude_matching', labelSourceText: 'source' })
    },
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
    visibility: { featureOverrides: identityRow({ featureVisibility: 'off' }) },
    labelText: { featureOverrides: identityRow({ labelText: 'renamed' }) },
    labelVisibility: { featureOverrides: identityRow({ labelVisibility: 'off' }) },
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
    featureOverrides: identityRow({ featureVisibility: 'off', labelText: 'renamed', labelVisibility: 'off' }),
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


test('decoration transforms mutate only their matched batch output; zero keeps EMPTY', () => {
  const catalogAdmission = { ...admission(), resultNames: ['a', 'b'] };
  const transform = () => {};
  const plan = compileDirectEditorMutationPlan({ catalogAdmission, resultTransforms: [null, transform] });
  assert.equal(plan.operationsByResult[0].callerTransforms.length, 0);
  assert.deepEqual(plan.operationsByResult[1].callerTransforms, [transform]);
  assert.equal(compileDirectEditorMutationPlan({ catalogAdmission, resultTransforms: [null, null] }).kind, 'EMPTY');
});

test('an edited Legend order compiles to one order operation; the default order compiles none (D-08)', async () => {
  const cds = { caption: 'CDS', originalCaption: 'CDS', color: '#111111', xPos: 22, yPos: 7 };
  const gc = { caption: 'GC percent', originalCaption: 'GC content', color: '#222222', xPos: 22, yPos: 31 };
  const manual = { caption: 'Manual', originalCaption: 'Manual', color: '#333333', xPos: 22, yPos: 55 };
  const compile = (legendEntries) => compileDirectEditorMutationPlan({
    catalogAdmission: admission(),
    legendEntries,
    originalLegendOrder: ['CDS', 'GC content']
  }).operationsByResult[0];

  const defaultOrder = compile([cds, gc, manual]);
  assert.equal(defaultOrder.legendOrder.length, 0);
  assert.deepEqual(defaultOrder.legendRenames.map(({ from, to, xPos, yPos }) => ({ from, to, xPos, yPos })), [
    { from: 'GC content', to: 'GC percent', xPos: 22, yPos: 31 }
  ]);

  const reordered = compile([gc, manual, cds]);
  assert.deepEqual(reordered.legendOrder.map(({ captions }) => [...captions]), [['GC percent', 'Manual', 'CDS']]);
  // The order owns the slots, so a renamed entry does not replay its old anchor.
  assert.deepEqual(reordered.legendRenames.map(({ xPos, yPos }) => [xPos, yPos]), [[null, null]]);

  const { installFakeSvgDom } = await import('./fake-svg-dom.mjs');
  const { applyEditorOperationsToMountedSvg } = await import('../../gbdraw/web/js/services/svg-result-ingestion.js');
  installFakeSvgDom();
  const entry = (caption, y) => `<g data-legend-key="${caption}"><path fill="#123456" transform="translate(0, ${y})"/><text transform="translate(22, ${y})"/></g>`;
  const svg = new DOMParser().parseFromString(
    `<svg viewBox="0 0 100 100"><g id="legend"><g id="feature_legend">${entry('CDS', 7)}${entry('GC content', 31)}${entry('Other', 55)}</g></g></svg>`
  ).documentElement;
  applyEditorOperationsToMountedSvg(svg, { ...reordered, legendAdds: [] });
  const placed = svg.querySelectorAll('g[data-legend-key]').map((group) => [
    group.getAttribute('data-legend-key'), group.querySelector('text').getAttribute('transform')
  ]);
  // Missing captions keep their place after the ordered ones.
  assert.deepEqual(placed, [
    ['GC percent', 'translate(22, 7)'],
    ['CDS', 'translate(22, 31)'],
    ['Other', 'translate(22, 55)']
  ]);
});

test('canvas padding reaches every candidate Result once and is idempotent (D-09)', async () => {
  const { installFakeSvgDom } = await import('./fake-svg-dom.mjs');
  const { applyCanvasPaddingToSvg } = await import('../../gbdraw/web/js/app/legend-layout/composition-actions.js');
  const { captureDecorationContinuity } = await import('../../gbdraw/web/js/app/legend-layout/decoration-continuity.js');
  installFakeSvgDom();
  const parse = () => new DOMParser().parseFromString('<svg viewBox="0 0 100 50" width="100px" height="50px"></svg>').documentElement;

  assert.equal(captureDecorationContinuity({ results: [], canvasPadding: { top: 0, right: 0, bottom: 0, left: 0 } }), null);
  const continuity = captureDecorationContinuity({ results: [], canvasPadding: { top: 0, right: 150, bottom: 0, left: 10 } });
  const transforms = continuity(null, { catalog: { items: [{}, {}] } });
  assert.equal(transforms.length, 2);
  const svgs = [parse(), parse()];
  transforms.forEach((transform, index) => transform(svgs[index]));
  svgs.forEach((svg) => {
    assert.equal(svg.getAttribute('viewBox'), '-10 0 260 50');
    assert.equal(svg.getAttribute('width'), '260px');
    assert.equal(svg.getAttribute('data-original-view-box'), '0 0 100 50');
  });
  assert.equal(applyCanvasPaddingToSvg(svgs[0], { right: 150, left: 10 }), false);
  assert.equal(svgs[0].getAttribute('viewBox'), '-10 0 260 50');
  assert.equal(applyCanvasPaddingToSvg(svgs[0], {}), true);
  assert.equal(svgs[0].getAttribute('viewBox'), '0 0 100 50');
  assert.equal(svgs[0].getAttribute('width'), '100px');
  // Zero padding leaves an unpadded Result byte-identical.
  const plain = parse();
  assert.equal(applyCanvasPaddingToSvg(plain, {}), false);
  assert.equal(plain.getAttribute('data-original-view-box'), null);
});

// Design Q4 6.2: an identity row reaches each Result through the rendered ID
// that Result draws it with, not through one rendered ID string.
test('an identity row projects onto every Result that draws its feature', () => {
  const plan = compileDirectEditorMutationPlan({
    catalogAdmission: {
      resultNames: ['one.svg', 'two.svg'],
      renderedTargetsByOverrideKey: new Map([[stableKey, [
        { resultIndex: 0, renderedId: 'f0001_record_1' },
        { resultIndex: 1, renderedId: 'f0001__instance_record_2_0123456789abcdef' }
      ]]]),
      resultIndexesByRenderedId: new Map()
    },
    featureOverrides: identityRow({ featureVisibility: 'off', labelText: 'renamed' })
  });
  assert.deepEqual(plan.operationsByResult.map((operations) => operations.featureVisibility), [
    [{ renderedId: 'f0001_record_1', mode: 'off' }],
    [{ renderedId: 'f0001__instance_record_2_0123456789abcdef', mode: 'off' }]
  ]);
  assert.deepEqual(plan.operationsByResult.map((operations) => operations.labelText.map(({ value }) => value)), [
    ['renamed'], ['renamed']
  ]);
});

// OV-45: a batch Result that draws none of a styled Legend row's features draws
// no row for it, so its Legend style may find the row absent; the Result that
// draws them still requires it.
test('a Legend style may miss its row only in a Result that draws none of its features', () => {
  const plan = compileDirectEditorMutationPlan({
    catalogAdmission: { ...admission(), resultNames: ['record-a.svg', 'record-b.svg'] },
    featureColorOverrides: { [stableKey]: { color: '#123456', caption: 'codon start two' } },
    manualSpecificRules: [{ feat: 'CDS', qual: 'hash', val: 'fef810304', color: '#123456', cap: 'codon start two' }],
    legendEntries: [{ caption: 'codon start two', originalCaption: 'codon start two', color: '#123456' }],
    originalLegendOrder: ['codon start two', 'CDS'],
    legendColorOverrides: { 'codon start two': '#123456' },
    legendStrokeOverrides: { 'codon start two': { strokeColor: '#445566', strokeWidth: 2 } }
  });
  assert.deepEqual(plan.operationsByResult.map(({ legendFills }) => legendFills), [
    [{ caption: 'codon start two', color: '#123456', allowMissing: false }],
    [{ caption: 'codon start two', color: '#123456', allowMissing: true }]
  ]);
  assert.deepEqual(plan.operationsByResult.map(({ legendStrokes }) => legendStrokes.map(
    ({ allowMissing, renderedIds }) => ({ allowMissing, renderedIds })
  )), [
    [{ allowMissing: false, renderedIds: ['f0001'] }],
    [{ allowMissing: true, renderedIds: [] }]
  ]);
});
