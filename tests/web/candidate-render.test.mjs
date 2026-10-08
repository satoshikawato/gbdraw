import assert from 'node:assert/strict';
import test from 'node:test';

import { compileDirectEditorMutationPlan } from '../../gbdraw/web/js/app/candidate-render.js';
import { legendRowFeatureIds } from '../../gbdraw/web/js/services/legend-svg.js';
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
  resultIndexesByRenderedId: new Map([['f0001', new Set([0])]]),
  // The rendered features of each Result with the renderer's fill.
  renderedFeaturesByResult: [new Map([['f0001', { fill_color: '#94a3b8' }]])]
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
  // The live Legend may draw the renamed row elsewhere (a longer caption wraps).
  const gc = { caption: 'GC percent', originalCaption: 'GC content', color: '#222222', xPos: 140, yPos: 7 };
  const manual = { caption: 'Manual', originalCaption: 'Manual', color: '#333333', xPos: 22, yPos: 55 };
  const compile = (legendEntries) => compileDirectEditorMutationPlan({
    catalogAdmission: admission(),
    legendEntries,
    originalLegendOrder: ['CDS', 'GC content']
  }).operationsByResult[0];

  const defaultOrder = compile([cds, gc, manual]);
  assert.equal(defaultOrder.legendOrder.length, 0);
  // A rename carries no position: the row keeps its place and the Legend
  // layout places it (OV-156).
  assert.deepEqual(defaultOrder.legendRenames, [{ from: 'GC content', to: 'GC percent', allowMissing: false }]);

  const reordered = compile([gc, manual, cds]);
  assert.deepEqual(reordered.legendOrder.map(({ captions }) => [...captions]), [['GC percent', 'Manual', 'CDS']]);
  assert.deepEqual(reordered.legendRenames, [{ from: 'GC content', to: 'GC percent', allowMissing: false }]);

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

  // OV-156: renaming a row that is not last, in the default order, keeps the
  // row in its place, even when the live Legend drew it elsewhere.
  const renamedInPlace = new DOMParser().parseFromString(
    `<svg viewBox="0 0 100 100"><g id="legend"><g id="feature_legend">${entry('CDS', 7)}${entry('GC content', 31)}${entry('Other', 55)}</g></g></svg>`
  ).documentElement;
  applyEditorOperationsToMountedSvg(renamedInPlace, { ...defaultOrder, legendAdds: [] });
  assert.deepEqual(renamedInPlace.querySelectorAll('g[data-legend-key]').map((group) => [
    group.getAttribute('data-legend-key'), group.querySelector('text').getAttribute('transform')
  ]), [
    ['CDS', 'translate(22, 7)'],
    ['GC percent', 'translate(22, 31)'],
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

// OV-46, OV-63: a Legend row with no known features is required by the plan in each
// Result. Admission lets Python's Legend row facts excuse the absence (suppressed by
// the draft, or drawn by another Result), so the plan carries no guess of its own.
test('a Legend row with no known features stays required in each Result of a batch', () => {
  const plan = compileDirectEditorMutationPlan({
    catalogAdmission: { ...admission(), resultNames: ['record-a.svg', 'record-b.svg'] },
    legendEntries: [{ caption: 'other proteins', originalCaption: 'other proteins', color: '#00aa00' }],
    originalLegendOrder: ['other proteins'],
    legendColorOverrides: { 'other proteins': '#00aa00' },
    legendStrokeOverrides: { 'other proteins': { strokeColor: '#445566', strokeWidth: 2 } }
  });
  assert.deepEqual(plan.operationsByResult.map(({ legendFills }) => legendFills), [
    [{ caption: 'other proteins', color: '#00aa00', allowMissing: false }],
    [{ caption: 'other proteins', color: '#00aa00', allowMissing: false }]
  ]);
  assert.deepEqual(plan.operationsByResult.map(({ legendStrokes }) => legendStrokes.map(
    ({ allowMissing, renderedIds }) => ({ allowMissing, renderedIds })
  )), [
    [{ allowMissing: false, renderedIds: [] }],
    [{ allowMissing: false, renderedIds: [] }]
  ]);
});

test('a Legend rename and style of a generated row are compiled the same way whatever the draft hides', () => {
  const compile = (extra = {}) => compileDirectEditorMutationPlan({
    catalogAdmission: admission(),
    legendEntries: [{ caption: 'Repeats', originalCaption: 'repeat_region', color: '#7b2cbf' }],
    originalLegendOrder: ['repeat_region'],
    legendColorOverrides: { Repeats: '#7b2cbf' },
    ...extra
  }).operationsByResult[0];
  const plain = compile();
  const hidden = compile({ featureOverrides: identityRow({ featureVisibility: 'off' }) });
  assert.deepEqual(hidden.legendFills, plain.legendFills);
  assert.deepEqual(hidden.legendRenames, plain.legendRenames);
  assert.deepEqual(plain.legendFills, [{ caption: 'repeat_region', color: '#7b2cbf', allowMissing: false }]);
  assert.equal(plain.legendRenames[0].allowMissing, false);
});

test('a Legend row in a single Result stays required in that Result', () => {
  const plan = compileDirectEditorMutationPlan({
    catalogAdmission: admission(),
    legendEntries: [{ caption: 'other proteins', originalCaption: 'other proteins', color: '#00aa00' }],
    originalLegendOrder: ['other proteins'],
    legendColorOverrides: { 'other proteins': '#00aa00' }
  });
  assert.deepEqual(plan.operationsByResult[0].legendFills, [
    { caption: 'other proteins', color: '#00aa00', allowMissing: false }
  ]);
});

// OV-81: Show Depth off sends Python no Depth source, so it cannot report the Depth
// row; the draft names it, and only that row may be absent (a Ghost row stays stale).
test('a Legend style on a Depth series the request left out may miss its row', () => {
  const compile = (unrequestedDepthCaptions) => compileDirectEditorMutationPlan({
    catalogAdmission: admission(),
    legendEntries: [
      { caption: 'depth', originalCaption: 'depth', color: '#7b2cbf' },
      { caption: 'Ghost', originalCaption: 'Ghost', color: '#123456' }
    ],
    originalLegendOrder: ['depth', 'Ghost'],
    legendColorOverrides: { depth: '#7b2cbf', Ghost: '#123456' },
    legendStrokeOverrides: { depth: { strokeColor: '#445566', strokeWidth: 2 } },
    unrequestedDepthCaptions
  }).operationsByResult[0];
  const requested = compile([]);
  assert.deepEqual(requested.legendFills.map(({ caption, allowMissing }) => [caption, allowMissing]), [
    ['depth', false], ['Ghost', false]
  ]);
  const hidden = compile(['depth']);
  assert.deepEqual(hidden.legendFills.map(({ caption, allowMissing }) => [caption, allowMissing]), [
    ['depth', true], ['Ghost', false]
  ]);
  assert.deepEqual(hidden.legendStrokes.map(({ caption, allowMissing }) => [caption, allowMissing]), [['depth', true]]);
});

// OV-88: one rule for the row Python is not sent: a Depth row renamed in the
// Legend is excused like an unrenamed one, its rename and the styles stored
// under the new name included.
test('a renamed Depth row the request left out may miss its row, rename and styles alike', () => {
  const compile = (unrequestedDepthCaptions) => compileDirectEditorMutationPlan({
    catalogAdmission: admission(),
    legendEntries: [
      { caption: 'Coverage', originalCaption: 'depth', color: '#7b2cbf' },
      { caption: 'Ghost name', originalCaption: 'Ghost', color: '#123456' }
    ],
    originalLegendOrder: ['depth', 'Ghost'],
    legendColorOverrides: { Coverage: '#7b2cbf', 'Ghost name': '#123456' },
    legendStrokeOverrides: { Coverage: { strokeColor: '#445566', strokeWidth: 2 } },
    unrequestedDepthCaptions
  }).operationsByResult[0];
  const pairs = (operations, key) => operations.map((operation) => [operation[key], operation.allowMissing]);
  const requested = compile([]);
  assert.deepEqual(pairs(requested.legendRenames, 'from'), [['depth', false], ['Ghost', false]]);
  assert.deepEqual(pairs(requested.legendFills, 'caption'), [['depth', false], ['Ghost', false]]);
  const hidden = compile(['depth']);
  assert.deepEqual(pairs(hidden.legendRenames, 'from'), [['depth', true], ['Ghost', false]]);
  assert.deepEqual(pairs(hidden.legendFills, 'caption'), [['depth', true], ['Ghost', false]]);
  assert.deepEqual(pairs(hidden.legendStrokes, 'caption'), [['depth', true]]);
});

// OV-123 (R3, PD-OI-066): one rule says which features a Legend row's stroke
// reaches. The live stroke reads the mounted Result and Generate the Result as
// it is drawn; both call it.
test('a Legend row stroke reaches the listed, named, and same-colored features without their own stroke', () => {
  const drawnFills = [['f1', '#54BCF8'], ['f1', '#54bcf8'], ['f2', '#d3d3d3'], ['f3', '#54bcf8'], ['f4', '#123456']];
  assert.deepEqual(legendRowFeatureIds({ color: '#54bcf8' }, { drawnFills }), ['f1', 'f3']);
  assert.deepEqual(legendRowFeatureIds({ color: '#54bcf8' }, { drawnFills, namedIds: ['f4'], ownStrokeIds: ['f3'] }), ['f1', 'f4']);
  assert.deepEqual(legendRowFeatureIds({ color: '#54bcf8', featureIds: ['f2', 'f9'] }, { drawnFills }), ['f2']);
  assert.deepEqual(legendRowFeatureIds({ color: 'none' }, { drawnFills }), []);
  assert.deepEqual(legendRowFeatureIds(undefined, { drawnFills, namedIds: ['f2'] }), ['f2']);
});

test('a stroke on a generated Legend row reaches the features each Result draws in its color', () => {
  const featureB = biologicalFeatureKey('record-a', 'feature-b');
  const plan = compileDirectEditorMutationPlan({
    catalogAdmission: {
      resultNames: ['record-a.svg', 'record-b.svg'],
      renderedTargetsByOverrideKey: new Map([
        [stableKey, [{ resultIndex: 0, renderedId: 'f0001' }]],
        [featureB, [{ resultIndex: 1, renderedId: 'f0002' }]]
      ]),
      resultIndexesByRenderedId: new Map([['f0001', new Set([0])], ['f0002', new Set([1])]]),
      renderedFeaturesByResult: [
        new Map([['f0001', { fill_color: '#54bcf8' }], ['f0003', { fill_color: '#54bcf8' }], ['f0004', { fill_color: '#d3d3d3' }]]),
        new Map([['f0002', { fill_color: '#d3d3d3' }], ['f0005', { fill_color: '#54bcf8' }]])
      ]
    },
    // f0001 keeps its own stroke; f0002 is drawn in the row's color by its fill edit.
    featureStrokeOverrides: { [stableKey]: { strokeColor: '#2a9d8f' } },
    featureColorOverrides: { [featureB]: { color: '#54bcf8' } },
    legendEntries: [{ caption: 'CDS', originalCaption: 'CDS', color: '#54bcf8' }],
    originalLegendOrder: ['CDS'],
    legendStrokeOverrides: { CDS: { strokeColor: '#e63946', strokeWidth: 3 } }
  });
  assert.deepEqual(plan.operationsByResult.map(({ legendStrokes }) => legendStrokes.map(
    ({ caption, allowMissing, renderedIds }) => ({ caption, allowMissing, renderedIds })
  )), [
    [{ caption: 'CDS', allowMissing: false, renderedIds: ['f0003'] }],
    [{ caption: 'CDS', allowMissing: false, renderedIds: ['f0002', 'f0005'] }]
  ]);
});

// OV-120: a renamed row that an earlier Generate did not draw (GC off, Show
// Depth off) waits in the drawing. Each later plan renames and styles the
// generated row where a Result draws it, and may miss it elsewhere; a row the
// drawing shows or deleted is not renamed twice. A shown row stays required.
test('a renamed row an earlier Generate hid is renamed and styled where a Result draws it (OV-120)', () => {
  const dormant = { caption: 'GC percent', originalCaption: 'GC content', color: '#7b2cbf' };
  const compile = (options) => compileDirectEditorMutationPlan({
    catalogAdmission: admission(),
    legendEntries: [{ caption: 'Ghost', originalCaption: 'Ghost', color: '#123456' }],
    originalLegendOrder: ['Ghost'],
    legendColorOverrides: { 'GC percent': '#7b2cbf', Ghost: '#123456' },
    legendStrokeOverrides: { 'GC percent': { strokeColor: '#445566', strokeWidth: 2 } },
    dormantLegendEntries: [dormant],
    ...options
  }).operationsByResult[0];
  const waiting = compile({});
  assert.deepEqual(waiting.legendRenames.map(({ from, to, allowMissing }) => [from, to, allowMissing]),
    [['GC content', 'GC percent', true]]);
  assert.deepEqual(waiting.legendFills.map(({ caption, allowMissing }) => [caption, allowMissing]), [
    ['GC content', true], ['Ghost', false]
  ]);
  assert.deepEqual(waiting.legendStrokes.map(({ caption, allowMissing }) => [caption, allowMissing]), [['GC content', true]]);
  // Without the waiting row, its style names a caption no Result draws.
  assert.deepEqual(compile({ dormantLegendEntries: [] }).legendRenames, []);
  // Shown again, or deleted, the row is not renamed from the waiting entry.
  const shown = compile({
    legendEntries: [{ caption: 'GC percent', originalCaption: 'GC content', color: '#7b2cbf' }],
    originalLegendOrder: ['GC content']
  });
  assert.deepEqual(shown.legendRenames.map(({ from, to, allowMissing }) => [from, to, allowMissing]),
    [['GC content', 'GC percent', false]]);
  const deleted = compile({
    deletedLegendEntries: [{ caption: 'GC content', originalCaption: 'GC content', color: '#a1a1a1' }],
    originalLegendOrder: ['Ghost', 'GC content']
  });
  assert.deepEqual(deleted.legendRenames, []);
});
