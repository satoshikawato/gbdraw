// CW-03 (docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md#computation-ownership):
// empty label and visibility overrides build no feature metadata or index.
import assert from 'node:assert/strict';
import test from 'node:test';

const { buildLabelOverrideRows, buildLabelOverrideTsv } = await import(
  '../../gbdraw/web/js/app/feature-editor/label-override-table.js'
);

const metrics = [];
globalThis.__GBDRAW_TEST_HOOKS__ = { onStructuralMetric: (metric) => metrics.push(metric) };
const buildCount = () => metrics.filter(({ name }) => name === 'labelOverrideTableBuildCount').length;

// Any read of the feature inputs throws, so an eager metadata build cannot pass.
const sentinel = () => new Proxy([], {
  get(target, key) {
    if (key === Symbol.iterator || key === 'length' || typeof key === 'string') {
      throw new Error('CW-03: feature metadata input was read');
    }
    return Reflect.get(target, key);
  }
});
const features = [
  { svg_id: 'f1', record_id: 'rec', feature_type: 'CDS', qualifiers: { gene: ['alpha'], locus_tag: ['L1'] } },
  { svg_id: 'f2', record_id: 'rec', feature_type: 'CDS', qualifiers: { gene: ['beta'], locus_tag: ['L2'] } }
];
const labels = [{ key: 'k1', featureId: 'f1', sourceText: 'alpha', text: 'alpha' }];

test('the metadata sentinel is armed', () => {
  assert.throws(() => buildLabelOverrideRows({ f1: 'renamed' }, {}, {
    extractedFeatures: sentinel(), editableLabels: []
  }), /CW-03: feature metadata input was read/);
});

test('empty overrides return an empty table without reading feature inputs', () => {
  metrics.length = 0;
  for (const visibilityOverrides of [undefined, {}]) {
    const result = buildLabelOverrideTsv({}, {}, {
      extractedFeatures: sentinel(), editableLabels: sentinel(), visibilityOverrides
    });
    assert.equal(result.tsv, '');
    assert.deepEqual(result.rows, []);
  }
  assert.equal(buildCount(), 0);
});

test('each non-empty override kind builds exactly once and emits its row', () => {
  const cases = [
    ['label text', { f1: 'renamed' }, {}, {}],
    ['visibility', {}, {}, { f2: 'off' }],
    ['bulk', {}, { alpha: 'ALPHA' }, {}]
  ];
  for (const [kind, featureOverrides, bulkOverrides, visibilityOverrides] of cases) {
    metrics.length = 0;
    const result = buildLabelOverrideTsv(featureOverrides, bulkOverrides, {
      extractedFeatures: features, editableLabels: labels, visibilityOverrides
    });
    // Positive control: the probe is live and records the real construction.
    assert.equal(buildCount(), 1, kind);
    assert.equal(metrics[0].featureCount, features.length, kind);
    assert.equal(result.rows.length > 0, true, kind);
    assert.match(result.tsv, /\t/, kind);
  }
});
