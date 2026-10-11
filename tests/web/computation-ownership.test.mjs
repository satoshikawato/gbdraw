// CW-03 (docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md#computation-ownership):
// without bulk label edits the label table projection reads no label or
// feature input. Per-feature edits are identity rows (design Q4), which the
// request carries without a table.
import assert from 'node:assert/strict';
import test from 'node:test';

const { buildBulkLabelProjection } = await import(
  '../../gbdraw/web/js/services/label-override-table.js'
);

const metrics = [];
globalThis.__GBDRAW_TEST_HOOKS__ = { onStructuralMetric: (metric) => metrics.push(metric) };
const buildCount = () => metrics.filter(({ name }) => name === 'labelOverrideTableBuildCount').length;

// Any read of the label inputs throws, so an eager build cannot pass.
const sentinel = () => new Proxy([], {
  get(target, key) {
    if (key === Symbol.iterator || key === 'length' || typeof key === 'string') {
      throw new Error('CW-03: label input was read');
    }
    return Reflect.get(target, key);
  }
});
const identity = (featureId) => JSON.stringify(['rec', featureId]);
const labelTargets = [{ identityKey: identity('f1'), sourceText: 'alpha' }];

test('the input sentinel is armed', () => {
  assert.throws(() => buildBulkLabelProjection({ alpha: 'ALPHA' }, {
    labelTargets: sentinel(), featureOverrides: {}
  }), /CW-03: label input was read/);
});

test('without bulk label edits the projection reads no input', () => {
  metrics.length = 0;
  const result = buildBulkLabelProjection({}, { labelTargets: sentinel(), featureOverrides: sentinel() });
  assert.deepEqual(result, { bulkLabelText: {}, rows: [] });
  assert.equal(buildCount(), 0);
});

test('a bulk edit builds once: known labels take its text, the rest a table rule', () => {
  metrics.length = 0;
  const result = buildBulkLabelProjection({ alpha: 'ALPHA', gamma: 'GAMMA', delta: '' }, {
    labelTargets,
    featureOverrides: {
      [identity('f2')]: { recordKey: 'rec', biologicalFeatureId: 'f2', labelSourceText: 'alpha' }
    }
  });
  assert.equal(buildCount(), 1);
  assert.equal(metrics[0].featureCount, labelTargets.length);
  assert.deepEqual(result.bulkLabelText, { [identity('f1')]: 'ALPHA', [identity('f2')]: 'ALPHA' });
  // A blank text hides labels; only the table rule says that.
  assert.deepEqual(result.rows, ['*\t*\tlabel\t^delta$\t', '*\t*\tlabel\t^gamma$\tGAMMA']);
});
