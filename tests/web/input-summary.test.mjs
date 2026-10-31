// The input summary Python sends beside each Result is admitted as sent, or the Result is refused.
import assert from 'node:assert/strict';
import test from 'node:test';

import { drawnInputSummaries, inputSummaryOf, validateInputSummaries } from '../../gbdraw/web/js/services/input-summary.js';
import { visibilityPaintState } from '../../gbdraw/web/js/app/result-paint-record.js';

const results = [{ name: 'a.svg' }, { name: 'b.svg' }];
const summary = (overrides = {}) => ({
  resultIndex: 1,
  resultName: 'b.svg',
  lines: ['Record EMPTY001.1 (empty.gb), record 1 of 1, 1,000 bp: drew 0 features.'],
  warnings: [{ code: 'no_drawable_features', message: 'Record EMPTY001.1 (empty.gb) has no feature to draw: it has no features.' }],
  ...overrides
});

test('a reply without input summaries has none, and a valid one is kept as sent', () => {
  assert.deepEqual(validateInputSummaries(undefined, results), []);
  const admitted = validateInputSummaries([summary()], results);
  assert.deepEqual(admitted, [summary()]);
  assert.ok(Object.isFrozen(admitted));
  assert.equal(inputSummaryOf(admitted, 1, 'b.svg'), admitted[0]);
  assert.equal(inputSummaryOf(admitted, 0, 'a.svg'), null);
  assert.equal(inputSummaryOf(admitted, 1, 'a.svg'), null, 'a Result is found by index and name');
});

test('a summary that does not name its Result, or has another shape, refuses the reply', () => {
  for (const bad of [
    null,
    {},
    [summary({ resultName: 'a.svg' })],
    [summary({ resultIndex: 2, resultName: undefined })],
    [summary({ resultIndex: -1 })],
    [summary({ lines: [''] })],
    [summary({ warnings: [{ code: 'No Features', message: 'x' }] })],
    [summary({ warnings: [{ code: 'no_drawable_features', message: 'x', extra: 1 }] })],
    [{ ...summary(), extra: true }]
  ]) {
    assert.throws(() => validateInputSummaries(bad, results), (error) => (
      error.code === 'RESULT_INVALID' && error.context.field === 'inputSummaries'
    ), JSON.stringify(bad));
  }
});

// A live hide changes what the Result draws without Python: the summary keeps
// the Feature visibility its render drew, and the page compares it with the
// drawing's (`visibilityPaintState`, the owner the paint record reads too).
test('each summary keeps the Feature visibility its render drew', () => {
  const drawing = { featureOverrides: {}, featureVisibilityManualRules: [] };
  const drawn = drawnInputSummaries(validateInputSummaries([summary()], results), visibilityPaintState(drawing));
  assert.ok(Object.isFrozen(drawn) && Object.isFrozen(drawn[0]));
  assert.deepEqual(drawn[0], { ...summary(), drawnVisibility: visibilityPaintState(drawing) });
  assert.equal(inputSummaryOf(drawn, 1, 'b.svg'), drawn[0]);
  const hidden = { ...drawing, featureOverrides: { k: { recordKey: 'r', biologicalFeatureId: 'f1', featureVisibility: 'off' } } };
  assert.notEqual(drawn[0].drawnVisibility, visibilityPaintState(hidden), 'a hide changes the drawn visibility');
});
