// The input summary Python sends beside each Result is admitted as sent, or the Result is refused.
import assert from 'node:assert/strict';
import test from 'node:test';

import { NOTICE_LINE_LIMIT, drawnInputSummaries, inputSummaryOf, validateInputSummaries } from '../../gbdraw/web/js/services/input-summary.js';
import { visibilityPaintState } from '../../gbdraw/web/js/app/result-paint-record.js';

const results = [{ name: 'a.svg' }, { name: 'b.svg' }];
const summary = (overrides = {}) => ({
  resultIndex: 1,
  resultName: 'b.svg',
  lines: ['Record EMPTY001.1 (empty.gb), record 1 of 1, 1,000 bp: drew 0 features.'],
  warnings: [{ code: 'no_drawable_features', message: 'Record EMPTY001.1 (empty.gb) has no feature to draw: it has no features.' }],
  records: [{ sourceLength: 1000, region: null }],
  linearWrap: null,
  ...overrides
});
// gbdraw/api/input_summary.py `LinearWrapPayload` of phage lambda at Auto (tests/test_web_wrap_diagnostics.py).
// Its records ride the summary (`records`), as they do for a Result that cut nothing.
const wrap = (overrides = {}) => ({
  rowLength: 10000, rowCount: 5, requested: 'auto', recordRows: [5], drawnAxisLength: 10000,
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

// Each drawn record: its source record's bases and the region of the source the request asked for (review 3 F2):
// a cropped record says which, so the Web compares it with the card. Its drawn bases are Python's to derive (D-71).
test('a summary keeps each record a cropped Result drew, as Python sent it', () => {
  const records = [{ sourceLength: 5000000, region: { start: 1, end: 200000 } }];
  assert.deepEqual(validateInputSummaries([summary({ records })], results)[0].records, records);
});

test('a summary keeps the rows its wrapped Result drew, as Python sent them', () => {
  for (const linearWrap of [wrap(), wrap({ requested: 20000, rowLength: 20000, rowCount: 3, recordRows: [3], drawnAxisLength: 20000 })]) {
    assert.deepEqual(validateInputSummaries([summary({ linearWrap })], results)[0].linearWrap, linearWrap);
  }
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
    [{ ...summary(), extra: true }],
    [(({ linearWrap, ...rest }) => rest)(summary())],
    [(({ records, ...rest }) => rest)(summary())],
    [{ ...summary(), recordLengths: [1000] }],
    [summary({ records: [] })],
    [summary({ records: [{ sourceLength: 0, region: null }] })],
    [summary({ records: [{ sourceLength: 1000.5, region: null }] })],
    [summary({ records: [{ region: null }] })],
    [summary({ records: [{ sourceLength: 1000 }] })],
    [summary({ records: [{ sourceLength: 1000, region: null, extra: 1 }] })],
    [summary({ records: [{ drawnLength: 1000, sourceLength: 1000, region: null }] })],
    [summary({ records: [{ sourceLength: 1000, region: { start: 1 } }] })],
    [summary({ records: [{ sourceLength: 1000, region: { start: 0, end: 1000 } }] })],
    [summary({ records: [{ sourceLength: 1000, region: [1, 1000] }] })],
    [summary({ linearWrap: {} })],
    [summary({ linearWrap: wrap({ requested: 'none' }) })],
    [summary({ linearWrap: wrap({ requested: '10kb' }) })],
    [summary({ linearWrap: wrap({ rowLength: 0 }) })],
    [summary({ linearWrap: wrap({ rowCount: 2.5 }) })],
    [summary({ linearWrap: wrap({ recordRows: [] }) })],
    [summary({ linearWrap: wrap({ recordRows: [5, 1] }) })],
    [summary({ linearWrap: wrap({ records: [] }) })],
    [summary({ linearWrap: wrap({ drawnAxisLength: undefined }) })],
    [summary({ linearWrap: { ...wrap(), extra: 1 } })]
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
  assert.deepEqual(drawn[0], {
    ...summary(), drawnVisibility: visibilityPaintState(drawing),
    warningGroups: [{ kind: 'left-out', shown: summary().warnings, rest: [] }]
  });
  assert.equal(inputSummaryOf(drawn, 1, 'b.svg'), drawn[0]);
  const hidden = { ...drawing, featureOverrides: { k: { recordKey: 'r', biologicalFeatureId: 'f1', featureVisibility: 'off' } } };
  assert.notEqual(drawn[0].drawnVisibility, visibilityPaintState(hidden), 'a hide changes the drawn visibility');
  const labelOnly = { ...drawing, featureOverrides: { k: { recordKey: 'r', biologicalFeatureId: 'f1', featureVisibility: null, labelVisibility: 'off' } } };
  assert.equal(visibilityPaintState(labelOnly), drawn[0].drawnVisibility, 'a label-only row sets no visibility');
});

// OV-2954 (Owner 2026-11-09 「Q2=C」): a comparison file whose record IDs repeat is drawn in full, so its warning
// alone does not say part of the input was not drawn (Python's INPUT_WARNINGS_OF_DRAWN_INPUT).
test('a warning of an input drawn in full does not say part of the input was not drawn', () => {
  const repeated = { code: 'comparison_record_id_repeated', message: 'Comparison sequence file a.fna: record ID X names more than one record.' };
  const [alone, both] = drawnInputSummaries(validateInputSummaries([
    summary({ resultIndex: 0, resultName: 'a.svg', warnings: [repeated] }),
    summary({ warnings: [repeated, ...summary().warnings] })
  ], results), '');
  assert.deepEqual(alone.warningGroups, [{ kind: 'drawn-in', shown: [repeated], rest: [] }]);
  // #1777 review 2, LOW 3: the "not drawn" heading sits over the warnings that leave input out only.
  assert.deepEqual(both.warningGroups, [
    { kind: 'left-out', shown: summary().warnings, rest: [] }, { kind: 'drawn-in', shown: [repeated], rest: [] }
  ]);
});

// OV-3495: a Linear LOSATP run of one record against 60 contigs warned of 59 empty record pairs, and the notice
// listed every one above the figure and pushed it out of view. Each group shows at most NOTICE_LINE_LIMIT lines:
// the first warnings and an "and {count} more" line that opens the rest.
test('a group of many warnings shows the first few and sets the rest apart', () => {
  const pairs = Array.from({ length: 59 }, (_, index) => ({
    code: 'comparison_table_empty',
    message: `Comparison of record #1 NZ_CP065372.1 and record #${index + 2} C${index}.1 has no rows, so no match is drawn.`
  }));
  const repeated = { code: 'comparison_record_id_repeated', message: 'Comparison sequence file a.fna: record ID X names more than one record.' };
  const [drawn] = drawnInputSummaries(validateInputSummaries([summary({ warnings: [...pairs, repeated] })], results), '');
  assert.ok(NOTICE_LINE_LIMIT >= 2 && NOTICE_LINE_LIMIT <= 5, 'a few lines, not the list');
  assert.deepEqual(drawn.warningGroups, [
    { kind: 'left-out', shown: pairs.slice(0, NOTICE_LINE_LIMIT - 1), rest: pairs.slice(NOTICE_LINE_LIMIT - 1) },
    { kind: 'drawn-in', shown: [repeated], rest: [] }
  ]);
  const fit = pairs.slice(0, NOTICE_LINE_LIMIT);
  const [fits] = drawnInputSummaries(validateInputSummaries([summary({ warnings: fit })], results), '');
  assert.deepEqual(fits.warningGroups, [{ kind: 'left-out', shown: fit, rest: [] }], 'a list that fits is shown whole');
});
