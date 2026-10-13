// The label fit report Python sends, the draft fields a row writes, and the
// dialog that shows the rows (LP-D43, LP-D45, LP-D48, OD-F14).
import assert from 'node:assert/strict';
import test from 'node:test';

globalThis.window = {
  Vue: {
    ref: (value) => ({ value }),
    reactive: (value) => value,
    computed: (getter) => ({ get value() { return getter(); } }),
    nextTick: async () => {}
  },
  DOMPurify: { sanitize: (value) => value }
};
globalThis.document = {};

const { state } = await import('../../gbdraw/web/js/state.js');
const { draftFieldsForConfigValues } = await import('../../gbdraw/web/js/services/session-request.js');
const { normalizeLabelFitFailure, normalizeLabelFitRow, normalizeResultLabelFits } = await import(
  '../../gbdraw/web/js/services/label-fit.js'
);
const { createLabelFitDialog } = await import('../../gbdraw/web/js/app/label-fit-dialog.js');

const Y_OFFSET = 'labels.unified_adjustment.outer_labels.y_radius_offset';
const row = (kind, status, extra = {}) => ({
  kind, title: kind.replace(/_/g, ' '), status, values: null, settings: [], canvasSize: null, reason: null, ...extra
});
const foundRow = (kind, value) => row(kind, 'found', {
  values: { [Y_OFFSET]: value }, settings: [{ name: 'Outer Y Offset', unit: '', value }],
  canvasSize: { width: 1708, height: 1130 }
});
const report = (rows, { drawn = false } = {}) => ({
  schema: 1, drawn, overlapPairs: drawn ? 0 : 3, crossingPairs: 5, definitionOverlaps: 0,
  records: [{ recordId: 'NC_1', labelKind: 'circular_outer', labelKindTitle: 'outer labels',
    labelCount: 80, fontSize: 18, overlapPairs: drawn ? 0 : 3, crossingPairs: 5, definitionOverlaps: 0,
    refused: false, counts: '3 label pairs overlap, 5 leader pairs cross' }],
  rows
});
const failure = (value, rowToken = null) => JSON.stringify({ report: value, rowToken });
const normalizeReport = (value) => normalizeLabelFitFailure(failure(value)).report;

test('the failure reader keeps Python\'s rows in order, the record counts, and the row token', () => {
  const rows = [foundRow('outward_vertical', 1.15), row('tighter_rows', 'not_found', { reason: 'No fit' }),
    row('combination', 'not_computed'),
    row('labels_out', 'found', { values: { 'labels.circular.scope': 'outer' },
      settings: [{ name: 'Label Mode', unit: '', value: 'Out' }], canvasSize: { width: 900, height: 900 } })];
  const { report: parsed, rowToken } = normalizeLabelFitFailure(failure(report(rows), 'circular:2:0'));
  assert.equal(rowToken, 'circular:2:0');
  assert.equal(parsed.records[0].counts, '3 label pairs overlap, 5 leader pairs cross');
  assert.deepEqual(parsed.rows[3].values, { 'labels.circular.scope': 'outer' });
  assert.equal(parsed.rows[3].settings[0].value, 'Out');
  assert.deepEqual(parsed.rows.map((item) => [item.kind, item.status]),
    [['outward_vertical', 'found'], ['tighter_rows', 'not_found'], ['combination', 'not_computed'], ['labels_out', 'found']]);
  assert.deepEqual(parsed.rows[0].values, { [Y_OFFSET]: 1.15 });
  assert.deepEqual(parsed.rows[0].canvasSize, { width: 1708, height: 1130 });
  assert.equal(parsed.rows[1].reason, 'No fit');
  assert.equal(normalizeLabelFitFailure(undefined), null);
  assert.ok(Object.isFrozen(parsed.rows[0].values));
});

test('a malformed report or row is RESULT_INVALID', () => {
  const invalid = (error) => error.code === 'RESULT_INVALID';
  assert.throws(() => normalizeLabelFitFailure('{'), invalid);
  assert.throws(() => normalizeLabelFitFailure(failure({ ...report([]), schema: 2 })), invalid);
  assert.throws(() => normalizeLabelFitFailure(failure(report([]), '')), invalid);
  const [record] = report([]).records;
  assert.throws(() => normalizeLabelFitFailure(failure({ ...report([]), records: [{ ...record, counts: 3 }] })), invalid);
  assert.throws(() => normalizeLabelFitRow(row('x', 'maybe')), invalid);
  // A found row has values and a canvas size.
  assert.throws(() => normalizeLabelFitRow(row('x', 'found')), invalid);
  assert.throws(() => normalizeLabelFitRow({ ...foundRow('x', 1), values: { [Y_OFFSET]: { wide: true } } }), invalid);
});

test('a Result\'s report names that Result by index and name', () => {
  const entry = { resultIndex: 0, resultName: 'fit.svg', report: report([], { drawn: true }), rowToken: 'circular:1:0' };
  const [parsed] = normalizeResultLabelFits([entry], [{ name: 'fit.svg' }]);
  assert.equal(parsed.rowToken, 'circular:1:0');
  assert.deepEqual(normalizeResultLabelFits(undefined, []), []);
  assert.throws(() => normalizeResultLabelFits([entry], [{ name: 'other.svg' }]), (error) => error.code === 'RESULT_INVALID');
});

test('a row\'s config paths become the draft fields the request mapping projects back', () => {
  const { adv } = state.drawings.circular;
  const draft = (values) => draftFieldsForConfigValues(values, 'circular', adv);
  assert.deepEqual(draft({ [Y_OFFSET]: 1.15 }), { form: {}, adv: { outer_label_y_offset: 1.15 } });
  assert.deepEqual(draft({ 'labels.font_size.short': 12, 'labels.font_size.long': 12 }),
    { form: {}, adv: { label_font_size: 12 } });
  assert.deepEqual(draft({ 'labels.spacing.circular': 0.5 }), { form: {}, adv: { circular_label_spacing: 0.5 } });
  // Word rows: every label outside, radial labels.
  assert.deepEqual(draft({ 'labels.circular.scope': 'outer' }), { form: { labels_mode: 'out' }, adv: {} });
  assert.deepEqual(draft({ 'labels.circular.placement': 'radial' }),
    { form: {}, adv: { circular_label_placement: 'radial' } });
  const unknown = (error) => error.code === 'INPUT_INVALID' && error.context.reason === 'UNKNOWN_CONFIG_PATH';
  assert.throws(() => draftFieldsForConfigValues({ 'labels.unknown': 1 }, 'circular', adv), unknown);
  assert.throws(() => draftFieldsForConfigValues({}, 'circular', adv), unknown);
});

const dialogHarness = ({ computeRow = async () => ({}) } = {}) => {
  const steps = [];
  const writes = [];
  const generates = [];
  const pending = { value: false };
  const closes = [];
  const dialog = createLabelFitDialog({
    reactive: (value) => value,
    computeRow,
    writeRowSettings: (values, mode) => writes.push({ values, mode }),
    withHistory: (label, handler, cancel) => async (choice) => {
      if (choice === 'cancel') return cancel();
      steps.push(label());
      pending.value = true;
      try { return handler(choice); } finally { pending.value = false; closes.splice(0).forEach((close) => close()); }
    },
    closeAfterChoice: (close) => closes.push(close),
    choicePending: () => pending.value,
    runAnalysis: async () => { generates.push(true); },
    errorSummary: () => 'The row could not be computed.'
  });
  return { dialog, steps, writes, generates };
};

test('Apply writes the chosen row as one History step, closes, and Generates; Cancel writes nothing', async () => {
  const rows = [row('outward_vertical', 'not_found', { reason: 'No fit' }), foundRow('smaller_font', 1.2)];
  const { dialog, steps, writes, generates } = dialogHarness();
  await dialog.openLabelFitDialog(normalizeReport(report(rows)), { mode: 'circular' });
  assert.equal(dialog.labelFitDialog.show, true);
  assert.equal(dialog.labelFitDialog.selectedKind, 'smaller_font');
  await dialog.resolveLabelFitChoice('apply');
  assert.deepEqual(steps, ['Fit labels: smaller font']);
  assert.deepEqual(writes, [{ values: { [Y_OFFSET]: 1.2 }, mode: 'circular' }]);
  assert.equal(generates.length, 1);
  assert.equal(dialog.labelFitDialog.show, false);

  await dialog.openLabelFitDialog(normalizeReport(report(rows)), { mode: 'circular' });
  await dialog.resolveLabelFitChoice('cancel');
  assert.equal(steps.length, 1);
  assert.equal(writes.length, 1);
  assert.equal(generates.length, 1);
  assert.equal(dialog.labelFitDialog.show, false);
});

test('rows Generate did not search are computed one after another when the dialog opens', async () => {
  const asked = [];
  const computeRow = async (token, kind) => {
    asked.push([token, kind]);
    if (kind === 'tighter_rows') throw new Error('worker');
    return { row: kind === 'outward_vertical' ? foundRow(kind, 1.3) : row(kind, 'not_found', { reason: 'No fit' }) };
  };
  const rows = [row('outward_vertical', 'not_computed'), row('tighter_rows', 'not_computed'),
    row('combination', 'not_computed')];
  const { dialog } = dialogHarness({ computeRow });
  await dialog.openLabelFitDialog(normalizeReport(report(rows, { drawn: true })),
    { rowToken: 'circular:1:0', mode: 'circular' });
  assert.deepEqual(asked, [['circular:1:0', 'outward_vertical'], ['circular:1:0', 'tighter_rows'],
    ['circular:1:0', 'combination']]);
  assert.deepEqual(dialog.labelFitDialog.rows.map((item) => item.status), ['found', 'failed', 'not_found']);
  assert.equal(dialog.labelFitDialog.rows[1].reason, 'The row could not be computed.');
  assert.equal(dialog.labelFitDialog.selectedKind, 'outward_vertical');
});

test('closing the dialog stops the rows still to compute', async () => {
  const asked = [];
  let release;
  const computeRow = (token, kind) => {
    asked.push(kind);
    return new Promise((resolve) => { release = () => resolve({ row: row(kind, 'not_found', { reason: 'No fit' }) }); });
  };
  const rows = [row('outward_vertical', 'not_computed'), row('combination', 'not_computed')];
  const { dialog } = dialogHarness({ computeRow });
  const computing = dialog.openLabelFitDialog(normalizeReport(report(rows, { drawn: true })),
    { rowToken: 'circular:1:0', mode: 'circular' });
  assert.equal(dialog.labelFitDialog.rows[0].status, 'computing');
  await dialog.resolveLabelFitChoice('cancel');
  release();
  await computing;
  assert.deepEqual(asked, ['outward_vertical']);
});
