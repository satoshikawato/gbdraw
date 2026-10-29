// The label fit report Python sends, the draft fields a row writes, and the
// dialog that shows the rows (LP-D43, LP-D45, LP-D48, OD-F14).
import assert from 'node:assert/strict';
import test from 'node:test';
import { modeDrawing } from './helpers/drawing-state.mjs';

import { t } from '../../gbdraw/web/js/utils/i18n.js';

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
const {
  labelCapacityNoticeText, labelCapacityText, normalizeLabelCapacityNotice, normalizeLabelFitFailure, normalizeLabelFitRow,
  normalizeResultLabelFits
} = await import('../../gbdraw/web/js/services/label-fit.js');
const { SOURCE_LOCALE, setUiLocale } = await import('../../gbdraw/web/js/utils/i18n.js');
const { createLabelFitDialog } = await import('../../gbdraw/web/js/app/label-fit-dialog.js');

const Y_OFFSET = 'labels.unified_adjustment.outer_labels.y_radius_offset';
const row = (kind, status, extra = {}) => ({
  kind, title: kind.replace(/_/g, ' '), status, values: null, trackOptions: {}, settings: [], canvasSize: null,
  reason: null, note: null, ...extra
});
const foundRow = (kind, value) => row(kind, 'found', {
  values: { [Y_OFFSET]: value }, settings: [{ name: 'Outer Y Offset', unit: '', value }],
  canvasSize: { width: 1708, height: 1130 }
});
const report = (rows, { drawn = false } = {}) => ({
  schema: 1, drawn, overlapPairs: drawn ? 0 : 3, crossingPairs: 5, definitionOverlaps: 0, regionOverlaps: 0,
  notice: drawn ? ['Leader lines cross (5 pairs).'] : [],
  records: [{ recordId: 'NC_1', recordIndex: 0, labelKind: 'circular_outer', labelKindTitle: 'outer labels',
    labelCount: 80, fontSize: 18, overlapPairs: drawn ? 0 : 3, crossingPairs: 5, definitionOverlaps: 0,
    regionOverlaps: 0, refused: false, capacityRatio: null, counts: '3 label pairs overlap, 5 leader pairs cross' }],
  rows,
  remedies: []
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
  const separate = normalizeLabelFitRow(row('separate_features', 'found', { values: { 'canvas.resolve_overlaps': true },
    settings: [{ name: 'Resolve Overlaps', unit: '', value: 'On' }], canvasSize: { width: 700, height: 700 } }));
  assert.deepEqual(separate.values, { 'canvas.resolve_overlaps': true });
  assert.equal(separate.settings[0].value, 'On');
  assert.deepEqual(parsed.rows[0].canvasSize, { width: 1708, height: 1130 });
  assert.equal(parsed.rows[1].reason, 'No fit');
  assert.equal(normalizeLabelFitFailure(undefined), null);
  assert.ok(Object.isFrozen(parsed.rows[0].values));
});

test('the notice bar shows Python\'s sentences and has no wording of its own', async () => {
  const { readFile } = await import('node:fs/promises');
  const parsed = normalizeReport(report([], { drawn: true }));
  assert.deepEqual(parsed.notice, ['Leader lines cross (5 pairs).']);
  assert.ok(Object.isFrozen(parsed.notice));
  const invalid = (error) => error.code === 'RESULT_INVALID';
  assert.throws(() => normalizeReport({ ...report([]), notice: [5] }), invalid);
  assert.throws(() => normalizeReport({ ...report([]), notice: undefined }), invalid);
  const html = await readFile(new URL('../../gbdraw/web/index.html', import.meta.url), 'utf8');
  const bar = html.slice(html.indexOf('data-testid="label-fit-notice"'), html.indexOf("{{ t('Choose a fix') }}</button>"));
  assert.match(bar, /labelFitNotice\.report\.notice/);
  assert.doesNotMatch(bar, /Leader lines cross|definition text|'pairs'/);
});

const noticeRow = () => row('separate_features', 'found_with_notice', { values: { 'canvas.resolve_overlaps': true },
  settings: [{ name: 'Resolve Overlaps', unit: '', value: 'On' }], canvasSize: { width: 1045, height: 873 },
  note: '14 leader pairs still cross' });

test('a row that leaves crossing leaders carries its settings and Python\'s note (LP-D57)', () => {
  const parsed = normalizeLabelFitRow(noticeRow());
  assert.equal(parsed.status, 'found_with_notice');
  assert.deepEqual(parsed.values, { 'canvas.resolve_overlaps': true });
  assert.deepEqual(parsed.canvasSize, { width: 1045, height: 873 });
  assert.equal(parsed.note, '14 leader pairs still cross');
  assert.equal(normalizeLabelFitRow(row('x', 'not_found')).note, null);
  const invalid = (error) => error.code === 'RESULT_INVALID';
  // Offered like a found row: it has values and a canvas size.
  assert.throws(() => normalizeLabelFitRow(row('x', 'found_with_notice')), invalid);
  assert.throws(() => normalizeLabelFitRow({ ...noticeRow(), note: 14 }), invalid);
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
  // A setting shows a number or the option word, never a raw switch.
  assert.throws(() => normalizeLabelFitRow({ ...foundRow('x', 1), settings: [{ name: 'Resolve Overlaps', unit: '', value: true }] }),
    invalid);
});

// LP-D68, LP-D82: a Linear figure drawn with region labels over other labels or
// lines warns with its count and offers no fix.
const regionReport = () => ({
  ...report([], { drawn: true }), crossingPairs: 0, regionOverlaps: 1,
  notice: ['Region labels overlap other labels or lines (1).'],
  records: [{ recordId: 'adjacent', recordIndex: null, labelKind: 'linear_region', labelKindTitle: 'region label', labelCount: 1,
    fontSize: 10, overlapPairs: 0, crossingPairs: 0, definitionOverlaps: 0, regionOverlaps: 1, refused: false,
    capacityRatio: null, counts: '1 region label overlap' }]
});

test('a region label report carries its count, and its notice offers no fix without rows (LP-D82)', async () => {
  const { readFile } = await import('node:fs/promises');
  const [entry] = normalizeResultLabelFits([{ resultIndex: 0, resultName: 'fit.svg', rowToken: null,
    report: regionReport() }], [{ name: 'fit.svg' }]);
  assert.equal(entry.report.regionOverlaps, 1);
  assert.equal(entry.report.records[0].regionOverlaps, 1);
  assert.deepEqual(entry.report.rows, []);
  assert.deepEqual(entry.report.notice, ['Region labels overlap other labels or lines (1).']);
  const invalid = (error) => error.code === 'RESULT_INVALID';
  assert.throws(() => normalizeReport({ ...regionReport(), regionOverlaps: -1 }), invalid);
  const html = await readFile(new URL('../../gbdraw/web/index.html', import.meta.url), 'utf8');
  assert.match(html, /<button v-if="labelFitNotice\.report\.rows\.length" [^>]*@click="openLabelFitNotice">\{\{ t\('Choose a fix'\) \}\}<\/button>/);
});

test('a Result\'s report names that Result by index and name', () => {
  const entry = { resultIndex: 0, resultName: 'fit.svg', report: report([], { drawn: true }), rowToken: 'circular:1:0' };
  const [parsed] = normalizeResultLabelFits([entry], [{ name: 'fit.svg' }]);
  assert.equal(parsed.rowToken, 'circular:1:0');
  assert.deepEqual(normalizeResultLabelFits(undefined, []), []);
  assert.throws(() => normalizeResultLabelFits([entry], [{ name: 'other.svg' }]), (error) => error.code === 'RESULT_INVALID');
});

test('a row\'s config paths become the draft fields the request mapping projects back', () => {
  const { adv } = modeDrawing(state, 'circular');
  const draft = (values) => draftFieldsForConfigValues(values, 'circular', adv);
  assert.deepEqual(draft({ [Y_OFFSET]: 1.15 }), { form: {}, adv: { outer_label_y_offset: 1.15 } });
  assert.deepEqual(draft({ 'labels.font_size.short': 12, 'labels.font_size.long': 12 }),
    { form: {}, adv: { label_font_size: 12 } });
  assert.deepEqual(draft({ 'labels.spacing.circular': 0.5 }), { form: {}, adv: { circular_label_spacing: 0.5 } });
  // Word rows: every label outside, radial labels.
  assert.deepEqual(draft({ 'labels.circular.scope': 'outer' }), { form: { labels_mode: 'out' }, adv: {} });
  assert.deepEqual(draft({ 'labels.circular.placement': 'radial' }),
    { form: {}, adv: { circular_label_placement: 'radial' } });
  // "Separate overlapping features" (LP-D52) switches the mode's Resolve Overlaps control on.
  assert.deepEqual(draft({ 'canvas.resolve_overlaps': true }), { form: {}, adv: { resolve_overlaps: true } });
  // "Smaller definition font" (D-12) sets Record Labels > Default font size.
  assert.deepEqual(draft({ 'objects.definition.circular.font_size': 16.5 }), { form: {}, adv: { def_font_size: 16.5 } });
  assert.deepEqual(draftFieldsForConfigValues({ 'canvas.resolve_overlaps': true }, 'linear', modeDrawing(state, 'linear').adv),
    { form: {}, adv: { resolve_overlaps: true } });
  // "Narrow the inside tracks" (PD-OI-097, N-a) sets the Circular Center Reserved Radius, a track option of the
  // request (by its request name, not a config path) that the request reads from the same draft field; no new
  // Session field.
  assert.deepEqual(draftFieldsForConfigValues({}, 'circular', adv, { centerReservedRadius: 231.99 }),
    { form: {}, adv: { center_reserved_radius: 231.99 } });
  const unknown = (error) => error.code === 'INPUT_INVALID' && error.context.reason === 'UNKNOWN_CONFIG_PATH';
  assert.throws(() => draftFieldsForConfigValues({}, 'circular', adv, { centerReservedRadius: '231.99' }), unknown);
  assert.throws(() => draftFieldsForConfigValues({}, 'linear', modeDrawing(state, 'linear').adv,
    { centerReservedRadius: 231.99 }), unknown);
  assert.throws(() => draftFieldsForConfigValues({}, 'circular', adv, { circularTrackAxisIndex: 1 }), unknown);
  // The old config-path spelling is no row value.
  assert.throws(() => draft({ 'tracks.center_reserved_radius': 231.99 }), unknown);
  assert.throws(() => draftFieldsForConfigValues({ 'labels.unknown': 1 }, 'circular', adv), unknown);
  assert.throws(() => draftFieldsForConfigValues({}, 'circular', adv), unknown);
});

test('a row may write exactly the config paths Python\'s rows send, in the mode that sends them', async () => {
  const { readFile } = await import('node:fs/promises');
  const { applyLabelFitValues } = await import('../../gbdraw/web/js/services/config.js');
  const fit = await readFile(new URL('../../gbdraw/labels/fit.py', import.meta.url), 'utf8');
  const constants = new Map([...fit.matchAll(/^([A-Z_]+_PATHS?) = (.+)$/gm)]
    .map(([, name, value]) => [name, [...value.matchAll(/"([^"]+)"/g)].map(([, path]) => path)]));
  const {
    SEPARATE_FEATURES_PATH: [separate], FONT_SIZE_PATHS: fontSizes, CENTER_RESERVED_RADIUS_PATH: [trackOption], ...single
  } = Object.fromEntries(constants);
  assert.equal(fontSizes.length, 2);
  const sample = { [separate]: true, 'labels.circular.scope': 'outer', 'labels.circular.placement': 'radial' };
  const unknown = (error) => error.code === 'INPUT_INVALID' && error.context.reason === 'UNKNOWN_CONFIG_PATH';
  // Python's row keeps the track option's own key; the Web payload sends it as a request track option instead.
  assert.throws(() => applyLabelFitValues(modeDrawing(state, 'circular'), { [trackOption]: 1.1 }, 'circular'), unknown);
  // Circular rows: every path Python names; Linear rows: separate overlapping features only.
  applyLabelFitValues(modeDrawing(state, 'circular'), Object.fromEntries(fontSizes.map((path) => [path, 12])), 'circular');
  for (const [path] of Object.values(single)) {
    applyLabelFitValues(modeDrawing(state, 'circular'), { [path]: sample[path] ?? 1.1 }, 'circular');
  }
  applyLabelFitValues(modeDrawing(state, 'linear'), { [separate]: true }, 'linear');
  // Label settings no row sends are not row writes.
  for (const [values, mode] of [
    [{ 'labels.unified_adjustment.inner_labels.x_radius_offset': 1.1 }, 'circular'],
    [{ 'labels.unified_adjustment.inner_labels.y_radius_offset': 1.1 }, 'circular'],
    [{ 'labels.spacing.linear': 0.5 }, 'linear'],
    [{ 'labels.linear.scope': 'out' }, 'linear']
  ]) {
    assert.throws(() => applyLabelFitValues(modeDrawing(state, mode), values, mode), unknown, JSON.stringify(values));
  }
});

const dialogHarness = ({ computeRow = async () => ({}), onWrite = () => {} } = {}) => {
  const steps = [];
  const writes = [];
  const generates = [];
  const shown = [];
  const pending = { value: false };
  const closes = [];
  const dialog = createLabelFitDialog({
    showRecordText: (recordId) => { shown.push(recordId); return true; },
    reactive: (value) => value,
    computeRow,
    writeRowSettings: (settings, mode) => { writes.push({ ...settings, mode }); onWrite(); },
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
  return { dialog, steps, writes, generates, shown, pending };
};

const definitionRecord = (recordId, recordIndex) => ({
  recordId, recordIndex, labelKind: 'circular_definition', labelKindTitle: 'definition lines', labelCount: 1, fontSize: 18,
  overlapPairs: 0, crossingPairs: 0, definitionOverlaps: 0, regionOverlaps: 0, refused: false, capacityRatio: null,
  counts: '1 line overlaps a track'
});

test('Edit Species / Strain closes the dialog, applies nothing, and shows that record\'s field by its index (E-b)', async () => {
  const { dialog, steps, writes, generates, shown } = dialogHarness();
  const drawn = report([foundRow('smaller_definition_font', 16.5)], { drawn: true });
  // Two records of a Multi-Record Canvas may share a record ID: each has its own button, by its drawn index.
  drawn.records = [...drawn.records, definitionRecord('NC_2', 1), definitionRecord('NC_2', 2)];
  await dialog.openLabelFitDialog(normalizeReport(drawn), { mode: 'circular' });
  // One button per record whose definition is over a track (not the outer-label record).
  assert.deepEqual(dialog.labelFitDialog.definitionRecords, [{ recordIndex: 1, recordId: 'NC_2' }, { recordIndex: 2, recordId: 'NC_2' }]);
  assert.equal(await dialog.editLabelFitRecordText(2), true);
  assert.deepEqual(shown, [2]);
  assert.equal(dialog.labelFitDialog.show, false);
  assert.deepEqual([steps, writes, generates], [[], [], []]);
  // A report without a definition record offers no such button.
  await dialog.openLabelFitDialog(normalizeReport(report([foundRow('smaller_font', 1.2)])), { mode: 'circular' });
  assert.deepEqual(dialog.labelFitDialog.definitionRecords, []);
  // Nor does a Generate that drew nothing: no drawn definition to edit (the preview may show an older figure).
  const failed = report([foundRow('smaller_font', 1.2)]);
  failed.records = [...failed.records, definitionRecord('NC_2', 1)];
  await dialog.openLabelFitDialog(normalizeReport(failed), { mode: 'circular' });
  assert.deepEqual(dialog.labelFitDialog.definitionRecords, []);
});

test('a report record\'s index is a drawn record index or null, never missing', () => {
  for (const recordIndex of [undefined, -1, 1.5, '0']) {
    const value = report([], { drawn: true });
    value.records = [{ ...value.records[0], recordIndex }];
    assert.throws(() => normalizeReport(value), (error) => error.code === 'RESULT_INVALID', String(recordIndex));
  }
  const linear = report([], { drawn: true });
  linear.records = [{ ...linear.records[0], recordIndex: null }];
  assert.equal(normalizeReport(linear).records[0].recordIndex, null);
});

test('over a drawn figure the dialog\'s Cancel reads Keep as drawn (K-a), and Edit Species / Strain is a button per record', async () => {
  const { readFile } = await import('node:fs/promises');
  const html = await readFile(new URL('../../gbdraw/web/index.html', import.meta.url), 'utf8');
  const start = html.indexOf('<choice-dialog v-if="labelFitDialog.show"');
  const markup = html.slice(start, html.indexOf('</choice-dialog>', start));
  assert.match(markup, /data-label-fit-keep @click="resolveLabelFitChoice\('cancel'\)"[^>]*>\{\{ labelFitDialog\.capacity \? t\('Close'\) : labelFitDialog\.drawn \? t\('Keep as drawn'\) : t\('Cancel'\) \}\}<\/button>/);
  assert.match(markup, /v-for="record in labelFitDialog\.definitionRecords" :key="record\.recordIndex"[^>]*:data-label-fit-edit-record="record\.recordIndex" @click="editLabelFitRecordText\(record\.recordIndex\)"/);
  const setup = await readFile(new URL('../../gbdraw/web/js/app/app-setup.js', import.meta.url), 'utf8');
  assert.match(setup, /showRecordText: \(recordIndex\) => previewItemPopup\.showRecordTextSettings\(recordIndex\)/);
});

test('while Apply writes, the dialog is busy and Escape (its Cancel) does nothing (PD-OI-088)', async () => {
  const { readFile } = await import('node:fs/promises');
  const escapes = [];
  let harness;
  harness = dialogHarness({ onWrite: () => {
    escapes.push(harness.dialog.resolveLabelFitChoice('cancel'));
    escapes.push(harness.dialog.labelFitDialog.show);
  } });
  const { dialog, steps, writes, generates } = harness;
  await dialog.openLabelFitDialog(normalizeReport(report([foundRow('smaller_font', 1.2)])), { mode: 'circular' });
  await dialog.resolveLabelFitChoice('apply');
  assert.equal(escapes.length, 2);
  assert.equal(await escapes[0], undefined);
  assert.equal(escapes[1], true, 'the Escape during Apply left the dialog open');
  // One step, the row's settings, one Generate: the Escape during Apply wrote and closed nothing.
  assert.deepEqual(steps, ['Fit labels: smaller font']);
  assert.equal(writes.length, 1);
  assert.equal(generates.length, 1);
  // The dialog is the shared choice dialog bound to the choice in progress: aria-busy,
  // its controls disabled, and Escape emits the Cancel that waits.
  const html = await readFile(new URL('../../gbdraw/web/index.html', import.meta.url), 'utf8');
  const open = html.match(/<choice-dialog v-if="labelFitDialog\.show"[^>]*>/s)?.[0] || '';
  assert.match(open, /:busy="dialogChoicePending \|\| labelFitDialog\.finishing"/);
  assert.match(open, /@cancel="resolveLabelFitChoice\('cancel'\)"/);
  const template = html.slice(html.indexOf('id="choice-dialog-template"'), html.indexOf('</script>', html.indexOf('id="choice-dialog-template"')));
  assert.match(template, /@keydown\.esc\.prevent\.stop="\$emit\('cancel'\)"/);
  assert.match(template, /:aria-busy="busy \? 'true' : null"/);
  assert.match(template, /<fieldset :disabled="busy"/);
  const setup = await readFile(new URL('../../gbdraw/web/js/app/app-setup.js', import.meta.url), 'utf8');
  assert.match(setup, /const dialogChoicePending = dialogChoice\.pending;/);
  assert.match(setup, /choicePending: \(\) => dialogChoice\.pending\.value,/);
});

test('only a Generate that drew nothing or the notice bar opens the dialog; a rerender failure does not', async () => {
  const { readFile, readdir } = await import('node:fs/promises');
  const root = new URL('../../gbdraw/web/js/', import.meta.url);
  const files = (await readdir(root, { recursive: true })).filter((name) => name.endsWith('.js'));
  const calls = [];
  for (const name of files) {
    const text = await readFile(new URL(name, root), 'utf8');
    for (const match of text.matchAll(/openLabelFitDialog\(/g)) calls.push(`${name}:${text.slice(0, match.index).split('\n').length}`);
  }
  assert.equal(calls.length, 2, calls.join(', '));
  const setup = await readFile(new URL('app/app-setup.js', root), 'utf8');
  // Generate's own outcome (afterGenerate), and the notice bar's "Choose a fix".
  assert.match(setup, /if \(result\?\.status === 'error' && result\.labelFit\) \{\n\s*void openLabelFitDialog\(result\.labelFit\.report/);
  assert.match(setup, /const openLabelFitNotice = \(\) => \{[^}]*return openLabelFitDialog\(entry\.report/s);
  // A label reflow (the rerender of an edit) that fails reports a live-edit failure, not a label fit.
  const run = await readFile(new URL('app/run-analysis.js', root), 'utf8');
  const reflowStart = run.indexOf('const runLabelReflowCandidate = async');
  const reflow = run.slice(reflowStart, run.indexOf('\n  };\n', reflowStart));
  // Both failures of a rerender report it: the engine's error and a live
  // setting value Generate would refuse (its projection).
  const reported = /\{\n\s*(?:logPostGbdrawTimings\(timingEntries\);\n\s*)?const error = formatError\([^;]*\);\n\s*labelReflowLastError\.value = liveEditFailure\(error\);\n\s*return \{ status: 'error', error \};/;
  assert.match(reflow, new RegExp(`if \\(execution\\.status === 'engine-error'\\) ${reported.source}`));
  assert.match(reflow, new RegExp(`projectLiveSettingWrites\\([^;]*;\\n\\s*\\} catch \\(e\\) ${reported.source}`));
  assert.doesNotMatch(reflow, /labelFit:/);
});

test('Apply writes the chosen row as one History step, closes, and Generates; Cancel writes nothing', async () => {
  const rows = [row('outward_vertical', 'not_found', { reason: 'No fit' }), foundRow('smaller_font', 1.2)];
  const { dialog, steps, writes, generates } = dialogHarness();
  await dialog.openLabelFitDialog(normalizeReport(report(rows)), { mode: 'circular' });
  assert.equal(dialog.labelFitDialog.show, true);
  assert.equal(dialog.labelFitDialog.selectedKind, 'smaller_font');
  await dialog.resolveLabelFitChoice('apply');
  assert.deepEqual(steps, ['Fit labels: smaller font']);
  assert.deepEqual(writes, [{ values: { [Y_OFFSET]: 1.2 }, trackOptions: {}, mode: 'circular' }]);
  assert.equal(generates.length, 1);
  assert.equal(dialog.labelFitDialog.show, false);

  await dialog.openLabelFitDialog(normalizeReport(report(rows)), { mode: 'circular' });
  await dialog.resolveLabelFitChoice('cancel');
  assert.equal(steps.length, 1);
  assert.equal(writes.length, 1);
  assert.equal(generates.length, 1);
  assert.equal(dialog.labelFitDialog.show, false);
});

test('a row that leaves crossing leaders is offered after the rows that fit, and Apply writes it', async () => {
  const { dialog, steps, writes, generates } = dialogHarness();
  await dialog.openLabelFitDialog(normalizeReport(report([row('radial', 'not_found', { reason: 'No fit' }), noticeRow()])),
    { mode: 'circular' });
  const rows = dialog.labelFitDialog.rows;
  assert.deepEqual(rows.map((item) => [item.kind, item.offered]), [['radial', false], ['separate_features', true]]);
  assert.equal(rows[1].note, '14 leader pairs still cross');
  // The only offered row is preselected.
  assert.equal(dialog.labelFitDialog.selectedKind, 'separate_features');
  await dialog.resolveLabelFitChoice('apply');
  assert.deepEqual(steps, ['Fit labels: separate features']);
  assert.deepEqual(writes, [{ values: { 'canvas.resolve_overlaps': true }, trackOptions: {}, mode: 'circular' }]);
  assert.equal(generates.length, 1);

  // A row that fits comes first in Python's order and is preselected.
  await dialog.openLabelFitDialog(normalizeReport(report([foundRow('smaller_font', 1.2), noticeRow()])),
    { mode: 'circular' });
  assert.equal(dialog.labelFitDialog.selectedKind, 'smaller_font');
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

// CAPACITY-STOP (PD-OI-097): labels that need much more room than the map has
// were not placed, so there are no rows. Python names the ways to make room; the
// Web words each with its own controls, and the dialog offers only Close.
const CAPACITY_REMEDIES = ['filter', 'label_text', 'font', 'linear'];
const capacityReport = (remedies = CAPACITY_REMEDIES, capacityRatio = 2.5) => {
  const value = report([]);
  return { ...value, overlapPairs: 0, crossingPairs: 0, remedies,
    records: [{ ...value.records[0], recordId: 'AP027078.1', labelCount: 562, fontSize: 8, overlapPairs: 0,
      crossingPairs: 0, capacityRatio, counts: 'need about 2.5 times the room around the map' }] };
};

test('a capacity stop carries each record\'s ratio and Python\'s remedies; an unknown remedy is RESULT_INVALID', () => {
  const parsed = normalizeReport(capacityReport());
  assert.deepEqual(parsed.records[0].capacity,
    { labelKind: 'circular_outer', labelCount: 562, fontSize: 8, capacityRatio: 2.5 });
  assert.ok(Object.isFrozen(parsed.records[0].capacity));
  assert.equal('capacityRatio' in parsed.records[0], false);
  assert.equal(parsed.records[0].counts, 'need about 2.5 times the room around the map');
  assert.deepEqual(parsed.remedies, CAPACITY_REMEDIES);
  assert.ok(Object.isFrozen(parsed.remedies));
  assert.deepEqual(parsed.rows, []);
  assert.equal(normalizeReport(report([])).records[0].capacity, null);
  assert.deepEqual(normalizeReport(report([])).remedies, []);
  const invalid = (error) => error.code === 'RESULT_INVALID';
  for (const remedies of [['shrink'], ['Filter'], 'filter', [3], undefined, null]) {
    assert.throws(() => normalizeReport({ ...capacityReport(), remedies }), invalid, JSON.stringify(remedies));
  }
  for (const ratio of [undefined, '2.5', -1, 0, true]) {
    const value = capacityReport();
    value.records = [{ ...value.records[0], capacityRatio: ratio }];
    assert.throws(() => normalizeReport(value), invalid, String(ratio));
  }
  // Only the kinds whose room Python reads carry a ratio (the Web words them, OV-2179).
  assert.equal(normalizeReport({ ...capacityReport(), records: [{ ...capacityReport().records[0], labelKind: 'circular_inner' }] })
    .records[0].capacity.labelKind, 'circular_inner');
  for (const labelKind of ['circular_embedded', 'linear', 'outer']) {
    const value = capacityReport();
    value.records = [{ ...value.records[0], labelKind }];
    assert.throws(() => normalizeReport(value), invalid, labelKind);
  }
});

// OV-2179: the room labels need is worded by the Web from Python's values, in the UI language; the record line of
// the dialog and the overlay's notice are whole sentences per side and count.
test('the room labels need reads in the UI language, from Python\'s values', async () => {
  const outer = { labelKind: 'circular_outer', labelCount: 562, fontSize: 8, capacityRatio: 2.5 };
  // `{record}` is the dialog's slot for the record ID in bold.
  assert.equal(labelCapacityText(outer), '{record}: 562 outer labels at 8 pt need about 2.5 times the room around the map');
  assert.equal(labelCapacityText({ ...outer, labelKind: 'circular_inner', labelCount: 1, fontSize: 7.5 }),
    '{record}: 1 inner label at 7.5 pt needs about 2.5 times the room around the map');
  const notice = normalizeLabelCapacityNotice(JSON.stringify({ recordId: 'AP027078.1', ...outer, fontSize: 6,
    capacityRatio: 1.8 }));
  assert.ok(Object.isFrozen(notice));
  assert.equal(labelCapacityNoticeText(/** @type {NonNullable<typeof notice>} */ (notice)), 'AP027078.1: 562 outer labels'
    + ' at 6 pt need about 1.8 times the room around the map, so they may not fit. Placing them can take minutes.');
  const { default: ja } = await import('../../gbdraw/web/js/utils/i18n-catalogs/ja.js');
  // This file's stand-in document has no element to take the page language.
  const page = globalThis.document;
  globalThis.document = undefined;
  try {
    await setUiLocale('ja', async () => ({ default: ja }));
    for (const text of [labelCapacityText(outer), labelCapacityNoticeText(/** @type {NonNullable<typeof notice>} */ (notice))]) {
      assert.match(text, /562/);
      assert.match(text, /^(\{record\}|AP027078\.1)：/, 'Japanese puts a full-width colon after the record');
      assert.doesNotMatch(text, /room around the map|outer labels|may not fit/, text);
    }
  } finally {
    await setUiLocale(SOURCE_LOCALE);
    globalThis.document = page;
  }
  // Text of another shape is no notice: Python's sentence, a missing value, another label kind.
  for (const value of ['AP027078.1: 562 outer labels', JSON.stringify({ ...outer, recordId: 7 }), JSON.stringify(outer),
    JSON.stringify({ recordId: 'r', ...outer, labelKind: 'circular_embedded' }),
    JSON.stringify({ recordId: 'r', ...outer, capacityRatio: 0 }), undefined, 3]) {
    assert.equal(normalizeLabelCapacityNotice(value), null, String(value));
  }
});

test('the Web words exactly the remedies Python can name, in Python\'s order', async () => {
  const { readFile } = await import('node:fs/promises');
  const { LABEL_CAPACITY_REMEDY_WORDS } = await import('../../gbdraw/web/js/services/label-fit.js');
  const fit = await readFile(new URL('../../gbdraw/labels/fit.py', import.meta.url), 'utf8');
  const literal = fit.match(/^LabelCapacityRemedy: TypeAlias = Literal\[([^\]]+)\]$/m)?.[1] || '';
  const ids = [...literal.matchAll(/"([a-z_]+)"/g)].map(([, id]) => id);
  assert.deepEqual(ids, CAPACITY_REMEDIES);
  assert.deepEqual(Object.keys(LABEL_CAPACITY_REMEDY_WORDS), ids);
  assert.ok(Object.isFrozen(LABEL_CAPACITY_REMEDY_WORDS));
});

test('the Web words exactly the label kinds whose room Python reads, in Python\'s order', async () => {
  const { readFile } = await import('node:fs/promises');
  const { LABEL_CAPACITY_KINDS } = await import('../../gbdraw/web/js/services/label-fit.js');
  const fit = await readFile(new URL('../../gbdraw/labels/fit.py', import.meta.url), 'utf8');
  const literal = fit.match(/^HorizontalLabelKind: TypeAlias = Literal\[([^\]]+)\]$/m)?.[1] || '';
  const kinds = [...literal.matchAll(/"([a-z_]+)"/g)].map(([, kind]) => kind);
  assert.deepEqual(kinds, ['circular_outer', 'circular_inner']);
  assert.deepEqual(LABEL_CAPACITY_KINDS, kinds);
  assert.ok(Object.isFrozen(LABEL_CAPACITY_KINDS));
});

test('the dialog over a capacity stop lists the remedies in Python\'s order, no rows, and only Close', async () => {
  const { readFile } = await import('node:fs/promises');
  const { dialog, steps, writes, generates } = dialogHarness();
  await dialog.openLabelFitDialog(normalizeReport(capacityReport(['filter', 'label_text', 'linear'])), { mode: 'circular' });
  const state = dialog.labelFitDialog;
  assert.equal(state.show, true);
  assert.equal(state.capacity, true);
  assert.deepEqual(state.remedies, [
    'Label fewer features: Blacklist or Whitelist under LABEL FILTERING.',
    'Shorten the label text: a priority rule, for example gene before product.',
    'Draw a Linear map.'
  ]);
  assert.deepEqual(state.rows, []);
  assert.deepEqual(state.records.map((record) => [record.recordId, record.capacity]),
    [['AP027078.1', { labelKind: 'circular_outer', labelCount: 562, fontSize: 8, capacityRatio: 2.5 }]]);
  assert.equal(state.selectedKind, '');
  await dialog.resolveLabelFitChoice('cancel');
  assert.equal(state.show, false);
  assert.deepEqual([steps, writes, generates], [[], [], []]);
  // A report with rows has no remedies.
  await dialog.openLabelFitDialog(normalizeReport(report([foundRow('smaller_font', 1.2)])), { mode: 'circular' });
  assert.equal(state.capacity, false);
  assert.deepEqual(state.remedies, []);
  assert.deepEqual(state.records.map((record) => record.capacity), [null]);

  const html = await readFile(new URL('../../gbdraw/web/index.html', import.meta.url), 'utf8');
  const start = html.indexOf('<choice-dialog v-if="labelFitDialog.show"');
  const markup = html.slice(start, html.indexOf('</choice-dialog>', start));
  assert.match(markup, /\{\{ labelFitDialog\.capacity\s+\? t\('Nothing was drawn\. These labels need much more room than the map has, so placing them was not tried\.'\)/);
  assert.match(markup, /<div v-if="labelFitDialog\.remedies\.length"[^>]*data-label-fit-remedies>\s*<p[^>]*>\{\{ t\('To make room:'\) \}\}<\/p>\s*<ul[^>]*>\s*<li v-for="remedy in labelFitDialog\.remedies" :key="remedy">\{\{ tm\(remedy\) \}\}<\/li>/);
  // No rows section and no Apply over a capacity stop; the one button reads Close.
  assert.match(markup, /<fieldset v-if="labelFitDialog\.rows\.length"/);
  assert.match(markup, /<p v-else-if="!labelFitDialog\.capacity"[^>]*>\{\{ t\('No setting in this dialog makes these labels fit\.'\) \}\}<\/p>/);
  assert.match(markup, /<button v-if="!labelFitDialog\.capacity" [^>]*@click="resolveLabelFitChoice\('apply'\)"/);
  assert.match(markup, />\{\{ labelFitDialog\.capacity \? t\('Close'\) : labelFitDialog\.drawn \? t\('Keep as drawn'\) : t\('Cancel'\) \}\}<\/button>/);
  // A record not placed for lack of room reads in the Web's words (OV-2179); other records keep their counts in brackets.
  assert.match(markup, /<i18n-text v-if="record\.capacity" :text="labelCapacityText\(record\.capacity\)"><template #record><span[^>]*>\{\{ record\.recordId \}\}<\/span><\/template><\/i18n-text>\s*<template v-else><span[^>]*>\{\{ record\.recordId \}\}<\/span>: \{\{ record\.labelCount \}\} \{\{ record\.labelKindTitle \}\} at \{\{ record\.fontSize \}\} pt \(\{\{ record\.counts \}\}\)<\/template>/);
});

test('when every row does not fit, the dialog also lists Python\'s ways to make room (NO-ROW-WAYS, OV-2180)', async () => {
  const notFound = (kind) => row(kind, 'not_found', { reason: 'Does not fit down to 6 pt' });
  const ways = ['filter', 'label_text', 'linear'];
  const words = [
    'Label fewer features: Blacklist or Whitelist under LABEL FILTERING.',
    'Shorten the label text: a priority rule, for example gene before product.',
    'Draw a Linear map.'
  ];
  const { dialog, steps, writes, generates } = dialogHarness();
  const state = dialog.labelFitDialog;
  await dialog.openLabelFitDialog(normalizeReport({ ...report([notFound('outward_both'), notFound('smaller_font')]),
    remedies: ways }), { mode: 'circular' });
  assert.equal(state.capacity, false);
  assert.deepEqual(state.rows.map((item) => [item.kind, item.status, item.offered]),
    [['outward_both', 'not_found', false], ['smaller_font', 'not_found', false]]);
  assert.deepEqual(state.remedies, words);
  assert.equal(state.selectedKind, '');
  await dialog.resolveLabelFitChoice('cancel');
  assert.equal(state.show, false);
  assert.deepEqual([steps, writes, generates], [[], [], []]);

  // Rows still to compute may fit: the ways show only once every row is known not to fit.
  const later = [];
  const computeRow = (token, kind) => new Promise((resolve) => later.push(() => resolve({ row: notFound(kind) })));
  const lazy = dialogHarness({ computeRow });
  const opened = lazy.dialog.openLabelFitDialog(normalizeReport({
    ...report([notFound('outward_both'), row('combination', 'not_computed')], { drawn: true }), remedies: ways
  }), { rowToken: 'circular:1:0', mode: 'circular' });
  assert.deepEqual(lazy.dialog.labelFitDialog.remedies, []);
  await Promise.resolve();
  later.splice(0).forEach((resolve) => resolve());
  await opened;
  assert.deepEqual(lazy.dialog.labelFitDialog.rows.map((item) => item.status), ['not_found', 'not_found']);
  assert.deepEqual(lazy.dialog.labelFitDialog.remedies, words);

  // A row that fits: no ways, even when the report names them.
  await dialog.openLabelFitDialog(normalizeReport({ ...report([notFound('outward_both'), foundRow('outward_vertical', 1.5)]),
    remedies: ways }), { mode: 'circular' });
  assert.deepEqual(state.remedies, []);
});

// The status line the shared choice dialog shows, from the label fit dialog's
// bindings and the choice dialog's template, for the harness's state.
const choiceDialogStatus = async () => {
  const { readFile } = await import('node:fs/promises');
  const html = await readFile(new URL('../../gbdraw/web/index.html', import.meta.url), 'utf8');
  const open = html.match(/<choice-dialog v-if="labelFitDialog\.show"[^>]*>/s)?.[0] || '';
  const busy = open.match(/ :busy="([^"]*)"/)?.[1] || 'false';
  const busyStatus = open.match(/ :busy-status="([^"]*)"/)?.[1] || 'null';
  const start = html.indexOf('id="choice-dialog-template"');
  const template = html.slice(start, html.indexOf('</script>', start));
  const lines = [...template.matchAll(/<p v-if="busy" role="status"[^>]*>\{\{ (.+?) \}\}<\/p>/g)].map((match) => match[1]);
  assert.equal(lines.length, 1, 'the choice dialog has one busy status line');
  const bindings = new Function('dialogChoicePending', 'labelFitDialog', 't', `return [${busy}, ${busyStatus}];`);
  const text = new Function('busyStatus', 't', `return ${lines[0]};`);
  return (/** @type {boolean} */ pending, /** @type {object} */ dialog) => {
    const [isBusy, status] = bindings(pending, dialog, t);
    return isBusy ? [text(status, t)] : [];
  };
};

test('Apply or Cancel during a row waits for that row, says so, and starts no further row (LP-D60)', async () => {
  const { readFile } = await import('node:fs/promises');
  const statusShown = await choiceDialogStatus();
  for (const choice of /** @type {const} */ (['cancel', 'apply'])) {
    const asked = [];
    let release;
    const computeRow = (token, kind) => {
      asked.push(kind);
      return new Promise((resolve) => { release = () => resolve({ row: row(kind, 'not_found', { reason: 'No value found within 30 seconds.' }) }); });
    };
    const rows = [row('outward_vertical', 'not_computed'), row('combination', 'not_computed'), foundRow('smaller_font', 1.2)];
    const applying = [];
    let harness;
    harness = dialogHarness({ computeRow, onWrite: () => {
      applying.push(statusShown(harness.pending.value, harness.dialog.labelFitDialog));
    } });
    const { dialog, steps, writes, generates, pending } = harness;
    const computing = dialog.openLabelFitDialog(normalizeReport(report(rows, { drawn: true })),
      { rowToken: 'circular:1:0', mode: 'circular' });
    assert.equal(dialog.labelFitDialog.rows[0].status, 'computing');
    assert.deepEqual(statusShown(pending.value, dialog.labelFitDialog), []);
    const chosen = dialog.resolveLabelFitChoice(choice);
    // The dialog stays open and busy, says only that it is finishing the row, and takes no second choice.
    assert.equal(dialog.labelFitDialog.finishing, true);
    assert.equal(dialog.labelFitDialog.show, true);
    assert.deepEqual(statusShown(pending.value, dialog.labelFitDialog), ['Finishing the choice being computed…'], choice);
    assert.equal(await dialog.resolveLabelFitChoice('apply'), undefined);
    assert.deepEqual([steps.length, writes.length, generates.length], [0, 0, 0]);
    release();
    await chosen;
    await computing;
    assert.deepEqual(asked, ['outward_vertical'], choice);
    assert.equal(dialog.labelFitDialog.finishing, false);
    assert.equal(dialog.labelFitDialog.show, false);
    // Apply: one History step and one Generate, "Applying an edit…" while it writes; Cancel: none (PD-OI-088).
    assert.deepEqual([steps.length, writes.length, generates.length], choice === 'apply' ? [1, 1, 1] : [0, 0, 0]);
    assert.deepEqual(applying, choice === 'apply' ? [['Applying an edit…']] : []);
  }
  // The label fit dialog's words; no status line of its own beside the shared one.
  const html = await readFile(new URL('../../gbdraw/web/index.html', import.meta.url), 'utf8');
  const body = html.slice(html.indexOf('<choice-dialog v-if="labelFitDialog.show"'), html.indexOf('</choice-dialog>', html.indexOf('data-label-fit-dialog')));
  assert.equal(body.match(/Finishing the choice being computed…/g)?.length, 1);
  assert.doesNotMatch(body.slice(body.indexOf('>') + 1), /Finishing/);
  const source = await readFile(new URL('../../gbdraw/web/js/app/label-fit-dialog.js', import.meta.url), 'utf8');
  assert.doesNotMatch(source, /Finishing/);
});

test('an exported interactive SVG opens the annotation of a region label drawn in front of the leader lines (OV-1289)', async () => {
  const { STANDALONE_INTERACTIVE_SCRIPT: script } = await import(
    '../../gbdraw/web/js/services/standalone-interactivity-assets.js'
  );
  const start = script.indexOf('function closestAnnotation(');
  const end = script.indexOf('\n  }\n', start);
  assert.ok(start >= 0 && end > start);
  const selector = '[data-gbdraw-interactive-annotation="true"][data-gbdraw-annotation-id]';
  const element = (attributes, parentNode, annotation = false) => ({
    parentNode,
    getAttribute: (name) => attributes[name] ?? null,
    matches: (query) => annotation && query === selector
  });
  const byId = new Map();
  const svg = {
    ownerDocument: { getElementById: (id) => byId.get(id) ?? null },
    contains: (node) => [...byId.values()].includes(node)
  };
  const record = element({}, svg);
  const annotationGroup = element({ id: 'annotation-r0' }, record, true);
  byId.set('annotation-r0', annotationGroup);
  const frontLayer = element({}, record);
  const frontItem = element({ 'data-gbdraw-region-label-for': 'annotation-r0' }, frontLayer);
  const frontText = element({}, frontItem);
  const stray = element({ 'data-gbdraw-region-label-for': 'missing' }, frontLayer);
  const closestAnnotation = new Function('svg', 'ANNOTATION_SELECTOR',
    `${script.slice(start, end + 4)}\nreturn closestAnnotation;`)(svg, selector);

  assert.equal(closestAnnotation(frontText), annotationGroup);
  assert.equal(closestAnnotation(element({}, annotationGroup)), annotationGroup);
  assert.equal(closestAnnotation(element({}, stray)), null);
});
