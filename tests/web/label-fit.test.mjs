// The label fit report Python sends, the draft fields a row writes, and the
// dialog that shows the rows (LP-D43, LP-D45, LP-D48, OD-F14).
import assert from 'node:assert/strict';
import test from 'node:test';

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
const { normalizeLabelFitFailure, normalizeLabelFitRow, normalizeResultLabelFits } = await import(
  '../../gbdraw/web/js/services/label-fit.js'
);
const { createLabelFitDialog } = await import('../../gbdraw/web/js/app/label-fit-dialog.js');

const Y_OFFSET = 'labels.unified_adjustment.outer_labels.y_radius_offset';
const row = (kind, status, extra = {}) => ({
  kind, title: kind.replace(/_/g, ' '), status, values: null, settings: [], canvasSize: null, reason: null, note: null,
  ...extra
});
const foundRow = (kind, value) => row(kind, 'found', {
  values: { [Y_OFFSET]: value }, settings: [{ name: 'Outer Y Offset', unit: '', value }],
  canvasSize: { width: 1708, height: 1130 }
});
const report = (rows, { drawn = false } = {}) => ({
  schema: 1, drawn, overlapPairs: drawn ? 0 : 3, crossingPairs: 5, definitionOverlaps: 0, regionOverlaps: 0,
  notice: drawn ? ['Leader lines cross (5 pairs).'] : [],
  records: [{ recordId: 'NC_1', labelKind: 'circular_outer', labelKindTitle: 'outer labels',
    labelCount: 80, fontSize: 18, overlapPairs: drawn ? 0 : 3, crossingPairs: 5, definitionOverlaps: 0,
    regionOverlaps: 0, refused: false, counts: '3 label pairs overlap, 5 leader pairs cross' }],
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
  records: [{ recordId: 'adjacent', labelKind: 'linear_region', labelKindTitle: 'region label', labelCount: 1,
    fontSize: 10, overlapPairs: 0, crossingPairs: 0, definitionOverlaps: 0, regionOverlaps: 1, refused: false,
    counts: '1 region label overlap' }]
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
  // "Separate overlapping features" (LP-D52) switches the mode's Resolve Overlaps control on.
  assert.deepEqual(draft({ 'canvas.resolve_overlaps': true }), { form: {}, adv: { resolve_overlaps: true } });
  assert.deepEqual(draftFieldsForConfigValues({ 'canvas.resolve_overlaps': true }, 'linear', state.drawings.linear.adv),
    { form: {}, adv: { resolve_overlaps: true } });
  const unknown = (error) => error.code === 'INPUT_INVALID' && error.context.reason === 'UNKNOWN_CONFIG_PATH';
  assert.throws(() => draftFieldsForConfigValues({ 'labels.unknown': 1 }, 'circular', adv), unknown);
  assert.throws(() => draftFieldsForConfigValues({}, 'circular', adv), unknown);
});

test('a row may write exactly the config paths Python\'s rows send, in the mode that sends them', async () => {
  const { readFile } = await import('node:fs/promises');
  const { applyLabelFitValues } = await import('../../gbdraw/web/js/services/config.js');
  const fit = await readFile(new URL('../../gbdraw/labels/fit.py', import.meta.url), 'utf8');
  const constants = new Map([...fit.matchAll(/^([A-Z_]+_PATHS?) = (.+)$/gm)]
    .map(([, name, value]) => [name, [...value.matchAll(/"([^"]+)"/g)].map(([, path]) => path)]));
  const { SEPARATE_FEATURES_PATH: [separate], FONT_SIZE_PATHS: fontSizes, ...single } = Object.fromEntries(constants);
  assert.equal(fontSizes.length, 2);
  const sample = { [separate]: true, 'labels.circular.scope': 'outer', 'labels.circular.placement': 'radial' };
  const unknown = (error) => error.code === 'INPUT_INVALID' && error.context.reason === 'UNKNOWN_CONFIG_PATH';
  // Circular rows: every path Python names; Linear rows: separate overlapping features only.
  applyLabelFitValues(state.drawings.circular, Object.fromEntries(fontSizes.map((path) => [path, 12])), 'circular');
  for (const [path] of Object.values(single)) {
    applyLabelFitValues(state.drawings.circular, { [path]: sample[path] ?? 1.1 }, 'circular');
  }
  applyLabelFitValues(state.drawings.linear, { [separate]: true }, 'linear');
  // Label settings no row sends are not row writes.
  for (const [values, mode] of [
    [{ 'labels.unified_adjustment.inner_labels.x_radius_offset': 1.1 }, 'circular'],
    [{ 'labels.unified_adjustment.inner_labels.y_radius_offset': 1.1 }, 'circular'],
    [{ 'labels.spacing.linear': 0.5 }, 'linear'],
    [{ 'labels.linear.scope': 'out' }, 'linear']
  ]) {
    assert.throws(() => applyLabelFitValues(state.drawings[mode], values, mode), unknown, JSON.stringify(values));
  }
});

const dialogHarness = ({ computeRow = async () => ({}), onWrite = () => {} } = {}) => {
  const steps = [];
  const writes = [];
  const generates = [];
  const pending = { value: false };
  const closes = [];
  const dialog = createLabelFitDialog({
    reactive: (value) => value,
    computeRow,
    writeRowSettings: (values, mode) => { writes.push({ values, mode }); onWrite(); },
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
  return { dialog, steps, writes, generates, pending };
};

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
  assert.deepEqual(writes, [{ values: { 'canvas.resolve_overlaps': true }, mode: 'circular' }]);
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
  const bindings = new Function('dialogChoicePending', 'labelFitDialog', `return [${busy}, ${busyStatus}];`);
  const text = new Function('busyStatus', 't', `return ${lines[0]};`);
  return (/** @type {boolean} */ pending, /** @type {object} */ dialog) => {
    const [isBusy, status] = bindings(pending, dialog);
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
