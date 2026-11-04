// The Web's Wrap Rows: the words of its refusals and of the export limit (one wording owner,
// utils/error-normalization.js), the form edit that switches what rows decide (services/linear-wrap-form.js
// through the layout transition), its blockers, and the length field's amount and unit. The value's request
// and Session channel is tested in session-request, mode-scoped-session and session50-load-resave.
import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import test from 'node:test';

import { normalizeUserFacingError, reasonText } from '../../gbdraw/web/js/utils/error-normalization.js';
import { createDefaultAdv, createDefaultForm } from '../../gbdraw/web/js/services/session-active-config-contract.js';
import {
  LINEAR_WRAP_FORM_FIELDS, exportLinearWrap, linearWrapAmountText, linearWrapBlockers, linearWrapFormPatch, linearWrapInputValue,
  linearWrapUnitOf, requestLinearWrapLeaf, requestLinearWrapLength, shownLinearWrap
} from '../../gbdraw/web/js/services/linear-wrap-form.js';
import {
  createFeaturePlacementActions, draftFeatureSlot, draftPlacementTargets
} from '../../gbdraw/web/js/app/feature-editor/placement-actions.js';
import { createLinearWrapField } from '../../gbdraw/web/js/app/linear-wrap-field.js';
import { movesARecordSideways } from '../../gbdraw/web/js/services/linear-record-layout.js';
import { createAutoValueDisplay } from '../../gbdraw/web/js/app/auto-value-display.js';
import { withDrawings } from './helpers/drawing-state.mjs';

// What Python sends for each wrap refusal (gbdraw/layout/linear_wrap.py; tests/test_web_wrap_diagnostics.py).
const WRAP_REFUSALS = Object.freeze({
  WRAP: ['COMPARISON_INPUT', 'Remove the comparison inputs (BLAST, comparison pairs, LOSAT, Similarity groups), or turn Wrap Rows off.'],
  WRAP_TRACK_LAYOUT: ['INPUT_INVALID', 'Set Track Layout to Above or Below, move the feature track off the axis in Custom Track Slots, or turn Wrap Rows off.'],
  WRAP_SCALE_STYLE: ['INPUT_INVALID', 'Set Scale Style to Ruler (Ticks), or turn Wrap Rows off.'],
  WRAP_NORMALIZE_LENGTH: ['INPUT_INVALID', 'Turn Normalize Record Lengths off, or turn Wrap Rows off.'],
  WRAP_GRID_ROW: ['INPUT_INVALID', 'Give each record its own row, or turn Wrap Rows off.'],
  WRAP_ANCHOR_FEATURE: ['INPUT_INVALID', 'Turn Wrap Rows off, then use Align at feature.'],
  WRAP_RECORD_POSITION: ['INPUT_INVALID', 'Undo the move, or turn Wrap Rows off.']
});

test('each wrap refusal names what rows cannot draw and the way out, the text its hint shows', () => {
  for (const [reason, [code, way]] of Object.entries(WRAP_REFUSALS)) {
    const text = reasonText(reason);
    assert.ok(text?.endsWith(way), reason);
    assert.doesNotMatch(text, /Preserved session settings|canvas\.linear\.wrap_length/, `${reason}: the retired advice`);
    const model = normalizeUserFacingError({ code, operation: 'generate', stage: 'render', context: { field: 'wrap_length', reason } });
    assert.deepEqual(model.context, { field: 'wrap_length', reason }, reason);
    assert.ok(model.summary.endsWith(`Field: Wrap Rows. ${text}`), `${reason}: ${model.summary}`);
  }
  assert.equal(reasonText('NOT_A_REASON'), null);
  assert.equal(reasonText('toString'), null, 'only the owner\'s own reasons');
});

test('an export over the format limit names the limit, the size and what would fit, and offers SVG', () => {
  const png = normalizeUserFacingError({ code: 'EXPORT_LIMIT', operation: 'export-png', stage: 'export-capture',
    context: { reason: 'PNG_SIZE', limit: 32767, width: 6250, height: 143647, rowLength: 500000, dpi: 68 } });
  assert.equal(png.summary, 'The figure is larger than this file format allows, so no file was made.'
    + ' PNG is limited to 32,767 px per side, and this figure is 6,250 x 143,647 px.'
    + ' A row length of 500,000 bp should fit (estimated from the drawn height).'
    + ' A DPI of 68 or lower would fit.'
    + ' Use SVG, which has no size limit.');
  assert.deepEqual(png.actions, ['use-svg']);
  const pdf = normalizeUserFacingError({ code: 'EXPORT_LIMIT', operation: 'export-pdf', stage: 'export-capture',
    context: { reason: 'PDF_SIZE', limit: 14400, width: 1500, height: 34475, rowLength: 200000 } });
  assert.equal(pdf.summary, 'The figure is larger than this file format allows, so no file was made.'
    + ' PDF is limited to 14,400 pt per side, and this figure is 1,500 x 34,475 pt.'
    + ' A row length of 200,000 bp should fit (estimated from the drawn height).'
    + ' Use SVG, which has no size limit.');
  // D-46: a row length or DPI is advised only when Python found one that fits. An unwrapped PDF (Wrap Rows off),
  // and a wrapped one no longer row length fits, name the limit and SVG only.
  for (const context of [{ limit: 14400, width: 1500, height: 20000 }, { limit: 14400, width: 20000, height: 900 }]) {
    const plain = normalizeUserFacingError({ code: 'EXPORT_LIMIT', operation: 'export-pdf', stage: 'export-capture',
      context: { reason: 'PDF_SIZE', ...context } });
    assert.equal(plain.summary, 'The figure is larger than this file format allows, so no file was made.'
      + ` PDF is limited to 14,400 pt per side, and this figure is ${context.width.toLocaleString('en-US')} x ${context.height.toLocaleString('en-US')} pt.`
      + ' Use SVG, which has no size limit.');
    assert.doesNotMatch(plain.summary, /row length|DPI/);
  }
  const unbounded = normalizeUserFacingError({ code: 'EXPORT_LIMIT', operation: 'export-png', stage: 'export-capture',
    context: { reason: 'PNG_SIZE', limit: 32767, width: -1, height: 1.5, dpi: 2e7 } });
  assert.deepEqual(unbounded.context, { reason: 'PNG_SIZE', limit: 32767 }, 'only bounded counts are kept');
  assert.equal(unbounded.summary, 'The figure is larger than this file format allows, so no file was made.'
    + ' Use SVG, which has no size limit.');
});

// The fields Wrap Rows switches, as a drawing's form holds them.
const switched = (form) => Object.fromEntries(LINEAR_WRAP_FORM_FIELDS.map((field) => [field, form[field]]));
const unwrappedForm = (overrides = {}) => ({ ...createDefaultForm(), linear_track_layout: 'middle', scale_style: 'bar',
  normalize_length: true, linear_ruler_on_axis: false, ...overrides });

test('turning Wrap Rows on switches what rows decide; turning it off writes only the length (case 1)', () => {
  const on = { linear_wrap_length: 'auto', linear_track_layout: 'above', scale_style: 'ruler', normalize_length: false,
    linear_ruler_on_axis: true };
  assert.deepEqual(linearWrapFormPatch(unwrappedForm(), 'auto'), on);
  // Above and Below are kept; a length switches the same fields (also when no record would be cut).
  assert.deepEqual(linearWrapFormPatch(unwrappedForm({ linear_track_layout: 'below' }), 20000),
    { ...on, linear_wrap_length: 20000, linear_track_layout: 'below' });
  assert.deepEqual(linearWrapFormPatch(unwrappedForm({ linear_track_layout: 'above' }), 'auto'), on);
  assert.deepEqual(linearWrapFormPatch(unwrappedForm(), null), { linear_wrap_length: null });
  // Every field a patch writes is a field the layout transition saves and restores.
  for (const value of ['auto', 20000, null]) {
    assert.deepEqual(Object.keys(linearWrapFormPatch(unwrappedForm(), value)).filter((key) => !LINEAR_WRAP_FORM_FIELDS.includes(key)), []);
  }
});

// A Linear drawing in simple mode, Features on axis and strands together: the feature slot splits into
// an Above and a Below lane.
const linearDrawing = () => ({ id: 'linear', mode: 'linear', form: unwrappedForm({ separate_strands: false }),
  adv: { ...createDefaultAdv('linear'), linear_track_slots_enabled: false }, featurePlacementOverrides: {},
  featureCatalog: { value: null } });
const placementStore = (drawing) => Object.assign(withDrawings(drawing), { mode: { value: 'linear' } });
const checkbox = (checked) => ({ target: { type: 'checkbox', checked, getAttribute: () => 'Wrap Rows' } });

test('Wrap Rows is one edit of the layout transition: Cancel restores every field, Reset is one step (case 2)', async () => {
  // No lane placement: the edit applies at once, all fields together, and the control's History adapter records
  // it as one step (no runUndoable step of its own).
  const plain = linearDrawing();
  const steps = [];
  const actions = createFeaturePlacementActions({ state: placementStore(plain), getCommittedRequest: () => null,
    isCurrentFeature: () => true, runUndoable: async (label, fn) => { steps.push(label); fn(); } });
  const before = switched(plain.form);
  assert.equal(actions.changeLayoutSetting(checkbox(true), 'linear_wrap_length', 'auto'), undefined);
  assert.deepEqual(switched(plain.form), linearWrapFormPatch(unwrappedForm(), 'auto'));
  assert.deepEqual(steps, []);
  assert.notDeepEqual(switched(plain.form), before);

  // A lane placement Features on axis -> Above would lose: the dialog asks, the form keeps every field.
  const laned = linearDrawing();
  laned.featurePlacementOverrides.l = { recordKey: 'r', biologicalFeatureId: 'cds', placement: { kind: 'lane', side: 'above', level: 1 } };
  assert.deepEqual(draftPlacementTargets(laned).map((target) => target.side ?? target.kind), ['main', 'above', 'below']);
  const ask = createFeaturePlacementActions({ state: placementStore(laned), getCommittedRequest: () => null,
    isCurrentFeature: () => true, runUndoable: async (label, fn) => { steps.push(label); fn(); } });
  const lanedBefore = JSON.stringify(laned.form);
  const control = checkbox(true);
  assert.equal(ask.changeLayoutSetting(control, 'linear_wrap_length', 'auto'), false);
  assert.equal(ask.layoutChange.open, true);
  assert.equal(JSON.stringify(laned.form), lanedBefore, 'Cancel restores every field the edit wrote');
  assert.equal(control.target.checked, false, 'the checkbox shows the kept value');
  await ask.resolveLayoutChange('cancel');
  assert.equal(JSON.stringify(laned.form), lanedBefore);
  assert.deepEqual(steps, []);
  ask.changeLayoutSetting(checkbox(true), 'linear_wrap_length', 'auto');
  await ask.resolveLayoutChange('reset');
  assert.deepEqual(steps, ['Change setting and reset Feature placements']);
  assert.deepEqual(switched(laned.form), linearWrapFormPatch(unwrappedForm(), 'auto'));
  assert.deepEqual(laned.featurePlacementOverrides, {});
});

// PD-OI-113 scenario revision 2: one rule for every source of a record's place. The blocker reads the stored
// translations (Align at feature, a Similarity alignment, a loaded CLI or Python Session alike); only x moves a
// record sideways, and whether a record is cut is Python's (WRAP_RECORD_POSITION).
// The x test repeats Python's `translation.x != 0` (gbdraw/layout/linear_wrap.py `refuse_undrawable_wrap`): a twin with
// no generator, held by these cases and tests/test_linear_wrap_plan.py's x 0 / x 40 pair.
test('a drawn record moved sideways is a translation the request carries with an x other than 0, whatever wrote it', () => {
  const drawn = [{ uid: 'record-1' }, { uid: 'record-2' }];
  assert.equal(movesARecordSideways([], drawn), false);
  assert.equal(movesARecordSideways([{ recordKey: 'record-1', x: 0, y: 0 }, { recordKey: 'record-2', x: 0, y: 12 }], drawn), false,
    'zero entries (a LOSAT request) and a vertical move only');
  assert.equal(movesARecordSideways([{ recordKey: 'record-1', x: 0, y: 0 }, { recordKey: 'record-2', x: -62.8, y: 0 }], drawn), true,
    'Align at feature on HmmtDNA and NC_002333.2');
  assert.equal(movesARecordSideways([{ recordKey: 'cli-record-1', x: 15, y: 0 }], [{ uid: 'cli-record-1' }]), true,
    'a loaded --anchor_feature Session');
  // Review 1 L1: the entry of a record no longer drawn (its File removed, the record turned off) is not read, as the
  // request leaves it out (`recordTranslationsFor`), so it blocks nothing Python would draw.
  assert.equal(movesARecordSideways([{ recordKey: 'record-1', x: 0, y: 0 }, { recordKey: 'record-2', x: -62.8, y: 0 }],
    [{ uid: 'record-1' }]), false);
});

test('each blocker keeps Wrap Rows from turning on, in its Generate refusal\'s words; Simple-mode Middle is none (case 3)', () => {
  const none = { comparisons: false, sharedRecordRow: false, featureAlignment: false, movedRecords: false, featureTrackOnAxis: false };
  assert.deepEqual(linearWrapBlockers(none), []);
  assert.deepEqual(linearWrapBlockers({ ...none, comparisons: true }), ['WRAP']);
  assert.deepEqual(linearWrapBlockers({ ...none, sharedRecordRow: true }), ['WRAP_GRID_ROW']);
  assert.deepEqual(linearWrapBlockers({ ...none, movedRecords: true }), ['WRAP_RECORD_POSITION']);
  // Review 1 M1 (D-65): Align at feature blocks also when its Apply moved no record; with a move, one hint.
  assert.deepEqual(linearWrapBlockers({ ...none, featureAlignment: true }), ['WRAP_RECORD_POSITION']);
  assert.deepEqual(linearWrapBlockers({ ...none, featureAlignment: true, movedRecords: true }), ['WRAP_RECORD_POSITION']);
  assert.deepEqual(linearWrapBlockers({ ...none, featureTrackOnAxis: true }), ['WRAP_TRACK_LAYOUT']);
  assert.deepEqual(linearWrapBlockers({ comparisons: true, sharedRecordRow: true, featureAlignment: true, movedRecords: true, featureTrackOnAxis: true }),
    ['WRAP', 'WRAP_GRID_ROW', 'WRAP_RECORD_POSITION', 'WRAP_TRACK_LAYOUT']);
  // The feature slot the blocker reads is the one the placement targets read: in simple mode, Middle puts the
  // features on the axis, which turning wrap on switches (no blocker); a custom stack's overlay slot is one.
  const drawing = linearDrawing();
  assert.equal(draftFeatureSlot(drawing).side, 'overlay');

  const store = withDrawings(linearDrawing());
  let inputs = { ...none, comparisons: true };
  const writes = [];
  const field = createLinearWrapField({ state: store, blockerInputs: () => inputs, watch: () => {}, shownWrap: () => ({ kind: 'off' }),
    changeLayoutSetting: (_event, fieldName, value) => { writes.push([fieldName, value]); } });
  assert.equal(field.linearWrapBlocked.value, true);
  assert.equal(field.linearWrapBlockedText.value, reasonText('WRAP'));
  inputs = { ...none, sharedRecordRow: true, featureTrackOnAxis: true };
  assert.equal(field.linearWrapBlockedText.value, `${reasonText('WRAP_GRID_ROW')} ${reasonText('WRAP_TRACK_LAYOUT')}`);
  // Turning it off is never blocked: the checkbox writes Off.
  field.changeLinearWrap(checkbox(false), 'on');
  assert.deepEqual(writes, [['linear_wrap_length', null]]);
});

test('the length field writes bp from an amount and a unit, and shows a value from elsewhere in its largest whole unit', () => {
  assert.equal(linearWrapInputValue('', 'kbp'), 'auto', 'empty is Auto');
  assert.equal(linearWrapInputValue('  ', 'Mbp'), 'auto');
  assert.equal(linearWrapInputValue('20', 'kbp'), 20000);
  assert.equal(linearWrapInputValue('1.1', 'kbp'), 1100);
  assert.equal(linearWrapInputValue('0.5', 'Mbp'), 500000);
  assert.equal(linearWrapInputValue('1.5', 'bp'), 1.5, 'Python refuses a length that is no whole bp');
  // The decimal amount is shifted by the unit exactly (D-57), not rounded: a non-whole bp is sent as typed.
  assert.equal(linearWrapInputValue('1.0000000001', 'kbp'), 1000.0000001, 'the CLI refuses --wrap_length 1.0000000001kb too');
  assert.equal(linearWrapInputValue('0.0000015', 'Mbp'), 1.5);
  assert.equal(linearWrapInputValue('1.10', 'kbp'), 1100);
  assert.equal(linearWrapInputValue('.5', 'kbp'), 500);
  assert.equal(linearWrapInputValue('20.', 'kbp'), 20000);
  assert.equal(linearWrapInputValue('0.000001', 'Mbp'), 1);
  assert.equal(linearWrapUnitOf(20000, 'bp'), 'kbp');
  assert.equal(linearWrapUnitOf(3_000_000, 'bp'), 'Mbp');
  assert.equal(linearWrapUnitOf(1500, 'Mbp'), 'bp');
  assert.equal(linearWrapUnitOf('auto', 'Mbp'), 'Mbp');
  assert.equal(linearWrapAmountText(20000, 'kbp'), '20');
  assert.equal(linearWrapAmountText('auto', 'kbp'), '');

  const drawing = linearDrawing();
  const store = withDrawings(drawing);
  const watchers = [];
  const field = createLinearWrapField({ state: store, watch: (_source, callback) => watchers.push(callback), shownWrap: () => ({ kind: 'unknown' }),
    blockerInputs: () => ({ comparisons: false, sharedRecordRow: false, featureAlignment: false, movedRecords: false, featureTrackOnAxis: false }),
    changeLayoutSetting: (_event, _field, value) => { drawing.form.linear_wrap_length = value; } });
  const changed = () => watchers.forEach((callback) => callback());
  field.changeLinearWrap({ target: { value: '20' } }, 'amount');
  changed();
  assert.equal(drawing.form.linear_wrap_length, 20000);
  assert.equal(field.linearWrapAmount.value, '20');
  // A unit change keeps the typed amount.
  field.changeLinearWrap({ target: { value: 'Mbp' } }, 'unit');
  changed();
  assert.equal(drawing.form.linear_wrap_length, 20_000_000);
  assert.equal(field.linearWrapUnit.value, 'Mbp');
  // An Undo or a Load writes 1500: it shows as 1500 bp.
  drawing.form.linear_wrap_length = 1500;
  changed();
  assert.equal(field.linearWrapUnit.value, 'bp');
  assert.equal(field.linearWrapAmount.value, '1500');
  field.changeLinearWrap({ target: { value: '' } }, 'amount');
  changed();
  assert.equal(drawing.form.linear_wrap_length, 'auto');
  assert.equal(field.linearWrapUnit.value, 'bp', 'Auto keeps the unit');
});

// Case 6 (summary part): what the shown Result drew with the drawing's value comes from Python's input summary
// and the committed request; the Web copies neither the Auto rule nor the drawn axis. "These records" compares each
// drawn record's source length and requested region, as Python sent them, with the card's (review 3 F2, D-53).
const LAMBDA = Object.freeze({ sourceLength: 48502, region: null });
const SUMMARY = Object.freeze({ records: [LAMBDA],
  linearWrap: { rowLength: 10000, rowCount: 5, requested: 'auto', recordRows: [5], drawnAxisLength: 10000 } });
const ONE_ROW = Object.freeze({ records: [{ sourceLength: 8000, region: null }], linearWrap: null });
/** @param {number} sourceLength @param {{ start: number, end: number } | null} [region] */
const card = (sourceLength, region = null) => ({ sourceLength, region });

test('the shown Result names its rows only when it drew this value over these records (case 6)', () => {
  const read = (inputs) => shownLinearWrap({ value: 'auto', requested: 'auto', summary: SUMMARY, records: [card(48502)], ...inputs });
  assert.deepEqual(read({}), { kind: 'rows', rowLength: 10000, rowCount: 5, drawnAxisLength: 10000 });
  assert.deepEqual(read({ value: null }), { kind: 'off' });
  // One row: the Result cut nothing over these records. A record added since (no Generate) is unknown: the next
  // Generate may cut it, so neither the text nor the Auto hints claim one row (review 2 F2).
  assert.deepEqual(read({ summary: ONE_ROW, records: [card(8000)] }), { kind: 'one-row' }, 'Auto with nothing cut');
  assert.deepEqual(read({ summary: ONE_ROW, records: [card(8000), card(5_000_000)] }), { kind: 'unknown' }, 'a record added since');
  assert.deepEqual(read({ summary: ONE_ROW, records: [card(5_000_000)] }), { kind: 'unknown' }, 'a record replaced since');
  assert.deepEqual(read({ summary: ONE_ROW, records: [] }), { kind: 'unknown' }, 'records not known');
  assert.deepEqual(read({ value: 20000, requested: 20000, summary: ONE_ROW, records: [card(8000), card(30000)] }), { kind: 'unknown' },
    'a typed length over the old longest record');
  for (const changed of [{ requested: undefined }, { requested: 20000 }, { summary: undefined }, { records: [card(48502), card(1000)] },
    { records: [] }, { value: 20000, requested: 20000 }]) {
    assert.deepEqual(read(changed), { kind: 'unknown' }, JSON.stringify(changed));
  }
  // The committed request's value: its override leaf, or a CLI request's full config; the leaf as stored for the
  // export check, which Python reads (review 3 F1, D-52).
  assert.equal(requestLinearWrapLength({ diagramOptions: { configOverrides: { 'canvas.linear.wrap_length': 'auto' } } }), 'auto');
  assert.equal(requestLinearWrapLength({ diagramOptions: { config: { canvas: { linear: { wrap_length: 20000 } } }, configOverrides: {} } }), 20000);
  assert.equal(requestLinearWrapLength({ diagramOptions: { configOverrides: {} } }), null);
  assert.equal(requestLinearWrapLeaf({ diagramOptions: { config: { canvas: { linear: { wrap_length: '5kb' } } }, configOverrides: {} } }), '5kb');
  assert.equal(requestLinearWrapLeaf({ diagramOptions: { configOverrides: {} } }), null);

  const drawing = linearDrawing();
  drawing.form.linear_wrap_length = 'auto';
  let reading = read({});
  const field = createLinearWrapField({ state: withDrawings(drawing), watch: () => {}, changeLayoutSetting: () => {},
    shownWrap: () => reading, blockerInputs: () => ({ comparisons: false, sharedRecordRow: false, featureAlignment: false, movedRecords: false, featureTrackOnAxis: false }) });
  assert.equal(field.linearWrapAutoText.value, '10 kbp per row, 5 rows (auto)');
  reading = { kind: 'one-row' };
  assert.equal(field.linearWrapAutoText.value, 'One row (auto)');
  reading = { kind: 'unknown' };
  assert.equal(field.linearWrapAutoText.value, 'Auto');
});

// Review 3 F2 (D-53): a cropped record. Python drew 200,000 bases of a 5 Mb record in 4 rows of 50 kbp; the card
// still asks for the same region of a record of the same length, so the field shows the rows, and a cropped record
// that is not cut shows one row. A changed region, or a source of another length, is unknown until Generate.
test('a cropped record shows the rows and one row Python drew while the card asks for the same region (review 3 F2)', () => {
  const region = { start: 1, end: 200_000 };
  const cropped = { records: [{ sourceLength: 5_000_000, region }],
    linearWrap: { rowLength: 50_000, rowCount: 4, requested: 'auto', recordRows: [4], drawnAxisLength: 50_000 } };
  const read = (inputs) => shownLinearWrap({ value: 'auto', requested: 'auto', summary: cropped, records: [card(5_000_000, region)], ...inputs });
  assert.deepEqual(read({}), { kind: 'rows', rowLength: 50_000, rowCount: 4, drawnAxisLength: 50_000 });
  const uncut = { records: [{ sourceLength: 5_000_000, region: { start: 1001, end: 9000 } }], linearWrap: null };
  assert.deepEqual(read({ summary: uncut, records: [card(5_000_000, { start: 1001, end: 9000 })] }), { kind: 'one-row' });
  assert.deepEqual(read({ records: [card(5_000_000, { start: 1, end: 300_000 })] }), { kind: 'unknown' }, 'another region');
  assert.deepEqual(read({ records: [card(5_000_000)] }), { kind: 'unknown' }, 'the region removed');
  assert.deepEqual(read({ records: [card(4_000_000, region)] }), { kind: 'unknown' }, 'a source of another length');
});

test('while wrapping, the Auto sizes and ticks follow the drawn axis Python sent, and the windows the whole record (case 6)', () => {
  const displayOf = (lengths, axis) => {
    const state = withDrawings({ mode: { value: 'linear' }, form: { scale_style: 'ruler' }, adv: {}, unmanagedConfigOverrides: {} });
    return createAutoValueDisplay(state, { linearRecordLengths: () => lengths, linearWrappedAxisLength: () => axis });
  };
  const wrapped = displayOf([48502], 10000);
  const oneRow = displayOf([10000], undefined);
  const whole = displayOf([48502], undefined);
  for (const key of ['scaleInterval', 'axisStrokeWidth', 'labelFontSize']) assert.equal(wrapped.autoValueText(key), oneRow.autoValueText(key), key);
  assert.notEqual(whole.autoValueText('scaleInterval'), oneRow.autoValueText('scaleInterval'), 'the fixture tells the axis from the record');
  assert.equal(wrapped.autoValueText('dinucleotideWindow'), whole.autoValueText('dinucleotideWindow'));
  // Not known: the length-dependent texts name their rule.
  const unknown = displayOf([48502], null);
  assert.equal(unknown.autoValueText('scaleInterval'), 'by record length (auto)');
  assert.equal(unknown.autoValueText('dinucleotideWindow'), whole.autoValueText('dinucleotideWindow'));
});

// Case 6 (port part): PNG and PDF ask the Worker's size check before the browser draws or downloads anything,
// with the SVG size and the PNG DPI; SVG and Interactive SVG never ask. The page has no DOM here: a minimal
// stand-in records what the export creates, so a refusal that reaches the canvas or the PDF staging fails.
const exportDom = (width, height) => {
  const created = [];
  const svg = {
    attributes: new Map([['width', String(width)], ['height', String(height)]]),
    getAttribute(name) { return this.attributes.get(name) ?? null; },
    setAttribute(name, value) { this.attributes.set(name, String(value)); },
    removeAttribute(name) { this.attributes.delete(name); },
    querySelectorAll: () => [],
    classList: { toggle() {}, remove() {}, length: 0 }
  };
  const saved = { document: globalThis.document, XMLSerializer: globalThis.XMLSerializer };
  globalThis.document = { createElement: (tag) => {
    created.push(tag);
    if (tag === 'template') return { set innerHTML(_markup) {}, content: { querySelector: () => svg } };
    throw new Error(`the export created a ${tag}`);
  } };
  globalThis.XMLSerializer = class { serializeToString() { return '<svg/>'; } };
  return { svg, created, restore: () => Object.assign(globalThis, saved) };
};

test('a PNG or PDF over its limit is refused by the Worker check before anything is drawn (case 6)', async () => {
  const { downloadPNG, downloadPDF, downloadSVG } = await import('../../gbdraw/web/js/services/export.js');
  const dom = exportDom(1500, 45966);
  try {
    const requests = [];
    const refusal = Object.assign(new Error('over the limit'), { code: 'EXPORT_LIMIT' });
    const checkExportLimit = async (request) => { requests.push(request); throw refusal; };
    const snapshot = { name: 'lambda.svg', svg: dom.svg, dpi: '300' };
    await assert.rejects(downloadPNG(snapshot, { checkExportLimit }), refusal);
    await assert.rejects(downloadPDF(snapshot, { checkExportLimit }), refusal);
    assert.deepEqual(requests, [
      { format: 'png', width: 1500, height: 45966, dpi: 300 },
      { format: 'pdf', width: 1500, height: 45966, dpi: null }
    ]);
    assert.deepEqual(dom.created.filter((tag) => tag !== 'template'), [], 'no canvas and no PDF staging');
    // A missing port is a wiring error, not a skipped check.
    await assert.rejects(downloadPNG(snapshot, {}), /size check is not available/);
    // SVG never asks: it has no size limit (the stand-in stops its download link).
    await assert.rejects(downloadSVG(snapshot, { checkExportLimit }), /created a a/);
    assert.equal(requests.length, 2, 'SVG asked no size check');
  } finally {
    dom.restore();
  }
});

// Review 2 F3 (D-56): while Wrap Rows is on, Normalize Record Lengths and Ruler on Axis are fixed only while they
// hold the value rows need, as Track Layout and Scale Style disable only the conflicting option. A loaded Session
// that holds the other value keeps the control editable toward it, so WRAP_NORMALIZE_LENGTH's first remedy works.
test('Normalize Record Lengths and Ruler on Axis are fixed only while they hold what rows need', () => {
  const drawing = linearDrawing();
  const field = createLinearWrapField({ state: withDrawings(drawing), watch: () => {}, changeLayoutSetting: () => {},
    shownWrap: () => ({ kind: 'unknown' }),
    blockerInputs: () => ({ comparisons: false, sharedRecordRow: false, featureAlignment: false, movedRecords: false, featureTrackOnAxis: false }) });
  const fixed = () => [field.linearWrapFixesNormalize.value, field.linearWrapFixesRulerOnAxis.value];
  drawing.form.linear_wrap_length = null;
  drawing.form.normalize_length = false;
  drawing.form.linear_ruler_on_axis = true;
  assert.deepEqual(fixed(), [false, false], 'wrap off: both free');
  // Right after turning on, both hold what rows need (the switch wrote them).
  Object.assign(drawing.form, linearWrapFormPatch(drawing.form, 'auto'));
  assert.deepEqual(fixed(), [true, true]);
  // A loaded `--normalize_length --wrap_length 50kb` Session: Normalize is on, so it stays editable.
  drawing.form.normalize_length = true;
  drawing.form.linear_ruler_on_axis = false;
  assert.deepEqual(fixed(), [false, false]);
});

// Review 2 F1 (D-55): the size check of an export reads the Result's committed wrap length and its records'
// lengths, and Python decides from them whether a record is cut (gbdraw/web_support/export_limit.py). A loaded
// Result has no input summary: the records the Web knows stand in, and a PDF the CLI refuses is refused.
const pythonExportCheck = (requests) => JSON.parse(execFileSync(process.env.PYTHON || 'python', ['-c', `
import json, sys
from gbdraw.config.models.canvas import linear_wrap_length_from_config
from gbdraw.exceptions import GbdrawError
from gbdraw.web_support.export_limit import check_web_export_limit_json
out = []
for request in json.load(sys.stdin):
    try:
        if "wrapLength" in request:
            out.append(linear_wrap_length_from_config(request["wrapLength"]))
        else:
            out.append(json.loads(check_web_export_limit_json(json.dumps(request))))
    except GbdrawError as error:
        out.append(error.diagnostic)
print(json.dumps(out))
`], { input: JSON.stringify(requests), encoding: 'utf8' }));

test('a PDF export of a wrapped Result without an input summary is refused when a record is cut, as on the CLI', async () => {
  // What the export sends: the committed request's wrap leaf as stored, with the records the summary names, or,
  // without a summary, the records the drawing's discovery reads, awaited at export (review 3 F1, D-52): each
  // record's source length and region, cropped or not (D-70); nothing while the records stay unknown, wrap is off,
  // or no Linear Result is shown.
  const whole = (length) => ({ sourceLength: length, region: null });
  const discovered = (lengths) => async () => lengths.map(whole);
  assert.deepEqual(await exportLinearWrap({ leaf: 10000, summary: undefined, discoverRecords: discovered([5_000_000]) }),
    { requested: 10000, records: [whole(5_000_000)] });
  assert.deepEqual(await exportLinearWrap({ leaf: 'auto', summary: SUMMARY, discoverRecords: discovered([]) }),
    { requested: 'auto', records: [whole(48502)] }, 'a Result drawn here: the records Python drew');
  assert.equal(await exportLinearWrap({ leaf: null, summary: undefined, discoverRecords: discovered([5_000_000]) }), null, 'wrap off');
  assert.equal(await exportLinearWrap({ leaf: undefined, summary: undefined, discoverRecords: discovered([5_000_000]) }), null,
    'no Linear request');
  assert.equal(await exportLinearWrap({ leaf: 'auto', summary: undefined, discoverRecords: discovered([]) }), null, 'records not known');

  // The state after a Load: no input summary, and the cards' record lists deferred until something reads them. The
  // export awaits the discovery as Generate does, so the 5 Mb record is known when the check is asked.
  const page = { deferred: true, records: [] };
  const discoverRecords = async () => {
    if (page.deferred) { await Promise.resolve(); page.deferred = false; page.records = [whole(5_000_000)]; }
    return page.records;
  };
  const loaded = (leaf) => exportLinearWrap({ leaf, summary: undefined, discoverRecords });
  const pdf = (linearWrap) => ({ format: 'pdf', width: 1500, height: 34_000, dpi: null, linearWrap });
  const [cut, unit, auto, none, uncut, plain, fraction] = pythonExportCheck([
    pdf(await loaded(10000)),
    // A CLI or Python API Session's leaf goes as stored: Python's validator reads it.
    pdf(await loaded('5kb')),
    pdf(await loaded('Auto')),
    pdf(await loaded('none')),
    pdf(await exportLinearWrap({ leaf: 'auto', summary: undefined, discoverRecords: discovered([8000]) })),
    pdf(await exportLinearWrap({ leaf: null, summary: undefined, discoverRecords: discovered([5_000_000]) })),
    { wrapLength: linearWrapInputValue('1.0000000001', 'kbp') }
  ]);
  for (const refusal of [cut, unit, auto]) {
    assert.equal(refusal.code, 'EXPORT_LIMIT', JSON.stringify(refusal));
    assert.equal(refusal.reason, 'PDF_SIZE');
  }
  assert.deepEqual(none, { ok: true }, '"none" is off');
  assert.deepEqual(uncut, { ok: true }, 'nothing is cut: an unwrapped PDF stays as it was');
  assert.deepEqual(plain, { ok: true });
  // D-57: a typed amount that is no whole bp reaches Python as typed and is refused with its diagnostic.
  assert.deepEqual(fraction, { code: 'INPUT_INVALID', configPath: 'canvas.linear.wrap_length' });
});

test('a loaded Result with a cropped record sends its region and Python checks the bases the drawing drew (D-70)', async () => {
  // The D-59 residual: a loaded Result whose card crops its record. The card's source length and region go as they
  // are; Python crops them with the drawing's rule, so 2 Mb of a 5 Mb record at 10 kbp rows is refused, and 8 kbp of
  // it, which no 10 kbp row cuts, is an unwrapped PDF, as the CLI checks them. The Web computes no cropped length.
  const cropped = (start, end) => async () => [{ sourceLength: 5_000_000, region: { start, end } }];
  const sent = await exportLinearWrap({ leaf: 10000, summary: undefined, discoverRecords: cropped(1_000_001, 3_000_000) });
  assert.deepEqual(sent, { requested: 10000, records: [{ sourceLength: 5_000_000, region: { start: 1_000_001, end: 3_000_000 } }] });
  const pdf = (linearWrap) => ({ format: 'pdf', width: 1500, height: 34_000, dpi: null, linearWrap });
  const [cut, short] = pythonExportCheck([
    pdf(sent),
    pdf(await exportLinearWrap({ leaf: 10000, summary: undefined, discoverRecords: cropped(1_001, 9_000) }))
  ]);
  assert.equal(cut.code, 'EXPORT_LIMIT', JSON.stringify(cut));
  assert.equal(cut.reason, 'PDF_SIZE');
  assert.deepEqual(short, { ok: true }, 'the crop cuts no row: an unwrapped PDF stays as it was');
});
