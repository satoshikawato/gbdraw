// The Web's Wrap Rows: the words of its refusals and of the export limit (one wording owner,
// utils/error-normalization.js), the form edit that switches what rows decide (services/linear-wrap-form.js
// through the layout transition), its blockers, and the length field's amount and unit. The value's request
// and Session channel is tested in session-request, mode-scoped-session and session50-load-resave.
import assert from 'node:assert/strict';
import test from 'node:test';

import { normalizeUserFacingError, reasonText } from '../../gbdraw/web/js/utils/error-normalization.js';
import { createDefaultAdv, createDefaultForm } from '../../gbdraw/web/js/services/session-active-config-contract.js';
import {
  LINEAR_WRAP_FORM_FIELDS, linearWrapAmountText, linearWrapBlockers, linearWrapFormPatch, linearWrapInputValue, linearWrapUnitOf
} from '../../gbdraw/web/js/services/linear-wrap-form.js';
import {
  createFeaturePlacementActions, draftFeatureSlot, draftPlacementTargets
} from '../../gbdraw/web/js/app/feature-editor/placement-actions.js';
import { createLinearWrapField } from '../../gbdraw/web/js/app/linear-wrap-field.js';
import { withDrawings } from './helpers/drawing-state.mjs';

// What Python sends for each wrap refusal (gbdraw/layout/linear_wrap.py; tests/test_web_wrap_diagnostics.py).
const WRAP_REFUSALS = Object.freeze({
  WRAP: ['COMPARISON_INPUT', 'Remove the comparison inputs (BLAST, comparison pairs, LOSAT, Similarity groups), or turn Wrap Rows off.'],
  WRAP_TRACK_LAYOUT: ['INPUT_INVALID', 'Set Track Layout to Above or Below, move the feature track off the axis in Custom Track Slots, or turn Wrap Rows off.'],
  WRAP_SCALE_STYLE: ['INPUT_INVALID', 'Set Scale Style to Ruler (Ticks), or turn Wrap Rows off.'],
  WRAP_NORMALIZE_LENGTH: ['INPUT_INVALID', 'Turn Normalize Record Lengths off, or turn Wrap Rows off.'],
  WRAP_GRID_ROW: ['INPUT_INVALID', 'Give each record its own row, or turn Wrap Rows off.'],
  WRAP_ANCHOR_FEATURE: ['INPUT_INVALID', 'Reset the alignment, or turn Wrap Rows off.']
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
    + ' Lower the PNG DPI, or use SVG, which has no size limit.');
  assert.deepEqual(png.actions, ['use-svg']);
  const pdf = normalizeUserFacingError({ code: 'EXPORT_LIMIT', operation: 'export-pdf', stage: 'export-capture',
    context: { reason: 'PDF_SIZE', limit: 14400, width: 1500, height: 34475 } });
  assert.equal(pdf.summary, 'The figure is larger than this file format allows, so no file was made.'
    + ' PDF is limited to 14,400 pt per side, and this figure is 1,500 x 34,475 pt.'
    + ' Set a longer row length in Wrap Rows, or use SVG, which has no size limit.');
  const unbounded = normalizeUserFacingError({ code: 'EXPORT_LIMIT', operation: 'export-png', stage: 'export-capture',
    context: { reason: 'PNG_SIZE', limit: 32767, width: -1, height: 1.5, dpi: 2e7 } });
  assert.deepEqual(unbounded.context, { reason: 'PNG_SIZE', limit: 32767 }, 'only bounded counts are kept');
  assert.equal(unbounded.summary, 'The figure is larger than this file format allows, so no file was made.'
    + ' Lower the PNG DPI, or use SVG, which has no size limit.');
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

test('each blocker keeps Wrap Rows from turning on, in its Generate refusal\'s words; Simple-mode Middle is none (case 3)', () => {
  const none = { comparisons: false, sharedRecordRow: false, featureAlignment: false, featureTrackOnAxis: false };
  assert.deepEqual(linearWrapBlockers(none), []);
  assert.deepEqual(linearWrapBlockers({ ...none, comparisons: true }), ['WRAP']);
  assert.deepEqual(linearWrapBlockers({ ...none, sharedRecordRow: true }), ['WRAP_GRID_ROW']);
  assert.deepEqual(linearWrapBlockers({ ...none, featureAlignment: true }), ['WRAP_ANCHOR_FEATURE']);
  assert.deepEqual(linearWrapBlockers({ ...none, featureTrackOnAxis: true }), ['WRAP_TRACK_LAYOUT']);
  assert.deepEqual(linearWrapBlockers({ comparisons: true, sharedRecordRow: true, featureAlignment: true, featureTrackOnAxis: true }),
    ['WRAP', 'WRAP_GRID_ROW', 'WRAP_ANCHOR_FEATURE', 'WRAP_TRACK_LAYOUT']);
  // The feature slot the blocker reads is the one the placement targets read: in simple mode, Middle puts the
  // features on the axis, which turning wrap on switches (no blocker); a custom stack's overlay slot is one.
  const drawing = linearDrawing();
  assert.equal(draftFeatureSlot(drawing).side, 'overlay');

  const store = withDrawings(linearDrawing());
  let inputs = { ...none, comparisons: true };
  const writes = [];
  const field = createLinearWrapField({ state: store, blockerInputs: () => inputs, watch: () => {},
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
  assert.equal(linearWrapUnitOf(20000, 'bp'), 'kbp');
  assert.equal(linearWrapUnitOf(3_000_000, 'bp'), 'Mbp');
  assert.equal(linearWrapUnitOf(1500, 'Mbp'), 'bp');
  assert.equal(linearWrapUnitOf('auto', 'Mbp'), 'Mbp');
  assert.equal(linearWrapAmountText(20000, 'kbp'), '20');
  assert.equal(linearWrapAmountText('auto', 'kbp'), '');

  const drawing = linearDrawing();
  const store = withDrawings(drawing);
  const watchers = [];
  const field = createLinearWrapField({ state: store, watch: (_source, callback) => watchers.push(callback),
    blockerInputs: () => ({ comparisons: false, sharedRecordRow: false, featureAlignment: false, featureTrackOnAxis: false }),
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
