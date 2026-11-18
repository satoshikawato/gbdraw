import assert from 'node:assert/strict';
import test from 'node:test';
import { CircularMeasureInput } from '../../gbdraw/web/js/app/circular-track-slots/measure-input.js';
import { setUiLocale } from '../../gbdraw/web/js/utils/i18n.js';
import { loadPageVue } from './helpers/vue-runtime.mjs';

const Vue = loadPageVue();
globalThis.window = { Vue };
const { reactive } = Vue;

const setup = (scalar = null, field = 'width') => {
  const props = reactive({ modelValue: scalar, slotId: 'gc_content', field, controlId: `slot-${field}`, autoText: '18 px (auto)' });
  const emissions = [];
  const controls = CircularMeasureInput.setup(props, { emit: (event, value) => {
    emissions.push({ event, value });
    props.modelValue = value;
  } });
  return { props, controls, emissions };
};

test('component reads manual fields and History restores without emitting or mirroring units', () => {
  const { props, controls, emissions } = setup({ value: '1.', unit: 'px' });
  assert.equal(controls.valueText.value, '1.');
  assert.equal(controls.selectedUnit.value, 'px');
  assert.equal(controls.valueName.value, 'Circular track slot gc_content Width value');
  assert.equal(controls.unitName.value, 'Circular track slot gc_content Width unit');
  assert.equal(controls.describedBy.value, 'circular-measure-help');
  controls.selectedUnit.value = 'px';
  assert.equal(emissions.length, 0);
  controls.selectedUnit.value = 'factor';
  assert.equal(emissions.length, 1);
  assert.deepEqual({ ...props.modelValue }, { value: '1.', unit: 'factor' });
  props.modelValue = { value: 20, unit: 'px' };
  assert.equal(controls.valueText.value, '20');
  assert.equal(controls.selectedUnit.value, 'px');
  assert.equal(emissions.length, 1);
});

test('Auto preference never emits, preserves saved blank, and resets with component lifetime', () => {
  const { props, controls, emissions } = setup('');
  assert.equal(controls.selectedUnit.value, 'px');
  controls.selectedUnit.value = 'factor';
  assert.equal(props.modelValue, '');
  assert.equal(emissions.length, 0);
  assert.equal(controls.describedBy.value, 'circular-measure-help slot-width-auto');
  controls.valueText.value = ' 1.5 ';
  assert.deepEqual({ ...props.modelValue }, { value: '1.5', unit: 'factor' });
  assert.equal(emissions.length, 1);
  controls.valueText.value = '';
  assert.equal(props.modelValue, null);
  assert.equal(setup(props.modelValue).controls.selectedUnit.value, 'px');
});

// OV-2314 (D-03): an Auto Width's hint is in px (the layout resolves a width
// in px) and an Auto Radius's in ×R. While a field is Auto its unit selector
// starts at the unit of its hint, from the same owner, so a number typed over
// the hint means what the hint says; a bare number in the slot stays a factor.
test('an Auto field starts in the unit of its Auto hint, and a number typed over it means that unit', async () => {
  const { CIRCULAR_MEASURE_AUTO_UNITS, formatCircularMeasureAuto } = await import('../../gbdraw/web/js/services/track-slot-display.js');
  assert.deepEqual({ ...CIRCULAR_MEASURE_AUTO_UNITS }, { width: 'px', radius: 'factor' });
  assert.equal(formatCircularMeasureAuto('width', 18), '18 px (auto)');
  assert.equal(formatCircularMeasureAuto('radius', 0.65), '0.65 R (auto)');
  const width = setup(null, 'width');
  assert.equal(width.controls.selectedUnit.value, 'px');
  width.controls.valueText.value = '18';
  assert.deepEqual({ ...width.props.modelValue }, { value: '18', unit: 'px' });
  const radius = setup(null, 'radius');
  assert.equal(radius.controls.selectedUnit.value, 'factor');
  radius.controls.valueText.value = '0.65';
  assert.deepEqual({ ...radius.props.modelValue }, { value: '0.65', unit: 'factor' });
  // A stored bare number is a factor in both fields, as the CLI reads it.
  assert.equal(setup('1.5', 'width').controls.selectedUnit.value, 'factor');
});

test('numeric edits use the shared suffix adapter and retain invalid drafts and units', () => {
  const { props, controls, emissions } = setup({ value: 20, unit: 'px' });
  controls.valueText.value = '65%';
  assert.deepEqual({ ...props.modelValue }, { value: '0.65', unit: 'factor' });
  controls.valueText.value = '20px';
  assert.deepEqual({ ...props.modelValue }, { value: '20', unit: 'px' });
  controls.valueText.value = '1epx';
  assert.equal(controls.valueText.value, '1epx');
  assert.equal(controls.error.value, 'Width must be a number greater than 0, in px or ×R.');
  assert.equal(controls.describedBy.value, 'circular-measure-help slot-width-error');
  controls.selectedUnit.value = 'factor';
  assert.deepEqual({ ...props.modelValue }, { value: '1epx', unit: 'factor' });
  assert.equal(emissions.length, 4);
});

// OV-1408: the label, the accessible names, and the invalid-value sentence
// (the validator's Width/Radius sentence) show in the UI language; English is
// unchanged.
test('the field shows its label, names, and invalid value in the UI language (OV-1408)', async () => {
  const shown = (field, scalar) => {
    const { controls } = setup(scalar, field);
    return [controls.label.value, controls.valueName.value, controls.unitName.value, controls.error.value];
  };
  assert.deepEqual(shown('width', '0'), ['Width', 'Circular track slot gc_content Width value',
    'Circular track slot gc_content Width unit', 'Width must be a number greater than 0, in px or ×R.']);
  assert.deepEqual(shown('radius', { value: 'bad', unit: 'px' }), ['Radius', 'Circular track slot gc_content Radius value',
    'Circular track slot gc_content Radius unit', 'Radius must be a number greater than 0, in px or ×R.']);
  assert.equal(shown('radius', '0.5')[3], '');
  await setUiLocale('ja');
  try {
    assert.deepEqual(shown('width', '0'), ['幅', '環状図のトラックスロットgc_content：幅の値',
      '環状図のトラックスロットgc_content：幅の単位', '幅は0より大きい数値を、pxまたは×Rで指定してください。']);
    assert.deepEqual(shown('radius', '-1'), ['半径', '環状図のトラックスロットgc_content：半径の値',
      '環状図のトラックスロットgc_content：半径の単位', '半径は0より大きい数値を、pxまたは×Rで指定してください。']);
  } finally {
    await setUiLocale('en');
  }
});
