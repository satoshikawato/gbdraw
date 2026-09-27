import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import vm from 'node:vm';
import test from 'node:test';
import { CircularMeasureInput } from '../../gbdraw/web/js/app/circular-track-slots/measure-input.js';

const context = vm.createContext({ console });
vm.runInContext(readFileSync('gbdraw/web/vendor/vue/vue.global.js', 'utf8'), context);
globalThis.window = { Vue: context.Vue };
const { reactive } = context.Vue;

const setup = (scalar = null) => {
  const props = reactive({ modelValue: scalar, slotId: 'gc_content', field: 'width', controlId: 'slot-width', autoText: '18 px (auto)' });
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
  assert.equal(controls.accessibleName.value, 'Circular track slot gc_content Width');
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
  assert.equal(controls.selectedUnit.value, 'factor');
  controls.selectedUnit.value = 'px';
  assert.equal(props.modelValue, '');
  assert.equal(emissions.length, 0);
  assert.equal(controls.describedBy.value, 'circular-measure-help slot-width-auto');
  controls.valueText.value = ' 1.5 ';
  assert.deepEqual({ ...props.modelValue }, { value: '1.5', unit: 'px' });
  assert.equal(emissions.length, 1);
  controls.valueText.value = '';
  assert.equal(props.modelValue, null);
  assert.equal(setup(props.modelValue).controls.selectedUnit.value, 'factor');
});

test('numeric edits use the shared suffix adapter and retain invalid drafts and units', () => {
  const { props, controls, emissions } = setup({ value: 20, unit: 'px' });
  controls.valueText.value = '65%';
  assert.deepEqual({ ...props.modelValue }, { value: '0.65', unit: 'factor' });
  controls.valueText.value = '20px';
  assert.deepEqual({ ...props.modelValue }, { value: '20', unit: 'px' });
  controls.valueText.value = '1epx';
  assert.equal(controls.valueText.value, '1epx');
  assert.ok(controls.view.value.error);
  assert.equal(controls.describedBy.value, 'circular-measure-help slot-width-error');
  controls.selectedUnit.value = 'factor';
  assert.deepEqual({ ...props.modelValue }, { value: '1epx', unit: 'factor' });
  assert.equal(emissions.length, 4);
});
