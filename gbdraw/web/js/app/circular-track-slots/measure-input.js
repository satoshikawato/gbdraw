// @ts-check
import {
  readCircularMeasure,
  writeCircularMeasureValue,
  changeCircularMeasureUnit
} from './measure-editor.js';

export const CircularMeasureInput = {
  template: '#circular-measure-input-template',
  props: {
    modelValue: { default: null },
    slotId: { type: String, required: true },
    field: { type: String, required: true },
    controlId: { type: String, required: true },
    autoText: { type: String, default: '' },
    disabled: { type: Boolean, default: false }
  },
  emits: ['update:modelValue'],
  setup(props, { emit }) {
    const { computed, ref } = window.Vue;
    // Only Auto has a local next-input preference. Manual units come from the slot.
    const autoUnit = ref('factor');
    const view = computed(() => readCircularMeasure(props.modelValue));
    const selectedUnit = computed({
      get: () => view.value.isAuto ? autoUnit.value : view.value.selectedUnit,
      set: (unit) => {
        if (view.value.isAuto) {
          autoUnit.value = unit;
          return;
        }
        const scalar = changeCircularMeasureUnit(props.modelValue, unit);
        if (scalar !== props.modelValue) emit('update:modelValue', scalar);
      }
    });
    const valueText = computed({
      get: () => view.value.valueText,
      set: (text) => emit('update:modelValue', writeCircularMeasureValue(text, selectedUnit.value))
    });
    const label = computed(() => props.field === 'width' ? 'Width' : 'Radius');
    const accessibleName = computed(() => `Circular track slot ${props.slotId} ${label.value}`);
    const describedBy = computed(() => [
      'circular-measure-help',
      view.value.isAuto && props.autoText ? `${props.controlId}-auto` : '',
      view.value.error ? `${props.controlId}-error` : ''
    ].filter(Boolean).join(' '));
    return { view, valueText, selectedUnit, label, accessibleName, describedBy };
  }
};
