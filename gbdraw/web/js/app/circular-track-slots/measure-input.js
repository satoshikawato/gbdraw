// @ts-check
import {
  readCircularMeasure,
  writeCircularMeasureValue,
  changeCircularMeasureUnit
} from '../../services/circular-track-measure.js';
import { circularMeasureInvalidParts } from '../../services/track-slot-validation.js';
import { errorText } from '../../utils/error-text.js';
import { msg, tm } from '../../utils/i18n.js';

/** @import { CircularMeasureField } from '../../services/track-slot-validation.js' */

/**
 * The texts of one field, as catalog keys: its label and the accessible names
 * of its value and unit controls (`{id}` is the slot id).
 * @typedef {object} CircularMeasureFieldTexts
 * @property {string} label
 * @property {string} valueName
 * @property {string} unitName
 */
/** @type {Readonly<Record<CircularMeasureField, CircularMeasureFieldTexts>>} */
const FIELD_TEXTS = Object.freeze({
  width: Object.freeze({
    label: msg('Width'),
    valueName: msg('Circular track slot {id} Width value'),
    unitName: msg('Circular track slot {id} Width unit')
  }),
  radius: Object.freeze({
    label: msg('Radius'),
    valueName: msg('Circular track slot {id} Radius value'),
    unitName: msg('Circular track slot {id} Radius unit')
  })
});

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
  /**
   * @param {{ modelValue: unknown, slotId: string, field: string, controlId: string, autoText: string, disabled: boolean }} props
   * @param {{ emit: (event: 'update:modelValue', value: unknown) => void }} context
   */
  setup(props, { emit }) {
    const { computed, ref } = window.Vue;
    // Only Auto has a local next-input preference. Manual units come from the slot.
    const autoUnit = ref('factor');
    /** @returns {CircularMeasureField} */
    const field = () => (props.field === 'width' ? 'width' : 'radius');
    const view = computed(() => readCircularMeasure(props.modelValue));
    const selectedUnit = computed({
      get: () => view.value.isAuto ? autoUnit.value : view.value.selectedUnit,
      set: (/** @type {string} */ unit) => {
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
      set: (/** @type {string} */ text) => emit('update:modelValue', writeCircularMeasureValue(text, selectedUnit.value))
    });
    const label = computed(() => tm(FIELD_TEXTS[field()].label));
    const valueName = computed(() => tm(FIELD_TEXTS[field()].valueName, { id: props.slotId }));
    const unitName = computed(() => tm(FIELD_TEXTS[field()].unitName, { id: props.slotId }));
    // The validator's sentence for this field (TK-12); empty while the value is valid.
    const error = computed(() => (view.value.invalid ? errorText(circularMeasureInvalidParts(field())) : ''));
    const describedBy = computed(() => [
      'circular-measure-help',
      view.value.isAuto && props.autoText ? `${props.controlId}-auto` : '',
      view.value.invalid ? `${props.controlId}-error` : ''
    ].filter(Boolean).join(' '));
    return { view, valueText, selectedUnit, label, valueName, unitName, error, describedBy };
  }
};
