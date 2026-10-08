// @ts-check
import { normalizeTrackIndex } from '../../services/circular-track-slot-model.js';

// The Circular Depth row's track index (TK-07). A committed value that is not
// a non-negative integer stays visible with a field error, and the row keeps
// its last valid index; the field never replaces the entry silently.
export const DepthTrackIndexInput = {
  template: '#depth-track-index-input-template',
  props: {
    modelValue: { default: null },
    slotId: { type: String, required: true },
    controlId: { type: String, required: true },
    helpId: { type: String, default: '' },
    disabled: { type: Boolean, default: false }
  },
  emits: ['update:modelValue'],
  /**
   * @param {{ modelValue: unknown, slotId: string, controlId: string, helpId: string, disabled: boolean }} props
   * @param {{ emit: (event: 'update:modelValue', value: number) => void }} context
   */
  setup(props, { emit }) {
    const { computed, ref, watch } = window.Vue;
    /** @type {{ value: string | null }} The rejected entry; null shows the row's index. */
    const rejected = ref(null);
    watch(() => props.modelValue, () => { rejected.value = null; });
    const text = computed(() => rejected.value ?? String(normalizeTrackIndex(props.modelValue) ?? ''));
    const error = computed(() => (
      rejected.value === null
        ? ''
        : `Enter a whole number of 0 or more. The row keeps track index ${normalizeTrackIndex(props.modelValue) ?? 0}.`
    ));
    const describedBy = computed(() => [props.helpId, error.value ? `${props.controlId}-error` : '']
      .filter(Boolean).join(' '));
    /** @param {string} value */
    const commit = (value) => {
      const index = normalizeTrackIndex(String(value).trim());
      if (index === null) {
        rejected.value = String(value);
        return;
      }
      rejected.value = null;
      if (index !== normalizeTrackIndex(props.modelValue)) emit('update:modelValue', index);
    };
    return { text, error, describedBy, commit };
  }
};
