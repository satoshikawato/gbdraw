// @ts-check
// `v-committed-value="x"`: the value a text, number or color field shows when
// it commits on change, blur or Enter, written in place of `:value="x"`.
// Vue patches a bound `value` on every render of the component that owns the
// field and puts the model's value back whenever the field's text differs, so
// text the user typed but had not committed was replaced by any render of the
// page (a title, Generate progress, the end of the boot, a search). This is
// the rule of Vue's own `v-model.lazy` instead: a render writes the model's
// value into the field unless the field has focus and the value is the one it
// already showed, and never during an IME composition. After the field loses
// focus it shows the model's value, so a value its commit refused or
// normalized does not stay on screen. tests/web/committed-value.test.mjs keeps
// `:value` off such fields in the templates.

/** @typedef {HTMLInputElement | HTMLTextAreaElement} CommittedField */
/** @typedef {{ value: unknown }} CommittedValueBinding */

/**
 * The model's value of one field as its text, and whether an IME composition
 * is open in it.
 * @typedef {object} CommittedFieldState
 * @property {string} model
 * @property {boolean} composing
 * @property {() => void} stop Removes the listeners.
 */

/** @type {WeakMap<CommittedField, CommittedFieldState>} */
const fieldStates = new WeakMap();

/**
 * The text a field shows for a model value: empty for none.
 * @param {unknown} value
 */
export const committedText = (value) => (value === null || value === undefined ? '' : String(value));

/**
 * Whether a render leaves the field's own text: while an IME composition is
 * open, and while the field has focus and the model's value has not changed.
 * @param {{ focused: boolean, composing: boolean, changed: boolean }} field
 */
export const keepsTypedText = ({ focused, composing, changed }) => composing || (focused && !changed);

/** @param {CommittedField} field @param {string} text */
const show = (field, text) => {
  if (field.value !== text) field.value = text;
};

/** @param {CommittedField} field */
const hasFocus = (field) => field.ownerDocument.activeElement === field;

export const committedValue = {
  /** @param {CommittedField} field @param {CommittedValueBinding} binding */
  beforeMount(field, { value }) {
    /** @type {CommittedFieldState} */
    const state = { model: committedText(value), composing: false, stop: () => {} };
    // The change event, and the commit its handler makes, come before blur; the
    // field shows the value that commit left once Vue has rendered it.
    const onBlur = () => {
      void window.Vue.nextTick(() => {
        if (field.isConnected && !hasFocus(field)) show(field, state.model);
      });
    };
    const onCompositionStart = () => { state.composing = true; };
    const onCompositionEnd = () => { state.composing = false; };
    field.addEventListener('blur', onBlur);
    field.addEventListener('compositionstart', onCompositionStart);
    field.addEventListener('compositionend', onCompositionEnd);
    state.stop = () => {
      field.removeEventListener('blur', onBlur);
      field.removeEventListener('compositionstart', onCompositionStart);
      field.removeEventListener('compositionend', onCompositionEnd);
    };
    fieldStates.set(field, state);
    show(field, state.model);
  },
  /** @param {CommittedField} field @param {CommittedValueBinding} binding */
  updated(field, { value }) {
    const state = fieldStates.get(field);
    if (!state) return;
    const model = committedText(value);
    const changed = model !== state.model;
    state.model = model;
    if (keepsTypedText({ focused: hasFocus(field), composing: state.composing, changed })) return;
    show(field, model);
  },
  /** @param {CommittedField} field */
  unmounted(field) {
    fieldStates.get(field)?.stop();
    fieldStates.delete(field);
  }
};
