// @ts-check
// `v-committed-value="x"`: the value a text, number or color field shows,
// written in place of `:value="x"` on every `<input>` and `<textarea>` the
// user types into. Vue patches a bound `value` on every render of the
// component that owns the field and puts the model's value back whenever the
// field's text differs, so text the user had typed but not committed, or text
// a commit refused, was replaced by any render of the page (a title, Generate
// progress, the end of the boot, a search). This directive is the one owner of
// what such a field shows:
// - A render writes the model's value into the field, unless an IME
//   composition is open, or the model's value is the one the field showed
//   before and the field has focus or a commit of it is still running (the
//   rule of Vue's own `v-model.lazy`, extended to a running commit).
// - When a commit ends (`change`, `blur`, or an Enter the field's own handler
//   took), the field shows the model's value once the commit has settled, with
//   or without focus: the typed text was either stored or refused. A handler whose
//   commit ends after an `await` hands its promise over with
//   `settleCommittedField`; a synchronous one needs nothing.
// tests/web/committed-value.test.mjs keeps `:value` off typed fields.

import { listHoldsFocus } from './list-focus.js';

/** @typedef {HTMLInputElement | HTMLTextAreaElement} CommittedField */
/** @typedef {{ value: unknown }} CommittedValueBinding */

/**
 * What the directive holds for one field.
 * @typedef {object} CommittedFieldState
 * @property {string} model The model's value as text, at the field's last render.
 * @property {boolean} composing Whether an IME composition is open.
 * @property {number} running The commits handed over that have not settled.
 * @property {Promise<void>} settled Resolves once every commit handed over so far has settled.
 * @property {number} turn The number of the latest commit end; an earlier one shows nothing.
 * @property {() => void} stop Removes the listeners.
 */

/** @type {WeakMap<EventTarget, CommittedFieldState>} */
const fieldStates = new WeakMap();

/**
 * The state of a field the directive has mounted. A field without it is not
 * bound with `v-committed-value`, a template or lifecycle error.
 * @param {EventTarget} field
 */
const stateOf = (field) => {
  const state = fieldStates.get(field);
  if (!state) throw new Error('This field is not bound with v-committed-value.');
  return state;
};

/**
 * The text a field shows for a model value: empty for none.
 * @param {unknown} value
 */
export const committedText = (value) => (value === null || value === undefined ? '' : String(value));

/**
 * Whether a render leaves the field's own text: while an IME composition is
 * open, and while the model's value has not changed and the field has focus or
 * a commit of it is running.
 * @param {{ focused: boolean, composing: boolean, committing: boolean, changed: boolean }} field
 */
export const keepsTypedText = ({ focused, composing, committing, changed }) => composing || ((focused || committing) && !changed);

/** @param {CommittedField} field @param {string} text */
const show = (field, text) => {
  if (field.value !== text) field.value = text;
};

// A field a keyed list moves has the focus again after the move (app/list-focus.js).
/** @param {CommittedField} field */
const hasFocus = (field) => field.ownerDocument.activeElement === field || listHoldsFocus(field);

/**
 * Hands the result of a commit of `field` to the directive, which shows the
 * model's value once it has settled; returns `commit`. A template calls it
 * around an asynchronous handler that is not given the field
 * (`settleCommittedField($event.target, rename(...))`), and a handler given the
 * field calls it around its own commit. Call it in the event that commits.
 * @template T
 * @param {EventTarget | null} field
 * @param {T} commit
 * @returns {T}
 */
export const settleCommittedField = (field, commit) => {
  if (!field) throw new Error('A commit was handed over without its field.');
  const state = stateOf(field);
  const done = Promise.resolve(commit).then(() => {}, () => {});
  state.running += 1;
  state.settled = Promise.all([state.settled, done]).then(() => {});
  void done.then(() => { state.running -= 1; });
  return commit;
};

/**
 * After a commit ends: Vue has run its handler; wait for the commits handed
 * over and for the render they cause, then show the model's value. A commit
 * end that came later shows instead, and text typed meanwhile into the
 * focused field stays.
 * @param {CommittedField} field
 * @param {CommittedFieldState} state
 */
const settle = async (field, state) => {
  state.turn += 1;
  const turn = state.turn;
  const text = field.value;
  await window.Vue.nextTick();
  /** @type {Promise<void>} */
  let settled;
  do {
    settled = state.settled;
    await settled;
    await window.Vue.nextTick();
  } while (settled !== state.settled);
  if (turn !== state.turn || !field.isConnected || state.composing) return;
  if (hasFocus(field) && field.value !== text) return;
  show(field, state.model);
};

export const committedValue = {
  /** @param {CommittedField} field @param {CommittedValueBinding} binding */
  beforeMount(field, { value }) {
    /** @type {CommittedFieldState} */
    const state = { model: committedText(value), composing: false, running: 0, settled: Promise.resolve(), turn: 0, stop: () => {} };
    // Vue adds the template's listeners before this hook, so a handler has run
    // (and handed its commit over) before these do.
    const onCommitEnd = () => { void settle(field, state); };
    // An Enter the template's handler took (`@keydown.enter.prevent`) ends its
    // commit here. Any other Enter is the browser's: it commits by the
    // `change` of the keypress that follows, and showing the model's value at
    // keydown would undo the typed text before that `change` reads it.
    /** @param {Event} event */
    const onKeydown = (event) => {
      const key = /** @type {KeyboardEvent} */ (event);
      if (key.key === 'Enter' && key.defaultPrevented && !key.isComposing) onCommitEnd();
    };
    const onCompositionStart = () => { state.composing = true; };
    const onCompositionEnd = () => { state.composing = false; };
    /** @type {[string, (event: Event) => void][]} */
    const listeners = [
      ['change', onCommitEnd], ['blur', onCommitEnd], ['keydown', onKeydown],
      ['compositionstart', onCompositionStart], ['compositionend', onCompositionEnd]
    ];
    for (const [type, listener] of listeners) field.addEventListener(type, listener);
    state.stop = () => {
      for (const [type, listener] of listeners) field.removeEventListener(type, listener);
    };
    fieldStates.set(field, state);
    show(field, state.model);
  },
  /** @param {CommittedField} field @param {CommittedValueBinding} binding */
  updated(field, { value }) {
    const state = stateOf(field);
    const model = committedText(value);
    const changed = model !== state.model;
    state.model = model;
    if (keepsTypedText({ focused: hasFocus(field), composing: state.composing, committing: state.running > 0, changed })) return;
    show(field, model);
  },
  /** @param {CommittedField} field */
  unmounted(field) {
    stateOf(field).stop();
    fieldStates.delete(field);
  }
};
