// @ts-check
// The Wrap Rows field of the Linear Layout card (PD-OI-113). The field writes
// the drawing's `form.linear_wrap_length` only through the layout transition
// (`changeLayoutSetting`), which also switches the settings rows decide
// (`linearWrapFormPatch`); the length's unit is view state (D-27).
/** @import { DrawingState } from '../state.js' */
/** @import { LinearWrapBlocker, LinearWrapBlockerInputs, LinearWrapUnit } from '../services/linear-wrap-form.js' */
import {
  isLinearWrapUnit, linearWrapAmountText, linearWrapBlockers, linearWrapInputValue, linearWrapUnitOf
} from '../services/linear-wrap-form.js';
import { reasonText } from '../utils/error-normalization.js';
import { joinSentences, tm } from '../utils/i18n.js';
import { computed, ref } from '../utils/reactivity.js';

/**
 * @typedef {object} LinearWrapFieldPorts
 * @property {{ activeDrawing: () => DrawingState }} state
 * @property {(event: Event, field: string, value?: unknown) => unknown} changeLayoutSetting The layout transition.
 * @property {() => LinearWrapBlockerInputs} blockerInputs What the shown Linear drawing holds that rows cannot draw with.
 * @property {(source: () => unknown, callback: () => void) => unknown} watch Vue `watch`: reads a value
 *   from outside the field (Load, Undo, Reset, a drawing switch).
 */

/** Which part of the field an edit comes from. @typedef {'on' | 'amount' | 'unit'} LinearWrapFieldPart */

/** @param {LinearWrapFieldPorts} ports */
export const createLinearWrapField = ({ state, changeLayoutSetting, blockerInputs, watch }) => {
  const form = () => state.activeDrawing().form;
  const linearWrapOn = computed(() => form().linear_wrap_length !== null);
  /** @type {{ value: LinearWrapUnit }} */
  const linearWrapUnit = ref(/** @type {LinearWrapUnit} */ ('kbp'));
  // The value the field wrote last: a value from anywhere else (Load, Undo,
  // Reset, another drawing) shows in the largest unit that holds it whole.
  /** @type {unknown} */
  let written;
  watch(() => [state.activeDrawing().id, form().linear_wrap_length], () => {
    const value = form().linear_wrap_length;
    if (value !== written) linearWrapUnit.value = linearWrapUnitOf(value, linearWrapUnit.value);
    written = undefined;
  });
  const linearWrapAmount = computed(() => linearWrapAmountText(form().linear_wrap_length, linearWrapUnit.value));
  /** @type {{ readonly value: LinearWrapBlocker[] }} */
  const blockers = computed(() => linearWrapBlockers(blockerInputs()));
  // D-29: a blocker keeps Wrap Rows from being turned on or its length changed, never from being turned off.
  const linearWrapBlocked = computed(() => blockers.value.length > 0);
  // One visible hint names every blocker, in the words of its Generate refusal.
  const linearWrapBlockedText = computed(() => joinSentences(...blockers.value.map((blocker) => tm(/** @type {string} */ (reasonText(blocker))))));
  /**
   * @param {Event} event A change of the checkbox, the length, or the unit.
   * @param {LinearWrapFieldPart} part
   */
  const changeLinearWrap = (event, part) => {
    const control = /** @type {HTMLInputElement | HTMLSelectElement} */ (event.target);
    // A unit change keeps the amount shown: 20 kbp becomes 20 Mbp.
    const amount = part === 'amount' ? control.value : linearWrapAmount.value;
    if (part === 'unit' && isLinearWrapUnit(control.value)) linearWrapUnit.value = control.value;
    const value = part === 'on'
      ? (/** @type {HTMLInputElement} */ (control).checked ? 'auto' : null)
      : linearWrapInputValue(amount, linearWrapUnit.value);
    written = value;
    return changeLayoutSetting(event, 'linear_wrap_length', value);
  };
  return { linearWrapOn, linearWrapUnit, linearWrapAmount, linearWrapBlocked, linearWrapBlockedText, changeLinearWrap };
};
