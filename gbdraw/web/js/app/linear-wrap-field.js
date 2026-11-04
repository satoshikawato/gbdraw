// @ts-check
// The Wrap Rows field of the Linear Layout card (PD-OI-113). The field writes
// the drawing's `form.linear_wrap_length` only through the layout transition
// (`changeLayoutSetting`), which also switches the settings rows decide
// (`linearWrapFormPatch`); the length's unit is view state (D-27).
/** @import { DrawingState } from '../state.js' */
/** @import { LinearWrapBlocker, LinearWrapBlockerInputs, LinearWrapUnit, ShownLinearWrap } from '../services/linear-wrap-form.js' */
import {
  isLinearWrapUnit, linearWrapAmountText, linearWrapBlockers, linearWrapInputValue, linearWrapUnitOf, offersRecordPositionReset
} from '../services/linear-wrap-form.js';
import { reasonText } from '../utils/error-normalization.js';
import { joinSentences, t, tm } from '../utils/i18n.js';
import { formatBp } from './auto-value-display.js';
import { computed, ref } from '../utils/reactivity.js';

/**
 * @typedef {object} LinearWrapFieldPorts
 * @property {{ activeDrawing: () => DrawingState }} state
 * @property {(event: Event, field: string, value?: unknown) => unknown} changeLayoutSetting The layout transition.
 * @property {() => LinearWrapBlockerInputs} blockerInputs What the shown Linear drawing holds that rows cannot draw with.
 * @property {() => ShownLinearWrap} shownWrap What the shown Result drew with the drawing's value.
 * @property {(source: () => unknown, callback: () => void) => unknown} watch Vue `watch`: reads a value
 *   from outside the field (Load, Undo, Reset, a drawing switch).
 */

/** Which part of the field an edit comes from. @typedef {'on' | 'amount' | 'unit'} LinearWrapFieldPart */

/** @param {LinearWrapFieldPorts} ports */
export const createLinearWrapField = ({ state, changeLayoutSetting, blockerInputs, shownWrap, watch }) => {
  const form = () => state.activeDrawing().form;
  const linearWrapOn = computed(() => form().linear_wrap_length !== null);
  // D-56: Normalize Record Lengths and Ruler on Axis are fixed only while they hold what rows need (off, on), as
  // Track Layout and Scale Style disable only the conflicting option; a loaded other value stays editable toward it.
  const linearWrapFixesNormalize = computed(() => linearWrapOn.value && !form().normalize_length);
  const linearWrapFixesRulerOnAxis = computed(() => linearWrapOn.value && Boolean(form().linear_ruler_on_axis));
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
  // Beside an empty length (Auto): what the shown Result drew with Auto, else just Auto.
  const linearWrapAutoText = computed(() => {
    const reading = shownWrap();
    if (form().linear_wrap_length !== 'auto') return t('Auto');
    if (reading.kind === 'rows') {
      return t('{length} per row, {rows} rows (auto)', { length: formatBp(reading.rowLength), rows: reading.rowCount });
    }
    return reading.kind === 'one-row' ? t('One row (auto)') : t('Auto');
  });
  /** @type {{ readonly value: LinearWrapBlocker[] }} */
  const blockers = computed(() => linearWrapBlockers(blockerInputs()));
  // D-29: a blocker keeps Wrap Rows from being turned on or its length changed, never from being turned off.
  const linearWrapBlocked = computed(() => blockers.value.length > 0);
  // One visible hint names every blocker, in the words of its Generate refusal.
  const linearWrapBlockedText = computed(() => joinSentences(...blockers.value.map((blocker) => tm(/** @type {string} */ (reasonText(blocker))))));
  // A drawn record moved sideways, or a Similarity alignment, offers Reset record positions under the hint, which names
  // it (PD-OI-113 scenario revision 2, Q5 A).
  const linearWrapResetsPositions = computed(() => offersRecordPositionReset(blockerInputs()));
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
  return {
    linearWrapOn, linearWrapFixesNormalize, linearWrapFixesRulerOnAxis, linearWrapUnit, linearWrapAmount, linearWrapAutoText,
    linearWrapBlocked, linearWrapBlockedText, linearWrapResetsPositions, changeLinearWrap
  };
};
