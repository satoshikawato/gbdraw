// @ts-check
// The Linear wrap length of a drawing's form (`form.linear_wrap_length`): what
// the form holds and how a stored value reads into it. The request leaf is
// `canvas.linear.wrap_length` (services/session-request.js compiles it);
// Python decides the row length and whether a record is cut
// (gbdraw/layout/linear_wrap.py), and parses a length typed with a unit.
/** @import { LinearWrapSummary } from './input-summary.js' */

/**
 * `null`: off; `'auto'`: Python's automatic row length; a number: a row length in bp, a positive safe integer.
 * @typedef {null | 'auto' | number} LinearWrapFormValue
 */
/**
 * A stored value read into the form: `lost` when the form cannot hold it (a length typed with a unit, such as
 * `"50kb"`, which only Python reads), and then the value is off.
 * @typedef {{ value: LinearWrapFormValue, lost: boolean }} LinearWrapReading
 */

/**
 * Reads a stored wrap length (a Session's form value, a request or config leaf, a preserved override): absent,
 * `null` and Python's `"none"` are off, `"auto"` is Auto, a positive safe integer is a length in bp.
 * @param {unknown} raw
 * @returns {LinearWrapReading}
 */
export const wrapFormValue = (raw) => {
  if (raw === undefined || raw === null || raw === 'none') return { value: null, lost: false };
  if (raw === 'auto') return { value: 'auto', lost: false };
  if (typeof raw === 'number' && Number.isSafeInteger(raw) && raw > 0) return { value: raw, lost: false };
  return { value: null, lost: true };
};

// Turning Wrap Rows on switches the settings rows cannot draw with, and the one rows decide (receipt PD-OI-113:
// each row draws its ruler on its axis), in the same edit: Features on axis to Above, Bar (Simple) to Ruler
// (Ticks), Normalize Record Lengths off, Ruler on Axis on. It does so also when no record is cut, so the figure is
// the form's visible values. Turning it off writes only the length; Undo restores the others.
/**
 * @typedef {{ linear_wrap_length: LinearWrapFormValue | number, linear_track_layout?: 'above' | 'below',
 *   scale_style?: 'ruler', normalize_length?: false, linear_ruler_on_axis?: true }} LinearWrapFormPatch
 *   A length that is no positive safe integer (a typed `1.5`) is sent as typed, and Python refuses it.
 */
/** The form fields a Wrap Rows edit writes; the layout transition saves and restores them. */
export const LINEAR_WRAP_FORM_FIELDS = Object.freeze(
  /** @type {const} */ (['linear_wrap_length', 'linear_track_layout', 'scale_style', 'normalize_length', 'linear_ruler_on_axis'])
);

/**
 * The fields one Wrap Rows edit writes.
 * @param {{ linear_track_layout: unknown }} form The drawing's form before the edit.
 * @param {LinearWrapFormValue | number} value
 * @returns {LinearWrapFormPatch}
 */
export const linearWrapFormPatch = (form, value) => (value === null
  ? { linear_wrap_length: null }
  : {
      linear_wrap_length: value,
      linear_track_layout: form.linear_track_layout === 'below' ? 'below' : 'above',
      scale_style: 'ruler',
      normalize_length: false,
      linear_ruler_on_axis: true
    });

/**
 * What the shown Linear drawing holds that rows cannot draw with, each read from the owner its request reads.
 * @typedef {object} LinearWrapBlockerInputs
 * @property {boolean} comparisons BLAST, comparison pairs, LOSAT or Similarity groups.
 * @property {boolean} sharedRecordRow Record Layout puts several records in one row.
 * @property {boolean} featureAlignment Records aligned at a feature (a similarity alignment).
 * @property {boolean} featureTrackOnAxis Custom Track Slots put the feature track on the axis (side overlay).
 */
/**
 * Python's refusal reason for each (gbdraw/layout/linear_wrap.py), so the field's hint and the Generate refusal
 * show one text (`reasonText`).
 * @typedef {'WRAP' | 'WRAP_GRID_ROW' | 'WRAP_ANCHOR_FEATURE' | 'WRAP_TRACK_LAYOUT'} LinearWrapBlocker
 */
/** @type {ReadonlyArray<[keyof LinearWrapBlockerInputs, LinearWrapBlocker]>} */
const BLOCKERS = Object.freeze([
  ['comparisons', 'WRAP'], ['sharedRecordRow', 'WRAP_GRID_ROW'],
  ['featureAlignment', 'WRAP_ANCHOR_FEATURE'], ['featureTrackOnAxis', 'WRAP_TRACK_LAYOUT']
]);
/**
 * What keeps Wrap Rows from being turned on or its length changed. Turning it off is never blocked, so a Session
 * that holds both can follow the refusal's "turn Wrap Rows off". Whether a record is cut is Python's.
 * @param {LinearWrapBlockerInputs} inputs
 * @returns {LinearWrapBlocker[]}
 */
export const linearWrapBlockers = (inputs) => BLOCKERS.filter(([input]) => inputs[input]).map(([, blocker]) => blocker);

// The length's unit is view state (D-27): the form and the Session hold bp.
export const LINEAR_WRAP_UNITS = Object.freeze({ bp: 1, kbp: 1_000, Mbp: 1_000_000 });
/** @typedef {keyof typeof LINEAR_WRAP_UNITS} LinearWrapUnit */
/** @type {readonly LinearWrapUnit[]} */
const LARGEST_FIRST = Object.freeze(['Mbp', 'kbp', 'bp']);

/** @param {unknown} unit @returns {unit is LinearWrapUnit} */
export const isLinearWrapUnit = (unit) => typeof unit === 'string' && Object.hasOwn(LINEAR_WRAP_UNITS, unit);

/**
 * The unit a length shows in when it comes from outside the field (Load, Undo, Reset): the largest one whose
 * amount is a whole number; Auto and Off keep `current`.
 * @param {unknown} value
 * @param {LinearWrapUnit} current
 * @returns {LinearWrapUnit}
 */
export const linearWrapUnitOf = (value, current) => (typeof value === 'number' && Number.isSafeInteger(value) && value > 0
  ? /** @type {LinearWrapUnit} */ (LARGEST_FIRST.find((unit) => value % LINEAR_WRAP_UNITS[unit] === 0))
  : current);

/**
 * The amount the length field shows: empty for Auto (D-28) and Off.
 * @param {unknown} value
 * @param {LinearWrapUnit} unit
 */
export const linearWrapAmountText = (value, unit) => (typeof value === 'number' ? String(value / LINEAR_WRAP_UNITS[unit]) : '');

/**
 * The value a typed amount writes: empty is Auto (D-28); otherwise amount x unit in bp, a whole number when the
 * decimal amount names one (`1.1` kbp is 1100 bp). Any other number is sent as typed and Python refuses it.
 * @param {string} text
 * @param {LinearWrapUnit} unit
 * @returns {'auto' | number}
 */
export const linearWrapInputValue = (text, unit) => {
  const trimmed = String(text).trim();
  const amount = Number(trimmed);
  if (!trimmed || !Number.isFinite(amount)) return 'auto';
  const product = amount * LINEAR_WRAP_UNITS[unit];
  const whole = Math.round(product);
  return Math.abs(product - whole) < 1e-6 ? whole : product;
};

/**
 * What the shown Result says about the drawing's Wrap Rows value. `off`: wrap is off. `rows`: the Result cut
 * records with this value over these records (its input summary's `linearWrap`). `one-row`: the Result asked for
 * this value and cut nothing, so the axis is the longest record. `unknown`: no Result drew this draft.
 * @typedef {{ kind: 'off' } | { kind: 'unknown' } | { kind: 'one-row' }
 *   | { kind: 'rows', rowLength: number, rowCount: number, drawnAxisLength: number }} ShownLinearWrap
 */
/**
 * @typedef {object} ShownLinearWrapInputs
 * @property {unknown} value The drawing's `form.linear_wrap_length`.
 * @property {LinearWrapFormValue | undefined} requested The shown Result's committed request's wrap length (read
 *   with `wrapFormValue`); undefined when no Result of this drawing is shown.
 * @property {LinearWrapSummary | null | undefined} summary The shown Result's input summary's `linearWrap`;
 *   undefined when the render sent no summary for it (a loaded Result).
 * @property {readonly number[]} recordLengths The lengths of the records the next Generate draws, [] while unknown.
 */
/**
 * Reads Python's summary only; the Auto rule, the cut test and the drawn axis stay Python's.
 * @param {ShownLinearWrapInputs} inputs
 * @returns {ShownLinearWrap}
 */
export const shownLinearWrap = ({ value, requested, summary, recordLengths }) => {
  if (value === null) return { kind: 'off' };
  if (requested === undefined || requested !== value || summary === undefined) return { kind: 'unknown' };
  if (summary === null) return { kind: 'one-row' };
  const sameRecords = summary.recordLengths.length === recordLengths.length
    && summary.recordLengths.every((length, index) => length === recordLengths[index]);
  return summary.requested === value && sameRecords
    ? { kind: 'rows', rowLength: summary.rowLength, rowCount: summary.rowCount, drawnAxisLength: summary.drawnAxisLength }
    : { kind: 'unknown' };
};
