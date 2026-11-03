// @ts-check
// The Linear wrap length of a drawing's form (`form.linear_wrap_length`): what
// the form holds and how a stored value reads into it. The request leaf is
// `canvas.linear.wrap_length` (services/session-request.js compiles it);
// Python decides the row length and whether a record is cut
// (gbdraw/layout/linear_wrap.py), and parses a length typed with a unit.

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
