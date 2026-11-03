// @ts-check
// The input summary Python sends beside each Result (`metadata.inputSummaries`,
// gbdraw/api/input_summary.py): what the Result read and drew, and the warnings
// for a record or a comparison table it draws nothing of. The CLI logs the same
// lines and warnings; Python owns every count and word, so the Web shows them
// as sent. Not saved in a Session.
//
// A live Feature visibility edit that no rerender draws (a hide) changes what
// the Result draws without Python, so the counts no longer describe it. Each
// entry keeps the Feature visibility of the render that sent it
// (`drawnVisibility`, `visibilityPaintState`); while the drawing's differs, the
// page says the summary describes the last render and keeps its warnings back.
import { diagnosticError } from '../utils/error-normalization.js';

/**
 * @typedef {object} InputWarning
 * @property {string} code A Python warning code, such as `no_drawable_features`.
 * @property {string} message
 */
/**
 * The rows a wrapped Linear Result drew (gbdraw/api/input_summary.py `LinearWrapPayload`).
 * @typedef {object} LinearWrapSummary
 * @property {number} rowLength The bases per row.
 * @property {number} rowCount
 * @property {'auto' | number} requested What the request asked for: `'auto'` or a length in bp.
 * @property {readonly number[]} recordRows The rows of each drawn record, in drawing order (the summary's
 *   `recordLengths` are their bases).
 * @property {number} drawnAxisLength The bases on the longest drawn axis, which sizes and ticks follow.
 */
/**
 * @typedef {object} InputSummary
 * @property {number} resultIndex
 * @property {string} resultName
 * @property {readonly string[]} lines
 * @property {readonly InputWarning[]} warnings
 * @property {readonly number[]} recordLengths The bases of each drawn record, in drawing order.
 * @property {LinearWrapSummary | null} linearWrap `null` when the Result cut no record.
 */
/**
 * @typedef {InputSummary & { readonly drawnVisibility: string }} DrawnInputSummary
 *   A Result's summary with the Feature visibility of the render that sent it (Web-owned).
 */

const SUMMARY_FIELDS = ['resultIndex', 'resultName', 'lines', 'warnings', 'recordLengths', 'linearWrap'];
const WARNING_FIELDS = ['code', 'message'];
const WRAP_FIELDS = ['rowLength', 'rowCount', 'requested', 'recordRows', 'drawnAxisLength'];

/** @param {unknown} value @param {readonly string[]} fields @returns {value is Record<string, unknown>} */
const hasExactly = (value, fields) => (
  value !== null && typeof value === 'object' && !Array.isArray(value)
  && Object.keys(value).length === fields.length
  && fields.every((field) => Object.prototype.hasOwnProperty.call(value, field))
);

/** @param {unknown} value @returns {value is string} */
const isText = (value) => typeof value === 'string' && value.length > 0;

/** @param {unknown} value @returns {value is InputWarning} */
const isWarning = (value) => hasExactly(value, WARNING_FIELDS)
  && isText(value.code) && /^[a-z_]+$/.test(value.code) && isText(value.message);

/** @param {unknown} value @returns {value is number} */
const isCount = (value) => Number.isSafeInteger(value) && Number(value) > 0;

/** @param {unknown} value @returns {value is readonly number[]} */
const isCounts = (value) => Array.isArray(value) && value.length > 0 && value.every(isCount);

/**
 * @param {unknown} value
 * @param {readonly number[]} recordLengths The summary's records: a wrapped Result names the rows of each.
 * @returns {value is LinearWrapSummary | null}
 */
const isLinearWrap = (value, recordLengths) => value === null || (hasExactly(value, WRAP_FIELDS)
  && isCount(value.rowLength) && isCount(value.rowCount) && isCount(value.drawnAxisLength)
  && (value.requested === 'auto' || isCount(value.requested))
  && isCounts(value.recordRows) && value.recordRows.length === recordLengths.length);

/**
 * The input summaries of a successful render, one per Result, each naming its Result by index and name.
 * @param {unknown} value
 * @param {{ readonly [index: number]: { readonly name?: unknown } | null | undefined }} results
 * @returns {readonly InputSummary[]}
 */
export const validateInputSummaries = (value, results) => {
  if (value === undefined) return Object.freeze([]);
  if (!Array.isArray(value) || value.some((entry) => (
    !hasExactly(entry, SUMMARY_FIELDS)
    || !Number.isSafeInteger(entry.resultIndex) || Number(entry.resultIndex) < 0
    || typeof entry.resultName !== 'string'
    || results[Number(entry.resultIndex)]?.name !== entry.resultName
    || !Array.isArray(entry.lines) || !entry.lines.every(isText)
    || !Array.isArray(entry.warnings) || !entry.warnings.every(isWarning)
    || !isCounts(entry.recordLengths) || !isLinearWrap(entry.linearWrap, entry.recordLengths)
  ))) throw diagnosticError('RESULT_INVALID', { field: 'inputSummaries' }, { stage: 'result-admission', operation: 'validateInputSummaries' });
  return Object.freeze(/** @type {InputSummary[]} */ (value));
};

/**
 * The summaries of a render, each with the Feature visibility the render drew.
 * @param {readonly InputSummary[]} summaries
 * @param {string} drawnVisibility `visibilityPaintState` of the drawing the render read.
 * @returns {readonly DrawnInputSummary[]}
 */
export const drawnInputSummaries = (summaries, drawnVisibility) => Object.freeze(
  summaries.map((entry) => Object.freeze({ ...entry, drawnVisibility }))
);

/**
 * The summary of the shown Result, or null when the render sent none for it.
 * @param {readonly DrawnInputSummary[]} summaries
 * @param {number} resultIndex
 * @param {string | undefined} resultName
 * @returns {DrawnInputSummary | null}
 */
export const inputSummaryOf = (summaries, resultIndex, resultName) => (
  summaries.find((entry) => entry.resultIndex === resultIndex && entry.resultName === resultName) ?? null
);
