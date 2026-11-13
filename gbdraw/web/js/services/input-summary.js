// @ts-check
// The input summary Python sends beside each Result (`metadata.inputSummaries`,
// gbdraw/api/input_summary.py): what the Result read and drew, and the warnings
// for a record or a comparison table it draws nothing of, or for a comparison
// file whose record IDs repeat (drawn in full, named ambiguously). The CLI logs the same
// lines and warnings; Python owns every count and word, so the Web shows them
// as sent. Not saved in a Session.
//
// A live Feature visibility edit that no rerender draws (a hide) changes what
// the Result draws without Python, so the counts no longer describe it. Each
// entry keeps the Feature visibility of the render that sent it
// (`drawnVisibility`, `visibilityPaintState`); while the drawing's differs, the
// page says the summary describes the last render and keeps its warnings back.
import { INPUT_WARNINGS_OF_DRAWN_INPUT } from '../contracts/schema.generated.js';
import { diagnosticError } from '../utils/error-normalization.js';
import { limitedList } from '../utils/i18n.js';

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
 *   `records`).
 * @property {number} drawnAxisLength The bases on the longest drawn axis, which sizes and ticks follow.
 */
/**
 * One drawn record (gbdraw/api/input_summary.py `RecordPayload`). Its drawn bases are Python's to derive from these
 * two facts with the drawing's crop rule (gbdraw/web_support/export_limit.py, D-70); the Web predicts none.
 * @typedef {object} RecordSummary
 * @property {number} sourceLength The bases of its source record.
 * @property {{ start: number, end: number } | null} region The region of the source the request asked for (1-based,
 *   both ends included, as asked: the end may run past the record); null when it asked for none.
 */
/**
 * @typedef {object} InputSummary
 * @property {number} resultIndex
 * @property {string} resultName
 * @property {readonly string[]} lines
 * @property {readonly InputWarning[]} warnings
 * @property {readonly RecordSummary[]} records Each drawn record, in drawing order.
 * @property {LinearWrapSummary | null} linearWrap `null` when the Result cut no record.
 */
/**
 * The warnings the preview notice shows in one group: `left-out` says part of the input was not drawn, `drawn-in`
 * warns of an input drawn in full. `shown` and `rest` are the group cut to `NOTICE_LINE_LIMIT` lines.
 * @typedef {{ readonly kind: 'left-out' | 'drawn-in' } & import('../utils/i18n.js').LimitedList<InputWarning>} InputWarningGroup
 */
/**
 * @typedef {InputSummary & {
 *   readonly drawnVisibility: string, readonly warningGroups: readonly InputWarningGroup[]
 * }} DrawnInputSummary
 *   A Result's summary with the Feature visibility of the render that sent it (Web-owned), and its warnings in the
 *   groups the notice shows, each with warnings, left-out first.
 */

/**
 * The lines a warning group of the preview notice takes at most, its "and {count} more" line included (OV-3495:
 * 59 empty record pairs of one LOSATP run pushed the figure out of view). The rest opens on request.
 */
export const NOTICE_LINE_LIMIT = 4;

const SUMMARY_FIELDS = ['resultIndex', 'resultName', 'lines', 'warnings', 'records', 'linearWrap'];
const WARNING_FIELDS = ['code', 'message'];
const RECORD_FIELDS = ['sourceLength', 'region'];
const REGION_FIELDS = ['start', 'end'];
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

/** @param {unknown} value @returns {value is RecordSummary} */
const isRecord = (value) => hasExactly(value, RECORD_FIELDS) && isCount(value.sourceLength)
  && (value.region === null || (hasExactly(value.region, REGION_FIELDS) && isCount(value.region.start) && isCount(value.region.end)));

/** @param {unknown} value @returns {value is readonly RecordSummary[]} */
const isRecords = (value) => Array.isArray(value) && value.length > 0 && value.every(isRecord);

/**
 * @param {unknown} value
 * @param {readonly RecordSummary[]} records The summary's records: a wrapped Result names the rows of each.
 * @returns {value is LinearWrapSummary | null}
 */
const isLinearWrap = (value, records) => value === null || (hasExactly(value, WRAP_FIELDS)
  && isCount(value.rowLength) && isCount(value.rowCount) && isCount(value.drawnAxisLength)
  && (value.requested === 'auto' || isCount(value.requested))
  && isCounts(value.recordRows) && value.recordRows.length === records.length);

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
    || !isRecords(entry.records) || !isLinearWrap(entry.linearWrap, entry.records)
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
  summaries.map((entry) => Object.freeze({
    ...entry,
    drawnVisibility,
    warningGroups: Object.freeze(/** @type {const} */ (['left-out', 'drawn-in']).flatMap((kind) => {
      // A warning of an input drawn in full (Python's `INPUT_WARNINGS_OF_DRAWN_INPUT`) leaves nothing out.
      const warnings = entry.warnings.filter((warning) => INPUT_WARNINGS_OF_DRAWN_INPUT.includes(warning.code) === (kind === 'drawn-in'));
      return warnings.length ? [Object.freeze({ kind, ...limitedList(warnings, NOTICE_LINE_LIMIT) })] : [];
    }))
  }))
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
