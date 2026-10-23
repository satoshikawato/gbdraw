// @ts-check
// The label fit report Python sends (gbdraw/web_support/label_fit.py): the one
// JavaScript reader of it. Python owns every value and word of the report:
// the counts and each record's counts text, the rows' titles, the settings a
// row changes by config path, the canvas sizes, and the reasons (LP-D43,
// LP-D48). Nothing here places labels or searches values; a row that was not
// computed is computed by the Worker (`computeLabelFitRow`, OD-F14).
import { diagnosticError } from '../utils/error-normalization.js';

/**
 * @typedef {'found' | 'found_with_notice' | 'not_found' | 'not_computed'} LabelFitRowStatus
 * `found_with_notice`: the row leaves no overlapping labels but leaves crossing
 * leaders or labels over the definition text, so the figure draws with the
 * notice bar (LP-D57); its `note` says what it leaves.
 */
/** @typedef {string & { readonly labelFitRowKind: unique symbol }} LabelFitRowKind A Python row kind, such as `outward_vertical`. */
/** @typedef {string & { readonly configPath: unique symbol }} ConfigPath A dotted config path, such as `labels.spacing.circular`. */
/** @typedef {string & { readonly labelFitRowToken: unique symbol }} LabelFitRowToken Names the drawn figure whose rows the Worker computes. */
/** @typedef {number | string | boolean} LabelFitValue A number, a word such as a label scope, or a switch. */
/**
 * @typedef {object} LabelFitSetting One control a row changes, as Python names it.
 * @property {string} name
 * @property {string} unit
 * @property {number | string} value A number, or the control's option label.
 */
/**
 * @typedef {object} LabelFitRow
 * @property {LabelFitRowKind} kind
 * @property {string} title
 * @property {LabelFitRowStatus} status
 * @property {Readonly<Record<ConfigPath, LabelFitValue>> | null} values The settings to write, by config path; set when offered.
 * @property {readonly LabelFitSetting[]} settings
 * @property {{ width: number, height: number } | null} canvasSize
 * @property {string | null} reason
 * @property {string | null} note What a `found_with_notice` row leaves, in Python's words.
 */
/**
 * @typedef {object} LabelFitRecord
 * @property {string} recordId
 * @property {string} labelKind
 * @property {string} labelKindTitle
 * @property {number} labelCount
 * @property {number} fontSize
 * @property {number} overlapPairs
 * @property {number} crossingPairs
 * @property {number} definitionOverlaps
 * @property {number} regionOverlaps Region labels over other labels or lines (Linear).
 * @property {boolean} refused Radial labels that could not be placed (counts may be 0).
 * @property {string} counts Python's counts, or why the labels were refused.
 */
/**
 * @typedef {object} LabelFitReport
 * @property {boolean} drawn True when the figure was drawn (crossing leaders only).
 * @property {number} overlapPairs
 * @property {number} crossingPairs
 * @property {number} definitionOverlaps Labels over the definition text.
 * @property {number} regionOverlaps Region labels over other labels or lines (Linear).
 * @property {readonly string[]} notice Python's sentences for the notice bar over a drawn figure.
 * @property {readonly LabelFitRecord[]} records
 * @property {readonly LabelFitRow[]} rows In the dialog order.
 */
/**
 * @typedef {object} ResultLabelFit The report of one drawn Result.
 * @property {number} resultIndex
 * @property {string} resultName
 * @property {LabelFitReport} report
 * @property {LabelFitRowToken | null} rowToken
 */
/**
 * @typedef {object} LabelFitFailure The report of a Generate that drew nothing.
 * @property {LabelFitReport} report
 * @property {LabelFitRowToken | null} rowToken Names the failure whose rows the Worker computes.
 */

const ROW_STATUSES = new Set(['found', 'found_with_notice', 'not_found', 'not_computed']);
const OFFERED_STATUSES = new Set(['found', 'found_with_notice']);

/**
 * Whether a row's status offers its settings to Apply.
 * @param {string} status
 * @returns {boolean}
 */
export const labelFitRowOffered = (status) => OFFERED_STATUSES.has(status);

const invalid = () => diagnosticError('RESULT_INVALID', {}, { stage: 'result-admission' });

/** @param {unknown} value @returns {value is Record<string, unknown>} */
const isObject = (value) => value !== null && typeof value === 'object' && !Array.isArray(value);
/** @param {unknown} value @returns {value is number} */
const isCount = (value) => Number.isSafeInteger(value) && Number(value) >= 0;
/** @param {unknown} value @returns {value is number} */
const isNumber = (value) => typeof value === 'number' && Number.isFinite(value);
/** @param {unknown} value @returns {value is number | string} */
const isShown = (value) => isNumber(value) || (typeof value === 'string' && value !== '');
/** @param {unknown} value @returns {value is LabelFitValue} */
const isValue = (value) => isShown(value) || typeof value === 'boolean';
/** @param {unknown} value @returns {value is LabelFitRowToken | null} */
const isToken = (value) => value === null || (typeof value === 'string' && value !== '');
/** @param {unknown} value @returns {value is string | null} */
const isOptionalText = (value) => value === null || value === undefined || typeof value === 'string';

/**
 * One row, as the report or the Worker's `computeLabelFitRow` gives it.
 * @param {unknown} value
 * @returns {LabelFitRow}
 */
export const normalizeLabelFitRow = (value) => {
  if (!isObject(value) || typeof value.kind !== 'string' || !value.kind
    || typeof value.title !== 'string' || typeof value.status !== 'string'
    || !ROW_STATUSES.has(value.status) || !Array.isArray(value.settings)
    || !isOptionalText(value.reason) || !isOptionalText(value.note)) throw invalid();
  const status = /** @type {LabelFitRowStatus} */ (value.status);
  const offered = labelFitRowOffered(status);
  /** @type {Record<ConfigPath, LabelFitValue> | null} */
  let values = null;
  if (offered) {
    if (!isObject(value.values) || !Object.values(value.values).every(isValue)) throw invalid();
    values = /** @type {Record<ConfigPath, LabelFitValue>} */ ({ ...value.values });
  }
  const size = value.canvasSize;
  if (offered && !(isObject(size) && isNumber(size.width) && isNumber(size.height))) throw invalid();
  const settings = value.settings.map((setting) => {
    if (!isObject(setting) || typeof setting.name !== 'string' || typeof setting.unit !== 'string'
      || !isShown(setting.value)) throw invalid();
    return Object.freeze({ name: setting.name, unit: setting.unit, value: setting.value });
  });
  return Object.freeze({
    kind: /** @type {LabelFitRowKind} */ (value.kind),
    title: value.title,
    status,
    values: values && Object.freeze(values),
    settings: Object.freeze(settings),
    canvasSize: offered && isObject(size)
      ? Object.freeze({ width: Number(size.width), height: Number(size.height) }) : null,
    reason: typeof value.reason === 'string' ? value.reason : null,
    note: typeof value.note === 'string' ? value.note : null
  });
};

/**
 * A report in a Result's metadata or a failure's.
 * @param {unknown} parsed
 * @returns {LabelFitReport}
 */
const normalizeLabelFitReport = (parsed) => {
  if (!isObject(parsed) || parsed.schema !== 1 || typeof parsed.drawn !== 'boolean'
    || !isCount(parsed.overlapPairs) || !isCount(parsed.crossingPairs) || !isCount(parsed.definitionOverlaps)
    || !isCount(parsed.regionOverlaps)
    || !Array.isArray(parsed.notice) || !parsed.notice.every((sentence) => typeof sentence === 'string')
    || !Array.isArray(parsed.records) || !Array.isArray(parsed.rows)) throw invalid();
  const records = parsed.records.map((record) => {
    if (!isObject(record) || typeof record.recordId !== 'string' || typeof record.labelKind !== 'string'
      || typeof record.labelKindTitle !== 'string' || !isCount(record.labelCount)
      || !isNumber(record.fontSize) || !isCount(record.overlapPairs) || !isCount(record.crossingPairs)
      || !isCount(record.definitionOverlaps) || !isCount(record.regionOverlaps) || typeof record.refused !== 'boolean' || typeof record.counts !== 'string') {
      throw invalid();
    }
    return Object.freeze({
      recordId: record.recordId,
      labelKind: record.labelKind,
      labelKindTitle: record.labelKindTitle,
      labelCount: record.labelCount,
      fontSize: record.fontSize,
      overlapPairs: record.overlapPairs,
      crossingPairs: record.crossingPairs,
      definitionOverlaps: record.definitionOverlaps,
      regionOverlaps: record.regionOverlaps,
      refused: record.refused,
      counts: record.counts
    });
  });
  return Object.freeze({
    drawn: parsed.drawn,
    overlapPairs: parsed.overlapPairs,
    crossingPairs: parsed.crossingPairs,
    definitionOverlaps: parsed.definitionOverlaps,
    regionOverlaps: parsed.regionOverlaps,
    notice: Object.freeze([...parsed.notice]),
    records: Object.freeze(records),
    rows: Object.freeze(parsed.rows.map(normalizeLabelFitRow))
  });
};

/**
 * The reports of a successful render's Results (`metadata.labelFitReports`).
 * Each names a Result of `results` by index and name.
 * @param {unknown} value
 * @param {readonly { name?: unknown }[]} results
 * @returns {readonly ResultLabelFit[]}
 */
export const normalizeResultLabelFits = (value, results) => {
  if (value === undefined) return Object.freeze([]);
  if (!Array.isArray(value)) throw invalid();
  return Object.freeze(value.map((entry) => {
    if (!isObject(entry) || !isCount(entry.resultIndex) || typeof entry.resultName !== 'string'
      || results[entry.resultIndex]?.name !== entry.resultName) {
      throw invalid();
    }
    if (!isToken(entry.rowToken)) throw invalid();
    const report = normalizeLabelFitReport(entry.report);
    return Object.freeze({
      resultIndex: entry.resultIndex,
      resultName: entry.resultName,
      report,
      rowToken: /** @type {LabelFitRowToken | null} */ (entry.rowToken)
    });
  }));
};

/**
 * The JSON text next to a Generate error whose labels do not fit: the report
 * and the token of its rows to compute. Absent is null.
 * @param {unknown} value
 * @returns {LabelFitFailure | null}
 */
export const normalizeLabelFitFailure = (value) => {
  if (value === null || value === undefined) return null;
  let parsed;
  try { parsed = JSON.parse(String(value)); } catch (_) { throw invalid(); }
  if (typeof value !== 'string' || !isObject(parsed) || !isToken(parsed.rowToken)) throw invalid();
  return Object.freeze({
    report: normalizeLabelFitReport(parsed.report),
    rowToken: /** @type {LabelFitRowToken | null} */ (parsed.rowToken)
  });
};
