// @ts-check
// The label fit report Python sends (gbdraw/web_support/label_fit.py): the one
// JavaScript reader of it. Python owns every value and word of the report:
// the counts, the rows' titles, the settings a row changes by config path,
// the canvas sizes, and the reasons (LP-D43, LP-D48). Nothing here places
// labels or searches values; a row that was not computed is computed by the
// Worker (`computeLabelFitRow`, OD-F14).
import { diagnosticError } from '../utils/error-normalization.js';

/** @typedef {'found' | 'not_found' | 'not_computed'} LabelFitRowStatus */
/** @typedef {string & { readonly labelFitRowKind: unique symbol }} LabelFitRowKind A Python row kind, such as `outward_vertical`. */
/** @typedef {string & { readonly configPath: unique symbol }} ConfigPath A dotted config path, such as `labels.spacing.circular`. */
/** @typedef {string & { readonly labelFitRowToken: unique symbol }} LabelFitRowToken Names the drawn figure whose rows the Worker computes. */
/**
 * @typedef {object} LabelFitSetting One control a row changes, as Python names it.
 * @property {string} name
 * @property {string} unit
 * @property {number} value
 */
/**
 * @typedef {object} LabelFitRow
 * @property {LabelFitRowKind} kind
 * @property {string} title
 * @property {LabelFitRowStatus} status
 * @property {Readonly<Record<ConfigPath, number>> | null} values The settings to write, by config path; set when `found`.
 * @property {readonly LabelFitSetting[]} settings
 * @property {{ width: number, height: number } | null} canvasSize
 * @property {string | null} reason
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
 */
/**
 * @typedef {object} LabelFitReport
 * @property {boolean} drawn True when the figure was drawn (crossing leaders only).
 * @property {number} overlapPairs
 * @property {number} crossingPairs
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

const ROW_STATUSES = new Set(['found', 'not_found', 'not_computed']);

const invalid = () => diagnosticError('RESULT_INVALID', {}, { stage: 'result-admission' });

/** @param {unknown} value @returns {value is Record<string, unknown>} */
const isObject = (value) => value !== null && typeof value === 'object' && !Array.isArray(value);
/** @param {unknown} value @returns {value is number} */
const isCount = (value) => Number.isSafeInteger(value) && Number(value) >= 0;
/** @param {unknown} value @returns {value is number} */
const isNumber = (value) => typeof value === 'number' && Number.isFinite(value);
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
    || !isOptionalText(value.reason)) throw invalid();
  const status = /** @type {LabelFitRowStatus} */ (value.status);
  /** @type {Record<ConfigPath, number> | null} */
  let values = null;
  if (status === 'found') {
    if (!isObject(value.values) || !Object.values(value.values).every(isNumber)) throw invalid();
    values = /** @type {Record<ConfigPath, number>} */ ({ ...value.values });
  }
  const size = value.canvasSize;
  if (status === 'found' && !(isObject(size) && isNumber(size.width) && isNumber(size.height))) throw invalid();
  const settings = value.settings.map((setting) => {
    if (!isObject(setting) || typeof setting.name !== 'string' || typeof setting.unit !== 'string'
      || !isNumber(setting.value)) throw invalid();
    return Object.freeze({ name: setting.name, unit: setting.unit, value: setting.value });
  });
  return Object.freeze({
    kind: /** @type {LabelFitRowKind} */ (value.kind),
    title: value.title,
    status,
    values: values && Object.freeze(values),
    settings: Object.freeze(settings),
    canvasSize: status === 'found' && isObject(size)
      ? Object.freeze({ width: Number(size.width), height: Number(size.height) }) : null,
    reason: typeof value.reason === 'string' ? value.reason : null
  });
};

/**
 * The report next to a Generate error, or in a Result's metadata. A JSON text
 * (the error path) is parsed first; an absent report is null.
 * @param {unknown} value
 * @returns {LabelFitReport | null}
 */
export const normalizeLabelFitReport = (value) => {
  if (value === null || value === undefined) return null;
  let parsed = value;
  if (typeof value === 'string') {
    try { parsed = JSON.parse(value); } catch (_) { throw invalid(); }
  }
  if (!isObject(parsed) || parsed.schema !== 1 || typeof parsed.drawn !== 'boolean'
    || !isCount(parsed.overlapPairs) || !isCount(parsed.crossingPairs)
    || !Array.isArray(parsed.records) || !Array.isArray(parsed.rows)) throw invalid();
  const records = parsed.records.map((record) => {
    if (!isObject(record) || typeof record.recordId !== 'string' || typeof record.labelKind !== 'string'
      || typeof record.labelKindTitle !== 'string' || !isCount(record.labelCount)
      || !isNumber(record.fontSize) || !isCount(record.overlapPairs) || !isCount(record.crossingPairs)) {
      throw invalid();
    }
    return Object.freeze({
      recordId: record.recordId,
      labelKind: record.labelKind,
      labelKindTitle: record.labelKindTitle,
      labelCount: record.labelCount,
      fontSize: record.fontSize,
      overlapPairs: record.overlapPairs,
      crossingPairs: record.crossingPairs
    });
  });
  return Object.freeze({
    drawn: parsed.drawn,
    overlapPairs: parsed.overlapPairs,
    crossingPairs: parsed.crossingPairs,
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
      || results[entry.resultIndex]?.name !== entry.resultName
      || !(entry.rowToken === null || (typeof entry.rowToken === 'string' && entry.rowToken))) {
      throw invalid();
    }
    const report = normalizeLabelFitReport(entry.report);
    if (!report) throw invalid();
    return Object.freeze({
      resultIndex: entry.resultIndex,
      resultName: entry.resultName,
      report,
      rowToken: /** @type {LabelFitRowToken | null} */ (entry.rowToken)
    });
  }));
};
