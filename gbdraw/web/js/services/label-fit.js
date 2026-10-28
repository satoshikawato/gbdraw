// @ts-check
// The label fit report Python sends (gbdraw/web_support/label_fit.py): the one
// JavaScript reader of it. Python owns every value and word of the report:
// the counts and each record's counts text, the rows' titles, the settings a
// row changes by config path and by request track option, the canvas sizes,
// and the reasons (LP-D43, LP-D48). Nothing here places labels or searches values; a row that was not
// computed is computed by the Worker (`computeLabelFitRow`, OD-F14).
// The exceptions are the Web's own words for labels that need much more room
// than the map has: the ways to make room (CAPACITY-STOP), which Python names
// by id while the Web names its controls, and the room such labels need
// (CAPACITY-STOP, CAPACITY-WARN), which Python sends as values (OV-2179).
import { diagnosticError } from '../utils/error-normalization.js';
import { msg, tn } from '../utils/i18n.js';
/** @import { CanonicalCircularTrackOptions } from '../contracts/request.generated.js' */

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
 * @typedef {Partial<Record<Extract<keyof CanonicalCircularTrackOptions, 'centerReservedRadius'>, number>>} LabelFitTrackOptions
 * The canonical render request's Circular track options (`diagramOptions.tracks`) a row sets, by their request
 * names: "narrow the inside tracks" sets `centerReservedRadius` (PD-OI-097).
 */
/**
 * @typedef {object} LabelFitRow
 * @property {LabelFitRowKind} kind
 * @property {string} title
 * @property {LabelFitRowStatus} status
 * @property {Readonly<Record<ConfigPath, LabelFitValue>> | null} values The config overrides to write, by config path; set when offered.
 * @property {Readonly<LabelFitTrackOptions>} trackOptions The request track options to write, by request name; empty
 *   when the row sets none.
 * @property {readonly LabelFitSetting[]} settings
 * @property {{ width: number, height: number } | null} canvasSize
 * @property {string | null} reason
 * @property {string | null} note What a `found_with_notice` row leaves, in Python's words.
 */
/**
 * @typedef {object} LabelFitRecord
 * @property {string} recordId
 * @property {number | null} recordIndex The index of the record's drawn groups (`data-gbdraw-record-index`), its
 *   identity where two records share an ID; null for a Linear record.
 * @property {string} labelKind
 * @property {string} labelKindTitle
 * @property {number} labelCount
 * @property {number} fontSize
 * @property {number} overlapPairs
 * @property {number} crossingPairs
 * @property {number} definitionOverlaps
 * @property {number} regionOverlaps Region labels over other labels or lines (Linear).
 * @property {boolean} refused Radial labels that could not be placed (counts may be 0).
 * @property {LabelCapacity | null} capacity The room the labels need, when they were not placed (CAPACITY-STOP; from
 *   Python's `capacityRatio` and the record's kind, count and font size), worded by `labelCapacityText`; null
 *   otherwise.
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
 * @property {readonly LabelCapacityRemedy[]} remedies The ways to make room for the records not placed, in Python's
 *   order; empty unless a record has a `capacity`.
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

/**
 * The ways to make room for labels that were not placed (CAPACITY-STOP), by
 * Python's id (`LabelCapacityRemedy` in gbdraw/labels/fit.py, in its order),
 * each worded with the Web's controls; the CLI words the same ids with its
 * options. A report naming another id is invalid.
 */
export const LABEL_CAPACITY_REMEDY_WORDS = Object.freeze({
  filter: msg('Label fewer features: Blacklist or Whitelist under LABEL FILTERING.'),
  label_text: msg('Shorten the label text: a priority rule, for example gene before product.'),
  font: msg('Use a smaller Label Font Size.'),
  linear: msg('Draw a Linear map.')
});
/** @typedef {keyof typeof LABEL_CAPACITY_REMEDY_WORDS} LabelCapacityRemedy */

/**
 * The Web's words for the room each kind of label needs (CAPACITY-STOP, CAPACITY-WARN), by Python's label kind id:
 * its keys are the kinds whose room Python reads before placing them (`HorizontalLabelKind` in
 * gbdraw/labels/fit.py, in its order). `line` is the dialog's record line, whose `{record}` stays in the
 * text for the dialog's `<i18n-text>` slot (the record ID in bold); `notice` is the processing overlay's sentence. A
 * report or notice naming another kind is invalid. Each form is a whole sentence per kind and count (OV-2179).
 */
const LABEL_CAPACITY_KIND_WORDS = Object.freeze({
  circular_outer: {
    /** @param {number} count @param {Record<string, string | number>} params */
    line: (count, params) => tn(count,
      '{record}: {count} outer label at {size} pt needs about {ratio} times the room around the map',
      '{record}: {count} outer labels at {size} pt need about {ratio} times the room around the map', params),
    /** @param {number} count @param {Record<string, string | number>} params */
    notice: (count, params) => tn(count,
      '{record}: {count} outer label at {size} pt needs about {ratio} times the room around the map, so it may not fit. Placing it can take minutes.',
      '{record}: {count} outer labels at {size} pt need about {ratio} times the room around the map, so they may not fit. Placing them can take minutes.',
      params)
  },
  circular_inner: {
    /** @param {number} count @param {Record<string, string | number>} params */
    line: (count, params) => tn(count,
      '{record}: {count} inner label at {size} pt needs about {ratio} times the room around the map',
      '{record}: {count} inner labels at {size} pt need about {ratio} times the room around the map', params),
    /** @param {number} count @param {Record<string, string | number>} params */
    notice: (count, params) => tn(count,
      '{record}: {count} inner label at {size} pt needs about {ratio} times the room around the map, so it may not fit. Placing it can take minutes.',
      '{record}: {count} inner labels at {size} pt need about {ratio} times the room around the map, so they may not fit. Placing them can take minutes.',
      params)
  }
});
/** @typedef {keyof typeof LABEL_CAPACITY_KIND_WORDS} LabelCapacityKind */
/** The label kinds whose room Python reads, in its order (`HorizontalLabelKind`). */
export const LABEL_CAPACITY_KINDS = Object.freeze(/** @type {LabelCapacityKind[]} */ (Object.keys(LABEL_CAPACITY_KIND_WORDS)));
/**
 * @typedef {object} LabelCapacity The room one side of a record's labels needs.
 * @property {LabelCapacityKind} labelKind
 * @property {number} labelCount
 * @property {number} fontSize
 * @property {number} capacityRatio As a multiple of the room around the map, to one decimal place as Python rounds it.
 */
/**
 * @typedef {LabelCapacity & { recordId: string }} LabelCapacityNotice A record's labels that may not fit, sent while
 * the render runs (CAPACITY-WARN; `LabelCapacityNoticePayload` in gbdraw/web_support/label_fit.py).
 */
/** @param {unknown} value @returns {value is LabelCapacityKind} */
const isCapacityKind = (value) => typeof value === 'string' && Object.hasOwn(LABEL_CAPACITY_KIND_WORDS, value);

/**
 * The dialog's line for a record's labels that were not placed (CAPACITY-STOP), with `{record}` left for its slot.
 * @param {LabelCapacity} capacity
 * @returns {string}
 */
export const labelCapacityText = ({ labelKind, labelCount, fontSize, capacityRatio }) => (
  LABEL_CAPACITY_KIND_WORDS[labelKind].line(labelCount, { size: fontSize, ratio: capacityRatio })
);

/**
 * The sentence the processing overlay shows for labels that may not fit (CAPACITY-WARN).
 * @param {LabelCapacityNotice} notice
 * @returns {string}
 */
export const labelCapacityNoticeText = ({ recordId, labelKind, labelCount, fontSize, capacityRatio }) => (
  LABEL_CAPACITY_KIND_WORDS[labelKind].notice(labelCount, { record: recordId, size: fontSize, ratio: capacityRatio })
);

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
/** @param {unknown} value @returns {value is LabelCapacityRemedy} */
const isRemedy = (value) => typeof value === 'string' && Object.hasOwn(LABEL_CAPACITY_REMEDY_WORDS, value);
/** @param {unknown} value @returns {value is number} */
const isRatio = (value) => isNumber(value) && value > 0;
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
    || !isObject(value.trackOptions) || !Object.values(value.trackOptions).every(isNumber)
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
    trackOptions: Object.freeze(/** @type {LabelFitTrackOptions} */ ({ ...value.trackOptions })),
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
    || !Array.isArray(parsed.records) || !Array.isArray(parsed.rows)
    || !Array.isArray(parsed.remedies) || !parsed.remedies.every(isRemedy)) throw invalid();
  const records = parsed.records.map((record) => {
    if (!isObject(record) || typeof record.recordId !== 'string'
      || !(record.recordIndex === null || isCount(record.recordIndex)) || typeof record.labelKind !== 'string'
      || typeof record.labelKindTitle !== 'string' || !isCount(record.labelCount)
      || !isNumber(record.fontSize) || !isCount(record.overlapPairs) || !isCount(record.crossingPairs)
      || !isCount(record.definitionOverlaps) || !isCount(record.regionOverlaps) || typeof record.refused !== 'boolean'
      || typeof record.counts !== 'string') {
      throw invalid();
    }
    // Only the kinds whose room Python reads carry a ratio.
    const { labelKind, capacityRatio } = record;
    /** @type {LabelCapacity | null} */
    let capacity = null;
    if (capacityRatio !== null) {
      if (!isRatio(capacityRatio) || !isCapacityKind(labelKind)) throw invalid();
      capacity = Object.freeze({ labelKind, labelCount: record.labelCount, fontSize: record.fontSize, capacityRatio });
    }
    return Object.freeze({
      recordId: record.recordId,
      recordIndex: record.recordIndex,
      labelKind: record.labelKind,
      labelKindTitle: record.labelKindTitle,
      labelCount: record.labelCount,
      fontSize: record.fontSize,
      overlapPairs: record.overlapPairs,
      crossingPairs: record.crossingPairs,
      definitionOverlaps: record.definitionOverlaps,
      regionOverlaps: record.regionOverlaps,
      refused: record.refused,
      capacity,
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
    rows: Object.freeze(parsed.rows.map(normalizeLabelFitRow)),
    remedies: Object.freeze([...parsed.remedies])
  });
};

/**
 * A notice the Worker relays while the render runs: the JSON text of Python's
 * `LabelCapacityNoticePayload`. Text of another shape is no notice (null).
 * @param {unknown} value
 * @returns {LabelCapacityNotice | null}
 */
export const normalizeLabelCapacityNotice = (value) => {
  if (typeof value !== 'string') return null;
  let parsed;
  try { parsed = JSON.parse(value); } catch (_) { return null; }
  if (!isObject(parsed) || typeof parsed.recordId !== 'string' || !isCapacityKind(parsed.labelKind)
    || !isCount(parsed.labelCount) || !isNumber(parsed.fontSize) || !isRatio(parsed.capacityRatio)) return null;
  return Object.freeze({
    recordId: parsed.recordId,
    labelKind: parsed.labelKind,
    labelCount: parsed.labelCount,
    fontSize: parsed.fontSize,
    capacityRatio: parsed.capacityRatio
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
