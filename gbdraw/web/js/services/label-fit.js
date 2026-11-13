// @ts-check
// The label fit report Python sends (gbdraw/web_support/label_fit.py): the one
// JavaScript reader of it. Python owns every value of the report and the words of
// each record's counts text, the rows' titles, and the settings a row changes by
// config path and by request track option, with the canvas sizes (LP-D43, LP-D48). Nothing here places labels or searches values; a row that was not
// computed is computed by the Worker (`computeLabelFitRow`, OD-F14).
// The exceptions are the Web's own words for labels that need much more room
// than the map has: the ways to make room (CAPACITY-STOP, NO-ROW-WAYS), which Python names
// by id while the Web names its controls, and the room such labels need
// (CAPACITY-STOP, CAPACITY-WARN), which Python sends as values (OV-2179).
// A row's title and the names and option words of its settings are English
// catalog keys: `tools/generate_web_contracts.py` lists them from Python for
// the catalogs (`tools/i18n/source-keys.mjs`), and the dialog shows them with `tm()`.
// Why a row has no value and what a row that draws with the notice leaves are
// Python's data, which the Web words in the shown language: each reason kind
// with its numbers, each count by its field name (OV-3109).
import { diagnosticError, isFigureLengthText } from '../utils/error-normalization.js';
import { listText, msg, numberText, t, tn } from '../utils/i18n.js';
import {
  LABEL_FIT_FIELDS, LABEL_FIT_LIMIT_REASON_KINDS, LABEL_FIT_OFFERED_ROW_STATUSES, LABEL_FIT_PLAIN_REASON_KINDS,
  LABEL_FIT_REGION_WARNING_FIELDS
} from '../contracts/schema.generated.js';
/** @import { CanonicalCircularTrackOptions } from '../contracts/request.generated.js' */

/**
 * @typedef {'found' | 'found_with_notice' | 'not_found' | 'not_computed'} LabelFitRowStatus
 * `found_with_notice`: the row leaves no overlapping labels but leaves crossing
 * leaders, labels over the definition text or region label warnings (LP-D57),
 * or a Circular label row fits the labels it moves while each record keeps
 * region label warnings it already had (OV-3105, D-112, D-115), so the figure
 * draws with the notice bar; its `remaining` says what the figure is drawn with.
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
/** @typedef {typeof LABEL_FIT_FIELDS[number]} LabelFitField A count of a label fit, by its Python field name. */
/**
 * @typedef {object} LabelFitCount One count of a label fit (Python's `LabelFitCountPayload`).
 * @property {LabelFitField} field
 * @property {number} count
 */
/**
 * @typedef {typeof LABEL_FIT_LIMIT_REASON_KINDS[number]} LabelFitLimitReasonKind A lever's limit (an offset, a font
 *   size in pt, a label spacing in px), or `search_bound`: the Worker's search for the row stopped after `limit` seconds.
 */
/**
 * @typedef {typeof LABEL_FIT_PLAIN_REASON_KINDS[number]} LabelFitPlainReasonKind A reason without numbers.
 */
/**
 * Why a row has no value (Python's `LabelFitRowReasonPayload`), by `kind`: a lever's `limit`; the Combination row's
 * `fontSize` and Y offset `limit`; the `counts` a switch row's assembly leaves (none to state: empty); or none.
 * @typedef {{ kind: LabelFitLimitReasonKind, limit: number }
 *   | { kind: 'combination_limit', fontSize: number, limit: number }
 *   | { kind: 'labels_do_not_fit', counts: readonly LabelFitCount[] }
 *   | { kind: LabelFitPlainReasonKind }} LabelFitRowReason
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
 * @property {LabelFitRowReason | null} reason Why the row has no value; worded by `labelFitRowReasonText`.
 * @property {readonly LabelFitCount[]} remaining What a `found_with_notice` row leaves, in its note's order (empty for
 *   another row); worded by `labelFitRowNoteText`.
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
 * @property {readonly LabelCapacityRemedy[]} remedies The ways to make room, in Python's order: for the records not
 *   placed (CAPACITY-STOP), or when no row has a value (NO-ROW-WAYS; the dialog lists them once every row is known
 *   not to fit).
 * @property {LargerFigure | null} largerFigure The requested sides the `figure_size` remedy names, each with Python's
 *   estimate (CAPACITY-STOP) or without one (true, NO-ROW-WAYS), a side left Auto null; null unless the remedy is
 *   offered.
 */
/**
 * A requested side of a figure to make larger: about how long it must be, as Python's message writes it (`829 mm`),
 * for labels not placed (CAPACITY-STOP); true when no row fits and Python cannot estimate it (NO-ROW-WAYS).
 * @typedef {string | true} LargerSide
 */
/**
 * The requested sides of a figure to make larger (Python's `LargerFigurePayload`): null for a side left Auto; at
 * least one side, and both requested sides have an estimate or neither has.
 * @typedef {{ width: LargerSide | null, height: LargerSide | null }} LargerFigure
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
  figure_size: msg('Make the figure larger (Figure size in Layout).'),
  linear: msg('Draw a Linear map.')
});
/** @typedef {keyof typeof LABEL_CAPACITY_REMEDY_WORDS} LabelCapacityRemedy */
// The `figure_size` remedy with Python's estimate, by the sides it names (D-14).
const LARGER_FIGURE_WORDS = Object.freeze({
  width: msg('Make the figure larger: Figure width about {width} (Figure size in Layout).'),
  height: msg('Make the figure larger: Figure height about {height} (Figure size in Layout).'),
  both: msg('Make the figure larger: Figure width about {width}, Figure height about {height} (Figure size in Layout).')
});
// The `figure_size` remedy without an estimate (NO-ROW-WAYS), by the sides it names (PD-OI-114 r3).
const LARGER_FIGURE_FIELD_WORDS = Object.freeze({
  width: msg('Make the figure larger: increase Figure width (Figure size in Layout).'),
  height: msg('Make the figure larger: increase Figure height (Figure size in Layout).'),
  both: msg('Make the figure larger: increase Figure width and Figure height (Figure size in Layout).')
});
/**
 * A way to make room as the dialog shows it: a catalog key and its values.
 * @typedef {{ text: string, params: Readonly<Record<string, string>> }} LabelCapacityRemedyText
 */
/**
 * The words of a way to make room; the `figure_size` remedy names each requested side, with Python's estimate when
 * it has one.
 * @param {LabelCapacityRemedy} remedy
 * @param {LargerFigure | null} largerFigure
 * @returns {LabelCapacityRemedyText}
 */
export const labelCapacityRemedyText = (remedy, largerFigure) => {
  if (remedy !== 'figure_size' || largerFigure === null) return { text: LABEL_CAPACITY_REMEDY_WORDS[remedy], params: {} };
  const { width, height } = largerFigure;
  const sides = width !== null && height !== null ? 'both' : width !== null ? 'width' : 'height';
  if (width === true || height === true) return { text: LARGER_FIGURE_FIELD_WORDS[sides], params: {} };
  return { text: LARGER_FIGURE_WORDS[sides], params: { ...(width !== null ? { width } : {}), ...(height !== null ? { height } : {}) } };
};
/** A side: null, true, or a length as Python's message writes it. @param {unknown} value @returns {value is LargerSide | null} */
const isLargerSide = (value) => value === null || value === true || isFigureLengthText(value);
/**
 * A report's `largerFigure`, absent (null) unless the `figure_size` remedy is offered.
 * @param {Record<string, unknown>} parsed
 * @returns {LargerFigure | null}
 */
const normalizeLargerFigure = (parsed) => {
  if (!Object.hasOwn(parsed, 'largerFigure')) return null;
  const value = parsed.largerFigure;
  if (!isObject(value) || Object.keys(value).sort().join(',') !== 'height,width' || !isLargerSide(value.width)
    || !isLargerSide(value.height) || (value.width === null && value.height === null)
    // Python estimates both requested sides or neither.
    || (value.width !== null && value.height !== null && (value.width === true) !== (value.height === true))
    || !Array.isArray(parsed.remedies) || !parsed.remedies.includes('figure_size')) throw invalid();
  return Object.freeze({ width: value.width, height: value.height });
};

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

/**
 * Each count of a label fit as one item of a list, by Python's field name (`LABEL_FIT_FIELDS`): the counts of a row
 * whose one assembly leaves labels that do not fit, and the region label warnings a row draws with.
 * @type {Readonly<Record<LabelFitField, (count: number) => string>>}
 */
const LABEL_FIT_COUNT_WORDS = Object.freeze({
  overlap_pairs: (count) => tn(count, '{count} label pair overlaps', '{count} label pairs overlap'),
  crossing_pairs: (count) => tn(count, '{count} leader pair crosses', '{count} leader pairs cross'),
  definition_overlaps: (count) => tn(count, '{count} label covers the definition text', '{count} labels cover the definition text'),
  region_overlaps: (count) => tn(count, '{count} region label overlap', '{count} region label overlaps'),
  region_cuts: (count) => tn(count, '{count} region label cut', '{count} region labels cut'),
  region_track_overlaps: (count) => tn(count, '{count} region label on another track', '{count} region labels on other tracks'),
  definition_track_overlaps: (count) => tn(count, '{count} line overlaps a track', '{count} lines overlap a track'),
  region_rooms_past_centre: (count) => tn(count, '{count} region label left out for want of room before the centre',
    '{count} region labels left out for want of room before the centre'),
  region_rooms_clipped: (count) => tn(count, '{count} region label left out in clipped slots',
    '{count} region labels left out in clipped slots'),
  region_cuts_left_out: (count) => tn(count, '{count} region label left out: no letter fits its circle',
    '{count} region labels left out: no letter fits their circles')
});
/**
 * What a row that draws with the notice still leaves, other than region label warnings, by field name.
 * @type {Readonly<Partial<Record<LabelFitField, (count: number) => string>>>}
 */
const LABEL_FIT_NOTE_WORDS = Object.freeze({
  crossing_pairs: (count) => tn(count, '{count} leader pair still crosses', '{count} leader pairs still cross'),
  definition_overlaps: (count) => tn(count, '{count} label still covers the definition text',
    '{count} labels still cover the definition text'),
  definition_track_overlaps: (count) => tn(count, '{count} definition line still overlaps a track',
    '{count} definition lines still overlap a track')
});
/** @type {ReadonlySet<string>} */
const REGION_WARNING_FIELDS = new Set(LABEL_FIT_REGION_WARNING_FIELDS);

/**
 * The note of a row that draws with the notice: what it still leaves, with the region label warnings the figure is
 * drawn with; null when it leaves nothing.
 * @param {readonly LabelFitCount[]} remaining
 * @returns {string | null}
 */
export const labelFitRowNoteText = (remaining) => {
  const regions = remaining.filter(({ field }) => REGION_WARNING_FIELDS.has(field));
  const clauses = remaining.flatMap(({ field, count }) => LABEL_FIT_NOTE_WORDS[field]?.(count) ?? []);
  if (regions.length) {
    clauses.push(t('the figure is drawn with {warnings}',
      { warnings: listText(regions.map(({ field, count }) => LABEL_FIT_COUNT_WORDS[field](count))) }));
  }
  return clauses.length ? listText(clauses) : null;
};

/**
 * Why a row has no value, in the shown language.
 * @param {LabelFitRowReason} reason
 * @returns {string}
 */
export const labelFitRowReasonText = (reason) => {
  switch (reason.kind) {
    case 'y_offset_limit': return t('Does not fit up to Y offset {limit}', { limit: numberText(reason.limit) });
    case 'x_offset_limit': return t('Does not fit up to X offset {limit}', { limit: numberText(reason.limit) });
    case 'xy_offset_limit': return t('Does not fit up to X and Y offset {limit}', { limit: numberText(reason.limit) });
    case 'font_limit': return t('Does not fit down to {size} pt', { size: numberText(reason.limit) });
    case 'spacing_limit': return t('Does not fit down to {spacing} px label spacing', { spacing: numberText(reason.limit) });
    case 'search_bound': return tn(reason.limit, 'No value found within {count} second.',
      'No value found within {count} seconds.', { count: numberText(reason.limit) });
    case 'combination_limit': return t('Does not fit at {size} pt up to Y offset {limit}',
      { size: numberText(reason.fontSize), limit: numberText(reason.limit) });
    case 'labels_do_not_fit': return reason.counts.length
      ? t('Labels still do not fit ({counts})',
        { counts: listText(reason.counts.map(({ field, count }) => LABEL_FIT_COUNT_WORDS[field](count))) })
      : t('Labels still do not fit');
    case 'outer_labels_only': return t('Does not move inner, embedded, or radial labels');
    case 'tracks_do_not_fit': return t('The tracks do not fit with this setting');
    case 'labels_not_placed': return t('The labels cannot be placed with this setting');
    case 'stale': return t('Generate again to compute this choice.');
  }
};

const ROW_STATUSES = new Set(['found', 'found_with_notice', 'not_found', 'not_computed']);
/** @type {ReadonlyMap<string, number>} Python's offered statuses, by preference (0 first). */
const OFFERED_STATUS_RANK = new Map(LABEL_FIT_OFFERED_ROW_STATUSES.map((status, rank) => [status, rank]));

/**
 * Whether a row's status offers its settings to Apply.
 * @param {string} status
 * @returns {boolean}
 */
export const labelFitRowOffered = (status) => OFFERED_STATUS_RANK.has(status);

/**
 * How Python prefers an offered row's status: 0 for a row that fits, then a row that draws with the notice.
 * @param {string} status
 * @returns {number}
 */
export const labelFitRowRank = (status) => OFFERED_STATUS_RANK.get(status) ?? OFFERED_STATUS_RANK.size;

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
/** @type {ReadonlySet<string>} */
const FIELDS = new Set(LABEL_FIT_FIELDS);
/** @type {ReadonlySet<string>} */
const LIMIT_REASON_KINDS = new Set(LABEL_FIT_LIMIT_REASON_KINDS);
/** @type {ReadonlySet<string>} */
const PLAIN_REASON_KINDS = new Set(LABEL_FIT_PLAIN_REASON_KINDS);

/**
 * Counts of a label fit, each a field Python names with its count.
 * @param {unknown} value
 * @returns {readonly LabelFitCount[]}
 */
const normalizeCounts = (value) => {
  if (!Array.isArray(value)) throw invalid();
  return Object.freeze(value.map((item) => {
    if (!isObject(item) || typeof item.field !== 'string' || !FIELDS.has(item.field) || !isCount(item.count)) throw invalid();
    return Object.freeze({ field: /** @type {LabelFitField} */ (item.field), count: item.count });
  }));
};

/**
 * A row's reason, with the fields its kind has.
 * @param {unknown} value
 * @returns {LabelFitRowReason | null}
 */
const normalizeReason = (value) => {
  if (value === null) return null;
  if (!isObject(value) || typeof value.kind !== 'string') throw invalid();
  const { kind } = value;
  if (LIMIT_REASON_KINDS.has(kind) && isNumber(value.limit)) {
    return Object.freeze({ kind: /** @type {LabelFitLimitReasonKind} */ (kind), limit: value.limit });
  }
  if (kind === 'combination_limit' && isNumber(value.fontSize) && isNumber(value.limit)) {
    return Object.freeze({ kind, fontSize: value.fontSize, limit: value.limit });
  }
  if (kind === 'labels_do_not_fit') return Object.freeze({ kind, counts: normalizeCounts(value.counts) });
  if (PLAIN_REASON_KINDS.has(kind)) return Object.freeze({ kind: /** @type {LabelFitPlainReasonKind} */ (kind) });
  throw invalid();
};

/**
 * One row, as the report or the Worker's `computeLabelFitRow` gives it.
 * @param {unknown} value
 * @returns {LabelFitRow}
 */
export const normalizeLabelFitRow = (value) => {
  if (!isObject(value) || typeof value.kind !== 'string' || !value.kind
    || typeof value.title !== 'string' || typeof value.status !== 'string'
    || !ROW_STATUSES.has(value.status) || !Array.isArray(value.settings)
    || !isObject(value.trackOptions) || !Object.values(value.trackOptions).every(isNumber)) throw invalid();
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
    reason: normalizeReason(value.reason),
    remaining: normalizeCounts(value.remaining)
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
    remedies: Object.freeze([...parsed.remedies]),
    largerFigure: normalizeLargerFigure(parsed)
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
