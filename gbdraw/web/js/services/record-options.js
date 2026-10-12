// @ts-check
import { circularRecordsToDraw } from './record-draw-selection.js';

/** @param {unknown} value */
const cleanText = (value) => String(value ?? '').trim();
const NULL_RECORD_SELECTOR_TOKENS = new Set([
  'none', 'null', 'jsnull', 'undefined', 'jsundefined', '-'
]);

/** @param {unknown} value */
export const isUnspecifiedRecordSelectorValue = (value) => {
  const text = cleanText(value);
  return !text || NULL_RECORD_SELECTOR_TOKENS.has(text.toLowerCase());
};

/** @param {unknown} value */
const isSafeRecordId = (value) => {
  const recordId = cleanText(value);
  return !isUnspecifiedRecordSelectorValue(recordId) && !recordId.startsWith('#');
};

/** @param {unknown} value */
export const formatRecordLength = (value) => {
  const numeric = Number(value);
  if (!Number.isInteger(numeric) || numeric <= 0) return 'length unavailable';
  return `${numeric.toLocaleString('en-US')} bp`;
};

/**
 * A discovered record as these functions read it; any other field is kept.
 * @typedef {{ sourceIndex?: unknown, selector?: unknown, recordId?: unknown }} RecordEntryInput
 */

/**
 * @template {RecordEntryInput} T
 * @param {readonly T[] | null | undefined} records
 */
export const buildDisambiguatedRecordEntries = (records) => {
  const normalized = (Array.isArray(records) ? records : []).map((record, index) => ({
    ...record,
    sourceIndex: Number.isInteger(record?.sourceIndex) ? record.sourceIndex : index,
    selector: cleanText(record?.selector) || `#${index + 1}`,
    recordId: cleanText(record?.recordId) || `Record_${index + 1}`
  }));
  const idCounts = new Map();
  normalized.forEach((record) => {
    idCounts.set(record.recordId, (idCounts.get(record.recordId) || 0) + 1);
  });
  return normalized.map((record) => {
    const duplicate = (idCounts.get(record.recordId) || 0) > 1;
    const usesIndex = duplicate || !isSafeRecordId(record.recordId);
    return {
      ...record,
      duplicate,
      usesIndex,
      value: usesIndex ? record.selector : record.recordId
    };
  });
};

/**
 * @template {RecordEntryInput} T
 * @param {readonly T[] | null | undefined} records
 * @param {unknown} requestedValue
 */
export const resolveDisambiguatedRecordSelection = (records, requestedValue) => {
  const entries = buildDisambiguatedRecordEntries(records);
  const requested = cleanText(requestedValue);
  if (!requested) return { status: 'unspecified', requested, entries, record: null };

  for (const field of ['value', 'selector', 'recordId']) {
    const matches = entries.filter((record) => record[field] === requested);
    if (matches.length === 1) {
      return { status: 'resolved', requested, entries, record: matches[0] };
    }
    if (matches.length > 1) {
      return { status: 'ambiguous', requested, entries, record: null };
    }
  }
  return { status: 'missing', requested, entries, record: null };
};

// The source records one Circular request draws: the selected record of a
// single presentation, otherwise every ON record (record selection; an
// explicit single choice draws its record ON or OFF). services/session-request.js
// builds the request from this set and the region annotation record catalog
// offers the same records. `presentedCount` is how many records the
// presentation offers (1 for a single choice), which decides single or batch;
// `omittedRecords` are the offered records that are OFF.
/**
 * @typedef {Object} CircularRequestRecordSetOptions
 * @property {readonly (RecordEntryInput & { record_id?: unknown })[] | null} [records] Discovered records (`record_id` or `recordId`, `recordLength`, `selector`).
 * @property {string} [selector]
 * @property {boolean} [multiRecordCanvas]
 * @property {string} [groupingIntent]
 * @property {readonly string[]} [recordsOff] The drawing's OFF records (`#N`).
 */

/** @param {CircularRequestRecordSetOptions} [options] */
export const resolveCircularRequestRecordSet = ({
  records,
  selector = '',
  multiRecordCanvas = false,
  groupingIntent = '',
  recordsOff = []
} = {}) => {
  const recordSelectors = buildDisambiguatedRecordEntries(
    (Array.isArray(records) ? records : []).map((record) => ({
      ...record,
      recordId: record?.record_id ?? record?.recordId
    }))
  );
  const requestedSelector = String(selector || '').trim();
  const selection = resolveDisambiguatedRecordSelection(recordSelectors, requestedSelector);
  const singlePresentation = !multiRecordCanvas && groupingIntent !== 'batch';
  // 'AMBIGUOUS' or 'NO_MATCH' (the RECORD_SELECTION reasons) when a single
  // presentation names a record the input does not resolve to exactly one.
  const selectionFailure = singlePresentation && requestedSelector && selection.status !== 'resolved'
    ? (selection.status === 'ambiguous' ? 'AMBIGUOUS' : 'NO_MATCH')
    : '';
  if (singlePresentation && selection.record) {
    return { recordSelectors, records: [selection.record], omittedRecords: [], presentedCount: 1,
      singlePresentation, selectionFailure };
  }
  const drawnRecords = circularRecordsToDraw(recordSelectors, recordsOff);
  const drawn = new Set(drawnRecords);
  return {
    recordSelectors,
    records: drawnRecords,
    omittedRecords: recordSelectors.filter((record) => !drawn.has(record)),
    presentedCount: recordSelectors.length,
    singlePresentation,
    selectionFailure: selectionFailure || (recordSelectors.length > 0 && drawnRecords.length === 0 ? 'NONE_DRAWN' : '')
  };
};

/**
 * What the drawn-record reader takes from a Circular drawing.
 * @typedef {object} CircularRecordDrawing
 * @property {{ circular_record_selector: string, multi_record_canvas: boolean,
 *   circular_region_start: unknown, circular_region_end: unknown }} form
 * @property {{ circular_grouping_intent: string }} adv
 * @property {readonly string[]} recordsOff
 */

/**
 * Whether a record set is one record offered alone: the journey that draws the
 * record's own presentation and can crop it to a region.
 * @param {{ presentedCount: number, singlePresentation: boolean }} recordSet
 */
export const isSingleCircularJourney = (recordSet) => recordSet.presentedCount <= 1 && recordSet.singlePresentation;

/**
 * The record lengths of each diagram the next Generate of a Circular drawing
 * draws, in drawing order: a Multi-Record Canvas is one diagram of every drawn
 * record, and any other record is a diagram of its own. Empty while a drawn
 * length is not known: no record is read or drawn, the selection fails, or
 * the one record is cropped to a region.
 * @param {readonly { selector?: string, record_id?: string, record_length?: number | null }[]} records
 * @param {CircularRecordDrawing} drawing
 * @returns {number[][]}
 */
export const circularDiagramRecordLengths = (records, drawing) => {
  const recordSet = resolveCircularRequestRecordSet({
    records: [...records],
    selector: drawing.form.circular_record_selector,
    multiRecordCanvas: drawing.form.multi_record_canvas,
    groupingIntent: drawing.adv.circular_grouping_intent,
    recordsOff: drawing.recordsOff
  });
  const cropped = isSingleCircularJourney(recordSet)
    && Boolean(cleanText(drawing.form.circular_region_start) || cleanText(drawing.form.circular_region_end));
  const lengths = recordSet.records.map((record) => Number(record.record_length));
  if (cropped || recordSet.selectionFailure || lengths.length === 0 || !lengths.every((length) => length > 0)) return [];
  return drawing.form.multi_record_canvas ? [lengths] : lengths.map((length) => [length]);
};
