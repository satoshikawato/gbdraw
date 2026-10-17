// @ts-check
import { circularRecordRequestKeys, circularRecordsToDraw } from './record-draw-selection.js';

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
 * @property {readonly (RecordEntryInput & { record_id?: unknown, record_length?: unknown })[] | null} [records] Discovered records (`record_id` or `recordId`, `record_length`, `selector`).
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
 * The region fields of the Circular form.
 * @typedef {object} CircularRegionForm
 * @property {unknown} circular_region_start
 * @property {unknown} circular_region_end
 */

/**
 * What the record-set readers take from a Circular drawing.
 * @typedef {object} CircularRecordDrawing
 * @property {CircularRegionForm & { circular_record_selector: string, multi_record_canvas: boolean }} form
 * @property {{ circular_grouping_intent: string }} adv
 * @property {readonly string[]} recordsOff
 */

/**
 * The source records the next Generate of a Circular drawing draws: the one
 * mapping from a drawing to the record selection. The request, the region
 * annotation record catalog, the Auto values and the track slot estimates
 * read it.
 * @param {CircularRequestRecordSetOptions['records']} records the records discovered in the input
 * @param {CircularRecordDrawing} drawing
 */
export const circularDrawingRecordSet = (records, drawing) => resolveCircularRequestRecordSet({
  records,
  selector: drawing.form.circular_record_selector,
  multiRecordCanvas: drawing.form.multi_record_canvas,
  groupingIntent: drawing.adv.circular_grouping_intent,
  recordsOff: drawing.recordsOff
});

/**
 * Whether a record set is one record offered alone: the journey that draws the
 * record's own presentation and can crop it to a region.
 * @param {{ presentedCount: number, singlePresentation: boolean }} recordSet
 */
export const isSingleCircularJourney = (recordSet) => recordSet.presentedCount <= 1 && recordSet.singlePresentation;

/**
 * The region the Circular form crops a record to (1-based, both ends
 * included), as its fields state it: null without one (both fields empty), and
 * `{ error }` while Generate cannot draw it. The one reader of those fields;
 * only a record offered alone is cropped (`isSingleCircularJourney`).
 * @param {CircularRegionForm} form
 * @param {{ recordLength?: unknown, record_length?: unknown } | null} [record] the record it crops, when known
 * @returns {{ start: number, end: number } | { error: string } | null}
 */
export const circularFormRegion = (form, record = null) => {
  const hasStart = cleanText(form.circular_region_start) !== '';
  const hasEnd = cleanText(form.circular_region_end) !== '';
  if (!hasStart && !hasEnd) return null;
  if (hasStart !== hasEnd) return { error: 'Circular region requires both Start and End coordinates.' };
  const start = Number(form.circular_region_start);
  const end = Number(form.circular_region_end);
  if (!Number.isInteger(start) || !Number.isInteger(end) || start < 1 || end < 1) {
    return { error: 'Circular region Start and End must be positive integers.' };
  }
  if (start > end) {
    return { error: 'Circular region Start must not exceed End. Use Reverse complement to change display orientation.' };
  }
  const recordLength = Number(record?.recordLength ?? record?.record_length);
  if (Number.isInteger(recordLength) && recordLength > 0 && end > recordLength) {
    return { error: `Circular region End (${end}) exceeds the selected record length (${recordLength}).` };
  }
  return { start, end };
};

/**
 * The region Generate crops a record to: the form's region when the drawing
 * offers that record alone and Generate can draw the region; null otherwise
 * (a region Generate refuses crops nothing).
 * @param {boolean} single whether the drawing offers one record alone (`isSingleCircularJourney`)
 * @param {CircularRegionForm} form
 * @param {{ recordLength?: unknown, record_length?: unknown } | null} record the record it crops, when known
 * @returns {{ start: number, end: number } | null}
 */
const circularCropRegion = (single, form, record) => {
  const region = single ? circularFormRegion(form, record) : null;
  return region !== null && !('error' in region) ? region : null;
};

/**
 * Whether the next Generate of a Circular drawing crops its record
 * (`circularCropRegion`).
 * @param {CircularRequestRecordSetOptions['records']} records the records discovered in the input
 * @param {CircularRecordDrawing} drawing
 */
export const circularDrawingCropped = (records, drawing) => {
  const recordSet = circularDrawingRecordSet(records, drawing);
  return circularCropRegion(isSingleCircularJourney(recordSet), drawing.form, recordSet.records[0] ?? null) !== null;
};

/**
 * One diagram of the next Circular Generate.
 * @typedef {object} CircularDiagramRecords
 * @property {string[]} recordKeys the request keys that can name its records
 * @property {number[]} lengths the drawn lengths of its records, where known
 */

/**
 * The diagrams the next Generate of a Circular drawing draws, in drawing
 * order: a Multi-Record Canvas is one diagram of every drawn record, and any
 * other record is a diagram of its own. A record is drawn at its own length,
 * or at the length of the region it is cropped to (Generate crops the record
 * before it sizes the diagram); a length that is not known is left out, and
 * so is a diagram with no known length. While the drawing draws no record,
 * the diagrams are those of the records it offers.
 * @param {CircularRequestRecordSetOptions['records']} records the records discovered in the input
 * @param {CircularRecordDrawing} drawing
 * @returns {CircularDiagramRecords[]}
 */
export const circularRequestDiagrams = (records, drawing) => {
  const recordSet = circularDrawingRecordSet(records, drawing);
  const single = isSingleCircularJourney(recordSet);
  const drawn = recordSet.records.length > 0 ? recordSet.records : recordSet.recordSelectors;
  // A discovered record keeps its fields through the record set.
  /** @param {{ record_length?: unknown, sourceIndex: number, recordKey?: unknown, recordId?: unknown, selector?: unknown }} record */
  const drawnLength = (record) => {
    const region = circularCropRegion(single, drawing.form, record);
    return region !== null ? region.end - region.start + 1 : Number(record.record_length);
  };
  return (drawing.form.multi_record_canvas ? [drawn] : drawn.map((record) => [record]))
    .map((diagram) => ({
      recordKeys: diagram.flatMap(circularRecordRequestKeys),
      lengths: diagram.map(drawnLength).filter((length) => length > 0)
    }))
    .filter((diagram) => diagram.lengths.length > 0);
};

/**
 * What the estimates read of the request that drew the shown Results.
 * @typedef {object} ShownRequest
 * @property {string} mode
 * @property {string} grouping `batch` draws one Result per record, in request order; any other grouping one Result
 * @property {readonly { recordKey: string }[]} records
 */

/**
 * The record lengths the Circular track slot estimates read: those of the
 * diagram of the next Generate that draws a record of the shown Result. The
 * record is found by its request key, not by the Result's position, since the
 * next Generate can draw other records. The first diagram when no Circular
 * Result is shown or the next Generate draws none of its records; empty while
 * no length is known.
 * @param {CircularDiagramRecords[]} diagrams
 * @param {ShownRequest | null} shownRequest
 * @param {number} resultIndex the shown Result
 * @returns {number[]}
 */
export const circularEstimateLengths = (diagrams, shownRequest, resultIndex) => {
  const circularRecords = shownRequest !== null && shownRequest.mode === 'circular' ? shownRequest.records : [];
  const shownRecords = shownRequest !== null && shownRequest.grouping === 'batch'
    ? circularRecords.slice(resultIndex, resultIndex + 1)
    : circularRecords;
  const shownKeys = new Set(shownRecords.map((record) => record.recordKey));
  const diagram = diagrams.find(({ recordKeys }) => recordKeys.some((key) => shownKeys.has(key))) || diagrams[0];
  return diagram ? diagram.lengths : [];
};
