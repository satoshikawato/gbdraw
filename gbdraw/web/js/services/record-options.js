// @ts-check
/** @import { LinearSeqWindow } from './drawing/inputs.js' */
import { circularRecordsToDraw } from './record-draw-selection.js';
/** @import { CircularListRecord } from './record-discovery.js' */
import { diagnosticError } from '../utils/error-normalization.js';
import { numberText } from '../utils/i18n.js';

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
  return `${numberText(numeric)} bp`;
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
export const resolveDisambiguatedRecordSelection = (records, requestedValue) => (
  selectAmongRecordEntries(buildDisambiguatedRecordEntries(records), requestedValue)
);

/**
 * The entry a requested value names among named entries: by `value`, then
 * `selector`, then `recordId`; a value two entries share is ambiguous.
 * @template {{ value: string, selector: string, recordId: string }} E
 * @param {readonly E[]} entries
 * @param {unknown} requestedValue
 * @returns {{ status: 'unspecified' | 'resolved' | 'ambiguous' | 'missing', requested: string, entries: readonly E[], record: E | null }}
 */
const selectAmongRecordEntries = (entries, requestedValue) => {
  const requested = cleanText(requestedValue);
  if (!requested) return { status: 'unspecified', requested, entries, record: null };

  for (const field of /** @type {const} */ (['value', 'selector', 'recordId'])) {
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

/**
 * The Circular records of a drawing (`circularRecordList`, in drawing order),
 * named as the Record select and the request name them: each is disambiguated
 * within its own card's file (a record ID twice in one file, or an unsafe one,
 * is named by its `#k`), `fileValue` is that name (the selector of the CLI's
 * form a request record of its file carries), and `value` is its identity,
 * which the Record select stores.
 * @template {Partial<CircularListRecord> & { recordId?: unknown, recordLength?: unknown }} T
 * @param {readonly T[] | null | undefined} records
 */
export const circularRecordEntries = (records) => {
  /** @type {Map<string, T[]>} */
  const byCard = new Map();
  (Array.isArray(records) ? records : []).forEach((record) => {
    const uid = cleanText(record?.uid);
    byCard.set(uid, [...(byCard.get(uid) || []), record]);
  });
  return [...byCard.values()].flatMap((list) => buildDisambiguatedRecordEntries(list.map((record, index) => ({
    ...record,
    sourceIndex: Number.isInteger(record?.k) ? Number(record.k) - 1 : index,
    recordId: record?.record_id ?? record?.recordId,
    recordLength: record?.record_length ?? record?.recordLength
  }))).map((entry) => ({ ...entry, fileValue: entry.value, value: cleanText(entry.identity) || entry.value })));
};

/**
 * The Circular record a Record select value names (`circularRecordEntries`):
 * its identity, or a selector or record ID that names one record.
 * @template {{ value: string, selector: string, recordId: string }} E
 * @param {readonly E[]} entries `circularRecordEntries`
 * @param {unknown} requestedValue
 */
export const circularRecordSelection = (entries, requestedValue) => selectAmongRecordEntries(entries, requestedValue);

// The source records one Circular request draws: the selected record of a
// single presentation, otherwise every ON record (record selection; an
// explicit single choice draws its record ON or OFF). services/session-request.js
// builds the request from this set and the region annotation record catalog
// offers the same records. `presentedCount` is how many records the
// presentation offers (1 for a single choice), which decides single or batch;
// `omittedRecords` are the offered records that are OFF.
/**
 * @typedef {Object} CircularRequestRecordSetOptions
 * @property {readonly (Partial<CircularListRecord> & { recordId?: unknown, recordLength?: unknown })[] | null} [records] Discovered records (`circularRecordList`).
 * @property {string} [selector]
 * @property {boolean} [multiRecordCanvas]
 * @property {string} [groupingIntent]
 * @property {readonly string[]} [recordsOff] The drawing's OFF records (their identities).
 */

/** @param {CircularRequestRecordSetOptions} [options] */
export const resolveCircularRequestRecordSet = ({
  records,
  selector = '',
  multiRecordCanvas = false,
  groupingIntent = '',
  recordsOff = []
} = {}) => {
  const recordSelectors = circularRecordEntries(records);
  const requestedSelector = String(selector || '').trim();
  const selection = circularRecordSelection(recordSelectors, requestedSelector);
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
 * The crop fields of the Circular form: the region its Start/End fields state, and the window list that holds what
 * those fields cannot (one window across the origin of a circular record, which a CLI or Python Session wrote).
 * @typedef {object} CircularRegionForm
 * @property {unknown} circular_region_start
 * @property {unknown} circular_region_end
 * @property {readonly LinearSeqWindow[]} circular_windows
 */

/**
 * The crop a CLI or Python Session gave one record of a Circular grid or batch (PD-OI-044 scenario revision 5),
 * which the drawing keeps while it draws that record in a grid or batch. `source` is the uid of the record's card;
 * `record` names the record in that card's file as the Session's request did (its identity, its record ID or `#k`,
 * or '' for the record of a file of one record), so it is read once the card's records are (`circularRecordRange`);
 * `windows` are its windows in source coordinates, each with its own orientation.
 * @typedef {{ source: string, record: string, windows: LinearSeqWindow[] }} CircularRecordRange
 */

/**
 * The range a drawing keeps for one of its read records (`circularRecordEntries`): the row of the record's card
 * that names it there (`circularRecordSelection` over that card's records), else null. The one reader of a
 * record's range: the request, the diagram lengths, the Record select and the record rows read it.
 * @template {{ uid?: unknown, value: string, selector: string, recordId: string }} E
 * @param {readonly CircularRecordRange[]} ranges
 * @param {readonly E[]} entries the drawing's read records, named
 * @param {E} entry
 * @returns {CircularRecordRange | null}
 */
export const circularRecordRange = (ranges, entries, entry) => {
  if (!ranges.length) return null;
  const uid = cleanText(entry.uid);
  const card = entries.filter((candidate) => cleanText(candidate.uid) === uid);
  return ranges.find((range) => range.source === uid
    && circularRecordSelection(card, range.record || (card.length === 1 ? card[0].value : '')).record?.value === entry.value)
    ?? null;
};

/**
 * The ranges of a drawing without those of the records of card `uid`: what a Remove of the card, or a replacement
 * of its file that retires its records' presentation, leaves.
 * @param {readonly CircularRecordRange[]} ranges @param {string} uid
 */
export const circularRangesWithoutCard = (ranges, uid) => ranges.filter((range) => range.source !== uid);

/**
 * What the record-set readers take from a Circular drawing.
 * @typedef {object} CircularRecordDrawing
 * @property {CircularRegionForm & { circular_record_selector: string, multi_record_canvas: boolean,
 *   circular_record_ranges: readonly CircularRecordRange[] }} form
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
 * A region of a record (1-based, both ends included): `start` and `end` as
 * typed, and `drawnEnd`, the End that Generate draws to. `recordLength` is set
 * when `drawnEnd` is shorter than `end`.
 * @typedef {{ start: number, end: number, drawnEnd: number, recordLength?: number }} RecordRegion
 */

/**
 * Where Generate's crop ends: Python owns the rule (`check_start_end_coords`
 * in gbdraw/crop_genbank.py) and this is the Web's one prediction of it. An End
 * past the end of a record of known length is drawn to the record end; a
 * region that starts after the record ends is refused. Without a known length
 * the region is drawn as typed, and Python applies the rule.
 * @param {number} start
 * @param {number} end not less than `start`
 * @param {unknown} recordLength
 * @returns {RecordRegion | null} null when the region starts after the record ends
 */
export const regionWithinRecord = (start, end, recordLength) => {
  const length = Number(recordLength);
  if (!Number.isInteger(length) || length <= 0 || end <= length) return { start, end, drawnEnd: end };
  return start > length ? null : { start, end, drawnEnd: length, recordLength: length };
};

/**
 * The bases a window draws, both ends included. A window whose start is after its end crosses the origin of a
 * circular record (D-34): its bases need the record's length, and are null until the record is read (Python's
 * `region_slice`). The one reader of a window's length for the window lists, a Linear card's window across the
 * origin, and the Circular crop.
 * @param {number} start @param {number} end
 * @param {number | null} recordLength the length of the record the window draws, once read
 */
export const windowBases = (start, end, recordLength) => (
  start <= end ? end - start + 1 : recordLength === null ? null : recordLength - start + 1 + end
);

/**
 * One row of the read-only window list of a record drawn in several windows (D-10): its number (1-based, in column
 * order), its span in source coordinates, its length, and whether it is drawn reverse-complemented.
 * @typedef {{ number: number, span: string, length: string, reverseComplement: boolean }} LinearRecordWindowRow
 */
/**
 * The rows of the read-only window list of the Circular Region block (a Linear card edits its windows,
 * app/linear-card-window-readers.js), in column order (D-08; windows without a column keep the request's
 * order after those with one, as the layout places them). A window whose start is after its end crosses the origin
 * of a circular record (D-34): its length needs the record's, and is empty until the record is read.
 * @param {readonly { start: number | null, end: number | null, reverseComplement: boolean, column: number | null }[]} windows
 * @param {number | null} recordLength the length of the record the card draws, once read
 * @returns {LinearRecordWindowRow[]}
 */
export const linearRecordWindowRows = (windows, recordLength) => windows
  .map((window, index) => ({ window, index }))
  .sort((a, b) => (a.window.column ?? Infinity) - (b.window.column ?? Infinity) || a.index - b.index)
  .map(({ window }, order) => {
    const { start, end } = window;
    const stated = start !== null && end !== null;
    const bases = stated ? windowBases(start, end, recordLength) : null;
    return {
      number: order + 1,
      span: stated ? `${numberText(start)}\u2013${numberText(end)}` : '',
      length: bases === null ? '' : `${numberText(bases)} bp`,
      reverseComplement: window.reverseComplement
    };
  });

/**
 * The region a Linear record card's fields state, for what Generate draws (the
 * analysis inputs, the match sequences, the note); the request projection
 * `linearRegionPayload` parses them itself and sends them as typed. null
 * without one (both fields empty); `{ error }`, a REGION_INVALID diagnostic
 * naming the card's 1-based `ordinal`, while Generate cannot draw it; else the
 * endpoints as typed (an End before the Start draws the reverse complement) and
 * `drawn`, the span Generate draws from a record of `recordLength`
 * (`regionWithinRecord`; null when it starts after the record ends).
 * @param {{ region_start?: unknown, region_end?: unknown }} seq
 * @param {number} ordinal
 * @param {unknown} [recordLength]
 * @returns {{ start: number, end: number, drawn: RecordRegion | null } | { error: Error } | null}
 */
export const linearCardRegion = (seq, ordinal, recordLength) => {
  /** @param {unknown} value */
  const stated = (value) => value !== null && value !== undefined && value !== '';
  /** @param {string} reason */
  const invalid = (reason) => ({ error: diagnosticError('REGION_INVALID', { inputOrdinal: ordinal, reason }) });
  const hasStart = stated(seq.region_start);
  if (hasStart !== stated(seq.region_end)) return invalid('BOTH_ENDPOINTS');
  if (!hasStart) return null;
  const start = Number(seq.region_start);
  const end = Number(seq.region_end);
  if (!Number.isFinite(start) || !Number.isFinite(end)) return invalid('FINITE');
  if (!Number.isInteger(start) || !Number.isInteger(end)) return invalid('INTEGER');
  if (start < 1 || end < 1) return invalid('POSITIVE_INTEGER');
  return { start, end, drawn: regionWithinRecord(Math.min(start, end), Math.max(start, end), recordLength) };
};

/**
 * The note beside the Region fields while Generate draws a shorter region than typed; '' otherwise.
 * @param {RecordRegion | null} region
 */
export const regionClampNote = (region) => (
  region?.recordLength === undefined ? '' : `Region reaches ${numberText(region.end)}, past the record end`
    + ` (${formatRecordLength(region.recordLength)}); Generate draws to ${numberText(region.drawnEnd)}.`
);

/**
 * The region the Circular form crops a record to, as its fields state it: null
 * without one (both fields empty), and `{ error }` while Generate cannot draw
 * it. The one reader of those fields; only a record offered alone is cropped
 * (`isSingleCircularJourney`).
 * @param {CircularRegionForm} form
 * @param {{ recordLength?: unknown, record_length?: unknown } | null} [record] the record it crops, when known
 * @returns {RecordRegion | { error: Error } | null}
 */
const circularFormRegion = (form, record = null) => {
  const hasStart = cleanText(form.circular_region_start) !== '';
  const hasEnd = cleanText(form.circular_region_end) !== '';
  if (!hasStart && !hasEnd) return null;
  /** @param {string} reason */
  const invalid = (reason) => ({ error: diagnosticError('REGION_INVALID', { field: 'region', reason }) });
  if (hasStart !== hasEnd) return invalid('BOTH_ENDPOINTS');
  const start = Number(form.circular_region_start);
  const end = Number(form.circular_region_end);
  if (!Number.isInteger(start) || !Number.isInteger(end) || start < 1 || end < 1) return invalid('POSITIVE_INTEGER');
  if (start > end) return invalid('ORDER');
  return regionWithinRecord(start, end, record?.recordLength ?? record?.record_length) ?? invalid('RECORD_BOUNDS');
};

/**
 * The crop of a Circular record as Generate sends it: its span and orientation, the bases it draws (null while the
 * length of a record whose window crosses the origin is not known), and the region of the Start/End fields it came
 * from (null for a listed window).
 * @typedef {{ start: number, end: number, reverseComplement: boolean, bases: number | null, region: RecordRegion | null }} CircularCrop
 */
/**
 * The crop of the record a Circular drawing offers alone: the window its list holds (one window across the origin,
 * D-45: the list is the mark, never the Start/End fields), else the region its fields state (`circularFormRegion`).
 * null without one, and `{ error }` while Generate cannot draw it. The one reader of the crop for the request, the
 * estimates and the notes.
 * @param {CircularRegionForm & { circular_reverse?: unknown }} form
 * @param {{ recordLength?: unknown, record_length?: unknown } | null} [record] the record it crops, when known
 * @returns {CircularCrop | { error: Error } | null}
 */
export const circularFormCrop = (form, record = null) => {
  const [listed] = form.circular_windows;
  if (listed) {
    const { start, end } = listed;
    if (start === null || end === null) {
      return { error: diagnosticError('REGION_INVALID', { field: 'region', reason: 'BOTH_ENDPOINTS' }) };
    }
    const length = Number(record?.recordLength ?? record?.record_length);
    const bases = windowBases(start, end, Number.isInteger(length) && length > 0 ? length : null);
    return { start, end, reverseComplement: listed.reverseComplement, bases, region: null };
  }
  const region = circularFormRegion(form, record);
  if (region === null || 'error' in region) return region;
  return {
    start: region.start, end: region.end, reverseComplement: Boolean(form.circular_reverse),
    bases: region.drawnEnd - region.start + 1, region
  };
};

/**
 * The crop Generate applies to a record: the form's crop when the drawing
 * offers that record alone and Generate can draw it; null otherwise
 * (a crop Generate refuses crops nothing).
 * @param {boolean} single whether the drawing offers one record alone (`isSingleCircularJourney`)
 * @param {CircularRegionForm} form
 * @param {{ recordLength?: unknown, record_length?: unknown } | null} record the record it crops, when known
 * @returns {CircularCrop | null}
 */
const circularCropRegion = (single, form, record) => {
  const crop = single ? circularFormCrop(form, record) : null;
  return crop !== null && !('error' in crop) ? crop : null;
};

/**
 * The record a Circular drawing offers alone (`isSingleCircularJourney`), which its next Generate draws, with the crop
 * Generate applies to it (`circularCropRegion`, null when it draws the whole record); null while the drawing offers
 * several records or draws none. The Topology rows read it (OV-3340).
 * @template {Partial<CircularListRecord> & { recordId?: unknown, recordLength?: unknown }} T
 * @param {readonly T[] | null | undefined} records the records discovered in the input
 * @param {CircularRecordDrawing} drawing
 */
export const circularDrawingAloneRecord = (records, drawing) => {
  const recordSet = circularDrawingRecordSet(records, drawing);
  const [record] = isSingleCircularJourney(recordSet) ? recordSet.records : [];
  return record === undefined ? null : { record, crop: circularCropRegion(true, drawing.form, record) };
};

/**
 * The records a Circular drawing draws cropped in a grid or batch, by identity: those its ranges keep
 * (`circularRecordRange`), each with the first window of its range (null for a range of no window); none while it
 * offers one record alone, whose crop is the form's (`circularDrawingAloneRecord`).
 * @param {CircularRequestRecordSetOptions['records']} records the records discovered in the input
 * @param {CircularRecordDrawing} drawing
 * @returns {Map<string, LinearSeqWindow | null>}
 */
export const circularDrawingCroppedRecords = (records, drawing) => {
  const recordSet = circularDrawingRecordSet(records, drawing);
  /** @type {Map<string, LinearSeqWindow | null>} */
  const cropped = new Map();
  if (isSingleCircularJourney(recordSet) || !drawing.form.circular_record_ranges.length) return cropped;
  recordSet.recordSelectors.forEach((entry) => {
    const range = circularRecordRange(drawing.form.circular_record_ranges, recordSet.recordSelectors, entry);
    if (range !== null) cropped.set(entry.value, range.windows[0] ?? null);
  });
  return cropped;
};

/**
 * The note beside the Circular Region fields (`regionClampNote`) for the
 * region the next Generate crops its record to.
 * @param {CircularRequestRecordSetOptions['records']} records the records discovered in the input
 * @param {CircularRecordDrawing} drawing
 */
export const circularDrawingRegionNote = (records, drawing) => {
  const recordSet = circularDrawingRecordSet(records, drawing);
  return regionClampNote(circularCropRegion(isSingleCircularJourney(recordSet), drawing.form, recordSet.records[0] ?? null)?.region ?? null);
};

/**
 * The read-only list of the window a Circular drawing's record keeps across the origin (`linearRecordWindowRows`),
 * with the length of the record it offers alone, once read; empty without one.
 * @param {CircularRequestRecordSetOptions['records']} records the records discovered in the input
 * @param {CircularRecordDrawing} drawing
 */
export const circularDrawingWindowRows = (records, drawing) => {
  if (!drawing.form.circular_windows.length) return [];
  const length = Number(circularDrawingRecordSet(records, drawing).records[0]?.record_length);
  return linearRecordWindowRows(drawing.form.circular_windows, Number.isInteger(length) && length > 0 ? length : null);
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
  /** @param {{ record_length?: unknown, identity?: unknown, uid?: unknown, value: string, selector: string, recordId: string }} record */
  const drawnLength = (record) => {
    const crop = circularCropRegion(single, drawing.form, record);
    if (crop !== null) return crop.bases ?? 0;
    const length = Number(record.record_length);
    const range = single ? null : circularRecordRange(drawing.form.circular_record_ranges, recordSet.recordSelectors, record);
    return range === null ? length : range.windows.reduce((sum, window) => (
      sum + (windowBases(Number(window.start), Number(window.end), Number.isInteger(length) && length > 0 ? length : null) ?? 0)
    ), 0);
  };
  return (drawing.form.multi_record_canvas ? [drawn] : drawn.map((record) => [record]))
    .map((diagram) => ({
      recordKeys: diagram.map((record) => cleanText(record.identity)),
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
