// @ts-check
// Sequence search (PD-OI-116): the data the panel sends to Python and reads
// back. Python reads the queries (`readSearchQueries`, the FASTA-only loader)
// and turns one search's LOSAT rows into an annotation set
// (`searchHitAnnotationSet`, gbdraw/web_support/sequence_search.py); this
// module types both boundaries, builds the searched sequences from the records
// the shown drawing draws, and names the set a search adds. It reads no
// sequence and decides no hit itself.
import { diagnosticError } from '../utils/error-normalization.js';
import { SEARCH_QUERY_COLORS } from '../utils/sequence-search-colors.generated.js';
import { sourceSpanOf } from './region-axis.js';

/** @import { RegionAxisRecord } from './region-axis.js' */
/** @import { RecordSpanInput, RecordSpanSequencesReply } from './diagram-worker-protocol.js' */
/** @import { AnnotationRecordSelector } from './annotation-state.js' */

/**
 * One query as Python read it (`SearchQueryRead`): its FASTA ID, its length in bases, and the FASTA text LOSAT reads.
 * @typedef {object} SearchQueryRead
 * @property {string} name
 * @property {number} length
 * @property {string} fasta
 */
/** A query of the panel: as read, with the hit colour the panel shows and the user may change. @typedef {SearchQueryRead & { color: string }} SearchQuery */

/**
 * One searched sequence (`SearchSubjectInput`): the ID LOSAT names it by (its FASTA header), the drawing's record
 * selector, the request's record key, and the source coordinate of its first base.
 * @typedef {object} SearchSubject
 * @property {string} subjectId
 * @property {AnnotationRecordSelector} record
 * @property {string} recordKey
 * @property {number} sourceStart
 */

/** @typedef {'query' | 'gradient'} SearchColorMode */
/** @typedef {{ evalue: number, bitscore: number, identity: number, alignmentLength: number }} SearchThresholds */

/**
 * The payload of `searchHitAnnotationSet` (`SequenceSearchInput`); every field is required.
 * @typedef {object} SequenceSearchInput
 * @property {string} table The raw outfmt 6 text LOSAT wrote.
 * @property {string} setId
 * @property {SearchSubject[]} subjects
 * @property {SearchThresholds} thresholds
 * @property {'blastn'} program
 * @property {{ name: string, length: number }[]} queries In panel order.
 * @property {SearchColorMode} colorMode
 * @property {Record<string, string>} queryColors
 */

/**
 * One hit of the reply (`SearchHitReply`).
 * @typedef {object} SearchHit
 * @property {string} annotationId
 * @property {string} query
 * @property {string} recordKey
 * @property {number} start
 * @property {number} end
 * @property {'+' | '-'} strand
 * @property {number} identity
 * @property {number} alignmentLength
 * @property {number} evalue
 * @property {number} bitscore
 * @property {number} qstart
 * @property {number} qend
 */
/** @typedef {{ name: string, threshold: number, rowsLeftOut: number, bestValue: number }} ThresholdTally */
/**
 * The `searchHitAnnotationSet` reply (`SequenceSearchReply`): the set as the canonical request writes it (null when no
 * hit passed the thresholds), its hits, and what the thresholds and IDs left out.
 * @typedef {object} SequenceSearchReply
 * @property {Record<string, unknown> | null} annotationSet
 * @property {SearchHit[]} hits
 * @property {number} rowsRead
 * @property {number} rowsKept
 * @property {ThresholdTally[]} thresholds
 * @property {number} unmatchedRows
 * @property {string[]} unmatchedSubjectIds
 */

/** @param {unknown} value @returns {value is Record<string, unknown>} */
const isRecord = (value) => Boolean(value) && typeof value === 'object' && !Array.isArray(value);
/** @param {Record<string, unknown>} value @param {readonly string[]} keys */
const strings = (value, keys) => keys.every((key) => typeof value[key] === 'string');
/** @param {Record<string, unknown>} value @param {readonly string[]} keys */
const numbers = (value, keys) => keys.every((key) => typeof value[key] === 'number' && Number.isFinite(value[key]));
/** @param {string} operation */
const invalid = (operation) => diagnosticError('RESULT_INVALID', {}, { stage: 'helper', operation });

/**
 * The `readSearchQueries` result, checked.
 * @param {unknown} result
 * @returns {SearchQueryRead[]}
 */
export const readSearchQueriesReply = (result) => {
  if (!isRecord(result) || !Array.isArray(result.queries)
    || !result.queries.every((query) => isRecord(query) && strings(query, ['name', 'fasta']) && Number.isSafeInteger(query.length))) {
    throw invalid('readSearchQueries');
  }
  return /** @type {SearchQueryRead[]} */ (result.queries);
};

/** @param {unknown} hit @returns {hit is SearchHit} */
const isHit = (hit) => isRecord(hit) && strings(hit, ['annotationId', 'query', 'recordKey'])
  && numbers(hit, ['start', 'end', 'identity', 'alignmentLength', 'evalue', 'bitscore', 'qstart', 'qend'])
  && (hit.strand === '+' || hit.strand === '-');

/**
 * The `searchHitAnnotationSet` result, checked.
 * @param {unknown} result
 * @returns {SequenceSearchReply}
 */
export const readSequenceSearchReply = (result) => {
  if (!isRecord(result) || !(result.annotationSet === null || isRecord(result.annotationSet))
    || !Array.isArray(result.hits) || !result.hits.every(isHit) || !numbers(result, ['rowsRead', 'rowsKept', 'unmatchedRows'])
    || !Array.isArray(result.thresholds) || !result.thresholds.every((item) => isRecord(item) && typeof item.name === 'string')
    || !Array.isArray(result.unmatchedSubjectIds)) {
    throw invalid('searchHitAnnotationSet');
  }
  return /** @type {SequenceSearchReply} */ (/** @type {unknown} */ (result));
};

/**
 * The panel's queries after a read: the queries as read, each with the colour of its place (or the colour it had
 * under the same name, so a reread keeps a changed colour).
 * @param {readonly SearchQueryRead[]} read
 * @param {readonly SearchQuery[]} [previous]
 * @returns {SearchQuery[]}
 */
export const panelQueries = (read, previous = []) => read.map((query, index) => ({
  ...query,
  color: previous.find((item) => item.name === query.name)?.color ?? SEARCH_QUERY_COLORS[index % SEARCH_QUERY_COLORS.length]
}));

/**
 * The spans a search reads: each record the shown drawing draws (every Result's), once, over the bases it draws, on
 * its source strand (a whole record from 1).
 * @param {readonly RegionAxisRecord[]} records
 * @returns {RecordSpanInput[]}
 */
export const searchSpans = (records) => {
  /** @type {Map<string, RecordSpanInput>} */
  const spans = new Map();
  for (const { recordKey, axis } of records) {
    if (spans.has(recordKey)) continue;
    const { start, end } = sourceSpanOf(axis.source, { start: 0, length: axis.source.drawnLength });
    spans.set(recordKey, { recordKey, start, end });
  }
  return [...spans.values()];
};

/**
 * The searched sequences of the spans Python read: LOSAT names each by its FASTA header (`recordId:start-end`), and
 * an annotation names its record by ID, or by its place when two drawn records share that ID.
 * @param {RecordSpanSequencesReply} read
 * @param {readonly RegionAxisRecord[]} records
 * @returns {SearchSubject[]}
 */
export const searchSubjects = (read, records) => read.spans.map((span) => {
  const shared = read.spans.filter((other) => other.recordId === span.recordId).length > 1;
  const index = records.find((record) => record.recordKey === span.recordKey)?.recordIndex ?? 0;
  return {
    subjectId: span.header,
    record: shared || !span.recordId ? { kind: 'recordIndex', index } : { kind: 'recordId', value: span.recordId },
    recordKey: span.recordKey,
    sourceStart: span.start
  };
});

/**
 * The id of the set the next search adds: "Search N", N the first number no set of the drawing has.
 * @param {readonly { id: string }[]} sets
 */
export const nextSearchSetId = (sets) => {
  const ids = new Set(sets.map(({ id }) => id));
  let n = 1;
  while (ids.has(`Search ${n}`)) n += 1;
  return `Search ${n}`;
};

/**
 * One changed value of a History step (`HistoryIntentChange`), as far as this module reads it.
 * @typedef {{ readonly path?: readonly (string | number)[], readonly before?: unknown, readonly after?: unknown }} SearchStepChange
 */
/** The History intent path of a drawing's annotation sets. */
const ANNOTATION_SETS_PATH = 'config.annotationSets';
/** @param {unknown} sets @returns {Map<string, string>} each set's id and its JSON */
const setsById = (sets) => new Map((Array.isArray(sets) ? sets : [])
  .filter((set) => isRecord(set) && typeof set.id === 'string')
  .map((set) => [/** @type {string} */ (set.id), JSON.stringify(set)]));

/**
 * The sets a History step added that a search added (the search's own step, its Undo, or its Redo): the ids the step
 * adds to the drawing's annotation sets whose set is the one the search wrote (`searched`, by id, the set's JSON).
 * The step's Undo removes them; a step that removed or renamed a set, or added a set no search wrote, names none.
 * @param {readonly SearchStepChange[] | null | undefined} changes
 * @param {ReadonlyMap<string, string>} searched
 * @returns {string[]}
 */
export const searchSetIdsOfStep = (changes, searched) => (Array.isArray(changes) ? changes : [])
  .filter((change) => (change.path ?? []).join('.') === ANNOTATION_SETS_PATH)
  .flatMap((change) => {
    const before = setsById(change.before);
    return [...setsById(change.after)].filter(([id, json]) => !before.has(id) && searched.get(id) === json).map(([id]) => id);
  });
