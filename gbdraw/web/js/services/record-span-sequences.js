// @ts-check
// The bases of spans of the drawn records, read by Python (`readRecordSpanSequences`,
// gbdraw/web_support/record_spans.py) from the request the shown Result was made from: Copy FASTA of Select region.
// The Web slices no sequence; this module sends the request and reads the reply.
/** @import { RecordSpanRead, RecordSpanSequence, RecordSpanSequencesPayload, RecordSpanSequencesReply } from './diagram-worker-protocol.js' */
/** @import { CommittedCanonicalSession } from './session-request.js' */
import { diagnosticError } from '../utils/error-normalization.js';

const FAILURE = /** @type {const} */ ({ stage: 'helper', operation: 'readRecordSpanSequences' });

/** @param {unknown} value @returns {value is Record<string, unknown>} */
const isRecord = (value) => Boolean(value) && typeof value === 'object' && !Array.isArray(value);

/** @param {unknown} value @returns {value is RecordSpanSequence} */
const isSpanSequence = (value) => isRecord(value)
  && ['recordKey', 'recordId', 'header', 'sequence'].every((key) => typeof value[key] === 'string')
  && ['start', 'end', 'length', 'recordLength'].every((key) => Number.isSafeInteger(value[key]))
  && typeof value.wrapsOrigin === 'boolean' && typeof value.originWindow === 'boolean';

/**
 * The helper's result, checked against `RecordSpanSequencesReply`: it reads every record asked and no other (a drawn
 * record may give two spans, its drawn bases and its origin window).
 * @param {unknown} result
 * @param {readonly RecordSpanRead[]} spans The spans asked.
 * @returns {RecordSpanSequencesReply}
 */
export const readRecordSpanSequencesReply = (result, spans) => {
  const asked = new Set(spans.map(({ recordKey }) => recordKey));
  if (!isRecord(result) || typeof result.fasta !== 'string' || !Array.isArray(result.spans)
    || !result.spans.every(isSpanSequence)
    || new Set(result.spans.map(({ recordKey }) => recordKey)).size !== asked.size
    || result.spans.some(({ recordKey }) => !asked.has(recordKey))) {
    throw diagnosticError('RESULT_INVALID', {}, FAILURE);
  }
  return /** @type {RecordSpanSequencesReply} */ (/** @type {unknown} */ (result));
};

/**
 * Reads spans of the records of the committed request, the one the shown Result was made from.
 * @param {CommittedCanonicalSession | null} committed
 * @param {readonly RecordSpanRead[]} spans
 * @param {(payload: RecordSpanSequencesPayload) => Promise<{ result?: unknown }>} runHelper Runs the Worker helper.
 * @returns {Promise<RecordSpanSequencesReply>}
 */
export const requestRecordSpanSequences = async (committed, spans, runHelper) => {
  // A Result with no committed request (none generated or loaded) has no input to read.
  if (!committed) throw diagnosticError('RESOURCE_INVALID', {}, FAILURE);
  const reply = await runHelper({ canonicalRequest: committed.renderRequest, resources: committed.resources, spans });
  return readRecordSpanSequencesReply(reply?.result, spans);
};
