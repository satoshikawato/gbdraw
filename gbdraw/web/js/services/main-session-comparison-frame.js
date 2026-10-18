// @ts-check
// origin/main Web Sessions (version 42 and older) stored the nucleotide
// comparison rows of a reverse-complemented Linear record after the reverse
// complement. Comparison tables now use the search frame (D-18, PD-OI-073),
// so Load rewrites those stored tables once; a re-saved Session keeps the
// search frame. The table is read and rewritten by Python (one table reader).
import { discoverGffFastaRecords, discoverSequenceRecords } from './record-discovery.js';
import { base64ToBytes, bytesToText, textToBase64, textToBytes } from './byte-utils.js';
import { isResourceBackedCanonicalComparison } from './canonical-comparisons.js';
import {
  DRAWN_FRAME_PROTEIN_RESULT_SESSION_VERSION, MAIN_DISPLAY_FRAME_SESSION_VERSION
} from '../contracts/schema.generated.js';

/** @import { SequenceRecord } from './record-discovery.js' */
/** @import { ResourceBackedComparison } from './canonical-comparisons.js' */
/** @import { RawSessionResourceDescriptor } from './session-resource-backing.js' */
/** @import { UnvalidatedFields } from './session-authority.js' */
/** @import { SessionChecked31, SessionChecked39, SessionChecked40, SessionChecked46 } from '../contracts/session.generated.js' */
/**
 * @import { CanonicalGenbankRecordSource, CanonicalGffFastaRecordSource, CanonicalLinearRenderRequest,
 *   CanonicalNucleotideBlastComparison, CanonicalRecord, CanonicalRecordIdSelector, CanonicalRecordIndexSelector,
 *   CanonicalRecordPresentation, CanonicalRecordRegion } from '../contracts/request.generated.js'
 */

// The conversion reads the stored request before `readSessionRequests` reads
// it, so each part is read through a view whose every field is unvalidated,
// named as the canonical request names it. A stored selector of an older
// request schema may name its record by `id`.
/** @typedef {UnvalidatedFields<keyof CanonicalRecord>} StoredRecord */
/** @typedef {UnvalidatedFields<keyof CanonicalGenbankRecordSource | keyof CanonicalGffFastaRecordSource>} StoredSource */
/** @typedef {UnvalidatedFields<keyof CanonicalRecordIdSelector | keyof CanonicalRecordIndexSelector | 'id'>} StoredSelector */
/** @typedef {UnvalidatedFields<keyof RawSessionResourceDescriptor>} StoredResource */
/** @typedef {{ length: number, reverse: boolean }} EndpointFrame */
/** @typedef {{ tableText: string, queryFrame: EndpointFrame, subjectFrame: EndpointFrame }} TableConversion */
/**
 * The record discovery the conversion reads a record length from.
 * @typedef {{
 *   sequence: (options: { file: File, format: string }) => Promise<SequenceRecord[]>,
 *   gffFasta: (options: { gffFile: File, fastaFile: File }) => Promise<SequenceRecord[]>
 * }} RecordDiscovery
 */

// The app injects the Worker table conversion, like the configuration validator.
/** @type {((request: TableConversion) => Promise<unknown>) | null} */
let tableConverter = null;
/** @param {typeof tableConverter} converter */
export const setMainSessionComparisonFrameConverter = (converter) => {
  tableConverter = typeof converter === 'function' ? converter : null;
};

/**
 * The entry of `list` a stored index names: the index is read as a property
 * key, as `list[index]` reads it, so `"1"` names entry 1 and `null` none.
 * @template E
 * @param {readonly E[]} list
 * @param {unknown} index
 * @returns {E | undefined}
 */
const entryAt = (list, index) => list.find((_, position) => String(position) === String(index));

/**
 * @param {Record<string, unknown>} resources `readSessionDocument` kept each descriptor it could read.
 * @param {unknown} resourceId
 * @returns {File | null} null for a file the Session does not hold
 */
const resourceFile = (resources, resourceId) => {
  const entry = /** @type {StoredResource} */ (resources[String(resourceId)]);
  if (!entry || typeof entry.data !== 'string') return null;
  return new File([base64ToBytes(entry.data)], String(entry.name || resourceId), {
    type: String(entry.type || 'text/plain')
  });
};

/**
 * @param {StoredRecord} record
 * @param {Record<string, unknown>} resources
 * @param {RecordDiscovery} discover
 * @returns {Promise<number | null>} null when the Session does not hold the record's files
 */
const selectedRecordLength = async (record, resources, discover) => {
  const source = /** @type {StoredSource} */ (record?.source);
  const gffFasta = source?.kind === 'gffFasta';
  const file = resourceFile(resources, gffFasta ? source.gffResourceId : source?.resourceId);
  const fastaFile = gffFasta ? resourceFile(resources, source.fastaResourceId) : null;
  if (!file || (gffFasta && !fastaFile)) return null;
  const catalog = fastaFile
    ? await discover.gffFasta({ gffFile: file, fastaFile })
    : await discover.sequence({ file, format: 'genbank' });
  const selector = /** @type {StoredSelector} */ (record?.selector || null);
  const selected = selector?.kind === 'recordIndex'
    ? entryAt(catalog, selector.index)
    : selector?.kind === 'recordId'
      ? catalog.find((entry) => entry.recordId === String(selector.value ?? selector.id ?? ''))
      : catalog.length === 1 ? catalog[0] : null;
  const length = Number(selected?.recordLength);
  if (!Number.isSafeInteger(length) || length < 1) {
    throw new Error('A reversed Linear record of this Session has no single record length to convert its comparison table.');
  }
  return length;
};

/**
 * @param {unknown} stored A record of the stored request, or nothing.
 * @param {Record<string, unknown>} resources
 * @param {RecordDiscovery} discover
 * @returns {Promise<EndpointFrame | null>} null when the length of a reversed record cannot be read: the
 *   Session does not hold its files
 */
const endpointFrame = async (stored, resources, discover) => {
  const record = /** @type {StoredRecord} */ (stored);
  const region = /** @type {UnvalidatedFields<keyof CanonicalRecordRegion>} */ (record?.region);
  if (!storedRecordReverses(record)) return { length: 0, reverse: false };
  const length = region
    ? Number(region.end) - Number(region.start) + 1
    : await selectedRecordLength(record, resources, discover);
  return length === null ? null : { length, reverse: true };
};

/**
 * @param {StoredRecord} record
 * @returns {boolean} whether the stored record is drawn reverse-complemented
 */
const storedRecordReverses = (record) => {
  const region = /** @type {UnvalidatedFields<keyof CanonicalRecordRegion>} */ (record?.region);
  const presentation = /** @type {UnvalidatedFields<keyof CanonicalRecordPresentation>} */ (record?.presentation);
  return Boolean(region?.reverseComplement || presentation?.reverseComplement);
};

/**
 * Whether a Session up to 47 stores protein results of a reverse-complemented
 * Linear record, which it stored in the drawn frame; Session 48 stores them in
 * the search frame and Python converts older ones (OV-1145). Both committed
 * requests are read: the shown one and the other mode's.
 * @param {{ version?: unknown, renderRequest?: unknown, otherModeResult?: unknown }} data
 * @returns {boolean}
 */
export const storesDrawnFrameProteinResults = (data) => {
  if (!(Number(data?.version) <= DRAWN_FRAME_PROTEIN_RESULT_SESSION_VERSION)) return false;
  const other = /** @type {UnvalidatedFields<'renderRequest'>} */ (data?.otherModeResult);
  return [data?.renderRequest, other?.renderRequest].some((stored) => {
    const request = /** @type {UnvalidatedFields<keyof CanonicalLinearRenderRequest>} */ (stored);
    if (request?.mode !== 'linear' || !Array.isArray(request.comparisons) || !Array.isArray(request.records)) return false;
    /** @type {readonly unknown[]} */
    const comparisons = request.comparisons;
    /** @type {readonly unknown[]} */
    const records = request.records;
    return comparisons.some((item) => isResourceBackedCanonicalComparison(
      /** @type {ResourceBackedComparison} */ (item)
    )) && records.some((record) => storedRecordReverses(/** @type {StoredRecord} */ (record)));
  });
};

/**
 * Return the Session with each stored nucleotide table of a reversed endpoint
 * in the search frame. Other Sessions are returned unchanged. A Session
 * before 31 holds no request (OV-862), so the caller converts only a
 * Session of a canonical band (`isCanonicalBand`).
 * @template {SessionChecked31 | SessionChecked39 | SessionChecked40 | SessionChecked46} T
 * @param {T} data The document `readSessionDocument` read.
 * @param {{ convertTable?: typeof tableConverter, discover?: RecordDiscovery }} [options]
 * @returns {Promise<T>}
 */
export const convertMainSessionComparisonFrames = async (
  data,
  {
    convertTable = tableConverter,
    discover = { sequence: discoverSequenceRecords, gffFasta: discoverGffFastaRecords }
  } = {}
) => {
  if (data.version > MAIN_DISPLAY_FRAME_SESSION_VERSION) return data;
  const request = /** @type {UnvalidatedFields<keyof CanonicalLinearRenderRequest>} */ (data.renderRequest);
  if (request?.mode !== 'linear' || !Array.isArray(request.comparisons)) return data;
  /** @type {readonly unknown[]} */
  const comparisons = request.comparisons;
  /** @type {readonly unknown[]} */
  const records = Array.isArray(request.records) ? request.records : [];
  let resources = data.resources;
  const converted = new Set();
  for (const stored of comparisons) {
    const comparison = /** @type {UnvalidatedFields<keyof CanonicalNucleotideBlastComparison>} */ (stored);
    const resourceId = comparison?.resourceId;
    if (comparison?.kind !== 'nucleotideBlast' || !resourceId || converted.has(resourceId)) continue;
    // A table the Session does not hold is the request reader's to drop, with its row.
    const table = /** @type {StoredResource} */ (resources[String(resourceId)]);
    if (!table) continue;
    const queryFrame = await endpointFrame(entryAt(records, comparison.queryRecordIndex), resources, discover);
    const subjectFrame = await endpointFrame(entryAt(records, comparison.subjectRecordIndex), resources, discover);
    // A record file the Session does not hold, like a table, is the request reader's to drop (OV-869).
    if (!queryFrame || !subjectFrame || (!queryFrame.reverse && !subjectFrame.reverse)) continue;
    if (!convertTable) throw new Error('Session comparison conversion service is unavailable.');
    const response = /** @type {UnvalidatedFields<'result', UnvalidatedFields<'tsv'>>} */ (await convertTable(
      { tableText: bytesToText(base64ToBytes(String(table.data)), { fatal: true }), queryFrame, subjectFrame }
    ));
    const tsv = String(response?.result?.tsv ?? '');
    resources = {
      ...resources,
      [String(resourceId)]: { ...table, encoding: 'base64', data: textToBase64(tsv), size: textToBytes(tsv).length }
    };
    converted.add(resourceId);
  }
  return converted.size > 0 ? /** @type {T} */ ({ ...data, resources }) : data;
};
