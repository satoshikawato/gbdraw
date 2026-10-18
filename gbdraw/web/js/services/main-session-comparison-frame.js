// @ts-check
// origin/main Web Sessions (version 42 and older) stored the nucleotide
// comparison rows of a reverse-complemented Linear record after the reverse
// complement. Comparison tables now use the search frame (D-18, PD-OI-073),
// so Load rewrites those stored tables once; a re-saved Session keeps the
// search frame. The table is read and rewritten by Python (one table reader).
// Sessions up to 47 stored the protein results of such a record in the drawn
// frame; Python converts them once the stored figure is shown (OV-1145, D-128).
import { discoverGffFastaRecords, discoverSequenceRecords } from './record-discovery.js';
import { base64ToBytes, bytesToText, sha256Hex, textToBase64, textToBytes } from './byte-utils.js';
import { isResourceBackedCanonicalComparison } from './canonical-comparisons.js';
import { diagnosticError } from '../utils/error-normalization.js';
import {
  DRAWN_FRAME_PROTEIN_RESULT_SESSION_VERSION, MAIN_DISPLAY_FRAME_SESSION_VERSION
} from '../contracts/schema.generated.js';

/** @import { SequenceRecord } from './record-discovery.js' */
/** @import { ResourceBackedComparison } from './canonical-comparisons.js' */
/** @import { RawSessionResourceDescriptor } from './session-resource-backing.js' */
/** @import { ResourceDescriptor } from './diagram-resource-staging.js' */
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
/** @typedef {Readonly<Record<string, Partial<ResourceDescriptor> | null | undefined>>} DrawingResources */
/**
 * One drawing's committed request and the resources it names, as the Worker
 * stages them (`convertStoredProteinResultFrames`).
 * @typedef {{ version: number, canonicalRequest: unknown, resources: DrawingResources }} ProteinResultConversion
 */
/** @typedef {{ renderRequest?: unknown, resources?: DrawingResources }} StoredDrawing A committed drawing. */
/** @typedef {{ data: string, size: number }} ConvertedResource The base64 bytes and size Python returned. */
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
// And the Worker protein result conversion (OV-1145).
/** @type {((request: ProteinResultConversion) => Promise<unknown>) | null} */
let proteinResultConverter = null;
/** @param {typeof proteinResultConverter} converter */
export const setStoredProteinResultFrameConverter = (converter) => {
  proteinResultConverter = typeof converter === 'function' ? converter : null;
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
 * @returns {StoredResource | null} null for a file the Session does not hold
 */
const heldResource = (resources, resourceId) => {
  const entry = /** @type {StoredResource} */ (resources[String(resourceId)]);
  return entry && typeof entry.data === 'string' ? entry : null;
};

/**
 * @param {Record<string, unknown>} resources
 * @param {unknown} resourceId
 * @returns {File | null} null for a file the Session does not hold
 */
const resourceFile = (resources, resourceId) => {
  const entry = heldResource(resources, resourceId);
  if (!entry) return null;
  return new File([base64ToBytes(String(entry.data))], String(entry.name || resourceId), {
    type: String(entry.type || 'text/plain')
  });
};

/**
 * @param {StoredRecord} record
 * @returns {{ primary: unknown, fasta: unknown } | { primary: unknown, fasta?: undefined }} the resource IDs of its
 *   source files: a GFF3 source names its FASTA file too
 */
const recordResourceIds = (record) => {
  const source = /** @type {StoredSource} */ (record?.source);
  return source?.kind === 'gffFasta'
    ? { primary: source.gffResourceId, fasta: source.fastaResourceId }
    : { primary: source?.resourceId };
};

/**
 * @param {StoredRecord} record
 * @param {Record<string, unknown>} resources
 * @param {RecordDiscovery} discover
 * @returns {Promise<number | null>} null when the Session does not hold the record's files
 */
const selectedRecordLength = async (record, resources, discover) => {
  const ids = recordResourceIds(record);
  const file = resourceFile(resources, ids.primary);
  const fastaFile = ids.fasta === undefined ? null : resourceFile(resources, ids.fasta);
  if (!file || (ids.fasta !== undefined && !fastaFile)) return null;
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
 * Whether a stored Linear request holds a protein result that touches one of
 * its reverse-complemented records. A precomputed table touches its two
 * records; a typed result the records its `generatedProteinComparison`
 * marker's pairs name, or every record when no marker names one. A reversed
 * record whose files the Session does not hold is the request reader's to
 * drop (OV-869). Indices are read as property keys, as `entryAt` reads them.
 * @param {unknown} stored
 * @param {Record<string, unknown>} resources
 * @returns {boolean}
 */
const storesReversedProteinResult = (stored, resources) => {
  const request = /** @type {UnvalidatedFields<keyof CanonicalLinearRenderRequest>} */ (stored);
  if (request?.mode !== 'linear' || !Array.isArray(request.comparisons) || !Array.isArray(request.records)) return false;
  /** @type {readonly unknown[]} */
  const comparisons = request.comparisons;
  /** @type {readonly unknown[]} */
  const records = request.records;
  const reversed = new Set(records.flatMap((item, position) => {
    const record = /** @type {StoredRecord} */ (item);
    const ids = recordResourceIds(record);
    const held = heldResource(resources, ids.primary) && (ids.fasta === undefined || heldResource(resources, ids.fasta));
    return storedRecordReverses(record) && held ? [String(position)] : [];
  }));
  /** @param {unknown} item */
  const touches = (item) => {
    const pair = /** @type {UnvalidatedFields<'queryRecordIndex' | 'subjectRecordIndex'>} */ (item);
    return reversed.has(String(pair?.queryRecordIndex)) || reversed.has(String(pair?.subjectRecordIndex));
  };
  const named = comparisons.flatMap((item) => {
    const marker = /** @type {UnvalidatedFields<'kind' | 'pairs'>} */ (item);
    return marker?.kind === 'generatedProteinComparison' && Array.isArray(marker.pairs)
      ? /** @type {readonly unknown[]} */ (marker.pairs) : [];
  });
  return reversed.size > 0 && comparisons.some((item) => {
    const comparison = /** @type {ResourceBackedComparison} */ (item);
    if (!isResourceBackedCanonicalComparison(comparison)) return false;
    return comparison.kind === 'precomputedProteinComparison' ? touches(comparison) : named.length === 0 || named.some(touches);
  });
};

/**
 * Whether a Session up to 47 stores protein results of a reverse-complemented
 * Linear record, which it stored in the drawn frame; Session 48 stores them in
 * the search frame and Python converts older ones (OV-1145). Both committed
 * requests are read: the shown one and the other mode's.
 * @param {{ version?: unknown, renderRequest?: unknown, otherModeResult?: unknown, resources?: unknown }} data
 * @returns {boolean}
 */
export const storesDrawnFrameProteinResults = (data) => {
  if (!(Number(data?.version) <= DRAWN_FRAME_PROTEIN_RESULT_SESSION_VERSION)) return false;
  const other = /** @type {UnvalidatedFields<'renderRequest'>} */ (data?.otherModeResult);
  const resources = /** @type {Record<string, unknown>} */ (data?.resources ?? {});
  return [data?.renderRequest, other?.renderRequest].some((stored) => storesReversedProteinResult(stored, resources));
};

const PROTEIN_RESULT_CONVERSION = /** @type {const} */ ({ stage: 'helper', operation: 'convertStoredProteinResultFrames' });

/**
 * B-LOAD (D-128): the protein results that each committed drawing of a
 * Session up to 47 stored for a reversed record, in the search frame. Python
 * converts each drawing whose results touch one (one converter, no JS copy);
 * the bytes come back by stored resource ID. A failure throws, so the Load
 * changes nothing (PD-OI-046).
 * @param {number} version
 * @param {ReadonlyArray<StoredDrawing | null | undefined>} drawings
 * @param {{ convert?: typeof proteinResultConverter }} [options]
 * @returns {Promise<Map<string, ConvertedResource>>} empty when no drawing needs it
 */
export const convertStoredProteinResultFrames = async (version, drawings, { convert = proteinResultConverter } = {}) => {
  /** @type {Map<string, ConvertedResource>} */
  const converted = new Map();
  if (!(version <= DRAWN_FRAME_PROTEIN_RESULT_SESSION_VERSION)) return converted;
  for (const drawing of drawings) {
    const resources = drawing?.resources ?? {};
    if (!drawing || !storesReversedProteinResult(drawing.renderRequest, resources)) continue;
    if (!convert) throw diagnosticError('HELPER_PROTOCOL', {}, PROTEIN_RESULT_CONVERSION);
    const response = /** @type {UnvalidatedFields<'result', UnvalidatedFields<'resources'>>} */ (await convert(
      { version, canonicalRequest: drawing.renderRequest, resources }
    ));
    const returned = response?.result?.resources;
    if (!returned || typeof returned !== 'object' || Array.isArray(returned)) {
      throw diagnosticError('RESULT_INVALID', {}, PROTEIN_RESULT_CONVERSION);
    }
    for (const [resourceId, value] of Object.entries(returned)) {
      const item = /** @type {UnvalidatedFields<keyof ConvertedResource>} */ (value);
      const previous = converted.get(resourceId);
      if (typeof item?.data !== 'string' || !Number.isSafeInteger(item.size) || Number(item.size) < 0
        || (previous && previous.data !== item.data)) {
        throw diagnosticError('RESULT_INVALID', {}, PROTEIN_RESULT_CONVERSION);
      }
      converted.set(resourceId, { data: item.data, size: Number(item.size) });
    }
  }
  return converted;
};

/**
 * `resources` with the converted bytes under their stored IDs: each keeps its
 * other fields, takes the new size, and a checksum it held is computed again.
 * @param {Record<string, unknown>} resources
 * @param {ReadonlyMap<string, ConvertedResource>} converted
 * @returns {Promise<Record<string, unknown>>} `resources` itself when it holds none of them
 */
export const withConvertedResources = async (resources, converted) => {
  const ids = Object.keys(resources).filter((resourceId) => converted.has(resourceId));
  if (ids.length === 0) return resources;
  const next = { ...resources };
  for (const resourceId of ids) {
    const { checksum, ...descriptor } = /** @type {RawSessionResourceDescriptor} */ (resources[resourceId]);
    const { data, size } = /** @type {ConvertedResource} */ (converted.get(resourceId));
    next[resourceId] = { ...descriptor, encoding: 'base64', data, size,
      ...(checksum ? { checksum: await sha256Hex(base64ToBytes(data)) } : {}) };
  }
  return next;
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
