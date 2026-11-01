// @ts-check
// origin/main Web Sessions (version 42 and older) stored the nucleotide
// comparison rows of a reverse-complemented Linear record after the reverse
// complement. Comparison tables now use the search frame (D-18, PD-OI-073),
// and Python rewrites those stored tables (one table reader); a re-saved
// Session keeps the search frame.
// Sessions up to 47 stored the protein results of such a record in the drawn
// frame; Python converts them once the stored figure is shown (OV-1145, D-128).
import { base64ToBytes, sha256Hex } from './byte-utils.js';
import { isResourceBackedCanonicalComparison } from './canonical-comparisons.js';
import { diagnosticError } from '../utils/error-normalization.js';
import { DRAWN_FRAME_PROTEIN_RESULT_SESSION_VERSION } from '../contracts/schema.generated.js';

/** @import { ResourceBackedComparison } from './canonical-comparisons.js' */
/** @import { RawSessionResourceDescriptor } from './session-resource-backing.js' */
/** @import { ResourceDescriptor } from './diagram-resource-staging.js' */
/** @import { UnvalidatedFields } from './session-authority.js' */
/**
 * @import { CanonicalGenbankRecordSource, CanonicalGffFastaRecordSource, CanonicalLinearRenderRequest,
 *   CanonicalRecord, CanonicalRecordPresentation, CanonicalRecordRegion } from '../contracts/request.generated.js'
 */

// The conversion reads the stored request before `readSessionRequests` reads
// it, so each part is read through a view whose every field is unvalidated,
// named as the canonical request names it.
/** @typedef {UnvalidatedFields<keyof CanonicalRecord>} StoredRecord */
/** @typedef {UnvalidatedFields<keyof CanonicalGenbankRecordSource | keyof CanonicalGffFastaRecordSource>} StoredSource */
/** @typedef {UnvalidatedFields<keyof RawSessionResourceDescriptor>} StoredResource */
/** @typedef {Readonly<Record<string, Partial<ResourceDescriptor> | null | undefined>>} DrawingResources */
/**
 * One drawing's committed request and the resources it names, as the Worker
 * stages them (`convertStoredProteinResultFrames`).
 * @typedef {{ version: number, canonicalRequest: unknown, resources: DrawingResources }} ProteinResultConversion
 */
/** @typedef {{ renderRequest?: unknown, resources?: DrawingResources }} StoredDrawing A committed drawing. */
/** @typedef {{ data: string, size: number }} ConvertedResource The base64 bytes and size Python returned. */
// The app injects the Worker protein result conversion (OV-1145), like the configuration validator.
/** @type {((request: ProteinResultConversion) => Promise<unknown>) | null} */
let proteinResultConverter = null;
/** @param {typeof proteinResultConverter} converter */
export const setStoredProteinResultFrameConverter = (converter) => {
  proteinResultConverter = typeof converter === 'function' ? converter : null;
};

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
 * drop (OV-869). Indices are read as property keys, as `list[index]` reads them.
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
 * The file a Load converts protein results of: its version, which Python's
 * converter reads, and whether it keeps them in the drawn frame
 * (`SourceFacts.proteinResultsInDrawnFrame`; a Session 27-44 does).
 * @typedef {{ version: number, inDrawnFrame: boolean }} ProteinFrameSource
 */
/**
 * B-LOAD (D-128): the protein results that each committed drawing of a
 * Session up to 47 stored for a reversed record, in the search frame. Python
 * converts each drawing whose results touch one (one converter, no JS copy);
 * the bytes come back by stored resource ID. A failure throws, so the Load
 * changes nothing (PD-OI-046).
 * @param {ProteinFrameSource} source
 * @param {ReadonlyArray<StoredDrawing | null | undefined>} drawings
 * @param {{ convert?: typeof proteinResultConverter }} [options]
 * @returns {Promise<Map<string, ConvertedResource>>} empty when no drawing needs it
 */
export const convertStoredProteinResultFrames = async ({ version, inDrawnFrame }, drawings, { convert = proteinResultConverter } = {}) => {
  /** @type {Map<string, ConvertedResource>} */
  const converted = new Map();
  if (!inDrawnFrame) return converted;
  for (const drawing of drawings) {
    const resources = drawing?.resources ?? {};
    if (!drawing || !storesReversedProteinResult(drawing.renderRequest, resources)) continue;
    if (!convert) throw diagnosticError('HELPER_PROTOCOL', {}, PROTEIN_RESULT_CONVERSION);
    const response = await convert({ version, canonicalRequest: drawing.renderRequest, resources });
    for (const [resourceId, item] of readConvertedResources(response, PROTEIN_RESULT_CONVERSION)) {
      const previous = converted.get(resourceId);
      if (previous && previous.data !== item.data) throw diagnosticError('RESULT_INVALID', {}, PROTEIN_RESULT_CONVERSION);
      converted.set(resourceId, item);
    }
  }
  return converted;
};

/**
 * The bytes by resource ID of a Worker conversion reply, `{ result: { resources: { [id]: { data, size } } } }`,
 * as Python returns them. A reply of another shape fails as RESULT_INVALID of `failure`.
 * @param {unknown} response
 * @param {{ stage: 'helper', operation: string }} failure
 * @returns {Map<string, ConvertedResource>}
 */
export const readConvertedResources = (response, failure) => {
  const returned = /** @type {UnvalidatedFields<'result', UnvalidatedFields<'resources'>>} */ (response)?.result?.resources;
  if (!returned || typeof returned !== 'object' || Array.isArray(returned)) {
    throw diagnosticError('RESULT_INVALID', {}, failure);
  }
  /** @type {Map<string, ConvertedResource>} */
  const converted = new Map();
  for (const [resourceId, value] of Object.entries(returned)) {
    const item = /** @type {UnvalidatedFields<keyof ConvertedResource>} */ (value);
    if (typeof item?.data !== 'string' || !Number.isSafeInteger(item.size) || Number(item.size) < 0) {
      throw diagnosticError('RESULT_INVALID', {}, failure);
    }
    converted.set(resourceId, { data: item.data, size: Number(item.size) });
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
