// @ts-check
// origin/main Web Sessions (version 42 and older) stored the nucleotide
// comparison rows of a reverse-complemented Linear record after the reverse
// complement. Comparison tables now use the search frame (D-18, PD-OI-073),
// so Load rewrites those stored tables once; a re-saved Session keeps the
// search frame. The table is read and rewritten by Python (one table reader).
import { discoverGffFastaRecords, discoverSequenceRecords } from './record-discovery.js';
import { base64ToBytes, bytesToText, textToBase64, textToBytes } from './byte-utils.js';

/** @import { SequenceRecord } from './record-discovery.js' */

/** @typedef {{ kind: 'gffFasta', gffResourceId: string, fastaResourceId: string } | { kind: 'genbank', resourceId: string }} FrameSource */
/** @typedef {{ kind: 'recordIndex', index: number } | { kind: 'recordId', value?: unknown, id?: unknown }} FrameSelector */
/** @typedef {{ kind?: string, resourceId?: string, queryRecordIndex: number, subjectRecordIndex: number }} FrameComparison */
/**
 * The part of a committed request record the frame conversion reads.
 * @typedef {{
 *   source: FrameSource,
 *   selector?: FrameSelector | null,
 *   region?: { start: number, end: number, reverseComplement?: unknown } | null,
 *   presentation?: { reverseComplement?: unknown } | null
 * }} FrameRecord
 */
/** @typedef {{ data?: string, name?: string, type?: string, [field: string]: unknown }} FrameResource */
/** @typedef {Record<string, FrameResource>} FrameResources A Session's encoded resources by id. */
/** @typedef {{ version?: number, renderRequest?: { mode?: unknown, comparisons?: FrameComparison[], records?: FrameRecord[] } | null, resources: FrameResources }} FrameSession */
/**
 * The record discovery the conversion reads a record length from.
 * @typedef {{
 *   sequence: (options: { file: File, format: string }) => Promise<SequenceRecord[]>,
 *   gffFasta: (options: { gffFile: File, fastaFile: File }) => Promise<SequenceRecord[]>
 * }} RecordDiscovery
 */

export const MAIN_DISPLAY_FRAME_SESSION_VERSION = 42;

// The app injects the Worker table conversion, like the configuration validator.
/** @type {((request: Record<string, any>) => Promise<any>) | null} */
let tableConverter = null;
/** @param {typeof tableConverter} converter */
export const setMainSessionComparisonFrameConverter = (converter) => {
  tableConverter = typeof converter === 'function' ? converter : null;
};

/**
 * @param {FrameResources | undefined} resources
 * @param {string} resourceId
 */
const resourceFile = (resources, resourceId) => {
  const entry = resources?.[resourceId];
  if (!entry || typeof entry.data !== 'string') {
    throw new Error(`Missing Session resource: ${resourceId}`);
  }
  return new File([base64ToBytes(entry.data)], String(entry.name || resourceId), {
    type: String(entry.type || 'text/plain')
  });
};

/**
 * @param {FrameRecord} record
 * @param {FrameResources | undefined} resources
 * @param {RecordDiscovery} discover
 */
const selectedRecordLength = async (record, resources, discover) => {
  const source = record?.source || {};
  const catalog = source.kind === 'gffFasta'
    ? await discover.gffFasta({
        gffFile: resourceFile(resources, source.gffResourceId),
        fastaFile: resourceFile(resources, source.fastaResourceId)
      })
    : await discover.sequence({ file: resourceFile(resources, source.resourceId), format: 'genbank' });
  const selector = record?.selector || null;
  const selected = selector?.kind === 'recordIndex'
    ? catalog[selector.index]
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
 * @param {FrameRecord | undefined} record
 * @param {FrameResources | undefined} resources
 * @param {RecordDiscovery} discover
 */
const endpointFrame = async (record, resources, discover) => {
  const reverse = Boolean(record?.region?.reverseComplement || record?.presentation?.reverseComplement);
  if (!reverse) return { length: 0, reverse: false };
  const region = record?.region;
  const length = region
    ? Number(region.end) - Number(region.start) + 1
    // `reverse` is true only when `record` holds a reverse complement, so the record exists.
    : await selectedRecordLength(/** @type {FrameRecord} */ (record), resources, discover);
  return { length, reverse: true };
};

/**
 * Return the Session with each stored nucleotide table of a reversed endpoint
 * in the search frame. Other Sessions are returned unchanged.
 * @template {FrameSession} T
 * @param {T} data
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
  const request = data?.renderRequest;
  // `Number.isInteger` does not narrow `version`; the second operand runs only for an integer.
  if (!Number.isInteger(data?.version) || /** @type {number} */ (data.version) > MAIN_DISPLAY_FRAME_SESSION_VERSION
    || request?.mode !== 'linear' || !Array.isArray(request.comparisons)) {
    return data;
  }
  const records = Array.isArray(request.records) ? request.records : [];
  let resources = data.resources;
  const converted = new Set();
  for (const comparison of request.comparisons) {
    const resourceId = comparison?.resourceId;
    if (comparison?.kind !== 'nucleotideBlast' || !resourceId || converted.has(resourceId)) continue;
    const queryFrame = await endpointFrame(records[comparison.queryRecordIndex], resources, discover);
    const subjectFrame = await endpointFrame(records[comparison.subjectRecordIndex], resources, discover);
    if (!queryFrame.reverse && !subjectFrame.reverse) continue;
    if (!convertTable) throw new Error('Session comparison conversion service is unavailable.');
    const response = await convertTable(
      { tableText: bytesToText(base64ToBytes(resources[resourceId].data), { fatal: true }), queryFrame, subjectFrame }
    );
    const tsv = String(response?.result?.tsv ?? '');
    resources = {
      ...resources,
      [resourceId]: { ...resources[resourceId], encoding: 'base64', data: textToBase64(tsv), size: textToBytes(tsv).length }
    };
    converted.add(resourceId);
  }
  return converted.size > 0 ? { ...data, resources } : data;
};
