// @ts-check
import { getSessionResourceSource } from './file-content-cache.js';
import { recordSourceOf } from './input-types.js';
/** @import { CanonicalRecordSource } from '../contracts/request.generated.js' */

/**
 * What owns a resource's payload bytes: the File the user chose, or another
 * object (the descriptor itself when no owner was set).
 * @typedef {WeakKey} ResourcePayloadOwner
 */

/** @type {WeakMap<object, ResourcePayloadOwner>} */
const payloadOwners = new WeakMap();

/** @param {unknown} value */
const objectKey = (value) => (
  value && (typeof value === 'object' || typeof value === 'function') ? value : null
);

/**
 * @template T
 * @param {T} descriptor
 * @param {unknown} owner
 * @returns {T}
 */
export const setResourcePayloadOwner = (descriptor, owner) => {
  const descriptorKey = objectKey(descriptor);
  const ownerKey = objectKey(owner);
  if (descriptorKey && ownerKey) payloadOwners.set(descriptorKey, ownerKey);
  return descriptor;
};

/**
 * @template T
 * @param {T} descriptor
 * @returns {T | ResourcePayloadOwner | null}
 */
export const getResourcePayloadOwner = (descriptor) => {
  const key = objectKey(descriptor);
  return (key && payloadOwners.get(key)) || key;
};

/**
 * The object whose bytes a File, a Session resource view or a descriptor
 * stands for: every view of one Session resource stands for its source
 * descriptor, anything else for its payload owner.
 * @param {object} value
 * @returns {object}
 */
export const getResourcePayloadIdentity = (value) => {
  const owner = /** @type {object} */ (getResourcePayloadOwner(value));
  return /** @type {{ descriptor?: object } | null} */ (getSessionResourceSource(owner))?.descriptor || owner;
};

/**
 * The id of the resource of `resources` whose bytes `file` stands for (the
 * File the user chose, or a Session resource view), of this kind, or null. A
 * request that names it stages the bytes the Worker holds under that id
 * (OV-1555).
 * @param {Readonly<Record<string, unknown>> | null | undefined} resources A committed Session's resources.
 * @param {object | null | undefined} file
 * @param {string} kind
 * @returns {string | null}
 */
export const resourceIdOwnedBy = (resources, file, kind) => {
  if (!resources || !file) return null;
  const identity = getResourcePayloadIdentity(file);
  return Object.keys(resources).find((resourceId) => {
    const descriptor = resources[resourceId];
    return typeof descriptor === 'object' && descriptor !== null
      && /** @type {{ kind?: unknown }} */ (descriptor).kind === kind && getResourcePayloadIdentity(descriptor) === identity;
  }) ?? null;
};

/**
 * The source a request names a record by when it draws these files: the files
 * of a record source of `sourceKind` (a GenBank file, a GFF3 file and its FASTA
 * file `fastaFile`, or a FASTA file alone), each by the id of the resource of
 * `resources` that holds it (`resourceIdOwnedBy`); null when one of them holds
 * no file. Both Generate's staging and the live record names read a card's
 * files through it.
 * @param {Readonly<Record<string, unknown>> | null | undefined} resources
 * @param {CanonicalRecordSource['kind']} sourceKind
 * @param {object | null | undefined} sourceFile
 * @param {object | null | undefined} fastaFile
 * @returns {CanonicalRecordSource | null}
 */
export const recordSourceOwnedBy = (resources, sourceKind, sourceFile, fastaFile) => recordSourceOf(
  sourceKind, ({ resource }, index) => resourceIdOwnedBy(resources, index === 0 ? sourceFile : fastaFile, resource)
);

// The resource ids of a record source's files, in source order (input-types.js owns the reading).
export { recordSourceResourceIds } from './input-types.js';
