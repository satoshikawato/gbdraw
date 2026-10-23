// @ts-check
import { base64ToBytesInTasks } from './byte-utils.js';
import { collectCanonicalResourceIds } from './canonical-resource-references.js';
import { decodeDepthText, isEncodedDepthFileEntry } from './depth-file-codec.js';
import { takeFileBytesForTransfer } from './file-content-cache.js';
import { getResourcePayloadIdentity, getResourcePayloadOwner } from './resource-payload-owner.js';
import {
  recordSessionLifecycleEvent,
  recordStructuralMetric
} from './runtime-test-hooks.js';

/** @import { DepthFilePayload } from './depth-file-codec.js' */
/** @import { ResourcePayloadOwner } from './resource-payload-owner.js' */
/** @import { FileLike } from './session-resource-backing.js' */

/**
 * A canonical render resource as a Session or request holds it: base64 text,
 * or an encoded depth table, with the file's metadata.
 * @typedef {{
 *   size: number, encoding?: unknown, data?: string | DepthFilePayload | null,
 *   name?: unknown, kind?: unknown, type?: unknown, lastModified?: unknown
 * }} ResourceDescriptor
 * @typedef {{ encoding?: unknown, data?: DepthFilePayload | null }} DepthResourceView
 *   A descriptor as the depth codec reads it; base64 `data` is a string, which never matches the codec's `schema`.
 */

// The base64 text of a resource is decoded by the one decoder a Session
// resource backing also uses: native chunks that yield to paint where the
// engine has them. `prepare` checks the decoded size against the descriptor.
/**
 * @param {unknown} encoded
 * @param {string} resourceId
 * @returns {Promise<ArrayBuffer>}
 */
const decodeBase64ForTransfer = async (encoded, resourceId) => {
  try {
    const bytes = await base64ToBytesInTasks(/** @type {string | null | undefined} */ (encoded));
    // The decoder writes into a buffer of exactly the decoded length.
    return /** @type {ArrayBuffer} */ (bytes.buffer);
  } catch (error) {
    throw new Error(`Canonical resource ${resourceId} contains invalid base64 data.`, { cause: error });
  }
};

/**
 * @param {string} resourceId
 * @param {Partial<ResourceDescriptor> | null | undefined} descriptor
 * @returns {ResourceDescriptor}
 */
const validateResourceDescriptor = (resourceId, descriptor) => {
  if (!descriptor || typeof descriptor !== 'object' || Array.isArray(descriptor)) {
    throw new Error(`Canonical render resource ${resourceId} must be an object.`);
  }
  if (!Number.isSafeInteger(descriptor.size) || /** @type {number} */ (descriptor.size) < 0) {
    throw new Error(`Canonical render resource ${resourceId} has an invalid byte size.`);
  }
  if (
    descriptor.encoding !== 'base64'
    && !isEncodedDepthFileEntry(/** @type {DepthResourceView} */ (descriptor))
  ) {
    throw new Error(`Canonical render resource ${resourceId} has an unsupported encoding.`);
  }
  // The checks above hold the object, its integer `size` and its `encoding`.
  return /** @type {ResourceDescriptor} */ (descriptor);
};

/**
 * @param {string} resourceId
 * @param {ResourceDescriptor} descriptor
 * @param {ResourcePayloadOwner | ResourceDescriptor | null} owner
 */
const descriptorTransferBuffer = async (resourceId, descriptor, owner) => {
  if (owner && owner !== descriptor) {
    // An owner other than the descriptor is the File the user chose.
    return takeFileBytesForTransfer(/** @type {FileLike} */ (owner));
  }
  if (isEncodedDepthFileEntry(/** @type {DepthResourceView} */ (descriptor))) {
    // The depth entry check above holds only for a `data` object.
    const bytes = new TextEncoder().encode(decodeDepthText(/** @type {DepthFilePayload} */ (descriptor.data)));
    recordStructuralMetric('decodedByteCount', bytes.byteLength, { resourceId });
    return bytes.buffer;
  }
  recordStructuralMetric('base64DecodeCount', 1, { resourceId });
  const buffer = await decodeBase64ForTransfer(descriptor.data, resourceId);
  recordStructuralMetric('decodedByteCount', buffer.byteLength, { resourceId });
  return buffer;
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

export const createDiagramResourceTransport = () => {
  let cachedResources = new Map();
  let nextCacheToken = 1;
  // A new Worker holds nothing: a helper staging prepared before the reset records nothing.
  let workerEpoch = 0;

  const reset = () => {
    cachedResources = new Map();
    nextCacheToken = 1;
    workerEpoch += 1;
  };

  // A render stages what its request names, and the Worker then holds only that.
  // A helper (`retain`) stages the resources it reads, and the Worker keeps the
  // rest for the next Generate (OV-1590), so the cache keeps them too.
  /**
   * @param {{
   *   request?: unknown,
   *   resources?: Readonly<Record<string, Partial<ResourceDescriptor> | null | undefined>> | null,
   *   retain?: boolean
   * }} [payload]
   */
  const prepare = async ({ request, resources, retain = false } = {}) => {
    const previousCachedResources = cachedResources;
    const preparedEpoch = workerEpoch;
    const referencedIds = Array.from(collectCanonicalResourceIds(request));
    const resourceManifest = [];
    const stagedResources = [];
    const nextCachedResources = new Map();
    let referencedDeclaredBytes = 0;
    let bytesTransferred = 0;
    let materializationCount = 0;

    recordSessionLifecycleEvent('resource-stage-start', {
      referencedResourceCount: referencedIds.length
    });
    for (const resourceId of referencedIds) {
      const descriptor = validateResourceDescriptor(resourceId, resources?.[resourceId]);
      const owner = getResourcePayloadOwner(descriptor);
      const payloadIdentity = typeof descriptor.data === 'string'
        ? descriptor.data
        : owner;
      const previous = previousCachedResources.get(resourceId);
      const cacheHit = Boolean(
        previous
        && previous.payloadIdentity === payloadIdentity
        && previous.size === descriptor.size
      );
      const cacheToken = cacheHit
        ? previous.cacheToken
        : `render-resource-${nextCacheToken++}`;
      referencedDeclaredBytes += descriptor.size;
      resourceManifest.push({
        resourceId,
        cacheToken,
        name: String(descriptor.name || `${resourceId}.bin`),
        kind: String(descriptor.kind || 'web-file'),
        type: String(descriptor.type || 'application/octet-stream'),
        size: descriptor.size,
        lastModified: Number(descriptor.lastModified) || 0
      });
      nextCachedResources.set(resourceId, {
        cacheToken,
        payloadIdentity,
        size: descriptor.size
      });
      if (cacheHit) {
        recordStructuralMetric('workerResourceCacheHitCount', 1, { resourceId });
        continue;
      }

      const bytes = await descriptorTransferBuffer(resourceId, descriptor, owner);
      if (!(bytes instanceof ArrayBuffer) || bytes.byteLength !== descriptor.size) {
        throw new Error(`Canonical render resource ${resourceId} byte size does not match.`);
      }
      materializationCount += 1;
      bytesTransferred += bytes.byteLength;
      recordStructuralMetric('resourceMaterializationCount', 1, { resourceId });
      recordStructuralMetric('transferredResourceBytes', bytes.byteLength, { resourceId });
      stagedResources.push({ resourceId, cacheToken, bytes });
    }
    const stagedResourceIds = stagedResources.map(({ resourceId }) => resourceId);
    recordSessionLifecycleEvent('resource-stage-end', {
      referencedResourceCount: referencedIds.length,
      referencedDeclaredBytes,
      bytesTransferred,
      materializationCount
    });

    return {
      resourceManifest,
      stagedResources,
      commit() {
        if (retain) {
          if (workerEpoch !== preparedEpoch) return false;
          stagedResourceIds.forEach((resourceId) => {
            cachedResources.set(resourceId, nextCachedResources.get(resourceId));
          });
          return true;
        }
        // A reset or another acknowledged request supersedes this preparation.
        if (cachedResources !== previousCachedResources) return false;
        cachedResources = nextCachedResources;
        return true;
      }
    };
  };

  return { prepare, reset };
};
