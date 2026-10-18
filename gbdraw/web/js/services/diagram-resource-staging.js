// @ts-check
import { collectCanonicalResourceIds } from './canonical-resource-references.js';
import { decodeDepthText, isEncodedDepthFileEntry } from './depth-file-codec.js';
import {
  readFileContentToken, releaseFileContentCache, resourceContentToken, takeFileBytesForTransfer
} from './file-content-cache.js';
import { getResourcePayloadOwner } from './resource-payload-owner.js';
import {
  recordSessionLifecycleEvent,
  recordStructuralMetric
} from './runtime-test-hooks.js';

/** @import { DepthFilePayload } from './depth-file-codec.js' */
/** @import { ResourceContentToken } from './file-content-cache.js' */
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

const BASE64_CHUNK_CHARACTERS = 1024 * 1024;

/**
 * @param {unknown} encoded
 * @param {number} declaredSize
 * @param {string} resourceId
 */
const decodeBase64ForTransfer = (encoded, declaredSize, resourceId) => {
  const source = String(encoded || '');
  if (source.length % 4 !== 0 || /\s/.test(source)) {
    throw new Error(`Canonical resource ${resourceId} contains invalid base64 data.`);
  }
  const bytes = new Uint8Array(declaredSize);
  let byteOffset = 0;
  for (let offset = 0; offset < source.length; offset += BASE64_CHUNK_CHARACTERS) {
    const chunk = atob(source.slice(offset, offset + BASE64_CHUNK_CHARACTERS));
    if (byteOffset + chunk.length > bytes.byteLength) {
      throw new Error(`Canonical resource ${resourceId} exceeds its declared byte size.`);
    }
    for (let index = 0; index < chunk.length; index += 1) {
      bytes[byteOffset + index] = chunk.charCodeAt(index);
    }
    byteOffset += chunk.length;
  }
  if (byteOffset !== declaredSize) {
    throw new Error(`Canonical resource ${resourceId} byte size does not match its descriptor.`);
  }
  return bytes.buffer;
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
  const buffer = decodeBase64ForTransfer(descriptor.data, descriptor.size, resourceId);
  recordStructuralMetric('decodedByteCount', buffer.byteLength, { resourceId });
  return buffer;
};

/**
 * A request a drawing has committed, with the resources it names.
 * @typedef {{ request: unknown, resources: Readonly<Record<string, unknown>> | null | undefined }} CommittedResourceSet
 */
/**
 * The committed request of every drawing; the Worker keeps the resources they
 * name. The composition root supplies it.
 * @callback ReadCommittedResourceSets
 * @returns {Iterable<CommittedResourceSet>}
 */

/**
 * What names a payload's bytes until they change: the base64 text, or the
 * object that owns them (the File the user chose, a Session resource view, or
 * the descriptor itself).
 * @param {ResourceDescriptor} descriptor
 * @returns {string | ResourcePayloadOwner}
 */
const payloadIdentityOf = (descriptor) => (
  typeof descriptor.data === 'string' ? descriptor.data : getResourcePayloadOwner(descriptor) ?? descriptor
);

/** @param {{ readCommittedResourceSets?: ReadCommittedResourceSets }} [ports] */
export const createDiagramResourceTransport = ({ readCommittedResourceSets = () => [] } = {}) => {
  // The tokens the Worker holds after the last acknowledged request.
  /** @type {ReadonlySet<ResourceContentToken>} */
  let heldTokens = new Set();
  // The token of each payload identity seen whose bytes the Worker holds.
  /** @type {Map<string | ResourcePayloadOwner, ResourceContentToken>} */
  let knownTokens = new Map();

  const reset = () => {
    heldTokens = new Set();
    knownTokens = new Map();
  };

  // The held tokens a drawing's committed request names. Only tokens already
  // known count: hashing a payload here would read bytes the Worker may not hold.
  const committedHeldTokens = () => {
    /** @type {Set<ResourceContentToken>} */
    const retained = new Set();
    for (const { request, resources } of readCommittedResourceSets()) {
      for (const resourceId of collectCanonicalResourceIds(request)) {
        const descriptor = resources && Object.hasOwn(resources, resourceId) ? resources[resourceId] : null;
        if (!descriptor || typeof descriptor !== 'object') continue;
        const token = knownTokens.get(payloadIdentityOf(/** @type {ResourceDescriptor} */ (descriptor)));
        if (token && heldTokens.has(token)) retained.add(token);
      }
    }
    return retained;
  };

  /** @param {{ request?: unknown, resources?: Readonly<Record<string, Partial<ResourceDescriptor> | null | undefined>> | null }} [payload] */
  const prepare = async ({ request, resources } = {}) => {
    const previousHeldTokens = heldTokens;
    const retained = committedHeldTokens();
    const referencedIds = Array.from(collectCanonicalResourceIds(request));
    const resourceManifest = [];
    const stagedResources = [];
    /** @type {Set<ResourceContentToken>} */
    const manifestTokens = new Set();
    let referencedDeclaredBytes = 0;
    let bytesTransferred = 0;
    let materializationCount = 0;

    recordSessionLifecycleEvent('resource-stage-start', {
      referencedResourceCount: referencedIds.length
    });
    for (const resourceId of referencedIds) {
      const descriptor = validateResourceDescriptor(resourceId, resources?.[resourceId]);
      const owner = getResourcePayloadOwner(descriptor);
      const identity = payloadIdentityOf(descriptor);
      /** @type {ArrayBufferLike | null} */
      let bytes = null;
      let cacheToken = knownTokens.get(identity);
      let hashedFile = false;
      if (!cacheToken) {
        if (owner && owner !== descriptor) {
          hashedFile = true;
          // The File's digest, cached on it: its bytes stay cached until the transfer below takes them.
          cacheToken = await readFileContentToken(/** @type {FileLike} */ (owner));
        } else {
          bytes = await descriptorTransferBuffer(resourceId, descriptor, owner);
          cacheToken = await resourceContentToken(new Uint8Array(bytes));
        }
        knownTokens.set(identity, cacheToken);
      }
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
      // Equal bytes under another resource ID of this request travel once.
      const cacheHit = previousHeldTokens.has(cacheToken) || manifestTokens.has(cacheToken);
      manifestTokens.add(cacheToken);
      if (cacheHit) {
        // A File read for its digest alone gives its bytes back, as a transfer would.
        if (hashedFile) releaseFileContentCache(/** @type {FileLike} */ (owner));
        recordStructuralMetric('workerResourceCacheHitCount', 1, { resourceId });
        continue;
      }

      bytes ??= await descriptorTransferBuffer(resourceId, descriptor, owner);
      if (!(bytes instanceof ArrayBuffer) || bytes.byteLength !== descriptor.size) {
        throw new Error(`Canonical render resource ${resourceId} byte size does not match.`);
      }
      materializationCount += 1;
      bytesTransferred += bytes.byteLength;
      recordStructuralMetric('resourceMaterializationCount', 1, { resourceId });
      recordStructuralMetric('transferredResourceBytes', bytes.byteLength, { resourceId });
      stagedResources.push({ resourceId, cacheToken, bytes });
    }
    recordSessionLifecycleEvent('resource-stage-end', {
      referencedResourceCount: referencedIds.length,
      referencedDeclaredBytes,
      bytesTransferred,
      materializationCount
    });

    return {
      resourceManifest,
      stagedResources,
      // The held tokens the Worker keeps besides this request's: it releases the rest.
      retainedResourceTokens: [...retained].filter((token) => !manifestTokens.has(token)),
      commit() {
        // A reset or another acknowledged request supersedes this preparation.
        if (heldTokens !== previousHeldTokens) return false;
        heldTokens = new Set([...manifestTokens, ...retained]);
        knownTokens = new Map([...knownTokens].filter(([, token]) => heldTokens.has(token)));
        return true;
      }
    };
  };

  return { prepare, reset };
};
