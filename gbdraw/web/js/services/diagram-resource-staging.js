// @ts-check
import { base64ToBytesInTasks } from './byte-utils.js';
import { collectCanonicalResourceIds } from './canonical-resource-references.js';
import { decodeDepthText, isEncodedDepthFileEntry } from './depth-file-codec.js';
import {
  readFileContentToken, releaseFileContentCache, resourceContentToken, takeFileBytesForTransfer
} from './file-content-cache.js';
import { getResourcePayloadIdentity, getResourcePayloadOwner } from './resource-payload-owner.js';
export { recordSourceOwnedBy, recordSourceResourceIds } from './resource-payload-owner.js';
import {
  recordSessionLifecycleEvent,
  recordStructuralMetric
} from './runtime-test-hooks.js';

/** @import { DepthFilePayload } from './depth-file-codec.js' */
/** @import { ResourceContentToken } from './file-content-cache.js' */
/** @import { ResourcePayloadOwner } from './resource-payload-owner.js' */
/** @import { FileLike, RawSessionResourceDescriptor } from './session-resource-backing.js' */

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
 * The files a drawing's draft holds now while its request omits them (a hidden
 * Depth file; `heldDepthFiles`), as the descriptors earlier requests carried.
 * @typedef {{ drawingId: string, files: readonly RawSessionResourceDescriptor[] }} DraftHeldFiles
 */
/**
 * Every drawing's draft held files, read at each staging, so a file leaves
 * the Worker once no draft holds it (gbdraw/web/CLAUDE.md, Runtime data flow).
 * The composition root supplies it.
 * @callback ReadDraftHeldFiles
 * @returns {Iterable<DraftHeldFiles>}
 */

/**
 * What names a payload's bytes until they change: its inline base64 text, else
 * the payload identity its owner gives (`getResourcePayloadIdentity`: the File
 * the user chose, or the source descriptor of every view of one Session
 * resource).
 * @param {ResourceDescriptor} descriptor
 * @returns {string | object}
 */
const stagingKeyOf = (descriptor) => (
  typeof descriptor.data === 'string' ? descriptor.data : getResourcePayloadIdentity(descriptor)
);

/** @param {{ readCommittedResourceSets?: ReadCommittedResourceSets, readDraftHeldFiles?: ReadDraftHeldFiles }} [ports] */
export const createDiagramResourceTransport = ({ readCommittedResourceSets = () => [], readDraftHeldFiles = () => [] } = {}) => {
  // The tokens the Worker holds after the last acknowledged request.
  /** @type {ReadonlySet<ResourceContentToken>} */
  let heldTokens = new Set();
  // The token of each payload identity seen whose bytes the Worker holds.
  /** @type {Map<string | object, ResourceContentToken>} */
  let knownTokens = new Map();
  // A new Worker holds nothing: a helper staging prepared before the reset records nothing.
  let workerEpoch = 0;

  const reset = () => {
    heldTokens = new Set();
    knownTokens = new Map();
    workerEpoch += 1;
  };

  // Adds the token of `descriptor` when the Worker holds its bytes. Only tokens
  // already known count: hashing a payload here would read bytes the Worker may not hold.
  /** @param {unknown} descriptor @param {Set<ResourceContentToken>} into */
  const addHeldToken = (descriptor, into) => {
    if (!descriptor || typeof descriptor !== 'object') return;
    const token = knownTokens.get(stagingKeyOf(/** @type {ResourceDescriptor} */ (descriptor)));
    if (token && heldTokens.has(token)) into.add(token);
  };

  // The held tokens a drawing's committed request names.
  const committedHeldTokens = () => {
    /** @type {Set<ResourceContentToken>} */
    const retained = new Set();
    for (const { request, resources } of readCommittedResourceSets()) {
      for (const resourceId of collectCanonicalResourceIds(request)) {
        addHeldToken(resources && Object.hasOwn(resources, resourceId) ? resources[resourceId] : null, retained);
      }
    }
    return retained;
  };

  // A render stages what its request names, and the Worker then holds only that.
  // A helper (`retain`) stages the resources it reads, and the Worker keeps the
  // rest for the next Generate (OV-1590), so the held tokens keep them too.
  // The files each drawing's draft holds now while its request omits them stay
  // too, and those of the drawing a render draws (`drawingId`) are its
  // `heldResourceTokens`; a file the Worker never received is not held, read or staged.
  /**
   * @param {{
   *   request?: unknown,
   *   resources?: Readonly<Record<string, Partial<ResourceDescriptor> | null | undefined>> | null,
   *   drawingId?: string | null,
   *   retain?: boolean
   * }} [payload]
   */
  const prepare = async ({ request, resources, drawingId = null, retain = false } = {}) => {
    const previousHeldTokens = heldTokens;
    const preparedEpoch = workerEpoch;
    const retained = committedHeldTokens();
    /** @type {Set<ResourceContentToken>} */
    const drawingHeldTokens = new Set();
    for (const draft of readDraftHeldFiles()) {
      for (const descriptor of draft.files) addHeldToken(descriptor, draft.drawingId === drawingId ? drawingHeldTokens : retained);
    }
    for (const token of drawingHeldTokens) retained.add(token);
    const referencedIds = Array.from(collectCanonicalResourceIds(request));
    const resourceManifest = [];
    /** @type {{ resourceId: string, cacheToken: ResourceContentToken, bytes: ArrayBuffer }[]} */
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
      const identity = stagingKeyOf(descriptor);
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
        // The input summary names the file the user chose (gbdraw/api/input_summary.py), not its Session name.
        name: owner instanceof File && owner.name ? owner.name : String(descriptor.name || `${resourceId}.bin`),
        kind: String(descriptor.kind || 'web-file'),
        type: String(descriptor.type || 'application/octet-stream'),
        size: descriptor.size,
        lastModified: Number(descriptor.lastModified) || 0
      });
      // Equal bytes under another resource ID of this request travel once, and bytes a helper
      // staged from a File stay held for the Generate that encodes that File (OV-1595).
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
      heldResourceTokens: [...drawingHeldTokens].filter((token) => !manifestTokens.has(token)),
      commit() {
        if (retain) {
          if (workerEpoch !== preparedEpoch) return false;
          heldTokens = new Set([...heldTokens, ...stagedResources.map(({ cacheToken }) => cacheToken)]);
          return true;
        }
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
