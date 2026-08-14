import { collectCanonicalResourceIds } from './canonical-resource-references.js';
import { decodeDepthText, isEncodedDepthFileEntry } from './depth-file-codec.js';
import { takeFileBytesForTransfer } from './file-content-cache.js';
import { getResourcePayloadOwner } from './resource-payload-owner.js';
import {
  recordSessionLifecycleEvent,
  recordStructuralMetric
} from './runtime-test-hooks.js';

const BASE64_CHUNK_CHARACTERS = 1024 * 1024;

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

const validateResourceDescriptor = (resourceId, descriptor) => {
  if (!descriptor || typeof descriptor !== 'object' || Array.isArray(descriptor)) {
    throw new Error(`Canonical render resource ${resourceId} must be an object.`);
  }
  if (!Number.isSafeInteger(descriptor.size) || descriptor.size < 0) {
    throw new Error(`Canonical render resource ${resourceId} has an invalid byte size.`);
  }
  if (
    descriptor.encoding !== 'base64'
    && !isEncodedDepthFileEntry(descriptor)
  ) {
    throw new Error(`Canonical render resource ${resourceId} has an unsupported encoding.`);
  }
  return descriptor;
};

const descriptorTransferBuffer = async (resourceId, descriptor, owner) => {
  if (owner && owner !== descriptor) {
    return takeFileBytesForTransfer(owner);
  }
  if (isEncodedDepthFileEntry(descriptor)) {
    const bytes = new TextEncoder().encode(decodeDepthText(descriptor.data));
    recordStructuralMetric('decodedByteCount', bytes.byteLength, { resourceId });
    return bytes.buffer;
  }
  recordStructuralMetric('base64DecodeCount', 1, { resourceId });
  const buffer = decodeBase64ForTransfer(descriptor.data, descriptor.size, resourceId);
  recordStructuralMetric('decodedByteCount', buffer.byteLength, { resourceId });
  return buffer;
};

export const createDiagramResourceTransport = () => {
  let cachedResources = new Map();
  let nextCacheToken = 1;

  const reset = () => {
    cachedResources = new Map();
    nextCacheToken = 1;
  };

  const prepare = async ({ request, resources } = {}) => {
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
      const payloadIdentity = owner === descriptor ? descriptor.data : owner;
      const previous = cachedResources.get(resourceId);
      const cacheHit = Boolean(
        previous
        && previous.owner === owner
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
        owner,
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
    recordSessionLifecycleEvent('resource-stage-end', {
      referencedResourceCount: referencedIds.length,
      referencedDeclaredBytes,
      bytesTransferred,
      materializationCount
    });

    return {
      resourceManifest,
      stagedResources,
      nextCachedResources,
      commit() {
        cachedResources = nextCachedResources;
      }
    };
  };

  return { prepare, reset };
};
