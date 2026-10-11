// @ts-check
import {
  bytesToBase64,
  getSessionResourceSource,
  readFileBytes,
  sha256Hex
} from './file-content-cache.js';
import {
  adoptCurrentSessionResources,
  createSessionResourceFileView,
  sessionResourceSource
} from './session-resource-backing.js';
import { hasBiologicalSessionInputs, isAdoptedCanonicalSession } from './session-authority.js';
import { recordStructuralMetric } from './runtime-test-hooks.js';
import {
  collectCanonicalResourceIds,
  isCanonicalResourceReferenceField
} from './canonical-resource-references.js';

/**
 * The shape `sessionResourceSource` returns: a single backing (`descriptor`,
 * `readBytes`) or a composite of backings (`descriptors`).
 * @typedef {{
 *   resourceId?: string,
 *   descriptor?: Record<string, any>,
 *   descriptors?: Record<string, any>[],
 *   readBytes?: () => any
 * }} SessionResourceSource
 */

const safeResourceLeaf = (value) => {
  // split() returns at least one element, so pop() is a string.
  const basename = /** @type {string} */ (String(value || 'resource.dat')
    .replace(/\\/g, '/')
    .split('/')
    .pop());
  const safe = basename
    .replace(/[^A-Za-z0-9._-]+/g, '_')
    .replace(/^[._]+|[._]+$/g, '');
  return safe || 'resource.dat';
};

const bindingForFile = (file, resourceId) => ({
  resourceId,
  name: file?.name === undefined ? 'file' : String(file.name),
  type: String(file?.type || ''),
  lastModified: Number(file?.lastModified) || 0
});

const sameEncodedPayload = (left, right) => (
  left === right || (
    left.encoding === right.encoding
    && left.data === right.data
    && Number(left.size) === Number(right.size)
  )
);

const rewriteResourceRefs = (value, aliases) => {
  if (Array.isArray(value)) {
    return value.map((item) => rewriteResourceRefs(item, aliases));
  }
  if (!value || typeof value !== 'object') return value;
  return Object.fromEntries(
    Object.entries(value).map(([key, item]) => [
      key,
      (
        isCanonicalResourceReferenceField(key)
        && typeof item === 'string'
        && aliases.has(item)
      )
        ? aliases.get(item)
        : rewriteResourceRefs(item, aliases)
    ])
  );
};

const rewriteOriginalNameHints = (webFiles, aliases) => {
  const identityAliases = [...aliases].every(([source, target]) => source === target);
  if (identityAliases) {
    const rewritten = { ...(webFiles || {}) };
    delete rewritten.bindings;
    if (Array.isArray(rewritten.linearRecordMetadata)) {
      rewritten.linearRecordMetadata = rewritten.linearRecordMetadata.map((entry) => {
        if (!entry || typeof entry !== 'object' || Array.isArray(entry)) return entry;
        const { losatFilename: _losatFilename, ...metadata } = entry;
        return metadata;
      });
    }
    return rewritten;
  }
  const rewritten = rewriteResourceRefs(webFiles || {}, aliases);
  delete rewritten.bindings;
  const hints = webFiles?.resourceOriginalNames;
  if (hints && typeof hints === 'object' && !Array.isArray(hints)) {
    rewritten.resourceOriginalNames = Object.fromEntries(
      Object.entries(hints)
        .filter(([resourceId]) => aliases.has(resourceId))
        .map(([resourceId, name]) => [aliases.get(resourceId), name])
    );
  }
  [
    'conservationLosatFastaSources',
    'conservationSequenceSources'
  ].forEach((field) => {
    if (!Array.isArray(webFiles?.[field])) return;
    rewritten[field] = webFiles[field].map((resourceId) => (
      resourceId && aliases.has(resourceId) ? aliases.get(resourceId) : null
    ));
  });
  if (Array.isArray(rewritten.linearRecordMetadata)) {
    rewritten.linearRecordMetadata = rewritten.linearRecordMetadata.map((entry) => {
      if (!entry || typeof entry !== 'object' || Array.isArray(entry)) return entry;
      const { losatFilename: _losatFilename, ...metadata } = entry;
      return metadata;
    });
  }
  return rewritten;
};

// The request metadata of `webFiles` that one committed request writes: Save
// keeps the other mode's (E1, `otherModeResult`) where the top-level request
// has none. Circular and Linear requests write different fields.
const REQUEST_WEB_FILE_FIELDS = [
  'linearRecordMetadata', 'conservationLosatFastaSources', 'conservationSequenceSources',
  'circularInputOriginalName', 'conservationBlastSource'
];
const mergeOtherRequestWebFiles = (webFiles, otherWebFiles, aliases) => {
  if (!otherWebFiles || typeof otherWebFiles !== 'object') return webFiles;
  const other = rewriteOriginalNameHints(otherWebFiles, aliases);
  const merged = { ...webFiles };
  REQUEST_WEB_FILE_FIELDS.forEach((field) => {
    if (!Object.hasOwn(merged, field) && Object.hasOwn(other, field)) merged[field] = other[field];
  });
  const names = { ...(other.resourceOriginalNames || {}), ...(webFiles.resourceOriginalNames || {}) };
  if (Object.keys(names).length > 0) merged.resourceOriginalNames = names;
  return merged;
};

/**
 * One resource table for the committed request and, when the other diagram
 * mode keeps a Result too (E1), that mode's committed request: both are
 * allocated through the same content-deduplicating allocator, so equal bytes
 * share one resource and a positional ID that names different bytes in the
 * two requests gets a new ID; each request's references are rewritten.
 * `state` holds the Files and the Linear records; `drawing` holds the Linear
 * comparison plan, whose uploaded comparison files the Session also binds.
 * @param {Record<string, any>} state
 * @param {Record<string, any> | null} committedRequest
 * @param {Record<string, any>} [drawing]
 * @param {Record<string, any> | null} [otherCommittedRequest]
 */
export const buildSessionResources = async (state, committedRequest, drawing, otherCommittedRequest = null) => {
  if (
    committedRequest === null
      ? hasBiologicalSessionInputs({ ...state?.files, linearSeqs: state?.linearSeqs })
        || Boolean(state?.results?.value?.length)
      : (!committedRequest
        || typeof committedRequest !== 'object'
        || !committedRequest.renderRequest
        || !committedRequest.resources)
  ) {
    throw new Error('A committed canonical render request is required to save a session.');
  }

  const resources = {};
  const aliases = new Map();
  const reuseEncodedResources = isAdoptedCanonicalSession(committedRequest);
  // A textual ID is canonical only within its own source table.
  const committedTable = adoptCurrentSessionResources(committedRequest?.resources || {});
  const candidatesBySize = new Map();
  const candidatesByDescriptor = new WeakMap();
  const usedNames = new Set();
  // Resources of an adopted table that nothing names yet: they keep their IDs
  // and names, and are written only when a binding or request uses them.
  /** @type {Map<string, Record<string, any>>} */
  const reserved = new Map();
  /** @param {string} id */
  const taken = (id) => Object.hasOwn(resources, id) || reserved.has(id);
  /** @param {string} id */
  const use = (id) => {
    const descriptor = reserved.get(id);
    if (descriptor) {
      resources[id] = descriptor;
      reserved.delete(id);
    }
    return id;
  };
  let nextResourceNumber = 1;
  let allocation = Promise.resolve();

  const nextResourceId = () => {
    let resourceId;
    do {
      resourceId = `resource-${String(nextResourceNumber++).padStart(4, '0')}`;
    } while (taken(resourceId));
    return resourceId;
  };
  const candidateFor = (source) => {
    const key = source.descriptor || source.bytes;
    if (candidatesByDescriptor.has(key)) {
      return candidatesByDescriptor.get(key);
    }
    const candidate = {
      ...source,
      size: source.descriptor?.size ?? source.bytes.byteLength,
      bytesPromise: null,
      identityPromise: null
    };
    candidatesByDescriptor.set(key, candidate);
    return candidate;
  };
  const bytesFor = candidate => {
    candidate.bytesPromise ??= Promise.resolve().then(() => candidate.bytes || candidate.readBytes());
    return candidate.bytesPromise;
  };
  const identityFor = candidate => {
    candidate.identityPromise ??= bytesFor(candidate).then(bytes => {
      // readBytes has already validated any declared checksum.
      if (candidate.descriptor?.checksum) {
        return candidate.descriptor.checksum.trim().toLowerCase().replace(/^sha256:/, '');
      }
      recordStructuralMetric('resourceIdentityHashCount', 1);
      return sha256Hex(bytes);
    });
    return candidate.identityPromise;
  };
  const register = (id, candidate) => {
    let bucket = candidatesBySize.get(candidate.size);
    if (!bucket) {
      bucket = { encoded: new Map(), identities: new Map(), pending: new Map() };
      candidatesBySize.set(candidate.size, bucket);
    }
    const descriptor = candidate.descriptor;
    if (descriptor) {
      let encoded = bucket.encoded.get(descriptor.encoding);
      if (!encoded) bucket.encoded.set(descriptor.encoding, encoded = new Map());
      if (!encoded.has(descriptor.data)) encoded.set(descriptor.data, id);
    }
    bucket.pending.set(id, candidate);
    usedNames.add((resources[id] || reserved.get(id)).name);
  };
  const equivalentId = async candidate => {
    const bucket = candidatesBySize.get(candidate.size);
    if (!bucket) return null;
    const descriptor = candidate.descriptor;
    const encodedId = descriptor && bucket.encoded.get(descriptor.encoding)?.get(descriptor.data);
    if (encodedId) return encodedId;
    const identity = await identityFor(candidate);
    if (bucket.identities.has(identity)) return bucket.identities.get(identity);
    for (const [id, existing] of bucket.pending) {
      const existingIdentity = await identityFor(existing);
      bucket.pending.delete(id);
      if (!bucket.identities.has(existingIdentity)) bucket.identities.set(existingIdentity, id);
      if (existingIdentity === identity) return id;
    }
    return null;
  };
  const allocate = (source, metadata, preferredId = '') => {
    // File arrays may resolve concurrently; allocation itself is one ordered
    // transaction so two equal sources cannot publish duplicate payloads.
    allocation = allocation.then(async () => {
      const candidate = candidateFor(source);
      const existing = preferredId && (resources[preferredId] || reserved.get(preferredId));
      if (source.descriptor?.checksum && source.descriptor !== existing
        && source.descriptor.checksum !== existing?.checksum) await bytesFor(candidate);
      if (existing && source.descriptor && sameEncodedPayload(existing, source.descriptor)) {
        return use(preferredId);
      }
      const equivalent = await equivalentId(candidate);
      if (equivalent) return use(equivalent);
      // Preserve the existing collision integrity check through the lazy backing.
      if (existing) await bytesFor(candidate);
      const descriptor = source.descriptor;
      const safePreferred = /^[a-z][a-z0-9]*(?:-[a-z0-9]+)*$/.test(preferredId)
        && !taken(preferredId)
        && descriptor?.name === safeResourceLeaf(descriptor.name)
        && !usedNames.has(descriptor.name);
      let id = safePreferred ? preferredId : nextResourceId();
      let name = safePreferred ? descriptor.name : `${id}-${safeResourceLeaf(metadata.name)}`;
      while (usedNames.has(name)) {
        id = nextResourceId();
        name = `${id}-${safeResourceLeaf(metadata.name)}`;
      }
      if (descriptor) {
        resources[id] = safePreferred ? descriptor : { ...descriptor, name };
      } else {
        const bytes = await bytesFor(candidate);
        recordStructuralMetric('base64EncodeCount', 1, { resourceName: metadata.name });
        recordStructuralMetric('encodedByteCount', bytes.byteLength, { resourceName: metadata.name });
        resources[id] = {
          kind: String(metadata.kind || 'web-file'), name,
          type: String(metadata.type || 'application/octet-stream'),
          size: bytes.byteLength, lastModified: Number(metadata.lastModified) || 0,
          encoding: 'base64', data: bytesToBase64(bytes)
        };
      }
      register(id, candidate);
      return id;
    });
    return allocation;
  };

  const committedResources = committedRequest?.resources || {};
  if (reuseEncodedResources) {
    // Only what the committed request and its Web metadata name is kept as
    // is; the bindings and the other mode's request (E1) are allocated below,
    // so a resource that nothing names any more is not written.
    const webFiles = committedRequest?.webFiles || {};
    const named = collectCanonicalResourceIds(committedRequest?.renderRequest);
    collectCanonicalResourceIds(webFiles, named);
    ['conservationLosatFastaSources', 'conservationSequenceSources'].forEach((field) => {
      (Array.isArray(webFiles[field]) ? webFiles[field] : [])
        .forEach((/** @type {unknown} */ id) => { if (typeof id === 'string' && id) named.add(id); });
    });
    Object.entries(committedResources).forEach(([id, descriptor]) => {
      if (named.has(id)) {
        resources[id] = descriptor;
        aliases.set(id, id);
      } else {
        reserved.set(id, descriptor);
      }
      register(id, candidateFor(sessionResourceSource(createSessionResourceFileView(committedTable, id))));
    });
  }
  for (const id of collectCanonicalResourceIds(committedRequest?.renderRequest)) {
    if (!Object.hasOwn(committedResources, id)) {
      throw new Error(`Committed render resource is missing: ${id}.`);
    }
    // A view of one resource id is a single backing, which always has readBytes.
    const source = /** @type {SessionResourceSource & { readBytes: () => any }} */ (
      sessionResourceSource(createSessionResourceFileView(committedTable, id))
    );
    if (!reuseEncodedResources) await source.readBytes();
    aliases.set(id, reuseEncodedResources ? id : await allocate(source, source.descriptor));
  }
  const otherAliases = new Map();
  if (otherCommittedRequest) {
    const reuseOther = isAdoptedCanonicalSession(otherCommittedRequest);
    const otherTable = otherCommittedRequest.resources === committedRequest?.resources
      ? committedTable
      : adoptCurrentSessionResources(otherCommittedRequest.resources || {});
    for (const id of collectCanonicalResourceIds(otherCommittedRequest.renderRequest)) {
      if (!Object.hasOwn(otherCommittedRequest.resources, id)) {
        throw new Error(`Committed render resource is missing: ${id}.`);
      }
      const source = /** @type {SessionResourceSource & { readBytes: () => any }} */ (
        sessionResourceSource(createSessionResourceFileView(otherTable, id))
      );
      if (!reuseOther) await source.readBytes();
      otherAliases.set(id, await allocate(source, source.descriptor, reuseOther ? id : ''));
    }
  }

  const bindSource = async (source, metadata) => bindingForFile(
    metadata, await allocate(source, metadata, source.resourceId)
  );
  const bindFile = async file => {
    if (!file) return null;
    const source = /** @type {SessionResourceSource | null} */ (getSessionResourceSource(file));
    if (Array.isArray(source?.descriptors)) {
      const components = [];
      for (const component of source.descriptors) components.push(await bindSource(component, component));
      const { resourceId: _resourceId, ...metadata } = bindingForFile(file, '');
      return { kind: 'composite', components, ...metadata };
    }
    if (source?.descriptor) return bindSource(source, file);
    return bindSource({ bytes: await readFileBytes(file) }, file);
  };

  const bindFileValue = async (value) => {
    if (Array.isArray(value)) {
      return Promise.all(value.map((item) => bindFileValue(item)));
    }
    return bindFile(value);
  };

  const files = state?.files || {};
  const linearSeqs = Array.isArray(state?.linearSeqs) ? state.linearSeqs : [];
  const linearComparisons = Array.isArray(drawing?.linearComparisonPlan?.edges)
    ? drawing.linearComparisonPlan.edges
    : [];

  const bindings = {
    schema: 2,
    c_gb: await bindFile(files.c_gb),
    c_gff: await bindFile(files.c_gff),
    c_fasta: await bindFile(files.c_fasta),
    c_depth: await bindFileValue(files.c_depth),
    c_conservation_blasts: await bindFileValue(files.c_conservation_blasts),
    c_conservation_blasts_source:
      files.c_conservation_blasts_source === 'losat-cache' ? 'losat-cache' : null,
    c_conservation_fastas: await bindFileValue(files.c_conservation_fastas),
    c_conservation_sequence_sources: await bindFileValue(
      files.c_conservation_sequence_sources
    ),
    d_color: await bindFile(files.d_color),
    t_color: await bindFile(files.t_color),
    blacklist: await bindFile(files.blacklist),
    whitelist: await bindFile(files.whitelist),
    qualifier_priority: await bindFile(files.qualifier_priority),
    linearSeqs: await Promise.all(linearSeqs.map(async (sequence) => ({
      uid: String(sequence?.uid || ''),
      gb: await bindFile(sequence?.gb),
      gff: await bindFile(sequence?.gff),
      fasta: await bindFile(sequence?.fasta),
      depth: await bindFileValue(sequence?.depth),
      losat_gencode: sequence?.losat_gencode ?? 1,
      definition: String(sequence?.definition || ''),
      record_subtitle: String(sequence?.record_subtitle || ''),
      file_definition: String(sequence?.file_definition || ''),
      file_subtitle: String(sequence?.file_subtitle || ''),
      inferred_definition: String(sequence?.inferred_definition || ''),
      region_record_id: String(sequence?.region_record_id || ''),
      region_start: sequence?.region_start ?? null,
      region_end: sequence?.region_end ?? null,
      region_reverse: Boolean(sequence?.region_reverse)
    }))),
    linearComparisons: await Promise.all(
      linearComparisons
        .filter((comparison) => comparison?.file)
        .map(async (comparison) => ({
          id: String(comparison?.id || ''),
          file: await bindFile(comparison.file)
        }))
    )
  };

  const identityAliases = [...aliases].every(([source, target]) => source === target);
  const otherIdentityAliases = [...otherAliases].every(([source, target]) => source === target);
  const webFiles = rewriteOriginalNameHints(committedRequest?.webFiles, aliases);

  return {
    renderRequest: identityAliases
      ? committedRequest?.renderRequest ?? null
      : rewriteResourceRefs(committedRequest?.renderRequest, aliases),
    ...(otherCommittedRequest ? {
      otherRenderRequest: otherIdentityAliases
        ? otherCommittedRequest.renderRequest
        : rewriteResourceRefs(otherCommittedRequest.renderRequest, otherAliases)
    } : {}),
    resources,
    webFiles: {
      ...(otherCommittedRequest
        ? mergeOtherRequestWebFiles(webFiles, otherCommittedRequest.webFiles, otherAliases)
        : webFiles),
      bindings
    }
  };
};
