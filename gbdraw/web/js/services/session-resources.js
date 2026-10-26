// @ts-check
import {
  bytesToBase64,
  encodedFilePayloadOf,
  getSessionResourceSource,
  readFileBytes,
  sha256Hex
} from './file-content-cache.js';
import {
  adoptCurrentSessionResources,
  createSessionResourceFileView,
  sessionResourceSource
} from './session-resource-backing.js';
import { isAdoptedCanonicalSession } from './session-authority.js';
import { recordStructuralMetric } from './runtime-test-hooks.js';
import {
  collectCanonicalResourceIds,
  isCanonicalResourceReferenceField
} from './canonical-resource-references.js';

/** @import { CanonicalRenderEnvelope } from './session-request.js' */
/** @import { EncodedFilePayload } from './file-content-cache.js' */
/** @import { RawSessionResourceDescriptor } from './session-resource-backing.js' */
/** @import { BoundFileValue, DrawingInputFiles, LinearSeq } from './drawing/inputs.js' */
/**
 * What a Session binds of a drawing (`webBindings`): a Circular drawing's input
 * fields, a Linear drawing's record rows and comparison plan.
 * @typedef {{
 *   readonly mode: string, readonly files: Readonly<Partial<DrawingInputFiles>>,
 *   readonly linearSeqs: readonly Readonly<LinearSeq>[],
 *   readonly linearComparisonPlan?: { readonly edges?: readonly LinearComparisonEdge[] } | null
 * }} SessionBindingDrawing
 */
/** @import { LinearComparisonEdge } from './linear-comparisons.js' */
/** @import { InputFileRef, InputFileStore } from './history-files.js' */

/**
 * What the allocator reads of a resource's source: a Session resource
 * (`descriptor`, `readBytes`), a File with the payload it was encoded as
 * (`encoded`, `readBytes`), or bytes read from a file.
 * @typedef {object} ResourceSource
 * @property {string} [resourceId]
 * @property {RawSessionResourceDescriptor} [descriptor]
 * @property {EncodedFilePayload} [encoded]
 * @property {Uint8Array} [bytes]
 * @property {() => Promise<Uint8Array>} [readBytes]
 */
/** @typedef {{ name?: unknown, kind?: unknown, type?: unknown, lastModified?: unknown }} ResourceMetadata */
/**
 * The shape `sessionResourceSource` returns: a single backing (`descriptor`,
 * `readBytes`) or a composite of backings (`descriptors`).
 * @typedef {ResourceSource & { descriptors?: (ResourceSource & ResourceMetadata)[] }} SessionResourceSource
 */
/** @typedef {ResourceSource & { size: unknown, bytesPromise: Promise<Uint8Array> | null, identityPromise: Promise<string> | null }} ResourceCandidate */
/** @typedef {{ encoded: Map<unknown, Map<unknown, string>>, identities: Map<string, string>, pending: Map<string, ResourceCandidate> }} SizeBucket */

/** @param {unknown} value */
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

/**
 * @param {{ name?: unknown, type?: unknown, lastModified?: unknown } | null | undefined} file
 * @param {string} resourceId
 */
const bindingForFile = (file, resourceId) => ({
  resourceId,
  name: file?.name === undefined ? 'file' : String(file.name),
  type: String(file?.type || ''),
  lastModified: Number(file?.lastModified) || 0
});

/**
 * @param {{ encoding?: unknown, data?: unknown, size?: unknown }} left
 * @param {{ encoding?: unknown, data?: unknown, size?: unknown }} right
 */
const sameEncodedPayload = (left, right) => (
  left === right || (
    left.encoding === right.encoding
    && left.data === right.data
    && Number(left.size) === Number(right.size)
  )
);

// A reference is replaced by another string id, so the rewritten value has the shape of `value`.
/**
 * @template T
 * @param {T} value
 * @param {Map<string, string>} aliases
 * @returns {T}
 */
const rewriteResourceRefs = (value, aliases) => {
  if (Array.isArray(value)) {
    return /** @type {T} */ (value.map((item) => rewriteResourceRefs(item, aliases)));
  }
  if (!value || typeof value !== 'object') return value;
  return /** @type {T} */ (Object.fromEntries(
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
  ));
};

/**
 * @param {CanonicalRenderEnvelope['webFiles'] | null | undefined} webFiles
 * @param {Map<string, string>} aliases
 * @returns {CanonicalRenderEnvelope['webFiles']}
 */
const rewriteOriginalNameHints = (webFiles, aliases) => {
  const identityAliases = [...aliases].every(([source, target]) => source === target);
  if (identityAliases) {
    const rewritten = { ...(webFiles || {}) };
    delete rewritten.bindings;
    if (Array.isArray(rewritten.linearRecordMetadata)) {
      rewritten.linearRecordMetadata = rewritten.linearRecordMetadata.map((/** @type {unknown} */ entry) => {
        if (!entry || typeof entry !== 'object' || Array.isArray(entry)) return entry;
        const { losatFilename: _losatFilename, ...metadata } = /** @type {Record<string, unknown>} */ (entry); // a non-array object
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
    rewritten.linearRecordMetadata = rewritten.linearRecordMetadata.map((/** @type {unknown} */ entry) => {
      if (!entry || typeof entry !== 'object' || Array.isArray(entry)) return entry;
      const { losatFilename: _losatFilename, ...metadata } = /** @type {Record<string, unknown>} */ (entry); // a non-array object
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
/**
 * @param {CanonicalRenderEnvelope['webFiles']} webFiles
 * @param {CanonicalRenderEnvelope['webFiles'] | null | undefined} otherWebFiles
 * @param {Map<string, string>} aliases
 */
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
 * What Save writes of one drawing: the drawing whose inputs its bindings name
 * (null: no bindings) and its committed request (null: none).
 * @typedef {{ drawing: SessionBindingDrawing | null, committed: CanonicalRenderEnvelope | null }} SessionResourceSet
 */
/**
 * One drawing's request and Web metadata, its resource references rewritten
 * to the Session's table, and the bindings of its own mode's inputs and the
 * file cards (a Session 49 drawing's `webFiles`).
 * @typedef {{
 *   renderRequest: CanonicalRenderEnvelope['renderRequest'] | null,
 *   webFiles: CanonicalRenderEnvelope['webFiles']
 * }} SessionResourceDrawing
 */
/**
 * One resource table for the committed requests and the input bindings of
 * every drawing Save writes: all are allocated through the same
 * content-deduplicating allocator, so equal bytes share one resource and a
 * positional ID that names different bytes in two requests gets a new ID;
 * each request's references are rewritten. The first set with a request is
 * the primary one: an adopted table's resources keep their IDs. `state` holds
 * the store. A Circular drawing binds its input fields, a Linear drawing its
 * records and comparison plan, whose uploaded comparison files the Session
 * also binds; every drawing binds its own five card files (R16-4).
 * @param {Record<string, any>} state
 * @param {readonly SessionResourceSet[]} sets
 * @returns {Promise<{ resources: Record<string, RawSessionResourceDescriptor>, written: SessionResourceDrawing[] }>} `written` in the order of `sets`.
 */
export const buildSessionDrawingResources = async (state, sets) => {
  const primary = sets.find((set) => set.committed !== null) ?? null;
  const committedRequest = primary?.committed ?? null;
  for (const { committed } of sets) {
    if (committed !== null && (typeof committed !== 'object' || !committed.renderRequest || !committed.resources)) {
      throw new Error('A committed canonical render request is required to save a session.');
    }
  }

  /** @type {Record<string, RawSessionResourceDescriptor>} */
  const resources = {};
  /** @type {Map<string, string>} */
  const aliases = new Map();
  const reuseEncodedResources = isAdoptedCanonicalSession(committedRequest);
  // A textual ID is canonical only within its own source table.
  const committedTable = adoptCurrentSessionResources(committedRequest?.resources || {});
  /** @type {Map<unknown, SizeBucket>} */
  const candidatesBySize = new Map();
  /** @type {WeakMap<WeakKey, ResourceCandidate>} */
  const candidatesByDescriptor = new WeakMap();
  /** @type {Set<unknown>} */
  const usedNames = new Set();
  // Resources of an adopted table that nothing names yet: they keep their IDs
  // and names, and are written only when a binding or request uses them.
  /** @type {Map<string, RawSessionResourceDescriptor>} */
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
  /** @type {Promise<unknown>} */
  let allocation = Promise.resolve();

  const nextResourceId = () => {
    /** @type {string} */
    let resourceId;
    do {
      resourceId = `resource-${String(nextResourceNumber++).padStart(4, '0')}`;
    } while (taken(resourceId));
    return resourceId;
  };
  /** @param {ResourceSource} source */
  const candidateFor = (source) => {
    // A source has a descriptor, an encoded payload, or bytes.
    const key = /** @type {WeakKey} */ (source.descriptor || source.encoded || source.bytes);
    if (candidatesByDescriptor.has(key)) {
      return /** @type {ResourceCandidate} */ (candidatesByDescriptor.get(key));
    }
    /** @type {ResourceCandidate} */
    const candidate = {
      ...source,
      size: source.descriptor?.size ?? source.encoded?.size ?? /** @type {Uint8Array} */ (source.bytes).byteLength,
      bytesPromise: null,
      identityPromise: null
    };
    candidatesByDescriptor.set(key, candidate);
    return candidate;
  };
  /** @param {ResourceCandidate} candidate */
  const bytesFor = candidate => {
    // A candidate without bytes is a single backing, which has readBytes.
    candidate.bytesPromise ??= Promise.resolve().then(
      () => candidate.bytes || /** @type {() => Promise<Uint8Array>} */ (candidate.readBytes)()
    );
    return candidate.bytesPromise;
  };
  /** @param {ResourceCandidate} candidate */
  const identityFor = candidate => {
    candidate.identityPromise ??= bytesFor(candidate).then(bytes => {
      // readBytes has already validated any declared checksum, which is a string.
      if (candidate.descriptor?.checksum) {
        return /** @type {string} */ (candidate.descriptor.checksum).trim().toLowerCase().replace(/^sha256:/, '');
      }
      recordStructuralMetric('resourceIdentityHashCount', 1);
      return sha256Hex(bytes);
    });
    return candidate.identityPromise;
  };
  /**
   * @param {string} id
   * @param {ResourceCandidate} candidate
   */
  const register = (id, candidate) => {
    let bucket = candidatesBySize.get(candidate.size);
    if (!bucket) {
      bucket = { encoded: new Map(), identities: new Map(), pending: new Map() };
      candidatesBySize.set(candidate.size, bucket);
    }
    const descriptor = candidate.descriptor || candidate.encoded;
    if (descriptor) {
      let encoded = bucket.encoded.get(descriptor.encoding);
      if (!encoded) bucket.encoded.set(descriptor.encoding, encoded = new Map());
      if (!encoded.has(descriptor.data)) encoded.set(descriptor.data, id);
    }
    bucket.pending.set(id, candidate);
    // `id` is in `resources` or `reserved`.
    usedNames.add((resources[id] || /** @type {RawSessionResourceDescriptor} */ (reserved.get(id))).name);
  };
  /** @param {ResourceCandidate} candidate */
  const equivalentId = async candidate => {
    const bucket = candidatesBySize.get(candidate.size);
    if (!bucket) return null;
    const descriptor = candidate.descriptor || candidate.encoded;
    const encodedId = descriptor && bucket.encoded.get(descriptor.encoding)?.get(descriptor.data);
    if (encodedId) return encodedId;
    const identity = await identityFor(candidate);
    if (bucket.identities.has(identity)) return /** @type {string} */ (bucket.identities.get(identity));
    for (const [id, existing] of bucket.pending) {
      const existingIdentity = await identityFor(existing);
      bucket.pending.delete(id);
      if (!bucket.identities.has(existingIdentity)) bucket.identities.set(existingIdentity, id);
      if (existingIdentity === identity) return id;
    }
    return null;
  };
  /**
   * @param {ResourceSource} source
   * @param {ResourceMetadata} metadata
   * @param {string} [preferredId]
   * @returns {Promise<string>}
   */
  const allocate = (source, metadata, preferredId = '') => {
    // File arrays may resolve concurrently; allocation itself is one ordered
    // transaction so two equal sources cannot publish duplicate payloads.
    allocation = allocation.then(async () => {
      const candidate = candidateFor(source);
      const existing = preferredId ? (resources[preferredId] || reserved.get(preferredId)) : undefined;
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
        // A source with a `resourceId` (the only way a preferred id arises) has a descriptor.
        && descriptor?.name === safeResourceLeaf(/** @type {RawSessionResourceDescriptor} */ (descriptor).name)
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
        // A File Generate already encoded keeps that text; other bytes are encoded here.
        /** @type {{ data: string, size: number }} */
        let encoded;
        if (source.encoded) {
          encoded = source.encoded;
        } else {
          const bytes = await bytesFor(candidate);
          recordStructuralMetric('base64EncodeCount', 1, { resourceName: metadata.name });
          recordStructuralMetric('encodedByteCount', bytes.byteLength, { resourceName: metadata.name });
          encoded = { data: bytesToBase64(bytes), size: bytes.byteLength };
        }
        resources[id] = {
          kind: String(metadata.kind || 'web-file'), name,
          type: String(metadata.type || 'application/octet-stream'),
          size: encoded.size, lastModified: Number(metadata.lastModified) || 0,
          encoding: 'base64', data: encoded.data
        };
      }
      register(id, candidate);
      return id;
    });
    // `allocation` is the promise just chained, which resolves to the allocated id.
    return /** @type {Promise<string>} */ (allocation);
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
      // A view of one resource id is a single backing.
      register(id, candidateFor(/** @type {ResourceSource} */ (
        sessionResourceSource(createSessionResourceFileView(committedTable, id))
      )));
    });
  }
  for (const id of collectCanonicalResourceIds(committedRequest?.renderRequest)) {
    if (!Object.hasOwn(committedResources, id)) {
      throw new Error(`Committed render resource is missing: ${id}.`);
    }
    // A view of one resource id is a single backing, which always has readBytes.
    const source = /** @type {SessionResourceSource & { descriptor: RawSessionResourceDescriptor, readBytes: () => Promise<Uint8Array> }} */ (
      sessionResourceSource(createSessionResourceFileView(committedTable, id))
    );
    if (!reuseEncodedResources) await source.readBytes();
    aliases.set(id, reuseEncodedResources ? id : await allocate(source, source.descriptor));
  }
  /** @type {Map<SessionResourceSet, Map<string, string>>} */
  const setAliases = new Map();
  if (primary) setAliases.set(primary, aliases);
  for (const set of sets) {
    const other = set.committed;
    if (!other || set === primary) continue;
    /** @type {Map<string, string>} */
    const otherAliases = new Map();
    setAliases.set(set, otherAliases);
    const reuseOther = isAdoptedCanonicalSession(other);
    const otherTable = other.resources === committedRequest?.resources
      ? committedTable
      : adoptCurrentSessionResources(other.resources || {});
    for (const id of collectCanonicalResourceIds(other.renderRequest)) {
      if (!Object.hasOwn(other.resources, id)) {
        throw new Error(`Committed render resource is missing: ${id}.`);
      }
      const source = /** @type {SessionResourceSource & { descriptor: RawSessionResourceDescriptor, readBytes: () => Promise<Uint8Array> }} */ (
        sessionResourceSource(createSessionResourceFileView(otherTable, id))
      );
      if (!reuseOther) await source.readBytes();
      otherAliases.set(id, await allocate(source, source.descriptor, reuseOther ? id : ''));
    }
  }

  /**
   * @param {ResourceSource} source
   * @param {ResourceMetadata} metadata
   */
  const bindSource = async (source, metadata) => bindingForFile(
    metadata, await allocate(source, metadata, source.resourceId)
  );
  /** @type {Pick<InputFileStore, 'fileOf'>} */
  const inputFiles = state.inputFiles;
  /**
   * The binding of one input file: a binding names it by ref (`state.inputFiles`).
   * @param {InputFileRef | null | undefined} ref
   */
  const bindFile = async (ref) => {
    const file = inputFiles.fileOf(ref);
    if (!file) return null;
    const source = /** @type {SessionResourceSource | null} */ (getSessionResourceSource(file));
    if (Array.isArray(source?.descriptors)) {
      const components = [];
      for (const component of source.descriptors) components.push(await bindSource(component, component));
      const { resourceId: _resourceId, ...metadata } = bindingForFile(file, '');
      return { kind: 'composite', components, ...metadata };
    }
    if (source?.descriptor) return bindSource(source, file);
    // The payload Generate encoded for this File matches the committed resource
    // by its text, so binding it reads, hashes and decodes nothing.
    const encoded = encodedFilePayloadOf(file);
    if (encoded) return bindSource({ encoded, readBytes: () => readFileBytes(file) }, file);
    return bindSource({ bytes: await readFileBytes(file) }, file);
  };

  /**
   * @param {BoundFileValue} value
   * @returns {Promise<unknown>}
   */
  const bindFileValue = async (value) => {
    if (Array.isArray(value)) {
      return Promise.all(value.map((item) => bindFileValue(item)));
    }
    // `Array.isArray` does not narrow a readonly list out of the union.
    return bindFile(/** @type {InputFileRef | null | undefined} */ (value));
  };

  /** @param {SessionBindingDrawing} drawing Its five card files (R16-4). */
  const cardBindings = async ({ files }) => ({
    d_color: await bindFile(files.d_color),
    t_color: await bindFile(files.t_color),
    blacklist: await bindFile(files.blacklist),
    whitelist: await bindFile(files.whitelist),
    qualifier_priority: await bindFile(files.qualifier_priority)
  });
  /** @param {SessionBindingDrawing} drawing The bindings of its mode's inputs. */
  const circularBindings = async ({ files }) => ({
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
    )
  });
  /** @param {SessionBindingDrawing} drawing */
  const linearBindings = async ({ linearSeqs, linearComparisonPlan }) => ({
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
      (linearComparisonPlan?.edges ?? [])
        .filter((/** @type {LinearComparisonEdge} */ comparison) => comparison?.file)
        .map(async (/** @type {LinearComparisonEdge} */ comparison) => ({
          id: String(comparison?.id || ''),
          file: await bindFile(comparison.file)
        }))
    )
  });
  /** @param {SessionBindingDrawing} drawing */
  const bindingsOf = async (drawing) => {
    // The cards first: a drawing's files are allocated in the order 48 wrote them.
    const cards = await cardBindings(drawing);
    return {
      schema: 2,
      ...(drawing.mode === 'linear' ? await linearBindings(drawing) : await circularBindings(drawing)),
      ...cards
    };
  };

  /** @type {SessionResourceDrawing[]} */
  const written = [];
  for (const set of sets) {
    const setAlias = setAliases.get(set) ?? new Map();
    const identity = [...setAlias].every(([source, target]) => source === target);
    const renderRequest = set.committed
      ? (identity ? set.committed.renderRequest : rewriteResourceRefs(set.committed.renderRequest, setAlias))
      : null;
    const webFiles = rewriteOriginalNameHints(set.committed?.webFiles, setAlias);
    written.push({ renderRequest, webFiles: set.drawing ? { ...webFiles, bindings: await bindingsOf(set.drawing) } : webFiles });
  }
  return { resources, written };
};

/**
 * `buildSessionDrawingResources` of one committed request, the other mode's
 * one (E1) and the bindings of `drawings` merged into one set, as the tests of
 * the allocator read it.
 * @param {Record<string, any>} state
 * @param {CanonicalRenderEnvelope | null} committedRequest
 * @param {readonly SessionBindingDrawing[]} [drawings] `state.drawingList()`
 * @param {CanonicalRenderEnvelope | null} [otherCommittedRequest]
 */
export const buildSessionResources = async (state, committedRequest, drawings = [], otherCommittedRequest = null) => {
  const built = await buildSessionDrawingResources(state, [
    { drawing: null, committed: committedRequest },
    ...(otherCommittedRequest ? [{ drawing: null, committed: otherCommittedRequest }] : []),
    // Without drawings, the bindings of an empty drawing of each mode.
    ...(drawings.length ? drawings : /** @type {SessionBindingDrawing[]} */ ([
      { mode: 'circular', files: {}, linearSeqs: [] }, { mode: 'linear', files: {}, linearSeqs: [] }
    ])).map((drawing) => ({ drawing, committed: null }))
  ]);
  const [top, other] = built.written;
  // One set of bindings: each mode's own fields, and the card files of the
  // first drawing (a test's fixture binds them in one drawing).
  const bindings = Object.assign({}, ...built.written.slice(otherCommittedRequest ? 2 : 1).reverse().map((drawing) => drawing.webFiles.bindings));
  return {
    renderRequest: top.renderRequest,
    ...(otherCommittedRequest ? { otherRenderRequest: other.renderRequest } : {}),
    resources: built.resources,
    webFiles: {
      ...(otherCommittedRequest ? mergeOtherRequestWebFiles(top.webFiles, other.webFiles, new Map()) : top.webFiles),
      bindings
    }
  };
};
