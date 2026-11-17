// @ts-check
import {
  contentFactsOfStream,
  contentTokenOfBytes,
  fileDescriptorOf,
  getSessionResourceSource,
  storedFileBytes
} from './file-content-cache.js';
import {
  adoptCurrentSessionResources,
  createSessionResourceFileView,
  fileBackedDescriptor,
  fileBackedReaderOf,
  renamedResourceDescriptor,
  resourcePayloadKey,
  sessionResourceSource,
  storedResourceBytes
} from './session-resource-backing.js';
import { isAdoptedCanonicalSession } from './session-authority.js';
import { recordStructuralMetric } from './runtime-test-hooks.js';
import { losatRecordGencode, ringGenomeIdLists } from './drawing/inputs.js';
import {
  collectCanonicalResourceIds,
  isCanonicalResourceReferenceField
} from './canonical-resource-references.js';

/** @import { CanonicalRenderEnvelope } from './session-request.js' */
/** @import { FileLike, RawSessionResourceDescriptor, StoredResourceBytes } from './session-resource-backing.js' */
/** @import { BoundFileValue, CircularSource, DrawingInputFiles, LinearSeq } from './drawing/inputs.js' */
/**
 * What a Session binds of a drawing (`webBindings`): a Circular drawing's File
 * cards and input fields, a Linear drawing's record rows and comparison plan.
 * @typedef {{
 *   readonly mode: string, readonly files: Readonly<Partial<DrawingInputFiles>>,
 *   readonly circularSources: readonly Readonly<CircularSource>[],
 *   readonly linearSeqs: readonly Readonly<LinearSeq>[],
 *   readonly linearComparisonPlan?: { readonly edges?: readonly LinearComparisonEdge[] } | null
 * }} SessionBindingDrawing
 */
/** @import { LinearComparisonEdge } from './linear-comparisons.js' */
/** @import { InputFileRef, InputFileStore } from './history-files.js' */

/**
 * What the allocator reads of a resource's source: a Session resource, or the
 * file-backed descriptor of a File (`descriptor`, `readBytes`).
 * @typedef {object} ResourceSource
 * @property {string} [resourceId]
 * @property {RawSessionResourceDescriptor} descriptor
 * @property {() => Promise<Uint8Array>} [readBytes]
 */
/** @typedef {{ name?: unknown, kind?: unknown, type?: unknown, lastModified?: unknown }} ResourceMetadata */
/**
 * The shape `sessionResourceSource` returns: a single backing (`descriptor`, `readBytes`).
 * @typedef {ResourceSource} SessionResourceSource
 */
/**
 * A source the allocator compares: what the Session file stores for it (`stored`, whose `size` it is
 * bucketed by) and the content token of those bytes once asked for.
 * @typedef {ResourceSource & { stored: StoredResourceBytes, size: number, identityPromise: Promise<string> | null }} ResourceCandidate
 */
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
 * @param {RawSessionResourceDescriptor} left
 * @param {RawSessionResourceDescriptor} right
 */
const sameEncodedPayload = (left, right) => (
  left === right || (
    left.encoding === right.encoding
    && resourcePayloadKey(left) === resourcePayloadKey(right)
    && Number(left.size) === Number(right.size)
  )
);

/**
 * A committed resource of a Save and the candidate that reads its bytes.
 * @typedef {{ id: string, descriptor: RawSessionResourceDescriptor, candidate: ResourceCandidate }} CommittedPayload
 */
/**
 * Reads an encoded source's bytes whole, which checks their size and any
 * declared checksum (`materializeBytes`): its base64 text holds them whole
 * already (W-14). A file-backed source is never read whole: the Session
 * writer checks its size as it streams.
 * @param {ResourceSource} source
 */
const checkEncodedPayload = async (source) => {
  if (!fileBackedReaderOf(source.descriptor) && source.readBytes) await source.readBytes();
};

/**
 * A File no Generate described, as a Save source: a file-backed descriptor of
 * what the Session file stores for it (`storedFileBytes`), its copy streamed
 * when the Session is written, with the size and token its first read found;
 * nothing is read here.
 * @param {FileLike & { name?: unknown, type?: unknown, lastModified?: unknown }} file
 * @returns {Promise<ResourceSource>}
 */
const storedFileSource = async (file) => {
  const stored = await storedFileBytes(file);
  const readStored = async () => new Uint8Array(await new Response(await stored.open()).arrayBuffer());
  const descriptor = fileBackedDescriptor({
    kind: 'web-file', name: String(file.name ?? 'file'), type: String(file.type || 'application/octet-stream'),
    size: stored.size, lastModified: Number(file.lastModified) || 0
  }, readStored, stored);
  return { descriptor, readBytes: readStored };
};

/**
 * The Save source of a File no Generate described (`storedFileSource`), made
 * once per Save.
 * @param {Map<object, Promise<ResourceSource>>} made The sources this Save made, by File.
 * @param {FileLike & { name?: unknown, type?: unknown, lastModified?: unknown }} file
 * @returns {Promise<ResourceSource>}
 */
const storedFileSourceOnce = (made, file) => {
  const known = made.get(file);
  if (known) return known;
  const source = storedFileSource(file);
  made.set(file, source);
  return source;
};

/**
 * The ID of the first of `earlier` whose bytes equal `entry`'s, or null (D-228).
 * Two encoded texts compare as text, with no read or hash; when either is
 * file-backed, the content identities of their bytes compare.
 * @param {readonly CommittedPayload[]} earlier
 * @param {CommittedPayload} entry
 * @param {(candidate: ResourceCandidate) => Promise<string>} identityOf
 * @returns {Promise<string | null>}
 */
const firstEqualPayload = async (earlier, entry, identityOf) => {
  for (const other of earlier) {
    const equal = fileBackedReaderOf(other.descriptor) || fileBackedReaderOf(entry.descriptor)
      ? await identityOf(other.candidate) === await identityOf(entry.candidate)
      : other.descriptor.encoding === entry.descriptor.encoding && other.descriptor.data === entry.descriptor.data;
    if (equal) return other.id;
  }
  return null;
};

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
    // One list of resource IDs per ring (`ringGenomeIdLists`).
    rewritten[field] = ringGenomeIdLists(webFiles[field]).map((ids) => ids.flatMap((resourceId) => {
      const alias = aliases.get(resourceId);
      return alias ? [alias] : [];
    }));
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
 * file cards (a Session 51 drawing's `webFiles`).
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
 * A File's bytes stay file-backed (`fileBackedDescriptor`): the Session
 * writer streams them into the Session file (`writeSessionDocumentFile` in config.js).
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
    const key = /** @type {WeakKey} */ (source.descriptor);
    const known = candidatesByDescriptor.get(key);
    if (known) return known;
    const stored = storedResourceBytes(source.descriptor);
    /** @type {ResourceCandidate} */
    const candidate = { ...source, stored, size: stored.size, identityPromise: null };
    candidatesByDescriptor.set(key, candidate);
    return candidate;
  };
  /**
   * The content token of what the Session file stores for a candidate, so
   * equal stored bytes are written once: known from a File's first read or a
   * Load's copy, else hashed from a stream of them, one piece at a time; an
   * encoded source's bytes, whole in its text already, are hashed whole.
   * @param {ResourceCandidate} candidate
   */
  const identityFor = candidate => {
    candidate.identityPromise ??= (async () => {
      const { descriptor, readBytes, stored } = candidate;
      if (!fileBackedReaderOf(descriptor) && readBytes) {
        const bytes = await readBytes();
        recordStructuralMetric('resourceIdentityHashCount', 1);
        return contentTokenOfBytes(bytes);
      }
      if (stored.token) return stored.token();
      recordStructuralMetric('resourceIdentityHashCount', 1);
      return (await contentFactsOfStream(await stored.open())).token;
    })();
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
    const descriptor = candidate.descriptor;
    if (descriptor) {
      let encoded = bucket.encoded.get(descriptor.encoding);
      if (!encoded) bucket.encoded.set(descriptor.encoding, encoded = new Map());
      const key = resourcePayloadKey(descriptor);
      if (!encoded.has(key)) encoded.set(key, id);
    }
    bucket.pending.set(id, candidate);
    // `id` is in `resources` or `reserved`.
    usedNames.add((resources[id] || /** @type {RawSessionResourceDescriptor} */ (reserved.get(id))).name);
  };
  /** @param {ResourceCandidate} candidate */
  const equivalentId = async candidate => {
    const bucket = candidatesBySize.get(candidate.size);
    if (!bucket) return null;
    const descriptor = candidate.descriptor;
    const encodedId = descriptor && bucket.encoded.get(descriptor.encoding)?.get(resourcePayloadKey(descriptor));
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
      if (source.descriptor.checksum && source.descriptor !== existing
        && source.descriptor.checksum !== existing?.checksum) await checkEncodedPayload(source);
      if (existing && sameEncodedPayload(existing, source.descriptor)) {
        return use(preferredId);
      }
      const equivalent = await equivalentId(candidate);
      if (equivalent) return use(equivalent);
      // Preserve the existing collision integrity check through the lazy backing.
      if (existing) await checkEncodedPayload(source);
      const descriptor = source.descriptor;
      const safePreferred = /^[a-z][a-z0-9]*(?:-[a-z0-9]+)*$/.test(preferredId)
        && !taken(preferredId)
        && descriptor.name === safeResourceLeaf(descriptor.name)
        && !usedNames.has(descriptor.name);
      let id = safePreferred ? preferredId : nextResourceId();
      let name = safePreferred ? descriptor.name : `${id}-${safeResourceLeaf(metadata.name)}`;
      while (usedNames.has(name)) {
        id = nextResourceId();
        name = `${id}-${safeResourceLeaf(metadata.name)}`;
      }
      resources[id] = safePreferred ? descriptor : renamedResourceDescriptor(descriptor, String(name));
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
    // so a resource that nothing names any more is not written. A resource
    // whose bytes, kind and original name equal an earlier one's is written
    // once (D-220, D-228): its id becomes an alias of the earlier id. Bytes are
    // compared only when kind, name and size match (`firstEqualPayload`).
    const webFiles = committedRequest?.webFiles || {};
    const named = collectCanonicalResourceIds(committedRequest?.renderRequest);
    collectCanonicalResourceIds(webFiles, named);
    ['conservationLosatFastaSources', 'conservationSequenceSources'].forEach((field) => {
      ringGenomeIdLists(webFiles[field]).flat().forEach((id) => named.add(id));
    });
    const hints = webFiles.resourceOriginalNames || {};
    /** @type {Map<string, CommittedPayload[]>} */
    const firstIds = new Map();
    for (const [id, descriptor] of Object.entries(committedResources)) {
      // A missing name equals only a missing name.
      const key = JSON.stringify([descriptor.kind, Object.hasOwn(hints, id) ? hints[id] : null, descriptor.size]);
      let earlier = firstIds.get(key);
      if (!earlier) firstIds.set(key, earlier = []);
      // A view of one resource id is a single backing.
      const entry = { id, descriptor, candidate: candidateFor(/** @type {ResourceSource} */ (
        sessionResourceSource(createSessionResourceFileView(committedTable, id))
      )) };
      const equal = await firstEqualPayload(earlier, entry, identityFor);
      if (equal) {
        aliases.set(id, named.has(id) ? use(equal) : equal);
        continue;
      }
      earlier.push(entry);
      if (named.has(id)) resources[id] = descriptor;
      else reserved.set(id, descriptor);
      aliases.set(id, id);
      register(id, entry.candidate);
    }
  }
  for (const id of collectCanonicalResourceIds(committedRequest?.renderRequest)) {
    if (!Object.hasOwn(committedResources, id)) {
      throw new Error(`Committed render resource is missing: ${id}.`);
    }
    // An adopted request's ids were kept or aliased above.
    if (reuseEncodedResources) continue;
    // A view of one resource id is a single backing, which always has readBytes.
    const source = /** @type {SessionResourceSource & { descriptor: RawSessionResourceDescriptor, readBytes: () => Promise<Uint8Array> }} */ (
      sessionResourceSource(createSessionResourceFileView(committedTable, id))
    );
    await checkEncodedPayload(source);
    aliases.set(id, await allocate(source, source.descriptor));
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
      if (!reuseOther) await checkEncodedPayload(source);
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
  /** @type {Map<object, Promise<ResourceSource>>} */
  const storedFileSources = new Map();
  /**
   * The binding of one input file: a binding names it by ref (`state.inputFiles`).
   * @param {InputFileRef | null | undefined} ref
   */
  const bindFile = async (ref) => {
    const file = inputFiles.fileOf(ref);
    if (!file) return null;
    const source = /** @type {SessionResourceSource | null} */ (getSessionResourceSource(file));
    if (source?.descriptor) return bindSource(source, file);
    // The descriptor Generate made of this File names the committed resource
    // by its reader, so binding it reads, hashes and decodes nothing.
    const described = fileDescriptorOf(file);
    const readFile = fileBackedReaderOf(described);
    if (described && readFile) {
      const descriptor = fileBackedDescriptor({
        kind: 'web-file', name: String(file.name ?? 'file'), type: String(file.type || 'application/octet-stream'),
        size: Number(described.size), lastModified: Number(file.lastModified) || 0
      }, readFile);
      return bindSource({ descriptor, readBytes: readFile }, file);
    }
    return bindSource(await storedFileSourceOnce(storedFileSources, file), file);
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
  /** @param {SessionBindingDrawing} drawing The bindings of its mode's inputs: one row per File card. */
  const circularBindings = async ({ files, circularSources }) => ({
    circularSources: await Promise.all(circularSources.map(async (card) => ({
      uid: card.uid,
      gb: await bindFile(card.gb),
      gff: await bindFile(card.gff),
      fasta: await bindFile(card.fasta)
    }))),
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
      losat_gencode: losatRecordGencode(sequence),
      definition: String(sequence?.definition || ''),
      record_subtitle: String(sequence?.record_subtitle || ''),
      file_definition: String(sequence?.file_definition || ''),
      file_subtitle: String(sequence?.file_subtitle || ''),
      inferred_definition: String(sequence?.inferred_definition || ''),
      region_record_id: String(sequence?.region_record_id || ''),
      region_start: sequence?.region_start ?? null, region_end: sequence?.region_end ?? null,
      region_reverse: Boolean(sequence?.region_reverse),
      windows: Array.isArray(sequence?.windows) ? sequence.windows.map((window) => ({ ...window })) : []
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
      { mode: 'circular', files: {}, linearSeqs: [], circularSources: [] },
      { mode: 'linear', files: {}, linearSeqs: [], circularSources: [] }
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
