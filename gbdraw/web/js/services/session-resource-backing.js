// @ts-check
import {
  base64ToBytesInTasks,
  bytesToText,
  checkGzipExpandedSize,
  expandGzipBytes,
  isGzipBytes,
  sha256Hex,
  textToBytes
} from './byte-utils.js';
import {
  decodeDepthText,
  isEncodedDepthFileEntry
} from './depth-file-codec.js';
import { CIRCULAR_SOURCE_FILE_SLOTS, SESSION_DRAWINGS_LAYOUT } from '../contracts/schema.generated.js';
import { recordStructuralMetric, testByteCap } from './runtime-test-hooks.js';
import { diagnosticError } from '../utils/error-normalization.js';
import { MAX_EXPANDED_INPUT_BYTES } from '../contracts/schema.generated.js';

/** @import { DepthFilePayload } from './depth-file-codec.js' */
/** @import { ResourceContentToken } from './file-content-cache.js' */

/**
 * A Session resource as the Session file stores it, before any check: every
 * field is unvalidated.
 * @typedef {{
 *   kind?: unknown, name?: unknown, type?: unknown, encoding?: unknown, data?: unknown,
 *   size?: unknown, lastModified?: unknown, checksum?: unknown
 * }} RawSessionResourceDescriptor
 */
/**
 * A Session resource after `validateDescriptor` that holds its payload: a
 * base64 text, or an encoded depth table (`encoding` names which `data` holds).
 * @typedef {{
 *   kind?: unknown, name?: unknown, type?: unknown, encoding: string, data: string | DepthFilePayload,
 *   size: number, lastModified?: number, checksum?: unknown
 * }} EncodedSessionResourceDescriptor
 */
/**
 * A Session resource after `validateDescriptor` whose reader gives its bytes
 * (`fileBackedDescriptor`): a File Generate described. It never leaves the page.
 * @typedef {{
 *   kind?: unknown, name?: unknown, type?: unknown, encoding?: undefined, data?: undefined,
 *   size: number, lastModified?: number, checksum?: unknown
 * }} FileBackedSessionResourceDescriptor
 */
/** @typedef {EncodedSessionResourceDescriptor | FileBackedSessionResourceDescriptor} SessionResourceDescriptor */
/**
 * The stand-in for a `File` that a Session resource backs.
 * @typedef {Readonly<{ name: string, type: string, size: number, lastModified: number }>} SessionResourceFileView
 */
/**
 * A `File` or a Session resource view: the bytes an input file names. An
 * `InputFileRef` is not one (`fileId` is never set), so a byte reader that is
 * handed a ref fails type checking; resolve it with `state.inputFiles.fileOf`.
 * @typedef {{
 *   name?: string, type?: string, size?: number, lastModified?: number,
 *   arrayBuffer?: () => Promise<ArrayBuffer>, text?: () => Promise<string>,
 *   fileId?: never
 * }} FileLike
 */
/**
 * @typedef {object} BackingFields
 * @property {string} resourceId
 * @property {string} name
 * @property {string} type
 * @property {number} size
 * @property {number} lastModified
 * @property {Promise<Uint8Array> | null} bytesPromise The stored bytes (the size the descriptor counts).
 * @property {Promise<Uint8Array> | null} [plainPromise] The plain bytes, kept unless the resource is a gzip input (`materializePlainBytes`).
 * @property {Promise<string> | null} textPromise The text, kept unless the resource is a gzip input (`materializeText`).
 * @property {null} filePromise
 * @property {boolean} dirty
 */
/**
 * The backing of one resource: its descriptor and the content it decoded.
 * @typedef {BackingFields & {
 *   descriptor: SessionResourceDescriptor, kind: string, encoding: string, checksum: string | null
 * }} LeafBacking
 */
/** @type {WeakMap<SessionResourceTable, Map<string, LeafBacking>>} */
const tableBackings = new WeakMap();
/** @type {WeakMap<WeakKey, LeafBacking>} */
const fileViewBackings = new WeakMap();

/**
 * The reader of each file-backed descriptor: a descriptor of the page, never
 * of a Session file, whose bytes a File (or a Session resource view) holds. It
 * has no `encoding` and no `data`; a Session writer encodes it
 * (`writeSessionDocumentFile`), so Generate keeps no base64 text of an input
 * (OV-2612). One reader per File: it is also the payload's key.
 * @typedef {() => Promise<Uint8Array>} FileBackedReader
 * @type {WeakMap<object, FileBackedReader>}
 */
const fileBackedReaders = new WeakMap();
/**
 * The opener of a stream of the bytes a Session file stores for a resource.
 * @typedef {() => Promise<ReadableStream<Uint8Array>>} FileBackedStreamOpener
 */
/**
 * What a Session file stores for a resource (its entry `resources/<id>`):
 * the size of the stored bytes, a stream of them, their content token
 * (`contentTokenOfBytes` in file-content-cache.js) when it is known without a
 * read, else null, and a Blob of them when they rest in a file the page does
 * not hold (an input's copy, PD-OI-121, the File a Generate described, or
 * the entry of the Session file a Load read when no copy could be written),
 * else null: Generate stages that Blob, and the page keeps no read of them
 * (OV-2753, OV-2754). A gzip input is stored as it is (receipt
 * session.file-container), so a File chosen as gzip stores its gzip bytes,
 * which the reader of its file-backed descriptor reads too and Generate sends.
 * @typedef {{
 *   size: number, open: FileBackedStreamOpener, token: (() => Promise<ResourceContentToken>) | null,
 *   blob: (() => Promise<Blob>) | null
 * }} StoredResourceBytes
 */
/**
 * The stored bytes of the resources a file-backed reader reads, where they
 * can be streamed (an input's copy, PD-OI-121, or an entry of the Session file
 * a Load read): a Session Save streams them and never holds them whole. Keyed
 * by the reader, so every descriptor that shares the reader shares them.
 * @type {WeakMap<FileBackedReader, StoredResourceBytes>}
 */
const fileBackedStored = new WeakMap();

/**
 * A file-backed descriptor of these fields, whose bytes (the `size` it
 * counts) `readBytes` reads, and which a Session file stores as `stored`
 * says; with no `stored`, it stores the bytes `readBytes` reads.
 * @template {{ kind?: unknown, name: string, type: string, size: number, lastModified: number }} F
 * @param {F} fields
 * @param {FileBackedReader} readBytes
 * @param {StoredResourceBytes | null} [stored]
 * @returns {F}
 */
export const fileBackedDescriptor = (fields, readBytes, stored = null) => {
  const descriptor = { ...fields };
  fileBackedReaders.set(descriptor, readBytes);
  if (stored) fileBackedStored.set(readBytes, stored);
  return descriptor;
};

/**
 * The reader of a file-backed descriptor, or null for any other value.
 * @param {unknown} descriptor
 * @returns {FileBackedReader | null}
 */
export const fileBackedReaderOf = (descriptor) => (
  descriptor && typeof descriptor === 'object' ? fileBackedReaders.get(descriptor) ?? null : null
);

/**
 * What a Session file stores for a checked Session resource descriptor, the
 * one owner of a saved entry's bytes and size: a file-backed descriptor's
 * stored bytes (`fileBackedDescriptor`), else the bytes its reader reads; an
 * encoded depth table's text; a base64 text's bytes (a resource of an older
 * JSON Session, decoded when the entry is written).
 * @param {RawSessionResourceDescriptor} descriptor
 * @returns {StoredResourceBytes}
 */
export const storedResourceBytes = (descriptor) => {
  const readFile = fileBackedReaderOf(descriptor);
  const stored = readFile ? fileBackedStored.get(readFile) : undefined;
  if (stored) return stored;
  /** @param {() => Promise<Uint8Array>} read */
  const whole = (read) => async () => new Blob([/** @type {Uint8Array<ArrayBuffer>} */ (await read())]).stream();
  /** @type {() => Promise<Uint8Array>} */
  let read;
  if (readFile) read = readFile;
  else if (isEncodedDepthFileEntry(/** @type {DepthEntryProbe} */ (descriptor))) {
    read = async () => textToBytes(decodeDepthText(/** @type {DepthFilePayload} */ (descriptor.data)));
  } else read = () => base64ToBytesInTasks(/** @type {string} */ (descriptor.data));
  return { size: /** @type {number} */ (descriptor.size), open: whole(read), token: null, blob: null };
};

/**
 * A stream of the bytes a Session file stores for a checked descriptor (`storedResourceBytes`).
 * @param {RawSessionResourceDescriptor} descriptor
 * @returns {Promise<ReadableStream<Uint8Array>>}
 */
export const sessionResourceStream = (descriptor) => storedResourceBytes(descriptor).open();

/**
 * What names a descriptor's payload until its bytes change: the reader of a
 * file-backed descriptor, else its `data` (base64 text, or a Depth payload).
 * @param {RawSessionResourceDescriptor} descriptor
 * @returns {unknown}
 */
export const resourcePayloadKey = (descriptor) => fileBackedReaderOf(descriptor) ?? descriptor.data;

/**
 * The same descriptor under another name: the one owner of a rename, so a
 * file-backed descriptor keeps its reader (a spread copy would lose it).
 * @template {RawSessionResourceDescriptor} D
 * @param {D} descriptor
 * @param {string} name
 * @returns {D}
 */
export const renamedResourceDescriptor = (descriptor, name) => {
  const renamed = { ...descriptor, name };
  const readFile = fileBackedReaderOf(descriptor);
  if (readFile) fileBackedReaders.set(renamed, readFile);
  return renamed;
};

const CHECKSUM_PATTERN = /^[0-9a-f]{64}$/;
/** @param {unknown} value A checked checksum: absent, or a SHA-256 digest text. */
const normalizedChecksum = (value) => (typeof value === 'string' && value !== ''
  ? value.trim().toLowerCase().replace(/^sha256:/, '') : null);

/**
 * One thing a stored resource descriptor holds that the app cannot read. A
 * fault of `checksum` or `lastModified` names a field the app recomputes; a
 * fault of no field (`null`) leaves the resource unreadable.
 * @typedef {{ field: 'checksum' | 'lastModified' | null }} ResourceDescriptorFault
 */
/**
 * The one check of a stored resource descriptor: Load drops what it names
 * (`readSessionResources`), the runtime refuses it (`validateDescriptor`).
 * @param {RawSessionResourceDescriptor | null | undefined} descriptor Unvalidated.
 * @returns {ResourceDescriptorFault[]}
 */
const resourceDescriptorFaults = (descriptor) => {
  if (!descriptor || typeof descriptor !== 'object' || Array.isArray(descriptor)) {
    return [{ field: null }];
  }
  /** @type {ResourceDescriptorFault[]} */
  const faults = [];
  const supportedPayload = fileBackedReaders.has(descriptor) || (descriptor.encoding === 'base64'
    ? typeof descriptor.data === 'string'
    // A `data` that is not an object (a string, or nothing) has no `schema`, so the guard reads false.
    : isEncodedDepthFileEntry(/** @type {{ encoding?: unknown, data?: { schema?: unknown } | null }} */ (descriptor)));
  if (!supportedPayload) faults.push({ field: null });
  if (typeof descriptor.size !== 'number' || !Number.isSafeInteger(descriptor.size) || descriptor.size < 0) {
    faults.push({ field: null });
  }
  if (
    descriptor.lastModified !== undefined
    && (typeof descriptor.lastModified !== 'number' || !Number.isFinite(descriptor.lastModified) || descriptor.lastModified < 0)
  ) {
    faults.push({ field: 'lastModified' });
  }
  const { checksum } = descriptor;
  if (checksum !== null && checksum !== undefined && checksum !== '') {
    if (typeof checksum !== 'string' || !CHECKSUM_PATTERN.test(normalizedChecksum(checksum) ?? '')) faults.push({ field: 'checksum' });
  }
  return faults;
};

/**
 * @param {string} resourceId
 * @param {RawSessionResourceDescriptor | null | undefined} descriptor Unvalidated; the first check proves an object.
 * @returns {LeafBacking}
 */
const validateDescriptor = (resourceId, descriptor) => {
  // R6: the typed Session-fields error; Load lists the fault's place before (`readSessionResources`).
  if (resourceDescriptorFaults(descriptor).length) throw diagnosticError('INPUT_INVALID', { field: 'schema', reason: 'FIELDS' });
  // `resourceDescriptorFaults` found none: an object with a supported payload and size.
  const checked = /** @type {SessionResourceDescriptor} */ (descriptor);
  return {
    descriptor: checked,
    resourceId,
    kind: String(checked.kind || 'web-file'),
    name: String(checked.name || 'file'),
    type: String(checked.type || 'application/octet-stream'),
    encoding: String(checked.encoding || ''),
    size: checked.size,
    lastModified: Number(checked.lastModified) || 0,
    checksum: normalizedChecksum(checked.checksum),
    bytesPromise: null,
    textPromise: null,
    filePromise: null,
    dirty: false
  };
};

/**
 * The id a resource table keys a stored resource by, or null when it is
 * empty or repeats an id `seen` holds (`seen` collects the ids).
 * @param {string} rawResourceId
 * @param {Set<string>} seen
 */
const tableResourceId = (rawResourceId, seen) => {
  const resourceId = String(rawResourceId || '').trim();
  if (!resourceId || seen.has(resourceId)) return null;
  seen.add(resourceId);
  return resourceId;
};

/**
 * Load: the stored resources the app can read (A8). A resource whose
 * descriptor cannot be read, or whose id is empty or repeats one, is dropped
 * and named in `dropped`; a `checksum` or `lastModified` the app recomputes is
 * dropped from its descriptor with no name. The table is copied where it
 * changed; the stored one is never changed.
 * @template T
 * @param {T} resources The Session's `resources`, unvalidated.
 * @returns {{ value: T, dropped: string[] }}
 */
export const readSessionResources = (resources) => {
  if (!resources || typeof resources !== 'object' || Array.isArray(resources)) return { value: resources, dropped: [] };
  /** @type {string[]} */
  const dropped = [];
  /** @type {Set<string>} */
  const seen = new Set();
  let changed = false;
  const kept = Object.entries(resources).flatMap(([rawResourceId, descriptor]) => {
    const resourceId = tableResourceId(rawResourceId, seen);
    const faults = resourceId === null ? [] : resourceDescriptorFaults(descriptor);
    if (resourceId === null || faults.some(({ field }) => field === null)) {
      dropped.push(rawResourceId);
      changed = true;
      return [];
    }
    if (!faults.length) return [[rawResourceId, descriptor]];
    changed = true;
    const fields = new Set(faults.map(({ field }) => field));
    return [[rawResourceId, Object.fromEntries(Object.entries(descriptor).filter(([field]) => !fields.has(/** @type {'checksum'} */ (field))))]];
  });
  return { value: changed ? /** @type {T} */ (Object.fromEntries(kept)) : resources, dropped };
};

/**
 * A file binding after `validateWebFileBindings`: an object with no own `kind`
 * or `components` that names a resource (schema 2 also proves the metadata).
 * @typedef {Readonly<{ resourceId: unknown, kind?: undefined, components?: undefined, [field: string]: unknown }>} WebFileLeafBinding
 */
/**
 * The `c_gb` of a Session 41-50 that bound two or more GenBank files as one:
 * Load makes one File card of each component (`circularSourcesOfModeScoped`).
 * @typedef {Readonly<{
 *   kind: 'composite', name: string, type: string, lastModified: number,
 *   components: Array<{ resourceId: string, name: string, type: string, lastModified: number }>
 * }>} WebFileCompositeBinding
 */
/**
 * One Circular File card of a Session 51 drawing (`circularSources` row,
 * Python's `CircularSourceCard`): a uid, distinct among the drawing's cards,
 * and each file slot (`CIRCULAR_SOURCE_FILE_SLOTS`) a leaf or null.
 * @typedef {Readonly<{ uid: string }
 *   & Record<(typeof CIRCULAR_SOURCE_FILE_SLOTS)[number], WebFileLeafBinding | null>>} WebFileCircularSourceBinding
 */
/**
 * A row of `linearSeqs` or of a comparison list: only its file slots are checked.
 * @typedef {Readonly<Record<string, unknown>> | null | undefined} WebFileBindingRow
 */
/**
 * `webFiles.bindings` after `validateWebFileBindings`. Each file slot holds null,
 * a leaf, or a list of them (a JSON list has no `kind`); only the `c_gb` of a
 * Session 41-50 may hold a composite. A Session 51 Circular drawing holds its
 * File cards in `circularSources`.
 * @typedef {Readonly<{
 *   schema: 1 | 2,
 *   c_gb?: WebFileCompositeBinding | WebFileLeafBinding | (ReadonlyArray<unknown> & { readonly kind?: undefined }) | null,
 *   circularSources?: ReadonlyArray<WebFileCircularSourceBinding>,
 *   [field: string]: unknown
 * }>} WebFileBindings
 */

// Binding grammar belongs beside the File backing; envelope/version and
// comparison-plan authority remain with their existing owners.
/**
 * One thing a Web input inventory holds that the app cannot read. `at` is
 * null for the inventory itself (no object, or a schema no writer wrote), or
 * the keys of one binding under `webFiles.bindings` (`['c_gb']`,
 * `['c_conservation_blasts', 1]`, `['c_depth', 0, 0]`, `['linearSeqs', 0, 'gb']`),
 * or of a field no writer wrote (`field`). `webFileBindingPath` writes it as text.
 * @typedef {{ at: ReadonlyArray<string | number> | null, field?: boolean }} WebFileBindingFault
 */
/**
 * The place of a binding as a row names it: `c_depth[0][0]`, `linearSeqs[0].gb`.
 * @param {ReadonlyArray<string | number>} at
 * @returns {string}
 */
export const webFileBindingPath = (at) => at.map((key, index) => (
  typeof key === 'number' ? `[${key}]` : `${index ? '.' : ''}${key}`)).join('');
// The fields of a Circular File card (`circularSources` row): its uid and its file slots.
const CIRCULAR_SOURCE_FIELDS = Object.freeze(['uid', ...CIRCULAR_SOURCE_FILE_SLOTS]);
/**
 * The binding lists whose index is a place, not an order: a Depth track
 * (`c_depth`, a Linear row's `depth`) and a Circular comparison series's
 * BLAST table and sequence (`c_conservation_blasts`, `c_conservation_fastas`,
 * `c_conservation_sequence_sources`, read by the series' `sourceIndex`).
 * Load leaves a dropped binding's place null, so no later file changes place:
 * the readers of a sequence or Depth list take it as a slot with no file, and
 * Load leaves out the series of an empty table place with its files (OV-893).
 * A binding of any other list is removed from it.
 * @type {ReadonlySet<string>}
 */
export const POSITIONAL_BINDING_LISTS = Object.freeze(new Set([
  'c_depth', 'depth', 'c_conservation_blasts', 'c_conservation_fastas', 'c_conservation_sequence_sources'
]));
/**
 * The input slot a schema-2 inventory holds (`WEB_BINDING_INPUT_SLOT`,
 * `gbdraw/session_io.py`): `c_gb` for a Session 41-50, which binds both
 * modes' inputs; a Session 51 drawing's own mode's (its one-drawing view
 * names the mode in `ui.mode`): `circularSources` or `linearSeqs`.
 * @typedef {'c_gb' | 'circularSources' | 'linearSeqs'} WebBindingInputSlot
 * @param {boolean} drawingView The inventory is a Session 51 drawing's.
 * @param {unknown} mode The view's `ui.mode`.
 * @returns {WebBindingInputSlot}
 */
export const webBindingInputSlot = (drawingView, mode) => (drawingView
  ? /** @type {WebBindingInputSlot} */ (SESSION_DRAWINGS_LAYOUT.bindingInputSlots[mode === 'linear' ? 'linear' : 'circular'])
  : 'c_gb');
/**
 * The one check of a Web input inventory: Load drops what it names
 * (`readWebFileBindings` in session-authority.js), Save and the runtime refuse
 * it (`validateWebFileBindings`).
 * @param {{ bindings?: unknown } | null | undefined} webFiles
 * @param {{ readonly [resourceId: string]: RawSessionResourceDescriptor | null | undefined } | null} [resources]
 *   The Session's `resources`, unvalidated: read by the id of each binding, and each entry read is checked.
 * @param {WebBindingInputSlot} [inputSlot] The input slot schema 2 holds: `c_gb` for a Session 41-50, the
 *   drawing's mode's for a Session 51 drawing's one-drawing view (`gbdraw/session_io.py`).
 * @returns {WebFileBindingFault[]}
 */
export const webFileBindingFaults = (webFiles, resources = {}, inputSlot = 'c_gb') => {
  if (!webFiles || !Object.hasOwn(webFiles, 'bindings')) return [];
  const bindings = webFiles.bindings;
  /**
   * @param {unknown} value
   * @returns {value is Record<string, unknown>}
   */
  const object = value => value !== null && typeof value === 'object' && !Array.isArray(value);
  // The inventory: an object of a schema a writer wrote; schema 2 holds its input slot.
  if (!object(bindings) || !(/** @type {unknown[]} */ ([1, 2])).includes(bindings.schema)) return [{ at: null }];
  const current = bindings.schema === 2;
  if (current && !Object.hasOwn(bindings, inputSlot)) return [{ at: null }];
  /** @param {Record<string, unknown>} value @returns {boolean} Its file metadata cannot be read. */
  const badMetadata = value => typeof value.name !== 'string' || typeof value.type !== 'string'
    || typeof value.lastModified !== 'number' || !Number.isFinite(value.lastModified) || value.lastModified < 0;
  /** @param {unknown} value @returns {boolean} The binding of one file cannot be read. */
  const badLeaf = value => {
    // A file binding is an object that is no composite and names a canonical resourceId.
    if (!object(value) || Object.hasOwn(value, 'kind') || Object.hasOwn(value, 'components')) return true;
    const id = current ? value.resourceId : String(value.resourceId || '').trim();
    if (typeof id !== 'string' || !id || id !== id.trim()) return true;
    // A `resources` that is null holds no resource, like any other value without this id (OV-636).
    if (resources === null || !Object.hasOwn(resources, id)) return true;
    if (!current) return false;
    return Object.keys(value).some(key => !['resourceId', 'name', 'type', 'lastModified'].includes(key))
      || badMetadata(value) || resourceDescriptorFaults(resources[id]).length > 0;
  };
  /** @param {unknown} binding @param {boolean} compositeAllowed @returns {boolean} */
  const badSingle = (binding, compositeAllowed) => {
    if (!(object(binding) && Object.hasOwn(binding, 'kind'))) return badLeaf(binding);
    // A composite is the schema 2 `c_gb` of a Session 41-50 only, with only its fields and two or more base64 components.
    if (!current || !compositeAllowed || binding.kind !== 'composite') return true;
    if (Object.keys(binding).some(key => !['kind', 'components', 'name', 'type', 'lastModified'].includes(key))) return true;
    const components = binding.components;
    if (!Array.isArray(components) || components.length < 2 || badMetadata(binding) || components.some(badLeaf)) return true;
    // `badLeaf` proved that `resources` holds each component's id and that its entry is an object.
    return components.some((/** @type {{ resourceId: string }} */ component) => (
      /** @type {Record<string, RawSessionResourceDescriptor>} */ (resources)[component.resourceId].encoding !== 'base64'));
  };
  /** @type {WebFileBindingFault[]} */
  const faults = [];
  /** @param {unknown} binding @param {ReadonlyArray<string | number>} at @param {boolean} [compositeAllowed] */
  const value = (binding, at, compositeAllowed = false) => {
    if (binding === null || binding === undefined) return;
    if (Array.isArray(binding)) {
      binding.forEach((item, index) => value(item, [...at, index]));
      return;
    }
    if (badSingle(binding, compositeAllowed)) faults.push({ at });
  };
  // A Session 41-50 binds the Circular input in c_gb, c_gff and c_fasta; a
  // Session 51 drawing in its File cards (`circularSources`).
  const modeScoped = inputSlot === 'c_gb';
  const slots = [
    ...(modeScoped ? SESSION_DRAWINGS_LAYOUT.modeScopedCircularInputSlots : []), 'c_depth', 'c_conservation_blasts',
    'c_conservation_fastas', 'c_conservation_sequence_sources', 'd_color',
    't_color', 'blacklist', 'whitelist', 'qualifier_priority'
  ];
  if (current) {
    Object.keys(bindings).filter(key => ![
      'schema', ...slots, 'c_conservation_blasts_source', 'linearSeqs', 'linearComparisons',
      ...(modeScoped ? [] : ['circularSources'])
    ].includes(key)).forEach(key => faults.push({ at: [key], field: true }));
  }
  slots.forEach(slot => value(bindings[slot], [slot], slot === 'c_gb'));
  if (!modeScoped && Object.hasOwn(bindings, 'circularSources')) {
    const cards = bindings.circularSources;
    // The File cards of a Circular drawing (`_validate_circular_source_rows`):
    // a card that is no row of exactly its fields or repeats a uid is dropped;
    // a file that cannot be read leaves its card without it.
    if (!Array.isArray(cards)) return [{ at: null }];
    /** @type {Set<string>} */
    const uids = new Set();
    cards.forEach((/** @type {unknown} */ card, index) => {
      const uid = object(card) ? card.uid : null;
      if (!object(card) || Object.keys(card).length !== CIRCULAR_SOURCE_FIELDS.length
          || !CIRCULAR_SOURCE_FIELDS.every(key => Object.hasOwn(card, key))
          || typeof uid !== 'string' || !uid || uid !== uid.trim() || uids.has(uid)) {
        faults.push({ at: ['circularSources', index] });
        return;
      }
      uids.add(uid);
      CIRCULAR_SOURCE_FILE_SLOTS.forEach(slot => {
        if (card[slot] !== null && badLeaf(card[slot])) faults.push({ at: ['circularSources', index, slot] });
      });
    });
  }
  (Array.isArray(bindings.linearSeqs) ? bindings.linearSeqs : []).forEach((/** @type {WebFileBindingRow} */ sequence, index) => {
    ['gb', 'gff', 'fasta', 'depth', 'blast'].forEach(slot => value(sequence?.[slot], ['linearSeqs', index, slot]));
  });
  ['linearComparisons', 'linearCanonicalComparisons'].forEach(slot => {
    (Array.isArray(bindings[slot]) ? bindings[slot] : []).forEach((/** @type {WebFileBindingRow} */ row, index) => (
      value(row?.file, [slot, index, 'file'])));
  });
  return faults;
};

/**
 * The inventory with no Depth file that names one of `resourceIds`: each keeps
 * its place in its list as an empty position.
 * @param {WebFileBindings} bindings
 * @param {ReadonlySet<string>} resourceIds
 * @returns {WebFileBindings}
 */
export const webFileBindingsWithoutDepth = (bindings, resourceIds) => {
  /** @param {unknown} value @returns {unknown} */
  const kept = (value) => {
    if (Array.isArray(value)) return value.map(kept);
    return resourceIds.has(String(/** @type {{ resourceId?: unknown } | null} */ (value)?.resourceId)) ? null : value;
  };
  return {
    ...bindings,
    ...(Object.hasOwn(bindings, 'c_depth') ? { c_depth: kept(bindings.c_depth) } : {}),
    ...(Array.isArray(bindings.linearSeqs) ? {
      linearSeqs: bindings.linearSeqs.map((/** @type {WebFileBindingRow} */ row) => (
        row && Object.hasOwn(row, 'depth') ? { ...row, depth: kept(row.depth) } : row))
    } : {})
  };
};

/**
 * @param {{ bindings?: unknown } | null | undefined} webFiles
 * @param {{ readonly [resourceId: string]: RawSessionResourceDescriptor | null | undefined } | null} [resources]
 * @param {WebBindingInputSlot} [inputSlot] As `webFileBindingFaults`.
 * @returns {WebFileBindings | undefined} the `webFiles.bindings` object, checked
 */
export const validateWebFileBindings = (webFiles, resources = {}, inputSlot = 'c_gb') => {
  // R6: the typed Session-fields error; Load lists each fault's place before (`readSessionInputs`).
  if (webFileBindingFaults(webFiles, resources, inputSlot).length) throw diagnosticError('INPUT_INVALID', { field: 'schema', reason: 'FIELDS' });
  // No fault: the checks prove the schema and each file slot.
  return /** @type {WebFileBindings | undefined} */ (webFiles?.bindings);
};

/**
 * @param {SessionResourceTable | null | undefined} table
 * @param {unknown} resourceId
 * @returns {LeafBacking}
 */
const requireBacking = (table, resourceId) => {
  const backings = table && tableBackings.get(table);
  if (!backings) {
    throw new TypeError('An adopted current-session resource table is required.');
  }
  const normalizedId = String(resourceId || '').trim();
  const backing = backings.get(normalizedId);
  if (!backing) throw new Error(`Missing canonical resource: ${normalizedId}`);
  return backing;
};

/**
 * @param {LeafBacking} backing
 * @returns {Promise<Uint8Array>}
 */
const decodedResourceBytes = async (backing) => {
  const readFile = fileBackedReaderOf(backing.descriptor);
  if (readFile) return readFile();
  // `validateDescriptor` paired `encoding` with `data`: a depth table when the guard holds, base64 text otherwise.
  if (isEncodedDepthFileEntry(/** @type {DepthEntryProbe} */ (backing.descriptor))) {
    return textToBytes(decodeDepthText(/** @type {DepthFilePayload} */ (backing.descriptor.data)));
  }
  recordStructuralMetric('base64DecodeCount', 1, {
    resourceId: backing.resourceId,
    resourceName: backing.name
  });
  return base64ToBytesInTasks(/** @type {string} */ (backing.descriptor.data));
};

/**
 * @param {LeafBacking} backing
 * @returns {Promise<Uint8Array>}
 */
const materializeBytes = (backing) => {
  if (!backing.bytesPromise) {
    recordStructuralMetric('resourceByteReadCount', 1, {
      resourceId: backing.resourceId,
      resourceName: backing.name
    });
    backing.bytesPromise = Promise.resolve().then(async () => {
      /** @type {Uint8Array} */
      let bytes;
      try {
        bytes = await decodedResourceBytes(backing);
      } catch (error) {
        // A File's read keeps its own error (a refused read is INPUT_UNREADABLE, PD-OI-121).
        if (fileBackedReaderOf(backing.descriptor)) throw error;
        throw new Error(
          `Session resource ${backing.resourceId} (${backing.name}) contains invalid encoded data.`,
          { cause: error }
        );
      }
      recordStructuralMetric('decodedByteCount', bytes.byteLength, {
        resourceId: backing.resourceId,
        resourceName: backing.name
      });
      if (bytes.byteLength !== backing.size) {
        throw new Error(
          `Session resource ${backing.resourceId} (${backing.name}) byte size does not match its descriptor.`
        );
      }
      if (backing.checksum) {
        const actual = await sha256Hex(bytes);
        if (actual !== backing.checksum) {
          throw new Error(
            `Session resource ${backing.resourceId} (${backing.name}) checksum does not match.`
          );
        }
      }
      return bytes;
    });
    keepWhileRunning(backing, 'bytesPromise');
  }
  return backing.bytesPromise;
};

/**
 * What a Session file stores for a leaf resource (`storedResourceBytes`).
 * @param {LeafBacking} backing
 */
const storedOf = (backing) => storedResourceBytes(backing.descriptor);

/**
 * A read of a leaf whose stored bytes rest in a file (`StoredResourceBytes.blob`)
 * is shared while it runs and kept by nothing once it ends, as a File's is
 * (`holdReadBytes` in file-content-cache.js): the next read reads the file
 * again, so a History step that names a view of it keeps no bytes (OV-2754).
 * Any other leaf keeps its reads until its content is released.
 * @param {LeafBacking} backing
 * @param {'bytesPromise' | 'plainPromise' | 'textPromise'} field The read just started.
 */
const keepWhileRunning = (backing, field) => {
  const pending = backing[field];
  if (!pending || !storedOf(backing).blob) return;
  const forget = () => {
    if (backing[field] === pending) backing[field] = null;
  };
  pending.then(forget, forget);
};

const MIB = 1024 * 1024;
/**
 * The plain bytes of an input's bytes, the one expansion of an input in the
 * page (a File's, `readFileBytes` in file-content-cache.js, and a Session
 * resource's), which only the page's readers of plain bytes or text call and
 * whose result nothing keeps: a gzip input streams through the browser's
 * gzip decoder to at most Python's `MAX_EXPANDED_INPUT_BYTES`, and past it
 * rejects with `INPUT_UNREADABLE` (`GZIP_EXPANDED_SIZE`), naming the file and
 * the cap in MiB (PD-OI-116), as Python refuses it at Generate; any other
 * bytes are given as they are.
 * @param {string} name The input's name.
 * @param {Uint8Array} bytes
 * @returns {Promise<Uint8Array>}
 */
export const expandInputBytes = (name, bytes) => {
  if (!isGzipBytes(bytes)) return Promise.resolve(bytes);
  recordStructuralMetric('inputExpansionCount', 1, { resourceName: name });
  const { maxBytes, overCap } = inputExpansionCap(name);
  return expandGzipBytes(/** @type {Uint8Array<ArrayBuffer>} */ (bytes), maxBytes, overCap);
};

/**
 * The cap of a gzip input's expanded size and its refusal, which names the
 * file and the cap in MiB (PD-OI-116); the one owner of both for the page.
 * @param {string} name The input's name.
 * @returns {{ maxBytes: number, overCap: () => Error }}
 */
const inputExpansionCap = (name) => {
  const maxBytes = testByteCap('maxExpandedInputBytes') ?? MAX_EXPANDED_INPUT_BYTES;
  return {
    maxBytes,
    overCap: () => diagnosticError('INPUT_UNREADABLE',
      { reason: 'GZIP_EXPANDED_SIZE', inputName: name, limit: Math.floor(maxBytes / MIB) }, { stage: 'read' })
  };
};

/**
 * Checks a chosen input against the cap of `expandInputBytes`, keeping nothing
 * expanded: a gzip input whose expanded size passes it rejects with the same
 * refusal, so a chosen file is refused right after it is chosen, whatever its
 * kind (PD-OI-116); any other bytes pass.
 * @param {string} name The input's name.
 * @param {Uint8Array} bytes Its bytes as chosen.
 * @returns {Promise<void>}
 */
export const checkInputExpansion = (name, bytes) => {
  if (!isGzipBytes(bytes)) return Promise.resolve();
  recordStructuralMetric('inputExpansionCheckCount', 1, { resourceName: name });
  const { maxBytes, overCap } = inputExpansionCap(name);
  return checkGzipExpandedSize(/** @type {Uint8Array<ArrayBuffer>} */ (bytes), maxBytes, overCap);
};

/**
 * Whether an error is the cap refusal of `expandInputBytes`: the answer Python gives the same input too, so a
 * browser fast path gives it instead of declining to the Worker (R4).
 * @param {unknown} error
 * @returns {boolean}
 */
export const isInputExpansionRefusal = (error) => {
  const { code, context } = /** @type {{ code?: unknown, context?: { reason?: unknown } }} */ (error ?? {});
  return code === 'INPUT_UNREADABLE' && context?.reason === 'GZIP_EXPANDED_SIZE';
};

/**
 * The plain bytes of a leaf resource: its stored bytes, a gzip input's
 * expanded (`expandInputBytes`). Those of a plain resource are kept until its
 * content is released, and the reads that run meanwhile share them; a gzip
 * input's expanded bytes are not kept once its read ends, so its next read
 * expands again.
 * @param {LeafBacking} backing
 * @returns {Promise<Uint8Array>}
 */
const materializePlainBytes = (backing) => {
  if (backing.plainPromise) return backing.plainPromise;
  const plain = materializeBytes(backing).then((bytes) => {
    if (!isGzipBytes(bytes)) return bytes;
    if (backing.plainPromise === plain) backing.plainPromise = null;
    return expandInputBytes(backing.name, bytes);
  });
  backing.plainPromise = plain;
  keepWhileRunning(backing, 'plainPromise');
  return plain;
};

/**
 * The text of a leaf resource's plain bytes (`materializePlainBytes`), kept as
 * they are: a gzip input's text is not.
 * @param {LeafBacking} backing
 * @returns {Promise<string>}
 */
const materializeText = (backing) => {
  if (backing.textPromise) return backing.textPromise;
  recordStructuralMetric('resourceTextReadCount', 1, {
    resourceId: backing.resourceId,
    resourceName: backing.name
  });
  // Started with the bytes' read, so the two share it; a failed read rejects the text.
  const plain = materializePlainBytes(backing);
  plain.catch(() => {});
  const text = materializeBytes(backing).then(async (bytes) => {
    if (isGzipBytes(bytes) && backing.textPromise === text) backing.textPromise = null;
    return bytesToText(await plain);
  });
  backing.textPromise = text;
  keepWhileRunning(backing, 'textPromise');
  return text;
};

/**
 * @param {Record<string, unknown> | null | undefined} resources Unvalidated; each descriptor is checked here.
 * @returns {SessionResourceTable}
 */
export const adoptCurrentSessionResources = (resources) => {
  if (!resources || typeof resources !== 'object' || Array.isArray(resources)) {
    throw new Error('Current session resources must be an object.');
  }
  /** @type {Map<string, LeafBacking>} */
  const backings = new Map();
  /** @type {Set<string>} */
  const seen = new Set();
  Object.entries(resources).forEach(([rawResourceId, descriptor]) => {
    const resourceId = tableResourceId(rawResourceId, seen);
    if (resourceId === null) {
      throw new Error('Current session resource IDs must be unique non-empty strings.');
    }
    backings.set(resourceId, validateDescriptor(resourceId, /** @type {RawSessionResourceDescriptor | null | undefined} */ (descriptor)));
  });
  /** @type {SessionResourceTable} */
  const table = Object.freeze({});
  tableBackings.set(table, backings);
  return table;
};

/**
 * @param {SessionResourceTable | null | undefined} table
 * @param {unknown} resourceId
 */
export const adoptedSessionResourceDescriptor = (table, resourceId) => (
  requireBacking(table, resourceId).descriptor
);

/**
 * @typedef {Readonly<Record<string, never>>} SessionResourceTable Opaque handle from adoptCurrentSessionResources
 * @typedef {{ name?: unknown, type?: unknown, lastModified?: unknown }} SessionResourceFileMetadata
 * @typedef {{ encoding?: unknown, data?: { schema?: unknown } | null }} DepthEntryProbe The part of a descriptor `isEncodedDepthFileEntry` reads; a string `data` has no `schema`.
 */

/**
 * @param {SessionResourceTable} table
 * @param {unknown} resourceId Checked against the table's own ids.
 * @param {SessionResourceFileMetadata} [metadata]
 * @returns {SessionResourceFileView}
 */
export const createSessionResourceFileView = (table, resourceId, metadata = {}) => {
  const backing = requireBacking(table, resourceId);
  const view = Object.freeze({
    name: metadata.name === undefined ? backing.name : String(metadata.name),
    type: metadata.type === undefined ? backing.type : String(metadata.type || ''),
    size: backing.size,
    lastModified: metadata.lastModified === undefined
      ? backing.lastModified
      : Number(metadata.lastModified) || 0
  });
  fileViewBackings.set(view, backing);
  return view;
};

/**
 * @param {unknown} file
 * @returns {file is SessionResourceFileView}
 */
export const isSessionResourceFileView = (file) => Boolean(
  file && (typeof file === 'object' || typeof file === 'function')
  && fileViewBackings.has(file)
);

/**
 * @param {FileLike | null | undefined} file
 * @param {SessionResourceDescriptor | null | undefined} descriptor the saved descriptor
 * @returns {boolean}
 */
export const matchesSessionResourceDescriptor = (file, descriptor) => {
  // A missing file is no key: `get` returns undefined for it.
  const backing = fileViewBackings.get(/** @type {FileLike} */ (file));
  // Native Files remain bound by the successful run's object identity.
  if (!backing || backing.dirty) return true;
  if (!descriptor) return false;
  return backing.descriptor.encoding === descriptor.encoding
    && resourcePayloadKey(backing.descriptor) === resourcePayloadKey(descriptor);
};

/** @param {unknown} file A file slot's value; only a clean Session resource view has a source. */
export const sessionResourceSource = (file) => {
  // The lookup is by identity. A value that is no key (a missing file, a
  // descriptor, a `File`) has no backing: `get` returns undefined for it.
  const backing = fileViewBackings.get(/** @type {FileLike} */ (file));
  if (!backing || backing.dirty) return null;
  return Object.freeze({
    resourceId: backing.resourceId,
    descriptor: backing.descriptor,
    readBytes: () => materializeBytes(backing)
  });
};

/**
 * What a read of `file` is kept under (`readOncePerFile`): a Session resource view by its resource, so every view of
 * one resource shares one read of it; any other file by itself. A resource's bytes never change once its table is
 * adopted (its descriptor is fixed then), so a read of one view is a read of every view; a new Load adopts a new
 * table, and a replaced file is a new File.
 * @param {FileLike} file
 * @returns {object}
 */
export const sessionResourceReadOwner = (file) => fileViewBackings.get(file) ?? file;

/**
 * The plain bytes of a Session resource view, which the page's readers read:
 * its stored bytes, a gzip input's expanded at every read
 * (`expandInputBytes`).
 * @param {FileLike | null | undefined} file
 * @returns {Promise<Uint8Array> | null}
 */
export const readSessionResourceBytes = (file) => {
  // A missing file is no key: `get` returns undefined for it.
  const backing = fileViewBackings.get(/** @type {FileLike} */ (file));
  if (!backing) return null;
  return materializePlainBytes(backing);
};

/**
 * The stored bytes of a Session resource view, the size its descriptor
 * counts (a gzip input as stored): what Generate transfers for a view and
 * what its content token hashes.
 * @param {FileLike | null | undefined} file
 * @returns {Promise<Uint8Array> | null}
 */
export const readSessionResourceStoredBytes = (file) => {
  // A missing file is no key: `get` returns undefined for it.
  const backing = fileViewBackings.get(/** @type {FileLike} */ (file));
  if (!backing) return null;
  return materializeBytes(backing);
};

/**
 * What a Session file stores for the resource of a Session resource view
 * (`storedResourceBytes`), or null for any other value.
 * @param {FileLike | null | undefined} file
 * @returns {StoredResourceBytes | null}
 */
export const sessionResourceStoredBytes = (file) => {
  // A missing file is no key: `get` returns undefined for it.
  const backing = fileViewBackings.get(/** @type {FileLike} */ (file));
  return backing ? storedOf(backing) : null;
};

/** @type {WeakMap<FileLike, FileBackedReader>} */
const viewReaders = new WeakMap();
/**
 * The one reader of a clean Session resource view's bytes, or null for any
 * other value. A file-backed descriptor Generate made of the view keeps it
 * (`fileBackedDescriptor`), so the view knows that descriptor by identity
 * (`matchesSessionResourceDescriptor`), with no read.
 * @param {unknown} file
 * @returns {FileBackedReader | null}
 */
export const sessionResourceViewReader = (file) => {
  const view = /** @type {FileLike} */ (file);
  const backing = fileViewBackings.get(view);
  if (!backing || backing.dirty) return null;
  let reader = viewReaders.get(view);
  if (!reader) {
    viewReaders.set(view, reader = () => /** @type {Promise<Uint8Array>} */ (readSessionResourceStoredBytes(view)));
    // A view stores what its resource stores (an entry a Load copied streams), so a Save after a Generate does too.
    if (fileBackedReaderOf(backing.descriptor)) fileBackedStored.set(reader, storedResourceBytes(backing.descriptor));
  }
  return reader;
};

/**
 * The text of a Session resource view: its plain bytes (a gzip input
 * expanded, `materializeText`).
 * @param {FileLike | null | undefined} file
 * @returns {Promise<string> | null}
 */
export const readSessionResourceText = (file) => {
  // A missing file is no key: `get` returns undefined for it.
  const backing = fileViewBackings.get(/** @type {FileLike} */ (file));
  if (!backing) return null;
  return materializeText(backing);
};

/** @param {FileLike | null | undefined} file */
export const releaseSessionResourceContent = (file) => {
  // A missing file is no key: `get` returns undefined for it.
  const backing = fileViewBackings.get(/** @type {FileLike} */ (file));
  if (!backing) return false;
  backing.bytesPromise = null;
  backing.plainPromise = null;
  backing.textPromise = null;
  backing.filePromise = null;
  return true;
};
