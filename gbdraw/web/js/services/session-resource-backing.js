// @ts-check
import {
  base64DecodedLastByte,
  base64ToBytesInTasks,
  bytesToText,
  sha256Hex,
  textToBytes
} from './byte-utils.js';
import {
  decodeDepthText,
  isEncodedDepthFileEntry
} from './depth-file-codec.js';
import { recordStructuralMetric } from './runtime-test-hooks.js';
import { diagnosticError } from '../utils/error-normalization.js';

/** @import { DepthFilePayload } from './depth-file-codec.js' */

/**
 * A Session resource as the Session file stores it, before any check: every
 * field is unvalidated.
 * @typedef {{
 *   kind?: unknown, name?: unknown, type?: unknown, encoding?: unknown, data?: unknown,
 *   size?: unknown, lastModified?: unknown, checksum?: unknown
 * }} RawSessionResourceDescriptor
 */
/**
 * A Session resource after `validateDescriptor`: a base64 text, or an encoded
 * depth table (`encoding` names which `data` holds).
 * @typedef {{
 *   kind?: unknown, name?: unknown, type?: unknown, encoding: string, data: string | DepthFilePayload,
 *   size: number, lastModified?: number, checksum?: unknown
 * }} SessionResourceDescriptor
 */
/**
 * The stand-in for a `File` that a Session resource backs.
 * @typedef {Readonly<{ name: string, type: string, size: number, lastModified: number }>} SessionResourceFileView
 */
/**
 * A `File` or a Session resource view, the two things a file slot holds.
 * @typedef {{
 *   name?: string, type?: string, size?: number, lastModified?: number,
 *   arrayBuffer?: () => Promise<ArrayBuffer>, text?: () => Promise<string>
 * }} FileLike
 */
/**
 * @typedef {object} BackingFields
 * @property {string} resourceId
 * @property {string} name
 * @property {string} type
 * @property {number} size
 * @property {number} lastModified
 * @property {Promise<Uint8Array> | null} bytesPromise
 * @property {Promise<string> | null} textPromise
 * @property {null} filePromise
 * @property {boolean} dirty
 */
/**
 * The backing of one resource: its descriptor and the content it decoded.
 * @typedef {BackingFields & {
 *   descriptor: SessionResourceDescriptor, kind: string, encoding: string, checksum: string | null,
 *   descriptors?: undefined, sourceBackings?: undefined, readBytes?: undefined, readText?: undefined
 * }} LeafBacking
 */
/**
 * One component of a combined view: the part's descriptor and a reader of its bytes.
 * @typedef {Readonly<{
 *   resourceId: string, descriptor: SessionResourceDescriptor, name: string, type: string,
 *   lastModified: number, readBytes: () => Promise<Uint8Array>
 * }>} CompositePart
 */
/**
 * The backing of a combined view of several base64 resources. Its `descriptor`
 * stays null until `matchesSessionResourceDescriptor` verifies a serialized form.
 * @typedef {BackingFields & {
 *   checksum: null, descriptor: SessionResourceDescriptor | null, descriptors: ReadonlyArray<CompositePart>,
 *   sourceBackings: LeafBacking[], readBytes?: () => Promise<Uint8Array>, readText?: () => Promise<string>
 * }} CompositeBacking
 */

/** @type {WeakMap<SessionResourceTable, Map<string, LeafBacking>>} */
const tableBackings = new WeakMap();
/** @type {WeakMap<WeakKey, LeafBacking | CompositeBacking>} */
const fileViewBackings = new WeakMap();
const LF_BYTE = 0x0A;

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
  const supportedPayload = descriptor.encoding === 'base64'
    ? typeof descriptor.data === 'string'
    // A `data` that is not an object (a string, or nothing) has no `schema`, so the guard reads false.
    : isEncodedDepthFileEntry(/** @type {{ encoding?: unknown, data?: { schema?: unknown } | null }} */ (descriptor));
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
 * The schema 2 `c_gb` binding that reads two or more base64 resources as one file.
 * @typedef {Readonly<{
 *   kind: 'composite', name: string, type: string, lastModified: number,
 *   components: Array<{ resourceId: string, name: string, type: string, lastModified: number }>
 * }>} WebFileCompositeBinding
 */
/**
 * A row of `linearSeqs` or of a comparison list: only its file slots are checked.
 * @typedef {Readonly<Record<string, unknown>> | null | undefined} WebFileBindingRow
 */
/**
 * `webFiles.bindings` after `validateWebFileBindings`. Each file slot holds null,
 * a leaf, or a list of them (a JSON list has no `kind`); only `c_gb` may hold a composite.
 * @typedef {Readonly<{
 *   schema: 1 | 2,
 *   c_gb?: WebFileCompositeBinding | WebFileLeafBinding | (ReadonlyArray<unknown> & { readonly kind?: undefined }) | null,
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
 * The one check of a Web input inventory: Load drops what it names
 * (`readWebFileBindings` in session-authority.js), Save and the runtime refuse
 * it (`validateWebFileBindings`).
 * @param {{ bindings?: unknown } | null | undefined} webFiles
 * @param {{ readonly [resourceId: string]: RawSessionResourceDescriptor | null | undefined } | null} [resources]
 *   The Session's `resources`, unvalidated: read by the id of each binding, and each entry read is checked.
 * @returns {WebFileBindingFault[]}
 */
export const webFileBindingFaults = (webFiles, resources = {}) => {
  if (!webFiles || !Object.hasOwn(webFiles, 'bindings')) return [];
  const bindings = webFiles.bindings;
  /**
   * @param {unknown} value
   * @returns {value is Record<string, unknown>}
   */
  const object = value => value !== null && typeof value === 'object' && !Array.isArray(value);
  // The inventory: an object of a schema a writer wrote; schema 2 holds `c_gb`.
  if (!object(bindings) || !(/** @type {unknown[]} */ ([1, 2])).includes(bindings.schema)) return [{ at: null }];
  const current = bindings.schema === 2;
  if (current && !Object.hasOwn(bindings, 'c_gb')) return [{ at: null }];
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
    // A composite is schema 2 `c_gb` only, with only its fields and two or more base64 components.
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
  const slots = [
    'c_gb', 'c_gff', 'c_fasta', 'c_depth', 'c_conservation_blasts',
    'c_conservation_fastas', 'c_conservation_sequence_sources', 'd_color',
    't_color', 'blacklist', 'whitelist', 'qualifier_priority'
  ];
  if (current) {
    Object.keys(bindings).filter(key => ![
      'schema', ...slots, 'c_conservation_blasts_source', 'linearSeqs', 'linearComparisons'
    ].includes(key)).forEach(key => faults.push({ at: [key], field: true }));
  }
  slots.forEach(slot => value(bindings[slot], [slot], slot === 'c_gb'));
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
 * @param {{ bindings?: unknown } | null | undefined} webFiles
 * @param {{ readonly [resourceId: string]: RawSessionResourceDescriptor | null | undefined } | null} [resources]
 * @returns {WebFileBindings | undefined} the `webFiles.bindings` object, checked
 */
export const validateWebFileBindings = (webFiles, resources = {}) => {
  // R6: the typed Session-fields error; Load lists each fault's place before (`readSessionInputs`).
  if (webFileBindingFaults(webFiles, resources).length) throw diagnosticError('INPUT_INVALID', { field: 'schema', reason: 'FIELDS' });
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
  }
  return backing.bytesPromise;
};

/**
 * @param {LeafBacking} backing
 * @returns {Promise<string>}
 */
const materializeText = (backing) => {
  if (!backing.textPromise) {
    recordStructuralMetric('resourceTextReadCount', 1, {
      resourceId: backing.resourceId,
      resourceName: backing.name
    });
    backing.textPromise = materializeBytes(backing).then((bytes) => bytesToText(bytes));
  }
  return backing.textPromise;
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
 * @param {SessionResourceTable} table
 * @param {Array<SessionResourceFileMetadata & { resourceId: string }>} components
 * @param {SessionResourceFileMetadata} [metadata]
 * @returns {SessionResourceFileView | null}
 */
export const createCombinedSessionResourceFileView = (
  table,
  components,
  metadata = {}
) => {
  const views = (Array.isArray(components) ? components : []).map(component => (
    createSessionResourceFileView(table, component.resourceId, component)
  ));
  // Each view above was created by `createSessionResourceFileView`, so it maps to a leaf backing.
  const backings = /** @type {LeafBacking[]} */ (views.map(view => fileViewBackings.get(view)));
  if (backings.length === 0) return null;
  if (backings.length === 1) {
    return createSessionResourceFileView(table, backings[0].resourceId, { ...components[0], ...metadata });
  }
  const appendedLineFeeds = backings.reduce((count, backing) => (
    count + Number(
      backing.size > 0
      // A combined view holds base64 resources only (`validateWebFileBindings`).
      && base64DecodedLastByte(/** @type {string} */ (backing.descriptor.data), backing.size) !== LF_BYTE
    )
  ), 0);
  /** @type {CompositeBacking} */
  const composite = {
    resourceId: backings.map(({ resourceId }) => resourceId).join('+'),
    name: metadata.name === undefined ? 'canonical-circular-records.gb' : String(metadata.name),
    type: metadata.type === undefined ? 'text/plain' : String(metadata.type),
    size: backings.reduce((total, backing) => total + backing.size, 0)
      + appendedLineFeeds,
    lastModified: metadata.lastModified === undefined
      ? Math.max(0, ...views.map(view => view.lastModified))
      : Number(metadata.lastModified) || 0,
    checksum: null,
    descriptor: null,
    descriptors: Object.freeze(backings.map((backing, index) => Object.freeze({
      resourceId: backing.resourceId,
      descriptor: backing.descriptor,
      name: views[index].name,
      type: views[index].type,
      lastModified: views[index].lastModified,
      readBytes: () => materializeBytes(backing)
    }))),
    sourceBackings: backings,
    bytesPromise: null,
    textPromise: null,
    filePromise: null,
    dirty: false
  };
  const readParts = () => Promise.all(backings.map((backing) => materializeBytes(backing)))
    .then((parts) => {
      const lengths = parts.map((part) => (
        part.byteLength > 0 && part[part.byteLength - 1] !== LF_BYTE
          ? part.byteLength + 1
          : part.byteLength
      ));
      const combined = new Uint8Array(lengths.reduce((total, length) => total + length, 0));
      let offset = 0;
      parts.forEach((part, index) => {
        combined.set(part, offset);
        offset += part.byteLength;
        if (lengths[index] > part.byteLength) combined[offset++] = LF_BYTE;
      });
      return combined;
    });
  composite.readBytes = () => {
    if (!composite.bytesPromise) composite.bytesPromise = readParts();
    return composite.bytesPromise;
  };
  composite.readText = () => {
    if (!composite.textPromise) {
      recordStructuralMetric('resourceTextReadCount', 1, {
        resourceId: composite.resourceId,
        resourceName: composite.name
      });
      // `readBytes` was assigned just above.
      composite.textPromise = /** @type {() => Promise<Uint8Array>} */ (composite.readBytes)().then((bytes) => bytesToText(bytes));
    }
    return composite.textPromise;
  };
  const view = Object.freeze({
    name: composite.name,
    type: composite.type,
    size: composite.size,
    lastModified: composite.lastModified
  });
  fileViewBackings.set(view, composite);
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
  // A combined view has no descriptor of its own until one is verified below: it equals none (OV-574).
  /** @param {SessionResourceDescriptor | null} entry */
  const matches = (entry) => entry !== null && entry.encoding === descriptor.encoding
    && entry.data === descriptor.data;
  if (matches(backing.descriptor)) return true;
  if (!backing.sourceBackings) return false;
  if (backing.sourceBackings.some((part) => matches(part.descriptor))) return true;
  if (descriptor.encoding !== 'base64' || descriptor.size !== backing.size) return false;

  /** @param {{ resourceId: string, name: string, descriptor: SessionResourceDescriptor }} entry */
  const decode = (entry) => {
    recordStructuralMetric('base64DecodeCount', 1, {
      resourceId: entry.resourceId, resourceName: entry.name
    });
    // The parts of a combined view are base64 text (`validateWebFileBindings`).
    const binary = atob(/** @type {string} */ (entry.descriptor.data));
    recordStructuralMetric('decodedByteCount', binary.length, {
      resourceId: entry.resourceId, resourceName: entry.name
    });
    return binary;
  };
  const combined = decode({ ...backing, descriptor });
  if (combined.length !== descriptor.size) return false;
  let offset = 0;
  for (const part of backing.sourceBackings) {
    const binary = decode(part);
    if (binary.length !== part.size || !combined.startsWith(binary, offset)) return false;
    offset += binary.length;
    if (binary.length && binary.charCodeAt(binary.length - 1) !== LF_BYTE) {
      if (combined.charCodeAt(offset++) !== LF_BYTE) return false;
    }
  }
  if (offset !== combined.length) return false;
  // The immutable composite's existing descriptor slot can now reference its
  // verified serialized representation, without retaining decoded content.
  backing.descriptor = descriptor;
  return true;
};

/** @param {unknown} file A file slot's value; only a clean Session resource view has a source. */
export const sessionResourceSource = (file) => {
  // The lookup is by identity. A value that is no key (a missing file, a
  // descriptor, a `File`) has no backing: `get` returns undefined for it.
  const backing = fileViewBackings.get(/** @type {FileLike} */ (file));
  if (!backing || backing.dirty) return null;
  if (Array.isArray(backing.descriptors)) {
    return Object.freeze({ descriptors: backing.descriptors });
  }
  return Object.freeze({
    resourceId: backing.resourceId,
    descriptor: backing.descriptor,
    // Only a combined view has `descriptors`, so this is a leaf backing.
    readBytes: () => materializeBytes(/** @type {LeafBacking} */ (backing))
  });
};

/**
 * @param {FileLike | null | undefined} file
 * @returns {Promise<Uint8Array> | null}
 */
export const readSessionResourceBytes = (file) => {
  // A missing file is no key: `get` returns undefined for it.
  const backing = fileViewBackings.get(/** @type {FileLike} */ (file));
  if (!backing) return null;
  if (typeof backing.readBytes === 'function') return backing.readBytes();
  // Every combined view defines `readBytes`, so this is a leaf backing.
  return materializeBytes(/** @type {LeafBacking} */ (backing));
};

/**
 * @param {FileLike | null | undefined} file
 * @returns {Promise<string> | null}
 */
export const readSessionResourceText = (file) => {
  // A missing file is no key: `get` returns undefined for it.
  const backing = fileViewBackings.get(/** @type {FileLike} */ (file));
  if (!backing) return null;
  if (typeof backing.readText === 'function') return backing.readText();
  // Every combined view defines `readText`, so this is a leaf backing.
  return materializeText(/** @type {LeafBacking} */ (backing));
};

/** @param {FileLike | null | undefined} file */
export const releaseSessionResourceContent = (file) => {
  // A missing file is no key: `get` returns undefined for it.
  const backing = fileViewBackings.get(/** @type {FileLike} */ (file));
  if (!backing) return false;
  /**
   * @param {LeafBacking | CompositeBacking} candidate
   * @returns {void}
   */
  const releaseBacking = (candidate) => {
    candidate.bytesPromise = null;
    candidate.textPromise = null;
    candidate.filePromise = null;
    candidate.sourceBackings?.forEach(releaseBacking);
  };
  releaseBacking(backing);
  return true;
};
