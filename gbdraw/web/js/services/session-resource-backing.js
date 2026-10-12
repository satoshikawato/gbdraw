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

/** @param {unknown} value */
const normalizedChecksum = (value) => {
  if (value === null || value === undefined || value === '') return null;
  if (typeof value !== 'string') {
    throw new Error('Session resource checksum must be a string.');
  }
  const normalized = value.trim().toLowerCase().replace(/^sha256:/, '');
  if (!/^[0-9a-f]{64}$/.test(normalized)) {
    throw new Error('Session resource checksum must be a SHA-256 digest.');
  }
  return normalized;
};

/**
 * @param {string} resourceId
 * @param {RawSessionResourceDescriptor} descriptor
 * @returns {LeafBacking}
 */
const validateDescriptor = (resourceId, descriptor) => {
  if (!descriptor || typeof descriptor !== 'object' || Array.isArray(descriptor)) {
    throw new Error(`Session resource ${resourceId} must be an object.`);
  }
  const supportedPayload = descriptor.encoding === 'base64'
    ? typeof descriptor.data === 'string'
    // A `data` that is not an object (a string, or nothing) has no `schema`, so the guard reads false.
    : isEncodedDepthFileEntry(/** @type {{ encoding?: unknown, data?: { schema?: unknown } | null }} */ (descriptor));
  if (!supportedPayload) {
    throw new Error(`Session resource ${resourceId} has an unsupported encoded payload.`);
  }
  if (typeof descriptor.size !== 'number' || !Number.isSafeInteger(descriptor.size) || descriptor.size < 0) {
    throw new Error(`Session resource ${resourceId} has an invalid declared byte size.`);
  }
  if (
    descriptor.lastModified !== undefined
    && (typeof descriptor.lastModified !== 'number' || !Number.isFinite(descriptor.lastModified) || descriptor.lastModified < 0)
  ) {
    throw new Error(`Session resource ${resourceId} has an invalid lastModified value.`);
  }
  return {
    // The checks above establish `encoding`, `data` and `size`.
    descriptor: /** @type {SessionResourceDescriptor} */ (descriptor),
    resourceId,
    kind: String(descriptor.kind || 'web-file'),
    name: String(descriptor.name || 'file'),
    type: String(descriptor.type || 'application/octet-stream'),
    encoding: String(descriptor.encoding || ''),
    size: descriptor.size,
    lastModified: Number(descriptor.lastModified) || 0,
    checksum: normalizedChecksum(descriptor.checksum),
    bytesPromise: null,
    textPromise: null,
    filePromise: null,
    dirty: false
  };
};

// Binding grammar belongs beside the File backing; envelope/version and
// comparison-plan authority remain with their existing owners.
/**
 * @param {{ bindings?: unknown } | null | undefined} webFiles
 * @param {Record<string, RawSessionResourceDescriptor>} [resources]
 * @returns {Record<string, unknown> | undefined} the `webFiles.bindings` object, checked
 */
export const validateWebFileBindings = (webFiles, resources = {}) => {
  if (!webFiles || !Object.hasOwn(webFiles, 'bindings')) return undefined;
  const bindings = webFiles.bindings;
  /**
   * @param {unknown} value
   * @returns {value is Record<string, unknown>}
   */
  const object = value => value !== null && typeof value === 'object' && !Array.isArray(value);
  if (!object(bindings)) throw new Error('Session webFiles.bindings must be an object.');
  if (!(/** @type {unknown[]} */ ([1, 2])).includes(bindings.schema)) throw new Error('Unsupported Web file binding schema.');
  const current = bindings.schema === 2;
  if (current && !Object.hasOwn(bindings, 'c_gb')) {
    throw new Error('Web binding schema 2 requires c_gb.');
  }
  /** @param {Record<string, unknown>} value */
  const metadata = value => {
    if (typeof value.name !== 'string' || typeof value.type !== 'string'
      || typeof value.lastModified !== 'number' || !Number.isFinite(value.lastModified) || value.lastModified < 0) {
      throw new Error('Invalid Web file binding metadata.');
    }
  };
  /** @param {unknown} value */
  const leaf = value => {
    if (!object(value)) throw new Error('Web file bindings must be objects or null.');
    if (Object.hasOwn(value, 'kind') || Object.hasOwn(value, 'components')) {
      throw new Error('A composite component must be an ordinary Web file binding.');
    }
    const id = current ? value.resourceId : String(value.resourceId || '').trim();
    if (typeof id !== 'string' || !id || id !== id.trim()) {
      throw new Error('A Web file binding requires a canonical resourceId.');
    }
    if (!Object.hasOwn(resources, id)) throw new Error(`Missing canonical resource: ${id}`);
    if (current) {
      if (Object.keys(value).some(key => !['resourceId', 'name', 'type', 'lastModified'].includes(key))) {
        throw new Error('Unknown Web file binding field.');
      }
      metadata(value);
      validateDescriptor(id, resources[id]);
    }
  };
  /**
   * @param {unknown} binding
   * @param {boolean} [compositeAllowed]
   * @returns {void}
   */
  const value = (binding, compositeAllowed = false) => {
    if (binding === null || binding === undefined) return;
    if (Array.isArray(binding)) {
      binding.forEach(item => value(item));
    } else if (object(binding) && Object.hasOwn(binding, 'kind')) {
      if (!current || !compositeAllowed || binding.kind !== 'composite') {
        throw new Error('Unsupported Web composite file binding.');
      }
      if (Object.keys(binding).some(key => !['kind', 'components', 'name', 'type', 'lastModified'].includes(key))) {
        throw new Error('Unknown or mixed composite binding fields.');
      }
      metadata(binding);
      if (!Array.isArray(binding.components) || binding.components.length < 2) {
        throw new Error('A composite binding requires at least two components.');
      }
      // `leaf` checks each component's `resourceId`, so it is a string here.
      binding.components.forEach((/** @type {{ resourceId: string }} */ component) => {
        leaf(component);
        if (resources[component.resourceId].encoding !== 'base64') {
          throw new Error('Composite components require base64 resources.');
        }
      });
    } else {
      leaf(binding);
    }
  };
  const slots = [
    'c_gb', 'c_gff', 'c_fasta', 'c_depth', 'c_conservation_blasts',
    'c_conservation_fastas', 'c_conservation_sequence_sources', 'd_color',
    't_color', 'blacklist', 'whitelist', 'qualifier_priority'
  ];
  if (current && Object.keys(bindings).some(key => ![
    'schema', ...slots, 'c_conservation_blasts_source', 'linearSeqs', 'linearComparisons'
  ].includes(key))) throw new Error('Unknown Web binding inventory field.');
  slots.forEach(slot => value(bindings[slot], slot === 'c_gb'));
  (Array.isArray(bindings.linearSeqs) ? bindings.linearSeqs : []).forEach((/** @type {Record<string, unknown> | null | undefined} */ sequence) => {
    ['gb', 'gff', 'fasta', 'depth', 'blast'].forEach(slot => value(sequence?.[slot]));
  });
  ['linearComparisons', 'linearCanonicalComparisons'].forEach(slot => {
    (Array.isArray(bindings[slot]) ? bindings[slot] : []).forEach(row => value(row?.file));
  });
  return bindings;
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
 * @param {Record<string, RawSessionResourceDescriptor> | null | undefined} resources
 * @returns {SessionResourceTable}
 */
export const adoptCurrentSessionResources = (resources) => {
  if (!resources || typeof resources !== 'object' || Array.isArray(resources)) {
    throw new Error('Current session resources must be an object.');
  }
  /** @type {Map<string, LeafBacking>} */
  const backings = new Map();
  Object.entries(resources).forEach(([rawResourceId, descriptor]) => {
    const resourceId = String(rawResourceId || '').trim();
    if (!resourceId || backings.has(resourceId)) {
      throw new Error('Current session resource IDs must be unique non-empty strings.');
    }
    backings.set(resourceId, validateDescriptor(resourceId, descriptor));
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
 * @param {string} resourceId
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
  /** @param {SessionResourceDescriptor | null} entry */
  const matches = (entry) => entry?.encoding === descriptor.encoding
    && entry?.data === descriptor.data;
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
