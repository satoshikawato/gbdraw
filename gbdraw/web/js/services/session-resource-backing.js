import {
  base64DecodedLastByte,
  base64ToBytes,
  bytesToText,
  sha256Hex,
  textToBytes
} from './byte-utils.js';
import {
  decodeDepthText,
  isEncodedDepthFileEntry
} from './depth-file-codec.js';
import { recordStructuralMetric } from './runtime-test-hooks.js';

const tableBackings = new WeakMap();
const fileViewBackings = new WeakMap();
const LF_BYTE = 0x0A;

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

const validateDescriptor = (resourceId, descriptor) => {
  if (!descriptor || typeof descriptor !== 'object' || Array.isArray(descriptor)) {
    throw new Error(`Session resource ${resourceId} must be an object.`);
  }
  const supportedPayload = descriptor.encoding === 'base64'
    ? typeof descriptor.data === 'string'
    : isEncodedDepthFileEntry(descriptor);
  if (!supportedPayload) {
    throw new Error(`Session resource ${resourceId} has an unsupported encoded payload.`);
  }
  if (!Number.isSafeInteger(descriptor.size) || descriptor.size < 0) {
    throw new Error(`Session resource ${resourceId} has an invalid declared byte size.`);
  }
  if (
    descriptor.lastModified !== undefined
    && (!Number.isFinite(descriptor.lastModified) || descriptor.lastModified < 0)
  ) {
    throw new Error(`Session resource ${resourceId} has an invalid lastModified value.`);
  }
  return {
    descriptor,
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
export const validateWebFileBindings = (webFiles, resources = {}) => {
  if (!webFiles || !Object.hasOwn(webFiles, 'bindings')) return undefined;
  const bindings = webFiles.bindings;
  const object = value => value !== null && typeof value === 'object' && !Array.isArray(value);
  if (!object(bindings)) throw new Error('Session webFiles.bindings must be an object.');
  if (![1, 2].includes(bindings.schema)) throw new Error('Unsupported Web file binding schema.');
  const current = bindings.schema === 2;
  if (current && !Object.hasOwn(bindings, 'c_gb')) {
    throw new Error('Web binding schema 2 requires c_gb.');
  }
  const metadata = value => {
    if (typeof value.name !== 'string' || typeof value.type !== 'string'
      || !Number.isFinite(value.lastModified) || value.lastModified < 0) {
      throw new Error('Invalid Web file binding metadata.');
    }
  };
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
      binding.components.forEach(component => {
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
  (Array.isArray(bindings.linearSeqs) ? bindings.linearSeqs : []).forEach(sequence => {
    ['gb', 'gff', 'fasta', 'depth', 'blast'].forEach(slot => value(sequence?.[slot]));
  });
  ['linearComparisons', 'linearCanonicalComparisons'].forEach(slot => {
    (Array.isArray(bindings[slot]) ? bindings[slot] : []).forEach(row => value(row?.file));
  });
  return bindings;
};

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

const decodedResourceBytes = (backing) => {
  if (isEncodedDepthFileEntry(backing.descriptor)) {
    return textToBytes(decodeDepthText(backing.descriptor.data));
  }
  recordStructuralMetric('base64DecodeCount', 1, {
    resourceId: backing.resourceId,
    resourceName: backing.name
  });
  return base64ToBytes(backing.descriptor.data);
};

const materializeBytes = (backing) => {
  if (!backing.bytesPromise) {
    recordStructuralMetric('resourceByteReadCount', 1, {
      resourceId: backing.resourceId,
      resourceName: backing.name
    });
    backing.bytesPromise = Promise.resolve().then(async () => {
      let bytes;
      try {
        bytes = decodedResourceBytes(backing);
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

export const adoptCurrentSessionResources = (resources) => {
  if (!resources || typeof resources !== 'object' || Array.isArray(resources)) {
    throw new Error('Current session resources must be an object.');
  }
  const backings = new Map();
  Object.entries(resources).forEach(([rawResourceId, descriptor]) => {
    const resourceId = String(rawResourceId || '').trim();
    if (!resourceId || backings.has(resourceId)) {
      throw new Error('Current session resource IDs must be unique non-empty strings.');
    }
    backings.set(resourceId, validateDescriptor(resourceId, descriptor));
  });
  const table = Object.freeze({});
  tableBackings.set(table, backings);
  return table;
};

export const adoptedSessionResourceDescriptor = (table, resourceId) => (
  requireBacking(table, resourceId).descriptor
);

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

export const createCombinedSessionResourceFileView = (
  table,
  components,
  metadata = {}
) => {
  const views = (Array.isArray(components) ? components : []).map(component => (
    createSessionResourceFileView(table, component.resourceId, component)
  ));
  const backings = views.map(view => fileViewBackings.get(view));
  if (backings.length === 0) return null;
  if (backings.length === 1) {
    return createSessionResourceFileView(table, backings[0].resourceId, { ...components[0], ...metadata });
  }
  const appendedLineFeeds = backings.reduce((count, backing) => (
    count + Number(
      backing.size > 0
      && base64DecodedLastByte(backing.descriptor.data, backing.size) !== LF_BYTE
    )
  ), 0);
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
      composite.textPromise = composite.readBytes().then((bytes) => bytesToText(bytes));
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

export const isSessionResourceFileView = (file) => Boolean(
  file && (typeof file === 'object' || typeof file === 'function')
  && fileViewBackings.has(file)
);

export const sessionResourceSource = (file) => {
  const backing = fileViewBackings.get(file);
  if (!backing || backing.dirty) return null;
  if (Array.isArray(backing.descriptors)) {
    return Object.freeze({ descriptors: backing.descriptors });
  }
  return Object.freeze({
    resourceId: backing.resourceId,
    descriptor: backing.descriptor,
    readBytes: () => materializeBytes(backing)
  });
};

export const readSessionResourceBytes = (file) => {
  const backing = fileViewBackings.get(file);
  if (!backing) return null;
  if (typeof backing.readBytes === 'function') return backing.readBytes();
  return materializeBytes(backing);
};

export const readSessionResourceText = (file) => {
  const backing = fileViewBackings.get(file);
  if (!backing) return null;
  if (typeof backing.readText === 'function') return backing.readText();
  return materializeText(backing);
};

export const releaseSessionResourceContent = (file) => {
  const backing = fileViewBackings.get(file);
  if (!backing) return false;
  const releaseBacking = (candidate) => {
    candidate.bytesPromise = null;
    candidate.textPromise = null;
    candidate.filePromise = null;
    candidate.sourceBackings?.forEach(releaseBacking);
  };
  releaseBacking(backing);
  return true;
};
