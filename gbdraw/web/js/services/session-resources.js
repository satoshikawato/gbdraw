import {
  base64ToBytes,
  bytesToText,
  bytesToBase64,
  getSessionResourceSource,
  readFileBytes,
  registerLazyFileContent,
  sha256Hex,
  textToBytes
} from './file-content-cache.js';
import {
  decodeDepthText,
  isEncodedDepthFileEntry
} from './depth-file-codec.js';
import { isAdoptedCanonicalSession } from './session-authority.js';

const resourceBytes = (entry) => {
  if (isEncodedDepthFileEntry(entry)) {
    return textToBytes(decodeDepthText(entry.data));
  }
  if (entry?.encoding !== 'base64' || typeof entry?.data !== 'string') {
    throw new Error('Session resources must contain supported embedded file data.');
  }
  return base64ToBytes(entry.data);
};

const adoptedResourceTables = new WeakSet();

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

const validateEncodedResourceEntry = (resourceId, entry) => {
  if (!entry || typeof entry !== 'object' || Array.isArray(entry)) {
    throw new Error(`Session resource ${resourceId} must be an object.`);
  }
  if (
    entry.encoding === 'base64'
      ? typeof entry.data !== 'string'
      : !isEncodedDepthFileEntry(entry)
  ) {
    throw new Error(`Session resource ${resourceId} has an unsupported encoded payload.`);
  }
  if (!Number.isSafeInteger(entry.size) || entry.size < 0) {
    throw new Error(`Session resource ${resourceId} has an invalid declared byte size.`);
  }
  if (
    entry.lastModified !== undefined
    && (!Number.isFinite(entry.lastModified) || entry.lastModified < 0)
  ) {
    throw new Error(`Session resource ${resourceId} has an invalid lastModified value.`);
  }
  return {
    resourceId,
    name: String(entry.name || 'file'),
    type: String(entry.type || 'application/octet-stream'),
    lastModified: Number(entry.lastModified) || 0,
    size: entry.size,
    checksum: normalizedChecksum(entry.checksum),
    entry,
    bytesPromise: null,
    textPromise: null,
    source: Object.freeze({ resourceId, descriptor: entry })
  };
};

const requireResourceBacking = (table, resourceId) => {
  if (!table || !adoptedResourceTables.has(table)) {
    throw new TypeError('An adopted current-session resource table is required.');
  }
  const normalizedId = String(resourceId || '').trim();
  const backing = table.backings.get(normalizedId);
  if (!backing) throw new Error(`Missing canonical resource: ${normalizedId}`);
  return backing;
};

const materializeResourceBytes = (backing) => {
  if (!backing.bytesPromise) {
    backing.bytesPromise = Promise.resolve().then(async () => {
      let bytes;
      try {
        bytes = resourceBytes(backing.entry);
      } catch (error) {
        throw new Error(
          `Session resource ${backing.resourceId} contains invalid encoded data.`,
          { cause: error }
        );
      }
      if (bytes.byteLength !== backing.size) {
        throw new Error(
          `Session resource ${backing.resourceId} byte size does not match its descriptor.`
        );
      }
      if (backing.checksum) {
        const actual = await sha256Hex(bytes);
        if (actual !== backing.checksum) {
          throw new Error(`Session resource ${backing.resourceId} checksum does not match.`);
        }
      }
      return bytes;
    });
  }
  return backing.bytesPromise;
};

const materializeResourceText = (backing) => {
  if (!backing.textPromise) {
    backing.textPromise = materializeResourceBytes(backing).then((bytes) => bytesToText(bytes));
  }
  return backing.textPromise;
};

export const adoptCurrentSessionResources = (resources) => {
  if (!resources || typeof resources !== 'object' || Array.isArray(resources)) {
    throw new Error('Current session resources must be an object.');
  }
  const table = { backings: new Map() };
  Object.entries(resources).forEach(([rawResourceId, entry]) => {
    const resourceId = String(rawResourceId || '').trim();
    if (!resourceId || table.backings.has(resourceId)) {
      throw new Error('Current session resource IDs must be unique non-empty strings.');
    }
    table.backings.set(
      resourceId,
      validateEncodedResourceEntry(resourceId, entry)
    );
  });
  adoptedResourceTables.add(table);
  return table;
};

export const getAdoptedSessionResourceDescriptor = (table, resourceId) => (
  requireResourceBacking(table, resourceId).entry
);

export const createSessionResourceFileView = (
  table,
  resourceId,
  metadata = {}
) => {
  const backing = requireResourceBacking(table, resourceId);
  const view = {
    name: String(metadata.name || backing.name),
    type: metadata.type === undefined ? backing.type : String(metadata.type || ''),
    size: backing.size,
    lastModified: metadata.lastModified === undefined
      ? backing.lastModified
      : Number(metadata.lastModified) || 0
  };
  registerLazyFileContent(view, {
    cacheKey: backing,
    readBytes: () => materializeResourceBytes(backing),
    readText: () => materializeResourceText(backing),
    sessionResource: backing.source
  });
  return Object.freeze(view);
};

export const createCombinedSessionResourceFileView = (
  table,
  resourceIds,
  metadata = {}
) => {
  const backings = (Array.isArray(resourceIds) ? resourceIds : [])
    .map((resourceId) => requireResourceBacking(table, resourceId));
  if (backings.length === 0) return null;
  if (backings.length === 1) {
    return createSessionResourceFileView(table, backings[0].resourceId, metadata);
  }
  const composite = {
    bytesPromise: null,
    textPromise: null
  };
  const readBytes = () => {
    if (!composite.bytesPromise) {
      composite.bytesPromise = Promise.all(
        backings.map((backing) => materializeResourceBytes(backing))
      ).then((parts) => {
        const lengths = parts.map((part) => (
          part.byteLength > 0 && part[part.byteLength - 1] !== 0x0A
            ? part.byteLength + 1
            : part.byteLength
        ));
        const combined = new Uint8Array(lengths.reduce((total, length) => total + length, 0));
        let offset = 0;
        parts.forEach((part, index) => {
          combined.set(part, offset);
          offset += part.byteLength;
          if (lengths[index] > part.byteLength) combined[offset++] = 0x0A;
        });
        return combined;
      });
    }
    return composite.bytesPromise;
  };
  const readText = () => {
    if (!composite.textPromise) composite.textPromise = readBytes().then(bytesToText);
    return composite.textPromise;
  };
  const declaredSize = backings.reduce((total, backing) => total + backing.size, 0)
    + backings.slice(0, -1).filter((backing) => backing.size > 0).length;
  const view = {
    name: String(metadata.name || 'canonical-circular-records.gb'),
    type: String(metadata.type || 'text/plain'),
    size: declaredSize,
    lastModified: metadata.lastModified === undefined
      ? Math.max(0, ...backings.map((backing) => backing.lastModified))
      : Number(metadata.lastModified) || 0
  };
  registerLazyFileContent(view, {
    cacheKey: composite,
    readBytes,
    readText,
    sessionResource: Object.freeze({
      descriptors: backings.map((backing) => backing.source)
    })
  });
  return Object.freeze(view);
};

const safeResourceLeaf = (value) => {
  const basename = String(value || 'resource.dat')
    .replace(/\\/g, '/')
    .split('/')
    .pop();
  const safe = basename
    .replace(/[^A-Za-z0-9._-]+/g, '_')
    .replace(/^[._]+|[._]+$/g, '');
  return safe || 'resource.dat';
};

const bindingForFile = (file, resourceId) => ({
  resourceId,
  name: String(file?.name || 'file'),
  type: String(file?.type || ''),
  lastModified: Number(file?.lastModified) || 0
});

const RESOURCE_REFERENCE_FIELDS = new Set([
  'resourceId',
  'gffResourceId',
  'fastaResourceId'
]);

const rewriteResourceRefs = (value, aliases) => {
  if (Array.isArray(value)) {
    return value.map((item) => rewriteResourceRefs(item, aliases));
  }
  if (!value || typeof value !== 'object') return value;
  return Object.fromEntries(
    Object.entries(value).map(([key, item]) => [
      key,
      (
        RESOURCE_REFERENCE_FIELDS.has(key)
        && typeof item === 'string'
        && aliases.has(item)
      )
        ? aliases.get(item)
        : rewriteResourceRefs(item, aliases)
    ])
  );
};

const collectResourceRefs = (value, target = new Set()) => {
  if (Array.isArray(value)) {
    value.forEach((item) => collectResourceRefs(item, target));
    return target;
  }
  if (!value || typeof value !== 'object') return target;
  Object.entries(value).forEach(([key, item]) => {
    if (
      RESOURCE_REFERENCE_FIELDS.has(key)
      && typeof item === 'string'
      && item.trim()
    ) {
      target.add(item.trim());
      return;
    }
    collectResourceRefs(item, target);
  });
  return target;
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

export const buildSessionResources = async (state, committedRequest) => {
  if (
    !committedRequest
    || typeof committedRequest !== 'object'
    || !committedRequest.renderRequest
    || !committedRequest.resources
  ) {
    throw new Error('A committed canonical render request is required to save a session.');
  }

  const resources = {};
  const aliases = new Map();
  const encodedDataToResourceId = new Map();
  const reuseEncodedResources = isAdoptedCanonicalSession(committedRequest);
  let nextResourceNumber = 1;

  const nextResourceId = () => {
    let resourceId;
    do {
      resourceId = `resource-${String(nextResourceNumber).padStart(4, '0')}`;
      nextResourceNumber += 1;
    } while (Object.prototype.hasOwnProperty.call(resources, resourceId));
    return resourceId;
  };

  const adoptEncodedResource = (resourceId, entry) => {
    const normalizedId = String(resourceId || '').trim();
    if (!normalizedId || !entry) {
      throw new Error('An adopted session resource requires an ID and descriptor.');
    }
    const existing = resources[normalizedId];
    if (existing && existing !== entry) {
      const samePayload = (
        existing.encoding === entry.encoding
        && existing.data === entry.data
        && Number(existing.size) === Number(entry.size)
      );
      if (!samePayload) {
        throw new Error(`Conflicting adopted session resource: ${normalizedId}.`);
      }
      resources[normalizedId] = entry;
    } else {
      resources[normalizedId] = entry;
    }
    aliases.set(normalizedId, normalizedId);
    if (entry.encoding === 'base64' && typeof entry.data === 'string') {
      encodedDataToResourceId.set(entry.data, normalizedId);
    }
    return normalizedId;
  };

  const addBytes = async (bytes, metadata = {}) => {
    const encoded = bytesToBase64(bytes);
    const existing = encodedDataToResourceId.get(encoded);
    if (existing) return existing;

    const resourceId = nextResourceId();
    encodedDataToResourceId.set(encoded, resourceId);
    resources[resourceId] = {
      kind: String(metadata.kind || 'web-file'),
      name: `${resourceId}-${safeResourceLeaf(metadata.name)}`,
      type: String(metadata.type || 'application/octet-stream'),
      size: bytes.byteLength,
      lastModified: Number(metadata.lastModified) || 0,
      encoding: 'base64',
      data: encoded
    };
    return resourceId;
  };

  const committedResourceIds = collectResourceRefs(committedRequest.renderRequest);
  for (const resourceId of committedResourceIds) {
    const entry = committedRequest.resources[resourceId];
    if (!entry) {
      throw new Error(`Committed render resource is missing: ${resourceId}.`);
    }
    const nextId = reuseEncodedResources
      ? adoptEncodedResource(resourceId, entry)
      : await addBytes(resourceBytes(entry), entry);
    aliases.set(resourceId, nextId);
  }

  const bindFile = async (file) => {
    if (!file) return null;
    const source = getSessionResourceSource(file);
    if (reuseEncodedResources && source?.resourceId && source?.descriptor) {
      const resourceId = adoptEncodedResource(source.resourceId, source.descriptor);
      return bindingForFile(file, resourceId);
    }
    if (reuseEncodedResources && Array.isArray(source?.descriptors)) {
      source.descriptors.forEach(({ resourceId, descriptor }) => {
        adoptEncodedResource(resourceId, descriptor);
      });
      return undefined;
    }
    const bytes = await readFileBytes(file);
    const resourceId = await addBytes(bytes, {
      kind: 'web-file',
      name: file.name,
      type: file.type,
      lastModified: file.lastModified
    });
    return bindingForFile(file, resourceId);
  };

  const bindFileValue = async (value) => {
    if (Array.isArray(value)) {
      return Promise.all(value.map((item) => bindFileValue(item)));
    }
    return bindFile(value);
  };

  const files = state?.files || {};
  const linearSeqs = Array.isArray(state?.linearSeqs) ? state.linearSeqs : [];
  const linearComparisons = Array.isArray(state?.linearComparisonPlan?.edges)
    ? state.linearComparisonPlan.edges
    : [];

  const bindings = {
    schema: 1,
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

  return {
    renderRequest: identityAliases
      ? committedRequest.renderRequest
      : rewriteResourceRefs(committedRequest.renderRequest, aliases),
    resources,
    webFiles: {
      ...rewriteOriginalNameHints(committedRequest.webFiles, aliases),
      bindings
    }
  };
};
