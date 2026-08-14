import {
  base64ToBytes,
  bytesToBase64,
  getSessionResourceSource,
  readFileBytes,
  sha256Hex,
  textToBytes
} from './file-content-cache.js';
import {
  decodeDepthText,
  isEncodedDepthFileEntry
} from './depth-file-codec.js';
import { isAdoptedCanonicalSession } from './session-authority.js';
import { recordStructuralMetric } from './runtime-test-hooks.js';

const resourceBytes = (entry) => {
  if (isEncodedDepthFileEntry(entry)) {
    return textToBytes(decodeDepthText(entry.data));
  }
  if (entry?.encoding !== 'base64' || typeof entry?.data !== 'string') {
    throw new Error('Session resources must contain supported embedded file data.');
  }
  return base64ToBytes(entry.data);
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
  const identityToResourceId = new Map();
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

  const adoptEncodedResource = (resourceId, descriptor) => {
    const normalizedId = String(resourceId || '').trim();
    if (!normalizedId || !descriptor) {
      throw new Error('An adopted session resource requires an ID and descriptor.');
    }
    const existing = resources[normalizedId];
    if (existing && existing !== descriptor) {
      const samePayload = (
        existing.encoding === descriptor.encoding
        && existing.data === descriptor.data
        && Number(existing.size) === Number(descriptor.size)
      );
      if (!samePayload) {
        throw new Error(`Conflicting adopted session resource: ${normalizedId}.`);
      }
    }
    resources[normalizedId] = descriptor;
    aliases.set(normalizedId, normalizedId);
    return normalizedId;
  };

  const addBytes = async (bytes, metadata = {}) => {
    const identity = `${bytes.byteLength}:${await sha256Hex(bytes)}`;
    const existing = identityToResourceId.get(identity);
    if (existing) return existing;

    recordStructuralMetric('base64EncodeCount', 1, {
      resourceName: String(metadata.name || 'file')
    });
    recordStructuralMetric('encodedByteCount', bytes.byteLength, {
      resourceName: String(metadata.name || 'file')
    });
    const encoded = bytesToBase64(bytes);
    const resourceId = nextResourceId();
    identityToResourceId.set(identity, resourceId);
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

  if (reuseEncodedResources) {
    Object.entries(committedRequest.resources).forEach(([resourceId, descriptor]) => {
      adoptEncodedResource(resourceId, descriptor);
    });
  }

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
