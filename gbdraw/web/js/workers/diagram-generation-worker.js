import { sendBoundedJson } from '../services/bounded-json-transport.js';
import { PYTHON_HELPERS } from '../app/python-helpers.js';
import {
  DIAGRAM_HELPER_OPERATIONS,
  readFeatureCatalogRevision,
  readHeldSegments,
  readSvgDocumentParts
} from '../services/diagram-worker-protocol.js';
import { normalizeUserFacingError } from '../utils/error-normalization.js';
import { ORTHOGROUP_INFERENCE_VERSION } from '../contracts/schema.generated.js';
import { createWorkerFontRegistry } from './worker-font-registry.js';

let runtimePromise = null;
let runtime = null;
let operationQueue = Promise.resolve();
let auxiliaryAcknowledgement = null;
const sendAuxiliaryResult = async (type, requestId, result) => {
  let wholeReply = false;
  try {
    await sendBoundedJson(result, (part, transfers = []) => {
      if (part.whole) {
        wholeReply = true;
        self.postMessage({ type, requestId, ok: true, result: part.value });
        return Promise.resolve();
      }
      return new Promise((resolve) => {
        auxiliaryAcknowledgement = { requestId, resolve };
        self.postMessage({ type, requestId, status: 'part', ...part }, transfers);
      });
    }, { consume: true });
    if (!wholeReply) self.postMessage({ type, requestId, ok: true });
  } finally {
    auxiliaryAcknowledgement = null;
  }
};
const RENDER_RESOURCE_CACHE = '/gbdraw-web-render-resource-cache';
const RENDER_WORKSPACE_MARKER = '.gbdraw-worker-render-workspace';
const RENDER_RESOURCE_ID_RE = /^[a-z][a-z0-9]*(?:-[a-z0-9]+)*$/;
const RENDER_RESOURCE_TOKEN_RE = /^render-resource-[1-9][0-9]*$/;
const cachedRenderResources = new Map();

export const buildPreparedResourceIdentityMap = (resourceManifest) => {
  if (!Array.isArray(resourceManifest)) {
    throw new TypeError('Prepared resource identities require a resource manifest.');
  }
  const identities = {};
  resourceManifest.forEach((metadata) => {
    const resourceId = String(metadata?.resourceId || '').trim();
    const cacheToken = String(metadata?.cacheToken || '');
    const size = Number(metadata?.size);
    if (!RENDER_RESOURCE_ID_RE.test(resourceId) || resourceId in identities) {
      throw new TypeError(`Invalid or duplicate prepared render resource '${resourceId}'.`);
    }
    if (!RENDER_RESOURCE_TOKEN_RE.test(cacheToken)) {
      throw new TypeError(`Invalid cache token for prepared render resource '${resourceId}'.`);
    }
    if (!Number.isSafeInteger(size) || size < 0) {
      throw new TypeError(`Invalid byte size for prepared render resource '${resourceId}'.`);
    }
    identities[resourceId] = { cacheToken, size };
  });
  return identities;
};

const emitTestLifecycle = (enabled, requestId, name, detail = {}) => {
  if (!enabled) return;
  self.postMessage({
    type: 'test-lifecycle',
    requestId,
    event: {
      name,
      timestamp: globalThis.performance?.now?.() ?? Date.now(),
      ...detail
    }
  });
};

export const collectGenerationResultTransferList = (payload) => {
  const buffers = new Set();
  const addBuffer = (value) => {
    if (value instanceof ArrayBuffer) buffers.add(value);
    else if (ArrayBuffer.isView(value)) buffers.add(value.buffer);
  };
  (Array.isArray(payload?.results) ? payload.results : []).forEach((result) => {
    addBuffer(result?.segments);
  });
  addBuffer(payload?.metadata);
  addBuffer(payload?.featureCatalog);
  return Array.from(buffers);
};

const bytesForDigest = (value) => {
  if (value instanceof ArrayBuffer) return new Uint8Array(value);
  if (ArrayBuffer.isView(value)) {
    return new Uint8Array(value.buffer, value.byteOffset, value.byteLength);
  }
  return null;
};

const hexDigest = async (bytes) => {
  const digest = new Uint8Array(await globalThis.crypto.subtle.digest('SHA-256', bytes));
  return Array.from(digest, (value) => value.toString(16).padStart(2, '0')).join('');
};

const invalidGenerationReply = () => serializeError({
  code: 'RESULT_INVALID', operation: 'generate', stage: 'result-admission'
});

/**
 * Check one Result's document against its `segments` buffer and the segments
 * the request named, and describe the whole document: its UTF-8 length, and a
 * manifest that names each children text by revision, held or sent.
 * @param {{ document?: unknown, segments?: unknown }} result
 * @param {ReadonlyMap<string, number>} heldBytes The request's held segments: revision -> bytes.
 * @param {TextEncoder} encoder
 */
const describeResultDocument = async (result, heldBytes, encoder) => {
  const segments = result.segments === undefined ? new Uint8Array(0) : bytesForDigest(result.segments);
  if (!segments) throw invalidGenerationReply();
  const parts = readSvgDocumentParts(result.document, segments.byteLength, invalidGenerationReply);
  let bytes = encoder.encode(parts.head).byteLength + encoder.encode(parts.tail).byteLength;
  const pieces = [];
  for (const piece of parts.pieces) {
    if (piece.kind === 'raw') {
      bytes += encoder.encode(piece.text).byteLength;
      pieces.push(piece);
      continue;
    }
    const { start, revision } = piece;
    if (piece.offset === undefined) {
      if (heldBytes.get(revision) !== piece.bytes) throw invalidGenerationReply();
    } else if (await hexDigest(segments.subarray(piece.offset, piece.offset + piece.bytes)) !== revision) {
      throw invalidGenerationReply();
    }
    bytes += encoder.encode(start).byteLength + '</g>'.length + piece.bytes;
    pieces.push({ kind: 'content', start, revision, bytes: piece.bytes });
  }
  return { bytes, manifest: { head: parts.head, tail: parts.tail, pieces } };
};

/**
 * Build a private identity from the buffers that are already being transferred.
 * SVG children, metadata and feature catalog payloads are never re-encoded or
 * concatenated here. Each sent children text is hashed and must equal its
 * revision; a held one is named by the revision the request held, so a
 * document has one fingerprint and one byte count whether its children were
 * held or sent. `featureCatalogDigest` is the digest of the catalog bytes sent,
 * or of the held catalog's bytes when Python reported it unchanged, so equal
 * artifacts get equal fingerprints either way.
 * @param {{ results?: unknown, metadata?: unknown }} payload The wrapper's reply.
 * @param {{ bytes: number, sha256: string } | null} [featureCatalogDigest]
 * @param {readonly { revision: string, bytes: number }[]} [heldSegments] The segments the request named.
 */
export const buildGeneratedArtifactTransportIdentity = async (
  payload,
  featureCatalogDigest = null,
  heldSegments = []
) => {
  const results = Array.isArray(payload?.results) ? payload.results : [];
  const encoder = new TextEncoder();
  const heldBytes = new Map(heldSegments.map(({ revision, bytes }) => [revision, bytes]));
  const resultManifest = [];
  let resultBytes = 0;
  let resultNameBytes = 0;

  for (let index = 0; index < results.length; index += 1) {
    const result = results[index] || {};
    const nameBytes = encoder.encode(String(result.name || ''));
    const document = await describeResultDocument(result, heldBytes, encoder);
    resultBytes += document.bytes;
    resultNameBytes += nameBytes.byteLength;
    resultManifest.push({
      index,
      nameBytes: nameBytes.byteLength,
      nameSha256: await hexDigest(nameBytes),
      documentBytes: document.bytes,
      documentSha256: await hexDigest(encoder.encode(JSON.stringify(document.manifest)))
    });
  }

  const metadataBytes = bytesForDigest(payload?.metadata);
  if (!metadataBytes) {
    throw new Error('Generated metadata was not encoded for Worker transport.');
  }
  const manifest = {
    schema: 3,
    results: resultManifest,
    metadata: {
      bytes: metadataBytes.byteLength,
      sha256: await hexDigest(metadataBytes)
    },
    featureCatalog: featureCatalogDigest
      ? { bytes: featureCatalogDigest.bytes, sha256: featureCatalogDigest.sha256 }
      : null
  };
  const catalogBytes = featureCatalogDigest?.bytes ?? 0;
  const fingerprint = await hexDigest(encoder.encode(JSON.stringify(manifest)));
  return Object.freeze({
    schema: 1,
    algorithm: 'SHA-256',
    fingerprint,
    resultBytes,
    resultNameBytes,
    metadataBytes: metadataBytes.byteLength + catalogBytes,
    // Parsed metadata retains strings and object structure in addition to the
    // transferred UTF-8 representation. Two bytes per source byte is a bounded,
    // deliberately conservative owner estimate; Result UTF-16 strings are added
    // separately by the main-thread handle. `resultBytes` counts the whole
    // admitted document, held or sent. The feature catalog counts as metadata
    // whether it was sent or kept (D-38 OQ-2).
    retainedBytes: resultBytes + resultNameBytes + (metadataBytes.byteLength + catalogBytes) * 2
  });
};

/**
 * Name the feature catalog of a Generate reply and fix its artifact identity.
 * A sent catalog is hashed; an "unchanged" reply reuses the digest of the held
 * catalog (`held`, read from the run message), so nothing is hashed again. The
 * page checks that "unchanged" names the revision it holds.
 * @param {{ featureCatalog?: unknown, featureCatalogRevision?: unknown, artifactIdentity?: unknown }} outcome
 * @param {{ revision: string, sha256: string, bytes: number } | null} held
 * @param {readonly { revision: string, bytes: number }[]} heldSegments The segments the request named.
 */
export const resolveGeneratedFeatureCatalog = async (outcome, held, heldSegments) => {
  const catalogBytes = bytesForDigest(outcome.featureCatalog);
  const digest = catalogBytes
    ? { bytes: catalogBytes.byteLength, sha256: await hexDigest(catalogBytes) }
    : held
      ? { bytes: held.bytes, sha256: held.sha256 }
      : null;
  const pythonRevision = typeof outcome.featureCatalogRevision === 'string' && outcome.featureCatalogRevision
    ? outcome.featureCatalogRevision
    : null;
  outcome.featureCatalogRevision = pythonRevision && digest ? { revision: pythonRevision, ...digest } : null;
  outcome.artifactIdentity = await buildGeneratedArtifactTransportIdentity(outcome, digest, heldSegments);
  return outcome;
};

export const serializeError = (error, options = {}) => normalizeUserFacingError(error || {}, options);

const attachCleanupDiagnostic = (primary, cleanupError) => {
  if (!primary || typeof primary !== 'object' || !cleanupError) return primary;
  primary.secondary = [...(Array.isArray(primary.secondary) ? primary.secondary : []),
    { code: 'CLEANUP_FAILED', stage: 'cleanup' }].slice(0, 2);
  return primary;
};

export const resolveGenerationCleanupOutcome = ({
  result = null,
  primaryError = null,
  destroyError = null,
  workspaceError = null
} = {}) => {
  const pythonError = (
    result?.error &&
    typeof result.error === 'object' &&
    !Array.isArray(result.error)
  )
    ? result.error
    : null;
  const primary = primaryError || pythonError || destroyError || workspaceError;
  if (!primaryError && !pythonError && primary) {
    Object.assign(primary, { code: 'CLEANUP_FAILED', stage: 'cleanup' });
  }
  if (primary) {
    if (primary !== destroyError && destroyError) {
      attachCleanupDiagnostic(
        primary,
        destroyError
      );
    }
    if (primary !== workspaceError && workspaceError) {
      attachCleanupDiagnostic(
        primary,
        workspaceError
      );
    }
  }
  if (primaryError) throw primaryError;
  if (pythonError) return result;
  if (destroyError) throw destroyError;
  if (workspaceError) throw workspaceError;
  return result;
};

const ensureLocalAsset = async (url, label) => {
  const response = await fetch(url, { method: 'HEAD', cache: 'no-store' });
  if (!response.ok) {
    throw new Error(`Missing packaged asset: ${label} (${response.status}) at ${url}`);
  }
  return url;
};

const readRuntimeCapabilities = (pyodide) => {
  const raw = pyodide.runPython(`
import json as _json
from gbdraw.api import get_web_runtime_capabilities as _get_web_runtime_capabilities
_json.dumps(_get_web_runtime_capabilities(), sort_keys=True)
  `);
  return JSON.parse(String(raw));
};

// The first-use preparation steps, reported as they start: a real step
// count, not a time estimate.
const RUNTIME_PREPARATION_STEPS = Object.freeze([
  'loading Python', 'loading the package installer', 'checking packages', 'installing packages', 'starting gbdraw'
]);

const initializeRuntime = async ({
  id = null,
  pyodideIndexUrl,
  pyodideModuleUrl,
  localWheelUrls = [],
  gbdrawWheelUrl
} = {}) => {
  if (runtime) return runtime;
  if (runtimePromise) return runtimePromise;
  if (!pyodideIndexUrl || !pyodideModuleUrl) {
    throw new Error('Diagram generation worker requires Pyodide asset URLs.');
  }
  if (!gbdrawWheelUrl) {
    throw new Error('Diagram generation worker requires the gbdraw wheel URL.');
  }

  /** @param {number} index */
  const reportStep = (index) => self.postMessage({
    id, type: 'init-progress', step: index + 1, total: RUNTIME_PREPARATION_STEPS.length, label: RUNTIME_PREPARATION_STEPS[index]
  });
  runtimePromise = (async () => {
    reportStep(0);
    const { loadPyodide } = await import(pyodideModuleUrl);
    const pyodide = await loadPyodide({
      indexURL: pyodideIndexUrl,
      packageBaseUrl: pyodideIndexUrl
    });
    reportStep(1);
    await pyodide.loadPackage('micropip');
    const micropip = pyodide.pyimport('micropip');

    reportStep(2);
    await Promise.all(
      localWheelUrls.map((url, index) =>
        ensureLocalAsset(url, `Pyodide dependency wheel #${index + 1}`)
      )
    );
    await ensureLocalAsset(gbdrawWheelUrl, 'gbdraw browser wheel');
    reportStep(3);
    await micropip.install(localWheelUrls);
    await micropip.install(gbdrawWheelUrl);
    reportStep(4);
    await pyodide.runPythonAsync(PYTHON_HELPERS);

    const capabilities = readRuntimeCapabilities(pyodide);
    runtime = { pyodide, capabilities };
    return runtime;
  })();

  try {
    return await runtimePromise;
  } catch (error) {
    runtimePromise = null;
    throw error;
  }
};

// The fonts registered in this Worker's Python (workers/worker-font-registry.js).
let fontRegistry = null;
const fontRegistryFor = (pyodide) => {
  if (!fontRegistry) {
    fontRegistry = createWorkerFontRegistry({
      fs: pyodide.FS,
      registerInPython: (entries) => callJsonHelper(pyodide, 'register_web_fonts_json', [JSON.stringify(entries)])
    });
  }
  return fontRegistry;
};

const HELPER_FILE_NAMES = Object.freeze({
  source: 'source.bin',
  gff: 'source.gff',
  fasta: 'source.fasta',
  pairs: 'pairs.json',
  visibility: 'feature-visibility.tsv',
  featureOverrides: 'feature-overrides.tsv',
  rawTsv: 'raw-losatp.tsv'
});

const requirePayloadObject = (payload, operation) => {
  if (!payload || typeof payload !== 'object' || Array.isArray(payload)) {
    throw new TypeError(`Diagram helper '${operation}' requires an object payload.`);
  }
  return payload;
};

const assertAllowedPayloadKeys = (payload, operation, allowedKeys) => {
  const allowed = new Set(allowedKeys);
  Object.keys(payload).forEach((key) => {
    if (!allowed.has(key)) {
      throw new TypeError(`Diagram helper '${operation}' does not accept '${key}'.`);
    }
  });
};

const removeWorkspace = (pyodide, path) => {
  if (!pyodide.FS.analyzePath(path).exists) return;
  pyodide.FS.readdir(path).forEach((entry) => {
    if (entry === '.' || entry === '..') return;
    const child = `${path}/${entry}`;
    const stat = pyodide.FS.stat(child);
    if (pyodide.FS.isDir(stat.mode)) removeWorkspace(pyodide, child);
    else pyodide.FS.unlink(child);
  });
  pyodide.FS.rmdir(path);
};

const withRequestWorkspace = async (pyodide, requestId, kind, callback) => {
  const workspace = `/gbdraw-web-${kind}-${Number(requestId) || 0}`;
  pyodide.FS.mkdir(workspace);
  let primaryError = null;
  try {
    return await callback(workspace);
  } catch (error) {
    primaryError = error;
    throw error;
  } finally {
    try {
      removeWorkspace(pyodide, workspace);
    } catch (cleanupError) {
      if (primaryError) {
        attachCleanupDiagnostic(
          primaryError,
          cleanupError
        );
      } else {
        throw { code: 'CLEANUP_FAILED', stage: 'cleanup' };
      }
    }
  }
};

const stageHelperFiles = (pyodide, workspace, files, allowedRoles) => {
  const allowed = new Set(allowedRoles);
  const paths = new Map();
  const entries = Array.isArray(files) ? files : [];
  entries.forEach((file) => {
    const role = String(file?.role || '').trim();
    if (!allowed.has(role)) {
      throw new TypeError(`Unexpected diagram helper file role '${role || '(blank)'}'.`);
    }
    if (paths.has(role)) {
      throw new TypeError(`Duplicate diagram helper file role '${role}'.`);
    }
    if (!(file?.bytes instanceof ArrayBuffer)) {
      throw new TypeError(`Diagram helper file '${role}' requires an ArrayBuffer.`);
    }
    const path = `${workspace}/${HELPER_FILE_NAMES[role]}`;
    pyodide.FS.writeFile(path, new Uint8Array(file.bytes));
    paths.set(role, path);
  });
  return paths;
};

const ensureRenderResourceCache = (pyodide) => {
  if (!pyodide.FS.analyzePath(RENDER_RESOURCE_CACHE).exists) {
    pyodide.FS.mkdir(RENDER_RESOURCE_CACHE);
  }
};

const stageRenderResources = async (
  pyodide,
  workspace,
  resourceManifest,
  stagedResources,
  testLifecycleEnabled,
  requestId
) => {
  if (!Array.isArray(resourceManifest) || !Array.isArray(stagedResources)) {
    throw new TypeError('Diagram generation requires a staged resource manifest.');
  }
  ensureRenderResourceCache(pyodide);
  const stagedById = new Map();
  stagedResources.forEach((entry) => {
    const resourceId = String(entry?.resourceId || '').trim();
    if (!RENDER_RESOURCE_ID_RE.test(resourceId) || stagedById.has(resourceId)) {
      throw new TypeError(`Invalid or duplicate staged render resource '${resourceId}'.`);
    }
    if (!RENDER_RESOURCE_TOKEN_RE.test(String(entry?.cacheToken || ''))) {
      throw new TypeError(`Invalid cache token for staged render resource '${resourceId}'.`);
    }
    if (!(entry?.bytes instanceof ArrayBuffer)) {
      throw new TypeError(`Staged render resource '${resourceId}' requires an ArrayBuffer.`);
    }
    stagedById.set(resourceId, entry);
  });

  const manifestIds = new Set();
  const resourcePaths = {};
  const resourcesDirectory = `${workspace}/resources`;
  pyodide.FS.mkdir(workspace);
  pyodide.FS.writeFile(`${workspace}/${RENDER_WORKSPACE_MARKER}`, new Uint8Array(0));
  pyodide.FS.mkdir(resourcesDirectory);
  for (let index = 0; index < resourceManifest.length; index += 1) {
    const metadata = resourceManifest[index];
    const resourceId = String(metadata?.resourceId || '').trim();
    const cacheToken = String(metadata?.cacheToken || '');
    const size = Number(metadata?.size);
    if (!RENDER_RESOURCE_ID_RE.test(resourceId) || manifestIds.has(resourceId)) {
      throw new TypeError(`Invalid or duplicate render resource '${resourceId}'.`);
    }
    if (!RENDER_RESOURCE_TOKEN_RE.test(cacheToken)) {
      throw new TypeError(`Invalid cache token for render resource '${resourceId}'.`);
    }
    if (!Number.isSafeInteger(size) || size < 0) {
      throw new TypeError(`Invalid byte size for render resource '${resourceId}'.`);
    }
    manifestIds.add(resourceId);
    let cached = cachedRenderResources.get(resourceId);
    const staged = stagedById.get(resourceId);
    if (staged) {
      if (staged.cacheToken !== cacheToken || staged.bytes.byteLength !== size) {
        throw new TypeError(`Staged render resource '${resourceId}' does not match its manifest.`);
      }
      emitTestLifecycle(testLifecycleEnabled, requestId, 'worker-resource-stage-start', {
        resourceId,
        bytes: size
      });
      if (cached && pyodide.FS.analyzePath(cached.path).exists) {
        pyodide.FS.unlink(cached.path);
      }
      const cachePath = `${RENDER_RESOURCE_CACHE}/${cacheToken}.bin`;
      if (pyodide.FS.analyzePath(cachePath).exists) pyodide.FS.unlink(cachePath);
      if (typeof pyodide.FS.createDataFile !== 'function') {
        throw new Error('Pyodide render resource ownership transfer is unavailable.');
      }
      pyodide.FS.createDataFile(
        RENDER_RESOURCE_CACHE,
        `${cacheToken}.bin`,
        new Uint8Array(staged.bytes),
        true,
        true,
        true
      );
      cached = { cacheToken, path: cachePath, size };
      cachedRenderResources.set(resourceId, cached);
      staged.bytes = null;
      stagedById.delete(resourceId);
      emitTestLifecycle(testLifecycleEnabled, requestId, 'worker-resource-stage-end', {
        resourceId,
        bytes: size
      });
    }
    if (
      !cached
      || cached.cacheToken !== cacheToken
      || cached.size !== size
      || !pyodide.FS.analyzePath(cached.path).exists
    ) {
      throw new Error(`Render resource cache miss for '${resourceId}'.`);
    }
    const requestPath = `${resourcesDirectory}/${String(index + 1).padStart(4, '0')}.bin`;
    if (typeof pyodide.FS.symlink !== 'function') {
      throw new Error('Pyodide render resource symbolic links are unavailable.');
    }
    pyodide.FS.symlink(cached.path, requestPath);
    resourcePaths[resourceId] = requestPath;
  }
  stagedById.forEach((_entry, resourceId) => {
    if (!manifestIds.has(resourceId)) {
      throw new TypeError(`Unexpected staged render resource '${resourceId}'.`);
    }
  });
  cachedRenderResources.forEach((cached, resourceId) => {
    if (manifestIds.has(resourceId)) return;
    if (pyodide.FS.analyzePath(cached.path).exists) pyodide.FS.unlink(cached.path);
    cachedRenderResources.delete(resourceId);
  });
  stagedResources.splice(0);
  return resourcePaths;
};

// Rendering and source-bound helpers enter the same canonical staging boundary.
const prepareCanonicalResources = (
  pyodide, workspace, { resourceManifest, stagedResources },
  { testLifecycleEnabled = false, requestId = 0 } = {}
) => stageRenderResources(
  pyodide, workspace, resourceManifest, stagedResources, testLifecycleEnabled, requestId
);

const requireHelperFile = (paths, role, operation) => {
  const path = paths.get(role);
  if (!path) throw new TypeError(`Diagram helper '${operation}' requires the '${role}' file.`);
  return path;
};

export const callJsonHelper = (pyodide, helperName, args) => {
  const helper = pyodide.globals.get('call_web_json_helper');
  let result = null;
  let primary = null;
  try {
    if (typeof helper !== 'function') throw { code: 'HELPER_PROTOCOL', stage: 'helper' };
    try {
      result = JSON.parse(String(helper(helperName, ...args) || 'null'));
    } catch (error) {
      throw normalizeUserFacingError(error, { code: 'RESULT_INVALID', stage: 'helper' });
    }
    if (result?.error) throw normalizeUserFacingError(result.error);
  } catch (error) {
    primary = serializeError(error, { stage: 'helper' });
  }
  try {
    helper?.destroy?.();
  } catch (error) {
    if (primary) attachCleanupDiagnostic(primary, error);
    else primary = { code: 'CLEANUP_FAILED', stage: 'cleanup' };
  }
  if (primary) throw primary;
  return result;
};

const jsonArgument = (value, fallback) => JSON.stringify(value ?? fallback);

const HELPER_OPERATION_SPECS = Object.freeze({
  [DIAGRAM_HELPER_OPERATIONS.RESOLVE_SIMILARITY_ALIGNMENT]: {
    keys: ['request', 'projection', 'resourceManifest', 'stagedResources'],
    fileRoles: [],
    run: async (pyodide, payload, _paths, _operation, workspace) => {
      let resourcePaths = {};
      if (payload.projection !== undefined) {
        resourcePaths = await prepareCanonicalResources(
          pyodide, `${workspace}/projection`, payload
        );
      } else if (payload.resourceManifest !== undefined || payload.stagedResources !== undefined) {
        throw new TypeError('Alignment resources require a projection context.');
      }
      return callJsonHelper(pyodide, 'resolve_similarity_alignment_json', [
        jsonArgument(payload.request, null), jsonArgument(payload.projection, null),
        jsonArgument(resourcePaths, {}), workspace
      ]);
    }
  },
  // Load Feature Edits TSV: Python reads the table against the records of the
  // committed request, staged as Generate stages them (design Q4 6.4, R4).
  [DIAGRAM_HELPER_OPERATIONS.READ_FEATURE_OVERRIDE_TABLE]: {
    keys: ['files', 'canonicalRequest', 'resourceManifest', 'stagedResources'],
    fileRoles: ['featureOverrides'],
    run: async (pyodide, payload, paths, operation, workspace) => callJsonHelper(
      pyodide, 'read_feature_override_table_json', [
        requireHelperFile(paths, 'featureOverrides', operation),
        jsonArgument(payload.canonicalRequest, null),
        jsonArgument(await prepareCanonicalResources(pyodide, `${workspace}/request`, payload), {}),
        `${workspace}/output`
      ]
    )
  },
  [DIAGRAM_HELPER_OPERATIONS.EVALUATE_RULES]: {
    keys: ['features', 'rules', 'kind'],
    fileRoles: [],
    run: (pyodide, payload) => callJsonHelper(pyodide, 'evaluate_rules_json', [
      jsonArgument(payload.features, []), jsonArgument(payload.rules, []), String(payload.kind || 'color')
    ])
  },
  // One row of the label fit dialog for a figure drawn with crossing leaders:
  // Python searches it on the figure it kept (OD-F14); nothing here places labels.
  [DIAGRAM_HELPER_OPERATIONS.COMPUTE_LABEL_FIT_ROW]: {
    keys: ['rowToken', 'kind'],
    fileRoles: [],
    run: (pyodide, payload) => callJsonHelper(pyodide, 'compute_label_fit_row_json', [
      String(payload.rowToken || ''), String(payload.kind || '')
    ])
  },
  [DIAGRAM_HELPER_OPERATIONS.REGISTER_FONTS]: {
    keys: ['fonts'],
    fileRoles: [],
    run: (pyodide, payload) => fontRegistryFor(pyodide).register(payload.fonts)
  },
  // The fonts an export carries (custom fonts D-05, D-13): the embedded-font
  // block of an SVG (`kind` svg or raster), and the TrueType faces of a PDF's
  // text runs (fonts by file name).
  [DIAGRAM_HELPER_OPERATIONS.EMBED_FONTS]: {
    keys: ['svg', 'kind'],
    fileRoles: [],
    run: (pyodide, payload) => callJsonHelper(pyodide, 'svg_font_embedding_json', [
      String(payload.svg || ''), String(payload.kind || 'svg')
    ])
  },
  [DIAGRAM_HELPER_OPERATIONS.PDF_FONTS]: {
    keys: ['runs', 'known'],
    fileRoles: [],
    run: (pyodide, payload) => callJsonHelper(pyodide, 'pdf_fonts_json', [
      jsonArgument(payload.runs, []), jsonArgument(payload.known, [])
    ])
  },
  [DIAGRAM_HELPER_OPERATIONS.VALIDATE_CONFIG_OVERRIDES]: {
    keys: [
      'mode',
      'config',
      'configOverrides',
      'managedPaths',
      'source'
    ],
    fileRoles: [],
    run: (pyodide, payload) => callJsonHelper(
      pyodide,
      'validate_web_config_overrides_json',
      [
        String(payload.mode || ''),
        jsonArgument(payload.config, null),
        jsonArgument(payload.configOverrides, {}),
        jsonArgument(payload.managedPaths, []),
        String(payload.source || '')
      ]
    )
  },
  [DIAGRAM_HELPER_OPERATIONS.EXTRACT_FIRST_FASTA]: {
    keys: ['files', 'format', 'regionSpec', 'recordSelector', 'reverseFlag'],
    fileRoles: ['source'],
    run: (pyodide, payload, paths, operation) => callJsonHelper(
      pyodide,
      'extract_first_fasta',
      [
        requireHelperFile(paths, 'source', operation),
        String(payload.format || 'genbank').trim().toLowerCase(),
        payload.regionSpec ?? null,
        payload.recordSelector ?? null,
        payload.reverseFlag ? '1' : '0'
      ]
    )
  },
  [DIAGRAM_HELPER_OPERATIONS.EXTRACT_CDS_PROTEIN_FASTA]: {
    keys: [
      'files',
      'format',
      'regionSpec',
      'recordSelector',
      'reverseFlag',
      'recordIndex',
      'recordInstanceKey',
      'featureOverrides'
    ],
    fileRoles: ['source', 'fasta', 'visibility'],
    run: (pyodide, payload, paths, operation) => {
      const format = String(payload.format || 'genbank').trim().toLowerCase();
      const fastaPath = paths.get('fasta') || null;
      if (format === 'gff' && !fastaPath) {
        throw new TypeError(`Diagram helper '${operation}' requires the 'fasta' file for GFF3.`);
      }
      return callJsonHelper(pyodide, 'extract_cds_protein_fasta', [
        requireHelperFile(paths, 'source', operation),
        format,
        fastaPath,
        payload.regionSpec ?? null,
        payload.recordSelector ?? null,
        payload.reverseFlag ? '1' : '0',
        payload.recordIndex ?? null,
        payload.recordInstanceKey ?? null,
        paths.get('visibility') || null,
        // This record's canonical featureOverrides rows (design Q4, 3.3).
        payload.featureOverrides == null ? null : JSON.stringify(payload.featureOverrides)
      ]);
    }
  },
  [DIAGRAM_HELPER_OPERATIONS.BUILD_PROTEIN_LOSAT_CACHE_KEYS]: {
    keys: ['identityManifest', 'pairs'],
    fileRoles: [],
    run: (pyodide, payload) => callJsonHelper(
      pyodide,
      'build_protein_losat_cache_keys_json',
      [
        jsonArgument(payload.identityManifest, {}),
        jsonArgument(payload.pairs, [])
      ]
    )
  },
  [DIAGRAM_HELPER_OPERATIONS.PROMOTE_LEGACY_LOSATP_CACHE]: {
    keys: [
      'candidates',
      'queryFasta',
      'subjectFasta',
      'queryProteinMap',
      'subjectProteinMap',
      'identityManifest',
      'expectedOptions'
    ],
    fileRoles: [],
    run: (pyodide, payload) => callJsonHelper(
      pyodide,
      'promote_legacy_losatp_cache_candidates',
      [
        jsonArgument(payload.candidates, []),
        String(payload.queryFasta || ''),
        String(payload.subjectFasta || ''),
        jsonArgument(payload.queryProteinMap, {}),
        jsonArgument(payload.subjectProteinMap, {}),
        jsonArgument(payload.identityManifest, {}),
        jsonArgument(payload.expectedOptions, {})
      ]
    )
  },
  [DIAGRAM_HELPER_OPERATIONS.RESOLVE_LEGACY_PROTEIN_REFERENCES]: {
    keys: ['proteinRecords', 'identityManifest', 'referenceIds'],
    fileRoles: [],
    run: (pyodide, payload) => callJsonHelper(
      pyodide,
      'resolve_legacy_protein_reference_map_json',
      [
        jsonArgument(payload.proteinRecords, []),
        jsonArgument(payload.identityManifest, {}),
        jsonArgument(payload.referenceIds, [])
      ]
    )
  },
  [DIAGRAM_HELPER_OPERATIONS.CONVERT_LOSATP_PAIRS_TO_GENOMIC_PAYLOAD]: {
    keys: [
      'files',
      'mode',
      'maxHits',
      'bitscore',
      'evalue',
      'identity',
      'alignmentLength',
      'collinearMinAnchors',
      'collinearMaxUnitGap',
      'collinearUnitMode',
      'collinearColorMode',
      'collinearAnchorMode',
      'collinearMaxDiagonalDrift',
      'collinearMaxConflictsInMergeGap',
      'collinearMaxParalogLinksPerOrthogroup',
      'collinearSearchScope',
      'collinearInferOrthogroups',
      'orthogroupMembershipMode',
      'orthogroupMemberMaxHits',
      'collinearMergeOrientation',
      'explicitDisplayPairs'
    ],
    fileRoles: ['pairs', 'rawTsv'],
    run: (pyodide, payload, paths, operation, workspace) => {
      const canonicalPath = `${workspace}/canonical-comparison.json`;
      const result = callJsonHelper(
      pyodide,
      'convert_losatp_blastp_pairs_to_genomic_payload',
      [
        requireHelperFile(paths, 'pairs', operation),
        requireHelperFile(paths, 'rawTsv', operation),
        payload.mode ?? 'pairwise',
        payload.maxHits ?? 5,
        payload.bitscore ?? 50,
        payload.evalue ?? '1e-2',
        payload.identity ?? 0,
        payload.alignmentLength ?? 0,
        payload.collinearMinAnchors ?? 1,
        payload.collinearMaxUnitGap ?? 0,
        payload.collinearUnitMode ?? 'auto',
        payload.collinearColorMode ?? 'orientation',
        payload.collinearAnchorMode ?? 'rbh',
        payload.collinearMaxDiagonalDrift ?? 0,
        payload.collinearMaxConflictsInMergeGap ?? 1,
        payload.collinearMaxParalogLinksPerOrthogroup ?? 2,
        payload.collinearSearchScope ?? 'adjacent',
        payload.orthogroupMembershipMode ?? ORTHOGROUP_INFERENCE_VERSION,
        payload.orthogroupMemberMaxHits ?? null,
        payload.collinearMergeOrientation ?? 'either',
        payload.collinearInferOrthogroups ?? true,
        canonicalPath,
        payload.explicitDisplayPairs === true
      ]
      );
      if (result?.canonicalResource) {
        const bytes = pyodide.FS.readFile(canonicalPath);
        if (bytes.byteLength !== result.canonicalResource.size) {
          throw new Error("Analysis canonical resource byte size does not match.");
        }
        result.canonicalResource.bytes = bytes;
      }
      return result;
    }
  },
  [DIAGRAM_HELPER_OPERATIONS.CONVERT_MAIN_SESSION_COMPARISON_FRAME]: {
    keys: ['tableText', 'queryFrame', 'subjectFrame'],
    fileRoles: [],
    run: (pyodide, payload) => callJsonHelper(
      pyodide,
      'main_session_table_text_to_search_frame_json',
      [
        String(payload.tableText || ''),
        jsonArgument(payload.queryFrame, {}),
        jsonArgument(payload.subjectFrame, {})
      ]
    )
  },
  [DIAGRAM_HELPER_OPERATIONS.HYDRATE_PROTEIN_LOSAT_TSV]: {
    keys: ['entry', 'identityManifest'],
    fileRoles: [],
    run: (pyodide, payload) => callJsonHelper(
      pyodide,
      'hydrate_protein_losat_tsv_json',
      [jsonArgument(payload.entry, {}), jsonArgument(payload.identityManifest, {})]
    )
  },
  [DIAGRAM_HELPER_OPERATIONS.LIST_SEQUENCE_RECORDS]: {
    keys: ['files', 'format'],
    fileRoles: ['source'],
    run: (pyodide, payload, paths, operation) => callJsonHelper(
      pyodide,
      'list_sequence_records',
      [
        requireHelperFile(paths, 'source', operation),
        String(payload.format || 'genbank').trim().toLowerCase()
      ]
    )
  },
  [DIAGRAM_HELPER_OPERATIONS.READ_COMPARISON_SEQUENCE]: {
    keys: ['files'],
    fileRoles: ['source'],
    run: (pyodide, _payload, paths, operation) => callJsonHelper(
      pyodide,
      'read_comparison_sequence_json',
      [requireHelperFile(paths, 'source', operation)]
    )
  },
  [DIAGRAM_HELPER_OPERATIONS.LIST_GFF_FASTA_RECORDS]: {
    keys: ['files'],
    fileRoles: ['gff', 'fasta'],
    run: (pyodide, _payload, paths, operation) => callJsonHelper(
      pyodide,
      'list_gff_fasta_records',
      [
        requireHelperFile(paths, 'gff', operation),
        requireHelperFile(paths, 'fasta', operation)
      ]
    )
  },
  [DIAGRAM_HELPER_OPERATIONS.MEASURE_LEGEND_TEXT]: {
    keys: ['caption', 'fontFamily', 'fontSize', 'config', 'configOverrides'],
    fileRoles: [],
    run: (pyodide, payload) => callJsonHelper(
      pyodide,
      'measure_legend_text_json',
      [
        String(payload.caption || ''),
        String(payload.fontFamily || 'Arial'),
        payload.fontSize ?? 14,
        jsonArgument(payload.config, null),
        jsonArgument(payload.configOverrides, {})
      ]
    )
  }
});

const runHelperOperation = async ({ operation, payload, requestId } = {}) => {
  if (!runtime?.pyodide) {
    throw serializeError({ code: 'WORKER_INIT', operation, stage: 'initialization' });
  }
  const normalizedOperation = String(operation || '').trim();
  const spec = Object.prototype.hasOwnProperty.call(
    HELPER_OPERATION_SPECS,
    normalizedOperation
  )
    ? HELPER_OPERATION_SPECS[normalizedOperation]
    : null;
  if (!spec) {
    throw serializeError({ code: 'HELPER_PROTOCOL', operation: normalizedOperation, stage: 'request-validation' });
  }
  let failureStage = 'request-validation';
  try {
    const normalizedPayload = requirePayloadObject(payload, normalizedOperation);
    assertAllowedPayloadKeys(normalizedPayload, normalizedOperation, [
      ...spec.keys
    ]);
    failureStage = 'resource-staging';
    return await withRequestWorkspace(
      runtime.pyodide,
      requestId,
      'helper',
      async (workspace) => {
        const paths = stageHelperFiles(
          runtime.pyodide,
          workspace,
          normalizedPayload.files,
          spec.fileRoles
        );
        failureStage = 'helper';
        return spec.run(runtime.pyodide, normalizedPayload, paths, normalizedOperation, workspace);
      }
    );
  } catch (error) {
    throw serializeError(error, {
      operation: normalizedOperation, stage: failureStage,
      code: failureStage === 'request-validation' ? 'HELPER_PROTOCOL' : 'UNKNOWN'
    });
  }
};

const runGeneration = async ({
  request,
  resourceManifest,
  stagedResources,
  fontIdentities = [],
  heldFeatureCatalogRevision = null,
  heldSegments: heldSegmentsInput = null,
  requestId,
  testLifecycleEnabled = false
} = {}) => {
  if (!runtime?.pyodide) {
    throw serializeError({ code: 'WORKER_INIT', operation: 'generate', stage: 'initialization' });
  }
  if (!request || typeof request !== 'object' || Array.isArray(request)) {
    throw serializeError({ code: 'HELPER_PROTOCOL', operation: 'generate', stage: 'request-validation' });
  }
  // The catalog revision the page holds for the shown drawing, with the
  // digest of its bytes (gbdraw/web/CLAUDE.md, Runtime data flow).
  const held = readFeatureCatalogRevision(
    heldFeatureCatalogRevision,
    () => serializeError({ code: 'HELPER_PROTOCOL', operation: 'generate', stage: 'request-validation' })
  );
  // The children texts the page holds (content revisions); Python omits their text.
  const heldSegments = readHeldSegments(
    heldSegmentsInput,
    () => serializeError({ code: 'HELPER_PROTOCOL', operation: 'generate', stage: 'request-validation' })
  );
  const { pyodide } = runtime;
  try {
    fontRegistryFor(pyodide).requireRegistered(fontIdentities);
  } catch (error) {
    throw serializeError(error, { operation: 'generate', stage: 'request-validation' });
  }
  const runWrapper = pyodide.globals.get('run_canonical_request_wrapper');
  const workspace = `/gbdraw-web-render-${Number(requestId) || 0}`;
  let result = null;
  let resultHandle = null;
  let primaryError = null;
  let pythonWrapperStarted = false;
  let failureStage = 'resource-staging';
  try {
    const newlyStagedResourceBytes = (
      Array.isArray(stagedResources) ? stagedResources : []
    ).reduce(
      (total, entry) => total + Number(entry?.bytes?.byteLength || 0),
      0
    );
    emitTestLifecycle(
      testLifecycleEnabled,
      requestId,
      'worker-workspace-preparation-start'
    );
    emitTestLifecycle(
      testLifecycleEnabled,
      requestId,
      'worker-resource-linking-start'
    );
    self.postMessage({ type: 'progress', requestId, stage: 'preparing-resources' });
    const resourcePaths = await prepareCanonicalResources(
      pyodide, workspace, { resourceManifest, stagedResources },
      { testLifecycleEnabled, requestId }
    );
    emitTestLifecycle(
      testLifecycleEnabled,
      requestId,
      'worker-resource-linking-end',
      {
        referencedResourceCount: Object.keys(resourcePaths).length,
        newlyStagedResourceBytes
      }
    );
    emitTestLifecycle(testLifecycleEnabled, requestId, 'worker-request-json-start');
    const requestJson = JSON.stringify(request);
    emitTestLifecycle(testLifecycleEnabled, requestId, 'worker-request-json-end', {
      characters: requestJson.length
    });
    emitTestLifecycle(testLifecycleEnabled, requestId, 'worker-resource-manifest-json-start');
    const resourcePathsJson = JSON.stringify(resourcePaths);
    const preparedResourceIdentitiesJson = JSON.stringify(
      buildPreparedResourceIdentityMap(resourceManifest)
    );
    emitTestLifecycle(testLifecycleEnabled, requestId, 'worker-resource-manifest-json-end', {
      characters: resourcePathsJson.length,
      preparedIdentityCharacters: preparedResourceIdentitiesJson.length
    });
    emitTestLifecycle(
      testLifecycleEnabled,
      requestId,
      'worker-workspace-preparation-end'
    );
    self.postMessage({ type: 'progress', requestId, stage: 'rendering' });
    emitTestLifecycle(testLifecycleEnabled, requestId, 'python-wrapper-start');
    pythonWrapperStarted = true;
    failureStage = 'render';
    resultHandle = runWrapper(
      requestJson,
      resourcePathsJson,
      workspace,
      testLifecycleEnabled,
      preparedResourceIdentitiesJson,
      held?.revision ?? null,
      JSON.stringify(heldSegments.map(({ revision }) => revision))
    );
    emitTestLifecycle(testLifecycleEnabled, requestId, 'python-wrapper-end');
    self.postMessage({ type: 'progress', requestId, stage: 'finalizing' });
    emitTestLifecycle(testLifecycleEnabled, requestId, 'result-object-conversion-start');
    failureStage = 'result-admission';
    result = typeof resultHandle?.toJs === 'function'
      ? resultHandle.toJs({ dict_converter: Object.fromEntries })
      : resultHandle;
    emitTestLifecycle(testLifecycleEnabled, requestId, 'result-object-conversion-end');
    const pythonDiagnostics = (
      result?._diagnostics
      && typeof result._diagnostics === 'object'
      && !Array.isArray(result._diagnostics)
    )
      ? result._diagnostics
      : null;
    if (result && typeof result === 'object') delete result._diagnostics;
    if (pythonDiagnostics) {
      emitTestLifecycle(testLifecycleEnabled, requestId, 'python-diagnostics', {
        timingsMs: pythonDiagnostics.timingsMs || {},
        metrics: pythonDiagnostics.metrics || {}
      });
    }
    emitTestLifecycle(testLifecycleEnabled, requestId, 'result-transport-ready', {
      transport: 'transferable-binary',
      featureCatalog: result?.featureCatalog != null ? 'sent' : result?.error ? 'none' : 'unchanged',
      heldFeatureCatalogRevision: Boolean(held),
      bytes: collectGenerationResultTransferList(result).reduce(
        (total, buffer) => total + buffer.byteLength,
        0
      )
    });
  } catch (error) {
    primaryError = serializeError(error, { operation: 'generate', stage: failureStage });
  }
  emitTestLifecycle(testLifecycleEnabled, requestId, 'worker-cleanup-start');
  let destroyError = null;
  try {
    resultHandle?.destroy?.();
  } catch (error) {
    destroyError = error;
  }
  try {
    runWrapper.destroy?.();
  } catch (error) {
    if (destroyError) {
      attachCleanupDiagnostic(
        destroyError,
        error
      );
    } else {
      destroyError = error;
    }
  }
  let workspaceError = null;
  try {
    let remainingEntries = [];
    if (pyodide.FS.analyzePath(workspace).exists) {
      remainingEntries = pyodide.FS.readdir(workspace).filter(
        (entry) => entry !== '.' && entry !== '..'
      );
      if (!pythonWrapperStarted) {
        removeWorkspace(pyodide, workspace);
      } else if (remainingEntries.length === 0) {
        // Pyodide's MEMFS can retain the now-empty top-level directory after
        // Python shutil.rmtree has removed its contents.
        pyodide.FS.rmdir(workspace);
      }
    }
    if (pyodide.FS.analyzePath(workspace).exists) {
      workspaceError = new Error(
        `Diagram render workspace cleanup invariant failed: ${workspace}` +
        (remainingEntries.length > 0
          ? ` (${remainingEntries.join(', ')})`
          : '')
      );
      try {
        removeWorkspace(pyodide, workspace);
      } catch (_cleanupError) {
        // Preserve the invariant error assembled above.
      }
    }
  } catch (error) {
    workspaceError = error;
  }
  emitTestLifecycle(testLifecycleEnabled, requestId, 'worker-cleanup-end');
  const outcome = resolveGenerationCleanupOutcome({
    result,
    primaryError,
    destroyError,
    workspaceError
  });
  if (outcome?.error) outcome.error = serializeError(outcome.error, { operation: 'generate' });
  if (!outcome?.error) {
    emitTestLifecycle(testLifecycleEnabled, requestId, 'worker-artifact-identity-start');
    await resolveGeneratedFeatureCatalog(outcome, held, heldSegments);
    emitTestLifecycle(testLifecycleEnabled, requestId, 'worker-artifact-identity-end', {
      fingerprint: outcome.artifactIdentity.fingerprint,
      retainedBytes: outcome.artifactIdentity.retainedBytes
    });
  }
  return outcome;
};

const runFeatureExtraction = async ({
  path,
  format = 'genbank',
  fastaPath = null,
  files = [],
  regionSpec = null,
  recordSelector = null,
  reverseFlag = false,
  selectedFeatures = null,
  featureVisibilityTablePath = null,
  includeBiologicalFeatures = false,
  requestId = 0
} = {}) => {
  if (!runtime?.pyodide) {
    throw serializeError({ code: 'WORKER_INIT', operation: 'feature-extraction', stage: 'initialization' });
  }
  const normalizedPath = String(path || '').trim();
  if (!normalizedPath) {
    throw serializeError({ code: 'HELPER_PROTOCOL', operation: 'feature-extraction', stage: 'request-validation', context: { field: 'input', reason: 'REQUIRED' } });
  }
  const normalizedFormat = String(format || 'genbank').trim().toLowerCase();
  const isGff = normalizedFormat === 'gff';
  const normalizedFastaPath = String(fastaPath || '').trim();
  if (isGff && !normalizedFastaPath) {
    throw serializeError({ code: 'FASTA_REQUIRED', operation: 'feature-extraction', stage: 'request-validation' });
  }

  const { pyodide } = runtime;
  const normalizedFiles = Array.isArray(files) ? files : [];
  const sourceFile = normalizedFiles.find(
    (file) => String(file?.path || '').trim() === normalizedPath
  );
  const fastaFile = isGff
    ? normalizedFiles.find((file) => String(file?.path || '').trim() === normalizedFastaPath)
    : null;
  const normalizedVisibilityPath = String(featureVisibilityTablePath || '').trim();
  const visibilityFile = normalizedVisibilityPath
    ? normalizedFiles.find((file) => String(file?.path || '').trim() === normalizedVisibilityPath)
    : null;
  const workspaceFiles = [
    { role: 'source', bytes: sourceFile?.bytes },
    ...(fastaFile ? [{ role: 'fasta', bytes: fastaFile.bytes }] : []),
    ...(visibilityFile ? [{ role: 'visibility', bytes: visibilityFile.bytes }] : [])
  ];

  let failureStage = 'resource-staging';
  try {
    return await withRequestWorkspace(pyodide, requestId, 'feature', async (workspace) => {
      const paths = stageHelperFiles(
        pyodide,
        workspace,
        workspaceFiles,
        ['source', 'fasta', 'visibility']
      );
      const sourcePath = requireHelperFile(paths, 'source', 'feature-extraction');
      const stagedFastaPath = isGff
        ? requireHelperFile(paths, 'fasta', 'feature-extraction')
        : null;
      const selectedFeaturesJson = Array.isArray(selectedFeatures) && selectedFeatures.length
        ? JSON.stringify(selectedFeatures)
        : null;
      failureStage = 'helper';
      return callJsonHelper(
        pyodide,
        isGff ? 'extract_features_from_gff_fasta' : 'extract_features_from_genbank',
        isGff
          ? [
            sourcePath,
            stagedFastaPath,
            regionSpec || null,
            recordSelector || null,
            reverseFlag ? '1' : '0',
            selectedFeaturesJson,
            paths.get('visibility') || null,
            Boolean(includeBiologicalFeatures)
          ]
          : [
            sourcePath,
            regionSpec || null,
            recordSelector || null,
            reverseFlag ? '1' : '0',
            selectedFeaturesJson,
            paths.get('visibility') || null,
            Boolean(includeBiologicalFeatures)
          ]
      );
    });
  } catch (error) {
    throw serializeError(error, { operation: 'feature-extraction', stage: failureStage });
  }
};

const handleWorkerMessage = async (data) => {
  const { id, requestId, type } = data;
  try {
    if (type === 'init') {
      const initialized = await initializeRuntime(data);
      // The replay of the page's fonts into a new Worker, font by font; the
      // reply names each font it lost (OV-686), and a Generate that names one
      // fails (FONT_NOT_LOADED).
      const lostFonts = await fontRegistryFor(initialized.pyodide).replay(data.fonts);
      self.postMessage({
        id,
        type: 'init',
        ok: true,
        capabilities: initialized.capabilities,
        lostFonts
      });
      return;
    }
    if (type === 'ping') {
      self.postMessage({ id, type: 'ping', ok: true });
      return;
    }
    if (type === 'feature-extraction') {
      const result = await runFeatureExtraction({
        ...(data.payload || {}),
        requestId
      });
      await sendAuxiliaryResult('feature-extraction', requestId, result);
      return;
    }
    if (type === 'helper') {
      const result = await runHelperOperation({
        operation: data.operation,
        payload: data.payload || {},
        requestId
      });
      await sendAuxiliaryResult('helper', requestId, result);
      return;
    }
    if (type !== 'run') {
      throw new Error(`Unsupported diagram generation worker message type '${type || '(blank)'}'.`);
    }

    emitTestLifecycle(
      data.testLifecycleEnabled,
      requestId,
      'run-message-received'
    );
    const results = await runGeneration({
      ...(data.payload || {}),
      requestId,
      testLifecycleEnabled: data.testLifecycleEnabled
    });
    self.postMessage(
      { requestId, type: 'run', ok: true, results },
      collectGenerationResultTransferList(results)
    );
  } catch (error) {
    self.postMessage({
      id,
      requestId,
      type: type || 'run',
      ok: false,
      error: serializeError(error, {
        operation: type === 'run' ? 'generate' : type === 'feature-extraction' ? 'feature-extraction' : data.operation,
        stage: type === 'init' ? 'initialization' : 'unknown',
        code: type === 'init' ? 'WORKER_INIT' : 'UNKNOWN'
      })
    });
  }
};

self.onmessage = (event) => {
  const data = event.data || {};
  if (data.type === 'auxiliary-ack') {
    if (auxiliaryAcknowledgement?.requestId === data.requestId) {
      const pending = auxiliaryAcknowledgement;
      auxiliaryAcknowledgement = null;
      pending.resolve();
    }
    return;
  }
  const scheduled = operationQueue.then(
    () => handleWorkerMessage(data),
    () => handleWorkerMessage(data)
  );
  operationQueue = scheduled.catch(() => {});
};
