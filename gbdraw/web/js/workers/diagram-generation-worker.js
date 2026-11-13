import { sendBoundedJson } from '../services/bounded-json-transport.js';
import { sha256Hex } from '../services/byte-utils.js';
import { PYTHON_HELPERS } from '../app/python-helpers.js';
import {
  DIAGRAM_HELPER_OPERATIONS,
  DIAGRAM_HELPER_STAGING_OPERATIONS,
  STAGED_RESOURCE_PAYLOAD_KEYS,
  readFeatureCatalogRevision,
  readHeldSegments,
  readSvgDocumentParts
} from '../services/diagram-worker-protocol.js';
import { diagnosticError, normalizeUserFacingError, OUT_OF_MEMORY_CODE } from '../utils/error-normalization.js';
import { msg } from '../utils/i18n.js';
import { ORTHOGROUP_INFERENCE_VERSION } from '../contracts/schema.generated.js';
import { recordSourceFilesOf } from '../services/input-types.js';
import { createWorkerFontRegistry } from './worker-font-registry.js';
import {
  fetchPyodideSnapshot, installRuntimePackages, pyodideSnapshotFailureDetail, pyodideSnapshotStep, runtimeWheelDownloads
} from './pyodide-snapshot.js';

let runtimePromise = null;
let runtime = null;
let operationQueue = Promise.resolve();
let auxiliaryAcknowledgement = null;
// The page's interrupt buffer (DiagramWorkerInitMessage); null where a Cancel terminates this Worker.
let interruptBuffer = null;
// The last request ID a Cancel named (DiagramWorkerCancelMessage): no reply of
// a request up to it is posted.
let canceledThrough = 0;
const isCanceled = (requestId) => typeof requestId === 'number' && requestId <= canceledThrough;
const sendAuxiliaryResult = async (type, requestId, result) => {
  if (isCanceled(requestId)) return;
  let wholeReply = false;
  try {
    await sendBoundedJson(result, (part, transfers = []) => {
      // A Cancel during a bounded reply: its remaining parts are dropped.
      if (isCanceled(requestId)) return Promise.resolve();
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
    if (!wholeReply && !isCanceled(requestId)) self.postMessage({ type, requestId, ok: true });
  } finally {
    auxiliaryAcknowledgement = null;
  }
};
const RENDER_RESOURCE_CACHE = '/gbdraw-web-render-resource-cache';
const RENDER_WORKSPACE_MARKER = '.gbdraw-worker-render-workspace';
const RENDER_RESOURCE_ID_RE = /^[a-z][a-z0-9]*(?:-[a-z0-9]+)*$/;
// A resource's token names its bytes (`contentTokenOfBytes` in services/file-content-cache.js).
const RENDER_RESOURCE_TOKEN_RE = /^[0-9a-f]{64}$/;
// The resources the Worker holds, by token: the current request's and those
// the page's drawings have committed (`retainedResourceTokens`).
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

/**
 * What Python's prepared input cache keeps besides this request's inputs
 * (A4b): the inputs of each drawing's last request, while the resources they
 * name are requested or retained, and, for this drawing, the entries of the
 * files its draft holds while the request omits them (`heldResourceTokens`,
 * which the Worker retains). No drawing: only this request's.
 */
export const buildPreparedInputRetention = (drawingId, retainedResourceTokens, heldResourceTokens) => {
  if (drawingId === null || drawingId === undefined) return null;
  // A drawing ID has the form of a resource ID.
  if (typeof drawingId !== 'string' || !RENDER_RESOURCE_ID_RE.test(drawingId)) {
    throw diagnosticError('HELPER_PROTOCOL', { reason: 'retention-drawing' });
  }
  // A render message always names the held tokens (empty when none); only helpers omit them.
  if (!Array.isArray(heldResourceTokens)) {
    throw diagnosticError('HELPER_PROTOCOL', { reason: 'retention-held' });
  }
  const resourceTokens = [...(retainedResourceTokens ?? [])];
  const heldTokens = [...heldResourceTokens];
  if (!heldTokens.every((token) => resourceTokens.includes(token))) {
    throw diagnosticError('HELPER_PROTOCOL', { reason: 'retention-held' });
  }
  return { drawing: drawingId, resourceTokens, heldTokens };
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

/**
 * The Pyodide file system calls a reply buffer uses.
 * @typedef {{ readFile(path: string): Uint8Array<ArrayBuffer>, unlink(path: string): void }} ReplyBufferFiles
 */

/**
 * The file each reply chunk crosses through. Python writes a chunk there and
 * the buffer takes it out at once, so one path serves every buffer of a Worker.
 */
const REPLY_CHUNK_PATH = '/gbdraw-web-reply-chunk';

/**
 * A reply buffer Python writes (gbdraw/web_support/reply_buffers.py): the
 * metadata or feature catalog JSON, or a Result's children texts, in UTF-8
 * chunks of about 1 MiB. Python writes each chunk to the file at `path`, and
 * `take` moves it out of the file system. A chunk crosses through a file
 * because Pyodide's buffer views read a heap address as a signed 32-bit number
 * and fail once the chunk sits above 2 GiB of the heap, while a file write
 * reads it as unsigned (OV-2800). Python returns the buffer itself; the Worker
 * joins its chunks after Python returns (`finishReplyBuffers`), so the reply
 * briefly takes about twice its size here: the chunk copies and the joined
 * buffer.
 */
export class ReplyBuffer {
  /** @type {ReplyBufferFiles} */
  #files;
  /** @type {Uint8Array<ArrayBuffer>[]} */
  #chunks = [];
  #byteLength = 0;

  /** @param {ReplyBufferFiles} files */
  constructor(files) {
    this.#files = files;
  }

  /** The file Python writes each chunk to before it calls `take`. */
  get path() {
    return REPLY_CHUNK_PATH;
  }

  /** Move the chunk Python wrote out of the file system into this buffer. */
  take() {
    const copy = this.#files.readFile(REPLY_CHUNK_PATH);
    this.#files.unlink(REPLY_CHUNK_PATH);
    this.#chunks.push(copy);
    this.#byteLength += copy.byteLength;
  }

  /** The chunks, in order, as one exact-length buffer; the chunk copies are dropped. */
  finish() {
    const bytes = new Uint8Array(this.#byteLength);
    let offset = 0;
    for (const copy of this.#chunks) {
      bytes.set(copy, offset);
      offset += copy.byteLength;
    }
    this.#chunks = [];
    return bytes.buffer;
  }
}

/**
 * A reply buffer over the Worker's file system; the Worker hands Python
 * `() => createReplyBuffer(pyodide.FS)` as the factory of the reply's buffers.
 * @param {ReplyBufferFiles} files
 */
export const createReplyBuffer = (files) => new ReplyBuffer(files);
/** @param {{ FS: ReplyBufferFiles }} pyodide The factory of reply buffers Python is handed. */
const replyBuffersOf = (pyodide) => () => createReplyBuffer(pyodide.FS);
const UTF8_DECODER = new TextDecoder();

/**
 * Join the reply buffers of a converted Generate result in place: its
 * `metadata`, its `featureCatalog` and each Result's `segments`. Called at
 * result admission, where a buffer too large to allocate fails as the
 * conversion of Python's bytes did before.
 * @param {unknown} result
 */
export const finishReplyBuffers = (result) => {
  if (!result || typeof result !== 'object') return;
  const { results } = /** @type {{ results?: unknown }} */ (result);
  for (const holder of [result, ...(Array.isArray(results) ? results : [])]) {
    if (!holder || typeof holder !== 'object') continue;
    const fields = /** @type {Record<string, unknown>} */ (holder);
    for (const key of ['metadata', 'featureCatalog', 'segments']) {
      const value = fields[key];
      if (value instanceof ReplyBuffer) fields[key] = value.finish();
    }
  }
};

/**
 * Python's Generate result as JavaScript values, with its reply buffers joined.
 * @param {{ toJs?: (options: { dict_converter: typeof Object.fromEntries }) => unknown } | null | undefined} resultHandle
 */
export const convertGenerationResult = (resultHandle) => {
  const result = typeof resultHandle?.toJs === 'function'
    ? resultHandle.toJs({ dict_converter: Object.fromEntries })
    : resultHandle;
  finishReplyBuffers(result);
  return result;
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
    } else if (await sha256Hex(segments.subarray(piece.offset, piece.offset + piece.bytes)) !== revision) {
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
      nameSha256: await sha256Hex(nameBytes),
      documentBytes: document.bytes,
      documentSha256: await sha256Hex(encoder.encode(JSON.stringify(document.manifest)))
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
      sha256: await sha256Hex(metadataBytes)
    },
    featureCatalog: featureCatalogDigest
      ? { bytes: featureCatalogDigest.bytes, sha256: featureCatalogDigest.sha256 }
      : null
  };
  const catalogBytes = featureCatalogDigest?.bytes ?? 0;
  const fingerprint = await sha256Hex(encoder.encode(JSON.stringify(manifest)));
  return Object.freeze({
    schema: 1,
    algorithm: 'SHA-256',
    fingerprint,
    resultBytes,
    resultNameBytes,
    metadataBytes: metadataBytes.byteLength + catalogBytes,
    // The page keeps the metadata and the feature catalog parsed, not their
    // UTF-8 bytes: one byte per source byte. Measured in Chromium (OV-3010),
    // the parsed E. coli O157 Sakai ring catalog (36.6 MB of JSON) holds 0.92
    // bytes per source byte. The page keeps each Result as its sanitized text
    // and its segments, so the main-thread handle adds `resultBytes` once more.
    // `resultBytes` counts the whole admitted document, held or sent. The
    // feature catalog counts as metadata whether it was sent or kept (D-38 OQ-2).
    retainedBytes: resultBytes + resultNameBytes + metadataBytes.byteLength + catalogBytes
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
    ? { bytes: catalogBytes.byteLength, sha256: await sha256Hex(catalogBytes) }
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

// The booted runtime with its capabilities. A runtime restored from a Pyodide
// snapshot first starts this Worker's own state, since the snapshot holds the
// state of the process that made it (OV-1000), and calls numpy and pandas, so
// a side module the restore linked wrongly fails here and the page boots cold.
const SNAPSHOT_BOOT_CHECK = `
_start_web_worker_state()
import numpy as _numpy, pandas as _pandas
if int(_numpy.arange(3).sum()) != 3 or _pandas.DataFrame({'a': [1]}).shape != (1, 1):
    raise RuntimeError('numpy or pandas returned a wrong result after the snapshot restore')
del _numpy, _pandas
`;
const bootedRuntime = (pyodide, boot) => {
  const raw = pyodide.runPython(`
${boot === 'snapshot' ? SNAPSHOT_BOOT_CHECK : ''}
import json as _json
from gbdraw.api import get_web_runtime_capabilities as _get_web_runtime_capabilities
_json.dumps(_get_web_runtime_capabilities(), sort_keys=True)
  `);
  return { pyodide, capabilities: JSON.parse(String(raw)), boot };
};

// The first-use preparation steps, reported as they start: a real step
// count, not a time estimate. The page translates the names where it shows them.
const RUNTIME_PREPARATION_STEPS = Object.freeze([
  msg('loading Python'), msg('loading the package installer'), msg('checking packages'), msg('installing packages'), msg('starting gbdraw')
]);

const initializeRuntime = async ({
  id = null,
  pyodideIndexUrl,
  pyodideModuleUrl,
  localWheelUrls = [],
  gbdrawWheelUrl,
  pyodideSnapshotUrl = null
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
    // The wheels download while Pyodide loads.
    const wheels = runtimeWheelDownloads({ localWheelUrls, gbdrawWheelUrl });
    // A test server may serve a snapshot of the runtime this function builds
    // (workers/pyodide-snapshot.js). An error in a step of its restore is a
    // snapshot failure (pyodideSnapshotStep); any other error is an ordinary
    // init error (OV-1021).
    const [{ loadPyodide }, snapshot] = await Promise.all([
      import(pyodideModuleUrl),
      pyodideSnapshotUrl
        ? pyodideSnapshotStep(() => fetchPyodideSnapshot({ folderUrl: pyodideSnapshotUrl, wheels }))
        : null
    ]);
    const load = () => loadPyodide({
      indexURL: pyodideIndexUrl,
      packageBaseUrl: pyodideIndexUrl,
      ...(snapshot && { _loadSnapshot: snapshot.memory, _snapshotDeserializer: snapshot.deserialize })
    });
    const pyodide = await (snapshot ? pyodideSnapshotStep(load) : load());
    if (snapshot) {
      reportStep(4);
      runtime = await pyodideSnapshotStep(async () => {
        await snapshot.complete(pyodide);
        return bootedRuntime(pyodide, 'snapshot');
      });
      return runtime;
    }
    reportStep(1);
    await installRuntimePackages(pyodide, {
      wheels,
      beforeInstall: async (installOrder) => {
        reportStep(2);
        await Promise.all(installOrder.map(({ bytes }) => bytes));
        reportStep(3);
      }
    });
    reportStep(4);
    await pyodide.runPythonAsync(PYTHON_HELPERS);
    runtime = bootedRuntime(pyodide, 'cold');
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
  rawTsv: 'raw-losatp.tsv',
  embedItem: 'embed-item.json'
});

const requirePayloadObject = (payload, operation) => {
  if (!payload || typeof payload !== 'object' || Array.isArray(payload)) {
    throw new TypeError(`Diagram helper '${operation}' requires an object payload.`);
  }
  return payload;
};

// The record's display topology an extraction reads (`display.isCircular` of its request record): true, false, or
// null for the topology its file states.
/** @param {unknown} value @param {string} operation @returns {boolean | null} */
const displayCircularArgument = (value, operation) => {
  if (value === null || typeof value === 'boolean') return value;
  throw new TypeError(`Diagram helper '${operation}' requires 'isCircular' to be true, false or null.`);
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

/**
 * The bytes of an input the page sent as a Blob (OV-2749): an input's copy, its
 * File, or bytes the page read. The Worker reads it whole, once, into memory it
 * then gives MEMFS: Python reads MEMFS line by line about 14 times faster than
 * a WORKERFS mount of the Blob, which reads the Blob anew on every call.
 * A read the browser refuses (the file was deleted, moved, or changed) and one
 * that finds no memory for the bytes both reject; a one-byte read tells them apart.
 * @param {Blob} blob
 * @returns {Promise<Uint8Array>}
 */
const readInputBlob = async (blob) => {
  try {
    return new Uint8Array(await blob.arrayBuffer());
  } catch {
    const readable = await blob.slice(0, 1).arrayBuffer().then(() => true, () => false);
    throw readable
      ? { code: 'RESOURCE_INVALID', stage: 'resource-staging', context: { reason: 'INPUT_MEMORY' } }
      : { code: 'INPUT_UNREADABLE', stage: 'resource-staging' };
  }
};

// A role named `<role>[]` takes several files in order (the files of one
// comparison genome); its paths are a list under `<role>`.
const stageHelperFiles = async (pyodide, workspace, files, allowedRoles) => {
  const allowed = new Set(allowedRoles);
  const paths = new Map();
  const entries = Array.isArray(files) ? files : [];
  for (const file of entries) {
    const role = String(file?.role || '').trim();
    const listed = allowed.has(`${role}[]`);
    if (!allowed.has(role) && !listed) {
      throw new TypeError(`Unexpected diagram helper file role '${role || '(blank)'}'.`);
    }
    if (!listed && paths.has(role)) {
      throw new TypeError(`Duplicate diagram helper file role '${role}'.`);
    }
    if (!(file?.blob instanceof Blob)) {
      throw new TypeError(`Diagram helper file '${role}' requires a Blob.`);
    }
    const placed = listed ? (paths.get(role) || []) : null;
    const path = placed
      ? `${workspace}/${placed.length + 1}-${HELPER_FILE_NAMES[role]}`
      : `${workspace}/${HELPER_FILE_NAMES[role]}`;
    paths.set(role, placed ? [...placed, path] : path);
    // MEMFS keeps the bytes it is given: the Worker read them for this file alone.
    pyodide.FS.writeFile(path, await readInputBlob(file.blob), { canOwn: true });
  }
  return paths;
};

// A helper reads each file where the Worker staged it, under its role's name
// (`HELPER_FILE_NAMES`). A refusal that names that file (a gzip input over its
// cap, PD-OI-116) names the file the user chose instead, the helper file's
// `name`, as Generate's refusal names its resource's (gbdraw/web_support/request_render.py).
export const nameStagedHelperInput = (error, files) => {
  const staged = error?.context?.inputName;
  if (typeof staged !== 'string' || !Object.values(HELPER_FILE_NAMES).includes(staged)) return error;
  const file = (Array.isArray(files) ? files : []).find((entry) => HELPER_FILE_NAMES[String(entry?.role || '')] === staged);
  const { inputName: _staged, ...context } = error.context;
  const name = typeof file?.name === 'string' ? file.name : '';
  return serializeError({ ...error, context: name ? { ...context, inputName: name } : context });
};

const ensureRenderResourceCache = (pyodide) => {
  if (!pyodide.FS.analyzePath(RENDER_RESOURCE_CACHE).exists) {
    pyodide.FS.mkdir(RENDER_RESOURCE_CACHE);
  }
};

/**
 * Link the manifest's resources into the request workspace, writing the staged
 * bytes into the cache first; then release every cached resource that neither
 * the manifest nor `retainedResourceTokens` names when `prune` is set (a render).
 */
const stageRenderResources = async (
  pyodide,
  workspace,
  resourceManifest,
  stagedResources,
  retainedResourceTokens,
  testLifecycleEnabled,
  requestId,
  prune
) => {
  if (!Array.isArray(resourceManifest) || !Array.isArray(stagedResources)) {
    throw new TypeError('Diagram generation requires a staged resource manifest.');
  }
  const retained = retainedResourceTokens === undefined ? [] : retainedResourceTokens;
  if (!Array.isArray(retained) || !retained.every((token) => RENDER_RESOURCE_TOKEN_RE.test(String(token)))) {
    throw new TypeError('Invalid retained render resource tokens.');
  }
  ensureRenderResourceCache(pyodide);
  const stagedByToken = new Map();
  stagedResources.forEach((entry) => {
    const cacheToken = String(entry?.cacheToken || '');
    if (!RENDER_RESOURCE_TOKEN_RE.test(cacheToken) || stagedByToken.has(cacheToken)) {
      throw new TypeError(`Invalid or duplicate cache token for staged render resource '${String(entry?.resourceId || '')}'.`);
    }
    if (!(entry?.blob instanceof Blob)) {
      throw new TypeError(`Staged render resource '${String(entry?.resourceId || '')}' requires a Blob.`);
    }
    stagedByToken.set(cacheToken, entry);
  });

  const manifestIds = new Set();
  const manifestTokens = new Set();
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
    manifestTokens.add(cacheToken);
    let cached = cachedRenderResources.get(cacheToken);
    const staged = stagedByToken.get(cacheToken);
    if (staged) {
      if (staged.blob.size !== size) {
        throw new TypeError(`Staged render resource '${resourceId}' does not match its manifest.`);
      }
      emitTestLifecycle(testLifecycleEnabled, requestId, 'worker-resource-stage-start', {
        resourceId,
        bytes: size
      });
      // The token names the bytes, so a file cached under it at this size already holds them. It
      // stays: records read from it keep reading it, and Python refuses a file replaced under them.
      if (!cached || cached.size !== size || !pyodide.FS.analyzePath(cached.path).exists) {
        const cachePath = `${RENDER_RESOURCE_CACHE}/${cacheToken}.bin`;
        if (typeof pyodide.FS.createDataFile !== 'function') {
          throw new Error('Pyodide render resource ownership transfer is unavailable.');
        }
        // The cache owns the bytes it read (canOwn): it holds them, not the Blob, so a copy
        // deleted later (PD-OI-121) leaves the resources a later Generate reuses intact.
        const bytes = await readInputBlob(staged.blob);
        if (pyodide.FS.analyzePath(cachePath).exists) pyodide.FS.unlink(cachePath);
        pyodide.FS.createDataFile(
          RENDER_RESOURCE_CACHE,
          `${cacheToken}.bin`,
          bytes,
          true,
          true,
          true
        );
        cached = { path: cachePath, size };
        cachedRenderResources.set(cacheToken, cached);
      }
      staged.blob = null;
      stagedByToken.delete(cacheToken);
      emitTestLifecycle(testLifecycleEnabled, requestId, 'worker-resource-stage-end', {
        resourceId,
        bytes: size
      });
    }
    if (
      !cached
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
  const [unexpected] = stagedByToken.values();
  if (unexpected) {
    throw new TypeError(`Unexpected staged render resource '${String(unexpected.resourceId || '')}'.`);
  }
  const kept = new Set([...manifestTokens, ...retained]);
  cachedRenderResources.forEach((cached, cacheToken) => {
    if (!prune || kept.has(cacheToken)) return;
    if (pyodide.FS.analyzePath(cached.path).exists) pyodide.FS.unlink(cached.path);
    cachedRenderResources.delete(cacheToken);
  });
  stagedResources.splice(0);
  return resourcePaths;
};

// Rendering and source-bound helpers enter the same canonical staging boundary.
// Only a render prunes the resources its manifest does not name: a helper reads
// the cache and leaves it to the next Generate (Worker cache rules, OV-1590).
const prepareCanonicalResources = (
  pyodide, workspace, { resourceManifest, stagedResources, retainedResourceTokens },
  { testLifecycleEnabled = false, requestId = 0, prune = false } = {}
) => stageRenderResources(
  pyodide, workspace, resourceManifest, stagedResources, retainedResourceTokens, testLifecycleEnabled, requestId, prune
);

// The staged paths of a canonical record source, in the order its kind names its
// files (recordSourceFilesOf, from the generated RECORD_SOURCE_FILES).
const stagedRecordSourcePaths = (recordSource, resourcePaths, operation) => {
  const files = recordSourceFilesOf(recordSource?.kind);
  const ids = files.map(({ field }) => recordSource[field]);
  const paths = ids.map((id) => (typeof id === 'string' && Object.hasOwn(resourcePaths, id) ? resourcePaths[id] : null));
  if (!paths.length || paths.some((path) => !path)) {
    throw new TypeError(`Diagram helper '${operation}' requires a staged record source.`);
  }
  return { source: paths[0], fasta: paths[1] ?? null };
};

const requireHelperFile = (paths, role, operation) => {
  const path = paths.get(role);
  if (!path) throw new TypeError(`Diagram helper '${operation}' requires the '${role}' file.`);
  return path;
};

// The files of a `<role>[]` role, at least one.
const requireHelperFiles = (paths, role, operation) => {
  const list = paths.get(role);
  if (!Array.isArray(list) || list.length === 0) {
    throw new TypeError(`Diagram helper '${operation}' requires the '${role}' files.`);
  }
  return list;
};

export const callJsonHelper = (pyodide, helperName, args) => {
  const helper = pyodide.globals.get('call_web_json_helper');
  let result = null;
  let primary = null;
  try {
    if (typeof helper !== 'function') throw { code: 'HELPER_PROTOCOL', stage: 'helper' };
    try {
      const reply = helper(helperName, ...args);
      // A helper that replies in UTF-8 (`_replied_in_utf8` in python-helpers.js) returns its reply buffer.
      result = JSON.parse(reply instanceof ReplyBuffer ? UTF8_DECODER.decode(reply.finish()) : String(reply || 'null'));
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
    keys: ['request', 'projection'],
    fileRoles: [],
    // The projection (the committed request and its orientations) is required: Python places each candidate as
    // the renderer does (D5-13).
    run: async (pyodide, payload, _paths, _operation, workspace) => callJsonHelper(
      pyodide, 'resolve_similarity_alignment_json', [
        jsonArgument(payload.request, null), jsonArgument(payload.projection, null),
        jsonArgument(await prepareCanonicalResources(pyodide, `${workspace}/projection`, payload), {}), workspace,
        jsonArgument(buildPreparedResourceIdentityMap(payload.resourceManifest), {})
      ]
    )
  },
  // Align at feature: Python resolves the anchor on the committed request,
  // staged as Generate stages it, and returns the record translations.
  [DIAGRAM_HELPER_OPERATIONS.RESOLVE_FEATURE_ANCHOR]: {
    keys: ['canonicalRequest', 'selector'],
    fileRoles: [],
    run: async (pyodide, payload, _paths, _operation, workspace) => callJsonHelper(
      pyodide, 'resolve_feature_anchor_json', [
        jsonArgument(payload.canonicalRequest, null),
        jsonArgument(payload.selector, null),
        jsonArgument(await prepareCanonicalResources(pyodide, `${workspace}/request`, payload), {}),
        `${workspace}/output`,
        jsonArgument(buildPreparedResourceIdentityMap(payload.resourceManifest), {})
      ]
    )
  },
  // Load Feature Edits TSV: Python reads the table against the records of the
  // committed request, staged as Generate stages them (design Q4 6.4, R4).
  [DIAGRAM_HELPER_OPERATIONS.READ_FEATURE_OVERRIDE_TABLE]: {
    keys: ['files', 'canonicalRequest'],
    fileRoles: ['featureOverrides'],
    run: async (pyodide, payload, paths, operation, workspace) => callJsonHelper(
      pyodide, 'read_feature_override_table_json', [
        requireHelperFile(paths, 'featureOverrides', operation),
        jsonArgument(payload.canonicalRequest, null),
        jsonArgument(await prepareCanonicalResources(pyodide, `${workspace}/request`, payload), {}),
        `${workspace}/output`,
        jsonArgument(buildPreparedResourceIdentityMap(payload.resourceManifest), {})
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
  // The qualifiers the candidate conditions of a color edit read (D5): the page sends the synthesis those only.
  [DIAGRAM_HELPER_OPERATIONS.SYNTHESIS_QUALIFIERS]: {
    keys: ['request'],
    fileRoles: [],
    run: (pyodide, payload) => callJsonHelper(pyodide, 'synthesis_qualifiers_json', [jsonArgument(payload.request, null)])
  },
  // The rules of a group color edit (OV-244): Python writes them over the whole catalog.
  [DIAGRAM_HELPER_OPERATIONS.SYNTHESIZE_COLOR_RULES]: {
    keys: ['request'],
    fileRoles: [],
    run: (pyodide, payload) => callJsonHelper(pyodide, 'synthesize_color_rules_json', [jsonArgument(payload.request, null)])
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
  // The size limit of a PNG or PDF (gbdraw/render/export_limits.py), checked before the browser draws it; for the
  // SVG of a Result drawn to a requested size, the root size its file states (gbdraw/render/document.py).
  [DIAGRAM_HELPER_OPERATIONS.CHECK_EXPORT_LIMIT]: {
    keys: ['format', 'width', 'height', 'dpi', 'linearWrap', 'figureSize'],
    fileRoles: [],
    run: (pyodide, payload) => callJsonHelper(pyodide, 'check_web_export_limit_json', [JSON.stringify({
      format: payload.format, width: payload.width, height: payload.height,
      dpi: payload.dpi ?? null, linearWrap: payload.linearWrap ?? null, figureSize: payload.figureSize ?? null
    })])
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
    keys: ['files', 'format', 'regionSpec', 'recordSelector', 'reverseFlag', 'isCircular'],
    fileRoles: ['source'],
    run: (pyodide, payload, paths, operation) => callJsonHelper(
      pyodide,
      'extract_first_fasta',
      [
        requireHelperFile(paths, 'source', operation),
        String(payload.format || 'genbank').trim().toLowerCase(),
        payload.regionSpec ?? null,
        payload.recordSelector ?? null,
        payload.reverseFlag ? '1' : '0',
        displayCircularArgument(payload.isCircular, operation)
      ]
    )
  },
  // Every record of a file as the CLI loader reads it, when the browser reader declines it (OV-1599).
  [DIAGRAM_HELPER_OPERATIONS.READ_ALL_SEQUENCE_RECORDS]: {
    keys: ['files', 'format'],
    fileRoles: ['source'],
    run: (pyodide, payload, paths, operation) => callJsonHelper(pyodide, 'read_all_sequence_records', [
      requireHelperFile(paths, 'source', operation),
      String(payload.format || 'genbank').trim().toLowerCase()
    ])
  },
  [DIAGRAM_HELPER_OPERATIONS.EXTRACT_CDS_PROTEIN_FASTA]: {
    keys: [
      'files',
      'format',
      'regionSpec',
      'recordSelector',
      'reverseFlag',
      'isCircular',
      'recordIndex',
      'recordInstanceKey',
      'featureOverrides',
      'windowOf',
      'recordSource'
    ],
    fileRoles: ['source', 'fasta', 'visibility'],
    // A record of the drawing is staged under its resource ids, so Python reads
    // the records the last Generate parsed (OV-1555); another file comes in `files`.
    run: async (pyodide, payload, paths, operation, workspace) => {
      const format = String(payload.format || 'genbank').trim().toLowerCase();
      const staged = payload.recordSource === undefined
        ? null
        : await prepareCanonicalResources(pyodide, `${workspace}/request`, payload);
      const { source, fasta: fastaPath } = staged
        ? stagedRecordSourcePaths(payload.recordSource, staged, operation)
        : { source: requireHelperFile(paths, 'source', operation), fasta: paths.get('fasta') || null };
      if (format === 'gff' && !fastaPath) {
        throw new TypeError(`Diagram helper '${operation}' requires the 'fasta' file for GFF3.`);
      }
      return callJsonHelper(pyodide, 'extract_cds_protein_fasta', [
        source,
        format,
        fastaPath,
        payload.regionSpec ?? null,
        payload.recordSelector ?? null,
        payload.reverseFlag ? '1' : '0',
        payload.recordIndex ?? null,
        payload.recordInstanceKey ?? null,
        paths.get('visibility') || null,
        // This record's canonical featureOverrides rows (design Q4, 3.3).
        payload.featureOverrides == null ? null : JSON.stringify(payload.featureOverrides),
        payload.windowOf ?? null,
        displayCircularArgument(payload.isCircular, operation),
        jsonArgument(staged, {}),
        jsonArgument(staged && buildPreparedResourceIdentityMap(payload.resourceManifest), {})
      ]);
    }
  },
  [DIAGRAM_HELPER_OPERATIONS.BUILD_PROTEIN_LOSAT_CACHE_KEYS]: {
    // `derived` (optional) asks for the key of the derived payload the raw results convert to.
    keys: ['identityManifest', 'pairs', 'derived'],
    fileRoles: [],
    run: (pyodide, payload) => callJsonHelper(
      pyodide,
      'build_protein_losat_cache_keys_json',
      [
        jsonArgument(payload.identityManifest, {}),
        jsonArgument(payload.pairs, []),
        jsonArgument(payload.derived, null)
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
        payload.collinearMergeOrientation,
        payload.collinearInferOrthogroups,
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
  // A Session up to 47 stored the protein results of a reversed Linear record
  // in the drawn frame: Python returns their bytes in the search frame by
  // stored resource ID (OV-1145), from the drawing's request and its staged files.
  [DIAGRAM_HELPER_OPERATIONS.CONVERT_STORED_PROTEIN_RESULT_FRAMES]: {
    keys: ['version', 'canonicalRequest'],
    fileRoles: [],
    run: async (pyodide, payload, _paths, _operation, workspace) => callJsonHelper(
      pyodide, 'stored_protein_results_to_search_frame_json', [
        jsonArgument({ version: payload.version, renderRequest: payload.canonicalRequest }, null),
        jsonArgument(await prepareCanonicalResources(pyodide, `${workspace}/request`, payload), {}),
        `${workspace}/output`,
        jsonArgument(buildPreparedResourceIdentityMap(payload.resourceManifest), {})
      ]
    )
  },
  // Load of a Session 27-44: Python returns its Session 51, the stored Results
  // a Session 51 cannot hold, and the Load repairs (DESIGN 4.3). The resources
  // travel inside the document, so the call names no file.
  [DIAGRAM_HELPER_OPERATIONS.UPGRADE_LEGACY_SESSION]: {
    keys: ['document', 'storedRequestSchema'],
    fileRoles: [],
    run: (pyodide, payload) => callJsonHelper(
      pyodide, 'legacy_session_to_drawings_json',
      [jsonArgument(payload.document, null), jsonArgument(payload.storedRequestSchema, null)]
    )
  },
  // Generate's ring search output in the search frame of the records its
  // request draws (PD-OI-073): Python converts it from the staged request.
  [DIAGRAM_HELPER_OPERATIONS.CONVERT_RING_TABLES_TO_SEARCH_FRAME]: {
    keys: ['resourceIds', 'canonicalRequest'],
    fileRoles: [],
    run: async (pyodide, payload, _paths, _operation, workspace) => callJsonHelper(
      pyodide, 'circular_ring_tables_in_search_frame_json', [
        jsonArgument(payload.resourceIds, []),
        jsonArgument(payload.canonicalRequest, null),
        jsonArgument(await prepareCanonicalResources(pyodide, `${workspace}/request`, payload), {}),
        `${workspace}/output`,
        jsonArgument(buildPreparedResourceIdentityMap(payload.resourceManifest), {})
      ]
    )
  },
  // The Orient records button: Python plans the last Generate's request with
  // auto_orient and returns each record's drawn orientation by record key.
  [DIAGRAM_HELPER_OPERATIONS.DECIDE_AUTO_ORIENTATION]: {
    keys: ['canonicalRequest'],
    fileRoles: [],
    run: async (pyodide, payload, _paths, _operation, workspace) => callJsonHelper(
      pyodide, 'decide_auto_orientation_json', [
        jsonArgument(payload.canonicalRequest, null),
        jsonArgument(await prepareCanonicalResources(pyodide, `${workspace}/request`, payload), {}),
        `${workspace}/output`,
        jsonArgument(buildPreparedResourceIdentityMap(payload.resourceManifest), {})
      ]
    )
  },
  // New drawing from a region: Python derives the drawing from the source
  // drawing's draft projection and returns it with its resource sources.
  [DIAGRAM_HELPER_OPERATIONS.DERIVE_REGION_DRAWING]: {
    keys: ['canonicalRequest', 'regions', 'adaptSizes'],
    fileRoles: [],
    run: async (pyodide, payload, _paths, _operation, workspace) => callJsonHelper(
      pyodide, 'derive_web_region_drawing_json', [
        jsonArgument(payload.canonicalRequest, null),
        jsonArgument(await prepareCanonicalResources(pyodide, `${workspace}/request`, payload), {}),
        `${workspace}/output`,
        jsonArgument(payload.regions, []),
        jsonArgument(payload.adaptSizes, true),
        jsonArgument(buildPreparedResourceIdentityMap(payload.resourceManifest), {})
      ]
    )
  },
  // Find in other files: one step (the reference, the search plan, the plan of its rows, the new drawing) over the
  // Files of the shown drawing's draft projection.
  [DIAGRAM_HELPER_OPERATIONS.FIND_IN_OTHER_FILES]: {
    keys: ['canonicalRequest', 'locate'],
    fileRoles: [],
    run: async (pyodide, payload, _paths, _operation, workspace) => callJsonHelper(
      pyodide, 'locate_web_json', [
        jsonArgument(payload.canonicalRequest, null),
        jsonArgument(await prepareCanonicalResources(pyodide, `${workspace}/request`, payload), {}),
        `${workspace}/output`,
        jsonArgument(payload.locate, null),
        jsonArgument(buildPreparedResourceIdentityMap(payload.resourceManifest), {})
      ]
    )
  },
  // Copy FASTA of Select region: Python slices the records of the request the shown Result was made from.
  [DIAGRAM_HELPER_OPERATIONS.READ_RECORD_SPAN_SEQUENCES]: {
    keys: ['canonicalRequest', 'spans'],
    fileRoles: [],
    run: async (pyodide, payload, _paths, _operation, workspace) => callJsonHelper(
      pyodide, 'read_record_span_sequences_json', [
        jsonArgument(payload.canonicalRequest, null),
        jsonArgument(await prepareCanonicalResources(pyodide, `${workspace}/request`, payload), {}),
        `${workspace}/output`,
        jsonArgument(payload.spans, []),
        jsonArgument(buildPreparedResourceIdentityMap(payload.resourceManifest), {})
      ]
    )
  },
  // The circled numbers and omitted lengths a Linear card shows beside its windows (PD-OI-124).
  [DIAGRAM_HELPER_OPERATIONS.PLAN_CARD_WINDOWS]: {
    keys: ['records'],
    fileRoles: [],
    run: (pyodide, payload) => callJsonHelper(pyodide, 'plan_card_windows_json', [jsonArgument(payload.records, [])])
  },
  // The columns a Source recipe's records table gives a placed file drawn whole, one per record (OV-3072).
  [DIAGRAM_HELPER_OPERATIONS.PLAN_RECORD_PLACES]: {
    keys: ['records'],
    fileRoles: [],
    run: (pyodide, payload) => callJsonHelper(pyodide, 'plan_record_places_json', [jsonArgument(payload.records, [])])
  },
  // Generate's LOSAT searches pass the E-value scale Python plans for the drawing's draft request, staged as Generate
  // stages its request (PD-OI-018 revision 6).
  [DIAGRAM_HELPER_OPERATIONS.PLAN_GENOME_SCALES]: {
    keys: ['canonicalRequest', 'program'],
    fileRoles: [],
    run: async (pyodide, payload, _paths, _operation, workspace) => callJsonHelper(
      pyodide, 'plan_genome_scales_json', [
        jsonArgument(payload.canonicalRequest, null),
        jsonArgument(await prepareCanonicalResources(pyodide, `${workspace}/request`, payload), {}),
        `${workspace}/output`,
        String(payload.program),
        // The file names a per-file warning names, as the input summary names them.
        jsonArgument(Object.fromEntries(payload.resourceManifest.map(({ resourceId, name }) => [resourceId, name])), {}),
        jsonArgument(buildPreparedResourceIdentityMap(payload.resourceManifest), {})
      ]
    )
  },
  // The sequences of features by identity (the Sequence tab, member Copy, sequence-field search), read by Python
  // from the request the shown Result was made from (D5).
  [DIAGRAM_HELPER_OPERATIONS.READ_FEATURE_SEQUENCES]: {
    keys: ['canonicalRequest', 'features'],
    fileRoles: [],
    run: async (pyodide, payload, _paths, _operation, workspace) => callJsonHelper(
      pyodide, 'read_feature_sequences_json', [
        replyBuffersOf(pyodide),
        jsonArgument(payload.canonicalRequest, null),
        jsonArgument(await prepareCanonicalResources(pyodide, `${workspace}/request`, payload), {}),
        `${workspace}/output`,
        jsonArgument(payload.features, []),
        jsonArgument(buildPreparedResourceIdentityMap(payload.resourceManifest), {})
      ]
    )
  },
  // A Circular Generate's ring searches pass the search space Python plans for the reference and each ring
  // (PD-OI-018 revision 6): the reference's record sources are staged under the ids the next Generate's request names
  // them by, so Python reads the parse that Generate renders; each ring comes as the ring reader's file lengths.
  [DIAGRAM_HELPER_OPERATIONS.PLAN_RING_SCALES]: {
    keys: ['program', 'referenceSources', 'rings'],
    fileRoles: [],
    run: async (pyodide, payload, _paths, _operation, workspace) => callJsonHelper(
      pyodide, 'plan_ring_scales_json', [
        String(payload.program),
        jsonArgument(payload.referenceSources, []),
        jsonArgument(await prepareCanonicalResources(pyodide, `${workspace}/request`, payload), {}),
        // The file names a per-file warning names, as the input summary names them.
        jsonArgument(Object.fromEntries(payload.resourceManifest.map(({ resourceId, name }) => [resourceId, name])), {}),
        jsonArgument(payload.rings, []),
        jsonArgument(buildPreparedResourceIdentityMap(payload.resourceManifest), {})
      ]
    )
  },
  // The bases of match spans of drawn or comparison records (the match popup's Copy and FASTA), on the displayed
  // strand, read by Python from the shown Result's request (D5).
  [DIAGRAM_HELPER_OPERATIONS.READ_MATCH_SEQUENCES]: {
    keys: ['canonicalRequest', 'spans'],
    fileRoles: [],
    run: async (pyodide, payload, _paths, _operation, workspace) => callJsonHelper(
      pyodide, 'read_match_sequences_json', [
        replyBuffersOf(pyodide),
        jsonArgument(payload.canonicalRequest, null),
        jsonArgument(await prepareCanonicalResources(pyodide, `${workspace}/request`, payload), {}),
        `${workspace}/output`,
        jsonArgument(payload.spans, []),
        jsonArgument(buildPreparedResourceIdentityMap(payload.resourceManifest), {})
      ]
    )
  },
  // One part of the stored catalog item of the shown Result with every sequence an exported Interactive SVG embeds,
  // in the export encoding, read from the shown Result's request (D5): the export reads it in parts of at most
  // FEATURE_READ_BASES bases of rows. The part's rows come as a file (`embedItem`), not as a string argument.
  [DIAGRAM_HELPER_OPERATIONS.READ_EMBED_SEQUENCES]: {
    keys: ['files', 'canonicalRequest', 'part'],
    fileRoles: ['embedItem'],
    run: async (pyodide, payload, paths, operation, workspace) => callJsonHelper(
      pyodide, 'read_embed_sequences_json', [
        replyBuffersOf(pyodide),
        jsonArgument(payload.canonicalRequest, null),
        jsonArgument(await prepareCanonicalResources(pyodide, `${workspace}/request`, payload), {}),
        `${workspace}/output`,
        requireHelperFile(paths, 'embedItem', operation),
        jsonArgument(payload.part, {}),
        jsonArgument(buildPreparedResourceIdentityMap(payload.resourceManifest), {})
      ]
    )
  },
  [DIAGRAM_HELPER_OPERATIONS.PLAN_LINEAR_UNITS]: {
    keys: ['records', 'pairs'],
    fileRoles: [],
    run: (pyodide, payload) => callJsonHelper(
      pyodide, 'plan_linear_units_json', [jsonArgument(payload.records, []), jsonArgument(payload.pairs, [])]
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
  // Sequence search: Python reads the pasted or chosen queries with the FASTA-only loader,
  // and turns one search's LOSAT rows into an annotation set (data, never SVG).
  [DIAGRAM_HELPER_OPERATIONS.READ_SEARCH_QUERIES]: {
    keys: ['files'],
    fileRoles: ['source'],
    run: (pyodide, _payload, paths, operation) => callJsonHelper(
      pyodide, 'read_search_queries_json', [requireHelperFile(paths, 'source', operation)]
    )
  },
  [DIAGRAM_HELPER_OPERATIONS.SEARCH_HIT_ANNOTATION_SET]: {
    keys: ['search'],
    fileRoles: [],
    run: (pyodide, payload) => callJsonHelper(
      pyodide, 'sequence_search_annotation_set_json', [jsonArgument(payload.search, null)]
    )
  },
  [DIAGRAM_HELPER_OPERATIONS.READ_COMPARISON_SEQUENCE]: {
    keys: ['files'],
    fileRoles: ['source[]'],
    run: (pyodide, _payload, paths, operation) => callJsonHelper(
      pyodide,
      'read_comparison_sequence_json',
      requireHelperFiles(paths, 'source', operation)
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
  }
});

/**
 * The payload keys a helper accepts: its own, and the staged resource keys of a
 * helper the page stages resources for. `null`: no such helper.
 * @param {string} operation
 * @returns {string[] | null}
 */
export const helperPayloadKeys = (operation) => {
  if (!Object.hasOwn(HELPER_OPERATION_SPECS, operation)) return null;
  const spec = HELPER_OPERATION_SPECS[/** @type {keyof typeof HELPER_OPERATION_SPECS} */ (operation)];
  return [...spec.keys, ...(DIAGRAM_HELPER_STAGING_OPERATIONS.includes(operation) ? STAGED_RESOURCE_PAYLOAD_KEYS : [])];
};

const runHelperOperation = async ({ operation, payload, requestId } = {}) => {
  if (!runtime?.pyodide) {
    throw serializeError({ code: 'WORKER_INIT', operation, stage: 'initialization' });
  }
  const normalizedOperation = String(operation || '').trim();
  const allowedKeys = helperPayloadKeys(normalizedOperation);
  const spec = allowedKeys ? HELPER_OPERATION_SPECS[normalizedOperation] : null;
  if (!spec || !allowedKeys) {
    throw serializeError({ code: 'HELPER_PROTOCOL', operation: normalizedOperation, stage: 'request-validation' });
  }
  let failureStage = 'request-validation';
  try {
    const normalizedPayload = requirePayloadObject(payload, normalizedOperation);
    assertAllowedPayloadKeys(normalizedPayload, normalizedOperation, allowedKeys);
    failureStage = 'resource-staging';
    return await withRequestWorkspace(
      runtime.pyodide,
      requestId,
      'helper',
      async (workspace) => {
        const paths = await stageHelperFiles(
          runtime.pyodide,
          workspace,
          normalizedPayload.files,
          spec.fileRoles
        );
        failureStage = 'helper';
        try {
          return await spec.run(runtime.pyodide, normalizedPayload, paths, normalizedOperation, workspace);
        } catch (error) {
          throw nameStagedHelperInput(error, normalizedPayload.files);
        }
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
  retainedResourceTokens,
  heldResourceTokens,
  drawingId = null,
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
      (total, entry) => total + Number(entry?.blob?.size || 0),
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
      pyodide, workspace, { resourceManifest, stagedResources, retainedResourceTokens },
      { testLifecycleEnabled, requestId, prune: true }
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
    const preparedInputRetentionJson = JSON.stringify(
      buildPreparedInputRetention(drawingId, retainedResourceTokens, heldResourceTokens)
    );
    // The file names the input summary gives the staged resources (gbdraw/api/input_summary.py).
    const resourceNamesJson = JSON.stringify(Object.fromEntries(
      resourceManifest.map(({ resourceId, name }) => [resourceId, name])
    ));
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
      JSON.stringify(heldSegments.map(({ revision }) => revision)),
      preparedInputRetentionJson,
      resourceNamesJson,
      // CAPACITY-WARN: Python calls this with the JSON text of labels that may not fit, before it places them;
      // the page words it while the render goes on. Pyodide passes the function as a JsProxy for the call.
      (/** @type {unknown} */ notice) => {
        self.postMessage({ type: 'progress', requestId, stage: 'rendering', notice: String(notice) });
      },
      replyBuffersOf(pyodide)
    );
    emitTestLifecycle(testLifecycleEnabled, requestId, 'python-wrapper-end');
    self.postMessage({ type: 'progress', requestId, stage: 'finalizing' });
    emitTestLifecycle(testLifecycleEnabled, requestId, 'result-object-conversion-start');
    failureStage = 'result-admission';
    result = convertGenerationResult(resultHandle);
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
      heldSegmentCount: heldSegments.length,
      bytes: collectGenerationResultTransferList(result).reduce(
        (total, buffer) => total + buffer.byteLength,
        0
      )
    });
  } catch (error) {
    // A Cancel stops Python with KeyboardInterrupt (the cancel message is read
    // only after this run); tests read where it stopped: the traceback's last line.
    const pythonError = /** @type {{ type?: unknown, message?: unknown }} */ (error);
    if (pythonError?.type === 'KeyboardInterrupt') {
      emitTestLifecycle(testLifecycleEnabled, requestId, 'python-keyboard-interrupt', {
        pythonError: String(pythonError.message ?? '').trim().split('\n').at(-1)
      });
    }
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
    { role: 'source', name: sourceFile?.name, blob: sourceFile?.blob },
    ...(fastaFile ? [{ role: 'fasta', name: fastaFile.name, blob: fastaFile.blob }] : []),
    ...(visibilityFile ? [{ role: 'visibility', name: visibilityFile.name, blob: visibilityFile.blob }] : [])
  ];

  let failureStage = 'resource-staging';
  try {
    return await withRequestWorkspace(pyodide, requestId, 'feature', async (workspace) => {
      const paths = await stageHelperFiles(
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
    throw nameStagedHelperInput(serializeError(error, { operation: 'feature-extraction', stage: failureStage }), workspaceFiles);
  }
};

const handleWorkerMessage = async (data) => {
  const { id, requestId, type } = data;
  try {
    if (type === 'init') {
      interruptBuffer = data.interruptBuffer instanceof Uint8Array ? data.interruptBuffer : null;
      let initialized;
      try {
        initialized = await initializeRuntime(data);
      } catch (error) {
        const snapshotError = pyodideSnapshotFailureDetail(error);
        if (snapshotError === null) throw error;
        // A snapshot restore error: the page replaces this Worker with one
        // that boots cold (services/diagram-generation.js).
        self.postMessage({ id, type: 'init', ok: false, snapshotError });
        return;
      }
      // One site for both boots: a Cancel's SIGINT reaches the handler that
      // _start_web_worker_state installed (D7).
      if (interruptBuffer) initialized.pyodide.setInterruptBuffer(interruptBuffer);
      // The replay of the page's fonts into a new Worker, font by font; the
      // reply names each font it lost (OV-686), and a Generate that names one
      // fails (FONT_NOT_LOADED).
      const lostFonts = await fontRegistryFor(initialized.pyodide).replay(data.fonts);
      self.postMessage({
        id,
        type: 'init',
        ok: true,
        capabilities: initialized.capabilities,
        lostFonts,
        boot: initialized.boot
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
    if (isCanceled(requestId)) {
      noteDroppedReply(results?.error);
      return;
    }
    self.postMessage(
      { requestId, type: 'run', ok: true, results },
      collectGenerationResultTransferList(results)
    );
  } catch (error) {
    const serialized = serializeError(error, {
      operation: type === 'run' ? 'generate' : type === 'feature-extraction' ? 'feature-extraction' : data.operation,
      stage: type === 'init' ? 'initialization' : 'unknown',
      code: type === 'init' ? 'WORKER_INIT' : 'UNKNOWN'
    });
    if (isCanceled(requestId)) {
      noteDroppedReply(serialized);
      return;
    }
    self.postMessage({ id, requestId, type: type || 'run', ok: false, error: serialized });
  }
};

/** @param {() => Promise<void> | void} step */
const enqueue = (step) => {
  const scheduled = operationQueue.then(step, step);
  operationQueue = scheduled.catch(() => {});
};

// The out-of-memory error of a canceled request's dropped reply (OV-3027): the
// drain answer carries it, so the page still replaces this engine.
/** @type {unknown} */
let droppedOutOfMemory = null;
/** @param {unknown} error The error of a reply this Worker drops. */
export const noteDroppedReply = (error) => {
  if (error && typeof error === 'object' && /** @type {{ code?: unknown }} */ (error).code === OUT_OF_MEMORY_CODE) {
    droppedOutOfMemory = error;
  }
};
/**
 * The answer to a Cancel's drain; it carries a dropped out-of-memory error once.
 * @param {string} id
 * @returns {import('../services/diagram-worker-protocol.js').DiagramWorkerCancelReply}
 */
export const cancelDrainReply = (id) => {
  const error = droppedOutOfMemory;
  droppedOutOfMemory = null;
  return error === null ? { type: 'cancel', id, ok: true } : { type: 'cancel', id, ok: true, error };
};

// The drain after a Cancel, once the operation in flight has ended: the
// interrupt and then Python's cancel flag are cleared, so the next request
// runs; the page sends none before this answer.
const drainCancel = (id) => {
  if (interruptBuffer) Atomics.store(interruptBuffer, 0, 0);
  if (runtime) callJsonHelper(runtime.pyodide, 'clear_web_cancel', []);
  self.postMessage(cancelDrainReply(id));
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
  if (data.type === 'cancel') {
    canceledThrough = Math.max(canceledThrough, Number(data.through) || 0);
    if (auxiliaryAcknowledgement && isCanceled(auxiliaryAcknowledgement.requestId)) {
      const pending = auxiliaryAcknowledgement;
      auxiliaryAcknowledgement = null;
      pending.resolve();
    }
    enqueue(() => drainCancel(data.id));
    return;
  }
  enqueue(() => handleWorkerMessage(data));
};
