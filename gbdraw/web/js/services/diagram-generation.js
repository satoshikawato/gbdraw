// @ts-check
import { createBoundedJsonReceiver } from './bounded-json-transport.js';
import { buildPyodideAssetManifest } from './pyodide-assets.js';
import { normalizeUserFacingError } from '../utils/error-normalization.js';
import { validateWebRuntimeCapabilities, DiagramRuntimeCompatibilityError } from './runtime-capabilities.js';
import {
  DIAGRAM_HELPER_OPERATION_NAMES,
  DIAGRAM_WORKER_INTERRUPT_SIGNAL,
  stagedHelperReferences,
  readFeatureCatalogRevision,
  readSvgDocumentParts
} from './diagram-worker-protocol.js';
import {
  recordSessionLifecycleEvent,
  recordStructuralMetric,
  runtimeTestHooksEnabled
} from './runtime-test-hooks.js';
import { createDiagramResourceTransport } from './diagram-resource-staging.js';
import { markCurrentWorkerGenerationResponse, SEGMENTS } from './current-worker-result-source.js';
import { normalizeLabelFitFailure } from './label-fit.js';
/** @import { LabelFitFailure } from './label-fit.js' */
import { requestFontFamily, webFontIdentitiesFor } from './font-store.js';
import { toRaw } from '../utils/vue-raw.js';

/** @import { UserFacingOptions } from '../utils/error-normalization.js' */
/** @import { ReadCommittedResourceSets, ResourceDescriptor } from './diagram-resource-staging.js' */
/** @import { CanonicalRecordSource } from '../contracts/request.generated.js' */
/** @import { BoundedJsonValue } from './bounded-json-transport.js' */
/** @import { WebFontEntry } from './font-store.js' */
/** @import { FeatureCatalog } from './feature-catalog.js' */
/** @import { FeatureCatalogRevision, WorkerHeldSegment } from './diagram-worker-protocol.js' */
/** @import { DiagramWorkerCancelMessage, DiagramWorkerCancelReply, DiagramWorkerInitMessage } from './diagram-worker-protocol.js' */
/** @import { DecodedContentPiece, DecodedSvgDocument, HeldSegment } from './current-worker-result-source.js' */

/**
 * The fingerprint of the artifact a Generate returned.
 * @typedef {object} GeneratedArtifactIdentity
 * @property {1} schema
 * @property {'SHA-256'} algorithm
 * @property {string} fingerprint
 * @property {number} retainedBytes
 * @property {number} resultBytes
 * @property {number} metadataBytes
 */

/**
 * What the Worker's Generate message carries once decoded.
 * @typedef {object} NormalizedGeneration
 * @property {unknown} results The decoded Result entries, or the engine's error payload.
 * @property {Record<string, unknown>} metadata
 * @property {GeneratedArtifactIdentity | null} artifactIdentity
 * @property {LabelFitFailure} [labelFit] With an engine error of labels that do not fit: Python's report and its row token.
 */

/** @typedef {NormalizedGeneration & { requestId: number }} GenerationResponse What a Generate request resolves with. */

/**
 * What the Worker sends in `results` of a Generate message, undecoded.
 * @typedef {object} WorkerGenerationPayload
 * @property {unknown} [error]
 * @property {unknown} [results]
 * @property {unknown} [metadata]
 * @property {unknown} [featureCatalog] The feature catalog's UTF-8 JSON bytes; absent when Python reports it unchanged.
 * @property {unknown} [featureCatalogRevision] The catalog revision with its digest, or null.
 * @property {unknown} [labelFit] The JSON text of a label fit report, next to an error.
 * @property {Partial<GeneratedArtifactIdentity> | null} [artifactIdentity]
 */

/**
 * A Result entry as the Worker sends it: its `document` pieces and the
 * `segments` buffer of the children texts it sent, or (a test, a caller of
 * the array shape) its `content` text.
 * @typedef {{ document?: unknown, segments?: unknown, content?: unknown, [key: string]: unknown }} WorkerResultEntry
 */

/**
 * A request for Generate: the canonical render request and the resources it names.
 * @typedef {object} GenerationPayload
 * @property {unknown} [request]
 * @property {Readonly<Record<string, Partial<ResourceDescriptor> | null | undefined>> | null} [resources]
 */

/**
 * A staged file the Worker receives by transfer.
 * @typedef {{ bytes?: unknown, [key: string]: unknown } | null | undefined} TransferFile
 */

/**
 * The payload of a feature-extraction or helper request.
 * @typedef {{ files?: TransferFile[], stagedResources?: TransferFile[], [key: string]: unknown }} AuxiliaryPayload
 */

/**
 * @typedef {object} HelperPayload What a helper request may carry for resource staging; its other keys are the operation's own.
 * @property {GenerationPayload['resources']} [resources]
 * @property {unknown} [canonicalRequest]
 * @property {{ canonicalRequest?: unknown } | null} [projection]
 * @property {CanonicalRecordSource} [recordSource] The record source of the drawing the blastp extraction reads.
 */

/**
 * A progress report of a request.
 * @typedef {{ requestId: number, stage: string } & Partial<RuntimePreparationStep>} RequestProgress
 */

/**
 * What a feature-extraction or helper request resolves with. The requests
 * relay the Worker's JSON through a promise that stays untyped: each operation
 * has its own result shape, which its callers read.
 * @typedef {{ requestId: number, result: BoundedJsonValue }} AuxiliaryResponse
 */

export { DIAGRAM_HELPER_OPERATIONS } from './diagram-worker-protocol.js';
export { resourceIdOwnedBy } from './diagram-resource-staging.js';

export class DiagramGenerationCanceledError extends Error {
  constructor(message = 'Diagram generation was canceled.') {
    super(message);
    this.name = 'DiagramGenerationCanceledError';
    this.canceled = true;
  }
}

/** @param {{ canceled?: unknown, name?: unknown } | null | undefined} error */
export const isDiagramGenerationCanceled = (error) =>
  Boolean(error?.canceled || error?.name === 'DiagramGenerationCanceledError');

/** @param {unknown} value */
const viewBytes = (value) => (
  value instanceof ArrayBuffer
    ? new Uint8Array(value)
    : ArrayBuffer.isView(value)
      ? new Uint8Array(value.buffer, value.byteOffset, value.byteLength)
      : null
);

/**
 * Decode one Result of a Generate reply into its admission candidate: the
 * document with each sent children text decoded once, and under `SEGMENTS`
 * the segment of the request's snapshot that each held piece names. A held
 * piece the snapshot lacks, or names with another length, is an invalid
 * reply. A Result with `content` text passes unchanged.
 * @param {WorkerResultEntry | string | number | boolean | null | undefined} result
 * @param {number} resultIndex
 * @param {ReadonlyMap<string, Readonly<HeldSegment>>} heldSegments The request's snapshot.
 */
const decodeGenerationResult = (result, resultIndex, heldSegments) => {
  if (!result || typeof result !== 'object' || !Object.hasOwn(result, 'document')) return result;
  const { document, segments, ...rest } = result;
  const buffer = segments === undefined ? new Uint8Array(0) : viewBytes(segments);
  if (!buffer) throw invalidGenerationResult();
  const parts = readSvgDocumentParts(document, buffer.byteLength, invalidGenerationResult);
  const decoder = new TextDecoder();
  /** @type {Map<string, Readonly<HeldSegment>>} */
  const held = new Map();
  if (segments !== undefined) {
    recordSessionLifecycleEvent('result-binary-decode-start', { resultIndex, bytes: buffer.byteLength });
  }
  const pieces = parts.pieces.map((piece) => {
    if (piece.kind === 'raw') return piece;
    const { start, revision, bytes } = piece;
    if (piece.offset === undefined) {
      const segment = heldSegments.get(revision);
      if (!segment || segment.bytes !== bytes) throw invalidGenerationResult();
      held.set(revision, segment);
      return /** @type {Readonly<DecodedContentPiece>} */ (Object.freeze({ kind: 'content', start, revision, bytes }));
    }
    const children = decoder.decode(buffer.subarray(piece.offset, piece.offset + bytes));
    return /** @type {Readonly<DecodedContentPiece>} */ (Object.freeze({ kind: 'content', start, revision, bytes, children }));
  });
  if (segments !== undefined) {
    recordStructuralMetric('resultBinaryDecodeCount', 1, { resultIndex });
    recordStructuralMetric('resultBinaryDecodedBytes', buffer.byteLength, { resultIndex });
    recordSessionLifecycleEvent('result-binary-decode-end', { resultIndex, bytes: buffer.byteLength });
  }
  /** @type {Readonly<DecodedSvgDocument>} */
  const decoded = Object.freeze({ head: parts.head, tail: parts.tail, pieces: Object.freeze(pieces) });
  return { ...rest, document: decoded, [SEGMENTS]: held };
};

/**
 * @param {unknown} results
 * @param {ReadonlyMap<string, Readonly<HeldSegment>>} heldSegments
 */
const decodeGenerationResults = (results, heldSegments) => (
  Array.isArray(results)
    ? results.map((result, resultIndex) => decodeGenerationResult(result, resultIndex, heldSegments))
    : results
);

/**
 * The two JSON payloads of a Generate reply and the names their decoding reports.
 * @typedef {{ event: string, countMetric: string, bytesMetric: string }} GenerationJsonPart
 */
/** @type {Readonly<Record<'metadata' | 'featureCatalog', GenerationJsonPart>>} */
const GENERATION_JSON_PARTS = Object.freeze({
  metadata: Object.freeze({
    event: 'result-metadata-binary-decode',
    countMetric: 'resultMetadataBinaryDecodeCount',
    bytesMetric: 'resultMetadataBinaryDecodedBytes'
  }),
  featureCatalog: Object.freeze({
    event: 'result-feature-catalog-binary-decode',
    countMetric: 'resultFeatureCatalogBinaryDecodeCount',
    bytesMetric: 'resultFeatureCatalogBinaryDecodedBytes'
  })
});

/**
 * @param {unknown} value The Worker's value: an object, or its UTF-8 JSON bytes.
 * @param {GenerationJsonPart} part
 * @returns {unknown}
 */
const decodeGenerationJson = (value, { event, countMetric, bytesMetric }) => {
  const bytes = viewBytes(value);
  if (!bytes) return value;
  recordSessionLifecycleEvent(`${event}-start`, {
    bytes: bytes.byteLength
  });
  const decoded = JSON.parse(new TextDecoder().decode(bytes));
  recordStructuralMetric(countMetric);
  recordStructuralMetric(bytesMetric, bytes.byteLength);
  recordSessionLifecycleEvent(`${event}-end`, {
    bytes: bytes.byteLength
  });
  return decoded;
};

/**
 * The revision each feature catalog this page received arrived with, and its
 * items' Result names then (before a rerender renames them, OV-136). Keyed by
 * the object the page owns: no copy, and the record dies with the object.
 * @type {WeakMap<object, { revision: Readonly<FeatureCatalogRevision>, resultNames: readonly string[] }>}
 */
const featureCatalogRevisions = new WeakMap();

/**
 * @param {FeatureCatalog} catalog
 * @returns {string[]}
 */
const itemResultNames = (catalog) => catalog.items.map((item) => String(item?.resultName ?? ''));

/**
 * A feature catalog the page received, with the revision it arrived with.
 * @typedef {{ catalog: FeatureCatalog, revision: Readonly<FeatureCatalogRevision> }} HeldFeatureCatalog
 */

/**
 * The catalog the page may name to the Worker for the drawing it renders: one
 * it received from this Worker, whose items still carry the Result names the
 * render expects. A rerender expects the names its Results keep (`resultNames`,
 * OV-136); a Generate expects the engine's names, recorded at receipt.
 * @param {FeatureCatalog | null | undefined} catalog
 * @param {readonly string[] | null} resultNames
 * @returns {HeldFeatureCatalog | null}
 */
const selectHeldFeatureCatalog = (catalog, resultNames) => {
  const record = catalog && typeof catalog === 'object' ? featureCatalogRevisions.get(catalog) : undefined;
  if (!catalog || !record) return null;
  const expected = resultNames ?? record.resultNames;
  const names = itemResultNames(catalog);
  if (names.length !== expected.length || names.some((name, index) => name !== expected[index])) return null;
  return { catalog, revision: record.revision };
};

/**
 * What the page holds for one render of the shown drawing, taken once when the
 * request starts: the feature catalog it may name, and the sanitized children
 * texts its committed Results keep. The reply is resolved against this
 * snapshot only, so a commit while the request runs changes nothing it names.
 * @typedef {object} HeldRenderState
 * @property {HeldFeatureCatalog | null} catalog
 * @property {ReadonlyMap<string, Readonly<HeldSegment>>} segments
 */

/** @type {Readonly<HeldRenderState>} */
const NOTHING_HELD = Object.freeze({ catalog: null, segments: new Map() });

/**
 * The held render state of the shown drawing: its catalog (as
 * `selectHeldFeatureCatalog` above) and the segments of its committed Results.
 * @param {FeatureCatalog | null | undefined} catalog
 * @param {readonly unknown[] | null | undefined} results The drawing's committed Results.
 * @param {readonly string[] | null} resultNames
 * @returns {Readonly<HeldRenderState>}
 */
export const selectHeldRenderState = (catalog, results, resultNames) => {
  /** @type {Map<string, Readonly<HeldSegment>>} */
  const segments = new Map();
  // A reactive Result would hand out a reactive Map for its symbol key.
  for (const result of toRaw(results) ?? []) {
    const raw = result && typeof result === 'object' ? toRaw(result) : null;
    /** @type {unknown} */
    const held = raw ? /** @type {{ [SEGMENTS]?: unknown }} */ (raw)[SEGMENTS] : undefined;
    if (held instanceof Map) held.forEach((segment, revision) => segments.set(revision, segment));
  }
  return Object.freeze({ catalog: selectHeldFeatureCatalog(catalog, resultNames), segments });
};

/**
 * The held segments a run message names: revisions and byte counts only, sorted.
 * @param {Readonly<HeldRenderState>} held
 * @returns {WorkerHeldSegment[]}
 */
const heldSegmentList = (held) => [...held.segments.values()]
  .map(({ revision, bytes }) => ({ revision, bytes }))
  .sort((left, right) => (left.revision < right.revision ? -1 : left.revision > right.revision ? 1 : 0));

const invalidGenerationResult = () => deserializeWorkerError({
  code: 'RESULT_INVALID', operation: 'generate', stage: 'result-admission'
});

/**
 * @param {Partial<GeneratedArtifactIdentity> | null | undefined} identity
 * @returns {Readonly<GeneratedArtifactIdentity> | null} null for an identity that fails the schema.
 */
const normalizeGeneratedArtifactIdentity = (identity) => {
  const fingerprint = String(identity?.fingerprint || '').toLowerCase();
  if (
    identity?.schema !== 1
    || identity?.algorithm !== 'SHA-256'
    || !/^[0-9a-f]{64}$/.test(fingerprint)
  ) return null;
  const retainedBytes = Number(identity.retainedBytes);
  const resultBytes = Number(identity.resultBytes);
  const metadataBytes = Number(identity.metadataBytes);
  if (
    !Number.isSafeInteger(retainedBytes) || retainedBytes < 0
    || !Number.isSafeInteger(resultBytes) || resultBytes < 0
    || !Number.isSafeInteger(metadataBytes) || metadataBytes < 0
  ) return null;
  return Object.freeze({
    schema: 1,
    algorithm: 'SHA-256',
    fingerprint,
    retainedBytes,
    resultBytes,
    metadataBytes
  });
};

/**
 * @param {unknown[] | WorkerGenerationPayload | string | number | boolean | null | undefined} payload
 * @param {Readonly<HeldRenderState>} [held] What the request named: an "unchanged" reply keeps its
 *   catalog, and a held piece takes its segment.
 * @returns {NormalizedGeneration}
 */
export const normalizeGenerationResponse = (payload, held = NOTHING_HELD) => {
  if (Array.isArray(payload)) {
    return markCurrentWorkerGenerationResponse({
      results: decodeGenerationResults(payload, held.segments),
      metadata: {},
      artifactIdentity: null
    });
  }
  if (!payload || typeof payload !== 'object') {
    return markCurrentWorkerGenerationResponse({
      results: [],
      metadata: {},
      artifactIdentity: null
    });
  }
  if (payload.error) {
    const labelFit = normalizeLabelFitFailure(payload.labelFit);
    return markCurrentWorkerGenerationResponse({
      results: payload,
      metadata: {},
      artifactIdentity: null,
      ...(labelFit ? { labelFit } : {})
    });
  }
  const results = Array.isArray(payload.results)
    ? decodeGenerationResults(payload.results, held.segments)
    : [];
  const decodedMetadata = decodeGenerationJson(payload.metadata, GENERATION_JSON_PARTS.metadata);
  // A copy: the catalog joins it below, and a decoded object is the caller's.
  /** @type {Record<string, unknown>} */
  const metadata = (
    decodedMetadata &&
    typeof decodedMetadata === 'object' &&
    !Array.isArray(decodedMetadata)
  )
    ? { ...decodedMetadata }
    : {};
  // The catalog travels beside the metadata, never inside it (step 3a).
  if (Object.hasOwn(metadata, 'featureCatalog')) throw invalidGenerationResult();
  const replyRevision = readFeatureCatalogRevision(payload.featureCatalogRevision, invalidGenerationResult);
  if (payload.featureCatalog == null) {
    // "Unchanged": only for the revision this request named.
    const heldCatalog = held.catalog;
    if (!replyRevision || !heldCatalog || heldCatalog.revision.revision !== replyRevision.revision) {
      throw invalidGenerationResult();
    }
    metadata.featureCatalog = heldCatalog.catalog;
    recordStructuralMetric('featureCatalogReuseCount');
    recordSessionLifecycleEvent('catalog.reused', {
      resultCount: Array.isArray(payload.results) ? payload.results.length : 0
    });
  } else {
    const catalog = decodeGenerationJson(payload.featureCatalog, GENERATION_JSON_PARTS.featureCatalog);
    if (replyRevision && catalog && typeof catalog === 'object' && Array.isArray(/** @type {{ items?: unknown }} */ (catalog).items)) {
      const received = /** @type {FeatureCatalog} */ (catalog);
      featureCatalogRevisions.set(received, {
        revision: replyRevision,
        resultNames: Object.freeze(itemResultNames(received))
      });
    }
    metadata.featureCatalog = catalog;
  }
  return markCurrentWorkerGenerationResponse({
    results,
    metadata,
    artifactIdentity: normalizeGeneratedArtifactIdentity(payload.artifactIdentity)
  });
};

const resolveWorkerUrl = () => new URL('../workers/diagram-generation-worker.js', import.meta.url).toString();

/** @type {Worker | null} */
let worker = null;
let workerInitialized = false;
/** @type {ReturnType<typeof validateWebRuntimeCapabilities> | null} */
let workerCapabilities = null;
/**
 * Why a Worker of this page could not restore the Pyodide snapshot; null
 * while none failed. The Workers started after a failure boot cold.
 * @type {string | null}
 */
let pyodideSnapshotError = null;
/**
 * The interrupt buffer of the current Worker (index 0, read by Pyodide), or
 * null on a page that is not cross-origin isolated, where a Cancel
 * terminates the Worker.
 * @type {Uint8Array | null}
 */
let workerInterruptBuffer = null;
/**
 * How long a Cancel waits for the kept Worker to answer (D7) before it
 * terminates it, as a Cancel without an interrupt buffer does.
 */
export const DIAGRAM_WORKER_CANCEL_TIMEOUT_MS = 3000;
/**
 * A request sent after a Cancel, waiting for the kept Worker's drain.
 * @typedef {{ resolve: () => void, reject: (error: unknown) => void }} DrainWaiter
 */
/**
 * The kept Worker's drain after a Cancel. It finishes when the Worker
 * answers, or when the Worker is terminated: its waiters then go on (a new
 * Worker boots for them), unless the termination carries an error (dispose),
 * which rejects them, as it rejects a pending initialization.
 * @typedef {object} WorkerDrain
 * @property {Set<DrainWaiter>} waiters
 * @property {(outcome: 'drained' | 'terminated', error?: unknown) => void} finish
 */
/** @type {WorkerDrain | null} */
let drainState = null;
/** @typedef {Readonly<{ step: number, total: number, label: string }>} RuntimePreparationStep */
/**
 * The pending Worker initialization: its promise, the reject that settles it,
 * the listener cleanup, and the preparation step listeners with the last
 * step reported. Null while no initialization is in flight.
 * @type {{ promise: Promise<Worker>, reject: (error: unknown) => void, cleanup: () => void,
 *   lastStep: RuntimePreparationStep | null, stepListeners: Set<(step: RuntimePreparationStep) => void> } | null}
 */
let initState = null;
/**
 * @typedef {object} GenerateRequest The one running Generate request.
 * @property {number} requestId
 * @property {boolean} settled
 * @property {(() => void) | null} cleanup Null until its Worker listeners exist.
 * @property {(value: GenerationResponse) => void} resolve
 * @property {(error: unknown) => void} reject
 */
/** @type {GenerateRequest | null} */
let activeRequest = null;
/**
 * @typedef {object} AuxiliaryRequest A feature-extraction or helper request
 *   that waits on the Worker.
 * @property {number} requestId
 * @property {(() => void) | null} cleanup Null until its listeners exist.
 * @property {(error: unknown) => void} reject
 */
/**
 * @typedef {object} PreparedResources What `resourceTransport.prepare` stages
 *   for the Worker.
 * @property {unknown[]} resourceManifest
 * @property {TransferFile[]} stagedResources
 * @property {readonly string[]} retainedResourceTokens The tokens the Worker keeps besides the manifest's.
 * @property {() => boolean} commit Acknowledges the staging; false when a reset
 *   or another acknowledged request superseded it.
 */
/** @type {Set<AuxiliaryRequest>} */
const activeFeatureRequests = new Set();
/** @type {Set<AuxiliaryRequest>} */
const activeHelperRequests = new Set();
let nextRequestId = 1;
/** @type {ReadCommittedResourceSets} */
const NO_COMMITTED_RESOURCES = () => [];
let readCommittedResourceSets = NO_COMMITTED_RESOURCES;
const resourceTransport = createDiagramResourceTransport({
  readCommittedResourceSets: () => readCommittedResourceSets()
});

/**
 * Name the committed request of every drawing, so the Worker keeps the
 * resources they name across Generates of other drawings.
 * @param {ReadCommittedResourceSets | null} read Null when the app unmounts.
 */
export const setDiagramResourceRetainSource = (read) => {
  readCommittedResourceSets = read || NO_COMMITTED_RESOURCES;
};

/** @type {ReadonlySet<string>} */
const diagramHelperOperationNames = new Set(DIAGRAM_HELPER_OPERATION_NAMES);

/**
 * A font the replay into a new Worker could not register (OV-686).
 * @typedef {{ identity: string, label: string }} LostWorkerFont
 */

/**
 * Where the Worker's fonts come from: the page's font store, and who is told
 * of a font a new Worker could not register again.
 * @typedef {object} DiagramWorkerFontSource
 * @property {() => readonly WebFontEntry[]} readFonts
 * @property {(lost: readonly LostWorkerFont[]) => void} onFontsLost
 */

/** @type {DiagramWorkerFontSource} */
const NO_FONT_SOURCE = Object.freeze({ readFonts: () => [], onFontsLost: () => {} });
/** @type {DiagramWorkerFontSource} */
let workerFontSource = NO_FONT_SOURCE;

/**
 * Name the fonts registered in the page (the font store's), so a new Worker
 * registers them all at start and a Generate names those its request uses.
 * @param {DiagramWorkerFontSource | null} source Null when the app unmounts.
 */
export const setDiagramWorkerFontSource = (source) => {
  workerFontSource = source || NO_FONT_SOURCE;
};

// The store's bytes: postMessage makes the one copy the Worker needs.
/** @param {readonly WebFontEntry[]} fonts */
const workerFontFiles = (fonts) => fonts.map((font) => ({ identity: font.identity, label: font.label, bytes: font.bytes }));

/** @param {unknown} lost @returns {LostWorkerFont[]} */
const lostWorkerFonts = (lost) => (Array.isArray(lost) ? lost : []).map((font) => ({
  identity: String(font?.identity || ''), label: String(font?.label || '')
}));

/**
 * The page's fonts a canonical request's text font names.
 * @param {unknown} request
 */
const requestFontIdentities = (request) => {
  const fontFamily = requestFontFamily(request);
  return fontFamily ? webFontIdentitiesFor(workerFontSource.readFonts(), fontFamily) : [];
};

// A Cancel can stop Python only through memory the page shares with the Worker.
const createInterruptBuffer = () => (
  globalThis.crossOriginIsolated === true && typeof SharedArrayBuffer === 'function'
    ? new Uint8Array(new SharedArrayBuffer(1))
    : null
);

/**
 * The init message of a new Worker, which takes the Worker's interrupt buffer.
 * @param {string} id
 * @returns {DiagramWorkerInitMessage}
 */
const buildInitPayload = (id) => {
  workerInterruptBuffer = createInterruptBuffer();
  return {
    type: 'init',
    id,
    fonts: workerFontFiles(workerFontSource.readFonts()),
    ...buildPyodideAssetManifest(),
    ...(pyodideSnapshotError === null ? {} : { pyodideSnapshotUrl: null }),
    interruptBuffer: workerInterruptBuffer
  };
};

/**
 * @param {string | Record<string, unknown> | null | undefined} serialized
 * @param {UserFacingOptions} [options]
 */
export const deserializeWorkerError = (
  serialized,
  options = {}
) => {
  const compatibility = serialized instanceof DiagramRuntimeCompatibilityError;
  // normalizeUserFacingError returns null only for a falsy value; the argument
  // here is `serialized || { code: 'UNKNOWN' }` or a literal, never falsy.
  const model = /** @type {NonNullable<ReturnType<typeof normalizeUserFacingError>>} */ (normalizeUserFacingError(compatibility ? { code: 'RUNTIME_INCOMPATIBLE', stage: 'initialization' } : serialized || { code: 'UNKNOWN' }, options));
  return Object.assign(compatibility ? new DiagramRuntimeCompatibilityError('Capability contract mismatch.') : new Error(model.summary), model, { message: model.summary, stack: model.summary });
};

const getWorker = () => {
  if (!worker) {
    worker = new Worker(resolveWorkerUrl(), { type: 'module' });
    recordStructuralMetric('workerConstructionCount');
    workerInitialized = false;
  }
  return worker;
};

const clearInitState = () => {
  if (initState?.cleanup) initState.cleanup();
  initState = null;
};

/** @param {unknown} error */
const rejectPendingInit = (error) => {
  if (!initState) return;
  const { reject } = initState;
  clearInitState();
  reject(error);
};

/**
 * Reject every feature-extraction and helper request that waits on the Worker.
 * @param {unknown} error
 */
const rejectWorkerRequests = (error) => {
  for (const requests of [activeFeatureRequests, activeHelperRequests]) {
    requests.forEach((request) => {
      request.cleanup?.();
      request.reject(error);
    });
    requests.clear();
  }
};

/** @param {unknown} [error] */
const terminateWorker = (error = null) => {
  if (error) rejectPendingInit(error);
  else clearInitState();
  rejectWorkerRequests(error || new Error('Diagram generation worker was terminated.'));
  if (worker) {
    worker.terminate();
    worker = null;
  }
  workerInitialized = false;
  workerCapabilities = null;
  workerInterruptBuffer = null;
  resourceTransport.reset();
  drainState?.finish('terminated', error);
};

/**
 * Stop the kept Worker's canceled requests (D7): write SIGINT into its
 * interrupt buffer, so Python stops at its next checkpoint, and ask it to
 * drain: it drops the replies of every request up to the last one sent and
 * answers once the operation in flight ends. Without an answer in time, it
 * is terminated.
 * @param {Worker} currentWorker
 * @param {Uint8Array} buffer
 */
const drainWorker = (currentWorker, buffer) => {
  const id = `diagram-cancel-${Date.now()}-${nextRequestId}`;
  /** @type {Set<DrainWaiter>} */
  const waiters = new Set();
  /** @param {MessageEvent} event */
  const handleMessage = (event) => {
    const data = /** @type {Partial<DiagramWorkerCancelReply>} */ (event.data || {});
    if (data.type === 'cancel' && data.id === id) finish('drained');
  };
  const terminate = () => {
    if (worker === currentWorker) terminateWorker();
  };
  const timer = setTimeout(terminate, DIAGRAM_WORKER_CANCEL_TIMEOUT_MS);
  /** @param {'drained' | 'terminated'} outcome @param {unknown} [error] */
  const finish = (outcome, error = null) => {
    clearTimeout(timer);
    currentWorker.removeEventListener('message', handleMessage);
    currentWorker.removeEventListener('error', terminate);
    drainState = null;
    recordStructuralMetric('workerCancelDrain', 1, { outcome });
    for (const waiter of waiters) {
      if (error) waiter.reject(error);
      else waiter.resolve();
    }
  };
  currentWorker.addEventListener('message', handleMessage);
  currentWorker.addEventListener('error', terminate);
  drainState = { waiters, finish };
  Atomics.store(buffer, 0, DIAGRAM_WORKER_INTERRUPT_SIGNAL);
  /** @type {DiagramWorkerCancelMessage} */
  const message = { type: 'cancel', id, through: nextRequestId - 1 };
  currentWorker.postMessage(message);
};

/**
 * Wait for the kept Worker's drain, if one runs (D7). A request checks that
 * it is still live after this wait, before it boots or uses a Worker.
 * @returns {Promise<void>}
 */
const awaitDrain = () => new Promise((resolve, reject) => {
  if (drainState) drainState.waiters.add({ resolve, reject });
  else resolve();
});

/**
 * Reject the requests that wait for the drain: a second Cancel cancels them
 * too, as a Cancel rejects a request that waits for an initialization.
 * @param {unknown} error
 */
const rejectDrainWaiters = (error) => {
  const waiters = [...(drainState?.waiters ?? [])];
  drainState?.waiters.clear();
  waiters.forEach((waiter) => waiter.reject(error));
};

/**
 * @param {((step: RuntimePreparationStep) => void) | null} [onStep]
 *   Each first-use preparation step as the Worker starts it; a caller that
 *   joins a running initialization first receives the current step.
 * @returns {Promise<Worker>}
 */
const ensureWorkerInitialized = (onStep = null) => {
  const currentWorker = getWorker();
  if (workerInitialized) return Promise.resolve(currentWorker);
  if (initState?.promise) {
    if (typeof onStep === 'function') {
      initState.stepListeners.add(onStep);
      if (initState.lastStep) onStep(initState.lastStep);
    }
    return initState.promise;
  }

  const id = `diagram-init-${Date.now()}-${nextRequestId}`;
  // The Promise executor runs synchronously and assigns both before anything reads them.
  /** @type {(worker: Worker) => void} */
  let resolveInit = () => {};
  /** @type {(error: unknown) => void} */
  let rejectInit = () => {};
  /** @type {Promise<Worker>} */
  const promise = new Promise((resolve, reject) => {
    resolveInit = resolve;
    rejectInit = reject;
  });

  const cleanup = () => {
    currentWorker.removeEventListener('message', handleMessage);
    currentWorker.removeEventListener('error', handleError);
    currentWorker.removeEventListener('messageerror', handleMessageError);
  };

  /** @param {unknown} error */
  const fail = (error) => {
    cleanup();
    initState = null;
    if (worker === currentWorker) terminateWorker();
    rejectInit(error);
  };

  /** @param {MessageEvent} event */
  function handleMessage(event) {
    const data = event.data || {};
    if (data.type === 'init-progress' && data.id === id && initState) {
      const step = Object.freeze({ step: Number(data.step) || 0, total: Number(data.total) || 0, label: String(data.label || '') });
      initState.lastStep = step;
      initState.stepListeners.forEach((listener) => listener(step));
      return;
    }
    if (data.type !== 'init' || data.id !== id) return;
    const stepListeners = initState?.stepListeners;
    cleanup();
    initState = null;
    if (!data.ok && typeof data.snapshotError === 'string') {
      // The snapshot restore failed: this Worker is discarded and a new one
      // boots cold, for this initialization and every later one.
      pyodideSnapshotError = data.snapshotError || 'unknown';
      if (worker === currentWorker) terminateWorker();
      ensureWorkerInitialized((step) => stepListeners?.forEach((listener) => listener(step)))
        .then(resolveInit, rejectInit);
      return;
    }
    if (data.ok) {
      try {
        workerCapabilities = validateWebRuntimeCapabilities(data.capabilities);
      } catch (error) {
        fail(error);
        return;
      }
      workerInitialized = true;
      recordStructuralMetric('workerBoot', 1, {
        boot: pyodideSnapshotError === null ? String(data.boot) : `cold after snapshot error: ${pyodideSnapshotError}`
      });
      const lost = lostWorkerFonts(data.lostFonts);
      if (lost.length > 0) workerFontSource.onFontsLost(lost);
      resolveInit(currentWorker);
      return;
    }
    fail(deserializeWorkerError(data.error, { code: 'WORKER_INIT', stage: 'initialization' }));
  }

  function handleError() {
    fail(deserializeWorkerError({ code: 'WORKER_INIT', stage: 'initialization' }));
  }

  function handleMessageError() {
    fail(deserializeWorkerError({ code: 'WORKER_INIT', stage: 'initialization' }));
  }

  currentWorker.addEventListener('message', handleMessage);
  currentWorker.addEventListener('error', handleError);
  currentWorker.addEventListener('messageerror', handleMessageError);
  initState = { promise, reject: rejectInit, cleanup, lastStep: null,
    stepListeners: new Set(typeof onStep === 'function' ? [onStep] : []) };
  recordStructuralMetric('workerInitializationCount');
  currentWorker.postMessage(buildInitPayload(id));

  return promise;
};

export const getDiagramGenerationWorkerCapabilities = () => workerCapabilities;

/** @param {{ files?: TransferFile[], stagedResources?: TransferFile[] } | null | undefined} payload */
const collectTransferList = (payload) => {
  /** @type {Set<ArrayBuffer>} */
  const buffers = new Set();
  [
    ...(payload?.files || []),
    ...(payload?.stagedResources || [])
  ].forEach((file) => {
    if (file?.bytes instanceof ArrayBuffer) {
      buffers.add(file.bytes);
    }
  });
  return Array.from(buffers);
};

/**
 * @param {GenerateRequest | null | undefined} request
 * @param {() => void} callback
 */
const settleActiveRequest = (request, callback) => {
  if (!request || request.settled) return;
  request.settled = true;
  if (activeRequest === request) activeRequest = null;
  request.cleanup?.();
  callback();
};

/**
 * @param {GenerationPayload} [payload]
 * @param {{ onProgress?: ((progress: RequestProgress) => void) | null, held?: Readonly<HeldRenderState>, drawingId?: string | null }} [options]
 *   `held`: what the page holds for the drawing it renders (selectHeldRenderState); the
 *   catalog's revision and each segment's revision and byte count are posted.
 *   `drawingId`: the drawing it renders; the Worker keeps each drawing's prepared inputs
 *   while the retain set names their resources (A4b).
 * @returns {Promise<GenerationResponse>}
 */
export const runDiagramGeneration = (payload = {}, { onProgress = null, held = NOTHING_HELD, drawingId = null } = {}) => {
  if (activeRequest) {
    return Promise.reject(new Error('A diagram generation request is already running.'));
  }

  const requestId = nextRequestId;
  nextRequestId += 1;

  // The Promise executor runs synchronously and assigns both before anything reads them.
  /** @type {(value: GenerationResponse) => void} */
  let resolveRequest = () => {};
  /** @type {(error: unknown) => void} */
  let rejectRequest = () => {};
  /** @type {Promise<GenerationResponse>} */
  const promise = new Promise((resolve, reject) => {
    resolveRequest = resolve;
    rejectRequest = reject;
  });
  /** @type {GenerateRequest} */
  const request = {
    requestId,
    settled: false,
    cleanup: null,
    resolve: resolveRequest,
    reject: rejectRequest
  };
  activeRequest = request;
  /** @param {string} stage @param {RuntimePreparationStep | null} [step] */
  const reportProgress = (stage, step = null) => {
    if (!request.settled && activeRequest === request && typeof onProgress === 'function') {
      onProgress({ requestId, stage, ...(step || {}) });
    }
  };

  (async () => {
    let failureStage = 'initialization';
    try {
      // A drain is not runtime preparation, but its wait is not a warm start either.
      const initializationWasWarm = workerInitialized && !drainState;
      if (!workerInitialized) reportProgress('preparing-runtime');
      const initializationStartedAt = globalThis.performance?.now?.() ?? Date.now();
      if (drainState) {
        await awaitDrain();
        if (request.settled || activeRequest !== request) return;
        // The drain ended by termination: a new Worker boots.
        if (!workerInitialized) reportProgress('preparing-runtime');
      }
      const currentWorker = await ensureWorkerInitialized((step) => reportProgress('preparing-runtime', step));
      recordSessionLifecycleEvent('worker-initialization-resolved', {
        reused: initializationWasWarm,
        durationMs: initializationWasWarm
          ? 0
          : (globalThis.performance?.now?.() ?? Date.now()) - initializationStartedAt
      });
      if (request.settled || activeRequest !== request) return;
      reportProgress('preparing-resources');
      failureStage = 'resource-staging';
      const preparedResources = await resourceTransport.prepare(payload);
      if (request.settled || activeRequest !== request) return;
      const workerPayload = {
        request: payload.request,
        resourceManifest: preparedResources.resourceManifest,
        stagedResources: preparedResources.stagedResources,
        retainedResourceTokens: preparedResources.retainedResourceTokens,
        drawingId,
        fontIdentities: requestFontIdentities(payload.request),
        heldFeatureCatalogRevision: held.catalog?.revision ?? null,
        heldSegments: heldSegmentList(held)
      };

      const cleanup = () => {
        currentWorker.removeEventListener('message', handleMessage);
        currentWorker.removeEventListener('error', handleError);
        currentWorker.removeEventListener('messageerror', handleMessageError);
      };

      /** @param {unknown} error */
      const fail = (error) => {
        settleActiveRequest(request, () => rejectRequest(error));
      };

      /** @param {unknown} error */
      const failWorker = (error) => {
        fail(error);
        if (worker === currentWorker) terminateWorker(error);
      };

      /** @param {MessageEvent} event */
      async function handleMessage(event) {
        const data = event.data || {};
        if (data.type === 'test-lifecycle' && data.requestId === requestId) {
          recordSessionLifecycleEvent(data.event?.name, data.event);
          return;
        }
        if (data.type === 'progress' && data.requestId === requestId) {
          reportProgress(data.stage);
          return;
        }
        if (data.type !== 'run' || data.requestId !== requestId) return;
        if (data.ok) {
          // A completed run acknowledges Worker resource staging even when
          // Python returns a structured render error. Artifact admission is
          // independent: the live Worker already owns this manifest's bytes.
          preparedResources.commit();
          const beforeResponse = runtimeTestHooksEnabled()
            ? globalThis.__GBDRAW_TEST_HOOKS__?.beforeDiagramGenerationResponse
            : null;
          if (typeof beforeResponse === 'function') {
            try {
              await beforeResponse({ requestId });
            } catch (error) {
              failWorker(error);
              return;
            }
            if (request.settled || activeRequest !== request) return;
          }
          /** @type {NormalizedGeneration} */
          let normalized;
          try {
            normalized = normalizeGenerationResponse(data.results, held);
          } catch (error) {
            failWorker(error);
            return;
          }
          settleActiveRequest(request, () => {
            resolveRequest(markCurrentWorkerGenerationResponse({
              requestId,
              ...normalized
            }));
          });
          return;
        }
        failWorker(deserializeWorkerError(data.error, { operation: 'generate' }));
      }

      function handleError() {
        failWorker(deserializeWorkerError({ code: 'UNKNOWN', operation: 'generate', stage: 'unknown' }));
      }

      function handleMessageError() {
        failWorker(deserializeWorkerError({ code: 'RESULT_INVALID', operation: 'generate', stage: 'result-admission' }));
      }

      request.cleanup = cleanup;
      currentWorker.addEventListener('message', handleMessage);
      currentWorker.addEventListener('error', handleError);
      currentWorker.addEventListener('messageerror', handleMessageError);
      recordSessionLifecycleEvent('worker-post-start', { requestId });
      currentWorker.postMessage(
        {
          type: 'run',
          requestId,
          payload: workerPayload,
          testLifecycleEnabled: runtimeTestHooksEnabled()
        },
        collectTransferList(workerPayload)
      );
      recordSessionLifecycleEvent('worker-post-end', { requestId });
    } catch (error) {
      if (request.settled || activeRequest !== request) return;
      settleActiveRequest(request, () => rejectRequest(deserializeWorkerError(error, { operation: 'generate', stage: failureStage,
        code: failureStage === 'initialization' ? 'WORKER_INIT' : 'RESOURCE_INVALID' })));
    }
  })();

  return promise;
};

/** @param {AuxiliaryPayload} [payload] */
export const runFeatureExtraction = (payload = {}) => {
  return runAuxiliaryWorkerRequest({
    type: 'feature-extraction',
    payload,
    activeRequests: activeFeatureRequests,
    operation: 'feature-extraction'
  });
};

/**
 * @param {{
 *   type: string,
 *   payload: AuxiliaryPayload,
 *   operation?: string | null,
 *   activeRequests: Set<AuxiliaryRequest>,
 *   prepareResources?: (() => Promise<PreparedResources>) | null,
 *   onProgress?: ((progress: RequestProgress) => void) | null
 * }} options
 */
const runAuxiliaryWorkerRequest = ({
  type,
  payload,
  operation = null,
  activeRequests,
  prepareResources = null,
  onProgress = null
}) => {
  const requestId = nextRequestId;
  nextRequestId += 1;

  // The Promise executor runs synchronously and assigns both before anything reads them.
  /** @type {(value: AuxiliaryResponse) => void} */
  let resolveRequest = () => {};
  /** @type {(error: unknown) => void} */
  let rejectRequest = () => {};
  const promise = new Promise((resolve, reject) => {
    resolveRequest = resolve;
    rejectRequest = reject;
  });

  (async () => {
    /** @type {AuxiliaryRequest | null} */
    let request = null;
    let failureStage = 'initialization';
    try {
      const reportStep = typeof onProgress === 'function'
        ? (/** @type {RuntimePreparationStep | null} */ step) => onProgress({ requestId, stage: 'preparing-runtime', ...(step || {}) }) : null;
      if (drainState) await awaitDrain();
      if (!workerInitialized) reportStep?.(null);
      const currentWorker = await ensureWorkerInitialized(reportStep);
      request = {
        requestId,
        cleanup: null,
        reject: rejectRequest
      };

      const cleanup = () => {
        currentWorker.removeEventListener('message', handleMessage);
        currentWorker.removeEventListener('error', handleError);
        currentWorker.removeEventListener('messageerror', handleMessageError);
        // `request` is assigned just above this closure exists; nothing calls
        // `cleanup` before then.
        activeRequests.delete(/** @type {AuxiliaryRequest} */ (request));
      };
      request.cleanup = cleanup;
      activeRequests.add(request);
      /** @type {PreparedResources | null} */
      let preparedResources = null;
      let workerPayload = payload;

      /** @param {unknown} error */
      const fail = (error) => {
        cleanup();
        rejectRequest(error);
      };

      const assembly = createBoundedJsonReceiver();
      /** @param {MessageEvent} event */
      function handleMessage(event) {
        const data = event.data || {};
        if (data.type !== type || data.requestId !== requestId) return;
        if (data.status === 'part') {
          try {
            assembly.receivePart(data);
            currentWorker.postMessage({ type: 'auxiliary-ack', requestId });
          } catch { handleMessageError(); }
          return;
        }
        cleanup();
        if (data.ok) {
          try {
            if (Object.hasOwn(data, 'result')) {
              assembly.receivePart({ kind: 'value', path: [], value: data.result });
            }
            const result = assembly.getValue();
            preparedResources?.commit();
            resolveRequest({ requestId, result });
          } catch { handleMessageError(); }
          return;
        }
        rejectRequest(deserializeWorkerError(data.error, { operation }));
      }

      function handleError() {
        fail(deserializeWorkerError({ code: 'UNKNOWN', operation, stage: 'unknown' }));
      }

      function handleMessageError() {
        const error = deserializeWorkerError({ code: 'RESULT_INVALID', operation, stage: 'result-admission' });
        fail(error);
        if (activeRequest) {
          const pending = activeRequest;
          settleActiveRequest(pending, () => pending.reject(error));
        }
        if (worker === currentWorker) terminateWorker(error);
      }

      if (prepareResources) {
        failureStage = 'resource-staging';
        preparedResources = await prepareResources();
        if (!activeRequests.has(request)) return;
        workerPayload = { ...payload, resourceManifest: preparedResources.resourceManifest,
          stagedResources: preparedResources.stagedResources,
          retainedResourceTokens: preparedResources.retainedResourceTokens };
      }
      currentWorker.addEventListener('message', handleMessage);
      currentWorker.addEventListener('error', handleError);
      currentWorker.addEventListener('messageerror', handleMessageError);
      currentWorker.postMessage(
        {
          type,
          requestId,
          ...(operation ? { operation } : {}),
          payload: workerPayload
        },
        collectTransferList(workerPayload)
      );
    } catch (error) {
      request?.cleanup?.();
      rejectRequest(isDiagramGenerationCanceled(error) ? error : deserializeWorkerError(error, { operation, stage: failureStage,
        code: failureStage === 'initialization' ? 'WORKER_INIT' : 'RESOURCE_INVALID' }));
    }
  })();

  return promise;
};

/**
 * @param {unknown} operation
 * @param {AuxiliaryPayload & HelperPayload} [payload]
 * @param {{ onProgress?: ((progress: RequestProgress) => void) | null }} [options]
 */
export const runDiagramHelperOperation = (operation, payload = {}, { onProgress = null } = {}) => {
  const normalizedOperation = String(operation || '').trim();
  if (!diagramHelperOperationNames.has(normalizedOperation)) {
    return Promise.reject(
      deserializeWorkerError({ code: 'HELPER_PROTOCOL', stage: 'request-validation' })
    );
  }
  if (!payload || typeof payload !== 'object' || Array.isArray(payload)) {
    return Promise.reject(deserializeWorkerError({ code: 'HELPER_PROTOCOL', operation: normalizedOperation, stage: 'request-validation' }));
  }
  recordStructuralMetric('diagramHelperRequest', 1, { operation: normalizedOperation });
  const references = stagedHelperReferences(normalizedOperation, payload);
  const { resources, ...helperPayload } = payload;
  return runAuxiliaryWorkerRequest({
    type: 'helper', operation: normalizedOperation,
    payload: references ? helperPayload : payload,
    prepareResources: references ? () => resourceTransport.prepare({
      request: references, resources, retain: true
    }) : null,
    activeRequests: activeHelperRequests,
    onProgress
  });
};

// Cancel ends only Worker-side work. A Cancel during JS-side preparation keeps
// the warm Worker; the caller's cancel check stops the pending dispatch.
// A booted Worker with an interrupt buffer is kept and drained (D7); a Cancel
// during its boot, or on a page that is not cross-origin isolated, terminates it.
export const cancelDiagramGeneration = () => {
  const request = activeRequest;
  if (!request && !initState && !activeFeatureRequests.size && !activeHelperRequests.size) {
    return false;
  }
  const error = new DiagramGenerationCanceledError();
  if (request) {
    settleActiveRequest(request, () => request.reject(error));
  }
  const currentWorker = worker;
  if (!currentWorker || !workerInitialized || !workerInterruptBuffer) {
    terminateWorker(error);
    return true;
  }
  rejectWorkerRequests(error);
  rejectDrainWaiters(error);
  // The Worker pruned its resources for the canceled request: send every byte again.
  resourceTransport.reset();
  if (!drainState) drainWorker(currentWorker, workerInterruptBuffer);
  return true;
};

export const disposeDiagramGenerationWorker = () => {
  const error = new DiagramGenerationCanceledError('Diagram generation worker was disposed.');
  const request = activeRequest;
  if (request) {
    settleActiveRequest(request, () => request.reject(error));
  }
  terminateWorker(error);
};
