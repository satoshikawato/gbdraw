// @ts-check
import { createBoundedJsonReceiver } from './bounded-json-transport.js';
import { buildPyodideAssetManifest } from './pyodide-assets.js';
import { normalizeUserFacingError } from '../utils/error-normalization.js';
import { validateWebRuntimeCapabilities, DiagramRuntimeCompatibilityError } from './runtime-capabilities.js';
import {
  DIAGRAM_HELPER_OPERATION_NAMES,
  DIAGRAM_HELPER_OPERATIONS
} from './diagram-worker-protocol.js';
import {
  recordSessionLifecycleEvent,
  recordStructuralMetric,
  runtimeTestHooksEnabled
} from './runtime-test-hooks.js';
import { createDiagramResourceTransport } from './diagram-resource-staging.js';
import { markCurrentWorkerGenerationResponse } from './current-worker-result-source.js';
import { normalizeLabelFitFailure } from './label-fit.js';
/** @import { LabelFitFailure } from './label-fit.js' */
import { webFontIdentitiesFor } from './font-store.js';

/** @import { UserFacingOptions } from '../utils/error-normalization.js' */
/** @import { ResourceDescriptor } from './diagram-resource-staging.js' */
/** @import { BoundedJsonValue } from './bounded-json-transport.js' */
/** @import { WebFontEntry } from './font-store.js' */

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
 * @property {unknown} [labelFit] The JSON text of a label fit report, next to an error.
 * @property {Partial<GeneratedArtifactIdentity> | null} [artifactIdentity]
 */

/**
 * A Result entry as the Worker sends it: its `content` is text or UTF-8 bytes.
 * @typedef {{ content?: unknown, [key: string]: unknown }} WorkerResultEntry
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

/**
 * @param {WorkerResultEntry | string | number | boolean | null | undefined} result
 * @param {number} resultIndex
 */
const decodeGenerationResultContent = (result, resultIndex) => {
  if (!result || typeof result !== 'object' || typeof result.content === 'string') {
    return result;
  }
  const bytes = result.content instanceof ArrayBuffer
    ? new Uint8Array(result.content)
    : ArrayBuffer.isView(result.content)
      ? new Uint8Array(
          result.content.buffer,
          result.content.byteOffset,
          result.content.byteLength
        )
      : null;
  if (!bytes) return result;
  recordSessionLifecycleEvent('result-binary-decode-start', {
    resultIndex,
    bytes: bytes.byteLength
  });
  const content = new TextDecoder().decode(bytes);
  recordStructuralMetric('resultBinaryDecodeCount', 1, { resultIndex });
  recordStructuralMetric('resultBinaryDecodedBytes', bytes.byteLength, { resultIndex });
  recordSessionLifecycleEvent('result-binary-decode-end', {
    resultIndex,
    bytes: bytes.byteLength,
    characters: content.length
  });
  return { ...result, content };
};

/** @param {unknown} results */
const decodeGenerationResults = (results) => (
  Array.isArray(results)
    ? results.map(decodeGenerationResultContent)
    : results
);

/**
 * @param {unknown} metadata The Worker's metadata: an object, or its UTF-8 JSON bytes.
 * @returns {unknown}
 */
const decodeGenerationMetadata = (metadata) => {
  const bytes = metadata instanceof ArrayBuffer
    ? new Uint8Array(metadata)
    : ArrayBuffer.isView(metadata)
      ? new Uint8Array(metadata.buffer, metadata.byteOffset, metadata.byteLength)
      : null;
  if (!bytes) return metadata;
  recordSessionLifecycleEvent('result-metadata-binary-decode-start', {
    bytes: bytes.byteLength
  });
  const decoded = JSON.parse(new TextDecoder().decode(bytes));
  recordStructuralMetric('resultMetadataBinaryDecodeCount');
  recordStructuralMetric('resultMetadataBinaryDecodedBytes', bytes.byteLength);
  recordSessionLifecycleEvent('result-metadata-binary-decode-end', {
    bytes: bytes.byteLength
  });
  return decoded;
};

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
 * @returns {NormalizedGeneration}
 */
export const normalizeGenerationResponse = (payload) => {
  if (Array.isArray(payload)) {
    return markCurrentWorkerGenerationResponse({
      results: decodeGenerationResults(payload),
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
    ? decodeGenerationResults(payload.results)
    : [];
  const decodedMetadata = decodeGenerationMetadata(payload.metadata);
  const metadata = (
    decodedMetadata &&
    typeof decodedMetadata === 'object' &&
    !Array.isArray(decodedMetadata)
  )
    // A non-array object, by the condition.
    ? /** @type {Record<string, unknown>} */ (decodedMetadata)
    : {};
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
 * @property {() => boolean} commit Acknowledges the staging; false when a reset
 *   or another acknowledged request superseded it.
 */
/** @type {Set<AuxiliaryRequest>} */
const activeFeatureRequests = new Set();
/** @type {Set<AuxiliaryRequest>} */
const activeHelperRequests = new Set();
let nextRequestId = 1;
const resourceTransport = createDiagramResourceTransport();

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
  const options = /** @type {{ diagramOptions?: { configOverrides?: Record<string, unknown>, config?: { objects?: { text?: { font_family?: unknown } } } } } | null} */ (request)?.diagramOptions;
  const fontFamily = options?.configOverrides?.['objects.text.font_family'] ?? options?.config?.objects?.text?.font_family;
  return fontFamily ? webFontIdentitiesFor(workerFontSource.readFonts(), fontFamily) : [];
};

/** @param {string} id */
const buildInitPayload = (id) => {
  return {
    type: 'init',
    id,
    fonts: workerFontFiles(workerFontSource.readFonts()),
    ...buildPyodideAssetManifest()
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

/** @param {unknown} [error] */
const terminateWorker = (error = null) => {
  if (error) rejectPendingInit(error);
  else clearInitState();
  activeFeatureRequests.forEach((request) => {
    request.cleanup?.();
    request.reject(error || new Error('Diagram generation worker was terminated.'));
  });
  activeFeatureRequests.clear();
  activeHelperRequests.forEach((request) => {
    request.cleanup?.();
    request.reject(error || new Error('Diagram generation worker was terminated.'));
  });
  activeHelperRequests.clear();
  if (worker) {
    worker.terminate();
    worker = null;
  }
  workerInitialized = false;
  workerCapabilities = null;
  resourceTransport.reset();
};

/**
 * @param {((step: RuntimePreparationStep) => void) | null} [onStep]
 *   Each first-use preparation step as the Worker starts it; a caller that
 *   joins a running initialization first receives the current step.
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
    cleanup();
    initState = null;
    if (data.ok) {
      try {
        workerCapabilities = validateWebRuntimeCapabilities(data.capabilities);
      } catch (error) {
        fail(error);
        return;
      }
      workerInitialized = true;
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
 * @param {{ onProgress?: ((progress: RequestProgress) => void) | null }} [options]
 * @returns {Promise<GenerationResponse>}
 */
export const runDiagramGeneration = (payload = {}, { onProgress = null } = {}) => {
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
      const initializationWasWarm = workerInitialized;
      if (!initializationWasWarm) reportProgress('preparing-runtime');
      const initializationStartedAt = globalThis.performance?.now?.() ?? Date.now();
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
        fontIdentities: requestFontIdentities(payload.request)
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
            normalized = normalizeGenerationResponse(data.results);
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
          stagedResources: preparedResources.stagedResources };
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
  // Helpers that read the records of a committed request stage its resources.
  const staging = normalizedOperation === DIAGRAM_HELPER_OPERATIONS.READ_FEATURE_OVERRIDE_TABLE ? payload
    : (normalizedOperation === DIAGRAM_HELPER_OPERATIONS.RESOLVE_SIMILARITY_ALIGNMENT ? payload.projection : undefined);
  const { resources, ...helperPayload } = payload;
  return runAuxiliaryWorkerRequest({
    type: 'helper', operation: normalizedOperation,
    payload: staging ? helperPayload : payload,
    prepareResources: staging ? () => resourceTransport.prepare({
      request: staging.canonicalRequest, resources
    }) : null,
    activeRequests: activeHelperRequests,
    onProgress
  });
};

// Cancel ends only Worker-side work. A Cancel during JS-side preparation keeps
// the warm Worker; the caller's cancel check stops the pending dispatch.
export const cancelDiagramGeneration = () => {
  const request = activeRequest;
  if (!request && !initState && !activeFeatureRequests.size && !activeHelperRequests.size) {
    return false;
  }
  const error = new DiagramGenerationCanceledError();
  if (request) {
    settleActiveRequest(request, () => request.reject(error));
  }
  terminateWorker(error);
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
