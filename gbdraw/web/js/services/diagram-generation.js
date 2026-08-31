import { buildPyodideAssetManifest } from './pyodide-assets.js';
import { normalizeUserFacingError } from './error-normalization.js';
import { validateWebRuntimeCapabilities } from './runtime-capabilities.js';
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

export { DIAGRAM_HELPER_OPERATIONS } from './diagram-worker-protocol.js';

export class DiagramGenerationCanceledError extends Error {
  constructor(message = 'Diagram generation was canceled.') {
    super(message);
    this.name = 'DiagramGenerationCanceledError';
    this.canceled = true;
  }
}

export const isDiagramGenerationCanceled = (error) =>
  Boolean(error?.canceled || error?.name === 'DiagramGenerationCanceledError');

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

const decodeGenerationResults = (results) => (
  Array.isArray(results)
    ? results.map(decodeGenerationResultContent)
    : results
);

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

export const normalizeGenerationResponse = (payload) => {
  if (Array.isArray(payload)) {
    return { results: decodeGenerationResults(payload), metadata: {}, artifactIdentity: null };
  }
  if (!payload || typeof payload !== 'object') {
    return { results: [], metadata: {}, artifactIdentity: null };
  }
  if (payload.error) {
    return { results: payload, metadata: {}, artifactIdentity: null };
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
    ? decodedMetadata
    : {};
  return {
    results,
    metadata,
    artifactIdentity: normalizeGeneratedArtifactIdentity(payload.artifactIdentity)
  };
};

const resolveWorkerUrl = () => new URL('../workers/diagram-generation-worker.js', import.meta.url).toString();

let worker = null;
let workerInitialized = false;
let workerCapabilities = null;
let initState = null;
let activeRequest = null;
const activeFeatureRequests = new Set();
const activeHelperRequests = new Set();
let nextRequestId = 1;
const resourceTransport = createDiagramResourceTransport();

const diagramHelperOperationNames = new Set(DIAGRAM_HELPER_OPERATION_NAMES);

const buildInitPayload = (id) => {
  return {
    type: 'init',
    id,
    ...buildPyodideAssetManifest()
  };
};

export const deserializeWorkerError = (
  serialized,
  fallbackMessage = 'Diagram generation worker failed'
) => {
  const message = serialized?.message ? String(serialized.message) : fallbackMessage;
  const error = new Error(message);
  error.name = serialized?.name ? String(serialized.name) : 'DiagramGenerationWorkerError';
  if (Array.isArray(serialized?.details)) {
    error.details = normalizeUserFacingError({
      message,
      details: serialized.details
    })?.details || [];
  }
  if (Array.isArray(serialized?.notes)) {
    error.notes = serialized.notes.map((note) => String(note));
  }
  if (serialized?.stack) error.stack = String(serialized.stack);
  return error;
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

const rejectPendingInit = (error) => {
  if (!initState) return;
  const { reject } = initState;
  clearInitState();
  reject(error);
};

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

const ensureWorkerInitialized = () => {
  const currentWorker = getWorker();
  if (workerInitialized) return Promise.resolve(currentWorker);
  if (initState?.promise) return initState.promise;

  const id = `diagram-init-${Date.now()}-${nextRequestId}`;
  let resolveInit;
  let rejectInit;
  const promise = new Promise((resolve, reject) => {
    resolveInit = resolve;
    rejectInit = reject;
  });

  const cleanup = () => {
    currentWorker.removeEventListener('message', handleMessage);
    currentWorker.removeEventListener('error', handleError);
    currentWorker.removeEventListener('messageerror', handleMessageError);
  };

  const fail = (error) => {
    cleanup();
    initState = null;
    if (worker === currentWorker) terminateWorker();
    rejectInit(error);
  };

  function handleMessage(event) {
    const data = event.data || {};
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
      resolveInit(currentWorker);
      return;
    }
    fail(deserializeWorkerError(data.error, 'Diagram generation worker initialization failed'));
  }

  function handleError(event) {
    fail(new Error(event.message || 'Diagram generation worker initialization error'));
  }

  function handleMessageError() {
    fail(new Error('Diagram generation worker initialization message could not be decoded'));
  }

  currentWorker.addEventListener('message', handleMessage);
  currentWorker.addEventListener('error', handleError);
  currentWorker.addEventListener('messageerror', handleMessageError);
  initState = { promise, reject: rejectInit, cleanup };
  recordStructuralMetric('workerInitializationCount');
  currentWorker.postMessage(buildInitPayload(id));

  return promise;
};

export const getDiagramGenerationWorkerCapabilities = () => workerCapabilities;

const collectTransferList = (payload) => {
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

const settleActiveRequest = (request, callback) => {
  if (!request || request.settled) return;
  request.settled = true;
  if (activeRequest === request) activeRequest = null;
  request.cleanup?.();
  callback();
};

export const runDiagramGeneration = (payload = {}) => {
  if (activeRequest) {
    return Promise.reject(new Error('A diagram generation request is already running.'));
  }

  const requestId = nextRequestId;
  nextRequestId += 1;

  let resolveRequest;
  let rejectRequest;
  const promise = new Promise((resolve, reject) => {
    resolveRequest = resolve;
    rejectRequest = reject;
  });
  const request = {
    requestId,
    settled: false,
    cleanup: null,
    resolve: resolveRequest,
    reject: rejectRequest
  };
  activeRequest = request;

  (async () => {
    try {
      const initializationWasWarm = workerInitialized;
      const initializationStartedAt = globalThis.performance?.now?.() ?? Date.now();
      const currentWorker = await ensureWorkerInitialized();
      recordSessionLifecycleEvent('worker-initialization-resolved', {
        reused: initializationWasWarm,
        durationMs: initializationWasWarm
          ? 0
          : (globalThis.performance?.now?.() ?? Date.now()) - initializationStartedAt
      });
      if (request.settled || activeRequest !== request) return;
      const preparedResources = await resourceTransport.prepare(payload);
      if (request.settled || activeRequest !== request) return;
      const workerPayload = {
        request: payload.request,
        resourceManifest: preparedResources.resourceManifest,
        stagedResources: preparedResources.stagedResources
      };

      const cleanup = () => {
        currentWorker.removeEventListener('message', handleMessage);
        currentWorker.removeEventListener('error', handleError);
        currentWorker.removeEventListener('messageerror', handleMessageError);
      };

      const fail = (error) => {
        settleActiveRequest(request, () => rejectRequest(error));
      };

      const failWorker = (error) => {
        fail(error);
        if (worker === currentWorker) terminateWorker(error);
      };

      async function handleMessage(event) {
        const data = event.data || {};
        if (data.type === 'test-lifecycle' && data.requestId === requestId) {
          recordSessionLifecycleEvent(data.event?.name, data.event);
          return;
        }
        if (data.type !== 'run' || data.requestId !== requestId) return;
        if (data.ok) {
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
          let normalized;
          try {
            normalized = normalizeGenerationResponse(data.results);
          } catch (error) {
            failWorker(error);
            return;
          }
          let resourcePromotionFinalized = false;
          const finalizeResourcePromotion = () => {
            if (resourcePromotionFinalized) return false;
            preparedResources.commit();
            resourcePromotionFinalized = true;
            return true;
          };
          settleActiveRequest(request, () => {
            resolveRequest({ requestId, ...normalized, finalizeResourcePromotion });
          });
          return;
        }
        failWorker(deserializeWorkerError(data.error, 'Diagram generation failed'));
      }

      function handleError(event) {
        failWorker(new Error(event.message || 'Diagram generation worker error'));
      }

      function handleMessageError() {
        failWorker(new Error('Diagram generation worker message could not be decoded'));
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
      settleActiveRequest(request, () => rejectRequest(error));
    }
  })();

  return promise;
};

export const runFeatureExtraction = (payload = {}) => {
  return runAuxiliaryWorkerRequest({
    type: 'feature-extraction',
    payload,
    activeRequests: activeFeatureRequests,
    fallbackMessage: 'Feature extraction failed'
  });
};

const runAuxiliaryWorkerRequest = ({
  type,
  payload,
  operation = null,
  activeRequests,
  fallbackMessage
}) => {
  const requestId = nextRequestId;
  nextRequestId += 1;

  let resolveRequest;
  let rejectRequest;
  const promise = new Promise((resolve, reject) => {
    resolveRequest = resolve;
    rejectRequest = reject;
  });

  (async () => {
    let request = null;
    try {
      const currentWorker = await ensureWorkerInitialized();
      request = {
        requestId,
        cleanup: null,
        reject: rejectRequest
      };

      const cleanup = () => {
        currentWorker.removeEventListener('message', handleMessage);
        currentWorker.removeEventListener('error', handleError);
        currentWorker.removeEventListener('messageerror', handleMessageError);
        activeRequests.delete(request);
      };
      request.cleanup = cleanup;
      activeRequests.add(request);

      const fail = (error) => {
        cleanup();
        rejectRequest(error);
      };

      function handleMessage(event) {
        const data = event.data || {};
        if (data.type !== type || data.requestId !== requestId) return;
        cleanup();
        if (data.ok) {
          resolveRequest({ requestId, result: data.result });
          return;
        }
        rejectRequest(deserializeWorkerError(data.error, fallbackMessage));
      }

      function handleError(event) {
        fail(new Error(event.message || 'Diagram generation worker error'));
      }

      function handleMessageError() {
        fail(new Error('Diagram generation worker message could not be decoded'));
      }

      currentWorker.addEventListener('message', handleMessage);
      currentWorker.addEventListener('error', handleError);
      currentWorker.addEventListener('messageerror', handleMessageError);
      currentWorker.postMessage(
        {
          type,
          requestId,
          ...(operation ? { operation } : {}),
          payload
        },
        collectTransferList(payload)
      );
    } catch (error) {
      request?.cleanup?.();
      rejectRequest(error);
    }
  })();

  return promise;
};

export const runDiagramHelperOperation = (operation, payload = {}) => {
  const normalizedOperation = String(operation || '').trim();
  if (!diagramHelperOperationNames.has(normalizedOperation)) {
    return Promise.reject(
      new TypeError(`Unsupported diagram helper operation '${normalizedOperation || '(blank)'}'.`)
    );
  }
  if (!payload || typeof payload !== 'object' || Array.isArray(payload)) {
    return Promise.reject(new TypeError('Diagram helper payload must be an object.'));
  }
  return runAuxiliaryWorkerRequest({
    type: 'helper',
    operation: normalizedOperation,
    payload,
    activeRequests: activeHelperRequests,
    fallbackMessage: `Diagram helper operation '${normalizedOperation}' failed`
  });
};

export const cancelDiagramGeneration = () => {
  const error = new DiagramGenerationCanceledError();
  const request = activeRequest;
  const hadActiveRequest = Boolean(
    request || initState || activeFeatureRequests.size || activeHelperRequests.size
  );
  if (request) {
    settleActiveRequest(request, () => request.reject(error));
  }
  terminateWorker(error);
  return hadActiveRequest;
};

export const disposeDiagramGenerationWorker = () => {
  const error = new DiagramGenerationCanceledError('Diagram generation worker was disposed.');
  const request = activeRequest;
  if (request) {
    settleActiveRequest(request, () => request.reject(error));
  }
  terminateWorker(error);
};
