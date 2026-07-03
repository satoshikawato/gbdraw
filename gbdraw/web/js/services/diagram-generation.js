import { buildPyodideAssetManifest } from './pyodide-assets.js';

export class DiagramGenerationCanceledError extends Error {
  constructor(message = 'Diagram generation was canceled.') {
    super(message);
    this.name = 'DiagramGenerationCanceledError';
    this.canceled = true;
  }
}

export const isDiagramGenerationCanceled = (error) =>
  Boolean(error?.canceled || error?.name === 'DiagramGenerationCanceledError');

export const normalizeGenerationResponse = (payload) => {
  if (Array.isArray(payload)) {
    return { results: payload, metadata: {} };
  }
  if (!payload || typeof payload !== 'object') {
    return { results: [], metadata: {} };
  }
  if (payload.error) {
    return { results: payload, metadata: {} };
  }
  const results = Array.isArray(payload.results) ? payload.results : [];
  const metadata = (
    payload.metadata &&
    typeof payload.metadata === 'object' &&
    !Array.isArray(payload.metadata)
  )
    ? payload.metadata
    : {};
  return { results, metadata };
};

const resolveWorkerUrl = () => new URL('../workers/diagram-generation-worker.js', import.meta.url).toString();

let worker = null;
let workerInitialized = false;
let initState = null;
let activeRequest = null;
const activeFeatureRequests = new Set();
let nextRequestId = 1;

const buildInitPayload = (id) => {
  return {
    type: 'init',
    id,
    ...buildPyodideAssetManifest()
  };
};

const deserializeWorkerError = (serialized, fallbackMessage = 'Diagram generation worker failed') => {
  const message = serialized?.message ? String(serialized.message) : fallbackMessage;
  const error = new Error(message);
  error.name = serialized?.name ? String(serialized.name) : 'DiagramGenerationWorkerError';
  if (serialized?.stack) error.stack = String(serialized.stack);
  return error;
};

const getWorker = () => {
  if (!worker) {
    worker = new Worker(resolveWorkerUrl(), { type: 'module' });
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
  if (worker) {
    worker.terminate();
    worker = null;
  }
  workerInitialized = false;
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
  currentWorker.postMessage(buildInitPayload(id));

  return promise;
};

export const preinitializeDiagramGenerationWorker = () => ensureWorkerInitialized();

const collectTransferList = (payload) => {
  const buffers = new Set();
  (payload?.files || []).forEach((file) => {
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
      const currentWorker = await ensureWorkerInitialized();
      if (request.settled || activeRequest !== request) return;

      const cleanup = () => {
        currentWorker.removeEventListener('message', handleMessage);
        currentWorker.removeEventListener('error', handleError);
        currentWorker.removeEventListener('messageerror', handleMessageError);
      };

      const fail = (error) => {
        settleActiveRequest(request, () => rejectRequest(error));
      };

      function handleMessage(event) {
        const data = event.data || {};
        if (data.type !== 'run' || data.requestId !== requestId) return;
        if (data.ok) {
          settleActiveRequest(request, () => {
            resolveRequest({ requestId, ...normalizeGenerationResponse(data.results) });
          });
          return;
        }
        fail(deserializeWorkerError(data.error, 'Diagram generation failed'));
      }

      function handleError(event) {
        fail(new Error(event.message || 'Diagram generation worker error'));
      }

      function handleMessageError() {
        fail(new Error('Diagram generation worker message could not be decoded'));
      }

      request.cleanup = cleanup;
      currentWorker.addEventListener('message', handleMessage);
      currentWorker.addEventListener('error', handleError);
      currentWorker.addEventListener('messageerror', handleMessageError);
      currentWorker.postMessage(
        {
          type: 'run',
          requestId,
          payload
        },
        collectTransferList(payload)
      );
    } catch (error) {
      if (request.settled || activeRequest !== request) return;
      settleActiveRequest(request, () => rejectRequest(error));
    }
  })();

  return promise;
};

export const runFeatureExtraction = (payload = {}) => {
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
        activeFeatureRequests.delete(request);
      };
      request.cleanup = cleanup;
      activeFeatureRequests.add(request);

      const fail = (error) => {
        cleanup();
        rejectRequest(error);
      };

      function handleMessage(event) {
        const data = event.data || {};
        if (data.type !== 'feature-extraction' || data.requestId !== requestId) return;
        cleanup();
        if (data.ok) {
          resolveRequest({ requestId, result: data.result });
          return;
        }
        rejectRequest(deserializeWorkerError(data.error, 'Feature extraction failed'));
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
          type: 'feature-extraction',
          requestId,
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

export const cancelDiagramGeneration = () => {
  const error = new DiagramGenerationCanceledError();
  const request = activeRequest;
  const hadActiveRequest = Boolean(request);
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
