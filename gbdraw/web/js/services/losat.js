import { LOSAT_THREADED_WASM_URL, WASI_SHIM_URL } from '../config.js';

const DEFAULT_WASM_PATH = './wasm/losat/losat.wasm';
const DEFAULT_THREADED_WASM_PATH = LOSAT_THREADED_WASM_URL || './wasm/losat/losat-threaded.wasm';
const DEFAULT_MAX_WORKERS = 4;
const DEFAULT_THREADED_MIN_FASTA_CHARS = 500000;
const DEFAULT_MAX_TOTAL_THREADS = 16;
const DEFAULT_MAX_THREADS_PER_JOB = 16;
const SUPPORTED_PROGRAMS = new Set(['blastn', 'tblastx', 'blastp']);

let wasiShimPromise = null;
const wasmModulePromises = new Map();
let directInstancePromise = null;
const threadedSupportPromises = new Map();

const resolveAssetUrl = (path) => new URL(path, window.location.href).toString();
const resolveWorkerUrl = () => new URL('../workers/losat-worker.js', import.meta.url).toString();
const resolveThreadedWorkerUrl = () => new URL('../workers/losat-threaded-worker.js', import.meta.url).toString();
const getNow = () => (globalThis.performance?.now ? performance.now() : Date.now());
const formatDuration = (startedAt) => `${((getNow() - startedAt) / 1000).toFixed(2)}s`;

const createAbortError = () => {
  const error = new Error('LOSAT run was canceled.');
  error.name = 'AbortError';
  error.canceled = true;
  return error;
};

const getAbortReason = (signal) => {
  const reason = signal?.reason;
  if (reason instanceof Error) return reason;
  if (reason !== undefined) {
    const error = new Error(String(reason));
    error.name = 'AbortError';
    error.canceled = true;
    return error;
  }
  return createAbortError();
};

const throwIfAborted = (signal) => {
  if (signal?.aborted) throw getAbortReason(signal);
};

const isAbortError = (error, signal) =>
  Boolean(signal?.aborted && (error === signal.reason || error?.canceled || error?.name === 'AbortError'));

const loadWasiShim = async () => {
  if (!wasiShimPromise) {
    const wasiShimUrl = resolveAssetUrl(WASI_SHIM_URL);
    wasiShimPromise = import(wasiShimUrl).catch((error) => {
      throw new Error(
        `Missing packaged asset: browser_wasi_shim module at ${wasiShimUrl}. ${error.message}`
      );
    });
  }
  return wasiShimPromise;
};

const loadLosatModule = async (wasmPath = DEFAULT_WASM_PATH) => {
  const key = String(wasmPath || DEFAULT_WASM_PATH);
  if (!wasmModulePromises.has(key)) {
    wasmModulePromises.set(key, (async () => {
      const resolvedWasmPath = resolveAssetUrl(wasmPath);
      const response = await fetch(resolvedWasmPath, { cache: 'no-store' });
      if (!response.ok) {
        throw new Error(`Missing packaged asset: LOSAT wasm not found at ${resolvedWasmPath}`);
      }
      const bytes = await response.arrayBuffer();
      return WebAssembly.compile(bytes);
    })());
  }
  return wasmModulePromises.get(key);
};

const concatUint8Arrays = (chunks) => {
  const total = chunks.reduce((sum, chunk) => sum + chunk.length, 0);
  const merged = new Uint8Array(total);
  let offset = 0;
  chunks.forEach((chunk) => {
    merged.set(chunk, offset);
    offset += chunk.length;
  });
  return merged;
};

const hasDirectLosatApi = (wasmModule) =>
  WebAssembly.Module.exports(wasmModule).some(
    (entry) => entry.kind === 'function' && entry.name === 'losat_web_run_pair'
  );

const moduleHasThreadedWasiShape = (wasmModule) => {
  const exports = WebAssembly.Module.exports(wasmModule);
  const imports = WebAssembly.Module.imports(wasmModule);
  return exports.some((entry) => entry.kind === 'function' && entry.name === '_start') &&
    exports.some((entry) => entry.kind === 'function' && entry.name === 'wasi_thread_start') &&
    imports.some((entry) => entry.module === 'env' && entry.name === 'memory' && entry.kind === 'memory') &&
    imports.some((entry) => entry.module === 'wasi' && entry.name === 'thread-spawn' && entry.kind === 'function');
};

const buildThreadingStatus = (state, message, details = {}) => ({
  state,
  message,
  ...details
});

export const getLosatThreadingSupport = async ({
  threadedWasmPath = DEFAULT_THREADED_WASM_PATH
} = {}) => {
  if (typeof Worker !== 'function') {
    return buildThreadingStatus('unavailable', 'Web Workers are unavailable in this browser.');
  }
  if (typeof SharedArrayBuffer !== 'function') {
    return buildThreadingStatus('unavailable', 'SharedArrayBuffer is unavailable.');
  }
  if (globalThis.crossOriginIsolated !== true) {
    return buildThreadingStatus('unavailable', 'Cross-origin isolation is not enabled.');
  }

  const key = String(threadedWasmPath || DEFAULT_THREADED_WASM_PATH);
  if (!threadedSupportPromises.has(key)) {
    threadedSupportPromises.set(key, (async () => {
      const wasmModule = await loadLosatModule(threadedWasmPath);
      if (!moduleHasThreadedWasiShape(wasmModule)) {
        return buildThreadingStatus(
          'unavailable',
          'Threaded LOSAT wasm is missing WASI thread imports or exports.',
          { wasmModule: null }
        );
      }
      return buildThreadingStatus(
        'available',
        'Threaded LOSAT is available.',
        { wasmModule }
      );
    })().catch((error) =>
      buildThreadingStatus(
        'unavailable',
        error?.message ? String(error.message) : String(error || 'Threaded LOSAT is unavailable.')
      )
    ));
  }
  return threadedSupportPromises.get(key);
};

const instantiateDirectLosat = async (
  { WASI, File, OpenFile, PreopenDirectory, ConsoleStdout },
  wasmModule
) => {
  if (!directInstancePromise) {
    directInstancePromise = (async () => {
      const stdout = new ConsoleStdout(() => {});
      const stderr = new ConsoleStdout(() => {});
      const stdin = new OpenFile(new File(new Uint8Array(), { readonly: true }));
      const preopen = new PreopenDirectory('.', new Map());
      const wasi = new WASI([], [], [stdin, stdout, stderr, preopen]);
      const instance = await WebAssembly.instantiate(wasmModule, {
        wasi_snapshot_preview1: wasi.wasiImport
      });
      if (typeof wasi.initialize === 'function') {
        wasi.initialize(instance);
      }
      return instance;
    })();
  }
  return directInstancePromise;
};

const copyBytesToWasm = (exports, bytes) => {
  if (bytes.length === 0) return { ptr: 0, len: 0 };
  const ptr = exports.losat_web_alloc(bytes.length);
  new Uint8Array(exports.memory.buffer).set(bytes, ptr);
  return { ptr, len: bytes.length };
};

const readWasmString = (exports, ptrFn, lenFn) => {
  const ptr = ptrFn();
  const len = lenFn();
  if (!ptr || !len) return '';
  const bytes = new Uint8Array(exports.memory.buffer, ptr, len).slice();
  return new TextDecoder().decode(bytes);
};

const runLosatPairDirect = async ({
  program,
  queryFasta,
  subjectFasta,
  outfmt = '6',
  extraArgs = [],
  wasiShim,
  wasmModule
}) => {
  const instance = await instantiateDirectLosat(wasiShim, wasmModule);
  const exports = instance.exports;
  const encoder = new TextEncoder();
  const argsText = Array.isArray(extraArgs) ? extraArgs.map(String).join('\0') : '';
  const allocations = [
    copyBytesToWasm(exports, encoder.encode(String(program))),
    copyBytesToWasm(exports, encoder.encode(String(queryFasta))),
    copyBytesToWasm(exports, encoder.encode(String(subjectFasta))),
    copyBytesToWasm(exports, encoder.encode(String(outfmt))),
    copyBytesToWasm(exports, encoder.encode(argsText))
  ];

  try {
    const exitCode = exports.losat_web_run_pair(
      allocations[0].ptr,
      allocations[0].len,
      allocations[1].ptr,
      allocations[1].len,
      allocations[2].ptr,
      allocations[2].len,
      allocations[3].ptr,
      allocations[3].len,
      allocations[4].ptr,
      allocations[4].len
    );
    if (exitCode !== 0) {
      const detail = readWasmString(
        exports,
        exports.losat_web_error_ptr,
        exports.losat_web_error_len
      );
      throw new Error(detail || `LOSAT exited with code ${exitCode}`);
    }
    return readWasmString(exports, exports.losat_web_result_ptr, exports.losat_web_result_len);
  } finally {
    allocations.forEach(({ ptr, len }) => {
      if (ptr && len) exports.losat_web_dealloc(ptr, len);
    });
    if (typeof exports.losat_web_clear === 'function') {
      exports.losat_web_clear();
    }
  }
};

export const runLosatPair = async ({
  program,
  queryFasta,
  subjectFasta,
  outfmt = '6',
  extraArgs = [],
  wasmPath = DEFAULT_WASM_PATH,
  signal
} = {}) => {
  throwIfAborted(signal);
  if (!program || !SUPPORTED_PROGRAMS.has(program)) {
    throw new Error('LOSAT program must be blastn, tblastx, or blastp.');
  }
  if (!queryFasta || !subjectFasta) {
    throw new Error('LOSAT requires both query and subject FASTA content.');
  }

  const [wasiShim, wasmModule] = await Promise.all([
    loadWasiShim(),
    loadLosatModule(wasmPath)
  ]);
  throwIfAborted(signal);
  const { WASI, File, OpenFile, PreopenDirectory, ConsoleStdout } = wasiShim;

  if (hasDirectLosatApi(wasmModule)) {
    const text = await runLosatPairDirect({
      program,
      queryFasta,
      subjectFasta,
      outfmt,
      extraArgs,
      wasiShim,
      wasmModule
    });
    throwIfAborted(signal);
    return text;
  }

  const encoder = new TextEncoder();
  const files = new Map([
    ['query.fa', new File(encoder.encode(queryFasta), { readonly: true })],
    ['subject.fa', new File(encoder.encode(subjectFasta), { readonly: true })]
  ]);
  const preopen = new PreopenDirectory('.', files);

  const stdoutChunks = [];
  const stderrChunks = [];
  const stdout = new ConsoleStdout((data) => stdoutChunks.push(data));
  const stderr = new ConsoleStdout((data) => stderrChunks.push(data));
  const stdin = new OpenFile(new File(new Uint8Array(), { readonly: true }));
  const fds = [stdin, stdout, stderr, preopen];

  const extraList = Array.isArray(extraArgs) ? extraArgs : [];
  const args = [
    'losat',
    program,
    '--query',
    'query.fa',
    '--subject',
    'subject.fa',
    '--outfmt',
    String(outfmt),
    ...extraList
  ];

  const wasi = new WASI(args, [], fds);
  const instance = await WebAssembly.instantiate(wasmModule, {
    wasi_snapshot_preview1: wasi.wasiImport
  });
  throwIfAborted(signal);
  const exitCode = wasi.start(instance);
  throwIfAborted(signal);

  const stdoutText = new TextDecoder().decode(concatUint8Arrays(stdoutChunks));
  const stderrText = new TextDecoder().decode(concatUint8Arrays(stderrChunks));

  if (exitCode !== 0) {
    const detail = stderrText.trim() || `LOSAT exited with code ${exitCode}`;
    throw new Error(detail);
  }

  return stdoutText;
};

const getDefaultConcurrency = (jobCount) => {
  if (!Number.isFinite(jobCount) || jobCount <= 0) return 0;
  const hardwareLimit = Math.max(1, (globalThis.navigator?.hardwareConcurrency || 4) - 1);
  return Math.min(jobCount, hardwareLimit, DEFAULT_MAX_WORKERS);
};

const getHardwareThreadBudget = () =>
  Math.max(1, Number(globalThis.navigator?.hardwareConcurrency || 4) || 4);

const normalizeTotalThreadBudget = (value) => {
  const hardwareBudget = getHardwareThreadBudget();
  const normalized = String(value ?? 'safe').trim().toLowerCase();
  if (!normalized || normalized === 'safe' || normalized === 'auto') {
    return Math.min(DEFAULT_MAX_TOTAL_THREADS, hardwareBudget);
  }
  if (normalized === 'available') return hardwareBudget;
  const parsed = Number(value);
  if (!Number.isInteger(parsed) || parsed < 1) {
    return Math.min(DEFAULT_MAX_TOTAL_THREADS, hardwareBudget);
  }
  return Math.max(1, Math.min(parsed, hardwareBudget));
};

const normalizeExecutionMode = (value) => {
  const mode = String(value || 'auto').trim().toLowerCase();
  return ['auto', 'serial', 'threaded'].includes(mode) ? mode : 'auto';
};

const normalizeThreadsPerJob = (value, { fallback = null, maxThreads = DEFAULT_MAX_THREADS_PER_JOB } = {}) => {
  if (value === undefined || value === null || String(value).trim().toLowerCase() === 'auto') {
    return fallback;
  }
  const parsed = Number(value);
  if (!Number.isInteger(parsed) || parsed < 1) return fallback;
  return Math.min(parsed, Math.max(1, Number(maxThreads) || DEFAULT_MAX_THREADS_PER_JOB));
};

const getAutoThreadsPerJob = (jobs) => {
  const budget = Math.min(DEFAULT_MAX_THREADS_PER_JOB, getHardwareThreadBudget());
  if (!Array.isArray(jobs) || jobs.length !== 1) return Math.max(1, Math.min(2, budget));
  return Math.max(2, Math.min(budget, getHardwareThreadBudget() - 1 || 1));
};

const getThreadedFastaCharCount = (job, sequenceStore) => {
  let total = 0;
  try {
    total += resolveJobSequence(sequenceStore, job.querySequenceKey, 'query').length;
    total += resolveJobSequence(sequenceStore, job.subjectSequenceKey, 'subject').length;
  } catch {
    return 0;
  }
  return total;
};

const shouldUseThreadedLosat = (jobs, sequenceStore, threadsPerJob) => {
  if (!Array.isArray(jobs) || jobs.length === 0 || threadsPerJob <= 1) return false;
  if (jobs.length === 1) {
    return getThreadedFastaCharCount(jobs[0], sequenceStore) >= DEFAULT_THREADED_MIN_FASTA_CHARS;
  }
  const totalChars = jobs.reduce((sum, job) => sum + getThreadedFastaCharCount(job, sequenceStore), 0);
  return totalChars >= DEFAULT_THREADED_MIN_FASTA_CHARS * 2;
};

const buildRuntimeStatus = (state, message, details = {}) => ({
  state,
  message,
  ...details
});

const notifyRuntimeStatus = (callback, status) => {
  if (typeof callback === 'function') callback(status);
};

const buildThreadedRuntimePlan = (jobs, options, sequenceStore) => {
  const hasExplicitTotalBudget = Number.isFinite(options.totalThreadBudget);
  const requestedPairWorkers = Number.isFinite(options.concurrency)
    ? Math.max(1, Math.floor(options.concurrency))
    : hasExplicitTotalBudget
      ? jobs.length
      : getDefaultConcurrency(jobs.length);
  const autoThreadsPerJob = getAutoThreadsPerJob(jobs);
  const totalBudget = normalizeTotalThreadBudget(options.totalThreadBudget);
  const requestedThreadsPerJob = normalizeThreadsPerJob(options.threadsPerJob, {
    fallback: autoThreadsPerJob,
    maxThreads: Math.max(1, totalBudget - 1)
  });
  const threadsPerJob = Math.max(1, Math.min(requestedThreadsPerJob, Math.max(1, totalBudget - 1)));
  const workersPerThreadedJob = Math.max(1, threadsPerJob + 1);
  const pairWorkers = Math.max(
    1,
    Math.min(jobs.length, requestedPairWorkers, Math.max(1, Math.floor(totalBudget / workersPerThreadedJob)))
  );
  return {
    pairWorkers,
    threadsPerJob,
    useful: shouldUseThreadedLosat(jobs, sequenceStore, threadsPerJob),
    totalBudget
  };
};

const formatPairErrorPrefix = (job) => {
  const pairNumber = Number.isInteger(job?.pairIndex) ? job.pairIndex + 1 : null;
  return pairNumber ? `LOSAT pair #${pairNumber}` : 'LOSAT pair';
};

const normalizeSequenceStore = (sequences) => {
  if (sequences instanceof Map) return sequences;
  if (Array.isArray(sequences)) {
    return new Map(
      sequences
        .filter((entry) => entry && entry.key !== undefined)
        .map((entry) => [String(entry.key), String(entry.fasta || '')])
    );
  }
  if (sequences && typeof sequences === 'object') {
    return new Map(Object.entries(sequences).map(([key, fasta]) => [String(key), String(fasta || '')]));
  }
  return new Map();
};

const resolveJobSequence = (sequenceStore, key, label) => {
  const normalizedKey = String(key || '');
  if (!normalizedKey || !sequenceStore.has(normalizedKey)) {
    throw new Error(`Missing LOSAT ${label} sequence for key '${normalizedKey || '(blank)'}'.`);
  }
  const fasta = sequenceStore.get(normalizedKey);
  if (!fasta) {
    throw new Error(`LOSAT ${label} sequence for key '${normalizedKey}' is empty.`);
  }
  return fasta;
};

const materializeJobSequences = (job, sequenceStore) => ({
  ...job,
  queryFasta: resolveJobSequence(sequenceStore, job.querySequenceKey, 'query'),
  subjectFasta: resolveJobSequence(sequenceStore, job.subjectSequenceKey, 'subject')
});

const buildJobSequencePayload = (job, sequenceStore, loadedKeys) => {
  const payload = [];
  const addSequence = (key, label) => {
    const normalizedKey = String(key || '');
    if (!normalizedKey || loadedKeys.has(normalizedKey)) return;
    payload.push({
      key: normalizedKey,
      fasta: resolveJobSequence(sequenceStore, normalizedKey, label)
    });
    loadedKeys.add(normalizedKey);
  };
  addSequence(job.querySequenceKey, 'query');
  addSequence(job.subjectSequenceKey, 'subject');
  return payload;
};

const runLosatPairsSequential = async (jobs, { onProgress, sequences, signal } = {}) => {
  throwIfAborted(signal);
  const startedAt = getNow();
  const sequenceStore = normalizeSequenceStore(sequences);
  console.info(`Running ${jobs.length} LOSAT pair${jobs.length === 1 ? '' : 's'} sequentially.`);
  const results = [];
  for (const job of jobs) {
    throwIfAborted(signal);
    try {
      const text = await runLosatPair({
        ...materializeJobSequences(job, sequenceStore),
        signal
      });
      throwIfAborted(signal);
      results.push({ ...job, text });
      if (typeof onProgress === 'function') {
        onProgress({ completed: results.length, total: jobs.length, job });
      }
    } catch (error) {
      if (isAbortError(error, signal)) throw getAbortReason(signal);
      const message = error?.message ? String(error.message) : String(error || 'Unknown LOSAT error');
      throw new Error(`${formatPairErrorPrefix(job)}: ${message}`);
    }
  }
  console.info(`Completed ${jobs.length} LOSAT pair${jobs.length === 1 ? '' : 's'} sequentially in ${formatDuration(startedAt)}.`);
  return results;
};

const initializeLosatWorker = (worker, payload, signal) =>
  new Promise((resolve, reject) => {
    let settled = false;
    const settle = (callback) => {
      if (settled) return;
      settled = true;
      cleanup();
      callback();
    };
    const handleAbort = () => {
      settle(() => reject(getAbortReason(signal)));
    };
    const cleanup = () => {
      worker.removeEventListener('message', handleMessage);
      worker.removeEventListener('error', handleError);
      worker.removeEventListener('messageerror', handleMessageError);
      signal?.removeEventListener?.('abort', handleAbort);
    };
    const handleMessage = (event) => {
      const data = event.data || {};
      if (data.id !== payload.id || data.type !== 'init') return;
      if (data.ok) {
        settle(() => resolve(worker));
        return;
      }
      settle(() => reject(new Error(data.error || 'LOSAT worker initialization failed')));
    };

    const handleError = (event) => {
      settle(() => reject(new Error(event.message || 'LOSAT worker initialization error')));
    };

    const handleMessageError = () => {
      settle(() => reject(new Error('LOSAT worker initialization message could not be decoded')));
    };

    if (signal?.aborted) {
      handleAbort();
      return;
    }
    signal?.addEventListener?.('abort', handleAbort, { once: true });
    worker.addEventListener('message', handleMessage);
    worker.addEventListener('error', handleError);
    worker.addEventListener('messageerror', handleMessageError);
    worker.postMessage(payload);
  });

const runLosatPairsWithWorkers = async (
  jobs,
  { concurrency, workerUrl, wasmPath, onProgress, sequences, signal } = {}
) => {
  throwIfAborted(signal);
  const sequenceStore = normalizeSequenceStore(sequences);
  const workerCount = Math.min(
    jobs.length,
    Math.max(1, Number.isFinite(concurrency) ? Math.floor(concurrency) : getDefaultConcurrency(jobs.length))
  );
  if (workerCount <= 0) return [];

  throwIfAborted(signal);
  const startedAt = getNow();
  const resolvedWorkerUrl = workerUrl || resolveWorkerUrl();
  const resolvedWasmUrl = resolveAssetUrl(wasmPath || DEFAULT_WASM_PATH);
  const resolvedWasiShimUrl = resolveAssetUrl(WASI_SHIM_URL);
  const wasmModule = await loadLosatModule(wasmPath || DEFAULT_WASM_PATH);
  throwIfAborted(signal);
  const results = new Array(jobs.length);
  const workers = [];
  const workerLoadedSequenceKeys = new Map();
  let nextJobIndex = 0;
  let completed = 0;
  let settled = false;
  let requestId = 0;

  try {
    for (let i = 0; i < workerCount; i += 1) {
      const worker = new Worker(resolvedWorkerUrl, { type: 'module' });
      workers.push(worker);
      workerLoadedSequenceKeys.set(worker, new Set());
    }
    await Promise.all(
      workers.map((worker, index) =>
        initializeLosatWorker(
          worker,
          {
            type: 'init',
            id: `init-${index}-${Date.now()}`,
            wasmModule,
            wasmUrl: resolvedWasmUrl,
            wasiShimUrl: resolvedWasiShimUrl
          },
          signal
        )
      )
    );
    throwIfAborted(signal);
  } catch (error) {
    workers.forEach((worker) => worker.terminate());
    throw error;
  }

  return new Promise((resolve, reject) => {
    let handleAbort = null;
    const cleanup = () => {
      if (handleAbort) signal?.removeEventListener?.('abort', handleAbort);
      workers.forEach((worker) => worker.terminate());
    };

    const fail = (error) => {
      if (settled) return;
      settled = true;
      cleanup();
      reject(error);
    };

    handleAbort = () => fail(getAbortReason(signal));
    if (signal?.aborted) {
      handleAbort();
      return;
    }
    signal?.addEventListener?.('abort', handleAbort, { once: true });

    const maybeResolve = () => {
      if (settled || completed < jobs.length) return;
      settled = true;
      cleanup();
      console.info(
        `Completed ${jobs.length} LOSAT pair${jobs.length === 1 ? '' : 's'} with ${workerCount} worker${workerCount === 1 ? '' : 's'} in ${formatDuration(startedAt)}.`
      );
      resolve(results);
    };

    const assignNext = (worker) => {
      if (settled) return;
      if (signal?.aborted) {
        fail(getAbortReason(signal));
        return;
      }
      if (nextJobIndex >= jobs.length) {
        maybeResolve();
        return;
      }
      const index = nextJobIndex;
      nextJobIndex += 1;
      const job = jobs[index];
      const id = `${Date.now()}-${requestId}`;
      requestId += 1;
      const loadedKeys = workerLoadedSequenceKeys.get(worker) || new Set();
      let sequencePayload = [];
      try {
        sequencePayload = buildJobSequencePayload(job, sequenceStore, loadedKeys);
        workerLoadedSequenceKeys.set(worker, loadedKeys);
      } catch (error) {
        fail(new Error(`${formatPairErrorPrefix(job)}: ${error?.message || error}`));
        return;
      }

      const handleMessage = (event) => {
        const data = event.data || {};
        if (data.id !== id) return;
        if (data.type && data.type !== 'run') return;
        worker.removeEventListener('message', handleMessage);
        worker.removeEventListener('error', handleError);
        worker.removeEventListener('messageerror', handleMessageError);

        if (!data.ok) {
          fail(new Error(`${formatPairErrorPrefix(job)}: ${data.error || 'LOSAT worker failed'}`));
          return;
        }
        if (signal?.aborted) {
          fail(getAbortReason(signal));
          return;
        }

        results[index] = { ...job, text: data.text || '' };
        completed += 1;
        if (typeof onProgress === 'function') {
          onProgress({ completed, total: jobs.length, job, index });
        }
        assignNext(worker);
      };

      const handleError = (event) => {
        worker.removeEventListener('message', handleMessage);
        worker.removeEventListener('error', handleError);
        worker.removeEventListener('messageerror', handleMessageError);
        fail(new Error(`${formatPairErrorPrefix(job)}: ${event.message || 'LOSAT worker error'}`));
      };

      const handleMessageError = () => {
        worker.removeEventListener('message', handleMessage);
        worker.removeEventListener('error', handleError);
        worker.removeEventListener('messageerror', handleMessageError);
        fail(new Error(`${formatPairErrorPrefix(job)}: LOSAT worker message could not be decoded`));
      };

      worker.addEventListener('message', handleMessage);
      worker.addEventListener('error', handleError);
      worker.addEventListener('messageerror', handleMessageError);
      worker.postMessage({
        type: 'run',
        id,
        job,
        sequences: sequencePayload
      });
    };

    try {
      console.info(`Running ${workerCount} LOSAT worker${workerCount === 1 ? '' : 's'}.`);
      workers.forEach(assignNext);
    } catch (error) {
      fail(error);
    }
  });
};

const runLosatPairsThreaded = async (
  jobs,
  {
    concurrency,
    onProgress,
    sequences,
    signal,
    threadedWasmPath = DEFAULT_THREADED_WASM_PATH,
    threadedWorkerUrl,
    threadsPerJob,
    wasmModule
  } = {}
) => {
  throwIfAborted(signal);
  const sequenceStore = normalizeSequenceStore(sequences);
  const workerCount = Math.min(
    jobs.length,
    Math.max(1, Number.isFinite(concurrency) ? Math.floor(concurrency) : 1)
  );
  if (workerCount <= 0) return [];

  const startedAt = getNow();
  const resolvedWorkerUrl = threadedWorkerUrl || resolveThreadedWorkerUrl();
  const resolvedWasmUrl = resolveAssetUrl(threadedWasmPath || DEFAULT_THREADED_WASM_PATH);
  const resolvedWasiShimUrl = resolveAssetUrl(WASI_SHIM_URL);
  const compiledModule = wasmModule || await loadLosatModule(threadedWasmPath || DEFAULT_THREADED_WASM_PATH);
  const effectiveThreads = Math.max(1, Number(threadsPerJob) || 1);
  const results = new Array(jobs.length);
  const activeWorkers = new Set();
  let nextJobIndex = 0;
  let completed = 0;
  let requestId = 0;
  let settled = false;

  return new Promise((resolve, reject) => {
    let handleAbort = null;
    const cleanup = () => {
      if (handleAbort) signal?.removeEventListener?.('abort', handleAbort);
      activeWorkers.forEach((worker) => worker.terminate());
      activeWorkers.clear();
    };
    const fail = (error) => {
      if (settled) return;
      settled = true;
      cleanup();
      reject(error);
    };
    const maybeResolve = () => {
      if (settled || completed < jobs.length) return;
      settled = true;
      cleanup();
      console.info(
        `Completed ${jobs.length} threaded LOSAT pair${jobs.length === 1 ? '' : 's'} with ${workerCount} pair worker${workerCount === 1 ? '' : 's'} and ${effectiveThreads} thread${effectiveThreads === 1 ? '' : 's'} per job in ${formatDuration(startedAt)}.`
      );
      resolve(results);
    };
    const launchNext = () => {
      if (settled) return;
      if (signal?.aborted) {
        fail(getAbortReason(signal));
        return;
      }
      if (nextJobIndex >= jobs.length) {
        maybeResolve();
        return;
      }

      const index = nextJobIndex;
      nextJobIndex += 1;
      const job = jobs[index];
      const id = `threaded-${Date.now()}-${requestId}`;
      requestId += 1;
      let payloadJob;
      try {
        payloadJob = materializeJobSequences(job, sequenceStore);
      } catch (error) {
        fail(new Error(`${formatPairErrorPrefix(job)}: ${error?.message || error}`));
        return;
      }

      let worker;
      try {
        worker = new Worker(resolvedWorkerUrl, { type: 'module' });
      } catch (error) {
        fail(new Error(`${formatPairErrorPrefix(job)}: ${error?.message || error}`));
        return;
      }
      activeWorkers.add(worker);

      const cleanupWorker = () => {
        worker.removeEventListener('message', handleMessage);
        worker.removeEventListener('error', handleError);
        worker.removeEventListener('messageerror', handleMessageError);
        worker.terminate();
        activeWorkers.delete(worker);
      };

      const handleMessage = (event) => {
        const data = event.data || {};
        if (data.id !== id) return;
        if (data.type && data.type !== 'run') return;
        cleanupWorker();

        if (!data.ok) {
          fail(new Error(`${formatPairErrorPrefix(job)}: ${data.error || 'Threaded LOSAT worker failed'}`));
          return;
        }
        if (signal?.aborted) {
          fail(getAbortReason(signal));
          return;
        }

        if (data.stderr) {
          console.info(`Threaded LOSAT stderr for ${formatPairErrorPrefix(job)}:\n${data.stderr}`);
        }
        results[index] = { ...job, text: data.text || '' };
        completed += 1;
        if (typeof onProgress === 'function') {
          onProgress({ completed, total: jobs.length, job, index, threaded: true, spawnCount: data.spawnCount || 0 });
        }
        launchNext();
        maybeResolve();
      };

      const handleError = (event) => {
        cleanupWorker();
        fail(new Error(`${formatPairErrorPrefix(job)}: ${event.message || 'Threaded LOSAT worker error'}`));
      };

      const handleMessageError = () => {
        cleanupWorker();
        fail(new Error(`${formatPairErrorPrefix(job)}: Threaded LOSAT worker message could not be decoded`));
      };

      worker.addEventListener('message', handleMessage);
      worker.addEventListener('error', handleError);
      worker.addEventListener('messageerror', handleMessageError);
      worker.postMessage({
        type: 'run',
        id,
        module: compiledModule,
        wasmUrl: resolvedWasmUrl,
        wasiShimUrl: resolvedWasiShimUrl,
        job,
        queryFasta: payloadJob.queryFasta,
        subjectFasta: payloadJob.subjectFasta,
        threadsPerJob: effectiveThreads
      });
    };

    handleAbort = () => fail(getAbortReason(signal));
    if (signal?.aborted) {
      handleAbort();
      return;
    }
    signal?.addEventListener?.('abort', handleAbort, { once: true });
    console.info(
      `Running threaded LOSAT with ${workerCount} pair worker${workerCount === 1 ? '' : 's'} and ${effectiveThreads} thread${effectiveThreads === 1 ? '' : 's'} per job.`
    );
    for (let i = 0; i < workerCount; i += 1) launchNext();
  });
};

export const runLosatPairsParallel = async (jobs, options = {}) => {
  const jobList = Array.isArray(jobs) ? jobs : [];
  if (jobList.length === 0) return [];
  throwIfAborted(options.signal);

  const executionMode = normalizeExecutionMode(options.executionMode);
  const sequenceStore = normalizeSequenceStore(options.sequences);
  let threadedFallbackReason = '';

  if (executionMode !== 'serial') {
    const plan = buildThreadedRuntimePlan(jobList, options, sequenceStore);
    const support = await getLosatThreadingSupport({
      threadedWasmPath: options.threadedWasmPath || DEFAULT_THREADED_WASM_PATH
    });
    throwIfAborted(options.signal);

    if (support.state === 'available' && (executionMode === 'threaded' || plan.useful)) {
      notifyRuntimeStatus(
        options.onRuntimeStatus,
        buildRuntimeStatus(
          'running',
          `Running threaded LOSAT: ${plan.threadsPerJob} thread${plan.threadsPerJob === 1 ? '' : 's'} per job, ${plan.pairWorkers} pair worker${plan.pairWorkers === 1 ? '' : 's'}.`,
          {
            mode: 'threaded',
            threadsPerJob: plan.threadsPerJob,
            pairWorkers: plan.pairWorkers,
            totalBudget: plan.totalBudget
          }
        )
      );
      try {
        const results = await runLosatPairsThreaded(jobList, {
          ...options,
          concurrency: plan.pairWorkers,
          sequences: options.sequences,
          threadedWasmPath: options.threadedWasmPath || DEFAULT_THREADED_WASM_PATH,
          threadsPerJob: plan.threadsPerJob,
          wasmModule: support.wasmModule
        });
        notifyRuntimeStatus(
          options.onRuntimeStatus,
          buildRuntimeStatus(
            'available',
            `Threaded LOSAT completed with ${plan.threadsPerJob} thread${plan.threadsPerJob === 1 ? '' : 's'} per job.`,
            { mode: 'threaded', threadsPerJob: plan.threadsPerJob, pairWorkers: plan.pairWorkers }
          )
        );
        return results;
      } catch (error) {
        if (isAbortError(error, options.signal)) throw getAbortReason(options.signal);
        if (executionMode === 'threaded') throw error;
        threadedFallbackReason = error?.message ? String(error.message) : String(error || 'Threaded LOSAT failed.');
        console.warn('Threaded LOSAT failed; falling back to serial browser execution.', error);
      }
    } else if (executionMode === 'threaded') {
      throw new Error(support.message || 'Threaded LOSAT is unavailable.');
    } else {
      threadedFallbackReason = support.state === 'available'
        ? 'Current LOSAT workload is below the threaded auto threshold.'
        : support.message || 'Threaded LOSAT is unavailable.';
    }

    if (threadedFallbackReason) {
      notifyRuntimeStatus(
        options.onRuntimeStatus,
        buildRuntimeStatus(
          'fallback',
          `Using serial LOSAT: ${threadedFallbackReason}`,
          { mode: 'serial', fallbackReason: threadedFallbackReason }
        )
      );
    }
  } else {
    notifyRuntimeStatus(
      options.onRuntimeStatus,
      buildRuntimeStatus('disabled', 'Using serial LOSAT by request.', { mode: 'serial' })
    );
  }

  try {
    return await runLosatPairsWithWorkers(jobList, options);
  } catch (error) {
    if (isAbortError(error, options.signal)) throw getAbortReason(options.signal);
    console.warn('LOSAT Worker pool failed; falling back to sequential execution.', error);
    return runLosatPairsSequential(jobList, options);
  }
};

export const prepareLosatRuntime = async ({
  wasmPath = DEFAULT_WASM_PATH,
  includeThreaded = false,
  threadedWasmPath = DEFAULT_THREADED_WASM_PATH
} = {}) => {
  const prepared = await Promise.all([
    loadWasiShim(),
    loadLosatModule(wasmPath)
  ]);
  if (!includeThreaded) return { wasiShim: prepared[0], wasmModule: prepared[1], threaded: null };
  const threaded = await getLosatThreadingSupport({ threadedWasmPath });
  return { wasiShim: prepared[0], wasmModule: prepared[1], threaded };
};
