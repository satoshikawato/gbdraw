import { WASI_SHIM_URL } from '../config.js';

const DEFAULT_WASM_PATH = './wasm/losat/losat.wasm';
const DEFAULT_MAX_WORKERS = 4;

let wasiShimPromise = null;
let wasmModulePromise = null;
let directInstancePromise = null;

const resolveAssetUrl = (path) => new URL(path, window.location.href).toString();
const resolveWorkerUrl = () => new URL('../workers/losat-worker.js', import.meta.url).toString();
const getNow = () => (globalThis.performance?.now ? performance.now() : Date.now());
const formatDuration = (startedAt) => `${((getNow() - startedAt) / 1000).toFixed(2)}s`;

export const createAbortError = () => {
  if (typeof DOMException === 'function') {
    return new DOMException('Generation cancelled', 'AbortError');
  }
  const error = new Error('Generation cancelled');
  error.name = 'AbortError';
  return error;
};

export const isAbortError = (error) => error?.name === 'AbortError';

const throwIfAborted = (signal) => {
  if (signal?.aborted) throw createAbortError();
};

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
  if (!wasmModulePromise) {
    wasmModulePromise = (async () => {
      const resolvedWasmPath = resolveAssetUrl(wasmPath);
      const response = await fetch(resolvedWasmPath, { cache: 'no-store' });
      if (!response.ok) {
        throw new Error(`Missing packaged asset: LOSAT wasm not found at ${resolvedWasmPath}`);
      }
      const bytes = await response.arrayBuffer();
      return WebAssembly.compile(bytes);
    })();
  }
  return wasmModulePromise;
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
  wasmPath = DEFAULT_WASM_PATH
} = {}) => {
  if (!program || (program !== 'blastn' && program !== 'tblastx')) {
    throw new Error('LOSAT program must be blastn or tblastx.');
  }
  if (!queryFasta || !subjectFasta) {
    throw new Error('LOSAT requires both query and subject FASTA content.');
  }

  const [wasiShim, wasmModule] = await Promise.all([
    loadWasiShim(),
    loadLosatModule(wasmPath)
  ]);
  const { WASI, File, OpenFile, PreopenDirectory, ConsoleStdout } = wasiShim;

  if (hasDirectLosatApi(wasmModule)) {
    return runLosatPairDirect({
      program,
      queryFasta,
      subjectFasta,
      outfmt,
      extraArgs,
      wasiShim,
      wasmModule
    });
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
  const exitCode = wasi.start(instance);

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

const formatPairErrorPrefix = (job) => {
  const pairNumber = Number.isInteger(job?.pairIndex) ? job.pairIndex + 1 : null;
  return pairNumber ? `LOSAT pair #${pairNumber}` : 'LOSAT pair';
};

const runLosatPairsSequential = async (jobs, {
  signal,
  onQueued,
  onStarted,
  onCompleted,
  onProgress
} = {}) => {
  const startedAt = getNow();
  console.info(`Running ${jobs.length} LOSAT pair${jobs.length === 1 ? '' : 's'} sequentially.`);
  const results = [];
  let completed = 0;
  onQueued?.({ jobs: [...jobs], activeJobs: [], completed, total: jobs.length });
  onProgress?.({ activeJobs: [], completed, total: jobs.length });
  for (const job of jobs) {
    throwIfAborted(signal);
    try {
      const activeJobs = [job];
      onStarted?.({ job, activeJobs, completed, total: jobs.length });
      onProgress?.({ activeJobs, completed, total: jobs.length });
      const text = await runLosatPair(job);
      throwIfAborted(signal);
      const result = { ...job, text };
      results.push(result);
      completed += 1;
      onCompleted?.({ result, activeJobs: [], completed, total: jobs.length });
      onProgress?.({ activeJobs: [], completed, total: jobs.length });
    } catch (error) {
      if (isAbortError(error)) throw error;
      const message = error?.message ? String(error.message) : String(error || 'Unknown LOSAT error');
      throw new Error(`${formatPairErrorPrefix(job)}: ${message}`);
    }
  }
  console.info(`Completed ${jobs.length} LOSAT pair${jobs.length === 1 ? '' : 's'} sequentially in ${formatDuration(startedAt)}.`);
  return results;
};

const runLosatPairsWithWorkers = async (jobs, {
  concurrency,
  workerUrl,
  wasmPath,
  signal,
  onQueued,
  onStarted,
  onCompleted,
  onProgress
} = {}) => {
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
  const activeJobs = new Map();
  let nextJobIndex = 0;
  let completed = 0;
  let settled = false;
  let requestId = 0;

  return new Promise((resolve, reject) => {
    const getActiveJobs = () => Array.from(activeJobs.values()).map((entry) => entry.job);
    const emitProgress = () => {
      onProgress?.({ activeJobs: getActiveJobs(), completed, total: jobs.length });
    };

    const cleanup = () => {
      signal?.removeEventListener('abort', handleAbort);
      workers.forEach((worker) => worker.terminate());
    };

    const fail = (error) => {
      if (settled) return;
      settled = true;
      cleanup();
      reject(error);
    };

    const handleAbort = () => {
      console.info(
        `Cancelling LOSAT worker pool after ${completed} / ${jobs.length} completed job${jobs.length === 1 ? '' : 's'}.`
      );
      fail(createAbortError());
    };

    if (signal?.aborted) {
      fail(createAbortError());
      return;
    }
    signal?.addEventListener('abort', handleAbort, { once: true });

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
        fail(createAbortError());
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
      activeJobs.set(worker, { job, index, id });
      onStarted?.({ job, activeJobs: getActiveJobs(), completed, total: jobs.length });
      emitProgress();

      const handleMessage = (event) => {
        const data = event.data || {};
        if (data.id !== id) return;
        if (settled) return;
        worker.removeEventListener('message', handleMessage);
        worker.removeEventListener('error', handleError);
        worker.removeEventListener('messageerror', handleMessageError);
        activeJobs.delete(worker);

        if (!data.ok) {
          fail(new Error(`${formatPairErrorPrefix(job)}: ${data.error || 'LOSAT worker failed'}`));
          return;
        }

        const result = { ...job, text: data.text || '' };
        results[index] = result;
        completed += 1;
        onCompleted?.({ result, activeJobs: getActiveJobs(), completed, total: jobs.length });
        emitProgress();
        assignNext(worker);
      };

      const handleError = (event) => {
        worker.removeEventListener('message', handleMessage);
        worker.removeEventListener('error', handleError);
        worker.removeEventListener('messageerror', handleMessageError);
        activeJobs.delete(worker);
        fail(new Error(`${formatPairErrorPrefix(job)}: ${event.message || 'LOSAT worker error'}`));
      };

      const handleMessageError = () => {
        worker.removeEventListener('message', handleMessage);
        worker.removeEventListener('error', handleError);
        worker.removeEventListener('messageerror', handleMessageError);
        activeJobs.delete(worker);
        fail(new Error(`${formatPairErrorPrefix(job)}: LOSAT worker message could not be decoded`));
      };

      worker.addEventListener('message', handleMessage);
      worker.addEventListener('error', handleError);
      worker.addEventListener('messageerror', handleMessageError);
      worker.postMessage({
        id,
        job: {
          ...job,
          wasmModule,
          wasmUrl: resolvedWasmUrl,
          wasiShimUrl: resolvedWasiShimUrl
        }
      });
    };

    try {
      for (let i = 0; i < workerCount; i += 1) {
        workers.push(new Worker(resolvedWorkerUrl, { type: 'module' }));
      }
      console.info(`Running ${workerCount} LOSAT worker${workerCount === 1 ? '' : 's'}.`);
      onQueued?.({ jobs: [...jobs], activeJobs: [], completed, total: jobs.length });
      emitProgress();
      workers.forEach(assignNext);
    } catch (error) {
      fail(error);
    }
  });
};

export const runLosatPairsParallel = async (jobs, options = {}) => {
  const jobList = Array.isArray(jobs) ? jobs : [];
  if (jobList.length === 0) return [];
  throwIfAborted(options.signal);

  try {
    return await runLosatPairsWithWorkers(jobList, options);
  } catch (error) {
    if (isAbortError(error)) throw error;
    console.warn('LOSAT Worker pool failed; falling back to sequential execution.', error);
    return runLosatPairsSequential(jobList, options);
  }
};

export const prepareLosatRuntime = async ({ wasmPath = DEFAULT_WASM_PATH } = {}) => {
  await Promise.all([
    loadWasiShim(),
    loadLosatModule(wasmPath)
  ]);
};
