import { WASI_SHIM_URL } from '../config.js';

const DEFAULT_WASM_PATH = './wasm/losat/losat.wasm';
const DEFAULT_MAX_WORKERS = 4;

let wasiShimPromise = null;
let wasmModulePromise = null;

const resolveAssetUrl = (path) => new URL(path, window.location.href).toString();
const resolveWorkerUrl = () => new URL('../workers/losat-worker.js', import.meta.url).toString();
const getNow = () => (globalThis.performance?.now ? performance.now() : Date.now());
const formatDuration = (startedAt) => `${((getNow() - startedAt) / 1000).toFixed(2)}s`;

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

  const [{ WASI, File, OpenFile, PreopenDirectory, ConsoleStdout }, wasmModule] = await Promise.all([
    loadWasiShim(),
    loadLosatModule(wasmPath)
  ]);

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

const runLosatPairsSequential = async (jobs) => {
  const startedAt = getNow();
  console.info(`Running ${jobs.length} LOSAT pair${jobs.length === 1 ? '' : 's'} sequentially.`);
  const results = [];
  for (const job of jobs) {
    try {
      const text = await runLosatPair(job);
      results.push({ ...job, text });
    } catch (error) {
      const message = error?.message ? String(error.message) : String(error || 'Unknown LOSAT error');
      throw new Error(`${formatPairErrorPrefix(job)}: ${message}`);
    }
  }
  console.info(`Completed ${jobs.length} LOSAT pair${jobs.length === 1 ? '' : 's'} sequentially in ${formatDuration(startedAt)}.`);
  return results;
};

const runLosatPairsWithWorkers = async (jobs, { concurrency, workerUrl, wasmPath } = {}) => {
  const workerCount = Math.min(
    jobs.length,
    Math.max(1, Number.isFinite(concurrency) ? Math.floor(concurrency) : getDefaultConcurrency(jobs.length))
  );
  if (workerCount <= 0) return [];

  const startedAt = getNow();
  const resolvedWorkerUrl = workerUrl || resolveWorkerUrl();
  const resolvedWasmUrl = resolveAssetUrl(wasmPath || DEFAULT_WASM_PATH);
  const resolvedWasiShimUrl = resolveAssetUrl(WASI_SHIM_URL);
  const wasmModule = await loadLosatModule(wasmPath || DEFAULT_WASM_PATH);
  const results = new Array(jobs.length);
  const workers = [];
  let nextJobIndex = 0;
  let completed = 0;
  let settled = false;
  let requestId = 0;

  return new Promise((resolve, reject) => {
    const cleanup = () => {
      workers.forEach((worker) => worker.terminate());
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
        `Completed ${jobs.length} LOSAT pair${jobs.length === 1 ? '' : 's'} with ${workerCount} worker${workerCount === 1 ? '' : 's'} in ${formatDuration(startedAt)}.`
      );
      resolve(results);
    };

    const assignNext = (worker) => {
      if (settled) return;
      if (nextJobIndex >= jobs.length) {
        maybeResolve();
        return;
      }
      const index = nextJobIndex;
      nextJobIndex += 1;
      const job = jobs[index];
      const id = `${Date.now()}-${requestId}`;
      requestId += 1;

      const handleMessage = (event) => {
        const data = event.data || {};
        if (data.id !== id) return;
        worker.removeEventListener('message', handleMessage);
        worker.removeEventListener('error', handleError);
        worker.removeEventListener('messageerror', handleMessageError);

        if (!data.ok) {
          fail(new Error(`${formatPairErrorPrefix(job)}: ${data.error || 'LOSAT worker failed'}`));
          return;
        }

        results[index] = { ...job, text: data.text || '' };
        completed += 1;
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
      workers.forEach(assignNext);
    } catch (error) {
      fail(error);
    }
  });
};

export const runLosatPairsParallel = async (jobs, options = {}) => {
  const jobList = Array.isArray(jobs) ? jobs : [];
  if (jobList.length === 0) return [];
  if (jobList.length === 1) return runLosatPairsSequential(jobList);

  try {
    return await runLosatPairsWithWorkers(jobList, options);
  } catch (error) {
    console.warn('LOSAT Worker pool failed; falling back to sequential execution.', error);
    return runLosatPairsSequential(jobList);
  }
};
