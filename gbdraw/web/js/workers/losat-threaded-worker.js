import { compileLosatWasm, concatUint8Arrays, losatCommandArgs } from '../services/losat-runtime.js';
import { LOSAT_WASM_BUILD } from '../services/losat-build.generated.js';

// The shared memory of a job: the initial pages the threaded file imports (recorded from the file), and at most
// 1 GiB (16384 pages) reserved, below the file's own maximum.
const MEMORY_INITIAL_PAGES = LOSAT_WASM_BUILD.threadedMemoryPages.initial;
const MEMORY_MAXIMUM_PAGES = Math.min(16384, LOSAT_WASM_BUILD.threadedMemoryPages.maximum);
const WORKER_PREPARE_TIMEOUT_MS = 120000;
const WORKER_START_TIMEOUT_MS = 30000;

const childWorkerUrl = new URL('./losat-wasi-thread-worker.js', import.meta.url).toString();

const compileLosatModule = async (module, wasmUrl) => {
  if (module instanceof WebAssembly.Module) return module;
  return compileLosatWasm(wasmUrl, 'threaded LOSAT wasm');
};

const getChildWorkerCount = (threadsPerJob) =>
  Math.max(0, Math.floor(Math.max(1, Number(threadsPerJob) || 1)) - 1);

const instantiateWithMemory = async ({ module, memory, wasi, spawnThread }) =>
  WebAssembly.instantiate(module, {
    env: { memory },
    wasi: { 'thread-spawn': spawnThread },
    wasi_snapshot_preview1: wasi.wasiImport
  });

const instantiateMain = async ({ module, wasi, spawnThread }) => {
  const memory = new WebAssembly.Memory({ initial: MEMORY_INITIAL_PAGES, maximum: MEMORY_MAXIMUM_PAGES, shared: true });
  const instance = await instantiateWithMemory({ module, memory, wasi, spawnThread });
  return { instance, memory };
};

const terminateWorkers = async (workers) => {
  const pending = [];
  for (const worker of workers) pending.push(worker.terminate());
  workers.clear();
  await Promise.allSettled(pending);
};

const attachThreadWorkerOutput = ({ worker, encoder, stdoutChunks, stderrChunks }) => {
  worker.addEventListener('message', (event) => {
    const data = event.data || {};
    if (data.stderr) stderrChunks.push(encoder.encode(String(data.stderr)));
    if (data.stdout) stdoutChunks.push(encoder.encode(String(data.stdout)));
    if (data.type === 'error' && data.error) {
      stderrChunks.push(encoder.encode(`${data.error}\n`));
    }
  });
  worker.addEventListener('error', (event) => {
    stderrChunks.push(encoder.encode(`${event.message || 'LOSAT WASI thread worker error'}\n`));
  });
};

const prepareThreadWorker = ({
  compiledModule,
  wasmUrl,
  mainMemory,
  args,
  env,
  wasiShimUrl,
  faultChannel,
  workers,
  encoder,
  stdoutChunks,
  stderrChunks,
  timeoutMs = WORKER_PREPARE_TIMEOUT_MS
}) =>
  new Promise((resolve, reject) => {
    let settled = false;
    let worker;
    const readyControl = new SharedArrayBuffer(Int32Array.BYTES_PER_ELEMENT);
    const readyView = new Int32Array(readyControl);
    const id = `prepare-${Date.now()}-${Math.random().toString(36).slice(2)}`;
    const settle = (callback) => {
      if (settled) return;
      settled = true;
      cleanup();
      callback();
    };
    const cleanup = () => {
      clearTimeout(timeout);
      if (!worker) return;
      worker.removeEventListener('message', handleMessage);
      worker.removeEventListener('error', handleError);
      worker.removeEventListener('messageerror', handleMessageError);
    };
    const handleMessage = (event) => {
      const data = event.data || {};
      if (data.id !== id || data.type !== 'prepare') return;
      if (data.ok) {
        settle(() => resolve({ worker, readyView }));
        return;
      }
      settle(() => {
        worker.terminate();
        workers.delete(worker);
        reject(new Error(data.error || 'LOSAT WASI thread worker preparation failed'));
      });
    };
    const handleError = (event) => {
      settle(() => {
        if (worker) {
          worker.terminate();
          workers.delete(worker);
        }
        reject(new Error(event.message || 'LOSAT WASI thread worker preparation error'));
      });
    };
    const handleMessageError = () => {
      settle(() => {
        if (worker) {
          worker.terminate();
          workers.delete(worker);
        }
        reject(new Error('LOSAT WASI thread worker preparation message could not be decoded'));
      });
    };
    const timeout = setTimeout(() => {
      settle(() => {
        if (worker) {
          worker.terminate();
          workers.delete(worker);
        }
        reject(new Error(`timed out preparing LOSAT WASI thread worker after ${Math.round(timeoutMs / 1000)}s`));
      });
    }, timeoutMs);

    try {
      worker = new Worker(childWorkerUrl, { type: 'module' });
      workers.add(worker);
      attachThreadWorkerOutput({ worker, encoder, stdoutChunks, stderrChunks });
      worker.addEventListener('message', handleMessage);
      worker.addEventListener('error', handleError);
      worker.addEventListener('messageerror', handleMessageError);
      worker.postMessage({
        type: 'prepare',
        id,
        readyControl,
        module: compiledModule,
        wasmUrl,
        memory: mainMemory,
        args,
        env,
        wasiShimUrl,
        faultChannel
      });
    } catch (error) {
      settle(() => {
        if (worker) {
          worker.terminate();
          workers.delete(worker);
        }
        reject(error);
      });
    }
  });

const runThreadedLosat = async ({
  module,
  wasmUrl,
  wasiShimUrl,
  job,
  queryFasta,
  subjectFasta,
  threadsPerJob,
  faultChannel
}) => {
  if (typeof SharedArrayBuffer !== 'function') {
    throw new Error('SharedArrayBuffer is unavailable; threaded LOSAT requires cross-origin isolation.');
  }
  const wasiShim = await import(wasiShimUrl);
  const { WASI, File, OpenFile, PreopenDirectory, ConsoleStdout } = wasiShim;
  const compiledModule = await compileLosatModule(module, wasmUrl);
  const exports = WebAssembly.Module.exports(compiledModule);
  const imports = WebAssembly.Module.imports(compiledModule);
  if (!exports.some((entry) => entry.name === '_start') ||
      !exports.some((entry) => entry.name === 'wasi_thread_start') ||
      !imports.some((entry) => entry.module === 'wasi' && entry.name === 'thread-spawn')) {
    throw new Error('Threaded LOSAT wasm is missing _start, wasi_thread_start, or wasi/thread-spawn.');
  }

  const encoder = new TextEncoder();
  const files = new Map([
    ['query.fa', new File(encoder.encode(queryFasta), { readonly: true })],
    ['subject.fa', new File(encoder.encode(subjectFasta), { readonly: true })]
  ]);
  const stdoutChunks = [];
  const stderrChunks = [];
  const stdout = new ConsoleStdout((data) => stdoutChunks.push(data));
  const stderr = new ConsoleStdout((data) => stderrChunks.push(data));
  const stdin = new OpenFile(new File(new Uint8Array(), { readonly: true }));
  const preopen = new PreopenDirectory('.', files);
  const effectiveThreads = Math.max(1, Number(threadsPerJob) || 1);
  const childWorkerCount = getChildWorkerCount(effectiveThreads);
  const args = losatCommandArgs({
    program: job.program,
    outfmt: job.outfmt,
    extraArgs: job.extraArgs,
    threads: effectiveThreads
  });
  const env = [
    `LOSAT_WASI_THREAD_CAP=${effectiveThreads}`,
    `RAYON_NUM_THREADS=${effectiveThreads}`,
    `LOSAT_WASM_THREADS_BROWSER=1`
  ];
  const wasi = new WASI(args, env, [stdin, stdout, stderr, preopen]);
  const workers = new Set();
  const threadSlots = [];
  let nextTid = 1;
  let spawnCount = 0;
  let mainMemory = null;

  const waitForThreadStart = ({ tid, worker, controlView, timeoutMs }) => {
    const wait = Atomics.wait(controlView, 0, 0, timeoutMs);
    if (wait === 'timed-out') {
      stderrChunks.push(encoder.encode(`timed out waiting for LOSAT WASI thread ${tid} to start after ${Math.round(timeoutMs / 1000)}s\n`));
      worker.terminate();
      workers.delete(worker);
      return -1;
    }
    if (Atomics.load(controlView, 0) !== 1) {
      worker.terminate();
      workers.delete(worker);
      return -1;
    }
    return tid;
  };

  const acquirePreparedThreadSlot = (timeoutMs) => {
    if (threadSlots.length === 0) return null;
    const deadline = performance.now() + timeoutMs;
    while (performance.now() < deadline) {
      for (const slot of threadSlots) {
        const state = Atomics.load(slot.readyView, 0);
        if (state === -1) {
          stderrChunks.push(encoder.encode('LOSAT WASI thread worker entered an error state\n'));
          return null;
        }
        if (state === 1 && Atomics.compareExchange(slot.readyView, 0, 1, 2) === 1) {
          return slot;
        }
      }
      const preparingSlot = threadSlots.find((slot) => Atomics.load(slot.readyView, 0) === 0);
      if (!preparingSlot) return null;
      const remaining = Math.max(1, Math.min(250, deadline - performance.now()));
      Atomics.wait(preparingSlot.readyView, 0, 0, remaining);
    }
    stderrChunks.push(encoder.encode(`timed out waiting for reusable LOSAT WASI thread worker slot after ${Math.round(timeoutMs / 1000)}s\n`));
    return null;
  };

  const spawnPreparedThread = ({ tid, startArg, control, controlView }) => {
    const slot = acquirePreparedThreadSlot(WORKER_PREPARE_TIMEOUT_MS);
    if (!slot) return null;
    slot.worker.postMessage({
      type: 'start',
      id: `start-${tid}-${Date.now()}`,
      tid,
      startArg,
      control
    });
    const startedTid = waitForThreadStart({ tid, worker: slot.worker, controlView, timeoutMs: WORKER_START_TIMEOUT_MS });
    if (startedTid === -1) {
      Atomics.store(slot.readyView, 0, -1);
      Atomics.notify(slot.readyView, 0, 1);
    }
    return startedTid;
  };

  const spawnColdThread = ({ tid, startArg, control, controlView }) => {
    let worker;
    try {
      worker = new Worker(childWorkerUrl, { type: 'module' });
    } catch (error) {
      stderrChunks.push(encoder.encode(`${error?.message || error}\n`));
      return -1;
    }
    workers.add(worker);
    attachThreadWorkerOutput({ worker, encoder, stdoutChunks, stderrChunks });
    worker.postMessage({
      type: 'run',
      module: compiledModule,
      wasmUrl,
      memory: mainMemory,
      startArg,
      tid,
      control,
      args,
      env,
      wasiShimUrl,
      faultChannel
    });
    return waitForThreadStart({ tid, worker, controlView, timeoutMs: WORKER_PREPARE_TIMEOUT_MS });
  };

  const spawnThread = (startArg) => {
    const tid = nextTid;
    nextTid += 1;
    spawnCount += 1;
    const control = new SharedArrayBuffer(Int32Array.BYTES_PER_ELEMENT);
    const controlView = new Int32Array(control);
    const preparedTid = spawnPreparedThread({ tid, startArg, control, controlView });
    if (preparedTid !== null) return preparedTid;
    return spawnColdThread({ tid, startArg, control, controlView });
  };

  try {
    const instantiated = await instantiateMain({
      module: compiledModule,
      wasi,
      spawnThread
    });
    mainMemory = instantiated.memory;
    const { instance } = instantiated;
    if (typeof instance.exports._start !== 'function') {
      throw new Error('Threaded LOSAT wasm does not export _start.');
    }
    if (childWorkerCount > 0) {
      const prepared = await Promise.all(
        Array.from({ length: childWorkerCount }, () =>
          prepareThreadWorker({
            compiledModule,
            wasmUrl,
            mainMemory,
            args,
            env,
            wasiShimUrl,
            faultChannel,
            workers,
            encoder,
            stdoutChunks,
            stderrChunks
          })
        )
      );
      threadSlots.push(...prepared);
    }
    const exitCode = wasi.start(instance);
    const stdoutText = new TextDecoder().decode(concatUint8Arrays(stdoutChunks));
    const stderrText = new TextDecoder().decode(concatUint8Arrays(stderrChunks));
    if (exitCode !== 0) {
      throw new Error(stderrText.trim() || `LOSAT exited with code ${exitCode}`);
    }
    return { text: stdoutText, stderr: stderrText, spawnCount };
  } finally {
    await terminateWorkers(workers);
  }
};

self.onmessage = async (event) => {
  const { id, type } = event.data || {};
  if (type !== 'run') {
    self.postMessage({ id, type: type || 'run', ok: false, error: `Unsupported threaded LOSAT worker message type '${type || '(blank)'}'.` });
    return;
  }

  try {
    const result = await runThreadedLosat(event.data || {});
    self.postMessage({
      id,
      type: 'run',
      ok: true,
      text: result.text,
      stderr: result.stderr,
      spawnCount: result.spawnCount
    });
  } catch (error) {
    self.postMessage({
      id,
      type: 'run',
      ok: false,
      error: error?.message ? String(error.message) : String(error || 'Threaded LOSAT worker failed'),
      trap: error instanceof WebAssembly.RuntimeError
    });
  }
};
