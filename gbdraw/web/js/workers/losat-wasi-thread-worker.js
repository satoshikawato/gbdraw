import { concatUint8Arrays } from '../services/losat-runtime.js';

const markWorkerStart = (control, state) => {
  const view = new Int32Array(control);
  Atomics.store(view, 0, state);
  Atomics.notify(view, 0, 1);
};

let preparedContext = null;
let preparedPayload = null;
let readyControl = null;

const markSlotState = (state) => {
  if (!readyControl) return;
  markWorkerStart(readyControl, state);
};

const instantiateThreadContext = async ({
  module,
  wasmUrl,
  memory,
  args,
  env,
  wasiShimUrl
}) => {
  if (!memory || !(memory.buffer instanceof SharedArrayBuffer)) {
    throw new Error('LOSAT WASI thread worker requires shared Wasm memory.');
  }
  const wasiShim = await import(wasiShimUrl);
  const { WASI, File, OpenFile, PreopenDirectory, ConsoleStdout } = wasiShim;
  const stdoutChunks = [];
  const stderrChunks = [];
  const stdout = new ConsoleStdout((data) => stdoutChunks.push(data));
  const stderr = new ConsoleStdout((data) => stderrChunks.push(data));
  const stdin = new OpenFile(new File(new Uint8Array(), { readonly: true }));
  const preopen = new PreopenDirectory('.', new Map());
  const wasi = new WASI(args || [], env || [], [stdin, stdout, stderr, preopen]);
  const compiledModule = module || await WebAssembly.compileStreaming(fetch(wasmUrl, { cache: 'no-store' }));
  const instance = await WebAssembly.instantiate(compiledModule, {
    env: { memory },
    wasi: { 'thread-spawn': () => -1 },
    wasi_snapshot_preview1: wasi.wasiImport
  });
  if (typeof wasi.initialize === 'function') {
    wasi.initialize(instance);
  }
  if (typeof instance.exports.wasi_thread_start !== 'function') {
    throw new Error('Threaded LOSAT wasm does not export wasi_thread_start.');
  }
  return { instance, stdoutChunks, stderrChunks };
};

const flushThreadOutput = (context, payload = {}) => {
  if (!context) return;
  const stdoutText = new TextDecoder().decode(concatUint8Arrays(context.stdoutChunks));
  const stderrText = new TextDecoder().decode(concatUint8Arrays(context.stderrChunks));
  self.postMessage({ ...payload, stdout: stdoutText, stderr: stderrText });
};

const startPreparedThread = ({ id, tid, startArg, control }) => {
  (async () => {
    try {
      if (!preparedContext) {
        throw new Error('LOSAT WASI thread worker was not prepared before start.');
      }
      markWorkerStart(control, 1);
      preparedContext.instance.exports.wasi_thread_start(tid, startArg);
      flushThreadOutput(preparedContext, { id, type: 'done', tid });
      preparedContext = null;
      markSlotState(0);
      if (!preparedPayload) {
        throw new Error('LOSAT WASI thread worker cannot reprepare without a saved payload.');
      }
      preparedContext = await instantiateThreadContext(preparedPayload);
      markSlotState(1);
    } catch (error) {
      preparedContext = null;
      markWorkerStart(control, -1);
      markSlotState(-1);
      self.postMessage({
        id,
        type: 'error',
        tid,
        error: error?.message ? String(error.message) : String(error || 'LOSAT WASI thread failed')
      });
    }
  })();
};

const runThreadFromCold = async (data) => {
  const {
    module,
    wasmUrl,
    memory,
    startArg,
    tid,
    control,
    args,
    env,
    wasiShimUrl
  } = data || {};

  try {
    const context = await instantiateThreadContext({ module, wasmUrl, memory, args, env, wasiShimUrl });
    markWorkerStart(control, 1);
    context.instance.exports.wasi_thread_start(tid, startArg);
    flushThreadOutput(context, { type: 'done', tid });
    self.close();
  } catch (error) {
    markWorkerStart(control, -1);
    self.postMessage({
      type: 'error',
      tid,
      error: error?.message ? String(error.message) : String(error || 'LOSAT WASI thread failed')
    });
    self.close();
  }
};

self.onmessage = async (event) => {
  const data = event.data || {};
  const { id, type } = data;

  if (type === 'prepare') {
    readyControl = data.readyControl || null;
    preparedPayload = {
      module: data.module,
      wasmUrl: data.wasmUrl,
      memory: data.memory,
      args: data.args,
      env: data.env,
      wasiShimUrl: data.wasiShimUrl
    };
    markSlotState(0);
    try {
      preparedContext = await instantiateThreadContext(preparedPayload);
      markSlotState(1);
      self.postMessage({ id, type: 'prepare', ok: true });
    } catch (error) {
      preparedContext = null;
      preparedPayload = null;
      markSlotState(-1);
      self.postMessage({
        id,
        type: 'prepare',
        ok: false,
        error: error?.message ? String(error.message) : String(error || 'LOSAT WASI thread prepare failed')
      });
    }
    return;
  }

  if (type === 'start') {
    startPreparedThread(data);
    return;
  }

  await runThreadFromCold(data);
};
