const markWorkerStart = (control, state) => {
  const view = new Int32Array(control);
  Atomics.store(view, 0, state);
  Atomics.notify(view, 0, 1);
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

self.onmessage = async (event) => {
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
  } = event.data || {};

  try {
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

    markWorkerStart(control, 1);
    instance.exports.wasi_thread_start(tid, startArg);

    const stdoutText = new TextDecoder().decode(concatUint8Arrays(stdoutChunks));
    const stderrText = new TextDecoder().decode(concatUint8Arrays(stderrChunks));
    self.postMessage({ type: 'done', tid, stdout: stdoutText, stderr: stderrText });
  } catch (error) {
    markWorkerStart(control, -1);
    self.postMessage({
      type: 'error',
      tid,
      error: error?.message ? String(error.message) : String(error || 'LOSAT WASI thread failed')
    });
  }
};
