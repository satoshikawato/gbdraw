// @ts-check
/**
 * The exports of the direct LOSAT wasm API (`hasDirectLosatApi` checks for it).
 * @typedef {object} LosatDirectExports
 * @property {WebAssembly.Memory} memory
 * @property {(len: number) => number} losat_web_alloc
 * @property {(ptr: number, len: number) => void} losat_web_dealloc
 * @property {(...args: number[]) => number} losat_web_run_pair
 * @property {() => number} losat_web_error_ptr
 * @property {() => number} losat_web_error_len
 * @property {() => number} losat_web_result_ptr
 * @property {() => number} losat_web_result_len
 * @property {() => void} [losat_web_clear]
 */

/**
 * What this module uses of the vendored WASI shim (`vendor/browser_wasi_shim`): five classes. A
 * file and a file descriptor are only handed back to the shim, so they are opaque here.
 * @typedef {object} LosatWasiShim
 * @property {new (args: string[], env: string[], fds: unknown[]) => {
 *   wasiImport: WebAssembly.ModuleImports,
 *   start: (instance: WebAssembly.Instance) => number,
 *   initialize?: (instance: WebAssembly.Instance) => void
 * }} WASI
 * @property {new (data: Uint8Array, options?: { readonly?: boolean }) => unknown} File
 * @property {new (file: unknown) => unknown} OpenFile
 * @property {new (name: string, contents: Map<string, unknown>) => unknown} PreopenDirectory
 * @property {new (write: (data: Uint8Array) => void) => unknown} ConsoleStdout
 */

/**
 * One pair run, by the main thread or a pool worker.
 * @typedef {object} LosatPairRun
 * @property {string} program 'blastn', 'tblastx', or 'blastp'
 * @property {string} queryFasta
 * @property {string} subjectFasta
 * @property {string} [outfmt]
 * @property {string[]} [extraArgs]
 * @property {LosatWasiShim} wasiShim
 * @property {WebAssembly.Module} wasmModule
 */

/** @type {Promise<WebAssembly.Instance> | null} */
let directInstancePromise = null;

/**
 * Fetch and compile a packaged LOSAT wasm. Its URL carries no version, so the
 * request always revalidates (`no-cache`): an unchanged file costs a 304, and a
 * file replaced at the same URL is never served from a stale cache.
 * @param {string} wasmUrl
 * @param {string} label the asset named when it is missing
 * @returns {Promise<WebAssembly.Module>}
 */
export const compileLosatWasm = async (wasmUrl, label) => {
  const response = await fetch(wasmUrl, { cache: 'no-cache' });
  if (!response.ok) {
    throw new Error(`Missing packaged asset: ${label} not found at ${wasmUrl}`);
  }
  return WebAssembly.compile(await response.arrayBuffer());
};

/** @param {Uint8Array[]} chunks */
export const concatUint8Arrays = (chunks) => {
  const total = chunks.reduce((sum, chunk) => sum + chunk.length, 0);
  const merged = new Uint8Array(total);
  let offset = 0;
  chunks.forEach((chunk) => {
    merged.set(chunk, offset);
    offset += chunk.length;
  });
  return merged;
};

/** @param {WebAssembly.Module} wasmModule */
export const hasDirectLosatApi = (wasmModule) =>
  WebAssembly.Module.exports(wasmModule).some(
    (entry) => entry.kind === 'function' && entry.name === 'losat_web_run_pair'
  );

/**
 * @param {LosatPairRun & { checkCanceled?: () => void }} run
 * @returns {Promise<string>}
 */
export const runLosatPairWasi = async ({
  program,
  queryFasta,
  subjectFasta,
  outfmt = '6',
  extraArgs = /** @type {string[]} */ ([]),
  wasiShim,
  wasmModule,
  checkCanceled = () => {}
}) => {
  const { WASI, File, OpenFile, PreopenDirectory, ConsoleStdout } = wasiShim;
  const encoder = new TextEncoder();
  const files = new Map([
    ['query.fa', new File(encoder.encode(queryFasta), { readonly: true })],
    ['subject.fa', new File(encoder.encode(subjectFasta), { readonly: true })]
  ]);
  const preopen = new PreopenDirectory('.', files);
  /** @type {Uint8Array[]} */
  const stdoutChunks = [];
  /** @type {Uint8Array[]} */
  const stderrChunks = [];
  const stdout = new ConsoleStdout((/** @type {Uint8Array} */ data) => stdoutChunks.push(data));
  const stderr = new ConsoleStdout((/** @type {Uint8Array} */ data) => stderrChunks.push(data));
  const stdin = new OpenFile(new File(new Uint8Array(), { readonly: true }));
  const args = [
    'losat',
    program,
    '--query',
    'query.fa',
    '--subject',
    'subject.fa',
    '--outfmt',
    String(outfmt),
    ...(Array.isArray(extraArgs) ? extraArgs : [])
  ];
  const wasi = new WASI(args, [], [stdin, stdout, stderr, preopen]);
  const instance = await WebAssembly.instantiate(wasmModule, {
    wasi_snapshot_preview1: wasi.wasiImport
  });
  checkCanceled();
  const exitCode = wasi.start(instance);
  checkCanceled();

  const stdoutText = new TextDecoder().decode(concatUint8Arrays(stdoutChunks));
  const stderrText = new TextDecoder().decode(concatUint8Arrays(stderrChunks));
  if (exitCode !== 0) {
    throw new Error(stderrText.trim() || `LOSAT exited with code ${exitCode}`);
  }
  return stdoutText;
};

/**
 * @param {LosatWasiShim} wasiShim
 * @param {WebAssembly.Module} wasmModule
 * @returns {Promise<WebAssembly.Instance>}
 */
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

/**
 * @param {LosatDirectExports} exports
 * @param {Uint8Array} bytes
 */
const copyBytesToWasm = (exports, bytes) => {
  if (bytes.length === 0) return { ptr: 0, len: 0 };
  const ptr = exports.losat_web_alloc(bytes.length);
  new Uint8Array(exports.memory.buffer).set(bytes, ptr);
  return { ptr, len: bytes.length };
};

/**
 * @param {LosatDirectExports} exports
 * @param {() => number} ptrFn
 * @param {() => number} lenFn
 */
const readWasmString = (exports, ptrFn, lenFn) => {
  const ptr = ptrFn();
  const len = lenFn();
  if (!ptr || !len) return '';
  const bytes = new Uint8Array(exports.memory.buffer, ptr, len).slice();
  return new TextDecoder().decode(bytes);
};

/**
 * @param {LosatPairRun} run
 * @returns {Promise<string>}
 */
export const runLosatPairDirect = async ({
  program,
  queryFasta,
  subjectFasta,
  outfmt = '6',
  extraArgs = /** @type {string[]} */ ([]),
  wasiShim,
  wasmModule
}) => {
  const instance = await instantiateDirectLosat(wasiShim, wasmModule);
  // hasDirectLosatApi admitted this module, so its exports are the direct API.
  const exports = /** @type {LosatDirectExports} */ (instance.exports);
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
