const DEFAULT_WASM_URL = '../../wasm/losat/losat.wasm';
const SUPPORTED_PROGRAMS = new Set(['blastn', 'tblastx', 'blastp']);

let wasiShimPromise = null;
let wasmModulePromise = null;
let wasmModuleUrl = null;
let directInstancePromise = null;

const loadWasiShim = async (wasiShimUrl) => {
  if (!wasiShimPromise) {
    wasiShimPromise = import(wasiShimUrl).catch((error) => {
      throw new Error(
        `Missing packaged asset: browser_wasi_shim module at ${wasiShimUrl}. ${error.message}`
      );
    });
  }
  return wasiShimPromise;
};

const loadLosatModule = async (wasmUrl = DEFAULT_WASM_URL) => {
  if (!wasmModulePromise || wasmModuleUrl !== wasmUrl) {
    wasmModuleUrl = wasmUrl;
    wasmModulePromise = (async () => {
      const response = await fetch(wasmUrl, { cache: 'no-store' });
      if (!response.ok) {
        throw new Error(`Missing packaged asset: LOSAT wasm not found at ${wasmUrl}`);
      }
      const bytes = await response.arrayBuffer();
      return WebAssembly.compile(bytes);
    })();
  }
  return wasmModulePromise;
};

const getLosatModule = async ({ wasmModule, wasmUrl } = {}) => {
  if (wasmModule instanceof WebAssembly.Module) return wasmModule;
  return loadLosatModule(wasmUrl);
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

const runLosatPair = async ({
  program,
  queryFasta,
  subjectFasta,
  outfmt = '6',
  extraArgs = [],
  wasmModule,
  wasmUrl,
  wasiShimUrl
} = {}) => {
  if (!program || !SUPPORTED_PROGRAMS.has(program)) {
    throw new Error('LOSAT program must be blastn, tblastx, or blastp.');
  }
  if (!queryFasta || !subjectFasta) {
    throw new Error('LOSAT requires both query and subject FASTA content.');
  }
  if (!wasiShimUrl) {
    throw new Error('LOSAT worker requires a browser_wasi_shim URL.');
  }

  const [wasiShim, compiledModule] = await Promise.all([
    loadWasiShim(wasiShimUrl),
    getLosatModule({ wasmModule, wasmUrl })
  ]);
  const { WASI, File, OpenFile, PreopenDirectory, ConsoleStdout } = wasiShim;

  if (hasDirectLosatApi(compiledModule)) {
    return runLosatPairDirect({
      program,
      queryFasta,
      subjectFasta,
      outfmt,
      extraArgs,
      wasiShim,
      wasmModule: compiledModule
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
  const instance = await WebAssembly.instantiate(compiledModule, {
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

self.onmessage = async (event) => {
  const { id, job } = event.data || {};
  try {
    const text = await runLosatPair(job);
    self.postMessage({ id, ok: true, text });
  } catch (error) {
    self.postMessage({
      id,
      ok: false,
      error: error?.message ? String(error.message) : String(error || 'Unknown LOSAT worker error')
    });
  }
};
