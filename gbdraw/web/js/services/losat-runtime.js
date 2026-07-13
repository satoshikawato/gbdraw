let directInstancePromise = null;

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

export const runLosatPairDirect = async ({
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
