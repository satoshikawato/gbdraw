const DEFAULT_WASM_PATH = './wasm/losat/losat.wasm';

let wasiShimPromise = null;
let wasmModulePromise = null;

const loadWasiShim = async () => {
  if (!wasiShimPromise) {
    wasiShimPromise = import(
      'https://unpkg.com/@bjorn3/browser_wasi_shim@0.4.2/dist/index.js'
    );
  }
  return wasiShimPromise;
};

const loadLosatModule = async (wasmPath = DEFAULT_WASM_PATH) => {
  if (!wasmModulePromise) {
    wasmModulePromise = (async () => {
      const response = await fetch(wasmPath);
      if (!response.ok) {
        throw new Error(`LOSAT wasm not found at ${wasmPath}`);
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
