const DEFAULT_WASM_URL = '../../wasm/losat/losat.wasm';

let wasiShimPromise = null;
let wasmModulePromise = null;
let wasmModuleUrl = null;

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
  if (!program || (program !== 'blastn' && program !== 'tblastx')) {
    throw new Error('LOSAT program must be blastn or tblastx.');
  }
  if (!queryFasta || !subjectFasta) {
    throw new Error('LOSAT requires both query and subject FASTA content.');
  }
  if (!wasiShimUrl) {
    throw new Error('LOSAT worker requires a browser_wasi_shim URL.');
  }

  const [{ WASI, File, OpenFile, PreopenDirectory, ConsoleStdout }, compiledModule] = await Promise.all([
    loadWasiShim(wasiShimUrl),
    getLosatModule({ wasmModule, wasmUrl })
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
