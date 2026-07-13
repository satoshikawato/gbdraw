import {
  hasDirectLosatApi,
  runLosatPairDirect,
  runLosatPairWasi
} from '../services/losat-runtime.js';

const DEFAULT_WASM_URL = '../../wasm/losat/losat.wasm';
const SUPPORTED_PROGRAMS = new Set(['blastn', 'tblastx', 'blastp']);

let wasiShimPromise = null;
let wasmModulePromise = null;
let wasmModuleUrl = null;
let runtime = null;
let sequenceStore = new Map();

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

  return runLosatPairWasi({
    program,
    queryFasta,
    subjectFasta,
    outfmt,
    extraArgs,
    wasiShim,
    wasmModule: compiledModule
  });
};

const initializeWorkerRuntime = async ({
  wasmModule,
  wasmUrl = DEFAULT_WASM_URL,
  wasiShimUrl,
  sequences = []
} = {}) => {
  if (!wasiShimUrl) {
    throw new Error('LOSAT worker requires a browser_wasi_shim URL.');
  }
  const entries = Array.isArray(sequences) ? sequences : [];
  sequenceStore = new Map(
    entries
      .filter((entry) => entry && entry.key !== undefined)
      .map((entry) => [String(entry.key), String(entry.fasta || '')])
  );
  runtime = {
    wasmModule,
    wasmUrl,
    wasiShimUrl
  };
  await Promise.all([
    loadWasiShim(wasiShimUrl),
    getLosatModule({ wasmModule, wasmUrl })
  ]);
};

const mergeSequencePayload = (sequences = []) => {
  if (!Array.isArray(sequences)) return;
  sequences.forEach((entry) => {
    if (!entry || entry.key === undefined) return;
    sequenceStore.set(String(entry.key), String(entry.fasta || ''));
  });
};

const resolveSequence = (key, label) => {
  const normalizedKey = String(key || '');
  if (!normalizedKey || !sequenceStore.has(normalizedKey)) {
    throw new Error(`Missing LOSAT ${label} sequence for key '${normalizedKey || '(blank)'}'.`);
  }
  const fasta = sequenceStore.get(normalizedKey);
  if (!fasta) {
    throw new Error(`LOSAT ${label} sequence for key '${normalizedKey}' is empty.`);
  }
  return fasta;
};

self.onmessage = async (event) => {
  const { id, type, job, sequences } = event.data || {};
  try {
    if (type === 'init') {
      await initializeWorkerRuntime(event.data || {});
      self.postMessage({ id, ok: true, type: 'init' });
      return;
    }
    if (type !== 'run') {
      throw new Error(`Unsupported LOSAT worker message type '${type || '(blank)'}'.`);
    }
    if (!runtime) {
      throw new Error('LOSAT worker has not been initialized.');
    }
    mergeSequencePayload(sequences);
    const text = await runLosatPair({
      ...job,
      queryFasta: resolveSequence(job?.querySequenceKey, 'query'),
      subjectFasta: resolveSequence(job?.subjectSequenceKey, 'subject'),
      wasmModule: runtime.wasmModule,
      wasmUrl: runtime.wasmUrl,
      wasiShimUrl: runtime.wasiShimUrl
    });
    self.postMessage({ id, ok: true, type: 'run', text });
  } catch (error) {
    self.postMessage({
      id,
      type: type || 'run',
      ok: false,
      error: error?.message ? String(error.message) : String(error || 'Unknown LOSAT worker error')
    });
  }
};
