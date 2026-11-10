// @ts-check
import { LOSAT_SEARCH_FLAGS } from '../contracts/schema.generated.js';
import { diagnosticError } from '../utils/error-normalization.js';
/** @import { DiagnosticError } from '../utils/error-normalization.js' */

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
 * @property {string} program one of `LOSAT_PROGRAMS`
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
 * The programs the app runs with LOSAT, in every execution mode (the main thread, the pair Worker and the threaded
 * Worker): Pairwise Comparisons' blastn, tblastx and blastp, and Sequence search's blastn and tblastn (PD-OI-117).
 * @type {ReadonlySet<string>}
 */
export const LOSAT_PROGRAMS = new Set(['blastn', 'tblastx', 'blastp', 'tblastn']);

/**
 * Refuses a pair run of a program the app does not run.
 * @param {unknown} program
 * @returns {asserts program is string}
 */
export function requireLosatProgram(program) {
  if (typeof program !== 'string' || !LOSAT_PROGRAMS.has(program)) {
    throw new Error(`LOSAT program must be one of ${[...LOSAT_PROGRAMS].join(', ')}.`);
  }
}

/**
 * The context of a LOSAT run that exited with an error: the code `LOSAT_RUNTIME` with the reason `FAILED`, as Python's
 * `run_losat_search` raises it for the CLI.
 * @typedef {object} LosatFailureContext
 * @property {'FAILED'} reason
 * @property {string} [losatMessage] The first line LOSAT wrote about the failure, without its `error:` prefix (for an
 *   argument LOSAT refuses, the option and the value). Absent when LOSAT wrote nothing.
 */

/** @typedef {DiagnosticError & { context: LosatFailureContext }} LosatRunFailure */

/**
 * The failure of a LOSAT run that exited with an error, for every engine (the serial reactor, a WASI command, the
 * threaded worker) and so for every caller (pair comparisons, Sequence search, Find in other files). A pool worker posts what LOSAT
 * reported (`reportedLosatFailure`) and the page raises the failure again here.
 * @param {string} reported LOSAT's error text: its stderr, or the reactor's error buffer; empty when it wrote none
 * @returns {LosatRunFailure}
 */
export const losatRunFailure = (reported) => {
  const line = String(reported).split(/\r?\n/).map((text) => text.trim()).find(Boolean) ?? '';
  const losatMessage = line.replace(/^error:\s*/i, '');
  return /** @type {LosatRunFailure} */ (diagnosticError('LOSAT_RUNTIME', {
    reason: 'FAILED', ...(losatMessage && { losatMessage })
  }, { stage: 'losat' }));
};

/**
 * What a LOSAT worker posts of a failure so the page can raise it again: what LOSAT reported (empty when it wrote
 * nothing) for a failure of `losatRunFailure`, `null` for any other failure.
 * @param {unknown} error
 * @returns {string | null}
 */
export const reportedLosatFailure = (error) => {
  if (!(error instanceof Error) || /** @type {Partial<LosatRunFailure>} */ (error).code !== 'LOSAT_RUNTIME') return null;
  const { context } = /** @type {LosatRunFailure} */ (error);
  return context.reason === 'FAILED' ? context.losatMessage ?? '' : null;
};

/**
 * The failure of a threaded LOSAT run whose host could not run LOSAT's threads: a thread worker that could not be made
 * or prepared, a thread that did not start in time, a thread worker that failed (OV-3316). The serial engine runs no
 * threads, so Auto continues with it and Threaded names it: `LOSAT_THREADING_UNAVAILABLE` with the reason
 * `THREADS_FAILED`.
 * @param {string} detail what failed in the host, kept as the message for the console; the normalizer words the failure
 * @returns {DiagnosticError}
 */
export const losatThreadsFailure = (detail) => Object.assign(
  diagnosticError('LOSAT_THREADING_UNAVAILABLE', { reason: 'THREADS_FAILED' }, { stage: 'losat' }),
  detail ? { message: detail } : {}
);

/**
 * Whether a failure is the threaded host failure (`losatThreadsFailure`), which the threaded job worker posts as
 * `hostFailure` so the page raises it again.
 * @param {unknown} error
 * @returns {boolean}
 */
export const isLosatThreadsFailure = (error) => error instanceof Error
  && /** @type {Partial<DiagnosticError>} */ (error).code === 'LOSAT_THREADING_UNAVAILABLE'
  && /** @type {DiagnosticError} */ (error).context.reason === 'THREADS_FAILED';

/** @type {Readonly<Record<string, string>>} */
const SEARCH_FLAGS = LOSAT_SEARCH_FLAGS;

/**
 * The words of a search's options as LOSAT reads them. Generate keeps the raw-cache spelling in its args
 * (`--task`, the CLI's cache spelling too); LOSAT reads the NCBI spelling (`-task`).
 * @param {readonly unknown[] | undefined} extraArgs
 * @returns {string[]}
 */
export const losatSearchArgs = (extraArgs) => (
  Array.isArray(extraArgs) ? extraArgs.map((arg) => SEARCH_FLAGS[String(arg)] ?? String(arg)) : []
);

/**
 * The argv of the LOSAT command for one pair whose inputs are `query.fa` and `subject.fa`.
 * @param {{ program: string, outfmt?: string, extraArgs?: readonly unknown[], threads?: number }} search
 * @returns {string[]}
 */
export const losatCommandArgs = ({ program, outfmt = '6', extraArgs, threads }) => {
  const args = [
    'losat', program, '-query', 'query.fa', '-subject', 'subject.fa', '-outfmt', String(outfmt || '6'),
    ...losatSearchArgs(extraArgs)
  ];
  if (threads !== undefined && !args.some((arg) => arg === '-num_threads' || arg.startsWith('-num_threads='))) {
    args.push('-num_threads', String(threads));
  }
  return args;
};

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
  const wasi = new WASI(losatCommandArgs({ program, outfmt, extraArgs }), [], [stdin, stdout, stderr, preopen]);
  const instance = await WebAssembly.instantiate(wasmModule, {
    wasi_snapshot_preview1: wasi.wasiImport
  });
  checkCanceled();
  const exitCode = wasi.start(instance);
  checkCanceled();

  const stdoutText = new TextDecoder().decode(concatUint8Arrays(stdoutChunks));
  const stderrText = new TextDecoder().decode(concatUint8Arrays(stderrChunks));
  if (exitCode !== 0) throw losatRunFailure(stderrText);
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
  const words = [program, queryFasta, subjectFasta, outfmt, losatSearchArgs(extraArgs).join('\0')];
  /** @type {{ ptr: number, len: number }[]} */
  const allocations = [];
  let exitCode;
  try {
    words.forEach((word) => allocations.push(copyBytesToWasm(exports, encoder.encode(String(word)))));
    exitCode = exports.losat_web_run_pair(...allocations.flatMap(({ ptr, len }) => [ptr, len]));
  } catch (error) {
    // A trap (LOSAT docs/web/abi_v1.md §9): the instance is never called again, not even to free or clear;
    // the next pair instantiates the module again.
    directInstancePromise = null;
    throw error;
  }

  try {
    if (exitCode !== 0) {
      const detail = readWasmString(
        exports,
        exports.losat_web_error_ptr,
        exports.losat_web_error_len
      );
      throw losatRunFailure(detail);
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
