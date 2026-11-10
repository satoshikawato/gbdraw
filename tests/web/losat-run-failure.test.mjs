// A LOSAT run that exits with an error (an argument LOSAT refuses, for example) reaches the user as the LOSAT
// failure Python raises for the CLI (`LOSAT_RUNTIME`, reason `FAILED`) with LOSAT's own first error line, not as an
// unrecognized failure (OV-3304). The one owner is `losatRunFailure` in services/losat-runtime.js; a pool worker
// posts what LOSAT reported and the page raises the same failure again.
import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import test from 'node:test';

const web = new URL('../../gbdraw/web/', import.meta.url);
const runtime = await import(new URL('js/services/losat-runtime.js', web).href);
const { normalizeUserFacingError } = await import(new URL('js/utils/error-normalization.js', web).href);
const { errorText } = await import(new URL('js/utils/error-text.js', web).href);
const wasiShimUrl = new URL('vendor/browser_wasi_shim/dist/index.js', web).href;

const REFUSED = ['--searchsp', '99999999999999999999'];
const pair = { program: 'blastn', queryFasta: '>q\nACGTACGTACGTACGTACGT\n', subjectFasta: '>s\nACGTACGTACGTACGTACGT\n' };

/** @param {unknown} error */
const isRefusal = (error) => error instanceof Error && error.code === 'LOSAT_RUNTIME' && error.stage === 'losat'
  && error.context.reason === 'FAILED' && /searchsp/.test(error.context.losatMessage);

test('a LOSAT failure names the search and shows what LOSAT reported, as text', () => {
  const failure = runtime.losatRunFailure("error: invalid value '9' for '-searchsp <INT_VALUE>': too large\n\nFor more information, try '--help'.\n");
  assert.deepEqual({ code: failure.code, stage: failure.stage, context: failure.context }, {
    code: 'LOSAT_RUNTIME', stage: 'losat',
    context: { reason: 'FAILED', losatMessage: "invalid value '9' for '-searchsp <INT_VALUE>': too large" }
  });
  const model = normalizeUserFacingError(failure, { operation: 'generate' });
  assert.equal(model.code, 'LOSAT_RUNTIME');
  assert.equal(model.stage, 'losat');
  // OV-3317: LOSAT said why it stopped, so Generate with the same settings stops the same way: the one action shows
  // the comparison settings.
  assert.deepEqual(model.actions, ['edit-comparison']);
  // The text says what the action leads to: Generate with the same settings stops the same way, and what LOSAT
  // reported names what to change.
  assert.equal(model.summary, 'LOSAT stopped the comparison search. Generating again with the same settings and inputs stops the same way; change what LOSAT names below, then Generate again.'
    + " LOSAT reported: invalid value '9' for '\u2011searchsp <INT_VALUE>': too large");
  // Sequence search runs from its own Search button: its text names that button, not Generate.
  assert.equal(normalizeUserFacingError(failure, { operation: 'searchHitAnnotationSet' }).summary,
    'The search could not run or returned unusable output. Check the search settings, then click Search again.'
    + ' The search runtime exited with an error; check its message and inputs.'
    + " LOSAT reported: invalid value '9' for '\u2011searchsp <INT_VALUE>': too large");
  // LOSAT's words are a value of the phrase, never a catalog key: a translation renders them as written, except that
  // a hyphen starting an option name does not end a line; the Diagnostics keep LOSAT's own text.
  assert.equal(errorText(model.summaryParts).endsWith("LOSAT reported: invalid value '9' for '\u2011searchsp <INT_VALUE>': too large"), true);
  assert.match(model.details[0].text, /losatMessage: invalid value '9' for '-searchsp <INT_VALUE>': too large/);
  assert.match(normalizeUserFacingError(runtime.losatRunFailure('unsupported argument: --searchsp, - x-y')).summary,
    /LOSAT reported: unsupported argument: \u2011\u2011searchsp, \u2011 x-y$/);
  // One line of at most 200 characters, control characters folded.
  const long = normalizeUserFacingError({ code: 'LOSAT_RUNTIME', stage: 'losat', context: { reason: 'FAILED', losatMessage: `a\u0000b\tc ${'x'.repeat(300)}` } });
  assert.equal(long.context.losatMessage, `a b c ${'x'.repeat(194)}`);
  // LOSAT wrote nothing: the failure keeps its code and reason, and no LOSAT line.
  const silent = runtime.losatRunFailure('');
  assert.deepEqual(silent.context, { reason: 'FAILED' });
  assert.doesNotMatch(normalizeUserFacingError(silent).summary, /LOSAT reported/);
  // Without LOSAT's words the cause is not known to repeat, so Retry stays.
  assert.deepEqual(normalizeUserFacingError(silent, { operation: 'generate' }).actions, ['edit-comparison', 'retry']);
  assert.equal(runtime.reportedLosatFailure(silent), '');
  assert.equal(runtime.reportedLosatFailure(failure), failure.context.losatMessage);
  assert.equal(runtime.reportedLosatFailure(new Error('LOSAT_RUNTIME')), null);
});

test('the serial reactor raises the LOSAT failure for an argument LOSAT refuses', async () => {
  const wasiShim = await import(wasiShimUrl);
  const wasmModule = await WebAssembly.compile(await readFile(new URL('wasm/losat/losat.wasm', web)));
  await assert.rejects(runtime.runLosatPairDirect({ ...pair, extraArgs: REFUSED, wasiShim, wasmModule }), isRefusal);
});

test('a WASI command that exits with an error raises the LOSAT failure with its first stderr line', async () => {
  const encoder = new TextEncoder();
  const wasiShim = {
    WASI: class {
      constructor(_args, _env, fds) { this.fds = fds; this.wasiImport = {}; }
      start() {
        this.fds[2].write(encoder.encode("\nerror: invalid value '99999999999999999999' for '-searchsp <INT_VALUE>'\n"));
        return 2;
      }
    },
    File: class {},
    OpenFile: class {},
    PreopenDirectory: class {},
    ConsoleStdout: class { constructor(write) { this.write = write; } }
  };
  const wasmModule = await WebAssembly.compile(new Uint8Array([0x00, 0x61, 0x73, 0x6d, 0x01, 0x00, 0x00, 0x00]));
  await assert.rejects(runtime.runLosatPairWasi({ ...pair, extraArgs: REFUSED, wasiShim, wasmModule }), isRefusal);
});

test('the serial pool worker posts what LOSAT reported', async (t) => {
  const previous = Object.getOwnPropertyDescriptor(globalThis, 'self');
  t.after(() => { if (previous) Object.defineProperty(globalThis, 'self', previous); else delete globalThis.self; });
  const posted = [];
  globalThis.self = { postMessage: (message) => posted.push(message) };
  await import(`${new URL('js/workers/losat-worker.js', web).href}?instance=${Date.now()}`);
  const wasmModule = await WebAssembly.compile(await readFile(new URL('wasm/losat/losat.wasm', web)));
  await globalThis.self.onmessage({ data: { type: 'init', id: 'i', wasmModule, wasiShimUrl } });
  await globalThis.self.onmessage({
    data: {
      type: 'run', id: 'r', job: { program: 'blastn', querySequenceKey: 'q', subjectSequenceKey: 's', extraArgs: REFUSED },
      sequences: [{ key: 'q', fasta: pair.queryFasta }, { key: 's', fasta: pair.subjectFasta }]
    }
  });
  assert.deepEqual(posted.map(({ id, ok }) => ({ id, ok })), [{ id: 'i', ok: true }, { id: 'r', ok: false }]);
  assert.match(posted[1].losatReported, /searchsp/);
  assert.equal(isRefusal(runtime.losatRunFailure(posted[1].losatReported)), true);
});

// OV-3316: the threaded engine's host could not run LOSAT's threads (a thread that did not start in time, a thread
// worker that could not be made). The failure says so and names the engines that can run, never an unknown failure.
test('a threaded host failure names the threaded engine and the engines that can run', () => {
  const failure = runtime.losatThreadsFailure('timed out waiting for LOSAT WASI thread 1 to start after 30s');
  assert.deepEqual({ code: failure.code, stage: failure.stage, context: failure.context },
    { code: 'LOSAT_THREADING_UNAVAILABLE', stage: 'losat', context: { reason: 'THREADS_FAILED' } });
  assert.match(failure.message, /thread 1 to start/);
  assert.equal(runtime.isLosatThreadsFailure(failure), true);
  assert.equal(runtime.isLosatThreadsFailure(runtime.losatRunFailure('error: x')), false);
  assert.equal(runtime.reportedLosatFailure(failure), null);
  const model = normalizeUserFacingError(failure, { operation: 'generate' });
  assert.equal(model.summary, 'Threaded LOSAT execution is unavailable in this browser environment. Select Serial or Auto execution, then Generate again.'
    + ' The threaded LOSAT engine could not start or keep its threads running in this browser.');
  assert.deepEqual(model.actions, ['edit-comparison', 'retry']);
  assert.doesNotMatch(model.summary, /thread 1/);
  // Sequence search and Find in other files run from their own Search button, not from Generate.
  for (const operation of ['searchHitAnnotationSet', 'findInOtherFiles']) {
    assert.equal(normalizeUserFacingError(failure, { operation }).summary,
      'Threaded LOSAT execution is unavailable in this browser environment. Select Serial or Auto execution, then click Search again.'
      + ' The threaded LOSAT engine could not start or keep its threads running in this browser.');
  }
});

// OV-3327: a trap of LOSAT's code (a LOSAT panic, a WebAssembly "unreachable") is the LOSAT failure with the reason
// TRAP on every engine, never an unknown failure. The same search traps the same way, so the panel offers the
// comparison settings, not Retry.
const isTrap = (error) => error instanceof Error && error.code === 'LOSAT_RUNTIME' && error.stage === 'losat'
  && error.context.reason === 'TRAP';
// LOSAT panics on a genetic code it does not know (OV-3326).
const TRAPPING = { ...pair, program: 'tblastx', extraArgs: ['--query-gencode', '99'] };

test('a LOSAT trap names the comparison search and leads to the comparison settings', () => {
  const failure = runtime.raisedLosatFailure(new WebAssembly.RuntimeError('unreachable'));
  assert.equal(isTrap(failure), true);
  assert.equal(failure.message, 'unreachable');
  assert.equal(runtime.isLosatTrapFailure(failure), true);
  assert.equal(runtime.reportedLosatFailure(failure), null);
  const plain = new Error('x');
  assert.equal(runtime.raisedLosatFailure(plain), plain);
  const model = normalizeUserFacingError(failure, { operation: 'generate' });
  assert.equal(model.summary, 'LOSAT stopped unexpectedly during the comparison search. Generating again with the same settings and inputs stops the same way; change the comparison settings or inputs, then Generate again.');
  assert.deepEqual(model.actions, ['edit-comparison']);
  assert.equal(normalizeUserFacingError(failure, { operation: 'findInOtherFiles' }).summary,
    'The search could not run or returned unusable output. Check the search settings, then click Search again. LOSAT stopped unexpectedly.');
});

test('the serial reactor raises the LOSAT trap for a genetic code LOSAT does not know', async () => {
  const wasiShim = await import(wasiShimUrl);
  const wasmModule = await WebAssembly.compile(await readFile(new URL('wasm/losat/losat.wasm', web)));
  await assert.rejects(runtime.runLosatPairDirect({ ...TRAPPING, wasiShim, wasmModule }), isTrap);
});

test('a WASI command that traps raises the LOSAT trap', async () => {
  const wasiShim = {
    WASI: class { constructor() { this.wasiImport = {}; } start() { throw new WebAssembly.RuntimeError('unreachable'); } },
    File: class {}, OpenFile: class {}, PreopenDirectory: class {}, ConsoleStdout: class {}
  };
  const wasmModule = await WebAssembly.compile(new Uint8Array([0x00, 0x61, 0x73, 0x6d, 0x01, 0x00, 0x00, 0x00]));
  await assert.rejects(runtime.runLosatPairWasi({ ...TRAPPING, wasiShim, wasmModule }), isTrap);
});

test('the serial pool worker posts a LOSAT trap as one', async (t) => {
  const previous = Object.getOwnPropertyDescriptor(globalThis, 'self');
  t.after(() => { if (previous) Object.defineProperty(globalThis, 'self', previous); else delete globalThis.self; });
  const posted = [];
  globalThis.self = { postMessage: (message) => posted.push(message) };
  await import(`${new URL('js/workers/losat-worker.js', web).href}?instance=trap-${Date.now()}`);
  const wasmModule = await WebAssembly.compile(await readFile(new URL('wasm/losat/losat.wasm', web)));
  await globalThis.self.onmessage({ data: { type: 'init', id: 'i', wasmModule, wasiShimUrl } });
  await globalThis.self.onmessage({
    data: {
      type: 'run', id: 'r', job: { program: 'tblastx', querySequenceKey: 'q', subjectSequenceKey: 's', extraArgs: TRAPPING.extraArgs },
      sequences: [{ key: 'q', fasta: pair.queryFasta }, { key: 's', fasta: pair.subjectFasta }]
    }
  });
  assert.deepEqual(posted.map(({ id, ok, losatReported, trap }) => ({ id, ok, losatReported, trap })),
    [{ id: 'i', ok: true, losatReported: undefined, trap: undefined }, { id: 'r', ok: false, losatReported: null, trap: true }]);
});
