const { test, expect } = require('@playwright/test');
const { evaluateWithRetainedPromise, routeCrossOriginIsolated } = require('./helpers/app-lifecycle.cjs');

// Stand-in for losat-threaded.wasm: _start spawns one WASI thread, then waits
// on a futex that only that thread could notify. The thread traps at once, as a
// LOSAT thread does when Chromium rejects memory.copy after another thread grew
// the shared memory. Without trap reporting, the job never settles.
const leb = (value) => {
  const bytes = [];
  do {
    let byte = value & 0x7f;
    value >>>= 7;
    if (value) byte |= 0x80;
    bytes.push(byte);
  } while (value);
  return bytes;
};
const name = (text) => [...leb(text.length), ...Buffer.from(text)];
const section = (id, entries) => {
  const body = [...leb(entries.length), ...entries.flat()];
  return [id, ...leb(body.length), ...body];
};
const I32 = 0x7f;
const threadTrapJobModule = Buffer.from([
  0x00, 0x61, 0x73, 0x6d, 0x01, 0x00, 0x00, 0x00,
  ...section(1, [[0x60, 1, I32, 1, I32], [0x60, 0, 0], [0x60, 2, I32, I32, 0]]),
  ...section(2, [
    [...name('wasi'), ...name('thread-spawn'), 0x00, 0],
    [...name('env'), ...name('memory'), 0x02, 0x03, ...leb(1), ...leb(16384)]
  ]),
  ...section(3, [[1], [2]]),
  ...section(7, [[...name('_start'), 0x00, 1], [...name('wasi_thread_start'), 0x00, 2]]),
  ...section(10, [
    [...leb(18), 0x00, 0x41, 0x00, 0x10, 0x00, 0x1a,
      0x41, 0x00, 0x41, 0x00, 0x42, 0x7f, 0xfe, 0x01, 0x02, 0x00, 0x1a, 0x0b],
    [...leb(3), 0x00, 0x00, 0x0b]
  ])
]);

test('a trapping WASI thread fails its threaded LOSAT job after bounded retries instead of stalling it', async ({ page }) => {
  test.setTimeout(60000);
  await routeCrossOriginIsolated(page.context(), {
    '/gbdraw/web/losat-thread-fault.html': { contentType: 'text/html', body: '<!doctype html><title>LOSAT thread fault</title>' },
    '/gbdraw/web/wasm/losat/thread-trap.wasm': { contentType: 'application/wasm', body: threadTrapJobModule }
  });
  await page.addInitScript(() => {
    Object.defineProperty(navigator, 'hardwareConcurrency', { get: () => 4 });
  });
  const threadWorkers = [];
  page.on('worker', (worker) => {
    if (!worker.url().includes('losat-wasi-thread-worker.js')) return;
    const entry = { closed: false };
    threadWorkers.push(entry);
    worker.on('close', () => { entry.closed = true; });
  });

  await page.goto('/gbdraw/web/losat-thread-fault.html');
  expect(await page.evaluate(() => crossOriginIsolated)).toBe(true);
  const outcome = await evaluateWithRetainedPromise(page, async () => {
    const { runLosatPairsParallel } = await import('/gbdraw/web/js/services/losat.js');
    const run = runLosatPairsParallel(
      [{ pairIndex: 0, program: 'blastp', querySequenceKey: 'q', subjectSequenceKey: 's' }],
      {
        executionMode: 'threaded', totalThreadBudget: 2, concurrency: 1, threadsPerJob: 2,
        threadedWasmPath: './wasm/losat/thread-trap.wasm',
        sequences: { q: '>q\nMKV', s: '>s\nMKV' }
      }
    ).then(
      () => ({ status: 'resolved' }),
      (error) => ({ status: 'rejected', message: String(error?.message || error) })
    );
    // Guard for the pre-fix behavior only: the job otherwise waits forever.
    const stalled = new Promise((resolve) => setTimeout(() => resolve({ status: 'stalled' }), 20000));
    return Promise.race([run, stalled]);
  });

  expect(outcome.status).toBe('rejected');
  expect(outcome.message).toMatch(/^LOSAT pair #1: LOSAT thread 1 failed: .*unreachable/);
  // The synthetic thread always traps: one attempt and two retries.
  expect(threadWorkers.length).toBe(3);
  // Ending each job also ends its WASI thread worker.
  await expect.poll(() => threadWorkers.every((entry) => entry.closed)).toBe(true);
});

// The shipped threaded file runs a pair in the browser: the job's shared memory meets the
// memory the file imports (its recorded page limits), and LOSAT reads the search's words.
test('the shipped threaded LOSAT wasm runs a blastp pair with two threads', async ({ page }) => {
  test.setTimeout(120000);
  await routeCrossOriginIsolated(page.context(), {
    '/gbdraw/web/losat-threaded-run.html': { contentType: 'text/html', body: '<!doctype html><title>LOSAT threaded run</title>' }
  });
  await page.goto('/gbdraw/web/losat-threaded-run.html');
  expect(await page.evaluate(() => crossOriginIsolated)).toBe(true);
  const protein = 'MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQ';
  const outcome = await evaluateWithRetainedPromise(page, async (sequence) => {
    const { runLosatPairsParallel } = await import('/gbdraw/web/js/services/losat.js');
    return runLosatPairsParallel(
      [{ pairIndex: 0, cacheKey: 'pair', program: 'blastp', querySequenceKey: 'q', subjectSequenceKey: 's',
        extraArgs: ['--max-hsps-per-subject', '1', '--max-target-seqs', '5'] }],
      {
        executionMode: 'threaded', totalThreadBudget: 2, concurrency: 1, threadsPerJob: 2,
        sequences: { q: `>q1\n${sequence}\n`, s: `>s1\n${sequence}\n` }
      }
    ).then(
      (results) => ({ status: 'resolved', text: results.map((result) => result.text).join('') }),
      (error) => ({ status: 'rejected', message: String(error?.message || error) })
    );
  }, protein);
  expect(outcome).toMatchObject({ status: 'resolved' });
  expect(outcome.text.trim().split('\t').slice(0, 4)).toEqual(['q1', 's1', '100.000', String(protein.length)]);
});
