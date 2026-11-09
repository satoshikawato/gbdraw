// The LOSAT jobs of a Find in other files search, by File (`runLosatJobsByFile` in services/losat.js, PD-OI-123): one
// job per File, one after another, the status "File k of N" before each; a File whose job fails keeps its reason and
// the next File still runs; Cancel stops the run and no File after it runs.
import assert from 'node:assert/strict';
import test from 'node:test';

const { runLosatJobsByFile } = await import('../../gbdraw/web/js/services/losat.js');

const jobs = ['a', 'b', 'c'].map((fileKey, pairIndex) => ({ pairIndex, fileKey, program: 'tblastn' }));

test('each File runs in turn, the status names it, and a failed File keeps its reason while the next still runs', async () => {
  const controller = new AbortController();
  const status = [];
  const ran = [];
  const results = await runLosatJobsByFile(jobs, async (job) => {
    ran.push([job.fileKey, status.at(-1)]);
    if (job.fileKey === 'b') throw new Error('LOSAT exited with status 1');
    return [{ ...job, text: `rows of ${job.fileKey}` }];
  }, { signal: controller.signal, status: (text) => status.push(text) });
  assert.deepEqual(ran, [
    ['a', 'Running LOSAT: File 1 of 3'], ['b', 'Running LOSAT: File 2 of 3'], ['c', 'Running LOSAT: File 3 of 3']
  ]);
  assert.deepEqual(results.map(({ fileKey, text, error }) => [fileKey, text, typeof error]), [
    ['a', 'rows of a', 'undefined'], ['b', '', 'string'], ['c', 'rows of c', 'undefined']
  ]);
  // The reason in the user-facing wording of the error owner (`normalizeCaughtError`, R6).
  assert.ok(results[1].error.length > 0);
});

test('Cancel during a File stops the run: the abort is thrown and no later File runs', async () => {
  const controller = new AbortController();
  const ran = [];
  const canceled = new Error('canceled');
  await assert.rejects(runLosatJobsByFile(jobs, async (job) => {
    ran.push(job.fileKey);
    if (job.fileKey === 'b') {
      controller.abort(canceled);
      throw canceled;
    }
    return [{ ...job, text: 'rows' }];
  }, { signal: controller.signal, status: () => {} }), (error) => error === canceled);
  assert.deepEqual(ran, ['a', 'b']);
});
