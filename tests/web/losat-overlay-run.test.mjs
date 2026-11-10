// A run under the full-window overlay (services/losat.js): its Cancel (`runUntilCanceled`: a step's reply or Cancel,
// whichever comes first; a reply after Cancel is discarded; no status after Cancel), and the LOSAT jobs of a Find in
// other files search by File (`runLosatJobsByFile`, PD-OI-123): one job per File, one after another, the status
// "File k of N" before each; a File whose job fails keeps its reason and the next File still runs; Cancel stops the
// run and no File after it runs.
import assert from 'node:assert/strict';
import test from 'node:test';

const { runLosatJobsByFile, runUntilCanceled } = await import('../../gbdraw/web/js/services/losat.js');
const { losatRunFailure } = await import('../../gbdraw/web/js/services/losat-runtime.js');

/** A step whose reply the test gives when it chooses. */
const held = () => {
  let reply = (/** @type {unknown} */ _value) => {};
  const step = new Promise((resolve) => { reply = resolve; });
  return { step, reply };
};

test('the steps of a run give their replies and the status line follows them', async () => {
  const written = [];
  const run = await runUntilCanceled(new AbortController().signal, (text) => written.push(text), async ({ until, status }) => {
    status('Reading');
    const files = await until(Promise.resolve(['a', 'b']));
    status('Planning');
    return { files, plan: await until(Promise.resolve('plan')) };
  });
  assert.deepEqual(run, { status: 'ok', value: { files: ['a', 'b'], plan: 'plan' } });
  assert.deepEqual(written, ['Reading', 'Planning']);
});

test('Cancel during a step ends the run at once; the step reply that comes later is discarded and nothing after it runs', async () => {
  const controller = new AbortController();
  const written = [];
  const { step, reply } = held();
  const after = [];
  const running = runUntilCanceled(controller.signal, (text) => written.push(text), async ({ until, status }) => {
    status('Reading');
    const files = await until(step);
    after.push(files);
    status('Planning');
    return files;
  });
  controller.abort(new Error('canceled'));
  // The run ends before the step replies.
  assert.deepEqual(await running, { status: 'canceled' });
  reply(['late']);
  await new Promise((resolve) => setTimeout(resolve, 0));
  assert.deepEqual(after, []);
  assert.deepEqual(written, ['Reading']);
});

test('a run that Cancel reaches after its last step, or whose task then fails, ends canceled', async () => {
  const controller = new AbortController();
  const done = await runUntilCanceled(controller.signal, () => {}, async ({ until, status }) => {
    const value = await until(Promise.resolve('reply'));
    controller.abort(new Error('canceled'));
    status('not written');
    return value;
  });
  assert.deepEqual(done, { status: 'canceled' });
  const aborted = new AbortController();
  aborted.abort(new Error('canceled'));
  const failed = await runUntilCanceled(aborted.signal, () => {}, async ({ until }) => until(Promise.resolve('x')));
  assert.deepEqual(failed, { status: 'canceled' });
});

test('an error of the task that is not Cancel reaches the caller', async () => {
  const boom = new Error('helper failed');
  await assert.rejects(
    runUntilCanceled(new AbortController().signal, () => {}, async ({ until }) => until(Promise.reject(boom))),
    (error) => error === boom
  );
});

const jobs = ['a', 'b', 'c'].map((fileKey, pairIndex) => ({ pairIndex, fileKey, program: 'tblastn' }));

test('each File runs in turn, the status names it, and a failed File keeps its reason while the next still runs', async () => {
  const controller = new AbortController();
  const status = [];
  const ran = [];
  const results = await runLosatJobsByFile(jobs, async (job) => {
    ran.push([job.fileKey, status.at(-1)]);
    if (job.fileKey === 'b') throw new Error('LOSAT exited with status 1');
    return [{ ...job, text: `rows of ${job.fileKey}` }];
  }, { signal: controller.signal, status: (text) => status.push(text) }, 'findInOtherFiles');
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
  }, { signal: controller.signal, status: () => {} }, 'findInOtherFiles'), (error) => error === canceled);
  assert.deepEqual(ran, ['a', 'b']);
});

test('a File whose LOSAT exits with an error says what LOSAT reported and names the Search, not Generate', async () => {
  const results = await runLosatJobsByFile(jobs.slice(0, 1), async () => {
    throw losatRunFailure("error: invalid value '9' for '-evalue <REAL>'\n");
  }, { signal: new AbortController().signal, status: () => {} }, 'findInOtherFiles');
  assert.match(results[0].error, /^The search could not run or returned unusable output\. Check the search settings, then click Search again\./);
  assert.match(results[0].error, /LOSAT reported: invalid value '9' for '\u2011evalue <REAL>'$/);
  assert.doesNotMatch(results[0].error, /Generate/);
});
