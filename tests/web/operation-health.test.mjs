import assert from 'node:assert/strict';
import { createRequire } from 'node:module';
import test from 'node:test';

// G-G(3) (Web GUI audit 2026-09-30): the shared Generate and Save helpers fail
// on an UNKNOWN diagnostic, an uncaught page error, or NaN/undefined in Run
// Info, and accept an opt-out only when it names an audit ID.
const require = createRequire(import.meta.url);
const { assertOperationHealth, assertOperationHealthStrict, takeOperationHealthProblems } = require('./helpers/app-lifecycle.cjs');

const fakePage = (app) => {
  const listeners = new Map();
  return {
    app,
    on: (type, listener) => listeners.set(type, listener),
    emit: (type, value) => listeners.get(type)?.(value),
    evaluate: async (callback) => {
      globalThis.window = { __GBDRAW_APP__: app };
      try {
        return callback();
      } finally {
        delete globalThis.window;
      }
    }
  };
};

test('a clean operation passes and each health failure is reported', async () => {
  const page = fakePage({ errorLog: null, lastRunInfo: { args: ['-w', '1000'], windowSize: 1000 } });
  assert.deepEqual(await assertOperationHealth(page), { errorCode: null, pageErrors: [], invalidRunInfo: [] });

  page.app.errorLog = { code: 'UNKNOWN', operation: 'generate', stage: 'request-validation', summary: 'failed' };
  await assert.rejects(assertOperationHealth(page), /must not report an UNKNOWN diagnostic/);
  page.app.errorLog = null;

  page.app.lastRunInfo = { args: ['-w', 'NaN'], stepSize: Number.NaN, nested: [{ label: 'undefined' }] };
  await assert.rejects(assertOperationHealth(page), (error) => {
    assert.match(error.message, /must not record NaN or undefined in Run Info/);
    assert.match(error.message, /lastRunInfo\.args\[1\]/);
    assert.match(error.message, /lastRunInfo\.stepSize=NaN/);
    assert.match(error.message, /lastRunInfo\.nested\[0\]\.label/);
    return true;
  });
  page.app.lastRunInfo = null;

  page.emit('pageerror', new Error('boom'));
  await assert.rejects(assertOperationHealth(page), /must not raise an uncaught page error/);
});

test('an opt-out names the audit ID of a known defect', async () => {
  const page = fakePage({
    errorLog: { code: 'UNKNOWN', operation: 'save-session', stage: 'save', summary: 'failed' },
    lastRunInfo: { windowSize: Number.NaN }
  });
  const health = await assertOperationHealth(page, { allowUnknown: 'IN-08', allowInvalidRunInfo: 'TR-04' });
  assert.equal(health.errorCode, 'UNKNOWN');
  assert.deepEqual(health.invalidRunInfo, ['lastRunInfo.windowSize=NaN']);
  for (const value of [true, 'unknown', 'IN08']) {
    await assert.rejects(assertOperationHealth(page, { allowUnknown: value }), /must name the audit ID/);
  }
});

test('an alert left by an earlier operation is not attributed to Save', async () => {
  const errorLog = { code: 'UNKNOWN', operation: 'generate', stage: 'render', summary: 'earlier' };
  const page = fakePage({ errorLog, lastRunInfo: null });
  const errorSignatureBefore = JSON.stringify([errorLog.code, errorLog.operation, errorLog.stage, errorLog.summary]);
  assert.equal((await assertOperationHealth(page, { errorSignatureBefore })).errorCode, null);
  page.app.errorLog = { ...errorLog, operation: 'save-session', summary: 'new' };
  await assert.rejects(assertOperationHealth(page, { errorSignatureBefore }), /UNKNOWN diagnostic/);
});

// The check the random walk runs after every step and the audit sweeps after
// every Load and operation: no opt-out, VALIDATION_UNCLASSIFIED fails as
// UNKNOWN does, console errors fail, and each call reports only the errors
// raised since the call before it.
test('the strict check reports every problem once and has no opt-out', async () => {
  const page = fakePage({ errorLog: null, lastRunInfo: null });
  assert.deepEqual(await takeOperationHealthProblems(page), []);
  await assertOperationHealthStrict(page, { operation: 'step 1' });

  page.app.errorLog = { code: 'VALIDATION_UNCLASSIFIED', operation: 'generate', stage: 'request-validation', summary: 'Input validation failed.' };
  page.app.lastRunInfo = { windowSize: Number.NaN };
  page.emit('pageerror', new Error('boom'));
  page.emit('console', { type: () => 'error', text: () => 'Failed to load resource' });
  page.emit('console', { type: () => 'warning', text: () => 'a warning is not an error' });
  assert.deepEqual(await takeOperationHealthProblems(page), [
    'VALIDATION_UNCLASSIFIED diagnostic: Input validation failed.',
    'uncaught page error: boom',
    'console error: Failed to load resource',
    'NaN or undefined in Run Info: lastRunInfo.windowSize=NaN'
  ]);
  // The errors were taken; the alert and Run Info are still there.
  assert.deepEqual(await takeOperationHealthProblems(page), [
    'VALIDATION_UNCLASSIFIED diagnostic: Input validation failed.',
    'NaN or undefined in Run Info: lastRunInfo.windowSize=NaN'
  ]);
  page.app.lastRunInfo = null;
  const errorSignatureBefore = JSON.stringify(['VALIDATION_UNCLASSIFIED', 'generate', 'request-validation', 'Input validation failed.']);
  assert.deepEqual(await takeOperationHealthProblems(page, { errorSignatureBefore }), [], 'an alert an earlier operation left');
  page.app.errorLog = { code: 'UNKNOWN', operation: 'generate', stage: 'render', summary: 'failed' };
  await assert.rejects(assertOperationHealthStrict(page, { operation: 'step 2', errorSignatureBefore }), /step 2 must leave no UNKNOWN/);
  page.app.errorLog = { code: 'LABELS_DO_NOT_FIT', operation: 'generate', stage: 'render', summary: 'Labels do not fit' };
  assert.deepEqual(await takeOperationHealthProblems(page), [], 'a classified alert is the operation\'s answer, not a health problem');
});
